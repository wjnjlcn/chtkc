#include <pthread.h>
#include <stdlib.h>

#include "check_all.h"
#include "../src/hash_map.h"

#define THREAD_COUNT 4


static KC__MemAllocator* ma;
static KC__HashMap* hm;
static size_t max_key_count;
static size_t unique_kmers_count;
static KC__count_t* count_array_in_hash;
static KC__count_t* count_array_out_hash;
static KC__unit_t* threads_kmers[THREAD_COUNT];
static size_t exported_count;


static void setup() {
    ma = KC__mem_allocator_create(1000000);
    hm = KC__hash_map_create(ma, 16, THREAD_COUNT);

    max_key_count = KC__hash_map_max_key_count(hm);
    unique_kmers_count = max_key_count * 2;

    count_array_in_hash = (KC__count_t*)malloc(sizeof(KC__count_t) * unique_kmers_count);
    count_array_out_hash = (KC__count_t*)malloc(sizeof(KC__count_t) * unique_kmers_count);
    for (size_t i = 0; i < unique_kmers_count; i++) {
        count_array_in_hash[i] = 0;
        count_array_out_hash[i] = 0;
    }

    for (size_t i = 0; i < THREAD_COUNT; i++) {
        threads_kmers[i] = (KC__unit_t*)malloc(sizeof(KC__unit_t) * unique_kmers_count);
        for (size_t n = 0; n < unique_kmers_count; n++) {
            threads_kmers[i][n] = n;
        }
    }

    exported_count = 0;
}

static void teardown() {
    KC__hash_map_free(ma, hm);
    KC__mem_allocator_free(ma);

    free(count_array_in_hash);
    free(count_array_out_hash);
    for (size_t i = 0; i < THREAD_COUNT; i++) {
        free(threads_kmers[i]);
    }
}

static void realloc_by_mem_limit(size_t mem_limit) {
    KC__hash_map_free(ma, hm);
    KC__mem_allocator_free(ma);

    ma = KC__mem_allocator_create(mem_limit);
    hm = KC__hash_map_create(ma, 16, THREAD_COUNT);

    max_key_count = KC__hash_map_max_key_count(hm);
    unique_kmers_count = max_key_count * 2;
}

static void randomize_thread_kmers(int n) {
    srandom((unsigned int)time(NULL) + n);
    for (size_t i = 0; i < THREAD_COUNT; i++) {
        for (size_t j = 0; j < unique_kmers_count; j++) {
            long rand_i = random() % THREAD_COUNT;
            long rand_j = random() % unique_kmers_count;
            KC__unit_t tmp = threads_kmers[i][j];
            threads_kmers[i][j] = threads_kmers[rand_i][rand_j];
            threads_kmers[rand_i][rand_j] = tmp;
        }
    }
}

static void* add_kmers(void* ptr) {
    size_t n = *((size_t *)ptr);

    KC__unit_t* kmers = threads_kmers[n];

    for (size_t m = 0; m < 2; m++) {
        for (size_t i = 0; i < unique_kmers_count; i++) {
            KC__unit_t kmer = kmers[i];

            if (!KC__hash_map_add_kmer(hm, n, &kmer)) {
                if (kmer >= unique_kmers_count) {
                    ck_abort();
                }
                __sync_fetch_and_add(&(count_array_out_hash[kmer]), 1);
            }
        }
    }

    KC__hash_map_finish_adding_kmers(hm, n);

    pthread_exit(NULL);
}

static void add_all_kmers() {
    pthread_t threads[THREAD_COUNT];
    size_t thread_ids[THREAD_COUNT];

    for (size_t i = 0; i < THREAD_COUNT; i++) {
        thread_ids[i] = i;
        pthread_create(&(threads[i]), NULL, add_kmers, &(thread_ids[i]));
    }

    for (size_t i = 0; i < THREAD_COUNT; i++) {
        pthread_join(threads[i], NULL);
    }
}

static void check_test_export_callback(const KC__unit_t* kmer, KC__count_t count, void* data) {
    if (data != NULL) {
        ck_abort();
    }

    KC__unit_t idx = kmer[0];
    if (idx >= unique_kmers_count) {
        ck_abort();
    }

    __sync_fetch_and_add(&(count_array_in_hash[idx]), count);
}

static void* export_kmers(void* ptr) {
    size_t thread_id = *((size_t*)ptr);

    size_t ec;
    KC__hash_map_export(hm, thread_id, check_test_export_callback, NULL, &ec);
    __sync_fetch_and_add(&exported_count, ec);

    pthread_exit(NULL);
}

static void export_all_kmers() {
    pthread_t threads[THREAD_COUNT];
    size_t thread_ids[THREAD_COUNT];

    for (size_t i = 0; i < THREAD_COUNT; i++) {
        thread_ids[i] = i;
        pthread_create(&(threads[i]), NULL, export_kmers, &(thread_ids[i]));
    }

    for (size_t i = 0; i < THREAD_COUNT; i++) {
        pthread_join(threads[i], NULL);
    }
}

static void check_results() {
    bool nodes_used_out;
    if (unique_kmers_count == max_key_count * 2) {
        nodes_used_out = true;
    } else if (unique_kmers_count == max_key_count / 2) {
        nodes_used_out = false;
    } else {
        ck_assert(false);
    }

    export_all_kmers();

    if (nodes_used_out) {
        ck_assert_msg(exported_count > max_key_count - THREAD_COUNT, "used nodes: %zu, max key count: %zu", exported_count, max_key_count);
    }
    ck_assert(exported_count <= max_key_count);

    for (size_t i = 0; i < unique_kmers_count; i++) {
        KC__count_t c1 = count_array_in_hash[i];
        KC__count_t c2 = count_array_out_hash[i];

        if (c1 + c2 != (THREAD_COUNT * 2)) {
            ck_abort_msg("%zu, in hash: %zu, out hash: %zu", i, c1, c2);
        }

        if (nodes_used_out) {
            if ((c1 != 0) && (c2 != 0))
                ck_abort_msg("%zu, in hash: %zu, out hash: %zu", i, c1, c2);

        } else {
            if (c2 != 0) {
                ck_abort_msg("%zu, in hash: %zu, out hash: %zu", i, c1, c2);
            }
        }
    }
}

START_TEST(test_table_capacity_one)
    {
        realloc_by_mem_limit(30000);

        // All K-mers will map to table position 0.
        KC__hash_map_set_table_capacity(hm, 1);

        randomize_thread_kmers(_i);

        add_all_kmers();
        check_results();
    }
END_TEST

START_TEST(test_normal_case)
    {
        randomize_thread_kmers(_i);

        for (size_t i = 0; i < 3; i++) {
            add_all_kmers();
            KC__hash_map_clear(hm);
        }
        for (size_t i = 0; i < unique_kmers_count; i++) {
            count_array_out_hash[i] = 0;
        }

        add_all_kmers();
        check_results();
    }
END_TEST

START_TEST(test_use_half_nodes)
    {
        unique_kmers_count = max_key_count / 2;
        randomize_thread_kmers(_i);

        add_all_kmers();
        check_results();
    }
END_TEST

START_TEST(test_rigorous)
    {
        add_all_kmers();
        check_results();
    }
END_TEST

static void test_export_count_callback(const KC__unit_t* kmer, KC__count_t count, void* data) {
    size_t* m = data;
    if (m == NULL) {
        ck_abort();
    }

    if ((kmer == NULL) || (count == 0)) {
        ck_abort();
    }

    (*m)++;
}

START_TEST(test_export_count)
    {
        randomize_thread_kmers(_i);
        add_all_kmers();

        size_t m;
        size_t exported_count;
        for (size_t i = 0; i < THREAD_COUNT; i++) {
            m = 0;
            KC__hash_map_export(hm, i, test_export_count_callback, &m, &exported_count);
            ck_assert(m == exported_count);
        }
    }
END_TEST


Suite* hash_map_suite() {
    TCase* tc_core = tcase_create("Core");
    tcase_add_checked_fixture(tc_core, setup, teardown);
    tcase_add_loop_test(tc_core, test_table_capacity_one, 0, 5);
    tcase_add_loop_test(tc_core, test_normal_case, 0, 5);
    tcase_add_loop_test(tc_core, test_use_half_nodes, 0, 5);
    tcase_add_loop_test(tc_core, test_export_count, 0, 5);

    // tcase_add_loop_test(tc_core, test_rigorous, 0, 1000);


    Suite* s = suite_create("Hash Map");
    suite_add_tcase(s, tc_core);

    return s;
}