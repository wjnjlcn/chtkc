#include <pthread.h>
#include <stdbool.h>
#include <stdlib.h>

#include "check_all.h"
#include "../src/buffer_queue.h"


static KC__MemAllocator* ma;
static KC__BufferQueue* bq;
static const char data[10] = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'};
static int produced_count;
static int consumed_data_count[10];
static int consumed_count;
static pthread_mutex_t produce_mtx;
static pthread_mutex_t consume_mtx;


static void setup() {
    ma = KC__mem_allocator_create(1000000);
    const size_t buffer_size = 5;
    const size_t buffers_count = 3;
    bq = KC__buffer_queue_create(ma, buffer_size, buffers_count);

    produced_count = 0;
    consumed_count = 0;

    for (int i = 0; i < 10; i++)
        consumed_data_count[i] = 0;

    pthread_mutex_init(&produce_mtx, NULL);
    pthread_mutex_init(&consume_mtx, NULL);
}

static void teardown() {
    KC__buffer_queue_free(ma, bq);
    pthread_mutex_destroy(&produce_mtx);
    pthread_mutex_destroy(&consume_mtx);
    KC__mem_allocator_free(ma);
}


static void* produce(void* ptr) {
    KC__BufferQueue* bq = ptr;

    while (true) {
        int i = -1;
        pthread_mutex_lock(&produce_mtx);
        if (produced_count != 10) {
            i = produced_count;
            produced_count++;
        }
        pthread_mutex_unlock(&produce_mtx);

        if (i == -1)
            break;

        KC__Buffer* buffer =  KC__buffer_queue_get_blank_buffer(bq);
        ck_assert(buffer->length == 0);
        ((char*)(buffer->data))[0] = data[i];
        buffer->length = 1;
        KC__buffer_queue_enqueue_filled_buffer(bq, buffer);
    }

    pthread_exit(NULL);
}

static void* consume(void* ptr) {
    KC__BufferQueue* bq = ptr;

    while (true) {
        KC__Buffer* buffer = KC__buffer_queue_dequeue_filled_buffer(bq);
        if (buffer == NULL)
            break;

        ck_assert(buffer->length == 1);

        char c = ((char*)(buffer->data))[0];
        for (int i = 0; i < 10; i++) {
            if (c == data[i]) {
                pthread_mutex_lock(&consume_mtx);
                consumed_data_count[i]++;
                pthread_mutex_unlock(&consume_mtx);
            }
        }

        pthread_mutex_lock(&consume_mtx);
        consumed_count++;
        pthread_mutex_unlock(&consume_mtx);

        KC__buffer_queue_recycle_blank_buffer(bq, buffer);
    }

    pthread_exit(NULL);
}

static void test_function(int p_count, int c_count) {
    pthread_t* p_threads = (pthread_t*)malloc(sizeof(pthread_t) * p_count);
    pthread_t* c_threads = (pthread_t*)malloc(sizeof(pthread_t) * c_count);

    KC__buffer_queue_start_input(bq);

    for (int i = 0; i < p_count; i++)
        pthread_create(&(p_threads[i]), NULL, produce, bq);
    for (int i = 0; i < c_count; i++)
        pthread_create(&(c_threads[i]), NULL, consume, bq);

    for (int i = 0; i < p_count; i++)
        pthread_join(p_threads[i], NULL);

    KC__buffer_queue_finish_input(bq);

    for (int i = 0; i < c_count; i++)
        pthread_join(c_threads[i], NULL);

    free(p_threads);
    free(c_threads);

    ck_assert(consumed_count == 10);

    for (int i = 0; i < 10; i++) {
        ck_assert(consumed_data_count[i] == 1);
    }
}

START_TEST(test_one_producer_multiple_consumers)
    {
        test_function(1, _i);
    }
END_TEST

START_TEST(test_multiple_producers_one_consumer)
    {
        test_function(_i, 1);
    }
END_TEST

START_TEST(test_multiple_producers_multiple_consumers)
    {
        test_function(_i, _i);
    }
END_TEST

Suite* buffer_queue_suite() {
    TCase* tc_core = tcase_create("Core");
    tcase_add_checked_fixture(tc_core, setup, teardown);
    tcase_add_loop_test(tc_core, test_one_producer_multiple_consumers, 1, 21);
    tcase_add_loop_test(tc_core, test_multiple_producers_one_consumer, 1, 21);
    tcase_add_loop_test(tc_core, test_multiple_producers_multiple_consumers, 1, 21);

    Suite* s = suite_create("Buffer Queue");
    suite_add_tcase(s, tc_core);

    return s;
}