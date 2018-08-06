#include <pthread.h>
#include <stdio.h>

#include "check_all.h"
#include "../src/file_writer.h"
#include "../src/header.h"


static KC__MemAllocator* ma;
static KC__FileWriter* fw;
static KC__BufferQueue* bq;
static const char* kmer_file_name = "../tests/test_files/test_write_kmers";
static const char* super_kmer_file_name = "../tests/test_files/test_write_super_kmers";

static void setup() {
    ma = KC__mem_allocator_create(1000000);

    fw = KC__file_writer_create(ma, kmer_file_name, NULL);

    bq = KC__buffer_queue_create(ma, 20, 10);
    KC__file_writer_link_modules(fw, bq);

    KC__file_writer_update_tmp_file(fw, super_kmer_file_name);
}

static void teardown() {
    if (fw != NULL) {
        KC__file_writer_free(ma, fw);
    }

    KC__buffer_queue_free(ma, bq);
    KC__mem_allocator_free(ma);
}

static inline void add_buffers() {
    KC__buffer_queue_start_input(bq);

    unsigned char contents[4][5] = {{0x0, 0x1, 0x2, 0x3},
                                    {0x4, 0x5, 0x6},
                                    {0x8, 0x9, 0xA, 0xB},
                                    {0x7, 0xC, 0xD, 0xE, 0xF}};

    KC__BufferType types[4] = {KC__BUFFER_TYPE_KMER,
                               KC__BUFFER_TYPE_SUPER_KMER,
                               KC__BUFFER_TYPE_KMER,
                               KC__BUFFER_TYPE_SUPER_KMER};

    uint32_t lengths[4] = {4, 3, 4, 5};

    for (size_t i = 0; i < 4; i++) {
        KC__Buffer* bf = KC__buffer_queue_get_blank_buffer(bq);
        bf->type = types[i];
        bf->length = lengths[i];
        memcpy(bf->data, contents[i], bf->length);
        KC__buffer_queue_enqueue_filled_buffer(bq, bf);
    }

    KC__buffer_queue_finish_input(bq);
}

static inline void check_files() {
    unsigned char contents[2][16] = {{0x0, 0x1, 0x2, 0x3, 0x8, 0x9, 0xA, 0xB},
                                     {0x3, 0x0, 0x0, 0x0, 0x4, 0x5, 0x6, 0x5, 0x0, 0x0, 0x0, 0x7, 0xC, 0xD, 0xE, 0xF}};

    size_t file_sizes[2] = {8, 16};

    const char* file_names[2] = {kmer_file_name, super_kmer_file_name};

    for (size_t i = 0; i < 2; i++) {
        FILE* fp = fopen(file_names[i], "rb");
        ck_assert(fp != NULL);

        fseek(fp, 0, SEEK_END);
        long file_size = ftell(fp);
        ck_assert(file_size == file_sizes[i]);

        rewind(fp);

        unsigned char read_content[50];
        size_t read_size = fread(read_content, 1, file_sizes[i], fp);
        ck_assert(read_size == file_sizes[i]);

        fclose(fp);

        ck_assert(memcmp(contents[i], read_content, file_sizes[i]) == 0);
    }
}

static inline void write_files() {
    pthread_t thread;
    pthread_create(&thread, NULL, KC__file_writer_work, fw);
    pthread_join(thread, NULL);
    size_t tmp_file_size = KC__file_writer_get_tmp_file_size(fw);
    ck_assert(tmp_file_size == 16);
}

START_TEST(test_function)
    {
        add_buffers();
        write_files();
        KC__file_writer_free(ma, fw);
        fw = NULL;
        check_files();
    }
END_TEST

Suite* file_writer_suite() {
    TCase* tc_core = tcase_create("Core");
    tcase_add_checked_fixture(tc_core, setup, teardown);
    tcase_add_test(tc_core, test_function);

    Suite* s = suite_create("File writer");
    suite_add_tcase(s, tc_core);

    return s;
}