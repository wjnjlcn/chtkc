#include <pthread.h>

#include "check_all.h"
#include "../src/file_reader.h"
#include "../src/buffer_queue.h"


static KC__MemAllocator* ma;
static KC__FileReader* fr;
static KC__BufferQueue* bq;


static void setup() {
    ma = KC__mem_allocator_create(1000000);
    // Create with type gzip ensures gzip buffer exists.
    fr = KC__file_reader_create(ma, 4, KC__FILE_COMPRESSION_TYPE_GZIP, 3);
    bq = KC__buffer_queue_create(ma, 20, 10);
    KC__file_reader_link_modules(fr, bq);
}

static void teardown() {
    KC__file_reader_free(ma, fr);
    KC__buffer_queue_free(ma, bq);
    KC__mem_allocator_free(ma);
}

static inline void read_files() {
    KC__buffer_queue_start_input(bq);

    pthread_t thread;
    pthread_create(&thread, NULL, KC__file_reader_work, fr);
    pthread_join(thread, NULL);

    KC__buffer_queue_finish_input(bq);
}

static void check_buffers(void (*check_buffer) (KC__Buffer*, size_t), size_t expected_buffers_count) {
    KC__Buffer* bf;
    size_t n = 0;
    while (1) {
        bf = KC__buffer_queue_dequeue_filled_buffer(bq);
        if (bf == NULL)
            break;

        check_buffer(bf, n);

        n++;

        KC__buffer_queue_recycle_blank_buffer(bq, bf);
    }

    ck_assert(n == expected_buffers_count);
}

static void check_fasta_buffer(KC__Buffer* bf, size_t i) {
    ck_assert(bf->type == KC__BUFFER_TYPE_FASTA);

    const char* content;
    switch (i) {
        case 0:
            content = ">1\nACGTA\n>2\nTCGAT\n";
            break;
        case 1:
            content = ">\nATCGATCG\nAACNCGNN\n";
            break;
        case 2:
            content = ">\nCCGGTT\n";
            break;
        case 3:
            content = ">\nNNNANNNC\nNNNGNNNN\n";
            break;
        case 4:
            content = ">\nACGT\n";
            break;
        default:
            ck_abort();
    }

    ck_assert(bf->length <= bf->size);
    ck_assert(bf->length == strlen(content));
    ck_assert(memcmp(content, bf->data, bf->length) == 0);
}

START_TEST(test_fasta)
    {
        char* file_names[] = {"../tests/test_files/test_fasta.fa"};
        KC__FileInputDescription input = {file_names, 1, KC__FILE_TYPE_FASTA, KC__FILE_COMPRESSION_TYPE_PLAIN};
        KC__file_reader_update_input(fr, input);
        read_files();

        check_buffers(check_fasta_buffer, 5);
    }
END_TEST

static void check_fastq_buffer(KC__Buffer* bf, size_t i) {
    ck_assert(bf->type == KC__BUFFER_TYPE_FASTQ);

    ck_assert_msg(bf->length == 17);

    const char* content;
    switch (i) {
        case 0:
            content = "@1\nACGTA\n+\n-----\n";
            break;
        case 1:
            content = "@2\nTGCAT\n+\n-----\n";
            break;
        case 2:
            content = "@3\nATCGA\n+\n-----\n";
            break;
        default:
            ck_abort();
    }

    ck_assert(memcmp(content, bf->data, bf->length) == 0);
}

START_TEST(test_fastq)
    {
        char* file_names[] = {"../tests/test_files/test_fastq_1.fq",
                              "../tests/test_files/test_fastq_2.fq"};
        KC__FileInputDescription input = {file_names, 2, KC__FILE_TYPE_FASTQ, KC__FILE_COMPRESSION_TYPE_PLAIN};
        KC__file_reader_update_input(fr, input);
        read_files();

        check_buffers(check_fastq_buffer, 3);
    }
END_TEST

START_TEST(test_gz)
    {
        char* file_names[] = {"../tests/test_files/test_fastq_1.fq.gz",
                              "../tests/test_files/test_fastq_2.fq.gz"};
        KC__FileInputDescription input = {file_names, 2, KC__FILE_TYPE_FASTQ, KC__FILE_COMPRESSION_TYPE_GZIP};
        KC__file_reader_update_input(fr, input);
        read_files();

        check_buffers(check_fastq_buffer, 3);
    }
END_TEST

START_TEST(test_cat_gz)
    {
        char* file_names[] = {"../tests/test_files/test_fastq_cat.fq.gz"};
        KC__FileInputDescription input = {file_names, 1, KC__FILE_TYPE_FASTQ, KC__FILE_COMPRESSION_TYPE_GZIP};
        KC__file_reader_update_input(fr, input);
        read_files();

        check_buffers(check_fastq_buffer, 3);
    }
END_TEST

static void check_super_kmer_buffer(KC__Buffer* bf, size_t i) {
    ck_assert(bf->type == KC__BUFFER_TYPE_SUPER_KMER);

    ck_assert_msg(bf->length == 4);

    unsigned char contents[4][4] = {{0x0, 0x1, 0x2, 0x3},
                                    {0x4, 0x5, 0x6, 0x7},
                                    {0x8, 0x9, 0xA, 0xB},
                                    {0xC, 0xD, 0xE, 0xF}};
    unsigned char* content;

    if (i < 4) {
        content = contents[i];
    } else {
        ck_abort();
    }

    ck_assert(memcmp(content, bf->data, bf->length) == 0);
}

START_TEST(test_super_kmer)
    {
        char* file_names[] = {"../tests/test_files/test_super_kmer"};
        KC__FileInputDescription input = {file_names, 1, KC__FILE_TYPE_SUPER_KMER, KC__FILE_COMPRESSION_TYPE_PLAIN};
        KC__file_reader_update_input(fr, input);
        read_files();

        check_buffers(check_super_kmer_buffer, 4);
    }
END_TEST

Suite* file_reader_suite() {
    TCase* tc_core = tcase_create("Core");
    tcase_add_checked_fixture(tc_core, setup, teardown);
    tcase_add_test(tc_core, test_fasta);
    tcase_add_test(tc_core, test_fastq);
    tcase_add_test(tc_core, test_gz);
    tcase_add_test(tc_core, test_cat_gz);
    tcase_add_test(tc_core, test_super_kmer);

    Suite* s = suite_create("File reader");
    suite_add_tcase(s, tc_core);

    return s;
}