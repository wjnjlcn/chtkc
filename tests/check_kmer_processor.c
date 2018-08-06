#include <stdlib.h>
#include <time.h>

#include "check_all.h"
#include "../src/kmer_processor.h"
#include "../src/utils.h"


#define TEST_KMER_CHAR_COUNT 24

static KC__MemAllocator *ma;
static KC__KmerProcessor *kp;
static size_t K;
static KC__Buffer buffer;
static size_t check_test_read_callback_called_times;
static size_t check_test_kmer_callback_called_times;
static char random_read[101];
static char random_read_rc[101];
static unsigned char test_canonical_kmer[TEST_KMER_CHAR_COUNT];

static KC__MemAllocator *ma2;
static size_t test_store_check_buffer_called_times;
static size_t test_export_kmer_buffers_count;

static KC__OutputParam output_param;


static void setup() {
    ma = KC__mem_allocator_create(1000000);
    buffer.data = (char *) malloc(1000);
    buffer.length = 0;

    check_test_read_callback_called_times = 0;
    check_test_kmer_callback_called_times = 0;

    ma2 = KC__mem_allocator_create(1000000);
    test_store_check_buffer_called_times = 0;
    test_export_kmer_buffers_count = 0;

    output_param.filter_min = 1;
    output_param.filter_max = KC__COUNT_MAX;
    output_param.count_max = KC__COUNT_MAX;
}

static void teardown() {
    free(buffer.data);
    KC__kmer_processor_free(ma, kp);
    KC__mem_allocator_free(ma);
    KC__mem_allocator_free(ma2);
}

static void init_kmer_processor_by_K() {
    kp = KC__kmer_processor_create(ma, 0, K, output_param);
}

static void init_kmer_processor_by_K_1() {
    kp = KC__kmer_processor_create(ma, 0, 1, output_param);
}

static void generate_random_read(int n) {
    srandom((unsigned int) time(NULL) + n);
    for (size_t i = 0; i < 100; i++) {
        char ch;
        char rc_ch;
        long r = random() % 4;
        switch (r) {
            case 0:
                ch = 'A';
                rc_ch = 'T';
                break;
            case 1:
                ch = 'C';
                rc_ch = 'G';
                break;
            case 2:
                ch = 'G';
                rc_ch = 'C';
                break;
            case 3:
                ch = 'T';
                rc_ch = 'A';
                break;
            default:
                ck_assert(false);
        }
        random_read[i] = ch;
        random_read_rc[99 - i] = rc_ch;
    }
    random_read[100] = '\0';
    random_read_rc[100] = '\0';
}

static inline void extract_kmer(unsigned char *kmer, const char *read, size_t start, size_t end) {
    unsigned char char_unit = 0;
    size_t m = 0;

    for (size_t i = end; i > start; i--) {
        size_t n = i - 1;

        char code;
        switch (read[n]) {
            case 'A':
            case 'a':
                code = 0x0;
                break;
            case 'C':
            case 'c':
                code = 0x1;
                break;
            case 'G':
            case 'g':
                code = 0x2;
                break;
            case 'T':
            case 't':
                code = 0x3;
                break;
            default:
                ck_assert(false);
        }
        char_unit |= (code << (m * 2));

        m++;
        if ((m == 4) || (n == start)) {
            m = 0;

            *kmer = char_unit;
            kmer++;

            char_unit = 0;
        }
    }
}

static void update_canonical_kmer_from_read(size_t n, size_t K) {
    unsigned char kmer[TEST_KMER_CHAR_COUNT];
    unsigned char rc_kmer[TEST_KMER_CHAR_COUNT];

    for (size_t i = 0; i < TEST_KMER_CHAR_COUNT; i++) {
        kmer[i] = 0;
        rc_kmer[i] = 0;
    }

    extract_kmer(kmer, random_read, n, n + K);
    extract_kmer(rc_kmer, random_read_rc, 100 - n - K, 100 - n);

    unsigned char *kmer_selected = NULL;
    for (int i = TEST_KMER_CHAR_COUNT - 1; i > -1; i--) {
        if (kmer[i] < rc_kmer[i]) {
            kmer_selected = kmer;
            break;
        } else if (kmer[i] > rc_kmer[i]) {
            kmer_selected = rc_kmer;
            break;
        }
    }
    if (kmer_selected == NULL) {
        kmer_selected = kmer;
    }

    for (size_t i = 0; i < TEST_KMER_CHAR_COUNT; i++) {
        test_canonical_kmer[i] = kmer_selected[i];
    }
}

static void copy_text_to_buffer(const char *text) {
    buffer.length = (uint32_t) strlen(text);
    memcpy(buffer.data, text, buffer.length);
}

static void check_test_read_callback(KC__KmerProcessor *kp, const char *read, size_t read_length) {
    ck_assert(kp != NULL);
    ck_assert(read != NULL);

    char r[100];
    memcpy(r, read, read_length);
    r[read_length] = '\0';

    switch (check_test_read_callback_called_times) {
        case 0:
            ck_assert(strcmp(r, "") == 0);
            break;
        case 1:
            ck_assert(strcmp(r, "AACCGGTT") == 0);
            break;
        case 2:
            ck_assert(strcmp(r, "ACGT") == 0);
            break;
        case 3:
            ck_assert(strcmp(r, "\rAA\nCC\r\nGG\rTT\n\n") == 0);
            break;
        case 4:
            ck_assert(strcmp(r, "\r\n\rACGT\r\n\n\rTGCA") == 0);
            break;
        default:
            break;
    }

    check_test_read_callback_called_times++;
}


static void check_test_kmer_callback_short_read(KC__KmerProcessor *kp, const KC__unit_t *canonical_kmer, size_t n, KC__unit_t last_code) {
    ck_assert(kp != NULL);
    ck_assert(K == 3);

    switch (check_test_kmer_callback_called_times) {
        case 0:
            ck_assert(n == 0);
            ck_assert(canonical_kmer[0] == 0x6);
            ck_assert(last_code == 0x2);
            break;
        case 1:
            ck_assert(n == 1);
            ck_assert(canonical_kmer[0] == 0x16);
            ck_assert(last_code == 0x2);
            break;
        case 2:
            ck_assert(n == 2);
            ck_assert(canonical_kmer[0] == 0x25);
            ck_assert(last_code == 0x1);
            break;
        case 3:
            ck_assert(n == 0);
            ck_assert(canonical_kmer[0] == 0x2C);
            ck_assert(last_code == 0x0);
            break;
        case 4:
            ck_assert(n == 1);
            ck_assert(canonical_kmer[0] == 0x30);
            ck_assert(last_code == 0x0);
            break;
        case 5:
            ck_assert(n == 0);
            ck_assert(canonical_kmer[0] == 0x5);
            ck_assert(last_code == 0x1);
            break;
        case 6:
            ck_assert(n == 1);
            ck_assert(canonical_kmer[0] == 0x16);
            ck_assert(last_code == 0x2);
            break;
        default:
            break;
    }

    check_test_kmer_callback_called_times++;
}

static void check_test_kmer_callback_long_read(KC__KmerProcessor *kp, const KC__unit_t *canonical_kmer, size_t n, KC__unit_t last_code) {
    ck_assert(kp != NULL);
    ck_assert(K == 67);

    switch (check_test_kmer_callback_called_times) {
        case 0:
            ck_assert(n == 0);
            ck_assert(canonical_kmer[0] == 0x2D4B91BC1B29B06B);
            ck_assert(canonical_kmer[1] == 0x6B1AC6FA16F4AD58);
            ck_assert(canonical_kmer[2] == 0x5);
            ck_assert(last_code == 0x3);
            break;
        case 1:
            ck_assert(n == 1);
            ck_assert(canonical_kmer[0] == 0x85E06B506C5B16AC);
            ck_assert(canonical_kmer[1] == 0xF1971BC1B91E87DA);
            ck_assert(canonical_kmer[2] == 0x16);
            ck_assert(last_code == 0x0);
            break;
        case 2:
            ck_assert(n == 2);
            ck_assert(canonical_kmer[0] == 0x1781AD41B16C5AB1);
            ck_assert(canonical_kmer[1] == 0xC65C6F06E47A1F6A);
            ck_assert(canonical_kmer[2] == 0x1B);
            ck_assert(last_code == 0x1);
            break;
        case 3:
            ck_assert(n == 3);
            ck_assert(canonical_kmer[0] == 0x60B52E46F06CA6C1);
            ck_assert(canonical_kmer[1] == 0x15AC6B1BE85BD2B5);
            ck_assert(canonical_kmer[2] == 0xB);
            ck_assert(last_code == 0x3);
            break;
        case 4:
            ck_assert(n == 4);
            ck_assert(canonical_kmer[0] == 0x582D4B91BC1B29B0);
            ck_assert(canonical_kmer[1] == 0xC56B1AC6FA16F4AD);
            ck_assert(canonical_kmer[2] == 0x12);
            ck_assert(last_code == 0x2);
            break;
        default:
            break;
    }

    check_test_kmer_callback_called_times++;
}

static void check_test_kmer_callback_random_read(KC__KmerProcessor *kp, const KC__unit_t *canonical_kmer, size_t n,
                                                 KC__unit_t last_code) {
    ck_assert(kp != NULL);
    ck_assert(last_code <= 0x3);

    size_t kmer_width = KC__calculate_kmer_width(K);

    update_canonical_kmer_from_read(n, K);
    KC__unit_t *test_kmer = (KC__unit_t *) test_canonical_kmer;

    for (size_t i = 0; i < kmer_width; i++) {
        if (canonical_kmer[i] != test_kmer[i]) {
            ck_abort_msg("Error checking read: %s", random_read);
        }
    }

    check_test_kmer_callback_called_times++;
}


START_TEST(test_handle_buffer_fasta)
    {
        buffer.type = KC__BUFFER_TYPE_FASTA;

        init_kmer_processor_by_K_1();
        KC__kmer_processor_set_read_callback(kp, check_test_read_callback);

        const char* fasta_text;
        switch (_i) {
            case 0:
                fasta_text = "---\n>\n\n" ">-->-\nAACCGGTT\n" ">\nACGT\n" ">\n" "\rAA\nCC\r\nGG\rTT\n\n" "\n>\n" "\r\n\rACGT\r\n\n\rTGCA";
                break;
            case 1:
                fasta_text = "---\r>\r\r" ">-->-\rAACCGGTT\r" ">\rACGT\r" ">\r" "\rAA\nCC\r\nGG\rTT\n\n" "\r>\r" "\r\n\rACGT\r\n\n\rTGCA" "\r";
                break;
            case 2:
                fasta_text = "---\r\n>\r\n\r\n" ">-->-\r\nAACCGGTT\r\n" ">\r\nACGT\r\n" ">\r\n" "\rAA\nCC\r\nGG\rTT\n\n" "\r\n>\r\n" "\r\n\rACGT\r\n\n\rTGCA" "\r\n";
                break;
            default:
                ck_abort();
        }
        copy_text_to_buffer(fasta_text);

        KC__kmer_processor_handle_buffer(kp, &buffer);

        ck_assert_msg(check_test_read_callback_called_times == 5);
    }
END_TEST

START_TEST(test_handle_buffer_fastq)
    {
        buffer.type = KC__BUFFER_TYPE_FASTQ;

        init_kmer_processor_by_K_1();
        KC__kmer_processor_set_read_callback(kp, check_test_read_callback);

        const char *fastq_text = "---\n@\n\n+\n+\n" "@--@-\nAACCGGTT\n+\n@---\n" "@\nACGT\n+";
        copy_text_to_buffer(fastq_text);

        KC__kmer_processor_handle_buffer(kp, &buffer);

        ck_assert(check_test_read_callback_called_times == 3);
    }
END_TEST

START_TEST(test_handle_buffer_super_kmer_1)
    {
        K = 3;
        init_kmer_processor_by_K();
        KC__kmer_processor_set_kmer_callback(kp, check_test_kmer_callback_short_read);

        const uint8_t content[11] = {0x3, 0x0, 0x0, 0x0, 0x2, 0xA4, 0x1, 0x1, 0xE, 0x1, 0x94};
        buffer.length = sizeof(uint8_t) * 11;
        memcpy(buffer.data, content, buffer.length);
        buffer.type = KC__BUFFER_TYPE_SUPER_KMER;

        KC__kmer_processor_handle_buffer(kp, &buffer);

        ck_assert(check_test_kmer_callback_called_times == 7);
    }
END_TEST

START_TEST(test_handle_buffer_super_kmer_2)
    {
        K = 67;
        init_kmer_processor_by_K();
        KC__kmer_processor_set_kmer_callback(kp, check_test_kmer_callback_long_read);

        const uint8_t content[23] = {0x1, 0x0, 0x0, 0x0, 0x4, 0x94, 0x4F, 0xD6, 0xE4, 0x43, 0x6E, 0xB4, 0xD2, 0xA7,
                                     0x52, 0x0B, 0xE9, 0x05, 0x39, 0xE5, 0x94, 0x3A, 0x2D};
        buffer.length = sizeof(uint8_t) * 23;
        memcpy(buffer.data, content, buffer.length);
        buffer.type = KC__BUFFER_TYPE_SUPER_KMER;

        KC__kmer_processor_handle_buffer(kp, &buffer);

        ck_assert(check_test_kmer_callback_called_times == 5);
    }
END_TEST

START_TEST(test_handle_read_short_read)
    {
        K = 3;
        init_kmer_processor_by_K();
        KC__kmer_processor_set_kmer_callback(kp, check_test_kmer_callback_short_read);

        const char *read = "NA\rC\r\nGGCNG\nCNGTAANNACCGNNN";
        size_t read_length = strlen(read);
        KC__kmer_processor_handle_read(kp, read, read_length);

        ck_assert(check_test_kmer_callback_called_times == 7);
    }
END_TEST

START_TEST(test_handle_read_long_read)
    {
        ck_assert(sizeof(KC__unit_t) == sizeof(uint64_t));

        K = 67;
        init_kmer_processor_by_K();
        KC__kmer_processor_set_kmer_callback(kp, check_test_kmer_callback_long_read);

        const char *read = "ACCG\rTTACG\r\nCCTACGTTAAC\nGTGCACTGGACT\n\r\nTCGGGACCTGAAC\n\nGGTCCAACGT\nACCGTACCGGGTACTG";
        size_t read_length = strlen(read);

        KC__kmer_processor_handle_read(kp, read, read_length);

        ck_assert(check_test_kmer_callback_called_times == 5);
    }
END_TEST

START_TEST(test_handle_read_random_read)
    {
        generate_random_read(_i);

        K = (size_t) _i;
        init_kmer_processor_by_K();
        KC__kmer_processor_set_kmer_callback(kp, check_test_kmer_callback_random_read);

        size_t read_length = strlen(random_read);
        KC__kmer_processor_handle_read(kp, random_read, read_length);

        ck_assert(check_test_kmer_callback_called_times == 101 - K);
    }
END_TEST

static void test_store_add_kmers(KC__HashMap *hm) {
    char long_read[285];
    for (size_t i = 0; i < 71; i++) {
        long_read[i * 4 + 0] = 'A';
        long_read[i * 4 + 1] = 'G';
        long_read[i * 4 + 2] = 'C';
        long_read[i * 4 + 3] = 'T';
    }
    long_read[284] = '\0';

    const char *reads[] = {"AAAAAA", "ACGTAAAAAACGTACG", long_read, "CGAT"};

    for (size_t i = 0; i < 4; i++) {
        KC__kmer_processor_handle_read(kp, reads[i], strlen(reads[i]));
        if (i == 0) {
            KC__hash_map_lock_keys(hm);
        }
    }

    KC__kmer_processor_finish(kp);
}

static KC__Buffer *test_store_alloc_buffer() {
    KC__Buffer *bf = (KC__Buffer *) KC__mem_alloc(ma2, sizeof(KC__Buffer), "test store buffer");
    bf->size = sizeof(uint8_t) * 77;
    bf->length = 0;
    bf->data = KC__mem_alloc(ma2, bf->size, "test store buffer data");
    return bf;
}

static void test_store_check_buffer(KC__Buffer *bf) {
    ck_assert(bf->type == KC__BUFFER_TYPE_SUPER_KMER);

    uint8_t data[77];
    const uint8_t tmp_1[12] = {0x3, 0x0, 0x0, 0x0, 0x3, 0xE4, 0x0, 0x5, 0x40, 0x4E, 0x2, 0xFF};
    const uint8_t tmp_2[5] = {0x1, 0x0, 0x0, 0x0, 0x18};
    const uint8_t tmp_3[6] = {0x1, 0x0, 0x0, 0x0, 0x0, 0xC9};

    switch (test_store_check_buffer_called_times) {
        case 0:
            ck_assert(bf->length == sizeof(uint8_t) * 77);
            memcpy(data, tmp_1, sizeof(uint8_t) * 12);
            for (size_t i = 12; i < 76; i++) {
                data[i] = 0xD8;
            }
            data[76] = 0x18;
            break;
        case 1:
            ck_assert(bf->length == sizeof(uint8_t) * 12);
            memcpy(data, tmp_2, sizeof(uint8_t) * 5);
            for (size_t i = 5; i < 12; i++) {
                data[i] = 0xD8;
            }
            break;
        case 2:
            ck_assert(bf->length == sizeof(uint8_t) * 6);
            memcpy(data, tmp_3, sizeof(uint8_t) * 6);
            break;
        default:
            ck_assert(false);
    }

    ck_assert(memcmp(data, bf->data, bf->length) == 0);

    KC__mem_free(ma2, bf->data);
    KC__mem_free(ma2, bf);

    test_store_check_buffer_called_times++;
}

START_TEST(test_store)
    {
        K = 4;
        init_kmer_processor_by_K();
        KC__kmer_processor_set_store_buffer_request_callback(kp, test_store_alloc_buffer);
        KC__kmer_processor_set_store_buffer_complete_callback(kp, test_store_check_buffer);

        KC__HashMap *hm = KC__hash_map_create(ma, K, 1);
        KC__kmer_processor_link_modules(kp, hm, NULL, NULL);

        test_store_add_kmers(hm);

        ck_assert(test_store_check_buffer_called_times == 3);

        KC__hash_map_free(ma, hm);
    }
END_TEST


static KC__Buffer* test_export_alloc_buffer() {
    KC__Buffer *bf = (KC__Buffer *) KC__mem_alloc(ma2, sizeof(KC__Buffer), "test store buffer");
    bf->size = 32;
    bf->length = 0;
    bf->data = KC__mem_alloc(ma2, bf->size, "test store buffer data");
    return bf;
}

static void test_export_check_buffer(KC__Buffer *bf) {
    ck_assert(bf->type == KC__BUFFER_TYPE_KMER);

    unsigned char data[32];
    for (size_t i = 0; i < 32; i++) {
        data[i] = 0;
    }

    size_t length;

    switch (test_export_kmer_buffers_count) {
        case 0:
            data[0] = 0x16;
            data[1] = 0x2;
            data[5] = 0x5A;
            data[6] = 0x3;
            data[10] = 0x1B;
            data[11] = 0x1;
            data[15] = 0x25;
            data[16] = 0x1;
            data[20] = 0x95;
            data[21] = 0x1;
            data[25] = 0x55;
            data[26] = 0x1;
            length = 30;
            break;
        case 1:
            data[0] = 0x56;
            data[1] = 0x2;
            data[5] = 0x36;
            data[6] = 0x1;
            length = 10;
            break;
        default:
            ck_abort();
    }

    test_export_kmer_buffers_count++;

    ck_assert(bf->length == length);
    ck_assert(memcmp(data, bf->data, bf->length) == 0);

    KC__mem_free(ma2, bf->data);
    KC__mem_free(ma2, bf);
}

START_TEST(test_export)
    {
        output_param.count_max = UINT32_MAX;

        K = 4;
        init_kmer_processor_by_K();
        KC__kmer_processor_set_store_buffer_request_callback(kp, test_export_alloc_buffer);
        KC__kmer_processor_set_store_buffer_complete_callback(kp, test_export_check_buffer);

        KC__HashMap *hm = KC__hash_map_create(ma, K, 1);
        KC__kmer_processor_link_modules(kp, hm, NULL, NULL);

        const char* reads[6] = {"ACCGG", "ACGT", "ACCGG", "AGCCCCGG", "CCCG", "ATCG"};
        for (size_t i = 0; i < 6; i++) {
            KC__kmer_processor_handle_read(kp, reads[i], strlen(reads[i]));
        }
        KC__kmer_processor_finish(kp);

        KC__kmer_processor_export_kmers(kp);

        ck_assert(test_export_kmer_buffers_count == 2);

        KC__hash_map_free(ma, hm);
    }
END_TEST

static void test_export_check_buffer_2(KC__Buffer *bf) {
    ck_assert(bf->type == KC__BUFFER_TYPE_KMER);

    unsigned char data[32];
    for (size_t i = 0; i < 32; i++) {
        data[i] = 0;
    }


    size_t length;

    unsigned char codes[9] = {0x69, 0x79, 0xE4, 0x06, 0x6F, 0x5C, 0xC6, 0x5B, 0x01};
    unsigned char codes_2[9] = {0x79, 0xE4, 0x06, 0x6F, 0x5C, 0xC6, 0x5B, 0xA5, 0x01};

    for (size_t i = 0; i < 9 ; i++) {
        data[i] = codes[i];
    }
    length = 9;

    switch (output_param.count_max) {
        case 255:
            data[9] = 0xFF;
            length += 1;
            break;
        case 300:
            data[9] = 0x2C;
            data[10] = 0x01;
            length += 2;
            break;
        case UINT32_MAX:
            data[9] = 0xF4;
            data[10] = 0x01;
            length += 4;
            break;
        default:
            ck_abort();
    }

    if (output_param.filter_min == 1) {
        for (size_t i = 0; i < 9; i++) {
            data[length + i] = codes_2[i];
        }
        length += 9;

        switch (output_param.count_max) {
            case 255:
                data[length] = 0x01;
                length += 1;
                break;
            case 300:
                data[length] = 0x01;
                length += 2;
                break;
            default:
                ck_abort();
        }

    }


    test_export_kmer_buffers_count++;

    ck_assert(bf->length == length);
    ck_assert(memcmp(data, bf->data, bf->length) == 0);

    KC__mem_free(ma2, bf->data);
    KC__mem_free(ma2, bf);
}

START_TEST(test_export_2)
    {
        switch (_i) {
            case 0:
                output_param.count_max = 255;
                break;
            case 1:
                output_param.count_max = 300;
                break;
            case 2:
                output_param.count_max = UINT32_MAX;
                output_param.filter_min = 2;
                break;
            default:
                ck_abort();
        }

        K = 33;
        init_kmer_processor_by_K();
        KC__kmer_processor_set_store_buffer_request_callback(kp, test_export_alloc_buffer);
        KC__kmer_processor_set_store_buffer_complete_callback(kp, test_export_check_buffer_2);

        KC__HashMap *hm = KC__hash_map_create(ma, K, 1);
        KC__kmer_processor_link_modules(kp, hm, NULL, NULL);

        const char* read = "CCCGTTACGCCTACGTTAACGTGCACTGCCGGC";
        for (size_t i = 0; i < 500; i++) {
            KC__kmer_processor_handle_read(kp, read, strlen(read));
        }
        read = "CGGCCCCGTTACGCCTACGTTAACGTGCACTGC";
        KC__kmer_processor_handle_read(kp, read, strlen(read));

        KC__kmer_processor_finish(kp);

        KC__kmer_processor_export_kmers(kp);

        ck_assert(test_export_kmer_buffers_count == 1);

        KC__hash_map_free(ma, hm);
    }
END_TEST


Suite *kmer_processor_suite() {
    TCase *tc_core = tcase_create("Core");
    tcase_add_checked_fixture(tc_core, setup, teardown);

    tcase_add_loop_test(tc_core, test_handle_buffer_fasta, 0, 3);
    tcase_add_test(tc_core, test_handle_buffer_fastq);

    tcase_add_test(tc_core, test_handle_read_short_read);
    tcase_add_test(tc_core, test_handle_read_long_read);
    tcase_add_loop_test(tc_core, test_handle_read_random_read, 1, 97);

    tcase_add_test(tc_core, test_handle_buffer_super_kmer_1);
    tcase_add_test(tc_core, test_handle_buffer_super_kmer_2);

    tcase_add_test(tc_core, test_store);
    tcase_add_test(tc_core, test_export);
    tcase_add_loop_test(tc_core, test_export_2, 0, 3);

    Suite *s = suite_create("Kmer Processor");
    suite_add_tcase(s, tc_core);

    return s;
}


/*
 * An example of K-mer extracting process.
 *
 * K = 8, W = 3, BIT = 6
 *
 *
 * d = [1 2 3 4 5 6 7 8]
 *
 * _0_1_2 _3_4_5 _6_7_8 -> kmer = _6_7_8 _3_4_5 _0_1_2
 * _0_8_7 _6_5_4 _3_2_1 -> rc_kmer = _3_2_1 _6_5_4 _0_8_7
 *
 *
 * d = [2 3 4 5 6 7 8 9]
 *
 * _0_2_3 _4_5_6 _7_8_9 -> kmer = _7_8_9 _4_5_6 _0_2_3
 * _0_9_8 _7_6_5 _4_3_2 -> rc_kmer = _4_3_2 _7_6_5 _0_9_8
 *
 *
 * d = [3 4 5 6 7 8 9 A]
 *
 * _0_3_4 _5_6_7 _8_9_A -> kmer = _8_9_A _5_6_7 _0_3_4
 * _0_A_9 _8_7_6 _5_4_3 -> rc_kmer = _5_4_3 _8_7_6 _0_A_9
 *
 */