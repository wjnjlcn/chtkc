/*
 * This file is part of CHTKC.
 *
 * CHTKC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CHTKC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CHTKC.  If not, see <https://www.gnu.org/licenses/>.
 *
 * Author: Jianan Wang
 */


#include <inttypes.h>
#include <string.h>
#include <pthread.h>
#include "kmer_processor.h"
#include "logging.h"
#include "assert.h"
#include "utils.h"
#include "param.h"


typedef struct {
    void* mem_block;

    size_t K;
    size_t W;

    size_t gen_w_init;
    size_t gen_s_init;
    size_t rc_gen_w_init;
    size_t rc_gen_s_init;

    size_t gen_w;
    size_t gen_s;
    size_t rc_gen_w;
    size_t rc_gen_s;

    KC__unit_t shift_mask;
    KC__unit_t rc_shift;

    KC__unit_t* kmer;
    KC__unit_t* rc_kmer;

} KC__KmerExtractUnit;


typedef enum {
    KC__KMER_STORE_ACTION_NEW = 0,
    KC__KMER_STORE_ACTION_EXPAND
} KC__KmerStoreAction;

typedef struct {
    KC__KmerStoreAction store_action;

    KC__Buffer* current_buffer;

    // Include expanded_bases_count and super K-mer codes.
    size_t super_kmer_info_max_size;

    uint32_t* super_kmers_count;

    // The length of super K-mer minus K, by base(A, C, G, T) count.
    uint8_t* expanded_bases_count;
    uint8_t* current_unit;
    // The count of bases in current_unit;
    size_t current_bases_count;

} KC__KmerStoreUnit;

typedef struct {
    KC__Buffer* buffer;
    size_t W;
    size_t C;
    size_t unit_size;

    KC__OutputParam output_param;
    size_t count_bit;

    size_t total_kmers_count;
    size_t unique_kmers_count;
    size_t exported_unique_kmers_count;
} KC__KmerExportUnit;


struct KC__KmerProcessor {
    size_t id;

    KC__KmerExtractUnit kmer_extract_unit;
    KC__KmerStoreUnit kmer_store_unit;
    KC__KmerExportUnit kmer_export_unit;

    void* tmp_kmers_mem;

    KC__HashMap* hash_map;
    KC__BufferQueue* read_buffer_queue;
    KC__BufferQueue* write_buffer_queue;

    KC__KmerProcessorReadCallback read_callback;
    KC__KmerProcessorKmerCallback kmer_callback;
    KC__KmerProcessorStoreBufferRequestCallback store_buffer_request_callback;
    KC__KmerProcessorStoreBufferCompleteCallback store_buffer_complete_callback;
};


static inline KC__unit_t KC__encode(char ch) {
    switch (ch) {
        case 'A':
        case 'a':
            return 0x0;
        case 'C':
        case 'c':
            return 0x1;
        case 'G':
        case 'g':
            return 0x2;
        case 'T':
        case 't':
            return 0x3;
        case '\n':
        case '\r':
            return 0x4;
        default:
            return 0x5;
    }
}

static inline KC__unit_t KC__get_rc_code(KC__unit_t code) {
    KC__ASSERT(code <= 0x3);
    return (KC__unit_t)0x3 - code;
}


static void KC__kmer_extract_unit_init(KC__KmerExtractUnit* keu, size_t K) {
    keu->K = K;
    keu->W = KC__calculate_kmer_width(K);

    size_t high_valid_bit = keu->K * 2 - KC__UNIT_BIT * (keu->W - 1);

    keu->gen_w_init = keu->W - 1;
    keu->gen_s_init = high_valid_bit - 2;
    keu->rc_gen_w_init = 0;
    keu->rc_gen_s_init = 0;

    keu->shift_mask = (high_valid_bit == KC__UNIT_BIT) ? KC__UNIT_MAX : (((size_t)1 << high_valid_bit) - 1);
    keu->rc_shift = keu->gen_s_init;
}

static inline void KC__kmer_extract_unit_generate_kmer(KC__KmerExtractUnit* keu, size_t n, KC__unit_t code) {
    if (n == 0) {
        keu->gen_w = keu->gen_w_init;
        keu->gen_s = keu->gen_s_init;
        keu->rc_gen_w = keu->rc_gen_w_init;
        keu->rc_gen_s = keu->rc_gen_s_init;

        for (size_t i = 0; i < keu->W; i++) {
            keu->kmer[i] = 0;
            keu->rc_kmer[i] = 0;
        }
    }

    keu->kmer[keu->gen_w] |= (code << keu->gen_s);
    keu->rc_kmer[keu->rc_gen_w] |= (KC__get_rc_code(code) << keu->rc_gen_s);

    if (keu->gen_s == 0) {
        keu->gen_w--;
        keu->gen_s = KC__UNIT_BIT - 2;
    } else {
        keu->gen_s -= 2;
    }

    if (keu->rc_gen_s == KC__UNIT_BIT - 2) {
        keu->rc_gen_w++;
        keu->rc_gen_s = 0;
    } else {
        keu->rc_gen_s += 2;
    }
}

static inline void KC__kmer_extract_unit_shift_kmer(KC__KmerExtractUnit* keu, KC__unit_t code) {
    size_t W = keu->W;

    KC__unit_t rc_code = KC__get_rc_code(code);

    KC__unit_t tmp_code;
    for (size_t i = 0; i < W - 1; i++) {
        tmp_code = (keu->kmer[i] >> (KC__UNIT_BIT - 2));
        keu->kmer[i] = (keu->kmer[i] << 2) | code;
        code = tmp_code;
    }
    keu->kmer[W - 1] = ((keu->kmer[W - 1] << 2) | code) & (keu->shift_mask);

    for (size_t i = W; i > 0; i--) {
        size_t shift = (i == W) ? keu->rc_shift : KC__UNIT_BIT - 2;
        size_t j = i - 1;
        tmp_code = (keu->rc_kmer[j] & 0x3);
        keu->rc_kmer[j] = ((keu->rc_kmer[j] >> 2) | (rc_code << shift));
        rc_code = tmp_code;
    }
}

static inline int KC__kmer_extract_unit_compare_kmers(KC__KmerExtractUnit* keu, const KC__unit_t* kmer_1, const KC__unit_t* kmer_2) {
    for (size_t i = keu->W; i > 0; i--) {
        size_t n = i - 1;
        if (kmer_1[n] < kmer_2[n]) {
            return -1;
        } else if (kmer_1[n] > kmer_2[n]) {
            return 1;
        }
    }
    return 0;
}

static void KC__kmer_store_unit_init(KC__KmerStoreUnit* ksu, size_t K) {
    ksu->store_action = KC__KMER_STORE_ACTION_NEW;
    ksu->current_buffer = NULL;

    size_t max_units_count = KC__calculate_kmer_width_by_unit_size(K + UINT8_MAX, sizeof(uint8_t));
    ksu->super_kmer_info_max_size = sizeof(uint8_t) * (max_units_count + 1);
}


KC__KmerProcessor* KC__kmer_processor_create(KC__MemAllocator* ma, size_t id, size_t K, KC__OutputParam output_param) {
    KC__KmerProcessor* kp = (KC__KmerProcessor*)KC__mem_aligned_alloc(ma, sizeof(KC__KmerProcessor), "kmer processor");

    kp->id = id;

    KC__kmer_extract_unit_init(&(kp->kmer_extract_unit), K);
    KC__kmer_store_unit_init(&(kp->kmer_store_unit), K);

    size_t kmer_size = KC__calculate_kmer_size(K);
    kp->tmp_kmers_mem = KC__mem_aligned_alloc(ma, kmer_size * 2, "kmer processor tmp kmers mem");
    char* mem = kp->tmp_kmers_mem;
    kp->kmer_extract_unit.kmer = (KC__unit_t*)(mem);
    kp->kmer_extract_unit.rc_kmer = (KC__unit_t*)(mem + kmer_size);

    kp->kmer_export_unit.output_param = output_param;


    kp->hash_map = NULL;
    kp->read_buffer_queue = NULL;
    kp->write_buffer_queue = NULL;

    KC__kmer_processor_set_read_callback(kp, KC__kmer_processor_handle_read);
    KC__kmer_processor_set_kmer_callback(kp, KC__kmer_processor_handle_kmer);

    kp->store_buffer_request_callback = NULL;
    kp->store_buffer_complete_callback = NULL;

    return kp;
}

void KC__kmer_processor_free(KC__MemAllocator* ma, KC__KmerProcessor* kp) {
    KC__mem_free(ma, kp->tmp_kmers_mem);
    KC__mem_free(ma, kp);
}

void KC__kmer_processor_link_modules(KC__KmerProcessor* kp, KC__HashMap* hash_map, KC__BufferQueue* read_buffer_queue, KC__BufferQueue* write_buffer_queue) {
    kp->hash_map = hash_map;
    kp->read_buffer_queue = read_buffer_queue;
    kp->write_buffer_queue = write_buffer_queue;
}

void KC__kmer_processor_set_read_callback(KC__KmerProcessor* kp, KC__KmerProcessorReadCallback read_callback) {
    kp->read_callback = read_callback;
}

void KC__kmer_processor_set_kmer_callback(KC__KmerProcessor* kp, KC__KmerProcessorKmerCallback kmer_callback) {
    kp->kmer_callback = kmer_callback;
}

void KC__kmer_processor_set_store_buffer_request_callback(KC__KmerProcessor* kp, KC__KmerProcessorStoreBufferRequestCallback request_callback) {
    kp->store_buffer_request_callback = request_callback;
}

void KC__kmer_processor_set_store_buffer_complete_callback(KC__KmerProcessor* kp, KC__KmerProcessorStoreBufferCompleteCallback complete_callback) {
    kp->store_buffer_complete_callback = complete_callback;
}

/**
 * Handle code extracted by read (or super K-mer), update K-mer extract unit.
 * If a canonical K-mer is ready, call the K-mer callback.
 * @param kp K-mer processor.
 * @param i The (i)th valid code of read (or super K-mer) is being handled.
 * @param code The code to be handled.
 */
static inline void KC__kmer_processor_handle_code(KC__KmerProcessor* kp, size_t i, KC__unit_t code) {
    KC__KmerExtractUnit* keu = &(kp->kmer_extract_unit);

    if (i < keu->K) {
        KC__kmer_extract_unit_generate_kmer(keu, i, code);

        if (i != keu->K - 1) {
            return;
        }

    } else {
        KC__kmer_extract_unit_shift_kmer(keu, code);
    }

    KC__unit_t *canonical_kmer;
    if (KC__kmer_extract_unit_compare_kmers(keu, keu->kmer, keu->rc_kmer) < 0)
        canonical_kmer = keu->kmer;
    else
        canonical_kmer = keu->rc_kmer;

    kp->kmer_callback(kp, canonical_kmer, i + 1 - keu->K, code);
}

static inline void KC__kmer_processor_handle_reads_buffer(KC__KmerProcessor* kp, const KC__Buffer* buffer) {
    const char* data = buffer->data;

    const char* prev_line = NULL;
    const char* next_line = NULL;

    const char* current_line = data;

    // Start index of a line;
    size_t line_start_idx = 0;

    size_t i = 0;
    while (true) {
        // Only analyze at the end of the line (or buffer).
        const bool end_of_buffer = (i == buffer->length);
        size_t line_end_idx = 0;

        bool end_of_line = false;

        if (end_of_buffer) {
            line_end_idx = i;

        } else {
            if (data[i] == '\n') {
                line_end_idx = i;
                end_of_line = true;

            } else if (data[i] == '\r') {
                line_end_idx = i;
                if ((i < buffer->length - 1) && (data[i + 1] == '\n')) {
                    i++;
                }
                end_of_line = true;
            }
        }

        if (end_of_buffer || end_of_line) {
            KC__ASSERT(line_end_idx >= line_start_idx);
            const size_t current_line_length = line_end_idx - line_start_idx;

            const size_t next_line_start_idx = i + 1;
            if (next_line_start_idx >= buffer->length) {
                next_line = NULL;
            } else {
                next_line = &(data[next_line_start_idx]);
            }

            bool current_line_is_read = false;
            bool update_current_line = false;

            switch (buffer->type) {
                case KC__BUFFER_TYPE_FASTA:
                    if ((prev_line != NULL) && (prev_line[0] == '>')) {
                        if ((next_line == NULL) || (next_line[0] == '>')) {
                            current_line_is_read = true;
                            update_current_line = true;
                        }
                    } else {
                        update_current_line = true;
                    }
                    break;
                case KC__BUFFER_TYPE_FASTQ:
                    // Prev line and next line at least has one character ('\n').
                    if ((prev_line != NULL) && (prev_line[0] == '@')) {
                        if ((next_line != NULL) && (next_line[0] == '+')) {
                            current_line_is_read = true;
                        }
                    }
                    update_current_line = true;
                    break;
                default:
                    KC__ASSERT(false);
                    break;
            }

            if (current_line_is_read) {
                kp->read_callback(kp, current_line, current_line_length);
            }

            if (end_of_buffer) {
                break;
            }

            if (update_current_line) {
                prev_line = current_line;
                line_start_idx = next_line_start_idx;
                current_line = data + line_start_idx;
            }
        }

        i++;
    }
}

static inline void KC__kmer_processor_handle_super_kmers_buffer(KC__KmerProcessor* kp, const KC__Buffer* buffer) {
    KC__ASSERT(buffer->type == KC__BUFFER_TYPE_SUPER_KMER);

    size_t super_kmers_count = *((uint32_t*)(buffer->data));
    uint8_t* p = (uint8_t*)((char*)(buffer->data) + sizeof(uint32_t));

    for (size_t n = 0; n < super_kmers_count; n++) {
        uint8_t expanded_bases_count = *p;
        size_t bases_count = kp->kmer_extract_unit.K + expanded_bases_count;

        uint8_t mask = 0x3;
        KC__unit_t code;

        p++;
        uint8_t unit = *p;
        size_t shift = 0;

        for (size_t i = 0; i < bases_count; i++) {
            if (shift == 8) {
                shift = 0;
                p++;
                unit = *p;
            }

            code = (unit >> shift) & mask;
            KC__kmer_processor_handle_code(kp, i, code);

            shift += 2;
        }

        p++;
    }

    KC__ASSERT((char*)p - (char*)(buffer->data) == buffer->length);
}

void KC__kmer_processor_handle_buffer(KC__KmerProcessor* kp, const KC__Buffer* buffer) {
    KC__ASSERT(buffer != NULL);

    if (buffer->length == 0) {
        return;
    }

    switch (buffer->type) {
        case KC__BUFFER_TYPE_FASTA:
        case KC__BUFFER_TYPE_FASTQ:
            KC__kmer_processor_handle_reads_buffer(kp, buffer);
            break;
        case KC__BUFFER_TYPE_SUPER_KMER:
            KC__kmer_processor_handle_super_kmers_buffer(kp, buffer);
            break;
        default:
            KC__ASSERT(false);
            break;
    }
}

/**
 * Handle read until it reaches the end of the read or unexpected characters (other than A, a, C, c, G, g, T, t).
 * @param kp K-mer processor.
 * @param read The read to be handled.
 * @param read_length The max count of characters can be handled.
 * @param handled_length The count of characters already been handled.
 */
static inline void KC__kmer_processor_handle_sub_read(KC__KmerProcessor* kp, const char* read, size_t read_length, size_t* handled_length) {
    KC__KmerExtractUnit* keu = &(kp->kmer_extract_unit);
    if (read_length < keu->K)
        return;

    size_t i = 0;
    size_t skipped_count = 0;
    while (true) {
        if (i == read_length)
            break;

        KC__unit_t code = KC__encode(read[i]);

        if (code == 0x5) {
            // Unexpected character.
            break;

        } else if (code == 0x4) {
            // Skipped character.
            skipped_count++;

        } else {
            // Valid character.
            size_t j = i - skipped_count;
            KC__kmer_processor_handle_code(kp, j, code);
        }

        i++;
    }

    *handled_length = i;
}

void KC__kmer_processor_handle_read(struct KC__KmerProcessor* kp, const char* read, size_t read_length) {
    while (true) {
        size_t handled_length = 0;
        KC__kmer_processor_handle_sub_read(kp, read, read_length, &handled_length);
        read += handled_length;
        read_length -= handled_length;

        if (read_length == 0)
            break;
        read++;
        read_length--;
    }
}

static inline void KC__kmer_store_unit_set_action(KC__KmerStoreUnit* ksu, KC__KmerStoreAction action) {
    ksu->store_action = action;
}

static inline void* KC__kmer_store_unit_mem_request(KC__KmerStoreUnit* ksu, size_t size) {
    KC__Buffer* buffer = ksu->current_buffer;
    char* mem = (char*)(buffer->data) + buffer->length;
    buffer->length += size;
    return mem;
}

static inline bool KC__kmer_store_unit_mem_sufficient(KC__KmerStoreUnit* ksu) {
    KC__Buffer* buffer = ksu->current_buffer;
    return (buffer->size - buffer->length >= ksu->super_kmer_info_max_size);
}

static inline void KC__kmer_store_unit_expand(KC__KmerStoreUnit* ksu, KC__unit_t code) {
    if ((ksu->current_unit != NULL) && (ksu->current_bases_count == 4)) {
        ksu->current_unit = NULL;
    }

    if (ksu->current_unit == NULL) {
        ksu->current_unit = (uint8_t*)KC__kmer_store_unit_mem_request(ksu, sizeof(uint8_t));
        *(ksu->current_unit) = 0;
        ksu->current_bases_count = 0;
    }

    *(ksu->current_unit) |= (code << (ksu->current_bases_count << 1));
    ksu->current_bases_count += 1;
}

static inline void KC__kmer_processor_store_buffer_request(KC__KmerProcessor* kp, KC__Buffer** buffer, KC__BufferType buffer_type) {
    KC__Buffer* bf;
    if (kp->write_buffer_queue != NULL) {
        bf = KC__buffer_queue_get_blank_buffer(kp->write_buffer_queue);
    } else {
        bf = kp->store_buffer_request_callback();
    }

    KC__ASSERT(bf != NULL);

    bf->type = buffer_type;
    *buffer = bf;
}

static inline void KC__kmer_processor_store_buffer_complete(KC__KmerProcessor* kp, KC__Buffer** buffer) {
    KC__Buffer* bf = *buffer;
    KC__ASSERT(bf != NULL);

    if (kp->write_buffer_queue != NULL) {
        KC__buffer_queue_enqueue_filled_buffer(kp->write_buffer_queue, bf);
    } else {
        kp->store_buffer_complete_callback(bf);
    }

    *buffer = NULL;
}

static inline void KC__kmer_processor_store_kmer(KC__KmerProcessor* kp, KC__unit_t last_code) {
    KC__KmerStoreUnit* ksu = &(kp->kmer_store_unit);
    KC__KmerExtractUnit* keu = &(kp->kmer_extract_unit);

    if (ksu->store_action == KC__KMER_STORE_ACTION_NEW) {
        if ((ksu->current_buffer != NULL) && (!KC__kmer_store_unit_mem_sufficient(ksu))) {
            KC__kmer_processor_store_buffer_complete(kp, &(ksu->current_buffer));
        }

        if (ksu->current_buffer == NULL) {
            KC__kmer_processor_store_buffer_request(kp, &(ksu->current_buffer), KC__BUFFER_TYPE_SUPER_KMER);
            KC__ASSERT(KC__kmer_store_unit_mem_sufficient(ksu));

            ksu->super_kmers_count = (uint32_t*)KC__kmer_store_unit_mem_request(ksu, sizeof(uint32_t));
            *(ksu->super_kmers_count) = 0;
        }

        *(ksu->super_kmers_count) += 1;
        ksu->expanded_bases_count = (uint8_t*)KC__kmer_store_unit_mem_request(ksu, sizeof(uint8_t));
        *(ksu->expanded_bases_count) = 0;
        ksu->current_unit = NULL;

        size_t w = keu->gen_w_init;
        size_t s = keu->gen_s_init;
        for (size_t i = 0; i < keu->K; i++) {
            KC__unit_t code = ((keu->kmer[w] >> s) & 0x3);
            KC__kmer_store_unit_expand(ksu, code);

            if (s == 0) {
                w--;
                s = KC__UNIT_BIT - 2;
            } else {
                s -= 2;
            }
        }

        KC__kmer_store_unit_set_action(ksu, KC__KMER_STORE_ACTION_EXPAND);


    } else if (ksu->store_action == KC__KMER_STORE_ACTION_EXPAND) {
        KC__kmer_store_unit_expand(ksu, last_code);

        *(ksu->expanded_bases_count) += 1;
        if (*(ksu->expanded_bases_count) == UINT8_MAX) {
            KC__kmer_store_unit_set_action(ksu, KC__KMER_STORE_ACTION_NEW);
        }

    } else {
        KC__ASSERT(false);
    }
}

void KC__kmer_processor_handle_kmer(KC__KmerProcessor* kp, const KC__unit_t* canonical_kmer, size_t n, KC__unit_t last_code) {
    KC__KmerStoreUnit* ksu = &(kp->kmer_store_unit);

    if (n == 0) {
        KC__kmer_store_unit_set_action(ksu, KC__KMER_STORE_ACTION_NEW);
    }

    if (KC__hash_map_add_kmer(kp->hash_map, kp->id, canonical_kmer)) {
        KC__kmer_store_unit_set_action(ksu, KC__KMER_STORE_ACTION_NEW);
        return;
    }

    KC__kmer_processor_store_kmer(kp, last_code);
}


static void KC__kmer_processor_export_kmers_callback(const KC__unit_t* kmer, KC__count_t count, void* data) {
    KC__KmerProcessor* kp = data;
    KC__KmerExportUnit* ktu = &(kp->kmer_export_unit);

    ktu->total_kmers_count += count;
    ktu->unique_kmers_count++;

    KC__OutputParam* p = &(ktu->output_param);

    if ((count < p->filter_min) || (count > p->filter_max)) {
        return;
    }
    if (count > p->count_max) {
        count = p->count_max;
    }

    ktu->exported_unique_kmers_count++;

    if (ktu->buffer == NULL) {
        KC__kmer_processor_store_buffer_request(kp, &(ktu->buffer), KC__BUFFER_TYPE_KMER);
    }

    char* current_location = (char*)(ktu->buffer->data) + ktu->buffer->length;

    for (size_t i = 0; i < ktu->W; i++) {
        size_t c = (i == ktu->W - 1) ? (ktu->C) : (KC__UNIT_BIT / 8);
        for (size_t j = 0; j < c; j++) {
            uint8_t* mem = (uint8_t*)current_location;
            *mem = (uint8_t)((kmer[i] >> (j * 8)) & 0xFF);
            current_location += sizeof(uint8_t);
        }
    }
    switch (ktu->count_bit) {
        case 8:
            *((uint8_t*)current_location) = (uint8_t)count;
            break;
        case 16:
            *((uint16_t*)current_location) = (uint16_t)count;
            break;
        case 32:
            *((uint32_t*)current_location) = (uint32_t)count;
            break;
        default:
            *((KC__count_t*)current_location) = (KC__count_t)count;
            break;
    }

    ktu->buffer->length += ktu->unit_size;

    if (ktu->unit_size > (ktu->buffer->size - ktu->buffer->length)) {
        KC__kmer_processor_store_buffer_complete(kp, &(ktu->buffer));
    }
}

void KC__kmer_processor_export_kmers(KC__KmerProcessor* kp) {
    KC__KmerExportUnit* ktu = &(kp->kmer_export_unit);

    size_t K = kp->kmer_extract_unit.K;
    size_t kmer_width = KC__calculate_kmer_width(K);
    size_t kmer_width_by_8 = KC__calculate_kmer_width_by_unit_size(K, sizeof(uint8_t));

    ktu->W = kmer_width;
    ktu->C = kmer_width_by_8 - (kmer_width - 1) * (KC__UNIT_BIT / 8);
    ktu->buffer = NULL;

    size_t count_bit;
    size_t count_size;
    KC__calculate_count_field(ktu->output_param.count_max, &count_bit, &count_size);

    ktu->count_bit = count_bit;
    ktu->unit_size = sizeof(uint8_t) * kmer_width_by_8 + count_size;

    ktu->total_kmers_count = 0;
    ktu->unique_kmers_count = 0;
    ktu->exported_unique_kmers_count = 0;

    KC__hash_map_export(kp->hash_map, kp->id, KC__kmer_processor_export_kmers_callback, kp, NULL);

    if (ktu->buffer != NULL) {
        KC__kmer_processor_store_buffer_complete(kp, &(ktu->buffer));
    }
}

void KC__kmer_processor_get_exported_kmers_stats(KC__KmerProcessor* kp, size_t* total_kmers_count, size_t* unique_kmers_count, size_t* exported_unique_kmers_count) {
    *total_kmers_count = kp->kmer_export_unit.total_kmers_count;
    *unique_kmers_count = kp->kmer_export_unit.unique_kmers_count;
    *exported_unique_kmers_count = kp->kmer_export_unit.exported_unique_kmers_count;
}


void KC__kmer_processor_finish(KC__KmerProcessor* kp) {
    KC__KmerStoreUnit* ksu = &(kp->kmer_store_unit);
    if (ksu->current_buffer != NULL) {
        KC__kmer_processor_store_buffer_complete(kp, &(ksu->current_buffer));
    }

    KC__hash_map_finish_adding_kmers(kp->hash_map, kp->id);
}

void* KC__kmer_processor_work_extract(void* ptr) {
    KC__KmerProcessor* kp = ptr;

    while (true) {
        KC__Buffer *buffer = KC__buffer_queue_dequeue_filled_buffer(kp->read_buffer_queue);
        if (buffer == NULL) {
            break;
        }

        KC__kmer_processor_handle_buffer(kp, buffer);
        KC__buffer_queue_recycle_blank_buffer(kp->read_buffer_queue, buffer);
    }

    KC__kmer_processor_finish(kp);

    pthread_exit(NULL);
}

void* KC__kmer_processor_work_export(void* ptr) {
    KC__KmerProcessor* kp = ptr;
    KC__kmer_processor_export_kmers(kp);
    pthread_exit(NULL);
}