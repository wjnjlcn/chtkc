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


#include <stdio.h>
#include <stdlib.h>
#include <argp.h>

#include "histo.h"
#include "utils.h"
#include "logging.h"
#include "types.h"
#include "header.h"


typedef struct {
    char* result_file_name;
    char* histo_file_name;
} KC__HistoParam;

typedef struct {
    size_t key;
    size_t value;
} KC__HistoItem;

static int KC__histo_items_compare(const void* item_1, const void* item_2) {
    size_t key_1 = ((KC__HistoItem*)item_1)->key;
    size_t key_2 = ((KC__HistoItem*)item_2)->key;
    if (key_1 < key_2)
        return -1;
    if (key_1 > key_2)
        return 1;
    return 0;
}

static error_t parse_histo_opt(int key, char* arg, struct argp_state* state) {
    KC__HistoParam *param = state->input;

    switch (key) {
        case 'o':
            param->histo_file_name = arg;
            break;
        case ARGP_KEY_ARG:
            if (state->arg_num == 0) {
                param->result_file_name = arg;
            }
            break;
        case ARGP_KEY_END:
            if (state->arg_num != 1)
                argp_usage(state);
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

void KC__histo(int argc, char** argv) {
    KC__HistoParam param;

    param.histo_file_name = "./KC__histo.txt";

    struct argp_option options[] = {
            {"out", 'o', "OUT", 0, "Output histo file path"},
            {0}
    };
    struct argp argp = {options, parse_histo_opt, "RESULT", "Generate histogram for k-mers."};
    argp_parse(&argp, argc, argv, 0, 0, &param);

    LOGGING_DEBUG("Result file name: %s", param.result_file_name);
    LOGGING_DEBUG("Histo file name: %s", param.histo_file_name);


    const char* file_name = param.result_file_name;
    FILE* fp = fopen(file_name, "rb");
    if (fp == NULL) {
        KC__file_error_exit(file_name, "Open", NULL);
    }

    KC__Header header;
    bool success = KC__read_header(&header, fp);
    if (!success) {
        KC__file_error_exit(file_name, "Read header", NULL);
    }

    LOGGING_DEBUG("K: %zu, count max: %zu, filter min: %zu, max: %zu", header.K, header.count_max, header.filter_min, header.filter_max);

    size_t count_bit;
    size_t count_size;
    KC__calculate_count_field(header.count_max, &count_bit, &count_size);

    const size_t kmer_size = KC__calculate_kmer_width_by_unit_size(header.K, sizeof(uint8_t)) * sizeof(uint8_t);
    const size_t mem_block_size_lim = 5000000;
    const size_t counts_array_length = 100000;

    void* mem_block = malloc(mem_block_size_lim);


    size_t* counts_array = (size_t*)malloc(sizeof(size_t) * counts_array_length);
    for (size_t i = 0; i < counts_array_length; i++) {
        counts_array[i] = 0;
    }

    size_t histo_items_array_capacity = 100;
    size_t histo_items_array_length = 0;
    KC__HistoItem* histo_items_array = (KC__HistoItem*)malloc(sizeof(KC__HistoItem) * histo_items_array_capacity);


    const size_t kmer_info_size = kmer_size + count_size;
    const size_t mem_block_size = mem_block_size_lim / kmer_info_size * kmer_info_size;

    while (true) {
        const size_t read_size = fread(mem_block, 1, mem_block_size, fp);

        if (ferror(fp)) {
            KC__file_error_exit(file_name, "Read", NULL);
        }

        if (read_size == 0) {
            break;
        }

        if (read_size % kmer_info_size != 0) {
            KC__file_error_exit(file_name, "Parse", "file is truncated");
        }

        const size_t read_kmer_info_count = read_size / kmer_info_size;
        char* p = mem_block;
        for (size_t i = 0; i < read_kmer_info_count; i++) {
            p += kmer_size;
            size_t count;

            switch (count_bit) {
                case 8:
                    count = *((uint8_t*)p);
                    break;
                case 16:
                    count = *((uint16_t*)p);
                    break;
                case 32:
                    count = *((uint32_t*)p);
                    break;
                default:
                    count = *((KC__count_t*)p);
                    break;
            }
            p += count_size;

            if (count < counts_array_length) {
                counts_array[count]++;
            } else {
                bool found = false;
                for (size_t n = 0; n < histo_items_array_length; n++) {
                    if (count == histo_items_array[n].key) {
                        histo_items_array[n].value++;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    if (histo_items_array_length == histo_items_array_capacity) {
                        histo_items_array_capacity *= 2;
                        histo_items_array = (KC__HistoItem*)realloc(histo_items_array, sizeof(KC__HistoItem) * histo_items_array_capacity);
                        LOGGING_DEBUG("Histo items array expanded capacity to %zu", histo_items_array_capacity);
                    }
                    histo_items_array[histo_items_array_length].key = count;
                    histo_items_array[histo_items_array_length].value = 1;
                    histo_items_array_length++;
                }
            }
        }
    }

    fclose(fp);

    qsort(histo_items_array, histo_items_array_length, sizeof(KC__HistoItem), KC__histo_items_compare);

    size_t total_kmers_count = 0;
    size_t unique_kmers_count = 0;

    fp = fopen(param.histo_file_name, "w");
    if (fp == NULL) {
        KC__file_error_exit(param.histo_file_name, "Open", NULL);
    }

    const char* fmt = "%zu\t%zu\n";
    for (size_t i = 0; i < counts_array_length; i++) {
        size_t c = counts_array[i];
        if (c != 0) {
            total_kmers_count += i * c;
            unique_kmers_count += c;
            fprintf(fp, fmt, i, c);
        }
    }

    for (size_t i = 0; i < histo_items_array_length; i++) {
        size_t n = histo_items_array[i].key;
        size_t c = histo_items_array[i].value;
        total_kmers_count += n * c;
        unique_kmers_count += c;
        fprintf(fp, fmt, n, c);
    }

    fclose(fp);

    LOGGING_DEBUG("Total K-mers count: %zu", total_kmers_count);
    LOGGING_DEBUG("Unique K-mers count: %zu", unique_kmers_count);

    free(mem_block);
    free(counts_array);
    free(histo_items_array);
}