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

#include "dump.h"
#include "utils.h"
#include "logging.h"
#include "types.h"
#include "header.h"
#include "assert.h"


typedef struct {
    char* result_file_name;
    char* dump_file_name;
} KC__DumpParam;


static error_t parse_dump_opt(int key, char* arg, struct argp_state* state) {
    KC__DumpParam *param = state->input;

    switch (key) {
        case 'o':
            param->dump_file_name = arg;
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

void KC__decode_kmers(char* mem_block, char* kmer_sequence, uint64_t k, size_t kmer_width) {
    size_t first_invalid_count = 4 - (k - (kmer_width - 1) * 4);

    size_t kmer_sequence_idx = 0;

    for (size_t i = 0; i < kmer_width; i++) {
        char* p = mem_block + (kmer_width - i - 1);
        uint8_t code = *((uint8_t*)p);

        for (size_t n = 0; n < 4; n++) {
            if (i == 0 && n < first_invalid_count) {
                continue;
            }
            uint8_t nt_code = (code >> (4 - n - 1) * 2) & 0x3;
            char nt;
            switch (nt_code) {
                case 0x0:
                    nt = 'A';
                    break;
                case 0x1:
                    nt = 'C';
                    break;
                case 0x2:
                    nt = 'G';
                    break;
                case 0x3:
                    nt = 'T';
                    break;
                default:
                    LOGGING_DEBUG("Unexpected nt code: %d", nt_code);
                    KC__ASSERT(false);
                    break;
            }
            kmer_sequence[kmer_sequence_idx] = nt;
            kmer_sequence_idx += 1;
        }
    }
}

void KC__dump(int argc, char** argv) {
    KC__DumpParam param;

    param.dump_file_name = "./KC__dump.txt";

    struct argp_option options[] = {
            {"out", 'o', "OUT", 0, "Output dump file path"},
            {0}
    };
    struct argp argp = {options, parse_dump_opt, "RESULT", "Dump the k-mers counting result."};
    argp_parse(&argp, argc, argv, 0, 0, &param);

    LOGGING_DEBUG("Result file name: %s", param.result_file_name);
    LOGGING_DEBUG("Dump file name: %s", param.dump_file_name);

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

    FILE* wfp = fopen(param.dump_file_name, "w");
    if (wfp == NULL) {
        KC__file_error_exit(param.dump_file_name, "Open", NULL);
    }

    size_t count_bit;
    size_t count_size;
    KC__calculate_count_field(header.count_max, &count_bit, &count_size);

    const size_t kmer_width = KC__calculate_kmer_width_by_unit_size(header.K, sizeof(uint8_t));
    const size_t kmer_size = kmer_width * sizeof(uint8_t);
    const size_t mem_block_size_lim = 5000000;

    void* mem_block = malloc(mem_block_size_lim);
    const size_t kmer_info_size = kmer_size + count_size;
    const size_t mem_block_size = mem_block_size_lim / kmer_info_size * kmer_info_size;

    char* kmer_sequence = malloc(header.K + 1);
    kmer_sequence[header.K] = '\0';

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
        char *p = mem_block;
        for (size_t i = 0; i < read_kmer_info_count; i++) {

            KC__decode_kmers(p, kmer_sequence, header.K, kmer_width);

            p += kmer_size;
            size_t count;

            switch (count_bit) {
                case 8:
                    count = *((uint8_t *) p);
                    break;
                case 16:
                    count = *((uint16_t *) p);
                    break;
                case 32:
                    count = *((uint32_t *) p);
                    break;
                default:
                    count = *((KC__count_t *) p);
                    break;
            }

            fprintf(wfp, "%s\t%zu\n", kmer_sequence, count);

            p += count_size;
        }
    }

    fclose(fp);
    fclose(wfp);

    free(kmer_sequence);
    free(mem_block);
}