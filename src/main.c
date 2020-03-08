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


#include <argp.h>
#include <string.h>
#include <time.h>
#include <stdio.h>

#include "kmer_counter.h"
#include "histo.h"
#include "dump.h"
#include "logging.h"
#include "assert.h"

const char* argp_program_version = "CHTKC 1.0.1";

static inline void KC__print_usage(const char* program_name) {
    printf("Usage: %s <CMD> [OPTION...] ARGS...\n"
           "  <CMD> is one of: count, histo, dump\n"
           "\n"
           "  -?, --help                 Give this help list\n"
           "  -V, --version              Print program version\n",
           program_name);
}

int main(int argc, char** argv) {
    KC__ASSERT(argp_program_version != NULL);

    if (argc < 1) {
        return 0;
    }
    const char* program_name = argv[0];

    KC__LOG_FILE = stderr;

    argc--;
    argv++;

    if (argc == 0) {
        KC__print_usage(program_name);
        return 0;
    }

    if (strcmp(argv[0], "count") == 0) {
        KC__Param param;
        KC__param_init(&param, argc, argv);

        time_t start_time = time(NULL);

        KC__MemAllocator *ma = KC__mem_allocator_create(param.mem_limit);

        KC__KmerCounter *kc = KC__kmer_counter_create(ma, &param);
        KC__kmer_counter_work(kc);
        KC__kmer_counter_free(ma, kc);

        KC__mem_allocator_free(ma);

        time_t end_time = time(NULL);
        LOGGING_INFO("Count running time: %zus", end_time - start_time);

        KC__param_destroy(&param);

    } else if (strcmp(argv[0], "histo") == 0) {
        KC__histo(argc, argv);

    } else if (strcmp(argv[0], "dump") == 0) {
        KC__dump(argc, argv);

    } else if ((strcmp(argv[0], "-?") == 0) || (strcmp(argv[0], "--help") == 0)) {
        KC__print_usage(program_name);

    } else if ((strcmp(argv[0], "-V") == 0) || (strcmp(argv[0], "--version") == 0)) {
        printf("%s\n", argp_program_version);

    } else {
        KC__print_usage(program_name);
    }

    return 0;
}