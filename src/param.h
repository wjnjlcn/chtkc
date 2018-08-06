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


#ifndef KC__PARAM_H
#define KC__PARAM_H

#include <stddef.h>
#include "types.h"


typedef struct {
    KC__count_t filter_min;
    KC__count_t filter_max;
    KC__count_t count_max;
} KC__OutputParam;

typedef struct {
    size_t K;

    size_t threads_count;
    size_t reading_threads_count;
    size_t kmer_processing_threads_count;

    char** input_file_names;
    size_t input_files_count;
    KC__FileType input_file_type;
    KC__FileCompressionType input_compression_type;

    const char* output_file_name;

    uint32_t read_buffer_size;
    size_t read_buffers_count;

    uint32_t write_buffer_size;
    size_t write_buffers_count;

    size_t mem_limit;

    KC__OutputParam output_param;

    const char* log_file_name;
} KC__Param;


void KC__param_init(KC__Param* param, int argc, char** argv);
void KC__param_destroy(KC__Param* param);

#endif
