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


#ifndef KC__FILE_READER_H
#define KC__FILE_READER_H

#include "types.h"
#include "mem_allocator.h"
#include "buffer_queue.h"

struct KC__FileReader;
typedef struct KC__FileReader KC__FileReader;

typedef struct {
    char** file_names;
    size_t files_count;
    KC__FileType file_type;
    KC__FileCompressionType compression_type;
} KC__FileInputDescription;

KC__FileReader* KC__file_reader_create(KC__MemAllocator* mem_allocator, size_t K, KC__FileCompressionType compression_type, size_t buffer_size);
void KC__file_reader_free(KC__MemAllocator* mem_allocator, KC__FileReader* file_reader);

void KC__file_reader_link_modules(KC__FileReader* file_reader, KC__BufferQueue* buffer_queue);
void KC__file_reader_update_input(KC__FileReader* file_reader, KC__FileInputDescription input);

void* KC__file_reader_work(void* ptr);

#endif
