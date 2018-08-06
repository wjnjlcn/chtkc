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


#ifndef KC__FILE_WRITER_H
#define KC__FILE_WRITER_H

#include "mem_allocator.h"
#include "buffer_queue.h"
#include "header.h"

struct KC__FileWriter;
typedef struct KC__FileWriter KC__FileWriter;

KC__FileWriter* KC__file_writer_create(KC__MemAllocator* mem_allocator, const char* output_file_name, const KC__Header* header);
void KC__file_writer_free(KC__MemAllocator* mem_allocator, KC__FileWriter* file_writer);

void KC__file_writer_link_modules(KC__FileWriter* file_writer, KC__BufferQueue* buffer_queue);
void KC__file_writer_update_tmp_file(KC__FileWriter* file_writer, const char* tmp_file_name);
size_t KC__file_writer_get_tmp_file_size(const KC__FileWriter* file_writer);

void* KC__file_writer_work(void* ptr);

#endif
