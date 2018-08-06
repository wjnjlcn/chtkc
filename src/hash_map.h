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


#ifndef KC__HASH_MAP_H
#define KC__HASH_MAP_H

#include <stdbool.h>
#include "types.h"
#include "mem_allocator.h"

struct KC__HashMap;
typedef struct KC__HashMap KC__HashMap;

typedef void (*KC__HashMapExportCallback) (const KC__unit_t* kmer, KC__count_t count, void* data);

KC__HashMap* KC__hash_map_create(KC__MemAllocator* mem_allocator, size_t K, size_t threads_count);
void KC__hash_map_free(KC__MemAllocator* mem_allocator, KC__HashMap* hash_map);

size_t KC__hash_map_max_key_count(const KC__HashMap* hash_map);
void KC__hash_map_set_table_capacity(KC__HashMap* hash_map, size_t capacity);
void KC__hash_map_lock_keys(KC__HashMap* hash_map);

void KC__hash_map_clear(KC__HashMap* hash_map);
bool KC__hash_map_add_kmer(KC__HashMap* hash_map, size_t thread_id, const KC__unit_t* kmer);
void KC__hash_map_finish_adding_kmers(KC__HashMap* hash_map, size_t thread_id);

void KC__hash_map_export(KC__HashMap* hash_map, size_t thread_id, KC__HashMapExportCallback callback, void* data, size_t* exported_count);

#endif
