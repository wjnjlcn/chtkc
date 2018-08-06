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


#ifndef KC__KMER_PROCESSOR_H
#define KC__KMER_PROCESSOR_H

#include "types.h"
#include "param.h"
#include "mem_allocator.h"
#include "buffer_queue.h"
#include "hash_map.h"


struct KC__KmerProcessor;
typedef struct KC__KmerProcessor KC__KmerProcessor;

typedef void (*KC__KmerProcessorReadCallback) (KC__KmerProcessor* kmer_processor, const char* read, size_t read_length);
typedef void (*KC__KmerProcessorKmerCallback) (KC__KmerProcessor* kmer_processor, const KC__unit_t* canonical_kmer, size_t n, KC__unit_t last_code);
typedef KC__Buffer* (*KC__KmerProcessorStoreBufferRequestCallback) (void);
typedef void (*KC__KmerProcessorStoreBufferCompleteCallback) (KC__Buffer* buffer);

KC__KmerProcessor* KC__kmer_processor_create(KC__MemAllocator* mem_allocator, size_t id, size_t K, KC__OutputParam output_param);
void KC__kmer_processor_free(KC__MemAllocator* mem_allocator, KC__KmerProcessor* kmer_processor);
void KC__kmer_processor_link_modules(KC__KmerProcessor* kmer_processor, KC__HashMap* hash_map, KC__BufferQueue* read_buffer_queue, KC__BufferQueue* write_buffer_queue);

void KC__kmer_processor_set_read_callback(KC__KmerProcessor* kmer_processor, KC__KmerProcessorReadCallback read_callback);
void KC__kmer_processor_set_kmer_callback(KC__KmerProcessor* kmer_processor, KC__KmerProcessorKmerCallback kmer_callback);
void KC__kmer_processor_set_store_buffer_request_callback(KC__KmerProcessor* kmer_processor, KC__KmerProcessorStoreBufferRequestCallback request_callback);
void KC__kmer_processor_set_store_buffer_complete_callback(KC__KmerProcessor* kmer_processor, KC__KmerProcessorStoreBufferCompleteCallback complete_callback);

void KC__kmer_processor_handle_buffer(KC__KmerProcessor* kmer_processor, const KC__Buffer* buffer);
void KC__kmer_processor_handle_read(KC__KmerProcessor* kmer_processor, const char* read, size_t read_length);
void KC__kmer_processor_handle_kmer(KC__KmerProcessor* kmer_processor, const KC__unit_t* kmer, size_t n, KC__unit_t last_code);
void KC__kmer_processor_finish(KC__KmerProcessor* kmer_processor);

void KC__kmer_processor_export_kmers(KC__KmerProcessor* kmer_processor);
void KC__kmer_processor_get_exported_kmers_stats(KC__KmerProcessor* kmer_processor, size_t* total_kmers_count, size_t* unique_kmers_count, size_t* exported_unique_kmers_count);

void* KC__kmer_processor_work_extract(void* ptr);
void* KC__kmer_processor_work_export(void* ptr);

#endif
