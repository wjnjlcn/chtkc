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


#ifndef KC__BUFFER_QUEUE_H
#define KC__BUFFER_QUEUE_H

#include <stdint.h>
#include "mem_allocator.h"


typedef enum {
    KC__BUFFER_TYPE_FASTA = 0,
    KC__BUFFER_TYPE_FASTQ,
    KC__BUFFER_TYPE_SUPER_KMER,
    KC__BUFFER_TYPE_KMER
} KC__BufferType;


typedef struct {
    void* data;
    KC__BufferType type;
    uint32_t size;
    uint32_t length;
} KC__Buffer;


struct KC__BufferQueue;
typedef struct KC__BufferQueue KC__BufferQueue;

KC__BufferQueue* KC__buffer_queue_create(KC__MemAllocator* mem_allocator, uint32_t buffer_size, size_t buffers_count);
void KC__buffer_queue_free(KC__MemAllocator* mem_allocator, KC__BufferQueue* buffer_queue);

/**
 * Inform queue to accept producing,
 * should be called before producers and consumers start running.
 * @param buffer_queue Buffer queue.
 */
void KC__buffer_queue_start_input(KC__BufferQueue* buffer_queue);

/**
 * Inform queue producing is finished,
 * should be called after producers stop, before consumers stop.
 * @param buffer_queue Buffer queue.
 */
void KC__buffer_queue_finish_input(KC__BufferQueue* buffer_queue);

/**
 * Get a blank buffer to produce, should be called by producer.
 * @param buffer_queue Buffer queue.
 * @return A blank buffer, producer can always get one.
 */
KC__Buffer* KC__buffer_queue_get_blank_buffer(KC__BufferQueue* buffer_queue);

/**
 * Enqueue a filled buffer, should be called by producer, always success.
 * @param buffer_queue Buffer queue.
 * @param filled_buffer Buffer filled by producer.
 */
void KC__buffer_queue_enqueue_filled_buffer(KC__BufferQueue* buffer_queue, KC__Buffer* filled_buffer);

/**
 * Dequeue a filled buffer, should be called by consumer.
 * @param bufferQueue Buffer queue.
 * @return Buffer to consume, if no buffers left to be consumed, return NULL.
 */
KC__Buffer* KC__buffer_queue_dequeue_filled_buffer(KC__BufferQueue* bufferQueue);

/**
 * Recycle a blank buffer, should be called by consumer, always success.
 * @param buffer_queue Buffer queue.
 * @param blank_buffer Blank buffer already used by consumer.
 */
void KC__buffer_queue_recycle_blank_buffer(KC__BufferQueue* buffer_queue, KC__Buffer* blank_buffer);

#endif
