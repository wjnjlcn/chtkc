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


#ifndef KC__QUEUE_H
#define KC__QUEUE_H

#include <stdbool.h>
#include "mem_allocator.h"


struct KC__Queue;
typedef struct KC__Queue KC__Queue;

KC__Queue* KC__queue_create(KC__MemAllocator* mem_allocator, size_t capacity);
void KC__queue_free(KC__MemAllocator* mem_allocator, KC__Queue* queue);

/**
 * Available items count of queue.
 * @param queue Queue.
 * @return The length of queue.
 */
size_t KC__queue_length(const KC__Queue* queue);

/**
 * The capacity of queue.
 * @param queue Queue.
 * @return The capacity of queue.
 */
size_t KC__queue_capacity(const KC__Queue* queue);

bool KC__queue_is_empty(const KC__Queue* queue);
bool KC__queue_is_full(const KC__Queue* queue);

void* KC__queue_front(const KC__Queue* queue);

/**
 * Enqueue an item into queue.
 * @param queue Queue.
 * @param item The item to enqueue.
 * @return If success return true, else return false.
 */
bool KC__queue_enqueue(KC__Queue* queue, void* item);

/**
 * Dequeue an item from queue.
 * @param queue Queue.
 * @return If success return the item, else return NULL.
 */
void* KC__queue_dequeue(KC__Queue* queue);

#endif
