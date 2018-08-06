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


#include "queue.h"
#include "assert.h"


struct KC__Queue {
    void** array;
    size_t capacity;
    size_t length;
    size_t front;
    size_t rear;
};


KC__Queue* KC__queue_create(KC__MemAllocator* ma, size_t capacity) {
    KC__ASSERT(capacity > 0);

    KC__Queue* queue = (KC__Queue*)KC__mem_alloc(ma, sizeof(KC__Queue), "queue");

    queue->capacity = capacity;
    queue->array = (void**)KC__mem_alloc(ma, sizeof(void*) * queue->capacity, "queue array");

    queue->front = 0;
    queue->rear = 0;
    queue->length = 0;
    return queue;
}

void KC__queue_free(KC__MemAllocator* ma, KC__Queue* queue) {
    KC__mem_free(ma, queue->array);
    KC__mem_free(ma, queue);
}

size_t KC__queue_length(const KC__Queue* queue) {
    return queue->length;
}

size_t KC__queue_capacity(const KC__Queue* queue) {
    return queue->capacity;
}

bool KC__queue_is_empty(const KC__Queue* queue) {
    return (queue->length == 0);
}

bool KC__queue_is_full(const KC__Queue* queue) {
    return (queue->length == queue->capacity);
}

void* KC__queue_front(const KC__Queue* queue) {
    if (KC__queue_is_empty(queue))
        return NULL;
    return queue->array[queue->front];
}

bool KC__queue_enqueue(KC__Queue* queue, void* item) {
    if (KC__queue_is_full(queue))
        return false;

    queue->array[queue->rear] = item;
    queue->rear++;
    if (queue->rear == queue->capacity)
        queue->rear = 0;
    queue->length++;
    return true;
}

void* KC__queue_dequeue(KC__Queue* queue) {
    if (KC__queue_is_empty(queue))
        return NULL;

    void * front = queue->array[queue->front];
    queue->front++;
    if (queue->front == queue->capacity)
        queue->front = 0;
    queue->length--;
    return front;
}