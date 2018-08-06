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


#include <pthread.h>
#include "buffer_queue.h"
#include "queue.h"
#include "logging.h"
#include "assert.h"


struct KC__BufferQueue {
    KC__Buffer* buffers;
    size_t buffers_count;

    KC__Queue* blank_buffers_queue;
    KC__Queue* filled_buffers_queue;

    bool input_finished;

    pthread_mutex_t blank_buffers_queue_mtx;
    pthread_mutex_t filled_buffers_queue_mtx;
    pthread_cond_t cv_has_blank_buffers;
    pthread_cond_t cv_has_filled_buffers;
};


KC__BufferQueue* KC__buffer_queue_create(KC__MemAllocator* ma, uint32_t buffer_size, size_t buffers_count) {
    KC__ASSERT(buffer_size > 0);
    KC__ASSERT(buffers_count > 0);

    KC__BufferQueue* bq = (KC__BufferQueue*)KC__mem_alloc(ma, sizeof(KC__BufferQueue), "buffer queue");

    bq->buffers_count = buffers_count;
    LOGGING_DEBUG("buffers count: %zu", bq->buffers_count);

    bq->buffers = (KC__Buffer*)KC__mem_alloc(ma, sizeof(KC__Buffer) * bq->buffers_count, "buffer queue buffers");

    for (size_t i = 0; i < bq->buffers_count; i++) {
        bq->buffers[i].data = KC__mem_alloc(ma, buffer_size, "buffer queue buffers data");
        bq->buffers[i].size = buffer_size;
    }

    bq->blank_buffers_queue = KC__queue_create(ma, bq->buffers_count);
    bq->filled_buffers_queue = KC__queue_create(ma, bq->buffers_count);

    for (size_t i = 0; i < bq->buffers_count; i++)
        KC__queue_enqueue(bq->blank_buffers_queue, &(bq->buffers[i]));

    bq->input_finished = true;

    pthread_mutex_init(&(bq->blank_buffers_queue_mtx), NULL);
    pthread_mutex_init(&(bq->filled_buffers_queue_mtx), NULL);
    pthread_cond_init(&(bq->cv_has_blank_buffers), NULL);
    pthread_cond_init(&(bq->cv_has_filled_buffers), NULL);

    return bq;
}

void KC__buffer_queue_free(KC__MemAllocator* ma, KC__BufferQueue* bq) {
    KC__ASSERT(KC__queue_is_full(bq->blank_buffers_queue));
    KC__ASSERT(KC__queue_is_empty(bq->filled_buffers_queue));

    KC__queue_free(ma, bq->blank_buffers_queue);
    KC__queue_free(ma, bq->filled_buffers_queue);

    for (size_t i = 0; i < bq->buffers_count; i++) {
        KC__mem_free(ma, bq->buffers[i].data);
    }
    KC__mem_free(ma, bq->buffers);

    pthread_mutex_destroy(&(bq->blank_buffers_queue_mtx));
    pthread_mutex_destroy(&(bq->filled_buffers_queue_mtx));
    pthread_cond_destroy(&(bq->cv_has_blank_buffers));
    pthread_cond_destroy(&(bq->cv_has_filled_buffers));

    KC__mem_free(ma, bq);
}

void KC__buffer_queue_start_input(KC__BufferQueue* bq) {
    pthread_mutex_lock(&(bq->filled_buffers_queue_mtx));
    bq->input_finished = false;
    pthread_mutex_unlock(&(bq->filled_buffers_queue_mtx));
}

void KC__buffer_queue_finish_input(KC__BufferQueue* bq) {
    pthread_mutex_lock(&(bq->filled_buffers_queue_mtx));
    bq->input_finished = true;
    pthread_cond_broadcast(&(bq->cv_has_filled_buffers));
    pthread_mutex_unlock(&(bq->filled_buffers_queue_mtx));
}

KC__Buffer* KC__buffer_queue_get_blank_buffer(KC__BufferQueue* bq) {
    KC__Buffer* blank_buffer;

    pthread_mutex_lock(&(bq->blank_buffers_queue_mtx));
    while (true) {
        blank_buffer = KC__queue_dequeue(bq->blank_buffers_queue);
        if (blank_buffer != NULL)
            break;
        pthread_cond_wait(&(bq->cv_has_blank_buffers), &(bq->blank_buffers_queue_mtx));
    }
    pthread_mutex_unlock(&(bq->blank_buffers_queue_mtx));

    blank_buffer->length = 0;
    return blank_buffer;
}

void KC__buffer_queue_enqueue_filled_buffer(KC__BufferQueue* bq, KC__Buffer* filled_buffer) {
    KC__ASSERT(filled_buffer->length <= filled_buffer->size);

    pthread_mutex_lock(&(bq->filled_buffers_queue_mtx));
    bool success = KC__queue_enqueue(bq->filled_buffers_queue, filled_buffer);
    KC__ASSERT(success);
    pthread_cond_signal(&(bq->cv_has_filled_buffers));
    pthread_mutex_unlock(&(bq->filled_buffers_queue_mtx));
}

KC__Buffer* KC__buffer_queue_dequeue_filled_buffer(KC__BufferQueue* bq) {
    KC__Buffer* filled_buffer;

    pthread_mutex_lock(&(bq->filled_buffers_queue_mtx));
    while (true) {
        filled_buffer = KC__queue_dequeue(bq->filled_buffers_queue);
        if (filled_buffer != NULL)
            break;
        if (bq->input_finished)
            break;
        pthread_cond_wait(&(bq->cv_has_filled_buffers), &(bq->filled_buffers_queue_mtx));
    }
    pthread_mutex_unlock(&(bq->filled_buffers_queue_mtx));

    return filled_buffer;
}

void KC__buffer_queue_recycle_blank_buffer(KC__BufferQueue* bq, KC__Buffer* blank_buffer) {
    pthread_mutex_lock(&(bq->blank_buffers_queue_mtx));
    bool success = KC__queue_enqueue(bq->blank_buffers_queue, blank_buffer);
    KC__ASSERT(success);
    pthread_cond_signal(&(bq->cv_has_blank_buffers));
    pthread_mutex_unlock(&(bq->blank_buffers_queue_mtx));
}