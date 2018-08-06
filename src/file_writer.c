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
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "file_writer.h"
#include "assert.h"
#include "logging.h"


struct KC__FileWriter {
    const char* output_file_name;
    FILE* output_file;

    const char* tmp_file_name;
    size_t tmp_file_size;

    KC__BufferQueue* buffer_queue;
};

KC__FileWriter* KC__file_writer_create(KC__MemAllocator* ma, const char* output_file_name, const KC__Header* header) {
    KC__FileWriter* fw = (KC__FileWriter*)KC__mem_alloc(ma, sizeof(KC__FileWriter), "file writer");

    fw->output_file_name = output_file_name;
    fw->output_file = fopen(fw->output_file_name, "wb");
    if (fw->output_file == NULL) {
        LOGGING_ERROR("Open output file error [%s]", fw->output_file_name);
        exit(EXIT_FAILURE);
    }

    if (header != NULL) {
        bool success = KC__write_header(header, fw->output_file);
        if (!success) {
            LOGGING_ERROR("Write header to file error [%s]", fw->output_file_name);
            exit(EXIT_FAILURE);
        }
    }

    fw->buffer_queue = NULL;

    return fw;
}

void KC__file_writer_free(KC__MemAllocator* ma, KC__FileWriter* fw) {
    KC__ASSERT(fw->output_file != NULL);
    fclose(fw->output_file);

    KC__mem_free(ma, fw);
}

void KC__file_writer_link_modules(KC__FileWriter* fw, KC__BufferQueue* buffer_queue) {
    fw->buffer_queue = buffer_queue;
}

void KC__file_writer_update_tmp_file(KC__FileWriter* fw, const char* tmp_file_name) {
    fw->tmp_file_name = tmp_file_name;
    fw->tmp_file_size = 0;
}

size_t KC__file_writer_get_tmp_file_size(const KC__FileWriter* fw) {
    return fw->tmp_file_size;
}

void* KC__file_writer_work(void* ptr) {
    KC__FileWriter* fw = ptr;

    FILE* tmp_file = NULL;

    if (fw->tmp_file_name != NULL) {
        tmp_file = fopen(fw->tmp_file_name, "wb");
        if (tmp_file == NULL) {
            LOGGING_ERROR("Open tmp file error [%s]", fw->tmp_file_name);
            exit(EXIT_FAILURE);
        }
    }

    while (true) {
        KC__Buffer* buffer = KC__buffer_queue_dequeue_filled_buffer(fw->buffer_queue);
        if (buffer == NULL) {
            break;
        }

        const char* file_name;
        FILE* file;
        bool write_buffer_length = false;

        switch (buffer->type) {
            case KC__BUFFER_TYPE_SUPER_KMER:
                file_name = fw->tmp_file_name;
                file = tmp_file;
                write_buffer_length = true;
                break;
            case KC__BUFFER_TYPE_KMER:
                file_name = fw->output_file_name;
                file = fw->output_file;
                break;
            default:
                KC__ASSERT(false);
                break;
        }

        if (write_buffer_length) {
            size_t write_size = fwrite(&(buffer->length), 1, sizeof(uint32_t), file);
            if (write_size < sizeof(uint32_t)) {
                LOGGING_ERROR("Write file error [%s]", file_name);
                exit(EXIT_FAILURE);
            }
        }

        size_t write_size = fwrite(buffer->data, 1, buffer->length, file);
        if (write_size < buffer->length) {
            LOGGING_ERROR("Write file error [%s]", file_name);
            exit(EXIT_FAILURE);
        }

        KC__buffer_queue_recycle_blank_buffer(fw->buffer_queue, buffer);
    }

    if (tmp_file != NULL) {
        long pos = ftell(tmp_file);
        if (pos >= 0) {
            fw->tmp_file_size = (size_t) pos;
        } else {
            LOGGING_ERROR("Getting tmp file size error [%s]", fw->tmp_file_name);
            exit(EXIT_FAILURE);
        }

        fclose(tmp_file);
    }

    pthread_exit(NULL);
}