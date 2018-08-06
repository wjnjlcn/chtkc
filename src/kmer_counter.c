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
#include <string.h>

#include "kmer_counter.h"
#include "file_reader.h"
#include "file_writer.h"
#include "kmer_processor.h"
#include "logging.h"
#include "param.h"
#include "header.h"


struct KC__KmerCounter {
    KC__Param* param;

    KC__FileReader** file_readers;
    size_t file_readers_count;

    KC__FileWriter* file_writer;

    KC__KmerProcessor** kmer_processors;
    size_t kmer_processors_count;

    KC__BufferQueue* read_buffer_queue;
    KC__BufferQueue* write_buffer_queue;

    KC__HashMap* hash_map;
};

KC__KmerCounter* KC__kmer_counter_create(KC__MemAllocator* ma, KC__Param* param) {
    KC__KmerCounter* kc = (KC__KmerCounter*)KC__mem_alloc(ma, sizeof(KC__KmerCounter), "kmer counter");
    kc->param = param;

    kc->file_readers_count = param->reading_threads_count;
    kc->file_readers = (KC__FileReader**)KC__mem_alloc(ma, sizeof(KC__FileReader*) * kc->file_readers_count, "kmer counter file readers");
    for (size_t i= 0; i < kc->file_readers_count; i++) {
        kc->file_readers[i] = KC__file_reader_create(ma, param->K, param->input_compression_type, param->read_buffer_size);
    }

    KC__Header header;
    header.K = param->K;
    header.count_max = param->output_param.count_max;
    header.filter_min = param->output_param.filter_min;
    header.filter_max = param->output_param.filter_max;

    kc->file_writer = KC__file_writer_create(ma, param->output_file_name, &header);

    kc->kmer_processors_count = param->kmer_processing_threads_count;
    kc->kmer_processors = (KC__KmerProcessor**)KC__mem_alloc(ma, sizeof(KC__KmerProcessor*) * kc->kmer_processors_count, "kmer counter kmer processors");
    for (size_t i = 0; i < kc->kmer_processors_count; i++) {
        kc->kmer_processors[i] = KC__kmer_processor_create(ma, i, param->K, param->output_param);
    }

    kc->read_buffer_queue = KC__buffer_queue_create(ma, param->read_buffer_size, param->read_buffers_count);
    kc->write_buffer_queue = KC__buffer_queue_create(ma, param->write_buffer_size, param->write_buffers_count);

    kc->hash_map = KC__hash_map_create(ma, param->K, kc->kmer_processors_count);

    for (size_t i= 0; i < kc->file_readers_count; i++) {
        KC__file_reader_link_modules(kc->file_readers[i], kc->read_buffer_queue);
    }
    KC__file_writer_link_modules(kc->file_writer, kc->write_buffer_queue);
    for (size_t i = 0; i < kc->kmer_processors_count; i++) {
        KC__kmer_processor_link_modules(kc->kmer_processors[i], kc->hash_map, kc->read_buffer_queue, kc->write_buffer_queue);
    }

    return kc;
}

void KC__kmer_counter_free(KC__MemAllocator* ma, KC__KmerCounter* kc) {
    for (size_t i= 0; i < kc->file_readers_count; i++) {
        KC__file_reader_free(ma, kc->file_readers[i]);
    }
    KC__mem_free(ma, kc->file_readers);
    KC__file_writer_free(ma, kc->file_writer);

    for (size_t i = 0; i < kc->kmer_processors_count; i++) {
        KC__kmer_processor_free(ma, kc->kmer_processors[i]);
    }
    KC__mem_free(ma, kc->kmer_processors);

    KC__buffer_queue_free(ma, kc->read_buffer_queue);
    KC__buffer_queue_free(ma, kc->write_buffer_queue);

    KC__hash_map_free(ma, kc->hash_map);

    KC__mem_free(ma, kc);
}

static inline void KC__kmer_counter_schedule_files(KC__FileInputDescription inputs[], size_t n, KC__Param* param) {
    size_t files_count_for_each = param->input_files_count / n;
    size_t remain_files_count = param->input_files_count % n;
    for (size_t i = 0; i < n; i++) {
        inputs[i].files_count = files_count_for_each;
    }
    for (size_t i = 0; i < remain_files_count; i++) {
        inputs[i].files_count += 1;
    }

    size_t offset = 0;
    for (size_t i = 0; i < n; i++) {
        inputs[i].file_names = &(param->input_file_names[offset]);
        offset += inputs[i].files_count;

        inputs[i].file_type = param->input_file_type;
        inputs[i].compression_type = param->input_compression_type;
    }
}

void KC__kmer_counter_work(KC__KmerCounter* kc) {
    KC__Param* param = kc->param;
    size_t inputs_count = kc->file_readers_count;

    size_t n = 0;

    KC__FileInputDescription inputs[inputs_count];
    KC__kmer_counter_schedule_files(inputs, inputs_count, param);

    size_t tmp_file_name_str_len = strlen(param->output_file_name) + strlen("_tmp_N") + 1;
    char tmp_file_1_name[tmp_file_name_str_len];
    strcpy(tmp_file_1_name, param->output_file_name);
    strcat(tmp_file_1_name, "_tmp_0");
    char tmp_file_2_name[tmp_file_name_str_len];
    strcpy(tmp_file_2_name, param->output_file_name);
    strcat(tmp_file_2_name, "_tmp_1");
    char* tmp_file_names[2] = {tmp_file_1_name, tmp_file_2_name};
    bool should_delete_tmp_files[2] = {true, false};
    int tmp_file_idx = 0;

    size_t total_kmers_count = 0;
    size_t unique_kmers_count = 0;
    size_t exported_unique_kmers_count = 0;


    while (true) {
        n++;

        LOGGING_INFO("Pass #%zu start.", n);

        pthread_t read_threads[inputs_count];
        pthread_t write_thread;
        pthread_t process_threads[kc->kmer_processors_count];

        // Inform buffer queues to start input.
        KC__buffer_queue_start_input(kc->read_buffer_queue);
        KC__buffer_queue_start_input(kc->write_buffer_queue);

        // Start reading thread.
        for (size_t i = 0; i < inputs_count; i++) {
            KC__file_reader_update_input(kc->file_readers[i], inputs[i]);
            pthread_create(&(read_threads[i]), NULL, KC__file_reader_work, kc->file_readers[i]);
        }

        // Start extracting threads.
        for (size_t i = 0; i < kc->kmer_processors_count; i++) {
            pthread_create(&(process_threads[i]), NULL, KC__kmer_processor_work_extract, kc->kmer_processors[i]);
        }

        // Start writing thread.
        const char* tmp_file_name = tmp_file_names[tmp_file_idx];
        KC__file_writer_update_tmp_file(kc->file_writer, tmp_file_name);
        pthread_create(&write_thread, NULL, KC__file_writer_work, kc->file_writer);

        // Reading thread finished, read buffer queue input finished.
        for (size_t i = 0; i < inputs_count; i++) {
            pthread_join(read_threads[i], NULL);
        }
        KC__buffer_queue_finish_input(kc->read_buffer_queue);

        // Extracting threads finished.
        for (size_t i = 0; i < kc->kmer_processors_count; i++) {
            pthread_join(process_threads[i], NULL);
        }

        // Start exporting threads.
        for (size_t i = 0; i < kc->kmer_processors_count; i++) {
            pthread_create(&(process_threads[i]), NULL, KC__kmer_processor_work_export, kc->kmer_processors[i]);
        }

        // Exporting threads finished, write buffer queue input finished.
        for (size_t i = 0; i < kc->kmer_processors_count; i++) {
            pthread_join(process_threads[i], NULL);
        }
        KC__buffer_queue_finish_input(kc->write_buffer_queue);

        // Writing thread finished.
        pthread_join(write_thread, NULL);


        for (size_t i = 0; i < kc->kmer_processors_count; i++) {
            size_t tc;
            size_t uc;
            size_t euc;
            KC__kmer_processor_get_exported_kmers_stats(kc->kmer_processors[i], &tc, &uc, &euc);
            total_kmers_count += tc;
            unique_kmers_count += uc;
            exported_unique_kmers_count += euc;
        }


        size_t tmp_file_size = KC__file_writer_get_tmp_file_size(kc->file_writer);
        LOGGING_DEBUG("Tmp file size: %zu", tmp_file_size);

        if (tmp_file_size == 0) {
            break;
        }

        inputs_count = 1;
        inputs[0].file_names = &(tmp_file_names[tmp_file_idx]);
        inputs[0].files_count = 1;
        inputs[0].file_type = KC__FILE_TYPE_SUPER_KMER;
        inputs[0].compression_type = KC__FILE_COMPRESSION_TYPE_PLAIN;

        tmp_file_idx = (tmp_file_idx + 1) % 2;
        should_delete_tmp_files[tmp_file_idx] = true;

        KC__hash_map_clear(kc->hash_map);
    }

    for (size_t i = 0; i < 2; i++) {
        if (should_delete_tmp_files[i]) {
            if (remove(tmp_file_names[i]) != 0) {
                LOGGING_WARNING("Delete file failed: %s", tmp_file_names[i]);
            }
        }
    }

    LOGGING_INFO("Total K-mers count: %zu", total_kmers_count);
    LOGGING_INFO("Unique K-mers count: %zu", unique_kmers_count);
    LOGGING_INFO("Exported unique K-mers count: %zu", exported_unique_kmers_count);
}