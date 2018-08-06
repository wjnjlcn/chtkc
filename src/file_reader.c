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
#include <string.h>
#include <zlib.h>

#include "file_reader.h"
#include "logging.h"
#include "assert.h"
#include "buffer_queue.h"


struct KC__FileReader {
    KC__FileInputDescription input;

    const char* file_name;

    KC__BufferQueue* buffer_queue;

    size_t K;

    z_stream gz_stream;
    void* gz_data;
    size_t gz_data_size;
};

KC__FileReader* KC__file_reader_create(KC__MemAllocator* ma, size_t K, KC__FileCompressionType compressionType, size_t buffer_size) {
    KC__FileReader* fr = (KC__FileReader*)KC__mem_alloc(ma, sizeof(KC__FileReader), "file reader");

    fr->K = K;
    fr->buffer_queue = NULL;

    fr->gz_data_size = buffer_size;
    fr->gz_data = NULL;

    switch (compressionType) {
        case KC__FILE_COMPRESSION_TYPE_PLAIN:
            break;
        case KC__FILE_COMPRESSION_TYPE_GZIP:
            fr->gz_data = KC__mem_alloc(ma, fr->gz_data_size, "file reader gz data");
            break;
        default:
            KC__ASSERT(false);
            break;
    }

    return fr;
}

void KC__file_reader_free(KC__MemAllocator* ma, KC__FileReader* fr) {
    if (fr->gz_data) {
        KC__mem_free(ma, fr->gz_data);
    }
    KC__mem_free(ma, fr);
}

void KC__file_reader_link_modules(KC__FileReader* fr, KC__BufferQueue* buffer_queue) {
    fr->buffer_queue = buffer_queue;
}

void KC__file_reader_update_input(KC__FileReader* fr, KC__FileInputDescription input) {
    fr->input = input;
    fr->file_name = NULL;
}

typedef enum {
    KC__FILE_READ_ERROR_OPEN = 0,
    KC__FILE_READ_ERROR_READ,
    KC__FILE_READ_ERROR_PARSE
} KC__FileReadError;

static inline void KC__file_reader_process_file_error_exit(KC__FileReader* fr, KC__FileReadError error, const char* msg) {
    char* s;
    switch (error) {
        case KC__FILE_READ_ERROR_OPEN:
            s = "Open";
            break;
        case KC__FILE_READ_ERROR_READ:
            s = "Read";
            break;
        case KC__FILE_READ_ERROR_PARSE:
            s = "Parse";
            break;
        default:
            KC__ASSERT(false);
            break;
    }

    if (msg == NULL) {
        LOGGING_ERROR("%s file error [%s]", s, fr->file_name);
    } else {
        LOGGING_ERROR("%s file error (%s) [%s]", s, msg, fr->file_name);
    }

    exit(EXIT_FAILURE);
}

static inline void KC__file_reader_request_buffer(KC__FileReader* fr, KC__Buffer** buffer) {
    KC__Buffer* bf = KC__buffer_queue_get_blank_buffer(fr->buffer_queue);
    KC__ASSERT(bf != NULL);

    KC__BufferType buffer_type;
    switch (fr->input.file_type) {
        case KC__FILE_TYPE_FASTA:
            buffer_type = KC__BUFFER_TYPE_FASTA;
            break;
        case KC__FILE_TYPE_FASTQ:
            buffer_type = KC__BUFFER_TYPE_FASTQ;
            break;
        case KC__FILE_TYPE_SUPER_KMER:
            buffer_type = KC__BUFFER_TYPE_SUPER_KMER;
            break;
        default:
            KC__ASSERT(false);
            break;
    }

    bf->type = buffer_type;

    *buffer = bf;
}

static inline void KC__file_reader_complete_buffer(KC__FileReader* fr, KC__Buffer** buffer) {
    KC__Buffer* bf = *buffer;
    KC__ASSERT(bf != NULL);

    KC__buffer_queue_enqueue_filled_buffer(fr->buffer_queue, bf);
    *buffer = NULL;
}

static inline void KC__file_reader_transfer_data(KC__FileReader* fr, KC__Buffer* current_buffer, KC__Buffer* extra_buffer, size_t extra_size) {
    KC__ASSERT(fr != NULL);
    current_buffer->length -= extra_size;
    const char *src = (char*)(current_buffer->data) + current_buffer->length;
    void* dest = extra_buffer->data;
    memcpy(dest, src, extra_size);
    extra_buffer->length = (uint32_t)extra_size;
}

static inline void KC__file_reader_modify_fasta_buffers(KC__FileReader* fr, KC__Buffer* current_buffer, KC__Buffer* extra_buffer) {
    char* data = current_buffer->data;

    size_t extra_size = 0;
    for (size_t i = current_buffer->length; i > 0; i--) {
        size_t n = i - 1;
        extra_size++;
        if (data[n] == '>')
            break;
    }

    if (extra_size < extra_buffer->size) {
        KC__file_reader_transfer_data(fr, current_buffer, extra_buffer, extra_size);

    } else {
        char* extra_data = extra_buffer->data;
        extra_data[0] = '>';
        extra_data[1] = '\n';

        size_t nt_count = 0;
        for (size_t i = current_buffer->length; i > 0; i--) {
            size_t n = i - 1;

            switch (data[n]) {
                case 'A':
                case 'a':
                case 'C':
                case 'c':
                case 'G':
                case 'g':
                case 'T':
                case 't':
                    extra_data[fr->K - nt_count] = data[n];
                    nt_count++;
                    break;
                default:
                    break;
            }

            if (nt_count == fr->K - 1) {
                break;
            }
        }

        if (nt_count != fr->K - 1) {
            KC__file_reader_process_file_error_exit(fr, KC__FILE_READ_ERROR_PARSE, "Too many unexpected characters");
        }

        extra_buffer->length = (uint32_t)(fr->K + 1);
    }
}

static inline void KC__file_reader_modify_fastq_buffers(KC__FileReader* fr, KC__Buffer* current_buffer, KC__Buffer* extra_buffer) {
    size_t extra_size = 0;

    char* data = current_buffer->data;

    for (size_t i = current_buffer->length; i > 0; i--) {
        size_t n = i - 1;
        extra_size++;
        if (data[n] == '@')
            break;
    }

    if (extra_size >= extra_buffer->size) {
        KC__file_reader_process_file_error_exit(fr, KC__FILE_READ_ERROR_PARSE, "Sequence may be too long");
    }

    KC__file_reader_transfer_data(fr, current_buffer, extra_buffer, extra_size);
}

static void KC__file_reader_read(KC__FileReader* fr, FILE* file, const size_t in_size, void* out, size_t* out_size, bool* end_of_file) {
    if (fr->input.compression_type == KC__FILE_COMPRESSION_TYPE_PLAIN) {
        const size_t read_size = fread(out, 1, in_size, file);
        if (ferror(file)) {
            KC__file_reader_process_file_error_exit(fr, KC__FILE_READ_ERROR_READ, NULL);
        }

        *out_size = read_size;
        *end_of_file = (read_size < in_size);

    } else if (fr->input.compression_type == KC__FILE_COMPRESSION_TYPE_GZIP) {
        size_t remain_size = in_size;
        size_t got_size = 0;
        bool eof = false;

        while (true) {
            if (fr->gz_stream.avail_in == 0) {
                const size_t read_size = fread(fr->gz_data, 1, fr->gz_data_size, file);

                if (ferror(file)) {
                    inflateEnd(&(fr->gz_stream));
                    KC__file_reader_process_file_error_exit(fr, KC__FILE_READ_ERROR_READ, NULL);
                }

                if (read_size == 0) {
                    eof = true;
                    break;
                }

                fr->gz_stream.avail_in = (uint)read_size;
                fr->gz_stream.next_in = fr->gz_data;
            }

            fr->gz_stream.avail_out = (uint)remain_size;
            fr->gz_stream.next_out = (void*)((char*)out + got_size);

            int ret = inflate(&(fr->gz_stream), Z_NO_FLUSH);
            if (ret != Z_OK) {
                bool occur_error = false;
                if (ret == Z_STREAM_END) {
                    // There may be another gz stream concatenated at the end.
                    if (inflateReset(&(fr->gz_stream)) != Z_OK) {
                        occur_error = true;
                    }
                } else {
                    occur_error = true;
                }
                if (occur_error) {
                    inflateEnd(&(fr->gz_stream));
                    KC__file_reader_process_file_error_exit(fr, KC__FILE_READ_ERROR_READ, "GZip file read error");
                }
            }

            got_size += remain_size - fr->gz_stream.avail_out;
            remain_size = fr->gz_stream.avail_out;

            if (remain_size == 0)
                break;
        }

        if (eof) {
            inflateEnd(&(fr->gz_stream));
        }

        *out_size = got_size;
        *end_of_file = eof;

    } else {
        KC__ASSERT(false);
    }
}

static void KC__file_reader_process_reads_file(KC__FileReader* fr) {
    FILE* file = fopen(fr->file_name, "r");
    if (file == NULL) {
        KC__file_reader_process_file_error_exit(fr, KC__FILE_READ_ERROR_OPEN, NULL);
    }

    switch (fr->input.compression_type) {
        case KC__FILE_COMPRESSION_TYPE_PLAIN:
            break;
        case KC__FILE_COMPRESSION_TYPE_GZIP:
            fr->gz_stream.zalloc = Z_NULL;
            fr->gz_stream.zfree = Z_NULL;
            fr->gz_stream.opaque = Z_NULL;
            fr->gz_stream.avail_in = 0;
            fr->gz_stream.next_in = Z_NULL;
            if (inflateInit2(&(fr->gz_stream), 31) != Z_OK) {
                LOGGING_ERROR("Init gz stream failed.");
                exit(EXIT_FAILURE);
            }
            break;
        default:
            KC__ASSERT(false);
            break;
    }

    KC__Buffer* current_buffer;
    KC__Buffer* extra_buffer;
    KC__file_reader_request_buffer(fr, &current_buffer);
    extra_buffer = NULL;

    while (true) {
        char* remain_location = (char*)(current_buffer->data) + current_buffer->length;
        const size_t remain_size = current_buffer->size - current_buffer->length;
        KC__ASSERT(remain_size != 0);

        size_t out_size;
        bool end_of_file;
        KC__file_reader_read(fr, file, remain_size, remain_location, &out_size, &end_of_file);

        current_buffer->length += (uint32_t)out_size;

        if (end_of_file) {
            break;
        }

        KC__ASSERT(extra_buffer == NULL);
        KC__file_reader_request_buffer(fr, &extra_buffer);

        switch (fr->input.file_type) {
            case KC__FILE_TYPE_FASTA:
                KC__file_reader_modify_fasta_buffers(fr, current_buffer, extra_buffer);
                break;
            case KC__FILE_TYPE_FASTQ:
                KC__file_reader_modify_fastq_buffers(fr, current_buffer, extra_buffer);
                break;
            default:
                KC__ASSERT(false);
                break;
        }

        KC__file_reader_complete_buffer(fr, &current_buffer);

        current_buffer = extra_buffer;
        extra_buffer = NULL;
    }

    KC__ASSERT(current_buffer != NULL);
    KC__ASSERT(extra_buffer == NULL);

    KC__file_reader_complete_buffer(fr, &current_buffer);

    fclose(file);
}

static void KC__file_reader_process_super_kmer_file(KC__FileReader* fr) {
    FILE* file = fopen(fr->file_name, "rb");
    if (file == NULL) {
        KC__file_reader_process_file_error_exit(fr, KC__FILE_READ_ERROR_OPEN, NULL);
    }

    while (true) {
        uint32_t buffer_length;

        size_t read_size = fread(&buffer_length, 1, sizeof(uint32_t), file);
        if (ferror(file)) {
            KC__file_reader_process_file_error_exit(fr, KC__FILE_READ_ERROR_READ, NULL);
        }

        if (read_size < sizeof(uint32_t)) {
            if (read_size == 0) {
                break;
            } else {
                KC__file_reader_process_file_error_exit(fr, KC__FILE_READ_ERROR_PARSE, "File is truncated");
            }
        }

        KC__Buffer* buffer;
        KC__file_reader_request_buffer(fr, &buffer);

        KC__ASSERT(buffer_length <= buffer->size);
        buffer->length = buffer_length;

        read_size = fread(buffer->data, 1, buffer_length, file);
        if (ferror(file)) {
            KC__file_reader_process_file_error_exit(fr, KC__FILE_READ_ERROR_READ, NULL);
        }

        if (read_size < buffer_length) {
            KC__file_reader_process_file_error_exit(fr, KC__FILE_READ_ERROR_PARSE, "File is truncated");
        }

        KC__file_reader_complete_buffer(fr, &buffer);
    }

    fclose(file);
}

void* KC__file_reader_work(void* ptr) {
    KC__FileReader* fr = ptr;

    for (size_t i = 0; i < fr->input.files_count; i++) {
        fr->file_name = fr->input.file_names[i];
        LOGGING_DEBUG("Start reading file %s", fr->file_name);

        switch (fr->input.file_type) {
            case KC__FILE_TYPE_FASTA:
            case KC__FILE_TYPE_FASTQ:
                KC__file_reader_process_reads_file(fr);
                break;
            case KC__FILE_TYPE_SUPER_KMER:
                KC__file_reader_process_super_kmer_file(fr);
                break;
            default:
                KC__ASSERT(false);
                break;
        }

        LOGGING_DEBUG("Finish reading file %s", fr->file_name);
    }

    fr->file_name = NULL;

    pthread_exit(NULL);
}