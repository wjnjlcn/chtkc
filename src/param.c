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


#include <argp.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "param.h"
#include "logging.h"


#define KC__OPT_FA 1
#define KC__OPT_FQ 2
#define KC__OPT_GZ 3
#define KC__OPT_COUNT_MAX 4
#define KC__OPT_FILTER_MIN 5
#define KC__OPT_FILTER_MAX 6
#define KC__OPT_BS 7
#define KC__OPT_RT 8
#define KC__OPT_LOG 9


static inline size_t KC__parse_number(struct argp_state* state, const char* arg, const char* info) {
    long n = strtol(arg, NULL, 0);
    if (n <= 0) {
        argp_error(state, "%s value invalid: %ld.", info, n);
    }
    return (size_t)n;
}

static error_t KC__parse_opt(int key, char* arg, struct argp_state* state) {
    KC__Param *param = state->input;
    size_t m;
    size_t n;

    switch (key) {
        case 'k':
            param->K = KC__parse_number(state, arg, "K-mer length");
            break;
        case 't':
            param->threads_count = KC__parse_number(state, arg, "Threads count");
            if (param->threads_count < 3) {
                argp_error(state, "Threads count cannot be less than 3.");
            }
            break;
        case 'm':
            m = strlen(arg);
            n = 0;
            if (m > 0) {
                switch (arg[m - 1]) {
                    case 'M':
                    case 'm':
                        n = 1000000;
                        break;
                    case 'G':
                    case 'g':
                        n = 1000000000;
                        break;
                    default:
                        break;
                }
            }
            if (n == 0) {
                argp_error(state, "Memory size not ends with M/G: %s.", arg);
            }
            param->mem_limit = KC__parse_number(state, arg, "Memory size") * n;
            break;
        case 'o':
            param->output_file_name = arg;
            break;
        case KC__OPT_FA:
            param->input_file_type = KC__FILE_TYPE_FASTA;
            break;
        case KC__OPT_FQ:
            param->input_file_type = KC__FILE_TYPE_FASTQ;
            break;
        case KC__OPT_GZ:
            param->input_compression_type = KC__FILE_COMPRESSION_TYPE_GZIP;
            break;
        case KC__OPT_COUNT_MAX:
            param->output_param.count_max = (KC__count_t)KC__parse_number(state, arg, "Count max");
            break;
        case KC__OPT_FILTER_MIN:
            param->output_param.filter_min = (KC__count_t)KC__parse_number(state, arg, "Filter min");
            break;
        case KC__OPT_FILTER_MAX:
            param->output_param.filter_max = (KC__count_t)KC__parse_number(state, arg, "Filter max");
            break;
        case KC__OPT_BS:
            param->read_buffer_size = (uint32_t)KC__parse_number(state, arg, "Buffer size");
            if (param->read_buffer_size == 0) {
                argp_error(state, "Buffer size cannot be 0.");
            }
            break;
        case KC__OPT_RT:
            param->reading_threads_count = KC__parse_number(state, arg, "Reading threads count");
            if (param->reading_threads_count < 1) {
                argp_error(state, "Reading threads count cannot be less than 1.");
            }
            break;
        case KC__OPT_LOG:
            param->log_file_name = arg;
            break;
        case ARGP_KEY_ARGS:
            param->input_file_names = (state->argv + state->next);
            param->input_files_count = (size_t)(state->argc - state->next);
            break;
        case ARGP_KEY_END:
            if (state->arg_num < 1)
                argp_usage(state);
            if (param->K == 0)
                argp_error(state, "K-mer length value must be provided.");
            if (param->mem_limit == 0)
                argp_error(state, "Memory size value must be provided.");
            if (param->input_file_type == KC__FILE_TYPE_UNKNOWN)
                argp_error(state, "Input file type (fa/fq) should be specified.");
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static inline size_t KC__get_processors_count() {
    long n_processors = sysconf(_SC_NPROCESSORS_ONLN);
    if (n_processors < 0) {
        LOGGING_ERROR("Getting system processors count failed.");
        exit(EXIT_FAILURE);
    } else {
        return (size_t)n_processors;
    }
}

void KC__param_init(KC__Param* param, int argc, char** argv) {
    param->K = 0;
    param->mem_limit = 0;

    param->threads_count = KC__get_processors_count();
    param->reading_threads_count = 0;

    param->input_file_type = KC__FILE_TYPE_UNKNOWN;
    param->input_compression_type = KC__FILE_COMPRESSION_TYPE_PLAIN;

    param->output_file_name = "./KC__output";
    param->log_file_name = NULL;

    param->read_buffer_size = 0;

    param->output_param.filter_min = 2;
    param->output_param.filter_max = KC__COUNT_MAX;
    param->output_param.count_max = 255;

    struct argp_option options[] = {
            {"kmer-len", 'k', "Length", 0, "Length of K-mer", 0},
            {"threads", 't', "N", 0, "Threads count", 0},
            {"mem", 'm', "M/G", 0, "Memory size", 0},
            {"out", 'o', "OUT", 0, "Output file path", 0},

            {"fa", KC__OPT_FA, 0, 0, "Input FASTA format", 1},
            {"fq", KC__OPT_FQ, 0, 0, "Input FASTQ format", 1},
            {"gz", KC__OPT_GZ, 0, 0, "Input GZIP compressed files", 1},

            {"count-max", KC__OPT_COUNT_MAX, "N", 0, "Max count value, default: 255", 2},
            {"filter-min", KC__OPT_FILTER_MIN, "N", 0, "Filter min value, default: 2", 2},
            {"filter-max", KC__OPT_FILTER_MAX, "N", 0, "Filter max value", 2},

            {"log", KC__OPT_LOG, "FILE", 0, "Log file", 3},

            {"bs", KC__OPT_BS, "SIZE", 0, "Buffer size", 4},
            {"rt", KC__OPT_RT, "N", 0, "Reading threads count", 4},
            {0}
    };
    struct argp argp = {options, KC__parse_opt, "FILE...", "Count k-mers."};
    argp_parse(&argp, argc, argv, 0, 0, param);


    param->kmer_processing_threads_count = param->threads_count - 2;

    bool reading_threads_count_provided = (param->reading_threads_count != 0);

    if (!reading_threads_count_provided) {
        param->reading_threads_count = 1;
        if (param->input_compression_type == KC__FILE_COMPRESSION_TYPE_GZIP) {
            param->reading_threads_count = (param->kmer_processing_threads_count + 7) / 8;
        }
    }

    if (param->input_files_count < param->reading_threads_count) {
        param->reading_threads_count = param->input_files_count;
        if (reading_threads_count_provided) {
            LOGGING_WARNING("Reduce reading threads count to number of files: %zu", param->reading_threads_count);
        }
    }


    param->write_buffer_size = 5000000;
    if (param->read_buffer_size == 0) {
        param->read_buffer_size = param->write_buffer_size;
        if (param->input_compression_type == KC__FILE_COMPRESSION_TYPE_GZIP) {
            param->read_buffer_size *= 4;
        }
    }

    param->read_buffers_count = param->kmer_processing_threads_count * 2;
    param->write_buffers_count = param->kmer_processing_threads_count * 2;


    if (param->log_file_name != NULL) {
        KC__LOG_FILE = fopen(param->log_file_name, "w");
        if (KC__LOG_FILE == NULL) {
            KC__LOG_FILE = stderr;
            LOGGING_ERROR("Open log file error: %s", param->log_file_name);
            exit(EXIT_FAILURE);
        }
    }


    LOGGING_DEBUG("K: %zu", param->K);
    LOGGING_DEBUG("Threads count(r/p): %zu(%zu/%zu)", param->threads_count, param->reading_threads_count, param->kmer_processing_threads_count);
    LOGGING_DEBUG("Memory limit: %zu", param->mem_limit);
    for (size_t i = 0; i < param->input_files_count; i++) {
        LOGGING_DEBUG("Input file #%zu: %s", i, param->input_file_names[i]);
    }
    LOGGING_DEBUG("Input file type: %d, compression type: %d", param->input_file_type, param->input_compression_type);
    LOGGING_DEBUG("Output files: %s", param->output_file_name);
    LOGGING_DEBUG("Buffer size(r/w): %zu/%zu, count(r/w): %zu/%zu", param->read_buffer_size, param->write_buffer_size, param->read_buffers_count, param->write_buffers_count);
    LOGGING_DEBUG("Count max: %zu, filter min: %zu, max: %zu", param->output_param.count_max, param->output_param.filter_min, param->output_param.filter_max);
}

void KC__param_destroy(KC__Param* param) {
    if (param->log_file_name != NULL) {
        fclose(KC__LOG_FILE);
        KC__LOG_FILE = stderr;
    }
}