#ifndef KC__CHECK_ALL_H
#define KC__CHECK_ALL_H

#include <check.h>

Suite* queue_suite();
Suite* buffer_queue_suite();
Suite* kmer_processor_suite();
Suite* hash_map_suite();
Suite* file_reader_suite();
Suite* file_writer_suite();

#endif
