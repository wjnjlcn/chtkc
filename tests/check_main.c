#include <stdlib.h>
#include "check_all.h"
#include "../src/logging.h"


int main() {
    KC__LOG_FILE = stderr;

    SRunner* sr = srunner_create(NULL);

    srunner_add_suite(sr, queue_suite());
    srunner_add_suite(sr, buffer_queue_suite());
    srunner_add_suite(sr, hash_map_suite());
    srunner_add_suite(sr, kmer_processor_suite());
    srunner_add_suite(sr, file_reader_suite());
    srunner_add_suite(sr, file_writer_suite());


    srunner_run_all(sr, CK_NORMAL);
    int number_failed = srunner_ntests_failed(sr);

    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}