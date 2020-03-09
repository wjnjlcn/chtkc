cmake_minimum_required(VERSION 3.0.2)
project(chtkc C)

set(CMAKE_C_STANDARD 99)


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
    add_definitions(-DDEBUG)
endif()

message(STATUS "Build: ${CMAKE_BUILD_TYPE}")


set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wpointer-arith")

set(SRC
        src/assert.h
        src/logging.h src/logging.c
        src/mem_allocator.h src/mem_allocator.c
        src/queue.h src/queue.c
        src/buffer_queue.h src/buffer_queue.c
        src/file_reader.h src/file_reader.c
        src/file_writer.h src/file_writer.c
        src/hash_map.h src/hash_map.c
        src/kmer_processor.h src/kmer_processor.c
        src/param.h src/param.c
        src/utils.h src/utils.c
        src/kmer_counter.h src/kmer_counter.c
        src/header.h src/header.c
        src/histo.h src/histo.c
        src/dump.h src/dump.c
        src/pthread_barrier.h src/pthread_barrier.c
        src/main.c)

add_executable(chtkc ${SRC})
target_link_libraries(chtkc pthread z)

add_executable(chtkco ${SRC})
target_compile_definitions(chtkco PRIVATE -DKC__MEM_OPT)
target_link_libraries(chtkco pthread z)

install(TARGETS chtkc chtkco DESTINATION bin)


if (KC__ENABLE_TESTS)
    message(STATUS "Enable tests")

    set(LIB_SRC ${SRC})
    list(REMOVE_ITEM LIB_SRC src/main.c)

    include(FindPkgConfig)
    pkg_search_module(CHECK check)

    set(TESTS_SRC
            ${LIB_SRC}
            tests/check_all.h
            tests/check_queue.c
            tests/check_buffer_queue.c
            tests/check_file_reader.c
            tests/check_file_writer.c
            tests/check_hash_map.c
            tests/check_kmer_processor.c
            tests/check_main.c)

    add_executable(check_chtkc ${TESTS_SRC})
    target_link_libraries(check_chtkc ${CHECK_STATIC_LDFLAGS} z)

    add_executable(check_chtkco ${TESTS_SRC})
    target_compile_definitions(check_chtkco PRIVATE -DKC__MEM_OPT)
    target_link_libraries(check_chtkco ${CHECK_STATIC_LDFLAGS} z)
endif()