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


#include <stdlib.h>
#include "mem_allocator.h"
#include "logging.h"
#include "assert.h"


struct KC__MemAllocator {
    size_t limit;
    size_t available;
    size_t allocated_count;
    size_t freed_count;
};


KC__MemAllocator* KC__mem_allocator_create(size_t mem_limit) {
    KC__MemAllocator* ma = (KC__MemAllocator*)malloc(sizeof(KC__MemAllocator));
    if (ma == NULL) {
        LOGGING_CRITICAL("Allocating memory for Mem Allocator failed.");
        exit(EXIT_FAILURE);
    }

    ma->limit = mem_limit;
    ma->available = ma->limit - sizeof(KC__MemAllocator);
    ma->allocated_count = 0;
    ma->freed_count = 0;

    return ma;
}

KC__MemAllocator* KC__mem_allocator_free(KC__MemAllocator* ma) {
    LOGGING_DEBUG("Mem           limit: %zu", ma->limit);
    LOGGING_DEBUG("Mem            used: %zu", ma->limit - ma->available);
    LOGGING_DEBUG("Mem allocated count: %zu", ma->allocated_count);
    LOGGING_DEBUG("Mem     freed count: %zu", ma->freed_count);
    KC__ASSERT(ma->allocated_count == ma->freed_count);
    free(ma);
}

size_t KC__mem_available(const KC__MemAllocator* ma) {
    return ma->available;
}

void* KC__mem_alloc(KC__MemAllocator* ma, size_t size, const char* name) {
    if (size <= ma->available) {
        void* mem = malloc(size);
        if (mem != NULL) {
            ma->allocated_count++;
            ma->available -= size;
            return mem;
        }
    }
    LOGGING_CRITICAL("Allocating memory for %s failed.", name);
    exit(EXIT_FAILURE);
}

void* KC__mem_aligned_alloc(KC__MemAllocator* ma, size_t size, const char* name) {
    if (size <= ma->available) {
        void* mem;
        int result = posix_memalign(&mem, 64, size);
        if (result == 0) {
            KC__ASSERT((size_t)mem % 64 == 0);
            ma->allocated_count++;
            ma->available -= size;
            return mem;
        }
    }
    LOGGING_CRITICAL("Allocating memory for %s failed.", name);
    exit(EXIT_FAILURE);
}

void KC__mem_free(KC__MemAllocator* ma, void* ptr) {
    free(ptr);
    ma->freed_count++;
}