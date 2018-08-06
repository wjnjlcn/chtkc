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


#ifndef KC__MEM_ALLOCATOR_H
#define KC__MEM_ALLOCATOR_H

#include <stddef.h>


struct KC__MemAllocator;
typedef struct KC__MemAllocator KC__MemAllocator;

KC__MemAllocator* KC__mem_allocator_create(size_t mem_limit);
KC__MemAllocator* KC__mem_allocator_free(KC__MemAllocator* mem_allocator);

size_t KC__mem_available(const KC__MemAllocator* mem_allocator);

void* KC__mem_alloc(KC__MemAllocator* mem_allocator, size_t size, const char* name);
void* KC__mem_aligned_alloc(KC__MemAllocator* mem_allocator, size_t size, const char* name);
void KC__mem_free(KC__MemAllocator* mem_allocator, void* ptr);

#endif
