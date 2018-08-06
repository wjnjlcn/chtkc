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


#ifndef KC__KMER_COUNTER_H
#define KC__KMER_COUNTER_H

#include "mem_allocator.h"
#include "param.h"

struct KC__KmerCounter;
typedef struct KC__KmerCounter KC__KmerCounter;

KC__KmerCounter* KC__kmer_counter_create(KC__MemAllocator* mem_allocator, KC__Param* param);
void KC__kmer_counter_free(KC__MemAllocator* mem_allocator, KC__KmerCounter* kmer_counter);

void KC__kmer_counter_work(KC__KmerCounter* kmer_counter);

#endif
