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


#ifndef KC__UTILS_H
#define KC__UTILS_H

#include <stddef.h>
#include <stdbool.h>

size_t KC__calculate_kmer_width_by_unit_size(size_t K, size_t unit_size);
size_t KC__calculate_kmer_width(size_t K);
size_t KC__calculate_kmer_size(size_t K);
size_t KC__max_prime_number(size_t limit);
void KC__calculate_count_field(size_t count_max, size_t* count_bit, size_t* count_size);
void KC__file_error_exit(const char* file_name, const char* action, const char* msg);

#endif
