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
#include <limits.h>

#include "utils.h"
#include "types.h"
#include "logging.h"


size_t KC__calculate_kmer_width_by_unit_size(size_t K, size_t unit_size) {
    size_t K_unit = unit_size * CHAR_BIT / 2;
    return (K + K_unit - 1) / K_unit;
}

size_t KC__calculate_kmer_width(size_t K) {
    return KC__calculate_kmer_width_by_unit_size(K, sizeof(KC__unit_t));
}

size_t KC__calculate_kmer_size(size_t K) {
    size_t kmer_width = KC__calculate_kmer_width(K);
    return sizeof(KC__unit_t) * kmer_width;
}


static inline bool KC__is_prime_number(size_t n) {
    if (n < 2)
        return false;
    for (size_t i = 2; i * i <= n; i++) {
        if (n % i == 0)
            return false;
    }
    return true;
}

size_t KC__max_prime_number(size_t limit) {
    size_t n = limit;
    while (!KC__is_prime_number(n) && (n > 0))
        n--;
    if (n == 0) {
        LOGGING_CRITICAL("Error getting prime number limited by %zu.", limit);
        exit(EXIT_FAILURE);
    }
    return n;
}

void KC__calculate_count_field(size_t count_max, size_t* count_bit, size_t* count_size) {
    if (count_max <= UINT8_MAX) {
        *count_bit = 8;
        *count_size = sizeof(uint8_t);
    } else if (count_max <= UINT16_MAX) {
        *count_bit = 16;
        *count_size = sizeof(uint16_t);
    } else if (count_max <= UINT32_MAX) {
        *count_bit = 32;
        *count_size = sizeof(uint32_t);
    } else {
        *count_bit = KC__COUNT_BIT;
        *count_size = sizeof(KC__count_t);
    }
}

void KC__file_error_exit(const char* file_name, const char* action, const char* msg) {
    if (msg == NULL) {
        LOGGING_ERROR("%s file error [%s]", action, file_name);
    } else {
        LOGGING_ERROR("%s file error (%s) [%s]", action, msg, file_name);
    }

    exit(EXIT_FAILURE);
}
