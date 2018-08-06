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


#include "header.h"


bool KC__write_header(const KC__Header* header, FILE* file) {
    uint64_t arr[4] = {header->K, header->count_max, header->filter_min, header->filter_max};
    for (size_t i = 0; i < 4; i++) {
        const size_t write_size = fwrite(&(arr[i]), 1, sizeof(uint64_t), file);
        if (write_size < sizeof(uint64_t)) {
            return false;
        }
    }
    return true;
}

bool KC__read_header(KC__Header* header, FILE* file) {
    uint64_t arr[4];
    for (size_t i = 0; i < 4; i++) {
        fread(&(arr[i]), 1, sizeof(uint64_t), file);
        if (ferror(file)) {
            return false;
        }
    }
    header->K = arr[0];
    header->count_max = arr[1];
    header->filter_min = arr[2];
    header->filter_max = arr[3];
    return true;
}