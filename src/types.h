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


#ifndef KC__TYPES_H
#define KC__TYPES_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>


#define KC__UNIT_BIT 64
typedef uint64_t KC__unit_t;
#define KC__UNIT_MAX UINT64_MAX

#define KC__COUNT_BIT 32
typedef uint32_t KC__count_t;
#define KC__COUNT_MAX UINT32_MAX

#ifdef KC__MEM_OPT
typedef uint32_t KC__node_id_t;
#define KC__NODE_ID_MAX UINT32_MAX
#else
typedef uint64_t KC__node_id_t;
#define KC__NODE_ID_MAX UINT64_MAX
#endif

#define KC__NODE_ID_NULL 0

typedef enum {
    KC__FILE_TYPE_FASTA = 0,
    KC__FILE_TYPE_FASTQ,
    KC__FILE_TYPE_SUPER_KMER,
    KC__FILE_TYPE_UNKNOWN
} KC__FileType;

typedef enum {
    KC__FILE_COMPRESSION_TYPE_PLAIN = 0,
    KC__FILE_COMPRESSION_TYPE_GZIP
} KC__FileCompressionType;

#endif //KC__TYPES_H
