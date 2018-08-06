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


#include <time.h>
#include <string.h>
#include "logging.h"


static char logging_time_str[255];

char* logging_time() {
    time_t t;
    time(&t);
    struct tm* local_t = localtime(&t);
    strftime(logging_time_str, 255, "%Y-%m-%d %H:%M:%S", local_t);
    return logging_time_str;
}

char* logging_src_name(char* src_path) {
    char* src_name = strrchr(src_path, '/');
    if (src_name)
        src_name += 1;
    else
        src_name = src_path;
    return src_name;
}