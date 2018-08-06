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


#ifndef LOGGING_H
#define LOGGING_H


#define LOGGING_LEVEL_DEBUG 1
#define LOGGING_LEVEL_INFO 2
#define LOGGING_LEVEL_WARNING 3
#define LOGGING_LEVEL_ERROR 4
#define LOGGING_LEVEL_CRITICAL 5
#define LOGGING_LEVEL_OFF 6


#ifdef DEBUG
#define LOGGING_LEVEL LOGGING_LEVEL_DEBUG
#else
#define LOGGING_LEVEL LOGGING_LEVEL_INFO
#endif

#ifdef LOGGING_OFF
#undef LOGGING_LEVEL
#define LOGGING_LEVEL LOGGING_LEVEL_OFF
#endif


#define LOGGING_OUTPUT(...) fprintf(KC__LOG_FILE, __VA_ARGS__)


#if LOGGING_LEVEL <= LOGGING_LEVEL_DEBUG
#include <unistd.h>
#include <sys/syscall.h>
#define LOGGING_FMT "%s | %-8s | %ld [%ld] | %s [%s(%d)] "
#define LOGGING_ARGS(LOGGING_TAG) logging_time(), LOGGING_TAG, syscall(SYS_getpid), syscall(SYS_gettid)-syscall(SYS_getpid), __FUNCTION__, logging_src_name(__FILE__), __LINE__
#else
#define LOGGING_FMT "%s | %-8s | "
#define LOGGING_ARGS(LOGGING_TAG) logging_time(), LOGGING_TAG
#endif


#if LOGGING_LEVEL <= LOGGING_LEVEL_DEBUG
#define LOGGING_DEBUG(MSG, ...) LOGGING_OUTPUT(LOGGING_FMT MSG "\n", LOGGING_ARGS("DEBUG"), ## __VA_ARGS__)
#else
#define LOGGING_DEBUG(MSG, ...)
#endif

#if LOGGING_LEVEL <= LOGGING_LEVEL_INFO
#define LOGGING_INFO(MSG, ...) LOGGING_OUTPUT(LOGGING_FMT MSG "\n", LOGGING_ARGS("INFO"), ## __VA_ARGS__)
#else
#define LOGGING_INFO(MSG, ...)
#endif

#if LOGGING_LEVEL <= LOGGING_LEVEL_WARNING
#define LOGGING_WARNING(MSG, ...) LOGGING_OUTPUT(LOGGING_FMT MSG "\n", LOGGING_ARGS("WARNING"), ## __VA_ARGS__)
#else
#define LOGGING_WARNING(MSG, ...)
#endif

#if LOGGING_LEVEL <= LOGGING_LEVEL_ERROR
#define LOGGING_ERROR(MSG, ...) LOGGING_OUTPUT(LOGGING_FMT MSG "\n", LOGGING_ARGS("ERROR"), ## __VA_ARGS__)
#else
#define LOGGING_ERROR(MSG, ...)
#endif

#if LOGGING_LEVEL <= LOGGING_LEVEL_CRITICAL
#define LOGGING_CRITICAL(MSG, ...) LOGGING_OUTPUT(LOGGING_FMT MSG "\n", LOGGING_ARGS("CRITICAL"), ## __VA_ARGS__)
#else
#define LOGGING_CRITICAL(MSG, ...)
#endif


#include <stdio.h>

FILE* KC__LOG_FILE;

char* logging_time();
char* logging_src_name(char* src_path);

#endif
