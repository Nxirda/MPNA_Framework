#pragma once 

#include <time.h>
#include "types.h"

//
typedef struct chrono_s {
	struct timespec ts_start;
	struct timespec ts_end;
} chrono_t;


// Could potentially reduce the sizes using unsigned of the most logical size (u16 / u32 ?)
typedef struct duration_s {
	i64 seconds;
	i32 nanoseconds;
} duration_t;

//
#define h hours
#define m minutes
#define s seconds
#define ms miliseconds
#define us microseconds
#define ns nanoseconds

//
static clockid_t __clock__ = CLOCK_MONOTONIC_RAW;

//
static inline i32 start_chrono(chrono_t *chrono)
{
	return clock_gettime(__clock__, &chrono->ts_start);
}

//same as start  (no need to pay for the function call)
static inline i32 stop_chrono(chrono_t *chrono)
{
	return clock_gettime(__clock__, &chrono->ts_end);
}

//
struct timespec get_elapsed_as_timespec(chrono_t self);

//
duration_t get_elapsed_time(chrono_t self);

f64 duration_as_s_f64(duration_t self);
f64 duration_as_ms_f64(duration_t self);
f64 duration_as_us_f64(duration_t self);
f64 duration_as_ns_f64(duration_t self);
