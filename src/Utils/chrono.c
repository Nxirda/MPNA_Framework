#include "chrono.h"

#include <stdio.h>

//
duration_t get_elapsed_time(chrono_t self)
{
	// Error handling
	i64 elapsed_s  = (self.ts_end.tv_sec  - self.ts_start.tv_sec );
	i64 elapsed_ns = (self.ts_end.tv_nsec - self.ts_start.tv_nsec);


	return (duration_t)
	{
		.s  =  elapsed_s,
		.ns =  elapsed_ns
	};
}

//
f64 duration_as_s_f64(duration_t self)
{
	return (f64)self.s + (f64)self.ns * 1.0e-9;
}

//
f64 duration_as_ms_f64(duration_t self)
{
	return (f64)self.s * 1.0e+3 + (f64)self.ns * 1.0e-6;
}

//
f64 duration_as_us_f64(duration_t self)
{
	return (f64)self.s * 1.0e+6 + (f64)self.ns * 1.0e-3;
}

//
f64 duration_as_ns_f64(duration_t self)
{
	return (f64)self.s * 1.0e+9 + (f64)self.ns;
}
