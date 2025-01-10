#include "chrono.h"

#include <stdio.h>

//
duration_t get_elapsed_time(chrono_t self)
{
	// Error handling
	i64 elapsed_s  = (self.ts_end.tv_sec  - self.ts_start.tv_sec );
	i64 elapsed_ns = (self.ts_end.tv_nsec - self.ts_start.tv_nsec);

	/* i64 hours = elapsed_s / 3600;
	elapsed_s %= 3600;

	i64 mins = elapsed_s / 60;
	elapsed_s %= 60; 

	i64 milis = elapsed_ns / 1000000;
	elapsed_ns %= 1000000;

	i64 micros = elapsed_ns / 1000;
	elapsed_ns %= 1000; */
	
	return (duration_t)
	{
		//.h  =  hours,
		//.m  =  mins,
		.s  =  elapsed_s,
		//.ms =  milis,
		//.us =  micros,
		.ns =  elapsed_ns
	};
}

// Need to add correct return things like : FIRST_IS_HIGHER instead of 1 or smth
i8 compare_durations(duration_t a, duration_t b)
{
	if(a.s > b.s)
		return 1;
	else if(a.s < b.s)
		return -1;

	if(a.ns > b.ns)
		return 1;
	else if(a.ns < b.ns)
		return -1;
	
	return 0;
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

// Might clean that up a bit and make it more appealing
void print_chrono(chrono_t self)
{
	// Error handling
	//duration_t elapsed = get_elapsed_time(self);
	
	printf("Measured time :\n");
	/* printf("%d h; %d m; %d s; %d ms; %d us; %d ns\n", elapsed.h , elapsed.m , elapsed.s , 
							  elapsed.ms, elapsed.us, elapsed.ns);	*/
}
