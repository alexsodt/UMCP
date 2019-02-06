#ifndef __randomglobalh__
#define __randomglobalh__

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

void init_random( int seed );

#ifdef __randomglobalc__

const gsl_rng_type *rng_T = gsl_rng_default;
gsl_rng *rng_x = NULL;

#else

extern const gsl_rng_type *rng_T;
extern gsl_rng *rng_x;

#endif


#endif
