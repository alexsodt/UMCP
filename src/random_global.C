#define __randomglobalc__

#include "gsl/gsl_randist.h"
#include "random_global.h"

void my_gsl_reseed( int seed )
{
	if( rng_x )
		gsl_rng_set( rng_x, seed );
}

void init_random( int seed )
{
        rng_x = gsl_rng_alloc(rng_T);
	gsl_rng_env_setup();

	my_gsl_reseed( seed);
}
