#define __gslrandomglobals__
#include "gsl/gsl_randist.h"
gsl_rng * r_gen_global = NULL;

void check_random_init( void )
{	
	if( !r_gen_global )
	{
		r_gen_global = gsl_rng_alloc (gsl_rng_taus);
		gsl_rng_set( r_gen_global, rand() );
	}
}
