#ifndef __lbfgsh__
#define __lbfgsh__
#define ENERGY_ERROR (1e37)

void l_bfgs_setup( int use_m, // the number of back vectors
		   int n_in,     // the size of the vector space on this node.
		   double *initial_x,
			double max_move,
		   double (*f)(double *p),
		   double (*fdf)(double *p, double *g)	
			 );

int l_bfgs_iteration( double *place );
void l_bfgs_clear( void );


#endif
