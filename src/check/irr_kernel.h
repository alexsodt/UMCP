#ifndef __irrkernelh__
#define __irrkernelh__

struct irr_kernel
{
	int val;
	int ndomains;	
	int ncoords_base;
	double *map; // dimension: 12 * ncoords_base * 3 * ndomains 
	void setup( int valence, int max_domain );
	int domain( double u, double v );
	double *get_map( double *u, double *v );
	void get_map_transform( double *u_u, double *u_v, double *v_u, double *v_v );
};

#endif
