#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//#define DEBUG_PRINT
static int limit_max_coord_move = 0;
static double max_move = -1;


#include "l-bfgs.h"
#include "parallel.h"

// the number of back vectors stored.
static int trigger_restart = 0;
static int m = 10;
static int n = 0;
static int cur_m = 0;
static double *back_s = NULL;
static double *back_y = NULL;
static double (*local_bfgs_f)(double *p) = NULL;
static double (*local_bfgs_fdf)(double *p, double *g) = NULL;
static double *q = NULL;
static double *z = NULL;
static double *alpha = NULL;
static double *beta  = NULL;
static double *rho   = NULL;
static double *cur_x     = NULL;
static double *last_g = NULL;
static double *last_x = NULL;
static int *vec_rotor = NULL;
static double initial_en;
static double *local_x_buffer = NULL;

void do_linesearch( double *search_step );

void l_bfgs_setup( int use_m, // the number of back vectors
		   int n_in,     // the size of the vector space on this node.
		   double *initial_x,
		   double max_move_in,
		   double (*f)(double *p),
		   double (*fdf)(double *p, double *g)	
			 )
{
	if( max_move_in > 0 )
	{
		limit_max_coord_move = 1;
//		printf("limiting the max move.\n");
		max_move = max_move_in;
	}

	m = use_m;
	n = n_in;

	back_s = (double *)malloc( sizeof(double) * n * use_m ); 
	back_y = (double *)malloc( sizeof(double) * n * use_m ); 

	local_bfgs_f = f;
	local_bfgs_fdf = fdf;
	q = (double *)malloc( sizeof(double) * n );
	z = (double *)malloc( sizeof(double) * n );

	alpha = (double *)malloc( sizeof(double) * use_m );
	beta  = (double *)malloc( sizeof(double) * use_m );
	rho   = (double *)malloc( sizeof(double) * use_m );
	last_g = (double *)malloc( sizeof(double) * n );
	last_x = (double *)malloc( sizeof(double) * n );
	cur_x  = (double *)malloc( sizeof(double) * n );
	
	vec_rotor = (int *)malloc( sizeof(int) * use_m );

	for( int i = 0; i < use_m; i++ )
		vec_rotor[i] = i;

	memcpy( cur_x, initial_x, sizeof(double) * n );

	local_x_buffer = (double *)malloc( sizeof(double) * n );
}

void l_bfgs_clear( void )
{
	cur_m = 0;
	free(back_s);
	free(back_y);
	free(q);
	free(z);
	free(alpha);
	free(beta);
	free(rho);
	free(last_g);
	free(last_x);
	free(cur_x);
	free(vec_rotor);
	free(local_x_buffer);
}

int l_bfgs_iteration( double *place )
{
	// evaluate gradient at (g_k)
	initial_en = local_bfgs_fdf( cur_x, q );	

//	printf("initial energy: %.14le\n", initial_en );

	double t;

	int vec_m_0 = vec_rotor[m-1];

	for( int i = m-1; i >= 0; i-- )
	{
		int mm1 = i-1;

		if( mm1 < 0 )
			mm1 = m-1;

		if( mm1 == m-1 )
			vec_rotor[i] = vec_m_0;
		else
			vec_rotor[i] = vec_rotor[mm1];
	}

	if( cur_m > 0 )
	{	
		for( int j = 0; j < n; j++ )
		{
			back_y[vec_rotor[0]*n+j] = q[j] - last_g[j];	 
			back_s[vec_rotor[0]*n+j] = cur_x[j] - last_x[j];
		}

		double dp = 0;

		for( int j = 0; j < n; j++ )
			dp += back_y[vec_rotor[0]*n+j] * back_s[vec_rotor[0]*n+j];


		if( fabs(dp) < 1e-40 )
		{
//			printf("Terminating: dp: %.14le\n", dp );

			return 0;
		}

		rho[vec_rotor[0]] = 1.0 / dp;	
	


//		printf("Computed rho: %le\n", rho[vec_rotor[0]] );

	}	

	memcpy( last_x, cur_x, sizeof(double) * n );
	memcpy( last_g, q, sizeof(double) * n );

	for( int i = 0; i < cur_m; i++ )
	{
		double local_alpha = 0;
		for( int j = 0; j < n; j++ )
			local_alpha += back_s[vec_rotor[i]*n+j] * q[j];
		local_alpha *= rho[vec_rotor[i]];

		for( int j = 0; j < n; j++ )
			q[j] -= local_alpha * back_y[vec_rotor[i]*n+j];

		alpha[vec_rotor[i]] = local_alpha;

//		printf("computed alpha k: %d %le\n", vec_rotor[i], alpha[vec_rotor[i]] );
	} 
	
	double hk0 = 1.0;

	if( cur_m > 0 )
	{
		double den = 0;
		for( int j = 0; j < n; j++ )
			den += back_y[vec_rotor[0]*n+j] * back_y[vec_rotor[0]*n+j];
		
		double num = 0;
		for( int j = 0; j < n; j++ )
			num += back_y[vec_rotor[0]*n+j] * back_s[vec_rotor[0]*n+j];
		hk0 = num / den;
			
	}
	
	for( int j = 0; j < n; j++ )
		z[j] = hk0 * q[j];

	for( int i = cur_m-1; i >= 0; i-- )
	{
		double beta = 0;

		for( int j = 0; j < n; j++ )
			beta += back_y[vec_rotor[i]*n+j] * z[j];		
		beta *= rho[vec_rotor[i]];

		for( int j = 0; j < n; j++ )
			z[j] += back_s[vec_rotor[i]*n+j] * (alpha[vec_rotor[i]] - beta); 
	}

	do_linesearch( z );

	cur_m++;

	if( cur_m > m )
		cur_m = m;

	memcpy( place, cur_x, sizeof(double) * n );

	return 1;
}


double evaluate( double alpha, double *cur_x, double *search_step )
{
//	printf("Line search evaluation.\n");

#ifdef PARALLEL
	MPI_Bcast( &alpha, 1, MPI_DOUBLE, BASE_TASK, MPI_COMM_WORLD );
#endif

	for( int j = 0; j < n; j++ )	
		local_x_buffer[j] = cur_x[j] - search_step[j] * alpha;

	return local_bfgs_f(local_x_buffer);
}

void do_linesearch( double *search_step )
{
	static int line_search_counter = 0;


	double max_alpha = 0.01;

	if( limit_max_coord_move )
	{
		double mm = 0;
		for( int i = 0; i < n; i++ )
		{
			if( fabs(search_step[i]) > mm )
				mm = fabs(search_step[i]);
		}
		
		max_alpha = fabs(max_move/mm);	
	}

#ifdef PARALLEL
	MPI_Bcast( &max_alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD ); 
	MPI_Bcast( &initial_en, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD ); 
#endif

	double alpha_low = 0;
	double f_low = initial_en;

	double alpha_high = max_alpha;

	double f_high = evaluate( alpha_high, cur_x, search_step );
	
	if( f_high < f_low )
	{
		for( int j = 0; j < n; j++ )
			cur_x[j] -= alpha_high * search_step[j];	
	
#ifdef DEBUG_PRINT
		printf("Linesearch: %le\n alpha %lf\n", f_high, alpha_high);
#endif
		return;
	}	
	
	int iter = 0;
	int n_max = 20;
	int n_cat = 300;

	if( f_low < ENERGY_ERROR)
	{
		while( f_high >= ENERGY_ERROR && iter < n_cat )
		{
			alpha_high /= 2;	
			f_high = evaluate( alpha_high, cur_x, search_step );
//			printf("alpha_high: %.14le f_hight: %.14le\n", alpha_high, f_high );
			n_cat++;
		}

		while( f_high < f_low && iter < n_cat)
		{
			alpha_high *= 1.5;	
			f_high = evaluate( alpha_high, cur_x, search_step );
//			printf("2 alpha_high: %.14le f_hight: %.14le\n", alpha_high, f_high);
			n_cat++;
		}
	}

#ifdef PRINT_LINE_SEARCH
//	if( line_search_counter % 3 == 0 )
	{
		for( double alpha = 0; alpha < 1e-6; alpha+= 1e-8 )
		{
			printf("line search %le %.14le\n", alpha, evaluate(alpha, cur_x, search_step ) );
		}
	}
	line_search_counter++;
#endif
//	printf("f_high: %le f_low: %le\n", f_high, f_low );

	if( f_high < f_low )
	{
		for( int j = 0; j < n; j++ )
			cur_x[j] -= alpha_high * search_step[j];	
	
//		printf("Linesearch: %le\n alpha %lf\n", f_high, alpha_high);
		return;
	}	
	
	
	int citer = 0;
	while( f_high < f_low && citer < n_cat )
	{
		alpha_high *=2;
		f_high = evaluate(alpha_high, cur_x, search_step);
		citer++;
	}

//	printf("max alpha: %.14le\n", max_alpha );
	if( alpha_high > max_alpha )
	{
		alpha_high = max_alpha;
		f_high = evaluate(alpha_high, cur_x, search_step);
	} 

	double alpha_mid = (alpha_high+alpha_low)/2;
	double f_mid = evaluate(alpha_mid, cur_x, search_step );
	
	while( f_mid > f_low && citer < n_cat )
	{
		alpha_mid = (alpha_mid+alpha_low)/2;
		f_mid = evaluate(alpha_mid, cur_x, search_step);
		citer++;
	}

	if( f_mid > f_high )
	{	
		printf("TRUNCATING linesearch, max move f_low: %.14le f_mid: %.14le f_high %.14le.\n", f_low, f_mid, f_high);
//		for( int j = 0; j < n; j++ )
//			cur_x[j] -= alpha_high * search_step[j];		
		return;	
	}

#ifdef DEBUG_PRINT
	printf("bracket: %lf %lf %lf, %lf %lf %lf\n", f_low, f_mid, f_high, alpha_low, alpha_mid, alpha_high );	
#endif
//	printf("Linesearch: %le\n alpha %lf\n", f_mid, alpha_mid);

//	return;


	while( iter < n_max)
	{	
		double alpha_try = (alpha_mid+alpha_low)/2;
		double f_try = evaluate(alpha_try,cur_x, search_step);

//		printf("f_try: %lf\n", f_try );

		if( f_try < f_mid )
		{
			alpha_high = alpha_mid;
			f_high = f_mid;

			alpha_mid = alpha_try;
			f_mid = f_try;
		}
		else
		{
			alpha_low = alpha_try;
			f_low = f_try;
		}
		
		alpha_try = (alpha_mid+alpha_high)/2;
		f_try = evaluate(alpha_try,cur_x, search_step);

		if( f_try < f_mid )
		{
			alpha_low = alpha_mid;
			f_low = f_mid;

			alpha_mid = alpha_try;
			f_mid = f_try;
		}
		else
		{
			alpha_high = alpha_try;
			f_high = f_try;
		}

		iter += 1; 
	}	 		

#ifdef DEBUG_PRINT
	printf("final alpha: %le\n", alpha_mid );
#endif	
	for( int j = 0; j < n; j++ )
		cur_x[j] -= alpha_mid * search_step[j];		

#ifdef DEBUG_PRINT
	printf("final_en: %lf\n", local_bfgs_f(cur_x) );
	fflush(stdout);	
#endif
}

int restart_triggered( void )
{
	return trigger_restart; 
}











