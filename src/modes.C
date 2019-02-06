#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "interp.h"
#include "mutil.h"
#include "lapack_we_use.h"
#include "gsl/gsl_sf.h"
#include "parallel.h"
#include "sparse.h"

static int maxv = 13;
static double rel_tol = 1e-20;

void surface::mode_perturb( double *r, double *ro, int qi, int qj, int nx, int ny, double Lx, double Ly, double mag, int do_cosine)
{
	double qx = qi * 2 * M_PI / Lx;
	double qy = qj * 2 * M_PI / Ly;

	if( qi > nx/2 ) qx = -(nx-qi) * 2 * M_PI / Lx;
	if( qj > ny/2 ) qy = -(ny-qj) * 2 * M_PI / Ly;
	double q2 = qx*qx+qy*qy;
//	double mode_mag = sqrt( mag / ( q2*q2 * Lx * Ly) ); 

	if( do_cosine )
	{
		for( int v = 0; v < nv; v++ ) 
			r[3*v+2] += mag * cos( qx * ro[3*v+0] + qy * ro[3*v+1] ); 
	}
	else
	{
		for( int v = 0; v < nv; v++ ) 
			r[3*v+2] += mag * sin( qx * ro[3*v+0] + qy * ro[3*v+1] ); 
	}
}


static double *d_hc_d_eps = NULL;
static double *d_hs_d_eps = NULL;
static double *d_r_d_eps = NULL;
static int do_spherical = 0;
static int do_set = 0;
static int n_modes = 0;

void surface::project_grad_to_modes( double *r, double *g, double *mg )
{
	if( do_set || do_spherical )
	{
		// change to dgemv to improve performance
		memset( mg, 0, sizeof(double) * n_modes );

		for( int m = 0; m < n_modes; m++ )
		for( int v = 0; v < 3*nv; v++ )
		{
			if( do_spherical ) 
				mg[m] += g[v] * d_r_d_eps[m*(3*nv)+v];
			else
				mg[m] += g[v] * d_hc_d_eps[m*(3*nv)+v];
		}	
	}
	else
	{
		mg[0] = 0;
		mg[1] = 0;

		for( int i = 0; i < 3*nv; i++ )
			mg[0] += d_hs_d_eps[i] * g[i];
		for( int i = 0; i < 3*nv; i++ )
			mg[1] += d_hc_d_eps[i] * g[i];
	}
}

int surface::setup_spherical_set( double *ro, int l_min, int l_max, double **output_qvals )
{
	do_set = 1;
	do_spherical = 1;

	double *Amat = (double *)malloc( sizeof(double) * (nv) * (nv) );
	memset( Amat, 0, sizeof(double) * (nv) * (nv ) );

	for( int v = 0; v < nv; v++ )
	{
		int val = theVertices[v].valence;

		Amat[(v)*nv+(v)] += 0.5;

		double w = 1.0 / (val * 2);

		for( int e = 0; e < val; e++ )
		{
			int j = theVertices[v].edges[e];

			Amat[(v)*nv+(j)] += w;
		}
	}

	int N = nv;
	int LDA = nv;
	int info = 0;

	double *oamat = (double *)malloc( sizeof(double) *  nv * nv );
	memcpy( oamat, Amat, sizeof(double) * nv * nv );
	
	int lwork = N*N; 
	double *work = (double *)malloc( sizeof(double) * lwork );
	int ipiv[N];
	
	
	int nset = 0;

	for( int l = l_min; l <= l_max; l++ )
	{
#ifdef M_MIN
			for( int m = M_MIN; m <= M_MAX; m++ )
#else
			for( int m = -l; m <= l; m++ )
#endif
			{
				nset++;
			}
	}
 
	*output_qvals = (double *)malloc( sizeof(double) * nset );
		
	nset = 0;

	for( int l = l_min; l <= l_max; l++ )
	{
#ifdef M_MIN
			for( int m = M_MIN; m <= M_MAX; m++ )
#else
			for( int m = -l; m <= l; m++ )
#endif
			{
				(*output_qvals)[nset] = l; // scaled later
				nset++;
			}
	} 
	n_modes = nset;

	double *r0_pos = (double *)malloc( sizeof(double) * 3 * nv );
	
	for( int v = 0; v < nv; v++ )
	{
		int val = theVertices[v].valence;

		r0_pos[3*v+0] = 0.5 * theVertices[v].r[0];
		r0_pos[3*v+1] = 0.5 * theVertices[v].r[1];
		r0_pos[3*v+2] = 0.5 * theVertices[v].r[2];

		double w = 1.0 / (val * 2);

		for( int e = 0; e < val; e++ )
		{
			int j = theVertices[v].edges[e];

			r0_pos[3*v+0] += w * theVertices[j].r[0];
			r0_pos[3*v+1] += w * theVertices[j].r[1];
			r0_pos[3*v+2] += w * theVertices[j].r[2];
		}
	}
	
	d_r_d_eps = (double *)malloc( sizeof(double) * 3 * nv * nset );
	memset( d_r_d_eps, 0, sizeof(double) * 3 * nv * nset );

#ifdef INVERT_DGETRF
	dgetrf( &N, &N, Amat, &LDA, ipiv, &info );
	if( info != 0 )
	{
		printf("info: %d Failure to invert A (dgetrf).\n", info);
		exit(1);
	} 
	dgetri( &N, Amat, &LDA, ipiv, work, &lwork, &info );
	
	if( info != 0 )
	{
		printf("info: %d Failure to invert A (dgetri).\n", info);
		exit(1);
	} 
#else
	
	double *mass_bvec = (double *)malloc( sizeof(double) * nset * 3 * nv ); // nset, and x,y,z.



	int m_cntr = 0;
	for( int ql = l_min; ql <= l_max; ql++ )
#ifdef M_MIN
	for( int qm = M_MIN; qm <= M_MAX; qm++, m_cntr++ )
#else
	for( int qm = -ql; qm <= ql; qm++, m_cntr++)
#endif
	for( int j = 0; j < nv; j++ )
	{		
		int l = ql;
		int m = qm;

		int am = m;
		if( m < 0 )
			am = -m;
	
		double *ro = r0_pos+3*j;
		double r = sqrt(ro[0]*ro[0]+ro[1]*ro[1]+ro[2]*ro[2]);
		double rn[3] = { ro[0]/r, ro[1]/r, ro[2]/r };
		double phi = acos( rn[2] );
		double theta = atan2( rn[1], rn[0] );
		double ylm = gsl_sf_legendre_sphPlm( l, am, rn[2] );

		if( m < 0 )
			ylm *= sqrt(2.0) * sin( am * theta );
		else if( m > 0 )
			ylm *= sqrt(2.0) * cos( m * theta );

		double dr[3] = { ylm * ro[0], ylm * ro[1], ylm * ro[2] };

		for( int v = 0; v < nv; v++ )
		{		
			// the parameter
			for( int vc = 0; vc < 3; vc++ )	
				mass_bvec[m_cntr*3*nv+vc*nv+j] = dr[vc];
		}
	}
	m_cntr=0;
	char uplo = 'U';
	int nrhs = 3*nset;
	dposv( &uplo, &N, &nrhs, oamat, &N, mass_bvec, &N, &info );
	printf("DPOSV: %d\n", info );
	//dgesv( &N, &nrhs, oamat, &N, ipiv, mass_bvec, &N, &info );
	//printf("DGESV: %d\n", info );
	

	for( int vec = 0; vec < nset; vec++ )
	{
		for( int vc = 0; vc < 3; vc++ )
		for( int v = 0; v < nv; v++ )
		{
			d_r_d_eps[vec*(3*nv)+3*v+vc] = mass_bvec[vec*3*nv+vc*nv+v];
		}
	}
	free(mass_bvec);
#endif
	


#ifdef INVERT_DGETRF

	int m_cntr = 0;
	for( int ql = l_min; ql <= l_max; ql++ )
#ifdef M_MIN
	for( int qm = M_MIN; qm <= M_MAX; qm++, m_cntr++ )
#else
	for( int qm = -ql; qm <= ql; qm++, m_cntr++)
#endif
	for( int j = 0; j < nv; j++ )
	{		
		int l = ql;
		int m = qm;

		int am = m;
		if( m < 0 )
			am = -m;
	
		double *ro = r0_pos+3*j;
		double r = sqrt(ro[0]*ro[0]+ro[1]*ro[1]+ro[2]*ro[2]);
		double rn[3] = { ro[0]/r, ro[1]/r, ro[2]/r };
		double phi = acos( rn[2] );
		double theta = atan2( rn[1], rn[0] );
		double ylm = gsl_sf_legendre_sphPlm( l, am, rn[2] );

		if( m < 0 )
			ylm *= sqrt(2.0) * sin( am * theta );
		else if( m > 0 )
			ylm *= sqrt(2.0) * cos( m * theta );

		double dr[3] = { ylm * ro[0], ylm * ro[1], ylm * ro[2] };

		for( int v = 0; v < nv; v++ )
		{		
			// the parameter
			for( int vc = 0; vc < 3; vc++ )
			{
					d_r_d_eps[m_cntr * 3 * nv + 3*v+vc] += Amat[(j)+(nv)*(v)] * dr[vc];
			}
		}
	}
#endif

	free(Amat);

	return n_modes;
}


void surface::setup_spherical_perturb( double *ro, int ql, int qm )
{
	do_spherical = 1;

	double *Amat = (double *)malloc( sizeof(double) * (3*nv) * (3*nv) );
	memset( Amat, 0, sizeof(double) * (3*nv) * (3*nv ) );

	for( int v = 0; v < nv; v++ )
	{
		int val = theVertices[v].valence;

		Amat[(3*v+0)*3*nv+(3*v+0)] += 0.5;
		Amat[(3*v+1)*3*nv+(3*v+1)] += 0.5;
		Amat[(3*v+2)*3*nv+(3*v+2)] += 0.5;

		double w = 1.0 / (val * 2);

		for( int e = 0; e < val; e++ )
		{
			int j = theVertices[v].edges[e];

			Amat[(3*v+0)*3*nv+(3*j+0)] += w;
			Amat[(3*v+1)*3*nv+(3*j+1)] += w;
			Amat[(3*v+2)*3*nv+(3*j+2)] += w;
		}
	}

#if 0	
	printf("A = {\n");
	for( int i = 0; i < 3*nv; i++ )
	{
		printf("{");
		for( int j = 0; j < 3 * nv; j++ )
		{
			printf("%lf", Amat[i*3*nv+j] );
			if( j != 3*nv-1 )
				printf(",");
		}
		printf("}");
		if( i != 3*nv-1 )
			printf(",\n");
		else
			printf("\n");
	}
	printf("};\n");
#endif

	int N = 3*nv;
	int LDA = 3*nv;
	int info = 0;

	double *oamat = (double *)malloc( sizeof(double) * 3 * nv * 3 * nv );
	memcpy( oamat, Amat, sizeof(double) * 3 * nv * 3 * nv );
	
	int lwork = N*N; 
	double *work = (double *)malloc( sizeof(double) * lwork );
	int ipiv[N];

	dgetrf( &N, &N, Amat, &LDA, ipiv, &info );
	if( info != 0 )
	{
		printf("info: %d Failure to invert A (dgetrf).\n", info);
		exit(1);
	} 
	dgetri( &N, Amat, &LDA, ipiv, work, &lwork, &info );
	
	if( info != 0 )
	{
		printf("info: %d Failure to invert A (dgetri).\n", info);
		exit(1);
	} 

#if 0 
	printf("Ainv = {\n");
	for( int i = 0; i < 3*nv; i++ )
	{
		printf("{");
		for( int j = 0; j < 3 * nv; j++ )
		{
			printf("%lf", Amat[i*3*nv+j] );
			if( j != 3*nv-1 )
				printf(",");
		}
		printf("}");
		if( i != 3*nv-1 )
			printf(",\n");
		else
			printf("\n");
	}
	printf("};\n");
#endif
	
	double *r0_pos = (double *)malloc( sizeof(double) * 3 * nv );
	
	for( int v = 0; v < nv; v++ )
	{
		int val = theVertices[v].valence;

		r0_pos[3*v+0] = 0.5 * theVertices[v].r[0];
		r0_pos[3*v+1] = 0.5 * theVertices[v].r[1];
		r0_pos[3*v+2] = 0.5 * theVertices[v].r[2];

		double w = 1.0 / (val * 2);

		for( int e = 0; e < val; e++ )
		{
			int j = theVertices[v].edges[e];

			r0_pos[3*v+0] += w * theVertices[j].r[0];
			r0_pos[3*v+1] += w * theVertices[j].r[1];
			r0_pos[3*v+2] += w * theVertices[j].r[2];
		}
	}

	d_r_d_eps = (double *)malloc( sizeof(double) * 3 * nv );
	memset( d_r_d_eps, 0, sizeof(double) * 3 * nv );
	for( int j = 0; j < nv; j++ )
	{		
		int l = ql;
		int m = qm;
		int am = m;
		if( m < 0 )
			am = -m;

	
		double *ro = r0_pos+3*j;
		double r = sqrt(ro[0]*ro[0]+ro[1]*ro[1]+ro[2]*ro[2]);
		double rn[3] = { ro[0]/r, ro[1]/r, ro[2]/r };
		double phi = acos( rn[2] );
		double theta = atan2( rn[1], rn[0] );
		double ylm = gsl_sf_legendre_sphPlm( l, am, rn[2] );

		if( m < 0 )
			ylm *= sqrt(2.0) * sin( am * theta );
		else if( m > 0 )
			ylm *= sqrt(2.0) * cos( m * theta );

		double dr[3] = { ylm * ro[0], ylm * ro[1], ylm * ro[2] };

		for( int v = 0; v < nv; v++ )
		{		
			// the parameter
			for( int vc = 0; vc < 3; vc++ )
			for( int jc = 0; jc < 3; jc++ ) 
			{
				d_r_d_eps[3*v+vc] += Amat[(3*j+jc)+(3*nv)*(3*v+vc)] * dr[jc];
			}
		}
	}	

	n_modes = 1;

#if 0
	for( int j = 0; j < nv; j++ )
	{
		int l = ql;
		int m = qm;
		int am = m;
		if( m < 0 )
			am = -m;

		double *ro = r0_pos+3*j;
		double r = sqrt(ro[0]*ro[0]+ro[1]*ro[1]+ro[2]*ro[2]);
		printf("ro: %lf %lf %lf\n", ro[0], ro[1], ro[2] );
		double rn[3] = { ro[0]/r, ro[1]/r, ro[2]/r };
		double phi = acos( rn[2] );
		double theta = atan2( rn[1], rn[0] );
		double ylm = gsl_sf_legendre_sphPlm( l, am, rn[2] );

		if( m < 0 )
			ylm *= sin( am * theta );
		else if( m > 0 )
			ylm *= cos( m * theta );

		double dr[3] = { ylm * ro[0], ylm * ro[1], ylm * ro[2] };
		// the spot

		double po[3] = {0,0,0};
		for( int k = 0; k < nv; k++ )
		{
			po[0] += oamat[(3*j+0)*(3*nv)+(3*k+0)] * d_r_d_eps[3*k+0];
			po[1] += oamat[(3*j+1)*(3*nv)+(3*k+1)] * d_r_d_eps[3*k+1];
			po[2] += oamat[(3*j+2)*(3*nv)+(3*k+2)] * d_r_d_eps[3*k+2];
		}
	
		printf("dr: %lf %lf %lf po: %lf %lf %lf\n", dr[0], dr[1], dr[2], po[0], po[1], po[2] );
	}

	exit(1);
#endif

	free(Amat);
}

int surface::setup_planar_set( double *ro, int l_min, int l_max, double **output_qvals )
{
	do_set = 1;
	double Lx = PBC_vec[0][0];
	double Ly = PBC_vec[1][1];

	double *Amat = (double *)malloc( sizeof(double) * (nv) * (nv) );
	memset( Amat, 0, sizeof(double) * (nv) * (nv ) );

	for( int v = 0; v < nv; v++ )
	{
		int val = theVertices[v].valence;

		Amat[(v)*nv+(v)] += 0.5;

		double w = 1.0 / (val * 2);

		for( int e = 0; e < val; e++ )
		{
			int j = theVertices[v].edges[e];

			Amat[(v)*nv+(j)] += w;
		}
	}

	int N = nv;
	int LDA =nv;
	int info = 0;

	double *oamat = (double *)malloc( sizeof(double) *  nv *  nv );
	memcpy( oamat, Amat, sizeof(double) *  nv *  nv );
	
	int lwork = N*N; 
	double *work = (double *)malloc( sizeof(double) * lwork );
	int ipiv[N];

#ifdef INVERT_DGETRF
	printf("DGETRF/DGETRI N: %d\n", N);
	dgetrf( &N, &N, Amat, &LDA, ipiv, &info );
	if( info != 0 )
	{
		printf("info: %d Failure to invert A (dgetrf).\n", info);
		exit(1);
	} 
	dgetri( &N, Amat, &LDA, ipiv, work, &lwork, &info );
	
	if( info != 0 )
	{
		printf("info: %d Failure to invert A (dgetri).\n", info);
		exit(1);
	} 
#endif
	double *r0_pos = (double *)malloc( sizeof(double) * 3 * nv );
	
	for( int v = 0; v < nv; v++ )
	{
		int val = theVertices[v].valence;

		r0_pos[3*v+0] = 0.5 * theVertices[v].r[0];
		r0_pos[3*v+1] = 0.5 * theVertices[v].r[1];
		r0_pos[3*v+2] = 0.5 * theVertices[v].r[2];

		double w = 1.0 / (val * 2);

		for( int e = 0; e < val; e++ )
		{
			int j = theVertices[v].edges[e];

			r0_pos[3*v+0] += w * theVertices[j].r[0];
			r0_pos[3*v+1] += w * theVertices[j].r[1];
			r0_pos[3*v+2] += w * theVertices[j].r[2];
		}
	}
	
	int nset = 0;
	
	for( int lx = 0; lx < (l_max+1); lx++ )
	for( int ly = 0; ly < (l_max+1); ly++ )
	for( int sc = 0; sc < 2; sc++ )
	{
		if( lx < l_min && ly < l_min )
			continue;

		int lx_mag = lx;
		int ly_mag = ly;
		
		if( lx_mag == 0 && ly_mag == 0 )
			continue;
		
		nset++;
	}
	n_modes = nset;
	
	*output_qvals = (double *)malloc( sizeof(double) * nset );

	d_hc_d_eps = (double *)malloc( sizeof(double) * 3 * nv * nset );
	d_hs_d_eps = NULL;

	memset( d_hc_d_eps, 0, sizeof(double) * 3 * nv * nset );

	int vec = 0;


#ifndef INVERT_DGETRF
	double *mass_bvec = (double *)malloc( sizeof(double) * nset * 3 * nv ); // nset, and x,y,z.

	for( int lx = 0; lx < (l_max+1); lx++ )
	for( int ly = 0; ly < (l_max+1); ly++ )
	for( int sc = 0; sc < 2; sc++ )
	{
		if( lx < l_min && ly < l_min )
			continue;

		int lx_mag = lx; 
		int ly_mag = ly; 
		
		if( lx_mag == 0 && ly_mag == 0 )
			continue;

		double qx = 2 * M_PI * lx_mag / Lx;
		double qy = 2 * M_PI * ly_mag / Ly;

		(*output_qvals)[vec] = sqrt(qx*qx+qy*qy);

		for( int j = 0; j < nv; j++ )
		{		
			// the spot
			double rc[3] = { 0,0, cos(qx*ro[3*j+0]+qy*ro[3*j+1]) };
			if( sc == 1 )
				rc[2] = sin(qx*ro[3*j+0]+qy*ro[3*j+1]);

			// h_q Exp( I q )
			// h_q exp( I (n-q)

			for( int vc = 0; vc < 3; vc++ )	
				mass_bvec[vec*3*nv+vc*nv+j] = rc[vc];
		}
		
		vec++;
	}

	char uplo = 'U';
	int nrhs = 3*nset;
	dposv( &uplo, &N, &nrhs, oamat, &N, mass_bvec, &N, &info );
	printf("DPOSV: %d\n", info );
	//dgesv( &N, &nrhs, oamat, &N, ipiv, mass_bvec, &N, &info );
	//printf("DGESV: %d\n", info );
	

	for( int vec = 0; vec < nset; vec++ )
	{
		for( int vc = 0; vc < 3; vc++ )
		for( int v = 0; v < nv; v++ )
		{
			d_hc_d_eps[vec*(3*nv)+3*v+vc] = mass_bvec[vec*3*nv+vc*nv+v];
		}
	}
	free(mass_bvec);
#endif

#ifdef INVERT_DGETRF
	vec =0;
	for( int lx = 0; lx < (l_max+1); lx++ )
	for( int ly = 0; ly < (l_max+1); ly++ )
	for( int sc = 0; sc < 2; sc++ )
	{
		if( lx < l_min && ly < l_min )
			continue;
		int lx_mag = lx; 
		int ly_mag = ly; 
		
		if( lx_mag == 0 && ly_mag == 0 )
			continue;

		double qx = 2 * M_PI * lx_mag / Lx;
		double qy = 2 * M_PI * ly_mag / Ly;

		(*output_qvals)[vec] = sqrt(qx*qx+qy*qy);

		for( int j = 0; j < nv; j++ )
		{		
			// the spot
			double rc[3] = { 0,0, cos(qx*ro[3*j+0]+qy*ro[3*j+1]) };
			if( sc == 1 )
				rc[2] = sin(qx*ro[3*j+0]+qy*ro[3*j+1]);

			// h_q Exp( I q )
			// h_q exp( I (n-q)
	
			for( int v = 0; v < nv; v++ )
			{		
				// the parameter
				for( int vc = 0; vc < 3; vc++ ) // should be diagonal but what the heck
					d_hc_d_eps[vec * (3*nv) + 3*v+vc] += Amat[(j)*(nv)+(v)] * rc[vc];
			}

		}
		
		vec++;
	}
#endif
	
	free(Amat);
	free(oamat);
	return nset;
}

void surface::setup_mode_perturb( double *ro, int qi, int qj, int nx, int ny, double Lx, double Ly )
{
	double qx = qi * 2 * M_PI / Lx;
	double qy = qj * 2 * M_PI / Ly;

	if( qi > nx/2 ) qx = -(nx-qi) * 2 * M_PI / Lx;
	if( qj > ny/2 ) qy = -(ny-qj) * 2 * M_PI / Ly;


	for( int v = 0; v < nv; v++ )
	{
		if( theVertices[v].valence != 6 )
		{
			// assumes it's regular so that Amat is symmetric.
			printf("Setup mode perturb called with irregular mesh.\n");
			exit(1);
		}
	}

	double *Amat = (double *)malloc( sizeof(double) * (3*nv) * (3*nv) );
	memset( Amat, 0, sizeof(double) * (3*nv) * (3*nv ) );

	for( int v = 0; v < nv; v++ )
	{
		int val = theVertices[v].valence;

		Amat[(3*v+0)*3*nv+(3*v+0)] += 0.5;
		Amat[(3*v+1)*3*nv+(3*v+1)] += 0.5;
		Amat[(3*v+2)*3*nv+(3*v+2)] += 0.5;

		double w = 1.0 / (val * 2);

		for( int e = 0; e < val; e++ )
		{
			int j = theVertices[v].edges[e];

			Amat[(3*v+0)*3*nv+(3*j+0)] += w;
			Amat[(3*v+1)*3*nv+(3*j+1)] += w;
			Amat[(3*v+2)*3*nv+(3*j+2)] += w;
		}
	}

#if 0	
	printf("A = {\n");
	for( int i = 0; i < 3*nv; i++ )
	{
		printf("{");
		for( int j = 0; j < 3 * nv; j++ )
		{
			printf("%lf", Amat[i*3*nv+j] );
			if( j != 3*nv-1 )
				printf(",");
		}
		printf("}");
		if( i != 3*nv-1 )
			printf(",\n");
		else
			printf("\n");
	}
	printf("};\n");
#endif
	char uplo = 'U';

	int N = 3*nv;
	int LDA = 3*nv;
	int info = 0;

	double *oamat = (double *)malloc( sizeof(double) * 3 * nv * 3 * nv );
	memcpy( oamat, Amat, sizeof(double) * 3 * nv * 3 * nv );

	dpotrf( &uplo, &N, Amat, &LDA, &info );
	if( info != 0 )
	{
		printf("info: %d Failure to invert A, is it not positive definite?\n", info);
		exit(1);
	} 
	dpotri( &uplo, &N, Amat, &LDA, &info );
	
	if( info != 0 )
	{
		printf("info: %d Failure to invert A, is it not positive definite?\n", info );
		exit(1);
	}
	
	for( int x = 0; x < 3*nv; x++ )
	{
		for( int y = x; y < 3*nv; y++ )
			Amat[x*(3*nv)+y] = Amat[y*(3*nv)+x]; 
	}
#if 0
	printf("Ainv = {\n");
	for( int i = 0; i < 3*nv; i++ )
	{
		printf("{");
		for( int j = 0; j < 3 * nv; j++ )
		{
			printf("%lf", Amat[i*3*nv+j] );
			if( j != 3*nv-1 )
				printf(",");
		}
		printf("}");
		if( i != 3*nv-1 )
			printf(",\n");
		else
			printf("\n");
	}
	printf("};\n");
#endif
	d_hc_d_eps = (double *)malloc( sizeof(double) * 3 * nv );
	d_hs_d_eps = (double *)malloc( sizeof(double) * 3 * nv );
	memset( d_hc_d_eps, 0, sizeof(double) * 3 * nv );
	memset( d_hs_d_eps, 0, sizeof(double) * 3 * nv );
	for( int j = 0; j < nv; j++ )
	{		
		// the spot
		double rc[3] = { 0,0,cos(qx*ro[3*j+0]+qy*ro[3*j+1]) };
		double rs[3] = { 0,0,sin(qx*ro[3*j+0]+qy*ro[3*j+1]) };

		for( int v = 0; v < nv; v++ )
		{		
			// the parameter
			for( int vc = 0; vc < 3; vc++ )
			for( int jc = 0; jc < 3; jc++ ) // should be diagonal but what the heck
			{
				d_hc_d_eps[3*v+vc] += Amat[(3*j+jc)*(3*nv)+(3*v+vc)] * rc[jc];
				d_hs_d_eps[3*v+vc] += Amat[(3*j+jc)*(3*nv)+(3*v+vc)] * rs[jc];
			}
		}
	}	

	n_modes = 2;

#if 0
	for( int j = 0; j < nv; j++ )
	{
		// the spot
		double rc[3] = { 0,0,cos(qx*ro[3*j+0]+qy*ro[3*j+1]) };
		double rs[3] = { 0,0,sin(qx*ro[3*j+0]+qy*ro[3*j+1]) };

		double po[3] = {0,0,0};
		for( int k = 0; k < nv; k++ )
		{
			po[0] += oamat[(3*j+0)*(3*nv)+3*k+0] * d_hc_d_eps[3*k+0];
			po[1] += oamat[(3*j+1)*(3*nv)+3*k+1] * d_hc_d_eps[3*k+1];
			po[2] += oamat[(3*j+2)*(3*nv)+3*k+2] * d_hc_d_eps[3*k+2];
		}
		printf("rc: %lf %lf %lf po: %lf %lf %lf\n", rc[0], rc[1], rc[2], po[0], po[1], po[2] );
	}

	exit(1);
#endif

	free(Amat);
}

void surface::mode_perturb( double *r,  double mag, int do_cosine ) // cosine is now which mode for spherical harmonics.
{
	if( do_spherical )
	{
		for( int v = 0; v < nv; v++ ) 	
		{
			r[3*v+0] += mag * d_r_d_eps[do_cosine*3*nv+3*v+0]; 
			r[3*v+1] += mag * d_r_d_eps[do_cosine*3*nv+3*v+1]; 
			r[3*v+2] += mag * d_r_d_eps[do_cosine*3*nv+3*v+2];
		} 
	}
	else if( do_set )
	{
		for( int v = 0; v < nv; v++ ) 	
		{
			r[3*v+0] += mag * d_hc_d_eps[do_cosine*3*nv+3*v+0]; 
			r[3*v+1] += mag * d_hc_d_eps[do_cosine*3*nv+3*v+1]; 
			r[3*v+2] += mag * d_hc_d_eps[do_cosine*3*nv+3*v+2];
		} 
	}
	else
	{
		if( do_cosine )
		{
			for( int v = 0; v < nv; v++ ) 	
			{
				r[3*v+0] += mag * d_hc_d_eps[3*v+0]; 
				r[3*v+1] += mag * d_hc_d_eps[3*v+1]; 
				r[3*v+2] += mag * d_hc_d_eps[3*v+2];
			} 
		}
		else
		{
			for( int v = 0; v < nv; v++ ) 
			{
				r[3*v+0] += mag * d_hs_d_eps[3*v+0]; 
				r[3*v+1] += mag * d_hs_d_eps[3*v+1]; 
				r[3*v+2] += mag * d_hs_d_eps[3*v+2];
			} 
		}
	}
}


void surface::load_least_squares_fitting( force_set *theForceSet )
{
	double *Amat = (double *)malloc( sizeof(double) * (nv) * (nv) );
	memset( Amat, 0, sizeof(double) * (nv) * (nv ) );


	theForceSet->frc_coef = (double *)malloc( sizeof(double) * maxv * theForceSet->npts );
	theForceSet->frc_coef_list = (int *)malloc( maxv * theForceSet->npts * sizeof(int) );
	theForceSet->frc_ncoef = (int *)malloc( theForceSet->npts  * sizeof(int) );
	
	
	for( int p = 0; p < theForceSet->npts; p++ )
	{
		/*
			for each point, the cross-correlation of each pair.
		*/

		int f = theForceSet->face_set[p];
		double u = theForceSet->uv_set[2*p+0];
		double v = theForceSet->uv_set[2*p+1];

		double coefs[MAX_VALENCE+6];
		int coef_list[MAX_VALENCE+6];
		int ncoef=0;
		get_pt_coeffs( f, u, v, coefs, coef_list, &ncoef ); 

		theForceSet->frc_ncoef[p] = ncoef;

		for( int c = 0; c < ncoef; c++ )
		{
			theForceSet->frc_coef[p*maxv+c] = coefs[c];
			theForceSet->frc_coef_list[p*maxv+c] = coef_list[c];
		}

		for( int c1 = 0; c1 < ncoef; c1++ )
		for( int c2 = 0; c2 < ncoef; c2++ )
			Amat[coef_list[c1]*nv+coef_list[c2]] += coefs[c1] * coefs[c2];			
	}


#if 0	
	printf("A = {\n");
	for( int i = 0; i < nv; i++ )
	{
		printf("{");
		for( int j = 0; j <  nv; j++ )
		{
			printf("%lf", Amat[i*nv+j] );
			if( j != nv-1 )
				printf(",");
		}
		printf("}");
		if( i != nv-1 )
			printf(",\n");
		else
			printf("\n");
	}
	printf("};\n");
#endif

	int N = nv;
	int LDA = nv;
	int info = 0;

	double *oamat = (double *)malloc( sizeof(double) * nv * nv );
	memcpy( oamat, Amat, sizeof(double) * nv * nv );
	
	int lwork = N*N; 
	double *work = (double *)malloc( sizeof(double) * lwork );
	int ipiv[N];

	dgetrf( &N, &N, Amat, &LDA, ipiv, &info );
	if( info != 0 )
	{
		printf("info: %d Failure to invert A (dgetrf).\n", info);
		exit(1);
	} 
	dgetri( &N, Amat, &LDA, ipiv, work, &lwork, &info );
	
	if( info != 0 )
	{
		printf("info: %d Failure to invert A (dgetri).\n", info);
		exit(1);
	} 

#if 0 
	printf("Ainv = {\n");
	for( int i = 0; i < nv; i++ )
	{
		printf("{");
		for( int j = 0; j < nv; j++ )
		{
			printf("%lf", Amat[i*nv+j] );
			if( j != nv-1 )
				printf(",");
		}
		printf("}");
		if( i != nv-1 )
			printf(",\n");
		else
			printf("\n");
	}
	printf("};\n");
#endif

	free(work);
	free(oamat);	

	theForceSet->frc_inverse = Amat;	
//	free(Amat);
}

void surface::test_force_set( force_set *theForceSet )
{
	double *r = (double *)malloc( sizeof(double ) * (3*nv+3) );
	r[3*nv+0] = 1.0;
	r[3*nv+1] = 1.0;
	r[3*nv+2] = 1.0;

	get(r);

	FILE *baseFile = fopen( "base.xyz", "w");
	fprintf(baseFile, "%d\n", theForceSet->npts );
	fprintf(baseFile, "base file\n" );
	
	FILE *targetFile = fopen( "target.xyz", "w");
	fprintf(targetFile, "%d\n", theForceSet->npts );
	fprintf(targetFile, "target file\n" );
	
	FILE *solutionFile = fopen( "solution.xyz", "w");
	fprintf(solutionFile, "%d\n", theForceSet->npts );
	fprintf(solutionFile, "solution file\n" );
	
	double *displ = (double *)malloc( sizeof(double) * 3 * nv );
	double *b = (double *)malloc( sizeof(double) *  3 * nv );
	memset( displ, 0, sizeof(double ) * 3 * nv );
	memset( b, 0, sizeof(double ) * 3 * nv );
	for( int x = 0; x < theForceSet->npts; x++ )
	{
		double rp[3];
		double np[3];
	
		int f = theForceSet->face_set[x];
		double u = theForceSet->uv_set[x*2+0];
		double v = theForceSet->uv_set[x*2+1];

		evaluateRNRM( f, u, v, rp, np, r );
		
		double tx = rp[0];
		double ty = rp[1];
		double tz = rp[2];

		double r = sqrt(tx*tx+ty*ty+tz*tz);

		double mag = 50;

//		double displ[3] = { mag*tx * tx * tz/r/r/r, -mag*tz*ty * tx/r/r/r, mag*tz * tz * tz/r/r/r };
		double displ[3] = { mag*(2*(double)rand()/(double)RAND_MAX-1), 
				    mag*(2*(double)rand()/(double)RAND_MAX-1),
				    mag*(2*(double)rand()/(double)RAND_MAX-1) };

		fprintf(baseFile, "C %lf %lf %lf\n", rp[0], rp[1], rp[2] );	
		fprintf(targetFile, "O %lf %lf %lf\n", rp[0]+displ[0], rp[1]+displ[1], rp[2]+displ[2] );	
		
		for( int c = 0; c < theForceSet->frc_ncoef[x]; c++ )
		{
			b[theForceSet->frc_coef_list[x*maxv+c]*3+0] += displ[0] * theForceSet->frc_coef[x*maxv+c];
			b[theForceSet->frc_coef_list[x*maxv+c]*3+1] += displ[1] * theForceSet->frc_coef[x*maxv+c];
			b[theForceSet->frc_coef_list[x*maxv+c]*3+2] += displ[2] * theForceSet->frc_coef[x*maxv+c];
		}		
	}	
		
	char trans = 'T';
	int incrxy = 3;
	double zero = 0.0;
	double one = 1.0;

	for( int c = 0; c < 3; c++ )
		dgemv( &trans, &nv, &nv, &one, theForceSet->frc_inverse, &nv, b+c, &incrxy, &zero, displ+c, &incrxy );

	for( int x = 0; x < nv; x++ )
	{
		r[3*x+0] += displ[3*x+0];
		r[3*x+1] += displ[3*x+1];
		r[3*x+2] += displ[3*x+2];
	}

	for( int x = 0; x < theForceSet->npts; x++ )
	{
		double rp[3];
		double np[3];
	
		int f = theForceSet->face_set[x];
		double u = theForceSet->uv_set[x*2+0];
		double v = theForceSet->uv_set[x*2+1];

		evaluateRNRM( f, u, v, rp, np, r );

		fprintf(solutionFile, "N %lf %lf %lf\n", rp[0], rp[1], rp[2] );	
	}	

	free(b);
}

void surface::evaluate_fp( force_set *theForceSet, double *xp, double *r )
{
	for( int x = 0; x < theForceSet->npts; x++ )
	{
		double np[3];
	
		int f = theForceSet->face_set[x];
		double u = theForceSet->uv_set[x*2+0];
		double v = theForceSet->uv_set[x*2+1];

		evaluateRNRM( f, u, v, xp+3*x, np, r );
	}	
}

void surface::debug_dKE_dx_and2( force_set *theForceSet, double *effM, double *vp, double *unit_f_vec, int f, double u, double v )
{
#if 0
	double dKEdx = 0;
	double d2KEdx2 = 0;

	dKE_dx_and2( theForceSet, effM, vp, unit_f_vec, f, u, v, &dKEdx, &d2KEdx2 );

	double T_cur = evaluate_T( theForceSet, vp, mesh_p );

	int neps = 100;
	double maxeps = 1;

	double *vcopy = (double *)malloc( sizeof(double) * 3 * nv );
	
	double *g = (double *)malloc( sizeof(double) * 3 * nv );

	for( int ieps = 0; ieps < neps; ieps++ )
	{
		memcpy( vcopy, vp, sizeof(double) * 3 * nv );

		double eps = ieps * (maxeps)/(double)neps;
	
		double frc[3] = { unit_f_vec[0] * eps,
				  unit_f_vec[1] * eps,
				  unit_f_vec[2] * eps};

		memset( g, 0, sizeof(double) * 3 * nv );			
		applyForceAtPoint( f, u, v, frc, g, theForceSet );

		// perturb the velocity

		for( int v1 = 0; v1 < nv; v1++ )
		for( int v2 = 0; v2 < nv; v2++ )
		{
			vcopy[3*v2+0] += (-g[3*v1+0]*effM[v2*nv+v1]);
			vcopy[3*v2+1] += (-g[3*v1+1]*effM[v2*nv+v1]);
			vcopy[3*v2+2] += (-g[3*v1+2]*effM[v2*nv+v1]);
		}
		
		double T_new = evaluate_T( theForceSet, vcopy );

		printf("%le T: %le approx %le\n", eps, T_new-T_cur, dKEdx * eps + d2KEdx2 * 0.5 * eps * eps  );
	}	

	free(vcopy);
	free(g);
#endif
}

void surface::dKE_dx_and2( force_set *theForceSet, double *effM, double *vp, double *unit_f_vec, int f, double u, double v, double *dKEdx, double *d2KEdx2 )
{
	// apply force at this point, what is the change in kinetic energy?

	// dv[v1] = \lambda effMv[v1,v2] * f[v2]

	double f_prop[3]={0,0,0};
	
	double coefs[MAX_VALENCE+6];
	int coef_list[MAX_VALENCE+6];
	int ncoef=0;
	get_pt_coeffs( f, u, v, coefs, coef_list, &ncoef ); 

	double *dv_coef = (double *)malloc( sizeof(double) * nv );
	memset( dv_coef, 0, sizeof(double) * nv );

	for( int v2 = 0; v2 < nv; v2++ )
	{
		for( int c1 = 0; c1 < ncoef; c1++ )
			dv_coef[v2] += effM[v2*nv+coef_list[c1]] * coefs[c1];

	}
	
	*dKEdx = 0;
	*d2KEdx2 = 0;

	double *mass = theForceSet->mass;

	for( int x = 0; x < theForceSet->npts; x++ )
	{
		// can skip out here if dv_coef is negligible.
		double vk[3] = { 0,0,0};

		// the virtual particle's velocity:

		for( int c = 0; c < theForceSet->frc_ncoef[x]; c++ )
		{
			vk[0] += theForceSet->frc_coef[x*maxv+c] * vp[theForceSet->frc_coef_list[x*maxv+c]*3+0];
			vk[1] += theForceSet->frc_coef[x*maxv+c] * vp[theForceSet->frc_coef_list[x*maxv+c]*3+1];
			vk[2] += theForceSet->frc_coef[x*maxv+c] * vp[theForceSet->frc_coef_list[x*maxv+c]*3+2];
		}	

		// the change in the virtual particle's velocity with the force magnitude:

		double sum = 0;
		for( int c1 = 0; c1 < theForceSet->frc_ncoef[x]; c1++ )
		{
			double coef = dv_coef[theForceSet->frc_coef_list[x*maxv+c1]] * theForceSet->frc_coef[x*maxv+c1];

			*dKEdx -= mass[x] * vk[0] * unit_f_vec[0] * coef; 
			*dKEdx -= mass[x] * vk[1] * unit_f_vec[1] * coef; 
			*dKEdx -= mass[x] * vk[2] * unit_f_vec[2] * coef; 
			
			sum += dv_coef[theForceSet->frc_coef_list[x*maxv+c1]] * theForceSet->frc_coef[x*maxv+c1];
		}

		*d2KEdx2 += mass[x] * unit_f_vec[0]*unit_f_vec[0]* sum*sum;
		*d2KEdx2 += mass[x] * unit_f_vec[1]*unit_f_vec[1]* sum*sum;
		*d2KEdx2 += mass[x] * unit_f_vec[2]*unit_f_vec[2]* sum*sum;
	}	

	free(dv_coef);
}

// the momentum

void surface::evaluate_momentum( force_set *theForceSet, double *vq, double *pmesh, double *pout )
{
	for( int x = 0; x < theForceSet->npts; x++ )
	{
		double vk[3] = { 0,0,0};

		for( int c = 0; c < theForceSet->frc_ncoef[x]; c++ )
		{
			vk[0] += theForceSet->frc_coef[x*maxv+c] * vq[theForceSet->frc_coef_list[x*maxv+c]*3+0];
			vk[1] += theForceSet->frc_coef[x*maxv+c] * vq[theForceSet->frc_coef_list[x*maxv+c]*3+1];
			vk[2] += theForceSet->frc_coef[x*maxv+c] * vq[theForceSet->frc_coef_list[x*maxv+c]*3+2];
		}	

		pout[0] += vk[0] * theForceSet->mass[x];
		pout[1] += vk[1] * theForceSet->mass[x];
		pout[2] += vk[2] * theForceSet->mass[x];

	}	
}

double surface::evaluate_T( double *vq, double *pmesh, double *vq_m1, double *pmesh_m1, double *alphas )
{
	double T = 0;
	

/*
	for( int x = 0; x < theForceSet->npts; x++ )
	{
		double vk[3] = { 0,0,0};

		for( int c = 0; c < theForceSet->frc_ncoef[x]; c++ )
		{
			vk[0] += theForceSet->frc_coef[x*maxv+c] * vq[theForceSet->frc_coef_list[x*maxv+c]*3+0];
			vk[1] += theForceSet->frc_coef[x*maxv+c] * vq[theForceSet->frc_coef_list[x*maxv+c]*3+1];
			vk[2] += theForceSet->frc_coef[x*maxv+c] * vq[theForceSet->frc_coef_list[x*maxv+c]*3+2];
		}	

		T += 0.5 * theForceSet->mass[x] * (vk[0]*vk[0]+vk[1]*vk[1]+vk[2]*vk[2]);		
	}	
*/

	double alt_T = 0;

	if( vq_m1 && pmesh_m1 )
	{
		for( int vx = 0; vx < par_info.nv; vx++ )
		{
			int v = par_info.verts[vx];

			alt_T += 0.5 * (pmesh[3*v+0]+pmesh_m1[3*v+0])/2 * (vq[3*v+0]+vq_m1[3*v+0])/2;
			alt_T += 0.5 * (pmesh[3*v+1]+pmesh_m1[3*v+1])/2 * (vq[3*v+1]+vq_m1[3*v+1])/2;
			alt_T += 0.5 * (pmesh[3*v+2]+pmesh_m1[3*v+2])/2 * (vq[3*v+2]+vq_m1[3*v+2])/2;
		}
	}
	else
	{
		for( int vx = 0; vx < par_info.nv; vx++ )
		{
			int v = par_info.verts[vx];
			alt_T += 0.5 * pmesh[3*v+0] * vq[3*v+0];
			alt_T += 0.5 * pmesh[3*v+1] * vq[3*v+1];
			alt_T += 0.5 * pmesh[3*v+2] * vq[3*v+2];
		}
	}
#ifdef PARALLEL
	ParallelSum(&alt_T,1);
#endif
//	printf("T: %.14le alt_T: %.14le\n", T, alt_T );

	return alt_T;
}

void surface::getSparseEffectiveMass( force_set * theForceSet, double *effective_mass, int *use_map, int *nuse, SparseMatrix **theMatrix )
{
	memset( effective_mass, 0, sizeof(double) * nv * nv );

	for( int x = 0; x < theForceSet->npts; x++ )
	{
		for( int c1 = 0; c1 < theForceSet->frc_ncoef[x]; c1++ )
		for( int c2 = 0; c2 < theForceSet->frc_ncoef[x]; c2++ )
			effective_mass[theForceSet->frc_coef_list[x*maxv+c1]*nv+theForceSet->frc_coef_list[x*maxv+c2]] += theForceSet->frc_coef[x*maxv+c1] * theForceSet->frc_coef[x*maxv+c2] * theForceSet->mass[x];
	}	

	// invert it.

	int N = nv;
	int LDA = nv;
	int info = 0;
	char uplo = 'U';
	
	dpotrf( &uplo, &N, effective_mass, &LDA, &info );

	if( info != 0 )
	{
		printf("info: %d Failure to invert effective mass matrix, is it not positive definite?\n", info);
		exit(1);
	} 

	dpotri( &uplo, &N, effective_mass, &LDA, &info );
	
	if( info != 0 )
	{
		printf("info: %d Failure to invert effective mass matrix, is it not positive definite?\n", info );
		exit(1);
	}

	for( int c1 = 0; c1 < nv; c1++ )
	for( int c2 = c1+1; c2 < nv; c2++ )
	{
		effective_mass[c1*nv+c2] = effective_mass[c2*nv+c1];	
	}

	
#ifdef PARALLEL
	ParallelBroadcast(effective_mass,nv*nv);
#endif


	(*theMatrix) = new SparseMatrix;

	(*theMatrix)->init( nv ); 

	for( int v2 = 0; v2 < nv; v2++ )
	{
		int use_it = 0;
		for( int v1x = 0; v1x < par_info.nv; v1x++ )
		{
			int v1 = par_info.verts[v1x];



			if( fabs(effective_mass[v1*nv+v2]) / sqrt(fabs(effective_mass[v1*nv+v1]*effective_mass[v2*nv+v2])) > rel_tol )
			{
				(*theMatrix)->coupleParameters( v1, v2, effective_mass[v1*nv+v2] );
				use_it=1;
			}
		}
				
		if( use_it )
		{
			use_map[*nuse] = v2;
			(*nuse)++;
		}
	}
	
	(*theMatrix)->setNeedSource();
	(*theMatrix)->compress();


	printf("nuse: %d nv: %d\n", *nuse, par_info.nv );

	double *temp_mat = (double *)malloc( sizeof(double) * (*nuse) * par_info.nv  );

	for( int v1x = 0; v1x < par_info.nv; v1x++ )
	{
		int v1 = par_info.verts[v1x];

		for( int v2x = 0; v2x < (*nuse); v2x++ )
		{
			int v2 = use_map[v2x];

			temp_mat[v1x*(*nuse)+v2x] = effective_mass[v1*nv+v2];
		}
	}	

	memcpy( effective_mass, temp_mat, sizeof(double) * (*nuse) * par_info.nv );

	free(temp_mat);

//#define DISTANCE_CHECK
#ifdef DISTANCE_CHECK

	double *rms_node_hop = ( double *) malloc( sizeof(double) * nv );
	double *nrms_node_hop = ( double *) malloc( sizeof(double) * nv );
	memset( rms_node_hop, 0, sizeof(double) * nv );
	memset( nrms_node_hop, 0, sizeof(double) * nv );
	
	int *hop_mat = (int *)malloc( sizeof(int) * nv * nv );


	for( int v = 0; v < nv; v++ )
	{
		int min_hops[nv];
	
		for( int v2 = 0; v2 < nv; v2++ )
			min_hops[v2] = nv;
		min_hops[v] = 0;

		int done = 0;

		while( !done )
		{
			done = 1;

			for( int v1 = 0; v1 < nv; v1++ )
			{
				int val = theVertices[v1].valence;

				for( int e = 0; e < val; e++ )
				{
					if( min_hops[theVertices[v1].edges[e]]+1 <
					    min_hops[v1] )
					{
						min_hops[v1] = min_hops[theVertices[v1].edges[e]]+1;
						done = 0;
					}
				}
			}
		}
	
		memcpy( hop_mat + v * nv, min_hops, sizeof(int) * nv ); 

	}

	int max_hop = 0;
	for( int x = 0; x < nv; x++ )
	for( int y = 0; y < nv; y++ )
	{
		if( hop_mat[x*nv+y] > max_hop ) max_hop = hop_mat[x*nv+y];
		rms_node_hop[hop_mat[x*nv+y]] += effective_mass[x*nv+y] * effective_mass[x*nv+y];
		nrms_node_hop[hop_mat[x*nv+y]] += 1;	
	}

	for( int e = 0; e < max_hop; e++ )
		printf("%d %le\n", e, sqrt(rms_node_hop[e]/nrms_node_hop[e]) );
	free(hop_mat);
	free(rms_node_hop);
	free(nrms_node_hop);
#endif

	
}

void surface::getEffectiveMass( force_set * theForceSet, double *effective_mass )
{
	memset( effective_mass, 0, sizeof(double) * nv * nv );

	for( int x = 0; x < theForceSet->npts; x++ )
	{
		for( int c1 = 0; c1 < theForceSet->frc_ncoef[x]; c1++ )
		for( int c2 = 0; c2 < theForceSet->frc_ncoef[x]; c2++ )
			effective_mass[theForceSet->frc_coef_list[x*maxv+c1]*nv+theForceSet->frc_coef_list[x*maxv+c2]] += theForceSet->frc_coef[x*maxv+c1] * theForceSet->frc_coef[x*maxv+c2] * theForceSet->mass[x];
	}	

	// invert it.

	int N = nv;
	int LDA = nv;
	int info = 0;
	char uplo = 'U';
	
	dpotrf( &uplo, &N, effective_mass, &LDA, &info );

	if( info != 0 )
	{
		printf("info: %d Failure to invert effective mass matrix, is it not positive definite?\n", info);
		exit(1);
	} 

	dpotri( &uplo, &N, effective_mass, &LDA, &info );
	
	if( info != 0 )
	{
		printf("info: %d Failure to invert effective mass matrix, is it not positive definite?\n", info );
		exit(1);
	}

	for( int c1 = 0; c1 < nv; c1++ )
	for( int c2 = c1+1; c2 < nv; c2++ )
	{
		effective_mass[c1*nv+c2] = effective_mass[c2*nv+c1];	
	}

#define DISTANCE_CHECK
#ifdef DISTANCE_CHECK

	double *rms_node_hop = ( double *) malloc( sizeof(double) * nv );
	double *nrms_node_hop = ( double *) malloc( sizeof(double) * nv );
	memset( rms_node_hop, 0, sizeof(double) * nv );
	memset( nrms_node_hop, 0, sizeof(double) * nv );
	
	int *hop_mat = (int *)malloc( sizeof(int) * nv * nv );


	for( int v = 0; v < nv; v++ )
	{
		int min_hops[nv];
	
		for( int v2 = 0; v2 < nv; v2++ )
			min_hops[v2] = nv;
		min_hops[v] = 0;

		int done = 0;

		while( !done )
		{
			done = 1;

			for( int v1 = 0; v1 < nv; v1++ )
			{
				int val = theVertices[v1].valence;

				for( int e = 0; e < val; e++ )
				{
					if( min_hops[theVertices[v1].edges[e]]+1 <
					    min_hops[v1] )
					{
						min_hops[v1] = min_hops[theVertices[v1].edges[e]]+1;
						done = 0;
					}
				}
			}
		}
	
		memcpy( hop_mat + v * nv, min_hops, sizeof(int) * nv ); 

	}

	int max_hop = 0;
	for( int x = 0; x < nv; x++ )
	for( int y = 0; y < nv; y++ )
	{
		if( hop_mat[x*nv+y] > max_hop ) max_hop = hop_mat[x*nv+y];
		rms_node_hop[hop_mat[x*nv+y]] += effective_mass[x*nv+y] * effective_mass[x*nv+y];
		nrms_node_hop[hop_mat[x*nv+y]] += 1;	
	}

	for( int e = 0; e < max_hop; e++ )
		printf("%d %le\n", e, sqrt(rms_node_hop[e]/(nrms_node_hop[e]+1e-20)) );
	free(hop_mat);
	free(rms_node_hop);
	free(nrms_node_hop);
#endif

	
}

void surface::applyForceAtPoint( int f, double u, double v, double *dp, double *force_vector, force_set *theForceSet )
{
	double coefs[MAX_VALENCE+6];
	int coef_list[MAX_VALENCE+6];
	int ncoef=0;
	get_pt_coeffs( f, u, v, coefs, coef_list, &ncoef ); 

	for( int c1 = 0; c1 < ncoef; c1++ )
	{
		force_vector[3*coef_list[c1]+0] += coefs[c1] * dp[0];
		force_vector[3*coef_list[c1]+1] += coefs[c1] * dp[1];
		force_vector[3*coef_list[c1]+2] += coefs[c1] * dp[2];
	}
	
}


void surface::testFAP( int f, double u, double v, double *dp, double *effM )
{
	double frc[3*nv];
	
	memset( frc, 0, sizeof(double) * 3 * nv );
	
	double coefs[MAX_VALENCE+6];
	int coef_list[MAX_VALENCE+6];
	int ncoef=0;
	get_pt_coeffs( f, u, v, coefs, coef_list, &ncoef ); 

	for( int c1 = 0; c1 < ncoef; c1++ )
	{
		frc[3*coef_list[c1]+0] += coefs[c1] * dp[0];
		frc[3*coef_list[c1]+1] += coefs[c1] * dp[1];
		frc[3*coef_list[c1]+2] += coefs[c1] * dp[2];
	}

	double dv_check[3*nv];
	memset( dv_check, 0, sizeof(double) * 3 * nv );

	for( int v1 = 0; v1 < nv; v1++ )
	for( int v2 = 0; v2 < nv; v2++ )
	{
		dv_check[3*v1+0] += effM[v1*nv+v2] * frc[v2*3+0];
		dv_check[3*v1+1] += effM[v1*nv+v2] * frc[v2*3+1];
		dv_check[3*v1+2] += effM[v1*nv+v2] * frc[v2*3+2];
	}
	
	double probe[3] = {0,0,0};

	for( int c1 = 0; c1 < ncoef; c1++ )
	{
		probe[0] += dv_check[3*coef_list[c1]+0] * coefs[c1];
		probe[1] += dv_check[3*coef_list[c1]+1] * coefs[c1];
		probe[2] += dv_check[3*coef_list[c1]+2] * coefs[c1];
	}

	printf("FRC: %le %le %le RES %le %le %le\n", dp[0], dp[1], dp[2], probe[0], probe[1], probe[2] );
}


void surface::velocityAtPoint( int f, double u, double v, double *vp, double *vel_out )
{
	double coefs[MAX_VALENCE+6];
	int coef_list[MAX_VALENCE+6];
	int ncoef=0;
	get_pt_coeffs( f, u, v, coefs, coef_list, &ncoef ); 

	vel_out[0] = 0;
	vel_out[1] = 0;
	vel_out[2] = 0;

	for( int c1 = 0; c1 < ncoef; c1++ )
	{
		vel_out[0] += coefs[c1] * vp[3*coef_list[c1]+0];
		vel_out[1] += coefs[c1] * vp[3*coef_list[c1]+1];
		vel_out[2] += coefs[c1] * vp[3*coef_list[c1]+2];
	}
}

void clearForceSet( force_set *theForceSet )
{
	free(theForceSet->face_set);
	free(theForceSet->uv_set);
	free(theForceSet->frc_inverse);
	free(theForceSet->frc_coef);
	free(theForceSet->mass);
	free(theForceSet->frc_coef_list);
	free(theForceSet->frc_ncoef);
}

void surface::debugDeformation( double *r )
{
	double AAPROX = PBC_vec[2][2]*PBC_vec[2][2];
	extern double kc;

	// 2 to 5, interesting
	for( int m = 1; m <= 1; m++ )
	{
		double q = 2 * m * M_PI / PBC_vec[2][2];
//		double mag = sqrt(1.0 / (q*q*q*q*kc*AAPROX)) * 2*(rand()-0.5) / (double)RAND_MAX;
		double mag = 0.0;
		for( int v = 0; v < nv; v++ )
		{
			double z = r[3*v+2];
			double rxy = sqrt(r[3*v+0]*r[3*v+0]+r[3*v+1]*r[3*v+1]);
			printf("rxy: %le\n", rxy );
			double th = atan2( r[3*v+1], r[3*v+0] );
	
		r[3*v+0] *= (1 + mag * sin(2*M_PI*z*m/PBC_vec[2][2])/rxy);
		r[3*v+1] *= (1 + mag * sin(2*M_PI*z*m/PBC_vec[2][2])/rxy);
//		r[3*v+0] *= (1 + mag * sin(th)/rxy);
//		r[3*v+1] *= (1 + mag * cos(th)/rxy);
	
			double tan_disp[3] = { r[3*v+1]/rxy, -r[3*v+0]/rxy, 0 };

//		r[3*v+0] += mag * 2*(rand()-0.5)/(double)RAND_MAX * tan_disp[0];
//		r[3*v+1] += mag * 2*(rand()-0.5)/(double)RAND_MAX * tan_disp[1];
//		r[3*v+2] += mag * 2*(rand()-0.5)/(double)RAND_MAX;

//		r[3*v+0] *= (1 + mag * 2*(rand()-0.5)/(double)RAND_MAX/rxy);
//		r[3*v+1] *= (1 + mag * 2*(rand()-0.5)/(double)RAND_MAX/rxy);
		}
	}	
}


