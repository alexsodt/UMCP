#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "lapack_we_use.h"
#include "irr_kernel.h"
//#include "interp.h"


void irr_kernel::setup( int valence, int max_domain )
{
	val = valence;
	ndomains = max_domain;
	int high_ev;

	double w = 3.0/(8*val);
	double w6 = 3.0/(8*6);
	ncoords_base  = 1 + val + 5; 
	// assemble the matrix:

	int ncoords_extra = 6;

	double *A = (double *)malloc( sizeof(double) * (ncoords_base+ncoords_extra) * ncoords_base );
	memset( A, 0, sizeof(double) * (ncoords_base+ncoords_extra)*ncoords_base );

	// point 0 (i), a carried vertex. 
	A[0*ncoords_base+0] = 1-val*w;
	for( int p = 1; p < 1 + val; p++ )
		A[0*ncoords_base+p] = w;
	
	A[1*ncoords_base+0]   = (3.0/8.0);	
	A[1*ncoords_base+1]   = (3.0/8.0);	
	A[1*ncoords_base+2]   = (1.0/8.0);	
	A[1*ncoords_base+val] = (1.0/8.0);	

	for( int p = 2; p < val; p++ )
	{
		A[p*ncoords_base+0]   = (3.0/8.0);
		A[p*ncoords_base+p]   = (3.0/8.0);
		A[p*ncoords_base+p-1] = (1.0/8.0);
		A[p*ncoords_base+p+1] = (1.0/8.0);
	}
	
	A[val*ncoords_base+0]     = (3.0/8.0);
	A[val*ncoords_base+val]   = (3.0/8.0);
	A[val*ncoords_base+val-1] = (1.0/8.0);
	A[val*ncoords_base+1]     = (1.0/8.0);


	A[(val+1)*ncoords_base+1] = (3.0/8.0);	
	A[(val+1)*ncoords_base+val] = (3.0/8.0);	
	A[(val+1)*ncoords_base+0] = (1.0/8.0);	
	A[(val+1)*ncoords_base+val+1] = (1.0/8.0);	

	// point 1 (k)			
	A[(val+2)*ncoords_base+1] = 1-6*w6;
		A[(val+2)*ncoords_base+(0)] = w6;
		A[(val+2)*ncoords_base+(val)] = w6;
		A[(val+2)*ncoords_base+(val+1)] = w6;
		A[(val+2)*ncoords_base+(val+2)] = w6;
		A[(val+2)*ncoords_base+(val+3)] = w6;
		A[(val+2)*ncoords_base+(2)] = w6;
	
	A[(val+3)*ncoords_base+1] = (3.0/8.0);	
	A[(val+3)*ncoords_base+2] = (3.0/8.0);	
	A[(val+3)*ncoords_base+0] = (1.0/8.0);	
	A[(val+3)*ncoords_base+(val+3)] = (1.0/8.0);	
	
	// point 2 (j)
	A[(val+4)*ncoords_base+2] = 1-6*w6;
		A[(val+4)*ncoords_base+(0)] = w6;
		A[(val+4)*ncoords_base+(1)] = w6;
		A[(val+4)*ncoords_base+(val+3)] = w6;
		A[(val+4)*ncoords_base+(val+4)] = w6;
		A[(val+4)*ncoords_base+(val+5)] = w6;
		A[(val+4)*ncoords_base+(3)] = w6;
	
	A[(val+5)*ncoords_base+2]       = (3.0/8.0);	
	A[(val+5)*ncoords_base+3]       = (3.0/8.0);	
	A[(val+5)*ncoords_base+0]       = (1.0/8.0);	
	A[(val+5)*ncoords_base+(val+5)] = (1.0/8.0);	


	// now, the extra coordinates.

	A[(val+6)*ncoords_base+1] = (3.0/8.0);
	A[(val+6)*ncoords_base+(val+1)] = (3.0/8.0); 
	A[(val+6)*ncoords_base+(val+2)] = (1.0/8.0); 
	A[(val+6)*ncoords_base+val] = (1.0/8.0); 
	
	A[(val+7)*ncoords_base+1] = (3.0/8.0);
	A[(val+7)*ncoords_base+(val+2)] = (3.0/8.0); 
	A[(val+7)*ncoords_base+(val+3)] = (1.0/8.0); 
	A[(val+7)*ncoords_base+(val+1)] = (1.0/8.0); 
	
	A[(val+8)*ncoords_base+1] = (3.0/8.0);
	A[(val+8)*ncoords_base+(val+3)] = (3.0/8.0); 
	A[(val+8)*ncoords_base+2] = (1.0/8.0); 
	A[(val+8)*ncoords_base+(val+2)] = (1.0/8.0); 
	
	A[(val+9)*ncoords_base+2] = (3.0/8.0);
	A[(val+9)*ncoords_base+(val+3)] = (3.0/8.0); 
	A[(val+9)*ncoords_base+1] = (1.0/8.0); 
	A[(val+9)*ncoords_base+(val+4)] = (1.0/8.0); 
	
	A[(val+10)*ncoords_base+2] = (3.0/8.0);
	A[(val+10)*ncoords_base+(val+4)] = (3.0/8.0); 
	A[(val+10)*ncoords_base+(val+5)] = (1.0/8.0); 
	A[(val+10)*ncoords_base+(val+3)] = (1.0/8.0); 
	
	A[(val+11)*ncoords_base+2] = (3.0/8.0);
	A[(val+11)*ncoords_base+(val+5)] = (3.0/8.0); 
	A[(val+11)*ncoords_base+3] = (1.0/8.0); 
	A[(val+11)*ncoords_base+(val+4)] = (1.0/8.0); 
		
	double *Aevec = (double *)malloc( sizeof(double) * ncoords_base * ncoords_base );
	double *Aevec_R = (double *)malloc( sizeof(double) * ncoords_base * ncoords_base );
	double *Aval = (double *)malloc( sizeof(double) * ncoords_base );
	double *a_copy = (double*)malloc(sizeof(double) * ncoords_base * ncoords_base );
	char jobvl = 'V';
	char jobvr = 'N';
	int N = ncoords_base;
	int LDA = N;
	double WR[ncoords_base];
	double WI[ncoords_base];
	int LDVL = N;
	int LDVR = N;
	int LWORK = 8*N + N*N;
	double *work = (double *)malloc( sizeof(double) * LWORK );
	int info;
	double *VR = NULL;
	

	for( int i = 0; i < ncoords_base; i++ )
	for( int j = 0; j < ncoords_base; j++ )
		a_copy[j*ncoords_base+i] = A[i*ncoords_base+j];

	dgeev( &jobvl, &jobvr, &N, a_copy, &LDA, WR, WI, Aevec, &LDVL, VR, &LDVR, work, &LWORK, &info );
	memcpy( Aevec_R, Aevec, sizeof(double) * ncoords_base * ncoords_base );
	int ipiv[N];
	dgetrf( &N, &N, Aevec_R, &N, ipiv, &info );
	dgetri( &N, Aevec_R, &N, ipiv, work, &LWORK, &info );
	memcpy( Aval, WR, sizeof(double) * N );


	for( int i = 0; i < ncoords_base; i++ )
	for( int j = 0; j < ncoords_base; j++ )
		a_copy[j*ncoords_base+i] = Aevec_R[i*ncoords_base+j];
	memcpy( Aevec_R, a_copy, sizeof(double) * ncoords_base * ncoords_base );

	double max_ev = -1e10;
	for( int x = 0; x < ncoords_base; x++ )
	{
		if( Aval[x] > max_ev )
		{
			max_ev = Aval[x];
			high_ev = x;
		}
	}


	if( info != 0 )
	{
		printf("DGEEV error valence %d.\n", val);
		exit(1);
	}

	Aval[high_ev] = 1;

	free(work);
	free(a_copy);

	// each irregular vertex will need a list.
/*
	for( int i = 0; i < nv; i++ )
	{
		int val = theVertices[i].valence;

		for( int e = 0; e < val; e++ )
		{
			int j = theVertices[i].edges[e];
			int ep1 = e+1;
			if( ep1 >= val ) ep1 -= val;

			int k = theVertices[i].edges[ep1];
				
			int valk = theVertices[k].valence;
			int valj = theVertices[j].valence;

			if( j < i || k < i )
				continue;
			
			int ej = theVertices[i].edge_rev[e];
			int ek = theVertices[i].edge_rev[ep1];

			if( val == 6 && ( valj != 6 || valk != 6 ) )
			{
				printf("Illegal mesh: must be subdivided once isolating low index on irregular mesh.\n");
				exit(1);
			}

			int ncoords_base  = 1 + val + 5; 
			int ncoords_extra = 6;
			// we also need to "bring along" coordinates for the other three triangular splines in the final pass.
			// but we only need to do the eigenvector analysis on the first pass.
			int point_list[ncoords_base+ncoords_extra];
					
			int em1 = e-1; if( em1 < 0 ) em1 += val;
			int em2 = e-2; if( em2 < 0 ) em2 += val;
			int em3 = e-3; if( em3 < 0 ) em3 += val;
			int em4 = e-4; if( em4 < 0 ) em4 += val;
			int em5 = e-5; if( em5 < 0 ) em5 += val;
			int em6 = e-6; if( em6 < 0 ) em6 += val;
			// ep1 already defined.
			int ep2 = e+2; if( ep2 >= val ) ep2 -= val;
			int ep3 = e+3; if( ep3 >= val ) ep3 -= val;
			
			int ejm1 = ej-1; if( ejm1 < 0 ) ejm1 += valj;
			int ejm2 = ej-2; if( ejm2 < 0 ) ejm2 += valj;

			int ejp2 = ej+2; if( ejp2 >= valj ) ejp2 -= valj;
			int ejp3 = ej+3; if( ejp3 >= valj ) ejp3 -= valj;
			
			int ekm1 = ek-1; if( ekm1 < 0 ) ekm1 += valk;
			int ekm2 = ek-2; if( ekm2 < 0 ) ekm2 += valk;
			int ekp2 = ek+2; if( ekp2 >= valk ) ekp2 -= valk;
			int ekp3 = ek+3; if( ekp3 >= valk ) ekp3 -= valk;

			point_list[0] = i;
			point_list[1] = k; // ep1
			point_list[2] = j; // e
			point_list[3] = theVertices[i].edges[em1];
			point_list[4] = theVertices[i].edges[em2];
			if( val >= 5 ) point_list[5] = theVertices[i].edges[em3];
			if( val >= 6 ) point_list[6] = theVertices[i].edges[em4];
			if( val >= 7 ) point_list[7] = theVertices[i].edges[em5];
			if( val >= 8 ) point_list[8] = theVertices[i].edges[em6];

			
			point_list[val+1] = theVertices[k].edges[ekm2];
			point_list[val+2] = theVertices[k].edges[ekp3];
			point_list[val+3] = theVertices[k].edges[ekp2];
			
			point_list[val+4] = theVertices[j].edges[ejp3];
			point_list[val+5] = theVertices[j].edges[ejp2];
*/
	// store this locally then copy to map.
	double coord_map[12 * ncoords_base];
	memset( coord_map, 0, sizeof(double) * 12 * ncoords_base );

	map = (double *)malloc( sizeof(double ) * 12 * ncoords_base * 3 * max_domain );

	for( int domain = 1; domain <= max_domain; domain++ )
	{
		for( int subd = 0; subd < 3; subd++ )
		{
			double pow2 = pow( 2.0, domain-1.0 );
			
//			scale_f_u *= pow2;
//			scale_f_v *= pow2;
	
//			fu *= pow2;
//			fv *= pow2;
	
			int pcycle[12];
	
			if( subd == 0 ) // fu > 0.5
			{
//				fv = 2 * fv;
//				fu = 2.0 * fu-1.0;
			
//				scale_f_u *= 2;
//				scale_f_v *= 2;
	
				pcycle[0] = 2;
				pcycle[1] = val+3;
				pcycle[2] = val+4;
				pcycle[3] = val+5;
				pcycle[4] = 3;
				pcycle[5] = 0;
				pcycle[6] = 1;
				pcycle[7] = val+2;
				pcycle[8] = val+8;
				pcycle[9] = val+9;
				pcycle[10] =val+10;
				pcycle[11] = val+11;
			
			}
			else if( subd == 1 ) // fv > 0.5 )
			{
//				fu = 2 * fu;
//				fv = 2.0 * fv-1.0;
	
//				scale_f_u *= 2;
//				scale_f_v *= 2;
			
				
				pcycle[0] = 1;
				pcycle[1] = val+2;
				pcycle[2] = val+3;
				pcycle[3] = 2;
				pcycle[4] = 0;
				pcycle[5] = val;
				pcycle[6] = val+1;
				pcycle[7] = val+6;
				pcycle[8] = val+7;
				pcycle[9] = val+8;
				pcycle[10] =  val+9;
				pcycle[11] =  val+4;
			}
			else // subd == 2
			{
//				fv = 1 - 2 * fv;
//				fu = 1 - 2 * fu;
				
//				scale_f_u *= -2;
//				scale_f_v *= -2;
				
				pcycle[0] = val+3; 
				pcycle[1] = 2;
				pcycle[2] = 1;
				pcycle[3] = val+2;
				pcycle[4] = val+8;
				pcycle[5] = val+9;
				pcycle[6] = val+4;
				pcycle[7] = val+5;
				pcycle[8] = 3;
				pcycle[9] = 0;
				pcycle[10] = val;
				pcycle[11] = val+1;
			}
				
			double A_prev[ncoords_base*ncoords_base];	
			
			for( int e = 0; e < ncoords_base; e++ )
			for( int j = 0; j < ncoords_base; j++ )
			{
				A_prev[e*ncoords_base+j] = 0;
	
				for( int i = 0; i < ncoords_base; i++ )
					A_prev[e*ncoords_base+j] += Aevec[i*ncoords_base+j] * pow( Aval[i], domain-1 ) * Aevec_R[i*ncoords_base+e];
			}
			
			double A_use[(ncoords_base+ncoords_extra)*ncoords_base];
			memset( A_use, 0, sizeof(double) * (ncoords_base+ncoords_extra)*ncoords_base );
			for( int f = 0; f < ncoords_base+ncoords_extra; f++ )
			{
				for( int i = 0; i < ncoords_base; i++ )
				{
					A_use[f*ncoords_base+i] = 0;
					
					for( int e = 0; e < ncoords_base; e++ )
						A_use[f*ncoords_base+i] += A[f*ncoords_base+e] * A_prev[e*ncoords_base+i];
				}
			}
	
			for( int y = 0; y < 12; y++ )
			{
				double sum = 0;
				
				for( int x = 0; x < ncoords_base; x++ )
				{
					// y, the spline base, x, the coordinate vertex base.
					sum += A_use[pcycle[y]*ncoords_base+x];
				}	
	
	
				for( int x = 0; x < ncoords_base; x++ )
				{
					// y, the spline base, x, the coordinate vertex base.
					coord_map[y*ncoords_base+x] = A_use[pcycle[y]*ncoords_base+x];
				}	
			}

			memcpy( map + (domain-1)* 3 * (12*ncoords_base) + subd * (12 * ncoords_base), coord_map, sizeof(double) * ncoords_base * 12 );
		}
	}
}

int irr_kernel::domain( double fu, double fv )
{
	double f_domain = -log10(fu+fv)/log10(2.0);

	if( !(f_domain < 1e300 ) || !(f_domain > -1e300) )
	{
		return -1;
	}
	if( fu+fv >= 1.0 -1e-9 )
		f_domain = 1;
	int domain = lround(ceil(f_domain));

	return domain;
}

void irr_kernel::get_map_transform( double *u_u, double *u_v, double *v_u, double *v_v )
{
	double fu = *u_u;
	double fv = *v_v;
	double f_domain = -log10(fu+fv)/log10(2.0);
	if( !(f_domain < 1e300 ) && !(f_domain > -1e300) )
	{
		*u_u=1;
		*v_v=1;
		return;
	}

	if( fu+fv >= 1.0 -1e-9 )
		f_domain = 1;
	int domain = lround(ceil(f_domain));
	
	if( domain-1 >= ndomains )
	{
		printf("ERROR: requested domain not provided by kernel.\n");
		exit(1);
	} 

	double pow2 = pow( 2.0, domain-1.0 );
	
	double scale_f_u = 1.0;
	double scale_f_v = 1.0;

	scale_f_u *= pow2;
	scale_f_v *= pow2;

	fu *= pow2;
	fv *= pow2;
		
	int subd = 0;
					
	if( fu > 0.5 )
	{
		subd = 0;
		fv = 2 * fv;
		fu = 2.0 * fu-1.0;

		scale_f_u *= 2.0;
		scale_f_v *= 2.0;

		*u_u = scale_f_u;
		*u_v = 0;
		*v_v = 0;
		*v_v = scale_f_v;
	}
	else if( fv > 0.5 )
	{
		subd = 1;
		fu = 2 * fu;
		fv = 2.0 * fv-1.0;

		scale_f_u *= 2;
		scale_f_v *= 2;
		
		*u_u = scale_f_u;
		*u_v = 0;
		*v_v = 0;
		*v_v = scale_f_v;
	}
	else
	{
		subd = 2;
		fv = 1 - 2 * fv;
		fu = 1 - 2 * fu;
		
		scale_f_u *= -2;
		scale_f_v *= -2;
		
		*u_u = scale_f_u;
		*u_v = 0;
		*v_v = 0;
		*v_v = scale_f_v;
	}
}

double *irr_kernel::get_map( double *u, double *v )
{
	double fu = *u;
	double fv = *v;
	double f_domain = -log10(fu+fv)/log10(2.0);
	if( !(f_domain < 1e300 ) && !(f_domain > -1e300) )
	{
		return NULL;
	}

	if( fu+fv >= 1.0 -1e-9 )
		f_domain = 1;
	int domain = lround(ceil(f_domain));
	
	if( domain-1 >= ndomains )
	{
		printf("ERROR: requested domain not provided by kernel.\n");
		exit(1);
	} 

	double pow2 = pow( 2.0, domain-1.0 );
	
	double scale_f_u = 1.0;
	double scale_f_v = 1.0;

	scale_f_u *= pow2;
	scale_f_v *= pow2;

	fu *= pow2;
	fv *= pow2;
		
	int subd = 0;
					
	if( fu > 0.5 )
	{
		subd = 0;
		fv = 2 * fv;
		fu = 2.0 * fu-1.0;

		scale_f_u *= 2.0;
		scale_f_v *= 2.0;
	}
	else if( fv > 0.5 )
	{
		subd = 1;
		fu = 2 * fu;
		fv = 2.0 * fv-1.0;

		scale_f_u *= 2;
		scale_f_v *= 2;
	}
	else
	{
		subd = 2;
		fv = 1 - 2 * fv;
		fu = 1 - 2 * fu;
		
		scale_f_u *= -2;
		scale_f_v *= -2;
	}

	*u = fu;
	*v = fv;

	return map + ((domain-1) * 3 +subd )* 12 * ncoords_base;  
	
}

