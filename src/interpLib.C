#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include <sys/time.h>
#include "interp.h"
#include "mutil.h"
#include "l-bfgs.h"
#include "simpleRandom.C"	
#include "m_triangles.h"
#include "config.h"
#include "gsl/gsl_sf.h"
#include "gauss.h"
#include "lapack_we_use.h"
#include "Bin.h"
#include "Point.h"
#include "fast_mm.h"
#include "parallel.h"


#define ATLAS

#define THRESH (1e-4)

//#define LOW_RES
#define MAX_INV_VALENCE 12

#define USE_G0
double eps_f_min = 5;

#ifdef LINUX
#include "clapack.h"
#endif
	
	double fix_edge_k = 0.0;
	double fix_edge_r0 = 1.0;
//#define USE_GPU
int activate_height_debug = 0;

	double *oarr = NULL;
#define FN_1
//#define FN_2


#ifdef LOW_RES
int n_v_finite = 1;
int n_u_gauss = 1;
#else
int n_v_finite = 4;
int n_u_gauss = 4;

#endif
int debug_bit = 0;
int do_face = 375;
int do_pt = 2;
int cnx = 0;
int cny = 0;
int cnx2 = 0;
int cny2 = 0;

double VPEN = 0;
double mag1 = 0;
double mag2 = 0;
double phase1 = 0;
double phase2 = 0;

//#define DEBUG_A
#ifdef LOW_RES
#define G_Q_P_1
#else
#define G_Q_P_3
//#define G_Q_P_15
#endif
//#define G_Q_P_3
//#define G_Q_P_4
//#define G_Q_P_PRINT
//#define G_Q_DEBUG
//#define G_Q_DEBUG_2
//#define DEBUG_C1
//#define DEBUG_NRMX
//#define PRINT_A
//#define AMAT_DEBUG

double cos_constraint = 1000000;
double tilt_modulus = 0.2 * 0.592;
double k_reg = 0;//1; // regularization elastic constant.
double kc = 14 ;;//0.01; // 14;//14;
double kc_irr = 14 ;//0.01; //14;//14;
double kg = 0;//-0.7*14;//-0.7*14 ;//-0.007; //-7;//-7;
double water_KV = 0.0056;
double vol0 = 508872847.505142*0.65;
double KA  = 0.0; //0.215;//0.215; // monolayer, doubled later
double KA4 = 0.0;
double micro_KA = 0.0; //0.215; // monolayer, doubled later
int do_print_area_ps = 0;
extern double KV_scale;
double Power( double a, double b);
double Sqrt( double a );
double tilt_cons = 1000;
double tilt_scale  = 1.0;
int nearMembraneWorker( double *r1, int nv1, double *r2, double **M, int mlow, int mhigh, int level, int max_level, double radius, 
		double factor, 
		double cur_u, double cur_v,
		double *puv_list, double *areas, int *npts );	

double ext_com_k = 0;
double force_face_0 = 0;
double ext_force_vec[3] = { 1/sqrt(3), 1/sqrt(3),1/sqrt(3) };


int compare_sense( int code1, int code2 )
{
	int code_matrix[3][3] =
	{
		{ -1, 1, -1 },
		{ 1, -1, 1  },
		{ -1, 1, -1 }
	};
		 
	return code_matrix[code1][code2];
}

void surface::generateVolumePlan( void )
{
	double *Aevec[1+MAX_INV_VALENCE];
	double *Aevec_R[1+MAX_INV_VALENCE];
	double *Aval[1+MAX_INV_VALENCE];
	double *Atot[1+MAX_INV_VALENCE];
	int high_ev[1+MAX_INV_VALENCE];			
	
	for( int val = 4; val <= MAX_INV_VALENCE; val++ )
	{
		double w = 3.0/(8*val);
		double w6 = 3.0/(8*6);
		int ncoords_base  = 1 + val + 5; 
		// assemble the matrix:

		int ncoords_extra = 6;
		Atot[val] = (double *)malloc( sizeof(double) * (ncoords_base+ncoords_extra) * ncoords_base ); 
	
		double *A = Atot[val];

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
			
		Aevec[val] = (double *)malloc( sizeof(double) * ncoords_base * ncoords_base );
		Aevec_R[val] = (double *)malloc( sizeof(double) * ncoords_base * ncoords_base );
		Aval[val] = (double *)malloc( sizeof(double) * ncoords_base );
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

		dgeev( &jobvl, &jobvr, &N, a_copy, &LDA, WR, WI, Aevec[val], &LDVL, VR, &LDVR, work, &LWORK, &info );
		memcpy( Aevec_R[val], Aevec[val], sizeof(double) * ncoords_base * ncoords_base );
		int ipiv[N];
		dgetrf( &N, &N, Aevec_R[val], &N, ipiv, &info );
		dgetri( &N, Aevec_R[val], &N, ipiv, work, &LWORK, &info );

		memcpy( Aval[val], WR, sizeof(double) * N );


		for( int i = 0; i < ncoords_base; i++ )
		for( int j = 0; j < ncoords_base; j++ )
			a_copy[j*ncoords_base+i] = Aevec_R[val][i*ncoords_base+j];
		memcpy( Aevec_R[val], a_copy, sizeof(double) * ncoords_base * ncoords_base );

		double max_ev = -1e10;
		for( int x = 0; x < ncoords_base; x++ )
		{
			if( Aval[val][x] > max_ev )
			{
				max_ev = Aval[val][x];
				high_ev[val] = x;
			}
		}


		if( info != 0 )
		{
			printf("DGEEV error valence %d.\n", val);
			exit(1);
		}
		else
		{
//#define EV_STRUCTURE_DEBUG
#ifdef EV_STRUCTURE_DEBUG
			printf("valence: %d\n", val );
			for( int i = 0; i < ncoords_base; i++ )
				printf("\teval %d %.14le %.14le\n", i, WR[i], WI[i] );
			printf("high ev vector:\n");

			for( int i = 0; i < ncoords_base; i++ )
				printf("\t%d %le\n", i, Aevec[val][high_ev[val]*ncoords_base+i] * Aevec_R[val][high_ev[val]*ncoords_base+i] );
#endif
		}

		Aval[val][high_ev[val]] = 1;

		free(work);
		free(a_copy);
	}

	int ni = 12;
	
	nfaces = 0;

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

			nfaces++;
		}
	}
	
	theVolumeFormulas = (volume_formula *)malloc( sizeof(volume_formula) * nfaces  );

	int trigger=0;

	int nv_space = ni * (ni+1)/2;

	int f = 0;
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
			printf("%d valence %d, %d valence %d, %d valence %d.\n",
				i, theVertices[i].valence,
				j, theVertices[j].valence,
				k, theVertices[k].valence );
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

			theVolumeFormulas[f].ni = ni;
			theVolumeFormulas[f].ncoor = ncoords_base;
			theVolumeFormulas[f].cp = (int *)malloc( sizeof(int) * ncoords_base ); 

			theVolumeFormulas[f].r_pbc = (double *)malloc( sizeof(double) * 3 * ncoords_base * nv_space );
			theVolumeFormulas[f].r_w = (double *)malloc( sizeof(double) * ncoords_base * nv_space ); 

			memcpy( theVolumeFormulas[f].cp, point_list, sizeof(int) * ncoords_base );
			memset( theVolumeFormulas[f].r_pbc, 0, sizeof(double) * 3 * ncoords_base );

			int ioff = 0;
			for( int iu = 0; iu < ni; iu++ )
			for( int iv = 0; iv < ni-iu; iv++, ioff++ )
			{
				double fu = iu / (double)(ni-1);
				double fv = iv / (double)(ni-1);
				double fw = 1-fu-fv;

				for( int p = 0; p < ncoords_base; p++ )
				{
					int *cp = point_list;

					double *r1 = theVertices[cp[p]].r;
					double *r0 = theVertices[cp[0]].r;
	
					double dr[3] = { r1[0] - r0[0], r1[1] - r0[1], r1[2] - r0[2] };
					double add[3]={0,0,0};

					MinImage3D( dr, PBC_vec, add );

					theVolumeFormulas[f].r_pbc[ioff*3*ncoords_base+3*p+0] = add[0]*PBC_vec[0][0] + add[1] * PBC_vec[1][0] + add[2] * PBC_vec[2][0];
					theVolumeFormulas[f].r_pbc[ioff*3*ncoords_base+3*p+1] = add[0]*PBC_vec[0][1] + add[1] * PBC_vec[1][1] + add[2] * PBC_vec[2][1];
					theVolumeFormulas[f].r_pbc[ioff*3*ncoords_base+3*p+2] = add[0]*PBC_vec[0][2] + add[1] * PBC_vec[1][2] + add[2] * PBC_vec[2][2];
				}

				double scale_f_u = 1.0;
				double scale_f_v = 1.0;

				if( iu == 0 && iv == 0 )
				{
					// the highest eigenvector is just copied in here.
					for( int x = 0; x < ncoords_base; x++ )
						theVolumeFormulas[f].r_w[ioff*ncoords_base+x] = Aevec[val][high_ev[val]*ncoords_base+x] * Aevec_R[val][high_ev[val]*ncoords_base+x];
				}
				else
				{
					double coord_map[12 * ncoords_base];
					memset( coord_map, 0, sizeof(double) * 12 * ncoords_base );

					if( val != 6 )
					{

						double f_domain = -log10(fu+fv)/log10(2.0);
						if( fu+fv >= 1.0 -1e-9 )
							f_domain = 1;
						int domain = lround(ceil(f_domain));
						double pow2 = pow( 2.0, domain-1.0 );
						
						scale_f_u *= pow2;
						scale_f_v *= pow2;

						fu *= pow2;
						fv *= pow2;

						int pcycle[12];

						if( fu > 0.5 )
						{
							fv = 2 * fv;
							fu = 2.0 * fu-1.0;
						
							scale_f_u *= 2;
							scale_f_v *= 2;

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
						else if( fv > 0.5 )
						{
							fu = 2 * fu;
							fv = 2.0 * fv-1.0;

							scale_f_u *= 2;
							scale_f_v *= 2;
						
							
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
						else
						{
							fv = 1 - 2 * fv;
							fu = 1 - 2 * fu;
							
							scale_f_u *= -2;
							scale_f_v *= -2;
							
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
								A_prev[e*ncoords_base+j] += Aevec[val][i*ncoords_base+j] * pow( Aval[val][i], domain-1 ) * Aevec_R[val][i*ncoords_base+e];
						}
						
						

						double A_use[(ncoords_base+ncoords_extra)*ncoords_base];
						memset( A_use, 0, sizeof(double) * (ncoords_base+ncoords_extra)*ncoords_base );
						for( int f = 0; f < ncoords_base+ncoords_extra; f++ )
						{
							for( int i = 0; i < ncoords_base; i++ )
							{
								A_use[f*ncoords_base+i] = 0;
								
								for( int e = 0; e < ncoords_base; e++ )
									A_use[f*ncoords_base+i] += Atot[val][f*ncoords_base+e] * A_prev[e*ncoords_base+i];
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
					}
					else
					{
						for( int i = 0; i < 12; i++ )
							coord_map[i*12+i] = 1;	
					}

					// get the coefficients.
						
						double u = fu;
						double v = fv;
						double w = 1 - u - v;
	
						double u2 = u*u;
						double u3 = u*u*u;
						double u4 = u*u*u*u;
						
						double v2 = v*v;
						double v3 = v*v*v;
						double v4 = v*v*v*v;
				
						double w2 = w*w;
						double w3 = w*w*w;
						double w4 = w*w*w*w;
						
					
					// 8 : 0
					// 7 : 1
					// 4 : 2
					// 5 : 3
					// 9 : 4
					// 12 : 5
					// 11 : 6
					// 10 : 7
					// 6 : 8
					// 3 : 9
					// 1 : 10
					// 2 : 11

	

					double n1 = (1.0/12.0)*(u4+2*u3*v); 
					double n2 = (1.0/12.0)*(u4+2*u3*w); 
					double n3 = (1.0/12.0)*(u4+2*u3*w+6*u3*v+6*u2*v*w+12*u2*v2+6*u*v2*w+6*u*v3+2*v3*w+v4); 
					double n4 = (1.0/12.0)*(6*u4+24*u3*w+24*u2*w2+8*u*w3+w4+24*u3*v+60*u2*v*w+36*u*v*w2+6*v*w3+24*u2*v2+36*u*v2*w+12*v2*w2+8*u*v3+6*v3*w+v4); 
					double n5 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+2*u3*v+6*u2*v*w+6*u*v*w2+2*v*w3); 
					double n6 = (1.0/12.0)*(2*u*v3+v4); 
					double n7 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+8*u3*v+36*u2*v*w+36*u*v*w2+8*v*w3+24*u2*v2+60*u*v2*w+24*v2*w2+24*u*v3+24*v3*w+6*v4); 
					double n8 = (1.0/12.0)*(u4+8*u3*w+24*u2*w2+24*u*w3+6*w4+6*u3*v+36*u2*v*w+60*u*v*w2+24*v*w3+12*u2*v2+36*u*v2*w+24*v2*w2+6*u*v3+8*v3*w+v4); 
					double n9 = (1.0/12.0)*(2*u*w3+w4); 
					double n10 = (1.0/12.0)*(2*v3*w+v4); 
					double n11 = (1.0/12.0)*(2*u*w3+w4+6*u*v*w2+6*v*w3+6*u*v2*w+12*v2*w2+2*u*v3+6*v3*w+v4); 
					double n12 = (1.0/12.0)*(w4+2*v*w3);

					double n1_alt = Power(u,4)/12. + (Power(u,3)*v)/6.;		
					double n2_alt = Power(u,4)/12. + (Power(u,3)*w)/6.;
					double n3_alt = Power(u,4)/12. + (Power(u,3)*v)/2. + Power(u,2)*Power(v,2) + (u*Power(v,3))/2. + Power(v,4)/12. + (Power(u,3)*w)/6. + (Power(u,2)*v*w)/2. + (u*Power(v,2)*w)/2. + (Power(v,3)*w)/6.;
					double n4_alt = Power(u,4)/2. + 2*Power(u,3)*v + 2*Power(u,2)*Power(v,2) + (2*u*Power(v,3))/3. + Power(v,4)/12. + 2*Power(u,3)*w + 5*Power(u,2)*v*w + 3*u*Power(v,2)*w + (Power(v,3)*w)/2. + 2*Power(u,2)*Power(w,2) + 3*u*v*Power(w,2) + Power(v,2)*Power(w,2) + (2*u*Power(w,3))/3. + (v*Power(w,3))/2. + Power(w,4)/12.;
					double n5_alt = Power(u,4)/12. + (Power(u,3)*v)/6. + (Power(u,3)*w)/2. + (Power(u,2)*v*w)/2. + Power(u,2)*Power(w,2) + (u*v*Power(w,2))/2. + (u*Power(w,3))/2. + (v*Power(w,3))/6. + Power(w,4)/12.;
					double n6_alt = (u*Power(v,3))/6. + Power(v,4)/12.;
					double n7_alt = Power(u,4)/12. + (2*Power(u,3)*v)/3. + 2*Power(u,2)*Power(v,2) + 2*u*Power(v,3) + Power(v,4)/2. + (Power(u,3)*w)/2. + 3*Power(u,2)*v*w + 5*u*Power(v,2)*w + 2*Power(v,3)*w + Power(u,2)*Power(w,2) + 3*u*v*Power(w,2) + 2*Power(v,2)*Power(w,2) + (u*Power(w,3))/2. + (2*v*Power(w,3))/3. + Power(w,4)/12.;
					double n8_alt = Power(u,4)/12. + (Power(u,3)*v)/2. + Power(u,2)*Power(v,2) + (u*Power(v,3))/2. + Power(v,4)/12. + (2*Power(u,3)*w)/3. + 3*Power(u,2)*v*w + 3*u*Power(v,2)*w + (2*Power(v,3)*w)/3. + 2*Power(u,2)*Power(w,2) + 5*u*v*Power(w,2) + 2*Power(v,2)*Power(w,2) + 2*u*Power(w,3) + 2*v*Power(w,3) + Power(w,4)/2.;
					double n9_alt = (u*Power(w,3))/6. + Power(w,4)/12.;
					double n10_alt = Power(v,4)/12. + (Power(v,3)*w)/6.;
					double n11_alt =(u*Power(v,3))/6. + Power(v,4)/12. + (u*Power(v,2)*w)/2. + (Power(v,3)*w)/2. + (u*v*Power(w,2))/2. + Power(v,2)*Power(w,2) + (u*Power(w,3))/6. + (v*Power(w,3))/2. + Power(w,4)/12.;
					double n12_alt = (v*Power(w,3))/6. + Power(w,4)/12.;
	
					double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };

					for( int x = 0; x < ncoords_base; x++ )
					{
						theVolumeFormulas[f].r_w[ioff*ncoords_base+x] = 0;
		
						for( int y = 0; y < 12; y++ )
							theVolumeFormulas[f].r_w[ioff*ncoords_base+x] += coord_map[y*ncoords_base+x] * ceff_map[y];

					}
				}
			}

			f++;
		}
	}
	
/*				
	FILE *theFile = fopen("plan.xyz", "w");
	fprintf(theFile, "%d\n", nf*2 );
	fprintf(theFile, "plan\n");

	for( int f = 0; f < nf; f++ )
	{
		double r[3] = {0,0,0};
		double ru[3] = {0,0,0};
		double rv[3] = {0,0,0};
		
		int np = theVolumeFormulas[f].ncoor;

		int *cp = theVolumeFormulas[f].cp;
		for( int p = 0; p < np; p++ )
		{
			r[0] += theVolumeFormulas[f].r_w[p] * (theVertices[cp[p]].r[0] + theVolumeFormulas[f].r_pbc[3*p+0]); 
			r[1] += theVolumeFormulas[f].r_w[p] * (theVertices[cp[p]].r[1] + theVolumeFormulas[f].r_pbc[3*p+1]); 
			r[2] += theVolumeFormulas[f].r_w[p] * (theVertices[cp[p]].r[2] + theVolumeFormulas[f].r_pbc[3*p+2]); 
			
			ru[0] += theVolumeFormulas[f].r_u[p] * (theVertices[cp[p]].r[0] + theVolumeFormulas[f].r_pbc[3*p+0]); 
			ru[1] += theVolumeFormulas[f].r_u[p] * (theVertices[cp[p]].r[1] + theVolumeFormulas[f].r_pbc[3*p+1]); 
			ru[2] += theVolumeFormulas[f].r_u[p] * (theVertices[cp[p]].r[2] + theVolumeFormulas[f].r_pbc[3*p+2]); 
			
			rv[0] += theVolumeFormulas[f].r_v[p] * (theVertices[cp[p]].r[0] + theVolumeFormulas[f].r_pbc[3*p+0]); 
			rv[1] += theVolumeFormulas[f].r_v[p] * (theVertices[cp[p]].r[1] + theVolumeFormulas[f].r_pbc[3*p+1]); 
			rv[2] += theVolumeFormulas[f].r_v[p] * (theVertices[cp[p]].r[2] + theVolumeFormulas[f].r_pbc[3*p+2]); 
		}

		normalize( ru );
		normalize( rv );

		double del = 0.1;
		double nrm[3];
		cross( ru, rv, nrm );
		normalize(nrm);
		fprintf(theFile, "C %lf %lf %lf\n", r[0], r[1], r[2] );
		fprintf(theFile, "O %lf %lf %lf\n", r[0]+del * nrm[0], r[1]+del * nrm[1], r[2]+del * nrm[2] );
	}
	fclose(theFile);
*/
}

double surface::dvolume( double *r, double *g, double scale )
{ // also returns volume.

	if( !theVolumeFormulas )
		generateVolumePlan();

	double vol = 0;
	
	for( int f = 0; f < nf_faces; f++ )
	{
		int ni = theVolumeFormulas[f].ni;

		int nv_space = ni*(ni+1)/2;

		double coords[3*ni*ni];

		int ioff=0;
		for( int iu = 0; iu < ni; iu++ )
		for( int iv = 0; iv < ni-iu; iv++, ioff++ )
		{
			double R[3] = {0,0,0};

			int *cp = theVolumeFormulas[f].cp;
			int np = theVolumeFormulas[f].ncoor;

			for( int p = 0; p < np; p++ )
			{
				R[0] += theVolumeFormulas[f].r_w[ioff*np+p] * (r[cp[p]*3+0] + theVolumeFormulas[f].r_pbc[ioff*np*3+3*p+0]); 
				R[1] += theVolumeFormulas[f].r_w[ioff*np+p] * (r[cp[p]*3+1] + theVolumeFormulas[f].r_pbc[ioff*np*3+3*p+1]); 
				R[2] += theVolumeFormulas[f].r_w[ioff*np+p] * (r[cp[p]*3+2] + theVolumeFormulas[f].r_pbc[ioff*np*3+3*p+2]); 
			}

			coords[(iu*ni+iv)*3+0] = R[0];
			coords[(iu*ni+iv)*3+1] = R[1];
			coords[(iu*ni+iv)*3+2] = R[2];
		}

		double dv_dr[(ni*ni)*3];
		memset( dv_dr, 0, sizeof(double) * ni*ni*3 );

		for( int i = 0; i < ni; i++ )
		for( int j = 0; j < ni-i; j++ )
		{
			for( int pass = 0; pass < 2; pass++ )
			{
				double *r1;
				double *r2;
				double *r3;
				int ir1, ir2, ir3;
				if( pass == 0 && (i + j + 1 < ni && i < ni-1 && j < ni-1) )
				{
					ir1 = i*ni+j;
					ir2 = (i+1)*ni+j;
					ir3 = i*ni+j+1;
					
					r1 = coords + ir1*3; 
					r2 = coords + ir2*3;
					r3 = coords + ir3*3;
				}
				else if( pass == 1 && (i-1 >= 0 && j+1 < ni && i+j+1 < ni) )
				{
					ir1 = i*ni+j;
					ir2 = i*ni+j+1;
					ir3 = (i-1)*ni+j+1;

					r1 = coords + ir1*3; 
					r2 = coords + ir2*3;
					r3 = coords + ir3*3;
				}
				else
					continue;
				

				double nrm[3];

				double dr1[3] = { r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]};
				double dr2[3] = { r3[0]-r1[0],r3[1]-r1[1],r3[2]-r1[2]};

				cross( dr1, dr2, nrm );

				double A = normalize(nrm)/2;

				double bc[3] = { (r1[0]+r2[0]+r3[0])/3,
						 (r1[1]+r2[1]+r3[1])/3,
						 (r1[2]+r2[2]+r3[2])/3 };

				vol += (1.0/3.0) * (bc[0] * nrm[0] + bc[1] * nrm[1] + bc[2] * nrm[2]) * A; 

				double r1x = r1[0];
				double r1y = r1[1];
				double r1z = r1[2];
				
				double r2x = r2[0];
				double r2y = r2[1];
				double r2z = r2[2];
				
				double r3x = r3[0];
				double r3y = r3[1];
				double r3z = r3[2];

				double d_nrmx_d_r1x, d_nrmx_d_r1y, d_nrmx_d_r1z;
				double d_nrmy_d_r1x, d_nrmy_d_r1y, d_nrmy_d_r1z;
				double d_nrmz_d_r1x, d_nrmz_d_r1y, d_nrmz_d_r1z;
				
				double d_nrmx_d_r2x, d_nrmx_d_r2y, d_nrmx_d_r2z;
				double d_nrmy_d_r2x, d_nrmy_d_r2y, d_nrmy_d_r2z;
				double d_nrmz_d_r2x, d_nrmz_d_r2y, d_nrmz_d_r2z;
				
				double d_nrmx_d_r3x, d_nrmx_d_r3y, d_nrmx_d_r3z;
				double d_nrmy_d_r3x, d_nrmy_d_r3y, d_nrmy_d_r3z;
				double d_nrmz_d_r3x, d_nrmz_d_r3y, d_nrmz_d_r3z;

				d_nrmx_d_r1x = -((-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)*(2*(r2y - r3y)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(-r2z + r3z)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5));
				d_nrmx_d_r1y = -((-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)*(2*(-r2x + r3x)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(r2z - r3z)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5)) + (r2z - r3z)/Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2));
				d_nrmx_d_r1z = -((-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)*(2*(r2x - r3x)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z) + 2*(-r2y + r3y)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5)) + (-r2y + r3y)/Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2));

				d_nrmy_d_r1x = -((r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z)*(2*(r2y - r3y)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(-r2z + r3z)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5)) + (-r2z + r3z)/Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2));
				d_nrmy_d_r1y = -((r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z)*(2*(-r2x + r3x)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(r2z - r3z)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5));
				d_nrmy_d_r1z = -((r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z)*(2*(r2x - r3x)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z) + 2*(-r2y + r3y)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5)) + (r2x - r3x)/Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2));
		
				d_nrmz_d_r1x = -((-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y)*(2*(r2y - r3y)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(-r2z + r3z)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5)) + (r2y - r3y)/Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2));

				d_nrmz_d_r1y = -((-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y)*(2*(-r2x + r3x)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(r2z - r3z)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5)) + (-r2x + r3x)/Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2));

				d_nrmz_d_r1z = -((-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y)*(2*(r2x - r3x)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z) + 2*(-r2y + r3y)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5));

				d_nrmx_d_r2x = -((-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)*(2*(-r1y + r3y)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(r1z - r3z)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5));
				d_nrmx_d_r2y = -((-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)*(2*(r1x - r3x)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(-r1z + r3z)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5)) + (-r1z + r3z)/Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2));
				d_nrmx_d_r2z = -((-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)*(2*(-r1x + r3x)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z) + 2*(r1y - r3y)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5)) + (r1y - r3y)/Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2));

				d_nrmy_d_r2x = -((r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z)*(2*(-r1y + r3y)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(r1z - r3z)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5)) + (r1z - r3z)/Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2));
				d_nrmy_d_r2y = -((r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z)*(2*(r1x - r3x)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(-r1z + r3z)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5));
				d_nrmy_d_r2z = -((r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z)*(2*(-r1x + r3x)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z) + 2*(r1y - r3y)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5)) + (-r1x + r3x)/Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2));
	
				d_nrmz_d_r2x = -((-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y)*(2*(-r1y + r3y)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(r1z - r3z)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5)) + (-r1y + r3y)/Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2));
				d_nrmz_d_r2y = -((-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y)*(2*(r1x - r3x)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(-r1z + r3z)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5)) + (r1x - r3x)/Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2));
				d_nrmz_d_r2z = -((-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y)*(2*(-r1x + r3x)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z) + 2*(r1y - r3y)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5));

				d_nrmx_d_r3x = -((-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)*(2*(r1y - r2y)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(-r1z + r2z)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5));
				d_nrmx_d_r3y = -((-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)*(2*(-r1x + r2x)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(r1z - r2z)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5)) + (r1z - r2z)/Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2));
				d_nrmx_d_r3z = -((-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)*(2*(r1x - r2x)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z) + 2*(-r1y + r2y)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5)) + (-r1y + r2y)/Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2));

				d_nrmy_d_r3x = -((r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z)*(2*(r1y - r2y)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(-r1z + r2z)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5)) + (-r1z + r2z)/Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2));
				d_nrmy_d_r3y = -((r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z)*(2*(-r1x + r2x)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(r1z - r2z)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5));
				d_nrmy_d_r3z = -((r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z)*(2*(r1x - r2x)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z) + 2*(-r1y + r2y)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5)) + (r1x - r2x)/Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2));

				d_nrmz_d_r3x = -((-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y)*(2*(r1y - r2y)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(-r1z + r2z)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5)) + (r1y - r2y)/Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2));
				d_nrmz_d_r3y = -((-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y)*(2*(-r1x + r2x)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(r1z - r2z)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5)) + (-r1x + r2x)/Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2));
				d_nrmz_d_r3z = -((-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y)*(2*(r1x - r2x)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z) + 2*(-r1y + r2y)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z)))/(2.*Power(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2),1.5));


				double d_A_d_r1x = (2*(r2y - r3y)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(-r2z + r3z)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z))/(4.*Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2)));
				double d_A_d_r1y = (2*(-r2x + r3x)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(r2z - r3z)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z))/(4.*Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2)));
				double d_A_d_r1z = (2*(r2x - r3x)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z) + 2*(-r2y + r3y)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z))/(4.*Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2)));
	
				double d_A_d_r2x = (2*(-r1y + r3y)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(r1z - r3z)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z))/(4.*Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2)));
				double d_A_d_r2y = (2*(r1x - r3x)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(-r1z + r3z)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z))/(4.*Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2)));
				double d_A_d_r2z = (2*(-r1x + r3x)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z) + 2*(r1y - r3y)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z))/(4.*Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2)));

				double d_A_d_r3x = (2*(r1y - r2y)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(-r1z + r2z)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z))/(4.*Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2)));
				double d_A_d_r3y = (2*(-r1x + r2x)*(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y) + 2*(r1z - r2z)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z))/(4.*Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2)));
				double d_A_d_r3z = (2*(r1x - r2x)*(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z) + 2*(-r1y + r2y)*(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z))/(4.*Sqrt(Power(-(r1y*r2x) + r1x*r2y + r1y*r3x - r2y*r3x - r1x*r3y + r2x*r3y,2) + Power(r1z*r2x - r1x*r2z - r1z*r3x + r2z*r3x + r1x*r3z - r2x*r3z,2) + Power(-(r1z*r2y) + r1y*r2z + r1z*r3y - r2z*r3y - r1y*r3z + r2y*r3z,2)));
			 				

				double d_v_d_nrmx = (1.0/3.0) * bc[0] * A;
				double d_v_d_nrmy = (1.0/3.0) * bc[1] * A;
				double d_v_d_nrmz = (1.0/3.0) * bc[2] * A;
				
				double d_v_d_A = (1.0/3.0) * (bc[0] * nrm[0] + bc[1] * nrm[1] + bc[2] * nrm[2]);

#ifdef DEBUG_A
				dv_dr[ir1*3+0] +=  0;
				dv_dr[ir1*3+1] +=  0;
				dv_dr[ir1*3+2] +=  0;
				
				dv_dr[ir2*3+0] +=  1;
				dv_dr[ir2*3+1] +=  0;
				dv_dr[ir2*3+2] +=  0;
				
				dv_dr[ir3*3+0] +=  0;
				dv_dr[ir3*3+1] +=  0;
				dv_dr[ir3*3+2] +=  0;
#else		
				dv_dr[ir1*3+0] += d_v_d_nrmx * d_nrmx_d_r1x;
				dv_dr[ir1*3+1] += d_v_d_nrmx * d_nrmx_d_r1y;
				dv_dr[ir1*3+2] += d_v_d_nrmx * d_nrmx_d_r1z;
				
				dv_dr[ir1*3+0] += d_v_d_nrmy * d_nrmy_d_r1x;
				dv_dr[ir1*3+1] += d_v_d_nrmy * d_nrmy_d_r1y;
				dv_dr[ir1*3+2] += d_v_d_nrmy * d_nrmy_d_r1z;
				
				dv_dr[ir1*3+0] += d_v_d_nrmz * d_nrmz_d_r1x;
				dv_dr[ir1*3+1] += d_v_d_nrmz * d_nrmz_d_r1y;
				dv_dr[ir1*3+2] += d_v_d_nrmz * d_nrmz_d_r1z;
				
				dv_dr[ir1*3+0] += d_v_d_A * d_A_d_r1x;
				dv_dr[ir1*3+1] += d_v_d_A * d_A_d_r1y;
				dv_dr[ir1*3+2] += d_v_d_A * d_A_d_r1z;
				
				dv_dr[ir2*3+0] += d_v_d_nrmx * d_nrmx_d_r2x;
				dv_dr[ir2*3+1] += d_v_d_nrmx * d_nrmx_d_r2y;
				dv_dr[ir2*3+2] += d_v_d_nrmx * d_nrmx_d_r2z;
				
				dv_dr[ir2*3+0] += d_v_d_nrmy * d_nrmy_d_r2x;
				dv_dr[ir2*3+1] += d_v_d_nrmy * d_nrmy_d_r2y;
				dv_dr[ir2*3+2] += d_v_d_nrmy * d_nrmy_d_r2z;
				
				dv_dr[ir2*3+0] += d_v_d_nrmz * d_nrmz_d_r2x;
				dv_dr[ir2*3+1] += d_v_d_nrmz * d_nrmz_d_r2y;
				dv_dr[ir2*3+2] += d_v_d_nrmz * d_nrmz_d_r2z;
				
				dv_dr[ir2*3+0] += d_v_d_A * d_A_d_r2x;
				dv_dr[ir2*3+1] += d_v_d_A * d_A_d_r2y;
				dv_dr[ir2*3+2] += d_v_d_A * d_A_d_r2z;
				
				dv_dr[ir3*3+0] += d_v_d_nrmx * d_nrmx_d_r3x;
				dv_dr[ir3*3+1] += d_v_d_nrmx * d_nrmx_d_r3y;
				dv_dr[ir3*3+2] += d_v_d_nrmx * d_nrmx_d_r3z;
				
				dv_dr[ir3*3+0] += d_v_d_nrmy * d_nrmy_d_r3x;
				dv_dr[ir3*3+1] += d_v_d_nrmy * d_nrmy_d_r3y;
				dv_dr[ir3*3+2] += d_v_d_nrmy * d_nrmy_d_r3z;
				
				dv_dr[ir3*3+0] += d_v_d_nrmz * d_nrmz_d_r3x;
				dv_dr[ir3*3+1] += d_v_d_nrmz * d_nrmz_d_r3y;
				dv_dr[ir3*3+2] += d_v_d_nrmz * d_nrmz_d_r3z;
				
				dv_dr[ir3*3+0] += d_v_d_A * d_A_d_r3x;
				dv_dr[ir3*3+1] += d_v_d_A * d_A_d_r3y;
				dv_dr[ir3*3+2] += d_v_d_A * d_A_d_r3z;
				
				// the barycenter:
				dv_dr[ir1*3+0] += (1.0/3.0) * nrm[0] * A * (1.0/3.0);
				dv_dr[ir1*3+1] += (1.0/3.0) * nrm[1] * A * (1.0/3.0);
				dv_dr[ir1*3+2] += (1.0/3.0) * nrm[2] * A * (1.0/3.0);
				
				dv_dr[ir2*3+0] += (1.0/3.0) * nrm[0] * A * (1.0/3.0);
				dv_dr[ir2*3+1] += (1.0/3.0) * nrm[1] * A * (1.0/3.0);
				dv_dr[ir2*3+2] += (1.0/3.0) * nrm[2] * A * (1.0/3.0);
				
				dv_dr[ir3*3+0] += (1.0/3.0) * nrm[0] * A * (1.0/3.0);
				dv_dr[ir3*3+1] += (1.0/3.0) * nrm[1] * A * (1.0/3.0);
				dv_dr[ir3*3+2] += (1.0/3.0) * nrm[2] * A * (1.0/3.0);
#endif
				if( my_isnan(d_v_d_A) || my_isnan(d_nrmx_d_r1x) || my_isnan(d_nrmx_d_r2x) || my_isnan(d_nrmx_d_r3x) )
				{
					printf("nan");
					exit(1);
				} 
			}

			
		}	 
		
		ioff=0;
		for( int iu = 0; iu < ni; iu++ )
		for( int iv = 0; iv < ni-iu; iv++, ioff++ )
		{
			double R[3] = {0,0,0};

			int *cp = theVolumeFormulas[f].cp;
			int np = theVolumeFormulas[f].ncoor;

			for( int p = 0; p < np; p++ )
			{
				g[cp[p]*3+0] += dv_dr[(iu*ni+iv)*3+0] * theVolumeFormulas[f].r_w[ioff*np+p] * scale;
				g[cp[p]*3+1] += dv_dr[(iu*ni+iv)*3+1] * theVolumeFormulas[f].r_w[ioff*np+p] * scale;
				g[cp[p]*3+2] += dv_dr[(iu*ni+iv)*3+2] * theVolumeFormulas[f].r_w[ioff*np+p] * scale;

			}
		}
	} 

	return vol;
	
}

double surface::volume( double *r)
{
	if( !theVolumeFormulas )
		generateVolumePlan();

	double vol = 0;	

	double rav = 0;
	double nav = 0;
	
	for( int f = 0; f < nfaces; f++ )
	{
		int ni = theVolumeFormulas[f].ni;

		int nv_space = ni*(ni+1)/2;

		double coords[3*ni*ni];

		int ioff=0;
		for( int iu = 0; iu < ni; iu++ )
		for( int iv = 0; iv < ni-iu; iv++, ioff++ )
		{
			double R[3] = {0,0,0};

			int *cp = theVolumeFormulas[f].cp;
			int np = theVolumeFormulas[f].ncoor;

			for( int p = 0; p < np; p++ )
			{
				R[0] += theVolumeFormulas[f].r_w[ioff*np+p] * (r[cp[p]*3+0] + theVolumeFormulas[f].r_pbc[ioff*np*3+3*p+0]); 
				R[1] += theVolumeFormulas[f].r_w[ioff*np+p] * (r[cp[p]*3+1] + theVolumeFormulas[f].r_pbc[ioff*np*3+3*p+1]); 
				R[2] += theVolumeFormulas[f].r_w[ioff*np+p] * (r[cp[p]*3+2] + theVolumeFormulas[f].r_pbc[ioff*np*3+3*p+2]); 
			}

			coords[(iu*ni+iv)*3+0] = R[0];
			coords[(iu*ni+iv)*3+1] = R[1];
			coords[(iu*ni+iv)*3+2] = R[2];
		}

		for( int i = 0; i < ni; i++ )
		for( int j = 0; j < ni-i; j++ )
		{
			
			for( int pass = 0; pass < 2; pass++ )
			{
				double *r1;
				double *r2;
				double *r3;

				if( pass == 0 && (i + j + 1 < ni && i < ni-1 && j < ni-1) )
				{
					r1 = coords + (i*ni+j)*3;
					r3 = coords + (i*ni+j+1)*3;
					r2 = coords + ((i+1)*ni+j)*3;
				}
				else if( pass == 1 && (i-1 >= 0 && j+1 < ni && i+j+1 < ni) )
				{
					r1 = coords + (i*ni+j)*3;
					r3 = coords + ((i-1)*ni+j+1)*3;
					r2 = coords + (i*ni+j+1)*3;
				}
				else
					continue;

				double nrm[3];

				double dr1[3] = { r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]};
				double dr2[3] = { r3[0]-r1[0],r3[1]-r1[1],r3[2]-r1[2]};

				cross( dr1, dr2, nrm );

				double A = normalize(nrm)/2;
				double bc[3] = { (r1[0]+r2[0]+r3[0])/3,
						 (r1[1]+r2[1]+r3[1])/3,
						 (r1[2]+r2[2]+r3[2])/3 };

#ifdef DEBUG_A
				vol += r2[0]; 
#else
				vol += (1.0/3.0) * (bc[0] * nrm[0] + bc[1] * nrm[1] + bc[2] * nrm[2]) * A; 

#endif
				rav += sqrt(bc[0]*bc[0]+bc[1]*bc[1]+bc[2]*bc[2]);
				nav += 1;

				if( !(vol < 0 || vol > -1 ) )
				{
					printf("nan.\n");
					exit(1);
				}
			}
		}	 
	} 

	return vol;
	
}

void surface::loadStartingConditions( double thick )
{
	// minimum distance from a fixed/protein point, and its height.
	double decay = 10.0;

	for( int i = 0; i < nv; i++ )
	{
		if( theVertices[i].protein_pt )
			continue;

		double min_height = thick;
		double min_r = 1e10;

		for( int j = 0; j < nv; j++ )
		{
			if( i == j ) continue;
			if( !theVertices[j].protein_pt ) continue;

			double dr[3] = { 
				theVertices[i].r[0] - theVertices[j].r[0],
				theVertices[i].r[1] - theVertices[j].r[1],
				theVertices[i].r[2] - theVertices[j].r[2] };

			double drn[3];
			MinImage3D( dr, PBC_vec, drn );		
			double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2] );

			if( r < min_r )
			{
				min_r = r;
				min_height = theVertices[j].r[2];
			}
		}		

		theVertices[i].r[2] = thick + (min_height - thick) * exp(-min_r/decay);
	}
	
}

int surface::loadAndCopyLattice( const char *fileName, surface *copyFrom )
{
	return loadLattice( fileName, 0.0, copyFrom );
}

int surface::loadLattice( const char *fileName, double noise, surface *copyFrom )
{
	on_surface = 0;
	disable_PBC = 0;
	cumulative_area = NULL;
	max_valence = 15; // points used.
	opencl_init = 0;
#ifdef FFTW
	h_in = NULL;
	h_out = NULL;
#endif
	bcs = NULL;
	nbc = 0;
	nfaces = 0;
	nfaces6 = 0;
	nfacesI = 0;

	ptBoxes = NULL;
	prBoxes = NULL;
	irr_formulas = NULL;
	single_formulas = NULL;

	edgeFormulas = NULL;
	theFormulas = NULL;
	theVolumeFormulas = NULL;
	nf_faces = 0;
	nf_irr_faces = 0;
	nf_irr_pts = 0;
	nf_g_q_p = 0;

	int n_rand = 1000;
	long *vals = (long *)malloc( sizeof(long) * n_rand);
	ran_start(1);
	ran_array(vals, n_rand );

	c0 = 0;
	PBC_vec[0][0] = 1e10;
	PBC_vec[0][1] = 0;
	PBC_vec[0][2] = 0;
	PBC_vec[1][0] = 0;
	PBC_vec[1][1] = 1e10;
	PBC_vec[1][2] = 0;
	PBC_vec[2][0] = 0;
	PBC_vec[2][1] = 0;
	PBC_vec[2][2] = 1e10;

	PBC[0] = 0;
	PBC[1] = 0;
	PBC[2] = 0;

	nv = 0;
	theVertices = NULL;

	char buffer[4096];
	FILE *theFile = fopen(fileName,"r");

	if( !theFile )
	{
		printf("Couldn't open lattice file '%s'.\n", fileName );
		exit(1);
	}

	getLine(theFile, buffer );
	int vec3 = 0;
	if( !strncasecmp( buffer, "3D", 2 ) ) vec3 = 1; 
	getLine(theFile, buffer );
		sscanf( buffer, "%lf %lf %lf", PBC_vec[0]+0, PBC_vec[0]+1, PBC_vec[0]+2 );
	getLine(theFile, buffer );
		sscanf( buffer, "%lf %lf %lf", PBC_vec[1]+0, PBC_vec[1]+1, PBC_vec[1]+2 );
	if( vec3 )
	{
		getLine(theFile, buffer );
		sscanf( buffer, "%lf %lf %lf", PBC_vec[2]+0, PBC_vec[2]+1, PBC_vec[2]+2 );
	}
	else
	{
		PBC_vec[2][0] = 0;
		PBC_vec[2][1] = 0;
		PBC_vec[2][2] = 1e6;
	}
	
	int fp = ftell(theFile);
	getLine( theFile, buffer );
	fix_sense = 1;
	if( !strncasecmp( buffer, "sense", 5) )
	{
		sscanf( buffer, "sense %d", &fix_sense );
		printf("Fixing sense at %d.\n", fix_sense );
	}
	else
		fseek( theFile, fp, SEEK_SET );

	nv = 0;
	int nvs = 10;
	theVertices = (vertex *)malloc( sizeof(vertex) * nvs );


	int line = 0;
	int total_valence = 0;
	while( !feof(theFile) )
	{
		getLine( theFile, buffer );
		line++;
		if( feof(theFile) ) break;
		int id;
		double r[3];
		int val;
		if( !strncasecmp( buffer, "ntri", 4 ) )
			break;
	
		int ns = sscanf( buffer, "P%d ", &id );

		double tan[3] = { 0,0,0}, der2[3] = {0,0,0}, midp[3] = {0,0,0}, midp_tan[3] = {0,0,0}, midp_der2[3] = { 0,0,0};

		int protein_pt = 0;
		int nfields = 0;
		int ext_prot = 0;
		if( ns == 1 )
		{
			nfields = sscanf( buffer, "P%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d", &id, r+0, r+1, r+2, 
					tan+0,tan+1,tan+2,			
					der2+0,der2+1,der2+2,
					midp+0,midp+1,midp+2,
					midp_tan+0,midp_tan+1,midp_tan+2,
					midp_der2+0,midp_der2+1,midp_der2+2,
					&val );
			protein_pt = 1;
			ext_prot = 1;
		}
		else
		{
			ns = sscanf( buffer, "Q%d ", &id );

			if( ns == 1 )
			{
				nfields = sscanf( buffer, "Q%d %lf %lf %lf %d", &id, r+0, r+1, r+2, &val );
				protein_pt = 1;
			}
			else
				nfields = sscanf( buffer, "%d %lf %lf %lf %d", &id, r+0, r+1, r+2, &val );
		}
		if( nfields >= 5 )
		{
			if( val < 5 || val > 7 )
			{
				printf("The valence of any point currently must be between 5 and 7.\n");
				exit(1);
			}
			if( val > MAX_VALENCE )
			{
				printf("Read point with valence %d.  MAX_VALENCE set to %d.\n",
					val, MAX_VALENCE );
			}
			if( nvs == nv )
			{
				nvs *= 2;
				theVertices = (vertex *)realloc( theVertices, nvs * sizeof(vertex) );
			}

			theVertices[nv].protein_pt = protein_pt;
			theVertices[nv].id = id;
			theVertices[nv].c0 = 0;
			theVertices[nv].valence = val;

			theVertices[nv].r[0] = r[0];
			theVertices[nv].r[1] = r[1];
			theVertices[nv].r[2] = r[2];

			total_valence += val;

			int p = goToField( buffer, nfields );

			if( p < 0 )
			{
				printf("Field read error at line %d.\n", line );
				return -1; 
			}

			readNInts( buffer+p, theVertices[nv].edges, MAX_VALENCE );	

			nv++;
		}
	}

	for( int i = 0; i < nv; i++ )
	{
		for( int ex = 0; ex < theVertices[i].valence; ex++ )
			theVertices[i].edgeIndices[ex] = -1;
	}


	int nload_tri = 0;
	int *load_tri = NULL;
	double *load_tilt = NULL;
	if( !strncasecmp( buffer, "ntri", 4 ) )
	{
		sscanf( buffer, "ntri %d", &nload_tri );
		
		load_tri = (int *)malloc( sizeof(int) * 3 * nload_tri );
		load_tilt = (double *)malloc( sizeof(double) * 3 * nload_tri );

		for( int t = 0; t < nload_tri; t++ )
		{
			getLine( theFile, buffer );				
			if( feof(theFile) )
			{
				printf("Error reading triangles.\n");
				exit(1);
			}
			int nr = sscanf(buffer, "%d %d %d %lf %lf %lf", load_tri+3*t+0, load_tri+3*t+1, load_tri+3*t+2, load_tilt+3*t+0, load_tilt+3*t+1, load_tilt+3*t+2 );
			if( nr < 6 )
			{
				load_tilt[3*t+0] = 0;
				load_tilt[3*t+1] = 0;
				load_tilt[3*t+2] = 1;
			}
		}
	}
	
	getLine( theFile, buffer );

	printf("Read %d lattice points.\n", nv);
	
	writeVertexXYZandPSFPeriodic( "check" );

	assignEdgePBC();

	nt = 0;
	nts = 10;

	theTriangles = (triangle *)malloc( sizeof(triangle) * nts );
	theEdges = (edge *)malloc( sizeof(edge) * (total_valence/2) ); 
	nedges = 0;
	nedgesSpace = total_valence/2;
	
	double com[3] = {0,0,0};
	for( int i = 0; i < nv; i++ )
	{
		theVertices[i].nfaces = 0;
		com[0] += theVertices[i].r[0];
		com[1] += theVertices[i].r[1];
		com[2] += theVertices[i].r[2];
	}

	com[0] /= nv;
	com[1] /= nv;
	com[2] /= nv;
	
	if( load_tri )
	{
		for( int t = 0; t < nload_tri; t++ )
		{
			int i = load_tri[3*t+0]; 
			int j = load_tri[3*t+1]; 
			int k = load_tri[3*t+2]; 

			if( j < i )
			{
				int t = i;
				i = j;
				j = t;
			}
			if( k < j )
			{
				int t = k;
				k = j;
				j = t;
			}
			if( j < i )
			{
				int t = i;
				i = j;
				j = t;
			}
					if( nt == nts )
					{
						nts *= 2;
						theTriangles = (triangle *)realloc( theTriangles, nts * sizeof(triangle) );
					}
	
					theTriangles[nt].ids[0] = i;
					theTriangles[nt].ids[1] = j;
					theTriangles[nt].ids[2] = k;
					theTriangles[nt].ft_vector = 0;

					theTriangles[nt].fixed_tilt[0] = load_tilt[3*t+0];
					theTriangles[nt].fixed_tilt[1] = load_tilt[3*t+1];
					theTriangles[nt].fixed_tilt[2] = load_tilt[3*t+2];
					int npp = 0;
					if( theVertices[i].protein_pt ) npp++;
					if( theVertices[j].protein_pt ) npp++;
					if( theVertices[k].protein_pt ) npp++;
					if( npp >= 2 )
					{
						if( fabs(theTriangles[nt].fixed_tilt[0]) > 0 || fabs(theTriangles[nt].fixed_tilt[1]) > 0 )
							theTriangles[nt].ft_vector = 1;
					}
					theVertices[i].faces[theVertices[i].nfaces++] = nt;
					theVertices[j].faces[theVertices[j].nfaces++] = nt;
					theVertices[k].faces[theVertices[k].nfaces++] = nt;

	
					int gotit[3] = { 0,0,0};
					for( int e = 0; e < nedges; e++ )
					{
						if( theEdges[e].vertices[0] == i && theEdges[e].vertices[1] == j )
							gotit[0] = 1+e; 
						if( theEdges[e].vertices[0] == i && theEdges[e].vertices[1] == k )
							gotit[1] = 1+e; 
						if( theEdges[e].vertices[0] == j && theEdges[e].vertices[1] == k )
							gotit[2] = 1+e; 
					}

					if( !gotit[0]  )
					{
						theEdges[nedges].vertices[0] = i;
						theEdges[nedges].vertices[1] = j;
						theEdges[nedges].faces[0] = nt;
						theEdges[nedges].faces[1] = -1;
						theEdges[nedges].faces[2] = -1;
						theEdges[nedges].code[0] = 2;
						theTriangles[nt].edges[0] = nedges;
						nedges++;
					}
					else
					{
						theEdges[gotit[0]-1].faces[1] = nt;
						theEdges[gotit[0]-1].code[1] = 2;
						theTriangles[nt].edges[0] = gotit[0]-1;
					}
					if( !gotit[1]  )
					{
						theEdges[nedges].vertices[0] = i;
						theEdges[nedges].vertices[1] = k;
						theEdges[nedges].faces[0] = nt;
						theEdges[nedges].faces[1] = -1;
						theEdges[nedges].faces[2] = -1;
						theEdges[nedges].code[0] = 1;
						theTriangles[nt].edges[2] = nedges;
						nedges++;
					}
					else
					{
						theEdges[gotit[1]-1].faces[1] = nt;
						theEdges[gotit[1]-1].code[1] = 1;
						theTriangles[nt].edges[2] = gotit[1]-1;
					}
					if( !gotit[2]  )
					{
						theEdges[nedges].vertices[0] = j;
						theEdges[nedges].vertices[1] = k;
						theEdges[nedges].faces[0] = nt;
						theEdges[nedges].faces[1] = -1;
						theEdges[nedges].faces[2] = -1;
						theEdges[nedges].code[0] = 0;
						theTriangles[nt].edges[1] = nedges;
						nedges++;
					}
					else
					{
						theEdges[gotit[2]-1].faces[1] = nt;
						theEdges[gotit[2]-1].code[1] = 0;
						theTriangles[nt].edges[1] = gotit[2]-1;
					}
					nt++;
		}

		for(int e = 0; e < nedges; e++ )
		{
			if( theEdges[e].faces[0] == -1 || theEdges[e].faces[1] == -1 )
			{
				printf("edge %d with vertices %d %d has incomplete faces.\n",
					e, theEdges[e].vertices[0], theEdges[e].vertices[1] );
			}
		}
	}
	else
	{
	for( int i = 0; i < nv; i++ )
	{
		for( int vi = 0; vi < theVertices[i].valence; vi++ )
		{
			int j = theVertices[i].edges[vi];

			if( j < i ) continue;


			for( int vj = 0; vj < theVertices[j].valence; vj++ )
			{
				int k = theVertices[j].edges[vj];

				if( k == i ) continue;
			
				if( k < j ) continue;
			
				int gotit = 0;

				for( int vk = 0; vk < theVertices[k].valence; vk++ )
				{
					int l = theVertices[k].edges[vk];

					if( l == i ) 
						gotit = 1;
				}

				if( gotit )
				{
					if( nt == nts )
					{
						nts *= 2;
						theTriangles = (triangle *)realloc( theTriangles, nts * sizeof(triangle) );
					}
	
					theTriangles[nt].ids[0] = i;
					theTriangles[nt].ids[1] = j;
					theTriangles[nt].ids[2] = k;

					theVertices[i].faces[theVertices[i].nfaces++] = nt;
					theVertices[j].faces[theVertices[j].nfaces++] = nt;
					theVertices[k].faces[theVertices[k].nfaces++] = nt;
	
					int gotit[3] = { 0,0,0};
					for( int e = 0; e < nedges; e++ )
					{
						if( theEdges[e].vertices[0] == i && theEdges[e].vertices[1] == j )
							gotit[0] = 1+e; 
						if( theEdges[e].vertices[0] == i && theEdges[e].vertices[1] == k )
							gotit[1] = 1+e; 
						if( theEdges[e].vertices[0] == j && theEdges[e].vertices[1] == k )
							gotit[2] = 1+e; 
					}

					int codes[3] = { 2, 1, 0 };
					int edge_code[3] = {0,2,1};
					int vput[3][2] = { {i,j}, {i,k}, {j,k} };
					for( int parse = 0; parse < 3; parse++ )
					{
						if( !gotit[parse] )
						{
							if( nedges == nedgesSpace )
							{
								printf("Something is wrong with the mesh. The program is trying to create more edges than the valence would have indicated.\n");
								exit(1);
							}

							theEdges[nedges].sense = 0;
							theEdges[nedges].fix_sense = 0;
							theEdges[nedges].vertices[0] = vput[parse][0];
							theEdges[nedges].vertices[1] = vput[parse][1];


	
							theEdges[nedges].faces[0] = nt;
							theEdges[nedges].faces[1] = -1;
							theEdges[nedges].faces[2] = -1;
							theEdges[nedges].code[0] = codes[parse];
							theTriangles[nt].edges[edge_code[parse]] = nedges;
							nedges++;
						}
						else
						{
							if( theEdges[gotit[parse]-1].faces[1]  == -1 )
							{
								theEdges[gotit[parse]-1].faces[1] = nt;
								theEdges[gotit[parse]-1].code[1] = codes[parse];
							}
							else
							{
								theEdges[gotit[parse]-1].faces[2] = nt;
								theEdges[gotit[parse]-1].code[2] = codes[parse];
							}
							theTriangles[nt].edges[edge_code[parse]] = gotit[parse]-1;
						}
					}
#if 0
					if( !gotit[0]  )
					{
						theEdges[nedges].vertices[0] = i;
						theEdges[nedges].vertices[1] = j;
						theEdges[nedges].faces[0] = nt;
						theEdges[nedges].faces[1] = -1;
						theEdges[nedges].faces[2] = -1;
						theEdges[nedges].code[0] = 2;
						theTriangles[nt].edges[0] = nedges;
						nedges++;
					}
					else
					{
						if( theEdges[gotit[0]-1].faces[1]  == -1 )
						{
							theEdges[gotit[0]-1].faces[1] = nt;
							theEdges[gotit[0]-1].code[1] = 2;
						}
						else
						{
							theEdges[gotit[0]-1].faces[2] = nt;
							theEdges[gotit[0]-1].code[2] = 2;
						}
						theTriangles[nt].edges[0] = gotit[0]-1;
					}
					if( !gotit[1]  )
					{
						theEdges[nedges].vertices[0] = i;
						theEdges[nedges].vertices[1] = k;
						theEdges[nedges].faces[0] = nt;
						theEdges[nedges].faces[1] = -1;
						theEdges[nedges].faces[2] = -1;
						theEdges[nedges].code[0] = 1;
						theTriangles[nt].edges[2] = nedges;
						nedges++;
					}
					else
					{
						if( theEdges[gotit[1]-1].faces[1]  == -1 )
						{
							theEdges[gotit[1]-1].faces[1] = nt;
							theEdges[gotit[1]-1].code[1] = 1;
						}
						else
						{
							theEdges[gotit[1]-1].faces[2] = nt;
							theEdges[gotit[1]-1].code[2] = 1;
						}
						theTriangles[nt].edges[2] = gotit[1]-1;
					}
					if( !gotit[2]  )
					{

						theEdges[nedges].vertices[0] = j;
						theEdges[nedges].vertices[1] = k;
						theEdges[nedges].faces[0] = nt;
						theEdges[nedges].faces[1] = -1;
						theEdges[nedges].faces[2] = -1;
						theEdges[nedges].code[0] = 0;
						theTriangles[nt].edges[1] = nedges;
						nedges++;
					}
					else
					{
						if( theEdges[gotit[2]-1].faces[1]  == -1 )
						{
							theEdges[gotit[2]-1].faces[1] = nt;
							theEdges[gotit[2]-1].code[1] = 0;
						}
						else
						{
							theEdges[gotit[2]-1].faces[2] = nt;
							theEdges[gotit[2]-1].code[2] = 0;
						}
						theTriangles[nt].edges[1] = gotit[2]-1;
					}
#endif
					nt++;
				}
			}
		}
	}
	}
	
	assignEdgeIndices();
	
	for( int i = 0; i < nv; i++ )
	{
		for( int ex = 0; ex < theVertices[i].valence; ex++ )
		{
			int j = theVertices[i].edges[ex];

			if( theVertices[i].edgeIndices[ex] == -1 )
			{
				theEdges[nedges].vertices[0] = i;
				theEdges[nedges].vertices[1] = j;

				theEdges[nedges].faces[0] = -1;
				theEdges[nedges].faces[1] = -1;
				theEdges[nedges].faces[2] = -1;

				for( int exj = 0; exj < theVertices[j].valence; exj++ )
				{
					if( theVertices[j].edges[exj] == i )
						theVertices[j].edgeIndices[exj] = nedges;
				}
				theVertices[i].edgeIndices[ex] = nedges;
				nedges++;
			}
		}
	}
	
	for( int t = 0; t < nt; t++ )
	{
		int i = theTriangles[t].ids[0];
		int j = theTriangles[t].ids[1];
		int k = theTriangles[t].ids[2];

		double r1[3] = { theVertices[i].r[0], theVertices[i].r[1], theVertices[i].r[2] };
		double r2[3] = { theVertices[j].r[0], theVertices[j].r[1], theVertices[j].r[2] };
		double r3[3] = { theVertices[k].r[0], theVertices[k].r[1], theVertices[k].r[2] };

		double dr1[3] = { r2[0] - r1[0], r2[1]-r1[1], r2[2] - r1[2] };
		double dr2[3] = { r3[0] - r1[0], r3[1]-r1[1], r3[2] - r1[2] };

		double put[3];

		MinImage3D( dr1, PBC_vec, theTriangles[t].pbc1 );
		MinImage3D( dr2, PBC_vec, theTriangles[t].pbc2 );
		double o[3] = {0,0,0};

//		dr1[2] = 0;
//		dr2[2] = 0;

		theTriangles[t].A0 = triangle_area( o, dr1, dr2 );
//		printf("triangle %d A0: %lf\n", t, theTriangles[t].A0  );
	}

	/*
		assign a sense to each edge as to how we are dividing up the tetrahedra.

	*/

	int alto = 0;
	for( int e = 0; e < nedges; e++ )
	{
		theEdges[e].sense = alto % 2;
		alto++;

		int p1 = theEdges[e].vertices[0]; 
		int p2 = theEdges[e].vertices[1]; 

		theEdges[e].fix_sense = 0;
		
	}


	int done = 0;
	int nswaps = 0;
	alto = 0;
	while( !done )
	{
		done = 1;

		for( int t = 0; t < nt; t++ )
		{
			// triangle has edge i to j, i to k, j to k.

			int sense[3] = {
				theEdges[theTriangles[t].edges[0]].sense,
				theEdges[theTriangles[t].edges[1]].sense,
				!theEdges[theTriangles[t].edges[2]].sense };

			if(
				(sense[0] && sense[1] && sense[2] ) ||	
				(!sense[0] && !sense[1] && !sense[2] ) )
			{
				int sw = ran_arr_cycle() % 3;
				alto++;
//				printf("sw: %d\n", sw );
				int mod = theTriangles[t].edges[sw];

				int p1 = theEdges[mod].vertices[0]; 
				int p2 = theEdges[mod].vertices[1]; 

				if( !theEdges[mod].fix_sense )
				{	
					theEdges[theTriangles[t].edges[sw]].sense = !theEdges[theTriangles[t].edges[sw]].sense;	
					nswaps++;	
				}
				done = 0;
			}
		}
		
	}



	for( int t = 0; t < nt; t++ )
	{
		int i = theTriangles[t].ids[0];
		int j = theTriangles[t].ids[1];
		int k = theTriangles[t].ids[2];
		
		double dr1[3] = { 
				theVertices[i].r[0] - theVertices[j].r[0],
				theVertices[i].r[1] - theVertices[j].r[1],
				theVertices[i].r[2] - theVertices[j].r[2] };
		double dr2[3] = { 
				theVertices[k].r[0] - theVertices[j].r[0],
				theVertices[k].r[1] - theVertices[j].r[1],
				theVertices[k].r[2] - theVertices[j].r[2] };

		double put[3] = {0,0,0};
		MinImage3D( dr1, PBC_vec, put );
		put[0] = 0;
		put[1] = 0;
		put[2] = 0;
		MinImage3D( dr2, PBC_vec, put );

		double fcom[3] = { 
			(3*theVertices[j].r[0] + dr1[0] + dr2[0] )/3,
			(3*theVertices[j].r[1] + dr1[1] + dr2[1] )/3,
			(3*theVertices[j].r[2] + dr1[2] + dr2[2] )/3 };

		normalize(fcom);
		double cp[3];

		cross( dr1, dr2, theTriangles[t].nrm);
		normalize(theTriangles[t].nrm);
		if( !(theTriangles[t].nrm[0] < 0 || theTriangles[t].nrm[0] > -2) )
		{
			printf("yikes.\n");
			exit(1);
		}
		theTriangles[t].dp  = fcom[0] * theTriangles[t].nrm[0];	
		theTriangles[t].dp += fcom[1] * theTriangles[t].nrm[1];	
		theTriangles[t].dp += fcom[2] * theTriangles[t].nrm[2];	

		theTriangles[t].sense = 0;
	}

	

	constructTriangles();
	int orientable_done = 0;
	int niter = 0;
	while( !orientable_done )
	{
		constructTriangles();

		for( int t = 0; t < nt; t++ )
			theTriangles[t].sense = 0;

		orientable_done = 1;
		theTriangles[0].sense = fix_sense;
	
		// set the sense of each triangle that we can.

		done = 0;
		while( !done )
		{
			done = 1;
			for( int e = 0; e < nedges; e++ )
			{
				int f1 = theEdges[e].faces[0];
				int f2 = theEdges[e].faces[1];
		
				if( f1 < 0 || f2 < 0 ) continue;
	
				if( theTriangles[f1].sense && theTriangles[f2].sense )
				{
					int target_sense = theTriangles[f1].sense * compare_sense( theEdges[e].code[0], theEdges[e].code[1] );
					if( theTriangles[f2].sense != target_sense )
					{
						printf("Surface is not orientable.\n");
						exit(1);
					}
				}

				if( theTriangles[f1].sense && !theTriangles[f2].sense )
				{
					theTriangles[f2].sense = theTriangles[f1].sense * compare_sense( theEdges[e].code[0], theEdges[e].code[1] );
					done = 0;
				}
				if( !theTriangles[f1].sense && theTriangles[f2].sense )
				{
					theTriangles[f1].sense = theTriangles[f2].sense * compare_sense( theEdges[e].code[0], theEdges[e].code[1] );
					done = 0;
				}
			}
		}
	

	
		int orientation_problem[nv];
		memset( orientation_problem, 0, sizeof(int) * nv );

		for( int i = 0; i < nv; i++ )
		{
			// we can tolerate two pockets of oriented points before there is a problem.

			int val = theVertices[i].valence;
			
			int *sorter = (int *)malloc( sizeof(int) * val );
			for( int p = 0; p < val; p++ )
				sorter[p] = -1;
			int nsorter = 0;
	
			for( int pocket = 0; pocket < 2; pocket++ )
			{
				int *local_sorter = (int *)malloc( sizeof(int) * val );
				for( int p = 0; p < val; p++ )
					local_sorter[p] = -1;
				int cur_spot = 1;
				int nrun = 0;

				for( int direction = 1; direction >= -1; direction -= 2 )
				{
					int to_place = 1;

					if( direction == 1)
					{		
						// find a place to start that hasn't been taken yet.

						int got_spot = 0;

						for( int j = 0; j < val; j++ )
						{
							int use = 1;

							for( int i = 0; i < nsorter; i++ )
							{
								if( sorter[i] == j )
									use = 0;	 
							}

							if( use )
							{
								got_spot = 1;
								local_sorter[0] = j;
								nrun = 1;
								break;
							}
						}

						if( !got_spot )
							break;

					}
					else
					{
						to_place = val-1;

						if( local_sorter[to_place] != -1 )
							break;
					}
		
					int last_vert = theVertices[i].edges[local_sorter[0]];
			
					int orientation_failure = 0;
		
					while( local_sorter[to_place] == -1 )
					{
						int x_cur_spot = to_place;
		
						for( int v = 0; v < theVertices[i].nfaces; v++ )
						{
							int tri = theVertices[i].faces[v];
	
							if( theTriangles[tri].sense == 0 )
								continue;
			
							int points[3] = { theTriangles[tri].ids[0], theTriangles[tri].ids[1], theTriangles[tri].ids[2] };
							
							if( points[0] == last_vert || points[1] == last_vert || points[2] == last_vert )
							{
								int next_pt = points[0];
								if( next_pt == i || next_pt == last_vert )
									next_pt = points[1];
								if( next_pt == i || next_pt == last_vert )
									next_pt = points[2];
				
								int type = -1;
								
								if( points[0] == i && points[1] == last_vert )
									type = 1;		
								if( points[1] == i && points[2] == last_vert )
									type = 1;		
								if( points[2] == i && points[0] == last_vert )
									type = 1;	
			
								if( type * theTriangles[tri].sense * direction < 0 )
								{
									for( int xv = 0; xv < val; xv++ )
									{
										if( theVertices[i].edges[xv] == next_pt )
											local_sorter[to_place] = xv;
									}
									last_vert = next_pt;

									if( direction > 0 )
										to_place++;
									else
										to_place--;
							
									if( to_place < 0 ) to_place += val;
									if( to_place >= val ) to_place -= val;
									nrun++;
									break;
								}	
							}

						}
						
						if( to_place == x_cur_spot )
							break;
					}
				}
				// we have a "run" of points.

				if( nrun + nsorter > val )
				{
					writeXYZandPSFPeriodic("bad");
//					writeXYZSurface("bad.xyz","bad.psf", this);
					printf("BAD ERROR.\n");
					exit(1);
				}

				int endp = val-1;
				if( nrun == val || local_sorter[endp] == -1 ) // take the whole thing
				{
					if( nrun == val && nsorter != 0 )
					{
						printf("here.\n");
					}
					memcpy( sorter+nsorter, local_sorter, sizeof(int) * nrun );
					nsorter += nrun;	
				}
				else // take a subset, account for loop-back.
				{
					int nforward = 0;
					int nrev     = 0;
					while( local_sorter[endp] != -1 && endp >= 0)
					{	
						nrev++;
						endp--;
					}		
					int firstp = 0;	
					while( local_sorter[firstp] != -1 && firstp < val)
					{	
						nforward++;
						firstp++;
					}	
					if( nrev > 0 )
					{
						memcpy( sorter+nsorter, local_sorter+val-nrev, sizeof(int)*nrev );
						nsorter+=nrev;
					}
					if( nforward > 0 )
					{
						memcpy( sorter+nsorter, local_sorter, sizeof(int)*nforward );
						nsorter+=nforward;
					}
				}

				free(local_sorter);

			}
			
			
			if( nsorter == val )
			{
				int oe[MAX_VALENCE];
				memcpy( oe, theVertices[i].edges, sizeof(int) * theVertices[i].valence );
				double opbc[3*MAX_VALENCE];
				memcpy( opbc, theVertices[i].edge_PBC, sizeof(double) * theVertices[i].valence*3 );
		
				for( int x = 0; x < theVertices[i].valence; x++ )
				{
					theVertices[i].edges[x] = oe[sorter[x]];
					theVertices[i].edge_PBC[3*x+0] = opbc[sorter[x]*3+0];
					theVertices[i].edge_PBC[3*x+1] = opbc[sorter[x]*3+1];
					theVertices[i].edge_PBC[3*x+2] = opbc[sorter[x]*3+2];
				}
			}
			else
			{
				orientation_problem[i] = 1;
				printf("Index %d could not be oriented!\n", i );
			}


			free(sorter);
		}
		
	
		setEdgeRev();

		assignEdgeIndices();
		// analyze for broken edges.

				
		int is_broken = 0;

		for( int i = 0; i < nv; i++ )
		{
			int val = theVertices[i].valence;

			for( int e = 0; e < val; e++ )
			{
				theVertices[i].broken[e] = 0;
				theVertices[i].broken_chain[e] = 0;
				theVertices[i].broken_looped[e] = 0;
		
				int ep1 = e+1;
				if( ep1 >= val ) ep1 -= val;
				int j = theVertices[i].edges[e];
				int k = theVertices[i].edges[ep1];
				int ej = theVertices[i].edge_rev[e];
				int ek = theVertices[i].edge_rev[ep1];
				int ekp1 = ek+1;
				if( ekp1 >= theVertices[k].valence )
					ekp1 -= theVertices[k].valence;
	
				if( theVertices[k].edges[ekp1] != j )
				{
					orientable_done = 0;

					if( orientation_problem[i] || orientation_problem[j] || orientation_problem[k] )
					{
				//		theVertices[i].broken_chain[e] = 1;
						theVertices[i].broken[e] = 1;
					}
					else
						theVertices[i].broken[e] = 1;
		
//					printf("broken: i: %d %d %d %d i: %d\n", i, j, k, theVertices[k].edges[ekp1], theVertices[k].edges[ek] );
					is_broken = 1;


				}
			}
		}	
		
		int max_loop = 2 * MAX_VALENCE;
		int *vert_loop = (int *)malloc( sizeof(int) * nv  );
		int *vert_loop_e = (int *)malloc( sizeof(int) * nv  );
		int *loops = (int *)malloc( sizeof(int) * nv * max_loop );
		int *loops_edge = (int *)malloc( sizeof(int) * nv * max_loop );
		int *loop_size = (int *)malloc( sizeof(int ) * nv );
		memset( loop_size, 0, sizeof(int) * nv );
		int *is_loop = (int *)malloc( sizeof(int) * nv );
		memset( is_loop, 0, sizeof(int) * nv );
		int nloops = 0;

		for( int i = 0; i < nv; i++ )
		{
			vert_loop[i] = -1;
			vert_loop_e[i] = -1;
		}

		for( int i = 0; i < nv; i++ )
		{
			if( orientation_problem[i] ) continue;

			int val = theVertices[i].valence;

			for( int e = 0; e < val; e++ )
			{
				if( theVertices[i].broken[e] && ! theVertices[i].broken_looped[e] && ! theVertices[i].broken_chain[e] )	
				{
					orientable_done = 0;
					int bad_loop = 0;

					loop_size[nloops] = 0;
					loops[nloops*max_loop+loop_size[nloops]] = i;
					loops_edge[nloops*max_loop+loop_size[nloops]] = e;
					loop_size[nloops] += 1;

					int cur_pt = i;
					int cur_edge = e;
					int next_pt = theVertices[i].edges[e];
					int term_with = -1;
					int term_with_e = -1;
					is_loop[nloops] = 1;
					while( next_pt != i )
					{
						int e_ind = theVertices[cur_pt].edgeIndices[cur_edge];

						int sense_defined = 0;

						for( int tx = 0; tx < 3; tx++ )
						{
							int tri1 = theEdges[e_ind].faces[tx];
							if( tri1 < 0 ) continue;
							if( theTriangles[tri1].sense )
								sense_defined = 1;
						}
				
						if( !sense_defined )
						{
							printf("loop");
							for( int x = 0; x < loop_size[nloops]; x++ )
								printf(" %d", loops[nloops*max_loop+x] );
							printf(" ending with undefined sense at index %d.\n", next_pt ); 
							term_with = next_pt;
							is_loop[nloops] = 0;
							bad_loop = 1;
							break;
						}


						int next_edge = theVertices[cur_pt].edge_rev[cur_edge]-1;
						if( next_edge < 0 )
							next_edge += theVertices[next_pt].valence;
						
						if( next_pt != i &&( !theVertices[next_pt].broken[next_edge] || theVertices[next_pt].broken_looped[next_edge]) )
						{
							if( theVertices[next_pt].broken_looped[next_edge] )
							{
								term_with = next_pt;
								term_with_e = next_edge;
							}
							is_loop[nloops] = 0;
							printf("ending with a problem.\n");
							break;
						} 
						
						loops[nloops*max_loop+loop_size[nloops]] = next_pt;
						loops_edge[nloops*max_loop+loop_size[nloops]] = next_edge;
						loop_size[nloops] += 1;

						if( orientation_problem[next_pt] )
						{
							is_loop[nloops] = 0;
							 break;
						}		
						theVertices[next_pt].broken_looped[next_edge] = 1;

						cur_pt = next_pt;
						cur_edge = next_edge;
						next_pt = theVertices[cur_pt].edges[cur_edge];
					}


			
					if( loop_size[nloops] >= 3 )
					{
						vert_loop[i] = nloops;	
						vert_loop_e[i] = e;

							printf("loop segment");
						for( int xx = 0; xx < loop_size[nloops]; xx++ )
							printf(" %d", loops[nloops*max_loop+xx] );
						printf("\n");

						nloops++;
					}
					else
					{
						printf("unused loop segment");
						for( int xx = 0; xx < loop_size[nloops]; xx++ )
						{
							int i = loops[nloops*max_loop+xx];
							int e = loops_edge[nloops*max_loop+xx];
							
							printf(" %d", loops[nloops*max_loop+xx] );

							theVertices[i].broken_looped[e] = 0;
						}
						printf("\n");
					}
				}
			}
		}
		
		assignEdgeIndices();

		int edge_created = 0;
		int edge_sealed = 0;
		for( int l = 0; l < nloops; l++ )
		{	
			if( loop_size[l] < 3 ) 
			{
				printf(" loop < 3, debug this.\n");
				continue;
//				exit(1);
			}

			validateCodes();

			// we are creating an edge between i and k, make sure ij and jk have triangles?

			int okay_s = 0;

			int valence[loop_size[l]];

			for( int x = 0; x < loop_size[l]; x++ )
				valence[x] = theVertices[loops[l*max_loop+x]].valence;

			for( int okay_s = 0; okay_s < loop_size[l]-2; okay_s++ )
			{
				int i = loops[l*max_loop+okay_s+0];
				int j = loops[l*max_loop+okay_s+1];
				int k = loops[l*max_loop+okay_s+2];

				if( loop_size[l] == 4 )
				{
					int vl = 0;
					if( okay_s == 0 )
						vl = loops[l*max_loop+3];
					else
						vl = loops[l*max_loop+0];

					if( theVertices[j].valence < 4 || theVertices[l].valence < 4 )
						continue;
					if( createEdge(i,j,k) )
					{
						edge_created = 1;
						break;
					}
				}	
				else if( loop_size[l] > 3 || !is_loop[l] )
				{
					if( createEdge( i, j, k ) )
					{
						edge_created = 1;
						break;
					}
				}
				else if ( is_loop[l] )
				{
					if( sealEdge( i, j, k ) )
					{
						edge_sealed = 1;
						break;
					}
	
				}
			}
			validateCodes();
		}	

		if( !edge_created && !orientable_done )
		{
			writeXYZandPSFPeriodic("bad");
			printf("Unable to create an edge.\n");
			exit(1);
		}	

		free(vert_loop);
		free(loops);
		free(loops_edge);
		free(loop_size);

		niter++;

		if( niter > 2000 )
		{
			writeXYZSurface("failed.xyz","failed.psf", this );
			
			exit(1);
		}
	}
			
	orientSurface();		

	
		
//	if( fix_sense == -1 )


	sortFaces();
//	writeXYZSurface("success.xyz","success.psf", this );
	
	for( int t = 0; t < nt; t++ )
	{
		int i = theTriangles[t].ids[0];
		int j = theTriangles[t].ids[1];
		int k = theTriangles[t].ids[2];

		double r1[3] = { theVertices[i].r[0], theVertices[i].r[1], theVertices[i].r[2] };
		double r2[3] = { theVertices[j].r[0], theVertices[j].r[1], theVertices[j].r[2] };
		double r3[3] = { theVertices[k].r[0], theVertices[k].r[1], theVertices[k].r[2] };

		double dr1[3] = { r2[0] - r1[0], r2[1]-r1[1], r2[2] - r1[2] };
		double dr2[3] = { r3[0] - r1[0], r3[1]-r1[1], r3[2] - r1[2] };

		double put[3];

		MinImage3D( dr1, PBC_vec, theTriangles[t].pbc1 );
		MinImage3D( dr2, PBC_vec, theTriangles[t].pbc2 );
		double o[3] = {0,0,0};


		theTriangles[t].A0 = triangle_area( o, dr1, dr2 );
	}

	
	int sense_zero = 0;
	for( int f = 0; f < nt; f++ )
	{
		if( theTriangles[f].sense == 0 )
		{
			printf("sense zero, indices %d %d %d\n", theTriangles[f].ids[0], theTriangles[f].ids[1], theTriangles[f].ids[2] );
			sense_zero = 1;
		}
	}
	if( sense_zero ) exit(1);

	for( int f = 0; f < nt; f++ )
	{
		theTriangles[f].nrm[0] *= theTriangles[f].sense;
		theTriangles[f].nrm[1] *= theTriangles[f].sense;
		theTriangles[f].nrm[2] *= theTriangles[f].sense;
	}

	for( int x = 0; x < nt; x++ )
	{
		if( theTriangles[x].sense == 0 )
		{
			printf("sense zero.\n");
			exit(1);
		}
		
		int t = x;

		int i = theTriangles[t].ids[0];
		int j = theTriangles[t].ids[1];
		int k = theTriangles[t].ids[2];
		
		double dr1[3] = { 
				theVertices[i].r[0] - theVertices[j].r[0],
				theVertices[i].r[1] - theVertices[j].r[1],
				theVertices[i].r[2] - theVertices[j].r[2] };
		double dr2[3] = { 
				theVertices[k].r[0] - theVertices[j].r[0],
				theVertices[k].r[1] - theVertices[j].r[1],
				theVertices[k].r[2] - theVertices[j].r[2] };

		double put[3];
		MinImage3D( dr1, PBC_vec, put );
		MinImage3D( dr2, PBC_vec, put );

		double fcom[3] = { 
			(3*theVertices[j].r[0] + dr1[0] + dr2[0] )/3,
			(3*theVertices[j].r[1] + dr1[1] + dr2[1] )/3,
			(3*theVertices[j].r[2] + dr1[2] + dr2[2] )/3 };

		normalize(fcom);
//		printf("%lf %lf %lf tri norm %lf %lf %lf\n", fcom[0], fcom[1], fcom[2], theTriangles[x].nrm[0], theTriangles[x].nrm[1], theTriangles[x].nrm[2] );
	}




	
	
	printf("Read %d triangles.\n", nt);

	for( int i = 0; i < nv; i++ )
	{
		for( int e = 0; e < theVertices[i].valence; e++ )
		{
			int j = theVertices[i].edges[e];

			for( int e2 = 0; e2 < theVertices[j].valence; e2++ )
			{
				if( theVertices[j].edges[e2] == i )
					theVertices[i].edge_rev[e] = e2;
			}
		}
	}

/*
	for( int f = 0; f < nt; f++ )
	{
		theTriangles[f].nrm[0] *= theTriangles[f].sense;
		theTriangles[f].nrm[1] *= theTriangles[f].sense;
		theTriangles[f].nrm[2] *= theTriangles[f].sense;
	}*/

	// compute vertex normals.

	for( int i = 0; i < nv; i++ )
	{
		{
			double atot = 0;
			double nrm[3] = { 0,0,0};
	
			for( int t = 0; t < theVertices[i].nfaces; t++ )
			{
				int f = theVertices[i].faces[t];
	
				atot += theTriangles[f].A0;
				
				nrm[0] += theTriangles[f].nrm[0] * theTriangles[f].A0; 
				nrm[1] += theTriangles[f].nrm[1] * theTriangles[f].A0; 
				nrm[2] += theTriangles[f].nrm[2] * theTriangles[f].A0; 
			}	
			

			theVertices[i].nrm[0] = nrm[0] / atot;
			theVertices[i].nrm[1] = nrm[1] / atot;
			theVertices[i].nrm[2] = nrm[2] / atot;
	
			
			

			normalize( theVertices[i].nrm );

			double dp = 
				theVertices[i].r[0] * theVertices[i].nrm[0] +
				theVertices[i].r[1] * theVertices[i].nrm[1] +
				theVertices[i].r[2] * theVertices[i].nrm[2];

//			printf("%lf %lf %lf nrm %lf %lf %lf dp: %lf\n", theVertices[i].r[0], theVertices[i].r[1], theVertices[i].r[2],
//				theVertices[i].nrm[0], theVertices[i].nrm[1], theVertices[i].nrm[2], dp ); 
		}

	}

	if( copyFrom )
	{
		for( int i = 0; i < nv; i++ )
		{
			memcpy( theVertices[i].broken, copyFrom->theVertices[i].broken, sizeof(char) * MAX_VALENCE );
			memcpy( theVertices[i].edges, copyFrom->theVertices[i].edges, sizeof(int) * MAX_VALENCE );
			memcpy( theVertices[i].edge_rev, copyFrom->theVertices[i].edge_rev, sizeof(int) * MAX_VALENCE );
			memcpy( theVertices[i].edge_PBC, copyFrom->theVertices[i].edge_PBC, sizeof(double) * MAX_VALENCE*3 );
			memcpy( theVertices[i].patch_A0, copyFrom->theVertices[i].patch_A0, sizeof(double) * MAX_VALENCE );
				
		}
	}
	else
	{

	// sort vertex faces and edges.

		for( int i = 0; i < nv; i++ )
		{
			int val = theVertices[i].valence;
	
			int sorter[val];
			for( int p = 0; p < val; p++ )
				sorter[p] = -1;
			int cur_spot = 1;
			sorter[0] = 0;
			
			int last_vert = theVertices[i].edges[0];
	
			while( cur_spot < val  )
			{
				for( int v = 0; v < val; v++ )
				{
					int tri = theVertices[i].faces[v];
	
					int points[3] = { theTriangles[tri].ids[0], theTriangles[tri].ids[1], theTriangles[tri].ids[2] };
					
					if( points[0] == last_vert || points[1] == last_vert || points[2] == last_vert )
					{
						int next_pt = points[0];
						if( next_pt == i || next_pt == last_vert )
							next_pt = points[1];
						if( next_pt == i || next_pt == last_vert )
							next_pt = points[2];
		
						int type = -1;
						
						if( points[0] == i && points[1] == last_vert )
							type = 1;		
						if( points[1] == i && points[2] == last_vert )
							type = 1;		
						if( points[2] == i && points[0] == last_vert )
							type = 1;	
	
						if( type * theTriangles[tri].sense < 0 )
						{
							for( int xv = 0; xv < val; xv++ )
							{
								if( theVertices[i].edges[xv] == next_pt )
									sorter[cur_spot] = xv;
							}
							last_vert = next_pt;
							cur_spot++;
							break;
						}	
					}
				}
			}
				
			int oe[MAX_VALENCE];
			memcpy( oe, theVertices[i].edges, sizeof(int) * theVertices[i].valence );
			double opbc[3*MAX_VALENCE];
			memcpy( opbc, theVertices[i].edge_PBC, sizeof(double) * theVertices[i].valence*3 );
			double otwist[3*MAX_VALENCE];
			double ob22[3*MAX_VALENCE];
			double oalt[3*MAX_VALENCE];
			double patch_area[3*MAX_VALENCE];
			double interpatch_area[3*MAX_VALENCE];
			char broken[MAX_VALENCE];
			memcpy( patch_area, theVertices[i].patch_A0, sizeof(double) *  theVertices[i].valence );	
			memcpy( broken, theVertices[i].broken, sizeof(char) * theVertices[i].valence );
			for( int x = 0; x < theVertices[i].valence; x++ )
			{
				theVertices[i].edges[x] = oe[sorter[x]];
				theVertices[i].edge_PBC[3*x+0] = opbc[sorter[x]*3+0];
				theVertices[i].edge_PBC[3*x+1] = opbc[sorter[x]*3+1];
				theVertices[i].edge_PBC[3*x+2] = opbc[sorter[x]*3+2];
				theVertices[i].patch_A0[x] = patch_area[sorter[x]];
				theVertices[i].broken[x] = broken[sorter[x]];
			}
		}

		if( 1 )
		{
			int t= 0;
		
			int i = theTriangles[t].ids[0];
			int j = theVertices[i].edges[0];
			int k = theVertices[i].edges[1];
		
			double r1[3] = { theVertices[i].r[0], theVertices[i].r[1], theVertices[i].r[2] };
			double r2[3] = { theVertices[j].r[0], theVertices[j].r[1], theVertices[j].r[2] };
			double r3[3] = { theVertices[k].r[0], theVertices[k].r[1], theVertices[k].r[2] };
		
			double dr1[3] = { r2[0] - r1[0], r2[1]-r1[1], r2[2] - r1[2] };
			double dr2[3] = { r3[0] - r1[0], r3[1]-r1[1], r3[2] - r1[2] };
		
			double put[3];
		
			MinImage3D( dr1, PBC_vec, put );
			MinImage3D( dr2, PBC_vec, put );
			double o[3] = {0,0,0};
		
			double cp[3];
			cross( dr2, dr1, cp );
		
			normalize(r1);
			normalize(cp);

			double dp = cp[0] * r1[0] + cp[1] * r1[1] + cp[2] * r1[2];

	
			if( dp > 0 && dp > 0.7)
			{
				printf("FLIPPING\n");
				flipOrientation();
				sortFaces();
			}
		}
	
		setEdgeRev();	

	}


	int ncr = 5;

	free(vals);
	
	int i = 0;		
	return 0;
}


void surface::subdivideSurface( surface *copy_surface )
{
/*
	This routine fills out the surface structure.

	It copies many of the variables from those from copy_surface. These are the initial lines of code.

*/

	on_surface = 0;
	disable_PBC = 0;
	cumulative_area = NULL;
	max_valence = 15; // points used.
	opencl_init = 0;
#ifdef FFTW
	h_in = NULL;
	h_out = NULL;
#endif
	bcs = NULL;
	nbc = 0;

	edgeFormulas = NULL;
	theFormulas = NULL;
	theVolumeFormulas = NULL;
	ptBoxes = NULL;
	prBoxes = NULL;
	nf_faces = 0;
	nf_g_q_p = 0;
	nf_irr_faces = 0;
	nf_irr_pts = 0;

	printf("Subdividing.\n");
	c0 = 0;
	PBC_vec[0][0] = copy_surface->PBC_vec[0][0];
	PBC_vec[0][1] =  copy_surface->PBC_vec[0][1];
	PBC_vec[0][2] =  copy_surface->PBC_vec[0][2];
	PBC_vec[1][0] =  copy_surface->PBC_vec[1][0];
	PBC_vec[1][1] =  copy_surface->PBC_vec[1][1];
	PBC_vec[1][2] =  copy_surface->PBC_vec[1][2];
	PBC_vec[2][0] =  copy_surface->PBC_vec[2][0];
	PBC_vec[2][1] =  copy_surface->PBC_vec[2][1];
	PBC_vec[2][2] =  copy_surface->PBC_vec[2][2];

	PBC[0] = copy_surface->PBC[0];
	PBC[1] = copy_surface->PBC[1];
	PBC[2] = copy_surface->PBC[2];

/*
	nvo is the Original Number of Vertices: the number from copy_surface.
*/

	int nvo = copy_surface->nv;

	nv = copy_surface->nv;
	
/*
	we calculate the new number of vertices so I can allocate arrays.
	I suppose I could calculate this by summing up all the valences instead.
*/

	for( int i = 0; i < nvo; i++ )
	{
		int val = copy_surface->theVertices[i].valence;

		for( int e = 0; e < val; e++ )
		{
			int j = copy_surface->theVertices[i].edges[e];
			int e2 = copy_surface->theVertices[i].edge_rev[e];

			if( j > i )
			{
				nv++;
			}
		}	
	}

	theVertices = (vertex *)malloc( sizeof(vertex) * nv );

	/*
		This code block computes the positions of the original indices in the next iteration of the subdivision.
		It uses the formula from Cirak.
	*/

	for( int v = 0; v < copy_surface->nv; v++ )
	{
		int val = copy_surface->theVertices[v].valence;

		double w = 3. / (8.*val);

		theVertices[v].r[0] = (1-val*w)*copy_surface->theVertices[v].r[0];
		theVertices[v].r[1] = (1-val*w)*copy_surface->theVertices[v].r[1];
		theVertices[v].r[2] = (1-val*w)*copy_surface->theVertices[v].r[2];

		for( int e = 0; e < val; e++ )
		{
			int j = copy_surface->theVertices[v].edges[e];

			double dr[3] = { 
				copy_surface->theVertices[j].r[0]- copy_surface->theVertices[v].r[0],
				copy_surface->theVertices[j].r[1]- copy_surface->theVertices[v].r[1],
				copy_surface->theVertices[j].r[2]- copy_surface->theVertices[v].r[2] };
			double add[3];
			MinImage3D(dr, PBC_vec, add );
			theVertices[v].r[0] += w*(copy_surface->theVertices[v].r[0]+dr[0]);
			theVertices[v].r[1] += w*(copy_surface->theVertices[v].r[1]+dr[1]);
			theVertices[v].r[2] += w*(copy_surface->theVertices[v].r[2]+dr[2]);
		}

		theVertices[v].valence = copy_surface->theVertices[v].valence;
	}

	// tnv is the index of the created indices. It starts at nvo, the index of the first new vertex.
	int tnv = nvo;	
	
	// this code block fills out some initial values for the new vertices.
	// *all* of the edge indexes for the original vertices change : they point to new subdivision points.
	// we fill these in as we create the vertices.

	for( int i = 0; i < nvo; i++ )
	{
		int val = copy_surface->theVertices[i].valence;

		for( int e = 0; e < val; e++ )
		{
			int em1 = e-1;
			int ep1 = e+1;

			// new variables are created for these values so that they can be wrapped:
			if( em1 < 0 ) em1 += val;
			if( ep1 >= val ) ep1 -= val;

			// These are computed using the convention that the edges are all stored with the same "clockwise" sense.
			// that is, if i and j are linked with edge e e2 (see j = ... below, and e2 = .. below)
			// then theVertices[i].edges[e-1] points to the same vertex as theVertices[j].edges[e2+1] because they form a triangle.
			int j = copy_surface->theVertices[i].edges[e];
			int k = copy_surface->theVertices[i].edges[em1];
			int l = copy_surface->theVertices[i].edges[ep1];
			int e2 = copy_surface->theVertices[i].edge_rev[e];
		
			if( j > i )
			{
				// the rule, from Cirak
				double fi = 3.0/8.0;
				double fj = 3.0/8.0;
				double fk = 1.0/8.0;
				double fl = 1.0/8.0;
			
				theVertices[tnv].r[0] = fi*copy_surface->theVertices[i].r[0];
                        	theVertices[tnv].r[1] = fi*copy_surface->theVertices[i].r[1];
                        	theVertices[tnv].r[2] = fi*copy_surface->theVertices[i].r[2];
				
				double dr[3] = { copy_surface->theVertices[j].r[0] - copy_surface->theVertices[i].r[0],
						 copy_surface->theVertices[j].r[1] - copy_surface->theVertices[i].r[1],
						 copy_surface->theVertices[j].r[2] - copy_surface->theVertices[i].r[2] };
				double add[3];
				MinImage3D(dr, PBC_vec, add );
				theVertices[tnv].r[0] += fj*(copy_surface->theVertices[i].r[0] +dr[0]);
                        	theVertices[tnv].r[1] += fj*(copy_surface->theVertices[i].r[1] +dr[1]);
                        	theVertices[tnv].r[2] += fj*(copy_surface->theVertices[i].r[2] +dr[2]);
				
				dr[0]= copy_surface->theVertices[k].r[0] - copy_surface->theVertices[i].r[0];
				dr[1]= copy_surface->theVertices[k].r[1] - copy_surface->theVertices[i].r[1];
				dr[2]= copy_surface->theVertices[k].r[2] - copy_surface->theVertices[i].r[2];
				MinImage3D(dr, PBC_vec, add );
				theVertices[tnv].r[0] += fk*(copy_surface->theVertices[i].r[0] +dr[0]);
                        	theVertices[tnv].r[1] += fk*(copy_surface->theVertices[i].r[1] +dr[1]);
                        	theVertices[tnv].r[2] += fk*(copy_surface->theVertices[i].r[2] +dr[2]);
				
				dr[0]= copy_surface->theVertices[l].r[0] - copy_surface->theVertices[i].r[0];
				dr[1]= copy_surface->theVertices[l].r[1] - copy_surface->theVertices[i].r[1];
				dr[2]= copy_surface->theVertices[l].r[2] - copy_surface->theVertices[i].r[2];
				MinImage3D(dr, PBC_vec, add );
				theVertices[tnv].r[0] += fl*(copy_surface->theVertices[i].r[0] +dr[0]);
                        	theVertices[tnv].r[1] += fl*(copy_surface->theVertices[i].r[1] +dr[1]);
                        	theVertices[tnv].r[2] += fl*(copy_surface->theVertices[i].r[2] +dr[2]);
				
				// this updates the edge index of the *original* vertices to point to the new indices.
				// the new indices have edges that are not yet set.
				theVertices[i].edges[e] = tnv;
				theVertices[j].edges[e2] = tnv;

				tnv++;
			}
		}	
	}

	tnv = nvo;
	
	// this code block sets the edges of the new indices.
	for( int i = 0; i < nvo; i++ )
	{
		int val = copy_surface->theVertices[i].valence;

		for( int e = 0; e < val; e++ )
		{
			int em1 = e-1;
			int ep1 = e+1;

			if( em1 < 0 ) em1 += val;
			if( ep1 >= val ) ep1 -= val;

			int j = copy_surface->theVertices[i].edges[e];
			int ej = copy_surface->theVertices[i].edge_rev[e];

			int ejm1 = ej-1;
			int ejp1 = ej+1;
			if( ejm1 < 0 ) ejm1 += copy_surface->theVertices[j].valence;
			if( ejp1 >= copy_surface->theVertices[j].valence ) ejp1 -= copy_surface->theVertices[j].valence;
			
			int k = copy_surface->theVertices[i].edges[em1];
			int l = copy_surface->theVertices[i].edges[ep1];

			if( j > i )
			{
				theVertices[tnv].valence = 6;
	
				// this is a counterclockwise loop around the point, i, i's point counterclockwise to it, j's point clockwise from j,
				// j, j's next counter-clockwise point, then the point clockwise from i.
				// (clockwise and counterclockwise can be interchanged here depending on how you're picturing it, as long as you're consistent).

				theVertices[tnv].edges[0] = i;
				theVertices[tnv].edges[1] = theVertices[i].edges[em1];
				theVertices[tnv].edges[2] = theVertices[j].edges[ejp1];
				theVertices[tnv].edges[3] = j;
				theVertices[tnv].edges[4] = theVertices[j].edges[ejm1];
				theVertices[tnv].edges[5] = theVertices[i].edges[ep1];
				tnv++;
			}
		}
		
	}
	
	theTriangles = NULL;
	theEdges = NULL;



	// here I fill out the indices that let us figure out edge indices of our partners	
	for( int i = 0; i < nv; i++ )
	{
		for( int e = 0; e < theVertices[i].valence; e++ )
		{
			int j = theVertices[i].edges[e];

			for( int e2 = 0; e2 < theVertices[j].valence; e2++ )
			{
				if( theVertices[j].edges[e2] == i )
					theVertices[i].edge_rev[e] = e2;
			}
		}
	}

#ifdef AMAT_DEBUG
	for( int i = 0; i < nv; i++ )
	{
		int val = theVertices[i].valence;

		if( theVertices[i].valence == 5)
		{
			int j = theVertices[i].edges[0];
			int k = theVertices[i].edges[1];
			printf("0: %.14le %.14le %.14le\n", 
				theVertices[i].r[0],
				theVertices[i].r[1],
				theVertices[i].r[2] );

			for( int v = 0; v < val; v++ )
			{
				int p = 1 - v;

				if( p < 0 ) p += val;

				int ind = theVertices[i].edges[p];

				printf("%d: %.14le %.14le %.14le\n", v+1,
					theVertices[ind].r[0],
					theVertices[ind].r[1],
					theVertices[ind].r[2] );
			} 

			int cp = 1+val;

			int ek = theVertices[i].edge_rev[1];
			int ej = theVertices[i].edge_rev[0];
			int ekm2 = ek-2; if( ekm2 < 0  ) ekm2 += theVertices[k].valence;
			int ekm3 = ek-3; if( ekm3 < 0  ) ekm3 += theVertices[k].valence;
			int ekm4 = ek-4; if( ekm4 < 0  ) ekm4 += theVertices[k].valence;

			int pr[3] = { ekm2, ekm3, ekm4 };
			for( int per = 0; per < 3; per++ )
			{
				int ind = theVertices[k].edges[pr[per]];

				printf("%d %.14le %.14le %.14le\n", cp, theVertices[ind].r[0], theVertices[ind].r[1], theVertices[ind].r[2] );
				cp++; 
			}
			
			int ejp2 = ej+2; if ( ejp2 >= theVertices[j].valence ) ejp2 -= theVertices[j].valence;
			int ejp3 = ej+3; if ( ejp3 >= theVertices[j].valence ) ejp3 -= theVertices[j].valence;
			pr[0] = ejp3;
			pr[1] = ejp2;
			for( int per = 0; per < 2; per++ )
			{
				int ind = theVertices[j].edges[pr[per]];

				printf("%d %.14le %.14le %.14le\n", cp, theVertices[ind].r[0], theVertices[ind].r[1], theVertices[ind].r[2] );
				cp++; 
			}
			break;
		}
	}
#endif
	for( int i = 0; i < nv; i++ )
		theVertices[i].c0 = 0;
	constructTriangles();
	assignEdgePBC();	
	sortFaces();	
}

void surface::assignEdgePBC( void )
{
	for( int i = 0; i < nv; i++ )
	{
		int val = theVertices[i].valence;

		for( int e = 0; e < val; e++ )
		{
			int k = theVertices[i].edges[e];

			double dr[3] = { theVertices[k].r[0] - theVertices[i].r[0],
					 theVertices[k].r[1] - theVertices[i].r[1],
					 theVertices[k].r[2] - theVertices[i].r[2] };
	
			MinImage3D( dr, PBC_vec, theVertices[i].edge_PBC+3*e );
		}
	}
}





void surface::put( double *r )
{
	for( int i = 0; i < nv; i++ )
	{
		theVertices[i].r[0]=r[3*i+0] ;
		theVertices[i].r[1]=r[3*i+1] ;
		theVertices[i].r[2]=r[3*i+2] ;
	}
}

void surface::get( double *r )
{
	for( int i = 0; i < nv; i++ )
	{
		r[3*i+0] = theVertices[i].r[0];
		r[3*i+1] = theVertices[i].r[1];
		r[3*i+2] = theVertices[i].r[2];
	}
}

void surface::setg0( double *r, double reset_factor )
{	
	if( !theFormulas )
		generatePlan();

	double avg = 0;
	double navg = 0;	
	double e = 0;
	double area = 0;
	double wgt = 0;
	for( int f = 0; f < nf_faces; f++ )
	{
		double A = 0;
		double c0 = 0;
		double c1 = 0;

		for( int p = 0; p < nf_g_q_p; p++ )
		{
			int frm = f*nf_g_q_p+p;
			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
			double nrm[3]={0,0,0}; 

			int *cp = theFormulas[f*nf_g_q_p+p].cp;
			int np = theFormulas[f*nf_g_q_p+p].ncoor;

			for( int p = 0; p < np; p++ )
			{
				R[0] += theFormulas[frm].r_w[p] * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				R[1] += theFormulas[frm].r_w[p] * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				R[2] += theFormulas[frm].r_w[p] * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				Ru[0] += theFormulas[frm].r_u[p] * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				Ru[1] += theFormulas[frm].r_u[p] * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				Ru[2] += theFormulas[frm].r_u[p] * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				Rv[0] += theFormulas[frm].r_v[p] * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				Rv[1] += theFormulas[frm].r_v[p] * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				Rv[2] += theFormulas[frm].r_v[p] * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuu[0] += theFormulas[frm].r_uu[p] * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuu[1] += theFormulas[frm].r_uu[p] * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuu[2] += theFormulas[frm].r_uu[p] * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuv[0] += theFormulas[frm].r_uv[p] * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuv[1] += theFormulas[frm].r_uv[p] * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuv[2] += theFormulas[frm].r_uv[p] * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSvv[0] += theFormulas[frm].r_vv[p] * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSvv[1] += theFormulas[frm].r_vv[p] * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSvv[2] += theFormulas[frm].r_vv[p] * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			}
		



		cross( Ru, Rv, nrm );
		normalize(nrm);

		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		double g = sqrt(RuRu*RvRv-RuRv*RuRv);
		avg += g;
		navg += 1;
		}

	} 

	avg /= navg;

	e = 0;
	area = 0;
	wgt = 0;

	double area_per_lipid = 65.;

	for( int f = 0; f < nf_faces; f++ )
	{
		double A = 0;
		double c0 = 0;
		double c1 = 0;
		double nlipids = 0;

		for( int p = 0; p < nf_g_q_p; p++ )
		{
			int frm = f*nf_g_q_p+p;
			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
			double nrm[3]={0,0,0}; 

			int *cp = theFormulas[f*nf_g_q_p+p].cp;
			int np = theFormulas[f*nf_g_q_p+p].ncoor;

			for( int p = 0; p < np; p++ )
			{
				R[0] += theFormulas[frm].r_w[p] * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				R[1] += theFormulas[frm].r_w[p] * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				R[2] += theFormulas[frm].r_w[p] * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				Ru[0] += theFormulas[frm].r_u[p] * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				Ru[1] += theFormulas[frm].r_u[p] * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				Ru[2] += theFormulas[frm].r_u[p] * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				Rv[0] += theFormulas[frm].r_v[p] * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				Rv[1] += theFormulas[frm].r_v[p] * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				Rv[2] += theFormulas[frm].r_v[p] * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuu[0] += theFormulas[frm].r_uu[p] * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuu[1] += theFormulas[frm].r_uu[p] * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuu[2] += theFormulas[frm].r_uu[p] * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuv[0] += theFormulas[frm].r_uv[p] * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuv[1] += theFormulas[frm].r_uv[p] * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuv[2] += theFormulas[frm].r_uv[p] * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSvv[0] += theFormulas[frm].r_vv[p] * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSvv[1] += theFormulas[frm].r_vv[p] * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSvv[2] += theFormulas[frm].r_vv[p] * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			}
		



		cross( Ru, Rv, nrm );
		normalize(nrm);

		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		double g = sqrt(RuRu*RvRv-RuRv*RuRv);

		if( reset_factor > 0 )
		{
			theFormulas[frm].g0 = g * reset_factor + (1-reset_factor) * theFormulas[frm].g0;
			theFormulas[frm].g0_base = theFormulas[frm].g0;
		}
		else
		{
			theFormulas[frm].g0 = g;
			theFormulas[frm].g0_base = g;
		}
		theFormulas[frm].RuRv0 = RuRv / sqrt(RuRu*RvRv);

		A += g * 0.5 * theFormulas[frm].weight;

		nlipids += g * 0.5 / area_per_lipid;
		}

		int t = theFormulas[f*nf_g_q_p].tri;

		theTriangles[t].nlipids = nlipids;
		theTriangles[t].f_lipids = 1.0;
		theTriangles[t].f_lipids_stashed = 1.0;
		

		if( theTriangles[t].composition.innerLeaflet &&
	            theTriangles[t].composition.outerLeaflet )
		{
			double i_area_unorm = 0;
			double o_area_unorm = 0;
	
			for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
			{
				i_area_unorm += theTriangles[t].composition.innerLeaflet[x] * bilayerComp.APL[x];
				o_area_unorm += theTriangles[t].composition.outerLeaflet[x] * bilayerComp.APL[x];
			}
	
			theTriangles[t].composition.A0 = A;
	
			double NFAC_i = A / i_area_unorm;
			double NFAC_o = A / i_area_unorm;
			
			for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
			{
				theTriangles[t].composition.innerLeaflet[x] *= NFAC_i;// / bilayerComp.APL[x];
				theTriangles[t].composition.outerLeaflet[x] *= NFAC_o;// / bilayerComp.APL[x];
			}
		}
	} 
	
	for( int f = 0; f < nf_irr_faces; f++ )
	{
		double A = 0;
		double c0 = 0;
		double c1 = 0;
		double nlipids = 0;

		for( int p = 0; p < nf_irr_pts; p++ )
		{
			int frm = f*nf_irr_pts+p;
			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
			double nrm[3]={0,0,0}; 

			int *cp = theIrregularFormulas[f*nf_irr_pts+p].cp;
			int np = theIrregularFormulas[f*nf_irr_pts+p].ncoor;

			for( int p = 0; p < np; p++ )
			{
				R[0] += theIrregularFormulas[frm].r_w[p] * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				R[1] += theIrregularFormulas[frm].r_w[p] * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				R[2] += theIrregularFormulas[frm].r_w[p] * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				Ru[0] += theIrregularFormulas[frm].r_u[p] * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				Ru[1] += theIrregularFormulas[frm].r_u[p] * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				Ru[2] += theIrregularFormulas[frm].r_u[p] * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				Rv[0] += theIrregularFormulas[frm].r_v[p] * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				Rv[1] += theIrregularFormulas[frm].r_v[p] * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				Rv[2] += theIrregularFormulas[frm].r_v[p] * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSuu[0] += theIrregularFormulas[frm].r_uu[p] * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSuu[1] += theIrregularFormulas[frm].r_uu[p] * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSuu[2] += theIrregularFormulas[frm].r_uu[p] * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSuv[0] += theIrregularFormulas[frm].r_uv[p] * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSuv[1] += theIrregularFormulas[frm].r_uv[p] * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSuv[2] += theIrregularFormulas[frm].r_uv[p] * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSvv[0] += theIrregularFormulas[frm].r_vv[p] * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSvv[1] += theIrregularFormulas[frm].r_vv[p] * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSvv[2] += theIrregularFormulas[frm].r_vv[p] * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
			}
		



		cross( Ru, Rv, nrm );
		normalize(nrm);

		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		double g = sqrt(RuRu*RvRv-RuRv*RuRv);
		theIrregularFormulas[frm].g0 = g;
		theIrregularFormulas[frm].g0_base = g;
		theIrregularFormulas[frm].RuRv0 = RuRv / sqrt(RuRu*RvRv);
		nlipids += g / area_per_lipid;

			A += g * theIrregularFormulas[frm].weight;
		}
		
		int t = theIrregularFormulas[f*nf_irr_pts].tri;

		theTriangles[t].nlipids = nlipids;
		theTriangles[t].f_lipids = 1.0;
		theTriangles[t].f_lipids_stashed = 1.0;

		if( theTriangles[t].composition.innerLeaflet &&
	            theTriangles[t].composition.outerLeaflet )
		{
			double i_area_unorm = 0;
			double o_area_unorm = 0;
	
			for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
			{
				i_area_unorm += theTriangles[t].composition.innerLeaflet[x] * bilayerComp.APL[x];
				o_area_unorm += theTriangles[t].composition.outerLeaflet[x] * bilayerComp.APL[x];
			}
			
			theTriangles[t].composition.A0 = A;
	
			double NFAC_i = A / i_area_unorm;
			double NFAC_o = A / i_area_unorm;
			
			for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
			{
				theTriangles[t].composition.innerLeaflet[x] *= NFAC_i;// / bilayerComp.APL[x];
				theTriangles[t].composition.outerLeaflet[x] *= NFAC_o;// / bilayerComp.APL[x];
			}
		}
	} 
}

double VA = 0,VC = 0, VA4=0, VR=0, AVC=0;
double av_h2 = 0, nav_h2 = 0;

/*
double surface::anglePenalty( double *r, int do_vertex )
{
	// theta = acos( drdu . drdv )
	// k(theta)^n ? 
}*/

double rho_scale = (1e-2);


double max_area_strain = 0;

double surface::faceEnergy( int f, double *r, double *p_uv, int doMonge )
{
	if( f < nf_faces )
	{
		if( doMonge )
			return fenergym( f, r, p_uv );
		else
			return fenergy( f, r, p_uv );
	}
	else
	{
		if( doMonge )
		{
			printf("Irregular Monge NYI.\n");
			exit(1);
		}
		else
			return ifenergy( f, r, p_uv );
	}
}

double surface::ifenergy( int f, double *r, double *p_uv )
{
	double alpha_x = r[3*nv];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];

	double e = 0;
		
	double face_energy_density = 0;
	double face_area = 0;

	for( int p = 0; p < nf_irr_pts; p++ )
	{
		double A0 = 0;
		double A = 0;
		double AP = 0;

		int frm = (f-nf_faces)*nf_irr_pts+p;
		double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
		double nrm[3]={0,0,0}; 

		int *cp = theIrregularFormulas[frm].cp;
		int np = theIrregularFormulas[frm].ncoor;
		int tri = theIrregularFormulas[frm].tri;

		for( int p = 0; p < np; p++ )
		{
			R[0] += theIrregularFormulas[frm].r_w[p] * alpha_x*(r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
			R[1] += theIrregularFormulas[frm].r_w[p] * alpha_y*(r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
			R[2] += theIrregularFormulas[frm].r_w[p] * alpha_z*(r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
			
			Ru[0] += theIrregularFormulas[frm].r_u[p] * alpha_x*(r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
			Ru[1] += theIrregularFormulas[frm].r_u[p] * alpha_y*(r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
			Ru[2] += theIrregularFormulas[frm].r_u[p] * alpha_z*(r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
			
			Rv[0] += theIrregularFormulas[frm].r_v[p] * alpha_x*(r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
			Rv[1] += theIrregularFormulas[frm].r_v[p] * alpha_y*(r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
			Rv[2] += theIrregularFormulas[frm].r_v[p] * alpha_z*(r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
			
			tSuu[0] += theIrregularFormulas[frm].r_uu[p] * alpha_x*(r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
			tSuu[1] += theIrregularFormulas[frm].r_uu[p] * alpha_y*(r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
			tSuu[2] += theIrregularFormulas[frm].r_uu[p] * alpha_z*(r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
			
			tSuv[0] += theIrregularFormulas[frm].r_uv[p] * alpha_x*(r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
			tSuv[1] += theIrregularFormulas[frm].r_uv[p] * alpha_y*(r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
			tSuv[2] += theIrregularFormulas[frm].r_uv[p] * alpha_z*(r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
			
			tSvv[0] += theIrregularFormulas[frm].r_vv[p] * alpha_x*(r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
			tSvv[1] += theIrregularFormulas[frm].r_vv[p] * alpha_y*(r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
			tSvv[2] += theIrregularFormulas[frm].r_vv[p] * alpha_z*(r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
		}
	
		cross( Ru, Rv, nrm );
		normalize(nrm);


		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		
		double RuPRuP = Ru[0] * Ru[0]/alpha_x/alpha_x + Ru[1] * Ru[1]/alpha_y/alpha_y + Ru[2]*Ru[2]/alpha_z/alpha_z;
		double RuPRvP = Ru[0] * Rv[0]/alpha_x/alpha_x + Ru[1] * Rv[1]/alpha_y/alpha_y + Ru[2]*Rv[2]/alpha_z/alpha_z;
		double RvPRvP = Rv[0] * Rv[0]/alpha_x/alpha_x + Rv[1] * Rv[1]/alpha_y/alpha_y + Rv[2]*Rv[2]/alpha_z/alpha_z;

		double g = sqrt(RuRu*RvRv-RuRv*RuRv);
//		if( f== 0 )//&& (cntr % 10000 == 0 ) )
//			printf("RuRv: %.14lf g: %le\n", RuRv/sqrt(RuRu)/sqrt(RvRv), g );
		double g_prior = sqrt(RuPRuP*RvPRvP-RuPRvP*RuPRvP);

		{
			double idet = 1.0 / (Ru[0] * Rv[1] - Ru[1] * Rv[0]); 
				
			double dudx = idet * Rv[1];
			double dudy =-idet * Rv[0];
			double dvdx =-idet * Ru[1];
			double dvdy = idet * Ru[0];

			double dhdx = Ru[2] * dudx + Rv[2] * dvdx;
			double dhdy = Ru[2] * dudy + Rv[2] * dvdy;
		
			av_h2 += (dhdx*dhdx+dhdy*dhdy)*g;
			nav_h2 += g;	
		}
		double nsuu = tSuu[0] * nrm[0] + tSuu[1] * nrm[1] + tSuu[2] * nrm[2];
		double nsuv = tSuv[0] * nrm[0] + tSuv[1] * nrm[1] + tSuv[2] * nrm[2];
		double nsvv = tSvv[0] * nrm[0] + tSvv[1] * nrm[1] + tSvv[2] * nrm[2];

		double Stot = (nsuu * RvRv + nsvv * RuRu -2*nsuv*RuRv)/(g*g);

//		Int[i*n+j] = 0.5 * kc * Stot * Stot * g; 
//		AInt[i*n+j] = g;

		double Sop[4] = { 1.0/(g*g) * (nsuu * RvRv  - nsuv * RuRv), (1.0/(g*g)) * ( nsuv*RvRv-nsvv*RuRv),
				  1.0/(g*g) * (nsuv * RuRu  - nsuu * RuRv), (1.0/(g*g)) * ( nsvv*RuRu-nsuv*RuRv) };

		double a = Sop[0];		
		double b = Sop[1];
		double c = Sop[2];
		double d = Sop[3];

		double c0 = theIrregularFormulas[frm].c0;
		double c1 = -0.5*(a+d-sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
		double c2 = -0.5*(a+d+sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
		
//		printf("%d %d IRREG %lf %lf %lf c %.14le %.14le g: %.14le w: %.14le\n", f, p, R[0], R[1], R[2], c1, c2, g, theIrregularFormulas[frm].weight );
			
		double en_tot = 0.5 * kc * (c1+c2-c0 ) * (c1+c2-c0) ;
		double en_tot_2 = 0.5 * kc * (c1+c2 ) * (c1+c2) ;
		double en_kg = kg * c1 * c2;

//		printf("e1: %lf e2: %lf\n e3: %lf\n", en_tot, en_tot_2, en_kg );

#ifdef LOCAL_LIPID_ENERGY
		double en = kg * c1 * c2;
#else
		double en = 0.5 * kc * (c1+c2-c0 ) * (c1+c2-c0) + kg * c1 * c2;
#endif
//		printf("e1: %lf e2: %lf\n e3: %lf\n en: %lf\n", en_tot, en_tot_2, en_kg, en );

		if( !( en >0 || en < 1 ) )
		{
			printf("nan en.\n");
			exit(1);
		}

//		printf("c1: %.14le c2: %.14le g: %.14le en: %.14le\n", c1, c2, g, en );
#ifdef FIXED_A
		double dA = theIrregularFormulas[frm].g0 * theIrregularFormulas[frm].weight;
#else
		double dA = g * theIrregularFormulas[frm].weight;
#endif	
		double dudv = 1.0; // this is merged into the weight, for irregular vertices.

		//Atot += dudv * dA;
		//Atot0 += dudv * theIrregularFormulas[frm].g0 * theIrregularFormulas[frm].weight;
		A += dudv * dA;
		A0 += dudv * theIrregularFormulas[frm].g0 * theIrregularFormulas[frm].weight;
		
		AP += dudv * g_prior * theIrregularFormulas[frm].weight;

		//area += dudv * dA;
	
		double RuRv0 = theIrregularFormulas[frm].RuRv0;
	
		double val = RuRv / sqrt(RuRu*RvRv);

		e += k_reg * (val-RuRv0) * (val - RuRv0);
		VR += k_reg * (val-RuRv0) * ( val - RuRv0);
//		VR += k_reg * (RuRv-RuRv0) * (RuRv-RuRv0);
//		if( g > 0 )
		{

#ifdef LOCAL_LIPID_ENERGY
			double atot_o=0, atot_i=0;
			for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
			{
				atot_o += bilayerComp.APL[x] * theTriangles[tri].composition.outerLeaflet[x];
				atot_i += bilayerComp.APL[x] * theTriangles[tri].composition.innerLeaflet[x];
			}
			for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
			{
				double f_o = bilayerComp.APL[x] * theTriangles[tri].composition.outerLeaflet[x] / atot_o;
				double f_i = bilayerComp.APL[x] * theTriangles[tri].composition.innerLeaflet[x] / atot_i;

				double dc_o = ( c1+c2 - bilayerComp.c0[x]);
				double dc_i = (-c1-c2 - bilayerComp.c0[x]);

				en += 0.5 * kc * dc_o*dc_o * f_o * 0.5;
				en += 0.5 * kc * dc_i*dc_i * f_i * 0.5;
			}
#endif
			e += dudv * dA * en; 			
			VC += dudv * dA * en;
			AVC += dudv *dA * (c1+c2);
//			face_energy_density += dudv * dA * 0.5 * kc * (c1+c2-c0 ) * (c1+c2-c0);
			face_energy_density += dudv * dA * en;
			face_area += dudv *dA ; 

			double a_strain = (A-A0)/A0;

			if( fabs(a_strain) > max_area_strain )
			{
				max_area_strain = fabs(a_strain);
			}
#ifdef MICRO_KA
			VA += 2 * 0.5 * A0 * micro_KA * ((AP-A0)/A0) * ((AP-A0)/A0);
			e += 2 * 0.5 * A0 * micro_KA * ((AP-A0)/A0) * ((AP-A0)/A0); 
#else
#ifndef GLOBAL_AREA
			VA += 2 * 0.5 * A0 * KA * ((A-A0)/A0) * ((A-A0)/A0);
			e += 2 * 0.5 * A0 * KA * ((A-A0)/A0) * ((A-A0)/A0); 
#endif
#endif
			 
		}
	}


	if( p_uv )
		e += penergy( f, r, p_uv, 0, face_energy_density/face_area );

	return e;
}

double surface::fenergy( int f, double *r, double *p_uv )
{
	double alpha_x = r[3*nv];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];

	double e = 0;
	
	double face_energy_density = 0;
	double face_area = 0;
	
	for( int p = 0; p < nf_g_q_p; p++ )
	{
		double A0 = 0;
		double A = 0;
		double AP = 0;

		if( debug_bit && (f != do_face || p != do_pt) )
			continue;
		int frm = f*nf_g_q_p+p;
		double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
		double nrm[3]={0,0,0}; 

		int *cp = theFormulas[f*nf_g_q_p+p].cp;
		int np = theFormulas[f*nf_g_q_p+p].ncoor;
		int tri = theFormulas[frm].tri;

		for( int p = 0; p < np; p++ )
		{

			R[0] += theFormulas[frm].r_w[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			R[1] += theFormulas[frm].r_w[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			R[2] += theFormulas[frm].r_w[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			Ru[0] += theFormulas[frm].r_u[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			Ru[1] += theFormulas[frm].r_u[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			Ru[2] += theFormulas[frm].r_u[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			Rv[0] += theFormulas[frm].r_v[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			Rv[1] += theFormulas[frm].r_v[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			Rv[2] += theFormulas[frm].r_v[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSuu[0] += theFormulas[frm].r_uu[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSuu[1] += theFormulas[frm].r_uu[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSuu[2] += theFormulas[frm].r_uu[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSuv[0] += theFormulas[frm].r_uv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSuv[1] += theFormulas[frm].r_uv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSuv[2] += theFormulas[frm].r_uv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSvv[0] += theFormulas[frm].r_vv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSvv[1] += theFormulas[frm].r_vv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSvv[2] += theFormulas[frm].r_vv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
		}
	

		cross( Ru, Rv, nrm );
		normalize(nrm);


		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		
			
		
		double RuPRuP = Ru[0] * Ru[0]/alpha_x/alpha_x + Ru[1] * Ru[1]/alpha_y/alpha_y + Ru[2]*Ru[2]/alpha_z/alpha_z;
		double RuPRvP = Ru[0] * Rv[0]/alpha_x/alpha_x + Ru[1] * Rv[1]/alpha_y/alpha_y + Ru[2]*Rv[2]/alpha_z/alpha_z;
		double RvPRvP = Rv[0] * Rv[0]/alpha_x/alpha_x + Rv[1] * Rv[1]/alpha_y/alpha_y + Rv[2]*Rv[2]/alpha_z/alpha_z;



		double g = sqrt(RuRu*RvRv-RuRv*RuRv);
//		if( f== 0 )//&& (cntr % 10000 == 0 ) )
//			printf("RuRv: %.14lf g: %le\n", RuRv/sqrt(RuRu)/sqrt(RvRv), g );
		double g_prior = sqrt(RuPRuP*RvPRvP-RuPRvP*RuPRvP);

		{
			double idet = 1.0 / (Ru[0] * Rv[1] - Ru[1] * Rv[0]); 
				
			double dudx = idet * Rv[1];
			double dudy =-idet * Rv[0];
			double dvdx =-idet * Ru[1];
			double dvdy = idet * Ru[0];

			double dhdx = Ru[2] * dudx + Rv[2] * dvdx;
			double dhdy = Ru[2] * dudy + Rv[2] * dvdy;
		
			av_h2 += (dhdx*dhdx+dhdy*dhdy)*g;
			nav_h2 += g;	
		}
		double nsuu = tSuu[0] * nrm[0] + tSuu[1] * nrm[1] + tSuu[2] * nrm[2];
		double nsuv = tSuv[0] * nrm[0] + tSuv[1] * nrm[1] + tSuv[2] * nrm[2];
		double nsvv = tSvv[0] * nrm[0] + tSvv[1] * nrm[1] + tSvv[2] * nrm[2];

		double Stot = (nsuu * RvRv + nsvv * RuRu -2*nsuv*RuRv)/(g*g);

//		Int[i*n+j] = 0.5 * kc * Stot * Stot * g; 
//		AInt[i*n+j] = g;

		double Sop[4] = { 1.0/(g*g) * (nsuu * RvRv  - nsuv * RuRv), (1.0/(g*g)) * ( nsuv*RvRv-nsvv*RuRv),
				  1.0/(g*g) * (nsuv * RuRu  - nsuu * RuRv), (1.0/(g*g)) * ( nsvv*RuRu-nsuv*RuRv) };

		double a = Sop[0];		
		double b = Sop[1];
		double c = Sop[2];
		double d = Sop[3];

		double c0 = theFormulas[frm].c0;
		double c1 = -0.5*(a+d-sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
		double c2 = -0.5*(a+d+sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
		
//		printf("%d %d REG %lf %lf %lf c %le %le g: %.14le\n", f, p, R[0], R[1], R[2], c1, c2, g  );
	
		double en_tot = 0.5 * kc * (c1+c2-c0 ) * (c1+c2-c0)  - 0.5 * kc * c0 * c0;
		double en_tot_2 = 0.5 * kc * (c1+c2 ) * (c1+c2) ;
		double en_kg = kg * c1 * c2;

//		printf("e1: %lf e2: %lf\n e3: %lf\n", en_tot, en_tot_2, en_kg );

#ifdef LOCAL_LIPID_ENERGY
		double en = kg * c1 * c2;
#else
		double en = 0.5 * kc * (c1+c2-c0 ) * (c1+c2-c0) + kg * c1 * c2;
#endif
//		printf("e1: %lf e2: %lf\n e3: %lf\n en: %lf\n", en_tot, en_tot_2, en_kg, en );

		if( !( en >0 || en < 1 ) )
		{
			printf("nan en.\n");
			exit(1);
		}

//		printf("c1: %.14le c2: %.14le g: %.14le en: %.14le\n", c1, c2, g, en );
#ifdef FIXED_A
		double dA = theFormulas[frm].g0 * theFormulas[frm].weight;
#else
		double dA = g * theFormulas[frm].weight;
#endif	
		double dudv = 0.5;


		//Atot += dudv * dA;
		//Atot0 += dudv * theFormulas[frm].g0 * theFormulas[frm].weight;
		A += dudv * dA;
		A0 += dudv * theFormulas[frm].g0 * theFormulas[frm].weight;
		
		AP += dudv * g_prior * theFormulas[frm].weight;

		//area += dudv * dA;

		double RuRv0 = theFormulas[frm].RuRv0;
		double val = RuRv / sqrt(RuRu*RvRv);
	
		//e += k_reg * (RuRv-RuRv0) * (RuRv-RuRv0);
		//VR += k_reg * (RuRv-RuRv0) * (RuRv-RuRv0);
		e += k_reg * (val-RuRv0) * (val - RuRv0);
		VR += k_reg * (val-RuRv0) * ( val - RuRv0);

//		if( g > 0 )
		{
			if( en > 1e8 )
			{
//				printf("Hm: e: %le %lf %lf\n", en, dudv, dA );
			}
#ifdef LOCAL_LIPID_ENERGY
			double atot_o=0, atot_i=0;
			for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
			{
				atot_o += bilayerComp.APL[x] * theTriangles[tri].composition.outerLeaflet[x];
				atot_i += bilayerComp.APL[x] * theTriangles[tri].composition.innerLeaflet[x];
			}
			for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
			{
				double f_o = bilayerComp.APL[x] * theTriangles[tri].composition.outerLeaflet[x] / atot_o;
				double f_i = bilayerComp.APL[x] * theTriangles[tri].composition.innerLeaflet[x] / atot_i;

				double dc_o = ( c1+c2 - bilayerComp.c0[x]);
				double dc_i = (-c1-c2 - bilayerComp.c0[x]);

				en += 0.5 * kc * dc_o*dc_o * f_o * 0.5;
				en += 0.5 * kc * dc_i*dc_i * f_i * 0.5;
			}
#endif
			e += dudv * dA * en; 			
			VC += dudv * dA * en;
//			if( f == 40 && p == 0 )
//			printf("e %d %d %le dA %le c1 %le c2 %le\n", f, p, dudv*dA*en,
//				dA, c1, c2 );
			AVC += dudv *dA * (c1+c2);
			face_energy_density += dudv * dA * 0.5 * kc * (c1+c2-c0 ) * (c1+c2-c0);
			face_area += dudv *dA; 

			double a_strain = (A-A0)/A0;

			if( fabs(a_strain) > max_area_strain )
			{
				max_area_strain = fabs(a_strain);
			}


#ifdef MICRO_KA
			VA += 2 * 0.5 * A0 * micro_KA * ((AP-A0)/A0) * ((AP-A0)/A0);
			e += 2 * 0.5 * A0 * micro_KA * ((AP-A0)/A0) * ((AP-A0)/A0); 
#else
#ifndef GLOBAL_AREA
			VA += 2 * 0.5 * A0 * KA * ((A-A0)/A0) * ((A-A0)/A0);
			e += 2 * 0.5 * A0 * KA * ((A-A0)/A0) * ((A-A0)/A0); 

			
#endif
#endif
			 
		}
	}


	if( p_uv )
		e += penergy( f, r, p_uv, 0, face_energy_density/face_area );

	return e;
}

double surface::fenergym( int f, double *r, double *p_uv )
{

  double alpha_x = r[3*nv];
  double alpha_y = r[3*nv+1];
  double alpha_z = r[3*nv+2];

  double e = 0;

  double face_energy_density = 0;
  double face_area = 0;

  for( int p = 0; p < nf_g_q_p; p++ )
    {
      double A0 = 0;
      double A = 0;
      double AP = 0;
      double Atot = 0;
      double Atot0 = 0;
      double area = 0;

      int frm = f*nf_g_q_p+p;
	double c0 = theFormulas[frm].c0;
      double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
      double nrm[3]={0,0,0};

      int *cp = theFormulas[f*nf_g_q_p+p].cp;
      int np = theFormulas[f*nf_g_q_p+p].ncoor;

      for( int p = 0; p < np; p++ )
	{
	  R[0] += theFormulas[frm].r_w[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]);
	  R[1] += theFormulas[frm].r_w[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]);
	  R[2] += theFormulas[frm].r_w[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]);

	  Ru[0] += theFormulas[frm].r_u[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]);
	  Ru[1] += theFormulas[frm].r_u[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]);
	  Ru[2] += theFormulas[frm].r_u[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]);

	  Rv[0] += theFormulas[frm].r_v[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]);
	  Rv[1] += theFormulas[frm].r_v[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]);
	  Rv[2] += theFormulas[frm].r_v[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]);

	  tSuu[0] += theFormulas[frm].r_uu[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]);
	  tSuu[1] += theFormulas[frm].r_uu[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]);
	  tSuu[2] += theFormulas[frm].r_uu[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]);

	  tSuv[0] += theFormulas[frm].r_uv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]);
	  tSuv[1] += theFormulas[frm].r_uv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]);
	  tSuv[2] += theFormulas[frm].r_uv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]);

	  tSvv[0] += theFormulas[frm].r_vv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]);
	  tSvv[1] += theFormulas[frm].r_vv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]);
	  tSvv[2] += theFormulas[frm].r_vv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]);
	}

      cross( Ru, Rv, nrm );
      normalize(nrm);


      double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
      double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
      double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

      double RuPRuP = Ru[0] * Ru[0]/alpha_x/alpha_x + Ru[1] * Ru[1]/alpha_y/alpha_y + Ru[2]*Ru[2]/alpha_z/alpha_z;
      double RuPRvP = Ru[0] * Rv[0]/alpha_x/alpha_x + Ru[1] * Rv[1]/alpha_y/alpha_y + Ru[2]*Rv[2]/alpha_z/alpha_z;
      double RvPRvP = Rv[0] * Rv[0]/alpha_x/alpha_x + Rv[1] * Rv[1]/alpha_y/alpha_y + Rv[2]*Rv[2]/alpha_z/alpha_z;

      double g = sqrt(RuRu*RvRv-RuRv*RuRv);

      double g_prior = sqrt(RuPRuP*RvPRvP-RuPRvP*RuPRvP);
      
      double idet = 1.0 / (Ru[0] * Rv[1] - Ru[1] * Rv[0]);

      double dudx = idet * Rv[1];
      double dudy =-idet * Rv[0];
      double dvdx =-idet * Ru[1];
      double dvdy = idet * Ru[0];

      double dhdx = Ru[2] * dudx + Rv[2] * dvdx;
      double dhdy = Ru[2] * dudy + Rv[2] * dvdy;
      
      double d2hdx2 = tSuu[2] * dudx * dudx + tSuv[2] * dudx * dvdx + tSuv[2] * dvdx * dudx + tSvv[2] * dvdx * dvdx;
      double d2hdy2 = tSuu[2] * dudy * dudy + tSuv[2] * dudy * dvdy + tSuv[2] * dvdy * dudy + tSvv[2] * dvdy * dvdy; 

	                	                
      av_h2 += (dhdx*dhdx+dhdy*dhdy)*g;
      nav_h2 += g;
      
      double en_tot = 0.5 * kc * (d2hdy2 + d2hdx2 - c0) * (d2hdy2 + d2hdx2 - c0);
      double en_tot_2 = 0.5 * kc * (d2hdy2 + d2hdx2) * (d2hdy2 + d2hdx2);

      double en = 0.5 * kc * (d2hdy2 + d2hdx2 - c0) * (d2hdy2 + d2hdx2 -c0) - 0.5 * kc * c0 * c0;
//      printf("etot: %lf etot2: %lf\n ekg: %lf\n en: %lf\n", en_tot, en_tot_2, en_kg, en );
      if( !( en >0 || en < 1 ) ) 
	{
	  printf("nan en.\n");
	  exit(1);
	}
      
#ifdef FIXED_A
      double dA = theFormulas[frm].g0 * theFormulas[frm].weight;
#else
      double dA = g * theFormulas[frm].weight;
#endif
      double dudv = 0.5;

      A += dudv * dA;
      Atot += dudv * dA;
      A0 += dudv * theFormulas[frm].g0 * theFormulas[frm].weight;
      Atot0 += dudv * theFormulas[frm].g0 * theFormulas[frm].weight;
      AP += dudv * g_prior * theFormulas[frm].weight;


      area += dudv * dA;

      {
	if( en > 1e8 )
	  {
//	    printf("Hm: e: %le %lf %lf\n", en, dudv, dA );
	  }
	e += dudv * dA * en;
	VC += dudv * dA * en;
	AVC += dudv *dA * (d2hdy2+d2hdx2);
	face_energy_density += dudv * dA * 0.5 * kc * (d2hdy2 + d2hdx2-c0 ) * (d2hdy2 + d2hdx2-c0);
	face_area += dudv *dA;  

	double a_strain = (A-A0)/A0;

	if( fabs(a_strain) > max_area_strain )
	  {
	    max_area_strain = fabs(a_strain);
	  }

#ifdef MICRO_KA
	VA += 2 * 0.5 * A0 * micro_KA * ((AP-A0)/A0) * ((AP-A0)/A0);
	e += 2 * 0.5 * A0 * micro_KA * ((AP-A0)/A0) * ((AP-A0)/A0);
#else
#ifndef GLOBAL_AREA
	VA += 2 * 0.5 * A0 * KA * ((A-A0)/A0) * ((A-A0)/A0);
	e += 2 * 0.5 * A0 * KA * ((A-A0)/A0) * ((A-A0)/A0);
#endif
#endif

      }

    }
 
	if( p_uv )
		e += penergy( f, r, p_uv, 1 /* do monge*/, face_energy_density/face_area );

  return e;
}

void surface::computeModifiedVertices( void )
{
	if( faces_for_vertex == NULL )
	{
		faces_for_vertex = (int**)malloc( sizeof(int*) * nv );
		nfaces_for_vertex = (int *)malloc( sizeof(int) * nv );
		nfacesSpace_for_vertex = (int *)malloc( sizeof(int) * nv );

		for( int v = 0; v < nv; v++ )
		{
			nfaces_for_vertex[v] = 0;
			nfacesSpace_for_vertex[v] = 10;
			faces_for_vertex[v] = (int *)malloc( sizeof(int) * nfacesSpace_for_vertex[v] );
		}

		for( int f = 0; f < nf_faces; f++ )
		{ // nf_g_q_p
			int frm = f * nf_g_q_p;
			for( int x = 0; x < theFormulas[frm].ncoor; x++ )
			{
				int i = theFormulas[frm].cp[x];

				if( nfacesSpace_for_vertex[i] == nfaces_for_vertex[i] )
				{
					nfacesSpace_for_vertex[i] *= 2;
					faces_for_vertex[i] = (int *)realloc( faces_for_vertex[i], sizeof(int) * nfacesSpace_for_vertex[i] );
				}

				faces_for_vertex[i][nfaces_for_vertex[i]] = f;
				nfaces_for_vertex[i]++;
			}	
		}
		
		for( int f = 0; f < nf_irr_faces; f++ )
		{ // nf_g_q_p
			int frm = f * nf_irr_pts;
			for( int x = 0; x < theIrregularFormulas[frm].ncoor; x++ )
			{
				int i = theIrregularFormulas[frm].cp[x];

				if( nfacesSpace_for_vertex[i] == nfaces_for_vertex[i] )
				{
					nfacesSpace_for_vertex[i] *= 2;
					faces_for_vertex[i] = (int *)realloc( faces_for_vertex[i], sizeof(int) * nfacesSpace_for_vertex[i] );
				}

				faces_for_vertex[i][nfaces_for_vertex[i]] = nf_faces+f;
				nfaces_for_vertex[i]++;
			}	
		}
	}
}

void surface::getModifiedFaces( int do_vertex, int *flist, int *nf, int *plist, int *np )
{
	if( !theFormulas )
		generatePlan();
	if( faces_for_vertex == NULL )
		computeModifiedVertices();	
/*
	if( faces_for_vertex == NULL )
	{
		faces_for_vertex = (int**)malloc( sizeof(int*) * nv );
		nfaces_for_vertex = (int *)malloc( sizeof(int) * nv );
		nfacesSpace_for_vertex = (int *)malloc( sizeof(int) * nv );

		for( int v = 0; v < nv; v++ )
		{
			nfaces_for_vertex[v] = 0;
			nfacesSpace_for_vertex[v] = 10;
			faces_for_vertex[v] = (int *)malloc( sizeof(int) * nfacesSpace_for_vertex[v] );
		}

		for( int f = 0; f < nf_faces; f++ )
		{ // nf_g_q_p
			int frm = f * nf_g_q_p;
			for( int x = 0; x < theFormulas[frm].ncoor; x++ )
			{
				int i = theFormulas[frm].cp[x];

				if( nfacesSpace_for_vertex[i] == nfaces_for_vertex[i] )
				{
					nfacesSpace_for_vertex[i] *= 2;
					faces_for_vertex[i] = (int *)realloc( faces_for_vertex[i], sizeof(int) * nfacesSpace_for_vertex[i] );
				}

				faces_for_vertex[i][nfaces_for_vertex[i]] = f;
				nfaces_for_vertex[i]++;
			}	
		}
		
		for( int f = 0; f < nf_irr_faces; f++ )
		{ // nf_g_q_p
			int frm = f * nf_irr_pts;
			for( int x = 0; x < theIrregularFormulas[frm].ncoor; x++ )
			{
				int i = theIrregularFormulas[frm].cp[x];

				if( nfacesSpace_for_vertex[i] == nfaces_for_vertex[i] )
				{
					nfacesSpace_for_vertex[i] *= 2;
					faces_for_vertex[i] = (int *)realloc( faces_for_vertex[i], sizeof(int) * nfacesSpace_for_vertex[i] );
				}

				faces_for_vertex[i][nfaces_for_vertex[i]] = nf_faces+f;
				nfaces_for_vertex[i]++;
			}	
		}
	}
*/
	
	// copy the modified faces into the requested array.
	memcpy( flist, faces_for_vertex[do_vertex], nfaces_for_vertex[do_vertex] * sizeof(int) );
	*nf = nfaces_for_vertex[do_vertex];
	*np = 0;

	for( int fx = 0; fx < *nf; fx++ )
	{
		int f = faces_for_vertex[do_vertex][fx];

		int t;

		if( f >= nf_faces )
			t = theIrregularFormulas[(f-nf_faces)*nf_irr_pts].tri;
		else
			t = theFormulas[f*nf_g_q_p].tri;

		for( int px = 0; px < theTriangles[t].np; px++ )
		{
			plist[*np] = theTriangles[t].plist[px];
			for( int i = 0; i < *np; i++ )
			{
				if( plist[i] == plist[*np] )
				{
					printf("Dupe somehow.\n");
				}
			}
			*np += 1;
		} 
	} 
}


double surface::energy( double *r, double *puv, int do_vertex, int *plist, int *np_found, int doMonge )
{
	if( !theFormulas )
		generatePlan();

	double alpha_x = r[3*nv];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];
	
	static int icntr = 0;
	
	VA = 0;
	VA4 = 0;
	VC = 0;
	AVC = 0;
	if( faces_for_vertex == NULL )
		computeModifiedVertices();	
/*	if( faces_for_vertex == NULL )
	{
		faces_for_vertex = (int**)malloc( sizeof(int*) * nv );
		nfaces_for_vertex = (int *)malloc( sizeof(int) * nv );
		nfacesSpace_for_vertex = (int *)malloc( sizeof(int) * nv );

		for( int v = 0; v < nv; v++ )
		{
			nfaces_for_vertex[v] = 0;
			nfacesSpace_for_vertex[v] = 10;
			faces_for_vertex[v] = (int *)malloc( sizeof(int) * nfacesSpace_for_vertex[v] );
		}

		for( int f = 0; f < nf_faces; f++ )
		{ // nf_g_q_p
			int frm = f * nf_g_q_p;
			for( int x = 0; x < theFormulas[frm].ncoor; x++ )
			{
				int i = theFormulas[frm].cp[x];

				if( nfacesSpace_for_vertex[i] == nfaces_for_vertex[i] )
				{
					nfacesSpace_for_vertex[i] *= 2;
					faces_for_vertex[i] = (int *)realloc( faces_for_vertex[i], sizeof(int) * nfacesSpace_for_vertex[i] );
				}

				faces_for_vertex[i][nfaces_for_vertex[i]] = f;
				nfaces_for_vertex[i]++;
			}	
	
		}
		
		for( int f = 0; f < nf_irr_faces; f++ )
		{ // nf_g_q_p
			int frm = f * nf_irr_pts;
			for( int x = 0; x < theIrregularFormulas[frm].ncoor; x++ )
			{
				int i = theIrregularFormulas[frm].cp[x];

				if( nfacesSpace_for_vertex[i] == nfaces_for_vertex[i] )
				{
					nfacesSpace_for_vertex[i] *= 2;
					faces_for_vertex[i] = (int *)realloc( faces_for_vertex[i], sizeof(int) * nfacesSpace_for_vertex[i] );
				}

				faces_for_vertex[i][nfaces_for_vertex[i]] = nf_faces+f;
				nfaces_for_vertex[i]++;
			}	
		}
	}*/

#ifdef GPU
	if( !oarr )
		oarr = (double *)malloc( sizeof(double) * nf_faces * nf_g_q_p );

	double alt_e;
	if( do_vertex < 0 )
	{
		alt_e = openCLEnergy( r, oarr );
		return alt_e;
	}
#endif

	double e = 0;
	double area = 0;
	double wgt = 0;

	int loop_start = 0;
	int loop_stop = nf_faces;

	if( do_vertex >= 0 )
	{
		if( par_info.my_id == BASE_TASK )
		{
			loop_stop = nfaces_for_vertex[do_vertex];
	
			av_h2 = 0;
			nav_h2 = 0;
			
			double Atot = 0;
			double Atot0 = 0;
		
			int n_p_found = 0;
		
			for( int fx = 0; fx < loop_stop; fx++ )
			{
				int f =fx;
				if( do_vertex >= 0 )
					f = faces_for_vertex[do_vertex][fx];
		
				if( debug_bit && (f != do_face) )
					continue;
		
				e += faceEnergy(f, r, puv, doMonge );
			
				if( do_vertex >= 0 )
				{
					int t = theFormulas[f*nf_g_q_p].tri;
		
					for( int px = 0; px < theTriangles[t].np; px++ )
					{
						plist[*np_found] = theTriangles[t].plist[px];
						for( int i = 0; i < *np_found; i++ )
						{
							if( plist[i] == plist[*np_found] )
							{
								printf("Dupe somehow.\n");
							}
						}
						*np_found += 1;
					} 
				}		
			}
		}
 
	}
	else
	{
		for( int fx = 0; fx < par_info.nf[surface_id]; fx++ )
		{
			int f = par_info.faces[surface_id][fx];
	
			e +=  faceEnergy( f, r, puv, doMonge ); 
		}
	}	
#if defined(GLOBAL_AREA) || defined(MICRO_KA)
	if( do_vertex < 0 && par_info.my_id == BASE_TASK )
	{
//		printf("Atot: %le Atot0: %le alphax: %lf alphay: %lf alphaz: %lf\n", Atot, Atot0, alpha_x, alpha_y, alpha_z  );
		e += 2 * 0.5 * Atot0 * KA * ((Atot-Atot0)/Atot0) * ((Atot-Atot0)/Atot0);
	}
#endif

//	av_h2 /= nav_h2;

//	printf("av_h2: %le\n", av_h2 );

//	if( do_vertex < 0 )
//		printf("e: %.14le alt_e: %.14le\n", e, alt_e );

//	exit(1);

	double ie = 0;

/*	DONE ABOVE NOW: DELETE THIS
	if( do_vertex < 0 )	
	{
		for( int f = 0; f < nf_irr_faces; f++ )
			e +=  faceEnergy( nf_faces + f, r, puv, doMonge ); 
//		ie = irregularEnergy( r );
	}
*/

	return e + ie;

}


void surface::area( double *r, int do_vertex, double *area, double *area0)
{
	if( !theFormulas )
		generatePlan();

	double alpha_x = r[3*nv];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];
	
	static int icntr = 0;
	
	VA = 0;
	VA4 = 0;
	VC = 0;
	AVC = 0;
	if( faces_for_vertex == NULL )
		computeModifiedVertices();	
/*	if( faces_for_vertex == NULL )
	{
		faces_for_vertex = (int**)malloc( sizeof(int*) * nv );
		nfaces_for_vertex = (int *)malloc( sizeof(int) * nv );
		nfacesSpace_for_vertex = (int *)malloc( sizeof(int) * nv );

		for( int v = 0; v < nv; v++ )
		{
			nfaces_for_vertex[v] = 0;
			nfacesSpace_for_vertex[v] = 10;
			faces_for_vertex[v] = (int *)malloc( sizeof(int) * nfacesSpace_for_vertex[v] );
		}

		for( int f = 0; f < nf_faces; f++ )
		{ // nf_g_q_p
			int frm = f * nf_g_q_p;
			for( int x = 0; x < theFormulas[frm].ncoor; x++ )
			{
				int i = theFormulas[frm].cp[x];

				if( nfacesSpace_for_vertex[i] == nfaces_for_vertex[i] )
				{
					nfacesSpace_for_vertex[i] *= 2;
					faces_for_vertex[i] = (int *)realloc( faces_for_vertex[i], sizeof(int) * nfacesSpace_for_vertex[i] );
				}

				faces_for_vertex[i][nfaces_for_vertex[i]] = f;
				nfaces_for_vertex[i]++;
			}	
	
		}
	}
*/
#ifdef GPU
	if( !oarr )
		oarr = (double *)malloc( sizeof(double) * nf_faces * nf_g_q_p );

	double alt_e;
	if( do_vertex < 0 )
	{
		alt_e = openCLEnergy( r, oarr );
		return alt_e;
	}
#endif

	double e = 0;
	double wgt = 0;

	int loop_start = 0;
	int loop_stop = nf_faces;

	if( do_vertex >= 0 )
		loop_stop = nfaces_for_vertex[do_vertex];

	av_h2 = 0;
	nav_h2 = 0;
	
	*area = 0;
	*area0 = 0;

	for( int fx = 0; fx < loop_stop; fx++ )
	{
		int f =fx;
		if( do_vertex >= 0 )
			f = faces_for_vertex[do_vertex][fx];

			if( debug_bit && (f != do_face) )
				continue;
		if( f < nf_faces )
		{
			for( int p = 0; p < nf_g_q_p; p++ )
			{
				double A0 = 0;
				double A = 0;
	
				if( debug_bit && (f != do_face || p != do_pt) )
					continue;
				int frm = f*nf_g_q_p+p;
				double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
				double nrm[3]={0,0,0}; 
	
				int *cp = theFormulas[f*nf_g_q_p+p].cp;
				int np = theFormulas[f*nf_g_q_p+p].ncoor;
	
				for( int p = 0; p < np; p++ )
				{
					R[0] += theFormulas[frm].r_w[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
					R[1] += theFormulas[frm].r_w[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
					R[2] += theFormulas[frm].r_w[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
					
					Ru[0] += theFormulas[frm].r_u[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
					Ru[1] += theFormulas[frm].r_u[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
					Ru[2] += theFormulas[frm].r_u[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
					
					Rv[0] += theFormulas[frm].r_v[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
					Rv[1] += theFormulas[frm].r_v[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
					Rv[2] += theFormulas[frm].r_v[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
					
					tSuu[0] += theFormulas[frm].r_uu[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
					tSuu[1] += theFormulas[frm].r_uu[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
					tSuu[2] += theFormulas[frm].r_uu[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
					
					tSuv[0] += theFormulas[frm].r_uv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
					tSuv[1] += theFormulas[frm].r_uv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
					tSuv[2] += theFormulas[frm].r_uv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
					
					tSvv[0] += theFormulas[frm].r_vv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
					tSvv[1] += theFormulas[frm].r_vv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
					tSvv[2] += theFormulas[frm].r_vv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				}
			
	
	
	
			cross( Ru, Rv, nrm );
			normalize(nrm);
	
			double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
			double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
			double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];
	
			double g = sqrt(RuRu*RvRv-RuRv*RuRv);
	
	//		printf("c1: %.14le c2: %.14le g: %.14le en: %.14le\n", c1, c2, g, en );
#ifdef FIXED_A
			double dA = theFormulas[frm].g0 * theFormulas[frm].weight;
#else
			double dA = g * theFormulas[frm].weight;
#endif	
			double dudv = 0.5;
	
			*area += dudv *dA;
			*area0 += dudv * theFormulas[frm].g0 * theFormulas[frm].weight;
	
			}
		}
		else
		{
		for( int p = 0; p < nf_irr_pts; p++ )
		{
			double A0 = 0;
			double A = 0;

			int frm = (f-nf_faces)*nf_irr_pts+p;

			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
			double nrm[3]={0,0,0}; 

			int *cp = theIrregularFormulas[frm].cp;
			int np = theIrregularFormulas[frm].ncoor;

			for( int p = 0; p < np; p++ )
			{
				R[0] += theIrregularFormulas[frm].r_w[p] * alpha_x*(r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				R[1] += theIrregularFormulas[frm].r_w[p] * alpha_y*(r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				R[2] += theIrregularFormulas[frm].r_w[p] * alpha_z*(r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				Ru[0] += theIrregularFormulas[frm].r_u[p] * alpha_x*(r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				Ru[1] += theIrregularFormulas[frm].r_u[p] * alpha_y*(r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				Ru[2] += theIrregularFormulas[frm].r_u[p] * alpha_z*(r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				Rv[0] += theIrregularFormulas[frm].r_v[p] * alpha_x*(r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				Rv[1] += theIrregularFormulas[frm].r_v[p] * alpha_y*(r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				Rv[2] += theIrregularFormulas[frm].r_v[p] * alpha_z*(r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
			}

		cross( Ru, Rv, nrm );
		normalize(nrm);

		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		double g = sqrt(RuRu*RvRv-RuRv*RuRv);

		double dA = g * theIrregularFormulas[frm].weight;
		double dudv = 1;
		*area += dudv * dA;
		*area0 += dudv * theIrregularFormulas[frm].g0 * theIrregularFormulas[frm].weight;
		}

		}		
	} 

	if( do_vertex < 0 ) {	
	for( int f = 0; f < nf_irr_faces; f++ )
	{
		double e_v[n_v_finite];
		memset( e_v, 0, sizeof(double) * n_v_finite );

		for( int p = 0; p < nf_irr_pts; p++ )
		{
			double A0 = 0;
			double A = 0;

			int frm = f*nf_irr_pts+p;

			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
			double nrm[3]={0,0,0}; 

			int *cp = theIrregularFormulas[f*nf_irr_pts+p].cp;
			int np = theIrregularFormulas[f*nf_irr_pts+p].ncoor;

			for( int p = 0; p < np; p++ )
			{
				R[0] += theIrregularFormulas[frm].r_w[p] * alpha_x*(r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				R[1] += theIrregularFormulas[frm].r_w[p] * alpha_y*(r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				R[2] += theIrregularFormulas[frm].r_w[p] * alpha_z*(r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				Ru[0] += theIrregularFormulas[frm].r_u[p] * alpha_x*(r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				Ru[1] += theIrregularFormulas[frm].r_u[p] * alpha_y*(r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				Ru[2] += theIrregularFormulas[frm].r_u[p] * alpha_z*(r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				Rv[0] += theIrregularFormulas[frm].r_v[p] * alpha_x*(r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				Rv[1] += theIrregularFormulas[frm].r_v[p] * alpha_y*(r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				Rv[2] += theIrregularFormulas[frm].r_v[p] * alpha_z*(r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
			}

		cross( Ru, Rv, nrm );
		normalize(nrm);

		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		double g = sqrt(RuRu*RvRv-RuRv*RuRv);

		double dA = g * theIrregularFormulas[frm].weight;
		double dudv = 1;
		*area += dudv * dA;
		*area0 += dudv * theIrregularFormulas[frm].g0 * theIrregularFormulas[frm].weight;

		}
		
	}
	} 

}

double surface::energyMonge( double *r, int do_vertex)
{
	if( !theFormulas )
		generatePlan();

	double alpha_x = r[3*nv];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];
	
	static int icntr = 0;
	VA = 0;
	VA4 = 0;
	VC = 0;
	AVC = 0;

	if( faces_for_vertex == NULL )
		computeModifiedVertices();
/*
	{
		faces_for_vertex = (int**)malloc( sizeof(int*) * nv );
		nfaces_for_vertex = (int *)malloc( sizeof(int) * nv );
		nfacesSpace_for_vertex = (int *)malloc( sizeof(int) * nv );

		for( int v = 0; v < nv; v++ )
		{
			nfaces_for_vertex[v] = 0;
			nfacesSpace_for_vertex[v] = 10;
			faces_for_vertex[v] = (int *)malloc( sizeof(int) * nfacesSpace_for_vertex[v] );
		}

		for( int f = 0; f < nf_faces; f++ )
		{ // nf_g_q_p
			int frm = f * nf_g_q_p;
			for( int x = 0; x < theFormulas[frm].ncoor; x++ )
			{
				int i = theFormulas[frm].cp[x];

				if( nfacesSpace_for_vertex[i] == nfaces_for_vertex[i] )
				{
					nfacesSpace_for_vertex[i] *= 2;
					faces_for_vertex[i] = (int *)realloc( faces_for_vertex[i], sizeof(int) * nfacesSpace_for_vertex[i] );
				}

				faces_for_vertex[i][nfaces_for_vertex[i]] = f;
				nfaces_for_vertex[i]++;
			}	
	
		}
	}
*/
#ifdef GPU
	if( !oarr )
		oarr = (double *)malloc( sizeof(double) * nf_faces * nf_g_q_p );

	double alt_e;
	if( do_vertex < 0 )
	{
		alt_e = openCLEnergy( r, oarr );
		return alt_e;
	}
#endif

	double e = 0;
	double area = 0;
	double wgt = 0;

	int loop_start = 0;
	int loop_stop = nf_faces;

	if( do_vertex >= 0 )
		loop_stop = nfaces_for_vertex[do_vertex];
	
	double Atot = 0;
	double Atot0 = 0;


	for( int fx = 0; fx < loop_stop; fx++ )
	{
		int f =fx;
		if( do_vertex >= 0 )
		f = faces_for_vertex[do_vertex][fx];

		if( debug_bit && (f != do_face) )
			continue;

		for( int p = 0; p < nf_g_q_p; p++ )
		{
			double A0 = 0;
			double A = 0;
			double AP = 0;

			if( debug_bit && (f != do_face || p != do_pt) )
				continue;
			int frm = f*nf_g_q_p+p;
			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
			double nrm[3]={0,0,0}; 

			int *cp = theFormulas[f*nf_g_q_p+p].cp;
			int np = theFormulas[f*nf_g_q_p+p].ncoor;

			for( int p = 0; p < np; p++ )
			{
				R[0] += theFormulas[frm].r_w[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				R[1] += theFormulas[frm].r_w[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				R[2] += theFormulas[frm].r_w[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				Ru[0] += theFormulas[frm].r_u[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				Ru[1] += theFormulas[frm].r_u[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				Ru[2] += theFormulas[frm].r_u[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				Rv[0] += theFormulas[frm].r_v[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				Rv[1] += theFormulas[frm].r_v[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				Rv[2] += theFormulas[frm].r_v[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuu[0] += theFormulas[frm].r_uu[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuu[1] += theFormulas[frm].r_uu[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuu[2] += theFormulas[frm].r_uu[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuv[0] += theFormulas[frm].r_uv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuv[1] += theFormulas[frm].r_uv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuv[2] += theFormulas[frm].r_uv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSvv[0] += theFormulas[frm].r_vv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSvv[1] += theFormulas[frm].r_vv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSvv[2] += theFormulas[frm].r_vv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			}

			// compute dudx, dudv
			
			double idet = 1.0 / (Ru[0] * Rv[1] - Ru[1] * Rv[0]); 
				
			double dudx = idet * Rv[1];
			double dudy =-idet * Rv[0];
			double dvdx =-idet * Ru[1];
			double dvdy = idet * Ru[0];

			double d2hdx2 = tSuu[2] * dudx * dudx +
					tSuv[2] * dudx * dvdx +
				        tSuv[2] * dvdx * dudx +
					tSvv[2] * dvdx * dvdx;
			
			double d2hdy2 = tSuu[2] * dudy * dudy +
					tSuv[2] * dudy * dvdy +
				        tSuv[2] * dvdy * dudy +
					tSvv[2] * dvdy * dvdy;

		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];
		
		double RuPRuP = Ru[0] * Ru[0]/alpha_x/alpha_x + Ru[1] * Ru[1]/alpha_y/alpha_y + Ru[2]*Ru[2]/alpha_z/alpha_z;
		double RuPRvP = Ru[0] * Rv[0]/alpha_x/alpha_x + Ru[1] * Rv[1]/alpha_y/alpha_y + Ru[2]*Rv[2]/alpha_z/alpha_z;
		double RvPRvP = Rv[0] * Rv[0]/alpha_x/alpha_x + Rv[1] * Rv[1]/alpha_y/alpha_y + Rv[2]*Rv[2]/alpha_z/alpha_z;

		double g_prior = sqrt(RuPRuP*RvPRvP-RuPRvP*RuPRvP);
		double g = sqrt(RuRu*RvRv-RuRv*RuRv);

		double en = 0.5 * kc * (d2hdy2 + d2hdx2) * (d2hdy2 + d2hdx2);

		if( !( en >0 || en < 1 ) )
		{
			printf("nan en.\n");
			exit(1);
		}

#ifdef FIXED_A
		double dA = theFormulas[frm].g0 * theFormulas[frm].weight;
#else
		double dA = g * theFormulas[frm].weight;
#endif	
	double dudv = 0.5;

		A += dudv * dA;
		Atot += dudv * dA;
		A0 += dudv * theFormulas[frm].g0 * theFormulas[frm].weight;
		Atot0 += dudv * theFormulas[frm].g0 * theFormulas[frm].weight;
		AP += dudv * g_prior * theFormulas[frm].weight;


		area += dudv * dA;
	
	//		if( g > 0 )
			{
				e += dudv * dA * en; 			
				VC += dudv * dA * en;
				AVC += dudv *dA * (d2hdx2+d2hdy2);

#ifdef MICRO_KA
				VA += 2 * 0.5 * A0 * micro_KA * ((AP-A0)/A0) * ((AP-A0)/A0);
				e += 2 * 0.5 * A0 * micro_KA * ((AP-A0)/A0) * ((AP-A0)/A0); 
#else
#ifndef GLOBAL_AREA
				VA += 2 * 0.5 * A0 * KA * ((A-A0)/A0) * ((A-A0)/A0);
				e += 2 * 0.5 * A0 * KA * ((A-A0)/A0) * ((A-A0)/A0); 
#endif
#endif
//			printf("E %d %.14le %.14le\n", f*nf_g_q_p+p, dudv * dA * en, 2 * 0.5 * A0 * KA * ((A-A0)/A0) * ((A-A0)/A0) );
//				printf("oarr: %.14le %.14le\n", dudv * dA * en, oarr[frm] );
//				printf("oarr: %.14le %.14le\n", tSvv[0], oarr[frm] );

			}

		}
			
	} 

#if defined(GLOBAL_AREA) || defined(MICRO_KA)
	if( do_vertex < 0 )
	{
//		printf("Atot: %le Atot0: %le\n", Atot, Atot0 );
		e += 2 * 0.5 * Atot0 * KA * ((Atot-Atot0)/Atot0) * ((Atot-Atot0)/Atot0);
	}
#endif


//	if( do_vertex < 0 )
//		printf("e: %.14le alt_e: %.14le\n", e, alt_e );

//	exit(1);
	return e;

}

double surface::irregularEnergy( double *r )
{
	if( !theFormulas )
		generatePlan();
	double alpha_x = r[3*nv];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];
//	printf("%lf %lf %lf\n", alpha_x, alpha_y, alpha_z );
	static int icntr = 0;
//	VA = 0;
//	VC = 0;

	double e = 0;
	double area = 0;
	double wgt = 0;

	for( int f = 0; f < nf_irr_faces; f++ )
	{
		double e_v[n_v_finite];
		memset( e_v, 0, sizeof(double) * n_v_finite );

		for( int p = 0; p < nf_irr_pts; p++ )
		{
			double A0 = 0;
			double A = 0;

			int frm = f*nf_irr_pts+p;

			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
			double nrm[3]={0,0,0}; 

			int *cp = theIrregularFormulas[f*nf_irr_pts+p].cp;
			int np = theIrregularFormulas[f*nf_irr_pts+p].ncoor;

			for( int p = 0; p < np; p++ )
			{
				R[0] += theIrregularFormulas[frm].r_w[p] * alpha_x*(r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				R[1] += theIrregularFormulas[frm].r_w[p] * alpha_y*(r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				R[2] += theIrregularFormulas[frm].r_w[p] * alpha_z*(r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				Ru[0] += theIrregularFormulas[frm].r_u[p] * alpha_x*(r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				Ru[1] += theIrregularFormulas[frm].r_u[p] * alpha_y*(r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				Ru[2] += theIrregularFormulas[frm].r_u[p] * alpha_z*(r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				Rv[0] += theIrregularFormulas[frm].r_v[p] * alpha_x*(r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				Rv[1] += theIrregularFormulas[frm].r_v[p] * alpha_y*(r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				Rv[2] += theIrregularFormulas[frm].r_v[p] * alpha_z*(r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSuu[0] += theIrregularFormulas[frm].r_uu[p] * alpha_x*(r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSuu[1] += theIrregularFormulas[frm].r_uu[p] * alpha_y*(r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSuu[2] += theIrregularFormulas[frm].r_uu[p] * alpha_z*(r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSuv[0] += theIrregularFormulas[frm].r_uv[p] * alpha_x*(r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSuv[1] += theIrregularFormulas[frm].r_uv[p] * alpha_y*(r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSuv[2] += theIrregularFormulas[frm].r_uv[p] * alpha_z*(r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSvv[0] += theIrregularFormulas[frm].r_vv[p] * alpha_x*(r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSvv[1] += theIrregularFormulas[frm].r_vv[p] * alpha_y*(r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSvv[2] += theIrregularFormulas[frm].r_vv[p] * alpha_z*(r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
			}

		cross( Ru, Rv, nrm );
		normalize(nrm);

		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		double g = sqrt(RuRu*RvRv-RuRv*RuRv);

		double nsuu = tSuu[0] * nrm[0] + tSuu[1] * nrm[1] + tSuu[2] * nrm[2];
		double nsuv = tSuv[0] * nrm[0] + tSuv[1] * nrm[1] + tSuv[2] * nrm[2];
		double nsvv = tSvv[0] * nrm[0] + tSvv[1] * nrm[1] + tSvv[2] * nrm[2];

		double Stot = (nsuu * RvRv + nsvv * RuRu -2*nsuv*RuRv)/(g*g);

//		Int[i*n+j] = 0.5 * kc * Stot * Stot * g; 
//		AInt[i*n+j] = g;

		double Sop[4] = { 1.0/(g*g) * (nsuu * RvRv  - nsuv * RuRv), (1.0/(g*g)) * ( nsuv*RvRv-nsvv*RuRv),
				  1.0/(g*g) * (nsuv * RuRu  - nsuu * RuRv), (1.0/(g*g)) * ( nsvv*RuRu-nsuv*RuRv) };

		double a = Sop[0];		
		double b = Sop[1];
		double c = Sop[2];
		double d = Sop[3];

		double c0 = theIrregularFormulas[frm].c0;
		double c1 = -0.5*(a+d-sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
		double c2 = -0.5*(a+d+sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
			
//		printf("e1: %lf e2: %lf\n", e1, e2 );
		double en = 0.5 * kc_irr * (c1+c2-c0 ) * (c1+c2-c0) + kg * c1 * c2;


		if( !( en >0 || en < 1 ) )
		{
			printf("nan en.\n");
			exit(1);
		}

//		printf("c1: %.14le c2: %.14le g: %.14le en: %.14le\n", c1, c2, g, en );
#ifdef FIXED_A
		double dA = theIrregularFormulas[frm].g0 * theIrregularFormulas[frm].weight;
#else
		double dA = g * theIrregularFormulas[frm].weight;
#endif	

		// factor of 1/2 is accounted for in the weights.
		double dudv = 1.0;

		A += dudv * dA;
		A0 += dudv * theIrregularFormulas[frm].g0 * theIrregularFormulas[frm].weight;

//		printf("A: %le en: %le running e: %le weight: %le g: %le\n", A, en, e, theIrregularFormulas[frm].weight, g );
		area += dudv * dA;
	
	//		if( g > 0 )
			{
				e += dudv * dA * en; 			
				VC += dudv * dA * en;
				AVC += dudv *dA * (c1+c2);
				VA += 2 * 0.5 * A0 * KA * ((A-A0)/A0) * ((A-A0)/A0);
				e += 2 * 0.5 * A0 * KA * ((A-A0)/A0) * ((A-A0)/A0); 
			}

//			printf("p: %d n_v_finite: %d en: %le dA: %le\n", p, n_v_finite, en, dA );
//			e_v[p/n_u_gauss] += dudv * dA * en;

/*			printf("u: %le v: %le g: %le c1: %le c2: %le ce: %le w: %le RuRu: %le RuRv: %le RvRv: %le ",
				theIrregularFormulas[frm].orig_u,
				theIrregularFormulas[frm].orig_v,
				g, c1, c2, en, theIrregularFormulas[frm].weight, RuRu, RuRv, RvRv ); 
			for( int p = 0; p < theIrregularFormulas[frm].ncoor; p++ )
				printf(" %le", theIrregularFormulas[frm].r_v[p] );
			printf("\n");
*/
							
	/*		else
			{
				printf("zero'ing g.\n");
				e += 1e5;
			}*/
		}
		
//		for( int v = 0; v < n_v_finite; v++ )
//			printf("%d %le\n", v, e_v[v] );
			
//		exit(1);
	} 



//	if( do_vertex < 0 )
//		printf("e: %.14le alt_e: %.14le\n", e, alt_e );

//	exit(1);
//	printf("irregular energy: %le\n", e );
	return e;

}

void surface::generatePlan( void )
{
	faces_for_vertex=NULL;
	nfaces_for_vertex=NULL;
	nfacesSpace_for_vertex=NULL;
	// each quadrature point has a formula for computing the position, tangents, and their derivatives wrt the control points.

	double *Aevec[1+MAX_INV_VALENCE];
	double *Aevec_R[1+MAX_INV_VALENCE];
	double *Aval[1+MAX_INV_VALENCE];
	double *Atot[1+MAX_INV_VALENCE];
	int high_ev[1+MAX_INV_VALENCE];			
	
	for( int val = 4; val <= MAX_INV_VALENCE; val++ )
	{
		double w = 3.0/(8*val);
		double w6 = 3.0/(8*6);
		int ncoords_base  = 1 + val + 5; 
		// assemble the matrix:

		int ncoords_extra = 6;
		Atot[val] = (double *)malloc( sizeof(double) * (ncoords_base+ncoords_extra) * ncoords_base ); 
	
		double *A = Atot[val];

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

		
/*
		A[2*ncoords_base+0]   = (3.0/8.0);	
		A[2*ncoords_base+2]   = (3.0/8.0);	
		A[2*ncoords_base+1]   = (1.0/8.0);	
		A[2*ncoords_base+3]   = (1.0/8.0);	

		// 0 3 2 4
		A[3*ncoords_base+0] = (3.0/8.0);	
		A[3*ncoords_base+3] = (3.0/8.0);	
		A[3*ncoords_base+2] = (1.0/8.0);	
		A[3*ncoords_base+4] = (1.0/8.0);	
		
		// 0 4 3 5
		A[4*ncoords_base+0] = (3.0/8.0);	
		A[4*ncoords_base+4] = (3.0/8.0);	
		A[4*ncoords_base+3] = (1.0/8.0);	
		A[4*ncoords_base+5] = (1.0/8.0);	
		
		if( val == 5 )
		{
			// 0 5 4 1
			A[5*ncoords_base+0] = (3.0/8.0);	
			A[5*ncoords_base+5] = (3.0/8.0);	
			A[5*ncoords_base+4] = (1.0/8.0);	
			A[5*ncoords_base+1] = (1.0/8.0);	
		}
		else if( val == 7 )
		{
			A[5*ncoords_base+0] = (3.0/8.0);	
			A[5*ncoords_base+5] = (3.0/8.0);	
			A[5*ncoords_base+4] = (1.0/8.0);	
			A[5*ncoords_base+6] = (1.0/8.0);	
			
			A[6*ncoords_base+0] = (3.0/8.0);	
			A[6*ncoords_base+6] = (3.0/8.0);	
			A[6*ncoords_base+5] = (1.0/8.0);	
			A[6*ncoords_base+7] = (1.0/8.0);	
			
			A[7*ncoords_base+0] = (3.0/8.0);	
			A[7*ncoords_base+7] = (3.0/8.0);	
			A[7*ncoords_base+6] = (1.0/8.0);	
			A[7*ncoords_base+1] = (1.0/8.0);	
		}
*/



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
			
#ifdef PRINT_A
			printf("A = {");
			for( int p = 0; p < ncoords_base; p++ )
			{
				printf("{");
				for( int q = 0; q < ncoords_base; q++ )
				{
					printf("%.14lf", A[p*ncoords_base+q] );
					if( q != ncoords_base-1)
						printf(",");
				}
				printf("}");
				if( p != ncoords_base )
					printf(",");
			}	
			printf("};");
#endif
		Aevec[val] = (double *)malloc( sizeof(double) * ncoords_base * ncoords_base );
		Aevec_R[val] = (double *)malloc( sizeof(double) * ncoords_base * ncoords_base );
		Aval[val] = (double *)malloc( sizeof(double) * ncoords_base );
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

		dgeev( &jobvl, &jobvr, &N, a_copy, &LDA, WR, WI, Aevec[val], &LDVL, VR, &LDVR, work, &LWORK, &info );
		memcpy( Aevec_R[val], Aevec[val], sizeof(double) * ncoords_base * ncoords_base );
		int ipiv[N];
		dgetrf( &N, &N, Aevec_R[val], &N, ipiv, &info );
		dgetri( &N, Aevec_R[val], &N, ipiv, work, &LWORK, &info );

/*
		for( int i = 0; i < ncoords_base; i++ )
		for( int j = 0; j < ncoords_base; j++ )
			a_copy[j*ncoords_base+i] = Aevec[val][i*ncoords_base+j];		
		memcpy( Aevec[val], a_copy, sizeof(double) * ncoords_base * ncoords_base );	
*/
		memcpy( Aval[val], WR, sizeof(double) * N );


		for( int i = 0; i < ncoords_base; i++ )
		for( int j = 0; j < ncoords_base; j++ )
			a_copy[j*ncoords_base+i] = Aevec_R[val][i*ncoords_base+j];
		memcpy( Aevec_R[val], a_copy, sizeof(double) * ncoords_base * ncoords_base );

		double max_ev = -1e10;
		for( int x = 0; x < ncoords_base; x++ )
		{
			if( Aval[val][x] > max_ev )
			{
				max_ev = Aval[val][x];
				high_ev[val] = x;
			}
		}


		if( info != 0 )
		{
			printf("DGEEV error valence %d.\n", val);
			exit(1);
		}
		else
		{
#ifdef EV_STRUCTURE_DEBUG
			printf("valence: %d\n", val );
			for( int i = 0; i < ncoords_base; i++ )
				printf("\teval %d %.14le %.14le\n", i, WR[i], WI[i] );
			printf("high ev vector:\n");

			for( int i = 0; i < ncoords_base; i++ )
				printf("\t%d %le\n", i, Aevec[val][high_ev[val]*ncoords_base+i] * Aevec_R[val][high_ev[val]*ncoords_base+i] );
#endif
		}

		Aval[val][high_ev[val]] = 1;

		free(work);
		free(a_copy);
	}

	int ni = 8;

	
	int nfaces = 0;

	
#ifdef G_Q_DEBUG
			int nuse = 10;
			int n_g_q_p = nuse *( nuse+1 ) /2;
	
			double **pts;

			pts = (double **)malloc( sizeof(double *) * n_g_q_p );
			double *weights = (double *)malloc( sizeof(double) * n_g_q_p );
			double eps = 1e-6;

			int pcntr = 0;
			for( int iu = 0; iu < nuse; iu++ )
			{
				for( int iv = 0; iv < nuse-iu; iv++ )
				{
					pts[pcntr] = (double *)malloc( sizeof(double) * 2 );

					pts[pcntr][0] = pow(2, -(nuse-iu)/2.0);// / (double)nuse;
					pts[pcntr][1] = pow(2, -(nuse-iv)/2.0);// / (double)nuse;
					//pts[pcntr][1] = eps*iv / (double)nuse;
					weights[pcntr] = 1.0 / (double)n_g_q_p;
					pcntr++;	
				}
			}

/*
			int n_g_q_p = 15;
			double pts[15][2] = {
				{0.1,0.1},
				{0.3,0.1},
				{0.5,0.1},
				{0.7,0.1},
				{0.9,0.1},
				{0.1,0.3},
				{0.3,0.3},
				{0.5,0.3},
				{0.7,0.3},
				{0.1,0.5},
				{0.3,0.5},
				{0.5,0.5},
				{0.1,0.7},
				{0.3,0.7},
				{0.1,0.9}
			};

			double del = 1.0 / 6.0;


			double weights[] = { 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
			};
*/
#elif defined(G_Q_DEBUG_2)
			int n_g_q_p = 2;
			double pts[][2] = {
				{1.0/1024.0, 1.0/1024.0},
				{1.0/1024.000001, 1.0/1024.000001}
			};
			double weights[] = { 0.5,0.5 };

#elif defined(G_Q_P_15)
			int n_g_q_p = 15;
			double pts[15][2] = {
				{0.1,0.1},
				{0.3,0.1},
				{0.5,0.1},
				{0.7,0.1},
				{0.9,0.1},
				{0.1,0.3},
				{0.3,0.3},
				{0.5,0.3},
				{0.7,0.3},
				{0.1,0.5},
				{0.3,0.5},
				{0.5,0.5},
				{0.1,0.7},
				{0.3,0.7},
				{0.1,0.9}
			};

			double del = 1.0 / 6.0;


			double weights[] = { 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
				1.0/15, 
			};
#elif defined(G_Q_P_PRINT)
			int n_g_q_p = 15;
			double pts[][2] = {
				{0, 0},
				{0.25, 0},
				{0.5, 0},
				{0.75, 0},
				{1.0, 0},
				{0, 0.25},
				{0.25, 0.25},
				{0.5, 0.25},
				{0.75, 0.25},
				{0, 0.5},
				{0.25, 0.5},
				{0.5, 0.5},
				{0, 0.75},
				{0.25,0.75},
				{0, 1.00}
			};
			double weights[] = { 0.0666666667,0.0666666667,0.0666666667,0.0666666667,0.0666666667,0.0666666667,0.0666666667,0.0666666667,0.0666666667,0.0666666667 , 0.0666666667, 0.0666666667, 0.0666666667, 0.066666667, 0.0666666667};
		
#elif defined(G_Q_P_1)
			int n_g_q_p = 1;
			double pts[][2] = {
				{1.0/3.0, 1.0/3.0}
			};
			double weights[] = { 1.0 };
#elif defined(G_Q_P_3)

			int n_g_q_p = 3;
			double pts[][2] = {
				{1.0/6.0, 1.0/6.0},
				{4.0/6.0, 1.0/6.0},
				{1.0/6.0, 4.0/6.0}
			};
			double weights[] = { 1.0/3.0, 1.0/3.0, 1.0/3.0 };

/*			int n_g_q_p = 3;
			double pts[][2] = {
				{1.0/4.0, 1.0/4.0},
				{1.0/2.0, 1.0/4.0},
				{1.0/4.0, 1.0/2.0}
			};
			double weights[] = { 1.0/3.0, 1.0/3.0, 1.0/3.0 };
*/
#else
			int n_g_q_p = 4;
			double pts[][2] = {
				{1.0/3.0, 1.0/3.0},
				{1.0/5.0, 1.0/5.0},
				{1.0/5.0, 3.0/5.0},
				{3.0/5.0, 1.0/5.0},
			};
			double weights[] = { -27.0/48.0, 25.0/48.0, 25.0/48.0, 25.0/48.0 };

#endif

	nf_irr_faces = 0;
	nf_irr_pts = (n_u_gauss) * n_v_finite; // one point extra for the triangular remainder region.

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

			if( theVertices[i].valence == 6 )
				nfaces++;
			else
				nf_irr_faces++;
		}
	}
	
	nf_faces = nfaces;
	nf_g_q_p = n_g_q_p;

	int nfi = 0;
	int nf = 0;
	int nfSpace = n_g_q_p * nfaces;
	int nef = 0;
	edgeFormulas = (formula *)malloc( sizeof(formula) * 2 * (nfaces+nf_irr_faces) );
	theFormulas = (formula *)malloc( sizeof(formula) * n_g_q_p * nfaces );
	theIrregularFormulas = (formula *)malloc( sizeof(formula) * nf_irr_pts * nf_irr_faces );

	disable_PBC = 1;

	int trigger=0;
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

			int tri = theVertices[i].faces[e];
			
			int ej = theVertices[i].edge_rev[e];
			int ek = theVertices[i].edge_rev[ep1];


			if( val == 6 && ( valj != 6 || valk != 6 ) )
			{

			printf("%d valence %d, %d valence %d, %d valence %d.\n",
				i, theVertices[i].valence,
				j, theVertices[j].valence,
				k, theVertices[k].valence );
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


			if( val != 6 ) // irregular, use custom plan.
			{
				int n_pt_quadrature = n_u_gauss;

				double k_vals[11] = 
				{
					0, 0, 0, 0, // unused
					0.375, // 4
					4.52254248593737e-01, // 5 
					0.5, // 6
					5.30872450464684e-01, // 7
					5.51776695296637e-01, // 8
					5.66511110779745e-01, // 9					
					5.77254248593736e-01 // 10	
				};


/*				double v_pts[] = { 0, 0.25, 0.5, 0.75 };

				// evaluation points.
				double u_pts[11][4][2] = 
				{
					{ },
					{ },
					{ },
					{ },
					{ },
					{ 
						{0.258707503209, 0.808802027429 }, // 0
						{0.170390766908, 0.599247199742 }, // 1
						{0.109372868997, 0.397287127896 }, // 2
						{0.0535585801263,    0.197828866407   }, // 3
					},
				};

				// Gauss weights				
				double w_pts[11][4][2] = 
				{
					{ },
					{ },
					{ },
					{ },
					{ },
					{ // 5
						{0.3461565082631408, 0.42928490375 }, // 0
						{0.302793914785,     0.342887865727    }, // 1
						{0.221382686717,     0.236850604628    }, // 2
						{0.118448886217,     0.121899207205    }, // 3
					}
				};
*/

				double k = k_vals[val];
				double v_pts[n_v_finite];
				double w_v_pts[n_v_finite];
	
				get_v_pts( v_pts, w_v_pts, n_v_finite, k );

				for( int v_pt = 0; v_pt < n_v_finite; v_pt++ )
				{
//					printf("v[%d]: %le w: %le\n", v_pt, v_pts[v_pt], w_v_pts[v_pt] );
//					double fv = v_pt / (double)n_v_finite;

					double fv = v_pts[v_pt];
					double v_factor = - (-1+pow(k, 4.0 * log2(fv))*fv*fv*fv*fv*fv) * log(2.0) / (log(32.0)+4.0*log(k));
					double v_weight = w_v_pts[v_pt] / v_factor; 

					double u_pts[n_u_gauss];
					double w_pts[n_u_gauss];

					//double ulim = 1-(fv+dv); 
					
					get_u_pts( u_pts, w_pts, n_u_gauss, k, fv );
			

//					printf("v: %lf\n", fv );
					for( int u_pt = 0; u_pt < n_u_gauss; u_pt++ )
					{
						fv = v_pts[v_pt];
//						printf("u: %lf w: %lf\n", u_pts[u_pt], w_pts[u_pt] );
						double fu = u_pts[u_pt];

						double factor = pow( 1.0 / (fu+fv) * pow( k, -log2( fu+fv) ), 4.0 );
						double u_weight = w_pts[u_pt] * factor; 
						double weight = v_weight * u_weight;
						double fw = 1-fu-fv;
		
						double o_fu = fu;
						double o_fv = fv;

						theIrregularFormulas[nfi].ncoor = ncoords_base;
						theIrregularFormulas[nfi].cp = (int *)malloc( sizeof(int) * ncoords_base ); 
		
						theIrregularFormulas[nfi].r_pbc = (double *)malloc( sizeof(double) * 3 * ncoords_base );
						theIrregularFormulas[nfi].r_w = (double *)malloc( sizeof(double) * ncoords_base ); 
						theIrregularFormulas[nfi].r_u = (double *)malloc( sizeof(double) * ncoords_base ); 
						theIrregularFormulas[nfi].r_v = (double *)malloc( sizeof(double) * ncoords_base ); 
						theIrregularFormulas[nfi].r_uu = (double *)malloc( sizeof(double) * ncoords_base ); 
						theIrregularFormulas[nfi].r_uv = (double *)malloc( sizeof(double) * ncoords_base ); 
						theIrregularFormulas[nfi].r_vv = (double *)malloc( sizeof(double) * ncoords_base ); 
						theIrregularFormulas[nfi].weight = weight;
//						printf("fu: %le fv: %le w: %le\n", fu, fv, weight );		

						memcpy( theIrregularFormulas[nfi].cp, point_list, sizeof(int) * ncoords_base );
						memset( theIrregularFormulas[nfi].r_pbc, 0, sizeof(double) * 3 * ncoords_base );
		
						for( int p = 0; p < ncoords_base; p++ )
						{
							int *cp = point_list;
		
							double *r1 = theVertices[cp[p]].r;
							double *r0 = theVertices[cp[0]].r;
			
							double dr[3] = { r1[0] - r0[0], r1[1] - r0[1], r1[2] - r0[2] };
							double add[3]={0,0,0};
		
							MinImage3D( dr, PBC_vec, add );
		
							if( add[0]*add[0]+add[1]*add[1]+add[2]*add[2] > 1e-9 ) disable_PBC = 0;

							theIrregularFormulas[nfi].r_pbc[3*p+0] = add[0]*PBC_vec[0][0] + add[1] * PBC_vec[1][0] + add[2] * PBC_vec[2][0];
							theIrregularFormulas[nfi].r_pbc[3*p+1] = add[0]*PBC_vec[0][1] + add[1] * PBC_vec[1][1] + add[2] * PBC_vec[2][1];
							theIrregularFormulas[nfi].r_pbc[3*p+2] = add[0]*PBC_vec[0][2] + add[1] * PBC_vec[1][2] + add[2] * PBC_vec[2][2];
						}
		
						double scale_f_u = 1.0;
						double scale_f_v = 1.0;
		
						double coord_map[12 * ncoords_base];
						memset( coord_map, 0, sizeof(double) * 12 * ncoords_base );
	
						double f_domain = -log10(fu+fv)/log10(2.0);
						if( fu+fv >= 1.0 -1e-9 )
							f_domain = 1;
						int domain = lround(ceil(f_domain));
						double pow2 = pow( 2.0, domain-1.0 );
						
						scale_f_u *= pow2;
						scale_f_v *= pow2;
	
						fu *= pow2;
						fv *= pow2;
	
						int pcycle[12];
	
						if( fu > 0.5 )
						{
							fv = 2 * fv;
							fu = 2.0 * fu-1.0;
						
							scale_f_u *= 2;
							scale_f_v *= 2;
	
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
						else if( fv > 0.5 )
						{
							fu = 2 * fu;
							fv = 2.0 * fv-1.0;
	
							scale_f_u *= 2;
							scale_f_v *= 2;
						
							
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
						else
						{
							fv = 1 - 2 * fv;
							fu = 1 - 2 * fu;
							
							scale_f_u *= -2;
							scale_f_v *= -2;
							
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
	
						
						for( int ex = 0; ex < ncoords_base; ex++ )
						for( int j = 0; j < ncoords_base; j++ )
						{
							A_prev[ex*ncoords_base+j] = 0;
	
							for( int i = 0; i < ncoords_base; i++ )
								A_prev[ex*ncoords_base+j] += Aevec[val][i*ncoords_base+j] * pow( Aval[val][i], domain-1 ) * Aevec_R[val][i*ncoords_base+ex];
						}
						
						
	
						double A_use[(ncoords_base+ncoords_extra)*ncoords_base];
						memset( A_use, 0, sizeof(double) * (ncoords_base+ncoords_extra)*ncoords_base );
						for( int f = 0; f < ncoords_base+ncoords_extra; f++ )
						{
							for( int ii = 0; ii < ncoords_base; ii++ )
							{
								A_use[f*ncoords_base+ii] = 0;
								
								for( int e = 0; e < ncoords_base; e++ )
									A_use[f*ncoords_base+ii] += Atot[val][f*ncoords_base+e] * A_prev[e*ncoords_base+ii];
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
	
					
	
					// get the coefficients.
						
						double u = fu;
						double v = fv;
						double w = 1 - u - v;
	
						double u2 = u*u;
						double u3 = u*u*u;
						double u4 = u*u*u*u;
						
						double v2 = v*v;
						double v3 = v*v*v;
						double v4 = v*v*v*v;
				
						double w2 = w*w;
						double w3 = w*w*w;
						double w4 = w*w*w*w;
							
						
						// 8 : 0
						// 7 : 1
						// 4 : 2
						// 5 : 3
						// 9 : 4
						// 12 : 5
						// 11 : 6
						// 10 : 7
						// 6 : 8
						// 3 : 9
						// 1 : 10
						// 2 : 11
	
		
	
						double n1 = (1.0/12.0)*(u4+2*u3*v); 
						double n2 = (1.0/12.0)*(u4+2*u3*w); 
						double n3 = (1.0/12.0)*(u4+2*u3*w+6*u3*v+6*u2*v*w+12*u2*v2+6*u*v2*w+6*u*v3+2*v3*w+v4); 
						double n4 = (1.0/12.0)*(6*u4+24*u3*w+24*u2*w2+8*u*w3+w4+24*u3*v+60*u2*v*w+36*u*v*w2+6*v*w3+24*u2*v2+36*u*v2*w+12*v2*w2+8*u*v3+6*v3*w+v4); 
						double n5 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+2*u3*v+6*u2*v*w+6*u*v*w2+2*v*w3); 
						double n6 = (1.0/12.0)*(2*u*v3+v4); 
						double n7 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+8*u3*v+36*u2*v*w+36*u*v*w2+8*v*w3+24*u2*v2+60*u*v2*w+24*v2*w2+24*u*v3+24*v3*w+6*v4); 
						double n8 = (1.0/12.0)*(u4+8*u3*w+24*u2*w2+24*u*w3+6*w4+6*u3*v+36*u2*v*w+60*u*v*w2+24*v*w3+12*u2*v2+36*u*v2*w+24*v2*w2+6*u*v3+8*v3*w+v4); 
						double n9 = (1.0/12.0)*(2*u*w3+w4); 
						double n10 = (1.0/12.0)*(2*v3*w+v4); 
						double n11 = (1.0/12.0)*(2*u*w3+w4+6*u*v*w2+6*v*w3+6*u*v2*w+12*v2*w2+2*u*v3+6*v3*w+v4); 
						double n12 = (1.0/12.0)*(w4+2*v*w3);
	
						double n1_alt = Power(u,4)/12. + (Power(u,3)*v)/6.;		
						double n2_alt = Power(u,4)/12. + (Power(u,3)*w)/6.;
						double n3_alt = Power(u,4)/12. + (Power(u,3)*v)/2. + Power(u,2)*Power(v,2) + (u*Power(v,3))/2. + Power(v,4)/12. + (Power(u,3)*w)/6. + (Power(u,2)*v*w)/2. + (u*Power(v,2)*w)/2. + (Power(v,3)*w)/6.;
						double n4_alt = Power(u,4)/2. + 2*Power(u,3)*v + 2*Power(u,2)*Power(v,2) + (2*u*Power(v,3))/3. + Power(v,4)/12. + 2*Power(u,3)*w + 5*Power(u,2)*v*w + 3*u*Power(v,2)*w + (Power(v,3)*w)/2. + 2*Power(u,2)*Power(w,2) + 3*u*v*Power(w,2) + Power(v,2)*Power(w,2) + (2*u*Power(w,3))/3. + (v*Power(w,3))/2. + Power(w,4)/12.;
						double n5_alt = Power(u,4)/12. + (Power(u,3)*v)/6. + (Power(u,3)*w)/2. + (Power(u,2)*v*w)/2. + Power(u,2)*Power(w,2) + (u*v*Power(w,2))/2. + (u*Power(w,3))/2. + (v*Power(w,3))/6. + Power(w,4)/12.;
						double n6_alt = (u*Power(v,3))/6. + Power(v,4)/12.;
						double n7_alt = Power(u,4)/12. + (2*Power(u,3)*v)/3. + 2*Power(u,2)*Power(v,2) + 2*u*Power(v,3) + Power(v,4)/2. + (Power(u,3)*w)/2. + 3*Power(u,2)*v*w + 5*u*Power(v,2)*w + 2*Power(v,3)*w + Power(u,2)*Power(w,2) + 3*u*v*Power(w,2) + 2*Power(v,2)*Power(w,2) + (u*Power(w,3))/2. + (2*v*Power(w,3))/3. + Power(w,4)/12.;
						double n8_alt = Power(u,4)/12. + (Power(u,3)*v)/2. + Power(u,2)*Power(v,2) + (u*Power(v,3))/2. + Power(v,4)/12. + (2*Power(u,3)*w)/3. + 3*Power(u,2)*v*w + 3*u*Power(v,2)*w + (2*Power(v,3)*w)/3. + 2*Power(u,2)*Power(w,2) + 5*u*v*Power(w,2) + 2*Power(v,2)*Power(w,2) + 2*u*Power(w,3) + 2*v*Power(w,3) + Power(w,4)/2.;
						double n9_alt = (u*Power(w,3))/6. + Power(w,4)/12.;
						double n10_alt = Power(v,4)/12. + (Power(v,3)*w)/6.;
						double n11_alt =(u*Power(v,3))/6. + Power(v,4)/12. + (u*Power(v,2)*w)/2. + (Power(v,3)*w)/2. + (u*v*Power(w,2))/2. + Power(v,2)*Power(w,2) + (u*Power(w,3))/6. + (v*Power(w,3))/2. + Power(w,4)/12.;
						double n12_alt = (v*Power(w,3))/6. + Power(w,4)/12.;
		
						double du_1 = Power(u,3)/3. + (Power(u,2)*v)/2.;
						double du_2 = Power(u,2)/2. - Power(u,3)/3. - (Power(u,2)*v)/2.; 	
						double du_3 =Power(u,2)/2. - Power(u,3)/3. + u*v - (Power(u,2)*v)/2. + Power(v,2)/2. - Power(v,3)/6.;
						double du_4 = 0.3333333333333333 + u - Power(u,2) - Power(u,3)/3. + v/2. - u*v - (Power(u,2)*v)/2. - Power(v,2) + Power(v,3)/3.;	
						double du_5 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. - v/2. + Power(u,2)*v + Power(v,2)/2. - Power(v,3)/6.;
						double du_6 = Power(v,3)/6.;	
						double du_7 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. + v/2. - 2*u*v + Power(u,2)*v - Power(v,2)/2. - Power(v,3)/6.;	
						double du_8 = -2*u + 2*Power(u,2) - Power(u,3)/3. - v + 2*u*v - (Power(u,2)*v)/2. + Power(v,2) - Power(v,3)/6.;	
						double du_9 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. + v/2. - (Power(u,2)*v)/2. - Power(v,2)/2. + Power(v,3)/6.;	
						double du_10 = -Power(v,3)/6.;
						double du_11 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. - v/2. + u*v - (Power(u,2)*v)/2. + Power(v,3)/3.;	
						double du_12 = 	-0.3333333333333333 + u - Power(u,2) + Power(u,3)/3. + v/2. - u*v + (Power(u,2)*v)/2. - Power(v,3)/6.;
	
						double dv_1 = Power(u,3)/6.;
						double dv_2 = -Power(u,3)/6.;
						double dv_3 = Power(u,2)/2. - Power(u,3)/6. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
						double dv_4 = 0.16666666666666666 + u/2. - Power(u,2)/2. - Power(u,3)/6. - 2*u*v - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
						double dv_5 = -0.16666666666666666 - u/2. + Power(u,3)/3. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
						double dv_6 = (u*Power(v,2))/2. + Power(v,3)/3.;
						double dv_7 = 0.3333333333333333 + u/2. - Power(u,2) + Power(u,3)/3. + v - u*v - Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
						double dv_8 = -u + Power(u,2) - Power(u,3)/6. - 2*v + 2*u*v + 2*Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
						double dv_9 = -0.3333333333333333 + u/2. - Power(u,3)/6. + v - u*v - Power(v,2) + (u*Power(v,2))/2. + Power(v,3)/3.;
						double dv_10 = Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
						double dv_11 = 0.16666666666666666 - u/2. + Power(u,2)/2. - Power(u,3)/6. - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
						double dv_12 = -0.16666666666666666 + u/2. - Power(u,2)/2. + Power(u,3)/6. + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		
						double d_uu_1 = Power(u,2) + u*v;
						double d_uu_2 = u - Power(u,2) - u*v;
						double d_uu_3 = u - Power(u,2) + v - u*v;
						double d_uu_4 = 1 - 2*u - Power(u,2) - v - u*v;
						double d_uu_5 = -2*u + 2*Power(u,2) + 2*u*v;
						double d_uu_6 = 0;
						double d_uu_7 = -2*u + 2*Power(u,2) - 2*v + 2*u*v;
						double d_uu_8 = -2 + 4*u - Power(u,2) + 2*v - u*v;
						double d_uu_9 = u - Power(u,2) - u*v;
						double d_uu_10 = 0;
						double d_uu_11 = u - Power(u,2) + v - u*v;
						double d_uu_12 = 1 - 2*u + Power(u,2) - v + u*v;
						
						double d_uv_1 = Power(u,2)/2.;
						double d_uv_2 = -Power(u,2)/2.;
						double d_uv_3 = u - Power(u,2)/2. + v - Power(v,2)/2.;
						double d_uv_4 = 0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
						double d_uv_5 = -0.5 + Power(u,2) + v - Power(v,2)/2.;
						double d_uv_6 = Power(v,2)/2.;
						double d_uv_7 = 0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
						double d_uv_8 = -1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
						double d_uv_9 = 0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
						double d_uv_10 = -Power(v,2)/2.;
						double d_uv_11 = -0.5 + u - Power(u,2)/2. + Power(v,2);
						double d_uv_12 = 0.5 - u + Power(u,2)/2. - Power(v,2)/2.;
						
						double d_vv_1 = 0;
						double d_vv_2 = 0;
						double d_vv_3 = u + v - u*v - Power(v,2);
						double d_vv_4 = -2*u - 2*v + 2*u*v + 2*Power(v,2);
						double d_vv_5 = u + v - u*v - Power(v,2);
						double d_vv_6 = u*v + Power(v,2);
						double d_vv_7 = 1 - u - 2*v - u*v - Power(v,2);
						double d_vv_8 = -2 + 2*u + 4*v - u*v - Power(v,2);
						double d_vv_9 = 1 - u - 2*v + u*v + Power(v,2);
						double d_vv_10 = v - u*v - Power(v,2);
						double d_vv_11 = -2*v + 2*u*v + 2*Power(v,2);
						double d_vv_12 = v - u*v - Power(v,2);
			
						double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
						double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
						double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
						
						double ceff_map_duu[12] = { d_uu_8, d_uu_7, d_uu_4, d_uu_5, d_uu_9, d_uu_12, d_uu_11, d_uu_10, d_uu_6, d_uu_3, d_uu_1, d_uu_2 };
						double ceff_map_duv[12] = { d_uv_8, d_uv_7, d_uv_4, d_uv_5, d_uv_9, d_uv_12, d_uv_11, d_uv_10, d_uv_6, d_uv_3, d_uv_1, d_uv_2 };
						double ceff_map_dvv[12] = { d_vv_8, d_vv_7, d_vv_4, d_vv_5, d_vv_9, d_vv_12, d_vv_11, d_vv_10, d_vv_6, d_vv_3, d_vv_1, d_vv_2 };
	
	
						for( int x = 0; x < ncoords_base; x++ )
						{
							theIrregularFormulas[nfi].r_w[x] = 0;
							theIrregularFormulas[nfi].r_u[x] = 0;
							theIrregularFormulas[nfi].r_v[x] = 0;
							theIrregularFormulas[nfi].r_uu[x] = 0;
							theIrregularFormulas[nfi].r_uv[x] = 0;
							theIrregularFormulas[nfi].r_vv[x] = 0;
			
							// ceff_map: the weight of this coordinate coord_map, val:6 y=x
							for( int y = 0; y < 12; y++ )
							{
								theIrregularFormulas[nfi].r_w[x] += coord_map[y*ncoords_base+x] * ceff_map[y];
								theIrregularFormulas[nfi].r_u[x] += coord_map[y*ncoords_base+x] * ceff_map_du[y] * scale_f_u;
								theIrregularFormulas[nfi].r_v[x] += coord_map[y*ncoords_base+x] * ceff_map_dv[y] * scale_f_v;
								
								theIrregularFormulas[nfi].r_uu[x] += coord_map[y*ncoords_base+x] * ceff_map_duu[y] * scale_f_u * scale_f_u;
								theIrregularFormulas[nfi].r_uv[x] += coord_map[y*ncoords_base+x] * ceff_map_duv[y] * scale_f_u * scale_f_v;
								theIrregularFormulas[nfi].r_vv[x] += coord_map[y*ncoords_base+x] * ceff_map_dvv[y] * scale_f_v * scale_f_v;
							}
	
							if( fabs(theIrregularFormulas[nfi].r_w[x]) > 1e10 )
							{
								printf("Err.\n");
								exit(1);
							}
						}
						
						theIrregularFormulas[nfi].c0 = theVertices[i].c0;
						theIrregularFormulas[nfi].vertex = i;
						theIrregularFormulas[nfi].edge = e;
		
						theIrregularFormulas[nfi].orig_u = o_fu; 
						theIrregularFormulas[nfi].orig_v = o_fv;
						theIrregularFormulas[nfi].tri = tri;
						theTriangles[tri].f = nf_faces+ (nfi / nf_irr_pts);

						nfi++;
					}
				}
			} 
			else
			{
				for( int g_q_p = 0; g_q_p < n_g_q_p; g_q_p++ )
				{
	
					double fu = pts[g_q_p][0];
					double fv = pts[g_q_p][1];
			
					if( fabs(fu) < 1e-8 && fabs(fv) < 1e-8 )
					{
						fu = 1e-8;
						fv = 1e-8;
					}
	
					double fw = 1-fu-fv;
	
					if( nfSpace == nf )
					{
						printf("Logical error calculating faces.\n");
						exit(1);
						nfSpace *= 2;
						theFormulas = (formula *)realloc( theFormulas, sizeof(formula) * nfSpace );
					}
	
					theFormulas[nf].ncoor = ncoords_base;
					theFormulas[nf].cp = (int *)malloc( sizeof(int) * ncoords_base ); 
	
					theFormulas[nf].tri = tri;
					theTriangles[tri].f = nf / nf_g_q_p;
					theFormulas[nf].r_pbc = (double *)malloc( sizeof(double) * 3 * ncoords_base );
					theFormulas[nf].r_w = (double *)malloc( sizeof(double) * ncoords_base ); 
					theFormulas[nf].r_u = (double *)malloc( sizeof(double) * ncoords_base ); 
					theFormulas[nf].r_v = (double *)malloc( sizeof(double) * ncoords_base ); 
					theFormulas[nf].r_uu = (double *)malloc( sizeof(double) * ncoords_base ); 
					theFormulas[nf].r_uv = (double *)malloc( sizeof(double) * ncoords_base ); 
					theFormulas[nf].r_vv = (double *)malloc( sizeof(double) * ncoords_base ); 
					theFormulas[nf].weight = weights[g_q_p]; 
					memcpy( theFormulas[nf].cp, point_list, sizeof(int) * ncoords_base );
					memset( theFormulas[nf].r_pbc, 0, sizeof(double) * 3 * ncoords_base );
	
					for( int p = 0; p < ncoords_base; p++ )
					{
						int *cp = point_list;
	
						double *r1 = theVertices[cp[p]].r;
						double *r0 = theVertices[cp[0]].r;
		
						double dr[3] = { r1[0] - r0[0], r1[1] - r0[1], r1[2] - r0[2] };
						double add[3]={0,0,0};
	
						MinImage3D( dr, PBC_vec, add );
							
						if( add[0]*add[0]+add[1]*add[1]+add[2]*add[2] > 1e-9 ) disable_PBC = 0;
	
						theFormulas[nf].r_pbc[3*p+0] = add[0]*PBC_vec[0][0] + add[1] * PBC_vec[1][0]+ add[2] * PBC_vec[2][0];
						theFormulas[nf].r_pbc[3*p+1] = add[0]*PBC_vec[0][1] + add[1] * PBC_vec[1][1]+ add[2] * PBC_vec[2][1];
						theFormulas[nf].r_pbc[3*p+2] = add[0]*PBC_vec[0][2] + add[1] * PBC_vec[1][2]+ add[2] * PBC_vec[2][2];
					}
	
					double scale_f_u = 1.0;
					double scale_f_v = 1.0;
	
					double coord_map[12 * ncoords_base];
					memset( coord_map, 0, sizeof(double) * 12 * ncoords_base );

					for( int ii = 0; ii < 12; ii++ )
						coord_map[ii*12+ii] = 1;	
	
					

					// get the coefficients.
						
						double u = fu;
						double v = fv;
						double w = 1 - u - v;
	
						double u2 = u*u;
						double u3 = u*u*u;
						double u4 = u*u*u*u;
						
						double v2 = v*v;
						double v3 = v*v*v;
						double v4 = v*v*v*v;
				
						double w2 = w*w;
						double w3 = w*w*w;
						double w4 = w*w*w*w;
						
					
					// 8 : 0
					// 7 : 1
					// 4 : 2
					// 5 : 3
					// 9 : 4
					// 12 : 5
					// 11 : 6
					// 10 : 7
					// 6 : 8
					// 3 : 9
					// 1 : 10
					// 2 : 11

	

					double n1 = (1.0/12.0)*(u4+2*u3*v); 
					double n2 = (1.0/12.0)*(u4+2*u3*w); 
					double n3 = (1.0/12.0)*(u4+2*u3*w+6*u3*v+6*u2*v*w+12*u2*v2+6*u*v2*w+6*u*v3+2*v3*w+v4); 
					double n4 = (1.0/12.0)*(6*u4+24*u3*w+24*u2*w2+8*u*w3+w4+24*u3*v+60*u2*v*w+36*u*v*w2+6*v*w3+24*u2*v2+36*u*v2*w+12*v2*w2+8*u*v3+6*v3*w+v4); 
					double n5 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+2*u3*v+6*u2*v*w+6*u*v*w2+2*v*w3); 
					double n6 = (1.0/12.0)*(2*u*v3+v4); 
					double n7 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+8*u3*v+36*u2*v*w+36*u*v*w2+8*v*w3+24*u2*v2+60*u*v2*w+24*v2*w2+24*u*v3+24*v3*w+6*v4); 
					double n8 = (1.0/12.0)*(u4+8*u3*w+24*u2*w2+24*u*w3+6*w4+6*u3*v+36*u2*v*w+60*u*v*w2+24*v*w3+12*u2*v2+36*u*v2*w+24*v2*w2+6*u*v3+8*v3*w+v4); 
					double n9 = (1.0/12.0)*(2*u*w3+w4); 
					double n10 = (1.0/12.0)*(2*v3*w+v4); 
					double n11 = (1.0/12.0)*(2*u*w3+w4+6*u*v*w2+6*v*w3+6*u*v2*w+12*v2*w2+2*u*v3+6*v3*w+v4); 
					double n12 = (1.0/12.0)*(w4+2*v*w3);

					double n1_alt = Power(u,4)/12. + (Power(u,3)*v)/6.;		
					double n2_alt = Power(u,4)/12. + (Power(u,3)*w)/6.;
					double n3_alt = Power(u,4)/12. + (Power(u,3)*v)/2. + Power(u,2)*Power(v,2) + (u*Power(v,3))/2. + Power(v,4)/12. + (Power(u,3)*w)/6. + (Power(u,2)*v*w)/2. + (u*Power(v,2)*w)/2. + (Power(v,3)*w)/6.;
					double n4_alt = Power(u,4)/2. + 2*Power(u,3)*v + 2*Power(u,2)*Power(v,2) + (2*u*Power(v,3))/3. + Power(v,4)/12. + 2*Power(u,3)*w + 5*Power(u,2)*v*w + 3*u*Power(v,2)*w + (Power(v,3)*w)/2. + 2*Power(u,2)*Power(w,2) + 3*u*v*Power(w,2) + Power(v,2)*Power(w,2) + (2*u*Power(w,3))/3. + (v*Power(w,3))/2. + Power(w,4)/12.;
					double n5_alt = Power(u,4)/12. + (Power(u,3)*v)/6. + (Power(u,3)*w)/2. + (Power(u,2)*v*w)/2. + Power(u,2)*Power(w,2) + (u*v*Power(w,2))/2. + (u*Power(w,3))/2. + (v*Power(w,3))/6. + Power(w,4)/12.;
					double n6_alt = (u*Power(v,3))/6. + Power(v,4)/12.;
					double n7_alt = Power(u,4)/12. + (2*Power(u,3)*v)/3. + 2*Power(u,2)*Power(v,2) + 2*u*Power(v,3) + Power(v,4)/2. + (Power(u,3)*w)/2. + 3*Power(u,2)*v*w + 5*u*Power(v,2)*w + 2*Power(v,3)*w + Power(u,2)*Power(w,2) + 3*u*v*Power(w,2) + 2*Power(v,2)*Power(w,2) + (u*Power(w,3))/2. + (2*v*Power(w,3))/3. + Power(w,4)/12.;
					double n8_alt = Power(u,4)/12. + (Power(u,3)*v)/2. + Power(u,2)*Power(v,2) + (u*Power(v,3))/2. + Power(v,4)/12. + (2*Power(u,3)*w)/3. + 3*Power(u,2)*v*w + 3*u*Power(v,2)*w + (2*Power(v,3)*w)/3. + 2*Power(u,2)*Power(w,2) + 5*u*v*Power(w,2) + 2*Power(v,2)*Power(w,2) + 2*u*Power(w,3) + 2*v*Power(w,3) + Power(w,4)/2.;
					double n9_alt = (u*Power(w,3))/6. + Power(w,4)/12.;
					double n10_alt = Power(v,4)/12. + (Power(v,3)*w)/6.;
					double n11_alt =(u*Power(v,3))/6. + Power(v,4)/12. + (u*Power(v,2)*w)/2. + (Power(v,3)*w)/2. + (u*v*Power(w,2))/2. + Power(v,2)*Power(w,2) + (u*Power(w,3))/6. + (v*Power(w,3))/2. + Power(w,4)/12.;
					double n12_alt = (v*Power(w,3))/6. + Power(w,4)/12.;
	
					double du_1 = Power(u,3)/3. + (Power(u,2)*v)/2.;
					double du_2 = Power(u,2)/2. - Power(u,3)/3. - (Power(u,2)*v)/2.; 	
					double du_3 =Power(u,2)/2. - Power(u,3)/3. + u*v - (Power(u,2)*v)/2. + Power(v,2)/2. - Power(v,3)/6.;
					double du_4 = 0.3333333333333333 + u - Power(u,2) - Power(u,3)/3. + v/2. - u*v - (Power(u,2)*v)/2. - Power(v,2) + Power(v,3)/3.;	
					double du_5 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. - v/2. + Power(u,2)*v + Power(v,2)/2. - Power(v,3)/6.;
					double du_6 = Power(v,3)/6.;	
					double du_7 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. + v/2. - 2*u*v + Power(u,2)*v - Power(v,2)/2. - Power(v,3)/6.;	
					double du_8 = -2*u + 2*Power(u,2) - Power(u,3)/3. - v + 2*u*v - (Power(u,2)*v)/2. + Power(v,2) - Power(v,3)/6.;	
					double du_9 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. + v/2. - (Power(u,2)*v)/2. - Power(v,2)/2. + Power(v,3)/6.;	
					double du_10 = -Power(v,3)/6.;
					double du_11 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. - v/2. + u*v - (Power(u,2)*v)/2. + Power(v,3)/3.;	
					double du_12 = 	-0.3333333333333333 + u - Power(u,2) + Power(u,3)/3. + v/2. - u*v + (Power(u,2)*v)/2. - Power(v,3)/6.;

					double dv_1 = Power(u,3)/6.;
					double dv_2 = -Power(u,3)/6.;
					double dv_3 = Power(u,2)/2. - Power(u,3)/6. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
					double dv_4 = 0.16666666666666666 + u/2. - Power(u,2)/2. - Power(u,3)/6. - 2*u*v - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
					double dv_5 = -0.16666666666666666 - u/2. + Power(u,3)/3. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
					double dv_6 = (u*Power(v,2))/2. + Power(v,3)/3.;
					double dv_7 = 0.3333333333333333 + u/2. - Power(u,2) + Power(u,3)/3. + v - u*v - Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
					double dv_8 = -u + Power(u,2) - Power(u,3)/6. - 2*v + 2*u*v + 2*Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
					double dv_9 = -0.3333333333333333 + u/2. - Power(u,3)/6. + v - u*v - Power(v,2) + (u*Power(v,2))/2. + Power(v,3)/3.;
					double dv_10 = Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
					double dv_11 = 0.16666666666666666 - u/2. + Power(u,2)/2. - Power(u,3)/6. - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
					double dv_12 = -0.16666666666666666 + u/2. - Power(u,2)/2. + Power(u,3)/6. + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
	
					double d_uu_1 = Power(u,2) + u*v;
					double d_uu_2 = u - Power(u,2) - u*v;
					double d_uu_3 = u - Power(u,2) + v - u*v;
					double d_uu_4 = 1 - 2*u - Power(u,2) - v - u*v;
					double d_uu_5 = -2*u + 2*Power(u,2) + 2*u*v;
					double d_uu_6 = 0;
					double d_uu_7 = -2*u + 2*Power(u,2) - 2*v + 2*u*v;
					double d_uu_8 = -2 + 4*u - Power(u,2) + 2*v - u*v;
					double d_uu_9 = u - Power(u,2) - u*v;
					double d_uu_10 = 0;
					double d_uu_11 = u - Power(u,2) + v - u*v;
					double d_uu_12 = 1 - 2*u + Power(u,2) - v + u*v;
					
					double d_uv_1 = Power(u,2)/2.;
					double d_uv_2 = -Power(u,2)/2.;
					double d_uv_3 = u - Power(u,2)/2. + v - Power(v,2)/2.;
					double d_uv_4 = 0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
					double d_uv_5 = -0.5 + Power(u,2) + v - Power(v,2)/2.;
					double d_uv_6 = Power(v,2)/2.;
					double d_uv_7 = 0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
					double d_uv_8 = -1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
					double d_uv_9 = 0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
					double d_uv_10 = -Power(v,2)/2.;
					double d_uv_11 = -0.5 + u - Power(u,2)/2. + Power(v,2);
					double d_uv_12 = 0.5 - u + Power(u,2)/2. - Power(v,2)/2.;
					
					double d_vv_1 = 0;
					double d_vv_2 = 0;
					double d_vv_3 = u + v - u*v - Power(v,2);
					double d_vv_4 = -2*u - 2*v + 2*u*v + 2*Power(v,2);
					double d_vv_5 = u + v - u*v - Power(v,2);
					double d_vv_6 = u*v + Power(v,2);
					double d_vv_7 = 1 - u - 2*v - u*v - Power(v,2);
					double d_vv_8 = -2 + 2*u + 4*v - u*v - Power(v,2);
					double d_vv_9 = 1 - u - 2*v + u*v + Power(v,2);
					double d_vv_10 = v - u*v - Power(v,2);
					double d_vv_11 = -2*v + 2*u*v + 2*Power(v,2);
					double d_vv_12 = v - u*v - Power(v,2);
		
					double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
					double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
					double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
					
					double ceff_map_duu[12] = { d_uu_8, d_uu_7, d_uu_4, d_uu_5, d_uu_9, d_uu_12, d_uu_11, d_uu_10, d_uu_6, d_uu_3, d_uu_1, d_uu_2 };
					double ceff_map_duv[12] = { d_uv_8, d_uv_7, d_uv_4, d_uv_5, d_uv_9, d_uv_12, d_uv_11, d_uv_10, d_uv_6, d_uv_3, d_uv_1, d_uv_2 };
					double ceff_map_dvv[12] = { d_vv_8, d_vv_7, d_vv_4, d_vv_5, d_vv_9, d_vv_12, d_vv_11, d_vv_10, d_vv_6, d_vv_3, d_vv_1, d_vv_2 };


					for( int x = 0; x < ncoords_base; x++ )
					{
						theFormulas[nf].r_w[x] = 0;
						theFormulas[nf].r_u[x] = 0;
						theFormulas[nf].r_v[x] = 0;
						theFormulas[nf].r_uu[x] = 0;
						theFormulas[nf].r_uv[x] = 0;
						theFormulas[nf].r_vv[x] = 0;
		
						// ceff_map: the weight of this coordinate coord_map, val:6 y=x
						for( int y = 0; y < 12; y++ )
						{
							theFormulas[nf].r_w[x] += coord_map[y*ncoords_base+x] * ceff_map[y];
							theFormulas[nf].r_u[x] += coord_map[y*ncoords_base+x] * ceff_map_du[y] * scale_f_u;
							theFormulas[nf].r_v[x] += coord_map[y*ncoords_base+x] * ceff_map_dv[y] * scale_f_v;
							
							theFormulas[nf].r_uu[x] += coord_map[y*ncoords_base+x] * ceff_map_duu[y] * scale_f_u * scale_f_u;
							theFormulas[nf].r_uv[x] += coord_map[y*ncoords_base+x] * ceff_map_duv[y] * scale_f_u * scale_f_v;
							theFormulas[nf].r_vv[x] += coord_map[y*ncoords_base+x] * ceff_map_dvv[y] * scale_f_v * scale_f_v;
						}

						if( fabs(theFormulas[nf].r_w[x]) > 1e10 )
						{
							printf("Err.\n");
							exit(1);
						}
					}
	
					theFormulas[nf].c0 = theVertices[i].c0;
					theFormulas[nf].vertex = i;
					theFormulas[nf].edge = e;
					theFormulas[nf].orig_u = pts[g_q_p][0]; 
					theFormulas[nf].orig_v = pts[g_q_p][1];
	
					nf++;
				}
			}
			// EDGE FORMULAS.
			
			for( int ep = 0; ep < 2; ep++ )
			{
				double fu = 0.5;
				double fv = 0;

				if( ep == 1 )
				{
					fu = 0;
					fv = 0.5;
				}
		
				double fw = 1-fu-fv;

				edgeFormulas[nef].ncoor = ncoords_base;
				edgeFormulas[nef].cp = (int *)malloc( sizeof(int) * ncoords_base ); 

				edgeFormulas[nef].r_pbc = (double *)malloc( sizeof(double) * 3 * ncoords_base );
				edgeFormulas[nef].r_w = (double *)malloc( sizeof(double) * ncoords_base ); 
				edgeFormulas[nef].r_u = (double *)malloc( sizeof(double) * ncoords_base ); 
				edgeFormulas[nef].r_v = (double *)malloc( sizeof(double) * ncoords_base ); 
				edgeFormulas[nef].r_uu = (double *)malloc( sizeof(double) * ncoords_base ); 
				edgeFormulas[nef].r_uv = (double *)malloc( sizeof(double) * ncoords_base ); 
				edgeFormulas[nef].r_vv = (double *)malloc( sizeof(double) * ncoords_base ); 
				edgeFormulas[nef].weight = 1.0; 
				memcpy( edgeFormulas[nef].cp, point_list, sizeof(int) * ncoords_base );
				memset( edgeFormulas[nef].r_pbc, 0, sizeof(double) * 3 * ncoords_base );

				for( int p = 0; p < ncoords_base; p++ )
				{
					int *cp = point_list;

					double *r1 = theVertices[cp[p]].r;
					double *r0 = theVertices[cp[0]].r;
	
					double dr[3] = { r1[0] - r0[0], r1[1] - r0[1], r1[2] - r0[2] };
					double add[3]={0,0,0};

					MinImage3D( dr, PBC_vec, add );
							
					if( add[0]*add[0]+add[1]*add[1]+add[2]*add[2] > 1e-9 ) disable_PBC = 0;

					edgeFormulas[nef].r_pbc[3*p+0] = add[0]*PBC_vec[0][0] + add[1] * PBC_vec[1][0]+ add[2] * PBC_vec[2][0];
					edgeFormulas[nef].r_pbc[3*p+1] = add[0]*PBC_vec[0][1] + add[1] * PBC_vec[1][1]+ add[2] * PBC_vec[2][1];
					edgeFormulas[nef].r_pbc[3*p+2] = add[0]*PBC_vec[0][2] + add[1] * PBC_vec[1][2]+ add[2] * PBC_vec[2][2];
				}

				double scale_f_u = 1.0;
				double scale_f_v = 1.0;

				double coord_map[12 * ncoords_base];
				memset( coord_map, 0, sizeof(double) * 12 * ncoords_base );

				if( val != 6 )
				{

					double f_domain = -log10(fu+fv)/log10(2.0);
					if( fu+fv >= 1.0 -1e-9 )
						f_domain = 1;
					int domain = lround(ceil(f_domain));
					double pow2 = pow( 2.0, domain-1.0 );
					
					scale_f_u *= pow2;
					scale_f_v *= pow2;

					fu *= pow2;
					fv *= pow2;

					int pcycle[12];

					if( fu > 0.5 )
					{
						fv = 2 * fv;
						fu = 2.0 * fu-1.0;
					
						scale_f_u *= 2;
						scale_f_v *= 2;

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
					else if( fv > 0.5 )
					{
						fu = 2 * fu;
						fv = 2.0 * fv-1.0;

						scale_f_u *= 2;
						scale_f_v *= 2;
					
						
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
					else
					{
						fv = 1 - 2 * fv;
						fu = 1 - 2 * fu;
						
						scale_f_u *= -2;
						scale_f_v *= -2;
						
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
							A_prev[e*ncoords_base+j] += Aevec[val][i*ncoords_base+j] * pow( Aval[val][i], domain-1 ) * Aevec_R[val][i*ncoords_base+e];
					}
					
					

					double A_use[(ncoords_base+ncoords_extra)*ncoords_base];
					memset( A_use, 0, sizeof(double) * (ncoords_base+ncoords_extra)*ncoords_base );
					for( int f = 0; f < ncoords_base+ncoords_extra; f++ )
					{
						for( int i = 0; i < ncoords_base; i++ )
						{
							A_use[f*ncoords_base+i] = 0;
							
							for( int e = 0; e < ncoords_base; e++ )
								A_use[f*ncoords_base+i] += Atot[val][f*ncoords_base+e] * A_prev[e*ncoords_base+i];
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
				}
				else
				{
					for( int i = 0; i < 12; i++ )
						coord_map[i*12+i] = 1;	
				}

				

				// get the coefficients.
					
					double u = fu;
					double v = fv;
					double w = 1 - u - v;

					double u2 = u*u;
					double u3 = u*u*u;
					double u4 = u*u*u*u;
					
					double v2 = v*v;
					double v3 = v*v*v;
					double v4 = v*v*v*v;
			
					double w2 = w*w;
					double w3 = w*w*w;
					double w4 = w*w*w*w;
					
				
				// 8 : 0
				// 7 : 1
				// 4 : 2
				// 5 : 3
				// 9 : 4
				// 12 : 5
				// 11 : 6
				// 10 : 7
				// 6 : 8
				// 3 : 9
				// 1 : 10
				// 2 : 11



				double n1 = (1.0/12.0)*(u4+2*u3*v); 
				double n2 = (1.0/12.0)*(u4+2*u3*w); 
				double n3 = (1.0/12.0)*(u4+2*u3*w+6*u3*v+6*u2*v*w+12*u2*v2+6*u*v2*w+6*u*v3+2*v3*w+v4); 
				double n4 = (1.0/12.0)*(6*u4+24*u3*w+24*u2*w2+8*u*w3+w4+24*u3*v+60*u2*v*w+36*u*v*w2+6*v*w3+24*u2*v2+36*u*v2*w+12*v2*w2+8*u*v3+6*v3*w+v4); 
				double n5 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+2*u3*v+6*u2*v*w+6*u*v*w2+2*v*w3); 
				double n6 = (1.0/12.0)*(2*u*v3+v4); 
				double n7 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+8*u3*v+36*u2*v*w+36*u*v*w2+8*v*w3+24*u2*v2+60*u*v2*w+24*v2*w2+24*u*v3+24*v3*w+6*v4); 
				double n8 = (1.0/12.0)*(u4+8*u3*w+24*u2*w2+24*u*w3+6*w4+6*u3*v+36*u2*v*w+60*u*v*w2+24*v*w3+12*u2*v2+36*u*v2*w+24*v2*w2+6*u*v3+8*v3*w+v4); 
				double n9 = (1.0/12.0)*(2*u*w3+w4); 
				double n10 = (1.0/12.0)*(2*v3*w+v4); 
				double n11 = (1.0/12.0)*(2*u*w3+w4+6*u*v*w2+6*v*w3+6*u*v2*w+12*v2*w2+2*u*v3+6*v3*w+v4); 
				double n12 = (1.0/12.0)*(w4+2*v*w3);

				double n1_alt = Power(u,4)/12. + (Power(u,3)*v)/6.;		
				double n2_alt = Power(u,4)/12. + (Power(u,3)*w)/6.;
				double n3_alt = Power(u,4)/12. + (Power(u,3)*v)/2. + Power(u,2)*Power(v,2) + (u*Power(v,3))/2. + Power(v,4)/12. + (Power(u,3)*w)/6. + (Power(u,2)*v*w)/2. + (u*Power(v,2)*w)/2. + (Power(v,3)*w)/6.;
				double n4_alt = Power(u,4)/2. + 2*Power(u,3)*v + 2*Power(u,2)*Power(v,2) + (2*u*Power(v,3))/3. + Power(v,4)/12. + 2*Power(u,3)*w + 5*Power(u,2)*v*w + 3*u*Power(v,2)*w + (Power(v,3)*w)/2. + 2*Power(u,2)*Power(w,2) + 3*u*v*Power(w,2) + Power(v,2)*Power(w,2) + (2*u*Power(w,3))/3. + (v*Power(w,3))/2. + Power(w,4)/12.;
				double n5_alt = Power(u,4)/12. + (Power(u,3)*v)/6. + (Power(u,3)*w)/2. + (Power(u,2)*v*w)/2. + Power(u,2)*Power(w,2) + (u*v*Power(w,2))/2. + (u*Power(w,3))/2. + (v*Power(w,3))/6. + Power(w,4)/12.;
				double n6_alt = (u*Power(v,3))/6. + Power(v,4)/12.;
				double n7_alt = Power(u,4)/12. + (2*Power(u,3)*v)/3. + 2*Power(u,2)*Power(v,2) + 2*u*Power(v,3) + Power(v,4)/2. + (Power(u,3)*w)/2. + 3*Power(u,2)*v*w + 5*u*Power(v,2)*w + 2*Power(v,3)*w + Power(u,2)*Power(w,2) + 3*u*v*Power(w,2) + 2*Power(v,2)*Power(w,2) + (u*Power(w,3))/2. + (2*v*Power(w,3))/3. + Power(w,4)/12.;
				double n8_alt = Power(u,4)/12. + (Power(u,3)*v)/2. + Power(u,2)*Power(v,2) + (u*Power(v,3))/2. + Power(v,4)/12. + (2*Power(u,3)*w)/3. + 3*Power(u,2)*v*w + 3*u*Power(v,2)*w + (2*Power(v,3)*w)/3. + 2*Power(u,2)*Power(w,2) + 5*u*v*Power(w,2) + 2*Power(v,2)*Power(w,2) + 2*u*Power(w,3) + 2*v*Power(w,3) + Power(w,4)/2.;
				double n9_alt = (u*Power(w,3))/6. + Power(w,4)/12.;
				double n10_alt = Power(v,4)/12. + (Power(v,3)*w)/6.;
				double n11_alt =(u*Power(v,3))/6. + Power(v,4)/12. + (u*Power(v,2)*w)/2. + (Power(v,3)*w)/2. + (u*v*Power(w,2))/2. + Power(v,2)*Power(w,2) + (u*Power(w,3))/6. + (v*Power(w,3))/2. + Power(w,4)/12.;
				double n12_alt = (v*Power(w,3))/6. + Power(w,4)/12.;

				double du_1 = Power(u,3)/3. + (Power(u,2)*v)/2.;
				double du_2 = Power(u,2)/2. - Power(u,3)/3. - (Power(u,2)*v)/2.; 	
				double du_3 =Power(u,2)/2. - Power(u,3)/3. + u*v - (Power(u,2)*v)/2. + Power(v,2)/2. - Power(v,3)/6.;
				double du_4 = 0.3333333333333333 + u - Power(u,2) - Power(u,3)/3. + v/2. - u*v - (Power(u,2)*v)/2. - Power(v,2) + Power(v,3)/3.;	
				double du_5 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. - v/2. + Power(u,2)*v + Power(v,2)/2. - Power(v,3)/6.;
				double du_6 = Power(v,3)/6.;	
				double du_7 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. + v/2. - 2*u*v + Power(u,2)*v - Power(v,2)/2. - Power(v,3)/6.;	
				double du_8 = -2*u + 2*Power(u,2) - Power(u,3)/3. - v + 2*u*v - (Power(u,2)*v)/2. + Power(v,2) - Power(v,3)/6.;	
				double du_9 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. + v/2. - (Power(u,2)*v)/2. - Power(v,2)/2. + Power(v,3)/6.;	
				double du_10 = -Power(v,3)/6.;
				double du_11 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. - v/2. + u*v - (Power(u,2)*v)/2. + Power(v,3)/3.;	
				double du_12 = 	-0.3333333333333333 + u - Power(u,2) + Power(u,3)/3. + v/2. - u*v + (Power(u,2)*v)/2. - Power(v,3)/6.;

				double dv_1 = Power(u,3)/6.;
				double dv_2 = -Power(u,3)/6.;
				double dv_3 = Power(u,2)/2. - Power(u,3)/6. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
				double dv_4 = 0.16666666666666666 + u/2. - Power(u,2)/2. - Power(u,3)/6. - 2*u*v - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
				double dv_5 = -0.16666666666666666 - u/2. + Power(u,3)/3. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
				double dv_6 = (u*Power(v,2))/2. + Power(v,3)/3.;
				double dv_7 = 0.3333333333333333 + u/2. - Power(u,2) + Power(u,3)/3. + v - u*v - Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
				double dv_8 = -u + Power(u,2) - Power(u,3)/6. - 2*v + 2*u*v + 2*Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
				double dv_9 = -0.3333333333333333 + u/2. - Power(u,3)/6. + v - u*v - Power(v,2) + (u*Power(v,2))/2. + Power(v,3)/3.;
				double dv_10 = Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
				double dv_11 = 0.16666666666666666 - u/2. + Power(u,2)/2. - Power(u,3)/6. - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
				double dv_12 = -0.16666666666666666 + u/2. - Power(u,2)/2. + Power(u,3)/6. + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;

				double d_uu_1 = Power(u,2) + u*v;
				double d_uu_2 = u - Power(u,2) - u*v;
				double d_uu_3 = u - Power(u,2) + v - u*v;
				double d_uu_4 = 1 - 2*u - Power(u,2) - v - u*v;
				double d_uu_5 = -2*u + 2*Power(u,2) + 2*u*v;
				double d_uu_6 = 0;
				double d_uu_7 = -2*u + 2*Power(u,2) - 2*v + 2*u*v;
				double d_uu_8 = -2 + 4*u - Power(u,2) + 2*v - u*v;
				double d_uu_9 = u - Power(u,2) - u*v;
				double d_uu_10 = 0;
				double d_uu_11 = u - Power(u,2) + v - u*v;
				double d_uu_12 = 1 - 2*u + Power(u,2) - v + u*v;
				
				double d_uv_1 = Power(u,2)/2.;
				double d_uv_2 = -Power(u,2)/2.;
				double d_uv_3 = u - Power(u,2)/2. + v - Power(v,2)/2.;
				double d_uv_4 = 0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
				double d_uv_5 = -0.5 + Power(u,2) + v - Power(v,2)/2.;
				double d_uv_6 = Power(v,2)/2.;
				double d_uv_7 = 0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
				double d_uv_8 = -1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
				double d_uv_9 = 0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
				double d_uv_10 = -Power(v,2)/2.;
				double d_uv_11 = -0.5 + u - Power(u,2)/2. + Power(v,2);
				double d_uv_12 = 0.5 - u + Power(u,2)/2. - Power(v,2)/2.;
				
				double d_vv_1 = 0;
				double d_vv_2 = 0;
				double d_vv_3 = u + v - u*v - Power(v,2);
				double d_vv_4 = -2*u - 2*v + 2*u*v + 2*Power(v,2);
				double d_vv_5 = u + v - u*v - Power(v,2);
				double d_vv_6 = u*v + Power(v,2);
				double d_vv_7 = 1 - u - 2*v - u*v - Power(v,2);
				double d_vv_8 = -2 + 2*u + 4*v - u*v - Power(v,2);
				double d_vv_9 = 1 - u - 2*v + u*v + Power(v,2);
				double d_vv_10 = v - u*v - Power(v,2);
				double d_vv_11 = -2*v + 2*u*v + 2*Power(v,2);
				double d_vv_12 = v - u*v - Power(v,2);
	
				double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
				double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
				double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
				
				double ceff_map_duu[12] = { d_uu_8, d_uu_7, d_uu_4, d_uu_5, d_uu_9, d_uu_12, d_uu_11, d_uu_10, d_uu_6, d_uu_3, d_uu_1, d_uu_2 };
				double ceff_map_duv[12] = { d_uv_8, d_uv_7, d_uv_4, d_uv_5, d_uv_9, d_uv_12, d_uv_11, d_uv_10, d_uv_6, d_uv_3, d_uv_1, d_uv_2 };
				double ceff_map_dvv[12] = { d_vv_8, d_vv_7, d_vv_4, d_vv_5, d_vv_9, d_vv_12, d_vv_11, d_vv_10, d_vv_6, d_vv_3, d_vv_1, d_vv_2 };


				for( int x = 0; x < ncoords_base; x++ )
				{
					edgeFormulas[nef].r_w[x] = 0;
					edgeFormulas[nef].r_u[x] = 0;
					edgeFormulas[nef].r_v[x] = 0;
					edgeFormulas[nef].r_uu[x] = 0;
					edgeFormulas[nef].r_uv[x] = 0;
					edgeFormulas[nef].r_vv[x] = 0;
	
					// ceff_map: the weight of this coordinate coord_map, val:6 y=x
					for( int y = 0; y < 12; y++ )
					{
						edgeFormulas[nef].r_w[x] += coord_map[y*ncoords_base+x] * ceff_map[y];
						edgeFormulas[nef].r_u[x] += coord_map[y*ncoords_base+x] * ceff_map_du[y] * scale_f_u;
						edgeFormulas[nef].r_v[x] += coord_map[y*ncoords_base+x] * ceff_map_dv[y] * scale_f_v;
						
						edgeFormulas[nef].r_uu[x] += coord_map[y*ncoords_base+x] * ceff_map_duu[y] * scale_f_u * scale_f_u;
						edgeFormulas[nef].r_uv[x] += coord_map[y*ncoords_base+x] * ceff_map_duv[y] * scale_f_u * scale_f_v;
						edgeFormulas[nef].r_vv[x] += coord_map[y*ncoords_base+x] * ceff_map_dvv[y] * scale_f_v * scale_f_v;
					}

					if( fabs(edgeFormulas[nef].r_w[x]) > 1e10 )
					{
						printf("Err.\n");
						exit(1);
					}
				}

				edgeFormulas[nef].c0 = theVertices[i].c0;

				if( ep == 0 )
				{
					edgeFormulas[nef].orig_u = 0.5; 
					edgeFormulas[nef].orig_v = 0;
				}
				else
				{
					edgeFormulas[nef].orig_u = 0.0; 
					edgeFormulas[nef].orig_v = 0.5;
				}
				nef++;
			}
		}
	}
	
	for( int f = 0; f < nf; f++ )
	{
		double r[3] = {0,0,0};
		
		int np = theFormulas[f].ncoor;

		int *cp = theFormulas[f].cp;
		for( int p = 0; p < np; p++ )
		{
			r[0] += theFormulas[f].r_w[p] *(theVertices[cp[p]].r[0] + theFormulas[f].r_pbc[3*p+0]); 
			r[1] += theFormulas[f].r_w[p] *(theVertices[cp[p]].r[1] + theFormulas[f].r_pbc[3*p+1]); 
			r[2] += theFormulas[f].r_w[p] *(theVertices[cp[p]].r[2] + theFormulas[f].r_pbc[3*p+2]); 
		}

		theFormulas[f].undeformed_r[0] = r[0];
		theFormulas[f].undeformed_r[1] = r[1];
		theFormulas[f].undeformed_r[2] = r[2];
	}
			
		

	double alpha_x = 1.0;
	double alpha_y = 1.0;
	double alpha_z = 1.0;
	
	double sum_reg_g = 0;

	for( int f = 0; f < nf; f++ )
	{
		double r[3] = {0,0,0};
		double ru[3] = {0,0,0};
		double rv[3] = {0,0,0};
		
		int np = theFormulas[f].ncoor;

		int *cp = theFormulas[f].cp;
		for( int p = 0; p < np; p++ )
		{
			r[0] += theFormulas[f].r_w[p] * alpha_x*(theVertices[cp[p]].r[0] + theFormulas[f].r_pbc[3*p+0]); 
			r[1] += theFormulas[f].r_w[p] * alpha_y*(theVertices[cp[p]].r[1] + theFormulas[f].r_pbc[3*p+1]); 
			r[2] += theFormulas[f].r_w[p] * alpha_z*(theVertices[cp[p]].r[2] + theFormulas[f].r_pbc[3*p+2]); 
			
			ru[0] += theFormulas[f].r_u[p] * alpha_x*(theVertices[cp[p]].r[0] + theFormulas[f].r_pbc[3*p+0]); 
			ru[1] += theFormulas[f].r_u[p] * alpha_y*(theVertices[cp[p]].r[1] + theFormulas[f].r_pbc[3*p+1]); 
			ru[2] += theFormulas[f].r_u[p] * alpha_z*(theVertices[cp[p]].r[2] + theFormulas[f].r_pbc[3*p+2]); 
			
			rv[0] += theFormulas[f].r_v[p] * alpha_x*(theVertices[cp[p]].r[0] + theFormulas[f].r_pbc[3*p+0]); 
			rv[1] += theFormulas[f].r_v[p] * alpha_y*(theVertices[cp[p]].r[1] + theFormulas[f].r_pbc[3*p+1]); 
			rv[2] += theFormulas[f].r_v[p] * alpha_z*(theVertices[cp[p]].r[2] + theFormulas[f].r_pbc[3*p+2]); 
		}
	
		double ruru = ru[0] * ru[0] + ru[1] * ru[1] + ru[2]*ru[2];
		double rurv = ru[0] * rv[0] + ru[1] * rv[1] + ru[2]*rv[2];
		double rvrv = rv[0] * rv[0] + rv[1] * rv[1] + rv[2]*rv[2];
	
		double g = sqrt(ruru*rvrv-rurv*rurv);

		sum_reg_g += g;

		normalize( ru );
		normalize( rv );

		double del = 2;
		double nrm[3];
		cross( ru, rv, nrm );
		normalize(nrm);
	}
		
	std_metric = sum_reg_g / nf;

	constructIrregularKernels();
}

void surface::updateFaceInfoForRandom( void )
{
	if( cumulative_area ) free(cumulative_area);
	if( cumulative_face ) free(cumulative_face);

	nfaces = nf_faces  + nf_irr_faces;
	cumulative_area = (double *)malloc( sizeof(double) * nfaces );
	cumulative_face = (int *)malloc( sizeof(int ) * nfaces );
	
	int xf = 0;
	
	double prev = 0;
	for( int f = 0; f < nf_faces; f++, xf++ )
	{
		double a = 0;
	
		for( int p = 0; p < nf_g_q_p; p++ )
			a += 0.5 * theFormulas[f*nf_g_q_p+p].g0 * theFormulas[f*nf_g_q_p+p].weight;
		prev += a;
		cumulative_area[xf] = prev;
	}
	
	for( int f = 0; f < nf_irr_faces; f++, xf++ )
	{
		double a = 0;
	
		for( int p = 0; p < nf_irr_pts; p++ )
			a += theIrregularFormulas[f*nf_irr_pts+p].g0 * theIrregularFormulas[f*nf_irr_pts+p].weight;
		prev += a;
		cumulative_area[xf] = prev;
	}
	
	for( int f = 0; f < nfaces; f++ )
		cumulative_area[f] /= prev;
	
}

void surface::randomPointOnSurface( int *face, double *u, double *v )
{
	if( !cumulative_area )
	{
		nfaces = nf_faces  + nf_irr_faces;
		cumulative_area = (double *)malloc( sizeof(double) * nfaces );
		cumulative_face = (int *)malloc( sizeof(int ) * nfaces );
	
		int xf = 0;
		
		double prev = 0;
		for( int f = 0; f < nf_faces; f++, xf++ )
		{
			double a = 0;
	
			for( int p = 0; p < nf_g_q_p; p++ )
				a += 0.5 * theFormulas[f*nf_g_q_p+p].g0 * theFormulas[f*nf_g_q_p+p].weight;
			prev += a;
			cumulative_area[xf] = prev;
		}
		
		for( int f = 0; f < nf_irr_faces; f++, xf++ )
		{
			double a = 0;
	
			for( int p = 0; p < nf_irr_pts; p++ )
				a += theIrregularFormulas[f*nf_irr_pts+p].g0 * theIrregularFormulas[f*nf_irr_pts+p].weight;
			prev += a;
			cumulative_area[xf] = prev;
		}
	
		for( int f = 0; f < nfaces; f++ )
			cumulative_area[f] /= prev;
	}

	double t1 = rand();

	double rnv = (t1) / (double)RAND_MAX;

	int min = 0;
	int max = nfaces-1;
	
	int done = 0;
	int cntr = 0;
			
	double lu = rand() / (double)RAND_MAX;	
	double lv = rand() / (double)RAND_MAX;	

	while( lu + lv > 1.0 )
	{
		lu = rand() / (double)RAND_MAX;	
		lv = rand() / (double)RAND_MAX;	
	}

	*u = lu;
	*v = lv;

	*face = -1;

	while(!done)
	{
		if( rnv < cumulative_area[(max+min)/2] )
			max = (max+min)/2;
		else
			min = (max+min)/2;

		if( max == min || max== 1 + min )
		{
			if( max >= nf_faces ) // an irrational point.
			{
				int tf = max - nf_faces;
				int val = theIrregularFormulas[tf*nf_irr_pts].ncoor - 6;
				int domain = kernels[val]->domain( lu, lv );

				while( domain < 1 || domain-1 >= MAX_DOMAIN )
				{
					double lu = rand() / (double)RAND_MAX;	
					double lv = rand() / (double)RAND_MAX;	
	
					while( lu + lv > 1.0 )
					{
						lu = rand() / (double)RAND_MAX;	
						lv = rand() / (double)RAND_MAX;	
					}

					*u = lu;
					*v = lv;
					
					domain = kernels[val]->domain( lu, lv );
				}	
			}

			*face = max;
			return;
		}

		if( cntr > 100 )
		{
			printf("logic problem.\n");	
			exit(1);
		}
	}
}


double surface::totalArea( void )
{
	if( !theFormulas )
		generatePlan();

	double dudv = 0.5;
	double area = 0;

	for( int f = 0; f < nf_faces; f++ )
		area += theFormulas[f].g0 * theFormulas[f].weight * dudv;

	return area;	
}	




int surface::sealEdge( int i, int j, int k )
{
	printf("sealing %d %d %d\n", i, j, k );
	int sorted[3] = { i,j,k};
	sort3(sorted);

	theTriangles[nt].ids[0] = sorted[0];
	theTriangles[nt].ids[1] = sorted[1];
	theTriangles[nt].ids[2] = sorted[2];

	theVertices[i].faces[theVertices[i].nfaces] = nt;
	theVertices[j].faces[theVertices[j].nfaces] = nt;
	theVertices[k].faces[theVertices[k].nfaces] = nt;

	theVertices[i].nfaces++;
        theVertices[j].nfaces++;
        theVertices[k].nfaces++;

	// ij
	for( int ex = 0; ex < theVertices[i].valence; ex++ )
	{
		if( theVertices[i].edges[ex] == j )
		{
			int e = theVertices[i].edgeIndices[ex];
			int the_code = 0;
			for( int px = 0; px < 3; px++ )	
			{
				if( theEdges[e].faces[px] == -1 )
				{
					theEdges[e].faces[px] = nt;
					the_code = theEdges[e].code[px] = getcode( theEdges[e].vertices, sorted );
					if( theEdges[e].code[px] < 0 )
					{
						printf("failed to assign code.\n");
						exit(1);	
					}
					break;
				}
			}
	
			if( the_code == 0 ) // jk
				theTriangles[nt].edges[1] = e; 
			else if( the_code == 1 ) // ik
				theTriangles[nt].edges[2] = e; 
			else   // ij
				theTriangles[nt].edges[0] = e; 
		}
	}
	
	// jk
	for( int ex = 0; ex < theVertices[j].valence; ex++ )
	{
		if( theVertices[j].edges[ex] == k )
		{
			int e = theVertices[j].edgeIndices[ex];
			int the_code = 0;
			for( int px = 0; px < 3; px++ )	
			{
				if( theEdges[e].faces[px] == -1 )
				{
					theEdges[e].faces[px] = nt;
					the_code = theEdges[e].code[px] = getcode( theEdges[e].vertices, sorted );
					if( theEdges[e].code[px] < 0 )
					{
						printf("failed to assign code.\n");
						exit(1);	
					}
					break;
				}
			}
	
			if( the_code == 0 ) // jk
				theTriangles[nt].edges[1] = e; 
			else if( the_code == 1 ) // ik
				theTriangles[nt].edges[2] = e; 
			else   // ij
				theTriangles[nt].edges[0] = e; 
		}
	}
	// ik
	for( int ex = 0; ex < theVertices[i].valence; ex++ )
	{
		if( theVertices[i].edges[ex] == k )
		{
			int e = theVertices[i].edgeIndices[ex];
			int the_code = 0;
			for( int px = 0; px < 3; px++ )	
			{
				if( theEdges[e].faces[px] == -1 )
				{
					theEdges[e].faces[px] = nt;
					the_code = theEdges[e].code[px] = getcode( theEdges[e].vertices, sorted );
					if( theEdges[e].code[px] < 0 )
					{
						printf("failed to assign code.\n");
						exit(1);	
					}
					break;
				}
			}
	
			if( the_code == 0 ) // jk
				theTriangles[nt].edges[1] = e; 
			else if( the_code == 1 ) // ik
				theTriangles[nt].edges[2] = e; 
			else   // ij
				theTriangles[nt].edges[0] = e; 
		}
	}

	theTriangles[nt].permutation = 0;
	theTriangles[nt].ft_vector = 0;
	theTriangles[nt].sense = 0;

	nt++;		

	return 1;
}

void surface::validateCodes( void )
{
	for( int e = 0; e < nedges; e++ )
	{
		int i = theEdges[e].vertices[0];
		int j = theEdges[e].vertices[1];
		
		int *verts = theEdges[e].vertices;

		for( int px = 0; px < 3; px++ )
		{
			if( theEdges[e].faces[px] >= 0 )
			{
				int check_code = getcode( verts, theTriangles[theEdges[e].faces[px]].ids );
				if( check_code != theEdges[e].code[px] )
				{
					printf("BAD CODE.\n");
					exit(1);
				}
			}
		}
	}
}


int surface::createEdge(int i, int j, int k )
{
	if( nt == nts )
	{
		nts *= 2;
		theTriangles = (triangle *)realloc( theTriangles, sizeof(triangle) * nts );
	}

	if( nedges == nedgesSpace )
	{
		nedgesSpace *= 2;
		theEdges = (edge *)realloc( theEdges, sizeof(edge) * nedgesSpace );
	}

	if( theVertices[j].valence < 4 )
	{
		printf("Attempting to create edge %d %d %d but rejected valence 3.\n", i, j, k );
		// refuse to possibly enclose a vertex of valence three.
		return 0;
	}

	if( theVertices[i].valence >= (MAX_INV_VALENCE-1) || theVertices[k].valence >= (MAX_INV_VALENCE-1) )
	{
		printf("Attempting to create edge %d %d %d but rejected valence %d.\n", i, j, k, MAX_INV_VALENCE );
		return 0;
	}
	printf("Creating edge %d %d %d\n", i, j, k );

	// NOT DONE YET.

	int sorted[3] = { i,j,k};
	int precode = sort3(sorted);

	int codemap[6] = { 1, 0, 2, 0, 2, 1 };   

	int link[2] = {i,k};
	if( k < i )
	{
		link[0] = k;
		link[1] = i;
	}
	theEdges[nedges].vertices[0] = link[0];
	theEdges[nedges].vertices[1] = link[1];
	theEdges[nedges].faces[0] = nt;
	theEdges[nedges].faces[1] = -1;
	theEdges[nedges].faces[2] = -1;
	theEdges[nedges].code[0] = getcode( theEdges[nedges].vertices, sorted );


	int val = theVertices[i].valence;
	theVertices[i].edges[val] = k;
	theVertices[i].edgeIndices[val] = nedges;
	double dr[3] = { theVertices[k].r[0] - theVertices[i].r[0],
			 theVertices[k].r[1] - theVertices[i].r[1],
			 theVertices[k].r[2] - theVertices[i].r[2] };

	MinImage3D( dr, PBC_vec, theVertices[i].edge_PBC+3*(val) );
	theVertices[i].valence++;

	double dr2[3] = { theVertices[i].r[0] - theVertices[k].r[0],
			 theVertices[i].r[1] - theVertices[k].r[1],
			 theVertices[i].r[2] - theVertices[k].r[2] };
	val = theVertices[k].valence;
	theVertices[k].edges[val] = i;
	theVertices[k].edgeIndices[val] = nedges;
	MinImage3D( dr2, PBC_vec, theVertices[k].edge_PBC+3*(val) );
	theVertices[k].valence++;

	theTriangles[nt].ids[0] = sorted[0];
	theTriangles[nt].ids[1] = sorted[1];
	theTriangles[nt].ids[2] = sorted[2];

	if( codemap[precode] == 0 ) // jk
		theTriangles[nt].edges[1] = nedges; 
	else if( codemap[precode] == 1 ) // ik
		theTriangles[nt].edges[2] = nedges; 
	else
		theTriangles[nt].edges[0] = nedges; 

	theTriangles[nt].permutation = 0;
	theTriangles[nt].ft_vector = 0;
	theTriangles[nt].sense = 0;


	for( int px = 0; px < 3; px++ )
	{
		int p = theTriangles[nt].ids[px];

		theVertices[p].faces[theVertices[p].nfaces] = nt;
		theVertices[p].nfaces++;
	}

	// link the other edges into this triangle .. ij and jk

	for( int ex = 0; ex < theVertices[i].valence; ex++ )
	{
		if( theVertices[i].edges[ex] == j )
		{
			int e = theVertices[i].edgeIndices[ex];
			int the_code = 0;
			for( int px = 0; px < 3; px++ )	
			{
				if( theEdges[e].faces[px] == -1 )
				{
					theEdges[e].faces[px] = nt;
					the_code = theEdges[e].code[px] = getcode( theEdges[e].vertices, sorted );
					if( theEdges[e].code[px] < 0 )
					{
						printf("failed to assign code.\n");
						exit(1);	
					}
					break;
				}
			}
	
			if( the_code == 0 ) // jk
				theTriangles[nt].edges[1] = e; 
			else if( the_code == 1 ) // ik
				theTriangles[nt].edges[2] = e; 
			else   // ij
				theTriangles[nt].edges[0] = e; 
		}
	}
	
	for( int ex = 0; ex < theVertices[j].valence; ex++ )
	{
		if( theVertices[j].edges[ex] == k )
		{
			int e = theVertices[j].edgeIndices[ex];
			int the_code = 0;
			for( int px = 0; px < 3; px++ )	
			{
				if( theEdges[e].faces[px] == -1 )
				{
					theEdges[e].faces[px] = nt;
					the_code = theEdges[e].code[px] = getcode( theEdges[e].vertices, sorted );
					if( theEdges[e].code[px] < 0 )
					{
						printf("failed to assign code.\n");
						exit(1);	
					}
					break;
				}
			}
			
			if( the_code == 0 ) // jk
				theTriangles[nt].edges[1] = e; 
			else if( the_code == 1 ) // ik
				theTriangles[nt].edges[2] = e; 
			else   // ij
				theTriangles[nt].edges[0] = e; 
		}
	}


	nedges++;
	nt++;
	
	return 1;	
}


void surface::assignEdgeIndices( void )
{
	for( int i = 0; i <nv; i++ )
	{
		for( int ex = 0; ex < theVertices[i].valence; ex++ )
		{
			if( theVertices[i].edgeIndices[ex] >= 0 )
				theVertices[i].edgeIndices[ex] = -2;
		}
	}

	for( int e = 0; e < nedges; e++ )
	{
		int i = theEdges[e].vertices[0];
		int j = theEdges[e].vertices[1];
		int val = theVertices[i].valence;
		for( int ex = 0; ex < val; ex++ )
		{
			if( theVertices[i].edges[ex] == j )
				theVertices[i].edgeIndices[ex] = e;
		}
		val = theVertices[j].valence;
		for( int ex = 0; ex < val; ex++ )
		{
			if( theVertices[j].edges[ex] == i )
				theVertices[j].edgeIndices[ex] = e;
		}
	}

	for( int i = 0; i < nv; i++ )
	{
		for( int ex = 0; ex < theVertices[i].valence; ex++ )
		{
			int j = theVertices[i].edges[ex];

			int e = theVertices[i].edgeIndices[ex];

			if( e == -1 ) continue;

			if( e == -2 )
			{
				printf("LACKS EDGE!!\n");
				exit(1);
			}

			if( !(i == theEdges[e].vertices[0] && j == theEdges[e].vertices[1]) && !(i == theEdges[e].vertices[1] && j == theEdges[e].vertices[0]) )	
			{
				printf("incomplete assignment.\n");
				exit(1);
			}
		}
	}
}

void surface::constructTriangles( void )
{
	nedges = 0;
	nt = 0;

	nts = 10;
	if( theTriangles )
	{
		for( int t = 0; t < nt; t++ )
		{
			if( theTriangles->npSpace > 0 )
				free(theTriangles->plist );
		}
		free(theTriangles);
	}
	if( theEdges )
		free(theEdges);
	theTriangles = (triangle *)malloc( sizeof(triangle) * nts );

	nedgesSpace = 0;

	for( int i = 0; i < nv; i++ )
		nedgesSpace += theVertices[i].valence;
	nedgesSpace /= 2;
	theEdges = (edge *)malloc( sizeof(edge) * nedgesSpace );

	for( int i = 0; i < nv; i++ )
	{
		theVertices[i].nfaces = 0;
		for( int e = 0; e < theVertices[i].valence; e++ )
			theVertices[i].edgeIndices[e] = -1;
	}
	for( int i = 0; i < nv; i++ )
	{

		for( int vi = 0; vi < theVertices[i].valence; vi++ )
		{
			int j = theVertices[i].edges[vi];

			if( j < i ) continue;

			for( int vj = 0; vj < theVertices[j].valence; vj++ )
			{
				int k = theVertices[j].edges[vj];

				if( k == i ) continue;
			
				if( k < j ) continue;
		
				int gotit = 0;

				for( int vk = 0; vk < theVertices[k].valence; vk++ )
				{
					int l = theVertices[k].edges[vk];

					if( l == i ) 
						gotit = 1;
				}

#ifdef DECLINE_SHARE // I take care of this below where I look for edges with three faces.
				int bad = 0;
		
				if( gotit )
				{
					// make sure they do not all share a vertex.

					for( int x = 0; x < theVertices[i].valence; x++ )
					{
						int l = theVertices[i].edges[x];
						
						for( int y = 0; y < theVertices[j].valence; y++ )
						{
							if( theVertices[j].edges[y] != l ) continue;
						
							for( int z = 0; z < theVertices[k].valence; z++ )
							{
								if( theVertices[k].edges[z] != l ) continue;
				
								printf("i: %d j: %d k: %d bad\n", i, j, k );						
								bad = 1;
							}
						}
					}
				}	
		
				if( bad ) continue;
#endif
				if( gotit )
				{
					if( nt == nts )
					{
						nts *= 2;
						theTriangles = (triangle *)realloc( theTriangles, nts * sizeof(triangle) );
					}
	
					theTriangles[nt].ids[0] = i;
					theTriangles[nt].ids[1] = j;
					theTriangles[nt].ids[2] = k;

					if( 
							(theVertices[i].nfaces >= 2*MAX_VALENCE ) ||
							(theVertices[j].nfaces >= 2*MAX_VALENCE ) ||
							(theVertices[k].nfaces >= 2*MAX_VALENCE )  )
					{
							printf("Too many faces.\n");
							exit(1);
					}

					theVertices[i].faces[theVertices[i].nfaces++] = nt;
					theVertices[j].faces[theVertices[j].nfaces++] = nt;
					theVertices[k].faces[theVertices[k].nfaces++] = nt;
	
					int gotit[3] = { 0,0,0};
					for( int e = 0; e < nedges; e++ )
					{
						if( theEdges[e].vertices[0] == i && theEdges[e].vertices[1] == j )
							gotit[0] = 1+e; 
						if( theEdges[e].vertices[0] == i && theEdges[e].vertices[1] == k )
							gotit[1] = 1+e; 
						if( theEdges[e].vertices[0] == j && theEdges[e].vertices[1] == k )
							gotit[2] = 1+e; 
					}

					if( !gotit[0]  )
					{
						theEdges[nedges].vertices[0] = i;
						theEdges[nedges].vertices[1] = j;
						theEdges[nedges].faces[0] = nt;
						theEdges[nedges].faces[1] = -1;
						theEdges[nedges].faces[2] = -1;
						theEdges[nedges].code[0] = 2;
						theTriangles[nt].edges[0] = nedges;
						nedges++;
					}
					else
					{
						if( theEdges[gotit[0]-1].faces[1]  == -1 )
						{
							theEdges[gotit[0]-1].faces[1] = nt;
							theEdges[gotit[0]-1].code[1] = 2;
						}
						else
						{
							theEdges[gotit[0]-1].faces[2] = nt;
							theEdges[gotit[0]-1].code[2] = 2;
						}
						theTriangles[nt].edges[0] = gotit[0]-1;
					}
					if( !gotit[1]  )
					{
						theEdges[nedges].vertices[0] = i;
						theEdges[nedges].vertices[1] = k;
						theEdges[nedges].faces[0] = nt;
						theEdges[nedges].faces[1] = -1;
						theEdges[nedges].faces[2] = -1;
						theEdges[nedges].code[0] = 1;
						theTriangles[nt].edges[2] = nedges;
						nedges++;
					}
					else
					{
						if( theEdges[gotit[1]-1].faces[1]  == -1 )
						{
							theEdges[gotit[1]-1].faces[1] = nt;
							theEdges[gotit[1]-1].code[1] = 1;
						}
						else
						{
							theEdges[gotit[1]-1].faces[2] = nt;
							theEdges[gotit[1]-1].code[2] = 1;
						}
						theTriangles[nt].edges[2] = gotit[1]-1;
					}
					if( !gotit[2]  )
					{
						theEdges[nedges].vertices[0] = j;
						theEdges[nedges].vertices[1] = k;
						theEdges[nedges].faces[0] = nt;
						theEdges[nedges].faces[1] = -1;
						theEdges[nedges].faces[2] = -1;
						theEdges[nedges].code[0] = 0;
						theTriangles[nt].edges[1] = nedges;
						nedges++;
					}
					else
					{
						if( theEdges[gotit[2]-1].faces[1]  == -1 )
						{
							theEdges[gotit[2]-1].faces[1] = nt;
							theEdges[gotit[2]-1].code[1] = 0;
						}
						else
						{
							theEdges[gotit[2]-1].faces[2] = nt;
							theEdges[gotit[2]-1].code[2] = 0;
						}
						theTriangles[nt].edges[1] = gotit[2]-1;
					}
					theTriangles[nt].np = 0;
					theTriangles[nt].npSpace = 0;
					theTriangles[nt].plist = NULL;
					theTriangles[nt].pc0 = NULL;
					theTriangles[nt].pa = NULL;
					nt++;
				}
			}
		}
	}
	
	assignEdgeIndices();
	
	for( int i = 0; i < nv; i++ )
	{
		for( int ex = 0; ex < theVertices[i].valence; ex++ )
		{
			int j = theVertices[i].edges[ex];

			if( theVertices[i].edgeIndices[ex] == -1 )
			{
				theEdges[nedges].vertices[0] = i;
				theEdges[nedges].vertices[1] = j;

				theEdges[nedges].faces[0] = -1;
				theEdges[nedges].faces[1] = -1;
				theEdges[nedges].faces[2] = -1;

				for( int exj = 0; exj < theVertices[j].valence; exj++ )
				{
					if( theVertices[j].edges[exj] == i )
						theVertices[j].edgeIndices[exj] = nedges;
				}
				theVertices[i].edgeIndices[ex] = nedges;
				nedges++;
			}
		}
	}
	
	for( int t = 0; t < nt; t++ )
	{
		int i = theTriangles[t].ids[0];
		int j = theTriangles[t].ids[1];
		int k = theTriangles[t].ids[2];

		double r1[3] = { theVertices[i].r[0], theVertices[i].r[1], theVertices[i].r[2] };
		double r2[3] = { theVertices[j].r[0], theVertices[j].r[1], theVertices[j].r[2] };
		double r3[3] = { theVertices[k].r[0], theVertices[k].r[1], theVertices[k].r[2] };

		double dr1[3] = { r2[0] - r1[0], r2[1]-r1[1], r2[2] - r1[2] };
		double dr2[3] = { r3[0] - r1[0], r3[1]-r1[1], r3[2] - r1[2] };

		double put[3];

		MinImage3D( dr1, PBC_vec, theTriangles[t].pbc1 );
		MinImage3D( dr2, PBC_vec, theTriangles[t].pbc2 );
		double o[3] = {0,0,0};

//		dr1[2] = 0;
//		dr2[2] = 0;

		theTriangles[t].A0 = triangle_area( o, dr1, dr2 );
//		printf("triangle %d A0: %lf\n", t, theTriangles[t].A0  );
	}
	
	if( nt > 0)
	{
		for( int t = 0; t < nt; t++ )	
		{

		int i = theTriangles[t].ids[0];
		int j = theTriangles[t].ids[1];
		int k = theTriangles[t].ids[2];
		
		double dr1[3] = { 
				theVertices[i].r[0] - theVertices[j].r[0],
				theVertices[i].r[1] - theVertices[j].r[1],
				theVertices[i].r[2] - theVertices[j].r[2] };
		double dr2[3] = { 
				theVertices[k].r[0] - theVertices[j].r[0],
				theVertices[k].r[1] - theVertices[j].r[1],
				theVertices[k].r[2] - theVertices[j].r[2] };

		double put[3];
		MinImage3D( dr1, PBC_vec, put );
		MinImage3D( dr2, PBC_vec, put );

		double fcom[3] = { 
			(3*theVertices[j].r[0] + dr1[0] + dr2[0] )/3,
			(3*theVertices[j].r[1] + dr1[1] + dr2[1] )/3,
			(3*theVertices[j].r[2] + dr1[2] + dr2[2] )/3 };

		normalize(fcom);
		double cp[3];

		cross( dr1, dr2, theTriangles[t].nrm);
		normalize(theTriangles[t].nrm);
		if( !(theTriangles[t].nrm[0] < 0 || theTriangles[t].nrm[0] > -2) )
		{
			printf("yikes.\n");
			exit(1);
		}
		theTriangles[t].dp  = fcom[0] * theTriangles[t].nrm[0];	
		theTriangles[t].dp += fcom[1] * theTriangles[t].nrm[1];	
		theTriangles[t].dp += fcom[2] * theTriangles[t].nrm[2];	

		theTriangles[t].sense = 0;
		}
	}
	else
	{
		printf("Construct triangles error.\n");
		exit(1);	
	}
	

	for( int t = 0; t < nt; t++ )
		theTriangles[t].draw = 1;

	for( int t = 0; t < nt; t++ )
	{
		theTriangles[t].composition.innerLeaflet = NULL;
		theTriangles[t].composition.outerLeaflet = NULL;
	}
}
	

void surface::orientSurface( void )
{

	constructTriangles();
//	writeXYZSurface("init.xyz","init.psf", this );
	

	int done = 0;
	int orientable_done = 0;
	int niter = 0;
	while( !orientable_done )
	{
		constructTriangles();

		for( int t = 0; t < nt; t++ )
			theTriangles[t].sense = 0;

		orientable_done = 1;
		theTriangles[0].sense = 1;
	
		// set the sense of each triangle that we can.

		done = 0;
		while( !done )
		{
			done = 1;
			for( int e = 0; e < nedges; e++ )
			{
				int f1 = theEdges[e].faces[0];
				int f2 = theEdges[e].faces[1];
		
				if( f1 < 0 || f2 < 0 ) continue;
	
				if( theTriangles[f1].sense && theTriangles[f2].sense )
				{
					int target_sense = theTriangles[f1].sense * compare_sense( theEdges[e].code[0], theEdges[e].code[1] );
					if( theTriangles[f2].sense != target_sense )
					{
						writeXYZandPSFPeriodic("bad");
						printf("Surface is not orientable.\n");
						exit(1);
					}
				}

				if( theTriangles[f1].sense && !theTriangles[f2].sense )
				{
					theTriangles[f2].sense = theTriangles[f1].sense * compare_sense( theEdges[e].code[0], theEdges[e].code[1] );
					done = 0;
				}
				if( !theTriangles[f1].sense && theTriangles[f2].sense )
				{
					theTriangles[f1].sense = theTriangles[f2].sense * compare_sense( theEdges[e].code[0], theEdges[e].code[1] );
					done = 0;
				}
			}
		}
	
		int orientation_problem[nv];
		memset( orientation_problem, 0, sizeof(int) * nv );

		for( int i = 0; i < nv; i++ )
		{

			// we can tolerate two pockets of oriented points before there is a problem.

			int val = theVertices[i].valence;
			
			int *sorter = (int *)malloc( sizeof(int) * val );
			for( int p = 0; p < val; p++ )
				sorter[p] = -1;
			int nsorter = 0;
	
			for( int pocket = 0; pocket < 2; pocket++ )
			{
				int *local_sorter = (int *)malloc( sizeof(int) * val );
				for( int p = 0; p < val; p++ )
					local_sorter[p] = -1;
				int cur_spot = 1;
				int nrun = 0;

				for( int direction = 1; direction >= -1; direction -= 2 )
				{
					int to_place = 1;

					if( direction == 1)
					{		
						// find a place to start that hasn't been taken yet.

						int got_spot = 0;

						for( int j = 0; j < val; j++ )
						{
							int use = 1;

							for( int i = 0; i < nsorter; i++ )
							{
								if( sorter[i] == j )
									use = 0;	 
							}

							if( use )
							{
								got_spot = 1;
								local_sorter[0] = j;
								nrun = 1;
								break;
							}
						}

						if( !got_spot )
							break;

					}
					else
					{
						to_place = val-1;

						if( local_sorter[to_place] != -1 )
							break;
					}
		
					int last_vert = theVertices[i].edges[local_sorter[0]];
			
					int orientation_failure = 0;
		
					while( local_sorter[to_place] == -1 )
					{
						int x_cur_spot = to_place;
		
						for( int v = 0; v < theVertices[i].nfaces; v++ )
						{
							int tri = theVertices[i].faces[v];
	
							if( theTriangles[tri].sense == 0 )
								continue;
			
							int points[3] = { theTriangles[tri].ids[0], theTriangles[tri].ids[1], theTriangles[tri].ids[2] };
							
							if( points[0] == last_vert || points[1] == last_vert || points[2] == last_vert )
							{
								int next_pt = points[0];
								if( next_pt == i || next_pt == last_vert )
									next_pt = points[1];
								if( next_pt == i || next_pt == last_vert )
									next_pt = points[2];
				
								int type = -1;
								
								if( points[0] == i && points[1] == last_vert )
									type = 1;		
								if( points[1] == i && points[2] == last_vert )
									type = 1;		
								if( points[2] == i && points[0] == last_vert )
									type = 1;	
			
								if( type * theTriangles[tri].sense * direction < 0 )
								{
									for( int xv = 0; xv < val; xv++ )
									{
										if( theVertices[i].edges[xv] == next_pt )
											local_sorter[to_place] = xv;
									}
									last_vert = next_pt;

									if( direction > 0 )
										to_place++;
									else
										to_place--;
							
									if( to_place < 0 ) to_place += val;
									if( to_place >= val ) to_place -= val;
									nrun++;
									break;
								}	
							}

						}
						
						if( to_place == x_cur_spot )
							break;
					}
				}
				// we have a "run" of points.

				if( nrun + nsorter > val )
				{
					//writeXYZSurface("bad.xyz","bad.psf", this);
					writeXYZandPSFPeriodic("bad");
					printf("BAD ERROR.\n");
					exit(1);
				}

				int endp = val-1;
				if( nrun == val || local_sorter[endp] == -1 ) // take the whole thing
				{
					if( nrun == val && nsorter != 0 )
					{
						printf("here.\n");
					}
					memcpy( sorter+nsorter, local_sorter, sizeof(int) * nrun );
					nsorter += nrun;	
				}
				else // take a subset, account for loop-back.
				{
					int nforward = 0;
					int nrev     = 0;
					while( local_sorter[endp] != -1 && endp >= 0)
					{	
						nrev++;
						endp--;
					}		
					int firstp = 0;	
					while( local_sorter[firstp] != -1 && firstp < val)
					{	
						nforward++;
						firstp++;
					}	
					if( nrev > 0 )
					{
						memcpy( sorter+nsorter, local_sorter+val-nrev, sizeof(int)*nrev );
						nsorter+=nrev;
					}
					if( nforward > 0 )
					{
						memcpy( sorter+nsorter, local_sorter, sizeof(int)*nforward );
						nsorter+=nforward;
					}
				}

				free(local_sorter);

			}
			
			
			if( nsorter == val )
			{
				int oe[MAX_VALENCE];
				memcpy( oe, theVertices[i].edges, sizeof(int) * theVertices[i].valence );
				double opbc[3*MAX_VALENCE];
				memcpy( opbc, theVertices[i].edge_PBC, sizeof(double) * theVertices[i].valence*3 );
		
				for( int x = 0; x < theVertices[i].valence; x++ )
				{
					theVertices[i].edges[x] = oe[sorter[x]];
					theVertices[i].edge_PBC[3*x+0] = opbc[sorter[x]*3+0];
					theVertices[i].edge_PBC[3*x+1] = opbc[sorter[x]*3+1];
					theVertices[i].edge_PBC[3*x+2] = opbc[sorter[x]*3+2];
				}
			}
			else
			{
				orientation_problem[i] = 1;
				printf("Index %d could not be oriented!\n", i );
			}

			free(sorter);
		}
		

		setEdgeRev();	

		assignEdgeIndices();
		// analyze for broken edges.

		niter++;

		if( niter > 2000 )
		{
			writeXYZSurface("failed.xyz","failed.psf", this );
			
			exit(1);
		}
	}
		
		
	setEdgeRev();	
	assignEdgeIndices();
	
}


void surface::nearestPointOnSurface( double *r, int *vertex  )
{
	int bestv = 0;
	double bestr2 = 1e10;

	for( int v = 0; v < nv; v++ )
	{
		double dr[3] = { r[0] - theVertices[v].r[0],
				 r[1] - theVertices[v].r[1],
				 r[2] - theVertices[v].r[2] };
		double put[3];
		MinImage3D( dr, PBC_vec, put );
		double r = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];	

		if( r < bestr2 )
		{
			bestr2 = r;
			bestv = v;
		}
	}	

	*vertex = bestv;

}



void surface::flipOrientation( void )
{
	for( int i = 0; i < nv; i++ )
	{
		int val = theVertices[i].valence;

		int oe[val];
		double opbc[val*3];
		memcpy( oe, theVertices[i].edges, sizeof(int) * val );
		memcpy( opbc, theVertices[i].edge_PBC, sizeof(double) * 3 * val );

		for( int x = 0; x < theVertices[i].valence; x++ )
		{
			theVertices[i].edges[x] = oe[val-x-1];
			theVertices[i].edge_PBC[3*x+0] = opbc[(val-x-1)*3+0];
			theVertices[i].edge_PBC[3*x+1] = opbc[(val-x-1)*3+1];
			theVertices[i].edge_PBC[3*x+2] = opbc[(val-x-1)*3+2];
		}
	}

	setEdgeRev();

}



void surface::addParticleToFace( int f, int pid, double p_c0, double p_area )
{
	int t = 0;

	if( f >= nf_faces )
		t = theIrregularFormulas[(f-nf_faces)*nf_irr_pts].tri;
	else
		t = theFormulas[f*nf_g_q_p].tri;

	if( theTriangles[t].npSpace == 0 )
	{
		theTriangles[t].npSpace = 1;
		theTriangles[t].plist = (int *)malloc( sizeof(int) * theTriangles[t].npSpace );
		theTriangles[t].pc0 = (double *)malloc( sizeof(double) * theTriangles[t].npSpace );
		theTriangles[t].pa = (double *)malloc( sizeof(double) * theTriangles[t].npSpace );
		theTriangles[t].np = 0;
	}
	else if( theTriangles[t].npSpace == theTriangles[t].np )
	{
		theTriangles[t].npSpace *= 2;
		theTriangles[t].plist = (int *)realloc( theTriangles[t].plist, sizeof(int) * theTriangles[t].npSpace );
		theTriangles[t].pc0 = (double *)realloc( theTriangles[t].pc0, sizeof(double) * theTriangles[t].npSpace );
		theTriangles[t].pa = (double *)realloc( theTriangles[t].pa, sizeof(double) * theTriangles[t].npSpace );
	}

	for( int tx = 0; tx < theTriangles[t].np; tx++ )
	{
		if( theTriangles[t].plist[tx] == pid )
		{
			printf("Duplicate point here!\n");
			exit(1);
		}
	}

	theTriangles[t].plist[theTriangles[t].np] = pid;
	theTriangles[t].pc0[theTriangles[t].np] = p_c0;
	theTriangles[t].pa[theTriangles[t].np] = p_area;
	theTriangles[t].np += 1; 

/* We are no longer averaging the spontaneous curvature into the mesh face.
   Now we compute the local curvature of each particle and subtract out the pure membrane energy.
   We do need the spontaneous curvature of the embedded particle.

	double a = 0;

	for( int p = 0; p < nf_g_q_p; p++ )
		a += theFormulas[f*nf_g_q_p+p].g0 * theFormulas[f*nf_g_q_p+p].weight * 0.5;

	for( int p = 0; p < nf_g_q_p; p++ )
		theFormulas[f*nf_g_q_p+p].c0 += c0_to_add / a;
*/
}

void surface::removeParticleFromFace( int f, int pid, double c0_to_add, double p_area )
{
	int t = 0;

	if( f >= nf_faces )
		t = theIrregularFormulas[(f-nf_faces)*nf_irr_pts].tri;
	else
		t = theFormulas[f*nf_g_q_p].tri;

	for( int px = 0; px < theTriangles[t].np; px++ )
	{
		if( theTriangles[t].plist[px] == pid )
		{
			theTriangles[t].pc0[px] = theTriangles[t].pc0[theTriangles[t].np-1];
			theTriangles[t].pa[px] = theTriangles[t].pa[theTriangles[t].np-1];
			theTriangles[t].plist[px] = theTriangles[t].plist[theTriangles[t].np-1];
			theTriangles[t].np--;
			break;
		}
	}

/* We are no longer averaging the spontaneous curvature into the mesh face.
   Now we compute the local curvature of each particle and subtract out the pure membrane energy.
   We do need the spontaneous curvature of the embedded particle.
	double a = 0;

	for( int p = 0; p < nf_g_q_p; p++ )
		a += theFormulas[f*nf_g_q_p+p].g0 * theFormulas[f*nf_g_q_p+p].weight * 0.5;

	for( int p = 0; p < nf_g_q_p; p++ )
		theFormulas[f*nf_g_q_p+p].c0 -= c0_to_add / a;
*/
}

void surface::evaluateRNRM( int f, double u, double v, double *rp, double *nrm, double *r )
{
	double alpha_x = r[3*nv+0];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];

	if( f >= nf_faces ) // an irrational face.
	{
		int frm = (f-nf_faces)*nf_irr_pts;

		int i = theIrregularFormulas[frm].vertex;
		int e = theIrregularFormulas[frm].edge;
		int val = theVertices[i].valence;

		irr_kernel *theKernel = kernels[val]; 
		int ncoords_base = val + 6;

		double fu = u;
		double fv = v;					

		int domain = theKernel->domain(fu,fv);
		if( domain-1 >= theKernel->ndomains )
		{
			printf("LOGICAL ERROR evaluating R,nrm on irregular vertex face.\n");
			exit(1);
		}	
		
		if( domain < 0 &&  fu < 1e-10 && fv < 1e-10  )
		{
			// special case.
			nrm[0] = nan("1");
			nrm[1] = nan("1");
			nrm[2] = nan("1"); // unknown for now;
			int v = i;

			int val = theVertices[i].valence;

			double lr[3] = { 0,0,0};
	
			double *ti = theVertices[i].r;
	
			lr[0] = 0.5 * theVertices[i].r[0]; 
			lr[1] = 0.5 * theVertices[i].r[1]; 
			lr[2] = 0.5 * theVertices[i].r[2]; 
	
			double w = 1.0 / (val * 2);
	
			for( int e = 0; e < val; e++ )
			{
				int j = theVertices[i].edges[e];
	
				double dr[3] = { theVertices[j].r[0] - ti[0], theVertices[j].r[1] - ti[1], theVertices[j].r[2] - ti[2] };
	
				double put[3];	

				if( !disable_PBC )
					MinImage3D( dr, PBC_vec, put );

				dr[0] += ti[0];
				dr[1] += ti[1];
				dr[2] += ti[2];
	
				lr[0] += dr[0]*w;
				lr[1] += dr[1]*w;
				lr[2] += dr[2]*w;
			}	
			rp[0] = lr[0] * alpha_x;	
			rp[1] = lr[1] * alpha_y;	
			rp[2] = lr[2] * alpha_z;	
			return;
		}
		double *theMap = theKernel->get_map( &fu, &fv );

		double u = fu;
		double v = fv;
		double w = 1 - u - v;


		if( u +v > 1.0 + THRESH || u < -THRESH || v < -THRESH )
		{
			printf("u: %.12le v: %.12le\n", u, v );
			printf("ERROR: outside domain.\n");
			exit(1);
		}

		double u2 = u*u;
		double u3 = u*u*u;
		double u4 = u*u*u*u;
		
		double v2 = v*v;
		double v3 = v*v*v;
		double v4 = v*v*v*v;
		
		double w2 = w*w;
		double w3 = w*w*w;
		double w4 = w*w*w*w;
			
		
		// 8 : 0
		// 7 : 1
		// 4 : 2
		// 5 : 3
		// 9 : 4
		// 12 : 5
		// 11 : 6
		// 10 : 7
		// 6 : 8
		// 3 : 9
		// 1 : 10
		// 2 : 11

	

		double n1 = (1.0/12.0)*(u4+2*u3*v); 
		double n2 = (1.0/12.0)*(u4+2*u3*w); 
		double n3 = (1.0/12.0)*(u4+2*u3*w+6*u3*v+6*u2*v*w+12*u2*v2+6*u*v2*w+6*u*v3+2*v3*w+v4); 
		double n4 = (1.0/12.0)*(6*u4+24*u3*w+24*u2*w2+8*u*w3+w4+24*u3*v+60*u2*v*w+36*u*v*w2+6*v*w3+24*u2*v2+36*u*v2*w+12*v2*w2+8*u*v3+6*v3*w+v4); 
		double n5 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+2*u3*v+6*u2*v*w+6*u*v*w2+2*v*w3); 
		double n6 = (1.0/12.0)*(2*u*v3+v4); 
		double n7 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+8*u3*v+36*u2*v*w+36*u*v*w2+8*v*w3+24*u2*v2+60*u*v2*w+24*v2*w2+24*u*v3+24*v3*w+6*v4); 
		double n8 = (1.0/12.0)*(u4+8*u3*w+24*u2*w2+24*u*w3+6*w4+6*u3*v+36*u2*v*w+60*u*v*w2+24*v*w3+12*u2*v2+36*u*v2*w+24*v2*w2+6*u*v3+8*v3*w+v4); 
		double n9 = (1.0/12.0)*(2*u*w3+w4); 
		double n10 = (1.0/12.0)*(2*v3*w+v4); 
		double n11 = (1.0/12.0)*(2*u*w3+w4+6*u*v*w2+6*v*w3+6*u*v2*w+12*v2*w2+2*u*v3+6*v3*w+v4); 
		double n12 = (1.0/12.0)*(w4+2*v*w3);

		double du_1 = Power(u,3)/3. + (Power(u,2)*v)/2.;
		double du_2 = Power(u,2)/2. - Power(u,3)/3. - (Power(u,2)*v)/2.; 	
		double du_3 =Power(u,2)/2. - Power(u,3)/3. + u*v - (Power(u,2)*v)/2. + Power(v,2)/2. - Power(v,3)/6.;
		double du_4 = 0.3333333333333333 + u - Power(u,2) - Power(u,3)/3. + v/2. - u*v - (Power(u,2)*v)/2. - Power(v,2) + Power(v,3)/3.;	
		double du_5 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. - v/2. + Power(u,2)*v + Power(v,2)/2. - Power(v,3)/6.;
		double du_6 = Power(v,3)/6.;	
		double du_7 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. + v/2. - 2*u*v + Power(u,2)*v - Power(v,2)/2. - Power(v,3)/6.;	
		double du_8 = -2*u + 2*Power(u,2) - Power(u,3)/3. - v + 2*u*v - (Power(u,2)*v)/2. + Power(v,2) - Power(v,3)/6.;	
		double du_9 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. + v/2. - (Power(u,2)*v)/2. - Power(v,2)/2. + Power(v,3)/6.;	
		double du_10 = -Power(v,3)/6.;
		double du_11 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. - v/2. + u*v - (Power(u,2)*v)/2. + Power(v,3)/3.;	
		double du_12 = 	-0.3333333333333333 + u - Power(u,2) + Power(u,3)/3. + v/2. - u*v + (Power(u,2)*v)/2. - Power(v,3)/6.;

		double dv_1 = Power(u,3)/6.;
		double dv_2 = -Power(u,3)/6.;
		double dv_3 = Power(u,2)/2. - Power(u,3)/6. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_4 = 0.16666666666666666 + u/2. - Power(u,2)/2. - Power(u,3)/6. - 2*u*v - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
		double dv_5 = -0.16666666666666666 - u/2. + Power(u,3)/3. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_6 = (u*Power(v,2))/2. + Power(v,3)/3.;
		double dv_7 = 0.3333333333333333 + u/2. - Power(u,2) + Power(u,3)/3. + v - u*v - Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_8 = -u + Power(u,2) - Power(u,3)/6. - 2*v + 2*u*v + 2*Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_9 = -0.3333333333333333 + u/2. - Power(u,3)/6. + v - u*v - Power(v,2) + (u*Power(v,2))/2. + Power(v,3)/3.;
		double dv_10 = Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_11 = 0.16666666666666666 - u/2. + Power(u,2)/2. - Power(u,3)/6. - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
		double dv_12 = -0.16666666666666666 + u/2. - Power(u,2)/2. + Power(u,3)/6. + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
	
		double d_uu_1 = Power(u,2) + u*v;
		double d_uu_2 = u - Power(u,2) - u*v;
		double d_uu_3 = u - Power(u,2) + v - u*v;
		double d_uu_4 = 1 - 2*u - Power(u,2) - v - u*v;
		double d_uu_5 = -2*u + 2*Power(u,2) + 2*u*v;
		double d_uu_6 = 0;
		double d_uu_7 = -2*u + 2*Power(u,2) - 2*v + 2*u*v;
		double d_uu_8 = -2 + 4*u - Power(u,2) + 2*v - u*v;
		double d_uu_9 = u - Power(u,2) - u*v;
		double d_uu_10 = 0;
		double d_uu_11 = u - Power(u,2) + v - u*v;
		double d_uu_12 = 1 - 2*u + Power(u,2) - v + u*v;
		
		double d_uv_1 = Power(u,2)/2.;
		double d_uv_2 = -Power(u,2)/2.;
		double d_uv_3 = u - Power(u,2)/2. + v - Power(v,2)/2.;
		double d_uv_4 = 0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
		double d_uv_5 = -0.5 + Power(u,2) + v - Power(v,2)/2.;
		double d_uv_6 = Power(v,2)/2.;
		double d_uv_7 = 0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
		double d_uv_8 = -1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
		double d_uv_9 = 0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
		double d_uv_10 = -Power(v,2)/2.;
		double d_uv_11 = -0.5 + u - Power(u,2)/2. + Power(v,2);
		double d_uv_12 = 0.5 - u + Power(u,2)/2. - Power(v,2)/2.;
		
		double d_vv_1 = 0;
		double d_vv_2 = 0;
		double d_vv_3 = u + v - u*v - Power(v,2);
		double d_vv_4 = -2*u - 2*v + 2*u*v + 2*Power(v,2);
		double d_vv_5 = u + v - u*v - Power(v,2);
		double d_vv_6 = u*v + Power(v,2);
		double d_vv_7 = 1 - u - 2*v - u*v - Power(v,2);
		double d_vv_8 = -2 + 2*u + 4*v - u*v - Power(v,2);
		double d_vv_9 = 1 - u - 2*v + u*v + Power(v,2);
		double d_vv_10 = v - u*v - Power(v,2);
		double d_vv_11 = -2*v + 2*u*v + 2*Power(v,2);
		double d_vv_12 = v - u*v - Power(v,2);
		
		double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
		double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
		double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
		
		double ceff_map_duu[12] = { d_uu_8, d_uu_7, d_uu_4, d_uu_5, d_uu_9, d_uu_12, d_uu_11, d_uu_10, d_uu_6, d_uu_3, d_uu_1, d_uu_2 };
		double ceff_map_duv[12] = { d_uv_8, d_uv_7, d_uv_4, d_uv_5, d_uv_9, d_uv_12, d_uv_11, d_uv_10, d_uv_6, d_uv_3, d_uv_1, d_uv_2 };
		double ceff_map_dvv[12] = { d_vv_8, d_vv_7, d_vv_4, d_vv_5, d_vv_9, d_vv_12, d_vv_11, d_vv_10, d_vv_6, d_vv_3, d_vv_1, d_vv_2 };
		
		rp[0] = 0;
		rp[1] = 0;
		rp[2] = 0;
		double lru[3] = { 0,0,0};
		double lrv[3] = { 0,0,0};
		for( int x = 0; x < ncoords_base; x++ )
		{
			int *cset = theVertices[i].irr_coord_set + e * ncoords_base;
			double *lr = r + cset[x]*3;

			for( int y = 0; y < 12; y++ )
			{
				rp[0] += (lr[0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map[y];
				rp[1] += (lr[1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map[y];
				rp[2] += (lr[2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map[y];
				
				lru[0] += (lr[0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_du[y];
				lru[1] += (lr[1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_du[y];
				lru[2] += (lr[2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_du[y];
				
				lrv[0] += (lr[0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_dv[y];
				lrv[1] += (lr[1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_dv[y];
				lrv[2] += (lr[2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_dv[y];
			}
		}
		
		rp[0] *= alpha_x;
		rp[1] *= alpha_y;
		rp[2] *= alpha_z;
		
		cross( lru, lrv, nrm );
		normalize(nrm);
	}
	else
	{
		int frm = f*nf_g_q_p;
	
		int *cp = theFormulas[frm].cp;
	
		double w = 1-u-v;
	
		double u2 = u*u;
		double u3 = u*u*u;
		double u4 = u*u*u*u;
		
		double v2 = v*v;
		double v3 = v*v*v;
		double v4 = v*v*v*v;
		
		double w2 = w*w;
		double w3 = w*w*w;
		double w4 = w*w*w*w;
			
		
		// 8 : 0
		// 7 : 1
		// 4 : 2
		// 5 : 3
		// 9 : 4
		// 12 : 5
		// 11 : 6
		// 10 : 7
		// 6 : 8
		// 3 : 9
		// 1 : 10
		// 2 : 11
	
		double n1 = (1.0/12.0)*(u4+2*u3*v); 
		double n2 = (1.0/12.0)*(u4+2*u3*w); 
		double n3 = (1.0/12.0)*(u4+2*u3*w+6*u3*v+6*u2*v*w+12*u2*v2+6*u*v2*w+6*u*v3+2*v3*w+v4); 
		double n4 = (1.0/12.0)*(6*u4+24*u3*w+24*u2*w2+8*u*w3+w4+24*u3*v+60*u2*v*w+36*u*v*w2+6*v*w3+24*u2*v2+36*u*v2*w+12*v2*w2+8*u*v3+6*v3*w+v4); 
		double n5 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+2*u3*v+6*u2*v*w+6*u*v*w2+2*v*w3); 
		double n6 = (1.0/12.0)*(2*u*v3+v4); 
		double n7 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+8*u3*v+36*u2*v*w+36*u*v*w2+8*v*w3+24*u2*v2+60*u*v2*w+24*v2*w2+24*u*v3+24*v3*w+6*v4); 
		double n8 = (1.0/12.0)*(u4+8*u3*w+24*u2*w2+24*u*w3+6*w4+6*u3*v+36*u2*v*w+60*u*v*w2+24*v*w3+12*u2*v2+36*u*v2*w+24*v2*w2+6*u*v3+8*v3*w+v4); 
		double n9 = (1.0/12.0)*(2*u*w3+w4); 
		double n10 = (1.0/12.0)*(2*v3*w+v4); 
		double n11 = (1.0/12.0)*(2*u*w3+w4+6*u*v*w2+6*v*w3+6*u*v2*w+12*v2*w2+2*u*v3+6*v3*w+v4); 
		double n12 = (1.0/12.0)*(w4+2*v*w3);
						
		double du_1 = Power(u,3)/3. + (Power(u,2)*v)/2.;
		double du_2 = Power(u,2)/2. - Power(u,3)/3. - (Power(u,2)*v)/2.; 	
		double du_3 =Power(u,2)/2. - Power(u,3)/3. + u*v - (Power(u,2)*v)/2. + Power(v,2)/2. - Power(v,3)/6.;
		double du_4 = 0.3333333333333333 + u - Power(u,2) - Power(u,3)/3. + v/2. - u*v - (Power(u,2)*v)/2. - Power(v,2) + Power(v,3)/3.;	
		double du_5 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. - v/2. + Power(u,2)*v + Power(v,2)/2. - Power(v,3)/6.;
		double du_6 = Power(v,3)/6.;	
		double du_7 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. + v/2. - 2*u*v + Power(u,2)*v - Power(v,2)/2. - Power(v,3)/6.;	
		double du_8 = -2*u + 2*Power(u,2) - Power(u,3)/3. - v + 2*u*v - (Power(u,2)*v)/2. + Power(v,2) - Power(v,3)/6.;	
		double du_9 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. + v/2. - (Power(u,2)*v)/2. - Power(v,2)/2. + Power(v,3)/6.;	
		double du_10 = -Power(v,3)/6.;
		double du_11 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. - v/2. + u*v - (Power(u,2)*v)/2. + Power(v,3)/3.;	
		double du_12 = 	-0.3333333333333333 + u - Power(u,2) + Power(u,3)/3. + v/2. - u*v + (Power(u,2)*v)/2. - Power(v,3)/6.;
	
		double dv_1 = Power(u,3)/6.;
		double dv_2 = -Power(u,3)/6.;
		double dv_3 = Power(u,2)/2. - Power(u,3)/6. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_4 = 0.16666666666666666 + u/2. - Power(u,2)/2. - Power(u,3)/6. - 2*u*v - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
		double dv_5 = -0.16666666666666666 - u/2. + Power(u,3)/3. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_6 = (u*Power(v,2))/2. + Power(v,3)/3.;
		double dv_7 = 0.3333333333333333 + u/2. - Power(u,2) + Power(u,3)/3. + v - u*v - Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_8 = -u + Power(u,2) - Power(u,3)/6. - 2*v + 2*u*v + 2*Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_9 = -0.3333333333333333 + u/2. - Power(u,3)/6. + v - u*v - Power(v,2) + (u*Power(v,2))/2. + Power(v,3)/3.;
		double dv_10 = Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_11 = 0.16666666666666666 - u/2. + Power(u,2)/2. - Power(u,3)/6. - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
		double dv_12 = -0.16666666666666666 + u/2. - Power(u,2)/2. + Power(u,3)/6. + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		
		double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
		double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
		double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
	
		rp[0] = 0;
		rp[1] = 0;
		rp[2] = 0;
	
		double lru[3] = {0,0,0};
		double lrv[3] = {0,0,0};
		int ncoords_base = theFormulas[frm].ncoor;
		for( int x = 0; x < ncoords_base; x++ )
		{
			double *lr = r + cp[x]*3;
	
			rp[0] += ceff_map[x] * (lr[0] + theFormulas[frm].r_pbc[3*x+0]);
			rp[1] += ceff_map[x] * (lr[1] + theFormulas[frm].r_pbc[3*x+1]);
			rp[2] += ceff_map[x] * (lr[2] + theFormulas[frm].r_pbc[3*x+2]);
			
			lru[0] += ceff_map_du[x] * (lr[0] + theFormulas[frm].r_pbc[3*x+0]);
			lru[1] += ceff_map_du[x] * (lr[1] + theFormulas[frm].r_pbc[3*x+1]);
			lru[2] += ceff_map_du[x] * (lr[2] + theFormulas[frm].r_pbc[3*x+2]);
			
			lrv[0] += ceff_map_dv[x] * (lr[0] + theFormulas[frm].r_pbc[3*x+0]);
			lrv[1] += ceff_map_dv[x] * (lr[1] + theFormulas[frm].r_pbc[3*x+1]);
			lrv[2] += ceff_map_dv[x] * (lr[2] + theFormulas[frm].r_pbc[3*x+2]);
		}
	
		rp[0] *= alpha_x;
		rp[1] *= alpha_y;
		rp[2] *= alpha_z;
		
		cross( lru, lrv, nrm );
		normalize(nrm);
	}	
}


void surface::sortFaces( void )
{
	// decide whether to flip the orientation. first triangle normal points away from origin.
	for( int i = 0; i < nv; i++ )
	{
		int val = theVertices[i].valence;

		if( val != theVertices[i].nfaces )
		{
			printf("surface::sortFaces called with valence != nfaces.\n");
			exit(1);
		}

		int sorted_faces[val];
		for( int e = 0; e < val; e++ )
		{
			int ep1 = e+1;
			if( ep1 >= val ) ep1 -= val;

			int j = theVertices[i].edges[e];
			int k = theVertices[i].edges[ep1];

			int arr[3] = {i,j,k};

			sort3(arr);

			int gotit = 0;
			for( int fx = 0; fx < val; fx++ )
			{
				int f = theVertices[i].faces[fx];
				if( theTriangles[f].ids[0] == arr[0] && theTriangles[f].ids[1] == arr[1] && theTriangles[f].ids[2] == arr[2] )
				{
					sorted_faces[e] = f;
					gotit = 1;
				}
			}
			if( !gotit )
			{
				printf("Couldn't find face for i,j,k triplet.\n");
				exit(1);
			}

		} 
		memcpy( theVertices[i].faces, sorted_faces, sizeof(int) * val );
	}	
	


	generateBorderMappings();
}

double surface::penergy( int f, double *r, double *p_uv, int do_monge, double face_energy_density )
{
	double alpha_x = r[3*nv];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];
		
	double en = 0;

	if( f >= nf_faces ) // an irregular face.
	{
		int frm = (f-nf_faces)*nf_irr_pts;

		int t = theIrregularFormulas[frm].tri;
		int i = theIrregularFormulas[frm].vertex;
		int e = theIrregularFormulas[frm].edge;

		int val = theVertices[i].valence;

		irr_kernel *theKernel = kernels[val]; 
		int ncoords_base = val + 6;
		
		for( int px = 0; px < theTriangles[t].np; px++ )
		{
			int pid = theTriangles[t].plist[px];
			double p_c0 = theTriangles[t].pc0[px];
				
//			NOTE 1: at some point I had a good reason not to calculate this.
//			if( fabs(p_c0 - theIrregularFormulas[frm].c0) < 1e-7 )
//				continue; 
			double p_area = theTriangles[t].pa[px];
	
			double fu = p_uv[2*pid+0];
			double fv = p_uv[2*pid+1];					
	
			int domain = theKernel->domain(fu,fv);
			if( domain-1 >= theKernel->ndomains )
			{
				printf("LOGICAL ERROR evaluating R,nrm on irregular vertex face.\n");
				exit(1);
			}	
			double *theMap = theKernel->get_map( &fu, &fv );
	
			double u = fu;
			double v = fv;
			double w = 1 - u - v;
	
			if( u +v > 1.0 + THRESH || u < -THRESH || v < -THRESH )
			{
				printf("ERROR: outside domain.\n");
				exit(1);
			}
	
			double u2 = u*u;
			double u3 = u*u*u;
			double u4 = u*u*u*u;
			
			double v2 = v*v;
			double v3 = v*v*v;
			double v4 = v*v*v*v;
			
			double w2 = w*w;
			double w3 = w*w*w;
			double w4 = w*w*w*w;
				
			
			// 8 : 0
			// 7 : 1
			// 4 : 2
			// 5 : 3
			// 9 : 4
			// 12 : 5
			// 11 : 6
			// 10 : 7
			// 6 : 8
			// 3 : 9
			// 1 : 10
			// 2 : 11
	
		
	
			double n1 = (1.0/12.0)*(u4+2*u3*v); 
			double n2 = (1.0/12.0)*(u4+2*u3*w); 
			double n3 = (1.0/12.0)*(u4+2*u3*w+6*u3*v+6*u2*v*w+12*u2*v2+6*u*v2*w+6*u*v3+2*v3*w+v4); 
			double n4 = (1.0/12.0)*(6*u4+24*u3*w+24*u2*w2+8*u*w3+w4+24*u3*v+60*u2*v*w+36*u*v*w2+6*v*w3+24*u2*v2+36*u*v2*w+12*v2*w2+8*u*v3+6*v3*w+v4); 
			double n5 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+2*u3*v+6*u2*v*w+6*u*v*w2+2*v*w3); 
			double n6 = (1.0/12.0)*(2*u*v3+v4); 
			double n7 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+8*u3*v+36*u2*v*w+36*u*v*w2+8*v*w3+24*u2*v2+60*u*v2*w+24*v2*w2+24*u*v3+24*v3*w+6*v4); 
			double n8 = (1.0/12.0)*(u4+8*u3*w+24*u2*w2+24*u*w3+6*w4+6*u3*v+36*u2*v*w+60*u*v*w2+24*v*w3+12*u2*v2+36*u*v2*w+24*v2*w2+6*u*v3+8*v3*w+v4); 
			double n9 = (1.0/12.0)*(2*u*w3+w4); 
			double n10 = (1.0/12.0)*(2*v3*w+v4); 
			double n11 = (1.0/12.0)*(2*u*w3+w4+6*u*v*w2+6*v*w3+6*u*v2*w+12*v2*w2+2*u*v3+6*v3*w+v4); 
			double n12 = (1.0/12.0)*(w4+2*v*w3);
	
			double du_1 = Power(u,3)/3. + (Power(u,2)*v)/2.;
			double du_2 = Power(u,2)/2. - Power(u,3)/3. - (Power(u,2)*v)/2.; 	
			double du_3 =Power(u,2)/2. - Power(u,3)/3. + u*v - (Power(u,2)*v)/2. + Power(v,2)/2. - Power(v,3)/6.;
			double du_4 = 0.3333333333333333 + u - Power(u,2) - Power(u,3)/3. + v/2. - u*v - (Power(u,2)*v)/2. - Power(v,2) + Power(v,3)/3.;	
			double du_5 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. - v/2. + Power(u,2)*v + Power(v,2)/2. - Power(v,3)/6.;
			double du_6 = Power(v,3)/6.;	
			double du_7 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. + v/2. - 2*u*v + Power(u,2)*v - Power(v,2)/2. - Power(v,3)/6.;	
			double du_8 = -2*u + 2*Power(u,2) - Power(u,3)/3. - v + 2*u*v - (Power(u,2)*v)/2. + Power(v,2) - Power(v,3)/6.;	
			double du_9 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. + v/2. - (Power(u,2)*v)/2. - Power(v,2)/2. + Power(v,3)/6.;	
			double du_10 = -Power(v,3)/6.;
			double du_11 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. - v/2. + u*v - (Power(u,2)*v)/2. + Power(v,3)/3.;	
			double du_12 = 	-0.3333333333333333 + u - Power(u,2) + Power(u,3)/3. + v/2. - u*v + (Power(u,2)*v)/2. - Power(v,3)/6.;
	
			double dv_1 = Power(u,3)/6.;
			double dv_2 = -Power(u,3)/6.;
			double dv_3 = Power(u,2)/2. - Power(u,3)/6. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_4 = 0.16666666666666666 + u/2. - Power(u,2)/2. - Power(u,3)/6. - 2*u*v - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
			double dv_5 = -0.16666666666666666 - u/2. + Power(u,3)/3. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_6 = (u*Power(v,2))/2. + Power(v,3)/3.;
			double dv_7 = 0.3333333333333333 + u/2. - Power(u,2) + Power(u,3)/3. + v - u*v - Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_8 = -u + Power(u,2) - Power(u,3)/6. - 2*v + 2*u*v + 2*Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_9 = -0.3333333333333333 + u/2. - Power(u,3)/6. + v - u*v - Power(v,2) + (u*Power(v,2))/2. + Power(v,3)/3.;
			double dv_10 = Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_11 = 0.16666666666666666 - u/2. + Power(u,2)/2. - Power(u,3)/6. - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
			double dv_12 = -0.16666666666666666 + u/2. - Power(u,2)/2. + Power(u,3)/6. + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		
			double d_uu_1 = Power(u,2) + u*v;
			double d_uu_2 = u - Power(u,2) - u*v;
			double d_uu_3 = u - Power(u,2) + v - u*v;
			double d_uu_4 = 1 - 2*u - Power(u,2) - v - u*v;
			double d_uu_5 = -2*u + 2*Power(u,2) + 2*u*v;
			double d_uu_6 = 0;
			double d_uu_7 = -2*u + 2*Power(u,2) - 2*v + 2*u*v;
			double d_uu_8 = -2 + 4*u - Power(u,2) + 2*v - u*v;
			double d_uu_9 = u - Power(u,2) - u*v;
			double d_uu_10 = 0;
			double d_uu_11 = u - Power(u,2) + v - u*v;
			double d_uu_12 = 1 - 2*u + Power(u,2) - v + u*v;
			
			double d_uv_1 = Power(u,2)/2.;
			double d_uv_2 = -Power(u,2)/2.;
			double d_uv_3 = u - Power(u,2)/2. + v - Power(v,2)/2.;
			double d_uv_4 = 0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
			double d_uv_5 = -0.5 + Power(u,2) + v - Power(v,2)/2.;
			double d_uv_6 = Power(v,2)/2.;
			double d_uv_7 = 0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
			double d_uv_8 = -1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
			double d_uv_9 = 0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
			double d_uv_10 = -Power(v,2)/2.;
			double d_uv_11 = -0.5 + u - Power(u,2)/2. + Power(v,2);
			double d_uv_12 = 0.5 - u + Power(u,2)/2. - Power(v,2)/2.;
			
			double d_vv_1 = 0;
			double d_vv_2 = 0;
			double d_vv_3 = u + v - u*v - Power(v,2);
			double d_vv_4 = -2*u - 2*v + 2*u*v + 2*Power(v,2);
			double d_vv_5 = u + v - u*v - Power(v,2);
			double d_vv_6 = u*v + Power(v,2);
			double d_vv_7 = 1 - u - 2*v - u*v - Power(v,2);
			double d_vv_8 = -2 + 2*u + 4*v - u*v - Power(v,2);
			double d_vv_9 = 1 - u - 2*v + u*v + Power(v,2);
			double d_vv_10 = v - u*v - Power(v,2);
			double d_vv_11 = -2*v + 2*u*v + 2*Power(v,2);
			double d_vv_12 = v - u*v - Power(v,2);
			
			double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
			double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
			double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
			
			double ceff_map_duu[12] = { d_uu_8, d_uu_7, d_uu_4, d_uu_5, d_uu_9, d_uu_12, d_uu_11, d_uu_10, d_uu_6, d_uu_3, d_uu_1, d_uu_2 };
			double ceff_map_duv[12] = { d_uv_8, d_uv_7, d_uv_4, d_uv_5, d_uv_9, d_uv_12, d_uv_11, d_uv_10, d_uv_6, d_uv_3, d_uv_1, d_uv_2 };
			double ceff_map_dvv[12] = { d_vv_8, d_vv_7, d_vv_4, d_vv_5, d_vv_9, d_vv_12, d_vv_11, d_vv_10, d_vv_6, d_vv_3, d_vv_1, d_vv_2 };
	
			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
			double nrm[3]={0,0,0}; 
	
			int *cset = theVertices[i].irr_coord_set + e * ncoords_base;
	
			for( int x = 0; x < ncoords_base; x++ )
			{
				for( int y = 0; y < 12; y++ )
				{
					R[0] += (r[3*cset[x]+0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map[y] * alpha_x;
					R[1] += (r[3*cset[x]+1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map[y] * alpha_y;
					R[2] += (r[3*cset[x]+2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map[y] * alpha_z;
					
					Ru[0] += (r[3*cset[x]+0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_du[y] * alpha_x;
					Ru[1] += (r[3*cset[x]+1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_du[y] * alpha_y;
					Ru[2] += (r[3*cset[x]+2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_du[y] * alpha_z;
					
					Rv[0] += (r[3*cset[x]+0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_dv[y] * alpha_x;
					Rv[1] += (r[3*cset[x]+1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_dv[y] * alpha_y;
					Rv[2] += (r[3*cset[x]+2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_dv[y] * alpha_z;
					
					tSuu[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_duu[y] * alpha_x;
					tSuu[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_duu[y] * alpha_y;
					tSuu[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_duu[y] * alpha_z;
					
					tSuv[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_duv[y] * alpha_x;
					tSuv[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_duv[y] * alpha_y;
					tSuv[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_duv[y] * alpha_z;
					
					tSvv[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_dvv[y] * alpha_x;
					tSvv[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_dvv[y] * alpha_y;
					tSvv[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_dvv[y] * alpha_z;
				}
			}
	
			cross( Ru, Rv, nrm );
			normalize(nrm);
		
		
			double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
			double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
			double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];
		
			double g = sqrt(RuRu*RvRv-RuRv*RuRv);
		
			double c1,c2;
		
			if( do_monge )
			{
			      double idet = 1.0 / (Ru[0] * Rv[1] - Ru[1] * Rv[0]);
		
				double dudx=idet*Rv[1];
				double dudy=-idet*Rv[0];
				double dvdx=-idet*Ru[1];
				double dvdy=idet*Ru[0];
		
				double dhdx=Ru[2]*dudx+Rv[2]*dvdx;
				double dhdy=Ru[2]*dudy+Rv[2]*dvdy;
		
				c1 = tSuu[2] * dudx * dudx + tSuv[2] * dudx * dvdx + tSuv[2] * dvdx * dudx + tSvv[2] * dvdx * dvdx;
				c2 = tSuu[2] * dudy * dudy + tSuv[2] * dudy * dvdy + tSuv[2] * dvdy * dudy + tSvv[2] * dvdy * dvdy; 
			}
			else
			{
				double nsuu = tSuu[0] * nrm[0] + tSuu[1] * nrm[1] + tSuu[2] * nrm[2];
				double nsuv = tSuv[0] * nrm[0] + tSuv[1] * nrm[1] + tSuv[2] * nrm[2];
				double nsvv = tSvv[0] * nrm[0] + tSvv[1] * nrm[1] + tSvv[2] * nrm[2];
		
				double Stot = (nsuu * RvRv + nsvv * RuRu -2*nsuv*RuRv)/(g*g);
		
				double Sop[4] = { 1.0/(g*g) * (nsuu * RvRv  - nsuv * RuRv), (1.0/(g*g)) * ( nsuv*RvRv-nsvv*RuRv),
						  1.0/(g*g) * (nsuv * RuRu  - nsuu * RuRv), (1.0/(g*g)) * ( nsvv*RuRu-nsuv*RuRv) };
		
				double a = Sop[0];		
				double b = Sop[1];
				double c = Sop[2];
				double d = Sop[3];
		
				double c0 = theIrregularFormulas[frm].c0;
				c1 = -0.5*(a+d-sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
				c2 = -0.5*(a+d+sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
			}
	
				
			en += 0.5*kc * p_area * ( c1 + c2 - p_c0 ) * ( c1 + c2 - p_c0 );

			//en -= 0.5*kc * p_area * ( c1 + c2 - c0 ) * ( c1 + c2 - c0 );

			en -= p_area * face_energy_density;
		}
	}
	else
	{
		int frm = f*nf_g_q_p;
		int t = theFormulas[frm].tri;
	
		for( int px = 0; px < theTriangles[t].np; px++ )
		{
			int pid = theTriangles[t].plist[px];
			double p_c0 = theTriangles[t].pc0[px];
//			NOTE 1: at some point I had a good reason not to calculate this.
//			if( fabs(p_c0 - theFormulas[frm].c0) < 1e-7 )
//				continue; 
			double p_area = theTriangles[t].pa[px];
	
			double u = p_uv[2*pid+0];
			double v = p_uv[2*pid+1];
	
			int *cp = theFormulas[frm].cp;
	
			double w = 1-u-v;
	
			double u2 = u*u;
			double u3 = u*u*u;
			double u4 = u*u*u*u;
			
			double v2 = v*v;
			double v3 = v*v*v;
			double v4 = v*v*v*v;
			
			double w2 = w*w;
			double w3 = w*w*w;
			double w4 = w*w*w*w;
	
			double n1 = (1.0/12.0)*(u4+2*u3*v); 
			double n2 = (1.0/12.0)*(u4+2*u3*w); 
			double n3 = (1.0/12.0)*(u4+2*u3*w+6*u3*v+6*u2*v*w+12*u2*v2+6*u*v2*w+6*u*v3+2*v3*w+v4); 
			double n4 = (1.0/12.0)*(6*u4+24*u3*w+24*u2*w2+8*u*w3+w4+24*u3*v+60*u2*v*w+36*u*v*w2+6*v*w3+24*u2*v2+36*u*v2*w+12*v2*w2+8*u*v3+6*v3*w+v4); 
			double n5 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+2*u3*v+6*u2*v*w+6*u*v*w2+2*v*w3); 
			double n6 = (1.0/12.0)*(2*u*v3+v4); 
			double n7 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+8*u3*v+36*u2*v*w+36*u*v*w2+8*v*w3+24*u2*v2+60*u*v2*w+24*v2*w2+24*u*v3+24*v3*w+6*v4); 
			double n8 = (1.0/12.0)*(u4+8*u3*w+24*u2*w2+24*u*w3+6*w4+6*u3*v+36*u2*v*w+60*u*v*w2+24*v*w3+12*u2*v2+36*u*v2*w+24*v2*w2+6*u*v3+8*v3*w+v4); 
			double n9 = (1.0/12.0)*(2*u*w3+w4); 
			double n10 = (1.0/12.0)*(2*v3*w+v4); 
			double n11 = (1.0/12.0)*(2*u*w3+w4+6*u*v*w2+6*v*w3+6*u*v2*w+12*v2*w2+2*u*v3+6*v3*w+v4); 
			double n12 = (1.0/12.0)*(w4+2*v*w3);
	
			double du_1 = Power(u,3)/3. + (Power(u,2)*v)/2.;
			double du_2 = Power(u,2)/2. - Power(u,3)/3. - (Power(u,2)*v)/2.; 	
			double du_3 = Power(u,2)/2. - Power(u,3)/3. + u*v - (Power(u,2)*v)/2. + Power(v,2)/2. - Power(v,3)/6.;
			double du_4 = 0.3333333333333333 + u - Power(u,2) - Power(u,3)/3. + v/2. - u*v - (Power(u,2)*v)/2. - Power(v,2) + Power(v,3)/3.;	
			double du_5 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. - v/2. + Power(u,2)*v + Power(v,2)/2. - Power(v,3)/6.;
			double du_6 = Power(v,3)/6.;	
			double du_7 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. + v/2. - 2*u*v + Power(u,2)*v - Power(v,2)/2. - Power(v,3)/6.;	
			double du_8 = -2*u + 2*Power(u,2) - Power(u,3)/3. - v + 2*u*v - (Power(u,2)*v)/2. + Power(v,2) - Power(v,3)/6.;	
			double du_9 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. + v/2. - (Power(u,2)*v)/2. - Power(v,2)/2. + Power(v,3)/6.;	
			double du_10 = -Power(v,3)/6.;
			double du_11 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. - v/2. + u*v - (Power(u,2)*v)/2. + Power(v,3)/3.;	
			double du_12 = 	-0.3333333333333333 + u - Power(u,2) + Power(u,3)/3. + v/2. - u*v + (Power(u,2)*v)/2. - Power(v,3)/6.;
	
			double dv_1 = Power(u,3)/6.;
			double dv_2 = -Power(u,3)/6.;
			double dv_3 = Power(u,2)/2. - Power(u,3)/6. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_4 = 0.16666666666666666 + u/2. - Power(u,2)/2. - Power(u,3)/6. - 2*u*v - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
			double dv_5 = -0.16666666666666666 - u/2. + Power(u,3)/3. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_6 = (u*Power(v,2))/2. + Power(v,3)/3.;
			double dv_7 = 0.3333333333333333 + u/2. - Power(u,2) + Power(u,3)/3. + v - u*v - Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_8 = -u + Power(u,2) - Power(u,3)/6. - 2*v + 2*u*v + 2*Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_9 = -0.3333333333333333 + u/2. - Power(u,3)/6. + v - u*v - Power(v,2) + (u*Power(v,2))/2. + Power(v,3)/3.;
			double dv_10 = Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_11 = 0.16666666666666666 - u/2. + Power(u,2)/2. - Power(u,3)/6. - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
			double dv_12 = -0.16666666666666666 + u/2. - Power(u,2)/2. + Power(u,3)/6. + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
	
			double d_uu_1 = Power(u,2) + u*v;
			double d_uu_2 = u - Power(u,2) - u*v;
			double d_uu_3 = u - Power(u,2) + v - u*v;
			double d_uu_4 = 1 - 2*u - Power(u,2) - v - u*v;
			double d_uu_5 = -2*u + 2*Power(u,2) + 2*u*v;
			double d_uu_6 = 0;
			double d_uu_7 = -2*u + 2*Power(u,2) - 2*v + 2*u*v;
			double d_uu_8 = -2 + 4*u - Power(u,2) + 2*v - u*v;
			double d_uu_9 = u - Power(u,2) - u*v;
			double d_uu_10 = 0;
			double d_uu_11 = u - Power(u,2) + v - u*v;
			double d_uu_12 = 1 - 2*u + Power(u,2) - v + u*v;
			
			double d_uv_1 = Power(u,2)/2.;
			double d_uv_2 = -Power(u,2)/2.;
			double d_uv_3 = u - Power(u,2)/2. + v - Power(v,2)/2.;
			double d_uv_4 = 0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
			double d_uv_5 = -0.5 + Power(u,2) + v - Power(v,2)/2.;
			double d_uv_6 = Power(v,2)/2.;
			double d_uv_7 = 0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
			double d_uv_8 = -1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
			double d_uv_9 = 0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
			double d_uv_10 = -Power(v,2)/2.;
			double d_uv_11 = -0.5 + u - Power(u,2)/2. + Power(v,2);
			double d_uv_12 = 0.5 - u + Power(u,2)/2. - Power(v,2)/2.;
			
			double d_vv_1 = 0;
			double d_vv_2 = 0;
			double d_vv_3 = u + v - u*v - Power(v,2);
			double d_vv_4 = -2*u - 2*v + 2*u*v + 2*Power(v,2);
			double d_vv_5 = u + v - u*v - Power(v,2);
			double d_vv_6 = u*v + Power(v,2);
			double d_vv_7 = 1 - u - 2*v - u*v - Power(v,2);
			double d_vv_8 = -2 + 2*u + 4*v - u*v - Power(v,2);
			double d_vv_9 = 1 - u - 2*v + u*v + Power(v,2);
			double d_vv_10 = v - u*v - Power(v,2);
			double d_vv_11 = -2*v + 2*u*v + 2*Power(v,2);
			double d_vv_12 = v - u*v - Power(v,2);
		
			double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
			double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
			double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
			
			double ceff_map_duu[12] = { d_uu_8, d_uu_7, d_uu_4, d_uu_5, d_uu_9, d_uu_12, d_uu_11, d_uu_10, d_uu_6, d_uu_3, d_uu_1, d_uu_2 };
			double ceff_map_duv[12] = { d_uv_8, d_uv_7, d_uv_4, d_uv_5, d_uv_9, d_uv_12, d_uv_11, d_uv_10, d_uv_6, d_uv_3, d_uv_1, d_uv_2 };
			double ceff_map_dvv[12] = { d_vv_8, d_vv_7, d_vv_4, d_vv_5, d_vv_9, d_vv_12, d_vv_11, d_vv_10, d_vv_6, d_vv_3, d_vv_1, d_vv_2 };
					
			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
			double nrm[3]={0,0,0}; 
	
			int np = theFormulas[frm].ncoor;
	
			for( int p = 0; p < np; p++ )
			{
	
				R[0] += ceff_map[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				R[1] += ceff_map[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				R[2] += ceff_map[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				Ru[0] += ceff_map_du[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				Ru[1] += ceff_map_du[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				Ru[2] += ceff_map_du[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				Rv[0] += ceff_map_dv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				Rv[1] += ceff_map_dv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				Rv[2] += ceff_map_dv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuu[0] += ceff_map_duu[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuu[1] += ceff_map_duu[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuu[2] += ceff_map_duu[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuv[0] += ceff_map_duv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuv[1] += ceff_map_duv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuv[2] += ceff_map_duv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSvv[0] += ceff_map_dvv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSvv[1] += ceff_map_dvv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSvv[2] += ceff_map_dvv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			}
			cross( Ru, Rv, nrm );
			normalize(nrm);
	
	
			double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
			double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
			double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];
	
			double g = sqrt(RuRu*RvRv-RuRv*RuRv);
	
			double c1,c2;
	      
			if( do_monge )
			{
			      double idet = 1.0 / (Ru[0] * Rv[1] - Ru[1] * Rv[0]);
	
				double dudx=idet*Rv[1];
				double dudy=-idet*Rv[0];
				double dvdx=-idet*Ru[1];
				double dvdy=idet*Ru[0];
	
				double dhdx=Ru[2]*dudx+Rv[2]*dvdx;
				double dhdy=Ru[2]*dudy+Rv[2]*dvdy;
	      
				c1 = tSuu[2] * dudx * dudx + tSuv[2] * dudx * dvdx + tSuv[2] * dvdx * dudx + tSvv[2] * dvdx * dvdx;
				c2 = tSuu[2] * dudy * dudy + tSuv[2] * dudy * dvdy + tSuv[2] * dvdy * dudy + tSvv[2] * dvdy * dvdy; 
			}
			else
			{
				double nsuu = tSuu[0] * nrm[0] + tSuu[1] * nrm[1] + tSuu[2] * nrm[2];
				double nsuv = tSuv[0] * nrm[0] + tSuv[1] * nrm[1] + tSuv[2] * nrm[2];
				double nsvv = tSvv[0] * nrm[0] + tSvv[1] * nrm[1] + tSvv[2] * nrm[2];
		
				double Stot = (nsuu * RvRv + nsvv * RuRu -2*nsuv*RuRv)/(g*g);
		
				double Sop[4] = { 1.0/(g*g) * (nsuu * RvRv  - nsuv * RuRv), (1.0/(g*g)) * ( nsuv*RvRv-nsvv*RuRv),
						  1.0/(g*g) * (nsuv * RuRu  - nsuu * RuRv), (1.0/(g*g)) * ( nsvv*RuRu-nsuv*RuRv) };
		
				double a = Sop[0];		
				double b = Sop[1];
				double c = Sop[2];
				double d = Sop[3];
		
				double c0 = theFormulas[frm].c0;
				c1 = -0.5*(a+d-sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
				c2 = -0.5*(a+d+sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
			}
	

			en += 0.5*kc * p_area * ( c1 + c2 - p_c0 ) * ( c1 + c2 - p_c0 );
//			en -= 0.5*kc * p_area * ( c1 + c2 - c0 ) * ( c1 + c2 - c0 );
			en -= p_area * face_energy_density;

		}
	}
	return en;
}

void surface::constructIrregularKernels( void )
{
	memset( kernels, 0, sizeof(irr_kernel *) * 20 );

	for( int i = 0; i < nv; i++ )
	{
		int val = theVertices[i].valence;

		if( val != 6 )
		{
			if( kernels[val] == NULL )
			{
				kernels[val] = (irr_kernel *)malloc( sizeof(irr_kernel) );
				kernels[val]->setup( val, MAX_DOMAIN ); 
			}

			int ncoords_base = kernels[val]->ncoords_base;

			theVertices[i].irr_coord_set = (int *)malloc( sizeof(int) * val * (ncoords_base) );

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
				printf("%d valence %d, %d valence %d, %d valence %d.\n",
				i, theVertices[i].valence,
				j, theVertices[j].valence,
				k, theVertices[k].valence );

					printf("Illegal mesh: must be subdivided once isolating low index on irregular mesh.\n");
					exit(1);
				}
	
				int ncoords_base  = 1 + val + 5; 
				// we also need to "bring along" coordinates for the other three triangular splines in the final pass.
				// but we only need to do the eigenvector analysis on the first pass.
				int point_list[ncoords_base];
						
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

				memcpy( theVertices[i].irr_coord_set+e*(ncoords_base), point_list, sizeof(int) * ncoords_base ); 

				// test it here.

				irr_kernel *theKernel = kernels[val];

				for( double tu = 0.0011568313814261763; tu <= 1.0; tu *= 1.8 )
				for( double tv = 0.000001; tv <= 1.0; tv *= 1.8 )
				{
					double fu = tu;
					double fv = tv;

					if( fu+fv > 1.0 ) 
						continue;
					int domain = theKernel->domain(fu,fv);
					if( domain-1 >= theKernel->ndomains )
						continue;
				
					double *theMap = theKernel->get_map( &fu, &fv );

					double u = fu;
					double v = fv;
					double w = 1 - u - v;

					if( u +v > 1.0 + THRESH || u < -THRESH || v < -THRESH  )
					{
						printf("ERROR: outside domain.\n");
						exit(1);
					}

					double u2 = u*u;
					double u3 = u*u*u;
					double u4 = u*u*u*u;
					
					double v2 = v*v;
					double v3 = v*v*v;
					double v4 = v*v*v*v;
			
					double w2 = w*w;
					double w3 = w*w*w;
					double w4 = w*w*w*w;
						
					
					// 8 : 0
					// 7 : 1
					// 4 : 2
					// 5 : 3
					// 9 : 4
					// 12 : 5
					// 11 : 6
					// 10 : 7
					// 6 : 8
					// 3 : 9
					// 1 : 10
					// 2 : 11

	

					double n1 = (1.0/12.0)*(u4+2*u3*v); 
					double n2 = (1.0/12.0)*(u4+2*u3*w); 
					double n3 = (1.0/12.0)*(u4+2*u3*w+6*u3*v+6*u2*v*w+12*u2*v2+6*u*v2*w+6*u*v3+2*v3*w+v4); 
					double n4 = (1.0/12.0)*(6*u4+24*u3*w+24*u2*w2+8*u*w3+w4+24*u3*v+60*u2*v*w+36*u*v*w2+6*v*w3+24*u2*v2+36*u*v2*w+12*v2*w2+8*u*v3+6*v3*w+v4); 
					double n5 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+2*u3*v+6*u2*v*w+6*u*v*w2+2*v*w3); 
					double n6 = (1.0/12.0)*(2*u*v3+v4); 
					double n7 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+8*u3*v+36*u2*v*w+36*u*v*w2+8*v*w3+24*u2*v2+60*u*v2*w+24*v2*w2+24*u*v3+24*v3*w+6*v4); 
					double n8 = (1.0/12.0)*(u4+8*u3*w+24*u2*w2+24*u*w3+6*w4+6*u3*v+36*u2*v*w+60*u*v*w2+24*v*w3+12*u2*v2+36*u*v2*w+24*v2*w2+6*u*v3+8*v3*w+v4); 
					double n9 = (1.0/12.0)*(2*u*w3+w4); 
					double n10 = (1.0/12.0)*(2*v3*w+v4); 
					double n11 = (1.0/12.0)*(2*u*w3+w4+6*u*v*w2+6*v*w3+6*u*v2*w+12*v2*w2+2*u*v3+6*v3*w+v4); 
					double n12 = (1.0/12.0)*(w4+2*v*w3);

					double n1_alt = Power(u,4)/12. + (Power(u,3)*v)/6.;		
					double n2_alt = Power(u,4)/12. + (Power(u,3)*w)/6.;
					double n3_alt = Power(u,4)/12. + (Power(u,3)*v)/2. + Power(u,2)*Power(v,2) + (u*Power(v,3))/2. + Power(v,4)/12. + (Power(u,3)*w)/6. + (Power(u,2)*v*w)/2. + (u*Power(v,2)*w)/2. + (Power(v,3)*w)/6.;
					double n4_alt = Power(u,4)/2. + 2*Power(u,3)*v + 2*Power(u,2)*Power(v,2) + (2*u*Power(v,3))/3. + Power(v,4)/12. + 2*Power(u,3)*w + 5*Power(u,2)*v*w + 3*u*Power(v,2)*w + (Power(v,3)*w)/2. + 2*Power(u,2)*Power(w,2) + 3*u*v*Power(w,2) + Power(v,2)*Power(w,2) + (2*u*Power(w,3))/3. + (v*Power(w,3))/2. + Power(w,4)/12.;
					double n5_alt = Power(u,4)/12. + (Power(u,3)*v)/6. + (Power(u,3)*w)/2. + (Power(u,2)*v*w)/2. + Power(u,2)*Power(w,2) + (u*v*Power(w,2))/2. + (u*Power(w,3))/2. + (v*Power(w,3))/6. + Power(w,4)/12.;
					double n6_alt = (u*Power(v,3))/6. + Power(v,4)/12.;
					double n7_alt = Power(u,4)/12. + (2*Power(u,3)*v)/3. + 2*Power(u,2)*Power(v,2) + 2*u*Power(v,3) + Power(v,4)/2. + (Power(u,3)*w)/2. + 3*Power(u,2)*v*w + 5*u*Power(v,2)*w + 2*Power(v,3)*w + Power(u,2)*Power(w,2) + 3*u*v*Power(w,2) + 2*Power(v,2)*Power(w,2) + (u*Power(w,3))/2. + (2*v*Power(w,3))/3. + Power(w,4)/12.;
					double n8_alt = Power(u,4)/12. + (Power(u,3)*v)/2. + Power(u,2)*Power(v,2) + (u*Power(v,3))/2. + Power(v,4)/12. + (2*Power(u,3)*w)/3. + 3*Power(u,2)*v*w + 3*u*Power(v,2)*w + (2*Power(v,3)*w)/3. + 2*Power(u,2)*Power(w,2) + 5*u*v*Power(w,2) + 2*Power(v,2)*Power(w,2) + 2*u*Power(w,3) + 2*v*Power(w,3) + Power(w,4)/2.;
					double n9_alt = (u*Power(w,3))/6. + Power(w,4)/12.;
					double n10_alt = Power(v,4)/12. + (Power(v,3)*w)/6.;
					double n11_alt =(u*Power(v,3))/6. + Power(v,4)/12. + (u*Power(v,2)*w)/2. + (Power(v,3)*w)/2. + (u*v*Power(w,2))/2. + Power(v,2)*Power(w,2) + (u*Power(w,3))/6. + (v*Power(w,3))/2. + Power(w,4)/12.;
					double n12_alt = (v*Power(w,3))/6. + Power(w,4)/12.;
	
					double du_1 = Power(u,3)/3. + (Power(u,2)*v)/2.;
					double du_2 = Power(u,2)/2. - Power(u,3)/3. - (Power(u,2)*v)/2.; 	
					double du_3 =Power(u,2)/2. - Power(u,3)/3. + u*v - (Power(u,2)*v)/2. + Power(v,2)/2. - Power(v,3)/6.;
					double du_4 = 0.3333333333333333 + u - Power(u,2) - Power(u,3)/3. + v/2. - u*v - (Power(u,2)*v)/2. - Power(v,2) + Power(v,3)/3.;	
					double du_5 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. - v/2. + Power(u,2)*v + Power(v,2)/2. - Power(v,3)/6.;
					double du_6 = Power(v,3)/6.;	
					double du_7 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. + v/2. - 2*u*v + Power(u,2)*v - Power(v,2)/2. - Power(v,3)/6.;	
					double du_8 = -2*u + 2*Power(u,2) - Power(u,3)/3. - v + 2*u*v - (Power(u,2)*v)/2. + Power(v,2) - Power(v,3)/6.;	
					double du_9 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. + v/2. - (Power(u,2)*v)/2. - Power(v,2)/2. + Power(v,3)/6.;	
					double du_10 = -Power(v,3)/6.;
					double du_11 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. - v/2. + u*v - (Power(u,2)*v)/2. + Power(v,3)/3.;	
					double du_12 = 	-0.3333333333333333 + u - Power(u,2) + Power(u,3)/3. + v/2. - u*v + (Power(u,2)*v)/2. - Power(v,3)/6.;

					double dv_1 = Power(u,3)/6.;
					double dv_2 = -Power(u,3)/6.;
					double dv_3 = Power(u,2)/2. - Power(u,3)/6. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
					double dv_4 = 0.16666666666666666 + u/2. - Power(u,2)/2. - Power(u,3)/6. - 2*u*v - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
					double dv_5 = -0.16666666666666666 - u/2. + Power(u,3)/3. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
					double dv_6 = (u*Power(v,2))/2. + Power(v,3)/3.;
					double dv_7 = 0.3333333333333333 + u/2. - Power(u,2) + Power(u,3)/3. + v - u*v - Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
					double dv_8 = -u + Power(u,2) - Power(u,3)/6. - 2*v + 2*u*v + 2*Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
					double dv_9 = -0.3333333333333333 + u/2. - Power(u,3)/6. + v - u*v - Power(v,2) + (u*Power(v,2))/2. + Power(v,3)/3.;
					double dv_10 = Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
					double dv_11 = 0.16666666666666666 - u/2. + Power(u,2)/2. - Power(u,3)/6. - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
					double dv_12 = -0.16666666666666666 + u/2. - Power(u,2)/2. + Power(u,3)/6. + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
	
					double d_uu_1 = Power(u,2) + u*v;
					double d_uu_2 = u - Power(u,2) - u*v;
					double d_uu_3 = u - Power(u,2) + v - u*v;
					double d_uu_4 = 1 - 2*u - Power(u,2) - v - u*v;
					double d_uu_5 = -2*u + 2*Power(u,2) + 2*u*v;
					double d_uu_6 = 0;
					double d_uu_7 = -2*u + 2*Power(u,2) - 2*v + 2*u*v;
					double d_uu_8 = -2 + 4*u - Power(u,2) + 2*v - u*v;
					double d_uu_9 = u - Power(u,2) - u*v;
					double d_uu_10 = 0;
					double d_uu_11 = u - Power(u,2) + v - u*v;
					double d_uu_12 = 1 - 2*u + Power(u,2) - v + u*v;
					
					double d_uv_1 = Power(u,2)/2.;
					double d_uv_2 = -Power(u,2)/2.;
					double d_uv_3 = u - Power(u,2)/2. + v - Power(v,2)/2.;
					double d_uv_4 = 0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
					double d_uv_5 = -0.5 + Power(u,2) + v - Power(v,2)/2.;
					double d_uv_6 = Power(v,2)/2.;
					double d_uv_7 = 0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
					double d_uv_8 = -1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
					double d_uv_9 = 0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
					double d_uv_10 = -Power(v,2)/2.;
					double d_uv_11 = -0.5 + u - Power(u,2)/2. + Power(v,2);
					double d_uv_12 = 0.5 - u + Power(u,2)/2. - Power(v,2)/2.;
					
					double d_vv_1 = 0;
					double d_vv_2 = 0;
					double d_vv_3 = u + v - u*v - Power(v,2);
					double d_vv_4 = -2*u - 2*v + 2*u*v + 2*Power(v,2);
					double d_vv_5 = u + v - u*v - Power(v,2);
					double d_vv_6 = u*v + Power(v,2);
					double d_vv_7 = 1 - u - 2*v - u*v - Power(v,2);
					double d_vv_8 = -2 + 2*u + 4*v - u*v - Power(v,2);
					double d_vv_9 = 1 - u - 2*v + u*v + Power(v,2);
					double d_vv_10 = v - u*v - Power(v,2);
					double d_vv_11 = -2*v + 2*u*v + 2*Power(v,2);
					double d_vv_12 = v - u*v - Power(v,2);
		
					double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
					double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
					double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
					
					double ceff_map_duu[12] = { d_uu_8, d_uu_7, d_uu_4, d_uu_5, d_uu_9, d_uu_12, d_uu_11, d_uu_10, d_uu_6, d_uu_3, d_uu_1, d_uu_2 };
					double ceff_map_duv[12] = { d_uv_8, d_uv_7, d_uv_4, d_uv_5, d_uv_9, d_uv_12, d_uv_11, d_uv_10, d_uv_6, d_uv_3, d_uv_1, d_uv_2 };
					double ceff_map_dvv[12] = { d_vv_8, d_vv_7, d_vv_4, d_vv_5, d_vv_9, d_vv_12, d_vv_11, d_vv_10, d_vv_6, d_vv_3, d_vv_1, d_vv_2 };

					double rp[3] = {0,0,0};

					for( int x = 0; x < ncoords_base; x++ )
					{
						int *cset = theVertices[i].irr_coord_set + e * ncoords_base;

						for( int y = 0; y < 12; y++ )
						{
							rp[0] += theVertices[cset[x]].r[0] * theMap[y*ncoords_base+x] * ceff_map[y];
							rp[1] += theVertices[cset[x]].r[1] * theMap[y*ncoords_base+x] * ceff_map[y];
							rp[2] += theVertices[cset[x]].r[2] * theMap[y*ncoords_base+x] * ceff_map[y];
						}
					}
		//			printf("DBG C %lf %lf %lf\n", rp[0], rp[1], rp[2] ); 
				}
			
			}
		}
	}	

}

void surface::duplicate( surface *copy_surface, int dupex, int dupey, int dupez )
{
	cumulative_area = NULL;
	max_valence = 15; // points used.
	opencl_init = 0;
#ifdef FFTW
	h_in = NULL;
	h_out = NULL;
#endif
	bcs = NULL;
	nbc = 0;

	edgeFormulas = NULL;
	theFormulas = NULL;
	theVolumeFormulas = NULL;
	ptBoxes = NULL;
	prBoxes = NULL;
	nf_faces = 0;
	nf_g_q_p = 0;
	nf_irr_faces = 0;
	nf_irr_pts = 0;
	
	double OLA = copy_surface->PBC_vec[0][0];
	double OLB = copy_surface->PBC_vec[1][1];
	double OLC = copy_surface->PBC_vec[2][2];

	printf("Duplicating.\n");
	c0 = 0;
	PBC_vec[0][0] = copy_surface->PBC_vec[0][0] * dupex;
	PBC_vec[0][1] =  copy_surface->PBC_vec[0][1];
	PBC_vec[0][2] =  copy_surface->PBC_vec[0][2];
	PBC_vec[1][0] =  copy_surface->PBC_vec[1][0];
	PBC_vec[1][1] =  copy_surface->PBC_vec[1][1] * dupey;
	PBC_vec[1][2] =  copy_surface->PBC_vec[1][2];
	PBC_vec[2][0] =  copy_surface->PBC_vec[2][0];
	PBC_vec[2][1] =  copy_surface->PBC_vec[2][1];
	PBC_vec[2][2] =  copy_surface->PBC_vec[2][2] * dupez;

	PBC[0] = copy_surface->PBC[0]*dupex;
	PBC[1] = copy_surface->PBC[1]*dupey;
	PBC[2] = copy_surface->PBC[2]*dupez;

	double LA = PBC_vec[0][0];
	double LB = PBC_vec[1][1];
	double LC = PBC_vec[2][2];

	int nvo = copy_surface->nv;

	nv = copy_surface->nv * dupex * dupey * dupez;
	
	theVertices = (vertex *)malloc( sizeof(vertex) * nv );

	int nbytes = sizeof(vertex)*nv;

	char *init = (char *)theVertices;

	for( int c = 0; c < nbytes; c++ )
		init[c] = 0xFA;

	int v2 = 0;

	for( int vx = 0; vx < dupex; vx++ )
	for( int vy = 0; vy < dupey; vy++ )
	for( int vz = 0; vz < dupez; vz++ )
	{
		for( int v = 0; v < copy_surface->nv; v++ )
		{
			int val = copy_surface->theVertices[v].valence;
	
	
			theVertices[v2].r[0] = copy_surface->theVertices[v].r[0] + vx * OLA;
			theVertices[v2].r[1] = copy_surface->theVertices[v].r[1] + vy * OLB;
			theVertices[v2].r[2] = copy_surface->theVertices[v].r[2] + vz * OLC;


			for( int xx = 0; xx < copy_surface->theVertices[v].valence; xx++ )
			{
				// does this cross a periodic threshold?
				int j = copy_surface->theVertices[v].edges[xx];

				double *r0 = copy_surface->theVertices[v].r;
				double *r1 = copy_surface->theVertices[j].r;
		
				double dr[3] = { r1[0] - r0[0], r1[1] - r0[1], r1[2] - r0[2] };
				double add[3]={0,0,0};

				MinImage3D( dr, copy_surface->PBC_vec, add );

				int dcx = (int)lround(add[0]);
				int dcy = (int)lround(add[1]);
				int dcz = (int)lround(add[2]);

				int celli[3] = { vx, vy, vz };
				int cellj[3] = { vx + dcx, vy + dcy, vz + dcz };

				// no loop.
				int alt_vert = copy_surface->theVertices[v].edges[xx] + (vx * dupey*dupez + vy * dupez + vz) * copy_surface->nv;

				if( dupex > 0 && cellj[0] > celli[0] )
				{
					alt_vert += dupey*dupez * copy_surface->nv;
					if( cellj[0] >= dupex )
						alt_vert -= nv;
				}
				else if( dupex > 0 && cellj[0] < celli[0] )
				{
					alt_vert -= dupey*dupez * copy_surface->nv;
					if( cellj[0] < 0 )
						alt_vert += nv;
				}
				
				if( dupey > 0 && cellj[1] > celli[1] )
				{
					alt_vert += dupez * copy_surface->nv;
					
					if( cellj[1] >= dupey )
						alt_vert -= dupey * dupez * copy_surface->nv;
				}
				else if( dupey > 0 && cellj[1] < celli[1] )
				{
					alt_vert -= dupez * copy_surface->nv;
					if(cellj[1] < 0 )
						alt_vert += dupey * dupez * copy_surface->nv;
				}
				
				if( dupez > 0 && cellj[2] > celli[2] )
				{
					alt_vert += copy_surface->nv;
					if( cellj[2] >= dupez )
						alt_vert -= dupez * copy_surface->nv;
				}
				else if( dupez > 0 && cellj[2] < celli[2] )
				{
					alt_vert -= copy_surface->nv;
					if( cellj[2] < 0 )
						alt_vert += dupez * copy_surface->nv;
				}

				if( alt_vert < 0 || alt_vert >= nv)
				{
					printf("ERROR.\n");
					exit(1);
				}

				theVertices[v2].edges[xx] = alt_vert; 
			}
			theVertices[v2].valence = copy_surface->theVertices[v].valence;

			v2++;
		}
	}

	for( int i = 0; i < nv; i++ )
	{
		for( int e = 0; e < theVertices[i].valence; e++ )
		{
			int j = theVertices[i].edges[e];

			int gotit = 0;

			for( int e2 = 0; e2 < theVertices[j].valence; e2++ )
			{
				if( theVertices[j].edges[e2] == i )
				{
					theVertices[i].edge_rev[e] = e2;
					gotit += 1;
				}
			}

			if( !gotit || gotit > 1 )
			{
				printf("ERROR reversing edges for index %d.\n", i );
				exit(1);
			}
		}
	}

	
	theTriangles = NULL;
	theEdges = NULL;

	for( int i = 0; i < nv; i++ )
		theVertices[i].c0 = 0;
	constructTriangles();
	sortFaces();	
	// must be called AFTER sort faces.
	assignEdgePBC();	
}
void surface::generateBorderMappings( void )
{
	printf("Generating border mappings.\n");
	// determine the mapping for u and v across edges of the triangle.
	
	for( int t = 0; t < nt; t++ )
	{
		int i = theTriangles[t].ids[0];
		int j = theTriangles[t].ids[1];
		int k = theTriangles[t].ids[2];

		int edge_i_j = -1;

		for( int e = 0; e < theVertices[i].valence; e++ )
		{
			if( theVertices[i].edges[e] == j )
				edge_i_j = e;
		}

		int edge_i_k = edge_i_j+1;
		if( edge_i_k >= theVertices[i].valence )
			edge_i_k -= theVertices[i].valence;
		
		if( theVertices[i].edges[edge_i_k] != k )
		{
			int t = j;	
			j = k;
			k = t;
		
			for( int e = 0; e < theVertices[i].valence; e++ )
			{
				if( theVertices[i].edges[e] == j )
					edge_i_j = e;
			}

			edge_i_k = edge_i_j+1;
			if( edge_i_k >= theVertices[i].valence )
				edge_i_k -= theVertices[i].valence;
		}

 
		if( edge_i_j == -1 )		
		{
			printf("Logical error configuring face uv mapping.\n");
			exit(1);
		}
	
		if( i > j || i > k )
		{	
			printf("Logical error. ids[0] should always be lowest.\n");
			exit(1);
		}
	
		int border_tri[3];
	
		int edge_j_i = theVertices[i].edge_rev[edge_i_j];
		int edge_j_l = edge_j_i-2;

		if( edge_j_l < 0 )
			edge_j_l += theVertices[j].valence;
	
		int l = theVertices[j].edges[edge_j_l];

		int edge_i_m = edge_i_j-1;
		if( edge_i_m < 0 )
			edge_i_m += theVertices[i].valence;
		int m = theVertices[i].edges[edge_i_m];

		int edge_i_n = edge_i_j+2;
		if( edge_i_n >= theVertices[i].valence )
			edge_i_n -= theVertices[i].valence;

		int n = theVertices[i].edges[edge_i_n];
		
		theTriangles[t].edge_type = 0; // i is the lowest of all indices, j is the lowest of j, k, l


		if( m < i )
		{
			int edge_m_j = theVertices[i].edge_rev[edge_i_m]-1;
			if( edge_m_j < 0 ) edge_m_j += theVertices[m].valence;

			border_tri[0] = theVertices[m].faces[edge_m_j];	
			theTriangles[t].edge_type += MASK_1 & ( (~theTriangles[t].edge_type) & (1<<0));
		}
		else
			border_tri[0] = theVertices[i].faces[edge_i_m];	
			
		int edge_k_j = theVertices[i].edge_rev[edge_i_k]+1;
		if( edge_k_j >= theVertices[k].valence )
			edge_k_j -= theVertices[k].valence;

	
		if( j < k && j < l )
		{
			int edge_j_k = theVertices[k].edge_rev[edge_k_j];
			int edge_j_l = edge_j_k-1;
			if( edge_j_l < 0 )
				edge_j_l += theVertices[j].valence;

			border_tri[1] = theVertices[j].faces[edge_j_l];	
		}
		else if( k < j && k < l )
		{
			
			border_tri[1] = theVertices[k].faces[edge_k_j];

			theTriangles[t].edge_type += MASK_2 & ( (~theTriangles[t].edge_type) & (1<<1));
		}
		else
		{
			int edge_k_l = edge_k_j+1;
			if( edge_k_l >= theVertices[k].valence )
				edge_k_l -= theVertices[k].valence;

			int edge_l_k = theVertices[k].edge_rev[edge_k_l];
			
			border_tri[1] = theVertices[l].faces[edge_l_k];
			theTriangles[t].edge_type += MASK_2 & ( (~theTriangles[t].edge_type) & (1<<2));
		}
	

		
		if( n < i )
		{
			int edge_i_n = edge_i_k+1;
			if( edge_i_n >= theVertices[i].valence )
				edge_i_n -= theVertices[i].valence;

			int edge_n_i = theVertices[i].edge_rev[edge_i_n];

			theTriangles[t].edge_type += MASK_3 & ( (~theTriangles[t].edge_type) & (1<<3));
			border_tri[2] = theVertices[n].faces[edge_n_i];
		}
		else
			border_tri[2] = theVertices[i].faces[edge_i_k];

		theTriangles[t].border_tri[0] = border_tri[0];
		theTriangles[t].border_tri[1] = border_tri[1];
		theTriangles[t].border_tri[2] = border_tri[2];


	}
}








void surface::moveParticleonSurface(int *f, double *u, double *v, double *p) {
	double *r_surface = (double *)malloc( sizeof(double) * (3 * nv + 3) );
        get(r_surface);
        r_surface[3*nv+0] = 1.0;
        r_surface[3*nv+1] = 1.0;
        r_surface[3*nv+2] = 1.0;
        double sigma = .036; // the standard deviation of the distance distribution
        double frc_duv[2] = {0,0}; // no force on the particle;                              
        double fstep[3]={0,0,0}; // mainly for debugging.
        double dt = 1; // the time step only needed if we are doing forces.
        int bisection_iter = 10; // the max iterations of bisection we will do. Often the correct step can be guessed right off the bat.
//	printf("%lf %lf \n", *u, *v);
        localMove( &*f, &*u, &*v, sigma, r_surface, frc_duv, dt, fstep, bisection_iter);
	double nrm[3];
	evaluateRNRM(*f, *u, *v, p, nrm, r_surface);

}



void surface::checkCurvature(FILE* file) {
        double *r = (double *)malloc( sizeof(double) * (3 * nv + 3) );
        get(r);
        r[3*nv+0] = 1.0;
        r[3*nv+1] = 1.0;
        r[3*nv+2] = 1.0;
	double eps = 1e-10;
	int f = 0;
	int nf = f;
	double u = 0;
	double v = eps; 
	double dt = 0.01;
	double du = 1 * dt;
	double dv = 0.0 * dt;
	int check = 0;

	for (int i = 0; i < 20; i++) {
		double tempdu = du;
		double tempdv = dv;	
		double tempu = u;
		double tempv = v;
		printf("both are same %d %d \n", f, nf);
		while (nf == f) {
			u = tempu;
			v = tempv;
			du = tempdu;
			dv = tempdv;
			if (u > .999 || v == 1.0) {
				u = .99;
				printf("got in here %lf %lf %lf %lf \n", u, v, du, dv);
				nf = nextFace(f, &u, &v, &du, &dv, r);
				printf("next face is %d \n", nf);
				if (nf > nf_faces) {
					printf("IRREGULAR\n");
				}
			} else {
				printf("in face %d: %lf %lf\n", f, u, v);
				double k;
				double c = surface::c(f, u, v, r,&k);
	
				double nrm[3];
				double p[3];
				evaluateRNRM(f, u, v, p, nrm, r);
				fprintf(file, "%d %lf %lf %lf %lf\n", check, c, nrm[0], nrm[1], nrm[2]);
		
				check++;
				tempu += tempdu;
				tempv += tempdv;
				nf = nextFace(f, &u, &v, &du, &dv, r);
			}
		}
		f = nf;
		double mag = sqrt(pow(du, 2) + pow(dv, 2));
		du = du/mag * dt;
		dv = dv/mag * dt;
		printf("next face is %d %lf %lf %lf %lf\n", nf, u, v, du, dv);
	}

/**
        int f = 0;
	int nf = f;
	double u = 0;
	double v = 0;
	double du = 0.01;
	double dv = 0;
	int check = 0;
	for (int i = 0; i < 100; i++) {
		if (dv == 0) {
			while (u < .99) {
				printf("in face %d %lf %lf\n", f, u, v);
				double c = surface::c(f, u, v, r);
				fprintf(file, "%d %lf\n", check, c);
				check++;
				u += du;
				v += dv;
			}
		} else if (du == 0) {
                        while (v > 0.01) {
                                printf("in face %d %lf %lf\n", f, u, v);
                                double c = surface::c(f, u, v, r);
                                fprintf(file, "%d %lf\n", check, c);
                                check++;
                                u += du;
                                v += dv;
                        }
		}
		
		nf = nextFace(f, &u, &v, &du, &dv, r);
		printf("next face is %d %lf %lf %lf %lf\n", nf, u, v, du, dv);
		double mag = sqrt(pow(du, 2) + pow(dv, 2));
		du = du/mag * 0.01;
		dv = 0;
		if (u == 0 && du < 0 && v > 0) {
			dv = du;
			du = 0;

		}
		printf("next face is %d %lf %lf %lf %lf\n", nf, u, v, du, dv);
		f = nf;
	}
*/






}

void surface::wrapPBC( double *dr1, double *alphas)
{
	double put[3];
	MinImage3D( dr1, PBC_vec, put, alphas );
	
}

double surface::cellVolume( void )
{
	double cp[3];
	double a1[3] = { PBC_vec[0][0], PBC_vec[0][1], PBC_vec[0][2] };
	double a2[3] = { PBC_vec[1][0], PBC_vec[1][1], PBC_vec[1][2] };
	double a3[3] = { PBC_vec[2][0], PBC_vec[2][1], PBC_vec[2][2] };

	cross( a2, a3, cp );
	double vol = a1[0] * cp[0] + a1[1] * cp[1] + a1[2] * cp[2];

	return vol;
}

	
void surface::setEdgeRev( void )
{	
	// resorted edges need to be reversed.
	for( int i = 0; i < nv; i++ )
	{
		for( int e = 0; e < theVertices[i].valence; e++ )
		{
			int j = theVertices[i].edges[e];

//			printf("i: %d e: %d j: %d\n", i, e, j );
			for( int e2 = 0; e2 < theVertices[j].valence; e2++ )
			{
				if( theVertices[j].edges[e2] == i )
					theVertices[i].edge_rev[e] = e2;
			}
		}
	}
}
