#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include <sys/time.h>
#include "interp.h"
#include "mutil.h"
#include "l-bfgs.h"
#include "m_triangles.h"
#include "config.h"
#include "gsl/gsl_sf.h"
#include "gauss.h"
#include "lapack_we_use.h"
#include "Bin.h"
#include "Point.h"
#include "fast_mm.h"
#include "parallel.h"
#define SKIP_SUPPORT_CHECK
#define MEGA_DEL (10.0)

#define USE_MAX_LEVEL 10
#define COLLISION_LEVEL 10 

#define FRACTIONAL_TOLERANCE

#ifdef FRACTIONAL_TOLERANCE
#define FRAC_TOL      0.0025
#define EXTREME_R	(1e-6)
#else
//#define RADIUS_TOL    (1)
#define RADIUS_TOL    (1e-3)
#endif

bool gjk_algorithm( double *r1, int nv1, double *r2, int nv2 );
bool linear_collision_worker(double* r1, int nv1, double* pt1, double *pt2, double *pt,  double radius);

int tr_tri_intersect3D (double *C1, double *P1, double *P2,
	     double *D1, double *Q1, double *Q2);

void convexHullBounds( 
	double *r1, int nv1, double **M, int mlow, int mhigh, int level, int target_level, double min[3], double max[3] ); 
int nearMembraneWorker( double *r1, int nv1, double *r2, double **M, int mlow, int mhigh, int level, int max_level, double radius, 
		double factor, 
		double cur_u, double cur_v,
		double *puv_list, double *areas, int *npts );	

int cmpfunc (const void* a, const void* b) {
	double x = *(double*)a;
	double y = *(double*)b;
	if (fabs(x) < fabs(y)) return -1;
	if (fabs(x) > fabs(y)) return 1;
	return 0;
}

int cmpfunci (const void* a, const void* b) {
	int x = *(int*)a;
	int y = *(int*)b;
	if (fabs(x) < fabs(y)) return -1;
	if (fabs(x) > fabs(y)) return 1;

	return 0;
}

int checkCollision( double *r1, int nv1, double *r2, int nv2, double **M, int mlow, int mhigh, int level, int max_level )
{
	int pts[3] = { 0, 1, 2};

#if 0
	double *C1 = r1+pts[0]*3;
	double *D1 = r2+pts[0]*3;
	double P1[3] = { r1[pts[1]*3+0] - r1[pts[0]*3+0], r1[pts[1]*3+1] - r1[pts[0]*3+1], r1[pts[1]*3+2] - r1[pts[0]*3+2] };
	double P2[3] = { r1[pts[2]*3+0] - r1[pts[0]*3+0], r1[pts[2]*3+1] - r1[pts[0]*3+1], r1[pts[2]*3+2] - r1[pts[0]*3+2] };
	double Q1[3] = { r2[pts[1]*3+0] - r2[pts[0]*3+0], r2[pts[1]*3+1] - r2[pts[0]*3+1], r2[pts[1]*3+2] - r2[pts[0]*3+2] };
	double Q2[3] = { r2[pts[2]*3+0] - r2[pts[0]*3+0], r2[pts[2]*3+1] - r2[pts[0]*3+1], r2[pts[2]*3+2] - r2[pts[0]*3+2] };

	int tt = (tr_tri_intersect3D( C1, P1, P2, D1, Q1, Q2 ) != 0);
#endif
	int gjk= gjk_algorithm( r1, nv1, r2, nv2 );

	if( !gjk ) {
		return 0;
	}
		

	if( level == max_level )
	{

		return 1;
	}

	int nvmax = 6+mhigh;



	for( int x1 = 0; x1 < 4; x1++ )
	{
		double sub1[nvmax*3];
		memset( sub1, 0, sizeof(double)*nvmax*3);

		int nv_use1 = 12;
		if( x1 == 0 ) nv_use1 = nv1;
	
		double *mat1 = M[nv1-6-mlow]+x1*(nv1 < 12 ? 12 : nv1)*nv1;

		for( int i = 0; i < nv_use1; i++ )
		{
			for( int j = 0; j < nv1; j++ )
			{
				sub1[i*3+0] += mat1[i*nv1+j] * r1[3*j+0];	
				sub1[i*3+1] += mat1[i*nv1+j] * r1[3*j+1];	
				sub1[i*3+2] += mat1[i*nv1+j] * r1[3*j+2];	
			}
		}

		for( int x2 = 0; x2 < 4; x2++ )
		{
			double sub2[nvmax*3];
			memset( sub2, 0, sizeof(double)*nvmax*3);

			int nv_use2 = 12;
			if( x2 == 0 ) nv_use2 = nv2;
			double *mat2 = M[nv2-6-mlow]+x2*(nv2 < 12 ? 12 : nv2)*nv2;

			for( int i = 0; i < nv_use2; i++ )
			{
				for( int j = 0; j < nv2; j++ )
				{
					sub2[i*3+0] += mat2[i*nv2+j] * r2[3*j+0];	
					sub2[i*3+1] += mat2[i*nv2+j] * r2[3*j+1];	
					sub2[i*3+2] += mat2[i*nv2+j] * r2[3*j+2];	
				}
			}	

			if( checkCollision( sub1, nv_use1, sub2, nv_use2, M, mlow, mhigh, level+1, max_level ) )
				return 1;	
		}
	}

	return 0;
} 

int collision( surface *surface1, surface *surface2, double **M, int mlow, int mhigh)
{
	/* for now, brute force comparison of all triangular pairs */


	for( int t = 0; t < surface1->nt; t++ )
	{
		triangle *tri1 = surface1->theTriangles+t;

		// f is the index into my structures for computing properties.
		int f1 = tri1->f;

		int *indices1;
		int val1 = 0;
		int np1 = 0;
		int base1 = 0;
		double *pbc1;

		if( f1 >= surface1->nf_faces )
		{
			f1 -= surface1->nf_faces;
			int formulas_per_face = surface1->nf_irr_pts;
			indices1 = surface1->theIrregularFormulas[f1*formulas_per_face].cp;
			pbc1 = surface1->theIrregularFormulas[f1*formulas_per_face].r_pbc;
			np1 = surface1->theIrregularFormulas[f1*formulas_per_face].ncoor;
			base1 = surface1->theIrregularFormulas[f1*formulas_per_face].vertex;
		}
		else
		{
			int formulas_per_face = surface1->nf_g_q_p;
			indices1 = surface1->theFormulas[f1*formulas_per_face].cp;
			pbc1 = surface1->theFormulas[f1*formulas_per_face].r_pbc;
			np1 = surface1->theFormulas[f1*formulas_per_face].ncoor;
			base1 = surface1->theFormulas[f1*formulas_per_face].vertex;
		}			

		double pts1[3*np1];
		for( int x = 0; x < np1; x++ )
		{
			pts1[3*x+0] = surface1->theVertices[indices1[x]].r[0] + pbc1[3*x+0]; 
			pts1[3*x+1] = surface1->theVertices[indices1[x]].r[1] + pbc1[3*x+1]; 
			pts1[3*x+2] = surface1->theVertices[indices1[x]].r[2] + pbc1[3*x+2];
		} 

		for( int t2 = 0; t2 < surface2->nt; t2++ )
		{
			triangle *tri2 = surface2->theTriangles+t2;

			// f is the index into my structures for computing properties.
			int f2 = tri2->f;
	
			int *indices2;
			int val2 = 0;
			int np2 = 0;
			int base2 = 0;
			double *pbc2;
	
			if( f2 >= surface2->nf_faces )
			{
				f2 -= surface2->nf_faces;
				int formulas_per_face = surface2->nf_irr_pts;
				indices2 = surface2->theIrregularFormulas[f2*formulas_per_face].cp;
				pbc2 = surface2->theIrregularFormulas[f2*formulas_per_face].r_pbc;
				np2 = surface2->theIrregularFormulas[f2*formulas_per_face].ncoor;
				base2 = surface2->theIrregularFormulas[f2*formulas_per_face].vertex;
			}
			else
			{
				int formulas_per_face = surface2->nf_g_q_p;
				indices2 = surface2->theFormulas[f2*formulas_per_face].cp;
				pbc2 = surface2->theFormulas[f2*formulas_per_face].r_pbc;
				np2 = surface2->theFormulas[f2*formulas_per_face].ncoor;
				base2 = surface2->theFormulas[f2*formulas_per_face].vertex;
			}			

			double pts2[3*np2];
			for( int x = 0; x < np2; x++ )
			{
				pts2[3*x+0] = surface2->theVertices[indices2[x]].r[0] + pbc2[3*x+0]; 
				pts2[3*x+1] = surface2->theVertices[indices2[x]].r[1] + pbc2[3*x+1]; 
				pts2[3*x+2] = surface2->theVertices[indices2[x]].r[2] + pbc2[3*x+2];
			} 
						
			if(  checkCollision( pts1, np1, pts2, np2, M, mlow, mhigh, 0, USE_MAX_LEVEL ) )
				return 1;
		}
		 
	}

	return 0;	
}


bool minimum_distance( double *r1, int nv1, double *r2, double radius, int do_quick_abort );
int checkSphereHullCollision( double *r1, int nv1, double *r2, double **M, int mlow, int mhigh, int level, int max_level, double radius, 
		double factor, 
		double cur_u, double cur_v,
		double *col_u, 	
		double *col_v )
{
	int pts[3] = { 0, 1, 2};

#if 0
	double *C1 = r1+pts[0]*3;
	double *D1 = r2+pts[0]*3;
	double P1[3] = { r1[pts[1]*3+0] - r1[pts[0]*3+0], r1[pts[1]*3+1] - r1[pts[0]*3+1], r1[pts[1]*3+2] - r1[pts[0]*3+2] };
	double P2[3] = { r1[pts[2]*3+0] - r1[pts[0]*3+0], r1[pts[2]*3+1] - r1[pts[0]*3+1], r1[pts[2]*3+2] - r1[pts[0]*3+2] };
	double Q1[3] = { r2[pts[1]*3+0] - r2[pts[0]*3+0], r2[pts[1]*3+1] - r2[pts[0]*3+1], r2[pts[1]*3+2] - r2[pts[0]*3+2] };
	double Q2[3] = { r2[pts[2]*3+0] - r2[pts[0]*3+0], r2[pts[2]*3+1] - r2[pts[0]*3+1], r2[pts[2]*3+2] - r2[pts[0]*3+2] };

	int tt = (tr_tri_intersect3D( C1, P1, P2, D1, Q1, Q2 ) != 0);
#endif
	bool min = minimum_distance( r1, nv1, r2, radius, level!=max_level );


	if( !min ) {
		return 0;
	}
		

	if( level == max_level )
	{
//		exit(1);
//	int min = minimum_distance( r1, nv1, r2, radius );
//	printf("here is it: %d", radius);
//		FILE *theFile = fopen("collision.xyz","w");
//		fprintf(theFile, "%d\n", nv1);
//		fprintf(theFile, "Collision points\n");
//		for( int x = 0; x < nv1; x++ )
//			fprintf(theFile, "C %.16le %.16le %.16le\n", r1[3*x+0], r1[3*x+1], r1[3*x+2] );
//		fclose(theFile);

	//	printf("Found collision at %le %le\n", cur_u, cur_v );
		*col_u = cur_u + factor * (1.0/3.0);
		*col_v = cur_v + factor * (1.0/3.0);

		return 1;
	}

	int nvmax = 6+mhigh;



	for( int x1 = 0; x1 < 4; x1++ )
	{
		double sub1[nvmax*3];
		memset( sub1, 0, sizeof(double)*nvmax*3);

		int nv_use1 = 12;
		if( x1 == 0 ) nv_use1 = nv1;
	
		double *mat1 = M[nv1-6-mlow]+x1*(nv1 < 12 ? 12 : nv1)*nv1;

		if( nv1 != 12 )
		{
//			printf("Using matrix %d (%d vertices).\n", nv1-6-mlow, nv1 );
		}
		double one = 1.0, zero = 0.0;
		int incrxy=3;

		double check1[nvmax*3];
		memset(check1, 0, sizeof(double) * 3 * nvmax );
		char trans = 'T';
		int nI = nv1;
		int nO = nv_use1;
#define USE_DGEMV
#ifdef USE_DGEMV
		if(  nv1 == 12  )
		{
			if( x1 == 0 )
				fast_sub6_0( r1, sub1 );
			else if( x1 == 1 )
				fast_sub6_1( r1, sub1 );
			else if( x1 == 2 )
				fast_sub6_2( r1, sub1 );
			else if( x1 == 3 )
				fast_sub6_3( r1, sub1 );
		}
		else
		{
			for( int c = 0; c < 3; c++ )
				dgemv( &trans, &nI, &nO, &one, mat1, &nv1, r1+c, &incrxy, &zero, sub1+c, &incrxy );
		}
#else
		for( int i = 0; i < nv_use1; i++ )
		{
			for( int j = 0; j < nv1; j++ )
			{
				sub1[i*3+0] += mat1[i*nv1+j] * r1[3*j+0];	
				sub1[i*3+1] += mat1[i*nv1+j] * r1[3*j+1];	
				sub1[i*3+2] += mat1[i*nv1+j] * r1[3*j+2];	
			}
		}
#endif

		double new_u = cur_u;
		double new_v = cur_v;
		
		double new_factor = factor;

		switch(x1)
		{
			case 0:
				new_factor *= 0.5;
				break;
			case 1:
				new_u += factor * 0.5;
				new_factor *= 0.5;
				break;
			case 2:
				new_v += factor * 0.5;
				new_factor *= 0.5;
				break;
			case 3:
				new_u += factor * 0.5;
				new_v += factor * 0.5;
				new_factor *= -0.5;
				break;
		
		}

			if( checkSphereHullCollision( sub1, nv_use1, r2, M, mlow, mhigh, level+1, max_level, radius,
				new_factor, new_u, new_v, col_u, col_v ) )
				return 1;	
	}

	return 0;
} 
double near_point_on_triangle( double v1x, double v1y, double v1z,
				double v2x, double v2y, double v2z,
				double v3x, double v3y, double v3z, int *type, int *which );

int checkLinearCollision( double *r1, int nv1, double *pt1, double *pt2, double *pt, double **M, int mlow, int mhigh, int level, int max_level, double radius, 
		double factor, 
		double cur_u, double cur_v,
		double *col_u, 	
		double *col_v )
{
	bool min = linear_collision_worker( r1, nv1, pt1, pt2, pt, radius );

//#define LINEAR_DEBUG

#ifdef LINEAR_DEBUG
	double triangle[9] = { r1[0], r1[1], r1[2],
			       r1[3], r1[4], r1[5],
			       r1[6], r1[7], r1[8] };
	
	double vec[3] = { pt2[0] - pt1[0],
			  pt2[1] - pt1[1],
			  pt2[2] - pt1[2] };

	double v1[3] = { triangle[0] - triangle[3],
			 triangle[1] - triangle[4],
			 triangle[2] - triangle[5] };	
	double v2[3] = { triangle[6] - triangle[3],
			 triangle[7] - triangle[4],
			 triangle[8] - triangle[5] };	
	double tnrm[3];
	cross( v1, v2, tnrm );
	normalize(tnrm);

	double t1 = pt1[0] * tnrm[0] + pt1[1] * tnrm[1] + pt1[2] * tnrm[2]; 
	double t2 = pt2[0] * tnrm[0] + pt2[1] * tnrm[1] + pt2[2] * tnrm[2]; 
	double k  = r1[0] * tnrm[0] + r1[1] * tnrm[1] + r1[2] * tnrm[2]; 

	

	double planar_point[3] = { pt1[0] + (k-t1)/(t2-t1) * (pt2[0]-pt1[0]),
			           pt1[1] + (k-t1)/(t2-t1) * (pt2[1]-pt1[1]),
			           pt1[2] + (k-t1)/(t2-t1) * (pt2[2]-pt1[2]) };
	int type,which;
	near_point_on_triangle( triangle[0]-planar_point[0],triangle[1]-planar_point[1],triangle[2]-planar_point[2],
				triangle[3]-planar_point[0],triangle[4]-planar_point[1],triangle[5]-planar_point[2],
				triangle[6]-planar_point[0],triangle[7]-planar_point[1],triangle[8]-planar_point[2],
					&type, &which );			
	

	if( (t1 < k && t2 < k) || ( t1 > k && t2 > k ) )
		type = 1;

	printf("line-tri collision: %d min: %d\n", (type == 0 ? 1 : 0), (min ? 1 : 0 ) );
	if( type == 0 && !min && level > 4 )
	{
		printf("GJK error I think.\n");
		bool min = linear_collision_worker( r1, nv1, pt1, pt2, pt, radius );
	near_point_on_triangle( triangle[0]-planar_point[0],triangle[1]-planar_point[1],triangle[2]-planar_point[2],
				triangle[3]-planar_point[0],triangle[4]-planar_point[1],triangle[5]-planar_point[2],
				triangle[6]-planar_point[0],triangle[7]-planar_point[1],triangle[8]-planar_point[2],
					&type, &which );			
		exit(1);
	}	
#endif
	if( !min ) {
		return 0;
	}
		

	if( level == max_level )
	{
		*col_u = cur_u + factor * (1.0/3.0);
		*col_v = cur_v + factor * (1.0/3.0);

		return 1;
	}

	int nvmax = 6+mhigh;

	for( int x1 = 0; x1 < 4; x1++ )
	{
		double sub1[nvmax*3];
		memset( sub1, 0, sizeof(double)*nvmax*3);

		int nv_use1 = 12;
		if( x1 == 0 ) nv_use1 = nv1;
	
		double *mat1 = M[nv1-6-mlow]+x1*(nv1 < 12 ? 12 : nv1)*nv1;

		if( nv1 != 12 )
		{
//			printf("Using matrix %d (%d vertices).\n", nv1-6-mlow, nv1 );
		}
		double one = 1.0, zero = 0.0;
		int incrxy=3;

		double check1[nvmax*3];
		memset(check1, 0, sizeof(double) * 3 * nvmax );
		char trans = 'T';
		int nI = nv1;
		int nO = nv_use1;
#define USE_DGEMV
#ifdef USE_DGEMV
		if(  nv1 == 12  )
		{
			if( x1 == 0 )
				fast_sub6_0( r1, sub1 );
			else if( x1 == 1 )
				fast_sub6_1( r1, sub1 );
			else if( x1 == 2 )
				fast_sub6_2( r1, sub1 );
			else if( x1 == 3 )
				fast_sub6_3( r1, sub1 );
		}
		else
		{
			for( int c = 0; c < 3; c++ )
				dgemv( &trans, &nI, &nO, &one, mat1, &nv1, r1+c, &incrxy, &zero, sub1+c, &incrxy );
		}
#else
		for( int i = 0; i < nv_use1; i++ )
		{
			for( int j = 0; j < nv1; j++ )
			{
				sub1[i*3+0] += mat1[i*nv1+j] * r1[3*j+0];	
				sub1[i*3+1] += mat1[i*nv1+j] * r1[3*j+1];	
				sub1[i*3+2] += mat1[i*nv1+j] * r1[3*j+2];	
			}
		}
#endif

		double new_u = cur_u;
		double new_v = cur_v;
		
		double new_factor = factor;

		switch(x1)
		{
			case 0:
				new_factor *= 0.5;
				break;
			case 1:
				new_u += factor * 0.5;
				new_factor *= 0.5;
				break;
			case 2:
				new_v += factor * 0.5;
				new_factor *= 0.5;
				break;
			case 3:
				new_u += factor * 0.5;
				new_v += factor * 0.5;
				new_factor *= -0.5;
				break;
		
		}

			if( checkLinearCollision( sub1, nv_use1, pt1, pt2, pt, M, mlow, mhigh, level+1, max_level, radius,
				new_factor, new_u, new_v, col_u, col_v ) )
				return 1;	
	}

	return 0;
} 

double surface::returnRadius (double *pt, int *col_f, double *col_u, double *col_v, double **M, int mlow, int mhigh, double distance, double currentR, double L, double defaultR, double *vertex_data, int *ptr_to_data, int *nump, int inside_outside, int disable_PBC_z, double reflecting_surface ) {

	// distance: the distance the particle moved since we last checked.	
	// currentR: the last particle/membrane distance
	// defaultR: the largest radius we care about.

	double check_fudge = 0.6;
	double alpha = 0.75;
	// checkR: currentR + alpha*distance, where we are *likely* to find the particle
	double checkR = currentR + alpha*distance;
	// there is no point in looking farther than defaultR since that's the largest distance we care about.
	if (checkR >= defaultR) {
		checkR = defaultR;
	}
	
	// if it moved farther than the current distance to the membrane, this is the distance from the membrane it could be.
	double checkClose = distance - currentR + check_fudge + reflecting_surface;

	if( checkClose < reflecting_surface + check_fudge )
		checkClose = reflecting_surface + check_fudge;

	if (checkR <= checkClose) {
		if (!withinRadius(pt, col_f, col_u, col_v, M, mlow, mhigh, L, checkClose, vertex_data, ptr_to_data, nump, disable_PBC_z )) {
			// it's not on the other side and we did a quick check to see that.
			return checkClose;
		} else {
			// if we moved far enough we need to check to see if we crossed to the other side.
			int f;
			double u, v, radius;
			// this checks to see which side we are on: checks the particles displacement from the surface relative to the surface normal.
			bool check = withinBoxedSurface( pt, &f, &u, &v, M, mlow, mhigh, &radius, L, inside_outside, disable_PBC_z, reflecting_surface );

			if (!check) {
				// not within the surface, outside it.
				*col_f = f;
				*col_u = u;
				*col_v = v;
				return -1;
			} else {
//				printf("withinBoxedSurface returns true, radius %le\n", radius);
				return radius;
			}
		}
	} else {
		bool did_collide = 1;
		double lowerR = 0;
		double upperR = checkR;
		double trial_R = upperR;
		while (did_collide) {
			if ((trial_R  > checkClose) && !withinRadius(pt, col_f, col_u, col_v, M, mlow, mhigh, L, trial_R, vertex_data, ptr_to_data, nump, disable_PBC_z )) {
				if (trial_R > defaultR) { 
//					printf("Trial R > defaultR.\n");
					return defaultR;
				}
//				printf("last distance %le, distance moved: %le trial_R: %le\n", currentR, distance, trial_R ); 
//				printf("Trial_R: %le\n", trial_R );
                        	return trial_R;
			} else if ((trial_R > checkClose) && withinRadius(pt, col_f, col_u, col_v, M, mlow, mhigh, L, trial_R, vertex_data, ptr_to_data, nump, disable_PBC_z )) {
				//printf("check2\n");
				if( trial_R < reflecting_surface )
					return -1;

				upperR = trial_R;
				trial_R = (lowerR + upperR)/2;
			} else if ((trial_R < checkClose) && !withinRadius(pt, col_f, col_u, col_v, M, mlow, mhigh, L, trial_R, vertex_data, ptr_to_data, nump, disable_PBC_z )) {
				//printf("check3\n");
				if (trial_R > defaultR) {
					printf("Trial_R > defaultR\n");
					return defaultR;
				}
				lowerR = trial_R;
				trial_R = (lowerR + upperR)/2;

			} else if ((trial_R < checkClose) && withinRadius(pt, col_f, col_u, col_v, M, mlow, mhigh, L, trial_R, vertex_data, ptr_to_data, nump, disable_PBC_z )) {
				
				if( trial_R < reflecting_surface )
					return -1;
				//printf("check4\n");
	        		int f;
				double u, v, radius;
		                bool check = withinBoxedSurface( pt, &f, &u, &v, M, mlow, mhigh, &radius, L, inside_outside, disable_PBC_z, reflecting_surface );
        	                if (!check) {
//					printf("Verifying, within boxed surface collision.\n");
                	                return -1;
                       		} else {
//					printf("Not verified, returning radius %le.\n", radius );
                                	return radius;
                        	}
			}

		}
	}


	return -1;
}
bool surface::withinRadius( double *pt_in, int *col_f, double *col_u, double *col_v, double **M, int mlow, int mhigh, double L, double radius,	
		double *vertex_data, int *ptr_to_data, int *nump, int disable_PBC_z  )
{
	 //for now, brute force comparison of all triangular pairs 

	double Lx = L;
	double Ly = L;
	double Lz = L;

	if( L< 0 )
	{
		Lx = PBC_vec[0][0];	
		Ly = PBC_vec[1][1];	
		Lz = PBC_vec[2][2];	
	}

	double pt[3] = { pt_in[0], pt_in[1], pt_in[2] };

	while( pt[0] < -Lx/2 ) pt[0] += Lx;
	while( pt[1] < -Ly/2 ) pt[1] += Ly;

	if( !disable_PBC_z )
		while( pt[2] < -Lz/2 ) pt[2] += Lz;
	
	while( pt[0] > Lx/2 ) pt[0] -= Lx;
	while( pt[1] > Ly/2 ) pt[1] -= Ly;

	if( !disable_PBC_z )
		while( pt[2] > Lz/2 ) pt[2] -= Lz;




	double trial_R = radius;

                                                struct timeval tp;
                                                gettimeofday(&tp, NULL);
                                                double time_1 = tp.tv_usec * (1e-6) + tp.tv_sec;
		// find boxes to search inside

		bool did_collide = 0;
	
		double tpt[3] = { pt[0] + Lx/2, pt[1] + Ly/2, pt[2] + Lz/2 };
	
		while( tpt[0] < 0 ) tpt[0] += Lx;
		while( tpt[1] < 0 ) tpt[1] += Ly;
		
		if( !disable_PBC_z )
			while( tpt[2] < 0 ) tpt[2] += Lz;
		
		while( tpt[0] >= Lx ) tpt[0] -= Lx;
		while( tpt[1] >= Ly ) tpt[1] -= Ly;

		if( !disable_PBC_z )
			while( tpt[2] >= Lz ) tpt[2] -= Lz;
			
		int cenx = xbox * (tpt[0]) / Lx; 
		int ceny = ybox * (tpt[1]) / Ly; 
		int cenz = zbox * (tpt[2]) / Lz;

		int low_bin_x = xbox * (tpt[0] - trial_R - max_width_x) / Lx;
		int low_bin_y = ybox * (tpt[1] - trial_R - max_width_y) / Ly;
		int low_bin_z = zbox * (tpt[2] - trial_R - max_width_z) / Lz;
		
		int high_bin_x = xbox * (tpt[0] + trial_R + max_width_x) / Lx;
		int high_bin_y = ybox * (tpt[1] + trial_R + max_width_y) / Ly;
		int high_bin_z = zbox * (tpt[2] + trial_R + max_width_z) / Lz;
		

		if( high_bin_x - low_bin_x >= xbox )
			high_bin_x = low_bin_x + xbox - 1;
		if( high_bin_y - low_bin_y >= ybox )
			high_bin_y = low_bin_y + ybox - 1;
		if( high_bin_z - low_bin_z >= zbox )
			high_bin_z = low_bin_z + zbox - 1;
	
		while( high_bin_x < 0 )
		{
			low_bin_x += xbox;
			high_bin_x += xbox;
		}
		while( high_bin_y < 0 )
		{
			low_bin_y += ybox;
			high_bin_y += ybox;
		}
		while( high_bin_z < 0 )
		{
			low_bin_z += zbox;
			high_bin_z += zbox;
		}

		int arrayX[xbox];
		int arrayY[ybox];
		int arrayZ[zbox];
		int sizeX = 0;
		int sizeY = 0;
		int sizeZ = 0;
		
		for( int t_xB = low_bin_x; t_xB <= high_bin_x; t_xB++, sizeX++ ) 
			arrayX[sizeX] = t_xB - cenx;
		qsort(arrayX, sizeX, sizeof(int), cmpfunci);
		for( int t_yB = low_bin_y; t_yB <= high_bin_y; t_yB++, sizeY++ ) 
			arrayY[sizeY] = t_yB - ceny;
		qsort(arrayY, sizeY, sizeof(int), cmpfunci);
		for( int t_zB = low_bin_z; t_zB <= high_bin_z; t_zB++, sizeZ++ ) 
			arrayZ[sizeZ] = t_zB - cenz;
		qsort(arrayZ, sizeZ, sizeof(int), cmpfunci);

	
/*	
		for( int t_xB = - binloop_x; t_xB <= binloop_x + fudge_x; t_xB++, sizeX++ ) 
			arrayX[sizeX] = t_xB;
		qsort(arrayX, sizeX, sizeof(int), cmpfunci);
		for( int t_yB = - binloop_y; t_yB <= binloop_y + fudge_y; t_yB++, sizeY++ ) 
			arrayY[sizeY] = t_yB;
		qsort(arrayY, sizeY, sizeof(int), cmpfunci);
		for( int t_zB = - binloop_z; t_zB <= binloop_z + fudge_z; t_zB++, sizeZ++ ) 
			arrayZ[sizeZ] = t_zB;
		qsort(arrayZ, sizeZ, sizeof(int), cmpfunci);
*/
		int numChecks = 0;
		int numTris = 0;
		for( int t_xB = 0; t_xB < sizeX && !did_collide; t_xB++ ) 
		{
			int xB =  cenx + arrayX[t_xB];

			if( xB < 0 ) xB += xbox;
			if( xB >= xbox ) xB -= xbox;

			for( int t_yB = 0; t_yB < sizeY && !did_collide; t_yB++ ) 
			{ 
				int yB =  ceny + arrayY[t_yB];
				
				if( yB < 0 ) yB += ybox;
				if( yB >= ybox ) yB -= ybox;

				for( int t_zB = 0; t_zB < sizeZ && !did_collide; t_zB++ )
				{
					int zB = cenz + arrayZ[t_zB];
				
					if( zB < 0 ) zB += zbox;
					if( zB >= zbox ) zB -= zbox;

					bool worthCheck = true;		
#ifndef SKIP_SUPPORT_CHECK
					if (box[xB][yB][zB].num_tri > 0) {
						double support = pt[0]*box[xB][yB][zB].search[0] + pt[1]*box[xB][yB][zB].search[1] + pt[2]*box[xB][yB][zB].search[2];
						double test = support + trial_R;
						double test2 = support - trial_R;
						double mindiff = box[xB][yB][zB].vmin;
						double maxdiff = box[xB][yB][zB].vmax;

						if ((test < mindiff && test2 < mindiff) || (test > maxdiff && test2 > maxdiff)) {
							worthCheck = false;
						} else {
							numChecks++;
						}
					}
#endif

					for (int ind = 0; ind < box[xB][yB][zB].num_tri && !did_collide && worthCheck; ind++) {
						numTris++;
						triangle *tri1 = theTriangles + box[xB][yB][zB].tri_list[ind];
                        			int f1 = tri1->f;
/*
						int *indices1;
                      				int val1 = 0;
                        			int np1 = 0;
                        			int base1 = 0;
						double *pbc1;
	
                        			if( f1 >= nf_faces )
                        			{

                                			f1 -= nf_faces;
                                			int formulas_per_face = nf_irr_pts;
                                			indices1 = theIrregularFormulas[f1*formulas_per_face].cp;
							pbc1 = theIrregularFormulas[f1*formulas_per_face].r_pbc;
                                			np1 = theIrregularFormulas[f1*formulas_per_face].ncoor;
                                			base1 = theIrregularFormulas[f1*formulas_per_face].vertex;
                        			}
                        			else
                        			{
                                			int formulas_per_face = nf_g_q_p;
                                			indices1 = theFormulas[f1*formulas_per_face].cp;
							pbc1 = theFormulas[f1*formulas_per_face].r_pbc;
                                			np1 = theFormulas[f1*formulas_per_face].ncoor;
                                			base1 = theFormulas[f1*formulas_per_face].vertex;
                        			}

                        			double pts1[3*np1];
                        			for( int x = 0; x < np1; x++ )
                        			{
                                			pts1[3*x+0] = theVertices[indices1[x]].r[0] + pbc1[3*x+0];
                                			pts1[3*x+1] = theVertices[indices1[x]].r[1] + pbc1[3*x+1];
                                			pts1[3*x+2] = theVertices[indices1[x]].r[2] + pbc1[3*x+2];
							
                        			}
					*/	

						int np1 = nump[tri1->f];
						double *pts1 = vertex_data + ptr_to_data[tri1->f]; 

						double usep[3] = { pt[0], pt[1], pt[2] };

						while( usep[0] - pts1[0] < -Lx/2 ) usep[0] += Lx;
						while( usep[1] - pts1[1] < -Ly/2 ) usep[1] += Ly;
						if( !disable_PBC_z )
							while( usep[2] - pts1[2] < -Lz/2 ) usep[2] += Lz;

						while( usep[0] - pts1[0] > Lx/2 ) usep[0] -= Lx;
						while( usep[1] - pts1[1] > Ly/2 ) usep[1] -= Ly;
						if( !disable_PBC_z )
							while( usep[2] - pts1[2] > Lz/2 ) usep[2] -= Lz;

                        			if(  checkSphereHullCollision( pts1, np1, usep, M, mlow, mhigh, 0, USE_MAX_LEVEL, trial_R,
                                			1.0, 1e-6, 1e-6, col_u, col_v ) )
                        			{
                                			*col_f =tri1->f;
                                			did_collide = 1;
		        			}
	
					} 
                		}
			}
		}
	return did_collide;
}

void surface::returnSizeofBox(int* x, int* y, int* z) {
	*x = xbox;
	*y = ybox; 
	*z = zbox;
}

void surface::box_system( double edge_length ) {
   
        double *r_surface = (double *)malloc( sizeof(double) * (3 * nv + 3) );
        get(r_surface);
	
	double av_edge_length = 0;
	double n_edge_length = 0;

	for ( int x = 0; x < nv; x++ )
	{
		int val = theVertices[x].valence;

		for( int e = 0; e < val; e++ )
		{
			int j = theVertices[x].edges[e];

			double *epbc = theVertices[x].edge_PBC+3*e;

			double dr[3] = { 
				r_surface[3*x+0] - r_surface[3*j+0] - (epbc[0]*PBC_vec[0][0] + epbc[1] * PBC_vec[1][0] + epbc[2]*PBC_vec[2][0]), 
				r_surface[3*x+1] - r_surface[3*j+1] - (epbc[0]*PBC_vec[0][1] + epbc[1] * PBC_vec[1][1] + epbc[2]*PBC_vec[2][1]), 
				r_surface[3*x+2] - r_surface[3*j+2] - (epbc[0]*PBC_vec[0][2] + epbc[1] * PBC_vec[1][2] + epbc[2]*PBC_vec[2][2]) };
			double r = normalize(dr);

			av_edge_length += r;
			n_edge_length += 1;
		}
	}

	av_edge_length /= n_edge_length;

	if( edge_length < 0 )
		edge_length = av_edge_length;

	double LA;
	double LB;
	double LC;

	LA = PBC_vec[0][0];
	LB = PBC_vec[1][1];
	LC = PBC_vec[2][2];
	

        double max_x=0, max_y=0, max_z=0;
	Point* averageTris = new Point[nt];
	
	double *M5 = (double *)malloc( sizeof(double) * 4 * 11 * 12 ); 
	double *M6 = (double *)malloc( sizeof(double) * 4 * 12 * 12 ); 
	double *M7 = (double *)malloc( sizeof(double) * 4 * 13 * 13 );

	int mlow = 5;
	int mhigh = 7;

	double *M[3] = { M5, M6, M7 };

	generateSubdivisionMatrices( M, mlow, mhigh );

                for( int t = 0; t < nt; t++ )
                {
                        triangle *tri1 = theTriangles+t;

                        int f1 = tri1->f;

                        int *indices1;
                        int val1 = 0;
                        int np1 = 0;
                        int base1 = 0;
			double *pbc1;
                        if( f1 >= nf_faces )
                        {
                                f1 -= nf_faces;
                                int formulas_per_face = nf_irr_pts;
                                indices1 = theIrregularFormulas[f1*formulas_per_face].cp;
				pbc1 = theIrregularFormulas[f1*formulas_per_face].r_pbc;
                                np1 = theIrregularFormulas[f1*formulas_per_face].ncoor;
                                base1 = theIrregularFormulas[f1*formulas_per_face].vertex;
                        }
                        else
                        {
                                int formulas_per_face = nf_g_q_p;
                                indices1 = theFormulas[f1*formulas_per_face].cp;
				pbc1 = theFormulas[f1*formulas_per_face].r_pbc;
                                np1 = theFormulas[f1*formulas_per_face].ncoor;
                                base1 = theFormulas[f1*formulas_per_face].vertex;
                        }

                        double pts1[3*np1];
                        for( int x = 0; x < np1; x++ )
                        {
                                pts1[3*x+0] = theVertices[indices1[x]].r[0] + pbc1[3*x+0];
                                pts1[3*x+1] = theVertices[indices1[x]].r[1] + pbc1[3*x+1];
                                pts1[3*x+2] = theVertices[indices1[x]].r[2] + pbc1[3*x+2];
                        }

                        Point hull_min = {pts1[0], pts1[1], pts1[2]};
                        Point hull_max = {pts1[0], pts1[1], pts1[2]};
#define BETTER_BOUNDS

#ifdef BETTER_BOUNDS
			double min[3] = { 1e10, 1e10, 1e10 };
			double max[3] = { -1e10, -1e10, -1e10 };

			convexHullBounds( pts1, np1, M, mlow, mhigh, 0, 3, min, max );

			hull_min.x = min[0];
			hull_min.y = min[1];
			hull_min.z = min[2];
			
			hull_max.x = max[0];
			hull_max.y = max[1];
			hull_max.z = max[2];
#else
                        for (int j = 1; j < np1; j++) {
                                if (pts1[3*j+0] > hull_max.x) {
                                        hull_max.x = pts1[3*j+0];
                                } else if (pts1[3*j+0] < hull_min.x) {
                                        hull_min.x = pts1[3*j+0];
                                }

                                if (pts1[3*j+1] > hull_max.y) {
                                        hull_max.y = pts1[3*j+1];
                                } else if (pts1[3*j+1] < hull_min.y) {
                                        hull_min.y = pts1[3*j+1];
                                }
                                
                                if (pts1[3*j+2] > hull_max.z) {
                                        hull_max.z = pts1[3*j+2];
                                } else if (pts1[3*j+2] < hull_min.z) {
                                        hull_min.z = pts1[3*j+2];
                                }
                        }
#endif

                        Point diff = hull_max - hull_min;
                        if (diff.x > max_x) {
                                max_x = diff.x;
                        }
                        if (diff.y > max_y) {
                                max_y = diff.y;
                        }
                        if (diff.z > max_z) {
                                max_z = diff.z;
                        }

			Point average = {(pts1[0] + pts1[3] + pts1[6])/3, (pts1[1] + pts1[4] + pts1[7])/3, (pts1[2] + pts1[5] + pts1[8])/3};
//			printf("avereage point is %d %d %d\n", average.x, average.y, average.z);
			averageTris[t] = average;
                }

	free(M5);
	free(M6);
	free(M7);

        // target about 20 x 20 x 20?

	double factor = 1;

	xlength = edge_length;
	ylength = edge_length;
	zlength = edge_length;

        xbox = factor * LA / xlength;
	ybox = factor * LB / ylength;
	zbox = factor * LC / zlength; 

//	ylength = LB / ybox;
//	zlength = LC / zbox;
	 
/*   
          xlength = 50;
          ylength = 50;
          zlength = 50;
*/
	double rmax = max_x;
	if( max_y > max_x )
		rmax = max_y;
	if( max_z > rmax )
		rmax = max_z;
	max_width_x = rmax; 
	max_width_y = rmax; 
	max_width_z = rmax; 
     

	printf("max_width_xyz: %le %le %le\n", max_width_x, max_width_y, max_width_z );

  /*     
        double txbox = LA/xlength;
        double tybox = LB/ylength;
        double tzbox = LC/zlength;

	xbox = (int)txbox;
	ybox = (int)tybox;
	zbox = (int)tzbox;
*/
//	printf("dims of box are %lf %lf %lf\n", max_x, max_y, max_z);
	printf("dims of box are %d %d %d\n", xbox, ybox, zbox);

        box = new Bin**[xbox];
        for (int i = 0; i < xbox; i++) {
                box[i] = new Bin*[ybox];
                for (int j = 0; j < ybox; j++) {
                        box[i][j] = new Bin[zbox];
			for (int k = 0; k < zbox; k++) {
				box[i][j][k].num_tri = 0;
				box[i][j][k].space = 10;
		//		box[i][j][k].tri_list = new int[box[i][j][k].space];
				box[i][j][k].tri_list =(int*) malloc(sizeof(int) * box[i][j][k].space);
				for (int size = 0; size < 10; size++) {
					box[i][j][k].tri_list[size] = 0;
				}
			}
                }
        }

	for (int i = 0; i < nt; i++) {

		double dx = averageTris[i].x + LA/2;
		while( dx < 0 ) dx += LA;
		while( dx >= LA ) dx -= LA;
		double dy = averageTris[i].y + LB/2;
		while( dy < 0 ) dy += LB;
		while( dy >= LB ) dy -= LB;
		double dz = averageTris[i].z + LC/2;
		while( dz < 0 ) dz += LC;
		while( dz >= LC ) dz -= LC;


		int x = (int)((dx) * xbox)/LA;
		int y = (int)((dy) * ybox)/LB;
		int z = (int)((dz) * zbox)/LC;

		

		if (box[x][y][z].num_tri == box[x][y][z].space) {
			box[x][y][z].space += 10;
			box[x][y][z].tri_list = (int*) realloc(box[x][y][z].tri_list, sizeof(int) * box[x][y][z].space);
			//delete[] box[x][y][z].tri_list;
			//box[x][y][z].tri_list = new int[box[x][y][z].space];
		}

		box[x][y][z].tri_list[box[x][y][z].num_tri] = i;
		box[x][y][z].num_tri++;

	}


	for (int i = 0; i < xbox; i++) {
		for (int j = 0; j < ybox; j++) {
			for (int k = 0; k < zbox; k++) {
//				printf("inside boxes %d\n", box[i][j][k].num_tri);
				double totNorm[3];
				totNorm[0] = 0;
				totNorm[1] = 0;
				totNorm[2] = 0;
				double vmin = 1e10;
				double vmax = -1e10;
				for (int ind = 0; ind < box[i][j][k].num_tri; ind++) {
					triangle *tri1 = theTriangles+box[i][j][k].tri_list[ind];
	                        	int f1 = tri1->f;
			               	int u = 0;
					int v = 0; 
					double rp[3];
			                double nrm[3];
                			r_surface[3*nv+0] = 1.0;
                			r_surface[3*nv+1] = 1.0;
                			r_surface[3*nv+2] = 1.0;
                			evaluateRNRM(f1, 1./3., 1./3., rp, nrm, r_surface );
				//printf("triangle normal is %lf %lf %lf\n", nrm[0], nrm[1], nrm[2]);
					for (int a = 0; a < 3; a++) {
						totNorm[a] = totNorm[a] + nrm[a];
					}
				}

				double magnitude = sqrt(pow(totNorm[0], 2) + pow(totNorm[1], 2) + pow(totNorm[2], 2));
				for (int a = 0; a < 3; a++) {
					totNorm[a] = totNorm[a]/magnitude;
				}					
				//printf("%lf %lf %lf\n", totNorm[0], totNorm[1], totNorm[2]);
				for (int ind = 0; ind < box[i][j][k].num_tri; ind++) {
                                        triangle *tri1 = theTriangles+box[i][j][k].tri_list[ind];
                                        int f1 = tri1->f;
        	                	int *indices1;
                	        	int val1 = 0;
                       			int np1 = 0;
	                        	int base1 = 0;
					double *pbc1;
		                        if( f1 >= nf_faces )
        		                {
                		                f1 -= nf_faces;
                        		        int formulas_per_face = nf_irr_pts;
                               			indices1 = theIrregularFormulas[f1*formulas_per_face].cp;
						pbc1 = theIrregularFormulas[f1*formulas_per_face].r_pbc;
	                                	np1 = theIrregularFormulas[f1*formulas_per_face].ncoor;
	        	                        base1 = theIrregularFormulas[f1*formulas_per_face].vertex;
        	        	        }
                	        	else
	                	        {
        	                	        int formulas_per_face = nf_g_q_p;
                	                	indices1 = theFormulas[f1*formulas_per_face].cp;
				pbc1 = theFormulas[f1*formulas_per_face].r_pbc;
	                        	        np1 = theFormulas[f1*formulas_per_face].ncoor;
        	                        	base1 = theFormulas[f1*formulas_per_face].vertex;
                	        	}
	
		                        double pts1[3*np1];
        		                for( int x = 0; x < np1; x++ )
                		        {
                        		        pts1[3*x+0] = theVertices[indices1[x]].r[0] + pbc1[3*x+0];
                                		pts1[3*x+1] = theVertices[indices1[x]].r[1] + pbc1[3*x+1];
                                		pts1[3*x+2] = theVertices[indices1[x]].r[2] + pbc1[3*x+2];
						
						double ptsupport = pts1[3*x+0] * totNorm[0] + pts1[3*x+1] * totNorm[1] + pts1[3*x+2] * totNorm[2];
//						printf("pt support is: %lf\n", ptsupport);
						if (ptsupport > vmax) {
							vmax = ptsupport;
						}
						if (ptsupport < vmin) {
							vmin = ptsupport;
						}

                        		}
				}
			//	box[i][j][k].search = new double[3];
				box[i][j][k].search[0] = totNorm[0];
				box[i][j][k].search[1] = totNorm[1];
				box[i][j][k].search[2] = totNorm[2];
				box[i][j][k].vmin = vmin;
				box[i][j][k].vmax = vmax;
			}
		}
	}
					
	free(r_surface);

	// print indeces in boxes
	for (int i = 0; i < xbox; i++) {
		for (int j = 0; j < ybox; j++) {
			for (int k = 0; k < zbox; k++) {
//				printf("%d %d %d: ", i, j, k);
				for (int z = 0; z < box[i][j][k].num_tri; z++) {

//				printf("%d ", box[i][j][k].tri_list[z]);

				}
				//printf("%lf %lf %lf", box[i][j][k].search[0], box[i][j][k].search[1], box[i][j][k].search[2]);
			//	printf("min: %lf max: %lf\n", box[i][j][k].vmin, box[i][j][k].vmax);
			//	printf("%lf ", sqrt(pow(box[i][j][k].search[0], 2) + pow(box[i][j][k].search[1], 2) + pow(box[i][j][k].search[2], 2)));
//				printf("\n");
			}
		}
	}


	delete [] averageTris;
}
	
bool surface::withinBoxedSurface(double* pt, int *f, double *u, double *v, double **M, int mlow, int mhigh, double *distance, double L, int inside_outside, int disable_PBC_z, double reflecting_surface ) {

 	if( inside_outside == 0 )
		return false;

                nearPointOnBoxedSurface( pt, f, u, v, M, mlow, mhigh, distance, disable_PBC_z);
                //nearPointOnSurface( pt, f, u, v, M5, M6, distance );

                double rp[3];
                double nrm[3];
                double *r_surface = (double *)malloc( sizeof(double) * (3 * nv + 3) );
                get(r_surface);
                r_surface[3*nv+0] = 1.0;
                r_surface[3*nv+1] = 1.0;
                r_surface[3*nv+2] = 1.0;
                //printf("col_f %d u: %le v: %le\n", *col_f, *col_u, *col_v );
                evaluateRNRM( *f, *u, *v, rp, nrm, r_surface );

                double pt_to_surface[3] = { pt[0] - rp[0], pt[1] - rp[1], pt[2] - rp[2] };
	
		if( !disable_PBC_z )
			wrapPBC( pt_to_surface, r_surface + 3 * nv  );

                //double rl = normalize(pt_to_surface);
                free(r_surface);
                //find dot product of normal and vector to point

		double dist = normalize(pt_to_surface);
//		printf("%le %le %le %d DOT NRM %le dist %le\n", pt[0], pt[1], pt[2], *f, pt_to_surface[0] * nrm[0] + pt_to_surface[1] * nrm[1] + pt_to_surface[2] * nrm[2], dist );
		if ((pt_to_surface[0] * nrm[0] + pt_to_surface[1] * nrm[1] + pt_to_surface[2] * nrm[2])*(inside_outside) > 0) {
			//not in surface

			
			return false;
		} else {
		//	printf("%lf %lf %lf, %lf %lf %lf\n", pt_to_surface[0], pt_to_surface[1], pt_to_surface[2], nrm[0], nrm[1], nrm[2]);

			if( dist > reflecting_surface )
			{
				if( *distance < reflecting_surface )
					return false;

				return true;
			}
			else
				return false;
		}
}



void surface::rebox_system( void )
{
	double LA = PBC_vec[0][0];
	double LB = PBC_vec[1][1];
	double LC = PBC_vec[2][2];

        double max_x=0, max_y=0, max_z=0;
	Point* averageTris = new Point[nt];
	
	double *M5 = (double *)malloc( sizeof(double) * 4 * 11 * 12 ); 
	double *M6 = (double *)malloc( sizeof(double) * 4 * 12 * 12 ); 
	double *M7 = (double *)malloc( sizeof(double) * 4 * 13 * 13 );

	int mlow = 5;
	int mhigh = 7;

	double *M[3] = { M5, M6, M7 };

	generateSubdivisionMatrices( M, mlow, mhigh );

                for( int t = 0; t < nt; t++ )
                {
                        triangle *tri1 = theTriangles+t;

                        int f1 = tri1->f;

                        int *indices1;
                        int val1 = 0;
                        int np1 = 0;
                        int base1 = 0;
			double *pbc1;
                        if( f1 >= nf_faces )
                        {
                                f1 -= nf_faces;
                                int formulas_per_face = nf_irr_pts;
                                indices1 = theIrregularFormulas[f1*formulas_per_face].cp;
				pbc1 = theIrregularFormulas[f1*formulas_per_face].r_pbc;
                                np1 = theIrregularFormulas[f1*formulas_per_face].ncoor;
                                base1 = theIrregularFormulas[f1*formulas_per_face].vertex;
                        }
                        else
                        {
                                int formulas_per_face = nf_g_q_p;
                                indices1 = theFormulas[f1*formulas_per_face].cp;
				pbc1 = theFormulas[f1*formulas_per_face].r_pbc;
                                np1 = theFormulas[f1*formulas_per_face].ncoor;
                                base1 = theFormulas[f1*formulas_per_face].vertex;
                        }

                        double pts1[3*np1];
                        for( int x = 0; x < np1; x++ )
                        {
                                pts1[3*x+0] = theVertices[indices1[x]].r[0] + pbc1[3*x+0];
                                pts1[3*x+1] = theVertices[indices1[x]].r[1] + pbc1[3*x+1];
                                pts1[3*x+2] = theVertices[indices1[x]].r[2] + pbc1[3*x+2];
                        }

                        Point hull_min = {pts1[0], pts1[1], pts1[2]};
                        Point hull_max = {pts1[0], pts1[1], pts1[2]};
#define BETTER_BOUNDS

#ifdef BETTER_BOUNDS
			double min[3] = { 1e10, 1e10, 1e10 };
			double max[3] = { -1e10, -1e10, -1e10 };

			convexHullBounds( pts1, np1, M, mlow, mhigh, 0, 1, min, max );

			hull_min.x = min[0];
			hull_min.y = min[1];
			hull_min.z = min[2];
			
			hull_max.x = max[0];
			hull_max.y = max[1];
			hull_max.z = max[2];
#else
                        for (int j = 1; j < np1; j++) {
                                if (pts1[3*j+0] > hull_max.x) {
                                        hull_max.x = pts1[3*j+0];
                                } else if (pts1[3*j+0] < hull_min.x) {
                                        hull_min.x = pts1[3*j+0];
                                }

                                if (pts1[3*j+1] > hull_max.y) {
                                        hull_max.y = pts1[3*j+1];
                                } else if (pts1[3*j+1] < hull_min.y) {
                                        hull_min.y = pts1[3*j+1];
                                }
                                
                                if (pts1[3*j+2] > hull_max.z) {
                                        hull_max.z = pts1[3*j+2];
                                } else if (pts1[3*j+2] < hull_min.z) {
                                        hull_min.z = pts1[3*j+2];
                                }
                        }
#endif

                        Point diff = hull_max - hull_min;
                        if (diff.x > max_x) {
                                max_x = diff.x;
                        }
                        if (diff.y > max_y) {
                                max_y = diff.y;
                        }
                        if (diff.z > max_z) {
                                max_z = diff.z;
                        }

			Point average = {(pts1[0] + pts1[3] + pts1[6])/3, (pts1[1] + pts1[4] + pts1[7])/3, (pts1[2] + pts1[5] + pts1[8])/3};
//			printf("avereage point is %d %d %d\n", average.x, average.y, average.z);
			averageTris[t] = average;
                }

	free(M5);
	free(M6);
	free(M7);

	double rmax = max_x;
	if( max_y > max_x )
		rmax = max_y;
	if( max_z > rmax )
		rmax = max_z;
	max_width_x = rmax; 
	max_width_y = rmax; 
	max_width_z = rmax; 
    
        for (int i = 0; i < xbox; i++) {
                for (int j = 0; j < ybox; j++) {
			for (int k = 0; k < zbox; k++) {
				box[i][j][k].num_tri = 0;
				for (int size = 0; size < box[i][j][k].space; size++) {
					box[i][j][k].tri_list[size] = 0;
				}
			}
                }
        }

	for (int i = 0; i < nt; i++) {

		double dx = averageTris[i].x + LA/2;
		while( dx < 0 ) dx += LA;
		while( dx >= LA ) dx -= LA;
		double dy = averageTris[i].y + LB/2;
		while( dy < 0 ) dy += LB;
		while( dy >= LB ) dy -= LB;
		double dz = averageTris[i].z + LC/2;
		while( dz < 0 ) dz += LC;
		while( dz >= LC ) dz -= LC;


		int x = (int)((dx) * xbox)/LA;
		int y = (int)((dy) * ybox)/LB;
		int z = (int)((dz) * zbox)/LC;

		

		if (box[x][y][z].num_tri == box[x][y][z].space) {
			box[x][y][z].space += 10;
			box[x][y][z].tri_list = (int*) realloc(box[x][y][z].tri_list, sizeof(int) * box[x][y][z].space);
			//delete[] box[x][y][z].tri_list;
			//box[x][y][z].tri_list = new int[box[x][y][z].space];
		}

		box[x][y][z].tri_list[box[x][y][z].num_tri] = i;
		box[x][y][z].num_tri++;

	}

        double *r_surface = (double *)malloc( sizeof(double) * (3 * nv + 3) );
        get(r_surface);

	for (int i = 0; i < xbox; i++) {
		for (int j = 0; j < ybox; j++) {
			for (int k = 0; k < zbox; k++) {
//				printf("inside boxes %d\n", box[i][j][k].num_tri);
				double totNorm[3];
				totNorm[0] = 0;
				totNorm[1] = 0;
				totNorm[2] = 0;
				double vmin = 1e10;
				double vmax = -1e10;
				for (int ind = 0; ind < box[i][j][k].num_tri; ind++) {
					triangle *tri1 = theTriangles+box[i][j][k].tri_list[ind];
	                        	int f1 = tri1->f;
			               	int u = 0;
					int v = 0; 
					double rp[3];
			                double nrm[3];
                			r_surface[3*nv+0] = 1.0;
                			r_surface[3*nv+1] = 1.0;
                			r_surface[3*nv+2] = 1.0;
                			evaluateRNRM(f1, 1./3., 1./3., rp, nrm, r_surface );
				//printf("triangle normal is %lf %lf %lf\n", nrm[0], nrm[1], nrm[2]);
					for (int a = 0; a < 3; a++) {
						totNorm[a] = totNorm[a] + nrm[a];
					}
				}

				double magnitude = sqrt(pow(totNorm[0], 2) + pow(totNorm[1], 2) + pow(totNorm[2], 2));
				for (int a = 0; a < 3; a++) {
					totNorm[a] = totNorm[a]/magnitude;
				}					
				//printf("%lf %lf %lf\n", totNorm[0], totNorm[1], totNorm[2]);
				for (int ind = 0; ind < box[i][j][k].num_tri; ind++) {
                                        triangle *tri1 = theTriangles+box[i][j][k].tri_list[ind];
                                        int f1 = tri1->f;
        	                	int *indices1;
                	        	int val1 = 0;
                       			int np1 = 0;
	                        	int base1 = 0;
					double *pbc1;
		                        if( f1 >= nf_faces )
        		                {
                		                f1 -= nf_faces;
                        		        int formulas_per_face = nf_irr_pts;
                               			indices1 = theIrregularFormulas[f1*formulas_per_face].cp;
						pbc1 = theIrregularFormulas[f1*formulas_per_face].r_pbc;
	                                	np1 = theIrregularFormulas[f1*formulas_per_face].ncoor;
	        	                        base1 = theIrregularFormulas[f1*formulas_per_face].vertex;
        	        	        }
                	        	else
	                	        {
        	                	        int formulas_per_face = nf_g_q_p;
                	                	indices1 = theFormulas[f1*formulas_per_face].cp;
				pbc1 = theFormulas[f1*formulas_per_face].r_pbc;
	                        	        np1 = theFormulas[f1*formulas_per_face].ncoor;
        	                        	base1 = theFormulas[f1*formulas_per_face].vertex;
                	        	}
	
		                        double pts1[3*np1];
        		                for( int x = 0; x < np1; x++ )
                		        {
                        		        pts1[3*x+0] = theVertices[indices1[x]].r[0] + pbc1[3*x+0];
                                		pts1[3*x+1] = theVertices[indices1[x]].r[1] + pbc1[3*x+1];
                                		pts1[3*x+2] = theVertices[indices1[x]].r[2] + pbc1[3*x+2];
						
						double ptsupport = pts1[3*x+0] * totNorm[0] + pts1[3*x+1] * totNorm[1] + pts1[3*x+2] * totNorm[2];
//						printf("pt support is: %lf\n", ptsupport);
						if (ptsupport > vmax) {
							vmax = ptsupport;
						}
						if (ptsupport < vmin) {
							vmin = ptsupport;
						}

                        		}
				}
			//	box[i][j][k].search = new double[3];
				box[i][j][k].search[0] = totNorm[0];
				box[i][j][k].search[1] = totNorm[1];
				box[i][j][k].search[2] = totNorm[2];
				box[i][j][k].vmin = vmin;
				box[i][j][k].vmax = vmax;
			}
		}
	}
					
	free(r_surface);

	// print indeces in boxes
	for (int i = 0; i < xbox; i++) {
		for (int j = 0; j < ybox; j++) {
			for (int k = 0; k < zbox; k++) {
//				printf("%d %d %d: ", i, j, k);
				for (int z = 0; z < box[i][j][k].num_tri; z++) {

//				printf("%d ", box[i][j][k].tri_list[z]);

				}
				//printf("%lf %lf %lf", box[i][j][k].search[0], box[i][j][k].search[1], box[i][j][k].search[2]);
			//	printf("min: %lf max: %lf\n", box[i][j][k].vmin, box[i][j][k].vmax);
			//	printf("%lf ", sqrt(pow(box[i][j][k].search[0], 2) + pow(box[i][j][k].search[1], 2) + pow(box[i][j][k].search[2], 2)));
//				printf("\n");
			}
		}
	}

	delete [] averageTris;

}

void convexHullBounds( 
	double *r1, int nv1, double **M, int mlow, int mhigh, int level, int target_level, double min[3], double max[3] ) 
{
	if( level == target_level )
	{
		for( int x = 0; x < nv1; x++ )
		{
			if( r1[3*x+0] < min[0] ) min[0] = r1[3*x+0];
			if( r1[3*x+1] < min[1] ) min[1] = r1[3*x+1];
			if( r1[3*x+2] < min[2] ) min[2] = r1[3*x+2];
			if( r1[3*x+0] > max[0] ) max[0] = r1[3*x+0];
			if( r1[3*x+1] > max[1] ) max[1] = r1[3*x+1];
			if( r1[3*x+2] > max[2] ) max[2] = r1[3*x+2];
		}

		return;
	}

	int nvmax = 6+mhigh;

	for( int x1 = 0; x1 < 4; x1++ )
	{
		double sub1[nvmax*3];
		memset( sub1, 0, sizeof(double)*nvmax*3);

		int nv_use1 = 12;
		if( x1 == 0 ) nv_use1 = nv1;
	
		double *mat1 = M[nv1-6-mlow]+x1*(nv1 < 12 ? 12 : nv1)*nv1;

		double one = 1.0, zero = 0.0;
		int incrxy=3;

		double check1[nvmax*3];
		memset(check1, 0, sizeof(double) * 3 * nvmax );
		char trans = 'T';
		int nI = nv1;
		int nO = nv_use1;
#define USE_DGEMV
#ifdef USE_DGEMV
		for( int c = 0; c < 3; c++ )
			dgemv( &trans, &nI, &nO, &one, mat1, &nv1, r1+c, &incrxy, &zero, sub1+c, &incrxy );
#else
		for( int i = 0; i < nv_use1; i++ )
		{
			for( int j = 0; j < nv1; j++ )
			{
				sub1[i*3+0] += mat1[i*nv1+j] * r1[3*j+0];	
				sub1[i*3+1] += mat1[i*nv1+j] * r1[3*j+1];	
				sub1[i*3+2] += mat1[i*nv1+j] * r1[3*j+2];	
			}
		}
#endif

		convexHullBounds( sub1, nv_use1, M, mlow, mhigh, level+1, target_level, min, max ); 
	}
} 

void surface::nearPointOnBoxedSurface( double *pt_in, int *col_f, double *col_u, double *col_v, double **M, int mlow, int mhigh, double *distance, double initial_distance, int disable_PBC_z)
{
	double pt[3] = { pt_in[0], pt_in[1], pt_in[2] };
	double Lx = PBC_vec[0][0];
	double Ly = PBC_vec[1][1];
	double Lz = PBC_vec[2][2];

/*	
	we do not need to loop back and check these boxes. 

*/
	int nbc = 0;
	int nbtot = xbox*ybox*zbox;
	int *box_checked_list = (int *)malloc( sizeof(int) * nbtot );
	char *box_blank = (char *)malloc( sizeof(char) * nbtot);
	memset( box_blank, 0, sizeof(char) * nbtot );
	int *uncheck_these = (int *)malloc( sizeof(int) * nbtot);
	int nuncheck = 0;



	double lower_R = 0;

	double *aPt = theVertices[0].r;

	// quick estimate of upper R...
	double upper_R = 500;
	if( initial_distance > 0 )
		upper_R = initial_distance;

	int bisection_done = 0;

#ifndef FRACTIONAL_TOLERANCE
	double TOL_BISECTION = RADIUS_TOL;
#endif
	int check = 0;
	bool initial = false;

	double nskipped = 0;
	double nnot_skipped = 0;

	while( pt[0] < -Lx/2 ) pt[0] += Lx;
	while( pt[0] >= Lx/2 ) pt[0] -= Lx;
	if( !disable_PBC_z )
		while( pt[1] < -Ly/2 ) pt[1] += Ly;

	while( pt[1] >= Ly/2 ) pt[1] -= Ly;
	while( pt[2] < -Lz/2 ) pt[2] += Lz;
	if( !disable_PBC_z )
		while( pt[2] >= Lz/2 ) pt[2] -= Lz;
	

	while( !bisection_done )
	{
		double trial_R = (lower_R+upper_R)/2;
//		printf("trial R is %lf\n", trial_R);
		if (initial == false) {
			trial_R = upper_R;
		}

                struct timeval tp;
                gettimeofday(&tp, NULL);
                double time_1 = tp.tv_usec * (1e-6) + tp.tv_sec;

		double tfudge = 0;

		pt[0] += Lx/2;
		pt[1] += Ly/2;
		pt[2] += Lz/2;

		double minx = pt[0] - trial_R - max_width_x - tfudge;
		double maxx = pt[0] + trial_R + max_width_x + tfudge;
		
		double miny = pt[1] - trial_R - max_width_x - tfudge;
		double maxy = pt[1] + trial_R + max_width_x + tfudge;
		
		double minz = pt[2] - trial_R - max_width_x - tfudge;
		double maxz = pt[2] + trial_R + max_width_x + tfudge;
		
#ifdef TOO_CONSERVATIVE
		double minx = pt[0] - trial_R - max_width_x - xlength - tfudge;
		double maxx = pt[0] + trial_R + max_width_x + xlength + tfudge;
		
		double miny = pt[1] - trial_R - max_width_x - xlength - tfudge;
		double maxy = pt[1] + trial_R + max_width_x + xlength + tfudge;
		
		double minz = pt[2] - trial_R - max_width_x - xlength - tfudge;
		double maxz = pt[2] + trial_R + max_width_x + xlength + tfudge;
#endif
		int cenx = xbox * pt[0] / Lx; 
		int ceny = ybox * pt[1] / Ly; 
		int cenz = zbox * pt[2] / Lz;
		
		pt[0] -= Lx/2;
		pt[1] -= Ly/2;
		pt[2] -= Lz/2;

		int binloop_x = xbox * (trial_R + max_width_x + xlength) / Lx;
		int binloop_y = ybox * (trial_R + max_width_y + ylength) / Ly;
		int binloop_z = zbox * (trial_R + max_width_z + zlength) / Lz;
	
		if( binloop_x >= xbox/2 ) binloop_x = xbox/2;
		if( binloop_y >= ybox/2 ) binloop_y = ybox/2;
		if( binloop_z >= zbox/2 ) binloop_z = zbox/2;
	
		int fudge_x = 0;
		int fudge_y = 0;
		int fudge_z = 0;

		if( binloop_x * 2 == xbox )
			fudge_x = -1;
		if( binloop_y * 2 == ybox )
			fudge_y = -1;
		if( binloop_z * 2 == zbox )
			fudge_z = -1;
	

					int did_collide = 0;
					int numChecks = 0;
					int numTris = 0;
//					printf("pulse\n");
			nbc = 0;
		
		int arrayX[xbox];
		int arrayY[ybox];
		int arrayZ[zbox];
		int sizeX = 0;
		int sizeY = 0;
		int sizeZ = 0;
		
		for( int t_xB = - binloop_x; t_xB <= binloop_x + fudge_x; t_xB++, sizeX++ ) 
			arrayX[sizeX] = t_xB;
		qsort(arrayX, sizeX, sizeof(int), cmpfunci);
		for( int t_yB = - binloop_y; t_yB <= binloop_y + fudge_y; t_yB++, sizeY++ ) 
			arrayY[sizeY] = t_yB;
		qsort(arrayY, sizeY, sizeof(int), cmpfunci);
		for( int t_zB = - binloop_z; t_zB <= binloop_z + fudge_z; t_zB++, sizeZ++ ) 
			arrayZ[sizeZ] = t_zB;
		qsort(arrayZ, sizeZ, sizeof(int), cmpfunci);

		for( int t_xB = 0; t_xB < sizeX && !did_collide; t_xB++ ) 
		{
			int xB =  cenx + arrayX[t_xB];

			if( xB < 0 ) xB += xbox;
			if( xB >= xbox ) xB -= xbox;

			for( int t_yB = 0; t_yB < sizeY && !did_collide; t_yB++ ) 
			{ 
				int yB =  ceny + arrayY[t_yB];
				
				if( yB < 0 ) yB += ybox;
				if( yB >= ybox ) yB -= ybox;

				for( int t_zB = 0; t_zB < sizeZ && !did_collide; t_zB++ )
				{
					int zB =  cenz + arrayZ[t_zB];
				
					if( zB < 0 ) zB += zbox;
					if( zB >= zbox ) zB -= zbox;

					bool worthCheck = true;	

					int box_id = (xB*ybox+yB)*zbox+zB;


					if( box_blank[box_id] ) 
					{
						nskipped++;
						continue;
					}
					else
					{
						nnot_skipped++;
					}

					if( nbc < nbtot )
					{
						box_checked_list[nbc] = box_id;
						nbc++;
					}
					else
					{
						printf("Looping over boxes twice?\n");
						exit(1);
					}
#ifndef SKIP_SUPPORT_CHECK	
					if (box[xB][yB][zB].num_tri > 0) {
						double support = pt[0]*box[xB][yB][zB].search[0] + pt[1]*box[xB][yB][zB].search[1] + pt[2]*box[xB][yB][zB].search[2];
						//double mindiff = abs(support - box[xB][yB][zB].vmin);
						//double maxdiff = abs(support - box[xB][yB][zB].vmax);

						double test = support + trial_R;
						double test2 = support - trial_R;
						double mindiff = box[xB][yB][zB].vmin;
						double maxdiff = box[xB][yB][zB].vmax;

/*						printf("support %le %le %le min %le max %le\n",
						 box[xB][yB][zB].search[0],
						 box[xB][yB][zB].search[1],
						 box[xB][yB][zB].search[2], mindiff, maxdiff );
*/
						if ((test < mindiff && test2 < mindiff) || (test > maxdiff && test2 > maxdiff)) {
							worthCheck = false;
						} else {
							numChecks++;
						}
					}
#endif

					for (int ind = 0; ind < box[xB][yB][zB].num_tri && !did_collide && worthCheck; ind++) {
						numTris++;
						//printf("cmae to check\n");
						//printf("index: %d\n", box[xB][yB][zB].tri_list[ind]);
						triangle *tri1 = theTriangles + box[xB][yB][zB].tri_list[ind];
                        			int f1 = tri1->f;
			//			printf("box: %d %d %d;  %d %lf\n", xB, yB, zB, box[xB][yB][zB].tri_list[ind], trial_R);                        			

						int *indices1;
                      				int val1 = 0;
                        			int np1 = 0;
                        			int base1 = 0;
						double *pbc1;
			//			printf("here1\n");

                        			if( f1 >= nf_faces )
                        			{
			//			printf("here2\n");

                                			f1 -= nf_faces;
                                			int formulas_per_face = nf_irr_pts;
                                			indices1 = theIrregularFormulas[f1*formulas_per_face].cp;
							pbc1 = theIrregularFormulas[f1*formulas_per_face].r_pbc;
                                			np1 = theIrregularFormulas[f1*formulas_per_face].ncoor;
                                			base1 = theIrregularFormulas[f1*formulas_per_face].vertex;
		//				printf("here3\n");
                        			}
                        			else
                        			{
		//				printf("here4\n");
                                			int formulas_per_face = nf_g_q_p;
                                			indices1 = theFormulas[f1*formulas_per_face].cp;
							pbc1 = theFormulas[f1*formulas_per_face].r_pbc;
                                			np1 = theFormulas[f1*formulas_per_face].ncoor;
                                			base1 = theFormulas[f1*formulas_per_face].vertex;
		//				printf("here5\n");
                        			}

                        			double pts1[3*np1];
                        			for( int x = 0; x < np1; x++ )
                        			{
                                			pts1[3*x+0] = theVertices[indices1[x]].r[0] + pbc1[3*x+0];
                                			pts1[3*x+1] = theVertices[indices1[x]].r[1] + pbc1[3*x+1];
                                			pts1[3*x+2] = theVertices[indices1[x]].r[2] + pbc1[3*x+2];
							
//							printf("point %d: %lf %lf %lf\n", x, pts1[3*x+0], pts1[3*x+1], pts1[3*x+2]);
                        			}
					
			//			printf("make it here\n");
						
						
						double usep[3] = { pt[0], pt[1], pt[2] };

						while( usep[0] - pts1[0] < -Lx/2 ) usep[0] += Lx;
						while( usep[1] - pts1[1] < -Ly/2 ) usep[1] += Ly;
						if( !disable_PBC_z ) 
							while( usep[2] - pts1[2] < -Lz/2 ) usep[2] += Lz;
						while( usep[0] - pts1[0] > Lx/2 ) usep[0] -= Lx;
						while( usep[1] - pts1[1] > Ly/2 ) usep[1] -= Ly;
						if( !disable_PBC_z ) 
							while( usep[2] - pts1[2] > Lz/2 ) usep[2] -= Lz;

                        			if(  checkSphereHullCollision( pts1, np1, usep, M, mlow, mhigh, 0, USE_MAX_LEVEL, trial_R,
                                			1.0, 1e-6, 1e-6, col_u, col_v ) )
                        			{
                                			*col_f =tri1->f;
                                			initial = true;
                                			did_collide = 1;

							for( int p = 0; p < nbc-1; p++ )
							{
								box_blank[box_checked_list[p]] = 1;
								uncheck_these[nuncheck] = box_checked_list[p];
								nuncheck++;
							}				
			

//							printf("new collided at %lf: %lf %lf\n", trial_R, lower_R, upper_R);
                		//			printf("COLLLLLLLLLIDDDDDEEEEEEEE %lf\n", trial_R);
		        			}
					} 
                		}
			}
		}
//		printf("WHATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
                                                gettimeofday(&tp, NULL);
                                                double time_2 = tp.tv_usec * (1e-6) + tp.tv_sec;
         //                                       printf("Radius %lf took %d checks, %d triangles. %le sec \n", trial_R, numChecks, numTris, time_2-time_1);
       		if (initial == false) {
       			lower_R = upper_R;
       			upper_R = upper_R * 2;
       			trial_R = (lower_R + upper_R)/2;
       		} else {
       			if( did_collide )
               			upper_R = trial_R;
       			else
               			lower_R = trial_R;
       		}

#ifdef FRACTIONAL_TOLERANCE
		double dr_f = 2*(upper_R-lower_R)/(upper_R+lower_R);

		if( dr_f < FRAC_TOL || upper_R - lower_R < EXTREME_R)
		{
			bisection_done = 1;
			*distance = lower_R;
		}
#else	
       		if( (upper_R - lower_R < TOL_BISECTION))
       		{
       			bisection_done = 1;
       			*distance = lower_R; // not trial_R I guess
       		}
#endif
	}

//	printf("nunchecked: %d\n", nuncheck );

//	printf("%lf skipped %lf not skipped.\n", nskipped, nnot_skipped );

	free( box_checked_list );
	free( box_blank );
	free( uncheck_these );

}

void surface::assembleNearList( double *pt_in, 
					int **f_list_io, double **puv_list_io, double **areas_io, int *npts, double **M, int mlow, int mhigh, double Rmax, int sub_limit, int disable_PBC )
{
	double pt[3] = { pt_in[0], pt_in[1], pt_in[2] };
	double Lx = PBC_vec[0][0];
	double Ly = PBC_vec[1][1];
	double Lz = PBC_vec[2][2];

	int nspace = 1;

	int *f_list = (int *)malloc( sizeof(int) * nspace );
	double *puv_list = (double *)malloc( sizeof(double) * 6 * nspace );
	double *areas    = (double *)malloc( sizeof(double) * nspace );
	*npts = 0;

/*	
	we do not need to loop back and check these boxes. 

*/
	int nbc = 0;
	int nbtot = xbox*ybox*zbox;

	if( ! disable_PBC )
	{
		while( pt[0] < -Lx/2 ) pt[0] += Lx;
		while( pt[0] >= Lx/2 ) pt[0] -= Lx;
		while( pt[1] < -Ly/2 ) pt[1] += Ly;
		while( pt[1] >= Ly/2 ) pt[1] -= Ly;
		while( pt[2] < -Lz/2 ) pt[2] += Lz;
		while( pt[2] >= Lz/2 ) pt[2] -= Lz;
	}

	double trial_R = Rmax;

	double tfudge = 0;

	double minx = pt[0] - trial_R - max_width_x - xlength - tfudge;
	double maxx = pt[0] + trial_R + max_width_x + xlength + tfudge;
	
	double miny = pt[1] - trial_R - max_width_x - xlength - tfudge;
	double maxy = pt[1] + trial_R + max_width_x + xlength + tfudge;
	
	double minz = pt[2] - trial_R - max_width_x - xlength - tfudge;
	double maxz = pt[2] + trial_R + max_width_x + xlength + tfudge;

	
	int cenx = xbox * pt[0] / Lx; 
	int ceny = ybox * pt[1] / Ly; 
	int cenz = zbox * pt[2] / Lz;

	int binloop_x = xbox * (trial_R + max_width_x + xlength) / Lx;
	int binloop_y = ybox * (trial_R + max_width_y + ylength) / Ly;
	int binloop_z = zbox * (trial_R + max_width_z + zlength) / Lz;

	if( binloop_x >= xbox/2 ) binloop_x = xbox/2;
	if( binloop_y >= ybox/2 ) binloop_y = ybox/2;
	if( binloop_z >= zbox/2 ) binloop_z = zbox/2;

	int fudge_x = 0;
	int fudge_y = 0;
	int fudge_z = 0;

	if( binloop_x * 2 == xbox )
		fudge_x = -1;
	if( binloop_y * 2 == ybox )
		fudge_y = -1;
	if( binloop_z * 2 == zbox )
		fudge_z = -1;


	int did_collide = 0;
	int numChecks = 0;
	int numTris = 0;
	nbc = 0;
		
	int arrayX[xbox];
	int arrayY[ybox];
	int arrayZ[zbox];
	int sizeX = 0;
	int sizeY = 0;
	int sizeZ = 0;
	
	for( int t_xB = - binloop_x; t_xB <= binloop_x + fudge_x; t_xB++, sizeX++ ) 
		arrayX[sizeX] = t_xB;
	qsort(arrayX, sizeX, sizeof(int), cmpfunci);
	for( int t_yB = - binloop_y; t_yB <= binloop_y + fudge_y; t_yB++, sizeY++ ) 
		arrayY[sizeY] = t_yB;
	qsort(arrayY, sizeY, sizeof(int), cmpfunci);
	for( int t_zB = - binloop_z; t_zB <= binloop_z + fudge_z; t_zB++, sizeZ++ ) 
		arrayZ[sizeZ] = t_zB;
	qsort(arrayZ, sizeZ, sizeof(int), cmpfunci);

	for( int t_xB = 0; t_xB < sizeX && !did_collide; t_xB++ ) 
	{
		int xB = xbox/2 + cenx + arrayX[t_xB];

		if( xB < 0 ) xB += xbox;
		if( xB >= xbox ) xB -= xbox;

		for( int t_yB = 0; t_yB < sizeY && !did_collide; t_yB++ ) 
		{ 
			int yB = ybox/2 + ceny + arrayY[t_yB];
			
			if( yB < 0 ) yB += ybox;
			if( yB >= ybox ) yB -= ybox;

			for( int t_zB = 0; t_zB < sizeZ && !did_collide; t_zB++ )
			{
				int zB = zbox/2 + cenz + arrayZ[t_zB];
			
				if( zB < 0 ) zB += zbox;
				if( zB >= zbox ) zB -= zbox;
				bool worthCheck = true;	

				int box_id = (xB*ybox+yB)*zbox+zB;

#ifndef SKIP_SUPPORT_CHECK	
				if (box[xB][yB][zB].num_tri > 0) 
				{
					double support = pt[0]*box[xB][yB][zB].search[0] + pt[1]*box[xB][yB][zB].search[1] + pt[2]*box[xB][yB][zB].search[2];
					//double mindiff = abs(support - box[xB][yB][zB].vmin);
					//double maxdiff = abs(support - box[xB][yB][zB].vmax);

					double test = support + trial_R;
					double test2 = support - trial_R;
					double mindiff = box[xB][yB][zB].vmin;
					double maxdiff = box[xB][yB][zB].vmax;

/*					printf("support %le %le %le min %le max %le\n",
					 box[xB][yB][zB].search[0],
					 box[xB][yB][zB].search[1],
					 box[xB][yB][zB].search[2], mindiff, maxdiff );
*/
					if ((test < mindiff && test2 < mindiff) || (test > maxdiff && test2 > maxdiff)) {
						worthCheck = false;
					} else {
						numChecks++;
					}
				}
#endif

				for (int ind = 0; ind < box[xB][yB][zB].num_tri && !did_collide && worthCheck; ind++) {
					numTris++;
					triangle *tri1 = theTriangles + box[xB][yB][zB].tri_list[ind];
                			int f1 = tri1->f;

					int *indices1;
              				int val1 = 0;
                			int np1 = 0;
                			int base1 = 0;
					double *pbc1;

                			if( f1 >= nf_faces )
                			{

                        			f1 -= nf_faces;
                        			int formulas_per_face = nf_irr_pts;
                        			indices1 = theIrregularFormulas[f1*formulas_per_face].cp;
						pbc1 = theIrregularFormulas[f1*formulas_per_face].r_pbc;
                        			np1 = theIrregularFormulas[f1*formulas_per_face].ncoor;
                        			base1 = theIrregularFormulas[f1*formulas_per_face].vertex;
                			}
                			else
                			{
                        			int formulas_per_face = nf_g_q_p;
                        			indices1 = theFormulas[f1*formulas_per_face].cp;
						pbc1 = theFormulas[f1*formulas_per_face].r_pbc;
                        			np1 = theFormulas[f1*formulas_per_face].ncoor;
                        			base1 = theFormulas[f1*formulas_per_face].vertex;
                			}

                			double pts1[3*np1];
                			for( int x = 0; x < np1; x++ )
                			{
                        			pts1[3*x+0] = theVertices[indices1[x]].r[0] + pbc1[3*x+0];
                        			pts1[3*x+1] = theVertices[indices1[x]].r[1] + pbc1[3*x+1];
                        			pts1[3*x+2] = theVertices[indices1[x]].r[2] + pbc1[3*x+2];
						
                			}
				
					if( *npts + pow( 4, sub_limit+1 ) >= nspace )
					{
						nspace += pow( 4, sub_limit+1 );
						areas = (double *)realloc( areas, sizeof(double) * nspace );
						puv_list = (double *)realloc( puv_list, sizeof(double) * 6 * nspace );
						f_list = (int*)realloc( f_list, sizeof(int) * nspace );
					}

					int cur_n_pts = *npts;

                			nearMembraneWorker( pts1, np1, pt, M, mlow, mhigh, 0, sub_limit, trial_R,
							1.0, 0, 0, puv_list, areas, npts ); 

					for( int x = cur_n_pts; x < *npts; x++ )
						f_list[x] = tri1->f;
				} 
        		}
		}
	}
	*f_list_io = f_list;
	*puv_list_io = puv_list;
	*areas_io = areas;

#ifdef NEAR_POINT_DEBUG
	printf("nearMembrane found %d points.\n", *npts );

	double *rv = (double *)malloc( sizeof(double) * 3 * (nv+1) );
	
	get(rv);

	rv[3*nv+0] = 1;
	rv[3*nv+1] = 1;
	rv[3*nv+2] = 1;
	

	for( int p = 0; p < *npts; p++ )
	{
		int f = (*f_list_io)[p];
		double u = ((*puv_list_io)[p*6+0] + (*puv_list_io)[p*6+2] + (*puv_list_io)[p*6+4])/3;
		double v = ((*puv_list_io)[p*6+1] + (*puv_list_io)[p*6+3] + (*puv_list_io)[p*6+5])/3;

		double rcheck[3], ncheck[3];

		evaluateRNRM( f, u, v, rcheck, ncheck, rv ); 

		double r = sqrt( 
				(pt[0]-rcheck[0])*(pt[0]-rcheck[0]) +
				(pt[1]-rcheck[1])*(pt[1]-rcheck[1]) +
				(pt[2]-rcheck[2])*(pt[2]-rcheck[2])  );

		printf("r: %le %le %le distance: %lf vs. cutoff %lf\n", rcheck[0], rcheck[1], rcheck[2], r, trial_R ); 
	}
#endif

}

int nearMembraneWorker( double *r1, int nv1, double *r2, double **M, int mlow, int mhigh, int level, int max_level, double radius, 
		double factor, 
		double cur_u, double cur_v,
		double *puv_list, double *areas, int *npts )	
{

	bool min = minimum_distance( r1, nv1, r2, radius, level!=max_level );

	if( !min ) {
		return 0;
	}

	if( level == max_level )
	{
/*		double eval_p[3];
		double eval_n[3];

		double avp[3] = { 
			(r1[0] + r1[3] + r1[6])/3,	
			(r1[1] + r1[4] + r1[7])/3,	
			(r1[2] + r1[5] + r1[8])/3 };
			
		double dr[3] = { r2[0] - avp[0], r2[1] - avp[1], r2[2] - avp[2] }; 

		double r2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];

		if( r2 < radius*radius )
		{ 
			printf("pt: %le %le %le r: %le radius: %le\n", avp[0], avp[1], avp[2], sqrt(r2), radius );

			areas[*npts] = triangle_area( r1, r1+3, r1+6 );
			// this is the middle of the triangle.
			puv_list[(*npts) *2 +0] = cur_u + factor * (1.0/3.0);
			puv_list[(*npts) *2 +1] = cur_v + factor * (1.0/3.0);
			(*npts) += 1;
		} 
*/
		// the triangle.

		puv_list[(*npts)*6+0] = cur_u;
		puv_list[(*npts)*6+1] = cur_v;

		puv_list[(*npts)*6+2] = cur_u + factor;
		puv_list[(*npts)*6+3] = cur_v;

		puv_list[(*npts)*6+4] = cur_u;
		puv_list[(*npts)*6+5] = cur_v + factor;
		areas[*npts] = 1.0;	
		(*npts) += 1;
	

		return 1;
	}

	int nvmax = 6+mhigh;

	for( int x1 = 0; x1 < 4; x1++ )
	{
		double sub1[nvmax*3];
		memset( sub1, 0, sizeof(double)*nvmax*3);

		int nv_use1 = 12;
		if( x1 == 0 ) nv_use1 = nv1;
	
		double *mat1 = M[nv1-6-mlow]+x1*(nv1 < 12 ? 12 : nv1)*nv1;

		if( nv1 != 12 )
		{
//			printf("Using matrix %d (%d vertices).\n", nv1-6-mlow, nv1 );
		}
		double one = 1.0, zero = 0.0;
		int incrxy=3;

		double check1[nvmax*3];
		memset(check1, 0, sizeof(double) * 3 * nvmax );
		char trans = 'T';
		int nI = nv1;
		int nO = nv_use1;
#define USE_DGEMV
#ifdef USE_DGEMV
		if(  nv1 == 12  )
		{
			if( x1 == 0 )
				fast_sub6_0( r1, sub1 );
			else if( x1 == 1 )
				fast_sub6_1( r1, sub1 );
			else if( x1 == 2 )
				fast_sub6_2( r1, sub1 );
			else if( x1 == 3 )
				fast_sub6_3( r1, sub1 );
		}
		else
		{
			for( int c = 0; c < 3; c++ )
				dgemv( &trans, &nI, &nO, &one, mat1, &nv1, r1+c, &incrxy, &zero, sub1+c, &incrxy );
		}
#else
		for( int i = 0; i < nv_use1; i++ )
		{
			for( int j = 0; j < nv1; j++ )
			{
				sub1[i*3+0] += mat1[i*nv1+j] * r1[3*j+0];	
				sub1[i*3+1] += mat1[i*nv1+j] * r1[3*j+1];	
				sub1[i*3+2] += mat1[i*nv1+j] * r1[3*j+2];	
			}
		}
#endif

		double new_u = cur_u;
		double new_v = cur_v;
		
		double new_factor = factor;

		switch(x1)
		{
			case 0:
				new_factor *= 0.5;
				break;
			case 1:
				new_u += factor * 0.5;
				new_factor *= 0.5;
				break;
			case 2:
				new_v += factor * 0.5;
				new_factor *= 0.5;
				break;
			case 3:
				new_u += factor * 0.5;
				new_v += factor * 0.5;
				new_factor *= -0.5;
				break;
		
		}
                			
		nearMembraneWorker( sub1, nv_use1, r2, M, mlow, mhigh, level+1, max_level, radius,
				new_factor, new_u, new_v, 
					puv_list, areas, npts ); 
	}


	return 0;
} 

// returns 0 for no collision.

int surface::linearCollisionPoint( double *pt1_in, double *pt2_in, int *col_f, double *col_u, double *col_v, double **M, int mlow, int mhigh, int disable_PBC_z )
{
	double Lx = PBC_vec[0][0];
	double Ly = PBC_vec[1][1];
	double Lz = PBC_vec[2][2];

	*col_f = -1;
	*col_u = 0;
	*col_v = 0;

	
	double pt2[3] = { pt2_in[0], pt2_in[1], pt2_in[2] };	
	double pt1[3] = { pt1_in[0], pt1_in[1], pt1_in[2] };	

	while( pt1[0] < -Lx/2 ) pt1[0] += Lx;
	while( pt1[0] >= Lx/2 ) pt1[0] -= Lx;
	while( pt1[1] < -Ly/2 ) pt1[1] += Ly;
	while( pt1[1] >= Ly/2 ) pt1[1] -= Ly;
	while( pt1[2] < -Lz/2 ) pt1[2] += Lz;
	while( pt1[2] >= Lz/2 ) pt1[2] -= Lz;

	double dr[3] = { pt2[0] - pt1[0], pt2[1] - pt1[1], pt2[2] - pt1[2] };

	while( dr[0] < -Lx/2 ) dr[0] += Lx;
	while( dr[0] >= Lx/2 ) dr[0] -= Lx;
	while( dr[1] < -Ly/2 ) dr[1] += Ly;
	while( dr[1] >= Ly/2 ) dr[1] -= Ly;
	while( dr[2] < -Lz/2 ) dr[2] += Lz;
	while( dr[2] >= Lz/2 ) dr[2] -= Lz;

	pt2[0] = pt1[0] + dr[0];
	pt2[1] = pt1[1] + dr[1];
	pt2[2] = pt1[2] + dr[2];

	double pt[3] = { pt1[0] + dr[0]/2, pt1[1] + dr[1]/2, pt1[2] + dr[2]/2 };
	double trial_R = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2])/2;

	double tfudge = 0;

	double minx = pt[0] - trial_R - max_width_x - xlength - tfudge;
	double maxx = pt[0] + trial_R + max_width_x + xlength + tfudge;
	
	double miny = pt[1] - trial_R - max_width_x - xlength - tfudge;
	double maxy = pt[1] + trial_R + max_width_x + xlength + tfudge;
	
	double minz = pt[2] - trial_R - max_width_x - xlength - tfudge;
	double maxz = pt[2] + trial_R + max_width_x + xlength + tfudge;

	
	int cenx = xbox * pt[0] / Lx; 
	int ceny = ybox * pt[1] / Ly; 
	int cenz = zbox * pt[2] / Lz;

	int binloop_x = xbox * (trial_R + max_width_x + xlength) / Lx;
	int binloop_y = ybox * (trial_R + max_width_y + ylength) / Ly;
	int binloop_z = zbox * (trial_R + max_width_z + zlength) / Lz;

	if( binloop_x >= xbox/2 ) binloop_x = xbox/2;
	if( binloop_y >= ybox/2 ) binloop_y = ybox/2;
	if( binloop_z >= zbox/2 ) binloop_z = zbox/2;

	int fudge_x = 0;
	int fudge_y = 0;
	int fudge_z = 0;

	if( binloop_x * 2 == xbox )
		fudge_x = -1;
	if( binloop_y * 2 == ybox )
		fudge_y = -1;
	if( binloop_z * 2 == zbox )
		fudge_z = -1;


	int found_pt = 0;
	nbc = 0;
	
	int arrayX[xbox];
	int arrayY[ybox];
	int arrayZ[zbox];
	int sizeX = 0;
	int sizeY = 0;
	int sizeZ = 0;
	
	for( int t_xB = - binloop_x; t_xB <= binloop_x + fudge_x; t_xB++, sizeX++ ) 
		arrayX[sizeX] = t_xB;
	qsort(arrayX, sizeX, sizeof(int), cmpfunci);
	for( int t_yB = - binloop_y; t_yB <= binloop_y + fudge_y; t_yB++, sizeY++ ) 
		arrayY[sizeY] = t_yB;
	qsort(arrayY, sizeY, sizeof(int), cmpfunci);
	for( int t_zB = - binloop_z; t_zB <= binloop_z + fudge_z; t_zB++, sizeZ++ ) 
		arrayZ[sizeZ] = t_zB;
	qsort(arrayZ, sizeZ, sizeof(int), cmpfunci);

//#define DEBUG_LINEAR_COLLISION

#ifdef DEBUG_LINEAR_COLLISION
	static int dbg_cntr = 0;
	double *rv = (double *)malloc( sizeof(double) * 3 * (nv+3) );
	rv[3*nv+0] = 1.0;
	rv[3*nv+1] = 1.0;
	rv[3*nv+2] = 1.0;
#endif

	for( int t_xB = 0; t_xB < sizeX && !found_pt; t_xB++ ) 
	{
		int xB = xbox/2 + cenx + arrayX[t_xB];

		if( xB < 0 ) xB += xbox;
		if( xB >= xbox ) xB -= xbox;

		for( int t_yB = 0; t_yB < sizeY && !found_pt; t_yB++ ) 
		{ 
			int yB = ybox/2 + ceny + arrayY[t_yB];
			
			if( yB < 0 ) yB += ybox;
			if( yB >= ybox ) yB -= ybox;

			for( int t_zB = 0; t_zB < sizeZ && !found_pt; t_zB++ )
			{
				int zB = zbox/2 + cenz + arrayZ[t_zB];
			
				if( zB < 0 ) zB += zbox;
				if( zB >= zbox ) zB -= zbox;

				bool worthCheck = true;	

				int box_id = (xB*ybox+yB)*zbox+zB;


#ifndef SKIP_SUPPORT_CHECK	
				if (box[xB][yB][zB].num_tri > 0) {
					double support = pt[0]*box[xB][yB][zB].search[0] + pt[1]*box[xB][yB][zB].search[1] + pt[2]*box[xB][yB][zB].search[2];
					//double mindiff = abs(support - box[xB][yB][zB].vmin);
					//double maxdiff = abs(support - box[xB][yB][zB].vmax);

					double test = support + trial_R;
					double test2 = support - trial_R;
					double mindiff = box[xB][yB][zB].vmin;
					double maxdiff = box[xB][yB][zB].vmax;

					if ((test < mindiff && test2 < mindiff) || (test > maxdiff && test2 > maxdiff)) 
						worthCheck = false;
				}
#endif

				for (int ind = 0; ind < box[xB][yB][zB].num_tri && !found_pt && worthCheck; ind++) {
					triangle *tri1 = theTriangles + box[xB][yB][zB].tri_list[ind];
                			int f1 = tri1->f;

					int *indices1;
              				int val1 = 0;
                			int np1 = 0;
                			int base1 = 0;
					double *pbc1;

                			if( f1 >= nf_faces )
                			{

                        			f1 -= nf_faces;
                        			int formulas_per_face = nf_irr_pts;
                        			indices1 = theIrregularFormulas[f1*formulas_per_face].cp;
						pbc1 = theIrregularFormulas[f1*formulas_per_face].r_pbc;
                        			np1 = theIrregularFormulas[f1*formulas_per_face].ncoor;
                        			base1 = theIrregularFormulas[f1*formulas_per_face].vertex;
                			}
                			else
                			{
                        			int formulas_per_face = nf_g_q_p;
                        			indices1 = theFormulas[f1*formulas_per_face].cp;
						pbc1 = theFormulas[f1*formulas_per_face].r_pbc;
                        			np1 = theFormulas[f1*formulas_per_face].ncoor;
                        			base1 = theFormulas[f1*formulas_per_face].vertex;
                			}

                			double pts1[3*np1];
                			for( int x = 0; x < np1; x++ )
                			{
                        			pts1[3*x+0] = theVertices[indices1[x]].r[0] + pbc1[3*x+0];;
                        			pts1[3*x+1] = theVertices[indices1[x]].r[1] + pbc1[3*x+1];;
                        			pts1[3*x+2] = theVertices[indices1[x]].r[2] + pbc1[3*x+2];;
						
                			}
				
					double usep1[3] = { pt1[0], pt1[1], pt1[2] };
					double usep2[3] = { pt2[0], pt2[1], pt2[2] };
					double usep[3] = { pt[0], pt[1], pt[2] };

					while( usep1[0] - pts1[0] < -Lx/2 ) { usep1[0] += Lx; usep2[0] += Lx; usep[0] += Lx; }
					while( usep1[1] - pts1[1] < -Ly/2 ) { usep1[1] += Ly; usep2[1] += Ly; usep[1] += Ly; }
					if( !disable_PBC_z ) 
						while( usep1[2] - pts1[2] < -Lz/2 ) { usep1[2] += Lz; usep2[2] += Lz; usep[2] += Lz; }
					while( usep1[0] - pts1[0] > Lx/2 ) { usep1[0] -= Lx; usep2[0] -= Lx; usep[0] -= Lx; }
					while( usep1[1] - pts1[1] > Ly/2 ) { usep1[1] -= Ly; usep2[1] -= Ly; usep[1] -= Ly; }
					if( !disable_PBC_z ) 
						while( usep1[2] - pts1[2] > Lz/2 ) {usep1[2] -= Lz; usep2[2] -= Lz; usep[2] -= Lz; }
					
	
                			if(  checkLinearCollision( pts1, np1, usep1, usep2, usep, M, mlow, mhigh, 0, COLLISION_LEVEL, trial_R,
                        			1.0, 0., 0., col_u, col_v ) )
                			{
                        			*col_f =tri1->f;
#ifdef DEBUG_LINEAR_COLLISION
						double rcol[3],ncol[3];
						evaluateRNRM( tri1->f, *col_u, *col_v, rcol, ncol, rv );
						printf("%d DEBUG collision at %le %le %le\n", dbg_cntr, rcol[0], rcol[1], rcol[2] );	
#else
                        			found_pt = 1;

#endif
	        			}
				} 
        		}
		}
	}

#ifdef DEBUG_LINEAR_COLLISION
	dbg_cntr++;
	free(rv);
#endif

	return found_pt;
}

void  surface::buildFaceData(  double **vertex_data_in, int *pointer_to_data, int *nump )
{
	int of = 0;
	int data_size = 0;

	double *data = NULL;

	for( int pass = 0; pass < 2; pass++ )
	{
		if( pass == 1 )
			data = (double *)malloc( sizeof(double) * data_size );

		data_size = 0;

		for( int f = 0; f < nt; f++ )
		{
			pointer_to_data[f] = data_size;
			int f1 = f;

			int *indices1;
			double *pbc1;
			int np1;
			int base1;

        	        if( f1 >= nf_faces )
        	        {	
	                	f1 -= nf_faces;
	                	int formulas_per_face = nf_irr_pts;
	                	indices1 = theIrregularFormulas[f1*formulas_per_face].cp;
				pbc1 = theIrregularFormulas[f1*formulas_per_face].r_pbc;
	                	np1 = theIrregularFormulas[f1*formulas_per_face].ncoor;
	                	base1 = theIrregularFormulas[f1*formulas_per_face].vertex;
	                }
	                else
	                {
	                	int formulas_per_face = nf_g_q_p;
	                	indices1 = theFormulas[f1*formulas_per_face].cp;
				pbc1 = theFormulas[f1*formulas_per_face].r_pbc;
	                	np1 = theFormulas[f1*formulas_per_face].ncoor;
	                	base1 = theFormulas[f1*formulas_per_face].vertex;
	                }

			if( pass == 1 )
			{
				for( int x = 0; x < np1; x++ )
				{
                			data[data_size+3*x+0] = theVertices[indices1[x]].r[0] + pbc1[3*x+0];;
                			data[data_size+3*x+1] = theVertices[indices1[x]].r[1] + pbc1[3*x+1];;
                			data[data_size+3*x+2] = theVertices[indices1[x]].r[2] + pbc1[3*x+2];;
				}
			}

			nump[f] = np1;
	
			data_size += 3 * np1;
		} 
	}

	*vertex_data_in = data;
}

int collisionForces( double *r1, int nv1, 
		     double *r2, int nv2, 
		     double **M, int mlow, int mhigh, 
		     int level, int max_level,
		     double *f1, double *f2, double cutoff, double alpha )
{
	int pts[3] = { 0, 1, 2};

	double min1[3]={r1[0],r1[1],r1[2]};
	double max1[3]={r1[0],r1[1],r1[2]};
	
	double min2[3]={r2[0],r2[1],r2[2]};
	double max2[3]={r2[0],r2[1],r2[2]};

	for( int c = 0; c < 3; c++ )
	{
		for( int v = 0; v < nv1; v++ )
		{
			if( r1[c] < min1[c] )
				min1[c] = r1[c];
			if( r1[c] > max1[c] )
				max1[c] = r1[c];
		}
		for( int v = 0; v < nv2; v++ )
		{
			if( r2[c] < min2[c] )
				min2[c] = r2[c];
			if( r2[c] > max2[c] )
				max2[c] = r2[c];
		}

		if( max1[c] + cutoff < min2[c] ||
		    max2[c] + cutoff < min1[c] )
			return 0;
	}	

	if( level == max_level )
	{
		double dr[3] = { 
			(r1[0]+r1[3]+r1[6])/3-(r2[0]+r2[3]+r2[6])/3,
			(r1[1]+r1[4]+r1[7])/3-(r2[1]+r2[4]+r2[7])/3,
			(r1[2]+r1[5]+r1[8])/3-(r2[2]+r2[5]+r2[8])/3 
				};
		double r = normalize(dr);
	
		double pref =  (1.0/alpha) * exp( -r / alpha );
	
		// -dvdr, negatives cancel.
		// center of geom of the first three points

		for( int p = 0; p < 3; p++ )
		{
			f1[p*3+0] += pref * dr[0];
			f1[p*3+1] += pref * dr[1];
			f1[p*3+2] += pref * dr[2];
			
			f2[p*3+0] -= pref * dr[0];
			f2[p*3+1] -= pref * dr[1];
			f2[p*3+2] -= pref * dr[2];
		}

		return 1;
	}

	int procd=0;
	
	int nvmax = 6+mhigh;

	for( int x1 = 0; x1 < 4; x1++ )
	{
		double sub1[nvmax*3];
		memset( sub1, 0, sizeof(double)*nvmax*3);
		double frc1[nvmax*3];
		memset( frc1, 0, sizeof(double) * nvmax*3);

		int do_proc1 = 0;
		int nv_use1 = 12;
		if( x1 == 0 ) nv_use1 = nv1;
	
		double *mat1 = M[nv1-6-mlow]+x1*(nv1 < 12 ? 12 : nv1)*nv1;

		for( int i = 0; i < nv_use1; i++ )
		{
			for( int j = 0; j < nv1; j++ )
			{
				sub1[i*3+0] += mat1[i*nv1+j] * r1[3*j+0];	
				sub1[i*3+1] += mat1[i*nv1+j] * r1[3*j+1];	
				sub1[i*3+2] += mat1[i*nv1+j] * r1[3*j+2];	
			}
		}

		for( int x2 = 0; x2 < 4; x2++ )
		{
			double sub2[nvmax*3];
			memset( sub2, 0, sizeof(double)*nvmax*3);
			double frc2[nvmax*3];
			memset( frc2, 0, sizeof(double) * nvmax*3);

			int nv_use2 = 12;
			if( x2 == 0 ) nv_use2 = nv2;
			double *mat2 = M[nv2-6-mlow]+x2*(nv2 < 12 ? 12 : nv2)*nv2;

			for( int i = 0; i < nv_use2; i++ )
			{
				for( int j = 0; j < nv2; j++ )
				{
					sub2[i*3+0] += mat2[i*nv2+j] * r2[3*j+0];	
					sub2[i*3+1] += mat2[i*nv2+j] * r2[3*j+1];	
					sub2[i*3+2] += mat2[i*nv2+j] * r2[3*j+2];	
				}
			}	

			int do_proc = collisionForces( 
					sub1, nv_use1, 
					sub2, nv_use2, 
					M, mlow, mhigh, 
					level+1, max_level,
					frc1, frc2, cutoff, alpha ); 

			if( do_proc )
			{
				do_proc1 = 1;

				// transform forces back to original
			
				for( int i = 0; i < nv_use2; i++ )
				{
					for( int j = 0; j < nv2; j++ )
					{
						f2[j*3+0] += mat2[i*nv2+j] * frc2[3*i+0];	
						f2[j*3+1] += mat2[i*nv2+j] * frc2[3*i+1];	
						f2[j*3+2] += mat2[i*nv2+j] * frc2[3*i+2];	
					}
				}	 
			}
		}
		
		if( do_proc1 )
		{
			procd = 1;

			for( int i = 0; i < nv_use1; i++ )
			{
				for( int j = 0; j < nv1; j++ )
				{
					f1[j*3+0] += mat1[i*nv1+j] * frc1[3*i+0];	
					f1[j*3+1] += mat1[i*nv1+j] * frc1[3*i+1];	
					f1[j*3+2] += mat1[i*nv1+j] * frc1[3*i+2];	
				}
			}
		}
	}

	return procd;
} 

void surfaceSurfaceCollisionForces( surface *surface1, surface *surface2, double *grad1, double *grad2, double alpha, double v0, double **M, int mlow, int mhigh )
{
	// V = rho1 rho2 v0 exp(-r/alpha)
	
	int eval_level = 1; // eval level one: triangle centers.
	double cutoff = 10 * alpha; // where we neglect exp(-r/alpha)

	// f_mesh: -dv dr_mesh, etc

	
	for( int t = 0; t < surface1->nt; t++ )
	{
		triangle *tri1 = surface1->theTriangles+t;

		// f is the index into my structures for computing properties.
		int f1 = tri1->f;

		int *indices1;
		int val1 = 0;
		int np1 = 0;
		int base1 = 0;
		double *pbc1;

		if( f1 >= surface1->nf_faces )
		{
			f1 -= surface1->nf_faces;
			int formulas_per_face = surface1->nf_irr_pts;
			indices1 = surface1->theIrregularFormulas[f1*formulas_per_face].cp;
			pbc1 = surface1->theIrregularFormulas[f1*formulas_per_face].r_pbc;
			np1 = surface1->theIrregularFormulas[f1*formulas_per_face].ncoor;
			base1 = surface1->theIrregularFormulas[f1*formulas_per_face].vertex;
		}
		else
		{
			int formulas_per_face = surface1->nf_g_q_p;
			indices1 = surface1->theFormulas[f1*formulas_per_face].cp;
			pbc1 = surface1->theFormulas[f1*formulas_per_face].r_pbc;
			np1 = surface1->theFormulas[f1*formulas_per_face].ncoor;
			base1 = surface1->theFormulas[f1*formulas_per_face].vertex;
		}			

		double pts1[3*np1];
		double frc1[3*np1];
		memset( frc1, 0, sizeof(double) * np1 );
		for( int x = 0; x < np1; x++ )
		{
			pts1[3*x+0] = surface1->theVertices[indices1[x]].r[0] + pbc1[3*x+0]; 
			pts1[3*x+1] = surface1->theVertices[indices1[x]].r[1] + pbc1[3*x+1]; 
			pts1[3*x+2] = surface1->theVertices[indices1[x]].r[2] + pbc1[3*x+2];
		} 

		for( int t2 = 0; t2 < surface2->nt; t2++ )
		{
			triangle *tri2 = surface2->theTriangles+t2;

			// f is the index into my structures for computing properties.
			int f2 = tri2->f;
	
			int *indices2;
			int val2 = 0;
			int np2 = 0;
			int base2 = 0;
			double *pbc2;
	
			if( f2 >= surface2->nf_faces )
			{
				f2 -= surface2->nf_faces;
				int formulas_per_face = surface2->nf_irr_pts;
				indices2 = surface2->theIrregularFormulas[f2*formulas_per_face].cp;
				pbc2 = surface2->theIrregularFormulas[f2*formulas_per_face].r_pbc;
				np2 = surface2->theIrregularFormulas[f2*formulas_per_face].ncoor;
				base2 = surface2->theIrregularFormulas[f2*formulas_per_face].vertex;
			}
			else
			{
				int formulas_per_face = surface2->nf_g_q_p;
				indices2 = surface2->theFormulas[f2*formulas_per_face].cp;
				pbc2 = surface2->theFormulas[f2*formulas_per_face].r_pbc;
				np2 = surface2->theFormulas[f2*formulas_per_face].ncoor;
				base2 = surface2->theFormulas[f2*formulas_per_face].vertex;
			}			

			double pts2[3*np2];
			double frc2[3*np2];
			memset( frc2, 0, sizeof(double) * np2 );
			for( int x = 0; x < np2; x++ )
			{
				pts2[3*x+0] = surface2->theVertices[indices2[x]].r[0] + pbc2[3*x+0]; 
				pts2[3*x+1] = surface2->theVertices[indices2[x]].r[1] + pbc2[3*x+1]; 
				pts2[3*x+2] = surface2->theVertices[indices2[x]].r[2] + pbc2[3*x+2];
			} 

			collisionForces( pts1, np1, pts2, np2, M, mlow, mhigh, 0, eval_level,
				frc1, frc2, cutoff, alpha );
			
			for( int x = 0; x < np2; x++ )
			{
				grad2[indices2[x]*3+0] += frc2[3*x+0];				
				grad2[indices2[x]*3+1] += frc2[3*x+1];				
				grad2[indices2[x]*3+2] += frc2[3*x+2];				
			} 
						
//			if(  checkCollision( pts1, np1, pts2, np2, M, mlow, mhigh, 0, USE_MAX_LEVEL ) )
//				return 1;
		}
			
		for( int x = 0; x < np1; x++ )
		{
			grad1[indices1[x]*3+0] += frc1[3*x+0];				
			grad1[indices1[x]*3+1] += frc1[3*x+1];				
			grad1[indices1[x]*3+2] += frc1[3*x+2];				
		} 
		 
	}
}
























