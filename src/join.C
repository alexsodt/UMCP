//#include <OpenCL/opencl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "util.h"
#include "interp.h"
#include "l-bfgs.h"
#include "m_triangles.h"
#include "parallel.h"
#include "mutil.h"
#include "mesh_utility.h"
#include "alignSet.h"
//#define SKIP_OPT

#define FIX_VOLUME

surface *minSurface = NULL;
double thick = 13.8;
extern int do_print_area_ps;
extern double av_vol_err;
extern double KV_scale;
extern double KA ;
void fdiff_check( double *coords );
int print_energy = 0;
extern double tilt_scale;
surface *upperSurface;
surface *lowerSurface;
volumes *theVolumes; 
 void amoeba(double **p, double y[], int ndim, double ftol,
               double (*funk)(double []), int *nfunk);
void amoebaOpt( double *r, int nv, double *rc, int nrc );
double surface_grad( double *r, double *g );
double surface_e( double *r );
	void curvatureGrad( double *g, surface *surface1, surface *surface2, double *r, double thick );
	double curvatureEnergy( surface *surface1, surface *surface2, double *r, double thick);
	void tiltGrad( double *g, surface *surface1, surface *surface2, double *r );
	void tiltGradExt( double *g, surface *surface1, surface *surface2, double *r );
	double tiltEnergy( surface *surface1, surface *surface2, double *r);
	double tiltEnergyExt( surface *surface1, surface *surface2, double *r, double *cons_energy);
double f(double *r);
double fdf( double *r, double *g);
void writeDualSurface( const char *fileNameXYZ, const char *fileNamePSF, surface *upperSurface, surface *lowerSurface, double *coords, double alpha );

void getAltCylinderMesh( double **rvals, int **edges, int **nedges, int *nv, int grain, double R, double LC);

int fix_alpha = 0;

void getVertPath( surface *theSurface, int vert, double join_radius, int **path_out, int *pathLen_out, double *r)
{
	// Find an ''edge path'' around object one that matches the radius.
	
	int n_do = 1000;

	int vert_path[n_do];

	double uvcen[2] = { 1.0/3, 1.0/3 };

	int f;

	if( theSurface->theVertices[vert].valence == 6 )
	{
		// start spot on the irregular vertex, gives the user a chance to pinpoint things.
		
		f = theSurface->theVertices[vert].faces[0];

		int *cp_use;
		if( f < theSurface->nf_faces )
			cp_use = theSurface->theFormulas[f*theSurface->nf_g_q_p].cp;
		else
			cp_use = theSurface->theIrregularFormulas[(f-theSurface->nf_faces)*theSurface->nf_irr_pts].cp;

		if( cp_use[0] == vert )
		{
			uvcen[0] = 0;
			uvcen[1] = 0;
		}
		else if( cp_use[1] == vert )
		{
			uvcen[0] = 0;
			uvcen[1] = 1;
		}
		else 
		{
			uvcen[0] = 1;
			uvcen[1] = 0;
		}
	}
	else
	{
		// start off the irregular vertex.

		f = theSurface->theTriangles[theSurface->theVertices[vert].faces[0]].f;
		uvcen[0] = 1.0/3;
		uvcen[1] = 1.0/3;
	}


	double cd1[2], cd2[2], c1,c2;
		

	theSurface->c( f, uvcen[0], uvcen[1], r, cd1, cd2, &c1, &c2 );
		
	double drdu[3];
	double drdv[3];

	theSurface->ru( f, uvcen[0], uvcen[1], r, drdu );	
	theSurface->rv( f, uvcen[0], uvcen[1], r, drdv );	

	double vec1[3] = { 
		drdu[0] * cd1[0] + drdv[0] * cd1[1],
		drdu[1] * cd1[0] + drdv[1] * cd1[1],
		drdu[2] * cd1[0] + drdv[2] * cd1[1] };
	
	double vec2[3] = { 
		drdu[0] * cd2[0] + drdv[0] * cd2[1],
		drdu[1] * cd2[0] + drdv[1] * cd2[1],
		drdu[2] * cd2[0] + drdv[2] * cd2[1] };
	double lv1 = normalize(vec1);
	double lv2 = normalize(vec2);



//	printf("%d\n", 2*n_do);
//	printf("circular path\n");

	// pick the sense to go along the normal properly

	double rcen[3], rnrm[3];

	theSurface->evaluateRNRM( f, uvcen[0], uvcen[1], rcen, rnrm, r );

	printf("r: %lf %lf %lf nrm: %lf %lf %lf\n", rcen[0], rcen[1], rcen[2], rnrm[0], rnrm[1], rnrm[2] );

	// the cd1 and cd2 directions are arbitrary and may require a flip

	double outp1[3], outp2[3];

	cross( vec1, vec2, outp1 );
	cross( drdu, drdv, outp2 );

	double dp = outp1[0]*outp2[0] + outp1[1] * outp2[1] + outp1[2]*outp2[2];
	double theta_sign = -1;

	if( dp < 0 )
		theta_sign *= -1;

	for( int i = 0; i < n_do; i++ )
	{
		double dx = cos( theta_sign * 2 * M_PI * i /(double)n_do );
		double dy = sin( theta_sign * 2 * M_PI * i /(double)n_do );

		double du = cd1[0] * join_radius * dx/lv1 + cd2[0] * join_radius * dy/lv2;
		double dv = cd1[1] * join_radius * dx/lv1 + cd2[1] * join_radius * dy/lv2;
		
		double u_cen = uvcen[0];
		double v_cen = uvcen[1];

		int f2=f, f_use;
		do 
		{	
			f_use = f2;

			f2 = theSurface->nextFace( f_use, &u_cen, &v_cen, &du, &dv, r );
		} while( f2 != f_use );
		
		double w = 1 - u_cen - v_cen;

		int *cp_use;
		if( f2 < theSurface->nf_faces )
			cp_use = theSurface->theFormulas[f2*theSurface->nf_g_q_p].cp;
		else
			cp_use = theSurface->theIrregularFormulas[(f2-theSurface->nf_faces)*theSurface->nf_irr_pts].cp;

		if( w > u_cen && w > v_cen) // u = small, v = small
			vert_path[i] = cp_use[0];
		else if( u_cen > v_cen ) // u large
			vert_path[i] = cp_use[2];
		else 
			vert_path[i] = cp_use[1];


		double r_eval[3], n_eval[3];

		theSurface->evaluateRNRM( f2, u_cen, v_cen, r_eval, n_eval, r );

//		printf("C %lf %lf %lf # vert %d\n", r[3*vert_path[i]+0], r[3*vert_path[i]+1], r[3*vert_path[i]+2], vert_path[i] );
//		printf("O %lf %lf %lf\n", r_eval[0], r_eval[1], r_eval[2] ); 
	}
	
	int mid_path[n_do];

	mid_path[0] = vert_path[0];
	int nmid = 1;
	int back_trip = 0;
	for( int i = 1; i < n_do; i++ )
	{
		if( vert_path[i] == mid_path[nmid-1] )
			continue; // same as previous point, move to next.

		if( vert_path[i] != mid_path[0] ) back_trip = 1;

		if( back_trip && vert_path[i] == mid_path[0] )
			break;

		mid_path[nmid] = vert_path[i];
		nmid++;	
	}	

//	printf("Mid-path, %d points.\n", nmid);	
	
	int final_path[nmid];
	int nfinal = 1;

	final_path[0] = mid_path[0];

	for( int i = 1; i < nmid; i++ )
	{
		// do we add i?

		int next_p = mid_path[0];

		if( i < nmid-1)
			next_p = mid_path[i+1];
		int here = final_path[nfinal-1];
		int skip = 0;
		for( int t = 0; t < theSurface->theVertices[here].valence; t++ )
		{
			if( theSurface->theVertices[here].edges[t] == next_p )
				skip = 1;	
		}
		
		if( skip ) continue;
		
		final_path[nfinal] = mid_path[i];
		nfinal++;
	}

	*path_out = (int*)malloc( sizeof(int) * nfinal );
	memcpy( *path_out, final_path, sizeof(int) * nfinal );
	*pathLen_out = nfinal;

//	printf("Pathlen: %d\n", nfinal );
}





int main( int argc, char **argv )
{
#ifdef PARALLEL
	int ierr;
	ierr = MPI_Init(&argc,&argv);
	if( ierr ) { exit(1); }
	quietParallel();
#endif
#ifdef FIXED_SEED
	srand(1);
#else
          struct timeval tp;
 
          gettimeofday( &tp, NULL );
 
          srand((long)tp.tv_usec);
#endif

	char buffer[4096];

	if( argc < 7 )
	{
		printf("Syntax: join mesh1 vert1 mesh2 vert2 radius cyl_len [auto]\n");
		printf("Auto toggles the radius until it can find a working mesh.\n");
		return -1;
	}
	
	int do_auto = 0;
	if( argc > 7 )
		do_auto = !strcasecmp( argv[7], "auto");

	surface *theSurface1 =(surface *)malloc( sizeof(surface) );
	theSurface1->loadLattice( argv[1], 0. );
	surface *theSurface2 =(surface *)malloc( sizeof(surface) );
	theSurface2->loadLattice( argv[3], 0. );

	theSurface1->generatePlan();
	theSurface2->generatePlan();

	double *r1 = (double *)malloc( sizeof(double) * (3 * theSurface1->nv+3) );
	double *r2 = (double *)malloc( sizeof(double) * (3 * theSurface2->nv+3) );

	int auto_done = 0;
	
	double dJoin = 1.0;
		
	double join_radius = atof(argv[5]);
	double join_len	   = atof(argv[6]);

	while( ! auto_done ) 
	{
		if( !do_auto ) 
			auto_done = 1;

		theSurface1->get(r1);
		theSurface1->setg0(r1);
		r1[theSurface1->nv*3] = 1.0;
		r1[theSurface1->nv*3+1] = 1.0;
		r1[theSurface1->nv*3+2] = 1.0;
		theSurface2->get(r2);
		theSurface2->setg0(r2);
		r2[theSurface2->nv*3] = 1.0;
		r2[theSurface2->nv*3+1] = 1.0;
		r2[theSurface2->nv*3+2] = 1.0;
	
		int mesh1_has_pbc = 0;
		int mesh2_has_pbc = 0;
	
		for( int v = 0; v < theSurface1->nv; v++ )
		{
			int val = theSurface1->theVertices[v].valence;
	
			for( int e = 0; e < val; e++ )
			{
				if( fabs(theSurface1->theVertices[v].edge_PBC[3*e+0]) > 1e-8 ) 
					mesh1_has_pbc=1;
				if( fabs(theSurface1->theVertices[v].edge_PBC[3*e+1]) > 1e-8 ) 
					mesh1_has_pbc=1;
				if( fabs(theSurface1->theVertices[v].edge_PBC[3*e+2]) > 1e-8 ) 
					mesh1_has_pbc=1;
			}
			if( mesh1_has_pbc ) break;
		}
	
		for( int v = 0; v < theSurface2->nv; v++ )
		{
			int val = theSurface2->theVertices[v].valence;
	
			for( int e = 0; e < val; e++ )
			{
				if( fabs(theSurface2->theVertices[v].edge_PBC[3*e+0]) > 1e-8 ) 
					mesh2_has_pbc=1;
				if( fabs(theSurface2->theVertices[v].edge_PBC[3*e+1]) > 1e-8 ) 
					mesh2_has_pbc=1;
				if( fabs(theSurface2->theVertices[v].edge_PBC[3*e+2]) > 1e-8 ) 
					mesh2_has_pbc=1;
			}
			if( mesh2_has_pbc ) break;
		}
	
	
		KA = 1.215;
	
		int vert1 = atoi(argv[2]);
	
		if( vert1 >= theSurface1->nv || vert1 < 0 )
		{
			printf("Invalid index for vert1.\n");
			exit(1);
	
		}
		int vert2 = atoi(argv[4]);
	
	
		int *path1;
		int pathLen1;
	
		void getVertPath( surface *theSurface, int vert, double join_radius, int **path_out, int *pathLen_out, double *r);
	
		printf("Getting the path to chop out of the first mesh.\n");
		getVertPath( theSurface1, vert1, join_radius, &path1, &pathLen1, r1); 
	
	#ifdef DEBUG_1
		printf("%d\n", theSurface1->nv );
		printf("okay\n");	
		for( int v = 0; v < theSurface1->nv; v++ )
		{
			int is_border = 0;
	
			for( int t = 0; t < pathLen1; t++ )
				if( path1[t] == v )
					is_border =1;
	
			if( is_border )
				printf("O %lf %lf %lf\n", r1[3*v+0], r1[3*v+1], r1[3*v+2] );
			else
				printf("C %lf %lf %lf\n", r1[3*v+0], r1[3*v+1], r1[3*v+2] );
		}
	#endif
	
		// now we need a path with the same length, approx, on the opposite mesh.
	
		int done = 0;
	
		double scale1 = 0.5;
		double scale2 = 1.5;
			
		int *path2;
		int pathLen2;
	
		r2[3*theSurface2->nv+0] = scale1;
		r2[3*theSurface2->nv+1] = scale1;
		r2[3*theSurface2->nv+2] = scale1;
	
	//	printf("Target is %d.\n", pathLen1 );
	
		printf("Getting the path to chop out of the second mesh and trying to get exactly the same length.\n");
		getVertPath( theSurface2, vert2, join_radius, &path2, &pathLen2, r2 );
	//	printf("scale1 path is %d.\n", pathLen2 );
	
		while( pathLen2 < pathLen1 )
		{
			free(path2);
			scale1 *= 0.95;
			r2[3*theSurface2->nv+0] = scale1;
			r2[3*theSurface2->nv+1] = scale1;
			r2[3*theSurface2->nv+2] = scale1;
	
			getVertPath( theSurface2, vert2, join_radius, &path2, &pathLen2, r2 );
	//		printf("new scale1 path is %d.\n", pathLen2 );
		}
	
		free(path2);
	
		r2[3*theSurface2->nv+0] = scale2;
		r2[3*theSurface2->nv+1] = scale2;
		r2[3*theSurface2->nv+2] = scale2;
		getVertPath( theSurface2, vert2, join_radius, &path2, &pathLen2, r2 );
	//		printf("scale2 path is %d.\n", pathLen2 );
	
		while( pathLen2 > pathLen1 )
		{
			free(path2);
			scale2 *= 1.05;
			r2[3*theSurface2->nv+0] = scale2;
			r2[3*theSurface2->nv+1] = scale2;
			r2[3*theSurface2->nv+2] = scale2;
	
			getVertPath( theSurface2, vert2, join_radius, &path2, &pathLen2, r2 );
	//		printf("new scale2 path is %d.\n", pathLen2 );
		}
		
		free(path2);
	
		double trial = (scale1+scale2)/2;
		int niter=0;
		while( !done && niter < 100 )
		{
			r2[3*theSurface2->nv+0] = trial;
			r2[3*theSurface2->nv+1] = trial;
			r2[3*theSurface2->nv+2] = trial;
	
			getVertPath( theSurface2, vert2, join_radius, &path2, &pathLen2, r2 );
	
			if( pathLen2 == pathLen1 )
				done = 1;
			else
			{
				free(path2);
					
			//	printf("trial: %lf scale1: %lf scale2: %lf\n", trial, scale1, scale2  );
				if( pathLen2 > pathLen1 )
					scale1 = trial;
				else
					scale2 = trial;
	
				trial = (scale1+scale2)/2;
			}	
	
			niter++;
		}	
	
		double mesh_scale_2 = trial;
		
		if( pathLen1 == pathLen2 )
			printf("Successfully found two paths of %d edges to chop out.\n", pathLen1 );
		else
		{
			printf("Couldn't find a good path for the second mesh.\n");
			exit(1);
		}	
		// mark ``interior'' vertices for deletion.
	
	
		int *vert_color1 = (int *)malloc( sizeof(int) * theSurface1->nv );
		int *vert_color2 = (int *)malloc( sizeof(int) * theSurface2->nv );
		int nv_final1 = 0;
		int nv_final2 = 0;
	
		for( int mesh = 0; mesh < 2; mesh ++ )
		{
			surface *colorSurface = theSurface1;
			int *color_array = vert_color1;
			int *path = path1;
			int pathLen = pathLen1;
			int centerVert = vert1;	
	
			done = 0;
	
			if( mesh == 1 )
			{
				centerVert = vert2;
				path = path2;
				pathLen = pathLen2;
				colorSurface = theSurface2;
				color_array = vert_color2;
			}
			int nv_final = colorSurface->nv;
		
			memset( color_array, 0, sizeof(int) * colorSurface->nv );
			for( int t = 0; t < pathLen; t++ )
				color_array[path[t]] = 2;
			color_array[centerVert] = 1;
			nv_final--;
			while(!done)
			{
				// color all the neighbors of the excised vertex, stopping at the ring.
				done = 1;
				for( int v = 0; v < colorSurface->nv; v++ )
				{
					if( color_array[v] == 1 )
					{
						for( int e = 0; e < colorSurface->theVertices[v].valence; e++ )
						{
							int v2 = colorSurface->theVertices[v].edges[e];
			
							if( color_array[v2] == 0 )
							{	
								nv_final--;
								color_array[v2] = 1;
								done = 0;
							}
						}
					}
				}
			}
	
			if( mesh == 0 )
				nv_final1 = nv_final;
			else 
				nv_final2 = nv_final;
		}
	
	
		double *rvals;
		int *edges;
		int *nedges;
		int nv;
		
		getAltCylinderMesh( &rvals, &edges, &nedges, &nv, pathLen1, join_radius, join_len ); 
	
		// discard edge-edges.
		
		// remove the edges of path1 pointing to path2 and path1
		for( int v = 0; v < pathLen1; v++ )
		{
			for( int t = 0; t < nedges[v]; )
			{
				if( edges[v*MAX_EDGES+t] >= nv-pathLen2 || edges[v*MAX_EDGES+t] < pathLen1 )
				{
					edges[v*MAX_EDGES+t] = edges[v*MAX_EDGES+nedges[v]-1];
					nedges[v]--;
				}	
				else
					t++;
			}
		}
		
		// remove the edges of path2 pointing to path1 and path2
		for( int v = nv-pathLen2; v < nv; v++ )
		{
			for( int t = 0; t < nedges[v]; )
			{
				if( edges[v*MAX_EDGES+t] < pathLen1 || edges[v*MAX_EDGES+t] >= nv-pathLen2 )
				{
					edges[v*MAX_EDGES+t] = edges[v*MAX_EDGES+nedges[v]-1];
					nedges[v]--;
				}	
				else
					t++;
			}
		}
	
		int *cyl_list = (int *)malloc(sizeof(int) * nv );
		for( int v = 0; v < nv; v++ )
			cyl_list[v] = v;
		
		// PBC wrap to the previous element of the path.
	
		for( int i = 1; i < pathLen1; i++ )
		{
			double dr[3] = { 
				r1[3*path1[i]+0] - r1[3*path1[i-1]+0],
				r1[3*path1[i]+1] - r1[3*path1[i-1]+1],
				r1[3*path1[i]+2] - r1[3*path1[i-1]+2] };
			while( dr[0] < -theSurface1->PBC_vec[0][0]/2 ) dr[0] += theSurface1->PBC_vec[0][0]; 
			while( dr[1] < -theSurface1->PBC_vec[1][1]/2 ) dr[1] += theSurface1->PBC_vec[1][1]; 
			while( dr[2] < -theSurface1->PBC_vec[2][2]/2 ) dr[2] += theSurface1->PBC_vec[2][2]; 
			
			while( dr[0] > theSurface1->PBC_vec[0][0]/2 ) dr[0] -= theSurface1->PBC_vec[0][0]; 
			while( dr[1] > theSurface1->PBC_vec[1][1]/2 ) dr[1] -= theSurface1->PBC_vec[1][1]; 
			while( dr[2] > theSurface1->PBC_vec[2][2]/2 ) dr[2] -= theSurface1->PBC_vec[2][2]; 
	
			r1[3*path1[i]+0] = r1[3*path1[i-1]+0] + dr[0];
			r1[3*path1[i]+1] = r1[3*path1[i-1]+1] + dr[1];
			r1[3*path1[i]+2] = r1[3*path1[i-1]+2] + dr[2];
		}
		
		for( int i = 1; i < pathLen2; i++ )
		{
			double dr[3] = { 
				r2[3*path2[i]+0] - r2[3*path2[i-1]+0],
				r2[3*path2[i]+1] - r2[3*path2[i-1]+1],
				r2[3*path2[i]+2] - r2[3*path2[i-1]+2] };
			while( dr[0] < -theSurface2->PBC_vec[0][0]/2 ) dr[0] += theSurface2->PBC_vec[0][0]; 
			while( dr[1] < -theSurface2->PBC_vec[1][1]/2 ) dr[1] += theSurface2->PBC_vec[1][1]; 
			while( dr[2] < -theSurface2->PBC_vec[2][2]/2 ) dr[2] += theSurface2->PBC_vec[2][2]; 
			
			while( dr[0] > theSurface2->PBC_vec[0][0]/2 ) dr[0] -= theSurface2->PBC_vec[0][0]; 
			while( dr[1] > theSurface2->PBC_vec[1][1]/2 ) dr[1] -= theSurface2->PBC_vec[1][1]; 
			while( dr[2] > theSurface2->PBC_vec[2][2]/2 ) dr[2] -= theSurface2->PBC_vec[2][2]; 
	
			r2[3*path2[i]+0] = r2[3*path2[i-1]+0] + dr[0];
			r2[3*path2[i]+1] = r2[3*path2[i-1]+1] + dr[1];
			r2[3*path2[i]+2] = r2[3*path2[i-1]+2] + dr[2];
		}
	
		alignStructuresOnAtomSet( r1, path1, rvals, cyl_list, pathLen1, nv );
		
		double cyl_com[3] = { 0,0,0};

		for( int c = 0; c < nv; c++ )
		{
			cyl_com[0] += rvals[3*c+0];
			cyl_com[1] += rvals[3*c+1];
			cyl_com[2] += rvals[3*c+2];
		}

		cyl_com[0] /= nv;
		cyl_com[1] /= nv;
		cyl_com[2] /= nv;
	
		int reverse_cyl_path2[pathLen2];
		int back_path2[nv];
	
		if( mesh2_has_pbc )
		{
			for( int p = 0; p < pathLen2; p++ )
			{
				reverse_cyl_path2[p] = nv-pathLen2+p;
				back_path2[reverse_cyl_path2[p]] = p;
			}
		}
		else
		{
			for( int p = 0; p < pathLen2; p++ )
			{
				reverse_cyl_path2[p] = nv-1-p;
				back_path2[reverse_cyl_path2[p]] = p;
			}
		}
	
		if( mesh2_has_pbc )
		{
			// find best cyclic permutation.
	
			int best_start = 0;
			double best_chi2 = 1e10;
	
			double com_path1[3]={0,0,0};
			double com_path2[3]={0,0,0};
	
			for( int v = 0; v < pathLen2; v++ )
			{
				com_path2[0] += r2[3*path2[v]+0];
				com_path2[1] += r2[3*path2[v]+1];
				com_path2[2] += r2[3*path2[v]+2];
				
				com_path1[0] += rvals[3*reverse_cyl_path2[v]+0];
				com_path1[1] += rvals[3*reverse_cyl_path2[v]+1];
				com_path1[2] += rvals[3*reverse_cyl_path2[v]+2];
			}
	
			com_path2[0] /= pathLen2;
			com_path2[1] /= pathLen2;
			com_path2[2] /= pathLen2;
			
			com_path1[0] /= pathLen2;
			com_path1[1] /= pathLen2;
			com_path1[2] /= pathLen2;
	
			for( int v = 0; v < theSurface2->nv; v++ )
			{
				r2[3*v+0] += com_path1[0] - com_path2[0];
				r2[3*v+1] += com_path1[1] - com_path2[1];
				r2[3*v+2] += com_path1[2] - com_path2[2];
			}
	
			for( int strt = 0; strt < pathLen2; strt++ )
			{
				double chi2 = 0;
				for( int p = 0; p < pathLen2; p++ )
				{
					int xp = p + strt;
					if( xp >= pathLen2 ) xp -= pathLen2;
	
					double dr[3] = { 
						rvals[3*reverse_cyl_path2[xp]+0] - (r2[3*path2[p]+0]),
						rvals[3*reverse_cyl_path2[xp]+1] - (r2[3*path2[p]+1]),
						rvals[3*reverse_cyl_path2[xp]+2] - (r2[3*path2[p]+2]) };
	
					while( dr[0] < -theSurface2->PBC_vec[0][0]/2 ) dr[0] += theSurface2->PBC_vec[0][0]; 
					while( dr[1] < -theSurface2->PBC_vec[1][1]/2 ) dr[1] += theSurface2->PBC_vec[1][1]; 
					while( dr[2] < -theSurface2->PBC_vec[2][2]/2 ) dr[2] += theSurface2->PBC_vec[2][2]; 
				
					while( dr[0] > theSurface2->PBC_vec[0][0]/2 ) dr[0] -= theSurface2->PBC_vec[0][0]; 
					while( dr[1] > theSurface2->PBC_vec[1][1]/2 ) dr[1] -= theSurface2->PBC_vec[1][1]; 
					while( dr[2] > theSurface2->PBC_vec[2][2]/2 ) dr[2] -= theSurface2->PBC_vec[2][2]; 
		
					double r2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
	
					chi2 += r2;		
				}
				printf("strt: %d chi2: %lf\n", strt, chi2 );
				if( chi2 < best_chi2 )
				{
					best_start = strt;
					best_chi2 = chi2;
				}
			}
	
			int copy[pathLen2];
			for( int p = 0; p < pathLen2; p++ )
			{
				int xp = p + best_start;
				if( xp >= pathLen2 ) xp -= pathLen2;
				copy[p]= reverse_cyl_path2[xp];
			}
			memcpy( reverse_cyl_path2, copy, sizeof(int) * pathLen2 );
	
		}
		else
		{
			alignStructuresOnAtomSet( rvals, reverse_cyl_path2, r2, path2, pathLen2, theSurface2->nv );
		}
	
	#define REDO_BACKPATH2
	
	#ifdef REDO_BACKPATH2
		if( mesh2_has_pbc )
		{
			for( int p = 0; p < pathLen2; p++ )
				back_path2[reverse_cyl_path2[p]] = p;
		}
		else
		{
			for( int p = 0; p < pathLen2; p++ )
				back_path2[reverse_cyl_path2[p]] = p;
		}
	#endif
	
	#ifdef DEBUG_ALL
		for( int i = 0; i < theSurface1->nv; i++ )
			printf("O %lf %lf %lf\n", r1[3*i+0], r1[3*i+1], r1[3*i+2] );
	
		for( int i = 0; i < nv; i++ )
		{
			printf("C %lf %lf %lf # %d", rvals[3*i+0], rvals[3*i+1], rvals[3*i+2], nedges[i] );
			for( int j = 0; j < nedges[i]; j++ )
				printf(" %d", edges[i*MAX_EDGES+j] );
			printf("\n");
		}
		
		for( int i = 0; i < theSurface2->nv; i++ )
			printf("N %lf %lf %lf\n", r2[3*i+0], r2[3*i+1], r2[3*i+2] );
	#endif
		FILE *joinFile = fopen("join.mesh","w");
		
		fprintf(joinFile, "3D saved surface\n");
		for( int t = 0; t < 3; t++ )
			fprintf(joinFile, "%lf %lf %lf\n", theSurface1->PBC_vec[t][0], theSurface1->PBC_vec[t][1], theSurface1->PBC_vec[t][2] );	
	
		int *map1 = (int *)malloc( sizeof(int) * theSurface1->nv );
		for( int t = 0; t < theSurface1->nv; t++ )
			map1[t] = -1; // just so I notice if something is wrong.
		int *map2 = (int *)malloc( sizeof(int) * theSurface2->nv );
		for( int t = 0; t < theSurface2->nv; t++ )
			map2[t] = -1; // just so I notice if something is wrong.
		int nmap1 = 0;
		int nmap2 = 0;
	
		for( int v = 0; v < theSurface1->nv; v++ )
		{
			if( vert_color1[v] == 0 || vert_color1[v] == 2 )	
			{
				map1[v] = nmap1;
				nmap1++;
			}
		}
		
		for( int v = 0; v < theSurface2->nv; v++ )
		{
			if( vert_color2[v] == 0 || vert_color2[v] == 2 )	
			{
				map2[v] = nmap2;
				nmap2++;
			}
		}
	
		int *cyl_map = (int *)malloc( sizeof(int) * nv );
		
		for( int x = 0; x < nv; x++ )
		{
			if( x < pathLen1 )
				cyl_map[x] = map1[path1[cyl_list[x]]];
			else if( x >= nv - pathLen2 )
				cyl_map[x] = map2[path2[back_path2[x]]] + nv - pathLen2;
				//cyl_map[x] = reverse_cyl_path2[(x-(nv-pathLen2))] + nv - pathLen1;
			else
				cyl_map[x] = x;
		}
	
		typedef struct
		{	
			int ind;
			double r[3];
			int valence;
			int edges[MAX_EDGES];
		} meshp;
	
		int write_to = 0;
		meshp *theMesh = (meshp *)malloc( sizeof(meshp) * (nmap1 + nmap2 + nv - pathLen1 - pathLen2 ) );
	
		printf("Writing join.mesh.\n");
	
		int cur_offset = 0;
		int total_verts = 0;
		for( int mesh = 0; mesh < 2; mesh++ )
		{
			surface *theSurface = theSurface1;
			int *color_array = vert_color1;
			int *path = path1;
			int pathLen = pathLen1;
			int *cylPath = cyl_list; 
			int *map = map1;
			int nmap = nmap1;
			double *vert_array = r1;
			double scale = 1.0;
			if( mesh == 1)
			{
				scale = mesh_scale_2;
				nmap = nmap2;
				cylPath = reverse_cyl_path2;
				theSurface = theSurface2;
				color_array = vert_color2;
				path = path2;
				pathLen = pathLen2;
				map = map2;
				vert_array = r2;
			}
	
			for( int v = 0; v < theSurface->nv; v++ )
			{
				int new_v = cur_offset + map[v];
				if( color_array[v] == 0 )	
				{
					theMesh[write_to].ind = new_v;
					theMesh[write_to].r[0] = vert_array[3*v+0]*scale;
					theMesh[write_to].r[1] = vert_array[3*v+1]*scale;
					theMesh[write_to].r[2] = vert_array[3*v+2]*scale;
					theMesh[write_to].valence = theSurface->theVertices[v].valence;
					for( int e = 0; e < theMesh[write_to].valence; e++ )
						theMesh[write_to].edges[e] = cur_offset + map[theSurface->theVertices[v].edges[e]];
					write_to++;
					
					total_verts++;	
	#if 0
					fprintf(joinFile, "%d %lf %lf %lf %d", new_v, vert_array[3*v+0]*scale, vert_array[3*v+1]*scale, vert_array[3*v+2]*scale, theSurface->theVertices[v].valence );
					for( int e = 0; e < theSurface->theVertices[v].valence; e++ )
						fprintf(joinFile, " %d", cur_offset + map[theSurface->theVertices[v].edges[e]] ); 
					fprintf(joinFile, "\n");
	#endif
				}
				else if( color_array[v] == 2 )
				{
					int new_edges[20];
					int nnew = 0;
		
					// the old vertices.
					for( int e = 0; e < theSurface->theVertices[v].valence; e++ )
					{
						int v2 = theSurface->theVertices[v].edges[e];
		
						if( color_array[v2] == 1 )
							continue;
		
	//					if( color_array[v2] == 0 || color_array[v2] == 2 )
							new_edges[nnew] = cur_offset + map[v2];
						nnew++;	
					}
	
					// the new vertices.
					
					int my_index = -1;
	
					for( int tv = 0; tv < pathLen; tv++ )
					{
						if( path[tv] == v )
							my_index = tv;
					}
						
	
					for( int e = 0; e < nedges[cylPath[my_index]]; e++ )
					{	// cylinder indices are offset by the first object, nmap1-pathLen1.
	
						new_edges[nnew] = nmap1 - pathLen1 + edges[cylPath[my_index]*MAX_EDGES+e];
						nnew++;
					}
					
					if( nnew < 5 || nnew > 7 )
					{
						printf("ERROR valence %d.\n", nnew );
	//					exit(1);
					}
	
					total_verts++;	
					
					theMesh[write_to].ind = new_v;
					theMesh[write_to].r[0] = vert_array[3*v+0]*scale;
					theMesh[write_to].r[1] = vert_array[3*v+1]*scale;
					theMesh[write_to].r[2] = vert_array[3*v+2]*scale;
					theMesh[write_to].valence = nnew;
					for( int e = 0; e < theMesh[write_to].valence; e++ )
					{
						theMesh[write_to].edges[e] = new_edges[e];
						
						
					}
					write_to++;
	#if 0 
					fprintf(joinFile, "%d %lf %lf %lf %d", new_v, vert_array[3*v+0], vert_array[3*v+1], vert_array[3*v+2], nnew );
					for( int e = 0; e < nnew; e++ )
						fprintf(joinFile, " %d",  new_edges[e] ); 
					fprintf(joinFile, "\n");
	#endif
				}
			}
	
			cur_offset += nmap; 
	
			if( mesh == 0 )
			{
				// after doing the first object, put down the cylinder.
	
				cur_offset -= pathLen1; //we lose the first pathLen1 of our points.
	
				for( int tv = pathLen1; tv < nv-pathLen2; tv++ )
				{
					total_verts++;	
					
					if( nedges[tv] < 5 || nedges[tv] > 7 )
					{
						printf("ERROR valence %d.\n", nedges[tv] );
					}
					
					theMesh[write_to].ind = cur_offset + tv;
					theMesh[write_to].r[0] = rvals[tv*3+0];
					theMesh[write_to].r[1] = rvals[tv*3+1];
					theMesh[write_to].r[2] = rvals[tv*3+2];
					theMesh[write_to].valence = nedges[tv];
	
					for( int e = 0; e < nedges[tv]; e++ )
					{
						if( edges[tv*MAX_EDGES+e] < pathLen1 )
							theMesh[write_to].edges[e] = cyl_map[edges[tv*MAX_EDGES+e]];
						else
							theMesh[write_to].edges[e] = cur_offset + cyl_map[edges[tv*MAX_EDGES+e]];
						
						if( theMesh[write_to].edges[e] == 264 && write_to == 197 )
						{
							printf("check.\n");
						}
					}
					write_to++;
					
	
	#if 0
					fprintf(joinFile, "%d %lf %lf %lf %d", cur_offset + tv, rvals[tv*3+0], rvals[3*tv+1], rvals[3*tv+2], nedges[tv] );
					
		
	
					for( int e = 0; e < nedges[tv]; e++ )
					{
						if( edges[tv*MAX_EDGES+e] < pathLen1 )
							fprintf(joinFile," %d", cyl_map[edges[tv*MAX_EDGES+e]] );
						else
							fprintf(joinFile," %d", cur_offset + cyl_map[edges[tv*MAX_EDGES+e]] );
					}
					fprintf(joinFile,"\n");
	#endif
				}
				
				cur_offset += nv - pathLen2; // advance to the end, with pathLen1 already sub'd
			}
		}
	
		int *sorter = (int *)malloc( sizeof(int) * write_to );
	
		for( int x = 0; x < write_to; x++ )
			sorter[x] = x;
	
	
		done = 0;
		while( !done )
		{
			done = 1;
			int prec[write_to];
			int npre = 0;
			for( int x = 0; x < write_to; x++ )
			{
				if( theMesh[sorter[x]].valence == 6 )
				{
					prec[npre] = x;
					npre++;
				} 
	
				if( theMesh[sorter[x]].valence != 6 && npre > 0 )
				{
					done = 0;
					int t = sorter[prec[0]];
					sorter[prec[0]] = sorter[x];
					sorter[x] = t;
					for( int p = 0; p < npre-1; p++ )	
						prec[p] = prec[p+1];
					npre--;
				}
			}
		}
		
		double com[3]={0,0,0};
		for( int v = 0; v < write_to; v++ )
		{
			com[0] += theMesh[sorter[v]].r[0];
			com[1] += theMesh[sorter[v]].r[1];
			com[2] += theMesh[sorter[v]].r[2];
		}
	
		com[0]/=write_to;
		com[1]/=write_to;
		com[2]/=write_to;
		
		double min[3]={1e10,1e10,1e10}, max[3]={-1e10,-1e10,-1e10};

		if( mesh1_has_pbc )
		{
			for( int v = 0; v < write_to; v++ )
			{
				theMesh[v].r[0] -= cyl_com[0];
				theMesh[v].r[1] -= cyl_com[1];
				theMesh[v].r[2] -= cyl_com[2];

				while( theMesh[v].r[0] < -theSurface1->PBC_vec[0][0]/2 ) theMesh[v].r[0] += theSurface1->PBC_vec[0][0];
				while( theMesh[v].r[1] < -theSurface1->PBC_vec[1][1]/2 ) theMesh[v].r[1] += theSurface1->PBC_vec[1][1];
				while( theMesh[v].r[2] < -theSurface1->PBC_vec[2][2]/2 ) theMesh[v].r[2] += theSurface1->PBC_vec[2][2];
				while( theMesh[v].r[0] > theSurface1->PBC_vec[0][0]/2 ) theMesh[v].r[0] -= theSurface1->PBC_vec[0][0];
				while( theMesh[v].r[1] > theSurface1->PBC_vec[1][1]/2 ) theMesh[v].r[1] -= theSurface1->PBC_vec[1][1];
				while( theMesh[v].r[2] > theSurface1->PBC_vec[2][2]/2 ) theMesh[v].r[2] -= theSurface1->PBC_vec[2][2];
			}	
		}
	

		for( int v = 0; v < write_to; v++ )
		{
//			 theMesh[sorter[v]].r[0]-=com[0];
//			 theMesh[sorter[v]].r[1]-=com[1];
//			 theMesh[sorter[v]].r[2]-=com[2];
	
			for( int c = 0; c < 3; c++ )
			{
				if( theMesh[sorter[v]].r[c] < min[c] )
					min[c] = theMesh[sorter[v]].r[c];
				if( theMesh[sorter[v]].r[c] > max[c] )
					max[c] = theMesh[sorter[v]].r[c];
			}
		}
	
//		printf("MIN: %le %le %le MAX %le %le %le\n", min[0], min[1], min[2], max[0], max[1], max[2] )

	
		for( int v = 0; v < write_to; v++ )
		{
			fprintf(joinFile, "%d %lf %lf %lf %d", v, 
				theMesh[sorter[v]].r[0],
				theMesh[sorter[v]].r[1],
				theMesh[sorter[v]].r[2],
				theMesh[sorter[v]].valence );
			for( int e = 0; e < theMesh[sorter[v]].valence; e++ )
				fprintf(joinFile, " %d", sorter[theMesh[sorter[v]].edges[e]] );	
			fprintf(joinFile,"\n");
		}
	
		fclose(joinFile);
	
		FILE *read = fopen("join.mesh", "r");
	
		int * edge_valence = (int *)malloc( sizeof(int) * (MAX_EDGES+1) * total_verts );
	
		getLine(read,buffer);
		getLine(read,buffer);
		getLine(read,buffer);
		getLine(read,buffer);
	
		while( !feof(read) )
		{
			getLine(read,buffer);
			if( feof(read) ) continue;
			int edges[20];
			int ind,val;
			double r[3];
			int nr = sscanf( buffer, "%d %lf %lf %lf %d %d %d %d %d %d %d %d %d",
				&ind, r+0,r+1,r+2, &val, 
					edges+0,
					edges+1,
					edges+2,
					edges+3,
					edges+4,
					edges+5,
					edges+6,
					edges+7 );
			edge_valence[ind*(MAX_EDGES+1)] = val;
			memcpy( edge_valence+ind*(MAX_EDGES+1)+1, edges, sizeof(int) * val );
		}
	
		int n_disasters = 0;
		int n_flaws = 0;
		for( int e = 0; e < total_verts; e++ )
		{
			if( edge_valence[e*(MAX_EDGES+1)+0] != 6 )
			{
				for( int ex = 0; ex < edge_valence[e*(MAX_EDGES+1)+0]; ex++ )
				{
					int e2 = edge_valence[e*(MAX_EDGES+1)+1+ex];
					if( e2 < e ) continue;
	
					if( edge_valence[e2*(MAX_EDGES+1)+0] != 6 )
						n_flaws++;
				}
			}
		}
	
		for( int e = 0; e < total_verts; e++ )
		{	
			int val = edge_valence[e*(MAX_EDGES+1)+0];
	
			for( int ex = 0; ex < edge_valence[e*(MAX_EDGES+1)+0]; ex++ )
			{
				int e2 = edge_valence[e*(MAX_EDGES+1)+1+ex];
			
				int gotit = 0;
				int val2 = edge_valence[e2*(MAX_EDGES+1)+0];
				for( int ex2 = 0; ex2 < edge_valence[e2*(MAX_EDGES+1)+0]; ex2++ )
				{
					if( edge_valence[e2*(MAX_EDGES+1)+1+ex2] == e )
						gotit=1;
				}
	
				if( !gotit )
				{
					printf("PROBLEM: vertex %d is bonded to %d but it is not reciprocated.\n", e, e2);
					n_disasters += 1;
				}
			}
		}

		free(rvals);
		free(edges);
		free(nedges);

		free(vert_color1);
		free(vert_color2);
		free(cyl_list);
		free(cyl_map);
		free(map1);
		free(map2);	
		free(edge_valence);
		free(sorter);
		free(theMesh);

		if( n_disasters > 0 )
			printf("There were %d disasters.\n", n_disasters );
		else if( n_flaws >= 1 )
			printf("There %s %d flawed %s with two irregular vertices. The mesh must be subdivided or the input to join adjusted.\n", (n_flaws > 1 ? "are" : "is" ), n_flaws, (n_flaws > 1 ? "triangles" : "triangle")  ); 
		else
			printf("This mesh can be used without subdivision.\n");	
	
		n_flaws += n_disasters;

		if( do_auto && n_flaws == 0 )
			auto_done = 1;

		if( !auto_done )
		{
			join_radius += dJoin;
			if( dJoin > 0 ) 
				dJoin += 1;
			else
				dJoin -= 1;
			dJoin *= -1;
			printf("dJoin: %lf join_radius: %lf\n", dJoin, join_radius );
		}
	}
	return 0;
}
	

