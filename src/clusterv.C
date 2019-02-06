#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mutil.h"

#include "clusterv.h"

typedef struct
{
	int nvertices;
	int site1, site2;
	int *vertices;
} ridge;
extern ridge *ridgeStorage;
extern double *vertStorage;
extern int nridges;
extern int nverts;

extern "C" int voronoi2D( double *pts, int npts, const char *unique );
extern "C" int voronoi( double *pts, int npts, const char *unique );

void get3DVoronoi( double *pts, int npts, const char *uniq, int *bonds, int *nbonds, int max_bonds, double PBC_vec[3][3], double max_cut )
{
	int npts_use = npts;
	double buffer = 20.0;

	double *temp_save = NULL; 
	int npts_orig = npts;

	double r_ext = 20.0;

	double maxr = 0;
	double midp[3] = { 0,0,0};
	for( int p = 0; p < npts; p++ )
	{
		midp[0] += pts[3*p+0];
		midp[1] += pts[3*p+1];
		midp[2] += pts[3*p+2];
	}

	midp[0] /= npts;
	midp[1] /= npts;
	midp[2] /= npts;

	for( int p = 0; p < npts; p++ )
	{
		double dr[3] = { pts[3*p+0] - midp[0], pts[3*p+1] - midp[1], pts[3*p+2] - midp[2] };
		double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

		if( r > maxr )
			maxr = r;
	}
	
	int nspace = npts + 0.2 * npts;

	double *pts_use = (double *)malloc( sizeof(double) * nspace * 3 );
	int *pt_map = (int *)malloc( sizeof(int) * nspace );
	
	for( int p = 0; p < npts; p++ )
	{
		pt_map[p] = p;

		pts_use[3*p+0] = pts[3*p+0];
		pts_use[3*p+1] = pts[3*p+1];
		pts_use[3*p+2] = pts[3*p+2];
	}
	
	npts_use = npts;

	double BUFFER = 15.0;
	
	for( int p = 0; p < npts; p++ )
	{
		nbonds[p] = 0;

		for( int dx = -1; dx <= 1; dx++ )
		for( int dy = -1; dy <= 1; dy++ )
		for( int dz = -1; dz <= 1; dz++ )
		{
			if( dx == 0 && dy == 0 && dz == 0 )
				continue;

			double newp[3] = { pts[3*p+0] + dx * PBC_vec[0][0] + dy * PBC_vec[1][0] + dz * PBC_vec[2][0],
					   pts[3*p+1] + dx * PBC_vec[0][1] + dy * PBC_vec[1][1] + dz * PBC_vec[2][1],
					   pts[3*p+2] + dx * PBC_vec[0][2] + dy * PBC_vec[1][2] + dz * PBC_vec[2][2] };
			// these points are outside of the domain.

			for( int p2 = 0; p2 < npts; p2++ )
			{
				double dr[3] = { newp[0] - pts[3*p2+0], newp[1] - pts[3*p2+1], newp[2] - pts[3*p2+2] };
				double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

				if( r < BUFFER )
				{
					if( npts_use == nspace )
					{
						nspace *= 2;
						pts_use = (double*)realloc( pts_use, sizeof(double) *3 * nspace );
						pt_map = (int*)realloc( pt_map, sizeof(int) * nspace );	
					}

					pts_use[3*npts_use+0] = newp[0];
					pts_use[3*npts_use+1] = newp[1];
					pts_use[3*npts_use+2] = newp[2];
					pt_map[npts_use] = p;

					npts_use+=1;
					break;
				}
			}
		}
	}	
	
	printf("%d\n", npts_use );
	printf("hi\n");
	for( int p = 0; p < npts_use; p++ )
	{
		if( p < npts )
			printf("C %lf %lf %lf\n", pts_use[3*p+0], pts_use[3*p+1], pts_use[3*p+2] );
		else
			printf("O %lf %lf %lf\n", pts_use[3*p+0], pts_use[3*p+1], pts_use[3*p+2] );
	}

	voronoi( pts_use, npts_use, uniq );	

	ridge *ridges = ridgeStorage;
	double *allVertices = vertStorage;
	
	 for( int x = 0; x < nridges; x++ )
	 {
	 	int nvertsl = ridges[x].nvertices;
	 	int s1 = ridges[x].site1;
	 	int s2 = ridges[x].site2;

		if( s1 >= npts && s2 >= npts )  continue;

		int p1 = pt_map[s1];
		int p2 = pt_map[s2];

		double dr[3] = { pts_use[3*s1+0] - pts_use[3*s2+0], pts_use[3*s1+1] - pts_use[3*s2+1], pts_use[3*s1+2] - pts_use[3*s2+2] };
		double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2] );

		if( r > max_cut ) continue;

		if( (p1 == 843 && p2 == 276) || ( p1 == 276 && p2 == 843) )
		{
			printf("hm.\n");
		}
		if( p1 == p2 ) continue;
		int gotit = 0;
		for( int px = 0; px < nbonds[p1]; px++ )
		{
			if( bonds[p1*max_bonds+px] == p2 ) 
				gotit =1;
		}
		if( !gotit && nbonds[p1] < max_bonds )
		{
			bonds[p1*max_bonds+nbonds[p1]] = p2;
			nbonds[p1] += 1;
		}
		gotit = 0;
		for( int px = 0; px < nbonds[p2]; px++ )
		{
			if( bonds[p2*max_bonds+px] == p1 ) 
				gotit =1;
		}
		if( !gotit && nbonds[p2] < max_bonds )
		{
			bonds[p2*max_bonds+nbonds[p2]] = p1;
			nbonds[p2] += 1;
		}

	}	

	free(pts_use); 
	free(pt_map );
}

void get2DHexVoronoiConnectivity( double *pts, int npts, const char *uniq, int *bonds, int *nbonds, int max_bonds, double BoxL[2] )
{
	int npts_use = npts;
	double buffer = 20.0;

	double *temp_save = NULL; 
	int npts_orig = npts;

	double r_ext = 20.0;
	double min[2] = { 1e10, 1e10 }, max[2] = { -1e10, -1e10 };

	double PBC_vec[2][3] = { 
			{ BoxL[0], 0, 0 },
			{ -BoxL[1]*cos(60.0*M_PI/180.0), BoxL[1]*sin(60.0*M_PI/180.0), 0 }
		};

	double maxr = 0;
	double midp[2] = { 0,0};
	for( int p = 0; p < npts; p++ )
	{
		midp[0] += pts[3*p+0];
		midp[1] += pts[3*p+1];
	}

	midp[0] /= npts;
	midp[1] /= npts;

	for( int p = 0; p < npts; p++ )
	{
		double dr[2] = { pts[3*p+0] - midp[0], pts[3*p+1] - midp[1] };
		double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]);

		if( r > maxr )
			maxr = r;
	}
	

	for( int p = 0; p < npts; p++ )
	{
		for( int dx = -1; dx <= 1; dx++ )
		for( int dy = -1; dy <= 1; dy++ )
		{
			if( dx == 0 && dy == 0 )
				continue;

			double newp[2] = { pts[3*p+0] + dx * PBC_vec[0][0] + dy * PBC_vec[1][0],
					   pts[3*p+1] + dx * PBC_vec[0][2] + dy * PBC_vec[1][1] };
			double dr[2] = { newp[0]-midp[0], newp[1]-midp[1] };
			double newr = sqrt(dr[0]*dr[0]+dr[1]*dr[1] );

			if( newr < maxr + r_ext )
				npts_use+=1;
		}
	}	
	
	double *pts_use = (double *)malloc( sizeof(double) * npts_use * 2 );
	int *pt_map = (int *)malloc( sizeof(int) * npts_use );
	
	for( int p = 0; p < npts; p++ )
	{
		pts_use[2*p+0] = pts[3*p+0];
		pts_use[2*p+1] = pts[3*p+1];
	}
	
	npts_use = npts;

	for( int p = 0; p < npts; p++ )
	{
		pt_map[p] = p;
		nbonds[p] = 0;
		for( int dx = -1; dx <= 1; dx++ )
		for( int dy = -1; dy <= 1; dy++ )
		{
			if( dx == 0 && dy == 0 )
				continue;
			double newp[2] = { pts[3*p+0] + dx * PBC_vec[0][0] + dy * PBC_vec[1][0],
					   pts[3*p+1] + dx * PBC_vec[0][2] + dy * PBC_vec[1][1] };
			double dr[2] = { newp[0]-midp[0], newp[1]-midp[1] };
			double newr = sqrt(dr[0]*dr[0]+dr[1]*dr[1] );

			if( newr < maxr + r_ext )
			{
				pt_map[npts_use] = p;
				pts_use[2*npts_use+0] = newp[0];
				pts_use[2*npts_use+1] = newp[1];
				npts_use += 1;
			}
		}
	}	

	printf("%d\n", npts_use );
	printf("hi\n");
	for( int p = 0; p < npts_use; p++ )
	{
		if( p < npts )
			printf("C %lf %lf 0.0\n", pts_use[2*p+0], pts_use[2*p+1] );
		else
			printf("O %lf %lf 0.0\n", pts_use[2*p+0], pts_use[2*p+1] );
	}

	voronoi2D( pts_use, npts_use, uniq );	

	ridge *ridges = ridgeStorage;
	double *allVertices = vertStorage;
	
	 for( int x = 0; x < nridges; x++ )
	 {
	 	int nvertsl = ridges[x].nvertices;
	 	int s1 = ridges[x].site1;
	 	int s2 = ridges[x].site2;

		if( s1 >= npts && s2 >= npts )  continue;

		int p1 = pt_map[s1];
		int p2 = pt_map[s2];


		int gotit = 0;
		for( int px = 0; px < nbonds[p1]; px++ )
		{
			if( bonds[p1*max_bonds+px] == p2 ) 
				gotit =1;
		}
		if( !gotit && nbonds[p1] < max_bonds )
		{
			bonds[p1*max_bonds+nbonds[p1]] = p2;
			nbonds[p1] += 1;
		}
		gotit = 0;
		for( int px = 0; px < nbonds[p2]; px++ )
		{
			if( bonds[p2*max_bonds+px] == p1 ) 
				gotit =1;
		}
		if( !gotit && nbonds[p2] < max_bonds )
		{
			bonds[p2*max_bonds+nbonds[p2]] = p1;
			nbonds[p2] += 1;
		}

	}	

	free(pts_use); 
	free(pt_map );
}

void get2DVoronoiConnectivity( double *pts_in, int npts, const char *uniq, int *bonds, int *nbonds, int max_bonds, double BoxL[2] )
{
	double *pts = (double *)malloc( sizeof(double) * 3 * npts );

	memcpy( pts, pts_in, sizeof(double ) * 3 * npts );

	int npts_use = npts;
	double bufferX = BoxL[0]/4;
	double bufferY = BoxL[1]/4;

	double *temp_save = NULL; 
	int npts_orig = npts;

	double min[2] = { 1e10, 1e10 }, max[2] = { -1e10, -1e10 };

	for( int p = 0; p < npts; p++ )
	{
		if( pts[3*p+0] < 0 ) pts[3*p+0] += BoxL[0];
		if( pts[3*p+1] < 0 ) pts[3*p+1] += BoxL[1];

		if( pts[3*p+0] > BoxL[0] ) pts[3*p+0] -= BoxL[0];
		if( pts[3*p+1] > BoxL[1] ) pts[3*p+1] -= BoxL[1];

		int bx = 0;
		int by = 0;

		if( pts[3*p+0] < bufferX )
			bx = 1;
		if( pts[3*p+1] < bufferY )
			by = 1;
		if( pts[3*p+0] > BoxL[0] - bufferX )
			bx = 1;
		if( pts[3*p+1] > BoxL[1] - bufferY )
			by = 1;

		if( bx && by )
			npts_use += 3;
		else if( bx )
			npts_use += 1;
		else if( by )
			npts_use += 1;
	}	
	
	double *pts_use = (double *)malloc( sizeof(double) * npts_use * 2 );
	int *pt_map = (int *)malloc( sizeof(int) * npts_use );
	
	for( int p = 0; p < npts; p++ )
	{
		pts_use[2*p+0] = pts[3*p+0];
		pts_use[2*p+1] = pts[3*p+1];
	}
	
	npts_use = npts;

	for( int p = 0; p < npts; p++ )
	{
		pt_map[p] = p;
		nbonds[p] = 0;
		int bx = 0;
		int by = 0;

		if( pts[3*p+0] < bufferX )
			bx = 1;
		if( pts[3*p+1] < bufferY )
			by = 1;
		if( pts[3*p+0] > BoxL[0] - bufferX )
			bx = -1;
		if( pts[3*p+1] > BoxL[1] - bufferY )
			by = -1;

		if( bx )
		{
			pt_map[npts_use] = p;
			pts_use[2*npts_use+0] = pts[3*p+0] + bx * BoxL[0];
			pts_use[2*npts_use+1] = pts[3*p+1];
			npts_use += 1;
		}

		if( by )
		{
			pt_map[npts_use] = p;
			pts_use[2*npts_use+0] = pts[3*p+0];
			pts_use[2*npts_use+1] = pts[3*p+1] + by * BoxL[1];
			npts_use += 1;
		}
		
		if( bx && by )
		{
			pt_map[npts_use] = p;
			pts_use[2*npts_use+0] = pts[p*3+0] + bx * BoxL[0];
			pts_use[2*npts_use+1] = pts[p*3+1] + by * BoxL[1];
			npts_use += 1;
		}

	}	

	voronoi2D( pts_use, npts_use, uniq );	

	ridge *ridges = ridgeStorage;
	double *allVertices = vertStorage;
	
	 for( int x = 0; x < nridges; x++ )
	 {
	 	int nvertsl = ridges[x].nvertices;
	 	int s1 = ridges[x].site1;
	 	int s2 = ridges[x].site2;

		if( s1 >= npts && s2 >= npts )  continue;

		int p1 = pt_map[s1];
		int p2 = pt_map[s2];


		int gotit = 0;
		for( int px = 0; px < nbonds[p1]; px++ )
		{
			if( bonds[p1*max_bonds+px] == p2 ) 
				gotit =1;
		}
		if( !gotit && nbonds[p1] < max_bonds )
		{
			bonds[p1*max_bonds+nbonds[p1]] = p2;
			nbonds[p1] += 1;
		}
		gotit = 0;
		for( int px = 0; px < nbonds[p2]; px++ )
		{
			if( bonds[p2*max_bonds+px] == p1 ) 
				gotit =1;
		}
		if( !gotit && nbonds[p2] < max_bonds )
		{
			bonds[p2*max_bonds+nbonds[p2]] = p1;
			nbonds[p2] += 1;
		}

	}	

	free(pts);
	free(pts_use); 
	free(pt_map );
}

void getSphericalVoronoiConnectivity( double *pts, int npts, const char *uniq, int *bonds, int *nbonds, int max_bonds )
{
	double *pts_use = (double *)malloc( sizeof(double) * 3 * (npts+1) );
	int npts_use = npts+1;

	memcpy( pts_use, pts, sizeof(double) * 3 * npts );

	pts_use[npts*3+0] = 0;
	pts_use[npts*3+1] = 0;
	pts_use[npts*3+2] = 0;

	voronoi( pts_use, npts_use, uniq );	

	ridge *ridges = ridgeStorage;
	double *allVertices = vertStorage;
	
	 for( int x = 0; x < nridges; x++ )
	 {
	 	int nvertsl = ridges[x].nvertices;
	 	int s1 = ridges[x].site1;
	 	int s2 = ridges[x].site2;

		if( s1 >= npts || s2 >= npts )  continue;

		int p1 = s1;
		int p2 = s2;

		int gotit = 0;
		for( int px = 0; px < nbonds[p1]; px++ )
		{
			if( bonds[p1*max_bonds+px] == p2 ) 
				gotit =1;
		}
		if( !gotit && nbonds[p1] < max_bonds )
		{
			bonds[p1*max_bonds+nbonds[p1]] = p2;
			nbonds[p1] += 1;
		}
		gotit = 0;
		for( int px = 0; px < nbonds[p2]; px++ )
		{
			if( bonds[p2*max_bonds+px] == p1 ) 
				gotit =1;
		}
		if( !gotit && nbonds[p2] < max_bonds )
		{
			bonds[p2*max_bonds+nbonds[p2]] = p1;
			nbonds[p2] += 1;
		}

	}	

	free(pts_use); 
}




double loopArea( double *pts, int npts )
{
	double area = 0;

	for( int x = npts-1; x >= 2; x-- )
		area += triangle_area( pts+0, pts+2*x, pts + 2*(x-1) );
	return area;
}

