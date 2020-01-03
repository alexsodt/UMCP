// by alex sodt
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#define TOL (1e-8)
#if 0
void rewrap( double *pts_in, double *required_pts, int *order1, int *order2, int *n1, int n2 )
{
	// required pts are inserted into the convex giftwrap.
		
	int ntot = *n1 + n2;
	double trial_pts[2*ntot];
	int order[ntot];

	for( int x = 0; x < *n1; x++ )
	{
		trial_pts[2*x+0] = pts_in[2*x+0];
		trial_pts[2*x+1] = pts_in[2*x+1];
	}

	for( int x = 0; x < n2; x++ )
	{
		trial_pts[2*(*n1+x)+0] = required_pts[2*x+0];
		trial_pts[2*(*n1+x)+1] = required_pts[2*x+1];
	}

	double phi[ntot];

	for( int x = 0; x < ntot; x++ )
		phi[x] = atan2( trial_pts[2*x+1], trial_pts[2*x+0] );

	int done = 0;

	while( !done )
	{
		int nconvex;

		giftwrap( trial_pts, order, ntot, &nconvex );

		// did this satisfy our requirement?

		int got_them[n2];
		memset( got_them, 0, sizeof(int) * n2 );

		for( int p = 0; p < nconvex; p++ )
		{
			if( order[p] >= n1 )
				got_them[order[p]-n1] = 1;
		}
	
		done = 1;

		for( int x = 0; x < n2; x++ )
		{
			if( got_them[x] == 0 )
			{
				done = 0;
			}
		}

		// remove points until we recover the points we want.
	}
	 
}
#endif

void giftwrap( double *pre_pts_in, int *ptsOrdered, int npts, int *nconvex, int expand )
{
	double pts_in[2*npts];
	int pdone[npts];
	memset( pdone, 0, sizeof(int) * npts );
	double eps = 1e-5;

	memcpy( pts_in, pre_pts_in, 2 * npts * sizeof(double) );

	if( expand )
	{
		// for points on a triangle this should give a circuit using all the points.
		
		double com[2] = {0,0};

		for( int p = 0; p < npts; p++ )
		{
			com[0] += pts_in[2*p+0]; 
			com[1] += pts_in[2*p+1]; 
		}

		com[0] /= npts;
		com[1] /= npts;
		
		for( int p = 0; p < npts; p++ )
		{
			pts_in[2*p+0] *= (1+eps); 
			pts_in[2*p+1] *= (1+eps); 
		}
	}

	if( npts == 1 )
	{
		printf("Giftwrap called with a single point.\n");
		exit(1);
	}
	double xmin = 1e10;

	int start_pt = -1;

	for( int p = 0; p < npts; p++ )
	{
		if( pts_in[2*p+0] < xmin )
		{
			xmin = pts_in[2*p+0];
			start_pt = p;
		}
	}

	if( npts == 2 )
	{
		ptsOrdered[0] = start_pt;
		ptsOrdered[1] = 1-start_pt;
		*nconvex = 2;
		return;
	}

	int newpts = 1;

	ptsOrdered[0] = start_pt;

	int use_pt = start_pt;
	double cur_pt[2] = { pts_in[2*start_pt+0], pts_in[2*start_pt+1] };
	int next_pt = -1;


	double prev_v[2] = {0,1};

	int pthit[npts];

	memset( pthit, 0, sizeof(int) * npts );

	pthit[start_pt] = 1;

	do {
		double min_dp= -1e10;

		next_pt = start_pt;

		for( int p = 0; p < npts; p++ )
		{
			if( pdone[p] ) continue;

			double dr[2] = { pts_in[2*p+0] - cur_pt[0],
					 pts_in[2*p+1] - cur_pt[1] };

			double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
			dr[0] /= r;
			dr[1] /= r;

			if( r < TOL )
				continue;

			double dp = prev_v[0] * dr[0] + prev_v[1] * dr[1]; 

			if( dp > min_dp )
			{
				min_dp = dp;
				next_pt = p;
			}
		}

		if( next_pt == start_pt )
			break;

		pdone[next_pt] = 1;	

		if( newpts >= npts )
		{
			printf("Giftwrap error 1.\n");
			exit(1);
		}
	

		prev_v[0] = pts_in[2*next_pt+0] - cur_pt[0];
		prev_v[1] = pts_in[2*next_pt+1] - cur_pt[1];

		double lpv = sqrt(prev_v[0]*prev_v[0]+prev_v[1]*prev_v[1]);
		prev_v[0] /= lpv;
		prev_v[1] /= lpv;

		if( next_pt != start_pt )
		{
			ptsOrdered[newpts] = next_pt;
			newpts++;
		}

		use_pt = next_pt;

		cur_pt[0] = pts_in[2*next_pt+0];
		cur_pt[1] = pts_in[2*next_pt+1];

		if( pthit[next_pt] && (start_pt!=next_pt) )
		{
			printf("Giftwrap error 2.\n");
			exit(1);
		}
		pthit[start_pt] = 1;

	} while (next_pt != start_pt );

	*nconvex = newpts; 
}

void angle_sort( double *pts_in, int *ptsOrdered, int npts, int *nconvex)
{
	double xmin = 1e10;

	int start_pt = -1;

	for( int p = 1; p < npts; p++ )
	{
		if( pts_in[2*p+0] < xmin )
		{
			xmin = pts_in[2*p+0];
			start_pt = p;
		}
	}

	double thetas[npts];
	for( int p = 1; p < npts; p++ )
		thetas[p-1] = atan2( pts_in[2*p+0] - pts_in[0], pts_in[2*p+1] - pts_in[1] );

	int sorter[npts];

	for( int p = 0; p < npts-1; p++ )
		sorter[p] = p;

	int done = 0;

	while( !done )
	{
		done = 1;
		
		for( int s = 0; s < npts-2; s++ )
		{
			if( thetas[sorter[s]] > thetas[sorter[s+1]] )
			{
				int t = sorter[s];
				sorter[s] = sorter[s+1];
				sorter[s+1] = t;
				done = 0;
			}
		}
	}

	ptsOrdered[0] = sorter[0]+1;

	int tso = 1;
	for( int x = 1; x < npts-1; x++ )
	{
		if( fabs(thetas[sorter[x]]-thetas[sorter[tso-1]]) < 1e-8 )
			continue;

		ptsOrdered[tso] = sorter[x]+1;
		tso++;
	}

	*nconvex = tso; 
}


