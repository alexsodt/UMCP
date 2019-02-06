#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "box.h"
#include "mutil.h"

#define DEL_2
//#define SPHERE_DEBUG
#define DEBUG


typedef struct 
{
	int np;
	int ns;
	int *plist;
} bin;

int placeOnAxis( int particle, bin *bins, int nx, int ny, double t1, double t2 );


static double cutoff = 1.0;
static double L = 1.0;
static int np = 1;
static int type = BOX_SPHERE;

static int nx = 1;
static int ny = 1;
static int nz = 1;

static int naxes = 1;

static double BoxLengthA;
static double BoxLengthB;
static double BoxLengthC;
static double R;
static double **axes;
static bin **bins;
static int *axis_for_particle;
static int *box_for_particle;

void initializeBox(  int type_in, int nparticles, double R_in, double Lx, double Ly, double Lz,  double cutoff_in )
{

	type = type_in;
	np = nparticles;	
	R = R_in;
	BoxLengthA = Lx;
	BoxLengthB = Ly;
	BoxLengthC = Lz;
	cutoff = cutoff_in;
	double cscale=1.0;

	if( type == BOX_SPHERE)
	{
		cscale = 1.25;
		double circ = 2 * M_PI * R;
		double d = 2 * R;

		printf("d: %lf\n", d );
		nx = d / (cscale*cutoff) - 1; //"z"
		ny = circ / (cscale*cutoff) - 1;
	
		if( nx < 1 ) nx = 1;
		if( ny < 1 ) ny = 1;
	}
	else if( type == BOX_CYLINDER )
	{
		double d = 2 * R;
		double circ = 2 * M_PI * R;

		nx = BoxLengthA / (cscale*cutoff) - 1;
		ny = circ / (cscale*cutoff) - 1;
	
		if( nx < 1 ) nx = 1;
		if( ny < 1 ) ny = 1;
	}
	else if( type == BOX_PLANE )
	{
		cscale = 1.0;
		nx = BoxLengthA / (cscale*cutoff) - 1;
		ny = BoxLengthB / (cscale*cutoff) - 1;
	
		if( nx < 1 ) nx = 1;
		if( ny < 1 ) ny = 1;
	}
	else if( type == BOX_3D )
	{
		cscale = 1.0;
		nx = BoxLengthA / (cscale*cutoff);
		ny = BoxLengthB / (cscale*cutoff);
		nz = BoxLengthC / (cscale*cutoff);
	
		if( nx < 1 ) nx = 1;
		if( ny < 1 ) ny = 1;
		if( nz < 1 ) nz = 1;
	}
	printf("nx: %d ny: %d nz: %d\n", nx, ny, nz );

	if( type == BOX_3D )
		naxes = 1;
	else if( type == BOX_SPHERE )
#ifdef SPHERE_DEBUG
		naxes = 1;
#else
		naxes = 6;
#endif		
	bins = (bin **)malloc( sizeof(bin *) * naxes );
	axes = (double **)malloc( sizeof(double *) * naxes ); 	
	
	for( int a = 0; a < naxes; a++ )
	{
		bins[a] = (bin *)malloc( sizeof( bin) * nx * ny * nz  );
		for( int b = 0; b < nx*ny*nz; b++ )
		{
			bins[a][b].np = 0;
			bins[a][b].ns = 10;
			bins[a][b].plist = (int *)malloc( sizeof(int) * bins[a][b].ns );
		}
		axes[a] = (double *)malloc( sizeof(double) * 9 ); 
	}

	if( type == BOX_SPHERE )
	{
		axes[0][0] = 1.0;
		axes[0][1] = 0.0;
		axes[0][2] = 0.0;
	
#ifndef SPHERE_DEBUG	
		axes[1][0] = 0.0;
		axes[1][1] = 1.0;
		axes[1][2] = 0.0;
		
		axes[2][0] = 0.0;
		axes[2][1] = 0.0;
		axes[2][2] = 1.0;
		
		axes[3][0] = 0.0;
		axes[3][1] = 1.0/sqrt(2.0);
		axes[3][2] = 1.0/sqrt(2.0);
		
		axes[4][0] = 1.0/sqrt(2.0);
		axes[4][1] = 0.0;
		axes[4][2] = 1.0/sqrt(2.0);
		
		axes[5][0] = 1.0/sqrt(2.0);
		axes[5][1] = 1.0/sqrt(2.0);
		axes[5][2] = 0.0;
#endif		
			printf("naxes: %d\n", naxes );
	}
	else
	{
		axes[0][0] = 0;
		axes[0][1] = 0;
		axes[0][2] = 1;
	}
	
	if( type == BOX_SPHERE || type == BOX_CYLINDER )
	{
		for( int a = 0; a < naxes; a++ )
		{
			double axis1[3] = { rand()/(double)RAND_MAX, rand()/(double)RAND_MAX, rand()/(double)RAND_MAX };

			double dp = axis1[0] * axes[a][0] + axis1[1] * axes[a][1] + axis1[2] * axes[a][2];

			axis1[0] -= dp * axes[a][0];
			axis1[1] -= dp * axes[a][1];
			axis1[2] -= dp * axes[a][2];

			double axis2[3];

			axis2[0] =  (axis1[1] * axes[a][2] - axis1[2] * axes[a][1]);
			axis2[1] = -(axis1[0] * axes[a][2] - axis1[2] * axes[a][0]);
			axis2[2] =  (axis1[0] * axes[a][1] - axis1[1] * axes[a][0]);

			normalize(axis1);
			normalize(axis2);

			memcpy( axes[a] + 3, axis1, sizeof(double) * 3 );
			memcpy( axes[a] + 6, axis2, sizeof(double) * 3 );
		}
	}

	axis_for_particle = (int *)malloc( sizeof(int) * nparticles );
	box_for_particle = (int *)malloc( sizeof(int) * nparticles * naxes );
	
	for( int p = 0; p < nparticles; p++ )
	{
		axis_for_particle[p] = -1;

		for( int a = 0; a < naxes; a++ )
			box_for_particle[p*naxes+a] = -1;
	}
	
}

void updateBox( int particle, double *r )
{
	removeBox( particle );
	place( particle, r );
}

void removeBox( int particle )
{
	for( int a = 0; a < naxes; a++ )
	{
		int b = box_for_particle[particle*naxes+a];

		int gotit = 0;

		for( int px = 0; px < bins[a][b].np; px++ )
		{
			if( bins[a][b].plist[px] == particle )
			{
				bins[a][b].plist[px] = bins[a][b].plist[bins[a][b].np-1];
				bins[a][b].np -= 1;
				gotit = 1;
				break;
			}
		}
#ifdef DEBUG
		if( !gotit )
		{
			printf("particle fetch error.\n");
			exit(1);
		}
#endif
	}
}

int placeOnAxis( int particle, bin *bins, int nx, int ny, int nz, double t1, double t2, double t3 )
{
	while( t1 < 0 )
		t1 += 1;
	while( t1 >= 1 )
		t1 -= 1;
	
	while( t2 < 0 )
		t2 += 1;
	while( t2 >= 1 )
		t2 -= 1;
	
	while( t3 < 0 )
		t3 += 1;
	while( t3 >= 1 )
		t3 -= 1;

	int bx = nx * t1;
	int by = ny * t2;
	int bz = nz * t3;

	while( bx < 0 ) bx += nx;
	while( by < 0 ) by += ny;
	while( bz < 0 ) bz += nz;
	while( bx >= nx ) bx -= nx;
	while( by >= ny ) by -= ny;
	while( bz >= nz ) bz -= nz;

	int bin = bx*ny*nz+by*nz+bz;
	if( bins[bin].ns == bins[bin].np )
	{
		bins[bin].ns *= 2;
		bins[bin].plist = (int *)realloc( bins[bin].plist, sizeof(int) * bins[bin].ns );
	}
	
	bins[bin].plist[bins[bin].np] = particle;
	bins[bin].np += 1;

	return bin;
}

int getNear( double *r, int *neighbors, int *del)
{
	int nn = 0;

		
	double t1 = r[0] / BoxLengthA;
	double t2 = r[1] / BoxLengthB;
	double t3 = r[2] / BoxLengthC;
			
	while( t1 < 0 )
		t1 += 1;
	while( t1 >= 1 )
		t1 -= 1;
	
	while( t2 < 0 )
		t2 += 1;
	while( t2 >= 1 )
		t2 -= 1;
	
	while( t3 < 0 )
		t3 += 1;
	while( t3 >= 1 )
		t3 -= 1;

	int bx = nx * t1;
	int by = ny * t2;
	int bz = nz * t3;

	while( bx < 0 ) bx += nx;
	while( by < 0 ) by += ny;
	while( bz < 0 ) bz += nz;
	while( bx >= nx ) bx -= nx;
	while( by >= ny ) by -= ny;
	while( bz >= nz ) bz -= nz;

	int b = bx*ny*nz+by*nz+bz;

	int dsearch = 1;

#ifdef DEL_2	
	int loopx[5] = { bx-2,bx-1, bx, bx+1,bx+2 };
	int loopy[5] = { by-2,by-1, by, by+1,by+2 };
	int loopz[5] = { bz-2,bz-1, bz, bz+1,bz+2 };
	int nlx = 5;
	int nly = 5;
	int nlz = 5;
#else
	int loopx[5] = { bx-1, bx, bx+1 };
	int loopy[5] = { by-1, by, by+1 };
	int loopz[5] = { bz-1, bz, bz+1 };
	int nlx = 3;
	int nly = 3;
	int nlz = 3;
#endif

	if( nx == 2 )
	{
		loopx[0] = bx;
		loopx[1] = bx+1;
		nlx = 2;
	}
	else if( nx == 1 )
	{
		loopx[0] = bx;
		nlx = 1;
	}
	
	if( ny == 2 )
	{
		loopy[0] = by;
		loopy[1] = by+1;
		nly = 2;
	}
	else if( ny == 1 )
	{
		loopy[0] = by;
		nly = 1;
	}
	
	if( nz == 2 )
	{
		loopz[0] = bz;
		loopz[1] = bz+1;
		nly = 2;
	}
	else if( nz == 1 )
	{
		loopy[0] = bz;
		nlz = 1;
	}

	for( int ix = 0; ix < nlx; ix++ )
	{
		while( loopx[ix] < 0 ) loopx[ix] += nx;
		while( loopx[ix] >= nx ) loopx[ix] -= nx;

		for( int ix2 = 0; ix2 < ix; ix2++ )
		{
			if( loopx[ix2] == loopx[ix] )
			{
				nlx = ix;
				break;
			}
		}
	}
	
	for( int iy = 0; iy < nly; iy++ )
	{
		while( loopy[iy] < 0 ) loopy[iy] += ny;
		while( loopy[iy] >= ny ) loopy[iy] -= ny;

		for( int iy2 = 0; iy2 < iy; iy2++ )
		{
			if( loopy[iy2] == loopy[iy] )
			{
				nly = iy;
				break;
			}
		}
	}
	
	for( int iz = 0; iz < nlz; iz++ )
	{
		while( loopz[iz] < 0 ) loopz[iz] += nz;
		while( loopz[iz] >= nz ) loopz[iz] -= nz;

		for( int iz2 = 0; iz2 < iz; iz2++ )
		{
			if( loopz[iz2] == loopz[iz] )
			{
				nlz = iz;
				break;
			}
		}
	}

	for( int ix = 0; ix < nlx; ix++ )
	for( int iy = 0; iy < nly; iy++ )
	for( int iz = 0; iz < nlz; iz++ )
	{
		int tx = loopx[ix];
		int ty = loopy[iy];
		int tz = loopz[iz];
		int cur_del[3] = {0,0,0};

		int b2 = tx * ny*nz + ty*nz +tz;

		for( int tp = 0; tp < bins[0][b2].np; tp++ )
		{
			int p2 = bins[0][b2].plist[tp];

			neighbors[nn] = p2;

			if( del )
			{
				del[3*nn+0] = cur_del[0];
				del[3*nn+1] = cur_del[1];
				del[3*nn+2] = cur_del[2];
			}

			nn++;
		}
	}

	return nn;
}

int getNear( int particle, int *neighbors )
{
	int nn = 0;
	int a = axis_for_particle[particle];
	int b = box_for_particle[particle*naxes+a];

	int bx = b / (ny*nz);
	int by = (b/nz) % ny;
	int bz = b % nz;

	int dsearch = 1;

#ifdef DEL_2	
	int loopx[5] = { bx-2,bx-1, bx, bx+1,bx+2 };
	int loopy[5] = { by-2,by-1, by, by+1,by+2 };
	int loopz[5] = { bz-2,bz-1, bz, bz+1,bz+2 };
	int nlx = 5;
	int nly = 5;
	int nlz = 5;
#else
	int loopx[5] = { bx-1, bx, bx+1 };
	int loopy[5] = { by-1, by, by+1 };
	int loopz[5] = { bz-1, bz, bz+1 };
	int nlx = 3;
	int nly = 3;
	int nlz = 3;
#endif
	if( ny > 10 && (type == BOX_SPHERE && (bx <= 0.8*(nx/2.) || bx >= (nx/2.+0.2*nx/2.) )) )
	{
		dsearch = 2;
		loopy[0] = by-2;
		loopy[1] = by-1;
		loopy[2] = by-0;
		loopy[3] = by+1;
		loopy[4] = by+2;
		nly = 5;
	}

	if( nx == 2 )
	{
		loopx[0] = bx;
		loopx[1] = bx+1;
		nlx = 2;
	}
	else if( nx == 1 )
	{
		loopx[0] = bx;
		nlx = 1;
	}
	
	if( ny == 2 )
	{
		loopy[0] = by;
		loopy[1] = by+1;
		nly = 2;
	}
	else if( ny == 1 )
	{
		loopy[0] = by;
		nly = 1;
	}
	
	if( nz == 2 )
	{
		loopz[0] = bz;
		loopz[1] = bz+1;
		nly = 2;
	}
	else if( nz == 1 )
	{
		loopy[0] = bz;
		nlz = 1;
	}

	for( int ix = 0; ix < nlx; ix++ )
	for( int iy = 0; iy < nly; iy++ )
	for( int iz = 0; iz < nlz; iz++ )
	{
		int tx = loopx[ix];
		int ty = loopy[iy];
		int tz = loopz[iz];

		if( type == BOX_SPHERE )
		{
			if( tx < 0 || tx >= nx ) continue;
			
			while( ty < 0 ) ty += ny;
			while( ty >= ny ) ty -= ny;
		}
		else if( type == BOX_CYLINDER || type == BOX_PLANE )
		{
			while( tx < 0 ) tx += nx;
			while( tx >= nx ) tx -= nx;
			while( ty < 0 ) ty += ny;
			while( ty >= ny ) ty -= ny;
		}
		else if( type == BOX_3D )
		{
			while( tx < 0 ) tx += nx;
			while( tx >= nx ) tx -= nx;
			while( ty < 0 ) ty += ny;
			while( ty >= ny ) ty -= ny;
			while( tz < 0 ) tz += nz;
			while( tz >= nz ) tz -= nz;
		}
		int b2 = tx * ny*nz + ty*nz +bz;

		for( int tp = 0; tp < bins[a][b2].np; tp++ )
		{
			int p2 = bins[a][b2].plist[tp];

			if( p2 != particle )
			{
				neighbors[nn] = p2;
				nn++;
			}
		}
	}

	return nn;
}

void place( int particle, double *r )
{
	double min_axis = 1e10;
	int best_axis = -1;

	if( type == BOX_3D )
	{
		double t1 = r[0] / BoxLengthA;
		double t2 = r[1] / BoxLengthB;
		double t3 = r[2] / BoxLengthC;
			
		axis_for_particle[particle] = 0;
		box_for_particle[particle] = placeOnAxis( particle, bins[0], nx, ny, nz, t1, t2,t3 );
	}/*
	else if( type == BOX_SPHERE || type == BOX_CYLINDER )
	{
		double aproj[naxes*3];

		for( int a = 0; a < naxes; a++ )
		{
			aproj[a*3+0] = r[0] * axes[a][0] + r[1] * axes[a][1] + r[2] * axes[a][2];
			aproj[a*3+1] = r[0] * axes[a][3] + r[1] * axes[a][4] + r[2] * axes[a][5];
			aproj[a*3+2] = r[0] * axes[a][6] + r[1] * axes[a][7] + r[2] * axes[a][8];
	
			if( fabs(aproj[3*a+0]) < min_axis )
			{
				min_axis = fabs(aproj[3*a+0]);
				best_axis = a;
			}
		}

		axis_for_particle[particle] = best_axis;
		
		for( int a = 0; a < naxes; a++ )
		{	
			double t1, t2, t3=0;

				if( type == BOX_SPHERE )
					t1 = aproj[a*3] / (2*R) + 0.5;
				else if( type == BOX_CYLINDER )
					t1 = aproj[a*3] / BoxLengthA;
				else if( type == BOX_3D )
					t1 = aproj[a*3]
			
				double t2 = (M_PI + atan2( aproj[a*3+2], aproj[a*3+1] )) / (2*M_PI);						
		
				if( !(t1 < 0 || t1 > -1 ) || !(t2 <0 || t2 > -1) )
				{
					printf("r: %lf %lf %lf\n", r[0], r[1], r[2] );
				}
			box_for_particle[particle*naxes+a] = placeOnAxis( particle, bins[a], nx, ny, t1, t2 );
//			printf("Placing particle %d on axis %d in box %d.\n",
//				particle, a, box_for_particle[particle*naxes+a] );
//			printf("box: %d\n", box_for_particle[particle*naxes+a] );
		}
	}
	else
	{
		double t1 = r[0] / BoxLengthA;
		double t2 = r[1] / BoxLengthB;
			
		axis_for_particle[particle] = 0;
		box_for_particle[particle] = placeOnAxis( particle, bins[0], nx, ny, t1, t2 );

	}*/
}

void report( int p1, int p2 )
{
/*	int a1 = axis_for_particle[p1];
	int a2 = axis_for_particle[p2];

	int bx11 = box_for_particle[p1*naxes+a1] / ny;
	int by11 = box_for_particle[p1*naxes+a1] % ny;
	int bx21 = box_for_particle[p2*naxes+a1] / ny;
	int by21 = box_for_particle[p2*naxes+a1] % ny;
	int bx12 = box_for_particle[p1*naxes+a2] / ny;
	int by12 = box_for_particle[p1*naxes+a2] % ny;
	int bx22 = box_for_particle[p2*naxes+a2] / ny;
	int by22 = box_for_particle[p2*naxes+a2] % ny;
	printf("p1: %d axis1: %d [%d %d] p2: %d [%d %d]\n",
		p1, axis_for_particle[p1], bx11, by11, p2,  bx21, by21 );
	printf("p1: %d axis2: %d [%d %d] p2: %d [%d %d]\n",
		p1, axis_for_particle[p2], bx12, by12, p2,  bx22, by22 );
*/
}






