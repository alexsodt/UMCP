#define __global_boxingc__

#include "global_boxing.h"
#include <stdlib.h>
#include <math.h>

//global_boxing *boxing = NULL;
//int global_boxing_init = 0;

void global_boxing::setup_boxing( double target_box_width, double PBC_in[3][3] )
{
	for( int i = 0; i < 3; i++ )
	for( int j = 0; j < 3; j++ )
		PBC[i][j] = PBC_in[i][j];
	
	nx = ceil(PBC[0][0] / target_box_width);	
	ny = ceil(PBC[1][1] / target_box_width);	
	nz = ceil(PBC[2][2] / target_box_width);	


	if( nx < 0 ) nx = 1;
	if( ny < 0 ) ny = 1;
	if( nz < 0 ) nz = 1;

	int tot_boxes = nx * ny * nz;

	if( tot_boxes > 1e6 )
	{
		if( nx > 100 ) nx = 100;
		if( ny > 100 ) ny = 100;
		if( nz > 100 ) nz = 100;
	}
	
	tot_boxes = nx * ny * nz;

	boxes = (box *)malloc( sizeof(box) * tot_boxes );

	for( int b = 0; b < tot_boxes; b++ )
	{
		boxes[b].npspace = 1;
		boxes[b].plist = (int *)malloc( sizeof(int) * boxes[b].npspace );
	}

	clearBoxing();
}

int global_boxing::getNearPts( double *r, int *plist, double rad_search )
{
	double tr[3] = { r[0], r[1],r[2] };

	while( tr[0] < 0 ) tr[0] += PBC[0][0]; 
	while( tr[1] < 0 ) tr[1] += PBC[1][1]; 
	while( tr[2] < 0 ) tr[2] += PBC[2][2]; 
	while( tr[0] > PBC[0][0] ) tr[0] -= PBC[0][0]; 
	while( tr[1] > PBC[1][1] ) tr[1] -= PBC[1][1]; 
	while( tr[2] > PBC[2][2] ) tr[2] -= PBC[2][2]; 

	int bx = nx * (tr[0] / PBC[0][0]);
	if( bx >= nx ) bx -= nx;
	int by = ny * (tr[1] / PBC[1][1]);
	if( by >= ny ) by -= ny;
	int bz = nz * (tr[2] / PBC[2][2]);
	if( bz >= nz ) bz -= nz;

	double boxes_per_x_search =  rad_search/(PBC[0][0] / nx);
	double boxes_per_y_search =  rad_search/(PBC[1][1] / ny);
	double boxes_per_z_search =  rad_search/(PBC[2][2] / nz);

	int dbx = 1 + floor(boxes_per_x_search);
	int dby = 1 + floor(boxes_per_y_search);
	int dbz = 1 + floor(boxes_per_z_search);

	int min_x = bx - dbx;
	int min_y = by - dby;
	int min_z = bz - dbz;

	int max_x = bx + dbx;
	int max_y = by + dby;
	int max_z = bz + dbz;

	if( max_x - min_x >= nx ) { min_x = 0; max_x = nx-1; }
	if( max_y - min_y >= ny ) { min_y = 0; max_y = ny-1; }
	if( max_z - min_z >= nz ) { min_z = 0; max_z = nz-1; }

	int np = 0;

	for( int ax = min_x; ax <= max_x; ax++ )
	for( int ay = min_y; ay <= max_y; ay++ )
	for( int az = min_z; az <= max_z; az++ )
	{
		int tx = ax; if( tx < 0 ) tx += nx; if( tx >= nx) tx -= nx;
		int ty = ay; if( ty < 0 ) ty += ny; if( ty >= ny) ty -= ny;
		int tz = az; if( tz < 0 ) tz += nz; if( tz >= nz) tz -= nz;

		int box = (tx * ny + ty) * nz + tz;	

		for( int tp = 0; tp <boxes[box].np; tp++ )
		{
			plist[np] = boxes[box].plist[tp];
			np++;
		}	
	}

	return np;
}

void global_boxing::setPBC( double PBC_vec[3][3], double *alphas )
{
	PBC[0][0] = PBC_vec[0][0] * alphas[0];
	PBC[0][1] = PBC_vec[0][1] * alphas[1];
	PBC[0][2] = PBC_vec[0][2] * alphas[2];
                                
	PBC[1][0] = PBC_vec[1][0] * alphas[0];
	PBC[1][1] = PBC_vec[1][1] * alphas[1];
	PBC[1][2] = PBC_vec[1][2] * alphas[2];
                                
	PBC[2][0] = PBC_vec[2][0] * alphas[0];
	PBC[2][1] = PBC_vec[2][1] * alphas[1];
	PBC[2][2] = PBC_vec[2][2] * alphas[2];
}

int global_boxing::addp( double *r, int id )
{
	double tr[3] = { r[0], r[1],r[2] };

	while( tr[0] < 0 ) tr[0] += PBC[0][0]; 
	while( tr[1] < 0 ) tr[1] += PBC[1][1]; 
	while( tr[2] < 0 ) tr[2] += PBC[2][2]; 
	while( tr[0] > PBC[0][0] ) tr[0] -= PBC[0][0]; 
	while( tr[1] > PBC[1][1] ) tr[1] -= PBC[1][1]; 
	while( tr[2] > PBC[2][2] ) tr[2] -= PBC[2][2]; 

	int bx = nx * (tr[0] / PBC[0][0]);
	if( bx >= nx ) bx -= nx;
	int by = ny * (tr[1] / PBC[1][1]);
	if( by >= ny ) by -= ny;
	int bz = nz * (tr[2] / PBC[2][2]);
	if( bz >= nz ) bz -= nz;

	int box = (bx * ny + by) * nz + bz;
	
	if( boxes[box].np >= boxes[box].npspace )
	{
		boxes[box].npspace *= 2;
		boxes[box].plist = (int *)realloc( boxes[box].plist, sizeof(int) * boxes[box].npspace );
	}
	
	boxes[box].plist[boxes[box].np] = id;
	boxes[box].np++;

	return box;
}



void global_boxing::clearBoxing( int *list, int nclear)
{
	if( list )
	{
		for( int x = 0; x < nclear; x++ )
			boxes[list[x]].np = 0;	
	}
	else
	{
	for( int bx =0; bx < nx; bx++ )
	for( int by =0; by < ny; by++ )
	for( int bz =0; bz < nz; bz++ )
		boxes[bx*ny*nz+by*nz+bz].np = 0;	
	}
}
