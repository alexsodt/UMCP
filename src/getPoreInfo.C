#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "util.h"

int main( int argc, char **argv )
{
	if( argc != 2 )
	{
		printf("Syntax: getPoreInfo file.rho\n");
		exit(1);
	}

	FILE *rhoFile = fopen(argv[1], "r");
	char *buffer = (char *)malloc( sizeof(char) * 100000 );
	if( !rhoFile )
	{
		printf("Cannot open density file '%s'.\n", argv[1] );
		exit(1);
	}

	double La,Lb,Lc;
	int nx, ny, nz;
	
	fscanf(rhoFile, "%lf %lf %lf\n", &La, &Lb, &Lc );
	fscanf(rhoFile, "%d %d %d\n", &nx, &ny, &nz );

	double *rho = (double *)malloc( sizeof(double) * nx * ny * nz );
	double *rho_xy = (double *)malloc( sizeof(double) * nx * ny * nz );
	double *rho_z = (double *)malloc( sizeof(double) * nz );

	for( int ix = 0; ix < nx; ix++ )
	for( int iy = 0; iy < ny; iy++ )
	{
		getLine( rhoFile, buffer );

		readNDoubles( buffer, rho+ix*ny*nz+iy*nz, nz );

		rho_xy[ix*ny+iy] = 0;
		for( int iz = 0; iz < nz; iz++ )
			 rho_xy[ix*ny+iy] += rho[ix*ny*nz+iy*nz+iz];
	}
	
	for( int iz = 0; iz < nz; iz++ )
	{
		for( int ix = 0; ix < nx; ix++ )
		for( int iy = 0; iy < ny; iy++ )
			rho_z[iz] += rho[ix*ny*nz+iy*nz+iz];
//		printf("%d %lf\n", iz, rho_z[iz] );
	}
	int i_max_z = 0;
	double max_z = rho_z[0];
	int i_min_z = 0;
	double min_z = rho_z[0];
	for( int iz = 0; iz  <nz; iz++ )
	{
		if( rho_z[iz] > max_z )
		{
			max_z = rho_z[iz];
			i_max_z = iz;
		}
		if( rho_z[iz] < min_z )
		{
			min_z = rho_z[iz];
			i_min_z = iz;
		}
	}

	int trigger = 0;
	double leaflets[4];
	double trigger_value = 0;
	int ileaflet = 0;

	for( int iz = 0; iz < nz; iz++ )
	{
		int use_z = iz + i_min_z;
		if( use_z >= nz )
			use_z -= nz;	
		// if we find a big peak count it until it comes down again.
		if( !trigger && rho_z[use_z] > 0.5 * max_z ) 
		{
			trigger = 1;
			trigger_value = rho_z[use_z];
			leaflets[ileaflet] = use_z;
		}

		if( trigger && rho_z[use_z] > trigger_value )
		{
			trigger_value = rho_z[use_z];
			leaflets[ileaflet] = use_z;
		}

		if( trigger && rho_z[use_z] < 0.5 * trigger_value )
		{
			trigger = 0;
			ileaflet++;
		}
	
		if( ileaflet == 4 ) break;
	}

	double lz[4] =
	{
		Lc * (leaflets[0]+0.5) / nz,
		Lc * (leaflets[1]+0.5) / nz,
		Lc * (leaflets[2]+0.5) / nz,
		Lc * (leaflets[3]+0.5) / nz
	};
	

	for( int dz = 1; dz < 4; dz++ )	
	{
		while( lz[dz] - lz[dz-1] < -Lc/2 )
		{
			for( int ddz = dz; ddz < 4; ddz++ )
				lz[ddz] += Lc;
		}
		
		while( lz[dz] - lz[dz-1] > Lc/2 )
		{
			for( int ddz = dz; ddz < 4; ddz++ )
				lz[ddz] -= Lc;
		}
	}	

//	printf("lz2: %le %le %le %le\n", lz[0], lz[1], lz[2], lz[3] );

	double thickness = fabs((lz[1]+lz[0])/2-(lz[3]+lz[2])/2);

	fclose(rhoFile);
	double best_chi2 = 1e10;
	double wrapto = 0;
	int nbins = nz;
	
	for( int zb = 0; zb < nz; zb++ )
	{
	        double zv = Lc * (zb+0.5) / (double)nz;
	
	         int zlow  = zb- nz/2;
	         int zhigh = zlow + nz;
	
	         double lchi2 = 0;
	         for( int iz = zlow; iz < zhigh; iz++ )
	         {
	                 double dz = Lc * (iz+0.5) / nz - zv;
	
	                 int iiz = iz;
	                 while( iiz < 0 ) iiz += nz;
	                 while( iiz >= nz ) iiz -= nz;
	
	                 lchi2 += rho_z[iiz] * (dz) * (dz);
	         }
	
	         if( lchi2 < best_chi2 )
	         {
	                 best_chi2 = lchi2;
	                 wrapto = zv;
	         }
	}
	
	double com[3] = {0,0,0};
	double ncom[3] = {0,0,0};

	com[2] = wrapto;
	ncom[2] = 1;

	// find the xy max
	
	double max_val = rho_xy[0];
	double min_val = rho_xy[0];
	double av_val  = 0;
	int min_xy = 0;
	int min_x = 0;
	int min_y = 0;

	for( int ix = 0; ix < nx; ix++ )
	for( int iy = 0; iy < ny; iy++ )
	{
		if( rho_xy[ix*ny+iy] > max_val )
			max_val = rho_xy[ix*ny+iy];
		if( rho_xy[ix*ny+iy] < min_val )
		{
			min_val = rho_xy[ix*ny+iy];
			min_x = ix;
			min_y = iy;
		}
		av_val += rho_xy[ix*ny+iy];
	}		

	com[0] = (min_x+0.5) * La / nx;
	com[1] = (min_y+0.5) * Lb / ny;
	ncom[0] += 1;
	ncom[1] += 1;

	av_val /= nx*ny;
	
	int *contig = (int *)malloc( sizeof(int) * nx * ny );
	memset( contig, 0, sizeof(int) * nx *ny );
	int *list = (int *)malloc( sizeof(int) * nx * ny );
	int nlist = 0;
	list[nlist++] = min_x*ny+min_y;
	contig[min_x*ny+min_y] = 1;

	double thresh = min_val + (max_val-min_val) * 0.025; 

	int done = 0;
	

	while( !done )
	{
		done = 1;
		
		for( int il = 0; il < nlist; il++ )
		{
			int t_x = list[il] / ny;
			int t_y = list[il] % ny;

			for( int dx = -1; dx <= 1; dx++ )
			for( int dy = -1; dy <= 1; dy++ )
			{
				int n_x = t_x+dx;
				int n_y = t_y+dy;

				if( n_x < 0 ) n_x += nx;
				if( n_x >= nx) n_x -= nx;
				if( n_y < 0 ) n_y += ny;
				if( n_y >= ny ) n_y -= ny;

				if( contig[n_x*ny+n_y] ) continue;

				if( rho_xy[n_x*ny+n_y] > thresh ) continue;

				list[nlist] = n_x*ny+n_y;
				nlist++;

				contig[n_x*ny+n_y] = 1;

				done = 0;

				int del[2] = { n_x - min_x, n_y - min_y };

				if( del[0] < -nx/2 ) del[0] += nx;
				if( del[0] >  nx/2 ) del[0] -= nx;
				if( del[1] < -ny/2 ) del[1] += ny;
				if( del[1] >  ny/2 ) del[1] -= ny;

				com[0] += del[0]*La/nx;
				com[1] += del[1]*Lb/ny;
	
				ncom[0] += 1;
				ncom[1] += 1;
			}
		}
	}

	double thresh_area = nlist * (La/nx) * (Lb/ny);
	
	double r_thresh = sqrt(thresh_area/M_PI);

	printf("xcom: %lf ycom: %lf zcom: %lf rthresh %lf thickness: %lf\n", com[0]/ncom[0], com[1]/ncom[1], com[2] / ncom[2], r_thresh, thickness);

	return -1;
}
