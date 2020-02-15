#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "interp.h"

void surface::setupBoxing( double *r, double *lcoords, int npts, double rad, int do_pr_boxing ) // setup surface or pt boxing.
{
	if( !theFormulas ) generatePlan();

	boxing *theBoxing = ptBoxes;
	if( do_pr_boxing ) theBoxing = prBoxes;

	if( !theBoxing  )
	{
		if( do_pr_boxing )
			theBoxing = prBoxes = (boxing *)malloc( sizeof(boxing) );
		else 
			theBoxing = ptBoxes = (boxing *)malloc( sizeof(boxing) );
	
		double l1 = sqrt(PBC_vec[0][0] * PBC_vec[0][0] + PBC_vec[0][1]  * PBC_vec[0][1] + PBC_vec[0][2]*PBC_vec[0][2]);
		double l2 = sqrt(PBC_vec[1][0] * PBC_vec[1][0] + PBC_vec[1][1]  * PBC_vec[1][1] + PBC_vec[1][2]*PBC_vec[1][2]);
		double l3 = sqrt(PBC_vec[2][0] * PBC_vec[2][0] + PBC_vec[2][1]  * PBC_vec[2][1] + PBC_vec[2][2]*PBC_vec[2][2]);

		int nx = 20; // for surface boxing
		
		if( do_pr_boxing )
		{	// for pt_boxing
			nx = 1 + PBC_vec[0][0] / rad;	

			if( nx > 50 ) nx = 100;
		}

		int ny = nx*ceil(l2/l1);
		int nz = nx*ceil(l3/l1);

		printf("nx: %d ny: %d nz: %d\n", nx, ny, nz );

		double boxQ = nx*ny*nz*sizeof(boxel) / (1024*1024);

		if( boxQ > 100 )
		{
			printf("Box size greater than 100 megabytes.\n");
			exit(1);
		}
	
	
		if( nx % 2 == 0 ) nx += 1;
		if( ny % 2 == 0 ) ny += 1;
		if( nz % 2 == 0 ) nz += 1;
		theBoxing->nx = nx;
		theBoxing->ny = ny;
		theBoxing->nz = nz;
		theBoxing->boxes = (boxel *)malloc( sizeof(boxel) * nx * ny * nz );

		for( int ix = 0; ix < nx; ix++ )
		for( int iy = 0; iy < ny; iy++ )
		for( int iz = 0; iz < nz; iz++ )
		{
			int b = ix*ny*nz+iy*nz+iz;

			theBoxing->boxes[b].np = 0;
			theBoxing->boxes[b].npSpace = 0;
			theBoxing->boxes[b].f = NULL;
			theBoxing->boxes[b].p = NULL;
		}	
	}

	int nx = theBoxing->nx;
	int ny = theBoxing->ny;
	int nz = theBoxing->nz;
	double La = PBC_vec[0][0];
	double Lb = PBC_vec[1][1];
	double Lc = PBC_vec[2][2];

	if( do_pr_boxing == 0 )
	{
		for( int f = 0; f < nf_faces; f++ )
		{
			for( int p = 0; p < nf_g_q_p; p++ )
			{
				int frm = f*nf_g_q_p+p;
				double R[3] = {0,0,0};
				double nrm[3]={0,0,0}; 
	
				int *cp = theFormulas[f*nf_g_q_p+p].cp;
				int np = theFormulas[f*nf_g_q_p+p].ncoor;
	
				for( int xp = 0; xp < np; xp++ )
				{
					R[0] += theFormulas[frm].r_w[xp] * (r[cp[xp]*3+0] + theFormulas[frm].r_pbc[3*xp+0]); 
					R[1] += theFormulas[frm].r_w[xp] * (r[cp[xp]*3+1] + theFormulas[frm].r_pbc[3*xp+1]); 
					R[2] += theFormulas[frm].r_w[xp] * (r[cp[xp]*3+2] + theFormulas[frm].r_pbc[3*xp+2]); 
				}
	
				while( R[0] < 0 )             R[0] += La;
				while( R[0] > La ) R[0] -= La;
				while( R[1] < 0 )             R[1] += Lb;
				while( R[1] > Lb ) R[1] -= Lb;
				while( R[2] < 0 )             R[2] += Lc;
				while( R[2] > Lc ) R[2] -= Lc;
	
				int bx = nx*R[0] / La;
				int by = ny*R[1] / Lb;
				int bz = nz*R[2] / Lc;
		
				lcoords[(f*nf_g_q_p+p)*3+0] = R[0];
				lcoords[(f*nf_g_q_p+p)*3+1] = R[1];
				lcoords[(f*nf_g_q_p+p)*3+2] = R[2];
	
				int b = bx * ny*nz+by*nz+bz;
				if( theBoxing->boxes[b].npSpace == theBoxing->boxes[b].np )
				{
					if( theBoxing->boxes[b].npSpace == 0 )
					{
						theBoxing->boxes[b].npSpace = 10;
						theBoxing->boxes[b].f = (int *)malloc( sizeof(int) * theBoxing->boxes[b].npSpace );
						theBoxing->boxes[b].p = (int *)malloc( sizeof(int) * theBoxing->boxes[b].npSpace );
					}
					else
					{
						theBoxing->boxes[b].npSpace *= 2;
						theBoxing->boxes[b].f = (int *)realloc( theBoxing->boxes[b].f, sizeof(int) * theBoxing->boxes[b].npSpace );
						theBoxing->boxes[b].p = (int *)realloc( theBoxing->boxes[b].p, sizeof(int) * theBoxing->boxes[b].npSpace );
					}
	
				} 
					
				int npb =theBoxing->boxes[b].np;
	
				theBoxing->boxes[b].f[npb] = f;
				theBoxing->boxes[b].p[npb] = p;
				theBoxing->boxes[b].np += 1;
			}
		}
		
		for( int f = 0; f < nf_irr_faces; f++ )
		{
			for( int fp = 0; fp < nf_irr_pts; fp++ )
			{
				int frm = f*nf_irr_pts+fp;
				double R[3] = {0,0,0};
				double nrm[3]={0,0,0}; 
	
				int *cp = theIrregularFormulas[frm].cp;
				int np  = theIrregularFormulas[frm].ncoor;
	
				for( int xp = 0; xp < np; xp++ )
				{
					R[0] += theIrregularFormulas[frm].r_w[xp] * (r[cp[xp]*3+0] + theIrregularFormulas[frm].r_pbc[3*xp+0]); 
					R[1] += theIrregularFormulas[frm].r_w[xp] * (r[cp[xp]*3+1] + theIrregularFormulas[frm].r_pbc[3*xp+1]); 
					R[2] += theIrregularFormulas[frm].r_w[xp] * (r[cp[xp]*3+2] + theIrregularFormulas[frm].r_pbc[3*xp+2]); 
				}
	
				while( R[0] < 0 )             R[0] += La;
				while( R[0] > La ) R[0] -= La;
				while( R[1] < 0 )             R[1] += Lb;
				while( R[1] > Lb ) R[1] -= Lb;
				while( R[2] < 0 )             R[2] += Lc;
				while( R[2] > Lc ) R[2] -= Lc;
				
				lcoords[nf_faces*nf_g_q_p*3+(f*nf_irr_pts+fp)*3+0] = R[0];
				lcoords[nf_faces*nf_g_q_p*3+(f*nf_irr_pts+fp)*3+1] = R[1];
				lcoords[nf_faces*nf_g_q_p*3+(f*nf_irr_pts+fp)*3+2] = R[2];
	
				int bx = nx*R[0] / La;
				int by = ny*R[1] / Lb;
				int bz = nz*R[2] / Lc;
		
				int b = bx * ny*nz+by*nz+bz;
				if( theBoxing->boxes[b].npSpace == theBoxing->boxes[b].np )
				{
					if( theBoxing->boxes[b].npSpace == 0 )
					{
						theBoxing->boxes[b].npSpace = 10;
						theBoxing->boxes[b].f = (int *)malloc( sizeof(int) * theBoxing->boxes[b].npSpace );
						theBoxing->boxes[b].p = (int *)malloc( sizeof(int) * theBoxing->boxes[b].npSpace );
					}
					else
					{
						theBoxing->boxes[b].npSpace *= 2;
						theBoxing->boxes[b].f = (int *)realloc( theBoxing->boxes[b].f, sizeof(int) * theBoxing->boxes[b].npSpace );
						theBoxing->boxes[b].p = (int *)realloc( theBoxing->boxes[b].p, sizeof(int) * theBoxing->boxes[b].npSpace );
					}
	
					int np =theBoxing->boxes[b].np;
	
					theBoxing->boxes[b].f[np] = nf_faces+f;
					theBoxing->boxes[b].p[np] = fp;
					theBoxing->boxes[b].np += 1;
				} 
			}
		}
		theBoxing->pbox = NULL;
	}
	else
	{
		theBoxing->pbox = (int *)malloc( sizeof(int) * npts ); 
	}
}

void surface::addParticle( double *p_in, int id, double alpha_x, double alpha_y, double alpha_z )
{
	double La = PBC_vec[0][0];
	double Lb = PBC_vec[1][1];
	double Lc = PBC_vec[2][2];
	int nx = prBoxes->nx;
	int ny = prBoxes->ny;
	int nz = prBoxes->nz;
	int bx = floor(nx*p_in[0] / alpha_x / La);
	int by = floor(ny*p_in[1] / alpha_y / Lb);
	int bz = floor(nz*p_in[2] / alpha_z / Lc);


	while( bx < 0 ) bx += nx;
	while( bx >= nx ) bx -= nx;
	while( by < 0 ) by += ny;
	while( by >= ny ) by -= ny;
	while( bz < 0 ) bz += nz;
	while( bz >= nz ) bz -= nz;

	int b = bx * ny * nz + by * nz + bz;

	int np = prBoxes->boxes[b].np ;

	if( np == prBoxes->boxes[b].npSpace )
	{	
		if( prBoxes->boxes[b].npSpace == 0 )
		{
			prBoxes->boxes[b].npSpace = 1;
			prBoxes->boxes[b].f = (int *)malloc( sizeof(int) * prBoxes->boxes[b].npSpace );
			
		}
		prBoxes->boxes[b].npSpace *= 2;
		prBoxes->boxes[b].f = (int *)realloc( prBoxes->boxes[b].f, sizeof(int) * prBoxes->boxes[b].npSpace );
	}

	prBoxes->pbox[id] = b;
	prBoxes->boxes[b].f[np] = id;
	prBoxes->boxes[b].np += 1;
}

void surface::updateParticle( double *p_in, int id, double alpha_x, double alpha_y, double alpha_z )
{
	double La = PBC_vec[0][0]*alpha_x;
	double Lb = PBC_vec[1][1]*alpha_y;
	double Lc = PBC_vec[2][2]*alpha_z;
	int nx = prBoxes->nx;
	int ny = prBoxes->ny;
	int nz = prBoxes->nz;
	int bx = floor(nx*p_in[0] / La);
	int by = floor(ny*p_in[1] / Lb);
	int bz = floor(nz*p_in[2] / Lc);


	while( bx < 0 ) bx += nx;
	while( bx >= nx ) bx -= nx;
	while( by < 0 ) by += ny;
	while( by >= ny ) by -= ny;
	while( bz < 0 ) bz += nz;
	while( bz >= nz ) bz -= nz;

	int b = prBoxes->pbox[id];
	int np = prBoxes->boxes[b].np;
 
	int gotit = 0;
	for( int tx = 0; tx < prBoxes->boxes[b].np; tx++ )
	{
		if( prBoxes->boxes[b].f[tx] == id )
		{
			gotit = 1;	
			// remove the pt.
			prBoxes->boxes[b].f[tx] = prBoxes->boxes[b].f[prBoxes->boxes[b].np-1];
			prBoxes->boxes[b].np--;
			break;
		}
	}

	if( !gotit )
	{
		printf("Boxing particle retrieval error.\n");
		exit(1);
	}
	
	b = bx * ny * nz + by * nz + bz;
	np = prBoxes->boxes[b].np;
	
	if( np == prBoxes->boxes[b].npSpace )
	{
		if( prBoxes->boxes[b].npSpace == 0 )
		{
			prBoxes->boxes[b].npSpace = 1;
			prBoxes->boxes[b].f = (int *)malloc( sizeof(int) * prBoxes->boxes[b].npSpace );
			
		}
		prBoxes->boxes[b].npSpace *= 2;
		prBoxes->boxes[b].f = (int *)realloc( prBoxes->boxes[b].f, sizeof(int) * prBoxes->boxes[b].npSpace );
	}

	prBoxes->pbox[id] = b;
	prBoxes->boxes[b].f[np] = id;
	prBoxes->boxes[b].np += 1;
}

int surface::getNear( double *p_in, int id_in, double *all_p, double cut, double alpha_x, double alpha_y, double alpha_z, int *plist, double *rads )
{
	double La = PBC_vec[0][0]*alpha_x;
	double Lb = PBC_vec[1][1]*alpha_y;
	double Lc = PBC_vec[2][2]*alpha_z;

	int nx = prBoxes->nx;
	int ny = prBoxes->ny;
	int nz = prBoxes->nz;

	int bx = floor(nx*p_in[0] / La);
	int by = floor(ny*p_in[1] / Lb);
	int bz = floor(nz*p_in[2] / Lc);

	double binw = La / prBoxes->nx;

	int del = ceil(cut / binw)+1;
	
	if( 2*del >= nx ) del = nx/2; 
	if( 2*del >= ny ) del = ny/2; 
	if( 2*del >= nz ) del = nz/2; 

	double minr2 = 1e10;

	int np = 0;

	for( int dx = -del; dx <= del; dx++ )
	{
		int tx = bx + dx;
		while( tx < 0   ) tx += nx;
		while( tx >= nx ) tx -= nx;

		for( int dy = -del; dy <= del; dy++ )
		{
			int ty = by + dy;
			while( ty < 0   ) ty += ny;
			while( ty >= ny ) ty -= ny;

			for( int dz = -del; dz <= del; dz++ )
			{
				int tz = bz + dz;
				while( tz < 0   ) tz += nz;
				while( tz >= nz ) tz -= nz;

				int tb = tx*ny*nz+ty*nz+tz;
				for( int lp = 0; lp < prBoxes->boxes[tb].np; lp++ )
				{
					int f = prBoxes->boxes[tb].f[lp];
//					printf("f: %d t: %d %d %d lp: %d np: %d npSpace %d\n",
//						f, tx, ty, tz, lp, prBoxes->boxes[tb].np, prBoxes->boxes[tb].npSpace );
	
					if( id_in == f ) continue;

					double *br = all_p + f*3;

					double dr[3] = { p_in[0] - br[0], p_in[1] - br[1], p_in[2] - br[2] };
					if( dr[0] < -La/2 ) dr[0] += La; 
					if( dr[1] < -Lb/2 ) dr[1] += Lb; 
					if( dr[2] < -Lc/2 ) dr[2] += Lc; 
					if( dr[0] >  La/2 ) dr[0] -= La; 
					if( dr[1] >  Lb/2 ) dr[1] -= Lb; 
					if( dr[2] >  Lc/2 ) dr[2] -= Lc; 
				
					double r2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];

			//		if( r2 < (rads[f]+rads[id_in])*(rads[f]+rads[id_in])  )
					if( r2 < cut*cut )
					{
						plist[np] = f;
						np++;
					}
				}
			}
		}
	}
	return np;	
}

void surface::nearPt( double *p_in, double *lcoords, double cut, double alpha_x, double alpha_y, double alpha_z, int *bf, int *bp )
{
	if( !ptBoxes )
	{
		printf("nearPt called without boxing setup.\n");
		exit(1);
	}
	double La = PBC_vec[0][0];
	double Lb = PBC_vec[1][1];
	double Lc = PBC_vec[2][2];
	int nx = ptBoxes->nx;
	int ny = ptBoxes->ny;
	int nz = ptBoxes->nz;
	int bx = floor(nx*p_in[0] / alpha_x / La);
	int by = floor(ny*p_in[1] / alpha_y / Lb);
	int bz = floor(nz*p_in[2] / alpha_z / Lc);

	double binw = La / ptBoxes->nx;

	int del = 1+cut / binw;
	
	if( 2*del >= nx ) del = nx/2; 
	if( 2*del >= ny ) del = ny/2; 
	if( 2*del >= nz ) del = nz/2; 

	double minr2 = 1e10;

	*bf = -1;

	for( int dx = -del; dx <= del; dx++ )
	{
		int tx = bx + dx;
		while( tx < 0   ) tx += nx;
		while( tx >= nx ) tx -= nx;

		for( int dy = -del; dy <= del; dy++ )
		{
			int ty = by + dy;
			while( ty < 0   ) ty += ny;
			while( ty >= ny ) ty -= ny;

			for( int dz = -del; dz <= del; dz++ )
			{
				int tz = bz + dz;
				while( tz < 0   ) tz += nz;
				while( tz >= nz ) tz -= nz;

				int tb = tx*ny*nz+ty*nz+tz;
				for( int lp = 0; lp < ptBoxes->boxes[tb].np; lp++ )
				{
					int f = ptBoxes->boxes[tb].f[lp];
					int p = ptBoxes->boxes[tb].p[lp];
	
					double *br;
					if( f < nf_faces )
						br = lcoords+(f*nf_g_q_p+p)*3;
					else
						br = lcoords+(nf_faces*nf_g_q_p+(f-nf_faces)*nf_irr_pts+p)*3;

					double dr[3] = { p_in[0] - br[0], p_in[1] - br[1], p_in[2] - br[2] };
					if( dr[0] < -La/2 ) dr[0] += La; 
					if( dr[1] < -Lb/2 ) dr[1] += Lb; 
					if( dr[2] < -Lc/2 ) dr[2] += Lc; 
					if( dr[0] >  La/2 ) dr[0] -= La; 
					if( dr[1] >  Lb/2 ) dr[1] -= Lb; 
					if( dr[2] >  Lc/2 ) dr[2] -= Lc; 
				
					double r2 = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

					if( r2 < minr2 )
					{
						minr2 = r2;
						*bf = f;
						*bp = p;	
					}	
				}
			}
		}
	}
}




