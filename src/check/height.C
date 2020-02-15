#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include <sys/time.h>
#include "interp.h"
#include "config.h"
#include "gsl/gsl_sf.h"
#include "gauss.h"
#include "mutil.h"

void surface::getHeight( double *h, double *r, int nx, int ny )
{
	double alpha_x = r[3*nv+0];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];
	if( !theFormulas )
		generatePlan();
	double *w = (double *)malloc( sizeof(double) * nx * ny );
	memset( w, 0, sizeof(double) * nx * ny );
	memset( h, 0, sizeof(double) * nx * ny );
 
	for( int f = 0; f < nf_faces; f++ )
	{
		for( int p = 0; p < nf_g_q_p; p++ )
		{
			int frm = f*nf_g_q_p+p;
			double R[3] = {0,0,0};
			double Ru[3] = { 0,0,0 };
			double Rv[3] = {0,0,0};
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
			}
		
			double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
			double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
			double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

			double g = sqrt(RuRu*RvRv-RuRv*RuRv);

			while( R[0] < 0 ) R[0] += PBC_vec[0][0]*alpha_x;
			while( R[0] >= PBC_vec[0][0]*alpha_x ) R[0] -= PBC_vec[0][0]*alpha_x;
			while( R[1] < 0 ) R[1] += PBC_vec[1][1]*alpha_y;
			while( R[1] >= PBC_vec[1][1]*alpha_y ) R[1] -= PBC_vec[1][1]*alpha_y;

			double fx = R[0] / (PBC_vec[0][0]*alpha_x);
			double fy = R[1] / (PBC_vec[1][1]*alpha_y);

			int ix = fx * nx;
			int iy = fy * ny;
			if( ix >= nx ) ix -= nx;
			if( iy >= ny ) iy -= ny;

			h[ix*ny+iy] += R[2] * theFormulas[frm].weight * theFormulas[frm].g0; 
			w[ix*ny+iy] += theFormulas[frm].weight * theFormulas[frm].g0;
		}		
	}

	for( int ix = 0; ix < nx; ix++ )
	for( int iy = 0; iy < ny; iy++ )
		h[ix*ny+iy] /= (1e-10 + w[ix*ny+iy]);
	free(w);
} 

void surface::getHeight2( double *h, double *r, int nx, int ny )
{
	double alpha_x = r[3*nv+0];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];
	if( !theVolumeFormulas )
		generateVolumePlan();
	double *w = (double *)malloc( sizeof(double) * nx * ny );
	memset( w, 0, sizeof(double) * nx * ny );
	memset( h, 0, sizeof(double) * nx * ny );

	double check_n = 0, check_d = 0;
	
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

				double R[3] = { (r1[0]+r2[0]+r3[0])/3,
						 (r1[1]+r2[1]+r3[1])/3,
						 (r1[2]+r2[2]+r3[2])/3 };

				while( R[0] < 0 ) R[0] += PBC_vec[0][0]*alpha_x;
				while( R[0] >= PBC_vec[0][0]*alpha_x ) R[0] -= PBC_vec[0][0]*alpha_x;
				while( R[1] < 0 ) R[1] += PBC_vec[1][1]*alpha_y;
				while( R[1] >= PBC_vec[1][1]*alpha_y ) R[1] -= PBC_vec[1][1]*alpha_y;
	
				double fx = R[0] / (PBC_vec[0][0]*alpha_x);
				double fy = R[1] / (PBC_vec[1][1]*alpha_y);
	
				int ix = fx * nx;
				int iy = fy * ny;
				if( ix >= nx ) ix -= nx;
				if( iy >= ny ) iy -= ny;


				h[ix*ny+iy] += R[2] * A; 
				w[ix*ny+iy] += A;
			}
		}	 
		
	} 
	
	for( int ix = 0; ix < nx; ix++ )
	for( int iy = 0; iy < ny; iy++ )
	{
		if( w[ix*ny+iy] < 1e-10 )
		{
			printf("BAD GRAIN.\n");
			exit(1);
		}
		h[ix*ny+iy] /= (1e-10 + w[ix*ny+iy]);
	}


	free(w);
} 

void surface::getHeight3( double *h, double *r, int nx, int ny )
{
	double alpha_x = r[3*nv+0];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];

	if( !theVolumeFormulas )
		generateVolumePlan();

	double check_n = 0, check_d = 0;

	for( int iy = 0; iy < ny; iy++ )
	for( int ix = 0; ix < nx; ix++ )
		h[iy*nx+ix] = 1e10;
	
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
			
				double xmm[2] = { 1e10,-1e10};	
				double ymm[2] = { 1e10,-1e10};	
		
				if( r1[0] < xmm[0] ) xmm[0] = r1[0];
				if( r1[0] > xmm[1] ) xmm[1] = r1[0];
				if( r1[1] < ymm[0] ) ymm[0] = r1[1];
				if( r1[1] > ymm[1] ) ymm[1] = r1[1];
				
				if( r2[0] < xmm[0] ) xmm[0] = r2[0];
				if( r2[0] > xmm[1] ) xmm[1] = r2[0];
				if( r2[1] < ymm[0] ) ymm[0] = r2[1];
				if( r2[1] > ymm[1] ) ymm[1] = r2[1];
				
				if( r3[0] < xmm[0] ) xmm[0] = r3[0];
				if( r3[0] > xmm[1] ) xmm[1] = r3[0];
				if( r3[1] < ymm[0] ) ymm[0] = r3[1];
				if( r3[1] > ymm[1] ) ymm[1] = r3[1];
	
				for( int dx = -1; dx <= 1; dx++ )
				for( int dy = -1; dy <= 1; dy++ )
				{
					double del[2] = { dx * PBC_vec[0][0], dy * PBC_vec[1][1] };

					int ileft  = nx*(xmm[0]+del[0])/PBC_vec[0][0];
					int iright = nx*(xmm[1]+del[0])/PBC_vec[0][0];
					
					int jleft  = ny*(ymm[0]+del[1])/PBC_vec[1][1];
					int jright = ny*(ymm[1]+del[1])/PBC_vec[1][1];

					if( xmm[0]+del[0] > 0 )
						ileft++;
					
					if( ymm[0]+del[1] > 0 )
						jleft++;

					if( ileft < 0 ) ileft = 0;
					if( iright >= nx ) iright = nx-1;
					
					if( jleft < 0 ) jleft = 0;
					if( jright >= ny ) jright = ny-1;

					for( int i = ileft; i <= iright; i++ )
					for( int j = jleft; j <= jright; j++ )
					{
						double p[3] = { i * PBC_vec[0][0] / nx-del[0], j * PBC_vec[1][1] / ny-del[1], 0};
				
							
						if( pointInTriangle( p, r1, r2, r3, 1e-8 ) )
						{
							double nrm[3];

							double dr1[3] = { r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]};
							double dr2[3] = { r3[0]-r1[0],r3[1]-r1[1],r3[2]-r1[2]};
		
							cross( dr1, dr2, nrm );
							normalize(nrm);
							double k = r1[0] * nrm[0] + r1[1] * nrm[1] + r1[2] * nrm[2];

							double pdp = p[0] * nrm[0] + p[1] * nrm[1];

							double t = (k - pdp) / nrm[2];
							// where in the triangle?
							h[j*nx+i] = t;
						}
					}	
			
				}	
			}
		}	 
		
	}

	for( int i = 0; i < ny; i++ )
	for( int j = 0; j < nx; j++ )
	{
		if( h[i*nx+j] > 1e9 )
		{		
			printf("missed point %d %d \n", i, j );
		}
	} 
} 

double surface::directFT( double *hq, double *hq2, double *r, int nx, int ny )
{
	double alpha_x = r[3*nv+0];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];

	if( !theVolumeFormulas )
		generateVolumePlan();

	double check_n = 0, check_d = 0;

	

	for( int iy = 0; iy < ny; iy++ )
	for( int ix = 0; ix < nx; ix++ )
	{
		hq[2*(iy*nx+ix)] = 0;
		hq[2*(iy*nx+ix)+1] = 0;
	}
	
	double Lx = PBC_vec[0][0]*alpha_x;
	double Ly = PBC_vec[1][1]*alpha_y;
	double Area = 0;
	double n = 0;

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
				R[0] += theVolumeFormulas[f].r_w[ioff*np+p] * (r[cp[p]*3+0] + theVolumeFormulas[f].r_pbc[ioff*np*3+3*p+0])*alpha_x; 
				R[1] += theVolumeFormulas[f].r_w[ioff*np+p] * (r[cp[p]*3+1] + theVolumeFormulas[f].r_pbc[ioff*np*3+3*p+1])*alpha_y; 
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
				
				double R[3] = { (r1[0]+r2[0]+r3[0])/3,
						 (r1[1]+r2[1]+r3[1])/3,
						 (r1[2]+r2[2]+r3[2])/3 };
				for( int qi = 0; qi < nx; qi++ )
				for( int qj = 0; qj < ny; qj++ )
				{
					double qx = 2 * M_PI * qi / Lx;
					double qy = 2 * M_PI * qj / Ly;
		
					if( qi > nx/2 )
						qx = -2 * M_PI * (nx-qi) / Lx; 
					if( qj > ny/2 )
						qy = -2 * M_PI * (ny-qj) / Ly; 

					hq[(qj*nx+qi)*2+0] += A * R[2] * cos( qx * R[0] + qy * R[1]  ); 
					hq[(qj*nx+qi)*2+1] += A * R[2] * sin( - qx * R[0] - qy * R[1]  ); 
				}

				Area += A;
			}
		}
	}

	
	int nq = 0;
	for( int qi = 0; qi < nx; qi++ )
	for( int qj = 0; qj < ny; qj++ )
	{
		double qx = 2 * M_PI * qi / Lx;
		double qy = 2 * M_PI * qj / Ly;

		if( qi > nx/2 )
			qx = -2 * M_PI * (nx-qi) / Lx; 
		if( qj > ny/2 )
			qy = -2 * M_PI * (ny-qj) / Ly; 

		hq2[2*nq+0] = sqrt(qx*qx+qy*qy);
		// could use the real area here.
		hq2[2*nq+1] = (1.0/(Lx*Ly)) * (hq[(qj*nx+qi)*2+0] * hq[(qj*nx+qi)*2+0] + hq[(qj*nx+qi)*2+1] * hq[(qj*nx+qi)*2+1]);
		nq++;
	}
	return hq[(0*nx+1)*2+0];
}
