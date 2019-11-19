#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include <sys/time.h>
#include "interp.h"
#include "3d_interp.h"
#include "fitRho.h"
#include "simulation.h"
#include "parallel.h"

// fitRho globals
double fitCoupling = 1.0;
double fitThickness = 15.0;

// fitRho statics
static double eps_f_min = 1.0;
int fitRho_activated = 0;
void Simulation::setupDensity( char *fileName )
{
	FILE *rhoFile = fopen(fileName, "r");
	char *buffer = (char *)malloc( sizeof(char) * 100000 );
	if( !rhoFile )
	{
		printf("Cannot open density file '%s'.\n", fileName );
		exit(1);
	}

	double La,Lb,Lc;
	int nx, ny, nz;
	
	fscanf(rhoFile, "%lf %lf %lf\n", &La, &Lb, &Lc );
	fscanf(rhoFile, "%d %d %d\n", &nx, &ny, &nz );

	double *rho = (double *)malloc( sizeof(double) * nx * ny * nz );

	for( int ix = 0; ix < nx; ix++ )
	for( int iy = 0; iy < ny; iy++ )
	{
		getLine( rhoFile, buffer );

		for( int iz = 0; iz < nz; iz++ )
			readNDoubles( buffer, rho+ix*ny*nz+iy*nz, nz );
	}

	fclose(rhoFile);
	
	double LXS = PBC_vec[0][0];
	double LYS = PBC_vec[1][1];
	double LZS = PBC_vec[2][2];

	double scale_x = La / LXS;
	double scale_y = Lb / LYS;
	double scale_z = Lc / LYS;

	PBC_vec[0][0] = La;
	PBC_vec[1][1] = Lb;
	PBC_vec[2][2] = Lc;

	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		for( int i = 0; i < sRec->theSurface->nv; i++ )
		{
			sRec->theSurface->theVertices[i].r[0] *= scale_x;
			sRec->theSurface->theVertices[i].r[1] *= scale_y;
			sRec->theSurface->theVertices[i].r[2] *= scale_z;
		}
	}
	double *rho_smooth = (double *)malloc( sizeof(double) * nx * ny * nz );
	
	double xbw = La/nx;
	double ybw = Lb/ny;
	double zbw = Lc/nz;

	double sigma =  sqrt(xbw*xbw+ybw*ybw+zbw*zbw)*3;

#define RUN_SMOOTHER

	#define N_SMOOTH	10
	printf("Smoothing.\n");
	for( int siter = 0; siter < N_SMOOTH; siter++ )
	{
		for( int ix = 0; ix < nx; ix++ )
		for( int iy = 0; iy < ny; iy++ )
		for( int iz = 0; iz < nz; iz++ )
		{
			double r = 0;

			double sum = 0;

			for( int dx = -1; dx <= 1; dx++ )
			for( int dy = -1; dy <= 1; dy++ )
			for( int dz = -1; dz <= 1; dz++ )
			{
				double dr[3] = { dx * xbw, dy * ybw, dz * zbw };
				double r2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
					
				int rx = ix+dx;
				int ry = iy+dy;
				int rz = iz+dz;

				if( rx < 0 ) rx += nx;
				if( rx >= nx ) rx -= nx;

				if( ry < 0 )   ry += ny;
				if( ry >= ny ) ry -= ny;
				
				if( rz < 0 )   rz += nz;
				if( rz >= nz ) rz -= nz;

				int alt_bin = rx*ny*nz+ry*nz+rz;

				sum += exp( -r2 / sigma );	
			}

			for( int dx = -1; dx <= 1; dx++ )
			for( int dy = -1; dy <= 1; dy++ )
			for( int dz = -1; dz <= 1; dz++ )
			{
				double dr[3] = { dx * xbw, dy * ybw, dz * zbw };
				double r2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
					
				int rx = ix+dx;
				int ry = iy+dy;
				int rz = iz+dz;

				if( rx < 0 ) rx += nx;
				if( rx >= nx ) rx -= nx;

				if( ry < 0 )   ry += ny;
				if( ry >= ny ) ry -= ny;
				
				if( rz < 0 )   rz += nz;
				if( rz >= nz ) rz -= nz;

				int alt_bin = rx*ny*nz+ry*nz+rz;

				r += rho[alt_bin] * exp( -r2 / sigma ) / sum;	
			}

			rho_smooth[ix*ny*nz+iy*nz+iz] = r;
		}
		memcpy( rho, rho_smooth, sizeof(double) * nx * ny *nz );
	}
	printf("Done smoothing.\n");
	setup_rho( rho, nx, ny, nz );
	free(buffer);	

	fitRho_activated = 1;
}

/* This returns the overlap energy and computes the overlap gradient of the surface with a 3d interpolated density of neutral surface atoms. */

double surface::rhoEnergy( double *r, double PBC_vec[3][3], double thickness_inner, double thickness_outer )
{
#ifdef PARALLEL
	if( par_info.my_id != BASE_TASK )
		return 0;
#endif
	if( ! fitRho_activated ) return 0;
	return rhoWorker( r, NULL, PBC_vec, 0, thickness_inner, thickness_outer, NULL, NULL );
}

double surface::rhoGrad( double *r, double *gr, double PBC_vec[3][3], double thickness_inner, double thickness_outer, double *tDerInner, double *tDerOuter )
{
#ifdef PARALLEL
	if( par_info.my_id != BASE_TASK )
		return 0;
#endif

	if( ! fitRho_activated ) return 0;

	return rhoWorker( r, gr, PBC_vec, 1, thickness_inner, thickness_outer, tDerInner, tDerOuter );
}

double surface::rhoWorker( double * r, double *gr, double PBC_vec[3][3], int do_grad, double thickness_inner, double thickness_outer, double *tDerInner, double *tDerOuter)
{
	if( !theFormulas )
		generatePlan();

	int static print_trigger = 0;

	double alpha_x = r[3*nv];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];

	double e = 0;

	int pgrid = 1;

	int use_max = 20;

	double *rGrad = (double *)malloc( sizeof(double) * 3 * use_max );
	double *nGrad = (double *)malloc( sizeof(double) * 3 * 3 * use_max );
	double *hGrad = (double *)malloc( sizeof(double) * 3 * use_max );
	double *kGrad = (double *)malloc( sizeof(double) * 3 * use_max );
	int *coor_list = (int *)malloc( sizeof(int) * use_max );
	int nCoor;

	double der_thickness[2] = { 0, 0 };

	int ngrid_pts = pgrid * (pgrid+1)/2;
	
	FILE *outFile = NULL;		
//	if( print_trigger % 10 == 0 && do_grad )
//		outFile = fopen("tempout.xyz","w");

	double ATOT = 0;
	for( int f = 0; f < nt; f++ )
	{
		if( f < nf_faces )
		{
			for( int q = 0; q < nf_g_q_p; q++ )
				ATOT += theFormulas[f*nf_g_q_p+q].g0*0.5;
		}
		else
		{
			for( int q = 0; q < nf_irr_pts; q++ )
				ATOT += theIrregularFormulas[(f-nf_faces)*nf_irr_pts+q].g0 * theIrregularFormulas[(f-nf_faces)*nf_irr_pts+q].weight;
		}
	}

	int nuv = pgrid*pgrid; 
	
	//
	double uv_array[2*nuv];
	double l = 1.0 / (double)pgrid;

	// advances per row:
	double d_row =  l * sin(60.0*M_PI/180.0);
	// advances per triangle;
	double d_row_half =  l * sin(30.0*M_PI/180.0);
	
	double uoff = d_row_half;
	int uv_pts = 0;
	for( int row = 0; row < pgrid; row++ )
	{
		double voff = d_row_half;

		for( int j = 0; j < pgrid-row; j++ )
		{
			uv_array[2*uv_pts+0] = row * l + (1.0/3.0) * l;
			uv_array[2*uv_pts+1] = j * l + (1.0/3.0) * l;
//			printf("%lf %lf\n", uv_array[2*uv_pts+0], uv_array[2*uv_pts+1]);
			uv_pts++;

			if( row > 0 )
			{
				uv_array[2*uv_pts+0] = row * l - (1.0/3.0) * l;
				uv_array[2*uv_pts+1] =   j * l + (2.0/3.0) * l;
//				printf("%lf %lf\n", uv_array[2*uv_pts+0], uv_array[2*uv_pts+1]);
				uv_pts++;
			}
		} 
	}


	int max_f;
	double max_u,max_v;
	double max_strain = 0;
	double max_c, max_k;

	for( int f = 0; f < nt; f++ )
	{
		double g0 = 0;
	
		if( f < nf_faces )
		{
			for( int q = 0; q < nf_g_q_p; q++ )
				g0 += theFormulas[f*nf_g_q_p+q].g0*0.5;

			double *rf = r + 3 * theFormulas[f*nf_g_q_p].cp[0];
		
			double dr[3] = { rf[0],rf[1],rf[2]-40 };

			double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

			if( r < 40.0 )
			{
//				printf("Break here.\n");
			}
		}
		else
		{
			for( int q = 0; q < nf_irr_pts; q++ )
				g0 += theIrregularFormulas[(f-nf_faces)*nf_irr_pts+q].g0 * theIrregularFormulas[(f-nf_faces)*nf_irr_pts+q].weight;
		}

		g0 /= ATOT;


		for( int tuv = 0; tuv < uv_pts; tuv++ )
		{
			double u = uv_array[2*tuv+0];
			double v = uv_array[2*tuv+1];

			double ctot=0;
			double k = 0;
			if( do_grad )
			{
				ctot = gradFetch( f, u, v, r,
				rGrad,
				nGrad,
				hGrad,
				kGrad,
				&nCoor,
				coor_list, &k );
			}
			else
				ctot = c( f, u, v, r, &k );
//#define FDIFF_CHECK		
#ifdef FDIFF_CHECK
{
			double eps = 1e-4;

			for( int cx = 0; cx < nCoor; cx++ )
			{
				for( int cart = 0; cart < 3; cart++ )
				{
					double dc[2];
					double dk[2];
					double dn[2][3];
					double dr[2][3];

					for( int pm = 0; pm < 2; pm++ )
					{
						r[3*coor_list[cx]+cart] += eps * (pm == 0 ? -1 : 1 );
						double k=0;
						dc[pm] = c(f,u,v,r,&k);
						dk[pm] = k;
						evaluateRNRM( f, u ,v, dr[pm], dn[pm], r ); 

						r[3*coor_list[cx]+cart] -= eps * (pm == 0 ? -1 : 1 );
					}
	
					double fd_k = (dk[1]-dk[0])/(2*eps);	
					double fd_c = (dc[1]-dc[0])/(2*eps);
					double fd_r[3] = { 
						(dr[1][0]-dr[0][0])/(2*eps),
						(dr[1][1]-dr[0][1])/(2*eps),
						(dr[1][2]-dr[0][2])/(2*eps) };
					double fd_n[3] = { 
						(dn[1][0]-dn[0][0])/(2*eps),
						(dn[1][1]-dn[0][1])/(2*eps),
						(dn[1][2]-dn[0][2])/(2*eps) };

					printf("f %d cart %d cx %d fdc: %.14le dc: %.14le fdk: %.14le dk: %.14le ac: %.14le fdr: %.14le ar: %.14le , fdn: %.14le %.14le %.14le an: %.14le %.14le %.14le\n",
						f,cart,cx,
						fd_c, 
						hGrad[3*cx+cart],
						fd_k,
						kGrad[3*cx+cart],
						fd_r[cart],
						rGrad[3*cx+cart],

						fd_n[0], fd_n[1], fd_n[2],
						nGrad[9*cx+cart+3*0], nGrad[9*cx+cart+3*1], nGrad[9*cx+cart+3*2] );
					
				}
			}
}			
#endif
			double rpt[3],nrm[3];

			evaluateRNRM( f, u, v, rpt, nrm, r );

			double use_thickness[2] = { thickness_inner, thickness_outer };
			
			double signs[2] = {-1,1};
			for( int leaflet = 0; leaflet < 2; leaflet++ )
			{
				double use_sign = signs[leaflet];

				double thickness = use_thickness[leaflet];
				double strain_f = exp( use_sign * thickness * ctot + thickness*thickness*k );

				double R[3] = { 
					rpt[0] + use_sign * nrm[0] * thickness* strain_f,
					rpt[1] + use_sign * nrm[1] * thickness* strain_f,
					rpt[2] + use_sign * nrm[2] * thickness* strain_f};

				double strain = strain_f - 1;

				if( strain_f > 200 )
				{
					e += 1e10;
					continue;
				}

				if( fabs(strain) > fabs(max_strain) )
				{
					max_f = f;
					max_u = u;
					max_v = v;
					max_c = ctot;
					max_k = k;
					max_strain = strain;
				}
				double d_R_d_thickness[3] = { 
						use_sign * nrm[0] * strain_f + use_sign * nrm[0] * thickness * strain_f * ( use_sign * ctot + 2 * thickness * k),
						use_sign * nrm[1] * strain_f + use_sign * nrm[1] * thickness * strain_f * ( use_sign * ctot + 2 * thickness * k),
						use_sign * nrm[2] * strain_f + use_sign * nrm[2] * thickness * strain_f * ( use_sign * ctot + 2 * thickness * k) };
//						use_sign * nrm[0] * (1+use_sign * thickness*ctot) + nrm[0] * ctot * thickness + use_sign * nrm[0] * 3 * thickness * thickness * k,
//						use_sign * nrm[1] * (1+use_sign * thickness*ctot) + nrm[1] * ctot * thickness + use_sign * nrm[1] * 3 * thickness * thickness * k,
//						use_sign * nrm[2] * (1+use_sign * thickness*ctot) + nrm[2] * ctot * thickness + use_sign * nrm[2] * 3 * thickness * thickness * k };

				double fractional[3] = { R[0] / PBC_vec[0][0], R[1] / PBC_vec[1][1], R[2] / PBC_vec[2][2] };

				e += -(g0/ngrid_pts) * fitCoupling * log( eps_f_min + eval_rho( fractional[0], fractional[1], fractional[2] ) );

				if( print_trigger % 100 == 0 && outFile )
				{
					fprintf(outFile, "%c %lf %lf %lf\n", (leaflet == 0 ? 'O' : 'C' ), R[0], R[1],R[2] );

				}
				if( do_grad )
				{
					double arg = eps_f_min + eval_rho( fractional[0], fractional[1], fractional[2] );

					double rho_g[3] = {0,0,0};
					eval_drho( fractional[0], fractional[1], fractional[2], rho_g );
	
					for( int cx = 0; cx < nCoor; cx++ )
					{
						double d_Rx_d_vx = rGrad[3*cx+0];
						double d_Rx_d_vy = 0;
						double d_Rx_d_vz = 0;
	
						double d_Ry_d_vx = 0;
						double d_Ry_d_vy = rGrad[3*cx+1];
						double d_Ry_d_vz = 0;
						
						double d_Rz_d_vx = 0;
						double d_Rz_d_vy = 0;
						double d_Rz_d_vz = rGrad[3*cx+2];
						d_Rx_d_vx += use_sign * nGrad[9*cx+0*3+0] * thickness * strain_f;
						d_Ry_d_vx += use_sign * nGrad[9*cx+1*3+0] * thickness * strain_f; 
						d_Rz_d_vx += use_sign * nGrad[9*cx+2*3+0] * thickness * strain_f; 
						
						d_Rx_d_vy += use_sign * nGrad[9*cx+0*3+1] * thickness * strain_f; 
						d_Ry_d_vy += use_sign * nGrad[9*cx+1*3+1] * thickness * strain_f; 
						d_Rz_d_vy += use_sign * nGrad[9*cx+2*3+1] * thickness * strain_f; 
						
						d_Rx_d_vz += use_sign * nGrad[9*cx+0*3+2] * thickness * strain_f; 
						d_Ry_d_vz += use_sign * nGrad[9*cx+1*3+2] * thickness * strain_f; 
						d_Rz_d_vz += use_sign * nGrad[9*cx+2*3+2] * thickness * strain_f; 
						
						d_Rx_d_vx += use_sign * nrm[0] * thickness * strain_f * thickness * use_sign  * hGrad[3*cx+0]; 
						d_Ry_d_vx += use_sign * nrm[1] * thickness * strain_f * thickness * use_sign  * hGrad[3*cx+0]; 
						d_Rz_d_vx += use_sign * nrm[2] * thickness * strain_f * thickness * use_sign  * hGrad[3*cx+0]; 
						                                                                             
						d_Rx_d_vy += use_sign * nrm[0] * thickness * strain_f * thickness * use_sign  * hGrad[3*cx+1]; 
						d_Ry_d_vy += use_sign * nrm[1] * thickness * strain_f * thickness * use_sign  * hGrad[3*cx+1]; 
						d_Rz_d_vy += use_sign * nrm[2] * thickness * strain_f * thickness * use_sign  * hGrad[3*cx+1]; 
						                                                                             
						d_Rx_d_vz += use_sign * nrm[0] * thickness * strain_f * thickness * use_sign  * hGrad[3*cx+2]; 
						d_Ry_d_vz += use_sign * nrm[1] * thickness * strain_f * thickness * use_sign  * hGrad[3*cx+2]; 
						d_Rz_d_vz += use_sign * nrm[2] * thickness * strain_f * thickness * use_sign  * hGrad[3*cx+2]; 
						
						d_Rx_d_vx += use_sign * nrm[0] * thickness * strain_f * thickness * thickness * kGrad[3*cx+0]; 
						d_Ry_d_vx += use_sign * nrm[1] * thickness * strain_f * thickness * thickness * kGrad[3*cx+0]; 
						d_Rz_d_vx += use_sign * nrm[2] * thickness * strain_f * thickness * thickness * kGrad[3*cx+0]; 
						                                                                                             
						d_Rx_d_vy += use_sign * nrm[0] * thickness * strain_f * thickness * thickness * kGrad[3*cx+1]; 
						d_Ry_d_vy += use_sign * nrm[1] * thickness * strain_f * thickness * thickness * kGrad[3*cx+1]; 
						d_Rz_d_vy += use_sign * nrm[2] * thickness * strain_f * thickness * thickness * kGrad[3*cx+1]; 
						                                                                                             
						d_Rx_d_vz += use_sign * nrm[0] * thickness * strain_f * thickness * thickness * kGrad[3*cx+2]; 
						d_Ry_d_vz += use_sign * nrm[1] * thickness * strain_f * thickness * thickness * kGrad[3*cx+2]; 
						d_Rz_d_vz += use_sign * nrm[2] * thickness * strain_f * thickness * thickness * kGrad[3*cx+2]; 

						gr[3*coor_list[cx]+0] -= fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[0] / PBC_vec[0][0] * alpha_x * d_Rx_d_vx;	
						gr[3*coor_list[cx]+1] -= fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[0] / PBC_vec[0][0] * alpha_x * d_Rx_d_vy;	
						gr[3*coor_list[cx]+2] -= fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[0] / PBC_vec[0][0] * alpha_x * d_Rx_d_vz;	
						
						gr[3*coor_list[cx]+0] -= fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[1] / PBC_vec[1][1] * alpha_y * d_Ry_d_vx;	
						gr[3*coor_list[cx]+1] -= fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[1] / PBC_vec[1][1] * alpha_y * d_Ry_d_vy;	
						gr[3*coor_list[cx]+2] -= fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[1] / PBC_vec[1][1] * alpha_y * d_Ry_d_vz;	
						
						gr[3*coor_list[cx]+0] -= fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[2] / PBC_vec[2][2] * alpha_z * d_Rz_d_vx;	
						gr[3*coor_list[cx]+1] -= fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[2] / PBC_vec[2][2] * alpha_z * d_Rz_d_vy;	
						gr[3*coor_list[cx]+2] -= fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[2] / PBC_vec[2][2] * alpha_z * d_Rz_d_vz;	

					}
						
					der_thickness[leaflet] -=  fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[0] / PBC_vec[0][0] * alpha_x  * d_R_d_thickness[0];
					der_thickness[leaflet] -=  fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[1] / PBC_vec[1][1] * alpha_y  * d_R_d_thickness[1];
					der_thickness[leaflet] -=  fitCoupling*(g0/ngrid_pts) * ( 1.0/arg) * rho_g[2] / PBC_vec[2][2] * alpha_z  * d_R_d_thickness[2];
				}
			}	
		}	
	} 

	free(rGrad);
	free(nGrad);
	free(hGrad);
	// thickness penalty to prohibit locking onto a single leaflet.

#if 1	
	double massive_k = 1e5;
	double min_thresh = 10.0;
	if( thickness_inner < min_thresh )
	{
		double dh = (thickness_inner-min_thresh);
		e += massive_k * (dh*dh);	
		if( do_grad )
			*tDerInner += 2 * massive_k * dh;
	}
	
	if( thickness_outer < min_thresh )
	{
		double dh = (thickness_outer-min_thresh);
		e += massive_k * (dh*dh);	

		if( do_grad )
			*tDerOuter += 2 * massive_k * dh;
	}
#endif
	if( do_grad )
	{
		printf("max strain %lf %lf %lf f %d u %lf v %lf\n", max_strain, max_c, max_k, max_f, max_u, max_v );
		*tDerInner += der_thickness[0];
		*tDerOuter += der_thickness[1];
	
		print_trigger++;

	}

	if( outFile )
	{
//		fclose(outFile);
	}

	return e;	
}

