#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "interp.h"
#include "mutil.h"
#include "l-bfgs.h"
#include "parallel.h"
#include "globals.h"

//#define DEBUG_NRMX
//#define DEBUG_G
double Sqrt( double );
double Power( double, double );
extern double kc;
extern double k_reg;
extern double KA;
extern double kg;
extern double k_strain;
extern int do_face, do_pt, debug_bit, on_surface;
static double THRESH = 1e-8;
double nrm_grad[6];

void surface::fdiff_check_grad( double *r )
{
	double *g = (double *)malloc( sizeof(double) * (3*nv+3) );

	memset( g, 0, sizeof(double) * (3*nv+3) );
	grad(r, g, NULL, NULL );

	double eps = 1e-5;
	
	const char *xyz = "xyz";

	for( int v = 0; v < nv; v++ )
	{
		for( int c = 0; c < 3; c++ )
		{
			r[3*v+c] += eps;

			double vp = energy( r, NULL );
			r[3*v+c] -= 2*eps;
			double vm = energy( r, NULL );


			double fd_g = (vp-vm) / (2*eps);

			double abs_diff = fabs(fd_g-g[3*v+c]);

			double perc_err = (2*fabs(fd_g-g[3*v+c]))/(fabs(fd_g)+fabs(g[3*v+c]));
	
			printf("vert %d cart %c fd %le g %le", v, xyz[c], fd_g, g[3*v+c] ); 
			if( perc_err < 1e-3 || abs_diff < 1e-7 )
				printf(" OK\n");
			else
				printf(" CHECK\n");
		}
	}

	free(g);

}

void surface::grad( double *r, double *gr, double *puv, double *pg )
{	
	if( !setup_for_parallel )
		setupParallel(this,NULL,0,-1);

	double alpha_x = r[3*nv];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];
	
	if( !theFormulas )
		generatePlan();

	double r_val = 0;
	double e = 0;
	double area = 0;
	double wgt = 0;
		double dudv = 0.5;

	for( int fx = 0; fx < par_info.nf; fx++ )
	{
		int f = par_info.faces[fx];

		if( f >= nf_faces )
			continue;

		double p_face_area = 0;

		int frm = f*nf_g_q_p;
		int t = theFormulas[frm].tri;
		//double g0 = 2*f_g0[f];
		for( int px = 0; px < theTriangles[t].np; px++ )
		{
			if( fabs(theTriangles[t].pc0[px] - theFormulas[frm].c0) < 1e-7 )
				continue; 
			p_face_area += theTriangles[t].pa[px];
		}
			if( debug_bit )
			{
			if( f != do_face )
				continue;
			}	
		double c0 = 0;
		double c1 = 0;

		double face_area = 0;
		double energy_density = 0;
	
		for( int p = 0; p < nf_g_q_p; p++ )
		{
			int frm = f*nf_g_q_p+p;
			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
			int *cp = theFormulas[f*nf_g_q_p+p].cp;
			int np = theFormulas[f*nf_g_q_p+p].ncoor;
			for( int p = 0; p < np; p++ )
			{
				R[0] += theFormulas[frm].r_w[p] * alpha_x * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				R[1] += theFormulas[frm].r_w[p] * alpha_y * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				R[2] += theFormulas[frm].r_w[p] * alpha_z * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				Ru[0] += theFormulas[frm].r_u[p] * alpha_x * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				Ru[1] += theFormulas[frm].r_u[p] * alpha_y * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				Ru[2] += theFormulas[frm].r_u[p] * alpha_z * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				Rv[0] += theFormulas[frm].r_v[p] * alpha_x * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				Rv[1] += theFormulas[frm].r_v[p] * alpha_y * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				Rv[2] += theFormulas[frm].r_v[p] * alpha_z * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuu[0] += theFormulas[frm].r_uu[p] * alpha_x * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuu[1] += theFormulas[frm].r_uu[p] * alpha_y * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuu[2] += theFormulas[frm].r_uu[p] * alpha_z * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuv[0] += theFormulas[frm].r_uv[p] * alpha_x * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuv[1] += theFormulas[frm].r_uv[p] * alpha_y * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuv[2] += theFormulas[frm].r_uv[p] * alpha_z * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSvv[0] += theFormulas[frm].r_vv[p] * alpha_x * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSvv[1] += theFormulas[frm].r_vv[p] * alpha_y * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSvv[2] += theFormulas[frm].r_vv[p] * alpha_z * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			}
		
			double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
			double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
			double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];
	
			double g = sqrt(RuRu*RvRv-RuRv*RuRv);
			double g0 = theFormulas[frm].g0;
			double nrm[3];		

			cross( Ru, Rv, nrm );
			normalize(nrm);


			double nsuu = tSuu[0] * nrm[0] + tSuu[1] * nrm[1] + tSuu[2] * nrm[2];
			double nsuv = tSuv[0] * nrm[0] + tSuv[1] * nrm[1] + tSuv[2] * nrm[2];
			double nsvv = tSvv[0] * nrm[0] + tSvv[1] * nrm[1] + tSvv[2] * nrm[2];
	
			double Stot = (nsuu * RvRv + nsvv * RuRu -2*nsuv*RuRv)/(g*g);

//		Int[i*n+j] = 0.5 * kc * Stot * Stot * g; 
//		AInt[i*n+j] = g;

			double Sop[4] = { 1.0/(g*g) * (nsuu * RvRv  - nsuv * RuRv), (1.0/(g*g)) * ( nsuv*RvRv-nsvv*RuRv),
					  1.0/(g*g) * (nsuv * RuRu  - nsuu * RuRv), (1.0/(g*g)) * ( nsvv*RuRu-nsuv*RuRv) };
	
			double a = Sop[0];		
			double b = Sop[1];
			double c = Sop[2];
			double d = Sop[3];
	
			double e1 = -0.5*(a+d-sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
			double e2 = -0.5*(a+d+sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
			
			double c0 = theFormulas[frm].c0;
//		printf("e1: %lf e2: %lf\n", e1, e2 );
			double en = 0.5 * kc * (e1 + e2 - c0 ) * (e1+e2 - c0);
			double dA = g * theFormulas[frm].weight;
			
	
			face_area  += g                   * theFormulas[frm].weight;	
			energy_density += 0.5 * kc * ( e1 + e2 - c0 ) * ( e1 + e2 - c0 ) * g * theFormulas[frm].weight;
		}

		for( int p = 0; p < nf_g_q_p; p++ )
		{
			double A = 0, A0 = 0;
			if( debug_bit )
			{
			if( f != do_face || p != do_pt )
				continue;
			}	
			int frm = f*nf_g_q_p+p;
			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
			double nrm[3]={0,0,0}; 
			double g0 = theFormulas[frm].g0;
			double RuRv0 = theFormulas[frm].RuRv0;

			int *cp = theFormulas[f*nf_g_q_p+p].cp;
			int np = theFormulas[f*nf_g_q_p+p].ncoor;
			double dedr[3*np];
			memset( dedr, 0, sizeof(double) * 3 *np );

			for( int p = 0; p < np; p++ )
			{
				R[0] += theFormulas[frm].r_w[p] * alpha_x * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				R[1] += theFormulas[frm].r_w[p] * alpha_y * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				R[2] += theFormulas[frm].r_w[p] * alpha_z * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				Ru[0] += theFormulas[frm].r_u[p] * alpha_x * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				Ru[1] += theFormulas[frm].r_u[p] * alpha_y * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				Ru[2] += theFormulas[frm].r_u[p] * alpha_z * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				Rv[0] += theFormulas[frm].r_v[p] * alpha_x * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				Rv[1] += theFormulas[frm].r_v[p] * alpha_y * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				Rv[2] += theFormulas[frm].r_v[p] * alpha_z * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuu[0] += theFormulas[frm].r_uu[p] * alpha_x * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuu[1] += theFormulas[frm].r_uu[p] * alpha_y * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuu[2] += theFormulas[frm].r_uu[p] * alpha_z * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuv[0] += theFormulas[frm].r_uv[p] * alpha_x * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuv[1] += theFormulas[frm].r_uv[p] * alpha_y * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuv[2] += theFormulas[frm].r_uv[p] * alpha_z * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSvv[0] += theFormulas[frm].r_vv[p] * alpha_x * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSvv[1] += theFormulas[frm].r_vv[p] * alpha_y * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSvv[2] += theFormulas[frm].r_vv[p] * alpha_z * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			}
		



		cross( Ru, Rv, nrm );
		normalize(nrm);

		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		double g = sqrt(RuRu*RvRv-RuRv*RuRv);

		double nsuu = tSuu[0] * nrm[0] + tSuu[1] * nrm[1] + tSuu[2] * nrm[2];
		double nsuv = tSuv[0] * nrm[0] + tSuv[1] * nrm[1] + tSuv[2] * nrm[2];
		double nsvv = tSvv[0] * nrm[0] + tSvv[1] * nrm[1] + tSvv[2] * nrm[2];

		double Stot = (nsuu * RvRv + nsvv * RuRu -2*nsuv*RuRv)/(g*g);

//		Int[i*n+j] = 0.5 * kc * Stot * Stot * g; 
//		AInt[i*n+j] = g;

		double Sop[4] = { 1.0/(g*g) * (nsuu * RvRv  - nsuv * RuRv), (1.0/(g*g)) * ( nsuv*RvRv-nsvv*RuRv),
				  1.0/(g*g) * (nsuv * RuRu  - nsuu * RuRv), (1.0/(g*g)) * ( nsvv*RuRu-nsuv*RuRv) };

		double a = Sop[0];		
		double b = Sop[1];
		double c = Sop[2];
		double d = Sop[3];

		double e1 = -0.5*(a+d-sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
		double e2 = -0.5*(a+d+sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
			
		double c0 = theFormulas[frm].c0;
//		printf("e1: %lf e2: %lf\n", e1, e2 );
		double en = 0.5 * kc * (e1 + e2 - c0 ) * (e1+e2 - c0);
		double dA = g * theFormulas[frm].weight;
		

		A += dudv * g * theFormulas[frm].weight;
		A0 += dudv * g0 * theFormulas[frm].weight;

		r_val += dudv * dA * (e1+e2)/2;
		area += dudv * dA;
		e += dudv * dA * en; 			

double d_nrmz_d_rux=0,d_nrmz_d_ruy=0,d_nrmz_d_ruz=0,d_nrmz_d_rvx=0,d_nrmz_d_rvy=0,d_nrmz_d_rvz=0,d_nrmy_d_rux=0,d_nrmy_d_ruy=0,d_nrmy_d_ruz=0,d_nrmy_d_rvx=0,d_nrmy_d_rvy=0,d_nrmy_d_rvz=0,d_nrmx_d_rux=0,d_nrmx_d_ruy=0,d_nrmx_d_ruz=0,d_nrmx_d_rvx=0,d_nrmx_d_rvy=0,d_nrmx_d_rvz=0,d_RuRu_d_rux=0,d_RuRu_d_ruy=0,d_RuRu_d_ruz=0,d_nsvv_d_rux=0,d_nsvv_d_ruy=0,d_nsvv_d_ruz=0,d_nsvv_d_rvx=0,d_nsvv_d_rvy=0,d_nsvv_d_rvz=0,d_nsvv_d_svvx=0,d_nsvv_d_svvy=0,d_nsvv_d_svvz=0,d_RuRv_d_rux=0,d_RuRv_d_ruy=0,d_RuRv_d_ruz=0,d_RuRv_d_rvx=0,d_RuRv_d_rvy=0,d_RuRv_d_rvz=0,d_nsuv_d_rux=0,d_nsuv_d_ruy=0,d_nsuv_d_ruz=0,d_nsuv_d_rvx=0,d_nsuv_d_rvy=0,d_nsuv_d_rvz=0,d_nsuv_d_suvx=0,d_nsuv_d_suvy=0,d_nsuv_d_suvz=0,d_RvRv_d_rvx=0,d_RvRv_d_rvy=0,d_RvRv_d_rvz=0,d_nsuu_d_rux=0,d_nsuu_d_ruy=0,d_nsuu_d_ruz=0,d_nsuu_d_rvx=0,d_nsuu_d_rvy=0,d_nsuu_d_rvz=0,d_nsuu_d_suux=0,d_nsuu_d_suuy=0,d_nsuu_d_suuz=0,d_g_d_rux=0,d_g_d_ruy=0,d_g_d_ruz=0,d_g_d_rvx=0,d_g_d_rvy=0,d_g_d_rvz=0,d_e_d_rux=0,d_e_d_ruy=0,d_e_d_ruz=0,d_e_d_rvx=0,d_e_d_rvy=0,d_e_d_rvz=0,d_e_d_suux=0,d_e_d_suuy=0,d_e_d_suuz=0,d_e_d_suvx=0,d_e_d_suvy=0,d_e_d_suvz=0,d_e_d_svvx=0,d_e_d_svvy=0,d_e_d_svvz=0,d_c2_d_rux=0,d_c2_d_ruy=0,d_c2_d_ruz=0,d_c2_d_rvx=0,d_c2_d_rvy=0,d_c2_d_rvz=0,d_c2_d_suux=0,d_c2_d_suuy=0,d_c2_d_suuz=0,d_c2_d_suvx=0,d_c2_d_suvy=0,d_c2_d_suvz=0,d_c2_d_svvx=0,d_c2_d_svvy=0,d_c2_d_svvz=0,d_d_d_rux=0,d_d_d_ruy=0,d_d_d_ruz=0,d_d_d_rvx=0,d_d_d_rvy=0,d_d_d_rvz=0,d_d_d_suvx=0,d_d_d_suvy=0,d_d_d_suvz=0,d_d_d_svvx=0,d_d_d_svvy=0,d_d_d_svvz=0,d_c_d_rux=0,d_c_d_ruy=0,d_c_d_ruz=0,d_c_d_rvx=0,d_c_d_rvy=0,d_c_d_rvz=0,d_c_d_suux=0,d_c_d_suuy=0,d_c_d_suuz=0,d_c_d_suvx=0,d_c_d_suvy=0,d_c_d_suvz=0,d_b_d_rux=0,d_b_d_ruy=0,d_b_d_ruz=0,d_b_d_rvx=0,d_b_d_rvy=0,d_b_d_rvz=0,d_b_d_suvx=0,d_b_d_suvy=0,d_b_d_suvz=0,d_b_d_svvx=0,d_b_d_svvy=0,d_b_d_svvz=0,d_a_d_rux=0,d_a_d_ruy=0,d_a_d_ruz=0,d_a_d_rvx=0,d_a_d_rvy=0,d_a_d_rvz=0,d_a_d_suux=0,d_a_d_suuy=0,d_a_d_suuz=0,d_a_d_suvx=0,d_a_d_suvy=0,d_a_d_suvz=0,d_c1_d_rux=0,d_c1_d_ruy=0,d_c1_d_ruz=0,d_c1_d_rvx=0,d_c1_d_rvy=0,d_c1_d_rvz=0,d_c1_d_suux=0,d_c1_d_suuy=0,d_c1_d_suuz=0,d_c1_d_suvx=0,d_c1_d_suvy=0,d_c1_d_suvz=0,d_c1_d_svvx=0,d_c1_d_svvy=0,d_c1_d_svvz=0,d_g_d_RvRv=0,d_g_d_RuRv=0,d_g_d_RuRu=0,d_nsvv_d_nrmz=0,d_nsvv_d_nrmy=0,d_nsvv_d_nrmx=0,d_nsuv_d_nrmz=0,d_nsuv_d_nrmy=0,d_nsuv_d_nrmx=0,d_nsuu_d_nrmz=0,d_nsuu_d_nrmy=0,d_nsuu_d_nrmx=0,d_d_d_RuRv=0,d_d_d_nsuv=0,d_d_d_RuRu=0,d_d_d_nsvv=0,d_d_d_g=0,d_c_d_RuRv=0,d_c_d_nsuv=0,d_c_d_RuRu=0,d_c_d_nsuu=0,d_c_d_g=0,d_b_d_RuRv=0,d_b_d_nsuv=0,d_b_d_RvRv=0,d_b_d_nsvv=0,d_b_d_g=0,d_a_d_RuRv=0,d_a_d_nsuv=0,d_a_d_RvRv=0,d_a_d_nsuu=0,d_a_d_g=0,d_e_d_c2=0,d_e_d_c1=0,d_e_d_g=0,d_c2_d_d=0,d_c2_d_c=0,d_c2_d_b=0,d_c2_d_a=0,d_c1_d_d=0,d_c1_d_c=0,d_c1_d_b=0,d_c1_d_a=0,junk;
			// the basic variables.

			d_g_d_RuRu = RvRv/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
			d_g_d_RuRv = -(RuRv/Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
			d_g_d_RvRv = RuRu/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  

			d_RuRu_d_rux = 2*Ru[0];
			d_RuRu_d_ruy = 2*Ru[1];
			d_RuRu_d_ruz = 2*Ru[2];
			
			d_RuRv_d_rux = Rv[0];
			d_RuRv_d_ruy = Rv[1];
			d_RuRv_d_ruz = Rv[2];
			
			d_RuRv_d_rvx = Ru[0];
			d_RuRv_d_rvy = Ru[1];
			d_RuRv_d_rvz = Ru[2];
			
			d_RvRv_d_rvx = 2*Rv[0];
			d_RvRv_d_rvy = 2*Rv[1];
			d_RvRv_d_rvz = 2*Rv[2];

			// intermediates.
			d_a_d_g = (-2*(-(nsuv*RuRv) + nsuu*RvRv))/Power(g,3);
			d_a_d_nsuu = RvRv/Power(g,2);
			d_a_d_RvRv = nsuu/Power(g,2);
			d_a_d_nsuv = -(RuRv/Power(g,2));
			d_a_d_RuRv = -(nsuv/Power(g,2));

			d_b_d_g = (-2*(-(nsvv*RuRv) + nsuv*RvRv))/Power(g,3);
			d_b_d_nsvv = -(RuRv/Power(g,2));
			d_b_d_RvRv = nsuv/Power(g,2);
			d_b_d_nsuv = RvRv/Power(g,2);
			d_b_d_RuRv = -(nsvv/Power(g,2));

			d_c_d_g = (-2*(nsuv*RuRu - nsuu*RuRv))/Power(g,3);
			d_c_d_nsuu = -(RuRv/Power(g,2));
			d_c_d_RuRu = nsuv/Power(g,2);
			d_c_d_nsuv = RuRu/Power(g,2);
			d_c_d_RuRv = -(nsuu/Power(g,2));

			d_d_d_g = (-2*(nsvv*RuRu - nsuv*RuRv))/Power(g,3);
			d_d_d_nsvv = RuRu/Power(g,2);
			d_d_d_RuRu = nsvv/Power(g,2);
			d_d_d_nsuv = -(RuRv/Power(g,2));
			d_d_d_RuRv = -(nsuv/Power(g,2));

			d_c1_d_a =-0.5*(1 - (2*a - 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
			d_c1_d_b =-(-1.*c)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2)) ;
			d_c1_d_c =-(-1.*b)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2)) ;
			d_c1_d_d =-0.5*(1 - (-2*a + 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
		
			d_c2_d_a = -0.5*(1 + (2*a - 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
			d_c2_d_b = -(1.*c)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2));
			d_c2_d_c = -(1.*b)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2));
			d_c2_d_d = -0.5*(1 + (-2*a + 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));

		/* here are the energy terms */

			d_e_d_g  = (0.5 * kc * ( e1 + e2 - c0) * ( e1 + e2 - c0) + kg * e1 * e2 ) * dudv * theFormulas[frm].weight;

			d_e_d_c1 = kc * (e1 + e2 - c0) * dudv * g * theFormulas[frm].weight;
			d_e_d_c2 = kc * (e1 + e2 - c0) * dudv * g * theFormulas[frm].weight;

			d_e_d_c1 += kg * e2 * dudv * g * theFormulas[frm].weight;
			d_e_d_c2 += kg * e1 * dudv * g * theFormulas[frm].weight;

			d_e_d_g += 2* KA * ((A-A0)/A0) * theFormulas[frm].weight * dudv;
			// particles need gradient.
			
			// derivative of dA in the numerator:
			d_e_d_g  += -p_face_area *  0.5 * kc * (e1+e2-c0)*(e1+e2-c0) * theFormulas[frm].weight / face_area;
			// derivative of dA in the denominator:
			d_e_d_g  += +p_face_area * energy_density * theFormulas[frm].weight / face_area / face_area;
			d_e_d_c1 += -p_face_area *  dA *  kc * (e1+e2-c0) / face_area;
			d_e_d_c2 += -p_face_area *  dA *  kc * (e1+e2-c0) / face_area;

		/* end energy terms */

			double fac = Power(Power(Ru[1]*Rv[0] - Ru[0]*Rv[1],2) + Power(Ru[2]*Rv[0] - Ru[0]*Rv[2],2) + Power(Ru[2]*Rv[1] - Ru[1]*Rv[2],2),1.5);

			d_nrmx_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac; 
			d_nrmx_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2)))/fac;
			d_nrmx_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac;
			
			d_nrmy_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2))))/fac;
			d_nrmy_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac; 
			d_nrmy_d_ruz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2])) + Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac;
	
			d_nrmz_d_rux = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2])/fac;
			d_nrmz_d_ruy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;
			d_nrmz_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2))) + (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;

			d_nrmx_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac; 
			d_nrmx_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2])))/fac;
			d_nrmx_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac;
			
			d_nrmy_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
			d_nrmy_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2])))/fac;
			d_nrmy_d_rvz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[0]*Ru[1]*Rv[0] - Power(Ru[0],2)*Rv[1] + Ru[2]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2])))/fac;
	
			d_nrmz_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;
			d_nrmz_d_rvy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[0]*Ru[2]*Rv[0] - Power(Ru[0],2)*Rv[2] + Ru[1]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
			d_nrmz_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;
			d_nsuu_d_nrmx = tSuu[0];
			d_nsuu_d_nrmy = tSuu[1];
			d_nsuu_d_nrmz = tSuu[2];
			
			d_nsuv_d_nrmx = tSuv[0];
			d_nsuv_d_nrmy = tSuv[1];
			d_nsuv_d_nrmz = tSuv[2];

			d_nsvv_d_nrmx = tSvv[0];
			d_nsvv_d_nrmy = tSvv[1];
			d_nsvv_d_nrmz = tSvv[2];

			d_nsuu_d_suux = nrm[0];
			d_nsuu_d_suuy = nrm[1];
			d_nsuu_d_suuz = nrm[2];
			
			d_nsuv_d_suvx = nrm[0];
			d_nsuv_d_suvy = nrm[1];
			d_nsuv_d_suvz = nrm[2];

			d_nsvv_d_svvx = nrm[0];
			d_nsvv_d_svvy = nrm[1];
			d_nsvv_d_svvz = nrm[2];



	d_nsvv_d_rux += d_nsvv_d_nrmz * d_nrmz_d_rux;
	d_nsvv_d_ruy += d_nsvv_d_nrmz * d_nrmz_d_ruy;
	d_nsvv_d_ruz += d_nsvv_d_nrmz * d_nrmz_d_ruz;
	d_nsvv_d_rvx += d_nsvv_d_nrmz * d_nrmz_d_rvx;
	d_nsvv_d_rvy += d_nsvv_d_nrmz * d_nrmz_d_rvy;
	d_nsvv_d_rvz += d_nsvv_d_nrmz * d_nrmz_d_rvz;
	d_nsvv_d_rux += d_nsvv_d_nrmy * d_nrmy_d_rux;
	d_nsvv_d_ruy += d_nsvv_d_nrmy * d_nrmy_d_ruy;
	d_nsvv_d_ruz += d_nsvv_d_nrmy * d_nrmy_d_ruz;
	d_nsvv_d_rvx += d_nsvv_d_nrmy * d_nrmy_d_rvx;
	d_nsvv_d_rvy += d_nsvv_d_nrmy * d_nrmy_d_rvy;
	d_nsvv_d_rvz += d_nsvv_d_nrmy * d_nrmy_d_rvz;
	d_nsvv_d_rux += d_nsvv_d_nrmx * d_nrmx_d_rux;
	d_nsvv_d_ruy += d_nsvv_d_nrmx * d_nrmx_d_ruy;
	d_nsvv_d_ruz += d_nsvv_d_nrmx * d_nrmx_d_ruz;
	d_nsvv_d_rvx += d_nsvv_d_nrmx * d_nrmx_d_rvx;
	d_nsvv_d_rvy += d_nsvv_d_nrmx * d_nrmx_d_rvy;
	d_nsvv_d_rvz += d_nsvv_d_nrmx * d_nrmx_d_rvz;
	d_nsuv_d_rux += d_nsuv_d_nrmz * d_nrmz_d_rux;
	d_nsuv_d_ruy += d_nsuv_d_nrmz * d_nrmz_d_ruy;
	d_nsuv_d_ruz += d_nsuv_d_nrmz * d_nrmz_d_ruz;
	d_nsuv_d_rvx += d_nsuv_d_nrmz * d_nrmz_d_rvx;
	d_nsuv_d_rvy += d_nsuv_d_nrmz * d_nrmz_d_rvy;
	d_nsuv_d_rvz += d_nsuv_d_nrmz * d_nrmz_d_rvz;
	d_nsuv_d_rux += d_nsuv_d_nrmy * d_nrmy_d_rux;
	d_nsuv_d_ruy += d_nsuv_d_nrmy * d_nrmy_d_ruy;
	d_nsuv_d_ruz += d_nsuv_d_nrmy * d_nrmy_d_ruz;
	d_nsuv_d_rvx += d_nsuv_d_nrmy * d_nrmy_d_rvx;
	d_nsuv_d_rvy += d_nsuv_d_nrmy * d_nrmy_d_rvy;
	d_nsuv_d_rvz += d_nsuv_d_nrmy * d_nrmy_d_rvz;
	d_nsuv_d_rux += d_nsuv_d_nrmx * d_nrmx_d_rux;
	d_nsuv_d_ruy += d_nsuv_d_nrmx * d_nrmx_d_ruy;
	d_nsuv_d_ruz += d_nsuv_d_nrmx * d_nrmx_d_ruz;
	d_nsuv_d_rvx += d_nsuv_d_nrmx * d_nrmx_d_rvx;
	d_nsuv_d_rvy += d_nsuv_d_nrmx * d_nrmx_d_rvy;
	d_nsuv_d_rvz += d_nsuv_d_nrmx * d_nrmx_d_rvz;
	d_nsuu_d_rux += d_nsuu_d_nrmz * d_nrmz_d_rux;
	d_nsuu_d_ruy += d_nsuu_d_nrmz * d_nrmz_d_ruy;
	d_nsuu_d_ruz += d_nsuu_d_nrmz * d_nrmz_d_ruz;
	d_nsuu_d_rvx += d_nsuu_d_nrmz * d_nrmz_d_rvx;
	d_nsuu_d_rvy += d_nsuu_d_nrmz * d_nrmz_d_rvy;
	d_nsuu_d_rvz += d_nsuu_d_nrmz * d_nrmz_d_rvz;
	d_nsuu_d_rux += d_nsuu_d_nrmy * d_nrmy_d_rux;
	d_nsuu_d_ruy += d_nsuu_d_nrmy * d_nrmy_d_ruy;
	d_nsuu_d_ruz += d_nsuu_d_nrmy * d_nrmy_d_ruz;
	d_nsuu_d_rvx += d_nsuu_d_nrmy * d_nrmy_d_rvx;
	d_nsuu_d_rvy += d_nsuu_d_nrmy * d_nrmy_d_rvy;
	d_nsuu_d_rvz += d_nsuu_d_nrmy * d_nrmy_d_rvz;
	d_nsuu_d_rux += d_nsuu_d_nrmx * d_nrmx_d_rux;
	d_nsuu_d_ruy += d_nsuu_d_nrmx * d_nrmx_d_ruy;
	d_nsuu_d_ruz += d_nsuu_d_nrmx * d_nrmx_d_ruz;
	d_nsuu_d_rvx += d_nsuu_d_nrmx * d_nrmx_d_rvx;
	d_nsuu_d_rvy += d_nsuu_d_nrmx * d_nrmx_d_rvy;
	d_nsuu_d_rvz += d_nsuu_d_nrmx * d_nrmx_d_rvz;
	d_g_d_rvx += d_g_d_RvRv * d_RvRv_d_rvx;
	d_g_d_rvy += d_g_d_RvRv * d_RvRv_d_rvy;
	d_g_d_rvz += d_g_d_RvRv * d_RvRv_d_rvz;
	d_g_d_rux += d_g_d_RuRv * d_RuRv_d_rux;
	d_g_d_ruy += d_g_d_RuRv * d_RuRv_d_ruy;
	d_g_d_ruz += d_g_d_RuRv * d_RuRv_d_ruz;
	d_g_d_rvx += d_g_d_RuRv * d_RuRv_d_rvx;
	d_g_d_rvy += d_g_d_RuRv * d_RuRv_d_rvy;
	d_g_d_rvz += d_g_d_RuRv * d_RuRv_d_rvz;
	d_g_d_rux += d_g_d_RuRu * d_RuRu_d_rux;
	d_g_d_ruy += d_g_d_RuRu * d_RuRu_d_ruy;
	d_g_d_ruz += d_g_d_RuRu * d_RuRu_d_ruz;
	d_d_d_rux += d_d_d_RuRv * d_RuRv_d_rux;
	d_d_d_ruy += d_d_d_RuRv * d_RuRv_d_ruy;
	d_d_d_ruz += d_d_d_RuRv * d_RuRv_d_ruz;
	d_d_d_rvx += d_d_d_RuRv * d_RuRv_d_rvx;
	d_d_d_rvy += d_d_d_RuRv * d_RuRv_d_rvy;
	d_d_d_rvz += d_d_d_RuRv * d_RuRv_d_rvz;
	d_d_d_rux += d_d_d_nsuv * d_nsuv_d_rux;
	d_d_d_ruy += d_d_d_nsuv * d_nsuv_d_ruy;
	d_d_d_ruz += d_d_d_nsuv * d_nsuv_d_ruz;
	d_d_d_rvx += d_d_d_nsuv * d_nsuv_d_rvx;
	d_d_d_rvy += d_d_d_nsuv * d_nsuv_d_rvy;
	d_d_d_rvz += d_d_d_nsuv * d_nsuv_d_rvz;
	d_d_d_suvx += d_d_d_nsuv * d_nsuv_d_suvx;
	d_d_d_suvy += d_d_d_nsuv * d_nsuv_d_suvy;
	d_d_d_suvz += d_d_d_nsuv * d_nsuv_d_suvz;
	d_d_d_rux += d_d_d_RuRu * d_RuRu_d_rux;
	d_d_d_ruy += d_d_d_RuRu * d_RuRu_d_ruy;
	d_d_d_ruz += d_d_d_RuRu * d_RuRu_d_ruz;
	d_d_d_rux += d_d_d_nsvv * d_nsvv_d_rux;
	d_d_d_ruy += d_d_d_nsvv * d_nsvv_d_ruy;
	d_d_d_ruz += d_d_d_nsvv * d_nsvv_d_ruz;
	d_d_d_rvx += d_d_d_nsvv * d_nsvv_d_rvx;
	d_d_d_rvy += d_d_d_nsvv * d_nsvv_d_rvy;
	d_d_d_rvz += d_d_d_nsvv * d_nsvv_d_rvz;
	d_d_d_svvx += d_d_d_nsvv * d_nsvv_d_svvx;
	d_d_d_svvy += d_d_d_nsvv * d_nsvv_d_svvy;
	d_d_d_svvz += d_d_d_nsvv * d_nsvv_d_svvz;
	d_d_d_rux += d_d_d_g * d_g_d_rux;
	d_d_d_ruy += d_d_d_g * d_g_d_ruy;
	d_d_d_ruz += d_d_d_g * d_g_d_ruz;
	d_d_d_rvx += d_d_d_g * d_g_d_rvx;
	d_d_d_rvy += d_d_d_g * d_g_d_rvy;
	d_d_d_rvz += d_d_d_g * d_g_d_rvz;
	d_c_d_rux += d_c_d_RuRv * d_RuRv_d_rux;
	d_c_d_ruy += d_c_d_RuRv * d_RuRv_d_ruy;
	d_c_d_ruz += d_c_d_RuRv * d_RuRv_d_ruz;
	d_c_d_rvx += d_c_d_RuRv * d_RuRv_d_rvx;
	d_c_d_rvy += d_c_d_RuRv * d_RuRv_d_rvy;
	d_c_d_rvz += d_c_d_RuRv * d_RuRv_d_rvz;
	d_c_d_rux += d_c_d_nsuv * d_nsuv_d_rux;
	d_c_d_ruy += d_c_d_nsuv * d_nsuv_d_ruy;
	d_c_d_ruz += d_c_d_nsuv * d_nsuv_d_ruz;
	d_c_d_rvx += d_c_d_nsuv * d_nsuv_d_rvx;
	d_c_d_rvy += d_c_d_nsuv * d_nsuv_d_rvy;
	d_c_d_rvz += d_c_d_nsuv * d_nsuv_d_rvz;
	d_c_d_suvx += d_c_d_nsuv * d_nsuv_d_suvx;
	d_c_d_suvy += d_c_d_nsuv * d_nsuv_d_suvy;
	d_c_d_suvz += d_c_d_nsuv * d_nsuv_d_suvz;
	d_c_d_rux += d_c_d_RuRu * d_RuRu_d_rux;
	d_c_d_ruy += d_c_d_RuRu * d_RuRu_d_ruy;
	d_c_d_ruz += d_c_d_RuRu * d_RuRu_d_ruz;
	d_c_d_rux += d_c_d_nsuu * d_nsuu_d_rux;
	d_c_d_ruy += d_c_d_nsuu * d_nsuu_d_ruy;
	d_c_d_ruz += d_c_d_nsuu * d_nsuu_d_ruz;
	d_c_d_rvx += d_c_d_nsuu * d_nsuu_d_rvx;
	d_c_d_rvy += d_c_d_nsuu * d_nsuu_d_rvy;
	d_c_d_rvz += d_c_d_nsuu * d_nsuu_d_rvz;
	d_c_d_suux += d_c_d_nsuu * d_nsuu_d_suux;
	d_c_d_suuy += d_c_d_nsuu * d_nsuu_d_suuy;
	d_c_d_suuz += d_c_d_nsuu * d_nsuu_d_suuz;
	d_c_d_rux += d_c_d_g * d_g_d_rux;
	d_c_d_ruy += d_c_d_g * d_g_d_ruy;
	d_c_d_ruz += d_c_d_g * d_g_d_ruz;
	d_c_d_rvx += d_c_d_g * d_g_d_rvx;
	d_c_d_rvy += d_c_d_g * d_g_d_rvy;
	d_c_d_rvz += d_c_d_g * d_g_d_rvz;
	d_b_d_rux += d_b_d_RuRv * d_RuRv_d_rux;
	d_b_d_ruy += d_b_d_RuRv * d_RuRv_d_ruy;
	d_b_d_ruz += d_b_d_RuRv * d_RuRv_d_ruz;
	d_b_d_rvx += d_b_d_RuRv * d_RuRv_d_rvx;
	d_b_d_rvy += d_b_d_RuRv * d_RuRv_d_rvy;
	d_b_d_rvz += d_b_d_RuRv * d_RuRv_d_rvz;
	d_b_d_rux += d_b_d_nsuv * d_nsuv_d_rux;
	d_b_d_ruy += d_b_d_nsuv * d_nsuv_d_ruy;
	d_b_d_ruz += d_b_d_nsuv * d_nsuv_d_ruz;
	d_b_d_rvx += d_b_d_nsuv * d_nsuv_d_rvx;
	d_b_d_rvy += d_b_d_nsuv * d_nsuv_d_rvy;
	d_b_d_rvz += d_b_d_nsuv * d_nsuv_d_rvz;
	d_b_d_suvx += d_b_d_nsuv * d_nsuv_d_suvx;
	d_b_d_suvy += d_b_d_nsuv * d_nsuv_d_suvy;
	d_b_d_suvz += d_b_d_nsuv * d_nsuv_d_suvz;
	d_b_d_rvx += d_b_d_RvRv * d_RvRv_d_rvx;
	d_b_d_rvy += d_b_d_RvRv * d_RvRv_d_rvy;
	d_b_d_rvz += d_b_d_RvRv * d_RvRv_d_rvz;
	d_b_d_rux += d_b_d_nsvv * d_nsvv_d_rux;
	d_b_d_ruy += d_b_d_nsvv * d_nsvv_d_ruy;
	d_b_d_ruz += d_b_d_nsvv * d_nsvv_d_ruz;
	d_b_d_rvx += d_b_d_nsvv * d_nsvv_d_rvx;
	d_b_d_rvy += d_b_d_nsvv * d_nsvv_d_rvy;
	d_b_d_rvz += d_b_d_nsvv * d_nsvv_d_rvz;
	d_b_d_svvx += d_b_d_nsvv * d_nsvv_d_svvx;
	d_b_d_svvy += d_b_d_nsvv * d_nsvv_d_svvy;
	d_b_d_svvz += d_b_d_nsvv * d_nsvv_d_svvz;
	d_b_d_rux += d_b_d_g * d_g_d_rux;
	d_b_d_ruy += d_b_d_g * d_g_d_ruy;
	d_b_d_ruz += d_b_d_g * d_g_d_ruz;
	d_b_d_rvx += d_b_d_g * d_g_d_rvx;
	d_b_d_rvy += d_b_d_g * d_g_d_rvy;
	d_b_d_rvz += d_b_d_g * d_g_d_rvz;
	d_a_d_rux += d_a_d_RuRv * d_RuRv_d_rux;
	d_a_d_ruy += d_a_d_RuRv * d_RuRv_d_ruy;
	d_a_d_ruz += d_a_d_RuRv * d_RuRv_d_ruz;
	d_a_d_rvx += d_a_d_RuRv * d_RuRv_d_rvx;
	d_a_d_rvy += d_a_d_RuRv * d_RuRv_d_rvy;
	d_a_d_rvz += d_a_d_RuRv * d_RuRv_d_rvz;
	d_a_d_rux += d_a_d_nsuv * d_nsuv_d_rux;
	d_a_d_ruy += d_a_d_nsuv * d_nsuv_d_ruy;
	d_a_d_ruz += d_a_d_nsuv * d_nsuv_d_ruz;
	d_a_d_rvx += d_a_d_nsuv * d_nsuv_d_rvx;
	d_a_d_rvy += d_a_d_nsuv * d_nsuv_d_rvy;
	d_a_d_rvz += d_a_d_nsuv * d_nsuv_d_rvz;
	d_a_d_suvx += d_a_d_nsuv * d_nsuv_d_suvx;
	d_a_d_suvy += d_a_d_nsuv * d_nsuv_d_suvy;
	d_a_d_suvz += d_a_d_nsuv * d_nsuv_d_suvz;
	d_a_d_rvx += d_a_d_RvRv * d_RvRv_d_rvx;
	d_a_d_rvy += d_a_d_RvRv * d_RvRv_d_rvy;
	d_a_d_rvz += d_a_d_RvRv * d_RvRv_d_rvz;
	d_a_d_rux += d_a_d_nsuu * d_nsuu_d_rux;
	d_a_d_ruy += d_a_d_nsuu * d_nsuu_d_ruy;
	d_a_d_ruz += d_a_d_nsuu * d_nsuu_d_ruz;
	d_a_d_rvx += d_a_d_nsuu * d_nsuu_d_rvx;
	d_a_d_rvy += d_a_d_nsuu * d_nsuu_d_rvy;
	d_a_d_rvz += d_a_d_nsuu * d_nsuu_d_rvz;
	d_a_d_suux += d_a_d_nsuu * d_nsuu_d_suux;
	d_a_d_suuy += d_a_d_nsuu * d_nsuu_d_suuy;
	d_a_d_suuz += d_a_d_nsuu * d_nsuu_d_suuz;
	d_a_d_rux += d_a_d_g * d_g_d_rux;
	d_a_d_ruy += d_a_d_g * d_g_d_ruy;
	d_a_d_ruz += d_a_d_g * d_g_d_ruz;
	d_a_d_rvx += d_a_d_g * d_g_d_rvx;
	d_a_d_rvy += d_a_d_g * d_g_d_rvy;
	d_a_d_rvz += d_a_d_g * d_g_d_rvz;
	d_c1_d_rux += d_c1_d_d * d_d_d_rux;
	d_c1_d_ruy += d_c1_d_d * d_d_d_ruy;
	d_c1_d_ruz += d_c1_d_d * d_d_d_ruz;
	d_c1_d_rvx += d_c1_d_d * d_d_d_rvx;
	d_c1_d_rvy += d_c1_d_d * d_d_d_rvy;
	d_c1_d_rvz += d_c1_d_d * d_d_d_rvz;
	d_c1_d_suvx += d_c1_d_d * d_d_d_suvx;
	d_c1_d_suvy += d_c1_d_d * d_d_d_suvy;
	d_c1_d_suvz += d_c1_d_d * d_d_d_suvz;
	d_c1_d_svvx += d_c1_d_d * d_d_d_svvx;
	d_c1_d_svvy += d_c1_d_d * d_d_d_svvy;
	d_c1_d_svvz += d_c1_d_d * d_d_d_svvz;
	d_c1_d_rux += d_c1_d_c * d_c_d_rux;
	d_c1_d_ruy += d_c1_d_c * d_c_d_ruy;
	d_c1_d_ruz += d_c1_d_c * d_c_d_ruz;
	d_c1_d_rvx += d_c1_d_c * d_c_d_rvx;
	d_c1_d_rvy += d_c1_d_c * d_c_d_rvy;
	d_c1_d_rvz += d_c1_d_c * d_c_d_rvz;
	d_c1_d_suux += d_c1_d_c * d_c_d_suux;
	d_c1_d_suuy += d_c1_d_c * d_c_d_suuy;
	d_c1_d_suuz += d_c1_d_c * d_c_d_suuz;
	d_c1_d_suvx += d_c1_d_c * d_c_d_suvx;
	d_c1_d_suvy += d_c1_d_c * d_c_d_suvy;
	d_c1_d_suvz += d_c1_d_c * d_c_d_suvz;
	d_c1_d_rux += d_c1_d_b * d_b_d_rux;
	d_c1_d_ruy += d_c1_d_b * d_b_d_ruy;
	d_c1_d_ruz += d_c1_d_b * d_b_d_ruz;
	d_c1_d_rvx += d_c1_d_b * d_b_d_rvx;
	d_c1_d_rvy += d_c1_d_b * d_b_d_rvy;
	d_c1_d_rvz += d_c1_d_b * d_b_d_rvz;
	d_c1_d_suvx += d_c1_d_b * d_b_d_suvx;
	d_c1_d_suvy += d_c1_d_b * d_b_d_suvy;
	d_c1_d_suvz += d_c1_d_b * d_b_d_suvz;
	d_c1_d_svvx += d_c1_d_b * d_b_d_svvx;
	d_c1_d_svvy += d_c1_d_b * d_b_d_svvy;
	d_c1_d_svvz += d_c1_d_b * d_b_d_svvz;
	d_c1_d_rux += d_c1_d_a * d_a_d_rux;
	d_c1_d_ruy += d_c1_d_a * d_a_d_ruy;
	d_c1_d_ruz += d_c1_d_a * d_a_d_ruz;
	d_c1_d_rvx += d_c1_d_a * d_a_d_rvx;
	d_c1_d_rvy += d_c1_d_a * d_a_d_rvy;
	d_c1_d_rvz += d_c1_d_a * d_a_d_rvz;
	d_c1_d_suux += d_c1_d_a * d_a_d_suux;
	d_c1_d_suuy += d_c1_d_a * d_a_d_suuy;
	d_c1_d_suuz += d_c1_d_a * d_a_d_suuz;
	d_c1_d_suvx += d_c1_d_a * d_a_d_suvx;
	d_c1_d_suvy += d_c1_d_a * d_a_d_suvy;
	d_c1_d_suvz += d_c1_d_a * d_a_d_suvz;
	d_c2_d_rux += d_c2_d_d * d_d_d_rux;
	d_c2_d_ruy += d_c2_d_d * d_d_d_ruy;
	d_c2_d_ruz += d_c2_d_d * d_d_d_ruz;
	d_c2_d_rvx += d_c2_d_d * d_d_d_rvx;
	d_c2_d_rvy += d_c2_d_d * d_d_d_rvy;
	d_c2_d_rvz += d_c2_d_d * d_d_d_rvz;
	d_c2_d_suvx += d_c2_d_d * d_d_d_suvx;
	d_c2_d_suvy += d_c2_d_d * d_d_d_suvy;
	d_c2_d_suvz += d_c2_d_d * d_d_d_suvz;
	d_c2_d_svvx += d_c2_d_d * d_d_d_svvx;
	d_c2_d_svvy += d_c2_d_d * d_d_d_svvy;
	d_c2_d_svvz += d_c2_d_d * d_d_d_svvz;
	d_c2_d_rux += d_c2_d_c * d_c_d_rux;
	d_c2_d_ruy += d_c2_d_c * d_c_d_ruy;
	d_c2_d_ruz += d_c2_d_c * d_c_d_ruz;
	d_c2_d_rvx += d_c2_d_c * d_c_d_rvx;
	d_c2_d_rvy += d_c2_d_c * d_c_d_rvy;
	d_c2_d_rvz += d_c2_d_c * d_c_d_rvz;
	d_c2_d_suux += d_c2_d_c * d_c_d_suux;
	d_c2_d_suuy += d_c2_d_c * d_c_d_suuy;
	d_c2_d_suuz += d_c2_d_c * d_c_d_suuz;
	d_c2_d_suvx += d_c2_d_c * d_c_d_suvx;
	d_c2_d_suvy += d_c2_d_c * d_c_d_suvy;
	d_c2_d_suvz += d_c2_d_c * d_c_d_suvz;
	d_c2_d_rux += d_c2_d_b * d_b_d_rux;
	d_c2_d_ruy += d_c2_d_b * d_b_d_ruy;
	d_c2_d_ruz += d_c2_d_b * d_b_d_ruz;
	d_c2_d_rvx += d_c2_d_b * d_b_d_rvx;
	d_c2_d_rvy += d_c2_d_b * d_b_d_rvy;
	d_c2_d_rvz += d_c2_d_b * d_b_d_rvz;
	d_c2_d_suvx += d_c2_d_b * d_b_d_suvx;
	d_c2_d_suvy += d_c2_d_b * d_b_d_suvy;
	d_c2_d_suvz += d_c2_d_b * d_b_d_suvz;
	d_c2_d_svvx += d_c2_d_b * d_b_d_svvx;
	d_c2_d_svvy += d_c2_d_b * d_b_d_svvy;
	d_c2_d_svvz += d_c2_d_b * d_b_d_svvz;
	d_c2_d_rux += d_c2_d_a * d_a_d_rux;
	d_c2_d_ruy += d_c2_d_a * d_a_d_ruy;
	d_c2_d_ruz += d_c2_d_a * d_a_d_ruz;
	d_c2_d_rvx += d_c2_d_a * d_a_d_rvx;
	d_c2_d_rvy += d_c2_d_a * d_a_d_rvy;
	d_c2_d_rvz += d_c2_d_a * d_a_d_rvz;
	d_c2_d_suux += d_c2_d_a * d_a_d_suux;
	d_c2_d_suuy += d_c2_d_a * d_a_d_suuy;
	d_c2_d_suuz += d_c2_d_a * d_a_d_suuz;
	d_c2_d_suvx += d_c2_d_a * d_a_d_suvx;
	d_c2_d_suvy += d_c2_d_a * d_a_d_suvy;
	d_c2_d_suvz += d_c2_d_a * d_a_d_suvz;
	d_e_d_rux += d_e_d_c2 * d_c2_d_rux;
	d_e_d_ruy += d_e_d_c2 * d_c2_d_ruy;
	d_e_d_ruz += d_e_d_c2 * d_c2_d_ruz;
	d_e_d_rvx += d_e_d_c2 * d_c2_d_rvx;
	d_e_d_rvy += d_e_d_c2 * d_c2_d_rvy;
	d_e_d_rvz += d_e_d_c2 * d_c2_d_rvz;
	d_e_d_suux += d_e_d_c2 * d_c2_d_suux;
	d_e_d_suuy += d_e_d_c2 * d_c2_d_suuy;
	d_e_d_suuz += d_e_d_c2 * d_c2_d_suuz;
	d_e_d_suvx += d_e_d_c2 * d_c2_d_suvx;
	d_e_d_suvy += d_e_d_c2 * d_c2_d_suvy;
	d_e_d_suvz += d_e_d_c2 * d_c2_d_suvz;
	d_e_d_svvx += d_e_d_c2 * d_c2_d_svvx;
	d_e_d_svvy += d_e_d_c2 * d_c2_d_svvy;
	d_e_d_svvz += d_e_d_c2 * d_c2_d_svvz;
	d_e_d_rux += d_e_d_c1 * d_c1_d_rux;
	d_e_d_ruy += d_e_d_c1 * d_c1_d_ruy;
	d_e_d_ruz += d_e_d_c1 * d_c1_d_ruz;
	d_e_d_rvx += d_e_d_c1 * d_c1_d_rvx;
	d_e_d_rvy += d_e_d_c1 * d_c1_d_rvy;
	d_e_d_rvz += d_e_d_c1 * d_c1_d_rvz;
	d_e_d_suux += d_e_d_c1 * d_c1_d_suux;
	d_e_d_suuy += d_e_d_c1 * d_c1_d_suuy;
	d_e_d_suuz += d_e_d_c1 * d_c1_d_suuz;
	d_e_d_suvx += d_e_d_c1 * d_c1_d_suvx;
	d_e_d_suvy += d_e_d_c1 * d_c1_d_suvy;
	d_e_d_suvz += d_e_d_c1 * d_c1_d_suvz;
	d_e_d_svvx += d_e_d_c1 * d_c1_d_svvx;
	d_e_d_svvy += d_e_d_c1 * d_c1_d_svvy;
	d_e_d_svvz += d_e_d_c1 * d_c1_d_svvz;
	d_e_d_rux += d_e_d_g * d_g_d_rux;
	d_e_d_ruy += d_e_d_g * d_g_d_ruy;
	d_e_d_ruz += d_e_d_g * d_g_d_ruz;
	d_e_d_rvx += d_e_d_g * d_g_d_rvx;
	d_e_d_rvy += d_e_d_g * d_g_d_rvy;
	d_e_d_rvz += d_e_d_g * d_g_d_rvz;
			
			double val = RuRv / sqrt(RuRu*RvRv);

			double d_val_d_rux = d_RuRv_d_rux / sqrt(RuRu*RvRv) -  d_RuRu_d_rux * RvRv*RuRv/(2*pow(RuRu*RvRv,3.0/2.0));
			double d_val_d_ruy = d_RuRv_d_ruy / sqrt(RuRu*RvRv) -  d_RuRu_d_ruy * RvRv*RuRv/(2*pow(RuRu*RvRv,3.0/2.0));
			double d_val_d_ruz = d_RuRv_d_ruz / sqrt(RuRu*RvRv) -  d_RuRu_d_ruz * RvRv*RuRv/(2*pow(RuRu*RvRv,3.0/2.0));
                                                                                                                                                                                      
			double d_val_d_rvx = d_RuRv_d_rvx / sqrt(RuRu*RvRv) -  d_RvRv_d_rvx * RuRu*RuRv/(2*pow(RuRu*RvRv,3.0/2.0));
			double d_val_d_rvy = d_RuRv_d_rvy / sqrt(RuRu*RvRv) -  d_RvRv_d_rvy * RuRu*RuRv/(2*pow(RuRu*RvRv,3.0/2.0));
			double d_val_d_rvz = d_RuRv_d_rvz / sqrt(RuRu*RvRv) -  d_RvRv_d_rvz * RuRu*RuRv/(2*pow(RuRu*RvRv,3.0/2.0));
			

			for( int p = 0; p < np; p++ )
			{
				gr[3*cp[p]+0] += d_e_d_rux * theFormulas[frm].r_u[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_ruy * theFormulas[frm].r_u[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_ruz * theFormulas[frm].r_u[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_rvx * theFormulas[frm].r_v[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_rvy * theFormulas[frm].r_v[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_rvz * theFormulas[frm].r_v[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_suux * theFormulas[frm].r_uu[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_suuy * theFormulas[frm].r_uu[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_suuz * theFormulas[frm].r_uu[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_suvx * theFormulas[frm].r_uv[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_suvy * theFormulas[frm].r_uv[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_suvz * theFormulas[frm].r_uv[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_svvx * theFormulas[frm].r_vv[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_svvy * theFormulas[frm].r_vv[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_svvz * theFormulas[frm].r_vv[p] * alpha_z; 
				

				gr[3*cp[p]+0] += k_reg * 2 * d_val_d_rux * (val-RuRv0) * theFormulas[frm].r_u[p] * alpha_x;
				gr[3*cp[p]+1] += k_reg * 2 * d_val_d_ruy * (val-RuRv0) * theFormulas[frm].r_u[p] * alpha_y;
				gr[3*cp[p]+2] += k_reg * 2 * d_val_d_ruz * (val-RuRv0) * theFormulas[frm].r_u[p] * alpha_z;
				                                         
				gr[3*cp[p]+0] += k_reg * 2 * d_val_d_rvx * (val-RuRv0) * theFormulas[frm].r_v[p] * alpha_x;
				gr[3*cp[p]+1] += k_reg * 2 * d_val_d_rvy * (val-RuRv0) * theFormulas[frm].r_v[p] * alpha_y;
				gr[3*cp[p]+2] += k_reg * 2 * d_val_d_rvz * (val-RuRv0) * theFormulas[frm].r_v[p] * alpha_z;
				
				dedr[3*p+0] += k_reg * 2 * d_val_d_rux * (val-RuRv0) * theFormulas[frm].r_u[p];
				dedr[3*p+1] += k_reg * 2 * d_val_d_ruy * (val-RuRv0) * theFormulas[frm].r_u[p];
				dedr[3*p+2] += k_reg * 2 * d_val_d_ruz * (val-RuRv0) * theFormulas[frm].r_u[p];
				                                       
				dedr[3*p+0] += k_reg * 2 * d_val_d_rvx * (val-RuRv0) * theFormulas[frm].r_v[p];
				dedr[3*p+1] += k_reg * 2 * d_val_d_rvy * (val-RuRv0) * theFormulas[frm].r_v[p];
				dedr[3*p+2] += k_reg * 2 * d_val_d_rvz * (val-RuRv0) * theFormulas[frm].r_v[p];

				dedr[3*p+0] += d_e_d_rux * theFormulas[frm].r_u[p]; 
				dedr[3*p+1] += d_e_d_ruy * theFormulas[frm].r_u[p]; 
				dedr[3*p+2] += d_e_d_ruz * theFormulas[frm].r_u[p]; 
				
				dedr[3*p+0] += d_e_d_rvx * theFormulas[frm].r_v[p]; 
				dedr[3*p+1] += d_e_d_rvy * theFormulas[frm].r_v[p]; 
				dedr[3*p+2] += d_e_d_rvz * theFormulas[frm].r_v[p]; 
				
				dedr[3*p+0] += d_e_d_suux * theFormulas[frm].r_uu[p]; 
				dedr[3*p+1] += d_e_d_suuy * theFormulas[frm].r_uu[p]; 
				dedr[3*p+2] += d_e_d_suuz * theFormulas[frm].r_uu[p]; 
				
				dedr[3*p+0] += d_e_d_suvx * theFormulas[frm].r_uv[p]; 
				dedr[3*p+1] += d_e_d_suvy * theFormulas[frm].r_uv[p]; 
				dedr[3*p+2] += d_e_d_suvz * theFormulas[frm].r_uv[p]; 
				
				dedr[3*p+0] += d_e_d_svvx * theFormulas[frm].r_vv[p]; 
				dedr[3*p+1] += d_e_d_svvy * theFormulas[frm].r_vv[p]; 
				dedr[3*p+2] += d_e_d_svvz * theFormulas[frm].r_vv[p]; 


			}
			
			for( int p = 0; p < np; p++ )
			{
				gr[3*nv+0] += dedr[3*p+0] * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				gr[3*nv+1] += dedr[3*p+1] * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				gr[3*nv+2] += dedr[3*p+2] * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			}

		}	

	} 

	r_val /= area;

	igrad( r, gr );

	pgrad( r, gr, puv, pg );

}
/*
void surface::igrad( double *r, double *gr )
{	
	double alpha_x = r[3*nv];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];
	if( !theIrregularFormulas )
		generatePlan();

	double r_val = 0;
	double e = 0;
	double area = 0;
	double wgt = 0;
		double dudv = 1.0;

//	printf("This gradient needs to be corrected to use the total area on an irregular face instead of the individual integration point area.\n");
//	exit(1);

	for( int f = 0; f < nf_irr_faces; f++ )
	{
		double p_face_area = 0;

		int frm = f*nf_irr_pts;
		int t = theIrregularFormulas[frm].tri;
		double g0 = theIrregularFormulas[frm].g0;//2*f_g0[nf_faces+f];
	
		for( int px = 0; px < theTriangles[t].np; px++ )
		{
			if( fabs(theTriangles[t].pc0[px] - theIrregularFormulas[frm].c0) < 1e-7 )
				continue; 
			p_face_area += theTriangles[t].pa[px];
		}
			if( debug_bit )
			{
			if( f != do_face )
				continue;
			}	
		double c0 = 0;
		double c1 = 0;
		double face_area = 0;
		double energy_density = 0;
		double A=0,A0=0;
		for( int p = 0; p < nf_irr_pts; p++ )
		{
			int frm = f*nf_irr_pts+p;
			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
			int *cp = theIrregularFormulas[f*nf_irr_pts+p].cp;
			int np = theIrregularFormulas[f*nf_irr_pts+p].ncoor;
			for( int p = 0; p < np; p++ )
			{

				R[0] += theIrregularFormulas[frm].r_w[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				R[1] += theIrregularFormulas[frm].r_w[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				R[2] += theIrregularFormulas[frm].r_w[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				Ru[0] += theIrregularFormulas[frm].r_u[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				Ru[1] += theIrregularFormulas[frm].r_u[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				Ru[2] += theIrregularFormulas[frm].r_u[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				Rv[0] += theIrregularFormulas[frm].r_v[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				Rv[1] += theIrregularFormulas[frm].r_v[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				Rv[2] += theIrregularFormulas[frm].r_v[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSuu[0] += theIrregularFormulas[frm].r_uu[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSuu[1] += theIrregularFormulas[frm].r_uu[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSuu[2] += theIrregularFormulas[frm].r_uu[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSuv[0] += theIrregularFormulas[frm].r_uv[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSuv[1] += theIrregularFormulas[frm].r_uv[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSuv[2] += theIrregularFormulas[frm].r_uv[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSvv[0] += theIrregularFormulas[frm].r_vv[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSvv[1] += theIrregularFormulas[frm].r_vv[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSvv[2] += theIrregularFormulas[frm].r_vv[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
			}

		double nrm[3];
		cross( Ru, Rv, nrm );
		normalize(nrm);

		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		double g = sqrt(RuRu*RvRv-RuRv*RuRv);

		double nsuu = tSuu[0] * nrm[0] + tSuu[1] * nrm[1] + tSuu[2] * nrm[2];
		double nsuv = tSuv[0] * nrm[0] + tSuv[1] * nrm[1] + tSuv[2] * nrm[2];
		double nsvv = tSvv[0] * nrm[0] + tSvv[1] * nrm[1] + tSvv[2] * nrm[2];

		double Stot = (nsuu * RvRv + nsvv * RuRu -2*nsuv*RuRv)/(g*g);

//		Int[i*n+j] = 0.5 * kc * Stot * Stot * g; 
//		AInt[i*n+j] = g;

		double Sop[4] = { 1.0/(g*g) * (nsuu * RvRv  - nsuv * RuRv), (1.0/(g*g)) * ( nsuv*RvRv-nsvv*RuRv),
				  1.0/(g*g) * (nsuv * RuRu  - nsuu * RuRv), (1.0/(g*g)) * ( nsvv*RuRu-nsuv*RuRv) };

		double a = Sop[0];		
		double b = Sop[1];
		double c = Sop[2];
		double d = Sop[3];

		double e1 = -0.5*(a+d-sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
		double e2 = -0.5*(a+d+sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
			
		double c0 = theIrregularFormulas[frm].c0;
//		printf("e1: %lf e2: %lf\n", e1, e2 );
		double en = 0.5 * kc * (e1 + e2 - c0 ) * (e1+e2 - c0);
		double dA = g * theIrregularFormulas[frm].weight;
		
			A +=  g * theIrregularFormulas[frm].weight;
			A0 += g0 * theIrregularFormulas[frm].weight;

			face_area  += dudv * g                   * theIrregularFormulas[frm].weight;	
			energy_density += 0.5 * kc * ( e1 + e2 - c0 ) * ( e1 + e2 - c0 ) * g * theIrregularFormulas[frm].weight;
		}

		

		for( int p = 0; p < nf_irr_pts; p++ )
		{
			if( debug_bit )
			{
			if( f != do_face || p != do_pt )
				continue;
			}	
			int frm = f*nf_irr_pts+p;
			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
			double nrm[3]={0,0,0}; 

			int *cp = theIrregularFormulas[f*nf_irr_pts+p].cp;
			int np = theIrregularFormulas[f*nf_irr_pts+p].ncoor;
			double dedr[3*np];
			memset( dedr, 0, sizeof(double) * 3 *np );

			for( int p = 0; p < np; p++ )
			{
				R[0] += theIrregularFormulas[frm].r_w[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				R[1] += theIrregularFormulas[frm].r_w[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				R[2] += theIrregularFormulas[frm].r_w[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				Ru[0] += theIrregularFormulas[frm].r_u[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				Ru[1] += theIrregularFormulas[frm].r_u[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				Ru[2] += theIrregularFormulas[frm].r_u[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				Rv[0] += theIrregularFormulas[frm].r_v[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				Rv[1] += theIrregularFormulas[frm].r_v[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				Rv[2] += theIrregularFormulas[frm].r_v[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSuu[0] += theIrregularFormulas[frm].r_uu[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSuu[1] += theIrregularFormulas[frm].r_uu[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSuu[2] += theIrregularFormulas[frm].r_uu[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSuv[0] += theIrregularFormulas[frm].r_uv[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSuv[1] += theIrregularFormulas[frm].r_uv[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSuv[2] += theIrregularFormulas[frm].r_uv[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSvv[0] += theIrregularFormulas[frm].r_vv[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSvv[1] += theIrregularFormulas[frm].r_vv[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSvv[2] += theIrregularFormulas[frm].r_vv[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
			}
		



		cross( Ru, Rv, nrm );
		normalize(nrm);

		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		double g = sqrt(RuRu*RvRv-RuRv*RuRv);

		double nsuu = tSuu[0] * nrm[0] + tSuu[1] * nrm[1] + tSuu[2] * nrm[2];
		double nsuv = tSuv[0] * nrm[0] + tSuv[1] * nrm[1] + tSuv[2] * nrm[2];
		double nsvv = tSvv[0] * nrm[0] + tSvv[1] * nrm[1] + tSvv[2] * nrm[2];

		double Stot = (nsuu * RvRv + nsvv * RuRu -2*nsuv*RuRv)/(g*g);

//		Int[i*n+j] = 0.5 * kc * Stot * Stot * g; 
//		AInt[i*n+j] = g;

		double Sop[4] = { 1.0/(g*g) * (nsuu * RvRv  - nsuv * RuRv), (1.0/(g*g)) * ( nsuv*RvRv-nsvv*RuRv),
				  1.0/(g*g) * (nsuv * RuRu  - nsuu * RuRv), (1.0/(g*g)) * ( nsvv*RuRu-nsuv*RuRv) };

		double a = Sop[0];		
		double b = Sop[1];
		double c = Sop[2];
		double d = Sop[3];

		double e1 = -0.5*(a+d-sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
		double e2 = -0.5*(a+d+sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
			
		double c0 = theIrregularFormulas[frm].c0;
//		printf("e1: %lf e2: %lf\n", e1, e2 );
		double en = 0.5 * kc * (e1 + e2 - c0 ) * (e1+e2 - c0);
		double dA = g * theIrregularFormulas[frm].weight;
		
//		A +=  g * theIrregularFormulas[frm].weight;
//		A0 += g0 * theIrregularFormulas[frm].weight;

		r_val += dudv * dA * (e1+e2)/2;
		area += dudv * dA;
		e += dudv * dA * en; 			
double d_nrmz_d_rux=0,d_nrmz_d_ruy=0,d_nrmz_d_ruz=0,d_nrmz_d_rvx=0,d_nrmz_d_rvy=0,d_nrmz_d_rvz=0,d_nrmy_d_rux=0,d_nrmy_d_ruy=0,d_nrmy_d_ruz=0,d_nrmy_d_rvx=0,d_nrmy_d_rvy=0,d_nrmy_d_rvz=0,d_nrmx_d_rux=0,d_nrmx_d_ruy=0,d_nrmx_d_ruz=0,d_nrmx_d_rvx=0,d_nrmx_d_rvy=0,d_nrmx_d_rvz=0,d_RuRu_d_rux=0,d_RuRu_d_ruy=0,d_RuRu_d_ruz=0,d_nsvv_d_rux=0,d_nsvv_d_ruy=0,d_nsvv_d_ruz=0,d_nsvv_d_rvx=0,d_nsvv_d_rvy=0,d_nsvv_d_rvz=0,d_nsvv_d_svvx=0,d_nsvv_d_svvy=0,d_nsvv_d_svvz=0,d_RuRv_d_rux=0,d_RuRv_d_ruy=0,d_RuRv_d_ruz=0,d_RuRv_d_rvx=0,d_RuRv_d_rvy=0,d_RuRv_d_rvz=0,d_nsuv_d_rux=0,d_nsuv_d_ruy=0,d_nsuv_d_ruz=0,d_nsuv_d_rvx=0,d_nsuv_d_rvy=0,d_nsuv_d_rvz=0,d_nsuv_d_suvx=0,d_nsuv_d_suvy=0,d_nsuv_d_suvz=0,d_RvRv_d_rvx=0,d_RvRv_d_rvy=0,d_RvRv_d_rvz=0,d_nsuu_d_rux=0,d_nsuu_d_ruy=0,d_nsuu_d_ruz=0,d_nsuu_d_rvx=0,d_nsuu_d_rvy=0,d_nsuu_d_rvz=0,d_nsuu_d_suux=0,d_nsuu_d_suuy=0,d_nsuu_d_suuz=0,d_g_d_rux=0,d_g_d_ruy=0,d_g_d_ruz=0,d_g_d_rvx=0,d_g_d_rvy=0,d_g_d_rvz=0,d_e_d_rux=0,d_e_d_ruy=0,d_e_d_ruz=0,d_e_d_rvx=0,d_e_d_rvy=0,d_e_d_rvz=0,d_e_d_suux=0,d_e_d_suuy=0,d_e_d_suuz=0,d_e_d_suvx=0,d_e_d_suvy=0,d_e_d_suvz=0,d_e_d_svvx=0,d_e_d_svvy=0,d_e_d_svvz=0,d_c2_d_rux=0,d_c2_d_ruy=0,d_c2_d_ruz=0,d_c2_d_rvx=0,d_c2_d_rvy=0,d_c2_d_rvz=0,d_c2_d_suux=0,d_c2_d_suuy=0,d_c2_d_suuz=0,d_c2_d_suvx=0,d_c2_d_suvy=0,d_c2_d_suvz=0,d_c2_d_svvx=0,d_c2_d_svvy=0,d_c2_d_svvz=0,d_d_d_rux=0,d_d_d_ruy=0,d_d_d_ruz=0,d_d_d_rvx=0,d_d_d_rvy=0,d_d_d_rvz=0,d_d_d_suvx=0,d_d_d_suvy=0,d_d_d_suvz=0,d_d_d_svvx=0,d_d_d_svvy=0,d_d_d_svvz=0,d_c_d_rux=0,d_c_d_ruy=0,d_c_d_ruz=0,d_c_d_rvx=0,d_c_d_rvy=0,d_c_d_rvz=0,d_c_d_suux=0,d_c_d_suuy=0,d_c_d_suuz=0,d_c_d_suvx=0,d_c_d_suvy=0,d_c_d_suvz=0,d_b_d_rux=0,d_b_d_ruy=0,d_b_d_ruz=0,d_b_d_rvx=0,d_b_d_rvy=0,d_b_d_rvz=0,d_b_d_suvx=0,d_b_d_suvy=0,d_b_d_suvz=0,d_b_d_svvx=0,d_b_d_svvy=0,d_b_d_svvz=0,d_a_d_rux=0,d_a_d_ruy=0,d_a_d_ruz=0,d_a_d_rvx=0,d_a_d_rvy=0,d_a_d_rvz=0,d_a_d_suux=0,d_a_d_suuy=0,d_a_d_suuz=0,d_a_d_suvx=0,d_a_d_suvy=0,d_a_d_suvz=0,d_c1_d_rux=0,d_c1_d_ruy=0,d_c1_d_ruz=0,d_c1_d_rvx=0,d_c1_d_rvy=0,d_c1_d_rvz=0,d_c1_d_suux=0,d_c1_d_suuy=0,d_c1_d_suuz=0,d_c1_d_suvx=0,d_c1_d_suvy=0,d_c1_d_suvz=0,d_c1_d_svvx=0,d_c1_d_svvy=0,d_c1_d_svvz=0,d_g_d_RvRv=0,d_g_d_RuRv=0,d_g_d_RuRu=0,d_nsvv_d_nrmz=0,d_nsvv_d_nrmy=0,d_nsvv_d_nrmx=0,d_nsuv_d_nrmz=0,d_nsuv_d_nrmy=0,d_nsuv_d_nrmx=0,d_nsuu_d_nrmz=0,d_nsuu_d_nrmy=0,d_nsuu_d_nrmx=0,d_d_d_RuRv=0,d_d_d_nsuv=0,d_d_d_RuRu=0,d_d_d_nsvv=0,d_d_d_g=0,d_c_d_RuRv=0,d_c_d_nsuv=0,d_c_d_RuRu=0,d_c_d_nsuu=0,d_c_d_g=0,d_b_d_RuRv=0,d_b_d_nsuv=0,d_b_d_RvRv=0,d_b_d_nsvv=0,d_b_d_g=0,d_a_d_RuRv=0,d_a_d_nsuv=0,d_a_d_RvRv=0,d_a_d_nsuu=0,d_a_d_g=0,d_e_d_c2=0,d_e_d_c1=0,d_e_d_g=0,d_c2_d_d=0,d_c2_d_c=0,d_c2_d_b=0,d_c2_d_a=0,d_c1_d_d=0,d_c1_d_c=0,d_c1_d_b=0,d_c1_d_a=0,junk;

			// the basic variables.

			d_g_d_RuRu = RvRv/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
			d_g_d_RuRv = -(RuRv/Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
			d_g_d_RvRv = RuRu/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  

			d_RuRu_d_rux = 2*Ru[0];
			d_RuRu_d_ruy = 2*Ru[1];
			d_RuRu_d_ruz = 2*Ru[2];
			
			d_RuRv_d_rux = Rv[0];
			d_RuRv_d_ruy = Rv[1];
			d_RuRv_d_ruz = Rv[2];
			
			d_RuRv_d_rvx = Ru[0];
			d_RuRv_d_rvy = Ru[1];
			d_RuRv_d_rvz = Ru[2];
			
			d_RvRv_d_rvx = 2*Rv[0];
			d_RvRv_d_rvy = 2*Rv[1];
			d_RvRv_d_rvz = 2*Rv[2];

			// intermediates.
			d_a_d_g = (-2*(-(nsuv*RuRv) + nsuu*RvRv))/Power(g,3);
			d_a_d_nsuu = RvRv/Power(g,2);
			d_a_d_RvRv = nsuu/Power(g,2);
			d_a_d_nsuv = -(RuRv/Power(g,2));
			d_a_d_RuRv = -(nsuv/Power(g,2));

			d_b_d_g = (-2*(-(nsvv*RuRv) + nsuv*RvRv))/Power(g,3);
			d_b_d_nsvv = -(RuRv/Power(g,2));
			d_b_d_RvRv = nsuv/Power(g,2);
			d_b_d_nsuv = RvRv/Power(g,2);
			d_b_d_RuRv = -(nsvv/Power(g,2));

			d_c_d_g = (-2*(nsuv*RuRu - nsuu*RuRv))/Power(g,3);
			d_c_d_nsuu = -(RuRv/Power(g,2));
			d_c_d_RuRu = nsuv/Power(g,2);
			d_c_d_nsuv = RuRu/Power(g,2);
			d_c_d_RuRv = -(nsuu/Power(g,2));

			d_d_d_g = (-2*(nsvv*RuRu - nsuv*RuRv))/Power(g,3);
			d_d_d_nsvv = RuRu/Power(g,2);
			d_d_d_RuRu = nsvv/Power(g,2);
			d_d_d_nsuv = -(RuRv/Power(g,2));
			d_d_d_RuRv = -(nsuv/Power(g,2));

			d_c1_d_a = -0.5*(1 - (2*a - 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
			d_c1_d_b =-(-1.*c)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2)) ;
			d_c1_d_c =-(-1.*b)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2)) ;
			d_c1_d_d =-0.5*(1 - (-2*a + 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
		
			d_c2_d_a = -0.5*(1 + (2*a - 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
			d_c2_d_b = -(1.*c)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2));
			d_c2_d_c = -(1.*b)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2));
			d_c2_d_d = -0.5*(1 + (-2*a + 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));

			//d_e_d_g  = 0.5 * kc * ( e1 + e2 - c0) * ( e1 + e2 - c0) * dudv * theIrregularFormulas[frm].weight;
#ifdef FIXED_A
			d_e_d_g = 0;
			d_e_d_c1 = kc * (e1 + e2 - c0) * dudv * g0 * theIrregularFormulas[frm].weight;
			d_e_d_c2 = kc * (e1 + e2 - c0) * dudv * g0 * theIrregularFormulas[frm].weight;
#else
			d_e_d_g  = (0.5 * kc * ( e1 + e2 - c0) * ( e1 + e2 - c0) + kg * e1 * e2 ) * dudv * theIrregularFormulas[frm].weight;
//			if( g < 0 )
//				d_e_d_g *= -1;

			d_e_d_c1 = kc * (e1 + e2 - c0) * dudv * g * theIrregularFormulas[frm].weight;
			d_e_d_c2 = kc * (e1 + e2 - c0) * dudv * g * theIrregularFormulas[frm].weight;

			d_e_d_c1 += kg * e2 * dudv * g * theIrregularFormulas[frm].weight;
			d_e_d_c2 += kg * e1 * dudv * g * theIrregularFormulas[frm].weight;
#endif

			d_e_d_g += 2* KA * ((A-A0)/A0) * theIrregularFormulas[frm].weight * dudv;

			
			// derivative of dA in the numerator:
			d_e_d_g  += -p_face_area *  0.5 * kc * (e1+e2-c0)*(e1+e2-c0) * theIrregularFormulas[frm].weight / face_area;
			// derivative of dA in the denominator:
			d_e_d_g  += +p_face_area * energy_density * theIrregularFormulas[frm].weight / face_area / face_area;
			d_e_d_c1 += -p_face_area *  dA *  kc * (e1+e2-c0) / face_area;
			d_e_d_c2 += -p_face_area *  dA *  kc * (e1+e2-c0) / face_area;

//			double fac = Power(Power(Ru[1]*Rv[0] - Ru[0]*Rv[1],2) + Power(Ru[2]*Rv[0] - Ru[0]*Rv[2],2) + Power(Ru[2]*Rv[1] - Ru[1]*Rv[2],2),1.5);
			double fac = Power(Power(Ru[1]*Rv[0] - Ru[0]*Rv[1],2) + Power(Ru[2]*Rv[0] - Ru[0]*Rv[2],2) + Power(Ru[2]*Rv[1] - Ru[1]*Rv[2],2),1.5);

			d_nrmx_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac; 
			d_nrmx_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2)))/fac;
			d_nrmx_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac;
			
			d_nrmy_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2))))/fac;
			d_nrmy_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac; 
			d_nrmy_d_ruz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2])) + Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac;
	
			d_nrmz_d_rux = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2])/fac;
			d_nrmz_d_ruy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;
			d_nrmz_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2))) + (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;

			d_nrmx_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac; 
			d_nrmx_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2])))/fac;
			d_nrmx_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac;
			
			d_nrmy_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
			d_nrmy_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2])))/fac;
			d_nrmy_d_rvz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[0]*Ru[1]*Rv[0] - Power(Ru[0],2)*Rv[1] + Ru[2]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2])))/fac;
	
			d_nrmz_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;
			d_nrmz_d_rvy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[0]*Ru[2]*Rv[0] - Power(Ru[0],2)*Rv[2] + Ru[1]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
			d_nrmz_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;
			d_nsuu_d_nrmx = tSuu[0];
			d_nsuu_d_nrmy = tSuu[1];
			d_nsuu_d_nrmz = tSuu[2];
			
			d_nsuv_d_nrmx = tSuv[0];
			d_nsuv_d_nrmy = tSuv[1];
			d_nsuv_d_nrmz = tSuv[2];

			d_nsvv_d_nrmx = tSvv[0];
			d_nsvv_d_nrmy = tSvv[1];
			d_nsvv_d_nrmz = tSvv[2];

			d_nsuu_d_suux = nrm[0];
			d_nsuu_d_suuy = nrm[1];
			d_nsuu_d_suuz = nrm[2];
			
			d_nsuv_d_suvx = nrm[0];
			d_nsuv_d_suvy = nrm[1];
			d_nsuv_d_suvz = nrm[2];

			d_nsvv_d_svvx = nrm[0];
			d_nsvv_d_svvy = nrm[1];
			d_nsvv_d_svvz = nrm[2];



	d_nsvv_d_rux += d_nsvv_d_nrmz * d_nrmz_d_rux;
	d_nsvv_d_ruy += d_nsvv_d_nrmz * d_nrmz_d_ruy;
	d_nsvv_d_ruz += d_nsvv_d_nrmz * d_nrmz_d_ruz;
	d_nsvv_d_rvx += d_nsvv_d_nrmz * d_nrmz_d_rvx;
	d_nsvv_d_rvy += d_nsvv_d_nrmz * d_nrmz_d_rvy;
	d_nsvv_d_rvz += d_nsvv_d_nrmz * d_nrmz_d_rvz;
	d_nsvv_d_rux += d_nsvv_d_nrmy * d_nrmy_d_rux;
	d_nsvv_d_ruy += d_nsvv_d_nrmy * d_nrmy_d_ruy;
	d_nsvv_d_ruz += d_nsvv_d_nrmy * d_nrmy_d_ruz;
	d_nsvv_d_rvx += d_nsvv_d_nrmy * d_nrmy_d_rvx;
	d_nsvv_d_rvy += d_nsvv_d_nrmy * d_nrmy_d_rvy;
	d_nsvv_d_rvz += d_nsvv_d_nrmy * d_nrmy_d_rvz;
	d_nsvv_d_rux += d_nsvv_d_nrmx * d_nrmx_d_rux;
	d_nsvv_d_ruy += d_nsvv_d_nrmx * d_nrmx_d_ruy;
	d_nsvv_d_ruz += d_nsvv_d_nrmx * d_nrmx_d_ruz;
	d_nsvv_d_rvx += d_nsvv_d_nrmx * d_nrmx_d_rvx;
	d_nsvv_d_rvy += d_nsvv_d_nrmx * d_nrmx_d_rvy;
	d_nsvv_d_rvz += d_nsvv_d_nrmx * d_nrmx_d_rvz;
	d_nsuv_d_rux += d_nsuv_d_nrmz * d_nrmz_d_rux;
	d_nsuv_d_ruy += d_nsuv_d_nrmz * d_nrmz_d_ruy;
	d_nsuv_d_ruz += d_nsuv_d_nrmz * d_nrmz_d_ruz;
	d_nsuv_d_rvx += d_nsuv_d_nrmz * d_nrmz_d_rvx;
	d_nsuv_d_rvy += d_nsuv_d_nrmz * d_nrmz_d_rvy;
	d_nsuv_d_rvz += d_nsuv_d_nrmz * d_nrmz_d_rvz;
	d_nsuv_d_rux += d_nsuv_d_nrmy * d_nrmy_d_rux;
	d_nsuv_d_ruy += d_nsuv_d_nrmy * d_nrmy_d_ruy;
	d_nsuv_d_ruz += d_nsuv_d_nrmy * d_nrmy_d_ruz;
	d_nsuv_d_rvx += d_nsuv_d_nrmy * d_nrmy_d_rvx;
	d_nsuv_d_rvy += d_nsuv_d_nrmy * d_nrmy_d_rvy;
	d_nsuv_d_rvz += d_nsuv_d_nrmy * d_nrmy_d_rvz;
	d_nsuv_d_rux += d_nsuv_d_nrmx * d_nrmx_d_rux;
	d_nsuv_d_ruy += d_nsuv_d_nrmx * d_nrmx_d_ruy;
	d_nsuv_d_ruz += d_nsuv_d_nrmx * d_nrmx_d_ruz;
	d_nsuv_d_rvx += d_nsuv_d_nrmx * d_nrmx_d_rvx;
	d_nsuv_d_rvy += d_nsuv_d_nrmx * d_nrmx_d_rvy;
	d_nsuv_d_rvz += d_nsuv_d_nrmx * d_nrmx_d_rvz;
	d_nsuu_d_rux += d_nsuu_d_nrmz * d_nrmz_d_rux;
	d_nsuu_d_ruy += d_nsuu_d_nrmz * d_nrmz_d_ruy;
	d_nsuu_d_ruz += d_nsuu_d_nrmz * d_nrmz_d_ruz;
	d_nsuu_d_rvx += d_nsuu_d_nrmz * d_nrmz_d_rvx;
	d_nsuu_d_rvy += d_nsuu_d_nrmz * d_nrmz_d_rvy;
	d_nsuu_d_rvz += d_nsuu_d_nrmz * d_nrmz_d_rvz;
	d_nsuu_d_rux += d_nsuu_d_nrmy * d_nrmy_d_rux;
	d_nsuu_d_ruy += d_nsuu_d_nrmy * d_nrmy_d_ruy;
	d_nsuu_d_ruz += d_nsuu_d_nrmy * d_nrmy_d_ruz;
	d_nsuu_d_rvx += d_nsuu_d_nrmy * d_nrmy_d_rvx;
	d_nsuu_d_rvy += d_nsuu_d_nrmy * d_nrmy_d_rvy;
	d_nsuu_d_rvz += d_nsuu_d_nrmy * d_nrmy_d_rvz;
	d_nsuu_d_rux += d_nsuu_d_nrmx * d_nrmx_d_rux;
	d_nsuu_d_ruy += d_nsuu_d_nrmx * d_nrmx_d_ruy;
	d_nsuu_d_ruz += d_nsuu_d_nrmx * d_nrmx_d_ruz;
	d_nsuu_d_rvx += d_nsuu_d_nrmx * d_nrmx_d_rvx;
	d_nsuu_d_rvy += d_nsuu_d_nrmx * d_nrmx_d_rvy;
	d_nsuu_d_rvz += d_nsuu_d_nrmx * d_nrmx_d_rvz;
	d_g_d_rvx += d_g_d_RvRv * d_RvRv_d_rvx;
	d_g_d_rvy += d_g_d_RvRv * d_RvRv_d_rvy;
	d_g_d_rvz += d_g_d_RvRv * d_RvRv_d_rvz;
	d_g_d_rux += d_g_d_RuRv * d_RuRv_d_rux;
	d_g_d_ruy += d_g_d_RuRv * d_RuRv_d_ruy;
	d_g_d_ruz += d_g_d_RuRv * d_RuRv_d_ruz;
	d_g_d_rvx += d_g_d_RuRv * d_RuRv_d_rvx;
	d_g_d_rvy += d_g_d_RuRv * d_RuRv_d_rvy;
	d_g_d_rvz += d_g_d_RuRv * d_RuRv_d_rvz;
	d_g_d_rux += d_g_d_RuRu * d_RuRu_d_rux;
	d_g_d_ruy += d_g_d_RuRu * d_RuRu_d_ruy;
	d_g_d_ruz += d_g_d_RuRu * d_RuRu_d_ruz;
	d_d_d_rux += d_d_d_RuRv * d_RuRv_d_rux;
	d_d_d_ruy += d_d_d_RuRv * d_RuRv_d_ruy;
	d_d_d_ruz += d_d_d_RuRv * d_RuRv_d_ruz;
	d_d_d_rvx += d_d_d_RuRv * d_RuRv_d_rvx;
	d_d_d_rvy += d_d_d_RuRv * d_RuRv_d_rvy;
	d_d_d_rvz += d_d_d_RuRv * d_RuRv_d_rvz;
	d_d_d_rux += d_d_d_nsuv * d_nsuv_d_rux;
	d_d_d_ruy += d_d_d_nsuv * d_nsuv_d_ruy;
	d_d_d_ruz += d_d_d_nsuv * d_nsuv_d_ruz;
	d_d_d_rvx += d_d_d_nsuv * d_nsuv_d_rvx;
	d_d_d_rvy += d_d_d_nsuv * d_nsuv_d_rvy;
	d_d_d_rvz += d_d_d_nsuv * d_nsuv_d_rvz;
	d_d_d_suvx += d_d_d_nsuv * d_nsuv_d_suvx;
	d_d_d_suvy += d_d_d_nsuv * d_nsuv_d_suvy;
	d_d_d_suvz += d_d_d_nsuv * d_nsuv_d_suvz;
	d_d_d_rux += d_d_d_RuRu * d_RuRu_d_rux;
	d_d_d_ruy += d_d_d_RuRu * d_RuRu_d_ruy;
	d_d_d_ruz += d_d_d_RuRu * d_RuRu_d_ruz;
	d_d_d_rux += d_d_d_nsvv * d_nsvv_d_rux;
	d_d_d_ruy += d_d_d_nsvv * d_nsvv_d_ruy;
	d_d_d_ruz += d_d_d_nsvv * d_nsvv_d_ruz;
	d_d_d_rvx += d_d_d_nsvv * d_nsvv_d_rvx;
	d_d_d_rvy += d_d_d_nsvv * d_nsvv_d_rvy;
	d_d_d_rvz += d_d_d_nsvv * d_nsvv_d_rvz;
	d_d_d_svvx += d_d_d_nsvv * d_nsvv_d_svvx;
	d_d_d_svvy += d_d_d_nsvv * d_nsvv_d_svvy;
	d_d_d_svvz += d_d_d_nsvv * d_nsvv_d_svvz;
	d_d_d_rux += d_d_d_g * d_g_d_rux;
	d_d_d_ruy += d_d_d_g * d_g_d_ruy;
	d_d_d_ruz += d_d_d_g * d_g_d_ruz;
	d_d_d_rvx += d_d_d_g * d_g_d_rvx;
	d_d_d_rvy += d_d_d_g * d_g_d_rvy;
	d_d_d_rvz += d_d_d_g * d_g_d_rvz;
	d_c_d_rux += d_c_d_RuRv * d_RuRv_d_rux;
	d_c_d_ruy += d_c_d_RuRv * d_RuRv_d_ruy;
	d_c_d_ruz += d_c_d_RuRv * d_RuRv_d_ruz;
	d_c_d_rvx += d_c_d_RuRv * d_RuRv_d_rvx;
	d_c_d_rvy += d_c_d_RuRv * d_RuRv_d_rvy;
	d_c_d_rvz += d_c_d_RuRv * d_RuRv_d_rvz;
	d_c_d_rux += d_c_d_nsuv * d_nsuv_d_rux;
	d_c_d_ruy += d_c_d_nsuv * d_nsuv_d_ruy;
	d_c_d_ruz += d_c_d_nsuv * d_nsuv_d_ruz;
	d_c_d_rvx += d_c_d_nsuv * d_nsuv_d_rvx;
	d_c_d_rvy += d_c_d_nsuv * d_nsuv_d_rvy;
	d_c_d_rvz += d_c_d_nsuv * d_nsuv_d_rvz;
	d_c_d_suvx += d_c_d_nsuv * d_nsuv_d_suvx;
	d_c_d_suvy += d_c_d_nsuv * d_nsuv_d_suvy;
	d_c_d_suvz += d_c_d_nsuv * d_nsuv_d_suvz;
	d_c_d_rux += d_c_d_RuRu * d_RuRu_d_rux;
	d_c_d_ruy += d_c_d_RuRu * d_RuRu_d_ruy;
	d_c_d_ruz += d_c_d_RuRu * d_RuRu_d_ruz;
	d_c_d_rux += d_c_d_nsuu * d_nsuu_d_rux;
	d_c_d_ruy += d_c_d_nsuu * d_nsuu_d_ruy;
	d_c_d_ruz += d_c_d_nsuu * d_nsuu_d_ruz;
	d_c_d_rvx += d_c_d_nsuu * d_nsuu_d_rvx;
	d_c_d_rvy += d_c_d_nsuu * d_nsuu_d_rvy;
	d_c_d_rvz += d_c_d_nsuu * d_nsuu_d_rvz;
	d_c_d_suux += d_c_d_nsuu * d_nsuu_d_suux;
	d_c_d_suuy += d_c_d_nsuu * d_nsuu_d_suuy;
	d_c_d_suuz += d_c_d_nsuu * d_nsuu_d_suuz;
	d_c_d_rux += d_c_d_g * d_g_d_rux;
	d_c_d_ruy += d_c_d_g * d_g_d_ruy;
	d_c_d_ruz += d_c_d_g * d_g_d_ruz;
	d_c_d_rvx += d_c_d_g * d_g_d_rvx;
	d_c_d_rvy += d_c_d_g * d_g_d_rvy;
	d_c_d_rvz += d_c_d_g * d_g_d_rvz;
	d_b_d_rux += d_b_d_RuRv * d_RuRv_d_rux;
	d_b_d_ruy += d_b_d_RuRv * d_RuRv_d_ruy;
	d_b_d_ruz += d_b_d_RuRv * d_RuRv_d_ruz;
	d_b_d_rvx += d_b_d_RuRv * d_RuRv_d_rvx;
	d_b_d_rvy += d_b_d_RuRv * d_RuRv_d_rvy;
	d_b_d_rvz += d_b_d_RuRv * d_RuRv_d_rvz;
	d_b_d_rux += d_b_d_nsuv * d_nsuv_d_rux;
	d_b_d_ruy += d_b_d_nsuv * d_nsuv_d_ruy;
	d_b_d_ruz += d_b_d_nsuv * d_nsuv_d_ruz;
	d_b_d_rvx += d_b_d_nsuv * d_nsuv_d_rvx;
	d_b_d_rvy += d_b_d_nsuv * d_nsuv_d_rvy;
	d_b_d_rvz += d_b_d_nsuv * d_nsuv_d_rvz;
	d_b_d_suvx += d_b_d_nsuv * d_nsuv_d_suvx;
	d_b_d_suvy += d_b_d_nsuv * d_nsuv_d_suvy;
	d_b_d_suvz += d_b_d_nsuv * d_nsuv_d_suvz;
	d_b_d_rvx += d_b_d_RvRv * d_RvRv_d_rvx;
	d_b_d_rvy += d_b_d_RvRv * d_RvRv_d_rvy;
	d_b_d_rvz += d_b_d_RvRv * d_RvRv_d_rvz;
	d_b_d_rux += d_b_d_nsvv * d_nsvv_d_rux;
	d_b_d_ruy += d_b_d_nsvv * d_nsvv_d_ruy;
	d_b_d_ruz += d_b_d_nsvv * d_nsvv_d_ruz;
	d_b_d_rvx += d_b_d_nsvv * d_nsvv_d_rvx;
	d_b_d_rvy += d_b_d_nsvv * d_nsvv_d_rvy;
	d_b_d_rvz += d_b_d_nsvv * d_nsvv_d_rvz;
	d_b_d_svvx += d_b_d_nsvv * d_nsvv_d_svvx;
	d_b_d_svvy += d_b_d_nsvv * d_nsvv_d_svvy;
	d_b_d_svvz += d_b_d_nsvv * d_nsvv_d_svvz;
	d_b_d_rux += d_b_d_g * d_g_d_rux;
	d_b_d_ruy += d_b_d_g * d_g_d_ruy;
	d_b_d_ruz += d_b_d_g * d_g_d_ruz;
	d_b_d_rvx += d_b_d_g * d_g_d_rvx;
	d_b_d_rvy += d_b_d_g * d_g_d_rvy;
	d_b_d_rvz += d_b_d_g * d_g_d_rvz;
	d_a_d_rux += d_a_d_RuRv * d_RuRv_d_rux;
	d_a_d_ruy += d_a_d_RuRv * d_RuRv_d_ruy;
	d_a_d_ruz += d_a_d_RuRv * d_RuRv_d_ruz;
	d_a_d_rvx += d_a_d_RuRv * d_RuRv_d_rvx;
	d_a_d_rvy += d_a_d_RuRv * d_RuRv_d_rvy;
	d_a_d_rvz += d_a_d_RuRv * d_RuRv_d_rvz;
	d_a_d_rux += d_a_d_nsuv * d_nsuv_d_rux;
	d_a_d_ruy += d_a_d_nsuv * d_nsuv_d_ruy;
	d_a_d_ruz += d_a_d_nsuv * d_nsuv_d_ruz;
	d_a_d_rvx += d_a_d_nsuv * d_nsuv_d_rvx;
	d_a_d_rvy += d_a_d_nsuv * d_nsuv_d_rvy;
	d_a_d_rvz += d_a_d_nsuv * d_nsuv_d_rvz;
	d_a_d_suvx += d_a_d_nsuv * d_nsuv_d_suvx;
	d_a_d_suvy += d_a_d_nsuv * d_nsuv_d_suvy;
	d_a_d_suvz += d_a_d_nsuv * d_nsuv_d_suvz;
	d_a_d_rvx += d_a_d_RvRv * d_RvRv_d_rvx;
	d_a_d_rvy += d_a_d_RvRv * d_RvRv_d_rvy;
	d_a_d_rvz += d_a_d_RvRv * d_RvRv_d_rvz;
	d_a_d_rux += d_a_d_nsuu * d_nsuu_d_rux;
	d_a_d_ruy += d_a_d_nsuu * d_nsuu_d_ruy;
	d_a_d_ruz += d_a_d_nsuu * d_nsuu_d_ruz;
	d_a_d_rvx += d_a_d_nsuu * d_nsuu_d_rvx;
	d_a_d_rvy += d_a_d_nsuu * d_nsuu_d_rvy;
	d_a_d_rvz += d_a_d_nsuu * d_nsuu_d_rvz;
	d_a_d_suux += d_a_d_nsuu * d_nsuu_d_suux;
	d_a_d_suuy += d_a_d_nsuu * d_nsuu_d_suuy;
	d_a_d_suuz += d_a_d_nsuu * d_nsuu_d_suuz;
	d_a_d_rux += d_a_d_g * d_g_d_rux;
	d_a_d_ruy += d_a_d_g * d_g_d_ruy;
	d_a_d_ruz += d_a_d_g * d_g_d_ruz;
	d_a_d_rvx += d_a_d_g * d_g_d_rvx;
	d_a_d_rvy += d_a_d_g * d_g_d_rvy;
	d_a_d_rvz += d_a_d_g * d_g_d_rvz;
	d_c1_d_rux += d_c1_d_d * d_d_d_rux;
	d_c1_d_ruy += d_c1_d_d * d_d_d_ruy;
	d_c1_d_ruz += d_c1_d_d * d_d_d_ruz;
	d_c1_d_rvx += d_c1_d_d * d_d_d_rvx;
	d_c1_d_rvy += d_c1_d_d * d_d_d_rvy;
	d_c1_d_rvz += d_c1_d_d * d_d_d_rvz;
	d_c1_d_suvx += d_c1_d_d * d_d_d_suvx;
	d_c1_d_suvy += d_c1_d_d * d_d_d_suvy;
	d_c1_d_suvz += d_c1_d_d * d_d_d_suvz;
	d_c1_d_svvx += d_c1_d_d * d_d_d_svvx;
	d_c1_d_svvy += d_c1_d_d * d_d_d_svvy;
	d_c1_d_svvz += d_c1_d_d * d_d_d_svvz;
	d_c1_d_rux += d_c1_d_c * d_c_d_rux;
	d_c1_d_ruy += d_c1_d_c * d_c_d_ruy;
	d_c1_d_ruz += d_c1_d_c * d_c_d_ruz;
	d_c1_d_rvx += d_c1_d_c * d_c_d_rvx;
	d_c1_d_rvy += d_c1_d_c * d_c_d_rvy;
	d_c1_d_rvz += d_c1_d_c * d_c_d_rvz;
	d_c1_d_suux += d_c1_d_c * d_c_d_suux;
	d_c1_d_suuy += d_c1_d_c * d_c_d_suuy;
	d_c1_d_suuz += d_c1_d_c * d_c_d_suuz;
	d_c1_d_suvx += d_c1_d_c * d_c_d_suvx;
	d_c1_d_suvy += d_c1_d_c * d_c_d_suvy;
	d_c1_d_suvz += d_c1_d_c * d_c_d_suvz;
	d_c1_d_rux += d_c1_d_b * d_b_d_rux;
	d_c1_d_ruy += d_c1_d_b * d_b_d_ruy;
	d_c1_d_ruz += d_c1_d_b * d_b_d_ruz;
	d_c1_d_rvx += d_c1_d_b * d_b_d_rvx;
	d_c1_d_rvy += d_c1_d_b * d_b_d_rvy;
	d_c1_d_rvz += d_c1_d_b * d_b_d_rvz;
	d_c1_d_suvx += d_c1_d_b * d_b_d_suvx;
	d_c1_d_suvy += d_c1_d_b * d_b_d_suvy;
	d_c1_d_suvz += d_c1_d_b * d_b_d_suvz;
	d_c1_d_svvx += d_c1_d_b * d_b_d_svvx;
	d_c1_d_svvy += d_c1_d_b * d_b_d_svvy;
	d_c1_d_svvz += d_c1_d_b * d_b_d_svvz;
	d_c1_d_rux += d_c1_d_a * d_a_d_rux;
	d_c1_d_ruy += d_c1_d_a * d_a_d_ruy;
	d_c1_d_ruz += d_c1_d_a * d_a_d_ruz;
	d_c1_d_rvx += d_c1_d_a * d_a_d_rvx;
	d_c1_d_rvy += d_c1_d_a * d_a_d_rvy;
	d_c1_d_rvz += d_c1_d_a * d_a_d_rvz;
	d_c1_d_suux += d_c1_d_a * d_a_d_suux;
	d_c1_d_suuy += d_c1_d_a * d_a_d_suuy;
	d_c1_d_suuz += d_c1_d_a * d_a_d_suuz;
	d_c1_d_suvx += d_c1_d_a * d_a_d_suvx;
	d_c1_d_suvy += d_c1_d_a * d_a_d_suvy;
	d_c1_d_suvz += d_c1_d_a * d_a_d_suvz;
	d_c2_d_rux += d_c2_d_d * d_d_d_rux;
	d_c2_d_ruy += d_c2_d_d * d_d_d_ruy;
	d_c2_d_ruz += d_c2_d_d * d_d_d_ruz;
	d_c2_d_rvx += d_c2_d_d * d_d_d_rvx;
	d_c2_d_rvy += d_c2_d_d * d_d_d_rvy;
	d_c2_d_rvz += d_c2_d_d * d_d_d_rvz;
	d_c2_d_suvx += d_c2_d_d * d_d_d_suvx;
	d_c2_d_suvy += d_c2_d_d * d_d_d_suvy;
	d_c2_d_suvz += d_c2_d_d * d_d_d_suvz;
	d_c2_d_svvx += d_c2_d_d * d_d_d_svvx;
	d_c2_d_svvy += d_c2_d_d * d_d_d_svvy;
	d_c2_d_svvz += d_c2_d_d * d_d_d_svvz;
	d_c2_d_rux += d_c2_d_c * d_c_d_rux;
	d_c2_d_ruy += d_c2_d_c * d_c_d_ruy;
	d_c2_d_ruz += d_c2_d_c * d_c_d_ruz;
	d_c2_d_rvx += d_c2_d_c * d_c_d_rvx;
	d_c2_d_rvy += d_c2_d_c * d_c_d_rvy;
	d_c2_d_rvz += d_c2_d_c * d_c_d_rvz;
	d_c2_d_suux += d_c2_d_c * d_c_d_suux;
	d_c2_d_suuy += d_c2_d_c * d_c_d_suuy;
	d_c2_d_suuz += d_c2_d_c * d_c_d_suuz;
	d_c2_d_suvx += d_c2_d_c * d_c_d_suvx;
	d_c2_d_suvy += d_c2_d_c * d_c_d_suvy;
	d_c2_d_suvz += d_c2_d_c * d_c_d_suvz;
	d_c2_d_rux += d_c2_d_b * d_b_d_rux;
	d_c2_d_ruy += d_c2_d_b * d_b_d_ruy;
	d_c2_d_ruz += d_c2_d_b * d_b_d_ruz;
	d_c2_d_rvx += d_c2_d_b * d_b_d_rvx;
	d_c2_d_rvy += d_c2_d_b * d_b_d_rvy;
	d_c2_d_rvz += d_c2_d_b * d_b_d_rvz;
	d_c2_d_suvx += d_c2_d_b * d_b_d_suvx;
	d_c2_d_suvy += d_c2_d_b * d_b_d_suvy;
	d_c2_d_suvz += d_c2_d_b * d_b_d_suvz;
	d_c2_d_svvx += d_c2_d_b * d_b_d_svvx;
	d_c2_d_svvy += d_c2_d_b * d_b_d_svvy;
	d_c2_d_svvz += d_c2_d_b * d_b_d_svvz;
	d_c2_d_rux += d_c2_d_a * d_a_d_rux;
	d_c2_d_ruy += d_c2_d_a * d_a_d_ruy;
	d_c2_d_ruz += d_c2_d_a * d_a_d_ruz;
	d_c2_d_rvx += d_c2_d_a * d_a_d_rvx;
	d_c2_d_rvy += d_c2_d_a * d_a_d_rvy;
	d_c2_d_rvz += d_c2_d_a * d_a_d_rvz;
	d_c2_d_suux += d_c2_d_a * d_a_d_suux;
	d_c2_d_suuy += d_c2_d_a * d_a_d_suuy;
	d_c2_d_suuz += d_c2_d_a * d_a_d_suuz;
	d_c2_d_suvx += d_c2_d_a * d_a_d_suvx;
	d_c2_d_suvy += d_c2_d_a * d_a_d_suvy;
	d_c2_d_suvz += d_c2_d_a * d_a_d_suvz;
	d_e_d_rux += d_e_d_c2 * d_c2_d_rux;
	d_e_d_ruy += d_e_d_c2 * d_c2_d_ruy;
	d_e_d_ruz += d_e_d_c2 * d_c2_d_ruz;
	d_e_d_rvx += d_e_d_c2 * d_c2_d_rvx;
	d_e_d_rvy += d_e_d_c2 * d_c2_d_rvy;
	d_e_d_rvz += d_e_d_c2 * d_c2_d_rvz;
	d_e_d_suux += d_e_d_c2 * d_c2_d_suux;
	d_e_d_suuy += d_e_d_c2 * d_c2_d_suuy;
	d_e_d_suuz += d_e_d_c2 * d_c2_d_suuz;
	d_e_d_suvx += d_e_d_c2 * d_c2_d_suvx;
	d_e_d_suvy += d_e_d_c2 * d_c2_d_suvy;
	d_e_d_suvz += d_e_d_c2 * d_c2_d_suvz;
	d_e_d_svvx += d_e_d_c2 * d_c2_d_svvx;
	d_e_d_svvy += d_e_d_c2 * d_c2_d_svvy;
	d_e_d_svvz += d_e_d_c2 * d_c2_d_svvz;
	d_e_d_rux += d_e_d_c1 * d_c1_d_rux;
	d_e_d_ruy += d_e_d_c1 * d_c1_d_ruy;
	d_e_d_ruz += d_e_d_c1 * d_c1_d_ruz;
	d_e_d_rvx += d_e_d_c1 * d_c1_d_rvx;
	d_e_d_rvy += d_e_d_c1 * d_c1_d_rvy;
	d_e_d_rvz += d_e_d_c1 * d_c1_d_rvz;
	d_e_d_suux += d_e_d_c1 * d_c1_d_suux;
	d_e_d_suuy += d_e_d_c1 * d_c1_d_suuy;
	d_e_d_suuz += d_e_d_c1 * d_c1_d_suuz;
	d_e_d_suvx += d_e_d_c1 * d_c1_d_suvx;
	d_e_d_suvy += d_e_d_c1 * d_c1_d_suvy;
	d_e_d_suvz += d_e_d_c1 * d_c1_d_suvz;
	d_e_d_svvx += d_e_d_c1 * d_c1_d_svvx;
	d_e_d_svvy += d_e_d_c1 * d_c1_d_svvy;
	d_e_d_svvz += d_e_d_c1 * d_c1_d_svvz;
	d_e_d_rux += d_e_d_g * d_g_d_rux;
	d_e_d_ruy += d_e_d_g * d_g_d_ruy;
	d_e_d_ruz += d_e_d_g * d_g_d_ruz;
	d_e_d_rvx += d_e_d_g * d_g_d_rvx;
	d_e_d_rvy += d_e_d_g * d_g_d_rvy;
	d_e_d_rvz += d_e_d_g * d_g_d_rvz;


#ifdef DEBUG_G
			for( int p = 0; p < np; p++ )
			{
				gr[3*cp[p]+0] += d_c2_d_rux * theIrregularFormulas[frm].r_u[p]; 
				gr[3*cp[p]+1] += d_c2_d_ruy * theIrregularFormulas[frm].r_u[p]; 
				gr[3*cp[p]+2] += d_c2_d_ruz * theIrregularFormulas[frm].r_u[p]; 
				
				gr[3*cp[p]+0] += d_c2_d_rvx * theIrregularFormulas[frm].r_v[p]; 
				gr[3*cp[p]+1] += d_c2_d_rvy * theIrregularFormulas[frm].r_v[p]; 
				gr[3*cp[p]+2] += d_c2_d_rvz * theIrregularFormulas[frm].r_v[p]; 
				
				gr[3*cp[p]+0] += d_c2_d_suux * theIrregularFormulas[frm].r_uu[p]; 
				gr[3*cp[p]+1] += d_c2_d_suuy * theIrregularFormulas[frm].r_uu[p]; 
				gr[3*cp[p]+2] += d_c2_d_suuz * theIrregularFormulas[frm].r_uu[p]; 
				
				gr[3*cp[p]+0] += d_c2_d_suvx * theIrregularFormulas[frm].r_uv[p]; 
				gr[3*cp[p]+1] += d_c2_d_suvy * theIrregularFormulas[frm].r_uv[p]; 
				gr[3*cp[p]+2] += d_c2_d_suvz * theIrregularFormulas[frm].r_uv[p]; 
				
				gr[3*cp[p]+0] += d_c2_d_svvx * theIrregularFormulas[frm].r_vv[p]; 
				gr[3*cp[p]+1] += d_c2_d_svvy * theIrregularFormulas[frm].r_vv[p]; 
				gr[3*cp[p]+2] += d_c2_d_svvz * theIrregularFormulas[frm].r_vv[p]; 
	
			}
#elif defined(DEBUG_NRMX)
			for( int p = 0; p < np; p++ )
			{
				gr[3*cp[p]+0] += d_nrmz_d_rux * theIrregularFormulas[frm].r_u[p]; 
				gr[3*cp[p]+1] += d_nrmz_d_ruy * theIrregularFormulas[frm].r_u[p]; 
				gr[3*cp[p]+2] += d_nrmz_d_ruz * theIrregularFormulas[frm].r_u[p]; 
				
				gr[3*cp[p]+0] += d_nrmz_d_rvx * theIrregularFormulas[frm].r_v[p]; 
				gr[3*cp[p]+1] += d_nrmz_d_rvy * theIrregularFormulas[frm].r_v[p]; 
				gr[3*cp[p]+2] += d_nrmz_d_rvz * theIrregularFormulas[frm].r_v[p]; 
				
//				gr[3*cp[p]+0] += d_nrmx_d_suux * theIrregularFormulas[frm].r_uu[p]; 
//				gr[3*cp[p]+1] += d_nrmx_d_suuy * theIrregularFormulas[frm].r_uu[p]; 
//				gr[3*cp[p]+2] += d_nrmx_d_suuz * theIrregularFormulas[frm].r_uu[p]; 
				
//				gr[3*cp[p]+0] += d_nsuu_d_suvx * theIrregularFormulas[frm].r_uv[p]; 
//				gr[3*cp[p]+1] += d_nsuu_d_suvy * theIrregularFormulas[frm].r_uv[p]; 
//				gr[3*cp[p]+2] += d_nsuu_d_suvz * theIrregularFormulas[frm].r_uv[p]; 
				
//				gr[3*cp[p]+0] += d_a_d_svvx * theIrregularFormulas[frm].r_vv[p]; 
//				gr[3*cp[p]+1] += d_a_d_svvy * theIrregularFormulas[frm].r_vv[p]; 
//				gr[3*cp[p]+2] += d_a_d_svvz * theIrregularFormulas[frm].r_vv[p]; 
			}

#else
			for( int p = 0; p < np; p++ )
			{
				gr[3*cp[p]+0] += d_e_d_rux * theIrregularFormulas[frm].r_u[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_ruy * theIrregularFormulas[frm].r_u[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_ruz * theIrregularFormulas[frm].r_u[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_rvx * theIrregularFormulas[frm].r_v[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_rvy * theIrregularFormulas[frm].r_v[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_rvz * theIrregularFormulas[frm].r_v[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_suux * theIrregularFormulas[frm].r_uu[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_suuy * theIrregularFormulas[frm].r_uu[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_suuz * theIrregularFormulas[frm].r_uu[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_suvx * theIrregularFormulas[frm].r_uv[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_suvy * theIrregularFormulas[frm].r_uv[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_suvz * theIrregularFormulas[frm].r_uv[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_svvx * theIrregularFormulas[frm].r_vv[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_svvy * theIrregularFormulas[frm].r_vv[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_svvz * theIrregularFormulas[frm].r_vv[p] * alpha_z; 
				
				dedr[3*p+0] += d_e_d_rux * theIrregularFormulas[frm].r_u[p]; 
				dedr[3*p+1] += d_e_d_ruy * theIrregularFormulas[frm].r_u[p]; 
				dedr[3*p+2] += d_e_d_ruz * theIrregularFormulas[frm].r_u[p]; 
				
				dedr[3*p+0] += d_e_d_rvx * theIrregularFormulas[frm].r_v[p]; 
				dedr[3*p+1] += d_e_d_rvy * theIrregularFormulas[frm].r_v[p]; 
				dedr[3*p+2] += d_e_d_rvz * theIrregularFormulas[frm].r_v[p]; 
				
				dedr[3*p+0] += d_e_d_suux * theIrregularFormulas[frm].r_uu[p]; 
				dedr[3*p+1] += d_e_d_suuy * theIrregularFormulas[frm].r_uu[p]; 
				dedr[3*p+2] += d_e_d_suuz * theIrregularFormulas[frm].r_uu[p]; 
				
				dedr[3*p+0] += d_e_d_suvx * theIrregularFormulas[frm].r_uv[p]; 
				dedr[3*p+1] += d_e_d_suvy * theIrregularFormulas[frm].r_uv[p]; 
				dedr[3*p+2] += d_e_d_suvz * theIrregularFormulas[frm].r_uv[p]; 
				
				dedr[3*p+0] += d_e_d_svvx * theIrregularFormulas[frm].r_vv[p]; 
				dedr[3*p+1] += d_e_d_svvy * theIrregularFormulas[frm].r_vv[p]; 
				dedr[3*p+2] += d_e_d_svvz * theIrregularFormulas[frm].r_vv[p]; 


			}
#endif
			for( int p = 0; p < np; p++ )
			{
				gr[3*nv+0] += dedr[3*p+0] * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				gr[3*nv+1] += dedr[3*p+1] * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				gr[3*nv+2] += dedr[3*p+2] * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
			}

		}	

	} 

	r_val /= area;

}
*/

// Particle helfrich gradient routine.

void surface::particle_H_grad( double *r, double *gr, int f, double u, double v, double p_area, double p_c0, double *p_uv_g )
{
	double alpha_x = r[3*nv];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];

	// gradient of the particle energy wrt the vertices.

	if( f < nf_faces )
	{
		int frm = nf_g_q_p*f;
	
		int *cp = theFormulas[frm].cp;

		double w = 1-u-v;

		double u2 = u*u;
		double u3 = u*u*u;
		double u4 = u*u*u*u;
		
		double v2 = v*v;
		double v3 = v*v*v;
		double v4 = v*v*v*v;
		
		double w2 = w*w;
		double w3 = w*w*w;
		double w4 = w*w*w*w;

		double n1 = (1.0/12.0)*(u4+2*u3*v); 
		double n2 = (1.0/12.0)*(u4+2*u3*w); 
		double n3 = (1.0/12.0)*(u4+2*u3*w+6*u3*v+6*u2*v*w+12*u2*v2+6*u*v2*w+6*u*v3+2*v3*w+v4); 
		double n4 = (1.0/12.0)*(6*u4+24*u3*w+24*u2*w2+8*u*w3+w4+24*u3*v+60*u2*v*w+36*u*v*w2+6*v*w3+24*u2*v2+36*u*v2*w+12*v2*w2+8*u*v3+6*v3*w+v4); 
		double n5 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+2*u3*v+6*u2*v*w+6*u*v*w2+2*v*w3); 
		double n6 = (1.0/12.0)*(2*u*v3+v4); 
		double n7 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+8*u3*v+36*u2*v*w+36*u*v*w2+8*v*w3+24*u2*v2+60*u*v2*w+24*v2*w2+24*u*v3+24*v3*w+6*v4); 
		double n8 = (1.0/12.0)*(u4+8*u3*w+24*u2*w2+24*u*w3+6*w4+6*u3*v+36*u2*v*w+60*u*v*w2+24*v*w3+12*u2*v2+36*u*v2*w+24*v2*w2+6*u*v3+8*v3*w+v4); 
		double n9 = (1.0/12.0)*(2*u*w3+w4); 
		double n10 = (1.0/12.0)*(2*v3*w+v4); 
		double n11 = (1.0/12.0)*(2*u*w3+w4+6*u*v*w2+6*v*w3+6*u*v2*w+12*v2*w2+2*u*v3+6*v3*w+v4); 
		double n12 = (1.0/12.0)*(w4+2*v*w3);

		double du_1 = Power(u,3)/3. + (Power(u,2)*v)/2.;
		double du_2 = Power(u,2)/2. - Power(u,3)/3. - (Power(u,2)*v)/2.; 	
		double du_3 = Power(u,2)/2. - Power(u,3)/3. + u*v - (Power(u,2)*v)/2. + Power(v,2)/2. - Power(v,3)/6.;
		double du_4 = 0.3333333333333333 + u - Power(u,2) - Power(u,3)/3. + v/2. - u*v - (Power(u,2)*v)/2. - Power(v,2) + Power(v,3)/3.;	
		double du_5 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. - v/2. + Power(u,2)*v + Power(v,2)/2. - Power(v,3)/6.;
		double du_6 = Power(v,3)/6.;	
		double du_7 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. + v/2. - 2*u*v + Power(u,2)*v - Power(v,2)/2. - Power(v,3)/6.;	
		double du_8 = -2*u + 2*Power(u,2) - Power(u,3)/3. - v + 2*u*v - (Power(u,2)*v)/2. + Power(v,2) - Power(v,3)/6.;	
		double du_9 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. + v/2. - (Power(u,2)*v)/2. - Power(v,2)/2. + Power(v,3)/6.;	
		double du_10 = -Power(v,3)/6.;
		double du_11 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. - v/2. + u*v - (Power(u,2)*v)/2. + Power(v,3)/3.;	
		double du_12 = 	-0.3333333333333333 + u - Power(u,2) + Power(u,3)/3. + v/2. - u*v + (Power(u,2)*v)/2. - Power(v,3)/6.;

		double dv_1 = Power(u,3)/6.;
		double dv_2 = -Power(u,3)/6.;
		double dv_3 = Power(u,2)/2. - Power(u,3)/6. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_4 = 0.16666666666666666 + u/2. - Power(u,2)/2. - Power(u,3)/6. - 2*u*v - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
		double dv_5 = -0.16666666666666666 - u/2. + Power(u,3)/3. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_6 = (u*Power(v,2))/2. + Power(v,3)/3.;
		double dv_7 = 0.3333333333333333 + u/2. - Power(u,2) + Power(u,3)/3. + v - u*v - Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_8 = -u + Power(u,2) - Power(u,3)/6. - 2*v + 2*u*v + 2*Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_9 = -0.3333333333333333 + u/2. - Power(u,3)/6. + v - u*v - Power(v,2) + (u*Power(v,2))/2. + Power(v,3)/3.;
		double dv_10 = Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_11 = 0.16666666666666666 - u/2. + Power(u,2)/2. - Power(u,3)/6. - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
		double dv_12 = -0.16666666666666666 + u/2. - Power(u,2)/2. + Power(u,3)/6. + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;

		double d_uu_1 = Power(u,2) + u*v;
		double d_uu_2 = u - Power(u,2) - u*v;
		double d_uu_3 = u - Power(u,2) + v - u*v;
		double d_uu_4 = 1 - 2*u - Power(u,2) - v - u*v;
		double d_uu_5 = -2*u + 2*Power(u,2) + 2*u*v;
		double d_uu_6 = 0;
		double d_uu_7 = -2*u + 2*Power(u,2) - 2*v + 2*u*v;
		double d_uu_8 = -2 + 4*u - Power(u,2) + 2*v - u*v;
		double d_uu_9 = u - Power(u,2) - u*v;
		double d_uu_10 = 0;
		double d_uu_11 = u - Power(u,2) + v - u*v;
		double d_uu_12 = 1 - 2*u + Power(u,2) - v + u*v;
		
		double d_uv_1 = Power(u,2)/2.;
		double d_uv_2 = -Power(u,2)/2.;
		double d_uv_3 = u - Power(u,2)/2. + v - Power(v,2)/2.;
		double d_uv_4 = 0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
		double d_uv_5 = -0.5 + Power(u,2) + v - Power(v,2)/2.;
		double d_uv_6 = Power(v,2)/2.;
		double d_uv_7 = 0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
		double d_uv_8 = -1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
		double d_uv_9 = 0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
		double d_uv_10 = -Power(v,2)/2.;
		double d_uv_11 = -0.5 + u - Power(u,2)/2. + Power(v,2);
		double d_uv_12 = 0.5 - u + Power(u,2)/2. - Power(v,2)/2.;
		
		double d_vv_1 = 0;
		double d_vv_2 = 0;
		double d_vv_3 = u + v - u*v - Power(v,2);
		double d_vv_4 = -2*u - 2*v + 2*u*v + 2*Power(v,2);
		double d_vv_5 = u + v - u*v - Power(v,2);
		double d_vv_6 = u*v + Power(v,2);
		double d_vv_7 = 1 - u - 2*v - u*v - Power(v,2);
		double d_vv_8 = -2 + 2*u + 4*v - u*v - Power(v,2);
		double d_vv_9 = 1 - u - 2*v + u*v + Power(v,2);
		double d_vv_10 = v - u*v - Power(v,2);
		double d_vv_11 = -2*v + 2*u*v + 2*Power(v,2);
		double d_vv_12 = v - u*v - Power(v,2);
		

		double d_uuu_1 = 2*u + v;
		double d_uuu_2 = 1-2*u-v; //u - Power(u,2) - u*v;
		double d_uuu_3 = 1-2*u-v;//u - Power(u,2) + v - u*v;
		double d_uuu_4 = -2-2*u-v;//1 - 2*u - Power(u,2) - v - u*v;
		double d_uuu_5 = -2+4*u+2*v;//-2*u + 2*Power(u,2) + 2*u*v;
		double d_uuu_6 = 0;
		double d_uuu_7 = -2+4*u+2*v;//-2*u + 2*Power(u,2) - 2*v + 2*u*v;
		double d_uuu_8 = 4-2*u-v;//-2 + 4*u - Power(u,2) + 2*v - u*v;
		double d_uuu_9 = 1-2*u-v;//u - Power(u,2) - u*v;
		double d_uuu_10 = 0;
		double d_uuu_11 = 1-2*u-v;//u - Power(u,2) + v - u*v;
		double d_uuu_12 = -2+2*u+v;//1 - 2*u + Power(u,2) - v + u*v;

		double d_uuv_1 = u;//Power(u,2)/2.;
		double d_uuv_2 = -u;//-Power(u,2)/2.;
		double d_uuv_3 = 1-u;//u - Power(u,2)/2. + v - Power(v,2)/2.;
		double d_uuv_4 = -1-u;//0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
		double d_uuv_5 = 2*u;//-0.5 + Power(u,2) + v - Power(v,2)/2.;
		double d_uuv_6 = 0;//Power(v,2)/2.;
		double d_uuv_7 = -2+2*u;//0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
		double d_uuv_8 = 2-u;//-1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
		double d_uuv_9 = -u;//0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
		double d_uuv_10 = 0;//-Power(v,2)/2.;
		double d_uuv_11 = 1-u;//-0.5 + u - Power(u,2)/2. + Power(v,2);
		double d_uuv_12 = -1+u;//0.5 - u + Power(u,2)/2. - Power(v,2)/2.;

		double d_uvv_1 = 0;
		double d_uvv_2 = 0;
		double d_uvv_3 = 1-v;//u + v - u*v - Power(v,2);
		double d_uvv_4 = -2+2*v;//-2*u - 2*v + 2*u*v + 2*Power(v,2);
		double d_uvv_5 = 1-v;//u + v - u*v - Power(v,2);
		double d_uvv_6 = v;//u*v + Power(v,2);
		double d_uvv_7 = -1-v;//1 - u - 2*v - u*v - Power(v,2);
		double d_uvv_8 = 2-v;//-2 + 2*u + 4*v - u*v - Power(v,2);
		double d_uvv_9 = -1+v;//1 - u - 2*v + u*v + Power(v,2);
		double d_uvv_10 = -v;//v - u*v - Power(v,2);
		double d_uvv_11 = 2*v;//-2*v + 2*u*v + 2*Power(v,2);
		double d_uvv_12 = -v;//v - u*v - Power(v,2);
		
		double d_vvv_1 = 0;
		double d_vvv_2 = 0;
		double d_vvv_3 = 1-u-2*v;//u + v - u*v - Power(v,2);
		double d_vvv_4 = -2+2*u+4*v;//-2*u - 2*v + 2*u*v + 2*Power(v,2);
		double d_vvv_5 = 1-u-2*v;//u + v - u*v - Power(v,2);
		double d_vvv_6 = u+2*v;//u*v + Power(v,2);
		double d_vvv_7 = -2-u-2*v;//1 - u - 2*v - u*v - Power(v,2);
		double d_vvv_8 = 4-u-2*v;//-2 + 2*u + 4*v - u*v - Power(v,2);
		double d_vvv_9 = -2+u+2*v;//1 - u - 2*v + u*v + Power(v,2);
		double d_vvv_10 = 1-u-2*v;//v - u*v - Power(v,2);
		double d_vvv_11 = -2+2*u+4*v;//-2*v + 2*u*v + 2*Power(v,2);
		double d_vvv_12 = 1-u-2*v;//v - u*v - Power(v,2);

		double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
		double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
		double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
		
		double ceff_map_duu[12] = { d_uu_8, d_uu_7, d_uu_4, d_uu_5, d_uu_9, d_uu_12, d_uu_11, d_uu_10, d_uu_6, d_uu_3, d_uu_1, d_uu_2 };
		double ceff_map_duv[12] = { d_uv_8, d_uv_7, d_uv_4, d_uv_5, d_uv_9, d_uv_12, d_uv_11, d_uv_10, d_uv_6, d_uv_3, d_uv_1, d_uv_2 };
		double ceff_map_dvv[12] = { d_vv_8, d_vv_7, d_vv_4, d_vv_5, d_vv_9, d_vv_12, d_vv_11, d_vv_10, d_vv_6, d_vv_3, d_vv_1, d_vv_2 };
		
		double ceff_map_duuu[12] = { d_uuu_8, d_uuu_7, d_uuu_4, d_uuu_5, d_uuu_9, d_uuu_12, d_uuu_11, d_uuu_10, d_uuu_6, d_uuu_3, d_uuu_1, d_uuu_2 };
		double ceff_map_duuv[12] = { d_uuv_8, d_uuv_7, d_uuv_4, d_uuv_5, d_uuv_9, d_uuv_12, d_uuv_11, d_uuv_10, d_uuv_6, d_uuv_3, d_uuv_1, d_uuv_2 };
		double ceff_map_duvv[12] = { d_uvv_8, d_uvv_7, d_uvv_4, d_uvv_5, d_uvv_9, d_uvv_12, d_uvv_11, d_uvv_10, d_uvv_6, d_uvv_3, d_uvv_1, d_uvv_2 };
		double ceff_map_dvvv[12] = { d_vvv_8, d_vvv_7, d_vvv_4, d_vvv_5, d_vvv_9, d_vvv_12, d_vvv_11, d_vvv_10, d_vvv_6, d_vvv_3, d_vvv_1, d_vvv_2 };
				
		double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};

		double tSuuu[3] = {0,0,0};
		double tSuuv[3] = {0,0,0};
		double tSuvv[3] = {0,0,0};
		double tSvvv[3] = {0,0,0};

		double nrm[3]={0,0,0}; 

		int np = theFormulas[frm].ncoor;
		double dedr[3*np];
		memset( dedr, 0, sizeof(double) * 3 *np );

		for( int p = 0; p < np; p++ )
		{

			R[0] += ceff_map[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			R[1] += ceff_map[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			R[2] += ceff_map[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			Ru[0] += ceff_map_du[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			Ru[1] += ceff_map_du[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			Ru[2] += ceff_map_du[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			Rv[0] += ceff_map_dv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			Rv[1] += ceff_map_dv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			Rv[2] += ceff_map_dv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSuu[0] += ceff_map_duu[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSuu[1] += ceff_map_duu[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSuu[2] += ceff_map_duu[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSuv[0] += ceff_map_duv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSuv[1] += ceff_map_duv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSuv[2] += ceff_map_duv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSvv[0] += ceff_map_dvv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSvv[1] += ceff_map_dvv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSvv[2] += ceff_map_dvv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSuuu[0] += ceff_map_duuu[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSuuu[1] += ceff_map_duuu[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSuuu[2] += ceff_map_duuu[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSuuv[0] += ceff_map_duuv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSuuv[1] += ceff_map_duuv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSuuv[2] += ceff_map_duuv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSuvv[0] += ceff_map_duvv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSuvv[1] += ceff_map_duvv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSuvv[2] += ceff_map_duvv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSvvv[0] += ceff_map_dvvv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSvvv[1] += ceff_map_dvvv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSvvv[2] += ceff_map_dvvv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
		}

		cross( Ru, Rv, nrm );
		normalize(nrm);


		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		double g = sqrt(RuRu*RvRv-RuRv*RuRv);

		double e1,e2;
      
		double nsuu = tSuu[0] * nrm[0] + tSuu[1] * nrm[1] + tSuu[2] * nrm[2];
		double nsuv = tSuv[0] * nrm[0] + tSuv[1] * nrm[1] + tSuv[2] * nrm[2];
		double nsvv = tSvv[0] * nrm[0] + tSvv[1] * nrm[1] + tSvv[2] * nrm[2];

		double Stot = (nsuu * RvRv + nsvv * RuRu -2*nsuv*RuRv)/(g*g);

		double Sop[4] = { 1.0/(g*g) * (nsuu * RvRv  - nsuv * RuRv), (1.0/(g*g)) * ( nsuv*RvRv-nsvv*RuRv),
				  1.0/(g*g) * (nsuv * RuRu  - nsuu * RuRv), (1.0/(g*g)) * ( nsvv*RuRu-nsuv*RuRv) };

		double a = Sop[0];		
		double b = Sop[1];
		double c = Sop[2];
		double d = Sop[3];

		double c0 = theFormulas[frm].c0;
		e1 = -0.5*(a+d-sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
		e2 = -0.5*(a+d+sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
	

double d_nrmz_d_rux=0,d_nrmz_d_ruy=0,d_nrmz_d_ruz=0,d_nrmz_d_rvx=0,d_nrmz_d_rvy=0,d_nrmz_d_rvz=0,d_nrmy_d_rux=0,d_nrmy_d_ruy=0,d_nrmy_d_ruz=0,d_nrmy_d_rvx=0,d_nrmy_d_rvy=0,d_nrmy_d_rvz=0,d_nrmx_d_rux=0,d_nrmx_d_ruy=0,d_nrmx_d_ruz=0,d_nrmx_d_rvx=0,d_nrmx_d_rvy=0,d_nrmx_d_rvz=0,d_RuRu_d_rux=0,d_RuRu_d_ruy=0,d_RuRu_d_ruz=0,d_nsvv_d_rux=0,d_nsvv_d_ruy=0,d_nsvv_d_ruz=0,d_nsvv_d_rvx=0,d_nsvv_d_rvy=0,d_nsvv_d_rvz=0,d_nsvv_d_svvx=0,d_nsvv_d_svvy=0,d_nsvv_d_svvz=0,d_RuRv_d_rux=0,d_RuRv_d_ruy=0,d_RuRv_d_ruz=0,d_RuRv_d_rvx=0,d_RuRv_d_rvy=0,d_RuRv_d_rvz=0,d_nsuv_d_rux=0,d_nsuv_d_ruy=0,d_nsuv_d_ruz=0,d_nsuv_d_rvx=0,d_nsuv_d_rvy=0,d_nsuv_d_rvz=0,d_nsuv_d_suvx=0,d_nsuv_d_suvy=0,d_nsuv_d_suvz=0,d_RvRv_d_rvx=0,d_RvRv_d_rvy=0,d_RvRv_d_rvz=0,d_nsuu_d_rux=0,d_nsuu_d_ruy=0,d_nsuu_d_ruz=0,d_nsuu_d_rvx=0,d_nsuu_d_rvy=0,d_nsuu_d_rvz=0,d_nsuu_d_suux=0,d_nsuu_d_suuy=0,d_nsuu_d_suuz=0,d_g_d_rux=0,d_g_d_ruy=0,d_g_d_ruz=0,d_g_d_rvx=0,d_g_d_rvy=0,d_g_d_rvz=0,d_e_d_rux=0,d_e_d_ruy=0,d_e_d_ruz=0,d_e_d_rvx=0,d_e_d_rvy=0,d_e_d_rvz=0,d_e_d_suux=0,d_e_d_suuy=0,d_e_d_suuz=0,d_e_d_suvx=0,d_e_d_suvy=0,d_e_d_suvz=0,d_e_d_svvx=0,d_e_d_svvy=0,d_e_d_svvz=0,d_c2_d_rux=0,d_c2_d_ruy=0,d_c2_d_ruz=0,d_c2_d_rvx=0,d_c2_d_rvy=0,d_c2_d_rvz=0,d_c2_d_suux=0,d_c2_d_suuy=0,d_c2_d_suuz=0,d_c2_d_suvx=0,d_c2_d_suvy=0,d_c2_d_suvz=0,d_c2_d_svvx=0,d_c2_d_svvy=0,d_c2_d_svvz=0,d_d_d_rux=0,d_d_d_ruy=0,d_d_d_ruz=0,d_d_d_rvx=0,d_d_d_rvy=0,d_d_d_rvz=0,d_d_d_suvx=0,d_d_d_suvy=0,d_d_d_suvz=0,d_d_d_svvx=0,d_d_d_svvy=0,d_d_d_svvz=0,d_c_d_rux=0,d_c_d_ruy=0,d_c_d_ruz=0,d_c_d_rvx=0,d_c_d_rvy=0,d_c_d_rvz=0,d_c_d_suux=0,d_c_d_suuy=0,d_c_d_suuz=0,d_c_d_suvx=0,d_c_d_suvy=0,d_c_d_suvz=0,d_b_d_rux=0,d_b_d_ruy=0,d_b_d_ruz=0,d_b_d_rvx=0,d_b_d_rvy=0,d_b_d_rvz=0,d_b_d_suvx=0,d_b_d_suvy=0,d_b_d_suvz=0,d_b_d_svvx=0,d_b_d_svvy=0,d_b_d_svvz=0,d_a_d_rux=0,d_a_d_ruy=0,d_a_d_ruz=0,d_a_d_rvx=0,d_a_d_rvy=0,d_a_d_rvz=0,d_a_d_suux=0,d_a_d_suuy=0,d_a_d_suuz=0,d_a_d_suvx=0,d_a_d_suvy=0,d_a_d_suvz=0,d_c1_d_rux=0,d_c1_d_ruy=0,d_c1_d_ruz=0,d_c1_d_rvx=0,d_c1_d_rvy=0,d_c1_d_rvz=0,d_c1_d_suux=0,d_c1_d_suuy=0,d_c1_d_suuz=0,d_c1_d_suvx=0,d_c1_d_suvy=0,d_c1_d_suvz=0,d_c1_d_svvx=0,d_c1_d_svvy=0,d_c1_d_svvz=0,d_g_d_RvRv=0,d_g_d_RuRv=0,d_g_d_RuRu=0,d_nsvv_d_nrmz=0,d_nsvv_d_nrmy=0,d_nsvv_d_nrmx=0,d_nsuv_d_nrmz=0,d_nsuv_d_nrmy=0,d_nsuv_d_nrmx=0,d_nsuu_d_nrmz=0,d_nsuu_d_nrmy=0,d_nsuu_d_nrmx=0,d_d_d_RuRv=0,d_d_d_nsuv=0,d_d_d_RuRu=0,d_d_d_nsvv=0,d_d_d_g=0,d_c_d_RuRv=0,d_c_d_nsuv=0,d_c_d_RuRu=0,d_c_d_nsuu=0,d_c_d_g=0,d_b_d_RuRv=0,d_b_d_nsuv=0,d_b_d_RvRv=0,d_b_d_nsvv=0,d_b_d_g=0,d_a_d_RuRv=0,d_a_d_nsuv=0,d_a_d_RvRv=0,d_a_d_nsuu=0,d_a_d_g=0,d_e_d_c2=0,d_e_d_c1=0,d_e_d_g=0,d_c2_d_d=0,d_c2_d_c=0,d_c2_d_b=0,d_c2_d_a=0,d_c1_d_d=0,d_c1_d_c=0,d_c1_d_b=0,d_c1_d_a=0,junk;
		// the basic variables.





		d_g_d_RuRu = RvRv/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
		d_g_d_RuRv = -(RuRv/Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
		d_g_d_RvRv = RuRu/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  

		d_RuRu_d_rux = 2*Ru[0];
		d_RuRu_d_ruy = 2*Ru[1];
		d_RuRu_d_ruz = 2*Ru[2];
		
		d_RuRv_d_rux = Rv[0];
		d_RuRv_d_ruy = Rv[1];
		d_RuRv_d_ruz = Rv[2];
		
		d_RuRv_d_rvx = Ru[0];
		d_RuRv_d_rvy = Ru[1];
		d_RuRv_d_rvz = Ru[2];
		
		d_RvRv_d_rvx = 2*Rv[0];
		d_RvRv_d_rvy = 2*Rv[1];
		d_RvRv_d_rvz = 2*Rv[2];

		// intermediates.
		d_a_d_g = (-2*(-(nsuv*RuRv) + nsuu*RvRv))/Power(g,3);
		d_a_d_nsuu = RvRv/Power(g,2);
		d_a_d_RvRv = nsuu/Power(g,2);
		d_a_d_nsuv = -(RuRv/Power(g,2));
		d_a_d_RuRv = -(nsuv/Power(g,2));

		d_b_d_g = (-2*(-(nsvv*RuRv) + nsuv*RvRv))/Power(g,3);
		d_b_d_nsvv = -(RuRv/Power(g,2));
		d_b_d_RvRv = nsuv/Power(g,2);
		d_b_d_nsuv = RvRv/Power(g,2);
		d_b_d_RuRv = -(nsvv/Power(g,2));

		d_c_d_g = (-2*(nsuv*RuRu - nsuu*RuRv))/Power(g,3);
		d_c_d_nsuu = -(RuRv/Power(g,2));
		d_c_d_RuRu = nsuv/Power(g,2);
		d_c_d_nsuv = RuRu/Power(g,2);
		d_c_d_RuRv = -(nsuu/Power(g,2));

		d_d_d_g = (-2*(nsvv*RuRu - nsuv*RuRv))/Power(g,3);
		d_d_d_nsvv = RuRu/Power(g,2);
		d_d_d_RuRu = nsvv/Power(g,2);
		d_d_d_nsuv = -(RuRv/Power(g,2));
		d_d_d_RuRv = -(nsuv/Power(g,2));

		d_c1_d_a =-0.5*(1 - (2*a - 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
		d_c1_d_b =-(-1.*c)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2)) ;
		d_c1_d_c =-(-1.*b)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2)) ;
		d_c1_d_d =-0.5*(1 - (-2*a + 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
	
		d_c2_d_a = -0.5*(1 + (2*a - 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
		d_c2_d_b = -(1.*c)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2));
		d_c2_d_c = -(1.*b)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2));
		d_c2_d_d = -0.5*(1 + (-2*a + 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));


		// this is the term we need to gradient:
		// en += 0.5*kc * p_area * ( c1 + c2 - p_c0 ) * ( c1 + c2 - p_c0 );
		d_e_d_c1 = kc * p_area * ( e1 + e2 - p_c0);
		d_e_d_c2 = kc * p_area * ( e1 + e2 - p_c0);
		d_e_d_g = 0;

		// that's it.

		double fac = Power(Power(Ru[1]*Rv[0] - Ru[0]*Rv[1],2) + Power(Ru[2]*Rv[0] - Ru[0]*Rv[2],2) + Power(Ru[2]*Rv[1] - Ru[1]*Rv[2],2),1.5);

		d_nrmx_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac; 
		d_nrmx_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2)))/fac;
		d_nrmx_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac;
		
		d_nrmy_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2))))/fac;
		d_nrmy_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac; 
		d_nrmy_d_ruz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2])) + Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac;

		d_nrmz_d_rux = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2])/fac;
		d_nrmz_d_ruy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;
		d_nrmz_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2))) + (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;

		d_nrmx_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac; 
		d_nrmx_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2])))/fac;
		d_nrmx_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac;
		
		d_nrmy_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
		d_nrmy_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2])))/fac;
		d_nrmy_d_rvz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[0]*Ru[1]*Rv[0] - Power(Ru[0],2)*Rv[1] + Ru[2]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2])))/fac;

		d_nrmz_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;
		d_nrmz_d_rvy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[0]*Ru[2]*Rv[0] - Power(Ru[0],2)*Rv[2] + Ru[1]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
		d_nrmz_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;

		double d_nrmx_d_u = d_nrmx_d_rux * tSuu[0] + d_nrmx_d_ruy * tSuu[1] + d_nrmx_d_ruz * tSuu[2] +
				    d_nrmx_d_rvx * tSuv[0] + d_nrmx_d_rvy * tSuv[1] + d_nrmx_d_rvz * tSuv[2];
		double d_nrmx_d_v = d_nrmx_d_rux * tSuv[0] + d_nrmx_d_ruy * tSuv[1] + d_nrmx_d_ruz * tSuv[2] +
				    d_nrmx_d_rvx * tSvv[0] + d_nrmx_d_rvy * tSvv[1] + d_nrmx_d_rvz * tSvv[2];
		double d_nrmy_d_u = d_nrmy_d_rux * tSuu[0] + d_nrmy_d_ruy * tSuu[1] + d_nrmy_d_ruz * tSuu[2] +
				    d_nrmy_d_rvx * tSuv[0] + d_nrmy_d_rvy * tSuv[1] + d_nrmy_d_rvz * tSuv[2];
		double d_nrmy_d_v = d_nrmy_d_rux * tSuv[0] + d_nrmy_d_ruy * tSuv[1] + d_nrmy_d_ruz * tSuv[2] +
				    d_nrmy_d_rvx * tSvv[0] + d_nrmy_d_rvy * tSvv[1] + d_nrmy_d_rvz * tSvv[2];
		double d_nrmz_d_u = d_nrmz_d_rux * tSuu[0] + d_nrmz_d_ruy * tSuu[1] + d_nrmz_d_ruz * tSuu[2] +
				    d_nrmz_d_rvx * tSuv[0] + d_nrmz_d_rvy * tSuv[1] + d_nrmz_d_rvz * tSuv[2];
		double d_nrmz_d_v = d_nrmz_d_rux * tSuv[0] + d_nrmz_d_ruy * tSuv[1] + d_nrmz_d_ruz * tSuv[2] +
				    d_nrmz_d_rvx * tSvv[0] + d_nrmz_d_rvy * tSvv[1] + d_nrmz_d_rvz * tSvv[2];

		double d_nsuu_d_u = tSuuu[0] * nrm[0] + tSuuu[1] * nrm[1] + tSuuu[2] * nrm[2];
		       d_nsuu_d_u += tSuu[0] * d_nrmx_d_u + tSuu[1] * d_nrmy_d_u + tSuu[2] * d_nrmz_d_u;
		double d_nsuu_d_v = tSuuv[0] * nrm[0] + tSuuv[1] * nrm[1] + tSuuv[2] * nrm[2];
		       d_nsuu_d_v += tSuu[0] * d_nrmx_d_v + tSuu[1] * d_nrmy_d_v + tSuu[2] * d_nrmz_d_v;
		
		double d_nsuv_d_u = tSuuv[0] * nrm[0] + tSuuv[1] * nrm[1] + tSuuv[2] * nrm[2];
		       d_nsuv_d_u += tSuv[0] * d_nrmx_d_u + tSuv[1] * d_nrmy_d_u + tSuv[2] * d_nrmz_d_u;
		double d_nsuv_d_v = tSuvv[0] * nrm[0] + tSuvv[1] * nrm[1] + tSuvv[2] * nrm[2];
		       d_nsuv_d_v += tSuv[0] * d_nrmx_d_v + tSuv[1] * d_nrmy_d_v + tSuv[2] * d_nrmz_d_v;

		double d_nsvv_d_u = tSuvv[0] * nrm[0] + tSuvv[1] * nrm[1] + tSuvv[2] * nrm[2];
		       d_nsvv_d_u += tSvv[0] * d_nrmx_d_u + tSvv[1] * d_nrmy_d_u + tSvv[2] * d_nrmz_d_u;
		double d_nsvv_d_v = tSvvv[0] * nrm[0] + tSvvv[1] * nrm[1] + tSvvv[2] * nrm[2];
		       d_nsvv_d_v += tSvv[0] * d_nrmx_d_v + tSvv[1] * d_nrmy_d_v + tSvv[2] * d_nrmz_d_v;

		double d_RuRu_d_u = 2*Ru[0] * tSuu[0] + 2*Ru[1] * tSuu[1] + 2 * Ru[2] * tSuu[2];
		double d_RuRu_d_v = 2*Ru[0] * tSuv[0] + 2*Ru[1] * tSuv[1] + 2 * Ru[2] * tSuv[2];

		double d_RuRv_d_u =  Ru[0] * tSuv[0] + Ru[1] * tSuv[1] + Ru[2] * tSuv[2];
		       d_RuRv_d_u += Rv[0] * tSuu[0] + Rv[1] * tSuu[1] + Rv[2] * tSuu[2];
		
		double d_RuRv_d_v =  Ru[0] * tSvv[0] + Ru[1] * tSvv[1] + Ru[2] * tSvv[2];
		       d_RuRv_d_v += Rv[0] * tSuv[0] + Rv[1] * tSuv[1] + Rv[2] * tSuv[2];

		double d_RvRv_d_u = 2*Rv[0] * tSuv[0] + 2*Rv[1] * tSuv[1] + 2 * Rv[2] * tSuv[2];
		double d_RvRv_d_v = 2*Rv[0] * tSvv[0] + 2*Rv[1] * tSvv[1] + 2 * Rv[2] * tSvv[2];

		double d_g_d_u = 0;
		double d_g_d_v = 0;

		d_g_d_u += d_g_d_RuRu * 2* Ru[0] * tSuu[0] + d_g_d_RuRu * 2* Ru[1] * tSuu[1] + d_g_d_RuRu * 2* Ru[2] * tSuu[2];
		d_g_d_v += d_g_d_RuRu * 2* Ru[0] * tSuv[0] + d_g_d_RuRu * 2* Ru[1] * tSuv[1] + d_g_d_RuRu * 2* Ru[2] * tSuv[2];

		// d_ru dru_du
		d_g_d_u += d_g_d_RuRv * Rv[0] * tSuu[0] + d_g_d_RuRv * Rv[1] * tSuu[1] + d_g_d_RuRv * Rv[2] * tSuu[2];
		// d_rv_drv_du
		d_g_d_u += d_g_d_RuRv * Ru[0] * tSuv[0] + d_g_d_RuRv * Ru[1] * tSuv[1] + d_g_d_RuRv * Ru[2] * tSuv[2];	
		// d_ru dru_dv
		d_g_d_v += d_g_d_RuRv * Rv[0] * tSuv[0] + d_g_d_RuRv * Rv[1] * tSuv[1] + d_g_d_RuRv * Rv[2] * tSuv[2];
		// d_rv drv_dv
		d_g_d_v += d_g_d_RuRv * Ru[0] * tSvv[0] + d_g_d_RuRv * Ru[1] * tSvv[1] + d_g_d_RuRv * Ru[2] * tSvv[2];
		
		d_g_d_u += d_g_d_RvRv * 2* Rv[0] * tSuv[0] + d_g_d_RvRv * 2* Rv[1] * tSuv[1] + d_g_d_RvRv * 2* Rv[2] * tSuv[2];
		d_g_d_v += d_g_d_RvRv * 2* Rv[0] * tSvv[0] + d_g_d_RvRv * 2* Rv[1] * tSvv[1] + d_g_d_RvRv * 2* Rv[2] * tSvv[2];

//			d_a_d_g = (-2*(-(nsuv*RuRv) + nsuu*RvRv))/Power(g,3);
//			d_a_d_nsuu = RvRv/Power(g,2);
//			d_a_d_RvRv = nsuu/Power(g,2);
//			d_a_d_nsuv = -(RuRv/Power(g,2));
//			d_a_d_RuRv = -(nsuv/Power(g,2));
		
		double d_a_d_u = d_a_d_nsuu * d_nsuu_d_u + d_a_d_nsuv * d_nsuv_d_u + d_a_d_RuRv * d_RuRv_d_u + d_a_d_RvRv * d_RvRv_d_u + d_a_d_g * d_g_d_u;
		double d_a_d_v = d_a_d_nsuu * d_nsuu_d_v + d_a_d_nsuv * d_nsuv_d_v + d_a_d_RuRv * d_RuRv_d_v + d_a_d_RvRv * d_RvRv_d_v + d_a_d_g * d_g_d_v;
		double d_b_d_u = d_b_d_nsvv * d_nsvv_d_u + d_b_d_nsuv * d_nsuv_d_u + d_b_d_RuRv * d_RuRv_d_u + d_b_d_RvRv * d_RvRv_d_u + d_b_d_g * d_g_d_u;
		double d_b_d_v = d_b_d_nsvv * d_nsvv_d_v + d_b_d_nsuv * d_nsuv_d_v + d_b_d_RuRv * d_RuRv_d_v + d_b_d_RvRv * d_RvRv_d_v + d_b_d_g * d_g_d_v;
		double d_c_d_u = d_c_d_nsuu * d_nsuu_d_u + d_c_d_nsuv * d_nsuv_d_u + d_c_d_RuRv * d_RuRv_d_u + d_c_d_RuRu * d_RuRu_d_u + d_c_d_g * d_g_d_u;
		double d_c_d_v = d_c_d_nsuu * d_nsuu_d_v + d_c_d_nsuv * d_nsuv_d_v + d_c_d_RuRv * d_RuRv_d_v + d_c_d_RuRu * d_RuRu_d_v + d_c_d_g * d_g_d_v;
		double d_d_d_u = d_d_d_nsvv * d_nsvv_d_u + d_d_d_nsuv * d_nsuv_d_u + d_d_d_RuRv * d_RuRv_d_u + d_d_d_RuRu * d_RuRu_d_u + d_d_d_g * d_g_d_u;
		double d_d_d_v = d_d_d_nsvv * d_nsvv_d_v + d_d_d_nsuv * d_nsuv_d_v + d_d_d_RuRv * d_RuRv_d_v + d_d_d_RuRu * d_RuRu_d_v + d_d_d_g * d_g_d_v;

		double d_c1_d_u = d_c1_d_a * d_a_d_u + d_c1_d_b * d_b_d_u + d_c1_d_c * d_c_d_u + d_c1_d_d * d_d_d_u;
		double d_c2_d_u = d_c2_d_a * d_a_d_u + d_c2_d_b * d_b_d_u + d_c2_d_c * d_c_d_u + d_c2_d_d * d_d_d_u;

		double d_c1_d_v = d_c1_d_a * d_a_d_v + d_c1_d_b * d_b_d_v + d_c1_d_c * d_c_d_v + d_c1_d_d * d_d_d_v;
		double d_c2_d_v = d_c2_d_a * d_a_d_v + d_c2_d_b * d_b_d_v + d_c2_d_c * d_c_d_v + d_c2_d_d * d_d_d_v;

		double d_e_d_u = d_e_d_c1 * d_c1_d_u + d_e_d_c2 * d_c2_d_u;
		double d_e_d_v = d_e_d_c1 * d_c1_d_v + d_e_d_c2 * d_c2_d_v;

		p_uv_g[0] += d_e_d_u;
		p_uv_g[1] += d_e_d_v;

//			printf("u %le v %le a: %.14le da_du: %le nrmx: %.14le dnrmx_du %.14le\n", u, v, b, d_b_d_u, nrm[0], d_nrmy_d_u );

		double e = 0.5* kc * p_area * ( e1 + e2 - p_c0 ) * ( e1 + e2 - p_c0 );
//			printf("e %le at uv: %le %le  pg: %le %le\n", e, p_uv[2*pid+0], p_uv[2*pid+1], d_e_d_u, d_e_d_v );
		
		d_nsuu_d_nrmx = tSuu[0];
		d_nsuu_d_nrmy = tSuu[1];
		d_nsuu_d_nrmz = tSuu[2];
		
		d_nsuv_d_nrmx = tSuv[0];
		d_nsuv_d_nrmy = tSuv[1];
		d_nsuv_d_nrmz = tSuv[2];

		d_nsvv_d_nrmx = tSvv[0];
		d_nsvv_d_nrmy = tSvv[1];
		d_nsvv_d_nrmz = tSvv[2];

		d_nsuu_d_suux = nrm[0];
		d_nsuu_d_suuy = nrm[1];
		d_nsuu_d_suuz = nrm[2];
		
		d_nsuv_d_suvx = nrm[0];
		d_nsuv_d_suvy = nrm[1];
		d_nsuv_d_suvz = nrm[2];

		d_nsvv_d_svvx = nrm[0];
		d_nsvv_d_svvy = nrm[1];
		d_nsvv_d_svvz = nrm[2];



d_nsvv_d_rux += d_nsvv_d_nrmz * d_nrmz_d_rux;
d_nsvv_d_ruy += d_nsvv_d_nrmz * d_nrmz_d_ruy;
d_nsvv_d_ruz += d_nsvv_d_nrmz * d_nrmz_d_ruz;
d_nsvv_d_rvx += d_nsvv_d_nrmz * d_nrmz_d_rvx;
d_nsvv_d_rvy += d_nsvv_d_nrmz * d_nrmz_d_rvy;
d_nsvv_d_rvz += d_nsvv_d_nrmz * d_nrmz_d_rvz;
d_nsvv_d_rux += d_nsvv_d_nrmy * d_nrmy_d_rux;
d_nsvv_d_ruy += d_nsvv_d_nrmy * d_nrmy_d_ruy;
d_nsvv_d_ruz += d_nsvv_d_nrmy * d_nrmy_d_ruz;
d_nsvv_d_rvx += d_nsvv_d_nrmy * d_nrmy_d_rvx;
d_nsvv_d_rvy += d_nsvv_d_nrmy * d_nrmy_d_rvy;
d_nsvv_d_rvz += d_nsvv_d_nrmy * d_nrmy_d_rvz;
d_nsvv_d_rux += d_nsvv_d_nrmx * d_nrmx_d_rux;
d_nsvv_d_ruy += d_nsvv_d_nrmx * d_nrmx_d_ruy;
d_nsvv_d_ruz += d_nsvv_d_nrmx * d_nrmx_d_ruz;
d_nsvv_d_rvx += d_nsvv_d_nrmx * d_nrmx_d_rvx;
d_nsvv_d_rvy += d_nsvv_d_nrmx * d_nrmx_d_rvy;
d_nsvv_d_rvz += d_nsvv_d_nrmx * d_nrmx_d_rvz;
d_nsuv_d_rux += d_nsuv_d_nrmz * d_nrmz_d_rux;
d_nsuv_d_ruy += d_nsuv_d_nrmz * d_nrmz_d_ruy;
d_nsuv_d_ruz += d_nsuv_d_nrmz * d_nrmz_d_ruz;
d_nsuv_d_rvx += d_nsuv_d_nrmz * d_nrmz_d_rvx;
d_nsuv_d_rvy += d_nsuv_d_nrmz * d_nrmz_d_rvy;
d_nsuv_d_rvz += d_nsuv_d_nrmz * d_nrmz_d_rvz;
d_nsuv_d_rux += d_nsuv_d_nrmy * d_nrmy_d_rux;
d_nsuv_d_ruy += d_nsuv_d_nrmy * d_nrmy_d_ruy;
d_nsuv_d_ruz += d_nsuv_d_nrmy * d_nrmy_d_ruz;
d_nsuv_d_rvx += d_nsuv_d_nrmy * d_nrmy_d_rvx;
d_nsuv_d_rvy += d_nsuv_d_nrmy * d_nrmy_d_rvy;
d_nsuv_d_rvz += d_nsuv_d_nrmy * d_nrmy_d_rvz;
d_nsuv_d_rux += d_nsuv_d_nrmx * d_nrmx_d_rux;
d_nsuv_d_ruy += d_nsuv_d_nrmx * d_nrmx_d_ruy;
d_nsuv_d_ruz += d_nsuv_d_nrmx * d_nrmx_d_ruz;
d_nsuv_d_rvx += d_nsuv_d_nrmx * d_nrmx_d_rvx;
d_nsuv_d_rvy += d_nsuv_d_nrmx * d_nrmx_d_rvy;
d_nsuv_d_rvz += d_nsuv_d_nrmx * d_nrmx_d_rvz;
d_nsuu_d_rux += d_nsuu_d_nrmz * d_nrmz_d_rux;
d_nsuu_d_ruy += d_nsuu_d_nrmz * d_nrmz_d_ruy;
d_nsuu_d_ruz += d_nsuu_d_nrmz * d_nrmz_d_ruz;
d_nsuu_d_rvx += d_nsuu_d_nrmz * d_nrmz_d_rvx;
d_nsuu_d_rvy += d_nsuu_d_nrmz * d_nrmz_d_rvy;
d_nsuu_d_rvz += d_nsuu_d_nrmz * d_nrmz_d_rvz;
d_nsuu_d_rux += d_nsuu_d_nrmy * d_nrmy_d_rux;
d_nsuu_d_ruy += d_nsuu_d_nrmy * d_nrmy_d_ruy;
d_nsuu_d_ruz += d_nsuu_d_nrmy * d_nrmy_d_ruz;
d_nsuu_d_rvx += d_nsuu_d_nrmy * d_nrmy_d_rvx;
d_nsuu_d_rvy += d_nsuu_d_nrmy * d_nrmy_d_rvy;
d_nsuu_d_rvz += d_nsuu_d_nrmy * d_nrmy_d_rvz;
d_nsuu_d_rux += d_nsuu_d_nrmx * d_nrmx_d_rux;
d_nsuu_d_ruy += d_nsuu_d_nrmx * d_nrmx_d_ruy;
d_nsuu_d_ruz += d_nsuu_d_nrmx * d_nrmx_d_ruz;
d_nsuu_d_rvx += d_nsuu_d_nrmx * d_nrmx_d_rvx;
d_nsuu_d_rvy += d_nsuu_d_nrmx * d_nrmx_d_rvy;
d_nsuu_d_rvz += d_nsuu_d_nrmx * d_nrmx_d_rvz;
d_g_d_rvx += d_g_d_RvRv * d_RvRv_d_rvx;
d_g_d_rvy += d_g_d_RvRv * d_RvRv_d_rvy;
d_g_d_rvz += d_g_d_RvRv * d_RvRv_d_rvz;
d_g_d_rux += d_g_d_RuRv * d_RuRv_d_rux;
d_g_d_ruy += d_g_d_RuRv * d_RuRv_d_ruy;
d_g_d_ruz += d_g_d_RuRv * d_RuRv_d_ruz;
d_g_d_rvx += d_g_d_RuRv * d_RuRv_d_rvx;
d_g_d_rvy += d_g_d_RuRv * d_RuRv_d_rvy;
d_g_d_rvz += d_g_d_RuRv * d_RuRv_d_rvz;
d_g_d_rux += d_g_d_RuRu * d_RuRu_d_rux;
d_g_d_ruy += d_g_d_RuRu * d_RuRu_d_ruy;
d_g_d_ruz += d_g_d_RuRu * d_RuRu_d_ruz;
d_d_d_rux += d_d_d_RuRv * d_RuRv_d_rux;
d_d_d_ruy += d_d_d_RuRv * d_RuRv_d_ruy;
d_d_d_ruz += d_d_d_RuRv * d_RuRv_d_ruz;
d_d_d_rvx += d_d_d_RuRv * d_RuRv_d_rvx;
d_d_d_rvy += d_d_d_RuRv * d_RuRv_d_rvy;
d_d_d_rvz += d_d_d_RuRv * d_RuRv_d_rvz;
d_d_d_rux += d_d_d_nsuv * d_nsuv_d_rux;
d_d_d_ruy += d_d_d_nsuv * d_nsuv_d_ruy;
d_d_d_ruz += d_d_d_nsuv * d_nsuv_d_ruz;
d_d_d_rvx += d_d_d_nsuv * d_nsuv_d_rvx;
d_d_d_rvy += d_d_d_nsuv * d_nsuv_d_rvy;
d_d_d_rvz += d_d_d_nsuv * d_nsuv_d_rvz;
d_d_d_suvx += d_d_d_nsuv * d_nsuv_d_suvx;
d_d_d_suvy += d_d_d_nsuv * d_nsuv_d_suvy;
d_d_d_suvz += d_d_d_nsuv * d_nsuv_d_suvz;
d_d_d_rux += d_d_d_RuRu * d_RuRu_d_rux;
d_d_d_ruy += d_d_d_RuRu * d_RuRu_d_ruy;
d_d_d_ruz += d_d_d_RuRu * d_RuRu_d_ruz;
d_d_d_rux += d_d_d_nsvv * d_nsvv_d_rux;
d_d_d_ruy += d_d_d_nsvv * d_nsvv_d_ruy;
d_d_d_ruz += d_d_d_nsvv * d_nsvv_d_ruz;
d_d_d_rvx += d_d_d_nsvv * d_nsvv_d_rvx;
d_d_d_rvy += d_d_d_nsvv * d_nsvv_d_rvy;
d_d_d_rvz += d_d_d_nsvv * d_nsvv_d_rvz;
d_d_d_svvx += d_d_d_nsvv * d_nsvv_d_svvx;
d_d_d_svvy += d_d_d_nsvv * d_nsvv_d_svvy;
d_d_d_svvz += d_d_d_nsvv * d_nsvv_d_svvz;
d_d_d_rux += d_d_d_g * d_g_d_rux;
d_d_d_ruy += d_d_d_g * d_g_d_ruy;
d_d_d_ruz += d_d_d_g * d_g_d_ruz;
d_d_d_rvx += d_d_d_g * d_g_d_rvx;
d_d_d_rvy += d_d_d_g * d_g_d_rvy;
d_d_d_rvz += d_d_d_g * d_g_d_rvz;
d_c_d_rux += d_c_d_RuRv * d_RuRv_d_rux;
d_c_d_ruy += d_c_d_RuRv * d_RuRv_d_ruy;
d_c_d_ruz += d_c_d_RuRv * d_RuRv_d_ruz;
d_c_d_rvx += d_c_d_RuRv * d_RuRv_d_rvx;
d_c_d_rvy += d_c_d_RuRv * d_RuRv_d_rvy;
d_c_d_rvz += d_c_d_RuRv * d_RuRv_d_rvz;
d_c_d_rux += d_c_d_nsuv * d_nsuv_d_rux;
d_c_d_ruy += d_c_d_nsuv * d_nsuv_d_ruy;
d_c_d_ruz += d_c_d_nsuv * d_nsuv_d_ruz;
d_c_d_rvx += d_c_d_nsuv * d_nsuv_d_rvx;
d_c_d_rvy += d_c_d_nsuv * d_nsuv_d_rvy;
d_c_d_rvz += d_c_d_nsuv * d_nsuv_d_rvz;
d_c_d_suvx += d_c_d_nsuv * d_nsuv_d_suvx;
d_c_d_suvy += d_c_d_nsuv * d_nsuv_d_suvy;
d_c_d_suvz += d_c_d_nsuv * d_nsuv_d_suvz;
d_c_d_rux += d_c_d_RuRu * d_RuRu_d_rux;
d_c_d_ruy += d_c_d_RuRu * d_RuRu_d_ruy;
d_c_d_ruz += d_c_d_RuRu * d_RuRu_d_ruz;
d_c_d_rux += d_c_d_nsuu * d_nsuu_d_rux;
d_c_d_ruy += d_c_d_nsuu * d_nsuu_d_ruy;
d_c_d_ruz += d_c_d_nsuu * d_nsuu_d_ruz;
d_c_d_rvx += d_c_d_nsuu * d_nsuu_d_rvx;
d_c_d_rvy += d_c_d_nsuu * d_nsuu_d_rvy;
d_c_d_rvz += d_c_d_nsuu * d_nsuu_d_rvz;
d_c_d_suux += d_c_d_nsuu * d_nsuu_d_suux;
d_c_d_suuy += d_c_d_nsuu * d_nsuu_d_suuy;
d_c_d_suuz += d_c_d_nsuu * d_nsuu_d_suuz;
d_c_d_rux += d_c_d_g * d_g_d_rux;
d_c_d_ruy += d_c_d_g * d_g_d_ruy;
d_c_d_ruz += d_c_d_g * d_g_d_ruz;
d_c_d_rvx += d_c_d_g * d_g_d_rvx;
d_c_d_rvy += d_c_d_g * d_g_d_rvy;
d_c_d_rvz += d_c_d_g * d_g_d_rvz;
d_b_d_rux += d_b_d_RuRv * d_RuRv_d_rux;
d_b_d_ruy += d_b_d_RuRv * d_RuRv_d_ruy;
d_b_d_ruz += d_b_d_RuRv * d_RuRv_d_ruz;
d_b_d_rvx += d_b_d_RuRv * d_RuRv_d_rvx;
d_b_d_rvy += d_b_d_RuRv * d_RuRv_d_rvy;
d_b_d_rvz += d_b_d_RuRv * d_RuRv_d_rvz;
d_b_d_rux += d_b_d_nsuv * d_nsuv_d_rux;
d_b_d_ruy += d_b_d_nsuv * d_nsuv_d_ruy;
d_b_d_ruz += d_b_d_nsuv * d_nsuv_d_ruz;
d_b_d_rvx += d_b_d_nsuv * d_nsuv_d_rvx;
d_b_d_rvy += d_b_d_nsuv * d_nsuv_d_rvy;
d_b_d_rvz += d_b_d_nsuv * d_nsuv_d_rvz;
d_b_d_suvx += d_b_d_nsuv * d_nsuv_d_suvx;
d_b_d_suvy += d_b_d_nsuv * d_nsuv_d_suvy;
d_b_d_suvz += d_b_d_nsuv * d_nsuv_d_suvz;
d_b_d_rvx += d_b_d_RvRv * d_RvRv_d_rvx;
d_b_d_rvy += d_b_d_RvRv * d_RvRv_d_rvy;
d_b_d_rvz += d_b_d_RvRv * d_RvRv_d_rvz;
d_b_d_rux += d_b_d_nsvv * d_nsvv_d_rux;
d_b_d_ruy += d_b_d_nsvv * d_nsvv_d_ruy;
d_b_d_ruz += d_b_d_nsvv * d_nsvv_d_ruz;
d_b_d_rvx += d_b_d_nsvv * d_nsvv_d_rvx;
d_b_d_rvy += d_b_d_nsvv * d_nsvv_d_rvy;
d_b_d_rvz += d_b_d_nsvv * d_nsvv_d_rvz;
d_b_d_svvx += d_b_d_nsvv * d_nsvv_d_svvx;
d_b_d_svvy += d_b_d_nsvv * d_nsvv_d_svvy;
d_b_d_svvz += d_b_d_nsvv * d_nsvv_d_svvz;
d_b_d_rux += d_b_d_g * d_g_d_rux;
d_b_d_ruy += d_b_d_g * d_g_d_ruy;
d_b_d_ruz += d_b_d_g * d_g_d_ruz;
d_b_d_rvx += d_b_d_g * d_g_d_rvx;
d_b_d_rvy += d_b_d_g * d_g_d_rvy;
d_b_d_rvz += d_b_d_g * d_g_d_rvz;
d_a_d_rux += d_a_d_RuRv * d_RuRv_d_rux;
d_a_d_ruy += d_a_d_RuRv * d_RuRv_d_ruy;
d_a_d_ruz += d_a_d_RuRv * d_RuRv_d_ruz;
d_a_d_rvx += d_a_d_RuRv * d_RuRv_d_rvx;
d_a_d_rvy += d_a_d_RuRv * d_RuRv_d_rvy;
d_a_d_rvz += d_a_d_RuRv * d_RuRv_d_rvz;
d_a_d_rux += d_a_d_nsuv * d_nsuv_d_rux;
d_a_d_ruy += d_a_d_nsuv * d_nsuv_d_ruy;
d_a_d_ruz += d_a_d_nsuv * d_nsuv_d_ruz;
d_a_d_rvx += d_a_d_nsuv * d_nsuv_d_rvx;
d_a_d_rvy += d_a_d_nsuv * d_nsuv_d_rvy;
d_a_d_rvz += d_a_d_nsuv * d_nsuv_d_rvz;
d_a_d_suvx += d_a_d_nsuv * d_nsuv_d_suvx;
d_a_d_suvy += d_a_d_nsuv * d_nsuv_d_suvy;
d_a_d_suvz += d_a_d_nsuv * d_nsuv_d_suvz;
d_a_d_rvx += d_a_d_RvRv * d_RvRv_d_rvx;
d_a_d_rvy += d_a_d_RvRv * d_RvRv_d_rvy;
d_a_d_rvz += d_a_d_RvRv * d_RvRv_d_rvz;
d_a_d_rux += d_a_d_nsuu * d_nsuu_d_rux;
d_a_d_ruy += d_a_d_nsuu * d_nsuu_d_ruy;
d_a_d_ruz += d_a_d_nsuu * d_nsuu_d_ruz;
d_a_d_rvx += d_a_d_nsuu * d_nsuu_d_rvx;
d_a_d_rvy += d_a_d_nsuu * d_nsuu_d_rvy;
d_a_d_rvz += d_a_d_nsuu * d_nsuu_d_rvz;
d_a_d_suux += d_a_d_nsuu * d_nsuu_d_suux;
d_a_d_suuy += d_a_d_nsuu * d_nsuu_d_suuy;
d_a_d_suuz += d_a_d_nsuu * d_nsuu_d_suuz;
d_a_d_rux += d_a_d_g * d_g_d_rux;
d_a_d_ruy += d_a_d_g * d_g_d_ruy;
d_a_d_ruz += d_a_d_g * d_g_d_ruz;
d_a_d_rvx += d_a_d_g * d_g_d_rvx;
d_a_d_rvy += d_a_d_g * d_g_d_rvy;
d_a_d_rvz += d_a_d_g * d_g_d_rvz;
d_c1_d_rux += d_c1_d_d * d_d_d_rux;
d_c1_d_ruy += d_c1_d_d * d_d_d_ruy;
d_c1_d_ruz += d_c1_d_d * d_d_d_ruz;
d_c1_d_rvx += d_c1_d_d * d_d_d_rvx;
d_c1_d_rvy += d_c1_d_d * d_d_d_rvy;
d_c1_d_rvz += d_c1_d_d * d_d_d_rvz;
d_c1_d_suvx += d_c1_d_d * d_d_d_suvx;
d_c1_d_suvy += d_c1_d_d * d_d_d_suvy;
d_c1_d_suvz += d_c1_d_d * d_d_d_suvz;
d_c1_d_svvx += d_c1_d_d * d_d_d_svvx;
d_c1_d_svvy += d_c1_d_d * d_d_d_svvy;
d_c1_d_svvz += d_c1_d_d * d_d_d_svvz;
d_c1_d_rux += d_c1_d_c * d_c_d_rux;
d_c1_d_ruy += d_c1_d_c * d_c_d_ruy;
d_c1_d_ruz += d_c1_d_c * d_c_d_ruz;
d_c1_d_rvx += d_c1_d_c * d_c_d_rvx;
d_c1_d_rvy += d_c1_d_c * d_c_d_rvy;
d_c1_d_rvz += d_c1_d_c * d_c_d_rvz;
d_c1_d_suux += d_c1_d_c * d_c_d_suux;
d_c1_d_suuy += d_c1_d_c * d_c_d_suuy;
d_c1_d_suuz += d_c1_d_c * d_c_d_suuz;
d_c1_d_suvx += d_c1_d_c * d_c_d_suvx;
d_c1_d_suvy += d_c1_d_c * d_c_d_suvy;
d_c1_d_suvz += d_c1_d_c * d_c_d_suvz;
d_c1_d_rux += d_c1_d_b * d_b_d_rux;
d_c1_d_ruy += d_c1_d_b * d_b_d_ruy;
d_c1_d_ruz += d_c1_d_b * d_b_d_ruz;
d_c1_d_rvx += d_c1_d_b * d_b_d_rvx;
d_c1_d_rvy += d_c1_d_b * d_b_d_rvy;
d_c1_d_rvz += d_c1_d_b * d_b_d_rvz;
d_c1_d_suvx += d_c1_d_b * d_b_d_suvx;
d_c1_d_suvy += d_c1_d_b * d_b_d_suvy;
d_c1_d_suvz += d_c1_d_b * d_b_d_suvz;
d_c1_d_svvx += d_c1_d_b * d_b_d_svvx;
d_c1_d_svvy += d_c1_d_b * d_b_d_svvy;
d_c1_d_svvz += d_c1_d_b * d_b_d_svvz;
d_c1_d_rux += d_c1_d_a * d_a_d_rux;
d_c1_d_ruy += d_c1_d_a * d_a_d_ruy;
d_c1_d_ruz += d_c1_d_a * d_a_d_ruz;
d_c1_d_rvx += d_c1_d_a * d_a_d_rvx;
d_c1_d_rvy += d_c1_d_a * d_a_d_rvy;
d_c1_d_rvz += d_c1_d_a * d_a_d_rvz;
d_c1_d_suux += d_c1_d_a * d_a_d_suux;
d_c1_d_suuy += d_c1_d_a * d_a_d_suuy;
d_c1_d_suuz += d_c1_d_a * d_a_d_suuz;
d_c1_d_suvx += d_c1_d_a * d_a_d_suvx;
d_c1_d_suvy += d_c1_d_a * d_a_d_suvy;
d_c1_d_suvz += d_c1_d_a * d_a_d_suvz;
d_c2_d_rux += d_c2_d_d * d_d_d_rux;
d_c2_d_ruy += d_c2_d_d * d_d_d_ruy;
d_c2_d_ruz += d_c2_d_d * d_d_d_ruz;
d_c2_d_rvx += d_c2_d_d * d_d_d_rvx;
d_c2_d_rvy += d_c2_d_d * d_d_d_rvy;
d_c2_d_rvz += d_c2_d_d * d_d_d_rvz;
d_c2_d_suvx += d_c2_d_d * d_d_d_suvx;
d_c2_d_suvy += d_c2_d_d * d_d_d_suvy;
d_c2_d_suvz += d_c2_d_d * d_d_d_suvz;
d_c2_d_svvx += d_c2_d_d * d_d_d_svvx;
d_c2_d_svvy += d_c2_d_d * d_d_d_svvy;
d_c2_d_svvz += d_c2_d_d * d_d_d_svvz;
d_c2_d_rux += d_c2_d_c * d_c_d_rux;
d_c2_d_ruy += d_c2_d_c * d_c_d_ruy;
d_c2_d_ruz += d_c2_d_c * d_c_d_ruz;
d_c2_d_rvx += d_c2_d_c * d_c_d_rvx;
d_c2_d_rvy += d_c2_d_c * d_c_d_rvy;
d_c2_d_rvz += d_c2_d_c * d_c_d_rvz;
d_c2_d_suux += d_c2_d_c * d_c_d_suux;
d_c2_d_suuy += d_c2_d_c * d_c_d_suuy;
d_c2_d_suuz += d_c2_d_c * d_c_d_suuz;
d_c2_d_suvx += d_c2_d_c * d_c_d_suvx;
d_c2_d_suvy += d_c2_d_c * d_c_d_suvy;
d_c2_d_suvz += d_c2_d_c * d_c_d_suvz;
d_c2_d_rux += d_c2_d_b * d_b_d_rux;
d_c2_d_ruy += d_c2_d_b * d_b_d_ruy;
d_c2_d_ruz += d_c2_d_b * d_b_d_ruz;
d_c2_d_rvx += d_c2_d_b * d_b_d_rvx;
d_c2_d_rvy += d_c2_d_b * d_b_d_rvy;
d_c2_d_rvz += d_c2_d_b * d_b_d_rvz;
d_c2_d_suvx += d_c2_d_b * d_b_d_suvx;
d_c2_d_suvy += d_c2_d_b * d_b_d_suvy;
d_c2_d_suvz += d_c2_d_b * d_b_d_suvz;
d_c2_d_svvx += d_c2_d_b * d_b_d_svvx;
d_c2_d_svvy += d_c2_d_b * d_b_d_svvy;
d_c2_d_svvz += d_c2_d_b * d_b_d_svvz;
d_c2_d_rux += d_c2_d_a * d_a_d_rux;
d_c2_d_ruy += d_c2_d_a * d_a_d_ruy;
d_c2_d_ruz += d_c2_d_a * d_a_d_ruz;
d_c2_d_rvx += d_c2_d_a * d_a_d_rvx;
d_c2_d_rvy += d_c2_d_a * d_a_d_rvy;
d_c2_d_rvz += d_c2_d_a * d_a_d_rvz;
d_c2_d_suux += d_c2_d_a * d_a_d_suux;
d_c2_d_suuy += d_c2_d_a * d_a_d_suuy;
d_c2_d_suuz += d_c2_d_a * d_a_d_suuz;
d_c2_d_suvx += d_c2_d_a * d_a_d_suvx;
d_c2_d_suvy += d_c2_d_a * d_a_d_suvy;
d_c2_d_suvz += d_c2_d_a * d_a_d_suvz;
d_e_d_rux += d_e_d_c2 * d_c2_d_rux;
d_e_d_ruy += d_e_d_c2 * d_c2_d_ruy;
d_e_d_ruz += d_e_d_c2 * d_c2_d_ruz;
d_e_d_rvx += d_e_d_c2 * d_c2_d_rvx;
d_e_d_rvy += d_e_d_c2 * d_c2_d_rvy;
d_e_d_rvz += d_e_d_c2 * d_c2_d_rvz;
d_e_d_suux += d_e_d_c2 * d_c2_d_suux;
d_e_d_suuy += d_e_d_c2 * d_c2_d_suuy;
d_e_d_suuz += d_e_d_c2 * d_c2_d_suuz;
d_e_d_suvx += d_e_d_c2 * d_c2_d_suvx;
d_e_d_suvy += d_e_d_c2 * d_c2_d_suvy;
d_e_d_suvz += d_e_d_c2 * d_c2_d_suvz;
d_e_d_svvx += d_e_d_c2 * d_c2_d_svvx;
d_e_d_svvy += d_e_d_c2 * d_c2_d_svvy;
d_e_d_svvz += d_e_d_c2 * d_c2_d_svvz;
d_e_d_rux += d_e_d_c1 * d_c1_d_rux;
d_e_d_ruy += d_e_d_c1 * d_c1_d_ruy;
d_e_d_ruz += d_e_d_c1 * d_c1_d_ruz;
d_e_d_rvx += d_e_d_c1 * d_c1_d_rvx;
d_e_d_rvy += d_e_d_c1 * d_c1_d_rvy;
d_e_d_rvz += d_e_d_c1 * d_c1_d_rvz;
d_e_d_suux += d_e_d_c1 * d_c1_d_suux;
d_e_d_suuy += d_e_d_c1 * d_c1_d_suuy;
d_e_d_suuz += d_e_d_c1 * d_c1_d_suuz;
d_e_d_suvx += d_e_d_c1 * d_c1_d_suvx;
d_e_d_suvy += d_e_d_c1 * d_c1_d_suvy;
d_e_d_suvz += d_e_d_c1 * d_c1_d_suvz;
d_e_d_svvx += d_e_d_c1 * d_c1_d_svvx;
d_e_d_svvy += d_e_d_c1 * d_c1_d_svvy;
d_e_d_svvz += d_e_d_c1 * d_c1_d_svvz;
d_e_d_rux += d_e_d_g * d_g_d_rux;
d_e_d_ruy += d_e_d_g * d_g_d_ruy;
d_e_d_ruz += d_e_d_g * d_g_d_ruz;
d_e_d_rvx += d_e_d_g * d_g_d_rvx;
d_e_d_rvy += d_e_d_g * d_g_d_rvy;
d_e_d_rvz += d_e_d_g * d_g_d_rvz;
		

		// dE/du 			


		for( int p = 0; p < np; p++ )
		{
			gr[3*cp[p]+0] += d_e_d_rux * ceff_map_du[p] * alpha_x; 
			gr[3*cp[p]+1] += d_e_d_ruy * ceff_map_du[p] * alpha_y; 
			gr[3*cp[p]+2] += d_e_d_ruz * ceff_map_du[p] * alpha_z; 
			
			gr[3*cp[p]+0] += d_e_d_rvx * ceff_map_dv[p] * alpha_x; 
			gr[3*cp[p]+1] += d_e_d_rvy * ceff_map_dv[p] * alpha_y; 
			gr[3*cp[p]+2] += d_e_d_rvz * ceff_map_dv[p] * alpha_z; 
			
			gr[3*cp[p]+0] += d_e_d_suux * ceff_map_duu[p] * alpha_x; 
			gr[3*cp[p]+1] += d_e_d_suuy * ceff_map_duu[p] * alpha_y; 
			gr[3*cp[p]+2] += d_e_d_suuz * ceff_map_duu[p] * alpha_z; 
			
			gr[3*cp[p]+0] += d_e_d_suvx * ceff_map_duv[p] * alpha_x; 
			gr[3*cp[p]+1] += d_e_d_suvy * ceff_map_duv[p] * alpha_y; 
			gr[3*cp[p]+2] += d_e_d_suvz * ceff_map_duv[p] * alpha_z; 
			
			gr[3*cp[p]+0] += d_e_d_svvx * ceff_map_dvv[p] * alpha_x; 
			gr[3*cp[p]+1] += d_e_d_svvy * ceff_map_dvv[p] * alpha_y; 
			gr[3*cp[p]+2] += d_e_d_svvz * ceff_map_dvv[p] * alpha_z; 
			
			dedr[3*p+0] += d_e_d_rux * ceff_map_du[p];
			dedr[3*p+1] += d_e_d_ruy * ceff_map_du[p];
			dedr[3*p+2] += d_e_d_ruz * ceff_map_du[p];
			                                         
			dedr[3*p+0] += d_e_d_rvx * ceff_map_dv[p];
			dedr[3*p+1] += d_e_d_rvy * ceff_map_dv[p];
			dedr[3*p+2] += d_e_d_rvz * ceff_map_dv[p];
			                                           
			dedr[3*p+0] += d_e_d_suux * ceff_map_duu[p]; 
			dedr[3*p+1] += d_e_d_suuy * ceff_map_duu[p]; 
			dedr[3*p+2] += d_e_d_suuz * ceff_map_duu[p]; 
			                                           
			dedr[3*p+0] += d_e_d_suvx * ceff_map_duv[p]; 
			dedr[3*p+1] += d_e_d_suvy * ceff_map_duv[p]; 
			dedr[3*p+2] += d_e_d_suvz * ceff_map_duv[p]; 
			                                           
			dedr[3*p+0] += d_e_d_svvx * ceff_map_dvv[p]; 
			dedr[3*p+1] += d_e_d_svvy * ceff_map_dvv[p]; 
			dedr[3*p+2] += d_e_d_svvz * ceff_map_dvv[p]; 


		}
		
		for( int p = 0; p < np; p++ )
		{
			gr[3*nv+0] += dedr[3*p+0] * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			gr[3*nv+1] += dedr[3*p+1] * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			gr[3*nv+2] += dedr[3*p+2] * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
		}
	}
	else
	{
		int frm = (f-nf_faces)*nf_irr_pts;

		int t = theIrregularFormulas[frm].tri;
		int i = theIrregularFormulas[frm].vertex;
		int ed = theIrregularFormulas[frm].edge;
		int *cp = theIrregularFormulas[frm].cp;
		int np = theIrregularFormulas[frm].ncoor;

		int val = theVertices[i].valence;

		irr_kernel *theKernel = kernels[val]; 
		int ncoords_base = val + 6;
		
		double fu = u;
		double fv = v;					

		int domain = theKernel->domain(fu,fv);
		if( domain-1 >= theKernel->ndomains )
		{
			printf("LOGICAL ERROR evaluating R,nrm on irregular vertex face.\n");
			exit(1);
		}	
		double *theMap = theKernel->get_map( &fu, &fv );
		double u_u = u, u_v = 0, v_u = 0, v_v = v;
		theKernel->get_map_transform( &u_u, &u_v, &v_u, &v_v );
		double domain_scale = pow(2, domain);

		double u = fu;
		double v = fv;
		double w = 1 - u - v;

		if( u +v > 1.0 || u < 0 || v < 0 )
		{
			printf("ERROR: outside domain.\n");
			exit(1);
		}

		double u2 = u*u;
		double u3 = u*u*u;
		double u4 = u*u*u*u;
		
		double v2 = v*v;
		double v3 = v*v*v;
		double v4 = v*v*v*v;
		
		double w2 = w*w;
		double w3 = w*w*w;
		double w4 = w*w*w*w;
			
		
		// 8 : 0
		// 7 : 1
		// 4 : 2
		// 5 : 3
		// 9 : 4
		// 12 : 5
		// 11 : 6
		// 10 : 7
		// 6 : 8
		// 3 : 9
		// 1 : 10
		// 2 : 11

	

		double n1 = (1.0/12.0)*(u4+2*u3*v); 
		double n2 = (1.0/12.0)*(u4+2*u3*w); 
		double n3 = (1.0/12.0)*(u4+2*u3*w+6*u3*v+6*u2*v*w+12*u2*v2+6*u*v2*w+6*u*v3+2*v3*w+v4); 
		double n4 = (1.0/12.0)*(6*u4+24*u3*w+24*u2*w2+8*u*w3+w4+24*u3*v+60*u2*v*w+36*u*v*w2+6*v*w3+24*u2*v2+36*u*v2*w+12*v2*w2+8*u*v3+6*v3*w+v4); 
		double n5 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+2*u3*v+6*u2*v*w+6*u*v*w2+2*v*w3); 
		double n6 = (1.0/12.0)*(2*u*v3+v4); 
		double n7 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+8*u3*v+36*u2*v*w+36*u*v*w2+8*v*w3+24*u2*v2+60*u*v2*w+24*v2*w2+24*u*v3+24*v3*w+6*v4); 
		double n8 = (1.0/12.0)*(u4+8*u3*w+24*u2*w2+24*u*w3+6*w4+6*u3*v+36*u2*v*w+60*u*v*w2+24*v*w3+12*u2*v2+36*u*v2*w+24*v2*w2+6*u*v3+8*v3*w+v4); 
		double n9 = (1.0/12.0)*(2*u*w3+w4); 
		double n10 = (1.0/12.0)*(2*v3*w+v4); 
		double n11 = (1.0/12.0)*(2*u*w3+w4+6*u*v*w2+6*v*w3+6*u*v2*w+12*v2*w2+2*u*v3+6*v3*w+v4); 
		double n12 = (1.0/12.0)*(w4+2*v*w3);

		double du_1 = Power(u,3)/3. + (Power(u,2)*v)/2.;
		double du_2 = Power(u,2)/2. - Power(u,3)/3. - (Power(u,2)*v)/2.; 	
		double du_3 =Power(u,2)/2. - Power(u,3)/3. + u*v - (Power(u,2)*v)/2. + Power(v,2)/2. - Power(v,3)/6.;
		double du_4 = 0.3333333333333333 + u - Power(u,2) - Power(u,3)/3. + v/2. - u*v - (Power(u,2)*v)/2. - Power(v,2) + Power(v,3)/3.;	
		double du_5 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. - v/2. + Power(u,2)*v + Power(v,2)/2. - Power(v,3)/6.;
		double du_6 = Power(v,3)/6.;	
		double du_7 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. + v/2. - 2*u*v + Power(u,2)*v - Power(v,2)/2. - Power(v,3)/6.;	
		double du_8 = -2*u + 2*Power(u,2) - Power(u,3)/3. - v + 2*u*v - (Power(u,2)*v)/2. + Power(v,2) - Power(v,3)/6.;	
		double du_9 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. + v/2. - (Power(u,2)*v)/2. - Power(v,2)/2. + Power(v,3)/6.;	
		double du_10 = -Power(v,3)/6.;
		double du_11 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. - v/2. + u*v - (Power(u,2)*v)/2. + Power(v,3)/3.;	
		double du_12 = 	-0.3333333333333333 + u - Power(u,2) + Power(u,3)/3. + v/2. - u*v + (Power(u,2)*v)/2. - Power(v,3)/6.;

		double dv_1 = Power(u,3)/6.;
		double dv_2 = -Power(u,3)/6.;
		double dv_3 = Power(u,2)/2. - Power(u,3)/6. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_4 = 0.16666666666666666 + u/2. - Power(u,2)/2. - Power(u,3)/6. - 2*u*v - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
		double dv_5 = -0.16666666666666666 - u/2. + Power(u,3)/3. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_6 = (u*Power(v,2))/2. + Power(v,3)/3.;
		double dv_7 = 0.3333333333333333 + u/2. - Power(u,2) + Power(u,3)/3. + v - u*v - Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_8 = -u + Power(u,2) - Power(u,3)/6. - 2*v + 2*u*v + 2*Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_9 = -0.3333333333333333 + u/2. - Power(u,3)/6. + v - u*v - Power(v,2) + (u*Power(v,2))/2. + Power(v,3)/3.;
		double dv_10 = Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_11 = 0.16666666666666666 - u/2. + Power(u,2)/2. - Power(u,3)/6. - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
		double dv_12 = -0.16666666666666666 + u/2. - Power(u,2)/2. + Power(u,3)/6. + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
	
		double d_uu_1 = Power(u,2) + u*v;
		double d_uu_2 = u - Power(u,2) - u*v;
		double d_uu_3 = u - Power(u,2) + v - u*v;
		double d_uu_4 = 1 - 2*u - Power(u,2) - v - u*v;
		double d_uu_5 = -2*u + 2*Power(u,2) + 2*u*v;
		double d_uu_6 = 0;
		double d_uu_7 = -2*u + 2*Power(u,2) - 2*v + 2*u*v;
		double d_uu_8 = -2 + 4*u - Power(u,2) + 2*v - u*v;
		double d_uu_9 = u - Power(u,2) - u*v;
		double d_uu_10 = 0;
		double d_uu_11 = u - Power(u,2) + v - u*v;
		double d_uu_12 = 1 - 2*u + Power(u,2) - v + u*v;
		
		double d_uv_1 = Power(u,2)/2.;
		double d_uv_2 = -Power(u,2)/2.;
		double d_uv_3 = u - Power(u,2)/2. + v - Power(v,2)/2.;
		double d_uv_4 = 0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
		double d_uv_5 = -0.5 + Power(u,2) + v - Power(v,2)/2.;
		double d_uv_6 = Power(v,2)/2.;
		double d_uv_7 = 0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
		double d_uv_8 = -1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
		double d_uv_9 = 0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
		double d_uv_10 = -Power(v,2)/2.;
		double d_uv_11 = -0.5 + u - Power(u,2)/2. + Power(v,2);
		double d_uv_12 = 0.5 - u + Power(u,2)/2. - Power(v,2)/2.;
		
		double d_vv_1 = 0;
		double d_vv_2 = 0;
		double d_vv_3 = u + v - u*v - Power(v,2);
		double d_vv_4 = -2*u - 2*v + 2*u*v + 2*Power(v,2);
		double d_vv_5 = u + v - u*v - Power(v,2);
		double d_vv_6 = u*v + Power(v,2);
		double d_vv_7 = 1 - u - 2*v - u*v - Power(v,2);
		double d_vv_8 = -2 + 2*u + 4*v - u*v - Power(v,2);
		double d_vv_9 = 1 - u - 2*v + u*v + Power(v,2);
		double d_vv_10 = v - u*v - Power(v,2);
		double d_vv_11 = -2*v + 2*u*v + 2*Power(v,2);
		double d_vv_12 = v - u*v - Power(v,2);

		double d_uuu_1 = 2*u + v;
		double d_uuu_2 = 1-2*u-v; //u - Power(u,2) - u*v;
		double d_uuu_3 = 1-2*u-v;//u - Power(u,2) + v - u*v;
		double d_uuu_4 = -2-2*u-v;//1 - 2*u - Power(u,2) - v - u*v;
		double d_uuu_5 = -2+4*u+2*v;//-2*u + 2*Power(u,2) + 2*u*v;
		double d_uuu_6 = 0;
		double d_uuu_7 = -2+4*u+2*v;//-2*u + 2*Power(u,2) - 2*v + 2*u*v;
		double d_uuu_8 = 4-2*u-v;//-2 + 4*u - Power(u,2) + 2*v - u*v;
		double d_uuu_9 = 1-2*u-v;//u - Power(u,2) - u*v;
		double d_uuu_10 = 0;
		double d_uuu_11 = 1-2*u-v;//u - Power(u,2) + v - u*v;
		double d_uuu_12 = -2+2*u+v;//1 - 2*u + Power(u,2) - v + u*v;

		double d_uuv_1 = u;//Power(u,2)/2.;
		double d_uuv_2 = -u;//-Power(u,2)/2.;
		double d_uuv_3 = 1-u;//u - Power(u,2)/2. + v - Power(v,2)/2.;
		double d_uuv_4 = -1-u;//0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
		double d_uuv_5 = 2*u;//-0.5 + Power(u,2) + v - Power(v,2)/2.;
		double d_uuv_6 = 0;//Power(v,2)/2.;
		double d_uuv_7 = -2+2*u;//0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
		double d_uuv_8 = 2-u;//-1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
		double d_uuv_9 = -u;//0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
		double d_uuv_10 = 0;//-Power(v,2)/2.;
		double d_uuv_11 = 1-u;//-0.5 + u - Power(u,2)/2. + Power(v,2);
		double d_uuv_12 = -1+u;//0.5 - u + Power(u,2)/2. - Power(v,2)/2.;

		double d_uvv_1 = 0;
		double d_uvv_2 = 0;
		double d_uvv_3 = 1-v;//u + v - u*v - Power(v,2);
		double d_uvv_4 = -2+2*v;//-2*u - 2*v + 2*u*v + 2*Power(v,2);
		double d_uvv_5 = 1-v;//u + v - u*v - Power(v,2);
		double d_uvv_6 = v;//u*v + Power(v,2);
		double d_uvv_7 = -1-v;//1 - u - 2*v - u*v - Power(v,2);
		double d_uvv_8 = 2-v;//-2 + 2*u + 4*v - u*v - Power(v,2);
		double d_uvv_9 = -1+v;//1 - u - 2*v + u*v + Power(v,2);
		double d_uvv_10 = -v;//v - u*v - Power(v,2);
		double d_uvv_11 = 2*v;//-2*v + 2*u*v + 2*Power(v,2);
		double d_uvv_12 = -v;//v - u*v - Power(v,2);
		
		double d_vvv_1 = 0;
		double d_vvv_2 = 0;
		double d_vvv_3 = 1-u-2*v;//u + v - u*v - Power(v,2);
		double d_vvv_4 = -2+2*u+4*v;//-2*u - 2*v + 2*u*v + 2*Power(v,2);
		double d_vvv_5 = 1-u-2*v;//u + v - u*v - Power(v,2);
		double d_vvv_6 = u+2*v;//u*v + Power(v,2);
		double d_vvv_7 = -2-u-2*v;//1 - u - 2*v - u*v - Power(v,2);
		double d_vvv_8 = 4-u-2*v;//-2 + 2*u + 4*v - u*v - Power(v,2);
		double d_vvv_9 = -2+u+2*v;//1 - u - 2*v + u*v + Power(v,2);
		double d_vvv_10 = 1-u-2*v;//v - u*v - Power(v,2);
		double d_vvv_11 = -2+2*u+4*v;//-2*v + 2*u*v + 2*Power(v,2);
		double d_vvv_12 = 1-u-2*v;//v - u*v - Power(v,2);
		
		double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
		double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
		double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
		
		double ceff_map_duu[12] = { d_uu_8, d_uu_7, d_uu_4, d_uu_5, d_uu_9, d_uu_12, d_uu_11, d_uu_10, d_uu_6, d_uu_3, d_uu_1, d_uu_2 };
		double ceff_map_duv[12] = { d_uv_8, d_uv_7, d_uv_4, d_uv_5, d_uv_9, d_uv_12, d_uv_11, d_uv_10, d_uv_6, d_uv_3, d_uv_1, d_uv_2 };
		double ceff_map_dvv[12] = { d_vv_8, d_vv_7, d_vv_4, d_vv_5, d_vv_9, d_vv_12, d_vv_11, d_vv_10, d_vv_6, d_vv_3, d_vv_1, d_vv_2 };
		
		double ceff_map_duuu[12] = { d_uuu_8, d_uuu_7, d_uuu_4, d_uuu_5, d_uuu_9, d_uuu_12, d_uuu_11, d_uuu_10, d_uuu_6, d_uuu_3, d_uuu_1, d_uuu_2 };
		double ceff_map_duuv[12] = { d_uuv_8, d_uuv_7, d_uuv_4, d_uuv_5, d_uuv_9, d_uuv_12, d_uuv_11, d_uuv_10, d_uuv_6, d_uuv_3, d_uuv_1, d_uuv_2 };
		double ceff_map_duvv[12] = { d_uvv_8, d_uvv_7, d_uvv_4, d_uvv_5, d_uvv_9, d_uvv_12, d_uvv_11, d_uvv_10, d_uvv_6, d_uvv_3, d_uvv_1, d_uvv_2 };
		double ceff_map_dvvv[12] = { d_vvv_8, d_vvv_7, d_vvv_4, d_vvv_5, d_vvv_9, d_vvv_12, d_vvv_11, d_vvv_10, d_vvv_6, d_vvv_3, d_vvv_1, d_vvv_2 };

		double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
		double nrm[3]={0,0,0}; 
		double tSuuu[3] = {0,0,0};
		double tSuuv[3] = {0,0,0};
		double tSuvv[3] = {0,0,0};
		double tSvvv[3] = {0,0,0};

		int *cset = theVertices[i].irr_coord_set + ed * ncoords_base;

		for( int x = 0; x < ncoords_base; x++ )
		{
			for( int y = 0; y < 12; y++ )
			{
				R[0] += (r[3*cset[x]+0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map[y] * alpha_x;
				R[1] += (r[3*cset[x]+1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map[y] * alpha_y;
				R[2] += (r[3*cset[x]+2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map[y] * alpha_z;
				
				Ru[0] += (r[3*cset[x]+0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_du[y] * alpha_x * u_u;
				Ru[1] += (r[3*cset[x]+1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_du[y] * alpha_y * u_u;
				Ru[2] += (r[3*cset[x]+2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_du[y] * alpha_z * u_u;
				
				Rv[0] += (r[3*cset[x]+0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_dv[y] * alpha_x * v_v;
				Rv[1] += (r[3*cset[x]+1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_dv[y] * alpha_y * v_v;
				Rv[2] += (r[3*cset[x]+2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_dv[y] * alpha_z * v_v;
				
				tSuu[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_duu[y] * alpha_x * u_u * u_u;
				tSuu[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_duu[y] * alpha_y * u_u * u_u;
				tSuu[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_duu[y] * alpha_z * u_u * u_u;
				
				tSuv[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_duv[y] * alpha_x * u_u * v_v;
				tSuv[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_duv[y] * alpha_y * u_u * v_v;
				tSuv[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_duv[y] * alpha_z * u_u * v_v;
				
				tSvv[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_dvv[y] * alpha_x * v_v * v_v;
				tSvv[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_dvv[y] * alpha_y * v_v * v_v;
				tSvv[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_dvv[y] * alpha_z * v_v * v_v;
		
				tSuuu[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_duuu[y] * alpha_x * u_u * u_u * u_u;
				tSuuu[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_duuu[y] * alpha_y * u_u * u_u * u_u;
				tSuuu[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_duuu[y] * alpha_z * u_u * u_u * u_u;
				                                                                                                                              
				tSuuv[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_duuv[y] * alpha_x * u_u * u_u * v_v;
				tSuuv[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_duuv[y] * alpha_y * u_u * u_u * v_v;
				tSuuv[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_duuv[y] * alpha_z * u_u * u_u * v_v;
				                                                                                                                              
				tSuvv[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_duvv[y] * alpha_x * u_u * v_v * v_v;
				tSuvv[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_duvv[y] * alpha_y * u_u * v_v * v_v;
				tSuvv[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_duvv[y] * alpha_z * u_u * v_v * v_v;
				
				tSvvv[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_dvvv[y] * alpha_x * v_v * v_v * v_v;
				tSvvv[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_dvvv[y] * alpha_y * v_v * v_v * v_v;
				tSvvv[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_dvvv[y] * alpha_z * v_v * v_v * v_v;

			}
		}
		
		cross( Ru, Rv, nrm );
		normalize(nrm);


		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		double g = sqrt(RuRu*RvRv-RuRv*RuRv);

		double e1,e2;
      
		double nsuu = tSuu[0] * nrm[0] + tSuu[1] * nrm[1] + tSuu[2] * nrm[2];
		double nsuv = tSuv[0] * nrm[0] + tSuv[1] * nrm[1] + tSuv[2] * nrm[2];
		double nsvv = tSvv[0] * nrm[0] + tSvv[1] * nrm[1] + tSvv[2] * nrm[2];

		double Stot = (nsuu * RvRv + nsvv * RuRu -2*nsuv*RuRv)/(g*g);

		double Sop[4] = { 1.0/(g*g) * (nsuu * RvRv  - nsuv * RuRv), (1.0/(g*g)) * ( nsuv*RvRv-nsvv*RuRv),
				  1.0/(g*g) * (nsuv * RuRu  - nsuu * RuRv), (1.0/(g*g)) * ( nsvv*RuRu-nsuv*RuRv) };

		double a = Sop[0];		
		double b = Sop[1];
		double c = Sop[2];
		double d = Sop[3];

		double c0 = theFormulas[frm].c0;
		e1 = -0.5*(a+d-sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
		e2 = -0.5*(a+d+sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
	

double d_nrmz_d_rux=0,d_nrmz_d_ruy=0,d_nrmz_d_ruz=0,d_nrmz_d_rvx=0,d_nrmz_d_rvy=0,d_nrmz_d_rvz=0,d_nrmy_d_rux=0,d_nrmy_d_ruy=0,d_nrmy_d_ruz=0,d_nrmy_d_rvx=0,d_nrmy_d_rvy=0,d_nrmy_d_rvz=0,d_nrmx_d_rux=0,d_nrmx_d_ruy=0,d_nrmx_d_ruz=0,d_nrmx_d_rvx=0,d_nrmx_d_rvy=0,d_nrmx_d_rvz=0,d_RuRu_d_rux=0,d_RuRu_d_ruy=0,d_RuRu_d_ruz=0,d_nsvv_d_rux=0,d_nsvv_d_ruy=0,d_nsvv_d_ruz=0,d_nsvv_d_rvx=0,d_nsvv_d_rvy=0,d_nsvv_d_rvz=0,d_nsvv_d_svvx=0,d_nsvv_d_svvy=0,d_nsvv_d_svvz=0,d_RuRv_d_rux=0,d_RuRv_d_ruy=0,d_RuRv_d_ruz=0,d_RuRv_d_rvx=0,d_RuRv_d_rvy=0,d_RuRv_d_rvz=0,d_nsuv_d_rux=0,d_nsuv_d_ruy=0,d_nsuv_d_ruz=0,d_nsuv_d_rvx=0,d_nsuv_d_rvy=0,d_nsuv_d_rvz=0,d_nsuv_d_suvx=0,d_nsuv_d_suvy=0,d_nsuv_d_suvz=0,d_RvRv_d_rvx=0,d_RvRv_d_rvy=0,d_RvRv_d_rvz=0,d_nsuu_d_rux=0,d_nsuu_d_ruy=0,d_nsuu_d_ruz=0,d_nsuu_d_rvx=0,d_nsuu_d_rvy=0,d_nsuu_d_rvz=0,d_nsuu_d_suux=0,d_nsuu_d_suuy=0,d_nsuu_d_suuz=0,d_g_d_rux=0,d_g_d_ruy=0,d_g_d_ruz=0,d_g_d_rvx=0,d_g_d_rvy=0,d_g_d_rvz=0,d_e_d_rux=0,d_e_d_ruy=0,d_e_d_ruz=0,d_e_d_rvx=0,d_e_d_rvy=0,d_e_d_rvz=0,d_e_d_suux=0,d_e_d_suuy=0,d_e_d_suuz=0,d_e_d_suvx=0,d_e_d_suvy=0,d_e_d_suvz=0,d_e_d_svvx=0,d_e_d_svvy=0,d_e_d_svvz=0,d_c2_d_rux=0,d_c2_d_ruy=0,d_c2_d_ruz=0,d_c2_d_rvx=0,d_c2_d_rvy=0,d_c2_d_rvz=0,d_c2_d_suux=0,d_c2_d_suuy=0,d_c2_d_suuz=0,d_c2_d_suvx=0,d_c2_d_suvy=0,d_c2_d_suvz=0,d_c2_d_svvx=0,d_c2_d_svvy=0,d_c2_d_svvz=0,d_d_d_rux=0,d_d_d_ruy=0,d_d_d_ruz=0,d_d_d_rvx=0,d_d_d_rvy=0,d_d_d_rvz=0,d_d_d_suvx=0,d_d_d_suvy=0,d_d_d_suvz=0,d_d_d_svvx=0,d_d_d_svvy=0,d_d_d_svvz=0,d_c_d_rux=0,d_c_d_ruy=0,d_c_d_ruz=0,d_c_d_rvx=0,d_c_d_rvy=0,d_c_d_rvz=0,d_c_d_suux=0,d_c_d_suuy=0,d_c_d_suuz=0,d_c_d_suvx=0,d_c_d_suvy=0,d_c_d_suvz=0,d_b_d_rux=0,d_b_d_ruy=0,d_b_d_ruz=0,d_b_d_rvx=0,d_b_d_rvy=0,d_b_d_rvz=0,d_b_d_suvx=0,d_b_d_suvy=0,d_b_d_suvz=0,d_b_d_svvx=0,d_b_d_svvy=0,d_b_d_svvz=0,d_a_d_rux=0,d_a_d_ruy=0,d_a_d_ruz=0,d_a_d_rvx=0,d_a_d_rvy=0,d_a_d_rvz=0,d_a_d_suux=0,d_a_d_suuy=0,d_a_d_suuz=0,d_a_d_suvx=0,d_a_d_suvy=0,d_a_d_suvz=0,d_c1_d_rux=0,d_c1_d_ruy=0,d_c1_d_ruz=0,d_c1_d_rvx=0,d_c1_d_rvy=0,d_c1_d_rvz=0,d_c1_d_suux=0,d_c1_d_suuy=0,d_c1_d_suuz=0,d_c1_d_suvx=0,d_c1_d_suvy=0,d_c1_d_suvz=0,d_c1_d_svvx=0,d_c1_d_svvy=0,d_c1_d_svvz=0,d_g_d_RvRv=0,d_g_d_RuRv=0,d_g_d_RuRu=0,d_nsvv_d_nrmz=0,d_nsvv_d_nrmy=0,d_nsvv_d_nrmx=0,d_nsuv_d_nrmz=0,d_nsuv_d_nrmy=0,d_nsuv_d_nrmx=0,d_nsuu_d_nrmz=0,d_nsuu_d_nrmy=0,d_nsuu_d_nrmx=0,d_d_d_RuRv=0,d_d_d_nsuv=0,d_d_d_RuRu=0,d_d_d_nsvv=0,d_d_d_g=0,d_c_d_RuRv=0,d_c_d_nsuv=0,d_c_d_RuRu=0,d_c_d_nsuu=0,d_c_d_g=0,d_b_d_RuRv=0,d_b_d_nsuv=0,d_b_d_RvRv=0,d_b_d_nsvv=0,d_b_d_g=0,d_a_d_RuRv=0,d_a_d_nsuv=0,d_a_d_RvRv=0,d_a_d_nsuu=0,d_a_d_g=0,d_e_d_c2=0,d_e_d_c1=0,d_e_d_g=0,d_c2_d_d=0,d_c2_d_c=0,d_c2_d_b=0,d_c2_d_a=0,d_c1_d_d=0,d_c1_d_c=0,d_c1_d_b=0,d_c1_d_a=0,junk;
		// the basic variables.





		d_g_d_RuRu = RvRv/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
		d_g_d_RuRv = -(RuRv/Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
		d_g_d_RvRv = RuRu/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  

		d_RuRu_d_rux = 2*Ru[0];
		d_RuRu_d_ruy = 2*Ru[1];
		d_RuRu_d_ruz = 2*Ru[2];
		
		d_RuRv_d_rux = Rv[0];
		d_RuRv_d_ruy = Rv[1];
		d_RuRv_d_ruz = Rv[2];
		
		d_RuRv_d_rvx = Ru[0];
		d_RuRv_d_rvy = Ru[1];
		d_RuRv_d_rvz = Ru[2];
		
		d_RvRv_d_rvx = 2*Rv[0];
		d_RvRv_d_rvy = 2*Rv[1];
		d_RvRv_d_rvz = 2*Rv[2];

		// intermediates.
		d_a_d_g = (-2*(-(nsuv*RuRv) + nsuu*RvRv))/Power(g,3);
		d_a_d_nsuu = RvRv/Power(g,2);
		d_a_d_RvRv = nsuu/Power(g,2);
		d_a_d_nsuv = -(RuRv/Power(g,2));
		d_a_d_RuRv = -(nsuv/Power(g,2));

		d_b_d_g = (-2*(-(nsvv*RuRv) + nsuv*RvRv))/Power(g,3);
		d_b_d_nsvv = -(RuRv/Power(g,2));
		d_b_d_RvRv = nsuv/Power(g,2);
		d_b_d_nsuv = RvRv/Power(g,2);
		d_b_d_RuRv = -(nsvv/Power(g,2));

		d_c_d_g = (-2*(nsuv*RuRu - nsuu*RuRv))/Power(g,3);
		d_c_d_nsuu = -(RuRv/Power(g,2));
		d_c_d_RuRu = nsuv/Power(g,2);
		d_c_d_nsuv = RuRu/Power(g,2);
		d_c_d_RuRv = -(nsuu/Power(g,2));

		d_d_d_g = (-2*(nsvv*RuRu - nsuv*RuRv))/Power(g,3);
		d_d_d_nsvv = RuRu/Power(g,2);
		d_d_d_RuRu = nsvv/Power(g,2);
		d_d_d_nsuv = -(RuRv/Power(g,2));
		d_d_d_RuRv = -(nsuv/Power(g,2));

		d_c1_d_a =-0.5*(1 - (2*a - 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
		d_c1_d_b =-(-1.*c)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2)) ;
		d_c1_d_c =-(-1.*b)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2)) ;
		d_c1_d_d =-0.5*(1 - (-2*a + 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
	
		d_c2_d_a = -0.5*(1 + (2*a - 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
		d_c2_d_b = -(1.*c)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2));
		d_c2_d_c = -(1.*b)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2));
		d_c2_d_d = -0.5*(1 + (-2*a + 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));


		// this is the term we need to gradient:
		// en += 0.5*kc * p_area * ( c1 + c2 - p_c0 ) * ( c1 + c2 - p_c0 );
		d_e_d_c1 = kc * p_area * ( e1 + e2 - p_c0);
		d_e_d_c2 = kc * p_area * ( e1 + e2 - p_c0);
		d_e_d_g = 0;

		// that's it.

		double fac = Power(Power(Ru[1]*Rv[0] - Ru[0]*Rv[1],2) + Power(Ru[2]*Rv[0] - Ru[0]*Rv[2],2) + Power(Ru[2]*Rv[1] - Ru[1]*Rv[2],2),1.5);

		d_nrmx_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac; 
		d_nrmx_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2)))/fac;
		d_nrmx_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac;
		
		d_nrmy_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2))))/fac;
		d_nrmy_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac; 
		d_nrmy_d_ruz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2])) + Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac;

		d_nrmz_d_rux = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2])/fac;
		d_nrmz_d_ruy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;
		d_nrmz_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2))) + (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;

		d_nrmx_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac; 
		d_nrmx_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2])))/fac;
		d_nrmx_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac;
		
		d_nrmy_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
		d_nrmy_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2])))/fac;
		d_nrmy_d_rvz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[0]*Ru[1]*Rv[0] - Power(Ru[0],2)*Rv[1] + Ru[2]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2])))/fac;

		d_nrmz_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;
		d_nrmz_d_rvy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[0]*Ru[2]*Rv[0] - Power(Ru[0],2)*Rv[2] + Ru[1]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
		d_nrmz_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;

		double d_nrmx_d_u = d_nrmx_d_rux * tSuu[0] + d_nrmx_d_ruy * tSuu[1] + d_nrmx_d_ruz * tSuu[2] +
				    d_nrmx_d_rvx * tSuv[0] + d_nrmx_d_rvy * tSuv[1] + d_nrmx_d_rvz * tSuv[2];
		double d_nrmx_d_v = d_nrmx_d_rux * tSuv[0] + d_nrmx_d_ruy * tSuv[1] + d_nrmx_d_ruz * tSuv[2] +
				    d_nrmx_d_rvx * tSvv[0] + d_nrmx_d_rvy * tSvv[1] + d_nrmx_d_rvz * tSvv[2];
		double d_nrmy_d_u = d_nrmy_d_rux * tSuu[0] + d_nrmy_d_ruy * tSuu[1] + d_nrmy_d_ruz * tSuu[2] +
				    d_nrmy_d_rvx * tSuv[0] + d_nrmy_d_rvy * tSuv[1] + d_nrmy_d_rvz * tSuv[2];
		double d_nrmy_d_v = d_nrmy_d_rux * tSuv[0] + d_nrmy_d_ruy * tSuv[1] + d_nrmy_d_ruz * tSuv[2] +
				    d_nrmy_d_rvx * tSvv[0] + d_nrmy_d_rvy * tSvv[1] + d_nrmy_d_rvz * tSvv[2];
		double d_nrmz_d_u = d_nrmz_d_rux * tSuu[0] + d_nrmz_d_ruy * tSuu[1] + d_nrmz_d_ruz * tSuu[2] +
				    d_nrmz_d_rvx * tSuv[0] + d_nrmz_d_rvy * tSuv[1] + d_nrmz_d_rvz * tSuv[2];
		double d_nrmz_d_v = d_nrmz_d_rux * tSuv[0] + d_nrmz_d_ruy * tSuv[1] + d_nrmz_d_ruz * tSuv[2] +
				    d_nrmz_d_rvx * tSvv[0] + d_nrmz_d_rvy * tSvv[1] + d_nrmz_d_rvz * tSvv[2];

		double d_nsuu_d_u = tSuuu[0] * nrm[0] + tSuuu[1] * nrm[1] + tSuuu[2] * nrm[2];
		       d_nsuu_d_u += tSuu[0] * d_nrmx_d_u + tSuu[1] * d_nrmy_d_u + tSuu[2] * d_nrmz_d_u;
		double d_nsuu_d_v = tSuuv[0] * nrm[0] + tSuuv[1] * nrm[1] + tSuuv[2] * nrm[2];
		       d_nsuu_d_v += tSuu[0] * d_nrmx_d_v + tSuu[1] * d_nrmy_d_v + tSuu[2] * d_nrmz_d_v;
		
		double d_nsuv_d_u = tSuuv[0] * nrm[0] + tSuuv[1] * nrm[1] + tSuuv[2] * nrm[2];
		       d_nsuv_d_u += tSuv[0] * d_nrmx_d_u + tSuv[1] * d_nrmy_d_u + tSuv[2] * d_nrmz_d_u;
		double d_nsuv_d_v = tSuvv[0] * nrm[0] + tSuvv[1] * nrm[1] + tSuvv[2] * nrm[2];
		       d_nsuv_d_v += tSuv[0] * d_nrmx_d_v + tSuv[1] * d_nrmy_d_v + tSuv[2] * d_nrmz_d_v;

		double d_nsvv_d_u = tSuvv[0] * nrm[0] + tSuvv[1] * nrm[1] + tSuvv[2] * nrm[2];
		       d_nsvv_d_u += tSvv[0] * d_nrmx_d_u + tSvv[1] * d_nrmy_d_u + tSvv[2] * d_nrmz_d_u;
		double d_nsvv_d_v = tSvvv[0] * nrm[0] + tSvvv[1] * nrm[1] + tSvvv[2] * nrm[2];
		       d_nsvv_d_v += tSvv[0] * d_nrmx_d_v + tSvv[1] * d_nrmy_d_v + tSvv[2] * d_nrmz_d_v;

		double d_RuRu_d_u = 2*Ru[0] * tSuu[0] + 2*Ru[1] * tSuu[1] + 2 * Ru[2] * tSuu[2];
		double d_RuRu_d_v = 2*Ru[0] * tSuv[0] + 2*Ru[1] * tSuv[1] + 2 * Ru[2] * tSuv[2];

		double d_RuRv_d_u =  Ru[0] * tSuv[0] + Ru[1] * tSuv[1] + Ru[2] * tSuv[2];
		       d_RuRv_d_u += Rv[0] * tSuu[0] + Rv[1] * tSuu[1] + Rv[2] * tSuu[2];
		
		double d_RuRv_d_v =  Ru[0] * tSvv[0] + Ru[1] * tSvv[1] + Ru[2] * tSvv[2];
		       d_RuRv_d_v += Rv[0] * tSuv[0] + Rv[1] * tSuv[1] + Rv[2] * tSuv[2];

		double d_RvRv_d_u = 2*Rv[0] * tSuv[0] + 2*Rv[1] * tSuv[1] + 2 * Rv[2] * tSuv[2];
		double d_RvRv_d_v = 2*Rv[0] * tSvv[0] + 2*Rv[1] * tSvv[1] + 2 * Rv[2] * tSvv[2];

		double d_g_d_u = 0;
		double d_g_d_v = 0;

		d_g_d_u += d_g_d_RuRu * 2* Ru[0] * tSuu[0] + d_g_d_RuRu * 2* Ru[1] * tSuu[1] + d_g_d_RuRu * 2* Ru[2] * tSuu[2];
		d_g_d_v += d_g_d_RuRu * 2* Ru[0] * tSuv[0] + d_g_d_RuRu * 2* Ru[1] * tSuv[1] + d_g_d_RuRu * 2* Ru[2] * tSuv[2];

		// d_ru dru_du
		d_g_d_u += d_g_d_RuRv * Rv[0] * tSuu[0] + d_g_d_RuRv * Rv[1] * tSuu[1] + d_g_d_RuRv * Rv[2] * tSuu[2];
		// d_rv_drv_du
		d_g_d_u += d_g_d_RuRv * Ru[0] * tSuv[0] + d_g_d_RuRv * Ru[1] * tSuv[1] + d_g_d_RuRv * Ru[2] * tSuv[2];	
		// d_ru dru_dv
		d_g_d_v += d_g_d_RuRv * Rv[0] * tSuv[0] + d_g_d_RuRv * Rv[1] * tSuv[1] + d_g_d_RuRv * Rv[2] * tSuv[2];
		// d_rv drv_dv
		d_g_d_v += d_g_d_RuRv * Ru[0] * tSvv[0] + d_g_d_RuRv * Ru[1] * tSvv[1] + d_g_d_RuRv * Ru[2] * tSvv[2];
		
		d_g_d_u += d_g_d_RvRv * 2* Rv[0] * tSuv[0] + d_g_d_RvRv * 2* Rv[1] * tSuv[1] + d_g_d_RvRv * 2* Rv[2] * tSuv[2];
		d_g_d_v += d_g_d_RvRv * 2* Rv[0] * tSvv[0] + d_g_d_RvRv * 2* Rv[1] * tSvv[1] + d_g_d_RvRv * 2* Rv[2] * tSvv[2];

//			d_a_d_g = (-2*(-(nsuv*RuRv) + nsuu*RvRv))/Power(g,3);
//			d_a_d_nsuu = RvRv/Power(g,2);
//			d_a_d_RvRv = nsuu/Power(g,2);
//			d_a_d_nsuv = -(RuRv/Power(g,2));
//			d_a_d_RuRv = -(nsuv/Power(g,2));
		
		double d_a_d_u = d_a_d_nsuu * d_nsuu_d_u + d_a_d_nsuv * d_nsuv_d_u + d_a_d_RuRv * d_RuRv_d_u + d_a_d_RvRv * d_RvRv_d_u + d_a_d_g * d_g_d_u;
		double d_a_d_v = d_a_d_nsuu * d_nsuu_d_v + d_a_d_nsuv * d_nsuv_d_v + d_a_d_RuRv * d_RuRv_d_v + d_a_d_RvRv * d_RvRv_d_v + d_a_d_g * d_g_d_v;
		double d_b_d_u = d_b_d_nsvv * d_nsvv_d_u + d_b_d_nsuv * d_nsuv_d_u + d_b_d_RuRv * d_RuRv_d_u + d_b_d_RvRv * d_RvRv_d_u + d_b_d_g * d_g_d_u;
		double d_b_d_v = d_b_d_nsvv * d_nsvv_d_v + d_b_d_nsuv * d_nsuv_d_v + d_b_d_RuRv * d_RuRv_d_v + d_b_d_RvRv * d_RvRv_d_v + d_b_d_g * d_g_d_v;
		double d_c_d_u = d_c_d_nsuu * d_nsuu_d_u + d_c_d_nsuv * d_nsuv_d_u + d_c_d_RuRv * d_RuRv_d_u + d_c_d_RuRu * d_RuRu_d_u + d_c_d_g * d_g_d_u;
		double d_c_d_v = d_c_d_nsuu * d_nsuu_d_v + d_c_d_nsuv * d_nsuv_d_v + d_c_d_RuRv * d_RuRv_d_v + d_c_d_RuRu * d_RuRu_d_v + d_c_d_g * d_g_d_v;
		double d_d_d_u = d_d_d_nsvv * d_nsvv_d_u + d_d_d_nsuv * d_nsuv_d_u + d_d_d_RuRv * d_RuRv_d_u + d_d_d_RuRu * d_RuRu_d_u + d_d_d_g * d_g_d_u;
		double d_d_d_v = d_d_d_nsvv * d_nsvv_d_v + d_d_d_nsuv * d_nsuv_d_v + d_d_d_RuRv * d_RuRv_d_v + d_d_d_RuRu * d_RuRu_d_v + d_d_d_g * d_g_d_v;

		double d_c1_d_u = d_c1_d_a * d_a_d_u + d_c1_d_b * d_b_d_u + d_c1_d_c * d_c_d_u + d_c1_d_d * d_d_d_u;
		double d_c2_d_u = d_c2_d_a * d_a_d_u + d_c2_d_b * d_b_d_u + d_c2_d_c * d_c_d_u + d_c2_d_d * d_d_d_u;

		double d_c1_d_v = d_c1_d_a * d_a_d_v + d_c1_d_b * d_b_d_v + d_c1_d_c * d_c_d_v + d_c1_d_d * d_d_d_v;
		double d_c2_d_v = d_c2_d_a * d_a_d_v + d_c2_d_b * d_b_d_v + d_c2_d_c * d_c_d_v + d_c2_d_d * d_d_d_v;
		double d_e_d_u = d_e_d_c1 * d_c1_d_u + d_e_d_c2 * d_c2_d_u;
		double d_e_d_v = d_e_d_c1 * d_c1_d_v + d_e_d_c2 * d_c2_d_v;

		p_uv_g[0] += d_e_d_u ;
		p_uv_g[1] += d_e_d_v ;

//			printf("u %le v %le a: %.14le da_du: %le nrmx: %.14le dnrmx_du %.14le\n", u, v, b, d_b_d_u, nrm[0], d_nrmy_d_u );

		double e = 0.5* kc * p_area * ( e1 + e2 - p_c0 ) * ( e1 + e2 - p_c0 );
//			printf("e %le at uv: %le %le  pg: %le %le\n", e, p_uv[2*pid+0], p_uv[2*pid+1], d_e_d_u, d_e_d_v );
		
		d_nsuu_d_nrmx = tSuu[0];
		d_nsuu_d_nrmy = tSuu[1];
		d_nsuu_d_nrmz = tSuu[2];
		
		d_nsuv_d_nrmx = tSuv[0];
		d_nsuv_d_nrmy = tSuv[1];
		d_nsuv_d_nrmz = tSuv[2];

		d_nsvv_d_nrmx = tSvv[0];
		d_nsvv_d_nrmy = tSvv[1];
		d_nsvv_d_nrmz = tSvv[2];

		d_nsuu_d_suux = nrm[0];
		d_nsuu_d_suuy = nrm[1];
		d_nsuu_d_suuz = nrm[2];
		
		d_nsuv_d_suvx = nrm[0];
		d_nsuv_d_suvy = nrm[1];
		d_nsuv_d_suvz = nrm[2];

		d_nsvv_d_svvx = nrm[0];
		d_nsvv_d_svvy = nrm[1];
		d_nsvv_d_svvz = nrm[2];



d_nsvv_d_rux += d_nsvv_d_nrmz * d_nrmz_d_rux;
d_nsvv_d_ruy += d_nsvv_d_nrmz * d_nrmz_d_ruy;
d_nsvv_d_ruz += d_nsvv_d_nrmz * d_nrmz_d_ruz;
d_nsvv_d_rvx += d_nsvv_d_nrmz * d_nrmz_d_rvx;
d_nsvv_d_rvy += d_nsvv_d_nrmz * d_nrmz_d_rvy;
d_nsvv_d_rvz += d_nsvv_d_nrmz * d_nrmz_d_rvz;
d_nsvv_d_rux += d_nsvv_d_nrmy * d_nrmy_d_rux;
d_nsvv_d_ruy += d_nsvv_d_nrmy * d_nrmy_d_ruy;
d_nsvv_d_ruz += d_nsvv_d_nrmy * d_nrmy_d_ruz;
d_nsvv_d_rvx += d_nsvv_d_nrmy * d_nrmy_d_rvx;
d_nsvv_d_rvy += d_nsvv_d_nrmy * d_nrmy_d_rvy;
d_nsvv_d_rvz += d_nsvv_d_nrmy * d_nrmy_d_rvz;
d_nsvv_d_rux += d_nsvv_d_nrmx * d_nrmx_d_rux;
d_nsvv_d_ruy += d_nsvv_d_nrmx * d_nrmx_d_ruy;
d_nsvv_d_ruz += d_nsvv_d_nrmx * d_nrmx_d_ruz;
d_nsvv_d_rvx += d_nsvv_d_nrmx * d_nrmx_d_rvx;
d_nsvv_d_rvy += d_nsvv_d_nrmx * d_nrmx_d_rvy;
d_nsvv_d_rvz += d_nsvv_d_nrmx * d_nrmx_d_rvz;
d_nsuv_d_rux += d_nsuv_d_nrmz * d_nrmz_d_rux;
d_nsuv_d_ruy += d_nsuv_d_nrmz * d_nrmz_d_ruy;
d_nsuv_d_ruz += d_nsuv_d_nrmz * d_nrmz_d_ruz;
d_nsuv_d_rvx += d_nsuv_d_nrmz * d_nrmz_d_rvx;
d_nsuv_d_rvy += d_nsuv_d_nrmz * d_nrmz_d_rvy;
d_nsuv_d_rvz += d_nsuv_d_nrmz * d_nrmz_d_rvz;
d_nsuv_d_rux += d_nsuv_d_nrmy * d_nrmy_d_rux;
d_nsuv_d_ruy += d_nsuv_d_nrmy * d_nrmy_d_ruy;
d_nsuv_d_ruz += d_nsuv_d_nrmy * d_nrmy_d_ruz;
d_nsuv_d_rvx += d_nsuv_d_nrmy * d_nrmy_d_rvx;
d_nsuv_d_rvy += d_nsuv_d_nrmy * d_nrmy_d_rvy;
d_nsuv_d_rvz += d_nsuv_d_nrmy * d_nrmy_d_rvz;
d_nsuv_d_rux += d_nsuv_d_nrmx * d_nrmx_d_rux;
d_nsuv_d_ruy += d_nsuv_d_nrmx * d_nrmx_d_ruy;
d_nsuv_d_ruz += d_nsuv_d_nrmx * d_nrmx_d_ruz;
d_nsuv_d_rvx += d_nsuv_d_nrmx * d_nrmx_d_rvx;
d_nsuv_d_rvy += d_nsuv_d_nrmx * d_nrmx_d_rvy;
d_nsuv_d_rvz += d_nsuv_d_nrmx * d_nrmx_d_rvz;
d_nsuu_d_rux += d_nsuu_d_nrmz * d_nrmz_d_rux;
d_nsuu_d_ruy += d_nsuu_d_nrmz * d_nrmz_d_ruy;
d_nsuu_d_ruz += d_nsuu_d_nrmz * d_nrmz_d_ruz;
d_nsuu_d_rvx += d_nsuu_d_nrmz * d_nrmz_d_rvx;
d_nsuu_d_rvy += d_nsuu_d_nrmz * d_nrmz_d_rvy;
d_nsuu_d_rvz += d_nsuu_d_nrmz * d_nrmz_d_rvz;
d_nsuu_d_rux += d_nsuu_d_nrmy * d_nrmy_d_rux;
d_nsuu_d_ruy += d_nsuu_d_nrmy * d_nrmy_d_ruy;
d_nsuu_d_ruz += d_nsuu_d_nrmy * d_nrmy_d_ruz;
d_nsuu_d_rvx += d_nsuu_d_nrmy * d_nrmy_d_rvx;
d_nsuu_d_rvy += d_nsuu_d_nrmy * d_nrmy_d_rvy;
d_nsuu_d_rvz += d_nsuu_d_nrmy * d_nrmy_d_rvz;
d_nsuu_d_rux += d_nsuu_d_nrmx * d_nrmx_d_rux;
d_nsuu_d_ruy += d_nsuu_d_nrmx * d_nrmx_d_ruy;
d_nsuu_d_ruz += d_nsuu_d_nrmx * d_nrmx_d_ruz;
d_nsuu_d_rvx += d_nsuu_d_nrmx * d_nrmx_d_rvx;
d_nsuu_d_rvy += d_nsuu_d_nrmx * d_nrmx_d_rvy;
d_nsuu_d_rvz += d_nsuu_d_nrmx * d_nrmx_d_rvz;
d_g_d_rvx += d_g_d_RvRv * d_RvRv_d_rvx;
d_g_d_rvy += d_g_d_RvRv * d_RvRv_d_rvy;
d_g_d_rvz += d_g_d_RvRv * d_RvRv_d_rvz;
d_g_d_rux += d_g_d_RuRv * d_RuRv_d_rux;
d_g_d_ruy += d_g_d_RuRv * d_RuRv_d_ruy;
d_g_d_ruz += d_g_d_RuRv * d_RuRv_d_ruz;
d_g_d_rvx += d_g_d_RuRv * d_RuRv_d_rvx;
d_g_d_rvy += d_g_d_RuRv * d_RuRv_d_rvy;
d_g_d_rvz += d_g_d_RuRv * d_RuRv_d_rvz;
d_g_d_rux += d_g_d_RuRu * d_RuRu_d_rux;
d_g_d_ruy += d_g_d_RuRu * d_RuRu_d_ruy;
d_g_d_ruz += d_g_d_RuRu * d_RuRu_d_ruz;
d_d_d_rux += d_d_d_RuRv * d_RuRv_d_rux;
d_d_d_ruy += d_d_d_RuRv * d_RuRv_d_ruy;
d_d_d_ruz += d_d_d_RuRv * d_RuRv_d_ruz;
d_d_d_rvx += d_d_d_RuRv * d_RuRv_d_rvx;
d_d_d_rvy += d_d_d_RuRv * d_RuRv_d_rvy;
d_d_d_rvz += d_d_d_RuRv * d_RuRv_d_rvz;
d_d_d_rux += d_d_d_nsuv * d_nsuv_d_rux;
d_d_d_ruy += d_d_d_nsuv * d_nsuv_d_ruy;
d_d_d_ruz += d_d_d_nsuv * d_nsuv_d_ruz;
d_d_d_rvx += d_d_d_nsuv * d_nsuv_d_rvx;
d_d_d_rvy += d_d_d_nsuv * d_nsuv_d_rvy;
d_d_d_rvz += d_d_d_nsuv * d_nsuv_d_rvz;
d_d_d_suvx += d_d_d_nsuv * d_nsuv_d_suvx;
d_d_d_suvy += d_d_d_nsuv * d_nsuv_d_suvy;
d_d_d_suvz += d_d_d_nsuv * d_nsuv_d_suvz;
d_d_d_rux += d_d_d_RuRu * d_RuRu_d_rux;
d_d_d_ruy += d_d_d_RuRu * d_RuRu_d_ruy;
d_d_d_ruz += d_d_d_RuRu * d_RuRu_d_ruz;
d_d_d_rux += d_d_d_nsvv * d_nsvv_d_rux;
d_d_d_ruy += d_d_d_nsvv * d_nsvv_d_ruy;
d_d_d_ruz += d_d_d_nsvv * d_nsvv_d_ruz;
d_d_d_rvx += d_d_d_nsvv * d_nsvv_d_rvx;
d_d_d_rvy += d_d_d_nsvv * d_nsvv_d_rvy;
d_d_d_rvz += d_d_d_nsvv * d_nsvv_d_rvz;
d_d_d_svvx += d_d_d_nsvv * d_nsvv_d_svvx;
d_d_d_svvy += d_d_d_nsvv * d_nsvv_d_svvy;
d_d_d_svvz += d_d_d_nsvv * d_nsvv_d_svvz;
d_d_d_rux += d_d_d_g * d_g_d_rux;
d_d_d_ruy += d_d_d_g * d_g_d_ruy;
d_d_d_ruz += d_d_d_g * d_g_d_ruz;
d_d_d_rvx += d_d_d_g * d_g_d_rvx;
d_d_d_rvy += d_d_d_g * d_g_d_rvy;
d_d_d_rvz += d_d_d_g * d_g_d_rvz;
d_c_d_rux += d_c_d_RuRv * d_RuRv_d_rux;
d_c_d_ruy += d_c_d_RuRv * d_RuRv_d_ruy;
d_c_d_ruz += d_c_d_RuRv * d_RuRv_d_ruz;
d_c_d_rvx += d_c_d_RuRv * d_RuRv_d_rvx;
d_c_d_rvy += d_c_d_RuRv * d_RuRv_d_rvy;
d_c_d_rvz += d_c_d_RuRv * d_RuRv_d_rvz;
d_c_d_rux += d_c_d_nsuv * d_nsuv_d_rux;
d_c_d_ruy += d_c_d_nsuv * d_nsuv_d_ruy;
d_c_d_ruz += d_c_d_nsuv * d_nsuv_d_ruz;
d_c_d_rvx += d_c_d_nsuv * d_nsuv_d_rvx;
d_c_d_rvy += d_c_d_nsuv * d_nsuv_d_rvy;
d_c_d_rvz += d_c_d_nsuv * d_nsuv_d_rvz;
d_c_d_suvx += d_c_d_nsuv * d_nsuv_d_suvx;
d_c_d_suvy += d_c_d_nsuv * d_nsuv_d_suvy;
d_c_d_suvz += d_c_d_nsuv * d_nsuv_d_suvz;
d_c_d_rux += d_c_d_RuRu * d_RuRu_d_rux;
d_c_d_ruy += d_c_d_RuRu * d_RuRu_d_ruy;
d_c_d_ruz += d_c_d_RuRu * d_RuRu_d_ruz;
d_c_d_rux += d_c_d_nsuu * d_nsuu_d_rux;
d_c_d_ruy += d_c_d_nsuu * d_nsuu_d_ruy;
d_c_d_ruz += d_c_d_nsuu * d_nsuu_d_ruz;
d_c_d_rvx += d_c_d_nsuu * d_nsuu_d_rvx;
d_c_d_rvy += d_c_d_nsuu * d_nsuu_d_rvy;
d_c_d_rvz += d_c_d_nsuu * d_nsuu_d_rvz;
d_c_d_suux += d_c_d_nsuu * d_nsuu_d_suux;
d_c_d_suuy += d_c_d_nsuu * d_nsuu_d_suuy;
d_c_d_suuz += d_c_d_nsuu * d_nsuu_d_suuz;
d_c_d_rux += d_c_d_g * d_g_d_rux;
d_c_d_ruy += d_c_d_g * d_g_d_ruy;
d_c_d_ruz += d_c_d_g * d_g_d_ruz;
d_c_d_rvx += d_c_d_g * d_g_d_rvx;
d_c_d_rvy += d_c_d_g * d_g_d_rvy;
d_c_d_rvz += d_c_d_g * d_g_d_rvz;
d_b_d_rux += d_b_d_RuRv * d_RuRv_d_rux;
d_b_d_ruy += d_b_d_RuRv * d_RuRv_d_ruy;
d_b_d_ruz += d_b_d_RuRv * d_RuRv_d_ruz;
d_b_d_rvx += d_b_d_RuRv * d_RuRv_d_rvx;
d_b_d_rvy += d_b_d_RuRv * d_RuRv_d_rvy;
d_b_d_rvz += d_b_d_RuRv * d_RuRv_d_rvz;
d_b_d_rux += d_b_d_nsuv * d_nsuv_d_rux;
d_b_d_ruy += d_b_d_nsuv * d_nsuv_d_ruy;
d_b_d_ruz += d_b_d_nsuv * d_nsuv_d_ruz;
d_b_d_rvx += d_b_d_nsuv * d_nsuv_d_rvx;
d_b_d_rvy += d_b_d_nsuv * d_nsuv_d_rvy;
d_b_d_rvz += d_b_d_nsuv * d_nsuv_d_rvz;
d_b_d_suvx += d_b_d_nsuv * d_nsuv_d_suvx;
d_b_d_suvy += d_b_d_nsuv * d_nsuv_d_suvy;
d_b_d_suvz += d_b_d_nsuv * d_nsuv_d_suvz;
d_b_d_rvx += d_b_d_RvRv * d_RvRv_d_rvx;
d_b_d_rvy += d_b_d_RvRv * d_RvRv_d_rvy;
d_b_d_rvz += d_b_d_RvRv * d_RvRv_d_rvz;
d_b_d_rux += d_b_d_nsvv * d_nsvv_d_rux;
d_b_d_ruy += d_b_d_nsvv * d_nsvv_d_ruy;
d_b_d_ruz += d_b_d_nsvv * d_nsvv_d_ruz;
d_b_d_rvx += d_b_d_nsvv * d_nsvv_d_rvx;
d_b_d_rvy += d_b_d_nsvv * d_nsvv_d_rvy;
d_b_d_rvz += d_b_d_nsvv * d_nsvv_d_rvz;
d_b_d_svvx += d_b_d_nsvv * d_nsvv_d_svvx;
d_b_d_svvy += d_b_d_nsvv * d_nsvv_d_svvy;
d_b_d_svvz += d_b_d_nsvv * d_nsvv_d_svvz;
d_b_d_rux += d_b_d_g * d_g_d_rux;
d_b_d_ruy += d_b_d_g * d_g_d_ruy;
d_b_d_ruz += d_b_d_g * d_g_d_ruz;
d_b_d_rvx += d_b_d_g * d_g_d_rvx;
d_b_d_rvy += d_b_d_g * d_g_d_rvy;
d_b_d_rvz += d_b_d_g * d_g_d_rvz;
d_a_d_rux += d_a_d_RuRv * d_RuRv_d_rux;
d_a_d_ruy += d_a_d_RuRv * d_RuRv_d_ruy;
d_a_d_ruz += d_a_d_RuRv * d_RuRv_d_ruz;
d_a_d_rvx += d_a_d_RuRv * d_RuRv_d_rvx;
d_a_d_rvy += d_a_d_RuRv * d_RuRv_d_rvy;
d_a_d_rvz += d_a_d_RuRv * d_RuRv_d_rvz;
d_a_d_rux += d_a_d_nsuv * d_nsuv_d_rux;
d_a_d_ruy += d_a_d_nsuv * d_nsuv_d_ruy;
d_a_d_ruz += d_a_d_nsuv * d_nsuv_d_ruz;
d_a_d_rvx += d_a_d_nsuv * d_nsuv_d_rvx;
d_a_d_rvy += d_a_d_nsuv * d_nsuv_d_rvy;
d_a_d_rvz += d_a_d_nsuv * d_nsuv_d_rvz;
d_a_d_suvx += d_a_d_nsuv * d_nsuv_d_suvx;
d_a_d_suvy += d_a_d_nsuv * d_nsuv_d_suvy;
d_a_d_suvz += d_a_d_nsuv * d_nsuv_d_suvz;
d_a_d_rvx += d_a_d_RvRv * d_RvRv_d_rvx;
d_a_d_rvy += d_a_d_RvRv * d_RvRv_d_rvy;
d_a_d_rvz += d_a_d_RvRv * d_RvRv_d_rvz;
d_a_d_rux += d_a_d_nsuu * d_nsuu_d_rux;
d_a_d_ruy += d_a_d_nsuu * d_nsuu_d_ruy;
d_a_d_ruz += d_a_d_nsuu * d_nsuu_d_ruz;
d_a_d_rvx += d_a_d_nsuu * d_nsuu_d_rvx;
d_a_d_rvy += d_a_d_nsuu * d_nsuu_d_rvy;
d_a_d_rvz += d_a_d_nsuu * d_nsuu_d_rvz;
d_a_d_suux += d_a_d_nsuu * d_nsuu_d_suux;
d_a_d_suuy += d_a_d_nsuu * d_nsuu_d_suuy;
d_a_d_suuz += d_a_d_nsuu * d_nsuu_d_suuz;
d_a_d_rux += d_a_d_g * d_g_d_rux;
d_a_d_ruy += d_a_d_g * d_g_d_ruy;
d_a_d_ruz += d_a_d_g * d_g_d_ruz;
d_a_d_rvx += d_a_d_g * d_g_d_rvx;
d_a_d_rvy += d_a_d_g * d_g_d_rvy;
d_a_d_rvz += d_a_d_g * d_g_d_rvz;
d_c1_d_rux += d_c1_d_d * d_d_d_rux;
d_c1_d_ruy += d_c1_d_d * d_d_d_ruy;
d_c1_d_ruz += d_c1_d_d * d_d_d_ruz;
d_c1_d_rvx += d_c1_d_d * d_d_d_rvx;
d_c1_d_rvy += d_c1_d_d * d_d_d_rvy;
d_c1_d_rvz += d_c1_d_d * d_d_d_rvz;
d_c1_d_suvx += d_c1_d_d * d_d_d_suvx;
d_c1_d_suvy += d_c1_d_d * d_d_d_suvy;
d_c1_d_suvz += d_c1_d_d * d_d_d_suvz;
d_c1_d_svvx += d_c1_d_d * d_d_d_svvx;
d_c1_d_svvy += d_c1_d_d * d_d_d_svvy;
d_c1_d_svvz += d_c1_d_d * d_d_d_svvz;
d_c1_d_rux += d_c1_d_c * d_c_d_rux;
d_c1_d_ruy += d_c1_d_c * d_c_d_ruy;
d_c1_d_ruz += d_c1_d_c * d_c_d_ruz;
d_c1_d_rvx += d_c1_d_c * d_c_d_rvx;
d_c1_d_rvy += d_c1_d_c * d_c_d_rvy;
d_c1_d_rvz += d_c1_d_c * d_c_d_rvz;
d_c1_d_suux += d_c1_d_c * d_c_d_suux;
d_c1_d_suuy += d_c1_d_c * d_c_d_suuy;
d_c1_d_suuz += d_c1_d_c * d_c_d_suuz;
d_c1_d_suvx += d_c1_d_c * d_c_d_suvx;
d_c1_d_suvy += d_c1_d_c * d_c_d_suvy;
d_c1_d_suvz += d_c1_d_c * d_c_d_suvz;
d_c1_d_rux += d_c1_d_b * d_b_d_rux;
d_c1_d_ruy += d_c1_d_b * d_b_d_ruy;
d_c1_d_ruz += d_c1_d_b * d_b_d_ruz;
d_c1_d_rvx += d_c1_d_b * d_b_d_rvx;
d_c1_d_rvy += d_c1_d_b * d_b_d_rvy;
d_c1_d_rvz += d_c1_d_b * d_b_d_rvz;
d_c1_d_suvx += d_c1_d_b * d_b_d_suvx;
d_c1_d_suvy += d_c1_d_b * d_b_d_suvy;
d_c1_d_suvz += d_c1_d_b * d_b_d_suvz;
d_c1_d_svvx += d_c1_d_b * d_b_d_svvx;
d_c1_d_svvy += d_c1_d_b * d_b_d_svvy;
d_c1_d_svvz += d_c1_d_b * d_b_d_svvz;
d_c1_d_rux += d_c1_d_a * d_a_d_rux;
d_c1_d_ruy += d_c1_d_a * d_a_d_ruy;
d_c1_d_ruz += d_c1_d_a * d_a_d_ruz;
d_c1_d_rvx += d_c1_d_a * d_a_d_rvx;
d_c1_d_rvy += d_c1_d_a * d_a_d_rvy;
d_c1_d_rvz += d_c1_d_a * d_a_d_rvz;
d_c1_d_suux += d_c1_d_a * d_a_d_suux;
d_c1_d_suuy += d_c1_d_a * d_a_d_suuy;
d_c1_d_suuz += d_c1_d_a * d_a_d_suuz;
d_c1_d_suvx += d_c1_d_a * d_a_d_suvx;
d_c1_d_suvy += d_c1_d_a * d_a_d_suvy;
d_c1_d_suvz += d_c1_d_a * d_a_d_suvz;
d_c2_d_rux += d_c2_d_d * d_d_d_rux;
d_c2_d_ruy += d_c2_d_d * d_d_d_ruy;
d_c2_d_ruz += d_c2_d_d * d_d_d_ruz;
d_c2_d_rvx += d_c2_d_d * d_d_d_rvx;
d_c2_d_rvy += d_c2_d_d * d_d_d_rvy;
d_c2_d_rvz += d_c2_d_d * d_d_d_rvz;
d_c2_d_suvx += d_c2_d_d * d_d_d_suvx;
d_c2_d_suvy += d_c2_d_d * d_d_d_suvy;
d_c2_d_suvz += d_c2_d_d * d_d_d_suvz;
d_c2_d_svvx += d_c2_d_d * d_d_d_svvx;
d_c2_d_svvy += d_c2_d_d * d_d_d_svvy;
d_c2_d_svvz += d_c2_d_d * d_d_d_svvz;
d_c2_d_rux += d_c2_d_c * d_c_d_rux;
d_c2_d_ruy += d_c2_d_c * d_c_d_ruy;
d_c2_d_ruz += d_c2_d_c * d_c_d_ruz;
d_c2_d_rvx += d_c2_d_c * d_c_d_rvx;
d_c2_d_rvy += d_c2_d_c * d_c_d_rvy;
d_c2_d_rvz += d_c2_d_c * d_c_d_rvz;
d_c2_d_suux += d_c2_d_c * d_c_d_suux;
d_c2_d_suuy += d_c2_d_c * d_c_d_suuy;
d_c2_d_suuz += d_c2_d_c * d_c_d_suuz;
d_c2_d_suvx += d_c2_d_c * d_c_d_suvx;
d_c2_d_suvy += d_c2_d_c * d_c_d_suvy;
d_c2_d_suvz += d_c2_d_c * d_c_d_suvz;
d_c2_d_rux += d_c2_d_b * d_b_d_rux;
d_c2_d_ruy += d_c2_d_b * d_b_d_ruy;
d_c2_d_ruz += d_c2_d_b * d_b_d_ruz;
d_c2_d_rvx += d_c2_d_b * d_b_d_rvx;
d_c2_d_rvy += d_c2_d_b * d_b_d_rvy;
d_c2_d_rvz += d_c2_d_b * d_b_d_rvz;
d_c2_d_suvx += d_c2_d_b * d_b_d_suvx;
d_c2_d_suvy += d_c2_d_b * d_b_d_suvy;
d_c2_d_suvz += d_c2_d_b * d_b_d_suvz;
d_c2_d_svvx += d_c2_d_b * d_b_d_svvx;
d_c2_d_svvy += d_c2_d_b * d_b_d_svvy;
d_c2_d_svvz += d_c2_d_b * d_b_d_svvz;
d_c2_d_rux += d_c2_d_a * d_a_d_rux;
d_c2_d_ruy += d_c2_d_a * d_a_d_ruy;
d_c2_d_ruz += d_c2_d_a * d_a_d_ruz;
d_c2_d_rvx += d_c2_d_a * d_a_d_rvx;
d_c2_d_rvy += d_c2_d_a * d_a_d_rvy;
d_c2_d_rvz += d_c2_d_a * d_a_d_rvz;
d_c2_d_suux += d_c2_d_a * d_a_d_suux;
d_c2_d_suuy += d_c2_d_a * d_a_d_suuy;
d_c2_d_suuz += d_c2_d_a * d_a_d_suuz;
d_c2_d_suvx += d_c2_d_a * d_a_d_suvx;
d_c2_d_suvy += d_c2_d_a * d_a_d_suvy;
d_c2_d_suvz += d_c2_d_a * d_a_d_suvz;
d_e_d_rux += d_e_d_c2 * d_c2_d_rux;
d_e_d_ruy += d_e_d_c2 * d_c2_d_ruy;
d_e_d_ruz += d_e_d_c2 * d_c2_d_ruz;
d_e_d_rvx += d_e_d_c2 * d_c2_d_rvx;
d_e_d_rvy += d_e_d_c2 * d_c2_d_rvy;
d_e_d_rvz += d_e_d_c2 * d_c2_d_rvz;
d_e_d_suux += d_e_d_c2 * d_c2_d_suux;
d_e_d_suuy += d_e_d_c2 * d_c2_d_suuy;
d_e_d_suuz += d_e_d_c2 * d_c2_d_suuz;
d_e_d_suvx += d_e_d_c2 * d_c2_d_suvx;
d_e_d_suvy += d_e_d_c2 * d_c2_d_suvy;
d_e_d_suvz += d_e_d_c2 * d_c2_d_suvz;
d_e_d_svvx += d_e_d_c2 * d_c2_d_svvx;
d_e_d_svvy += d_e_d_c2 * d_c2_d_svvy;
d_e_d_svvz += d_e_d_c2 * d_c2_d_svvz;
d_e_d_rux += d_e_d_c1 * d_c1_d_rux;
d_e_d_ruy += d_e_d_c1 * d_c1_d_ruy;
d_e_d_ruz += d_e_d_c1 * d_c1_d_ruz;
d_e_d_rvx += d_e_d_c1 * d_c1_d_rvx;
d_e_d_rvy += d_e_d_c1 * d_c1_d_rvy;
d_e_d_rvz += d_e_d_c1 * d_c1_d_rvz;
d_e_d_suux += d_e_d_c1 * d_c1_d_suux;
d_e_d_suuy += d_e_d_c1 * d_c1_d_suuy;
d_e_d_suuz += d_e_d_c1 * d_c1_d_suuz;
d_e_d_suvx += d_e_d_c1 * d_c1_d_suvx;
d_e_d_suvy += d_e_d_c1 * d_c1_d_suvy;
d_e_d_suvz += d_e_d_c1 * d_c1_d_suvz;
d_e_d_svvx += d_e_d_c1 * d_c1_d_svvx;
d_e_d_svvy += d_e_d_c1 * d_c1_d_svvy;
d_e_d_svvz += d_e_d_c1 * d_c1_d_svvz;
d_e_d_rux += d_e_d_g * d_g_d_rux;
d_e_d_ruy += d_e_d_g * d_g_d_ruy;
d_e_d_ruz += d_e_d_g * d_g_d_ruz;
d_e_d_rvx += d_e_d_g * d_g_d_rvx;
d_e_d_rvy += d_e_d_g * d_g_d_rvy;
d_e_d_rvz += d_e_d_g * d_g_d_rvz;
		

		// dE/du 			

		double dedr[3*np];
		memset( dedr, 0, sizeof(double) * 3 *np );
		for( int p = 0; p < np; p++ )
		for( int y = 0; y < 12; y++ )
		{
			gr[3*cp[p]+0] += d_e_d_rux *  ceff_map_du[y] * theMap[y*ncoords_base+p] * alpha_x; 
			gr[3*cp[p]+1] += d_e_d_ruy *  ceff_map_du[y] * theMap[y*ncoords_base+p] * alpha_y; 
			gr[3*cp[p]+2] += d_e_d_ruz *  ceff_map_du[y] * theMap[y*ncoords_base+p] * alpha_z; 
			
			gr[3*cp[p]+0] += d_e_d_rvx *  ceff_map_dv[y] * theMap[y*ncoords_base+p] * alpha_x; 
			gr[3*cp[p]+1] += d_e_d_rvy *  ceff_map_dv[y] * theMap[y*ncoords_base+p] * alpha_y; 
			gr[3*cp[p]+2] += d_e_d_rvz *  ceff_map_dv[y] * theMap[y*ncoords_base+p] * alpha_z; 
			
			gr[3*cp[p]+0] += d_e_d_suux * ceff_map_duu[y] * theMap[y*ncoords_base+p] * alpha_x; 
			gr[3*cp[p]+1] += d_e_d_suuy * ceff_map_duu[y] * theMap[y*ncoords_base+p] * alpha_y; 
			gr[3*cp[p]+2] += d_e_d_suuz * ceff_map_duu[y] * theMap[y*ncoords_base+p] * alpha_z; 
			
			gr[3*cp[p]+0] += d_e_d_suvx * ceff_map_duv[y] * theMap[y*ncoords_base+p] * alpha_x; 
			gr[3*cp[p]+1] += d_e_d_suvy * ceff_map_duv[y] * theMap[y*ncoords_base+p] * alpha_y; 
			gr[3*cp[p]+2] += d_e_d_suvz * ceff_map_duv[y] * theMap[y*ncoords_base+p] * alpha_z; 
			
			gr[3*cp[p]+0] += d_e_d_svvx * ceff_map_dvv[y] * theMap[y*ncoords_base+p] * alpha_x; 
			gr[3*cp[p]+1] += d_e_d_svvy * ceff_map_dvv[y] * theMap[y*ncoords_base+p] * alpha_y; 
			gr[3*cp[p]+2] += d_e_d_svvz * ceff_map_dvv[y] * theMap[y*ncoords_base+p] * alpha_z; 
			
			dedr[3*p+0] += d_e_d_rux *  ceff_map_du[y] * theMap[y*ncoords_base+p];
			dedr[3*p+1] += d_e_d_ruy *  ceff_map_du[y] * theMap[y*ncoords_base+p];
			dedr[3*p+2] += d_e_d_ruz *  ceff_map_du[y] * theMap[y*ncoords_base+p];
			                                                                     
			dedr[3*p+0] += d_e_d_rvx *  ceff_map_dv[y] * theMap[y*ncoords_base+p];
			dedr[3*p+1] += d_e_d_rvy *  ceff_map_dv[y] * theMap[y*ncoords_base+p];
			dedr[3*p+2] += d_e_d_rvz *  ceff_map_dv[y] * theMap[y*ncoords_base+p];
			                                                                       
			dedr[3*p+0] += d_e_d_suux * ceff_map_duu[y] * theMap[y*ncoords_base+p]; 
			dedr[3*p+1] += d_e_d_suuy * ceff_map_duu[y] * theMap[y*ncoords_base+p]; 
			dedr[3*p+2] += d_e_d_suuz * ceff_map_duu[y] * theMap[y*ncoords_base+p]; 
			                                                                      
			dedr[3*p+0] += d_e_d_suvx * ceff_map_duv[y] * theMap[y*ncoords_base+p]; 
			dedr[3*p+1] += d_e_d_suvy * ceff_map_duv[y] * theMap[y*ncoords_base+p]; 
			dedr[3*p+2] += d_e_d_suvz * ceff_map_duv[y] * theMap[y*ncoords_base+p]; 
			                                                                      
			dedr[3*p+0] += d_e_d_svvx * ceff_map_dvv[y] * theMap[y*ncoords_base+p]; 
			dedr[3*p+1] += d_e_d_svvy * ceff_map_dvv[y] * theMap[y*ncoords_base+p]; 
			dedr[3*p+2] += d_e_d_svvz * ceff_map_dvv[y] * theMap[y*ncoords_base+p]; 
		}
		for( int p = 0; p < np; p++ )
		{
			gr[3*nv+0] += dedr[3*p+0] * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
			gr[3*nv+1] += dedr[3*p+1] * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
			gr[3*nv+2] += dedr[3*p+2] * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
		}
	}	
}
void surface::pgrad( double *r, double *gr, double *p_uv, double *pg )
{
	if( !setup_for_parallel )
		setupParallel(this,NULL,0,-1);

	double alpha_x = r[3*nv];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];
	// gradient of the particle energy wrt the vertices.

	for( int fx = 0; fx < par_info.nf; fx++ )
	{
		int f = par_info.faces[fx];

		if( f >= nf_faces )
			continue;

		int frm = f*nf_g_q_p;
		int t = theFormulas[frm].tri;
	
		for( int px = 0; px < theTriangles[t].np; px++ )
		{
			int pid = theTriangles[t].plist[px];
			double p_c0 = theTriangles[t].pc0[px];

//		 	NOTE 1: at some point I had a good reason to make this exception.	
//			if( fabs(p_c0 - theFormulas[frm].c0) < 1e-7 )
//				continue; 

			double p_area = theTriangles[t].pa[px];
	
			double u = p_uv[2*pid+0];
			double v = p_uv[2*pid+1];

			int *cp = theFormulas[frm].cp;
	
			double w = 1-u-v;
	
			double u2 = u*u;
			double u3 = u*u*u;
			double u4 = u*u*u*u;
			
			double v2 = v*v;
			double v3 = v*v*v;
			double v4 = v*v*v*v;
			
			double w2 = w*w;
			double w3 = w*w*w;
			double w4 = w*w*w*w;
	
			double n1 = (1.0/12.0)*(u4+2*u3*v); 
			double n2 = (1.0/12.0)*(u4+2*u3*w); 
			double n3 = (1.0/12.0)*(u4+2*u3*w+6*u3*v+6*u2*v*w+12*u2*v2+6*u*v2*w+6*u*v3+2*v3*w+v4); 
			double n4 = (1.0/12.0)*(6*u4+24*u3*w+24*u2*w2+8*u*w3+w4+24*u3*v+60*u2*v*w+36*u*v*w2+6*v*w3+24*u2*v2+36*u*v2*w+12*v2*w2+8*u*v3+6*v3*w+v4); 
			double n5 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+2*u3*v+6*u2*v*w+6*u*v*w2+2*v*w3); 
			double n6 = (1.0/12.0)*(2*u*v3+v4); 
			double n7 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+8*u3*v+36*u2*v*w+36*u*v*w2+8*v*w3+24*u2*v2+60*u*v2*w+24*v2*w2+24*u*v3+24*v3*w+6*v4); 
			double n8 = (1.0/12.0)*(u4+8*u3*w+24*u2*w2+24*u*w3+6*w4+6*u3*v+36*u2*v*w+60*u*v*w2+24*v*w3+12*u2*v2+36*u*v2*w+24*v2*w2+6*u*v3+8*v3*w+v4); 
			double n9 = (1.0/12.0)*(2*u*w3+w4); 
			double n10 = (1.0/12.0)*(2*v3*w+v4); 
			double n11 = (1.0/12.0)*(2*u*w3+w4+6*u*v*w2+6*v*w3+6*u*v2*w+12*v2*w2+2*u*v3+6*v3*w+v4); 
			double n12 = (1.0/12.0)*(w4+2*v*w3);
	
			double du_1 = Power(u,3)/3. + (Power(u,2)*v)/2.;
			double du_2 = Power(u,2)/2. - Power(u,3)/3. - (Power(u,2)*v)/2.; 	
			double du_3 = Power(u,2)/2. - Power(u,3)/3. + u*v - (Power(u,2)*v)/2. + Power(v,2)/2. - Power(v,3)/6.;
			double du_4 = 0.3333333333333333 + u - Power(u,2) - Power(u,3)/3. + v/2. - u*v - (Power(u,2)*v)/2. - Power(v,2) + Power(v,3)/3.;	
			double du_5 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. - v/2. + Power(u,2)*v + Power(v,2)/2. - Power(v,3)/6.;
			double du_6 = Power(v,3)/6.;	
			double du_7 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. + v/2. - 2*u*v + Power(u,2)*v - Power(v,2)/2. - Power(v,3)/6.;	
			double du_8 = -2*u + 2*Power(u,2) - Power(u,3)/3. - v + 2*u*v - (Power(u,2)*v)/2. + Power(v,2) - Power(v,3)/6.;	
			double du_9 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. + v/2. - (Power(u,2)*v)/2. - Power(v,2)/2. + Power(v,3)/6.;	
			double du_10 = -Power(v,3)/6.;
			double du_11 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. - v/2. + u*v - (Power(u,2)*v)/2. + Power(v,3)/3.;	
			double du_12 = 	-0.3333333333333333 + u - Power(u,2) + Power(u,3)/3. + v/2. - u*v + (Power(u,2)*v)/2. - Power(v,3)/6.;
	
			double dv_1 = Power(u,3)/6.;
			double dv_2 = -Power(u,3)/6.;
			double dv_3 = Power(u,2)/2. - Power(u,3)/6. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_4 = 0.16666666666666666 + u/2. - Power(u,2)/2. - Power(u,3)/6. - 2*u*v - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
			double dv_5 = -0.16666666666666666 - u/2. + Power(u,3)/3. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_6 = (u*Power(v,2))/2. + Power(v,3)/3.;
			double dv_7 = 0.3333333333333333 + u/2. - Power(u,2) + Power(u,3)/3. + v - u*v - Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_8 = -u + Power(u,2) - Power(u,3)/6. - 2*v + 2*u*v + 2*Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_9 = -0.3333333333333333 + u/2. - Power(u,3)/6. + v - u*v - Power(v,2) + (u*Power(v,2))/2. + Power(v,3)/3.;
			double dv_10 = Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_11 = 0.16666666666666666 - u/2. + Power(u,2)/2. - Power(u,3)/6. - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
			double dv_12 = -0.16666666666666666 + u/2. - Power(u,2)/2. + Power(u,3)/6. + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
	
			double d_uu_1 = Power(u,2) + u*v;
			double d_uu_2 = u - Power(u,2) - u*v;
			double d_uu_3 = u - Power(u,2) + v - u*v;
			double d_uu_4 = 1 - 2*u - Power(u,2) - v - u*v;
			double d_uu_5 = -2*u + 2*Power(u,2) + 2*u*v;
			double d_uu_6 = 0;
			double d_uu_7 = -2*u + 2*Power(u,2) - 2*v + 2*u*v;
			double d_uu_8 = -2 + 4*u - Power(u,2) + 2*v - u*v;
			double d_uu_9 = u - Power(u,2) - u*v;
			double d_uu_10 = 0;
			double d_uu_11 = u - Power(u,2) + v - u*v;
			double d_uu_12 = 1 - 2*u + Power(u,2) - v + u*v;
			
			double d_uv_1 = Power(u,2)/2.;
			double d_uv_2 = -Power(u,2)/2.;
			double d_uv_3 = u - Power(u,2)/2. + v - Power(v,2)/2.;
			double d_uv_4 = 0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
			double d_uv_5 = -0.5 + Power(u,2) + v - Power(v,2)/2.;
			double d_uv_6 = Power(v,2)/2.;
			double d_uv_7 = 0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
			double d_uv_8 = -1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
			double d_uv_9 = 0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
			double d_uv_10 = -Power(v,2)/2.;
			double d_uv_11 = -0.5 + u - Power(u,2)/2. + Power(v,2);
			double d_uv_12 = 0.5 - u + Power(u,2)/2. - Power(v,2)/2.;
			
			double d_vv_1 = 0;
			double d_vv_2 = 0;
			double d_vv_3 = u + v - u*v - Power(v,2);
			double d_vv_4 = -2*u - 2*v + 2*u*v + 2*Power(v,2);
			double d_vv_5 = u + v - u*v - Power(v,2);
			double d_vv_6 = u*v + Power(v,2);
			double d_vv_7 = 1 - u - 2*v - u*v - Power(v,2);
			double d_vv_8 = -2 + 2*u + 4*v - u*v - Power(v,2);
			double d_vv_9 = 1 - u - 2*v + u*v + Power(v,2);
			double d_vv_10 = v - u*v - Power(v,2);
			double d_vv_11 = -2*v + 2*u*v + 2*Power(v,2);
			double d_vv_12 = v - u*v - Power(v,2);
			

			double d_uuu_1 = 2*u + v;
			double d_uuu_2 = 1-2*u-v; //u - Power(u,2) - u*v;
			double d_uuu_3 = 1-2*u-v;//u - Power(u,2) + v - u*v;
			double d_uuu_4 = -2-2*u-v;//1 - 2*u - Power(u,2) - v - u*v;
			double d_uuu_5 = -2+4*u+2*v;//-2*u + 2*Power(u,2) + 2*u*v;
			double d_uuu_6 = 0;
			double d_uuu_7 = -2+4*u+2*v;//-2*u + 2*Power(u,2) - 2*v + 2*u*v;
			double d_uuu_8 = 4-2*u-v;//-2 + 4*u - Power(u,2) + 2*v - u*v;
			double d_uuu_9 = 1-2*u-v;//u - Power(u,2) - u*v;
			double d_uuu_10 = 0;
			double d_uuu_11 = 1-2*u-v;//u - Power(u,2) + v - u*v;
			double d_uuu_12 = -2+2*u+v;//1 - 2*u + Power(u,2) - v + u*v;
	
			double d_uuv_1 = u;//Power(u,2)/2.;
			double d_uuv_2 = -u;//-Power(u,2)/2.;
			double d_uuv_3 = 1-u;//u - Power(u,2)/2. + v - Power(v,2)/2.;
			double d_uuv_4 = -1-u;//0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
			double d_uuv_5 = 2*u;//-0.5 + Power(u,2) + v - Power(v,2)/2.;
			double d_uuv_6 = 0;//Power(v,2)/2.;
			double d_uuv_7 = -2+2*u;//0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
			double d_uuv_8 = 2-u;//-1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
			double d_uuv_9 = -u;//0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
			double d_uuv_10 = 0;//-Power(v,2)/2.;
			double d_uuv_11 = 1-u;//-0.5 + u - Power(u,2)/2. + Power(v,2);
			double d_uuv_12 = -1+u;//0.5 - u + Power(u,2)/2. - Power(v,2)/2.;

			double d_uvv_1 = 0;
			double d_uvv_2 = 0;
			double d_uvv_3 = 1-v;//u + v - u*v - Power(v,2);
			double d_uvv_4 = -2+2*v;//-2*u - 2*v + 2*u*v + 2*Power(v,2);
			double d_uvv_5 = 1-v;//u + v - u*v - Power(v,2);
			double d_uvv_6 = v;//u*v + Power(v,2);
			double d_uvv_7 = -1-v;//1 - u - 2*v - u*v - Power(v,2);
			double d_uvv_8 = 2-v;//-2 + 2*u + 4*v - u*v - Power(v,2);
			double d_uvv_9 = -1+v;//1 - u - 2*v + u*v + Power(v,2);
			double d_uvv_10 = -v;//v - u*v - Power(v,2);
			double d_uvv_11 = 2*v;//-2*v + 2*u*v + 2*Power(v,2);
			double d_uvv_12 = -v;//v - u*v - Power(v,2);
			
			double d_vvv_1 = 0;
			double d_vvv_2 = 0;
			double d_vvv_3 = 1-u-2*v;//u + v - u*v - Power(v,2);
			double d_vvv_4 = -2+2*u+4*v;//-2*u - 2*v + 2*u*v + 2*Power(v,2);
			double d_vvv_5 = 1-u-2*v;//u + v - u*v - Power(v,2);
			double d_vvv_6 = u+2*v;//u*v + Power(v,2);
			double d_vvv_7 = -2-u-2*v;//1 - u - 2*v - u*v - Power(v,2);
			double d_vvv_8 = 4-u-2*v;//-2 + 2*u + 4*v - u*v - Power(v,2);
			double d_vvv_9 = -2+u+2*v;//1 - u - 2*v + u*v + Power(v,2);
			double d_vvv_10 = 1-u-2*v;//v - u*v - Power(v,2);
			double d_vvv_11 = -2+2*u+4*v;//-2*v + 2*u*v + 2*Power(v,2);
			double d_vvv_12 = 1-u-2*v;//v - u*v - Power(v,2);
	
			double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
			double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
			double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
			
			double ceff_map_duu[12] = { d_uu_8, d_uu_7, d_uu_4, d_uu_5, d_uu_9, d_uu_12, d_uu_11, d_uu_10, d_uu_6, d_uu_3, d_uu_1, d_uu_2 };
			double ceff_map_duv[12] = { d_uv_8, d_uv_7, d_uv_4, d_uv_5, d_uv_9, d_uv_12, d_uv_11, d_uv_10, d_uv_6, d_uv_3, d_uv_1, d_uv_2 };
			double ceff_map_dvv[12] = { d_vv_8, d_vv_7, d_vv_4, d_vv_5, d_vv_9, d_vv_12, d_vv_11, d_vv_10, d_vv_6, d_vv_3, d_vv_1, d_vv_2 };
			
			double ceff_map_duuu[12] = { d_uuu_8, d_uuu_7, d_uuu_4, d_uuu_5, d_uuu_9, d_uuu_12, d_uuu_11, d_uuu_10, d_uuu_6, d_uuu_3, d_uuu_1, d_uuu_2 };
			double ceff_map_duuv[12] = { d_uuv_8, d_uuv_7, d_uuv_4, d_uuv_5, d_uuv_9, d_uuv_12, d_uuv_11, d_uuv_10, d_uuv_6, d_uuv_3, d_uuv_1, d_uuv_2 };
			double ceff_map_duvv[12] = { d_uvv_8, d_uvv_7, d_uvv_4, d_uvv_5, d_uvv_9, d_uvv_12, d_uvv_11, d_uvv_10, d_uvv_6, d_uvv_3, d_uvv_1, d_uvv_2 };
			double ceff_map_dvvv[12] = { d_vvv_8, d_vvv_7, d_vvv_4, d_vvv_5, d_vvv_9, d_vvv_12, d_vvv_11, d_vvv_10, d_vvv_6, d_vvv_3, d_vvv_1, d_vvv_2 };
					
			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};

			double tSuuu[3] = {0,0,0};
			double tSuuv[3] = {0,0,0};
			double tSuvv[3] = {0,0,0};
			double tSvvv[3] = {0,0,0};
	
			double nrm[3]={0,0,0}; 
	
			int np = theFormulas[frm].ncoor;
			double dedr[3*np];
			memset( dedr, 0, sizeof(double) * 3 *np );
	
			for( int p = 0; p < np; p++ )
			{
	
				R[0] += ceff_map[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				R[1] += ceff_map[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				R[2] += ceff_map[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				Ru[0] += ceff_map_du[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				Ru[1] += ceff_map_du[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				Ru[2] += ceff_map_du[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				Rv[0] += ceff_map_dv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				Rv[1] += ceff_map_dv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				Rv[2] += ceff_map_dv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuu[0] += ceff_map_duu[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuu[1] += ceff_map_duu[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuu[2] += ceff_map_duu[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuv[0] += ceff_map_duv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuv[1] += ceff_map_duv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuv[2] += ceff_map_duv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSvv[0] += ceff_map_dvv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSvv[1] += ceff_map_dvv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSvv[2] += ceff_map_dvv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuuu[0] += ceff_map_duuu[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuuu[1] += ceff_map_duuu[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuuu[2] += ceff_map_duuu[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuuv[0] += ceff_map_duuv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuuv[1] += ceff_map_duuv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuuv[2] += ceff_map_duuv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuvv[0] += ceff_map_duvv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuvv[1] += ceff_map_duvv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuvv[2] += ceff_map_duvv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSvvv[0] += ceff_map_dvvv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSvvv[1] += ceff_map_dvvv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSvvv[2] += ceff_map_dvvv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			}

			cross( Ru, Rv, nrm );
			normalize(nrm);
	
	
			double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
			double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
			double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];
	
			double g = sqrt(RuRu*RvRv-RuRv*RuRv);
	
			double e1,e2;
	      
			double nsuu = tSuu[0] * nrm[0] + tSuu[1] * nrm[1] + tSuu[2] * nrm[2];
			double nsuv = tSuv[0] * nrm[0] + tSuv[1] * nrm[1] + tSuv[2] * nrm[2];
			double nsvv = tSvv[0] * nrm[0] + tSvv[1] * nrm[1] + tSvv[2] * nrm[2];
	
			double Stot = (nsuu * RvRv + nsvv * RuRu -2*nsuv*RuRv)/(g*g);
	
			double Sop[4] = { 1.0/(g*g) * (nsuu * RvRv  - nsuv * RuRv), (1.0/(g*g)) * ( nsuv*RvRv-nsvv*RuRv),
					  1.0/(g*g) * (nsuv * RuRu  - nsuu * RuRv), (1.0/(g*g)) * ( nsvv*RuRu-nsuv*RuRv) };
	
			double a = Sop[0];		
			double b = Sop[1];
			double c = Sop[2];
			double d = Sop[3];
	
			double c0 = theFormulas[frm].c0;
			e1 = -0.5*(a+d-sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
			e2 = -0.5*(a+d+sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
		

double d_nrmz_d_rux=0,d_nrmz_d_ruy=0,d_nrmz_d_ruz=0,d_nrmz_d_rvx=0,d_nrmz_d_rvy=0,d_nrmz_d_rvz=0,d_nrmy_d_rux=0,d_nrmy_d_ruy=0,d_nrmy_d_ruz=0,d_nrmy_d_rvx=0,d_nrmy_d_rvy=0,d_nrmy_d_rvz=0,d_nrmx_d_rux=0,d_nrmx_d_ruy=0,d_nrmx_d_ruz=0,d_nrmx_d_rvx=0,d_nrmx_d_rvy=0,d_nrmx_d_rvz=0,d_RuRu_d_rux=0,d_RuRu_d_ruy=0,d_RuRu_d_ruz=0,d_nsvv_d_rux=0,d_nsvv_d_ruy=0,d_nsvv_d_ruz=0,d_nsvv_d_rvx=0,d_nsvv_d_rvy=0,d_nsvv_d_rvz=0,d_nsvv_d_svvx=0,d_nsvv_d_svvy=0,d_nsvv_d_svvz=0,d_RuRv_d_rux=0,d_RuRv_d_ruy=0,d_RuRv_d_ruz=0,d_RuRv_d_rvx=0,d_RuRv_d_rvy=0,d_RuRv_d_rvz=0,d_nsuv_d_rux=0,d_nsuv_d_ruy=0,d_nsuv_d_ruz=0,d_nsuv_d_rvx=0,d_nsuv_d_rvy=0,d_nsuv_d_rvz=0,d_nsuv_d_suvx=0,d_nsuv_d_suvy=0,d_nsuv_d_suvz=0,d_RvRv_d_rvx=0,d_RvRv_d_rvy=0,d_RvRv_d_rvz=0,d_nsuu_d_rux=0,d_nsuu_d_ruy=0,d_nsuu_d_ruz=0,d_nsuu_d_rvx=0,d_nsuu_d_rvy=0,d_nsuu_d_rvz=0,d_nsuu_d_suux=0,d_nsuu_d_suuy=0,d_nsuu_d_suuz=0,d_g_d_rux=0,d_g_d_ruy=0,d_g_d_ruz=0,d_g_d_rvx=0,d_g_d_rvy=0,d_g_d_rvz=0,d_e_d_rux=0,d_e_d_ruy=0,d_e_d_ruz=0,d_e_d_rvx=0,d_e_d_rvy=0,d_e_d_rvz=0,d_e_d_suux=0,d_e_d_suuy=0,d_e_d_suuz=0,d_e_d_suvx=0,d_e_d_suvy=0,d_e_d_suvz=0,d_e_d_svvx=0,d_e_d_svvy=0,d_e_d_svvz=0,d_c2_d_rux=0,d_c2_d_ruy=0,d_c2_d_ruz=0,d_c2_d_rvx=0,d_c2_d_rvy=0,d_c2_d_rvz=0,d_c2_d_suux=0,d_c2_d_suuy=0,d_c2_d_suuz=0,d_c2_d_suvx=0,d_c2_d_suvy=0,d_c2_d_suvz=0,d_c2_d_svvx=0,d_c2_d_svvy=0,d_c2_d_svvz=0,d_d_d_rux=0,d_d_d_ruy=0,d_d_d_ruz=0,d_d_d_rvx=0,d_d_d_rvy=0,d_d_d_rvz=0,d_d_d_suvx=0,d_d_d_suvy=0,d_d_d_suvz=0,d_d_d_svvx=0,d_d_d_svvy=0,d_d_d_svvz=0,d_c_d_rux=0,d_c_d_ruy=0,d_c_d_ruz=0,d_c_d_rvx=0,d_c_d_rvy=0,d_c_d_rvz=0,d_c_d_suux=0,d_c_d_suuy=0,d_c_d_suuz=0,d_c_d_suvx=0,d_c_d_suvy=0,d_c_d_suvz=0,d_b_d_rux=0,d_b_d_ruy=0,d_b_d_ruz=0,d_b_d_rvx=0,d_b_d_rvy=0,d_b_d_rvz=0,d_b_d_suvx=0,d_b_d_suvy=0,d_b_d_suvz=0,d_b_d_svvx=0,d_b_d_svvy=0,d_b_d_svvz=0,d_a_d_rux=0,d_a_d_ruy=0,d_a_d_ruz=0,d_a_d_rvx=0,d_a_d_rvy=0,d_a_d_rvz=0,d_a_d_suux=0,d_a_d_suuy=0,d_a_d_suuz=0,d_a_d_suvx=0,d_a_d_suvy=0,d_a_d_suvz=0,d_c1_d_rux=0,d_c1_d_ruy=0,d_c1_d_ruz=0,d_c1_d_rvx=0,d_c1_d_rvy=0,d_c1_d_rvz=0,d_c1_d_suux=0,d_c1_d_suuy=0,d_c1_d_suuz=0,d_c1_d_suvx=0,d_c1_d_suvy=0,d_c1_d_suvz=0,d_c1_d_svvx=0,d_c1_d_svvy=0,d_c1_d_svvz=0,d_g_d_RvRv=0,d_g_d_RuRv=0,d_g_d_RuRu=0,d_nsvv_d_nrmz=0,d_nsvv_d_nrmy=0,d_nsvv_d_nrmx=0,d_nsuv_d_nrmz=0,d_nsuv_d_nrmy=0,d_nsuv_d_nrmx=0,d_nsuu_d_nrmz=0,d_nsuu_d_nrmy=0,d_nsuu_d_nrmx=0,d_d_d_RuRv=0,d_d_d_nsuv=0,d_d_d_RuRu=0,d_d_d_nsvv=0,d_d_d_g=0,d_c_d_RuRv=0,d_c_d_nsuv=0,d_c_d_RuRu=0,d_c_d_nsuu=0,d_c_d_g=0,d_b_d_RuRv=0,d_b_d_nsuv=0,d_b_d_RvRv=0,d_b_d_nsvv=0,d_b_d_g=0,d_a_d_RuRv=0,d_a_d_nsuv=0,d_a_d_RvRv=0,d_a_d_nsuu=0,d_a_d_g=0,d_e_d_c2=0,d_e_d_c1=0,d_e_d_g=0,d_c2_d_d=0,d_c2_d_c=0,d_c2_d_b=0,d_c2_d_a=0,d_c1_d_d=0,d_c1_d_c=0,d_c1_d_b=0,d_c1_d_a=0,junk;
			// the basic variables.





			d_g_d_RuRu = RvRv/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
			d_g_d_RuRv = -(RuRv/Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
			d_g_d_RvRv = RuRu/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  

			d_RuRu_d_rux = 2*Ru[0];
			d_RuRu_d_ruy = 2*Ru[1];
			d_RuRu_d_ruz = 2*Ru[2];
			
			d_RuRv_d_rux = Rv[0];
			d_RuRv_d_ruy = Rv[1];
			d_RuRv_d_ruz = Rv[2];
			
			d_RuRv_d_rvx = Ru[0];
			d_RuRv_d_rvy = Ru[1];
			d_RuRv_d_rvz = Ru[2];
			
			d_RvRv_d_rvx = 2*Rv[0];
			d_RvRv_d_rvy = 2*Rv[1];
			d_RvRv_d_rvz = 2*Rv[2];

			// intermediates.
			d_a_d_g = (-2*(-(nsuv*RuRv) + nsuu*RvRv))/Power(g,3);
			d_a_d_nsuu = RvRv/Power(g,2);
			d_a_d_RvRv = nsuu/Power(g,2);
			d_a_d_nsuv = -(RuRv/Power(g,2));
			d_a_d_RuRv = -(nsuv/Power(g,2));

			d_b_d_g = (-2*(-(nsvv*RuRv) + nsuv*RvRv))/Power(g,3);
			d_b_d_nsvv = -(RuRv/Power(g,2));
			d_b_d_RvRv = nsuv/Power(g,2);
			d_b_d_nsuv = RvRv/Power(g,2);
			d_b_d_RuRv = -(nsvv/Power(g,2));

			d_c_d_g = (-2*(nsuv*RuRu - nsuu*RuRv))/Power(g,3);
			d_c_d_nsuu = -(RuRv/Power(g,2));
			d_c_d_RuRu = nsuv/Power(g,2);
			d_c_d_nsuv = RuRu/Power(g,2);
			d_c_d_RuRv = -(nsuu/Power(g,2));

			d_d_d_g = (-2*(nsvv*RuRu - nsuv*RuRv))/Power(g,3);
			d_d_d_nsvv = RuRu/Power(g,2);
			d_d_d_RuRu = nsvv/Power(g,2);
			d_d_d_nsuv = -(RuRv/Power(g,2));
			d_d_d_RuRv = -(nsuv/Power(g,2));

			d_c1_d_a =-0.5*(1 - (2*a - 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
			d_c1_d_b =-(-1.*c)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2)) ;
			d_c1_d_c =-(-1.*b)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2)) ;
			d_c1_d_d =-0.5*(1 - (-2*a + 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
		
			d_c2_d_a = -0.5*(1 + (2*a - 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
			d_c2_d_b = -(1.*c)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2));
			d_c2_d_c = -(1.*b)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2));
			d_c2_d_d = -0.5*(1 + (-2*a + 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));


			// this is the term we need to gradient:
			// en += 0.5*kc * p_area * ( c1 + c2 - p_c0 ) * ( c1 + c2 - p_c0 );
			d_e_d_c1 = kc * p_area * ( e1 + e2 - p_c0);
			d_e_d_c2 = kc * p_area * ( e1 + e2 - p_c0);
			d_e_d_g = 0;

			// that's it.

			double fac = Power(Power(Ru[1]*Rv[0] - Ru[0]*Rv[1],2) + Power(Ru[2]*Rv[0] - Ru[0]*Rv[2],2) + Power(Ru[2]*Rv[1] - Ru[1]*Rv[2],2),1.5);

			d_nrmx_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac; 
			d_nrmx_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2)))/fac;
			d_nrmx_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac;
			
			d_nrmy_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2))))/fac;
			d_nrmy_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac; 
			d_nrmy_d_ruz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2])) + Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac;
	
			d_nrmz_d_rux = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2])/fac;
			d_nrmz_d_ruy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;
			d_nrmz_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2))) + (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;

			d_nrmx_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac; 
			d_nrmx_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2])))/fac;
			d_nrmx_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac;
			
			d_nrmy_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
			d_nrmy_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2])))/fac;
			d_nrmy_d_rvz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[0]*Ru[1]*Rv[0] - Power(Ru[0],2)*Rv[1] + Ru[2]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2])))/fac;
	
			d_nrmz_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;
			d_nrmz_d_rvy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[0]*Ru[2]*Rv[0] - Power(Ru[0],2)*Rv[2] + Ru[1]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
			d_nrmz_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;

			double d_nrmx_d_u = d_nrmx_d_rux * tSuu[0] + d_nrmx_d_ruy * tSuu[1] + d_nrmx_d_ruz * tSuu[2] +
					    d_nrmx_d_rvx * tSuv[0] + d_nrmx_d_rvy * tSuv[1] + d_nrmx_d_rvz * tSuv[2];
			double d_nrmx_d_v = d_nrmx_d_rux * tSuv[0] + d_nrmx_d_ruy * tSuv[1] + d_nrmx_d_ruz * tSuv[2] +
					    d_nrmx_d_rvx * tSvv[0] + d_nrmx_d_rvy * tSvv[1] + d_nrmx_d_rvz * tSvv[2];
			double d_nrmy_d_u = d_nrmy_d_rux * tSuu[0] + d_nrmy_d_ruy * tSuu[1] + d_nrmy_d_ruz * tSuu[2] +
					    d_nrmy_d_rvx * tSuv[0] + d_nrmy_d_rvy * tSuv[1] + d_nrmy_d_rvz * tSuv[2];
			double d_nrmy_d_v = d_nrmy_d_rux * tSuv[0] + d_nrmy_d_ruy * tSuv[1] + d_nrmy_d_ruz * tSuv[2] +
					    d_nrmy_d_rvx * tSvv[0] + d_nrmy_d_rvy * tSvv[1] + d_nrmy_d_rvz * tSvv[2];
			double d_nrmz_d_u = d_nrmz_d_rux * tSuu[0] + d_nrmz_d_ruy * tSuu[1] + d_nrmz_d_ruz * tSuu[2] +
					    d_nrmz_d_rvx * tSuv[0] + d_nrmz_d_rvy * tSuv[1] + d_nrmz_d_rvz * tSuv[2];
			double d_nrmz_d_v = d_nrmz_d_rux * tSuv[0] + d_nrmz_d_ruy * tSuv[1] + d_nrmz_d_ruz * tSuv[2] +
					    d_nrmz_d_rvx * tSvv[0] + d_nrmz_d_rvy * tSvv[1] + d_nrmz_d_rvz * tSvv[2];

			double d_nsuu_d_u = tSuuu[0] * nrm[0] + tSuuu[1] * nrm[1] + tSuuu[2] * nrm[2];
			       d_nsuu_d_u += tSuu[0] * d_nrmx_d_u + tSuu[1] * d_nrmy_d_u + tSuu[2] * d_nrmz_d_u;
			double d_nsuu_d_v = tSuuv[0] * nrm[0] + tSuuv[1] * nrm[1] + tSuuv[2] * nrm[2];
			       d_nsuu_d_v += tSuu[0] * d_nrmx_d_v + tSuu[1] * d_nrmy_d_v + tSuu[2] * d_nrmz_d_v;
			
			double d_nsuv_d_u = tSuuv[0] * nrm[0] + tSuuv[1] * nrm[1] + tSuuv[2] * nrm[2];
			       d_nsuv_d_u += tSuv[0] * d_nrmx_d_u + tSuv[1] * d_nrmy_d_u + tSuv[2] * d_nrmz_d_u;
			double d_nsuv_d_v = tSuvv[0] * nrm[0] + tSuvv[1] * nrm[1] + tSuvv[2] * nrm[2];
			       d_nsuv_d_v += tSuv[0] * d_nrmx_d_v + tSuv[1] * d_nrmy_d_v + tSuv[2] * d_nrmz_d_v;

			double d_nsvv_d_u = tSuvv[0] * nrm[0] + tSuvv[1] * nrm[1] + tSuvv[2] * nrm[2];
			       d_nsvv_d_u += tSvv[0] * d_nrmx_d_u + tSvv[1] * d_nrmy_d_u + tSvv[2] * d_nrmz_d_u;
			double d_nsvv_d_v = tSvvv[0] * nrm[0] + tSvvv[1] * nrm[1] + tSvvv[2] * nrm[2];
			       d_nsvv_d_v += tSvv[0] * d_nrmx_d_v + tSvv[1] * d_nrmy_d_v + tSvv[2] * d_nrmz_d_v;

			double d_RuRu_d_u = 2*Ru[0] * tSuu[0] + 2*Ru[1] * tSuu[1] + 2 * Ru[2] * tSuu[2];
			double d_RuRu_d_v = 2*Ru[0] * tSuv[0] + 2*Ru[1] * tSuv[1] + 2 * Ru[2] * tSuv[2];

			double d_RuRv_d_u =  Ru[0] * tSuv[0] + Ru[1] * tSuv[1] + Ru[2] * tSuv[2];
			       d_RuRv_d_u += Rv[0] * tSuu[0] + Rv[1] * tSuu[1] + Rv[2] * tSuu[2];
			
			double d_RuRv_d_v =  Ru[0] * tSvv[0] + Ru[1] * tSvv[1] + Ru[2] * tSvv[2];
			       d_RuRv_d_v += Rv[0] * tSuv[0] + Rv[1] * tSuv[1] + Rv[2] * tSuv[2];

			double d_RvRv_d_u = 2*Rv[0] * tSuv[0] + 2*Rv[1] * tSuv[1] + 2 * Rv[2] * tSuv[2];
			double d_RvRv_d_v = 2*Rv[0] * tSvv[0] + 2*Rv[1] * tSvv[1] + 2 * Rv[2] * tSvv[2];

			double d_g_d_u = 0;
			double d_g_d_v = 0;

			d_g_d_u += d_g_d_RuRu * 2* Ru[0] * tSuu[0] + d_g_d_RuRu * 2* Ru[1] * tSuu[1] + d_g_d_RuRu * 2* Ru[2] * tSuu[2];
			d_g_d_v += d_g_d_RuRu * 2* Ru[0] * tSuv[0] + d_g_d_RuRu * 2* Ru[1] * tSuv[1] + d_g_d_RuRu * 2* Ru[2] * tSuv[2];

			// d_ru dru_du
			d_g_d_u += d_g_d_RuRv * Rv[0] * tSuu[0] + d_g_d_RuRv * Rv[1] * tSuu[1] + d_g_d_RuRv * Rv[2] * tSuu[2];
			// d_rv_drv_du
			d_g_d_u += d_g_d_RuRv * Ru[0] * tSuv[0] + d_g_d_RuRv * Ru[1] * tSuv[1] + d_g_d_RuRv * Ru[2] * tSuv[2];	
			// d_ru dru_dv
			d_g_d_v += d_g_d_RuRv * Rv[0] * tSuv[0] + d_g_d_RuRv * Rv[1] * tSuv[1] + d_g_d_RuRv * Rv[2] * tSuv[2];
			// d_rv drv_dv
			d_g_d_v += d_g_d_RuRv * Ru[0] * tSvv[0] + d_g_d_RuRv * Ru[1] * tSvv[1] + d_g_d_RuRv * Ru[2] * tSvv[2];
			
			d_g_d_u += d_g_d_RvRv * 2* Rv[0] * tSuv[0] + d_g_d_RvRv * 2* Rv[1] * tSuv[1] + d_g_d_RvRv * 2* Rv[2] * tSuv[2];
			d_g_d_v += d_g_d_RvRv * 2* Rv[0] * tSvv[0] + d_g_d_RvRv * 2* Rv[1] * tSvv[1] + d_g_d_RvRv * 2* Rv[2] * tSvv[2];

//			d_a_d_g = (-2*(-(nsuv*RuRv) + nsuu*RvRv))/Power(g,3);
//			d_a_d_nsuu = RvRv/Power(g,2);
//			d_a_d_RvRv = nsuu/Power(g,2);
//			d_a_d_nsuv = -(RuRv/Power(g,2));
//			d_a_d_RuRv = -(nsuv/Power(g,2));
			
			double d_a_d_u = d_a_d_nsuu * d_nsuu_d_u + d_a_d_nsuv * d_nsuv_d_u + d_a_d_RuRv * d_RuRv_d_u + d_a_d_RvRv * d_RvRv_d_u + d_a_d_g * d_g_d_u;
			double d_a_d_v = d_a_d_nsuu * d_nsuu_d_v + d_a_d_nsuv * d_nsuv_d_v + d_a_d_RuRv * d_RuRv_d_v + d_a_d_RvRv * d_RvRv_d_v + d_a_d_g * d_g_d_v;
			double d_b_d_u = d_b_d_nsvv * d_nsvv_d_u + d_b_d_nsuv * d_nsuv_d_u + d_b_d_RuRv * d_RuRv_d_u + d_b_d_RvRv * d_RvRv_d_u + d_b_d_g * d_g_d_u;
			double d_b_d_v = d_b_d_nsvv * d_nsvv_d_v + d_b_d_nsuv * d_nsuv_d_v + d_b_d_RuRv * d_RuRv_d_v + d_b_d_RvRv * d_RvRv_d_v + d_b_d_g * d_g_d_v;
			double d_c_d_u = d_c_d_nsuu * d_nsuu_d_u + d_c_d_nsuv * d_nsuv_d_u + d_c_d_RuRv * d_RuRv_d_u + d_c_d_RuRu * d_RuRu_d_u + d_c_d_g * d_g_d_u;
			double d_c_d_v = d_c_d_nsuu * d_nsuu_d_v + d_c_d_nsuv * d_nsuv_d_v + d_c_d_RuRv * d_RuRv_d_v + d_c_d_RuRu * d_RuRu_d_v + d_c_d_g * d_g_d_v;
			double d_d_d_u = d_d_d_nsvv * d_nsvv_d_u + d_d_d_nsuv * d_nsuv_d_u + d_d_d_RuRv * d_RuRv_d_u + d_d_d_RuRu * d_RuRu_d_u + d_d_d_g * d_g_d_u;
			double d_d_d_v = d_d_d_nsvv * d_nsvv_d_v + d_d_d_nsuv * d_nsuv_d_v + d_d_d_RuRv * d_RuRv_d_v + d_d_d_RuRu * d_RuRu_d_v + d_d_d_g * d_g_d_v;

			double d_c1_d_u = d_c1_d_a * d_a_d_u + d_c1_d_b * d_b_d_u + d_c1_d_c * d_c_d_u + d_c1_d_d * d_d_d_u;
			double d_c2_d_u = d_c2_d_a * d_a_d_u + d_c2_d_b * d_b_d_u + d_c2_d_c * d_c_d_u + d_c2_d_d * d_d_d_u;
	
			double d_c1_d_v = d_c1_d_a * d_a_d_v + d_c1_d_b * d_b_d_v + d_c1_d_c * d_c_d_v + d_c1_d_d * d_d_d_v;
			double d_c2_d_v = d_c2_d_a * d_a_d_v + d_c2_d_b * d_b_d_v + d_c2_d_c * d_c_d_v + d_c2_d_d * d_d_d_v;

			double d_e_d_u = d_e_d_c1 * d_c1_d_u + d_e_d_c2 * d_c2_d_u;
			double d_e_d_v = d_e_d_c1 * d_c1_d_v + d_e_d_c2 * d_c2_d_v;

			pg[2*pid+0] += d_e_d_u;
			pg[2*pid+1] += d_e_d_v;
	
//			printf("u %le v %le a: %.14le da_du: %le nrmx: %.14le dnrmx_du %.14le\n", u, v, b, d_b_d_u, nrm[0], d_nrmy_d_u );

			double e = 0.5* kc * p_area * ( e1 + e2 - p_c0 ) * ( e1 + e2 - p_c0 );
//			printf("e %le at uv: %le %le  pg: %le %le\n", e, p_uv[2*pid+0], p_uv[2*pid+1], d_e_d_u, d_e_d_v );
			
			d_nsuu_d_nrmx = tSuu[0];
			d_nsuu_d_nrmy = tSuu[1];
			d_nsuu_d_nrmz = tSuu[2];
			
			d_nsuv_d_nrmx = tSuv[0];
			d_nsuv_d_nrmy = tSuv[1];
			d_nsuv_d_nrmz = tSuv[2];

			d_nsvv_d_nrmx = tSvv[0];
			d_nsvv_d_nrmy = tSvv[1];
			d_nsvv_d_nrmz = tSvv[2];

			d_nsuu_d_suux = nrm[0];
			d_nsuu_d_suuy = nrm[1];
			d_nsuu_d_suuz = nrm[2];
			
			d_nsuv_d_suvx = nrm[0];
			d_nsuv_d_suvy = nrm[1];
			d_nsuv_d_suvz = nrm[2];

			d_nsvv_d_svvx = nrm[0];
			d_nsvv_d_svvy = nrm[1];
			d_nsvv_d_svvz = nrm[2];



	d_nsvv_d_rux += d_nsvv_d_nrmz * d_nrmz_d_rux;
	d_nsvv_d_ruy += d_nsvv_d_nrmz * d_nrmz_d_ruy;
	d_nsvv_d_ruz += d_nsvv_d_nrmz * d_nrmz_d_ruz;
	d_nsvv_d_rvx += d_nsvv_d_nrmz * d_nrmz_d_rvx;
	d_nsvv_d_rvy += d_nsvv_d_nrmz * d_nrmz_d_rvy;
	d_nsvv_d_rvz += d_nsvv_d_nrmz * d_nrmz_d_rvz;
	d_nsvv_d_rux += d_nsvv_d_nrmy * d_nrmy_d_rux;
	d_nsvv_d_ruy += d_nsvv_d_nrmy * d_nrmy_d_ruy;
	d_nsvv_d_ruz += d_nsvv_d_nrmy * d_nrmy_d_ruz;
	d_nsvv_d_rvx += d_nsvv_d_nrmy * d_nrmy_d_rvx;
	d_nsvv_d_rvy += d_nsvv_d_nrmy * d_nrmy_d_rvy;
	d_nsvv_d_rvz += d_nsvv_d_nrmy * d_nrmy_d_rvz;
	d_nsvv_d_rux += d_nsvv_d_nrmx * d_nrmx_d_rux;
	d_nsvv_d_ruy += d_nsvv_d_nrmx * d_nrmx_d_ruy;
	d_nsvv_d_ruz += d_nsvv_d_nrmx * d_nrmx_d_ruz;
	d_nsvv_d_rvx += d_nsvv_d_nrmx * d_nrmx_d_rvx;
	d_nsvv_d_rvy += d_nsvv_d_nrmx * d_nrmx_d_rvy;
	d_nsvv_d_rvz += d_nsvv_d_nrmx * d_nrmx_d_rvz;
	d_nsuv_d_rux += d_nsuv_d_nrmz * d_nrmz_d_rux;
	d_nsuv_d_ruy += d_nsuv_d_nrmz * d_nrmz_d_ruy;
	d_nsuv_d_ruz += d_nsuv_d_nrmz * d_nrmz_d_ruz;
	d_nsuv_d_rvx += d_nsuv_d_nrmz * d_nrmz_d_rvx;
	d_nsuv_d_rvy += d_nsuv_d_nrmz * d_nrmz_d_rvy;
	d_nsuv_d_rvz += d_nsuv_d_nrmz * d_nrmz_d_rvz;
	d_nsuv_d_rux += d_nsuv_d_nrmy * d_nrmy_d_rux;
	d_nsuv_d_ruy += d_nsuv_d_nrmy * d_nrmy_d_ruy;
	d_nsuv_d_ruz += d_nsuv_d_nrmy * d_nrmy_d_ruz;
	d_nsuv_d_rvx += d_nsuv_d_nrmy * d_nrmy_d_rvx;
	d_nsuv_d_rvy += d_nsuv_d_nrmy * d_nrmy_d_rvy;
	d_nsuv_d_rvz += d_nsuv_d_nrmy * d_nrmy_d_rvz;
	d_nsuv_d_rux += d_nsuv_d_nrmx * d_nrmx_d_rux;
	d_nsuv_d_ruy += d_nsuv_d_nrmx * d_nrmx_d_ruy;
	d_nsuv_d_ruz += d_nsuv_d_nrmx * d_nrmx_d_ruz;
	d_nsuv_d_rvx += d_nsuv_d_nrmx * d_nrmx_d_rvx;
	d_nsuv_d_rvy += d_nsuv_d_nrmx * d_nrmx_d_rvy;
	d_nsuv_d_rvz += d_nsuv_d_nrmx * d_nrmx_d_rvz;
	d_nsuu_d_rux += d_nsuu_d_nrmz * d_nrmz_d_rux;
	d_nsuu_d_ruy += d_nsuu_d_nrmz * d_nrmz_d_ruy;
	d_nsuu_d_ruz += d_nsuu_d_nrmz * d_nrmz_d_ruz;
	d_nsuu_d_rvx += d_nsuu_d_nrmz * d_nrmz_d_rvx;
	d_nsuu_d_rvy += d_nsuu_d_nrmz * d_nrmz_d_rvy;
	d_nsuu_d_rvz += d_nsuu_d_nrmz * d_nrmz_d_rvz;
	d_nsuu_d_rux += d_nsuu_d_nrmy * d_nrmy_d_rux;
	d_nsuu_d_ruy += d_nsuu_d_nrmy * d_nrmy_d_ruy;
	d_nsuu_d_ruz += d_nsuu_d_nrmy * d_nrmy_d_ruz;
	d_nsuu_d_rvx += d_nsuu_d_nrmy * d_nrmy_d_rvx;
	d_nsuu_d_rvy += d_nsuu_d_nrmy * d_nrmy_d_rvy;
	d_nsuu_d_rvz += d_nsuu_d_nrmy * d_nrmy_d_rvz;
	d_nsuu_d_rux += d_nsuu_d_nrmx * d_nrmx_d_rux;
	d_nsuu_d_ruy += d_nsuu_d_nrmx * d_nrmx_d_ruy;
	d_nsuu_d_ruz += d_nsuu_d_nrmx * d_nrmx_d_ruz;
	d_nsuu_d_rvx += d_nsuu_d_nrmx * d_nrmx_d_rvx;
	d_nsuu_d_rvy += d_nsuu_d_nrmx * d_nrmx_d_rvy;
	d_nsuu_d_rvz += d_nsuu_d_nrmx * d_nrmx_d_rvz;
	d_g_d_rvx += d_g_d_RvRv * d_RvRv_d_rvx;
	d_g_d_rvy += d_g_d_RvRv * d_RvRv_d_rvy;
	d_g_d_rvz += d_g_d_RvRv * d_RvRv_d_rvz;
	d_g_d_rux += d_g_d_RuRv * d_RuRv_d_rux;
	d_g_d_ruy += d_g_d_RuRv * d_RuRv_d_ruy;
	d_g_d_ruz += d_g_d_RuRv * d_RuRv_d_ruz;
	d_g_d_rvx += d_g_d_RuRv * d_RuRv_d_rvx;
	d_g_d_rvy += d_g_d_RuRv * d_RuRv_d_rvy;
	d_g_d_rvz += d_g_d_RuRv * d_RuRv_d_rvz;
	d_g_d_rux += d_g_d_RuRu * d_RuRu_d_rux;
	d_g_d_ruy += d_g_d_RuRu * d_RuRu_d_ruy;
	d_g_d_ruz += d_g_d_RuRu * d_RuRu_d_ruz;
	d_d_d_rux += d_d_d_RuRv * d_RuRv_d_rux;
	d_d_d_ruy += d_d_d_RuRv * d_RuRv_d_ruy;
	d_d_d_ruz += d_d_d_RuRv * d_RuRv_d_ruz;
	d_d_d_rvx += d_d_d_RuRv * d_RuRv_d_rvx;
	d_d_d_rvy += d_d_d_RuRv * d_RuRv_d_rvy;
	d_d_d_rvz += d_d_d_RuRv * d_RuRv_d_rvz;
	d_d_d_rux += d_d_d_nsuv * d_nsuv_d_rux;
	d_d_d_ruy += d_d_d_nsuv * d_nsuv_d_ruy;
	d_d_d_ruz += d_d_d_nsuv * d_nsuv_d_ruz;
	d_d_d_rvx += d_d_d_nsuv * d_nsuv_d_rvx;
	d_d_d_rvy += d_d_d_nsuv * d_nsuv_d_rvy;
	d_d_d_rvz += d_d_d_nsuv * d_nsuv_d_rvz;
	d_d_d_suvx += d_d_d_nsuv * d_nsuv_d_suvx;
	d_d_d_suvy += d_d_d_nsuv * d_nsuv_d_suvy;
	d_d_d_suvz += d_d_d_nsuv * d_nsuv_d_suvz;
	d_d_d_rux += d_d_d_RuRu * d_RuRu_d_rux;
	d_d_d_ruy += d_d_d_RuRu * d_RuRu_d_ruy;
	d_d_d_ruz += d_d_d_RuRu * d_RuRu_d_ruz;
	d_d_d_rux += d_d_d_nsvv * d_nsvv_d_rux;
	d_d_d_ruy += d_d_d_nsvv * d_nsvv_d_ruy;
	d_d_d_ruz += d_d_d_nsvv * d_nsvv_d_ruz;
	d_d_d_rvx += d_d_d_nsvv * d_nsvv_d_rvx;
	d_d_d_rvy += d_d_d_nsvv * d_nsvv_d_rvy;
	d_d_d_rvz += d_d_d_nsvv * d_nsvv_d_rvz;
	d_d_d_svvx += d_d_d_nsvv * d_nsvv_d_svvx;
	d_d_d_svvy += d_d_d_nsvv * d_nsvv_d_svvy;
	d_d_d_svvz += d_d_d_nsvv * d_nsvv_d_svvz;
	d_d_d_rux += d_d_d_g * d_g_d_rux;
	d_d_d_ruy += d_d_d_g * d_g_d_ruy;
	d_d_d_ruz += d_d_d_g * d_g_d_ruz;
	d_d_d_rvx += d_d_d_g * d_g_d_rvx;
	d_d_d_rvy += d_d_d_g * d_g_d_rvy;
	d_d_d_rvz += d_d_d_g * d_g_d_rvz;
	d_c_d_rux += d_c_d_RuRv * d_RuRv_d_rux;
	d_c_d_ruy += d_c_d_RuRv * d_RuRv_d_ruy;
	d_c_d_ruz += d_c_d_RuRv * d_RuRv_d_ruz;
	d_c_d_rvx += d_c_d_RuRv * d_RuRv_d_rvx;
	d_c_d_rvy += d_c_d_RuRv * d_RuRv_d_rvy;
	d_c_d_rvz += d_c_d_RuRv * d_RuRv_d_rvz;
	d_c_d_rux += d_c_d_nsuv * d_nsuv_d_rux;
	d_c_d_ruy += d_c_d_nsuv * d_nsuv_d_ruy;
	d_c_d_ruz += d_c_d_nsuv * d_nsuv_d_ruz;
	d_c_d_rvx += d_c_d_nsuv * d_nsuv_d_rvx;
	d_c_d_rvy += d_c_d_nsuv * d_nsuv_d_rvy;
	d_c_d_rvz += d_c_d_nsuv * d_nsuv_d_rvz;
	d_c_d_suvx += d_c_d_nsuv * d_nsuv_d_suvx;
	d_c_d_suvy += d_c_d_nsuv * d_nsuv_d_suvy;
	d_c_d_suvz += d_c_d_nsuv * d_nsuv_d_suvz;
	d_c_d_rux += d_c_d_RuRu * d_RuRu_d_rux;
	d_c_d_ruy += d_c_d_RuRu * d_RuRu_d_ruy;
	d_c_d_ruz += d_c_d_RuRu * d_RuRu_d_ruz;
	d_c_d_rux += d_c_d_nsuu * d_nsuu_d_rux;
	d_c_d_ruy += d_c_d_nsuu * d_nsuu_d_ruy;
	d_c_d_ruz += d_c_d_nsuu * d_nsuu_d_ruz;
	d_c_d_rvx += d_c_d_nsuu * d_nsuu_d_rvx;
	d_c_d_rvy += d_c_d_nsuu * d_nsuu_d_rvy;
	d_c_d_rvz += d_c_d_nsuu * d_nsuu_d_rvz;
	d_c_d_suux += d_c_d_nsuu * d_nsuu_d_suux;
	d_c_d_suuy += d_c_d_nsuu * d_nsuu_d_suuy;
	d_c_d_suuz += d_c_d_nsuu * d_nsuu_d_suuz;
	d_c_d_rux += d_c_d_g * d_g_d_rux;
	d_c_d_ruy += d_c_d_g * d_g_d_ruy;
	d_c_d_ruz += d_c_d_g * d_g_d_ruz;
	d_c_d_rvx += d_c_d_g * d_g_d_rvx;
	d_c_d_rvy += d_c_d_g * d_g_d_rvy;
	d_c_d_rvz += d_c_d_g * d_g_d_rvz;
	d_b_d_rux += d_b_d_RuRv * d_RuRv_d_rux;
	d_b_d_ruy += d_b_d_RuRv * d_RuRv_d_ruy;
	d_b_d_ruz += d_b_d_RuRv * d_RuRv_d_ruz;
	d_b_d_rvx += d_b_d_RuRv * d_RuRv_d_rvx;
	d_b_d_rvy += d_b_d_RuRv * d_RuRv_d_rvy;
	d_b_d_rvz += d_b_d_RuRv * d_RuRv_d_rvz;
	d_b_d_rux += d_b_d_nsuv * d_nsuv_d_rux;
	d_b_d_ruy += d_b_d_nsuv * d_nsuv_d_ruy;
	d_b_d_ruz += d_b_d_nsuv * d_nsuv_d_ruz;
	d_b_d_rvx += d_b_d_nsuv * d_nsuv_d_rvx;
	d_b_d_rvy += d_b_d_nsuv * d_nsuv_d_rvy;
	d_b_d_rvz += d_b_d_nsuv * d_nsuv_d_rvz;
	d_b_d_suvx += d_b_d_nsuv * d_nsuv_d_suvx;
	d_b_d_suvy += d_b_d_nsuv * d_nsuv_d_suvy;
	d_b_d_suvz += d_b_d_nsuv * d_nsuv_d_suvz;
	d_b_d_rvx += d_b_d_RvRv * d_RvRv_d_rvx;
	d_b_d_rvy += d_b_d_RvRv * d_RvRv_d_rvy;
	d_b_d_rvz += d_b_d_RvRv * d_RvRv_d_rvz;
	d_b_d_rux += d_b_d_nsvv * d_nsvv_d_rux;
	d_b_d_ruy += d_b_d_nsvv * d_nsvv_d_ruy;
	d_b_d_ruz += d_b_d_nsvv * d_nsvv_d_ruz;
	d_b_d_rvx += d_b_d_nsvv * d_nsvv_d_rvx;
	d_b_d_rvy += d_b_d_nsvv * d_nsvv_d_rvy;
	d_b_d_rvz += d_b_d_nsvv * d_nsvv_d_rvz;
	d_b_d_svvx += d_b_d_nsvv * d_nsvv_d_svvx;
	d_b_d_svvy += d_b_d_nsvv * d_nsvv_d_svvy;
	d_b_d_svvz += d_b_d_nsvv * d_nsvv_d_svvz;
	d_b_d_rux += d_b_d_g * d_g_d_rux;
	d_b_d_ruy += d_b_d_g * d_g_d_ruy;
	d_b_d_ruz += d_b_d_g * d_g_d_ruz;
	d_b_d_rvx += d_b_d_g * d_g_d_rvx;
	d_b_d_rvy += d_b_d_g * d_g_d_rvy;
	d_b_d_rvz += d_b_d_g * d_g_d_rvz;
	d_a_d_rux += d_a_d_RuRv * d_RuRv_d_rux;
	d_a_d_ruy += d_a_d_RuRv * d_RuRv_d_ruy;
	d_a_d_ruz += d_a_d_RuRv * d_RuRv_d_ruz;
	d_a_d_rvx += d_a_d_RuRv * d_RuRv_d_rvx;
	d_a_d_rvy += d_a_d_RuRv * d_RuRv_d_rvy;
	d_a_d_rvz += d_a_d_RuRv * d_RuRv_d_rvz;
	d_a_d_rux += d_a_d_nsuv * d_nsuv_d_rux;
	d_a_d_ruy += d_a_d_nsuv * d_nsuv_d_ruy;
	d_a_d_ruz += d_a_d_nsuv * d_nsuv_d_ruz;
	d_a_d_rvx += d_a_d_nsuv * d_nsuv_d_rvx;
	d_a_d_rvy += d_a_d_nsuv * d_nsuv_d_rvy;
	d_a_d_rvz += d_a_d_nsuv * d_nsuv_d_rvz;
	d_a_d_suvx += d_a_d_nsuv * d_nsuv_d_suvx;
	d_a_d_suvy += d_a_d_nsuv * d_nsuv_d_suvy;
	d_a_d_suvz += d_a_d_nsuv * d_nsuv_d_suvz;
	d_a_d_rvx += d_a_d_RvRv * d_RvRv_d_rvx;
	d_a_d_rvy += d_a_d_RvRv * d_RvRv_d_rvy;
	d_a_d_rvz += d_a_d_RvRv * d_RvRv_d_rvz;
	d_a_d_rux += d_a_d_nsuu * d_nsuu_d_rux;
	d_a_d_ruy += d_a_d_nsuu * d_nsuu_d_ruy;
	d_a_d_ruz += d_a_d_nsuu * d_nsuu_d_ruz;
	d_a_d_rvx += d_a_d_nsuu * d_nsuu_d_rvx;
	d_a_d_rvy += d_a_d_nsuu * d_nsuu_d_rvy;
	d_a_d_rvz += d_a_d_nsuu * d_nsuu_d_rvz;
	d_a_d_suux += d_a_d_nsuu * d_nsuu_d_suux;
	d_a_d_suuy += d_a_d_nsuu * d_nsuu_d_suuy;
	d_a_d_suuz += d_a_d_nsuu * d_nsuu_d_suuz;
	d_a_d_rux += d_a_d_g * d_g_d_rux;
	d_a_d_ruy += d_a_d_g * d_g_d_ruy;
	d_a_d_ruz += d_a_d_g * d_g_d_ruz;
	d_a_d_rvx += d_a_d_g * d_g_d_rvx;
	d_a_d_rvy += d_a_d_g * d_g_d_rvy;
	d_a_d_rvz += d_a_d_g * d_g_d_rvz;
	d_c1_d_rux += d_c1_d_d * d_d_d_rux;
	d_c1_d_ruy += d_c1_d_d * d_d_d_ruy;
	d_c1_d_ruz += d_c1_d_d * d_d_d_ruz;
	d_c1_d_rvx += d_c1_d_d * d_d_d_rvx;
	d_c1_d_rvy += d_c1_d_d * d_d_d_rvy;
	d_c1_d_rvz += d_c1_d_d * d_d_d_rvz;
	d_c1_d_suvx += d_c1_d_d * d_d_d_suvx;
	d_c1_d_suvy += d_c1_d_d * d_d_d_suvy;
	d_c1_d_suvz += d_c1_d_d * d_d_d_suvz;
	d_c1_d_svvx += d_c1_d_d * d_d_d_svvx;
	d_c1_d_svvy += d_c1_d_d * d_d_d_svvy;
	d_c1_d_svvz += d_c1_d_d * d_d_d_svvz;
	d_c1_d_rux += d_c1_d_c * d_c_d_rux;
	d_c1_d_ruy += d_c1_d_c * d_c_d_ruy;
	d_c1_d_ruz += d_c1_d_c * d_c_d_ruz;
	d_c1_d_rvx += d_c1_d_c * d_c_d_rvx;
	d_c1_d_rvy += d_c1_d_c * d_c_d_rvy;
	d_c1_d_rvz += d_c1_d_c * d_c_d_rvz;
	d_c1_d_suux += d_c1_d_c * d_c_d_suux;
	d_c1_d_suuy += d_c1_d_c * d_c_d_suuy;
	d_c1_d_suuz += d_c1_d_c * d_c_d_suuz;
	d_c1_d_suvx += d_c1_d_c * d_c_d_suvx;
	d_c1_d_suvy += d_c1_d_c * d_c_d_suvy;
	d_c1_d_suvz += d_c1_d_c * d_c_d_suvz;
	d_c1_d_rux += d_c1_d_b * d_b_d_rux;
	d_c1_d_ruy += d_c1_d_b * d_b_d_ruy;
	d_c1_d_ruz += d_c1_d_b * d_b_d_ruz;
	d_c1_d_rvx += d_c1_d_b * d_b_d_rvx;
	d_c1_d_rvy += d_c1_d_b * d_b_d_rvy;
	d_c1_d_rvz += d_c1_d_b * d_b_d_rvz;
	d_c1_d_suvx += d_c1_d_b * d_b_d_suvx;
	d_c1_d_suvy += d_c1_d_b * d_b_d_suvy;
	d_c1_d_suvz += d_c1_d_b * d_b_d_suvz;
	d_c1_d_svvx += d_c1_d_b * d_b_d_svvx;
	d_c1_d_svvy += d_c1_d_b * d_b_d_svvy;
	d_c1_d_svvz += d_c1_d_b * d_b_d_svvz;
	d_c1_d_rux += d_c1_d_a * d_a_d_rux;
	d_c1_d_ruy += d_c1_d_a * d_a_d_ruy;
	d_c1_d_ruz += d_c1_d_a * d_a_d_ruz;
	d_c1_d_rvx += d_c1_d_a * d_a_d_rvx;
	d_c1_d_rvy += d_c1_d_a * d_a_d_rvy;
	d_c1_d_rvz += d_c1_d_a * d_a_d_rvz;
	d_c1_d_suux += d_c1_d_a * d_a_d_suux;
	d_c1_d_suuy += d_c1_d_a * d_a_d_suuy;
	d_c1_d_suuz += d_c1_d_a * d_a_d_suuz;
	d_c1_d_suvx += d_c1_d_a * d_a_d_suvx;
	d_c1_d_suvy += d_c1_d_a * d_a_d_suvy;
	d_c1_d_suvz += d_c1_d_a * d_a_d_suvz;
	d_c2_d_rux += d_c2_d_d * d_d_d_rux;
	d_c2_d_ruy += d_c2_d_d * d_d_d_ruy;
	d_c2_d_ruz += d_c2_d_d * d_d_d_ruz;
	d_c2_d_rvx += d_c2_d_d * d_d_d_rvx;
	d_c2_d_rvy += d_c2_d_d * d_d_d_rvy;
	d_c2_d_rvz += d_c2_d_d * d_d_d_rvz;
	d_c2_d_suvx += d_c2_d_d * d_d_d_suvx;
	d_c2_d_suvy += d_c2_d_d * d_d_d_suvy;
	d_c2_d_suvz += d_c2_d_d * d_d_d_suvz;
	d_c2_d_svvx += d_c2_d_d * d_d_d_svvx;
	d_c2_d_svvy += d_c2_d_d * d_d_d_svvy;
	d_c2_d_svvz += d_c2_d_d * d_d_d_svvz;
	d_c2_d_rux += d_c2_d_c * d_c_d_rux;
	d_c2_d_ruy += d_c2_d_c * d_c_d_ruy;
	d_c2_d_ruz += d_c2_d_c * d_c_d_ruz;
	d_c2_d_rvx += d_c2_d_c * d_c_d_rvx;
	d_c2_d_rvy += d_c2_d_c * d_c_d_rvy;
	d_c2_d_rvz += d_c2_d_c * d_c_d_rvz;
	d_c2_d_suux += d_c2_d_c * d_c_d_suux;
	d_c2_d_suuy += d_c2_d_c * d_c_d_suuy;
	d_c2_d_suuz += d_c2_d_c * d_c_d_suuz;
	d_c2_d_suvx += d_c2_d_c * d_c_d_suvx;
	d_c2_d_suvy += d_c2_d_c * d_c_d_suvy;
	d_c2_d_suvz += d_c2_d_c * d_c_d_suvz;
	d_c2_d_rux += d_c2_d_b * d_b_d_rux;
	d_c2_d_ruy += d_c2_d_b * d_b_d_ruy;
	d_c2_d_ruz += d_c2_d_b * d_b_d_ruz;
	d_c2_d_rvx += d_c2_d_b * d_b_d_rvx;
	d_c2_d_rvy += d_c2_d_b * d_b_d_rvy;
	d_c2_d_rvz += d_c2_d_b * d_b_d_rvz;
	d_c2_d_suvx += d_c2_d_b * d_b_d_suvx;
	d_c2_d_suvy += d_c2_d_b * d_b_d_suvy;
	d_c2_d_suvz += d_c2_d_b * d_b_d_suvz;
	d_c2_d_svvx += d_c2_d_b * d_b_d_svvx;
	d_c2_d_svvy += d_c2_d_b * d_b_d_svvy;
	d_c2_d_svvz += d_c2_d_b * d_b_d_svvz;
	d_c2_d_rux += d_c2_d_a * d_a_d_rux;
	d_c2_d_ruy += d_c2_d_a * d_a_d_ruy;
	d_c2_d_ruz += d_c2_d_a * d_a_d_ruz;
	d_c2_d_rvx += d_c2_d_a * d_a_d_rvx;
	d_c2_d_rvy += d_c2_d_a * d_a_d_rvy;
	d_c2_d_rvz += d_c2_d_a * d_a_d_rvz;
	d_c2_d_suux += d_c2_d_a * d_a_d_suux;
	d_c2_d_suuy += d_c2_d_a * d_a_d_suuy;
	d_c2_d_suuz += d_c2_d_a * d_a_d_suuz;
	d_c2_d_suvx += d_c2_d_a * d_a_d_suvx;
	d_c2_d_suvy += d_c2_d_a * d_a_d_suvy;
	d_c2_d_suvz += d_c2_d_a * d_a_d_suvz;
	d_e_d_rux += d_e_d_c2 * d_c2_d_rux;
	d_e_d_ruy += d_e_d_c2 * d_c2_d_ruy;
	d_e_d_ruz += d_e_d_c2 * d_c2_d_ruz;
	d_e_d_rvx += d_e_d_c2 * d_c2_d_rvx;
	d_e_d_rvy += d_e_d_c2 * d_c2_d_rvy;
	d_e_d_rvz += d_e_d_c2 * d_c2_d_rvz;
	d_e_d_suux += d_e_d_c2 * d_c2_d_suux;
	d_e_d_suuy += d_e_d_c2 * d_c2_d_suuy;
	d_e_d_suuz += d_e_d_c2 * d_c2_d_suuz;
	d_e_d_suvx += d_e_d_c2 * d_c2_d_suvx;
	d_e_d_suvy += d_e_d_c2 * d_c2_d_suvy;
	d_e_d_suvz += d_e_d_c2 * d_c2_d_suvz;
	d_e_d_svvx += d_e_d_c2 * d_c2_d_svvx;
	d_e_d_svvy += d_e_d_c2 * d_c2_d_svvy;
	d_e_d_svvz += d_e_d_c2 * d_c2_d_svvz;
	d_e_d_rux += d_e_d_c1 * d_c1_d_rux;
	d_e_d_ruy += d_e_d_c1 * d_c1_d_ruy;
	d_e_d_ruz += d_e_d_c1 * d_c1_d_ruz;
	d_e_d_rvx += d_e_d_c1 * d_c1_d_rvx;
	d_e_d_rvy += d_e_d_c1 * d_c1_d_rvy;
	d_e_d_rvz += d_e_d_c1 * d_c1_d_rvz;
	d_e_d_suux += d_e_d_c1 * d_c1_d_suux;
	d_e_d_suuy += d_e_d_c1 * d_c1_d_suuy;
	d_e_d_suuz += d_e_d_c1 * d_c1_d_suuz;
	d_e_d_suvx += d_e_d_c1 * d_c1_d_suvx;
	d_e_d_suvy += d_e_d_c1 * d_c1_d_suvy;
	d_e_d_suvz += d_e_d_c1 * d_c1_d_suvz;
	d_e_d_svvx += d_e_d_c1 * d_c1_d_svvx;
	d_e_d_svvy += d_e_d_c1 * d_c1_d_svvy;
	d_e_d_svvz += d_e_d_c1 * d_c1_d_svvz;
	d_e_d_rux += d_e_d_g * d_g_d_rux;
	d_e_d_ruy += d_e_d_g * d_g_d_ruy;
	d_e_d_ruz += d_e_d_g * d_g_d_ruz;
	d_e_d_rvx += d_e_d_g * d_g_d_rvx;
	d_e_d_rvy += d_e_d_g * d_g_d_rvy;
	d_e_d_rvz += d_e_d_g * d_g_d_rvz;
			

			// dE/du 			

	
			for( int p = 0; p < np; p++ )
			{
				gr[3*cp[p]+0] += d_e_d_rux * ceff_map_du[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_ruy * ceff_map_du[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_ruz * ceff_map_du[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_rvx * ceff_map_dv[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_rvy * ceff_map_dv[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_rvz * ceff_map_dv[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_suux * ceff_map_duu[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_suuy * ceff_map_duu[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_suuz * ceff_map_duu[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_suvx * ceff_map_duv[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_suvy * ceff_map_duv[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_suvz * ceff_map_duv[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_svvx * ceff_map_dvv[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_svvy * ceff_map_dvv[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_svvz * ceff_map_dvv[p] * alpha_z; 
				
				dedr[3*p+0] += d_e_d_rux * ceff_map_du[p];
				dedr[3*p+1] += d_e_d_ruy * ceff_map_du[p];
				dedr[3*p+2] += d_e_d_ruz * ceff_map_du[p];
				                                         
				dedr[3*p+0] += d_e_d_rvx * ceff_map_dv[p];
				dedr[3*p+1] += d_e_d_rvy * ceff_map_dv[p];
				dedr[3*p+2] += d_e_d_rvz * ceff_map_dv[p];
				                                           
				dedr[3*p+0] += d_e_d_suux * ceff_map_duu[p]; 
				dedr[3*p+1] += d_e_d_suuy * ceff_map_duu[p]; 
				dedr[3*p+2] += d_e_d_suuz * ceff_map_duu[p]; 
				                                           
				dedr[3*p+0] += d_e_d_suvx * ceff_map_duv[p]; 
				dedr[3*p+1] += d_e_d_suvy * ceff_map_duv[p]; 
				dedr[3*p+2] += d_e_d_suvz * ceff_map_duv[p]; 
				                                           
				dedr[3*p+0] += d_e_d_svvx * ceff_map_dvv[p]; 
				dedr[3*p+1] += d_e_d_svvy * ceff_map_dvv[p]; 
				dedr[3*p+2] += d_e_d_svvz * ceff_map_dvv[p]; 


			}
			
			for( int p = 0; p < np; p++ )
			{
				gr[3*nv+0] += dedr[3*p+0] * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				gr[3*nv+1] += dedr[3*p+1] * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				gr[3*nv+2] += dedr[3*p+2] * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			}

		}	

	}


	for( int fx = 0; fx < par_info.nf; fx++ )
	{
		int tf = par_info.faces[fx];

		if( tf < nf_faces )
			continue;

		int f = tf - nf_faces;

		int frm = f*nf_irr_pts;

		int t = theIrregularFormulas[frm].tri;
		int i = theIrregularFormulas[frm].vertex;
		int e = theIrregularFormulas[frm].edge;
		int *cp = theIrregularFormulas[frm].cp;
		int np = theIrregularFormulas[frm].ncoor;

		int val = theVertices[i].valence;

		irr_kernel *theKernel = kernels[val]; 
		int ncoords_base = val + 6;
		
		for( int px = 0; px < theTriangles[t].np; px++ )
		{
			int pid = theTriangles[t].plist[px];
			double p_c0 = theTriangles[t].pc0[px];

			double p_area = theTriangles[t].pa[px];
	
			double fu = p_uv[2*pid+0];
			double fv = p_uv[2*pid+1];					
	
			int domain = theKernel->domain(fu,fv);
			if( domain-1 >= theKernel->ndomains )
			{
				printf("LOGICAL ERROR evaluating R,nrm on irregular vertex face.\n");
				exit(1);
			}	
			double *theMap = theKernel->get_map( &fu, &fv );
			double u_u = p_uv[2*pid+0], u_v = 0, v_u = 0, v_v = p_uv[2*pid+1];
			theKernel->get_map_transform( &u_u, &u_v, &v_u, &v_v );
	
			double u = fu;
			double v = fv;
			double w = 1 - u - v;
	
			if( u +v > 1.0 || u < 0 || v < 0 )
			{
				printf("ERROR: outside domain.\n");
				exit(1);
			}
	
			double u2 = u*u;
			double u3 = u*u*u;
			double u4 = u*u*u*u;
			
			double v2 = v*v;
			double v3 = v*v*v;
			double v4 = v*v*v*v;
			
			double w2 = w*w;
			double w3 = w*w*w;
			double w4 = w*w*w*w;
				
			
			// 8 : 0
			// 7 : 1
			// 4 : 2
			// 5 : 3
			// 9 : 4
			// 12 : 5
			// 11 : 6
			// 10 : 7
			// 6 : 8
			// 3 : 9
			// 1 : 10
			// 2 : 11
	
		
	
			double n1 = (1.0/12.0)*(u4+2*u3*v); 
			double n2 = (1.0/12.0)*(u4+2*u3*w); 
			double n3 = (1.0/12.0)*(u4+2*u3*w+6*u3*v+6*u2*v*w+12*u2*v2+6*u*v2*w+6*u*v3+2*v3*w+v4); 
			double n4 = (1.0/12.0)*(6*u4+24*u3*w+24*u2*w2+8*u*w3+w4+24*u3*v+60*u2*v*w+36*u*v*w2+6*v*w3+24*u2*v2+36*u*v2*w+12*v2*w2+8*u*v3+6*v3*w+v4); 
			double n5 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+2*u3*v+6*u2*v*w+6*u*v*w2+2*v*w3); 
			double n6 = (1.0/12.0)*(2*u*v3+v4); 
			double n7 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+8*u3*v+36*u2*v*w+36*u*v*w2+8*v*w3+24*u2*v2+60*u*v2*w+24*v2*w2+24*u*v3+24*v3*w+6*v4); 
			double n8 = (1.0/12.0)*(u4+8*u3*w+24*u2*w2+24*u*w3+6*w4+6*u3*v+36*u2*v*w+60*u*v*w2+24*v*w3+12*u2*v2+36*u*v2*w+24*v2*w2+6*u*v3+8*v3*w+v4); 
			double n9 = (1.0/12.0)*(2*u*w3+w4); 
			double n10 = (1.0/12.0)*(2*v3*w+v4); 
			double n11 = (1.0/12.0)*(2*u*w3+w4+6*u*v*w2+6*v*w3+6*u*v2*w+12*v2*w2+2*u*v3+6*v3*w+v4); 
			double n12 = (1.0/12.0)*(w4+2*v*w3);
	
			double du_1 = Power(u,3)/3. + (Power(u,2)*v)/2.;
			double du_2 = Power(u,2)/2. - Power(u,3)/3. - (Power(u,2)*v)/2.; 	
			double du_3 =Power(u,2)/2. - Power(u,3)/3. + u*v - (Power(u,2)*v)/2. + Power(v,2)/2. - Power(v,3)/6.;
			double du_4 = 0.3333333333333333 + u - Power(u,2) - Power(u,3)/3. + v/2. - u*v - (Power(u,2)*v)/2. - Power(v,2) + Power(v,3)/3.;	
			double du_5 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. - v/2. + Power(u,2)*v + Power(v,2)/2. - Power(v,3)/6.;
			double du_6 = Power(v,3)/6.;	
			double du_7 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. + v/2. - 2*u*v + Power(u,2)*v - Power(v,2)/2. - Power(v,3)/6.;	
			double du_8 = -2*u + 2*Power(u,2) - Power(u,3)/3. - v + 2*u*v - (Power(u,2)*v)/2. + Power(v,2) - Power(v,3)/6.;	
			double du_9 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. + v/2. - (Power(u,2)*v)/2. - Power(v,2)/2. + Power(v,3)/6.;	
			double du_10 = -Power(v,3)/6.;
			double du_11 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. - v/2. + u*v - (Power(u,2)*v)/2. + Power(v,3)/3.;	
			double du_12 = 	-0.3333333333333333 + u - Power(u,2) + Power(u,3)/3. + v/2. - u*v + (Power(u,2)*v)/2. - Power(v,3)/6.;
	
			double dv_1 = Power(u,3)/6.;
			double dv_2 = -Power(u,3)/6.;
			double dv_3 = Power(u,2)/2. - Power(u,3)/6. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_4 = 0.16666666666666666 + u/2. - Power(u,2)/2. - Power(u,3)/6. - 2*u*v - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
			double dv_5 = -0.16666666666666666 - u/2. + Power(u,3)/3. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_6 = (u*Power(v,2))/2. + Power(v,3)/3.;
			double dv_7 = 0.3333333333333333 + u/2. - Power(u,2) + Power(u,3)/3. + v - u*v - Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_8 = -u + Power(u,2) - Power(u,3)/6. - 2*v + 2*u*v + 2*Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_9 = -0.3333333333333333 + u/2. - Power(u,3)/6. + v - u*v - Power(v,2) + (u*Power(v,2))/2. + Power(v,3)/3.;
			double dv_10 = Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
			double dv_11 = 0.16666666666666666 - u/2. + Power(u,2)/2. - Power(u,3)/6. - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
			double dv_12 = -0.16666666666666666 + u/2. - Power(u,2)/2. + Power(u,3)/6. + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		
			double d_uu_1 = Power(u,2) + u*v;
			double d_uu_2 = u - Power(u,2) - u*v;
			double d_uu_3 = u - Power(u,2) + v - u*v;
			double d_uu_4 = 1 - 2*u - Power(u,2) - v - u*v;
			double d_uu_5 = -2*u + 2*Power(u,2) + 2*u*v;
			double d_uu_6 = 0;
			double d_uu_7 = -2*u + 2*Power(u,2) - 2*v + 2*u*v;
			double d_uu_8 = -2 + 4*u - Power(u,2) + 2*v - u*v;
			double d_uu_9 = u - Power(u,2) - u*v;
			double d_uu_10 = 0;
			double d_uu_11 = u - Power(u,2) + v - u*v;
			double d_uu_12 = 1 - 2*u + Power(u,2) - v + u*v;
			
			double d_uv_1 = Power(u,2)/2.;
			double d_uv_2 = -Power(u,2)/2.;
			double d_uv_3 = u - Power(u,2)/2. + v - Power(v,2)/2.;
			double d_uv_4 = 0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
			double d_uv_5 = -0.5 + Power(u,2) + v - Power(v,2)/2.;
			double d_uv_6 = Power(v,2)/2.;
			double d_uv_7 = 0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
			double d_uv_8 = -1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
			double d_uv_9 = 0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
			double d_uv_10 = -Power(v,2)/2.;
			double d_uv_11 = -0.5 + u - Power(u,2)/2. + Power(v,2);
			double d_uv_12 = 0.5 - u + Power(u,2)/2. - Power(v,2)/2.;
			
			double d_vv_1 = 0;
			double d_vv_2 = 0;
			double d_vv_3 = u + v - u*v - Power(v,2);
			double d_vv_4 = -2*u - 2*v + 2*u*v + 2*Power(v,2);
			double d_vv_5 = u + v - u*v - Power(v,2);
			double d_vv_6 = u*v + Power(v,2);
			double d_vv_7 = 1 - u - 2*v - u*v - Power(v,2);
			double d_vv_8 = -2 + 2*u + 4*v - u*v - Power(v,2);
			double d_vv_9 = 1 - u - 2*v + u*v + Power(v,2);
			double d_vv_10 = v - u*v - Power(v,2);
			double d_vv_11 = -2*v + 2*u*v + 2*Power(v,2);
			double d_vv_12 = v - u*v - Power(v,2);

			double d_uuu_1 = 2*u + v;
			double d_uuu_2 = 1-2*u-v; //u - Power(u,2) - u*v;
			double d_uuu_3 = 1-2*u-v;//u - Power(u,2) + v - u*v;
			double d_uuu_4 = -2-2*u-v;//1 - 2*u - Power(u,2) - v - u*v;
			double d_uuu_5 = -2+4*u+2*v;//-2*u + 2*Power(u,2) + 2*u*v;
			double d_uuu_6 = 0;
			double d_uuu_7 = -2+4*u+2*v;//-2*u + 2*Power(u,2) - 2*v + 2*u*v;
			double d_uuu_8 = 4-2*u-v;//-2 + 4*u - Power(u,2) + 2*v - u*v;
			double d_uuu_9 = 1-2*u-v;//u - Power(u,2) - u*v;
			double d_uuu_10 = 0;
			double d_uuu_11 = 1-2*u-v;//u - Power(u,2) + v - u*v;
			double d_uuu_12 = -2+2*u+v;//1 - 2*u + Power(u,2) - v + u*v;
	
			double d_uuv_1 = u;//Power(u,2)/2.;
			double d_uuv_2 = -u;//-Power(u,2)/2.;
			double d_uuv_3 = 1-u;//u - Power(u,2)/2. + v - Power(v,2)/2.;
			double d_uuv_4 = -1-u;//0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
			double d_uuv_5 = 2*u;//-0.5 + Power(u,2) + v - Power(v,2)/2.;
			double d_uuv_6 = 0;//Power(v,2)/2.;
			double d_uuv_7 = -2+2*u;//0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
			double d_uuv_8 = 2-u;//-1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
			double d_uuv_9 = -u;//0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
			double d_uuv_10 = 0;//-Power(v,2)/2.;
			double d_uuv_11 = 1-u;//-0.5 + u - Power(u,2)/2. + Power(v,2);
			double d_uuv_12 = -1+u;//0.5 - u + Power(u,2)/2. - Power(v,2)/2.;

			double d_uvv_1 = 0;
			double d_uvv_2 = 0;
			double d_uvv_3 = 1-v;//u + v - u*v - Power(v,2);
			double d_uvv_4 = -2+2*v;//-2*u - 2*v + 2*u*v + 2*Power(v,2);
			double d_uvv_5 = 1-v;//u + v - u*v - Power(v,2);
			double d_uvv_6 = v;//u*v + Power(v,2);
			double d_uvv_7 = -1-v;//1 - u - 2*v - u*v - Power(v,2);
			double d_uvv_8 = 2-v;//-2 + 2*u + 4*v - u*v - Power(v,2);
			double d_uvv_9 = -1+v;//1 - u - 2*v + u*v + Power(v,2);
			double d_uvv_10 = -v;//v - u*v - Power(v,2);
			double d_uvv_11 = 2*v;//-2*v + 2*u*v + 2*Power(v,2);
			double d_uvv_12 = -v;//v - u*v - Power(v,2);
			
			double d_vvv_1 = 0;
			double d_vvv_2 = 0;
			double d_vvv_3 = 1-u-2*v;//u + v - u*v - Power(v,2);
			double d_vvv_4 = -2+2*u+4*v;//-2*u - 2*v + 2*u*v + 2*Power(v,2);
			double d_vvv_5 = 1-u-2*v;//u + v - u*v - Power(v,2);
			double d_vvv_6 = u+2*v;//u*v + Power(v,2);
			double d_vvv_7 = -2-u-2*v;//1 - u - 2*v - u*v - Power(v,2);
			double d_vvv_8 = 4-u-2*v;//-2 + 2*u + 4*v - u*v - Power(v,2);
			double d_vvv_9 = -2+u+2*v;//1 - u - 2*v + u*v + Power(v,2);
			double d_vvv_10 = 1-u-2*v;//v - u*v - Power(v,2);
			double d_vvv_11 = -2+2*u+4*v;//-2*v + 2*u*v + 2*Power(v,2);
			double d_vvv_12 = 1-u-2*v;//v - u*v - Power(v,2);
			
			double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
			double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
			double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
			
			double ceff_map_duu[12] = { d_uu_8, d_uu_7, d_uu_4, d_uu_5, d_uu_9, d_uu_12, d_uu_11, d_uu_10, d_uu_6, d_uu_3, d_uu_1, d_uu_2 };
			double ceff_map_duv[12] = { d_uv_8, d_uv_7, d_uv_4, d_uv_5, d_uv_9, d_uv_12, d_uv_11, d_uv_10, d_uv_6, d_uv_3, d_uv_1, d_uv_2 };
			double ceff_map_dvv[12] = { d_vv_8, d_vv_7, d_vv_4, d_vv_5, d_vv_9, d_vv_12, d_vv_11, d_vv_10, d_vv_6, d_vv_3, d_vv_1, d_vv_2 };
			
			double ceff_map_duuu[12] = { d_uuu_8, d_uuu_7, d_uuu_4, d_uuu_5, d_uuu_9, d_uuu_12, d_uuu_11, d_uuu_10, d_uuu_6, d_uuu_3, d_uuu_1, d_uuu_2 };
			double ceff_map_duuv[12] = { d_uuv_8, d_uuv_7, d_uuv_4, d_uuv_5, d_uuv_9, d_uuv_12, d_uuv_11, d_uuv_10, d_uuv_6, d_uuv_3, d_uuv_1, d_uuv_2 };
			double ceff_map_duvv[12] = { d_uvv_8, d_uvv_7, d_uvv_4, d_uvv_5, d_uvv_9, d_uvv_12, d_uvv_11, d_uvv_10, d_uvv_6, d_uvv_3, d_uvv_1, d_uvv_2 };
			double ceff_map_dvvv[12] = { d_vvv_8, d_vvv_7, d_vvv_4, d_vvv_5, d_vvv_9, d_vvv_12, d_vvv_11, d_vvv_10, d_vvv_6, d_vvv_3, d_vvv_1, d_vvv_2 };
	
			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
			double nrm[3]={0,0,0}; 
			double tSuuu[3] = {0,0,0};
			double tSuuv[3] = {0,0,0};
			double tSuvv[3] = {0,0,0};
			double tSvvv[3] = {0,0,0};
	
			int *cset = theVertices[i].irr_coord_set + e * ncoords_base;
	
			for( int x = 0; x < ncoords_base; x++ )
			{
				for( int y = 0; y < 12; y++ )
				{
					R[0] += (r[3*cset[x]+0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map[y] * alpha_x;
					R[1] += (r[3*cset[x]+1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map[y] * alpha_y;
					R[2] += (r[3*cset[x]+2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map[y] * alpha_z;
					
					Ru[0] += (r[3*cset[x]+0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_du[y] * alpha_x * u_u;
					Ru[1] += (r[3*cset[x]+1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_du[y] * alpha_y * u_u;
					Ru[2] += (r[3*cset[x]+2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_du[y] * alpha_z * u_u;
					
					Rv[0] += (r[3*cset[x]+0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_dv[y] * alpha_x * v_v;
					Rv[1] += (r[3*cset[x]+1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_dv[y] * alpha_y * v_v;
					Rv[2] += (r[3*cset[x]+2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_dv[y] * alpha_z * v_v;
					
					tSuu[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_duu[y] * alpha_x * u_u * u_u;
					tSuu[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_duu[y] * alpha_y * u_u * u_u;
					tSuu[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_duu[y] * alpha_z * u_u * u_u;
					
					tSuv[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_duv[y] * alpha_x * u_u * v_v;
					tSuv[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_duv[y] * alpha_y * u_u * v_v;
					tSuv[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_duv[y] * alpha_z * u_u * v_v;
					
					tSvv[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_dvv[y] * alpha_x * v_v * v_v;
					tSvv[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_dvv[y] * alpha_y * v_v * v_v;
					tSvv[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_dvv[y] * alpha_z * v_v * v_v;
			
					tSuuu[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_duuu[y] * alpha_x * u_u * u_u * u_u;
					tSuuu[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_duuu[y] * alpha_y * u_u * u_u * u_u;
					tSuuu[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_duuu[y] * alpha_z * u_u * u_u * u_u;
					
					tSuuv[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_duuv[y] * alpha_x * u_u * u_u * v_v;
					tSuuv[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_duuv[y] * alpha_y * u_u * u_u * v_v;
					tSuuv[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_duuv[y] * alpha_z * u_u * u_u * v_v;
					
					tSuvv[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_duvv[y] * alpha_x * u_u * v_v * v_v;
					tSuvv[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_duvv[y] * alpha_y * u_u * v_v * v_v;
					tSuvv[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_duvv[y] * alpha_z * u_u * v_v * v_v;
					
					tSvvv[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_dvvv[y] * alpha_x * v_v * v_v * v_v;
					tSvvv[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_dvvv[y] * alpha_y * v_v * v_v * v_v;
					tSvvv[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_dvvv[y] * alpha_z * v_v * v_v * v_v;
	
				}
			}
			
			cross( Ru, Rv, nrm );
			normalize(nrm);
	
	
			double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
			double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
			double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];
	
			double g = sqrt(RuRu*RvRv-RuRv*RuRv);
	
			double e1,e2;
	      
			double nsuu = tSuu[0] * nrm[0] + tSuu[1] * nrm[1] + tSuu[2] * nrm[2];
			double nsuv = tSuv[0] * nrm[0] + tSuv[1] * nrm[1] + tSuv[2] * nrm[2];
			double nsvv = tSvv[0] * nrm[0] + tSvv[1] * nrm[1] + tSvv[2] * nrm[2];
	
			double Stot = (nsuu * RvRv + nsvv * RuRu -2*nsuv*RuRv)/(g*g);
	
			double Sop[4] = { 1.0/(g*g) * (nsuu * RvRv  - nsuv * RuRv), (1.0/(g*g)) * ( nsuv*RvRv-nsvv*RuRv),
					  1.0/(g*g) * (nsuv * RuRu  - nsuu * RuRv), (1.0/(g*g)) * ( nsvv*RuRu-nsuv*RuRv) };
	
			double a = Sop[0];		
			double b = Sop[1];
			double c = Sop[2];
			double d = Sop[3];
	
			double c0 = theFormulas[frm].c0;
			e1 = -0.5*(a+d-sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
			e2 = -0.5*(a+d+sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
		

double d_nrmz_d_rux=0,d_nrmz_d_ruy=0,d_nrmz_d_ruz=0,d_nrmz_d_rvx=0,d_nrmz_d_rvy=0,d_nrmz_d_rvz=0,d_nrmy_d_rux=0,d_nrmy_d_ruy=0,d_nrmy_d_ruz=0,d_nrmy_d_rvx=0,d_nrmy_d_rvy=0,d_nrmy_d_rvz=0,d_nrmx_d_rux=0,d_nrmx_d_ruy=0,d_nrmx_d_ruz=0,d_nrmx_d_rvx=0,d_nrmx_d_rvy=0,d_nrmx_d_rvz=0,d_RuRu_d_rux=0,d_RuRu_d_ruy=0,d_RuRu_d_ruz=0,d_nsvv_d_rux=0,d_nsvv_d_ruy=0,d_nsvv_d_ruz=0,d_nsvv_d_rvx=0,d_nsvv_d_rvy=0,d_nsvv_d_rvz=0,d_nsvv_d_svvx=0,d_nsvv_d_svvy=0,d_nsvv_d_svvz=0,d_RuRv_d_rux=0,d_RuRv_d_ruy=0,d_RuRv_d_ruz=0,d_RuRv_d_rvx=0,d_RuRv_d_rvy=0,d_RuRv_d_rvz=0,d_nsuv_d_rux=0,d_nsuv_d_ruy=0,d_nsuv_d_ruz=0,d_nsuv_d_rvx=0,d_nsuv_d_rvy=0,d_nsuv_d_rvz=0,d_nsuv_d_suvx=0,d_nsuv_d_suvy=0,d_nsuv_d_suvz=0,d_RvRv_d_rvx=0,d_RvRv_d_rvy=0,d_RvRv_d_rvz=0,d_nsuu_d_rux=0,d_nsuu_d_ruy=0,d_nsuu_d_ruz=0,d_nsuu_d_rvx=0,d_nsuu_d_rvy=0,d_nsuu_d_rvz=0,d_nsuu_d_suux=0,d_nsuu_d_suuy=0,d_nsuu_d_suuz=0,d_g_d_rux=0,d_g_d_ruy=0,d_g_d_ruz=0,d_g_d_rvx=0,d_g_d_rvy=0,d_g_d_rvz=0,d_e_d_rux=0,d_e_d_ruy=0,d_e_d_ruz=0,d_e_d_rvx=0,d_e_d_rvy=0,d_e_d_rvz=0,d_e_d_suux=0,d_e_d_suuy=0,d_e_d_suuz=0,d_e_d_suvx=0,d_e_d_suvy=0,d_e_d_suvz=0,d_e_d_svvx=0,d_e_d_svvy=0,d_e_d_svvz=0,d_c2_d_rux=0,d_c2_d_ruy=0,d_c2_d_ruz=0,d_c2_d_rvx=0,d_c2_d_rvy=0,d_c2_d_rvz=0,d_c2_d_suux=0,d_c2_d_suuy=0,d_c2_d_suuz=0,d_c2_d_suvx=0,d_c2_d_suvy=0,d_c2_d_suvz=0,d_c2_d_svvx=0,d_c2_d_svvy=0,d_c2_d_svvz=0,d_d_d_rux=0,d_d_d_ruy=0,d_d_d_ruz=0,d_d_d_rvx=0,d_d_d_rvy=0,d_d_d_rvz=0,d_d_d_suvx=0,d_d_d_suvy=0,d_d_d_suvz=0,d_d_d_svvx=0,d_d_d_svvy=0,d_d_d_svvz=0,d_c_d_rux=0,d_c_d_ruy=0,d_c_d_ruz=0,d_c_d_rvx=0,d_c_d_rvy=0,d_c_d_rvz=0,d_c_d_suux=0,d_c_d_suuy=0,d_c_d_suuz=0,d_c_d_suvx=0,d_c_d_suvy=0,d_c_d_suvz=0,d_b_d_rux=0,d_b_d_ruy=0,d_b_d_ruz=0,d_b_d_rvx=0,d_b_d_rvy=0,d_b_d_rvz=0,d_b_d_suvx=0,d_b_d_suvy=0,d_b_d_suvz=0,d_b_d_svvx=0,d_b_d_svvy=0,d_b_d_svvz=0,d_a_d_rux=0,d_a_d_ruy=0,d_a_d_ruz=0,d_a_d_rvx=0,d_a_d_rvy=0,d_a_d_rvz=0,d_a_d_suux=0,d_a_d_suuy=0,d_a_d_suuz=0,d_a_d_suvx=0,d_a_d_suvy=0,d_a_d_suvz=0,d_c1_d_rux=0,d_c1_d_ruy=0,d_c1_d_ruz=0,d_c1_d_rvx=0,d_c1_d_rvy=0,d_c1_d_rvz=0,d_c1_d_suux=0,d_c1_d_suuy=0,d_c1_d_suuz=0,d_c1_d_suvx=0,d_c1_d_suvy=0,d_c1_d_suvz=0,d_c1_d_svvx=0,d_c1_d_svvy=0,d_c1_d_svvz=0,d_g_d_RvRv=0,d_g_d_RuRv=0,d_g_d_RuRu=0,d_nsvv_d_nrmz=0,d_nsvv_d_nrmy=0,d_nsvv_d_nrmx=0,d_nsuv_d_nrmz=0,d_nsuv_d_nrmy=0,d_nsuv_d_nrmx=0,d_nsuu_d_nrmz=0,d_nsuu_d_nrmy=0,d_nsuu_d_nrmx=0,d_d_d_RuRv=0,d_d_d_nsuv=0,d_d_d_RuRu=0,d_d_d_nsvv=0,d_d_d_g=0,d_c_d_RuRv=0,d_c_d_nsuv=0,d_c_d_RuRu=0,d_c_d_nsuu=0,d_c_d_g=0,d_b_d_RuRv=0,d_b_d_nsuv=0,d_b_d_RvRv=0,d_b_d_nsvv=0,d_b_d_g=0,d_a_d_RuRv=0,d_a_d_nsuv=0,d_a_d_RvRv=0,d_a_d_nsuu=0,d_a_d_g=0,d_e_d_c2=0,d_e_d_c1=0,d_e_d_g=0,d_c2_d_d=0,d_c2_d_c=0,d_c2_d_b=0,d_c2_d_a=0,d_c1_d_d=0,d_c1_d_c=0,d_c1_d_b=0,d_c1_d_a=0,junk;
			// the basic variables.





			d_g_d_RuRu = RvRv/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
			d_g_d_RuRv = -(RuRv/Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
			d_g_d_RvRv = RuRu/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  

			d_RuRu_d_rux = 2*Ru[0];
			d_RuRu_d_ruy = 2*Ru[1];
			d_RuRu_d_ruz = 2*Ru[2];
			
			d_RuRv_d_rux = Rv[0];
			d_RuRv_d_ruy = Rv[1];
			d_RuRv_d_ruz = Rv[2];
			
			d_RuRv_d_rvx = Ru[0];
			d_RuRv_d_rvy = Ru[1];
			d_RuRv_d_rvz = Ru[2];
			
			d_RvRv_d_rvx = 2*Rv[0];
			d_RvRv_d_rvy = 2*Rv[1];
			d_RvRv_d_rvz = 2*Rv[2];

			// intermediates.
			d_a_d_g = (-2*(-(nsuv*RuRv) + nsuu*RvRv))/Power(g,3);
			d_a_d_nsuu = RvRv/Power(g,2);
			d_a_d_RvRv = nsuu/Power(g,2);
			d_a_d_nsuv = -(RuRv/Power(g,2));
			d_a_d_RuRv = -(nsuv/Power(g,2));

			d_b_d_g = (-2*(-(nsvv*RuRv) + nsuv*RvRv))/Power(g,3);
			d_b_d_nsvv = -(RuRv/Power(g,2));
			d_b_d_RvRv = nsuv/Power(g,2);
			d_b_d_nsuv = RvRv/Power(g,2);
			d_b_d_RuRv = -(nsvv/Power(g,2));

			d_c_d_g = (-2*(nsuv*RuRu - nsuu*RuRv))/Power(g,3);
			d_c_d_nsuu = -(RuRv/Power(g,2));
			d_c_d_RuRu = nsuv/Power(g,2);
			d_c_d_nsuv = RuRu/Power(g,2);
			d_c_d_RuRv = -(nsuu/Power(g,2));

			d_d_d_g = (-2*(nsvv*RuRu - nsuv*RuRv))/Power(g,3);
			d_d_d_nsvv = RuRu/Power(g,2);
			d_d_d_RuRu = nsvv/Power(g,2);
			d_d_d_nsuv = -(RuRv/Power(g,2));
			d_d_d_RuRv = -(nsuv/Power(g,2));

			d_c1_d_a =-0.5*(1 - (2*a - 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
			d_c1_d_b =-(-1.*c)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2)) ;
			d_c1_d_c =-(-1.*b)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2)) ;
			d_c1_d_d =-0.5*(1 - (-2*a + 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
		
			d_c2_d_a = -0.5*(1 + (2*a - 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
			d_c2_d_b = -(1.*c)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2));
			d_c2_d_c = -(1.*b)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2));
			d_c2_d_d = -0.5*(1 + (-2*a + 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));


			// this is the term we need to gradient:
			// en += 0.5*kc * p_area * ( c1 + c2 - p_c0 ) * ( c1 + c2 - p_c0 );
			d_e_d_c1 = kc * p_area * ( e1 + e2 - p_c0);
			d_e_d_c2 = kc * p_area * ( e1 + e2 - p_c0);
			d_e_d_g = 0;

			// that's it.

			double fac = Power(Power(Ru[1]*Rv[0] - Ru[0]*Rv[1],2) + Power(Ru[2]*Rv[0] - Ru[0]*Rv[2],2) + Power(Ru[2]*Rv[1] - Ru[1]*Rv[2],2),1.5);

			d_nrmx_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac; 
			d_nrmx_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2)))/fac;
			d_nrmx_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac;
			
			d_nrmy_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2))))/fac;
			d_nrmy_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac; 
			d_nrmy_d_ruz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2])) + Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac;
	
			d_nrmz_d_rux = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2])/fac;
			d_nrmz_d_ruy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;
			d_nrmz_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2))) + (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;

			d_nrmx_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac; 
			d_nrmx_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2])))/fac;
			d_nrmx_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac;
			
			d_nrmy_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
			d_nrmy_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2])))/fac;
			d_nrmy_d_rvz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[0]*Ru[1]*Rv[0] - Power(Ru[0],2)*Rv[1] + Ru[2]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2])))/fac;
	
			d_nrmz_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;
			d_nrmz_d_rvy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[0]*Ru[2]*Rv[0] - Power(Ru[0],2)*Rv[2] + Ru[1]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
			d_nrmz_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;

			double d_nrmx_d_u = d_nrmx_d_rux * tSuu[0] + d_nrmx_d_ruy * tSuu[1] + d_nrmx_d_ruz * tSuu[2] +
					    d_nrmx_d_rvx * tSuv[0] + d_nrmx_d_rvy * tSuv[1] + d_nrmx_d_rvz * tSuv[2];
			double d_nrmx_d_v = d_nrmx_d_rux * tSuv[0] + d_nrmx_d_ruy * tSuv[1] + d_nrmx_d_ruz * tSuv[2] +
					    d_nrmx_d_rvx * tSvv[0] + d_nrmx_d_rvy * tSvv[1] + d_nrmx_d_rvz * tSvv[2];
			double d_nrmy_d_u = d_nrmy_d_rux * tSuu[0] + d_nrmy_d_ruy * tSuu[1] + d_nrmy_d_ruz * tSuu[2] +
					    d_nrmy_d_rvx * tSuv[0] + d_nrmy_d_rvy * tSuv[1] + d_nrmy_d_rvz * tSuv[2];
			double d_nrmy_d_v = d_nrmy_d_rux * tSuv[0] + d_nrmy_d_ruy * tSuv[1] + d_nrmy_d_ruz * tSuv[2] +
					    d_nrmy_d_rvx * tSvv[0] + d_nrmy_d_rvy * tSvv[1] + d_nrmy_d_rvz * tSvv[2];
			double d_nrmz_d_u = d_nrmz_d_rux * tSuu[0] + d_nrmz_d_ruy * tSuu[1] + d_nrmz_d_ruz * tSuu[2] +
					    d_nrmz_d_rvx * tSuv[0] + d_nrmz_d_rvy * tSuv[1] + d_nrmz_d_rvz * tSuv[2];
			double d_nrmz_d_v = d_nrmz_d_rux * tSuv[0] + d_nrmz_d_ruy * tSuv[1] + d_nrmz_d_ruz * tSuv[2] +
					    d_nrmz_d_rvx * tSvv[0] + d_nrmz_d_rvy * tSvv[1] + d_nrmz_d_rvz * tSvv[2];

			double d_nsuu_d_u = tSuuu[0] * nrm[0] + tSuuu[1] * nrm[1] + tSuuu[2] * nrm[2];
			       d_nsuu_d_u += tSuu[0] * d_nrmx_d_u + tSuu[1] * d_nrmy_d_u + tSuu[2] * d_nrmz_d_u;
			double d_nsuu_d_v = tSuuv[0] * nrm[0] + tSuuv[1] * nrm[1] + tSuuv[2] * nrm[2];
			       d_nsuu_d_v += tSuu[0] * d_nrmx_d_v + tSuu[1] * d_nrmy_d_v + tSuu[2] * d_nrmz_d_v;
			
			double d_nsuv_d_u = tSuuv[0] * nrm[0] + tSuuv[1] * nrm[1] + tSuuv[2] * nrm[2];
			       d_nsuv_d_u += tSuv[0] * d_nrmx_d_u + tSuv[1] * d_nrmy_d_u + tSuv[2] * d_nrmz_d_u;
			double d_nsuv_d_v = tSuvv[0] * nrm[0] + tSuvv[1] * nrm[1] + tSuvv[2] * nrm[2];
			       d_nsuv_d_v += tSuv[0] * d_nrmx_d_v + tSuv[1] * d_nrmy_d_v + tSuv[2] * d_nrmz_d_v;

			double d_nsvv_d_u = tSuvv[0] * nrm[0] + tSuvv[1] * nrm[1] + tSuvv[2] * nrm[2];
			       d_nsvv_d_u += tSvv[0] * d_nrmx_d_u + tSvv[1] * d_nrmy_d_u + tSvv[2] * d_nrmz_d_u;
			double d_nsvv_d_v = tSvvv[0] * nrm[0] + tSvvv[1] * nrm[1] + tSvvv[2] * nrm[2];
			       d_nsvv_d_v += tSvv[0] * d_nrmx_d_v + tSvv[1] * d_nrmy_d_v + tSvv[2] * d_nrmz_d_v;

			double d_RuRu_d_u = 2*Ru[0] * tSuu[0] + 2*Ru[1] * tSuu[1] + 2 * Ru[2] * tSuu[2];
			double d_RuRu_d_v = 2*Ru[0] * tSuv[0] + 2*Ru[1] * tSuv[1] + 2 * Ru[2] * tSuv[2];

			double d_RuRv_d_u =  Ru[0] * tSuv[0] + Ru[1] * tSuv[1] + Ru[2] * tSuv[2];
			       d_RuRv_d_u += Rv[0] * tSuu[0] + Rv[1] * tSuu[1] + Rv[2] * tSuu[2];
			
			double d_RuRv_d_v =  Ru[0] * tSvv[0] + Ru[1] * tSvv[1] + Ru[2] * tSvv[2];
			       d_RuRv_d_v += Rv[0] * tSuv[0] + Rv[1] * tSuv[1] + Rv[2] * tSuv[2];

			double d_RvRv_d_u = 2*Rv[0] * tSuv[0] + 2*Rv[1] * tSuv[1] + 2 * Rv[2] * tSuv[2];
			double d_RvRv_d_v = 2*Rv[0] * tSvv[0] + 2*Rv[1] * tSvv[1] + 2 * Rv[2] * tSvv[2];

			double d_g_d_u = 0;
			double d_g_d_v = 0;

			d_g_d_u += d_g_d_RuRu * 2* Ru[0] * tSuu[0] + d_g_d_RuRu * 2* Ru[1] * tSuu[1] + d_g_d_RuRu * 2* Ru[2] * tSuu[2];
			d_g_d_v += d_g_d_RuRu * 2* Ru[0] * tSuv[0] + d_g_d_RuRu * 2* Ru[1] * tSuv[1] + d_g_d_RuRu * 2* Ru[2] * tSuv[2];

			// d_ru dru_du
			d_g_d_u += d_g_d_RuRv * Rv[0] * tSuu[0] + d_g_d_RuRv * Rv[1] * tSuu[1] + d_g_d_RuRv * Rv[2] * tSuu[2];
			// d_rv_drv_du
			d_g_d_u += d_g_d_RuRv * Ru[0] * tSuv[0] + d_g_d_RuRv * Ru[1] * tSuv[1] + d_g_d_RuRv * Ru[2] * tSuv[2];	
			// d_ru dru_dv
			d_g_d_v += d_g_d_RuRv * Rv[0] * tSuv[0] + d_g_d_RuRv * Rv[1] * tSuv[1] + d_g_d_RuRv * Rv[2] * tSuv[2];
			// d_rv drv_dv
			d_g_d_v += d_g_d_RuRv * Ru[0] * tSvv[0] + d_g_d_RuRv * Ru[1] * tSvv[1] + d_g_d_RuRv * Ru[2] * tSvv[2];
			
			d_g_d_u += d_g_d_RvRv * 2* Rv[0] * tSuv[0] + d_g_d_RvRv * 2* Rv[1] * tSuv[1] + d_g_d_RvRv * 2* Rv[2] * tSuv[2];
			d_g_d_v += d_g_d_RvRv * 2* Rv[0] * tSvv[0] + d_g_d_RvRv * 2* Rv[1] * tSvv[1] + d_g_d_RvRv * 2* Rv[2] * tSvv[2];

//			d_a_d_g = (-2*(-(nsuv*RuRv) + nsuu*RvRv))/Power(g,3);
//			d_a_d_nsuu = RvRv/Power(g,2);
//			d_a_d_RvRv = nsuu/Power(g,2);
//			d_a_d_nsuv = -(RuRv/Power(g,2));
//			d_a_d_RuRv = -(nsuv/Power(g,2));
			
			double d_a_d_u = d_a_d_nsuu * d_nsuu_d_u + d_a_d_nsuv * d_nsuv_d_u + d_a_d_RuRv * d_RuRv_d_u + d_a_d_RvRv * d_RvRv_d_u + d_a_d_g * d_g_d_u;
			double d_a_d_v = d_a_d_nsuu * d_nsuu_d_v + d_a_d_nsuv * d_nsuv_d_v + d_a_d_RuRv * d_RuRv_d_v + d_a_d_RvRv * d_RvRv_d_v + d_a_d_g * d_g_d_v;
			double d_b_d_u = d_b_d_nsvv * d_nsvv_d_u + d_b_d_nsuv * d_nsuv_d_u + d_b_d_RuRv * d_RuRv_d_u + d_b_d_RvRv * d_RvRv_d_u + d_b_d_g * d_g_d_u;
			double d_b_d_v = d_b_d_nsvv * d_nsvv_d_v + d_b_d_nsuv * d_nsuv_d_v + d_b_d_RuRv * d_RuRv_d_v + d_b_d_RvRv * d_RvRv_d_v + d_b_d_g * d_g_d_v;
			double d_c_d_u = d_c_d_nsuu * d_nsuu_d_u + d_c_d_nsuv * d_nsuv_d_u + d_c_d_RuRv * d_RuRv_d_u + d_c_d_RuRu * d_RuRu_d_u + d_c_d_g * d_g_d_u;
			double d_c_d_v = d_c_d_nsuu * d_nsuu_d_v + d_c_d_nsuv * d_nsuv_d_v + d_c_d_RuRv * d_RuRv_d_v + d_c_d_RuRu * d_RuRu_d_v + d_c_d_g * d_g_d_v;
			double d_d_d_u = d_d_d_nsvv * d_nsvv_d_u + d_d_d_nsuv * d_nsuv_d_u + d_d_d_RuRv * d_RuRv_d_u + d_d_d_RuRu * d_RuRu_d_u + d_d_d_g * d_g_d_u;
			double d_d_d_v = d_d_d_nsvv * d_nsvv_d_v + d_d_d_nsuv * d_nsuv_d_v + d_d_d_RuRv * d_RuRv_d_v + d_d_d_RuRu * d_RuRu_d_v + d_d_d_g * d_g_d_v;

			double d_c1_d_u = d_c1_d_a * d_a_d_u + d_c1_d_b * d_b_d_u + d_c1_d_c * d_c_d_u + d_c1_d_d * d_d_d_u;
			double d_c2_d_u = d_c2_d_a * d_a_d_u + d_c2_d_b * d_b_d_u + d_c2_d_c * d_c_d_u + d_c2_d_d * d_d_d_u;
	
			double d_c1_d_v = d_c1_d_a * d_a_d_v + d_c1_d_b * d_b_d_v + d_c1_d_c * d_c_d_v + d_c1_d_d * d_d_d_v;
			double d_c2_d_v = d_c2_d_a * d_a_d_v + d_c2_d_b * d_b_d_v + d_c2_d_c * d_c_d_v + d_c2_d_d * d_d_d_v;

			double d_e_d_u = d_e_d_c1 * d_c1_d_u + d_e_d_c2 * d_c2_d_u;
			double d_e_d_v = d_e_d_c1 * d_c1_d_v + d_e_d_c2 * d_c2_d_v;

			pg[2*pid+0] += d_e_d_u;
			pg[2*pid+1] += d_e_d_v;
	
//			printf("u %le v %le a: %.14le da_du: %le nrmx: %.14le dnrmx_du %.14le\n", u, v, b, d_b_d_u, nrm[0], d_nrmy_d_u );

			double e = 0.5* kc * p_area * ( e1 + e2 - p_c0 ) * ( e1 + e2 - p_c0 );
//			printf("e %le at uv: %le %le  pg: %le %le\n", e, p_uv[2*pid+0], p_uv[2*pid+1], d_e_d_u, d_e_d_v );
			
			d_nsuu_d_nrmx = tSuu[0];
			d_nsuu_d_nrmy = tSuu[1];
			d_nsuu_d_nrmz = tSuu[2];
			
			d_nsuv_d_nrmx = tSuv[0];
			d_nsuv_d_nrmy = tSuv[1];
			d_nsuv_d_nrmz = tSuv[2];

			d_nsvv_d_nrmx = tSvv[0];
			d_nsvv_d_nrmy = tSvv[1];
			d_nsvv_d_nrmz = tSvv[2];

			d_nsuu_d_suux = nrm[0];
			d_nsuu_d_suuy = nrm[1];
			d_nsuu_d_suuz = nrm[2];
			
			d_nsuv_d_suvx = nrm[0];
			d_nsuv_d_suvy = nrm[1];
			d_nsuv_d_suvz = nrm[2];

			d_nsvv_d_svvx = nrm[0];
			d_nsvv_d_svvy = nrm[1];
			d_nsvv_d_svvz = nrm[2];



	d_nsvv_d_rux += d_nsvv_d_nrmz * d_nrmz_d_rux;
	d_nsvv_d_ruy += d_nsvv_d_nrmz * d_nrmz_d_ruy;
	d_nsvv_d_ruz += d_nsvv_d_nrmz * d_nrmz_d_ruz;
	d_nsvv_d_rvx += d_nsvv_d_nrmz * d_nrmz_d_rvx;
	d_nsvv_d_rvy += d_nsvv_d_nrmz * d_nrmz_d_rvy;
	d_nsvv_d_rvz += d_nsvv_d_nrmz * d_nrmz_d_rvz;
	d_nsvv_d_rux += d_nsvv_d_nrmy * d_nrmy_d_rux;
	d_nsvv_d_ruy += d_nsvv_d_nrmy * d_nrmy_d_ruy;
	d_nsvv_d_ruz += d_nsvv_d_nrmy * d_nrmy_d_ruz;
	d_nsvv_d_rvx += d_nsvv_d_nrmy * d_nrmy_d_rvx;
	d_nsvv_d_rvy += d_nsvv_d_nrmy * d_nrmy_d_rvy;
	d_nsvv_d_rvz += d_nsvv_d_nrmy * d_nrmy_d_rvz;
	d_nsvv_d_rux += d_nsvv_d_nrmx * d_nrmx_d_rux;
	d_nsvv_d_ruy += d_nsvv_d_nrmx * d_nrmx_d_ruy;
	d_nsvv_d_ruz += d_nsvv_d_nrmx * d_nrmx_d_ruz;
	d_nsvv_d_rvx += d_nsvv_d_nrmx * d_nrmx_d_rvx;
	d_nsvv_d_rvy += d_nsvv_d_nrmx * d_nrmx_d_rvy;
	d_nsvv_d_rvz += d_nsvv_d_nrmx * d_nrmx_d_rvz;
	d_nsuv_d_rux += d_nsuv_d_nrmz * d_nrmz_d_rux;
	d_nsuv_d_ruy += d_nsuv_d_nrmz * d_nrmz_d_ruy;
	d_nsuv_d_ruz += d_nsuv_d_nrmz * d_nrmz_d_ruz;
	d_nsuv_d_rvx += d_nsuv_d_nrmz * d_nrmz_d_rvx;
	d_nsuv_d_rvy += d_nsuv_d_nrmz * d_nrmz_d_rvy;
	d_nsuv_d_rvz += d_nsuv_d_nrmz * d_nrmz_d_rvz;
	d_nsuv_d_rux += d_nsuv_d_nrmy * d_nrmy_d_rux;
	d_nsuv_d_ruy += d_nsuv_d_nrmy * d_nrmy_d_ruy;
	d_nsuv_d_ruz += d_nsuv_d_nrmy * d_nrmy_d_ruz;
	d_nsuv_d_rvx += d_nsuv_d_nrmy * d_nrmy_d_rvx;
	d_nsuv_d_rvy += d_nsuv_d_nrmy * d_nrmy_d_rvy;
	d_nsuv_d_rvz += d_nsuv_d_nrmy * d_nrmy_d_rvz;
	d_nsuv_d_rux += d_nsuv_d_nrmx * d_nrmx_d_rux;
	d_nsuv_d_ruy += d_nsuv_d_nrmx * d_nrmx_d_ruy;
	d_nsuv_d_ruz += d_nsuv_d_nrmx * d_nrmx_d_ruz;
	d_nsuv_d_rvx += d_nsuv_d_nrmx * d_nrmx_d_rvx;
	d_nsuv_d_rvy += d_nsuv_d_nrmx * d_nrmx_d_rvy;
	d_nsuv_d_rvz += d_nsuv_d_nrmx * d_nrmx_d_rvz;
	d_nsuu_d_rux += d_nsuu_d_nrmz * d_nrmz_d_rux;
	d_nsuu_d_ruy += d_nsuu_d_nrmz * d_nrmz_d_ruy;
	d_nsuu_d_ruz += d_nsuu_d_nrmz * d_nrmz_d_ruz;
	d_nsuu_d_rvx += d_nsuu_d_nrmz * d_nrmz_d_rvx;
	d_nsuu_d_rvy += d_nsuu_d_nrmz * d_nrmz_d_rvy;
	d_nsuu_d_rvz += d_nsuu_d_nrmz * d_nrmz_d_rvz;
	d_nsuu_d_rux += d_nsuu_d_nrmy * d_nrmy_d_rux;
	d_nsuu_d_ruy += d_nsuu_d_nrmy * d_nrmy_d_ruy;
	d_nsuu_d_ruz += d_nsuu_d_nrmy * d_nrmy_d_ruz;
	d_nsuu_d_rvx += d_nsuu_d_nrmy * d_nrmy_d_rvx;
	d_nsuu_d_rvy += d_nsuu_d_nrmy * d_nrmy_d_rvy;
	d_nsuu_d_rvz += d_nsuu_d_nrmy * d_nrmy_d_rvz;
	d_nsuu_d_rux += d_nsuu_d_nrmx * d_nrmx_d_rux;
	d_nsuu_d_ruy += d_nsuu_d_nrmx * d_nrmx_d_ruy;
	d_nsuu_d_ruz += d_nsuu_d_nrmx * d_nrmx_d_ruz;
	d_nsuu_d_rvx += d_nsuu_d_nrmx * d_nrmx_d_rvx;
	d_nsuu_d_rvy += d_nsuu_d_nrmx * d_nrmx_d_rvy;
	d_nsuu_d_rvz += d_nsuu_d_nrmx * d_nrmx_d_rvz;
	d_g_d_rvx += d_g_d_RvRv * d_RvRv_d_rvx;
	d_g_d_rvy += d_g_d_RvRv * d_RvRv_d_rvy;
	d_g_d_rvz += d_g_d_RvRv * d_RvRv_d_rvz;
	d_g_d_rux += d_g_d_RuRv * d_RuRv_d_rux;
	d_g_d_ruy += d_g_d_RuRv * d_RuRv_d_ruy;
	d_g_d_ruz += d_g_d_RuRv * d_RuRv_d_ruz;
	d_g_d_rvx += d_g_d_RuRv * d_RuRv_d_rvx;
	d_g_d_rvy += d_g_d_RuRv * d_RuRv_d_rvy;
	d_g_d_rvz += d_g_d_RuRv * d_RuRv_d_rvz;
	d_g_d_rux += d_g_d_RuRu * d_RuRu_d_rux;
	d_g_d_ruy += d_g_d_RuRu * d_RuRu_d_ruy;
	d_g_d_ruz += d_g_d_RuRu * d_RuRu_d_ruz;
	d_d_d_rux += d_d_d_RuRv * d_RuRv_d_rux;
	d_d_d_ruy += d_d_d_RuRv * d_RuRv_d_ruy;
	d_d_d_ruz += d_d_d_RuRv * d_RuRv_d_ruz;
	d_d_d_rvx += d_d_d_RuRv * d_RuRv_d_rvx;
	d_d_d_rvy += d_d_d_RuRv * d_RuRv_d_rvy;
	d_d_d_rvz += d_d_d_RuRv * d_RuRv_d_rvz;
	d_d_d_rux += d_d_d_nsuv * d_nsuv_d_rux;
	d_d_d_ruy += d_d_d_nsuv * d_nsuv_d_ruy;
	d_d_d_ruz += d_d_d_nsuv * d_nsuv_d_ruz;
	d_d_d_rvx += d_d_d_nsuv * d_nsuv_d_rvx;
	d_d_d_rvy += d_d_d_nsuv * d_nsuv_d_rvy;
	d_d_d_rvz += d_d_d_nsuv * d_nsuv_d_rvz;
	d_d_d_suvx += d_d_d_nsuv * d_nsuv_d_suvx;
	d_d_d_suvy += d_d_d_nsuv * d_nsuv_d_suvy;
	d_d_d_suvz += d_d_d_nsuv * d_nsuv_d_suvz;
	d_d_d_rux += d_d_d_RuRu * d_RuRu_d_rux;
	d_d_d_ruy += d_d_d_RuRu * d_RuRu_d_ruy;
	d_d_d_ruz += d_d_d_RuRu * d_RuRu_d_ruz;
	d_d_d_rux += d_d_d_nsvv * d_nsvv_d_rux;
	d_d_d_ruy += d_d_d_nsvv * d_nsvv_d_ruy;
	d_d_d_ruz += d_d_d_nsvv * d_nsvv_d_ruz;
	d_d_d_rvx += d_d_d_nsvv * d_nsvv_d_rvx;
	d_d_d_rvy += d_d_d_nsvv * d_nsvv_d_rvy;
	d_d_d_rvz += d_d_d_nsvv * d_nsvv_d_rvz;
	d_d_d_svvx += d_d_d_nsvv * d_nsvv_d_svvx;
	d_d_d_svvy += d_d_d_nsvv * d_nsvv_d_svvy;
	d_d_d_svvz += d_d_d_nsvv * d_nsvv_d_svvz;
	d_d_d_rux += d_d_d_g * d_g_d_rux;
	d_d_d_ruy += d_d_d_g * d_g_d_ruy;
	d_d_d_ruz += d_d_d_g * d_g_d_ruz;
	d_d_d_rvx += d_d_d_g * d_g_d_rvx;
	d_d_d_rvy += d_d_d_g * d_g_d_rvy;
	d_d_d_rvz += d_d_d_g * d_g_d_rvz;
	d_c_d_rux += d_c_d_RuRv * d_RuRv_d_rux;
	d_c_d_ruy += d_c_d_RuRv * d_RuRv_d_ruy;
	d_c_d_ruz += d_c_d_RuRv * d_RuRv_d_ruz;
	d_c_d_rvx += d_c_d_RuRv * d_RuRv_d_rvx;
	d_c_d_rvy += d_c_d_RuRv * d_RuRv_d_rvy;
	d_c_d_rvz += d_c_d_RuRv * d_RuRv_d_rvz;
	d_c_d_rux += d_c_d_nsuv * d_nsuv_d_rux;
	d_c_d_ruy += d_c_d_nsuv * d_nsuv_d_ruy;
	d_c_d_ruz += d_c_d_nsuv * d_nsuv_d_ruz;
	d_c_d_rvx += d_c_d_nsuv * d_nsuv_d_rvx;
	d_c_d_rvy += d_c_d_nsuv * d_nsuv_d_rvy;
	d_c_d_rvz += d_c_d_nsuv * d_nsuv_d_rvz;
	d_c_d_suvx += d_c_d_nsuv * d_nsuv_d_suvx;
	d_c_d_suvy += d_c_d_nsuv * d_nsuv_d_suvy;
	d_c_d_suvz += d_c_d_nsuv * d_nsuv_d_suvz;
	d_c_d_rux += d_c_d_RuRu * d_RuRu_d_rux;
	d_c_d_ruy += d_c_d_RuRu * d_RuRu_d_ruy;
	d_c_d_ruz += d_c_d_RuRu * d_RuRu_d_ruz;
	d_c_d_rux += d_c_d_nsuu * d_nsuu_d_rux;
	d_c_d_ruy += d_c_d_nsuu * d_nsuu_d_ruy;
	d_c_d_ruz += d_c_d_nsuu * d_nsuu_d_ruz;
	d_c_d_rvx += d_c_d_nsuu * d_nsuu_d_rvx;
	d_c_d_rvy += d_c_d_nsuu * d_nsuu_d_rvy;
	d_c_d_rvz += d_c_d_nsuu * d_nsuu_d_rvz;
	d_c_d_suux += d_c_d_nsuu * d_nsuu_d_suux;
	d_c_d_suuy += d_c_d_nsuu * d_nsuu_d_suuy;
	d_c_d_suuz += d_c_d_nsuu * d_nsuu_d_suuz;
	d_c_d_rux += d_c_d_g * d_g_d_rux;
	d_c_d_ruy += d_c_d_g * d_g_d_ruy;
	d_c_d_ruz += d_c_d_g * d_g_d_ruz;
	d_c_d_rvx += d_c_d_g * d_g_d_rvx;
	d_c_d_rvy += d_c_d_g * d_g_d_rvy;
	d_c_d_rvz += d_c_d_g * d_g_d_rvz;
	d_b_d_rux += d_b_d_RuRv * d_RuRv_d_rux;
	d_b_d_ruy += d_b_d_RuRv * d_RuRv_d_ruy;
	d_b_d_ruz += d_b_d_RuRv * d_RuRv_d_ruz;
	d_b_d_rvx += d_b_d_RuRv * d_RuRv_d_rvx;
	d_b_d_rvy += d_b_d_RuRv * d_RuRv_d_rvy;
	d_b_d_rvz += d_b_d_RuRv * d_RuRv_d_rvz;
	d_b_d_rux += d_b_d_nsuv * d_nsuv_d_rux;
	d_b_d_ruy += d_b_d_nsuv * d_nsuv_d_ruy;
	d_b_d_ruz += d_b_d_nsuv * d_nsuv_d_ruz;
	d_b_d_rvx += d_b_d_nsuv * d_nsuv_d_rvx;
	d_b_d_rvy += d_b_d_nsuv * d_nsuv_d_rvy;
	d_b_d_rvz += d_b_d_nsuv * d_nsuv_d_rvz;
	d_b_d_suvx += d_b_d_nsuv * d_nsuv_d_suvx;
	d_b_d_suvy += d_b_d_nsuv * d_nsuv_d_suvy;
	d_b_d_suvz += d_b_d_nsuv * d_nsuv_d_suvz;
	d_b_d_rvx += d_b_d_RvRv * d_RvRv_d_rvx;
	d_b_d_rvy += d_b_d_RvRv * d_RvRv_d_rvy;
	d_b_d_rvz += d_b_d_RvRv * d_RvRv_d_rvz;
	d_b_d_rux += d_b_d_nsvv * d_nsvv_d_rux;
	d_b_d_ruy += d_b_d_nsvv * d_nsvv_d_ruy;
	d_b_d_ruz += d_b_d_nsvv * d_nsvv_d_ruz;
	d_b_d_rvx += d_b_d_nsvv * d_nsvv_d_rvx;
	d_b_d_rvy += d_b_d_nsvv * d_nsvv_d_rvy;
	d_b_d_rvz += d_b_d_nsvv * d_nsvv_d_rvz;
	d_b_d_svvx += d_b_d_nsvv * d_nsvv_d_svvx;
	d_b_d_svvy += d_b_d_nsvv * d_nsvv_d_svvy;
	d_b_d_svvz += d_b_d_nsvv * d_nsvv_d_svvz;
	d_b_d_rux += d_b_d_g * d_g_d_rux;
	d_b_d_ruy += d_b_d_g * d_g_d_ruy;
	d_b_d_ruz += d_b_d_g * d_g_d_ruz;
	d_b_d_rvx += d_b_d_g * d_g_d_rvx;
	d_b_d_rvy += d_b_d_g * d_g_d_rvy;
	d_b_d_rvz += d_b_d_g * d_g_d_rvz;
	d_a_d_rux += d_a_d_RuRv * d_RuRv_d_rux;
	d_a_d_ruy += d_a_d_RuRv * d_RuRv_d_ruy;
	d_a_d_ruz += d_a_d_RuRv * d_RuRv_d_ruz;
	d_a_d_rvx += d_a_d_RuRv * d_RuRv_d_rvx;
	d_a_d_rvy += d_a_d_RuRv * d_RuRv_d_rvy;
	d_a_d_rvz += d_a_d_RuRv * d_RuRv_d_rvz;
	d_a_d_rux += d_a_d_nsuv * d_nsuv_d_rux;
	d_a_d_ruy += d_a_d_nsuv * d_nsuv_d_ruy;
	d_a_d_ruz += d_a_d_nsuv * d_nsuv_d_ruz;
	d_a_d_rvx += d_a_d_nsuv * d_nsuv_d_rvx;
	d_a_d_rvy += d_a_d_nsuv * d_nsuv_d_rvy;
	d_a_d_rvz += d_a_d_nsuv * d_nsuv_d_rvz;
	d_a_d_suvx += d_a_d_nsuv * d_nsuv_d_suvx;
	d_a_d_suvy += d_a_d_nsuv * d_nsuv_d_suvy;
	d_a_d_suvz += d_a_d_nsuv * d_nsuv_d_suvz;
	d_a_d_rvx += d_a_d_RvRv * d_RvRv_d_rvx;
	d_a_d_rvy += d_a_d_RvRv * d_RvRv_d_rvy;
	d_a_d_rvz += d_a_d_RvRv * d_RvRv_d_rvz;
	d_a_d_rux += d_a_d_nsuu * d_nsuu_d_rux;
	d_a_d_ruy += d_a_d_nsuu * d_nsuu_d_ruy;
	d_a_d_ruz += d_a_d_nsuu * d_nsuu_d_ruz;
	d_a_d_rvx += d_a_d_nsuu * d_nsuu_d_rvx;
	d_a_d_rvy += d_a_d_nsuu * d_nsuu_d_rvy;
	d_a_d_rvz += d_a_d_nsuu * d_nsuu_d_rvz;
	d_a_d_suux += d_a_d_nsuu * d_nsuu_d_suux;
	d_a_d_suuy += d_a_d_nsuu * d_nsuu_d_suuy;
	d_a_d_suuz += d_a_d_nsuu * d_nsuu_d_suuz;
	d_a_d_rux += d_a_d_g * d_g_d_rux;
	d_a_d_ruy += d_a_d_g * d_g_d_ruy;
	d_a_d_ruz += d_a_d_g * d_g_d_ruz;
	d_a_d_rvx += d_a_d_g * d_g_d_rvx;
	d_a_d_rvy += d_a_d_g * d_g_d_rvy;
	d_a_d_rvz += d_a_d_g * d_g_d_rvz;
	d_c1_d_rux += d_c1_d_d * d_d_d_rux;
	d_c1_d_ruy += d_c1_d_d * d_d_d_ruy;
	d_c1_d_ruz += d_c1_d_d * d_d_d_ruz;
	d_c1_d_rvx += d_c1_d_d * d_d_d_rvx;
	d_c1_d_rvy += d_c1_d_d * d_d_d_rvy;
	d_c1_d_rvz += d_c1_d_d * d_d_d_rvz;
	d_c1_d_suvx += d_c1_d_d * d_d_d_suvx;
	d_c1_d_suvy += d_c1_d_d * d_d_d_suvy;
	d_c1_d_suvz += d_c1_d_d * d_d_d_suvz;
	d_c1_d_svvx += d_c1_d_d * d_d_d_svvx;
	d_c1_d_svvy += d_c1_d_d * d_d_d_svvy;
	d_c1_d_svvz += d_c1_d_d * d_d_d_svvz;
	d_c1_d_rux += d_c1_d_c * d_c_d_rux;
	d_c1_d_ruy += d_c1_d_c * d_c_d_ruy;
	d_c1_d_ruz += d_c1_d_c * d_c_d_ruz;
	d_c1_d_rvx += d_c1_d_c * d_c_d_rvx;
	d_c1_d_rvy += d_c1_d_c * d_c_d_rvy;
	d_c1_d_rvz += d_c1_d_c * d_c_d_rvz;
	d_c1_d_suux += d_c1_d_c * d_c_d_suux;
	d_c1_d_suuy += d_c1_d_c * d_c_d_suuy;
	d_c1_d_suuz += d_c1_d_c * d_c_d_suuz;
	d_c1_d_suvx += d_c1_d_c * d_c_d_suvx;
	d_c1_d_suvy += d_c1_d_c * d_c_d_suvy;
	d_c1_d_suvz += d_c1_d_c * d_c_d_suvz;
	d_c1_d_rux += d_c1_d_b * d_b_d_rux;
	d_c1_d_ruy += d_c1_d_b * d_b_d_ruy;
	d_c1_d_ruz += d_c1_d_b * d_b_d_ruz;
	d_c1_d_rvx += d_c1_d_b * d_b_d_rvx;
	d_c1_d_rvy += d_c1_d_b * d_b_d_rvy;
	d_c1_d_rvz += d_c1_d_b * d_b_d_rvz;
	d_c1_d_suvx += d_c1_d_b * d_b_d_suvx;
	d_c1_d_suvy += d_c1_d_b * d_b_d_suvy;
	d_c1_d_suvz += d_c1_d_b * d_b_d_suvz;
	d_c1_d_svvx += d_c1_d_b * d_b_d_svvx;
	d_c1_d_svvy += d_c1_d_b * d_b_d_svvy;
	d_c1_d_svvz += d_c1_d_b * d_b_d_svvz;
	d_c1_d_rux += d_c1_d_a * d_a_d_rux;
	d_c1_d_ruy += d_c1_d_a * d_a_d_ruy;
	d_c1_d_ruz += d_c1_d_a * d_a_d_ruz;
	d_c1_d_rvx += d_c1_d_a * d_a_d_rvx;
	d_c1_d_rvy += d_c1_d_a * d_a_d_rvy;
	d_c1_d_rvz += d_c1_d_a * d_a_d_rvz;
	d_c1_d_suux += d_c1_d_a * d_a_d_suux;
	d_c1_d_suuy += d_c1_d_a * d_a_d_suuy;
	d_c1_d_suuz += d_c1_d_a * d_a_d_suuz;
	d_c1_d_suvx += d_c1_d_a * d_a_d_suvx;
	d_c1_d_suvy += d_c1_d_a * d_a_d_suvy;
	d_c1_d_suvz += d_c1_d_a * d_a_d_suvz;
	d_c2_d_rux += d_c2_d_d * d_d_d_rux;
	d_c2_d_ruy += d_c2_d_d * d_d_d_ruy;
	d_c2_d_ruz += d_c2_d_d * d_d_d_ruz;
	d_c2_d_rvx += d_c2_d_d * d_d_d_rvx;
	d_c2_d_rvy += d_c2_d_d * d_d_d_rvy;
	d_c2_d_rvz += d_c2_d_d * d_d_d_rvz;
	d_c2_d_suvx += d_c2_d_d * d_d_d_suvx;
	d_c2_d_suvy += d_c2_d_d * d_d_d_suvy;
	d_c2_d_suvz += d_c2_d_d * d_d_d_suvz;
	d_c2_d_svvx += d_c2_d_d * d_d_d_svvx;
	d_c2_d_svvy += d_c2_d_d * d_d_d_svvy;
	d_c2_d_svvz += d_c2_d_d * d_d_d_svvz;
	d_c2_d_rux += d_c2_d_c * d_c_d_rux;
	d_c2_d_ruy += d_c2_d_c * d_c_d_ruy;
	d_c2_d_ruz += d_c2_d_c * d_c_d_ruz;
	d_c2_d_rvx += d_c2_d_c * d_c_d_rvx;
	d_c2_d_rvy += d_c2_d_c * d_c_d_rvy;
	d_c2_d_rvz += d_c2_d_c * d_c_d_rvz;
	d_c2_d_suux += d_c2_d_c * d_c_d_suux;
	d_c2_d_suuy += d_c2_d_c * d_c_d_suuy;
	d_c2_d_suuz += d_c2_d_c * d_c_d_suuz;
	d_c2_d_suvx += d_c2_d_c * d_c_d_suvx;
	d_c2_d_suvy += d_c2_d_c * d_c_d_suvy;
	d_c2_d_suvz += d_c2_d_c * d_c_d_suvz;
	d_c2_d_rux += d_c2_d_b * d_b_d_rux;
	d_c2_d_ruy += d_c2_d_b * d_b_d_ruy;
	d_c2_d_ruz += d_c2_d_b * d_b_d_ruz;
	d_c2_d_rvx += d_c2_d_b * d_b_d_rvx;
	d_c2_d_rvy += d_c2_d_b * d_b_d_rvy;
	d_c2_d_rvz += d_c2_d_b * d_b_d_rvz;
	d_c2_d_suvx += d_c2_d_b * d_b_d_suvx;
	d_c2_d_suvy += d_c2_d_b * d_b_d_suvy;
	d_c2_d_suvz += d_c2_d_b * d_b_d_suvz;
	d_c2_d_svvx += d_c2_d_b * d_b_d_svvx;
	d_c2_d_svvy += d_c2_d_b * d_b_d_svvy;
	d_c2_d_svvz += d_c2_d_b * d_b_d_svvz;
	d_c2_d_rux += d_c2_d_a * d_a_d_rux;
	d_c2_d_ruy += d_c2_d_a * d_a_d_ruy;
	d_c2_d_ruz += d_c2_d_a * d_a_d_ruz;
	d_c2_d_rvx += d_c2_d_a * d_a_d_rvx;
	d_c2_d_rvy += d_c2_d_a * d_a_d_rvy;
	d_c2_d_rvz += d_c2_d_a * d_a_d_rvz;
	d_c2_d_suux += d_c2_d_a * d_a_d_suux;
	d_c2_d_suuy += d_c2_d_a * d_a_d_suuy;
	d_c2_d_suuz += d_c2_d_a * d_a_d_suuz;
	d_c2_d_suvx += d_c2_d_a * d_a_d_suvx;
	d_c2_d_suvy += d_c2_d_a * d_a_d_suvy;
	d_c2_d_suvz += d_c2_d_a * d_a_d_suvz;
	d_e_d_rux += d_e_d_c2 * d_c2_d_rux;
	d_e_d_ruy += d_e_d_c2 * d_c2_d_ruy;
	d_e_d_ruz += d_e_d_c2 * d_c2_d_ruz;
	d_e_d_rvx += d_e_d_c2 * d_c2_d_rvx;
	d_e_d_rvy += d_e_d_c2 * d_c2_d_rvy;
	d_e_d_rvz += d_e_d_c2 * d_c2_d_rvz;
	d_e_d_suux += d_e_d_c2 * d_c2_d_suux;
	d_e_d_suuy += d_e_d_c2 * d_c2_d_suuy;
	d_e_d_suuz += d_e_d_c2 * d_c2_d_suuz;
	d_e_d_suvx += d_e_d_c2 * d_c2_d_suvx;
	d_e_d_suvy += d_e_d_c2 * d_c2_d_suvy;
	d_e_d_suvz += d_e_d_c2 * d_c2_d_suvz;
	d_e_d_svvx += d_e_d_c2 * d_c2_d_svvx;
	d_e_d_svvy += d_e_d_c2 * d_c2_d_svvy;
	d_e_d_svvz += d_e_d_c2 * d_c2_d_svvz;
	d_e_d_rux += d_e_d_c1 * d_c1_d_rux;
	d_e_d_ruy += d_e_d_c1 * d_c1_d_ruy;
	d_e_d_ruz += d_e_d_c1 * d_c1_d_ruz;
	d_e_d_rvx += d_e_d_c1 * d_c1_d_rvx;
	d_e_d_rvy += d_e_d_c1 * d_c1_d_rvy;
	d_e_d_rvz += d_e_d_c1 * d_c1_d_rvz;
	d_e_d_suux += d_e_d_c1 * d_c1_d_suux;
	d_e_d_suuy += d_e_d_c1 * d_c1_d_suuy;
	d_e_d_suuz += d_e_d_c1 * d_c1_d_suuz;
	d_e_d_suvx += d_e_d_c1 * d_c1_d_suvx;
	d_e_d_suvy += d_e_d_c1 * d_c1_d_suvy;
	d_e_d_suvz += d_e_d_c1 * d_c1_d_suvz;
	d_e_d_svvx += d_e_d_c1 * d_c1_d_svvx;
	d_e_d_svvy += d_e_d_c1 * d_c1_d_svvy;
	d_e_d_svvz += d_e_d_c1 * d_c1_d_svvz;
	d_e_d_rux += d_e_d_g * d_g_d_rux;
	d_e_d_ruy += d_e_d_g * d_g_d_ruy;
	d_e_d_ruz += d_e_d_g * d_g_d_ruz;
	d_e_d_rvx += d_e_d_g * d_g_d_rvx;
	d_e_d_rvy += d_e_d_g * d_g_d_rvy;
	d_e_d_rvz += d_e_d_g * d_g_d_rvz;
			

			// dE/du 			

			double dedr[3*np];
			memset( dedr, 0, sizeof(double) * 3 *np );
			for( int p = 0; p < np; p++ )
			for( int y = 0; y < 12; y++ )
			{
				gr[3*cp[p]+0] += d_e_d_rux *  ceff_map_du[y] * theMap[y*ncoords_base+p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_ruy *  ceff_map_du[y] * theMap[y*ncoords_base+p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_ruz *  ceff_map_du[y] * theMap[y*ncoords_base+p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_rvx *  ceff_map_dv[y] * theMap[y*ncoords_base+p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_rvy *  ceff_map_dv[y] * theMap[y*ncoords_base+p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_rvz *  ceff_map_dv[y] * theMap[y*ncoords_base+p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_suux * ceff_map_duu[y] * theMap[y*ncoords_base+p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_suuy * ceff_map_duu[y] * theMap[y*ncoords_base+p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_suuz * ceff_map_duu[y] * theMap[y*ncoords_base+p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_suvx * ceff_map_duv[y] * theMap[y*ncoords_base+p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_suvy * ceff_map_duv[y] * theMap[y*ncoords_base+p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_suvz * ceff_map_duv[y] * theMap[y*ncoords_base+p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_svvx * ceff_map_dvv[y] * theMap[y*ncoords_base+p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_svvy * ceff_map_dvv[y] * theMap[y*ncoords_base+p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_svvz * ceff_map_dvv[y] * theMap[y*ncoords_base+p] * alpha_z; 
				
				dedr[3*p+0] += d_e_d_rux *  ceff_map_du[y] * theMap[y*ncoords_base+p];
				dedr[3*p+1] += d_e_d_ruy *  ceff_map_du[y] * theMap[y*ncoords_base+p];
				dedr[3*p+2] += d_e_d_ruz *  ceff_map_du[y] * theMap[y*ncoords_base+p];
				                                                                     
				dedr[3*p+0] += d_e_d_rvx *  ceff_map_dv[y] * theMap[y*ncoords_base+p];
				dedr[3*p+1] += d_e_d_rvy *  ceff_map_dv[y] * theMap[y*ncoords_base+p];
				dedr[3*p+2] += d_e_d_rvz *  ceff_map_dv[y] * theMap[y*ncoords_base+p];
				                                                                       
				dedr[3*p+0] += d_e_d_suux * ceff_map_duu[y] * theMap[y*ncoords_base+p]; 
				dedr[3*p+1] += d_e_d_suuy * ceff_map_duu[y] * theMap[y*ncoords_base+p]; 
				dedr[3*p+2] += d_e_d_suuz * ceff_map_duu[y] * theMap[y*ncoords_base+p]; 
				                                                                      
				dedr[3*p+0] += d_e_d_suvx * ceff_map_duv[y] * theMap[y*ncoords_base+p]; 
				dedr[3*p+1] += d_e_d_suvy * ceff_map_duv[y] * theMap[y*ncoords_base+p]; 
				dedr[3*p+2] += d_e_d_suvz * ceff_map_duv[y] * theMap[y*ncoords_base+p]; 
				                                                                      
				dedr[3*p+0] += d_e_d_svvx * ceff_map_dvv[y] * theMap[y*ncoords_base+p]; 
				dedr[3*p+1] += d_e_d_svvy * ceff_map_dvv[y] * theMap[y*ncoords_base+p]; 
				dedr[3*p+2] += d_e_d_svvz * ceff_map_dvv[y] * theMap[y*ncoords_base+p]; 
			}

/*	
			for( int p = 0; p < np; p++ )
			{
				gr[3*cp[p]+0] += d_e_d_rux * theIrregularFormulas[frm].r_u[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_ruy * theIrregularFormulas[frm].r_u[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_ruz * theIrregularFormulas[frm].r_u[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_rvx * theIrregularFormulas[frm].r_v[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_rvy * theIrregularFormulas[frm].r_v[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_rvz * theIrregularFormulas[frm].r_v[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_suux * theIrregularFormulas[frm].r_uu[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_suuy * theIrregularFormulas[frm].r_uu[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_suuz * theIrregularFormulas[frm].r_uu[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_suvx * theIrregularFormulas[frm].r_uv[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_suvy * theIrregularFormulas[frm].r_uv[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_suvz * theIrregularFormulas[frm].r_uv[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_svvx * theIrregularFormulas[frm].r_vv[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_svvy * theIrregularFormulas[frm].r_vv[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_svvz * theIrregularFormulas[frm].r_vv[p] * alpha_z; 
				
				dedr[3*p+0] += d_e_d_rux * theIrregularFormulas[frm].r_u[p]; 
				dedr[3*p+1] += d_e_d_ruy * theIrregularFormulas[frm].r_u[p]; 
				dedr[3*p+2] += d_e_d_ruz * theIrregularFormulas[frm].r_u[p]; 
				
				dedr[3*p+0] += d_e_d_rvx * theIrregularFormulas[frm].r_v[p]; 
				dedr[3*p+1] += d_e_d_rvy * theIrregularFormulas[frm].r_v[p]; 
				dedr[3*p+2] += d_e_d_rvz * theIrregularFormulas[frm].r_v[p]; 
				
				dedr[3*p+0] += d_e_d_suux * theIrregularFormulas[frm].r_uu[p]; 
				dedr[3*p+1] += d_e_d_suuy * theIrregularFormulas[frm].r_uu[p]; 
				dedr[3*p+2] += d_e_d_suuz * theIrregularFormulas[frm].r_uu[p]; 
				
				dedr[3*p+0] += d_e_d_suvx * theIrregularFormulas[frm].r_uv[p]; 
				dedr[3*p+1] += d_e_d_suvy * theIrregularFormulas[frm].r_uv[p]; 
				dedr[3*p+2] += d_e_d_suvz * theIrregularFormulas[frm].r_uv[p]; 
				
				dedr[3*p+0] += d_e_d_svvx * theIrregularFormulas[frm].r_vv[p]; 
				dedr[3*p+1] += d_e_d_svvy * theIrregularFormulas[frm].r_vv[p]; 
				dedr[3*p+2] += d_e_d_svvz * theIrregularFormulas[frm].r_vv[p]; 


			}
*/
			for( int p = 0; p < np; p++ )
			{
				gr[3*nv+0] += dedr[3*p+0] * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				gr[3*nv+1] += dedr[3*p+1] * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				gr[3*nv+2] += dedr[3*p+2] * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
			}
		}	
	}

}

#if 0
void surface::strainGrad( double *r, double *gr, int vary_g  )
{
	// computes the contribution of the lambda-strain, including g values (if requested) to g.
	double *strains = (double *)malloc( sizeof(double) * nt );
	double *d_e_d_g = (double *)malloc( sizeof(double) * nt );
	double *d_e_d_g0 = (double *)malloc( sizeof(double) * nt );

	if( vary_g )
		memcpy( f_g0, r +3 * nv + 3, sizeof(double) * nt );

	normalize_g0();
	
	memset( d_e_d_g,  0, sizeof(double) * nt );
	memset( d_e_d_g0, 0, sizeof(double) * nt );

	for( int t = 0; t < nt; t++ )
		strains[t] = fstrain(t,r);
	
	for( int t = 0; t < nt; t++ )
	{
		// d_e_d_g:
	
		for( int fx = 0; fx < 3; fx++ )
		{ 
			int fn = fneighbors[3*t+fx];

			// derivative wrt total face g (which is average of gs).	
			d_e_d_g[t]  += k_strain * (strains[t] - strains[fn]) / (f_g0[t]);
			// derivative wrt normalized g0s:
			d_e_d_g0[t] += k_strain * (strains[t]-strains[fn]) * (-strains[t]/f_g0[t]-1.0/f_g0[t]);
		}
	}	


	// now compute vertex derivatives.
	// regular faces first.
	
	double alpha_x = r[3*nv];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];
	if( !theFormulas )
		generatePlan();

	double r_val = 0;
	double e = 0;
	double area = 0;
	double wgt = 0;
	double dudv = 0.5;

	for( int f = 0; f < nf_faces; f++ )
	{
		double p_face_area = 0;
		double g0 = 2*f_g0[f];

		int frm = f*nf_g_q_p;
		int t = theFormulas[frm].tri;
	
		for( int p = 0; p < nf_g_q_p; p++ )
		{
			double A = 0, A0 = 0;
			int frm = f*nf_g_q_p+p;
			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
			double nrm[3]={0,0,0}; 

			int *cp = theFormulas[f*nf_g_q_p+p].cp;
			int np = theFormulas[f*nf_g_q_p+p].ncoor;
			double dedr[3*np];
			memset( dedr, 0, sizeof(double) * 3 *np );

			for( int p = 0; p < np; p++ )
			{
				R[0] += theFormulas[frm].r_w[p] * alpha_x * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				R[1] += theFormulas[frm].r_w[p] * alpha_y * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				R[2] += theFormulas[frm].r_w[p] * alpha_z * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				Ru[0] += theFormulas[frm].r_u[p] * alpha_x * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				Ru[1] += theFormulas[frm].r_u[p] * alpha_y * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				Ru[2] += theFormulas[frm].r_u[p] * alpha_z * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				Rv[0] += theFormulas[frm].r_v[p] * alpha_x * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				Rv[1] += theFormulas[frm].r_v[p] * alpha_y * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				Rv[2] += theFormulas[frm].r_v[p] * alpha_z * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuu[0] += theFormulas[frm].r_uu[p] * alpha_x * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuu[1] += theFormulas[frm].r_uu[p] * alpha_y * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuu[2] += theFormulas[frm].r_uu[p] * alpha_z * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSuv[0] += theFormulas[frm].r_uv[p] * alpha_x * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSuv[1] += theFormulas[frm].r_uv[p] * alpha_y * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSuv[2] += theFormulas[frm].r_uv[p] * alpha_z * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				tSvv[0] += theFormulas[frm].r_vv[p] * alpha_x * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				tSvv[1] += theFormulas[frm].r_vv[p] * alpha_y * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				tSvv[2] += theFormulas[frm].r_vv[p] * alpha_z * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			}
		



		cross( Ru, Rv, nrm );
		normalize(nrm);

		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		double gv = sqrt(RuRu*RvRv-RuRv*RuRv);

double d_nrmz_d_rux=0,d_nrmz_d_ruy=0,d_nrmz_d_ruz=0,d_nrmz_d_rvx=0,d_nrmz_d_rvy=0,d_nrmz_d_rvz=0,d_nrmy_d_rux=0,d_nrmy_d_ruy=0,d_nrmy_d_ruz=0,d_nrmy_d_rvx=0,d_nrmy_d_rvy=0,d_nrmy_d_rvz=0,d_nrmx_d_rux=0,d_nrmx_d_ruy=0,d_nrmx_d_ruz=0,d_nrmx_d_rvx=0,d_nrmx_d_rvy=0,d_nrmx_d_rvz=0,d_RuRu_d_rux=0,d_RuRu_d_ruy=0,d_RuRu_d_ruz=0,d_nsvv_d_rux=0,d_nsvv_d_ruy=0,d_nsvv_d_ruz=0,d_nsvv_d_rvx=0,d_nsvv_d_rvy=0,d_nsvv_d_rvz=0,d_nsvv_d_svvx=0,d_nsvv_d_svvy=0,d_nsvv_d_svvz=0,d_RuRv_d_rux=0,d_RuRv_d_ruy=0,d_RuRv_d_ruz=0,d_RuRv_d_rvx=0,d_RuRv_d_rvy=0,d_RuRv_d_rvz=0,d_nsuv_d_rux=0,d_nsuv_d_ruy=0,d_nsuv_d_ruz=0,d_nsuv_d_rvx=0,d_nsuv_d_rvy=0,d_nsuv_d_rvz=0,d_nsuv_d_suvx=0,d_nsuv_d_suvy=0,d_nsuv_d_suvz=0,d_RvRv_d_rvx=0,d_RvRv_d_rvy=0,d_RvRv_d_rvz=0,d_nsuu_d_rux=0,d_nsuu_d_ruy=0,d_nsuu_d_ruz=0,d_nsuu_d_rvx=0,d_nsuu_d_rvy=0,d_nsuu_d_rvz=0,d_nsuu_d_suux=0,d_nsuu_d_suuy=0,d_nsuu_d_suuz=0,d_g_d_rux=0,d_g_d_ruy=0,d_g_d_ruz=0,d_g_d_rvx=0,d_g_d_rvy=0,d_g_d_rvz=0,d_e_d_rux=0,d_e_d_ruy=0,d_e_d_ruz=0,d_e_d_rvx=0,d_e_d_rvy=0,d_e_d_rvz=0,d_e_d_suux=0,d_e_d_suuy=0,d_e_d_suuz=0,d_e_d_suvx=0,d_e_d_suvy=0,d_e_d_suvz=0,d_e_d_svvx=0,d_e_d_svvy=0,d_e_d_svvz=0,d_c2_d_rux=0,d_c2_d_ruy=0,d_c2_d_ruz=0,d_c2_d_rvx=0,d_c2_d_rvy=0,d_c2_d_rvz=0,d_c2_d_suux=0,d_c2_d_suuy=0,d_c2_d_suuz=0,d_c2_d_suvx=0,d_c2_d_suvy=0,d_c2_d_suvz=0,d_c2_d_svvx=0,d_c2_d_svvy=0,d_c2_d_svvz=0,d_d_d_rux=0,d_d_d_ruy=0,d_d_d_ruz=0,d_d_d_rvx=0,d_d_d_rvy=0,d_d_d_rvz=0,d_d_d_suvx=0,d_d_d_suvy=0,d_d_d_suvz=0,d_d_d_svvx=0,d_d_d_svvy=0,d_d_d_svvz=0,d_c_d_rux=0,d_c_d_ruy=0,d_c_d_ruz=0,d_c_d_rvx=0,d_c_d_rvy=0,d_c_d_rvz=0,d_c_d_suux=0,d_c_d_suuy=0,d_c_d_suuz=0,d_c_d_suvx=0,d_c_d_suvy=0,d_c_d_suvz=0,d_b_d_rux=0,d_b_d_ruy=0,d_b_d_ruz=0,d_b_d_rvx=0,d_b_d_rvy=0,d_b_d_rvz=0,d_b_d_suvx=0,d_b_d_suvy=0,d_b_d_suvz=0,d_b_d_svvx=0,d_b_d_svvy=0,d_b_d_svvz=0,d_a_d_rux=0,d_a_d_ruy=0,d_a_d_ruz=0,d_a_d_rvx=0,d_a_d_rvy=0,d_a_d_rvz=0,d_a_d_suux=0,d_a_d_suuy=0,d_a_d_suuz=0,d_a_d_suvx=0,d_a_d_suvy=0,d_a_d_suvz=0,d_c1_d_rux=0,d_c1_d_ruy=0,d_c1_d_ruz=0,d_c1_d_rvx=0,d_c1_d_rvy=0,d_c1_d_rvz=0,d_c1_d_suux=0,d_c1_d_suuy=0,d_c1_d_suuz=0,d_c1_d_suvx=0,d_c1_d_suvy=0,d_c1_d_suvz=0,d_c1_d_svvx=0,d_c1_d_svvy=0,d_c1_d_svvz=0,d_g_d_RvRv=0,d_g_d_RuRv=0,d_g_d_RuRu=0,d_nsvv_d_nrmz=0,d_nsvv_d_nrmy=0,d_nsvv_d_nrmx=0,d_nsuv_d_nrmz=0,d_nsuv_d_nrmy=0,d_nsuv_d_nrmx=0,d_nsuu_d_nrmz=0,d_nsuu_d_nrmy=0,d_nsuu_d_nrmx=0,d_d_d_RuRv=0,d_d_d_nsuv=0,d_d_d_RuRu=0,d_d_d_nsvv=0,d_d_d_g=0,d_c_d_RuRv=0,d_c_d_nsuv=0,d_c_d_RuRu=0,d_c_d_nsuu=0,d_c_d_g=0,d_b_d_RuRv=0,d_b_d_nsuv=0,d_b_d_RvRv=0,d_b_d_nsvv=0,d_b_d_g=0,d_a_d_RuRv=0,d_a_d_nsuv=0,d_a_d_RvRv=0,d_a_d_nsuu=0,d_a_d_g=0,d_e_d_c2=0,d_e_d_c1=0,d_c2_d_d=0,d_c2_d_c=0,d_c2_d_b=0,d_c2_d_a=0,d_c1_d_d=0,d_c1_d_c=0,d_c1_d_b=0,d_c1_d_a=0,junk;
			// the basic variables.

			d_g_d_RuRu = RvRv/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
			d_g_d_RuRv = -(RuRv/Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
			d_g_d_RvRv = RuRu/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  

			d_RuRu_d_rux = 2*Ru[0];
			d_RuRu_d_ruy = 2*Ru[1];
			d_RuRu_d_ruz = 2*Ru[2];
			
			d_RuRv_d_rux = Rv[0];
			d_RuRv_d_ruy = Rv[1];
			d_RuRv_d_ruz = Rv[2];
			
			d_RuRv_d_rvx = Ru[0];
			d_RuRv_d_rvy = Ru[1];
			d_RuRv_d_rvz = Ru[2];
			
			d_RvRv_d_rvx = 2*Rv[0];
			d_RvRv_d_rvy = 2*Rv[1];
			d_RvRv_d_rvz = 2*Rv[2];

	d_g_d_rvx += d_g_d_RvRv * d_RvRv_d_rvx;
	d_g_d_rvy += d_g_d_RvRv * d_RvRv_d_rvy;
	d_g_d_rvz += d_g_d_RvRv * d_RvRv_d_rvz;
	d_g_d_rux += d_g_d_RuRv * d_RuRv_d_rux;
	d_g_d_ruy += d_g_d_RuRv * d_RuRv_d_ruy;
	d_g_d_ruz += d_g_d_RuRv * d_RuRv_d_ruz;
	d_g_d_rvx += d_g_d_RuRv * d_RuRv_d_rvx;
	d_g_d_rvy += d_g_d_RuRv * d_RuRv_d_rvy;
	d_g_d_rvz += d_g_d_RuRv * d_RuRv_d_rvz;
	d_g_d_rux += d_g_d_RuRu * d_RuRu_d_rux;
	d_g_d_ruy += d_g_d_RuRu * d_RuRu_d_ruy;
	d_g_d_ruz += d_g_d_RuRu * d_RuRu_d_ruz;

	// KA A0 * ( A - A0) / A0 
	
	// need KA 

	
	
	d_e_d_rux += 0.5*d_e_d_g[f] * d_g_d_rux * theFormulas[frm].weight;
	d_e_d_ruy += 0.5*d_e_d_g[f] * d_g_d_ruy * theFormulas[frm].weight;
	d_e_d_ruz += 0.5*d_e_d_g[f] * d_g_d_ruz * theFormulas[frm].weight;
	d_e_d_rvx += 0.5*d_e_d_g[f] * d_g_d_rvx * theFormulas[frm].weight;
	d_e_d_rvy += 0.5*d_e_d_g[f] * d_g_d_rvy * theFormulas[frm].weight;
	d_e_d_rvz += 0.5*d_e_d_g[f] * d_g_d_rvz * theFormulas[frm].weight;
			
			d_e_d_g0[f] += 2 * KA * (1 - gv*gv / (g0*g0) ) * dudv * theFormulas[frm].weight; 

//			if( fabs(g/g0-1)>1e-5 )
//			{
//				printf("f: %d g: %le g0: %le\n", f, g, g0 );
//				exit(1);
//			}

			for( int p = 0; p < np; p++ )
			{
				gr[3*cp[p]+0] += d_e_d_rux * theFormulas[frm].r_u[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_ruy * theFormulas[frm].r_u[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_ruz * theFormulas[frm].r_u[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_rvx * theFormulas[frm].r_v[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_rvy * theFormulas[frm].r_v[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_rvz * theFormulas[frm].r_v[p] * alpha_z; 
				
				
				dedr[3*p+0] += d_e_d_rux * theFormulas[frm].r_u[p]; 
				dedr[3*p+1] += d_e_d_ruy * theFormulas[frm].r_u[p]; 
				dedr[3*p+2] += d_e_d_ruz * theFormulas[frm].r_u[p]; 
				
				dedr[3*p+0] += d_e_d_rvx * theFormulas[frm].r_v[p]; 
				dedr[3*p+1] += d_e_d_rvy * theFormulas[frm].r_v[p]; 
				dedr[3*p+2] += d_e_d_rvz * theFormulas[frm].r_v[p]; 
			}
			
			for( int p = 0; p < np; p++ )
			{
				gr[3*nv+0] += dedr[3*p+0] * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				gr[3*nv+1] += dedr[3*p+1] * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				gr[3*nv+2] += dedr[3*p+2] * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			}
		}	
	} 

	for( int f = 0; f < nf_irr_faces; f++ )
	{

		int frm = f*nf_irr_pts;
		int t = theIrregularFormulas[frm].tri;
		//double g0 = f_g0[f+nf_faces];
		//double g0 = theIrregularFormulas[frm].g0;
		double g0 = theIrregularFormulas[frm].g0;//2*f_g0[nf_faces+f];
		double RuRv0 = theIrregularFormulas[frm].RuRv0;

		for( int p = 0; p < nf_irr_pts; p++ )
		{
			int frm = f*nf_irr_pts+p;
			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
			double nrm[3]={0,0,0}; 

			int *cp = theIrregularFormulas[f*nf_irr_pts+p].cp;
			int np = theIrregularFormulas[f*nf_irr_pts+p].ncoor;
			double dedr[3*np];
			memset( dedr, 0, sizeof(double) * 3 *np );

			for( int p = 0; p < np; p++ )
			{
				R[0] += theIrregularFormulas[frm].r_w[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				R[1] += theIrregularFormulas[frm].r_w[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				R[2] += theIrregularFormulas[frm].r_w[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				Ru[0] += theIrregularFormulas[frm].r_u[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				Ru[1] += theIrregularFormulas[frm].r_u[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				Ru[2] += theIrregularFormulas[frm].r_u[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				Rv[0] += theIrregularFormulas[frm].r_v[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				Rv[1] += theIrregularFormulas[frm].r_v[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				Rv[2] += theIrregularFormulas[frm].r_v[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSuu[0] += theIrregularFormulas[frm].r_uu[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSuu[1] += theIrregularFormulas[frm].r_uu[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSuu[2] += theIrregularFormulas[frm].r_uu[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSuv[0] += theIrregularFormulas[frm].r_uv[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSuv[1] += theIrregularFormulas[frm].r_uv[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSuv[2] += theIrregularFormulas[frm].r_uv[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSvv[0] += theIrregularFormulas[frm].r_vv[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSvv[1] += theIrregularFormulas[frm].r_vv[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSvv[2] += theIrregularFormulas[frm].r_vv[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
			}
		



		cross( Ru, Rv, nrm );
		normalize(nrm);

		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		double g = sqrt(RuRu*RvRv-RuRv*RuRv);

double d_nrmz_d_rux=0,d_nrmz_d_ruy=0,d_nrmz_d_ruz=0,d_nrmz_d_rvx=0,d_nrmz_d_rvy=0,d_nrmz_d_rvz=0,d_nrmy_d_rux=0,d_nrmy_d_ruy=0,d_nrmy_d_ruz=0,d_nrmy_d_rvx=0,d_nrmy_d_rvy=0,d_nrmy_d_rvz=0,d_nrmx_d_rux=0,d_nrmx_d_ruy=0,d_nrmx_d_ruz=0,d_nrmx_d_rvx=0,d_nrmx_d_rvy=0,d_nrmx_d_rvz=0,d_RuRu_d_rux=0,d_RuRu_d_ruy=0,d_RuRu_d_ruz=0,d_nsvv_d_rux=0,d_nsvv_d_ruy=0,d_nsvv_d_ruz=0,d_nsvv_d_rvx=0,d_nsvv_d_rvy=0,d_nsvv_d_rvz=0,d_nsvv_d_svvx=0,d_nsvv_d_svvy=0,d_nsvv_d_svvz=0,d_RuRv_d_rux=0,d_RuRv_d_ruy=0,d_RuRv_d_ruz=0,d_RuRv_d_rvx=0,d_RuRv_d_rvy=0,d_RuRv_d_rvz=0,d_nsuv_d_rux=0,d_nsuv_d_ruy=0,d_nsuv_d_ruz=0,d_nsuv_d_rvx=0,d_nsuv_d_rvy=0,d_nsuv_d_rvz=0,d_nsuv_d_suvx=0,d_nsuv_d_suvy=0,d_nsuv_d_suvz=0,d_RvRv_d_rvx=0,d_RvRv_d_rvy=0,d_RvRv_d_rvz=0,d_nsuu_d_rux=0,d_nsuu_d_ruy=0,d_nsuu_d_ruz=0,d_nsuu_d_rvx=0,d_nsuu_d_rvy=0,d_nsuu_d_rvz=0,d_nsuu_d_suux=0,d_nsuu_d_suuy=0,d_nsuu_d_suuz=0,d_g_d_rux=0,d_g_d_ruy=0,d_g_d_ruz=0,d_g_d_rvx=0,d_g_d_rvy=0,d_g_d_rvz=0,d_e_d_rux=0,d_e_d_ruy=0,d_e_d_ruz=0,d_e_d_rvx=0,d_e_d_rvy=0,d_e_d_rvz=0,d_e_d_suux=0,d_e_d_suuy=0,d_e_d_suuz=0,d_e_d_suvx=0,d_e_d_suvy=0,d_e_d_suvz=0,d_e_d_svvx=0,d_e_d_svvy=0,d_e_d_svvz=0,d_c2_d_rux=0,d_c2_d_ruy=0,d_c2_d_ruz=0,d_c2_d_rvx=0,d_c2_d_rvy=0,d_c2_d_rvz=0,d_c2_d_suux=0,d_c2_d_suuy=0,d_c2_d_suuz=0,d_c2_d_suvx=0,d_c2_d_suvy=0,d_c2_d_suvz=0,d_c2_d_svvx=0,d_c2_d_svvy=0,d_c2_d_svvz=0,d_d_d_rux=0,d_d_d_ruy=0,d_d_d_ruz=0,d_d_d_rvx=0,d_d_d_rvy=0,d_d_d_rvz=0,d_d_d_suvx=0,d_d_d_suvy=0,d_d_d_suvz=0,d_d_d_svvx=0,d_d_d_svvy=0,d_d_d_svvz=0,d_c_d_rux=0,d_c_d_ruy=0,d_c_d_ruz=0,d_c_d_rvx=0,d_c_d_rvy=0,d_c_d_rvz=0,d_c_d_suux=0,d_c_d_suuy=0,d_c_d_suuz=0,d_c_d_suvx=0,d_c_d_suvy=0,d_c_d_suvz=0,d_b_d_rux=0,d_b_d_ruy=0,d_b_d_ruz=0,d_b_d_rvx=0,d_b_d_rvy=0,d_b_d_rvz=0,d_b_d_suvx=0,d_b_d_suvy=0,d_b_d_suvz=0,d_b_d_svvx=0,d_b_d_svvy=0,d_b_d_svvz=0,d_a_d_rux=0,d_a_d_ruy=0,d_a_d_ruz=0,d_a_d_rvx=0,d_a_d_rvy=0,d_a_d_rvz=0,d_a_d_suux=0,d_a_d_suuy=0,d_a_d_suuz=0,d_a_d_suvx=0,d_a_d_suvy=0,d_a_d_suvz=0,d_c1_d_rux=0,d_c1_d_ruy=0,d_c1_d_ruz=0,d_c1_d_rvx=0,d_c1_d_rvy=0,d_c1_d_rvz=0,d_c1_d_suux=0,d_c1_d_suuy=0,d_c1_d_suuz=0,d_c1_d_suvx=0,d_c1_d_suvy=0,d_c1_d_suvz=0,d_c1_d_svvx=0,d_c1_d_svvy=0,d_c1_d_svvz=0,d_g_d_RvRv=0,d_g_d_RuRv=0,d_g_d_RuRu=0,d_nsvv_d_nrmz=0,d_nsvv_d_nrmy=0,d_nsvv_d_nrmx=0,d_nsuv_d_nrmz=0,d_nsuv_d_nrmy=0,d_nsuv_d_nrmx=0,d_nsuu_d_nrmz=0,d_nsuu_d_nrmy=0,d_nsuu_d_nrmx=0,d_d_d_RuRv=0,d_d_d_nsuv=0,d_d_d_RuRu=0,d_d_d_nsvv=0,d_d_d_g=0,d_c_d_RuRv=0,d_c_d_nsuv=0,d_c_d_RuRu=0,d_c_d_nsuu=0,d_c_d_g=0,d_b_d_RuRv=0,d_b_d_nsuv=0,d_b_d_RvRv=0,d_b_d_nsvv=0,d_b_d_g=0,d_a_d_RuRv=0,d_a_d_nsuv=0,d_a_d_RvRv=0,d_a_d_nsuu=0,d_a_d_g=0,d_e_d_c2=0,d_e_d_c1=0,d_c2_d_d=0,d_c2_d_c=0,d_c2_d_b=0,d_c2_d_a=0,d_c1_d_d=0,d_c1_d_c=0,d_c1_d_b=0,d_c1_d_a=0,junk;

			// the basic variables.

			d_g_d_RuRu = RvRv/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
			d_g_d_RuRv = -(RuRv/Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
			d_g_d_RvRv = RuRu/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  

			d_RuRu_d_rux = 2*Ru[0];
			d_RuRu_d_ruy = 2*Ru[1];
			d_RuRu_d_ruz = 2*Ru[2];
			
			d_RuRv_d_rux = Rv[0];
			d_RuRv_d_ruy = Rv[1];
			d_RuRv_d_ruz = Rv[2];
			
			d_RuRv_d_rvx = Ru[0];
			d_RuRv_d_rvy = Ru[1];
			d_RuRv_d_rvz = Ru[2];
			
			d_RvRv_d_rvx = 2*Rv[0];
			d_RvRv_d_rvy = 2*Rv[1];
			d_RvRv_d_rvz = 2*Rv[2];

	d_g_d_rvx += d_g_d_RvRv * d_RvRv_d_rvx;
	d_g_d_rvy += d_g_d_RvRv * d_RvRv_d_rvy;
	d_g_d_rvz += d_g_d_RvRv * d_RvRv_d_rvz;
	d_g_d_rux += d_g_d_RuRv * d_RuRv_d_rux;
	d_g_d_ruy += d_g_d_RuRv * d_RuRv_d_ruy;
	d_g_d_ruz += d_g_d_RuRv * d_RuRv_d_ruz;
	d_g_d_rvx += d_g_d_RuRv * d_RuRv_d_rvx;
	d_g_d_rvy += d_g_d_RuRv * d_RuRv_d_rvy;
	d_g_d_rvz += d_g_d_RuRv * d_RuRv_d_rvz;
	d_g_d_rux += d_g_d_RuRu * d_RuRu_d_rux;
	d_g_d_ruy += d_g_d_RuRu * d_RuRu_d_ruy;
	d_g_d_ruz += d_g_d_RuRu * d_RuRu_d_ruz;

	d_e_d_rux += d_e_d_g[f+nf_faces] * d_g_d_rux * theIrregularFormulas[frm].weight;
	d_e_d_ruy += d_e_d_g[f+nf_faces] * d_g_d_ruy * theIrregularFormulas[frm].weight;
	d_e_d_ruz += d_e_d_g[f+nf_faces] * d_g_d_ruz * theIrregularFormulas[frm].weight;
	d_e_d_rvx += d_e_d_g[f+nf_faces] * d_g_d_rvx * theIrregularFormulas[frm].weight;
	d_e_d_rvy += d_e_d_g[f+nf_faces] * d_g_d_rvy * theIrregularFormulas[frm].weight;
	d_e_d_rvz += d_e_d_g[f+nf_faces] * d_g_d_rvz * theIrregularFormulas[frm].weight;
			
			d_e_d_g0[f+nf_faces] +=  KA * (1 - g*g / (g0*g0) ) *  theIrregularFormulas[frm].weight; 

			for( int p = 0; p < np; p++ )
			{
				gr[3*cp[p]+0] += d_e_d_rux * theIrregularFormulas[frm].r_u[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_ruy * theIrregularFormulas[frm].r_u[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_ruz * theIrregularFormulas[frm].r_u[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_rvx * theIrregularFormulas[frm].r_v[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_rvy * theIrregularFormulas[frm].r_v[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_rvz * theIrregularFormulas[frm].r_v[p] * alpha_z; 
				
				dedr[3*p+0] += d_e_d_rux * theIrregularFormulas[frm].r_u[p]; 
				dedr[3*p+1] += d_e_d_ruy * theIrregularFormulas[frm].r_u[p]; 
				dedr[3*p+2] += d_e_d_ruz * theIrregularFormulas[frm].r_u[p]; 
				
				dedr[3*p+0] += d_e_d_rvx * theIrregularFormulas[frm].r_v[p]; 
				dedr[3*p+1] += d_e_d_rvy * theIrregularFormulas[frm].r_v[p]; 
				dedr[3*p+2] += d_e_d_rvz * theIrregularFormulas[frm].r_v[p]; 

			}
			for( int p = 0; p < np; p++ )
			{
				gr[3*nv+0] += dedr[3*p+0] * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				gr[3*nv+1] += dedr[3*p+1] * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				gr[3*nv+2] += dedr[3*p+2] * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
			}

		}	

	} 

	r_val /= area;
	
	if( vary_g )
	{
		// derivative of the energy with respect to g0
		double rsum[f_nregions];
		memset( rsum, 0, sizeof(double) * f_nregions );
			
		for( int f2 = 0; f2 < nt; f2++ )
			rsum[f_g_region[f2]] += d_e_d_g0[f2] * f_g0[f2];

		for( int f1 = 0; f1 < nt; f1++ )
		{
				// dE_dg0_f1
			gr[3*nv+3+f1] += rsum[f_g_region[f1]] * (-1.0/f_g_sum0[f_g_region[f1]]); 
			gr[3*nv+3+f1] += d_e_d_g0[f1]; 
		}	
	}


	// loop over frames, computing dg/dr everything.

	free( strains );
	free( d_e_d_g );
	free( d_e_d_g0 );
}
#endif

void surface::pointGradient( int face, double u, double v, double *r, double *gr, double *g_puv, double *de_dr, double *de_dnrm, double frac_mult )
{
	// updates the gradient wrt to vertices as well as the point (g_puv).
	// input: de_dr, the derivative wrt this point.
	//	  de_dnrm, the derivative wrt this norm.

	// need, derivative of the point d_rx wrt vertex positions	
	
	double alpha_x = r[3*nv];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];
	
	if( !theFormulas )
		generatePlan();

	if( face < nf_faces )
	{
		int f = face;
		int frm = f*nf_g_q_p;
		int t = theFormulas[frm].tri;

		int *cp = theFormulas[frm].cp;

		double w = 1-u-v;

		double u2 = u*u;
		double u3 = u*u*u;
		double u4 = u*u*u*u;
		
		double v2 = v*v;
		double v3 = v*v*v;
		double v4 = v*v*v*v;
		
		double w2 = w*w;
		double w3 = w*w*w;
		double w4 = w*w*w*w;

		double n1 = (1.0/12.0)*(u4+2*u3*v); 
		double n2 = (1.0/12.0)*(u4+2*u3*w); 
		double n3 = (1.0/12.0)*(u4+2*u3*w+6*u3*v+6*u2*v*w+12*u2*v2+6*u*v2*w+6*u*v3+2*v3*w+v4); 
		double n4 = (1.0/12.0)*(6*u4+24*u3*w+24*u2*w2+8*u*w3+w4+24*u3*v+60*u2*v*w+36*u*v*w2+6*v*w3+24*u2*v2+36*u*v2*w+12*v2*w2+8*u*v3+6*v3*w+v4); 
		double n5 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+2*u3*v+6*u2*v*w+6*u*v*w2+2*v*w3); 
		double n6 = (1.0/12.0)*(2*u*v3+v4); 
		double n7 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+8*u3*v+36*u2*v*w+36*u*v*w2+8*v*w3+24*u2*v2+60*u*v2*w+24*v2*w2+24*u*v3+24*v3*w+6*v4); 
		double n8 = (1.0/12.0)*(u4+8*u3*w+24*u2*w2+24*u*w3+6*w4+6*u3*v+36*u2*v*w+60*u*v*w2+24*v*w3+12*u2*v2+36*u*v2*w+24*v2*w2+6*u*v3+8*v3*w+v4); 
		double n9 = (1.0/12.0)*(2*u*w3+w4); 
		double n10 = (1.0/12.0)*(2*v3*w+v4); 
		double n11 = (1.0/12.0)*(2*u*w3+w4+6*u*v*w2+6*v*w3+6*u*v2*w+12*v2*w2+2*u*v3+6*v3*w+v4); 
		double n12 = (1.0/12.0)*(w4+2*v*w3);

		double du_1 = Power(u,3)/3. + (Power(u,2)*v)/2.;
		double du_2 = Power(u,2)/2. - Power(u,3)/3. - (Power(u,2)*v)/2.; 	
		double du_3 = Power(u,2)/2. - Power(u,3)/3. + u*v - (Power(u,2)*v)/2. + Power(v,2)/2. - Power(v,3)/6.;
		double du_4 = 0.3333333333333333 + u - Power(u,2) - Power(u,3)/3. + v/2. - u*v - (Power(u,2)*v)/2. - Power(v,2) + Power(v,3)/3.;	
		double du_5 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. - v/2. + Power(u,2)*v + Power(v,2)/2. - Power(v,3)/6.;
		double du_6 = Power(v,3)/6.;	
		double du_7 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. + v/2. - 2*u*v + Power(u,2)*v - Power(v,2)/2. - Power(v,3)/6.;	
		double du_8 = -2*u + 2*Power(u,2) - Power(u,3)/3. - v + 2*u*v - (Power(u,2)*v)/2. + Power(v,2) - Power(v,3)/6.;	
		double du_9 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. + v/2. - (Power(u,2)*v)/2. - Power(v,2)/2. + Power(v,3)/6.;	
		double du_10 = -Power(v,3)/6.;
		double du_11 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. - v/2. + u*v - (Power(u,2)*v)/2. + Power(v,3)/3.;	
		double du_12 = 	-0.3333333333333333 + u - Power(u,2) + Power(u,3)/3. + v/2. - u*v + (Power(u,2)*v)/2. - Power(v,3)/6.;

		double dv_1 = Power(u,3)/6.;
		double dv_2 = -Power(u,3)/6.;
		double dv_3 = Power(u,2)/2. - Power(u,3)/6. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_4 = 0.16666666666666666 + u/2. - Power(u,2)/2. - Power(u,3)/6. - 2*u*v - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
		double dv_5 = -0.16666666666666666 - u/2. + Power(u,3)/3. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_6 = (u*Power(v,2))/2. + Power(v,3)/3.;
		double dv_7 = 0.3333333333333333 + u/2. - Power(u,2) + Power(u,3)/3. + v - u*v - Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_8 = -u + Power(u,2) - Power(u,3)/6. - 2*v + 2*u*v + 2*Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_9 = -0.3333333333333333 + u/2. - Power(u,3)/6. + v - u*v - Power(v,2) + (u*Power(v,2))/2. + Power(v,3)/3.;
		double dv_10 = Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_11 = 0.16666666666666666 - u/2. + Power(u,2)/2. - Power(u,3)/6. - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
		double dv_12 = -0.16666666666666666 + u/2. - Power(u,2)/2. + Power(u,3)/6. + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;

		double d_uu_1 = Power(u,2) + u*v;
		double d_uu_2 = u - Power(u,2) - u*v;
		double d_uu_3 = u - Power(u,2) + v - u*v;
		double d_uu_4 = 1 - 2*u - Power(u,2) - v - u*v;
		double d_uu_5 = -2*u + 2*Power(u,2) + 2*u*v;
		double d_uu_6 = 0;
		double d_uu_7 = -2*u + 2*Power(u,2) - 2*v + 2*u*v;
		double d_uu_8 = -2 + 4*u - Power(u,2) + 2*v - u*v;
		double d_uu_9 = u - Power(u,2) - u*v;
		double d_uu_10 = 0;
		double d_uu_11 = u - Power(u,2) + v - u*v;
		double d_uu_12 = 1 - 2*u + Power(u,2) - v + u*v;
		
		double d_uv_1 = Power(u,2)/2.;
		double d_uv_2 = -Power(u,2)/2.;
		double d_uv_3 = u - Power(u,2)/2. + v - Power(v,2)/2.;
		double d_uv_4 = 0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
		double d_uv_5 = -0.5 + Power(u,2) + v - Power(v,2)/2.;
		double d_uv_6 = Power(v,2)/2.;
		double d_uv_7 = 0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
		double d_uv_8 = -1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
		double d_uv_9 = 0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
		double d_uv_10 = -Power(v,2)/2.;
		double d_uv_11 = -0.5 + u - Power(u,2)/2. + Power(v,2);
		double d_uv_12 = 0.5 - u + Power(u,2)/2. - Power(v,2)/2.;
		
		double d_vv_1 = 0;
		double d_vv_2 = 0;
		double d_vv_3 = u + v - u*v - Power(v,2);
		double d_vv_4 = -2*u - 2*v + 2*u*v + 2*Power(v,2);
		double d_vv_5 = u + v - u*v - Power(v,2);
		double d_vv_6 = u*v + Power(v,2);
		double d_vv_7 = 1 - u - 2*v - u*v - Power(v,2);
		double d_vv_8 = -2 + 2*u + 4*v - u*v - Power(v,2);
		double d_vv_9 = 1 - u - 2*v + u*v + Power(v,2);
		double d_vv_10 = v - u*v - Power(v,2);
		double d_vv_11 = -2*v + 2*u*v + 2*Power(v,2);
		double d_vv_12 = v - u*v - Power(v,2);
	
		double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
		double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
		double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
		
		double ceff_map_duu[12] = { d_uu_8, d_uu_7, d_uu_4, d_uu_5, d_uu_9, d_uu_12, d_uu_11, d_uu_10, d_uu_6, d_uu_3, d_uu_1, d_uu_2 };
		double ceff_map_duv[12] = { d_uv_8, d_uv_7, d_uv_4, d_uv_5, d_uv_9, d_uv_12, d_uv_11, d_uv_10, d_uv_6, d_uv_3, d_uv_1, d_uv_2 };
		double ceff_map_dvv[12] = { d_vv_8, d_vv_7, d_vv_4, d_vv_5, d_vv_9, d_vv_12, d_vv_11, d_vv_10, d_vv_6, d_vv_3, d_vv_1, d_vv_2 };
				
		double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
		double nrm[3]={0,0,0}; 

		int np = theFormulas[frm].ncoor;
		double dedr[3*np];
		memset( dedr, 0, sizeof(double) * 3 *np );

		for( int p = 0; p < np; p++ )
		{

			R[0] += ceff_map[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			R[1] += ceff_map[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			R[2] += ceff_map[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			Ru[0] += ceff_map_du[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			Ru[1] += ceff_map_du[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			Ru[2] += ceff_map_du[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			Rv[0] += ceff_map_dv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			Rv[1] += ceff_map_dv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			Rv[2] += ceff_map_dv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
			tSuu[0] += ceff_map_duu[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSuu[1] += ceff_map_duu[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSuu[2] += ceff_map_duu[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSuv[0] += ceff_map_duv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSuv[1] += ceff_map_duv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSuv[2] += ceff_map_duv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSvv[0] += ceff_map_dvv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSvv[1] += ceff_map_dvv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSvv[2] += ceff_map_dvv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
		}
		cross( Ru, Rv, nrm );
		normalize(nrm);

		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];
		
		double g = sqrt(RuRu*RvRv-RuRv*RuRv);


double d_nrmz_d_rux=0,d_nrmz_d_ruy=0,d_nrmz_d_ruz=0,d_nrmz_d_rvx=0,d_nrmz_d_rvy=0,d_nrmz_d_rvz=0,d_nrmy_d_rux=0,d_nrmy_d_ruy=0,d_nrmy_d_ruz=0,d_nrmy_d_rvx=0,d_nrmy_d_rvy=0,d_nrmy_d_rvz=0,d_nrmx_d_rux=0,d_nrmx_d_ruy=0,d_nrmx_d_ruz=0,d_nrmx_d_rvx=0,d_nrmx_d_rvy=0,d_nrmx_d_rvz=0,d_RuRu_d_rux=0,d_RuRu_d_ruy=0,d_RuRu_d_ruz=0,d_nsvv_d_rux=0,d_nsvv_d_ruy=0,d_nsvv_d_ruz=0,d_nsvv_d_rvx=0,d_nsvv_d_rvy=0,d_nsvv_d_rvz=0,d_nsvv_d_svvx=0,d_nsvv_d_svvy=0,d_nsvv_d_svvz=0,d_RuRv_d_rux=0,d_RuRv_d_ruy=0,d_RuRv_d_ruz=0,d_RuRv_d_rvx=0,d_RuRv_d_rvy=0,d_RuRv_d_rvz=0,d_nsuv_d_rux=0,d_nsuv_d_ruy=0,d_nsuv_d_ruz=0,d_nsuv_d_rvx=0,d_nsuv_d_rvy=0,d_nsuv_d_rvz=0,d_nsuv_d_suvx=0,d_nsuv_d_suvy=0,d_nsuv_d_suvz=0,d_RvRv_d_rvx=0,d_RvRv_d_rvy=0,d_RvRv_d_rvz=0,d_nsuu_d_rux=0,d_nsuu_d_ruy=0,d_nsuu_d_ruz=0,d_nsuu_d_rvx=0,d_nsuu_d_rvy=0,d_nsuu_d_rvz=0,d_nsuu_d_suux=0,d_nsuu_d_suuy=0,d_nsuu_d_suuz=0,d_g_d_rux=0,d_g_d_ruy=0,d_g_d_ruz=0,d_g_d_rvx=0,d_g_d_rvy=0,d_g_d_rvz=0,d_e_d_rux=0,d_e_d_ruy=0,d_e_d_ruz=0,d_e_d_rvx=0,d_e_d_rvy=0,d_e_d_rvz=0,d_e_d_suux=0,d_e_d_suuy=0,d_e_d_suuz=0,d_e_d_suvx=0,d_e_d_suvy=0,d_e_d_suvz=0,d_e_d_svvx=0,d_e_d_svvy=0,d_e_d_svvz=0,d_c2_d_rux=0,d_c2_d_ruy=0,d_c2_d_ruz=0,d_c2_d_rvx=0,d_c2_d_rvy=0,d_c2_d_rvz=0,d_c2_d_suux=0,d_c2_d_suuy=0,d_c2_d_suuz=0,d_c2_d_suvx=0,d_c2_d_suvy=0,d_c2_d_suvz=0,d_c2_d_svvx=0,d_c2_d_svvy=0,d_c2_d_svvz=0,d_d_d_rux=0,d_d_d_ruy=0,d_d_d_ruz=0,d_d_d_rvx=0,d_d_d_rvy=0,d_d_d_rvz=0,d_d_d_suvx=0,d_d_d_suvy=0,d_d_d_suvz=0,d_d_d_svvx=0,d_d_d_svvy=0,d_d_d_svvz=0,d_c_d_rux=0,d_c_d_ruy=0,d_c_d_ruz=0,d_c_d_rvx=0,d_c_d_rvy=0,d_c_d_rvz=0,d_c_d_suux=0,d_c_d_suuy=0,d_c_d_suuz=0,d_c_d_suvx=0,d_c_d_suvy=0,d_c_d_suvz=0,d_b_d_rux=0,d_b_d_ruy=0,d_b_d_ruz=0,d_b_d_rvx=0,d_b_d_rvy=0,d_b_d_rvz=0,d_b_d_suvx=0,d_b_d_suvy=0,d_b_d_suvz=0,d_b_d_svvx=0,d_b_d_svvy=0,d_b_d_svvz=0,d_a_d_rux=0,d_a_d_ruy=0,d_a_d_ruz=0,d_a_d_rvx=0,d_a_d_rvy=0,d_a_d_rvz=0,d_a_d_suux=0,d_a_d_suuy=0,d_a_d_suuz=0,d_a_d_suvx=0,d_a_d_suvy=0,d_a_d_suvz=0,d_c1_d_rux=0,d_c1_d_ruy=0,d_c1_d_ruz=0,d_c1_d_rvx=0,d_c1_d_rvy=0,d_c1_d_rvz=0,d_c1_d_suux=0,d_c1_d_suuy=0,d_c1_d_suuz=0,d_c1_d_suvx=0,d_c1_d_suvy=0,d_c1_d_suvz=0,d_c1_d_svvx=0,d_c1_d_svvy=0,d_c1_d_svvz=0,d_g_d_RvRv=0,d_g_d_RuRv=0,d_g_d_RuRu=0,d_nsvv_d_nrmz=0,d_nsvv_d_nrmy=0,d_nsvv_d_nrmx=0,d_nsuv_d_nrmz=0,d_nsuv_d_nrmy=0,d_nsuv_d_nrmx=0,d_nsuu_d_nrmz=0,d_nsuu_d_nrmy=0,d_nsuu_d_nrmx=0,d_d_d_RuRv=0,d_d_d_nsuv=0,d_d_d_RuRu=0,d_d_d_nsvv=0,d_d_d_g=0,d_c_d_RuRv=0,d_c_d_nsuv=0,d_c_d_RuRu=0,d_c_d_nsuu=0,d_c_d_g=0,d_b_d_RuRv=0,d_b_d_nsuv=0,d_b_d_RvRv=0,d_b_d_nsvv=0,d_b_d_g=0,d_a_d_RuRv=0,d_a_d_nsuv=0,d_a_d_RvRv=0,d_a_d_nsuu=0,d_a_d_g=0,d_e_d_c2=0,d_e_d_c1=0,d_e_d_g=0,d_c2_d_d=0,d_c2_d_c=0,d_c2_d_b=0,d_c2_d_a=0,d_c1_d_d=0,d_c1_d_c=0,d_c1_d_b=0,d_c1_d_a=0,junk;
		// the basic variables.

		d_g_d_RuRu = RvRv/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
		d_g_d_RuRv = -(RuRv/Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
		d_g_d_RvRv = RuRu/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
	
		d_RuRu_d_rux = 2*Ru[0];
		d_RuRu_d_ruy = 2*Ru[1];
		d_RuRu_d_ruz = 2*Ru[2];
		
		d_RuRv_d_rux = Rv[0];
		d_RuRv_d_ruy = Rv[1];
		d_RuRv_d_ruz = Rv[2];
		
		d_RuRv_d_rvx = Ru[0];
		d_RuRv_d_rvy = Ru[1];
		d_RuRv_d_rvz = Ru[2];
		
		d_RvRv_d_rvx = 2*Rv[0];
		d_RvRv_d_rvy = 2*Rv[1];
		d_RvRv_d_rvz = 2*Rv[2];
		
		d_g_d_rvx += d_g_d_RvRv * d_RvRv_d_rvx;
		d_g_d_rvy += d_g_d_RvRv * d_RvRv_d_rvy;
		d_g_d_rvz += d_g_d_RvRv * d_RvRv_d_rvz;
		d_g_d_rux += d_g_d_RuRv * d_RuRv_d_rux;
		d_g_d_ruy += d_g_d_RuRv * d_RuRv_d_ruy;
		d_g_d_ruz += d_g_d_RuRv * d_RuRv_d_ruz;
		d_g_d_rvx += d_g_d_RuRv * d_RuRv_d_rvx;
		d_g_d_rvy += d_g_d_RuRv * d_RuRv_d_rvy;
		d_g_d_rvz += d_g_d_RuRv * d_RuRv_d_rvz;
		d_g_d_rux += d_g_d_RuRu * d_RuRu_d_rux;
		d_g_d_ruy += d_g_d_RuRu * d_RuRu_d_ruy;
		d_g_d_ruz += d_g_d_RuRu * d_RuRu_d_ruz;

		// that's it.

		double fac = Power(Power(Ru[1]*Rv[0] - Ru[0]*Rv[1],2) + Power(Ru[2]*Rv[0] - Ru[0]*Rv[2],2) + Power(Ru[2]*Rv[1] - Ru[1]*Rv[2],2),1.5);

		d_nrmx_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac; 
		d_nrmx_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2)))/fac;
		d_nrmx_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac;
		
		d_nrmy_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2))))/fac;
		d_nrmy_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac; 
		d_nrmy_d_ruz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2])) + Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac;

		d_nrmz_d_rux = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2])/fac;
		d_nrmz_d_ruy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;
		d_nrmz_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2))) + (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;

		d_nrmx_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac; 
		d_nrmx_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2])))/fac;
		d_nrmx_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac;
		
		d_nrmy_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
		d_nrmy_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2])))/fac;
		d_nrmy_d_rvz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[0]*Ru[1]*Rv[0] - Power(Ru[0],2)*Rv[1] + Ru[2]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2])))/fac;

		d_nrmz_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;
		d_nrmz_d_rvy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[0]*Ru[2]*Rv[0] - Power(Ru[0],2)*Rv[2] + Ru[1]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
		d_nrmz_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;
		
		d_e_d_rux = (de_dnrm[0] * d_nrmx_d_rux + de_dnrm[1] * d_nrmy_d_rux + de_dnrm[2] * d_nrmz_d_rux)*frac_mult;
		d_e_d_ruy = (de_dnrm[0] * d_nrmx_d_ruy + de_dnrm[1] * d_nrmy_d_ruy + de_dnrm[2] * d_nrmz_d_ruy)*frac_mult;
		d_e_d_ruz = (de_dnrm[0] * d_nrmx_d_ruz + de_dnrm[1] * d_nrmy_d_ruz + de_dnrm[2] * d_nrmz_d_ruz)*frac_mult;
		
		d_e_d_rvx = (de_dnrm[0] * d_nrmx_d_rvx + de_dnrm[1] * d_nrmy_d_rvx + de_dnrm[2] * d_nrmz_d_rvx)*frac_mult;
		d_e_d_rvy = (de_dnrm[0] * d_nrmx_d_rvy + de_dnrm[1] * d_nrmy_d_rvy + de_dnrm[2] * d_nrmz_d_rvy)*frac_mult;
		d_e_d_rvz = (de_dnrm[0] * d_nrmx_d_rvz + de_dnrm[1] * d_nrmy_d_rvz + de_dnrm[2] * d_nrmz_d_rvz)*frac_mult;

		double d_e_d_rx = de_dr[0] * frac_mult;
		double d_e_d_ry = de_dr[1] * frac_mult;
		double d_e_d_rz = de_dr[2] * frac_mult;

		// dn/du dn/dv
		double Nu[3]={0,0,0},Nv[3]={0,0,0};

		Nu[0] += d_nrmx_d_rux * tSuu[0] + d_nrmx_d_rvx * tSuv[0];
		Nu[1] += d_nrmy_d_rux * tSuu[0] + d_nrmy_d_rvx * tSuv[0];
		Nu[2] += d_nrmz_d_rux * tSuu[0] + d_nrmz_d_rvx * tSuv[0];
		                                           
		Nu[0] += d_nrmx_d_ruy * tSuu[1] + d_nrmx_d_rvy * tSuv[1];
		Nu[1] += d_nrmy_d_ruy * tSuu[1] + d_nrmy_d_rvy * tSuv[1];
		Nu[2] += d_nrmz_d_ruy * tSuu[1] + d_nrmz_d_rvy * tSuv[1];
                                                           
		Nu[0] += d_nrmx_d_ruz * tSuu[2] + d_nrmx_d_rvz * tSuv[2];
		Nu[1] += d_nrmy_d_ruz * tSuu[2] + d_nrmy_d_rvz * tSuv[2];
		Nu[2] += d_nrmz_d_ruz * tSuu[2] + d_nrmz_d_rvz * tSuv[2];
                                                           
		Nv[0] += d_nrmx_d_rux * tSuv[0] + d_nrmx_d_rvx * tSvv[0];
		Nv[1] += d_nrmy_d_rux * tSuv[0] + d_nrmy_d_rvx * tSvv[0];
		Nv[2] += d_nrmz_d_rux * tSuv[0] + d_nrmz_d_rvx * tSvv[0];
		                                           
		Nv[0] += d_nrmx_d_ruy * tSuv[1] + d_nrmx_d_rvy * tSvv[1];
		Nv[1] += d_nrmy_d_ruy * tSuv[1] + d_nrmy_d_rvy * tSvv[1];
		Nv[2] += d_nrmz_d_ruy * tSuv[1] + d_nrmz_d_rvy * tSvv[1];
                                                           
		Nv[0] += d_nrmx_d_ruz * tSuv[2] + d_nrmx_d_rvz * tSvv[2];
		Nv[1] += d_nrmy_d_ruz * tSuv[2] + d_nrmy_d_rvz * tSvv[2];
		Nv[2] += d_nrmz_d_ruz * tSuv[2] + d_nrmz_d_rvz * tSvv[2];

		nrm_grad[0] = Nu[0];
		nrm_grad[1] = Nu[1];
		nrm_grad[2] = Nu[2];
		
		nrm_grad[3] = Nv[0];
		nrm_grad[4] = Nv[1];
		nrm_grad[5] = Nv[2];

		// dE/du dE/dv


		g_puv[0] += de_dr[0] * Ru[0] + de_dr[1] * Ru[1] + de_dr[2] * Ru[2];
		g_puv[1] += de_dr[0] * Rv[0] + de_dr[1] * Rv[1] + de_dr[2] * Rv[2];
		
		g_puv[0] += de_dnrm[0] * Nu[0] + de_dnrm[1] * Nu[1] + de_dnrm[2] * Nu[2];
		g_puv[1] += de_dnrm[0] * Nv[0] + de_dnrm[1] * Nv[1] + de_dnrm[2] * Nv[2];

	

		// normal is constant wrt to u and v at first order.
	
		for( int p = 0; p < np; p++ )
		{
			gr[3*cp[p]+0] += d_e_d_rx * ceff_map[p] * alpha_x;
			gr[3*cp[p]+1] += d_e_d_ry * ceff_map[p] * alpha_y;
			gr[3*cp[p]+2] += d_e_d_rz * ceff_map[p] * alpha_z;

			gr[3*cp[p]+0] += d_e_d_rux * ceff_map_du[p] * alpha_x; 
			gr[3*cp[p]+1] += d_e_d_ruy * ceff_map_du[p] * alpha_y; 
			gr[3*cp[p]+2] += d_e_d_ruz * ceff_map_du[p] * alpha_z; 
			
			gr[3*cp[p]+0] += d_e_d_rvx * ceff_map_dv[p] * alpha_x; 
			gr[3*cp[p]+1] += d_e_d_rvy * ceff_map_dv[p] * alpha_y; 
			gr[3*cp[p]+2] += d_e_d_rvz * ceff_map_dv[p] * alpha_z; 
			
			dedr[3*p+0] += d_e_d_rx * ceff_map[p];
			dedr[3*p+1] += d_e_d_ry * ceff_map[p];
			dedr[3*p+2] += d_e_d_rz * ceff_map[p];
			
			dedr[3*p+0] += d_e_d_rux * ceff_map_du[p];
			dedr[3*p+1] += d_e_d_ruy * ceff_map_du[p];
			dedr[3*p+2] += d_e_d_ruz * ceff_map_du[p];
			                                         
			dedr[3*p+0] += d_e_d_rvx * ceff_map_dv[p];
			dedr[3*p+1] += d_e_d_rvy * ceff_map_dv[p];
			dedr[3*p+2] += d_e_d_rvz * ceff_map_dv[p];
		}
		
		if ( on_surface )
		  {


		    for( int p = 0; p < np; p++ )
		      {
			
			gr[3*cp[p]+0] -= (1/g) * kT * d_g_d_rux * ceff_map_du[p] * alpha_x;
			gr[3*cp[p]+1] -= (1/g) * kT * d_g_d_ruy * ceff_map_du[p] * alpha_y;
			gr[3*cp[p]+2] -= (1/g) * kT * d_g_d_ruz * ceff_map_du[p] * alpha_z;

			gr[3*cp[p]+0] -= (1/g) * kT * d_g_d_rvx * ceff_map_dv[p] * alpha_x;
			gr[3*cp[p]+1] -= (1/g) * kT * d_g_d_rvy * ceff_map_dv[p] * alpha_y;
			gr[3*cp[p]+2] -= (1/g) * kT * d_g_d_rvz * ceff_map_dv[p] * alpha_z;

			dedr[3*p+0] -= (1/g) * kT * d_g_d_rux * ceff_map_du[p];
			dedr[3*p+1] -= (1/g) * kT * d_g_d_ruy * ceff_map_du[p];
			dedr[3*p+2] -= (1/g) * kT * d_g_d_ruz * ceff_map_du[p];

			dedr[3*p+0] -= (1/g) * kT * d_g_d_rvx * ceff_map_dv[p];
			dedr[3*p+1] -= (1/g) * kT * d_g_d_rvy * ceff_map_dv[p];
			dedr[3*p+2] -= (1/g) * kT * d_g_d_rvz * ceff_map_dv[p];

		      }
		  }

		for( int p = 0; p < np; p++ )
		{
			gr[3*nv+0] += dedr[3*p+0] * (r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			gr[3*nv+1] += dedr[3*p+1] * (r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			gr[3*nv+2] += dedr[3*p+2] * (r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
		}

	}
	else
	{
		int f = face - nf_faces;
		int frm = f*nf_irr_pts;
		int t = theIrregularFormulas[frm].tri;
		int *cp = theIrregularFormulas[frm].cp;
		int i = theIrregularFormulas[frm].vertex;
		int e = theIrregularFormulas[frm].edge;
		int val = theVertices[i].valence;
		int np = theIrregularFormulas[frm].ncoor;

		irr_kernel *theKernel = kernels[val]; 
		int ncoords_base = val + 6;

		double fu = u;
		double fv = v;					

		int domain = theKernel->domain(fu,fv);

		if( domain-1 >= theKernel->ndomains )
		{
			printf("LOGICAL ERROR evaluating R,nrm on irregular vertex face.\n");
			exit(1);
		}	
		
		if( domain < 0 &&  fu < 1e-10 && fv < 1e-10  )
		{
		}

		double u_u = fu, u_v=0, v_u=0, v_v = fv;
		theKernel->get_map_transform( &u_u, &u_v, &v_u, &v_v );
		double *theMap = theKernel->get_map( &fu, &fv );
		double domain_scale = pow(2,domain);

		double u = fu;
		double v = fv;
		double w = 1 - u - v;


		if( u +v > 1.0 + THRESH || u < -THRESH || v < -THRESH )
		{
			printf("u: %.12le v: %.12le\n", u, v );
			printf("ERROR: outside domain.\n");
			exit(1);
		}

		double u2 = u*u;
		double u3 = u*u*u;
		double u4 = u*u*u*u;
		
		double v2 = v*v;
		double v3 = v*v*v;
		double v4 = v*v*v*v;
		
		double w2 = w*w;
		double w3 = w*w*w;
		double w4 = w*w*w*w;
			
		
		// 8 : 0
		// 7 : 1
		// 4 : 2
		// 5 : 3
		// 9 : 4
		// 12 : 5
		// 11 : 6
		// 10 : 7
		// 6 : 8
		// 3 : 9
		// 1 : 10
		// 2 : 11

	

		double n1 = (1.0/12.0)*(u4+2*u3*v); 
		double n2 = (1.0/12.0)*(u4+2*u3*w); 
		double n3 = (1.0/12.0)*(u4+2*u3*w+6*u3*v+6*u2*v*w+12*u2*v2+6*u*v2*w+6*u*v3+2*v3*w+v4); 
		double n4 = (1.0/12.0)*(6*u4+24*u3*w+24*u2*w2+8*u*w3+w4+24*u3*v+60*u2*v*w+36*u*v*w2+6*v*w3+24*u2*v2+36*u*v2*w+12*v2*w2+8*u*v3+6*v3*w+v4); 
		double n5 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+2*u3*v+6*u2*v*w+6*u*v*w2+2*v*w3); 
		double n6 = (1.0/12.0)*(2*u*v3+v4); 
		double n7 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+8*u3*v+36*u2*v*w+36*u*v*w2+8*v*w3+24*u2*v2+60*u*v2*w+24*v2*w2+24*u*v3+24*v3*w+6*v4); 
		double n8 = (1.0/12.0)*(u4+8*u3*w+24*u2*w2+24*u*w3+6*w4+6*u3*v+36*u2*v*w+60*u*v*w2+24*v*w3+12*u2*v2+36*u*v2*w+24*v2*w2+6*u*v3+8*v3*w+v4); 
		double n9 = (1.0/12.0)*(2*u*w3+w4); 
		double n10 = (1.0/12.0)*(2*v3*w+v4); 
		double n11 = (1.0/12.0)*(2*u*w3+w4+6*u*v*w2+6*v*w3+6*u*v2*w+12*v2*w2+2*u*v3+6*v3*w+v4); 
		double n12 = (1.0/12.0)*(w4+2*v*w3);

		double du_1 = Power(u,3)/3. + (Power(u,2)*v)/2.;
		double du_2 = Power(u,2)/2. - Power(u,3)/3. - (Power(u,2)*v)/2.; 	
		double du_3 =Power(u,2)/2. - Power(u,3)/3. + u*v - (Power(u,2)*v)/2. + Power(v,2)/2. - Power(v,3)/6.;
		double du_4 = 0.3333333333333333 + u - Power(u,2) - Power(u,3)/3. + v/2. - u*v - (Power(u,2)*v)/2. - Power(v,2) + Power(v,3)/3.;	
		double du_5 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. - v/2. + Power(u,2)*v + Power(v,2)/2. - Power(v,3)/6.;
		double du_6 = Power(v,3)/6.;	
		double du_7 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. + v/2. - 2*u*v + Power(u,2)*v - Power(v,2)/2. - Power(v,3)/6.;	
		double du_8 = -2*u + 2*Power(u,2) - Power(u,3)/3. - v + 2*u*v - (Power(u,2)*v)/2. + Power(v,2) - Power(v,3)/6.;	
		double du_9 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. + v/2. - (Power(u,2)*v)/2. - Power(v,2)/2. + Power(v,3)/6.;	
		double du_10 = -Power(v,3)/6.;
		double du_11 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. - v/2. + u*v - (Power(u,2)*v)/2. + Power(v,3)/3.;	
		double du_12 = 	-0.3333333333333333 + u - Power(u,2) + Power(u,3)/3. + v/2. - u*v + (Power(u,2)*v)/2. - Power(v,3)/6.;

		double dv_1 = Power(u,3)/6.;
		double dv_2 = -Power(u,3)/6.;
		double dv_3 = Power(u,2)/2. - Power(u,3)/6. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_4 = 0.16666666666666666 + u/2. - Power(u,2)/2. - Power(u,3)/6. - 2*u*v - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
		double dv_5 = -0.16666666666666666 - u/2. + Power(u,3)/3. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_6 = (u*Power(v,2))/2. + Power(v,3)/3.;
		double dv_7 = 0.3333333333333333 + u/2. - Power(u,2) + Power(u,3)/3. + v - u*v - Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_8 = -u + Power(u,2) - Power(u,3)/6. - 2*v + 2*u*v + 2*Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_9 = -0.3333333333333333 + u/2. - Power(u,3)/6. + v - u*v - Power(v,2) + (u*Power(v,2))/2. + Power(v,3)/3.;
		double dv_10 = Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_11 = 0.16666666666666666 - u/2. + Power(u,2)/2. - Power(u,3)/6. - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
		double dv_12 = -0.16666666666666666 + u/2. - Power(u,2)/2. + Power(u,3)/6. + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
	
		double d_uu_1 = Power(u,2) + u*v;
		double d_uu_2 = u - Power(u,2) - u*v;
		double d_uu_3 = u - Power(u,2) + v - u*v;
		double d_uu_4 = 1 - 2*u - Power(u,2) - v - u*v;
		double d_uu_5 = -2*u + 2*Power(u,2) + 2*u*v;
		double d_uu_6 = 0;
		double d_uu_7 = -2*u + 2*Power(u,2) - 2*v + 2*u*v;
		double d_uu_8 = -2 + 4*u - Power(u,2) + 2*v - u*v;
		double d_uu_9 = u - Power(u,2) - u*v;
		double d_uu_10 = 0;
		double d_uu_11 = u - Power(u,2) + v - u*v;
		double d_uu_12 = 1 - 2*u + Power(u,2) - v + u*v;
		
		double d_uv_1 = Power(u,2)/2.;
		double d_uv_2 = -Power(u,2)/2.;
		double d_uv_3 = u - Power(u,2)/2. + v - Power(v,2)/2.;
		double d_uv_4 = 0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
		double d_uv_5 = -0.5 + Power(u,2) + v - Power(v,2)/2.;
		double d_uv_6 = Power(v,2)/2.;
		double d_uv_7 = 0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
		double d_uv_8 = -1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
		double d_uv_9 = 0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
		double d_uv_10 = -Power(v,2)/2.;
		double d_uv_11 = -0.5 + u - Power(u,2)/2. + Power(v,2);
		double d_uv_12 = 0.5 - u + Power(u,2)/2. - Power(v,2)/2.;
		
		double d_vv_1 = 0;
		double d_vv_2 = 0;
		double d_vv_3 = u + v - u*v - Power(v,2);
		double d_vv_4 = -2*u - 2*v + 2*u*v + 2*Power(v,2);
		double d_vv_5 = u + v - u*v - Power(v,2);
		double d_vv_6 = u*v + Power(v,2);
		double d_vv_7 = 1 - u - 2*v - u*v - Power(v,2);
		double d_vv_8 = -2 + 2*u + 4*v - u*v - Power(v,2);
		double d_vv_9 = 1 - u - 2*v + u*v + Power(v,2);
		double d_vv_10 = v - u*v - Power(v,2);
		double d_vv_11 = -2*v + 2*u*v + 2*Power(v,2);
		double d_vv_12 = v - u*v - Power(v,2);
		
		double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
		double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
		double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
		
		double ceff_map_duu[12] = { d_uu_8, d_uu_7, d_uu_4, d_uu_5, d_uu_9, d_uu_12, d_uu_11, d_uu_10, d_uu_6, d_uu_3, d_uu_1, d_uu_2 };
		double ceff_map_duv[12] = { d_uv_8, d_uv_7, d_uv_4, d_uv_5, d_uv_9, d_uv_12, d_uv_11, d_uv_10, d_uv_6, d_uv_3, d_uv_1, d_uv_2 };
		double ceff_map_dvv[12] = { d_vv_8, d_vv_7, d_vv_4, d_vv_5, d_vv_9, d_vv_12, d_vv_11, d_vv_10, d_vv_6, d_vv_3, d_vv_1, d_vv_2 };

		double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
		double nrm[3]={0,0,0}; 
		
		for( int x = 0; x < ncoords_base; x++ )
		{
			int *cset = theVertices[i].irr_coord_set + e * ncoords_base;
			double *lr = r + cset[x]*3;

			for( int y = 0; y < 12; y++ )
			{
				R[0] += (lr[0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map[y];
				R[1] += (lr[1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map[y];
				R[2] += (lr[2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map[y];
				
				Ru[0] += (lr[0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_du[y] * u_u;
				Ru[1] += (lr[1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_du[y] * u_u;
				Ru[2] += (lr[2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_du[y] * u_u;
				
				Rv[0] += (lr[0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_dv[y] * v_v;
				Rv[1] += (lr[1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_dv[y] * v_v;
				Rv[2] += (lr[2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_dv[y] * v_v;
				
				tSuu[0] += (lr[0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_duu[y] * u_u *u_u;
				tSuu[1] += (lr[1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_duu[y] * u_u *u_u;
				tSuu[2] += (lr[2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_duu[y] * u_u *u_u;
				                                                                                                                
				tSuv[0] += (lr[0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_duv[y] * u_u *v_v;
				tSuv[1] += (lr[1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_duv[y] * u_u *v_v;
				tSuv[2] += (lr[2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_duv[y] * u_u *v_v;
				                                                                                                                
				tSvv[0] += (lr[0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_dvv[y] * v_v *v_v;
				tSvv[1] += (lr[1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_dvv[y] * v_v *v_v;
				tSvv[2] += (lr[2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_dvv[y] * v_v *v_v;
				
			}
		}

		double dedr[3*np];
		memset( dedr, 0, sizeof(double) * 3 *np );
		
		cross( Ru, Rv, nrm );
		normalize(nrm);

		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		double g = sqrt(RuRu*RvRv-RuRv*RuRv);

double d_nrmz_d_rux=0,d_nrmz_d_ruy=0,d_nrmz_d_ruz=0,d_nrmz_d_rvx=0,d_nrmz_d_rvy=0,d_nrmz_d_rvz=0,d_nrmy_d_rux=0,d_nrmy_d_ruy=0,d_nrmy_d_ruz=0,d_nrmy_d_rvx=0,d_nrmy_d_rvy=0,d_nrmy_d_rvz=0,d_nrmx_d_rux=0,d_nrmx_d_ruy=0,d_nrmx_d_ruz=0,d_nrmx_d_rvx=0,d_nrmx_d_rvy=0,d_nrmx_d_rvz=0,d_RuRu_d_rux=0,d_RuRu_d_ruy=0,d_RuRu_d_ruz=0,d_nsvv_d_rux=0,d_nsvv_d_ruy=0,d_nsvv_d_ruz=0,d_nsvv_d_rvx=0,d_nsvv_d_rvy=0,d_nsvv_d_rvz=0,d_nsvv_d_svvx=0,d_nsvv_d_svvy=0,d_nsvv_d_svvz=0,d_RuRv_d_rux=0,d_RuRv_d_ruy=0,d_RuRv_d_ruz=0,d_RuRv_d_rvx=0,d_RuRv_d_rvy=0,d_RuRv_d_rvz=0,d_nsuv_d_rux=0,d_nsuv_d_ruy=0,d_nsuv_d_ruz=0,d_nsuv_d_rvx=0,d_nsuv_d_rvy=0,d_nsuv_d_rvz=0,d_nsuv_d_suvx=0,d_nsuv_d_suvy=0,d_nsuv_d_suvz=0,d_RvRv_d_rvx=0,d_RvRv_d_rvy=0,d_RvRv_d_rvz=0,d_nsuu_d_rux=0,d_nsuu_d_ruy=0,d_nsuu_d_ruz=0,d_nsuu_d_rvx=0,d_nsuu_d_rvy=0,d_nsuu_d_rvz=0,d_nsuu_d_suux=0,d_nsuu_d_suuy=0,d_nsuu_d_suuz=0,d_g_d_rux=0,d_g_d_ruy=0,d_g_d_ruz=0,d_g_d_rvx=0,d_g_d_rvy=0,d_g_d_rvz=0,d_e_d_rux=0,d_e_d_ruy=0,d_e_d_ruz=0,d_e_d_rvx=0,d_e_d_rvy=0,d_e_d_rvz=0,d_e_d_suux=0,d_e_d_suuy=0,d_e_d_suuz=0,d_e_d_suvx=0,d_e_d_suvy=0,d_e_d_suvz=0,d_e_d_svvx=0,d_e_d_svvy=0,d_e_d_svvz=0,d_c2_d_rux=0,d_c2_d_ruy=0,d_c2_d_ruz=0,d_c2_d_rvx=0,d_c2_d_rvy=0,d_c2_d_rvz=0,d_c2_d_suux=0,d_c2_d_suuy=0,d_c2_d_suuz=0,d_c2_d_suvx=0,d_c2_d_suvy=0,d_c2_d_suvz=0,d_c2_d_svvx=0,d_c2_d_svvy=0,d_c2_d_svvz=0,d_d_d_rux=0,d_d_d_ruy=0,d_d_d_ruz=0,d_d_d_rvx=0,d_d_d_rvy=0,d_d_d_rvz=0,d_d_d_suvx=0,d_d_d_suvy=0,d_d_d_suvz=0,d_d_d_svvx=0,d_d_d_svvy=0,d_d_d_svvz=0,d_c_d_rux=0,d_c_d_ruy=0,d_c_d_ruz=0,d_c_d_rvx=0,d_c_d_rvy=0,d_c_d_rvz=0,d_c_d_suux=0,d_c_d_suuy=0,d_c_d_suuz=0,d_c_d_suvx=0,d_c_d_suvy=0,d_c_d_suvz=0,d_b_d_rux=0,d_b_d_ruy=0,d_b_d_ruz=0,d_b_d_rvx=0,d_b_d_rvy=0,d_b_d_rvz=0,d_b_d_suvx=0,d_b_d_suvy=0,d_b_d_suvz=0,d_b_d_svvx=0,d_b_d_svvy=0,d_b_d_svvz=0,d_a_d_rux=0,d_a_d_ruy=0,d_a_d_ruz=0,d_a_d_rvx=0,d_a_d_rvy=0,d_a_d_rvz=0,d_a_d_suux=0,d_a_d_suuy=0,d_a_d_suuz=0,d_a_d_suvx=0,d_a_d_suvy=0,d_a_d_suvz=0,d_c1_d_rux=0,d_c1_d_ruy=0,d_c1_d_ruz=0,d_c1_d_rvx=0,d_c1_d_rvy=0,d_c1_d_rvz=0,d_c1_d_suux=0,d_c1_d_suuy=0,d_c1_d_suuz=0,d_c1_d_suvx=0,d_c1_d_suvy=0,d_c1_d_suvz=0,d_c1_d_svvx=0,d_c1_d_svvy=0,d_c1_d_svvz=0,d_g_d_RvRv=0,d_g_d_RuRv=0,d_g_d_RuRu=0,d_nsvv_d_nrmz=0,d_nsvv_d_nrmy=0,d_nsvv_d_nrmx=0,d_nsuv_d_nrmz=0,d_nsuv_d_nrmy=0,d_nsuv_d_nrmx=0,d_nsuu_d_nrmz=0,d_nsuu_d_nrmy=0,d_nsuu_d_nrmx=0,d_d_d_RuRv=0,d_d_d_nsuv=0,d_d_d_RuRu=0,d_d_d_nsvv=0,d_d_d_g=0,d_c_d_RuRv=0,d_c_d_nsuv=0,d_c_d_RuRu=0,d_c_d_nsuu=0,d_c_d_g=0,d_b_d_RuRv=0,d_b_d_nsuv=0,d_b_d_RvRv=0,d_b_d_nsvv=0,d_b_d_g=0,d_a_d_RuRv=0,d_a_d_nsuv=0,d_a_d_RvRv=0,d_a_d_nsuu=0,d_a_d_g=0,d_e_d_c2=0,d_e_d_c1=0,d_e_d_g=0,d_c2_d_d=0,d_c2_d_c=0,d_c2_d_b=0,d_c2_d_a=0,d_c1_d_d=0,d_c1_d_c=0,d_c1_d_b=0,d_c1_d_a=0,junk;
		// the basic variables.

		d_g_d_RuRu = RvRv/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
		d_g_d_RuRv = -(RuRv/Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
		d_g_d_RvRv = RuRu/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
		

		d_RuRu_d_rux = 2*Ru[0];
		d_RuRu_d_ruy = 2*Ru[1];
		d_RuRu_d_ruz = 2*Ru[2];
		
		d_RuRv_d_rux = Rv[0];
		d_RuRv_d_ruy = Rv[1];
		d_RuRv_d_ruz = Rv[2];
		
		d_RuRv_d_rvx = Ru[0];
		d_RuRv_d_rvy = Ru[1];
		d_RuRv_d_rvz = Ru[2];
		
		d_RvRv_d_rvx = 2*Rv[0];
		d_RvRv_d_rvy = 2*Rv[1];
		d_RvRv_d_rvz = 2*Rv[2];
		
		d_g_d_rvx += d_g_d_RvRv * d_RvRv_d_rvx;
		d_g_d_rvy += d_g_d_RvRv * d_RvRv_d_rvy;
		d_g_d_rvz += d_g_d_RvRv * d_RvRv_d_rvz;
		d_g_d_rux += d_g_d_RuRv * d_RuRv_d_rux;
		d_g_d_ruy += d_g_d_RuRv * d_RuRv_d_ruy;
		d_g_d_ruz += d_g_d_RuRv * d_RuRv_d_ruz;
		d_g_d_rvx += d_g_d_RuRv * d_RuRv_d_rvx;
		d_g_d_rvy += d_g_d_RuRv * d_RuRv_d_rvy;
		d_g_d_rvz += d_g_d_RuRv * d_RuRv_d_rvz;
		d_g_d_rux += d_g_d_RuRu * d_RuRu_d_rux;
		d_g_d_ruy += d_g_d_RuRu * d_RuRu_d_ruy;
		d_g_d_ruz += d_g_d_RuRu * d_RuRu_d_ruz;

		// that's it.

		double fac = Power(Power(Ru[1]*Rv[0] - Ru[0]*Rv[1],2) + Power(Ru[2]*Rv[0] - Ru[0]*Rv[2],2) + Power(Ru[2]*Rv[1] - Ru[1]*Rv[2],2),1.5);

		d_nrmx_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac; 
		d_nrmx_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2)))/fac;
		d_nrmx_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac;
		
		d_nrmy_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2))))/fac;
		d_nrmy_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac; 
		d_nrmy_d_ruz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2])) + Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac;

		d_nrmz_d_rux = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2])/fac;
		d_nrmz_d_ruy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;
		d_nrmz_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2))) + (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;

		d_nrmx_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac; 
		d_nrmx_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2])))/fac;
		d_nrmx_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac;
		
		d_nrmy_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
		d_nrmy_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2])))/fac;
		d_nrmy_d_rvz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[0]*Ru[1]*Rv[0] - Power(Ru[0],2)*Rv[1] + Ru[2]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2])))/fac;

		d_nrmz_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;
		d_nrmz_d_rvy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[0]*Ru[2]*Rv[0] - Power(Ru[0],2)*Rv[2] + Ru[1]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
		d_nrmz_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;
		
		d_e_d_rux = (de_dnrm[0] * d_nrmx_d_rux + de_dnrm[1] * d_nrmy_d_rux + de_dnrm[2] * d_nrmz_d_rux)*frac_mult;
		d_e_d_ruy = (de_dnrm[0] * d_nrmx_d_ruy + de_dnrm[1] * d_nrmy_d_ruy + de_dnrm[2] * d_nrmz_d_ruy)*frac_mult;
		d_e_d_ruz = (de_dnrm[0] * d_nrmx_d_ruz + de_dnrm[1] * d_nrmy_d_ruz + de_dnrm[2] * d_nrmz_d_ruz)*frac_mult;
		
		d_e_d_rvx = (de_dnrm[0] * d_nrmx_d_rvx + de_dnrm[1] * d_nrmy_d_rvx + de_dnrm[2] * d_nrmz_d_rvx)*frac_mult;
		d_e_d_rvy = (de_dnrm[0] * d_nrmx_d_rvy + de_dnrm[1] * d_nrmy_d_rvy + de_dnrm[2] * d_nrmz_d_rvy)*frac_mult;
		d_e_d_rvz = (de_dnrm[0] * d_nrmx_d_rvz + de_dnrm[1] * d_nrmy_d_rvz + de_dnrm[2] * d_nrmz_d_rvz)*frac_mult;

		double d_e_d_rx = de_dr[0]*frac_mult;
		double d_e_d_ry = de_dr[1]*frac_mult;
		double d_e_d_rz = de_dr[2]*frac_mult;

		// dn/du dn/dv
		double Nu[3]={0,0,0},Nv[3]={0,0,0};

		Nu[0] += d_nrmx_d_rux * tSuu[0] + d_nrmx_d_rvx * tSuv[0];
		Nu[1] += d_nrmy_d_rux * tSuu[0] + d_nrmy_d_rvx * tSuv[0];
		Nu[2] += d_nrmz_d_rux * tSuu[0] + d_nrmz_d_rvx * tSuv[0];
		                                           
		Nu[0] += d_nrmx_d_ruy * tSuu[1] + d_nrmx_d_rvy * tSuv[1];
		Nu[1] += d_nrmy_d_ruy * tSuu[1] + d_nrmy_d_rvy * tSuv[1];
		Nu[2] += d_nrmz_d_ruy * tSuu[1] + d_nrmz_d_rvy * tSuv[1];
                                                           
		Nu[0] += d_nrmx_d_ruz * tSuu[2] + d_nrmx_d_rvz * tSuv[2];
		Nu[1] += d_nrmy_d_ruz * tSuu[2] + d_nrmy_d_rvz * tSuv[2];
		Nu[2] += d_nrmz_d_ruz * tSuu[2] + d_nrmz_d_rvz * tSuv[2];
                                                           
		Nv[0] += d_nrmx_d_rux * tSuv[0] + d_nrmx_d_rvx * tSvv[0];
		Nv[1] += d_nrmy_d_rux * tSuv[0] + d_nrmy_d_rvx * tSvv[0];
		Nv[2] += d_nrmz_d_rux * tSuv[0] + d_nrmz_d_rvx * tSvv[0];
		                                           
		Nv[0] += d_nrmx_d_ruy * tSuv[1] + d_nrmx_d_rvy * tSvv[1];
		Nv[1] += d_nrmy_d_ruy * tSuv[1] + d_nrmy_d_rvy * tSvv[1];
		Nv[2] += d_nrmz_d_ruy * tSuv[1] + d_nrmz_d_rvy * tSvv[1];
                                                           
		Nv[0] += d_nrmx_d_ruz * tSuv[2] + d_nrmx_d_rvz * tSvv[2];
		Nv[1] += d_nrmy_d_ruz * tSuv[2] + d_nrmy_d_rvz * tSvv[2];
		Nv[2] += d_nrmz_d_ruz * tSuv[2] + d_nrmz_d_rvz * tSvv[2];

		nrm_grad[0] = Nu[0];
		nrm_grad[1] = Nu[1];
		nrm_grad[2] = Nu[2];
		
		nrm_grad[3] = Nv[0];
		nrm_grad[4] = Nv[1];
		nrm_grad[5] = Nv[2];

		// dE/du dE/dv
		g_puv[0] += de_dr[0] * Ru[0] + de_dr[1] * Ru[1] + de_dr[2] * Ru[2];
		g_puv[1] += de_dr[0] * Rv[0] + de_dr[1] * Rv[1] + de_dr[2] * Rv[2];
		
		g_puv[0] += de_dnrm[0] * Nu[0] + de_dnrm[1] * Nu[1] + de_dnrm[2] * Nu[2];
		g_puv[1] += de_dnrm[0] * Nv[0] + de_dnrm[1] * Nv[1] + de_dnrm[2] * Nv[2];

	

		// normal is constant wrt to u and v at first order.
	
		for( int p = 0; p < np; p++ )
		for( int y = 0; y < 12; y++ )
		{
			gr[3*cp[p]+0] += d_e_d_rx * ceff_map[y] * theMap[y*ncoords_base+p]* alpha_x;
			gr[3*cp[p]+1] += d_e_d_ry * ceff_map[y] * theMap[y*ncoords_base+p]* alpha_y;
			gr[3*cp[p]+2] += d_e_d_rz * ceff_map[y] * theMap[y*ncoords_base+p]* alpha_z;

			gr[3*cp[p]+0] += d_e_d_rux * ceff_map_du[y] * theMap[y*ncoords_base+p]* alpha_x * u_u; 
			gr[3*cp[p]+1] += d_e_d_ruy * ceff_map_du[y] * theMap[y*ncoords_base+p]* alpha_y * u_u; 
			gr[3*cp[p]+2] += d_e_d_ruz * ceff_map_du[y] * theMap[y*ncoords_base+p]* alpha_z * u_u; 
			
			gr[3*cp[p]+0] += d_e_d_rvx * ceff_map_dv[y] * theMap[y*ncoords_base+p]* alpha_x * v_v; 
			gr[3*cp[p]+1] += d_e_d_rvy * ceff_map_dv[y] * theMap[y*ncoords_base+p]* alpha_y * v_v; 
			gr[3*cp[p]+2] += d_e_d_rvz * ceff_map_dv[y] * theMap[y*ncoords_base+p]* alpha_z * v_v; 
			
			dedr[3*p+0] += d_e_d_rx * ceff_map[y]* theMap[y*ncoords_base+p];
			dedr[3*p+1] += d_e_d_ry * ceff_map[y]* theMap[y*ncoords_base+p];
			dedr[3*p+2] += d_e_d_rz * ceff_map[y]* theMap[y*ncoords_base+p];
			
			dedr[3*p+0] += d_e_d_rux * ceff_map_du[y]* theMap[y*ncoords_base+p] * u_u;
			dedr[3*p+1] += d_e_d_ruy * ceff_map_du[y]* theMap[y*ncoords_base+p] * u_u;
			dedr[3*p+2] += d_e_d_ruz * ceff_map_du[y]* theMap[y*ncoords_base+p] * u_u;
			                                                                         
			dedr[3*p+0] += d_e_d_rvx * ceff_map_dv[y]* theMap[y*ncoords_base+p] * v_v;
			dedr[3*p+1] += d_e_d_rvy * ceff_map_dv[y]* theMap[y*ncoords_base+p] * v_v;
			dedr[3*p+2] += d_e_d_rvz * ceff_map_dv[y]* theMap[y*ncoords_base+p] * v_v;
		}
		
		if ( on_surface )
                  {
                    for( int p = 0; p < np; p++ )
		    for( int y = 0; y < 12; y++ )
                      {

                        gr[3*cp[p]+0] -= (1/g) * kT * d_g_d_rux * ceff_map_du[y] * theMap[y*ncoords_base+p]* alpha_x * u_u;
                        gr[3*cp[p]+1] -= (1/g) * kT * d_g_d_ruy * ceff_map_du[y] * theMap[y*ncoords_base+p]* alpha_y * u_u;
                        gr[3*cp[p]+2] -= (1/g) * kT * d_g_d_ruz * ceff_map_du[y] * theMap[y*ncoords_base+p]* alpha_z * u_u;

                        gr[3*cp[p]+0] -= (1/g) * kT * d_g_d_rvx * ceff_map_dv[y] * theMap[y*ncoords_base+p]* alpha_x * v_v;
                        gr[3*cp[p]+1] -= (1/g) * kT * d_g_d_rvy * ceff_map_dv[y] * theMap[y*ncoords_base+p]* alpha_y * v_v;
                        gr[3*cp[p]+2] -= (1/g) * kT * d_g_d_rvz * ceff_map_dv[y] * theMap[y*ncoords_base+p]* alpha_z * v_v;

                        dedr[3*p+0] -= (1/g) * kT * d_g_d_rux * ceff_map_du[y]* theMap[y*ncoords_base+p] * u_u;
                        dedr[3*p+1] -= (1/g) * kT * d_g_d_ruy * ceff_map_du[y]* theMap[y*ncoords_base+p] * u_u;
                        dedr[3*p+2] -= (1/g) * kT * d_g_d_ruz * ceff_map_du[y]* theMap[y*ncoords_base+p] * u_u;

                        dedr[3*p+0] -= (1/g) * kT * d_g_d_rvx * ceff_map_dv[y]* theMap[y*ncoords_base+p] * v_v;
                        dedr[3*p+1] -= (1/g) * kT * d_g_d_rvy * ceff_map_dv[y]* theMap[y*ncoords_base+p] * v_v;
                        dedr[3*p+2] -= (1/g) * kT * d_g_d_rvz * ceff_map_dv[y]* theMap[y*ncoords_base+p] * v_v;

                      }
                  }

		for( int p = 0; p < np; p++ )
		{
			gr[3*nv+0] += dedr[3*p+0] * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
			gr[3*nv+1] += dedr[3*p+1] * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
			gr[3*nv+2] += dedr[3*p+2] * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
		}
	}

}

void surface::igrad( double *r, double *gr )
{	
	if( !setup_for_parallel )
		setupParallel(this,NULL,0,-1);

	double alpha_x = r[3*nv];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];
	if( !theIrregularFormulas )
		generatePlan();

	double r_val = 0;
	double e = 0;
	double area = 0;
	double wgt = 0;
		double dudv = 1.0;

//	printf("This gradient needs to be corrected to use the total area on an irregular face instead of the individual integration point area.\n");
//	exit(1);

	for( int fx = 0; fx < par_info.nf; fx++ )
	{
		int tf = par_info.faces[fx];

		if( tf < nf_faces )
			continue;

		int f = tf - nf_faces;


		double p_face_area = 0;

		int frm = f*nf_irr_pts;
		int t = theIrregularFormulas[frm].tri;
	
		for( int px = 0; px < theTriangles[t].np; px++ )
		{
			if( fabs(theTriangles[t].pc0[px] - theIrregularFormulas[frm].c0) < 1e-7 )
				continue; 
			p_face_area += theTriangles[t].pa[px];
		}
			if( debug_bit )
			{
			if( f != do_face )
				continue;
			}	
		double c0 = 0;
		double c1 = 0;
		double face_area = 0;
		double energy_density = 0;
//		double A=0,A0=0; would prefer to do it this way.

		for( int p = 0; p < nf_irr_pts; p++ )
		{
			int frm = f*nf_irr_pts+p;
			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
			int *cp = theIrregularFormulas[f*nf_irr_pts+p].cp;
			int np = theIrregularFormulas[f*nf_irr_pts+p].ncoor;
			double g0 = theIrregularFormulas[frm].g0;//2*f_g0[nf_faces+f];
			for( int p = 0; p < np; p++ )
			{

				R[0] += theIrregularFormulas[frm].r_w[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				R[1] += theIrregularFormulas[frm].r_w[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				R[2] += theIrregularFormulas[frm].r_w[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				Ru[0] += theIrregularFormulas[frm].r_u[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				Ru[1] += theIrregularFormulas[frm].r_u[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				Ru[2] += theIrregularFormulas[frm].r_u[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				Rv[0] += theIrregularFormulas[frm].r_v[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				Rv[1] += theIrregularFormulas[frm].r_v[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				Rv[2] += theIrregularFormulas[frm].r_v[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSuu[0] += theIrregularFormulas[frm].r_uu[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSuu[1] += theIrregularFormulas[frm].r_uu[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSuu[2] += theIrregularFormulas[frm].r_uu[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSuv[0] += theIrregularFormulas[frm].r_uv[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSuv[1] += theIrregularFormulas[frm].r_uv[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSuv[2] += theIrregularFormulas[frm].r_uv[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSvv[0] += theIrregularFormulas[frm].r_vv[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSvv[1] += theIrregularFormulas[frm].r_vv[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSvv[2] += theIrregularFormulas[frm].r_vv[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
			}
		
		double nrm[3];
		cross( Ru, Rv, nrm );
		normalize(nrm);

		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		double g = sqrt(RuRu*RvRv-RuRv*RuRv);

		double nsuu = tSuu[0] * nrm[0] + tSuu[1] * nrm[1] + tSuu[2] * nrm[2];
		double nsuv = tSuv[0] * nrm[0] + tSuv[1] * nrm[1] + tSuv[2] * nrm[2];
		double nsvv = tSvv[0] * nrm[0] + tSvv[1] * nrm[1] + tSvv[2] * nrm[2];

		double Stot = (nsuu * RvRv + nsvv * RuRu -2*nsuv*RuRv)/(g*g);

//		Int[i*n+j] = 0.5 * kc * Stot * Stot * g; 
//		AInt[i*n+j] = g;

		double Sop[4] = { 1.0/(g*g) * (nsuu * RvRv  - nsuv * RuRv), (1.0/(g*g)) * ( nsuv*RvRv-nsvv*RuRv),
				  1.0/(g*g) * (nsuv * RuRu  - nsuu * RuRv), (1.0/(g*g)) * ( nsvv*RuRu-nsuv*RuRv) };

		double a = Sop[0];		
		double b = Sop[1];
		double c = Sop[2];
		double d = Sop[3];

		double e1 = -0.5*(a+d-sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
		double e2 = -0.5*(a+d+sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
			
		double c0 = theIrregularFormulas[frm].c0;
//		printf("e1: %lf e2: %lf\n", e1, e2 );
		double en = 0.5 * kc * (e1 + e2 - c0 ) * (e1+e2 - c0);
		double dA = g * theIrregularFormulas[frm].weight;

//		would prefer to do it this way, especially for irregular faces!!		
//			A +=  g * theIrregularFormulas[frm].weight;
//			A0 += g0 * theIrregularFormulas[frm].weight;

			face_area  += dudv * g                   * theIrregularFormulas[frm].weight;	
			energy_density += 0.5 * kc * ( e1 + e2 - c0 ) * ( e1 + e2 - c0 ) * g * theIrregularFormulas[frm].weight;
		}

		

		for( int p = 0; p < nf_irr_pts; p++ )
		{
			if( debug_bit )
			{
			if( f != do_face || p != do_pt )
				continue;
			}	
			int frm = f*nf_irr_pts+p;
			double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
			double nrm[3]={0,0,0}; 

			int *cp = theIrregularFormulas[f*nf_irr_pts+p].cp;
			int np = theIrregularFormulas[f*nf_irr_pts+p].ncoor;
			double dedr[3*np];
			memset( dedr, 0, sizeof(double) * 3 *np );
			double A=0, A0=0;
			double g0 = theIrregularFormulas[frm].g0;//2*f_g0[nf_faces+f];
			double RuRv0 = theIrregularFormulas[frm].RuRv0;
			for( int p = 0; p < np; p++ )
			{
				R[0] += theIrregularFormulas[frm].r_w[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				R[1] += theIrregularFormulas[frm].r_w[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				R[2] += theIrregularFormulas[frm].r_w[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				Ru[0] += theIrregularFormulas[frm].r_u[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				Ru[1] += theIrregularFormulas[frm].r_u[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				Ru[2] += theIrregularFormulas[frm].r_u[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				Rv[0] += theIrregularFormulas[frm].r_v[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				Rv[1] += theIrregularFormulas[frm].r_v[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				Rv[2] += theIrregularFormulas[frm].r_v[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSuu[0] += theIrregularFormulas[frm].r_uu[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSuu[1] += theIrregularFormulas[frm].r_uu[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSuu[2] += theIrregularFormulas[frm].r_uu[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSuv[0] += theIrregularFormulas[frm].r_uv[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSuv[1] += theIrregularFormulas[frm].r_uv[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSuv[2] += theIrregularFormulas[frm].r_uv[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
				
				tSvv[0] += theIrregularFormulas[frm].r_vv[p] * alpha_x * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				tSvv[1] += theIrregularFormulas[frm].r_vv[p] * alpha_y * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				tSvv[2] += theIrregularFormulas[frm].r_vv[p] * alpha_z * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
			}
		



		cross( Ru, Rv, nrm );
		normalize(nrm);

		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		double g = sqrt(RuRu*RvRv-RuRv*RuRv);

		double nsuu = tSuu[0] * nrm[0] + tSuu[1] * nrm[1] + tSuu[2] * nrm[2];
		double nsuv = tSuv[0] * nrm[0] + tSuv[1] * nrm[1] + tSuv[2] * nrm[2];
		double nsvv = tSvv[0] * nrm[0] + tSvv[1] * nrm[1] + tSvv[2] * nrm[2];

		double Stot = (nsuu * RvRv + nsvv * RuRu -2*nsuv*RuRv)/(g*g);

//		Int[i*n+j] = 0.5 * kc * Stot * Stot * g; 
//		AInt[i*n+j] = g;

		double Sop[4] = { 1.0/(g*g) * (nsuu * RvRv  - nsuv * RuRv), (1.0/(g*g)) * ( nsuv*RvRv-nsvv*RuRv),
				  1.0/(g*g) * (nsuv * RuRu  - nsuu * RuRv), (1.0/(g*g)) * ( nsvv*RuRu-nsuv*RuRv) };

		double a = Sop[0];		
		double b = Sop[1];
		double c = Sop[2];
		double d = Sop[3];

		double e1 = -0.5*(a+d-sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
		double e2 = -0.5*(a+d+sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
			
		double c0 = theIrregularFormulas[frm].c0;
//		printf("e1: %lf e2: %lf\n", e1, e2 );
		double en = 0.5 * kc * (e1 + e2 - c0 ) * (e1+e2 - c0);
		double dA = g * theIrregularFormulas[frm].weight;
		
		A =  g * theIrregularFormulas[frm].weight;
		A0 = g0 * theIrregularFormulas[frm].weight;

		r_val += dudv * dA * (e1+e2)/2;
		area += dudv * dA;
		e += dudv * dA * en; 			
double d_nrmz_d_rux=0,d_nrmz_d_ruy=0,d_nrmz_d_ruz=0,d_nrmz_d_rvx=0,d_nrmz_d_rvy=0,d_nrmz_d_rvz=0,d_nrmy_d_rux=0,d_nrmy_d_ruy=0,d_nrmy_d_ruz=0,d_nrmy_d_rvx=0,d_nrmy_d_rvy=0,d_nrmy_d_rvz=0,d_nrmx_d_rux=0,d_nrmx_d_ruy=0,d_nrmx_d_ruz=0,d_nrmx_d_rvx=0,d_nrmx_d_rvy=0,d_nrmx_d_rvz=0,d_RuRu_d_rux=0,d_RuRu_d_ruy=0,d_RuRu_d_ruz=0,d_nsvv_d_rux=0,d_nsvv_d_ruy=0,d_nsvv_d_ruz=0,d_nsvv_d_rvx=0,d_nsvv_d_rvy=0,d_nsvv_d_rvz=0,d_nsvv_d_svvx=0,d_nsvv_d_svvy=0,d_nsvv_d_svvz=0,d_RuRv_d_rux=0,d_RuRv_d_ruy=0,d_RuRv_d_ruz=0,d_RuRv_d_rvx=0,d_RuRv_d_rvy=0,d_RuRv_d_rvz=0,d_nsuv_d_rux=0,d_nsuv_d_ruy=0,d_nsuv_d_ruz=0,d_nsuv_d_rvx=0,d_nsuv_d_rvy=0,d_nsuv_d_rvz=0,d_nsuv_d_suvx=0,d_nsuv_d_suvy=0,d_nsuv_d_suvz=0,d_RvRv_d_rvx=0,d_RvRv_d_rvy=0,d_RvRv_d_rvz=0,d_nsuu_d_rux=0,d_nsuu_d_ruy=0,d_nsuu_d_ruz=0,d_nsuu_d_rvx=0,d_nsuu_d_rvy=0,d_nsuu_d_rvz=0,d_nsuu_d_suux=0,d_nsuu_d_suuy=0,d_nsuu_d_suuz=0,d_g_d_rux=0,d_g_d_ruy=0,d_g_d_ruz=0,d_g_d_rvx=0,d_g_d_rvy=0,d_g_d_rvz=0,d_e_d_rux=0,d_e_d_ruy=0,d_e_d_ruz=0,d_e_d_rvx=0,d_e_d_rvy=0,d_e_d_rvz=0,d_e_d_suux=0,d_e_d_suuy=0,d_e_d_suuz=0,d_e_d_suvx=0,d_e_d_suvy=0,d_e_d_suvz=0,d_e_d_svvx=0,d_e_d_svvy=0,d_e_d_svvz=0,d_c2_d_rux=0,d_c2_d_ruy=0,d_c2_d_ruz=0,d_c2_d_rvx=0,d_c2_d_rvy=0,d_c2_d_rvz=0,d_c2_d_suux=0,d_c2_d_suuy=0,d_c2_d_suuz=0,d_c2_d_suvx=0,d_c2_d_suvy=0,d_c2_d_suvz=0,d_c2_d_svvx=0,d_c2_d_svvy=0,d_c2_d_svvz=0,d_d_d_rux=0,d_d_d_ruy=0,d_d_d_ruz=0,d_d_d_rvx=0,d_d_d_rvy=0,d_d_d_rvz=0,d_d_d_suvx=0,d_d_d_suvy=0,d_d_d_suvz=0,d_d_d_svvx=0,d_d_d_svvy=0,d_d_d_svvz=0,d_c_d_rux=0,d_c_d_ruy=0,d_c_d_ruz=0,d_c_d_rvx=0,d_c_d_rvy=0,d_c_d_rvz=0,d_c_d_suux=0,d_c_d_suuy=0,d_c_d_suuz=0,d_c_d_suvx=0,d_c_d_suvy=0,d_c_d_suvz=0,d_b_d_rux=0,d_b_d_ruy=0,d_b_d_ruz=0,d_b_d_rvx=0,d_b_d_rvy=0,d_b_d_rvz=0,d_b_d_suvx=0,d_b_d_suvy=0,d_b_d_suvz=0,d_b_d_svvx=0,d_b_d_svvy=0,d_b_d_svvz=0,d_a_d_rux=0,d_a_d_ruy=0,d_a_d_ruz=0,d_a_d_rvx=0,d_a_d_rvy=0,d_a_d_rvz=0,d_a_d_suux=0,d_a_d_suuy=0,d_a_d_suuz=0,d_a_d_suvx=0,d_a_d_suvy=0,d_a_d_suvz=0,d_c1_d_rux=0,d_c1_d_ruy=0,d_c1_d_ruz=0,d_c1_d_rvx=0,d_c1_d_rvy=0,d_c1_d_rvz=0,d_c1_d_suux=0,d_c1_d_suuy=0,d_c1_d_suuz=0,d_c1_d_suvx=0,d_c1_d_suvy=0,d_c1_d_suvz=0,d_c1_d_svvx=0,d_c1_d_svvy=0,d_c1_d_svvz=0,d_g_d_RvRv=0,d_g_d_RuRv=0,d_g_d_RuRu=0,d_nsvv_d_nrmz=0,d_nsvv_d_nrmy=0,d_nsvv_d_nrmx=0,d_nsuv_d_nrmz=0,d_nsuv_d_nrmy=0,d_nsuv_d_nrmx=0,d_nsuu_d_nrmz=0,d_nsuu_d_nrmy=0,d_nsuu_d_nrmx=0,d_d_d_RuRv=0,d_d_d_nsuv=0,d_d_d_RuRu=0,d_d_d_nsvv=0,d_d_d_g=0,d_c_d_RuRv=0,d_c_d_nsuv=0,d_c_d_RuRu=0,d_c_d_nsuu=0,d_c_d_g=0,d_b_d_RuRv=0,d_b_d_nsuv=0,d_b_d_RvRv=0,d_b_d_nsvv=0,d_b_d_g=0,d_a_d_RuRv=0,d_a_d_nsuv=0,d_a_d_RvRv=0,d_a_d_nsuu=0,d_a_d_g=0,d_e_d_c2=0,d_e_d_c1=0,d_e_d_g=0,d_c2_d_d=0,d_c2_d_c=0,d_c2_d_b=0,d_c2_d_a=0,d_c1_d_d=0,d_c1_d_c=0,d_c1_d_b=0,d_c1_d_a=0,junk;

			// the basic variables.

			d_g_d_RuRu = RvRv/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
			d_g_d_RuRv = -(RuRv/Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
			d_g_d_RvRv = RuRu/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  

			d_RuRu_d_rux = 2*Ru[0];
			d_RuRu_d_ruy = 2*Ru[1];
			d_RuRu_d_ruz = 2*Ru[2];
			
			d_RuRv_d_rux = Rv[0];
			d_RuRv_d_ruy = Rv[1];
			d_RuRv_d_ruz = Rv[2];
			
			d_RuRv_d_rvx = Ru[0];
			d_RuRv_d_rvy = Ru[1];
			d_RuRv_d_rvz = Ru[2];
			
			d_RvRv_d_rvx = 2*Rv[0];
			d_RvRv_d_rvy = 2*Rv[1];
			d_RvRv_d_rvz = 2*Rv[2];

			// intermediates.
			d_a_d_g = (-2*(-(nsuv*RuRv) + nsuu*RvRv))/Power(g,3);
			d_a_d_nsuu = RvRv/Power(g,2);
			d_a_d_RvRv = nsuu/Power(g,2);
			d_a_d_nsuv = -(RuRv/Power(g,2));
			d_a_d_RuRv = -(nsuv/Power(g,2));

			d_b_d_g = (-2*(-(nsvv*RuRv) + nsuv*RvRv))/Power(g,3);
			d_b_d_nsvv = -(RuRv/Power(g,2));
			d_b_d_RvRv = nsuv/Power(g,2);
			d_b_d_nsuv = RvRv/Power(g,2);
			d_b_d_RuRv = -(nsvv/Power(g,2));

			d_c_d_g = (-2*(nsuv*RuRu - nsuu*RuRv))/Power(g,3);
			d_c_d_nsuu = -(RuRv/Power(g,2));
			d_c_d_RuRu = nsuv/Power(g,2);
			d_c_d_nsuv = RuRu/Power(g,2);
			d_c_d_RuRv = -(nsuu/Power(g,2));

			d_d_d_g = (-2*(nsvv*RuRu - nsuv*RuRv))/Power(g,3);
			d_d_d_nsvv = RuRu/Power(g,2);
			d_d_d_RuRu = nsvv/Power(g,2);
			d_d_d_nsuv = -(RuRv/Power(g,2));
			d_d_d_RuRv = -(nsuv/Power(g,2));

			d_c1_d_a = -0.5*(1 - (2*a - 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
			d_c1_d_b =-(-1.*c)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2)) ;
			d_c1_d_c =-(-1.*b)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2)) ;
			d_c1_d_d =-0.5*(1 - (-2*a + 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
		
			d_c2_d_a = -0.5*(1 + (2*a - 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));
			d_c2_d_b = -(1.*c)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2));
			d_c2_d_c = -(1.*b)/Sqrt(SAFETY+1e-100+Power(a,2) + 4*b*c - 2*a*d + Power(d,2));
			d_c2_d_d = -0.5*(1 + (-2*a + 2*d)/(1e-100+2.*Sqrt(SAFETY+Power(a,2) + 4*b*c - 2*a*d + Power(d,2))));

			//d_e_d_g  = 0.5 * kc * ( e1 + e2 - c0) * ( e1 + e2 - c0) * dudv * theIrregularFormulas[frm].weight;
#ifdef FIXED_A
			d_e_d_g = 0;
			d_e_d_c1 = kc * (e1 + e2 - c0) * dudv * g0 * theIrregularFormulas[frm].weight;
			d_e_d_c2 = kc * (e1 + e2 - c0) * dudv * g0 * theIrregularFormulas[frm].weight;
#else
			d_e_d_g  = (0.5 * kc * ( e1 + e2 - c0) * ( e1 + e2 - c0) + kg * e1 * e2 ) * dudv * theIrregularFormulas[frm].weight;
//			if( g < 0 )
//				d_e_d_g *= -1;

			d_e_d_c1 = kc * (e1 + e2 - c0) * dudv * g * theIrregularFormulas[frm].weight;
			d_e_d_c2 = kc * (e1 + e2 - c0) * dudv * g * theIrregularFormulas[frm].weight;

			d_e_d_c1 += kg * e2 * dudv * g * theIrregularFormulas[frm].weight;
			d_e_d_c2 += kg * e1 * dudv * g * theIrregularFormulas[frm].weight;
#endif

/*			if( g < 0 )
			{
				printf("zero-ing curvature.\n");
				d_e_d_g = 0;
				d_e_d_c1 = 0;
				d_e_d_c2 = 0;
			}
*/			
			d_e_d_g += 2* KA * ((A-A0)/A0) * theIrregularFormulas[frm].weight * dudv;

			
			// derivative of dA in the numerator:
			d_e_d_g  += -p_face_area *  0.5 * kc * (e1+e2-c0)*(e1+e2-c0) * theIrregularFormulas[frm].weight / face_area;
			// derivative of dA in the denominator:
			d_e_d_g  += +p_face_area * energy_density * theIrregularFormulas[frm].weight / face_area / face_area;
			d_e_d_c1 += -p_face_area *  dA *  kc * (e1+e2-c0) / face_area;
			d_e_d_c2 += -p_face_area *  dA *  kc * (e1+e2-c0) / face_area;

//			double fac = Power(Power(Ru[1]*Rv[0] - Ru[0]*Rv[1],2) + Power(Ru[2]*Rv[0] - Ru[0]*Rv[2],2) + Power(Ru[2]*Rv[1] - Ru[1]*Rv[2],2),1.5);
			double fac = Power(Power(Ru[1]*Rv[0] - Ru[0]*Rv[1],2) + Power(Ru[2]*Rv[0] - Ru[0]*Rv[2],2) + Power(Ru[2]*Rv[1] - Ru[1]*Rv[2],2),1.5);

			d_nrmx_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac; 
			d_nrmx_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2)))/fac;
			d_nrmx_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac;
			
			d_nrmy_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2))))/fac;
			d_nrmy_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac; 
			d_nrmy_d_ruz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2])) + Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac;
	
			d_nrmz_d_rux = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2])/fac;
			d_nrmz_d_ruy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;
			d_nrmz_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2))) + (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;

			d_nrmx_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac; 
			d_nrmx_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2])))/fac;
			d_nrmx_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac;
			
			d_nrmy_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
			d_nrmy_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2])))/fac;
			d_nrmy_d_rvz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[0]*Ru[1]*Rv[0] - Power(Ru[0],2)*Rv[1] + Ru[2]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2])))/fac;
	
			d_nrmz_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;
			d_nrmz_d_rvy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[0]*Ru[2]*Rv[0] - Power(Ru[0],2)*Rv[2] + Ru[1]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
			d_nrmz_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;
/*

			d_nrmx_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac;
			d_nrmx_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2)))/fac; 
			d_nrmx_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[1]*Rv[0]*Rv[1] + Ru[2]*Rv[0]*Rv[2] - Ru[0]*(Power(Rv[1],2) + Power(Rv[2],2))))/fac;
			
			d_nrmy_d_rux = -((Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2))))/fac;
			d_nrmy_d_ruy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2]) - Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac;
			d_nrmy_d_ruz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Rv[1]*(Ru[0]*Rv[0] + Ru[2]*Rv[2])) + Ru[1]*(Power(Rv[0],2) + Power(Rv[2],2)))/fac;
	
			d_nrmz_d_rux = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2])/fac;
			d_nrmz_d_ruy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2)) - (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;
			d_nrmz_d_ruz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[2]*(Power(Rv[0],2) + Power(Rv[1],2))) + (Ru[0]*Rv[0] + Ru[1]*Rv[1])*Rv[2]))/fac;

			d_nrmx_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac;
			d_nrmx_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2])))/fac;
			d_nrmx_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Power(Ru[1],2)*Rv[0] - Ru[0]*Ru[1]*Rv[1] + Ru[2]*(Ru[2]*Rv[0] - Ru[0]*Rv[2]))/fac;
			
			d_nrmy_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
			d_nrmy_d_rvy = -((Ru[2]*Rv[0] - Ru[0]*Rv[2])*(-(Ru[0]*Ru[1]*Rv[0]) + Power(Ru[0],2)*Rv[1] + Ru[2]*(Ru[2]*Rv[1] - Ru[1]*Rv[2])))/fac;
			d_nrmy_d_rvz = -((Ru[1]*Rv[0] - Ru[0]*Rv[1])*(Ru[0]*Ru[1]*Rv[0] - Power(Ru[0],2)*Rv[1] + Ru[2]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2])))/fac;
	
			d_nrmz_d_rvx = (Ru[2]*Rv[1] - Ru[1]*Rv[2])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;
			d_nrmz_d_rvy = (Ru[2]*Rv[0] - Ru[0]*Rv[2])*(Ru[0]*Ru[2]*Rv[0] - Power(Ru[0],2)*Rv[2] + Ru[1]*(Ru[2]*Rv[1] - Ru[1]*Rv[2]))/fac;
			d_nrmz_d_rvz = (Ru[1]*Rv[0] - Ru[0]*Rv[1])*(-(Ru[0]*Ru[2]*Rv[0]) + Power(Ru[0],2)*Rv[2] + Ru[1]*(-(Ru[2]*Rv[1]) + Ru[1]*Rv[2]))/fac;
*/
			d_nsuu_d_nrmx = tSuu[0];
			d_nsuu_d_nrmy = tSuu[1];
			d_nsuu_d_nrmz = tSuu[2];
			
			d_nsuv_d_nrmx = tSuv[0];
			d_nsuv_d_nrmy = tSuv[1];
			d_nsuv_d_nrmz = tSuv[2];

			d_nsvv_d_nrmx = tSvv[0];
			d_nsvv_d_nrmy = tSvv[1];
			d_nsvv_d_nrmz = tSvv[2];

			d_nsuu_d_suux = nrm[0];
			d_nsuu_d_suuy = nrm[1];
			d_nsuu_d_suuz = nrm[2];
			
			d_nsuv_d_suvx = nrm[0];
			d_nsuv_d_suvy = nrm[1];
			d_nsuv_d_suvz = nrm[2];

			d_nsvv_d_svvx = nrm[0];
			d_nsvv_d_svvy = nrm[1];
			d_nsvv_d_svvz = nrm[2];



	d_nsvv_d_rux += d_nsvv_d_nrmz * d_nrmz_d_rux;
	d_nsvv_d_ruy += d_nsvv_d_nrmz * d_nrmz_d_ruy;
	d_nsvv_d_ruz += d_nsvv_d_nrmz * d_nrmz_d_ruz;
	d_nsvv_d_rvx += d_nsvv_d_nrmz * d_nrmz_d_rvx;
	d_nsvv_d_rvy += d_nsvv_d_nrmz * d_nrmz_d_rvy;
	d_nsvv_d_rvz += d_nsvv_d_nrmz * d_nrmz_d_rvz;
	d_nsvv_d_rux += d_nsvv_d_nrmy * d_nrmy_d_rux;
	d_nsvv_d_ruy += d_nsvv_d_nrmy * d_nrmy_d_ruy;
	d_nsvv_d_ruz += d_nsvv_d_nrmy * d_nrmy_d_ruz;
	d_nsvv_d_rvx += d_nsvv_d_nrmy * d_nrmy_d_rvx;
	d_nsvv_d_rvy += d_nsvv_d_nrmy * d_nrmy_d_rvy;
	d_nsvv_d_rvz += d_nsvv_d_nrmy * d_nrmy_d_rvz;
	d_nsvv_d_rux += d_nsvv_d_nrmx * d_nrmx_d_rux;
	d_nsvv_d_ruy += d_nsvv_d_nrmx * d_nrmx_d_ruy;
	d_nsvv_d_ruz += d_nsvv_d_nrmx * d_nrmx_d_ruz;
	d_nsvv_d_rvx += d_nsvv_d_nrmx * d_nrmx_d_rvx;
	d_nsvv_d_rvy += d_nsvv_d_nrmx * d_nrmx_d_rvy;
	d_nsvv_d_rvz += d_nsvv_d_nrmx * d_nrmx_d_rvz;
	d_nsuv_d_rux += d_nsuv_d_nrmz * d_nrmz_d_rux;
	d_nsuv_d_ruy += d_nsuv_d_nrmz * d_nrmz_d_ruy;
	d_nsuv_d_ruz += d_nsuv_d_nrmz * d_nrmz_d_ruz;
	d_nsuv_d_rvx += d_nsuv_d_nrmz * d_nrmz_d_rvx;
	d_nsuv_d_rvy += d_nsuv_d_nrmz * d_nrmz_d_rvy;
	d_nsuv_d_rvz += d_nsuv_d_nrmz * d_nrmz_d_rvz;
	d_nsuv_d_rux += d_nsuv_d_nrmy * d_nrmy_d_rux;
	d_nsuv_d_ruy += d_nsuv_d_nrmy * d_nrmy_d_ruy;
	d_nsuv_d_ruz += d_nsuv_d_nrmy * d_nrmy_d_ruz;
	d_nsuv_d_rvx += d_nsuv_d_nrmy * d_nrmy_d_rvx;
	d_nsuv_d_rvy += d_nsuv_d_nrmy * d_nrmy_d_rvy;
	d_nsuv_d_rvz += d_nsuv_d_nrmy * d_nrmy_d_rvz;
	d_nsuv_d_rux += d_nsuv_d_nrmx * d_nrmx_d_rux;
	d_nsuv_d_ruy += d_nsuv_d_nrmx * d_nrmx_d_ruy;
	d_nsuv_d_ruz += d_nsuv_d_nrmx * d_nrmx_d_ruz;
	d_nsuv_d_rvx += d_nsuv_d_nrmx * d_nrmx_d_rvx;
	d_nsuv_d_rvy += d_nsuv_d_nrmx * d_nrmx_d_rvy;
	d_nsuv_d_rvz += d_nsuv_d_nrmx * d_nrmx_d_rvz;
	d_nsuu_d_rux += d_nsuu_d_nrmz * d_nrmz_d_rux;
	d_nsuu_d_ruy += d_nsuu_d_nrmz * d_nrmz_d_ruy;
	d_nsuu_d_ruz += d_nsuu_d_nrmz * d_nrmz_d_ruz;
	d_nsuu_d_rvx += d_nsuu_d_nrmz * d_nrmz_d_rvx;
	d_nsuu_d_rvy += d_nsuu_d_nrmz * d_nrmz_d_rvy;
	d_nsuu_d_rvz += d_nsuu_d_nrmz * d_nrmz_d_rvz;
	d_nsuu_d_rux += d_nsuu_d_nrmy * d_nrmy_d_rux;
	d_nsuu_d_ruy += d_nsuu_d_nrmy * d_nrmy_d_ruy;
	d_nsuu_d_ruz += d_nsuu_d_nrmy * d_nrmy_d_ruz;
	d_nsuu_d_rvx += d_nsuu_d_nrmy * d_nrmy_d_rvx;
	d_nsuu_d_rvy += d_nsuu_d_nrmy * d_nrmy_d_rvy;
	d_nsuu_d_rvz += d_nsuu_d_nrmy * d_nrmy_d_rvz;
	d_nsuu_d_rux += d_nsuu_d_nrmx * d_nrmx_d_rux;
	d_nsuu_d_ruy += d_nsuu_d_nrmx * d_nrmx_d_ruy;
	d_nsuu_d_ruz += d_nsuu_d_nrmx * d_nrmx_d_ruz;
	d_nsuu_d_rvx += d_nsuu_d_nrmx * d_nrmx_d_rvx;
	d_nsuu_d_rvy += d_nsuu_d_nrmx * d_nrmx_d_rvy;
	d_nsuu_d_rvz += d_nsuu_d_nrmx * d_nrmx_d_rvz;
	d_g_d_rvx += d_g_d_RvRv * d_RvRv_d_rvx;
	d_g_d_rvy += d_g_d_RvRv * d_RvRv_d_rvy;
	d_g_d_rvz += d_g_d_RvRv * d_RvRv_d_rvz;
	d_g_d_rux += d_g_d_RuRv * d_RuRv_d_rux;
	d_g_d_ruy += d_g_d_RuRv * d_RuRv_d_ruy;
	d_g_d_ruz += d_g_d_RuRv * d_RuRv_d_ruz;
	d_g_d_rvx += d_g_d_RuRv * d_RuRv_d_rvx;
	d_g_d_rvy += d_g_d_RuRv * d_RuRv_d_rvy;
	d_g_d_rvz += d_g_d_RuRv * d_RuRv_d_rvz;
	d_g_d_rux += d_g_d_RuRu * d_RuRu_d_rux;
	d_g_d_ruy += d_g_d_RuRu * d_RuRu_d_ruy;
	d_g_d_ruz += d_g_d_RuRu * d_RuRu_d_ruz;
	d_d_d_rux += d_d_d_RuRv * d_RuRv_d_rux;
	d_d_d_ruy += d_d_d_RuRv * d_RuRv_d_ruy;
	d_d_d_ruz += d_d_d_RuRv * d_RuRv_d_ruz;
	d_d_d_rvx += d_d_d_RuRv * d_RuRv_d_rvx;
	d_d_d_rvy += d_d_d_RuRv * d_RuRv_d_rvy;
	d_d_d_rvz += d_d_d_RuRv * d_RuRv_d_rvz;
	d_d_d_rux += d_d_d_nsuv * d_nsuv_d_rux;
	d_d_d_ruy += d_d_d_nsuv * d_nsuv_d_ruy;
	d_d_d_ruz += d_d_d_nsuv * d_nsuv_d_ruz;
	d_d_d_rvx += d_d_d_nsuv * d_nsuv_d_rvx;
	d_d_d_rvy += d_d_d_nsuv * d_nsuv_d_rvy;
	d_d_d_rvz += d_d_d_nsuv * d_nsuv_d_rvz;
	d_d_d_suvx += d_d_d_nsuv * d_nsuv_d_suvx;
	d_d_d_suvy += d_d_d_nsuv * d_nsuv_d_suvy;
	d_d_d_suvz += d_d_d_nsuv * d_nsuv_d_suvz;
	d_d_d_rux += d_d_d_RuRu * d_RuRu_d_rux;
	d_d_d_ruy += d_d_d_RuRu * d_RuRu_d_ruy;
	d_d_d_ruz += d_d_d_RuRu * d_RuRu_d_ruz;
	d_d_d_rux += d_d_d_nsvv * d_nsvv_d_rux;
	d_d_d_ruy += d_d_d_nsvv * d_nsvv_d_ruy;
	d_d_d_ruz += d_d_d_nsvv * d_nsvv_d_ruz;
	d_d_d_rvx += d_d_d_nsvv * d_nsvv_d_rvx;
	d_d_d_rvy += d_d_d_nsvv * d_nsvv_d_rvy;
	d_d_d_rvz += d_d_d_nsvv * d_nsvv_d_rvz;
	d_d_d_svvx += d_d_d_nsvv * d_nsvv_d_svvx;
	d_d_d_svvy += d_d_d_nsvv * d_nsvv_d_svvy;
	d_d_d_svvz += d_d_d_nsvv * d_nsvv_d_svvz;
	d_d_d_rux += d_d_d_g * d_g_d_rux;
	d_d_d_ruy += d_d_d_g * d_g_d_ruy;
	d_d_d_ruz += d_d_d_g * d_g_d_ruz;
	d_d_d_rvx += d_d_d_g * d_g_d_rvx;
	d_d_d_rvy += d_d_d_g * d_g_d_rvy;
	d_d_d_rvz += d_d_d_g * d_g_d_rvz;
	d_c_d_rux += d_c_d_RuRv * d_RuRv_d_rux;
	d_c_d_ruy += d_c_d_RuRv * d_RuRv_d_ruy;
	d_c_d_ruz += d_c_d_RuRv * d_RuRv_d_ruz;
	d_c_d_rvx += d_c_d_RuRv * d_RuRv_d_rvx;
	d_c_d_rvy += d_c_d_RuRv * d_RuRv_d_rvy;
	d_c_d_rvz += d_c_d_RuRv * d_RuRv_d_rvz;
	d_c_d_rux += d_c_d_nsuv * d_nsuv_d_rux;
	d_c_d_ruy += d_c_d_nsuv * d_nsuv_d_ruy;
	d_c_d_ruz += d_c_d_nsuv * d_nsuv_d_ruz;
	d_c_d_rvx += d_c_d_nsuv * d_nsuv_d_rvx;
	d_c_d_rvy += d_c_d_nsuv * d_nsuv_d_rvy;
	d_c_d_rvz += d_c_d_nsuv * d_nsuv_d_rvz;
	d_c_d_suvx += d_c_d_nsuv * d_nsuv_d_suvx;
	d_c_d_suvy += d_c_d_nsuv * d_nsuv_d_suvy;
	d_c_d_suvz += d_c_d_nsuv * d_nsuv_d_suvz;
	d_c_d_rux += d_c_d_RuRu * d_RuRu_d_rux;
	d_c_d_ruy += d_c_d_RuRu * d_RuRu_d_ruy;
	d_c_d_ruz += d_c_d_RuRu * d_RuRu_d_ruz;
	d_c_d_rux += d_c_d_nsuu * d_nsuu_d_rux;
	d_c_d_ruy += d_c_d_nsuu * d_nsuu_d_ruy;
	d_c_d_ruz += d_c_d_nsuu * d_nsuu_d_ruz;
	d_c_d_rvx += d_c_d_nsuu * d_nsuu_d_rvx;
	d_c_d_rvy += d_c_d_nsuu * d_nsuu_d_rvy;
	d_c_d_rvz += d_c_d_nsuu * d_nsuu_d_rvz;
	d_c_d_suux += d_c_d_nsuu * d_nsuu_d_suux;
	d_c_d_suuy += d_c_d_nsuu * d_nsuu_d_suuy;
	d_c_d_suuz += d_c_d_nsuu * d_nsuu_d_suuz;
	d_c_d_rux += d_c_d_g * d_g_d_rux;
	d_c_d_ruy += d_c_d_g * d_g_d_ruy;
	d_c_d_ruz += d_c_d_g * d_g_d_ruz;
	d_c_d_rvx += d_c_d_g * d_g_d_rvx;
	d_c_d_rvy += d_c_d_g * d_g_d_rvy;
	d_c_d_rvz += d_c_d_g * d_g_d_rvz;
	d_b_d_rux += d_b_d_RuRv * d_RuRv_d_rux;
	d_b_d_ruy += d_b_d_RuRv * d_RuRv_d_ruy;
	d_b_d_ruz += d_b_d_RuRv * d_RuRv_d_ruz;
	d_b_d_rvx += d_b_d_RuRv * d_RuRv_d_rvx;
	d_b_d_rvy += d_b_d_RuRv * d_RuRv_d_rvy;
	d_b_d_rvz += d_b_d_RuRv * d_RuRv_d_rvz;
	d_b_d_rux += d_b_d_nsuv * d_nsuv_d_rux;
	d_b_d_ruy += d_b_d_nsuv * d_nsuv_d_ruy;
	d_b_d_ruz += d_b_d_nsuv * d_nsuv_d_ruz;
	d_b_d_rvx += d_b_d_nsuv * d_nsuv_d_rvx;
	d_b_d_rvy += d_b_d_nsuv * d_nsuv_d_rvy;
	d_b_d_rvz += d_b_d_nsuv * d_nsuv_d_rvz;
	d_b_d_suvx += d_b_d_nsuv * d_nsuv_d_suvx;
	d_b_d_suvy += d_b_d_nsuv * d_nsuv_d_suvy;
	d_b_d_suvz += d_b_d_nsuv * d_nsuv_d_suvz;
	d_b_d_rvx += d_b_d_RvRv * d_RvRv_d_rvx;
	d_b_d_rvy += d_b_d_RvRv * d_RvRv_d_rvy;
	d_b_d_rvz += d_b_d_RvRv * d_RvRv_d_rvz;
	d_b_d_rux += d_b_d_nsvv * d_nsvv_d_rux;
	d_b_d_ruy += d_b_d_nsvv * d_nsvv_d_ruy;
	d_b_d_ruz += d_b_d_nsvv * d_nsvv_d_ruz;
	d_b_d_rvx += d_b_d_nsvv * d_nsvv_d_rvx;
	d_b_d_rvy += d_b_d_nsvv * d_nsvv_d_rvy;
	d_b_d_rvz += d_b_d_nsvv * d_nsvv_d_rvz;
	d_b_d_svvx += d_b_d_nsvv * d_nsvv_d_svvx;
	d_b_d_svvy += d_b_d_nsvv * d_nsvv_d_svvy;
	d_b_d_svvz += d_b_d_nsvv * d_nsvv_d_svvz;
	d_b_d_rux += d_b_d_g * d_g_d_rux;
	d_b_d_ruy += d_b_d_g * d_g_d_ruy;
	d_b_d_ruz += d_b_d_g * d_g_d_ruz;
	d_b_d_rvx += d_b_d_g * d_g_d_rvx;
	d_b_d_rvy += d_b_d_g * d_g_d_rvy;
	d_b_d_rvz += d_b_d_g * d_g_d_rvz;
	d_a_d_rux += d_a_d_RuRv * d_RuRv_d_rux;
	d_a_d_ruy += d_a_d_RuRv * d_RuRv_d_ruy;
	d_a_d_ruz += d_a_d_RuRv * d_RuRv_d_ruz;
	d_a_d_rvx += d_a_d_RuRv * d_RuRv_d_rvx;
	d_a_d_rvy += d_a_d_RuRv * d_RuRv_d_rvy;
	d_a_d_rvz += d_a_d_RuRv * d_RuRv_d_rvz;
	d_a_d_rux += d_a_d_nsuv * d_nsuv_d_rux;
	d_a_d_ruy += d_a_d_nsuv * d_nsuv_d_ruy;
	d_a_d_ruz += d_a_d_nsuv * d_nsuv_d_ruz;
	d_a_d_rvx += d_a_d_nsuv * d_nsuv_d_rvx;
	d_a_d_rvy += d_a_d_nsuv * d_nsuv_d_rvy;
	d_a_d_rvz += d_a_d_nsuv * d_nsuv_d_rvz;
	d_a_d_suvx += d_a_d_nsuv * d_nsuv_d_suvx;
	d_a_d_suvy += d_a_d_nsuv * d_nsuv_d_suvy;
	d_a_d_suvz += d_a_d_nsuv * d_nsuv_d_suvz;
	d_a_d_rvx += d_a_d_RvRv * d_RvRv_d_rvx;
	d_a_d_rvy += d_a_d_RvRv * d_RvRv_d_rvy;
	d_a_d_rvz += d_a_d_RvRv * d_RvRv_d_rvz;
	d_a_d_rux += d_a_d_nsuu * d_nsuu_d_rux;
	d_a_d_ruy += d_a_d_nsuu * d_nsuu_d_ruy;
	d_a_d_ruz += d_a_d_nsuu * d_nsuu_d_ruz;
	d_a_d_rvx += d_a_d_nsuu * d_nsuu_d_rvx;
	d_a_d_rvy += d_a_d_nsuu * d_nsuu_d_rvy;
	d_a_d_rvz += d_a_d_nsuu * d_nsuu_d_rvz;
	d_a_d_suux += d_a_d_nsuu * d_nsuu_d_suux;
	d_a_d_suuy += d_a_d_nsuu * d_nsuu_d_suuy;
	d_a_d_suuz += d_a_d_nsuu * d_nsuu_d_suuz;
	d_a_d_rux += d_a_d_g * d_g_d_rux;
	d_a_d_ruy += d_a_d_g * d_g_d_ruy;
	d_a_d_ruz += d_a_d_g * d_g_d_ruz;
	d_a_d_rvx += d_a_d_g * d_g_d_rvx;
	d_a_d_rvy += d_a_d_g * d_g_d_rvy;
	d_a_d_rvz += d_a_d_g * d_g_d_rvz;
	d_c1_d_rux += d_c1_d_d * d_d_d_rux;
	d_c1_d_ruy += d_c1_d_d * d_d_d_ruy;
	d_c1_d_ruz += d_c1_d_d * d_d_d_ruz;
	d_c1_d_rvx += d_c1_d_d * d_d_d_rvx;
	d_c1_d_rvy += d_c1_d_d * d_d_d_rvy;
	d_c1_d_rvz += d_c1_d_d * d_d_d_rvz;
	d_c1_d_suvx += d_c1_d_d * d_d_d_suvx;
	d_c1_d_suvy += d_c1_d_d * d_d_d_suvy;
	d_c1_d_suvz += d_c1_d_d * d_d_d_suvz;
	d_c1_d_svvx += d_c1_d_d * d_d_d_svvx;
	d_c1_d_svvy += d_c1_d_d * d_d_d_svvy;
	d_c1_d_svvz += d_c1_d_d * d_d_d_svvz;
	d_c1_d_rux += d_c1_d_c * d_c_d_rux;
	d_c1_d_ruy += d_c1_d_c * d_c_d_ruy;
	d_c1_d_ruz += d_c1_d_c * d_c_d_ruz;
	d_c1_d_rvx += d_c1_d_c * d_c_d_rvx;
	d_c1_d_rvy += d_c1_d_c * d_c_d_rvy;
	d_c1_d_rvz += d_c1_d_c * d_c_d_rvz;
	d_c1_d_suux += d_c1_d_c * d_c_d_suux;
	d_c1_d_suuy += d_c1_d_c * d_c_d_suuy;
	d_c1_d_suuz += d_c1_d_c * d_c_d_suuz;
	d_c1_d_suvx += d_c1_d_c * d_c_d_suvx;
	d_c1_d_suvy += d_c1_d_c * d_c_d_suvy;
	d_c1_d_suvz += d_c1_d_c * d_c_d_suvz;
	d_c1_d_rux += d_c1_d_b * d_b_d_rux;
	d_c1_d_ruy += d_c1_d_b * d_b_d_ruy;
	d_c1_d_ruz += d_c1_d_b * d_b_d_ruz;
	d_c1_d_rvx += d_c1_d_b * d_b_d_rvx;
	d_c1_d_rvy += d_c1_d_b * d_b_d_rvy;
	d_c1_d_rvz += d_c1_d_b * d_b_d_rvz;
	d_c1_d_suvx += d_c1_d_b * d_b_d_suvx;
	d_c1_d_suvy += d_c1_d_b * d_b_d_suvy;
	d_c1_d_suvz += d_c1_d_b * d_b_d_suvz;
	d_c1_d_svvx += d_c1_d_b * d_b_d_svvx;
	d_c1_d_svvy += d_c1_d_b * d_b_d_svvy;
	d_c1_d_svvz += d_c1_d_b * d_b_d_svvz;
	d_c1_d_rux += d_c1_d_a * d_a_d_rux;
	d_c1_d_ruy += d_c1_d_a * d_a_d_ruy;
	d_c1_d_ruz += d_c1_d_a * d_a_d_ruz;
	d_c1_d_rvx += d_c1_d_a * d_a_d_rvx;
	d_c1_d_rvy += d_c1_d_a * d_a_d_rvy;
	d_c1_d_rvz += d_c1_d_a * d_a_d_rvz;
	d_c1_d_suux += d_c1_d_a * d_a_d_suux;
	d_c1_d_suuy += d_c1_d_a * d_a_d_suuy;
	d_c1_d_suuz += d_c1_d_a * d_a_d_suuz;
	d_c1_d_suvx += d_c1_d_a * d_a_d_suvx;
	d_c1_d_suvy += d_c1_d_a * d_a_d_suvy;
	d_c1_d_suvz += d_c1_d_a * d_a_d_suvz;
	d_c2_d_rux += d_c2_d_d * d_d_d_rux;
	d_c2_d_ruy += d_c2_d_d * d_d_d_ruy;
	d_c2_d_ruz += d_c2_d_d * d_d_d_ruz;
	d_c2_d_rvx += d_c2_d_d * d_d_d_rvx;
	d_c2_d_rvy += d_c2_d_d * d_d_d_rvy;
	d_c2_d_rvz += d_c2_d_d * d_d_d_rvz;
	d_c2_d_suvx += d_c2_d_d * d_d_d_suvx;
	d_c2_d_suvy += d_c2_d_d * d_d_d_suvy;
	d_c2_d_suvz += d_c2_d_d * d_d_d_suvz;
	d_c2_d_svvx += d_c2_d_d * d_d_d_svvx;
	d_c2_d_svvy += d_c2_d_d * d_d_d_svvy;
	d_c2_d_svvz += d_c2_d_d * d_d_d_svvz;
	d_c2_d_rux += d_c2_d_c * d_c_d_rux;
	d_c2_d_ruy += d_c2_d_c * d_c_d_ruy;
	d_c2_d_ruz += d_c2_d_c * d_c_d_ruz;
	d_c2_d_rvx += d_c2_d_c * d_c_d_rvx;
	d_c2_d_rvy += d_c2_d_c * d_c_d_rvy;
	d_c2_d_rvz += d_c2_d_c * d_c_d_rvz;
	d_c2_d_suux += d_c2_d_c * d_c_d_suux;
	d_c2_d_suuy += d_c2_d_c * d_c_d_suuy;
	d_c2_d_suuz += d_c2_d_c * d_c_d_suuz;
	d_c2_d_suvx += d_c2_d_c * d_c_d_suvx;
	d_c2_d_suvy += d_c2_d_c * d_c_d_suvy;
	d_c2_d_suvz += d_c2_d_c * d_c_d_suvz;
	d_c2_d_rux += d_c2_d_b * d_b_d_rux;
	d_c2_d_ruy += d_c2_d_b * d_b_d_ruy;
	d_c2_d_ruz += d_c2_d_b * d_b_d_ruz;
	d_c2_d_rvx += d_c2_d_b * d_b_d_rvx;
	d_c2_d_rvy += d_c2_d_b * d_b_d_rvy;
	d_c2_d_rvz += d_c2_d_b * d_b_d_rvz;
	d_c2_d_suvx += d_c2_d_b * d_b_d_suvx;
	d_c2_d_suvy += d_c2_d_b * d_b_d_suvy;
	d_c2_d_suvz += d_c2_d_b * d_b_d_suvz;
	d_c2_d_svvx += d_c2_d_b * d_b_d_svvx;
	d_c2_d_svvy += d_c2_d_b * d_b_d_svvy;
	d_c2_d_svvz += d_c2_d_b * d_b_d_svvz;
	d_c2_d_rux += d_c2_d_a * d_a_d_rux;
	d_c2_d_ruy += d_c2_d_a * d_a_d_ruy;
	d_c2_d_ruz += d_c2_d_a * d_a_d_ruz;
	d_c2_d_rvx += d_c2_d_a * d_a_d_rvx;
	d_c2_d_rvy += d_c2_d_a * d_a_d_rvy;
	d_c2_d_rvz += d_c2_d_a * d_a_d_rvz;
	d_c2_d_suux += d_c2_d_a * d_a_d_suux;
	d_c2_d_suuy += d_c2_d_a * d_a_d_suuy;
	d_c2_d_suuz += d_c2_d_a * d_a_d_suuz;
	d_c2_d_suvx += d_c2_d_a * d_a_d_suvx;
	d_c2_d_suvy += d_c2_d_a * d_a_d_suvy;
	d_c2_d_suvz += d_c2_d_a * d_a_d_suvz;
	d_e_d_rux += d_e_d_c2 * d_c2_d_rux;
	d_e_d_ruy += d_e_d_c2 * d_c2_d_ruy;
	d_e_d_ruz += d_e_d_c2 * d_c2_d_ruz;
	d_e_d_rvx += d_e_d_c2 * d_c2_d_rvx;
	d_e_d_rvy += d_e_d_c2 * d_c2_d_rvy;
	d_e_d_rvz += d_e_d_c2 * d_c2_d_rvz;
	d_e_d_suux += d_e_d_c2 * d_c2_d_suux;
	d_e_d_suuy += d_e_d_c2 * d_c2_d_suuy;
	d_e_d_suuz += d_e_d_c2 * d_c2_d_suuz;
	d_e_d_suvx += d_e_d_c2 * d_c2_d_suvx;
	d_e_d_suvy += d_e_d_c2 * d_c2_d_suvy;
	d_e_d_suvz += d_e_d_c2 * d_c2_d_suvz;
	d_e_d_svvx += d_e_d_c2 * d_c2_d_svvx;
	d_e_d_svvy += d_e_d_c2 * d_c2_d_svvy;
	d_e_d_svvz += d_e_d_c2 * d_c2_d_svvz;
	d_e_d_rux += d_e_d_c1 * d_c1_d_rux;
	d_e_d_ruy += d_e_d_c1 * d_c1_d_ruy;
	d_e_d_ruz += d_e_d_c1 * d_c1_d_ruz;
	d_e_d_rvx += d_e_d_c1 * d_c1_d_rvx;
	d_e_d_rvy += d_e_d_c1 * d_c1_d_rvy;
	d_e_d_rvz += d_e_d_c1 * d_c1_d_rvz;
	d_e_d_suux += d_e_d_c1 * d_c1_d_suux;
	d_e_d_suuy += d_e_d_c1 * d_c1_d_suuy;
	d_e_d_suuz += d_e_d_c1 * d_c1_d_suuz;
	d_e_d_suvx += d_e_d_c1 * d_c1_d_suvx;
	d_e_d_suvy += d_e_d_c1 * d_c1_d_suvy;
	d_e_d_suvz += d_e_d_c1 * d_c1_d_suvz;
	d_e_d_svvx += d_e_d_c1 * d_c1_d_svvx;
	d_e_d_svvy += d_e_d_c1 * d_c1_d_svvy;
	d_e_d_svvz += d_e_d_c1 * d_c1_d_svvz;
	d_e_d_rux += d_e_d_g * d_g_d_rux;
	d_e_d_ruy += d_e_d_g * d_g_d_ruy;
	d_e_d_ruz += d_e_d_g * d_g_d_ruz;
	d_e_d_rvx += d_e_d_g * d_g_d_rvx;
	d_e_d_rvy += d_e_d_g * d_g_d_rvy;
	d_e_d_rvz += d_e_d_g * d_g_d_rvz;

/*
if( fabs(d_nsvv_d_nrmz) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_nsvv_d_nrmy) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_nsvv_d_nrmx) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_nsuv_d_nrmz) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_nsuv_d_nrmy) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_nsuv_d_nrmx) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_nsuu_d_nrmz) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_nsuu_d_nrmy) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_nsuu_d_nrmx) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_g_d_RvRv) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_g_d_RuRv) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_g_d_RuRu) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_d_d_RuRv) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_d_d_nsuv) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_d_d_RuRu) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_d_d_nsvv) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_d_d_g) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_c_d_RuRv) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_c_d_nsuv) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_c_d_RuRu) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_c_d_nsuu) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_c_d_g) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_b_d_RuRv) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_b_d_nsuv) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_b_d_RvRv) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_b_d_nsvv) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_b_d_g) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_a_d_RuRv) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_a_d_nsuv) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_a_d_RvRv) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_a_d_nsuu) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_a_d_g) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_c1_d_d) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_c1_d_c) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_c1_d_b) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_c1_d_a) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_c2_d_d) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_c2_d_c) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_c2_d_b) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_c2_d_a) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_e_d_c2) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_e_d_c1) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
if( fabs(d_e_d_g) < 1e-10 ) { printf("Possible error here.\n"); exit(1); }
*/
			double val = RuRv / sqrt(RuRu*RvRv);
			double d_val_d_rux = d_RuRv_d_rux / sqrt(RuRu*RvRv) - RuRv*RvRv/(2*pow(RuRu*RvRv,3.0/2.0)) * d_RuRu_d_rux;
			double d_val_d_ruy = d_RuRv_d_ruy / sqrt(RuRu*RvRv) - RuRv*RvRv/(2*pow(RuRu*RvRv,3.0/2.0)) * d_RuRu_d_ruy;
			double d_val_d_ruz = d_RuRv_d_ruz / sqrt(RuRu*RvRv) - RuRv*RvRv/(2*pow(RuRu*RvRv,3.0/2.0)) * d_RuRu_d_ruz;
                                                                                                                                                                                      
			double d_val_d_rvx = d_RuRv_d_rvx / sqrt(RuRu*RvRv) - RuRu*RuRv/(2*pow(RuRu*RvRv,3.0/2.0)) * d_RvRv_d_rvx;
			double d_val_d_rvy = d_RuRv_d_rvy / sqrt(RuRu*RvRv) - RuRu*RuRv/(2*pow(RuRu*RvRv,3.0/2.0)) * d_RvRv_d_rvy;
			double d_val_d_rvz = d_RuRv_d_rvz / sqrt(RuRu*RvRv) - RuRu*RuRv/(2*pow(RuRu*RvRv,3.0/2.0)) * d_RvRv_d_rvz;


#ifdef DEBUG_G
			for( int p = 0; p < np; p++ )
			{
				gr[3*cp[p]+0] += d_c2_d_rux * theIrregularFormulas[frm].r_u[p]; 
				gr[3*cp[p]+1] += d_c2_d_ruy * theIrregularFormulas[frm].r_u[p]; 
				gr[3*cp[p]+2] += d_c2_d_ruz * theIrregularFormulas[frm].r_u[p]; 
				
				gr[3*cp[p]+0] += d_c2_d_rvx * theIrregularFormulas[frm].r_v[p]; 
				gr[3*cp[p]+1] += d_c2_d_rvy * theIrregularFormulas[frm].r_v[p]; 
				gr[3*cp[p]+2] += d_c2_d_rvz * theIrregularFormulas[frm].r_v[p]; 
				
				gr[3*cp[p]+0] += d_c2_d_suux * theIrregularFormulas[frm].r_uu[p]; 
				gr[3*cp[p]+1] += d_c2_d_suuy * theIrregularFormulas[frm].r_uu[p]; 
				gr[3*cp[p]+2] += d_c2_d_suuz * theIrregularFormulas[frm].r_uu[p]; 
				
				gr[3*cp[p]+0] += d_c2_d_suvx * theIrregularFormulas[frm].r_uv[p]; 
				gr[3*cp[p]+1] += d_c2_d_suvy * theIrregularFormulas[frm].r_uv[p]; 
				gr[3*cp[p]+2] += d_c2_d_suvz * theIrregularFormulas[frm].r_uv[p]; 
				
				gr[3*cp[p]+0] += d_c2_d_svvx * theIrregularFormulas[frm].r_vv[p]; 
				gr[3*cp[p]+1] += d_c2_d_svvy * theIrregularFormulas[frm].r_vv[p]; 
				gr[3*cp[p]+2] += d_c2_d_svvz * theIrregularFormulas[frm].r_vv[p]; 
	
			}
#elif defined(DEBUG_NRMX)
			for( int p = 0; p < np; p++ )
			{
				gr[3*cp[p]+0] += d_nrmz_d_rux * theIrregularFormulas[frm].r_u[p]; 
				gr[3*cp[p]+1] += d_nrmz_d_ruy * theIrregularFormulas[frm].r_u[p]; 
				gr[3*cp[p]+2] += d_nrmz_d_ruz * theIrregularFormulas[frm].r_u[p]; 
				
				gr[3*cp[p]+0] += d_nrmz_d_rvx * theIrregularFormulas[frm].r_v[p]; 
				gr[3*cp[p]+1] += d_nrmz_d_rvy * theIrregularFormulas[frm].r_v[p]; 
				gr[3*cp[p]+2] += d_nrmz_d_rvz * theIrregularFormulas[frm].r_v[p]; 
				
//				gr[3*cp[p]+0] += d_nrmx_d_suux * theIrregularFormulas[frm].r_uu[p]; 
//				gr[3*cp[p]+1] += d_nrmx_d_suuy * theIrregularFormulas[frm].r_uu[p]; 
//				gr[3*cp[p]+2] += d_nrmx_d_suuz * theIrregularFormulas[frm].r_uu[p]; 
				
//				gr[3*cp[p]+0] += d_nsuu_d_suvx * theIrregularFormulas[frm].r_uv[p]; 
//				gr[3*cp[p]+1] += d_nsuu_d_suvy * theIrregularFormulas[frm].r_uv[p]; 
//				gr[3*cp[p]+2] += d_nsuu_d_suvz * theIrregularFormulas[frm].r_uv[p]; 
				
//				gr[3*cp[p]+0] += d_a_d_svvx * theIrregularFormulas[frm].r_vv[p]; 
//				gr[3*cp[p]+1] += d_a_d_svvy * theIrregularFormulas[frm].r_vv[p]; 
//				gr[3*cp[p]+2] += d_a_d_svvz * theIrregularFormulas[frm].r_vv[p]; 
			}

#else
			for( int p = 0; p < np; p++ )
			{
				gr[3*cp[p]+0] += d_e_d_rux * theIrregularFormulas[frm].r_u[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_ruy * theIrregularFormulas[frm].r_u[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_ruz * theIrregularFormulas[frm].r_u[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_rvx * theIrregularFormulas[frm].r_v[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_rvy * theIrregularFormulas[frm].r_v[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_rvz * theIrregularFormulas[frm].r_v[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_suux * theIrregularFormulas[frm].r_uu[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_suuy * theIrregularFormulas[frm].r_uu[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_suuz * theIrregularFormulas[frm].r_uu[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_suvx * theIrregularFormulas[frm].r_uv[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_suvy * theIrregularFormulas[frm].r_uv[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_suvz * theIrregularFormulas[frm].r_uv[p] * alpha_z; 
				
				gr[3*cp[p]+0] += d_e_d_svvx * theIrregularFormulas[frm].r_vv[p] * alpha_x; 
				gr[3*cp[p]+1] += d_e_d_svvy * theIrregularFormulas[frm].r_vv[p] * alpha_y; 
				gr[3*cp[p]+2] += d_e_d_svvz * theIrregularFormulas[frm].r_vv[p] * alpha_z; 
				
				gr[3*cp[p]+0] += k_reg * 2 * d_val_d_rux * (val-RuRv0) * theIrregularFormulas[frm].r_u[p] * alpha_x;
				gr[3*cp[p]+1] += k_reg * 2 * d_val_d_ruy * (val-RuRv0) * theIrregularFormulas[frm].r_u[p] * alpha_y;
				gr[3*cp[p]+2] += k_reg * 2 * d_val_d_ruz * (val-RuRv0) * theIrregularFormulas[frm].r_u[p] * alpha_z;
				                                                   
				gr[3*cp[p]+0] += k_reg * 2 * d_val_d_rvx * (val-RuRv0) * theIrregularFormulas[frm].r_v[p] * alpha_x;
				gr[3*cp[p]+1] += k_reg * 2 * d_val_d_rvy * (val-RuRv0) * theIrregularFormulas[frm].r_v[p] * alpha_y;
				gr[3*cp[p]+2] += k_reg * 2 * d_val_d_rvz * (val-RuRv0) * theIrregularFormulas[frm].r_v[p] * alpha_z;
				
				dedr[3*p+0] += k_reg * 2 * d_val_d_rux * (val-RuRv0) * theIrregularFormulas[frm].r_u[p] * alpha_x;
				dedr[3*p+1] += k_reg * 2 * d_val_d_ruy * (val-RuRv0) * theIrregularFormulas[frm].r_u[p] * alpha_y;
				dedr[3*p+2] += k_reg * 2 * d_val_d_ruz * (val-RuRv0) * theIrregularFormulas[frm].r_u[p] * alpha_z;
				                                                 
				dedr[3*p+0] += k_reg * 2 * d_val_d_rvx * (val-RuRv0) * theIrregularFormulas[frm].r_v[p] * alpha_x;
				dedr[3*p+1] += k_reg * 2 * d_val_d_rvy * (val-RuRv0) * theIrregularFormulas[frm].r_v[p] * alpha_y;
				dedr[3*p+2] += k_reg * 2 * d_val_d_rvz * (val-RuRv0) * theIrregularFormulas[frm].r_v[p] * alpha_z;
				
				dedr[3*p+0] += d_e_d_rux * theIrregularFormulas[frm].r_u[p]; 
				dedr[3*p+1] += d_e_d_ruy * theIrregularFormulas[frm].r_u[p]; 
				dedr[3*p+2] += d_e_d_ruz * theIrregularFormulas[frm].r_u[p]; 
				
				dedr[3*p+0] += d_e_d_rvx * theIrregularFormulas[frm].r_v[p]; 
				dedr[3*p+1] += d_e_d_rvy * theIrregularFormulas[frm].r_v[p]; 
				dedr[3*p+2] += d_e_d_rvz * theIrregularFormulas[frm].r_v[p]; 
				
				dedr[3*p+0] += d_e_d_suux * theIrregularFormulas[frm].r_uu[p]; 
				dedr[3*p+1] += d_e_d_suuy * theIrregularFormulas[frm].r_uu[p]; 
				dedr[3*p+2] += d_e_d_suuz * theIrregularFormulas[frm].r_uu[p]; 
				
				dedr[3*p+0] += d_e_d_suvx * theIrregularFormulas[frm].r_uv[p]; 
				dedr[3*p+1] += d_e_d_suvy * theIrregularFormulas[frm].r_uv[p]; 
				dedr[3*p+2] += d_e_d_suvz * theIrregularFormulas[frm].r_uv[p]; 
				
				dedr[3*p+0] += d_e_d_svvx * theIrregularFormulas[frm].r_vv[p]; 
				dedr[3*p+1] += d_e_d_svvy * theIrregularFormulas[frm].r_vv[p]; 
				dedr[3*p+2] += d_e_d_svvz * theIrregularFormulas[frm].r_vv[p]; 


			}
#endif
			for( int p = 0; p < np; p++ )
			{
				gr[3*nv+0] += dedr[3*p+0] * (r[cp[p]*3+0] + theIrregularFormulas[frm].r_pbc[3*p+0]); 
				gr[3*nv+1] += dedr[3*p+1] * (r[cp[p]*3+1] + theIrregularFormulas[frm].r_pbc[3*p+1]); 
				gr[3*nv+2] += dedr[3*p+2] * (r[cp[p]*3+2] + theIrregularFormulas[frm].r_pbc[3*p+2]); 
			}

		}	

	} 

	r_val /= area;
}

// g0, dGdu, dGdv, dGduu, dGduv, dGdvv
//

double surface::dG( int f, double u, double v, double grad[6], double *r)
{
	double alpha_x = r[3*nv];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];
	double g = 1;

	if( f < nf_faces )
	{
		int frm = f*nf_g_q_p;
		int t = theFormulas[frm].tri;
		int *cp = theFormulas[frm].cp;

		double w = 1-u-v;

		double u2 = u*u;
		double u3 = u*u*u;
		double u4 = u*u*u*u;
		
		double v2 = v*v;
		double v3 = v*v*v;
		double v4 = v*v*v*v;
		
		double w2 = w*w;
		double w3 = w*w*w;
		double w4 = w*w*w*w;

		double n1 = (1.0/12.0)*(u4+2*u3*v); 
		double n2 = (1.0/12.0)*(u4+2*u3*w); 
		double n3 = (1.0/12.0)*(u4+2*u3*w+6*u3*v+6*u2*v*w+12*u2*v2+6*u*v2*w+6*u*v3+2*v3*w+v4); 
		double n4 = (1.0/12.0)*(6*u4+24*u3*w+24*u2*w2+8*u*w3+w4+24*u3*v+60*u2*v*w+36*u*v*w2+6*v*w3+24*u2*v2+36*u*v2*w+12*v2*w2+8*u*v3+6*v3*w+v4); 
		double n5 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+2*u3*v+6*u2*v*w+6*u*v*w2+2*v*w3); 
		double n6 = (1.0/12.0)*(2*u*v3+v4); 
		double n7 = (1.0/12.0)*(u4+6*u3*w+12*u2*w2+6*u*w3+w4+8*u3*v+36*u2*v*w+36*u*v*w2+8*v*w3+24*u2*v2+60*u*v2*w+24*v2*w2+24*u*v3+24*v3*w+6*v4); 
		double n8 = (1.0/12.0)*(u4+8*u3*w+24*u2*w2+24*u*w3+6*w4+6*u3*v+36*u2*v*w+60*u*v*w2+24*v*w3+12*u2*v2+36*u*v2*w+24*v2*w2+6*u*v3+8*v3*w+v4); 
		double n9 = (1.0/12.0)*(2*u*w3+w4); 
		double n10 = (1.0/12.0)*(2*v3*w+v4); 
		double n11 = (1.0/12.0)*(2*u*w3+w4+6*u*v*w2+6*v*w3+6*u*v2*w+12*v2*w2+2*u*v3+6*v3*w+v4); 
		double n12 = (1.0/12.0)*(w4+2*v*w3);

		double du_1 = Power(u,3)/3. + (Power(u,2)*v)/2.;
		double du_2 = Power(u,2)/2. - Power(u,3)/3. - (Power(u,2)*v)/2.; 	
		double du_3 = Power(u,2)/2. - Power(u,3)/3. + u*v - (Power(u,2)*v)/2. + Power(v,2)/2. - Power(v,3)/6.;
		double du_4 = 0.3333333333333333 + u - Power(u,2) - Power(u,3)/3. + v/2. - u*v - (Power(u,2)*v)/2. - Power(v,2) + Power(v,3)/3.;	
		double du_5 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. - v/2. + Power(u,2)*v + Power(v,2)/2. - Power(v,3)/6.;
		double du_6 = Power(v,3)/6.;	
		double du_7 = 0.16666666666666666 - Power(u,2) + (2*Power(u,3))/3. + v/2. - 2*u*v + Power(u,2)*v - Power(v,2)/2. - Power(v,3)/6.;	
		double du_8 = -2*u + 2*Power(u,2) - Power(u,3)/3. - v + 2*u*v - (Power(u,2)*v)/2. + Power(v,2) - Power(v,3)/6.;	
		double du_9 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. + v/2. - (Power(u,2)*v)/2. - Power(v,2)/2. + Power(v,3)/6.;	
		double du_10 = -Power(v,3)/6.;
		double du_11 = -0.16666666666666666 + Power(u,2)/2. - Power(u,3)/3. - v/2. + u*v - (Power(u,2)*v)/2. + Power(v,3)/3.;	
		double du_12 = 	-0.3333333333333333 + u - Power(u,2) + Power(u,3)/3. + v/2. - u*v + (Power(u,2)*v)/2. - Power(v,3)/6.;

		double dv_1 = Power(u,3)/6.;
		double dv_2 = -Power(u,3)/6.;
		double dv_3 = Power(u,2)/2. - Power(u,3)/6. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_4 = 0.16666666666666666 + u/2. - Power(u,2)/2. - Power(u,3)/6. - 2*u*v - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
		double dv_5 = -0.16666666666666666 - u/2. + Power(u,3)/3. + u*v + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_6 = (u*Power(v,2))/2. + Power(v,3)/3.;
		double dv_7 = 0.3333333333333333 + u/2. - Power(u,2) + Power(u,3)/3. + v - u*v - Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_8 = -u + Power(u,2) - Power(u,3)/6. - 2*v + 2*u*v + 2*Power(v,2) - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_9 = -0.3333333333333333 + u/2. - Power(u,3)/6. + v - u*v - Power(v,2) + (u*Power(v,2))/2. + Power(v,3)/3.;
		double dv_10 = Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;
		double dv_11 = 0.16666666666666666 - u/2. + Power(u,2)/2. - Power(u,3)/6. - Power(v,2) + u*Power(v,2) + (2*Power(v,3))/3.;
		double dv_12 = -0.16666666666666666 + u/2. - Power(u,2)/2. + Power(u,3)/6. + Power(v,2)/2. - (u*Power(v,2))/2. - Power(v,3)/3.;

		double d_uu_1 = Power(u,2) + u*v;
		double d_uu_2 = u - Power(u,2) - u*v;
		double d_uu_3 = u - Power(u,2) + v - u*v;
		double d_uu_4 = 1 - 2*u - Power(u,2) - v - u*v;
		double d_uu_5 = -2*u + 2*Power(u,2) + 2*u*v;
		double d_uu_6 = 0;
		double d_uu_7 = -2*u + 2*Power(u,2) - 2*v + 2*u*v;
		double d_uu_8 = -2 + 4*u - Power(u,2) + 2*v - u*v;
		double d_uu_9 = u - Power(u,2) - u*v;
		double d_uu_10 = 0;
		double d_uu_11 = u - Power(u,2) + v - u*v;
		double d_uu_12 = 1 - 2*u + Power(u,2) - v + u*v;
		
		double d_uv_1 = Power(u,2)/2.;
		double d_uv_2 = -Power(u,2)/2.;
		double d_uv_3 = u - Power(u,2)/2. + v - Power(v,2)/2.;
		double d_uv_4 = 0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
		double d_uv_5 = -0.5 + Power(u,2) + v - Power(v,2)/2.;
		double d_uv_6 = Power(v,2)/2.;
		double d_uv_7 = 0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
		double d_uv_8 = -1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
		double d_uv_9 = 0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
		double d_uv_10 = -Power(v,2)/2.;
		double d_uv_11 = -0.5 + u - Power(u,2)/2. + Power(v,2);
		double d_uv_12 = 0.5 - u + Power(u,2)/2. - Power(v,2)/2.;
		
		double d_vv_1 = 0;
		double d_vv_2 = 0;
		double d_vv_3 = u + v - u*v - Power(v,2);
		double d_vv_4 = -2*u - 2*v + 2*u*v + 2*Power(v,2);
		double d_vv_5 = u + v - u*v - Power(v,2);
		double d_vv_6 = u*v + Power(v,2);
		double d_vv_7 = 1 - u - 2*v - u*v - Power(v,2);
		double d_vv_8 = -2 + 2*u + 4*v - u*v - Power(v,2);
		double d_vv_9 = 1 - u - 2*v + u*v + Power(v,2);
		double d_vv_10 = v - u*v - Power(v,2);
		double d_vv_11 = -2*v + 2*u*v + 2*Power(v,2);
		double d_vv_12 = v - u*v - Power(v,2);
		

		double d_uuu_1 = 2*u + v;
		double d_uuu_2 = 1-2*u-v; //u - Power(u,2) - u*v;
		double d_uuu_3 = 1-2*u-v;//u - Power(u,2) + v - u*v;
		double d_uuu_4 = -2-2*u-v;//1 - 2*u - Power(u,2) - v - u*v;
		double d_uuu_5 = -2+4*u+2*v;//-2*u + 2*Power(u,2) + 2*u*v;
		double d_uuu_6 = 0;
		double d_uuu_7 = -2+4*u+2*v;//-2*u + 2*Power(u,2) - 2*v + 2*u*v;
		double d_uuu_8 = 4-2*u-v;//-2 + 4*u - Power(u,2) + 2*v - u*v;
		double d_uuu_9 = 1-2*u-v;//u - Power(u,2) - u*v;
		double d_uuu_10 = 0;
		double d_uuu_11 = 1-2*u-v;//u - Power(u,2) + v - u*v;
		double d_uuu_12 = -2+2*u+v;//1 - 2*u + Power(u,2) - v + u*v;

		double d_uuv_1 = u;//Power(u,2)/2.;
		double d_uuv_2 = -u;//-Power(u,2)/2.;
		double d_uuv_3 = 1-u;//u - Power(u,2)/2. + v - Power(v,2)/2.;
		double d_uuv_4 = -1-u;//0.5 - u - Power(u,2)/2. - 2*v + Power(v,2);
		double d_uuv_5 = 2*u;//-0.5 + Power(u,2) + v - Power(v,2)/2.;
		double d_uuv_6 = 0;//Power(v,2)/2.;
		double d_uuv_7 = -2+2*u;//0.5 - 2*u + Power(u,2) - v - Power(v,2)/2.;
		double d_uuv_8 = 2-u;//-1 + 2*u - Power(u,2)/2. + 2*v - Power(v,2)/2.;
		double d_uuv_9 = -u;//0.5 - Power(u,2)/2. - v + Power(v,2)/2.;
		double d_uuv_10 = 0;//-Power(v,2)/2.;
		double d_uuv_11 = 1-u;//-0.5 + u - Power(u,2)/2. + Power(v,2);
		double d_uuv_12 = -1+u;//0.5 - u + Power(u,2)/2. - Power(v,2)/2.;

		double d_uvv_1 = 0;
		double d_uvv_2 = 0;
		double d_uvv_3 = 1-v;//u + v - u*v - Power(v,2);
		double d_uvv_4 = -2+2*v;//-2*u - 2*v + 2*u*v + 2*Power(v,2);
		double d_uvv_5 = 1-v;//u + v - u*v - Power(v,2);
		double d_uvv_6 = v;//u*v + Power(v,2);
		double d_uvv_7 = -1-v;//1 - u - 2*v - u*v - Power(v,2);
		double d_uvv_8 = 2-v;//-2 + 2*u + 4*v - u*v - Power(v,2);
		double d_uvv_9 = -1+v;//1 - u - 2*v + u*v + Power(v,2);
		double d_uvv_10 = -v;//v - u*v - Power(v,2);
		double d_uvv_11 = 2*v;//-2*v + 2*u*v + 2*Power(v,2);
		double d_uvv_12 = -v;//v - u*v - Power(v,2);
		
		double d_vvv_1 = 0;
		double d_vvv_2 = 0;
		double d_vvv_3 = 1-u-2*v;//u + v - u*v - Power(v,2);
		double d_vvv_4 = -2+2*u+4*v;//-2*u - 2*v + 2*u*v + 2*Power(v,2);
		double d_vvv_5 = 1-u-2*v;//u + v - u*v - Power(v,2);
		double d_vvv_6 = u+2*v;//u*v + Power(v,2);
		double d_vvv_7 = -2-u-2*v;//1 - u - 2*v - u*v - Power(v,2);
		double d_vvv_8 = 4-u-2*v;//-2 + 2*u + 4*v - u*v - Power(v,2);
		double d_vvv_9 = -2+u+2*v;//1 - u - 2*v + u*v + Power(v,2);
		double d_vvv_10 = 1-u-2*v;//v - u*v - Power(v,2);
		double d_vvv_11 = -2+2*u+4*v;//-2*v + 2*u*v + 2*Power(v,2);
		double d_vvv_12 = 1-u-2*v;//v - u*v - Power(v,2);

		double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
		double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
		double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
		
		double ceff_map_duu[12] = { d_uu_8, d_uu_7, d_uu_4, d_uu_5, d_uu_9, d_uu_12, d_uu_11, d_uu_10, d_uu_6, d_uu_3, d_uu_1, d_uu_2 };
		double ceff_map_duv[12] = { d_uv_8, d_uv_7, d_uv_4, d_uv_5, d_uv_9, d_uv_12, d_uv_11, d_uv_10, d_uv_6, d_uv_3, d_uv_1, d_uv_2 };
		double ceff_map_dvv[12] = { d_vv_8, d_vv_7, d_vv_4, d_vv_5, d_vv_9, d_vv_12, d_vv_11, d_vv_10, d_vv_6, d_vv_3, d_vv_1, d_vv_2 };
		
		double ceff_map_duuu[12] = { d_uuu_8, d_uuu_7, d_uuu_4, d_uuu_5, d_uuu_9, d_uuu_12, d_uuu_11, d_uuu_10, d_uuu_6, d_uuu_3, d_uuu_1, d_uuu_2 };
		double ceff_map_duuv[12] = { d_uuv_8, d_uuv_7, d_uuv_4, d_uuv_5, d_uuv_9, d_uuv_12, d_uuv_11, d_uuv_10, d_uuv_6, d_uuv_3, d_uuv_1, d_uuv_2 };
		double ceff_map_duvv[12] = { d_uvv_8, d_uvv_7, d_uvv_4, d_uvv_5, d_uvv_9, d_uvv_12, d_uvv_11, d_uvv_10, d_uvv_6, d_uvv_3, d_uvv_1, d_uvv_2 };
		double ceff_map_dvvv[12] = { d_vvv_8, d_vvv_7, d_vvv_4, d_vvv_5, d_vvv_9, d_vvv_12, d_vvv_11, d_vvv_10, d_vvv_6, d_vvv_3, d_vvv_1, d_vvv_2 };
				
		double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};

		double tSuuu[3] = {0,0,0};
		double tSuuv[3] = {0,0,0};
		double tSuvv[3] = {0,0,0};
		double tSvvv[3] = {0,0,0};

		double nrm[3]={0,0,0}; 

		int np = theFormulas[frm].ncoor;

		for( int p = 0; p < np; p++ )
		{
			R[0] += ceff_map[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			R[1] += ceff_map[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			R[2] += ceff_map[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			Ru[0] += ceff_map_du[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			Ru[1] += ceff_map_du[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			Ru[2] += ceff_map_du[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			Rv[0] += ceff_map_dv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			Rv[1] += ceff_map_dv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			Rv[2] += ceff_map_dv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSuu[0] += ceff_map_duu[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSuu[1] += ceff_map_duu[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSuu[2] += ceff_map_duu[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSuv[0] += ceff_map_duv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSuv[1] += ceff_map_duv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSuv[2] += ceff_map_duv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSvv[0] += ceff_map_dvv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSvv[1] += ceff_map_dvv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSvv[2] += ceff_map_dvv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSuuu[0] += ceff_map_duuu[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSuuu[1] += ceff_map_duuu[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSuuu[2] += ceff_map_duuu[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSuuv[0] += ceff_map_duuv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSuuv[1] += ceff_map_duuv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSuuv[2] += ceff_map_duuv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSuvv[0] += ceff_map_duvv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSuvv[1] += ceff_map_duvv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSuvv[2] += ceff_map_duvv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
			tSvvv[0] += ceff_map_dvvv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
			tSvvv[1] += ceff_map_dvvv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
			tSvvv[2] += ceff_map_dvvv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
		}

		cross( Ru, Rv, nrm );
		normalize(nrm);


		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];

		g = sqrt(RuRu*RvRv-RuRv*RuRv);

double d_nrmz_d_rux=0,d_nrmz_d_ruy=0,d_nrmz_d_ruz=0,d_nrmz_d_rvx=0,d_nrmz_d_rvy=0,d_nrmz_d_rvz=0,d_nrmy_d_rux=0,d_nrmy_d_ruy=0,d_nrmy_d_ruz=0,d_nrmy_d_rvx=0,d_nrmy_d_rvy=0,d_nrmy_d_rvz=0,d_nrmx_d_rux=0,d_nrmx_d_ruy=0,d_nrmx_d_ruz=0,d_nrmx_d_rvx=0,d_nrmx_d_rvy=0,d_nrmx_d_rvz=0,d_RuRu_d_rux=0,d_RuRu_d_ruy=0,d_RuRu_d_ruz=0,d_nsvv_d_rux=0,d_nsvv_d_ruy=0,d_nsvv_d_ruz=0,d_nsvv_d_rvx=0,d_nsvv_d_rvy=0,d_nsvv_d_rvz=0,d_nsvv_d_svvx=0,d_nsvv_d_svvy=0,d_nsvv_d_svvz=0,d_RuRv_d_rux=0,d_RuRv_d_ruy=0,d_RuRv_d_ruz=0,d_RuRv_d_rvx=0,d_RuRv_d_rvy=0,d_RuRv_d_rvz=0,d_nsuv_d_rux=0,d_nsuv_d_ruy=0,d_nsuv_d_ruz=0,d_nsuv_d_rvx=0,d_nsuv_d_rvy=0,d_nsuv_d_rvz=0,d_nsuv_d_suvx=0,d_nsuv_d_suvy=0,d_nsuv_d_suvz=0,d_RvRv_d_rvx=0,d_RvRv_d_rvy=0,d_RvRv_d_rvz=0,d_nsuu_d_rux=0,d_nsuu_d_ruy=0,d_nsuu_d_ruz=0,d_nsuu_d_rvx=0,d_nsuu_d_rvy=0,d_nsuu_d_rvz=0,d_nsuu_d_suux=0,d_nsuu_d_suuy=0,d_nsuu_d_suuz=0,d_g_d_rux=0,d_g_d_ruy=0,d_g_d_ruz=0,d_g_d_rvx=0,d_g_d_rvy=0,d_g_d_rvz=0,d_e_d_rux=0,d_e_d_ruy=0,d_e_d_ruz=0,d_e_d_rvx=0,d_e_d_rvy=0,d_e_d_rvz=0,d_e_d_suux=0,d_e_d_suuy=0,d_e_d_suuz=0,d_e_d_suvx=0,d_e_d_suvy=0,d_e_d_suvz=0,d_e_d_svvx=0,d_e_d_svvy=0,d_e_d_svvz=0,d_c2_d_rux=0,d_c2_d_ruy=0,d_c2_d_ruz=0,d_c2_d_rvx=0,d_c2_d_rvy=0,d_c2_d_rvz=0,d_c2_d_suux=0,d_c2_d_suuy=0,d_c2_d_suuz=0,d_c2_d_suvx=0,d_c2_d_suvy=0,d_c2_d_suvz=0,d_c2_d_svvx=0,d_c2_d_svvy=0,d_c2_d_svvz=0,d_d_d_rux=0,d_d_d_ruy=0,d_d_d_ruz=0,d_d_d_rvx=0,d_d_d_rvy=0,d_d_d_rvz=0,d_d_d_suvx=0,d_d_d_suvy=0,d_d_d_suvz=0,d_d_d_svvx=0,d_d_d_svvy=0,d_d_d_svvz=0,d_c_d_rux=0,d_c_d_ruy=0,d_c_d_ruz=0,d_c_d_rvx=0,d_c_d_rvy=0,d_c_d_rvz=0,d_c_d_suux=0,d_c_d_suuy=0,d_c_d_suuz=0,d_c_d_suvx=0,d_c_d_suvy=0,d_c_d_suvz=0,d_b_d_rux=0,d_b_d_ruy=0,d_b_d_ruz=0,d_b_d_rvx=0,d_b_d_rvy=0,d_b_d_rvz=0,d_b_d_suvx=0,d_b_d_suvy=0,d_b_d_suvz=0,d_b_d_svvx=0,d_b_d_svvy=0,d_b_d_svvz=0,d_a_d_rux=0,d_a_d_ruy=0,d_a_d_ruz=0,d_a_d_rvx=0,d_a_d_rvy=0,d_a_d_rvz=0,d_a_d_suux=0,d_a_d_suuy=0,d_a_d_suuz=0,d_a_d_suvx=0,d_a_d_suvy=0,d_a_d_suvz=0,d_c1_d_rux=0,d_c1_d_ruy=0,d_c1_d_ruz=0,d_c1_d_rvx=0,d_c1_d_rvy=0,d_c1_d_rvz=0,d_c1_d_suux=0,d_c1_d_suuy=0,d_c1_d_suuz=0,d_c1_d_suvx=0,d_c1_d_suvy=0,d_c1_d_suvz=0,d_c1_d_svvx=0,d_c1_d_svvy=0,d_c1_d_svvz=0,d_g_d_RvRv=0,d_g_d_RuRv=0,d_g_d_RuRu=0,d_nsvv_d_nrmz=0,d_nsvv_d_nrmy=0,d_nsvv_d_nrmx=0,d_nsuv_d_nrmz=0,d_nsuv_d_nrmy=0,d_nsuv_d_nrmx=0,d_nsuu_d_nrmz=0,d_nsuu_d_nrmy=0,d_nsuu_d_nrmx=0,d_d_d_RuRv=0,d_d_d_nsuv=0,d_d_d_RuRu=0,d_d_d_nsvv=0,d_d_d_g=0,d_c_d_RuRv=0,d_c_d_nsuv=0,d_c_d_RuRu=0,d_c_d_nsuu=0,d_c_d_g=0,d_b_d_RuRv=0,d_b_d_nsuv=0,d_b_d_RvRv=0,d_b_d_nsvv=0,d_b_d_g=0,d_a_d_RuRv=0,d_a_d_nsuv=0,d_a_d_RvRv=0,d_a_d_nsuu=0,d_a_d_g=0,d_e_d_c2=0,d_e_d_c1=0,d_e_d_g=0,d_c2_d_d=0,d_c2_d_c=0,d_c2_d_b=0,d_c2_d_a=0,d_c1_d_d=0,d_c1_d_c=0,d_c1_d_b=0,d_c1_d_a=0,junk;
		// the basic variables.





		d_g_d_RuRu = RvRv/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
		d_g_d_RuRv = -(RuRv/Sqrt(-Power(RuRv,2) + RuRu*RvRv));  
		d_g_d_RvRv = RuRu/(2.*Sqrt(-Power(RuRv,2) + RuRu*RvRv));  

		double d_RuRu_d_u = 2*Ru[0] * tSuu[0] + 2*Ru[1] * tSuu[1] + 2 * Ru[2] * tSuu[2];
		double d_RuRu_d_v = 2*Ru[0] * tSuv[0] + 2*Ru[1] * tSuv[1] + 2 * Ru[2] * tSuv[2];

		double d_RuRv_d_u =  Ru[0] * tSuv[0] + Ru[1] * tSuv[1] + Ru[2] * tSuv[2];
		       d_RuRv_d_u += Rv[0] * tSuu[0] + Rv[1] * tSuu[1] + Rv[2] * tSuu[2];
		
		double d_RuRv_d_v =  Ru[0] * tSvv[0] + Ru[1] * tSvv[1] + Ru[2] * tSvv[2];
		       d_RuRv_d_v += Rv[0] * tSuv[0] + Rv[1] * tSuv[1] + Rv[2] * tSuv[2];

		double d_RvRv_d_u = 2*Rv[0] * tSuv[0] + 2*Rv[1] * tSuv[1] + 2 * Rv[2] * tSuv[2];
		double d_RvRv_d_v = 2*Rv[0] * tSvv[0] + 2*Rv[1] * tSvv[1] + 2 * Rv[2] * tSvv[2];
		double d_g_d_u = 0;
		double d_g_d_v = 0;

		d_g_d_u += d_g_d_RuRu * 2* Ru[0] * tSuu[0] + d_g_d_RuRu * 2* Ru[1] * tSuu[1] + d_g_d_RuRu * 2* Ru[2] * tSuu[2];
		d_g_d_v += d_g_d_RuRu * 2* Ru[0] * tSuv[0] + d_g_d_RuRu * 2* Ru[1] * tSuv[1] + d_g_d_RuRu * 2* Ru[2] * tSuv[2];

		// d_ru dru_du
		d_g_d_u += d_g_d_RuRv * Rv[0] * tSuu[0] + d_g_d_RuRv * Rv[1] * tSuu[1] + d_g_d_RuRv * Rv[2] * tSuu[2];
		// d_rv_drv_du
		d_g_d_u += d_g_d_RuRv * Ru[0] * tSuv[0] + d_g_d_RuRv * Ru[1] * tSuv[1] + d_g_d_RuRv * Ru[2] * tSuv[2];	
		// d_ru dru_dv
		d_g_d_v += d_g_d_RuRv * Rv[0] * tSuv[0] + d_g_d_RuRv * Rv[1] * tSuv[1] + d_g_d_RuRv * Rv[2] * tSuv[2];
		// d_rv drv_dv
		d_g_d_v += d_g_d_RuRv * Ru[0] * tSvv[0] + d_g_d_RuRv * Ru[1] * tSvv[1] + d_g_d_RuRv * Ru[2] * tSvv[2];
		
		d_g_d_u += d_g_d_RvRv * 2* Rv[0] * tSuv[0] + d_g_d_RvRv * 2* Rv[1] * tSuv[1] + d_g_d_RvRv * 2* Rv[2] * tSuv[2];
		d_g_d_v += d_g_d_RvRv * 2* Rv[0] * tSvv[0] + d_g_d_RvRv * 2* Rv[1] * tSvv[1] + d_g_d_RvRv * 2* Rv[2] * tSvv[2];

		// second derivatives.

		double den = pow( RuRu * RvRv - RuRv * RuRv, 1.5 );
		// diagonal
		double d2_g_d2_RuRu = - RvRv * RvRv / (4 * den );
		double d2_g_d2_RuRv = - RuRu * RvRv / den;
		double d2_g_d2_RvRv = - RuRu * RuRu / (4 * den ); 
		// cross
		double d2_g_d_RuRu_RuRv = RuRv * RvRv / (2*den);
		double d2_g_d_RuRu_RvRv = (-2 * RuRv * RuRv + RvRv * RuRu) / (4*den); 
		double d2_g_d_RuRv_RvRv = RuRu * RuRv / (2*den);


		grad[0] = d_g_d_u;
		grad[1] = d_g_d_v;

		// g_uu
		grad[2] = d2_g_d2_RuRu * d_RuRu_d_u * d_RuRu_d_u + 
			  d2_g_d2_RuRv * d_RuRv_d_u * d_RuRv_d_u +
			  d2_g_d2_RvRv * d_RvRv_d_u * d_RvRv_d_u +
			  2 * d2_g_d_RuRu_RuRv * d_RuRu_d_u * d_RuRv_d_u +
			  2 * d2_g_d_RuRu_RvRv * d_RuRu_d_u * d_RvRv_d_u +
			  2 * d2_g_d_RuRv_RvRv * d_RuRv_d_u * d_RvRv_d_u;
		
		// g_uv
		grad[3] = d2_g_d2_RuRu * d_RuRu_d_u * d_RuRu_d_v + 
			  d2_g_d2_RuRv * d_RuRv_d_u * d_RuRv_d_v +
			  d2_g_d2_RvRv * d_RvRv_d_u * d_RvRv_d_v +
			  d2_g_d_RuRu_RuRv * d_RuRu_d_v * d_RuRv_d_u +
			  d2_g_d_RuRu_RvRv * d_RuRu_d_v * d_RvRv_d_u +
			  d2_g_d_RuRv_RvRv * d_RuRv_d_v * d_RvRv_d_u + 
			  d2_g_d_RuRu_RuRv * d_RuRu_d_u * d_RuRv_d_v +
			  d2_g_d_RuRu_RvRv * d_RuRu_d_u * d_RvRv_d_v +
			  d2_g_d_RuRv_RvRv * d_RuRv_d_u * d_RvRv_d_v;
		
		// g_vv
		grad[4] = d2_g_d2_RuRu * d_RuRu_d_v * d_RuRu_d_v + 
			  d2_g_d2_RuRv * d_RuRv_d_v * d_RuRv_d_v +
			  d2_g_d2_RvRv * d_RvRv_d_v * d_RvRv_d_v +
			  2 * d2_g_d_RuRu_RuRv * d_RuRu_d_v * d_RuRv_d_v +
			  2 * d2_g_d_RuRu_RvRv * d_RuRu_d_v * d_RvRv_d_v +
			  2 * d2_g_d_RuRv_RvRv * d_RuRv_d_v * d_RvRv_d_v;

		double d2_RuRu_d2_u = 
			2*(tSuu[0] * tSuu[0] + tSuu[1] * tSuu[1] + tSuu[2] * tSuu[2]) +
			2*(tSuuu[0] * Ru[0] + tSuuu[1] * Ru[1] + tSuuu[2] * Ru[2]);
		
		double d2_RuRu_d_uv = 
			2*(tSuu[0] * tSuv[0] + tSuu[1] * tSuv[1] + tSuu[2] * tSuv[2]) +
			2*(tSuuv[0] * Ru[0] + tSuuv[1] * Ru[1] + tSuuv[2] * Ru[2]);
		
		double d2_RuRu_d2_v = 
			2*(tSuv[0] * tSuv[0] + tSuv[1] * tSuv[1] + tSuv[2] * tSuv[2]) +
			2*(tSuvv[0] * Ru[0] + tSuvv[1] * Ru[1] + tSuvv[2] * Ru[2]);
		
		double d2_RuRv_d2_u = 
			2*(tSuu[0] * tSuv[0] + tSuu[1] * tSuv[1] + tSuu[2] * tSuv[2]) +
			  (tSuuu[0] * Rv[0] + tSuuu[1] * Rv[1] + tSuuu[2] * Rv[2]) +
			  (Ru[0] * tSuuv[0] + Ru[1] * tSuuv[1] + Ru[2] * tSuuv[2]);
		
		double d2_RuRv_d_uv = 
			 tSuuv[0] * Rv[0]+tSuuv[1] * Rv[1]+tSuuv[2] * Rv[2] +
			 tSuu [0]* tSvv[0]+tSuu [1]* tSvv[1]+tSuu [2]* tSvv[2] +	
			 tSuv [0]* tSuv[0]+tSuv [1]* tSuv[1]+tSuv [2]* tSuv[2] +
			 Ru  [0]*tSuvv[0]+Ru  [1]*tSuvv[1]+Ru  [2]*tSuvv[2];		
		
		double d2_RuRv_d2_v = 
			2*(tSvv[0] * tSuv[0] + tSvv[1] * tSuv[1] + tSvv[2] * tSuv[2]) +
			  (tSvvv[0] * Ru[0] + tSvvv[1] * Ru[1] + tSvvv[2] * Ru[2]) +
			  (Rv[0] * tSuvv[0] + Rv[1] * tSuvv[1] + Rv[2] * tSuvv[2]);
		
		double d2_RvRv_d2_u = 
			2*(tSuv[0] * tSuv[0] + tSuv[1] * tSuv[1] + tSuv[2] * tSuv[2]) +
			2*(tSuuv[0] * Rv[0] + tSuuv[1] * Rv[1] + tSuuv[2] * Rv[2]);
		
		double d2_RvRv_d_uv = 
			2*(tSvv[0] * tSuv[0] + tSvv[1] * tSuv[1] + tSvv[2] * tSuv[2]) +
			2*(tSuvv[0] * Rv[0] + tSuvv[1] * Rv[1] + tSuvv[2] * Rv[2]);
		
		double d2_RvRv_d2_v = 
			2*(tSvv[0] * tSvv[0] + tSvv[1] * tSvv[1] + tSvv[2] * tSvv[2]) +
			2*(tSvvv[0] * Rv[0] + tSvvv[1] * Rv[1] + tSvvv[2] * Rv[2]);
		

		grad[2] += d_g_d_RuRu * d2_RuRu_d2_u; 
		grad[2] += d_g_d_RuRv * d2_RuRv_d2_u; 
		grad[2] += d_g_d_RvRv * d2_RvRv_d2_u; 
		
		grad[3] += d_g_d_RuRu * d2_RuRu_d_uv; 
		grad[3] += d_g_d_RuRv * d2_RuRv_d_uv; 
		grad[3] += d_g_d_RvRv * d2_RvRv_d_uv; 
		
		grad[4] += d_g_d_RuRu * d2_RuRu_d2_v; 
		grad[4] += d_g_d_RuRv * d2_RuRv_d2_v; 
		grad[4] += d_g_d_RvRv * d2_RvRv_d2_v; 
	}
	else
	{	
		grad[0] = -1e100;
		grad[1] = -1e100;
		grad[2] = -1e100;
		grad[3] = -1e100;
		grad[4] = -1e100;
	}

	return g;
}


