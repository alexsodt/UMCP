#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "interp.h"
#include "mutil.h"
#include <math.h>
#include "gsl/gsl_randist.h"
#include "gsl/gsl_rng.h"
#include "random_global.h"
#define MOVE_CHECK

static double THRESH = 1e-8;

void surface::get_pt_coeffs( int f, double u, double v, double *coeffs, int *coord_list, int *ncoords )
{
	if( f >= nf_faces ) // an irrational face.
	{
		int frm = (f-nf_faces)*nf_irr_pts;

		int i = theIrregularFormulas[frm].vertex;
		int e = theIrregularFormulas[frm].edge;
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
		
		if( domain < 0 &&  fu < 1e-10 && fv < 1e-10  )
		{
			// special case.
			int v = i;

			int val = theVertices[i].valence;

			double lr[3] = { 0,0,0};
	
			double *ti = theVertices[i].r;
	
			lr[0] = 0.5 * theVertices[i].r[0]; 
			lr[1] = 0.5 * theVertices[i].r[1]; 
			lr[2] = 0.5 * theVertices[i].r[2]; 
	
			double w = 1.0 / (val * 2);
	
			*ncoords = 1 + theVertices[i].valence;

			coord_list[0] = i;
			coeffs[0] = 0.5;

			for( int e = 0; e < val; e++ )
			{
				int j = theVertices[i].edges[e];

				coord_list[1+e] = j;
				coeffs[1+e] = w;
			}	

			return;

		}
		double *theMap = theKernel->get_map( &fu, &fv );

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

		*ncoords = ncoords_base;
		
		double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
		for( int x = 0; x < ncoords_base; x++ )
		{
			int *cset = theVertices[i].irr_coord_set + e * ncoords_base;

			coeffs[x] = 0;
			coord_list[x] = cset[x];

			for( int y = 0; y < 12; y++ )
				coeffs[x] += theMap[y*ncoords_base+x] * ceff_map[y];
		}
	}
	else
	{
		int frm = f*nf_g_q_p;
	
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
						
		double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
	
		int ncoords_base = theFormulas[frm].ncoor;

		*ncoords = ncoords_base;

		for( int x = 0; x < ncoords_base; x++ )
		{
			coeffs[x] = ceff_map[x];
			coord_list[x] = cp[x];
		}
	
	}	
}

// gets the coefficients, as well as their derivatives.

void surface::get_pt_dcoeffs( int f, double u, double v, double *coeffs, int *coord_list, int *ncoords )
{
	if( f >= nf_faces ) // an irrational face.
	{
		int frm = (f-nf_faces)*nf_irr_pts;

		int i = theIrregularFormulas[frm].vertex;
		int e = theIrregularFormulas[frm].edge;
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
		
		if( domain < 0 &&  fu < 1e-10 && fv < 1e-10  )
		{
			// special case.
			int v = i;

			int val = theVertices[i].valence;

			double lr[3] = { 0,0,0};
	
			double *ti = theVertices[i].r;
	
			lr[0] = 0.5 * theVertices[i].r[0]; 
			lr[1] = 0.5 * theVertices[i].r[1]; 
			lr[2] = 0.5 * theVertices[i].r[2]; 
	
			double w = 1.0 / (val * 2);
	
			*ncoords = 1 + theVertices[i].valence;

			coord_list[0] = i;
			coeffs[0] = 0.5;

			for( int e = 0; e < val; e++ )
			{
				int j = theVertices[i].edges[e];

				coord_list[1+e] = j;
				coeffs[1+e] = w;
			}	

			return;

		}
		double *theMap = theKernel->get_map( &fu, &fv );
		double u_u = u, u_v=0, v_u=0, v_v = v;
		
		theKernel->get_map_transform( &u_u, &u_v, &v_u, &v_v );

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
		
		double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
		double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
		double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };

		*ncoords = ncoords_base;
		
		for( int x = 0; x < ncoords_base; x++ )
		{
			int *cset = theVertices[i].irr_coord_set + e * ncoords_base;

			coord_list[x] = cset[x];
			coeffs[x] = 0;
			coeffs[ncoords_base+x] = 0;
			coeffs[2*ncoords_base+x] = 0;
			for( int y = 0; y < 12; y++ )
			{
				coeffs[x] += theMap[y*ncoords_base+x] * ceff_map[y];
				coeffs[ncoords_base+x] += theMap[y*ncoords_base+x] * ceff_map_du[y] * u_u;
				coeffs[2*ncoords_base+x] += theMap[y*ncoords_base+x] * ceff_map_dv[y] * v_v;
			}
		}
	}
	else
	{
		int frm = f*nf_g_q_p;
	
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
		
		double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
		double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
		double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
	
		int ncoords_base = theFormulas[frm].ncoor;

		*ncoords = ncoords_base;

		for( int x = 0; x < ncoords_base; x++ )
		{
			coeffs[x]                = ceff_map[x];
			coeffs[ncoords_base+x]   = ceff_map_du[x];
			coeffs[2*ncoords_base+x] = ceff_map_dv[x];
			coord_list[x] = cp[x];
		}
	
	}	
}


double surface::g( int f, double u, double v, double *r )
{
	double alpha_x = r[3*nv+0];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];

	if( f >= nf_faces ) // an irrational face.
	{
		int frm = (f-nf_faces)*nf_irr_pts;

		int i = theIrregularFormulas[frm].vertex;
		int e = theIrregularFormulas[frm].edge;
		int val = theVertices[i].valence;

		irr_kernel *theKernel = kernels[val]; 
		int ncoords_base = val + 6;

		double fu = u;
		double fv = v;					

		int domain = theKernel->domain(fu,fv);
		if( domain-1 >= theKernel->ndomains )
			return 1e-100;

		double *theMap = theKernel->get_map( &fu, &fv );

		double u = fu;
		double v = fv;
		double w = 1 - u - v;


		if( u +v > 1.0+THRESH || u < -THRESH || v < -THRESH )
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
		
		double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
		double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
		double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
		
		double lru[3] = { 0,0,0};
		double lrv[3] = { 0,0,0};
		for( int x = 0; x < ncoords_base; x++ )
		{
			int *cset = theVertices[i].irr_coord_set + e * ncoords_base;
			double *lr = r + cset[x]*3;

			for( int y = 0; y < 12; y++ )
			{
				lru[0] += (lr[0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_du[y];
				lru[1] += (lr[1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_du[y];
				lru[2] += (lr[2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_du[y];
				
				lrv[0] += (lr[0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_dv[y];
				lrv[1] += (lr[1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_dv[y];
				lrv[2] += (lr[2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_dv[y];
			}
		}
	
		double nrm[3];	
		cross( lru, lrv, nrm );

		return pow( 4., domain) * normalize(nrm);
	}
	else
	{
		int frm = f*nf_g_q_p;
	
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
		
		double ceff_map[12] = { n8, n7, n4, n5, n9, n12, n11, n10, n6, n3, n1, n2 };
		double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
		double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
	
		double lru[3] = {0,0,0};
		double lrv[3] = {0,0,0};
		int ncoords_base = theFormulas[frm].ncoor;
		for( int x = 0; x < ncoords_base; x++ )
		{
			double *lr = r + cp[x]*3;
	
			lru[0] += ceff_map_du[x] * (lr[0] + theFormulas[frm].r_pbc[3*x+0]);
			lru[1] += ceff_map_du[x] * (lr[1] + theFormulas[frm].r_pbc[3*x+1]);
			lru[2] += ceff_map_du[x] * (lr[2] + theFormulas[frm].r_pbc[3*x+2]);
			
			lrv[0] += ceff_map_dv[x] * (lr[0] + theFormulas[frm].r_pbc[3*x+0]);
			lrv[1] += ceff_map_dv[x] * (lr[1] + theFormulas[frm].r_pbc[3*x+1]);
			lrv[2] += ceff_map_dv[x] * (lr[2] + theFormulas[frm].r_pbc[3*x+2]);
		}
	
		double nrm[3];
		cross( lru, lrv, nrm );
		return normalize(nrm);
	}	

	return 1e-100;
}

void surface::r2der( int f, double u, double v, double *r, double *dr_uu, double *dr_uv, double *dr_vv )
{
	double alpha_x = r[3*nv+0];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];

	memset( dr_uu, 0, sizeof(double) * 3 );
	memset( dr_uv, 0, sizeof(double) * 3 );
	memset( dr_vv, 0, sizeof(double) * 3 );

	if( f >= nf_faces ) // an irrational face.
	{
		int frm = (f-nf_faces)*nf_irr_pts;

		int i = theIrregularFormulas[frm].vertex;
		int e = theIrregularFormulas[frm].edge;
		int val = theVertices[i].valence;

		irr_kernel *theKernel = kernels[val]; 
		int ncoords_base = val + 6;

		double fu = u;
		double fv = v;					

		int domain = theKernel->domain(fu,fv);
		if( domain-1 >= theKernel->ndomains )
		{
			printf("ERROR trying to evaluate ru outside of irregular domain.\n");
			exit(1);
		}
		double u_u = u, u_v=0, v_u=0, v_v = v;
		theKernel->get_map_transform( &u_u, &u_v, &v_u, &v_v );
		double *theMap = theKernel->get_map( &fu, &fv );

		double u = fu;
		double v = fv;
		double w = 1 - u - v;

		if( u +v > 1.0+THRESH || u < -THRESH || v < -THRESH )
		{
			printf("u: %lf v: %lf\n", u, v );
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
		
		
		double ceff_map_duu[12] = { d_uu_8, d_uu_7, d_uu_4, d_uu_5, d_uu_9, d_uu_12, d_uu_11, d_uu_10, d_uu_6, d_uu_3, d_uu_1, d_uu_2 };
		double ceff_map_duv[12] = { d_uv_8, d_uv_7, d_uv_4, d_uv_5, d_uv_9, d_uv_12, d_uv_11, d_uv_10, d_uv_6, d_uv_3, d_uv_1, d_uv_2 };
		double ceff_map_dvv[12] = { d_vv_8, d_vv_7, d_vv_4, d_vv_5, d_vv_9, d_vv_12, d_vv_11, d_vv_10, d_vv_6, d_vv_3, d_vv_1, d_vv_2 };
		
		double tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
	
		int *cset = theVertices[i].irr_coord_set + e * ncoords_base;
		
	
		for( int x = 0; x < ncoords_base; x++ )
		{
			for( int y = 0; y < 12; y++ )
			{
				dr_uu[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_duu[y] * alpha_x * u_u * u_u;
				dr_uu[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_duu[y] * alpha_y * u_u * u_u;
				dr_uu[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_duu[y] * alpha_z * u_u * u_u;
				
				dr_uv[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_duv[y] * alpha_x * u_u * v_v;
				dr_uv[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_duv[y] * alpha_y * u_u * v_v;
				dr_uv[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_duv[y] * alpha_z * u_u * v_v;
				
				dr_vv[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_dvv[y] * alpha_x * v_v * v_v;
				dr_vv[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_dvv[y] * alpha_y * v_v * v_v;
				dr_vv[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_dvv[y] * alpha_z * v_v * v_v;
			}
		}
	}
	else
	{
		int frm = f*nf_g_q_p;
	
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
		
		double ceff_map_duu[12] = { d_uu_8, d_uu_7, d_uu_4, d_uu_5, d_uu_9, d_uu_12, d_uu_11, d_uu_10, d_uu_6, d_uu_3, d_uu_1, d_uu_2 };
		double ceff_map_duv[12] = { d_uv_8, d_uv_7, d_uv_4, d_uv_5, d_uv_9, d_uv_12, d_uv_11, d_uv_10, d_uv_6, d_uv_3, d_uv_1, d_uv_2 };
		double ceff_map_dvv[12] = { d_vv_8, d_vv_7, d_vv_4, d_vv_5, d_vv_9, d_vv_12, d_vv_11, d_vv_10, d_vv_6, d_vv_3, d_vv_1, d_vv_2 };
	
		int ncoords_base = theFormulas[frm].ncoor;
			
		double nrm[3]={0,0,0}; 
		for( int p = 0; p < ncoords_base; p++ )
		{
			double *lr = r + cp[p]*3;
				
				
				dr_uu[0] += ceff_map_duu[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				dr_uu[1] += ceff_map_duu[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				dr_uu[2] += ceff_map_duu[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				dr_uv[0] += ceff_map_duv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				dr_uv[1] += ceff_map_duv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				dr_uv[2] += ceff_map_duv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
				
				dr_vv[0] += ceff_map_dvv[p] * alpha_x*(r[cp[p]*3+0] + theFormulas[frm].r_pbc[3*p+0]); 
				dr_vv[1] += ceff_map_dvv[p] * alpha_y*(r[cp[p]*3+1] + theFormulas[frm].r_pbc[3*p+1]); 
				dr_vv[2] += ceff_map_dvv[p] * alpha_z*(r[cp[p]*3+2] + theFormulas[frm].r_pbc[3*p+2]); 
			
		}	
	}	
}

void surface::ru( int f, double u, double v, double *r, double *dr_u )
{
	double *alphas = r+3*nv;
	dr_u[0] = 0;
	dr_u[1] = 0;
	dr_u[2] = 0;
	double alpha_x = r[3*nv+0];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];

	if( f >= nf_faces ) // an irrational face.
	{
		int frm = (f-nf_faces)*nf_irr_pts;

		int i = theIrregularFormulas[frm].vertex;
		int e = theIrregularFormulas[frm].edge;
		int val = theVertices[i].valence;

		irr_kernel *theKernel = kernels[val]; 
		int ncoords_base = val + 6;

		double u_u = u, u_v=0, v_u=0, v_v = v;
		
		theKernel->get_map_transform( &u_u, &u_v, &v_u, &v_v );

		double fu = u;
		double fv = v;					

		int domain = theKernel->domain(fu,fv);
		if( domain-1 >= theKernel->ndomains )
		{
			printf("ERROR trying to evaluate ru outside of irregular domain.\n");
			exit(1);
		}
		double *theMap = theKernel->get_map( &fu, &fv );

		double u = fu;
		double v = fv;
		double w = 1 - u - v;

		if( u +v > 1.0+THRESH || u < -THRESH || v < -THRESH )
		{
			printf("u: %lf v: %lf\n", u, v );
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

		
		double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
	
		dr_u[0] = 0;	
		dr_u[1] = 0;	
		dr_u[2] = 0;	


		for( int x = 0; x < ncoords_base; x++ )
		{
			int *cset = theVertices[i].irr_coord_set + e * ncoords_base;
			double *lr = r + cset[x]*3;

			for( int y = 0; y < 12; y++ )
			{
				dr_u[0] += (lr[0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_du[y] * u_u*alphas[0];
				dr_u[1] += (lr[1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_du[y] * u_u*alphas[1];
				dr_u[2] += (lr[2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_du[y] * u_u*alphas[2];
			}
		}
	}
	else
	{
		int frm = f*nf_g_q_p;
	
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
	
		
		double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
	
		dr_u[0] = 0;
		dr_u[1] = 0;
		dr_u[2] = 0;
		int ncoords_base = theFormulas[frm].ncoor;
		for( int x = 0; x < ncoords_base; x++ )
		{
			double *lr = r + cp[x]*3;
			
			dr_u[0] += ceff_map_du[x] * (lr[0] + theFormulas[frm].r_pbc[3*x+0])*alphas[0];
			dr_u[1] += ceff_map_du[x] * (lr[1] + theFormulas[frm].r_pbc[3*x+1])*alphas[1];
			dr_u[2] += ceff_map_du[x] * (lr[2] + theFormulas[frm].r_pbc[3*x+2])*alphas[2];
		}	
	}	
}

void surface::rv( int f, double u, double v, double *r, double *dr_v )
{
	double *alphas = r+3*nv;
	dr_v[0] = 0;
	dr_v[1] = 0;
	dr_v[2] = 0;
	double alpha_x = r[3*nv+0];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];

	if( f >= nf_faces ) // an irrational face.
	{
		int frm = (f-nf_faces)*nf_irr_pts;

		int i = theIrregularFormulas[frm].vertex;
		int e = theIrregularFormulas[frm].edge;
		int val = theVertices[i].valence;

		irr_kernel *theKernel = kernels[val]; 
		int ncoords_base = val + 6;
		double u_u = u, u_v=0, v_u=0, v_v = v;
		
		theKernel->get_map_transform( &u_u, &u_v, &v_u, &v_v );

		double fu = u;
		double fv = v;					

		int domain = theKernel->domain(fu,fv);
		if( domain-1 >= theKernel->ndomains )
		{
			printf("ERROR trying to evaluate ru outside of irregular domain.\n");
			exit(1);
		}

		double *theMap = theKernel->get_map( &fu, &fv );

		double u = fu;
		double v = fv;
		double w = 1 - u - v;

		if( u +v > 1.0+THRESH || u < -THRESH || v < -THRESH )
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
		
		double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
	
		dr_v[0] = 0;	
		dr_v[1] = 0;	
		dr_v[2] = 0;	
		for( int x = 0; x < ncoords_base; x++ )
		{
			int *cset = theVertices[i].irr_coord_set + e * ncoords_base;
			double *lr = r + cset[x]*3;

			for( int y = 0; y < 12; y++ )
			{
				dr_v[0] += (lr[0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_dv[y] * v_v*alphas[0];
				dr_v[1] += (lr[1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_dv[y] * v_v*alphas[1];
				dr_v[2] += (lr[2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_dv[y] * v_v*alphas[2];
			}
		}
	}
	else
	{
		int frm = f*nf_g_q_p;
	
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
		
		double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
	
		dr_v[0] = 0;
		dr_v[1] = 0;
		dr_v[2] = 0;
		int ncoords_base = theFormulas[frm].ncoor;
		for( int x = 0; x < ncoords_base; x++ )
		{
			double *lr = r + cp[x]*3;
			
			dr_v[0] += ceff_map_dv[x] * (lr[0] + theFormulas[frm].r_pbc[3*x+0])*alphas[0];
			dr_v[1] += ceff_map_dv[x] * (lr[1] + theFormulas[frm].r_pbc[3*x+1])*alphas[1];
			dr_v[2] += ceff_map_dv[x] * (lr[2] + theFormulas[frm].r_pbc[3*x+2])*alphas[2];
		}	
	}	

}

int surface::nextFace( int f, double *u_in, double *v_in, double *du_in, double *dv_in, double *r, double *mom, double *coord_transform )
{
	if( *du_in == 0 && *dv_in == 0 )
		return f;

	int t_cur;

	if( f < nf_faces )
		t_cur = theFormulas[f*nf_g_q_p+0].tri;
	else
		t_cur = theIrregularFormulas[(f-nf_faces)*nf_irr_pts+0].tri;

	// I am going to return the point on the "next face" where we end up, and the leftover du, dv.

	double u = *u_in;
	double v = *v_in;
	double du = *du_in;
	double dv = *dv_in;

	if( u + v > 1 + THRESH || u < -THRESH || v < -THRESH )
	{
		printf("ERROR: pt not on face.\n");
		exit(1);
	}

	// which face do we leave out of ?

	double nu = u + du;
	double nv = v + dv;

	if( nu > 0 && nv > 0 && nu+nv < 1 )
	{
		*u_in += *du_in;
		*v_in += *dv_in;
		*du_in = 0;
		*dv_in = 0;

		return f;

	}

	// compute t parameter for leaving each face.

	double t_1 = 1e10;
	double t_2 = 1e10;
	double t_3 = 1e10;

	if( dv < 0 )
		t_1 = -v / dv;
	if( du+dv > 0 )
		t_2 = (1-u-v)/(du+dv);
	if( du < 0 )
		t_3 = -u / du;

	int face;
			
	double ru_old[3];
	double rv_old[3];

	double rold[3];
	double rnew[3];
	double nrm[3];

	int edge_type, sub_type=0;

	if( t_1 < t_2 && t_1 < t_3 )
	{
		edge_type = 0;
		u += t_1 * du;
		v += t_1 * dv;
		du -= t_1 * du;
		dv -= t_1 * dv;

#ifdef MOVE_CHECK		
		evaluateRNRM( f, u, v, rold, nrm, r );
#endif
		ru(f, u, v, r, ru_old );
		rv(f, u, v, r, rv_old ); 			

		// leaves by border zero

		face = theTriangles[theTriangles[t_cur].border_tri[0]].f;

		if( theTriangles[t_cur].edge_type & MASK_1 )
		{
			sub_type = 1;
			*u_in = u;
			*v_in = 1- *u_in;
	
		}	
		else
		{
			*u_in = 0;
			*v_in = u;
		}
		
	}
	else if( t_2 < t_3 )
	{
		edge_type = 1;
		u += t_2 * du;
		v += t_2 * dv;
		du -= t_2 * du;
		dv -= t_2 * dv;

#ifdef MOVE_CHECK		
		evaluateRNRM( f, u, v, rold, nrm, r );
#endif
		ru(f, u, v, r, ru_old );
		rv(f, u, v, r, rv_old ); 			
		// leaves by border one
		
		face = theTriangles[theTriangles[t_cur].border_tri[1]].f;
		
		if( !(theTriangles[t_cur].edge_type & MASK_2) ) // j is lowest.
		{
			*u_in = 0;
			*v_in = v;

		}	
		else if ( theTriangles[t_cur].edge_type & (1<<1) ) // k is lowest
		{
			sub_type = 1;
			*u_in = u;
			*v_in = 0;
		}
		else // l is lowest. 
		{
			sub_type = 2;
			*v_in = u;
			*u_in = v;
		}
	}
	else
	{	
		edge_type = 2;
		u += t_3 * du;
		v += t_3 * dv;
		du -= t_3 * du;
		dv -= t_3 * dv;
		
#ifdef MOVE_CHECK		
		evaluateRNRM( f, u, v, rold, nrm, r );
#endif
		ru(f, u, v, r, ru_old );
		rv(f, u, v, r, rv_old ); 			

		// leaves by border two
		face = theTriangles[theTriangles[t_cur].border_tri[2]].f;
		
		if ( theTriangles[t_cur].edge_type & (MASK_3) ) // n is lowest
		{
			sub_type = 1;
			*v_in = v;
			*u_in = 1-*v_in;	
		}
		else
		{
			*u_in = v;
			*v_in = 0;
		}

	}

	double ru_new[3];
	double rv_new[3];

#ifdef MOVE_CHECK		
	evaluateRNRM( face, *u_in, *v_in, rnew, nrm, r );

	double dr[3] = { rnew[0] - rold[0], rnew[1] - rold[1], rnew[2] - rold[2] };
	double put[3];
	MinImage3D( dr, PBC_vec, put, r+3*this->nv );

	double dr_check = sqrt( dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2] );

	if( dr_check > 1e-10 )
	{
		printf("BORDER MOVE ERROR %le. edgetype: %d subtype: %d\n", dr_check, edge_type, sub_type );
		exit(1);
	}
#endif
	ru(face, *u_in, *v_in, r, ru_new );
	rv(face, *u_in, *v_in, r, rv_new ); 			
	
	double ruo1 = ru_old[0] * ru_old[0] + ru_old[1] * ru_old[1] + ru_old[2]*ru_old[2];
	double run1 = ru_new[0] * ru_old[0] + ru_new[1] * ru_old[1] + ru_new[2]*ru_old[2];
	double ruo2 = ru_old[0] * rv_old[0] + ru_old[1] * rv_old[1] + ru_old[2]*rv_old[2];
	double run2 = ru_new[0] * rv_old[0] + ru_new[1] * rv_old[1] + ru_new[2]*rv_old[2];
	
	double rvo1 = rv_old[0] * ru_old[0] + rv_old[1] * ru_old[1] + rv_old[2]*ru_old[2];
	double rvn1 = rv_new[0] * ru_old[0] + rv_new[1] * ru_old[1] + rv_new[2]*ru_old[2];
	double rvo2 = rv_old[0] * rv_old[0] + rv_old[1] * rv_old[1] + rv_old[2]*rv_old[2];
	double rvn2 = rv_new[0] * rv_old[0] + rv_new[1] * rv_old[1] + rv_new[2]*rv_old[2];

//	(du r_u + dv r_v == du' r_u' + dv' r_v') . ru 
//	(du r_u + dv r_v == du' r_u' + dv' r_v') . rv 

		

	double dup = -(-(du*ruo2*rvn1) + du*ruo1*rvn2 + dv*rvn2*rvo1 - dv*rvn1*rvo2)/(run2*rvn1 - run1*rvn2);
	double dvp = -(-(du*run2*ruo1) + du*run1*ruo2 - dv*run2*rvo1 + dv*run1*rvo2)/(run2*rvn1 - run1*rvn2);

/*
	printf("tangent: %lf %lf %lf dup: %lf dvp: %lf du: %lf dv: %lf\n", 
		du * ru_old[0] + dv * rv_old[0],
		du * ru_old[1] + dv * rv_old[1],
		du * ru_old[2] + dv * rv_old[2],
			dup, dvp, du, dv );
*/
	if( mom )
	{
		double rotor[4] = 
		{
			(ruo2*rvn1 - ruo1*rvn2) / (run2*rvn1-run1*rvn2),  -(rvn2*rvo1-rvn1*rvo2) / (run2*rvn1-run1*rvn2),
			(run2*ruo1 - run1*ruo2) / (run2*rvn1-run1*rvn2),   (run2*rvo1-run1*rvo2) / (run2*rvn1-run1*rvn2)
		};

		if( coord_transform )
		{
			double t[4];
	
			//t[0] = coord_transform[0] * rotor[0] + coord_transform[1] * rotor[2];
			//t[1] = coord_transform[0] * rotor[1] + coord_transform[1] * rotor[3];
			//t[2] = coord_transform[2] * rotor[0] + coord_transform[3] * rotor[2];
			//t[3] = coord_transform[2] * rotor[1] + coord_transform[3] * rotor[3];
			t[0] = rotor[0]*coord_transform[0] +  rotor[1]*coord_transform[2];
			t[1] = rotor[0]*coord_transform[1] +  rotor[1]*coord_transform[3];
			t[2] = rotor[2]*coord_transform[0] +  rotor[3]*coord_transform[2];
			t[3] = rotor[2]*coord_transform[1] +  rotor[3]*coord_transform[3];
			memcpy( coord_transform, t, sizeof(double) * 4 );
		}

		double mu = mom[0] * rotor[0] + mom[1] * rotor[1];
		double mv = mom[0] * rotor[2] + mom[1] * rotor[3];
	

//		printf("du_old %le dv_old %le du_new %le dv_new %le mu_old %le mv_old %le mu_new %le mv_new %le\n",
//			du, dv, dup, dvp, mom[0], mom[1], mu, mv );
		mom[0] = mu;
		mom[1] = mv;
	}

	*du_in = dup;
	*dv_in = dvp;

	return face;
}

double surface::trialMove( int *f_in, double *u_in, double *v_in, double duv[2], double lambda, double *r, double *base_position )
{
	double u = *u_in;
	double v = *v_in;
	double du = duv[0] * lambda;
	double dv = duv[1] * lambda;

	int f_1 = *f_in, nf = *f_in;
	do {
		f_1 = nf;
		nf = nextFace( f_1, &u, &v, &du, &dv, r );
	} while( nf != f_1 );
	
	u += du;
	v += dv;
	*f_in = nf;
	*u_in = u;
	*v_in = v;

	double rval[3];
	double nval[3];

	evaluateRNRM( *f_in, *u_in, *v_in, rval, nval, r);
 	
	double dr[3] = { rval[0] - base_position[0],
			 rval[1] - base_position[1], 
			 rval[2] - base_position[2] };
	
	double l = sqrt( dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

	return l;
}


void surface::localMove( int *f_in, double *u_in, double *v_in, double sigma, double *r, double *frc_duv, double dt, double *fstep, int max_correct_iterations)
{
	int f = *f_in;
	double u = *u_in;
	double v = *v_in;
	double drdu[3];
	double drdv[3];

	ru(f, u, v, r, drdu );
	rv(f, u, v, r, drdv );

	double lu = sqrt(drdu[0]*drdu[0]+drdu[1]*drdu[1]+drdu[2]*drdu[2]);
	double lv = sqrt(drdv[0]*drdv[0]+drdv[1]*drdv[1]+drdv[2]*drdv[2]);

	double U[3] = { drdu[0],drdu[1],drdu[2]};
	double V[3] = {0,0,0};
	double proj;

	proj = (drdv[0] * U[0] + drdv[1] * U[1] + drdv[2] * U[2])/(U[0] * U[0] + U[1] * U[1] + U[2] * U[2]);

	V[0] += drdv[0] - proj * U[0];
	V[1] += drdv[1] - proj * U[1];
	V[2] += drdv[2] - proj * U[2];

	double nrm_U = normalize(U);
	double nrm_V = normalize(V);
	
	if( ! rng_x )
		init_random(0);

	double theta = 2 * M_PI * rand() / (double)RAND_MAX;
	double l = gsl_ran_gaussian(rng_x, sqrt(2.0) * sigma);

	double cvec[3] = { U[0] * cos(theta) + V[0] * sin(theta),
			   U[1] * cos(theta) + V[1] * sin(theta),
			   U[2] * cos(theta) + V[2] * sin(theta) };
	cvec[0] *= l;
	cvec[1] *= l;
	cvec[2] *= l;

	double cvu = cvec[0] * drdu[0] + cvec[1] * drdu[1] + cvec[2] * drdu[2];
	double cvv = cvec[0] * drdv[0] + cvec[1] * drdv[1] + cvec[2] * drdv[2];
	
	double coeff_cross = (cvv - cvu * proj)/nrm_V/nrm_V;
	
	double du = cvu/lu/lu - coeff_cross * proj;
	double dv = coeff_cross;  

//	double du = (cvec[0] * drdu[0] + cvec[1] * drdu[1] + cvec[2] * drdu[2])/lu/lu;
//	double dv = (drdv[0] * cvec[0] + drdv[1] * cvec[1] + drdv[2] * cvec[2])/nrm_V/nrm_V;
//	dv -= 

	

/*
	// dudr * lam == vec
	double du = (cvec[0] * drdu[0] + cvec[1] * drdu[1] + cvec[2] * drdu[2])/lu;
	cvec[0] -= du * drdu[0]/lu;
	cvec[1] -= du * drdu[1]/lu;
	cvec[2] -= du * drdu[2]/lu;
	double dv = (cvec[0] * drdv[0] + cvec[1] * drdv[1] + cvec[2] * drdv[2])/lv;

	cvec[0] -= dv * drdv[0]/lv;
	cvec[1] -= dv * drdv[1]/lv;
	cvec[2] -= dv * drdv[2]/lv;
*/
	cvec[0] -= du * drdu[0] + dv * drdv[0];
	cvec[1] -= du * drdu[1] + dv * drdv[1];
	cvec[2] -= du * drdu[2] + dv * drdv[2];

//	printf("resid: %le\n", cvec[0]*cvec[0]+cvec[1]*cvec[1]+cvec[2]*cvec[2] );



	// x1 r1 + x2 r2 == z

	/*	double theta = 0;
	double dx,dy;
       	while( theta <= 2 * M_PI )
	{
		dx = (U[0] * cos(theta) + V[0] * sin(theta));
		dy = (U[1] * cos(theta) + V[1] * sin(theta));
		printf("%lf %lf\n", dx, dy);
		theta += M_PI/16;
		}*/

        int i;
	
//	double du, dv;


	// sigma r 

//	du = gsl_ran_gaussian(rng_x, sigma/nrm_U); // dr * / (dr/du) == du
//	dv = gsl_ran_gaussian(rng_x, sigma/nrm_V);
//	printf("%lf %lf %lf %lf\n", du, dv, sigma, nrm_U);

//        gsl_rng_free(x);
	du /= sqrt(dt);
	dv /= sqrt(dt);

#ifdef DBG_PRINT
	printf("drdu: %le %le %le\n", drdu[0], drdu[1], drdu[2] );
	printf("drdv: %le %le %le\n", drdv[0], drdv[1], drdv[2] );
	printf("U: %le %le %le\n", U[0], U[1], U[2] );
	printf("V: %le %le %le\n", V[0], V[1], V[2] );
#endif
	if( frc_duv )
	{
		// key: U and V are orthonormal.  find the inverse transformation to convert tangent space vectors to u,v:
		double ut1 = U[0] * drdu[0] + U[1] * drdu[1] + U[2] * drdu[2];
		double ut2 = V[0] * drdu[0] + V[1] * drdu[1] + V[2] * drdu[2]; 
		double vt1 = U[0] * drdv[0] + U[1] * drdv[1] + U[2] * drdv[2];
		double vt2 = V[0] * drdv[0] + V[1] * drdv[1] + V[2] * drdv[2];
			      
		double idet = 1.0 / (ut1 * vt2 - ut2 * vt1 );

		// invert this matrix.
		 
		double dudt1 = idet * vt2;
		//double dudt2 =-idet * ut2;
		//double dvdt1 =-idet * vt1;
		double dudt2 =-idet * vt1;
		double dvdt1 =-idet * ut2;
		double dvdt2 = idet * ut1;

		// the gradient in uv space.
		double dE_du = frc_duv[0];
		double dE_dv = frc_duv[1];

		
		// the gradient in the arbitrary tangent space.
		double G_t1 = dE_du * dudt1 + dE_dv * dvdt1;
		double G_t2 = dE_du * dudt2 + dE_dv * dvdt2;


		double fac = 1e8 / 0.592;

		// given this, does it take the right step size???
	
		// the step in uv space
			
		double step_u = G_t1 * dudt1 + G_t2 * dudt2;
		double step_v = G_t1 * dvdt1 + G_t2 * dvdt2; 

		if( dt*sqrt(frc_duv[0]*frc_duv[0]+frc_duv[1]*frc_duv[1]) > 1e5 )
		{
			printf("f %d u %le v %le real step: %le %le uv step %lf %lf l1 l2 %lf %lf\n", f, u, v, dt*G_t1, dt*G_t2, dt*step_u, dt*step_v, nrm_U, nrm_V ); 
			printf("duv before: %le %le after %le %le\n", du, dv, du+step_u, dv+step_v );
		}

		if( fstep )
		{
			int tf = f;
			double tu = u;
			double tv = v;
			double tdu = step_u*dt;
			double tdv = step_v*dt;
	
			double r_in[3];
			double njunk[3];
	
			evaluateRNRM( f, u, v, r_in, njunk, r );
			{
				int nf = tf;
				do {
					tf = nf;
					nf = nextFace( tf, &tu, &tv, &tdu, &tdv, r );
				} while( nf != tf );
			}
			
			double r_out[3];
			evaluateRNRM( tf, tu, tv, r_out, njunk, r );
	
			fstep[0] = r_out[0] - r_in[0];
			fstep[1] = r_out[1] - r_in[1];
			fstep[2] = r_out[2] - r_in[2];
		}
		du += step_u;
		dv += step_v;
	}

	du *= dt;
	dv *= dt;
	
	if( max_correct_iterations > 0 && f >= nf_faces ) // only on irregular faces.
	{
		// what is the cartesian length we expect to take? 
	
		double expec[3] = { du * drdu[0] + dv * drdv[0],
				    du * drdu[1] + dv * drdv[1],	
				    du * drdu[2] + dv * drdv[2] };	
		double basep[3], basen[3];
		evaluateRNRM( f,u,v,basep,basen,r);
		double le = sqrt( expec[0]*expec[0]+expec[1]*expec[1]+expec[2]*expec[2]);
			
		double lambda_low = 0;
		double lambda_mid = 1.0;
		double lambda_high = 2.0;		

		if(  f >= nf_faces )
		{
			lambda_mid = 0.5;
			lambda_high = 1.0;
		}

		double rlow = 0, rmid, rhigh;
		int tf = f;
		double tu = u;
		double tv = v;	
		double duv[2] = { du, dv };
		tf = f; tu = u; tv = v; rmid  = trialMove( &tf, &tu, &tv, duv, lambda_mid, r, basep );

		int bestf = tf;
		double bestu = tu;
		double bestv = tv; 

		tf = f; tu = u; tv = v; rhigh = trialMove( &tf, &tu, &tv, duv, lambda_high, r, basep );
//		printf("Start: rlow: %le rmid: %le rhigh: %le\n", rlow, rmid, rhigh );
		while( rhigh < le  )
		{
			lambda_high *= 1.25;
			tf = f; tu = u; tv = v; rhigh = trialMove( &tf, &tu, &tv, duv, lambda_high, r, basep );
//			printf("iterating lambda_high up to %le got %le trying for %le\n", lambda_high, rhigh, le );
		}
	
		for( int ir = 0; ir < max_correct_iterations; ir++ )
		{
			double rtest;	
			double ltest = (lambda_mid+lambda_high)/2;
			tf = f; tu = u; tv = v; rtest  = trialMove( &tf, &tu, &tv, duv, ltest, r, basep );
			
			if( rtest < le )
			{
				rlow = rmid;
				lambda_low = lambda_mid;
				rmid = rtest;
				lambda_mid = ltest;
				bestf = tf;
				bestu = tu;
				bestv = tv;
//				printf("testing mid-high and we are on the left, moving bracket up to %le %le %le\n", lambda_low, lambda_mid, lambda_high ); 
			}
			else
			{	
				rhigh = rtest;
				lambda_high = ltest;
//				printf("testing mid-high and we are on the right, moving bracket down a bit to %le %le %le\n", lambda_low, lambda_mid, lambda_high ); 
			}
			
			ltest = (lambda_mid+lambda_low)/2;
			tf = f; tu = u; tv = v; rtest  = trialMove( &tf, &tu, &tv, duv, ltest, r, basep );
			
			if( rtest > le )
			{
				rhigh = rmid;
				lambda_high= lambda_mid;
				rmid = rtest;
				lambda_mid = ltest;
				bestf = tf;
				bestu = tu;
				bestv = tv;
//				printf("testing mid-low and we are on the right, moving bracket down to to %le %le %le\n", lambda_low, lambda_mid, lambda_high ); 
			}
			else
			{	
				rlow = rtest;
				lambda_low = ltest;
//				printf("testing mid-low and we are on the left, moving bracket up a bit to %le %le %le\n", lambda_low, lambda_mid, lambda_high ); 
			}

			if( fabs(le-rmid) < 1e-6 )
				break;
		}

//#define DEBUG_IRR
#ifdef DEBUG_IRR
		if( f >= nf_faces )
			printf(" irr lambda_scale %le Found step %le trying for %le del %le\n", lambda_mid, rmid, le, rmid-le ); 
		else
			printf(" reg lambda_scale %le Found step %le trying for %le del %le\n", lambda_mid, rmid, le, rmid-le ); 
#endif
		*f_in = bestf;
		*u_in = bestu;
		*v_in = bestv; 
	}
	else
	{
		int f_1 = f, nf = f;
		do {
			f_1 = nf;
			nf = nextFace( f_1, &u, &v, &du, &dv, r );
		} while( nf != f_1 );
	
		u += du;
		v += dv;
		*f_in = nf;
		*u_in = u;
		*v_in = v;
	}
}

int surface::neighborList( int f, int *neighbors )
{
	return 0;
}

int surface::map( int face_from, int face_to, double u, double v )
{
	return 0;
}


double surface::c( int f, double u, double v, double *r )
{
	double alpha_x = r[3*nv+0];
	double alpha_y = r[3*nv+1];
	double alpha_z = r[3*nv+2];

	int do_monge = 0;

	if( f >= nf_faces ) // an irrational face.
	{
		int frm = (f-nf_faces)*nf_irr_pts;

		int i = theIrregularFormulas[frm].vertex;
		int e = theIrregularFormulas[frm].edge;
		int val = theVertices[i].valence;

		irr_kernel *theKernel = kernels[val]; 
		int ncoords_base = val + 6;

		double fu = u;
		double fv = v;					

		int domain = theKernel->domain(fu,fv);
		if( domain-1 >= theKernel->ndomains )
		{
			printf("ERROR trying to evaluate ru outside of irregular domain.\n");
			exit(1);
		}

		double *theMap = theKernel->get_map( &fu, &fv );

		double u = fu;
		double v = fv;
		double w = 1 - u - v;

		if( u +v > 1.0+THRESH || u < -THRESH || v < -THRESH )
		{
			printf("u: %lf v: %lf\n", u, v );
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
		
		double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
		double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
		
		double ceff_map_duu[12] = { d_uu_8, d_uu_7, d_uu_4, d_uu_5, d_uu_9, d_uu_12, d_uu_11, d_uu_10, d_uu_6, d_uu_3, d_uu_1, d_uu_2 };
		double ceff_map_duv[12] = { d_uv_8, d_uv_7, d_uv_4, d_uv_5, d_uv_9, d_uv_12, d_uv_11, d_uv_10, d_uv_6, d_uv_3, d_uv_1, d_uv_2 };
		double ceff_map_dvv[12] = { d_vv_8, d_vv_7, d_vv_4, d_vv_5, d_vv_9, d_vv_12, d_vv_11, d_vv_10, d_vv_6, d_vv_3, d_vv_1, d_vv_2 };
		
		double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
		double nrm[3]={0,0,0}; 
	
		int *cset = theVertices[i].irr_coord_set + e * ncoords_base;
	
		for( int x = 0; x < ncoords_base; x++ )
		{
			for( int y = 0; y < 12; y++ )
			{
				Ru[0] += (r[3*cset[x]+0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_du[y] * alpha_x;
				Ru[1] += (r[3*cset[x]+1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_du[y] * alpha_y;
				Ru[2] += (r[3*cset[x]+2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_du[y] * alpha_z;
				
				Rv[0] += (r[3*cset[x]+0] + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_dv[y] * alpha_x;
				Rv[1] += (r[3*cset[x]+1] + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_dv[y] * alpha_y;
				Rv[2] += (r[3*cset[x]+2] + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_dv[y] * alpha_z;
				
				tSuu[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_duu[y] * alpha_x;
				tSuu[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_duu[y] * alpha_y;
				tSuu[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_duu[y] * alpha_z;
				
				tSuv[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_duv[y] * alpha_x;
				tSuv[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_duv[y] * alpha_y;
				tSuv[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_duv[y] * alpha_z;
				
				tSvv[0] += (r[3*cset[x]+0]  + theIrregularFormulas[frm].r_pbc[3*x+0]) * theMap[y*ncoords_base+x] * ceff_map_dvv[y] * alpha_x;
				tSvv[1] += (r[3*cset[x]+1]  + theIrregularFormulas[frm].r_pbc[3*x+1]) * theMap[y*ncoords_base+x] * ceff_map_dvv[y] * alpha_y;
				tSvv[2] += (r[3*cset[x]+2]  + theIrregularFormulas[frm].r_pbc[3*x+2]) * theMap[y*ncoords_base+x] * ceff_map_dvv[y] * alpha_z;
			}
		}
		cross( Ru, Rv, nrm );
		normalize(nrm);
	
		/*
			irregular vertex domain scaling (domain^2) is absent here. it cancels.
		*/

	
		double RuRu = Ru[0] * Ru[0] + Ru[1] * Ru[1] + Ru[2]*Ru[2];
		double RuRv = Ru[0] * Rv[0] + Ru[1] * Rv[1] + Ru[2]*Rv[2];
		double RvRv = Rv[0] * Rv[0] + Rv[1] * Rv[1] + Rv[2]*Rv[2];
	
		double g = sqrt(RuRu*RvRv-RuRv*RuRv);
	
		double c1,c2;

		double csum = 0;
	
		if( do_monge )
		{
		      double idet = 1.0 / (Ru[0] * Rv[1] - Ru[1] * Rv[0]);
	
			double dudx=idet*Rv[1];
			double dudy=-idet*Rv[0];
			double dvdx=-idet*Ru[1];
			double dvdy=idet*Ru[0];
	
			double dhdx=Ru[2]*dudx+Rv[2]*dvdx;
			double dhdy=Ru[2]*dudy+Rv[2]*dvdy;
	
			c1 = tSuu[2] * dudx * dudx + tSuv[2] * dudx * dvdx + tSuv[2] * dvdx * dudx + tSvv[2] * dvdx * dvdx;
			c2 = tSuu[2] * dudy * dudy + tSuv[2] * dudy * dvdy + tSuv[2] * dvdy * dudy + tSvv[2] * dvdy * dvdy; 

			csum = c1+c2;
		}
		else
		{
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
			c1 = -0.5*(a+d-sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
			c2 = -0.5*(a+d+sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
			
			csum = c1+c2;
		}


		return csum;
	}
	else
	{
		int frm = f*nf_g_q_p;
	
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
		
		double ceff_map_du[12] = { du_8, du_7, du_4, du_5, du_9, du_12, du_11, du_10, du_6, du_3, du_1, du_2 };
		double ceff_map_dv[12] = { dv_8, dv_7, dv_4, dv_5, dv_9, dv_12, dv_11, dv_10, dv_6, dv_3, dv_1, dv_2 };
		
		double ceff_map_duu[12] = { d_uu_8, d_uu_7, d_uu_4, d_uu_5, d_uu_9, d_uu_12, d_uu_11, d_uu_10, d_uu_6, d_uu_3, d_uu_1, d_uu_2 };
		double ceff_map_duv[12] = { d_uv_8, d_uv_7, d_uv_4, d_uv_5, d_uv_9, d_uv_12, d_uv_11, d_uv_10, d_uv_6, d_uv_3, d_uv_1, d_uv_2 };
		double ceff_map_dvv[12] = { d_vv_8, d_vv_7, d_vv_4, d_vv_5, d_vv_9, d_vv_12, d_vv_11, d_vv_10, d_vv_6, d_vv_3, d_vv_1, d_vv_2 };
	
		int ncoords_base = theFormulas[frm].ncoor;
			
		double R[3] = {0,0,0}, Ru[3]={0,0,0}, Rv[3]={0,0,0}, tSuu[3]={0,0,0}, tSuv[3]={0,0,0}, tSvv[3]={0,0,0};
		double nrm[3]={0,0,0}; 
		for( int p = 0; p < ncoords_base; p++ )
		{
			double *lr = r + cp[p]*3;
				
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

		double c1,c2;
      
		if( do_monge )
		{
		      double idet = 1.0 / (Ru[0] * Rv[1] - Ru[1] * Rv[0]);

			double dudx=idet*Rv[1];
			double dudy=-idet*Rv[0];
			double dvdx=-idet*Ru[1];
			double dvdy=idet*Ru[0];

			double dhdx=Ru[2]*dudx+Rv[2]*dvdx;
			double dhdy=Ru[2]*dudy+Rv[2]*dvdy;
      
			c1 = tSuu[2] * dudx * dudx + tSuv[2] * dudx * dvdx + tSuv[2] * dvdx * dudx + tSvv[2] * dvdx * dvdx;
			c2 = tSuu[2] * dudy * dudy + tSuv[2] * dudy * dvdy + tSuv[2] * dvdy * dudy + tSvv[2] * dvdy * dvdy; 
		}
		else
		{
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
			c1 = -0.5*(a+d-sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
			c2 = -0.5*(a+d+sqrt(SAFETY+a*a+4*b*c-2*a*d+d*d));
		}
		return c1 + c2;
	}
	
	return 0;
}


void surface::fetchPuv( int f, double u, double v, double *P_uv, double *rsurf )
{
	double drdu[3];
	double drdv[3];
	
	ru(f, u, v, rsurf, drdu );
	rv(f, u, v, rsurf, drdv );

	P_uv[0] = drdu[0] * drdu[0] + drdu[1] * drdu[1] + drdu[2] * drdu[2];
	P_uv[1] = drdu[0] * drdv[0] + drdu[1] * drdv[1] + drdu[2] * drdv[2];
	P_uv[2] = P_uv[1];
	P_uv[3] = drdv[0] * drdv[0] + drdv[1] * drdv[1] + drdv[2] * drdv[2];
}

int surface::fetchCiu( int f, double u, double v, double **C_iu, double *rsurf )
{
	double *alphas = rsurf + 3 * nv;
	double coefs[MAX_VALENCE+6];
	int coef_list[MAX_VALENCE+6];
	int ncoef=0;
	get_pt_coeffs( f, u, v, coefs, coef_list, &ncoef ); 

	// this is the coupling matrix between the mesh and the particle

	(*C_iu) = (double *)malloc(sizeof(double) * (ncoef*3*2) ); 

	double drdu[3];
	double drdv[3];
	
	ru(f, u, v, rsurf, drdu );
	rv(f, u, v, rsurf, drdv );

	double *mat = *C_iu;

	for( int c = 0; c < ncoef; c++ )
	{	// dr_dqx * drx_du 
		mat[c*2*3+0*3+0] = alphas[0] * coefs[c] * drdu[0];  
		mat[c*2*3+0*3+1] = alphas[1] * coefs[c] * drdu[1];  
		mat[c*2*3+0*3+2] = alphas[2] * coefs[c] * drdu[2];  
		
		mat[c*2*3+1*3+0] = alphas[0] * coefs[c] * drdv[0];  
		mat[c*2*3+1*3+1] = alphas[1] * coefs[c] * drdv[1];  
		mat[c*2*3+1*3+2] = alphas[2] * coefs[c] * drdv[2];  
	}

	return ncoef;
}



int surface::fetchdP_duv( int f, double u, double v, double **dP_duv, double *rsurf )
{
	double *alphas = rsurf + 3 * nv;
	double coefs[3*(MAX_VALENCE+6)];
	int coef_list[MAX_VALENCE+6];
	int ncoef=0;
	get_pt_dcoeffs( f, u, v, coefs, coef_list, &ncoef ); 

	// this is the derivative of the PUV matrix wrt u and v, as well as wrt each coordinate.

	(*dP_duv) = (double *)malloc(sizeof(double) * (4*2+ncoef*3*4) ); 

	double drdu[3];
	double drdv[3];
	
	double druu[3]={0,0,0};
	double druv[3]={0,0,0};
	double drvv[3]={0,0,0};
	
	ru(f, u, v, rsurf, drdu );
	rv(f, u, v, rsurf, drdv );

	r2der( f, u, v, rsurf, druu, druv, drvv );

	double *mat = *dP_duv;

	// derivative of the matrix wrt u
	mat[0] = 2*(druu[0] * drdu[0] + druu[1] * drdu[1] + druu[2] * drdu[2]);
	mat[1] =    druu[0] * drdv[0] + druu[1] * drdv[1] + druu[2] * drdv[2] + 
		    druv[0] * drdu[0] + druv[1] * drdu[1] + druv[2] * drdu[2];
	mat[2] = mat[1];
	mat[3] = 2*(druv[0] * drdv[0] + druv[1] * drdv[1] + druv[2] * drdv[2]);
	
	// derivative of the matrix wrt v
	mat[4] = 2*(druv[0] * drdu[0] + druv[1] * drdu[1] + druv[2] * drdu[2]);
	mat[5] =    druv[0] * drdv[0] + druv[1] * drdv[1] + druv[2] * drdv[2] + 
		    drvv[0] * drdu[0] + drvv[1] * drdu[1] + drvv[2] * drdu[2];
	mat[6] = mat[5];
	mat[7] = 2*(drvv[0] * drdv[0] + drvv[1] * drdv[1] + drvv[2] * drdv[2]);

	// derivative of the matrix wrt coordinates.

	for( int c = 0; c < ncoef; c++ )
	{	// dr_dqx * drx_du 

		// this is the derivative of the mesh coupling matrix wrt u. 

		// iu,du
		mat[8+c*4*3+0*3+0] = alphas[0] * coefs[c] * druu[0] + alphas[0] * coefs[ncoef+c] * drdu[0];  
		mat[8+c*4*3+0*3+1] = alphas[1] * coefs[c] * druu[1] + alphas[1] * coefs[ncoef+c] * drdu[1];  
		mat[8+c*4*3+0*3+2] = alphas[2] * coefs[c] * druu[2] + alphas[2] * coefs[ncoef+c] * drdu[2];  
		
		// iu,dv
		mat[8+c*4*3+1*3+0] = alphas[0] * coefs[c] * druv[0] + alphas[0] * coefs[2*ncoef+c] * drdu[0];  
		mat[8+c*4*3+1*3+1] = alphas[1] * coefs[c] * druv[1] + alphas[1] * coefs[2*ncoef+c] * drdu[1];  
		mat[8+c*4*3+1*3+2] = alphas[2] * coefs[c] * druv[2] + alphas[2] * coefs[2*ncoef+c] * drdu[2];  
		
		// iv,du
		mat[8+c*4*3+2*3+0] = alphas[0] * coefs[c] * druv[0] + alphas[0] * coefs[ncoef+c] * drdv[0];  
		mat[8+c*4*3+2*3+1] = alphas[1] * coefs[c] * druv[1] + alphas[1] * coefs[ncoef+c] * drdv[1];  
		mat[8+c*4*3+2*3+2] = alphas[2] * coefs[c] * druv[2] + alphas[2] * coefs[ncoef+c] * drdv[2];  
		
		// iv,dv
		mat[8+c*4*3+3*3+0] = alphas[0] * coefs[c] * drvv[0] + alphas[0] * coefs[2*ncoef+c] * drdv[0];  
		mat[8+c*4*3+3*3+1] = alphas[1] * coefs[c] * drvv[1] + alphas[1] * coefs[2*ncoef+c] * drdv[1];  
		mat[8+c*4*3+3*3+2] = alphas[2] * coefs[c] * drvv[2] + alphas[2] * coefs[2*ncoef+c] * drdv[2];  
	}
		

	return ncoef;
}




