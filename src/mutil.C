#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "mutil.h"
#include <math.h>
#include "lapack_we_use.h"
double Power( double a, double b) { return pow(a,b); }
double Sqrt( double a) { return sqrt(a); }
double dot( double *a, double *b )
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double cross( double *dr1, double *dr2, double *cp )
{
	cp[0] = (dr1[1] * dr2[2] - dr1[2] * dr2[1]);
	cp[1] =-(dr1[0] * dr2[2] - dr1[2] * dr2[0]);
	cp[2] = (dr1[0] * dr2[1] - dr1[1] * dr2[0]);

	double l = sqrt( cp[0]*cp[0]+cp[1]*cp[1]+cp[2]*cp[2]);

	return l;
}


double length3(double*dr)
{
	double lr = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
	return lr;
}
double normalize( double *dr )
{
	double lr = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

	if( fabs(lr) > 0  )
	{
		dr[0] /= lr;
		dr[1] /= lr;
		dr[2] /= lr;
	}

	return lr;
}

int line_segment_triangle_intersection( double *r1, double *r2, double *v1, double *v2, double *v3, double fudge )
{
	double plane_normal[3];

	double dr1[3] = { v1[0]-v2[0],v1[1]-v2[1],v1[2]-v2[2] };
	double dr2[3] = { v3[0]-v2[0],v3[1]-v2[1],v3[2]-v2[2] };
	double dr3[3] = { v3[0]-v1[0],v3[1]-v1[1],v3[2]-v1[2] };

	cross(dr1,dr2,plane_normal);
	
	double k = plane_normal[0] * v2[0] + plane_normal[1] * v2[1] + plane_normal[2] * v2[2];

	// solve	r1 + t (r2-r1) . plane_normal == k

	double dr[3] = {r2[0]-r1[0], r2[1]-r1[1],r2[2]-r1[2]};
	
	double dp_dr = dr[0] * plane_normal[0] + dr[1] * plane_normal[1] + dr[2] * plane_normal[2];
	double dp_1 = r1[0] * plane_normal[0] + r1[1] * plane_normal[1] + r1[2] * plane_normal[2];

	double t = (k - dp_1)/(dp_dr);

	if( t < -fudge || t > 1+fudge ) 
		return 0;

	double p_proj[3] = {
			r1[0] + t * dr[0],
			r1[1] + t * dr[1],
			r1[2] + t * dr[2] };

	double cross_vec1[3];
	double cross_vec2[3];
	double cross_vec3[3];
	
	// v1 - v2
	cross( dr1, plane_normal, cross_vec1 ); 
	// v3 - v2
	cross( dr2, plane_normal, cross_vec2 ); 
	// v3 - v1
	cross( dr3, plane_normal, cross_vec3 ); 

	double vec1A[3] = { p_proj[0] - v2[0], p_proj[1] - v2[1], p_proj[2] - v2[2] };
	double dp_1A = vec1A[0] * cross_vec1[0] + vec1A[1] * cross_vec1[1] + vec1A[2] * cross_vec1[2];	 
	double dp_1B =  dr2[0] * cross_vec1[0] + dr2[1] * cross_vec1[1] + dr2[2] * cross_vec1[2];	 
	double vec2A[3] = { p_proj[0] - v2[0], p_proj[1] - v2[1], p_proj[2] - v2[2] };
	double dp_2A = vec2A[0] * cross_vec2[0] + vec2A[1] * cross_vec2[1] + vec2A[2] * cross_vec2[2];	 
	double dp_2B = dr1[0] * cross_vec2[0] + dr1[1] * cross_vec2[1] + dr1[2] * cross_vec2[2];	 
	double vec3A[3] = { p_proj[0] - v1[0], p_proj[1] - v1[1], p_proj[2] - v1[2] };
	double dp_3A = vec3A[0] * cross_vec3[0] + vec3A[1] * cross_vec3[1] + vec3A[2] * cross_vec3[2];	 
	double dp_3B = -dr1[0] * cross_vec3[0] - dr1[1] * cross_vec3[1] - dr1[2] * cross_vec3[2];	 

	if( dp_1A * dp_1B >= 0 && dp_2A * dp_2B >= 0 && dp_3A * dp_3B >= 0 )
	{
		return 1;
	}
	return 0;
}

double triangle_area( double *pt1, double *pt2, double *pt3 )
{
	double v1x = pt1[0];
	double v1y = pt1[1];
	double v1z = pt1[2];
	
	double v2x = pt2[0];
	double v2y = pt2[1];
	double v2z = pt2[2];
	
	double v3x = pt3[0];
	double v3y = pt3[1];
	double v3z = pt3[2];

	return 0.5 * Sqrt(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
    Power(-(v1z*v2x) + v1x*v2z + v1z*v3x - v2z*v3x - v1x*v3z + v2x*v3z,2) + 
    Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2));

}

void d_triangle_area( double *der, double *pt1, double *pt2, double *pt3 )
{
	double v1x = pt1[0];
	double v1y = pt1[1];
	double v1z = pt1[2];
	
	double v2x = pt2[0];
	double v2y = pt2[1];
	double v2z = pt2[2];
	
	double v3x = pt3[0];
	double v3y = pt3[1];
	double v3z = pt3[2];

	der[0] = 0.5*(2*(v2y - v3y)*(-(v2y*v3x) + v1y*(-v2x + v3x) + v1x*(v2y - v3y) + v2x*v3y) + 
     2*(v2z - v3z)*(-(v2z*v3x) + v1z*(-v2x + v3x) + v1x*(v2z - v3z) + v2x*v3z))/
   (2.*Sqrt(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2)));

	der[1] = 0.5*(2*(v2x - v3x)*(v1y*(v2x - v3x) + v2y*v3x - v2x*v3y + v1x*(-v2y + v3y)) + 
     2*(v2z - v3z)*(-(v2z*v3y) + v1z*(-v2y + v3y) + v1y*(v2z - v3z) + v2y*v3z))/
   (2.*Sqrt(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2)));

	der[2] = 0.5*(2*(v2x - v3x)*(v1z*(v2x - v3x) + v2z*v3x - v2x*v3z + v1x*(-v2z + v3z)) + 
     2*(v2y - v3y)*(v1z*(v2y - v3y) + v2z*v3y - v2y*v3z + v1y*(-v2z + v3z)))/
   (2.*Sqrt(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2)));

	der[3] = 0.5*(2*(v1y - v3y)*(v1y*(v2x - v3x) + v2y*v3x - v2x*v3y + v1x*(-v2y + v3y)) + 
     2*(v1z - v3z)*(v1z*(v2x - v3x) + v2z*v3x - v2x*v3z + v1x*(-v2z + v3z)))/
   (2.*Sqrt(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2)));

	der[4] = 0.5*(2*(v1x - v3x)*(-(v2y*v3x) + v1y*(-v2x + v3x) + v1x*(v2y - v3y) + v2x*v3y) + 
     2*(v1z - v3z)*(v1z*(v2y - v3y) + v2z*v3y - v2y*v3z + v1y*(-v2z + v3z)))/
   (2.*Sqrt(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2)));

	der[5] = 0.5*(2*(v1x - v3x)*(-(v2z*v3x) + v1z*(-v2x + v3x) + v1x*(v2z - v3z) + v2x*v3z) + 
     2*(v1y - v3y)*(-(v2z*v3y) + v1z*(-v2y + v3y) + v1y*(v2z - v3z) + v2y*v3z))/
   (2.*Sqrt(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2)));

	der[6] = 0.5*(2*(-v1y + v2y)*(v1y*(v2x - v3x) + v2y*v3x - v2x*v3y + v1x*(-v2y + v3y)) + 
     2*(v1z - v2z)*(-(v2z*v3x) + v1z*(-v2x + v3x) + v1x*(v2z - v3z) + v2x*v3z))/
   (2.*Sqrt(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2)));

	der[7] = 0.5*(2*(v1x - v2x)*(v1y*(v2x - v3x) + v2y*v3x - v2x*v3y + v1x*(-v2y + v3y)) + 
     2*(-v1z + v2z)*(v1z*(v2y - v3y) + v2z*v3y - v2y*v3z + v1y*(-v2z + v3z)))/
   (2.*Sqrt(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2)));

	der[8] = 0.5*(2*(-v1x + v2x)*(-(v2z*v3x) + v1z*(-v2x + v3x) + v1x*(v2z - v3z) + v2x*v3z) + 
     2*(v1y - v2y)*(v1z*(v2y - v3y) + v2z*v3y - v2y*v3z + v1y*(-v2z + v3z)))/
   (2.*Sqrt(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2)));

}

void dnrm( double *g, double *pt1, double *pt2, double *pt3 )
{
	double v1x = pt1[0];
	double v1y = pt1[1];
	double v1z = pt1[2];
	
	double v2x = pt2[0];
	double v2y = pt2[1];
	double v2z = pt2[2];
	
	double v3x = pt3[0];
	double v3y = pt3[1];
	double v3z = pt3[2];
	
	// d nrm_x, v1x
	g[0] = -((-(v2z*v3y) + v1z*(-v2y + v3y) + v1y*(v2z - v3z) + v2y*v3z)*
      (2*(v2y - v3y)*(-(v2y*v3x) + v1y*(-v2x + v3x) + v1x*(v2y - v3y) + v2x*v3y) + 
        2*(v2z - v3z)*(-(v2z*v3x) + v1z*(-v2x + v3x) + v1x*(v2z - v3z) + v2x*v3z)))/
   (2.*Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5));
	g[1] = ((v1z*(v2x - v3x) + v2z*v3x - v2x*v3z + v1x*(-v2z + v3z))*
     (v1z*v2x*v2z + Power(v2y,2)*v3x - v1z*v2z*v3x + Power(v2z,2)*v3x + 
       v1y*(v2x - v3x)*(v2y - v3y) - v2x*v2y*v3y - v2y*v3x*v3y + v2x*Power(v3y,2) - 
       v1z*v2x*v3z - v2x*v2z*v3z + v1z*v3x*v3z - v2z*v3x*v3z + v2x*Power(v3z,2) - 
       v1x*(Power(v2y,2) + Power(v2z,2) - 2*v2y*v3y + Power(v3y,2) - 2*v2z*v3z + 
          Power(v3z,2))))/
   Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
     Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
     Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5);
	g[2] = -(((v1y*(v2x - v3x) + v2y*v3x - v2x*v3y + v1x*(-v2y + v3y))*
       (v1z*v2x*v2z + Power(v2y,2)*v3x - v1z*v2z*v3x + Power(v2z,2)*v3x + 
         v1y*(v2x - v3x)*(v2y - v3y) - v2x*v2y*v3y - v2y*v3x*v3y + v2x*Power(v3y,2) - 
         v1z*v2x*v3z - v2x*v2z*v3z + v1z*v3x*v3z - v2z*v3x*v3z + v2x*Power(v3z,2) - 
         v1x*(Power(v2y,2) + Power(v2z,2) - 2*v2y*v3y + Power(v3y,2) - 2*v2z*v3z + 
            Power(v3z,2))))/
     Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5));
	g[3] = -((-(v2z*v3y) + v1z*(-v2y + v3y) + v1y*(v2z - v3z) + v2y*v3z)*
      (2*(v1y - v3y)*(v1y*(v2x - v3x) + v2y*v3x - v2x*v3y + v1x*(-v2y + v3y)) + 
        2*(v1z - v3z)*(v1z*(v2x - v3x) + v2z*v3x - v2x*v3z + v1x*(-v2z + v3z))))/
   (2.*Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5));
	g[4] = -(((v1z*(v2x - v3x) + v2z*v3x - v2x*v3z + v1x*(-v2z + v3z))*
       (Power(v1y,2)*(v2x - v3x) + Power(v1z,2)*(v2x - v3x) + v1x*v2y*v3y - 
         v2y*v3x*v3y - v1x*Power(v3y,2) + v2x*Power(v3y,2) + 
         v1y*(v2y*v3x - 2*v2x*v3y + v3x*v3y + v1x*(-v2y + v3y)) + v1x*v2z*v3z - 
         v2z*v3x*v3z - v1x*Power(v3z,2) + v2x*Power(v3z,2) + 
         v1z*(-(v1x*v2z) + v2z*v3x + v1x*v3z - 2*v2x*v3z + v3x*v3z)))/
     Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5));
	g[5] = ((v1y*(v2x - v3x) + v2y*v3x - v2x*v3y + v1x*(-v2y + v3y))*
     (Power(v1y,2)*(v2x - v3x) + Power(v1z,2)*(v2x - v3x) + v1x*v2y*v3y - v2y*v3x*v3y - 
       v1x*Power(v3y,2) + v2x*Power(v3y,2) + 
       v1y*(v2y*v3x - 2*v2x*v3y + v3x*v3y + v1x*(-v2y + v3y)) + v1x*v2z*v3z - 
       v2z*v3x*v3z - v1x*Power(v3z,2) + v2x*Power(v3z,2) + 
       v1z*(-(v1x*v2z) + v2z*v3x + v1x*v3z - 2*v2x*v3z + v3x*v3z)))/
   Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
     Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
     Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5);
	g[6] = -((-(v2z*v3y) + v1z*(-v2y + v3y) + v1y*(v2z - v3z) + v2y*v3z)*
      (2*(v1y - v2y)*(-(v2y*v3x) + v1y*(-v2x + v3x) + v1x*(v2y - v3y) + v2x*v3y) + 
        2*(-v1z + v2z)*(v1z*(v2x - v3x) + v2z*v3x - v2x*v3z + v1x*(-v2z + v3z))))/
   (2.*Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5));
	g[7] = ((v1z*(v2x - v3x) + v2z*v3x - v2x*v3z + v1x*(-v2z + v3z))*
     (v1x*Power(v2y,2) + v1x*Power(v2z,2) + Power(v1y,2)*(v2x - v3x) + 
       Power(v1z,2)*(v2x - v3x) - Power(v2y,2)*v3x - Power(v2z,2)*v3x - v1x*v2y*v3y + 
       v2x*v2y*v3y - v1y*(-2*v2y*v3x + v1x*(v2y - v3y) + v2x*(v2y + v3y)) - 
       v1x*v2z*v3z + v2x*v2z*v3z - v1z*(-2*v2z*v3x + v1x*(v2z - v3z) + v2x*(v2z + v3z))))
    /Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
     Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
     Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5);
	g[8] = ((v1y*(v2x - v3x) + v2y*v3x - v2x*v3y + v1x*(-v2y + v3y))*
     (-(v1x*Power(v2y,2)) - v1x*Power(v2z,2) + Power(v2y,2)*v3x + Power(v2z,2)*v3x + 
       Power(v1y,2)*(-v2x + v3x) + Power(v1z,2)*(-v2x + v3x) + v1x*v2y*v3y - 
       v2x*v2y*v3y + v1y*(-2*v2y*v3x + v1x*(v2y - v3y) + v2x*(v2y + v3y)) + 
       v1x*v2z*v3z - v2x*v2z*v3z + v1z*(-2*v2z*v3x + v1x*(v2z - v3z) + v2x*(v2z + v3z))))
    /Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
     Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
     Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5);


	g[9] = -(((v1z*(v2y - v3y) + v2z*v3y - v2y*v3z + v1y*(-v2z + v3z))*
       (v1z*v2y*v2z - v2x*v2y*v3x + v2y*Power(v3x,2) + v1x*(v2x - v3x)*(v2y - v3y) + 
         Power(v2x,2)*v3y - v1z*v2z*v3y + Power(v2z,2)*v3y - v2x*v3x*v3y - v1z*v2y*v3z - 
         v2y*v2z*v3z + v1z*v3y*v3z - v2z*v3y*v3z + v2y*Power(v3z,2) - 
         v1y*(Power(v2x,2) + Power(v2z,2) - 2*v2x*v3x + Power(v3x,2) - 2*v2z*v3z + 
            Power(v3z,2))))/
     Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5));
	g[10] = -((v1z*(v2x - v3x) + v2z*v3x - v2x*v3z + v1x*(-v2z + v3z))*
      (2*(v2x - v3x)*(v1y*(v2x - v3x) + v2y*v3x - v2x*v3y + v1x*(-v2y + v3y)) + 
        2*(v2z - v3z)*(-(v2z*v3y) + v1z*(-v2y + v3y) + v1y*(v2z - v3z) + v2y*v3z)))/
   (2.*Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5));
	g[11] = ((v1y*(v2x - v3x) + v2y*v3x - v2x*v3y + v1x*(-v2y + v3y))*
     (-(v1z*v2y*v2z) + v2x*v2y*v3x - v2y*Power(v3x,2) - v1x*(v2x - v3x)*(v2y - v3y) - 
       Power(v2x,2)*v3y + v1z*v2z*v3y - Power(v2z,2)*v3y + v2x*v3x*v3y + v1z*v2y*v3z + 
       v2y*v2z*v3z - v1z*v3y*v3z + v2z*v3y*v3z - v2y*Power(v3z,2) + 
       v1y*(Power(v2x,2) + Power(v2z,2) - 2*v2x*v3x + Power(v3x,2) - 2*v2z*v3z + 
          Power(v3z,2))))/
   Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
     Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
     Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5);
	g[12] = ((v1z*(v2y - v3y) + v2z*v3y - v2y*v3z + v1y*(-v2z + v3z))*
     (v1y*v2x*v3x - v1y*Power(v3x,2) + v2y*Power(v3x,2) + Power(v1x,2)*(v2y - v3y) + 
       Power(v1z,2)*(v2y - v3y) - v2x*v3x*v3y + 
       v1x*(-2*v2y*v3x + v1y*(-v2x + v3x) + (v2x + v3x)*v3y) + v1y*v2z*v3z - 
       v2z*v3y*v3z - v1y*Power(v3z,2) + v2y*Power(v3z,2) + 
       v1z*(v2z*v3y + (-2*v2y + v3y)*v3z + v1y*(-v2z + v3z))))/
   Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
     Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
     Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5);
	g[13] = -((v1z*(v2x - v3x) + v2z*v3x - v2x*v3z + v1x*(-v2z + v3z))*
      (2*(v1x - v3x)*(-(v2y*v3x) + v1y*(-v2x + v3x) + v1x*(v2y - v3y) + v2x*v3y) + 
        2*(v1z - v3z)*(v1z*(v2y - v3y) + v2z*v3y - v2y*v3z + v1y*(-v2z + v3z))))/
   (2.*Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5));
	g[14] = ((v1y*(v2x - v3x) + v2y*v3x - v2x*v3y + v1x*(-v2y + v3y))*
     (v1y*v2x*v3x - v1y*Power(v3x,2) + v2y*Power(v3x,2) + Power(v1x,2)*(v2y - v3y) + 
       Power(v1z,2)*(v2y - v3y) - v2x*v3x*v3y + 
       v1x*(-2*v2y*v3x + v1y*(-v2x + v3x) + (v2x + v3x)*v3y) + v1y*v2z*v3z - 
       v2z*v3y*v3z - v1y*Power(v3z,2) + v2y*Power(v3z,2) + 
       v1z*(v2z*v3y + (-2*v2y + v3y)*v3z + v1y*(-v2z + v3z))))/
   Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
     Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
     Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5);
	g[15] = -(((Power(v1z,2)*v2y - v1z*v2y*v2z + v2x*v2y*v3x + Power(v1x,2)*(v2y - v3y) - 
         Power(v1z,2)*v3y - Power(v2x,2)*v3y + 2*v1z*v2z*v3y - Power(v2z,2)*v3y - 
         v1x*(v2x*v2y + v1y*(v2x - v3x) + v2y*v3x - 2*v2x*v3y) + 
         v1y*(Power(v2x,2) - v2x*v3x - (v1z - v2z)*(v2z - v3z)) - v1z*v2y*v3z + 
         v2y*v2z*v3z)*(v1z*(v2y - v3y) + v2z*v3y - v2y*v3z + v1y*(-v2z + v3z)))/
     Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5));
	g[16] = -((v1z*(v2x - v3x) + v2z*v3x - v2x*v3z + v1x*(-v2z + v3z))*
      (2*(-v1x + v2x)*(-(v2y*v3x) + v1y*(-v2x + v3x) + v1x*(v2y - v3y) + v2x*v3y) + 
        2*(v1z - v2z)*(-(v2z*v3y) + v1z*(-v2y + v3y) + v1y*(v2z - v3z) + v2y*v3z)))/
   (2.*Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5));
	g[17] = -(((v1y*(v2x - v3x) + v2y*v3x - v2x*v3y + v1x*(-v2y + v3y))*
       (Power(v1z,2)*v2y - v1z*v2y*v2z + v2x*v2y*v3x + Power(v1x,2)*(v2y - v3y) - 
         Power(v1z,2)*v3y - Power(v2x,2)*v3y + 2*v1z*v2z*v3y - Power(v2z,2)*v3y - 
         v1x*(v2x*v2y + v1y*(v2x - v3x) + v2y*v3x - 2*v2x*v3y) + 
         v1y*(Power(v2x,2) - v2x*v3x - (v1z - v2z)*(v2z - v3z)) - v1z*v2y*v3z + 
         v2y*v2z*v3z))/
     Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5));
	// d nrm_z /d v1x
	g[18] = ((-(v1y*v2y*v2z) + v2x*v2z*v3x - v2z*Power(v3x,2) + v1y*v2z*v3y + v2y*v2z*v3y - 
       v2z*Power(v3y,2) + v1z*(Power(v2x,2) + Power(v2y,2) - 2*v2x*v3x + Power(v3x,2) - 
          2*v2y*v3y + Power(v3y,2)) - v1x*(v2x - v3x)*(v2z - v3z) - Power(v2x,2)*v3z + 
       v1y*v2y*v3z - Power(v2y,2)*v3z + v2x*v3x*v3z - v1y*v3y*v3z + v2y*v3y*v3z)*
     (v1z*(v2y - v3y) + v2z*v3y - v2y*v3z + v1y*(-v2z + v3z)))/
   Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
     Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
     Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5);
 	g[19] = -(((-(v1y*v2y*v2z) + v2x*v2z*v3x - v2z*Power(v3x,2) + v1y*v2z*v3y + v2y*v2z*v3y - 
         v2z*Power(v3y,2) + v1z*(Power(v2x,2) + Power(v2y,2) - 2*v2x*v3x + 
            Power(v3x,2) - 2*v2y*v3y + Power(v3y,2)) - v1x*(v2x - v3x)*(v2z - v3z) - 
         Power(v2x,2)*v3z + v1y*v2y*v3z - Power(v2y,2)*v3z + v2x*v3x*v3z - v1y*v3y*v3z + 
         v2y*v3y*v3z)*(v1z*(v2x - v3x) + v2z*v3x - v2x*v3z + v1x*(-v2z + v3z)))/
     Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5));
	g[20] = -((-(v2y*v3x) + v1y*(-v2x + v3x) + v1x*(v2y - v3y) + v2x*v3y)*
      (2*(v2x - v3x)*(v1z*(v2x - v3x) + v2z*v3x - v2x*v3z + v1x*(-v2z + v3z)) + 
        2*(v2y - v3y)*(v1z*(v2y - v3y) + v2z*v3y - v2y*v3z + v1y*(-v2z + v3z))))/
   (2.*Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5));
	g[21] = ((v1z*(v2y - v3y) + v2z*v3y - v2y*v3z + v1y*(-v2z + v3z))*
     (v1z*v2x*v3x - v1z*Power(v3x,2) + v2z*Power(v3x,2) + v1z*v2y*v3y - 
       v1z*Power(v3y,2) + v2z*Power(v3y,2) + Power(v1x,2)*(v2z - v3z) + 
       Power(v1y,2)*(v2z - v3z) - v2x*v3x*v3z - v2y*v3y*v3z + 
       v1x*(-2*v2z*v3x + v1z*(-v2x + v3x) + (v2x + v3x)*v3z) + 
       v1y*(-2*v2z*v3y + v1z*(-v2y + v3y) + (v2y + v3y)*v3z)))/
   Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
     Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
     Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5);
	g[22] = ((v1z*(v2x - v3x) + v2z*v3x - v2x*v3z + v1x*(-v2z + v3z))*
     (-(v1z*v2x*v3x) + v1z*Power(v3x,2) - v2z*Power(v3x,2) - v1z*v2y*v3y + 
       v1z*Power(v3y,2) - v2z*Power(v3y,2) + v2x*v3x*v3z + v2y*v3y*v3z + 
       Power(v1x,2)*(-v2z + v3z) + Power(v1y,2)*(-v2z + v3z) + 
       v1x*(v1z*(v2x - v3x) + 2*v2z*v3x - (v2x + v3x)*v3z) + 
       v1y*(v1z*(v2y - v3y) + 2*v2z*v3y - (v2y + v3y)*v3z)))/
   Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
     Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
     Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5);
	g[23] = -((-(v2y*v3x) + v1y*(-v2x + v3x) + v1x*(v2y - v3y) + v2x*v3y)*
      (2*(v1x - v3x)*(-(v2z*v3x) + v1z*(-v2x + v3x) + v1x*(v2z - v3z) + v2x*v3z) + 
        2*(v1y - v3y)*(-(v2z*v3y) + v1z*(-v2y + v3y) + v1y*(v2z - v3z) + v2y*v3z)))/
   (2.*Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5));
	g[24] = ((-(v2z*v3y) + v1z*(-v2y + v3y) + v1y*(v2z - v3z) + v2y*v3z)*
     (Power(v1y,2)*v2z - v1y*v2y*v2z + v2x*v2z*v3x + 
       v1z*(Power(v2x,2) - v2x*v3x - (v1y - v2y)*(v2y - v3y)) - v1y*v2z*v3y + 
       v2y*v2z*v3y + Power(v1x,2)*(v2z - v3z) - Power(v1y,2)*v3z - Power(v2x,2)*v3z + 
       2*v1y*v2y*v3z - Power(v2y,2)*v3z - 
       v1x*(v2x*v2z + v1z*(v2x - v3x) + v2z*v3x - 2*v2x*v3z)))/
   Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
     Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
     Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5);
	g[25] = ((v1z*(v2x - v3x) + v2z*v3x - v2x*v3z + v1x*(-v2z + v3z))*
     (Power(v1y,2)*v2z - v1y*v2y*v2z + v2x*v2z*v3x + 
       v1z*(Power(v2x,2) - v2x*v3x - (v1y - v2y)*(v2y - v3y)) - v1y*v2z*v3y + 
       v2y*v2z*v3y + Power(v1x,2)*(v2z - v3z) - Power(v1y,2)*v3z - Power(v2x,2)*v3z + 
       2*v1y*v2y*v3z - Power(v2y,2)*v3z - 
       v1x*(v2x*v2z + v1z*(v2x - v3x) + v2z*v3x - 2*v2x*v3z)))/
   Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
     Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
     Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5);
	g[26] = -((-(v2y*v3x) + v1y*(-v2x + v3x) + v1x*(v2y - v3y) + v2x*v3y)*
      (2*(-v1y + v2y)*(-(v2z*v3y) + v1z*(-v2y + v3y) + v1y*(v2z - v3z) + v2y*v3z) + 
        2*(v1x - v2x)*(v1z*(v2x - v3x) + v2z*v3x - v2x*v3z + v1x*(-v2z + v3z))))/
   (2.*Power(Power(v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y,2) + 
       Power(v1z*v2x - v1x*v2z - v1z*v3x + v2z*v3x + v1x*v3z - v2x*v3z,2) + 
       Power(v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z,2),1.5));
}

double tvol( double *pt1, double *pt2, double *pt3, double *pt4 )
{
	double v1x = pt1[0];
	double v1y = pt1[1];
	double v1z = pt1[2];
	
	double v2x = pt2[0];
	double v2y = pt2[1];
	double v2z = pt2[2];
	
	double v3x = pt3[0];
	double v3y = pt3[1];
	double v3z = pt3[2];
	
	double v4x = pt4[0];
	double v4y = pt4[1];
	double v4z = pt4[2];

	double vol = ((-(v2y*v3x) + v2x*v3y + v2y*v4x - v3y*v4x - v2x*v4y + v3x*v4y)*(v1z - v4z) + 
     (v1y - v4y)*(v2z*v3x - v2x*v3z - v2z*v4x + v3z*v4x + v2x*v4z - v3x*v4z) + 
     (v1x - v4x)*(-(v2z*v3y) + v2y*v3z + v2z*v4y - v3z*v4y - v2y*v4z + v3y*v4z))/6.;

	return vol;
}

void d_tvol( double *der, double *pt1, double *pt2, double *pt3, double *pt4 )
{
	double v1x = pt1[0];
	double v1y = pt1[1];
	double v1z = pt1[2];
	
	double v2x = pt2[0];
	double v2y = pt2[1];
	double v2z = pt2[2];
	
	double v3x = pt3[0];
	double v3y = pt3[1];
	double v3z = pt3[2];
	
	double v4x = pt4[0];
	double v4y = pt4[1];
	double v4z = pt4[2];
	double vol = ((-(v2y*v3x) + v2x*v3y + v2y*v4x - v3y*v4x - v2x*v4y + v3x*v4y)*(v1z - v4z) + 
     (v1y - v4y)*(v2z*v3x - v2x*v3z - v2z*v4x + v3z*v4x + v2x*v4z - v3x*v4z) + 
     (v1x - v4x)*(-(v2z*v3y) + v2y*v3z + v2z*v4y - v3z*v4y - v2y*v4z + v3y*v4z))/6.;


	double sgn = 1;

//	if( vol < 0 )
//		sgn = -1;
	if( fabs(vol) < 1e-5 )
	{
//		printf("Very small vol.\n");
//		exit(1);
	}
	der[0] = sgn * ((-(v2z*v3y) + v2y*v3z + v2z*v4y - v3z*v4y - v2y*v4z + v3y*v4z)/6. );
	der[1] = sgn * ((v2z*v3x - v2x*v3z - v2z*v4x + v3z*v4x + v2x*v4z - v3x*v4z)/6. );
	der[2] = sgn * ((-(v2y*v3x) + v2x*v3y + v2y*v4x - v3y*v4x - v2x*v4y + v3x*v4y)/6. );

	der[3] = sgn * (((v3y - v4y)*(v1z - v4z) + (v1y - v4y)*(-v3z + v4z))/6. );
	der[4] = sgn * (((-v3x + v4x)*(v1z - v4z) + (v1x - v4x)*(v3z - v4z))/6. );
	der[5] = sgn * (((v3x - v4x)*(v1y - v4y) + (v1x - v4x)*(-v3y + v4y))/6. );

	der[6] = sgn * (((-v2y + v4y)*(v1z - v4z) + (v1y - v4y)*(v2z - v4z))/6. );
	der[7] = sgn * (((v2x - v4x)*(v1z - v4z) + (v1x - v4x)*(-v2z + v4z))/6. );
	der[8] = sgn * (((-v2x + v4x)*(v1y - v4y) + (v1x - v4x)*(v2y - v4y))/6. );

	der[9] = sgn * ((v2z*v3y - v2y*v3z + (-v2z + v3z)*(v1y - v4y) - v2z*v4y + v3z*v4y + 
     (v2y - v3y)*(v1z - v4z) + v2y*v4z - v3y*v4z)/6. );
	der[10] = sgn * ((-(v2z*v3x) + v2x*v3z + (v2z - v3z)*(v1x - v4x) + v2z*v4x - v3z*v4x + 
     (-v2x + v3x)*(v1z - v4z) - v2x*v4z + v3x*v4z)/6. );
	der[11] = sgn * ((v2y*v3x - v2x*v3y + (-v2y + v3y)*(v1x - v4x) - v2y*v4x + v3y*v4x + 
     (v2x - v3x)*(v1y - v4y) + v2x*v4y - v3x*v4y)/6. );

}

int pointInTriangle( double *pt, double *pt1, double *pt2, double *pt3, double tol  )
{

	double t1x = pt[0];
	double t1y = pt[1];
	double t1z = 0;
	
	double p1x = pt1[0];
	double p1y = pt1[1];
	double p1z = 0;

	double p2x = pt2[0];
	double p2y = pt2[1];
	double p2z = 0;

	double p3x = pt3[0];
	double p3y = pt3[1];
	double p3z = 0;

double sgn1 = (p1y*(p2x - t1x) + p2y*t1x - p2x*t1y + p1x*(-p2y + t1y))*(p1y*(p3x - t1x) + p3y*t1x - p3x*t1y + p1x*(-p3y + t1y)) + 
   (-(p2z*t1x) + p1z*(-p2x + t1x) + p1x*(p2z - t1z) + p2x*t1z)*
    (-(p3z*t1x) + p1z*(-p3x + t1x) + p1x*(p3z - t1z) + p3x*t1z) + 
   (p1z*(p2y - t1y) + p2z*t1y - p2y*t1z + p1y*(-p2z + t1z))*(p1z*(p3y - t1y) + p3z*t1y - p3y*t1z + p1y*(-p3z + t1z));
		
double sgn2 = (-(p3y*t1x) + p2y*(-p3x + t1x) + p2x*(p3y - t1y) + p3x*t1y)*(p1y*(p2x - t1x) + p2y*t1x - p2x*t1y + p1x*(-p2y + t1y)) + 
   (-(p3z*t1x) + p2z*(-p3x + t1x) + p2x*(p3z - t1z) + p3x*t1z)*
    (p1z*(p2x - t1x) + p2z*t1x - p2x*t1z + p1x*(-p2z + t1z)) + 
   (-(p3z*t1y) + p2z*(-p3y + t1y) + p2y*(p3z - t1z) + p3y*t1z)*(p1z*(p2y - t1y) + p2z*t1y - p2y*t1z + p1y*(-p2z + t1z));

double sgn3 = (p1y*(p3x - t1x) + p3y*t1x - p3x*t1y + p1x*(-p3y + t1y))*(p2y*(p3x - t1x) + p3y*t1x - p3x*t1y + p2x*(-p3y + t1y)) + 
   (p1z*(p3x - t1x) + p3z*t1x - p3x*t1z + p1x*(-p3z + t1z))*(p2z*(p3x - t1x) + p3z*t1x - p3x*t1z + p2x*(-p3z + t1z)) + 
   (p1z*(p3y - t1y) + p3z*t1y - p3y*t1z + p1y*(-p3z + t1z))*(p2z*(p3y - t1y) + p3z*t1y - p3y*t1z + p2y*(-p3z + t1z));	

	if( sgn1 < tol && sgn2 < tol && sgn3 < tol )
		return 1;
	return 0;
}


void MinImage3D( double *dr1, double PBC_vec[3][3], double *put_vec, double *alphas )
{
	int done = 0;

	put_vec[0] = 0;
	put_vec[1] = 0;
	put_vec[2] = 0;

	double use_alpha[3] = {1,1,1};

	if( alphas )
	{
		use_alpha[0] = alphas[0];
		use_alpha[1] = alphas[1];
		use_alpha[2] = alphas[2];
	}

	double tr[3] = { dr1[0], dr1[1], dr1[2] };

	while( !done )
	{
		done = 1;

		double r2 = tr[0]*tr[0]+tr[1]*tr[1]+tr[2]*tr[2];

		for( int dx = -1; dx <= 1; dx++ )
		for( int dy = -1; dy <= 1; dy++ )
		for( int dz = -1; dz <= 1; dz++ )
		{
			tr[0] += (dx * PBC_vec[0][0] + dy * PBC_vec[1][0] + dz * PBC_vec[2][0])*use_alpha[0];
			tr[1] += (dx * PBC_vec[0][1] + dy * PBC_vec[1][1] + dz * PBC_vec[2][1])*use_alpha[1];
			tr[2] += (dx * PBC_vec[0][2] + dy * PBC_vec[1][2] + dz * PBC_vec[2][2])*use_alpha[2];

			double nr2 = tr[0]*tr[0]+tr[1]*tr[1]+tr[2]*tr[2];

			if( nr2 >= r2 - 1e-7 )
			{
				tr[0] -= (dx * PBC_vec[0][0] + dy * PBC_vec[1][0] + dz * PBC_vec[2][0])*use_alpha[0];
				tr[1] -= (dx * PBC_vec[0][1] + dy * PBC_vec[1][1] + dz * PBC_vec[2][1])*use_alpha[1];
				tr[2] -= (dx * PBC_vec[0][2] + dy * PBC_vec[1][2] + dz * PBC_vec[2][2])*use_alpha[2];

			}
			else	
			{
				r2 = nr2;

				put_vec[0] += dx;
				put_vec[1] += dy;
				put_vec[2] += dz;
				
				done = 0;
			}
		}
	}

	dr1[0] = tr[0];
	dr1[1] = tr[1];
	dr1[2] = tr[2];
}

void fprintPSFAtomExt( FILE *theFile, int atn, int res, const char *atname, const char *resname)
{
	char tbuf[256];

	sprintf(tbuf, "%d", atn );

	for( int d = 0; d < 10 - strlen(tbuf); d++ )
		fprintf(theFile, " ");
	fprintf(theFile, "%s", tbuf );

	fprintf(theFile, " SEGN     ");
	
	sprintf(tbuf, "%d", res );

	fprintf(theFile, "%s", tbuf );
	for( int d = 0; d < 9 - strlen(tbuf); d++ )
		fprintf(theFile, " ");

	//fprintf(theFile, "  1    ");

	fprintf(theFile, "%s", resname );

	fprintf(theFile, "        %c        %c        0.000000        1.0000           0   0.00000      0.000000E-01\n", atname[0], atname[0] );
}

void writePSF( FILE *thePSF, int nv, char *atomNames, int *bond_list, int nbonds )
{
	fprintf(thePSF, "PSF EXT\n");
	fprintf(thePSF, "\n");
	
	char tbuf[256];
	sprintf(tbuf, "%d", nv );
	
	for( int x = 0; x < 8 - strlen(tbuf); x++ )
		fprintf(thePSF, " ");
	fprintf(thePSF, "%s", tbuf );	
	
	fprintf(thePSF, " !NATOM\n");

	for( int i = 0; i < nv; i++ )
	{
		if( atomNames )
		{
			char name[2] = { atomNames[i], '\0' }; 
			fprintPSFAtomExt( thePSF, 1+i, 1+i, name, name );
		}
		else
			fprintPSFAtomExt( thePSF, 1+i, 1+i, "C", "C" );
	}
	sprintf(tbuf, "%d", nbonds );
	
	for( int x = 0; x < 8 - strlen(tbuf); x++ )
		fprintf(thePSF, " ");
	fprintf(thePSF, "%s", tbuf );	
	fprintf(thePSF, " !NBOND: nbonds\n");
	
	int cr = 0;
	for( int b = 0; b < nbonds; b++ )
	{
		sprintf(tbuf, "%d", 1+bond_list[2*b+0] );
	
		for( int x = 0; x < 10 - strlen(tbuf); x++ )
			fprintf(thePSF, " ");
		fprintf(thePSF, "%s", tbuf );	
	
		sprintf(tbuf, "%d", 1+bond_list[2*b+1]);
	
		for( int x = 0; x < 10 - strlen(tbuf); x++ )
			fprintf(thePSF, " ");
		fprintf(thePSF, "%s", tbuf );	
	
		cr++;
	
		if( cr % 4 == 0 )
			fprintf(thePSF, "\n");
	}

}


int sort3( int *ind )
{ // first to third

	int code = 0; // i j k -> ik
	if( ind[0] > ind[1] )
	{
		int t = ind[0];	
		ind[0] = ind[1];
		ind[1] = t;
		code = 1;
		// j i k  -> jk
	}
	if( ind[1] > ind[2] )
	{
		int t = ind[1];	
		ind[1] = ind[2];
		ind[2] = t;

		// code 2: i k j (2)
		// code 3: j k i (0)

		code += 2;
	}
	if( ind[0] > ind[1] )
	{
		int t = ind[0];	
		ind[0] = ind[1];
		ind[1] = t;
		// code 4: k i j (2)
		// code 5: k j i (1)
		code += 2;
	}

	return code;
}

int getcode( int *verts, int *sorted )
{
	int i = verts[0];
	int j = verts[1];

	if( j < i )
	{
		int t = i;
		i = j;
		j = t;
	}

	if( i == sorted[0] && j == sorted[1] )
		return 2;
	if( i == sorted[0] && j == sorted[2] )
		return 1;
	if( i == sorted[1] && j == sorted[2] )
		return 0;
	return -1;
}

int iabs( int del )
{
	if( del < 0 )
		return -del;
	return del;
}

double Power( double a, int i )
{
	return pow(a,(double)i);
}

double Mul22( double *A, double *B, double *out )
{
	double O11 = A[0] * B[0] + A[1] * B[2];
	double O12 = A[0] * B[1] + A[1] * B[3];
	double O21 = A[2] * B[0] + A[3] * B[2];
	double O22 = A[2] * B[1] + A[3] * B[3];

	out[0] = O11;
	out[1] = O12;
	out[2] = O21;
	out[3] = O22;

	return 0;
}

void MatVec( double *a, double *b, double *c, int n, int m)
{
	// A:  m fast, n slow
	// B: n
	// C: m
	
	char transy='N';	
	double one =1.0;
	double zero = 0.0;
	int incr = 1;
	dgemv( &transy, &m, &n, &one, a, &m, b, &incr, &zero, c, &incr );
}

int nearInteriorPointOnTriangle( double *test_pt, double *vert1, double *vert2, double *vert3, double *output)
{
	double vec1[3] = { vert1[0] - vert2[0],
			   vert1[1] - vert2[1],
			   vert1[2] - vert2[2] };
	
	double vec2[3] = { vert3[0] - vert2[0],
			   vert3[1] - vert2[1],
			   vert3[2] - vert2[2] };
	
	double vec3[3] = { vert3[0] - vert1[0],
			   vert3[1] - vert1[1],
			   vert3[2] - vert1[2] };
	
	// triangle is convex so this must be interior point.
	double interior_pt[3] = { 
		(vert1[0]+vert2[0]+vert3[0])/3,
		(vert1[1]+vert2[1]+vert3[1])/3,
		(vert1[2]+vert2[2]+vert3[2])/3 };

	double nrm[3];
	cross( vec1, vec2, nrm );
	normalize(nrm);

	double dp = nrm[0] * vert2[0] + nrm[1] * vert2[1] + nrm[2] * vert2[2];
	// defines the tangent plane.
	
	double dpp = nrm[0] * test_pt[0] + nrm[1] * test_pt[1] + nrm[2] * test_pt[2];
	double dist_to_plane = fabs(dpp-dp);
	
	output[0] = test_pt[0] + (dp-dpp) * nrm[0];
	output[1] = test_pt[1] + (dp-dpp) * nrm[1];
	output[2] = test_pt[2] + (dp-dpp) * nrm[2];
	
	// is the output interior?

	double check_vec1[3];
	cross( vec1, nrm, check_vec1 );
	double tvec1[3] = { interior_pt[0] - vert1[0], interior_pt[1] - vert1[1], interior_pt[2] - vert1[2] };
	double cvec1[3] = { output[0] - vert1[0], output[1] - vert1[1], output[2] - vert1[2] };
	double sign1 = tvec1[0] * check_vec1[0] + tvec1[1] * check_vec1[1] + tvec1[2] * check_vec1[2];
	double csign1 = cvec1[0] * check_vec1[0] + cvec1[1] * check_vec1[1] + cvec1[2] * check_vec1[2];
	if( csign1 * sign1 < 0 ) return 0;
	
	double check_vec2[3];
	cross( vec2, nrm, check_vec2 );
	double tvec2[3] = { interior_pt[0] - vert2[0], interior_pt[1] - vert2[1], interior_pt[2] - vert2[2] };
	double cvec2[3] = { output[0] - vert2[0], output[1] - vert2[1], output[2] - vert2[2] };
	double sign2 = tvec2[0] * check_vec2[0] + tvec2[1] * check_vec2[1] + tvec2[2] * check_vec2[2];
	double csign2 = cvec2[0] * check_vec2[0] + cvec2[1] * check_vec2[1] + cvec2[2] * check_vec2[2];
	if( csign2 * sign2 < 0 ) return 0;
	
	double check_vec3[3];
	cross( vec3, nrm, check_vec3 );
	double tvec3[3] = { interior_pt[0] - vert3[0], interior_pt[1] - vert3[1], interior_pt[2] - vert3[2] };
	double cvec3[3] = { output[0] - vert3[0], output[1] - vert3[1], output[2] - vert3[2] };
	double sign3 = tvec3[0] * check_vec3[0] + tvec3[1] * check_vec3[1] + tvec3[2] * check_vec3[2];
	double csign3 = cvec3[0] * check_vec3[0] + cvec3[1] * check_vec3[1] + cvec3[2] * check_vec3[2];
	if( csign3 * sign3 < 0 ) return 0;

	return 1;
}

double segmentSegmentDist( double *r1A, double *r1B, double *r2A, double *r2B, double *t1_out, double *t2_out)
{
	double r1Ax = r1A[0];	
	double r1Ay = r1A[1];	
	double r1Az = r1A[2];	
	
	double r1Bx = r1B[0];	
	double r1By = r1B[1];	
	double r1Bz = r1B[2];	
	
	double r2Ax = r2A[0];	
	double r2Ay = r2A[1];	
	double r2Az = r2A[2];	
	
	double r2Bx = r2B[0];	
	double r2By = r2B[1];	
	double r2Bz = r2B[2];	

	double t1 = -(((2*r1Ax*(r2Ax - r2Bx) - 2*r2Ax*(r2Ax - r2Bx) + 2*r1Ay*(r2Ay - r2By) - 
          2*r2Ay*(r2Ay - r2By) + 2*r1Az*(r2Az - r2Bz) - 2*r2Az*(r2Az - r2Bz))*
        (-2*(r1Ax - r1Bx)*r2Ax - 2*(r1Ay - r1By)*r2Ay - 2*(r1Az - r1Bz)*r2Az + 
          2*(r1Ax - r1Bx)*r2Bx + 2*(r1Ay - r1By)*r2By + 2*(r1Az - r1Bz)*r2Bz) - 
       (-2*r1Ax*(r1Ax - r1Bx) - 2*r1Ay*(r1Ay - r1By) - 2*r1Az*(r1Az - r1Bz) + 
          2*(r1Ax - r1Bx)*r2Ax + 2*(r1Ay - r1By)*r2Ay + 2*(r1Az - r1Bz)*r2Az)*
        (2*r2Ax*(r2Ax - r2Bx) - 2*(r2Ax - r2Bx)*r2Bx + 2*r2Ay*(r2Ay - r2By) - 
          2*(r2Ay - r2By)*r2By + 2*r2Az*(r2Az - r2Bz) - 2*(r2Az - r2Bz)*r2Bz))/
     ((-2*r1Ax*(r2Ax - r2Bx) + 2*r1Bx*(r2Ax - r2Bx) - 2*r1Ay*(r2Ay - r2By) + 
          2*r1By*(r2Ay - r2By) - 2*r1Az*(r2Az - r2Bz) + 2*r1Bz*(r2Az - r2Bz))*
        (-2*(r1Ax - r1Bx)*r2Ax - 2*(r1Ay - r1By)*r2Ay - 2*(r1Az - r1Bz)*r2Az + 
          2*(r1Ax - r1Bx)*r2Bx + 2*(r1Ay - r1By)*r2By + 2*(r1Az - r1Bz)*r2Bz) - 
       (2*r1Ax*(r1Ax - r1Bx) - 2*(r1Ax - r1Bx)*r1Bx + 2*r1Ay*(r1Ay - r1By) - 
          2*(r1Ay - r1By)*r1By + 2*r1Az*(r1Az - r1Bz) - 2*(r1Az - r1Bz)*r1Bz)*
        (2*r2Ax*(r2Ax - r2Bx) - 2*(r2Ax - r2Bx)*r2Bx + 2*r2Ay*(r2Ay - r2By) - 
          2*(r2Ay - r2By)*r2By + 2*r2Az*(r2Az - r2Bz) - 2*(r2Az - r2Bz)*r2Bz)));

	double expr1 = Power(r2Ax,2) + Power(r2Ay,2) + Power(r2Az,2) - 2*r2Ax*r2Bx + Power(r2Bx,2) - 
   2*r2Ay*r2By + Power(r2By,2) - 2*r2Az*r2Bz + Power(r2Bz,2);

	double t2 = -((r1Ax*r2Ax)/expr1) + Power(r2Ax,2)/expr1 - (r1Ay*r2Ay)/expr1 + Power(r2Ay,2)/expr1 - 
   (r1Az*r2Az)/expr1 + Power(r2Az,2)/expr1 + (r1Ax*r2Bx)/expr1 - (r2Ax*r2Bx)/expr1 + 
   (r1Ay*r2By)/expr1 - (r2Ay*r2By)/expr1 + (r1Az*r2Bz)/expr1 - (r2Az*r2Bz)/expr1 + 
   (r1Ax*r2Ax*t1)/expr1 - (r1Bx*r2Ax*t1)/expr1 + (r1Ay*r2Ay*t1)/expr1 - 
   (r1By*r2Ay*t1)/expr1 + (r1Az*r2Az*t1)/expr1 - (r1Bz*r2Az*t1)/expr1 - 
   (r1Ax*r2Bx*t1)/expr1 + (r1Bx*r2Bx*t1)/expr1 - (r1Ay*r2By*t1)/expr1 + 
   (r1By*r2By*t1)/expr1 - (r1Az*r2Bz*t1)/expr1 + (r1Bz*r2Bz*t1)/expr1;

	double p1[3] = { 
		r1Ax + t1 * (r1Bx-r1Ax),
		r1Ay + t1 * (r1By-r1Ay),
		r1Az + t1 * (r1Bz-r1Az) };
	double p2[3] = { 
		r2Ax + t1 * (r2Bx-r2Ax),
		r2Ay + t1 * (r2By-r2Ay),
		r2Az + t1 * (r2Bz-r2Az) };
	double dr[3] = { p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2] };

	*t1_out = t1;
	*t2_out = t2;

	return normalize(dr);
}


double dihe( double *r1, double *r2, double *r3, double *r4 ) 
{
        double b1[3] = {r2[0] - r1[0], r2[1] -r1[1], r2[2] - r1[2] };
        double b2[3] = {r3[0] - r2[0], r3[1] - r2[1], r3[2] - r2[2] };
        double b3[3] = {r4[0] - r3[0], r4[1] - r3[1], r4[2] - r3[2] };

        double lb2 = sqrt(b2[0]*b2[0]+b2[1]*b2[1]+b2[2]*b2[2]);
        double b2xb3[3] = {  b2[1] * b3[2] - b2[2] * b3[1],
                           -(b2[0] * b3[2] - b2[2] * b3[0]),
                             b2[0] * b3[1] - b2[1] * b3[0] };
        double b1xb2[3] = {  b1[1] * b2[2] - b1[2] * b2[1],
                           -(b1[0] * b2[2] - b1[2] * b2[0]),
                             b1[0] * b2[1] - b1[1] * b2[0] };

        double b1db2xb3 = b1[0] * b2xb3[0] + b1[1] * b2xb3[1] + b1[2] * b2xb3[2];
        double xx = b1xb2[0] * b2xb3[0] + b1xb2[1] * b2xb3[1] + b1xb2[2] * b2xb3[2];

        double arg1= lb2 * b1db2xb3;
        double arg2 = xx; 

        if( fabs( arg1 ) < 1e-10 )
                arg1 = 1e-10;
        if( fabs( arg2 ) < 1e-10 )
                arg2 = 1e-10;   
 
              double phi = atan2( arg1, arg2 );
                      return phi;
}
