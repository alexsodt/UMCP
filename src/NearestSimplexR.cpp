#include "functions.h"
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define TYPE_INTERIOR   0
#define TYPE_EDGE       1
#define TYPE_POINT      2

Check Functions::NearestSimplex()
{
	Point vec10 = simplex[0] - simplex[1];
	Point vec20 = simplex[0] - simplex[2];
	Point vec21 = simplex[1] - simplex[2];
	Point vec31 = simplex[1] - simplex[3];
	Point vec32 = simplex[2] - simplex[3];
	Point vec02 = simplex[2] - simplex[0];
	Point vec03 = simplex[3] - simplex[0];
	Point vec13 = simplex[3] - simplex[1];
	Point n1 = cross( vec10, vec20 );
	Point n2 = cross( vec21, vec31 );
	Point n3 = cross( vec32, vec02 );
	Point n4 = cross( vec03, vec13 );
	
	double k1 = dot(simplex[0],n1);
	double k2 = dot(simplex[1],n2);
	double k3 = dot(simplex[2],n3);
	double k4 = dot(simplex[3],n4);
	
	double l1 = dot(simplex[3],n1);
	double l2 = dot(simplex[0],n2);
	double l3 = dot(simplex[1],n3);
	double l4 = dot(simplex[2],n4);

	int con1 = (k1 > 0 && l1 < k1) || ( k1 < 0 && l1 > k1 );
	int con2 = (k2 > 0 && l2 < k2) || ( k2 < 0 && l2 > k2 );
	int con3 = (k3 > 0 && l3 < k3) || ( k3 < 0 && l3 > k3 );
	int con4 = (k4 > 0 && l4 < k4) || ( k4 < 0 && l4 > k4 );

	if( con1 && con2 && con3 && con4 ) {
		Check new_check = {true, vec10};
		return new_check;
	}
 
	if( elements != 4 )
	{
		printf("NearestSimplex must be called with a tetrahedron currently.\n");
		exit(1);
	}

	int which1,type1;	
	int which2,type2;	
	int which3,type3;	
	int which4,type4;	

	// edges: 01 02 03 12 13 23
	//         0  1  2  3  4  5
	// edges:
	// 12 02 01	
	int emap1[3] = { 3, 1, 0 };	
	int pmap1[3] = { 0, 1, 2 };
	double d1 = near_point_on_triangle( 
		simplex[0].x, simplex[0].y, simplex[0].z, // u = 1, v = 0 
		simplex[1].x, simplex[1].y, simplex[1].z, // u = 0, v = 1
		simplex[2].x, simplex[2].y, simplex[2].z, // u = 0, v = 0
		&type1, &which1 );
	if( type1 == TYPE_EDGE) which1 = emap1[which1];
	if( type1 == TYPE_POINT) which1 = pmap1[which1];
	// 01 13 03
	int emap2[3] = { 0, 4, 2 };
	int pmap2[3] = { 3, 0, 1 };
	double d2 = near_point_on_triangle( 
		simplex[3].x, simplex[3].y, simplex[3].z, 
		simplex[0].x, simplex[0].y, simplex[0].z, 
		simplex[1].x, simplex[1].y, simplex[1].z, &type2, &which2 );
	if( type2 == TYPE_EDGE) which2 = emap2[which2];
	if( type2 == TYPE_POINT) which2 = pmap2[which2];
	// 03 02 23
	int emap3[3] = { 2, 1, 5 };
	int pmap3[3] = { 0, 2, 3 };
	double d3 = near_point_on_triangle(                    
		simplex[2].x, simplex[2].y, simplex[2].z, 
		simplex[3].x, simplex[3].y, simplex[3].z, 
		simplex[0].x, simplex[0].y, simplex[0].z, &type3, &which3 );
	if( type3 == TYPE_EDGE) which3 = emap3[which3];
	if( type3 == TYPE_POINT) which3 = pmap3[which3];
	// 23 13 12
	int emap4[3] = { 5, 4, 3 };
	int pmap4[3] = { 1, 2, 3 };
	double d4 = near_point_on_triangle(                    
		simplex[1].x, simplex[1].y, simplex[1].z, 
		simplex[2].x, simplex[2].y, simplex[2].z, 
		simplex[3].x, simplex[3].y, simplex[3].z, &type4, &which4 );
	if( type4 == TYPE_EDGE) which4 = emap4[which4];
	if( type4 == TYPE_POINT) which4 = pmap4[which4];

	int ninterior = 0;
	if( type1 == TYPE_INTERIOR ) ninterior++;
	if( type2 == TYPE_INTERIOR ) ninterior++;
	if( type3 == TYPE_INTERIOR ) ninterior++;
	if( type4 == TYPE_INTERIOR ) ninterior++;
	
	
	int npoints = 0;
	if( type1 == TYPE_POINT ) npoints++;
	if( type2 == TYPE_POINT ) npoints++;
	if( type3 == TYPE_POINT ) npoints++;
	if( type4 == TYPE_POINT ) npoints++;

//	printf("%d %d %d %d which %d %d %d %d\n", type1, type2, type3, type4, which1, which2, which3, which4 );

	int np[4] = {0,0,0,0};

	if( type1 == TYPE_POINT ) np[which1] += 1;
	if( type2 == TYPE_POINT ) np[which2] += 1;
	if( type3 == TYPE_POINT ) np[which3] += 1;
	if( type4 == TYPE_POINT ) np[which4] += 1;

	if( type1 == TYPE_EDGE && type2 == TYPE_EDGE && which1 == which2 )
	{	
		elements = 2;
		// EDGE 0, 1
	}
	else if( type2 == TYPE_EDGE && type3 == TYPE_EDGE && which2 == which3 )
	{	
		elements = 2;
		simplex[1] = simplex[3];
		// EDGE 0, 3
	}
	else if( type3 == TYPE_EDGE && type4 == TYPE_EDGE && which3 == which4 )
	{	
		elements = 2;
		simplex[0] = simplex[2];
		simplex[1] = simplex[3];
		// EDGE 2, 3
	}
	else if( type4 == TYPE_EDGE && type1 == TYPE_EDGE && which1 == which4 )
	{	
		elements = 2;
		simplex[0] = simplex[2];
		// EDGE 1, 2
	}
	else if( type1 == TYPE_EDGE && type3 == TYPE_EDGE && which1 == which3 )
	{	
		elements = 2;
		simplex[1] = simplex[2];
		// EDGE 0, 2
	}
	else if( type4 == TYPE_EDGE && type2 == TYPE_EDGE && which2 == which4 )
	{	
		elements = 2;
		simplex[0] = simplex[3];
		// EDGE 1, 3
	}
	else if( np[3] == 3 )	
	{
		elements = 1;
		simplex[0] = simplex[3];
	}
	else if( np[2] == 3  )	
	{
		elements = 1;
		simplex[0] = simplex[2];
	}
	else if( np[1] == 3 )	
	{
		elements = 1;
		simplex[0] = simplex[1];
	}
	else if( np[0] == 3  )	
	{
		elements = 1;
	}
	else
	{	// it's one of the interior points.
			
		elements = 3;
		
		double cur_dist = 0;
		if( d1 < d2 && d1 < d3 && d1 < d4 )
		{
			cur_dist = sqrt(d1);
		}
		else if( d2 < d3 && d2 < d4 )
		{
			simplex[2] = simplex[3];	
			cur_dist = sqrt(d2);
		}
		else if( d3 < d4  )
		{
			simplex[1] = simplex[3];	
			cur_dist = sqrt(d3);
		}
		else
		{
			simplex[0] = simplex[3];	
			cur_dist = sqrt(d4);
		}

//		printf("Cur dist: %le\n", cur_dist );
		
	}

	Point D;	
	Point v10 = simplex[0] - simplex[1];
	Point v20 = simplex[0] - simplex[2];
	Point neg1 = negative(simplex[1]);
	Point tempCross = cross(v10, neg1);
	switch( elements )
	{
		case 1:
			D = negative(simplex[0]);
			break;
		case 2:
    			D = negative(cross(tempCross, v10));
			break;
		case 3:
			D = cross( v10, v20 );
			break;
	}
		
	double k = dot( simplex[0], D);
	if( k > 0 )
		D = negative(D);

	Check new_check = {false, D};
	return new_check;
}

inline double Power2( double u )
{
        return u*u;
}

double Functions::near_point_on_triangle( double v1x, double v1y, double v1z,
				double v2x, double v2y, double v2z,
				double v3x, double v3y, double v3z, int *type, int *which )
{
/*
 * MATHEMATICA solution to:
  
p = u {v1x, v1y, v1z} + v {v2x, v2y, v2z} + (1 - u - v) {v3x, v3y, v3z}
Solve[D[p.p, u] == 0 && D[p.p, v] == 0, {u, v}]

interior points have u>0 && v >0 && u+v < 1

 * */
	double u = -((-4*(Power2(v2x - v3x) + Power2(v2y - v3y) + Power2(v2z - v3z))*
        ((v1x - v3x)*v3x + (v1y - v3y)*v3y + (v1z - v3z)*v3z) + 
       4*((v1x - v3x)*(v2x - v3x) + (v1y - v3y)*(v2y - v3y) + (v1z - v3z)*(v2z - v3z))*
        ((v2x - v3x)*v3x + (v2y - v3y)*v3y + (v2z - v3z)*v3z))/
     (4*Power2((v1x - v3x)*(v2x - v3x) + (v1y - v3y)*(v2y - v3y) + 
          (v1z - v3z)*(v2z - v3z)) - 
       4*(Power2(v1x - v3x) + Power2(v1y - v3y) + Power2(v1z - v3z))*
        (Power2(v2x - v3x) + Power2(v2y - v3y) + Power2(v2z - v3z))));

	double v = (-(u*v1x*v2x) - u*v1y*v2y - u*v1z*v2z + u*v1x*v3x - v2x*v3x + u*v2x*v3x + Power2(v3x) - 
     u*Power2(v3x) + u*v1y*v3y - v2y*v3y + u*v2y*v3y + Power2(v3y) - u*Power2(v3y) + 
     u*v1z*v3z - v2z*v3z + u*v2z*v3z + Power2(v3z) - u*Power2(v3z))/
   (Power2(v2x) + Power2(v2y) + Power2(v2z) - 2*v2x*v3x + Power2(v3x) - 2*v2y*v3y + 
     Power2(v3y) - 2*v2z*v3z + Power2(v3z));

	double dr[3] = { 
		u * v1x + v * v2x + (1-u-v) * v3x,
		u * v1y + v * v2y + (1-u-v) * v3y,
		u * v1z + v * v2z + (1-u-v) * v3z };
	
	double d= dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];		
	
	double r1 = v1x*v1x+v1y*v1y+v1z*v1z;
	double r2 = v2x*v2x+v2y*v2y+v2z*v2z;
	double r3 = v3x*v3x+v3y*v3y+v3z*v3z;

	*which = 0;	

	if( u > 0 && v > 0 && u + v < 1 )
	{
		*type = TYPE_INTERIOR;
	} 
	else
	{
		int cnd1 = u <= 0;
		int cnd2 = v <= 0;
		int cnd3 = u+v >= 1.0;

		int e_check1 = 0;
		int e_check2 = 0;
		int e_check3 = 0;


		if( cnd1 )
		{
			double v_edge_0 = (-(v2x*v3x) + Power2(v3x) - v2y*v3y + Power2(v3y) - v2z*v3z + Power2(v3z))/ (Power2(v2x - v3x) + Power2(v2y - v3y) + Power2(v2z - v3z));
			if( v_edge_0 >= 0 && v_edge_0 <= 1.0 )
				e_check1 = 1;
		}
		if( cnd2 )
		{
			double u_edge_0 = (-(v1x*v3x) + Power2(v3x) - v1y*v3y + Power2(v3y) - v1z*v3z + Power2(v3z))/ (Power2(v1x - v3x) + Power2(v1y - v3y) + Power2(v1z - v3z));
			if( u_edge_0 >= 0 && u_edge_0 <= 1.0 )
				e_check2 = 1;
		}
		if( cnd3 )
		{
			double v_opp = (Power2(v1x) + Power2(v1y) + Power2(v1z) - v1x*v2x - v1y*v2y - v1z*v2z)/ (Power2(v1x) + Power2(v1y) + Power2(v1z) - 2*v1x*v2x + Power2(v2x) - 2*v1y*v2y +      Power2(v2y) - 2*v1z*v2z + Power2(v2z));
			double u_opp = 1-v_opp;
			
			if( u_opp >= 0 && v_opp >= 0 )
				e_check3 = 1;
		}

		if( e_check1 )
		{
			*type = TYPE_EDGE;
			*which = 0; 
		}
		else if( e_check2 )
		{
			*type = TYPE_EDGE;
			*which = 1; 
		}
		else if( e_check3 )
		{
			*type = TYPE_EDGE;
			*which = 2; 
		}
		else
		{
			*type = TYPE_POINT;
			if( r1 < r2 && r1 < r3 )
				*which = 0;
			else if( r2 < r3 )
				*which = 1;
			else
				*which = 2;
								
		}
	}
	if( u < 0 || v < 0 || u+v >= 1.0 )
		d = 1e300;
//#define DEBUG
#ifdef DEBUG
/*
 * DEBUG
 *
 * */

	// line between 0,0 and 0,1
	double v_edge_0 = (-(v2x*v3x) + Power2(v3x) - v2y*v3y + Power2(v3y) - v2z*v3z + Power2(v3z))/ (Power2(v2x - v3x) + Power2(v2y - v3y) + Power2(v2z - v3z));
	double u_edge_0 = (-(v1x*v3x) + Power2(v3x) - v1y*v3y + Power2(v3y) - v1z*v3z + Power2(v3z))/ (Power2(v1x - v3x) + Power2(v1y - v3y) + Power2(v1z - v3z));
	double v_opp = (Power2(v1x) + Power2(v1y) + Power2(v1z) - v1x*v2x - v1y*v2y - v1z*v2z)/ (Power2(v1x) + Power2(v1y) + Power2(v1z) - 2*v1x*v2x + Power2(v2x) - 2*v1y*v2y +      Power2(v2y) - 2*v1z*v2z + Power2(v2z));
	double u_opp = 1-v_opp;
	
	double r_edge_0[3] = { 
		v_edge_0 * v2x + (1-v_edge_0) * v3x,
		v_edge_0 * v2y + (1-v_edge_0) * v3y,
		v_edge_0 * v2z + (1-v_edge_0) * v3z };
	double r_edge_1[3] = { 
		u_edge_0 * v1x + (1-u_edge_0) * v3x,
		u_edge_0 * v1y + (1-u_edge_0) * v3y,
		u_edge_0 * v1z + (1-u_edge_0) * v3z };
	double r_edge_2[3] = { 
		u_opp * v1x + v_opp * v2x,
		u_opp * v1y + v_opp * v2y,
		u_opp * v1z + v_opp * v2z };
	
	double r_int[3] = {
		u * v1x + v* v2x + (1-u-v)*v3x, 
		u * v1y + v* v2y + (1-u-v)*v3y, 
		u * v1z + v* v2z + (1-u-v)*v3z };
	
	double r2_e1 = r_edge_0[0]*r_edge_0[0] +  r_edge_0[1]*r_edge_0[1] + r_edge_0[2]*r_edge_0[2];
	double r2_e2 = r_edge_1[0]*r_edge_1[0] +  r_edge_1[1]*r_edge_1[1] + r_edge_1[2]*r_edge_1[2];
	double r2_e3 = r_edge_2[0]*r_edge_2[0] +  r_edge_2[1]*r_edge_2[1] + r_edge_2[2]*r_edge_2[2];

	if( v_edge_0 < 0 || v_edge_0 >= 1 )
		r2_e1 = 1e10;
	if( u_edge_0 < 0 || u_edge_0 >= 1 )
		r2_e2 = 1e10;
	if( u_opp < 0 || u_opp >= 1.0 || v_opp < 0 || v_opp >= 1.0 )
		r2_e3 = 1e10;

	double rs[7] = { r2_e1, r2_e2, r2_e3, r1, r2, r3, d };
	int types[7] = { TYPE_EDGE, TYPE_EDGE, TYPE_EDGE, TYPE_POINT, TYPE_POINT, TYPE_POINT, TYPE_INTERIOR };
	int whichs[7] = { 0, 1, 2, 0, 1, 2, 0 };
	int sorter[7] = {0,1,2,3,4,5,6};

	int done = 0;
	while( !done )
	{
		done = 1;
		for( int x = 0; x < 6; x++ )
		{
			if( rs[sorter[x]] > rs[sorter[x+1]] )
			{
				int t = sorter[x];
				sorter[x] = sorter[x+1];
				sorter[x+1] = t;
				done = 0;
			}
		}
	}

#ifdef DEBUG_PRINT
	printf("DBG EDGE 0 %le r2 %le\n", v_edge_0, r2_e1 );
	printf("DBG EDGE 1 %le r2 %le\n", u_edge_0, r2_e2 );
	printf("DBG EDGE 2 %le %le r2 %le\n", u_opp, v_opp, r2_e3 );
	printf("DBG PT 0 r2 %le\n", r1 );
	printf("DBG PT 1 r2 %le\n", r2 );
	printf("DBG PT 2 r2 %le\n", r3 );
	printf("DBG INT %le %le r2 %le\n", u, v, r_int[0]*r_int[0] +  r_int[1]*r_int[1] + r_int[2]*r_int[2] );

	printf("FUNCTION RETURNS TYPE %d WHICH %d r2 %le\n", *type, *which, d );
#endif
	static int iseq = 0;

	if( *type != types[sorter[0]] || *which != whichs[sorter[0]] )
	{
		printf("%d LOGICAL ERROR.\n", iseq);

		FILE *dbgFile = fopen("dbg.xyz","w");
	
		fprintf(dbgFile,"3002\n");
		fprintf(dbgFile,"DBG\n");
		fprintf(dbgFile,"He 0.0 0.0 0.0\n");
		fprintf(dbgFile,"F %lf %lf %lf\n", r_int[0], r_int[1], r_int[2] );
		// v = 0
		for( int iv = 0; iv < 1000; iv++ )
		{
			double u = -5 + 10 * iv / 1000.0;
			double v = 0;

	double r_plot[3] = {
		u * v1x + v* v2x + (1-u-v)*v3x, 
		u * v1y + v* v2y + (1-u-v)*v3y, 
		u * v1z + v* v2z + (1-u-v)*v3z };
			fprintf(dbgFile,"C %lf %lf %lf\n", r_plot[0], r_plot[1], r_plot[2] );			
		}
		// u = 0
		for( int iv = 0; iv < 1000; iv++ )
		{
			double v = -5 + 10 * iv / 1000.0;
			double u = 0;

	double r_plot[3] = {
		u * v1x + v* v2x + (1-u-v)*v3x, 
		u * v1y + v* v2y + (1-u-v)*v3y, 
		u * v1z + v* v2z + (1-u-v)*v3z };
			fprintf(dbgFile,"N %lf %lf %lf\n", r_plot[0], r_plot[1], r_plot[2] );			
		}
		// u+v = 1
		for( int iv = 0; iv < 1000; iv++ )
		{
			double u = -5 + 10 * iv / 1000.0;
			double v = 1-u;

	double r_plot[3] = {
		u * v1x + v* v2x + (1-u-v)*v3x, 
		u * v1y + v* v2y + (1-u-v)*v3y, 
		u * v1z + v* v2z + (1-u-v)*v3z };
			fprintf(dbgFile,"O %lf %lf %lf\n", r_plot[0], r_plot[1], r_plot[2] );			
		}


		exit(1);
	}
	d = rs[sorter[0]];
	iseq++;
#endif

	return d;
}


