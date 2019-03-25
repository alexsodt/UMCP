#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "interp.h"
#include "util.h"
#include "pdb.h"
#include "mutil.h"
#include <math.h>
#include "alignSet.h"
#include "pcomplex.h"
#include "mutil.h"
#include "input.h"

#define SMOOTH_TRIANGLES
// routine to write the current structure to a tachyon input file, including the option to linearly interpolate frames for smoothing.

static double *interp_data = NULL;
static int np_prev = -1; // number of points at previous request for interpolation
static double scale = -1;
static int frame_number = 0;
static int has_pbc = -1;
static int do_cyl = 0;
// Tachyon parameters
static int resolution = 1024;
static double color[3] = { 0, 1, 0 };
static double fixed_updir[3] = {0,0,1};
static double fixed_center[3] = {0,0,0};
static double fixed_view[3] = {0,0,1};
static int fixed_views_set = 0;

int surface::writeTachyon( const char *name,
				int grid_uv, // resolution, how many points we do on one side of a triangle 
				int nint, // number of interpolation frames
				double *r, // current coordinates				
				pcomplex **allComplexes,
				int ncomplex,
				parameterBlock *params,
				srd_integrator *srd,
				int face_center // put the origin on this face.
			)
{
	int overlay_mesh = params->tachyon_overlay_mesh;
	int np = nv;
	int nsites_total = nv+1;
	for( int c = 0; c < ncomplex; c++ )
		nsites_total += allComplexes[c]->nsites;

	int n_frames_to_write = nint;
	if( nint < 1 )
		n_frames_to_write = 1;

	if( n_frames_to_write > 1 )
	{
		if( np_prev != np )
		{
			free( interp_data );
			np_prev = -1;
			scale = 0;
		}

		if( !interp_data )
		{
			interp_data = (double *)malloc( sizeof(double) * 3 * nsites_total );			
			memcpy( interp_data, r, sizeof(double) * 3 * nsites_total  );
			return 0; // no frames currently and interpolation requested.
		} 
	}

	if( scale <= 0 )
	{
		double min = 1e10;
		double max = -1e10;
		
		for( int x = 0; x < nv; x++ )
		for( int xc = 0; xc < 3; xc++ )
		{
			if( r[3*x+xc] < min ) min = r[3*x+xc];
			if( r[3*x+xc] > max ) max = r[3*x+xc];
		}
	
		scale = 0.9 / (max-min);	

		printf("scale: %lf\n", scale );
	}



	if( has_pbc < 0 )	
	{
		// determine if it has PBC. If it doesn't, we will subtract off the center of mass

		has_pbc = 0;
		for( int x = 0; x < np && !has_pbc; x++ )
		{
			int val = theVertices[x].valence;

			for( int y = 0; y < val; y++ )
			for( int xc = 0; xc < 3; xc++ )
			{
				if( fabs(theVertices[x].edge_PBC[3*y+xc]) > 1e-7 )
				{
					if( !has_pbc && xc == 2 )
						has_pbc = 2;
					else if( xc != 2 )
						has_pbc = 1; 
				}
			}
		}
	}

	double alpha[3] = { r[3*nv+0], r[3*nv+1], r[3*nv+2] };

	double Lx = PBC_vec[0][0] * alpha[0];
	double Ly = PBC_vec[1][1] * alpha[1];
	double Lz = PBC_vec[2][2] * alpha[2];

	do_cyl = 0;
	int do_planar = 0;
	if( has_pbc )
		do_planar = 1;
	if( has_pbc == 2 )
	{
		do_planar = 0;
		do_cyl = 1;
	}
	double t = 0;

	double *use_frame = (double *)malloc( sizeof(double) * 3 * (nsites_total) );


	use_frame[3*nv+0] = alpha[0];	
	use_frame[3*nv+1] = alpha[1];	
	use_frame[3*nv+2] = alpha[2];	

//	use_frame[3*np+0] = 1.0;
//	use_frame[3*np+1] = 1.0;
//	use_frame[3*np+2] = 1.0;

	if( n_frames_to_write == 1 )
	{
		memcpy( use_frame, r, sizeof(double) * 3 * (nv+1) );
		
		int off = nv+1;	
		for( int c = 0; c < ncomplex; c++ )
		{
			memcpy( use_frame+off*3, allComplexes[c]->rall, sizeof(double) * 3 * allComplexes[c]->nsites );
			off += allComplexes[c]->nsites;
		}
	}

        int *tris = (int *)malloc( sizeof(int) * (grid_uv-1)*(grid_uv-1)*4 );
        int ntri = 0;

        int map[grid_uv*grid_uv];

        int iseq = 0;
        for( int iu = 0; iu < grid_uv; iu++ )
        for( int iv = 0; iv < grid_uv - iu; iv++, iseq++ )
                map[iu*grid_uv+iv] = iseq;
      
	for( int iu = 0; iu < grid_uv; iu++ )
        for( int iv = 0; iv < grid_uv - iu; iv++ )
        {
                if( iu < grid_uv-1 && iv < grid_uv-1 && iu+iv+1 < grid_uv )
                {
                        tris[ntri*4+0] = map[iu*grid_uv+iv]; // base. 
                        tris[ntri*4+1] = map[(iu+1)*grid_uv+iv]; // j
                        tris[ntri*4+2] = map[iu*grid_uv+iv+1];  // k
                        tris[ntri*4+3] = 1;
                        ntri++;
                }

                if( iu > 0 && iv < grid_uv-1 && iu+iv+1 < grid_uv )
                {
                        tris[ntri*4+0] = map[iu*grid_uv+iv];  // base
                        tris[ntri*4+1] = map[(iu)*grid_uv+iv+1];  //k
                        tris[ntri*4+2] = map[(iu-1)*grid_uv+iv+1];  // less u, greater v
                        tris[ntri*4+3] = 0;
                        ntri++;
                }
	}

	int n_frames_written = 0;
	
	if( interp_data && n_frames_to_write > 1 )
	{
		int off = nv+1;
		for( int c = 0; c < ncomplex; c++ )
		{
			for( int p = 0; p < allComplexes[c]->nsites; p++ )
			{
				double dr[3] = {
					interp_data[3*(off+p)+0] - allComplexes[c]->rall[3*p+0],
					interp_data[3*(off+p)+1] - allComplexes[c]->rall[3*p+1],
					interp_data[3*(off+p)+2] - allComplexes[c]->rall[3*p+2] };
				double put[3];
				MinImage3D( dr, PBC_vec, put );

				interp_data[3*(off+p)+0] = allComplexes[c]->rall[3*p+0] + dr[0];
				interp_data[3*(off+p)+1] = allComplexes[c]->rall[3*p+1] + dr[1];
				interp_data[3*(off+p)+2] = allComplexes[c]->rall[3*p+2] + dr[2];
			}

			off += allComplexes[c]->nsites;
		}
	}

	for( int i = 0; i < n_frames_to_write; i++, t += 1.0 /(double)(n_frames_to_write), n_frames_written++ )
	{
		if( n_frames_to_write > 1 )
		{
			for( int x = 0; x < (nsites_total)*3; x++ )
				use_frame[x] = (1-t) * interp_data[x] + t * r[x]; 
		}

		char fileName[strlen(name) + 16 ];
		char num5[256];
		print5(frame_number, num5 );
		sprintf(fileName, "%s%s.objs", name, num5 );		
		FILE *theFile = fopen(fileName,"w");

	double center[3] = { 0, 0, -2 };
	double view[3]   = { 0 ,0, 1 };
	double updir[3] = { 0, 1, 0 };
	int resolutiony = resolution;

	double zoom = 0.9;

	if( do_planar )
	{
	
			view[0] = 0;
			view[1] = -0.8;
			view[2] = -0.2;
			
			center[0] = -2*view[0];
			center[1] = -2*view[1];
			center[2] = -2*view[2];
			
			updir[0] = 0;
			updir[1] = -1;
			updir[2] = 0;
			resolutiony = ceil(resolution*1.0);
			if( resolutiony % 2 == 1 ) resolutiony++;		
	//	zoom = 2.0;
	}
	
	double bg[3] = {1,1,1};

	if( do_cyl )
	{
#define VIEW_Z
#ifdef VIEW_Z
		view[0] = 0;
		view[1] = 0;
		view[2] = 1;
		
		center[0] = 0;
		center[1] = 0;
		center[2] = -1;
		
		updir[0] = 0;
		updir[1] = 1;
		updir[2] = 0;
		zoom = 0.5;

		bg[0] = 1;
		bg[1] = 1;
		bg[2] = 1;

#else	
		view[0] = -0.5;
		view[1] = -0.5;
		view[2] = 2;
		
		center[0] = 0.25;
		center[1] = 0.25;
		center[2] = -1;
		
		updir[0] = 0;
		updir[1] = 1;
		updir[2] = 1;
		zoom = 0.5;
#endif
	}
	
	
	if( face_center >= 0 )	                 
	{
		double rcen[3];
		double ncen[3];       	

		evaluateRNRM( face_center, 1.0/3.0, 1.0/3.0, rcen, ncen, use_frame );


		if( !fixed_views_set )
		{
			center[0] = rcen[0] * scale;
			center[1] = rcen[1] * scale;
			center[2] = rcen[2] * scale;

			updir[0] = ncen[0];
			updir[1] = ncen[1];
			updir[2] = ncen[2];
		
			// designed for a cylinder in z but should work for a sphere too I think.
			double xydir[3] = {0, 0, 1};
	
			double t2[3];
			double t1[3];

			cross(  xydir, ncen, t2 );
			cross(  ncen, t2, t1 );

			view[0] = t2[0];
			view[1] = t2[1];
			view[2] = t2[2];
			// view is now side-on
			// rotate 45 back into the normal dir
			
			double origin[3] = {0,0,0};
			rotateArbitrary( view, t1, origin,   1, M_PI/6 );   
			rotateArbitrary( view, ncen, origin, 1, M_PI/4 );   
			// still side on

/*	
			double angle = 85;
			view[0] = cos(angle*M_PI/180.0) * view90[0] + sin(angle*M_PI/180.0) * ncen[0]; 
			view[1] = cos(angle*M_PI/180.0) * view90[1] + sin(angle*M_PI/180.0) * ncen[1]; 
			view[2] = cos(angle*M_PI/180.0) * view90[2] + sin(angle*M_PI/180.0) * ncen[2]; 
*/	
			normalize(view);
			normalize(updir);
			
			normalize(updir);
			
			fixed_view[0] = view[0];
			fixed_view[1] = view[1];
			fixed_view[2] = view[2];

			fixed_updir[0] = updir[0];
			fixed_updir[1] = updir[1];
			fixed_updir[2] = updir[2];

			fixed_center[0] = center[0] ;
			fixed_center[1] = center[1] ;
			fixed_center[2] = center[2] ;
			fixed_views_set = 1;
		}
		else
		{
			view[0] = fixed_view[0];
			view[1] = fixed_view[1];
			view[2] = fixed_view[2];

			updir[0] = fixed_updir[0];
			updir[1] = fixed_updir[1];
			updir[2] = fixed_updir[2];

			center[0] = fixed_center[0];
			center[1] = fixed_center[1];
			center[2] = fixed_center[2];
		}
			
		center[0] -= view[0];
		center[1] -= view[1];
		center[2] -= view[2];
	}

	fprintf(theFile,
"Begin_Scene\n"
"Resolution %d %d\n"
"Shader_Mode Medium\n"
"  Trans_VMD\n"
"End_Shader_Mode\n"
"Camera\n"
"  Projection Orthographic\n"
"  Zoom %lf\n"
"  Aspectratio 1\n"
"  Antialiasing 4\n"
"  Raydepth 10\n"
"  Center  %lf %lf %lf \n"
"  Viewdir %lf %lf %lf \n"
"  Updir   %lf %lf %lf \n"
"End_Camera\n"
"Directional_Light Direction %lf %lf %lf Color 1.0 1.0 1.0\n"
"Directional_Light Direction 0.1 -0.1 %lf Color 1.0 1.0 1.0\n"
"Directional_Light Direction -1 -2 %lf Color 0.2 0.2 0.2\n"
"Background %lf %lf %lf\n",
	resolution, resolutiony, zoom, center[0], center[1], center[2], view[0], view[1], view[2], updir[0], updir[1], updir[2], -updir[0], -updir[1], -updir[2], (do_planar ? -1.0 : 1.0), (do_planar ? -0.5 : 0.5), bg[0], bg[1], bg[2] );
	
		int num = grid_uv * (grid_uv+1)/2;

		int dx_min_lim=0, dx_max_lim=0;
		int dy_min_lim=0, dy_max_lim=0;
		int dz_min_lim=0, dz_max_lim=0;
		if( do_planar && 0)
		{
			dx_min_lim=-1;
			dx_max_lim=1;
			dy_min_lim=-1;
			dy_max_lim=1;
		}

		if( do_cyl )
		{
			dz_min_lim = -1;
			dz_max_lim = 1;
		}

		for( int dx = dx_min_lim; dx <= dx_max_lim; dx++ )
		for( int dy = dy_min_lim; dy <= dy_max_lim; dy++ )
		for( int dz = dz_min_lim; dz <= dz_max_lim; dz++ )
		{
			for( int f = 0; f < nf_faces + nf_irr_faces; f++ )
			{	
				double Rvals[3*num];
				double Nvals[3*num];
		
				iseq = 0;

/*
				int tri;
				if( f < nf_faces )
					tri = theFormulas[f*nf_g_q_p].tri;
				else
					tri = theIrregularFormulas[(f-nf_faces)*nf_irr_pts].tri;

				double fl = theTriangles[tri].f_lipids;

				double color[3] = { 1.0, 1.0, 1.0 };
				double color_thresh = 0.2;
				if( fl < 1 )
				{
					color[1] -= (1-fl)/color_thresh;
					color[2] -= (1-fl)/color_thresh;
					if( fl < 1-color_thresh )
					{
						color[1] = 0;
						color[2] = 0;
					}
				}
				else
				{
					color[0] -= (fl-1)/color_thresh;
					color[1] -= (fl-1)/color_thresh;
					if( fl > 1 + color_thresh )
					{
						color[0] = 0;
						color[1] = 0;
					}
				}
*/
		                for( int iu = 0; iu < grid_uv; iu++ )
		                for( int iv = 0; iv < grid_uv - iu; iv++, iseq++ )
		                {
		                        // evaluate R, N, c.
		
		                        double u = iu/(double)(grid_uv-1);
		                        double v = iv/(double)(grid_uv-1);

					if( f >= nf_faces && iu == 0 && iv == 0 )
					{
		                        	evaluateRNRM( f, u+1e-6, v+1e-6, Rvals+iseq*3, Nvals+iseq*3, use_frame );
								/*	
						int v = theIrregularFormulas[(f-nf_faces)*nf_irr_pts+0].vertex;
						int e = theIrregularFormulas[(f-nf_faces)*nf_irr_pts+0].edge;
				
						int val = theVertices[v].valence;	
						double lr[3] = { 0,0,0};

						double *ti = theVertices[v].r;

						lr[0] = 0.5 * use_frame[v*3+0]; 
						lr[1] = 0.5 * use_frame[v*3+1]; 
						lr[2] = 0.5 * use_frame[v*3+2]; 

						double w = 1.0 / (val * 2);

						for( int e = 0; e < val; e++ )
						{
							int j = theVertices[v].edges[e];

							double dr[3] = { use_frame[j*3+0] - ti[0], use_frame[j*3+1] - ti[1], use_frame[j*3+2] - ti[2] };

							double put[3];	
							MinImage3D( dr, PBC_vec, put );
							dr[0] += ti[0];
							dr[1] += ti[1];
							dr[2] += ti[2];
				
							lr[0] += dr[0]*w;
							lr[1] += dr[1]*w;
							lr[2] += dr[2]*w;
						}

						Rvals[iseq*3+0] = lr[0];
						Rvals[iseq*3+1] = lr[1];
						Rvals[iseq*3+2] = lr[2];*/
					} 	
					else
			
		                        	evaluateRNRM( f, u, v, Rvals+iseq*3, Nvals+iseq*3, use_frame );
		                }

				for( int x1 = 0; x1 < ntri; x1++ )
				{       
				        int ir1 = tris[x1*4+0];
				        int ir2 = tris[x1*4+1];
				        int ir3 = tris[x1*4+2];
					double nfac = 1.0;
					if( do_planar )
						nfac = 1;
#ifdef SMOOTH_TRIANGLES
					fprintf(theFile, "STri\n");
#else
					fprintf(theFile, "Tri\n");
#endif
					fprintf(theFile, "V0 %lf %lf %lf\n", Rvals[ir1*3+0]*scale + dx*Lx*scale, Rvals[ir1*3+1]*scale+ dy*Ly*scale, Rvals[ir1*3+2]*scale+ dz*Lz*scale );
					fprintf(theFile, "V1 %lf %lf %lf\n", Rvals[ir2*3+0]*scale + dx*Lx*scale, Rvals[ir2*3+1]*scale+ dy*Ly*scale, Rvals[ir2*3+2]*scale+ dz*Lz*scale );
					fprintf(theFile, "V2 %lf %lf %lf\n", Rvals[ir3*3+0]*scale + dx*Lx*scale, Rvals[ir3*3+1]*scale+ dy*Ly*scale, Rvals[ir3*3+2]*scale+ dz*Lz*scale );
#ifdef SMOOTH_TRIANGLES
					fprintf(theFile, "N0 %lf %lf %lf\n", nfac*Nvals[ir1*3+0]*scale, nfac*Nvals[ir1*3+1]*scale, nfac*Nvals[ir1*3+2]*scale );
					fprintf(theFile, "N1 %lf %lf %lf\n", nfac*Nvals[ir2*3+0]*scale, nfac*Nvals[ir2*3+1]*scale, nfac*Nvals[ir2*3+2]*scale );
					fprintf(theFile, "N2 %lf %lf %lf\n", nfac*Nvals[ir3*3+0]*scale, nfac*Nvals[ir3*3+1]*scale, nfac*Nvals[ir3*3+2]*scale );
#endif
					fprintf(theFile, "Texture\n"
							"Ambient 0 Diffuse 0.45 Specular 0.2 Opacity 1.0\n"
							 "Phong Metal 0.5 Phong_size 40 Color %lf %lf %lf TexFunc 0\n", color[0], color[1], color[2] );
				}
			} 	
		
			if( overlay_mesh )
			{
				for( int i = 0; i < nv; i++ )
				{
					fprintf(theFile, "Sphere\n");
					fprintf(theFile, "Center %lf %lf %lf\n", alpha[0]*use_frame[3*i+0]*scale +dx*Lx*scale, alpha[1]*use_frame[3*i+1]*scale +dy*Ly*scale, alpha[2]*use_frame[3*i+2]*scale + dz * Lz*scale );
					fprintf(theFile, "Rad %lf\n", 9.0 * scale );
					fprintf(theFile, "Texture\n"
							"Ambient 0.1 Diffuse 0.4 Specular 0.0 Opacity 1.0\n"
							 "Phong Metal 0.5 Phong_size 40 Color 1.0 1.0 1.0 TexFunc 0\n");
					for( int e = 0; e < theVertices[i].valence; e++ )
					{
						double sgn=1;
						int j = theVertices[i].edges[e];	
						double rj[3] = { 
							(use_frame[3*j+0]+sgn*theVertices[i].edge_PBC[3*e+0] * PBC_vec[0][0]+sgn*theVertices[i].edge_PBC[3*e+1] * PBC_vec[1][0]+sgn*theVertices[i].edge_PBC[3*e+2] * PBC_vec[2][0] )*scale + dx*Lx*scale,
							(use_frame[3*j+1]+sgn*theVertices[i].edge_PBC[3*e+1] * PBC_vec[0][1]+sgn*theVertices[i].edge_PBC[3*e+1] * PBC_vec[1][1]+sgn*theVertices[i].edge_PBC[3*e+2] * PBC_vec[2][1] )*scale + dy*Ly*scale, 
							(use_frame[3*j+2]+sgn*theVertices[i].edge_PBC[3*e+2] * PBC_vec[0][2]+sgn*theVertices[i].edge_PBC[3*e+1] * PBC_vec[1][2]+sgn*theVertices[i].edge_PBC[3*e+2] * PBC_vec[2][2] )*scale + dz*Lz*scale };
						fprintf(theFile, "FCylinder\n");
						fprintf(theFile, "Base %lf %lf %lf\n",  alpha[0]*(use_frame[3*i+0]*scale +dx*Lx*scale), alpha[1]*(use_frame[3*i+1]*scale +dy*Ly*scale), alpha[2]*(use_frame[3*i+2]*scale + dz * Lz *scale) );
						fprintf(theFile, "Apex %lf %lf %lf\n", alpha[0]*rj[0], alpha[1]*rj[1], alpha[2]*rj[2] );
						fprintf(theFile, "Rad %lf\n", 5.0 * scale ); // six angstroms
						fprintf(theFile, "Texture\n"
							"Ambient 0.1 Diffuse 0.4 Specular 0.0 Opacity 1\n"
							 "Phong Metal 0.5 Phong_size 40 Color 0 0 1.0 TexFunc 0\n");
					}
				}

			}	
			int off = nv+1;
			for( int c = 0; c < ncomplex; c++ )
			{
				double use_p[3*allComplexes[c]->nsites];
				for( int p = 0; p < allComplexes[c]->nsites; p++ )
				{

					double p_rad = allComplexes[c]->sigma[p];
					double pcolor[3] = { 1, 0, 1 };

					if( allComplexes[c]->att_sigma[p] > 0 )
					{
						p_rad = allComplexes[c]->att_sigma[p];
						pcolor[0] = 0;
						pcolor[1] = 0;
						pcolor[2] = 1;
					}

						// wrap PBC
	
						double pbc_w[3] = { use_frame[3*(off+p)+0], 
								    use_frame[3*(off+p)+1],						
								    use_frame[3*(off+p)+2] };

//						while( pbc_w[0] < 0  ) pbc_w[0] += Lx;						
//						while( pbc_w[0] >  Lx ) pbc_w[0] -= Lx;						
						
//						while( pbc_w[1] < 0 ) pbc_w[1] += Ly;						
//						while( pbc_w[1] >  Ly ) pbc_w[1] -= Ly;						
						
						while( pbc_w[2] < 0 ) pbc_w[2] += Lz;						
						while( pbc_w[2] >  Lz ) pbc_w[2] -= Lz;						

						use_p[3*p+0] = pbc_w[0];
						use_p[3*p+1] = pbc_w[1];
						use_p[3*p+2] = pbc_w[2];

						fprintf(theFile, "Sphere\n");
						fprintf(theFile, "Center %lf %lf %lf\n",  
	pbc_w[0]*scale +dx*Lx*scale, pbc_w[1]*scale +dy*Ly*scale, pbc_w[2]*scale + dz * Lz*scale );
						fprintf(theFile, "Rad %lf\n", p_rad * scale );
						fprintf(theFile, "Texture\n"
							"Ambient 0.1 Diffuse 0.4 Specular 0.0 Opacity 1.0\n"
							 "Phong Metal 0.5 Phong_size 40 Color %lf %lf %lf TexFunc 0\n", pcolor[0], pcolor[1], pcolor[2] );
				}

				int nbonds = allComplexes[c]->getNBonds();
				int bonds[2*nbonds];
				allComplexes[c]->putBonds( bonds );
	
				for( int tb = 0; tb < nbonds; tb++ )
				{
					double p1[3] = {
					use_p[3*bonds[2*tb]+0],
					use_p[3*bonds[2*tb]+1],
					use_p[3*bonds[2*tb]+2] };
					double p2[3] = {
					use_p[3*bonds[2*tb+1]+0],
					use_p[3*bonds[2*tb+1]+1],
					use_p[3*bonds[2*tb+1]+2] };

					while( p2[0] - p1[0] > Lx /2 ) p2[0] -= Lx;
					while( p2[0] - p1[0]< -Lx /2 ) p2[0] += Lx;
					while( p2[1] - p1[1] > Ly /2 ) p2[1] -= Ly;
					while( p2[1] - p1[1]< -Ly /2 ) p2[1] += Ly;
					while( p2[2] - p1[2] > Lz /2 ) p2[2] -= Lz;
					while( p2[2] - p1[2]< -Lz /2 ) p2[2] += Lz;
						fprintf(theFile, "FCylinder\n");
						fprintf(theFile, "Base %lf %lf %lf\n",  p1[0]*scale +dx*Lx*scale, p1[1]*scale +dy*Ly*scale, p1[2]*scale + dz * Lz *scale );
						fprintf(theFile, "Apex %lf %lf %lf\n",  p2[0]*scale +dx*Lx*scale, p2[1]*scale +dy*Ly*scale, p2[2]*scale + dz * Lz *scale );
						fprintf(theFile, "Rad %lf\n", 10.0 * scale ); // one nm for now
						fprintf(theFile, "Texture\n"
							"Ambient 0.1 Diffuse 0.4 Specular 0.0 Opacity 1\n"
							 "Phong Metal 0.5 Phong_size 40 Color 1.0 1.0 1.0 TexFunc 0\n");
				}

				off += allComplexes[c]->nsites;
			}
		}

		

		fprintf(theFile, "End_Scene\n");

		fclose(theFile );
	
		frame_number++;
	}

	if( interp_data && n_frames_to_write > 1 )
	{	
		memcpy( interp_data, r, sizeof(double) * 3 * (nsites_total) );
	}
	free(use_frame);
	free(tris);

	return n_frames_written;
}

#if 0
static int nhel = 0;
static int *helices = NULL;
static double *rvals = NULL;
static double *aligned_rvals = NULL;
static double centers[9];
static double normals[9];
static int data_per_helix = 6;

void loadBar( const char *proteinFileName, int c1, char chain1, int c2, char chain2, int c3, char chain3 )
{
	int nhel_space = 1;
	
	helices = (int *)malloc( sizeof(int) * nhel_space * data_per_helix );
	rvals = (double *)malloc( sizeof(double) * 6 * nhel_space );
	nhel = 0;

	FILE *theFile = fopen(proteinFileName,"r");

	if( !theFile )
	{
		printf("Couldn't open file '%s'.\n", proteinFileName );
		exit(1);
	}

	char buffer[4096];

	memset( centers, 0, sizeof(double) * 9 );
	while( !feof(theFile) ) 
	{	
		getLine( theFile, buffer );
		if( !strncasecmp( buffer, "HELIX", 5 ) )
		{
			int num;
			char resName1[256], resName2[256];;
			int r1,r2;
			char chain1,chain2;
			int nr = sscanf( buffer, "HELIX %d %d %s %c %d %s %c %d",
				&num, &num, resName1, &chain1, &r1, resName2, &chain2, &r2 );
			printf("nr: %d buffer %s\n", nr, buffer );
			if( nr == 8 )
			{
				if( nhel == nhel_space )
				{	
					nhel_space *= 2;
					helices = (int *)realloc( helices, sizeof(int) * nhel_space * data_per_helix );
					rvals = (double *)realloc( rvals, sizeof(double) * 6 * nhel_space );
				}
	
				helices[nhel*data_per_helix+0] = r1;
				helices[nhel*data_per_helix+1] = r2;
				helices[nhel*data_per_helix+2] = chain1;
				helices[nhel*data_per_helix+3] = chain2;
	
				nhel++;
			}
		}
		else if( !strncasecmp( buffer, "ATOM", 4 ) )
		{
			struct atom_rec at;
			readATOM( buffer, &at );

			if( strcasecmp( at.atname, "CA" ) ) continue;

			if( at.res == c1 && at.chain == chain1 )
			{
				printf("Setting center 1 to %lf %lf %lf\n", at.x, at.y, at.z );
				centers[0] = at.x;
				centers[1] = at.y;
				centers[2] = at.z;
			}
			else if( at.res == c2 && at.chain == chain2 )
			{
				printf("Setting center 2 to %lf %lf %lf\n", at.x, at.y, at.z );
				centers[3] = at.x;
				centers[4] = at.y;
				centers[5] = at.z;
			}
			else if( at.res == c3 && at.chain == chain3 )
			{
				printf("Setting center 3 to %lf %lf %lf\n", at.x, at.y, at.z );
				centers[6] = at.x;
				centers[7] = at.y;
				centers[8] = at.z;
			}
			for( int h = 0; h < nhel; h++ )
			{
				if( helices[h*data_per_helix+0] == at.res && helices[h*data_per_helix+2] == at.chain )
				{
					rvals[h*6+0] = at.x;
					rvals[h*6+1] = at.y;
					rvals[h*6+2] = at.z;
					helices[h*data_per_helix+0] = -1;
				}
				else if( helices[h*data_per_helix+1] == at.res && helices[h*data_per_helix+3] == at.chain )
				{
					rvals[h*6+3] = at.x;
					rvals[h*6+4] = at.y;
					rvals[h*6+5] = at.z;
					helices[h*data_per_helix+1] = -1;
				}
			}
			at.zap();						
		}
	}

	printf("Found %d helices.\n", nhel );

	int ok = 1;
	int ngot = 0;
	for( int h = 0; h < nhel; h++ )
	{
		if( helices[h*data_per_helix+0] != -1 || helices[h*data_per_helix+1] != -1 )
			ok = 0;
		else ngot++;
	}
	printf("Accounted for %d/%d helices.\n", ngot, nhel );

	// set the normals used for the orientation.

	double dr1[3] = { centers[0] - centers[3], centers[1] - centers[4], centers[2] - centers[5] };
	double dr2[3] = { centers[6] - centers[3], centers[7] - centers[4], centers[8] - centers[5] };

	normalize(dr1);
	normalize(dr2);
	
	double off[3];
		
	// arbitrary, the "y" direction
	cross( dr1, dr2, off );

	double averagex[3] = { dr1[0] - dr2[0], dr1[1] - dr2[1], dr1[2] - dr2[2] };
	normalize(averagex);

	// the center norm:
	cross( averagex, off, normals+3 );
	normalize(normals+3);

	// center 1 normal:
	cross( dr1, off, normals );
	normalize(normals);
	
	// negate dr2 by switching the cross product order:
	cross( off, dr2, normals+6 );
	normalize(normals);

	// partition helices between the centers.

	for( int h = 0; h < nhel; h++ )
	for( int comp = 0; comp < 2; comp++ )
	{
		double *com = rvals+h*6+comp*3;	
		double mind = 1e10;

		for( int c = 0; c < 3; c++ )
		{
			double dr[3] = { com[0] - centers[c*3+0], com[1] - centers[c*3+1], com[2] - centers[c*3+2] };

			double r=sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

			if( r < mind )
			{
				helices[data_per_helix*h+4+comp] = c;
				mind = r;
			}
		}
		printf("assigning helix %d component %d to center %d.\n", h, comp, helices[h*data_per_helix+4+comp] );
	}

	aligned_rvals = (double *)malloc( sizeof(double) * 6 * nhel );

}

void drawBarDomain( FILE *theFile, double scale, double *r, double *n, double PBC_vec[3][3], int dx, int dy, int dz)
{
	double Lx = PBC_vec[0][0];
	double Ly = PBC_vec[1][1];
	double Lz = PBC_vec[2][2];

	if( nhel == 0 )
		loadBar( "4i1q.pdb", 156, 'A', 73, 'A', 156, 'B' ); 

	// align dr and normals to center dr and normals.

	int align_2[4] = { 0, 1, 2, 3 };
	double ralign1[12+nhel*6];
	int align_1[4+nhel*2];
	align_1[0]=0;
	align_1[1]=1;
	align_1[2]=2;
	align_1[3]=3;
	// center 0
	for( int center = 0; center < 3; center++ )
	{
		int nalign_1 = 3;
		int nalign_2 = 3;
		double ralign2[12] = {	
			// the center.
			r[center*3+0], r[center*3+1], r[center*3+2],
			// the center + the normal
			r[center*3+0]+n[center*3+0], r[center*3+1]+n[center*3+1], r[center*3+2]+n[center*3+2],
			// the middle bead.
			r[3], r[4], r[5],
			0,0,0
		};
		
		ralign1[0] = centers[center*3+0];		
		ralign1[1] = centers[center*3+1];		
		ralign1[2] = centers[center*3+2];		
		
		ralign1[3] = centers[center*3+0]+normals[center*3+0];		
		ralign1[4] = centers[center*3+1]+normals[center*3+1];		
		ralign1[5] = centers[center*3+2]+normals[center*3+2];		
		
		ralign1[6] = centers[3];		
		ralign1[7] = centers[4];		
		ralign1[8] = centers[5];		
			

		if( center == 1 )
		{
			ralign2[6] = r[0];
			ralign2[7] = r[1];
			ralign2[8] = r[2];
			
			ralign2[9]  = r[6];
			ralign2[10] = r[7];
			ralign2[11] = r[8];
		
			ralign1[6] = centers[0];		
			ralign1[7] = centers[1];		
			ralign1[8] = centers[2];		
			
			ralign1[9]  = centers[6];		
			ralign1[10] = centers[7];		
			ralign1[11] = centers[8];

			nalign_2 = 4;
			nalign_1 = 4;		

		}
		
		int align_inst[3][3] = {
			2,0,1,
			0,2,3,
			2,0,1
		};
			
		for( int alx = 0; alx < 2; alx++ )
		{
			int c_al_to = align_inst[center][0];
			int c_al   = align_inst[center][1+alx];

			double put[3];
			double dr[3] = { 
				ralign2[c_al*3+0] - ralign2[c_al_to*3+0], 
				ralign2[c_al*3+1] - ralign2[c_al_to*3+1], 
				ralign2[c_al*3+2] - ralign2[c_al_to*3+2] }; 

			MinImage3D( dr, PBC_vec,put);

			ralign2[c_al*3+0] = ralign2[c_al_to*3+0] + dr[0];
			ralign2[c_al*3+1] = ralign2[c_al_to*3+1] + dr[1];
			ralign2[c_al*3+2] = ralign2[c_al_to*3+2] + dr[2];
		}
	
		int nalign_base = nalign_2;
	
		int nalign = nalign_1;

		for( int h = 0; h < nhel; h++ )
		for( int comp = 0; comp < 2; comp++ )
		{
			if( helices[h*data_per_helix+4] == center )
			{
				ralign1[nalign_1*3+0] = rvals[h*6+comp*3+0];
				ralign1[nalign_1*3+1] = rvals[h*6+comp*3+1];
				ralign1[nalign_1*3+2] = rvals[h*6+comp*3+2];
				align_1[nalign_1] = nalign_1;
				nalign_1++;
			}
		}

#if 0 
		printf("align p %lf %lf %lf and %lf %lf %lf and %lf %lf %lf\n",
			ralign1[0], ralign1[1], ralign1[2],
			ralign1[3], ralign1[4], ralign1[5],
			ralign1[6], ralign1[7], ralign1[8] );
		
		printf("align p2 %lf %lf %lf and %lf %lf %lf and %lf %lf %lf\n",
			ralign2[0], ralign2[1], ralign2[2],
			ralign2[3], ralign2[4], ralign2[5],
			ralign2[6], ralign2[7], ralign2[8] );
#endif
		alignStructuresOnAtomSet( ralign2,  align_2, ralign1, align_1, nalign_2, nalign_1 );

		nalign_1 = nalign_base;

		for( int h = 0; h < nhel; h++ )
		for( int comp = 0; comp < 2; comp++ )
		{
			if( helices[h*data_per_helix+4] == center )
			{
				aligned_rvals[h*6+comp*3+0]=ralign1[nalign_1*3+0];
				aligned_rvals[h*6+comp*3+1]=ralign1[nalign_1*3+1];
				aligned_rvals[h*6+comp*3+2]=ralign1[nalign_1*3+2];

				nalign_1++;
			}
		}
	}
	
	for( int h = 0; h < nhel; h++ )
	{
		double end1[3] = { aligned_rvals[h*6+0], aligned_rvals[h*6+1], aligned_rvals[h*6+2] };
		double end2[3] = { aligned_rvals[h*6+3], aligned_rvals[h*6+4], aligned_rvals[h*6+5] };
	
		fprintf(theFile, "FCylinder\n");
		fprintf(theFile, "Base %lf %lf %lf\n", end1[0]*scale+dx*Lx*scale, end1[1]*scale+dy*Ly*scale, end1[2]*scale +dz*Lz*scale);
		fprintf(theFile, "Apex %lf %lf %lf\n", end2[0]*scale+dx*Lx*scale, end2[1]*scale+dy*Ly*scale, end2[2]*scale +dz*Lz*scale);
		fprintf(theFile, "Rad %lf\n", 6.0 * scale ); // six angstroms
		fprintf(theFile, "Texture\n"
				"Ambient 0 Diffuse 0.45 Specular 0.4 Opacity 1\n"
				 "Phong Metal 0.5 Phong_size 40 Color 1 0 0 TexFunc 0\n");
		
		fprintf(theFile, "Ring\n");
		fprintf(theFile, "Center %lf %lf %lf\n", end1[0]*scale+dx*Lx*scale, end1[1]*scale+dy*Ly*scale, end1[2]*scale + dz*Lz*scale );
		double dr[3] = { end2[0]-end1[0], end2[1]-end1[1], end2[2]-end1[2] };
		normalize(dr);
		fprintf(theFile, "Normal %lf %lf %lf\n", dr[0],dr[1],dr[2] );
		fprintf(theFile, "Inner 0.0 Outer %lf\n", 6.0 * scale ); // six angstroms
		fprintf(theFile, "Texture\n"
				"Ambient 0 Diffuse 0.45 Specular 0.4 Opacity 1\n"
				 "Phong Metal 0.5 Phong_size 40 Color 1 0 0 TexFunc 0\n");
		
		fprintf(theFile, "Ring\n");
		fprintf(theFile, "Center %lf %lf %lf\n", end2[0]*scale+dx*Lx*scale, end2[1]*scale+dy*Ly*scale, end2[2]*scale + dz*Lz*scale );
		fprintf(theFile, "Normal %lf %lf %lf\n", dr[0],dr[1],dr[2] );
		fprintf(theFile, "Inner 0.0 Outer %lf\n", 6.0 * scale ); // six angstroms
		fprintf(theFile, "Texture\n"
				"Ambient 0 Diffuse 0.45 Specular 0.4 Opacity 1\n"
				 "Phong Metal 0.5 Phong_size 40 Color 1 0 0 TexFunc 0\n");
		
	}

}
#endif











