#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "interp.h"
#include "mutil.h"
#include "l-bfgs.h"
#include "pcomplex.h"
#include "mutil.h"
	
static int plim = 1;

void surface::writeSurface( const char *fileName)
{
	int total_valence = 0;

	for( int i = 0; i < nv; i++ )
	{
		total_valence += theVertices[i].valence;
		if( theVertices[i].broken[theVertices[i].valence-1] )
			total_valence -= 1;
	}

	int npr = 4;
	int npts_print = nv + total_valence * 4;

	FILE *surfaceXYZ = fopen(fileName,"w");

	fprintf(surfaceXYZ, "%d\n", npts_print );
	fprintf(surfaceXYZ, "simple\n");

	for( int i = 0; i < nv; i++ )
	{
		int val = theVertices[i].valence;
			
		fprintf(surfaceXYZ, "N %lf %lf %lf\n", theVertices[i].r[0], theVertices[i].r[1], theVertices[i].r[2] );

	}
	
	for( int i = 0; i < nv; i++ )
	{
		int val = theVertices[i].valence;
			
		for( int e = 0; e < val; e++ )
		{
			int j = theVertices[i].edges[e];

			double dr[3] = { 
				theVertices[j].r[0] - theVertices[i].r[0] + theVertices[i].edge_PBC[3*e+0]*PBC_vec[0][0] + theVertices[i].edge_PBC[3*e+1] * PBC_vec[1][0]+ theVertices[i].edge_PBC[3*e+2] * PBC_vec[2][0],
				theVertices[j].r[1] - theVertices[i].r[1] + theVertices[i].edge_PBC[3*e+0]*PBC_vec[0][1] + theVertices[i].edge_PBC[3*e+1] * PBC_vec[1][1]+ theVertices[i].edge_PBC[3*e+2] * PBC_vec[2][1],
				theVertices[j].r[2] - theVertices[i].r[2] + theVertices[i].edge_PBC[3*e+0]*PBC_vec[0][2] + theVertices[i].edge_PBC[3*e+1] * PBC_vec[1][2]+ theVertices[i].edge_PBC[3*e+2] * PBC_vec[2][2] };
		
			for( int it = 0; it < npr; it++ )
			{ 
				double t = (it+0.5)/(double)npr;
				double tr[3] = { 
					theVertices[i].r[0] + dr[0] * t,
					theVertices[i].r[1] + dr[1] * t,
					theVertices[i].r[2] + dr[2] * t };
				
				fprintf(surfaceXYZ, "O %lf %lf %lf # %d %d\n", tr[0], tr[1], tr[2], i , j );
			}
		}
	}

	fclose(surfaceXYZ);
	
	void debugDSduv( double *b );


//	debugDSduv( theVertices[0].btriangle );

}








static	int narray[3] = { 1, 4, 16 };
static	double use_array[3][16][2] =
	{
		{	{ 1.0/3.0, 1.0/3.0 } },

		{	
			{1.0/3.0,               1.0/3.0},

			{1.0/6.0,               1.0/6.0},
			{ 0.5 + 1.0/6.0,        1.0/6.0 },
			{ 1.0/6.0,        0.5 + 1.0/6.0 }
		}
	};

static	int which = 1;



void surface::writeToSurfaceGen( surface *lowerSurface, double thickness, const char *job_name)
{

const char *header2 =
"%%EndComments\n"
"\n"
"/arrowdict 14 dict def \n"
"arrowdict begin\n"
" /mtrx matrix def \n"
"end\n"
"\n"
"/arrow\n"
" { arrowdict begin\n"
"   /headlength exch def \n"
"   /halfheadthickness exch 2 div def \n"
"   /halfthickness exch 2 div def \n"
"   /tipy exch def /tipx exch def \n"
"   /taily exch def /tailx exch def \n"
"\n"
"   /dx tipx tailx sub def \n"
"   /dy tipy taily sub def \n"
"   /arrowlength dx dx mul dy dy mul add sqrt def \n"
"   /angle dy dx atan def \n"
"   /base arrowlength headlength sub def \n"
"   /savematrix mtrx currentmatrix def \n"
" \n"
"   tailx taily translate\n"
"   angle rotate\n"
"   0 halfthickness neg moveto\n"
"   base halfthickness neg lineto\n"
"   base halfheadthickness neg lineto\n"
"   arrowlength 0 lineto\n"
"   base halfheadthickness lineto\n"
"   base halfthickness lineto\n"
"   0 halfthickness lineto\n"
"   closepath\n"
"   savematrix setmatrix\n"
"  end\n"
"} def\n";

const char *ps_header = 
"%%!PS-Adobe-3.0 EPSF-3.0\n"
"%%%%BoundingBox: 0 0 %lf %lf\n"
"%%%%Orientation: Portrait\n"
"%%%%Pages: 1\n"
"%%%%EndComments\n"
"%%%%Page: 1 1\n";

	int nim = 7;
	int do_hex_BC = 1;
	int nleaflets = 1;
	
	if( !do_hex_BC )
		nim = 1;

	const char *thename = "test";

	double avL = 1.0;
	double target_h = thickness;
			
	char fileName[256];
		
	sprintf(fileName, "%sHeight.ps", job_name  );

	FILE *psFile = fopen(fileName,"w");
	
	double width = 800;
	double height = 800;
	
	fprintf(psFile, ps_header, width, height );
	fprintf(psFile, "%s", header2 );
	fprintf(psFile, "%lf %lf translate\n", width/2, height/2 );
	
	double scale = 3;
	for( int ix = -1; ix <= 1; ix++ )
	for( int iy = -1; iy <= 1; iy++ )
	{
		double shift[3] = { ix * PBC_vec[0][0] + iy * PBC_vec[1][0],
				    ix * PBC_vec[0][1] + iy * PBC_vec[1][1], 0 }; 
	for( int t = 0; t < nt; t++ )
	{
		int i = theTriangles[t].ids[0];
		int j = theTriangles[t].ids[1];
		int k = theTriangles[t].ids[2];

		vertex *theVerticesB = lowerSurface->theVertices;
		double *pbc1 = theTriangles[t].pbc1;
		double *pbc2 = theTriangles[t].pbc2;
		double alpha = 1.0;
		double p1[3] = { theVertices[i].r[0], theVertices[i].r[1], theVertices[i].r[2] };

		double p2[3] = {theVertices[j].r[0]	+pbc1[0] * PBC_vec[0][0]*alpha + pbc1[1] * PBC_vec[1][0],
		theVertices[j].r[1]	+pbc1[0] * PBC_vec[0][1]*alpha + pbc1[1] * PBC_vec[1][1],
		theVertices[j].r[2]	+pbc1[0] * PBC_vec[0][2]*alpha + pbc1[1] * PBC_vec[1][2] };
	                                                        
		double p3[3] = {theVertices[k].r[0]	+pbc2[0] * PBC_vec[0][0]*alpha + pbc2[1] * PBC_vec[1][0],
		theVertices[k].r[1]	+pbc2[0] * PBC_vec[0][1]*alpha + pbc2[1] * PBC_vec[1][1],
		theVertices[k].r[2]	+pbc2[0] * PBC_vec[0][2]*alpha + pbc2[1] * PBC_vec[1][2] };
		
		double p1b[3] = { theVerticesB[i].r[0], theVerticesB[i].r[1], theVerticesB[i].r[2] };

		double p2b[3] = {theVerticesB[j].r[0]	+pbc1[0] * PBC_vec[0][0]*alpha + pbc1[1] * PBC_vec[1][0],
		theVerticesB[j].r[1]	+pbc1[0] * PBC_vec[0][1]*alpha + pbc1[1] * PBC_vec[1][1],
		theVerticesB[j].r[2]	+pbc1[0] * PBC_vec[0][2]*alpha + pbc1[1] * PBC_vec[1][2] };
	                                                        
		double p3b[3] = {theVerticesB[k].r[0]	+pbc2[0] * PBC_vec[0][0]*alpha + pbc2[1] * PBC_vec[1][0],
		theVerticesB[k].r[1]	+pbc2[0] * PBC_vec[0][1]*alpha + pbc2[1] * PBC_vec[1][1],
		theVerticesB[k].r[2]	+pbc2[0] * PBC_vec[0][2]*alpha + pbc2[1] * PBC_vec[1][2] };

		double del = 1.0/9.0;

//		for( double i = del; i < 1.0-del+1e-10; i += del )
//		for( double j = del; j < 1.0-del+1e-10; j += del )

		double vertices[9] = { p1[0],p1[1],p1[2],p2[0],p2[1],p2[2],p3[0],p3[1],p3[2] };
		double fi = 1.0 / 3.0;
		double fj = 1.0 / 3.0;
		double color_min = thickness - 2;
		double color_max = thickness+2;
		double color1[3] = {0,0,1};
		double color2[3] = {1,0,0};
		{
			double k = 1-fi-fj;

			if( k < 0 ) continue;

			double p[3] = { p1[0]*fi+p2[0]*fj+p3[0]*k,
					p1[1]*fi+p2[1]*fj+p3[1]*k,
					p1[2]*fi+p2[2]*fj+p3[2]*k }; 
			double pb[3] = { p1b[0]*fi+p2b[0]*fj+p3b[0]*k,
					 p1b[1]*fi+p2b[1]*fj+p3b[1]*k,
					 p1b[2]*fi+p2b[2]*fj+p3b[2]*k }; 
				
			double dr[3] = { p[0] - pb[0], p[1] - pb[1], p[2] - pb[2] };	

			double h = p[2];
			double f = (h - color_min)/(color_max-color_min); 


			if( f < 0 ) f = 0;
			if( f > 1 ) f = 1;
	
			double color[3] = { (1-f)*color1[0] + f * color2[0],
					    (1-f)*color1[1] + f * color2[1],
					    (1-f)*color1[2] + f * color2[2] };
	
			if( f < 0.5 )
			{
				f*=2;
				color[0] = (1-f)*color1[0];
				color[1] = (1-f)*color1[1];
				color[2] = (1-f)*color1[2];
			}
			else
			{
				f-=0.5;
				f*=2;
				color[0] = f*color2[0];
				color[1] = f*color2[1];
				color[2] = f*color2[2];
			}
		
			fprintf(psFile, "%lf %lf %lf setrgbcolor\n",
				color[0], color[1], color[2] );

			for( int pass = 0; pass < 2; pass++ )
			{
				fprintf(psFile, "newpath\n");
				fprintf(psFile, "%lf %lf moveto\n", (vertices[0]+shift[0])*scale, (vertices[1]+shift[1])*scale );
				fprintf(psFile, "%lf %lf lineto\n", (vertices[3]+shift[0])*scale, (vertices[4]+shift[1])*scale );
				fprintf(psFile, "%lf %lf lineto\n", (vertices[6]+shift[0])*scale, (vertices[7]+shift[1])*scale );
				fprintf(psFile, "closepath\n");
				if( pass == 1 )
					fprintf(psFile, "fill\n");
				else
					fprintf(psFile, "stroke\n");
			}
		}
	}
	}	
}

void surface::writeXYZSurface( const char *fileNameXYZ, const char *fileNamePSF, surface *theSurface )
{
	int nv = theSurface->nv;
	
	int maxbonds = theSurface->nv*MAX_VALENCE + theSurface->nt;

	int *bond_list = (int*)malloc( sizeof(int) * maxbonds *2 );
	int nbonds =0;

	FILE *theXYZ = fopen(fileNameXYZ, "w");
	FILE *thePSF = fopen(fileNamePSF,"w");

	int p_base = 0;
	int nbx = 0;
	
	fprintf(thePSF, "PSF EXT\n");
	fprintf(thePSF, "\n");

	
	char tbuf[256];
	sprintf(tbuf, "%d", nv );
	
	for( int x = 0; x < 8 - strlen(tbuf); x++ )
		fprintf(thePSF, " ");
	fprintf(thePSF, "%s", tbuf );	
	
	fprintf(thePSF, " !NATOM\n");

	
	


	fprintf(theXYZ, "%d\n", nv );
	fprintf(theXYZ, "dual surface\n");
	for( int i = 0; i < theSurface->nv; i++ )
	{
		fprintf(theXYZ, "C %lf %lf %lf # upper i: %d\n", 
			theSurface->theVertices[i].r[0], 
			theSurface->theVertices[i].r[1], 
			theSurface->theVertices[i].r[2], i );
		fprintPSFAtomExt( thePSF, 1+i, 1+i, "C", "C" );

		int val = theSurface->theVertices[i].valence;
//		printf("%d valence: %d\n", i, val );

		for( int e = 0; e < val; e++ )
		{
			int j = theSurface->theVertices[i].edges[e];


			if( j > i && 
			fabs(theSurface->theVertices[i].edge_PBC[3*e+0]) < 1e-5 &&
			fabs(theSurface->theVertices[i].edge_PBC[3*e+2]) < 1e-5 &&
			fabs(theSurface->theVertices[i].edge_PBC[3*e+1]) < 1e-5  )
//			if( j > i )
			{
				bond_list[2*nbonds+0] = i;
				bond_list[2*nbonds+1] = j;
				nbonds++;
			}
		}
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

	free(bond_list);
	fclose(theXYZ);
	fclose(thePSF);
}

void surface::saveLimitingSurface( const char *fileName )
{
	FILE *theFile = fopen(fileName, "w");

	fprintf(theFile, "R 1\n");
	fprintf(theFile, "%lf %lf %lf\n", PBC_vec[0][0], PBC_vec[0][1], PBC_vec[0][2] );
	fprintf(theFile, "%lf %lf %lf\n", PBC_vec[1][0], PBC_vec[1][1], PBC_vec[1][2] );
	fprintf(theFile, "%lf %lf %lf\n", PBC_vec[2][0], PBC_vec[2][1], PBC_vec[2][2] );
	
	for( int v = 0; v < nv; v++ )
	{
		int val = theVertices[v].valence;

		
		double lr[3] = { 0,0,0};

		lr[0] = 0.5 * theVertices[v].r[0]; 
		lr[1] = 0.5 * theVertices[v].r[1]; 
		lr[2] = 0.5 * theVertices[v].r[2]; 

		double w = 1.0 / (val * 2);

		for( int e = 0; e < val; e++ )
		{
			int j = theVertices[v].edges[e];

			lr[0] += theVertices[j].r[0]*w;
			lr[1] += theVertices[j].r[1]*w;
			lr[2] += theVertices[j].r[2]*w;
		}
	


		fprintf(theFile, "%d %lf %lf %lf %d", v, lr[0], lr[1], lr[2], 
			theVertices[v].valence );

		for( int e = 0; e < theVertices[v].valence; e++ )
			fprintf(theFile, " %d", theVertices[v].edges[e] );
		fprintf(theFile, "\n");
	}

	int ntri = 0;
	
	for( int i = 0; i < nv; i++ )
	{
		int val = theVertices[i].valence;
		
		for( int e = 0; e < val; e++ )
		{
			int j = theVertices[i].edges[e];

			int ep1 = e+1;
			if( ep1 >= val ) ep1 -= val;

			int k = theVertices[i].edges[ep1];

			if( i < j && i < k )
				ntri++;
		}
	}

	fprintf(theFile, "ntri %d\n", ntri );
	
	for( int i = 0; i < nv; i++ )
	{
		int val = theVertices[i].valence;
		
		for( int e = 0; e < val; e++ )
		{
			int j = theVertices[i].edges[e];

			int ep1 = e+1;
			if( ep1 >= val ) ep1 -= val;

			int k = theVertices[i].edges[ep1];

			if( i < j && i < k )
			{
				int pj = j;
				int pk = k;
				if( j > k )
				{
					pj = k;
					pk = j;
				}
				fprintf(theFile, "%d %d %d\n", i, pj, pk );
			}
		}
	}

	
}


void surface::saveSurface( const char *fileName )
{
	FILE *theFile = fopen(fileName, "w");

	fprintf(theFile, "3D saved surface\n");
	fprintf(theFile, "%lf %lf %lf\n", PBC_vec[0][0], PBC_vec[0][1], PBC_vec[0][2] );
	fprintf(theFile, "%lf %lf %lf\n", PBC_vec[1][0], PBC_vec[1][1], PBC_vec[1][2] );
	fprintf(theFile, "%lf %lf %lf\n", PBC_vec[2][0], PBC_vec[2][1], PBC_vec[2][2] );
	
	for( int v = 0; v < nv; v++ )
	{
		fprintf(theFile, "%d %lf %lf %lf %d", v, theVertices[v].r[0], theVertices[v].r[1], theVertices[v].r[2], 
			theVertices[v].valence );

		for( int e = 0; e < theVertices[v].valence; e++ )
			fprintf(theFile, " %d", theVertices[v].edges[e] );
		fprintf(theFile, "\n");
	}

	int ntri = 0;
	
	for( int i = 0; i < nv; i++ )
	{
		int val = theVertices[i].valence;
		
		for( int e = 0; e < val; e++ )
		{
			int j = theVertices[i].edges[e];

			int ep1 = e+1;
			if( ep1 >= val ) ep1 -= val;

			int k = theVertices[i].edges[ep1];

			if( i < j && i < k )
				ntri++;
		}
	}

	fprintf(theFile, "ntri %d\n", ntri );
	
	for( int i = 0; i < nv; i++ )
	{
		int val = theVertices[i].valence;
		
		for( int e = 0; e < val; e++ )
		{
			int j = theVertices[i].edges[e];

			int ep1 = e+1;
			if( ep1 >= val ) ep1 -= val;

			int k = theVertices[i].edges[ep1];

			if( i < j && i < k )
			{
				int pj = j;
				int pk = k;
				if( j > k )
				{
					pj = k;
					pk = j;
				}
				fprintf(theFile, "%d %d %d\n", i, pj, pk );
			}
		}
	}

	
}

void surface::writeLimitingSurfacePSF( FILE *theFile, pcomplex **allComplexes, int ncomplexes)
{
	int npts = 0;

	int ppd = 0;


	for( int fi = 0; fi <= plim; fi++ )
	for( int fj = 0; fj <= plim-fi; fj++ )
	{
		double f1 = fi / (double)(plim+1e-14);
		double f2 = fj / (double)(plim+1e-14);
		
		ppd++;
	}

	npts = ppd * nt;

	int nbonds = 0;	
		
	int index[plim+1][plim+1];
	for( int fi = 0; fi <= plim; fi++)
	for( int fj = 0; fj <= plim; fj++)
		index[fi][fj] = -1;

	int t_index = 0;
	for( int fi = 0; fi <= plim; fi++ )
	for( int fj = 0; fj <= plim-fi; fj++)
	{
		index[fi][fj] = t_index;
		t_index++;
	}

	int *the_bonds = NULL;


	for( int pass = 0; pass < 2; pass++ )
	{
		if( pass == 1 )
			the_bonds = (int *)malloc( sizeof(int) * 2 * nbonds );

		nbonds=0;

		int i_tot = 0;
		for( int t = 0; t < nt; t++ )
		{	
			int base_index = i_tot;

			for( int fi = 0; fi <= plim; fi++ )
			for( int fj = 0; fj <= plim-fi; fj++ )
			{
				double f1 = fi / (double)plim;
				double f2 = fj / (double)plim;
			
				int bonds[3][2] = {
					{fi+1, fj-1},
					{fi+1, fj},
					{ fi, fj+1}
				};
	
				for( int b = 0; b < 3; b++ )
				{
					if( bonds[b][0] >= 0 && bonds[b][1] >= 0 && bonds[b][0]+bonds[b][1] <= plim )
					{
						int i1 = base_index + index[fi][fj];
						int i2 = base_index + index[bonds[b][0]][bonds[b][1]];
						
						if( pass == 1 )
						{
							the_bonds[2*nbonds+0] = i1; 
							the_bonds[2*nbonds+1] = i2; 
						}
					
						nbonds++;
					}
				}
		
			}

			i_tot += ppd; 
		}

		int curp = npts;

		for( int c = 0; c < ncomplexes; c++ )
		{
			int c_nbonds = allComplexes[c]->getNBonds();
			
			if( pass == 1 )
			{
				allComplexes[c]->putBonds( the_bonds + nbonds*2 );
	
				for( int x = 0; x < c_nbonds; x++ )
				{
					the_bonds[2*nbonds+2*x] += curp;
					the_bonds[2*nbonds+2*x+1] += curp;
	
				}
				nbonds += c_nbonds;
		
			}
			else
				nbonds += c_nbonds;
			
			curp += allComplexes[c]->nsites;
		}
	}

	int cpts = 0;
	for( int c = 0; c < ncomplexes; c++ )
		cpts += allComplexes[c]->nsites;


	writePSF( theFile, npts+cpts, NULL, the_bonds, nbonds );

	free(the_bonds);
}

void surface::writeLimitingSurface( FILE *theFile, pcomplex **allComplexes, int ncomplexes, double *alphas )
{
	double *r = (double *)malloc( sizeof(double) * 3 * (nv+1) );

	get(r);

	if( !alphas )
	{
		r[3*nv+0] = 1.0;
		r[3*nv+1] = 1.0;
		r[3*nv+2] = 1.0;
	}
	else
	{
		r[3*nv+0] = alphas[0];
		r[3*nv+1] = alphas[1];
		r[3*nv+2] = alphas[2];
	}
	int npts = 0;

	int ppd = 0;


	for( int fi = 0; fi <= plim; fi++ )
	for( int fj = 0; fj <= plim-fi; fj++ )
	{
		double f1 = fi / (double)plim;
		double f2 = fj / (double)plim;
		
		ppd++;
	}

	npts = ppd * nt;
	
	for( int c = 0; c < ncomplexes; c++ )
		npts += allComplexes[c]->nsites;

	fprintf(theFile, "%d\n", npts );
	fprintf(theFile, "limiting surface\n");
	
	for( int t = 0; t < nt; t++ )
	{	
		for( int fi = 0; fi <= plim; fi++ )
		for( int fj = 0; fj <= plim-fi; fj++ )
		{
			double f1 = fi / (double)plim;
			double f2 = fj / (double)plim;
			
			double rl[3],nrm[3];
		
			evaluateRNRM( t, f1, f2, rl, nrm, r );
	
	
			fprintf(theFile, "C %lf %lf %lf\n", rl[0], rl[1], rl[2] );
		}
	}

	free(r);

	for( int c = 0; c < ncomplexes; c++ )
	{
		for( int p = 0; p < allComplexes[c]->nsites; p++ )
			fprintf(theFile, "O %lf %lf %lf\n", allComplexes[c]->rall[3*p+0], allComplexes[c]->rall[3*p+1], allComplexes[c]->rall[3*p+2] );
	}
	fflush(theFile);
}


void surface::writeLimitTriangles( FILE *theFile )
{
	double *rl = (double *)malloc( sizeof(double) * 3 * nv );

	for( int v = 0; v < nv; v++ )
	{
		int val = theVertices[v].valence;

		double lr[3] = { 0,0,0};

		double *ti = theVertices[v].r;

		lr[0] = 0.5 * theVertices[v].r[0]; 
		lr[1] = 0.5 * theVertices[v].r[1]; 
		lr[2] = 0.5 * theVertices[v].r[2]; 

		double w = 1.0 / (val * 2);

		for( int e = 0; e < val; e++ )
		{
			int j = theVertices[v].edges[e];

			double dr[3] = { theVertices[j].r[0] - ti[0], theVertices[j].r[1] - ti[1], theVertices[j].r[2] - ti[2] };

			double put[3];	
			MinImage3D( dr, PBC_vec, put );
			dr[0] += ti[0];
			dr[1] += ti[1];
			dr[2] += ti[2];

			lr[0] += dr[0]*w;
			lr[1] += dr[1]*w;
			lr[2] += dr[2]*w;
		}

		rl[3*v+0] = lr[0];
		rl[3*v+1] = lr[1];
		rl[3*v+2] = lr[2];
	}

	double min_val = 1e10;
	double max_val = -1e10;

	double avn = 0;

	for( int t = 0; t < nt; t++ )
	{
			
	}

	for( int t = 0; t < nt; t++ )
	{
		int i = theTriangles[t].ids[0];
		int j = theTriangles[t].ids[1];
		int k = theTriangles[t].ids[2];

		double v1 = (theVertices[i].store_val+theVertices[j].store_val+theVertices[k].store_val);
		double v2 = (theVertices[i].store_val2+theVertices[j].store_val2+theVertices[k].store_val2);
	
		double vi = theVertices[i].store_val / (theVertices[i].store_val2+1e-10);
		double vj = theVertices[j].store_val / (theVertices[j].store_val2+1e-10);
		double vk = theVertices[k].store_val / (theVertices[k].store_val2+1e-10);
	
		double val = (vi+vj+vk)/3;
		if( val < min_val )
			min_val = val;
		if( val > max_val )
			max_val = val;
	}

	int color_by_value = 0;

	if( fabs(max_val-min_val) > 1e-4 )
	{
		min_val = -0.03;
		max_val = 0.03;
		color_by_value = 1;
	}
	int range = 300;
	int offset = 17;
	if( !color_by_value )
		fprintf( theFile, "draw color blue\n" );
	else 
	{
		// load index colors.
		double color1[3] = { 0, 0, 1 };
		double color2[3] = { 1, 0, 0 };
		for( int i = 0; i < range; i++ )
		{	
			if( i < range/2 )
			{
				double f = (range/2.0-i)/(range/2.0);
				fprintf( theFile, "color change rgb %d %lf %lf %lf\n", i+offset,
					f * color1[0] + (1-f), f* color1[1] + (1-f), f*color1[2] + (1-f) );
			}
			else 
			{
				double f = (i-range/2.0)/(range/2.0);
				fprintf( theFile, "color change rgb %d %lf %lf %lf\n", i+offset,
					f * color2[0] + 1-f, f* color2[1] + 1-f, f*color2[2] + 1-f );
			}
		}
	}


	for( int t = 0; t < nt; t++ )
	{
		int i = theTriangles[t].ids[0];
		int j = theTriangles[t].ids[1];
		int k = theTriangles[t].ids[2];

		if( !theTriangles[t].draw )
		{
			continue;
		}
	
		double dr1[3] = { rl[3*j+0] - rl[3*i+0], rl[3*j+1] - rl[3*i+1], rl[3*j+2] - rl[3*i+2] };
		double dr2[3] = { rl[3*k+0] - rl[3*i+0], rl[3*k+1] - rl[3*i+1], rl[3*k+2] - rl[3*i+2] };

		double put1[3];
		double put2[3];

		MinImage3D( dr1, PBC_vec, put1 );
		MinImage3D( dr2, PBC_vec, put2 );

		double v1 = fabs(put1[0]) + fabs(put1[1]) + fabs(put1[2]);
		double v2 = fabs(put2[0]) + fabs(put2[1]) + fabs(put2[2]);

		double color_below = 0.2;
		double color_above = 0.6;
		if( v1 < 1e-3 && v2 < 1e-3 )
		{
			if( color_by_value )
			{
				double vi = theVertices[i].store_val / (theVertices[i].store_val2+1e-10);
				double vj = theVertices[j].store_val / (theVertices[j].store_val2+1e-10);
				double vk = theVertices[k].store_val / (theVertices[k].store_val2+1e-10);
	
				double val = (vi+vj+vk)/3;

				double f = (val-min_val) / (max_val-min_val);
				//double scale_f = (f-color_below)/(color_above-color_below);
				//f = scale_f;
				if( f > 0.999 ) f = 0.999;
				if( f < 0 ) f = 0;
	
				int color = offset + f * range;
	
				fprintf( theFile, "draw color %d\n", color );				
			}

			fprintf( theFile, "draw triangle { %lf %lf %lf } { %lf %lf %lf } { %lf %lf %lf }\n",
				rl[3*i+0], rl[3*i+1], rl[3*i+2],
				rl[3*j+0], rl[3*j+1], rl[3*j+2],
				rl[3*k+0], rl[3*k+1], rl[3*k+2] );
		}
	}	

	free(rl);
}

void surface::writeXYZandPSFPeriodic( const char *baseName )
{
	char fileName[256];

	sprintf(fileName, "%s.psf", baseName );
	FILE *thePSF = fopen(fileName, "w");
	sprintf(fileName, "%s.xyz", baseName );
	FILE *theXYZ = fopen(fileName, "w");

	double totv = 0;
	for( int x = 0; x < nv; x++ )
		totv += theVertices[x].valence;
	totv /= 2;

	int *all_bonds = (int *)malloc( sizeof(int) * totv *3 );

	int b = 0;

	int nvSpace = nv;
	double *coord_space = (double *)malloc( sizeof(double) * nvSpace * 3 );
	int nv_write = nv;

	for( int i = 0; i < nv; i++ )
	{
		coord_space[3*i+0] = theVertices[i].r[0];
		coord_space[3*i+1] = theVertices[i].r[1];
		coord_space[3*i+2] = theVertices[i].r[2];

		for( int e = 0; e < theVertices[i].valence; e++ )
		{
			int j = theVertices[i].edges[e];
			if( j < i ) continue;

			double dr[3] = { 
				theVertices[j].r[0] - theVertices[i].r[0],
				theVertices[j].r[1] - theVertices[i].r[1],
				theVertices[j].r[2] - theVertices[i].r[2] };
			double put[3];
			MinImage3D( dr, PBC_vec, put );

			if( put[0]*put[0]+put[1]*put[1]+put[2]*put[2] < 1e-3 )
			{
				all_bonds[2*b+0] = i;
				all_bonds[2*b+1] = j;
				b++;
			}
			else
			{
				if( nvSpace == nv_write )
				{
					nvSpace *= 2;
					coord_space = (double *)realloc( coord_space, sizeof(double) * 3 * nvSpace );	
				}

				coord_space[nv_write*3+0] = theVertices[i].r[0] + dr[0];
				coord_space[nv_write*3+1] = theVertices[i].r[1] + dr[1];
				coord_space[nv_write*3+2] = theVertices[i].r[2] + dr[2];
				
				all_bonds[2*b+0] = i;
				all_bonds[2*b+1] = nv_write;

				nv_write++;
				b++;
			}
		} 
	} 

	writePSF( thePSF, nv_write, NULL, all_bonds, b );

	fprintf(theXYZ, "%d\n", nv_write );
	fprintf(theXYZ, "PSF and XYZ, periodic\n");
	for( int i = 0; i < nv_write; i++ )
		fprintf(theXYZ, "C %lf %lf %lf\n", coord_space[3*i+0], coord_space[3*i+1], coord_space[3*i+2] );
	fclose(theXYZ);
	fclose(thePSF);

	free(all_bonds);
	free(coord_space);
}

void surface::writeLimitPSF( FILE *theFile, int np, double dist_nrm)
{
	double totv = 0;
	for( int x = 0; x < nv; x++ )
		totv += theVertices[x].valence;
	totv /= 2;

	int nbonds = totv;
	if( dist_nrm > 0 )
		nbonds += np;
		

	int *all_bonds = (int *)malloc( sizeof(int) * nbonds *2 );
	
	int b = 0;
	for( int i = 0; i < nv; i++ )
	{
		for( int e = 0; e < theVertices[i].valence; e++ )
		{
			int j = theVertices[i].edges[e];
			if( j < i ) continue;

			double dr[3] = { 
				theVertices[i].r[0] - theVertices[j].r[0],
				theVertices[i].r[1] - theVertices[j].r[1],
				theVertices[i].r[2] - theVertices[j].r[2] };
			double put[3];
			MinImage3D( dr, PBC_vec, put );

			if( put[0]*put[0]+put[1]*put[1]+put[2]*put[2] < 1e-3 )
			{
				all_bonds[2*b+0] = i;
				all_bonds[2*b+1] = j;
				b++;
			}
		} 
	} 

	if( dist_nrm > 0 )
	{
		for( int px = 0; px < np; px++ )	
		{
			all_bonds[2*b+0] = nv + 2*px;
			all_bonds[2*b+1] = nv + 2*px+1;
			b++;
		}
	}
	if( np > 0 )
	{
		int np_print = np;
		if( dist_nrm > 0 )
			np_print += np;
		char *atomNames = (char *)malloc( sizeof(char) * (1+nv+np_print ) );
		for( int p = 0; p < nv; p++ )
			atomNames[p] = 'C';
		if( dist_nrm > 0 )
		{
			for( int p = nv; p < nv+np*2; p += 2 )
			{
				atomNames[p] = 'O';	
				atomNames[p+1] = 'N';	
			}
		}
		else
		{	
			for( int p = nv; p < nv+np; p++ )
				atomNames[p] = 'O';	
		}
		writePSF( theFile, nv+np_print, atomNames, all_bonds, b );
		free(atomNames);
	}
	else
		writePSF( theFile, nv, NULL, all_bonds, b );

	free(all_bonds);
}

void surface::writeBARLimitStructure( FILE *theFile, int *pf, double *puv, int np )
{
	fprintf(theFile, "%d\n", nv+np );
	fprintf(theFile, "BAR code limit structure\n");

	for( int v = 0; v < nv; v++ )
	{
		int val = theVertices[v].valence;

		double lr[3] = { 0,0,0};

		double *ti = theVertices[v].r;

		lr[0] = 0.5 * theVertices[v].r[0]; 
		lr[1] = 0.5 * theVertices[v].r[1]; 
		lr[2] = 0.5 * theVertices[v].r[2]; 

		double w = 1.0 / (val * 2);

		for( int e = 0; e < val; e++ )
		{
			int j = theVertices[v].edges[e];

			double dr[3] = { theVertices[j].r[0] - ti[0], theVertices[j].r[1] - ti[1], theVertices[j].r[2] - ti[2] };

			double put[3];	
			MinImage3D( dr, PBC_vec, put );
			dr[0] += ti[0];
			dr[1] += ti[1];
			dr[2] += ti[2];

			lr[0] += dr[0]*w;
			lr[1] += dr[1]*w;
			lr[2] += dr[2]*w;
		}

		fprintf(theFile, "C %lf %lf %lf\n", lr[0], lr[1], lr[2] );	
	}
	
	if( np > 0 ) {
	
		double * r = (double *)malloc( sizeof(double) * 3 * (nv +1 ) );
		get(r);
		r[3*nv+0] = 1.0;
		r[3*nv+1] = 1.0;
		r[3*nv+2] = 1.0;

		double prev[3] = { 0,0,0};
		for( int p = 0; p < np; p++ )
		{
			int frm = pf[p]*nf_g_q_p;
	
			int *cp = theFormulas[frm].cp;
	
			double u = puv[2*p+0];
			double v = puv[2*p+1];
	
			
			double rp[3], nrm[3];
			evaluateRNRM( pf[p], u, v, rp, nrm,  r ); 

			if( p % 3 != 0)
			{
				double dr[3] = { rp[0] - prev[0], rp[1] - prev[1], rp[2] - prev[2] };
	
				double put[3];
				MinImage3D( dr, PBC_vec, put );

				rp[0] = prev[0] + dr[0];
				rp[1] = prev[1] + dr[1];
				rp[2] = prev[2] + dr[2];
			}

			if( p % 3 == 1 )	
				fprintf(theFile, "O %lf %lf %lf\n", rp[0], rp[1], rp[2] );	
			else
				fprintf(theFile, "N %lf %lf %lf\n", rp[0], rp[1], rp[2] );	

			prev[0] = rp[0];
			prev[1] = rp[1];
			prev[2] = rp[2];
		}

		free(r);
	}


}

void surface::writeLimitStructure( FILE *theFile, int *pf, double *puv, int np, double dist_nrm)
{
	int extra_np = (dist_nrm > 0 ? np : 0 );
	fprintf(theFile, "%d\n", nv+np+extra_np );
	fprintf(theFile, "limit structure\n");

	for( int v = 0; v < nv; v++ )
	{
		int val = theVertices[v].valence;

		double lr[3] = { 0,0,0};

		double *ti = theVertices[v].r;

		lr[0] = 0.5 * theVertices[v].r[0]; 
		lr[1] = 0.5 * theVertices[v].r[1]; 
		lr[2] = 0.5 * theVertices[v].r[2]; 

		double w = 1.0 / (val * 2);

		for( int e = 0; e < val; e++ )
		{
			int j = theVertices[v].edges[e];

			double dr[3] = { theVertices[j].r[0] - ti[0], theVertices[j].r[1] - ti[1], theVertices[j].r[2] - ti[2] };

			double put[3];	
			MinImage3D( dr, PBC_vec, put );
			dr[0] += ti[0];
			dr[1] += ti[1];
			dr[2] += ti[2];

			lr[0] += dr[0]*w;
			lr[1] += dr[1]*w;
			lr[2] += dr[2]*w;
		}

		fprintf(theFile, "C %lf %lf %lf\n", lr[0], lr[1], lr[2] );	
	}
	
	if( np > 0 ) {
	
		double * r = (double *)malloc( sizeof(double) * 3 * (nv +1 ) );
		get(r);
		r[3*nv+0] = 1.0;
		r[3*nv+1] = 1.0;
		r[3*nv+2] = 1.0;
		for( int p = 0; p < np; p++ )
		{
			int frm = pf[p]*nf_g_q_p;
	
			int *cp = theFormulas[frm].cp;
	
			double u = puv[2*p+0];
			double v = puv[2*p+1];
	
			
			double rp[3], nrm[3];
			evaluateRNRM( pf[p], u, v, rp, nrm,  r ); 
	
			fprintf(theFile, "O %lf %lf %lf\n", rp[0], rp[1], rp[2] );	
			if( dist_nrm > 0 )
			{
				fprintf(theFile, "N %lf %lf %lf\n", 
					rp[0] + dist_nrm * nrm[0],
					rp[1] + dist_nrm * nrm[1],
					rp[2] + dist_nrm * nrm[2] );
			}
		}

		free(r);
	}


}
	
void surface::writeStructure( FILE *theFile )
{
	double alpha_x = 1.0;
	double alpha_y = 1.0;
	double alpha_z = 1.0;
	int nf = nf_faces*nf_g_q_p;
	fprintf(theFile, "%d\n", nf + nf_irr_faces*nf_irr_pts );
	fprintf(theFile, "structure\n");
	
	for( int f = 0; f < nf_irr_faces*nf_irr_pts; f++ )
	{
		double r[3] = {0,0,0};
		double ru[3] = {0,0,0};
		double rv[3] = {0,0,0};
		
		int np = theIrregularFormulas[f].ncoor;

		int *cp = theIrregularFormulas[f].cp;
		for( int p = 0; p < np; p++ )
		{
			r[0] += theIrregularFormulas[f].r_w[p] * alpha_x*(theVertices[cp[p]].r[0] + theIrregularFormulas[f].r_pbc[3*p+0]); 
			r[1] += theIrregularFormulas[f].r_w[p] * alpha_y*(theVertices[cp[p]].r[1] + theIrregularFormulas[f].r_pbc[3*p+1]); 
			r[2] += theIrregularFormulas[f].r_w[p] * alpha_z*(theVertices[cp[p]].r[2] + theIrregularFormulas[f].r_pbc[3*p+2]); 
			
			ru[0] += theIrregularFormulas[f].r_u[p] * alpha_x*(theVertices[cp[p]].r[0] + theIrregularFormulas[f].r_pbc[3*p+0]); 
			ru[1] += theIrregularFormulas[f].r_u[p] * alpha_y*(theVertices[cp[p]].r[1] + theIrregularFormulas[f].r_pbc[3*p+1]); 
			ru[2] += theIrregularFormulas[f].r_u[p] * alpha_z*(theVertices[cp[p]].r[2] + theIrregularFormulas[f].r_pbc[3*p+2]); 
			
			rv[0] += theIrregularFormulas[f].r_v[p] * alpha_x*(theVertices[cp[p]].r[0] + theIrregularFormulas[f].r_pbc[3*p+0]); 
			rv[1] += theIrregularFormulas[f].r_v[p] * alpha_y*(theVertices[cp[p]].r[1] + theIrregularFormulas[f].r_pbc[3*p+1]); 
			rv[2] += theIrregularFormulas[f].r_v[p] * alpha_z*(theVertices[cp[p]].r[2] + theIrregularFormulas[f].r_pbc[3*p+2]); 
		}

		normalize( ru );
		normalize( rv );

		double del = 0.1;
		double nrm[3];
		cross( ru, rv, nrm );
		normalize(nrm);
		fprintf(theFile, "N %lf %lf %lf\n", r[0], r[1], r[2] );
//		fprintf(theFile, "O %lf %lf %lf\n", r[0]+del * nrm[0], r[1]+del * nrm[1], r[2]+del * nrm[2] );
	}

	for( int f = 0; f < nf; f++ )
	{
		double r[3] = {0,0,0};
		double ru[3] = {0,0,0};
		double rv[3] = {0,0,0};
		
		int np = theFormulas[f].ncoor;

		int *cp = theFormulas[f].cp;
		for( int p = 0; p < np; p++ )
		{
			r[0] += theFormulas[f].r_w[p] * alpha_x*(theVertices[cp[p]].r[0] + theFormulas[f].r_pbc[3*p+0]); 
			r[1] += theFormulas[f].r_w[p] * alpha_y*(theVertices[cp[p]].r[1] + theFormulas[f].r_pbc[3*p+1]); 
			r[2] += theFormulas[f].r_w[p] * alpha_z*(theVertices[cp[p]].r[2] + theFormulas[f].r_pbc[3*p+2]); 
			
			ru[0] += theFormulas[f].r_u[p] * alpha_x*(theVertices[cp[p]].r[0] + theFormulas[f].r_pbc[3*p+0]); 
			ru[1] += theFormulas[f].r_u[p] * alpha_y*(theVertices[cp[p]].r[1] + theFormulas[f].r_pbc[3*p+1]); 
			ru[2] += theFormulas[f].r_u[p] * alpha_z*(theVertices[cp[p]].r[2] + theFormulas[f].r_pbc[3*p+2]); 
			
			rv[0] += theFormulas[f].r_v[p] * alpha_x*(theVertices[cp[p]].r[0] + theFormulas[f].r_pbc[3*p+0]); 
			rv[1] += theFormulas[f].r_v[p] * alpha_y*(theVertices[cp[p]].r[1] + theFormulas[f].r_pbc[3*p+1]); 
			rv[2] += theFormulas[f].r_v[p] * alpha_z*(theVertices[cp[p]].r[2] + theFormulas[f].r_pbc[3*p+2]); 
		}

		normalize( ru );
		normalize( rv );

		double del = 0.1;
		double nrm[3];
		cross( ru, rv, nrm );
		normalize(nrm);
		fprintf(theFile, "C %lf %lf %lf\n", r[0], r[1], r[2] );
//		fprintf(theFile, "O %lf %lf %lf\n", r[0]+del * nrm[0], r[1]+del * nrm[1], r[2]+del * nrm[2] );
	}
}


void surface::writeVertexXYZandPSFPeriodic( const char *baseName )
{
	char fileName[256];

	sprintf(fileName, "%s.psf", baseName );
	FILE *thePSF = fopen(fileName, "w");
	sprintf(fileName, "%s.xyz", baseName );
	FILE *theXYZ = fopen(fileName, "w");

	double totv = 0;
	for( int x = 0; x < nv; x++ )
		totv += theVertices[x].valence;
	totv /= 2;

	int *all_bonds = (int *)malloc( sizeof(int) * totv *3 );

	int b = 0;

	int nvSpace = nv;
	double *coord_space = (double *)malloc( sizeof(double) * nvSpace * 3 );
	int nv_write = nv;

	for( int i = 0; i < nv; i++ )
	{
		coord_space[3*i+0] = theVertices[i].r[0];
		coord_space[3*i+1] = theVertices[i].r[1];
		coord_space[3*i+2] = theVertices[i].r[2];
		
		for( int e = 0; e < theVertices[i].valence; e++ )
		{
			int j = theVertices[i].edges[e];
	//		if( j < i ) continue;
			if( i < j ) continue;
			all_bonds[2*b+0] = i;
			all_bonds[2*b+1] = j;
			b++;
		}
	} 

	writePSF( thePSF, nv_write, NULL, all_bonds, b );

	fprintf(theXYZ, "%d\n", nv_write );
	fprintf(theXYZ, "PSF and XYZ, periodic\n");
	for( int i = 0; i < nv_write; i++ )
		fprintf(theXYZ, "C %lf %lf %lf\n", coord_space[3*i+0], coord_space[3*i+1], coord_space[3*i+2] );
	fclose(theXYZ);
	fclose(thePSF);

	free(all_bonds);
	free(coord_space);
}

void Simulation::writeLimitingSurface( FILE *theFile )
{
	int npts = 0;
	int ppd = 0;

	for( int fi = 0; fi <= plim; fi++ )
	for( int fj = 0; fj <= plim-fi; fj++ )
	{
		double f1 = fi / (double)plim;
		double f2 = fj / (double)plim;
		
		ppd++;
	}

	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
		npts += ppd * sRec->theSurface->nt;
	
	npts += nsites_at_psfwrite + visualization_cache;

	fprintf(theFile, "%d\n", npts );
	fprintf(theFile, "limiting surface\n");

	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		surface *theSurface = sRec->theSurface;
		int nv = theSurface->nv;
		int nt = theSurface->nt;
		double *r = (double *)malloc( sizeof(double) * 3 * (nv+1) );
	
		r[3*nv+0] = alpha[0];
		r[3*nv+1] = alpha[1];
		r[3*nv+2] = alpha[2];
		theSurface->get(r);
		
		for( int t = 0; t < nt; t++ )
		{	
			for( int fi = 0; fi <= plim; fi++ )
			for( int fj = 0; fj <= plim-fi; fj++ )
			{
				double f1 = fi / (double)plim;
				double f2 = fj / (double)plim;
				
				double rl[3],nrm[3];
			
				theSurface->evaluateRNRM( t, f1, f2, rl, nrm, r );
		
		
				fprintf(theFile, "C %lf %lf %lf\n", rl[0], rl[1], rl[2] );
			}
		}

		free(r);
	}

#define RD_HACK_REMOVE_ME
#ifdef RD_HACK_REMOVE_ME
	// first pass, oxygen
	
	int extra = 0;
	for( int pass = 0; pass < 2; pass++ )
	{
		int nsites_written = extra;

		for( int c = 0; c < ncomplex; c++ )
		{
			if( allComplexes[c]->disabled ) continue;
	
			if( pass == 0 && !strcasecmp( allComplexes[c]->complex_name, "simpleLipid" ) )
			{
				for( int p = 0; p < allComplexes[c]->nsites; p++ )
				{
					fprintf(theFile, "O %lf %lf %lf\n", allComplexes[c]->rall[3*p+0], allComplexes[c]->rall[3*p+1], allComplexes[c]->rall[3*p+2] );
					nsites_written++;
				}
			} 
			else if( pass == 1 && strcasecmp( allComplexes[c]->complex_name, "simpleLipid" ) )
			{
				for( int p = 0; p < allComplexes[c]->nsites; p++ )
				{
					fprintf(theFile, "N %lf %lf %lf\n", allComplexes[c]->rall[3*p+0], allComplexes[c]->rall[3*p+1], allComplexes[c]->rall[3*p+2] );
					nsites_written++;
				}
			} 
		}

		if( pass == 0 )
		{
			if( nsites_written > nsites_at_psfwrite )
				extra = nsites_written - nsites_at_psfwrite;

			for( int tx = nsites_written; tx < nsites_at_psfwrite; tx++ )
				fprintf(theFile, "O %lf %lf %lf\n", 
					PBC_vec[0][0] * alpha[0],
					PBC_vec[1][1] * alpha[1],
					PBC_vec[2][2] * alpha[2] );
		}
		else
		{
			for( int tx = nsites_written; tx < visualization_cache; tx++ )
				fprintf(theFile, "N %lf %lf %lf\n", 
					PBC_vec[0][0] * alpha[0],
					PBC_vec[1][1] * alpha[1],
					PBC_vec[2][2] * alpha[2] );
		}
	}
#else	
	int nsites_written = 0;
	for( int c = 0; c < ncomplex; c++ )
	{
		if( allComplexes[c]->disabled ) continue;

		for( int p = 0; p < allComplexes[c]->nsites; p++ )
		{
			fprintf(theFile, "O %lf %lf %lf\n", allComplexes[c]->rall[3*p+0], allComplexes[c]->rall[3*p+1], allComplexes[c]->rall[3*p+2] );
			nsites_written++;
		}
	}

	for( int tx = nsites_written; tx < nsites_at_psfwrite + visualization_cache; tx++ )
		fprintf(theFile, "O %lf %lf %lf\n", 
			PBC_vec[0][0] * alpha[0],
			PBC_vec[1][1] * alpha[1],
			PBC_vec[2][2] * alpha[2] );
#endif
	fflush(theFile);

}
        		
void Simulation::writeLimitingSurfacePSF(FILE *theFile )
{
	int ppd = 0;

	for( int fi = 0; fi <= plim; fi++ )
	for( int fj = 0; fj <= plim-fi; fj++ )
	{
		double f1 = fi / (double)(plim+1e-14);
		double f2 = fj / (double)(plim+1e-14);
		
		ppd++;
	}
	
	int npts = 0;
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
		npts += ppd * sRec->theSurface->nt;
	
	int nbonds = 0;	
		
	int index[plim+1][plim+1];
	for( int fi = 0; fi <= plim; fi++)
	for( int fj = 0; fj <= plim; fj++)
		index[fi][fj] = -1;
	
	int t_index = 0;
	for( int fi = 0; fi <= plim; fi++ )
	for( int fj = 0; fj <= plim-fi; fj++)
	{
		index[fi][fj] = t_index;
		t_index++;
	}

	int *the_bonds = NULL;
	
	int disable_complex_bonds = visualization_cache > 0;

	for( int pass = 0; pass < 2; pass++ )
	{
		if( pass == 1 )
			the_bonds = (int *)malloc( sizeof(int) * 2 * nbonds );
		nbonds=0;
		int i_tot = 0;
		for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
		{
	
			for( int t = 0; t < sRec->theSurface->nt; t++ )
			{	
				int base_index = i_tot;
	
				for( int fi = 0; fi <= plim; fi++ )
				for( int fj = 0; fj <= plim-fi; fj++ )
				{
					double f1 = fi / (double)plim;
					double f2 = fj / (double)plim;
				
					int bonds[3][2] = {
						{fi+1, fj-1},
						{fi+1, fj},
						{ fi, fj+1}
					};
		
					for( int b = 0; b < 3; b++ )
					{
						if( bonds[b][0] >= 0 && bonds[b][1] >= 0 && bonds[b][0]+bonds[b][1] <= plim )
						{
							int i1 = base_index + index[fi][fj];
							int i2 = base_index + index[bonds[b][0]][bonds[b][1]];
							
							if( pass == 1 )
							{
								the_bonds[2*nbonds+0] = i1; 
								the_bonds[2*nbonds+1] = i2; 
							}
						
							nbonds++;
						}
					}
			
				}
	
				i_tot += ppd; 
			}
		}	
		int curp = npts;

		if( !disable_complex_bonds )
		{
			for( int c = 0; c < ncomplex; c++ )
			{
				int c_nbonds = allComplexes[c]->getNBonds();
				
				if( pass == 1 )
				{
					allComplexes[c]->putBonds( the_bonds + nbonds*2 );
	
					for( int x = 0; x < c_nbonds; x++ )
					{
						the_bonds[2*nbonds+2*x] += curp;
						the_bonds[2*nbonds+2*x+1] += curp;
	
					}
					nbonds += c_nbonds;
			
				}
				else
					nbonds += c_nbonds;
				
				curp += allComplexes[c]->nsites;
			}
		}
	}

	int cpts = 0;
	for( int c = 0; c < ncomplex; c++ )
		cpts += allComplexes[c]->nsites;
	
	nsites_at_psfwrite = cpts;
	cpts += visualization_cache;

	writePSF( theFile, npts+cpts, NULL, the_bonds, nbonds );

	free(the_bonds);

}


