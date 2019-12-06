#define __lipidcompositionc__
#include "interp.h"
#include "util.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "library.h"
#include "parallel.h"

static int nlipid_space = 0;
static int nlipids = 0;

int nLipidsInLibrary(void)
{
	return nlipids;
}

const fixed_lipid library_fixed[] =
{
	// Chen/Rand 1997
	// PO values copied from DO.
	{ "ZERO",         0, 68.3,  "C21" },
	{ "POPC", -1.0/87.3, 68.3,  "C21" },
	{ "DOPE", -1.0/29.4, 63.4,  "C21" },
	{ "POPE", -1.0/29.4, 58.8,  "C21" },
	{ "DOPC", -1.0/87.3, 69.7,  "C21" }, 
	{ "CHOL", -1.0/27.2, 31.0,  "C21" }, 
	{ "CHL1", -1.0/27.2, 31.0,  "O3" }, 
	{ "DOG",  -1.0/11.5, 59.0,  "C21" },
	{ "DOPS",  1.0/144.0, 71.6, "C21" }, // Fuller 2003, area from charmm f.f.
	{ "CER180",	  0, 59.0,  "C1F" },
	{ "SSM",  	0,   59.0,  "C1F" },
	{ "PSM",  	0,   59.0,  "C1F" }
};

void dumpLibrary( void )
{
	printf("Library:\n");
	printf("--------\n");
	for( int lx = 0; lx < nlipids; lx++ )
		printf(" %s c0: %lf APL: %lf\n", 
			lipidLibrary[lx].name, lipidLibrary[lx].c0, lipidLibrary[lx].APL );	
	printf("--------\n");
}


void surface::readLipidComposition( FILE *inputFile )
{
	nlipid_space = sizeof(library_fixed)/sizeof(fixed_lipid);
	lipidLibrary = (lipid *)malloc( sizeof(lipid) * nlipid_space );
	
	for( int i = 0; i < nlipid_space; i++ )
	{
		lipidLibrary[i].c0 = library_fixed[i].c0;
		lipidLibrary[i].APL = library_fixed[i].APL;
		lipidLibrary[i].name = (char *)malloc( sizeof(char)* (1+strlen(library_fixed[i].name)) );
		strcpy( lipidLibrary[i].name, library_fixed[i].name );
		lipidLibrary[i].ns_atom = (char *)malloc( sizeof(char)* (1+strlen(library_fixed[i].ns_atom)) );
		strcpy( lipidLibrary[i].ns_atom, library_fixed[i].ns_atom );
	}
	nlipids = nlipid_space;

#if 0
	if( !inputFile )
	{	// load all POPC (eventually). now load it in as a composition.
		bilayerComp.nlipidTypes = 1;	
		bilayerComp.c0 = (double *)malloc( sizeof(double) );
		bilayerComp.APL = (double *)malloc( sizeof(double) );
		bilayerComp.num_C = (double *)malloc( sizeof(double) );
		bilayerComp.sum_C = (double *)malloc( sizeof(double) );
		memset( bilayerComp.num_C, 0, sizeof(double) );
		memset( bilayerComp.sum_C, 0, sizeof(double) );
		bilayerComp.names = (char **)malloc( sizeof(char *) );
		bilayerComp.names[0] = (char *)malloc( sizeof(char) * (1 + strlen(library_fixed[0].name) ) );
		strcpy( bilayerComp.names[0], library_fixed[0].name);
				
		for( int t = 0; t < nt; t++ )
		{
			theTriangles[t].composition.innerLeaflet = (double *)malloc( sizeof(double) * bilayerComp.nlipidTypes );
			theTriangles[t].composition.outerLeaflet = (double *)malloc( sizeof(double) * bilayerComp.nlipidTypes );
			memset( theTriangles[t].composition.innerLeaflet, 0, sizeof(double) * bilayerComp.nlipidTypes );
			memset( theTriangles[t].composition.outerLeaflet, 0, sizeof(double) * bilayerComp.nlipidTypes );
		}
	}
	else
	{
#endif
		char *buffer = (char *)malloc( sizeof(char) * 4096 );
	
	
		
		
		int nlipid_types_used = 0;
	
	
		for( int pass = 0; pass < 3; pass++ )
		{
			if( pass == 1 )
			{
				// the library is established.
				bilayerComp.nlipidTypes = 0;
				bilayerComp.c0  = (double *)malloc( sizeof(double) * nlipids ); 
				bilayerComp.APL = (double *)malloc( sizeof(double) * nlipids ); 
				bilayerComp.names = (char **)malloc( sizeof(char *) * nlipids );
		
				bilayerComp.num_C = (double *)malloc( sizeof(double) *nlipids );
				bilayerComp.sum_C = (double *)malloc( sizeof(double) *nlipids );
				memset( bilayerComp.num_C, 0, sizeof(double) *nlipids );
				memset( bilayerComp.sum_C, 0, sizeof(double) *nlipids );
			}
			if( pass == 2 )
			{
				printf("nlipidtypes: %d\n", bilayerComp.nlipidTypes );
				if( bilayerComp.nlipidTypes == 0  )
				{
					bilayerComp.nlipidTypes = 1;	
					bilayerComp.c0 = (double *)malloc( sizeof(double) );
					bilayerComp.APL = (double *)malloc( sizeof(double) );
					bilayerComp.names = (char **)malloc( sizeof(char *) );
					bilayerComp.names[0] = (char *)malloc( sizeof(char) * (1 + strlen(library_fixed[0].name) ) );
					strcpy( bilayerComp.names[0], library_fixed[0].name);
					bilayerComp.c0[0] = library_fixed[0].c0;
					bilayerComp.APL[0] = library_fixed[0].APL;
				
					bilayerComp.num_C = (double *)malloc( sizeof(double) );
					bilayerComp.sum_C = (double *)malloc( sizeof(double) );
					memset( bilayerComp.num_C, 0, sizeof(double)  );
					memset( bilayerComp.sum_C, 0, sizeof(double)  );
				}
				for( int t = 0; t < nt; t++ )
				{
					theTriangles[t].composition.innerLeaflet = (double *)malloc( sizeof(double) * bilayerComp.nlipidTypes );
					theTriangles[t].composition.outerLeaflet = (double *)malloc( sizeof(double) * bilayerComp.nlipidTypes );
					memset( theTriangles[t].composition.innerLeaflet, 0, sizeof(double) * bilayerComp.nlipidTypes );
					memset( theTriangles[t].composition.outerLeaflet, 0, sizeof(double) * bilayerComp.nlipidTypes );
				}
					
				bilayerComp.input_innerLeaflet = (double *)malloc( sizeof(double) * bilayerComp.nlipidTypes );
				bilayerComp.input_outerLeaflet = (double *)malloc( sizeof(double) * bilayerComp.nlipidTypes );
				memset( bilayerComp.input_innerLeaflet, 0, sizeof(double) * bilayerComp.nlipidTypes );
				memset( bilayerComp.input_outerLeaflet, 0, sizeof(double) * bilayerComp.nlipidTypes );
			}
			if( inputFile )
				rewind( inputFile );
			while( inputFile && !feof(inputFile) )
			{
				getLine( inputFile, buffer );
		
				if( feof(inputFile) ) break;
		
				if( !strncasecmp( buffer, "lipid ", strlen("lipid ") )  || !strncasecmp( buffer, "lipid\t", strlen("lipid\t") ) ) 
				{
					char junk[256];
					char command[256];
					int nr = sscanf( buffer, "lipid %s", command );
					if( !strcasecmp(command,"library") && pass == 0 )
					{
						char lipidName[256];
						double APL;
						double c0;
		
						char nsName[256];
						int nr = sscanf( buffer, "%s %s %s %lf %lf %s", 
							junk, command, lipidName, &APL, &c0, nsName );
						if( nr >= 5 )
						{
							int found = -1;

							for( int l = 0; l < nlipids; l++ )
							{
								if( !strcasecmp( lipidLibrary[l].name, lipidName ) )
									found = l;
							}
			
							if( found == -1 )
							{
								printf("Adding lipid %s with c0 %lf and area-per-lipid %lf to the library.\n", lipidName, c0, APL );
								if( nlipids == nlipid_space )
								{
									nlipid_space *= 2;
					
									lipidLibrary = (lipid *)realloc( lipidLibrary, sizeof(lipid) * nlipid_space );
								}								
		
								lipidLibrary[nlipids].name = (char *)malloc( sizeof(char) * ( 1 + strlen(lipidName)) );
								strcpy( lipidLibrary[nlipids].name, lipidName );
								lipidLibrary[nlipids].c0 = c0;
								lipidLibrary[nlipids].APL = APL;
								lipidLibrary[nlipids].ns_atom = NULL;
								found = nlipids;
								nlipids++;
							}
							else
							{
								printf("Replacing lipid %s with c0 %lf and area-per-lipid %lf in the library.\n", lipidName, c0, APL );
								lipidLibrary[found].c0 = c0;
								lipidLibrary[found].APL = APL;
							}

							if( nr >= 6 )
							{
								if( lipidLibrary[found].ns_atom ) 
									free(lipidLibrary[found].ns_atom );
								lipidLibrary[found].ns_atom = (char *)malloc( sizeof(char) * (1+strlen(nsName) ) );
								strcpy( lipidLibrary[found].ns_atom, nsName );
							}
						}	
						else
						{
							printf("Error reading library directive '%s'. Expected format 'lipid library lipidName APL c0'\n", buffer );
							exit(1);
						}
					}
					else if( (!strcasecmp( command, "inner" ) || !strcasecmp( command, "outer" ))  && pass > 0)
					{
						char lipidName[256];
						double parts;
						int nr = sscanf( buffer, "%s %s %s %lf", junk, junk, lipidName, &parts );
						// find this in the library.
				
						if( nr != 4 )
						{
							printf("Error reading command '%s'. Expected format 'lipid %s lipidName parts\n",
									command, command );
							exit(1);
						}
	
						int lpos = -1;		
						for( int lx = 0; lx < nlipids; lx++ )
						{
							if( !strcasecmp( lipidLibrary[lx].name, lipidName ) )
				 				lpos = lx;
						}
	
						if( lpos == -1 )
						{
							printf("Couldn't find lipid '%s' in the library.\n", lipidName );
		
							dumpLibrary();
						}
							
						int gotit = -1;
						for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
						{
							if( !strcasecmp( bilayerComp.names[x], lipidName ) )
								gotit = x;
						}	
	
						if( pass == 1 )
						{
	
							if( gotit == -1 )
							{
								int ind = bilayerComp.nlipidTypes;
								bilayerComp.names[ind] = (char *)malloc( sizeof(char) * ( 1 + strlen(lipidLibrary[lpos].name) ) );
								strcpy( bilayerComp.names[ind], lipidLibrary[lpos].name);
								bilayerComp.c0[ind] = lipidLibrary[lpos].c0;
								bilayerComp.APL[ind] = lipidLibrary[lpos].APL;
								bilayerComp.nlipidTypes++;
							}
						}
						else if( pass == 2 )
						{
							if( !strcasecmp( command, "inner" ) )
							{
								for( int t = 0; t < nt; t++ )
									theTriangles[t].composition.innerLeaflet[gotit] += parts;	
								bilayerComp.input_innerLeaflet[gotit] += parts;
							}
							else
							{
								for( int t = 0; t < nt; t++ )
									theTriangles[t].composition.outerLeaflet[gotit] += parts;	
								bilayerComp.input_outerLeaflet[gotit] += parts;
							}
						}
					}
				}
			}
		}
		free(buffer);
#if 0
	}
#endif
	for( int t = 0; t < nt; t++ )
	{
		double po=0,pi=0;
		for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
		{
			po += theTriangles[t].composition.outerLeaflet[x];
			pi += theTriangles[t].composition.innerLeaflet[x];
		}		

		if( po == 0 )
			theTriangles[t].composition.outerLeaflet[0] = 1;
		if( pi == 0 )
			theTriangles[t].composition.innerLeaflet[0] = 1;
		if( po == 0 )
			bilayerComp.input_outerLeaflet[0] = 1;
		if( pi == 0 )
			bilayerComp.input_innerLeaflet[0] = 1;
	}

	printf("--------------\n");
	printf("Inner leaflet:\n");
	double po=0,pi=0;
	int t = 0;
	for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
	{
		po += theTriangles[t].composition.outerLeaflet[x];
		pi += theTriangles[t].composition.innerLeaflet[x];
	}		
	
	for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
	{
		if( theTriangles[t].composition.innerLeaflet[x] > 0 )
			printf("%lf%% %s, c0 %lf APL %lf\n",
				100 * theTriangles[t].composition.innerLeaflet[x] / pi,
				bilayerComp.names[x], bilayerComp.c0[x], bilayerComp.APL[x] );
	}
	printf("Outer leaflet:\n");		
	for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
	{
		if( theTriangles[t].composition.outerLeaflet[x] > 0 )
			printf("%lf%% %s, c0 %lf APL %lf\n",
				100 * theTriangles[t].composition.outerLeaflet[x] / pi,
				bilayerComp.names[x], bilayerComp.c0[x], bilayerComp.APL[x] );
	}		
	printf("--------------\n");
	
}

void surface::getVPassData( double **verts, int *nvert, int **tris, int ** eft_tris, int *ntri, int **edges, int *nedges_in, int edge_dim )
{
	*nvert = nv;

	(*verts) = (double *)malloc( sizeof(double) * 3 * nv );

	for( int v = 0; v < nv; v++ )
	{
		(*verts)[3*v+0] = theVertices[v].r[0];
		(*verts)[3*v+1] = theVertices[v].r[1];
		(*verts)[3*v+2] = theVertices[v].r[2];
	}

	(*tris) = (int *)malloc( sizeof(int) * 3 * nt );
	(*eft_tris) = (int *)malloc( sizeof(int) * 3 * nt );

	*ntri = nt;

	for( int t = 0; t < nt; t++ )
	{
		(*tris)[3*t+0] = theTriangles[t].ids[0];
		(*tris)[3*t+1] = theTriangles[t].ids[1];
		(*tris)[3*t+2] = theTriangles[t].ids[2];

		(*eft_tris)[3*t+0] = theTriangles[t].edges[0];
		(*eft_tris)[3*t+1] = theTriangles[t].edges[1];
		(*eft_tris)[3*t+2] = theTriangles[t].edges[2];
	}

	*(nedges_in) = nedges;
	
	(*edges) = (int *)malloc( sizeof(int) * (2+edge_dim) * nedges );	

	for( int e = 0; e < nedges; e++ )
	{
		(*edges)[e*(2+edge_dim)+0] = theEdges[e].vertices[0];
		(*edges)[e*(2+edge_dim)+1] = theEdges[e].vertices[1];
		(*edges)[e*(2+edge_dim)+2] = theEdges[e].faces[0];
		(*edges)[e*(2+edge_dim)+3] = theEdges[e].faces[1];
	}
}

void surface::lipidSync(void)
{
#ifdef PARALLEL
	int *list = (int *)malloc( sizeof(int) * nt );
	double *fl = (double *)malloc( sizeof(double) * nt * 2 * bilayerComp.nlipidTypes);

	for( int p = 0; p < par_info.nprocs; p++ )
	{
		if( p == par_info.my_id )
		{
			for( int fx = 0; fx < par_info.nf[surface_id]; fx++ )
			{
				int f = par_info.faces[surface_id][fx];
				int t;
				if( f < nf_faces )
					t = theFormulas[f*nf_g_q_p].tri;
				else
					t = theIrregularFormulas[(f-nf_faces)*nf_irr_pts].tri;

				list[fx] = t;

				for( int lx = 0; lx < bilayerComp.nlipidTypes; lx++ )
				{
					fl[fx*bilayerComp.nlipidTypes*2+lx] = theTriangles[t].composition.innerLeaflet[lx];
					fl[fx*bilayerComp.nlipidTypes*2+bilayerComp.nlipidTypes+lx] = theTriangles[t].composition.outerLeaflet[lx];
				}
			}
		}

		int nvals = par_info.nf[surface_id];
		
		MPI_Bcast( &nvals, 1, MPI_INT, p, MPI_COMM_WORLD );		
		MPI_Bcast( list, nvals, MPI_INT, p, MPI_COMM_WORLD ); 
		MPI_Bcast( fl, nvals*bilayerComp.nlipidTypes*2, MPI_DOUBLE, p, MPI_COMM_WORLD ); 
	
		for( int tx = 0; tx < nvals; tx++ )
		{
			int t = list[tx];

			for( int lx = 0; lx < bilayerComp.nlipidTypes; lx++ )
			{
				theTriangles[t].composition.innerLeaflet[lx]	=	fl[tx*bilayerComp.nlipidTypes*2+lx];
				theTriangles[t].composition.outerLeaflet[lx]	=	fl[tx*bilayerComp.nlipidTypes*2+bilayerComp.nlipidTypes+lx];
			}
		}
	}
	free(list);
	free(fl);
#endif
}

void surface::lipidBroadcast(void)
{
#ifdef PARALLEL
	double *fl = (double *)malloc( sizeof(double) * nt * 2 * bilayerComp.nlipidTypes );

	if( par_info.my_id == BASE_TASK )
	{
		for( int t = 0; t < nt; t++ )
		{
			for( int lx = 0; lx < bilayerComp.nlipidTypes; lx++ )
			{	
				fl[t*2*bilayerComp.nlipidTypes+lx]                         = theTriangles[t].composition.innerLeaflet[lx];
				fl[t*2*bilayerComp.nlipidTypes+bilayerComp.nlipidTypes+lx] = theTriangles[t].composition.outerLeaflet[lx];
			}
		}	
	}

	ParallelBroadcast( fl, nt*2*bilayerComp.nlipidTypes );

	for( int t = 0; t < nt; t++ )
	{
		for( int lx = 0; lx < bilayerComp.nlipidTypes; lx++ )
		{	
			theTriangles[t].composition.innerLeaflet[lx]=fl[t*2*bilayerComp.nlipidTypes+lx];
			theTriangles[t].composition.outerLeaflet[lx]=fl[t*2*bilayerComp.nlipidTypes+bilayerComp.nlipidTypes+lx];
		}
	}
	
	free(fl);
#endif
}
