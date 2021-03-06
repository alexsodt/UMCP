#include "interp.h"
#include "util.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

struct fixed_lipid
{	
	const char *name;
	double c0;
	double APL;
};

struct lipid
{
	char *name;
	double c0;
	double APL;
};

static lipid *lipidLibrary = NULL;
static int nlipid_space = 0;
static int nlipids = 0;

const fixed_lipid library_fixed[] =
{
	// Chen/Rand 1997
	// PO values copied from DO.
	{ "POPC", -1.0/87.3, 68.3 },
	{ "DOPE", -1.0/29.4, 63.4 },
	{ "POPE", -1.0/29.4, 58.8 },
	{ "DOPC", -1.0/87.3, 69.7 }, 
	{ "CHOL", -1.0/27.2, 31.0 }, 
	{ "DOG",  -1.0/11.5, 59.0 },
	{ "DOPS",  1.0/144.0, 71.6 } // Fuller 2003, area from charmm f.f.
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
	}
	nlipids = nlipid_space;

	if( !inputFile )
	{	// load all POPC (eventually). now load it in as a composition.
		bilayerComp.nlipidTypes = 1;	
		bilayerComp.c0 = (double *)malloc( sizeof(double) );
		bilayerComp.APL = (double *)malloc( sizeof(double) );
		bilayerComp.names = (char **)malloc( sizeof(char *) );
		bilayerComp.names[0] = (char *)malloc( sizeof(char) * (1 + strlen(library_fixed[0].name) ) );
		strcpy( bilayerComp.names[0], library_fixed[0].name);
	}
	else
	{
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
			}
			if( pass == 2 )
			{
				printf("nlipidtypes: %d\n", bilayerComp.nlipidTypes );
				if( bilayerComp.nlipidTypes == 0 )
				{
					bilayerComp.nlipidTypes = 1;	
					bilayerComp.c0 = (double *)malloc( sizeof(double) );
					bilayerComp.APL = (double *)malloc( sizeof(double) );
					bilayerComp.names = (char **)malloc( sizeof(char *) );
					bilayerComp.names[0] = (char *)malloc( sizeof(char) * (1 + strlen(library_fixed[0].name) ) );
					strcpy( bilayerComp.names[0], library_fixed[0].name);
					bilayerComp.c0[0] = library_fixed[0].c0;
					bilayerComp.APL[0] = library_fixed[0].APL;
				}
				for( int t = 0; t < nt; t++ )
				{
					theTriangles[t].composition.innerLeaflet = (double *)malloc( sizeof(double) * bilayerComp.nlipidTypes );
					theTriangles[t].composition.outerLeaflet = (double *)malloc( sizeof(double) * bilayerComp.nlipidTypes );
					memset( theTriangles[t].composition.innerLeaflet, 0, sizeof(double) * bilayerComp.nlipidTypes );
					memset( theTriangles[t].composition.outerLeaflet, 0, sizeof(double) * bilayerComp.nlipidTypes );
				}
			}
			rewind( inputFile );
			while( !feof(inputFile) )
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
		
						int nr = sscanf( buffer, "%s %s %s %lf %lf", 
							junk, command, lipidName, &APL, &c0 );
						if( nr == 5 )
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
							nlipids++;
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
							}
							else
							{
								for( int t = 0; t < nt; t++ )
									theTriangles[t].composition.outerLeaflet[gotit] += parts;	
							}
						}
					}
				}
			}
		}
		free(buffer);
	}

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
