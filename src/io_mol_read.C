#define __io_mol_readc__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"
#include "library.h"
#include "io_mol_read.h"
#include "mutil.h"

// trajectory reading provided by either LOOS or some legacy functions I wrote -AJS
// header file is the same, regardless.

//void io_readStructureFile( char *fileName );
//void io_readFrame( char *fileName );

#include "pdb.h"

#ifdef HAVE_LIBLOOS
#include "loos.hpp"
#else
#include "dcd.h"
#endif

typedef struct
{
	FILE *file_handle;
	struct atom_rec *cur_frame;
	struct atom_rec *last_frame;
	int natoms;
	double cur_PBC[3][3];	
	int nframes;
	int populated;
	int read_init;
	int n_ns_atoms;
	int *ns_atoms;
	void align_on_ns();
} open_frame_file;

static open_frame_file *theFrame= NULL;

// only handles psf/pdb for structure and dcd/pdb for coords
// not particularly robust or well-tested.
// LOOS up to date and tested by broader community.

int io_nframes( void )
{
	if( !theFrame )
	{
		printf("Internal error: Read hasn't started for trajectory file yet.\n");
		exit(1);
	}

	return theFrame->nframes;
}

void io_readStructureFile( char *fileName )
{
	int natoms = 0;
#ifdef HAVE_LIBLOOS

#else 
	if( strlen(fileName) < 4 || (strcasecmp(fileName+strlen(fileName)-3,"pdb") && strcasecmp(fileName+strlen(fileName)-3,"psf")) )
	{	
		printf("File must have either pdb or psf extension.\n");
		printf("For a wide variety of I/O, build with LOOS.\n");
		exit(1);
	}

	FILE *theFile = fopen(fileName,"r");

	if( !strcasecmp(fileName+strlen(fileName)-3,"pdb") )
		loadPSFfromPDB(theFile);
	else
		loadPSF(theFile);
	natoms = curNAtoms();
	fclose(theFile);
#endif
	if( !theFrame )
	{
		theFrame = (open_frame_file *)malloc( sizeof(open_frame_file) );
		theFrame->file_handle = NULL;
		theFrame->natoms = natoms;
		theFrame->cur_frame = (struct atom_rec *)malloc( sizeof(struct atom_rec) * theFrame->natoms );
		theFrame->last_frame = NULL;
		
		theFrame->n_ns_atoms = 0;	
		theFrame->ns_atoms = NULL;

		theFrame->cur_PBC[0][0] = 0;
		theFrame->cur_PBC[0][1] = 0;
		theFrame->cur_PBC[0][2] = 0;
		
		theFrame->cur_PBC[1][0] = 0;
		theFrame->cur_PBC[1][1] = 0;
		theFrame->cur_PBC[1][2] = 0;
		
		theFrame->cur_PBC[2][0] = 0;
		theFrame->cur_PBC[2][1] = 0;
		theFrame->cur_PBC[2][2] = 0;

		theFrame->read_init = 0;
		theFrame->populated = 0;
		theFrame->nframes = 0;
	}	
}

void io_initialize_read( char *fileName )
{
#ifndef HAVE_LIBLOOS
	if( !theFrame->file_handle )
	{
		if( strlen(fileName) < 4 || (strcasecmp(fileName+strlen(fileName)-3,"pdb") && strcasecmp(fileName+strlen(fileName)-3,"dcd")) )
		{	
			printf("File must have either pdb or psf extension.\n");
			printf("For a wide variety of I/O, build with LOOS.\n");
			exit(1);
		}

		theFrame->file_handle = fopen(fileName,"r");

		if( !strcasecmp(fileName+strlen(fileName)-3,"pdb") )
			theFrame->nframes = 1;
		else
		{
			readDCDHeader( theFrame->file_handle );	
			setAligned();
			theFrame->nframes = curNFrames();
		}
	}
	theFrame->read_init = 1;
#endif	
}

void io_readFrame( char *fileName )
{
	if( !theFrame )
	{
		printf("Internal error: Cannot read a trajectory frame without structural information.\n");
		exit(1);
	}	
	
	if( !theFrame->read_init )
		io_initialize_read(fileName);

	int init = theFrame->populated;

	if( theFrame->populated )
	{	// if the last frame is populated, swap.
		if( !theFrame->last_frame )
			theFrame->last_frame = (struct atom_rec *)malloc( sizeof(struct atom_rec) * theFrame->natoms );
		else
		{
			for( int a = 0; a < theFrame->natoms; a++ )
				theFrame->last_frame[a].zap();
		}

		struct atom_rec *temp = theFrame->cur_frame;	
		theFrame->cur_frame = theFrame->last_frame;
		theFrame->last_frame = temp;	
	}

#ifdef HAVE_LIBLOOS

#else
	if( !strcasecmp(fileName+strlen(fileName)-3,"dcd") )
		loadFrame(theFrame->file_handle, theFrame->cur_frame);  
	else
		loadPDB(theFrame->file_handle, theFrame->cur_frame, theFrame->natoms);
	double La,Lb,Lc,alpha,beta,gamma;
	PBCD(&La,&Lb,&Lc,&alpha,&beta,&gamma);
	theFrame->cur_PBC[0][0] = La;
	theFrame->cur_PBC[1][1] = Lb;
	theFrame->cur_PBC[2][2] = Lc;
#endif

	if( !init )
	{
		int nlibrary = nLipidsInLibrary();
		int nspace = 10;
		theFrame->ns_atoms = (int *)malloc( sizeof(int) * nspace );
		int pres = -1;
		int added = 0;
		char pseg[256];
		pseg[0] = '\0';

		struct atom_rec *at = theFrame->cur_frame;

		for( int a = 0; a < theFrame->natoms; a++ )
		{
			if( (at[a].segid && !strcasecmp( pseg, at[a].segid )) || at[a].res != pres || !added )
			{
				if( at[a].res != pres )
					added = 0;

				for( int xl = 0; xl < nlibrary; xl++ )
				{
					if( !lipidLibrary[xl].ns_atom ) continue;

					if( !strcasecmp( lipidLibrary[xl].name, at[a].resname) &&
					    !strcasecmp( lipidLibrary[xl].ns_atom, at[a].atname) )
					{
						added = 1;
						if( theFrame->n_ns_atoms == nspace )
						{
							nspace *= 2;	
							theFrame->ns_atoms = (int *)realloc( theFrame->ns_atoms, sizeof(int) * nspace );
						}
						theFrame->ns_atoms[theFrame->n_ns_atoms] = a;
						theFrame->n_ns_atoms+=1;
					} 
				}
			}
		
			if( strlen(at[a].segid) < 255 )
				strcpy( pseg, at[a].segid );
			else
				pseg[0] = '\0';

			pres = at[a].res;		
		}
	}

	theFrame->populated = 1;
}

int io_nNSAtoms( void )
{
	if( !theFrame )
	{
		printf("Internal error: Read hasn't started for trajectory file yet.\n");
		exit(1);
	}

	return theFrame->n_ns_atoms;	
}

int io_getNSAtom( int index )
{
	if( index < 0 || index >= theFrame->n_ns_atoms )
	{
		printf("ERROR: requested atom index out of bounds.\n");
		exit(1);
	}

	return theFrame->ns_atoms[index];
}

struct atom_rec *io_getCurrentFrame( void )
{
	return theFrame->cur_frame;
}

// align on the neutral surface atoms.

void io_align( void )
{
	theFrame->align_on_ns();
}

void open_frame_file::align_on_ns( void )
{
	struct atom_rec *at = cur_frame;
	if( n_ns_atoms == 0 ) return;
	if( last_frame == NULL ) return;

	double com_prev[3] = { 0,0,0};
	for( int t = 0; t < n_ns_atoms; t++ )
	{
		com_prev[0] += at[ns_atoms[t]].x;
		com_prev[1] += at[ns_atoms[t]].y;
		com_prev[2] += at[ns_atoms[t]].z;
	}	

	com_prev[0] /= n_ns_atoms;
	com_prev[1] /= n_ns_atoms;
	com_prev[2] /= n_ns_atoms;

	double cur_com[3] = { 0,0,0};
	
	for( int t = 0; t < n_ns_atoms; t++ )
	{
		double dr[3] = { 
			at[ns_atoms[t]].x - last_frame[ns_atoms[t]].x,
			at[ns_atoms[t]].y - last_frame[ns_atoms[t]].y,
			at[ns_atoms[t]].z - last_frame[ns_atoms[t]].z };
		double put[3];
		MinImage3D( dr, cur_PBC, put );

		cur_com[0] += last_frame[ns_atoms[t]].x + dr[0];
		cur_com[1] += last_frame[ns_atoms[t]].y + dr[1];
		cur_com[2] += last_frame[ns_atoms[t]].z + dr[2];
	}	

	cur_com[0] /= n_ns_atoms;
	cur_com[1] /= n_ns_atoms;
	cur_com[2] /= n_ns_atoms;

	for( int t = 0; t < n_ns_atoms; t++ )
	{
		cur_frame[t].x -= (cur_com[0] - com_prev[0]);
		cur_frame[t].y -= (cur_com[1] - com_prev[1]);
		cur_frame[t].z -= (cur_com[2] - com_prev[2]);
	}
}



