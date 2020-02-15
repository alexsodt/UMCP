#ifndef __io_mol_readh__
#define __io_mol_readh__
#include "pdb.h"
int io_nframes( void);
void io_readStructureFile( char *fileName );
void io_readFrame( char *fileName );
struct atom_rec *io_getCurrentFrame();
int io_nNSAtoms( void );
int io_getNSAtom( int index );
void io_align( void );
void io_initialize_read( char *fileName );
#endif
