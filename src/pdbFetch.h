#ifndef __pdbfetchh__
#define __pdbfetchh__

#include "pdb.h"
#include "dcd.h"

void pdbFetch( struct atom_rec **out_pdb, int *nout, const char *dir, const char *file );

#endif
