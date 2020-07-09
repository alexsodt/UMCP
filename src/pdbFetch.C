#include "pdb.h"
#include "dcd.h"

const char *libDir = "/Users/sodtaj/git_projects/HD/lib";

void pdbFetch( struct atom_rec **out_pdb, int *nout, const char *dir, const char *file )
{
	int lenLIB = strlen(libDir);
	int lenDir = strlen(dir);
	int lenFile = strlen(file);
	int lenTot = lenLIB + lenDir + lenFile + 256;
	char fileName[lenTot];
	sprintf(fileName, "%s/%s/%s.pdb", libDir, dir, file );

	FILE *theFile = fopen(fileName,"r");

	if( !theFile )
	{
		printf("Couldn't open file '%s' from pdbFetch.\n", fileName );
		exit(1);
	}

	loadPSFfromPDB( theFile );
	int nat = curNAtoms();

	(*out_pdb) = (struct atom_rec *)malloc( sizeof(struct atom_rec) * nat );
	*nout = nat;

	rewind(theFile);

	loadPDB( theFile, *out_pdb );	

	fclose(theFile);
}
