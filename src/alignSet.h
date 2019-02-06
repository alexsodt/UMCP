#ifndef __ALIGNSET__
#define __ALIGNSET__

void rotateArbitrary( double *xyz_rot, double *axis, double *origin, int nat, double val );
double alignStructuresOnAtomSet( double *xyz1_in, 
			int *atoms_align1,
			double *xyz2_in,
			int *atoms_align2,
			int nat,
			int nat_tot );
double alignStructuresOnAtomSetXY( double *xyz1_in, 
			int *atoms_align1,
			double *xyz2_in,
			int *atoms_align2,
			int nat,
			int nat_tot );
void displacePlanes( double *xyz_to_displace, int nat, double fac );
void rotatePlanar( double *xyz_rot, double pt[3], int nat, double val );
void flipCoords( double *xyz, int nat, double *porig, double *n );
void rotateAxial( double *xyz_rot, int at1, int at2, int nat, double val );
double getChi2( double *xyz1, 
			double *xyz2,
			int nat );

#endif
