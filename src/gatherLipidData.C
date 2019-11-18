//#include <OpenCL/opencl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/time.h>
#include <math.h>
#include "util.h"
#include "interp.h"
#include "l-bfgs.h"
#include "m_triangles.h"
#include "pdb.h"
#include "dcd.h"
#include "clusterv.h"
#include "mutil.h"
#include "alignSet.h"
#include "amoeba_code/amoeba.h"

//#define DEBUG

	#define MAX_SHELL	10
#if !defined(PLOT_CHL) && !defined(GET_DIRECTORS) && !defined(GET_THICKNESS)
#define PLOT_CHL
//#define GET_DIRECTORS
#endif

#define DO_EXPAND_XY
#define FIXED_SEED
#define FORCE_PASS 	0
#define TOLERANT_PASS 	1
#define RESTRICTED_PASS 2
#define RANDOM_PASS	3
//#define SKIP_OPT

//#define APPLY_EXTERNAL


	double *M[3];
	
	int mlow = 5;
	int mhigh = 7;


double CMIN = -0.05;
double CMAX =  0.05;
#define N_CDIST_BINS 100
typedef struct lipid_record {
	const char *resname;
	const char *surfaceName;
	const char *tailName1;
	const char *tailName2;
	double c0;
	double cdist[N_CDIST_BINS];
	double sum_c_obs_upper;
	double n_obs_upper;
	double sum_c_obs_lower;
	double n_obs_lower;
} lipid_record;

typedef struct process_record 
{
	int rec;
	int atn;
	int res;
	int surface;
	int spot;
} process_record;

lipid_record lipids[] =
{
	{"DOPE", "C21", "C218", "C318", -1.0/30.0  },
	{"DOPS", "C21", "C218", "C318",  1.0/140.0 },
	{"DOPC", "C21", "C218", "C318", -1.0/200.0 }
};
	





surface *minSurface = NULL;
double loose_normal_tol = -0.6;
double normal_tol = -0.2;
double thick = 13.8;
extern int do_print_area_ps;
extern double av_vol_err;
extern double KV_scale;
void fdiff_check( double *coords );
int print_energy = 0;
extern double tilt_scale;
surface *upperSurface;
surface *lowerSurface;
volumes *theVolumes; 
 void amoeba(double **p, double y[], int ndim, double ftol,
               double (*funk)(double []), int *nfunk);
void amoebaOpt( double *r, int nv, double *rc, int nrc );
double surface_grad( double *r, double *g );
double surface_e( double *r );
	void curvatureGrad( double *g, surface *surface1, surface *surface2, double *r, double thick );
	double curvatureEnergy( surface *surface1, surface *surface2, double *r, double thick);
	void tiltGrad( double *g, surface *surface1, surface *surface2, double *r );
	void tiltGradExt( double *g, surface *surface1, surface *surface2, double *r );
	double tiltEnergy( surface *surface1, surface *surface2, double *r);
	double tiltEnergyExt( surface *surface1, surface *surface2, double *r, double *cons_energy);
double f(double *r);
double fdf( double *r, double *g);
void writeDualSurface( const char *fileNameXYZ, const char *fileNamePSF, surface *upperSurface, surface *lowerSurface, double *coords, double alpha );
void label_points( int *edges, int *edges_for_point, int *nedges_for_point, int *tris, int xp1, int xp2, int xp3, int KMED, int *labels, int *pt_elim );

int isTriangle( int *tri, int ntri, int *tri_active, int p1, int p2, int p3 );
int extract( double *rp, int nuse, int *bonds, int *nbonds, int max_bonds, double rcut, double PBC_vec[3][3], int do_write_struct );
int fix_alpha = 0;
static int max_edges = 20;

int swapEdge( int e, int *edges, int *edges_for_point, int *nedges_for_point, int *tri, int *eft_tris, int loose );
void cleanEdges( int v, int *edges, int *edges_for_point, int *nedges_for_point, int *tri, int *eft_tris, int nedges, int *pt_elim  );
int hardDebugEdge( int e, int *edges, int *edges_for_point, int *nedges_for_point, int *tri, int *eft_tris, int loose );

double opt_shift( surface *surface1, surface *surface2, double *pt_set, int npts, double La, double Lb, double Lc, double shift[3], int just_report );
	
struct region_site
{
	double curve_sampled;
	double curve_sampled_c;
	double curve_sampled_l;
	double chl_orient[3];
	double lip_orient[3];
	double tilt_sampled;
	double tilt_sampled_c;
	double tilt_sampled_l;
	double h_sampled;
	double h_sampled_c;
	double h_sampled_l;
	double *n_chl_b;
	double *n_lip_b;
	double n_chl;
	double n_lip;
	double phi_err;
	double phi;
	double R;
	double theta;
	int nblks;

	void setNBlks( int nblks );
	void clear(void);
	void computePhi(void);
void  addChl( double c, double h, double tilt, double R, double ltheta, int blk, double tx, double ty, double tz);
void  addLip( double c, double h, double tilt, double R, double ltheta, int blk, double tx, double ty, double tz );
};

void region_site::computePhi( void )
{
	double av = 0;
	double av2 = 0;

	for( int b = 0; b < nblks; b++ )
	{
		double lav = (n_chl_b[b] ) / ( n_chl_b[b] + n_lip_b[b] + 1e-10 );
		printf("blk: %d lav: %lf\n", b, lav );
		av += lav;
		av2 += lav*lav; 
	}
	
	av /= nblks;
	av2 /= nblks;

	phi = av;
	phi_err = sqrt(av2-av*av) / sqrt(nblks);
}

void region_site::clear( void )
{
	chl_orient[0] = 0;
	chl_orient[1] = 0;
	chl_orient[2] = 0;
	
	lip_orient[0] = 0;
	lip_orient[1] = 0;
	lip_orient[2] = 0;

	tilt_sampled   = 0;
	tilt_sampled_c = 0;
	tilt_sampled_l = 0;
	curve_sampled   = 0;
	curve_sampled_c = 0;
	curve_sampled_l = 0;
	h_sampled   = 0;
	h_sampled_c = 0;
	h_sampled_l = 0;
	n_chl = 0;
	n_lip = 0;
	R = 0;
	theta = 0;
	nblks = 0; 
}

void region_site::setNBlks( int nb_in )
{
	nblks = nb_in;
	n_chl_b = (double *)malloc( sizeof(double) * nb_in );
	n_lip_b = (double *)malloc( sizeof(double) * nb_in );

	memset( n_chl_b, 0, sizeof(double) * nb_in );	
	memset( n_lip_b, 0, sizeof(double) * nb_in );	
}

void region_site::addChl( double c, double h, double t, double lR, double ltheta, int blk, double tx, double ty, double tz)
{
	if( n_lip+n_chl > 0 )
	{
		double cur_theta =theta / (n_lip+n_chl);
		while( ltheta -cur_theta > M_PI )
			ltheta -= 2*M_PI;
		while( ltheta -cur_theta < -M_PI )
			ltheta += 2*M_PI;
	} 

	chl_orient[0] += tx;
	chl_orient[1] += ty;
	chl_orient[2] += tz;

	theta += ltheta;
	tilt_sampled += t;
	tilt_sampled_c += t;
	curve_sampled += c;
	curve_sampled_c += c;
	h_sampled += h;
	h_sampled_c += h;
	n_chl += 1;
	R += lR;

	n_chl_b[blk] += 1;
}

void  region_site::addLip( double c, double h, double t, double lR, double ltheta, int blk, double tx, double ty, double tz)
{
	if( n_lip+n_chl > 0 )
	{
		double cur_theta =theta / (n_lip+n_chl);
		while( ltheta -cur_theta > M_PI )
			ltheta -= 2*M_PI;
		while( ltheta -cur_theta <- M_PI )
			ltheta += 2*M_PI;
	} 
	
	lip_orient[0] += tx;
	lip_orient[1] += ty;
	lip_orient[2] += tz;

	theta += ltheta;
	curve_sampled += c;
	curve_sampled_l += c;
	tilt_sampled += t;
	tilt_sampled_l += t;
	h_sampled += h;
	h_sampled_l += h;
	n_lip += 1;
	R += lR;
	
	n_lip_b[blk] += 1;
}

int main( int argc, char **argv )
{
	double bin_width = 3.0;

#ifdef FIXED_SEED
	srand(1);
#else
          struct timeval tp;
 
          gettimeofday( &tp, NULL );
 
          srand((long)tp.tv_usec);
#endif

	char buffer[4096];

	/*
		write out upper rho, lower rho, midplane rho, protein rho
		// dimensions:
	*/

	if( argc < 7 )
	{
		printf("Syntax: gatherLipidData.exe psf align.pdb upper.surface midplane.surface lower.surface dcd1 [dcd2 ...]\n");
		return -1;
	}
	       
	FILE *psfFile = fopen(argv[1],"r");
	
	if( !psfFile )
	{
		printf("Couldn't open PSF file '%s'.\n", argv[1] );
		exit(1);
	}

	if( !strcasecmp( argv[1] + strlen(argv[1])-3, "pdb" ) ) 
		loadPSFfromPDB( psfFile );    
        else
		loadPSF( psfFile );

	int nat = curNAtoms();

	struct atom_rec *align_at = (struct atom_rec *)malloc( sizeof(atom_rec) * nat );
	FILE *alignFile = fopen(argv[2],"r");

	if( !alignFile )
	{
		printf("Could not open alignment file '%s'.\n", argv[2] );
		return -1;
	}
	int nat_align = loadPDB( alignFile, align_at, nat );
	int nuse = 0;
	FILE *dcdFile = fopen( argv[6], "r");
	
	if( !dcdFile )
	{
		printf("Couldn't open dcd file '%s'.\n", argv[6] );
		exit(1);
	}

        readDCDHeader(dcdFile);
	setAligned();

	int nframes = curNFrames();
	double nframes_tot = 0;
	struct atom_rec *at = (struct atom_rec*)malloc( sizeof(struct atom_rec) * nat );
	
	loadFrame( dcdFile, at );

	double last_PBC[3];
	double cur_PBC[3];
	double La,Lb,Lc,alpha,beta,gamma;	
	PBCD( &La, &Lb, &Lc, &alpha, &beta, &gamma );

	if( !DCDsuccess() )
	{
		printf("Could not load intended frames\n");
		exit(1);
	}

	fclose(dcdFile);
	
	int *main_set  = (int *)malloc( sizeof(int) * nat );
	int nalign = 0;
	
	int *align_set = (int *)malloc( sizeof(int) * nat_align );
	for( int am = 0; am < nat; am++ )
	{
		if( strncasecmp( at[am].segid, "PRO", 3 ) && strncasecmp( at[am].atname, "BB", 2 ) && strncasecmp( at[am].atname, "SC", 2 ) )
			continue;
		
		int gotit = -1;
		for( int aa = 0; aa < nat_align && gotit == -1; aa++ )
		{
			
			if( at[am].res == align_at[aa].res &&
			   !strcasecmp( at[am].atname, align_at[aa].atname) && 
			   !strcasecmp( at[am].resname, align_at[aa].resname) )
			{
				gotit = aa;	
			}		
		}

		if( gotit >= 0 )
		{
			align_set[nalign] = gotit;
			main_set[nalign] = am;
			nalign++;
		}
	}

	cur_PBC[0] = La;
	cur_PBC[1] = Lb;
	cur_PBC[2] = Lc;

	double *prot_align = (double*)malloc( sizeof(double) * curNAtoms() * 3 );

	loadFrame( dcdFile, at );

	int *prot_atoms = (int *)malloc( sizeof(int) * curNAtoms() );
	int np=0;
	


	double *align_array = (double *)malloc( sizeof(double) * nalign * 3 );
	int *align_set_ind = (int *)malloc( sizeof(int) * nalign );
	for( int c = 0; c < nalign; c++ )
	{
		align_set_ind[c] = c;
		align_array[3*c+0] = align_at[align_set[c]].x;
		align_array[3*c+1] = align_at[align_set[c]].y;
		align_array[3*c+2] = align_at[align_set[c]].z;
	}
	
	// remove center from alignment.

	double acom[3] = { 0,0,0};
	for( int c = 0; c < nalign; c++ )
	{
		acom[0] += align_array[3*c+0];
		acom[1] += align_array[3*c+1];
		acom[2] += align_array[3*c+2];
	}

	acom[0] /= nalign;
	acom[1] /= nalign;
	acom[2] /= nalign;

	for( int a = 0; a < nat; a++ )
		at[a].zap();
			
	surface *upperSurface = (surface *)malloc( sizeof(surface) );
	surface *midplaneSurface = (surface *)malloc( sizeof(surface) );
	surface *lowerSurface = (surface *)malloc( sizeof(surface) );

	upperSurface->loadLattice( argv[3], 0. );
	midplaneSurface->loadLattice( argv[4], 0. );
	lowerSurface->loadLattice( argv[5], 0. );

	surface *theSurfaces[3] = { upperSurface, midplaneSurface, lowerSurface };

#if 0
	for( int ss = 0; ss < 3; ss++ )
	{
		surface *mod = theSurfaces[ss];
		double scomz=0;
		double del = 0;	
		
		for( int v = 0; v < mod->nv; v++ )
			scomz += mod->theVertices[v].r[2] / mod->nv;
		while( scomz+del - acom[2] > mod->PBC_vec[2][2]/2 )
			del -= mod->PBC_vec[2][2];
		while( scomz+del - acom[2] < -mod->PBC_vec[2][2]/2 )
			del += mod->PBC_vec[2][2];
		for( int v = 0; v < mod->nv; v++ )
			mod->theVertices[v].r[2] += del;
	}
#endif

	
	double *ru = (double *)malloc(sizeof(double) * (upperSurface->nv*3+3));
	double *rm = (double *)malloc(sizeof(double) * (midplaneSurface->nv*3+3));
	double *rl = (double *)malloc(sizeof(double) * (lowerSurface->nv*3+3));
	upperSurface->get(ru);
	midplaneSurface->get(rm);
	lowerSurface->get(rl);

	ru[3*upperSurface->nv] = 1.0;
	ru[3*upperSurface->nv+1] = 1.0;
	ru[3*upperSurface->nv+2] = 1.0;
	rm[3*upperSurface->nv] = 1.0;
	rm[3*upperSurface->nv+1] = 1.0;
	rm[3*upperSurface->nv+2] = 1.0;
	rl[3*upperSurface->nv] = 1.0;
	rl[3*upperSurface->nv+1] = 1.0;
	rl[3*upperSurface->nv+2] = 1.0;

	upperSurface->generatePlan();
	midplaneSurface->generatePlan();
	lowerSurface->generatePlan();
	
	upperSurface->box_system();
	lowerSurface->box_system();
	
	// this has all the coordinates we will need for processing. they all get aligned however we wish.

	double *main_arr = NULL;
	// number of frames we've processed.
	int ntraj = 0;

	int nproc = 0;
	int nproc_size = 10;
	process_record *proc = (process_record *)malloc( sizeof(process_record) * nproc_size ); 
				
	int cntr = nalign;
	
	void generateSubdivisionMatrices( double **M, int mlow, int mhigh );

	double *M5 = (double *)malloc( sizeof(double) * 4 * 11 * 12 ); 
	double *M6 = (double *)malloc( sizeof(double) * 4 * 12 * 12 ); 
	double *M7 = (double *)malloc( sizeof(double) * 4 * 13 * 13 );
	
	M[0] = M5;
	M[1] = M6;
	M[2] = M7;
	upperSurface->generateSubdivisionMatrices( M, mlow, mhigh );


#ifdef DEBUG
	FILE *debugFile = fopen("debug.vmd", "w");
#endif

	int n_pt_set = 0;
	double *pt_set = NULL;
	int *get_set = NULL;
	

	double sign_flip_upper = 1;
	double sign_flip_lower = 1;


	double r_test[3], n_test[3];
	upperSurface->evaluateRNRM( 0, 0.33, 0.33, r_test, n_test, ru );
	if( n_test[2] < 0 )
		sign_flip_upper = -1;
	lowerSurface->evaluateRNRM( 0, 0.33, 0.33, r_test, n_test, rl );
	if( n_test[2] > 0 )
		sign_flip_lower = -1;

	printf("sign_flips: %lf %lf\n", sign_flip_lower, sign_flip_upper );
				
	int nrecs = sizeof(lipids)/sizeof(lipid_record);

	for( int r = 0; r < nrecs; r++ )
	{
		memset( lipids[r].cdist, 0, sizeof(double) * N_CDIST_BINS );
		lipids[r].sum_c_obs_upper = 0;
		lipids[r].sum_c_obs_lower = 0;
		lipids[r].n_obs_upper = 0;	
		lipids[r].n_obs_lower = 0;	
	}

	for( int c = 6; c < argc; c++ )
	{
		FILE *dcdFile = fopen(argv[c], "r");

        	readDCDHeader(dcdFile);
		setAligned();

		int nframes = curNFrames();
		
		int trip = 1;
		for( int f = 0; f < nframes; f++, nframes_tot++)
		{
			loadFrame( dcdFile, at );

			double La,Lb,Lc,alpha,beta,gamma;	
			PBCD( &La, &Lb, &Lc, &alpha, &beta, &gamma );
			
			if( !DCDsuccess() )
			{
				printf("Could not load intended frames.\n");
				exit(1);
			}

			if( trip )
			{
				// set up arrays to process.
				for( int a = 0; a < curNAtoms(); a++ )
				{
					int rec = -1;
					int create = 0;
					int do_surface = 0;
	
					for( int r = 0; r < nrecs; r++ )
					{
						if( !strcasecmp( at[a].resname, lipids[r].resname ) )
						{
							if( !strcasecmp( at[a].atname, lipids[r].surfaceName ) )
							{
								create = 1;
								do_surface = 1;
								rec = r;
							}
							if( !strcasecmp( at[a].atname, lipids[r].tailName1 ) || !strcasecmp( at[a].atname, lipids[r].tailName2 ) )
							{
								create = 1;
								rec = r;
							}
						}
					}

					if( create )
					{
						if( nproc == nproc_size )
						{
							nproc_size *= 2;
							proc = (process_record *)realloc( proc, sizeof(process_record) * nproc_size );
						}

						proc[nproc].atn = a;
						proc[nproc].res = at[a].res;
						proc[nproc].rec = rec;
						proc[nproc].surface = do_surface;
						proc[nproc].spot = cntr;

						// store all 9 X/Y PBC images.
						cntr += 9;

						nproc++;
					}
				}
				
				int frac = 10;
				pt_set = (double *)malloc( sizeof(double) * 3 * (nproc/(double)frac +1) );
				get_set = (int *)malloc( sizeof(int) * (nproc/(double)frac+1) );
				int off = 0;

				for( int x = 0; x < nproc; x += frac )
				{
					if( proc[x].surface )
					{
						get_set[off] = proc[x].spot + 4; // the central pbc image.	
						off++;
					}
				}
				n_pt_set = off;

				main_arr = (double *)malloc( sizeof(double) * 3 * cntr ); 
			}

			// remove protein cross-PBC discontinuities:
#if 1
			// center in z around the protein.
			double align_com[3] = { 0, 0, 0 };
			double cur_align[3] = { 0,0,0};
			char pseg[256];
			int pres = -1;
			int pactive =0;
			for( int a = 0; a < curNAtoms(); a++ )
			{
				if( !strcasecmp( at[a].segid, pseg) || ((at[a].res != pres) && !( (!strncasecmp(at[a].segid, "PRO",3) || !strcasecmp( at[a].atname, "BB") || !strncasecmp( at[a].atname, "SC",2)  )  && pactive))  ) 
				{
					cur_align[0] = align_com[0];
					cur_align[1] = align_com[1];
					cur_align[2] = align_com[2];
				}
	
				while( at[a].x - cur_align[0] < -La/2 ) at[a].x += La;
				while( at[a].y - cur_align[1] < -Lb/2 ) at[a].y += Lb;
				while( at[a].z - cur_align[2] < -Lc/2 ) at[a].z += Lc;
				while( at[a].x - cur_align[0] >  La/2 ) at[a].x -= La;
				while( at[a].y - cur_align[1] >  Lb/2 ) at[a].y -= Lb;
				while( at[a].z - cur_align[2] >  Lc/2 ) at[a].z -= Lc;
	
				if( !strncasecmp( at[a].segid, "PRO",3) || !strcasecmp( at[a].atname, "BB") || !strncasecmp( at[a].atname, "SC", 2) ) 
					pactive = 1;
				else
					pactive = 0;

				cur_align[0] = at[a].x;
				cur_align[1] = at[a].y;
				cur_align[2] = at[a].z;
				pres = at[a].res;
				strcpy( pseg, at[a].segid );
			}
#endif




			align_com[0] = 0;
			align_com[1] = 0;
			align_com[2] = 0;

			for( int x = 0; x < nalign; x++ )
			{	
				align_com[0] += at[main_set[x]].x;
				align_com[1] += at[main_set[x]].y;
				align_com[2] += at[main_set[x]].z;
			}
	
			align_com[0] /= nalign;	
			align_com[1] /= nalign;	
			align_com[2] /= nalign;	

			int curp = 0;
			for( int x = 0; x < nalign; x++ )
			{
				main_arr[curp*3+0] = at[main_set[x]].x;
				main_arr[curp*3+1] = at[main_set[x]].y;
				main_arr[curp*3+2] = at[main_set[x]].z;
				curp++;
			}

			for( int p = 0; p < nproc; p++ )
			{
				int pbc_spot = 0;
				for( int dx = -1; dx <= 1; dx++ )
				for( int dy = -1; dy <= 1; dy++, pbc_spot++ )
				{
					main_arr[3*(proc[p].spot+pbc_spot)+0] = at[proc[p].atn].x + dx * La;
					main_arr[3*(proc[p].spot+pbc_spot)+1] = at[proc[p].atn].y + dy * Lb;
					main_arr[3*(proc[p].spot+pbc_spot)+2] = at[proc[p].atn].z;
				}
			}

			double rmsd = alignStructuresOnAtomSetXY( align_array, align_set_ind, main_arr, align_set_ind, nalign, cntr );	

			double pcen[3] = { 0,0,0};

			for( int x = 0; x < nalign; x++ )
			{
				pcen[0] += main_arr[3*x+0];
				pcen[1] += main_arr[3*x+1];
				pcen[2] += main_arr[3*x+2];
			}
			pcen[0] /= nalign;
			pcen[1] /= nalign;
			pcen[2] /= nalign;

//			printf("pcen: %le %le %le\n", pcen[0], pcen[1], pcen[2] );
			// wrap z to a midplane atom near the protein.


			double wrapto = 0;
/*			double minr2 = 1e10;

			for( int p = 0; p < nproc; p++ )
			{
				if( proc[p].surface ) continue;

				double dr[3] = { 
					main_arr[3*(proc[p].spot+4)+0] - pcen[0],	
					main_arr[3*(proc[p].spot+4)+1] - pcen[1],	
					main_arr[3*(proc[p].spot+4)+2] - pcen[2] };
				while( dr[0] < -La/2 ) dr[0] += La;	
				while( dr[0] >  La/2 ) dr[0] -= La;	
				while( dr[1] < -Lb/2 ) dr[1] += Lb;	
				while( dr[1] >  Lb/2 ) dr[1] -= Lb;	
				while( dr[2] < -Lc/2 ) dr[2] += Lc;	
				while( dr[2] >  Lc/2 ) dr[2] -= Lc;	
				
				double r2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];

				if( r2 < minr2 )
				{
					minr2 = r2;
					wrapto = dr[2] + pcen[2];
				}
			}
*/

#define N_BINS_MOLDIST 100

                         double best_chi2 = 1e10;
 
                         int nbins = N_BINS_MOLDIST;
                         double moldist[N_BINS_MOLDIST];
                         memset( moldist, 0, sizeof(double) * N_BINS_MOLDIST );

			for( int p = 0; p < nproc; p++ )
			{
				if( proc[p].surface ) continue;

                                      double tz = main_arr[proc[p].spot*3+2];

                                      while( tz < 0 ) tz += Lc;
                                      while( tz >= Lc ) tz -= Lc;

                                      int zb = N_BINS_MOLDIST * tz / Lc; // this is right
                                      if( zb < 0 ) zb = 0;
                                      if( zb >= N_BINS_MOLDIST ) zb = N_BINS_MOLDIST-1;
                                      moldist[zb] += 1;
			}

                         for( int zb = 0; zb < nbins; zb++ )
                         {
                                 double zv = Lc * (zb+0.5) / (double)N_BINS_MOLDIST;
 
                                  int zlow  = zb- nbins/2;
                                  int zhigh = zlow + nbins;
 
                                  double lchi2 = 0;
                                  for( int iz = zlow; iz < zhigh; iz++ )
                                  {
                                          double dz = (iz+0.5) / nbins - zv;
 
                                          int iiz = iz;
                                          while( iiz < 0 ) iiz += nbins;
                                          while( iiz >= nbins ) iiz -= nbins;
 
                                          lchi2 += moldist[iiz] * (dz) * (dz);
                                  }
 
                                  if( lchi2 < best_chi2 )
                                  {
                                          best_chi2 = lchi2;
                                          wrapto = zv;
                                  }
                         }

			while( wrapto - pcen[2] < -Lc/2 ) wrapto += Lc;
			while( wrapto - pcen[2] >  Lc/2 ) wrapto -= Lc;
		
//			printf("wrapto: %le\n", wrapto );

			double lipid_com_z = 0;

			for( int x = nalign; x < cntr; x++ )
			{
				while( main_arr[3*x+2] - wrapto < -Lc/2 )
					main_arr[3*x+2] += Lc;
				while( main_arr[3*x+2] - wrapto >  Lc/2 )
					main_arr[3*x+2] -= Lc;


				lipid_com_z += main_arr[3*x+2];
			}

			lipid_com_z /= (1e-15 + (cntr-nalign) );

			// match to midplane cntrl pt z com.
	
			double surface_midz = 0;
			for( int x = 0; x < midplaneSurface->nv; x++ )
				surface_midz += rm[3*x+2];
			surface_midz /= midplaneSurface->nv;

//			printf("surface mid %le lipid_com_z %le\n", surface_midz, lipid_com_z );

			for( int x = 0; x < cntr; x++ )
				main_arr[3*x+2] += (surface_midz - lipid_com_z);

/*
			for( int x = 0; x < curp; x++ )
			{
				main_arr[3*x+0] -= pcen[0];
				main_arr[3*x+1] -= pcen[1];
			}
*/

			for( int x = 0; x < n_pt_set; x++ )
			{
				pt_set[3*x+0] = main_arr[3*get_set[x]+0];
				pt_set[3*x+1] = main_arr[3*get_set[x]+1];
				pt_set[3*x+2] = main_arr[3*get_set[x]+2];
			}

			double out_shift[3] = {0,0,0};
//			double chi2 = 0;
			double chi2 = opt_shift( upperSurface, lowerSurface, pt_set, n_pt_set, La, Lb, Lc, out_shift, 1 /* just report */ );
//			printf("out_shift: %lf %lf %lf\n", out_shift[0], out_shift[1], out_shift[2] );
	
			for( int x = 0; x < cntr; x++ )
			{
				main_arr[3*x+0] += out_shift[0];
				main_arr[3*x+1] += out_shift[1];
				main_arr[3*x+2] += out_shift[2];
			}
			for( int p = 0; p < nproc; p++ )
			{
//				if( at[proc[p].atn].res != 1176 ) continue;
				if( proc[p].surface )
				{
				int pbc_spot = 0;

				double close_lower = 1e10;
				double close_upper = 1e10;
				double c_at_u = 0;
				double c_at_l = 0;
				int atn = proc[p].atn;
				double close_point_lower[3] = { 0,0,0};	
				double close_point_upper[3] = { 0,0,0};	
				
				double source_point_lower[3] = { 0,0,0};	
				double source_point_upper[3] = { 0,0,0};	

//				for( int dx = -1; dx <= 1; dx++ )
//				for( int dy = -1; dy <= 1; dy++, pbc_spot++ )
				int dx = 0;
				int dy = 0;
				pbc_spot = 4;
				{
					double *rp = main_arr+(proc[p].spot+pbc_spot)*3;
			
					if( p == 6 && dx == 0 && dy == -1 )
					{
						
					}
					
			//		void surface::nearPointOnSurface( double *pt, int *col_f, double *col_u, double *col_v, double **M, int mlow, int mhigh, double *distance )
						double dist_u;
						int f_u;
						double col_u_u, col_v_u;
						//upperSurface->nearPointOnSurface( rp, &f_u, &col_u_u, &col_v_u, M, mlow, mhigh, &dist_u );					
						upperSurface->nearPointOnBoxedSurface( rp, &f_u, &col_u_u, &col_v_u, M, mlow, mhigh, &dist_u );					
						double pt_upper[3], nrm_upper[3];

						upperSurface->evaluateRNRM( f_u, col_u_u, col_v_u, pt_upper, nrm_upper, ru );
						double dist_l;
						int f_l;
						double col_u_l, col_v_l;
						lowerSurface->nearPointOnBoxedSurface( rp, &f_l, &col_u_l, &col_v_l, M, mlow, mhigh, &dist_l );					

						double pt_lower[3], nrm_lower[3];

						lowerSurface->evaluateRNRM( f_l, col_u_l, col_v_l, pt_lower, nrm_lower, rl );


						double dr_u[3] = { rp[0] - pt_upper[0], rp[1] - pt_upper[1], rp[2] - pt_upper[2] };
						while( dr_u[0] < -La/2 ) dr_u[0] += La;
						while( dr_u[1] < -Lb/2 ) dr_u[1] += Lb;
						while( dr_u[2] < -Lc/2 ) dr_u[2] += Lc;
						while( dr_u[0] > La/2 ) dr_u[0] -= La;
						while( dr_u[1] > Lb/2 ) dr_u[1] -= Lb;
						while( dr_u[2] > Lc/2 ) dr_u[2] -= Lc;
						double dr_l[3] = { rp[0] - pt_lower[0], rp[1] - pt_lower[1], rp[2] - pt_lower[2] };
						while( dr_l[0] < -La/2 ) dr_l[0] += La;
						while( dr_l[1] < -Lb/2 ) dr_l[1] += Lb;
						while( dr_l[2] < -Lc/2 ) dr_l[2] += Lc;
						while( dr_l[0] > La/2 ) dr_l[0] -= La;
						while( dr_l[1] > Lb/2 ) dr_l[1] -= Lb;
						while( dr_l[2] > Lc/2 ) dr_l[2] -= Lc;
					
						double r2_u = sqrt(dr_u[0]*dr_u[0] + dr_u[1]*dr_u[1] +dr_u[2]*dr_u[2]);
						double r2_l = sqrt(dr_l[0]*dr_l[0] + dr_l[1]*dr_l[1] +dr_l[2]*dr_l[2]);
	
/*						printf("at %s res %d name %s pt %le %le %le upper: %le (%le %le %le) lower: %le (%le %le %le)\n", 
						at[atn].atname, at[atn].res, at[atn].resname, rp[0], rp[1], rp[2], 
						dist_u, pt_upper[0], pt_upper[1], pt_upper[2],
						dist_l, pt_lower[0], pt_lower[1], pt_lower[2] );
*/

						if( r2_u < close_upper )
						{
//							printf("CLOSE UPPER %le\n", r2_u );
							close_upper = r2_u;
							double k;
							c_at_u = sign_flip_upper * upperSurface->c( f_u, col_u_u, col_v_u, ru, &k );
							close_point_upper[0] = pt_upper[0];
							close_point_upper[1] = pt_upper[1];
							close_point_upper[2] = pt_upper[2];
							source_point_upper[0] = rp[0];
							source_point_upper[1] = rp[1];
							source_point_upper[2] = rp[2];
						}
						if( r2_l < close_lower )
						{
							close_lower = r2_l;
							c_at_l =  sign_flip_lower * lowerSurface->c( f_l, col_u_l, col_v_l, rl );
							close_point_lower[0] = pt_lower[0];
							close_point_lower[1] = pt_lower[1];
							close_point_lower[2] = pt_lower[2];
							source_point_lower[0] = rp[0];
							source_point_lower[1] = rp[1];
							source_point_lower[2] = rp[2];
						}
				}
		

				if( close_upper < close_lower )
				{
					lipids[proc[p].rec].sum_c_obs_upper += c_at_u;
					lipids[proc[p].rec].n_obs_upper += 1;
					int cbin = N_CDIST_BINS * (c_at_u-CMIN)/(CMAX-CMIN);
					if( cbin < 0 ) cbin = 0;
					if( cbin >= N_CDIST_BINS ) cbin = N_CDIST_BINS-1;
					lipids[proc[p].rec].cdist[cbin] += 1;
				
#ifdef DEBUG
					fprintf(debugFile, "draw line { %lf %lf %lf } { %lf %lf %lf }\n", source_point_upper[0], source_point_upper[1], source_point_upper[2], close_point_upper[0], close_point_upper[1], close_point_upper[2] );
#endif
				//	printf("%s UPPER %le c %le\n", at[atn].resname, close_upper, c_at_u );
				}
				else
				{	
					lipids[proc[p].rec].sum_c_obs_lower += c_at_l;
					lipids[proc[p].rec].n_obs_lower += 1;
					int cbin = N_CDIST_BINS * (c_at_l-CMIN)/(CMAX-CMIN);
					if( cbin < 0 ) cbin = 0;
					if( cbin >= N_CDIST_BINS ) cbin = N_CDIST_BINS-1;
					lipids[proc[p].rec].cdist[cbin] += 1;
#ifdef DEBUG
					fprintf(debugFile, "draw line { %lf %lf %lf } { %lf %lf %lf }\n", source_point_lower[0], source_point_lower[1], source_point_lower[2], close_point_lower[0], close_point_lower[1], close_point_lower[2] );
#endif
				//	printf("%s LOWER %le c %le\n", at[atn].resname, close_lower, c_at_l );
				}
				}
				fflush(stdout);	
			}
			
			printf("LIPIDS UPPER");
			for( int r = 0; r < nrecs; r++ )
			{
				if( lipids[r].n_obs_upper > 0 )
					printf(" %s %le", lipids[r].resname, lipids[r].sum_c_obs_upper / lipids[r].n_obs_upper );
				lipids[r].sum_c_obs_upper = 0;
				lipids[r].n_obs_upper = 0;
			}
			printf(" LOWER");
			for( int r = 0; r < nrecs; r++ )
			{
				if( lipids[r].n_obs_lower > 0 )
					printf(" %s %le", lipids[r].resname, lipids[r].sum_c_obs_lower / lipids[r].n_obs_lower );
				lipids[r].sum_c_obs_lower = 0;
				lipids[r].n_obs_lower = 0;
			}
			printf(" CHI2: %le", chi2 );
			printf("\n");	
			ntraj++;

			for( int a = 0; a < nat; a++ )
				at[a].zap();

			// we can know if we've already passed through the main loop.
			trip = 0;

#ifdef DEBUG
			fclose(debugFile);
			exit(1);
#endif
		}
	}
}

static int opt_n_pts = 0;
static surface *opt_surface1=NULL;
static surface *opt_surface2=NULL;
static double opt_La = 0;
static double opt_Lb = 0;
static double opt_Lc = 0;
static double *r_set = NULL;
static double *opt_ru = NULL;
static double *opt_rl = NULL;

double pt_amoeba( double *shifts )
{
	double shift_dx = shifts[1];
	double shift_dy = shifts[2];
	double shift_dz = shifts[3];

	while( shift_dx < -opt_La/2 ) shift_dx += opt_La;
	while( shift_dy < -opt_Lb/2 ) shift_dy += opt_Lb;
	while( shift_dz < -opt_Lc/2 ) shift_dz += opt_Lc;
	while( shift_dx >  opt_La/2 ) shift_dx -= opt_La;
	while( shift_dy >  opt_Lb/2 ) shift_dy -= opt_Lb;
	while( shift_dz >  opt_Lc/2 ) shift_dz -= opt_Lc;

	double chi2 = 0;

	for( int xx = 0; xx < opt_n_pts; xx++ )
	{
		int pbc_spot = 0;

		double close_lower = 1e10;
		double close_upper = 1e10;
		double close_point_lower[3] = { 0,0,0};	
		double close_point_upper[3] = { 0,0,0};	

		for( int dx = -1; dx <= 1; dx++ )
		for( int dy = -1; dy <= 1; dy++ )
		for( int dz = -1; dz <= 1; dz++ )
		{
			double rp[3] = { r_set[3*xx+0] + dx * opt_La + shift_dx, r_set[3*xx+1] + dy * opt_Lb + shift_dy, r_set[3*xx+2] + dz * opt_Lc + shift_dz };
		
			if( rp[0] < -opt_La * 0.65 ) continue;
			if( rp[0] >  opt_La * 0.65 ) continue;
			if( rp[1] < -opt_Lb * 0.65 ) continue;
			if( rp[1] >  opt_Lb * 0.65 ) continue;
			if( rp[2] < -opt_Lc * 0.65 ) continue;
			if( rp[2] >  opt_Lc * 0.65 ) continue;

			double dist_u;
			int f_u;
			double col_u_u, col_v_u;
			//opt_surface1->nearPointOnSurface( rp, &f_u, &col_u_u, &col_v_u, M, mlow, mhigh, &dist_u );					
			opt_surface1->nearPointOnBoxedSurface( rp, &f_u, &col_u_u, &col_v_u, M, mlow, mhigh, &dist_u );					
			double pt_upper[3], nrm_upper[3];

			opt_surface1->evaluateRNRM( f_u, col_u_u, col_v_u, pt_upper, nrm_upper, opt_ru );
			double dist_l;
			int f_l;
			double col_u_l, col_v_l;
			opt_surface2->nearPointOnBoxedSurface( rp, &f_l, &col_u_l, &col_v_l, M, mlow, mhigh, &dist_l );					

			double pt_lower[3], nrm_lower[3];

			opt_surface2->evaluateRNRM( f_l, col_u_l, col_v_l, pt_lower, nrm_lower, opt_rl );


			double dr_u[3] = { rp[0] - pt_upper[0], rp[1] - pt_upper[1], rp[2] - pt_upper[2] };
			while( dr_u[0] < -opt_La/2 ) dr_u[0] += opt_La;
			while( dr_u[1] < -opt_Lb/2 ) dr_u[1] += opt_Lb;
			while( dr_u[2] < -opt_Lc/2 ) dr_u[2] += opt_Lc;
			while( dr_u[0] > opt_La/2 ) dr_u[0] -= opt_La;
			while( dr_u[1] > opt_Lb/2 ) dr_u[1] -= opt_Lb;
			while( dr_u[2] > opt_Lc/2 ) dr_u[2] -= opt_Lc;
			double dr_l[3] = { rp[0] - pt_lower[0], rp[1] - pt_lower[1], rp[2] - pt_lower[2] };
			while( dr_l[0] < -opt_La/2 ) dr_l[0] += opt_La;
			while( dr_l[1] < -opt_Lb/2 ) dr_l[1] += opt_Lb;
			while( dr_l[2] < -opt_Lc/2 ) dr_l[2] += opt_Lc;
			while( dr_l[0] > opt_La/2 ) dr_l[0] -= opt_La;
			while( dr_l[1] > opt_Lb/2 ) dr_l[1] -= opt_Lb;
			while( dr_l[2] > opt_Lc/2 ) dr_l[2] -= opt_Lc;
			
			double r2_u = sqrt(dr_u[0]*dr_u[0] + dr_u[1]*dr_u[1] +dr_u[2]*dr_u[2]);
			double r2_l = sqrt(dr_l[0]*dr_l[0] + dr_l[1]*dr_l[1] +dr_l[2]*dr_l[2]);
	

			if( r2_u < close_upper )
			{
				close_upper = r2_u;
				close_point_upper[0] = pt_upper[0];
				close_point_upper[1] = pt_upper[1];
				close_point_upper[2] = pt_upper[2];
			}
			if( r2_l < close_lower )
			{
				close_lower = r2_l;
				close_point_lower[0] = pt_lower[0];
				close_point_lower[1] = pt_lower[1];
				close_point_lower[2] = pt_lower[2];
			}
		}
	
		if( close_upper < close_lower )
		{
			chi2 += close_upper*close_upper;
		}
		else
		{	
			chi2 += close_lower*close_lower;
		}

//		printf("close: %le %le chi2: %le shift %le %le %le\n", close_upper, close_lower, chi2, shift_dx, shift_dy, shift_dz );
	}

//	printf("chi2: %le\n", chi2 );
	return chi2;
}

double opt_shift( surface *surface1, surface *surface2, double *pt_set, int npts, double La, double Lb, double Lc, double shift[3], int just_report )
{
	opt_ru = (double *)malloc( sizeof(double) * 3 * (surface1->nv+1) );
	opt_rl = (double *)malloc( sizeof(double) * 3 * (surface2->nv+1) );

	opt_ru[3*surface1->nv+0] = 1.0;
	opt_ru[3*surface1->nv+1] = 1.0;
	opt_ru[3*surface1->nv+2] = 1.0;
	
	opt_rl[3*surface2->nv+0] = 1.0;
	opt_rl[3*surface2->nv+1] = 1.0;
	opt_rl[3*surface2->nv+2] = 1.0;

	surface1->get(opt_ru);
	surface2->get(opt_rl);

	opt_surface1 = surface1;
	opt_surface2 = surface2;
	opt_La = La;
	opt_Lb = Lb;
	opt_Lc = Lc;
	r_set = pt_set;
	opt_n_pts = npts;


	int nparams = 3;
	
	double **p = (double **)malloc( sizeof(double*) * (2+nparams) );
	for( int i = 0; i <= 1+nparams; i++ )
	        p[i] = (double *)malloc( sizeof(double) * (1+nparams) );
	double *y = (double *)malloc( sizeof(double) * (nparams+2) );

	double vals[4] = { 0, 0, 0, 0 };	
	
	for( int i = 1; i <= 1+nparams; i++ )
	for( int j = 1; j <= nparams; j++ )
	{
	        if( i == 1 )
	                p[i][j] = vals[j];
	        else if( i-1 == j )
	                p[i][j] = vals[j] + 1;
	        else
	                p[i][j] = vals[j];
	}
	
	
	for( int i = 1; i <= nparams+1; i++ )
	{
	        y[i] = pt_amoeba(p[i]);
		if( just_report ) break;
	}
	int nfunk = 0;
	if( !just_report )		
		amoeba( p, y, nparams, 1e-3, pt_amoeba, &nfunk );


	shift[0] = p[1][1];
	shift[1] = p[1][2];
	shift[2] = p[1][3];

	while( shift[0] < -opt_La/2 ) shift[0] += opt_La;
	while( shift[1] < -opt_Lb/2 ) shift[1] += opt_Lb;
	while( shift[2] < -opt_Lc/2 ) shift[2] += opt_Lc;
	
	while( shift[0] >  opt_La/2 ) shift[0] -= opt_La;
	while( shift[1] >  opt_Lb/2 ) shift[1] -= opt_Lb;
	while( shift[2] >  opt_Lc/2 ) shift[2] -= opt_Lc;

	double chi2 = y[1];
	free(y);
	for( int i = 0; i <= 1 + nparams; i++ )
		free(p[i]);
	free(p);
	
	free(opt_ru);
	free(opt_rl);	

	return chi2;
}





