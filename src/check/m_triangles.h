#ifndef __m_trianglesh__
#define __m_trianglesh__

#define MAX_NEIGHBORS 20
struct m_triangle;

struct triHeader
{	
	double use_La;
	int nsub;	
	int level;
	int hex_bc;
};

struct triElement
{
	double vertices[9];
	double cen[3];
	double avh;
	int navh;
	double vector[3];
	int nvec;
	int vertices_setup;
};

typedef struct m_triangle
{
	int unique_index[3];
	double cen[3];
	double vertices[9];
	m_triangle *sub_m_triangles;
	int nsub;
	int vertices_setup;
	double avh;
	double navh;

	double vector[3];
	int nvec;

	void init( void );
	void subdivide( int level );
	struct m_triangle *nearm_triangle( double *r ); 
	void addValue( double h, double *r, double *PBC_vecs );
	void addVector( double *vec, double *r, double *PBC_vecs );
	void fillm_triangle( FILE *theFile, double color_min, double color_max, double *color1, double *color2, double scale, double *shift, int by_density=0 );
	void fillDensity( FILE *theFile, double color_min, double color_max, double *color1, double *color2, double scale, double *shift );
	void drawBorder( FILE *theFile, double scale, double shift[3]);
	void drawArrow( FILE *theFile, double scale, double *shift );
	void printValue( void );
	int countm_triangles(void);
	void writem_triangles( FILE *theFile );
	void readm_triangles( FILE *theFile );
	void writeToFile( FILE *theFile ,double use_La, int level, int do_hex_bc );
	void readFromFile( FILE *theFile ,double *use_La, int *do_hex_bc );
	int nmax( void );
	void clearUnique( void );
	void assignUnique( m_triangle *header, int *cur_index, double * PBC_vecs );
	int assignUniqueIndices( double *PBC_vecs );

	int matchUnique(  double *r, double *PBC_vecs );
	int printPosition( int index );
	void testThreshold( int *assigned, double cut );
	void getNeighbors( int *nneighbors, int *neighbors );
	void writeLatticePoints( double *r, double *w, double *PBC_vec );
	int countAssignedm_triangles( int *assigned );
	void printAssignedm_triangles( FILE * theFile, int *assigned, int *map );
} m_triangle;

#endif
