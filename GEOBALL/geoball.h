/*****************************************************************************
 * GeoBall.c  : Computes the volume and surface area of a union of balls
 *             using the Alpha Shape theory introduced by Herbert
 *             Edelsbrunner.
 *             Also computes the derivatives of the volume and surface
 *             areas with respect to the coordinates of the centers of
 *	       the balls, if needed.
 *
 * Information on the Alpha Shape method, and on the calculationof the derivatives 
 * can be found in:
 *             
 *		1. H. Edelsbrunner, "The Union of balls and its dual shape", Discrete 
 *		   Comput. Geom., 13:415-440 (1995)
 *
 *		2. H. Edelsbrunner and E.P. Mucke, "Three-dimensional alpha shapes",
 *                 ACM Trans. Graphics, 13:43-72 (1994)
 *
 *              3. J. Liang, H. Edelsbrunner, P. Fu, P.V. Sudhakar and S. Subramaniam,
 *                 "Analytical shape computation of macromolecules: I. Molecular area",
 *                 Proteins: Struct. Func. Genet., 33:1-17 (1998)
 *
 *		4. H. Edelsbrunner and P. Koehl, "The weighted volume derivative of a
 *		   space filling diagram", Proc. Natl. Acad. Sci. (USA), 100:2203-2208,
 *		   (2003)
 *
 *		5. R. Bryant, H. Edelsbrunner, P. Koehl and M. Levitt, "The area
 *		   derivative of a space-filling diagram", Discrete Comput. Geom.,
 *                 32:293-308 (2004)
 *
 *		6. H. Edelsbrunner and P. Koehl, "The geometry of biomolecular solvation",
 *		   Combinatorial and Discrete Geometry (MSRI Series), to appear (2005)
 *
 * Calculations in the program GeoBall are performed first in floating point; if the
 * results are inconclusive, the programs switches to arbitrary precision arithmetics,
 * using the package GMP. 
 *
 * In cases of degeneracies (i.e. a geometric predicate has a value of 0), the
 * program applies the procedure defined as "Simulation of Simplicity"
 *
 *		1. H. Edelsbrunner and E.P. Mucke, "Simulation of Simplicity: a technique
 *		   to cope with degenerate cases in geometric algorithms", ACM trans.
 *		   Graphics, 9:66-104 (1990)
 *
 * For historical reasons, the program is written as a mixture of Fortran and C.
 *
 * The whole suite of program is distributed under the GNU Lesser general Public License.
 *
 * Copyright (C) 2007 Patrice Koehl
 *
 */

#include "defines.h"
#include "stdio.h"
#include "time.h"


/*****************************************************************************************
 * Naming convention between C and Fortran 
 *
 * Let us consider two Fortran subroutines : foo, and foo_with_underscore
 *
 * Case 1: the name of the Fortran subroutine (foo) does not contain an underscore (_)
 *	   In the C program, we must add ONE underscore to the Fortran name:
 *	   the C program refers to the subroutine as:
 *		foo_(...)
 *	   while the Fortran code writes it as:
 *		foo(...)
 *	   This is independent of the compiler pair (at least for gcc/f77 and 
 *	   Intel icc/ifort)
 *
 * Case 2: the name of the Fortran subroutine (foo_with_underscore) contains at least one 
 * underscore (_)
 *	   Treatment of case 2 is compiler dependent!
 *	   - The Intel compiler treats this case as if it was case 1, i.e. requires 
 *	     ONE underscore at the end of the Fortran name in the C program
 *	   - the gnu f77 however requires TWO underscores at the end of the Fortran name 
 *	     in the C program.
 *
 * I solve this by introducing two functions, FTName1 and FTName2, where FTName2 is 
 * compiler dependent
 */

#define F77Name1(x) x##_   /* Need to add one underscore to Fortran program */

#if defined intel || defined(_MSC_VER)
#define F77Name2(x) x##_   /* Need to add one underscore to Fortran program for 
				Intel compiler */
#else
#define F77Name2(x) x##__   /* Need to add two underscores to Fortran program for GNU 
				f77 compiler */
#endif

/**************************************************************************************
* local common variables
*/

static char input_file[FSIZE], output_file[FSIZE], alpha_file[FSIZE];
static char delaunay_file[FSIZE];
#ifndef WIN32
static int const YES = 1;
static int const NO = 0;
static int const ERROR = 0;
static int const SAFE = 1;
static int const ONE = 1;
#else
#define YES 1
#define NO 0
#define ERROR 0
#define SAFE 1
#define ONE 1
#endif

static int AlphaFileName = 0;
static int DelaunayFileName = 0;

static int natom;
static int nat_orig;
static int nredundant;
static int nredalpha;
static double vol,surf,wvol,wsurf;

static int list_redundant[MAX_ATOM];
static int list_redalpha[MAX_ATOM];

/* Variables defining the vertices:
	xcoord, ycoord, zcoord : 3D coordinates of each vertex
	radius		       : radius for weighted Delaunay
	coord		       : array transferred to Fortran
								*/

static double xcoord[MAX_ATOM],ycoord[MAX_ATOM],zcoord[MAX_ATOM],radius[MAX_ATOM];
static double coord[3*MAX_ATOM];

/* Variables for volume and surface calculation:
	coef:		: we compute a weighted sum of the individual surface 
			  and volume; coef is the array containing the weights
	volume, surface	: contribution of each atom to the total surface and volume of the molecule
	deriv_vol
	deriv_surf	: derivatives of surface and volume
								*/
static double coef[MAX_ATOM];
static double volume[MAX_ATOM], surface[MAX_ATOM];
static double deriv_vol[3*MAX_ATOM], deriv_surf[3*MAX_ATOM];


/**************************************************************************************/
/* program parameters */

static char usage[] = "\n\
usage:    \n\
	%s DATA  -s <switch> -o <res_file> \n\
with:     \n\
	DATA ............... path name of input data file \n\
	-s <switch> ........ switch = 0: computes surface area and volume only, \n\
				      1: computes area, volume and derivatives  \n\
	-o <res_file> ...... path name of output result file (default DATA.vol) \n\
        -d <delaunay_file>.. path name of file for simplices of Delaunay (optional) \n\
        -a <alpha_file>..... path name of file for simplices of dual complex (optional) \n\
";

/**************************************************************************************/
/* Common structures between Fortran and C                                             */

/* Common blocks in Fortran are equivalent to data structures in C:

if in Fortran we have:

	real*8	a,b,c
	integer i,j,k

	common /values/ a,b,c,i,j,k

We use in C:

	extern struct {
		double a,b,c;
		int i,j,k;
	} values_

Some remarks:

	- the name of the common black (here "values") should be the same as the name of
	  the struct in C, with the addition of one (or two!) underscore

	- the underscore in the name of the structure follows the same rules as for
	  procedure (see above the definition of F77Name1 and F77Name2)

	- use "extern" if the common block / structure is first referenced in the Fortran
	  program; if the first reference is in the C program, remove extern

	- put variables in the same order in the structure and in the common block!!

For the Delaunay triangulation and the alpha complex, there are a series of variables that can be used indifferently in the C code and the Fortran code:

1) Variables that describe each tetrahedron of the (weighted) Delaunay triangulation:

	ntetra		: number of tetrehedrons
	tetra		: indices of the 4 vertices defining the tetrahedron
	tetra_neighbour	: the four neighbours of the tetrahedron. If the tetrahedron is (a,b,c,d),
			  the four neighbors are given in that order:
				tetraneighbor[i][0] is the neighbor sharing face b,c,d
				tetraneighbor[i][1] is the neighbor sharing face a,c,d
				tetraneighbor[i][2] is the neighbor sharing face a,b,d
				tetraneighbor[i][3] is the neighbor sharing face a,b,c
	tetra_nindex	: for each neighbor of the tetrahedron:
			  the neighbor of a,b,c,d sharing face b,c,d is a tetrahedron A,B,C,D
			  that in fact corresponds to a permutation of (b,c,d,e) where e is the
			  fourth vertex. tetra_nindex[][0] gives the position of e in (A,B,C,D)
			  (1, 2 3 or 4)
	tetra_link	: indices (in triangle list), of the 4 triangles (b,c,d), (a,c,d),
			  (a,b,d) and (a,b,c)
	tetra_status	: flag for each tetrahedron:
				1 if the tetrahedron is active
				0 if the tetrahedron is "dead" (i.e. does not belong to Delaunay
								triangulation)
	tetra_orient	: flag for orientation of the tetrahedron:
				1 if (a,b,c,d) in ccw order
				-1 otherwise
	tetra_hull	: flag for position of the tetrahedron:
				1 if (at least) one face is on the convex hull
				0 if tetrehedron is interior
	order		: Order of the tetrahedrons, to fit with a first depth search

2) Variables that describe each triangle of the (weighted) Delaunay triangulation:

	ntrig		: number of triangles
	trig		: indices of the 3 vertices defining the triangle
	trig_alp	: flag for each triangle:
				1 if the triangle is active
				0 if the triangle is "dead" (i.e. does not belong to Alpha complex)
	trig_link	: indices (in edge list), of the 3 edges (b,c), (a,c),
			  and (b,d)
	trig_hull	: flag for position of the triangle:
				1 if triangle is on the convex hull
				0 if triangle is interior
	trig_coef	: flag for each triangle (used for volume calculation only)
	
3) Variables that describe the edges of the (weighted) Delaunay triangulation:

	nedge		: number of edges
	edge		: indices of the 2 vertices defining the edge
	edge_hull	: flag for position of the edge:
				1 if both vertices are on convex hull
				0 otherwise
	edge_alp	: flag for each edge:
				1 if the edge is active
				0 if the edge is "dead" (i.e. does not belong to Alpha complex)

4) Variables that link vertices to tetrahedra

	vertex_nlink	: number of tetrahedra attached to each vertex
	vertex_link	: list of tetrahedra attached to each vertex
												*/
/* Now list all Fortran common blocks / C data structures that lies at the interface */


/**************************************************************************************/
/* local procedures with common variables */

int command_line(char *argv[], int argc);
int read_data();
int alpha_setup();
int alpha_geom(int switchval);
void get_volume();
void get_deriv();
void print_volume(double vol, double surf);
void write_delaunay();
void write_alphacx();

extern void clear_all_gmp_arrays(int *nvertex);
extern void clear_alf_gmp();
extern void F77Name2(adjust_nsphere)(double coord[3* MAX_ATOM], double radius[MAX_ATOM], int *natom);
extern void F77Name1(setup)(double coord[3* MAX_ATOM], double radius[MAX_ATOM], int *natom);
extern void F77Name1(regular3)(int *nredundant, int list_redundant[MAX_ATOM]);
extern void F77Name2(define_triangles)();
extern void F77Name2(define_edges)();
extern void F77Name2(remove_inf)();
extern void F77Name1(alfcx)(double *alpha_val, int *nred, int list_redalpha[MAX_ATOM]);
extern void F77Name2(readjust_nsphere)(int *natom, int *nred, int list_redalpha[MAX_ATOM]);
extern void F77Name2(compute_vol)(double coef[MAX_ATOM], int *option);
extern void F77Name2(alpha_result)(double *vol, double *surf, double *wvol, double *wsurf, double volume[MAX_ATOM], double surface[MAX_ATOM]);
extern void F77Name2(alpha_deriv)(double deriv_vol[3*MAX_ATOM],double deriv_surf[3*MAX_ATOM]);


extern struct{
	int ntetra;
	int tetra[MAX_TETRA][4];
	int tetra_neighbor[MAX_TETRA][4];
	} F77Name2(tetra_zone);

extern struct{
	char tetra_status[MAX_TETRA];
	char tetra_orient[MAX_TETRA];
	char tetra_nindex[MAX_TETRA][4];
	} F77Name2(tetra_stat);

extern struct{
	int tetra_position[MAX_TETRA];
	} F77Name2(tetra_print);

extern struct{
	int ntrig;
	int trig[MAX_TRIG][3];
	} F77Name2(trig_zone);

extern struct{
	char trig_coef[MAX_TETRA];
	} F77Name2(trig_fact);

extern struct{
	int nedge;
	int edge[MAX_EDGE][2];
	} F77Name2(edge_zone);

extern struct{
	int npoints;
	int nvertex;
	int redinfo[MAX_POINT];
	} F77Name2(vertex_zone);

extern struct{
	int tetra_link[MAX_TETRA][4];
	int trig_link[MAX_TRIG][3];
	} F77Name1(links);

extern struct {
	char tetra_hull[MAX_TETRA];
	char trig_hull[MAX_TRIG];
	char edge_hull[MAX_EDGE];
	char vertex_hull[MAX_POINT];
	} F77Name1(hull);

extern struct {
	int vertex_nlink;
	int vertex_link[MAX_LINK];
	} F77Name2(vertex_tet);

extern struct {
	char tetra_alp[MAX_TETRA];
	char trig_alp[MAX_TRIG];
	char edge_alp[MAX_EDGE];
	char vertex_alp[MAX_POINT];
	} F77Name2(alpha_complex);

int GeoBallMain(int argc, char **argv);