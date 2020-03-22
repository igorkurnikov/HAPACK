/* $ID: g94_globals.h
 *
 * Global variables of Gaussian-94
 * for HARLEM and other program using Gaussian
 */

#ifndef G94_GLOBALS_H
#define G94_GLOBALS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "f2c.h"

//#define GAUSSVER 98
#undef GAUSSVER 

/* GAUSSIAN Chk-point files numbers:  Description is in L1.F */

#ifdef  DEFMAXATM
const int MAXATM  = DEFMAXATM;
#else
const int MAXATM = 20000;
#endif

const int MAXAT1 = MAXATM+1;

#ifdef  DEFMAXSHL
  const int MAXSHL  = DEFMAXSHL;
#else
  const int MAXSHL = 20000;
#endif

const int MAXPRM = 3*MAXSHL;
const int MAXSH1 = MAXSHL+1;

#ifdef  DEFMAXNZ
const int MAXNZ = DEFMAXNZ;
#else
  /* default for G98  */
  const int MAXNZ = 20000;
#endif

const int  MAXSUB = 40;
const int  MAXDIR = 80*sizeof(int);

/* GAUSSIAN COMMON DEFINITIONS */

typedef struct {             // description in L1.F
  int    natoms;         //!<  THE NUMBER OF ATOMS IN THE MOLECULE.
  int    icharge;        //!< THE TOTAL ELECTRIC CHARGE ON THE MOLECULE
  int    multip;         //!< THE MOLECULE'S SPIN MULTIPLICITY
  int    nae;            //!< THE NUMBER OF ELECTRONS OF ALPHA SPIN
  int    nbe;            //!< THE NUMBER OF ELECTRONS OF BETA  SPIN
  int    ne;             //!< THE NUMBER OF ELECTRONS
  int    nbasis;         //!< THE NUMBER OF BASIS FUNCTIONS
  int    ian[MAXATM];    //!< THE ATOMIC NUMBER OF ATOM I.
  double atmchg[MAXATM]; //!< THE NUCLEAR CHARGE OF ATOM I.
  double c[MAXATM][3];   //!< THE X, Y, AND Z, COORDINATES OF THE ATOMS IN
                             //!< ATOMIC UNITS.  THE COORDINATES ARE STORED IN
                             //!< THE ORDER X1, Y1, Z1, X2, Y2, Z2, X3, ETC.
  int    iattyp[MAXATM];  
  double atchmm[MAXATM];
} mol_type;

struct b_type{
  double exx[MAXPRM];
  double c1[MAXPRM];
  double c2[MAXPRM];
  double c3[MAXPRM];
  double x[MAXSHL];
  double y[MAXSHL];
  double z[MAXSHL];
  int    jan[MAXSHL];
  int    shella[MAXSHL];
  int    shelln[MAXSHL];
  int    shellt[MAXSHL];
  int    shellc[MAXSHL];
  int    aos[MAXSHL];
  int    aon[MAXSHL];
  int    nshell;
  int    maxtyp;
};

typedef struct
{
	int maxm[2];
} maxmem_type;


typedef struct {
  double scale[MAXSHL];
} scale_type;

//      COMMON /B/

//      EXX    ... CONTAINS THE GUASSIAN EXPONENTS FOR ALL THE
//                 THE PRIMITIVE SHELLS.  THE ARRAY SHELLA
//                 CONTAINS POINTERS INTO EXX FOR THE VARIOUS
//                 PRIMITIVE SHELLS.
//      C1     ... CONTAINS THE S COEFFICIENTS FOR ALL THE PRIMITIVE
//                 SHELLS.  INDEXED BY SHELLA.
//      C2     ... CONTAINS THE P COEFFICIENTS FOR THE PRIMITIVE
//                 SHELLS.  INDEXED BY SHELLA.
//      C3     ... THIS ARRAY REALLY CONSISTS OF 3 SUB-ARRAYS
//                 (EACH 400 WORDS LONG) THAT CONTAIN THE
//                 D AND F COEFFICIENTS PLUS THE D AND F
//                 POINTER TABLE FOR D AND F TYPE SHELLS.
//                 THESE ARRAYS ARE C3, C4 AND SHLADF, RESPECTIVELY.
//      X      ... THE X-CARTESIAN COORDINATE FOR EACH OF THE
//                 PRIMITIVE SHELLS.
//      Y      ... THE Y-CARTESIAN COORDINATE FOR EACH OF THE
//                 PRIMITIVE SHELLS.
//      Z      ... THE Z-CARTESIAN COORDINATE FOR EACH OF THE
//                 PRIMITIVE SHELLS.
//      JAN    ... Number of the center associated with each shell.
//      SHELLA ... SHELLA(I) CONTAINS THE STARTING LOCATION
//                 WITHIN (EXX,C1,C2) OF THE DATA
//                 (EXPONENTS, S-COEFFICIENTS, P-COEFFICIENTS)
//                 FOR THE I-TH PRIMITIVE SHELL.
//      SHELLN ... SHELLN(I) CONTAINS THE NUMBER OF PRIMITIVE
//                 GAUSSIANS IN THE I-TH PRIMITIVE SHELL.
//      SHELLT ... CONTAINS THE MAXIMUM ANGULAR QUANTUM NUMBER
//                 OF THE I-TH SHELL.
//      SHELLC ... CONTAINS THE SHELL CONSTRAINT FOR THE I-TH
//                 SHELL (SEE TABLE BELOW).
//      AOS    ... AOS(I) GIVES THE STARTING ATOMIC ORBITAL BASIS
//                 FUNCTION NUMBER (IE NUMBER WITHIN THE LIST
//                 OF ATOMIC ORBITAL BASIS FUNCTIONS) OF THE
//                 I-TH SHELL.  NOTE THAT AOS IS ALWAYS FILLED
//                 AS THOUGH THE SHELL CONTAINED ALL POSSIBLE
//                 LOWER ANGULAR MOMENTUM FUNCTIONS.  SEE ROUTINE
//                 GenBas FOR DETAILS ON HOW AOS IS FILLED.
//      AON    ... Name of the shell (e.g. 2SI, 5D, 10F).
//      NSHELL ... THE NUMBER OF PRIMITIVE SHELLS.
//      MAXTYP ... THE HIGHEST ANGULAR QUANTUM NUMBER PRESENT
//                 IN THE BASIS.

//      THE FOLLOWING TABLE SUMMARIZES THE RELATIONSHIP BETWEEN
//      SHELLT AND SHELLC.

//      =========================================
//      TYPE   FUNCTIONS          SHELLT   SHELLC
//      =========================================
//        S     S                     0        0
//       SP     S PX PY PZ            1        0
//      SPD     S PX PY PZ            2        0
//              XX YY ZZ XY XZ YZ
//        P     PX PY PZ              1        1
//        D     XX YY ZZ XY XZ YZ     2        2
//        F     XXX YYY ZZZ XYY       3        2
//              XXY XXZ XZZ YZZ
//              YYZ XYZ
//      =========================================



typedef struct {
  int nchain;
  int idum;
} nchain_type;

typedef struct {
  int iunit[20]; // see description in L1.F
} munit_type;


//       COMMON/MUNIT/ IUNIT(20)

//       DEFINES THE UNIT NUMBERS FOR ALL EXTERNAL DATA FILES NEEDED
//       BY GAUSSIAN.  IUNIT IS INITIALIZED BY SUBROUTINE DEFUNT IN
//       LINK 1.
//       THE PURPOSE OF THIS COMMON IS TO CENTRALIZE THE DEFINITIONS
//       OF THE FORTRAN LOGICAL UNITS REQUIRED IN THE GAUSSIAN SYSTEM.
//       THUS, IF IT IS EVER NECESSARY TO CHANGE TO DIFFERENT LOGICAL
//       UNITS ON ANOTHER MACHINE, ONE SHOULD BE ABLE TO MERELY
//       CHANGE THE DEFINITIONS IN DEFUNT.


//       =================================================================
//       IUNIT( )            USE
//       =================================================================
//          1           NOT USED.
//          2           PRIMARY FORTRAN INPUT (internal input file).
//          3           PRIMARY FORTRAN PRINTED OUTPUT.
//          4           PRIMARY FORTRAN PUNCHED OUTPUT.
//          5           NOT USED.
//          6           NOT USED.
//          7           NTRAN unit for the read-write file.
//                      Effectively hard-wired to 1 everywhere.
//          8           Fortran I/O unit for DRT information.  See L806
//                      for details.  Not used in standard SDGuga routes.
//          9           CHECK-POINT FILE.
//          12          INTEGRALS.
//          13          Sorted AO integrals.
//          16          D2E UNIT.
//          17          Output unit for the archive entry.  Not used under
//                      VAX/VMS or Cray/COS.
//          18          ECP formula tape unit.
//          19          Polyatom integral tape unit.
//      ==================================================================


typedef struct {
  double  phycon[30];
} phycon_type;

//       COMMON /PHYCON/ PHYCON(30)
//
//       THIS COMMON BLOCK HOLDS THE VALUES OF THE PHYSICAL CONSTANTS
//       WHICH ARE USED IN THE PROGRAM.  IT IS INITIALIZED BY SUBROUTINE
//       PHYFIL IN LINK1 WHICH ALSO DOCUMENTS THE SOURCES FOR THE NUMBERS.
//
//       PHYCON( 1) ... TOANG,  ANGSTROMS PER BOHR
//       PHYCON( 2) ... TOKG,   KILOGRAMS PER AMU
//       PHYCON( 3) ... TOE,    ESU PER ELECTRON CHARGE
//       PHYCON( 4) ... PLANK,  PLANK CONSTANT
//       PHYCON( 5) ... AVOG,   AVOGADRO NUMBER
//       PHYCON( 6) ... JPCAL,  JOULES PER CALORIE
//       PHYCON( 7) ... TOMET,  METRES PER BOHR
//       PHYCON( 8) ... HARTRE, JOULES PER HARTRE
//       PHYCON( 9) ... SLIGHT, SPEED OF LIGHT
//       PHYCON(10) ... BOLTZ,  BOLTZMAN CONSTANT
//       PHYCON(11) ... Fine structure constant.
//       PHYCON(12) ... EMKG, Electron mass in KG.
//       PHYCON(13) ... NOT USED AT PRESENT


typedef struct {
  int  iop[50];
} iop_type;

typedef struct {
  int  nsubst;
  int  ipdsub;
  int  linkn  [MAXSUB];
  int  lendir [MAXSUB];
  int  linkd  [MAXSUB][MAXDIR];
} substn_type;

typedef struct {
  double tstart[3];
  double tstop[3];
  double elapsd[3];
 } clcks_type;

typedef struct {
  int killnk;
  int kilcnt;
  int refcnt;
  int kill;
  int numprc;
  int klpd1[2];
  double jobtim;
 } kjob_type;

typedef struct {
   logical debug;
   int ntrout;
} ntr002_type;

typedef struct {
  int ispect;
  int lspect;
  int nrorb;
  int noa;
  int nva;
  int nob;
  int nvb;
  int noaob;
  int noava;
  int noavb;
  int nobva;
  int nobvb;
  int nvavb;
  int noa2;
  int noa3;
  int nob2;
  int nob3;
  int nva2;
  int nva3;
  int nvb2;
  int nvb3;
  int novaa;
  int novab;
  int novbb;
  int maxbuc;
  int ieval;
  int ioab;
  int loab;
} orb_type;

typedef struct {
  int idb1;
  int idb2;
  int idb3;
  int idb4;
  int idb5;
  int idb6;
  int idb7;
  int idb8;
  int idb9;
  int idb10;
  int idb11;
  int idb12;
  int idb13;
  int idb14;
  int idb15;
  int idb16;
  int idb17;
  int idb18;
  int idb19;
  int idb20;
  int idb21;
  int iad1;
  int iad2;
  int iad3;
  int ias1;
  int ias2;
  int iwd1;
  int iwd2;
  int iwd3;
  int iws1;
  int iws2;
  int iscr1;
  int iscr2;
  int iscr3;
  int iscr4;
  int iscr5;
  int iscr6;
  int iscr7;
  int iscr8;
  int iscr9;
  int iscr10;
  int iscr11;
  int iscr12;
  int iscr13;
  int iscr14;
  int iscr15;
  int iscr16;
  int iscr17;
  int iscr18;
  int iscr19;
} bucknr_type;

typedef struct {
  int in;   
  int iout;
  int ipunch;
} io_type;

typedef struct {
  int nsymop;
  int nreps;
  int lblrep[32];
  double chrtbl[16][10];
  double symops[10][9];
  int iprmut[10][MAXATM];
} repcom_type;

typedef struct {
  int ismode;
  int mode;
  int istat;
  int last;
  int ntx;
  int iux[5];
  int icon;
  int nrpext;
  int knit1;
  int knit2;
  int ibase;
  int ibasd[2];
  int dbase;
  int dbasd[2];
  int ireset[2];
  int iq;
  int ifil;
  int intcnt;
  int itotal;
  int limint;
  int nwpi;
  int nwiib;
  int isym2e;
} ibf_type;

#ifdef GAU_MAIN

int      f77argc  = 0;     /*   command line arguments          */
char*    f77argv  =NULL;   /*   for fortran getarg and iargc    */

#else

extern int     f77argc ;     /*   command line arguments          */
extern char*   f77argv ;     /*   for fortran getarg and iargc    */

#endif

#if defined(_MSC_VER)
#define DllImport __declspec( dllimport ) 
#else
#define DllImport extern
#endif

#if GAUSSVER == 98

/* Gaussian COMMONS */

DllImport mol_type     mol_;
DllImport int      info_[20];  // see description in L1.F
DllImport nchain_type  nchain_;
DllImport munit_type   munit_;
DllImport phycon_type  phycon_;
DllImport iop_type     iop_;
DllImport substn_type  substn_;
DllImport clcks_type   clcks_;
DllImport kjob_type    kjob_;
DllImport ntr002_type  ntr002_;
DllImport b_type       b_;
DllImport b_type       b2_;
DllImport orb_type     orb_;
DllImport bucknr_type  bucknr_;
DllImport  io_type      io_;
DllImport repcom_type  repcom_;
DllImport ibf_type     ibf_;
DllImport scale_type   scale_;
DllImport maxmem_type  maxmem_;

#endif

#ifdef __cplusplus
        }
#endif


#endif  /* End G94_GLOBALS_H */
