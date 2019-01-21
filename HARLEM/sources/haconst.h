/*! \file haconst.h
 
   Application-wide Constants and Macros in HARLEM
 
   Partially Derived from:

   rasmol.h
   RasMol2 Molecular Graphics
   Roger Sayle, August 1995
   Version 2.6

  \author Igor Kurnikov
  \date 1999-2002
 
*/

#if !defined(HACONST_H)
#define HACONST_H

#if !defined(SWIG)
#ifndef True
const int True = 1;
const int False = 0;
#endif
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE  0
#endif

#ifndef PI   /* Avoid Linux Warnings! */
const double PI = 3.14159265358979323846;
#endif
#if !defined(CONSTANTS_H)
const double BOHR_TO_ANG       = 0.529177249;        //!< Bohrs to Angstroms
#endif

const double ANG_TO_BOHR       = 1.0/BOHR_TO_ANG;    //!< Angstroms to Bohrs
const double HARTREE_TO_EV     = 27.21138386;        //!< HARTREE TO eV  +/- 0.00000068 eV

#if !defined(CONSTANTS_H)
const double HARTREE_TO_KT     = 1059.642;       //!< HARTREE TO kT at 298K
#endif
const double KT_TO_KJ_PER_MOL     = 2.47638;       //!< kT at 298K to kJ/mol
const double HARTREE_TO_JOULES = 4.35974394e-18;  //!< Joules(J) in Hartree - exact
const double AU_TO_KG          = 1.6605402e-27;  //!< Kg per atomic unit mass  

const double EV_TO_KT       = HARTREE_TO_KT/HARTREE_TO_EV;    //!< EV to kT at 298K
const double AVOG_NUM       = 6.02214179e+23;  //!< mol^-1 Avogadro constant +/- 0.00000030 e+23 
const double CAL_TO_JOULES  = 4.184;          //!< Joules per calorie (thermochemical) exact
const double PLANCK         = 6.6260755e-34;  //!< Planck constant (Joules * sec)  
const double LIGHT_SPEED_CM = 2.99792458e+10; //!< Speed of light (cm/s)
const double BOLTZ          = 1.3806504e-23;   //!< Boltzmann constant (Joules per Kelvin) +/- 0.0000024e-23
const double R_HALF_KCAL    = BOLTZ*AVOG_NUM*0.5*1.0e-3/CAL_TO_JOULES;  //!< 1/2 R  - multiplied by T gives kinetic energy per degree of freedom

const double HARTREE_TO_KCAL=  AVOG_NUM*HARTREE_TO_JOULES/(1000.0*CAL_TO_JOULES);
const double KT_300_KCAL = 300.0 * BOLTZ * AVOG_NUM/(1000.0*CAL_TO_JOULES); //!< kT in kcal/mol for 300K 

const double RAD_TO_DEG = 180.0/PI;
const double DEG_TO_RAD = PI/180.0;

const double EL_ANG_TO_DEBYE = 4.8033324;  //!< conversion (electron * ang to debye )
const double AU_TO_DEBYE     = BOHR_TO_ANG * EL_ANG_TO_DEBYE;   //!< conversion coef ( Atomic units ( electron * bohr) to debye )

typedef int int_4;
typedef unsigned int uint_4;
typedef char int_1;
typedef unsigned char uint_1;
typedef short int_2;
  
// TNT MACROS:

#define TNT_BASE_OFFSET (1)

//---------------------------------------------------------------------
// Define this macro if you want  TNT to ensure all refernces
// are within the bounds of the array.  This encurs a run-time
// overhead, of course, but is recommended while developing
// code.  It can be turned off for production runs.
// 
//       #define TNT_BOUNDS_CHECK
//---------------------------------------------------------------------
//
#define TNT_BOUNDS_CHECK
#ifdef TNT_NO_BOUNDS_CHECK
#undef TNT_BOUNDS_CHECK
#endif

//---------------------------------------------------------------------
// Define this macro if you want to utilize matrix and vector
// regions.  This is typically on, but you can save some
// compilation time by turning it off.  If you do this and
// attempt to use regions you will get an error message.
//
//       #define TNT_USE_REGIONS
//---------------------------------------------------------------------
//
#define TNT_USE_REGIONS

#define AbsFun(a)    (((a)<0)? -(a) : (a))
#define MinFun(a,b)  (((a)<(b))? (a) : (b) )
#define MaxFun(a,b)  (((a)>(b))? (a) : (b) )


#if defined(_MSC_VER) || defined(__DECCXX) || defined(GNU)
#define _fmalloc   malloc
#define _ffree     free
#define _fstrcmp   strcmp
#define _fmemset   memset
#endif

#if defined(_MSC_VER)
#define _fstrnicmp strnicmp
#endif

enum RunMode{ RUN_FOREGROUND=0, RUN_BACKGROUND=1 };  //!< Type of the submit of an external program

const int ItemCount = 8;
const int AdvPickAtom = 0;
const int AdvPickNumber = 1;
const int AdvSelectCount = 2;
const int AdvName        = 3;
const int AdvIdent       = 4;
const int AdvClass       = 5;
const int AdvImage       = 6;
const int AdvPickCoord   = 7;

#undef HARLEM_PYTHON_NO
#undef HARLEM_WX_PLPLOT_NO

#endif /* !defined(HACONST_H) */
