/***************************************************************************
 *   Nikolay Simakov                                                       *
 *   nsimakov@andrew.cmu.edu                                               *
 *   Maria Kurnikova Research Group                                        *
 *   http://crete.chem.cmu.edu/                                            *
 *   Carnegie Mellon University, 2006                                      *
 ***************************************************************************/

#ifndef APBSMOD_H
#define APBSMOD_H

#include "hastl.h"
#include "vec3d.h"
#include "hastring.h"
#include "hacompmod.h"
#include "hasurface.h"

class HaSurface;
class MolSet;
class AtomGroup;
class HaResidue;
class AltChemState;


class APBSMod : public HaCompMod
{
public:
  APBSMod(MolSet* new_phost_mset = NULL);
  virtual ~APBSMod();

  double E;
  int nx, ny, nz;            //!< dimensions of the grid, default = 65 n=c*2^(l+1)+1, C!=0, l=nlev, nlev=4=> n=c*32+1=33,65,97,129,161
  float GridScale;//!<Grid Scale
  float epsi, epsout;       //!< internal and external dielectric constant
  float rionst;             //!< ionic strength , default = 0.0
  float exrad, radprb;      //!< Ion exclusion radius, and probe radius, default exrad = 2.0, radprb = 1.4 
  int boundary;              //!< Boundary: 0=Zero,1=Apprx Coul,2=Focus,3=Coul , default = 0
  int iper[3];               //!< periodic boundary flags
  int nlit;                  //!< number of linear Poisson-Boltz. iterations, default = -1
  int nnit;                  //!< number of non-linear Poisson-Boltz. iterations, default = 0
  char chgm[5];//!< way to distrebute charge: spl0 - 8-nodes, spl2 - Cubic B-spline discretization
  char bcfl[6];//!< zero/sdh/mdh/focus
  char PBEType[5];//!< lpbe or npbe
  static double d;
  int nlev;
  int Run();
protected:
  int SavePQRFile(const char* filename);

};

#endif
