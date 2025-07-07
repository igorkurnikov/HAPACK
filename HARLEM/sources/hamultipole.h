/*! \file hamultipole.h
    
    interface to classes defining 1e operators
 
    \author  Igor Kurnikov  
    \date 1998-2004

*/   
#ifndef HAMULTIPOLE_H
#define HAMULTIPOLE_H

#include "halinalg.h"
#include "hastl.h"

class HaQCMod;
class HaBasisSet;
class GauBasisSet;

class HaOperR 
//! Electric Dipole moment operator
{
public:
  
  HaOperR();
  virtual ~HaOperR(void);

  int FillMat(HaBasisSet* pbset,HaMat_doubleArr& rmats);
  int EvalGauBasisSet(GauBasisSet* pbset,HaMat_doubleArr& rmats);

  protected:
	
};


class HaOperRDelt 
//! Operator R X Grad ( -i Magnetic moment operator)
{
public:
  
  HaOperRDelt(HaQCMod* ptr_qc_mod_new);
  virtual ~HaOperRDelt(void);
  
  int FillMat(HaBasisSet* pbset, HaMat_doubleArr& rmats);
  int EvalGauBasisSet(GauBasisSet* pbset,HaMat_doubleArr& rmats);

  bool RecalcFromHr(); //!< calculate Rx Grad matricies using Grad=[H,r]
  bool RecalcFromHr2(); //!< calculate Rx Grad matricies using Grad=[H,r]
  bool RecalcLondon(GauBasisSet* pbset); //!< set London Fock matrix for magnetic operator 
  bool LondonDaltonCalc(); //!< Assign London orbital magnetic operator 
                           //!< using DALTON 

  void SetLondonAO(){ i_lond=1; }

  HaQCMod* ptr_qc_mod;
  HaMat_doubleArr data;
  int i_lond;
	         
};


class HaOperGrad
//! Operator Delta ( Gradient)
{
public:
  HaOperGrad();
  virtual ~HaOperGrad(void);

  int FillMat(HaBasisSet* pbset, HaMat_doubleArr& rmats);
  int EvalGauBasisSet(GauBasisSet* pbset,HaMat_doubleArr& fmats);
	         
};

class HaOperKinEner 
//! Kinetic Energy Operator
{
public:
	HaOperKinEner();
	virtual ~HaOperKinEner();

//	int FillMat(HaBasisSet* pbset, HaMat_doubleArr& rmats);
	int EvalGauBasisSet(GauBasisSet* pbset, HaMat_double& fmat);
	
};


#endif /* !HAMULTIPOLE_H */
