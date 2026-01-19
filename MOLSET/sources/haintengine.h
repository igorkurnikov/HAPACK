/*! \file haintengine.h
 
    Classes to calculate elemenatal integrals on basis function in HARLEM
    and make summations over densities 
 
    \author Igor Kurnikov 
 
    \date 1998-2004
*/

#ifndef HAINTENGINE_H
#define HAINTENGINE_H

#include "halinalg.h"
#include "haatombasis.h"

class HaQCMod;
class HaBasisSet;

enum TwoElComb{REG=0,RAFF1=11,RAFF2=12,RAFF3=13}; 
// type of 2-e combinations 
// REG - regular 2-e integrals (ia|jb) in MO or (ks|lt) AO notation
// RAFF1 - Raffenetti 1 combination: 
//      in MO notation: (A(1)-B)_ia,jb  4(ia|jb) - [(ij|ab) + (ib|ja)]     
//      in AO notation: (A(1)-B)'_ks,lt  4(ks|lt) - [(kl|st) + (kt|sl)]
//
//
// 


class HaIntEngine
//! Class to generate 1-e and 2-e integrals on basis functions
{
public:

	HaIntEngine();
	HaIntEngine(const HaQCMod & host);
	virtual ~HaIntEngine();
	
//	static int CalcOvlpMat(HaBasisSet* pbas1, HaBasisSet* pbas2, HaMat_double& ovlp_mat);
//	static int Eval1eOp(HaBasisSet* pbas1, HaBasisSet* pbas2,HaMat_double& ovlp_mat, int oper_type);

	int Gen2eInt(HaMat_double & R1, HaMat_double & R2, HaMat_double & R3,
		const TwoElComb &Type2e) const; // load 2e integrals

	static int IntCanonToSq(double *pcan, HaMat_double & SqInt, int N);
	// convert 2e integral from canonic linear form to square (ij|kl) form

	const HaQCMod* ptr_qc_mod;
};






#endif /* !HAINTENGINE_H */
