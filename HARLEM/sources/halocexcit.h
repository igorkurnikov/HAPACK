//  halocexcit.h
//
//  Classes to define local excitation in HARLEM
//
//
//  Igor Kurnikov , University of Pittsburgh 
//
//  Revisions: January 19 1997
//
#ifndef HALOCEXCIT_H
#define HALOCEXCIT_H

#include "halocorb.h"
#include "haatgroup.h"

class HaQCMod;

class HaLocExcit
// 
// class to define Local Excitation ( particle-hole pair )
// in the given Atomic Basis
// 
{
public:
  HaLocExcit();
  HaLocExcit(ArrayOrb3D & p1, ArrayOrb3D & p2);
  virtual ~HaLocExcit();

// Expand Local Excitation to transition density matrix
  int Expand(double* dmatr, HaQCMod & host) const;
 
protected:
  ArrayOrb3D* PartOrb;
  ArrayOrb3D* HoleOrb;

};

enum CISex_store_mode {AO_EXCIT,MO_EXCIT,AO_EXCIT_T};
//  Store mode of Excitation array starage
//  AO_EXCIT -   Array of transition densities 
//               premultiplied by S^(-1) (standard in G94)
//  AO_EXCIT_T - True AO transition density as obtained 
//               from transformation of CIS vector in MO
//  MO_EXCIT - array of amplitudes of occ_MO -> vac_MO 
//             excitations
//



#endif /* !HALOCEXCIT_H */
