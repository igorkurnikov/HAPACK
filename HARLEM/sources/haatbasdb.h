/*! \file haatbasdb.h
    
    interface to HaAtBas and HaAtBasDB class
    Classes to manage Atomic basis set DataBase.

    \author Igor Kurnikov , University of Pittsburgh 
    \date 1998-2002

*/
#ifndef HAATBASDB_H
#define HAATBASDB_H

#include "hastl.h"
#include "haatombasis.h"


class HaAtBasDB 
//!< Class to define database of quantum chemical gaussian basis sets in HARLEM 
{
public:
  HaAtBasDB(void);
  virtual ~HaAtBasDB(void);
 
  bool Init();
  bool AddAtomBasis(const GauAtomBasis& atbas); 
  GauAtomBasis* Extract(const std::string & BName, const std::string & Atomlbl);

  // Axxiliary functions:
  bool Init321G();        //!< Initialize 3-21G basis set entry 
  bool Init321PPG();      //!< Initialize 3-21++G basis set entry 
 
  bool InitHay1_dz();     //!< Initialize 3-21G for Light Atoms and Hay-Wadt dz basis for Heavy Atoms
                          //!< with Pseudo-Potential 
  bool InitHay1_dzPP();   //!< Initialize 3-21++G for Light Atoms and Hay-Wadt dz basis for Heavy Atoms
                          //!< with Pseudo-Potential 

  bool SetFromDaltonFile(const std::string & BName, const std::string & Atomlbl);
 

protected:
  StrPtrMap bas_map;
};



#endif /* !HAATBASDB_H */
