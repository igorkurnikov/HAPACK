//  stmmod.h
//
//  HARLEM
//
//  Classes 
//  to model STM current  
//
//  definitions
//
//  Igor Kurnikov , University of Pittsburgh 
//  Jianping Lin, Duke University
//
//  Created:   October 26 2001
//

#ifndef STMMOD_H
#define STMMOD_H

#include "hastl.h"
#include "command.h"
#include "hacompmod.h"


class StmMod : public HaCompMod
{
public:
	StmMod(HaMolSet* new_phost_mset = NULL);
	virtual ~StmMod();

	void SetStdParams();
	
	int CalcTMatr1();

	int num_mol_ao;
	double tmat_el;
	double wave_ve_1;
	double wave_ve_2;
	double wave_ve_3;
	
protected:

};

#endif // end if !defined(STMMOD_H) 
