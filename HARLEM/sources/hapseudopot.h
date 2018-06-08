// HARLEM
// Igor Kurnikov
// University of Pittsburgh
//
// hapseudopot.h: interface for the classes to manage core pseudopotentials
// 
// Created: December 1998
// Revised: December 1998
//////////////////////////////////////////////////////////////////////

#ifndef HAPSEUDOPOT_H
#define HAPSEUDOPOT_H

#include "hastl.h"
#include "hastring.h"


class PseudoTerm
{
public:

	PseudoTerm();
	PseudoTerm(int new_npower, double new_expon, double new_coef);
	~PseudoTerm();

	int npower;
	double expon;
	double coef;
};

class PseudoBlock
{
public:
	PseudoBlock();
	~PseudoBlock();

	std::string description;
	vector<PseudoTerm> terms;

};

class HaPseudoPot
//! class to define Core Pseudopotentials 
{
public:
  HaPseudoPot();
  virtual ~HaPseudoPot();

  bool SaveGaussInp(std::ostream& os) const;
  int GetNCoreEl() const { return ncore_el; } 
  
  friend class HaPseudoPotDB;
protected:

  std::string pot_name;
  std::string at_symbol;
  unsigned int max_angm;
  unsigned int ncore_el;
     
  vector<PseudoBlock> blocks;
};

class HaPseudoPotRef
{
public:
	HaPseudoPotRef();
	HaPseudoPotRef(const std::string& new_pot_name, const std::string& new_at_label);
	virtual ~HaPseudoPotRef();
	
	friend class HaPseudoPotDB;

  bool operator == (const HaPseudoPotRef &  rhs) const ;
  bool operator <  (const HaPseudoPotRef &  rhs) const ;

protected:

	std::string pot_name;
	std::string at_label;

};


class HaPseudoPotDB
{
public:
	HaPseudoPotDB();
	~HaPseudoPotDB();

	HaPseudoPot* Extract(const std::string& pot_name, const std::string& at_label); 
	bool Init();

	
protected:
	map<HaPseudoPotRef,HaPseudoPot*, less<HaPseudoPotRef> > dat; // Atomic Orbital symbol

	bool InitHayWadt_1(); // Initialize Hay and Wadt Pseudopotential 1 (only valence electrons)
};


#endif /* !HAPSEUDOPOT_H */
