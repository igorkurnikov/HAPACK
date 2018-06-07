/*!  \file harpaham.h
 
    Classes to define RPA Hamiltonian in HARLEM.

    \author Igor Kurnikov  
    \date 1998-2002

*/
#ifndef HARPAHAM_H
#define HARPAHAM_H

#include "harpavec.h"


class HaRPAHam 
//! Class to define CIS or RPA Hamiltonian 
{
public:
	HaRPAHam();
	virtual ~HaRPAHam();

	bool SetOpMode(int imode); //!< Set Mode of Hamiltonian action
	bool SetEnergy(const double NewEne); //!< Set Energy if Calc (E-H)


bool Apply(vector<HaRPAvec> & RPAvec); //!< Apply RPA Hamiltonian to the array
                                       //!<  of RPA vectors defined as 

bool Apply_init(vector<HaRPAvec> & RPAvec) const; 
      // Apply zero order RPA Hamiltonian 


vector<HaRPAvec> 
operator *(const vector<HaRPAvec> & RPAvec); //!< Apply Hamiltonian

protected: 
	double Ene;
	enum OperMode {FULL=0,MO_DIAG=1,E_MIN_H=2,E_MIN_H0=3} opmode; 
};


class HaRPAResolv 
//! Class to define CIS or RPA Resolvent operator 
{
public:
	HaRPAResolv();
	virtual ~HaRPAResolv();


bool SetOpMode(const int imode);

bool SetEnergy(const double NewEne); 

bool SetImag(const bool new_imag);


bool Apply_Init(vector<HaRPAvec> & RPAvec) const; //!< Apply Zero order approximation
                                                  //!< (diagonal MO) RPA resolvent

bool Apply(vector<HaRPAvec> & RPAvec) const; //!< Apply RPA Resolvent 
                                             //!< to the array of RPA vectors defined as 

vector<HaRPAvec> solve(const vector<HaRPAvec> & CISvec) const; //!< Apply RPA Resolvent 
                                                               //!< in the given approximation

protected: 
	bool imag;    //!< flag to indicate that RPAvec is imaginary
	double Ene;   //!< Energy for the Green Functions
	enum OperMode {FULL=0, MO_DIAG} opmode; 
	                      //!< Code of approximation method used 
	                      //!< to calculate resolvent matrix elements
};



#endif /* !HARPAHAM_H */
