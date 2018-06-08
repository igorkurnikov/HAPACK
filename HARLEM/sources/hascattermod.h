/*! \file hascattermod.h

   Classes to model electronic scattering in HARLEM  

   \author Igor Kurnikov  
   \date 2000-2002
*/

#ifndef HASCATTERMOD_H
#define HASCATTERMOD_H

#include "hastl.h"
#include "command.h"
#include "hacompmod.h"

class HaNonLocField3D;


class HaScatterMod : public HaCompMod
{
public:
	HaScatterMod(HaMolSet* new_phost_mset = NULL);
	virtual ~HaScatterMod();

	void SetStdParams();

	int FindGridHamEigVec(int num_vec); // find lowest num_vec 
	                       // eigen values of the grid hamiltonian matrix
	
	int CalcPseudoPotGrid();
	int CalcPseudoPotAO();
	double GetPsPGridElem(  int nr_1, int nc_1, int nl_1, 
						  int nr_2, int nc_2, int nl_2); 
	double  GetPsP_xyz(double x1, double y1, double z1,
						 double x2, double y2, double z2);
	double GetPsP_Gauss_xyz(double x1, double y1, double z1,
						    double x2, double y2, double z2);

	double GetKinEne_Gauss_xyz(double x1, double y1, double z1,
						    double x2, double y2, double z2);

	int SetCoreHamOnGrid();

	HaNonLocField3D* pot_grid;
	double psp_thresh_val; // Threshold value for the matrix elements of the 
	                       // pseudopotential below which this element is discarded
	double pt_exp;         // Exponent for the gaussian basis function located on the grid (bohr^-2)
	int calc_kin_ene_flag; // calculate kinetic energy instead of pseudopotential  
	std::string psp_fname;
	
protected:
	HaMat_double pseudo_pot_ao;	// matrix of the pseudo potential 
	                            // in the Gaussian Atomic basis set
};

#endif // end if !defined(HASCATTERMOD_H) 
