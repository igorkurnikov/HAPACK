/*! \file halocorb.h

    Classes to define Local Orbital object in HARLEM.

    \author Igor Kurnikov  
    \date 1997-2004  
*/
#ifndef HALOCORB_H
#define HALOCORB_H

#include "hastl.h"
#include "haatombasis.h"

class HaQCMod;

class LinCombOrb3D:  public ArrayOrb3D
//! An Array of 3D orbitals as linear combination of basis set functions 
{
public:
    LinCombOrb3D();
    LinCombOrb3D(ArrayOrb3D* new_bas, const double* coef_new, int nv);
	virtual ~LinCombOrb3D();

	void SetStdParams();  //!< Set Standard Parameters to Data

	virtual std::string GetClassName() const { return "LinCombOrb3D"; }

	bool IsEmpty();  //!< Check if the collection doesn't contain any orbitals
	void Clear();  //!< Delete All Orbitals in the collection

	int GetOrbIdxByID(const char* id); //!< Return zero based index of orbital by ID
	virtual std::string GetLabel(int idx);  //!< Get Label of Orbital with given index
    virtual Vec3D* GetHostPt(int idx);  //!< Get 3D point the function is associated with
	virtual const Vec3D* GetHostPt(int idx) const;  //!< Get 3D point the function is associated with
	virtual int TransferBetweenAtoms(PtrPtrMap& pt_corr_map); //!< Set new origin points for the basis set functions 
    virtual TiXmlElement* AddXml(TiXmlElement* parent_element,const char* name = "", int option = 0) const; //!< Add XML description of the orbitals as a child element of parent_element 
    virtual int LoadXml(const TiXmlElement* xml_element, int option = 0); 

	int CreateEmptyOrbs(int n_orb,ArrayOrb3D* new_bas); //!< Create n empty orbitals in the array and set basis set
	int AddOrbs(LinCombOrb3D& orbs);  //!< Add a linear combination of orbitals to the array
	void SetOrbLabel(int idx, const std::string& orb_lbl); //!< Set Label of the orbital in the array
	
    int ProjectToBasis(HaMat_double& coef_new, const ArrayOrb3D* basis_new); //!< project orbitals to a new basis set

	int GetNOrbs() const;    //!< Get Number of orbitals
	virtual int GetNBfunc() const{ return GetNOrbs(); }  //!< Get the number of orbitals in the set

    int TrCoefRot(const HaMat_double& rot_mat); //!< Transform expansion coefficients for atom rotation 

	int internal_basis;  //!< if TRUE delete basis on destruction (default = FALSE)

	HaMat_double coef;  //!< expansion coefficients - number of columns = number of vectors in the array
	ArrayOrb3D* bas;    //!< the basis set for orbitals 

	StrVec ids;    //!< String identificators of orbitals
	VecPtr at_ptr; //!< Pointers to atoms that orbitals are located on

//! Compute overlap matrix between to arrays of orbitals 
	static int CalcOvlpMat(LinCombOrb3D* pbas1, LinCombOrb3D* pbas2, HaMat_double& ovlp_mat);
//! Evaluate 1-e operator on arrays of orbitals using operator matrix on basis functions
	static int Eval1eOp(LinCombOrb3D* pbas1, LinCombOrb3D* pbas2, const HaMat_double& bop_mat, HaMat_double& op_mat);

};

#endif /* !HALOCORB_H */
