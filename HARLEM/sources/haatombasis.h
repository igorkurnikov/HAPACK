/*!  \file haatombasis.h

    Classes to define Atomic Basis object in HARLEM.

    \author Igor Kurnikov  

    \date 1997-2003
*/
#ifndef HAATOMBASIS_H
#define HAATOMBASIS_H

#include "hastring.h"
#include "haio.h"
#include "halinalg.h"
#include "haatom.h"

class InternalBasis;
class AtomContainer;

struct b_type;

class HaBasisSet
//! Abstract class to describe a basis set of orbitals or configurations
{
public:
    virtual int GetNBfunc() const = 0;       //!< Return the number of basis functions in the set
	virtual std::string GetClassName() const = 0; //!< Return Class Name (Type)  
	virtual std::string GetLabel(int idx) = 0;  //!< Get label of the basis function ( 0-base index)
	virtual Vec3D* GetHostPt(int idx) = 0;   //!< Get 3D point the function is associated with ( 0-base index)
	virtual const Vec3D* GetHostPt(int idx) const = 0;   //!< Get 3D point the function is associated with ( 0-base index)
	virtual int TransferBetweenAtoms(PtrPtrMap& pt_corr_map) = 0; //!< Set new origin points for the basis set functions 

	static int CalcOvlpMat(const HaBasisSet* pbas1, const HaBasisSet* pbas2, HaMat_double& ovlp_mat); //!< compute overlap matrix for basis sets
    
	static std::string GetID(const HaBasisSet* pbas); //!< Generate String ID for the basis set
	static int RemoveCachedMatForBasis(const HaBasisSet* pbas); //!< Remove cached Matrices for the basis set
	static int ClearMatCache(); //!< Clear Cache of Overlap Matricies
	static HaMat_double* GetCachedOvlpMat(const HaBasisSet* pbas1, const HaBasisSet* pbas2); //!< Get Cached Ovelap Matrix between bases sets
    static int SaveInCacheOverlapMap(const HaBasisSet* pbas1, const HaBasisSet* pbas2, const HaMat_double& smat); //!< Put Overlap Matrix of two basis set to cache

protected:

	static StrPtrMap ovlp_map_cache;
};

class ArrayOrb3D : public HaBasisSet
//!< An abstract class to describe an array of orbitals in 3D space
{
public:
	int ExpandInBas(HaMat_double& coef, HaBasisSet& bset); //!< Find expansion coefficients of Array of orbitals in a given basis set

	virtual int GetTransfMat(HaMat_double& trans_mat, const HaMat_double& rot_mat); //!< Set Transformation matrix for coefficients of vectors expanded in the basis when the basis rotated by rot_mat
	HaMat_double GetTransfMat(const HaMat_double& rot_mat ); //!< Get Transformation matrix for expansion coefficients upon rotation
//	virtual double EvalLinCombInPoint(double x, double y, double z, const double* cf ) const = 0 ; //!< Evaluate a linear combination of function in the point
//	virtual double GetExtent(int i, double tol) const = 0; //!< Get function space extent for tolerance tol

	virtual TiXmlElement* AddXml(TiXmlElement* parent_element, const char* name = "", int option = 0) const = 0; //!< Add XML description of the orbitals as a child element of parent_element 
    virtual int SaveXML(FILE* file_out, int option=0) const;
	virtual int LoadXmlFile(FILE* file_inp, const char* tag_name = "", int option = 0);
	virtual int LoadXml(const TiXmlElement* xml_element, int option=0 ) = 0;
	static ArrayOrb3D* CreateObjectWithType(const char* type); //!< Create Empty Basis of a given type
};

class GauShell: public ArrayOrb3D 
//! Class to represent a group of Gaussian Basis functions with the same radial functions 
//! and different angular dependence 
{
public:
  GauShell();
  GauShell(int new_l_ang, const int NGauss);
  GauShell(int new_l_ang, const HaMat_double& cf_new);
//  GauShell(GauShell & OldShell);
  virtual ~GauShell();

  int GetL() const { return l_ang;} //!< Get Angular momentum of the shell
    
  const char* GetShellSymbol() const;        //!< Get a symbol corresponding to a shell type 
  const char* GetShellFunSymbol(int ifun); //!< Get Symbol of function ifun of the shell (0-based)

  virtual std::string GetLabel(int idx);  //!< Get Orbital Text Label (0-based idx)
  virtual Vec3D* GetHostPt(int idx){ return NULL; }  //!< Get 3D point the function is associated with
  virtual const Vec3D* GetHostPt(int idx) const{ return NULL; }  //!< Get 3D point the function is associated with

  virtual std::string GetClassName() const { return "GauShell"; } 
  virtual int GetNBfunc() const;     //!< return a number of basis functions in the shell

  virtual int GetTransfMat(HaMat_double& trans_mat, const HaMat_double& rot_mat); //!< Set Transformation matrix for coefficients of vectors expanded in the basis when the basis rotated by rot_mat
  virtual int TransferBetweenAtoms( PtrPtrMap& pt_corr_map); //!< Set new origin points for the basis set functions 
  virtual TiXmlElement* AddXml(TiXmlElement* parent_element,const char* name = "", int option = 0) const; //!< Add XML description of the orbitals as a child element of parent_element 
  virtual int LoadXml(const TiXmlElement* gau_shell_element, int option=0 );

  int GetNBfuncCart() const;         //!< Compute the number basis functions if converted to Cartesian Basis Functions 

  double EvalLinCombInPoint( double x, double y, double z, const double* cf) const; //!< Evaluate a linear combination of function in the point
  double GetExtent(int i, double tol) const; //!< Get function space extent for tolerance tol

  int GetNumGauss(void) const;             //!< return the number of contracted Gaussians in the shell
  bool SetNumGauss(const int NGauss);      //!< set the number of Contracted gaussians 
  
  bool SetCoef(const double* NewCoef);           //!< set gaussian exponents and contraction coefficients
  bool SetCoef(const HaMat_double& new_cf_mat);  //!< set gaussian exponents and contraction coefficients
  const HaMat_double& GetCoef() const { return coef;} //!< get gaussian exponents and contraction coefficients

  double GetExp( int ig) { return coef(1,ig); }   //!< Get Exponent of a Gaussian component ( 1-based index)
  double GetCoef( int ig ) { return coef(2,ig); } //!< Get Expansion coef of the Gaussian component ( 1-based index)

  void SaveGaussianInp(std::ostream& os) const;        //!< Save Shell description as an input for Gaussian 
  // Comparisons operators (for STD library)
  bool operator == (const GauShell & rhs) const;
  bool operator  < (const GauShell & rhs) const;

  int Normalize(); //!< Normalize shell coefficients

protected:

  int spherical; //!<  flag to use spherical ( D_x2-y2, D_z2 etc ) or coordinate Gaussian Basis functions
                         
  int l_ang;          //!< Angular momentum of the shell
  int NumGauss;       //!< number of elemental Gaussian Functions in Shell
  bool DestroyCoef();
  HaMat_double coef;  //!< an array of gaussian exponents and contraction coefficients coef(2,NumGauss) exp and expansion coef
                  
};

typedef std::vector<GauShell> ShellsType;

class HaPseudoPot;

class GauAtomBasis: public ArrayOrb3D 
//! Class for Gaussian Atomic Basis set 
//! currently identified by the name of the basis set (like "6-31G")
//! and the nuclear charge of the atom. 
{

  friend class HaQCMod;

public:
  GauAtomBasis();
  GauAtomBasis(const std::string& NewBasName, const std::string& NewAtomType);
  GauAtomBasis(const std::string& NewBasName, HaAtom* aptr );
  GauAtomBasis(const GauAtomBasis & ref);
  virtual ~GauAtomBasis();

  void SetDefaultParams();

  GauAtomBasis& copy_from(const GauAtomBasis & ref);    //!< Copy content from a reference
  int SetForAtom(const char* BasName, HaAtom* aptr ); //!< Set Basis set on an atom 

  const char* GetBasName() const;            //!< Get Atomic Basis Name
  bool SetBasName(const std::string & name);    //!< Set Atomic Basis Name

  bool           SetAtHost(HaAtom* new_host_atom);
        HaAtom*  GetAtHost();
  const HaAtom*  GetAtHost() const;

  virtual std::string GetLabel(int idx); //!< Get Orbital label
  virtual Vec3D* GetHostPt(int idx){ return host_atom; }  //!< Get 3D point the function is associated with
  virtual const Vec3D* GetHostPt(int idx) const { return host_atom; }  //!< Get 3D point the function is associated with
  virtual int TransferBetweenAtoms(PtrPtrMap& pt_corr_map);    //!< Set new origin points for the basis set functions 
  virtual TiXmlElement* AddXml(TiXmlElement* parent_element, const char* name = "",int option = 0) const; //!< Add XML description of the orbitals as a child element of parent_element 
  virtual int LoadXml(const TiXmlElement* gau_shell_element, int option=0 );

  bool SetFromGaussianInp(std::istream& is);
  void SaveGaussianInp(std::ostream& p_stream_out) const;
  
  std::string GetAtomType() const;
  bool SetAtomType(const std::string& atype); 
  virtual int GetNBfunc() const;      //!<  return the number of Basis functions in that the AtomBasis
  int GetNBfuncCart() const;  //!< Compute the number basis functions if converted to Cartesian Basis Functions 

  virtual std::string GetClassName() const { return "GauAtomBasis"; } //!< Return Class Name 

  bool AddShell(GauShell & shl);         //!< Add Gaussian Shell to the Basis
//  int Normalize(); //!< Normalize coefficients of all shells of the Atomic Basis
  
  void Clear();                         //!< delete all class data
  void ClearCoef() { Shells.clear(); }  //!< delete Coef matricies

  void SetPseudoPotName(const std::string& new_pot_name);
  bool SetPseudoPotFromName(); //!< Set PseudoPotential From DataBase using 
                               // Name of the PseudoPotential and the name of the atom
  void SetPseudoPotPtr(const HaPseudoPot* new_ppot);
  const HaPseudoPot* GetPseudoPot() const { return ppot; } 
  bool IsSetPseudoPot() const { return (ppot != NULL); }

  int GetNumElectr() const; //!<  Get Number of active electrons of the atom of the basis 
                           //!< taking into account effective core models 

  // Comparisons operators (for STL ): 
  bool operator == (const GauAtomBasis & rhs) const;
  bool operator <  (const GauAtomBasis & rhs) const;

  bool Print_info(std::ostream &sout, const int level) const;

protected:
  HaAtom* host_atom;       //!< Atom - Basis Set reside on	 
  int internal_atom_flag;  //!< if TRUE host_atom is internal to the Atom Basis and will be destroyed upon Atom Basis set destruction
  std::string BasName;        //!< atomic basis set name 
  std::string AtomType;       //!< Atom Std Label
  std::string PseudoName;     //!< name of the Pseudopotential associated with the basis
  const HaPseudoPot* ppot; //!< Pointer to the Pseudo Potential class 
public:
  std::vector<GauShell> Shells; //!< list of Gaussian shells forming the basis
};


class GauBasisSet:  public ArrayOrb3D
//! Class to represent Gaussian Basis Set in HARLEM 
{
public:
	GauBasisSet();
	virtual ~GauBasisSet(); 

	void Clear();

	typedef std::vector<GauAtomBasis>::iterator AtomBasIterator;

	std::string GetName() const; //!< Get Name of the Gaussian Basis Set  ("GEN" if not described by a single name)
	bool IsGeneric() const;  //!< Check if basis set should be considered generic
	void SetGeneric();       //!< Set basis set to be generic

	int LoadToGaussianBas(b_type& gaub) const; //!< Load Basis info to Gaussian basis structure
	int LoadToGaussianBCommon() const;  //!< Load /B/ GAUSSIAN Common with Basis Set info 
	int LoadToGaussianB2Common() const; //!< Load /B2/ GAUSSIAN Common with Basis Set info 
	InternalBasis* CreateIPackBas();   //!< Create and return a pointer to IPACK InternalBasis object with Basis Set info

  	virtual int GetNBfunc() const;     //!< Compute the number of basis functions in the basis set
	virtual std::string GetClassName() const { return "GauBasisSet"; } 
    int GetNBfuncCart() const;         //!< Compute the number basis functions if converted to Cartesian Basis Functions 
	int pure_fun_flag;                 //!< if !=0 pure functions are assumed otherwise cartesian function are used 

	virtual std::string GetLabel(int idx); //!< Get Orbital Text Label
	virtual Vec3D* GetHostPt(int idx);  //!< Get 3D point the function is associated with
	virtual const Vec3D* GetHostPt(int idx) const;  //!< Get 3D point the function is associated with

    virtual int GetTransfMat(HaMat_double& transf_mat, const HaMat_double& rot_mat); //!< Set Transformation matrix for coefficients of vectors expanded in the basis when the basis rotated by rot_mat
    virtual int TransferBetweenAtoms(PtrPtrMap& pt_corr_map); //!< Set new origin points for the basis set functions 
    virtual TiXmlElement* AddXml(TiXmlElement* parent_element, const char* name = "", int option = 0) const; //!< Add XML description of the orbitals as a child element of parent_element 
    virtual int LoadXml(const TiXmlElement* xml_element, int option=0 );
	
	int InitForMolSet(const char* bname, MolSet* pmset); //!< Init Basis for all atoms of the Molecular Set
    int InitForAtoms(const char* bname, AtomContainer* at_coll); //!< Init Basis for all atoms of the Atom Collection
	GauAtomBasis* AddBasisToAtom(const char* bas_name, HaAtom* aptr); //!< Add BasisSet To Atom

	GauAtomBasis& GetAtBasByIdx(int i);  //!< Get atomic basis set by index
	int GetNumAtBas() const;             //!< Get the number of atomic basis set making given Gaussian Basis set

	static int CalcOvlpMat(GauBasisSet* pbas1, GauBasisSet* pbas2, HaMat_double& ovlp_mat); //!< Compute overlap matrix between two basis sets (may be the same)

	int GetNumCnt() const;                      //!< GetNumber of Centers
	int GetCntCoord(HaMat_double& coord) const; //!< Get array of coordinates of orbitals centers
	int GetCntCoordArr(Vec3DValArray& crd_arr) const; //!< Get array of coordinates of orbitals centers
	int GetNumElectr() const;                   //!< Get the number of electrons assuming tot mol charge=0, 

	int GetAtBasIdxForOrb(int i);     //!< Get an index of atom basis for a given basis function (0-based)
	int RecompFstBasVec();            //!< Recompute an array of first basis functions (fst_bas_fun_idx) of atomic basis sets 

    int MatchBasisSet(const GauBasisSet* basis_frag, IntIntMap& frag_bas_fun_map, HaVec_double* bas_pert_vec = NULL); //!< Map Basis functions of the fragment to the current basis set, bas_pert_vec - perturbatins factors for fragment atom basis

	std::vector<GauAtomBasis> at_bas_vec;           //!< vector of atom bases    
	HaVec_int atom_bas_idx;                    //!< Atom basis indexes of basis functions 
	HaVec_int fst_bas_fun_idx;                 //!< Index of the first orbitals of the atomic basis sets in the GauBasSet
	StrVec bf_lbls;                            //!< Labels of basis functions

protected:
	std::string bset_name; 

};


//class SlaterAtomBasis
//! Class for Slater Atomic Basis set 
//{
//public:
//	SlaterAtomBasis();
//	virtual ~SlaterAtomBasis();

//	bool          SetHostPt(Vec3D* new_host);
//        Vec3D*  GetHostPt();
//	const Vec3D*  GetHostPt() const;

//    int GetNBfunc() const;      //!<  return the number of Basis functions in that the Atom Basis

//	int norb;
//	int maxl;
//	int nel;
//	double zetas;
//	double zetap;
//	HaVec_double zetad;
//	HaVec_double cfd;

//	Vec3D* pt_host;
//};

//class SlaterBasisSet:  public ArrayOrb3D
//! Class to represent Slater Orbital Basis Set in HARLEM
//{
//public:
//	SlaterBasisSet();
//	virtual ~SlaterBasisSet(); 

//	vector<SlaterAtomBasis> at_bas_vec;           //!< vector of atom bases    

//	static int CalcOvlpMat(SlaterBasisSet* pbas1, SlaterBasisSet* pbas2, HaMat_double& ovlp_mat); //!< Compute overlap matrix between two basis sets (may be the same)
	
//  int GetNBfunc() const;      //!<  return the number of Basis functions in that the basis set

//};

#endif /* !HAATOMBASIS_H */
