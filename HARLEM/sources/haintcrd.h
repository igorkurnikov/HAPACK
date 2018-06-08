/*! \file haintcrd.h

    Classes to define internal coordinates manipulations in HARLEM

    \author Igor Kurnikov 
 
    \date 2002

*/
#ifndef HAINTCRD_H
#define HAINTCRD_H

#include "haatgroup.h"
#include "hacoord.h"

class CrdAssignRule
//! Base class to define rules for setting atoms coordinates  
{
public:
	CrdAssignRule();
	virtual ~CrdAssignRule() {}

	int GetPriority() const;       //!< Get relative priority of execution for the rule ( smallest executed first)
	void SetPriority( int priority_new ); //!< Get relative priority of execution for the rule 

	virtual int SaveXMLToStream(std::ostream& os, const harlem::HashMap* popt = NULL ) const { return TRUE; } 

	virtual int SetManagedAtomCrd() = 0; //!< Set coordinates of the managed atom(s) according to the rule 

protected:

	int priority; //!< Priority to apply the synchronization rule ( smallest executed first )
};

class SingleAtomCrdRule : public CrdAssignRule
{
//!< Rule to set coordinates of a single atom
public:
	SingleAtomCrdRule( HaAtom* p_mng_atom_par );                
	virtual ~SingleAtomCrdRule();

	HaAtom* GetManagedAtom(); //!< Get atom managed by the rule
	int SetManagedAtom( HaAtom* aptr ); //!< Set atom managed by the rule

	virtual int SaveXMLToStream(std::ostream& os, const harlem::HashMap* popt = NULL ) const; 
	virtual int ReplaceRefPt( Vec3D* p_ref_old, Vec3D* p_ref_new) = 0; //!< Replace one of the reference points

protected:
	HaAtom* p_mng_atom;
};


class FixedCrdRule : public SingleAtomCrdRule
//!< Cartesian Coordinates of the Managed Atom are set to fixed values 
{
public:
	FixedCrdRule(HaAtom* p_mng_atom_par );
	virtual ~FixedCrdRule();

	virtual int SetManagedAtomCrd(); 

	virtual int SaveXMLToStream(std::ostream& os, const harlem::HashMap* popt = NULL ) const; 
	virtual int ReplaceRefPt( Vec3D* p_ref_old, Vec3D* p_ref_new);

	void SetCrd( const Vec3D& crd); //! Set XYZ coordinate values
	void SetParFromCurPos(); //!< Set coordinate values of the rule from current managed atom position

	void SetX( double val );
	void SetY( double val );
	void SetZ( double val );
	double GetX() const;
	double GetY() const;
	double GetZ() const;

private:
	Vec3D crd; 
};

class SameAtomCrdRule : public SingleAtomCrdRule
//!< Atom coordinate synchronization rule: Managed Atom coordinates are set to those of the reference atom  
{
public:
	SameAtomCrdRule(HaAtom* p_mng_atom, HaAtom* p_ref_atom);
	virtual ~SameAtomCrdRule();

	virtual int SetManagedAtomCrd();

	void SetX( double val );
	void SetY( double val );
	void SetZ( double val );

	virtual int SaveXMLToStream(std::ostream& os, const harlem::HashMap* popt = NULL ) const; 

	HaAtom* GetRefAtom() { return p_ref_atom; } //!< Get Reference Atom of rule

	virtual int ReplaceRefPt( Vec3D* p_ref_old, Vec3D* p_ref_new); //!< Replace one of the reference points

private:
	HaAtom* p_ref_atom; //!< Reference Atom to synchronize coordinates
};

class Pt2CrdRule : public SingleAtomCrdRule
//!< Atom coordinate synchronization rule: Managed Atom coordinates are set to lay on a line between two reference atoms  
{
public:
	Pt2CrdRule(HaAtom* p_mng_atom, HaAtom* p_ref_1, HaAtom* p_ref_2);
	virtual ~Pt2CrdRule();

	virtual int SetManagedAtomCrd();

	HaAtom* GetRefAtom1() { return p_ref_1; } //!< Get Reference Atom 1
	HaAtom* GetRefAtom2() { return p_ref_2; } //!< Get Reference Atom 2

	virtual int ReplaceRefPt( Vec3D* p_ref_old, Vec3D* p_ref_new); //!< Replace one of the reference points

	void SetBondLen(double bond_len_new); //!< Set Distance between p_host_atom and p_ref_1

	virtual int SaveXMLToStream(std::ostream& os, const harlem::HashMap* popt = NULL ) const; 

private:
	HaAtom* p_ref_1; //!< Reference Atom 1 
	HaAtom* p_ref_2; //!< Reference Atom 2

	double bond_len; //!< Distance between p_mng_atom and p_ref_1
};

class Pt3CrdRule : public SingleAtomCrdRule
//!< Atom coordinates synchronization rule: Managed Atom coordinates are set from 3 points positions, a bond length, a val angle and a torsional angle 
{
public:
	Pt3CrdRule(HaAtom* p_mng_atom, Vec3D* p_ref_1, Vec3D* p_ref_2, Vec3D* p_ref_3);
	virtual ~Pt3CrdRule();

	virtual int SetManagedAtomCrd();

	int SetParFromPtPos(const Vec3D* pt_mng, const Vec3D* pt_ref_1, const Vec3D* pt_ref_2,const Vec3D* pt_ref_3); //!< Set Params of the rule (bond_len,..etc from point coordinates) 
	int SetParFromCurPos(); //!< Set Params of the rule (bond_len,..etc from current reference points coordinates) 
	int SetBondLen(double bond_len); //!< Set Distance between p_host_atom and p_ref_1
	int SetValAng (double val_ang);  //!< Set Valence angle between p_host_atom, p_ref_1 and p_ref_2
	int SetDihAng (double dih_ang);  //!< Set Dihedral angle betweem p_host_atom, p_ref_1, p_ref_2 and p_ref_3

	double GetBondLen() const;   //!< Get Distance between p_host_atom and p_ref_1
	double GetValAng()  const;   //!< Get Valence angle between p_host_atom, p_ref_1 and p_ref_2
	double GetDihAng()  const;   //!< Get Dihedral angle betweem p_host_atom, p_ref_1, p_ref_2 and p_ref_3

	Vec3D* GetRefPt1();  //!< Get Reference point 1 ( defining bond length) 
	Vec3D* GetRefPt2();  //!< Get Reference point 2 ( defining valence angle ) 
	Vec3D* GetRefPt3();  //!< Get Reference point 3 ( defining dihedral angle ) 

	virtual int ReplaceRefPt( Vec3D* p_ref_old, Vec3D* p_ref_new); //!< Replace one of the reference points

	virtual int SaveXMLToStream(std::ostream& os, const harlem::HashMap* popt = NULL ) const; 

private:

	Vec3D* p_ref_1;
	Vec3D* p_ref_2;
	Vec3D* p_ref_3;

	double bond_len; //!< Distance between p_mng_atom and p_ref_1
	double val_ang;  //!< Valence angle between p_mng_atom, p_ref_1 and p_ref_2
	double dih_ang; //!< Dihedral angle betweem p_mng_atom, p_ref_1, p_ref_2 and p_ref_3
};

class DihedralAngleCoord : public harlem::Coord
//! Internal coordinate node to specify dihedral angle  
{
public:
	
	DihedralAngleCoord( HaAtom* new_aptr1, HaAtom* new_aptr2, HaAtom* new_aptr3, HaAtom* aptr4);
    virtual ~DihedralAngleCoord();

// Virtual from HaCoord:
	virtual std::string GetClassName() const { return "DihedralAngleCoord";}  //!< Get Class Name of the Coordinate

	virtual harlem::Coord* clone();                     //!< Get a copy of coordinates
	virtual HaVec_double AsVecDouble() const;     //!< Transform Coordinates to a Vector of double values 
	virtual int SetFrom(const harlem::Coord* pcrd);     //!< Set Coordinates from the other coordinate object   
	virtual int SetFromVecDouble(const HaVec_double& dbl_vec);      //!< Set Coordinates from a vector of double values
	virtual int LoadFromStream(std::istream& is, const harlem::HashMap* popt = NULL );       //!< Read Coordinates from stream
	virtual int SaveToStream(std::ostream&  os,  const harlem::HashMap* popt = NULL ) const; //!< Write Coordinates to stream

	HaAtom* aptr1; //!<  1-st atom of the dihedral 
	HaAtom* aptr2; //!<  2-nd atom of the dihedral ( or 1st atom of the bond or valence angle)
	HaAtom* aptr3; //!<  3-rd atom of the dihedral ( or 2nd atom of the bond or valence angle)
	HaAtom* aptr4; //!<  4-th atom of the dihedral ( or 3nd atom of the valence angle)

	double GetDihVal();                //!< Compute Current Dihedral value 
	int SetDihVal(double new_dih_val); //!< Set new dihedral value for an internal coordinate

	int FindMovingAtoms();  //!< find atoms that move when internal coordinates change

	AtomGroup moving_atoms;   //!< Atoms that move when internal coordinates change

};

enum ElemCrdType { UNDEF_ELEM_CRD, X_ELEM_CRD, Y_ELEM_CRD, Z_ELEM_CRD, 
                   LEN_ELEM_CRD, ANG_ELEM_CRD, DIH_ELEM_CRD };
  
class ElemCrd
{
public:
	ElemCrd();
	ElemCrd( ElemCrdType type, CrdAssignRule* p_rule );
	virtual ~ElemCrd();
	
	ElemCrdType GetType() const { return type; }
	CrdAssignRule* GetCrdRule() { return p_rule; }
	
	void   SetValue( double val); //!< Set Value of the coordinate in storage units (Ang for X,Y,Z crd and bonds, radians for angles)
	double GetValue() const;  //!< Get Value of the coordinate in storage units (Ang for X,Y,Z crd and bonds, radians for angles)

	void SetDisplayValue( double val ); //!< Get Value of the coordinate in display units (Ang for X,Y,Z crd and bonds, degrees for angles)
	double GetDisplayValue() const; //!< Get Value of the coordinate in display units (Ang for X,Y,Z crd and bonds, degrees for angles)
	
	HaAtom* GetManagedAtom(); //!< Get Atom associated with elemental coordinate
	CrdAssignRule* GetCrdAssignRule(); //!< Get Atom Coordinate assing rule associated with the coordinate

	void SetTag( const std::string& tag_new); //!< Set Tag for the coordinate  
	std::string GetTag() const;               //!< get string tag of the coordinate

	void SetFrozen( bool set_flag = true ); //!< Set Coordinate Frozen flag on/off 
	bool IsFrozen() const;                //!< Check if coordinate if frozen

protected:
	bool frozen_flag;
	std::string tag;
	CrdAssignRule* p_rule;
	ElemCrdType type;
};

class ZMatLoadOptions : public harlem::HashMap
{
public:
	ZMatLoadOptions();
	ZMatLoadOptions( const ZMatLoadOptions& ref );
	virtual ~ZMatLoadOptions();

	virtual void Copy( const harlem::HashMap& ref );
	virtual harlem::HashMap* clone() const; 

	void SetStdOptions(); //!< Set standard options values

	void SetLoadAtID(bool set_flag = true); //!< Set Z-matrix lines expected to start with atom ID (sequence number or atom string reference )  
	bool ToLoadAtID() const;                  //!<  Check if Z-matrix lines expected to start with atom ID (sequence number or atom string reference )   

	void SetLoadAtElem(bool set_flag = true); //!< Set Z-matrix lines expected to contain atom element ( element number or Standard symbol )  
	bool ToLoadAtElem() const;                  //!< Check if Z-matrix lines expected to contain atom element ( element number or Standard symbol )

protected:

	bool load_at_id;  
	bool load_at_elem;

};

class ZMatSaveOptions : public harlem::HashMap
//!< Class with options to write Z-matrix 
{
public:
	ZMatSaveOptions(); 
	ZMatSaveOptions( const ZMatSaveOptions& ref );
	virtual ~ZMatSaveOptions();

	virtual void Copy( const harlem::HashMap& ref );
	virtual harlem::HashMap* clone() const; 

	void SetStdOptions(); //!< Set standard options values

	void SetSaveAtSeqNum(bool set_flag = true); //!< Set writing atom sequential number at the beginning of atom lines on/off  
	bool ToSaveAtSeqNum() const;                  //!< Check if to write atom sequential number at line start 

	void SetSaveAtElem(bool set_flag = true); //!< Set write atom element number on/off 
	bool ToSaveAtElem() const;                  //!< Check if to write atom element number

	void SetSaveAtSymbol(bool set_flag = true); //!< Set write atom symbol on/off 
	bool ToSaveAtSymbol() const;                  //!< Check if to write atom symbol

	void SetSaveTags(bool set_flag = true); //!< Write tags in place of tagged elemental coordinates
	bool ToSaveTags() const; //!< Check if Write tags for elemental tagged variables

protected:

	bool save_at_seq_num;  
	bool save_at_elem;
	bool save_at_symbol;
	bool save_tags;

};

class ZMatCrd 
{
public:
	ZMatCrd(HaMolSet* pmset_new);
	virtual ~ZMatCrd();

	void Clear(); //!< Clear content
	bool IsEmpty() const; //!< Check if empty
	
	void InitStdZMat();  //!< Init Standard Z-matrix for atoms 
	void InitAllXYZ();   //!< Init all XYZ coordinates for atoms

	int OnDelAtoms( AtomContainer& del_atoms ); //!< Modify Z-matrix upon atoms deletion

	int LoadFromString( const std::string& str, const harlem::HashMap* popt = NULL ); //!< Init Z-matrix from string
	int LoadFromStream(std::istream& is, const harlem::HashMap* popt = NULL ); //!< Init Z-matrix from stream
	
	int SetAtomCrd(); //!< Set Cartesian atom coordinates of the molecular set associated with Z-matrix 
	int GetCrdSnapshot( HaVec_double& crd_arr ); //!< Get array of cartesian coordinates of atoms corresponding to values of Z-matrix coordinates
	int SetFromAtomCrd(); //!< Set Z-matrix ( bond lengths, angles values etc) from current atom coordinates of the associated Moleculat Set 
	int SetFromCrdSnapshot( const HaVec_double& crd_arr ); //!<  Set Z-matrix from an array of cartesian atom coordinates 
	int GetElemCrdVal( HaVec_double& elem_crd_val_arr , bool unfrozen = true ) const; //!< Get Array of internal (elemental) coordinates values ( only unfrozen if unfrozen = true )   
	int SetFromElemCrdVal( const HaVec_double& elem_crd_val_arr , bool unfrozen = true ); //!< //!< Set Array of internal (elemental) coordinates values ( only unfrozen if unfrozen = true ) 
	int TransDerivToIntCrd( const HaVec_double& deriv_cart_crd, HaVec_double& deriv_int_crd, bool unfrozen = true ); //!< Transform Array of derivative in cartesian coordinates to internal coordinates 

	int SaveToStream(std::ostream& os, const harlem::HashMap* popt = NULL ); //!< Write Z-matrix to stream  
	std::string SaveToString(const harlem::HashMap* popt = NULL); //!< Write Z-matrix to string

	int SaveXMLToStream(std::ostream& os, const harlem::HashMap* popt = NULL ); //!< Save Z-matrix in XML format
	int LoadXMLNode( rapidxml::xml_node<>* node_zmat, const harlem::HashMap* popt = NULL );  //!< Load Z-matrix from XML node

	int GetNZ() const; //!< Get total number of atom centers used in Z-Matrix including dummy atoms 
	int GetNCrd() const; //!< Get total number of elemental coordinates in Z-matrix
	int GetNCrdUnFrozen() const; //!< Get the number of unfrozen elemental coordinates
	void FreezeCrdAll();   //!< Freeze all elemental coordinates 
	void UnFreezeCrdAll(); //!< Unfreeze all elemental coordinates
	int CalcBGMatr(); //!< Calculate Wilson B and G matricies for use in optimization calculations
	
	bool SetXYZCrd( int iat, double x, double y, double z ); //!< Set Cartesian coordinates for atom with index iat (0-based)

	void SetCrdDesc(int      iat, int iat_r,      int iat_ang,      int iat_dih); //!< Set Internal coordinate description for atom with index iat (0-based)
	void SetCrdDesc(HaAtom* aptr, HaAtom*  aptr_r, HaAtom* aptr_ang, HaAtom* aptr_dih); //!< Set Internal coordinate description for atom aptr
	void SetCrdDesc(const std::string& at_nm, const std::string& at_r_nm, const std::string& at_a_nm, 
		            const std::string& at_d_nm); //!< Set Internal coordinate description using atom references

	ElemCrd* GetRCrdForAtom  ( HaAtom* aptr ); //!< Get Distance elemental coordinate for atom  
	ElemCrd* GetAngCrdForAtom( HaAtom* aptr ); //!< Get Valence Angle elemental coordinate for atom 
	ElemCrd* GetDihCrdForAtom( HaAtom* aptr ); //!< Get Dihedral Angle elemental coordinate for atom

	ElemCrd* GetCrdByTag( const std::string& tag); //!< Get Elemental coordinate by tag
	ElemCrd* GetCrdByIdx( int idx ) { return elem_crds[idx]; } //!<  Get Elemental coordinate by Sequential Number (0-based index) 

	void SetTagRCrd  ( HaAtom* aptr, const std::string& tag ); //!< Set tag for distance coordinate for a given atom  
	void SetTagAngCrd( HaAtom* aptr, const std::string& tag ); //!< Set tag for angle    coordinate for a given atom  
	void SetTagDihCrd( HaAtom* aptr, const std::string& tag ); //!< Set tag for dihedral coordinate for a given atom
	
	void SetRVal(HaAtom* aptr, double rval, CoordUnits units = ANGSTROM_U); //!< Set distance coordinate for atom 
	void SetRVal(int iat,      double rval, CoordUnits units = ANGSTROM_U); //!< Set distance coordinate for atom idx (0-based)
	void SetRVal(const std::string& at_nm,  double rval, CoordUnits units = ANGSTROM_U); //!< Set distance coordinate for atom string reference (at_nm)

	void SetAngVal(HaAtom* aptr, double aval, AngleUnits units = DEGREE_U); //!< Set valence angle coordinate for atom 
	void SetAngVal(int iat,      double aval, AngleUnits units = DEGREE_U); //!< Set valence angle coordinate for atom idx (0-based)
	void SetAngVal(const std::string& at_nm, double aval, AngleUnits units = DEGREE_U); //!< Set valence angle coordinate for atom atom string reference (at_nm)

	void SetDihVal(HaAtom* aptr, double dval, AngleUnits units = DEGREE_U); //!< Set dihedral angle coordinate for atom 
	void SetDihVal(int iat,      double dval, AngleUnits units = DEGREE_U); //!< Set dihedral angle coordinate for atom idx (0-based)
	void SetDihVal(const std::string& at_nm, double dval, AngleUnits units = DEGREE_U); //!< Set dihedral angle coordinate for atom atom string reference (at_nm)

	HaAtom* AddDummyAtom(const std::string& at_name = "X"); //!< Add Dummy Atom to the Z-matrix 

protected:
	HaMolSet* pmset;

	std::vector< ElemCrd* >   elem_crds; 
	std::vector< SingleAtomCrdRule* > rules;

	std::vector<HaAtom*> dummy_atoms;

	static ZMatLoadOptions  opt_read_default; 
	static ZMatSaveOptions opt_save_default;

	HaMat_int    ib_mat;
	HaMat_double  b_mat;
	HaMat_double  g_mat; 
};



#endif /* !HAINTCRD_H */