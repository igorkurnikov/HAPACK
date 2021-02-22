/*!  \file mm_params.h

    Classes to define different parameters of Molecular Mechanics Model and computations

    \author Igor Kurnikov 
    \date 2010-

*/
#ifndef MM_PARAMS_H
#define MM_PARAMS_H

#include "hastring.h"
#include "hatypes.h"

//template<typename T> 
//class HaEnum2 : public HaEnum1
//{
//public:
//	HaEnum2() {}
//	virtual ~HaEnum2();
//	
//	operator int() const { return v_; }
//	bool operator==(const T& val) const { return v_ == val; }
//	bool operator!=(const T& val) const { return v_ != val; } 
//
//	virtual IntStrMap& GetLabelsMap() const { return labels; }
//	virtual int& value() { return (int &) v_; }
//	virtual const char* label() const { return labels[v_].c_str();}
//	virtual int SetWithValue(int val) { v_ = (T) val; return TRUE; }
//	virtual int SetWithLabel(const char* label)
//	{
//		IntStrMap::iterator itr;
//		for(itr = labels.begin(); itr != labels.end(); itr++)
//		{
//			if( (*itr).second == label_set ) 
//			{
//				v_ = (T) (*itr).first;
//			}
//			return TRUE;
//		}
//		return FALSE;
//	}
//
//	HaEnum2(T v): v_(v) {}
//private:
//	T v_;
//	static IntStrMap labels;
//};

class HaMolMechMod;
class MolMechModel;

class MMRunType : public HaEnum1
//! Enum class for Molecular Mechanics simulation type 
{
public:
	MMRunType();
	virtual ~MMRunType();

	enum Value { MD_RUN = 0, MIN_RUN = 1, ENER_RUN = 2};
    
	MMRunType& operator=( int value) { SetWithValue(value); return (*this); }
	operator int() const { return v_; }
	bool operator==(const Value& val) const { return v_ == val; }
	bool operator!=(const Value& val) const { return v_ != val; } 

	virtual IntStrMap& GetLabelsMap() const { return labels; }
	virtual int& value() { return (int &) v_; }
	virtual const char* label() const { return labels[v_].c_str();}
	virtual int SetWithValue(int val);
	virtual int SetWithLabel(const char* label);

	MMRunType(Value v): v_(v) {}
private:
	Value v_;
	static IntStrMap labels;
};

namespace swig {
	const MMRunType MD_RUN   = MMRunType::MD_RUN;
	const MMRunType MIN_RUN  = MMRunType::MIN_RUN;
	const MMRunType ENER_RUN = MMRunType::ENER_RUN;
}

class TempCtrlMethod : public HaEnum1
//! Extended Enum class for Molecular Dynamics Temperature Control Method 
{
public:
	TempCtrlMethod();
	virtual ~TempCtrlMethod();

	enum Value { CONST_ENE_MD = 0, CONST_TEMP_BERENDSEN = 1, CONST_TEMP_RANDOMIZED = 2, 
		         CONST_TEMP_LANGEVIN  = 3 };

	operator int() const { return v_; }
	bool operator==(const Value& val) const { return v_ == val; }
	bool operator!=(const Value& val) const { return v_ != val; } 

	virtual IntStrMap& GetLabelsMap() const { return labels; }
	virtual int& value() { return (int &) v_; }
	virtual const char* label() const { return labels[v_].c_str();}
	virtual int SetWithValue(int val);
	virtual int SetWithLabel(const char* label);

	TempCtrlMethod(Value v): v_(v) {}
private:
	Value v_;
	static IntStrMap labels;
	
};

namespace swig {
	const TempCtrlMethod CONST_ENE_MD          = TempCtrlMethod::CONST_ENE_MD;
	const TempCtrlMethod CONST_TEMP_BERENDSEN  = TempCtrlMethod::CONST_TEMP_BERENDSEN;
	const TempCtrlMethod CONST_TEMP_RANDOMIZED = TempCtrlMethod::CONST_TEMP_RANDOMIZED;
	const TempCtrlMethod CONST_TEMP_LANGEVIN   = TempCtrlMethod::CONST_TEMP_LANGEVIN;
}

class EneMinMethod : public HaEnum1
//! Enum class for Energy Minimization Method 
{
public:
	EneMinMethod();
	virtual ~EneMinMethod();

	enum Value { CONJ_GRAD = 0, SD_AND_CG = 1, STEEPEST_DESCENT = 2 };

	operator int() const { return v_; }
	bool operator==(const Value& val) const { return v_ == val; }
	bool operator!=(const Value& val) const { return v_ != val; } 

	virtual IntStrMap& GetLabelsMap() const { return labels; }
	virtual int& value() { return (int &) v_; }
	virtual const char* label() const { return labels[v_].c_str();}
	virtual int SetWithValue(int val);
	virtual int SetWithLabel(const char* label);

	EneMinMethod(Value v): v_(v) {}
private:
	Value v_;
	static IntStrMap labels;
};

namespace swig {
	const EneMinMethod CONJ_GRAD = EneMinMethod::CONJ_GRAD;
	const EneMinMethod SD_AND_CG = EneMinMethod::SD_AND_CG;
	const EneMinMethod STEEPEST_DESCENT = EneMinMethod::STEEPEST_DESCENT;
}

class MMReadInitCrdType : public HaEnum1
//! Enum class for type of inital coordinates read  
{
public:
	MMReadInitCrdType();
	virtual ~MMReadInitCrdType();

	enum Value { READ_X_FORM = 1, READ_X_BIN =2, READ_XV_BIN = 4, 
		         READ_XV_FORM = 5, READ_XVBOX_BIN=6, READ_XVBOX_FORM = 7 };

	operator int() const { return v_; }
	bool operator==(const Value& val) const { return v_ == val; }
	bool operator!=(const Value& val) const { return v_ != val; } 

	virtual IntStrMap& GetLabelsMap() const { return labels; }
	virtual int& value() { return (int &) v_; }
	virtual const char* label() const { return labels[v_].c_str();}
	virtual int SetWithValue(int val);
	virtual int SetWithLabel(const char* label);

	MMReadInitCrdType(Value v): v_(v) {}
private:
	Value v_;
	static IntStrMap labels;
};

namespace swig {
	const MMReadInitCrdType READ_X_FORM  = MMReadInitCrdType::READ_X_FORM;
	const MMReadInitCrdType READ_X_BIN   = MMReadInitCrdType::READ_X_BIN;
	const MMReadInitCrdType READ_XV_BIN  = MMReadInitCrdType::READ_XV_BIN;
	const MMReadInitCrdType READ_XV_FORM    = MMReadInitCrdType::READ_XV_FORM;
	const MMReadInitCrdType READ_XVBOX_BIN  = MMReadInitCrdType::READ_XVBOX_BIN;
	const MMReadInitCrdType READ_XVBOX_FORM = MMReadInitCrdType::READ_XVBOX_FORM;
}

class CrdFormatParam : public HaEnum1
//! Enum class for format of saved coordinates   
{
public:
	CrdFormatParam();
	virtual ~CrdFormatParam();

	enum Value { FORMATTED = 0, BINARY = 1 };

	operator int() const { return v_; }
	bool operator==(const Value& val) const { return v_ == val; }
	bool operator!=(const Value& val) const { return v_ != val; } 

	virtual IntStrMap& GetLabelsMap() const { return labels; }
	virtual int& value() { return (int &) v_; }
	virtual const char* label() const { return labels[v_].c_str();}
	virtual int SetWithValue(int val);
	virtual int SetWithLabel(const char* label);

	CrdFormatParam(Value v): v_(v) {}
private:
	Value v_;
	static IntStrMap labels;
};


class PerBoundaryCondType : public HaEnum1
//! Enum class for Periodical Boundary conditions type 
{
public:
	PerBoundaryCondType();
	virtual ~PerBoundaryCondType();

	void SetMolMechMod(HaMolMechMod* p_mm_mod ); //!< Set Molecular Mechanics module associated with the variable

	enum Value { NO_PERIODICITY = 0, CONST_VOL = 1, CONST_PRES = 2};

	operator int() const { return v_; }
	bool operator==(const Value& val) const { return v_ == val; }
	bool operator!=(const Value& val) const { return v_ != val; } 

	virtual IntStrMap& GetLabelsMap() const { return labels; }
	virtual int& value() { return (int &) v_; }
	virtual const char* label() const { return labels[v_].c_str();}
	virtual int SetWithValue(int val);
	virtual int SetWithLabel(const char* label);

	virtual std::vector<std::string> GetActiveLabels(); //!< Get labels compatible with current state of molecular mechanics module
	virtual void SetCompatValue(); //!< Set value of the variable compatible with the current context (state of Molecular Mechanics Module)

	PerBoundaryCondType(Value v): v_(v) {}

protected:
	Value v_;
	static IntStrMap labels;
	HaMolMechMod* p_mm_mod; 
};

namespace swig {
	const PerBoundaryCondType NO_PERIODICITY  = PerBoundaryCondType::NO_PERIODICITY;
	const PerBoundaryCondType CONST_VOL       = PerBoundaryCondType::CONST_VOL;
	const PerBoundaryCondType CONST_PRES      = PerBoundaryCondType::CONST_PRES;
}


class PressureRegMethod : public HaEnum1
//! Enum class for Pressure Regulation method in Molecular Mechanics simulation 
{
public:
	enum Value { NO_CRD_SCALING=0, ISOTROP_CRD_SCALING=1, ANISOTROP_CRD_SCALING=2, 
				 CRD_SCALING_XY_AND_Z = 3, CRD_SCALING_ONLY_Z = 4, CRD_SCALING_XZ_AND_Y = 5, 
				 CRD_SCALING_YZ_AND_X = 6 };

	PressureRegMethod();
	PressureRegMethod(Value v): v_(v) {}
	virtual ~PressureRegMethod();

	operator int() const { return v_; }
	bool operator==(const Value& val) const { return v_ == val; }
	bool operator!=(const Value& val) const { return v_ != val; } 

	virtual IntStrMap& GetLabelsMap() const { return labels; }
	virtual int& value() { return (int &) v_; }
	virtual const char* label() const { return labels[v_].c_str();}
	virtual int SetWithValue(int val);
	virtual int SetWithLabel(const char* label);

private:
	Value v_;
	static IntStrMap labels;
};

namespace swig {
	const PressureRegMethod NO_CRD_SCALING         = PressureRegMethod::NO_CRD_SCALING;
	const PressureRegMethod ISOTROP_CRD_SCALING    = PressureRegMethod::ISOTROP_CRD_SCALING;
	const PressureRegMethod ANISOTROP_CRD_SCALING  = PressureRegMethod::ANISOTROP_CRD_SCALING;
	const PressureRegMethod CRD_SCALING_XY_AND_Z   = PressureRegMethod::CRD_SCALING_XY_AND_Z;
	const PressureRegMethod CRD_SCALING_ONLY_Z     = PressureRegMethod::CRD_SCALING_ONLY_Z;
	const PressureRegMethod CRD_SCALING_XZ_AND_Y   = PressureRegMethod::CRD_SCALING_XZ_AND_Y;
	const PressureRegMethod CRD_SCALING_YZ_AND_X   = PressureRegMethod::CRD_SCALING_YZ_AND_X;
}

class MMElectrMethod : public HaEnum1
//! Enum class for method of modeling of electrostatic interactions in MM calculations
{
public:
	enum Value { DIST_DEP_DIEL = 0, PME_METHOD = 1, GEN_BORN = 2, SCREENED_COULOMB = 3};

	MMElectrMethod();
	MMElectrMethod(Value v): v_(v) {}
	virtual ~MMElectrMethod();

	void SetMMModel(MolMechModel* p_mm_model ); //!< Set Molecular Mechanics model associated with the variable
	
	operator int() const { return v_; }
	bool operator==(const Value& val) const { return v_ == val; }
	bool operator!=(const Value& val) const { return v_ != val; } 

	virtual IntStrMap& GetLabelsMap() const { return labels; }
	virtual int& value() { return (int &) v_; }
	virtual const char* label() const { return labels[v_].c_str();}
	virtual int SetWithValue(int val);
	virtual int SetWithLabel(const char* label);

	virtual std::vector<std::string> GetActiveLabels(); //!< Get labels compatible with current state of molecular mechanics module
	virtual void SetCompatValue(); //!< Set value of the variable compatible with the current context (state of Molecular Mechanics Module)

protected:
	Value v_;
	static IntStrMap labels;
	MolMechModel* p_mm_model;
};

namespace swig {
	const MMElectrMethod DIST_DEP_DIEL    = MMElectrMethod::DIST_DEP_DIEL;
	const MMElectrMethod PME_METHOD       = MMElectrMethod::PME_METHOD;
	const MMElectrMethod GEN_BORN         = MMElectrMethod::GEN_BORN;
	const MMElectrMethod SCREENED_COULOMB = MMElectrMethod::SCREENED_COULOMB;
}


class StartVelMethod : public HaEnum1
//! Enum class for Method to generate starting velocities 
{
public:
	StartVelMethod();
	virtual ~StartVelMethod();

	enum Value {MAXWELL_START_VEL = 3, READ_START_VEL = 4};

	operator int() const { return v_; }
	bool operator==(const Value& val) const { return v_ == val; }
	bool operator!=(const Value& val) const { return v_ != val; } 

	virtual IntStrMap& GetLabelsMap() const { return labels; }
	virtual int& value() { return (int &) v_; }
	virtual const char* label() const { return labels[v_].c_str();}
	virtual int SetWithValue(int val);
	virtual int SetWithLabel(const char* label);

	StartVelMethod(Value v): v_(v) {}
private:
	Value v_;
	static IntStrMap labels;
};

namespace swig {
	const StartVelMethod MAXWELL_START_VEL = StartVelMethod::MAXWELL_START_VEL;
	const StartVelMethod READ_START_VEL    = StartVelMethod::READ_START_VEL;
}


class OmitInteractionsParam : public HaEnum1
//! Enum class for a parameter to omit certain interaction terms in MM Model 
{
public:
	OmitInteractionsParam();
	virtual ~OmitInteractionsParam();

	enum Value { CALC_ALL_INTER = 1, OMIT_BONDS_H = 2, OMIT_BONDS = 3 ,
                 OMIT_BONDS_VANG_H = 4, OMIT_BONDS_VANG = 5,           
		         OMIT_BONDS_VANG_DIH_H = 6, OMIT_BONDS_VANG_DIH = 7 };

	operator int() const { return v_; }
	bool operator==(const Value& val) const { return v_ == val; }
	bool operator!=(const Value& val) const { return v_ != val; } 

	virtual IntStrMap& GetLabelsMap() const { return labels; }
	virtual int& value() { return (int &) v_; }
	virtual const char* label() const { return labels[v_].c_str();}
	virtual int SetWithValue(int val);
	virtual int SetWithLabel(const char* label);

	OmitInteractionsParam(Value v): v_(v) {}
private:
	Value v_;
	static IntStrMap labels;
};

namespace swig {
	const OmitInteractionsParam CALC_ALL_INTER = OmitInteractionsParam::CALC_ALL_INTER;
	const OmitInteractionsParam OMIT_BONDS_H   = OmitInteractionsParam::OMIT_BONDS_H;
	const OmitInteractionsParam OMIT_BONDS     = OmitInteractionsParam::OMIT_BONDS;
	const OmitInteractionsParam OMIT_BONDS_VANG_H      = OmitInteractionsParam::OMIT_BONDS_VANG_H;
	const OmitInteractionsParam OMIT_BONDS_VANG_DIH_H  = OmitInteractionsParam::OMIT_BONDS_VANG_DIH_H;
	const OmitInteractionsParam OMIT_BONDS_VANG_DIH    = OmitInteractionsParam::OMIT_BONDS_VANG_DIH;
}

class MMShakeParam : public HaEnum1
//! Enum class for Molecular Mechanics Shake Algorithm to freeze bonds 
{
public:
	MMShakeParam();
	virtual ~MMShakeParam();

	enum Value {NO_SHAKE = 1, H_ATOM_SHAKE=2, ALL_BOND_SHAKE=3};

	operator int() const { return v_; }
	bool operator==(const Value& val) const { return v_ == val; }
	bool operator!=(const Value& val) const { return v_ != val; } 

	virtual IntStrMap& GetLabelsMap() const { return labels; }
	virtual int& value() { return (int &) v_; }
	virtual const char* label() const { return labels[v_].c_str();}
	virtual int SetWithValue(int val);
	virtual int SetWithLabel(const char* label);

	MMShakeParam(Value v): v_(v) {}
private:
	Value v_;
	static IntStrMap labels;
};

namespace swig {
	const MMShakeParam NO_SHAKE       = MMShakeParam::NO_SHAKE;
	const MMShakeParam H_ATOM_SHAKE   = MMShakeParam::H_ATOM_SHAKE;
	const MMShakeParam ALL_BOND_SHAKE = MMShakeParam::ALL_BOND_SHAKE;
}

class MMExternalProg : public HaEnum1
//! Enum class for parameter defining external Molecular Mechanics simulation program
{
public:
	MMExternalProg();
	virtual ~MMExternalProg();

	enum Value { /* PMEMD_9 = 0, SANDER_9 = 1, PMEMD_10 = 2, SANDER_10 = 3, 
		         PMEMD_12 = 4, */ SANDER_12 = 5, PMEMD_12 = 6, PMEMD_18 = 7, TINKER_51 = 8, GROMACS_51 = 9 };

	operator int() const { return v_; }
	bool operator==(const Value& val) const { return v_ == val; }
	bool operator!=(const Value& val) const { return v_ != val; } 

	virtual IntStrMap& GetLabelsMap() const { return labels; }
	virtual int& value() { return (int &) v_; }
	virtual const char* label() const { return labels[v_].c_str();}
	virtual int SetWithValue(int val);
	virtual int SetWithLabel(const char* label);

	MMExternalProg(Value v): v_(v) {}
private:
	Value v_;
	static IntStrMap labels;
};

namespace swig {
	/* const MMExternalProg PMEMD_9    = MMExternalProg::PMEMD_9;
	const MMExternalProg SANDER_9   = MMExternalProg::SANDER_9;
	const MMExternalProg PMEMD_10   = MMExternalProg::PMEMD_10;
	const MMExternalProg SANDER_10  = MMExternalProg::SANDER_10; */
	const MMExternalProg PMEMD_12   = MMExternalProg::PMEMD_12;
	const MMExternalProg SANDER_12  = MMExternalProg::SANDER_12;
	const MMExternalProg PMEMD_18   = MMExternalProg::PMEMD_18;
	const MMExternalProg TINKER_51  = MMExternalProg::TINKER_51;
	const MMExternalProg GROMACS_51 = MMExternalProg::GROMACS_51;
}

class ForceFieldType : public HaEnum1
//! Enum class for Molecular Mechanics Force Field Type 
{
public:
	ForceFieldType();
	virtual ~ForceFieldType();

	enum Value { AMBER_94 = 0, AMBER_99_SB = 1, AMBER_99_BSC0 = 2, AMBER_03 = 3, AMBER_10 = 4, 
		         AMOEBA = 5, ARROW_5_14_CT=6, ARROW_2_0 = 7, UNKNOWN_FF = 8};
    
	ForceFieldType& operator=( int value) { SetWithValue(value); return (*this); }
	operator int() const { return v_; }
	bool operator==(const Value& val) const { return v_ == val; }
	bool operator!=(const Value& val) const { return v_ != val; } 

	virtual IntStrMap& GetLabelsMap() const { return labels; }
	virtual int& value() { return (int &) v_; }
	virtual const char* label() const { return labels[v_].c_str();}
	virtual int SetWithValue(int val);
	virtual int SetWithLabel(const char* label);

	ForceFieldType(Value v): v_(v) {}
private:
	Value v_;
	static IntStrMap labels;
};

namespace swig {
	const ForceFieldType AMBER_94   = ForceFieldType::AMBER_94;
	const ForceFieldType AMBER_99_SB   = ForceFieldType::AMBER_99_SB;
	const ForceFieldType AMBER_99_BSC0 = ForceFieldType::AMBER_99_BSC0;
	const ForceFieldType AMBER_03   = ForceFieldType::AMBER_03;
	const ForceFieldType AMBER_10   = ForceFieldType::AMBER_10;
	const ForceFieldType AMOEBA     = ForceFieldType::AMOEBA;
	const ForceFieldType ARROW_5_14_CT = ForceFieldType::ARROW_5_14_CT;
	const ForceFieldType ARROW_2_0     = ForceFieldType::ARROW_2_0;
	const ForceFieldType UNKNOWN_FF = ForceFieldType::UNKNOWN_FF;
}

class AtomContactType : public HaEnum1
//! Enum class for Type of Atom Atom Distance Constraint Type 
{
public:
	AtomContactType();
	virtual ~AtomContactType();

	enum Value { HARMONIC_CNT = 0, VDW_CNT_6_12 = 1, VDW_CNT_10_12 = 2, VDW_CNT_6_12_NO_REP = 3 };
    
	AtomContactType& operator=( int value) { SetWithValue(value); return (*this); }
	operator int() const { return v_; }
	bool operator==(const Value& val) const { return v_ == val; }
	bool operator!=(const Value& val) const { return v_ != val; } 

	virtual IntStrMap& GetLabelsMap() const { return labels; }
	virtual int& value() { return (int &) v_; }
	virtual const char* label() const { return labels[v_].c_str();}
	virtual int SetWithValue(int val);
	virtual int SetWithLabel(const char* label);

	AtomContactType(Value v): v_(v) {}
private:
	Value v_;
	static IntStrMap labels;
};

namespace swig {
	const AtomContactType VDW_CNT_6_12          = AtomContactType::VDW_CNT_6_12;
	const AtomContactType VDW_CNT_10_12  = AtomContactType::VDW_CNT_10_12;
	const AtomContactType HARMONIC_CNT          = AtomContactType::HARMONIC_CNT;
}


#endif // End #define MM_PARAMS_H
