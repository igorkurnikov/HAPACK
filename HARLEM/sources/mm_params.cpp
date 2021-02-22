/*!  \file mm_params.cpp

    Classes to define different parameters of Molecular Mechanics Model and computations

    \author Igor Kurnikov 
    \date 2010-

*/

#include <mpi.h>

#include "haconst.h"
#include "hastl.h"
#include "haio.h"
#include "hamolset.h"
#include "mm_params.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "hamolmech.h"

#define MM_PARAMS_CPP

IntStrMap SetMMRunTypeLbls()
{
	IntStrMap lmap;
	lmap[MMRunType::MD_RUN]   = "Molecular Dynamics";
	lmap[MMRunType::MIN_RUN]  = "Energy Minimization";
	lmap[MMRunType::ENER_RUN] = "Single Energy";
	return lmap;
}

IntStrMap MMRunType::labels = SetMMRunTypeLbls();  

MMRunType::MMRunType()
{
	v_= MD_RUN;	
}
	
MMRunType::~MMRunType()
{

}

int MMRunType::SetWithValue(int val)
{
	v_ = (Value) val;
	return TRUE;
}

int MMRunType::SetWithLabel(const char* label_set)
{
	IntStrMap::iterator itr;
	for(itr = labels.begin(); itr != labels.end(); itr++)
	{
		if( (*itr).second == label_set ) 
		{
			v_ = (Value) (*itr).first;
			return TRUE;
		}
	}
	return FALSE;
}

IntStrMap SetTempCtrlMethodLbls()
{
	IntStrMap lmap;
	lmap[TempCtrlMethod::CONST_ENE_MD]          = "Const Energy";
	lmap[TempCtrlMethod::CONST_TEMP_BERENDSEN]    = "Const Temperature Berendsen";
	lmap[TempCtrlMethod::CONST_TEMP_RANDOMIZED] = "Const Temperature Randomized";
	lmap[TempCtrlMethod::CONST_TEMP_LANGEVIN]   = "Langevin Dynamics";
	return lmap;
}

IntStrMap TempCtrlMethod::labels = SetTempCtrlMethodLbls();  

TempCtrlMethod::TempCtrlMethod()
{
	v_= CONST_TEMP_LANGEVIN;	
}
	
TempCtrlMethod::~TempCtrlMethod()
{

}

int TempCtrlMethod::SetWithValue(int val)
{
	v_ = (Value) val;
	return TRUE;
}


int TempCtrlMethod::SetWithLabel(const char* label_set)
{
	IntStrMap::iterator itr;
	for(itr = labels.begin(); itr != labels.end(); itr++)
	{
		if( (*itr).second == label_set ) 
		{
			v_ = (Value) (*itr).first;
			return TRUE;
		}
	}
	return FALSE;
}


IntStrMap SetEneMinMethodLbls()
{
	IntStrMap lmap;
	
	lmap[EneMinMethod::CONJ_GRAD] =  "Conjugate Gradients";
	lmap[EneMinMethod::SD_AND_CG] =  "Steepest Decent and Conjugate Gradients";
	lmap[EneMinMethod::STEEPEST_DESCENT] = "Steepest Decent";
	return lmap;
}

IntStrMap EneMinMethod::labels = SetEneMinMethodLbls();  

EneMinMethod::EneMinMethod()
{
	v_= STEEPEST_DESCENT;	
}
	
EneMinMethod::~EneMinMethod()
{

}

int EneMinMethod::SetWithValue(int val)
{
	v_ = (Value) val;
	return TRUE;
}

int EneMinMethod::SetWithLabel(const char* label_set)
{
	IntStrMap::iterator itr;
	for(itr = labels.begin(); itr != labels.end(); itr++)
	{
		if( (*itr).second == label_set ) 
		{
			v_ = (Value) (*itr).first;
			return TRUE;
		}
	}
	return FALSE;
}

IntStrMap SetMMReadInitCrdTypeLbls()
{
	IntStrMap lmap;
	lmap[MMReadInitCrdType::READ_X_FORM]   = "Read Coords Formatted";
	lmap[MMReadInitCrdType::READ_X_BIN]   = "Read Coords Binary";
	lmap[MMReadInitCrdType::READ_XV_BIN]      = "Read Coords and Vels Binary";
	lmap[MMReadInitCrdType::READ_XV_FORM]     = "Read Coords and Vels Formatted";
	lmap[MMReadInitCrdType::READ_XVBOX_BIN]   = "Read Coords, Vels and Box Binary";
	lmap[MMReadInitCrdType::READ_XVBOX_FORM]   = "Read Coords, Vels and Box Formatted";

	return lmap;
}

IntStrMap MMReadInitCrdType::labels = SetMMReadInitCrdTypeLbls();  

MMReadInitCrdType::MMReadInitCrdType()
{
	v_= READ_X_FORM;	
}
	
MMReadInitCrdType::~MMReadInitCrdType()
{

}

int MMReadInitCrdType::SetWithValue(int val)
{
	v_ = (Value) val;
	return TRUE;
}

int MMReadInitCrdType::SetWithLabel(const char* label_set)
{
	IntStrMap::iterator itr;
	for(itr = labels.begin(); itr != labels.end(); itr++)
	{
		if( (*itr).second == label_set ) 
		{
			v_ = (Value) (*itr).first;
			return TRUE;
		}
	}
	return FALSE;
}

IntStrMap SetCrdFormatParamLbls()
{
	IntStrMap lmap;
	lmap[CrdFormatParam::BINARY]     = "Binary";
	lmap[CrdFormatParam::FORMATTED]  = "Formatted";
	return lmap;
}

IntStrMap CrdFormatParam::labels = SetCrdFormatParamLbls();  

CrdFormatParam::CrdFormatParam()
{
	v_= FORMATTED;	
}
	
CrdFormatParam::~CrdFormatParam()
{

}

int CrdFormatParam::SetWithValue(int val)
{
	v_ = (Value) val;
	return TRUE;
}

int CrdFormatParam::SetWithLabel(const char* label_set)
{
	IntStrMap::iterator itr;
	for(itr = labels.begin(); itr != labels.end(); itr++)
	{
		if( (*itr).second == label_set ) 
		{
			v_ = (Value) (*itr).first;
			return TRUE;
		}
	}
	return FALSE;
}

IntStrMap SetPerBoundaryCondTypeLbls()
{
	IntStrMap lmap;
	lmap[PerBoundaryCondType::NO_PERIODICITY]   = "No Periodic Box";
	lmap[PerBoundaryCondType::CONST_VOL]        = "Const Volume";
	lmap[PerBoundaryCondType::CONST_PRES]       = "Const Pressure";
	return lmap;
}

IntStrMap PerBoundaryCondType::labels = SetPerBoundaryCondTypeLbls();  

PerBoundaryCondType::PerBoundaryCondType()
{
	v_= NO_PERIODICITY;	
	p_mm_mod = NULL;
}
	
PerBoundaryCondType::~PerBoundaryCondType()
{

}

void PerBoundaryCondType::SetMolMechMod(HaMolMechMod* p_mm_mod_par )
{
	p_mm_mod = p_mm_mod_par;
}

int PerBoundaryCondType::SetWithValue(int val)
{
	v_ = (Value) val;
	return TRUE;
}

int PerBoundaryCondType::SetWithLabel(const char* label_set)
{
	IntStrMap::iterator itr;
	for(itr = labels.begin(); itr != labels.end(); itr++)
	{
		if( (*itr).second == label_set ) 
		{
			v_ = (Value) (*itr).first;
			return TRUE;
		}
	}
	return FALSE;
}

std::vector<std::string> PerBoundaryCondType::GetActiveLabels()
{
	if( p_mm_mod == NULL ) return GetAllLabels();
	MolSet* pmset = p_mm_mod->GetMolSet();
	IntStrMap& lbl_map = GetLabelsMap();
	std::vector<std::string> labels;
	IntStrMap::iterator itr;
	for(itr = lbl_map.begin(); itr != lbl_map.end(); itr++)
	{
		int ival = (*itr).first;
		if( pmset->per_bc->IsSet() )
		{
			if( ival == NO_PERIODICITY ) continue;
		}
		else
		{
			if( ival == CONST_VOL ) continue;
			if( ival == CONST_PRES ) continue;
		}
		labels.push_back((*itr).second);
	}
	SetCompatValue();
	return labels;
}

void PerBoundaryCondType::SetCompatValue()
{
	if( p_mm_mod == NULL ) return;
	MolSet* pmset = p_mm_mod->GetMolSet();
	if( pmset->per_bc->IsSet() )
	{
		if( v_ == NO_PERIODICITY ) 
		{
			if( p_mm_mod->run_type == p_mm_mod->run_type.MD_RUN )
			{
				v_ = CONST_PRES;
			}
			else
			{
				v_ = CONST_VOL;
			}
		}
		if( v_ == CONST_PRES )
		{
			if( p_mm_mod->run_type == p_mm_mod->run_type.MIN_RUN || 
				p_mm_mod->run_type == p_mm_mod->run_type.ENER_RUN )
			{
				v_ = CONST_VOL;
			}
		}
	}
	else
	{
		if( v_ == CONST_VOL ) v_ = NO_PERIODICITY;
		if( v_ == CONST_PRES ) v_ = NO_PERIODICITY;
	}
}


IntStrMap SetPressureRegMethodLbls()
{
	IntStrMap lmap;
	lmap[PressureRegMethod::NO_CRD_SCALING]        = "No Coord Scaling";
	lmap[PressureRegMethod::ISOTROP_CRD_SCALING]   = "Isotropic Coord Scaling";
	lmap[PressureRegMethod::ANISOTROP_CRD_SCALING] = "Anisotropic Coord Scaling";
	lmap[PressureRegMethod::CRD_SCALING_XY_AND_Z] = "Coord Scaling XY and Z";
	lmap[PressureRegMethod::CRD_SCALING_ONLY_Z]   = "Coord Scaling Only Z";
	lmap[PressureRegMethod::CRD_SCALING_XZ_AND_Y] = "Coord Scaling XZ and Y";
	lmap[PressureRegMethod::CRD_SCALING_YZ_AND_X] = "Coord Scaling YZ and X";
	return lmap;
}

IntStrMap PressureRegMethod::labels = SetPressureRegMethodLbls();  

PressureRegMethod::PressureRegMethod()
{
	v_= ISOTROP_CRD_SCALING;	
}
	
PressureRegMethod::~PressureRegMethod()
{

}

int PressureRegMethod::SetWithValue(int val)
{
	v_ = (Value) val;
	return TRUE;
}

int PressureRegMethod::SetWithLabel(const char* label_set)
{
	IntStrMap::iterator itr;
	for(itr = labels.begin(); itr != labels.end(); itr++)
	{
		if( (*itr).second == label_set ) 
		{
			v_ = (Value) (*itr).first;
			return TRUE;
		}
	}
	return FALSE;
}

IntStrMap SetMMElectrMethodLbls()
{
	IntStrMap lmap;
	lmap[MMElectrMethod::DIST_DEP_DIEL]    = "Distant Dependent Dielectric";
	lmap[MMElectrMethod::PME_METHOD]       = "Particle Mesh Ewald(PME)";
	lmap[MMElectrMethod::GEN_BORN]         = "Generalized Born";
	lmap[MMElectrMethod::SCREENED_COULOMB] = "Screened Coulomb";
	return lmap;
}

IntStrMap MMElectrMethod::labels = SetMMElectrMethodLbls();  

MMElectrMethod::MMElectrMethod()
{
	v_= PME_METHOD;	
	p_mm_model = NULL;
}
	
MMElectrMethod::~MMElectrMethod()
{

}

void MMElectrMethod::SetMMModel(MolMechModel* p_mm_model_par )
{
	p_mm_model = p_mm_model_par;
}

int MMElectrMethod::SetWithValue(int val)
{
	v_ = (Value) val;
	return TRUE;
}

int MMElectrMethod::SetWithLabel(const char* label_set)
{
	IntStrMap::iterator itr;
	for(itr = labels.begin(); itr != labels.end(); itr++)
	{
		if( (*itr).second == label_set ) 
		{
			v_ = (Value) (*itr).first;
			return TRUE;
		}
	}
	return FALSE;
}

std::vector<std::string> MMElectrMethod::GetActiveLabels()
{
	if( p_mm_model == NULL ) return GetAllLabels();
	MolSet* pmset = p_mm_model->GetMolSet();
	IntStrMap& lbl_map = GetLabelsMap();
	std::vector<std::string> labels;
	IntStrMap::iterator itr;
	for(itr = lbl_map.begin(); itr != lbl_map.end(); itr++)
	{
		int ival = (*itr).first;
		if( pmset->per_bc->IsSet() )
		{
			if( ival == GEN_BORN ) continue;
		}
		else
		{
			if( ival == PME_METHOD ) continue;
		}
		labels.push_back((*itr).second);
	}
	SetCompatValue();
	return labels;
}

void MMElectrMethod::SetCompatValue()
{
	if( p_mm_model == NULL ) return;
	MolSet* pmset = p_mm_model->GetMolSet();
	if( pmset->per_bc->IsSet() )
	{
		if( v_ == GEN_BORN ) v_ = PME_METHOD;
	}
	else
	{
		if( v_ == PME_METHOD ) v_ = GEN_BORN;
	}
}

IntStrMap SetStartVelMethodLbls()
{
	IntStrMap lmap;
	lmap[StartVelMethod::MAXWELL_START_VEL] = "Maxwell Start Velocities";
	lmap[StartVelMethod::READ_START_VEL]    = "Read Start Velocities";
	return lmap;
}

IntStrMap StartVelMethod::labels = SetStartVelMethodLbls();  

StartVelMethod::StartVelMethod()
{
	v_= MAXWELL_START_VEL;	
}
	
StartVelMethod::~StartVelMethod()
{

}

int StartVelMethod::SetWithValue(int val)
{
	v_ = (Value) val;
	return TRUE;
}

int StartVelMethod::SetWithLabel(const char* label_set)
{
	IntStrMap::iterator itr;
	for(itr = labels.begin(); itr != labels.end(); itr++)
	{
		if( (*itr).second == label_set ) 
		{
			v_ = (Value) (*itr).first;
			return TRUE;
		}
	}
	return FALSE;
}

IntStrMap SetOmitInteractionsParamLbls()
{
	IntStrMap lmap;
	lmap[OmitInteractionsParam::CALC_ALL_INTER]   = "Calculate All interactions";
	lmap[OmitInteractionsParam::OMIT_BONDS_H]    = "Omit Bonds with H";
	lmap[OmitInteractionsParam::OMIT_BONDS] = "Omit Bonds";
	lmap[OmitInteractionsParam::OMIT_BONDS_VANG_H] = "Omit Bonds and Val Angles with H";
	lmap[OmitInteractionsParam::OMIT_BONDS_VANG] = "Omit Bonds and Val Angles";
	lmap[OmitInteractionsParam::OMIT_BONDS_VANG_DIH_H] = "Omit Bonds, Angles and Dihedrals with H";
	lmap[OmitInteractionsParam::OMIT_BONDS_VANG_DIH] = "Omit Bonds, Angles and Dihedrals";
	return lmap;
}

IntStrMap OmitInteractionsParam::labels = SetOmitInteractionsParamLbls();  

OmitInteractionsParam::OmitInteractionsParam()
{
	v_= CALC_ALL_INTER;	
}
	
OmitInteractionsParam::~OmitInteractionsParam()
{

}

int OmitInteractionsParam::SetWithValue(int val)
{
	v_ = (Value) val;
	return TRUE;
}

int OmitInteractionsParam::SetWithLabel(const char* label_set)
{
	IntStrMap::iterator itr;
	for(itr = labels.begin(); itr != labels.end(); itr++)
	{
		if( (*itr).second == label_set ) 
		{
			v_ = (Value) (*itr).first;
			return TRUE;
		}
	}
	return FALSE;
}

IntStrMap SetMMShakeParamLbls()
{
	IntStrMap lmap;
	lmap[MMShakeParam::NO_SHAKE]       = "No SHAKE";
	lmap[MMShakeParam::H_ATOM_SHAKE]   = "SHAKE for bonds with H";
	lmap[MMShakeParam::ALL_BOND_SHAKE] = "SHAKE for all bonds";
	return lmap;
}

IntStrMap MMShakeParam::labels = SetMMShakeParamLbls();  

MMShakeParam::MMShakeParam()
{
	v_= NO_SHAKE;	
}
	
MMShakeParam::~MMShakeParam()
{

}

int MMShakeParam::SetWithValue(int val)
{
	v_ = (Value) val;
	return TRUE;
}

int MMShakeParam::SetWithLabel(const char* label_set)
{
	IntStrMap::iterator itr;
	for(itr = labels.begin(); itr != labels.end(); itr++)
	{
		if( (*itr).second == label_set ) 
		{
			v_ = (Value) (*itr).first;
			return TRUE;
		}
	}
	return FALSE;
}

IntStrMap SetMMExternalProgLbls()
{
	IntStrMap lmap;
//	lmap[MMExternalProg::PMEMD_9  ]  = "AMBER 9: PMEMD";
//	lmap[MMExternalProg::SANDER_9 ]  = "AMBER 9: SANDER";
//	lmap[MMExternalProg::PMEMD_10]   = "AMBER 10: PMEMD";
//	lmap[MMExternalProg::SANDER_10 ] = "AMBER 10: SANDER";
//	lmap[MMExternalProg::PMEMD_12]   = "AMBER 12: PMEMD";
//	lmap[MMExternalProg::SANDER_12 ] = "AMBER 12: SANDER";
	lmap[MMExternalProg::PMEMD_18]   = "AMBER 18: PMEMD";
	lmap[MMExternalProg::TINKER_51]  = "TINKER 5.1";
	lmap[MMExternalProg::GROMACS_51] = "GROMACS 5.1";
	
	return lmap;
}

IntStrMap MMExternalProg::labels = SetMMExternalProgLbls();  

MMExternalProg::MMExternalProg()
{
	v_= PMEMD_18;	
}
	
MMExternalProg::~MMExternalProg()
{

}

int MMExternalProg::SetWithValue(int val)
{
	v_ = (Value) val;
	return TRUE;
}

int MMExternalProg::SetWithLabel(const char* label_set)
{
	IntStrMap::iterator itr;
	for(itr = labels.begin(); itr != labels.end(); itr++)
	{
		if( (*itr).second == label_set ) 
		{
			v_ = (Value) (*itr).first;
			return TRUE;
		}
	}
	return FALSE;
}




IntStrMap SetForceFieldTypeLbls()
{
	IntStrMap lmap;
	lmap[ForceFieldType::AMBER_94  ]     = "AMBER_94";
	lmap[ForceFieldType::AMBER_99_SB ]   = "AMBER_99_SB";
	lmap[ForceFieldType::AMBER_99_BSC0 ] = "AMBER_99_BSC0";
	lmap[ForceFieldType::AMBER_03]   = "AMBER_03";
	lmap[ForceFieldType::AMBER_10]   = "AMBER_10";
	lmap[ForceFieldType::AMOEBA]     = "AMOEBA";
	lmap[ForceFieldType::ARROW_5_14_CT]     = "ARROW_5.14_CT";
	lmap[ForceFieldType::ARROW_2_0] = "ARROW_2.0";
	lmap[ForceFieldType::UNKNOWN_FF] = "UNKNOWN_FF";
	
	return lmap;
}

IntStrMap ForceFieldType::labels = SetForceFieldTypeLbls();  

ForceFieldType::ForceFieldType()
{
	v_= AMBER_94;	
}
	
ForceFieldType::~ForceFieldType()
{

}

int ForceFieldType::SetWithValue(int val)
{
	v_ = (Value) val;
	return TRUE;
}

int ForceFieldType::SetWithLabel(const char* label_set)
{
	IntStrMap::iterator itr;
	for(itr = labels.begin(); itr != labels.end(); itr++)
	{
		if( (*itr).second == label_set ) 
		{
			v_ = (Value) (*itr).first;
			return TRUE;
		}
	}
	return FALSE;
}

IntStrMap SetAtomAtomContactTypeLbls()
{
	IntStrMap lmap;
	lmap[AtomContactType::HARMONIC_CNT]    = "HARMONIC";
	lmap[AtomContactType::VDW_CNT_6_12]    = "VDW_6_12";
	lmap[AtomContactType::VDW_CNT_10_12]   = "VDW_10_12";
	lmap[AtomContactType::VDW_CNT_6_12_NO_REP]   = "VDW_6_12_NO_REP";
	return lmap;
}

IntStrMap AtomContactType::labels = SetAtomAtomContactTypeLbls();  

AtomContactType::AtomContactType()
{
	v_= HARMONIC_CNT;	
}
	
AtomContactType::~AtomContactType()
{

}

int AtomContactType::SetWithValue(int val)
{
	v_ = (Value) val;
	return TRUE;
}

int AtomContactType::SetWithLabel(const char* label_set)
{
	IntStrMap::iterator itr;
	for(itr = labels.begin(); itr != labels.end(); itr++)
	{
		if( (*itr).second == label_set ) 
		{
			v_ = (Value) (*itr).first;
			return TRUE;
		}
	}
	return FALSE;
}

