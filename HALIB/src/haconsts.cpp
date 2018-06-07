//
// C++ Implementation: haconst
//
// Description:
//
//
// Author: mikola <mikola@linux>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "haconsts.h"

double HaConsts::e  = 1.60217653141414E-19; //!< elementary charge [C]
double HaConsts::Na = 6.02214151010101E23;  //!< Avogadro constant
double HaConsts::R  = 8.31447215151515;     //!< molar gas constant [J/(mol*K)]
double HaConsts::k  = 1.38065052424243E-23; //!< Boltzmann constant R/NA [J/K]
//Non-SI units accepted for use with the SI
double HaConsts::eV = 1.60217653141414E-19; //!< electron volt: (e/C) J
//Atomic and Nuclear
double HaConsts::a0    = 0.52917721081818E-10; //!< Bohr radius α/4πR∞ = 4πε0ħ2/mee2 [m]
double HaConsts::Eh    = 4.35974417757576E-18; //!< Hartree energy e2/4πε0a0 = 2R∞hc = α2mec2 [J]
double HaConsts::Eh_eV = 4.35974417757576E-18; //!< Hartree energy e2/4πε0a0 = 2R∞hc = α2mec2 [eV]

//Some units
double HaConsts::Tstd  = 298.15;               //! Standart Temperature [K]
//Convertional Coefition
//const double J_to_K= 7.24296313131313E22; //!< E = mc2 = hc/λ = hv = kT,
double HaConsts::J_to_eV = 6.24150947535353E18; //!< E = mc2 = hc/λ = hv = kT,
double HaConsts::J_to_Eh = 2.29371257393939E17; //!< J to Hartree
double HaConsts::kT_to_J = k*Tstd;//!< kT to J
double HaConsts::J_to_kT = 1.0/kT_to_J;

double HaConsts::kT_to_kJ_mol = R*Tstd*0.001;//!< kT to kJ/mol (actuall RT)
double HaConsts::kJ_mol_to_kT = 1.0/kT_to_kJ_mol;//!< kJ/mol to kT
double HaConsts::kT_to_Eh = kT_to_J*J_to_Eh;//!< kT to Hartree
double HaConsts::Bohr_to_A = a0*1.0E10;    //!< Bohr to Angstroms
double HaConsts::A_to_Bohr = 1.0/Bohr_to_A;    //!< Angstroms to Bohr
//float values have suffix f
float HaConsts::ef    = 1.602177E-19f; //!< elementary charge [C]
float HaConsts::Naf   = 6.022142E23f;  //!< Avogadro constant

double HaConsts::calorie_to_J=4.184; //!< = 4.184;             calorie in J
double HaConsts::kT_to_kcal_mol=HaConsts::kT_to_kJ_mol/HaConsts::calorie_to_J;


HaConsts::HaConsts()
{
}
HaConsts::~HaConsts()
{
}
