//
// C++ Interface: haconst
//
// Description: 
//
//
// Author: mikola <mikola@linux>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef HALIBHACONST_H
#define HALIBHACONST_H
/*! Class which contain physical constants and convertional coeficients
Most Constand based on CRC handbook of chemistry and physics, 87th Edition 2006-2007
if you want to add smth, but not sure about precise value go to http://www.hbcpnetbase.com/.
To know: float - 4bytes  3.4e +/- 38 (7 digits), double - 8bytes 1.7e +/- 308 (15 digits)
*/

#if !defined(HACONSTS_H_MIKOLA)
#define HACONSTS_H_MIKOLA

class HaConsts
{
  public:
    HaConsts();
    ~HaConsts();
    //An Abbreviated List of the CODATA Recommended Values of the Fundamental Constants of Physics and Chemistry Based on the 2002 Adjustment
    static double e;   //!<    = 1.60217653141414E-19;  elementary charge [C]
    static double Na;  //!<    = 6.02214151010101E23;   Avogadro constant
    static double R;   //!<    = 8.31447215151515;      molar gas constant [J/(mol*K)]
    static double k;   //!<  = 1.38065052424243E-23;    Boltzmann constant R/NA [J/K]
    //Non-SI units accepted for use with the SI
    static double eV;  //!<  = 1.60217653141414E-19;    electron volt: (e/C) J
    //Atomic and Nuclear
    static double a0;    //!<  = 0.52917721081818E-10;  Bohr radius α/4πR∞ = 4πε0ħ2/mee2 [m]
    static double Eh;    //!<  = 4.35974417757576E-18;  Hartree energy e2/4πε0a0 = 2R∞hc = α2mec2 [J]
    static double Eh_eV; //!< = 4.35974417757576E-18;   Hartree energy e2/4πε0a0 = 2R∞hc = α2mec2 [eV]

    //Some units
    static double Tstd;  //!<  = 298.15;                Standart Temperature [K]
    //Convertional Coefition
    //static const double J_to_K= 7.24296313131313E22;   E = mc2 = hc/λ = hv = kT,
    static double J_to_eV;  //!< = 6.24150947535353E18;  E = mc2 = hc/λ = hv = kT,
    static double J_to_Eh;  //!< = 2.29371257393939E17;  J to Hartree
    static double kT_to_J;  //!< = k*Tstd;   kT to J
    static double J_to_kT;  //!< = 1.0/kT_to_J;
    
    static double kT_to_kJ_mol; //!<  = R*Tstd*0.001;     kT to kJ/mol (actuall RT)
    static double kJ_mol_to_kT; //!< = 1.0/kT_to_kJ_mol;  kJ/mol to kT
    static double kT_to_Eh;     //!< = kT_to_J*J_to_Eh;   kT to Hartree
    static double Bohr_to_A;    //!< = a0*1.0E10;         Bohr to Angstroms
    static double A_to_Bohr;    //!< = 1.0/Bohr_to_A;     Angstroms to Bohr
    //float values have suffix f
    static float ef;            //!< = 1.602177E-19f;     elementary charge [C]
    static float Naf;           //!< = 6.022142E23f;      Avogadro constant
		
		static double calorie_to_J; //!< = 4.184;             calorie in J
		static double kT_to_kcal_mol;//! = kT_to_kJ_mol/calorie_to_J kT to kcal/mol
};

typedef HaConsts HaC;

#endif

#endif //! defined (HACONSTS_H_MIKOLA)=======

