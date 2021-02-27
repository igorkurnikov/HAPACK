/*! \file haqchem.cpp

     Basic Classes 
     to perform Quantum Chemical Calculations in HARLEM.

    \author Igor Kurnikov , University of Pittsburgh 
    \date 1998-2002 

*/

#define HAQCHEM_CPP

#pragma warning (disable:4786)

#include <mpi.h>
#include <boost/algorithm/string.hpp>
#include <wx/thread.h>

#include "errorip.h"
#include "randomip.h"
#include "hamolecule.h"
#include "hamolset.h"
#include "haatbasdb.h"
#include "hapseudopot.h"
#include "haqchem.h"
#include "hamultipole.h"
#include "tokens.h"
#include "hasurface.h"
#include "harlemapp.h"
#include "hamatdb.h"
#include "gaufile.h"
#include "hacoord.h"
#include "hagaussian.h"


#ifdef USE_IPACK  // IPACK headers

#ifdef max // to resolve a conflict with IPACK functions
#undef max
#endif

#ifdef min
#undef min
#endif

#include "const.h"
#include "basis.h"
#include "operators.h"
#include "parallel.h"
#include "f.h"
#include "randomip.h"
#endif 

using namespace harlem;

RunOptions HaQCMod::run_opt_default;

/*! \page page_qchem Quantum Chemical calculations with HARLEM

  Basic setup of a quantum chemical model of a molecular system in HARLEM is done in
HaQCMod class. The class is derived from HaCompMod class and HaTextCmdTarget
class, thus it has a basic functionality of a computational module and can accept text
commands. Three dialog classes are directly associated with HaQCMod: QChemParDlg,
LoadQCDatDlg, WaveFunAnalDlg. QChemParDlg manages setup of quantum chemical
model, prepare input files and run external quantum chemical programs. LoadQCDatDlg
dialog class manages load of the quantum chemical information such as molecular orbital
coe.cients and energies obtained by the external program. WafeFunAnalDlg dialog class
allows interactive graphical analysis of the quantum chemical data loaded in to the HaQC-
Mod class such as display of molecular orbitals with isolevel 3D contours.
HaQCMod specify gaussian basis set of the molecular system. Basis set information is
handled by HaAtBasDB class. Several basis sets are stored inside the HaAtBasDB code,
other basis set can be loaded from the files with basis set description in DALTON format
stored in a specific directory of HARLEM installation. In addition to the regular gaussian
basis set HaQCMod maintains active basis set which is a subset of the gaussian basis set
of the system and currently can correspond to valence atomic orbitals ( setid option of the
SetLocOrb function equal "VALENCE AO") or can coincide with a full basis set (setid
equal "FULL BASIS" ).
122
There are several classes which calculate matricies of one-electronic operators such as
matricies of electrical and magnetical dipole operators (classes HaOperR and HaOperRDelt
correspondingly), operators of kinetic energy (HaOperKinEner class). 
All these operators are derived from the class HaOper1e (C++ abstruction
for a general one-electronic operator) have the same member function calling interface and
calculate the their matricies by the call to DENBAS function of the GAUSSIAN utility
library, which is linked to the program. For the DENBAS function to work HaQCMod
static member function LoadGauCom fills GAUSSIAN common blocks
specifing basis set of the system. We plan to substitute GAUSSIAN libraries by freeware
gaussian integral package to remove copyright restrictions for the distribution of HARLEM.
For a divide-and-conquer analysis HaQCMod member functions InsertLocOrbSubMat and
ExtractLocOrbSubMat are especially important. These functions allow to manipulate ma-
trix of any operator in the basis of active(valence) orbitals in terms of Atomic Group
submatricies. Using these functions, for example, the matrix of the eective hamiltonian
for the donor/acceptor electronic coupling calculations is divided into atomic group-group
submatrices or assembled from group-group submatricies in divide-and-conquer calcula-
tions.
To perform calculations for optical response properties the most important classes are
HaRPAvec(an abstraction for one electron pertubation of the Hartree-Fock ground state),
HaRPAHam(an abstraction for RPA (Random Phase Approximation) for electronic Hamil-
tonian ), HaRPAResolv(an abstraction for electronic Hamiltonian resolvent in the RPA
approximation). HaRPAResolv class for any instance of HaRPAVec class specing RPA
vector X calculate RPA vector jZi = (E jZi jZi routine from GAUSSIAN utility libraries linked to HARLEM (see appropriate equations
in chapter 8. FORDIR perform appropriate convolutions of two-electron integrals with
transition density of the RPA vector. Scalar product of Z with another RPA vector Y
which is acomplished by a static function SProd of the HaRPAvec class Green Function
matrix element between X and Y, which are related to optical observables as discussed in

 */

//QCIntEngineType HaQCMod::int_engine = QCIntEngineType::INT_ENGINE_IPACK;
QCIntEngineType HaQCMod::int_engine = QCIntEngineType::INT_ENGINE_GAUSS;
int HaQCMod::max_gauss_mem = 12000000;


HaQCMod::HaQCMod(MolSet* new_phost_mset):
HaCompMod(COMP_MOD_QCHEM,new_phost_mset)
{
	SetStdParams();
}

HaQCMod::~HaQCMod()
{
	if( p_ndo_pars_db != NULL) delete p_ndo_pars_db;
	p_ndo_pars_db = NULL;
}

void HaQCMod::SetStdParams()
{
	charge=0;
	mult=1;
	SetHF();

	load_mo_flag = true;
	m_bas_name = "3-21G";
	m_loc_orb_set_id = "FULL_ATOMIC_BASIS";
	
	SetUseSymmetry(TRUE);
	SetCalcPolar(false);

	qpole.resize(6); qpole = 0.0;

	MO_coef.newsize(0,0);
	MOene.newsize(0); 

	stop_calc_flag = FALSE;

	m_grid_size = 21;
	distr_ext_chrg = 0.0;

	max_scf_iter = 200;
	temp0_fermi  = 0.0;
	iter_temp = 20;

	max_it_avg = 20;    
	max_it_noavg = 0;
    iuhf = FALSE;
	conv_dm = 1.0E-7;

    guess_only = FALSE;

	ActBas = NULL;
	allocated_act_basis = TRUE;

	SetExtChCrdOffset( 0.01 );


#if defined(GAUSSVER) 
	kjob_.numprc=1; // For GAUSSIAN set number of processors to 1, is it necessary?
#endif
}

void HaQCMod::set_max_gauss_mem(int new_max_mem)
{
	max_gauss_mem = new_max_mem;
}


bool HaQCMod::Print_info(ostream &sout, const int level) const
// Output information about the Quantum Chemical Model
{
	AtBasisType::const_iterator itrb=AtBasis.at_bas_vec.begin();
	for(; itrb != AtBasis.at_bas_vec.end(); itrb++)
	{
		char buf[256];
		const HaAtom* aptr= (*itrb).GetAtHost();
		if(aptr == NULL)
		{
			strcpy(buf,"NULL_ATOM");
		}
		else
		{
			aptr->FillRef(buf);
		}
		sout << "Basis set for Atom " << buf << endl;
        (*itrb).Print_info(sout,level);
		sout << endl;
	}
	return true;
}

int HaQCMod::SaveXMLToStream(std::ostream& os, const harlem::SaveOptions* popt ) const
{
	os << "<qchem> mset=\"" << GetMolSet()->GetName() << "\">" << std::endl;
	os << "<atbas name=\""<< GetBasName() << "\"/>" << std::endl;
	os << "<loc_orb_bas id=\"" << GetLocOrbSetID() << "\"/>" << std::endl;
	os << "</qchem>" << std::endl;
	return TRUE;
}

int HaQCMod::GetNumCnt() const
{
	if(phost_mset == NULL) return 0;

	return(phost_mset->GetNAtoms());

}

bool HaQCMod::GetCntCharges(HaVec_double& charges) const
{
	if(phost_mset == NULL) 
	{
		cerr << " Error in HaQCMod::GetCntCharges " << endl;
		cerr << " phost_mset == NULL " << endl;
		return false;
	}
	int nn= phost_mset->GetNAtoms();
	charges.newsize(nn);
	const HaAtom* aptr;
	int i=1;
	
	AtomIteratorMolSet aitr(phost_mset);
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		charges(i)=aptr->GetElemNo();
		i++;
	}
	return true;
}

QCIntEngineType HaQCMod::GetQCIntEngine()
{
	return int_engine;
}
	
void HaQCMod::SetQCIntEngine(const QCIntEngineType& int_engine_new)
{
	int_engine = int_engine_new;
}



void
HaQCMod::SetCharge(int NewCharge) 
{
	this->charge = NewCharge;
}

void
HaQCMod::SetMult(int NewMult) 
{
	this->mult=NewMult ;
}


int HaQCMod::GetCharge() const
{
	return(this->charge);
}

int HaQCMod::GetMult() const
{
	return(this->mult);
}

int HaQCMod::GetNelectr() const
{
    int ne=0;
	AtBasisType::const_iterator bitr;

    for(bitr= AtBasis.at_bas_vec.begin(); bitr != AtBasis.at_bas_vec.end(); bitr++)
	{
		ne += (*bitr).GetNumElectr();
    }
	ne = ne - GetCharge();
    return(ne);
}

void HaQCMod::SetField(double el_x, double el_y, double el_z)
{
	el_field[0] = el_x;
	el_field[1] = el_y;
	el_field[2] = el_z;
}

void HaQCMod::SetFieldAU(double el_x, double el_y, double el_z)
{
	el_field[0] = el_x/(BOHR_TO_ANG*BOHR_TO_ANG);
	el_field[1] = el_y/(BOHR_TO_ANG*BOHR_TO_ANG);
	el_field[2] = el_z/(BOHR_TO_ANG*BOHR_TO_ANG);
}

int HaQCMod::GetNumMO() const
{
	return ( MOene.size());
}

bool HaQCMod::InitBasis(const std::string& bname)  
// initialize basis set with a given name on all atoms of the system
{
	AtBasis.Clear();
	ActBas = NULL;
	ovlp_mat.clear();
	MO_coef.clear();
	MOene.clear();

	if( ActBas && allocated_act_basis)
	{
		delete ActBas;
		ActBas = NULL;
	}
	m_bas_name = bname;
	
	AtBasis.InitForMolSet(bname.c_str(), phost_mset); 

	HaMat_double ss_test;
	HaBasisSet::CalcOvlpMat( &AtBasis, &AtBasis, ss_test);
//	PrintLog(" HaQCMod::InitBasis()  ss.r0(0,100) = %12.6f \n", ss_test.r0(0,100) );

	return true;
}

bool HaQCMod::InitBasis(GauFile& gfile)
{
#if defined(GAUSSVER)

	int iunit = -GauFile::IO_basis; 
	int iaxx=2+7*MAXSHL;
    int len  = ((15*MAXSHL) + (int)intowp_(&iaxx));
    int offset = 0;
	GauFile gbas;
	void* pdata= (void*) &gbas;
	int res= gfile.fileio("read",iunit,len,pdata,offset);
	return true;
#else
	PrintLog("HaQCMod::InitBasis(GauFile& gfile) \n");
	PrintLog("GAUSSIAN Library is not present \n");
	return false;
#endif 
}

int HaQCMod::FBasFunPos(const HaAtom* ref_aptr) const
{
	if(ref_aptr == NULL)
	{
		ErrorInMod("HaQCMod::FBasFunPos()",
			       " Atom Pointer ref_aptr is NULL ");
		return(-1);

	}
	const HaAtom* aptr;

	AtBasisType::const_iterator bsitr;
	int ifst=1;

	for(bsitr=AtBasis.at_bas_vec.begin() ; bsitr != AtBasis.at_bas_vec.end(); bsitr++)
	{
		aptr= (*bsitr).GetAtHost();
		if(aptr == ref_aptr)
			break;
		ifst += (*bsitr).GetNBfunc(); 
	}


	if(bsitr == AtBasis.at_bas_vec.end()) 
	{
		cerr << " Error in HaQCMod::FBasFunPos() : " << endl ;
		char buf[256];
		ref_aptr->FillRef(buf);
		cerr << " Atom Basis of atom " << buf << endl;
		cerr << " is not found " << endl;
		return (-1);
	}

	return ifst;
}

const HaAtom* HaQCMod::GetAtomOfAO(const int idx_AO) const
{
	assert(idx_AO > 0);
	if(AtBasis.at_bas_vec.empty()) return NULL;

	AtBasisType::const_iterator itr;

	int nbas=0;
	
	for( itr= AtBasis.at_bas_vec.begin(); itr != AtBasis.at_bas_vec.end(); itr++)
	{
		nbas += (*itr).GetNBfunc();	
		if(nbas >= idx_AO) 
			return ((*itr).GetAtHost());
	}
	return NULL;
	
}

int HaQCMod::GetNBfunc() const
// return the number of basis functions:
{
	if(AtBasis.at_bas_vec.empty()) return 0;

	AtBasisType::const_iterator itr;

	int nbas=0;
	
	for( itr=AtBasis.at_bas_vec.begin(); itr != AtBasis.at_bas_vec.end(); itr++)
	{
		int n1= (*itr).GetNBfunc();
		nbas += n1;	
	}
	return nbas;
}

int HaQCMod::ProjMatToActBas(HaMat_double& fmat, HaMat_double& fmat_lb)
{
	if(ovlp_mat.num_rows() == 0) InitBasOvlp();
	
	int nb  = GetNBfunc();
	int nab = GetNActiveOrb();
	if( nab <= 0 || (nab > nb)) return FALSE;
	
	HaMat_double cl(nb,nab);
	
	if((void*)ActBas == (void*) &AtBasis )
	{
		fmat_lb = fmat;
		return TRUE;
	}

	std::string bas_type_1 = ActBas->GetClassName();

	if(bas_type_1 == "LinCombOrb3D")
	{
		LinCombOrb3D* lcmb = (LinCombOrb3D*)ActBas;
        LinCombOrb3D::Eval1eOp(lcmb,lcmb,fmat,fmat_lb);
	}

	return true;
}

int HaQCMod::CalcEPfromMO(HaMat_double& gm, double ene)
{
	int nb= GetNBfunc();
	
	if( nb <= 0 ) 
	{
		PrintLog(" HaQCMod::CalcEPfromMO \n ");
		PrintLog(" Basis set is not set for the host molecule \n");
		return FALSE;
	}

	if( MO_coef.num_cols() == 0 || MOene.size() == 0)
	{
		PrintLog(" Error in: HaElPropHF::Recalc() : \n");
		PrintLog(" Molecular orbitals or/and their energies are not set \n");  
		return false;
	}

	gm.newsize(nb,nb);
	gm = 0.0;

	int j,k,l;

	int num_mo = MOene.size();

	// Build Electronic Propagator matrix from Molecular Orbitals
	for(j=1; j <= num_mo; j++)
	{
		double r_et_min_ej= 1.0/(ene - MOene(j) );
		for(k=1; k <= nb; k++)
		{
			for( l=1; l <= nb; l++)
			{
				gm(k,l)+= MO_coef(k,j) * MO_coef(l,j) * r_et_min_ej;
			}
		}
	}

	HaMat_double& ss = GetOvlpMat();

	HaMat_double scr;

	matmult(scr,ss,gm);
	matmult(gm,scr,ss); // gm = S * ( sum_j  C_kj*C_lj/(E_T - E_j) ) * S  

	return TRUE;
}

bool HaQCMod::SetWaveFunType(const char* wave_fun_type_cs)
{
	std::string str_wave_fun_type = wave_fun_type_cs;
	boost::to_upper(str_wave_fun_type);

	if(str_wave_fun_type == "HARTREE-FOCK" )
	{
		wave_fun_type = harlem::qc::HARTREE_FOCK;
		return true;
	}
	if(str_wave_fun_type == "CNDO/2")
	{
		wave_fun_type=harlem::qc::NDO;
		ndo_method = harlem::qc::CNDO_2;
		return true;
	}
	if(str_wave_fun_type == "INDO/2")
	{
		wave_fun_type= harlem::qc::NDO;
		ndo_method = harlem::qc::INDO_2;
		return true;
	}
	if( str_wave_fun_type == "ZINDO/1")
	{
		wave_fun_type=harlem::qc::NDO;
		ndo_method = harlem::qc::ZINDO_1;
		return true;
	}
	if( str_wave_fun_type == "ZINDO/S")
	{
		wave_fun_type=harlem::qc::NDO;
		ndo_method = harlem::qc::ZINDO_S;
		return true;
	}
	if( str_wave_fun_type == "EXTENDED_HUCKEL" )
	{
		wave_fun_type= harlem::qc::EXTENDED_HUCKEL;
		return true;
	}
	if( str_wave_fun_type == "MP2")
	{
		wave_fun_type= harlem::qc::MP2;
		return true;
	}
	if( str_wave_fun_type == "DFT_B3LYP" ) 
	{
		wave_fun_type= harlem::qc::DFT_B3LYP;
		return true;
	}
	PrintLog(" Error in HaQCMod::SetWaveFunType() \n");
	PrintLog(" Unknown Wave Function type: %s \n", str_wave_fun_type.c_str()); 
	PrintLog(" Wave function type set to Hartree-Fock \n");

	wave_fun_type = harlem::qc::HARTREE_FOCK;

	return false;
}

int HaQCMod::RunExtHuckel()
{
    InitBasis("MINB6G");
    InitLocOrb("FULL_ATOMIC_BASIS");
	HaMat_double& ss = GetOvlpMat();
	InitHuckHam(huck_ham,ss,AtBasis);
	HaMat_double::DiagMat(huck_ham,ss,MO_coef,MOene);
	return TRUE;
}

void HaQCMod::SetSinglePtCalc()
{
	calc_type = SINGLE_PT_CALC;
}

void HaQCMod::SetEneMinCalc()
{
	calc_type = ENERGY_MIN_CALC;
}

void HaQCMod::SetTransStateCalc()
{
	calc_type = TRANS_STATE_CALC;
}

bool HaQCMod::IsSinglePtCalc() const
{
	return (calc_type == SINGLE_PT_CALC);
}

bool HaQCMod::IsEneMinCalc() const
{
	return (calc_type == ENERGY_MIN_CALC);
}

bool HaQCMod::IsTransStateCalc() const
{
	return (calc_type == TRANS_STATE_CALC);
}

int HaQCMod::CalcEnergy()
{
	SetSinglePtCalc();
	int ires = Run();
	return ires;
}

int HaQCMod::RunMinEne()
{
	SetEneMinCalc();	
	int ires = Run();
	return ires;
}

bool HaQCMod::IsGenBasisSet() const
{
	return AtBasis.IsGeneric();
}
	
void HaQCMod::SetBasisSetGen()
{
	AtBasis.SetGeneric();
}

bool HaQCMod::IsUsingSymmetry() const
{
	return (use_symmetry != 0);
}

void HaQCMod::SetUseSymmetry( int use_symmetry_par )
{
	use_symmetry = use_symmetry_par;
}

void HaQCMod::SetHF()
{
	wave_fun_type = harlem::qc::HARTREE_FOCK;
}

void HaQCMod::SetMP2()
{
	wave_fun_type = harlem::qc::MP2;
}

void HaQCMod::SetDFT()
{
	wave_fun_type = harlem::qc::DFT_B3LYP;
}

void HaQCMod::SetB3LYP()
{
	wave_fun_type = harlem::qc::DFT_B3LYP;
}

void HaQCMod::SetCCSD_T()
{
	wave_fun_type = harlem::qc::CCSD_T;
}

bool HaQCMod::IsHF() const
{
	return (wave_fun_type == harlem::qc::HARTREE_FOCK );
}

bool HaQCMod::IsDFT() const
{
	return (wave_fun_type == harlem::qc::DFT_B3LYP );
}

bool HaQCMod::IsB3LYP() const
{
	return (wave_fun_type == harlem::qc::DFT_B3LYP );
}
	
bool HaQCMod::IsMP2() const
{
	return (wave_fun_type == harlem::qc::MP2 );
}

bool HaQCMod::IsCCSD_T() const
{
	return ( wave_fun_type == harlem::qc::CCSD_T );
}

bool HaQCMod::IsUsingSCRF() const
{
	return ( scrf_method != harlem::qc::SCRF_NO );
}

void HaQCMod::SetSCRF()
{
	scrf_method = harlem::qc::SCRF_PCM;
}

void HaQCMod::SetSCRFMethod(int set_method_par )
{
	if( set_method_par == (int) harlem::qc::SCRF_PCM )
	{
		scrf_method = harlem::qc::SCRF_PCM;
	}
	else
	{
		scrf_method = harlem::qc::SCRF_NO;
	}
}

ZMatCrd* HaQCMod::GetZMat() 
{ 
	return GetMolSet()->GetZMat(); 
}

void HaQCMod::SetPrefix( const std::string& prefix_par )
{
	prefix = prefix_par;
}

std::string HaQCMod::GetPrefix() const
{
	std::string prefix_loc = prefix;
	if( prefix_loc.empty() ) prefix_loc = GetMolSet()->GetName();
	return prefix_loc;
}

bool HaQCMod::UsePseudoPot() const
{
	AtBasisType::const_iterator bitr;
	for(bitr= AtBasis.at_bas_vec.begin(); bitr != AtBasis.at_bas_vec.end(); bitr++)
	{
		if((*bitr).IsSetPseudoPot())
			return true;
	}
	return false;
}

bool HaQCMod::InitLocOrb(const char* setid)
// Initiate local orbitals with a method given by setid idenfificator
{	
	HaAtom* aptr;

	if(ActBas && allocated_act_basis) 
	{
		delete ActBas;
		ActBas = NULL;
	}
	m_loc_orb_set_id=setid;
	
	if( IsLocOrbFullBasis() ) // If Active basis coincide with Full basis set 
	{
		allocated_act_basis = FALSE;
		ActBas = &AtBasis;
		return true;
	}
	
	if( m_loc_orb_set_id == "VALENCE_AO")
	{
		ActBas = new GauBasisSet();
		((GauBasisSet*)ActBas)->InitForMolSet("MINB6G",GetMolSet());
		return true;
	}
	return false;

	int nb = AtBasis.GetNBfunc();

	HaMat_double cf_axx(nb,nb,0.0); // orbital expansion coefficients
	StrVec orb_lbs;  // orbital labels

	int nabf = 0;  // the number of active basis functions
	int i_fst = 0;

	VecPtr atoms_orb;
	map<std::string, HaMat_double, less<std::string> > tr_cf;

	AtBasisType::iterator bitr = AtBasis.at_bas_vec.begin();
	for(; bitr != AtBasis.at_bas_vec.end(); bitr++)
	{
		aptr= (HaAtom*) (*bitr).GetAtHost();
		if(aptr == NULL )
		{
			cerr << "Error in HaQCMod::InitLocOrb() " << endl;
			cerr << "Cannot find host Atom of the Atom Basis Set " << endl;
			return false;
		}

		std::string bname ( (*bitr).GetBasName() );
		std::string atsymb= aptr->GetStdSymbol();

		int ic = tr_cf.count(atsymb);
		int nf, nv;

		int idx_lbl = -1;
		int elem = aptr->GetElemNo();
		if(elem <= 2) idx_lbl = 0;
		if(elem > 2 && elem <= 10) idx_lbl = 1;
		if(elem > 10 && elem <= 18) idx_lbl = 2;
		if(elem > 18 && elem <= 20) idx_lbl = 3;
		if(elem > 20 && elem <= 36) idx_lbl = 4;
		if(elem > 36 && elem <= 48) idx_lbl = 5;

		if(elem > 48) 
		{
			PrintLog(" Can't set valence orbitals for element > %d \n", elem);
			continue;
		}

		HaMat_double cf_mb;
//		StrVec lbls;
		if(ic == 0)
		{
			GauBasisSet bset_min;
			GauAtomBasis min_bas("MINB6G",aptr);
			bset_min.at_bas_vec.push_back(min_bas);

			GauBasisSet bset_at;
			GauAtomBasis atbas;
			atbas.copy_from(*bitr);
			atbas.SetAtHost(aptr);
			bset_at.at_bas_vec.push_back(atbas);

			HaMat_double s_bb,s_bm,tmp;
			HaBasisSet::CalcOvlpMat(&bset_at,&bset_at,s_bb);
			tmp = s_bb;
			HaBasisSet::CalcOvlpMat(&bset_at,&bset_min,s_bm);

			HaMat_double::solv_lin_syst_1(tmp,s_bm);

			matmult(tmp,s_bb,s_bm);
			matmult_T1(s_bb,s_bm,tmp);

			nv = s_bm.num_cols();
			nf = s_bm.num_rows();

			int i,j;
			for( i = 0; i < nv; i++)
			{
				double scale = s_bb.GetVal_idx0(i,i);
				for(j = 0; j < nf; j++)
				{
					s_bm.SetVal_idx0(j,i, scale*s_bm.GetVal_idx0(j,i));
				}
//				lbls.push_back(min_bas.GetLabel(i));
			}
			tr_cf[atsymb] = s_bm;
			cf_mb = s_bm;
		}
		else
		{
			cf_mb = tr_cf[atsymb];
//			lbls = lbl_vecs[atsymb];
			nf = cf_mb.num_rows(); 
			nv = cf_mb.num_cols();
		}

		int i,j;
		for( i = 0; i < nv; i++)
		{
			atoms_orb.push_back(aptr);
			for(j =0; j < nf; j++)
			{
				cf_axx.SetVal_idx0( i_fst + j,nabf,cf_mb.GetVal_idx0(j,i));
			}
//			std::string lbl = min_bas.GetLabel(i);
//			orb_lbs.push_back(lbl);
			nabf++;
		}
		i_fst += (*bitr).GetNBfunc();
	}

	ActBas = new LinCombOrb3D(&AtBasis, cf_axx.begin(), nabf);
	((LinCombOrb3D*)ActBas)->ids = orb_lbs;
    ((LinCombOrb3D*)ActBas)->at_ptr = atoms_orb;

	return true;
}

bool HaQCMod::IsLocOrbFullBasis()
{
	return( m_loc_orb_set_id == "FULL_ATOMIC_BASIS");
}

int HaQCMod::GetNActiveOrb() const
{
	if(ActBas == NULL) return 0;
	return(ActBas->GetNBfunc());	
}

int HaQCMod::GetLocOrbIdxOfGrp(const std::string& gid, HaVec_int & lorb_idx ) const
{
	ChemGroup* gptr=phost_mset->GetChemGroupByID( gid );
	if(gptr == NULL)
	{
		cerr << " Error in HaQCMod::GetLocOrbIdxOfGrp() " << endl;
		cerr << " Can't find group with ID: " << gid << endl;
		return False;
	}

	int nab=ActBas->GetNBfunc();

	int i;
	list<int> tmp_list;
	for( i=1; i <= nab; i++)
	{
		HaAtom* aptr = (HaAtom*) ActBas->GetHostPt(i-1);
		if( gptr->IsMember(aptr) )
		{
			tmp_list.push_back(i);
		}
	}

	int ns=tmp_list.size();
	lorb_idx.newsize(ns);

	list<int>::iterator itr;
	i=1;
	for( itr= tmp_list.begin(); itr != tmp_list.end(); itr++ )
	{
		lorb_idx(i)=*itr;
		i++;
	}

	return true;

}

bool HaQCMod::ExtractLocOrbSubMat(const std::string & gid1, const std::string & gid2,
		                    const HaMat_double & ActOrbMat, 
							HaMat_double & ActOrbSubMat) const
{
	int nactorb= GetNActiveOrb();
	if( ActOrbMat.num_rows() != nactorb || ActOrbMat.num_cols() != nactorb) 
	{
		PrintLog(" Error in HaQCMod::ExtractLocOrbSubMat() \n");
		PrintLog(" the dimensions of ActOrbMat do not coincide with the number of active orbitals");
        return false;
	}
	
	HaVec_int ilgr1, ilgr2; // indexes of active orbitals of two groups
	GetLocOrbIdxOfGrp(gid1, ilgr1 );
	GetLocOrbIdxOfGrp(gid2, ilgr2 );
	int n1=ilgr1.size();
	int n2=ilgr2.size();
	
	ActOrbSubMat.newsize(n1,n2);
	for(int i=1; i <= n1; i++)
	{
		for(int j=1; j <= n2; j++)  
		{
			ActOrbSubMat(i,j)= ActOrbMat( ilgr1(i), ilgr2(j) );
		}
	}
	return true;
}


bool HaQCMod::InsertLocOrbSubMat(const std::string & gid1, const std::string & gid2,
		                    HaMat_double & ActOrbMat, 
							const HaMat_double & ActOrbSubMat) const 
{
	int nab= GetNActiveOrb();
	if( ActOrbMat.num_rows() != nab || ActOrbMat.num_cols() != nab) 
	{
		PrintLog(" Error in HaQCMod::InsertLocOrbSubMat() \n");
		PrintLog(" the dimensions of ActOrbMat do not coincide with the number of active orbitals");
		return false;
	}

	HaVec_int ilgr1, ilgr2; // indexes of active orbitals of two groups
	GetLocOrbIdxOfGrp(gid1, ilgr1 );
	GetLocOrbIdxOfGrp(gid2, ilgr2 );
	int n1=ilgr1.size();
	int n2=ilgr2.size();

	if( ActOrbSubMat.num_rows() != n1 )
	{
		cerr << " Error in HaQCMod::InsertLocOrbSubMat() " << endl;
		cerr << " Number of local orbitals n1= " << n1 << " in the group Id= " << gid1
			 <<  endl;
		cerr << " Doesn't equal to the number of rows in ActOrbSubMat " 
			 << ActOrbSubMat.num_rows() << endl;
		cerr << endl;
		return false;
	}

	if( ActOrbSubMat.num_cols() != n2 )
	{
		cerr << " Error in HaQCMod::InsertLocOrbSubMat() " << endl;
		cerr << " Number of local orbitals n2= " << n2 << " in the group Id= " << gid2
			 <<  endl;
		cerr << " Doesn't equal to the number of coloumns in ActOrbSubMat " 
			 << ActOrbSubMat.num_cols() << endl;
		cerr << endl;
		return false;
	}
	
	for(int i=1; i <= n1; i++)
	{
		for(int j=1; j <= n2; j++)  
		{
			ActOrbMat( ilgr1(i), ilgr2(j) ) = ActOrbSubMat(i,j);
		}
	}
	return true;
}

bool HaQCMod::EvalLinCombOnGrid(const HaVec_double& orb_coef, ArrayOrb3D& bas_set, HaField3D& mo_field)
{
	if(bas_set.GetNBfunc() != orb_coef.size() )
	{
		cerr << " Error in HaQCMod::BuildMOgrid() " << endl;
		cerr << " The size of Orbital coef matrix " << orb_coef.size() 
			 << " is not equal to the number of basis functions " << bas_set.GetNBfunc() << endl;
		return false;
	}
	
	mo_field.FillZeros();

	int icur=1;
	int i,j,k;
	int nx,ny,nz;

	nx= mo_field.GetNx();
	ny= mo_field.GetNy();
	nz= mo_field.GetNz();

	std::string bas_type = bas_set.GetClassName();
	
	if( bas_type == "GauBasisSet" )
	{
		GauBasisSet& gau_bas = (GauBasisSet&)bas_set;
		AtBasisType::iterator bitr;
		for(bitr = gau_bas.at_bas_vec.begin(); bitr != gau_bas.at_bas_vec.end(); bitr++)
		{
			ShellsType::iterator shitr;
			const HaAtom* aptr= (*bitr).GetAtHost();
			for(shitr = (*bitr).Shells.begin(); shitr != (*bitr).Shells.end(); shitr++)
			{
				int nf= (*shitr).GetNBfunc(); 

				for( k=0; k < nz; k++)
				{
					for ( j = 0; j < ny; j++)
					{
						for( i=0; i < nx; i++)
						{
							float x,y,z;
							bool result= mo_field.GetXYZ(x,y,z,i,j,k);
							float* pval = mo_field.GetValPtr(i,j,k);
							if(pval != NULL)
							{
								float fv = (float)(*shitr).EvalLinCombInPoint( x- aptr->GetX(), y - aptr->GetY(), z - aptr->GetZ(),
									&(orb_coef(icur)));
								(*pval)+= fv;
							}
						}
					}
				}
				icur+= nf;
			}
		}
	}
	else if( bas_type == "LinCombOrb3D" )
	{
		LinCombOrb3D& lcmb =  (LinCombOrb3D&) bas_set;
		int norb = lcmb.GetNOrbs();
		int nb_in = lcmb.bas->GetNBfunc();
		HaVec_double cf_red(nb_in,0.0);

		int i,j;
		for(i=0; i < nb_in; i++)
		{
			for( j = 0; j < norb; j++)
			{
				cf_red[i] += orb_coef[j]*lcmb.coef.GetVal_idx0(i,j);
			}
		}
		EvalLinCombOnGrid(cf_red,*(lcmb.bas), mo_field);
	}

	if(False) // Print 3D field
	{
		char buf[256];
		
		cout << " Field 3D: " << endl;
		for( k = 0; k < nz; k++ )
		{
			cout << " Layer #= " << k << endl;
			for( j = 0; j < ny; j++ ) 
			{
				for (i = 0; i < nx ; i++ )
				{	
					sprintf(buf,"%12.6e  ", *(mo_field.GetValPtr(i,j,k)) );
					cout << buf;
				}
				cout << endl;
			}
		}
	}
	return true;
}
	
bool HaQCMod::CreateMOcontour(const int imo, const double mo_isolvl, const int ngrid)
{
	bool result=false;
	char buf[256];
	if( MO_coef.num_cols() >= imo && MOene.size() > 0 && imo <= MOene.size() ) 
	{
		int nb = GetNBfunc();
		HaVec_double orb_coef( nb );
		if( nb != MO_coef.num_rows() )
		{
			cerr << " Error in HaQCMod::CreateMOcontour() " << endl;
			cerr << " The number of basis functions " << nb << 
				    " is not equal to the size of MOs " <<  MO_coef.num_rows() << endl;
			return false;
		}
                int i;
		for(i=1; i <= nb; i++)
		{
			orb_coef(i)= MO_coef(i,imo);
		}
		VecPtr surfs_ptr = CreateOrbContour(orb_coef, AtBasis, mo_isolvl , ngrid);
		int nsurf = surfs_ptr.size();
		for(i = 0; i < nsurf; i++)
		{
			result = true;
			sprintf(buf,"%d",imo);
			std::string str_idx = buf;
			boost::trim(str_idx);
		
			std::string cname = "MO_" + str_idx;
			if( i == 0) cname += "_POS";
			if( i == 1) cname += "_NEG";
			HaDisplayedSurface* pcnt = (HaDisplayedSurface*)surfs_ptr[i];
			pcnt->SetObjName(cname.c_str());
		}

	}
	return result;
}

VecPtr HaQCMod::CreateOrbContour(const HaVec_double& orb_coef, ArrayOrb3D& bas_set, const double mo_isolvl, const int ngrid)
{
	MolSet* pmset = GetMolSet();
	double xmin, xmax, ymin, ymax, zmin, zmax;
	pmset->GetMinMaxCrd( xmin, ymin, zmin, xmax, ymax, zmax);
	HaField3D mo_field;
	mo_field.SetDimensions(ngrid, ngrid, ngrid);
	mo_field.SetGridCornersCoord(xmin - 3.0, ymin - 3.0, zmin - 3.0, xmax + 3.0, ymax + 3.0, zmax + 3.0 );

	EvalLinCombOnGrid(orb_coef, bas_set, mo_field);
	
	VecPtr surfs;

	HaDisplayedSurface* sptr;

// Positive Potential Isosurface
	sptr = new HaDisplayedSurface();
	if(sptr == NULL)
		return surfs;

	pmset->AddObject3D(sptr);
	surfs.push_back((void*)sptr);

	bool result;
	result = sptr->calc_isosurf(&mo_field, mo_isolvl);
	
	sptr->ColourUniform(255,0,0);

// Negative Potential Isosurface

	sptr = new HaDisplayedSurface();
	if(sptr == NULL)
		return surfs;

	pmset->AddObject3D(sptr);
	surfs.push_back((void*)sptr);

	result = sptr->calc_isosurf(&mo_field, -mo_isolvl);
	sptr->ColourUniform(0,0,255);

	return surfs;
}

bool HaQCMod::BuildFockMatFromMOs(HaMat_double& fock_matrix, double cut_ene)
{
	if(debug_level > 5) 
		PrintLog(" Enter HaQCMod::BuildFockMatFromMOs() \n");

	int nb = GetNBfunc();
	if( nb <= 0 )
	{
		ErrorInMod("HaQCMod::BuildFockMatFromMOs()",
	    " Build of Fockian failed. Number of basis function <= 0");
		return false;
	}
	
	if( MO_coef.num_cols() == 0 || MOene.size() == 0)
	{
		ErrorInMod("HaQCMod::BuildFockMatFromMOs()",
	    " Molecular orbitals or/and their energies are not set ");
		return false;
	}
	
	fock_matrix.newsize(nb,nb);
	fock_matrix = 0.0;
	
	int num_mo = MOene.size();
	int k,l,j;

	// Build Electronic Propagator matrix from Molecular Orbitals
	for(j=1; j <= num_mo; j++)
	{
		double ej= MOene(j) ;
		if(ej < cut_ene)
			continue;
		for(k=1; k <= nb; k++)
		{
			for( l=1; l <= nb; l++)
			{
				fock_matrix(k,l)+= MO_coef(k,j) * MO_coef(l,j) * ej;
			}
		}
	}

	if(wave_fun_type == harlem::qc::NDO) return true;

	HaMat_double scr;
	HaMat_double& ss = GetOvlpMat();
	matmult(scr,ss,fock_matrix);
	matmult(fock_matrix,scr,ss); // fock_mat = S * ( sum_j  C_kj*C_lj*E_j)  * S  

	return true;
}

// static StrDoubleMap huck_pars_vela;
static StrDoubleMap huck_pars_std;

int HaQCMod::InitHuckParsStd()
{
	huck_pars_std["H_1S"]  = 13.6/HARTREE_TO_EV;  // HOFFMANN_JCP_63_39_1397
	huck_pars_std["HE_1S"] = 60.0/HARTREE_TO_EV; // Igor addition for completeness 

    huck_pars_std["LI_2S"] = 5.4/HARTREE_TO_EV;  // ANDERSON_HOFFMANN_JCP(74)_60_4271
	huck_pars_std["LI_2P"] = 3.5/HARTREE_TO_EV;  // ANDERSON_HOFFMANN_JCP(74)_60_4271
    huck_pars_std["BE_2S"] = 10.0/HARTREE_TO_EV; // CANADELL_INORG_CHEM(83)_22_3856
	huck_pars_std["BE_2P"] =  6.0/HARTREE_TO_EV; // CANADELL_INORG_CHEM(83)_22_3856
    huck_pars_std["B_2S"]  = 15.2/HARTREE_TO_EV; // ANDERSON_HOFFMANN_JCP(74)_60_4271
	huck_pars_std["B_2P"]  =  8.5/HARTREE_TO_EV; // ANDERSON_HOFFMANN_JCP(74)_60_4271
    huck_pars_std["C_2S"]  = 21.4/HARTREE_TO_EV; // HOFFMANN_JCP(63)_39_1397
	huck_pars_std["C_2P"]  = 11.4/HARTREE_TO_EV; // HOFFMANN_JCP(63)_39_1397
    huck_pars_std["N_2S"]  = 26.0/HARTREE_TO_EV; // HOFFMANN_JCP(63)_39_1397
	huck_pars_std["N_2P"]  = 13.4/HARTREE_TO_EV; // HOFFMANN_JCP(63)_39_1397
    huck_pars_std["O_2S"]  = 32.3/HARTREE_TO_EV; // HOFFMANN_JCP(63)_39_1397
	huck_pars_std["O_2P"]  = 14.8/HARTREE_TO_EV; // HOFFMANN_JCP(63)_39_1397
    huck_pars_std["F_2S"]  = 40.0/HARTREE_TO_EV; // ANDERSON_HOFFMANN_JCP(74)_60_4271
	huck_pars_std["F_2P"]  = 18.1/HARTREE_TO_EV; // ANDERSON_HOFFMANN_JCP(74)_60_4271
    huck_pars_std["NE_2S"]  = 60.0/HARTREE_TO_EV; // Igor addition for completeness
	huck_pars_std["NE_2P"]  = 60.0/HARTREE_TO_EV; // Igor addition for completeness

    huck_pars_std["NA_3S"] = 5.1/HARTREE_TO_EV;   // ANDERSON_HOFFMANN_JCP(74)_60_4271
	huck_pars_std["NA_3P"] = 3.0/HARTREE_TO_EV;   // ANDERSON_HOFFMANN_JCP(74)_60_4271
    huck_pars_std["MG_3S"] = 8.201/HARTREE_TO_EV;  // FROM VELA_JPC(88)_92_5688 - For completeness
	huck_pars_std["MG_3P"] = 4.641/HARTREE_TO_EV;  // FROM VELA_JPC(88)_92_5688 - For completeness
    huck_pars_std["AL_3S"] = 12.3/HARTREE_TO_EV;   // ANDERSON_HOFFMANN_JCP(74)_60_4271
	huck_pars_std["AL_3P"] = 6.5/HARTREE_TO_EV;    // ANDERSON_HOFFMANN_JCP(74)_60_4271
    huck_pars_std["SI_3S"] = 17.3/HARTREE_TO_EV;   // TRONG_HOFFMANN_JACS(78)_100_110   
	huck_pars_std["SI_3P"] = 9.2/HARTREE_TO_EV;    // TRONG_HOFFMANN_JACS(78)_100_110
    huck_pars_std["P_3S"]  = 18.6/HARTREE_TO_EV;   // SUMMERVILLE_HOFFMAN_JACS(76)_98_7240
	huck_pars_std["P_3P"]  = 14.0/HARTREE_TO_EV;   // SUMMERVILLE_HOFFMAN_JACS(76)_98_7240
    huck_pars_std["S_3S"]  = 20.0/HARTREE_TO_EV;    // CHEN_HOFFMANN_JACS(76)_98_1647
	huck_pars_std["S_3P"]  = 11.0/HARTREE_TO_EV;    // CHEN_HOFFMANN_JACS(76)_98_1647
    huck_pars_std["CL_3S"] = 26.3/HARTREE_TO_EV;   // SUMMERVILLE_HOFFMAN_JACS(76)_98_7240
	huck_pars_std["CL_3P"] = 14.2/HARTREE_TO_EV;   // SUMMERVILLE_HOFFMAN_JACS(76)_98_7240

    huck_pars_std["FE_4S"] =  7.6/HARTREE_TO_EV;   // SAILLARD_HOFFMAN_JACS(84)_106_2006
	huck_pars_std["FE_4P"] =  3.8/HARTREE_TO_EV;   // SAILLARD_HOFFMAN_JACS(84)_106_2006
	huck_pars_std["FE_3D"] =  9.2/HARTREE_TO_EV;   // SAILLARD_HOFFMAN_JACS(84)_106_2006
    huck_pars_std["CO_4S"] =  7.8/HARTREE_TO_EV;   // SAILLARD_HOFFMAN_JACS(84)_106_2006
	huck_pars_std["CO_4P"] =  3.8/HARTREE_TO_EV;   // SAILLARD_HOFFMAN_JACS(84)_106_2006
	huck_pars_std["CO_3D"] =  9.7/HARTREE_TO_EV;   // SAILLARD_HOFFMAN_JACS(84)_106_2006
    huck_pars_std["NI_4S"] =  7.8/HARTREE_TO_EV;   // SAILLARD_HOFFMAN_JACS(84)_106_2006
	huck_pars_std["NI_4P"] =  3.7/HARTREE_TO_EV;   // SAILLARD_HOFFMAN_JACS(84)_106_2006
	huck_pars_std["NI_3D"] =  9.9/HARTREE_TO_EV;   // SAILLARD_HOFFMAN_JACS(84)_106_2006
    huck_pars_std["CU_4S"] =  11.4/HARTREE_TO_EV;   // HAY_HOFFMAN_JACS(75)_97_4884
	huck_pars_std["CU_4P"] =  6.06/HARTREE_TO_EV;   // HAY_HOFFMAN_JACS(75)_97_4884
	huck_pars_std["CU_3D"] =  14.0/HARTREE_TO_EV;   // HAY_HOFFMAN_JACS(75)_97_4884
    huck_pars_std["ZN_4S"] =  12.41/HARTREE_TO_EV;   // SILVESTRE_ISR_J_CHEM(83)_23_139
	huck_pars_std["ZN_4P"] =  6.53/HARTREE_TO_EV;    // SILVESTRE_ISR_J_CHEM(83)_23_139
	huck_pars_std["ZN_3D"] =  19.792/HARTREE_TO_EV;  // FROM VELA_JPC(88)_92_5688 - For completeness

    huck_pars_std["RU_5S"] =  7.73/HARTREE_TO_EV;    // TATSUMI_HOFFMANN_JACS(81)_103_3340
	huck_pars_std["RU_5P"] =  4.44/HARTREE_TO_EV;    // TATSUMI_HOFFMANN_JACS(81)_103_3340
	huck_pars_std["RU_4D"] =  11.23/HARTREE_TO_EV;   // TATSUMI_HOFFMANN_JACS(81)_103_3340

	return TRUE;
}

int HaQCMod::InitHuckParsVela()
{
	huck_pars_std["H_1S"]  = 12.565/HARTREE_TO_EV;  
	huck_pars_std["HE_1S"] = 60.0  /HARTREE_TO_EV; 

    huck_pars_std["LI_2S"] =  5.439/HARTREE_TO_EV;  
	huck_pars_std["LI_2P"] =  3.61 /HARTREE_TO_EV;  
    huck_pars_std["BE_2S"] =  9.663/HARTREE_TO_EV; 
	huck_pars_std["BE_2P"] =  6.072/HARTREE_TO_EV; 
    huck_pars_std["B_2S"]  = 14.440/HARTREE_TO_EV; 
	huck_pars_std["B_2P"]  =  8.438/HARTREE_TO_EV; 
    huck_pars_std["C_2S"]  = 19.654/HARTREE_TO_EV; 
	huck_pars_std["C_2P"]  = 11.129/HARTREE_TO_EV; 
    huck_pars_std["N_2S"]  = 25.366/HARTREE_TO_EV; 
	huck_pars_std["N_2P"]  = 13.900/HARTREE_TO_EV; 
    huck_pars_std["O_2S"]  = 31.600/HARTREE_TO_EV; 
	huck_pars_std["O_2P"]  = 16.776/HARTREE_TO_EV; 
    huck_pars_std["F_2S"]  = 38.373/HARTREE_TO_EV; 
	huck_pars_std["F_2P"]  = 19.769/HARTREE_TO_EV; 
    huck_pars_std["NE_2S"]  = 60.0/HARTREE_TO_EV; 
	huck_pars_std["NE_2P"]  = 60.0/HARTREE_TO_EV; 

    huck_pars_std["NA_3S"] =  5.331/HARTREE_TO_EV;   
	huck_pars_std["NA_3P"] =  3.178/HARTREE_TO_EV;   
    huck_pars_std["MG_3S"] =  8.201/HARTREE_TO_EV;  
	huck_pars_std["MG_3P"] =  4.641/HARTREE_TO_EV;  
    huck_pars_std["AL_3S"] = 11.794/HARTREE_TO_EV;   
	huck_pars_std["AL_3P"] =  5.976/HARTREE_TO_EV;   
    huck_pars_std["SI_3S"] = 15.370/HARTREE_TO_EV;   
	huck_pars_std["SI_3P"] =  7.990/HARTREE_TO_EV;    
    huck_pars_std["P_3S"]  = 19.013/HARTREE_TO_EV;   
	huck_pars_std["P_3P"]  = 10.011/HARTREE_TO_EV;   
    huck_pars_std["S_3S"]  = 22.761/HARTREE_TO_EV;    
	huck_pars_std["S_3P"]  = 12.081/HARTREE_TO_EV;    
    huck_pars_std["CL_3S"] = 26.636/HARTREE_TO_EV;   
	huck_pars_std["CL_3P"] = 14.217/HARTREE_TO_EV;   

    huck_pars_std["FE_4S"] =  7.659/HARTREE_TO_EV;  
	huck_pars_std["FE_4P"] =  4.104/HARTREE_TO_EV;   
	huck_pars_std["FE_3D"] = 10.822/HARTREE_TO_EV;   
    huck_pars_std["CO_4S"] =  7.897/HARTREE_TO_EV;   
	huck_pars_std["CO_4P"] =  4.148/HARTREE_TO_EV;   
	huck_pars_std["CO_3D"] = 11.625/HARTREE_TO_EV;   
    huck_pars_std["NI_4S"] =  8.125/HARTREE_TO_EV;   
	huck_pars_std["NI_4P"] =  4.184/HARTREE_TO_EV;   
	huck_pars_std["NI_3D"] = 12.404/HARTREE_TO_EV;   
    huck_pars_std["CU_4S"] =  8.345/HARTREE_TO_EV;   
	huck_pars_std["CU_4P"] =  4.216/HARTREE_TO_EV;   
//	huck_pars_std["CU_3D"] = 13.162/HARTREE_TO_EV;   
	huck_pars_std["CU_3D"] = 11.500/HARTREE_TO_EV; // Changed by Igor to move Cu orbitals to the MO gap region   
	
    huck_pars_std["ZN_4S"] = 10.141/HARTREE_TO_EV;   
	huck_pars_std["ZN_4P"] =  5.100/HARTREE_TO_EV;    
	huck_pars_std["ZN_3D"] = 19.792/HARTREE_TO_EV;  

    huck_pars_std["RU_5S"] =  7.293/HARTREE_TO_EV;    
	huck_pars_std["RU_5P"] =  4.022/HARTREE_TO_EV;    
	huck_pars_std["RU_4D"] = 11.533/HARTREE_TO_EV;   

	return TRUE;
}


int HaQCMod::InitHuckHam(HaMat_double& hmat, HaMat_double& ss, ArrayOrb3D& bas)
{
	int nb = bas.GetNBfunc();
	if(ss.num_cols() != nb || ss.num_rows() != nb)
	{
		HaBasisSet::CalcOvlpMat(&bas,&bas,ss);
	}

//	if(huck_pars_std.size() == 0) InitHuckParsStd();
	if(huck_pars_std.size() == 0) InitHuckParsVela();

	hmat.newsize(nb,nb);
	hmat = 0.0;

	std::string bas_type = bas.GetClassName();

	if( bas_type == "GauBasisSet" )
	{
		GauBasisSet* pgbas = (GauBasisSet*) &bas;
		int i;
		for(i = 0; i < nb; i++)
		{
			std::string lbl_orb;
			if( i < pgbas->bf_lbls.size()) lbl_orb = pgbas->bf_lbls[i] ;
			HaAtom* aptr = (HaAtom*) pgbas->GetHostPt(i);
			std::string at_symb = HaAtom::GetStdSymbolElem(aptr->GetElemNo());

			std::string key_str = at_symb + "_";
			if(lbl_orb.size() > 0) key_str += lbl_orb[0];
			if(lbl_orb.size() > 1) key_str += lbl_orb[1];
			boost::to_upper(key_str);

			int icount = huck_pars_std.count(key_str.c_str());
			if(icount > 0) 
			{
				hmat.SetVal_idx0(i,i,-huck_pars_std[key_str]);
			}
			else
			{
				hmat.SetVal_idx0(i,i,-2.0);
			}
		}
		for( i = 0; i < nb; i++)
		{
			double hii = hmat.GetVal_idx0(i,i);
			int j;
			for( j = 0; i < nb; i++)
			{
				if( i == j) continue;
				
				double hjj = hmat.GetVal_idx0(j,j);
				double sij = ss.GetVal_idx0(i,j);
				hmat.SetVal_idx0(i,j, 1.75*(hii+ hjj)*0.5*sij);
			}
		}
		return TRUE;
	}
	return FALSE;
}


bool HaQCMod::InitBasOvlp()
{
	if(AtBasis.GetNBfunc() != ovlp_mat.num_cols() ) 
	{
		GauBasisSet::CalcOvlpMat(&AtBasis,&AtBasis,ovlp_mat);
	}
	return true;
}

int HaQCMod::GetNumAlphaEl(int active_bas) const
{
	int nel=0;

	if(active_bas && ActBas != NULL)
	{
		std::string bas_type = ActBas->GetClassName();
		if( bas_type == "GauBasisSet")
		{
			nel = ((GauBasisSet*)ActBas)->GetNumElectr();
		}
	}
	else
	{
		nel = AtBasis.GetNumElectr();
	}
	nel = nel - charge;
	nel += (mult-1);
	nel = nel/2;
	return nel;
}

int HaQCMod::GetNumBetaEl(int active_bas) const
{
	int nel=0;

	std::string bas_type = ActBas->GetClassName();
	if(active_bas && ActBas != NULL)
	{
		if( bas_type == "GauBasisSet")
		{
			nel = ((GauBasisSet*)ActBas)->GetNumElectr();
		}
	}
	else
	{
		nel = AtBasis.GetNumElectr();
	}
	nel = nel - charge;
	nel -= (mult-1);
	nel = nel/2;
	return nel;
}


HaMat_double& HaQCMod::GetOvlpMat()
{
	if( AtBasis.GetNBfunc() != ovlp_mat.num_cols() ) InitBasOvlp();
	return ovlp_mat;
}

bool HaQCMod::Init1eDens(GauFile& gfile)
{
	return true;
}

bool HaQCMod::InitMOs(GauFile& gfile)
{
	if(MOene.size() != 0)
	{
		MOene.newsize(0);
	}

	int iunit,len,offset;
	void* pdata;
	int res;

    mol_type gmol;

	iunit= -GauFile::IO_mol;
	len = (int)sizeof(gmol) / ((int)sizeof(double));
	offset = 0;
	pdata = &gmol;
	gfile.fileio("read",iunit,len,pdata,offset);

    int nmo = gmol.nbasis;

	MO_coef.newsize(nmo,nmo);
	iunit = -GauFile::IO_ca;
	len = nmo*nmo;
	offset = 0;
	pdata = MO_coef.begin();

    res= gfile.fileio("read",iunit,len,pdata,offset);
	
	MOene.newsize(nmo);

	iunit= -GauFile::IO_emo;
    len= nmo;
    offset= 0;
    pdata= MOene.begin();
  
    res= gfile.fileio("read",iunit,len,pdata,offset);
    if(res != 0 )
		return false;
	
	return true;
}


static bool fill_float_vec(double* fmat, FILE* fp, const int nn)
// Fill float vector from fchk file
{
	int i,j, nl, ic;
	char buf[256];
	char* str;
	
	ic=0;
	nl= nn/5;

	for( i=0; i < nl ; i++)
	{	
		str=fgets(buf,256,fp);
		if(str == NULL)
		{
			cerr << " Error in fill_float_mat() " << endl;
			cerr << " Error in reading line from the file " << endl;
			return false;
		}
		istrstream is(buf);

		for(j =0; j < 5; j++)
		{
			is >> fmat[ic];
			ic++;
		}
	}
	
	int nr= nn%5;
	
	if(nr != 0)
	{
		str=fgets(buf,256,fp);
		if(str == NULL)
		{
			cerr << " Error in fill_float_mat() " << endl;
			cerr << " Error in reading line from the file " << endl;
			return false;
		}
		istrstream is(buf);
		
		for(j= 1; j <= nr; j++)
		{
			is >> fmat[ic];
			ic++;
		}
	}
	return true;
}

bool HaQCMod::LoadDataFromFChk(const char* fname)
{
	char buf[256];
	bool result;
	FILE* fp= fopen(fname,"r");
	if(fp == NULL)
	{
		cerr << " Error in HaQCMod::LoadDataFromFChk() " << endl;
		cerr << " Error opening formatted checkpoint file " << fname << endl;
		return false;
	}

	result= find_line_in_file(fp, "Number of atoms", buf, 256);
	if(!result)
	{
		cerr << " Error in HaQCMod::LoadDataFromFChk() " << endl;
		cerr << " Cannot find a number of Atoms in Fchk file " << endl;
		return false;
	}
	int na;
	sscanf(buf+44,"%d",&na);

	MolSet* pmset = GetMolSet();
	int nat_real = pmset->GetNAtoms() - pmset->GetNDumAtoms();

	if( na != nat_real )
	{
		cerr << " Error in HaQCMod::LoadDataFromFChk() " << endl;
		cerr << " Number of atoms in the host molecule " << GetMolSet()->GetNAtoms() << endl;
		cerr << " Is not equal to the number of atoms in fchk file " << na << endl;
		return false;
	}

	result= find_line_in_file(fp, "Number of basis functions", buf, 256);
	if(!result)
	{
		cerr << " Error in HaQCMod::LoadDataFromFChk() " << endl;
		cerr << " Cannot find a number of Basis Functions in Fchk file " << endl;
		return false;
	}

	int nb ;
	sscanf(buf+44,"%d",&nb);
	if( nb != GetNBfunc() )
	{
		cerr << " Error in HaQCMod::LoadDataFromFChk() " << endl;
		cerr << " The number of basis functions " << GetNBfunc() << endl;
		cerr << " Is not equal to the number of basis functions in the fchk file " << nb << endl;
		return false;
	}

	int num_mo=nb;
	result= find_line_in_file(fp, "Number of independent functions", buf, 256);
	if(result)
	{
		sscanf(buf+44,"%d",&num_mo);
	}
	else
	{
		cerr << " Warning in HaQCMod::LoadDataFromFChk() " << endl;
		cerr << " Cannot find the number of independent basis functions in Fchk file " << endl;
		cerr << " The number of MOs is assumed to be equal to the number of basis functions " << endl;
		rewind(fp);
	}

	if(load_mo_flag)
	{
		result= find_line_in_file(fp, "Alpha Orbital Energies", buf, 256);
		if(!result)
		{
			cerr << " Error in HaQCMod::LoadDataFromFChk() " << endl;
			cerr << " Cannot find Molecular Orbital energies in Fchk file " << endl;
			return false;
		}
		
		MOene.newsize(num_mo);
		
		fill_float_vec(MOene.begin(), fp, num_mo);

		result= find_line_in_file(fp, "Alpha MO coefficients", buf, 256);
		if(!result)
		{
			cerr << " Error in HaQCMod::LoadDataFromFChk() " << endl;
			cerr << " Cannot find Molecular Orbital Coeffcients in Fchk file " << endl;
			return false;
		}

		MO_coef.newsize(nb,num_mo);
		fill_float_vec(MO_coef.begin(), fp, nb*num_mo);
	}	
	return true;

}

int
HaQCMod::PrepGauss()
{
#if defined(GAUSSVER) 
	maxmem_.maxm[0] = 32000000;
	maxmem_.maxm[1] = 32000000;	
#endif
	return TRUE;
}

double 
calc_chem_pot(double* ene_mo, int nmo, int nel, double beta)
{
	double mu_1 = ene_mo[0];
	double mu_2 = ene_mo[nmo-1];

	double mu;
	int i,j;
	for( i = 0; i < 100; i++)
	{
	   mu = (mu_1 + mu_2)/2.0;
	   double tot_el = 0.0;
	   for(j = 0; j < nmo; j++)
	   {
		  tot_el +=1.0/(exp(beta*(ene_mo[j] - mu)) + 1);
	   }
	   double delt_nel = tot_el - (double)nel;
	   if( fabs(delt_nel) < 1.0e-5) break;
       if( delt_nel > 0.0) 
	   {
		   mu_2 = mu;
	   }
	   else
	   {
           mu_1 = mu;
	   }
	}
	return mu;
}

int  HaQCMod::FormDenMat(const HaMat_double& cmo, double* pa, int nel, double* pene_mo, double temp_fermi)
{
//
// Can probably be more effective if first transpose cmo and the multiply
// making dot-product of cols of transpose matrix, see form2 of Gaussian
//
// Assume lower triangle storage of pa ( should change to symmetric matrix when realized)
//
	int nb = cmo.num_rows();
	int nmo = cmo.num_cols();
	int i,j,k;
	int ntt = nb*(nb+1)/2;
	for(i=0; i < ntt; i++)
	{
		pa[i] = 0.0;
	}
    
	double mu = 0.0;
	HaVec_double orb_pop;
	if(temp_fermi > 1.0e-5 && pene_mo != NULL)
	{
		double beta = HARTREE_TO_KT*298.0/temp_fermi;
		mu = calc_chem_pot(pene_mo,nmo,nel,beta);
		orb_pop.newsize(nmo);
		for(j = 0; j < nmo; j++)
		{
            orb_pop[j] = 1.0/(exp(beta*(pene_mo[j] - mu)) + 1);
		}
	}

	int idx = 0;
	for(i=0; i < nb; i++)
	{
		for(j=0; j <= i; j++)
		{
			if( orb_pop.size() == nmo)
			{
			   for(k= 0; k < nmo; k++)
			   {
				  pa[idx] += orb_pop[k]*cmo.r0(i,k)* cmo.r0(j,k);
			   }               
			}
			else
			{
			   for(k= 0; k < nel; k++)
			   {
				  pa[idx] += cmo.r0(i,k)* cmo.r0(j,k);
			   }
			}
			idx++;
		}
	}
	return TRUE;
}	

class CNDOThread: public wxThread
{
public:
	CNDOThread(HaQCMod* ptr_qc_mod_new) { ptr_qc_mod = ptr_qc_mod_new; }
	virtual ExitCode Entry()
	{
		ptr_qc_mod->RunCNDOThread();
		return 0;
	}
	
	HaQCMod* ptr_qc_mod;
};

int HaQCMod::RunCNDO()
{
	wxThread* cndo_run_thread = new CNDOThread(this);
	cndo_run_thread->Create();
	cndo_run_thread->Run();
	return TRUE;
}

int HaQCMod::Run(const RunOptions* popt_par )
{
	const RunOptions* popt = popt_par;
	if( popt == NULL ) popt = &run_opt_default;

	MolSet* pmset = GetMolSet();
	HaGaussMod* p_gauss = pmset->GetGaussMod(true);
	p_gauss->SaveInpFile();
	p_gauss->Run();
	return TRUE;
}

int HaQCMod::StopCalc()
{
	stop_calc_flag = TRUE;
	return True;
}

void HaQCMod::SetEne( double ene_par )
{
	ene = ene_par;
}

void HaQCMod::SetHFEne( double ene_hf_par )
{	
	ene_hf = ene_hf_par;
}

void HaQCMod::SetDFTEne( double ene_dft_par )
{
	ene_dft = ene_dft_par;
}

double HaQCMod::GetEne() const
{
	return ene;
}

double HaQCMod::GetHFEne() const
{
	return ene_hf;
}

double HaQCMod::GetDFTEne() const
{
	return ene_dft;
}

int HaQCMod::SetExtCharge( HaAtom* aptr, double ch)
{
	if( aptr == NULL ) return FALSE;
	ext_chrg[aptr] = ch;
	return TRUE;
}

int HaQCMod::SetExtCharge( const std::string& at_ref, double ch)
{
	HaAtom* aptr = GetMolSet()->GetAtomByRef( at_ref.c_str());
	if( aptr == NULL ) return FALSE;
	return SetExtCharge(aptr,ch);
}

void HaQCMod::SetExtChCrdOffset( double crd_offset )
{
	ext_ch_crd_offset = crd_offset;
}

#ifdef USE_IPACK

int HaQCMod::InitIPack()
{
	if(timer == NULL)
	{
		timer = new Timer;
		init_timer();
		timer -> start(t_main);
		timer -> start(t_sp);
		init_ftable();
	}

	if( pm == NULL)
	{
		int iargc = 0; 
		char* arg = "";
		char** argv_p = &arg;
		
		pm = new ParMachine(&iargc, &argv_p);
	}

	return TRUE;
}

int HaQCMod::TestIPack1() 
{
	InitIPack();

	InternalBasisList iblist;

	InternalBasis* pibas = AtBasis.CreateIPackBas();
	
	iblist.add(pibas);

	int same = True;
	int nb = pibas->no_function();

	ARRAY<Mat> ovlp_arr;

	overlap(iblist  ,ovlp_arr);
	normalize(iblist,ovlp_arr);
	ovlp_arr[0].write("SMATRIX");

//	ARRAY<Mat> ovlp_arr(1);
//	ovlp_arr[0].reset(nb,nb);

//	compute_overlap(*pibas,*pibas, ovlp_arr, same);

	cout << ovlp_arr[0];

	timer -> stop(t_main);
	timer -> print();
//    delete timer;
//    timer = 0;

	return True;
}

int HaQCMod::TestRandomGen()
{
	Random rand_num_gen(5);

	int i;
	
	ofstream frand("Random_num.dat");

	for(i = 0; i < 10000; i++)
	{
		frand << rand_num_gen() << "  ";
		rand_num_gen();
		frand << rand_num_gen() << "  " << endl;
		rand_num_gen();
	}

	return TRUE;
}

int HaQCMod::TestIPack2()
{
	InitIPack();

    InternalBasisList iblist;

	InternalBasis* pibas = AtBasis.CreateIPackBas();

	iblist.add(pibas);

	int same = True;
	int nb = pibas->no_function();

	ARRAY<Mat> nuc_attr_arr;
	Molecule ipack_mol;

	MolSet* pmset = GetMolSet();

	AtomIteratorMolSet aitr(pmset);
	HaAtom* aptr;

	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		int elem = aptr->GetElemNo();
		Location at_loc( aptr->GetX(), aptr->GetY(), aptr->GetZ() );
		ipack_mol.add(elem,at_loc,"test");
	}

	ARRAY<Mat> ovlp_arr;
	overlap(iblist  ,ovlp_arr);

	nuclear(iblist,ipack_mol,nuc_attr_arr);
	normalize(iblist,nuc_attr_arr);

	nuc_attr_arr[0].write("NUC_ATTR_INT");

	cout << nuc_attr_arr[0];

	timer -> stop(t_main);
	timer -> print();
    delete timer;
    timer = 0;

	return True;
}

#endif