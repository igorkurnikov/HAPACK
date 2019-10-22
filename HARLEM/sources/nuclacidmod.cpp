/*!  \file nuclacidmod.cpp

    Classes to model Nucleic Acid program from HARLEM  

    \author Igor Kurnikov 
    \date  2002

*/

#define NUCLACIDMOD_CPP

#include "haio.h"

#include <boost/algorithm/string.hpp>

#include "nuclacidmod.h"
#include "hamolecule.h"
#include "moleditor.h"
#include "haresdb.h"
#include "hamolview.h"
#include "harlemapp.h"
#include "math.h"

NuclAcidMod::NuclAcidMod(MolSet* new_phost_mset):
HaCompMod(COMP_MOD_NUCL_ACID,new_phost_mset)
{
	SetStdParams();
}

NuclAcidMod::~NuclAcidMod()
{
 
}

int NuclAcidMod::SetStdParams()
{	
	p_dna_mol = FindDNAMol();

	ChainIteratorMolecule chitr(p_dna_mol);
	HaChain* ch1;
	HaChain* ch2;

//	seq[0] = "GCGCGCGCGCGCGCGC";
//	seq[1] = "CGCGCGCGCGCGCGCG";
	seq[0] = "";
	seq[1] = "";
	seq[2] = "";
	seq[3] = "";

	if(p_dna_mol == NULL)
	{
		seq[0] = "GC";
		seq[1] = "CG";
	}
	else
	{
		ch1 = chitr.GetFirstChain();
		ch2 = chitr.GetNextChain();
		HaResidue* pres;
		std::string tmp;
		ResidueIteratorChain ritr_ch_1(ch1);
		for(pres = ritr_ch_1.GetFirstRes(); pres; pres = ritr_ch_1.GetNextRes())
		{
			tmp = pres->GetName();
			seq[0] += tmp[0];
		}

		ResidueIteratorChain ritr_ch_2(ch2);
		for(pres = ritr_ch_2.GetFirstRes(); pres; pres = ritr_ch_2.GetNextRes())
		{
			tmp = pres->GetName();
			seq[1] += tmp[0];
		}
	}
 
	nsym_unit.newsize(4);
	nsym_unit = 0;
	homon_symm_flags.newsize(3);
    homon_symm_flags = FALSE;
	homon_symm_offs.newsize(3);
	homon_symm_offs = -1;
	 	
	nbreak_1 = 0;
	nbreak_2 = 0;

	dir_strand.newsize(4);

	dir_strand[0] = 1;
	dir_strand[1] = -1;
	dir_strand[2] = 1;
	dir_strand[3] = -1;

	init_flag = FALSE;
	tot_energy = 0.0;

	sup_helix_rad = 0.0;
	sup_helix_pit = 0.0;
	fst_twist = 0.0;

	update_var_flag = FALSE;

	ene_per_unit_flag = FALSE;

	max_iter = 2000;
	SetFFtype("Flex");
	SetDielSlope(0.16);
	SetPhosChrg(-0.5);
	std::string nucl_dir = pApp->harlem_home_dir + "residues_db/";

#if defined(INT_JUMNA)
	strcpy_to_fort(&jmdirs_.libn[0], nucl_dir.c_str(),80);

	int i;
	for(i=0; i< N6A_J; i++)
	{
		mnn_.lock[i] = 0;
	}
#endif

	return TRUE;
}


int
NuclAcidMod::SetSeq(const char* seq_str)
{
	seq[0] = seq_str;
    int res = GenComplStrand();
	return res;
}

int NuclAcidMod::GenComplStrand()
{
	boost::trim(seq[0]);
	boost::to_upper(seq[0]);
	int len = seq[0].size();
	int i;
	seq[1] = "";
	for(i=0; i < len; i++)
	{
		if(seq[0][i] == 'A') seq[1] += 'T';
		else if(seq[0][i] == 'T') seq[1] += 'A';
		else if(seq[0][i] == 'G') seq[1] += 'C';
		else if(seq[0][i] == 'C') seq[1] += 'G';
		else
		{
			PrintLog(" Error NuclAcidMod::GenComplStrand() \n");
			PrintLog(" Unknown nucleotide symbol: %c \n",seq[0][i]);
			return FALSE;
		}
	}
	return TRUE;	
}

int NuclAcidMod::SetFFtype(const char* new_ff_type)
{
	int invalid = TRUE;
	std::string ff_name = new_ff_type;
	boost::trim(ff_name);
	boost::to_upper(ff_name);
	if(ff_name == "FLEX")
	{
		invalid = FALSE;
	    force_field = FLEX_FF;
#if defined(INT_JUMNA)
		strcpy_to_fort(&cha_.parm[0],"Flex",32);
#endif
	}
	else if(ff_name == "AMBER91")
	{
		invalid = FALSE;
		force_field = AMBER91_FF;
#if defined(INT_JUMNA)
		strcpy_to_fort(&cha_.parm[0],"Amber91",32);
#endif
	}
	else if(ff_name == "AMBER94")
	{
		invalid = FALSE;
		force_field = AMBER94_FF;
#if defined(INT_JUMNA)
		strcpy_to_fort(&cha_.parm[0],"Amber94",32);
#endif
	}

	if(!invalid) return TRUE;

	PrintLog(" Error in NuclAcidMod::SetFFtype() \n");
	PrintLog(" Invalid Force Field Type %s \n",new_ff_type);
	return FALSE;

}

int 
NuclAcidMod::SetFFtypeIdx(const int i_ff_type)
{
	if(i_ff_type == 1)
		SetFFtype("FLEX");
	else if(i_ff_type == 2)
		SetFFtype("AMBER91");
	else if(i_ff_type == 2)
		SetFFtype("AMBER94");
	else 
		SetFFtype("FLEX");

	return TRUE;
}

int
NuclAcidMod::SetDielSlope( double slope_new)
{
	diel_slope = slope_new;
#if defined(INT_JUMNA)
	datjm_.slope = slope_new;
#endif
	return TRUE;
}

int
NuclAcidMod::SetPhosChrg( double phos_chrg_new)
{
	phos_chrg = phos_chrg_new;
#if defined(INT_JUMNA)
	datjm_.phos = phos_chrg_new;
#endif
	return TRUE;
}

int NuclAcidMod::SetHelCoord(int i_strand, int i_res, int i_crd, double crd_val)
{
	if(p_dna_mol == NULL) return FALSE;
	int istr;
    
	int nres = 0;
	ChainIteratorMolecule chitr(p_dna_mol);
	HaChain* chain = chitr.GetFirstChain();
	for( istr = 1; istr < i_strand; istr++)
	{
		if(chain == NULL) return FALSE;
		nres += chain->GetNRes();
		chain = chitr.GetNextChain();
	}
	if(chain == NULL) return FALSE;
	
	nres += i_res;

	if(i_crd > 6) return FALSE;
    hel_crd(6*(nres-1) + i_crd) = crd_val; 
	
	return TRUE;
}

int NuclAcidMod::SetBBCoord(int i_strand, int i_res, int i_crd, double crd_val)
{

	return TRUE;
}

int 
NuclAcidMod::LockBBCoord(int i_strand, int i_res, int i_crd, int do_lock)
{

	return TRUE;	
}


int 
NuclAcidMod::SetShlxTwist(double shlx_tw)
{
	if(p_dna_mol == NULL) return FALSE;
    if(fabs(sup_helix_rad) < 0.000001)
	{
		PrintLog(" Error in SetShlxTwist()\n");
		PrintLog(" SuperHelix Radius is not set \n");
		return FALSE;
	}
//	double old_tw = hel_crd(6);
//	double delt_tw = shlx_tw - old_tw;
	SetHelCoord(1, 1, 6, shlx_tw);
	SetHelCoord(2, 1, 6, shlx_tw);
	
//	int nbbv   = NumBBCoord();
//	int nbbv_f = NumFreeBBCoord();

//	int i;
//	int k = nbbv;
//	int ki = nbbv_f;
//	int locked;

//	int ir_abs = 0;
//	int istrand = 0;
//	HaChain* chain = p_dna_mol->GetFirstChain();
//	if(chain == NULL) return FALSE;
//	for(; chain; chain = p_dna_mol->GetNextChain())
//	{
//		istrand++;
//		int nres = chain->GetNRes();
//		int ir;
		
//		for(ir = 1; ir <= nres; ir++)
//		{	
//			ir_abs++;
//			for(i=0; i < 5; i++)
//			{
//				locked = mnn_.lock[k];
//				if(ir == 1 && i == 2 && istrand == 1) continue;
//				if(!locked) ki++;
//			}
//			double tw = strjm_.hel[5*N2_J + (ir_abs - 1)] + delt_tw;
//			if( tw > 360.0) tw= tw - 360.0;
//			
//			strjm_.hel[5*N2_J + (ir_abs - 1)] = tw;
//			hel_crd(6*ir) = tw;
//			locked = mnn_.lock[k];
//			if(!locked) 
//			{
//				mnn_.var[ki] = tw;
//				ki++;
//			}
//		}
//	}
	
	return TRUE;
}

int NuclAcidMod::LockHelCoord(int i_strand, int i_res, int i_crd, int do_lock )
//! i_strand - strand index (1-based)
//! i_res    - index of residue in the strand( 1-based)
//! i_crd    - helical coordinate index (1-based)
{
	if(p_dna_mol == NULL) return FALSE;
	int istr;
    
	int nres = 0;
	ChainIteratorMolecule chitr(p_dna_mol);
	HaChain* chain = chitr.GetFirstChain();
	for( istr = 1; istr < i_strand; istr++)
	{
		if(chain == NULL) return FALSE;
		nres += chain->GetNRes();
		chain = chitr.GetNextChain();
	}
	if(chain == NULL) return FALSE;
	
	nres += i_res;

	int ioff = NumBBCoord();

	if(do_lock)
	{
		lock_hel[ (nres - 1)*6 + (i_crd-1) ] = 1;
#if defined(INT_JUMNA)
		mnn_.lock[ioff + (nres - 1)*6 + i_crd - 1] = 1;
#endif
	}
	else
	{
		lock_hel[ (nres - 1)*6 + (i_crd-1) ] = 0;
#if defined(INT_JUMNA)
		mnn_.lock[ioff + (nres - 1)*6 + i_crd - 1] = 0;
#endif
	}

	update_var_flag = TRUE;

//	strcpy_to_fort(&strjm_.code[(nres-1)*8],jumna_code,8);

	return TRUE;
}

int 
NuclAcidMod::IsHelCoordLocked(int ir, int i_crd)
{
	return lock_hel[ir*6 + i_crd];
}

int 
NuclAcidMod::BuildNuclAcid()
{
#if defined(INT_JUMNA)
//	strcpy_to_fort(jfinp_.finp,"gc_build_3.inp",80);
//  put to default:
	strcpy_to_fort(&cha_.out[0],"jumnatst",32);
	strcpy_to_fort(&cha_.pdb[0],"jumnatst",32);

	int nstrands = 0;     // Number of strands in the system
	int len_strand[4] = {0,0,0,0};  // Number of nucleotides in the strands  (positive means 5'->3' direction, negative 3'->5') 

	int is;
	for(is = 1; is < 4; is++)
	{
		boost::trim(seq[is]); 
		boost::to_upper(seq[is]); 
	}

	int nlen = seq[0].size();
	len_strand[0] = nlen;
	if(nlen == 0 || nlen > 40 )
	{
		ErrorInMod("NuclAcidMod::BuildNuclAcid()",
			       "Length of the nucleotide sequence should be > 0 and <= 40 \n");
		return FALSE;
	}
    nstrands = 1;

	for(is = 1; is < 4; is++)
	{
		if( seq[is].size() > nlen)
		{
			ErrorInMod("NuclAcidMod::BuildNuclAcid()",
				" First chain should be the longest");
			return FALSE;
		}
		if(seq[is].size() == 0) break;
		
		nstrands++;
		len_strand[is] = seq[is].size()*dir_strand[is];
		int ndiff = seq[is].size() - nlen;
		for(int j=0; j < ndiff; j++)
			seq[is] += "-";
	}

	std::string tot_seq;
	for(is = 0; is < nstrands; is++)
	{
		tot_seq += seq[is];
	}

	strcpy_to_fort(strjm_.seq,tot_seq.c_str(),120);
    strjm_.nst = nstrands;

	int i,ir;
	for(i=0; i < 4; i++)
	{
		strjm_.idr[i] = len_strand[i];
	}
	symjm_.ksym[0] = nsym_unit[0];
    symjm_.ksym[1] = nsym_unit[1];
    symjm_.ksym[2] = nsym_unit[2];
    symjm_.nbrk[0]  = nbreak_1;
    symjm_.nbrk[1]  = nbreak_2;

	datjm_.ecen   = ene_per_unit_flag;
	datjm_.homo   = homon_symm_flags[0];
	datjm_.homo2  = homon_symm_flags[1];
	datjm_.homo3  = homon_symm_flags[2];
	datjm_.mhomo  = homon_symm_offs[0];
	datjm_.mhomo2 = homon_symm_offs[1];
	datjm_.mhomo3 = homon_symm_offs[2];
	
	datjm_.rad   = sup_helix_rad;
	datjm_.pit   = sup_helix_pit;

	int j;
// Set default B-conformation:
	int ir_abs = -1;
	for(is=0; is < 4; is++)
	{
		for(ir = 0; ir < abs(len_strand[is]);ir++)
		{
			ir_abs++;
			for(j = 0; j < 6; j++)
			{
				slg_.hst[ j*N2_J + ir_abs] = 0;
			}
			for(j = 0; j < N8_J; j++)
			{
				slg_.bst[ j*N2_J + ir_abs] = 0;
			}
		
			double rise; 
			double twist; 

			rise = 3.4;

			strcpy_to_fort(&strjm_.code[ir_abs*8],"--------",8);
			twist = 36.0;
		
			if( ir == 0 )
			{
				rise = 0.0;
				twist = 0.0;
			}
			strjm_.hel[ 2*N2_J + ir_abs]= rise;
			slg_.hst[ 2*N2_J + ir_abs] = 1;
            strjm_.hel[ 5*N2_J + ir_abs]= twist;
			slg_.hst[ 2*N2_J + ir_abs] = 1;
		}
	}
	
	datjm_.opt  = 0;
	jumcall_();
	
//	datjm_.homo = TRUE;
//	datjm_.ecen = TRUE;

	int nr = ir_abs + 1;

	hel_crd.newsize(6*nr);
	lock_hel.newsize(6*nr);
	lock_hel = 0;

	bb_crd.newsize(12*nr);
	lock_bb.newsize(12*nr);
	lock_bb = 0;

	CreateMolFromJumna();
	if(p_dna_mol == NULL) return FALSE;

    LockHelCoord(1,1,3,1);
	if(!IsSupHlxConstr())LockHelCoord(1,1,6,1);

	int nr1 = p_dna_mol->GetNRes();

	init_flag = TRUE;
#endif
	return TRUE;
}

int 
NuclAcidMod::UpdateXYZ()
{
#if defined(INT_JUMNA)
	microb_();
	SetCoordsFromJumna();
#endif
	return TRUE;
}

int 
NuclAcidMod::MinEne()
{
#if defined(INT_JUMNA)
	int_4 ipl[N9_J];
	int_4 iabas = 0;
	int i;
	for(i=0; i< N9_J; i++)
		ipl[i] = 0;

	if(!init_flag)
	{
		BuildNuclAcid();
	}

	datjm_.opt = 1;
	datjm_.maxn = max_iter;
	minim_();
	helix_();
	backbo_();
	penalty_();
	disth_();

	if(force_field == FLEX_FF)
	{
		ecomp_();
	}
	else if( force_field == AMBER91_FF)
	{
		ecomp91_();
	}
	else if( force_field == AMBER94_FF)
	{
		ecomp94_();
	}

	tot_energy = enf_.ener - enf_.epen;

	PrintLog(" Tot ene=%9.3f(wo/pen),%9.3f(w/pen) penalty= %9.3f \n",
		       tot_energy, enf_.ener, enf_.epen);
	PrintLog(" Vdw = %9.3f  Electr= %9.3f \n",
		     (enf_.repl + enf_.disp),enf_.elec);
	PrintLog(" VAng ene = %9.3f  Torsional = %9.3f \n",
		       enf_.eang, enf_.etog);

#endif
	return TRUE;
}

int 
NuclAcidMod::CalcEne()
{
#if defined(INT_JUMNA) 

	int_4 ipl[N9_J];
	int_4 iabas = 0;
	int i;
	for(i=0; i< N9_J; i++)
		ipl[i] = 0;

    datjm_.opt = -1;
	if(!init_flag)
	{
		BuildNuclAcid();
	}
	
	microb_();
	helix_();
	backbo_();
	penalty_();
	if(force_field == FLEX_FF)
	{
		ecomp_();
	}
	else if( force_field == AMBER91_FF)
	{
		ecomp91_();
	}
	else if( force_field == AMBER94_FF)
	{
		ecomp94_();
	}
	//		closejm_(&ipl[0],&iabas);

	PrintLog(" Tot ene=%9.3f(w/pen),%9.3f(wo/pen) penalty= %9.3f \n",
		       enf_.ener, (enf_.ener - enf_.epen), enf_.epen);
	PrintLog(" Vdw = %9.3f  Electr= %9.3f \n",
		     (enf_.repl + enf_.disp),enf_.elec);
	PrintLog(" VAng ene = %9.3f  Torsional = %9.3f \n",
		       enf_.eang, enf_.etog);

#endif

	return TRUE;
}

int
NuclAcidMod::UpdateVarCoord()
{
#if defined(INT_JUMNA)
    setvar_();
    int_4 nspc,nbcc,nhlc,nklc,nlgc;
	equivjm_(&nspc,&nbcc,&nhlc,&nklc,&nlgc);
	update_var_flag = FALSE;
#endif
	return TRUE;
}


int NuclAcidMod::SetCoordsFromJumna()
{
#if defined(INT_JUMNA)

    int na = mrc_.kam;
	MolSet* pmset = GetMolSet();
	HaMolecule* pMol= pmset->GetMolByName("DNA");
	if(pMol == NULL)
	{
		ErrorInMod("NuclAcidMod::SetCoordsFromJumna()",
			       " DNA Molecule is not set");
		return FALSE;
	}
	if(pMol->GetNAtoms() != na)
	{
		PrintLog("Error in NuclAcidMod::SetCoordsFromJumna()\n");
		PrintLog("The number of atoms in the DNA molecule %d \n do not equal atoms in JUMNA %d \n",
			      pMol->GetNAtoms(),na);
		return FALSE;
	}

	AtomIteratorMolecule aitr(pMol);
	HaAtom* aptr;

	int ia = -1;

	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		ia++;
		aptr->SetX(mrc_.corm[ia]);
		aptr->SetY(mrc_.corm[N1_J+ia]);
        aptr->SetZ(mrc_.corm[N1_J*2+ia]);
	}

    int nr = pMol->GetNRes();
	int nrj = strjm_.nto;

	if( nr != nrj)
	{
		PrintLog("Error in NuclAcidMod::SetCoordsFromJumna()\n");
		PrintLog("The number of residues in the DNA molecule %d \n do not equal to the number residues in JUMNA %d \n",
			      nr,nrj);
		return FALSE;
	}

	if(hel_crd.size() != nr * 6 )
	{
		hel_crd.newsize(nr * 6);
	}

	int ir;
    
	ChainIteratorMolecule chitr(p_dna_mol);
    HaChain* chain1 = chitr.GetFirstChain();

	int nr1 = chain1->GetNRes();
	int iv;

	for(ir = 0; ir < nr; ir++)  // Set Backbone variables
	{
		int niv = NumIndBBVarRes(ir);
		for( iv = 0; iv < niv; iv++)
		{
			int koff = BBOffset(ir);
			bb_crd[12*ir + iv] = flx_.sap[koff + iv];
		}
	}

	int nt = NumBBCoord() - 1;	
	for(ir = 0; ir < nr; ir++)  // Set Thymine methyl torsion
	{
		bb_crd[12*ir + 11] = 0.0;
		if( IsThymineRes(ir))
		{
			nt++;
			bb_crd[12*ir + 11] = flx_.sap[nt];
		}
	}

    double zdel,wdel,zglo,wglo;

	for(ir = 0; ir < nr; ir++)  // Set Helical variables
	{
		if( ir < nr1 )
		{
			double sign = 1.0;
			zglo = strjm_.hel[2*N2_J + ir];
			wglo = strjm_.hel[5*N2_J + ir];
			if(IsSupHlxConstr() && ir > 0) wglo = strjm_.hel[5*N2_J + ir] - strjm_.hel[5*N2_J + ir - 1];
			if(wglo >  180.0) wglo=wglo- 360.0;
			if(wglo < -180.0) wglo=wglo + 360.0;
			hel_crd[6*ir] = -1.0*strjm_.hel[0*N2_J +ir];	
			hel_crd[6*ir+1] = strjm_.hel[1*N2_J +ir];	
			hel_crd[6*ir+2] = zglo;
			hel_crd[6*ir+3] = strjm_.hel[3*N2_J +ir];	
			hel_crd[6*ir+4] = -1.0*strjm_.hel[4*N2_J +ir];
			hel_crd[6*ir+5] = wglo;
		}
		else
		{
			int i1 = IdxResInChain(ir);

			if( IsFstResInChain(ir))
			{
				zdel = strjm_.hel[2*N2_J + ir];
				wdel = strjm_.hel[5*N2_J + ir];
				zglo = strjm_.hel[2*N2_J + i1] + zdel;
				wglo = strjm_.hel[5*N2_J + i1] + wdel;
				if(IsSupHlxConstr()) wglo = wdel;
			}
			else
			{
				if(!symjm_.ihl[ir] )
				{
					zglo = strjm_.hel[2*N2_J + ir];
					wglo = strjm_.hel[5*N2_J + ir];
					zdel = zdel + zglo - strjm_.hel[2*N2_J + i1];
					wdel = wdel + wglo - strjm_.hel[5*N2_J + i1];
				}
				else
				{
					zglo = (strjm_.hel[2*N2_J + i1] + strjm_.hel[2*N2_J + ir]) - zdel;
					wglo = (strjm_.hel[5*N2_J + i1] + strjm_.hel[5*N2_J + ir]) - wdel;
					zdel = strjm_.hel[2*N2_J + ir];
					wdel = strjm_.hel[5*N2_J + ir];
				}
				if(IsSupHlxConstr()) wglo = strjm_.hel[5*N2_J + ir] - strjm_.hel[5*N2_J + ir - 1];
			}
			if(wglo >  180.0) wglo=wglo-360.0;
			if(wglo < -180.0) wglo=wglo + 360.0;
			hel_crd[6*ir] = -1.0*strjm_.hel[0*N2_J +ir];	
			hel_crd[6*ir+1] = strjm_.hel[1*N2_J +ir];	
			hel_crd[6*ir+2] = zglo;
			hel_crd[6*ir+3] = strjm_.hel[3*N2_J +ir];	
			hel_crd[6*ir+4] = -1.0*strjm_.hel[4*N2_J +ir];
			hel_crd[6*ir+5] = wglo;
		}
	}

#endif
	return TRUE;
}


int 
NuclAcidMod::CalcLocCrdSys()
{
	if(p_dna_mol == NULL) 
	{
		HaMolecule* pmol = FindDNAMol();
		if(pmol == NULL) return FALSE;
	}

	HaMolecule::ResidueIterator ritr(p_dna_mol);

	HaResidue* pres;
	for(pres = ritr.GetFirstRes(); pres; pres = ritr.GetNextRes())
	{
		pres->CalcStdCrdSys();
	}
	return TRUE;
}

int NuclAcidMod::CalcAxisPar1(HaMat_double& ax_hlxc, double& sum, HaVec_double& gra, HaVec_double& scp)
{
	if(p_dna_mol == NULL) return FALSE;

	CalcLocCrdSys();

	HaChain* ch;
	ChainIteratorMolecule chitr(p_dna_mol);
	HaChain* ch1 = chitr.GetFirstChain();
	HaChain* ch2 = chitr.GetNextChain();

	if( ch1 == NULL || ch2 == NULL) return FALSE;

	int nr = ch1->GetNRes();

	Vec3DValArray b,s,up,us,um,q;

	u.resize(nr);
	h.resize(nr);
	o.resize(nr);
	b.resize(nr*4);
	s.resize(nr);
	up.resize(nr);
	us.resize(nr);
	um.resize(nr);
	q.resize(nr);

	// from /der/
	Vec3DValArray ud,pd;
	ud.resize(nr*4);
	pd.resize(nr*4);
	
	int ir = -1;
	int j;
	ResidueIteratorChain ritr_ch_1(ch1);
	ResidueIteratorChain ritr_ch_2(ch2);

	HaResidue* pres1 = ritr_ch_1.GetFirstRes();
	HaResidue* pres2 = ritr_ch_2.GetFirstRes();
	 
	for(; pres1 != NULL; pres1 = ritr_ch_1.GetNextRes(), pres2 = ritr_ch_2.GetNextRes() )
	{
		ir++;
		double xdi = ax_hlxc(1,ir+1);
		double ydi = ax_hlxc(2,ir+1);
		double cln = ax_hlxc(3,ir+1);
		double tip = ax_hlxc(4,ir+1);
		double ca  = cos(DEG_TO_RAD*(-tip));
        double sa  = sin(DEG_TO_RAD*(-tip));

		Vec3D r,v,w,t,d;

		r = pres1->std_crd_sys[1];
		v = pres1->std_crd_sys[2];
		t = v;
		t.Rotate(r,ca,sa);
		Vec3D::VecProduct(d,r,t);
		ca  = cos(DEG_TO_RAD*(-cln));
		sa  = sin(DEG_TO_RAD*(-cln));

		r = d;
		v = t;
		u[ir] = v;
		u[ir].Rotate(r,ca,sa);
		
		Vec3D::VecProduct(w,u[ir],d);

		h[ir]= pres1->std_crd_sys[3] - xdi*d - ydi*w;
		o[ir] = h[ir];

		double dot = 0.0;
		Vec3D::diff(v,pres2->std_crd_sys[3],h[ir]);

		dot += Vec3D::DotProduct(u[ir],v);
		dot = dot/2;
		
		o[ir] = h[ir] + dot*u[ir];

		Vec3D::diff(b[ir],o[ir],pres1->std_crd_sys[3]);
		Vec3D::diff(b[nr+ir],o[ir],pres2->std_crd_sys[3]);
        	
		if(ir > 0) s[ir] = o[ir] - o[ir-1];
		
//		PrintLog( "for res %d \n", (ir+1));
//		PrintLog( "s = (%9.3f,%9.3f,%9.3f) , o =(%9.3f,%9.3f%9.3f) \n",
//			      s[ir][0],s[ir][1],s[ir][2],o[ir][0],o[ir][1],o[ir][2]);
	}

	HaMat_double qr(nr,3*4);
	HaMat_double qp(nr,3*4);
	
	dif.newsize(nr);
	dif = 0.0;

	HaResidue* prev;
	HaResidue* pres;

	int k;
	for(k = 0; k < 2; k++)
	{
		if( k == 0) ch = ch1;
		if( k == 1) ch = ch2;

		int ir = 0;
		ResidueIteratorChain ritr_ch(ch);
		prev = ritr_ch.GetFirstRes();
		pres = ritr_ch.GetNextRes();
		for(; pres != NULL; prev = pres,pres = ritr_ch.GetNextRes())
		{
			ir++;

			Vec3D r1,r2;
			double d1,c1,d2,c2;
			for(j = 0; j < 3; j++)
			{
				r1 = prev->std_crd_sys[j];
				d1 = Vec3D::DotProduct(u[ir-1],r1);
				c1 = Vec3D::DotProduct(b[nr*k + ir-1],r1);

				r2 = pres->std_crd_sys[j];
				d2 = Vec3D::DotProduct(u[ir],r2);
				c2 = Vec3D::DotProduct(b[nr*k+ir],r2);
				
//				PrintLog(" k ir j c1 d1 c2 d2 = %3d %3d %3d %12.6f %12.6f %12.6f %12.6f \n",
//					       k,ir,j,c1,d1,c2,d2);

				qr(ir+1,k*3+j+1)= d2 - d1;
				qp(ir+1,k*3+j+1)= c2 - c1;

				double si1 = (d2 - d1)*(d2 - d1);
				dif[ir-1] += si1*10;
				scp[0] += si1;

				double si2 = (c2 - c1)*(c2 - c1);
				dif[ir-1] += si2;
				scp[1] += si2;
			}

			double qi1 = 0.0;
			for(j = 0; j < 3; j++)
			{
			   up[ir][j] = u[ir][j] - u[ir-1][j];
			   qi1 += up[ir][j]*up[ir][j];
			}
			dif[ir-1] += qi1*10;
			scp[2] += qi1;

			double dot = 0.0;
			for(j = 0; j < 3; j++)
			{
			   us[ir][j] = u[ir][j] + u[ir-1][j];
			   dot += us[ir][j]*us[ir][j];
			}
			for(j = 0; j < 3; j++)
			{
				um[ir][j] = us[ir][j]/dot;
			}
			dot = Vec3D::DotProduct(us[ir],s[ir]);
			double qi2 = 0.0;
			for(j = 0; j < 3; j++)
			{
				q[ir][j] = s[ir][j] - um[ir][j]*dot;
				qi2 +=  q[ir][j]*q[ir][j];
			}
			dif[ir-1] += qi2;
			scp[3] += qi2;
		}
	}
	sum = 10.0*(scp[0] + scp[2]) + scp[1] + scp[3];

// Checked until here:

	Vec3D zero_v;
	zero_v.SetX(0.0);
	zero_v.SetY(0.0);
	zero_v.SetZ(0.0);

// Compute derivative of u and p

	pres1 = ritr_ch_1.GetFirstRes();
    pres2 = ritr_ch_2.GetFirstRes();
	ir = -1;	
	for(; pres1; pres1 = ritr_ch_1.GetNextRes(),pres2 = ritr_ch_2.GetNextRes() )
	{
		ir++;
		Vec3D d,e,v;
			
		d = pres1->std_crd_sys[1];
		Vec3D::diff(e,pres1->std_crd_sys[3],h[ir]);
		Vec3D::VecProduct(v,d,u[ir]);
		
		ud[3*nr + ir] = -1.0*v;
		v.normalize();
		
		Vec3D::VecProduct(ud[2*nr + ir],u[ir],v);
		ud[1*nr + ir] = zero_v;
		ud[0*nr + ir] = zero_v;

		pd[0*nr + ir] = -1.0*v;
		pd[1*nr + ir] = -1.0*ud[2*nr + ir];
		Vec3D::VecProduct(pd[2*nr + ir],v,e);
		Vec3D::VecProduct(pd[3*nr + ir],d,e);

		int l;
		for(l = 0; l < 4; l++)
		{
			Vec3D vv = pd[l*nr + ir];
//      	PrintLog("PD1: (ir = %3d l=%3d) = %12.6f %12.6f %12.6f \n",ir,l,vv[0],vv[1],vv[2]);
		}
// Checked to here

		HaVec_double dot3(4);
		for(l= 0; l < 4; l++)
		{
			dot3[l] = Vec3D::DotProduct(u[ir],pd[l*nr + ir]);
		}	
		
		v = pres2->std_crd_sys[3] - h[ir];
		double dot1 = Vec3D::DotProduct(u[ir],v);

		for(l= 0; l < 4; l++)
		{
			double dot2 = Vec3D::DotProduct(ud[l*nr + ir],v);
//			PrintLog("PD1: ir,l, dot1,dot2,dot3 %3d %3d %12.6f %12.6f %12.6f \n",
//				      ir,l,dot1,dot2,dot3[l]);
			pd[l*nr + ir] = pd[l*nr + ir] + 0.5*dot1*ud[l*nr + ir] + 0.5*(dot2 - dot3[l])*u[ir];
		}
	}

	for(ir = 0; ir < nr; ir++)
	{
		for(int l = 0; l < 4; l++)
		{
//			Vec3D vv = pd[l*nr + ir];
//			PrintLog("PD (ir = %3d l=%3d) = %12.6f %12.6f %12.6f \n",ir,l,vv[0],vv[1],vv[2]);
//			Vec3D vv = ud[l*nr + ir];
//			PrintLog("UD (ir = %3d l=%3d) = %12.6f %12.6f %12.6f \n",ir,l,vv[0],vv[1],vv[2]);
		}
	}


	gra.newsize(4*nr);
	gra = 0.0;
	
	for(k = 0; k < 2; k++)
	{
	   if( k == 0) ch = ch1;
	   if( k == 1) ch = ch2;
	   int j;
	   for(j=0; j <3; j++)
	   {
		  HaVec_double du1(4,0.0);
		  HaVec_double du2(4,0.0);
		  HaVec_double dv1(4,0.0);
		  HaVec_double dv2(4,0.0);

		  ResidueIteratorChain ritr_ch(ch);

		  HaResidue* prev = ritr_ch.GetFirstRes();
	      HaResidue* pres = ritr_ch.GetNextRes();

	      int ir = 0;
		  int l,m;
		  double dot,dot1,dot2,dot3,dot4;
		  double ds1,ds2,dq1,dq2;
	      for(; pres; prev= pres,pres = ritr_ch.GetNextRes())
		  {
		     ir++;
		     Vec3D r1,r2;
		     r2 = pres->std_crd_sys[j];
		     r1 = prev->std_crd_sys[j];
		     
		     for(l = 0; l < 4; l++)
			 {
				m = (ir -1)*4 + l;	
				dot = Vec3D::DotProduct(ud[l*nr+ ir-1],r1);
				ds1 = -20.0* qr(ir+1,k*3+j+1)*dot;
				dot = Vec3D::DotProduct(pd[l*nr+ ir-1],r1);
				ds2 = -2.0* qp(ir+1,k*3+j+1)*dot;
				gra[m] += ds1 + ds2 + du1[l] + du2[l];
				dot = Vec3D::DotProduct(ud[l*nr+ ir],r2);
				du1[l] = 20.0*qr(ir+1,k*3+j+1)*dot;
				dot = Vec3D::DotProduct(pd[l*nr+ ir],r2);
				du2[l] = 2.0*qp(ir+1,k*3+j+1)*dot;
//				PrintLog("G1:k,j,ir,l,ds1,ds2,du1[l],du2[l]= %3d %3d %3d %3d %12.6f %12.6f %12.6f %12.6f \n",
//					      k,j,ir,l,ds1,ds2,du1[l],du2[l]);
			 }
		  }
          if(j == 0)
		  {
		     ResidueIteratorChain ritr_ch(ch);
			 HaResidue* prev = ritr_ch.GetFirstRes();
	         HaResidue* pres = ritr_ch.GetNextRes();
	         int ir = 0;
			 Vec3D dq;
	         for(; pres; prev= pres,pres = ritr_ch.GetNextRes())
			 {
		        ir++;
				for(l= 0 ; l < 4; l++)
				{
				  m = (ir -1)*4 + l;	
				  dq1 = -20.0*Vec3D::DotProduct(up[ir],ud[l*nr+ir-1]);
				  dot1 = Vec3D::DotProduct(ud[l*nr+ir-1],s[ir]);
				  dot2 = Vec3D::DotProduct(us[ir],pd[l*nr+ir-1]);
				  dot3 = Vec3D::DotProduct(um[ir],s[ir]);
				  dot4 = Vec3D::DotProduct(us[ir],ud[l*nr+ir-1]);
				  dq = -1.0*pd[l*nr+ir-1] - (dot1-dot2)*um[ir] 
				       -dot3* ud[l*nr+ir-1] + 2.0*dot3*dot4*um[ir];
				  dq2 = 2.0*Vec3D::DotProduct(q[ir],dq);
				  gra[m] += dq1 + dq2 + dv1[l] + dv2[l];
				  dv1[l]= 20.0*Vec3D::DotProduct(up[ir],ud[l*nr+ir]);
				  dot1 = Vec3D::DotProduct(ud[l*nr+ir],s[ir]);
				  dot2 = Vec3D::DotProduct(us[ir],pd[l*nr+ir]);
				  dot4 = Vec3D::DotProduct(us[ir],ud[l*nr+ir]);
				  dq = pd[l*nr+ir] - (dot1+dot2)*um[ir] 
				       -dot3* ud[l*nr+ir] + 2.0*dot3*dot4*um[ir];
				  dv2[l] = 2.0*Vec3D::DotProduct(q[ir],dq);
//				PrintLog("G2:k,j,ir,l,dq1,dq2,dv1[l],dv2[l]= %3d %3d %3d %3d %12.6f %12.6f %12.6f %12.6f \n",
//					      k,j,ir,l,dq1,dq2,dv1[l],dv2[l]);
				}
			 }
		  }
		  m = (nr - 1)*4;
		  for(l =0; l <4; l++)
		  {
			  gra[m+l] += du1[l] + du2[l] + dv1[l] + dv2[l];
		  }
	   }
	}

	int m= -1;
	for(ir = 0 ; ir < nr; ir++)
	{
		for(j=0; j < 4; j++)
		{
			m++;
//			PrintLog(" gra(%3d)= %12.6f \n",m,gra[m]);
		}
	}
	return TRUE;
}

int NuclAcidMod::CalcBend()
{
	if(p_dna_mol == NULL) 
	{
		HaMolecule* pmol = FindDNAMol();
		if(pmol == NULL) return FALSE;
	}

	HaChain* ch;
	ChainIteratorMolecule chitr(p_dna_mol);
	HaChain* ch1 = chitr.GetFirstChain();
	HaChain* ch2 = chitr.GetNextChain();

	int nr = ch1->GetNRes();

	if( ch1 == NULL || ch2 == NULL) return FALSE;

	HaResidue* pres1 = ch1->GetFirstRes();
	HaResidue* pres2 = ch2->GetFirstRes();

	hel.newsize(6,nr*4);
	hold.newsize(2,nr*4);
	vold.newsize(4,nr);
	bend.newsize(nr,4);
	vkin.newsize(7,nr*4);

	uho.resize(nr*4);
	hho.resize(nr*4);
	ul.resize(nr*4);

	Vec3DValArray p(nr*4),w(nr*4),v(nr*4);
	Vec3DValArray va(nr);

	int k,ir;
    int m = -1;
	for(k=0; k < 2; k++)
	{
		for(ir = 0; ir < nr; ir++)
		{		
			m++;
			uho[m] = u[ir];
			hho[m] = o[ir];
		}
	}

	for(k=0; k < 2; k++)
	{
		if(k==0) ch = ch1;
		if(k==1) ch = ch2;

// local axis system

		ResidueIteratorChain ritr_ch(ch);
		HaResidue* pres = ritr_ch.GetFirstRes();
		ir = -1;
		hel(3,nr*k+1) = 0.0;
		hel(6,nr*k+1) = 0.0;
		for( ; pres; pres = ritr_ch.GetNextRes())
		{
			ir++;
			if(k > 0) ul[k*nr + ir] = -1.0* ul[ir];
			else 
				ul[ir] = uho[ir];

			Vec3D vc,da,e,f,t,q;
			vc = pres->std_crd_sys[3] - hho[nr*k+ir];
			double dot= Vec3D::DotProduct(vc,ul[nr*k + ir]);
			p[nr*k+ir] = hho[nr*k+ir] + dot*ul[nr*k + ir];
			da = pres->std_crd_sys[1];
			dot = Vec3D::DotProduct(ul[nr*k + ir],da);
			vc = da - dot*ul[nr*k + ir];
			w[nr*k + ir] = vc;
			w[nr*k + ir].normalize();
		    if(k == 0 && ir > 0)
			{
				double dit = Vec3D::DotProduct(w[k*nr+ir-1],w[k*nr+ir]);
				if(dit < 0.0) w[k*nr+ir].Scale(-1.0);
			}
			Vec3D::VecProduct(v[k*nr+ir],w[k*nr+ir],ul[k*nr+ir]);

            e = pres->std_crd_sys[3] - p[nr*k+ir];
			f = pres->std_crd_sys[2];

			hel(1,nr*k+ir+1) = Vec3D::DotProduct(e,v[k*nr+ir]);
			hel(2,nr*k+ir+1) = Vec3D::DotProduct(e,w[k*nr+ir]);
			dot = Vec3D::DotProduct(w[k*nr+ir],da);
			if( fabs(dot) > 1.0) dot = (dot < 0.0)? -1.0 : 1.0;
			hel(4,nr*k+ir+1) = acos(dot)*RAD_TO_DEG;
			Vec3D::VecProduct(t,w[k*nr+ir],da);
			if(Vec3D::DotProduct(t,v[k*nr+ir]) < 0.0) hel(4,nr*k+ir+1) *= (-1.0);
			Vec3D::VecProduct(q,v[k*nr+ir],da);
			double rq = q.length();
			dot = Vec3D::DotProduct(q,f)/rq;
			if( fabs(dot) > 1.0) dot = (dot < 0.0)? -1.0 : 1.0;
			hel(5,nr*k+ir+1) = acos(dot)*RAD_TO_DEG;
			Vec3D::VecProduct(t,q,f);
			if(Vec3D::DotProduct(t,da) < 0.0) hel(5,nr*k+ir+1) *= (-1.0);
		}

// Overall bend

		if( k == 0 )
		{
			double dot;
			Vec3D d,t,r,c,v1,n1;
			int ic;
			for(ic = 1; ic < (nr-1); ic++)
			{
				d = o[nr-1] - o[0];
				t = o[ic] - o[0];
				d.normalize();
				r = ul[k*nr+ic];
				double drd = Vec3D::DotProduct(r,d);
				double drt = Vec3D::DotProduct(r,t);
				c = (drt/drd)*d - t;
				double rc = c.length();
				v1 = v[k*nr+ic];
				dot = Vec3D::DotProduct(c,v1)/rc;
				if(fabs(dot) > 1.0) dot = (dot < 0.0)? -1.0 : 1.0;
				bend(ic+1,k+1) = acos(dot) * RAD_TO_DEG;
				Vec3D::VecProduct(n1,v1,c);
				if(Vec3D::DotProduct(r,n1) < 0.0) bend(ic+1,k+1) *= (-1.0);
			}
			dot = Vec3D::DotProduct(ul[k*nr],ul[k*nr+nr-1]);
			if(fabs(dot) > 1.0) dot = (dot < 0.0)? -1.0 : 1.0;
			bend(1,k+1) = acos(dot) *RAD_TO_DEG;
			d = o[nr-1] - o[nr-2];
			t = o[1] - o[0];
			double rd = d.length();
			double rt = t.length();
			dot = Vec3D::DotProduct(d,t)/(rd*rt);
			if(fabs(dot) > 1.0) dot = (dot < 0.0)? -1.0 : 1.0;
			bend(nr,k+1) = acos(dot)*RAD_TO_DEG;
		}
	}
//	for(ir = 1; ir <= nr; ir++)
//	{
//		PrintLog("BEND1: %3d %12.6f\n",ir,bend(ir,1));
//	}

	for(ir = 0; ir < nr; ir++)
	{
		Vec3D vc;
		vc[0] = 0.0; vc[1] = 0.0; vc[2] = 0.0;

		double dot;
		for(k = 0; k < 2; k++)
		{
			if(k == 0) dot = -1.0;
			else
				dot = Vec3D::DotProduct(ul[ir],ul[nr*k + ir]);
			double sig = (dot > 0)? -1.0 : 1.0;
			vc[0] += v[k*nr+ir][0]*sig;
			vc[1] += v[k*nr+ir][1];
			vc[2] += v[k*nr+ir][2];
		}
		double r = vc.length();
		va[ir] = (1.0/r)*vc;
	}

//    for(k =0; k < 2; k++)
//	{
//		for(ir = 0; ir < nr; ir++)
//		{
//			PrintLog("UL: %3d %3d %12.6f %12.6f %12.6f \n",
//				k,ir,ul[nr*k + ir][0],ul[nr*k + ir][1],ul[nr*k + ir][2]);
//		}
//	}


	for(k=0; k < 2; k++)
	{
		Vec3D n,vc,q,d,f,pl,pu,t,r,fp;
		for(ir = 1; ir < nr; ir++)
		{
//          mean plane axis system
			
			double dot;

			n = ul[ir-1] + ul[ir];
			n.normalize();
			vc = o[ir-1] - o[ir];
			vkin(7,ir+1) = vc.length();
			q= 0.5*(o[ir-1] + o[ir]);
			vc = va[ir-1] + va[ir];
			dot = Vec3D::DotProduct(n,vc);
			d = vc - dot*n;
			d.normalize();
			Vec3D::VecProduct(f,n,d);
			int ka = 0;
			int km = k;
			int kp = k;
// u vector intersection
			vc = q - p[km*nr + ir-1];
			double dot1 = Vec3D::DotProduct(n,vc);
			double dot2 = Vec3D::DotProduct(n,ul[ka*nr+ir-1]);
			double dl = dot1/dot2;
			pl = p[km*nr + ir-1] + dl*ul[ka*nr+ir-1];
			vc = p[km*nr + ir] - q;
			dot1 = Vec3D::DotProduct(n,vc);
			dot2 = Vec3D::DotProduct(n,ul[ka*nr+ir]);
			double du = dot1/dot2;
			pu = p[km*nr + ir] - du*ul[ka*nr+ir];
			hel(3,nr*k+ir+1) = dl + du;
// shift parameter
			vc = pu - pl;
			vkin(1,k*nr+ir+1) = Vec3D::DotProduct(d,vc);
			int idir = 1;
			if(k==1) idir = -1;
			vkin(2,k*nr+ir+1) = idir*Vec3D::DotProduct(f,vc);
			vkin(5,k*nr+ir+1) = vc.length();
			dot = Vec3D::DotProduct(ul[ka*nr+ir-1],ul[ka*nr+ir]);
			if(fabs(dot) > 1.0) dot = (dot < 0.0)? -1.0 : 1.0;
			vkin(6,k*nr+ir+1) = acos(dot)*RAD_TO_DEG;
// kink pars
			Vec3D::VecProduct(t,ul[ka*nr+ir],d);
			double rt = t.length();
			dot = Vec3D::DotProduct(f,t)/rt;
			if(fabs(dot) > 1.0) dot = (dot < 0.0)? -1.0 : 1.0;
			double cln = acos(dot)*RAD_TO_DEG;
			Vec3D::VecProduct(vc,f,t);
			if(Vec3D::DotProduct(vc,d) < 0.0) cln = -cln;
			Vec3D::VecProduct(r,d,t);
			double rr = r.length();
			dot = Vec3D::DotProduct(ul[ka*nr+ir],r)/rr;
			if(fabs(dot) > 1.0) dot = (dot < 0.0)? -1.0 : 1.0;
			double tip = acos(dot)*RAD_TO_DEG;
			Vec3D::VecProduct(vc,r,ul[ka*nr+ir]);
			if(Vec3D::DotProduct(vc,t) < 0.0) tip = -tip;
			vkin(3,k*nr+ir+1) = 2.0* cln;
			vkin(4,k*nr+ir+1) = 2.0* tip * idir;
			hel(6,nr*k+ir+1)=0.0;
			int l;
			for(l = ir-1; l <= ir; l++)
			{
				double sa = sin(DEG_TO_RAD*(-cln));
				double ca = cos(DEG_TO_RAD*cln);
				if(l == ir) sa = sin(DEG_TO_RAD*cln);
				fp = f;
				fp.Rotate(d,ca,sa);
				dot = Vec3D::DotProduct(fp,w[k*nr+l]);
				if(fabs(dot) > 1.0) dot = (dot < 0.0)? -1.0 : 1.0;
				double wdg = acos(dot)*RAD_TO_DEG;
				Vec3D::VecProduct(vc,fp,w[k*nr+l]);
				dot = Vec3D::DotProduct(vc,ul[ka*nr+l]);
				if(l == (ir -1 ) && dot > 0) wdg = - wdg;
				if(l == ir       && dot < 0) wdg = - wdg;
				hel(6,nr*k+ir +1) += wdg;
			}
			hel(6,nr*k+ir+1) -= 360.0*floor( hel(6,nr*k+ir+1)/360.0);
			if(fabs(hel(6,nr*k+ir+1)) > 180.0) 
				hel(6,nr*k+ir+1) = (hel(6,nr*k+ir+1) > 0)? hel(6,nr*k+ir+1) - 360.0: hel(6,nr*k+ir+1) + 360.0;

// Jumna style paramters
			if(k == 0)
			{
				Vec3D ps;
				d = v[ir-1];
				f = w[ir-1];
				vc = p[ir] - p[ir-1];
				double dot1 = Vec3D::DotProduct(ul[ir-1],vc);
				double dot2 = Vec3D::DotProduct(ul[ir-1],ul[ir]);
				double dd = dot1/dot2;
				hold(1,ir+1) = dd;
				ps = p[ir] - dd*ul[ir];
				vc = ps - p[ir-1];
				vold(1,ir+1) = Vec3D::DotProduct(d,vc);
                vold(2,ir+1) = Vec3D::DotProduct(f,vc);
				Vec3D::VecProduct(t,ul[ir],d);
				double rt = t.length();
				double dot = Vec3D::DotProduct(f,t)/rt;
				if(fabs(dot) > 1.0) dot = (dot < 0.0)? -1.0 : 1.0;
			    double cln = acos(dot)*RAD_TO_DEG;
				Vec3D::VecProduct(vc,f,t);
				if(Vec3D::DotProduct(vc,d) < 0.0) cln = -cln;
				Vec3D::VecProduct(r,d,t);
				double rr = r.length();
				dot = Vec3D::DotProduct(ul[ir],r)/rr;
				if(fabs(dot) > 1.0) dot = (dot < 0.0)? -1.0 : 1.0;
			    double tip = acos(dot)*RAD_TO_DEG;
				Vec3D::VecProduct(vc,r,ul[ir]);
				if(Vec3D::DotProduct(vc,t) < 0.0) tip = -tip;
				vold(3,ir+1) = cln;
				vold(4,ir+1) = tip;
				double ca = cos(DEG_TO_RAD*cln);
				double sa = sin(DEG_TO_RAD*cln);
				fp = f;
				fp.Rotate(d,ca,sa);
				dot = Vec3D::DotProduct(fp,w[ir]);
				if(fabs(dot) > 1.0) dot = (dot < 0.0)? -1.0 : 1.0;
			    double wdg = acos(dot)*RAD_TO_DEG;
				Vec3D::VecProduct(vc,fp,w[ir]);
				if(Vec3D::DotProduct(vc,ul[ir]) < 0)wdg = -wdg;
				hold(2,ir+1) = wdg;
			}
		}
	}
	for(k=1; k < 2; k++)
	{
		double dzl = 0.0;
		double dwl = 0.0;
		for(ir = 0; ir < nr; ir++)
		{
			Vec3D vc;
			vc = p[k*nr+ir] - p[ir];
			double dzu = vc.length();
			if(Vec3D::DotProduct(vc,ul[ir]) < 0) dzu = -dzu;
			double dot = - Vec3D::DotProduct(w[ir],w[k*nr+ir]);
			if(fabs(dot) > 1.0) dot = (dot < 0.0)? -1.0 : 1.0;
			double dwu = acos(dot)*RAD_TO_DEG;
			Vec3D::VecProduct(vc,w[ir],w[k*nr+ir]);
			if(Vec3D::DotProduct(vc,ul[ir]) > 0) dwu = -dwu;
			hold(1,k*nr+ir+1) = hold(1,ir+1) + dzu - dzl;
			hold(2,k*nr+ir+1) = hold(2,ir+1) + dwu - dwl;
			if(ir == 0 && ir >= 0) 
			{
				hel(3,nr*k+ir+1) += dzu;
				hel(6,nr*k+ir+1) += dwu;
			}
			dzl = dzu;
			dwl = dwu;
		}
	}

//	for(k = 0; k < 2; k++)
//	{
//		for(ir = 0; ir < nr; ir++)
//		{
//			int m = nr*k+ir+1;
//			PrintLog(" HEL2: %3d %3d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f b\n",
//				      k,ir,hel(1,m),hel(2,m),hel(3,m),hel(4,m),hel(5,m),hel(6,m));
//		}
//	}

//	for(k = 0; k < 2; k++)
//	{
//		for(ir = 0; ir < nr; ir++)
//		{
//			int m = nr*k+ir+1;
//			PrintLog(" HOLD: %3d %3d %10.4f %10.4f \n",
//				      k,ir,hold(1,m),hold(2,m));
//		}
//	}

//	for(ir = 0; ir < nr; ir++)
//	{
//		int m = ir+1;
//		PrintLog(" VOLD: %3d %10.4f %10.4f %10.4f %10.4f \n",
//				  ir,vold(1,m),vold(2,m),vold(3,m),vold(4,m));
//	}

//	for(k = 0; k < 2; k++)
//	{
//		for(ir = 0; ir < nr; ir++)
//		{
//			int m = nr*k+ir+1;
//			PrintLog(" VKIN: %3d %3d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",
//				      k,ir,vkin(1,m),vkin(2,m),vkin(3,m),vkin(4,m),vkin(5,m),vkin(6,m),vkin(7,m));
//		}
//	}

	return TRUE;
}

static double av_ang(double a1,double a2)
{
	double ca = cos(DEG_TO_RAD*a1);
	double sa = sin(DEG_TO_RAD*a1);
	double cb = cos(DEG_TO_RAD*a2);
	double sb = sin(DEG_TO_RAD*a2);
	double diff = acos(ca*cb + sa*sb)*RAD_TO_DEG;
	if((cb*sa - sb*ca) < 0.0) diff = -diff;
	double ff = a2 + diff/2.0;
	return ff;
}

static double diff_ang(double a1,double a2)
{
	double ca = cos(DEG_TO_RAD*a1);
	double sa = sin(DEG_TO_RAD*a1);
	double cb = cos(DEG_TO_RAD*a2);
	double sb = sin(DEG_TO_RAD*a2);
	double diff = acos(ca*cb + sa*sb)*RAD_TO_DEG;
	if((cb*sa - sb*ca) < 0.0) diff = -diff;
	return diff;
}

int 
NuclAcidMod::CalcGlobHlxCrd()
{
	if(p_dna_mol == NULL)
	{
		HaMolecule* pmol = FindDNAMol();
		if(pmol == NULL) return FALSE;
	}

	if(hel.num_rows() == 0) 
	{
		CalcAxis();
	    CalcBend();
	}
	
	HaChain* ch;
	ChainIteratorMolecule chitr(p_dna_mol);
	HaChain* ch1 = chitr.GetFirstChain();
	HaChain* ch2 = chitr.GetNextChain();

	if( ch1 == NULL || ch2 == NULL) return FALSE;

	int nr = ch1->GetNRes();

	StrVec lbl(nr*4);

	HaResidue* pres;

	int m = -1;
	int k;
	for(k = 0; k < 2; k++)
	{
		if(k == 0) ch = ch1;
		if(k == 1) ch = ch2;
		ResidueIteratorChain ritr_ch(ch);
		pres = ritr_ch.GetFirstRes();
		for(; pres; pres = ritr_ch.GetNextRes())
		{
			m++;
			char buf[256];
			pres->FillRef(buf,1);
			lbl[m] = buf;
		}
	}

	PrintLog(" B: Global Base-Axis Parameters: \n");
	PrintLog("         Xdisp    Ydisp    Inclin      Tip  \n");
	PrintLog("          (dx)     (dy)     (eta)    (theta)  \n\n");

	int ir;
	for(k = 0; k < 2; k++)
	{
		for(ir = 0; ir < nr; ir++)
		{
			m = k*nr+ir+1;
			PrintLog("%s %9.3f %9.3f %9.3f %9.3f \n",lbl[m-1].c_str(),hel(1,m),hel(2,m),hel(4,m),hel(5,m));
		}
		PrintLog("\n");
	}
	PrintLog("\n\n");

	HaVec_double hav(6);

	PrintLog(" Global Base Pair - Axis Parameters: \n");
	PrintLog("               Xdisp     Ydisp   Inclin       Tip  \n");
	PrintLog("               (dx)       (dy)    (eta)     (theta)  \n\n");
	
	hav = 0.0;
	
	for(ir = 0; ir < nr; ir++)
	{
		double xdi = (hel(1,ir+1) + hel(1,nr+ir+1))/2.0;
		double ydi = (hel(2,ir+1) - hel(2,nr+ir+1))/2.0;
		double cln = av_ang(hel(4,ir+1),hel(4,nr+ir+1));
		double tip = av_ang(hel(5,ir+1),-hel(5,nr+ir+1));
		if(fabs(cln) > 180.0) cln = (cln > 0)? cln - 360.0: cln + 360.0;
		if(fabs(tip) > 180.0) tip = (tip > 0)? tip - 360.0: tip + 360.0;
		
		hav(1) += xdi;
		hav(2) += ydi;
		hav(4) += cln;
		hav(5) += tip;

		PrintLog("%s-%s %9.3f %9.3f %9.3f %9.3f \n",
			      lbl[ir].c_str(),lbl[nr+ir].c_str(),xdi,ydi,cln,tip);
	}
	PrintLog("\n");
    PrintLog("Averages:   %9.3f %9.3f %9.3f %9.3f \n",
			             hav(1)/nr,hav(2)/nr,hav(4)/nr,hav(5)/nr);

	
	PrintLog("\n\n");

	bs_bs_pars.newsize(6,nr);

	PrintLog(" D: Global Base - Base Parameters: \n");
	PrintLog("              Shear     Stretch   Stagger   Buckle    Propel    Opening '\n");
	PrintLog("               (Sx)       (Sy)     (Sz)     (kappa)   (omega)   (sigma)  \n\n");
	
	hav = 0.0;
	double stg = 0.0;
	double opn = 0.0;

	for(ir = 1; ir <= nr; ir++)
	{
		double str = hel(2,ir) + hel(2,nr+ir);
		double pro = diff_ang(hel(5,ir),-hel(5,nr+ir));
		if(fabs(pro) > 180.0) pro = (pro > 0)? pro - 360.0: pro + 360.0;
		stg += (hel(3,ir) - hel(3,nr+ir));
		opn += (hel(6,ir) - hel(6,nr+ir));

		double shr = hel(1,ir) - hel(1,nr+ir);
		double buc = diff_ang(hel(4,ir),hel(4,nr+ir));
		
		bs_bs_pars(1,ir) = shr;
		bs_bs_pars(2,ir) = str;
		bs_bs_pars(3,ir) = stg;
		bs_bs_pars(4,ir) = buc;
		bs_bs_pars(5,ir) = pro;
		bs_bs_pars(6,ir) = opn;

		hav(1) += shr;
		hav(2) += str;
		hav(3) += stg;
		hav(4) += buc;
		hav(5) += pro;
		hav(6) += opn;

		PrintLog("%s-%s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
			      lbl[ir-1].c_str(),lbl[nr+ir-1].c_str(),shr,str,stg,buc,pro,opn);
	}
	PrintLog("\n");
    PrintLog("Averages:   %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f \n",
			             hav(1)/nr,hav(2)/nr,hav(3)/nr,hav(4)/nr,hav(5)/nr,hav(6)/nr);

	PrintLog("\n\n");

	PrintLog(" E: Global Inter-Base Parameters: \n");
	PrintLog("              Shift    Slide       Rise      Tilt     Roll     Twist   \n");
	PrintLog("               (Dx)     (Dy)       (Dz)     (tau)     (rho)   (Omega) \n\n");

	for(k = 0; k < 2; k++)
	{
		int idir = 1;
		if(k==1) idir = -1;
		for(ir = 2; ir <= nr; ir++)
		{
			double shift = hel(1,nr*k+ir) + vkin(1,ir) - hel(1,nr*k+ir-1);
			double slide = hel(2,nr*k+ir) + vkin(2,ir)*idir - hel(2,nr*k+ir-1);
			double rise  = hel(3,nr*k+ir);
			double tilt  = hel(4,nr*k+ir) + vkin(3,ir) - hel(4,nr*k+ir-1);
			double roll  = hel(5,nr*k+ir) + vkin(4,ir)*idir - hel(5,nr*k+ir-1);
			double twist = hel(6,nr*k+ir);
			if(fabs(roll) > 180.0) roll = (roll > 0)? roll - 360.0: roll + 360.0;

			m = nr*k + ir - 1;
			PrintLog("%s-%s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f \n",
				lbl[m-1].c_str(),lbl[m].c_str(),shift,slide,rise,tilt,roll,twist);
		}
		PrintLog("\n");
	}
	PrintLog("\n\n");

	PrintLog(" F: Global Inter-Base pair Parameters: \n");
	PrintLog("                        Shift    Slide    Rise     Tilt    Roll    Twist   \n");
	PrintLog("                         (Dx)     (Dy)    (Dz)     (tau)   (rho)   (Omega) \n\n");

	hav = 0.0;
	int nav = 0;
	k = 1;
	for(ir = 2; ir <= nr; ir++)
	{
		nav++;
		double xs = (hel(1,ir) + hel(1,nr*k+ir))/2.0;
		double xm = (hel(1,ir-1) + hel(1,nr*k+ir-1))/2.0;
		double ys = (hel(2,ir) - hel(2,nr*k+ir))/2.0;
		double ym = (hel(2,ir-1) - hel(2,nr*k+ir-1))/2.0;
		double ts = av_ang(hel(4,ir),hel(4,nr*k+ir));
		double tm = av_ang(hel(4,ir-1),hel(4,nr*k+ir-1));
		double ps = av_ang(hel(5,ir),-hel(5,nr*k+ir));
		double pm = av_ang(hel(5,ir-1),-hel(5,nr*k+ir-1));


		double shift = xs + vkin(1,ir) - xm;
		double slide = ys + vkin(2,ir) - ym;
		double rise  = (hel(3,ir)+hel(3,nr*k+ir))/2.0;
		double tilt  = ts + vkin(3,ir) - tm;
		double roll  = ps + vkin(4,ir) - pm;
		double twist = (hel(6,ir) + hel(6,nr*k+ir))/2.0;
		if(fabs(roll) > 180.0) roll = (roll > 0)? roll - 360.0: roll + 360.0;

		hav(1) += shift;
		hav(2) += slide;
		hav(3) += rise;
		hav(4) += tilt;
		hav(5) += roll;
		hav(6) += twist;

		int m1 = ir - 1;
		m = nr*k + ir - 1;
		PrintLog("%s-%s/%s-%s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n",
				lbl[m1-1].c_str(),lbl[m-1].c_str(),lbl[m1].c_str(),lbl[m].c_str(),
				shift,slide,rise,tilt,roll,twist);
	}
	PrintLog("\n");
    PrintLog("Averages:               %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n",
			             hav(1)/(nr-1),hav(2)/(nr-1),hav(3)/(nr-1),
						 hav(4)/(nr-1),hav(5)/(nr-1),hav(6)/(nr-1));
	
	PrintLog("\n\n");


	PrintLog(" I: Global Axis Curvature \n");
	PrintLog("              Ax      Ay     Ainc     Atip      Adis     Angle     Path  \n");
	PrintLog(" \n");

	for(ir = 2; ir <= nr; ir++)
	{
		nav++;
		double dsx = vkin(1,ir);
		double dsy = vkin(2,ir);
		double ktl = vkin(3,ir);
		double kpr = vkin(4,ir);
		double dis = vkin(5,ir);
		double the = vkin(6,ir);
		double pat = vkin(7,ir);
		
		PrintLog("%s-%s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n",
				lbl[ir-2].c_str(),lbl[ir-1].c_str(),
				dsx,dsy,ktl,kpr,dis,the,pat);
	}
	
	PrintLog("\n\n");
	PrintLog(" Overall Axis bend ... UU= %8.3f  PP= %8.3f  \n",bend(1,1),bend(nr,1));   
	PrintLog("\n\n");

	PrintLog(" Offset  L.Dir   wrt end-to-end vector \n");
	PrintLog(" \n");

	Vec3D v0 = hho[0];
	Vec3D n = hho[nr-1] - v0;
	double rn = n.length();
	n.normalize();
	double pl = 0.0;
	for(ir = 0; ir < nr; ir++)
	{
		if(ir > 0) pl += vkin(7,ir+1);
		Vec3D v = hho[ir] - v0;
		double dot = Vec3D::DotProduct(v,n);
		Vec3D d;
		d = v - dot*n;
		double rd = d.length();
		PrintLog("%s  %9.3f  %9.3f \n",lbl[ir].c_str(),rd,bend(ir+1,1));
	}
	PrintLog("\n\n");
	PrintLog(" Path length = %9.3f , End-to-End= %9.3f  Shortening = %9.3f percent \n",
               pl,rn,100.0*(1.0 - rn/pl));

	PrintLog("\n\n");

	return TRUE;
}

int NuclAcidMod::CalcBBCrd()
{
	if(p_dna_mol == NULL)
	{
		HaMolecule* pmol = FindDNAMol();
		if(pmol == NULL) return FALSE;
	}

	HaChain* ch;
	ChainIteratorMolecule chitr(p_dna_mol);
	HaChain* ch1= chitr.GetFirstChain();
	HaChain* ch2= chitr.GetNextChain();

	int nr = ch1->GetNRes();

	StrVec lbl(nr*4);

	HaResidue* pres;

	int m = -1;
	int k;
	for(k = 0; k < 2; k++)
	{
		if(k == 0) ch = ch1;
		if(k == 1) ch = ch2;
		ResidueIteratorChain ritr_ch(ch);
		pres = ritr_ch.GetFirstRes();
		for(; pres; pres = ritr_ch.GetNextRes())
		{
			m++;
			char buf[256];
			pres->FillRef(buf,1);
			lbl[m] = buf;
		}
	}

	int ir;

	tor.newsize(13,4*nr);
	sug.newsize(2,4*nr);

	for(k=0; k <2; k++)
	{
		if(k==0) ch = ch1;
		if(k==1) ch = ch2;
		
		HaResidue* pres;
		HaResidue* prev;
		HaResidue* pnext;

		HaVec_double store(16,0.0);

		ResidueIteratorChain ritr_ch(ch);
		pres = ritr_ch.GetFirstRes();
		ir = 0;
		for( ;pres; pres = ritr_ch.GetNextRes())
		{
			ir++;
			HaAtom* at_c8  =NULL; // 1
			HaAtom* at_n9  =NULL; // 2
			HaAtom* at_c1p =NULL; // 3
			HaAtom* at_c2p =NULL; // 4
			HaAtom* at_c3p =NULL; // 5
			HaAtom* at_c4p =NULL; // 6
			HaAtom* at_o1p =NULL; // 7
			HaAtom* at_o3p =NULL; // 8
			HaAtom* at_p   =NULL; // 9
			HaAtom* at_o5p =NULL; // 10
			HaAtom* at_c5p_nxt =NULL; // 11
			HaAtom* at_c4p_nxt =NULL; // 12
			HaAtom* at_c5p     =NULL; // 13
			HaAtom* at_o5p_prv =NULL; // 14
			HaAtom* at_c4      =NULL; // 15

			prev = pres->GetPrevResInChain();
			pnext = pres->GetNextResInChain();
 
			if(k == 1)
			{
				HaResidue* paxx = pnext;
				pnext = prev;
				prev  = paxx;
			}

			if(pres->IsPurine())
				at_c8 = pres->GetAtomByName("C8");
			else if(pres->IsPyrimidine())
				at_c8 = pres->GetAtomByName("C6");
			
			if(pres->IsPurine())
				at_n9  = pres->GetAtomByName("N9");
			else if(pres->IsPyrimidine())
				at_n9  = pres->GetAtomByName("N1");

			at_c1p = pres->GetAtomByName("C1X");
			if(at_c1p == NULL)at_c1p = pres->GetAtomByName("C1'");
			if(at_c1p == NULL)at_c1p = pres->GetAtomByName("C1*");

			at_c2p = pres->GetAtomByName("C2X");
			if(at_c2p == NULL)at_c2p = pres->GetAtomByName("C2'");
			if(at_c2p == NULL)at_c2p = pres->GetAtomByName("C2*");

			at_c3p = pres->GetAtomByName("C3X");
			if(at_c3p == NULL)at_c3p = pres->GetAtomByName("C3'");
			if(at_c3p == NULL)at_c3p = pres->GetAtomByName("C3*");

			at_c4p = pres->GetAtomByName("C4X");
			if(at_c4p == NULL)at_c4p = pres->GetAtomByName("C4'");
			if(at_c4p == NULL)at_c4p = pres->GetAtomByName("C4*");

			at_o1p = pres->GetAtomByName("O1X");
			if(at_o1p == NULL)at_o1p = pres->GetAtomByName("O4X");
			if(at_o1p == NULL)at_o1p = pres->GetAtomByName("O1'");
			if(at_o1p == NULL)at_o1p = pres->GetAtomByName("O1*");
			if(at_o1p == NULL)at_o1p = pres->GetAtomByName("O4'");
			if(at_o1p == NULL)at_o1p = pres->GetAtomByName("O4*");

			at_o3p = pres->GetAtomByName("O3X");
			if(at_o3p == NULL)at_o3p = pres->GetAtomByName("O3'");
			if(at_o3p == NULL)at_o3p = pres->GetAtomByName("O3*");

			if(pnext != NULL)
			{
				at_p = pres->GetAtomByName("P");
				if(at_p == NULL)at_p = pres->GetAtomByName("P");
				if(at_p == NULL)at_p = pres->GetAtomByName("P");

				at_o5p = pres->GetAtomByName("O5X");
				if(at_o5p == NULL)at_o5p = pres->GetAtomByName("O5'");
				if(at_o5p == NULL)at_o5p = pres->GetAtomByName("O5*");

				at_c5p_nxt = pnext->GetAtomByName("C5X");
				if(at_c5p_nxt == NULL)at_c5p_nxt = pnext->GetAtomByName("C5'");
				if(at_c5p_nxt == NULL)at_c5p_nxt = pnext->GetAtomByName("C5*");

				at_c4p_nxt = pnext->GetAtomByName("C4X");
				if(at_c4p_nxt == NULL)at_c4p_nxt = pnext->GetAtomByName("C4'");
				if(at_c4p_nxt == NULL)at_c4p_nxt = pnext->GetAtomByName("C4*");
			}

			at_c5p = pres->GetAtomByName("C5X");
			if(at_c5p == NULL)at_c5p = pres->GetAtomByName("C5'");
			if(at_c5p == NULL)at_c5p = pres->GetAtomByName("C5*");

			if(prev != NULL)
			{
				at_o5p_prv = prev->GetAtomByName("O5X");
			    if(at_o5p_prv == NULL)at_o5p_prv = prev->GetAtomByName("O5'");
			    if(at_o5p_prv == NULL)at_o5p_prv = prev->GetAtomByName("O5*");
			}

			if(pres->IsPurine())
				at_c4 = pres->GetAtomByName("C4");
			else if(pres->IsPyrimidine())
				at_c4 = pres->GetAtomByName("C2");

			store = 999.0;

			store(1) = Vec3D::CalcAngle(at_o1p,at_c1p,at_c2p)*RAD_TO_DEG;  
			store(2) = Vec3D::CalcAngle(at_c1p,at_c2p,at_c3p)*RAD_TO_DEG; 
	        store(3) = Vec3D::CalcAngle(at_c2p,at_c3p,at_c4p)*RAD_TO_DEG;  
	        store(4) = Vec3D::CalcTorsion(at_c8,at_n9,at_c1p,at_o1p)*RAD_TO_DEG;  
	        store(5) = Vec3D::CalcTorsion(at_o1p,at_c1p,at_c2p,at_c3p)*RAD_TO_DEG;  
	        store(6) = Vec3D::CalcTorsion(at_c1p,at_c2p,at_c3p,at_c4p)*RAD_TO_DEG;
            store(7) = Vec3D::CalcTorsion(at_c2p,at_c3p,at_c4p,at_o1p)*RAD_TO_DEG;
	        store(8) = Vec3D::CalcTorsion(at_c3p,at_c4p,at_o1p,at_c1p)*RAD_TO_DEG;
			store(9)  = Vec3D::CalcTorsion(at_c4p,at_o1p,at_c1p,at_c2p)*RAD_TO_DEG;

			if(at_c5p_nxt != NULL)
			{
				store(10) = Vec3D::CalcTorsion(at_c4p,at_c3p,at_o3p,at_p)*RAD_TO_DEG;
			    store(11) = Vec3D::CalcTorsion(at_c3p,at_o3p,at_p,at_o5p)*RAD_TO_DEG;
			}
			if(at_o5p_prv != NULL )
			{
				store(12) = Vec3D::CalcTorsion(at_o5p_prv,at_c5p,at_c4p,at_c3p)*RAD_TO_DEG;
			}
	        store(13) = Vec3D::CalcTorsion(at_c5p,at_c4p,at_c3p,at_o3p)*RAD_TO_DEG;
			if(at_c5p_nxt != NULL)
			{
				store(14) = Vec3D::CalcTorsion(at_o3p,at_p,at_o5p,at_c5p_nxt)*RAD_TO_DEG;
				store(15) = Vec3D::CalcTorsion(at_p,at_o5p,at_c5p_nxt,at_c4p_nxt)*RAD_TO_DEG;
			}
	        store(16) = Vec3D::CalcTorsion(at_c4,at_n9,at_c1p,at_o1p)*RAD_TO_DEG;

			int l;
//			PrintLog("AN: %3d%3d",k,ir);
//			for(l=1; l <= 16; l++)
//			{
//				PrintLog("%8.3f",store(l));
//			}
//			PrintLog("\n");
			
			for(l=1; l <=6; l++)
			{
				tor(l,nr*k+ir+1) = store(l);
			}
			for(l=10; l <=16; l++)
			{
				tor(l-3,nr*k+ir+1) = store(l);
			}

			double a = 0.0;
			double b = 0.0;
			
			for(l = 1; l <=5; l++)
			{
				int j = l+5;
				if(l == 5) j = 5;
				a += store(j)*cos(DEG_TO_RAD*(144.0*(l-1)));
				b += store(j)*sin(DEG_TO_RAD*(144.0*(l-1)));
			}
			a *= 0.4;
			b *= (-0.4);
			double amp = sqrt(a*a + b*b);
			double phase = 0.0;
			if(amp > 0.0)
			{
				double cp = a/amp;
				double sp = b/amp;
				if( fabs(cp) > 1.0) cp = (cp < 0.0)? -1.0 : 1.0;
				phase = acos(cp)*RAD_TO_DEG;
				if( sp < 0.0) phase = 360.0 - phase;
			}
			sug(1,nr*k+ir+1) = amp;
			sug(2,nr*k+ir+1) = phase;
			PrintLog("SUG: %s %8.3f %8.3f \n",lbl[nr*k+ir].c_str(),amp,phase);
		}
	}

	return TRUE;
}


class MinAx1 : public HaMinimizer
{
public:
	MinAx1(NuclAcidMod* new_pna_mod) { pna_mod = new_pna_mod;}
	virtual ~MinAx1() {}

	NuclAcidMod* pna_mod;

	virtual int CalcValGrad(HaVec_double& x, double& val, HaVec_double& grad)
	{
		HaVec_double scp(4),gra1;

		int ns = x.size();
		int nr = x.size()/4;
		
		HaMat_double hel(4,nr,x.begin());

		pna_mod->CalcAxisPar1(hel, val, gra1, scp);

		grad.newsize(ns);

		int ir,j,m;
		m = -1;
		for(ir = 0; ir < nr; ir++)
		{
			for(j = 0; j < 4; j++)
			{
				m++;
				if(j < 2)
					grad[m] = gra1[m];
				else
					grad[m] = gra1[m]*DEG_TO_RAD;

			}
		}

		nitr++;

		if(nitr == 1)
		{
			glast.newsize(ns);
		}
	
		double df = 0.0;
		if( nitr > 1)
		{
			df = fabs( val - flast);
		}

		int i;
		double gmax = 0.0;
		for(i=0; i < ns; i++)
		{
			gmax = MaxFun(fabs(grad[i]),gmax);
		}
	
		PrintLog(" iter= %5d  f= %12.6e  df = %12.6e  gmax = %12.6e \n",
			       nitr, val,df,gmax);
		
		flast = val;
		for(i=0; i < ns; i++)
			glast[i] = grad[i];

		return TRUE;

	}


};

int NuclAcidMod::CalcAxis()
{
	int lin_ax     = 0;
	int ibreak_ax  = 0;

	if(p_dna_mol == NULL)
	{
		HaMolecule* pmol = FindDNAMol();
		if(pmol == NULL) return FALSE;
	}
	
	ChainIteratorMolecule chitr(p_dna_mol);
	HaChain* ch1 = chitr.GetFirstChain();
	int nr = ch1->GetNRes();

	HaMat_double ax_hlxc(4,nr,0.0);
	double sum;
    HaVec_double gra;
	HaVec_double scp(4,0.0);

    CalcAxisPar1(ax_hlxc, sum, gra, scp);
	PrintLog("sum = %12.6f scp = %12.6f %12.6f %12.6f %12.6f \n",
	          sum,scp[0],scp[1],scp[2],scp[3]);


    MinAx1 min_ax(this);

	min_ax.SetNVar(4*nr);
	min_ax.scale.newsize(4*nr);

	int ir,j;
	int m = -1;
	for(ir = 0; ir < nr; ir++)
	{
		for(j=0; j < 4; j++)
		{
			m++;
			if(j < 2)
				min_ax.scale[m] = 0.5;
			else
				min_ax.scale[m] = 1.5;
		}	
	}


	HaVec_double v_init(4*nr,ax_hlxc.begin());
	min_ax.SetInitPoint(v_init);

	min_ax.Minimize();
	for(ir = 0; ir < nr; ir++)
	{
		PrintLog("%3d) U: %9.3f %9.3f %9.3f P: %9.3f %9.3f %9.3f",
			         (ir+1),u[ir][0],u[ir][1],u[ir][2],o[ir][0],o[ir][1],o[ir][2]);
		if( ir < (nr -1)) PrintLog(" D: %9.3f ",dif[ir]);
		PrintLog("\n");
	}

	return TRUE;
}

HaMolecule*
NuclAcidMod::FindDNAMol()
{
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	if(nmol == 0) return NULL;
	
	int i;
	for(i = 0; i < nmol; i++)
	{
		HaMolecule* pmol = pmset->GetMolByIdx(i);
		HaMolecule::ResidueIterator ritr(pmol);
		HaResidue* pres = ritr.GetFirstRes();
		if(pres == NULL) continue;
		if(pres->IsNucleo()) 
		{
			p_dna_mol = pmol;
			return pmol;
		}
	}

	return NULL;
}


int NuclAcidMod::CalcLocHlxCrd(int for_bp)
{
	int ires = CalcLocCrdSys();
	if(!ires) return FALSE;
	
	ChainIteratorMolecule chitr(p_dna_mol);
	HaChain* ch1 = chitr.GetFirstChain();
	HaChain* ch2 = chitr.GetNextChain();
	HaChain* ch;
	ResidueIteratorChain ritr_ch_1(ch1);
	ResidueIteratorChain ritr_ch_2(ch2);


	if( ch1 == NULL || ch2 == NULL) return FALSE;

	int nr = ch1->GetNRes();

	StrVec lbl(nr*4);

	HaResidue* pres;

	int m = -1;
	int k;
	for(k = 0; k < 2; k++)
	{
		if(k == 0) ch = ch1;
		if(k == 1) ch = ch2;
		ResidueIteratorChain ritr_ch(ch);
		pres = ritr_ch.GetFirstRes();
		for(; pres; pres = ritr_ch.GetNextRes())
		{
			m++;
			char buf[256];
			pres->FillRef(buf,1);
			lbl[m] = buf;
		}
	}


	int i,ir;
	char buf[256],buf2[256],buf3[256],buf4[256];

	int idir;

	if(!for_bp)
	{
		for(i=0; i < 2; i++)
		{
			if(i == 0) ch  =  ch1;
			if(i == 1) ch  =  ch2;
            if(i == 0) idir =  1;
			if(i == 1) idir = -1;

			ResidueIteratorChain ritr_ch(ch);
			HaResidue* prev = ritr_ch.GetFirstRes();
			HaResidue* pres = ritr_ch.GetNextRes();

			PrintLog(" Local Base helical parmeters: \n\n");
	
            ir = 0;
			for(; pres; prev= pres,pres = ritr_ch.GetNextRes())
			{
                Vec3D c0_1,v1_1,v2_1,v3_1;
				Vec3D c0_2,v1_2,v2_2,v3_2;
				ir++;

				if(pres->std_crd_sys.empty() || prev->std_crd_sys.empty())
				{
					PrintLog(" Std coord system are not set for some residues \n");
					continue;
				}

				c0_1 = prev->std_crd_sys[3];
				c0_2 = pres->std_crd_sys[3];
				v1_1 = prev->std_crd_sys[0];
				v1_2 = pres->std_crd_sys[0];
				v2_1 = prev->std_crd_sys[1];
				v2_2 = pres->std_crd_sys[1];
				v3_1 = prev->std_crd_sys[2];
				v3_2 = pres->std_crd_sys[2];
				
				double shift,slide,rise,tilt,roll,twist;

				Vec3D::CalcHlxParams(c0_1,v1_1,v2_1,v3_1,c0_2,v1_2,v2_2,v3_2,
					                 shift,slide,rise,tilt,roll,twist,idir);

				pres->FillRef(buf);
				prev->FillRef(buf2);

				PrintLog(" %s-%s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n", 
					      lbl[nr*i+ir-1].c_str(),lbl[nr*i+ir].c_str(),shift,slide,rise,tilt,roll,twist); 

			}
			PrintLog("\n\n");
		}
	}
	else
	{
		HaResidue* prev1 = ritr_ch_1.GetFirstRes();
		HaResidue* pres1 = ritr_ch_1.GetNextRes();
		HaResidue* prev2 = ritr_ch_2.GetFirstRes();
		HaResidue* pres2 = ritr_ch_2.GetNextRes();
		
		PrintLog(" Local Base-Pair helical parmeters: \n\n");

		ir = 0;
		for(; pres1; prev1= pres1,pres1 = ritr_ch_1.GetNextRes(),prev2= pres2,pres2 = ritr_ch_2.GetNextRes())
		{
			ir++;
			Vec3D c0_1,v1_1,v2_1,v3_1;
			Vec3D c0_2,v1_2,v2_2,v3_2;
			
			if(pres1->std_crd_sys.empty() || prev1->std_crd_sys.empty() || 
				pres2->std_crd_sys.empty() || prev2->std_crd_sys.empty() )
			{
				PrintLog(" Std coord system are not set for some residues \n");
				continue;
			}
			
			Vec3D::sum(c0_1, prev1->std_crd_sys[3], prev2->std_crd_sys[3]);
			c0_1.Scale(0.5);
			Vec3D::sum(c0_2, pres1->std_crd_sys[3], pres2->std_crd_sys[3]);
			c0_2.Scale(0.5);
			Vec3D::diff(v3_1,prev1->std_crd_sys[2],prev2->std_crd_sys[2]);
			Vec3D::diff(v3_2,pres1->std_crd_sys[2],pres2->std_crd_sys[2]);
	
			v3_1.normalize();
			v3_2.normalize();
			
			Vec3D::sum(v1_1, prev1->std_crd_sys[0],prev2->std_crd_sys[0]);
			Vec3D::sum(v1_2, pres1->std_crd_sys[0],pres2->std_crd_sys[0]);
			
			v1_1.normalize();
			v1_2.normalize();
			
			Vec3D::VecProduct(v2_1,v3_1,v1_1);
			Vec3D::VecProduct(v2_2,v3_2,v1_2);
			
			double shift,slide,rise,tilt,roll,twist;
			
			idir = 1;

			Vec3D::CalcHlxParams(c0_1,v1_1,v2_1,v3_1,c0_2,v1_2,v2_2,v3_2,
				shift,slide,rise,tilt,roll,twist,idir);
			
			prev1->FillRef(buf);
			pres1->FillRef(buf2);
			prev2->FillRef(buf3);
			pres2->FillRef(buf4);
					
			PrintLog(" %s-%s/%s-%s/:%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n", 
				lbl[ir-1].c_str(),lbl[nr+ir-1].c_str(),lbl[ir].c_str(),lbl[nr+ir].c_str(),
				shift,slide,rise,tilt,roll,twist); 
		}
	}
	return TRUE;
}

int NuclAcidMod::SetIntCoordsToJumna()
{
#if defined(INT_JUMNA)
    int na = mrc_.kam;
	
	HaMolecule* pMol= p_dna_mol;
	if(pMol == NULL)
	{
		ErrorInMod("NuclAcidMod::SetIntCoordToJumna()",
			       " DNA Molecule is not set");
		return FALSE;
	}
	if(pMol->GetNAtoms() != na)
	{
		PrintLog("Error in NuclAcidMod::SetIntCoordToJumna()\n");
		PrintLog("The number of atoms in the DNA molecule %d \n do not equal atoms in JUMNA %d \n",
			      pMol->GetNAtoms(),na);
		return FALSE;
	}

    int nr = pMol->GetNRes();
	int nrj = strjm_.nto;

	if( nr != nrj)
	{
		PrintLog("Error in NuclAcidMod::SetCoordToJumna()\n");
		PrintLog("The number of residues in the DNA molecule %d \n do not equal to the number residues in JUMNA %d \n",
			      nr,nrj);
		return FALSE;
	}

	if(hel_crd.size() != nr * 6)
	{
		PrintLog("Error in NuclAcidMod::SetCoordToJumna()\n");
		PrintLog(" Internal Coordinates Array does not have proper dimensions");
		return FALSE;
	}

	HaVec_double hel_axx(nr*6);

	int ir,i;
	for(ir = 0; ir < nr; ir++)
	{
		for(i = 0; i < 6; i++)
		{
			double sign = 1.0;
			if(i == 0 || i == 4) sign = -1.0; 
            hel_axx[6*ir + i] = sign * hel_crd[6*ir + i];
		}
		for(i = 0; i < 11; i++)
		{
//			if( i == 11 && Is3endRes(ir)) // the nucleotide is on 3' end  
//			{
//				strjm_.set[7*N2_J + ir] = bb_crd[12*ir + i];
//			}
//			else
//			{
				strjm_.set[i*N2_J + ir] = bb_crd[12*ir + i];
//			}
		}
	}

	HaChain* chain;
	HaChain* chain2;
	ChainIteratorMolecule chitr(p_dna_mol);
	if(IsSupHlxConstr() ) // Transform to global twist variables
	{
		ir = -1;
		for( chain = chitr.GetFirstChain(); chain; chain = chitr.GetNextChain() )
        {
			double tw = 0.0;
			HaResidue* pres;
			ResidueIteratorChain ritr_ch(chain);
			for(pres = ritr_ch.GetFirstRes(); pres; pres = ritr_ch.GetNextRes())
			{
				ir++;
				tw += hel_axx[6*ir + 5];
				if(tw >= 360.0) tw -= 360.0;
				hel_axx[6*ir + 5] = tw;
			}
		}
	}

    chain = chitr.GetFirstChain();
	nr = chain->GetNRes();
	ir = nr - 1;
	for( chain2 = chitr.GetNextChain(); chain; chain = chitr.GetNextChain())
	{
		double zloc = 0.0;
		double wloc = 0.0;
		for(int ir1 = 0; ir1 < nr; ir1++)
		{
			ir++;
			zloc = zloc + hel_axx[6*ir + 2] - hel_axx[6*ir1 + 2];
			wloc = wloc + hel_axx[6*ir + 5] - hel_axx[6*ir1 + 5];

			if(symjm_.ihl[ir])
			{
				hel_axx[6*ir + 2] = zloc;
				if(!IsSupHlxConstr()) 
				{
					hel_axx[6*ir + 5] = wloc;
				}
			}
		}
	}
	int ki = -1;
	int iv;
	nr = p_dna_mol->GetNRes();
	for(ir = 0; ir < nr; ir++)  // Set Backbone variables
	{
		int niv = NumIndBBVarRes(ir);
		for( iv = 0; iv < niv; iv++)
		{
			if( lock_bb[12*ir + iv] == 0 )
			{
				ki++;
				mnn_.var[ki] = bb_crd[12*ir + iv];
			}
		}
	}

	for(ir = 0; ir < nr; ir++)  // Set Helical variables
	{
		for( iv = 0; iv < 6; iv++)
		{
			if(lock_hel[6*ir + iv] == 0)
			{
				ki++;
				mnn_.var[ki] = hel_axx[6*ir + iv];
			}
//			else
//			{
//				strjm_.hel[iv*N2_J + ir] = hel_axx[6*ir + iv];
//			}
		}
	}

	for(ir = 0; ir < nr; ir++) // Set thymine methyl torsional variable
	{
		if( IsThymineRes(ir) && lock_bb[12*ir + 11] == 0 )
		{
			ki++;
			mnn_.var[ki] = bb_crd[12*ir + 11];
		}
	}

	if(IsSupHlxConstr()) 
	{	
		int i = symjm_.isur;
		if(i > 0) mnn_.var[i - 1] = sup_helix_rad;
		i = symjm_.isup;
		if(i > 0) mnn_.var[i - 1] = sup_helix_pit;
	}

	setbac_();
	microb_();

#endif

	return TRUE;
}

int
NuclAcidMod::BBOffset(int ir)
{
#if defined(INT_JUMNA)
	return ind_.kapt[ir];
#else
	return 1;
#endif
}

int
NuclAcidMod::IsRibRes(int ir)
{
#if defined(INT_JUMNA)
	return extjm_.ribose[ir];
#else
	return 1;
#endif
}

int 
NuclAcidMod::IsThymineRes(int ir)
{
#if defined(INT_JUMNA)
	return strjm_.lthy[ir];
#else
	return 1;
#endif
}

int 
NuclAcidMod::NumIndBBVarRes(int ir)
{
#if defined(INT_JUMNA)
	int in = strjm_.itr[ir];
	return flx_.kap[ 3*(in - 1) ];
#else
	return 1;
#endif
}

int 
NuclAcidMod::NumAllBBVarRes(int ir)
{
#if defined(INT_JUMNA)
	int in = strjm_.itr[ir];
	return flx_.kap[ 3*(in-1) + 1 ];
#else
	return 1;
#endif
}

int
NuclAcidMod::Is3endRes(int ir)
{
#if defined(INT_JUMNA)
	return AbsFun(ind_.ise[ir]);
#else
	return 1;
#endif
}

int 
NuclAcidMod::IsFstResInChain(int ir) 
{
	int ires = 0;
#if defined(INT_JUMNA)
	if( ind_.ise[ir] < 0) ires = 1;
#endif
	return ires;
}

int 
NuclAcidMod::IdxResInChain(int ir) 
{
#if defined(INT_JUMNA)
	return (strjm_.ilq[ir] - 1);
#else
	return 1;
#endif
}

int 
NuclAcidMod::ChainIdxOfRes(int ir)
{
#if defined(INT_JUMNA)
	return strjm_.ilq[ir + N2_J];
#else
	return 1;
#endif
}

int NuclAcidMod::IsSupHlxConstr()
{
	return( (sup_helix_rad > fabs(0.001)) );	
}

int NuclAcidMod::CalcPdistToShlxCnt()
{
	Vec3D cnt_shlx;
	cnt_shlx.SetX(sup_helix_rad);
    cnt_shlx.SetY(0.0);
	cnt_shlx.SetZ(0.0);
	
	if(p_dna_mol == NULL) return FALSE;
	
	HaResidue* pres;
	HaChain*   chain;

	ChainIteratorMolecule chitr(p_dna_mol);

	PrintLog(" P distances to Superhelix center \n"); 
	for(chain = chitr.GetFirstChain();chain; chain = chitr.GetNextChain())
	{
		ResidueIteratorChain ritr_ch(chain);
		for(pres = ritr_ch.GetFirstRes();pres; pres = ritr_ch.GetNextRes())
		{
			HaAtom* aptr = pres->GetAtomByName("P");
			if(aptr != NULL)
			{
				double dist = Vec3D::CalcDistance(aptr,&cnt_shlx,ANGSTROM_U);
				PrintLog("%s %9.4f \n",(aptr->GetRef(HaAtom::ATOMREF_NO_MOL)).c_str(),dist);
			}
		}
	}
	PrintLog(" \n"); 

	return TRUE;
}

int NuclAcidMod::GetNRes() const
{
	if(p_dna_mol == NULL) return 0;
	return( p_dna_mol->GetNRes());
}

int NuclAcidMod::CreateMolFromJumna()
{
#if defined(INT_JUMNA)

	char buf[256];
	int na = mrc_.kam;
	PrintLog(" na = %5d \n",na);
	int i,j;
	MolSet* pmset = GetMolSet();
	MolEditor* p_mol_editor = pmset->GetMolEditor(true);
	HaMolecule* pMol= pmset->CreateMolecule();
	if(pMol == NULL) return FALSE;

	p_dna_mol = pMol;
	pMol->SetObjName("DNA");

	int nchain = strjm_.nst;
	HaVec_int res_in_chain(nchain);
	for(i=0; i < nchain; i++)
		res_in_chain[i] = AbsFun(strjm_.idr[i]);

	int nres_old = -100; 
	HaResidue* pres;
	
	int ich = -1;

	char ch_symb = 'A' - 1;
	int mres_loc = 0;

	HaChain* chain;
	for(i=0; i < na; i++)
	{
		int nres = mrc_.nunit[i];

		if(nres != nres_old)
		{
			mres_loc++;
			if( nres == 1 || mres_loc > res_in_chain[ich])
			{
				ch_symb++;
                ich++;
				chain = pMol->AddChain(ch_symb);
				mres_loc = 1;
			}
			pres = chain->AddResidue(mres_loc);
			if(pres == NULL) return FALSE;
//			int if_bsat = mrc_.nuc[nres-1] + extjm_.iofs[nres];
//			strncpy(buf,&mrc_.munit[4*if_bsat],4);
//			buf[4] = 0;
//			std::string res_name = buf;
			std::string res_name = "";
			res_name += seq[ich][mres_loc-1];

			boost::trim(res_name);
			pres->SetName(res_name.c_str());
			nres_old = nres;
		}
		HaAtom* aptr = pres->AddNewAtom();
		aptr->SetElemNo(mrc_.imch[i]);
		strncpy(buf,&mrc_.mnam[4*i],4);
		buf[4] = 0;
		for(j=0; j < 4; j++)
		{
			if(buf[j] == '\'') buf[j] = 'X';
			if(buf[j] == '\"') buf[j] = 'Y';
		}
		std::string atname = buf;
		boost::trim(atname);
		aptr->SetName(atname.c_str());
		aptr->SetX(mrc_.corm[i] );
		aptr->SetY(mrc_.corm[N1_J+i] );
        aptr->SetZ(mrc_.corm[N1_J*2+i] );
		aptr->SetCharge(mrc_.dmon[i]);
	}
	SetCoordsFromJumna();
	p_mol_editor->CreateCovBonds(pMol);
	HaMolView* pView= pmset->GetActiveMolView();
	if(!pView)
		return TRUE;
	pView->ReDrawFlag |= RFInitial;
	pView->InitialTransform();
	pView->DefaultRepresentation();	
	pmset->RefreshAllViews();

#endif
	return TRUE;
}

int
NuclAcidMod::SaveConfig()
{
	int_4 zero = 0;
#if defined(INT_JUMNA)
	putbac_(&zero);
#endif
	return TRUE;
}

int
NuclAcidMod::RestoreConfig()
{
	int_4 one = 1;
#if defined(INT_JUMNA)
	putbac_(&one);
#endif
	return TRUE;
}

int 
NuclAcidMod::NumBBCoord()
{
	int nv = 0;

	if(p_dna_mol == NULL) return nv;
	int nr = p_dna_mol->GetNRes();

//	int ir;

//	for(ir=0; ir < nr; ir++)
//	{
//        int in = strjm_.itr[ir];
//		nv += flx_.kap[3*(in-1)];
//	}
//	PrintLog(" Computed number of Indep Backbone vars = %d, mnn_.nbac = %d \n",
//		       nv, mnn_.nbac);

#if defined(INT_JUMNA)
	nv = mnn_.ntba;
#endif
	return nv;
}

double
NuclAcidMod::GetAtomCrd(int at_num, int crd_num)
{
#if defined(INT_JUMNA)
	if(at_num > mrc_.kam || crd_num > 3) return 0.0;
	return mrc_.corm[ (crd_num - 1)*N1_J + (at_num - 1) ];
#else
	return 0.0;
#endif

}

double
NuclAcidMod::GetVarNum(int crd_num)
{
#if defined(INT_JUMNA)
	return mnn_.var[crd_num];
#else
	return 0.0;
#endif

}


int 
NuclAcidMod::NumFreeBBCoord()
{
	int nv = 0;
#if defined(INT_JUMNA)
	nv = mnn_.nbac;
#endif
	return nv;
}

int 
NuclAcidMod::IdxLastHelCoord()
{
	int nv = 0;
#if defined(INT_JUMNA)
	nv = mnn_.nthe;
#endif
	return nv;
}

int 
NuclAcidMod::IdxLastFreeHelCoord()
{
	int nv = 0;
#if defined(INT_JUMNA)
	nv = mnn_.nhel;
#endif
	return nv;
}

int 
NuclAcidMod::ReadAXEfile(const char* axe_file)
{
	char fname[40];

	strcpy_to_fort(fname,axe_file,40);
#if defined(INT_JUMNA)
	reset_(fname,40);
#endif
	return TRUE;
}


int 
NuclAcidMod::SaveAXEfile(const char* axe_file)
{
	char fname[40];
	strcpy_to_fort(fname,axe_file,40);
#if defined(INT_JUMNA)
	axeout_(fname,40);
#endif
	return TRUE;	
}
