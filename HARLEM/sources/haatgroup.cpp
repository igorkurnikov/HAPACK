/*! \file haatgroup.cpp

    Classes to define atomic group and residues in HARLEM.

    \author Igor Kurnikov  
    \date 1998-2002
*/

#define HAATGROUP_CPP

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdexcept>

#include <mpi.h>

#include <boost/algorithm/string.hpp>

#include "rapidxml.hpp"

#include "haatgroup.h"
#include "hamolecule.h"
#include "haresdb.h"
#include "hamolset.h"
#include "abstree.h"

#include "mm_force_field.h"

int HaResidue::RegisterResName( const std::string& res_name)
{
	StrIntMap::iterator itr = res_name_refno_map.find(res_name);
	if( itr != res_name_refno_map.end() ) return (*itr).second;
	ResNames.push_back(res_name);
	int refno = ResNames.size() - 1;
	res_name_refno_map[res_name] = refno;
	return refno;
}

void HaResidue::InitStdResNames()
{
	if (ResNames.size() > 20) return;
	ResNames.reserve(200);

	RegisterResName("ALA");  // 0 
	RegisterResName("GLY");  // 1
	RegisterResName("LEU");
	RegisterResName("SER");
	RegisterResName("VAL");
	RegisterResName("THR");  // 5
	RegisterResName("LYS");
	RegisterResName("ASP");
	RegisterResName("ILE");
	RegisterResName("ASN");
	RegisterResName("GLU");  // 10
	RegisterResName("PRO");
	RegisterResName("ARG");
	RegisterResName("PHE");
	RegisterResName("GLN");
	RegisterResName("TYR");  // 15
	RegisterResName("HIS");
	RegisterResName("CYS");
	RegisterResName("MET");
	RegisterResName("TRP");

	RegisterResName("ASX");  // 20
	RegisterResName("GLX");
	RegisterResName("PCA");
	RegisterResName("HYP");

	//  DNA Nucleotides
	RegisterResName("A");    // 24
	RegisterResName("C");    // 25
	RegisterResName("G");    // 26
	RegisterResName("T");    // 27

	// RNA Nucleotides  
	RegisterResName("U");    // 28
	RegisterResName("+U");
	RegisterResName("I");    // 30
	RegisterResName("1MA");
	RegisterResName("5MC");
	RegisterResName("OMC");
	RegisterResName("1MG");
	RegisterResName("2MG");  // 35
	RegisterResName("M2G");
	RegisterResName("7MG");
	RegisterResName("OMG");
	RegisterResName("YG");
	RegisterResName("H2U");  // 40
	RegisterResName("5MU");
	RegisterResName("PSU");

	// Lipids
	RegisterResName("HPC");  // 43
	//  Miscellaneous  
	RegisterResName("UNK");  // 44
	RegisterResName("ACE");  // 45
	RegisterResName("FOR");  // 46
	RegisterResName("HOH");
	RegisterResName("DOD");
	RegisterResName("SO4");
	RegisterResName("PO4");  // 50
	RegisterResName("NAD");
	RegisterResName("COA");
	RegisterResName("NAP");
	RegisterResName("TRM");
}

void HaResidue::InitResSynonym()
{
	if (ResSynonym_to_std.size() > 10) return;
	ResSynonym_to_std["ADE"] = "A";    /*  Adenosine   */
	ResSynonym_to_std["GUA"] = "G";    /*  Guanosine   */
	ResSynonym_to_std["THY"] = "T#D";    /*  Thymidine   */
	ResSynonym_to_std["URI"] = "U";    /*  Uridine     */
	ResSynonym_to_std["CYT"] = "C";    /*  Cytosine    */

	ResSynonym_to_std["A5"] = "A#R5";
	ResSynonym_to_std["A3"] = "A#R3";
	ResSynonym_to_std["AN"] = "A#RN";
	ResSynonym_to_std["DA"] = "A#D";
	ResSynonym_to_std["DA3"] = "A#D3";
	ResSynonym_to_std["DA5"] = "A#D5";
	ResSynonym_to_std["DAN"] = "A#DN";

	ResSynonym_to_std["G5"] = "G#R5";
	ResSynonym_to_std["G3"] = "G#R3";
	ResSynonym_to_std["GN"] = "G#RN";
	ResSynonym_to_std["DG"] = "G#D";
	ResSynonym_to_std["DG3"] = "G#D3";
	ResSynonym_to_std["DG5"] = "G#D5";
	ResSynonym_to_std["DGN"] = "G#DN";

	ResSynonym_to_std["C5"] = "C#R5";
	ResSynonym_to_std["C3"] = "C#R3";
	ResSynonym_to_std["CN"] = "C#RN";
	ResSynonym_to_std["DC"] = "C#D";
	ResSynonym_to_std["DC3"] = "C#D3";
	ResSynonym_to_std["DC5"] = "C#D5";
	ResSynonym_to_std["DCN"] = "C#DN";

	ResSynonym_to_std["DT"]  = "T#D";
	ResSynonym_to_std["DT3"] = "T#D3";
	ResSynonym_to_std["DT5"] = "T#D5";
	ResSynonym_to_std["DTN"] = "T#DN";

	ResSynonym_to_std["U5"] = "U#R5";
	ResSynonym_to_std["U3"] = "U#R3";
	ResSynonym_to_std["UN"] = "U#RN";
	
	ResSynonym_to_std["NALA"] = "ALA#NT";  /*  N-terminal ALA  */
	ResSynonym_to_std["NARG"] = "ARG#NT";  /*  N-terminal ARG  */
	ResSynonym_to_std["NASN"] = "ASN#NT";  /*  N-terminal ASN  */
	ResSynonym_to_std["NGLN"] = "GLN#NT";  /*  N-terminal GLN  */
	ResSynonym_to_std["NASP"] = "ASP#NT";  /*  N-terminal ASP  */
	ResSynonym_to_std["NHIS"] = "HIS#NT";  /*  N-terminal HIS DELTA  */
	ResSynonym_to_std["NHIE"] = "HIS#EPSILON_NT";  /*  N-terminal HIS EPSILON  */
	ResSynonym_to_std["NHID"] = "HIS#NT";  /*  N-terminal HIS DELTA  */
	ResSynonym_to_std["NHSE"] = "HIS#EPSILON_NT";  /*  N-terminal HIS EPSILON  */
	ResSynonym_to_std["NHSD"] = "HIS#NT";  /*  N-terminal HIS DELTA  */
	ResSynonym_to_std["NHIP"] = "HIS#PROT_NT";  /*  N-terminal HIS PROTONATED  */
	ResSynonym_to_std["NPRO"] = "PRO#NT";  /*  N-terminal PRO   */
	ResSynonym_to_std["NTRP"] = "TRP#NT";  /*  N-terminal TRP   */
	ResSynonym_to_std["NCYS"] = "CYS#NT";  /*  N-terminal CYS   */
	ResSynonym_to_std["NCYX"] = "CYS#UNPROT_NT";  /*  N-terminal CYS UNPROTONATED   */
	ResSynonym_to_std["NGLY"] = "GLY#NT";  /*  N-terminal GLY   */
	ResSynonym_to_std["NILE"] = "ILE#NT";  /*  N-terminal ILE   */
	ResSynonym_to_std["NLEU"] = "LEU#NT";  /*  N-terminal LEU   */
	ResSynonym_to_std["NLYS"] = "LYS#NT";  /*  N-terminal LYS   */
	ResSynonym_to_std["NMET"] = "MET#NT";  /*  N-terminal MET   */
	ResSynonym_to_std["NPHE"] = "PHE#NT";  /*  N-terminal PHE   */
	ResSynonym_to_std["NSER"] = "SER#NT";  /*  N-terminal SER   */
	ResSynonym_to_std["NTYR"] = "TYR#NT";  /*  N-terminal TYR   */
	ResSynonym_to_std["NVAL"] = "VAL#NT";  /*  N-terminal VAL   */
	ResSynonym_to_std["NTHR"] = "THR#NT";  /*  N-terminal THR   */
	ResSynonym_to_std["NGLU"] = "GLU#NT";  /*  N-terminal THR   */

	ResSynonym_to_std["CALA"] = "ALA#CT";  /*  C-terminal ALA  */
	ResSynonym_to_std["CARG"] = "ARG#CT";  /*  C-terminal ARG  */
	ResSynonym_to_std["CASN"] = "ASN#CT";  /*  C-terminal ASN  */
	ResSynonym_to_std["CGLN"] = "GLN#CT";  /*  C-terminal GLN  */
	ResSynonym_to_std["CASP"] = "ASP#CT";  /*  C-terminal ASP  */
	ResSynonym_to_std["CHIS"] = "HIS#CT";  /*  C-terminal HIS DELTA  */
	ResSynonym_to_std["CHIE"] = "HIS#EPSILON_CT";  /*  C-terminal HIS EPSILON  */
	ResSynonym_to_std["CHID"] = "HIS#CT";  /*  C-terminal HIS DELTA  */
	ResSynonym_to_std["CHSE"] = "HIS#EPSILON_CT";  /*  C-terminal HIS EPSILON  */
	ResSynonym_to_std["CHSD"] = "HIS#CT";  /*  C-terminal HIS DELTA  */
	ResSynonym_to_std["CHIP"] = "HIS#PROT_CT";  /*  C-terminal HIS PROTONATED  */
	ResSynonym_to_std["CPRO"] = "PRO#CT";  /*  C-terminal PRO   */
	ResSynonym_to_std["CTRP"] = "TRP#CT";  /*  C-terminal TRP   */
	ResSynonym_to_std["CCYS"] = "CYS#CT";  /*  C-terminal CYS   */
	ResSynonym_to_std["CCYX"] = "CYS#UNPROT_CT";  /*  C-terminal CYS UNPROTONATED   */
	ResSynonym_to_std["CGLY"] = "GLY#CT";  /*  C-terminal GLY   */
	ResSynonym_to_std["CILE"] = "ILE#CT";  /*  C-terminal ILE   */
	ResSynonym_to_std["CLEU"] = "LEU#CT";  /*  C-terminal LEU   */
	ResSynonym_to_std["CLYS"] = "LYS#CT";  /*  C-terminal LYS   */
	ResSynonym_to_std["CMET"] = "MET#CT";  /*  C-terminal MET   */
	ResSynonym_to_std["CPHE"] = "PHE#CT";  /*  C-terminal PHE   */
	ResSynonym_to_std["CSER"] = "SER#CT";  /*  C-terminal SER   */
	ResSynonym_to_std["CTYR"] = "TYR#CT";  /*  C-terminal TYR   */
	ResSynonym_to_std["CVAL"] = "VAL#CT";  /*  C-terminal VAL   */
	ResSynonym_to_std["CTHR"] = "THR#CT";  /*  C-terminal THR   */
	ResSynonym_to_std["CGLU"] = "GLU#CT";  /*  C-terminal THR   */

	ResSynonym_to_std["HID"] = "HIS";  /*  histidine - hydrogen at delta position  */
	ResSynonym_to_std["HIE"] = "HIS#EPSILON";  /*  histidine - hydrogen at epsilon position  */
	ResSynonym_to_std["HIP"] = "HIS#PROT";  /*  histidine - protonated  */
	ResSynonym_to_std["CYX"] = "CYS#UNPROT";  /*  Cystine unprotonated*/

	// Make a reverse map from standard to AMBER residue names:
	for ( const auto& kv : ResSynonym_to_std ) 
	{
		ResSynonym_std_to_AMBER[kv.second] = kv.first;
	}

	// Next conversions are not one to one

	ResSynonym_to_std["HSD"] = "HIS";  /*  histidine - hydrogen at delta position  */
	ResSynonym_to_std["HSE"] = "HIS#EPSILON";  /*  histidine - hydrogen at epsilon position  */
	ResSynonym_to_std["HSP"] = "HIS#PROT";  /*  histidine - protonated  */

	ResSynonym_to_std["CPR"] = "PRO";  /*  Cis-proline */
	ResSynonym_to_std["TRY"] = "TRP";  /*  Cis-proline */
	ResSynonym_to_std["CSH"] = "CYS";  /*  Cystine     */
	ResSynonym_to_std["CSM"] = "CYS";  /*  Cystine     */
	ResSynonym_to_std["CYH"] = "CYS";  /*  Cystine     */

	ResSynonym_to_std["D2O"] = "DOD";  /*  Heavy Water */
	ResSynonym_to_std["H2O"] = "HOH";  /*  Water       */
	ResSynonym_to_std["SOL"] = "HOH";  /*  Water       */
	ResSynonym_to_std["TIP"] = "HOH";  /*  Water       */
	ResSynonym_to_std["WAT"] = "HOH";  /*  Water       */

	ResSynonym_std_to_AMBER["HOH"] = "WAT";

}

std::vector<std::string> HaResidue::ResNames;
StrStrMap  HaResidue::ResSynonym_to_std;
StrStrMap  HaResidue::ResSynonym_std_to_AMBER;

StrIntMap  HaResidue::res_name_refno_map;

AtomGroup::AtomGroup(AtomExpr* expr, MolSet* pmset)
{ 
	SetFromExpr(expr, pmset);	
}

AtomGroup::AtomGroup(const AtomGroup& ref_atset): vector<HaAtom*>(ref_atset)
{
	id = ref_atset.id;
}


AtomGroup::~AtomGroup()
{
	
}


bool AtomGroup::InsertAtom(HaAtom* aptr)
{
	if(aptr == NULL)
		return false;
	this->push_back(aptr); 
	return true;	
}


bool AtomGroup::DeleteAtom(HaAtom* aptr)
{
	if(aptr == NULL)
		return false;

	vector<HaAtom*>::iterator aitr;

	for(aitr = this->begin(); aitr != this->end();)
	{
		if( (*aitr) == aptr)
		{
			aitr = this->erase(aitr);
		}
		else
		{
			aitr++;
		}
	} 
	return true;
}

HaAtom* AtomGroup::GetAtomByName(const std::string& at_id)
{   
   AtomIteratorAtomGroup aitr(this);
   HaAtom* aptr;
   for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
   {
	   if(aptr->GetName() == boost::trim_copy(boost::to_upper_copy(at_id)) ) return aptr;
   }
   return NULL;
}

int AtomGroup::DelSelAtoms()
{
	vector<HaAtom*>::iterator aitr;

	int ndel= 0; 

	for(aitr = this->begin(); aitr != this->end();)
	{
		if( (*aitr)->Selected())
		{
			aitr = this->erase(aitr);
			ndel++;
		}
		else
		{
			aitr++;
		}
	} 
	return ndel;
}

int AtomGroup::DeleteAtoms(const PtrSet& ptr_set)
{
	vector<HaAtom*>::iterator aitr;

	int ndel= 0; 

	for(aitr = this->begin(); aitr != this->end();)
	{
		if( ptr_set.IsMember(*aitr))
		{
			aitr = this->erase(aitr);
			ndel++;
		}
		else
		{
			aitr++;
		}
	} 
	return ndel;	
}



int AtomGroup::IsMember(const HaAtom* aptr) const
{
	vector<HaAtom*>::const_iterator itr;
	
	for(itr= begin(); itr != end(); itr++)
	{
		if((*itr) == aptr)
			return( True);
	}
    return( False );
}

void AtomGroup::SetFromExpr(AtomExpr* expr, MolSet* pmset)
{	
	clear();
	HaAtom* aptr;
    if( pmset )
	{
	    AtomIteratorMolSet aitr(pmset);
	    for(aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom())
		{
			if( expr->EvaluateExprFor(aptr) )
				this->push_back(aptr);
		}
	}
}

void AtomGroup::SetFromExprStr(const char* expr_str, MolSet* pmset )
{
	clear();
	AtomExpr* p_expr;
	if( (p_expr = AtomExpr::ParseExpression(expr_str,pmset)) != NULL )
	{   
		SetFromExpr(p_expr, pmset);
		delete p_expr;
	}
}

void AtomGroup::AddFromExpr(AtomExpr* expr, MolSet* pmset)
{	
	if( !pmset) return;
	HaAtom* aptr;
	std::set<HaAtom*,less<HaAtom*> > old_atoms;
	AtomIteratorAtomGroup aitr_g(this);

	for(aptr = aitr_g.GetFirstAtom(); aptr; aptr = aitr_g.GetNextAtom())
	{
		old_atoms.insert(aptr);
	}

	AtomIteratorMolSet aitr(pmset);
	for(aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom())
	{
		if( old_atoms.count(aptr) > 0 ) continue;
		if( expr->EvaluateExprFor(aptr) ) this->push_back(aptr);
	}
}

void AtomGroup::AddFromExprStr(const char* expr_str, MolSet* pmset )
{
	AtomExpr* p_expr;
	if( (p_expr = AtomExpr::ParseExpression(expr_str,pmset)) != NULL )
	{      
		AddFromExpr(p_expr, pmset);
		delete p_expr;
	}
}

void AtomGroup::DeleteAtomsExpr(AtomExpr* expr, MolSet* pmset)
{	
	if( !pmset) return;
	HaAtom* aptr;
	std::set<HaAtom*,less<HaAtom*> > old_atoms;
	AtomIteratorAtomGroup aitr_g(this);

	for(aptr = aitr_g.GetFirstAtom(); aptr; aptr = aitr_g.GetNextAtom())
	{
		old_atoms.insert(aptr);
	}
	clear();

	AtomIteratorMolSet aitr(pmset);
	for(aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom())
	{
		if( old_atoms.count(aptr) > 0 && !expr->EvaluateExprFor(aptr) )
			this->push_back(aptr);
	}
}

void AtomGroup::DeleteAtomsExprStr(const char* expr_str, MolSet* pmset )
{
	AtomExpr* p_expr;
	if( (p_expr = AtomExpr::ParseExpression(expr_str,pmset)) != NULL )
	{   
		DeleteAtomsExpr(p_expr, pmset);
		delete p_expr;
	}
}

void AtomGroup::KeepOnlyAtomsExpr(AtomExpr* expr, MolSet* pmset)
{	
	if( !pmset) return;
	HaAtom* aptr;
	std::set<HaAtom*,less<HaAtom*> > old_atoms;
	AtomIteratorAtomGroup aitr_g(this);

	for(aptr = aitr_g.GetFirstAtom(); aptr; aptr = aitr_g.GetNextAtom())
	{
		old_atoms.insert(aptr);
	}
	clear();

	AtomIteratorMolSet aitr(pmset);
	for(aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom())
	{
		if( old_atoms.count(aptr) > 0 && expr->EvaluateExprFor(aptr) )
			this->push_back(aptr);
	}
}

void AtomGroup::KeepOnlyAtomsExprStr(const char* expr_str, MolSet* pmset )
{
	AtomExpr* p_expr;
	if( (p_expr = AtomExpr::ParseExpression(expr_str,pmset)) != NULL )
	{      
		KeepOnlyAtomsExpr(p_expr, pmset);
		delete p_expr;
	}
}


AtomIterator* AtomGroup::GetAtomIteratorPtr() 
{ 
	return new AtomIteratorAtomGroup(this); 
}

int
AtomGroup::GetNAtoms() const
{
	return size();
}

PointIterator* 
AtomGroup::GetPointIteratorPtr() 
{ 
	return new AtomIteratorAtomGroup(this); 
}

PointIterator_const* AtomGroup::GetPointIteratorPtr() const 
{ 
	return new AtomIteratorAtomGroup_const(this); 
}

void AtomGroup::SetID(const std::string& new_id)
{
	id = new_id;
	boost::trim(id);
	boost::to_upper(id);
}

const char* AtomGroup::GetID() const
{
	return(id.c_str());
}


ChemGroup::ChemGroup()
{
	phost_mset= NULL;
	id="";
	protect=1.0;
}


ChemGroup::ChemGroup(MolSet* new_phost_mset, const char* new_id)
{
	phost_mset= new_phost_mset;
	id=new_id;
	protect=1.0;
}


ChemGroup::~ChemGroup()
{

}

 
bool ChemGroup::FillRef(char* buf) const
{
	sprintf(buf,"chemgrp=%s",id.c_str());
	return true;
}

std::string ChemGroup::GetIDFromRef(const std::string& buf)
{
	std::string tmp_id=buf.substr(8);
	boost::trim(tmp_id);
	boost::to_upper(tmp_id);

	return tmp_id;
}
 
bool ChemGroup::Print_info(ostream &sout, const int level) const
{
	sout << " Atom Group with id " << id << endl;
	sout << " contains  " << size() << "  atoms" << endl;
	vector<HaAtom*>::const_iterator caitr;
	for ( caitr= begin(); caitr != end(); caitr++);
	{
		char buf[256];
		(*caitr)->FillRef(buf);
		sout << buf << endl;
	}
	return true;
}

double ChemGroup::GetProtect() const
// return screening factor 
{ 
	return protect;
}

bool ChemGroup::SetProtect(const double new_protect) 
		  // return screening factor 
{ 
	if( (new_protect >= 0.0) && (new_protect <= 1.0) )
	{
		protect=new_protect;
		return true;
	}
	return false;	
}



HaResidue::HaResidue()
{
	Clear();
}


//HaResidue::HaResidue(const HaResidue& ref_res):
//AtomGroup(ref_res)
//{
//	ResetData();
//	SetParamFrom(ref_res);
//	phost_ch = ref_res.phost_ch;
//}



HaResidue::HaResidue(HaChain* new_phost)
{
	Clear();
	phost_ch=new_phost;
//	this->reserve(10);
}


HaResidue::~HaResidue()
{

}

void HaResidue::Clear()
{
	phost_ch=NULL;                        
    serno=0;                
    width=0.0;                 
    col1=0;                 
    col2=0;                 
	insert=' ';                
	refno=0;                
	struc=0;                
	flag=0; 
}

bool HaResidue::SetParamFrom(const HaResidue& res_ref)
{
    serno=  res_ref.serno;                
    width=  res_ref.width;                 
    col1 =  res_ref.col1;                 
    col2 =  res_ref.col2;                 
	insert= res_ref.insert;                
	refno = res_ref.refno;                
	struc = res_ref.struc;                
	flag  = res_ref.flag;
    SetNameModifier(res_ref.NameModifier.c_str());
    
	return true;
}

HaMolecule* HaResidue::GetHostMol()
{
	return phost_ch->GetHostMol();
}

const HaMolecule* HaResidue::GetHostMol() const
{
	return phost_ch->GetHostMol();
}

MolSet* HaResidue::GetHostMolSet()
{
	return (phost_ch->GetHostMol())->GetHostMolSet();
}

const MolSet* HaResidue::GetHostMolSet() const
{
	return (phost_ch->GetHostMol())->GetHostMolSet();
}

int HaResidue::GetNAtomsNonProxy() const
{
	int na = this->size();
	int i;
	int n_nproxy = 0;
	for(i = 0; i < na; i++)
	{
		const HaAtom* aptr = (*this)[i];
		if( !aptr->IsProxy() ) n_nproxy++;
	}
	return n_nproxy;
}


HaResidue* HaResidue::GetPrevResInChain()
{
	std::multimap<int, HaResidue*, less<int> >::iterator ritr;

	ritr = phost_ch->res_map.find( this->GetSerNo() );
	if(ritr == phost_ch->res_map.begin() || ritr == phost_ch->res_map.end()) 
		return NULL;
	ritr--;
	return( (*ritr).second );
}

HaResidue* HaResidue::GetPrevResInChain() const
{
	std::multimap<int, HaResidue*, less<int> >::const_iterator ritr;

	ritr = phost_ch->res_map.find( this->GetSerNo() );
	if(ritr == phost_ch->res_map.begin() || ritr == phost_ch->res_map.end()) 
		return NULL;
	ritr--;
	return( (*ritr).second );
}

HaResidue* HaResidue::GetNextResInChain()
{
	std::multimap<int, HaResidue*, less<int> >::iterator ritr;
	ritr = phost_ch->res_map.find( this->GetSerNo() );

	if(ritr == phost_ch->res_map.end()) return NULL;
	ritr++;
	if(ritr == phost_ch->res_map.end()) return NULL;

	return( (*ritr).second );
}

HaResidue* HaResidue::GetNextResInChain() const
{
	std::multimap<int, HaResidue*, less<int> >::const_iterator ritr;
	ritr = phost_ch->res_map.find( this->GetSerNo() );

	if(ritr == phost_ch->res_map.end()) return NULL;
	ritr++;
	if(ritr == phost_ch->res_map.end()) return NULL;

	return( (*ritr).second );
}


bool HaResidue::IsBonded(HaResidue* res2) 
{
	if(res2 == NULL) return false;
	vector<HaAtom*>::const_iterator aitr1, aitr2;
	for(aitr1= this->begin(); aitr1 != this->end(); aitr1++)
	{
		for(aitr2= res2->begin(); aitr2 != res2->end(); aitr2++)
		{
			if( (*aitr1)->IsBonded(*(*aitr2)) ) return true;
		}	
	}
	return false;
}


int HaResidue::HasBackBHBond(HaResidue* res2)
{
	HaAtom* pn1 = GetAtomByName("N");
	HaAtom* po1 = GetAtomByName("O");
	HaAtom* pn2 = res2->GetAtomByName("N");
	HaAtom* po2 = res2->GetAtomByName("O");

	if(pn1 == NULL || po1 == NULL || pn2 == NULL || po2 == NULL)
	{
		return FALSE;
	}
	
	MolSet* pmset = this->GetHostMolSet();
	if(pmset->AreHBonded(pn1,po2) ) return TRUE;
	if(pmset->AreHBonded(pn2,po1) ) return TRUE;

#if 0 // Check for hydrogen bond without attaching hydrogens
    int i;
	for(i = 0; i < 2; i++)
	{
	   Vec3D hn;
       HaAtom* pn = pn1;
	   HaAtom* po = po2;
	   if(i == 1) 
	   {
		   pn = pn2;
		   po = po1;
	   }

	   int hn_set = FALSE;
       HaAtom* pca = NULL;
  	   HaAtom* pc  = NULL;
       
       AtomGroup bnd_at;
       pn->GetBondedAtoms(bnd_at);

	   HaAtom* aptr;
	   AtomIteratorAtomGroup aitr(&bnd_at);
	   for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	   {
          if( !strcmp(aptr->GetName(),"H"))
		  {
             hn.SetX(aptr->GetX()); hn.SetY(aptr->GetY()); hn.SetZ(aptr->GetZ());
             hn_set = TRUE;
		     break;
		  }
		  else if(!strcmp(aptr->GetName(),"CA"))
		  {
             pca = aptr;
		  }
		  else
		  {
             pc = aptr;
		  }
	   }
	   if(!hn_set && (pca != NULL && pc != NULL) )
	   {
		  Vec3D::SetAtomPos( &hn, pn, pca, pc, 1.08, PI*0.6667, PI);
		  hn_set = TRUE;
	   }

	   if(hn_set)
	   {
		   double d_no    = Vec3D::CalcDistance(pn,po,ANGSTROM_U);
		   double ang_nho  = Vec3D::CalcAngle(pn,&hn,po)*RAD_TO_DEG;
		  if( d_no < 4.2 && ang_nho > 100.0)
		  {
			 return TRUE;
		  }
	   }
	}
#endif

	return FALSE;
}



bool HaResidue::FillRef(char* buf,int mode) const
{
	const char* str;
	str=GetName();
		
	int i,j;

	j=0;

	const HaMolecule* phost_mol= GetHostMol();

	if( mode == 0)
	{
		j+= sprintf(buf+j,"$%s$",phost_mol->GetObjName());
	}

	for(i=0; i<4; i++)
	{
		if(str[i] == 0) break;
		if( str[i]!=' ') j += sprintf(buf+j,"%c", str[i]);
	}
	j += sprintf(buf+j,"%d",serno);

	if( phost_ch->ident!=' ' )
	{   
	    j += sprintf(buf+j, ":%c",phost_ch->ident);
	}
	return true;
}

std::string HaResidue::GetRef() const
{
	char buf[256];
	FillRef(buf);
	return buf;
}


std::string HaResidue::GetFullName() const
{
	std::string name(GetName());
	if(!NameModifier.empty()) name+= "#" + NameModifier;
	return name;
}


std::string HaResidue::GetResNameFromFullName(const char* res_full_name)
{
	std::string fname_str = res_full_name;
	size_t sep_pos = fname_str.find_first_of("#");
	if( sep_pos == std::string::npos ) { return res_full_name; }
	return fname_str.substr(0,sep_pos);
}

bool HaResidue::SetUniqueAtomNames()
{
	HaAtom* aptr;
	bool uniq_names = true;

	using std::string;
	using std::set;
	using std::multimap;
	using std::map;

	int iadd=1;
	multimap<string, HaAtom*, less<string> > atname_at_map;
	multimap<string, HaAtom*, less<string> >::iterator aa_map_itr;

	set<string, less<string> > atname_not_unique;
	set<string, less<string> > atname_all;

	AtomIteratorAtomGroup aitr(this);
	for(aptr= aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		aa_map_itr= atname_at_map.find( aptr->GetName() );
		if(aa_map_itr != atname_at_map.end())
		{
			uniq_names= false;
			atname_not_unique.insert( aptr->GetName() );
		}
		atname_all.insert( aptr->GetName() );
		atname_at_map.insert( pair<string,HaAtom*>( aptr->GetName(), aptr) );
	}
	if(uniq_names) return true;
	
	char buf[256];
	FillRef(buf);
	
	cerr << " Warning: in HaResidue::SetUniqueAtomNames() " << endl;
	cerr << " Atom Names in the Residue " << buf << " are not unique " << endl;
	cerr << " and will be changed " << endl;

	char lbi[4];
	set<string, less<string> >::iterator atname_itr;

	multimap<string,HaAtom*, less<string> >::iterator aa_map_itr1;
	multimap<string,HaAtom*, less<string> >::iterator aa_map_itr2;

	for(atname_itr = atname_not_unique.begin(); atname_itr != atname_not_unique.end(); atname_itr++)
	{
		aa_map_itr1 = atname_at_map.lower_bound(*atname_itr);
		aa_map_itr2 = atname_at_map.upper_bound(*atname_itr);
		aa_map_itr1++;
		if( aa_map_itr1 == aa_map_itr2) continue;
		
		for( aa_map_itr = aa_map_itr1; aa_map_itr != aa_map_itr2; aa_map_itr++)
		{
			std::string el_name = ((*aa_map_itr).second)->GetStdSymbol();
			std::string name= ((*aa_map_itr).second)->GetStdSymbol();
			for( ;iadd < 1000; iadd++)
			{
				if(iadd < 10)
					sprintf(lbi,"%1d",iadd);
				else if (iadd >= 10 && iadd < 100)
					sprintf(lbi,"%2d",iadd);
//				else if (iadd >= 100 && iadd < 1000 && (name.size() < 2) )
//					sprintf(lbi,"%3d",iadd);
				else
				{
					return false;
				}
				
				std::string name = el_name + lbi;
				if( atname_all.find(name) == atname_all.end())
				{
					((*aa_map_itr).second)->SetName(name);
					atname_all.insert(name);
					break;
				}			
			}
		}
	}
	return true;
}

std::string HaResidue::GetUniqueAtomName(int elem_no)
{
	set<std::string> at_names;
	HaAtom* aptr;
	std::string atname("");
	AtomIteratorAtomGroup aitr(this);
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		at_names.insert(aptr->GetName());
	}

    std::string at_symb = HaAtom::GetStdSymbolElem(elem_no);

	int max = 1000;
	if(at_symb.length() == 2)
		max = 100;

	int i;
	char buf[10];
	for(i=1; i < max; i++)
	{
        if(i<10)
		{
			sprintf(buf,"%s%1d",at_symb.c_str(),i);
		}
        if(i<100 && i > 9)
		{
			sprintf(buf,"%s%2d",at_symb.c_str(),i);
		}
        if(i<1000 && i > 99)
		{
			sprintf(buf,"%s%3d",at_symb.c_str(),i);
		}

        atname = buf;

		if(at_names.find(atname) == at_names.end())
			return atname;
	}

	return "";
}


const char* HaResidue::GetResNameInTable(const int j)
{
	if(j >= ResNames.size()) return "";
	return( ResNames[j].c_str());
}


const char* HaResidue::GetName() const
{
	return(ResNames[refno].c_str());
}

void HaResidue::SetName(const std::string& res_name_par, int convert_res_names )
{
	std::string new_name = boost::trim_copy(res_name_par);
	boost::to_upper( new_name );

	if( convert_res_names && ResSynonym_to_std.count(new_name.c_str()) > 0 )
	{
		std::string fname = ResSynonym_to_std.GetVal( new_name.c_str() );
		
		size_t isep = fname.find('#');
		if( isep == std::string::npos )
		{
			new_name = fname;
		}
		else
		{
			new_name = fname.substr(0,isep);
            NameModifier = fname.substr(isep+1);
		}
	}

	bool contain_digit=false;

	int i;
	for( i=0; i < new_name.size(); i++)
	{
		if( isdigit( new_name[i] ) )
		{
			contain_digit= true;
			new_name[i]= 'A' + (new_name[i] - '0');
		}
	}
	
	if(contain_digit)
	{
		PrintLog(" Warning in HaResidue::SetName() \n"); 
		PrintLog(" Residue name %s contain a digit \n", res_name_par.c_str() );
		PrintLog(" Has been changed to %s \n", new_name.c_str()); 
	}
 
	if( new_name.empty()) new_name = "MOL";
	refno = RegisterResName(new_name);
}

bool HaResidue::SplitResidue()
{
	int idx=0;
	bool is_large = false;
	vector<HaAtom*>::iterator aitr;
	for(aitr= begin(); aitr != end(); aitr++)
	{
		idx++;
		if(idx >= 100) 
		{	
			is_large = true;
			break;
		}
	}

	if(!is_large) 
	{
		return true;
	}

	bool create_new_res = true;
	HaResidue* cur_res = NULL;
	for(; aitr != end(); )
	{
		if(create_new_res)
		{
			idx = 1;
			int res_serno= phost_ch->GetUniqResSerNo(0);
			cur_res=phost_ch->AddResidue(res_serno);
			if(cur_res == NULL) return false;
			cur_res->SetName(this->GetName());
			create_new_res = false;
		}	
		cur_res->InsertAtom(*aitr);
		(*aitr)->SetHostRes(cur_res);

		aitr = this->erase(aitr);
		idx++;
		if(idx >= 100) 
		{
			create_new_res = true;
		}
	}
	return true;
}

bool AtomGroup::SwapAtoms(int i1, int i2)
{
	int na = this->GetNAtoms();
	if (i1 < 0 || i2 < 0 || i1 >= na || i2 >= na) return false;
	HaAtom* pat1 = (*this)[i1];
	(*this)[i1] = (*this)[i2];
	(*this)[i2] = pat1;
	return true;
}

int AtomGroup::HasSelectedAtoms()
{
	AtomIteratorAtomGroup aitr(this);
	HaAtom* aptr;
	
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if(aptr->Selected()) return TRUE;
	}
	return FALSE;
}

void AtomGroup::SelectAtomsAll()
{
	AtomIteratorAtomGroup aitr(this);
	HaAtom* aptr;	
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		aptr->Select();
	}
}

AtomIntMap HaResidue::GetAtomSeqNumMap()
{
	AtomIntMap at_seq_num_map;

	HaAtom* aptr = nullptr;
	AtomIteratorAtomGroup aitr(this);

	int i = 0;

	for (aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		at_seq_num_map[aptr] = i;
		i++;
	}
	return at_seq_num_map;
}

CAtomIntMap HaResidue::GetAtomSeqNumMap() const
{
	CAtomIntMap at_seq_num_map_loc;

	const HaAtom* aptr = nullptr;
	AtomIteratorAtomGroup_const aitr(this);
	int i = 0;

	for (aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		at_seq_num_map_loc[aptr] = i;
		i++;
	}
	return at_seq_num_map_loc;
}


bool HaResidue::IsAmino() const
{
	return( refno <= 23 );
}

bool HaResidue::IsAminoNucleo() const
{
	return( refno <= 42);
}

bool HaResidue::IsNucleo() const
{
	return( refno >= 24 && refno <=42 );
}

bool HaResidue::IsProtein() const
{
	return( refno <= 23  || ( refno >=43  && refno <=46 ) );
}

bool HaResidue::IsDNA() const
{
	return( refno >= 24 && refno <= 27 );
}

bool HaResidue::IsSolvent() const
{
	return( refno >=46 && refno <=49 );
}

bool HaResidue::IsWater() const
{
	return( refno == 46  || refno == 47 );
}

bool HaResidue::IsIon() const
{
	return( refno ==48  || refno == 49);
}

bool HaResidue::IsPyrimidine() const  
{
	return( IsCytosine() || IsThymine());
}

bool HaResidue::IsPurine() const
{
	return(IsAdenine() || IsGuanine());
}

bool HaResidue::IsRNA() const
{
         return(IsNucleo() && !IsThymine());
}

bool HaResidue::IsProline() const
{
	return(refno == 11);
}

bool HaResidue::IsHistidine() const
{
    return( refno == 16 );
}

bool HaResidue::IsCysteine() const
{
    return( refno ==17 );
}

bool HaResidue::IsAdenine() const
{
	return( refno == 24);
}

bool HaResidue::IsCytosine() const
{
    return ( refno ==25 );
}

bool HaResidue::IsGuanine() const
{
	return( refno == 26);
}

bool HaResidue::IsThymine() const
{
	return( refno==27);
}

bool HaResidue::IsCoenzyme() const
{
	return( refno >= 50 && refno <= 52 );
}

bool HaResidue::IsTerm() const
{
	return( stricmp_trunc( GetName(), "TRM") == 0 );
}

int HaResidue::CalcStdCrdSys(int fit_std_geom)
{
	if(this->IsNucleo())
	{
		HaAtom* ptr_c1p = NULL;
			
		ptr_c1p = GetAtomByName("C1X");
        if(ptr_c1p == NULL) ptr_c1p = GetAtomByName("C1'");
        if(ptr_c1p == NULL) ptr_c1p = GetAtomByName("C1*");
		if(ptr_c1p == NULL) return FALSE;
		
		HaAtom* ptr1 = NULL;
		if( this->IsPurine()) ptr1 = GetAtomByName("N9");
        if( this->IsPyrimidine()) ptr1 = GetAtomByName("N1");
		if(ptr1 == NULL) return FALSE;

		HaAtom* ptr2 = NULL;
		if( this->IsPurine()) ptr2 = GetAtomByName("C4");
        if( this->IsPyrimidine()) ptr2 = GetAtomByName("C2");
		if(ptr2 == NULL) return FALSE;

		HaAtom* ptr3 = NULL;
		if( this->IsPurine()) ptr3 = GetAtomByName("C8");
        if( this->IsPyrimidine()) ptr3 = GetAtomByName("C6");
		if(ptr3 == NULL) return FALSE;

		Vec3D a,c,nrm;

		Vec3D::diff(a,*ptr_c1p,*ptr1);
		Vec3D::diff(c,*ptr2,*ptr1);

		Vec3D::VecProduct(nrm,a,c);
		nrm.normalize();

		std_crd_sys.resize(4);

		std_crd_sys[2] = nrm;

        a.normalize();
		Vec3D cnt = a;
		cnt.Scale(4.5033);

		double ca = cos(DEG_TO_RAD*132.193);
		double sa = sin(DEG_TO_RAD*132.193);
		
		cnt.Rotate(nrm,ca,sa);

		Vec3D::sum(std_crd_sys[3],(*ptr1),cnt);

		Vec3D v2 = a;

		ca = cos(-DEG_TO_RAD*54.512);
		sa = sin(-DEG_TO_RAD*54.512);

		v2.Rotate(nrm,ca,sa);

		std_crd_sys[1] = v2;

		Vec3D v1; 
		Vec3D::VecProduct(v1,v2,nrm);

		std_crd_sys[0] = v1;

		return TRUE;
	}

	return FALSE;
}



/* Note: curr == prev->GetNextResInChain()! */
double HaResidue::CalcPhiAngle(const HaResidue* prev, const HaResidue* curr) 
{
    const HaAtom* prevc;
    const HaAtom* currca;
    const HaAtom* currc;
    const HaAtom* currn;

    if( !(prevc  = prev->GetAtomByName("C" )) ) return( 360.0 );
    if( !(currca = curr->GetAtomByName("CA")) ) return( 360.0 );
    if( !(currc  = curr->GetAtomByName("C" )) ) return( 360.0 );
    if( !(currn  = curr->GetAtomByName("N" )) ) return( 360.0 );

    return( RAD_TO_DEG*HaAtom::CalcDihedral(prevc,currn,currca,currc) );
}


/* Note: next == curr->GetNextResInChain()! */
double  HaResidue::CalcPsiAngle(const HaResidue* curr, const HaResidue* next ) 
{
    const HaAtom* nextn;
    const HaAtom* currca;
    const HaAtom* currc;
    const HaAtom* currn;

    if( !(nextn  = next->GetAtomByName("N" )) ) return( 360.0 );
    if( !(currca = curr->GetAtomByName("CA")) ) return( 360.0 );
    if( !(currc  = curr->GetAtomByName("C" )) ) return( 360.0 );
    if( !(currn  = curr->GetAtomByName("N" )) ) return( 360.0 );

    return( RAD_TO_DEG*HaAtom::CalcDihedral(currn,currca,currc,nextn) );
}

int HaResidue::InterpolResParams(const char* res_fname_1, const char* res_fname_2, double weight_1)
{
//!  Set parameters of the residue intermediate between template residues  res_fname_1 and res_fname_2
//!  Atom parameters(charges) the same in two templates are not modified (provision if they are altered from the standard
//!  for example if there are more than one alternative proteonation states for a residue)
//!
    HaResidue* templ_res_1;
    HaResidue* templ_res_2;

	HaResDB* p_res_db = HaResDB::GetDefaultResDB();	
	templ_res_1 = p_res_db->GetTemplateForResidue(res_fname_1);

	if( templ_res_1 == NULL)
	{
        PrintLog("Error in HaResidue::InterpolResParams() \n");
		PrintLog("Can't find a residue template %s in the DB \n",res_fname_1);
		return FALSE;
	}

	bool do_weight = true;
	if( fabs(weight_1 - 1.0) < 1.0E-6 ) do_weight = false;

	double weight_2 = 1.0 - weight_1;

    templ_res_2 = p_res_db->GetTemplateForResidue(res_fname_2);

	if(do_weight && templ_res_2 == NULL)
	{
        PrintLog("Error in HaResidue::InterpolResParams() \n");
		PrintLog("Can't find a template for the residue %s \n",res_fname_2);
		return FALSE;
	}

	AtomIteratorResidue aitr(this);
	AtomIteratorResidue aitr1(templ_res_1);
	AtomIteratorResidue aitr2(templ_res_2);

    HaAtom* aptr;
	for( aptr= aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
	{
		HaAtom* aptr1 = NULL; 
		HaAtom* aptr2 = NULL;
		std::string at_name = aptr->GetName();
		if( templ_res_1 ) aptr1 = templ_res_1->GetAtomByName(at_name.c_str());
        if( templ_res_2 ) aptr2 = templ_res_2->GetAtomByName(at_name.c_str());

		double ch = 0.0;
		double ch1 = aptr1->GetCharge();
        double ch2 = aptr2->GetCharge();

		if( fabs(ch1 - ch2) > 1.0E-6)
		{
		    if( aptr1 ) ch += aptr1->GetCharge() * weight_1 ;
		    if( aptr2 ) ch += aptr2->GetCharge() * weight_2 ;
		    aptr->SetCharge(ch);
		}
	}
        
    return TRUE;
}

HaResidue* HaResidue::GetTemplate()
{
	HaResidue* res_templ;
	std::string res_full_name = GetFullName();
	HaResDB* p_res_db = HaResDB::GetDefaultResDB();	
	res_templ = p_res_db->GetTemplateForResidue(res_full_name.c_str());
	return res_templ;
}

int HaResidue::CheckStruct()
{
	if( MMForceField::ff_type_default == ForceFieldType::AMOEBA)
	{
		return CheckStructMortLib( MMForceField::ff_type_default );
	}

	try
	{
		std::string res_ref = GetRef();
		std::string res_fname = GetFullName();

		HaResidue* res_templ = GetTemplate();
		if(res_templ == NULL) throw std::runtime_error( " Residue " + res_ref + " does not have a DB template " );

		HaAtom* aptr;
		HaAtom* aptr_templ;

		AtomIteratorAtomGroup aitr(this);
		int n_extra = 0;
		int n_miss = 0;
		int i_res_valid = TRUE;
		for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
			std::string atname = aptr->GetName();
			aptr_templ = res_templ->GetAtomByName(atname.c_str());
			if( aptr_templ == NULL)
			{
				if(n_extra == 0)
				{
					PrintLog(" Residue %s (%s) have atoms not found in template: \n", 
						res_ref.c_str(), res_fname.c_str());
				}	 
				n_extra++;
				i_res_valid = FALSE;
				PrintLog( " %s ", atname.c_str());
			}
		}
		if(n_extra > 0) PrintLog( " \n ");

		AtomIteratorAtomGroup aitr_t(res_templ);
		for(aptr_templ = aitr_t.GetFirstAtom(); aptr_templ; aptr_templ = aitr_t.GetNextAtom())
		{
			if( aptr_templ->IsProxy() ) continue;
			std::string atname = aptr_templ->GetName();
			aptr = this->GetAtomByName(atname.c_str());
			if( aptr == NULL)
			{
				if(n_miss == 0)
				{
					PrintLog(" Residue %s (%s) does not have atoms found in template: \n", 
						res_ref.c_str(), res_fname.c_str());
				}
				n_miss++;
				i_res_valid = FALSE;
				PrintLog( " %s ", atname.c_str());
			}
		}

		if(n_miss > 0) throw std::runtime_error(" "); 

		if( i_res_valid ) 
		{
			aptr = aitr.GetFirstAtom();
			aptr_templ = aitr_t.GetFirstAtom();
			for( ; aptr; aptr = aitr.GetNextAtom(),aptr_templ = aitr_t.GetNextAtom() )
			{
				std::string at_name = aptr->GetName();
				std::string at_name_templ = aptr_templ->GetName();
				if( at_name != at_name_templ ) throw std::runtime_error(" Residue " + res_ref +  " has atom order different from that in the template " + res_fname);
			}
		}
	}
	catch( std::exception& ex )
	{
		PrintLog("%s \n",ex.what());
		return FALSE;
	}
	return TRUE;
}


int HaResidue::SetStdCharges()
{
	HaResidue* res_templ = GetTemplate();
	if(res_templ == NULL) return FALSE;

	AtomIteratorResidue aitr(this);
	HaAtom* aptr;
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		std::string at_name = aptr->GetName();
		HaAtom* aptr_templ = res_templ->GetAtomByName(at_name.c_str());
		if(aptr_templ == NULL) continue;
		double ch = aptr_templ->GetCharge();
		aptr->SetCharge(ch);
	}
	return TRUE;
}


AtomIteratorGen AtomContainer::__iter__()
{
	return AtomIteratorGen(this);
}

AtomIteratorGen AtomContainer::GetAtomIterator()
{
	return AtomIteratorGen(this);
}

AtomIteratorGen::AtomIteratorGen( AtomContainer* pat_cont_new)
{
	pat_cont = pat_cont_new;
	pat_itr  = pat_cont->GetAtomIteratorPtr();
	first_called = 0;

}

AtomIteratorGen::AtomIteratorGen( const AtomIteratorGen& ref)
{
	this->pat_cont = ref.pat_cont;
	pat_itr  = this->pat_cont->GetAtomIteratorPtr();
	first_called = 0;
}


AtomIteratorGen::~AtomIteratorGen()
{
	delete pat_itr;
}

AtomIteratorGen AtomIteratorGen::__iter__() const
{
	return (*this);
}

HaAtom* AtomIteratorGen::next()
{
	HaAtom* aptr;
	if(first_called) 
	{
		aptr = this->GetNextAtom();
	}
	else
	{
		aptr = GetFirstAtom();
		first_called = 1;
	}
	if( aptr == NULL)
	{
		throw std::out_of_range("Stop Atom Iterations");	
	}
	return aptr;
}

HaAtom* AtomIteratorGen::__next__()
{
	return next();
}




AtomIteratorAtomGroup::AtomIteratorAtomGroup(AtomGroup* new_p_at_group)
{
	p_at_group = new_p_at_group;
	if(p_at_group == NULL) { return; }
	
    aitrm = p_at_group->begin(); 
}

AtomIteratorAtomGroup::~AtomIteratorAtomGroup()
{
	
}

HaAtom* AtomIteratorAtomGroup::GetFirstAtom()
{
	if(p_at_group == NULL) { return NULL; }
	aitrm = p_at_group->begin();   
	if( aitrm == p_at_group->end() )
	{
		return NULL;
	}   
	return (*aitrm);
}

HaAtom* AtomIteratorAtomGroup::GetNextAtom()
{
	if(p_at_group == NULL) { return NULL; }
	aitrm++;
	if( aitrm == p_at_group->end() )
	{
		return NULL;
	}
	return (*aitrm);
}

AtomIteratorAtomGroup_const::AtomIteratorAtomGroup_const(const AtomGroup* new_p_at_group)
{
	p_at_group = new_p_at_group;
	if(p_at_group == NULL) { return; }
	
    aitrm = p_at_group->begin(); 
}

AtomIteratorAtomGroup_const::~AtomIteratorAtomGroup_const()
{

}

const HaAtom* AtomIteratorAtomGroup_const::GetFirstAtom() 
{
	if(p_at_group == NULL) { return NULL; }
	aitrm = p_at_group->begin();   
	if( aitrm == p_at_group->end() )
	{
		return NULL;
	}   
	return (*aitrm);
}

const HaAtom* AtomIteratorAtomGroup_const::GetNextAtom() 
{
	if(p_at_group == NULL) { return NULL; }
	aitrm++;
	if( aitrm == p_at_group->end() )
	{
		return NULL;
	}
 
	return (*aitrm);
}

AtomIteratorChain::AtomIteratorChain()
{
	chain = NULL;
}

AtomIteratorChain::AtomIteratorChain(HaChain* new_chain)
{
	chain = new_chain;
	if(chain == NULL) return;

	GetFirstAtom();
}

AtomIteratorChain::~AtomIteratorChain()
{

}

int AtomIteratorChain::SetForChain(HaChain* new_chain)
{
	chain = new_chain;
	if(chain == NULL) return FALSE;

	GetFirstAtom();
	return TRUE;
}

HaAtom* AtomIteratorChain::GetFirstAtom()
{
	for( ritrm = chain->res_arr.begin(); ritrm != chain->res_arr.end(); ritrm++)
	{
		if( (*ritrm)->GetNAtoms() > 0) break;
	}

	if( ritrm == chain->res_arr.end() ) return NULL;

	aitr = (*ritrm)->begin();
	return (*aitr);
}

HaAtom* AtomIteratorChain::GetNextAtom()
{
	aitr++;
	if( aitr != (*ritrm)->end() ) return (*aitr);

	ritrm++;
	for(; ritrm != chain->res_arr.end(); ritrm++)
	{
		if( (*ritrm)->GetNAtoms() > 0) break;
	}
	if( ritrm == chain->res_arr.end() ) return NULL;

	aitr = (*ritrm)->begin();
	return (*aitr);
}

AtomIteratorChain_const::AtomIteratorChain_const()
{
	chain = NULL;
}

AtomIteratorChain_const::AtomIteratorChain_const(const HaChain* new_chain)
{
	chain = new_chain;
	if(chain == NULL) return;

	GetFirstAtom();
}

AtomIteratorChain_const::~AtomIteratorChain_const()
{

}

int AtomIteratorChain_const::SetForChain(const HaChain* new_chain)
{
	chain = new_chain;
	if(chain == NULL) return FALSE;

	GetFirstAtom();
	return TRUE;
}

const HaAtom* AtomIteratorChain_const::GetFirstAtom()
{
	for( ritrm = chain->res_arr.begin(); ritrm != chain->res_arr.end(); ritrm++)
	{
		if( (*ritrm)->GetNAtoms() > 0) break;
	}

	if( ritrm == chain->res_arr.end() ) return NULL;

	aitr = (*ritrm)->begin();
	return (*aitr);
}

const HaAtom* AtomIteratorChain_const::GetNextAtom()
{
	aitr++;
	if( aitr != (*ritrm)->end() ) return (*aitr);

	ritrm++;
	for(; ritrm != chain->res_arr.end(); ritrm++)
	{
		if( (*ritrm)->GetNAtoms() > 0) break;
	}
	if( ritrm == chain->res_arr.end() ) return NULL;

	aitr = (*ritrm)->begin();
	return (*aitr);
}

ResidueIteratorChain::ResidueIteratorChain(HaChain* new_chain)
{
	chain = new_chain;
	if(chain == NULL) { return; }
	
	GetFirstRes();
}

ResidueIteratorChain::ResidueIteratorChain(const ResidueIteratorChain& ritr_ref)
{
	chain = ritr_ref.chain;
	res_itr = ritr_ref.res_itr;
}


ResidueIteratorChain::~ResidueIteratorChain()
{

}

HaResidue* ResidueIteratorChain::GetFirstRes()
{
   if(chain == nullptr) return nullptr; 
   res_itr = chain->res_arr.begin();
   if(res_itr == chain->res_arr.end()) return nullptr;
   return (*res_itr);
}

HaResidue* ResidueIteratorChain::GetNextRes()
{
  if(chain == nullptr) { return nullptr; }

  res_itr++;

  if( res_itr != chain->res_arr.end() )
  {
	   return (*res_itr);
  }
  return nullptr;
}

HaChain::HaChain()
{
	SetDefaultParam();
}

HaChain::HaChain(HaMolecule* new_phost_mol, const char new_ident)
{
	SetDefaultParam();
	phost_mol=new_phost_mol;
	ident = new_ident;
}


HaChain::~HaChain()
{
	int i, nr = res_arr.size();
	for(i = 0; i < nr; i++)
	{
		HaResidue* pres = res_arr[i];
		delete pres;
	}
	res_arr.clear();
	res_map.clear();
}

int HaChain::GetNAtoms() const
{
	int na = 0;
	int i, nr = res_arr.size();
	for (i = 0; i < nr; i++)
	{
		HaResidue* pres = res_arr[i];
		na += pres->GetNAtoms();
	}
	return na;
}

int HaChain::IsMember(const HaAtom* aptr) const
{
	if( aptr == NULL) return FALSE;
	if( aptr->GetHostChain() == this) return TRUE;

	return FALSE;
}

int HaChain::GetNumPt() const
{
	return GetNAtoms();
}


void HaChain::SetDefaultParam()
{
	ident = ' ';
	phost_mol=NULL;
}

bool HaChain::SetParamFrom(const HaChain& chain_ref)
{
	ident = chain_ref.ident;
	return true;
}

PointIterator* HaChain::GetPointIteratorPtr()
{
	return GetAtomIteratorPtr();
}

PointIterator_const* HaChain::GetPointIteratorPtr() const
{
	return GetAtomIteratorPtr();
}

HaResidue* HaChain::AddResidue(int res_ser_no)
{
	if( res_map.count(res_ser_no) > 0 )
	{  
	   PrintLogCount(1, " Warning in HaChain::AddResidue() \n  Residue Number %d is not unique \n", res_ser_no );  // Limit number of prints 
	}
	HaResidue* pres = new HaResidue(this);
	res_arr.push_back(pres);

	pres->serno = res_ser_no;
	std::pair<int,HaResidue*> ir_pair(res_ser_no,pres);
	res_map.insert(ir_pair);
	return(pres);
}

static const int MIN_MAIN_RES_SERNO=0;
static const int MIN_TERM_RES_SERNO=5000;

int HaChain::GetUniqResSerNo(int term_res_flag) const
{
	std::multimap<int, HaResidue*, less<int> >::const_iterator ritr;
	int new_serno;
	
	if(term_res_flag)
		new_serno=MIN_TERM_RES_SERNO;
	else
		new_serno=MIN_MAIN_RES_SERNO;

	for(ritr = res_map.begin(); ritr != res_map.end(); ritr++)
	{
		if((*ritr).second->GetSerNo() > new_serno) new_serno= (*ritr).second->GetSerNo();
	}
	new_serno++;
	return(new_serno);
}

bool HaChain::SetUniqueResNo()
{
	std::multimap<int, HaResidue*, less<int> >::iterator ritr;

	int last_serno = -999999;
	int found_same_serno = FALSE;
	for(ritr = res_map.begin(); ritr != res_map.end(); ritr++ )
	{
		int cur_serno= (*ritr).second->GetSerNo();
		if( cur_serno == last_serno )
		{
			found_same_serno = TRUE;
			PrintLog(" Found Residue with the same serial numbers in the chain \n");
			PrintLog(" Will Renumber \n");
			break;
		}
		last_serno = cur_serno;
	}
	if( !found_same_serno ) return true;

	size_t const nr = res_arr.size();
	res_map.clear();

	last_serno = -999999;
	int ir;
	for( ir = 0; ir < nr; ir++ )
	{
		HaResidue* pres = res_arr[ir];
		int cur_serno = pres->GetSerNo();
		if( cur_serno <= last_serno ) 
		{
			cur_serno = last_serno + 1;
		}
		pres->serno = cur_serno;
		std::pair<int,HaResidue*> ir_pair(cur_serno,pres);
		res_map.insert(ir_pair);
		last_serno = cur_serno;
	}
	return true;
}

HaResidue* HaChain::GetFirstRes()
{
	if( res_map.empty() ) return NULL;
	HaResidue* pres = (*res_arr.begin());
	return pres;
}

HaResidue* HaChain::GetResBySerNo(const int res_ser_no)
{
	std::multimap<int,HaResidue*, less<int> >::iterator ritr = res_map.find(res_ser_no);
	if( ritr != res_map.end() ) return (*ritr).second;
	return NULL;
}

HaAtom* HaResidue::GetAtomByName(const std::string& at_ref )
{
	vector<HaAtom*>::iterator aitr;
    for( aitr=begin(); aitr != end(); aitr++ )
	{
        if( !stricmp_trunc( (*aitr)->GetName(), at_ref.c_str() )) return( *aitr );
	}

	int ipos = at_ref.find('.');
	if( ipos != std::string::npos)
	{
		std::string res_prefix = at_ref.substr(0,ipos);
		std::string atn = at_ref.substr(ipos+1);
		boost::to_upper(res_prefix);
		boost::to_upper(atn);
		HaResidue* pres_n = NULL;
		if( res_prefix == "PREV" )
		{
			pres_n = this->GetPrevResInChain();
		}
		if( res_prefix == "NEXT" )
		{
			pres_n = this->GetNextResInChain();
		}
		if( pres_n != NULL )
		{
			HaAtom* aptr = pres_n->GetAtomByName(atn.c_str());
			return aptr;
		}
	}
    return( NULL );	
}

const HaAtom* HaResidue::GetAtomByName(const std::string& at_ref) const
{
    vector<HaAtom*>::const_iterator aitr;
    for( aitr=begin(); aitr != end(); aitr++ )
	{
        if( !stricmp_trunc( (*aitr)->GetName(), at_ref.c_str() )) return( *aitr );
	}
	
	int ipos = at_ref.find('.');
	if( ipos != std::string::npos)
	{
		std::string res_prefix = at_ref.substr(0,ipos);
		std::string atn = at_ref.substr(ipos+1);
		boost::to_upper(res_prefix);
		boost::to_upper(atn);
		HaResidue* pres_n = NULL;
		if( res_prefix == "PREV" )
		{
			pres_n = this->GetPrevResInChain();
		}
		if( res_prefix == "NEXT" )
		{
			pres_n = this->GetNextResInChain();
		}
		if( pres_n != NULL )
		{
			HaAtom* aptr = pres_n->GetAtomByName(atn.c_str());
			return aptr;
		}
	}
    return( NULL );	
}

HaAtom* HaResidue::AddNewAtom()
{
	HaAtom* aptr = new HaAtom;
	aptr->SetHostRes(this);

	this->push_back(aptr);
 	
    if( this->IsSolvent() ) aptr->flag |= HeteroFlag;

	return( aptr );
}

static const int ALL_NEIB_RESOLVED   = 0x0001;
static const int HAS_NO_RESOLVED_NEIB = 0x0002;
static const int HAS_NO_NEIB = 0x0004;
static const int IS_RESOLVED = 0x0008; 
static const int IS_EXTERNAL_AT = 0x0010;


class AtomNode 
//! Axxiliary class for attaching missing atoms in a residue using a residue template
{
public:
	AtomNode() {}
	AtomNode(const AtomNode& ref_node)
	{
		at = ref_node.at;
		ar = ref_node.ar;
		ntype = ref_node.ntype;
		nresolv_neib = ref_node.nresolv_neib;

		bonded_nodes = ref_node.bonded_nodes;
	}
	virtual ~AtomNode() {}

	HaAtom* at; //!< Atom of the residue template 
	HaAtom* ar; //!< Corresponding Atom of the current residue
	            // if ar == NULL - the atom is missing	
	int ntype;  //!< Indicator of resolved state of the atom and his neighbors
	int nresolv_neib; //!< the numbor of resolved neighbor atoms
	vector<int> bonded_nodes; //!< the node indexes the template atoms(nodes) connected 
	                          //!< to the template atom in the template residue
};


static std::vector<AtomNode> nodes;

int HaResidue::AddMissingAtoms( ADD_ATOM_TYPE atom_type )
{
	char buf[256];char buf2[256];

	std::string fres_name = this->GetFullName();
	std::string res_ref   = this->GetRef();

	try
	{	
		HaResDB* p_res_db = HaResDB::GetDefaultResDB();	
		
		HaMolecule* res_templ = p_res_db->GetMolTemplForRes( fres_name.c_str() );
		HaResidue* prtempl = p_res_db->GetTemplateForResidue( fres_name.c_str() );
		if( prtempl == NULL || res_templ == NULL ) throw runtime_error(" Can't find residue template " + fres_name + " in the database" );

		if( this->IsWater())
		{
			return AddWaterHydrogens();
		}

		HaAtom* aptr; HaAtom* atempl;

		nodes.clear();
		AtomIntMap at_nodes_map; // map connecting atom pointer to their positions in nodes vector; 

		std::set<std::string, less<std::string> > at_names;

		AtomIteratorAtomGroup aitr_this(this);

		// Fill Atom<->Atom template nodes array for atoms with templates that can be resolved with atom names 

		for(aptr= aitr_this.GetFirstAtom(); aptr; aptr= aitr_this.GetNextAtom())
		{
			std::string at_name = aptr->GetName();
			if( at_names.find(aptr->GetName()) != at_names.end() ) throw runtime_error(" Residue has multiple atoms with the name " + at_name ); 

			atempl= prtempl->GetAtomByName(aptr->GetName());
			if(atempl == NULL) throw runtime_error(" Residue Template doesn't have atom " + at_name + " found in current residue ");

			AtomNode node;
			node.at = atempl;
			node.ar = aptr;
			node.ntype = IS_RESOLVED;
			node.nresolv_neib = 0; 
			nodes.push_back(node);
			at_nodes_map[atempl] = nodes.size() - 1;
		}

		int residue_complete = TRUE;


		// Fill nodes array for Non-proxy atoms that can not be resolved by atom names 

		AtomIteratorMolecule aitr_templ(res_templ);
		for(atempl = aitr_templ.GetFirstAtom(); atempl; atempl = aitr_templ.GetNextAtom())
		{
			if( atempl->IsProxy() ) continue; 
			aptr= this->GetAtomByName(atempl->GetName());
			if(aptr == NULL)
			{
				residue_complete = FALSE;
				AtomNode node;
				node.at = atempl;
				node.ar = NULL;
				node.ntype = 0;
				node.nresolv_neib = 0;
				nodes.push_back(node);
				at_nodes_map[atempl] = nodes.size() - 1;
			}
		}
		if(residue_complete) return TRUE;
			
		int i;
		int nsize = nodes.size();

 // Match to proxy atoms of the residue template to atoms of residues connected to the given residue: 

		for( i = 0; i < nsize; i++)
		{
			atempl = nodes[i].at;
			AtomGroup bonded_templ_atoms;
			atempl->GetBondedAtoms(bonded_templ_atoms);
			AtomIteratorAtomGroup aitr_bonded_templ_atoms(&bonded_templ_atoms);

			if(nodes[i].ar == NULL) continue; // template atom is not matched yet...
			aptr = nodes[i].ar;
			AtomGroup bonded_atoms;

			int found_boundary_atom = FALSE;
			aptr->GetBondedAtoms(bonded_atoms);

			AtomIteratorAtomGroup aitr_bonded_atoms(&bonded_atoms);

			HaAtom* bond_at_templ;
			for( bond_at_templ = aitr_bonded_templ_atoms.GetFirstAtom(); bond_at_templ; bond_at_templ = aitr_bonded_templ_atoms.GetNextAtom())
			{
				if(bond_at_templ->GetHostRes() != prtempl || bond_at_templ->IsProxy() ) // proxy atom of the template (does not belong to the residue template) 
				{
					HaAtom* bond_at;
					for(bond_at = aitr_bonded_atoms.GetFirstAtom(); bond_at; bond_at= aitr_bonded_atoms.GetNextAtom())
					{
						if( bond_at->GetHostRes() != this )
						{
							if( bond_at_templ->IsProxy() && !bond_at->IsMatch(bond_at_templ) ) continue; // Check if proxy atom of the residue template match by name the atom of the connected residue 
							AtomNode node;
							node.at = bond_at_templ;
							node.ar = bond_at;
							node.ntype = IS_RESOLVED | IS_EXTERNAL_AT;
							node.nresolv_neib = 1;
							nodes.push_back(node);
							at_nodes_map[bond_at_templ] = nodes.size() - 1;
							found_boundary_atom = TRUE;
							break;
						}
					}
				}
				if(found_boundary_atom) break;
			}
		}

		// Fill Bonding information of nodes    

		for(i = 0; i < nodes.size(); i++ )
		{
			atempl= nodes[i].at;
			AtomGroup bonded_templ_atoms;
			int ibt= atempl->GetBondedAtoms(bonded_templ_atoms);
			// PrintLog("Atom template %s has %d bonded atoms \n", atempl->GetRef().c_str(), bonded_templ_atoms.size());
			AtomIteratorAtomGroup aitr_bonded_templ_atoms(&bonded_templ_atoms);

			int found_boundary_atom = FALSE;

			HaAtom* bond_at_templ;
			for( bond_at_templ = aitr_bonded_templ_atoms.GetFirstAtom(); bond_at_templ; bond_at_templ = aitr_bonded_templ_atoms.GetNextAtom())
			{
				if(at_nodes_map.find( bond_at_templ ) != at_nodes_map.end())
				{
					int i_bond_node = at_nodes_map[bond_at_templ];
					nodes[i].bonded_nodes.push_back( i_bond_node );
					if( nodes[i_bond_node].ntype | IS_RESOLVED )  nodes[i].nresolv_neib++;
				}
			}
			// PrintLog("Atom template %s has %d bonded_nodes \n", atempl->GetRef().c_str(), nodes[i].bonded_nodes.size());
		}

		if( !AddMissingAtoms_2(prtempl, atom_type) ) throw std::runtime_error(" Error in AddMissingAtoms_2() ");
	}
	catch( std::exception& ex )
	{	
		PrintLog(" Error in HaResidue::AddMissingAtom() for residue %s \n",res_ref.c_str() );
		PrintLog("%s\n",ex.what());
		return FALSE;
	}
	return TRUE;
}

int HaResidue::AddWaterHydrogens()
{
	MolSet* pmset = GetHostMolSet();
	HaAtom* aptr_o = this->GetAtomByName("O");
	if(aptr_o == NULL)
	{
		PrintLog(" Error adding Hydrogens to the water molecule: No O atom \n");
		return FALSE;
	}
	HaAtom* aptr_h1 = this->GetAtomByName("H1");
	HaAtom* aptr_h2 = this->GetAtomByName("H2");
	
	if(aptr_h1 == NULL && aptr_h2 == NULL)
	{
		aptr_h1 = this->AddNewAtom();
		if(aptr_h1 == NULL)
		{
			PrintLog(" Error adding hydrogen to atom water molecule \n");
			return FALSE;
		}
		aptr_h1->SetName("H1");
		aptr_h1->SetElemNo(1);

		aptr_h1->SetX( aptr_o->GetX() + 0.957);
		aptr_h1->SetY( aptr_o->GetY() );
		aptr_h1->SetZ( aptr_o->GetZ() );

		HaBond* bptr= pmset->AddBond(aptr_h1, aptr_o );
		if(bptr) bptr->DrawWire();
		
		aptr_h2 = this->AddNewAtom();
		if( aptr_h2 == NULL)
		{
			PrintLog(" Error adding hydrogen atom to water molecule \n");
			return FALSE;
		}
		aptr_h2->SetName("H2");
		aptr_h2->SetElemNo(1);

		aptr_h2->SetX( aptr_o->GetX() - 0.24);
		aptr_h2->SetY( aptr_o->GetY() + 0.927);
		aptr_h2->SetZ( aptr_o->GetZ() );

		bptr = pmset->AddBond(aptr_h2, aptr_o );
		if(bptr) bptr->DrawWire();
	}
	return TRUE;
}
	
int HaResidue::AddMissingAtoms_2( HaResidue* prtempl, ADD_ATOM_TYPE atom_type)
// Adding atoms to the residue using the filled nodes array mapping atoms of the residue to the residue template
{ 
	MolSet* pmset = GetHostMolSet();
	HaAtom* atempl;

    int i;
	for(;;)
	{
		int residue_complete= TRUE;

		for(i = 0 ; i < nodes.size(); i++)
		{
			HaAtom* aptr2 = NULL;   HaAtom* aptr3  = NULL; HaAtom* aptr4   = NULL;
			HaAtom* atempl2 = NULL; HaAtom* atempl3= NULL; HaAtom* atempl4 = NULL;

			aptr2   = nodes[i].ar; // to this atom missing atoms will be attached
			atempl2 = nodes[i].at; // atom corresponding to aptr2 in the residue template

			if( aptr2 == NULL) continue; // if atom is not resolved go to the next node
				
			if(nodes[i].ntype & ALL_NEIB_RESOLVED ) continue;
				 
			if(nodes[i].nresolv_neib == 0 ) continue;
				
			if( (atom_type & ADD_POLAR_HYDROGENS) && !aptr2->IsHBDonor()) continue;
				
			// find atoms aptr3 and aptr4 and corresponding atoms atempl3, atempl4 
			// in the template residue which will be used to 
			// attach new atom aptr in the same relative position as its template atempl

			int j, iat; int iat3 = -1;
			for( j =0; j < nodes[i].bonded_nodes.size(); j++)
			{
				iat = nodes[i].bonded_nodes[j];
				if( nodes[iat].ar == NULL) continue;

				if(aptr3 == NULL)
				{
					aptr3   = nodes[iat].ar;
					atempl3 = nodes[iat].at;
					iat3 = iat;
					continue;
				}
				else 
				{
					if (nodes[iat].ar != aptr3)
					{
						aptr4 = nodes[iat].ar;
						atempl4 = nodes[iat].at;
						break;
					}
				}
			}
			if(aptr3 == NULL) continue;

			if(aptr4 == NULL)
			{
				for(j = 0; j < nodes[iat3].bonded_nodes.size(); j++)
				{
					iat = nodes[iat3].bonded_nodes[j];
					if(nodes[iat].ar == NULL)  continue;
					if(nodes[iat].ar == aptr2) continue;
					aptr4   = nodes[iat].ar;
					atempl4 = nodes[iat].at;
					// PrintLog("Added aptr4: %s as a neighbor of aptr3: %s \n", aptr4->GetRef().c_str(), aptr3->GetRef().c_str());
					break;
				}
			}
			
			for(j =0; j < nodes[i].bonded_nodes.size(); j++) 
			{
				iat = nodes[i].bonded_nodes[j];
				if(nodes[iat].ar == NULL) // if atom is missing in the structure 
				{
					atempl = nodes[iat].at; // the template of the missing atom
					if( (atom_type & (ADD_HYDROGENS | ADD_POLAR_HYDROGENS)) && 
						!atempl->IsHydrogen()) continue;

					if( (atom_type & ADD_HEAVY_ATOMS) && atempl->IsHydrogen())
						continue;

					std::string ns = "NULL";
					// PrintLog("Adding Atom for the template: %s \n", atempl->GetRef().c_str() );
					// PrintLog("atempl2: %s  atempl3: %s atempl4: %s \n", 
					//	atempl2 == NULL ? ns.c_str() : atempl2->GetRef().c_str(),   
					//	atempl3 == NULL ? ns.c_str() : atempl3->GetRef().c_str(),
					//	atempl4 == NULL ? ns.c_str() : atempl4->GetRef().c_str()
					//);
					// PrintLog("at2: %s  at3: %s at4: %s \n", 
					//	aptr2 == NULL ? ns.c_str() : aptr2->GetRef().c_str(),
					//	aptr3 == NULL ? ns.c_str() : aptr3->GetRef().c_str(),
					//	aptr4 == NULL ? ns.c_str() : aptr4->GetRef().c_str());
					// PrintLog("at2 crd : %8.3f %8.3f %8.3f \n at3 crd: %8.3f %8.3f %8.3f \n", 
					//	aptr2->GetX(), aptr2->GetY(), aptr2->GetZ(),
					//	aptr3->GetX(), aptr3->GetY(), aptr3->GetZ());
					HaAtom* new_atom= HaAtom::AddAtomFromTempl(aptr2,aptr3,aptr4, atempl, atempl2, atempl3, atempl4);
 					if(new_atom == NULL) continue;
					nodes[iat].ar = new_atom;

					int nbond = nodes[iat].bonded_nodes.size();
					int ib;
					for( ib = 0; ib < nbond; ib++ )
					{
						int iat2 = nodes[iat].bonded_nodes[ib];
						HaAtom* aptr_p = nodes[iat2].ar;
						if(aptr_p != NULL )
						{	
							HaBond* new_bond = pmset->AddBond(new_atom, aptr_p );
							new_bond->DrawWire();
						}
					}
					if(aptr4 == NULL ) 
					{
						aptr4 = new_atom;
						atempl4 = atempl;
					}
				}
			}
			nodes[i].ntype |= ALL_NEIB_RESOLVED;
		}
		if(residue_complete) return TRUE;
	}
	return TRUE;
}

HaAtom* HaAtom::AddAtomFromTempl( HaAtom* aptr2, HaAtom* aptr3, HaAtom* aptr4, 
		                     const HaAtom* aptr_templ, const HaAtom* aptr_templ_2, const HaAtom* aptr_templ_3, 
							 const HaAtom* aptr_templ_4) 
// add atom from a template
{
	HaAtom* new_atom = NULL;
	try 
	{
		if( (aptr2 == NULL) || (aptr3 == NULL) || (aptr4 == NULL) || (aptr_templ == NULL) || 
			(aptr_templ_2 == NULL) || (aptr_templ_3 == NULL) || (aptr_templ_4 == NULL) ) 
		{
			throw std::runtime_error("One of that atom pointers is NULL");
		}

		if (aptr2 == aptr3)
		{
			throw::runtime_error("Reference atoms aptr2 and patr3 are the same ");
		}

		if (aptr2->GetHostMolSet() != aptr3->GetHostMolSet() )
		{
			throw::runtime_error("Reference atoms aptr2 and aptr3 do not belong to the same Molecular Set");
		}

		double dist      = Vec3D::CalcDistance( aptr_templ, aptr_templ_2);
		double val_angle = Vec3D::CalcAngle   ( aptr_templ, aptr_templ_2, aptr_templ_3);

		double dih_angle;  

		if( aptr4 != NULL )
			dih_angle= Vec3D::CalcTorsion ( aptr_templ, aptr_templ_2, aptr_templ_3, aptr_templ_4);
		else
			dih_angle= 0.0;

		MolSet* pmset = aptr2->GetHostMolSet();

		if(aptr3->GetHostMolSet() != pmset || aptr4->GetHostMolSet() != pmset) 
		{
			throw::runtime_error("Atoms to attach do not belong to the same Molecular Set");
		}

		HaResidue* pres = aptr2->GetHostRes();

		new_atom = pres->AddNewAtom();
		if(new_atom == NULL)  throw::runtime_error("Can not add new atom");

		new_atom->SetParamFrom(*aptr_templ);

		Vec3D aptr4_2;
		if(aptr4 != NULL)
		{
			aptr4_2.SetX(aptr4->GetX()); 
			aptr4_2.SetY(aptr4->GetY()); 
			aptr4_2.SetZ(aptr4->GetZ()); 
		}
		else
		{
			aptr4_2.SetX(aptr3->GetX() + 1.0 ); 
			aptr4_2.SetY(aptr3->GetY() + 1.0 ); 
			aptr4_2.SetZ(aptr3->GetZ() + 1.0 ); 
		}

		Vec3D::SetAtomPos(new_atom, aptr2, aptr3, &aptr4_2, dist, val_angle, dih_angle);
	}
	catch(std::exception& ex )
	{
		PrintLog("Error in HaResidue::AddAtomFromTempl(): \n");
		if (aptr_templ != NULL) PrintLog("Adding atom for the template %s \n", aptr_templ->GetRef().c_str());
		PrintLog("%s\n",ex.what());
	}
	return new_atom;
}

PeriodicUnitInfo::PeriodicUnitInfo()
{
	spacegroup = "";
	orthogonal_flag = -1;
	octahedron_flag = false;

	ucell.resize(3);
	recip_ucell.resize(3);

	pbc_a = 0.0;
	pbc_b = 0.0;
	pbc_c = 0.0;

	alpha = 90.0*DEG_TO_RAD;
	beta  = 90.0*DEG_TO_RAD;
	gamma = 90.0*DEG_TO_RAD;

	cut_factor[0] = 0.0;
	cut_factor[1] = 0.0;
	cut_factor[2] = 0.0;

    reclng[0] = 0.0;
	reclng[1] = 0.0;
	reclng[2] = 0.0;
}

PeriodicUnitInfo::PeriodicUnitInfo(const PeriodicUnitInfo& ref)
{
	spacegroup = ref.spacegroup;
	orthogonal_flag = ref.orthogonal_flag;

	ucell = ref.ucell;
	recip_ucell = ref.recip_ucell;

	pbc_a = ref.pbc_a;
	pbc_b = ref.pbc_b;
	pbc_c = ref.pbc_c;

	alpha = ref.alpha;
	beta  = ref.beta;
	gamma = ref.gamma;

	cut_factor[0] = ref.cut_factor[0];
	cut_factor[1] = ref.cut_factor[1];
	cut_factor[2] = ref.cut_factor[2];

    reclng[0] = ref.reclng[0];
	reclng[1] = ref.reclng[1];
	reclng[2] = ref.reclng[2];
}

PeriodicUnitInfo::~PeriodicUnitInfo()
{

}

int PeriodicUnitInfo::IsSet() const 
{ 
	return !( orthogonal_flag == -1); 
}
	
int PeriodicUnitInfo::IsOrthogonal() const 
{ 
	return (orthogonal_flag == 1); 
}

bool PeriodicUnitInfo::IsOctahedron() const
{
	return octahedron_flag;
}

int PeriodicUnitInfo::IsValid() const 
{ 
	if( GetA() < 0.01 ) return FALSE;
	if( GetB() < 0.01 ) return FALSE;
	if( GetC() < 0.01 ) return FALSE;
	if( GetAlpha() < 0.01 ) return FALSE;
	if( GetBeta()  < 0.01 ) return FALSE;
	if( GetGamma() < 0.01 ) return FALSE;
	
	return TRUE;
}

void PeriodicUnitInfo::Set( bool to_set)
{
	if(to_set) 
	{
		if ( fabs( alpha - 90.0*DEG_TO_RAD ) < 0.01 && fabs( beta  - 90.0*DEG_TO_RAD ) < 0.01 && fabs( gamma - 90.0*DEG_TO_RAD ) < 0.01 )
		{
			orthogonal_flag = 1;
		}
		else
		{
			orthogonal_flag = 0;
		}
	}
	else
	{
		orthogonal_flag = -1;
	}
}

void PeriodicUnitInfo::SetOctahedron( bool to_set )
{
	if(to_set) 
	{
		octahedron_flag = true;
		orthogonal_flag = 0;
		alpha = 109.4712190 * DEG_TO_RAD;
		beta  = 109.4712190 * DEG_TO_RAD;
		gamma = 109.4712190 * DEG_TO_RAD;
	}
	else
	{
		octahedron_flag = false;
	}
}

int PeriodicUnitInfo::SetStdBox(AtomContainer* p_at_coll)
{
	if( p_at_coll == NULL) return FALSE;
	if( p_at_coll->GetNAtoms() == 0) return TRUE;

	double xmin =  1.0E+10;
	double ymin =  1.0E+10;
	double zmin =  1.0E+10;
	double xmax =  -1.0E+10;
	double ymax =  -1.0E+10;
	double zmax =  -1.0E+10;

	AtomIteratorGen aitr(p_at_coll);
	HaAtom* aptr;
	for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom())
	{
		double dval;
		double rad;
		rad = HaAtom::ElemVDWRadius(aptr->GetElemNo());
		dval = aptr->GetX() - rad;
		if( dval < xmin ) xmin = dval;
		dval = aptr->GetX() + rad;
		if( dval > xmax ) xmax = dval;
		dval = aptr->GetY() - rad;
		if( dval < ymin ) ymin = dval;
		dval = aptr->GetY() + rad;
		if( dval > ymax ) ymax = dval;
		dval = aptr->GetZ() - rad;
		if( dval < zmin ) zmin = dval;
		dval = aptr->GetZ() + rad;
		if( dval > zmax ) zmax = dval;
	}
	SetBox(xmax - xmin, ymax - ymin, zmax - zmin);

	return TRUE;
}


int PeriodicUnitInfo::SetBox(double a, double b, double c, double alpha_n, double beta_n, double gamma_n)
{
	if( a < 0.0001) a = 1.0;
	if( b < 0.0001) b = 1.0;
	if( c < 0.0001) c = 1.0;

	this->pbc_a = a;
	this->pbc_b = b;
	this->pbc_c = c;
	this->alpha = alpha_n;
	this->beta  = beta_n;
	this->gamma = gamma_n;


  if ( fabs( alpha - 90.0*DEG_TO_RAD ) < 0.01 && fabs( beta  - 90.0*DEG_TO_RAD ) < 0.01 && fabs( gamma - 90.0*DEG_TO_RAD ) < 0.01 )
  {
	  this->orthogonal_flag = 1;
  }
  else
  {
	  this->orthogonal_flag = 0;
  }

  int i;

  if (orthogonal_flag == 1)
  {
    for( i = 0; i < ucell.size(); i++)
		ucell[i].SetZeros();
    ucell[0][0] = a;
    ucell[1][1] = b;
    ucell[2][2] = c;
	for( i = 0; i < 3; i++)
		cut_factor[i] = 1.0;
  }
  else
  {
	  ucell[0][0] = a;
	  ucell[0][1] = 0.0;
	  ucell[0][2] = 0.0;
	  ucell[1][0] = b * cos(gamma);
	  ucell[1][1] = b * sin(gamma);
	  ucell[1][2] = 0.0; 
	  ucell[2][0] = c * cos(beta);
	  ucell[2][1] = ( b * c * cos(alpha) - ucell[2][0] * ucell[1][0] )/ucell[1][1] ;
	  ucell[2][2] = sqrt( c * c - ucell[2][0] * ucell[2][0] - ucell[2][1] * ucell[2][1] );

    //  Cut factors to correct for "spherical cutoff protrusion" for non-orhtogonal cells

      cut_factor[0] = 1.0 / (sin(beta) * sin(gamma))  ;
      cut_factor[1] = 1.0 / (sin(alpha) * sin(gamma)) ;
      cut_factor[2] = 1.0 / (sin(alpha) * sin(beta))  ;
  }

  if (orthogonal_flag == 1) 
  {
	for( i = 0; i < ucell.size(); i++)
		recip_ucell[i].SetZeros();
   
	recip_ucell[0][0] = 1.0/a;
    recip_ucell[1][1] = 1.0/b;
    recip_ucell[2][2] = 1.0/c;
    
    reclng[0] = a;
    reclng[1] = b;
    reclng[2] = c;
    ucell_vol = a * b * c;
  }
  else
  {
	  Vec3D u23,u31,u12;
	  Vec3D::VecProduct( u23, ucell[1], ucell[2]);
	  Vec3D::VecProduct( u31, ucell[2], ucell[0]);
	  Vec3D::VecProduct( u12, ucell[0], ucell[1]);

	  ucell_vol = Vec3D::DotProduct( ucell[0], u23);
    
	  for(i = 0; i < 3; i++)
	  {
		recip_ucell[0][i] = u23[i]/ucell_vol;
		recip_ucell[1][i] = u31[i]/ucell_vol;
		recip_ucell[2][i] = u12[i]/ucell_vol;
	  }

	  reclng[0] = 1.0/recip_ucell[0].length();  
	  reclng[1] = 1.0/recip_ucell[1].length();
      reclng[2] = 1.0/recip_ucell[2].length();
  }

  if (orthogonal_flag == 1 ) 
  {
	  ucell_sph = a;
	  if( b < ucell_sph) ucell_sph = b;
	  if( c < ucell_sph) ucell_sph = c;
	  ucell_sph = 0.5 * ucell_sph;
  }
  else
  {
     ucell_sph = a + b + c;
     for( i=0; i < 3; i++)
	 {
		 double dp = Vec3D::DotProduct(recip_ucell[i], ucell[i]);
         double dist = dp * reclng[i];
         if (dist < ucell_sph) ucell_sph = dist;
	 }
     ucell_sph = 0.5 * ucell_sph;
  }
  return TRUE;
}

int AtomContainer::GetStdPosition(HaMat_double& rot_std, Vec3D& trans_std)
//!
//! Determine standard rotation matrix and standard translation vector 
//! for a given orientation of the molecule
//!
{
	this->GetAverageCoord(trans_std[0], trans_std[1], trans_std[2]);
	if(!this->GetStdRotMat(rot_std))
	{
       return FALSE;
	}
	return TRUE;
}

int AtomContainer::GetStdPositionMomInertia(HaMat_double& rot_std, Vec3D& trans_std)
//!
//! Determine standard rotation matrix from the moment of inertia tensor 
//! and standard translation vector 
//! for a given orientation of the molecule
//!
{
	this->GetAverageCoord(trans_std[0], trans_std[1], trans_std[2]);
	if(!this->GetStdMomInertRotMat(rot_std))
	{
       return FALSE;
	}
	return TRUE;
}

int AtomContainer::SetPosition(const HaMat_double& rot_new, const Vec3D& trans_new)
//!
//! Rotate and translate the molecule according to standardized rotation
//! matrix rot_new and translation vector trans_new
//!
{
    Vec3D trans_old;
	HaMat_double rot_old(3,3), rt_diff(3,3);

	if(!this->GetStdPosition(rot_old,trans_old))
		return FALSE;
	 
	matmult_T2(rt_diff, rot_new, rot_old);

	this->RotateAtoms(rt_diff, trans_old);

	Vec3D diff = trans_new - trans_old;

	this->TranslateAtoms( diff );

	return TRUE;
}


int AtomContainer::SetPositionMomInertia(const HaMat_double& rot_new, const Vec3D& trans_new)
//!
//! Rotate and translate the molecule according to standardtized rotation
//! matrix rot_new (using principal axes of inertia) and translation vector trans_new
//!
{
    Vec3D trans_old;
	HaMat_double rot_old(3,3), rt_diff(3,3);
	Quaternion q;

	if(!this->GetStdPositionMomInertia(rot_old,trans_old))
		return FALSE;
//	 Quaternion::RotMatToQuaternion(rot_old, q);
	 //   PrintLog("rt_diff %2.3f %2.3f %2.3f\n",rot_old(1,1), rot_old(1,2), rot_old(1,3));
	//	PrintLog("rt_diff %2.3f %2.3f %2.3f\n",rot_old(2,1), rot_old(2,2), rot_old(2,3));
   	//	PrintLog("rt_diff %2.3f %2.3f %2.3f\n",rot_old(3,1), rot_old(3,2), rot_old(3,3));

//	 Quaternion::QuaternionToRotMat(q, rot_old);
 //  		PrintLog("rot_old %2.3f %2.3f %2.3f\n",rot_old(1,1), rot_old(1,2), rot_old(1,3));
	//	PrintLog("rot_old %2.3f %2.3f %2.3f\n",rot_old(2,1), rot_old(2,2), rot_old(2,3));
   	//	PrintLog("rot_old %2.3f %2.3f %2.3f\n",rot_old(3,1), rot_old(3,2), rot_old(3,3));
    //-------
	matmult_T2(rt_diff, rot_new, rot_old);
	//PrintLog("rt_diff %2.3f %2.3f %2.3f\n",rt_diff(1,1), rt_diff(1,2), rt_diff(1,3));
	//PrintLog("rt_diff %2.3f %2.3f %2.3f\n",rt_diff(2,1), rt_diff(2,2), rt_diff(2,3));
	//PrintLog("rt_diff %2.3f %2.3f %2.3f\n",rt_diff(3,1), rt_diff(3,2), rt_diff(3,3));
	this->RotateAtoms(rt_diff, trans_old);

	Vec3D diff = trans_new - trans_old;

	this->TranslateAtoms( diff );

	return TRUE;
}

void AtomContainer::GetPosEulerTrans( double& phi, double& cos_theta, double& psi, Vec3D& trans)
{
    HaMat_double rot(3,3); 
	GetStdPosition(rot,trans);
	Rot3D::RotMatToEuler(rot,phi,cos_theta,psi);
}

int AtomContainer::SetPosEulerTrans(double phi, double cos_theta, double psi, const Vec3D& trans)
//!
//! Rotate and translate the molecule according to standard rotation
//! matrix determined by euler angles and translation vectors transx, transy, transz
//!
{
    HaMat_double rot(3,3); 
	Rot3D::EulerToRotMat(phi, cos_theta, psi, rot );
    return SetPosition(rot,trans);
}

void AtomContainer::GetQuaternionTrans(Quaternion& q, Vec3D& trans)
{
    HaMat_double rot(3,3);
	
	GetStdPositionMomInertia(rot,trans);
	//PrintLog("Get %2.3f %2.3f %2.3f \n", rot(1,1), rot(1,2), rot(1,3) );
	//PrintLog("Get %2.3f %2.3f %2.3f \n", rot(2,1), rot(2,2), rot(2,3) );
	//PrintLog("Get %2.3f %2.3f %2.3f \n", rot(3,1), rot(3,2), rot(3,3) );
	Quaternion::RotMatToQuaternion(rot, q);
//PrintLog("GET");
//q.PrintOn();

}

int AtomContainer::SetQuaternionTrans(const Quaternion& q, const Vec3D& trans)
//!
//! Rotate and translate the molecule according to quaternion rotation
//! and translation vector transx, transy, transz
//!
{
    HaMat_double rot(3,3); 
	Quaternion::QuaternionToRotMat(q, rot);

//PrintLog("SET");
//Quaternion r = q;
//r.PrintOn();
 // 	PrintLog("set rot_mat %2.3f %2.3f %2.3f\n",rot(1,1), rot(1,2), rot(1,3));
	//PrintLog("set rot_mat %2.3f %2.3f %2.3f\n",rot(2,1), rot(2,2), rot(2,3));
 // 	PrintLog("set rot_mat %2.3f %2.3f %2.3f\n",rot(3,1), rot(3,2), rot(3,3));
    
	return SetPositionMomInertia(rot,trans);
}

bool AtomContainer::GetStdRotMat(HaMat_double& rot_mat)
{
   rot_mat.newsize(3,3);

   rot_mat(1,1) = 1.0; rot_mat(1,2) = 0.0; rot_mat(1,3) = 0.0;
   rot_mat(2,1) = 0.0; rot_mat(2,2) = 1.0; rot_mat(2,3) = 0.0;
   rot_mat(3,1) = 0.0; rot_mat(3,2) = 0.0; rot_mat(3,3) = 1.0;


   int na = GetNAtoms();
   if( na < 2) return false;

   HaAtom* aptr;

   double x1,y1,z1;
   double x2,y2,z2;
   Vec3D v1,v2,v3;
   
   AtomIteratorGen aitr(this);

   aptr = aitr.GetFirstAtom();

   x1= aptr->GetX();
   y1= aptr->GetY();
   z1= aptr->GetZ();

   aptr = aitr.GetNextAtom();

   x2= aptr->GetX();
   y2= aptr->GetY();
   z2= aptr->GetZ();

   double xc, yc, zc, vlen;
   GetAverageCoord(xc,yc,zc);

   v1.SetX( x1 - xc);
   v1.SetY( y1 - yc);
   v1.SetZ( z1 - zc);

   vlen = v1.length();

   if(vlen <= DBL_EPSILON) return false;

   v1.Scale( 1.0/vlen);

   v2.SetX( x2 - xc);
   v2.SetY( y2 - yc);
   v2.SetZ( z2 - zc);

   if(v2.length() <= DBL_EPSILON) return false;

   Vec3D::VecProduct(v3,v1,v2);

   vlen = v3.length();

   if( vlen <= DBL_EPSILON) return false;

   v3.Scale( 1.0/vlen);

   Vec3D::VecProduct(v2,v3,v1);

   rot_mat(1,1) = v1.GetX(); rot_mat(1,2) = v2.GetX(); rot_mat(1,3) = v3.GetX();
   rot_mat(2,1) = v1.GetY(); rot_mat(2,2) = v2.GetY(); rot_mat(2,3) = v3.GetY();
   rot_mat(3,1) = v1.GetZ(); rot_mat(3,2) = v2.GetZ(); rot_mat(3,3) = v3.GetZ();

   return true;
}

bool AtomContainer::GetStdMomInertRotMat(HaMat_double& rot_mat)
{
   HaMat_double cc(3,3);
   HaVec_double eig_val(3);

   rot_mat.newsize(3,3);

   rot_mat(1,1) = 1.0; rot_mat(1,2) = 0.0; rot_mat(1,3) = 0.0;
   rot_mat(2,1) = 0.0; rot_mat(2,2) = 1.0; rot_mat(2,3) = 0.0;
   rot_mat(3,1) = 0.0; rot_mat(3,2) = 0.0; rot_mat(3,3) = 1.0;

   HaAtom* aptr;

   double x1,y1,z1;
  
   double xc, yc, zc;
   GetAverageCoord(xc,yc,zc);

   double c11= 0;
   double c22= 0;
   double c33= 0;
   double c12= 0;
   double c13= 0;
   double c23= 0;

   AtomIteratorGen aitr(this);
   for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
   {
	   //			name = aptr->GetName();
	   //			if(name == "C" || name == "CA" || name == "N")
	   //			{
	   x1 = aptr->GetX() - xc;
	   y1=  aptr->GetY() - yc;
	   z1 = aptr->GetZ() - zc;
	   
	   c11+= (y1*y1+ z1*z1)*(aptr->GetMass());
	   c22+= (x1*x1+ z1*z1)*(aptr->GetMass());
	   c33+= (x1*x1+ y1*y1)*(aptr->GetMass());
	   c12 -= x1*y1*(aptr->GetMass());
	   c13 -= x1*z1*(aptr->GetMass());
	   c23 -= y1*z1*(aptr->GetMass());
	   //			}
   }
   
   cc.SetVal_idx0(0,0, c11);
   cc.SetVal_idx0(0,1, c12);
   cc.SetVal_idx0(0,2, c13);
   cc.SetVal_idx0(1,0, c12);
   cc.SetVal_idx0(1,1, c22);
   cc.SetVal_idx0(1,2, c23);
   cc.SetVal_idx0(2,0, c13);
   cc.SetVal_idx0(2,1, c23);
   cc.SetVal_idx0(2,2, c33);
   HaMat_double::mat_sdiag(cc, rot_mat, eig_val);
//   PrintLog("eig_val %2.3f %2.3f %2.3f\n", eig_val[0], eig_val[1], eig_val[2]);
   Vec3D vec1;
   Vec3D vec2;
   Vec3D vec3;
   
   aptr = aitr.GetFirstAtom();
   aptr = aitr.GetNextAtom();
   x1 = aptr->GetX() - xc;
   y1=  aptr->GetY() - yc;
   z1 = aptr->GetZ() - zc;
   
   vec1.SetX(x1);
   vec1.SetY(y1);
   vec1.SetZ(z1);
//   double vlen = vec1.length();
//   vec1.Scale( 1.0/vlen);
   vec2.SetX(vec1[0]);
   vec2.SetY(vec1[1]);
   vec2.SetZ(vec1[2]);
   
   double proj= rot_mat(1,1)*vec1[0] + rot_mat(2,1)*vec1[1] + rot_mat(3,1)*vec1[2];
   if(proj < 0 )
   {
	   rot_mat(1,1) = -rot_mat(1,1);
	   rot_mat(2,1) = -rot_mat(2,1);
	   rot_mat(3,1) = -rot_mat(3,1);
//	   PrintLog("proj %2.3f \n", proj);
   }
   vec1[0] = rot_mat(1,1); 
   vec1[1] = rot_mat(2,1); 
   vec1[2] = rot_mat(3,1);
   proj= rot_mat(1,2)*vec2[0] + rot_mat(2,2)*vec2[1] + rot_mat(3,2)*vec2[2];
   if(proj > 0 )
   {
	   rot_mat(1,2) = -rot_mat(1,2);
	   rot_mat(2,2) = -rot_mat(2,2);
	   rot_mat(3,2) = -rot_mat(3,2);
//	   PrintLog("proj2 %2.3f \n", proj);
   }
   vec2[0] = rot_mat(1,2); vec2[1] = rot_mat(2,2); vec2[2] = rot_mat(3,2);
   Vec3D::VecProduct(vec3,vec1,vec2);   
   rot_mat(1,3) = vec3[0];
   rot_mat(2,3) = vec3[1];
   rot_mat(3,3) = vec3[2];



		
	//	PrintLog("Rrot_mat %2.3f %2.3f %2.3f\n",rot_mat(2,1), rot_mat(2,2), rot_mat(2,3));
   	//	PrintLog("Rrot_mat %2.3f %2.3f %2.3f\n",rot_mat(3,1), rot_mat(3,2), rot_mat(3,3));
//		PrintLog("cc %2.3f %2.3f %2.3f\n",cc(1,1), cc(2,1), cc(3,1));

   return true;
}

int AtomContainer::RotateAtoms( const HaMat_double& rot_mat, const Vec3D& cnt)  
{	
	HaAtom* aptr;
	AtomIteratorGen aitr(this);

	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		double x = aptr->GetX() - cnt.GetX(); 
		double y = aptr->GetY() - cnt.GetY();
		double z = aptr->GetZ() - cnt.GetZ();
		//PrintLog("Xini= %2.3f ", aptr->GetY());
		aptr->SetX(cnt.GetX() + rot_mat(1,1) * x + rot_mat(1,2) * y + rot_mat(1,3) * z);
		aptr->SetY(cnt.GetY() + rot_mat(2,1) * x + rot_mat(2,2) * y + rot_mat(2,3) * z);
		aptr->SetZ(cnt.GetZ() + rot_mat(3,1) * x + rot_mat(3,2) * y + rot_mat(3,3) * z);
		//PrintLog("Xafter = %2.3f\n", aptr->GetY());
	}
	return TRUE;
}

int AtomContainer::TranslateAtoms(const Vec3D& tr_vec )
{
	bool add_x = ( fabs(tr_vec[0]) > DBL_EPSILON );
	bool add_y = ( fabs(tr_vec[1]) > DBL_EPSILON );
	bool add_z = ( fabs(tr_vec[2]) > DBL_EPSILON );
	
	AtomIteratorGen aitr(this);
	HaAtom* aptr;

	for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom())
	{
		if( add_x ) aptr->SetX(aptr->GetX() + tr_vec[0]);
		if( add_y ) aptr->SetY(aptr->GetY() + tr_vec[1]);
		if( add_z ) aptr->SetZ(aptr->GetZ() + tr_vec[2]);
	}
	
	return True;
}

AtomLoadOptions::AtomLoadOptions()
{
	SetStdOptions();
}

AtomLoadOptions::AtomLoadOptions( const AtomLoadOptions& ref)
{
	SetStdOptions(); 
	Copy(ref);
}

AtomLoadOptions::~AtomLoadOptions()
{
	
}

void AtomLoadOptions::SetStdOptions()
{
	SetDefaultMolName("MOL");
	SetCalcBonds(true);
	SetUniqueAtNames(true);
}

void AtomLoadOptions::Copy( const harlem::HashMap &ref)
{
	const AtomLoadOptions* pref = dynamic_cast<const AtomLoadOptions*>(&ref);
	if( pref != NULL )
	{
		calc_bonds        = pref->calc_bonds;
		unique_atom_names = pref->unique_atom_names;
		mol_name_default  = pref->mol_name_default;
	}
}

harlem::HashMap* AtomLoadOptions::clone() const
{
	AtomLoadOptions* p_opt = new AtomLoadOptions(*this);
	return p_opt;
}


AtomSaveOptions::AtomSaveOptions()
{
	SetStdOptions();
}

AtomSaveOptions::AtomSaveOptions( const AtomSaveOptions& ref)
{
	SetStdOptions();
	Copy(ref);
}

AtomSaveOptions::~AtomSaveOptions()
{

}

void AtomSaveOptions::SetStdOptions()
{
	save_selected   = FALSE;
	save_transform  = FALSE;
	save_connect    = TRUE;
	save_atom_ref   = TRUE; 
	save_amber_pdb  = TRUE;
	save_sep_wat_mol = FALSE;

	at_ref_type = HaAtom::ATOMREF_ELEM_NO;
}

void AtomSaveOptions::Copy( const harlem::HashMap& ref )
{	
	using namespace harlem;

	SaveOptions::Copy(ref);

	const AtomSaveOptions* pref = dynamic_cast<const AtomSaveOptions*>(&ref);
	if( pref != NULL )
	{
		save_selected  = pref->save_selected;
		save_connect   = pref->save_connect;
		save_transform = pref->save_transform;
		save_atom_ref  = pref->save_atom_ref;
		save_amber_pdb = pref->save_amber_pdb;

		at_ref_type    = pref->at_ref_type;
	}
}

harlem::HashMap* AtomSaveOptions::clone() const
{
	AtomSaveOptions* p_opt = new AtomSaveOptions(*this);
	return p_opt;
}

CrdSnapshot::CrdSnapshot(AtomContainer* p_at_cont_new, const std::string& name_new )
{
	p_at_cont = p_at_cont_new;
	name = name_new;
	SaveCurrentAtomCrd();
}

CrdSnapshot::~CrdSnapshot()
{

}

void CrdSnapshot::Clear()
{
	name.clear();
	crd.clear();
	pbox.clear();
}

int CrdSnapshot::SaveCrd( const HaVec_double& crd_par )
{
	if( p_at_cont == NULL) return FALSE;
	int na = p_at_cont->GetNAtoms();

	int ncrd = 3*na;

	if( crd_par.size() != ncrd )
	{
		PrintLog("Error in CrdSnapshot::SaveCrd() \n");
		PrintLog(" Number of atoms  %d in atom collection \n does not match the size of the coordinate array = %d \n",
			       na, crd_par.size());
		return FALSE;
	}

	crd = crd_par;
	return TRUE;
}

int CrdSnapshot::SavePBox(const HaVec_double& pbox_par )
{
	if( pbox_par.size() != 0 && pbox_par.size() != 3 && pbox_par.size() != 6 )
	{
		PrintLog(" Error in CrdSnapshot::SavePBox() \n");
		PrintLog(" Invalid size size of periodical box = %d \n", pbox_par.size() );
		return FALSE;
	}
	pbox = pbox_par;
	return TRUE;
}

int CrdSnapshot::SaveCurrentAtomCrd()
{
	if( p_at_cont == NULL) return FALSE;
	
	int ncrd = p_at_cont->GetNAtoms()*3;
	
	crd.resize(ncrd);

	std::auto_ptr<AtomIterator> p_aitr( p_at_cont->GetAtomIteratorPtr() );
	int j = 0;
	HaAtom* aptr;
	for( aptr = p_aitr->GetFirstAtom(); aptr; aptr = p_aitr->GetNextAtom() )
	{
		crd[j] = aptr->GetX(); j++;
		crd[j] = aptr->GetY(); j++;
		crd[j] = aptr->GetZ(); j++;
	}

	MolSet* pmset = dynamic_cast<MolSet*>( p_at_cont );
	if( pmset != NULL && pmset->per_bc->IsSet() ) 
	{
		pbox.resize(6);
		pbox[0] = pmset->per_bc->GetA();
		pbox[1] = pmset->per_bc->GetB();
		pbox[2] = pmset->per_bc->GetC();
		pbox[3] = pmset->per_bc->GetAlpha();
		pbox[4] = pmset->per_bc->GetBeta();
		pbox[5] = pmset->per_bc->GetGamma();
	}
	else
	{
		pbox.clear();
	}
	return TRUE;
}

int CrdSnapshot::SetAtomCrd()
{
	if( !IsValid() )
	{
		PrintLog(" Error in CrdSnapshot::SetAtomCrd() \n");
		PrintLog(" Snapshot is not valid \n");
		return FALSE;
	}

	std::auto_ptr<AtomIterator> p_aitr( p_at_cont->GetAtomIteratorPtr() );
	int j = 0;
	HaAtom* aptr;
	for( aptr = p_aitr->GetFirstAtom(); aptr; aptr = p_aitr->GetNextAtom() )
	{
		aptr->SetX(crd[j]); j++;
		aptr->SetY(crd[j]); j++;
		aptr->SetZ(crd[j]); j++;
	}

	bool set_box = false; 
	MolSet* pmset = dynamic_cast<MolSet*>( p_at_cont );
	if( pmset != NULL && pmset->per_bc->IsSet() ) set_box = true;

	if( pbox.size() == 0 ) set_box = false;

	if( set_box )
	{
		if( pbox.size() == 3 )
		{
			pmset->per_bc->SetBox(pbox[0], pbox[1], pbox[2]);
		}
		else if( pbox.size() == 6 )
		{
			pmset->per_bc->SetBox(pbox[0], pbox[1], pbox[2], pbox[3], pbox[4], pbox[5]);  
		}
	}	
	return TRUE; 
}

int CrdSnapshot::IsValid() const
{
	if( p_at_cont == NULL) return FALSE;

	int na = p_at_cont->GetNAtoms();
	int ncrd = 3*na;

	if( ncrd != crd.size() ) return FALSE;
	return TRUE;
}

int CrdSnapshot::HasPBox() const
{
	if( pbox.size() == 0 ) return FALSE;
	return TRUE;
}

std::string CrdSnapshot::GetName() const
{
	return name;
}

void CrdSnapshot::SetName(const std::string& name_new)
{
	name = name_new;
}

std::string CrdSnapshot::GetDesc() const
{
	return desc;
}

void CrdSnapshot::SetDesc(const std::string& desc_new)
{
	desc = desc_new;
}

int CrdSnapshot::LoadXMLNode( rapidxml::xml_node<>* node, const harlem::HashMap* popt_par )
{
	using namespace rapidxml;

	std::auto_ptr<harlem::HashMap> popt_auto(  (popt_par == NULL) ? new harlem::HashMap() : popt_par->clone() );
	harlem::HashMap* popt = popt_auto.get();

	MolSet* pmset = NULL;
	AtomGroup* p_at_grp = NULL;

	try
	{
		if( popt->has_a("MSET_PTR") ) 
		{
			pmset = boost::any_cast<MolSet*> (popt->get_a("MSET_PTR"));
		}
		if( popt->has_a("ATGRP_PTR") ) 
		{
			p_at_grp = boost::any_cast<AtomGroup*> (popt->get_a("ATGRP_PTR"));
		}

		std::string snap_name;
		std::string snap_desc;

		xml_attribute<>* attr;
		for( attr = node->first_attribute(); attr; attr = attr->next_attribute() )
		{	
			std::string tag_attr = attr->name();
			if( boost::iequals(tag_attr, "name") ) 
			{
				snap_name = attr->value();
				this->SetName(snap_name);
			}
			else if( boost::iequals(tag_attr, "atgrp") ) 
			{
				std::string grp_id = attr->value();
				if( p_at_grp == NULL )
				{
					if( pmset == NULL ) throw std::runtime_error(" molset is not set ");
					p_at_grp = pmset->GetAtomGroupByID( grp_id.c_str() );
					if( !p_at_grp ) throw std::runtime_error(" No atom Group " + grp_id);
				}
			}
		}

		std::string val_str;
		std::vector<std::string> str_vec;

		HaVec_double pbox_snap;
		xml_node<> *node2 = node->first_node();
		for( ; node2; node2 = node2->next_sibling() )
		{
			std::string tag_node2 = node2->name();
			if( boost::iequals(tag_node2, "pbox") )
			{
				val_str = node2->value();
				boost::trim(val_str);
				boost::split(str_vec,val_str,boost::is_any_of(" ,\r\n"),boost::token_compress_on);

				int n = str_vec.size();
				pbox_snap.resize(n);
				for( int i = 0; i < n ; i++ )
				{
					pbox_snap[i] = atof( str_vec[i].c_str() );
					if( i >= 3 && i <= 5 )
					{
						pbox_snap[i] *= DEG_TO_RAD;
					}
				}
			}
			else if( boost::iequals(tag_node2, "desc") )
			{
				snap_desc = node2->value();
				SetDesc( snap_desc );
			}
		}

		val_str = node->value();
		boost::trim(val_str);
		boost::split(str_vec,val_str,boost::is_any_of(" ,\r\n"),boost::token_compress_on);

		int n = str_vec.size();
		int ncrd = 0;
		if( p_at_grp != NULL )
		{
			ncrd = p_at_grp->GetNAtoms()*3;
		}
		else
		{
			ncrd = pmset->GetNAtoms()*3;
		}

		if( ncrd != n ) throw std::runtime_error(" Invalid number of coordinates in xml node "); 

		if( p_at_grp != NULL )
		{
			p_at_cont = p_at_grp;
		}
		else
		{
			p_at_cont = pmset;
		}

		this->crd.resize(ncrd);

		int i;
		for( i = 0; i < ncrd; i++)
		{
			this->crd[i] = atof( str_vec[i].c_str() );
		}
		if( p_at_grp == NULL && (pbox_snap.size() > 0) ) SavePBox( pbox_snap ); 
	}
	catch( const std::exception& ex )
	{
		PrintLog("CrdSnapshot::LoadXMLNode() \n");
		PrintLog("%s\n",ex.what());
		return FALSE;
	}
	return TRUE;
}

int CrdSnapshot::OnDelAtoms( AtomContainer& at_del_par )
{
	AtomIteratorGen aitr_del(&at_del_par );
	HaAtom* aptr;

    std::set<HaAtom*> del_atoms;

	for(aptr = aitr_del.GetFirstAtom(); aptr; aptr = aitr_del.GetNextAtom())
	{
		del_atoms.insert(aptr);
	}

	if( !IsValid() ) return FALSE;

	HaVec_double crd_new;
	crd_new.reserve(crd.size());
	AtomIteratorGen aitr(p_at_cont);
	int i = 0;
	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
	{
		if( del_atoms.count(aptr) == 0 ) 
		{
			for( int k=0; k < 3; k++)
			{
				crd_new.push_back(crd[i]); i++;
			}
		}
		else
		{
			i += 3;
		}
	}
	crd = crd_new;
	return TRUE;
}

CrdSnapshotIterator::CrdSnapshotIterator( MolSet* pmset )
{
	first_call = true;
	itr_curr = pmset->crd_snapshots.begin();
	itr_end  = pmset->crd_snapshots.end();
}

CrdSnapshotIterator::CrdSnapshotIterator(const CrdSnapshotIterator& ref)
{
	first_call = ref.first_call;
	itr_curr   = ref.itr_curr;
	itr_end    = ref.itr_end;
}

CrdSnapshotIterator::~CrdSnapshotIterator()
{

}

CrdSnapshotIterator& CrdSnapshotIterator::operator=( const CrdSnapshotIterator& ref)
{
	first_call = ref.first_call;
	itr_curr   = ref.itr_curr;
	itr_end    = ref.itr_end;

	return(*this);
}

CrdSnapshot* CrdSnapshotIterator::next()
{
	if(!first_call )
	{
		itr_curr++;
	}
	else
	{
		first_call = false;
	}
	if( itr_curr == itr_end ) 
	{
		throw std::out_of_range("Stop Coordinate Snapshot Iterations");
		return NULL;
	}
	return (*itr_curr);
}

CrdSnapshotIterator CrdSnapshotIterator::__iter__() const
{
	return (*this);
}


int AtomLoadOptions::ConvertResNames()
{
	if (!this->has_i("CONVERT_RES_NAMES")) return 0;
	return this->get_i("CONVERT_RES_NAMES");
}

void AtomLoadOptions::SetConvertResNames( int convert_res_names_opt )
{
	if (convert_res_names_opt < 0 || convert_res_names_opt > 1) this->set_i("CONVERT_RES_NAMES", 0);
	this->set_i("CONVERT_RES_NAMES", convert_res_names_opt);
}


MutationMap::MutationMap(MolSet* pmset_par)
{
	pmset = pmset_par;
}

MutationMap::~MutationMap()
{

}

void MutationMap::Clear() noexcept
{
	this->pmol1 = NULL;
	this->pmol2 = NULL;
	this->atom_atom_map.clear();
}

bool MutationMap::IsValid()
{
	try
	{
		if (pmol1 == NULL) throw std::runtime_error(" First molecule is not defined ");
		if (pmol2 == NULL) throw std::runtime_error(" Second molecule is not defined ");
		if ( atom_atom_map.empty() ) throw std::runtime_error(" Atom-Atom correspondence map is empty ");
	}
	catch (const std::exception& ex)
	{
		PrintLog("Mutation Map is no valid %s \n", ex.what());
		return false;
	}
	return true;
}

int MutationMap::LoadArbalestMutMap(std::string fname) 
{
	if(!(pmset->GetNMol() == 2))
	{
		PrintLog("Error in MutationMap::LoadArbalestMutMap() \n");
		PrintLog("Number of Molecules in Molecular Set = %d  is not equal 2", pmset->GetNMol());
		return FALSE;
	}
	this->Clear();

	pmol1 = pmset->GetMolByIdx(0);
	pmol2 = pmset->GetMolByIdx(1);
	
	try
	{ 
		TiXmlDocument doc;
		bool bres = doc.LoadFile(fname.c_str());

		const TiXmlElement* root_element;
		const TiXmlElement* data_element;

		if (!bres) throw std::runtime_error("Error to load an XML file with a mutation map ");
	
		root_element = doc.FirstChildElement();
		if (root_element == NULL)  throw std::runtime_error("No ROOT Element");
		
		const TiXmlElement* resA_elem = root_element->FirstChildElement("ResA");
		if( resA_elem == NULL ) throw std::runtime_error("No ResA Element");
		std::string res_name_a = resA_elem->GetText();

		const TiXmlElement* resB_elem = root_element->FirstChildElement("ResB");
		if (resB_elem == NULL) throw std::runtime_error("No ResB Element");
		std::string res_name_b = resB_elem->GetText();
		
		HaResidue* pres = NULL;
		HaResidue* pres_a = NULL;
		HaResidue* pres_b = NULL;
		
		ResidueIteratorMolecule ritr_a(pmol1);
		for (pres = ritr_a.GetFirstRes(); pres; pres = ritr_a.GetNextRes())
		{
			std::string res_name = pres->GetName();
			if (res_name == res_name_a)
			{
				pres_a = pres;
				break;
			}
		}
		if( pres_a == NULL ) throw std::runtime_error( "Can not find residue in MOL1 the name = " + res_name_a  );

		ResidueIteratorMolecule ritr_b(pmol2);
		for (pres = ritr_b.GetFirstRes(); pres; pres = ritr_b.GetNextRes())
		{
			std::string res_name = pres->GetName();
			if (res_name == res_name_b  )
			{
				pres_b = pres;
				break;
			}
		}
		if ( pres_b == NULL ) throw std::runtime_error("Can not find residue in MOL2 with the name " + res_name_b);

		const TiXmlElement* atoms_mapping_elem = root_element->FirstChildElement("AtomsMapping");
		if (atoms_mapping_elem == NULL ) throw std::runtime_error("Can not find AtomsMapping element");

		const TiXmlElement* atoms_elem = atoms_mapping_elem->FirstChildElement();
		for (; atoms_elem != NULL; atoms_elem = atoms_elem->NextSiblingElement())
		{
			const char* sattr_a = atoms_elem->Attribute("A");
			const char* sattr_b = atoms_elem->Attribute("B");
			if (sattr_a != NULL && sattr_b != NULL)
			{
				HaAtom* aptr_a = pres_a->GetAtomByName(sattr_a);
				HaAtom* aptr_b = pres_b->GetAtomByName(sattr_b);
				if (aptr_a == NULL) PrintLog(" Can not find Atom Name %s in ResA ", sattr_a);
				if (aptr_b == NULL) PrintLog(" Can not find Atom Name %s in ResB ", sattr_b);
				if (aptr_a != NULL && aptr_b != NULL)
				{
					this->atom_atom_map[aptr_a] = aptr_b;
				}
			}
		}
	}
	catch(const std::exception & ex)
	{
		PrintLog(" Error in MutationMap::LoadArbalestMutMap \n");
		PrintLog(" Failed to Load Mutation Map XML File %s \n", fname.c_str());
		PrintLog("%s\n", ex.what());
		return FALSE;
	}
	return TRUE;

}

int MutationMap::SaveArbalestMutMap(std::string fname)
{
	if (!this->IsValid()) return FALSE;
	std::ofstream os(fname);
	if (os.fail()) return FALSE;
	os << "<?xml version=\"1.0\" encoding=\"windows - 1251\"?>" << std::endl;
	std::string title = pmol1->GetObjName();
	title += "-}";
	title += pmol2->GetObjName();
	os << "<MutationMap Title=\"" + title + "\">" << std::endl;

	HaResidue* pres_1 = pmol1->GetFirstChain()->GetFirstRes();
	std::string res_name_1 = pres_1->GetName();
	HaResidue* pres_2 = pmol2->GetFirstChain()->GetFirstRes();
	std::string res_name_2 = pres_2->GetName();

	os << "    <ResA>" << res_name_1 << "</ResA>" << std::endl;
	os << "    <ResB>" << res_name_2 << "</ResB>" << std::endl;
	os << "    <AtomsMapping>" << std::endl;
	
	auto mitr = atom_atom_map.begin();
	for (; mitr != atom_atom_map.end(); mitr++)
	{
		HaAtom* aptr1 = (*mitr).first;
		HaAtom* aptr2 = (*mitr).second;
		os << "        <Atoms A=\"" << aptr1->GetName() << "\"  B=\"" << aptr2->GetName() << "\"/> " << std::endl;
	}
	os << "    </AtomsMapping>" << std::endl;
	os << "</MutationMap>" << std::endl;
	return TRUE;
}