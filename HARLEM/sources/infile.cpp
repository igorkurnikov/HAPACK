/*! \file infile.cpp
 
   Functions to read molecular structure in different formats
  
   \author Igor Kurnikov 
   \date 1998-2002
 
   RasMol2 Molecular Graphics
   Roger Sayle, August 1995
   Version 2.6
 */
#include "haconst.h"

#include <mpi.h>

#include <malloc.h>

#ifndef sun386
#include <stdlib.h>
#endif

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include <string>
#include <ctype.h>
#include <stdio.h>
#include <math.h>

#include "rapidxml.hpp"

#define INFILE
#include "hastl.h"
#include "haatom.h"
#include "habond.h"
#include "hamolecule.h"
#include "moleditor.h"
#include "haqchem.h"
#include "etcoupl.h"
#include "hamolview.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "hamolmech.h"
#include "gaufile.h"
#include "abstree.h"
#include "command.h"


#if !defined(_MSC_VER)
#  include <sys/types.h>
#  include <sys/time.h>
#else
#  include <time.h>
#endif

typedef struct {
	  char src[4];
	  char dst[4];
	  } ConvTable;

#define MAXALCATOM   5
static ConvTable AlcAtomTable[MAXALCATOM] = {
    { { 'S', 'O', '2', ' ' }, { ' ', 'S', '2', ' ' } },  /*  1 */
    { { 'C', 'A', 'R', ' ' }, { ' ', 'C', ' ', ' ' } },  /*  2 */
    { { 'N', 'A', 'R', ' ' }, { ' ', 'N', ' ', ' ' } },  /*  3 */
    { { 'N', 'A', 'M', ' ' }, { ' ', 'N', ' ', ' ' } },  /*  4 */
    { { 'N', 'P', 'L', '3' }, { ' ', 'N', '3', ' ' } },  /*  5 */
				 };

static char PDBInsert;

const  int  MAXBUF=5000;
static char buf[MAXBUF];

void HaMolecule::UpdateFeature(SecStructElement* ptr,int mask )
{
    HaChain  *chain;
    HaResidue  *group;
	
    ChainIteratorMolecule ch_itr(this);
	for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
	{
        if( chain->ident == ptr->chain )
        {
			ResidueIteratorChain ritr_ch(chain);
			for(group = ritr_ch.GetFirstRes(); group; group = ritr_ch.GetNextRes())
			{
				if( (ptr->init <= group->serno)  && (group->serno <= ptr->term))
					group->struc |= mask;
			}
			continue;
        }
	}
}
 
void HaMolecule::ProcessFeatures()
{
	list<SecStructElement>::iterator fitr;

    structsource = SourcePDB;
	
	for(fitr= Features.begin(); fitr != Features.end(); fitr++)
	{
		if( (*fitr).type == FeatHelix )
		{   
			UpdateFeature( &(*fitr), HelixFlag );
		} 
		else if( (*fitr).type == FeatSheet )
		{   
			UpdateFeature( &(*fitr), SheetFlag );
		} 
		else /* FeatTurn */
		{   
			UpdateFeature( &(*fitr), TurnFlag );
		}
    }
}
 

/*==============================*/
/* Molecule File Format Parsing */
/*==============================*/

static int res_ser_no_old = 0;

void MolSet::ProcessPDBAtom(const std::string& line, int heta, IntPtrMap& id_at_map, 
							  HaMolecule* pMol, HaChain* &pch_cur, HaResidue* &pres_cur)
{
	HaBond  *bptr;
	HaAtom  *ptr;
	double dx,dy,dz;
	int temp,serno;

	using boost::lexical_cast;
	using boost::bad_lexical_cast;
	using boost::trim_copy;
	using std::string;

	try
	{
		try 
		{
			dx = lexical_cast<double>( trim_copy(line.substr(30,8)) );
			dy = lexical_cast<double>( trim_copy(line.substr(38,8)) );
			dz = lexical_cast<double>( trim_copy(line.substr(46,8)) );
		}
		catch(bad_lexical_cast&) 
		{
			throw std::runtime_error("Invalid coordinates");
		}

		if( (fabs(dx - 9999.0) < 0.0001) && (fabs(dy - 9999.0) < 0.0001) && (fabs(dz - 9999.0) < 0.0001) )
		{
			throw std::runtime_error(" XPLOR Pseudo atom record is ignored "); 
		}

//		if( (fabs(dx - 0.0) < 0.0001) && (fabs(dy - 0.0) < 0.0001) && (fabs(dz - 0.0) < 0.0001) )
//		{
//			throw std::runtime_error(" zero coordinate pseudo atom record is ignored "); 
//		}

		string serno_str = line.substr(22,4); // Reading Residue ID
		string temp_str;
		int i;
		for( i = 0; i < serno_str.size(); i++)
		{
			if( isdigit(serno_str[i]) || serno_str[i] == '-' ) temp_str += serno_str[i]; 
		}
		serno_str = temp_str;

		try { serno = lexical_cast<int>( serno_str ); }
		catch( bad_lexical_cast&){ serno = 0;}

//		if(serno == 0) 
//		{
//			if(pres_cur == NULL)
//			{
//				serno = 1;
//			}
//			else
//			{
//				serno = pres_cur->serno;
//			}
//		}

		if( !pres_cur                       // no active residue 
			|| (res_ser_no_old != serno )    // new residue number in the stream 
			|| (pch_cur->ident != line[21])
			|| (PDBInsert != line[26]) )
		{
			PDBInsert = line[26];
			if( !pch_cur || (pch_cur->ident != line[21]) )
			{   
				pch_cur = pMol->AddChain( line[21] );
			}

			res_ser_no_old = serno;
			pres_cur = pch_cur->AddResidue(serno);

			if( pres_cur == NULL )
			{
				PrintLog(" Try to find unique number for a residue \n" );
				serno    = pch_cur->GetUniqResSerNo(false);
				pres_cur = pch_cur->AddResidue(serno);
			}

			pres_cur->SetName( line.substr(17,4) );
		}

		int at_serno;

		try { 
			at_serno= lexical_cast<int>( trim_copy(line.substr(6,5)) ); 
		}
		catch( bad_lexical_cast&){ 
			throw std::runtime_error(" Invalid serial atom number ");
		}

		string tmp = line.substr(12,4);

		// Deal with alternative atom positions
		// Will Load only the first alternative atom position

		if( line[16] != ' ' && line[16] != 'A')
		{
			throw std::runtime_error(" Omit Alternative Atom position ");
		}

		ptr = pres_cur->AddNewAtom();

		ptr->SetName(tmp.c_str());
		string at_name = ptr->GetName();
		if( !at_name.empty() )
		{
			ptr->SetElemNo( HaAtom::GetElemNoFromChar(at_name[0]) );
		}
		//	ptr->SetElemNo( HaAtom::GetElemNoFromName( ptr->GetName(), pres_cur ) );

		id_at_map[at_serno] = (void*)ptr;

		ptr->SetX_Ang(dx);
		ptr->SetY_Ang(dy);
		ptr->SetZ_Ang(dz);
		
		try {
			ptr->tempf = lexical_cast<double>( trim_copy(line.substr(60,7)) );
		}
		catch( bad_lexical_cast&){ 
			ptr->tempf = 0.0;
//			throw std::runtime_error(" Invalid Temp field ");
		}

		/* Handle CONCORD PDB Files */
		if( (line[12]==' ') && islower(line[14]) && !trim_copy(line.substr(15,7)).empty() )
		{
			ptr->SetName( line.substr(13,4).c_str() );
			ptr->SetElemNo( HaAtom::GetElemNoFromName( ptr->GetName(), pres_cur ) );
		}

		if( heta ) ptr->flag |= HeteroFlag;
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Warning in MolSet::ProcessPDBAtom()  Reading line: \n" );
		PrintLog("%s \n",line.c_str());
		PrintLog("%s\n",ex.what());
	}
}
 
int HaMolecule::SetTermResNames()
{
	HaChain* chain;
	ChainIteratorMolecule ch_itr(this);
	for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
	{
		ResidueIteratorChain ritr_ch(chain);
		HaResidue* pres = ritr_ch.GetFirstRes();
		if(!pres) continue;
		if(!pres->IsAmino()) continue;
		pres->SetNameModifier("NT");
		HaResidue* plres = NULL;
		for( ; pres != NULL; pres = ritr_ch.GetNextRes())
		{
			if(pres->IsAmino()) plres = pres;
		}
		if( plres != NULL && plres != pres)
		{
			plres->SetNameModifier("CT");
		}	
	}
	return TRUE;
}

int HaMolecule::SetHISNames()
{
	HaChain* chain;
	ChainIteratorMolecule ch_itr(this);
	for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
	{
		ResidueIteratorChain ritr_ch(chain);
		HaResidue* pres = ritr_ch.GetFirstRes();
		for( ; pres; pres = ritr_ch.GetNextRes() )
		{
			if(!pres) continue;
			if(!pres->IsAmino()) continue;
			std::string resname = pres->GetName();
			std::string name_mod = pres->GetNameModifier();
			if( resname != "HIS" ) continue;
			if( name_mod == "EPSILON" || name_mod == "EPSILON_NT" || name_mod == "EPSILON_CT" ) continue;
			if( name_mod == "PROT"    || name_mod == "PROT_NT"    || name_mod == "PROT_CT") continue;

			HaAtom* phe2 = pres->GetAtomByName("HE2");
			HaAtom* phd1 = pres->GetAtomByName("HD1");

			if( phe2 && !phd1 )
			{
				if( name_mod == "NT") pres->SetNameModifier("EPSILON_NT");
				else if( name_mod == "CT") pres->SetNameModifier("EPSILON_CT");
				else pres->SetNameModifier("EPSILON");
			}

			if( phe2 && phd1 )
			{
				if( name_mod == "NT") pres->SetNameModifier("PROT_NT");
				else if( name_mod == "CT") pres->SetNameModifier("PROT_CT");
				else pres->SetNameModifier("PROT");
			}
		}
	}
	return TRUE;
}



int HaMolecule::SetCysBridgeNames()
{
	AtomIteratorMolecule aitr(this);
	HaAtom* aptr;
	for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom() )
	{
		if( strcmp(aptr->GetName(),"SG") != 0 ) continue;
		HaAtom::BondIterator bitr = aptr->Bonds_begin();
		for( ; bitr != aptr->Bonds_end(); ++bitr)
		{
			HaBond* bptr = (*bitr);
			if( bptr->GetFirstAtom() > bptr->GetSecondAtom() ) continue;

			HaAtom* aptr2 = bptr->GetFirstAtom() == aptr ? bptr->GetSecondAtom() : bptr->GetFirstAtom();
			if( strcmp(aptr2->GetName(),"SG") != 0 ) continue;

			HaResidue* pres1 = bptr->srcatom->GetHostRes(); 
			HaResidue* pres2 = bptr->dstatom->GetHostRes(); 

			if(!pres1->IsCysteine())continue;
			if(!pres2->IsCysteine())continue;

			std::string fname = pres1->GetFullName();

			if(fname == "CYS")
			{
				pres1->SetNameModifier("UNPROT");
			}
			else if(fname == "CYS#CT")
			{
				pres1->SetNameModifier("UNPROT_CT");
			}
			else if(fname == "CYS#NT")
			{
				pres1->SetNameModifier("UNPROT_NT");
			}

			fname = pres2->GetFullName();

			if(fname == "CYS")
			{
				pres2->SetNameModifier("UNPROT");
			}
			else if(fname == "CYS#CT")
			{
				pres2->SetNameModifier("UNPROT_CT");
			}
			else if(fname == "CYS#NT")
			{
				pres2->SetNameModifier("UNPROT_NT");
			}
		}
	}

	return TRUE;
}

int MolSet::LoadPDBFile(const char* fname , int flag )
{
	using boost::lexical_cast;
	using boost::trim_copy;

	p_mol_editor = this->GetMolEditor(true);
	std::ifstream is(fname);

	if( is.fail() )
	{   
        PrintLog("\nError: File  %s not found!\n",fname );
        return( FALSE );
	}

    SecStructElement  *ptr;
    int srcatm, dstatm;
    int i,ignore;
	std::string str,line;

	HaMolecule* pMol = AddNewMolecule();
	HaChain*    pch_cur  = NULL;
	HaResidue*  pres_cur = NULL;

	std::string molname = harlem::GetPrefixFromFullName(fname);
	pMol->SetObjName(molname.c_str());

    ignore = False;

	IntPtrMap id_at_map;

    for(;;)
    {  
		try
		{
			std::getline(is,line);
			line.resize(81,' ');
			if( is.fail() ) break;

			if( line[0] == 'A' )
			{   
				if( !ignore && !strncmp("ATOM",line.c_str(),4) )
					ProcessPDBAtom( line, False, id_at_map, pMol, pch_cur, pres_cur);
			} 
			else 
			{
				switch( line[0] )
				{   
				case('C'):    
					if( line.size() >= 4 && !strncmp("CONE",line.c_str(),4) )
					{   
						if( ignore || flag ) continue;
						std::string atom_id_str = trim_copy(line.substr(6, 5));
						if (atom_id_str.empty()) continue;
						srcatm = lexical_cast<int>( atom_id_str );
						if( srcatm )
						{
							for( i=11; i<=36 && line[i]; i+=5 )
							{   
								atom_id_str = trim_copy(line.substr(i, 5));
								if (atom_id_str.empty()) continue;
								dstatm = lexical_cast<int>( atom_id_str );
								if( dstatm && (dstatm>srcatm) )
								{
									IntPtrMap::iterator mitr;
									mitr = id_at_map.find(srcatm);
									if(mitr == id_at_map.end())
										continue;
									HaAtom* aptr1 = (HaAtom*)(*mitr).second;
									mitr = id_at_map.find(dstatm);
									if(mitr == id_at_map.end())
										continue;
									HaAtom* aptr2 = (HaAtom*)(*mitr).second;
									AddBond(aptr1,aptr2);
								}
							}
						}
					}  
					else if( line.size() >= 4 && !strncmp("CRYS",line.c_str(),4) )
					{   
						per_bc->spacegroup = trim_copy( line.substr(55,11) );

						double px = lexical_cast<double>( trim_copy(line.substr(6,9)) );
						double py = lexical_cast<double>( trim_copy(line.substr(15,9)) );
						double pz = lexical_cast<double>( trim_copy(line.substr(24,9)) );

						double p_alpha = DEG_TO_RAD* lexical_cast<double>( trim_copy(line.substr(33,7)) );
						double p_beta  = DEG_TO_RAD* lexical_cast<double>( trim_copy(line.substr(40,7)) );
						double p_gamma = DEG_TO_RAD* lexical_cast<double>( trim_copy(line.substr(47,7)) );

						per_bc->SetBox(px,py,pz,p_alpha,p_beta,p_gamma);
					} 

				case('E'):    
					if( line.size() >= 4 && !strncmp("ENDM",line.c_str(),4) )
					{   /* break after single model??? */
						if( flag )
						{   
							pres_cur = NULL;
							pch_cur  = NULL;
						} 
						else 
							ignore = True;

					} 
					else if( line.size() >= 3 && !strncmp("END",line.c_str(),3) )
					{
				// Treat END same a TER 
						pres_cur =  NULL;
						pch_cur  =  NULL;
						break;
					}
				case('H'):    
					if( line.size() >= 4 && !strncmp("HETA",line.c_str(),4) )
					{   
						if( !ignore ) ProcessPDBAtom(line,True, id_at_map, pMol, pch_cur, pres_cur);
					} 
					else if( line.size() >= 4 && !strncmp("HELI",line.c_str(),4) )
					{   
						if( ignore ) continue;

						/* Remaining HELIX record fields   */
						/* 38-39 .... Helix Classification */
						/* 31 ....... Same Chain as 19?    */
						ptr = pMol->AddFeature();
						ptr->type = FeatHelix;
						ptr->chain = line[19];
						ptr->init = lexical_cast<int>( trim_copy(line.substr(21,4)) );
						ptr->term = lexical_cast<int>( trim_copy(line.substr(33,4)) );

					} 
					else if( line.size() >= 4 && !strncmp("HEAD",line.c_str(),4) )
					{   
						pMol->classification = trim_copy( line.substr(10,40) );
						pMol->identcode      = trim_copy( line.substr(62,4) );
					}
					break;

				case('S'):    
					if( line.size() >= 4 && !strncmp("SHEE",line.c_str(),4) )
					{   
						if( ignore ) break;
						/* Remaining SHEET record fields   */
						/* 38-39 .... Strand Parallelism   */
						/* 32 ....... Same Chain as 21?    */
						ptr = pMol->AddFeature();
						ptr->type = FeatSheet;
						ptr->chain = line[21];
						ptr->init = lexical_cast<int>( trim_copy(line.substr(22,4)) );
						ptr->term = lexical_cast<int>( trim_copy(line.substr(33,4)) );
					}
					break;

				case('T'):    
					if( line.size() >= 4 && !strncmp("TURN",line.c_str(),4) )
					{   
						if( ignore ) continue;

						ptr = pMol->AddFeature();
						ptr->type = FeatTurn;
						ptr->chain = line[19];
						ptr->init = lexical_cast<int>( trim_copy(line.substr(20,4)) );
						ptr->term = lexical_cast<int>( trim_copy(line.substr(31,4)) );
					} 
					else if( line.size() >= 3 && !strncmp("TER",line.c_str(),3) )
					{   
						pres_cur = NULL;
						pch_cur  = NULL;					
					}
					break;
				}
			}
		}
		catch( const std::exception& ex )
		{
			PrintLog(" Warning in MolSet::LoadPDBFile() \n");
			PrintLog(" Reading Line %s\n",line.c_str());
			PrintLog("%s\n",ex.what());
		}
	}
	pMol->FixChainsIdent();
	pMol->SetTermResNames();
	//	if( p_load_opt_default->ToCalcBonds() )
//	{
		p_mol_editor->CreateCovBonds(pMol);
//	}

	pMol->SetCysBridgeNames();
	pMol->SetHISNames();

    if( !pMol->Features.empty() ) pMol->ProcessFeatures();

	if( p_load_opt_default->UniqueAtNames() )
	{
		pMol->SetUniqueAtomNames();
	}

    return( TRUE );
}

bool HaMolecule::FixChainsIdent()
{
	if(GetNChains() <= 1) return true;
	
	list<HaChain>::iterator citr1;
	list<HaChain>::iterator citr2;

	set<char, less<char> > not_used_id;
	set<char, less<char> > used_id;

	set<char, less<char> >::iterator iditr;

	char id;

	not_used_id.insert(' ');
	for(id ='A'; id <= 'Z'; id++)
	{
		not_used_id.insert(id);
	}

	for(citr1 = Chains.begin(); citr1 != Chains.end(); citr1++)
	{
		id = (*citr1).ident;
		if( used_id.find( id ) ==  used_id.end() ) // unique id
		{	
			iditr = not_used_id.find( id );
			if( iditr != not_used_id.end() )
				not_used_id.erase(iditr);
			used_id.insert(id);
		}
		else             // not unique id
		{
			if( !not_used_id.empty() )
			{
				iditr = not_used_id.begin();
				id = *iditr;
				(*citr1).ident = id;
				not_used_id.erase(iditr);
			}
		}
	}
	return true;
}


int MolSet::LoadMDLFile(const char* fname )
{
	using boost::trim_copy;
	using boost::lexical_cast;

	MolEditor* p_mol_editor = this->GetMolEditor(true);

	std::ifstream is(fname);

	if( is.fail() )
	{   
        PrintLog("\nError: File  %s not found!\n",fname );
        return( FALSE );
	}

    HaAtom  *ptr;
	
    int i,type;
    int atoms, bonds;
    int srcatm,dstatm;
    char *cptr;
	std::string str,line;
	
	HaMolecule* pMol=AddNewMolecule();
	HaChain*    pch_cur  = NULL;
	HaResidue*  pres_cur = NULL;

	std::string molname = harlem::GetPrefixFromFullName(fname);
	pMol->SetObjName(molname.c_str());
	
	std::getline(is,line);	
	boost::trim(line);
    pMol->SetObjName(line.c_str());
	
    std::getline(is,line); /* Program Details */
    std::getline(is,line); /* Comments */
	
    std::getline(is,line);
	try {
		atoms = lexical_cast<int>( trim_copy(line.substr(0,3)) );
		bonds = lexical_cast<int>( trim_copy(line.substr(3,3)) );
	}
	catch(boost::bad_lexical_cast&) 
	{
		atoms = 0;
	}
	
	IntPtrMap id_at_map;

    if( !atoms ) 
	{
		PrintLog(" Num atoms = 0 \n");
		return( False );
	}
	
    pres_cur = pMol->AddChainAndResidue();
    for( i=1; i<=atoms; i++ )
    {   
		std::getline(is,line);
		if( is.fail()) break;
		
		try
		{
			ptr = pres_cur->AddNewAtom();

			int idx = 31;
			while( line[idx] == ' ' ) idx++;
			std::string at_name = line.substr(idx,4);

			ptr->SetName(at_name.c_str());
			ptr->SetElemNo( HaAtom::GetElemNoFromName( ptr->GetName(), pres_cur ) );

			idx = 0;
			idx = lexical_cast<int>( trim_copy(line.substr(36,3)) );

			switch( idx )
			{   
			case(1):  ptr->tempf =  3.0;  break;
			case(2):  ptr->tempf =  2.0;  break;
			case(3):  ptr->tempf =  1.0;  break;
			case(5):  ptr->tempf = -1.0;  break;
			case(6):  ptr->tempf = -2.0;  break;
			case(7):  ptr->tempf = -3.0;  break;
			default:  ptr->tempf = 0.0;
			}
			id_at_map[i] = (void*)ptr;

			ptr->SetX_Ang( lexical_cast<double>( trim_copy(line.substr( 0,10)) )/10000.0  );
			ptr->SetY_Ang( lexical_cast<double>( trim_copy(line.substr( 10,10)) )/10000.0 );
			ptr->SetZ_Ang( lexical_cast<double>( trim_copy(line.substr( 20,10)) )/10000.0 );
		}
		catch( const std::exception& ex )
		{
			PrintLog(" Error in MolSet::LoadMDLFile()   line:\n");
			PrintLog(" %s\n", line.c_str());
			PrintLog(" %s\b", ex.what());
		}
	}

	for( i=0; i<bonds; i++ )
	{   
		std::getline(is,line);
		if( is.fail()) break;

		try
		{
			srcatm = lexical_cast<int>( trim_copy(line.substr(0,3)) );
			dstatm = lexical_cast<int>( trim_copy(line.substr(3,3)) );
			type =   lexical_cast<int>( trim_copy(line.substr(6,3)) );
			HaBond* pbnd = AddBond((HaAtom*)id_at_map[srcatm],(HaAtom*)id_at_map[dstatm]);

			if( type==2 ) pbnd->SetDouble();            /* DOUBLE */
			else if( type==3 ) pbnd->SetTriple();       /* TRIPLE */
			else if( type==4 ) pbnd->SetAromatic();     /* AROMATIC */
		}
		catch( const std::exception& ex )
		{
			PrintLog(" Error in MolSet::LoadMDLFile()   line:\n");
			PrintLog(" %s\n", line.c_str());
			PrintLog(" %s\b", ex.what());
		}
	}
	

	if( p_load_opt_default->ToCalcBonds() )
	{
		p_mol_editor->CreateCovBonds(pMol);
	}
	if( p_load_opt_default->UniqueAtNames() )
	{
		pMol->SetUniqueAtomNames();
	}
	
	return( True );
}

int MolSet::SetCoordFromFile(const char* fname, int iform)
{
	MolSet axx_mset;
	int ires;
	if(iform == FormatGUESS || iform == FormatPDB)
	{
		ires = axx_mset.LoadPDBFile(fname);
	}
	else if(iform == FormatHarlem)
	{
        ires = axx_mset.LoadHarlemFile(fname);
	}
	else if(iform == FormatXYZ)
	{
        ires = axx_mset.LoadXYZFile(fname);
	}
	else if(iform == FormatHIN)
	{
        ires = axx_mset.LoadHINFile(fname);
	}

	if( !ires)
	{
		PrintLog("ERROR loading file %s \n",fname);
        return FALSE;
	}

	AtomIteratorMolSet aitr1(this);
	AtomIteratorMolSet  aitr2(&axx_mset);

	HaAtom* aptr1 = aitr1.GetFirstAtom();
	HaAtom* aptr2 = aitr2.GetFirstAtom();

	int na1 = GetNAtoms();
	int na2 = axx_mset.GetNAtoms();

	if( na1 != na2)
	{
		PrintLog("WARNING: MolSet::SetCoordFromFile() \n");
		PrintLog("Number of atoms in the molecule read %d \n",na2);
		PrintLog("Is not equal to the number of atoms in current mol set %d\n",na1);
	}
	while( aptr1 != NULL && aptr2 != NULL)
	{
	   aptr1->SetX(aptr2->GetX());
	   aptr1->SetY(aptr2->GetY());
	   aptr1->SetZ(aptr2->GetZ());

	   aptr1 = aitr1.GetNextAtom();
	   aptr2 = aitr2.GetNextAtom();
	}
	AnnounceGeomChange();
	return TRUE;
}

int MolSet::SetCrdFromArray( const HaVec_double& crd_arr )
{
	int na = GetNAtoms();
	if( crd_arr.size() != 3*na )
	{
		PrintLog(" Error in MolSet::SetCrdFromArray() \n");
		PrintLog(" size of coordinate array = %d  not equal to 3* natoms = %d \n", crd_arr.size(), 3*na );
		return FALSE;
	}
	int i = 0;
	AtomIteratorMolSet aitr(this);
	HaAtom* aptr;
	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
	{
		aptr->SetX( crd_arr[i] ); i++;
		aptr->SetY( crd_arr[i] ); i++;
		aptr->SetZ( crd_arr[i] ); i++;
	}
	return TRUE;
}

int MolSet::LoadXYZFile(const char* fname, const AtomLoadOptions* popt_arg )
{
	ifstream is_f(fname);
	if( !is_f.good() )
	{
		PrintLog(" Error in MolSet::LoadXYZFile() \n");
		PrintLog(" Error to open file %s for reading \n",fname);
		return FALSE;
	}
	std::string mol_name = harlem::GetPrefixFromFullName(fname);

	std::auto_ptr<AtomLoadOptions> popt_auto( popt_arg != NULL ? (AtomLoadOptions*) popt_arg->clone() : (AtomLoadOptions*) p_load_opt_default->clone() );
	AtomLoadOptions* popt = popt_auto.get();

	popt->SetDefaultMolName(mol_name);

	int ires = LoadXYZStream( is_f, popt );

	return ires;
}

int MolSet::LoadHINFile(const char* fname, const AtomLoadOptions* popt_arg )
{
	ifstream is_f(fname);
	if( !is_f.good() )
	{
		PrintLog(" Error in MolSet::LoadHINFile() \n");
		PrintLog(" Error to open file %s for reading \n",fname);
		return FALSE;
	}
	std::string mol_name = harlem::GetPrefixFromFullName(fname);

	std::auto_ptr<AtomLoadOptions> popt_auto( popt_arg != NULL ? (AtomLoadOptions*) popt_arg->clone() : (AtomLoadOptions*) p_load_opt_default->clone() );
	AtomLoadOptions* popt = popt_auto.get();

	popt->SetDefaultMolName(mol_name);

	int ires = LoadHINStream( is_f, popt );

	return ires;
}

int MolSet::LoadXYZStream( std::istream& is_arg, const AtomLoadOptions* p_opt_arg )
{
	std::auto_ptr<AtomLoadOptions> p_opt_auto( p_opt_arg == NULL ? (AtomLoadOptions*) p_load_opt_default->clone() : (AtomLoadOptions*) p_opt_arg->clone() );
	AtomLoadOptions* p_opt = p_opt_auto.get();
	
	MolEditor* p_mol_editor = this->GetMolEditor(true);
 
	if( !is_arg.good() )
	{   
        PrintLog("Error in MolSet::LoadXYZStream() \n");
        PrintLog("Can not read from input stream \n");
        return( FALSE );
	}

	using harlem::IsFloat;
	using harlem::IsInt;
	using boost::lexical_cast;
	using std::string;

    int na;
    HaAtom* aptr;
    int i;
	std::string line;

	HaMolecule* pMol = AddNewMolecule();
	pMol->SetObjName(p_opt->GetDefaultMolName().c_str() );
	HaChain*   pch_cur  = NULL;
	HaResidue* pres_cur = NULL; 

	AtomGroup atoms_new;
	std::vector<std::vector<int> > bonds; 
	std::string msg;

	try
	{
		pres_cur = pMol->AddChainAndResidue();
		pch_cur  = pres_cur->GetHostChain();
		int nat_res = 0;
		int res_idx = 1;
		int iat = 0;
		int idx_line = -1;
		int idx_chk  = -1;
		int na_max = 10000000;

		std::vector<string> str_arr;

		while(1)
		{   
			std::getline(is_arg,line);
			if( is_arg.fail() ) break;
			idx_line++;

			boost::trim(line);
			boost::split(str_arr, line, boost::is_any_of("\t "),boost::token_compress_on);
			
			if( str_arr.size() == 0 ) continue;
			if( str_arr.size() < 3  ) 
			{
				if( idx_line == 0 && harlem::IsInt(str_arr[0])) 
				{
					na = lexical_cast<int>( str_arr[0]);
				}
				continue;
			}

			int idx_check;
			int elem = -1;
			int ff_num = -1;
			string ats;
			double xpos = -9999.00;
			double ypos = -9999.00;
			double zpos = -9999.00;
			
			if( str_arr.size() == 3 )
			{
				if( IsFloat(str_arr[0]) && IsFloat(str_arr[1]) && IsFloat(str_arr[2]) )
				{
					iat++;
					elem = 6;
					ats = "C" + lexical_cast<string>(iat);
					xpos = lexical_cast<double>(str_arr[0]);
					ypos = lexical_cast<double>(str_arr[1]);
					zpos = lexical_cast<double>(str_arr[2]);
				}
				else
				{
					continue;
				}
			}
			else if( str_arr.size() >= 4 )
			{
				if( IsFloat(str_arr[1]) && IsFloat(str_arr[2]) && IsFloat(str_arr[3]) )
				{	
					iat++;
					xpos = lexical_cast<double>(str_arr[1]);
					ypos = lexical_cast<double>(str_arr[2]);
					zpos = lexical_cast<double>(str_arr[3]);
					if( IsInt(str_arr[0]) )  
					{
						elem = lexical_cast<int>(str_arr[0]);
						ats = HaAtom::GetStdSymbolElem(elem) + lexical_cast<string>(iat);
					}
					else
					{
						ats = str_arr[0];
						if( ats.find_first_of("0123456789") == string::npos )
						{
							ats = ats + lexical_cast<string>(iat);
						}
						elem = HaAtom::GetElemNoFromChar(ats[0]);
					}
					if( elem < 1 ) elem = 6;
				}
				else if( IsFloat(str_arr[0]) && IsFloat(str_arr[1]) && IsFloat(str_arr[2]) )
				{
					iat++;
					xpos = lexical_cast<double>(str_arr[0]);
					ypos = lexical_cast<double>(str_arr[1]);
					zpos = lexical_cast<double>(str_arr[2]);
					if( IsInt(str_arr[3]) ) 
					{
						elem = lexical_cast<int>(str_arr[3]);
						ats = HaAtom::GetStdSymbolElem(elem) + lexical_cast<string>(iat);
					}
					else
					{
						ats = str_arr[3];
						if( ats.find_first_of("0123456789") == string::npos )
						{
							ats = ats + lexical_cast<string>(iat);
						}
						elem = HaAtom::GetElemNoFromChar(ats[0]);
					}
					if( elem < 1 ) elem = 6;
				}
				else if( str_arr.size() >= 5 && !IsFloat(str_arr[0]) && !IsFloat(str_arr[1]) )
				{
					if( IsFloat( str_arr[2]) && IsFloat(str_arr[3]) && IsFloat(str_arr[4]) ) 
					{
						iat++;
						xpos = lexical_cast<double>(str_arr[2]);
						ypos = lexical_cast<double>(str_arr[3]);
						zpos = lexical_cast<double>(str_arr[4]);
						if( IsInt( str_arr[0] ) )
						{
							idx_chk = lexical_cast<int>( str_arr[0] );
							if( IsInt(str_arr[1]) ) 
							{
								elem = lexical_cast<int>( str_arr[1] );
								ats = HaAtom::GetStdSymbolElem(elem) + lexical_cast<string>(iat);
							}
						}
						if( str_arr.size() > 5 && IsInt(str_arr[5]) )
						{
							ff_num = lexical_cast<int>( str_arr[5] );
						}
					}
					else if( str_arr.size() >= 6 && !IsFloat(str_arr[2]))
					{
						if( IsFloat( str_arr[3]) && IsFloat(str_arr[4]) && IsFloat(str_arr[5]) )
						{
							iat++;
							xpos = lexical_cast<double>(str_arr[3]);
							ypos = lexical_cast<double>(str_arr[4]);
							zpos = lexical_cast<double>(str_arr[5]);
							if( IsInt( str_arr[0] ) )
							{
								idx_chk = lexical_cast<int>( str_arr[0] );
								if( IsInt(str_arr[1]) )
								{
									elem = lexical_cast<int>( str_arr[1] );
									ats = HaAtom::GetStdSymbolElem(elem) + lexical_cast<string>(iat);
								}
							}
						}
					}
				}
			}

			if( elem < 0 ) continue;
			nat_res++;

			if(nat_res >= 100) 
			{
				res_idx++;
				pres_cur = pch_cur->AddResidue(res_idx);
				pres_cur->SetName( "MOL" );
				nat_res = 1;
			}

			aptr = pres_cur->AddNewAtom();
			atoms_new.push_back(aptr);

			aptr->SetElemNo(elem);
			aptr->SetName(ats.c_str());
			aptr->SetX_Ang( xpos );
			aptr->SetY_Ang( ypos );
			aptr->SetZ_Ang( zpos );
			
			if( ff_num > 0 )
			{
				aptr->SetFFSymbol( lexical_cast<string>(ff_num) );
			}
				
	//		int idx_b;
	//		for(;;)
	//		{
	//			is_line >> idx_b;
	//			if( is_line.fail() ) break;
	//			if( bonds.empty() ) bonds.resize(na);
	//			bonds[i].push_back(idx_b);
	//		}
		}
	}
	catch( std::exception ex )
	{
		PrintLog(" Error in MolSet::LoadXYZFile() \n");
		PrintLog( ex.what() );
	}

	//int j;
	//if( !bonds.empty() )
	//{
	//	for(i = 0; i < na; i++)
	//	{
	//		for(j = 0; j < bonds[i].size(); j++)
	//		{
	//			int idx2 = bonds[i][j];
	//			if( idx2 > na || idx2 < 1 ) continue;
	//			HaAtom* aptr1 = atoms_new[i];
	//			HaAtom* aptr2 = atoms_new[idx2-1];
	//			if( !aptr1->IsBonded(*aptr2) )
	//			{
	//				AddBond(aptr1,aptr2);
	//			}
	//		}
	//	}
	//}
	
//	if( p_opt->ToCalcBonds() && bonds.empty() )
//	{
		p_mol_editor->CreateCovBonds(pMol);
//	}
	if( p_opt->UniqueAtNames() )
	{
		pMol->SetUniqueAtomNames();
	}
	return( TRUE );
}

int MolSet::LoadHINStream( std::istream& is_arg, const AtomLoadOptions* p_opt_arg )
{
	std::auto_ptr<AtomLoadOptions> p_opt_auto( p_opt_arg == NULL ? (AtomLoadOptions*) p_load_opt_default->clone() : (AtomLoadOptions*) p_opt_arg->clone() );
	AtomLoadOptions* p_opt = p_opt_auto.get();
	
	MolEditor* p_mol_editor = this->GetMolEditor(true);
 
	if( !is_arg.good() )
	{   
        PrintLog("Error in MolSet::LoadXYZStream() \n");
        PrintLog("Can not read from input stream \n");
        return( FALSE );
	}

	using harlem::IsFloat;
	using harlem::IsInt;
	using boost::lexical_cast;
	using std::string;

    int na;
    HaAtom* aptr;
    int i;
	std::string line;

	HaMolecule* pMol = NULL;
	
	HaChain*   pch_cur  = NULL;
	HaResidue* pres_cur = NULL;

	std::map<HaAtom*,int> atom_idx_map;
	std::map<int,HaAtom*> idx_atom_map;  // indexes of atoms in the molecule

	std::string bond_types_str = "sdtav";

//	std::map<HaAtom*,std::vector<std::string>> bonds_map;

//	AtomGroup atoms_new;
//	std::vector<std::vector<int> > bonds; 
	std::string msg;

	try
	{
		int nat_res = 0;
		int res_idx = 1;
		int iat = 0;
		int idx_line = -1;
		int idx_chk  = -1;
		int na_max = 10000000;
		bool read_mol_prop = false;
		bool read_mset_prop_1 = true;
		bool read_mset_prop_2 = false;
		HaAtom* pat = NULL;

		std::vector<string> str_arr;

		std::string mol_name = "MOL";
		while(1)
		{   
			std::getline(is_arg,line);
			if( is_arg.fail() ) break;
			idx_line++;

			boost::trim(line);
			boost::split(str_arr, line, boost::is_any_of("\t "),boost::token_compress_on);
			
			if( str_arr.size() == 0 ) continue;
			
			if( line[0] == ';' ) // process comment line
			{
				std::string cmnt = line.substr(1);
				if (read_mset_prop_1)
				{
					this->comments1.push_back(cmnt);
					continue;
				}
				if (read_mset_prop_2)
				{
					this->comments2.push_back(cmnt);
					continue;
				}


				if (read_mol_prop && pMol)
				{
					pMol->comments.push_back(cmnt);
					continue;
				}
				if (pat)
				{
					pat->comments.push_back(cmnt);
				}
				continue;
			}

			if (read_mol_prop && pMol)
			{
				if (str_arr[0] == "charge")
				{
					if (str_arr.size() > 1)
					{
						try
						{
							pMol->charge = std::stoi(str_arr[1]);
						}
						catch (const std::out_of_range& ex)
						{

						}
					}
				}
			}

			if( str_arr[0] == "mol"  ) 
			{
				pMol = AddNewMolecule();
				
				if( str_arr.size() > 2 )  mol_name = str_arr[2];
				else if( str_arr.size() > 1) mol_name = str_arr[1];

				if( mol_name[0] == '\"') mol_name.erase(0,1);
				if( mol_name.size() > 1 && mol_name[ mol_name.size() - 1] == '\"') mol_name.erase(mol_name.size() - 1,1);

				if( str_arr.size() > 2) pMol->SetObjName( mol_name.c_str());

				idx_atom_map.clear();
				read_mol_prop = true;
				read_mset_prop_1 = false;
				read_mset_prop_2 = false;

				continue;
			}
			if( str_arr[0] == "endmol"  ) 
			{
				pMol = NULL;
				if( pres_cur != NULL && pres_cur->IsAmino() )  pres_cur->SetNameModifier("CT");
				pch_cur  = NULL;
				pres_cur = NULL;
				pat = NULL;
				read_mset_prop_2 = true;
				continue;
			}

			if( str_arr[0] == "res"  ) 
			{
				if( pMol == NULL ) throw std::runtime_error("No starting mol line");

				char chain_id_res = ' ';
				if( str_arr.size() > 5 )  chain_id_res = str_arr[5][0];

				bool chain_started = false;
				if( pch_cur == NULL )  
				{
					pch_cur = pMol->AddChain(chain_id_res);
					chain_started = true;
				}

				if( pch_cur->ident != chain_id_res )
				{
					if( pres_cur != NULL && pres_cur->IsAmino() )  pres_cur->SetNameModifier("CT");
					pch_cur = pMol->AddChain( chain_id_res );
					chain_started = true;
				}

				int res_ser_no = 1;
				if( str_arr.size() > 3 && harlem::IsInt(str_arr[3]) ) 
				{
					res_ser_no = boost::lexical_cast<int>(str_arr[3]);
				}
				else if( str_arr.size() > 1 && harlem::IsInt(str_arr[1]) )
				{
					res_ser_no = boost::lexical_cast<int>(str_arr[1]);
				}
				else
				{
					throw std::runtime_error("No residue number in res line");
				}
				pres_cur = pch_cur->AddResidue( res_ser_no );
				
				if( str_arr.size() > 2 ) pres_cur->SetName( str_arr[2] );
				if( pres_cur->IsAmino() && chain_started ) pres_cur->SetNameModifier("NT");

				pat = NULL;
				read_mol_prop = false;

				continue;
			}			

			if( str_arr[0] == "endres"  ) 
			{
				pat = NULL;
				continue;
			}

			if( str_arr[0] == "atom" ) 
			{
				if (str_arr.size() < 10)
				{
					PrintLog("atom line has too few fields\n%s\n", line.c_str());
					continue;
				}

				if( pch_cur == NULL )  pch_cur = pMol->AddChain(' ');
				if (pres_cur == NULL)
				{
					pres_cur = pch_cur->AddResidue(1);
					std::string res_name = mol_name;
					if (harlem::IsInt(mol_name)) res_name = "RES";
					pres_cur->SetName(res_name);
				}
				
				if( !harlem::IsInt(str_arr[1])) throw std::runtime_error( "atom index is not integer");
				int idx_at = boost::lexical_cast<int>(str_arr[1]);
				if (idx_atom_map.count(idx_at))
				{
					throw std::runtime_error(" not unique atom index in the molecule ");
					PrintLog(" not unique atom index in the molecule \n%s\n", line.c_str());
					continue;
				}
				pat = pres_cur->AddNewAtom(); 

				atom_idx_map[pat] = idx_at;
				idx_atom_map[idx_at] = pat;

				read_mol_prop = false;
				
				if (str_arr[3] != "-")
				{
					int elem_std = HaAtom::GetElemNoFromName(str_arr[3], pres_cur);
					pat->SetElemNo(elem_std);
				}

				std::string at_name = str_arr[2];

				if( isdigit(at_name[0]) ) 
				{
					at_name.erase(0,1);
					at_name.push_back( str_arr[2][0]);
				}

				if (at_name == "-")
				{
					at_name = pat->GetStdSymbol() + harlem::ToString(idx_at);
				}
				pat->SetName( at_name );
				
				if( str_arr[4] != "-")
				{
					pat->SetFFSymbol( str_arr[4] );
				}

				if( harlem::IsInt(str_arr[6]) || harlem::IsFloat(str_arr[6]) )
				{
					double at_ch = boost::lexical_cast<double>(str_arr[6]); 
					pat->SetCharge( at_ch );
				}

				double x = std::stod(str_arr[7]);
				double y = std::stod(str_arr[8]);
				double z = std::stod(str_arr[9]);

//				if( !harlem::IsFloat(str_arr[7]) ) throw std::runtime_error( "atom X coordinate is invalid");
//				if( !harlem::IsFloat(str_arr[8]) ) throw std::runtime_error( "atom Y coordinate is invalid");
//				if( !harlem::IsFloat(str_arr[9]) ) throw std::runtime_error( "atom Z coordinate is invalid");

				pat->SetX_Ang( x );
				pat->SetY_Ang( y );
				pat->SetZ_Ang( z );

				if (str_arr.size() < 11) continue;

				if (!harlem::IsInt(str_arr[10]))
				{
					PrintLog("Warning in MolSet::LoadHINStream() in line \n%s\n", line.c_str());
					PrintLog("Number of bonds field in atom record is not integer");
					continue;
				}

				int nb = boost::lexical_cast<int>(str_arr[10]);
				if ( str_arr.size() < (11 + 2 * nb))
				{
					PrintLog("Warning in MolSet::LoadHINStream()\n");
					PrintLog("The Number atom bond descriptors do not match the number of bonds\n %s \n - adjusting expected number of bonds \n", line.c_str());
					nb = (str_arr.size() - 11) / 2;
					if (nb < 0) nb = 0;
				}
				int i;
				for( i = 0; i < nb; i++)
				{
					if( !harlem::IsInt( str_arr[11 + i*2] ) )
					{
						PrintLog("Invalid atom index of bonded atom in %s \n",line.c_str());
						break;
					}
					int iat_bnd_idx = boost::lexical_cast<int>(str_arr[11 + i*2]);
					if( idx_atom_map.count( iat_bnd_idx ) == 0 ) continue;
					HaAtom* bnd_atom = idx_atom_map[ iat_bnd_idx ];
					HaBond* pbnd = AddBond( pat, bnd_atom );

					std::string btype_str = str_arr[11 + i*2 + 1];
					
					if( btype_str.size() != 1 || bond_types_str.find( str_arr[11 + i*2 + 1][0] ) == std::string::npos )
					{
						PrintLog("Invalid bond type of bonded atom in %s \n",line.c_str());
						break;
					}
					if( btype_str[0] == 'd' ) 
					{
						pbnd->SetDouble();
					}
					if( btype_str[0] == 'a' ) 
					{
						pbnd->SetAromatic();
					}
					if( btype_str[0] == 'v' ) 
					{
						pbnd->SetVirtual();
					}
				}
			}	// atom line processing
		} // while(1) 
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in MolSet::LoadHINStream()  Reading line: \n" );
		PrintLog("%s \n",line.c_str());
		PrintLog("%s\n",ex.what());
		return FALSE;
	}

	if( p_opt->ToCalcBonds() )
	{
		p_mol_editor->CreateCovBonds(pMol);
	}
	if( p_opt->UniqueAtNames() )
	{
		pMol->SetUniqueAtomNames();
	}
	return( TRUE );
}


int MolSet::LoadXMLStream (std::istream& is, const AtomLoadOptions* popt_arg )
{ 
	using namespace rapidxml;

	char* buffer = NULL;
	try
	{
		clock_t init_time = clock();

		std::auto_ptr<AtomLoadOptions> popt_auto( popt_arg == NULL ? (AtomLoadOptions*) p_load_opt_default->clone() : (AtomLoadOptions*) popt_arg->clone() );
		AtomLoadOptions* popt = popt_auto.get();

		std::streampos length;
		is.seekg (0, std::ios::end);
		length = is.tellg();
//		PrintLog("MolSet::LoadXMLStream() length = %d \n", (int)length);
		is.seekg (0, std::ios::beg);
		buffer = new char[(size_t)length+1];
		is.read(buffer,length);
//		std::streampos nread = is.gcount();
//		PrintLog("MolSet::LoadXMLStream() nread = %d \n", (int) nread);
		buffer[(size_t)length] = 0;
		  
//		clock_t time_1 = clock();
//		PrintLog(" LoadXMLStream() xml buffer load time = %12.6f sec \n", (time_1 - init_time)/((double)CLOCKS_PER_SEC));

		xml_document<> doc;
		doc.parse<0>(&buffer[0]);  
		xml_node<>* node0 = doc.first_node();

		if( node0 == NULL ) throw std::runtime_error(" No Head node in the XML stream ");
		std::string tag = node0->name();
		if( !boost::iequals(tag, "HARLEM_DATA") ) throw std::runtime_error(" Head node in the XML stream is not HARLEM_DATA ");

		xml_node<>* node1 = node0->first_node();
		if( node1 == NULL ) throw std::runtime_error(" HARLEM_DATA does not have any XML nodes");
		for( ; node1; node1 = node1->next_sibling() )
		{
			tag = node1->name(); 

			//		clock_t time_2 = clock();
			//		PrintLog(" LoadXMLStream() xml parsing time = %12.6f sec \n", (time_2 - time_1)/((double)CLOCKS_PER_SEC));

			if( boost::iequals(tag, "molset") )
			{
				this->LoadXMLNode(node1, popt );
			}
			if( boost::iequals(tag, "crd_snap") )
			{
				CrdSnapshot* psnap = new CrdSnapshot(this);
				harlem::HashMap opt;
				opt.set_a("MSET_PTR",this);
				int ires = psnap->LoadXMLNode(node1,&opt); 
				if( ires )
				{
					crd_snapshots.push_back( psnap );
				}
				else
				{
					delete psnap;
				}
			}
		}
		clock_t time_3 = clock();
		//		PrintLog(" LoadXMLStream() xml processing time = %12.6f sec \n", (time_3 - time_2)/((double)CLOCKS_PER_SEC));
	}
	catch( const std::exception& ex )
	{
		PrintLog("Error in MolSet::LoadXMLStream() \n" );
		PrintLog("%s\n",ex.what());
		if( buffer) delete buffer;
		return FALSE;
	}
	if( buffer ) delete[] buffer;
	return TRUE;
}
 
 
int MolSet::LoadMol2File(const char* fname )
{
	MolEditor* p_mol_editor = this->GetMolEditor(true);
    std::ifstream is(fname);

	if( is.fail() )
	{   
        PrintLog("\nError: File  %s not found!\n",fname );
        return( FALSE );
	}

    double xpos, ypos, zpos;
    long features, sets, serno;
    long atoms, bonds, structs;
    long srcatm, dstatm;
	std::string str,line;
 
	HaMolecule* pMol    = AddNewMolecule();
	HaChain*    pch_cur  = NULL;
	HaResidue*  pres_cur = NULL; 

	std::string molname = harlem::GetPrefixFromFullName(fname);
	pMol->SetObjName(molname.c_str());

    char name[8];
    char type[8];
 
    HaAtom  *ptr;
    int i;
 	
	IntPtrMap id_at_map;

    for(;;)
    {   
		std::getline(is,line);
		if( is.fail() ) break;

		try
		{
			if( !line[0] || line[0] == '#' )
				continue;

			if( !strncmp("@<TRIPOS>MOLECULE",line.c_str(),17) || !strncmp("@MOLECULE",line.c_str(),9) )
			{   		
				std::getline(is,line);  /* Molecule Name */
				boost::trim(line);

				pMol->SetObjName(line.c_str());

				std::getline(is,line);
				std::istringstream iss(line);
				atoms = bonds = structs = features = sets = 0;
				iss >> atoms;
				iss >> bonds;
				iss >> structs;
				iss >> features;
				iss >> sets;

				std::getline(is,line);  /* Molecule Type  */
				std::getline(is,line);  /* Charge Type    */

			} 
			else if( !strncmp("@<TRIPOS>ATOM",line.c_str(),13) || !strncmp("@ATOM",line.c_str(),5) )
			{   
				if( !atoms ) continue;

				pMol->AddChainAndResidue();
				for( i=0; i<atoms; i++ )
				{    
					std::getline(is,line);
					ptr = pres_cur->AddNewAtom();

					std::istringstream iss(line);
					iss >> serno;
					iss >> name;
					iss >> xpos;
					iss >> ypos;
					iss >> zpos;
					iss >> type;

					ptr->SetName( type );
					ptr->SetElemNo( HaAtom::GetElemNoFromName( ptr->GetName(), pres_cur ) );
					id_at_map[serno]= (void*)ptr;
					/* ptr->serno = i; */

					ptr->SetX_Ang( xpos );
					ptr->SetY_Ang( ypos );
					ptr->SetZ_Ang( zpos );
				}
			} 
			else if( !strncmp("@<TRIPOS>BOND",line.c_str(),13) || !strncmp("@BOND",line.c_str(),5) )
			{
				for( i=0; i<bonds; i++ )
				{   
					std::getline(is,line);
					std::istringstream iss(line);
					iss >> serno;
					iss >> srcatm;
					iss >> dstatm;
					iss >> type;

					HaBond* pbnd = AddBond((HaAtom*)id_at_map[srcatm],(HaAtom*)id_at_map[dstatm]);

					if( !strncmp(type,"ar",2) ) pbnd->SetAromatic();  /* AROMATIC */          
					else if( *type == '2' ) pbnd->SetDouble();        /* DOUBLE */
					else if( *type == '3' ) pbnd->SetTriple();       /* TRIPLE */
				}
			}
		}
		catch( const std::exception& ex )
		{
			PrintLog(" Error in MolSet::LoadMol2File()  line:\n");
			PrintLog(" %s\n", line.c_str());
			PrintLog(" %s\n", ex.what());
		}
	}

	if( p_load_opt_default->ToCalcBonds() )
	{
		p_mol_editor->CreateCovBonds(pMol);
	}
	if( p_load_opt_default->UniqueAtNames() )
	{
		pMol->SetUniqueAtomNames();
	}

	return( True );
}

  
/*=================================*/
/* Molecule File Format Generation */
/*=================================*/
 
int MolSet::SaveXYZRadFile(const char* filename )
{
	FILE* fout = fopen(filename,"w");
	if(fout == NULL) return FALSE;
	AtomIteratorMolSet aitr(this);
	HaAtom* aptr;
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if( !p_save_opt_default->save_selected || aptr->Selected())
		{
			fprintf(fout," %12.6f %12.6f %12.6f %12.6f \n", aptr->GetX_Ang(),aptr->GetY_Ang(),aptr->GetZ_Ang(),
					aptr->radius);
		}
	}
	fclose(fout);
    return( True );
}

int MolSet::SaveDimerXYZFile(const char* prefix )
{
	if (!this->IsDimer())
	{
		PrintLog("Molecular Set is not a dimer \n Dimer XYZ files are not saved \n");
		return FALSE;
	}

	std::vector<std::string> suffix = {"_A","_B"};
	std::vector<double> chm(2, 0.0);
	std::vector<int> nam(2, 0);

	double ch=0.0;
	int i;

	std::vector< std::vector<HaAtom*> > submol(2);
	submol[0] = this->GetAtomsSubMol(0);
	submol[1] = this->GetAtomsSubMol(1);

	int nat = this->GetNAtoms();

	for (i = 0; i < 2; i++)
	{
		nam[i] = submol[i].size();
		for (HaAtom* aptr : submol[i])
		{
			chm[i] += aptr->GetCharge();
		}
	}
	if (nat != nam[0] + nam[1])
	{
		PrintLog("Error in %s\n,Number of atoms in the dimer %d  is not equal to the sum of atom numbers in submolecules %d %d\n",
			__func__, nat, nam[0], nam[1]);
		return FALSE;
	}

	std::string fname    = (std::string)prefix + ".xyz";
	std::string fname_m1 = (std::string)prefix + "_A";
	std::string fname_m2 = (std::string)prefix + "_B";

	ofstream os(fname);
	ofstream os1(fname_m1);
	ofstream os2(fname_m2);

	if (os.fail() || os1.fail() || os2.fail())
	{
		PrintLog("Error in %s \n Error to open xyz files for writing \n", __func__);
		return FALSE;
	}

	std::string mol_name = prefix;
	
	size_t ip1 = mol_name.find("/");
	if (ip1 != std::string::npos) mol_name = mol_name.substr(ip1+1);
	size_t ip2 = mol_name.find("\\");
	if (ip2 != std::string::npos) mol_name = mol_name.substr(ip2+1);

	os << boost::format("%3d\n") % nat;
	os << boost::format("%s CH1=%d CH2=%d \n") % mol_name % (int)chm[0] % (int)chm[1];

	os1 << boost::format("%3d\n") % nam[0];
	os1 << boost::format("%s CH1=%d CH2=%d \n") % mol_name % (int)chm[0] % (int)chm[1];

	os2 << boost::format("%3d\n") % nam[1];
	os2 << boost::format("%s CH1=%d CH2=%d \n") % mol_name % (int)chm[0] % (int)chm[1];

	i = 1;
	for (HaAtom* aptr : submol[0])
	{
		std::string selem = aptr->GetStdSymbol();
		std::string ss = selem + "1" + std::to_string(i);

		os1 << boost::format("%8s %16.6f %16.6f %16.6f\n") % ss % aptr->GetX() % aptr->GetY() % aptr->GetZ();
		os <<  boost::format("%8s %16.6f %16.6f %16.6f\n") % ss % aptr->GetX() % aptr->GetY() % aptr->GetZ();
		i++;
	}

	i = 1;
	for (HaAtom* aptr : submol[1])
	{
		std::string selem = aptr->GetStdSymbol();
		std::string ss = selem + "2" + std::to_string(i);

		os2 << boost::format("%8s %16.6f %16.6f %16.6f\n") % ss % aptr->GetX() % aptr->GetY() % aptr->GetZ();
		os <<  boost::format("%8s %16.6f %16.6f %16.6f\n") % ss % aptr->GetX() % aptr->GetY() % aptr->GetZ();
		i++;
	}
	
	return TRUE; 
}

static int set_string_from_istream(std::string& str, istream& is);

int MolSet::LoadHarlemFile (const char* fname, const AtomLoadOptions* p_opt )
{
//	PrintLog(" MolSet::LoadHarlemFile() pt 1 \n");

//	clock_t init_time = clock();

	std::ifstream is(fname);
	if(is.fail()) 
	{
		PrintLog(" Error in MolSet::LoadHarlemFile() \n");
		PrintLog(" Error to open file %s\n",fname);
		return FALSE;
	}
	
	int ires;

	std::string line;
	std::getline(is,line);
	if( line.find("<?xml") == std::string::npos ) 
	{
//		PrintLog(" MolSet::LoadHarlemFile():  Attempt to open as OLD HARLEM file \n");
		is.close();
		FILE* fp = fopen(fname,"r");
		if( fp == NULL ) return FALSE;
		ires =  LoadOldHarlemFile(fp);
	}
	else
	{
		is.close();
		is.open(fname,ios::binary);
//	    PrintLog(" MolSet::LoadHarlemFile() load file in new XML format \n");
		ires = LoadXMLStream (is, p_opt );
	}
	
//	clock_t fin_time = clock();
//	PrintLog(" MolSet::LoadHarlemFile time = %12.6f sec \n", (fin_time - init_time)/((double)CLOCKS_PER_SEC));

	return ires;
}

int MolSet::LoadOldHarlemFile(FILE* fp )
{
	HaMolecule* pMol     = NULL; 
	HaChain*    pch_cur  = NULL;
	HaResidue*  pres_cur = NULL;
	HaAtom* aptr = NULL;

	enum READ_MODE {FIND_SECTION, READ_ATOMS, READ_BONDS, READ_UNIT_CELL, READ_GROUPS, READ_ATOM_LISTS, 
		            READ_CHARGE_MAPS, READ_MOLSET_PAR, READ_MOL_PAR, READ_QCHEM_MOD, READ_ET_COUPL_MOD, 
					READ_MOLMECH_MOD} read_mode;

	enum ET_READ_SUBMODE { ET_FIND_SECTION, READ_DONOR, READ_ACCEPTOR} et_read_submode;
	enum MMECH_READ_SUBMODE { MM_FIND_SECTION, READ_IMPR_ANGLES, READ_SPECIAL_VAL_BONDS, READ_SPECIAL_VAL_ANGLES } mm_read_submode;
	
	std::string atname_old="";
	std::string res_name_old="";
	int res_serno_old=-1;

	int ires;

	char* token;
	char chain_id_old;

	HaQCMod* ptr_qc_mod = NULL;
	ETCouplMod* ptr_et_coupl_mod = NULL;
	HaMolMechMod* ptr_mmech_mod = NULL;
	AtomDoubleMap* charge_map = NULL;

	AtomGroup* pdon = NULL;
	AtomGroup* pacc = NULL;
	ChemGroup* gptr = NULL;
	AtomGroup*  atlist = NULL;
	int nat_to_read = 0;

	read_mode = FIND_SECTION;

	map<int,HaAtom*, less<int> > id_atom_map;
	map<int,HaAtom*, less<int> >::iterator mitr1,mitr2;
	// Load Molecule Parameters:
		for(;;)
		{
			
			char* result=fgets(buf,MAXBUF,fp);
			if(!result || !strncmp(buf,"#END MOL_SET",12))
				break;
			
			if(!strncmp(buf,"#MOLECULE",9))
			{
				pMol=AddNewMolecule();
				read_mode= READ_MOL_PAR;
				chain_id_old='$'; // To reset old chain name
				id_atom_map.clear();
				continue;
			}
			else if(!strncmp(buf,"#ATOMS",6))
			{
				read_mode = READ_ATOMS;
				continue;
			}
			else if(!strncmp(buf,"#BONDS",6))
			{
				read_mode = READ_BONDS;
				continue;
			}
			else if(!strncmp(buf,"#UNIT CELL",10))
			{
				read_mode = READ_UNIT_CELL;
				continue;
			}
			else if(!strncmp(buf,"#GROUPS",7))
			{
				read_mode = READ_GROUPS;
				continue;
			}
			else if(!strncmp(buf,"#ATOM LISTS",11))
			{
				read_mode = READ_ATOM_LISTS;
				continue;
			}
			else if(!strncmp(buf,"#CHARGE MAPS",12))
			{
				read_mode = READ_CHARGE_MAPS;
				continue;
			}
			else if(!strncmp(buf,"#QCHEM MODULE",13))
			{
				read_mode = READ_QCHEM_MOD;
				ptr_qc_mod = GetQCMod(true);
				continue;
			}
			else if(!strncmp(buf,"#ET COUPLING MODULE",19))
			{
				read_mode = READ_ET_COUPL_MOD;
				et_read_submode = ET_FIND_SECTION;
				ptr_et_coupl_mod = GetETCouplMod(true);
				continue;
			}
			else if(!strncmp(buf,"#MOLECULAR MECHANICS MODULE",27))
			{
				read_mode = READ_MOLMECH_MOD;
				mm_read_submode = MM_FIND_SECTION;
				ptr_mmech_mod = GetMolMechMod(true);
				continue;
			}
			
			switch (read_mode)
			{
			case(READ_MOL_PAR):
				
				if(!strncmp(buf,"MOLNAME=",8))
				{
					token= strtok(buf+8," \n\t");
					if(token)
					{
						pMol->SetObjName(token);
//						PrintLog("MolSet::LoadHarlemFile()  Molecule Name : %s \n", token);
					}
					else
					{
						PrintLog("MolSet::LoadHarlemFile() : \n");
						PrintLog("Error Reading Molecule Name \n");
					}
					break;
				}
				break;
				
				// Load Atoms:
			case(READ_ATOMS):
				if(!strncmp(buf,"#END",4))
				{
					read_mode = FIND_SECTION;
					break;
				}
				else
				{
					if(pMol== NULL) break;
					istrstream is(buf);
					int atid;
					int elno;
					std::string atname;
					std::string chain_id_str;
					char chain_id;
					std::string res_name;
					std::string res_name_mod;
					int semicol_pos;
					int res_serno;
					
					double xx, yy, zz;
					double at_charge = 0.0;
					double at_mass = 0.0;
					is >> atid;
					is >> elno;
					set_string_from_istream(atname,is);
					is >> xx;
					is >> yy;
					is >> zz;
					set_string_from_istream(chain_id_str,is);
					chain_id=*(chain_id_str.c_str());
					set_string_from_istream(res_name,is);
					semicol_pos = res_name.find_first_of('#');
					if( semicol_pos != -1 )
					{
						res_name_mod = res_name.substr( semicol_pos + 1);
						res_name = res_name.substr(0, semicol_pos);
					}
					is >> res_serno;
					is >> at_charge;
					if(is.fail())
						at_charge =0.0;

					std::string hybrid_str;				
					ires = set_string_from_istream(hybrid_str,is);
					if(ires == FALSE)
						hybrid_str = "";

					std::string hb_status_str;
				    ires = set_string_from_istream(hb_status_str,is);
					if(ires == FALSE)
						hb_status_str = "";

					std::string ff_symbol_str;
				    ires = set_string_from_istream(ff_symbol_str,is);
					if(ires == FALSE)
						ff_symbol_str = "";

					is >> at_mass; 
					if( is.fail() )
						at_mass = 1.0;
					
					if(chain_id != chain_id_old)
					{
						pch_cur = pMol->AddChain(chain_id);
						chain_id_old=chain_id;
                        res_name_old = "";
                        res_serno_old = -11111;
					}
					if(res_name != res_name_old || res_serno != res_serno_old)
					{
						res_name_old=res_name;
						res_serno_old=res_serno;
						pres_cur = pch_cur->AddResidue(res_serno);
						pres_cur->SetName(res_name.c_str());
						if( !res_name_mod.empty()) pres_cur->SetNameModifier(res_name_mod.c_str());
					}
					aptr = pres_cur->AddNewAtom();
					id_atom_map[atid] = aptr;	

					char atn[5];
					atn[4] =0;
					if(atname.size() < 4) // start atom name from the second atom 
					{
						atn[0] = ' ';
						strncpy(&atn[1],atname.c_str(),3);
					}
					else
					{
						strncpy(atn,atname.c_str(),4);
					}
					aptr->SetElemNo(elno);
					aptr->SetName(atn);
					aptr->SetX_Ang(xx);
					aptr->SetY_Ang(yy);
					aptr->SetZ_Ang(zz);
					aptr->SetCharge(at_charge);
					aptr->SetHybrid(hybrid_str);
					aptr->SetHBStatus(hb_status_str.c_str());
					aptr->SetFFSymbol(ff_symbol_str);
					aptr->SetMass(at_mass);
				}
				break;
				
				
				// Load Bonds
			case(READ_BONDS):
				if(!pMol || !strncmp(buf,"#END",4))
				{
					read_mode = FIND_SECTION;
					break;
				}
				else
				{
					istrstream is(buf);
					int src_id, dst_id;
					std::string btype_str;
					is >> src_id;
					is >> dst_id;
					is >> btype_str;
					
					mitr1 = id_atom_map.find(src_id);
					mitr2 = id_atom_map.find(dst_id);
					if(mitr1 != id_atom_map.end() && mitr2 != id_atom_map.end() )
					{
						HaAtom* aptr1 = (*mitr1).second;
						HaAtom* aptr2 = (*mitr2).second;
						HaBond* bptr = AddBond(aptr1,aptr2);
						if( btype_str == "DOUBLE") bptr->SetDouble();
						if( btype_str == "AROMATIC") bptr->SetAromatic();
						if( btype_str == "VIRTUAL") bptr->SetVirtual();
					}
				}
				break;
				// Read Unit Cell parameters
			case(READ_UNIT_CELL):
				if(!pMol || !strncmp(buf,"#END",4))
				{
					read_mode = FIND_SECTION;
					break;
				}
				else
				{
					istrstream is(buf);
					double temp;
					is >> temp;
					double px = temp;
					is >> temp;
					double py = temp;
					is >> temp;
					double pz = temp;
					
					is >> temp;
					double p_alpha = temp * DEG_TO_RAD;
					is >> temp;
					double p_beta = temp * DEG_TO_RAD;
					is >> temp;
					double p_gamma = temp * DEG_TO_RAD;
					this->per_bc->SetBox(px,py,pz,p_alpha,p_beta,p_gamma);
				}
				break;
		        // Load Charge Maps:
			case(READ_CHARGE_MAPS):
				if(!strncmp(buf,"#END CHARGE MAPS",16))
				{
					read_mode = FIND_SECTION;
					break;
				}
				else
				{
					if(!strncmp(buf,"#BEGIN MAP",10))
					{
						std::string map_name( &buf[10]);
						boost::trim(map_name);
                        charge_map = CreateChargeMap( map_name.c_str() );
						break;
					}
					if(!strncmp(buf,"#END MAP",8))
					{
                        charge_map = NULL;
						break;
					}
					istrstream is(buf);
					std::string atom_ref;
					double ch;
					is >> atom_ref;
					is >> ch;

					HaAtom* aptr= GetAtomByRef(atom_ref.c_str());
					charge_map->SetValue(aptr,ch);
				}
				break;

				// Load Atomic Structural Groups:
			case(READ_GROUPS):
			case(READ_ATOM_LISTS):
				if(!strncmp(buf,"#END",4))
				{
					read_mode = FIND_SECTION;
                    nat_to_read = 0;
					gptr = NULL;
					atlist = NULL;
					break;
				}
				else
				{
					istrstream is(buf);
					if( nat_to_read < 1)
					{
						std::string gid;
						is >> gid;
						if( read_mode == READ_GROUPS)
						{
							gptr=AddBlankChemGroup(gid);
							atlist = gptr;
						}
						else
						{
							atlist = AddAtomGroup(gid.c_str());
						}

						if( read_mode == READ_GROUPS)
						{
							double protect;
							is >> protect;
							gptr->SetProtect(protect);
						}
						is >> nat_to_read;
						if( nat_to_read > GetNAtoms())
						{
							PrintLog("Error In LoadHarlemFile() \n");
							PrintLog("Number of atoms in the groups is too large \n");
						}
					}
						
					while( nat_to_read > 0)
					{
						int atid;
						std::string atom_ref;
						is >> atom_ref;
						boost::trim(atom_ref);

						if(atom_ref == "###")
						{
                            break;
						}
						if( is.fail())
						{
							PrintLog(" Error Reading atom group description from HARLEM File %s \n",
								     " End of line before the group ends ");
							break;  
						}

						atid = atoi(atom_ref.c_str());
						if(atid != 0 && pMol != NULL)
						{
							mitr1 = id_atom_map.find(atid);
							if(mitr1 != id_atom_map.end() )
							{
								aptr= (*mitr1).second;
								atlist->InsertAtom(aptr);
							}
						}
						else
						{
							aptr= GetAtomByRef(atom_ref.c_str());
							if(aptr != NULL)
							{
								atlist->InsertAtom(aptr);
							}
						}
						nat_to_read--;
					}
						
					if( nat_to_read == 0)
					{
                       atlist = NULL;
					   gptr = NULL;
					}
				}
				break;
				// Load QChem Module Parameters:
			case(READ_QCHEM_MOD):
				if(!strncmp(buf,"#END",4))
				{
					read_mode = FIND_SECTION;
					break;
				}
				else
				{
					//if(!strncmp(buf,"BASIS=",6))
					//{
					//	
					//	token= strtok(buf+6," \n\t");
					//	if(token)
					//	{
					//		if(ptr_qc_mod != NULL) ptr_qc_mod->InitBasis(token);
					//	}
					//	else
					//	{
					//		cout << "MolSet::LoadHarlemFile() : " << endl;
					//		cout << " Error Reading Basis Name " << endl;
					//	}
					//	break;
					//}
					//else if(!strncmp(buf,"LOC_ORB_BASIS=",14))
					//{
					//	token= strtok(buf+14," \n\t");
					//	if(token)
					//	{
					//		if(ptr_qc_mod != NULL) ptr_qc_mod->InitLocOrb(token);
					//	}
					//	else
					//	{
					//		cout << "MolSet::LoadHarlemFile() : " << endl;
					//		cout << " Error Reading Basis Name " << endl;
					//	}
					//	break;
					//}
				}
				break;
				// Load ET Coupling mode parameters:
			case(READ_ET_COUPL_MOD):
				if(!strncmp(buf,"#END",4))
				{
					read_mode = FIND_SECTION;
					et_read_submode = ET_FIND_SECTION;
					break;
				}
				else
				{
					if( et_read_submode == ET_FIND_SECTION)
					{
						if( !strncmp(buf,"DONOR",5) )
						{
							et_read_submode = READ_DONOR;
							pdon = GetAtomGroupByID("DONOR");
							if(pdon == NULL) pdon = AddAtomGroup("DONOR");
						}
						if( !strncmp(buf,"ACCEPTOR",8) )
						{
							et_read_submode = READ_ACCEPTOR;
                            pacc = GetAtomGroupByID("ACCEPTOR");
							if(pacc == NULL) pacc = AddAtomGroup("ACCEPTOR");
						}
						if( !strncmp(buf,"END_DONOR",9) )
						{
							et_read_submode = ET_FIND_SECTION;
						}
						if( !strncmp(buf,"END_ACCEPTOR",12) )
						{
							et_read_submode = ET_FIND_SECTION;
						}
					}
					else if( et_read_submode == READ_DONOR || et_read_submode == READ_ACCEPTOR)
					{
						if( !strncmp(buf,"END_DONOR",9) )
						{
							et_read_submode = ET_FIND_SECTION;
						}
						else if( !strncmp(buf,"END_ACCEPTOR",12) )
						{
							et_read_submode = ET_FIND_SECTION;
						}  
						else
						{
							std::string atom_ref(buf);
							boost::trim(atom_ref);
							HaAtom* aptr= GetAtomByRef(atom_ref.c_str());
							if(aptr != NULL)
							{
								if(et_read_submode == READ_DONOR)
									if(pdon) pdon->InsertAtom(aptr);
								if(et_read_submode == READ_ACCEPTOR)
									if(pacc) pacc->InsertAtom(aptr);
							}
						}
					}
				}
			case(READ_MOLMECH_MOD):
				if(!strncmp(buf,"#END",4))
				{
					read_mode = FIND_SECTION;
					mm_read_submode = MM_FIND_SECTION;
					break;
				}
				else
				{
					if( mm_read_submode == MM_FIND_SECTION)
					{
						if( !strncmp(buf,"IMPROPER ANGLES",15) )
						{
							mm_read_submode = READ_IMPR_ANGLES;
						}
						if( !strncmp(buf,"SPECIAL VALENCE BONDS",21) )
						{
							mm_read_submode = READ_SPECIAL_VAL_BONDS;
						}
						if( !strncmp(buf,"SPECIAL VALENCE ANGLES",22) )
						{
							mm_read_submode = READ_SPECIAL_VAL_ANGLES;
						}
					}
					else if( mm_read_submode == READ_IMPR_ANGLES )
					{
						if( !strncmp(buf,"END IMPROPER ANGLES",19) )
						{
							mm_read_submode = MM_FIND_SECTION;
						}
						else 
						{
							std::string aref1, aref2, aref3, aref4;
							istrstream is(buf);
							is >> aref1;
							is >> aref2;
							is >> aref3;
							is >> aref4;
	
							HaAtom* aptr1 = GetAtomByRef(aref1.c_str());
							HaAtom* aptr2 = GetAtomByRef(aref2.c_str());
							HaAtom* aptr3 = GetAtomByRef(aref3.c_str());
							HaAtom* aptr4 = GetAtomByRef(aref4.c_str());
							
							if( aptr1 == NULL || aptr2 == NULL || aptr3 == NULL || aptr4 == NULL )
							{
								PrintLog(" Error reading Atom References in dihedral angle %s - %s - %s - %s \n",
									      aref1.c_str(), aref2.c_str(), aref3.c_str(), aref4.c_str() );
								continue;
							}
							if( ptr_mmech_mod != NULL)
							{
								MolMechModel* p_mm_model = ptr_mmech_mod->GetMolMechModel();
								p_mm_model->AddImprDihedral(aptr1,  aptr2, aptr3, aptr4);

								if(debug_flag & file_reading_debug)
								{
									PrintLog(" Add Improper Angle between atoms %s %s %s %s \n",
										      aref1.c_str(), aref2.c_str(), aref3.c_str(), aref4.c_str() ); 
								}
							}
						}
					}
					else if( mm_read_submode == READ_SPECIAL_VAL_BONDS )
					{
						if( !strncmp(buf,"END SPECIAL VALENCE BONDS",25) )
						{
							mm_read_submode = MM_FIND_SECTION;
						}
						else 
						{
							std::string aref1, aref2;
							double r0,fc;
							istrstream is(buf);
							is >> aref1;
							is >> aref2;
							is >> r0;
							is >> fc;
							
							HaAtom* aptr1 = GetAtomByRef(aref1.c_str());
							HaAtom* aptr2 = GetAtomByRef(aref2.c_str());
							
							if( aptr1 == NULL || aptr2 == NULL )
							{
								PrintLog(" Error reading Atom References in special val bond %s - %s \n",
									      aref1.c_str(), aref2.c_str() );
								continue;
							}
							if( ptr_mmech_mod != NULL)
							{
								MolMechModel* p_mm_model = ptr_mmech_mod->GetMolMechModel();
								p_mm_model->SetMMBond(aptr1,aptr2,r0,fc,MolMechModel::SET_SPEC);

								if(debug_flag & file_reading_debug)
								{
									PrintLog(" Add Special Bond between atoms %s %s with r0= %12.6f fc= %12.6f \n",
										      aref1.c_str(), aref2.c_str(), r0,fc ); 
								}
							}
						}
					}
					else if( mm_read_submode == READ_SPECIAL_VAL_ANGLES )
					{
						if( !strncmp(buf,"END SPECIAL VALENCE ANGLES",26) )
						{
							mm_read_submode = MM_FIND_SECTION;
						}
						else 
						{
							std::string aref1, aref2, aref3;
							double a0,fc;
							istrstream is(buf);
							is >> aref1;
							is >> aref2;
							is >> aref3;
							is >> a0;
							is >> fc;
							
							HaAtom* aptr1 = GetAtomByRef(aref1.c_str());
							HaAtom* aptr2 = GetAtomByRef(aref2.c_str());
							HaAtom* aptr3 = GetAtomByRef(aref3.c_str());
							
							if( aptr1 == NULL || aptr2 == NULL || aptr3 == NULL )
							{
								PrintLog(" Error reading Atom References in special val angle %s - %s - %s \n",
									      aref1.c_str(), aref2.c_str(), aref3.c_str() );
								continue;
							}
							if( ptr_mmech_mod != NULL)
							{
								MolMechModel* p_mm_model = ptr_mmech_mod->GetMolMechModel();
								p_mm_model->SetValAngle(aptr1,aptr2,aptr3,a0,fc,MolMechModel::SET_SPEC);

								if(debug_flag & file_reading_debug)
								{
									PrintLog(" Add Special Angle between atoms %s %s %s with a0= %12.6f fc= %12.6f \n",
										      aref1.c_str(), aref2.c_str(), aref3.c_str(), a0,fc ); 
								}
							}
						}
					}



				}
		}		
    }
	
	if(pMol && pMol->GetNAtoms() == 0) 
	{
		DeleteMol(pMol);
	}

	fclose(fp);
	return True;
}

int
MolSet::LoadAmberPrepFile(const char* fname)
{
	MolEditor* p_mol_editor = this->GetMolEditor(true);
    FILE* fp = fopen(fname,"r");

	if( !fp )
	{   
        PrintLog("\nError: File '");
        PrintLog(fname);
        PrintLog("' not found!\n\n");
        return( False );
	}

	HaMolecule* pMol     = NULL; 
	HaChain*    pch_cur  = NULL;
	HaResidue*  pres_cur = NULL;
	HaAtom* aptr = NULL;
	char buf[256];

	HaMolMechMod* ptr_mm_mod = GetMolMechMod(true);

	int i;
	char* result;
	for(i= 0; i < 2; i++)
	{
		result=fgets(buf,255,fp);
		if(result == NULL)
			return TRUE;
	}

	for(;;)
	{
		for(i= 0; i < 3; i++)
		{
			result=fgets(buf,255,fp);
			if(result == NULL)
				return TRUE;
		}
		istrstream is(buf);
		std::string mol_name;
		is >> mol_name;
		pMol=AddNewMolecule();
		pch_cur =pMol->AddChain(' ');
		pMol->SetObjName(mol_name.c_str());

		bool build_bonds = false;

		pres_cur = pch_cur->AddResidue(1);
		pres_cur->SetName(mol_name.c_str());


		for(i= 0; i < 2; i++)
		{
			result=fgets(buf,255,fp);
			if(result == NULL)
				return TRUE;
		}
		map<int, Vec3D*, less<int> > num_aptr_map;

		Vec3D ref_pt1; 
		Vec3D ref_pt2; 
		Vec3D ref_pt3;
		
		ref_pt1.SetX(0.0);  ref_pt1.SetY(1.0); ref_pt1.SetZ(0.0); 
		ref_pt2.SetX(0.0);  ref_pt2.SetY(0.0); ref_pt2.SetZ(0.0);
		ref_pt3.SetX(1.0);  ref_pt2.SetY(0.0); ref_pt2.SetZ(0.0);
	
		num_aptr_map[1] = &ref_pt1;
		num_aptr_map[2] = &ref_pt2;
		num_aptr_map[3] = &ref_pt3;

		for(;;)
		{
			result=fgets(buf,255,fp);
			istrstream is2(buf);
			int iat1;
			is2 >> iat1;
			if(is2.fail())
				break;
			std::string atname;
			std::string ff_symbol;
			std::string branch_symbol;

			std::string numstr;

			int iat2, iat3, iat4;
			double bond_len;
			double val_angle;
			double dih_angle;
			double at_ch;
			double x_coord, y_coord, z_coord;

			is2 >> atname;
			is2 >> ff_symbol;
			is2 >> branch_symbol;
			is2 >> numstr;
			
			bool internal_coord= false;

			int dot_pos = numstr.find_first_of('.');
			if( dot_pos == -1)
			{
				internal_coord = true;
			}
			else
				build_bonds = true;

			
			if( internal_coord)
			{
				iat2 = atoi(numstr.c_str());
				is2 >> iat3 >> iat4;
				is2 >> bond_len >> val_angle >> dih_angle;
				is2 >> at_ch;
			}
			else
			{
				x_coord = atof(numstr.c_str());
				is2 >> y_coord >> z_coord;
				is2 >> at_ch;
			}

			if(iat1 > 3)
			{
				char atn[5];
				aptr= pres_cur->AddNewAtom();
				strncpy(atn,atname.c_str(),5);
				num_aptr_map[iat1]= aptr;
				aptr->SetName(atn);
				aptr->SetElemNo(HaAtom::GetElemNoFromName(atname, pres_cur) );
			}
			else
				continue;
			if(aptr == NULL)
				continue;

			int err_ref_atom = FALSE;
			
			if(internal_coord)
			{
				if(num_aptr_map.find(iat2) == num_aptr_map.end() || 
					num_aptr_map.find(iat3) == num_aptr_map.end() ||
					num_aptr_map.find(iat4) == num_aptr_map.end() )
				{
					ErrorInMod( "MolSet::LoadAmberPrepFile()",
						"Can't find reference points In the map");
					aptr->SetX(0.0); aptr->SetY(0.0); aptr->SetZ(0.0);
				}
				else
				{
					Vec3D* ptr2 = num_aptr_map[iat2]; 
					Vec3D* ptr3 = num_aptr_map[iat3];
					Vec3D* ptr4 = num_aptr_map[iat4];
					val_angle *= DEG_TO_RAD;
					dih_angle *= DEG_TO_RAD;
					Vec3D::SetAtomPos(aptr, ptr2, ptr3, ptr4, 
						bond_len, val_angle, -dih_angle);  
					if(iat2 > 3) AddBond(aptr, (HaAtom*) ptr2 );  
				}
			}
			else
			{
				aptr->SetX_Ang(x_coord);
				aptr->SetY_Ang(y_coord);
				aptr->SetZ_Ang(z_coord);
			}

			aptr->SetCharge(at_ch);
			aptr->SetFFSymbol(ff_symbol);
			aptr->SetMass(HaAtom::StdElemMass(aptr->GetElemNo()) );
		}

		
		for(;;)
		{
			result=fgets(buf,255,fp);
			if(result == NULL)
				return TRUE;
			if( strncmp(buf,"LOOP",4) == 0)
			{
				for(;;)
				{
					result=fgets(buf,255,fp);
					if(result == NULL) 
						return TRUE;
					istrstream is3(buf);
					std::string atname_1;
					std::string atname_2;
					is3 >> atname_1 >> atname_2;
					if(is3.fail())
						break;
					HaAtom* at1 = pres_cur->GetAtomByName(atname_1.c_str());
					HaAtom* at2 = pres_cur->GetAtomByName(atname_2.c_str());
					if(at1 == NULL || at2 == NULL)
						continue;
					AddBond(at1,at2);
				}
			}
			if( strncmp(buf,"IMPROPER",8) == 0)
			{
				for(;;)
				{
					result=fgets(buf,255,fp);
					if(result == NULL) 
						return TRUE;
					istrstream is3(buf);
					std::string atname_1;
					std::string atname_2;
					std::string atname_3;
					std::string atname_4;

					is3 >> atname_1 >> atname_2 >> atname_3 >> atname_4;
					if(is3.fail())
						break;
					HaAtom* at1 = pres_cur->GetAtomByName(atname_1.c_str());
					HaAtom* at2 = pres_cur->GetAtomByName(atname_2.c_str());
					HaAtom* at3 = pres_cur->GetAtomByName(atname_3.c_str());
					HaAtom* at4 = pres_cur->GetAtomByName(atname_4.c_str());

					if(at1 == NULL || at2 == NULL || at3 == NULL || at4 == NULL)
						continue;
					if(ptr_mm_mod != NULL) 
					{
						MolMechModel* p_mm_model = ptr_mm_mod->GetMolMechModel();
						p_mm_model->AddImprDihedral(at1,at2,at3,at4);
					}
				}
			}
			if( strncmp(buf,"DONE",4) == 0)
			{
				/*  Add Terminal Atoms: */
				
				if(pres_cur->IsProtein())
				{
					char atn[4];
					HaAtom* aptr_n = pres_cur->GetAtomByName("N");
					HaAtom* aptr_ca = pres_cur->GetAtomByName("CA");
					HaAtom* aptr_h =  pres_cur->GetAtomByName("H");
					HaAtom* aptr_c =  pres_cur->GetAtomByName("C");
					HaAtom* aptr_o =  pres_cur->GetAtomByName("O");
					if(stricmp_loc(pres_cur->GetName(), "PRO") == 0 )
						aptr_h = pres_cur->GetAtomByName("CG");
					
					if( aptr_n != NULL && aptr_ca != NULL && aptr_h != NULL &&
						aptr_c != NULL && aptr_o != NULL)
					{	
						pres_cur = pch_cur->AddResidue(2);
						
						pres_cur->SetName("HED");
						aptr= pres_cur->AddNewAtom();
						strncpy(atn,"HT",4);
						aptr->SetName(atn);
						aptr->SetElemNo(1);
						double bond_len =1.08;
						double val_angle = 115.0*DEG_TO_RAD;
						double dih_angle = 180.0 *DEG_TO_RAD;
						if(aptr != NULL)
						{
							Vec3D::SetAtomPos(aptr, aptr_n, aptr_ca, aptr_h, 
								bond_len, val_angle, dih_angle);  
							AddBond(aptr, aptr_n );
							if(ptr_mm_mod != NULL) 
							{
								MolMechModel* p_mm_model = ptr_mm_mod->GetMolMechModel();
								p_mm_model->AddImprDihedral(aptr, aptr_ca, aptr_n, aptr_h);
							}
						}
						pres_cur = pch_cur->AddResidue(3);
						pres_cur->SetName("TAL");
						aptr = pres_cur->AddNewAtom();
						strncpy(atn,"HT",4);
						aptr->SetElemNo(1);
						aptr->SetName(atn);
						
						val_angle = 123.5*DEG_TO_RAD;
						if(aptr != NULL)
						{
							Vec3D::SetAtomPos(aptr, aptr_c, aptr_ca, aptr_o, 
								bond_len, val_angle, dih_angle);  
							AddBond(aptr, aptr_c);
							if(ptr_mm_mod != NULL) 
							{
								MolMechModel* p_mm_model = ptr_mm_mod->GetMolMechModel();
								p_mm_model->AddImprDihedral(aptr_ca, aptr, aptr_c, aptr_o);
							}
						}
					}
				}
				break;
			}
		}
		if(build_bonds) p_mol_editor->CreateCovBonds(pMol);
	}

	fclose(fp);
	return TRUE;
}

int MolSet::LoadAmberTopFile(const char* fname)
{
	MolEditor* p_mol_editor = this->GetMolEditor(true);
    FILE* fp = fopen(fname,"r");

	if( !fp )
	{   
        PrintLog("\nError: File '");
        PrintLog(fname);
        PrintLog("' not found!\n\n");
        return( False );
	}

	HaMolecule* pMol = NULL; 
	HaChain*   pch_cur  = NULL;
	HaResidue* pres_cur = NULL;
	HaAtom* aptr = NULL;

	char buf[256];

	HaMolMechMod* ptr_mm_mod = GetMolMechMod(true);

	int i;
	char* result;
	for(i= 0; i < 2; i++)
	{
		result=fgets(buf,255,fp);
		if(result == NULL)
			return TRUE;
	}

	for(;;)
	{
		for(i= 0; i < 3; i++)
		{
			result=fgets(buf,255,fp);
			if(result == NULL)
				return TRUE;
		}
		istrstream is(buf);
		std::string mol_name;
		is >> mol_name;
		pMol=AddNewMolecule();
		pch_cur =pMol->AddChain(' ');
		pMol->SetObjName(mol_name.c_str());

		bool build_bonds = false;

		pres_cur = pch_cur->AddResidue(1);
		pres_cur->SetName(mol_name.c_str());


		for(i= 0; i < 2; i++)
		{
			result=fgets(buf,255,fp);
			if(result == NULL)
				return TRUE;
		}
		map<int, Vec3D*, less<int> > num_aptr_map;

		Vec3D ref_pt1; 
		Vec3D ref_pt2; 
		Vec3D ref_pt3;
		
		ref_pt1.SetX(0.0);  ref_pt1.SetY(1.0); ref_pt1.SetZ(0.0); 
		ref_pt2.SetX(0.0);  ref_pt2.SetY(0.0); ref_pt2.SetZ(0.0);
		ref_pt3.SetX(1.0);  ref_pt2.SetY(0.0); ref_pt2.SetZ(0.0);
	
		num_aptr_map[1] = &ref_pt1;
		num_aptr_map[2] = &ref_pt2;
		num_aptr_map[3] = &ref_pt3;

		for(;;)
		{
			result=fgets(buf,255,fp);
			istrstream is2(buf);
			int iat1;
			is2 >> iat1;
			if(is2.fail())
				break;
			std::string atname;
			std::string ff_symbol;
			std::string branch_symbol;

			std::string numstr;

			int iat2, iat3, iat4;
			double bond_len;
			double val_angle;
			double dih_angle;
			double at_ch;
			double x_coord, y_coord, z_coord;

			is2 >> atname;
			is2 >> ff_symbol;
			is2 >> branch_symbol;
			is2 >> numstr;
			
			bool internal_coord= false;

			int dot_pos = numstr.find_first_of('.');
			if( dot_pos == -1)
			{
				internal_coord = true;
			}
			else
				build_bonds = true;

			
			if( internal_coord)
			{
				iat2 = atoi(numstr.c_str());
				is2 >> iat3 >> iat4;
				is2 >> bond_len >> val_angle >> dih_angle;
				is2 >> at_ch;
			}
			else
			{
				x_coord = atof(numstr.c_str());
				is2 >> y_coord >> z_coord;
				is2 >> at_ch;
			}

			if(iat1 > 3)
			{
				char atn[5];
				aptr=pres_cur->AddNewAtom();
				strncpy(atn,atname.c_str(),5);
				num_aptr_map[iat1]= aptr;
				aptr->SetName(atn);
				aptr->SetElemNo(HaAtom::GetElemNoFromName(atname, aptr->GetHostRes()) );
			}
			else
				continue;
			if(aptr == NULL)
				continue;

			int err_ref_atom = FALSE;
			
			if(internal_coord)
			{
				if(num_aptr_map.find(iat2) == num_aptr_map.end() || 
					num_aptr_map.find(iat3) == num_aptr_map.end() ||
					num_aptr_map.find(iat4) == num_aptr_map.end() )
				{
					ErrorInMod( "MolSet::LoadAmberTopFile()",
						"Can't find reference points In the map");
					aptr->SetX(0.0); aptr->SetY(0.0); aptr->SetZ(0.0);
				}
				else
				{
					Vec3D* ptr2 = num_aptr_map[iat2]; 
					Vec3D* ptr3 = num_aptr_map[iat3];
					Vec3D* ptr4 = num_aptr_map[iat4];
					val_angle *= DEG_TO_RAD;
					dih_angle *= DEG_TO_RAD;
					Vec3D::SetAtomPos(aptr, ptr2, ptr3, ptr4, 
						bond_len, val_angle, -dih_angle);  
					if(iat2 > 3) AddBond(aptr, (HaAtom*) ptr2 );  
				}
			}
			else
			{
				aptr->SetX_Ang(x_coord);
				aptr->SetY_Ang(y_coord);
				aptr->SetZ_Ang(z_coord);
			}

			aptr->SetCharge(at_ch);
			aptr->SetFFSymbol(ff_symbol);
			aptr->SetMass(HaAtom::StdElemMass(aptr->GetElemNo()) );
		}

		
		for(;;)
		{
			result=fgets(buf,255,fp);
			if(result == NULL)
				return TRUE;
			if( strncmp(buf,"LOOP",4) == 0)
			{
				for(;;)
				{
					result=fgets(buf,255,fp);
					if(result == NULL) 
						return TRUE;
					istrstream is3(buf);
					std::string atname_1;
					std::string atname_2;
					is3 >> atname_1 >> atname_2;
					if(is3.fail())
						break;
					HaAtom* at1 = pres_cur->GetAtomByName(atname_1.c_str());
					HaAtom* at2 = pres_cur->GetAtomByName(atname_2.c_str());
					if(at1 == NULL || at2 == NULL)
						continue;
					AddBond(at1,at2);
				}
			}
			if( strncmp(buf,"IMPROPER",8) == 0)
			{
				for(;;)
				{
					result=fgets(buf,255,fp);
					if(result == NULL) 
						return TRUE;
					istrstream is3(buf);
					std::string atname_1;
					std::string atname_2;
					std::string atname_3;
					std::string atname_4;

					is3 >> atname_1 >> atname_2 >> atname_3 >> atname_4;
					if(is3.fail())
						break;
					HaAtom* at1 = pres_cur->GetAtomByName(atname_1.c_str());
					HaAtom* at2 = pres_cur->GetAtomByName(atname_2.c_str());
					HaAtom* at3 = pres_cur->GetAtomByName(atname_3.c_str());
					HaAtom* at4 = pres_cur->GetAtomByName(atname_4.c_str());

					if(at1 == NULL || at2 == NULL || at3 == NULL || at4 == NULL)
						continue;
					if(ptr_mm_mod != NULL) 
					{
						MolMechModel* p_mm_model = ptr_mm_mod->GetMolMechModel();
						p_mm_model->AddImprDihedral(at1,at2,at3,at4);
					}
				}
			}
			if( strncmp(buf,"DONE",4) == 0)
			{
				/*  Add Terminal Atoms: */
				
				if(pres_cur->IsProtein())
				{
					char atn[4];
					HaAtom* aptr_n  = pres_cur->GetAtomByName("N");
					HaAtom* aptr_ca = pres_cur->GetAtomByName("CA");
					HaAtom* aptr_h  = pres_cur->GetAtomByName("H");
					HaAtom* aptr_c  = pres_cur->GetAtomByName("C");
					HaAtom* aptr_o  = pres_cur->GetAtomByName("O");
					if(stricmp_loc(pres_cur->GetName(), "PRO") == 0 )
						aptr_h = pres_cur->GetAtomByName("CG");
					

					if( aptr_n != NULL && aptr_ca != NULL && aptr_h != NULL &&
						aptr_c != NULL && aptr_o != NULL)
					{
						
						pres_cur = pch_cur->AddResidue(2);
						
						pres_cur->SetName("HED");
						aptr= pres_cur->AddNewAtom();
						strncpy(atn,"HT",4);
						aptr->SetName(atn);
						aptr->SetElemNo(1);
						double bond_len =1.08;
						double val_angle = 115.0*DEG_TO_RAD;
						double dih_angle = 180.0*DEG_TO_RAD;
						if(aptr != NULL)
						{
							Vec3D::SetAtomPos(aptr, aptr_n, aptr_ca, aptr_h, 
								bond_len, val_angle, dih_angle);  
							AddBond(aptr, aptr_n );
							if(ptr_mm_mod != NULL) 
							{
								MolMechModel* p_mm_model = ptr_mm_mod->GetMolMechModel();
								p_mm_model->AddImprDihedral(aptr, aptr_ca, aptr_n, aptr_h);
							}
						}
						pres_cur = pch_cur->AddResidue(3);
						pres_cur->SetName("TAL");
						aptr = pres_cur->AddNewAtom();
						strncpy(atn,"HT",4);
						aptr->SetElemNo(1);
						aptr->SetName(atn);
						
						val_angle = 123.5*DEG_TO_RAD;
						if(aptr != NULL)
						{
							Vec3D::SetAtomPos(aptr, aptr_c, aptr_ca, aptr_o, 
								bond_len, val_angle, dih_angle);  
							AddBond(aptr, aptr_c );
							if(ptr_mm_mod != NULL) 
							{
								MolMechModel* p_mm_model = ptr_mm_mod->GetMolMechModel();
								p_mm_model->AddImprDihedral(aptr_ca, aptr, aptr_c, aptr_o);
							}
						}
					}
				}
				break;
			}
		}
		if(build_bonds) p_mol_editor->CreateCovBonds(pMol);
	}

	fclose(fp);
	return TRUE;
}


static int set_string_from_istream(std::string& str, istream& is)
{
	char ch;
	str.erase();
	int inside_quotes=0;

	for(;;) // find beginning of the string
	{
		is.get(ch);
		if(is.fail() ) return FALSE;
		if(isspace(ch))continue;
		if(ch =='"')
		{
			is.get(ch);
			inside_quotes=1;
			break;
		}
		break;
	}
	for(;;) // read the string in 
	{
		if(is.fail()) return TRUE;
		if(ch == '"') return TRUE;
		if(isspace(ch) && !inside_quotes) return TRUE;
		str+=ch;
		is.get(ch);
	}
}

int MolSet::LoadRWFMolecule(const char* fname)
{
	PrintLog(" This is LoadRWFMolecule \n ");

	HaMolecule* pMol=AddNewMolecule();	

	GauFile  gfile;
	gfile.set_file_name(fname);
	gfile.set_gau_iunit(1);
	gfile.set_open_mode("unknown");
    gfile.open();

	std::string molname = harlem::GetPrefixFromFullName(fname);
	pMol->SetObjName(molname.c_str());

	pMol->InitAtoms(gfile); 
	
	HaQCMod* ptr_qc_mod = GetQCMod(true); 

	ptr_qc_mod->InitBasis(gfile); // Load GauBas Class from Gaussian File 
 	
//    std::string bname("3-21G");
  
//    InitBasis(bname);
//    InitLocOrb("FULL_ATOMIC_BASIS"); 

//    InitMOs(gfile);

//    InitBasOvlp();

	gfile.close();

	return (True);
}
