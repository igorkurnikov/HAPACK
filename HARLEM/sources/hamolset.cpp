/*! \file hamolset.cpp

    Set of Molecules class in HARLEM.
 
    \author Igor Kurnikov
    \date 1999-2002

*/  

#define HAMOLSET_CPP

#if !defined(HARLEM_PYTHON_NO)
#include "Python.h"
#endif

#include  <float.h>
#include  <math.h>

#include <mpi.h>

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

#include "rapidxml.hpp"

#include "harlemapp.h"
#include "haio.h"
#include "atom_mapping.h"
#include "haatgroup.h"
#include "hamolset.h"
#include "command.h"
#include "hamolecule.h"
#include "tokens.h"
#include "abstree.h"
#include "hamolview.h"
#include "haqchem.h"
#include "etcoupl.h"
#include "hacompmod.h"
#include "hagaussian.h"
#include "hadalton.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "hamolmech.h"
#include "haintermol.h"
#include "haresdb.h"
#include "hasurface.h"
#include "electrostmod.h"
#include "protonredox.h"
#include "moleditor.h"

#include "argraph.h"
#include "argedit.h"
#include "vf2_sub_state.h"
#include "match.h"

#include "randomip.h"
#include <tinyxml.h>
#include "haobject.h"
#include "apbsmod.h"

#include "wx/log.h"
#include "hawx_add.h"
#include "molset_evt_handler.h"
#include "hamatdb.h"



AtomLoadOptions* MolSet::p_load_opt_default = new AtomLoadOptions();
AtomSaveOptions* MolSet::p_save_opt_default = new AtomSaveOptions();

MolSet* MolSet::CurMolSet = NULL;

void SetCurMolSet(MolSet* pmset) 
{ 
	if( pmset == NULL)
	{
		MolSet::CurMolSet  = NULL;
		CurMolView = NULL;
	}
	MolSet::CurMolSet = pmset;
#if !defined(HARLEM_PYTHON_NO)
	PyGILState_STATE gstate;
	gstate = PyGILState_Ensure();
	int ires = PyRun_SimpleString("import molset");
	ires = PyRun_SimpleString("mset_c = molset.GetCurMolSet\(\) ");
	PyGILState_Release(gstate);
	if(!ires) return;
	if( pmset != NULL )
	{
		gstate = PyGILState_Ensure();
		ires = PyRun_SimpleString(
			"mmod_c  = mset_c.GetMolMechMod(False)\n"
			"elmod_c = mset_c.GetElectrostMod(False)\n"
			"qcmod_c = mset_c.GetQCMod(False)\n"
		);
		PyGILState_Release(gstate);
	}
#endif
}

MolSet::MolSet()
{	
	CurMolSet=this;
	
	debug_flag = 0;

	mset_pview = NULL;
	canvas_wx = NULL;

	per_bc = new PeriodicUnitInfo();

	SSBonds_found  = false;
	HBonds_found   = false;
	to_find_backb  = true;

	p_zmat = new ZMatCrd(this);
	p_mol_editor  = new MolEditor();
	p_evt_h       = new MolSetEvtHandler(this);

	if( pApp ) pApp->AddMolSet(this);

	this->SetName("MOLSET");
}

MolSet::~MolSet()
{
	if(CurMolSet == this) CurMolSet = NULL;
	DeleteAll();
	if(p_mol_editor)  delete p_mol_editor;
	if(per_bc)  delete per_bc;
	if(p_evt_h) delete p_evt_h;
	if(p_zmat)  delete p_zmat;
	
	if( pApp ) pApp->DeleteMolSet(this);
}


void MolSet::DeleteAll()
{	
	ReleaseAllFragments();

	vector<HaCompMod*>::iterator mitr;
	for(mitr = CompModules.begin(); mitr != CompModules.end(); mitr++)
	{
		HaCompMod* comp_mod = *mitr;
		if( comp_mod != NULL ) delete (comp_mod);
	}
	CompModules.clear();

	list<Object3D*>::iterator oitr;
	for(oitr = ViewObjects.begin(); oitr != ViewObjects.end(); oitr++)
	{
		if( (*oitr) != NULL && (*oitr)->GetObjType() != OBJ3D_MOLECULE )
			delete (*oitr);
	}
	ViewObjects.clear();

	AtomIteratorMolSet aitr(this);
	HaAtom* aptr;
	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if( aptr != NULL) delete aptr;
	}

	MoleculesType::iterator mol_itr;
	for( mol_itr=HostMolecules.begin(); mol_itr != HostMolecules.end(); mol_itr++)
	{
		if((*mol_itr) != NULL) delete *mol_itr;
	}
	HostMolecules.clear();
	ClearBackbone();
	
	DeleteCrdSnapshots();
} 



int MolSet::SavePDBFile(const char* filename)
{
	double x, y, z;
    
	HaChain  *chain;
	HaResidue  *group;
	HaAtom  *aptr;
	int count;

	FILE* DataFile = fopen( filename, "w" );
	if( !DataFile )
	{
		PrintLog("\n");
		return( False );
	}

	count = 1;
	MoleculesType::iterator mol_itr;
	for( mol_itr=HostMolecules.begin(); mol_itr != HostMolecules.end(); mol_itr++)
	{
		ChainIteratorMolecule ch_itr(*mol_itr);
		for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
		{
			ResidueIteratorChain ritr_ch(chain);
			for(group = ritr_ch.GetFirstRes(); group; group = ritr_ch.GetNextRes())
			{
				std::string res_name = group->GetName();
				if( p_save_opt_default->save_amber_pdb )
				{
					res_name = MMForceField::GetAmberResName( group->GetFullName() );
				}
					
				AtomIteratorAtomGroup aitr_group(group);
				for(aptr= aitr_group.GetFirstAtom(); aptr; aptr = aitr_group.GetNextAtom())
				{
					if( !p_save_opt_default->save_selected || aptr->Selected())
					{
//						if( prev && (chain->ident!=ch) )
//							fprintf( DataFile, "TER   %5d      %.3s %c%4d \n",
//							count++, prev->GetName(), ch, prev->serno);
						
						if( aptr->flag&HeteroFlag )
						{
							fputs("HETATM",DataFile);
						}
						else
						{
							fputs("ATOM  ",DataFile);
						}
						
						std::string atname = aptr->GetName();
						if( p_save_opt_default->save_amber_pdb )
						{
							atname = MMForceField::GetAmberAtName( atname, group->GetFullName() );
						}

						if(atname.size() < 4) atname.insert(0," ");

						int k;
						for(k=0; k < 4; k++)
						{
							if(atname.size() < 4)
								atname+= " ";
						}
						
						if(atname.size() > 4)
							atname = atname.substr(0,4);
						
						for(k=0; k < 4; k++)
						{
//							if(res_name.size() < 4) res_name.insert(0," ");
							if(res_name.size() < 4) res_name+= " ";
						}
						
						if(res_name.size() > 4) res_name = res_name.substr(0,4);
						
						fprintf( DataFile, "%5d %.4s %.4s%c%4d    ", count++, atname.c_str(), res_name.c_str(),
																      chain->ident, group->serno );
						
						HaMolView* pView = GetActiveMolView();
						if( p_save_opt_default->save_transform && pView != NULL)
						{
							pView->GetTransfCoord(aptr->GetX_Ang(), aptr->GetY_Ang(), aptr->GetZ_Ang(),x,y,z);
						}
						else
						{
							x = aptr->GetX_Ang();
							y = aptr->GetY_Ang();
							z = aptr->GetZ_Ang();
						}
						fprintf(DataFile,"%8.3f%8.3f%8.3f",x,y,z);
						fprintf(DataFile,"  1.00%6.2f\n",aptr->tempf);
					}
				}
			}
			fprintf( DataFile, "TER  \n");
		}
	}
	fputs("END   \n",DataFile);
	fclose( DataFile );
	return( True );
}

int MolSet::SavePQRFile(const char* filename, bool SaveChainLetter)
{
	double x, y, z;
	HaChain  *chain;
	HaResidue  *group;
	HaAtom  *aptr;
	int count;
	
	FILE* DataFile = fopen( filename, "w" );
	if( !DataFile )
	{
			PrintLog("\n");
			return( False );
	}
	
	count = 1;
	MoleculesType::iterator mol_itr;
	for( mol_itr=HostMolecules.begin(); mol_itr != HostMolecules.end(); mol_itr++)
	{
		ChainIteratorMolecule ch_itr(*mol_itr);
		for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
		{
			ResidueIteratorChain ritr_ch(chain);
			for(group = ritr_ch.GetFirstRes(); group; group = ritr_ch.GetNextRes())
			{
				AtomIteratorAtomGroup aitr_group(group);
				for(aptr= aitr_group.GetFirstAtom(); aptr; aptr = aitr_group.GetNextAtom())
				{
					if( !p_save_opt_default->save_selected || aptr->Selected())
					{
//						if( prev && (chain->ident!=ch) )
//							fprintf( DataFile, "TER   %5d      %.3s %c%4d \n",
//							count++, prev->GetName(), ch, prev->serno);
						
						if( aptr->flag&HeteroFlag )
						{
							fputs("HETATM",DataFile);
						}
						else
							fputs("ATOM  ",DataFile);
						
						std::string atname(aptr->GetName());
						if(atname.size() < 4)
							atname.insert(0," ");

						int k;
						for(k=0; k < 3; k++)
						{
							if(atname.size() < 4)
								atname+= " ";
						}
						
						if(atname.size() > 4)
							atname = atname.substr(0,3);
						
						std::string res_name(group->GetName());

						for(k=0; k < 3; k++)
						{
							if(res_name.size() < 3)
								res_name.insert(0," ");
						}
						
						if(res_name.size() > 3)
							res_name = res_name.substr(0,2);
						
						if(SaveChainLetter)
							fprintf( DataFile, "%5d %.4s %.3s %c%4d    ",
								count++, atname.c_str(), res_name.c_str(),
								chain->ident, group->serno );
						else
							fprintf( DataFile, "%5d %.4s %.3s  %4d    ",
								 count++, atname.c_str(), res_name.c_str(),
								 group->serno );
						
						HaMolView* pView = GetActiveMolView();
						if(p_save_opt_default->save_transform && pView != NULL)
						{
							pView->GetTransfCoord(aptr->GetX_Ang(),aptr->GetY_Ang(),aptr->GetZ_Ang(),x,y,z);
						}
						else
						{
							x = aptr->GetX_Ang();
							y = aptr->GetY_Ang();
							z = aptr->GetZ_Ang();
						}
						
						fprintf(DataFile,"%8.3f%8.3f%8.3f",x,y,z);
						fprintf(DataFile,"%8.4f%7.4f\n",aptr->GetCharge(),aptr->radius);
						
					}
				}
			}
			fprintf( DataFile, "TER  \n");
		}
	}
    fputs("END   \n",DataFile);
    fclose( DataFile );
    return( True );
}

int MolSet::SavePQRFreeFile(const char* filename)
{
	double x, y, z;
	HaChain  *chain;
	HaResidue  *group;
	HaAtom  *aptr;
	int count;
	
	FILE* DataFile = fopen( filename, "w" );
	if( !DataFile )
	{
		PrintLog("\n");
		return( False );
	}
	
	count = 1;
	MoleculesType::iterator mol_itr;
	for( mol_itr=HostMolecules.begin(); mol_itr != HostMolecules.end(); mol_itr++)
	{
		ChainIteratorMolecule ch_itr(*mol_itr);
		for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
		{
			ResidueIteratorChain ritr_ch(chain);
			for(group = ritr_ch.GetFirstRes(); group; group = ritr_ch.GetNextRes())
			{
				AtomIteratorAtomGroup aitr_group(group);
				for(aptr= aitr_group.GetFirstAtom(); aptr; aptr = aitr_group.GetNextAtom())
				{
					if( !p_save_opt_default->save_selected || aptr->Selected())
					{
//						if( prev && (chain->ident!=ch) )
//							fprintf( DataFile, "TER   %5d      %.3s %c%4d \n",
//							count++, prev->GetName(), ch, prev->serno);
						
						if( aptr->flag&HeteroFlag )
						{
							fputs("HETATM",DataFile);
						}
						else
							fputs("ATOM  ",DataFile);
						
						std::string atname(aptr->GetName());
						if(atname.size() < 4)
							atname.insert(0," ");

						int k;
						for(k=0; k < 3; k++)
						{
							if(atname.size() < 4)
								atname+= " ";
						}
						
						if(atname.size() > 4)
							atname = atname.substr(0,3);
						
						std::string res_name(group->GetName());

						for(k=0; k < 3; k++)
						{
							if(res_name.size() < 3)
								res_name.insert(0," ");
						}
						
						if(res_name.size() > 3)
							res_name = res_name.substr(0,2);
						
						
						fprintf( DataFile, "%5d %.4s %.3s %c%4d    ",
										 count++, atname.c_str(), res_name.c_str(),
																					 chain->ident, group->serno );
						
						HaMolView* pView = GetActiveMolView();
						if(p_save_opt_default->save_transform && pView != NULL)
						{
							pView->GetTransfCoord(aptr->GetX_Ang(),aptr->GetY_Ang(),aptr->GetZ_Ang(),x,y,z);
						}
						else
						{
							x = aptr->GetX_Ang();
							y = aptr->GetY_Ang();
							z = aptr->GetZ_Ang();
						}
						
						fprintf(DataFile,"%16.9f%16.9f%16.9f",x,y,z);
						fprintf(DataFile,"%16.9f%16.9f\n",aptr->GetCharge(),aptr->radius);
						
					}
				}
			}
			fprintf( DataFile, "TER  \n");
		}
	}
	fputs("END   \n",DataFile);
	fclose( DataFile );
	return( True );
}

int MolSet::SaveHINToStream(std::ostream& os ) const
{      
	if( os.fail() ) return FALSE;
	char buf[128];

	HaChain* chain;
	HaResidue* pres;
	
	MoleculesType::const_iterator mol_itr;
	
	int ires = 0;
	int imol = 0;
	for (mol_itr = HostMolecules.begin(); mol_itr != HostMolecules.end(); mol_itr++)
	{
		imol++;
		std::string mol_name_full = (*mol_itr)->GetObjName();
		std::string mol_name = mol_name_full;

		while (mol_name.size() > 0)
		{
			if (std::isdigit(mol_name[mol_name.size() - 1]) || mol_name[mol_name.size() - 1] == '_')
			{
				mol_name.erase(mol_name.size() - 1);
				continue;
			}
			break;
		}

	
		const HaMolecule* pmol_c = *mol_itr;
		CAtomIntMap at_seqn_map = pmol_c->GetAtomSeqNumMap();
		ChainIteratorMolecule ch_itr(*mol_itr);
		int iat = 0;

		HaResidue* pres_fst_ch0 = NULL;
		HaChain* chain_fst = ch_itr.GetFirstChain();
		if (chain_fst) pres_fst_ch0 = chain_fst->GetFirstRes();

		if(pres_fst_ch0 != NULL && pres_fst_ch0->IsWater() && this->p_save_opt_default->save_sep_wat_mol )
		{ 
			os << "mol " << imol << " " << '\"' << "HOH" << '\"' << std::endl;
		}
		else
		{
			os << "mol " << imol << " " << '\"' << mol_name << '\"' << std::endl;
		}


		for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
		{
			ResidueIteratorChain ritr_ch(chain);
			char id_chain = chain->ident;
			HaResidue* pres_fst = ritr_ch.GetFirstRes();
			for(pres = pres_fst; pres; pres = ritr_ch.GetNextRes())
			{
				ires++;
				bool save_res_as_mol = false;
				bool save_res_info = true;

				if (this->p_save_opt_default->save_sep_wat_mol)
				{
					save_res_as_mol = true;
					save_res_info = false;
				}
				if( save_res_as_mol  && !(pres == pres_fst_ch0) && pres->IsWater() )
				{
					if(!(pres == pres_fst) ) os << "endmol " << imol << std::endl;
					imol++;
					os << "mol " << imol << " " << '\"' << pres->GetName() << '\"' << std::endl;
				}
				if (save_res_as_mol)
				{
					iat = 0;
					at_seqn_map = ((const HaResidue*)pres)->GetAtomSeqNumMap();
				}
				
				std::string res_name = pres->GetName(); 
				std::string name_mod = pres->GetNameModifier();
				std::string res_name_save = res_name;

				if( res_name == "HIS")
				{
					res_name_save = "HID";
					if( name_mod == "EPSILON") res_name_save = "HIE";
					if( name_mod == "PROT")    res_name_save = "HIP";
				}

				if( res_name == "CYS")
				{
					if( name_mod == "UNPROT") res_name_save = "CYX";
				}

				if(save_res_info) os << "res " << pres->GetSerNo() << "  " << res_name_save << "  " << pres->GetSerNo() << " - " << id_chain << std::endl;
				const HaAtom* aptr;
				std::vector<HaBond*> bonds;
				AtomIteratorAtomGroup aitr_group(pres);
				for( aptr = aitr_group.GetFirstAtom(); aptr; aptr = aitr_group.GetNextAtom())
				{
//					if( !p_save_opt_default->save_selected || aptr->Selected())
//					{	
						iat++;
						os << "atom " << iat << " " << aptr->GetName() << " " << aptr->GetStdSymbol();
						std::string ff_symbol = aptr->GetFFSymbol();	
						boost::trim(ff_symbol);
						if( ff_symbol.empty() ) ff_symbol = aptr->GetStdSymbol();
						os << " " << ff_symbol;
						os << " s ";
						double at_ch = aptr->GetCharge();
						if( fabs(at_ch) < 0.000001 || fabs(at_ch - 1.0) < 0.000001 || fabs(at_ch + 1.0) < 0.000001 || fabs(at_ch - 2.0) < 0.000001 || fabs(at_ch + 2.0) < 0.000001 )
						{
							sprintf(buf,"%2.0f",at_ch);
						}
						else if( fabs(at_ch - 0.5) < 0.000001 || fabs(at_ch + 0.5) < 0.000001 || fabs(at_ch - 1.5) < 0.000001 || fabs(at_ch + 1.5) < 0.000001 )
						{
							sprintf(buf,"%3.1f",at_ch);
						}
						else
						{
							sprintf(buf,"%9.6f",at_ch);
						}
						os << buf;
						
						double x,y,z;

//						HaMolView* pView = GetActiveMolView();
//						if(p_save_opt_default->save_transform && pView != NULL)
//						{
//							pView->GetTransfCoord(aptr->GetX_Ang(),aptr->GetY_Ang(),aptr->GetZ_Ang(),x,y,z);
//						}
//						else
//						{
							x = aptr->GetX_Ang();
							y = aptr->GetY_Ang();
							z = aptr->GetZ_Ang();
//						}

						sprintf(buf,"%16.9f %16.9f %16.9f",x,y,z);
						os << buf;
//						os << boost::str(boost::format("%16.9f %16.9f %16.9f") % (x,y,z));
						int nb = aptr->GetNBonds();
						os << " " << nb;

						HaAtom::BondIterator_const bitr = aptr->Bonds_begin();

						for( ; bitr != aptr->Bonds_end(); bitr++)
						{
							const HaBond* pbond = *bitr;
							const HaAtom* aptr_b = pbond->GetFirstAtom();
							if( aptr_b == aptr ) aptr_b = pbond->GetSecondAtom();
							std::string bond_type_str = "s";
							if( pbond->IsDouble() ) bond_type_str = "d";
							if( pbond->IsTriple() ) bond_type_str = "t";
							if( pbond->IsAromatic() ) bond_type_str = "a";
							if( pbond->IsVirtual() ) bond_type_str = "v";
							int i_bat = at_seqn_map[aptr_b];
							i_bat++;
							os << " " << i_bat << " " << bond_type_str;
						}
//					} // if selected 
					os << std::endl;	
				} // end atom
				if( save_res_info )  os << "endres " << pres->GetSerNo() << std::endl;
//				if( save_res_as_mol && !(pres == pres_fst_ch0) ) os << "endmol " << imol << std::endl;
			} //  end res
		} // end chain
		os << "endmol " << imol << std::endl;
	} // end mol

	return TRUE;
}

int MolSet::SaveXMLToStream(std::ostream& os, const AtomSaveOptions* popt_arg ) const
{
	if( os.fail() ) return FALSE;

	std::auto_ptr<AtomSaveOptions>popt_auto( popt_arg == NULL ? (AtomSaveOptions*) p_save_opt_default->clone() : (AtomSaveOptions*) popt_arg->clone() );
	AtomSaveOptions* popt = popt_auto.get();

	bool save_header = popt->ToSaveHeader();
	bool save_footer = popt->ToSaveFooter();

	popt->SetSaveHeader(false);
	popt->SetSaveFooter(false);

	char buf[256];

	HaChain* chain;
	HaResidue* pres;
	vector<HaAtom*>::iterator paitr;
	
	CAtomIntMap at_seqn_map  = GetAtomSeqNumMap();

	if( save_header )
	{
		os << harlem::StdXMLHeader()     << std::endl;
		os << harlem::HarlemDataHeader() << std::endl;
	}

	os << "<molset ";
	os << "name=\"" << GetName() << "\" ";
	os << ">" << std::endl;

	int seq_n = 0;

	MoleculesType::const_iterator mol_itr;
	for( mol_itr=HostMolecules.begin(); mol_itr != HostMolecules.end(); mol_itr++)
	{
		os << "<molecule ";
		os << "name=\"" << (*mol_itr)->GetObjName() << "\" ";
		os << "/>" << std::endl;
		
		ChainIteratorMolecule ch_itr(*mol_itr);
		for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
		{
			os << "<chain ";
	        os << " id=\"" << chain->ident << "\" ";
			os << "/>" << std::endl;
				
			ResidueIteratorChain ritr_ch(chain);
			for(pres = ritr_ch.GetFirstRes(); pres; pres = ritr_ch.GetNextRes())
			{
				os << "<res ";
				std::string res_name = pres->GetName();
				os << "name=\"" << res_name << "\" ";
				
				std::string name_mod = pres->GetNameModifier();
				if( !name_mod.empty() ) os << "name_mod=\"" << name_mod << "\" ";

				sprintf(buf,"%d",pres->GetSerNo());
				std::string res_ser_no = buf;
				if( !res_ser_no.empty() ) os << " no=\"" << res_ser_no << "\" ";
				
				os << "/> " << std::endl;

				const HaAtom* aptr;
				
				int i_in_res = 0;
				AtomIteratorResidue aitr(pres);
				for( aptr = aitr.GetFirstAtom() ; aptr ; aptr = aitr.GetNextAtom() )
				{
					i_in_res++;
					seq_n++;

					os << "  <atom";
					
					if( aptr->GetHybrid() != NO_HYBRID) 
					{
						os << " hybr=\"";
						if ( aptr->GetHybrid() == SP_HYBRID ) os << "SP\" ";
						else if (  aptr->GetHybrid() == SP2_HYBRID ) os << "SP2\" ";
						else if (  aptr->GetHybrid() == SP3_HYBRID ) os << "SP3\" ";
					}

					std::string da_symbol = "";
					if(aptr->IsHBDonor())    da_symbol+= "D";
					if(aptr->IsHBAcceptor()) da_symbol+= "A";
					if( !da_symbol.empty() ) os << " hb_st=\"" << da_symbol << "\" ";
					
					std::string ff_symbol = aptr->GetFFSymbol();
					if( !ff_symbol.empty() ) os << " ff_s=\"" << ff_symbol << "\" ";

					double chrg = aptr->GetCharge();
					sprintf(buf," chrg=\"%10.5f\" ",chrg);
					if( fabs(chrg) > 1.0e-8 ) os << buf;

					double mass = aptr->GetMass();
					if( fabs(mass - aptr->GetStdMass()) > 0.001 ) 
					{
						sprintf(buf," mass=\"%10.5f\" ", mass );   
						os << buf;
					}

					os << "> ";

					double x, y, z;
					const HaMolView* pView = GetActiveMolView();
					if( p_save_opt_default->save_transform && pView != NULL )
					{
						pView->GetTransfCoord(aptr->GetX(), aptr->GetY(), aptr->GetZ(),x,y,z);
					}
					else
					{
						x = aptr->GetX();
						y = aptr->GetY();
						z = aptr->GetZ();
					}

					std::string at_name = aptr->GetName();
					if( at_name.empty() ) 
					{
						at_name = aptr->GetStdSymbol();
						boost::to_upper( at_name );
						at_name += harlem::ToString(i_in_res);
					}
					
					sprintf(buf,"%4d %5s  %d %16.9f %16.9f %16.9f ",seq_n, at_name.c_str(), aptr->GetElemNo(), x,y,z);
					os << buf;

					if( aptr->GetNBonds() > 0 )
					{	
						os << "<bonds> ";

						HaAtom::BondIterator_const bitr = aptr->Bonds_begin();
						for( ; bitr !=  aptr->Bonds_end(); bitr++ )
						{
							const HaBond* pb = *bitr;
							const HaAtom* aptr_b = pb->GetFirstAtom(); 
							if( aptr_b == aptr ) aptr_b = pb->GetSecondAtom(); 

							if( at_seqn_map.count(aptr_b) > 0 )
							{
								int seq_n_b = at_seqn_map[aptr_b] + 1; 
								os << seq_n_b;
								if( pb->IsAromatic() ) os << "@AR";
								else if( pb->IsDouble() ) os << "@D";
								else if( pb->IsTriple() ) os << "@T";
								os << "  ";
							}
						}
						
						os << "</bonds>";
					}
					if( fabs(aptr->tempf) > DBL_EPSILON )
					{
						sprintf(buf,"%9.4f",aptr->tempf);
						os << "<tempf>" << buf << "</tempf>";
					}
					os << "</atom>" << std::endl;
				}
			}
		}
	}
	
	if( per_bc->IsSet() )
	{
		os << "<ucell> ";
		sprintf(buf," %16.9f %16.9f %16.9f %16.9f %16.9f %16.9f ",per_bc->GetA(),per_bc->GetB(),per_bc->GetC(),
			per_bc->GetAlpha()*RAD_TO_DEG,per_bc->GetBeta()*RAD_TO_DEG, per_bc->GetGamma()*RAD_TO_DEG ) ;
		os << buf << "</ucell>" << std::endl;
	}

	if( !p_zmat->IsEmpty() )
	{
		p_zmat->SaveXMLToStream(os, popt);
	}

	const ChemGroup* gptr;
	ChemGroupsType::const_iterator gitr = ChemGroups.begin();
	for(; gitr != ChemGroups.end() ; gitr++ )
	{
		gptr = &(*gitr);
		int ng=gptr->size();
		std::string id = gptr->GetID();
		double protect = gptr->GetProtect();
		os << "<chem_grp ";
		if( !name_mset.empty() ) os << "id=\"" << id << "\"";
		sprintf(buf," prot=\"%8.5f\"  size=\"%.1d\" ",gptr->GetProtect(),ng);
		os << buf << ">" << std::endl;
		
		const HaAtom* aptr;

		int num_in_line = 0;
		AtomIteratorAtomGroup_const aitr_chem_group(gptr);
		for( aptr= aitr_chem_group.GetFirstAtom(); aptr; aptr= aitr_chem_group.GetNextAtom())
		{
			aptr->FillRef(buf);
			os << " " << buf << " ";
			num_in_line++;
			if(num_in_line == 10 )
			{
				os << "  " << std::endl; 
				num_in_line = 0;
			}	
		}
		if( num_in_line != 0 ) os << " " << std::endl;

		os << "</chem_grp>" << std::endl;
	}

	const AtomGroup* lptr;
	AtomGroupIteratorMolSet_const litr(this);
	for(lptr = litr.GetFirst(); lptr; lptr = litr.GetNext())
	{
		os << "<atgrp id=\"" << lptr->GetID() << "\"";
		int ng=lptr->size();
		sprintf(buf," size=\"%2d\" ",ng);
		os << buf;
		os << ">" << std::endl;

		int num_in_line = 0;
		const HaAtom* aptr;

		const HaMolecule* pmol = NULL; 

		AtomIteratorAtomGroup_const aitr(lptr);
		for(aptr= aitr.GetFirstAtom(); aptr; aptr=aitr.GetNextAtom())
		{
			if( aptr->GetHostMol() != pmol )
			{
				pmol = aptr->GetHostMol();
				if( num_in_line != 0) os << std::endl;
				os << "  molname=" << pmol->GetObjName() << " " << std::endl;	
			}

			aptr->FillRef(buf,HaAtom::ATOMREF_NO_MOL);
			os << " " << buf << " ";
			num_in_line++;

			if(num_in_line == 10 )
			{
				os << " " << std::endl;
				num_in_line = 0;
			}
		}
		os << " " << std::endl;
		os << "</atgrp>" << std::endl;
	}
	if(!ChargeMaps.empty())
	{
		int nm = ChargeMaps.size();
		int i;
		for(i = 0; i< nm; i++)
		{
			os << "<chrg_map ";
			std::string name = ChargeMaps[i].GetName();
			if( !name.empty() ) os << "name=\"" << name << "\" ";
			
			os << ">" << std::endl;

			AtomDoubleMap::const_iterator mitr;
			mitr = ChargeMaps[i].begin();
			for(; mitr != ChargeMaps[i].end(); mitr++)
			{
				HaAtom* aptr = (HaAtom*) (*mitr).first;
				double ch    = (*mitr).second;
					
				aptr->FillRef(buf);
				os << " " << buf;
				sprintf(buf," %16.9f ",ch);
				os << buf << std::endl;
			}
			os << "</chrg_map>" << std::endl;
		}
	}

	SaveCrdSnapshots(os);

	os << "</molset>" << std::endl;
	
	int nm = CompModules.size();
	int im;
	for( im = 0; im < nm; im++)
	{
		CompModules[im]->SaveXMLToStream(os, popt);
	}

	if( save_footer )
	{
		os << harlem::HarlemDataFooter() << std::endl;
	}
	
	return FALSE;
}



int MolSet::SaveCrdSnapshots(std::ostream& os, const harlem::HashMap* popt_par ) const
{
	std::auto_ptr<harlem::HashMap> popt_auto( (popt_par == NULL) ? new harlem::HashMap() : popt_par->clone() );
	harlem::HashMap* popt = popt_auto.get();

	char buf[256];

	if(!crd_snapshots.empty())
	{
		int save_header = popt->get_i("SAVE_HEADER");
		int save_footer = popt->get_i("SAVE_FOOTER");
		int save_mset_name = popt->get_i("SAVE_MSET_NAME");

		if( save_header )
		{
			os << harlem::StdXMLHeader()     << std::endl;
			os << harlem::HarlemDataHeader() << std::endl;
		}

		int ns = crd_snapshots.size();
		int i;
		for(i = 0; i< ns; i++)
		{
			os << "<crd_snap ";
			std::string name = crd_snapshots[i]->GetName();
			if( !name.empty() ) os << "name=\"" << name << "\" ";
			if( save_mset_name )
			{
				os << "mset_name=\"" << this->GetName() << "\" ";
			}

			AtomGroup* p_at_grp = dynamic_cast<AtomGroup*>(crd_snapshots[i]->p_at_cont);
			if( p_at_grp != NULL )
			{
				os << "atgrp=\"" << p_at_grp->GetID() << "\" ";
			}
			
			os << ">" << std::endl;

			std::string desc = crd_snapshots[i]->GetDesc();
			if( !desc.empty() )
			{
				os << "<desc>" << desc << "</desc>" << std::endl;
			}

			HaVec_double crd_snap = crd_snapshots[i]->GetCrd();
			
			write_double_array_chuncks(os,crd_snap,6,FLOAT_F12_7);
			
//			if( crd_snap.size() % 6 != 0 ) os << std::endl; 

			if( crd_snapshots[i]->HasPBox() )
			{
				HaVec_double pbox_snap = crd_snapshots[i]->GetPBox();
				if( pbox_snap.size() >= 6 )
				{
					pbox_snap[3] *= RAD_TO_DEG;
					pbox_snap[4] *= RAD_TO_DEG;
					pbox_snap[5] *= RAD_TO_DEG;
				}
				os << "<pbox>";
				int i;
				for( i = 0; i < pbox_snap.size(); i++ )
				{
					sprintf(buf," %12.6f",pbox_snap[i]); os << buf;
				}
				os << "</pbox>" << std::endl;
			}
			os << "</crd_snap>" << std::endl;
		}
		if( save_footer )
		{
			os << harlem::HarlemDataFooter() << std::endl;
		}
	}	
	return TRUE;
}

int MolSet::SaveCrdSnapshots(const std::string& fname, const harlem::HashMap* popt_par ) const
{
	std::auto_ptr<harlem::HashMap> popt_auto( (popt_par == NULL) ? new harlem::HashMap() : popt_par->clone() );
	harlem::HashMap* popt = popt_auto.get();

	if( popt_par == NULL )
	{
		popt->set_i("SAVE_HEADER",1);
		popt->set_i("SAVE_MSET_NAME",1);
		popt->set_i("SAVE_FOOTER",1);
	}
	
	std::ofstream os(fname.c_str() );
	if( os.fail() )
	{
		PrintLog(" Error in MolSet::SaveCrdSnapshots() \n");
		PrintLog(" Error to open file %s\n", fname.c_str());
		return FALSE;
	}
	int ires = SaveCrdSnapshots(os,popt);
	return ires;
}

int  MolSet::LoadCrdSnapshots(const std::string& fname, const harlem::HashMap* popt )
{
	using namespace rapidxml;
	std::ifstream is(fname.c_str(), ios::binary);
	if(is.fail()) 
	{
		PrintLog(" Error in MolSet::LoadCrdSnapshots \n");
		PrintLog(" Error to open file %s\n",fname.c_str());
		return FALSE;
	}
	int ires = LoadXMLStream (is, NULL );

	return FALSE;

}


int MolSet::LoadXMLNode( rapidxml::xml_node<>* node_mset, const AtomLoadOptions* popt_arg )
{
	using namespace rapidxml;

//	PrintLog( "MolSet::LoadXMLNode() pt 1 \n");

	std::auto_ptr<AtomLoadOptions> popt_auto( popt_arg == NULL ? (AtomLoadOptions*) p_load_opt_default->clone() : (AtomLoadOptions*) popt_arg->clone() );
	AtomLoadOptions* popt = popt_auto.get();
 
	try
	{
//		PrintLog( "MolSet::LoadXMLNode() pt 2 \n");

		char buf[256];
		char attr_buf[256];

		std::string tag = node_mset->name();
		if( !boost::iequals(tag, "molset") ) throw std::runtime_error(" Name of the node is not molset ");

//		PrintLog( "MolSet::LoadXMLNode() pt 3 \n");

		xml_attribute<>* attr;
		for( attr = node_mset->first_attribute(); attr; attr = attr->next_attribute() )
		{
			std::string tag_attr = attr->name();
			if( boost::iequals(tag_attr, "name") )
			{
				this->SetName(attr->value());
			}
		}

//		PrintLog( "MolSet::LoadXMLNode() pt 4 \n");

		HaMolecule* pmol   = NULL;
		HaChain*    pchain = NULL;
		HaResidue*  pres   = NULL;
		int res_ser_no = 0;

		std::map<std::string, HaAtom* > id_at_map; // map of atom IDs to atom pointers
		std::map<HaAtom*, std::vector<std::string> > at_bond_str_map;

		xml_node<>* node1 = node_mset->first_node();
		for( ; node1; node1 = node1->next_sibling() )
		{
			tag = node1->name();
			boost::to_lower(tag);
//			PrintLog( "MolSet::LoadXMLNode() pt 5   name = %s \n", tag.c_str());
			
			if( boost::equals(tag, "atom")  )
			{
				if( pmol == NULL     ) pmol = AddNewMolecule();
				if( pchain == NULL  ) 
				{
					pchain = pmol->AddChain(' ');
					res_ser_no = 0;
				}
				if( pres == NULL ) 
				{
					res_ser_no++;
					pres = pchain->AddResidue( res_ser_no ); 
				}
				HaAtom* aptr = pres->AddNewAtom();

//				for( int i = 0; i < 100; i++)
//				{
				for( attr = node1->first_attribute(); attr; attr = attr->next_attribute() )
				{		
//					strcpy(attr_buf,attr->name());
//					char* ch_ptr;
//					ch_ptr = &attr_buf[0];
//					for(; *ch_ptr != 0; ch_ptr++)
//					{
//						*ch_ptr = (char)::tolower((int)*ch_ptr);
//					}
					std::string tag_attr = attr->name();
					boost::to_lower(tag_attr);
					
					if( boost::equals(tag_attr, "chrg") ) 
// 		            if( attr_buf[0] == 'c' )
					{
						std::string chrg_str = attr->value();
						aptr->SetCharge( atof(chrg_str.c_str()));
					}
					else if( boost::equals(tag_attr, "ff_s") )
//					else if( attr_buf[0] == 'f'  )
					{
						aptr->SetFFSymbol( attr->value() );
					}
					else if( boost::equals(tag_attr, "hb_st") )
//					else if( attr_buf[0] == 'h'  )
					{
						aptr->SetHBStatus( attr->value() );
					}
					else if( boost::equals(tag_attr, "mass") ) 
//					else if( attr_buf[0] == 'm'  )
					{
						std::string mass_str = attr->value();
						aptr->SetMass( atof(mass_str.c_str()));
					}
					else if( boost::equals(tag_attr, "name") )
//					else if( attr_buf[0] == 'n'  )
					{
						aptr->SetNameFast( attr->value() );
					}
					else if( boost::equals(tag_attr, "hybr") ) 
//					else if( attr_buf[0] == 'h'  )
					{
						aptr->SetHybrid( attr->value() );
					}

					else if( boost::equals(tag_attr, "id") )
//					else if( attr_buf[0] == 'i'  )
					{
						std::string id_str = attr->value();
						if( id_at_map.count(id_str) != 0 ) throw std::runtime_error(" atom_id = " + id_str + " is not unique ");
						if( !id_str.empty() ) id_at_map[id_str] = aptr;
					}
				}
//				}

				std::string val_str = node1->value();
				boost::trim(val_str);
				std::vector<std::string> str_vec;
				boost::split(str_vec,val_str,boost::is_any_of(" ,"),boost::token_compress_on);
				if( str_vec.size() == 6 )
				{
					std::string id_str = str_vec[0];
					if( id_at_map.count(id_str) != 0 ) throw std::runtime_error(" atom_id = " + id_str + " is not unique ");
					id_at_map[id_str] = aptr;
					
					aptr->SetNameFast(str_vec[1]);

					std::string elem_str = str_vec[2];
					if( harlem::IsInt( elem_str) ) 
					{
						aptr->SetElemNo( atoi(elem_str.c_str()) );
					}
					else 
					{
						int elem_no = HaAtom::GetElemNoFromName(elem_str);
						aptr->SetElemNo( elem_no );
					}
					
					aptr->SetX( atof(str_vec[3].c_str()) );
					aptr->SetY( atof(str_vec[4].c_str()) );
					aptr->SetZ( atof(str_vec[5].c_str()) );
					
				}
				else
				{
					throw std::runtime_error("value of <atom> element should contain 6 strings \n" + val_str );
				}

				xml_node<> *node2 = node1->first_node();
				for( ; node2; node2 = node2->next_sibling() )
				{
				//	std::string tag_node2 = node2->name();
				//	if( boost::iequals(tag_node2, "bonds") )
					if( node2->name()[0] == 'b')
					{
						std::string bonds_str = node2->value();
					
						/*boost::tokenizer<> tok(bonds_str);
						boost::tokenizer<>::iterator titr;
						for( titr = tok.begin(); titr!=tok.end();++titr)
						{
							if( id_at_map.count(*titr) == 0) continue;
							HaAtom* aptr2 = id_at_map[*titr];
							if( aptr->IsBonded(*aptr2) ) continue;
							HaBond* pbnd = this->AddBond(aptr,aptr2);
						}*/
						
						boost::trim( bonds_str );
						std::vector<std::string> str_vec;
						str_vec.reserve(4);
						boost::split(str_vec,bonds_str,boost::is_any_of(" ,"),boost::token_compress_on);
						at_bond_str_map[aptr] = str_vec;
					}
					if( node2->name()[0] == 't')
					{
						if( node2->name()[1] != 'e') continue;
						if( node2->name()[2] != 'm') continue;
						if( node2->name()[3] != 'p') continue;

						std::string tempf_str = node2->value();

					
						/*boost::tokenizer<> tok(bonds_str);
						boost::tokenizer<>::iterator titr;
						for( titr = tok.begin(); titr!=tok.end();++titr)
						{
							if( id_at_map.count(*titr) == 0) continue;
							HaAtom* aptr2 = id_at_map[*titr];
							if( aptr->IsBonded(*aptr2) ) continue;
							HaBond* pbnd = this->AddBond(aptr,aptr2);
						}*/
						
						boost::trim( tempf_str );
						double tempf_loc = harlem::ToDouble(tempf_str);
						aptr->tempf = tempf_loc;
					}
				}
				std::string at_ref_str = aptr->GetRef();
				if( id_at_map.count(at_ref_str) != 0 )
				{
					PrintLog("Warning in MolSet::LoadXML() : Atom Reference %s is not unique \n", at_ref_str.c_str());
				}
				else
				{
					id_at_map[at_ref_str] = aptr;
				}
				continue;
			}  // end if(tag == "atom") 

			if( boost::equals(tag, "res")  )
			{
				if( pmol == NULL ) pmol = AddNewMolecule();
				if( pchain == NULL ) pchain = pmol->AddChain(' ');

				res_ser_no++;
				std::string res_name = "RES";
				std::string name_mod = "";

				for( attr = node1->first_attribute(); attr; attr = attr->next_attribute() )
				{		
					std::string tag_attr = attr->name();
					if( boost::iequals(tag_attr, "name") ) res_name = attr->value();
					if( boost::iequals(tag_attr, "name_mod") ) name_mod = attr->value();
					if( boost::iequals(tag_attr, "no") ) 
					{
						std::string ser_no_str = attr->value();
						if( harlem::IsInt(ser_no_str) ) res_ser_no = atoi(ser_no_str.c_str());
					}
				}
				pres = pchain->AddResidue( res_ser_no );
				pres->SetName( res_name.c_str());
				if( !name_mod.empty() )pres->SetNameModifier(name_mod.c_str());
				continue;
			} // end if(tag == "res")

			if( boost::equals(tag, "chain")  )
			{
				if( pmol == NULL ) pmol = AddNewMolecule();
				char id_chain = ' ';
				for( attr = node1->first_attribute(); attr; attr = attr->next_attribute() )
				{		
					std::string tag_attr = attr->name();
					if( boost::iequals(tag_attr, "id") ) 
					{
						std::string id_str = attr->value();
						if( !id_str.empty() ) id_chain = id_str[0];
					}
				}
				pchain = pmol->AddChain( id_chain );
				res_ser_no = 0;
				continue;
			} // end if( tag == "chain")

			if( boost::equals(tag, "molecule")  )
			{
				pmol = AddNewMolecule();
				for( attr = node1->first_attribute(); attr; attr = attr->next_attribute() )
				{		
					std::string tag_attr = attr->name();
					if( boost::iequals(tag_attr, "name") )
					{
						pmol->SetObjName(attr->value());
					}
				}
				continue;
			} // end if(tag == "molecule")

			if( boost::equals(tag, "ucell")  )
			{
				std::string val_str = node1->value();
				boost::trim(val_str);
				std::vector<std::string> str_vec;
				boost::split(str_vec,val_str,boost::is_any_of(" ,"),boost::token_compress_on);
				
				double px = 1.0;
				double py = 1.0;
				double pz = 1.0;
				
				double p_alpha = DEG_TO_RAD * 90.0;
				double p_beta  = DEG_TO_RAD * 90.0;
				double p_gamma = DEG_TO_RAD * 90.0;
				
				if( str_vec.size() < 3 )
				{
					PrintLog(" Element <ucell> contains less than 3 strings: periodical boundary conditions will not be set ");
					continue;
				}
				
				px = atof( str_vec[0].c_str() );
				py = atof( str_vec[1].c_str() );
				pz = atof( str_vec[2].c_str() );
				
				if( str_vec.size() > 5 )
				{
					p_alpha = DEG_TO_RAD * atof( str_vec[3].c_str() );
					p_beta  = DEG_TO_RAD * atof( str_vec[4].c_str() );
					p_gamma = DEG_TO_RAD * atof( str_vec[5].c_str() );
				}

				per_bc->SetBox(px,py,pz,p_alpha,p_beta,p_gamma);
				continue;
			} // if( tag == "ucell" )

			if( tag == "zmat" )
			{
				ZMatCrd* pzm = this->GetZMat();
				pzm->Clear();
				int ires = pzm->LoadXMLNode(node1,popt);
				if(!ires)
				{
					delete pzm;
					pzm = NULL;
				}
				continue;
			} // if ( tag == "zmat" )

			if( tag == "atgrp" )
			{
				std::string grp_id = "GRP";
				int nsize = 10;
				for( attr = node1->first_attribute(); attr; attr = attr->next_attribute() )
				{		
					std::string tag_attr = attr->name();
					if( boost::iequals(tag_attr, "id") ) grp_id = attr->value();
					else if( boost::iequals(tag_attr, "size") ) 
					{
						std::string size_str = attr->value();
						nsize = atoi( size_str.c_str() );
					}
				}
				AtomGroup* pgrp = AddAtomGroup( grp_id.c_str() ) ;
				pgrp->reserve(nsize);

				std::string val_str = node1->value();
				boost::trim(val_str);
				std::vector<std::string> str_vec;
				boost::split(str_vec,val_str,boost::is_any_of(" ,\r\n"),boost::token_compress_on);

				int n = str_vec.size();
				int i;

				std::string mol_name_def = "";

				for( i = 0; i < n; i++ )
				{
					HaAtom* aptr_g = NULL;
					std::string token = str_vec[i];

					size_t eq_pos = token.find("molname=");
					if( eq_pos != std::string::npos )
					{
						mol_name_def = boost::to_upper_copy(token.substr(8));
						continue;
					}

					std::string at_id_str = str_vec[i];
					std::string at_id_full;

					if( at_id_str.find('$') != std::string::npos || mol_name_def.empty() )
					{
						at_id_full = at_id_str;
					}
					else
					{
						at_id_full = "$" + mol_name_def + "$" + at_id_str;
					}
					
					if( id_at_map.count(at_id_full) != 0 )
					{
						aptr_g = id_at_map[at_id_full];
					}
					else
					{
						aptr_g = GetAtomByRef( at_id_full.c_str() );
					}
					if( aptr_g != NULL ) pgrp->InsertAtom(aptr_g);
				} // end for(i = 0; i < n; i++ )
				continue;
			} // end if( tag == "atgrp") 

			if( tag == "crd_snap" )
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
			} // end if( tag == "crd_snap") 
		} 		
		
		std::map<HaAtom*, std::vector<std::string> >::iterator bitr;
		for( bitr = at_bond_str_map.begin(); bitr != at_bond_str_map.end(); bitr++ )
		{
			HaAtom* aptr1 = (*bitr).first;
			std::vector<std::string>& bnd_at_ids = (*bitr).second;
			int nb = bnd_at_ids.size();
			for( int ib = 0; ib < nb; ib++)
			{
				std::string& bstr = bnd_at_ids[ib];
				size_t sep_pos = bstr.find('@');
				
				std::string at_id_str = bstr.substr(0,sep_pos);
				if( id_at_map.count(at_id_str) == 0) 
				{
					PrintLog(" Unknown atom ID in bond description %s \n", at_id_str.c_str() );
					continue;
				}

				HaAtom* aptr2 = id_at_map[at_id_str];
				if( aptr1->IsBonded(*aptr2) ) continue;

				HaBond* pbnd = this->AddBond(aptr1,aptr2);
				if( pbnd == NULL ) continue;
				if( sep_pos != std::string::npos )
				{
					std::string bord_str = bstr.substr(sep_pos+1);
					if( boost::iequals(bord_str,"AR")) pbnd->SetAromatic();
					if( boost::iequals(bord_str,"D")) pbnd->SetDouble();
					if( boost::iequals(bord_str,"T")) pbnd->SetTriple();
				}
			}
		}
	}
	catch( const std::exception& ex)
	{
		PrintLog(" Error in MolSet::LoadXMLNode() \n");
		PrintLog("%s\n",ex.what());
		return FALSE;
	}
	return TRUE;
}

ZMatCrd* MolSet::GetZMat( const harlem::HashMap* popt )
{
	return p_zmat;
}

int MolSet::SaveOldHarlemStream(std::ostream& os)
{
	if( os.fail() ) return FALSE;

	HaAtom  *aptr;
    HaBond  *bptr;
    const char *ptr;

	double x, y, z;

	char buf[256];

	HaChain* chain;
	HaResidue* pres;
	vector<HaAtom*>::iterator paitr;
	AtomIntMap at_id_map;
	AtomIntMap::iterator mitr; 

	os << "#BEGIN MOL_SET " << std::endl;
	if(!name_mset.empty())
	{
		os << "MOL_SET_NAME=" << name_mset << std::endl;
	}

	MoleculesType::iterator mol_itr;
	for( mol_itr=HostMolecules.begin(); mol_itr != HostMolecules.end(); mol_itr++)
	{
		os << "#MOLECULE " << std::endl;
		std::string mol_name = (*mol_itr)->GetObjName();
		at_id_map.clear();
		int atid = 1;
		if(!mol_name.empty())
		{
			os << "MOLNAME=" << mol_name << std::endl;
		}
		
		os << "#ATOMS " << std::endl;

		ChainIteratorMolecule ch_itr(*mol_itr);
		for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
		{
			ResidueIteratorChain ritr_ch(chain);
			for(pres = ritr_ch.GetFirstRes(); pres; pres = ritr_ch.GetNextRes())
			{
				std::string full_res_name = pres->GetFullName();
				for(paitr=pres->begin(), aptr=*paitr; paitr != pres->end(); paitr++, aptr=*paitr)
				{
					sprintf(buf,"%5d ",atid); os << buf;
					sprintf(buf,"%3d ",aptr->GetElemNo()); os << buf;
					sprintf(buf,"\"%.4s\" ",aptr->GetName()); os << buf;

					at_id_map[aptr] = atid;
					atid++;

					HaMolView* pView = GetActiveMolView();
					if( p_save_opt_default->save_transform && pView != NULL )
					{
						pView->GetTransfCoord(aptr->GetX(), aptr->GetY(), aptr->GetZ(),x,y,z);
					}
					else
					{
						x = aptr->GetX();
						y = aptr->GetY();
						z = aptr->GetZ();
					}
					
					sprintf(buf,"%16.9f %16.9f %16.9f",x,y,z);
					os << buf;
					sprintf(buf," \"%.1s\" ",&chain->ident);
					os << buf;
					
					sprintf(buf,"\"%s\"  %6d  %10.5f ",full_res_name.c_str(),pres->serno,aptr->GetCharge());
					os << buf;
					if( aptr->GetHybrid() == NO_HYBRID) 
					{
						sprintf(buf,"\"%s\" ","NO_HYBRID");
						os << buf;
					}
					else if ( aptr->GetHybrid() == SP_HYBRID )
					{
						sprintf(buf,"\"%s\" ","SP");
						os << buf;
					}
					else if (  aptr->GetHybrid() == SP2_HYBRID )
					{
						sprintf(buf,"\"%s\" ","SP2");
						os << buf;
					}
					else if (  aptr->GetHybrid() == SP3_HYBRID )
					{
						sprintf(buf,"\"%s\" ","SP3");
						os << buf;
					}

					std::string da_symbol;
					if(aptr->IsHBDonor())    da_symbol+= "D";
					if(aptr->IsHBAcceptor()) da_symbol+= "A";
					
					sprintf(buf,"\"%s\" ",da_symbol.c_str()   ); os << buf;
					sprintf(buf,"\"%s\" ",aptr->GetFFSymbol() ); os << buf;
					sprintf(buf," %10.5f ",aptr->GetMass() );    os << buf;
					os << std::endl;
				}
			}
		}
		os << "#END ATOMS " << std::endl;
		
		os << "#BONDS " << std::endl;
		
		AtomIteratorMolecule aitr_m(*mol_itr);
		for( aptr = aitr_m.GetFirstAtom(); aptr; aptr = aitr_m.GetNextAtom() )
		{
			HaAtom::BondIterator bitr = aptr->Bonds_begin();
			for(; bitr != aptr->Bonds_end(); bitr++)
			{
				bptr = (*bitr);
				mitr = at_id_map.find(bptr->srcatom);
				if(mitr == at_id_map.end() ) continue;
				int iat1 = (*mitr).second;
				mitr = at_id_map.find(bptr->dstatom);
				if(mitr == at_id_map.end() ) continue;
				int iat2 = (*mitr).second;
				if( iat2 < iat1 ) continue;

				sprintf(buf,"%5d %5d  ",iat1,iat2);
				os << buf;

				std::string bond_flag_str = "SINGLE";
				if( bptr->flag & AromBondFlag )
				{
					bond_flag_str = "AROMATIC";
				}
				else if( bptr->flag & TripBondFlag )
				{
					bond_flag_str = "TRIPLE";
				}
				else if( bptr->flag & DoubBondFlag )
				{
					bond_flag_str = "DOUBLE";
				}

				os << bond_flag_str << std::endl;
			}
		}
		os << "#END BONDS " << std::endl;
		os << "#END MOLECULE " << std::endl;
	}
	
	if( per_bc->IsSet() )
	{
		os << "#UNIT CELL " << std::endl;
		sprintf(buf," %16.9f %16.9f %16.9f %16.9f %16.9f %16.9f ",per_bc->GetA(),per_bc->GetB(),per_bc->GetC(),
			per_bc->GetAlpha()*RAD_TO_DEG,per_bc->GetBeta()*RAD_TO_DEG, per_bc->GetGamma()*RAD_TO_DEG ) ;
		os << buf << std::endl;
		os << "#END UNIT CELL " << std::endl;
	}

	if(!ChargeMaps.empty())
	{
		os << "#CHARGE MAPS " << std::endl;
		int nm = ChargeMaps.size();
		int i;
		for(i = 0; i< nm; i++)
		{
			os << "#BEGIN MAP " << ChargeMaps[i].GetName() << std::endl;
			AtomDoubleMap::iterator mitr;
			mitr = ChargeMaps[i].begin();
			for(; mitr != ChargeMaps[i].end(); mitr++)
			{
				HaAtom* aptr = (HaAtom*) (*mitr).first;
				double ch    = (*mitr).second;
					
				aptr->FillRef(buf);
				os << " " << buf;
				sprintf(buf," %16.9f ",ch);
				os << buf << std::endl;
			}
			os << "#END MAP " << std::endl;
		}
		os << "#END CHARGE MAPS " << std::endl;
	}

	os << "#GROUPS " << std::endl;
	
	ChemGroup* gptr;
	ChemGroupIterator gitr(this);
	for(gptr = gitr.GetFirst(); gptr; gptr = gitr.GetNext())
	{
		int ng=gptr->size();
		sprintf(buf," %s  %8.5f  %5d ",gptr->GetID(),gptr->GetProtect(),ng);
		os << buf;
		
		int num_in_line = 0;
		AtomIteratorAtomGroup aitr_chem_group(gptr);
		for(aptr= aitr_chem_group.GetFirstAtom(); aptr; aptr= aitr_chem_group.GetNextAtom())
		{
			aptr->FillRef(buf);
			os << " " << buf << " ";
			num_in_line++;
			if(num_in_line == 30 )
			{
				os << " ### " << std::endl << "         "; 
				num_in_line = 0;
			}
		}
		os << " " << std::endl;
	}
	os << "#END GROUPS " << std::endl;

	os << "#ATOM LISTS " << std::endl;
	AtomGroup* lptr;
	AtomGroupIteratorMolSet litr(this);
	for(lptr = litr.GetFirst(); lptr; lptr = litr.GetNext())
	{
		os << " " << lptr->GetID() << " ";
		int ng=lptr->size();
		sprintf(buf," %6d ",ng);
		os << buf;
		int num_in_line = 0;
		AtomIteratorAtomGroup aitr(lptr);
		for(aptr= aitr.GetFirstAtom(); aptr; aptr=aitr.GetNextAtom())
		{
			if(num_in_line == 30 )
			{
				os << " ### " << std::endl << "         ";
				num_in_line = 0;
			}
			aptr->FillRef(buf);
			os << " " << buf << " ";
			num_in_line++;
		}
		os << " " << std::endl;
	}
	os << "#END ATOM LISTS " << std::endl;

	HaQCMod* ptr_qc_mod = GetQCMod(false);
	if(ptr_qc_mod != NULL)
	{
		os << "#QCHEM MODULE " << std::endl;
		os << "BASIS=" << ptr_qc_mod->GetBasName() << std::endl;
		os << "LOC_ORB_BASIS=" << ptr_qc_mod->GetLocOrbSetID() << std::endl;
		os << "#END QCHEM MODULE " << std::endl;
	}

//	ETCouplMod* ptr_et_coupl = GetETCouplMod(false);
//	if( ptr_et_coupl != NULL)
//	{
//		fprintf(DataFile,"#ET COUPLING MODULE \n");
//		fprintf(DataFile,"#END ET COUPLING MODULE \n");
//	}

	HaMolMechMod* ptr_mm_mod = GetMolMechMod(false);
	if( ptr_mm_mod != NULL )
	{
		MolMechModel* p_mm_model = ptr_mm_mod->GetMolMechModel();
		os << "#MOLECULAR MECHANICS MODULE " << std::endl;
		
		vector<MMDihedral>::iterator iditr;
		
		if( !p_mm_model->ImprDihedrals.empty())
		{
			os << "IMPROPER ANGLES" << std::endl;
			for(iditr = p_mm_model->ImprDihedrals.begin(); iditr != p_mm_model->ImprDihedrals.end(); iditr++)
			{
				MMDihedral& impr_dihedral = *iditr;
				if( impr_dihedral.pt1 == NULL || impr_dihedral.pt2 == NULL || 
					impr_dihedral.pt3 == NULL || impr_dihedral.pt4 == NULL )
				{
					continue;
				}

				aptr = (HaAtom*) impr_dihedral.pt1;
				aptr->FillRef(buf);
				os << buf << "  ";

				aptr = (HaAtom*) impr_dihedral.pt2;
				aptr->FillRef(buf);
				os << buf << "  ";

				aptr = (HaAtom*) impr_dihedral.pt3;
				aptr->FillRef(buf);
				os << buf << "  ";

				aptr = (HaAtom*) impr_dihedral.pt4;
				aptr->FillRef(buf);
				os << buf << "  " << std::endl;
			}
			os << "END IMPROPER ANGLES" << std::endl;
		}

		int nv = 0;
		set<MMBond, less<MMBond> >::iterator mbitr = p_mm_model->MBonds.begin();

		for(; mbitr != p_mm_model->MBonds.end(); mbitr++)
		{
			MMBond& bnd = (MMBond&) *mbitr;
			if( bnd.set_type == MolMechModel::SET_SPEC)
			{
				if( nv == 0) os << "SPECIAL VALENCE BONDS" << std::endl;
				nv++;
				
				aptr = (HaAtom*) bnd.pt1;
				aptr->FillRef(buf);
				os << buf << "  ";

				aptr = (HaAtom*) bnd.pt2;
				aptr->FillRef(buf);
				os << buf << "  ";
				sprintf(buf,"%12.6f  %12.6f ",bnd.r0,bnd.fc);
				os << buf << std::endl;
			}
		}
		if( nv > 0) os << "END SPECIAL VALENCE BONDS" << std::endl;

		nv = 0;

		set<MMValAngle, less<MMValAngle> >::iterator vaitr = p_mm_model->ValAngles.begin();

		for(; vaitr != p_mm_model->ValAngles.end(); vaitr++)
		{
			MMValAngle& vang = (MMValAngle&) *vaitr;
			if( vang.set_type == MolMechModel::SET_SPEC)
			{
				if( nv == 0) os << "SPECIAL VALENCE ANGLES" << std::endl;
				nv++;
				 
				aptr = (HaAtom*) vang.pt1;
				aptr->FillRef(buf);
				os << buf << "  ";

				aptr = (HaAtom*) vang.pt2;
				aptr->FillRef(buf);
				os << buf << "  ";

				aptr = (HaAtom*) vang.pt3;
				aptr->FillRef(buf);
				os << buf << "  ";

				sprintf(buf,"%12.6f  %12.6f ",vang.a0,vang.fc);
				os << buf << std::endl;
			}  
		}
		if( nv > 0) os << "END SPECIAL VALENCE ANGLES" << std::endl;
		os << "#END MOLECULAR MECHANICS MODULE " << std::endl;
	}

	os << "#END MOL_SET" << std::endl;
	return TRUE;
}

int MolSet::SaveHarlemFile(const char* filename, const AtomSaveOptions* popt_arg )
{
	std::auto_ptr<AtomSaveOptions> popt_auto( popt_arg == NULL ? (AtomSaveOptions*) p_save_opt_default->clone() : (AtomSaveOptions*) popt_arg->clone() );
	AtomSaveOptions* popt = popt_auto.get();

	ofstream fout(filename);
	if( fout.fail())
	{
		PrintLog(" Error in MolSet::SaveHarlemFile()  opening file %s\n",filename);
		return FALSE;
	}	
	popt->SetSaveHeader(true);
	popt->SetSaveFooter(true);
	
	int ires = SaveXMLToStream(fout, popt);

	return ires;
}


int MolSet::SaveHINFile(const char* filename )
{
	ofstream fout(filename);
	if( fout.fail())
	{
		PrintLog(" Error in MolSet::SaveHINFile()  opening file %s\n",filename);
		return FALSE;
	}	
	
	int ires = SaveHINToStream(fout);
	return ires;
}

int MolSet::SaveOldHarlemFile(const char* filename )
{
	ofstream fout(filename);
	if( fout.fail())
	{
		PrintLog(" Error in MolSet::SaveOldHarlemFile()   opening file %s\n",filename);
		return FALSE;
	}
	int ires = SaveOldHarlemStream(fout);
	return ires;
}

bool MolSet::DeleteMol(HaMolecule* pMol)
{
	return( DeleteAtoms(*pMol) );
}

double MolSet::FindClosestContact(HaAtom* atc1,HaAtom* atc2)
{
	double dist_min = 999999.0;
	atc1 = NULL;
	atc2 = NULL;
	AtomIteratorMolSet aitr1(this);
	AtomIteratorMolSet aitr2(this);
	HaAtom* aptr1;
	HaAtom* aptr2;

	for(aptr1 = aitr1.GetFirstAtom(); aptr1; aptr1 = aitr1.GetNextAtom())
	{
		for(aptr2 = aitr2.GetFirstAtom(); aptr2; aptr2 = aitr2.GetNextAtom())
		{
			if(aptr1 == aptr2) continue;
			double dist = Vec3D::CalcDistance(aptr1,aptr2);
			if(dist < dist_min)
			{
				dist_min = dist;
				atc1 = aptr1;
				atc2 = aptr2;
			}
		}
	}
	if(atc1 != NULL && atc2 != NULL)
	{
		char buf1[256],buf2[256];
		atc1->FillRef(buf1);
		atc2->FillRef(buf2);
		PrintLog(" Minimal atom-atom distance is %12.6f (Ang) between Atoms %s and %s \n",
			       dist_min, buf1,buf2);
	}
	return dist_min;
}

int MolSet::ExecuteCommand(CmdParser& cmd_pr)
{
	char buf[256];
	int option;

	cmd_pr.ResetCursorPosition();

    if( !cmd_pr.FetchToken() ) return FALSE;

	HaMolView* pView=NULL;
	pView= GetActiveMolView();
	
    switch( cmd_pr.CurToken )
    {			
	case(SelectTok): 
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken )
		{   
			option = NormAtomFlag;
			if( pView && pView->HetaGroups ) option |= HeteroFlag;
			if( pView && pView->Hydrogens )  option |= HydrogenFlag;
			SelectAtomsMask(option);
		} 
		else if( cmd_pr.CurToken==AllTok )
		{   
			SelectAtomsAll();
			DisplaySelectCount();
		} 
		else if( cmd_pr.CurToken==NoneTok )
		{   
			UnSelectAtomsAll();
			DisplaySelectCount();
		} 
		else
		{
			AtomExpr* p_expr;
			if( (p_expr = cmd_pr.ParseExpression(0,this)) != NULL )
			{   
				if( !cmd_pr.CurToken )
				{   
					SelectAtomsExprObj(p_expr);
				} 
				else 
				{
					PrintLog("Invalid command syntax\n");
				}
				delete p_expr;
			}
		}
		break;

	case(RestrictTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken )
		{   
			option = NormAtomFlag;
			if( pView && pView->HetaGroups ) 
				option |= HeteroFlag;
			if( pView && pView->Hydrogens )  
				option |= HydrogenFlag;
			SelectAtomsMask(option);
			if(pView) pView->RestrictSelected();
			if(pView) pView->ReDrawFlag |= RFRefresh;
		} 
		else if( cmd_pr.CurToken==AllTok )
        {   
			SelectAtomsAll();
			pView->RestrictSelected();
            pView->ReDrawFlag |= RFRefresh;
        } 
		else if( cmd_pr.CurToken==NoneTok )
        {   
			UnSelectAtomsAll();
			pView->RestrictSelected();
            pView->ReDrawFlag |= RFRefresh;
        } 
		else
		{
			AtomExpr* p_expr;
			if( (p_expr = cmd_pr.ParseExpression(0,this)) != NULL )
			{   
				if( !cmd_pr.CurToken )
				{   
					SelectAtomsExprObj(p_expr);
					pView->RestrictSelected();
					pView->ReDrawFlag |= RFRefresh;
				} 
				else 
				{
					PrintLog("Invalid command syntax\n");
				}
				delete p_expr;
			} 
		}
		break;
	case(OverlapMolTok):
		try
		{
			cmd_pr.FetchToken(); // this worked
			if( cmd_pr.CurToken != IdentTok ) throw std::runtime_error(" First argument should be an atom group name ");
				
			AtomGroup* group1;
			AtomGroup* group2;
			group1 = GetAtomGroupByID(cmd_pr.TokenIdent.c_str());
			if( group1 == NULL ) throw std::runtime_error(" First argument should be an atom group name ");

			cmd_pr.FetchToken();
			if( cmd_pr.CurToken != IdentTok ) throw std::runtime_error(" Second argument should be an atom group name ");
						
			group2 = GetAtomGroupByID(cmd_pr.TokenIdent.c_str());
			if( group2 == NULL ) throw std::runtime_error(" Second argument should be an atom group name ");
					
			OverlapMol(*group1,*group2); //modifies the coordinates of one molecule to overlap with the other
			RefreshAllViews(RFRefresh | RFApply); //plots the molecule
		}
		catch(std::exception& ex)
		{
			PrintLog("Invalid command syntax\n");
			PrintLog("%s\n",ex.what());
		}
		break;

	case(AlignOverlapMolTok):
		cmd_pr.FetchToken(); //this worked
		if(!cmd_pr.CurToken)
		{
			PrintLog("Invalid command syntax\n");
		}
		else
		{
			AtomExpr* p_expr;
			if( (p_expr=cmd_pr.ParseExpression(0,this)) != NULL) //read an expression from command line
			{
				AtomGroup firstatset(p_expr,this); //constructs atom set out of expression
				delete p_expr;
				if( (p_expr = cmd_pr.ParseExpression(0,this)) != NULL )
				{
					AtomGroup secatset(p_expr,this); // reads the second expression separated by space
					delete p_expr;
					if( !cmd_pr.CurToken )
					{
						if( secatset.size() != 0)
						{
							HaAtom* aptr = secatset[0];
							AlignOverlapMol(firstatset,aptr->GetHostMol()); //modifies the coordinates of one molecule to overlap with the other
							RefreshAllViews(RFRefresh | RFApply); //plots the molecules
						}
					}
					else
					{
						PrintLog("Invalid command syntax\n");
					}
				}		
			}
			else
			{
				PrintLog("Invalid command syntax\n");
			}
		}
		break;

    case(DefineTok):  
		cmd_pr.FetchToken();
		if( !(cmd_pr.CurToken == IdentTok)) 
		{   
			PrintLog(" Missing or invalid (reserved word) atom group name" );
			break;
		}
		{
			std::string group_name = cmd_pr.TokenIdent;

			AtomGroup* pgroup;
			pgroup = GetAtomGroupByID( group_name.c_str());

			if(pgroup == NULL) 
			{
				pgroup = AddAtomGroup(group_name.c_str());
			}
				
			if( cmd_pr.FetchToken() )
			{   
				AtomExpr* p_expr;
				if( (p_expr = cmd_pr.ParseExpression(0,this)) != NULL )
				{   
					pgroup->SetFromExpr(p_expr,this);
					delete p_expr;
				} 
			}
			if(pgroup)
			{
				sprintf(buf," Atom group %s defined with %d atoms",
					group_name.c_str(),pgroup->GetNAtoms());
				PrintMessage(buf);
			}
		}
		break;

	case(RefreshTok):
		RefreshAllViews();
		break;
		
	case(ShowTok):
		ExecuteShowCommand(cmd_pr);
		break;

	}
		
	if(pView && pView->ReDrawFlag)
		RefreshAllViews();

	return True;
}

int MolSet::ExecuteShowCommand(CmdParser& cmd_pr)
{
    double temp;
	
	if(GetNMol() < 1 )
		return False;

	HaMolecule* pMol = HostMolecules[0];

    switch( cmd_pr.FetchToken() )
    {
	case(InfoTok):
		pMol->DescribeMolecule();
		break;
		
	case(SequenceTok):
		pMol->DescribeSequence();
		break;
		
	case(SymmetryTok):
		PrintLog("\n");
		
		if( !per_bc->spacegroup.empty() )
		{
			PrintLog("Space Group ...... %s\n", per_bc->spacegroup.c_str());
						
			PrintLog("Unit cell A ...... %g\n", per_bc->GetA());
			PrintLog("Unit cell B ...... %g\n", per_bc->GetB());
			PrintLog("Unit cell C ...... %g\n", per_bc->GetC());
			
			temp = RAD_TO_DEG*per_bc->GetAlpha();
			PrintLog("Unit cell alpha .. %g\n", temp);
			
			temp = RAD_TO_DEG*per_bc->GetBeta();
			PrintLog("Unit cell beta ... %g\n", temp);
			
			temp = RAD_TO_DEG*per_bc->GetGamma();
			PrintLog("Unit cell gamma .. %g\n", temp);
		}
		else
			PrintLog("No crystal symmetry data!\n");
		PrintLog("\n");
		break;
		
	default:
		PrintLog("Invalid command argument\n");
    }
	return True;
}

bool MolSet::DeleteAtomWithRef(const char* atref)
{
	HaAtom* aptr = GetAtomByRef(atref);
	if(aptr == NULL) return false;

	return(DeleteAtom(aptr));

}

bool MolSet::DeleteAtoms(AtomContainer& atcoll)
{
	int nmol;

	AtomIteratorGen aitr(&atcoll);
	HaAtom* aptr;

    PtrSet del_atoms;

	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if( aptr->GetHostMolSet() == this && !del_atoms.IsMember(aptr)) 
		{
			del_atoms.insert(aptr);
		}
	}

	if( del_atoms.size() > 0 ) OnAtomSeqChange();

	vector<HaCompMod*>::iterator mitr;
	for(mitr = CompModules.begin(); mitr != CompModules.end(); mitr++)
	{
		(*mitr)->OnDelAtoms(atcoll);
	}

	p_zmat->OnDelAtoms(atcoll);
	
	std::vector<CrdSnapshot*>::iterator snap_itr;
	for(snap_itr = crd_snapshots.begin(); snap_itr != crd_snapshots.end(); )
	{
		CrdSnapshot* psnap = (*snap_itr);
		if( !psnap->IsValid() ) 
		{
			PrintLog(" Deleting invalid snapshot %s \n",psnap->GetName().c_str());
			delete psnap;
			snap_itr = crd_snapshots.erase(snap_itr);
			continue;
		}
		psnap->OnDelAtoms(atcoll);
		if( psnap->GetCrd().size() == 0 )
		{
			delete psnap;
			snap_itr = crd_snapshots.erase(snap_itr);
			continue;
		}
		snap_itr++;
	}

	nmol = GetNMol();
	vector<HaMolecule*>::iterator mol_itr;

	std::vector<HaBond*>::iterator bitr1 = Bonds.begin();
	std::vector<HaBond*>::iterator bitr2 = Bonds.begin();
	
	for( ; bitr1 != Bonds.end(); ++bitr1 )
	{
		 HaBond* bptr = *bitr1;
		 if( del_atoms.IsMember( bptr->GetFirstAtom() ) || del_atoms.IsMember( bptr->GetSecondAtom() ) ) 
		 {
			 delete bptr;
			 continue;
		 }
		 (*bitr2) = (*bitr1);
		 ++bitr2;
	}

	if( bitr2 != Bonds.end() ) Bonds.erase(bitr2,Bonds.end());

	set<HaHBond, less<HaHBond> >::iterator hbitr,hbitr2;
	for( hbitr= HBonds.begin(); hbitr != HBonds.end(); )
	{
		HaAtom* at1 = (*hbitr).src;
		HaAtom* at2 = (*hbitr).dst;
		if( del_atoms.IsMember(at1) || del_atoms.IsMember(at2) )
		{
			hbitr2 = hbitr; hbitr++;
			HBonds.erase(hbitr2);
		}
		else
		{
			hbitr++;
		}
	}

	std::vector<HaBond*>::iterator bbitr;
	for( bbitr = BackboneBonds.begin(); bbitr != BackboneBonds.end(); )
	{
		HaAtom* at1 = (*bbitr)->srcatom;
		HaAtom* at2 = (*bbitr)->dstatom;
		if( del_atoms.IsMember(at1) || del_atoms.IsMember(at2) )
		{
			HaBond* bbptr = (*bbitr);
			delete bbptr;
			bbitr = BackboneBonds.erase(bbitr);
		}
		else
		{
			bbitr++;
		}
	}
	
	for(mol_itr = HostMolecules.begin(); mol_itr != HostMolecules.end(); )
	{
		HaMolecule* pMol = (*mol_itr);

		list<HaChain>::iterator chain_itr;
		multimap<int, HaResidue*, less<int> >::iterator res_itr, res_itr2;
		
		for( chain_itr = pMol->Chains.begin(); chain_itr != pMol->Chains.end(); )
		{
			std::set<HaResidue*> del_res;
            for( res_itr = (*chain_itr).res_map.begin(); res_itr != (*chain_itr).res_map.end();)
			{
                HaResidue* pres = (*res_itr).second;
				pres->DeleteAtoms( del_atoms );
				if( (*res_itr).second->size() == 0)
				{
                    res_itr2 = res_itr;
                    res_itr++;
					del_res.insert(pres);
					delete pres;
					(*chain_itr).res_map.erase(res_itr2);
				}
				else
				{
					res_itr++;
				}
			}
			auto ritr_arr = (*chain_itr).res_arr.begin();
			for(; ritr_arr != (*chain_itr).res_arr.end();)
			{
				HaResidue* pres = *ritr_arr;
				if (del_res.count(*ritr_arr) > 0)
				{
					ritr_arr = (*chain_itr).res_arr.erase(ritr_arr);
				}
				else
				{
					ritr_arr++;
				}
			}
			if( (*chain_itr).GetNRes() == 0 )
			{
				chain_itr = pMol->Chains.erase(chain_itr);
			}
			else
			{
			    chain_itr++;
			}
		}
		
		if(pMol->GetNAtoms() == 0)
		{
			list<Object3D*>::iterator oitr;
			for(oitr = ViewObjects.begin(); oitr != ViewObjects.end(); oitr++)
			{
				if( (*oitr) == pMol)
				{
					delete pMol;
					ViewObjects.erase(oitr);
					break;
				}
			}
			mol_itr = HostMolecules.erase(mol_itr);
		}
		else
		{
			mol_itr++;
		}
	}
	
	list<ChemGroup>::iterator cg_itr;
	for( cg_itr = ChemGroups.begin(); cg_itr != ChemGroups.end(); )
	{
		(*cg_itr).DeleteAtoms( del_atoms );
		if( (*cg_itr).size() == 0)
		{
			cg_itr = ChemGroups.erase(cg_itr);
		}
		else
		{
			cg_itr++;
		}
	}

	list<AtomGroup>::iterator atg_itr;
	for( atg_itr = NamedAtomGroups.begin(); atg_itr != NamedAtomGroups.end(); )
	{
		(*atg_itr).DeleteAtoms( del_atoms );
		if( (*atg_itr).size() == 0)
		{
			atg_itr = NamedAtomGroups.erase(atg_itr);
		}
		else
		{
			atg_itr++;
		}
	}  

	PtrSet::iterator aitr_d;
	for(aitr_d = del_atoms.begin(); aitr_d != del_atoms.end(); aitr_d++)
	{
		delete (HaAtom*)(*aitr_d);
	}
	return true;
}


bool MolSet::DeleteAtom(HaAtom* aptr)
{
	if(aptr == NULL) return false;
	if( aptr->GetHostMolSet() != this)
		return false;

	AtomGroup atlist;
	atlist.InsertAtom(aptr);
	bool bres = DeleteAtoms(atlist);
	return bres;
}

void MolSet::OnAtomSeqChange()
{
	at_seq_num_map.clear();
	at_array.clear();
}

void MolSet::OnChangePeriodicity()
{
	HaMolMechMod* p_mm_mod = GetMolMechMod();
	if( p_mm_mod ) p_mm_mod->OnChangePeriodicity();
}

PointIterator* MolSet::GetPointIteratorPtr()
{
	return GetAtomIteratorPtr();
}

PointIterator_const* MolSet::GetPointIteratorPtr() const
{
	return GetAtomIteratorPtr();
}

int MolSet::GetNumPt() const
{
	return GetNAtoms();
}

AtomIterator* MolSet::GetAtomIteratorPtr()
{
	AtomIteratorMolSet* p_aitr = new AtomIteratorMolSet(this);
	return p_aitr;
}

AtomIterator_const* MolSet::GetAtomIteratorPtr() const
{
	AtomIteratorMolSet_const* p_aitr = new AtomIteratorMolSet_const(this);
	return p_aitr;
}

int MolSet::IsMember(const HaAtom* aptr) const 
{
	AtomIteratorMolSet_const aitr(this);
	const HaAtom* aptr1 = aitr.GetFirstAtom();
	for(; aptr1; aptr1 = aitr.GetNextAtom())
	{
		if( aptr1 == aptr ) return TRUE;
	}
	return FALSE;
}

int MolSet::GetNMol() const
{
	return(HostMolecules.size());
}


int MolSet::GetNAtoms() const
{
	int nat=0;
	for(int i=0; i < HostMolecules.size(); i++)
	{
		nat+= HostMolecules[i]->GetNAtoms();
	}
	return(nat);
}

int MolSet::GetNDumAtoms() const
{
    AtomIteratorMolSet_const aitr(this);
	const HaAtom* aptr;
	int n_dummy_atoms = 0;
	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
       if( aptr->IsDummy() ) n_dummy_atoms++;
	}
	return n_dummy_atoms;
}

int MolSet::GetNBonds() const
{
    return(Bonds.size());
}

int MolSet::GetNHBonds() const
{
    return(HBonds.size());
}

int MolSet::GetNSSBonds() const
{
	int nss = 0;

	std::vector<HaBond*>::const_iterator bitr;
	
	for( bitr = Bonds.begin() ; bitr != Bonds.end(); bitr++ )
	{
		const HaBond* bptr = *bitr;
		int elem1 = bptr->GetFirstAtom()->GetElemNo();
		int elem2 = bptr->GetSecondAtom()->GetElemNo();
		if( elem1 == 16 && elem2 == 16 ) nss++;
	}
    return( nss );
}

int MolSet::GetNBackbBonds() const
{
    return(BackboneBonds.size());
}

int MolSet::GetNRes() const
{
	int nres=0;
	for(int i=0; i < HostMolecules.size(); i++)
	{
		nres+= HostMolecules[i]->GetNRes();
	}
	return(nres);
}

int
MolSet::GetNChains() const
{
	int nchains=0;
	for(int i=0; i < HostMolecules.size(); i++)
	{
		nchains+= HostMolecules[i]->GetNChains();
	}
	return(nchains);
}


HaMolView* MolSet::GetActiveMolView()
{
	return mset_pview;
}

const HaMolView* MolSet::GetActiveMolView() const
{
	return mset_pview;
}


void MolSet::SelectAtomsAll()
{  
	AtomIteratorMolSet aitr(this);
	HaAtom* aptr;

	for( aptr= aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		aptr->Select();
	}
}

void MolSet::SelectAtoms(AtomContainer* pat_cont)
{
	if(pat_cont == NULL) 
		return;
	AtomIteratorMolSet aitr(this);
	HaAtom* aptr;
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		aptr->UnSelect();
	}
	
	AtomIteratorGen aitr_s(pat_cont);
	for( aptr = aitr_s.GetFirstAtom(); aptr; aptr = aitr_s.GetNextAtom())
	{
		aptr->Select();
	}
}


void MolSet::UnSelectAtomsAll()
{
	AtomIteratorMolSet aitr(this);
	HaAtom* aptr;

	for( aptr= aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		aptr->UnSelect();
	}
}

void MolSet::RevertAtomSelection()
{
	AtomIteratorMolSet aitr(this);
	HaAtom* aptr;

	for( aptr= aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if( aptr->Selected()) 
		{
			aptr->UnSelect();
		}
		else
		{
			aptr->Select();
		}
	}
}

int MolSet::AreHBonded(HaAtom* src, HaAtom* dst) const
{
	HaHBond hbond(src,dst);
	
	set<HaHBond, less<HaHBond> >::const_iterator hbitr= HBonds.find(hbond);
	if(hbitr != HBonds.end()) return TRUE;

	HaHBond hbond2(dst,src);
	hbitr= HBonds.find(hbond2);
	if(hbitr != HBonds.end()) return TRUE;

	return FALSE;
}

HaResidue* MolSet::GetResByRef(const char* res_ref)
{
	if(HostMolecules.empty())
		return NULL;

	std::string res_ref_str(res_ref);

	HaMolecule* pMol= HostMolecules[0];
	int ibeg_res_ref_in_mol=0;
	if(res_ref_str[0] == '$')
	{
		int iend_mol_name = res_ref_str.find('$',1);
		if( iend_mol_name != -1)
		{
			std::string mol_name= res_ref_str.substr(1,iend_mol_name-1);
			pMol= GetMolByName(mol_name.c_str());
			ibeg_res_ref_in_mol = iend_mol_name+1;
		}
	}
	if(pMol != NULL)
	{
		return(pMol->GetResByRef( res_ref_str.substr(ibeg_res_ref_in_mol).c_str() ) );
	}
	return NULL;	
}

HaAtom* MolSet::GetAtomByRef(const std::string& at_ref_str)
{
	if(HostMolecules.empty()) return NULL;

	HaMolecule* pMol= HostMolecules[0];
	int ibeg_at_ref_in_mol=0;
	if(at_ref_str[0] == '$')
	{
		int iend_mol_name = at_ref_str.find('$',1);
		if( iend_mol_name != -1)
		{
			std::string mol_name= at_ref_str.substr(1,iend_mol_name-1);
			pMol= GetMolByName(mol_name.c_str());
			ibeg_at_ref_in_mol = iend_mol_name+1;
		}
	}
	if(pMol != NULL)
	{
		return(pMol->GetAtomByRef( at_ref_str.substr(ibeg_at_ref_in_mol).c_str() ) );
	}
	return NULL;	
}
	

AtomIntMap MolSet::GetAtomSeqNumMap()
{
	at_seq_num_map.clear();

	HaAtom* aptr;
	AtomIteratorMolSet aitr(this);

	int i = 0;

	for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom())
	{
		at_seq_num_map[aptr] = i;
		i++;
	}
	return at_seq_num_map;
}

CAtomIntMap MolSet::GetAtomSeqNumMap() const
{
	CAtomIntMap at_seq_num_map_loc;

	const HaAtom* aptr;
	AtomIteratorMolSet_const aitr(this);
	int i = 0;

	for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom())
	{
		at_seq_num_map_loc[aptr] = i;
		i++;
	}
	return at_seq_num_map_loc;
}


bool MolSet::GetAtomsByRef(const char* at_ref, AtomGroup& at_set)
{
	at_set.clear();
	if(HostMolecules.empty())
		return false;

	HaMolecule* pMol= HostMolecules[0];

	std::string at_ref_str(at_ref);

	HaAtom* aptr;

	int ibeg_at_ref_in_mol=0;
	if(at_ref_str[0] == '$')
	{
		int iend_mol_ref = at_ref_str.find('$',1);
		if( iend_mol_ref != -1)
		{
			std::string mol_ref= at_ref_str.substr(1,iend_mol_ref-1);
			pMol= GetMolByRef(mol_ref.c_str());
			ibeg_at_ref_in_mol = iend_mol_ref+1;
		}
	}
	else if(at_ref_str.length() > 8 && at_ref_str.substr(0,8) == "chemgrp=")
	{
		std::string chem_grp_id = at_ref_str.substr(8);
		ChemGroup* chem_grp = GetChemGroupByID(chem_grp_id);
		if(chem_grp != NULL)
		{
			AtomIteratorAtomGroup aitr(chem_grp);
			for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
			{
                at_set.InsertAtom(aptr);
			}
		}
		return true;
	}

	if(ibeg_at_ref_in_mol >= at_ref_str.length() || at_ref_str[ibeg_at_ref_in_mol] == '*')
	{
		AtomIteratorMolecule aitr(pMol);
		for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
		{
			at_set.InsertAtom(aptr);
		}
		return true;
	}
	
	int idx_pt = at_ref_str.find_first_of('.',ibeg_at_ref_in_mol);
	if(idx_pt != -1 )
	{
		aptr = pMol->GetAtomByRef( at_ref_str.substr(ibeg_at_ref_in_mol).c_str() );
		if( aptr == NULL)
			return false;
		else
		{
			at_set.InsertAtom(aptr);
			return true;
		}
	}
	else
	{
		HaResidue* pres= pMol->GetResByRef(	at_ref_str.substr(ibeg_at_ref_in_mol) );
		if( pres == NULL )
			return false;
		else
		{
			AtomIteratorResidue aitr2(pres);
			aptr;
			for(aptr = aitr2.GetFirstAtom();aptr; aptr = aitr2.GetNextAtom())
			{
				at_set.InsertAtom(aptr);
			}
			return true;
		}	

	}
	return false;
}

HaAtom* MolSet::GetAtomBySeqNum(int seq_num)
{
	if( at_array.empty() ) InitAtomIdx();
	return at_array[seq_num];
}

int MolSet::GetSeqNumForAtom( HaAtom* aptr)
{
	if( at_seq_num_map.empty() ) InitAtomIdx();
	AtomIntMap::const_iterator mitr = at_seq_num_map.find(aptr);
	if( mitr == at_seq_num_map.end() ) return -1;
	return (*mitr).second;
}

void MolSet::InitAtomIdx()
{
	at_seq_num_map.clear();
	at_array.clear();

	AtomIteratorMolSet aitr(this);
	HaAtom* aptr;
	int idx = 0;
	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
	{
		at_array.push_back(aptr);
		at_seq_num_map[aptr]=idx;
	}
}

HaMolecule* MolSet::GetMolByName(const char* mol_name)
{
	MoleculesType::iterator mol_itr;
	for( mol_itr=HostMolecules.begin(); mol_itr != HostMolecules.end(); mol_itr++)
	{
		if( !strcmp((*mol_itr)->GetObjName(), mol_name ) )
		{
			return (*mol_itr);
		}
	}
	return NULL;
}

const HaMolecule* MolSet::GetMolByName(const char* mol_name) const
{
	MoleculesType::const_iterator mol_itr;
	for( mol_itr=HostMolecules.begin(); mol_itr != HostMolecules.end(); mol_itr++)
	{
		if( !strcmp((*mol_itr)->GetObjName(), mol_name ) )
		{
			return (*mol_itr);
		}
	}
	return NULL;
}

const HaMolecule* MolSet::GetMolByRef(const char* mol_ref_par ) const
{
	std::string mol_ref = mol_ref_par;
	bool has_serno = false;
	int ser_no = -1;
	size_t br_beg = mol_ref.find('[');
	size_t br_end = mol_ref.find(']');
	if (br_beg != std::string::npos && br_end != std::string::npos && br_beg < br_end)
	{
		std::string serno_str = mol_ref.substr(br_beg + 1, br_end - br_beg - 1);
		if (harlem::IsInt(serno_str))
		{
			has_serno = true;
			ser_no = std::stoi(serno_str);
		}
	}
	std::string mol_name = mol_ref;
	if (br_beg != std::string::npos) mol_name = mol_ref.substr(br_beg);

	if (!has_serno)
	{
		if (harlem::IsInt(mol_name))
		{
			ser_no = std::stoi(mol_name);
			mol_name = "";
		}
	}
	const HaMolecule* pMol = NULL;
	if (has_serno)
	{
		if (serno_mol_map.count(ser_no) > 0)
		{
			bool name_match = false;
			for (auto itr = serno_mol_map.find(ser_no); itr != serno_mol_map.end(); itr++)
			{
				pMol = itr->second;
				if (mol_name.empty() || mol_name == pMol->GetName())
				{
					name_match = true;
					break;
				}
				if (!name_match)
				{
					PrintLog("Warning in __FUNCTION__ :\n");
					PrintLog("Molecule name is not matched  in molecule Reference %s \n", mol_ref_par );
				}
			}
		}
	}
	else
	{
		if (name_mol_map.count(mol_name) > 0)
		{
			pMol = name_mol_map.find(mol_name)->second;
		}
	}
	return pMol;
}

HaMolecule* MolSet::GetMolByRef(const char* mol_ref) 
{ 
	return (HaMolecule*) ((const MolSet*)(this))->GetMolByRef( mol_ref);
}


const char* MolSet::GetName() const
{
	return name_mset.c_str();
}


ChemGroupIterator::ChemGroupIterator(MolSet* new_pmset)
{
	pmset = new_pmset;
	if(pmset == NULL) return;
    CurGroupItr = pmset->ChemGroups.end();
}
	

ChemGroupIterator::~ChemGroupIterator()
{

}

ChemGroup* ChemGroupIterator::GetFirst()
{
	if(!pmset || pmset->ChemGroups.empty()) return(NULL);
	CurGroupItr = pmset->ChemGroups.begin();
	return( &(*CurGroupItr) );
}

ChemGroup* ChemGroupIterator::GetNext()
{
	if(CurGroupItr == pmset->ChemGroups.end() ) return(NULL);
	CurGroupItr++;
	if(CurGroupItr == pmset->ChemGroups.end()) return(NULL);
	return( &(*CurGroupItr) );
}

AtomGroupIteratorMolSet::AtomGroupIteratorMolSet(MolSet* new_pmset)
{
	pmset = new_pmset;
	if(pmset == NULL) return;
    CurListItr = pmset->NamedAtomGroups.end();
}
	

AtomGroupIteratorMolSet::~AtomGroupIteratorMolSet()
{

}

AtomGroup* AtomGroupIteratorMolSet::GetFirst()
{
	if(!pmset || pmset->NamedAtomGroups.empty()) return(NULL);
	CurListItr = pmset->NamedAtomGroups.begin();
	return( &(*CurListItr) );
}

AtomGroup* AtomGroupIteratorMolSet::GetNext()
{
	if(CurListItr == pmset->NamedAtomGroups.end() ) return(NULL);
	CurListItr++;
	if(CurListItr == pmset->NamedAtomGroups.end()) return(NULL);
	return( &(*CurListItr) );
}

AtomGroupIteratorMolSet_const::AtomGroupIteratorMolSet_const(const MolSet* new_pmset)
{
	pmset = new_pmset;
	if(pmset == NULL) return;
    CurListItr = pmset->NamedAtomGroups.end();
}
	

AtomGroupIteratorMolSet_const::~AtomGroupIteratorMolSet_const()
{

}

const AtomGroup* AtomGroupIteratorMolSet_const::GetFirst() 
{
	if(!pmset || pmset->NamedAtomGroups.empty()) return(NULL);
	CurListItr = pmset->NamedAtomGroups.begin();
	return( &(*CurListItr) );
}

const AtomGroup* AtomGroupIteratorMolSet_const::GetNext() 
{
	if(CurListItr == pmset->NamedAtomGroups.end() ) return(NULL);
	CurListItr++;
	if(CurListItr == pmset->NamedAtomGroups.end()) return(NULL);
	return( &(*CurListItr) );
}


int MolSet::WrapAndCenter(const std::string& grp_name, const Vec3D& cnt_crd )  
{
	try
	{
		AtomGroup* p_grp = GetAtomGroupByID( grp_name.c_str() );
		if( p_grp == NULL ) throw std::runtime_error(" No Atom Group with name " + grp_name );

		HaAtom* aptr;
		AtomIteratorMolSet aitr(this);
		Vec3D cur_cnt;
		cur_cnt = p_grp->GetAverageCoord();
		Vec3D trans = cnt_crd - cur_cnt;
		TranslateAtoms(trans);
		WrapToUnitCell();
	}
	catch( const std::exception& ex )
	{
		PrintLog("Error in MolSet::WrapAndCenter() \n");
		PrintLog("%s\n",ex.what());
		return FALSE;
	}
	return TRUE;
}

int MolSet::WrapToUnitCell()
{
	if( !per_bc->IsSet() ) return FALSE; 
	MolEditor* p_mol_ed = GetMolEditor(true);
	return p_mol_ed->WrapToUnitCell(this,per_bc);
}

ResidueIteratorMolSet MolSet::GetResidueIterator()
{
	ResidueIteratorMolSet ritr(this);
	return ritr;
}

int MolSet::GetNChemGroups() const
{
	return(ChemGroups.size());
}
	

ChemGroup* MolSet::AddBlankChemGroup(const std::string& gid )
{
	std::string gid_loc=gid;
	if(gid_loc == "BUF")  // For Buffer region Atomic Groups
	{
		gid_loc=GetUniqChemGrpID(1);
	}
	else if(gid_loc == "" || !CheckUniqChemGrpID(gid_loc))
	{
		gid_loc=GetUniqChemGrpID(0);
	}
	ChemGroups.push_back(ChemGroup(this,gid_loc.c_str()));
    chemg_idx.clear();
	return &(ChemGroups.back());
}

bool MolSet::DeleteChemGroup(const std::string& gid )
{
	ChemGroupsType::iterator gitr;
	for(gitr=ChemGroups.begin(); gitr != ChemGroups.end(); gitr++)
	{
		if(gid == (*gitr).GetID())
		{
			ChemGroups.erase(gitr);
			chemg_idx.clear();
			return true;
		}
	}
	return false;
}

bool MolSet::DeleteChemGroupPtr( ChemGroup* grp_ptr  )
{
	if( grp_ptr == NULL) return false;
	ChemGroupsType::iterator gitr;
	for(gitr=ChemGroups.begin(); gitr != ChemGroups.end(); gitr++)
	{
		if(grp_ptr == &(*gitr))
		{
			ChemGroups.erase(gitr);
			chemg_idx.clear();
			return true;
		}
	}
	return false;
}

bool MolSet::DeleteAtomGroup(const char* id )
{
	NamedAtomGroupsType::iterator gitr;
	std::string gid = id;
	for(gitr= NamedAtomGroups.begin(); gitr != NamedAtomGroups.end(); gitr++)
	{
		if(gid == (*gitr).GetID())
		{
			NamedAtomGroups.erase(gitr);
			return true;
		}
	}
	return false;
}

bool MolSet::DeleteAtomGroupPtr( AtomGroup* alist_ptr  )
{
	if( alist_ptr == NULL) return false;
	NamedAtomGroupsType::iterator gitr;
	for(gitr= NamedAtomGroups.begin(); gitr != NamedAtomGroups.end(); gitr++)
	{
		if(alist_ptr == &(*gitr))
		{
			NamedAtomGroups.erase(gitr);
			return true;
		}
	}
	return false;
}

int MolSet::CreateAxxMol(const char* mol_name, const char* id)
{
	AtomGroup* pgrp = GetAtomGroupByID(id);
	if(pgrp == NULL)
	{
		PrintLog("Error in MolSet::CreateAxxMol() \n");
	    PrintLog("No atom group with the name %s \n",id);
		return FALSE;
	}

	HaMolecule* pmol = GetMolByName(mol_name);
	if(pmol != NULL)
	{
		DeleteMol(pmol);
	}
	
	pmol = AddNewMolecule();
	pmol->SetObjName(mol_name);

	AtomIteratorAtomGroup aitr(pgrp);
	HaAtom* aptr;
	HaResidue* pres_old = NULL;
	HaChain* pchain_old = NULL;

	HaChain* chain = NULL;
	HaResidue* pres = NULL;

	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		HaResidue* pres_ref   = aptr->GetHostRes();
		HaChain*   pchain_ref = aptr->GetHostChain();

		if(pchain_ref != pchain_old)
		{
			chain = pmol->AddChain(pchain_ref->ident); 
			pchain_old = pchain_ref;
		}

		if(pres_ref != pres_old)
		{
			int new_resno = pres_ref->GetSerNo();
			pres= chain->AddResidue(new_resno);
			if( pres == NULL)
			{
				int new_resno = chain->GetUniqResSerNo(false);
				pres = chain->AddResidue(new_resno);
			}
			std::string ref_name = pres_ref->GetName();
			std::string mod_name = ref_name;
			pres->SetName(mod_name.c_str());
			pres_old = pres_ref;
		}

		HaAtom* aptr_new  = pres->AddNewAtom();
		aptr_new->SetParamFrom(*aptr);
		aptr_new->SetX( aptr_new->GetX() + 0.000001); // shift atom position to avoid infinite reuplsion energies
		aptr_new->SetDummy();
	}
	return TRUE;
}


bool MolSet::SetChemGrpSelected(const std::string& gid)
{
	ChemGroup* gptr=GetChemGroupByID(gid);
	if(!gptr)
	{
		PrintLog(" Error in MolSet::SetChemGrpSelected \n");
		PrintLog(" there is no group %s in the molecular set \n", gptr->GetID());
		return false;
	}

	gptr->clear();

	HaAtom* aptr;
	
	AtomIteratorMolSet aitr(this);
	for(aptr= aitr.GetFirstAtom(); aptr; aptr=aitr.GetNextAtom())
	{
		if(aptr->Selected())
			gptr->InsertAtom(aptr);
	}
	return true;
}

bool MolSet::SetStdChemGroups()
{
	ChemGroups.clear();
	HaResidue* group;
	HaResidue *res_next;
	HaResidue *res_prev;
	HaChain* chain;
	ChemGroup* gcur_ptr;
	HaAtom* aptr;

	ResidueIteratorMolSet ritr(this);
	for(group = ritr.GetFirstRes(); group; group = ritr.GetNextRes())
	{
		// Peptide backbone
		cout << " Calc residue " << group->GetSerNo() << endl;
		if(group->IsProtein())
		{
			res_next=group->GetNextResInChain();
			bool result1=group->IsBonded(res_next);
			if( (!res_next) || (!res_next->IsProtein()) || (!result1) )
			{
				gcur_ptr=AddBlankChemGroup();
				gcur_ptr->InsertAtom(group->GetAtomByName("C"));
				gcur_ptr->InsertAtom(group->GetAtomByName("O"));
				if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
			}
			gcur_ptr=AddBlankChemGroup();
			aptr=group->GetAtomByName("N");
			gcur_ptr->InsertAtom(aptr);
			gcur_ptr->InsertAtom(group->GetAtomByName("H"));

			res_prev= group->GetPrevResInChain();
			bool result2=false;
			if(res_prev != NULL) result2=group->IsBonded(res_prev);
 			if(res_prev && res_prev->IsProtein()   && result2 )
			{
				gcur_ptr->InsertAtom(res_prev->GetAtomByName("C"));
				gcur_ptr->InsertAtom(res_prev->GetAtomByName("O"));
			}
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}

		if(group->IsTerm())
		{
			gcur_ptr=AddBlankChemGroup("BUF");
			vector<HaAtom*>::iterator aitr;
			for(aitr=group->begin(); aitr != group->end(); aitr++)
			{
				gcur_ptr->InsertAtom(*aitr);
			}
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
			continue;
		}


		if( !strcmp_trunc(group->GetName(),"ALA") )   // ALA
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}

		else if( !strcmp_trunc(group->GetName(),"GLY" ) ) // GLY
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}

		else if( !strcmp_trunc(group->GetName(),"LEU") ) // LEU
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CG"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HG"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom( group->GetAtomByName("CD1")  );
			gcur_ptr->InsertAtom( group->GetAtomByName("HD11") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HD12") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HD13") );
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom( group->GetAtomByName("CD2")  );
			gcur_ptr->InsertAtom( group->GetAtomByName("HD21") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HD22") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HD23") );
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}

		else if(!strcmp_trunc(group->GetName(),"SER" ) ) // SER
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom( group->GetAtomByName("OG") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HG") );
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}

		else if(!strcmp_trunc(group->GetName(),"VAL"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom( group->GetAtomByName("CG1")  );
			gcur_ptr->InsertAtom( group->GetAtomByName("HG11") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HG12") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HG13") );
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom( group->GetAtomByName("CG2")  );
			gcur_ptr->InsertAtom( group->GetAtomByName("HG21") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HG22") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HG23") );
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

		}

		else if(!strcmp_trunc(group->GetName(),"THR"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom( group->GetAtomByName("OG1")  );
			gcur_ptr->InsertAtom( group->GetAtomByName("HG1") );
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom( group->GetAtomByName("CG2")  );
			gcur_ptr->InsertAtom( group->GetAtomByName("HG21") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HG22") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HG23") );
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}

		else if(!strcmp_trunc(group->GetName(),"LYS"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom( group->GetAtomByName("CG")  );
			gcur_ptr->InsertAtom( group->GetAtomByName("HG2") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HG3") );
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom( group->GetAtomByName("CD")  );
			gcur_ptr->InsertAtom( group->GetAtomByName("HD2") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HD3") );
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom( group->GetAtomByName("CE")  );
			gcur_ptr->InsertAtom( group->GetAtomByName("HE2") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HE3") );
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom( group->GetAtomByName("NZ")  );
			gcur_ptr->InsertAtom( group->GetAtomByName("HZ1") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HZ2") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HZ3") );
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}

		else if(!strcmp_trunc(group->GetName(),"ASP"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom( group->GetAtomByName("CG")  );
			gcur_ptr->InsertAtom( group->GetAtomByName("OD1") );
			gcur_ptr->InsertAtom( group->GetAtomByName("OD2") );
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}

		else if(!strcmp_trunc(group->GetName(),"ILE"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom( group->GetAtomByName("CG1")  );
			gcur_ptr->InsertAtom( group->GetAtomByName("HG12") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HG13") );
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom( group->GetAtomByName("CG2")  );
			gcur_ptr->InsertAtom( group->GetAtomByName("HG21") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HG22") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HG23") );
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom( group->GetAtomByName("CD1")  );
			gcur_ptr->InsertAtom( group->GetAtomByName("HD11") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HD12") );
			gcur_ptr->InsertAtom( group->GetAtomByName("HD13") );
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}

		else if(!strcmp_trunc(group->GetName(),"ASN"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom( group->GetAtomByName("CG")  );
			gcur_ptr->InsertAtom( group->GetAtomByName("OD1") );
			gcur_ptr->InsertAtom( group->GetAtomByName("ND2"));
			gcur_ptr->InsertAtom( group->GetAtomByName("HD21"));
			gcur_ptr->InsertAtom( group->GetAtomByName("HD22"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}

		else if(!strcmp_trunc(group->GetName(),"GLU"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CG"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HG2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HG3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
			
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CD"));
			gcur_ptr->InsertAtom(group->GetAtomByName("OE1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("OE2"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}


		else if(!strcmp_trunc(group->GetName(),"PRO"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CG"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HG2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HG3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CD"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HD2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HD3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

		}

		else if(!strcmp_trunc(group->GetName(),"ARG"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CG"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HG2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HG3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CD"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HD2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HD3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("NE"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HE"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CZ"));
			gcur_ptr->InsertAtom(group->GetAtomByName("NH1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HH11"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HH12"));
			gcur_ptr->InsertAtom(group->GetAtomByName("NH2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HH21"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HH22"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}

		else if(!strcmp_trunc(group->GetName(),"PHE"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CG"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CD1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HD1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CD2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HD2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CE1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HE1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CE2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HE2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CZ"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HZ"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}

		else if(!strcmp_trunc(group->GetName(),"TYR"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CG"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CD1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HD1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CD2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HD2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CE1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HE1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CE2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HE2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CZ"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
			
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("OH"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HH"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}

		else if(!strcmp_trunc(group->GetName(),"GLN"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CG"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HG2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HG3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom( group->GetAtomByName("CD")  );
			gcur_ptr->InsertAtom( group->GetAtomByName("OE1") );
			gcur_ptr->InsertAtom( group->GetAtomByName("NE2"));
			gcur_ptr->InsertAtom( group->GetAtomByName("HE21"));
			gcur_ptr->InsertAtom( group->GetAtomByName("HE22"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}

		else if(group->IsHistidine() ) // HIS
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CG"));
			gcur_ptr->InsertAtom(group->GetAtomByName("ND1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HD1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CD2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HD2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CE1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HE1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("NE2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HE2"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}

		else if(group->IsCysteine() ) // CYS
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("SG"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HG"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}

		else if(!strcmp_trunc(group->GetName(),"MET"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CG"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HG2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HG3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("SD"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CE"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HE1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HE2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HE3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}
		else if(!strcmp_trunc(group->GetName(),"TRP"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HA"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CG"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CD1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HD1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CD2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("NE1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HE1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CE2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CE3"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HE3"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CZ2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HZ2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CZ3"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HZ3"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CH2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HH2"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}
		else if(!strcmp_trunc(group->GetName(),"HEM"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("FE"));
			gcur_ptr->InsertAtom(group->GetAtomByName("NA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("C1A"));
			gcur_ptr->InsertAtom(group->GetAtomByName("C2A"));
			gcur_ptr->InsertAtom(group->GetAtomByName("C3A"));
			gcur_ptr->InsertAtom(group->GetAtomByName("C4A"));
			gcur_ptr->InsertAtom(group->GetAtomByName("NB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("C1B"));
			gcur_ptr->InsertAtom(group->GetAtomByName("C2B"));
			gcur_ptr->InsertAtom(group->GetAtomByName("C3B"));
			gcur_ptr->InsertAtom(group->GetAtomByName("C4B"));
			gcur_ptr->InsertAtom(group->GetAtomByName("NC"));
			gcur_ptr->InsertAtom(group->GetAtomByName("C1C"));
			gcur_ptr->InsertAtom(group->GetAtomByName("C2C"));
			gcur_ptr->InsertAtom(group->GetAtomByName("C3C"));
			gcur_ptr->InsertAtom(group->GetAtomByName("C4C"));
			gcur_ptr->InsertAtom(group->GetAtomByName("ND"));
			gcur_ptr->InsertAtom(group->GetAtomByName("C1D"));
			gcur_ptr->InsertAtom(group->GetAtomByName("C2D"));
			gcur_ptr->InsertAtom(group->GetAtomByName("C3D"));
			gcur_ptr->InsertAtom(group->GetAtomByName("C4D"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CHA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HGM"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CHB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HDM"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CHC"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HAM"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CHD"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HBM"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CAA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HP71"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HP72"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CBA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HP73"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HP74"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CGA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("O1A"));
			gcur_ptr->InsertAtom(group->GetAtomByName("O2A"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CMA"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HM81"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HM82"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HM83"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CMB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HM11"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HM12"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HM13"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CAB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HV21"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CBB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HV22"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HV23"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
			
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CMC"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HM31"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HM32"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HM33"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CAC"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HV41"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CBC"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HV42"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HV43"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CMD"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HM51"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HM52"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HM53"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
			
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CAD"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HP61"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HP62"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CBD"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HP63"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HP64"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CGD"));
			gcur_ptr->InsertAtom(group->GetAtomByName("O1D"));
			gcur_ptr->InsertAtom(group->GetAtomByName("O2D"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		
		}
		else if(!strcmp_trunc(group->GetName(),"RBP"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("RU"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CG"));
			gcur_ptr->InsertAtom(group->GetAtomByName("ND1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HD1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CD2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HD2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CE1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HE1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("NE2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HE2"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CB"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HB3"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("NL1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CL2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CL3"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HL3"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CL4"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HL4"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CL5"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HL5"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CL6"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HL6"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("NM1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CM2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CM3"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HM3"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CM4"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HM4"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CM5"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HM5"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CM6"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HM6"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("NN1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CN2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CN3"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HN3"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CN4"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HN4"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CN5"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HN5"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CN6"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HN6"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("NP1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CP2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CP3"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HP3"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CP4"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HP4"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CP5"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HP5"));
			gcur_ptr->InsertAtom(group->GetAtomByName("CP6"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HP6"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}
		else if(!strcmp_trunc(group->GetName(),"RUA"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("RU"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("N1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("H11"));
			gcur_ptr->InsertAtom(group->GetAtomByName("H12"));
			gcur_ptr->InsertAtom(group->GetAtomByName("H13"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("N2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("H21"));
			gcur_ptr->InsertAtom(group->GetAtomByName("H22"));
			gcur_ptr->InsertAtom(group->GetAtomByName("H23"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("N3"));
			gcur_ptr->InsertAtom(group->GetAtomByName("H31"));
			gcur_ptr->InsertAtom(group->GetAtomByName("H32"));
			gcur_ptr->InsertAtom(group->GetAtomByName("H33"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("N4"));
			gcur_ptr->InsertAtom(group->GetAtomByName("H41"));
			gcur_ptr->InsertAtom(group->GetAtomByName("H42"));
			gcur_ptr->InsertAtom(group->GetAtomByName("H43"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());

			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("N5"));
			gcur_ptr->InsertAtom(group->GetAtomByName("H51"));
			gcur_ptr->InsertAtom(group->GetAtomByName("H52"));
			gcur_ptr->InsertAtom(group->GetAtomByName("H53"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}
		else if(!strcmp_trunc(group->GetName(),"FSF") ||
				!strcmp_trunc(group->GetName(),"FSA") ||
				!strcmp_trunc(group->GetName(),"FSB"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("FE1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("S1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("FE2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("S2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("FE3"));
			gcur_ptr->InsertAtom(group->GetAtomByName("S3"));
			gcur_ptr->InsertAtom(group->GetAtomByName("FE4"));
			gcur_ptr->InsertAtom(group->GetAtomByName("S4"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}
		else if(!strcmp_trunc(group->GetName(),"PCL") ||
				!strcmp_trunc(group->GetName(),"PCA") ||
				!strcmp_trunc(group->GetName(),"PCB"))
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("FE1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("FE2"));
			gcur_ptr->InsertAtom(group->GetAtomByName("FE3"));
			gcur_ptr->InsertAtom(group->GetAtomByName("FE4"));
			gcur_ptr->InsertAtom(group->GetAtomByName("FE5"));
			gcur_ptr->InsertAtom(group->GetAtomByName("FE6"));
			gcur_ptr->InsertAtom(group->GetAtomByName("FE7"));
			gcur_ptr->InsertAtom(group->GetAtomByName("FE8"));
			gcur_ptr->InsertAtom(group->GetAtomByName("S1A"));
			gcur_ptr->InsertAtom(group->GetAtomByName("S2A"));
			gcur_ptr->InsertAtom(group->GetAtomByName("S3A"));
			gcur_ptr->InsertAtom(group->GetAtomByName("S4A"));
			gcur_ptr->InsertAtom(group->GetAtomByName("S2B"));
			gcur_ptr->InsertAtom(group->GetAtomByName("S3B"));
			gcur_ptr->InsertAtom(group->GetAtomByName("S4B"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}
		else if( !strcmp_trunc(group->GetName(),"MTL" ) ) // Methylene
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("CT"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HT1"));
			gcur_ptr->InsertAtom(group->GetAtomByName("HT2"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}
		else if( !strcmp_trunc(group->GetName(),"HTM" ) ) // Methylene
		{
			gcur_ptr=AddBlankChemGroup();
			gcur_ptr->InsertAtom(group->GetAtomByName("HT"));
			if(gcur_ptr->size() == 0) DeleteChemGroup(gcur_ptr->GetID());
		}
	}
	return true;

}

bool MolSet::CheckUniqChemGrpID(const std::string& gid)
{
	ChemGroupsType::iterator itr;
	for(itr=ChemGroups.begin(); itr != ChemGroups.end(); itr++)
	{
		if(gid == (*itr).GetID()) return false;
	}
	return true;
}

static const int MIN_MAIN_GRP_ID=10;
static const int MAX_MAIN_GRP_ID=99999;
static const int MIN_BUF_GRP_ID=100000;
static const int MAX_BUF_GRP_ID=200000;

std::string MolSet::GetUniqChemGrpID(int buf_reg_flag)
// but_reg_flag - TRUE if ChemGroup belongs to buffer Region
{
	int low=MIN_MAIN_GRP_ID;
	int high=MAX_MAIN_GRP_ID;
	if(buf_reg_flag)
	{
		low=MIN_BUF_GRP_ID;
		high=MAX_BUF_GRP_ID;
	}
	int gidmax=low;
	ChemGroupsType::iterator itr;
	for(itr=ChemGroups.begin(); itr != ChemGroups.end(); itr++)
	{
		std::string cur_gid=(*itr).GetID();
		int igid = atoi(cur_gid.c_str());
		if( (igid > gidmax) && (igid <= high)) gidmax= igid;
	}
	gidmax++;
	if(gidmax > high)
	{
		PrintLog(" Error in MolSet::GetUniqChemGroupID() \n");
        PrintLog(" gidmax > high \n");
	}
	char buf[256];
	sprintf(buf,"%d",gidmax);
	return buf;
}

ChemGroup& MolSet::GetChemGroupByIdx(const int index)
{
	if(chemg_idx.size() == 0 )
	{
		int n = ChemGroups.size();
		chemg_idx.resize(n);
		ChemGroupsType::iterator gitr;
		int i = 0;
		for(gitr = ChemGroups.begin(); gitr != ChemGroups.end(); gitr++)
		{
			 chemg_idx[i] = &(*gitr);
			 i++;
		}
	}
	return(*((ChemGroup*)chemg_idx[index-1]));	
}

//const ChemGroup&
//MolSet::GetChemGroupByIdx(const int index) const
//{
//	return(ChemGroups[index-1]);	
//}


ChemGroup* MolSet::GetChemGroupByID(const std::string& gid)
{
	ChemGroupsType::iterator gitr;
	
	for(gitr=ChemGroups.begin(); gitr != ChemGroups.end(); gitr++)
	{
		if(gid == (*gitr).GetID() )
		{
			return( &(*gitr));
		}
	}
	return NULL;
}

ChemGroup*
MolSet::GetChemGroupByAtom(const HaAtom* aptr)
{
	ChemGroupsType::iterator gitr;
	
	for(gitr=ChemGroups.begin(); gitr != ChemGroups.end(); gitr++)
	{
		if((*gitr).IsMember(aptr) )
		{
			return( &(*gitr));
		}
	}
	return NULL;
}

AtomGroup* MolSet::SetAtomGroupFromSelection( const char* id)
{
	AtomGroup* atlist = GetAtomGroupByID(id);
	if( atlist == NULL) 
		atlist = AddAtomGroup(id);
	if(atlist == NULL) return atlist;
	atlist->clear();
	AtomIteratorMolSet aitr(this);
	HaAtom* aptr;
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if(aptr->Selected()) atlist->InsertAtom(aptr);
	}
	return atlist;
}

int MolSet::AssociateFragment(MolSet* frag)
{
	if(frag == NULL) return FALSE;

	int nf= Fragments.size();
	int i;
	for(i=0; i < nf; i++)
	{
		if( Fragments[i] == (void*)frag) return TRUE;	
	}

	Fragments.push_back(frag);

	AtomMapping* frag_atom_maping = new AtomMapping(this,frag);
	frag_atom_maps[frag] = frag_atom_maping;
	
	frag_atom_maping->Map2to1ByAtomDistance();
	
	return TRUE;
}

int MolSet::ReleaseFragment(MolSet* frag)
{
	if(frag == NULL) return FALSE;

	if( frag_atom_maps.count(frag) > 0)
	{
		AtomMapping* at_map_ptr = (AtomMapping*) frag_atom_maps[frag];
		delete at_map_ptr;
		frag_atom_maps.erase(frag);
	}
	
	vector<MolSet*>::iterator fitr;

	for(fitr = Fragments.begin(); fitr != Fragments.end(); fitr++)
	{
		if( (*fitr) == (void*)frag) 
		{
			Fragments.erase(fitr);
			return TRUE;
		}
	}

	return FALSE;
}

int MolSet::DeleteFragment(MolSet* frag)
{
	int ires = ReleaseFragment(frag);

	if( !ires) return FALSE;

	delete frag;
	return TRUE;
}

int MolSet::ReleaseAllFragments()
{
	PtrPtrMap_itr mitr = frag_atom_maps.begin();

	for(; mitr != frag_atom_maps.end(); mitr++ )
	{
		AtomMapping* at_map_ptr = (AtomMapping*) mitr.GetVal();
		delete at_map_ptr;
    }
	frag_atom_maps.clear();

	vector<MolSet*>::iterator fitr;

	for(fitr = Fragments.begin(); fitr != Fragments.end(); )
	{
		fitr = Fragments.erase(fitr);
	}

	return TRUE;	
}

int MolSet::DeleteAllFragments()
{
	PtrPtrMap::iterator mitr = frag_atom_maps.begin();
	PtrPtrMap::iterator mitr_end = frag_atom_maps.end();

	for(; mitr != mitr_end; mitr++ )
	{
		AtomMapping* at_map_ptr = (AtomMapping*) mitr.GetVal();
		delete at_map_ptr;
    }
	frag_atom_maps.clear();
	
	vector<MolSet*>::iterator fitr;
	for(fitr = Fragments.begin(); fitr != Fragments.end(); )
	{
		MolSet* frag = (MolSet*)(*fitr);
		delete frag;
		fitr = Fragments.erase(fitr);
	}

	return TRUE;	
}


int MolSet::BuildFragmentAtomMap(MolSet* frag, AtomAtomMap& frag_atom_map)
{
	bool debug = false;
	frag_atom_map.clear();

	HaAtom* aptr;
	HaAtom* aptr_f;
	AtomIteratorMolSet aitr_f(frag);

// Set Atom-Atom maps for atoms with the same references:
	for( aptr_f = aitr_f.GetFirstAtom(); aptr_f; aptr_f = aitr_f.GetNextAtom() )
	{
		std::string at_ref_f = aptr_f->GetRef();
		aptr = this->GetAtomByRef( at_ref_f );
		if( aptr != NULL )
		{
			if( debug ) PrintLog(" MolSet::BuildFragmentAtomMap for atom %s match is found ", at_ref_f.c_str() ); 
			frag_atom_map[aptr] = aptr_f;
			frag_atom_map[aptr_f] = aptr;
		}
		else
		{
			if( debug ) PrintLog(" MolSet::BuildFragmentAtomMap for atom %s match is not found ", at_ref_f.c_str() ); 
		}
	}

// Set Atom-Atom map for coinciding points

	BoxPartition part_table;
	double xmin,ymin,zmin,xmax,ymax,zmax;

	frag->GetMinMaxCrd(xmin,ymin,zmin,xmax,ymax,zmax);
	
	part_table.SetBoundaries(xmin-0.5,ymin-0.5,zmin-0.5,xmax+0.5,ymax+0.5,zmax+0.5);
	part_table.SetRegionRad(0.0);

	part_table.DistributePointsToCells(*frag);
	
	AtomIteratorMolSet aitr(this);

	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		int idx = part_table.GetPointCellIdx(aptr);
		if( idx >= 0)
		{
			int np = part_table[idx].size();
			int i;
			for( i = 0 ; i < np; i++)
			{
				aptr_f = (HaAtom*) part_table[idx][i];
				if( Vec3D::CalcDistance(aptr,aptr_f) < 0.002)
				{
					if( frag_atom_map.count(aptr) > 0 &&  frag_atom_map[aptr] != aptr_f ) 
					{
						PrintLog(" Warning in MolSet::BuildFragmentAtomMap() atom %s is already mapped to different atom remapping by distance " );
					}
					if( frag_atom_map.count(aptr_f) > 0 &&  frag_atom_map[aptr_f] != aptr ) 
					{
						PrintLog(" Warning in MolSet::BuildFragmentAtomMap() atom %s is already mapped to different atom remapping by distance " );
					}
					frag_atom_map[aptr] = aptr_f;
					frag_atom_map[aptr_f] = aptr;
				}	
			}
		}
	}
	return TRUE;
}

int MolSet::SelectAtomsMatchingFragment(MolSet* frag)
{
	int nf = frag_atom_maps.size();

	AtomMapping* at_map_ptr = (AtomMapping*) frag_atom_maps[frag];
	if( at_map_ptr->p_ac_2 == frag )
	{
		std::map< HaAtom*, HaAtom*, less<HaAtom*> >::iterator mitr;
		for( mitr = at_map_ptr->atmap_2to1.begin(); mitr != at_map_ptr->atmap_2to1.end(); mitr++)
		{
			HaAtom* p_ref_atom = (*mitr).second;
			p_ref_atom->Select();
		}
	}
	return TRUE;
}


int MolSet::IsFragment(const MolSet* pmset)
{
	if(pmset == NULL) return FALSE;
	int nf = Fragments.size();
	int i;
	for( i = 0; i < nf; i++)
	{
		if( Fragments[i] == (const void*) pmset ) return TRUE;
	}
	return FALSE;
}

int MolSet::FragmentIdx(const MolSet* pmset)
{
	if(pmset == NULL) return -1;
	int nf = Fragments.size();
	
	int i;
	for( i = 0; i < nf; i++)
	{
		if( Fragments[i] == (const void*) pmset ) 
		{
			return i;
		}
	}
	return -1;
}
	

MolSet* MolSet::CreateFragmentFromSelection(const char* frag_name, FragmentCreatePars* params )
{
	HaAtom *aptr, *faptr;
	HaMolecule* pMol;
	HaMolecule*	pfMol;
	HaChain *chain;
	HaResidue *group, *fres;
	HaBond* bptr;

	HaChain*   pch_cur;
	HaResidue* pres_cur;

	if( this->GetNChemGroups() == 0) this->SetStdChemGroups();

	MolSet* pfrag= new MolSet;
	pfrag->SetName( frag_name );
	
	int cur_serno=0;
	vector<HaAtom*>::iterator paitr;

	int mol_new=0;
	int ch_new=0;
	int res_new=0;

// First cycle on Atomic Groups and include into selected all atoms
// in partially selected groups
	
	ChemGroupsType::iterator gitr;
	for(gitr=ChemGroups.begin(); gitr != ChemGroups.end();gitr++)
	{
		int include=0;
		AtomIteratorAtomGroup aitr_chem_group(&(*gitr));
		for( aptr= aitr_chem_group.GetFirstAtom(); aptr; aptr= aitr_chem_group.GetNextAtom())
		{
			if(aptr->Selected())
			{
				include=1;
				break;
			}
		}
		if(!include) continue;

		for( aptr= aitr_chem_group.GetFirstAtom(); aptr; aptr= aitr_chem_group.GetNextAtom())
		{
			aptr->Select();
		}
	}

// cycle on all heavy atoms and include all hydrogen bound to these atoms:

	AtomIteratorMolSet aitr(this);
	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if( !aptr->Selected() ) continue;
		if( aptr->IsHydrogen() ) continue;
		AtomGroup bonded_atoms;
		aptr->GetBondedAtoms( bonded_atoms );
		AtomIteratorAtomGroup aitr_b( &bonded_atoms );
		HaAtom* aptr_b;
		for( aptr_b = aitr_b.GetFirstAtom(); aptr_b; aptr_b = aitr_b.GetNextAtom())
		{
			if( aptr_b->IsHydrogen()) aptr_b->Select();
		}
	}

	AtomAtomMap frag_at_map;
	AtomAtomMap::iterator mitr;

// Copy Selected atoms from parent molecule to the fragment

	MoleculesType::iterator mol_itr;
	for( mol_itr=HostMolecules.begin(); mol_itr != HostMolecules.end(); mol_itr++)
	{
		mol_new=1;
		ChainIteratorMolecule chitr(*mol_itr);
		for( chain= chitr.GetFirstChain();chain; chain= chitr.GetNextChain())
		{
			ch_new=1;
			ResidueIteratorChain ritr_ch(chain);
			for(group = ritr_ch.GetFirstRes();group;group = ritr_ch.GetNextRes())
			{
				res_new=1;
				for(paitr=group->begin(); paitr != group->end(); paitr++)
				{
					if((*paitr)->Selected())
					{
						if(mol_new)
						{
							pfMol=pfrag->AddNewMolecule();
							pfMol->SetObjName((*mol_itr)->GetObjName());
							mol_new = 0;
						}
						if(ch_new)
						{
							pch_cur = pfMol->AddChain(chain->ident);
							pch_cur->SetParamFrom(*chain);
							ch_new = 0;
						}
						if(res_new)
						{
							fres = pch_cur->AddResidue(group->GetSerNo());
							fres->SetParamFrom(*group);
							res_new = 0;
						}
						faptr = fres->AddNewAtom();
						faptr->SetParamFrom(*(*paitr));
						frag_at_map[*paitr] = faptr ;
					}
				}
			}
		}
	}

	AssociateFragment(pfrag);

	HaAtom *faptr1, *faptr2;
	HaBond *fbptr;

	// Copy Atomic Groups:

	for(gitr=ChemGroups.begin(); gitr != ChemGroups.end();gitr++)
	{
		AtomIteratorAtomGroup aitr_chem_group(&(*gitr));
		aptr = aitr_chem_group.GetFirstAtom();
		if(aptr == NULL || !aptr->Selected()) continue;

		ChemGroup* fpgrp= pfrag->AddBlankChemGroup((*gitr).GetID());
		fpgrp->SetProtect((*gitr).GetProtect());

		for( aptr= aitr_chem_group.GetFirstAtom(); aptr; aptr= aitr_chem_group.GetNextAtom())
		{
			pMol=aptr->GetHostMol();
			std::string mol_name = pMol->GetObjName();
			pfMol= pfrag->GetMolByName(mol_name.c_str());
			if( pMol == NULL || pfMol == NULL) 
			{
				PrintLog("MolSet::CreateFragementFromSelection() pMol == NULL || pfMol == NULL \n");
				continue;
			}
			mitr = frag_at_map.find(aptr);
			if(mitr != frag_at_map.end() )
			faptr1= (*mitr).second;
			fpgrp->InsertAtom(faptr1);
		}
	}

	int igid_trm=MIN_BUF_GRP_ID+1;

	BondIteratorMolSet bitr(this);
	for(bptr=bitr.GetFirstBond(); bptr; bptr = bitr.GetNextBond())
	{
		if( bptr->srcatom->Selected() && bptr->dstatom->Selected())
		// if both bonded atoms belong to the selected region just
		// include this bond into the fragment
		{
			pMol  = bptr->srcatom->GetHostMol();
			pfMol = pfrag->GetMolByName(pMol->GetObjName());
			if( pfMol == NULL ) 
			{
				PrintLog("Error in MolSet::CreateFragmenFromSelection() \n");
				continue;
			}
			
			mitr = frag_at_map.find(bptr->srcatom);
			if(mitr == frag_at_map.end())
				continue;
			faptr1= (*mitr).second;
			mitr = frag_at_map.find(bptr->dstatom);
			if(mitr == frag_at_map.end())
				continue;
			faptr2=(*mitr).second;
			fbptr=pfrag->AddBond(faptr1,faptr2);
			fbptr->SetParamFrom(*bptr);
			continue;
		}

		if( bptr->srcatom->Selected() || bptr->dstatom->Selected())
		// if only one of the bonded atoms belongs to the fragment
		// generate truncated hydrogen atom instead of another one
		{
			HaAtom* aptr1=NULL;
			HaAtom* aptr2=NULL;
			if(bptr->srcatom->Selected())
			{
				aptr1=bptr->srcatom;
				aptr2=bptr->dstatom;
			}
			else
			{
				aptr1=bptr->dstatom;
				aptr2=bptr->srcatom;
			}
			pMol  = bptr->srcatom->GetHostMol();
			pfMol = pfrag->GetMolByName(pMol->GetObjName());
			
			mitr = frag_at_map.find(aptr1);
			if(mitr == frag_at_map.end() )
				continue;

			faptr1 = (*mitr).second;
			pch_cur = faptr1->GetHostChain();
		
			int res_serno= pch_cur->GetUniqResSerNo(1);
			HaResidue* pfres = pch_cur->AddResidue(res_serno);
			pfres->SetName("TRM");

			faptr2 = pfres->AddNewAtom();
         
			char buf[256];
	        sprintf(buf,"%d",igid_trm);

			ChemGroup* fpgrp= pfrag->AddBlankChemGroup(buf);
			igid_trm++;
		    fpgrp->SetProtect(1.0);
			fpgrp->InsertAtom(faptr2);

			HaAtom::SetCoordSubstH(aptr1,aptr2,faptr2);
			pfrag->AddBond(faptr1,faptr2);
		}	
	}

	// Set QChem Parameters: Basis Set and Local Active Orb Basis Set:
	HaQCMod* ptr_qc_mod= GetQCMod(false);

	if(ptr_qc_mod != NULL)
	{
		HaQCMod* ptr_qc_mod_frag= pfrag->GetQCMod(true);

		const std::string bname(ptr_qc_mod->GetBasName());
		if(!bname.empty())
		{
			ptr_qc_mod_frag->InitBasis( bname.c_str() );
		}
		const std::string lo_set_id(ptr_qc_mod->GetLocOrbSetID());
		if(!lo_set_id.empty())
		{
			ptr_qc_mod_frag->InitLocOrb(lo_set_id.c_str());
		}
	}

	return pfrag;
}

int MolSet::SyncFragmentCoord(MolSet* frag)
{
	if(!this->IsFragment(frag) )
	{
		PrintLog("Error in MolSet::SyncFragmentCoord() \n");
        PrintLog("Argument is not a fragment of the current Molecular Set \n");
		return FALSE;
	}
	if( frag_atom_maps.count(frag) == 0 )
	{
		PrintLog("Error in MolSet::SyncFragmentCoord() \n");
		PrintLog("AtomMapping is not set for the fragment \n");
		return FALSE;
	}

	AtomMapping* at_map_ptr = (AtomMapping*) frag_atom_maps[frag];
	int ires = at_map_ptr->SyncAtomCrd2From1();

	return ires;
}

bool MolSet::CalcDipole()
{
	int na=GetNAtoms();
	
	const HaAtom* aptr;
	double chrg=0.0;
	HaVec_double dipm(3,0.0);
	HaVec_double center_pt(3,0.0);

	int nn=0;
	
	AtomIteratorMolSet aitr(this);
	for(aptr= aitr.GetFirstAtom(); aptr; aptr=aitr.GetNextAtom())
	{
		if(aptr != NULL && aptr->Selected())
		{
			double ach = aptr->GetCharge();
			chrg+=ach;
			nn++;
			dipm(1) += ach*(aptr->GetX());
			dipm(2) += ach*(aptr->GetY());
			dipm(3) += ach*(aptr->GetZ());
				
			center_pt(1)+= aptr->GetX();
			center_pt(2)+= aptr->GetY();
			center_pt(3)+= aptr->GetZ();
		}		
	}
	dipm(1) -= chrg* ( center_pt(1)/nn );
	dipm(2) -= chrg* ( center_pt(2)/nn );
	dipm(3) -= chrg* ( center_pt(3)/nn );

	PrintLog(" Charge of the selected molecule(s) as a sum of partial atomic charges is: \n");
	PrintLog("%16.9f e \n\n",chrg);
	PrintLog(" Dipole moment of the molecule(s) from atomic partial charges is: \n");
	PrintLog(" dx= %16.9f  dy=%16.9f dz= %16.9f   (in e*Ang ) \n ",dipm(1),dipm(2),dipm(3));

	return true;
}

bool MolSet::SetVdwRadii()
{
	HaAtom* aptr;

	AtomIteratorMolSet aitr(this);

	for(aptr= aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		int elno = aptr->GetElemNo();
		aptr->radius = HaAtom::ElemVDWRadius(elno);
	}
	return true;
}

bool MolSet::SetParseRadii()
{
	HaAtom* aptr;

	AtomIteratorMolSet aitr(this);

	for(aptr= aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		int elno = aptr->GetElemNo();
		if( elno == 1 )
			aptr->radius = 1.0;
		else if( elno == 6 )
			aptr->radius = 1.7;
		else if( elno == 7 )
			aptr->radius = 1.5;
		else if( elno == 8 )
			aptr->radius = 1.6;
		else if( elno == 11 ) // Na
			aptr->radius = 1.85;
		else if( elno == 15 ) // P
			aptr->radius = 2.0;
		else if( elno == 16 ) // S
			aptr->radius = 1.9;
		else if( elno == 19 ) // K
			aptr->radius = 2.17;
		else if( elno == 26 ) // Fe
			aptr->radius = 2.0;
		else if( elno == 44 ) // Ru 
			aptr->radius = 2.0;
		else
			aptr->radius = HaAtom::ElemVDWRadius(elno);
	}
	return true;
}

bool MolSet::SetHPPRadii()
{
	HaAtom* aptr;

	AtomIteratorMolSet aitr(this);

	for(aptr= aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		int elno = aptr->GetElemNo();
		if( elno == 1 )
			aptr->radius = 1.2;  // from h++
		else if( elno == 6 )
			aptr->radius = 1.7;  // from h++
		else if( elno == 7 )
			aptr->radius = 1.55; // from h++
		else if( elno == 8 )
			aptr->radius = 1.5;  // from h++
		else if( elno == 11 ) // Na
			aptr->radius = 1.85;  // PARSE
		else if( elno == 15 ) // P
			aptr->radius = 2.0;   // PARSE
		else if( elno == 16 ) // S 
			aptr->radius = 1.8;   // from h++
		else if( elno == 19 ) // K
			aptr->radius = 2.17;  // PARSE
		else if( elno == 26 ) // Fe
			aptr->radius = 2.0;   // PARSE
		else if( elno == 44 ) // Ru 
			aptr->radius = 2.0;   // PARSE
		else
			aptr->radius = HaAtom::ElemVDWRadius(elno);
	}
	return true;
}

double MolSet::CalculatePotential( double x, double y, double z )
{
    HaAtom  *aptr;
    double dx,dy,dz;
    double result;
    double dist;
    double max;

    /* 8.0 Angstrom Cut Off */
    max = 8.0*8.0;

    result = 0.0;
   
	AtomIteratorMolSet aitr(this);

    for(aptr=aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
    {
		dx = aptr->GetX()-x;
        if( (dist=dx*dx) < max )
        {
			dy = aptr->GetY() - y;
            if( (dist+=dy*dy) < max )
            {
				dz = aptr->GetZ() - z;
                if( (dist+=dz*dz) < max )
                    result += aptr->GetCharge() / sqrt(dist);
            }
        }
    }

	/* Dielectric constant 10.0 */
    result = (HARTREE_TO_KT)*result/10.0;
    return( result );
}

AtomDoubleMap* MolSet::GetChargeMapByName(const char* map_name)
{
	int nn= ChargeMaps.size();
	int i;
	for(i = 0; i < nn; i++)
	{
		if( stricmp_loc(ChargeMaps[i].GetName(),map_name) == 0 )
		{
			return &ChargeMaps[i];
		}
	}
	return (AtomDoubleMap*) NULL;
}

AtomDoubleMap* MolSet::CreateChargeMap(const char* map_name)
{
	ChargeMaps.push_back(AtomDoubleMap(map_name));
	return &ChargeMaps.back();

}

int 
MolSet::SetChargeMapByCurrentCharges(const char* map_name)
{
	AtomDoubleMap* charge_map = GetChargeMapByName(map_name);
    if(charge_map == NULL )
	{
		charge_map = CreateChargeMap(map_name);
	}
	charge_map->clear();

	AtomIteratorMolSet aitr(this);
	HaAtom* aptr;
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
	{
		charge_map->insert(AtomDoubleMap::value_type(aptr,aptr->GetCharge()) );
	}
	return True;
}

int 
MolSet::SetChargesFromChargeMap(AtomDoubleMap* charge_map)
{
	if(charge_map == NULL) 
		return False;

	AtomIteratorMolSet aitr(this);
	HaAtom* aptr;
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
	{
		aptr->SetCharge(charge_map->GetValue(aptr));
	}

	return True;
	
}

bool
MolSet::Print_info(ostream &sout, const int level)
{
	MoleculesType::iterator mol_itr;
	for( mol_itr=HostMolecules.begin(); mol_itr != HostMolecules.end(); mol_itr++)
	{
		(*mol_itr)->Print_info(sout,level);
	}

	HaQCMod* ptr_qc_mod = GetQCMod(false);
	if(ptr_qc_mod != NULL)
		ptr_qc_mod->Print_info(sout, level);
	return true;
}

HaCompMod* MolSet::GetCompModule(int mtype, bool create_module )
{
	vector<HaCompMod*>::iterator mitr;
	for(mitr = CompModules.begin(); mitr != CompModules.end(); mitr++)
	{
		if( (*mitr) == NULL) continue;
		if( (*mitr)->GetType() == mtype)
			return (*mitr);
	}
	HaCompMod* pmod = NULL;;
	if(create_module)
	{
		pmod= HaCompMod::CreateCompMod(mtype,this);
		CompModules.push_back(pmod);
	}
	return pmod;
}

const HaCompMod* MolSet::GetCompModule( int mtype ) const  
{
	vector<HaCompMod*>::const_iterator mitr;
	for(mitr = CompModules.begin(); mitr != CompModules.end(); mitr++)
	{
		if( (*mitr) == NULL) continue;
		if( (*mitr)->GetType() == mtype) return (*mitr);
	}
	return NULL;
}

ETCouplMod* MolSet::GetETCouplMod( const bool create_module )
{
	return( (ETCouplMod*) GetCompModule( COMP_MOD_ET_COUPL ,create_module) );
}

HaQCMod* MolSet::GetQCMod( const bool create_module )
{
	return( (HaQCMod*) GetCompModule( COMP_MOD_QCHEM ,create_module) );
}

const HaQCMod* MolSet::GetQCMod() const
{
	return( (const HaQCMod*) GetCompModule( COMP_MOD_QCHEM ) );
}



HaGaussMod* MolSet::GetGaussMod( const bool create_module )
{
	return( (HaGaussMod*) GetCompModule( COMP_MOD_GAUSSIAN ,create_module) );
}

HaDaltonMod* MolSet::GetDaltonMod( const bool create_module )
{
	return( (HaDaltonMod*) GetCompModule( COMP_MOD_DALTON ,create_module) );
}

ElectrostMod* MolSet::GetElectrostMod( const bool create_module )
{
	if( ElectrostMod::ActiveElectrMod == COMP_MOD_ELECTROST)
	{
		return( (ElectrostMod*) GetCompModule( COMP_MOD_ELECTROST ,create_module) );
	}
	else
	{
		printf("Will Use ElMod instead of ElectrostMod\n");
		return( (ElectrostMod*) GetCompModule( COMP_MOD_EL ,create_module) );
	}
}
#ifdef ELMOD_COMPILE
ElMod* MolSet::GetElMod( const bool create_module )
{
	return( (ElMod*) GetCompModule( COMP_MOD_EL ,create_module) );
}
#endif
pKaCalcMod* MolSet::GetpKaCalcMod( const bool create_module )
{
	pKaCalcMod* m_pKaCalcMod=(pKaCalcMod*) GetCompModule( COMP_MOD_PKA_CALC ,create_module);
	if(m_pKaCalcMod==NULL)
		printf("Cannot create pKaCalcMod\n");
	return m_pKaCalcMod;
}

PNPMod* MolSet::GetPNPMod( const bool create_module )
{
	return( (PNPMod*) GetCompModule( COMP_MOD_PNP ,create_module) );
}
APBSMod* MolSet::GetAPBSMod( const bool create_module )
{
  return( (APBSMod*) GetCompModule( COMP_MOD_APBS ,create_module) );
}
HaInterMolMod* MolSet::GetInterMolMod( const bool create_module )
{
	return( (HaInterMolMod*) GetCompModule( COMP_MOD_INTERMOL ,create_module) );
}

HaMolMechMod* MolSet::GetMolMechMod( const bool create_module )
{
	return( (HaMolMechMod*) GetCompModule( COMP_MOD_MOLMECH ,create_module) );
}

const HaMolMechMod* MolSet::GetMolMechMod() const 
{
	return( (const HaMolMechMod*) GetCompModule( COMP_MOD_MOLMECH ) );
}

MDTrajAnalMod*  MolSet::GetTrajAnalMod( bool create_module)
{
	HaMolMechMod* p_mm_mod = this->GetMolMechMod( create_module );
	if( p_mm_mod == NULL ) return NULL;
	return p_mm_mod->GetTrajAnalMod();
}

HaScatterMod* MolSet::GetScatterMod( const bool create_module )
{
	HaScatterMod* pMod = (HaScatterMod*) GetCompModule( COMP_MOD_SCATTER ,create_module );
	return( pMod);
}

StmMod* MolSet::GetSTMMod( const bool create_module )
{
	StmMod* pMod = (StmMod*) GetCompModule( COMP_MOD_STM ,create_module );
	return( pMod);
}

NuclAcidMod* MolSet::GetNuclAcidMod( const bool create_module )
{
	NuclAcidMod* pMod = (NuclAcidMod*) GetCompModule( COMP_MOD_NUCL_ACID ,create_module );
	return( pMod);
}

HaZindoMod* MolSet::GetZindoMod( const bool create_module )
{
	return( (HaZindoMod*) GetCompModule( COMP_MOD_ZINDO ,create_module) );
}

ProtonRedoxMod* MolSet::GetProtonRedoxMod( const bool create_module )
{
	return( (ProtonRedoxMod*) GetCompModule( COMP_MOD_PROTON_REDOX ,create_module) );
}

HaEmpiricalMod* MolSet::GetEmpiricalMod( const bool create_module )
{
	return( (HaEmpiricalMod*) GetCompModule( COMP_MOD_EMPIRICAL ,create_module) );
}

HaMolMembraneMod* MolSet::GetMolMembraneMod(const bool create_module) // jose
{
	return( (HaMolMembraneMod*) GetCompModule( COMP_MOD_MEMBRANE , create_module) ); 
}

HaFlexMod* MolSet::GetFlexMod(const bool create_module) 
{
	return( (HaFlexMod*) GetCompModule( COMP_MOD_FLEX , create_module) ); 
}

CollectCrdAnalMod*  MolSet::GetCollectCrdAnalMod(const bool create_module )
{
	return( (CollectCrdAnalMod*) MolSet::GetCompModule( COMP_MOD_CLUSTER_ANAL , create_module) ); 
}

MolEditor* MolSet::GetMolEditor(const bool create_module) 
{
	if( p_mol_editor != NULL) return p_mol_editor;
	if( !create_module ) return NULL;
	p_mol_editor = new MolEditor();
	return( p_mol_editor ); 
}


void MolSet::RenumberGrp()
{
	PrintLog(" MolSet::OnRenumberGrp() \n"); 
	ChemGroup* gptr;
	ChemGroupIterator gitr(this);
	int igid=10;
	for(gptr= gitr.GetFirst(); gptr; gptr= gitr.GetNext())
	{
	    char buf[256];
	    sprintf(buf,"%d",igid);
		gptr->id = buf;
		igid++;
	}
}



bool MolSet::AddObject3D(Object3D* new_view_object)
{
	if(new_view_object == NULL)
		return false;
	ViewObjects.push_back(new_view_object);
	return true;
}


bool MolSet::DeleteObject3D(Object3D* pobj)
{
	if(pobj != NULL)
		return( DeleteObject3D(pobj->GetObjName()) );
	else 
		return false;
}


bool MolSet::DeleteObject3D(const std::string& obj_name)
{
	list<Object3D*>::iterator oitr;
	for(oitr = ViewObjects.begin(); oitr != ViewObjects.end(); oitr++)
	{
		if( !strcmp((*oitr)->GetObjName(), obj_name.c_str()) )
		{
			if((*oitr)->GetObjType() == OBJ3D_MOLECULE)
			{
				MoleculesType::iterator mol_itr;
				for( mol_itr=HostMolecules.begin(); mol_itr != HostMolecules.end(); mol_itr++)
				{
					HaMolecule* pMol = (HaMolecule*) (*mol_itr);
					if( !strcmp( pMol->GetObjName(),obj_name.c_str() ) )
					{
						DeleteAtoms(*pMol);
						return true;
					}
				}
			}
			else
			{
				delete(*oitr);
				ViewObjects.erase(oitr);
			}
			return true;
		}
	}
	return false;
}

AtomGroup* MolSet::AddAtomGroup(const char* id)
{
	NamedAtomGroups.push_back(AtomGroup());
	AtomGroup& at_arr= NamedAtomGroups.back();
    at_arr.SetID(id);
	return &at_arr;
}

AtomGroup* MolSet::GetAtomGroupByID( const char* id)
{
	NamedAtomGroupsType::iterator gitr;
	for( gitr = NamedAtomGroups.begin(); gitr != NamedAtomGroups.end(); gitr++)
	{
		if( stricmp_loc((*gitr).GetID(),id) == 0 ) return &(*gitr);
	}
	return NULL;
}

HaDisplayedSurface* MolSet::CalcMolSurfDens()
{
	HaField3D mol_dens;

	double MinX, MinY, MinZ, MaxX, MaxY, MaxZ;
	GetMinMaxCrd(MinX, MinY, MinZ, MaxX, MaxY, MaxZ);
	
	mol_dens.SetGridCornersCoord( MinX - 10.0, MinY - 10.0, MinZ - 10.0, 
		                          MaxX + 10.0, MaxY + 10.0, MaxZ + 10.0);

	int nx, ny, nz;
	double grid_space = 0.5;

	nx = (int)((mol_dens.GetXmax() - mol_dens.GetXmin())/grid_space) + 1;
	ny = (int)((mol_dens.GetYmax() - mol_dens.GetYmin())/grid_space) + 1;
	nz = (int)((mol_dens.GetZmax() - mol_dens.GetZmin())/grid_space) + 1;

	mol_dens.SetDimensions(nx, ny, nz);

	double stepx = mol_dens.GetXstep();
	double stepy = mol_dens.GetYstep();
	double stepz = mol_dens.GetZstep();

	mol_dens.SetDimensions(nx, ny, nz);
	mol_dens.FillZeros();

	double solv_rad = 3.0;

	AtomIteratorMolSet aitr(this);
	HaAtom* aptr;

	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
	{
		double xa = aptr->GetX();
		double ya = aptr->GetY();
		double za = aptr->GetZ();
		double frad = aptr->radius + solv_rad;
		double frad2 = frad*frad;
		
		int mrad = (int)(frad/grid_space) + 2;
		
		int ix, iy, iz;
		
		mol_dens.GetClosestGridPoint(xa,ya,za, ix,iy,iz);
		
		int lx = MaxFun(0, ix - mrad); 
		int mx = MinFun(nx-1, ix + mrad);
		int ly = MaxFun(0, iy - mrad); 
		int my = MinFun(ny-1, iy + mrad);
		int lz = MaxFun(0, iz - mrad); 
		int mz = MinFun(nz-1, iz + mrad);
		
		for( ix = lx ; ix <= nx ; ix++ )
		{
			for( iy = ly ; iy <= ny ; iy++ )
			{
				for( iz = lz ; iz <= nz ; iz++ )
				{
					float xp, yp, zp;
					mol_dens.GetXYZ(xp, yp , zp, ix, iy, iz);
					double dist = ( xp - xa)*(xp - xa);
					dist+= ( yp - ya)*(yp - ya);
					dist+= ( zp - za)*(zp - za);
					if( dist < frad2)
					{
						float* fld_ptr = mol_dens.GetValPtr(ix,iy,iz);
						(*fld_ptr) = 1.0; ;
					}
				}
			}
		}
	}
	
	HaDisplayedSurface* sptr;

// Positive Potential Isosurface
	sptr = new HaDisplayedSurface();
	if(sptr == NULL)
		return NULL;

	AddObject3D(sptr);
	
	bool result;
	result = sptr->calc_isosurf(&mol_dens, 0.5);

	sptr->SetObjName("VdW_surface");

	return sptr;
}

HaDisplayedSurface* MolSet::GetMolSurface(int create_flag)
{
	HaDisplayedSurface* sptr = NULL;
	list<Object3D*>::iterator oitr;
	for(oitr = this->ViewObjects.begin(); oitr != this->ViewObjects.end(); oitr++)
	{
		if( (*oitr)->GetObjType() == OBJ3D_SURFACE )
		{
			std::string obj_name = (*oitr)->GetObjName();
			if( obj_name == "VdW_surface") sptr = (HaDisplayedSurface*) (*oitr);
		}
	}
	if( sptr != NULL) return sptr;
	if(!create_flag)  return sptr;

	sptr = new HaDisplayedSurface();
	if(sptr == NULL) return NULL;

	sptr->SetObjName("VdW_surface");
	AddObject3D(sptr);
	return sptr;
}

HaDisplayedSurface* MolSet::CalcMolSurface(int surf_type)
{
	HaDisplayedSurface* sptr = GetMolSurface(TRUE);
	if(sptr == NULL) return sptr;

	AtomGroup at_arr;
	HaAtom* aptr;
	
	AtomIteratorMolSet aitr(this);
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if(aptr->Selected()) at_arr.push_back(aptr);
	}

	int result;
	double solv_rad = 1.4;

	result = HaSurface::CalcMolSurf(sptr, surf_type, solv_rad, at_arr);

	return sptr;
}

int MolSet::CalcSolventAccessArea()
{
	AtomGroup at_arr;
	HaAtom* aptr;
	
	AtomIteratorMolSet aitr(this);
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if(aptr->Selected()) at_arr.push_back(aptr);
	}

	int nat = at_arr.size();

	HaMat_double coord(3,nat);
	HaVec_double rad(nat);
	
	int i;
	for(i = 0; i < nat; i++)
	{
		coord.r0(0,i) = at_arr[i]->GetX();
		coord.r0(1,i) = at_arr[i]->GetY();
		coord.r0(2,i) = at_arr[i]->GetZ();
		rad[i] =  at_arr[i]->radius;
	}

	HaSurface surf;

	double solv_rad = 1.4;

	int ires = surf.CalcMolSurfAlpha(FALSE,solv_rad,coord,rad);

	for(i = 0; i < nat; i++)
	{
		at_arr[i]->solv_access_area = surf.surface_alpha[i];
	}

	return TRUE;
}

int MolSet::SaveCrdExclVolArb()
{
	HaMolecule* pmol = this->GetMolByName("EXCLUDED_VOLUME");
	if( pmol == NULL )
	{
		PrintLog("No EXCLUDED_VOLUME Molecule ");
		return FALSE;
	}
	if( pmol->GetNAtoms() == 0 )
	{
		PrintLog("EXCLUDED_VOLUME Molecule Is Empty");
		return FALSE;
	}

	char buf[128];
	std::string fname = "excluded_vol.xml";

	std::ofstream os(fname.c_str());
	os << "<\?xml version=\"1.0\" encoding=\"utf-8\"\?>" << std::endl;
	os << "<ExcludedVolume>" << std::endl;
	if( os.fail() )
	{
		PrintLog("Failed to Write to %s \n",fname.c_str());
		return FALSE;
	}

	AtomIteratorMolecule aitr(pmol);
	HaAtom* aptr;
	for( aptr = aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		sprintf(buf,"<Point X=\"%7.2f\" Y=\"%7.2f\" Z=\"%7.2f\"/>",aptr->GetX_Ang(), aptr->GetY_Ang(),aptr->GetZ_Ang());
		os << buf << std::endl;
	}	
	os << "</ExcludedVolume>" << std::endl;

	return TRUE;
}


int MolSet::CreateExcludedVolumeMol()
{
	char buf[128];
	HaAtom* aptr;
	AtomGroup selected_atoms;
	AtomIteratorMolSet aitr(this);
	
	for( aptr = aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		selected_atoms.InsertAtom(aptr);
	}
	if( selected_atoms.empty() ) 
	{
		PrintLog("No Atoms Selected() \n");
		return FALSE;
	}
	
	double MinX, MinY, MinZ, MaxX, MaxY, MaxZ;
	selected_atoms.GetMinMaxCrd(MinX, MinY, MinZ, MaxX, MaxY, MaxZ);
	HaField3D mol_dens;

	double grid_space = 0.5;

	MinX = MinX - 5.0;
	MinY = MinY - 5.0;
	MinZ = MinZ - 5.0;

	MinX = ((int)(MinX/grid_space) - 1)*grid_space;
	MinY = ((int)(MinY/grid_space) - 1)*grid_space;
	MinZ = ((int)(MinZ/grid_space) - 1)*grid_space;

	MaxX = MaxX + 5.0;
	MaxY = MaxY + 5.0;
	MaxZ = MaxZ + 5.0;

	int nx, ny, nz;
	int ix, iy, iz;

	nx = (int)((MaxX - MinX)/grid_space) + 1;
	ny = (int)((MaxY - MinY)/grid_space) + 1;
	nz = (int)((MaxZ - MinZ)/grid_space) + 1;

	MaxX = MinX + grid_space*(nx-1);
	MaxY = MinY + grid_space*(ny-1);
	MaxZ = MinZ + grid_space*(nz-1);

	mol_dens.SetGridCornersCoord( MinX, MinY, MinZ, MaxX, MaxY, MaxZ);
	mol_dens.SetDimensions(nx, ny, nz);

	double stepx = mol_dens.GetXstep();
	double stepy = mol_dens.GetYstep();
	double stepz = mol_dens.GetZstep();

	mol_dens.FillZeros();

	double solv_rad = 1.51;
	int mrad;
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
	{
		double xa = aptr->GetX();
		double ya = aptr->GetY();
		double za = aptr->GetZ();
		double frad = aptr->radius + solv_rad;
		double frad2 = frad*frad;
		
		mrad = (int)(frad/grid_space) + 2;
		
		mol_dens.GetClosestGridPoint(xa,ya,za, ix,iy,iz);
		
		int lx = MaxFun(0, ix - mrad); 
		int mx = MinFun(nx-1, ix + mrad);
		int ly = MaxFun(0, iy - mrad); 
		int my = MinFun(ny-1, iy + mrad);
		int lz = MaxFun(0, iz - mrad); 
		int mz = MinFun(nz-1, iz + mrad);
		
		for( ix = lx ; ix <  mx ; ix++ )
		{
			for( iy = ly ; iy < my ; iy++ )
			{
				for( iz = lz ; iz <  mz ; iz++ )
				{
					float xp, yp, zp;
					mol_dens.GetXYZ(xp, yp , zp, ix, iy, iz);
					double dist = ( xp - xa)*(xp - xa);
					dist+= ( yp - ya)*(yp - ya);
					dist+= ( zp - za)*(zp - za);
					if( dist < frad2)
					{
						float* fld_ptr = mol_dens.GetValPtr(ix,iy,iz);
						(*fld_ptr) = 1.0; 
					}
				}
			}
		}
	}

	mrad = (int)(solv_rad/grid_space) + 1;
	int mrad2 = (int)((solv_rad/grid_space)*(solv_rad/grid_space));

	vector< vector<int> > update_vec;

	for( ix = 0; ix < nx; ix++)
	{
		for( iy = 0; iy < ny; iy++)
		{
			for( iz = 0; iz < nz; iz++)
			{
				if( mol_dens.GetValue(ix,iy,iz) > 0.5 ) continue;

				int lx = MaxFun(0, ix - mrad); 
				int mx = MinFun(nx, ix + mrad);
				int ly = MaxFun(0, iy - mrad); 
				int my = MinFun(ny, iy + mrad);
				int lz = MaxFun(0, iz - mrad); 
				int mz = MinFun(nz, iz + mrad);

				for( int ix2 = lx; ix2 < mx; ix2++)
				{
					for( int iy2 = ly; iy2 < my; iy2++)
					{
						for( int iz2 = lz; iz2 < mz; iz2++)
						{
							if( ((ix2 - ix)*(ix2 - ix) + (iy2 - iy)*(iy2 - iy) + (iz2 - iz)*(iz2 - iz) ) > mrad2 ) continue;
							if( mol_dens.GetValue(ix2,iy2,iz2) < 0.5 ) continue;

							std::vector<int> idx_vec;
							idx_vec.push_back(ix2);
							idx_vec.push_back(iy2);
							idx_vec.push_back(iz2);
							update_vec.push_back(idx_vec);
						}
					}
				}
			}
		}
	}

	int np = update_vec.size();
	for(int i = 0; i < np; i++)
	{
		ix = update_vec[i][0];
		iy = update_vec[i][1];
		iz = update_vec[i][2];
		float* fld_ptr = mol_dens.GetValPtr(ix,iy,iz);
		(*fld_ptr) = 0.0; 
	}

	HaMolecule* pmol = this->GetMolByName("EXCLUDED_VOLUME");
	if( pmol != NULL) this->DeleteMol(pmol);

	pmol = this->AddNewMolecule();
	pmol->SetObjName("EXCLUDED_VOLUME");

	HaChain* pch= pmol->AddChain('A');
	HaResidue* pres = pch->AddResidue(1);
	pres->SetName("MOL");
	int ires = 1;
	int iat = 0;

	map< HaAtom*,vector<int> > aptr_idx_map;
	map< vector<int>,HaAtom* > idx_aptr_map;

	for( ix = 0; ix < mol_dens.GetNx(); ix++)
	{
		for( iy = 0; iy < mol_dens.GetNy(); iy++)
		{
			for( iz = 0; iz < mol_dens.GetNz(); iz++)
			{
				if( mol_dens.GetValue(ix,iy,iz) > 0.5 ) 
				{
					iat++;
					if( iat % 1000 == 0)
					{
						ires++;
						pres = pch->AddResidue(ires);
						iat = 1;
					}
					sprintf(buf,"%d",iat);
					std::string at_name = buf;
					boost::trim(at_name);
					if( iat < 10 ) at_name = (std::string)"0" + at_name;
					if( iat < 100 ) at_name = (std::string)"0" + at_name;
					at_name = (std::string) "A" + at_name;
					aptr = pres->AddNewAtom();
					aptr->SetName(at_name);
					aptr->SetElemNo(1);
					
					Vec3D crd_pt = mol_dens.GetGridPtCrd(ix,iy,iz);
					aptr->SetCoordFrom(crd_pt);

					std::vector<int> idx_grid;
					idx_grid.push_back(ix);
					idx_grid.push_back(iy);
					idx_grid.push_back(iz);
					aptr_idx_map[ aptr ] = idx_grid;
					idx_aptr_map[ idx_grid ] = aptr;
				}
			}
		}
	}

	AtomIteratorMolecule aitr_mol(pmol);
	for( aptr = aitr_mol.GetFirstAtom(); aptr; aptr = aitr_mol.GetNextAtom())
	{
		vector<int> idx_grid = aptr_idx_map[ aptr ];
		for( int i = 0; i < 3; i++)
		{
			vector<int> idx_grid_2 = idx_grid;
			idx_grid_2[i] = idx_grid[i] - 1;
			if( idx_aptr_map.count(idx_grid_2) > 0)
			{
				HaAtom* aptr2 = idx_aptr_map[idx_grid_2];
				HaAtom::CreateBond(aptr,aptr2);
			}
			idx_grid_2[i] = idx_grid[i] + 1;
			if( idx_aptr_map.count(idx_grid_2) > 0)
			{
				HaAtom* aptr2 = idx_aptr_map[idx_grid_2];
				HaAtom::CreateBond(aptr,aptr2);
			}
		}
	}
	
	HaMolView* pview = this->GetActiveMolView();
	if(pview)
	{
		pview->InitialTransform();
	    pview->DefaultRepresentation();
	}
	this->RefreshAllViews( RFRefresh | RFColour | RFApply );
	this->AnnounceGeomChange();


	return TRUE;
}

std::string MolSet::GetUniqueMolName(const std::string& suggest_name)
{
	vector<HaMolecule*>::iterator mol_itr;
	std::string trial_name= suggest_name;
	bool name_found = false;

	int idx = 0;
	char chidx[10];

	while (name_found == false)
	{
		name_found = true;
		for( mol_itr=HostMolecules.begin(); mol_itr != HostMolecules.end(); mol_itr++)
		{
			std::string mname= (*mol_itr)->GetObjName();
			if(mname == trial_name)
			{
				name_found = false;
				idx++;
				if(idx < 10)
					sprintf(chidx,"%1d",idx);
				else if( idx < 100)
					sprintf(chidx,"%2d",idx);
				else if( idx < 1000)
					sprintf(chidx,"%3d",idx);
				else if( idx < 10000)
					sprintf(chidx,"%4d",idx);

				trial_name= suggest_name + chidx;
				break;
			}
		}
	}
	return trial_name;
}

HaMolecule* MolSet::AddNewMolecule( int mol_ser_no )
{
	int last_id = 0;
	if (!serno_mol_map.empty()) last_id = serno_mol_map.rbegin()->first;
	if (mol_ser_no < 0)
	{
		if (HostMolecules.empty())
			mol_ser_no = 1;
		else
			mol_ser_no = last_id + 1;
	}

	if ( serno_mol_map.count(mol_ser_no) > 0)
	{
		PrintLog(" Warning in __FUNCTION__  \n");
		PrintLog(" Molecule Number %d is not unique \n", mol_ser_no);
	}

	HaMolecule* pmol = new HaMolecule(this);
	HostMolecules.push_back(pmol);

	pmol->serno = mol_ser_no;

	std::pair<int, HaMolecule*> im_pair(mol_ser_no, pmol);
	serno_mol_map.insert(im_pair);

	std::pair<std::string, HaMolecule*> nm_pair(pmol->GetName(), pmol);
	name_mol_map.insert(nm_pair);
	
	return(pmol);
} 

HaMolecule* MolSet::GetFirstMolecule()
{
	if(HostMolecules.empty())
		return NULL;
	return HostMolecules[0];
}

HaMolecule* MolSet::GetMolByIdx(int imol)
{
	if(HostMolecules.empty() || imol >= HostMolecules.size())
		return NULL;
	return HostMolecules[imol];
}

HaMolecule* MolSet::GetMolByIdx0(int imol)
{
	if (HostMolecules.empty() || imol >= HostMolecules.size())
		return NULL;
	return HostMolecules[imol];
}

HaMolecule* MolSet::GetMolByIdx1(int imol)
{
	if (HostMolecules.empty() || imol > HostMolecules.size())
		return NULL;
	return HostMolecules[imol-1];
}

double MolSet::OverlapMol(AtomGroup& firstatset, AtomGroup& secatset)
{
//! Superimpose the two structures such that chosen set of atoms
//! \param firstatset atom group of the first molecule
//! \param secatset   atom group of the second molecule corresponding to firstatset of the first molecule
//! \return eps(Ang) - atom RMS deviation in superimposed structures
//! \return -1.0 on error

// check that two sets have the same number of atoms and are not empty

//	PrintLog(" MolSet::OverlapMol() pt 1 \n");
	int nat;
    nat = firstatset.size();
	if (nat == 0)
	{
		PrintMessage("ERROR: First Atom Set is empty \n");
		return -1.0;
	}

	if (secatset.size() != nat)
	{
		PrintMessage("ERROR: number of atoms in two sets must be the same");
		return -1.0;
	}

// check that these sets belong to different molecules
	HaAtom *aptr=0;
	AtomIteratorAtomGroup aitr_fst_atset(&firstatset);
	AtomIteratorAtomGroup aitr_sec_atset(&secatset);

    aptr= aitr_fst_atset.GetFirstAtom();
// which molecule firstatset belongs
	HaMolecule *pMol1=0;
	HaMolecule *pMol2=0;
//check the first set
	pMol1 = aptr->GetHostMol();
    for(aptr= aitr_fst_atset.GetFirstAtom();aptr;aptr = aitr_fst_atset.GetNextAtom())
    {
		if(aptr->GetHostMol() != pMol1)
		{
			PrintMessage("Atoms in set 1 do not belong to the same molecule");
			return -1.0;
		}
	}

// check that second set starts in the second molecule

	aptr = aitr_sec_atset.GetFirstAtom();
	pMol2 = aptr->GetHostMol();

	if (pMol1 == pMol2 )
	{
		PrintMessage(" Sets must be from different molecules");
		return -1.0;
	}

//check the second set
   for(aptr=aitr_sec_atset.GetFirstAtom();aptr;aptr =aitr_sec_atset.GetNextAtom())
    {
		if(aptr->GetHostMol() != pMol2)
		{
			PrintMessage("Atoms in set 2 do not belong to the same molecule");
			return -1.0;
		}
	}

// create variables needed for fortran subroutine and call it	
	HaMat_double rot_mat(3,3,0.0);
    rot_mat(1,1) = 1.0;  rot_mat(2,2) = 1.0; rot_mat(3,3) = 1.0; 
	HaVec_double trans_vec (3,0.0);//fortran style 1D array
    double eps;

	int ires = PointContainer::GetSuperimposeMat( firstatset, secatset, rot_mat,  trans_vec, eps);

    double rms = eps;

	PrintLog("RMS of the fit (Ang)%f \n",rms);

// do coordinate transformation for the second molecule;
	 aptr = aitr_sec_atset.GetFirstAtom();
	 pMol2 = aptr->GetHostMol();

	 double x,y,z, xnew, ynew,znew;

	 AtomIteratorMolecule aitr_secmol(pMol2);

	 for (aptr= aitr_secmol.GetFirstAtom();aptr;aptr= aitr_secmol.GetNextAtom())
	 {
		x= aptr->GetX();
       	y = aptr->GetY();
		z = aptr->GetZ();

		xnew = x*rot_mat(1,1)+y*rot_mat(1,2)+z*rot_mat(1,3)+trans_vec(1);
		ynew = x*rot_mat(2,1)+y*rot_mat(2,2)+z*rot_mat(2,3)+trans_vec(2);
		znew = x*rot_mat(3,1)+y*rot_mat(3,2)+z*rot_mat(3,3)+trans_vec(3);

		aptr->SetX(xnew);
		aptr->SetY(ynew);
		aptr->SetZ(znew);
	 }
	 return rms;
}

typedef struct{
	int nc;
	double eps;
	PtrPtrMap* p_best_match;
	VecPtr* p_atvec1;
	VecPtr* p_atvec2;
} MatchInfo1;


class AtomDestroyer: public AttrDestroyer
//! Don't do anything atoms are allocated externally
{ 
public:
   virtual void destroy(void *p)
   { 
//	   delete p;      
   }
};

class AtomComparator: public AttrComparator
{
public:
	AtomComparator(){}

	virtual bool compatible(void *pa, void *pb)
    {
		HaAtom* at1 = (HaAtom*) pa;
		HaAtom* at2 = (HaAtom*) pb;
		int elno1 = at1->GetElemNo();
		int elno2 = at2->GetElemNo();
		return (elno1 == elno2);
	}
};



bool visit_f1(int n, node_id ni1[], node_id ni2[], void *usr_data)
{
	MatchInfo1* pinfo = (MatchInfo1*) usr_data;

	pinfo->nc++;

	AtomGroup ats1, ats2;

	int i;
	HaAtom* aptr;
	for( i = 0 ; i < n; i++)
	{
		int idx1 = ni1[i];
		aptr = (HaAtom*)(*(pinfo->p_atvec1))[idx1];
		ats1.InsertAtom(aptr);
//		PrintLog(" At %d - %9.3f %9.3f %9.3f \n", 
//        idx1, aptr->GetX_Ang(), aptr->GetY_Ang(), aptr->GetZ_Ang());
		int idx2 = ni2[i];
		aptr = (HaAtom*)(*(pinfo->p_atvec2))[idx2];
		ats2.InsertAtom(aptr);
//		PrintLog(" At %d - %9.3f %9.3f %9.3f \n", 
//        idx2, aptr->GetX_Ang(), aptr->GetY_Ang(), aptr->GetZ_Ang() );
//		PrintLog("\n");
	}	

// create variables needed for fortran subroutine and call it	
	HaMat_double rot_mat(3,3,0.0);
    rot_mat(1,1) = 1.0;  rot_mat(2,2) = 1.0; rot_mat(3,3) = 1.0; 
	HaVec_double trans_vec (3,0.0);//fortran style 1D array
    double eps;

	int ires = PointContainer::GetSuperimposeMat( ats1, ats2, rot_mat,  trans_vec, eps);

	if( ires )
	{
//		PrintLog("\nMatch number %d \n",pinfo->nc);
//		PrintLog("RMS of the fit (Ang)%f \n",eps);
		if( eps < pinfo->eps )
		{
            pinfo->eps = eps;
	        for(i = 0 ; i < n; i++)
			{
				int idx1 = ni1[i];
				HaAtom* aptr1 = (HaAtom*)(*(pinfo->p_atvec1))[idx1];
				
				int idx2 = ni2[i];
				HaAtom* aptr2 = (HaAtom*)(*(pinfo->p_atvec2))[idx2];

				(*(pinfo->p_best_match))[aptr1] = aptr2;  
			}
			PrintLog("\n");
		}
	}
	return false;
}


double
MolSet::AlignOverlapMol(AtomGroup& atset1, HaMolecule* pMol2, PtrPtrMap* fit, HaVec_double* p_trans, HaMat_double* p_rot)
//! Input Parameters:
//! atset1 - list of atoms of the molecule1 that will be superimposed by the corresponding 
//!          atoms of the molecule 2 
//! Corresponding atoms of two molecules are found using graph matching algorithm 
//! finding the best match between two molecules
//! molecules 1 and 2 should have the same number of atoms.
//!
//! Return:   eps - best RMS to fit atset1 by corresponding atoms of molecule 2    
//!                                
//! Optional parameters: fit - map of atom pointers of two molecules in the best fit
//!                      trans,rot - translation vector and rotational matrix 
//!                                   to superimpose molecule 2 to molecule 1
//!                            
{
	int nat;
    nat = atset1.size();
	if (nat == 0)
	{
		PrintMessage("ERROR: First Atom Set is empty \n");
		return FALSE;
	}

	char buf1[256]; 
	char buf2[256];
// check that these sets belong to different molecules
	HaAtom* aptr;
    HaMolecule* pMol1;
// which molecule atset1 belongs
//check the first set
	
	AtomIteratorAtomGroup aitr1(&atset1);
	aptr = aitr1.GetFirstAtom();
	pMol1 = aptr->GetHostMol();

    for(aptr=aitr1.GetFirstAtom();aptr;aptr = aitr1.GetNextAtom())
    {
		if(aptr->GetHostMol() != pMol1)
		{
			PrintMessage("Atoms in set 1 do not belong to the same molecule");
			return FALSE;
		}
	}

	if (pMol1 == pMol2 )
	{
		PrintMessage(" Sets must be from different molecules");
		return FALSE;
	}

	int nm1 = pMol1->GetNAtoms();

	if( nm1 != pMol2->GetNAtoms())
	{
		PrintLog(" Molecules 1 and 2 has different number of atoms \n");
		PrintLog(" They can not be aligned \n");
		return FALSE;
	}

    ARGEdit ed1,ed2;  // Graph bulders for molecules 1 and 2

	AtomIteratorMolecule aitr_m1(pMol1),aitr_m2(pMol2);

	HaBond* bptr;

    PtrIntMap atmap1;
	PtrIntMap atmap2;

	VecPtr atvec1(nm1);
	VecPtr atvec2(nm1);

// Set Graph for the second molecule

	int i = 0;
	for(aptr = aitr_m1.GetFirstAtom(); aptr; aptr = aitr_m1.GetNextAtom())
	{
		ed1.InsertNode(aptr);
		atmap1[aptr] = i;
        atvec1[i] = aptr;
		i++;
	}

	for(aptr = aitr_m1.GetFirstAtom(); aptr; aptr = aitr_m1.GetNextAtom())
	{
		HaAtom::BondIterator bitr = aptr->Bonds_begin();
		for( ; bitr != aptr->Bonds_end(); ++bitr )
		{
			HaBond* bptr = (*bitr);
			HaAtom* at1 = bptr->srcatom;
			HaAtom* at2 = bptr->dstatom;
			if( at2 < at1 ) continue;
			int idx1,idx2;

			if( atmap1.find(at1) != atmap1.end() && atmap1.find(at2) != atmap1.end())
			{
				idx1 = atmap1[at1];
				idx2 = atmap1[at2];
				ed1.InsertEdge(idx1, idx2, NULL);
				ed1.InsertEdge(idx2, idx1, NULL);
			}
		}
	}

// Set Graph for the second molecule

	i = 0;
	for(aptr = aitr_m2.GetFirstAtom(); aptr; aptr = aitr_m2.GetNextAtom())
	{
		ed2.InsertNode(aptr);
		atmap2[aptr] = i;
		atvec2[i] = aptr;
		i++;
	}
	
	for(aptr = aitr_m2.GetFirstAtom(); aptr; aptr = aitr_m2.GetNextAtom())
	{
		HaAtom::BondIterator bitr = aptr->Bonds_begin();
		for( ; bitr != aptr->Bonds_end(); ++bitr )
		{
			HaBond* bptr = (*bitr);
			HaAtom* at1 = bptr->srcatom;
			HaAtom* at2 = bptr->dstatom;
			if( at2 < at1 ) continue;
			int idx1,idx2;

			if( atmap2.find(at1) != atmap2.end() && atmap2.find(at2) != atmap2.end())
			{
				idx1 = atmap2[at1];
				idx2 = atmap2[at2];
				ed2.InsertEdge(idx1, idx2, NULL);
				ed2.InsertEdge(idx2, idx1, NULL);
			}
		}
	}

// Now Graph can be constructed...
//  Graph g1(&ed1);
//	Graph g2(&ed2);

	ARGraph<HaAtom,void> g1(&ed1);
	ARGraph<HaAtom,void> g2(&ed2);

    g1.SetNodeDestroyer(new AtomDestroyer());
    g2.SetNodeDestroyer(new AtomDestroyer());

    g1.SetNodeComparator(new AtomComparator());

// Set inital state for VF2 method:

    VF2SubState s1(&g1, &g2);

	int nc = 0;

	PtrPtrMap best_match;
	MatchInfo1 info1;

	info1.nc = 0;
	info1.eps = 10000.0;
	info1.p_best_match = &best_match;
	info1.p_atvec1 = &atvec1;
	info1.p_atvec2 = &atvec2;

	if(!match(&s1,visit_f1, &info1))
	{
		PrintLog("No match between molecules has been found \n");
		return FALSE;
	}

	PrintLog("Found %d matching between molecules \n", info1.nc);
	PrintLog("RMS of the best fit of whole molecule %f \n",info1.eps);

    if(fit != NULL) fit->clear();

	PrintLog(" Matching atoms of two molecules: \n");
	for( aptr = aitr_m1.GetFirstAtom(); aptr; aptr = aitr_m1.GetNextAtom())
	{
		HaAtom* aptr2 = (HaAtom*)(*(info1.p_best_match))[aptr];
		int idx1 = atmap1[aptr]+1;
		int idx2 = atmap2[aptr2]+1;
		aptr->FillRef(buf1, HaAtom::ATOMREF_NO_MOL);
		aptr2->FillRef(buf2,HaAtom::ATOMREF_NO_MOL);

		if( fit != NULL) fit->SetVal((void*)aptr,(void*)aptr2);

		PrintLog(" %4d %4d        %15s %15s  \n", idx1,idx2,buf1,buf2); 
	}
	PrintLog("\n\n");

	AtomGroup atset2; // AtomGroup of the second molecule corresponding to atset1

//	PrintLog("Atoms of the two sets to superimpose: \n");

	for(aptr = aitr1.GetFirstAtom(); aptr; aptr = aitr1.GetNextAtom() )
	{
		HaAtom* aptr2 = (HaAtom*)(*(info1.p_best_match))[aptr];

		atset2.InsertAtom(aptr2);
	}

	HaMat_double rot_mat(3,3,0.0);
    rot_mat(1,1) = 1.0; rot_mat(2,2) = 1.0; rot_mat(3,3) = 1.0; 
	HaVec_double trans_vec (3,0.0);//fortran style 1D array
	double eps_set;

	int ires = PointContainer::GetSuperimposeMat( atset1, atset2, rot_mat,  trans_vec, eps_set);

	if(!ires)
	{
		ErrorInMod(" MolSet::AlignOverlapMol()","Failed to superimpose atom sets");
		return -1.0;
	}
	if(p_rot != NULL) (*p_rot) = rot_mat;
	if(p_trans != NULL) (*p_trans) = trans_vec;

	PrintLog("RMS of the best fit of selected atom set %f \n", eps_set);

// do coordinate transformation for the second molecule;

	 for (aptr= aitr_m2.GetFirstAtom();aptr;aptr= aitr_m2.GetNextAtom())
	 {
		double x= aptr->GetX();
       	double y = aptr->GetY();
		double z = aptr->GetZ();
		
		double xnew = x*rot_mat(1,1)+y*rot_mat(1,2)+z*rot_mat(1,3)+trans_vec(1);
		double ynew = x*rot_mat(2,1)+y*rot_mat(2,2)+z*rot_mat(2,3)+trans_vec(2);
		double znew = x*rot_mat(3,1)+y*rot_mat(3,2)+z*rot_mat(3,3)+trans_vec(3);

		aptr->SetX(xnew);
		aptr->SetY(ynew);
		aptr->SetZ(znew);
	 }

// Compute RMS of the fit for the whole molecule

	 double rms = 0.0;
	 for(aptr = aitr_m1.GetFirstAtom(); aptr; aptr = aitr_m1.GetNextAtom())
	 {
		 HaAtom* aptr2 = (HaAtom*)(*(info1.p_best_match))[aptr];
		 int idx1 = atmap1[aptr];
		 int idx2 = atmap2[aptr2];
		 rms += (aptr->GetX() - aptr2->GetX())*(aptr->GetX() - aptr2->GetX());
		 rms += (aptr->GetY() - aptr2->GetY())*(aptr->GetY() - aptr2->GetY());
		 rms += (aptr->GetZ() - aptr2->GetZ())*(aptr->GetZ() - aptr2->GetZ());
	 }

	 rms = sqrt(rms /nm1);

	 PrintLog(" RMS of atoms of all atoms of two molecules for a transformation \n");
	 PrintLog(" that superimpose two subsets of the molecules: \n");
	 PrintLog(" RMS_3 = %12.6f \n", rms ); 

	 return eps_set;		
}


AtomIteratorMolSet::AtomIteratorMolSet(MolSet* pmset_new)
{
	pritr = new ResidueIteratorMolSet(pmset_new);	
	GetFirstAtom();
	first_called = false;
}

AtomIteratorMolSet::AtomIteratorMolSet(const AtomIteratorMolSet& ref)
{
	pritr = new ResidueIteratorMolSet(*ref.pritr);
	aitr_res      = ref.aitr_res;
	aitr_res_end  = ref.aitr_res_end;
	first_called  = ref.first_called;
}


AtomIteratorMolSet::~AtomIteratorMolSet()
{
	delete pritr;
}

PointIterator* AtomIteratorMolSet::clone() const
{
	AtomIteratorMolSet* p_aitr = new AtomIteratorMolSet(*this);
	return p_aitr;
}

HaAtom* AtomIteratorMolSet::GetFirstAtom()
{
   HaResidue* pres;
   for( pres = pritr->GetFirstRes(); pres; pres = pritr->GetNextRes() )
   {
	   if(pres->empty()) continue;
	   aitr_res     = pres->begin();
	   aitr_res_end = pres->end(); 
	   return (*aitr_res);
   }
   return NULL;
}

inline HaAtom* AtomIteratorMolSet::GetNextAtom()
{
	aitr_res++;
	if( aitr_res != aitr_res_end ) return (*aitr_res);

	HaResidue* pres = pritr->GetNextRes();
	for(; pres; pres = pritr->GetNextRes())
	{
		aitr_res     = pres->begin();
		aitr_res_end = pres->end();
		if( aitr_res != aitr_res_end ) return (*aitr_res);
	}
	return NULL;
}

AtomIteratorMolSet AtomIteratorMolSet::__iter__() const
{
	return (*this);
}

HaAtom* AtomIteratorMolSet::next()
{
	HaAtom* aptr;
	if(first_called) 
	{
		aptr = GetNextAtom();
	}
	else
	{
		aptr = GetFirstAtom();
		first_called = true;
	}
	if( aptr == NULL)
	{
		throw std::out_of_range("Stop Atom Iterations");	
	}
	return aptr;
}

AtomIteratorMolSet_const::AtomIteratorMolSet_const(const MolSet* pmset_new )
{
	pritr = new ResidueIteratorMolSet_const(pmset_new);   
	GetFirstAtom();
}

AtomIteratorMolSet_const::AtomIteratorMolSet_const(const AtomIteratorMolSet_const& ref)
{
	pritr = new ResidueIteratorMolSet_const(*ref.pritr);
	aitr_res      = ref.aitr_res;
	aitr_res_end  = ref.aitr_res_end;
}

AtomIteratorMolSet_const::~AtomIteratorMolSet_const()
{
	delete pritr;
}

PointIterator_const* AtomIteratorMolSet_const::clone() const
{
	AtomIteratorMolSet_const* p_aitr = new AtomIteratorMolSet_const(*this);
	return p_aitr;
}

const HaAtom* AtomIteratorMolSet_const::GetFirstAtom()
{
   const HaResidue* pres;
   for( pres = pritr->GetFirstRes(); pres; pres = pritr->GetNextRes() )
   {
	   aitr_res     = pres->begin();
	   aitr_res_end = pres->end();
	   if( aitr_res != aitr_res_end ) return (*aitr_res);
   }
   return NULL;
}

const HaAtom* AtomIteratorMolSet_const::GetNextAtom()
{
	aitr_res++;
	if( aitr_res != aitr_res_end ) return (*aitr_res);

	const HaResidue* pres = pritr->GetNextRes();

	for(; pres; pres = pritr->GetNextRes())
	{
	   aitr_res     = pres->begin();
	   aitr_res_end = pres->end();
	   if( aitr_res != aitr_res_end ) return (*aitr_res);
	}
	return NULL;
}

ResidueIteratorMolSet::ResidueIteratorMolSet(MolSet* pmset)
{
	mol_itr_begin = pmset->HostMolecules.begin();
	mol_itr_end   = pmset->HostMolecules.end();

	mol_itr = mol_itr_begin;
	
	if( mol_itr != mol_itr_end )
	{
        ch_itr = (*mol_itr)->Chains.begin(); 
		if( ch_itr != (*mol_itr)->Chains.end() )
		{
           res_itr = (*ch_itr).res_arr.begin();
		}
	}
	first_called = 0;
}

ResidueIteratorMolSet::ResidueIteratorMolSet(const ResidueIteratorMolSet& ref)
{
	res_itr     = ref.res_itr;
	ch_itr      = ref.ch_itr;
	mol_itr     = ref.mol_itr;

	mol_itr_begin = ref.mol_itr_begin;
	mol_itr_end   = ref.mol_itr_end;

	first_called = 0;
}

ResidueIteratorMolSet::~ResidueIteratorMolSet()
{

}

HaResidue* ResidueIteratorMolSet::GetFirstRes()
{	
   mol_itr = mol_itr_begin;
   if( mol_itr == mol_itr_end )
   {
       return NULL;
   }

   while( mol_itr != mol_itr_end )
   {
       if( (*mol_itr)->Chains.empty() ) 
	   { 
		   mol_itr++ ;
	   }
	   else
	   {
		   ch_itr = (*mol_itr)->Chains.begin();
		   break;
	   }
   }

   if( mol_itr == mol_itr_end )
   {
	   return NULL;
   }

   while( ch_itr != (*mol_itr)->Chains.end() )
   {
		if( (*ch_itr).res_arr.empty() )
		{
			ch_itr++;
		}
		else
		{
            res_itr = (*ch_itr).res_arr.begin();
			return (*res_itr);
		}
   }

   return NULL;
}

HaResidue* ResidueIteratorMolSet::GetNextRes()
{
  res_itr++;

  if( res_itr != (*ch_itr).res_arr.end())
  {
	   return (*res_itr);
  }

  ch_itr++;
  while(1)
  {
	  while(ch_itr != (*mol_itr)->Chains.end())
	  {
		  if( (*ch_itr).res_arr.empty() ) 
		  { 
			  ch_itr++;
		  }
		  else
		  {
			  res_itr = (*ch_itr).res_arr.begin();
			  return (*res_itr);
		  }
	  }
	  mol_itr++;
	  if(mol_itr == mol_itr_end) break;
	  ch_itr = (*mol_itr)->Chains.begin();
  }

  return NULL;
}

HaResidue* ResidueIteratorMolSet::GetCurrRes()
{
	if(mol_itr == mol_itr_end) return NULL;
	return( (HaResidue*)&(*res_itr));
}

ResidueIteratorMolSet ResidueIteratorMolSet::__iter__() const
{
	return (*this);
}

HaResidue* ResidueIteratorMolSet::next()
{
	HaResidue* rptr;
	if (first_called)
	{
		rptr = this->GetNextRes();
	}
	else
	{
		rptr = GetFirstRes();
		first_called = 1;
	}
	if (rptr == NULL)
	{
		throw std::out_of_range("Stop Residue Iterations");
	}
	return rptr;
}

HaResidue* ResidueIteratorMolSet::__next__()
{
	return next();
}


ResidueIteratorMolSet_const::ResidueIteratorMolSet_const(const MolSet* pmset)
{
	mol_itr_begin = pmset->HostMolecules.begin();
	mol_itr_end   = pmset->HostMolecules.end();
	
	mol_itr = mol_itr_begin;

	if( mol_itr != mol_itr_end )
	{
        ch_itr = (*mol_itr)->Chains.begin(); 
		if( ch_itr != (*mol_itr)->Chains.end() )
		{
           res_itr = (*ch_itr).res_arr.begin();
		}
	}
}

ResidueIteratorMolSet_const::~ResidueIteratorMolSet_const()
{

}

const HaResidue* ResidueIteratorMolSet_const::GetFirstRes()
{
   mol_itr = mol_itr_begin;
   
   if( mol_itr == mol_itr_end )
   {
       return NULL;
   }

   while( mol_itr != mol_itr_end )
   {
       if( (*mol_itr)->Chains.empty() ) 
	   { 
		   mol_itr++ ;
	   }
	   else
	   {
		   ch_itr = (*mol_itr)->Chains.begin();
		   break;
	   }
   }

   if( mol_itr == mol_itr_end )
   {
	   return NULL;
   }

   while( ch_itr != (*mol_itr)->Chains.end() )
   {
		if( (*ch_itr).res_arr.empty() )
		{
			ch_itr++;
		}
		else
		{
            res_itr = (*ch_itr).res_arr.begin();
			return (*res_itr);
		}
   }

   return NULL;

}

const HaResidue* ResidueIteratorMolSet_const::GetNextRes()
{
  res_itr++;
  if( res_itr != (*ch_itr).res_arr.end())
  {
	   return (*res_itr);
  }

  ch_itr++;
  while(1)
  {
	  while(ch_itr != (*mol_itr)->Chains.end())
	  {
		  if( (*ch_itr).res_arr.empty() ) 
		  { 
			  ch_itr++;
		  }
		  else
		  {
			  res_itr = (*ch_itr).res_arr.begin();
			  return (*res_itr);
		  }
	  }
	  mol_itr++;
	  if(mol_itr == mol_itr_end )
		  break;
	  ch_itr = (*mol_itr)->Chains.begin();
  }

  return NULL;
}

ChainIteratorMolSet::ChainIteratorMolSet(MolSet* new_pmset)
{
	pmset = new_pmset;
	if(pmset == NULL) { return; }

	mol_itr = pmset->HostMolecules.begin();
	if( mol_itr != pmset->HostMolecules.end() )
	{
        ch_itr = (*mol_itr)->Chains.begin(); 
	}
}

ChainIteratorMolSet::~ChainIteratorMolSet()
{

}

HaChain* ChainIteratorMolSet::GetFirstChain()
{
   if(pmset == NULL) { return NULL; }
	
   mol_itr = pmset->HostMolecules.begin();
   if( mol_itr == pmset->HostMolecules.end() )
   {
       return NULL;
   }

   while( mol_itr != pmset->HostMolecules.end())
   {
       if( (*mol_itr)->Chains.empty() ) 
	   { 
		   mol_itr++ ;
	   }
	   else
	   {
		   ch_itr = (*mol_itr)->Chains.begin();
		   return &(*ch_itr);
	   }
   }
   return NULL;
}

HaChain* ChainIteratorMolSet::GetNextChain()
{
  if(pmset == NULL) { return NULL; }
  if(mol_itr == pmset->HostMolecules.end() ) return NULL;

  ch_itr++;

  if( ch_itr != (*mol_itr)->Chains.end()) return &(*ch_itr);

  mol_itr++;

  while( mol_itr != pmset->HostMolecules.end())
  {
      if( (*mol_itr)->Chains.empty() ) 
	  { 
		  mol_itr++ ;
	  }
	  else
	  {
		  ch_itr = (*mol_itr)->Chains.begin();
		  return &(*ch_itr);
	  }
  }  
  return NULL;
}


BondIteratorMolSet::BondIteratorMolSet(MolSet* new_pmset)
{
	pmset = new_pmset;
	if(pmset == NULL) { return; }
	
	bitrm = pmset->Bonds.begin();
}

BondIteratorMolSet::~BondIteratorMolSet()
{

}

HaBond* BondIteratorMolSet::GetFirstBond()
{
   if(pmset == NULL) { return NULL; }
	
   bitrm = pmset->Bonds.begin();

   if( bitrm == pmset->Bonds.end() ) return NULL;
   return (*bitrm);
}

HaBond* BondIteratorMolSet::GetNextBond()
{
  if(pmset == NULL) { return NULL; }
  bitrm++;
  if( bitrm != pmset->Bonds.end())
  {
	   return (*bitrm);
  }

  return NULL;
}

HBondIteratorMolSet::HBondIteratorMolSet(MolSet* new_pmset)
{
	pmset = new_pmset;
	bitrm = pmset->HBonds.begin(); 
}

HBondIteratorMolSet::~HBondIteratorMolSet()
{

}

HaHBond* HBondIteratorMolSet::GetFirstBond()
{
   if(pmset == NULL) { return NULL; }
	
   bitrm = pmset->HBonds.begin();	
   if( bitrm == pmset->HBonds.end() ) return NULL;
   
   return (HaHBond*) &(*bitrm);
}

HaHBond* HBondIteratorMolSet::GetNextBond()
{
   if(pmset == NULL) { return NULL; }
	
   bitrm++;
   if( bitrm == pmset->HBonds.end() ) return NULL;
   
   return (HaHBond*) &(*bitrm);
}

void 
MolSet::ClearPickedAtoms()
{
	picked_atoms.clear();
}

void 
MolSet::DisplaySelectCount()
{
    char buffer[40];
	
	int NumSelected = 0;
	
	HaAtom* aptr;
	
	AtomIteratorMolSet aitr(this);
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if(aptr->Selected())
			NumSelected++;
	}
	
	if( NumSelected==0 )
	{   
		PrintMessage("No atoms selected!");
	} 
	else if( NumSelected > 1 )
	{   
		sprintf(buffer,"%d atoms selected!", NumSelected);
		PrintMessage(buffer);
	} 
	else 
		PrintMessage("1 atom selected!");
}

void MolSet::SelectAtomsInBoundaryBox()
{
	if(!per_bc->IsSet()) return;

	vector<AtomGroup> mols;

	this->p_mol_editor->SplitToMolecules(this,mols);

	Vec3D f;
	int i,j,k;
	HaAtom* aptr;

	int nmol = mols.size();

	for(i = 0; i < nmol; i++)
	{
		int select_mol = TRUE;
		int na = mols[i].size();
		for(j = 0; j < na; j++)
		{
			aptr = mols[i][j];

			for(k = 0; k < 3; k++)
			{
				double f = Vec3D::DotProduct(*aptr,per_bc->recip_ucell[k]) - 0.5;
				int n = (int)f;
				double dn = (double) n; 
				if( (dn < f) && ((f - dn) > 0.5) ) dn = dn + 1.0; 
				if( (dn > f) && ((dn - f) > 0.5) ) dn = dn - 1.0;
				if(fabs(dn) > 0.0001 )
				{
					select_mol = FALSE;
				}
			}
		}
		for(j = 0 ; j < na; j++)
		{
			aptr = mols[i][j];
			if( select_mol )
			{
				aptr->Select();
			}
			else
			{
				aptr->UnSelect();
			}
		}
	}
}

void MolSet::SelectAtomsMask( int mask )
{
    HaBond  *bptr;
    HaAtom  *aptr;

	AtomIteratorMolSet aitr(this);
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
    {
        if( aptr->flag & mask )
        {   
			aptr->Select();
        } 
		else 
			aptr->UnSelect();
	}
    DisplaySelectCount();

	BondIteratorMolSet bitr(this);
    if( HaMolView::ZoneBoth )
    {   
		for(bptr = bitr.GetFirstBond();bptr;bptr = bitr.GetNextBond())
		{
			if( bptr->srcatom->Selected() && bptr->dstatom->Selected() )
			{   
				bptr->Select();
			} 
			else 
				bptr->UnSelect();
		}
    } 
	else
	{
        for(bptr = bitr.GetFirstBond();bptr;bptr = bitr.GetNextBond())
		{
           if( bptr->srcatom->Selected() || bptr->dstatom->Selected() )
           {   
			   bptr->Select();
           } 
		   else 
			   bptr->UnSelect();
		}
	}
}

void MolSet::ClearBackbone()
{
	int nb = BackboneBonds.size();
	int i;
	for( i = 0; i < nb; i++)
	{
		delete BackboneBonds[i];
	}
	BackboneBonds.clear();
}


void MolSet::SelectAtomsExprObj( AtomExpr* expr )
{
    HaBond* bptr;
	HaAtom* aptr;

	AtomIteratorMolSet aitr(this);
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if( expr->EvaluateExprFor(aptr) )
		{   
			aptr->Select();
		} 
		else
		{
			aptr->UnSelect();
		}
	}

    DisplaySelectCount();
	BondIteratorMolSet bitr(this);

    if( HaMolView::ZoneBoth )
    {   
        for(bptr = bitr.GetFirstBond();bptr;bptr = bitr.GetNextBond())
		{
           if( bptr->srcatom->Selected() && bptr->dstatom->Selected() )
           {   
			   bptr->Select();
           } 
		   else bptr->UnSelect();
		}
    } 
	else
	{
        for(bptr = bitr.GetFirstBond();bptr;bptr = bitr.GetNextBond())
		{
           if( bptr->srcatom->Selected() || bptr->dstatom->Selected() )
           {   
			   bptr->Select();
           } 
		   else 
			   bptr->UnSelect();
		}
	}
}

void MolSet::SelectAtomsExpr( const char* expr_str)
{
	CmdParser cmd_pr;
	cmd_pr.SetCmdLine(expr_str);
	AtomExpr* p_expr;
	if( (p_expr = cmd_pr.ParseExpression(0,this)) != NULL )
	{   
		if( !cmd_pr.CurToken )
		{   
			SelectAtomsExprObj(p_expr);
		} 
		else 
		{
			PrintLog("Invalid selection string syntax\n");
		}
		delete p_expr;
	}
}

int MolSet::DescribeSecStruct()
{
	int nmol = GetNMol();
	int i;

	PrintLog("\n Secondary Structure Description: \n\n"); 
	for(i = 0; i < nmol; i++)
	{
		HaMolecule* pMol = HostMolecules[i];
		
		PrintLog(" Molecule  %s :\n\n",pMol->GetObjName());
		
		list<SecStructElement>::iterator fitr;

		PrintLog(" Alpha helicies: \n\n");

		for(fitr = pMol->Features.begin(); fitr != pMol->Features.end(); fitr++)
		{
			if( (*fitr).type == FeatHelix )
			{
				PrintLog(" chain %c   res %d - %d  \n", (*fitr).chain, (*fitr).init, (*fitr).term);
			}
		}

		PrintLog("\n Beta sheets: \n\n");

		for(fitr = pMol->Features.begin(); fitr != pMol->Features.end(); fitr++)
		{
			if( (*fitr).type == FeatSheet )
			{
				PrintLog(" chain %c   res %d - %d  \n", (*fitr).chain, (*fitr).init, (*fitr).term);
			}
		}

		PrintLog("\n Turns: \n\n");

		for(fitr = pMol->Features.begin(); fitr != pMol->Features.end(); fitr++)
		{
			if( (*fitr).type == FeatTurn )
			{
				PrintLog(" chain %c   res %d - %d  \n", (*fitr).chain, (*fitr).init, (*fitr).term);
			}
		}
	}

	return TRUE;
}

int MolSet::PrintHBonds()
{
	int nmol = GetNMol();

	set<HaHBond, less<HaHBond> >::iterator hitr;

	int nhb = GetNHBonds();

	PrintLog("\n HBonds Description \n\n");

	PrintLog(" Total Number of H-Bonds = %d \n\n", nhb);
	
	for( hitr = HBonds.begin(); hitr!= HBonds.end(); hitr++) 
	{
		HaAtom* src = (*hitr).src;
		HaAtom* dst = (*hitr).dst;

		std::string src_id = src->GetRef();
		std::string dst_id = dst->GetRef();

		//PrintLog(" %s -> %s \n", src_id.c_str(), dst_id.c_str());
		PrintLog("%d %d", src->GetSerNo(),dst->GetSerNo()); //<< jose addition 
		fstream hbonds;
		hbonds.open("hbonds.dat",ios::out | ios::app);
		hbonds << src->GetSerNo() << " " << dst->GetSerNo() << "\n";
		//>> jose addition
	}
	PrintLog("\n");

	return TRUE;
}

int MolSet::RenumberSelectedRes(int start_num)
{
	HaChain* p_chain = NULL;
	HaResidue* pres;
	ResidueIteratorMolSet ritr(this);

	std::vector<HaResidue*> sel_res;
	std::set<HaResidue*, less<HaResidue*> > sel_res_set;

	for(pres = ritr.GetFirstRes(); pres; pres = ritr.GetNextRes() )
	{
		if(!pres->HasSelectedAtoms()) continue;
		if( p_chain == NULL ) p_chain = pres->GetHostChain();
		if( pres->GetHostChain() != p_chain)
		{
			PrintLog(" Error in MolSet::RenumberSelectedRes() \n");
			PrintLog(" Selected Residues do not belong to the same chain \n"); 
			return FALSE;
		}
		sel_res.push_back(pres);
		sel_res_set.insert(pres);
	}	

	int nr = sel_res.size();
	if( nr == 0) 
	{
		PrintLog(" MolSet::RenumberSelectedRes() \n");
		PrintLog(" No Residues Selected \n");
		return FALSE;
	}

	std::set<int, less<int> > used_res_num;
	ResidueIteratorChain ritr_ch(p_chain);
	for(pres = ritr_ch.GetFirstRes(); pres; pres = ritr_ch.GetNextRes() )
	{
		if( !pres->HasSelectedAtoms() ) used_res_num.insert( pres->GetSerNo() );
	}
	
	HaVec_int res_num_new;
	res_num_new.resize(nr);
	
	int ir;
	for( ir = 0; ir < nr; ir++)
	{
		int ser_num_new = start_num + ir;
		if( used_res_num.count(ser_num_new) > 0  )
		{
			PrintLog(" Error in MolSet::RenumberSelectedRes() \n"); 
			PrintLog(" New residue number %d  will conflict with an exisitng residue in the same chain \n", ser_num_new);
			PrintLog(" Residue Renumbering will not be performed \n");
			return FALSE;
		}
		res_num_new[ir] = ser_num_new;
	}

	std::multimap<int, HaResidue*, less<int> >::iterator ritr_m1,ritr_m2;

	for(ritr_m1 = p_chain->res_map.begin(); ritr_m1 != p_chain->res_map.end(); )
	{
		HaResidue* pres = (*ritr_m1).second;
		if( sel_res_set.count(pres) > 0 )
		{
			ritr_m2 = ritr_m1;
			ritr_m1++;
			p_chain->res_map.erase(ritr_m2);
		}
		else
		{
			ritr_m1++;
		}
	}

	for( ir = 0; ir < nr; ir++)
	{
		HaResidue* pres = sel_res[ir];
		pres->serno = start_num + ir;
		std::pair<int,HaResidue*> ir_pair(pres->GetSerNo(),pres);
		p_chain->res_map.insert(ir_pair);
	}
	return TRUE;
}

int MolSet::AnnounceGeomChange()
{
	RefreshAllViews(RFApply | RFRefresh);
	return TRUE;
}

TiXmlElement* MolSet::AddXml(TiXmlElement* parent_element, const char* name, int option) const
{
	AtomIteratorMolSet_const aitr(this);
	const HaAtom* aptr;

	TiXmlElement* molset_element = new TiXmlElement("ATOMS");

	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		aptr->AddXml(molset_element);
	}
	parent_element->LinkEndChild(molset_element);

	return molset_element;
}

int MolSet::SaveXML(FILE* file_out, int option) const
{
	TiXmlDocument doc;
	TiXmlDeclaration * decl = new TiXmlDeclaration( "1.0", "", "" );
	doc.LinkEndChild( decl );

	TiXmlElement* root_element = new TiXmlElement("HARLEM_DATA");
		
	doc.LinkEndChild(root_element);

	this->AddXml( root_element );

	bool bres = doc.SaveFile( file_out );
	
	if( !bres ) return FALSE;

	return TRUE;
}

FragmentCreatePars::FragmentCreatePars()
{
	add_hydr = 1;
}

CrdSnapshotIterator MolSet::GetCrdSnapshots()
{
	CrdSnapshotIterator sitr(this);
	return sitr;
}

void MolSet::DeleteCrdSnapshots()
{
	if( !crd_snapshots.empty() )
	{
		int ns = crd_snapshots.size();
		int i;
		for( i = 0; i < ns; i++ )
		{ 
			delete crd_snapshots[i];
		}
	}
	crd_snapshots.clear();
}

CrdSnapshot* MolSet::AddCrdSnapshot(const std::string& snap_name_new )
{
	std::string snap_name = "SNAPSHOT";
	if( !snap_name_new.empty() ) snap_name = snap_name_new;
	CrdSnapshot* p_snap = new CrdSnapshot( this, snap_name );
	crd_snapshots.push_back(p_snap);
	return p_snap;
}

CrdSnapshot* MolSet::AddCrdSnapshotForGroup(const std::string& grp_name, const std::string& snap_name_new )
{
	AtomGroup* p_at_grp = this->GetAtomGroupByID( grp_name.c_str() );
	if( p_at_grp == NULL )
	{
		PrintLog(" Error in MolSet::AddCrdSnapshotForGroup() \n");
		PrintLog(" No atom group %s \n", grp_name.c_str() );
		return NULL;
	}
	std::string snap_name = "GRP_SNAPSHOT";
	if( !snap_name_new.empty() ) snap_name = snap_name_new;
	CrdSnapshot* p_snap = new CrdSnapshot( p_at_grp, snap_name );
	crd_snapshots.push_back(p_snap);
	return p_snap;
}

CrdSnapshot* MolSet::GetCrdSnapshotByName(const char* snp_name, bool create )
{
	int i;
	for(i = 0; i < crd_snapshots.size(); i++)
	{
		if( stricmp_loc(crd_snapshots[i]->GetName().c_str(),snp_name) == 0) return crd_snapshots[i];
	}
	if( create ) return AddCrdSnapshot ( snp_name );
	return NULL;
}

int MolSet::DeleteCrdSnapshot( CrdSnapshot* psnap )
{
	vector<CrdSnapshot*>::iterator itr;
	for( itr = crd_snapshots.begin(); itr != crd_snapshots.end(); )
	{
		if( *itr == psnap ) 
		{
			itr = crd_snapshots.erase(itr);
			delete psnap;
			return TRUE;
		}
		itr++;
	}
	return FALSE;
}

int MolSet::SetCrdFromSnapshot( CrdSnapshot* psnap )
{
	int ires;
	try
	{
		if( psnap == NULL ) throw std::runtime_error("psnap == NULL");
		ires = psnap->SetAtomCrd();
	}
	catch( const std::exception& ex )
	{
		PrintLog("MolSet::SetCrdFromSnapshot() \n");
		PrintLog("%s\n",ex.what());
		return FALSE;
	}
	return ires;
}

int MolSet::SetCrdFromSnapshot( const std::string& snap_name )
{
	CrdSnapshot* psnap = GetCrdSnapshotByName( snap_name.c_str() );
	if( psnap == NULL ) 
	{
		PrintLog("Error in MolSet::SetCrdFromSnapshot() \n");
		PrintLog("No CrdSnapshot with name %s \n",snap_name.c_str());
		return FALSE;
	}
	return SetCrdFromSnapshot( psnap );
}


///////////////////////////////////////////////////////////////////////////////
PyAccMolSetProp::PyAccMolSetProp(MolSet* new_pmset)
{
	this->pmset=new_pmset;
}
PyAccMolSetProp::~PyAccMolSetProp()
{
	//aptr->GetName());
			//Result->append(aptr->GetSerNo())
			//HaResidue* pres= aptr->GetHostRes(); 
			//pres->GetName(),4);
			//pres->GetSerNo() );
			//HaChain* chain= aptr->GetHostChain();
			//chain->ident);
			//aptr->charge
}
std::vector<int>* PyAccMolSetProp::GetAtomsSerNoAsVec()
{
	//Calculate serial number instead aptr->GetSerNo() since it too long
	std::vector<int> *Result=new std::vector<int>();
	Result->reserve(this->pmset->GetNAtoms());
	int ser_no=1;
	if(this->pmset != NULL) 
	{
		HaAtom* aptr;
	    AtomIteratorMolSet aitr(this->pmset);

		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if( this->pmset->p_save_opt_default->save_selected && !aptr->Selected())
				continue;
			Result->push_back(ser_no);
			ser_no++;
		}
	}
	return Result;
}
std::vector<int>* PyAccMolSetProp::GetResidueSerNoAsVec()
{
	std::vector<int> *Result=new std::vector<int>();
	Result->reserve(this->pmset->GetNAtoms());
	if(this->pmset != NULL) 
	{
		HaAtom* aptr;
	    AtomIteratorMolSet aitr(this->pmset);

		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if( this->pmset->p_save_opt_default->save_selected && !aptr->Selected())
				continue;
			HaResidue* pres= aptr->GetHostRes(); 
			Result->push_back(pres->GetSerNo());
		}
	}
	return Result;
}
std::vector<double>* PyAccMolSetProp::GetAtomsChargeAsVec()
{
	std::vector<double> *Result=new std::vector<double>();
	Result->reserve(this->pmset->GetNAtoms());
	if(this->pmset != NULL) 
	{
		HaAtom* aptr;
	    AtomIteratorMolSet aitr(this->pmset);

		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if( this->pmset->p_save_opt_default->save_selected && !aptr->Selected())
				continue;
			Result->push_back(aptr->charge);
		}
	}
	return Result;
}
std::vector<double>* PyAccMolSetProp::GetAtomsRadiusAsVec()
{
	std::vector<double> *Result=new std::vector<double>();
	Result->reserve(this->pmset->GetNAtoms());
	if(this->pmset != NULL) 
	{
		HaAtom* aptr;
	    AtomIteratorMolSet aitr(this->pmset);

		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if( this->pmset->p_save_opt_default->save_selected && !aptr->Selected())
				continue;
			Result->push_back(aptr->radius);
		}
	}
	return Result;
}
std::vector<std::string>* PyAccMolSetProp::GetAtomsNameAsVec()
{
	std::vector<std::string> *Result=new std::vector<std::string>();
	Result->reserve(this->pmset->GetNAtoms());
	if(this->pmset != NULL) 
	{
		HaAtom* aptr;
	    AtomIteratorMolSet aitr(this->pmset);

		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if( this->pmset->p_save_opt_default->save_selected && !aptr->Selected())
				continue;
			Result->push_back(std::string(aptr->GetName()));
		}
	}
	return Result;
}
std::vector<std::string>* PyAccMolSetProp::GetResidueNameAsVec()
{
	std::vector<std::string> *Result=new std::vector<std::string>();
	Result->reserve(this->pmset->GetNAtoms());
	if(this->pmset != NULL) 
	{
		HaAtom* aptr;
	    AtomIteratorMolSet aitr(this->pmset);

		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if( this->pmset->p_save_opt_default->save_selected && !aptr->Selected())
				continue;
			HaResidue* pres= aptr->GetHostRes(); 
			Result->push_back(std::string(pres->GetName()));
		}
	}
	return Result;
}
std::vector<double>* PyAccMolSetProp::GetAtomsCoorXAsVec()
{
	std::vector<double> *Result=new std::vector<double>();
	Result->reserve(this->pmset->GetNAtoms());
	if(this->pmset != NULL) 
	{
		HaAtom* aptr;
	    AtomIteratorMolSet aitr(this->pmset);

		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if( this->pmset->p_save_opt_default->save_selected && !aptr->Selected())
				continue;
			Result->push_back(aptr->GetX());
		}
	}
	return Result;
}
std::vector<double>* PyAccMolSetProp::GetAtomsCoorYAsVec()
{
	std::vector<double> *Result=new std::vector<double>();
	Result->reserve(this->pmset->GetNAtoms());
	if(this->pmset != NULL) 
	{
		HaAtom* aptr;
	    AtomIteratorMolSet aitr(this->pmset);

		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if( this->pmset->p_save_opt_default->save_selected && !aptr->Selected())
				continue;
			Result->push_back(aptr->GetY());
		}
	}
	return Result;
}
std::vector<double>* PyAccMolSetProp::GetAtomsCoorZAsVec()
{
	std::vector<double> *Result=new std::vector<double>();
	Result->reserve(this->pmset->GetNAtoms());
	if(this->pmset != NULL) 
	{
		HaAtom* aptr;
	    AtomIteratorMolSet aitr(this->pmset);

		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if( this->pmset->p_save_opt_default->save_selected && !aptr->Selected())
				continue;
			Result->push_back(aptr->GetZ());
		}
	}
	return Result;
}
std::vector<double>* PyAccMolSetProp::GetAtomsIonExcludedRadiusAsVec(double Rion)
{
	std::vector<double> *Result=new std::vector<double>();
	Result->reserve(this->pmset->GetNAtoms());
	if(this->pmset != NULL) 
	{
		HaAtom* aptr;
	    AtomIteratorMolSet aitr(this->pmset);

		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if( this->pmset->p_save_opt_default->save_selected && !aptr->Selected())
				continue;
			Result->push_back(aptr->radius+Rion);
		}
	}
	return Result;
}
void PyAccMolSetProp::WriteAtomParamFileForPNP(const char *filename,\
	std::vector<int>* ResidueSerNo, std::vector<std::string>* ResidueName, std::vector<int>* AtomsSerNo, std::vector<std::string>* AtomsName, \
	std::vector<double>* AtomsCoorX, std::vector<double>* AtomsCoorY, std::vector<double>* AtomsCoorZ,\
	std::vector<double>* AtomsCharge, std::vector<double>* AtomsRadius,  std::vector<double>* AtomsIER1, std::vector<double>* AtomsIER2,\
	std::vector<float>* SR_A_K, std::vector<float>* SR_N_K, std::vector<float>* SR_A_Cl, std::vector<float>* SR_N_Cl)
{
	FILE *fout;
	fout=fopen(filename,"w");
	fprintf(fout,"@ResNum    ResName    AtmNum  AtmName               x               y               z               Q           Rdiel           Rier1           Rier2           SR-A1           SR-N1           SR-A2           SR-N2\n");
	int i;
	int Natoms=AtomsSerNo->size();
	for(i=0;i<Natoms;i++)
	{
		//std::string resname=;
		fprintf(fout,"%9d %8s %9d %8s ",(*ResidueSerNo)[i],(*ResidueName)[i].c_str(),(*AtomsSerNo)[i],(*AtomsName)[i].c_str());
		fprintf(fout,"%15.4f %15.4f %15.4f ",(*AtomsCoorX)[i],(*AtomsCoorY)[i],(*AtomsCoorZ)[i]);
		fprintf(fout,"%15.4f %15.4f %15.4f %15.4f ",(*AtomsCharge)[i],(*AtomsRadius)[i],(*AtomsIER1)[i],(*AtomsIER2)[i]);
		fprintf(fout,"%15.4f %15.4f %15.4f %15.4f\n",(*SR_A_K)[i],(*SR_N_K)[i],(*SR_A_Cl)[i],(*SR_N_Cl)[i]);
	}
    fclose(fout);
}
///////////////////////////////////////////////////////////////////////////////

#if defined(HA_NOGUI)
void MolSet::RefreshAllViews(long lHint)
{

}
void MolSet::SetName(const char* new_name)
{
    name_mset = new_name;
}
#endif
