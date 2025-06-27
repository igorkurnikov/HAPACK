/*! \file haatom.cpp

   Classes to define Atom object in HARLEM.

  \author Igor Kurnikov 
  \date 1997-2002

*/

#define HAATOM_CPP

#include <assert.h>
#include <math.h>

#include <mpi.h>

#include "haio.h"

#include <boost/algorithm/string.hpp>

#include "tinyxml.h"

#include "haatom.h"
#include "haatgroup.h"
#include "haatombasis.h"
#include "g94_globals.h"
#include "hamolecule.h"
#include "hamolset.h"
#include "haresdb.h"

#include "mm_force_field.h"
#include "mm_elements.h"

ElemStruct ElementArr[MAXELEMNO] =  {
    { { ' ', ' ' }, 0.72,  1.44,  1.000,   12, ""             },  /*   0 */
    { { 'H', ' ' }, 0.32,  1.1 ,  1.0079,    4, "HYDROGEN"     },  /*   1 */
    { { 'H', 'e' }, 1.6,   2.2,   4.0026,    5, "HELIUM"       },  /*   2 */
    { { 'L', 'i' }, 0.68,  1.22,  6.9400,   14, "LITHIUM"      },  /*   3 */
    { { 'B', 'e' }, 0.352, 0.628, 9.0122,  12, "BERYLLIUM"    },  /*   4 */
    { { 'B', ' ' }, 0.832, 1.548, 10.810,  13, "BORON"        },  /*   5 */
    { { 'C', ' ' }, 0.72 , 1.548, 12.011,   0, "CARBON"       },  /*   6 */
    { { 'N', ' ' }, 0.68 , 1.4  , 14.0067,  1, "NITROGEN"     },  /*   7 */
    { { 'O', ' ' }, 0.68 , 1.348, 15.9994,  2, "OXYGEN"       },  /*   8 */
    { { 'F', ' ' }, 0.64 , 1.3  , 18.9964,  6, "FLUORINE"     },  /*   9 */
    { { 'N', 'e' }, 1.12 , 2.02 , 20.17,   12, "NEON"         },  /*  10 */
    { { 'N', 'a' }, 0.972, 2.2  , 22.9898,  7, "SODIUM"       },  /*  11 */
    { { 'M', 'g' }, 1.1  , 1.5  , 24.305,  15, "MAGNESIUM"    },  /*  12 */
    { { 'A', 'l' }, 1.352, 1.5  , 26.9815,  9, "ALUMINIUM"    },  /*  13 */
    { { 'S', 'i' }, 1.2  , 2.2  , 26.086,   6, "SILICON"      },  /*  14 */
    { { 'P', ' ' }, 1.036, 1.88 , 30.9738,  8, "PHOSPHORUS"   },  /*  15 */  /* 262? */
    { { 'S', ' ' }, 1.02 , 1.808, 32.0600,  3, "SULFUR"       },  /*  16 */
    { { 'C', 'l' }, 1    , 1.748, 35.453,  13, "CHLORINE"     },  /*  17 */
    { { 'A', 'r' }, 1.568, 2.768, 39.948,  12, "ARGON"        },  /*  18 */
    { { 'K', ' ' }, 1.328, 2.388, 39.090,  12, "POTASSIUM"    },  /*  19 */
    { { 'C', 'a' }, 0.992, 1.948, 40.08,    9, "CALCIUM"      },  /*  20 */
    { { 'S', 'c' }, 1.44 , 1.32 , 44.9559, 12, "SCANDIUM"     },  /*  21 */
    { { 'T', 'i' }, 1.472, 1.948, 47.9000,  9, "TITANIUM"     },  /*  22 */
    { { 'V', ' ' }, 1.328, 1.06 , 50.941,  12, "VANADIUM"     },  /*  23 */
    { { 'C', 'r' }, 1.352, 1.128, 51.966,   9, "CHROMIUM"     },  /*  24 */
    { { 'M', 'n' }, 1.352, 1.188, 54.938,   9, "MANGANESE"    },  /*  25 */
    { { 'F', 'e' }, 1.34 , 1.948, 55.847,   8, "IRON"         },  /*  26 */
    { { 'C', 'o' }, 1.328, 1.128, 58.9332, 12, "COBALT"       },  /*  27 */
    { { 'N', 'i' }, 1.62 , 1.24 , 58.71,   10, "NICKEL"       },  /*  28 */  /* >375! */
    { { 'C', 'u' }, 1.52 , 1.148, 63.546,  10, "COPPER"       },  /*  29 */
    { { 'Z', 'n' }, 1.448, 1.148, 65.38,   10, "ZINC"         },  /*  30 */
    { { 'G', 'a' }, 1.22 , 1.548, 69.72,   12, "GALLIUM"      },  /*  31 */
    { { 'G', 'e' }, 1.168, 3.996, 72.59,   12, "GERMANIUM"    },  /*  32 */  /* 1225? */
    { { 'A', 's' }, 1.208, 0.828, 74.9216, 12, "ARSENIC"      },  /*  33 */
    { { 'S', 'e' }, 1.22 , 0.9  , 78.96,   12, "SELENIUM"     },  /*  34 */
    { { 'B', 'r' }, 1.208, 1.748, 79.904,  10, "BROMINE"      },  /*  35 */
    { { 'K', 'r' }, 1.6  , 1.9  , 83.80,   12, "KRYPTON"      },  /*  36 */
    { { 'R', 'b' }, 1.472, 2.648, 85.467,  12, "RUBIDIUM"     },  /*  37 */
    { { 'S', 'r' }, 1.12 , 2.02 , 87.62,   12, "STRONTIUM"    },  /*  38 */
    { { 'Y', ' ' }, 1.78 , 1.608, 88.9059, 12, "YTTRIUM"      },  /*  39 */
    { { 'Z', 'r' }, 1.56 , 1.42 , 91.22,   12, "ZIRCONIUM"    },  /*  40 */
    { { 'N', 'b' }, 1.48 , 1.328, 92.9064, 12, "NIOBIUM"      },  /*  41 */
    { { 'M', 'o' }, 1.472, 1.748, 95.94,   12, "MOLYBDENUM"   },  /*  42 */
    { { 'T', 'c' }, 1.352, 1.8  , 98.9062, 12, "TECHNETIUM"   },  /*  43 */
    { { 'R', 'u' }, 1.6  , 2.0  , 101.07,  12, "RUTHENIUM"    },  /*  44  Cov & Vdw was 1.4, 1.2 Ang*/
    { { 'R', 'h' }, 1.448, 1.22 , 102.9055,12, "RHODIUM"      },  /*  45 */
    { { 'P', 'd' }, 1.5  , 1.44 , 106.4,   12, "PALLADIUM"    },  /*  46 */
    { { 'A', 'g' }, 1.592, 1.548, 107.868,  9, "SILVER"       },  /*  47 */
    { { 'C', 'd' }, 1.688, 1.748, 112.40,  12, "CADMIUM"      },  /*  48 */
    { { 'I', 'n' }, 1.632, 1.448, 114.82,  12, "INDIUM"       },  /*  49 */
    { { 'S', 'n' }, 1.46 , 1.668, 118.69,  12, "TIN",         },  /*  50 */
    { { 'S', 'b' }, 1.46 , 1.12 , 121.75,  12, "ANTIMONY"     },  /*  51 */
    { { 'T', 'e' }, 1.472, 1.26 , 127.60,  12, "TELLURIUM"    },  /*  52 */
    { { 'I', ' ' }, 1.4  , 1.748, 126.9045,11, "IODINE"       },  /*  53 */
    { { 'X', 'e' }, 1.7  , 2.1  , 131.30,  12, "XENON"        },  /*  54 */
    { { 'C', 's' }, 1.672, 3.008, 132.9054,12, "CAESIUM"      },  /*  55 */
    { { 'B', 'a' }, 1.34 , 2.408, 137.34,   8, "BARIUM"       },  /*  56 */
    { { 'L', 'a' }, 1.872, 1.828, 138.9055,12, "LANTHANUM"    },  /*  57 */
    { { 'C', 'e' }, 1.832, 1.86 , 140.12,  12, "CERIUM"       },  /*  58 */
    { { 'P', 'r' }, 1.82 , 1.62 , 140.9077,12, "PRASEODYMIUM" },  /*  59 */
    { { 'N', 'd' }, 1.808, 1.788, 144.24,  12, "NEODYMIUM"    },  /*  60 */
    { { 'P', 'm' }, 1.8  , 1.76 , 145.0,   12, "PROMETHIUM"   },  /*  61 */
    { { 'S', 'm' }, 1.8  , 1.74 , 150.4,   12, "SAMARIUM"     },  /*  62 */
    { { 'E', 'u' }, 1.992, 1.96 , 151.96,  12, "EUROPIUM"     },  /*  63 */
    { { 'G', 'd' }, 1.792, 1.688, 157.25,  12, "GADOLINIUM"   },  /*  64 */
    { { 'T', 'b' }, 1.76 , 1.66 , 158.9254,12, "TERBIUM"      },  /*  65 */
    { { 'D', 'y' }, 1.752, 1.628, 162.50,  12, "DYSPROSIUM"   },  /*  66 */
    { { 'H', 'o' }, 1.74 , 1.608, 164.9304,12, "HOLMIUM"      },  /*  67 */
    { { 'E', 'r' }, 1.728, 1.588, 167.26,  12, "ERBIUM"       },  /*  68 */
    { { 'T', 'm' }, 1.72 , 1.568, 168.9342,12, "THULIUM"      },  /*  69 */
    { { 'Y', 'b' }, 1.94 , 1.54 , 173.04,  12, "YTTERBIUM"    },  /*  70 */
    { { 'L', 'u' }, 1.72 , 1.528, 174.97,  12, "LUTETIUM"     },  /*  71 */
    { { 'H', 'f' }, 1.568, 1.4  , 178.49,  12, "HAFNIUM"      },  /*  72 */
    { { 'T', 'a' }, 1.432, 1.22 , 280.947, 12, "TANTALUM"     },  /*  73 */
    { { 'W', ' ' }, 1.368, 1.26 , 183.85,  12, "TUNGSTEN"     },  /*  74 */
    { { 'R', 'e' }, 1.352, 1.3  , 186.2,   12, "RHENIUM"      },  /*  75 */
    { { 'O', 's' }, 1.368, 1.58 , 190.2,   12, "OSMIUM"       },  /*  76 */
    { { 'I', 'r' }, 1.32 , 1.22 , 192.22,  12, "IRIDIUM"      },  /*  77 */
    { { 'P', 't' }, 1.5  , 1.548, 195.09,  12, "PLATINUM"     },  /*  78 */
    { { 'A', 'u' }, 1.5  , 1.448, 196.9665, 6, "GOLD"         },  /*  79 */
    { { 'H', 'g' }, 1.7  , 1.98 , 200.59,  12, "MERCURY"      },  /*  80 */
    { { 'T', 'l' }, 1.552, 1.708, 204.37,  12, "THALLIUM"     },  /*  81 */
    { { 'P', 'b' }, 1.54 , 2.16 , 207.200, 12, "LEAD"         },  /*  82 */
    { { 'B', 'i' }, 1.54 , 1.728, 208.9808,12, "BISMUTH"      },  /*  83 */
    { { 'P', 'o' }, 1.68 , 1.208, 209.0,   12, "POLONIUM"     },  /*  84 */
    { { 'A', 't' }, 1.208, 1.12 , 210.0,   12, "ASTATINE"     },  /*  85 */
    { { 'R', 'n' }, 1.9  , 2.3  , 222.0,   12, "RADON"        },  /*  86 */
    { { 'F', 'r' }, 1.8  , 3.24 , 223.0,   12, "FRANCIUM"     },  /*  87 */
    { { 'R', 'a' }, 1.432, 2.568, 226.0254,12, "RADIUM"       },  /*  88 */
    { { 'A', 'c' }, 1.18 , 2.12 , 227.00,  12, "ACTINIUM"     },  /*  89 */
    { { 'T', 'h' }, 1.02 , 1.84 , 232.0381,12, "THORIUM"      },  /*  90 */
    { { 'P', 'a' }, 0.888, 1.6  , 231.0359,12, "PROTACTINIUM" },  /*  91 */
    { { 'U', ' ' }, 0.968, 1.748, 238.029, 12, "URANIUM"      },  /*  92 */
    { { 'N', 'p' }, 0.952, 1.708, 237.0482,12, "NEPTUNIUM"    },  /*  93 */
    { { 'P', 'u' }, 0.928, 1.668, 244.0,   12, "PLUTONIUM"    },  /*  94 */
    { { 'A', 'm' }, 0.92 , 1.66 , 243.0,   12, "AMERICIUM"    },  /*  95 */
    { { 'C', 'm' }, 0.912, 1.648, 247.0,   12, "CURIUM"       },  /*  96 */
    { { 'B', 'k' }, 0.9  , 1.64 , 247.0,   12, "BERKELIUM"    },  /*  97 */
    { { 'C', 'f' }, 0.888, 1.628, 251.0,   12, "CALIFORNIUM"  },  /*  98 */
    { { 'E', 's' }, 0.88 , 1.62 , 254.0,   12, "EINSTEINIUM"  },  /*  99 */
    { { 'F', 'm' }, 0.872, 1.608, 257.0,   12, "FERMIUM"      },  /* 100 */
    { { 'M', 'd' }, 0.86 , 1.6  , 258.0,   12, "MENDELEVIUM"  },  /* 101 */
    { { 'N', 'o' }, 0.848, 1.588, 259.0,   12, "NOBELIUM"     },  /* 102 */
    { { 'L', 'r' }, 0.84 , 1.58 , 260.0,   12, "LAWRENCIUM"   }   /* 103 */ /* Lw? */
        };

std::string HaAtom::GetRef(int ref_type) const
{
	char buf[256];
	FillRef(buf, ref_type);
	return buf;	
}


HaAtom::HaAtom()
{
	ps_ff_par = std::make_shared<AtomFFParam>();
	Clear();
}

HaAtom::HaAtom(const HaAtom& atom_ref)
{
	ps_ff_par = std::make_shared<AtomFFParam>();
	this->SetParamFrom(atom_ref);
}

HaAtom::~HaAtom()
{

}

void HaAtom::Clear()
{
	x=y=z=0;
	
	tempf=0.0;
	col=0;
	elemno = DUMMY_ELEM;
	proxy_flag = false;
	radius = 1.0;
	image_radius = 1.0;
	refno=0;
	flag=0;
	Select();
	flag |= NonBondFlag;
	flag |= NormAtomFlag;
	SetDisplayed(true);

	irad=0;
	solv_access_area = 0.0;

	bonds.clear();

	hybrid = NO_HYBRID;

	phost_res=nullptr;
}


bool HaAtom::SetParamFrom( const HaAtom& atom_ref)
{
	pos[0] =   atom_ref.pos[0];
	pos[1] =   atom_ref.pos[1];
	pos[2] =   atom_ref.pos[2];

	*ps_ff_par = *(atom_ref.ps_ff_par);
	
	x=atom_ref.x;
	y=atom_ref.y;
	z=atom_ref.z;
	
	tempf=atom_ref.tempf;
	col=atom_ref.col;
	label=atom_ref.label;
	elemno=atom_ref.elemno;
	radius = atom_ref.radius;
    image_radius = atom_ref.image_radius;
	refno=atom_ref.refno;
	flag=atom_ref.flag;
	irad=atom_ref.irad;
	hybrid = atom_ref.hybrid;

	return true;	
}

int HaAtom::GetElemNo(void) const
{
  return elemno;
}

void HaAtom::SetElemNo(const int new_elemno)
{
	elemno = new_elemno;

	radius = HaAtom::ElemVDWRadius(elemno);
	this->SetVdWRad( HaAtom::ElemVDWRadius(elemno) );
	image_radius = HaAtom::ElemVDWRadius(elemno);
	this->SetMass( HaAtom::StdElemMass(elemno));
	if(elemno == 1) 
	{
		this->flag |= HydrogenFlag;
	}
	else
	{
		this->flag &= ~HydrogenFlag;
	}
}
void HaAtom::SetDummy()
{
	SetElemNo( DUMMY_ELEM );
}

bool HaAtom::IsDummy() const
{
	return( GetElemNo() == DUMMY_ELEM );
}

bool HaAtom::IsProxy() const
{
	return proxy_flag;
}
  
void HaAtom::SetProxy(bool proxy_flag_new)
{
	proxy_flag = proxy_flag_new;
}

std::string HaAtom::GetReplacedAtName() const
{
	return replaced_atom_name;
}
  
void HaAtom::SetReplacedAtName( const std::string& repl_atom_name )
{
	replaced_atom_name = repl_atom_name;
}

std::string HaAtom::GetDescription() const
{
	return desc;
}
  
void HaAtom::SetDescription( const std::string& desc_new )
{
	desc = desc_new;
}

void HaAtom::SetNameFast(const std::string& atname_par )
{
	refno = RegisterAtName(atname_par);
}

void HaAtom::SetName(const std::string& atname_par )
{
	std::string atname = boost::trim_copy(atname_par);
	boost::to_upper(atname);

	if(atname.size() == 0) { atname = "C"; }

// Handle AMBER PDB atom names (with a digit in front) 
	if( isdigit(atname[0]) && atname.size() > 1)
	{
		std::string tmp =  atname.substr(1); 
		tmp += atname[0];
		atname = tmp;
	}
	HaResidue* pres = this->GetHostRes();

	if( pres ) atname = MMForceField::GetAtNameFromAmber( atname, pres->GetFullName() ); 

	std::string modified_atname;
	for(int i = 0; i < atname.size(); i++)
	{
		if(isspace(atname[i])) continue;
		if(atname[i] == '*') 
			modified_atname += "X";
		else
			modified_atname += atname[i];
	}

	refno = RegisterAtName(modified_atname);

//  if(GetElemNo() == 0) SetElemNo( GetElemNoFromName(ElemDesc[refno_new].c_str(), this->GetHostRes() ) );

}

const char* HaAtom::GetName() const
{
	return( ElemDesc[this->refno].c_str() );
}

int HaAtom::GetNBonds() const
{
	return bonds.size();
}

std::string HaAtom::GetStdSymbol() const
{
	return(GetStdSymbolElem(elemno));
}



HaResidue* HaAtom::GetHostRes()
{
	return phost_res;
}

const HaResidue* HaAtom::GetHostRes() const
{
	return phost_res;
}



HaChain* HaAtom::GetHostChain()
{
	return phost_res->GetHostChain();
}

const HaChain* HaAtom::GetHostChain() const
{
	return phost_res->GetHostChain();
}

HaMolecule* HaAtom::GetHostMol()
{
	HaChain* phost_ch=GetHostChain();
	assert(phost_ch != NULL);
	return phost_ch->GetHostMol();
}

const HaMolecule* HaAtom::GetHostMol() const
{
	const HaChain* phost_ch=GetHostChain();
	assert(phost_ch != NULL);
	return phost_ch->GetHostMol();
}

MolSet* HaAtom::GetHostMolSet()
{
	HaMolecule* phost_mol=GetHostMol();
	assert(phost_mol != NULL);
	return phost_mol->GetHostMolSet();
}

int HaAtom::GetSerNo() const
{
    const MolSet* pmset = GetHostMolSet();
	const HaAtom* aptr; 
	int ser_no = 0;
    AtomIteratorMolSet_const aitr(pmset);
	for(aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom())
	{
		ser_no++;
		if(aptr == this) { return ser_no; }
	}
	return 0;
}

const MolSet* HaAtom::GetHostMolSet() const
{
	const HaMolecule* phost_mol=GetHostMol();
	assert(phost_mol != NULL);
	return phost_mol->GetHostMolSet();
}


ChemGroup* HaAtom::GetHostChemGroup()
{
	MolSet* pmset=GetHostMolSet();

	ChemGroup* gptr;
	ChemGroupIterator gitr(pmset);
	
	for(gptr=gitr.GetFirst();gptr;gptr=gitr.GetNext())
	{
		if(gptr->HasAtom(this)) return gptr;
	}
	return NULL;
}

void HaAtom::SetHostRes(HaResidue* new_phost_res)
{
	phost_res=new_phost_res;
}

bool HaAtom::CreateBond(HaAtom* aptr1, HaAtom* aptr2)
{
	MolSet* pmset = aptr1->GetHostMolSet();
	if( aptr2->GetHostMolSet() != pmset ) return false;

	HaBond* bptr= pmset->AddBond(aptr1,aptr2 );
	if(bptr != NULL) bptr->DrawWire();

	return false;
}

bool HaAtom::DeleteBond(HaAtom* aptr1, HaAtom* aptr2)
{
	MolSet* pmset = aptr1->GetHostMolSet();
	pmset->DeleteBond(aptr1,aptr2);
	return true;
}

int HaAtom::GetBondedAtoms(AtomGroup& bonded_atoms_out)
{
	bonded_atoms_out.clear();
	auto bitr = bonds.begin();
	for(; bitr != bonds.end(); bitr++)
	{
		if( (*bitr)->srcatom != this ) bonded_atoms_out.push_back( (*bitr)->srcatom );
		else bonded_atoms_out.push_back( (*bitr)->dstatom );
	}
	return( bonded_atoms_out.size() );
}

AtomGroup HaAtom::GetBondedAtoms()
{
	AtomGroup group;
	this->GetBondedAtoms(group);
	return group;
}

//int HaAtom::GetBonds(std::vector<HaBond*>& bonds_out )
//{
//	bonds_out = (*p_bonds);
//	return( bonds_out.size() );
//}
//
//int HaAtom::GetBonds(std::vector<const HaBond*>& bonds_out ) const 
//{
//	int nb = p_bonds->size();
//	int i;
//	bonds_out.resize(nb);
//	for( i = 0; i < nb; i++)
//	{
//		bonds_out[i] = p_bonds->at(i);
//	}
//	return( bonds_out.size() );
//}

int HaAtom::GetHBondAcc(AtomGroup& bonded_atoms_out)
{
	bonded_atoms_out.clear();
	
	MolSet* pmset = this->GetHostMolSet();

	HaHBond hbond_start(this, (HaAtom*)1);

	std::set<HaHBond>::iterator ritr;

	ritr= pmset->HBonds.lower_bound(hbond_start);

	int nb = 0;

	for( ; ritr != pmset->HBonds.end(); ritr++)
	{
		if( (*ritr).src != this) return nb;
		HaAtom* dst = (*ritr).dst;
		bonded_atoms_out.InsertAtom( (*ritr).dst );
		nb++;
	}

	return nb;
}

bool HaAtom::IsBonded(HaAtom& atm2)
{
	if( &atm2 == this ) return NULL;
	
	auto bitr = bonds.begin();
	for( ; bitr != bonds.end(); bitr++ )
	{
		if( (*bitr)->srcatom == &atm2 || (*bitr)->dstatom == &atm2 ) return true; 		
	}
	return false;
}

void HaAtom::RemoveBond(HaBond* pb)
{
	for( auto bitr = bonds.begin(); bitr != bonds.end(); )
	{
		if( (*bitr).get() == pb ) 
		{
			bitr = bonds.erase(bitr);
		}
		else
		{
			bitr++;
		}
	}
}

int HaAtom::GetReachableAtoms(AtomGroup& block_atoms, HaAtom* aptr2, AtomGroup& reach_atoms, int& loop)
{
	loop = FALSE;
    reach_atoms.clear();

	std::set<HaAtom*> ratom_set;
	std::stack<HaAtom*> pending_atoms;
	AtomGroup bonded_atoms;
	HaAtom* aptr;
	pending_atoms.push(aptr2);
	ratom_set.insert(aptr2);
	aptr2->GetBondedAtoms(bonded_atoms);
	AtomIteratorAtomGroup at_itr(&bonded_atoms);
	for( aptr = at_itr.GetFirstAtom(); aptr; aptr = at_itr.GetNextAtom())
	{
		if(!block_atoms.HasAtom(aptr)) 
		{
			pending_atoms.push(aptr);
		}
	}

	while( !pending_atoms.empty())
    {
		HaAtom* aptr = pending_atoms.top();
		pending_atoms.pop();
		if( ratom_set.count(aptr) == 0)
		{
			if(block_atoms.HasAtom(aptr)) // loop!
			{
				loop = TRUE;
				return FALSE;
			}
			ratom_set.insert(aptr);
			aptr->GetBondedAtoms(bonded_atoms);
			AtomIteratorAtomGroup aitr(&bonded_atoms);
			for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
			{
				pending_atoms.push(aptr);
			}
		}
	}

	std::set<HaAtom*>::iterator sa_itr, send;

	sa_itr = ratom_set.begin();
	send   = ratom_set.end();
	
	for( ; sa_itr != send; sa_itr++)
	{
		reach_atoms.InsertAtom(*sa_itr);
	}

	return TRUE;
}

static const int NATMAX=110;

static std::string StdAtmLbl[NATMAX]= {  "X",                        // 0 -Dummy
   "H",                               "He",                       // 1-2
   "Li","Be","B", "C", "N", "O", "F", "Ne",                       // 3-10
   "Na","Mg","Al","Si","P", "S", "Cl","Ar",                       // 11-18
   "K", "Ca",                                                     // 19-20
   "Sc","Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn",             // 21-30
             "Ga","Ge","As","Se","Br","Kr",                       // 31-36
   "Rb","Sr",                                                     // 37-38
   "Y" ,"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",             // 39-48
             "In","Sn","Sb","Te","I" ,"Xe",                       // 49-54
   "Cs","Ba","La",                                                // 55-57
   "Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu", // 58-71
        "Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg",             // 72-80
		     "Tl","Pb","Bi","Po","At","Rn",                       // 81-86
   "Fr","Ra","Ac",                                                // 87-89
   "Th","Pa","U" ,"Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr", // 90-103
        "Rf","Ha","Sg","Ns","Hs","Mt"                             // 104-109
};


int HaAtom::AtTypeFromLbl(const std::string & Label)
{
	for(int i=1; i < NATMAX; i++)
	{
		if(Label == StdAtmLbl[i])
			return i;
	}
	return 0;
}


bool HaAtom::Print_info(std::ostream &sout, const int level) const
{
	assert(level >= 0);
	sout << "element #= " << GetElemNo() << " Name = " << GetName() << "\n";

	sout << " coord= " << pos[0] << "  "
		               << pos[1] << "  "
					   << pos[2] << "\n";
    return true;
}

bool HaAtom::FillRef(char* buf, int ref_type) const
{
	int j,i;

	const HaChain*    phost_chain=GetHostChain();
	const HaMolecule* phost_mol= GetHostMol();
	const MolSet*   pmset = GetHostMolSet();

	const char* res_name = phost_res->GetName();
	
	j=0;

	if( ref_type == ATOMREF_FULL || (ref_type == ATOMREF_STD && pmset->GetNMol() > 1) )
	{
		j+= sprintf(buf+j,"$%s$",phost_mol->GetObjName());
	}
	
	if( ref_type != ATOMREF_NO_RES )
	{
		for(i=0; i<4; i++)
		{
			if( res_name[i] == 0) break;
			if( res_name[i] != ' ')j += sprintf(buf+j,"%c", res_name[i]);
		}
		j += sprintf(buf+j,"%d",phost_res->GetSerNo());
		
		if( phost_chain->ident != ' ' )
		{
			j += sprintf(buf+j, ":%c",phost_chain->ident);
		}

		j += sprintf(buf+j,".%s", GetName() );
	}
	else
	{
        j += sprintf(buf+j,"%s",GetName() );
	}
	
	*(buf+j)=0;
	
	return true;	
}

StrVec HaAtom::ElemDesc;
StrIntMap HaAtom::at_name_refno_map;

int HaAtom::RegisterAtName(const std::string& at_name)
{
	StrIntMap::iterator itr = at_name_refno_map.find(at_name);
	if( itr != at_name_refno_map.end() ) return (*itr).second;
	ElemDesc.push_back(at_name);
	int refno = ElemDesc.size() - 1;
	at_name_refno_map[at_name] = refno;
	return refno;
}

void HaAtom::FillStdAtomTypes()
{
	if (ElemDesc.size() > 10) return;
	ElemDesc.reserve(200);
	RegisterAtName("N");    // 0
	RegisterAtName("CA");   // 1
	RegisterAtName("C");    // 2
	RegisterAtName("O");    // 3   0-3   Amino Acid Backbone  
	RegisterAtName("C\\");   // 4
	RegisterAtName("OT");   // 5
	RegisterAtName("S");    // 6
	RegisterAtName("P");    // 7   4-7   Shapely Amino Backbone
	RegisterAtName("O1P");  // 8 
	RegisterAtName("O2P");  // 9 
	RegisterAtName("O5X");  // 10 
	RegisterAtName("C5");   // 11 
	RegisterAtName("C4");   // 12
	RegisterAtName("O4");   // 13
    RegisterAtName("C3");   // 14
	RegisterAtName("O3");   // 15
	RegisterAtName("C2");   // 16
    RegisterAtName("O2");   // 17
	RegisterAtName("C1");   // 18  7-18  Nucleic Acid Backbone
	RegisterAtName("CA2");  // 19  Shapely Special 
	RegisterAtName("SG");   // 20  Cysteine Sulfur 
	RegisterAtName("N1");   // 21
	RegisterAtName("N2");   // 22
	RegisterAtName("N3");   // 23
	RegisterAtName("N4");   // 24
	RegisterAtName("N6");   // 25
	RegisterAtName("O2");   // 26
	RegisterAtName("O4");   // 27
	RegisterAtName("O6");   // 28 21-28 Nucleic Acid H-Bonding
	RegisterAtName("HTM1");   // Terminal Haydrogens
	RegisterAtName("HTM2");   // Terminal Haydrogens
	RegisterAtName("HTM3");   // Terminal Haydrogens
	RegisterAtName("HTM4");   // Terminal Haydrogens
	RegisterAtName("HTM5");   // Terminal Haydrogens
}

bool HaAtom::IsSameName(const HaAtom* aptr_ref) const
{
	return( this->refno == aptr_ref->refno);
}

bool HaAtom::IsMatch(const HaAtom* atempl) const
{
	if( atempl == NULL ) return false;
	if( atempl->IsProxy() ) 
	{
		if( !atempl->replaced_atom_name.empty() )
		{
			if( stricmp_trunc( this->GetName(), replaced_atom_name.c_str() ) != 0 ) return false;
		}
	}
	return true;
}

// Flag Manipulation Utilities

void HaAtom::Select()    
{ 
	flag |= SelectFlag; 
}

void HaAtom::UnSelect()  
{ 
	flag &= ~SelectFlag; 
}

int HaAtom::Selected() const
{
	return (flag & SelectFlag);
}

void HaAtom::SelectBondedHydrogens()
{
	AtomGroup bonded_atoms = this->GetBondedAtoms();
	for (HaAtom* aptr : bonded_atoms)
		if (aptr->IsHydrogen()) aptr->Select();
}

bool HaAtom::IsDrawSphere() const
{
	return ( (flag & SphereFlag) != 0);
}

void HaAtom::SetDrawSphere(bool set_mode)
{
	if(set_mode)
		flag |= SphereFlag;
	else
		flag &= ~SphereFlag;
}

bool
HaAtom::IsDisplayed() const
{
	return ( (flag & DisplayedFlag) != 0);
}

void
HaAtom::SetDisplayed(bool set_mode)
{
	if(set_mode)
		flag |= DisplayedFlag;
	else
		flag &= ~DisplayedFlag;
}

bool
HaAtom::IsHydrogen() const
{
	return ( GetElemNo() == 1 );
}

bool
HaAtom::IsHBDonor() const
{
	return( (flag & HBDonorFlag) != 0 );
}

bool
HaAtom::IsHBAcceptor() const
{
	return( (flag & HBAcceptorFlag) != 0 );
}

void
HaAtom::SetHBDonor(bool set_mode)
{
	if(set_mode)
		flag |= HBDonorFlag;
	else
		flag &= ~HBDonorFlag;
}

void
HaAtom::SetHBAcceptor(bool set_mode)
{
	if(set_mode)
		flag |= HBAcceptorFlag;
	else
		flag &= ~HBAcceptorFlag;
}

bool HaAtom::IsAlphaCarbon() const
{
	return( refno == 1);
}

bool HaAtom::IsSugarPhosphate() const
{
	return( refno == 7);
}

bool HaAtom::IsAminoBackbone() const
{
	return( refno <=3 );
}

bool
HaAtom::IsShapelyBackbone() const
{
	return( refno <= 7);
}

bool
HaAtom::IsNucleicBackbone() const
{
	return( refno >=7 &&  refno<=18 );
}

bool
HaAtom::IsShapelySpecial() const
{
	return( refno ==19 );
}

bool
HaAtom::IsCysteineSulfur() const
{
	return((x)==20);
}

std::string HaAtom::GetStdSymbolElem(int elem)
{
    if(elem < 1 || elem > 109)
	{
       return "";
	}	
	return(StdAtmLbl[elem].c_str());
}

int HaAtom::GetElemNoFromName(const std::string& at_name, const HaResidue* pres)
{
    char ch1,ch2;
    const char *ptr;
	int is_coenzyme = FALSE;
	int is_protein  = FALSE;
	int is_nucleo   = FALSE;

    ptr = at_name.c_str();
	std::string res_name="";

	if(pres != NULL)
	{
		if(pres->IsCoenzyme())
			is_coenzyme = TRUE;
		if(pres->IsProtein())
			is_protein = TRUE;
		if(pres->IsNucleo())
			is_nucleo = TRUE;

		res_name = pres->GetName();
	}

	if(*ptr == ' ') ptr++;
	if(*ptr == ' ') ptr++;
	if( is_coenzyme )
	{   /* Exceptions to Brookhaven Atom Naming! */
		ch1 = ' ';
	}
	else
	{
		ch1 = ptr[0];
	}
	ch2 = ptr[1];
	
	ch1 = toupper(ch1);
	ch2 = toupper(ch2);

	// /* Handle HG, HD , CE, CD etc.. in Amino Acids! */
	if( is_protein )
	{
         if(ch1 != ' ')   ch2=' ';
	}
	else if( is_nucleo )
	{
		if(ch1 != ' ')   ch2=' ';
	}
	else if(!strcmp_trunc(res_name.c_str(),"RBP") )
	{
		if(ch1 != ' ' && ch1 != 'R') ch2=' ';
	}
	else if(!strcmp_trunc(res_name.c_str(),"HEM") )
	{
		if(ch1 != ' ' && ch1 != 'F') ch2=' ';
	}
	else if(!strcmp_trunc(res_name.c_str(),"FSF")  || !strcmp_trunc(res_name.c_str(),"FSA")
		    || !strcmp_trunc(res_name.c_str(),"FSB"))
	{
		if(ch1 != ' ' && ch1 != 'F') ch2=' ';
	}
	else if(!strcmp_trunc(res_name.c_str(),"PCL")  || !strcmp_trunc(res_name.c_str(),"PCA")
		    || !strcmp_trunc(res_name.c_str(),"PCB"))
	{
		if(ch1 != ' ' && ch1 != 'F') ch2=' ';
	}

	int elem;

    switch( ch1 )
    {
		case(' '): elem = GetElemNoFromChar(ch2);
			       return elem;
		
        case('A'):  switch( ch2 )
                    {   case('C'):  return( 89 );
                        case('G'):  return( 47 );
                        case('L'):  return( 13 );
                        case('M'):  return( 95 );
                        case('R'):  return( 18 );
                        case('S'):  return( 33 );
                        case('T'):  return( 85 );
                        case('U'):  return( 79 );
                    }
                    break;

        case('B'):  switch( ch2 )
                    {   case(' '):  return(  5 );
                        case('A'):  return( 56 );
                        case('E'):  return(  4 );
                        case('I'):  return( 83 );
                        case('K'):  return( 97 );
                        case('R'):  return( 35 );
                    }
                    break;

        case('C'):
			switch( ch2 )
			{
				case(' '):  return(  6 );
				case('A'):  return( 20 );
				case('D'):  return( 48 );
				case('E'):  return( 58 );
				case('F'):  return( 98 );
				case('L'):  return( 17 );
				case('M'):  return( 96 );
				case('O'):  return( 27 );
				case('R'):  return( 24 );
				case('S'):  return( 55 );
				case('U'):  return( 29 );
			}
			break;

		case('D'):
			if( ch2==' ' )
			{
				return(  1 );
			}
			else if( ch2=='Y' )
				return( 66 );
			break;

        case('E'):
			if( ch2=='R' )
			{
				return( 68 );
			}
			else if( ch2=='S' )
			{
				return( 99 );
			}
			else if( ch2=='U' )
				return( 63 );
			break;

        case('F'):
			if( ch2==' ' )
			{
				return(   9 );
			}
			else if( ch2=='E' )
			{
				return(  26 );
			}
			else if( ch2=='M' )
			{
				return( 100 );
			}
			else if( ch2=='R' )
				return(  87 );
			break;

        case('G'):
			if( ch2=='A' )
			{
				return( 31 );
			}
			else if( ch2=='D' )
			{
				return( 64 );
			}
			else if( ch2=='E' )
				return( 32 );
			break;

        case('H'):
			if( ch2==' ' )
			{
				return(  1 );
			}
			else if( ch2=='E' )
			{
				return(  2 );
			}
			else if( ch2=='F' )
			{
				return( 72 );
			}
			else if( ch2=='G' )
			{
				return( 80 );
			}
			else if( ch2=='O' )
				return( 67 );
			break;

        case('I'):
			if( ch2==' ' )
			{
				return( 53 );
			}
			else if( ch2=='N' )
			{
				return( 49 );
			}
			else if( ch2=='R' )
				return( 77 );
			break;

        case('K'):  if( ch2==' ' )
                    {   return( 19 );
                    } else if( ch2=='R' )
                        return( 36 );
                    break;

        case('L'):  if( ch2==' ' )
                    {   return(   1 );
                    } else if( ch2=='A' )
                    {   return(  57 );
                    } else if( ch2=='I' )
                    {   return(   3 );
                    } else if( (ch2=='R') || (ch2=='W') )
                    {   return( 103 );
                    } else if( ch2=='U' )
                        return(  71 );
                    break;

        case('M'):  if( ch2=='D' )
                    {   return( 101 );
                    } else if( ch2=='G' )
                    {   return(  12 );
                    } else if( ch2=='N' )
                    {   return(  25 );
                    } else if( ch2=='O' )
                        return(  42 );
                    break;

        case('N'):  switch( ch2 )
                    {   case(' '):  return(   7 );
                        case('A'):  return(  11 );
                        case('B'):  return(  41 );
                        case('D'):  return(  60 );
                        case('E'):  return(  10 );
                        case('I'):  return(  28 );
                        case('O'):  return( 102 );
                        case('P'):  return(  93 );
                    }
                    break;

        case('O'):  if( ch2==' ' )
                    {   return(  8 );
                    } else if( ch2=='S' )
                        return( 76 );
                    break;

        case('P'):  switch( ch2 )
                    {   case(' '):  return( 15 );
                        case('A'):  return( 91 );
                        case('B'):  return( 82 );
                        case('D'):  return( 46 );
                        case('M'):  return( 61 );
                        case('O'):  return( 84 );
                        case('R'):  return( 59 );
                        case('T'):  return( 78 );
                        case('U'):  return( 94 );
                    }
                    break;

        case('R'):  switch( ch2 )
                    {   case('A'):  return( 88 );
                        case('B'):  return( 37 );
                        case('E'):  return( 75 );
                        case('H'):  return( 45 );
                        case('N'):  return( 86 );
                        case('U'):  return( 44 );
                    }
                    break;

        case('S'):  switch( ch2 )
                    {   case(' '):  return( 16 );
                        case('B'):  return( 51 );
                        case('C'):  return( 21 );
                        case('E'):  return( 34 );
                        case('I'):  return( 14 );
                        case('M'):  return( 62 );
                        case('N'):  return( 50 );
                        case('R'):  return( 38 );
                    }
                    break;

        case('T'):  switch( ch2 )
                    {   case('A'):  return( 73 );
                        case('B'):  return( 65 );
                        case('C'):  return( 43 );
                        case('E'):  return( 52 );
                        case('H'):  return( 90 );
                        case('I'):  return( 22 );
                        case('L'):  return( 81 );
                        case('M'):  return( 69 );
                    }
                    break;

        case('U'):  if( ch2==' ' )
                        return( 92 );
                    break;

        case('V'):  if( ch2==' ' )
                        return( 23 );
                    break;

        case('W'):  if( ch2==' ' )
                        return( 74 );
                    break;

        case('X'):  if( ch2=='E' )
                        return( 54 );
                    break;

        case('Y'):  if( ch2==' ' )
                    {   return( 39 );
                    } else if( ch2=='B' )
                        return( 70 );
                    break;

        case('Z'):  if( ch2=='N' )
                    {   return( 30 );
                    } else if( ch2=='R' )
                        return( 40 );
                    break;
    }

    if( (ch1>='0') && (ch1<='9') )
        if( (ch2=='H') || (ch2=='D') )
            return( 1 ); /* Hydrogen */

    /* If all else fails! */
    return GetElemNoFromChar(ch1);
}

int HaAtom::GetElemNoFromChar(char ch)
{
	switch( ch )
    {
	    case('B'):  return(  5 );
        case('C'):  return(  6 );
        case('D'):  return(  1 );
        case('F'):  return(  9 );
        case('H'):  return(  1 );
        case('I'):  return( 53 );
        case('K'):  return( 19 );
        case('L'):  return(  1 );
        case('N'):  return(  7 );
        case('O'):  return(  8 );
        case('P'):  return( 15 );
        case('S'):  return( 16 );
        case('U'):  return( 92 );
        case('V'):  return( 23 );
        case('W'):  return( 74 );
        case('Y'):  return( 39 );
    }
	return 0;
}



static const double VDWCarbon   = 1.87;
static const double VDWNitrogen = 1.5;
static const double VDWOxygen   = 1.4;
static const double VDWSulfur  = 1.85;


double HaAtom::ElemVDWRadius(int elem, bool united_atom_flag)
{
    if( united_atom_flag )
        switch( elem )
        {
			case(  6 ):  return( VDWCarbon );
            case(  7 ):  return( VDWNitrogen );
            case(  8 ):  return( VDWOxygen );
            case( 16 ):  return( VDWSulfur );
        }
    return( ElementArr[elem].vdwrad );
}

double HaAtom::ElemDuttonRadius(int elem)
{
        switch( elem )
        {
			case(  1 ):  return( 0.0 );
			case(  6 ):  return( 2.0 );
            case(  7 ):  return( 2.0 );
            case(  8 ):  return( 1.6 );
			case( 15 ):  return( 1.9 );
            case( 16 ):  return( 1.8 );
			case( 26 ):  return( 2.8 );
        }

	    return( ElementArr[elem].vdwrad );
}


double HaAtom::StdElemMass(int elem)
{
	if( elem < 0 || elem > 103 ) 
	{
		return 1.0;
	}
	return( ElementArr[elem].mass );
}

double HaAtom::GetStdMass() const
{
	return StdElemMass(this->GetElemNo() );
}

int HaAtom::SetHybrid(const std::string& hybrid_str_par )
{
	std::string hybrid_str = hybrid_str_par;
	boost::trim(hybrid_str);
	boost::to_upper(hybrid_str);

	if(hybrid_str == "SP" )
	{
		hybrid = SP_HYBRID;
		return TRUE;
	}
	else if(hybrid_str == "SP2" )  
	{
		hybrid = SP2_HYBRID;
		return TRUE;
	}
	else if(hybrid_str == "SP3" )  
	{
		hybrid = SP3_HYBRID;
		return TRUE;
	}
	hybrid = NO_HYBRID;
	return FALSE;
}

int HaAtom::GetHybridTextStr(std::string& hybrid_str) const
{
	if(hybrid == NO_HYBRID)
	{
		hybrid_str = "";
	}
	else if( hybrid == SP_HYBRID)
	{
		hybrid_str = "SP";
	}
	else if( hybrid == SP2_HYBRID)
	{
		hybrid_str = "SP2";
	}
	else if( hybrid == SP3_HYBRID)
	{
		hybrid_str = "SP3";
	}
	return TRUE;
}

int HaAtom::SetHBStatus(const char* hb_status_c_str)
{
	std::string hb_status_str(hb_status_c_str);
	boost::trim(hb_status_str);
	boost::to_upper(hb_status_str);
	
	SetHBDonor(false);
	SetHBAcceptor(false);
	int i;
	for(i = 0; i < hb_status_str.size(); i++)
	{
		if( hb_status_str[i] == 'D')
			SetHBDonor(true);
		else if( hb_status_str[i] == 'A')
			SetHBAcceptor(true);
	}
	return TRUE;
}

int HaAtom::GetHBStatusTextStr(std::string& hb_status_str)
{
	hb_status_str = "";
	if( IsHBDonor() )
		hb_status_str+= "D";
	if( IsHBAcceptor() )
		hb_status_str+= "A";
	return TRUE;
}

double HaAtom::StdBondLen(const HaAtom* aptr1, const HaAtom* aptr2)
{
	if(aptr1 == NULL || aptr2 == NULL)
	{
		ErrorInMod( "HaAtom::StdBondLen()", 
			        "Atom pointer is NULL");
		return 0.0;
	}
	int elem1 = aptr1->GetElemNo();
	int elem2 = aptr2->GetElemNo();

	HYBRIDIZATION hb1 = aptr1->GetHybrid();
    HYBRIDIZATION hb2 = aptr2->GetHybrid();

	if(hb1 == NO_HYBRID) hb1 = SP3_HYBRID;
    if(hb2 == NO_HYBRID) hb2 = SP3_HYBRID;

// From J. March "Advanced Organic Chemistry":
// Huheey "Inorganic Chemistry"

	if( elem1 > elem2)
	{
		int axx = elem1;
		elem1 = elem2;
		elem2 = axx;
		HYBRIDIZATION haxx = hb1;
		hb1 = hb2;
		hb2 = hb1;
	}

	double len = 1.5;

	switch(elem1)
	{
	   case 1: 
		   if(elem2 == 6)
		   {
			  if(hb2 == SP3_HYBRID)       len = 1.09;
			  else if( hb2 == SP2_HYBRID) len = 1.08;
			  else if( hb2 == SP_HYBRID)  len = 1.08;
		   }
		   else if( elem2 == 7)    len = 1.01;   // N - Huheey
		   else if( elem2 == 8)    len = 0.96;   // O - Huheey
           else if( elem2 == 9)    len = 0.918;  // F - Huheey 
           else if( elem2 == 17)   len = 1.274;  // Cl - Huheey 
		   else if( elem2 == 35)   len = 1.408;  // Br - Huheey 
		   else if( elem2 == 53)   len = 1.608;  // I  - Huheey 

		   break;
	   case 6:
		   if(elem2 == 6)
		   {
			  if(hb1 == SP3_HYBRID)
			  {
				  if(hb2 == SP3_HYBRID)       len = 1.53; // March // 1.54 - Huheey
				  else if( hb2 == SP2_HYBRID) len = 1.51; // March // 1.34 - Huheey
	              else if( hb2 == SP_HYBRID)  len = 1.47; // March // 1.20 - Huheey
			  }
			  else if( hb1 == SP2_HYBRID)
			  {
				  if(hb2 == SP3_HYBRID)       len = 1.51;
				  else if( hb2 == SP2_HYBRID) len = 1.48;
	              else if( hb2 == SP_HYBRID)  len = 1.43;
			  }
			  else if( hb1 == SP_HYBRID)
			  {
				  if(hb2 == SP3_HYBRID)       len = 1.47;
				  else if( hb2 == SP2_HYBRID) len = 1.43;
	              else if( hb2 == SP_HYBRID)  len = 1.38;
			  }
		   }
		   else if( elem2 == 7) // N
		   {
			  if(hb1 == SP3_HYBRID)          len = 1.47;
			  else if( hb1 == SP2_HYBRID)    len = 1.38;
		   }
		   else if( elem2 == 8) // O
		   {
			  if(hb1 == SP3_HYBRID)          len = 1.43;
			  else if( hb1 == SP2_HYBRID)    len = 1.34;
		   }
		   else if( elem2 == 9) // F
		   {
			  if(hb1 == SP3_HYBRID)          len = 1.40;
			  else if( hb1 == SP2_HYBRID)    len = 1.34;
              else if( hb1 == SP_HYBRID)     len = 1.27;
		   }
		   else if( elem2 == 16) // S 
		   {
			  if(hb1 == SP3_HYBRID)          len = 1.82;
			  else if( hb1 == SP2_HYBRID)    len = 1.75;
			  else if( hb1 == SP_HYBRID)     len = 1.68;
		   }
		   else if( elem2 == 17) // Cl 
		   {
			  if(hb1 == SP3_HYBRID)          len = 1.79;
			  else if( hb1 == SP2_HYBRID)    len = 1.73;
			  else if( hb1 == SP_HYBRID)     len = 1.63;
		   }
		   else if( elem2 == 35) // Br 
		   {
			  if(hb1 == SP3_HYBRID)          len = 1.97;
			  else if( hb1 == SP2_HYBRID)    len = 1.88;
			  else if( hb1 == SP_HYBRID)     len = 1.79;
		   }
		   else if( elem2 == 53) // I 
		   {
			  if(hb1 == SP3_HYBRID)          len = 2.16;
			  else if( hb1 == SP2_HYBRID)    len = 2.10;
			  else if( hb1 == SP_HYBRID)     len = 1.99;
		   }
		   break;
	   case 7:
           if( elem2 == 7)          len = 1.45; // N-N - Huheey
		   else if( elem2 == 8)     len = 1.40; // O   - Huheey
		   else if( elem2 == 9)     len = 1.36; // F   - Huheey
		   else if( elem2 == 17)    len = 1.75; // Cl  - Huheey
           else if( elem2 == 26)    len = 2.00; // Fe  - me 
           else if( elem2 == 44)    len = 2.00; // Ru  - me 
		   else if( elem2 == 76)    len = 2.00; // Os  - me 
           break;
	}
	return len;
}

bool HaAtom::SetCoordSubstH(const HaAtom* aptr1, const HaAtom* aptr2,
		                   HaAtom* haptr)
{
	double h_dist=1.08;
	
	if( aptr1 == NULL || aptr2 == NULL || haptr == NULL )
	{
		ErrorInMod("HaAtom::SetCoordSubstH()",
		  "Some of atom pointers are NULL ");
		return false;
	}

	double x1=aptr1->GetX();
	double y1=aptr1->GetY();
	double z1=aptr1->GetZ();
	double x2=aptr2->GetX();
	double y2=aptr2->GetY();
	double z2=aptr2->GetZ();

	double ss=sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
	if(ss < 0.000001)
	{
		PrintLog(" HaAtom::GetCoordSubstH \n");
		PrintLog(" distance between atoms is equal zero \n");
		return false;
	}
	haptr->SetX( x1+ h_dist*(x2-x1)/ss);
	haptr->SetY( y1+ h_dist*(y2-y1)/ss);
	haptr->SetZ( z1+ h_dist*(z2-z1)/ss);

	haptr->SetElemNo(1);
	std::string at_name = "H";
	haptr->SetName(at_name.c_str());
	
	return true;

}

double HaAtom::GetMass() const 
{ 
	return ps_ff_par->mass;
}              

void HaAtom::SetMass(double new_mass) 
{ 
	ps_ff_par->mass = new_mass;
}

double HaAtom::GetVdWRad() const
{
	return ps_ff_par->rad_vdw;
}

void HaAtom::SetVdWRad(double new_vdw_rad)
{
	ps_ff_par->rad_vdw = new_vdw_rad;
}

double HaAtom::GetVdWEne() const
{
	return ps_ff_par->ene_vdw;
}

void HaAtom::SetVdWEne(double new_ene_vdw)
{
	ps_ff_par->ene_vdw = new_ene_vdw;
}

const std::string HaAtom::GetFFSymbol() const 
{
	if (ps_ff_par) return ps_ff_par->ff_symbol;
	return {};
} 

void HaAtom::SetFFSymbol( const std::string& new_ff_symbol)
{
	if (!ps_ff_par) ps_ff_par = std::make_shared<AtomFFParam>();
	ps_ff_par->ff_symbol = new_ff_symbol; 
}

void HaAtom::SetCharge(double new_charge)
{ 
	if (!ps_ff_par) ps_ff_par = std::make_shared<AtomFFParam>();
	ps_ff_par->charge = new_charge;
} 

double HaAtom::GetCharge() const
{ 
	if (ps_ff_par) return ps_ff_par->GetCharge();
	return 0.0;
} 

TiXmlElement* HaAtom::AddXml(TiXmlElement* parent_element, const char* name, int option) const
{
	if( parent_element == NULL) return NULL;
	TiXmlElement* atom_element;
	
	if(strlen(name) == 0)
	{
		atom_element = new TiXmlElement("HaAtom");
	}
	else
	{
		atom_element = new TiXmlElement(name);
		atom_element->SetAttribute("TYPE","HaAtom");
	}	
	atom_element->SetAttribute("elemno", elemno);
	atom_element->SetAttribute("name", this->GetName());
	atom_element->SetDoubleAttribute("x", this->GetX());
	atom_element->SetDoubleAttribute("y", this->GetY());
	atom_element->SetDoubleAttribute("z", this->GetZ());
    atom_element->SetAttribute("idi", reinterpret_cast<std::uintptr_t>(this) );

	parent_element->LinkEndChild(atom_element);

	return atom_element;
}

int HaAtom::LoadXml(const TiXmlElement* xml_element, int option)
{
	try
	{
		if( xml_element == NULL) throw std::runtime_error("xml_element == NULL");
  
		int itmp, ires;
		double dtmp;

		HaResidue* pres = this->GetHostRes(); 

		std::string type = "normal";
		if( xml_element->CStrAttribute("type") ) type = xml_element->CStrAttribute("type");
		if( type == "proxy" ) this->SetProxy();
		if( xml_element->CStrAttribute("replaced_at") )
		{
			std::string repl_atname = xml_element->CStrAttribute("replaced_at");
			this->SetReplacedAtName(repl_atname);
		}

		if( xml_element->CStrAttribute("bond_to") ) 
		{
			std::vector<std::string> str_vec;
			std::string text = xml_element->CStrAttribute("bond_to");
			boost::split(str_vec,text,boost::is_any_of("; "),boost::token_compress_on);
			int i;

			for(i = 0; i < str_vec.size(); i++ )
			{
				std::string bat_name = str_vec[i];
				HaAtom* bat = pres->GetAtomByName( bat_name.c_str() );
				if( pres && bat && bat != this )
				{
					HaAtom::CreateBond(this,bat);
				}
			}
		}

		if( xml_element->CStrAttribute("desc") ) this->SetDescription( xml_element->CStrAttribute("desc") );

		ires = xml_element->QueryIntAttribute("elemno",&itmp);
		if( ires == TIXML_SUCCESS ) this->SetElemNo(itmp);

		if( xml_element->CStrAttribute("name") ) this->SetName( xml_element->CStrAttribute("name") );

		ires = xml_element->QueryDoubleAttribute("x",&dtmp);
		if(ires == 0) this->SetX(dtmp);
		ires = xml_element->QueryDoubleAttribute("y",&dtmp);
		if(ires == 0) this->SetY(dtmp);
		ires = xml_element->QueryDoubleAttribute("z",&dtmp);
		if(ires == 0) this->SetZ(dtmp);

		const TiXmlElement* data_element = xml_element->FirstChildElement();
		while( data_element )
		{
			std::string data_elem_tag = data_element->Value();
			if( data_elem_tag == "intcrd") // read internal coordinates
			{
				std::vector<std::string> str_arr;
				std::string int_crd_str = data_element->GetText();
				std::istrstream is( int_crd_str.c_str()  );
				std::string atn1,atn2,atn3;
				double dist, angle, torsion;
				
				is >> atn1;
				is >> dist;
				is >> atn2;
				is >> angle;
				is >> atn3;
				is >> torsion;
				if( is.fail() ) throw std::runtime_error(" error reading internal coordinates: " + int_crd_str);
				HaAtom* aptr1 = pres->GetAtomByName(atn1.c_str());
				if( aptr1 == NULL) throw std::runtime_error(" error reading internal coordinates: not found atom " + atn1 );
				HaAtom* aptr2 = pres->GetAtomByName(atn2.c_str());
				if( aptr2 == NULL) throw std::runtime_error(" error reading internal coordinates: not found atom " + atn2 );
				HaAtom* aptr3 = pres->GetAtomByName(atn3.c_str());
				if( aptr3 == NULL) throw std::runtime_error(" error reading internal coordinates: not found atom " + atn3 );
				int ires = Vec3D::SetAtomPos(this,aptr1,aptr2,aptr3,dist,angle,torsion);
				if(!ires ) throw std::runtime_error(" error reading internal coordinates: not able to set atom coordinates ");
			}
			data_element = data_element->NextSiblingElement();
		}
	}
	catch(std::exception& ex)
	{
		PrintLog(" Error in HaAtom::LoadXml() \n");
		PrintLog(" Loading info for atom %s\n",GetRef().c_str() );
		return FALSE;
	}
	return TRUE;
}


AtomDoubleMap::AtomDoubleMap(const char* new_name)
{
	SetName(new_name);
}

AtomDoubleMap::~AtomDoubleMap()
{

}

double AtomDoubleMap::GetValue(HaAtom* aptr) const
{
	map<HaAtom*,double>::const_iterator itr;
	itr = find(aptr);
	if( itr == this->end()) return 0.0;

	return (*itr).second;
}

int AtomDoubleMap::SetValue(HaAtom* aptr, double new_val)
{
	map<HaAtom*,double>::iterator itr;
	itr = find(aptr);
	if( itr == this->end())
	{
		this->insert(AtomDoubleMap::value_type(aptr,new_val));
		return True;
	}
	
	(*itr).second = new_val;

	return True;
}

AtomAtomMap::AtomAtomMap()
{

}

AtomAtomMap::~AtomAtomMap()
{

}

HaAtom* AtomAtomMap::GetValue(HaAtom* aptr)
{
	return this->at(aptr);
}


void AtomAtomMap::SetValue(HaAtom* aptr, HaAtom* val)
{
	(*this)[aptr] = val;
}
