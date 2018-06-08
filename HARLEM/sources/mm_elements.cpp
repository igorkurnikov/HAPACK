/*!  \file mm_elements.cpp

    Molecular Mechanics Model elements 

    \author Igor Kurnikov 
    \date 2010-

*/

#define MM_ELEMENTS_CPP

#include "mpi.h"

#include <float.h>
#include <math.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "tinyxml.h"
#include "hamolmech.h"
#include "mm_elements.h"
#include "mm_model.h"

AtomFFParam::AtomFFParam()
{
	Clear();
}

AtomFFParam::AtomFFParam(HaAtom* aptr_ref_new)
{
	Clear();
	aptr_ref = aptr_ref_new;
	if( aptr_ref != NULL )  at_name = aptr_ref->GetName();
}

AtomFFParam::~AtomFFParam()
{
	Clear();
}

void AtomFFParam::Clear()
{
	polar  = 0.0;
	hpolar  = 0.0;
	charge = 0.0;
	dipole.clear();
	qpole.clear();
	bisect_flag = FALSE;
	at_name.clear();
	ff_symbol.clear();
	ff_polar_symbol.clear();
	frame_atom_names.clear();
	frame_atoms.clear();

	damp_polar_strength = 0.0;
	damp_polar_sensitivity = 0.0;
	damp_polar_rad = 0.0;

	screen_polar = 0.0;
	
	aptr_ref = NULL;
}

int AtomFFParam::LoadXml(const TiXmlElement* xml_element, int option )
{
	char buf[256];
	if( xml_element == NULL) return FALSE;
	std::string elem_name = xml_element->ValueStr();
	try
	{
		int ival;
		double dval;

		if( elem_name == "dummy_at" )
		{
			if( xml_element->CStrAttribute("name") ) ff_polar_symbol = xml_element->CStrAttribute("pol_type_id");
		}

		if( aptr_ref == NULL ) throw std::runtime_error(" Atom Reference pointer is NULL");
		if( xml_element->QueryDoubleAttribute("chg",&dval ) == TIXML_SUCCESS )
		{
			charge = dval;
		}

		if( xml_element->CStrAttribute("pol_type_id") ) ff_polar_symbol = xml_element->CStrAttribute("pol_type_id");
		if( xml_element->CStrAttribute("ff_type") )     ff_symbol       = xml_element->CStrAttribute("ff_type");

		const TiXmlElement* data_element;

		data_element = xml_element->FirstChildElement();
		for( ; data_element;  data_element = data_element->NextSiblingElement())
		{
			std::string name = data_element->ValueStr();
			std::string text;
			std::vector<std::string> str_vec;
			if( name == "frame" )
			{
				if( data_element->CStrAttribute("bisect") )
				{
					std::string bisect_val = data_element->CStrAttribute("bisect");
					if( bisect_val == "1" ) bisect_flag = TRUE;
				}
				text = data_element->GetText();
				boost::split(str_vec,text,boost::is_any_of(" "),boost::token_compress_on);
				if( str_vec.size() == 0 ) throw std::runtime_error("Empty frame atom names list");
				frame_atom_names = str_vec;
				
			}
			else if( name == "dipole" )
			{
				text = data_element->GetText();
				boost::split(str_vec,text,boost::is_any_of(" "),boost::token_compress_on);
				if( str_vec.size() != 3 ) throw std::runtime_error("The number of dipole components is not equal to 3");
				dipole.resize(3);
				for( int i = 0; i < 3; i++ )
				{
					try 
					{
						dipole[i] = boost::lexical_cast<double>(str_vec[i]);
					} 
					catch( boost::bad_lexical_cast& )
					{
						throw std::runtime_error(" Error to convert dipole component string to double ");
					}
				}
			}
			else if( name == "qpole")
			{
				text = data_element->GetText();
				boost::split(str_vec,text,boost::is_any_of(" "),boost::token_compress_on);
				if( str_vec.size() != 6 ) throw std::runtime_error("The number of quadrupole components is not equal to 6");
				qpole.resize(6);
				for( int i = 0; i < 6; i++ )
				{
					try 
					{
						qpole[i] = boost::lexical_cast<double>(str_vec[i]);
					} 
					catch( boost::bad_lexical_cast& )
					{
						throw std::runtime_error(" Error to convert quadrupole component string to double ");
					}
				}
			}
			else if( name == "polar" )
			{
				text = data_element->GetText();
				boost::split(str_vec,text,boost::is_any_of(" "),boost::token_compress_on);
				polar = 0.0;
				if( str_vec.size() != 1 ) throw std::runtime_error("The number of polarization components is not equal to 1");
				try 
				{
					polar = boost::lexical_cast<double>(str_vec[0]);
				} 
				catch( boost::bad_lexical_cast& )
				{
					throw std::runtime_error(" Error to convert polarization component string to double ");
				}
			}
			else if( name == "hpolar" )
			{
				text = data_element->GetText();
				boost::split(str_vec,text,boost::is_any_of(" "),boost::token_compress_on);
				hpolar = 0.0;
				if( str_vec.size() != 1 ) throw std::runtime_error("The number of hyperpolarization components is not equal to 1");
				try 
				{
					hpolar = boost::lexical_cast<double>(str_vec[0]);
				} 
				catch( boost::bad_lexical_cast& )
				{
					throw std::runtime_error(" Error to convert hyperpolarization component string to double ");
				}
			}
			else if( name == "damp_polar" )
			{
				text = data_element->GetText();
				boost::split(str_vec,text,boost::is_any_of(" "),boost::token_compress_on);
				damp_polar_strength = 0.0;
				damp_polar_sensitivity = 0.0;
				damp_polar_rad = 0.0;
				if( str_vec.size() != 3 ) throw std::runtime_error("The number of damping polarization components is not equal to 3");
				try 
				{
					damp_polar_strength    = boost::lexical_cast<double>(str_vec[0]);
					damp_polar_sensitivity = boost::lexical_cast<double>(str_vec[1]);
					damp_polar_rad         = boost::lexical_cast<double>(str_vec[2]);
				} 
				catch( boost::bad_lexical_cast& )
				{
					throw std::runtime_error(" Error to convert damping polarization component string to double ");
				}
			}
			else if( name == "screen_polar" )
			{
				text = data_element->GetText();
				boost::split(str_vec,text,boost::is_any_of(" "),boost::token_compress_on);
				screen_polar = 0.0;
				if( str_vec.size() != 1 ) throw std::runtime_error("The number of polarization screening components is not equal to 1");
				try 
				{
					screen_polar = boost::lexical_cast<double>(str_vec[0]);
				} 
				catch( boost::bad_lexical_cast& )
				{
					throw std::runtime_error(" Error to convert polarization screen string to double ");
				}
			}
		}
	}
	catch( const std::exception& ex )
	{
		PrintLog("Error in AtomFFParam::LoadXml atom template: %s \n", at_name.c_str() );
		PrintLog("Error: %s \n",ex.what());
		return FALSE;
	}
	return TRUE;

}

int AtomFFParam::HasDipole() const
{
	return (dipole.size() > 2);
}
	
int AtomFFParam::HasQPole()  const
{
	return (qpole.size() > 0);
}
	
int AtomFFParam::HasPolar()  const
{
	return (polar > 0.0);	
};    
	
int AtomFFParam::HasFrameAtomNames() const
{
	return (frame_atom_names.size() > 0);
}
	
int AtomFFParam::IsFrameSet()  const
{
	return frame_atoms.size() > 0;
}

bool AtomFFParam::AlterPolarizability() const
{
	return ( damp_polar_strength > DBL_EPSILON ); 
}

bool AtomFFParam::IsPolarPositionDep() const
{
	return ( damp_polar_sensitivity > DBL_EPSILON ); 
}

bool AtomFFParam::HasScreenPolar() const
{
	return ( screen_polar > DBL_EPSILON ); 
}

bool AtomFFParam::HasHPolar() const
{
	return ( hpolar > DBL_EPSILON ); 
}

double AtomFFParam::GetCharge() const
{
	return charge;
}
	
HaVec_double AtomFFParam::GetDipole() const
{
	return dipole;
}
	
HaVec_double AtomFFParam::GetQPole() const
{
	return qpole;
}
	
int AtomFFParam::IsChiralFrame() const
{
	return (frame_atom_names.size() == 4);
}

int AtomFFParam::IsBisectFrame() const
{
	return bisect_flag;
}
	
int AtomFFParam::SetFrameFromAtomNames()
{
	char buf[256];
	try
	{
		if( aptr_ref == NULL ) throw std::runtime_error(" Pointer to the reference atom is NULL \n");
		HaResidue* pres = aptr_ref->GetHostRes();
		if( pres == NULL ) throw std::runtime_error(" Pointer to reference residue is NULL \n");
		if( frame_atom_names.size() < 3 ) throw std::runtime_error( " There are less than 3 atom names specified \n ");
		frame_atoms.clear();
		int na = frame_atom_names.size();
		int i;
		for( i = 0; i < na; i++)
		{
			std::string at_frame_name = frame_atom_names[i];
			HaAtom* aptr = pres->GetAtomByName(at_frame_name.c_str());
			if( aptr == NULL ) 
			{
				sprintf(buf," No atom with name %s in the residue \n ", at_frame_name.c_str());
				throw std::runtime_error(buf);
			}
			frame_atoms.push_back(aptr);
		}
	}
	catch( std::exception& ex )
	{
		PrintLog(" Error in AtomFFParam::SetFrameFromAtomNames()  Atom: %s \n",at_name.c_str() );
		PrintLog(" Error: %s \n",ex.what());
		return FALSE;
	}
	return TRUE;
}
 
ResFFTemplate::ResFFTemplate(HaResidue* p_res_templ_new)
{
	p_res_templ = p_res_templ_new;
	if( p_res_templ != NULL ) res_name = p_res_templ->GetFullName(); 
}

ResFFTemplate::~ResFFTemplate()
{
	Clear();
}

void ResFFTemplate::Clear()
{
	if( atom_params.size() > 0 )
	{
		int n = atom_params.size();
		int i;
		for( i = 0; i < n; i++ )
		{
			delete atom_params[i];
		}
		atom_params.clear();
		at_name_ff_param_map.clear();
	}

	res_name.clear();
	res_ff_version.clear();
	p_res_templ = NULL;
}

int ResFFTemplate::SetResFFVersion(const std::string& res_ff_version_new)
{
	res_ff_version = res_ff_version_new;
	return TRUE;
}

std::string ResFFTemplate::GetFullName()
{
	std::string full_name = res_name;
	if( !res_ff_version.empty() ) 
	{
		full_name += ":";
		full_name += res_ff_version;
	}
	return res_ff_version;
}

AtomFFParam* ResFFTemplate::GetAtomFFParam(const std::string& at_name)
{
	std::string at_name_loc = at_name;
	boost::trim(at_name_loc);
	boost::to_upper(at_name_loc);

	if( at_name_ff_param_map.count(at_name_loc) == 0) return NULL;

	return at_name_ff_param_map[at_name_loc];
}

int ResFFTemplate::SetAtomFFParam(const std::string& at_name, AtomFFParam* p_at_ff_param)
{
	try
	{
		if( p_res_templ == NULL ) throw std::runtime_error(" p_res_templ == NULL ");
		if( p_res_templ->GetAtomByName(at_name.c_str()) == NULL) throw std::runtime_error(" No atom with this name in residue template ");
		if( at_name_ff_param_map.count(at_name) > 0 ) throw std::runtime_error(" Atom FF Params already exist for this atom name ");
		
		atom_params.push_back( p_at_ff_param );
		at_name_ff_param_map[at_name] = p_at_ff_param;
	}
	catch( std::exception& ex )
	{
		PrintLog(" Error Setting Atom FF param %s for residue FF template: %s \n",
			       at_name.c_str(),GetFullName().c_str());
		PrintLog(" Error: %s \n",ex.what());
		return FALSE;
	}
	return TRUE;
}

int ResFFTemplate::LoadXml(const TiXmlElement* xml_element, int option )
{
	char buf[256];
	if( xml_element == NULL) return FALSE;
	try
	{
		if( p_res_templ == NULL ) throw std::runtime_error(" Residue Template pointer is NULL");
		const TiXmlElement* data_element;
		data_element = xml_element->FirstChildElement();
		for( ; data_element;  data_element = data_element->NextSiblingElement() )
		{
			std::string element_type = data_element->ValueStr();
			if(element_type == "atomff")
			{
				std::string at_name = data_element->Attribute("name");
				HaAtom* aptr = p_res_templ->GetAtomByName(at_name.c_str());
				if( aptr == NULL)
				{
					sprintf(buf,"No atom name %s in the residue %s \n",
					        at_name.c_str(), p_res_templ->GetFullName().c_str());
					throw std::runtime_error(buf);
				}
				AtomFFParam* p_at_ff = new AtomFFParam(aptr);
				int ires = p_at_ff->LoadXml(data_element);
				if( !ires )
				{
					delete p_at_ff;
					sprintf(buf,"Error Loading FF parameters for atom %s\n",at_name.c_str());
					throw std::runtime_error(buf);
				}
				if( at_name_ff_param_map.count(at_name) > 0 )
				{
					delete p_at_ff;
					sprintf(buf,"Atom FF parameters for atom %s are already set\n",at_name.c_str());
					throw std::runtime_error(buf);
				}
				atom_params.push_back(p_at_ff);
				at_name_ff_param_map[at_name] = p_at_ff;
			}
			if(element_type == "improper_dihedrals")
			{ 
				std::string text = data_element->GetText();
				std::vector<std::string> str_vec;
				boost::split(str_vec,text,boost::is_any_of(" "),boost::token_compress_on);
				
				if( str_vec.size() % 4 != 0 ) throw std::runtime_error("The number of Improper angle atom symbols is not mutiple of 4");
				int nimpr = str_vec.size()/4;
				int impr;
				for( impr = 0; impr < nimpr; impr++ )
				{
					std::vector<std::string> impr_at_str;
					for(int j = 0; j < 4; j++)
					{
						impr_at_str.push_back( str_vec[4*impr+j] );
					}
					improper_dihedrals.push_back(impr_at_str);
				}
			}
		}
	}
	catch( std::exception& ex )
	{
		PrintLog("Error in ResFFTemplate::LoadXml() residue template: %s  version: %s\n", 
			     res_name.c_str(),this->res_ff_version.c_str());
		PrintLog("Error: %s \n",ex.what());
		return FALSE;
	}
	return TRUE;
}

HaResidue* ResFFTemplate::GetResTemplate()
{
	return p_res_templ;
}

MMBond::MMBond()
{
  pt1 = pt2 = NULL;
  r0 = 0.0;
  fc = 0.0;
  set_type = MolMechModel::NOT_SET;
}

MMBond::MMBond(HaAtom* new_pt1, HaAtom* new_pt2)
{
	if(new_pt1 == NULL  ||  new_pt2 == NULL )
	{
		ErrorInMod("MMBond::MMBond()", 
		           " One of the HaAtom pointers is NULL ");
		pt1 = pt2 = NULL;
		return;
	}
	if(new_pt1 < new_pt2)
	{
		pt1= new_pt1;
		pt2= new_pt2;
	}
	else
	{
		pt1 = new_pt2;
		pt2 = new_pt1;
	}
	r0 = 0.0;
	fc = 0.0;
	set_type = MolMechModel::NOT_SET;
}

MMBond::MMBond(const MMBond& bnd_ref)
{
	pt1 = bnd_ref.pt1;
	pt2 = bnd_ref.pt2;
	r0  = bnd_ref.r0;
	fc  = bnd_ref.fc;
	set_type = bnd_ref.set_type;
}

MMBond::~MMBond()
{

}

bool MMBond::operator==(const MMBond& rhs) const
{
	return( pt1 == rhs.pt1 && pt2 == rhs.pt2);
}

bool MMBond::operator < (const MMBond& rhs) const
{
	if( pt1 < rhs.pt1) 
		return true;
	if( pt1 > rhs.pt1)
		return false;
	return( pt2 < rhs.pt2);
}


MMValAngle::MMValAngle()
{
	pt1 = NULL;
	pt2 = NULL;
	pt3 = NULL;
	a0 = 0.0;
	fc = 0.0;
	set_type = MolMechModel::NOT_SET;
}



MMValAngle::MMValAngle(HaAtom* new_pt1, HaAtom* new_pt2, HaAtom* new_pt3):
pt1(new_pt1),pt2(new_pt2),pt3(new_pt3)
{
	if(new_pt1 == NULL || new_pt2 == NULL  || new_pt3 == NULL )
	{
		ErrorInMod("MMValAngle::MMValAngle()", 
		           " One of the HaAtom pointers is NULL "); 
		pt1 = pt2 = pt3 = NULL;
	}
	a0 = 0.0;
	fc = 0.0;
	set_type = MolMechModel::NOT_SET;
}
	
MMValAngle::~MMValAngle()
{

}

bool MMValAngle::operator < (const MMValAngle & rhs) const
{
	if(pt1 > rhs.pt1) return false;
	if(pt1 < rhs.pt1) return true;
	if(pt2 > rhs.pt2) return false;
	if(pt2 < rhs.pt2) return true;
	if(pt3 > rhs.pt3) return false;
	if(pt3 < rhs.pt3) return true;
	return false;
}

bool MMValAngle::operator == (const MMValAngle & rhs) const
{
	if(pt1 != rhs.pt1) return false;
	if(pt2 != rhs.pt2) return false;
	if(pt3 != rhs.pt3) return false;
	return true;
}

int MMDihedral::AddTerm(double pn_new, double phase_new, double pk_new, double idivf_new)
{
	if( fabs( idivf_new ) < DBL_EPSILON )
	{
		ErrorInMod("MMDihedral::AddTerm()"," idivf == 0 ");
		return False;
	}

	pn.push_back(pn_new);
	phase.push_back(phase_new);
	pk.push_back(pk_new/idivf_new);
	idivf.push_back(idivf_new);
	return TRUE;
}


MMDihedral::MMDihedral()
{
  pt1 = pt2 = pt3 = pt4 = NULL;

  improper = false;
  calc_14 = true;
  set_type = MolMechModel::NOT_SET;
}

MMDihedral::MMDihedral(HaAtom* new_pt1, HaAtom* new_pt2, 
					   HaAtom* new_pt3, HaAtom* new_pt4, bool improper_flag)
{
	if(new_pt1 == NULL || new_pt2 == NULL ||
	   new_pt3 == NULL || new_pt4 == NULL ) 
	   {

		ErrorInMod("MMDihedral::MMDihedral()", 
		           " One of the HaAtom pointers is NULL ");
		return;
	}
	pt1= new_pt1;
	pt2= new_pt2;
	pt3= new_pt3;
	pt4= new_pt4;
	
	improper = improper_flag;
	calc_14 = true;
	set_type = MolMechModel::NOT_SET;
}
	
MMDihedral::~MMDihedral()
{

}

bool MMDihedral::operator == ( const MMDihedral& rhs) const
{
	if( pt1 != rhs.pt1 || pt2 != rhs.pt2 || pt3 != rhs.pt3 || pt4 != rhs.pt4)
		return false;
	return true;
}

bool MMDihedral::operator < ( const MMDihedral& rhs) const
{
	if( pt1 < rhs.pt1)
		return true;
	else if( pt1 > rhs.pt1)
		return true;

	if( pt2 < rhs.pt2)
		return true;
	else if( pt2 > rhs.pt2)
		return true;

	if( pt3 < rhs.pt3)
		return true;
	else if( pt3 > rhs.pt3)
		return true;

	return (pt4 < rhs.pt4);
}

AtomContact::AtomContact()
{
  pt1 = pt2 = NULL;
  cf.resize(2);

  cnt_type = AtomContactType::HARMONIC_CNT;
  set_type = MolMechModel::NOT_SET;
}

AtomContact::AtomContact(HaAtom* new_pt1, HaAtom* new_pt2, const AtomContactType& cnt_type_new)
{
	cnt_type = cnt_type_new;
	if(new_pt1 == NULL  ||  new_pt2 == NULL )
	{
		ErrorInMod("AtomContact::AtomContact()", 
		           " One of the HaAtom pointers is NULL ");
		pt1 = pt2 = NULL;
		return;
	}
	if(new_pt1 < new_pt2)
	{
		pt1 = new_pt1;
		pt2= new_pt2;
	}
	else
	{
		pt1 = new_pt2;
		pt2 = new_pt1;
	}

	cf.resize(2);
	cf[0] = 0.0;
	cf[1] = 0.0;

	set_type = MolMechModel::NOT_SET;
}


AtomContact::AtomContact(const AtomContact& cnt_ref)
{
	pt1      = cnt_ref.pt1;
	pt2      = cnt_ref.pt2;
	cf       = cnt_ref.cf;

	cnt_type = cnt_ref.cnt_type;
	set_type = cnt_ref.set_type;
}

AtomContact::~AtomContact()
{

}

bool AtomContact::operator==(const AtomContact& rhs) const
{
	if( this->cnt_type != rhs.cnt_type) return false;
	return( pt1 == rhs.pt1 && pt2 == rhs.pt2);
}

bool AtomContact::operator < (const AtomContact& rhs) const
{
	if( pt1 < rhs.pt1) 
		return true;
	if( pt1 > rhs.pt1)
		return false;
	return( pt2 < rhs.pt2);
}

int AtomContact::SetParamsEneR(double ene_min, double rmin)
{
	ene_min = fabs(ene_min);
	if( cnt_type == cnt_type.HARMONIC_CNT )
	{
		if( cf.size() != 2 ) cf.resize(2);
		cf[0] = rmin;
		cf[1] = ene_min;
	}
    else if(cnt_type == cnt_type.VDW_CNT_6_12 || cnt_type == cnt_type.VDW_CNT_6_12_NO_REP )
	{
		if(cf.size() != 2 ) cf.resize(2);	
		double r6 = rmin*rmin*rmin;
		r6 = r6*r6;
	
		cf[0] = ene_min*r6*r6;  // A - coef at R^-12
		cf[1] = 2.0*ene_min*r6; // B - coef at R^-6
	}
	else if( cnt_type == cnt_type.VDW_CNT_10_12 )
	{
		if(cf.size() != 2 ) cf.resize(2);
		double r2  = rmin*rmin;
		double r10 = r2 * r2 * r2 * r2 * r2;
	
		cf[0] = 5.0*ene_min*r10*r2; // A - coef at R^-12
		cf[1] = 6.0*ene_min*r10;    // B - coef at R^-10
		
	}
	return TRUE;
}

double AtomContact::GetRMin() const 
{ 
	double rmin = -1.0;
	if( cnt_type == cnt_type.HARMONIC_CNT )
	{
		if(cf.size() != 2 ) return rmin;
		return cf[0];
	}
	else if( cnt_type == cnt_type.VDW_CNT_6_12 || cnt_type == cnt_type.VDW_CNT_6_12_NO_REP )
	{
		if(cf.size() != 2 ) return rmin;
		if(cf[1] < 0.1e-10) return rmin;
		rmin = pow(2.0*cf[0]/cf[1],1.0/6.0);
		return rmin;
	}
	else if( cnt_type == cnt_type.VDW_CNT_10_12 )
	{
		if(cf.size() != 2 ) return rmin;
		if(cf[1] < 0.1e-10) return rmin;
		rmin = sqrt((6.0*cf[0])/(5.0*cf[1]));
		return rmin;
	}
	return rmin;
}

double AtomContact::GetEneMin() const 
{ 
	double emin = 1.0e10;
	if( cnt_type == cnt_type.HARMONIC_CNT )
	{
		if( cf.size() != 2 ) return cf[1];
	}
	else if( cnt_type == cnt_type.VDW_CNT_10_12 )
	{
		if(cf.size() != 2 ) return emin;
		if(cf[1] < 0.1e-10) return emin;
		double rmin6 = 2.0*cf[0]/cf[1];
		emin = cf[0]/(rmin6*rmin6) - cf[1]/rmin6;
		return emin;
	}
	else if( cnt_type == cnt_type.VDW_CNT_10_12 )
	{
		if(cf.size() != 2 ) return emin;
		if(cf[1] < 0.1e-10) return emin;
		double rmin2 = (6.0*cf[0])/(5.0*cf[1]);
		double rmin10 = rmin2*rmin2*rmin2*rmin2*rmin2;
		emin = cf[0]/(rmin10*rmin2) - cf[1]/rmin10;
		return emin;
	}
	else if( cnt_type == cnt_type.VDW_CNT_6_12_NO_REP )
	{
		return 0.0;
	}
	return emin;	
}

double AtomContact::GetHarmForceConst() const
{
	if( cnt_type != cnt_type.HARMONIC_CNT || cf.size() < 2 ) return -1.0;
	return cf[1];
}