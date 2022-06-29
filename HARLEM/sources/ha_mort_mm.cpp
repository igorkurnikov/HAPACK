/*! \file ha_mort_mm.cpp

   Classes and Function to define interface bewtween Molecular Mechanics parts MORT library and HARLEM.

  \author Igor Kurnikov 
  \date 2010-

*/

#include "mpi.h"

#include "boost/format.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "object.hpp"
#include "guilib.hpp"
#include "format.hpp"
#include "ambfmt.hpp"
#include "pdbent.hpp"

#include "harlemapp.h"
#include "halinalg.h"
#include "hamolset.h"
#include "hamolecule.h"
#include "haresdb.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "mm_force_field.h"
#include "mm_driver_amber.h"
#include "ha_mort_mm.h"

//static mort::database_t* CreateMortDB()
//{
//	
//}

MortForceField::MortForceField()
{
	p_mdb    = new mort::database_t();
	p_atomff = new mort::molecule_t(); 
	p_poleff = new mort::molecule_t();
}

MortForceField::~MortForceField()
{
	delete p_mdb;
	delete p_atomff;
	delete p_poleff;
}

int MMForceField::InitMortFF()
{
	PrintLog(" MMForceField::InitMortFF() pt 1 \n");
	if( p_mort_ff != NULL) 
	{
		delete p_mort_ff;
	}
	p_mort_ff = new MortForceField();

    mort::database_t* p_mdb = p_mort_ff->p_mdb;
	
	std::vector<std::string> res_files;
	std::vector<std::string> parm_files;

	if( ff_type == ForceFieldType::AMOEBA ) 
	{
		res_files.push_back("amoeba_amino.off");
		res_files.push_back("amoeba_aminont.off");
		res_files.push_back("amoeba_aminoct.off");
//		res_files.push_back("amoeba_watbox.off");
		res_files.push_back("amoeba_wat.off");
		res_files.push_back("amoeba_add_1.off");
//		res_files.push_back("amoeba_nmet.off");

		parm_files.push_back("amoebapro.prm");
//		parm_files.push_back("amoebapro09v4.prm");
	}
	else if( ff_type == ForceFieldType::AMBER_94 )
	{
//		res_files.push_back("all_nucleic94.lib");
		res_files.push_back("all_amino94.lib");
		res_files.push_back("all_aminoct94.lib");
		res_files.push_back("all_aminont94.lib");
//		res_files.push_back("ions94.lib");
//		res_files.push_back("solvents.lib");

		parm_files.push_back("parm94.dat");
		parm_files.push_back("gaff.dat");
	}
	else if( ff_type == ForceFieldType::AMBER_99_SB )
	{
		res_files.push_back("all_nucleic94.lib");
		res_files.push_back("all_amino94.lib");
		res_files.push_back("all_aminoct94.lib");
		res_files.push_back("all_aminont94.lib");
		res_files.push_back("ions94.lib");
		res_files.push_back("solvents.lib");

		parm_files.push_back("parm99.dat");
		parm_files.push_back("frcmod.ff99SB");
		parm_files.push_back("gaff.dat");
	}
	else if( ff_type == ForceFieldType::AMBER_99_BSC0 )
	{
		res_files.push_back("all_nucleic94.lib");
		res_files.push_back("all_amino94.lib");
		res_files.push_back("all_aminoct94.lib");
		res_files.push_back("all_aminont94.lib");
		res_files.push_back("DNA_CI.lib");
		res_files.push_back("ions94.lib");
		res_files.push_back("solvents.lib");

		parm_files.push_back("parm99.dat");
		parm_files.push_back("frcmod.ff99SB");
		parm_files.push_back("frcmod.parmbsc0");
		parm_files.push_back("gaff.dat");
	}
	else if( ff_type == ForceFieldType::AMBER_03 )
	{
		res_files.push_back("all_nucleic94.lib");
		res_files.push_back("all_amino03.lib");
		res_files.push_back("all_aminoct03.lib");
		res_files.push_back("all_aminont03.lib");
		res_files.push_back("DNA_CI.lib");
		res_files.push_back("ions94.lib");
		res_files.push_back("solvents.lib");

		parm_files.push_back("parm99.dat");
		parm_files.push_back("frcmod.ff03");
		parm_files.push_back("gaff.dat");
	}
	else if( ff_type == ForceFieldType::AMBER_10 )  //   ff10 =  ff99SB for proteins; ff99bsc0 for DNA; ff99sbsc_chiOL3 for RNA
	{
		res_files.push_back("amino10.lib");
		res_files.push_back("aminoct10.lib");
		res_files.push_back("aminont10.lib");
		res_files.push_back("nucleic10.lib");
		res_files.push_back("ions08.lib");
		res_files.push_back("solvents.lib");

		parm_files.push_back("parm10.dat");
		parm_files.push_back("frcmod.ionsjc_tip3p");
		parm_files.push_back("gaff.dat");
	}

	int i;
	try 
	{
		std::string amber_lib_dir;
		std::string amber_parm_dir;
		amber_lib_dir  = pApp->res_db_dir + "amber_lib"  + path_sep;
		amber_parm_dir = pApp->res_db_dir + "amber_parm" + path_sep;

		for( i = 0; i < res_files.size(); i++)
		{
			std::string res_fname = amber_lib_dir + res_files[i];
			std::ifstream stream( res_fname.c_str() );
			if( !stream.is_open() )
			{
				std::string msg = "Can not open residue file ";
				msg += res_files[i];
				msg += "\n";
				throw std::runtime_error(msg);
			}

			PrintLog(" Loading residue File %s \n", res_files[i].c_str());

			while( stream )
			{
				boost::shared_ptr< mort::molecule_t > pmol( new mort::molecule_t() );    
				mort::read_off( stream, *pmol );
				p_mdb->set( pmol->get_s(mort::NAME), pmol );
			} 
		}

		if( ff_type == ForceFieldType::AMOEBA )
		{
//			for( i = 0; i < parm_files.size(); i++)
			for( i = 0; i < tinker_param_files.size(); i++)
			{
//				std::string parm_path = amber_parm_dir + parm_files[i];
				std::string parm_path = tinker_param_files[i];
				PrintLog(" Loading Parameter File %s \n", parm_path.c_str());
				std::ifstream stream_parm( parm_path.c_str() );
				mort::read_amoeba_frc( stream_parm, (*p_mort_ff->p_atomff), (*p_mort_ff->p_poleff ) );
			}
		}
		else if ( ff_type == ForceFieldType::AMBER_94 || ff_type == ForceFieldType::AMBER_99_SB || ff_type == ForceFieldType::AMBER_99_BSC0
			      || ff_type == ForceFieldType::AMBER_03 || ff_type == ForceFieldType::AMBER_10 )
		{
			for( i = 0; i < parm_files.size(); i++)
			{
				std::string parm_path = amber_parm_dir + parm_files[i];
				PrintLog(" Loading Parameter File %s \n", parm_path.c_str());
				std::ifstream stream_parm( parm_path.c_str() );
				mort::read_frc( stream_parm, (*p_mort_ff->p_atomff) );
			}	
			p_mort_ff->p_atomff->set_s(mort::NAME,"_amberffp");
		}
	}
	catch (const std::runtime_error& ex )
	{
		PrintLog(" Error in MMForceField::InitMortFF() \n %s \n",ex.what());
		return FALSE;
	}
	PrintLog(" MMForceField::InitMortFF() pt end \n");
	return TRUE;
}

int MolMechModel::ClearMortModel()
{
	if( p_mort_model != NULL) delete p_mort_model;
	p_mort_model = new mort::molecule_t;
	return TRUE;
}

static int get_iparm( mort::morf_t& p, const mort::hashid_t& parmid )
{
	int ivalue;
	if( p.get_i(parmid, ivalue) )
	{
		return p.get_i(parmid);
	}
	else
	{
		std::string str = p.get_s(parmid);
		try { ivalue = boost::lexical_cast<int>(str);
		}
		catch(boost::bad_lexical_cast&)
		{
			PrintLog(" Error in get_iparm() converting %s \n",str.c_str());
			return 999999;
		}
	}
	return ivalue;
}

int AmberMMModel::InitAmberModelAmoeba()
{
	PrintLog(" AmberMMModel::InitAmberModelAmoeba() pt 1 \n");
	try
	{
		if( p_mm_model->ff_type != ForceFieldType::AMOEBA ) throw std::runtime_error(" Amoeba Force Field is not set for the molecular mechanics model ");
	
		mort::molecule_t& mol = *p_mm_model->p_mort_model;
		MMForceField* p_ff = MMForceField::GetMMForceField( ForceFieldType::AMOEBA, TRUE );

		if( p_ff == NULL )  throw std::runtime_error(" Can not load AMOEBA force field ");
		
		AtomIntMap atoms_idx = pmset->GetAtomSeqNumMap();

		//	int mol_sol_par = (*p_mm_model->p_mort_model).get_i(mort::SOLUTE);

		mort::molecule_t* aff = p_ff->p_mort_ff->p_atomff;
		mort::molecule_t* pff = p_ff->p_mort_ff->p_poleff; 

		mort::energee_t e( mol );

		mort::excl_t& excl = e.m_excl;   // Build excluded atom list ?
		exclude( mol, excl, 4 );
		mort::prmtop::amoeba::setup_adjust( mol, excl, *pff );

		e.assignparm( *aff, mort::AMOEBA );
		mort::parmset_t& params = e.m_parmset;
		const mort::nabparm_t& prm = e.getnabparm();

		natom  = prm.Natom;
		ntypes = prm.Ntypes;
		nttyp  = ntypes * (ntypes + 1) / 2;
		nbonh  = prm.Nbonh;
		ntheth = prm.Ntheth;
		nphih  = prm.Nphih;

		next   = prm.Nnb;
		nres   = prm.Nres;
		nbona  = prm.Nbona;
		ntheta = prm.Ntheta;
		nphia  = prm.Nphia;
		numbnd = prm.Numbnd;
		numang = prm.Numang;
		nptra  = prm.Nptra;
		nphb   = prm.Nphb;

		if( prm.IfBox == 1 )
		{
			pmset->per_bc->Set(true);
			pmset->per_bc->orthogonal_flag = true;
		}
		else if( prm.IfBox == 2 )
		{
			pmset->per_bc->Set(true);
			pmset->per_bc->orthogonal_flag = false;
			pmset->per_bc->octahedron_flag = true;
		}
		else if( prm.IfBox == 3 )
		{
			pmset->per_bc->Set(true);
			pmset->per_bc->orthogonal_flag = false;
		}
		else
		{
			pmset->per_bc->Set(false);
		}
		max_res_size = prm.Nmxrs;

		int ifcap_dummy  = prm.IfCap;
		int numextra_dummy = prm.Numextra;

		bonda_idx  = nbonh + 1;
		anglea_idx = ntheth + 1;
		diheda_idx = nphih + 1;  

		gbl_bond_allocsize  = nbonh + nbona;
		if( nbona == 0 ) gbl_bond_allocsize++; 
		gbl_angle_allocsize = ntheth + ntheta;
		if( ntheta == 0 ) gbl_angle_allocsize++; 
		gbl_dihed_allocsize = nphih + nphia;
		if( nphia == 0 ) gbl_dihed_allocsize++;

		int i;
		atm_igraph.resize(natom);
		for(i = 0; i < natom; i++)
		{
			const std::string str(prm.AtomNames);
			atm_igraph[i] = str.substr(4*i,4);
		}

		atm_charge.resize(natom);
		for(i = 0; i < natom; i++)
		{
			atm_charge[i] = prm.Charges[i]/p_mm_model->diel_const;

			if( p_amber_driver->emulate_ext_amber )
			{
				MMDriverAmber::ModifyFormatVal(atm_charge[i],FLOAT_E16_8);
			}
		}

		atm_mass.resize(natom);
		atm_mass_inv.resize(natom);
		tmass = 0.0;
		for(i = 0; i < natom; i++)
		{
			atm_mass[i] = prm.Masses[i];
			tmass += atm_mass[i];
			atm_mass_inv[i] = 1.0/atm_mass[i];
			if( p_amber_driver->emulate_ext_amber )
			{
				MMDriverAmber::ModifyFormatVal(atm_mass_inv[i],FLOAT_E16_8);
			}
		}

		atm_iac.resize(natom); // IAC: indexes of atoms in the array of atom types (Fortran based 1)
		for(i = 0; i < natom; i++)
		{
			atm_iac[i] = prm.Iac[i];
		}

		atm_numex.resize(natom);  // NUMBER_EXCLUDED_ATOMS
		for(i = 0; i < natom; i++)
		{
			atm_numex[i] = prm.Iblo[i];
		}

		typ_ico.resize(ntypes*ntypes);  // NONBONDED_PARM_INDEX: indexes of pairs of atom types in the arrays of nonbonded parameters
		for(i = 0; i < ntypes*ntypes; i++)
		{
			typ_ico[i] = prm.Cno[i];
		}

		res_labels.resize(nres);  // Residue Labels:
		for(i = 0; i < nres; i++)
		{
			const std::string str(prm.ResNames);
			res_labels[i] = str.substr(4*i,4);
		}

		gbl_res_atms.resize(nres+1); // RESIDUE_POINTER: indexes of first atoms of residues
		for(i = 0; i < nres; i++)
		{
			gbl_res_atms[i] = prm.Ipres[i];
		}
		gbl_res_atms[nres] = natom + 1;

		gbl_rk.resize(numbnd); // BOND_FORCE_CONSTANT:  RK : force constant for the bonds of each type, kcal/mol
		for(i = 0; i < numbnd; i++)
		{
			gbl_rk[i] = prm.Rk[i];
			if( p_amber_driver->emulate_ext_amber )
			{
				MMDriverAmber::ModifyFormatVal(gbl_rk[i],FLOAT_E16_8);
			}
		}

		gbl_req.resize(numbnd); // BOND_EQUIL_VALUE:  REQ : equilibrium bond length for the bonds of each type, angstroms
		for(i = 0; i < numbnd; i++)
		{
			gbl_req[i] = prm.Req[i];
			if( p_amber_driver->emulate_ext_amber )
			{
				MMDriverAmber::ModifyFormatVal(gbl_req[i],FLOAT_E16_8);
			}
		}

		gbl_tk.resize(numang);  // ANGLE_FORCE_CONSTANT: TK : force constant for the angles of each type, kcal/mol A**2
		for(i = 0; i < numang; i++)
		{
			gbl_tk[i] = prm.Tk[i];
			if( p_amber_driver->emulate_ext_amber )
			{
				MMDriverAmber::ModifyFormatVal(gbl_tk[i],FLOAT_E16_8);
			}
		}

		gbl_teq.resize(numang);  // ANGLE_EQUIL_VALUE: TEQ : the equilibrium angle for the angles of each type, degrees
		for(i = 0; i < numang; i++)
		{
			gbl_teq[i] = prm.Teq[i];
			if( p_amber_driver->emulate_ext_amber )
			{
				MMDriverAmber::ModifyFormatVal(gbl_teq[i],FLOAT_E16_8);
			}
		}

		gbl_pk.resize(nptra);  // DIHEDRAL_FORCE_CONSTANT: PK : force constant for the dihedrals of each type, kcal/mol
		for(i = 0; i < nptra; i++)
		{
			gbl_pk[i] = prm.Pk[i];
			if( p_amber_driver->emulate_ext_amber )
			{
				MMDriverAmber::ModifyFormatVal(gbl_pk[i],FLOAT_E16_8);
			}
		}

		gbl_pn.resize(nptra);  // DIHEDRAL_PERIODICITY: PN : periodicity of the dihedral of a given type
		for(i = 0; i < nptra; i++)
		{
			gbl_pn[i] = prm.Pn[i];
		}

		gbl_phase.resize(nptra);  // DIHEDRAL_PHASE: PHASE : phase of the dihedral of a given type
		for(i = 0; i < nptra; i++)
		{
			gbl_phase[i] = prm.Phase[i];
			if( p_amber_driver->emulate_ext_amber )
			{
				MMDriverAmber::ModifyFormatVal(gbl_phase[i],FLOAT_E16_8);
			}
		}

		CalcAddDihParams(); // Compute additional dihedral angles parameters for vectorized processing

		// scee = *prm.Scee;  // TEMPORAL FIX IGOR 
		// scnb = *prm.Scnb;  // TEMPORAL FIX IGOR

		gbl_cn1.resize(nttyp);   // LENNARD_JONES_ACOEF :  // CN1 : Lennard Jones r**12 terms for all possible atom type
		for(i = 0; i < nttyp; i++)		                   // interactions, indexed by ICO and IAC; for atom i and j
		{			                                       // where i < j, the index into this array is as follows
			gbl_cn1[i] = prm.Cn1[i];				       // (assuming the value of ICO(INDEX) is positive):
								                           // CN1(ICO(NTYPES*(IAC(i)-1)+IAC(j))).
			if( p_amber_driver->emulate_ext_amber )
			{
				MMDriverAmber::ModifyFormatVal(gbl_cn1[i],FLOAT_E16_8);
			}
		}
		gbl_cn2.resize(nttyp);   // LENNARD_JONES_BCOEF :  // CN2 : Lennard Jones r**6 terms for all possible atom type interactions
		for(i = 0; i < nttyp; i++)		                   
		{			                                       
			gbl_cn2[i] = prm.Cn2[i];	
			if( p_amber_driver->emulate_ext_amber )
			{
				MMDriverAmber::ModifyFormatVal(gbl_cn2[i],FLOAT_E16_8);
			}
		}

		gbl_bond.resize(3*gbl_bond_allocsize);
		for(i = 0; i < nbonh; i++) // BONDS_INC_HYDROGEN
		{
			gbl_bond[3*i]   = prm.BondHAt1[i]/3 + 1;
			gbl_bond[3*i+1] = prm.BondHAt2[i]/3 + 1;
			gbl_bond[3*i+2] = prm.BondHNum[i];
		}

		for(i = 0; i < nbona; i++) // BONDS_WITHOUT_HYDROGEN
		{
			gbl_bond[3*(nbonh + i)  ]  = prm.BondAt1[i]/3 + 1;
			gbl_bond[3*(nbonh + i)+1]  = prm.BondAt2[i]/3 + 1;
			gbl_bond[3*(nbonh + i)+2]  = prm.BondNum[i];
		}

		gbl_angle.resize(4*gbl_angle_allocsize);
		for(i = 0; i < ntheth; i++) // ANGLES_INC_HYDROGEN
		{
			gbl_angle[4*i  ] = prm.AngleHAt1[i]/3 + 1;
			gbl_angle[4*i+1] = prm.AngleHAt2[i]/3 + 1;
			gbl_angle[4*i+2] = prm.AngleHAt3[i]/3 + 1;
			gbl_angle[4*i+3] = prm.AngleHNum[i];
		}

		for(i = 0; i < ntheta; i++) // ANGLES_WITHOUT_HYDROGEN
		{
			gbl_angle[4*(ntheth + i)]   = prm.AngleAt1[i]/3 + 1;
			gbl_angle[4*(ntheth + i)+1] = prm.AngleAt2[i]/3 + 1;
			gbl_angle[4*(ntheth + i)+2] = prm.AngleAt3[i]/3 + 1;
			gbl_angle[4*(ntheth + i)+3] = prm.AngleNum[i];
		}

		gbl_dihed.resize(5*gbl_dihed_allocsize);
		for(i = 0; i < nphih; i++) // DIHEDRALS_INC_HYDROGEN
		{
			gbl_dihed[5*i]    = prm.DihHAt1[i]/3 + 1;
			gbl_dihed[5*i+1]  = prm.DihHAt2[i]/3 + 1;
			gbl_dihed[5*i+2]  = abs(prm.DihHAt3[i])/3 + 1;
			if( prm.DihHAt3[i] < 0 ) gbl_dihed[5*i+2] = -gbl_dihed[5*i+2];
			gbl_dihed[5*i+3]  = abs(prm.DihHAt4[i])/3 + 1; 
			if( prm.DihHAt4[i] < 0 ) gbl_dihed[5*i+3] = -gbl_dihed[5*i+3];
			gbl_dihed[5*i+4]  = prm.DihHNum[i];
		}

		for(i = 0; i < nphia; i++) // DIHEDRALS_WITHOUT_HYDROGEN
		{
			gbl_dihed[5*(nphih + i)]    = prm.DihAt1[i]/3 + 1;
			gbl_dihed[5*(nphih + i)+1]  = prm.DihAt2[i]/3 + 1;
			gbl_dihed[5*(nphih + i)+2]  = abs(prm.DihAt3[i])/3 + 1;
			if( prm.DihAt3[i] < 0 ) gbl_dihed[5*(nphih + i)+2] = -gbl_dihed[5*(nphih + i)+2];
			gbl_dihed[5*(nphih + i)+3]  = abs(prm.DihAt4[i])/3 + 1; 
			if( prm.DihAt4[i] < 0 ) gbl_dihed[5*(nphih + i)+3] = -gbl_dihed[5*(nphih + i)+3];
			gbl_dihed[5*(nphih + i)+4]  = prm.DihNum[i];
		}

		gbl_natex.resize(next);  // EXCLUDED_ATOMS_LIST
		// NATEX : the excluded atom list. To get the excluded list for atom 
		// "i" you need to traverse the NUMEX list, adding up all
		// the previous NUMEX values, since NUMEX(i) holds the number
		// of excluded atoms for atom "i", not the index into the 
		// NATEX list. Let IEXCL = SUM(NUMEX(j), j=1,i-1), then
		// excluded atoms are NATEX(IEXCL) to NATEX(IEXCL+NUMEX(i)).
		for(i = 0; i < next; i++)
		{
			gbl_natex[i] = prm.ExclAt[i];
		}

		atm_isymbl.resize(natom); //  AMBER ATOM FORCE FIELD SYMBOLS
		for(i = 0; i < natom; i++)
		{
			const std::string& str(prm.AtomSym);
			atm_isymbl[i] = str.substr(4*i,4);
		}

		atm_itree.resize(natom); // TREE_CHAIN_CLASSIFICATION
		for(i = 0; i < natom; i++)
		{
			const std::string& str(prm.AtomTree);
			atm_itree[i] = str.substr(4*i,4);
		}

		// write_iarray( os, "JOIN_ARRAY",         	prm.TreeJoin, prm.Natom );
		// write_iarray( os, "IROTAT",             	prm.AtomRes,  prm.Natom );

		if( prm.IfBox ) // Additional info for periodic boundary conditions and solvent info
		{
			// SOLVENT_POINTERS
			n_solute_res =  prm.Iptres;   
			nspm         =  prm.Nspm;           
			n_solute_mol =  prm.Nspsol - 1;

			atm_nsp.resize(nspm);
			for(i = 0; i < nspm; i++)
			{
				atm_nsp[i] = prm.Boundary[i];
			}

			pmset->per_bc->SetBox(prm.Box[0],prm.Box[1],prm.Box[2],90.0*DEG_TO_RAD,prm.Box[3]*DEG_TO_RAD);
		}

		//   write_string( os, "RADIUS_SET", GBPARM_NAME[2] );

		atm_gb_radii.resize(natom);
		for(i = 0; i < natom; i++)
		{
			atm_gb_radii[i] = prm.Rborn[i];
			if( p_amber_driver->emulate_ext_amber )
			{
				MMDriverAmber::ModifyFormatVal(atm_gb_radii[i],FLOAT_E16_8);
			}
		}

		atm_gb_fs.resize(natom);
		for(i = 0; i < natom; i++)
		{
			atm_gb_fs[i] = prm.Fs[i];
			if( p_amber_driver->emulate_ext_amber )
			{
				MMDriverAmber::ModifyFormatVal(atm_gb_fs[i],FLOAT_E16_8);
			}

		}

		// AMOEBA parameters

		atm_poltype.resize(natom);
		atm_element.resize(natom);
		atm_class_idx.resize(natom);

		mort::atomiter_t maitr = mol.atom_begin();
		i  = 0;
		for( ; maitr != mol.atom_end(); ++maitr ) 
		{
			atm_poltype[i]   = get_iparm((*maitr),mort::POLTYPE);
			atm_element[i]   = get_iparm((*maitr),mort::ELEMENT);
			atm_class_idx[i] = (*maitr).get_s(mort::TYPE);
			i++;
		}

		n_bond_amoeba = mol.nbond();
		gbl_bond_amoeba.resize( 3*n_bond_amoeba );
		i = 0;
		mort::bonditer_t mbitr = mol.bond_begin(); // Bond iterator of MORT molecule 
		for( ; mbitr != mol.bond_end(); mbitr++ )
		{
			gbl_bond_amoeba[3*i]    = mort::atom_1st(*mbitr).get_i(mort::ID);
			gbl_bond_amoeba[3*i+1]  = mort::atom_2nd(*mbitr).get_i(mort::ID);
			gbl_bond_amoeba[3*i+2]  = (*mbitr).get_i(mort::TYPEID);
			i++;
		}

		bond_amoeba_params.resize(2);
		bond_amoeba_params[0] = params.bond[0];
		bond_amoeba_params[1] = params.bond[1];
		n_bond_amoeba_params  = bond_amoeba_params[0].size();

		bond_amoeba_ftab_degree = 4;
		bond_amoeba_ftab_coef.resize(5);
		bond_amoeba_ftab_coef[0] = 0.0;
		bond_amoeba_ftab_coef[1] = 0.0;
		bond_amoeba_ftab_coef[2] = 1.0;
		bond_amoeba_ftab_coef[3] = aff->get_d("bond-cubic");
		bond_amoeba_ftab_coef[4] = aff->get_d("bond-quartic");

		n_urey_bond = 0;
		mort::angliter_t mang_itr = mol.angl_begin();
		for( ; mang_itr != mol.angl_end(); mang_itr++)
		{
			mort::numvec urey = (*mang_itr).get_v(mort::UREY);
			if( urey[0] == 0.0 ) continue;
			n_urey_bond++;
		}
		gbl_bond_urey.resize(3*n_urey_bond);

		i = 0;
		for( mang_itr = mol.angl_begin(); mang_itr != mol.angl_end(); mang_itr++)
		{
			mort::numvec urey = (*mang_itr).get_v(mort::UREY);
			if( urey[0] == 0.0 ) continue;

			gbl_bond_urey[3*i  ] = mort::atom_1st(*mang_itr).get_i(mort::ID);
			gbl_bond_urey[3*i+1] = mort::atom_3rd(*mang_itr).get_i(mort::ID);
			gbl_bond_urey[3*i+2] = (*mang_itr).get_i(mort::UREYID);
			i++;
		}

		bond_urey_params.resize(2);
		bond_urey_params[0] = params.urey[0];
		bond_urey_params[1] = params.urey[1];
		n_urey_bond_params = bond_urey_params[0].size();

		bond_urey_ftab_degree = 2;
		bond_urey_ftab_coef.resize(3);
		bond_urey_ftab_coef[0] = 0.0;
		bond_urey_ftab_coef[1] = 0.0;
		bond_urey_ftab_coef[2] = 1.0;

		n_angle_amoeba = 0;
		for( mang_itr = mol.angl_begin(); mang_itr != mol.angl_end(); mang_itr++)
		{
			if( (*mang_itr).noops() != 0 ) continue;
			n_angle_amoeba++;
		}

		gbl_angle_amoeba_reg.resize(4*n_angle_amoeba);
		i = 0;
		for( mang_itr = mol.angl_begin(); mang_itr != mol.angl_end(); mang_itr++)
		{
			if( (*mang_itr).noops() != 0 ) continue;

			gbl_angle_amoeba_reg[4*i  ] = mort::atom_1st(*mang_itr).get_i(mort::ID);
			gbl_angle_amoeba_reg[4*i+1] = mort::atom_2nd(*mang_itr).get_i(mort::ID);
			gbl_angle_amoeba_reg[4*i+2] = mort::atom_3rd(*mang_itr).get_i(mort::ID);
			gbl_angle_amoeba_reg[4*i+3] = (*mang_itr).get_i(mort::TYPEID);
			i++;
		}

		angle_amoeba_params.resize(2);
		angle_amoeba_params[0] = params.angl[0];
		angle_amoeba_params[1] = params.angl[1];
		n_angle_amoeba_params = angle_amoeba_params[0].size();

		angle_amoeba_ftab_degree = 6;
		angle_amoeba_ftab_coef.resize(7);
		angle_amoeba_ftab_coef[0] = 0.0;
		angle_amoeba_ftab_coef[1] = 0.0;
		angle_amoeba_ftab_coef[2] = 1.0;
		angle_amoeba_ftab_coef[3] = aff->get_d("angle-cubic");
		angle_amoeba_ftab_coef[4] = aff->get_d("angle-quartic");
		angle_amoeba_ftab_coef[5] = aff->get_d("angle-pentic");
		angle_amoeba_ftab_coef[6] = aff->get_d("angle-sextic");

		n_trig_angles = mol.noops();
		gbl_angle_amoeba_trig.resize(5*n_trig_angles);

		mort::impriter_t oops = mol.impr_begin();
		i = 0;

		for( ; oops != mol.impr_end(); oops++)
		{
			gbl_angle_amoeba_trig[5*i  ]  = mort::atom_1st(*oops).get_i(mort::ID);
			gbl_angle_amoeba_trig[5*i+1]  = mort::atom_3rd(*oops).get_i(mort::ID);
			gbl_angle_amoeba_trig[5*i+2]  = mort::atom_2nd(*oops).get_i(mort::ID);
			gbl_angle_amoeba_trig[5*i+3]  = mort::atom_4th(*oops).get_i(mort::ID);
			gbl_angle_amoeba_trig[5*i+4]  = mort::angl_t::get( mort::atom_1st(*oops), mort::atom_3rd(*oops), mort::atom_2nd(*oops) ).get_i(mort::TYPEID);
			i++;
		}

		n_opbend_angles = mol.noops();
		gbl_opbend_angle.resize(5*n_opbend_angles);

		i = 0;

		for( oops = mol.impr_begin(); oops != mol.impr_end(); oops++)
		{
			gbl_opbend_angle[5*i  ] = mort::atom_1st(*oops).get_i(mort::ID);
			gbl_opbend_angle[5*i+1] = mort::atom_3rd(*oops).get_i(mort::ID);
			gbl_opbend_angle[5*i+2] = mort::atom_2nd(*oops).get_i(mort::ID);
			gbl_opbend_angle[5*i+3] = mort::atom_4th(*oops).get_i(mort::ID);
			gbl_opbend_angle[5*i+4] = (*oops).get_i(mort::TYPEID);
			i++;
		}

		opbend_angle_params = params.oops[0];
		n_opbend_angles_params = opbend_angle_params.size();

		n_tors_amoeba = 0;
		mort::diheiter_t tors_itr = mol.dihe_begin();
		for( ; tors_itr != mol.dihe_end(); tors_itr++)
		{
			double fc = (*tors_itr).get_d(mort::FORCE);
			if( fc < 1e-6 && fc > -1e-6 ) continue;
			n_tors_amoeba++;
		}

		gbl_amoeba_tors_angle.resize(5*n_tors_amoeba);
		i = 0;
		for( tors_itr = mol.dihe_begin(); tors_itr != mol.dihe_end(); tors_itr++)
		{
			double fc = (*tors_itr).get_d(mort::FORCE);
			if( fc < 1e-6 && fc > -1e-6 ) continue; 

			gbl_amoeba_tors_angle[5*i  ] = mort::atom_1st(*tors_itr).get_i(mort::ID);
			gbl_amoeba_tors_angle[5*i+1] = mort::atom_2nd(*tors_itr).get_i(mort::ID);
			gbl_amoeba_tors_angle[5*i+2] = mort::atom_3rd(*tors_itr).get_i(mort::ID);
			gbl_amoeba_tors_angle[5*i+3] = mort::atom_4th(*tors_itr).get_i(mort::ID);
			gbl_amoeba_tors_angle[5*i+4] = (*tors_itr).get_i(mort::TYPEID);
			i++;
		}

		tors_amoeba_params.resize(3);
		tors_amoeba_params[0] = params.tors[0];
		tors_amoeba_params[1] = params.tors[1];
		tors_amoeba_params[2] = params.tors[2];
		n_tors_amoeba_params = tors_amoeba_params[0].size();

		n_pi_torsions = mol.nptor();
		gbl_pi_tors_angle.resize(7*n_pi_torsions);

		mort::mobjiter_t pi_tors_itr = mol.ptor_begin();
		i = 0;
		for( ; pi_tors_itr != mol.ptor_end(); pi_tors_itr++)
		{
			gbl_pi_tors_angle[7*i  ] = mort::atom_1st(*pi_tors_itr).get_i(mort::ID);
			gbl_pi_tors_angle[7*i+1] = mort::atom_2nd(*pi_tors_itr).get_i(mort::ID);
			gbl_pi_tors_angle[7*i+2] = mort::atom_3rd(*pi_tors_itr).get_i(mort::ID);
			gbl_pi_tors_angle[7*i+3] = mort::atom_4th(*pi_tors_itr).get_i(mort::ID);
			mort::morf_t atm5 = *((*pi_tors_itr).atom_begin()+4);
			gbl_pi_tors_angle[7*i+4] = atm5.get_i(mort::ID);
			mort::morf_t atm6 = *((*pi_tors_itr).atom_begin()+5);
			gbl_pi_tors_angle[7*i+5] = atm6.get_i(mort::ID);
			gbl_pi_tors_angle[7*i+6] = (*pi_tors_itr).get_i(mort::TYPEID);
			i++;
		}

		pi_tors_params.resize(3);
		pi_tors_params[0] = params.ptor[0];
		pi_tors_params[1] = params.ptor[1];
		pi_tors_params[2] = params.ptor[2];
		n_pi_torsions_params = pi_tors_params[0].size();

		n_stretch_bend = 0;
		for( mang_itr = mol.angl_begin(); mang_itr != mol.angl_end(); mang_itr++)
		{
			mort::numvec strbnd = (*mang_itr).get_v(mort::STRBND);
			if( strbnd.size() == 0 ) continue;
			if( fabs(strbnd[0]) < 1e-8  ) continue;
			n_stretch_bend++;
		}
		gbl_str_bend_angle.resize(4*n_stretch_bend);

		i = 0;
		for( mang_itr = mol.angl_begin(); mang_itr != mol.angl_end(); mang_itr++)
		{
			mort::numvec strbnd = (*mang_itr).get_v(mort::STRBND);
			if( strbnd.size() == 0 ) continue;
			if( fabs(strbnd[0]) < 1e-8  ) continue;

			gbl_str_bend_angle[4*i  ] = mort::atom_1st(*mang_itr).get_i(mort::ID);
			gbl_str_bend_angle[4*i+1] = mort::atom_2nd(*mang_itr).get_i(mort::ID);
			gbl_str_bend_angle[4*i+2] = mort::atom_3rd(*mang_itr).get_i(mort::ID);
			gbl_str_bend_angle[4*i+3] = (*mang_itr).get_i((int)mort::STRBNDID);     // IGOR TEMPORAL FIX  compensating error in MORT library 
			i++;
		}

		str_bend_params.resize(4);
		str_bend_params[0] = params.strbnd[0];
		str_bend_params[1] = params.strbnd[1];
		str_bend_params[2] = params.strbnd[2];
		str_bend_params[3] = params.strbnd[3];
		n_stretch_bend_params = str_bend_params[0].size();

		n_tors_tors = mol.ntor2();
		gbl_tors_tors.resize(6*n_tors_tors);

		mort::mobjiter_t tors_tors_itr = mol.tor2_begin();
		i = 0;
		for( ; tors_tors_itr != mol.tor2_end(); tors_tors_itr++)
		{
			gbl_tors_tors[6*i  ] = mort::atom_1st(*tors_tors_itr).get_i(mort::ID);
			gbl_tors_tors[6*i+1] = mort::atom_2nd(*tors_tors_itr).get_i(mort::ID);
			gbl_tors_tors[6*i+2] = mort::atom_3rd(*tors_tors_itr).get_i(mort::ID);
			gbl_tors_tors[6*i+3] = mort::atom_4th(*tors_tors_itr).get_i(mort::ID);
			mort::morf_t atm5 = *((*tors_tors_itr).atom_begin()+4);
			gbl_tors_tors[6*i+4] = atm5.get_i(mort::ID);
			gbl_tors_tors[6*i+5] = (*tors_tors_itr).get_i(mort::TYPEID);
			i++;
		}

		n_tors_tors_params = aff->ntor2();
		tors_tors_id_params.resize(n_tors_tors_params);
		tors_tors_params.resize(n_tors_tors_params);

		mort::mobjiter_t parm = aff->tor2_begin();
		i = 0;
		for( ; parm != aff->tor2_end(); ++parm )
		{
			tors_tors_params[i].resize(6);
			tors_tors_params[i][0].resize(25);
			tors_tors_params[i][1].resize(25);
			tors_tors_params[i][2].resize(25*25);
			tors_tors_params[i][3].resize(25*25);
			tors_tors_params[i][4].resize(25*25);
			tors_tors_params[i][5].resize(25*25);

			tors_tors_id_params[i] = parm->absid() + 1;

			double start = -180.0;
			for( int j=0; j < 25; j++ )
			{
				tors_tors_params[i][0][j] = start + j*15.0;
				tors_tors_params[i][1][j] = start + j*15.0;
			}

			typedef boost::shared_ptr< std::vector<double> > dvec_ptr;
			dvec_ptr func = boost::any_cast< dvec_ptr >( parm->get_a(mort::TOR2FUNC) );

			for( int j=0; j < 625; ++j )
			{
				tors_tors_params[i][2][j] = func->at(6*j+2);
				tors_tors_params[i][3][j] = func->at(6*j+3);
				tors_tors_params[i][4][j] = func->at(6*j+4);
				tors_tors_params[i][5][j] = func->at(6*j+5);
			}
			i++;
		}

		atm_amoeba_vdw_type.resize(natom);
		atm_parent_id.resize(natom);
		atm_parent_weight.resize(natom);

		maitr = mol.atom_begin();
		i  = 0;
		for( maitr = mol.atom_begin(); maitr != mol.atom_end(); ++maitr ) 
		{
			atm_amoeba_vdw_type[i]   = (*maitr).get_i(mort::TYPEID);
			atm_parent_weight[i]     = (*maitr).get_d(mort::VBUFF);
			if( maitr->get_i(mort::ELEMENT) == mort::HYDROGEN )
			{
				atm_parent_id[i] = maitr->atom_begin()->get_i(mort::ID);
			}
			else
			{
				atm_parent_id[i] = maitr->get_i(mort::ID);
			}
			i++;
		}

		vdw_buffer_delta = 0.07; 
		vdw_buffer_gamma = 0.12;

		n_vdw_params =  params.vdw[0].size();
		amoeba_vdw_rstars.resize(n_vdw_params*n_vdw_params);
		amoeba_vdw_depths.resize(n_vdw_params*n_vdw_params);

		int ij = 0;
		for( i=0; i < n_vdw_params; ++i )
		{
			for( int j=0; j < n_vdw_params; ++j )
			{
				double ri = params.vdw[1][i];
				double rj = params.vdw[1][j];
				amoeba_vdw_rstars[ij] = (ri*ri*ri + rj*rj*rj)/(ri*ri + rj*rj);

				double ei = params.vdw[0][i];
				double ej = params.vdw[0][j];
				amoeba_vdw_depths[ij] = (4*ei*ej)/pow( sqrt(ei)+sqrt(ej), 2 );

				if( p_amber_driver->emulate_ext_amber )
				{
					MMDriverAmber::ModifyFormatVal(amoeba_vdw_rstars[ij],FLOAT_E16_8);
					MMDriverAmber::ModifyFormatVal(amoeba_vdw_depths[ij],FLOAT_E16_8);
				}

				ij++;
			}
		}

		num_adjust_list = excl.full_size();
		atm_adjust_list.resize( num_adjust_list*3);

		ij = 0;
		for( i=0; i < (int)excl.list.size(); ++i )
		{
			for( int j=0; j < (int)excl.list[i].size(); ++j )
			{
				atm_adjust_list[3*ij  ] = i+1;
				atm_adjust_list[3*ij+1] = excl.list[i][j];
				atm_adjust_list[3*ij+2] = excl.dist[i][j];
				ij++;
			}
		}

		adjust_vdw_weights.resize(9);

		adjust_vdw_weights[0] = 0.0; adjust_vdw_weights[1] = 0.0; adjust_vdw_weights[2] = 1.0; 
		adjust_vdw_weights[3] = 1.0; adjust_vdw_weights[4] = 0.0; adjust_vdw_weights[5] = 0.0; 
		adjust_vdw_weights[6] = 1.0; adjust_vdw_weights[7] = 1.0; adjust_vdw_weights[8] = 1.0; 

		adjust_mpole_weights.resize(9);

		adjust_mpole_weights[0] = 0.0; adjust_mpole_weights[1] = 0.0; adjust_mpole_weights[2] = 0.4; 
		adjust_mpole_weights[3] = 0.8; adjust_mpole_weights[4] = 0.0; adjust_mpole_weights[5] = 0.0; 
		adjust_mpole_weights[6] = 0.4; adjust_mpole_weights[7] = 0.8; adjust_mpole_weights[8] = 1.0; 

		adjust_direct_weights.resize(9);

		adjust_direct_weights[0] = 1.0; adjust_direct_weights[1] = 1.0; adjust_direct_weights[2] = 1.0; 
		adjust_direct_weights[3] = 1.0; adjust_direct_weights[4] = 0.0; adjust_direct_weights[5] = 0.0; 
		adjust_direct_weights[6] = 0.0; adjust_direct_weights[7] = 0.0; adjust_direct_weights[8] = 0.0; 

		adjust_polar_weights.resize(9);

		adjust_polar_weights[0] = 0.0; adjust_polar_weights[1] = 0.0; adjust_polar_weights[2] = 1.0; 
		adjust_polar_weights[3] = 1.0; adjust_polar_weights[4] = 0.0; adjust_polar_weights[5] = 0.0; 
		adjust_polar_weights[6] = 0.5; adjust_polar_weights[7] = 1.0; adjust_polar_weights[8] = 1.0; 

		adjust_mutual_weights.resize(9);

		adjust_mutual_weights[0] = 1.0; adjust_mutual_weights[1] = 1.0; adjust_mutual_weights[2] = 1.0; 
		adjust_mutual_weights[3] = 1.0; adjust_mutual_weights[4] = 1.0; adjust_mutual_weights[5] = 1.0; 
		adjust_mutual_weights[6] = 1.0; adjust_mutual_weights[7] = 1.0; adjust_mutual_weights[8] = 1.0; 

		num_local_multipoles = natom;
		atm_multipoles.resize(10*num_local_multipoles);  
		atm_multipoles = 0.0;
		atm_polar.resize(natom);
		atm_hpolar.resize(natom);
		atm_screen_polar.resize(natom);
		atm_screen_polar = sqrt( 0.39 );
		damp_polar_strength.resize(natom);
		damp_polar_strength = 0.0;
		damp_polar_sensitivity.resize(natom);
		damp_polar_sensitivity = 0.0;
		damp_polar_rad.resize(natom);
		damp_polar_rad = 1.0;
	
		if( p_mm_model->setup_params_from_mort_flag )
		{
			vector< vector<int> > chirials;
			vector< vector<int> > regulars;

			maitr = mol.atom_begin();
			for( ; maitr != mol.atom_end(); ++maitr )
			{
				mort::prmtop::amoeba::create_atomic_frame( *maitr, *pff, chirials, regulars );
			}
			
			num_chiral_frames = chirials.size();
			atm_chiral_frames.resize( 3*num_chiral_frames );
			for( i=0; i < num_chiral_frames; ++i )
			{
				atm_chiral_frames[3*i  ] = chirials[i][0];
				atm_chiral_frames[3*i+1] = chirials[i][1];
				atm_chiral_frames[3*i+2] = chirials[i][2];
			}

			num_reg_frames = regulars.size();
			atm_reg_frames.resize( 5*num_reg_frames );
			for( i=0; i < num_reg_frames; ++i )
			{
				atm_reg_frames[5*i  ] = regulars[i][0];
				atm_reg_frames[5*i+1] = regulars[i][1];
				atm_reg_frames[5*i+2] = regulars[i][2];
				atm_reg_frames[5*i+3] = regulars[i][3];
				atm_reg_frames[5*i+4] = regulars[i][4];
			}

			maitr = mol.atom_begin();
			i = 0;
			for( ; maitr != mol.atom_end(); ++maitr )
			{
				atm_multipoles[10*i  ] = maitr->get_d(mort::PCHG); 
				if( maitr->has_v(mort::POLE) )
				{
					mort::numvec pole = maitr->get_v(mort::POLE);

					atm_multipoles[10*i+1] = pole[0]; 
					atm_multipoles[10*i+2] = pole[1]; 
					atm_multipoles[10*i+3] = pole[2]; 
					atm_multipoles[10*i+4] = pole[3]*0.5; 
					atm_multipoles[10*i+5] = pole[5]*0.5; 
					atm_multipoles[10*i+6] = pole[8]*0.5; 
					atm_multipoles[10*i+7] = pole[4]; 
					atm_multipoles[10*i+8] = pole[6]; 
					atm_multipoles[10*i+9] = pole[7];
				}
				i++;
			}

			i  = 0;
			for( maitr = mol.atom_begin(); maitr != mol.atom_end(); ++maitr ) 
			{
				std::string polt = (*maitr).get_s(mort::POLTYPE);
			 	mort::morf_t pole( *pff, mort::ATOM, atoi(polt.c_str())-1 );
				atm_polar[i]   = pole.get_d(mort::POLAR);
				
	//			PrintLog(" InitAmberModelAmoeba() pt 1  i = %4d  atm_polar[i] = %s \n", i, polt.c_str());
				if( atoi(polt.c_str()) == 202 ) // water Oxygen - polarizability is affected by presence of Magnesium and Calcium 
				{	
					damp_polar_strength[i] = 0.0;
		//			damp_polar_sensitivity[i] = 0.5;
					damp_polar_sensitivity[i] = 0.0;
					damp_polar_sensitivity[i] = 0.0;
					damp_polar_rad[i] = 0.0;
					atm_screen_polar[i] = sqrt(0.39);
		//			atm_screen_polar[i] = 0.095/sqrt(0.39);
		//			atm_screen_polar[i] = 0.05/sqrt(0.39);
				}
				else if( atoi(polt.c_str()) == 210 ) // Magnesium is affecting polarization of water oxygen
				{	
//					damp_polar_strength[i] = 1.0;
					damp_polar_strength[i] = 0.0;
					damp_polar_sensitivity[i] = 0.0;
//					damp_polar_rad[i] = 1.0;
					damp_polar_rad[i] = 0.0;
					atm_screen_polar[i] = 0.095/sqrt(0.39);
		//			atm_screen_polar[i] = 0.05/sqrt(0.39);
				}
				else if( atoi(polt.c_str()) == 211 ) // Calcium 
				{	
					damp_polar_strength[i] = 1.0;
					damp_polar_sensitivity[i] = 0.1;
					damp_polar_rad[i] = 1.5;
					atm_screen_polar[i] = 0.159/sqrt(0.39);
				}
				else
				{
					damp_polar_strength[i] = 0.0;
					damp_polar_sensitivity[i] = 0.0;
					damp_polar_rad[i] = 0.0;
					atm_screen_polar[i] = sqrt(0.39);
				}
				i++;
			} 
		}
		else
		{
			i = 0;
			AtomIntMap& atoms_idx = p_mm_model->GetAtIdxMap(TRUE);
			HaResDB* p_res_db = HaResDB::GetDefaultResDB();
			ResidueIteratorMolSet ritr(this->pmset);
			HaResidue* pres;
			for( pres = ritr.GetFirstRes(); pres; pres = ritr.GetNextRes())
			{
				std::string res_fname = pres->GetFullName();
				ResFFTemplate* p_res_ff_templ = p_ff->GetResidueTemplate(res_fname);
				if( p_res_ff_templ == NULL ) throw std::runtime_error(" No FF template for residue " + res_fname );

				StrAtomMap templ_atname_to_res_map;
				int ires = p_res_db->GetTemplResAtNameAtomMap( pres, templ_atname_to_res_map ); 
				if(!ires) throw std::runtime_error(" Can not set residue template atom name -atom maps for residue " + res_fname );
				
				AtomIteratorAtomGroup aitr_res(pres);
				HaAtom* aptr;
				int j;
				for( aptr = aitr_res.GetFirstAtom(); aptr; aptr = aitr_res.GetNextAtom() )
				{
					std::string at_name = aptr->GetName();
					AtomFFParam* p_at_ff = p_res_ff_templ->GetAtomFFParam(at_name);
					if( p_at_ff == NULL ) throw std::runtime_error(" No FF template for atom name " + at_name + "In residue " + res_fname );
					
					AddAtomFrames(aptr, p_at_ff, &templ_atname_to_res_map );

					for(j =0; j < 10; j++)
					{
						atm_multipoles[10*i + j ] = 0.0;
					}

					atm_multipoles[10*i  ] = p_at_ff->GetCharge();
					if( p_at_ff->HasDipole() )
					{
						HaVec_double dipole = p_at_ff->GetDipole();
						atm_multipoles[10*i+1] = dipole[0]; 
						atm_multipoles[10*i+2] = dipole[1]; 
						atm_multipoles[10*i+3] = dipole[2];
					}
					
					if( p_at_ff->HasQPole() )
					{
						HaVec_double qpole = p_at_ff->GetQPole();
						atm_multipoles[10*i+4] = qpole[0]*0.5; 
						atm_multipoles[10*i+5] = qpole[2]*0.5; 
						atm_multipoles[10*i+6] = qpole[5]*0.5;
						atm_multipoles[10*i+7] = qpole[1]; 
						atm_multipoles[10*i+8] = qpole[3]; 
						atm_multipoles[10*i+9] = qpole[4];
					}

					if( p_at_ff->HasPolar() )
					{
						atm_polar[i] = p_at_ff->polar;
					}
					else
					{
						atm_polar[i] = 0.0;
					}

					if( p_at_ff->HasHPolar() )
					{
						atm_hpolar[i] = p_at_ff->hpolar;
					}
					else
					{
						atm_hpolar[i] = 0.0;
					}

					if( p_at_ff->HasScreenPolar() )
					{
						atm_screen_polar[i] = p_at_ff->screen_polar/sqrt(0.39);
					}
					else
					{
						atm_screen_polar[i] = sqrt(0.39);
					}

					damp_polar_strength[i]    = p_at_ff->damp_polar_strength;
					damp_polar_sensitivity[i] = p_at_ff->damp_polar_sensitivity;
					damp_polar_rad[i] = p_at_ff->damp_polar_rad;

		//			PrintLog(" InitAmberModelAmoeba() pt 2  i = %4d  elem_no = %d \n", i,aptr->GetElemNo());
					
					//if( aptr->GetElemNo() == 8 ) // water Oxygen - polarizability is affected by the presence of Magnesium and Calcium 
					//{	
					//	damp_polar_strength[i] = 0.0;
					//	damp_polar_sensitivity[i] = 0.1;
					//	damp_polar_rad[i] = 0.0;
					//}
					//else if( aptr->GetElemNo() == 12 ) // Magnesium
					//{
					//	damp_polar_strength[i] = 1.0; 
					//	damp_polar_sensitivity[i] = 0.0; 
					//	damp_polar_rad[i] = 1.0; 
					//}
					//else if( aptr->GetElemNo() == 20 ) // Calcium
					//{
					//	damp_polar_strength[i] = 1.0; 
					//	damp_polar_sensitivity[i] = 0.0; 
					//	damp_polar_rad[i] = 1.5; 
					//}
					//else
					//{
					//	damp_polar_strength[i] = 0.0; 
					//	damp_polar_sensitivity[i] = 0.0; 
					//	damp_polar_rad[i] = 0.0; 
					//}
					i++;
				}
			}
		}

		num_reg_frames = atm_reg_frames.size()/5;
		num_chiral_frames = atm_chiral_frames.size()/3;

//		PrintLog(" AmberMMModel::InitAmberModelAmoeba():  Chiral Frames: \n"); 
//		for( i = 0; i < num_chiral_frames; i++)
//		{
//			PrintLog(" %d    %5d %5d %5d \n",i,atm_chiral_frames[3*i],atm_chiral_frames[3*i+1],atm_chiral_frames[3*i+2]);
//		}
//		PrintLog(" AmberMMModel::InitAmberModelAmoeba():  Regular Frames: \n"); 
//		for( i = 0; i < num_reg_frames; i++)
//		{ 
//			PrintLog(" %d    %5d %5d %5d %5d %5d \n",i,atm_reg_frames[5*i],atm_reg_frames[5*i+1],atm_reg_frames[5*i+2],
//				                                       atm_reg_frames[5*i+3],atm_reg_frames[5*i+4]);
//		}
//		PrintLog(" AmberMMModel::InitAmberModelAmoeba():  Atomic multipoles: \n"); 
//		for( i = 0; i < 3; i++)
//		{
//			PrintLog(" %d    %12.6f \n ",i,atm_multipoles[10*i]);
//			PrintLog("       %12.6f %12.6f %12.6f \n ",atm_multipoles[10*i+1],atm_multipoles[10*i+2],atm_multipoles[10*i+3]);
//			PrintLog("       %12.6f \n  ", atm_multipoles[10*i+4]*2.0);
//			PrintLog("       %12.6f %12.6f \n  ", atm_multipoles[10*i+7],atm_multipoles[10*i+5]*2.0);
//			PrintLog("       %12.6f %12.6f %12.6f \n \n", atm_multipoles[10*i+8],atm_multipoles[10*i+9],atm_multipoles[10*i+6]*2.0);
//		}

		is_polarizable.resize(natom);
		atm_qterm.resize(natom);
		double polmin = 1.0e-12; 
		for( i = 0; i < natom; i++)
		{
			if( atm_polar[i] > polmin )
			{
				is_polarizable[i] = TRUE;
				atm_qterm[i] = 1.0 / sqrt( atm_polar[i] );
//				atm_qterm[i] = 1.0;   // IGOR TMP universal atm_qterm to check variable polarization code
				if( p_amber_driver->emulate_ext_amber )
				{
					MMDriverAmber::ModifyFormatVal(atm_qterm[i],FLOAT_E16_8);
				}
			}
			else
			{
				is_polarizable[i] = FALSE;
				atm_qterm[i] = 1.0;
			}	
		}
	}
	catch( std::exception& ex)
	{
		PrintLog("Error in AmmberMMModel::InitAmberModelAmoeba() \n");
		PrintLog("%s\n",ex.what());
		return FALSE;
	}

	// morf_t pole = find_pole( poleff, *atom );
     //               format( pole.get_d(POLAR) );

					//string polt = atom.get_s(POLTYPE);
     //           return morf_t( poleff, ATOM, atoi(polt.c_str())-1 );


	iamoeba = 1;
	PrintLog(" AmberMMModel::InitAmberModelAmoeba() pt end \n");
	return TRUE;
}

int MMForceField::SaveResFFTemplatesFromMort(const char* fname, MolSet* pmset_res )
{
	ofstream os(fname);
	if(!os.good()) return FALSE;
	
	char buf[256];

//	std::map<std::string, HaResidue*, less<std::string> > str_templ_map;
	HaResDB* p_res_db = HaResDB::GetDefaultResDB();

//	int nm = p_res_db->HostMolecules.size();
	int i;
//	for(i = 0; i < nm; i++)
//	{ 
//		HaMolecule* p_mol_templ = p_res_db->HostMolecules[i];
//		std::string templ_name = p_mol_templ->GetObjName();
//		HaResidue*  p_res_templ = p_res_db->GetTemplateForResidue(templ_name.c_str());
//		std::string amber_res_name = MMForceField::GetAmberResName(templ_name);
//		str_templ_map[amber_res_name] = p_res_templ;
//	}

	os << harlem::StdXMLHeader()     << std::endl;
	os << harlem::HarlemDataHeader() << std::endl;

	mort::molecule_t mort_res_mol; // MORT Molecule corresponding res_mol
	try
	{
		pmset_res->SetMortMol( mort_res_mol, this->ff_type );
		mort::parmset_t params;
		mort::parametrize( mort_res_mol, *(this->p_mort_ff->p_atomff), params );
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in MMForceField::SaveResFFTemplatesFromMort() \n");
		PrintLog(" Error to parametrize molecular set using MORT \n");
		PrintLog(" %s \n",ex.what() );
		return FALSE;
	}

	ResidueIteratorMolSet ritr_hlm(pmset_res);
	HaResidue* pres;

//	  mort::database_t::iterator mitr;
 //   mitr = p_mort_ff->p_mdb->begin();
 //   for( ; mitr != p_mort_ff->p_mdb->end(); ++mitr )
//  mort::resditer_t ritr_mort = mitr->resd_begin();
	mort::resditer_t ritr_mort = mort_res_mol.resd_begin();
	mort::resditer_t ritr_mort_prev = mort_res_mol.resd_end();
	
	for( pres = ritr_hlm.GetFirstRes(); pres; pres = ritr_hlm.GetNextRes(), ++ritr_mort )
    {   
		mort::resditer_t ritr_mort_next = ritr_mort;
		++ritr_mort_next;

		std::string res_full_name = pres->GetFullName();
		std::string amber_res_name = MMForceField::GetAmberResName(res_full_name);

		if( !p_mort_ff->p_mdb->has( amber_res_name ) )
		{
			PrintLog(" Can't find MORT residue template %s  corresponding to HARLEM residue %s \n",
				       amber_res_name.c_str(),res_full_name.c_str());
			continue;
		}
//		std::string amber_res_name = mitr->first;
//		boost::shared_ptr< mort::molecule_t > pmol = boost::dynamic_pointer_cast< mort::molecule_t >( mitr->second );
		boost::shared_ptr< mort::molecule_t > pmol = boost::dynamic_pointer_cast< mort::molecule_t >( p_mort_ff->p_mdb->get_mol(amber_res_name)); 
		
//		HaResidue*  p_res_templ = NULL;
//		PrintLog(" Residue Name = %s\n",amber_res_name.c_str());
//		if( str_templ_map.count(amber_res_name) > 0 )
//		{
//			p_res_templ = str_templ_map[amber_res_name];
//		}
		HaResidue*  p_res_templ = p_res_db->GetTemplateForResidue( res_full_name.c_str() );
		if( p_res_templ == NULL )
		{
			PrintLog(" Can not find HARLEM residue template %s for AMBER residue template %s \n", res_full_name.c_str(), amber_res_name.c_str()); 
			continue;
		}
		
		vector< vector<int> > chirials;
		vector< vector<int> > regulars;

		if( !pmol ) continue;
		try 
		{
			if( ff_type == ForceFieldType::AMOEBA )
			{
				switch_amoeba_type( *pmol,*(p_mort_ff->p_poleff) );
			}
			//			mort::parametrize( *pmol, *(p_mort_ff->p_atomff), parmset );  needed to use parametrization file from MORT FF it will fail for pure residues 

			if( ff_type == ForceFieldType::AMOEBA )
			{
				mort::atomiter_t aitr = pmol->atom_begin();
				for( ; aitr != pmol->atom_end(); ++aitr )
				{
					try
					{
						mort::prmtop::amoeba::create_atomic_frame( *aitr, *(p_mort_ff->p_poleff), chirials, regulars );
					}
					catch( std::exception& ex1)
					{
						std::string at_name = aitr->get_s(mort::NAME);
						int at_id = aitr->get_i(mort::POLTYPEID);
						PrintLog(" Error:  %s \n",ex1.what() );
						PrintLog(" Cannot make a frame for atom = %s  pol_id = %d of residue %s \n",
							at_name.c_str(), at_id, (*ritr_mort).name().c_str() );
					}
				}
			}
		}
		catch( std::exception& ex )
		{
			PrintLog(" Error to parametrize: residue %s \n",(*ritr_mort).name().c_str());
			PrintLog(" Error:  %s \n",ex.what() );
			continue;
		}

		std::string res_name_templ = "###";
		if( p_res_templ ) res_name_templ = p_res_templ->GetHostMol()->GetObjName();
		os << "<resff rname=\"" << res_name_templ << "\" ffname=\"" << ff_type.label() << "\" >" << std::endl;

		AtomIteratorAtomGroup aitr_hlm(p_res_templ);
		HaAtom* aptr = aitr_hlm.GetFirstAtom();

		int at_idx = 0;
		
		mort::atomiter_t aitr = pmol->atom_begin();
		for( ; aitr != pmol->atom_end(); ++aitr )
		{
			int at_id = aitr->get_i(mort::ID);
			std::string at_name_amber;
			std::string at_type;
			double      at_chg   = -10000.0;
			int         pol_type_id = -100;
			int         elem     = -100;
			double      polar    = -1000.0;
			mort::numvec pole;

			int bisect = FALSE;
			int at_id_2 = -100;
			int at_id_3 = -100;
			int at_id_4 = -100;
			std::string at_name_2;
			std::string at_name_3;
			std::string at_name_4;

			int i;
			int nreg  = 0;
			int nchir = 0;
			for( i = 0; i < regulars.size(); i++)
			{
				if( regulars[i][0] != at_id ) continue;
				if( nreg == 0 ) 
				{
					at_id_2 = regulars[i][3];
					if( regulars[i][4] == 2) bisect = TRUE;
				}
				if( nreg == 1 ) at_id_3 = regulars[i][3];
				nreg++;
			}
			for( i = 0; i < chirials.size(); i++)
			{
				if( chirials[i][0] != at_id ) continue;
				if( nchir == 0 ) at_id_4 = chirials[i][1];
				nchir++;
			}

			mort::atomiter_t aitr2 = pmol->atom_begin();
			for( ; aitr2 != pmol->atom_end(); ++aitr2 )
			{
				int id_m = aitr2->get_i(mort::ID);
				if( at_id_2 == id_m ) at_name_2 = aitr2->get_s(mort::NAME);
				if( at_id_3 == id_m ) at_name_3 = aitr2->get_s(mort::NAME);
				if( at_id_4 == id_m ) at_name_4 = aitr2->get_s(mort::NAME);
			}

			if( !at_name_2.empty() ) at_name_2 = MMForceField::GetAtNameFromAmber( at_name_2, res_name_templ);
			if( !at_name_3.empty() ) at_name_3 = MMForceField::GetAtNameFromAmber( at_name_3, res_name_templ);
			if( !at_name_4.empty() ) at_name_4 = MMForceField::GetAtNameFromAmber( at_name_4, res_name_templ);

			if( (*aitr).has_s(mort::NAME) ) at_name_amber  =  (*aitr).get_s(mort::NAME);
			if( (*aitr).has_i(mort::POLTYPEID) ) pol_type_id =  (*aitr).get_i(mort::POLTYPEID); 
			if( (*aitr).has_s(mort::TYPE) ) at_type  =  (*aitr).get_s(mort::TYPE);
			if( (*aitr).has_d(mort::PCHG) ) at_chg   =  (*aitr).get_d(mort::PCHG);
			if( (*aitr).has_i(mort::ELEMENT) ) elem     =  (*aitr).get_i(mort::ELEMENT); 
			if( (*aitr).has_v(mort::POLE) ) pole = aitr->get_v(mort::POLE);

			if( ff_type == ForceFieldType::AMOEBA )
			{
				mort::morf_t pol_atm_templ( *(p_mort_ff->p_poleff), mort::ATOM, pol_type_id-1 );
				polar   = pol_atm_templ.get_d(mort::POLAR);
			}

			std::string harlem_at_name = "###";
			std::string harlem_at_name_amber = "###";

			if( aptr ) 
			{
				harlem_at_name = aptr->GetName();
			}

			std::string at_name_conv_from_amber = MMForceField::GetAtNameFromAmber( at_name_amber, res_name_templ );

			if( harlem_at_name != at_name_conv_from_amber )
			{
				PrintLog("Amber residue: %s  harlem_res_name_templ = %s \n",amber_res_name.c_str(), res_name_templ.c_str());
				PrintLog(" at_idx= %d  at_name_conv_from_amber = %s  is not equal to harlem_at_name = %s \n",
					at_idx, at_name_conv_from_amber.c_str(), harlem_at_name.c_str()); 
			}

			os << "  <atomff " ;
			os << " name=\""   << harlem_at_name << "\""; 
			//				if( elem     > -100 )    os << " elem=\""     << elem << "\""; 
			if( pol_type_id > -100 ) os << " pol_type_id=\"" << pol_type_id << "\""; 
			os << " ff_type=\""     << at_type << "\""; 
			if( at_chg   > -1000.0 ) os << " chg=\""      << at_chg << "\"";
			os << " >";

			if( nreg > 0 )
			{
				os << std::endl;
				os << "      <frame ";
				if( bisect )    os << "bisect=\"1\" ";
				if( nchir > 0 ) os << "chiral=\"1\" ";
				os << "> ";
				os << harlem_at_name << "  " << at_name_2 << " " << at_name_3 << " " ;
				if( nchir > 0 ) os << at_name_4 << " ";
				os << "</frame>" << std::endl;
			}
			if( pole.size() > 0 ) 
			{
				os << "      <dipole>";
				sprintf(buf," %10.5f %10.5f %10.5f ",pole[0],pole[1],pole[2]);
				os << buf << "</dipole>" << std::endl;
				os << "      <qpole> " << std::endl;
				sprintf(buf,"       %10.5f \n       %10.5f %10.5f \n       %10.5f %10.5f %10.5f \n",
					pole[3],pole[4],pole[5],pole[6],pole[7],pole[8]);	    
				os << buf;
				os << "      </qpole>" << std::endl;
			}
			if( polar > -100.0 ) 
			{
				sprintf(buf,"      <polar> %10.5f </polar>",polar);
				os << buf << std::endl;
			}
			os << "  </atomff>" << std::endl;
			aptr = aitr_hlm.GetNextAtom();
			at_idx++;
		}
		
		std::set<int, less<int> > res_at_id;
		mort::atomiter_t aitr_rm = ritr_mort->atom_begin();
		for( ; aitr_rm !=  ritr_mort->atom_end(); ++aitr_rm )
		{
			int at_id = aitr_rm->get_i(mort::ID);
			res_at_id.insert(at_id);
		}

		os << "<improper_dihedrals>" << std::endl;
		mort::impriter_t impr_itr = mort_res_mol.impr_begin();
		for( ; impr_itr != mort_res_mol.impr_end(); ++impr_itr )
		{
			mort::atom_range ar = impr_itr->atoms();
			if( ar.size() != 4 )
			{
				PrintLog(" Invalid improper angle \n");
				continue;
			}
			int n_in_res = 0;
			StrVec at_str(4);
			for( i = 0; i < 4; i++)
			{
				int at_id = ar[i].get_i(mort::ID);
				if( res_at_id.count(at_id) > 0 ) 
				{
					n_in_res++;
					at_str[i] = ar[i].name();
					at_str[i] = MMForceField::GetAtNameFromAmber( at_str[i], res_name_templ);
				}
				else
				{
					std::string res_name_cnt = ar[i].resd().name();
					if( ritr_mort_prev != mort_res_mol.resd_end() && ar[i].resd() == (*ritr_mort_prev) )
					{
						res_name_cnt = "PREV";
					}
					if( ritr_mort_next != mort_res_mol.resd_end() && ar[i].resd() == (*ritr_mort_next) )
					{
						res_name_cnt = "NEXT";
					}
					at_str[i] = res_name_cnt + "." + ar[i].name();
				}
			}
			if( n_in_res < 2 ) continue;
			os << "   " << at_str[0] << "  " << at_str[1] << "   " << at_str[2]; 
			os << "   " << at_str[3] << std::endl;
		}
		os << "</improper_dihedrals>" << std::endl;

		os << "</resff>" << std::endl;
		ritr_mort_prev = ritr_mort;
	}

	os << harlem::HarlemDataFooter() << endl;
	return TRUE;
}
