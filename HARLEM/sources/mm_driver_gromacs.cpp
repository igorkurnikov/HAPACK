/*!  \file mm_driver_gromacs.cpp

    Molecular Mechanics simulations using GROMACS package  

    \author Igor Kurnikov 
    \date 2010-

*/
#include <mpi.h>
#include <math.h>
#include <tinyxml.h>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include "harlemapp.h"
#include "hamolecule.h"
#include "hamolset.h"
#include "hamolmech.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "mm_driver_gromacs.h"

MMDriverGromacs::MMDriverGromacs(HaMolMechMod* p_mm_mod_new)
{
	p_mm_mod = p_mm_mod_new;
	p_mm_model    = p_mm_mod->p_mm_model;
	pmset = p_mm_mod->GetMolSet();
	std::string prefix = pmset->GetName() + std::string("_gmx");
	this->SetFileNamesWithPrefix(prefix);

	to_save_input_files = TRUE;
}

MMDriverGromacs::~MMDriverGromacs()
{

}

std::set<std::string> MMDriverGromacs::std_gmx_mols = { "HOH","WAT","SOL","NA+","CL-","MG2+" };

void MMDriverGromacs::PartitionAtomsToMolecules()
{
	gmx_mol_partition.clear();

	ChainIteratorMolSet ch_itr(this->pmset);
	for (HaChain* p_ch = ch_itr.GetFirstChain(); p_ch; p_ch = ch_itr.GetNextChain())
	{
		HaMolecule* p_mol = p_ch->GetHostMol();

		std::string ch_str(1, p_ch->ident);
		if (ch_str == " ") ch_str = "0";
		std::string mol_name_gmx = p_mol->GetName() + std::string("_") + ch_str;

		gmx_mol_partition.push_back(AtomGroup());
		AtomGroup* p_gmx_mol = &gmx_mol_partition.back();
		p_gmx_mol->SetID(mol_name_gmx);

		ResidueIteratorChain ritr_ch(p_ch);
		for (HaResidue* pres = ritr_ch.GetFirstRes(); pres; pres = ritr_ch.GetNextRes())
		{	
			for (HaAtom* aptr: *pres)
			{
				p_gmx_mol->push_back(aptr);
			}

			std::string res_name = pres->GetName();
			if (std_gmx_mols.count(res_name) > 1)
			{
				p_gmx_mol->SetID(res_name);
				gmx_mol_partition.push_back(AtomGroup());
				p_gmx_mol = &gmx_mol_partition.back();
				p_gmx_mol->SetID(mol_name_gmx);
			}
		}

		if (p_gmx_mol->size() == 0) gmx_mol_partition.pop_back();
	}
}

void MMDriverGromacs::SetFileNamesWithPrefix(std::string prefix)
{
	inp_fname = prefix + ".mdp";
	top_fname = prefix + ".top";
	init_crd_fname = prefix + "_init.gro";
	restr_crd_fname = prefix + "_restr.gro";
	trj_fname = prefix + ".trr";
	ene_fname = prefix + ".ene";
}

int MMDriverGromacs::SaveAllInpFiles()
{	
	PrintLog("Save GROMACS mdp file %s\n", inp_fname.c_str());
	SaveMdpFile();
	PrintLog("Save GROMACS top file %s\n", top_fname.c_str());
	SaveGromacsTopFile();
	PrintLog("Save GROMACS Init Crd file %s \n", init_crd_fname.c_str());
	pmset->SaveGROFile(init_crd_fname.c_str());
	PrintLog("Save GROMACS Restr Crd file %s \n", restr_crd_fname.c_str());
	pmset->SaveGROFile(restr_crd_fname.c_str());
	to_save_input_files = FALSE;
	return TRUE;
}

int MMDriverGromacs::SaveMdpFile()
{
	if (p_mm_model->to_init_mm_model) p_mm_mod->InitMolMechModel();
	std::ofstream os(this->inp_fname);
	if (os.fail())
	{
		PrintLog("Error in MMDriverGromacs::SaveMdpFile() \n");
		PrintLog("Can't create file %s \n", this->inp_fname.c_str());
		return FALSE;
	}
	return SaveMdpToStream(os);
}

int MMDriverGromacs::SaveGromacsTopFile()
{
	if(p_mm_model->to_init_mm_model) p_mm_mod->InitMolMechModel();

	PartitionAtomsToMolecules();

	std::ofstream os( this->top_fname );
	if(os.fail()) 
	{
		PrintLog("Error in MMDriverGromacs::SaveTopFile() \n");
		PrintLog("Can't create file %s \n", this->top_fname.c_str() );
		return FALSE;
	}
	return SaveGromacsTopToStream(os);
}

int MMDriverGromacs::SaveMdpToStream(std::ostream& os)
{
	int i;
	char buf[256];

	try
	{
		os << ";  GROMACS mdp file created by HARLEM" << std::endl << std::endl;

		if (p_mm_mod->run_type == p_mm_mod->run_type.MIN_RUN)
		{
			if (p_mm_mod->min_type.value() == p_mm_mod->min_type.STEEPEST_DESCENT)
			{
				os << " integrator=steep          ; Steepest Decent minimization " << std::endl;
			}
			else if (p_mm_mod->min_type.value() == p_mm_mod->min_type.CONJ_GRAD)
			{
				os << " integrator=cg            ; Conjugate Gradient minimization " << std::endl;
				os << " nstcgsteep= " << p_mm_mod->num_steep_descent_steps << "        ; frequency of performing 1 steepest descent step while doing conjugate gradient energy minimization " << std::endl;

			}
			else
			{
				os << " integrator=steep " << std::endl;
			}
			os << " emstep= " << p_mm_mod->init_min_step << "       ; initial step-size  " << std::endl;
			os << " emtol= " << p_mm_mod->grad_cnvrg_val << "       ; Gradient convergence tolerance   " << std::endl;
			os << " nsteps= " << p_mm_mod->max_num_minim_steps << "       ; Maximum number of minimization steps " << std::endl;
		}
		else if (p_mm_mod->run_type == p_mm_mod->run_type.MD_RUN)
		{
			if (p_mm_mod->temp_control_method.CONST_TEMP_LANGEVIN)
			{
				os << " integrator=sd            ; Molecular Dynamics with Langevin Integrator " << std::endl;
			}
			else
			{
				os << " integrator=md            ; Molecular Dynamics with leap-frog allgorithm " << std::endl;
			}

			os << " dt=" << p_mm_mod->md_time_step << "         ; time step for integration " << std::endl;
			if (fabs(p_mm_mod->start_time) > 0.0000001)
			{
				os << " tinit=" << p_mm_mod->start_time << "        ; starting time for your run (ps)  " << std::endl;
			}
			os << " nsteps=" << p_mm_mod->num_md_steps  << "        ; maximum number of steps to integrate " << std::endl;
			if (p_mm_mod->remove_rb_motion_freq != 100)
			{
				os << " nstcomm=" << p_mm_mod->remove_rb_motion_freq << "       ; frequency for center of mass motion removal " << std::endl;
			}
		}
		else
		{
			throw std::runtime_error((std::string)("Run type") + p_mm_mod->run_type.label() + " is not supported in GROMACS ");
		}

		// Holonomic bond and angle constraints:

		if (p_mm_mod->shake_constr == p_mm_mod->shake_constr.NO_SHAKE)
		{
			os << " constraints=none " << "   ;   No bonds converted to constraints. " << std::endl;
		}
		else if (p_mm_mod->shake_constr == p_mm_mod->shake_constr.H_ATOM_SHAKE)
		{
			os << " constraints=h-bonds " << "   ;  Convert the bonds with H-atoms to constraints. " << std::endl;
			os << " constraint-algorithm=LINCS " << "   ;  LINear Constraint Solver. The accuracy in set with lincs-order,  " << std::endl;
			os << " lincs-order=6 " << "   ; Highest order in the expansion of the constraint coupling matrix " << std::endl;
		}
		else if (p_mm_mod->shake_constr == p_mm_mod->shake_constr.ALL_BOND_SHAKE)
		{
			os << " constraints=all-bonds " << "   ; Convert all bonds to constraints. " << std::endl;
			os << " constraint-algorithm=LINCS " << "   ;  LINear Constraint Solver. The accuracy in set with lincs-order,  " << std::endl;
		}
		else
		{
			throw std::runtime_error((std::string)("Constraint method ") + p_mm_mod->shake_constr.label() + " is not supported in GROMACS ");
		}
		// fprintf(fp, " tol=%12.7f, ", p_mm_mod->shake_tol );


		// Print out control:
		os << std::endl << "; MM Run Output control " << std::endl << std::endl;

		if (p_mm_mod->wrt_log_freq != 1000)
		{
			os << " nstlog= " << p_mm_mod->wrt_log_freq << "        ;  number of steps that elapse between writing energies to the log file  " << std::endl;
		}
		if (p_mm_mod->traj_wrt_format = p_mm_mod->traj_wrt_format.XTC)
		{
			if (p_mm_mod->wrt_coord_freq > 0)
			{
				os << " nstxout-compressed= " << p_mm_mod->wrt_coord_freq << "        ;  number of steps that elapse between writing position coordinates using lossy compression (xtc file) " << std::endl;
			}
		}
		else
		{
			if (p_mm_mod->wrt_coord_freq > 0)
			{
				os << " nstxout= " << p_mm_mod->wrt_coord_freq << "        ;  number of steps that elapse between writing coordinates to the output trajectory file (trr file) " << std::endl;
			}
		}
		if (p_mm_mod->wrt_vel_freq > 0)
		{
			os << " nstvout= " << p_mm_mod->wrt_vel_freq << "        ;  number of steps that elapse between writing velocities to the output trajectory file  " << std::endl;
		}
		if( p_mm_mod->wrt_ener_freq != 1000 )
		{ 
			os << " nstenergy= " << p_mm_mod->wrt_ener_freq   << "        ;  number of steps that elapse between writing energies to energy file  " << std::endl;
		}

		// fprintf(fp, " ntwr=%d, ", p_mm_mod->wrt_rstrt_freq);
		// fprintf(fp, " ntwprt=%d, ", p_mm_mod->limit_wrt_atoms);
		// fprintf(fp, " ntf=%d, ", p_mm_model->omit_interactions.value());

		os << std::endl << "; Pressure Control " << std::endl << std::endl;

		if ( p_mm_mod->period_bcond == p_mm_mod->period_bcond.NO_PERIODICITY )
		{
			os << " pbc=no " << "        ;  Use no periodic boundary conditions, ignore the box  " << std::endl;
		}
		else if (p_mm_mod->period_bcond == p_mm_mod->period_bcond.CONST_PRES || p_mm_mod->period_bcond == p_mm_mod->period_bcond.CONST_VOL)
		{
			os << " pbc=xyz " << "        ;  Use periodic boundary conditions in all directions  " << std::endl;
			if (p_mm_mod->run_type == p_mm_mod->run_type.MD_RUN)
			{
				if (p_mm_mod->period_bcond == p_mm_mod->period_bcond.CONST_VOL || p_mm_mod->pressure_reg_method == p_mm_mod->pressure_reg_method.NO_CRD_SCALING)
				{
					os << " pcoupl=no " << "        ;  No pressure coupling. This means a fixed box size  " << std::endl;
				}
				else if (p_mm_mod->period_bcond == p_mm_mod->period_bcond.CONST_PRES)
				{
					// os << " pcoupl=Berendsen " << "        ;  Exponential relaxation pressure coupling  " << std::endl;
					// pressure regulation
					if (p_mm_mod->pressure_reg_method == p_mm_mod->pressure_reg_method.ISOTROP_CRD_SCALING)
					{
						// os << " pcoupltype=isotropic " << "    ;  Isotropic pressure coupling  " << std::endl;
						os << " ref-p= " << p_mm_mod->ref_pressure << "    ;  The reference pressure for coupling  " << std::endl;
						os << " tau-p= " << p_mm_mod->press_relax_time << "    ;  The time constant for pressure coupling  " << std::endl;
						if (fabs(p_mm_mod->compressibility - 44.6) > 1.0)
						{
							os << " compressibility= " << p_mm_mod->compressibility * 1.0E-6 << "    ;  The compressibility ( bar^-1)  " << std::endl;
						}
					}
					else
					{
						throw std::runtime_error((std::string)("Pressure control method ") + p_mm_mod->pressure_reg_method.label() + " is not supported in GROMACS ");
					}
				}
				else
				{
					throw std::runtime_error((std::string)("Periodic boundary condition: ") + p_mm_mod->period_bcond.label() + " is not supported in GROMACS ");
				}
			}
		}

		os << std::endl << "; Non-bonded interactions " << std::endl << std::endl;
 		if (p_mm_model->electr_method == p_mm_model->electr_method.PME_METHOD )
		{
			os << " coulombtype=PME " << "            ;  Fast Smooth Particle-Mesh Ewald (SPME) electrostatics  " << std::endl;
		}
		else if (p_mm_model->electr_method == p_mm_model->electr_method.SCREENED_COULOMB )
		{
			os << " coulombtype=Cut-off " << "        ; Plain cut-off with pair list radius rlist and Coulomb cut-off rcoulomb  " << std::endl;
			if (fabs(p_mm_model->diel_const - 1.0) > 0.00001)
			{
				os << " epsilon-r= " << p_mm_model->diel_const << "        ; The relative dielectric constant  " << std::endl;
			}
		}
		else
		{
			throw std::runtime_error((std::string)("Electrostatics method ") + p_mm_model->electr_method.label() + " is not supported in GROMACS ");
		}

		os << " cutoff-scheme=Verlet " << "        ;  Generate a pair list with buffering " << std::endl;
		os << " rlist= " << boost::format("%8.3f") % (p_mm_model->nb_cut_dist * 0.1) << "        ;  Cut-off distance for the short-range neighbor list (nm). " << std::endl;

		os << " rcoulomb= " << (p_mm_model->nb_cut_dist * 0.1) << "        ;  distance for the Coulomb cut-off (nm) " << std::endl;
		os << " vdwtype=Cut-off " << "        ;  Plain cut-off  " << std::endl;
		os << " rvdw= " << (p_mm_model->nb_cut_dist * 0.1) << "        ;  distance for the LJ cut-off (nm) " << std::endl;

		// temperature regulation:
		if (p_mm_mod->run_type == p_mm_mod->run_type.MD_RUN)
		{
			os << std::endl << "; Temperature Control " << std::endl << std::endl;

			if (p_mm_mod->temp_control_method == p_mm_mod->temp_control_method.CONST_ENE_MD)
			{
				os << " tcoupl=no " << "        ;  No temperature coupling." << std::endl;
			}
			else if (p_mm_mod->temp_control_method == p_mm_mod->temp_control_method.CONST_TEMP_LANGEVIN || p_mm_mod->temp_control_method == p_mm_mod->temp_control_method.CONST_TEMP_BERENDSEN)
			{
				if (p_mm_mod->temp_control_method == p_mm_mod->temp_control_method.CONST_TEMP_BERENDSEN)
				{
					os << " tcoupl=berendsen " << "        ;  Temperature coupling with a Berendsen thermostat." << std::endl;
				}
				os << " tc-grps= system" << "              ; groups to couple to separate temperature baths  " << std::endl;
				os << " ref-t= " << p_mm_mod->ref_temp << "        ;  reference temperature for coupling (K) " << std::endl;
				os << " tau-t= " << p_mm_mod->langevin_dump_const << "        ;  time constant for coupling (ps) " << std::endl;
			}
			else
			{
				throw std::runtime_error((std::string)("Temperature control method ") + p_mm_mod->temp_control_method.label() + " is not supported in GROMACS ");
			}

			if (fabs(p_mm_mod->init_temp - 0.0) > 0.00001)
			{
				os << " gen-vel=yes " << "   ;  Generate velocities in gmx grompp according to a Maxwell distribution " << std::endl;
				os << " gen-temp= " << p_mm_mod->init_temp << std::endl;
			}
		}


		// Frozen and restrained atoms: To be expanded here

		// fprintf(fp, "ibelly=%d, ", p_amber_model->ibelly);
		// int ntr_loc = 0; if (p_amber_model->natc > 0) ntr_loc = 1;
		// fprintf(fp, "ntr=%d, ", ntr_loc);
		// fprintf(fp, "\n");

	}
	catch (const std::exception& ex)
	{
		PrintLog(" Error in MMDriverGromacs::SaveMdpToStream() \n");
		PrintLog("%s\n", ex.what());
		return FALSE;
	}
	return FALSE;
}

int MMDriverGromacs::SaveGromacsTopToStream(std::ostream& os)
{
	int i;
	char buf[256];
	
	std::string ffdb_prefix;

	if (p_mm_model->ff_type == p_mm_model->ff_type.AMBER_94)  ffdb_prefix = "amber94.ff";
	else if (p_mm_model->ff_type == p_mm_model->ff_type.AMBER_99_SB) ffdb_prefix = "amber99sb.ff";
	else if (p_mm_model->ff_type == p_mm_model->ff_type.AMBER_03)  ffdb_prefix = "amber03.ff";
	else if (p_mm_model->ff_type == p_mm_model->ff_type.AMBER_10)  ffdb_prefix = "amber03.ff";  // TMP IGOR
	else ffdb_prefix = "amber99sb.ff";

	std::set<std::string> gmx_mol_defined;

	for (std::string ss : this->std_gmx_mols)
		gmx_mol_defined.insert(ss);

	std::list<std::string> system_mol_list;  // list of molecules of the system  and their counts  

	try
	{
		os << "; GROMACS top created by HARLEM" << std::endl;
		os << std::endl;

		os << "#include \"" << ffdb_prefix << "/forcefield.itp\" " << std::endl;
		
		std::string last_group = "";
		int last_group_count = 0;

		for (AtomGroup& group : this->gmx_mol_partition)
		{
			std::string grp_id = group.GetID();
			if (gmx_mol_defined.count(grp_id) > 0)
			{
				if (grp_id == last_group)
				{
					last_group_count++;
					continue;
				}
				continue;
			}
			if (!last_group.empty())
				system_mol_list.push_back((boost::format("%12s %6d") % last_group % last_group_count).str());
			last_group = grp_id;
			last_group_count = 1;

			// save GROMACS molecule definition
			
			AtomIntMap at_idx_map;
			int idx = 0;
			for (HaAtom* aptr : group)
			{
				at_idx_map[aptr] = idx;
				++idx;
			}

			os << "  " << std::endl;
			os << "[ moleculetype ]" << std::endl;
			os << ";name            nrexcl" << std::endl;
			os << " " << grp_id << "           3 " << std::endl;
			os << "  " << std::endl;

			SaveAtomsToStream(os, group, at_idx_map);
			SaveBondsToStream(os, group, at_idx_map);
			Save14PairsToStream(os, group, at_idx_map);
			SaveAnglesToStream(os, group, at_idx_map);
			SaveDihedralsToStream(os, group, at_idx_map);

			gmx_mol_defined.insert(grp_id);

			last_group = grp_id;
			last_group_count = 1;
		}

		if (!last_group.empty()) 
			system_mol_list.push_back((boost::format("%12s %6d") % last_group % last_group_count).str());

		os << std::endl << "; Include water topology " << std::endl;
		os << "#include \"" << ffdb_prefix << "/tip3p.itp\" " << std::endl;

		os << std::endl << "; Include topology for ions" << std::endl;
		os << "#include \"" << ffdb_prefix << "/ions.itp\" " << std::endl;

		os << std::endl << "[ system ]" << std::endl;
		os << "; Name " << std::endl;
		os << "System" << std::endl;

		os << std::endl << "[ molecules ]" << std::endl;
		os << "; Compound        #mols" << std::endl;
		for (std::string mol_str : system_mol_list)
			os << mol_str << std::endl;
		os << "  ";

	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in MMDriverGromacs::SaveTopToStream() \n");
		PrintLog("%s\n",ex.what());
		return FALSE;
	}
	return FALSE;
}

int MMDriverGromacs::SaveAtomTypesToStream(std::ostream& os)
{
		char buf[256];
		std::set<std::string> at_ff_names;

		os << "[ atomtypes ]" << std::endl;
		os << ";name   bond_type     mass     charge   ptype   sigma         epsilon       Amb" << std::endl;
//		os << ";   c        c           0.00000  0.00000   A     3.39967e-01   3.59824e-01 ; 1.91  0.0860" << std::endl;

		int i;
		for( i = 0; i < p_mm_model->Atoms.size(); i++ )
		{
			std::string ff_s = p_mm_model->Atoms[i]->GetFFSymbol();
			boost::trim(ff_s);
			if( ff_s.empty() ) continue;
			if( at_ff_names.count(ff_s) > 0 ) continue;
			at_ff_names.insert(ff_s);

			double vdw_rad = p_mm_model->Atoms[i]->vdw_rad;
			double ew = p_mm_model->Atoms[i]->ew;

			double sigma   = 2*vdw_rad*0.1/pow(2.0,1.0/6.0);
			double epsilon = 4.184*ew;
			
			sprintf( buf, "%4s %4s    0.00000    0.00000  A  %14.5e %14.5e  ;  %9.5f  %9.6f ",
				           ff_s.c_str(),ff_s.c_str(), sigma, epsilon, vdw_rad, ew );
			os << buf << std::endl;
		}
		os << "  " << std::endl;
		return TRUE;
}

int MMDriverGromacs::SaveAtomsToStream(std::ostream& os, AtomGroup& group, AtomIntMap& at_idx_map)
{
	char buf[256];
		
	os << "[ atoms ]" << std::endl;
	os << ";   nr  type   resi  res  atom    cgnr     charge      mass       ; qtot   bond_type" << std::endl;
	double qtot = 0.0;
	for( HaAtom* aptr: group )
	{
		std::string atn = aptr->GetName();
		double ch = aptr->GetCharge();
		double mass = aptr->GetMass();
		HaResidue* pres = aptr->GetHostRes();
		std::string ff_s = aptr->GetFFSymbol();
		std::string res_name = pres->GetName();
		int resi = pres->GetSerNo();
		qtot += ch;
		int idx = at_idx_map[aptr] + 1;
		sprintf(buf,"%6d %5s%6d %5s %5s %6d %12.6f %12.6f ; qtot %7.3f ",idx,ff_s.c_str(),resi,res_name.c_str(),atn.c_str(),idx,ch,mass,qtot);
		os << buf << std::endl;
	}
	os << "  " << std::endl;

	return TRUE;
}

int MMDriverGromacs::SaveBondsToStream(std::ostream& os, AtomGroup& group, AtomIntMap& at_idx_map)
{
	char buf[256];
	os << "[ bonds ]" << std::endl;
	os << ";   ai     aj funct   r             k   " << std::endl;

	set<MMBond, less<MMBond> >::iterator mbitr = p_mm_model->MBonds.begin();
	for(; mbitr != p_mm_model->MBonds.end(); mbitr++)
	{
		MMBond& bnd = (MMBond&) *mbitr;		
		HaAtom* aptr1 = (HaAtom*) bnd.pt1;
		HaAtom* aptr2 = (HaAtom*) bnd.pt2;

		if( at_idx_map.count(aptr1) == 0 ) continue;
		if( at_idx_map.count(aptr2) == 0 ) continue;

		int idx1 = at_idx_map[aptr1]+1;  // convert to 1-based index
		int idx2 = at_idx_map[aptr2]+1;  // convert to 1-based index

		double r0 = bnd.r0 * 0.1; // Ang -> nm
		double fc = bnd.fc*4.184*4.184*100; // kcal/mol/Ang^2  ->kJ/nm^2

		std::string at_lbl_1 = aptr1->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_2 = aptr2->GetRef(HaAtom::ATOMREF_STD);

		sprintf(buf,"%6d  %6d   1 %14.4e%14.4e ; %s - %s ",idx1,idx2, r0,fc, at_lbl_1.c_str(), at_lbl_2.c_str());
		os << buf << std::endl;
	}
	os << "  " << std::endl;
	return TRUE;
}

int MMDriverGromacs::Save14PairsToStream(std::ostream& os, AtomGroup& group, AtomIntMap& at_idx_map)
{
	char buf[256];
	os << "[ pairs ]" << std::endl;
	os << ";   ai     aj    funct" << std::endl;

	vector<MMDihedral>::iterator ditr;
	for( ditr = p_mm_model->Dihedrals.begin(); ditr != p_mm_model->Dihedrals.end(); ditr++)
	{
		HaAtom* aptr1 = (HaAtom*) (*ditr).pt1;
		HaAtom* aptr2 = (HaAtom*) (*ditr).pt4;

		if( at_idx_map.count(aptr1) == 0 ) continue;
		if( at_idx_map.count(aptr2) == 0 ) continue;

		int idx1 = at_idx_map[aptr1]+1;  // convert to 1-based index
		int idx2 = at_idx_map[aptr2]+1;  // convert to 1-based index

		std::string at_lbl_1 = aptr1->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_2 = aptr2->GetRef(HaAtom::ATOMREF_STD);

		sprintf(buf,"%6d  %6d  1 ; %s - %s ",idx1,idx2, at_lbl_1.c_str(), at_lbl_2.c_str());
		os << buf << std::endl;
	}
	os << "  " << std::endl;
	return TRUE;
}

int MMDriverGromacs::SaveAnglesToStream(std::ostream& os, AtomGroup& group, AtomIntMap& at_idx_map)
{
	char buf[256];
	os << "[ angles ]" << std::endl;
	os << ";   ai     aj     ak    funct   theta         cth " << std::endl;

	set<MMValAngle, less<MMValAngle> >::iterator vaitr = p_mm_model->ValAngles.begin();
	for(; vaitr != p_mm_model->ValAngles.end(); vaitr++)
	{
		MMValAngle& ang = (MMValAngle&)(*vaitr);
		
		double a0 = ang.a0;
		double fc = ang.fc * 4.184;

		HaAtom* aptr1 = (HaAtom*) ang.pt1;
		HaAtom* aptr2 = (HaAtom*) ang.pt2;
		HaAtom* aptr3 = (HaAtom*) ang.pt3;

		if( at_idx_map.count(aptr1) == 0 ) continue;
		if( at_idx_map.count(aptr2) == 0 ) continue;
		if( at_idx_map.count(aptr3) == 0 ) continue;

		int idx1 = at_idx_map[aptr1]+1;  // convert to 1-based index
		int idx2 = at_idx_map[aptr2]+1;  // convert to 1-based index
		int idx3 = at_idx_map[aptr3]+1;  // convert to 1-based index

		std::string at_lbl_1 = aptr1->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_2 = aptr2->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_3 = aptr3->GetRef(HaAtom::ATOMREF_STD);

		sprintf(buf,"%6d %6d %6d  1  %14.4e%14.4e; %s - %s - %s ",idx1,idx2,idx3, a0,fc,at_lbl_1.c_str(), at_lbl_2.c_str(),at_lbl_3.c_str());
		os << buf << std::endl;
	}
	os << "  " << std::endl;
	return TRUE;
}

int MMDriverGromacs::SaveDihedralsToStream(std::ostream& os, AtomGroup& group, AtomIntMap& at_idx_map)
{
	char buf[256];
	
	os << "[ dihedrals ] ; proper " << std::endl;
	os << ";    i      j      k      l   func   phase     kd      pn" << std::endl;

	int i;
	std::vector<MMDihedral>::iterator ditr;		
	for(ditr = p_mm_model->Dihedrals.begin(); ditr != p_mm_model->Dihedrals.end(); ditr++)
	{
		MMDihedral& dih = (MMDihedral&)(*ditr);
		
		std::vector<double> pk = dih.pk;
		std::vector<double> pn = dih.pn;
		std::vector<double> phase = dih.phase;

		HaAtom* aptr1 = (HaAtom*) dih.pt1;
		HaAtom* aptr2 = (HaAtom*) dih.pt2;
		HaAtom* aptr3 = (HaAtom*) dih.pt3;
		HaAtom* aptr4 = (HaAtom*) dih.pt4;

		if( at_idx_map.count(aptr1) == 0 ) continue;
		if( at_idx_map.count(aptr2) == 0 ) continue;
		if( at_idx_map.count(aptr3) == 0 ) continue;
		if( at_idx_map.count(aptr4) == 0 ) continue;

		int idx1 = at_idx_map[aptr1]+1;  // convert to 1-based index
		int idx2 = at_idx_map[aptr2]+1;  // convert to 1-based index
		int idx3 = at_idx_map[aptr3]+1;  // convert to 1-based index
		int idx4 = at_idx_map[aptr4]+1;  // convert to 1-based index

		std::string at_lbl_1 = aptr1->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_2 = aptr2->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_3 = aptr3->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_4 = aptr4->GetRef(HaAtom::ATOMREF_STD);

		for( i = 0; i < pk.size(); i++)
		{
			pk[i] *= 4.184;
			sprintf(buf,"%6d %6d %6d %6d   1  %14.4e%14.4e %2.0f  ; %s - %s - %s - %s ",idx1,idx2,idx3,idx4, phase[i], pk[i], pn[i], 
				at_lbl_1.c_str(), at_lbl_2.c_str(),at_lbl_3.c_str(),at_lbl_4.c_str());
			os << buf << std::endl;
		}
	}
	os << "  " << std::endl;

	os << "[ dihedrals ] ; improper " << std::endl;		
	os << ";    i      j      k      l   func   phase     kd      pn" << std::endl;
	
	for(ditr = p_mm_model->ImprDihedrals.begin(); ditr != p_mm_model->ImprDihedrals.end(); ditr++)
	{
		MMDihedral& dih = (MMDihedral&)(*ditr);
		
		std::vector<double> pk = dih.pk;
		std::vector<double> pn = dih.pn;
		std::vector<double> phase = dih.phase;

		HaAtom* aptr1 = (HaAtom*) dih.pt1;
		HaAtom* aptr2 = (HaAtom*) dih.pt2;
		HaAtom* aptr3 = (HaAtom*) dih.pt3;
		HaAtom* aptr4 = (HaAtom*) dih.pt4;

		if( at_idx_map.count(aptr1) == 0 ) continue;
		if( at_idx_map.count(aptr2) == 0 ) continue;
		if( at_idx_map.count(aptr3) == 0 ) continue;
		if( at_idx_map.count(aptr4) == 0 ) continue;

		int idx1 = at_idx_map[aptr1]+1;  // convert to 1-based index
		int idx2 = at_idx_map[aptr2]+1;  // convert to 1-based index
		int idx3 = at_idx_map[aptr3]+1;  // convert to 1-based index
		int idx4 = at_idx_map[aptr4]+1;  // convert to 1-based index

		std::string at_lbl_1 = aptr1->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_2 = aptr2->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_3 = aptr3->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_4 = aptr4->GetRef(HaAtom::ATOMREF_STD);

		for( i = 0; i < pk.size(); i++)
		{
			pk[i] *= 4.184;
			sprintf(buf,"%6d %6d %6d %6d   1  %14.4e%14.4e %2.0f  ; %s - %s - %s - %s ",idx1,idx2,idx3,idx4, phase[i], pk[i], pn[i], 
				at_lbl_1.c_str(), at_lbl_2.c_str(),at_lbl_3.c_str(),at_lbl_4.c_str());
			os << buf << std::endl;
		}	
	}
			
	os << "  " << std::endl;
	return TRUE;
}
