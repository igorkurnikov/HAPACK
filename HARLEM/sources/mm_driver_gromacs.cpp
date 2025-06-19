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

std::set<std::string> MMDriverGromacs::std_gmx_mols = { "HOH","WAT","SOL","NA","CL","MG","NA+","CL-","MG2+"};

void MMDriverGromacs::PartitionAtomsToMolecules()
{
	gmx_mol_partition.clear();

	ChainIteratorMolSet ch_itr(this->pmset);
	for (HaChain* p_ch = ch_itr.GetFirstChain(); p_ch; p_ch = ch_itr.GetNextChain())
	{
		HaMolecule* p_mol = p_ch->GetHostMol();

		std::string ch_str(1, p_ch->ident);
		std::string mol_name_gmx = p_mol->GetName();
		if (p_mol->GetNChains() > 0)
		{
			if (ch_str == " ") ch_str = "0";
			mol_name_gmx = p_mol->GetName() + std::string("_") + ch_str;
		}

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
			if (std_gmx_mols.count(res_name) > 0)
			{
				p_gmx_mol->SetID(res_name);
				gmx_mol_partition.push_back(AtomGroup());
				p_gmx_mol = &gmx_mol_partition.back();
				p_gmx_mol->SetID(mol_name_gmx);
			}
		}

		if (p_gmx_mol->size() == 0) gmx_mol_partition.pop_back();
	}
	//PrintLog(" System partition into GROMACS molecules: \n");
	//for (AtomGroup& group : gmx_mol_partition)
	//	PrintLog(" %s : %6d atoms \n", group.GetID(), group.size() );

}

void MMDriverGromacs::SetFileNamesWithPrefix(std::string prefix)
{
	inp_fname = prefix + ".mdp";
	top_fname = prefix + ".top";
	init_crd_fname = prefix + "_init.gro";
	restr_crd_fname = prefix + "_restr.gro";
	trj_fname = prefix + ".trr";
	ene_fname = prefix + ".ene";
	run_fname = prefix + ".sh";
	tpr_fname = prefix + ".tpr";
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
	PrintLog("Save GROMACS run file %s\n", run_fname.c_str());
	SaveRunFile();
	to_save_input_files = FALSE;
	return TRUE;
}

int MMDriverGromacs::SaveMdpFile()
{
	if (p_mm_model->to_init_mm_model) p_mm_mod->InitMolMechModel();
	std::ofstream os(this->inp_fname, std::ios::binary);
	if (os.fail())
	{
		PrintLog("Error in MMDriverGromacs::SaveMdpFile() \n");
		PrintLog("Can't create file %s \n", this->inp_fname.c_str());
		return FALSE;
	}
	return SaveMdpToStream(os);
}

int MMDriverGromacs::SaveRunFile()
{
	if (p_mm_model->to_init_mm_model) p_mm_mod->InitMolMechModel();
	std::ofstream os(this->run_fname, std::ios::binary);
	if (os.fail())
	{
		PrintLog("Error in MMDriverGromacs::SaveRunFile() \n");
		PrintLog("Can't create file %s \n", this->run_fname.c_str());
		return FALSE;
	}
	os << "#!/bin/sh -f \n";
	os << boost::format("gmx grompp -f %s -p %s -c %s -r %s -o %s -maxwarn 2 ") %
		inp_fname % top_fname % init_crd_fname % restr_crd_fname % tpr_fname << " \n";
	os << boost::format("gmx mdrun -s %s ") % tpr_fname << " \n";
	
	return TRUE;
}

int MMDriverGromacs::SaveGromacsTopFile()
{
	if(p_mm_model->to_init_mm_model) p_mm_mod->InitMolMechModel();

	PartitionAtomsToMolecules();

	std::ofstream os( this->top_fname, std::ios::binary);
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
		os << ";  GROMACS mdp file created by HARLEM\n\n";

		if (p_mm_mod->run_type == p_mm_mod->run_type.MIN_RUN)
		{
			if (p_mm_mod->min_type.value() == p_mm_mod->min_type.STEEPEST_DESCENT)
			{
				os << " integrator=steep          ; Steepest Decent minimization \n";
			}
			else if (p_mm_mod->min_type.value() == p_mm_mod->min_type.CONJ_GRAD)
			{
				os << " integrator=cg            ; Conjugate Gradient minimization \n";
				os << " nstcgsteep= " << p_mm_mod->num_steep_descent_steps << "        ; frequency of performing 1 steepest descent step while doing conjugate gradient energy minimization \n";

			}
			else
			{
				os << " integrator=steep \n";
			}
			os << " emstep= " << p_mm_mod->init_min_step << "       ; initial step-size  \n";
			os << " emtol= " << p_mm_mod->grad_cnvrg_val << "       ; Gradient convergence tolerance   \n";
			os << " nsteps= " << p_mm_mod->max_num_minim_steps << "       ; Maximum number of minimization steps \n";
		}
		else if (p_mm_mod->run_type == p_mm_mod->run_type.MD_RUN)
		{
			if (p_mm_mod->temp_control_method.CONST_TEMP_LANGEVIN)
			{
				os << " integrator=sd            ; Molecular Dynamics with Langevin Integrator \n";
			}
			else
			{
				os << " integrator=md            ; Molecular Dynamics with leap-frog allgorithm \n";
			}

			os << " dt=" << p_mm_mod->md_time_step << "         ; time step for integration \n";
			if (fabs(p_mm_mod->start_time) > 0.0000001)
			{
				os << " tinit=" << p_mm_mod->start_time << "        ; starting time for your run (ps)  \n";
			}
			os << " nsteps=" << p_mm_mod->num_md_steps  << "        ; maximum number of steps to integrate \n";
			if (p_mm_mod->remove_rb_motion_freq != 100)
			{
				os << " nstcomm=" << p_mm_mod->remove_rb_motion_freq << "       ; frequency for center of mass motion removal \n";
			}
		}
		else
		{
			throw std::runtime_error((std::string)("Run type") + p_mm_mod->run_type.label() + " is not supported in GROMACS ");
		}

		// Holonomic bond and angle constraints:

		if (p_mm_mod->shake_constr == p_mm_mod->shake_constr.NO_SHAKE)
		{
			os << " constraints=none " << "   ;   No bonds converted to constraints. \n";
		}
		else if (p_mm_mod->shake_constr == p_mm_mod->shake_constr.H_ATOM_SHAKE)
		{
			os << " constraints=h-bonds " << "   ;  Convert the bonds with H-atoms to constraints. \n";
			os << " constraint-algorithm=LINCS " << "   ;  LINear Constraint Solver. The accuracy in set with lincs-order,  \n";
			os << " lincs-order=6 " << "   ; Highest order in the expansion of the constraint coupling matrix \n";
		}
		else if (p_mm_mod->shake_constr == p_mm_mod->shake_constr.ALL_BOND_SHAKE)
		{
			os << " constraints=all-bonds " << "   ; Convert all bonds to constraints. " << std::endl;
			os << " constraint-algorithm=LINCS " << "   ;  LINear Constraint Solver. The accuracy in set with lincs-order,  \n";
		}
		else
		{
			throw std::runtime_error((std::string)("Constraint method ") + p_mm_mod->shake_constr.label() + " is not supported in GROMACS ");
		}
		// fprintf(fp, " tol=%12.7f, ", p_mm_mod->shake_tol );


		// Print out control:
		os << std::endl << "; MM Run Output control \n\n";

		if (p_mm_mod->wrt_log_freq != 1000)
		{
			os << " nstlog= " << p_mm_mod->wrt_log_freq << "        ;  number of steps that elapse between writing energies to the log file  \n";
		}
		if (p_mm_mod->traj_wrt_format = p_mm_mod->traj_wrt_format.XTC)
		{
			if (p_mm_mod->wrt_coord_freq > 0)
			{
				os << " nstxout-compressed= " << p_mm_mod->wrt_coord_freq << "        ;  number of steps that elapse between writing position coordinates using lossy compression (xtc file) \n";
			}
		}
		else
		{
			if (p_mm_mod->wrt_coord_freq > 0)
			{
				os << " nstxout= " << p_mm_mod->wrt_coord_freq << "        ;  number of steps that elapse between writing coordinates to the output trajectory file (trr file) \n";
			}
		}
		if (p_mm_mod->wrt_vel_freq > 0)
		{
			os << " nstvout= " << p_mm_mod->wrt_vel_freq << "        ;  number of steps that elapse between writing velocities to the output trajectory file  \n";
		}
		if( p_mm_mod->wrt_ener_freq != 1000 )
		{ 
			os << " nstenergy= " << p_mm_mod->wrt_ener_freq   << "        ;  number of steps that elapse between writing energies to energy file  \n";
		}

		// fprintf(fp, " ntwr=%d, ", p_mm_mod->wrt_rstrt_freq);
		// fprintf(fp, " ntwprt=%d, ", p_mm_mod->limit_wrt_atoms);
		// fprintf(fp, " ntf=%d, ", p_mm_model->omit_interactions.value());

		os << std::endl << "; Pressure Control \n\n";

		if ( p_mm_mod->period_bcond == p_mm_mod->period_bcond.NO_PERIODICITY )
		{
			os << " pbc=no " << "        ;  Use no periodic boundary conditions, ignore the box  \n";
		}
		else if (p_mm_mod->period_bcond == p_mm_mod->period_bcond.CONST_PRES || p_mm_mod->period_bcond == p_mm_mod->period_bcond.CONST_VOL)
		{
			os << " pbc=xyz " << "        ;  Use periodic boundary conditions in all directions  \n";
			if (p_mm_mod->run_type == p_mm_mod->run_type.MD_RUN)
			{
				if (p_mm_mod->period_bcond == p_mm_mod->period_bcond.CONST_VOL || p_mm_mod->pressure_reg_method == p_mm_mod->pressure_reg_method.NO_CRD_SCALING)
				{
					os << " pcoupl=no " << "        ;  No pressure coupling. This means a fixed box size  \n";
				}
				else if (p_mm_mod->period_bcond == p_mm_mod->period_bcond.CONST_PRES)
				{
					// os << " pcoupl=Berendsen " << "        ;  Exponential relaxation pressure coupling  " << std::endl;
					// pressure regulation
					if (p_mm_mod->pressure_reg_method == p_mm_mod->pressure_reg_method.ISOTROP_CRD_SCALING)
					{
						// os << " pcoupltype=isotropic " << "    ;  Isotropic pressure coupling  " << std::endl;
						os << " ref-p= " << p_mm_mod->ref_pressure << "    ;  The reference pressure for coupling  \n";
						os << " tau-p= " << p_mm_mod->press_relax_time << "    ;  The time constant for pressure coupling  \n";
						if (fabs(p_mm_mod->compressibility - 44.6) > 1.0)
						{
							os << " compressibility= " << p_mm_mod->compressibility * 1.0E-6 << "    ;  The compressibility ( bar^-1)  \n";
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

		os << std::endl << "; Non-bonded interactions \n\n";
 		if (p_mm_model->electr_method == p_mm_model->electr_method.PME_METHOD )
		{
			os << " coulombtype=PME " << "            ;  Fast Smooth Particle-Mesh Ewald (SPME) electrostatics  \n";
		}
		else if (p_mm_model->electr_method == p_mm_model->electr_method.SCREENED_COULOMB )
		{
			os << " coulombtype=Cut-off " << "        ; Plain cut-off with pair list radius rlist and Coulomb cut-off rcoulomb  \n";
			if (fabs(p_mm_model->diel_const - 1.0) > 0.00001)
			{
				os << " epsilon-r= " << p_mm_model->diel_const << "        ; The relative dielectric constant  \n";
			}
		}
		else
		{
			throw std::runtime_error((std::string)("Electrostatics method ") + p_mm_model->electr_method.label() + " is not supported in GROMACS ");
		}

		os << " cutoff-scheme=Verlet " << "        ;  Generate a pair list with buffering \n";
		os << " rlist= " << boost::format("%8.3f") % (p_mm_model->nb_cut_dist * 0.1) << "        ;  Cut-off distance for the short-range neighbor list (nm). \n";

		os << " rcoulomb= " << (p_mm_model->nb_cut_dist * 0.1) << "        ;  distance for the Coulomb cut-off (nm) \n";
		os << " vdwtype=Cut-off " << "        ;  Plain cut-off  " << std::endl;
		os << " rvdw= " << (p_mm_model->nb_cut_dist * 0.1) << "        ;  distance for the LJ cut-off (nm) \n";

		// temperature regulation:
		if (p_mm_mod->run_type == p_mm_mod->run_type.MD_RUN)
		{
			os << "\n; Temperature Control \n\n";

			if (p_mm_mod->temp_control_method == p_mm_mod->temp_control_method.CONST_ENE_MD)
			{
				os << " tcoupl=no " << "        ;  No temperature coupling.\n";
			}
			else if (p_mm_mod->temp_control_method == p_mm_mod->temp_control_method.CONST_TEMP_LANGEVIN || p_mm_mod->temp_control_method == p_mm_mod->temp_control_method.CONST_TEMP_BERENDSEN)
			{
				if (p_mm_mod->temp_control_method == p_mm_mod->temp_control_method.CONST_TEMP_BERENDSEN)
				{
					os << " tcoupl=berendsen " << "        ;  Temperature coupling with a Berendsen thermostat.\n";
				}
				os << " tc-grps= system" << "              ; groups to couple to separate temperature baths  \n";
				os << " ref-t= " << p_mm_mod->ref_temp << "        ;  reference temperature for coupling (K) \n";
				os << " tau-t= " << p_mm_mod->langevin_dump_const << "        ;  time constant for coupling (ps) \n";
			}
			else
			{
				throw std::runtime_error((std::string)("Temperature control method ") + p_mm_mod->temp_control_method.label() + " is not supported in GROMACS ");
			}

			if (fabs(p_mm_mod->init_temp - 0.0) > 0.00001)
			{
				os << " gen-vel=yes " << "   ;  Generate velocities in gmx grompp according to a Maxwell distribution \n";
				os << " gen-temp= " << p_mm_mod->init_temp << "\n";
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

	AtomGroup* p_restr_atoms = p_mm_model->GetRestrAtoms();
	std::set<HaAtom*> restr_atoms_set;
	if (p_restr_atoms)
		for (HaAtom* aptr : *p_restr_atoms)
			restr_atoms_set.insert(aptr);

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
			if (grp_id == "HOH" || grp_id == "WAT") grp_id = "SOL";

			if (!last_group.empty() && grp_id == last_group)
			{
				last_group_count++;
				continue;
			}
			
			if (gmx_mol_defined.count(grp_id) == 0)  // Save description for a new GROMACS molecule 
			{
				AtomIntMap at_idx_map;
				int idx = 0;
				bool has_restr_atoms = false;
				for (HaAtom* aptr : group)
				{
					if (restr_atoms_set.count(aptr) > 0) has_restr_atoms = true;
					at_idx_map[aptr] = idx;
					++idx;
				}

				os << "  \n";
				os << "[ moleculetype ]\n";
				os << ";name            nrexcl \n";
				os << " " << grp_id << "           3 \n";
				os << "  \n";

				SaveAtomsToStream(os, group, at_idx_map);
				SaveBondsToStream(os, group, at_idx_map);
				Save14PairsToStream(os, group, at_idx_map);
				SaveAnglesToStream(os, group, at_idx_map);
				SaveDihedralsToStream(os, group, at_idx_map);

				if (has_restr_atoms)
				{
					std::string fname_restr = grp_id + "_restr.itp";
					PrintLog("Save GROMACS position restraints file %s \n", fname_restr.c_str());

					os << "#include \"" << fname_restr << "\" \n\n";

					std::ofstream os_restr(fname_restr, std::ios::binary);
					if (os_restr.fail())
					{
						PrintLog("Error in %s %d \n", __FILE__,__LINE__);
						PrintLog("Can't create file %s \n", fname_restr.c_str());
						return FALSE;
					}
					os_restr << "; position restraints for Molecule " << grp_id << "\n\n";
					os_restr << "[position_restraints] \n";
					os_restr << ";  i funct       fcx        fcy        fcz  \n";
					int funct_restr = 1;
					double k_restr = 1000.0;
					for (HaAtom* aptr : group)
					{
						if (restr_atoms_set.count(aptr) > 0)
						{
							int idx_at = at_idx_map[aptr] + 1;
							os_restr << boost::format(" %6d %4d  %10.3f %10.3f %10.3f ") % idx_at % funct_restr % k_restr % k_restr % k_restr << "\n";
						}	
					}
					os_restr << " \n";
				}

				gmx_mol_defined.insert(grp_id);
			}

			if (!last_group.empty())
				system_mol_list.push_back((boost::format("%-30s %-6d") % last_group % last_group_count).str());
			last_group = grp_id;
			last_group_count = 1;
		}

		if (!last_group.empty()) 
			system_mol_list.push_back((boost::format("%-30s %-6d") % last_group % last_group_count).str());

		os << "\n" << "; Include water topology \n";
		os << "#include \"" << ffdb_prefix << "/tip3p.itp\" \n";

		os << "\n" << "; Include topology for ions \n";
		os << "#include \"" << ffdb_prefix << "/ions.itp\" \n";

		os << "\n" << "[ system ]\n";
		os << "; Name \n";
		os << pmset->GetName() << "  \n";

		os << std::endl << "[ molecules ] \n";
		os << "; Compound        #mols \n" ;
		for (std::string mol_str : system_mol_list)
			os << mol_str << "  \n";
		os << "  \n";
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

bool MMDriverGromacs::SaveAtomsToStream(std::ostream& os, AtomGroup& group, AtomIntMap& at_idx_map)
{		
	bool has_mut_atoms = p_mm_model->has_mut_atoms_in_group(group);

	os << "[ atoms ]\n";
	if (has_mut_atoms)
	{
		os << ";   nr       type  resnr residue  atom   cgnr     charge       mass              typeB    chargeB     massB \n";
	}
	else
	{
		os << ";   nr  type   resi  res  atom    cgnr     charge      mass       ; qtot   bond_type \n";
	}

	double qtot = 0.0;
	double qtot_mut = 0.0;
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

		std::string ff_s_mut = ff_s;
		double ch_mut = ch;
		double mass_mut = mass;
		
		if (p_mm_model->atom_mut_params.count(aptr) > 0)
		{
			ff_s_mut = p_mm_model->atom_mut_params[aptr]->ff_symbol;
			ch_mut   = p_mm_model->atom_mut_params[aptr]->GetCharge();
			//mass_mut = p_mm_model->atom_mut_params[aptr]  //  IGOR_TMP FIX!
		}
		qtot_mut += ch_mut;

		int idx = at_idx_map[aptr] + 1;
		os << boost::format("%6d %5s %6d %5s %5s %6d %12.6f %12.6f ") % idx % ff_s % resi % res_name % atn % idx % ch % mass;
		if (has_mut_atoms)
		{
			os << boost::format("%5s %12.6f %12.6f ") % ff_s_mut % ch_mut % mass_mut;
		}
		os << boost::format("; qtot %7.3f ") % qtot;
		if (has_mut_atoms)
		{
			os << boost::format("; qtot_mut %7.3f ") % qtot_mut;
		}
		os << "\n";
	}
	os << "  " << "\n";
	return TRUE;
}

bool MMDriverGromacs::SaveBondsToStream(std::ostream& os, AtomGroup& group, AtomIntMap& at_idx_map)
{
	bool has_mut_atoms = p_mm_model->has_mut_atoms_in_group(group);

	os << "[ bonds ]\n";
	if(has_mut_atoms)
	{ 
		os << ";   ai     aj funct   r1            k1    r2            k2\n";
	}
	else
	{
		os << ";   ai     aj funct   r             k   \n";
	}

	set<MMBond>::iterator mbitr = p_mm_model->MBonds.begin();
	for(; mbitr != p_mm_model->MBonds.end(); mbitr++)
	{
		MMBond& bnd = (MMBond&) *mbitr;		
		HaAtom* aptr1 = (HaAtom*) bnd.pt1;
		HaAtom* aptr2 = (HaAtom*) bnd.pt2;

		if( at_idx_map.count(aptr1) == 0 ) continue;
		if( at_idx_map.count(aptr2) == 0 ) continue;

		const int idx1 = at_idx_map[aptr1]+1;  // convert to 1-based index
		const int idx2 = at_idx_map[aptr2]+1;  // convert to 1-based index

		const double r0 = bnd.r0 * 0.1; // Ang -> nm
		const double fc = bnd.fc*4.184*4.184*100; // kcal/mol/Ang^2  ->kJ/nm^2

		std::string at_lbl_1 = aptr1->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_2 = aptr2->GetRef(HaAtom::ATOMREF_STD);

		os << boost::format("%6d  %6d   1 %14.4e%14.4e ") % idx1 % idx2 % r0 % fc;
		if (has_mut_atoms)
		{
			auto mbitr_mut = p_mm_model->MBonds_mut.find(bnd);
			if (mbitr_mut != p_mm_model->MBonds_mut.end())
			{
				const double r0_m = mbitr_mut->r0 * 0.1; // Ang -> nm
				const double fc_m = mbitr_mut->fc * 4.184 * 4.184 * 100; // kcal/mol/Ang^2  ->kJ/nm^2
				os << boost::format(" %14.4e%14.4e ") % r0_m % fc_m;
			}
			else
			{
				os << boost::format(" %14.4e%14.4e ") % r0 % fc;
			}
		}
		os << boost::format("; %s - %s \n") % at_lbl_1, at_lbl_2;
	}

	if (has_mut_atoms)
	{
		for (; mbitr != p_mm_model->MBonds_mut.end(); mbitr++)  // Mutated bonds not found in p_mm_model->MBonds
		{
			MMBond& bnd_m = (MMBond&)*mbitr;
			HaAtom* aptr1 = (HaAtom*)bnd_m.pt1;
			HaAtom* aptr2 = (HaAtom*)bnd_m.pt2;

			if (at_idx_map.count(aptr1) == 0) continue;
			if (at_idx_map.count(aptr2) == 0) continue;

			const int idx1 = at_idx_map[aptr1] + 1;  // convert to 1-based index
			const int idx2 = at_idx_map[aptr2] + 1;  // convert to 1-based index

			const double r0 = bnd_m.r0 * 0.1; // Ang -> nm
			const double fc = bnd_m.fc * 4.184 * 4.184 * 100; // kcal/mol/Ang^2  ->kJ/nm^2

			std::string at_lbl_1 = aptr1->GetRef(HaAtom::ATOMREF_STD);
			std::string at_lbl_2 = aptr2->GetRef(HaAtom::ATOMREF_STD);

			os << boost::format("%6d  %6d   1 %14.4e%14.4e %14.4e%14.4e") % idx1 % idx2 % r0 % fc % r0 % fc;
			os << boost::format("; %s - %s \n") % at_lbl_1, at_lbl_2;
		}
	}

	os << "  \n";
	return true;
}

bool MMDriverGromacs::Save14PairsToStream(std::ostream& os, AtomGroup& group, AtomIntMap& at_idx_map)
{
	bool has_mut_atoms = p_mm_model->has_mut_atoms_in_group(group);
	int im = -1;
	for (auto* p_dih_list : {&(p_mm_model->Dihedrals), &(p_mm_model->Dihedrals_mut)} )
	{
		im++;
		if (im == 1 && !has_mut_atoms) continue;

		os << "[ pairs ]\n";
		os << ";   ai     aj    funct\n";

		for (shared_ptr<MMDihedral> ditr : *p_dih_list )
		{
			HaAtom* aptr1 = ditr->pt1;
			HaAtom* aptr2 = ditr->pt4;

			if (at_idx_map.count(aptr1) == 0) continue;
			if (at_idx_map.count(aptr2) == 0) continue;

			int idx1 = at_idx_map[aptr1] + 1;  // convert to 1-based index
			int idx2 = at_idx_map[aptr2] + 1;  // convert to 1-based index

			std::string at_lbl_1 = aptr1->GetRef(HaAtom::ATOMREF_STD);
			std::string at_lbl_2 = aptr2->GetRef(HaAtom::ATOMREF_STD);

			os << boost::format("%6d  %6d  1 ; %s - %s \n") % idx1 % idx2 % at_lbl_1 % at_lbl_2;
		}
		os << " \n";
	}
	return true;
}

bool MMDriverGromacs::SaveAnglesToStream(std::ostream& os, AtomGroup& group, AtomIntMap& at_idx_map)
{
	bool has_mut_atoms = p_mm_model->has_mut_atoms_in_group(group);

	os << "[ angles ] \n";
	if (has_mut_atoms)
	{
		os << ";   ai     aj     ak    funct   theta1        k1    theta2    k2\n";
	}
	else
	{
		os << ";   ai     aj     ak    funct   theta         k \n";
	}

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

		os << boost::format("%6d %6d %6d  1 %14.4e%14.4e ") % idx1 % idx2 % idx3 % a0 % fc;
		if (has_mut_atoms)
		{
			auto vitr_mut = p_mm_model->ValAngles_mut.find(ang);
			if (vitr_mut != p_mm_model->ValAngles_mut.end())
			{
				const double a0_m = vitr_mut->a0; // degrees
				const double fc_m = vitr_mut->fc * 4.184; // kJ/rad^2  ??
				os << boost::format(" %14.4e%14.4e ") % a0_m % fc_m;
			}
			else
			{
				os << boost::format(" %14.4e%14.4e ") % a0 % fc;
			}
		}
		os << boost::format("; %s - %s - %s \n") % at_lbl_1, at_lbl_2, at_lbl_3;
	}

	if (has_mut_atoms)
	{
		for (; vaitr != p_mm_model->ValAngles_mut.end(); vaitr++)  // Mutated Valence Angles not found in p_mm_model->ValAngles
		{
			MMValAngle& ang_m = (MMValAngle&)(*vaitr);

			double a0 = ang_m.a0;
			double fc = ang_m.fc * 4.184;

			HaAtom* aptr1 = (HaAtom*)ang_m.pt1;
			HaAtom* aptr2 = (HaAtom*)ang_m.pt2;
			HaAtom* aptr3 = (HaAtom*)ang_m.pt3;

			if (at_idx_map.count(aptr1) == 0) continue;
			if (at_idx_map.count(aptr2) == 0) continue;
			if (at_idx_map.count(aptr3) == 0) continue;

			const int idx1 = at_idx_map[aptr1] + 1;  // convert to 1-based index
			const int idx2 = at_idx_map[aptr2] + 1;  // convert to 1-based index
			const int idx3 = at_idx_map[aptr3] + 1;  // convert to 1-based index

			std::string at_lbl_1 = aptr1->GetRef(HaAtom::ATOMREF_STD);
			std::string at_lbl_2 = aptr2->GetRef(HaAtom::ATOMREF_STD);
			std::string at_lbl_3 = aptr3->GetRef(HaAtom::ATOMREF_STD);

			os << boost::format("%6d %6d %6d  1 %14.4e%14.4e %14.4e%14.4e") % idx1 % idx2 % idx3 % a0 % fc % a0 % fc; 
			os << boost::format("; %s - %s - %s \n") % at_lbl_1, at_lbl_2, at_lbl_3;
		}
	}
	os << "  \n";

	return true;
}

bool MMDriverGromacs::SaveDihedralsToStream(std::ostream& os, AtomGroup& group, AtomIntMap& at_idx_map)
{
	bool has_mut_atoms = p_mm_model->has_mut_atoms_in_group(group);

	os << "[ dihedrals ] ; proper dihedrals \n";
	if (has_mut_atoms)
	{
		os << ";    i      j      k      l   func   phase     kd      pn   phase2     kd2      pn2";
	}
	else
	{
		os << ";    i      j      k      l   func   phase     kd      pn \n";
	}

	bool mutated_state = false;
	for (auto* p_dih_list : { &(p_mm_model->Dihedrals),&(p_mm_model->Dihedrals_mut) })
	{
		for (shared_ptr<MMDihedral> ditr : *(p_dih_list))
		{
			MMDihedral& dih = (MMDihedral&)(*ditr);

			std::vector<double> pk = dih.pk;
			std::vector<double> pn = dih.pn;
			std::vector<double> phase = dih.phase;

			HaAtom* aptr1 = (HaAtom*)dih.pt1;
			HaAtom* aptr2 = (HaAtom*)dih.pt2;
			HaAtom* aptr3 = (HaAtom*)dih.pt3;
			HaAtom* aptr4 = (HaAtom*)dih.pt4;

			if (at_idx_map.count(aptr1) == 0) continue;
			if (at_idx_map.count(aptr2) == 0) continue;
			if (at_idx_map.count(aptr3) == 0) continue;
			if (at_idx_map.count(aptr4) == 0) continue;

			int idx1 = at_idx_map[aptr1] + 1;  // convert to 1-based index
			int idx2 = at_idx_map[aptr2] + 1;  // convert to 1-based index
			int idx3 = at_idx_map[aptr3] + 1;  // convert to 1-based index
			int idx4 = at_idx_map[aptr4] + 1;  // convert to 1-based index

			std::string at_lbl_1 = aptr1->GetRef(HaAtom::ATOMREF_STD);
			std::string at_lbl_2 = aptr2->GetRef(HaAtom::ATOMREF_STD);
			std::string at_lbl_3 = aptr3->GetRef(HaAtom::ATOMREF_STD);
			std::string at_lbl_4 = aptr4->GetRef(HaAtom::ATOMREF_STD);

			for (int i = 0; i < pk.size(); i++)
			{
				pk[i] *= 4.184;

				os << boost::format("%6d %6d %6d %6d   1 ") % idx1 % idx2 % idx3 % idx4;
				if (mutated_state)
				{
					double phase1 = 0.0;
					double pk1 = 0.0;
					double pn1 = 1.0;
					os << boost::format(" %14.4e %14.4e %2.0f") % phase1 % pk1 % pn1;
					if (has_mut_atoms)
					{
						os << boost::format(" %14.4e %14.4e %2.0f") % phase[i] % pk[i] % pn[i];
					}
				}
				else
				{
					os << boost::format( " %14.4e %14.4e %2.0f") % phase[i] % pk[i] % pn[i];
					if (has_mut_atoms)
					{
						double phase_mut = 0.0;
						double pk_mut = 0.0;
						double pn_mut = 1.0;
						os << boost::format(" %14.4e %14.4e %2.0f") % phase_mut % pk_mut % pn_mut;
					}				}
				os << boost::format("; %s - %s - %s - %s \n") % at_lbl_1 % at_lbl_2 % at_lbl_3 % at_lbl_4;
			}
		}
		mutated_state = true;
	}
	os << "  " << std::endl;

	os << "[ dihedrals ] ; improper dihedrals \n";
	if (has_mut_atoms)
	{
		os << ";    i      j      k      l   func   phase     kd      pn  phase2     kd2      pn2 \n";
	}
	else
	{
		os << ";    i      j      k      l   func   phase     kd      pn \n";
	}
	
	mutated_state = false;
	for (auto* p_dih_list : { &(p_mm_model->ImprDihedrals),&(p_mm_model->ImprDihedrals_mut) })
	{
		for (shared_ptr<MMDihedral> ditr : (*p_dih_list))
		{
			MMDihedral& dih = (MMDihedral&)(*ditr);

			std::vector<double> pk = dih.pk;
			std::vector<double> pn = dih.pn;
			std::vector<double> phase = dih.phase;

			HaAtom* aptr1 = (HaAtom*)dih.pt1;
			HaAtom* aptr2 = (HaAtom*)dih.pt2;
			HaAtom* aptr3 = (HaAtom*)dih.pt3;
			HaAtom* aptr4 = (HaAtom*)dih.pt4;

			if (at_idx_map.count(aptr1) == 0) continue;
			if (at_idx_map.count(aptr2) == 0) continue;
			if (at_idx_map.count(aptr3) == 0) continue;
			if (at_idx_map.count(aptr4) == 0) continue;

			int idx1 = at_idx_map[aptr1] + 1;  // convert to 1-based index
			int idx2 = at_idx_map[aptr2] + 1;  // convert to 1-based index
			int idx3 = at_idx_map[aptr3] + 1;  // convert to 1-based index
			int idx4 = at_idx_map[aptr4] + 1;  // convert to 1-based index

			std::string at_lbl_1 = aptr1->GetRef(HaAtom::ATOMREF_STD);
			std::string at_lbl_2 = aptr2->GetRef(HaAtom::ATOMREF_STD);
			std::string at_lbl_3 = aptr3->GetRef(HaAtom::ATOMREF_STD);
			std::string at_lbl_4 = aptr4->GetRef(HaAtom::ATOMREF_STD);

			for (int i = 0; i < pk.size(); i++)
			{
				pk[i] *= 4.184;

				os << boost::format("%6d %6d %6d %6d   1 ") % idx1 % idx2 % idx3 % idx4;
				if (mutated_state)
				{
					double phase1 = 0.0;
					double pk1 = 0.0;
					double pn1 = 1.0;
					os << boost::format(" %14.4e %14.4e %2.0f") % phase1 % pk1 % pn1;
					if (has_mut_atoms)
					{
						os << boost::format(" %14.4e %14.4e %2.0f") % phase[i] % pk[i] % pn[i];
					}
				}
				else
				{
					os << boost::format(" %14.4e %14.4e %2.0f") % phase[i] % pk[i] % pn[i];
					if (has_mut_atoms)
					{
						double phase_mut = 0.0;
						double pk_mut = 0.0;
						double pn_mut = 1.0;
						os << boost::format(" %14.4e %14.4e %2.0f") % phase_mut % pk_mut % pn_mut;
					}
				}
				os << boost::format("; %s - %s - %s - %s \n") % at_lbl_1 % at_lbl_2 % at_lbl_3 % at_lbl_4;
			}
		}
		mutated_state = true;
	}
	os << "  " << std::endl;
	return TRUE;
}
