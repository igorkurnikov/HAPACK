/*!  \file mm_driver_arbalest.cpp

    Molecular Mechanics simulations using ARBALEST package  

    \author Igor Kurnikov 
    \date 2025-

*/
#include <mpi.h>
#include <math.h>
#include <tinyxml.h>

#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

namespace fs = std::filesystem;

#include "harlemapp.h"
#include "hamolecule.h"
#include "hamolset.h"
#include "hamolmech.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "mm_driver_arbalest.h"

MMDriverArbalest::MMDriverArbalest(HaMolMechMod* p_mm_mod_new)
{
	p_mm_mod = p_mm_mod_new;
	p_mm_model    = p_mm_mod->p_mm_model;
	pmset = p_mm_mod->GetMolSet();
	
	p_mm_mod->traj_wrt_format = p_mm_mod->traj_wrt_format.TRR;
	
	this->arbalest_exe = "ARBALEST";
	std::string prefix = pmset->GetName();
	this->SetFileNamesWithPrefix(prefix);

	to_save_input_files = TRUE;
}

MMDriverArbalest::~MMDriverArbalest()
{

}

void MMDriverArbalest::SetFileNamesWithPrefix(std::string prefix)
{
	prefix += "_arb";
	if (p_mm_mod->run_type == p_mm_mod->run_type.MIN_RUN)
		prefix += "_min";
	if (p_mm_mod->run_type == p_mm_mod->run_type.MD_RUN)
		prefix += "_md";

	config_fname = prefix + ".xml";
	trj_fname = prefix + ".trr";
	ene_fname = prefix + ".ene";
	run_fname = prefix + ".sh";
}

int MMDriverArbalest::SaveAllInpFiles()
{	
	PrintLog("Save ARBALEST config file %s\n", config_fname);
	SaveConfigFile();
	PrintLog("Save ARBALEST run file %s\n", run_fname.c_str());
	SaveRunFiles();
	to_save_input_files = FALSE;
	return TRUE;
}

bool MMDriverArbalest::SaveConfigFile()
{
	if (p_mm_model->to_init_mm_model) p_mm_mod->InitMolMechModel();

	bool bres = false;
	if (p_mm_mod->run_ti)
	{
		fs::path p(this->config_fname);
		std::string prefix = p.stem().string();
		for (int ilmb = 0; ilmb < p_mm_mod->lambda_ti_v.size(); ilmb++)
		{
			std::string config_fname_lmb = prefix + "_L" + std::to_string(ilmb) + ".xml";
			p_mm_mod->idx_lambda_ti = ilmb;
			p_mm_mod->lambda_ti = p_mm_mod->lambda_ti_v[ilmb];
			PrintLog("Save ARBALEST config file %s \n", config_fname_lmb);
			std::ofstream os(config_fname_lmb, std::ios::binary);
			if (os.fail())
			{
				PrintLog("Error in MMDriverArbalest::SaveConfigFile() \n");
				PrintLog("Can't create file %s \n", config_fname_lmb);
				return false;
			}
			bres = SaveConfigToStream(os);
		}
	}
	else
	{
		PrintLog("Save ARBALEST config file %s \n", config_fname);
		std::ofstream os(this->config_fname, std::ios::binary);
		if (os.fail())
		{
			PrintLog("Error in MMDriverArbalest::SaveConfigFile() \n");
			PrintLog("Can't create file %s \n", this->config_fname);
			return FALSE;
		}
		bres = SaveConfigToStream(os);
	}
	return bres;
}

bool MMDriverArbalest::SaveRunFiles()
{
	if (p_mm_model->to_init_mm_model) p_mm_mod->InitMolMechModel();
	fs::path p(this->config_fname);
	std::string prefix = p.stem().string();

	if (p_mm_mod->run_ti)
	{
		std::string run_all_fname = prefix + "_all.sh";
		std::ofstream os_run_all(run_all_fname, std::ios::binary);

		int job_prefix_size = std::min<size_t>(prefix.size(), 3);
		std::string job_prefix = prefix.substr(0,job_prefix_size);

		for (int ilmb = 0; ilmb < p_mm_mod->lambda_ti_v.size(); ilmb++)
		{
			std::string prefix_lmb = prefix + "_L" + std::to_string(ilmb);
			std::string config_fname_lmb = prefix_lmb + ".xml";
			std::string init_crd_fname_lmb = prefix_lmb + "_init.gro";
			std::string run_fname_lmb = prefix_lmb + ".sh";

			p_mm_mod->idx_lambda_ti = ilmb;
			p_mm_mod->lambda_ti = p_mm_mod->lambda_ti_v[ilmb];

			PrintLog("Save ARBALEST run script %s \n", run_fname_lmb);
			std::ofstream os(run_fname_lmb, std::ios::binary);
			if (os.fail())
			{
				PrintLog("Error in MMDriverArbalest::SaveRunFiles() \n");
				PrintLog("Can't create file %s \n", run_fname_lmb);
				return false;
			}
			os << "#!/bin/bash \n";
			os << "#SBATCH --job-name=" << job_prefix << "_L" << ilmb << "       # job name \n";
			if (this->IsUsingGPU()) os << "#SBATCH --partition=gpu-part    # partition(queue) \n";
			else os << "#SBATCH --partition=cpu-part          # partition(queue) \n";

			os << "#SBATCH -t 5-24:00                    # time limit: (D-HH:MM) \n";
			os << "#SBATCH  --ntasks-per-node=" << this->GetNumCpu() << "         # number of cpu cores \n";
			os << "#SBATCH --gpus-per-node=1             # number of GPU(s) per node \n";
			
			if (ilmb == 0)
			{
				os_run_all << "#!/bin/bash \n";
				os_run_all << "#SBATCH --job-name=" << job_prefix << "       # job name \n";
				if (this->IsUsingGPU()) os_run_all << "#SBATCH --partition=gpu-part    # partition(queue) \n";
				else os_run_all << "#SBATCH --partition=cpu-part          # partition(queue) \n";
				os_run_all << "#SBATCH -t 5-24:00                    # time limit: (D-HH:MM) \n";
				os_run_all << "#SBATCH  --ntasks-per-node=" << this->GetNumCpu() << "         # number of cpu cores \n";
				os_run_all << "#SBATCH --gpus-per-node=1             # number of GPU(s) per node \n";
			}

			os << boost::format("%s --config %s --omp %d ") % arbalest_exe % config_fname_lmb % this->GetNumCpu();
			if (this->IsUsingGPU()) os << " --gpu 1  --gpudeviceid " << this->GetGPUID();
			os << " \n";

			os_run_all << boost::format("%s --config %s --omp %d ") % arbalest_exe % config_fname_lmb % this->GetNumCpu();
			if (this->IsUsingGPU()) os_run_all << " --gpu 1  --gpudeviceid " << this->GetGPUID();
			os_run_all << " \n";
		}
	}
	else
	{
		PrintLog("Save ARBALEST run script %s \n", this->run_fname);
		std::ofstream os(this->run_fname, std::ios::binary);

		if (os.fail())
		{
			PrintLog("Error in MMDriverArbalest::SaveRunFiles() \n");
			PrintLog("Can't create file %s \n", this->run_fname);
			return false;
		}
		os << "#!/bin/sh -f -x \n";
		os << boost::format("%s --config %s --omp %d ") % arbalest_exe % config_fname % this->GetNumCpu();
		if (this->IsUsingGPU()) os << " --gpu 1  --gpudeviceid " << this->GetGPUID();
		os << " \n";
	}
	return true;
}

int MMDriverArbalest::SaveConfigToStream(std::ostream& os)
{
	os << "<Configuration> \n";
	os << "  <ConfigureInteractions> \n";
	os << "    <ForceField Title =\"1018\"> \n";
	os << "      <FFConfig>Input/FFConfig.xml</FFConfig> \n";
	os << "      <Settings> \n";
	os << "        <Param Title=\"RTabulated\">10.0</Param> \n";
	os << "        <Param Title=\"NTabulated\">10000</Param> \n";
	os << "        <Param Title=\"EqTolerance\">0.0001</Param> \n";
	os << "        <Param Title=\"PolarConvergenceTol\">0.0001</Param> \n";
	os << "        <Param Title=\"PolarConvergenceTolPME\">0.0</Param> \n";
	os << "        <Param Title=\"PolarEquilibrationMaxSteps\">30</Param> \n";
	os << "      </Settings> \n";
	os << "    </ForceField> \n";

	os << "    <Volume> \n";
	os << boost::format("      <Lx>%f</Lx> \n") % pmset->per_bc->GetA();
	os << boost::format("      <Ly>%f</Ly> \n") % pmset->per_bc->GetB();
	os << boost::format("      <Lz>%f</Lz> \n") % pmset->per_bc->GetC();
	os << "    </Volume> \n";
	os << "    <Boundary> \n";
	os << "      <BoundX>PERIODIC</BoundX> \n";
	os << "      <BoundY>PERIODIC</BoundY> \n";
	os << "      <BoundZ>PERIODIC</BoundZ> \n";
	os << "    </Boundary> \n";

	os << "    <LongRange> \n";
	os << boost::format("      <RCutOff>%f</RCutOff> \n") % p_mm_model->nb_cut_dist;
	os << boost::format("      <RSwitch>%f</RSwitch> \n") % p_mm_model->nb_cut_dist;
	os << boost::format("      <RCutOffVdW>%f</RCutOffVdW> \n") % p_mm_model->nb_cut_dist;
	os << boost::format("      <RSwitchVdW>%f</RSwitchVdW> \n") % p_mm_model->nb_cut_dist;
	os << boost::format("      <dRCutOff>2.0</dRCutOff> \n");
	os << boost::format("      <NLstRefresh>100</NLstRefresh> \n");

	os << boost::format("      <LRangeCoulomb Type=\"PME\"> \n");
	os << boost::format("        <Settings> \n");
	os << boost::format("          <Param Title=\"GridFragmentationX\">40</Param> \n");
	os << boost::format("          <Param Title=\"GridFragmentationY\">40</Param> \n");
	os << boost::format("          <Param Title=\"GridFragmentationZ\">40</Param> \n");
	os << boost::format("          <Param Title=\"BSplineOrder\">5</Param> \n");
	os << boost::format("          <Param Title=\"Tolerance\">5</Param> \n");
	os << boost::format("          <Param Title=\"DielConstRF\">1000000000000.0</Param> \n");
	os << boost::format("        </Settings> \n");
	os << boost::format("      </LRangeCoulomb> \n");

	os << boost::format("      <LRangeVdW> \n");
	os << boost::format("        <CorrectEnergy>true</CorrectEnergy> \n");
	os << boost::format("        <CorrectPressure>true</CorrectPressure> \n");
	os << boost::format("      </LRangeVdW> \n");
	os << "    </LongRange> \n";

	os << boost::format("  <Constraints> \n");
	os << boost::format("    <Tolerance>0.001</Tolerance> \n");
	os << boost::format("    <ConstrGeom>EquilibriumGeom</ConstrGeom> \n");
	os << boost::format("  </Constraints> \n");
	os << "  </ConfigureInteractions> \n";

	os << "  <Topology> \n";

	std::set<std::string> mol_defined;
	std::string pos_restr_list;

	for (int idx_mol = 0; idx_mol < pmset->HostMolecules.size(); idx_mol++)
	{
		HaMolecule* pmol = pmset->HostMolecules[idx_mol];
		std::string mol_name = pmol->GetName();

		std::set<HaAtom*> saved_atoms;
		AtomIteratorMolecule aitr_m(pmol);
		for (HaAtom* aptr = aitr_m.GetFirstAtom(); aptr; aptr = aitr_m.GetNextAtom())
		{
			if (aptr->IsDummy()) continue; // Do not apply restraints to Dummy atoms
			saved_atoms.insert(aptr);
		}

		std::vector<std::string> pos_restr_desc_and_list = GetPosRestraintsDescAndList(saved_atoms);
		pos_restr_list += pos_restr_desc_and_list[1];

		SaveMolDefToStream(os, pmol, mol_defined, pos_restr_desc_and_list[0]);
	}

	os << "  </Topology> \n";
	os << "  <MolecularSystem ForcePlacement=\"true\"> \n";

	os << "    <SystemState Title=\"SystemStateA\"> \n";
	os << "      <Load> \n";
	os << "        <Molecules> \n";
	for (std::string mol_name : mol_defined)
	{
		os << boost::format("          <Molecule Instance=\"StateA\" Title=\"%s\" />\n") % mol_name;
	}
	os << "        </Molecules> \n";

	std::string system_hin_fname = pmset->GetName() + std::string("_SYS.hin");
	
	AtomSaveOptions opt;
	opt.save_sep_solv_mol = TRUE;
	opt.alchemical_state = AlchemicalState::STATE_A;
	pmset->SaveHINFile(system_hin_fname, opt);
	os << "        <StructureFiles> \n";
	os << boost::format("          <File Type=\"HIN\">%s</File>\n") % system_hin_fname;
	os << "        </StructureFiles> \n";

	if (pos_restr_list.size() > 0)
	{
		os << "        <Restraints> \n";
		os << pos_restr_list;
		os << "        </Restraints> \n";
	}
	os << "      </Load> \n";
	os << "    </SystemState> \n";

	if (pmset->has_mut_atoms())
	{
		os << "    <SystemState Title=\"SystemStateB\" TransitionFrom=\"SystemStateA\"> \n";
		os << "      <Transform> \n";
		os << "        <Molecules> \n";

		for (HaMolecule* pmol : pmset->HostMolecules)
		{
			std::string mol_name = pmol->GetName();
			if (pmol->has_mut_atoms())
			{
				os << boost::format("          <Molecule Title=\"%s\" TargetInstance=\"StateB\" /> \n") % mol_name;
			}
		}
		os << "        </Molecules> \n";
		os << "      </Transform> \n";
		os << "    </SystemState> \n";
	}
	os << "  </MolecularSystem> \n";

	if (pmset->has_mut_atoms())
	{
		if (p_mm_mod->run_ti)
		{
			os << "  <ConfigureTI> \n";
			os << "    <Transition Title=\"ProteinABMutation\" SeqID=\"1\"> \n";
			os << "      <Settings> \n";
			os << "        <Param Title=\"StateI\">SystemStateA</Param> \n";
			os << "        <Param Title=\"StateF\">SystemStateB</Param> \n";
			std::string lambdas_str;
			for (double lmb : p_mm_mod->lambda_ti_v)
				lambdas_str += (boost::format("%d ") % lmb).str();
			os << "        <Param Title = \"LambdaValues\">" << lambdas_str << "</Param> \n";
			os << "        <Param Title=\"InitialConformation\">true</Param> \n";
			os << "        <Param Title=\"ScaleFactorPwr\">2.0</Param> \n";
			os << "        <Param Title=\"SoftCoringRadius\">1.5</Param> \n";
			os << "        <Param Title=\"SoftCoringPwr\">1.0</Param> \n";
			os << "      </Settings> \n";
			os << "    </Transition> \n";
			os << "  </ConfigureTI> \n";
		}
	}

	os << "  <GroupsDefinition> \n";
	os << "    <Group Title=\"Water\" Type=\"UNITED\"> \n";
	os << "      <Set> \n";
	os << "        <LogicalItems> \n";
	os << "          <Item>HOH:</Item> \n";
	os << "        </LogicalItems> \n";
	os << "      </Set> \n";
	os << "    </Group> \n";
	os << "  </GroupsDefinition> \n";

	os << "  <TaskSequence> \n";

	if (p_mm_mod->run_ti)
	{
		os << "    <Task Type=\"SetTIPoint\" SeqID=\"1\"> \n";
		os << "      <Settings> \n";
		os << "        <Param Title=\"Transition\">ProteinABMutation</Param> \n";
		os << "        <Param Title=\"TIPoint\">" << p_mm_mod->idx_lambda_ti << "</Param> \n";
		os << "      </Settings> \n";
		os << "      <ConfigureAlgorithms /> \n";
		os << "    </Task> \n";
	}

	os << "    <Task Type=\"MinEnergy\" SeqID=\"3\"> \n";
	os << "      <Settings> \n";
	os << "        <Param Title=\"Method\">SteepestDescent</Param> \n";
	os << "        <Param Title=\"ConvergenceAccuracy\">0.01</Param> \n";
	os << "        <Param Title=\"NumIterations\">1000</Param> \n";
	os << "      </Settings> \n";
	os << "      <ConfigureAlgorithms /> \n";

	os << "      <ConsoleOutput Do=\"True\"> \n";
	os << "        <Output Screen=\"SYSTEM\" Frequency=\"100\" Mode=\"WATERFALL\" NRows=\"100\"> \n";
	os << "          <Item Title=\"EnrgPot\" Param=\"EnrgPot\" FunctionalGroup=\"SYSTEM\" /> \n";
	os << "          <Item FunctionalGroup=\"SYSTEM\" Param=\"EnrgRestraint\" Title=\"EnrgRestraint\" /> \n";
	os << "        </Output> \n";
	os << "      </ConsoleOutput> \n";

	os << "       <FileOutput OutputFolder=\"./Output\"> \n";
	os << "         <Output DataType=\"GRO\" Frequency=\"50\"> \n";
	os << "           <Settings /> \n";
	os << "           <FunctionalGroups> \n";
	os << "             <Group>SYSTEM</Group> \n";
	os << "           </FunctionalGroups> \n";
	os << "         </Output> \n";
	os << "       </FileOutput> \n";
	os << "    </Task> \n";

	if (p_mm_mod->run_type == p_mm_mod->run_type.MD_RUN)
	{
		os << "    <Task Type=\"MD\" SeqID=\"5\"> \n";
		os << "      <Settings> \n";
		os << "        <Param Title=\"TimeStep\">" << p_mm_mod->md_time_step << "</Param> \n";
		os << "        <Param Title=\"NumIterations\">" << p_mm_mod->num_md_steps << "</Param> \n";
		os << "        <Param Title=\"FastForceStepDivider\">8</Param> \n";
		if (p_mm_mod->run_ti)
		{
			os << "        <Param Title=\"CalcLambdaNeighbors\">1</Param> \n";
		}
		os << "      </Settings> \n";
		os << "      <ConfigureAlgorithms> \n";
		os << "        <Algorithm Type=\"INIT_VELOCITY\"> \n";
		os << "          <Settings> \n";
		os << "            <Param Title=\"RefTemp\">" << p_mm_mod->ref_temp << "</Param> \n";
		os << "            <Param Title=\"Seed\">2025</Param> \n";
		os << "          </Settings> \n";
		os << "          <FunctionalGroups> \n";
		os << "            <Group>SYSTEM</Group> \n";
		os << "          </FunctionalGroups> \n";
		os << "        </Algorithm> \n";

		os << "        <Algorithm Type=\"THERMOSTAT\"> \n";
		os << "          <Settings> \n";
		os << "            <Param Title=\"Type\">Langevin</Param> \n";
		os << "            <Param Title=\"RefTemp\">" << p_mm_mod->ref_temp << "</Param> \n";
		os << "            <Param Title=\"RelaxTime\">0.1</Param> \n";
		// os << "            <Param Title=\"NChain\">6</Param> \n";  // Only for Nose-Hoover
		os << "          </Settings> \n";
		os << "          <FunctionalGroups> \n";
		os << "            <Group>SYSTEM</Group> \n";
		os << "          </FunctionalGroups> \n";
		os << "        </Algorithm> \n";

		if (p_mm_mod->period_bcond == p_mm_mod->period_bcond.CONST_PRES)
		{
			os << "        <Algorithm Type=\"BAROSTAT\"> \n";
			os << "          <Settings> \n";
			os << "            <Param Title=\"Type\">Berendsen</Param> \n";
			os << "            <Param Title=\"RefPressure\">1.0</Param> \n";
			os << "            <Param Title=\"RelaxTime\">0.5</Param> \n";
			os << "            <Param Title=\"Compressibility\">4.5e-5</Param> \n";
			os << "          </Settings> \n";
			os << "          <FunctionalGroups> \n";
			os << "            <Group>SYSTEM</Group> \n";
			os << "          </FunctionalGroups> \n";
			os << "        </Algorithm> \n";
		}
		os << "      </ConfigureAlgorithms> \n";

		os << "      <ConsoleOutput Do=\"True\"> \n";
		os << "        <Output Screen=\"SYSTEM\" Frequency=\"100\" Mode=\"WATERFALL\" NRows=\"100\"> \n";
		os << "          <Item Title=\"Temp\" Param=\"Temp\" FunctionalGroup=\"SYSTEM\" /> \n";
		os << "          <Item Title=\"EnrgPot\" Param=\"EnrgPot\" FunctionalGroup=\"SYSTEM\" /> \n";
		os << "          <Item FunctionalGroup=\"SYSTEM\" Param=\"EnrgRestraint\" Title=\"EnrgRestraint\" /> \n";
		os << "          <Item FunctionalGroup=\"SYSTEM\" Param=\"Pressure\" Title=\"Pressure\" /> \n";
		os << "          <Item FunctionalGroup=\"SYSTEM\" Param=\"Density\" Title=\"Density\" /> \n";
		os << "        </Output> \n";
		os << "      </ConsoleOutput> \n";

		os << "       <FileOutput OutputFolder=\"./Output\"> \n";
		os << "          <Output DataType=\"ENE\" Frequency=\"50\"> \n";
		os << "            <Settings> \n";
		os << "              <Param Title=\"Time\">Time</Param> \n";
		os << "              <Param Title=\"Temp\">Temp</Param> \n";
		os << "              <Param Title=\"EnrgPot\">EnrgPot</Param> \n";
		os << "              <Param Title=\"EnrgKin\">EnrgKin</Param> \n";
		os << "              <Param Title=\"Pressure\">Pressure</Param> \n";
		os << "              <Param Title=\"EnrgRestraint\">EnrgRestraint</Param> \n";
		os << "              <Param Title=\"EnrgCoul\">EnrgCoul</Param> \n";
		os << "              <Param Title=\"EnrgVdW\">EnrgVdW</Param> \n";
		os << "              <Param Title=\"EnrgPol\">EnrgPol</Param> \n";
		os << "              <Param Title=\"EnrgBond\">EnrgBond</Param> \n";
		os << "              <Param Title=\"EnrgAngle\">EnrgAngle</Param> \n";
		os << "              <Param Title=\"EnrgStretch\">EnrgStretch</Param> \n";
		os << "              <Param Title=\"EnrgTorsion\">EnrgTorsion</Param> \n";
		os << "              <Param Title=\"EnrgOOP\">EnrgOOP</Param> \n";
		os << "              <Param Title=\"EnrgShell\">EnrgShell</Param> \n";
		os << "              <Param Title=\"EnrgVdW_LRCor\">EnrgVdW_LRCor</Param> \n";
		os << "              <Param Title=\"dHdL_Coul\">dHdL_Coul</Param> \n";
		os << "              <Param Title=\"dHdL_Coul_SC\">dHdL_Coul_SC</Param> \n";
		os << "              <Param Title=\"dHdL_VdW\">dHdL_VdW</Param> \n";
		os << "              <Param Title=\"dHdL_VdW_SC\">dHdL_VdW_SC</Param> \n";
		os << "              <Param Title=\"dHdL_Pol\">dHdL_Pol</Param> \n";
		os << "              <Param Title=\"dHdL_Bonds\">dHdL_Bonds</Param> \n";
		os << "              <Param Title=\"dHdL_Angles\">dHdL_Angles</Param> \n";
		os << "              <Param Title=\"dHdL_StretchBend\">dHdL_StretchBend</Param> \n";
		os << "              <Param Title=\"dHdL_Tors\">dHdL_Tors</Param> \n";
		os << "              <Param Title=\"dHdL_OutPlane\">dHdL_OutPlane</Param> \n";
		os << "              <Param Title=\"dHdL_LRCor\">dHdL_LRCor</Param> \n";
		os << "              <Param Title=\"dHdL\">dHdL</Param> \n";
		os << "            </Settings> \n";
		os << "            <FunctionalGroups> \n";
		os << "              <Group>SYSTEM</Group> \n";
		os << "            </FunctionalGroups>\n";
		os << "          </Output> \n";

		if (p_mm_mod->run_ti)
		{
			os << "          <Output DataType=\"BAR\" Frequency=\"50\"> \n";
			os << "            <Settings> \n";
			os << "              <Param Title=\"Time\">Time</Param> \n";
			os << "              <Param Title=\"EnrgPot\">EnrgPot</Param> \n";
			os << "              <Param Title=\"EnrgCoul\">EnrgCoul</Param> \n";
			os << "              <Param Title=\"EnrgVdW\">EnrgVdW</Param> \n";
			os << "              <Param Title=\"EnrgPol\">EnrgPol</Param> \n";
			os << "              <Param Title=\"EnrgBond\">EnrgBond</Param> \n";
			os << "              <Param Title=\"EnrgAngle\">EnrgAngle</Param> \n";
			os << "              <Param Title=\"EnrgStretch\">EnrgStretch</Param> \n";
			os << "              <Param Title=\"EnrgTorsion\">EnrgTorsion</Param> \n";
			os << "              <Param Title=\"EnrgOOP\">EnrgOOP</Param> \n";
			os << "              <Param Title=\"EnrgShell\">EnrgShell</Param> \n";
			os << "              <Param Title=\"EnrgVdW_LRCor\">EnrgVdW_LRCor</Param> \n";
			os << "            </Settings> \n";
			os << "            <FunctionalGroups> \n";
			os << "              <Group>SYSTEM</Group> \n";
			os << "            </FunctionalGroups>\n";
			os << "          </Output> \n";
		}

		os << "          <Output DataType=\"XTC\" Frequency=\"5000\"> \n";
		os << "            <Settings> \n";
		os << "              <Param Title=\"RecordGRO\">true</Param> \n";
		os << "              <Param Title=\"RecordVEL\">false</Param> \n";
		os << "              <Param Title=\"RecordFORCE\">false</Param> \n";
		os << "            </Settings> \n";
		os << "            <FunctionalGroups> \n";
		os << "              <Group>SYSTEM</Group> \n";
		os << "            </FunctionalGroups>\n";
		os << "          </Output> \n";
		os << "       </FileOutput> \n";
		os << "    </Task> \n";
	}

	os << "  </TaskSequence> \n";
	os << "</Configuration> \n";

	return true;
}

bool MMDriverArbalest::SaveMolDefToStream(std::ostream& os, HaMolecule* pmol, std::set<std::string>& mol_defined, std::string pos_restr_desc )
{
	std::string mol_name = pmol->GetName();
	std::string mol_a_hin_fname = mol_name + ".hin";

	if( mol_defined.count(mol_name) > 0) return true;

	if (pmol->IsSolvent())
	{
		std::set<std::string> solv_ion_names = pmol->GetResidueNames();
		for (std::string solv_ion_res_name : solv_ion_names)
		{
			if (mol_defined.count(solv_ion_res_name) > 0) continue;
			SaveStdMolDefToStream(os, solv_ion_res_name);
			mol_defined.insert(solv_ion_res_name);
		}
		return true;
	}
	
	AtomSaveOptions opt;
	opt.SetSavedAtoms(*pmol);
	opt.alchemical_state = AlchemicalState::STATE_A;
	pmset->SaveHINFile(mol_a_hin_fname, opt);

	mol_defined.insert(mol_name);
	

	os << boost::format("    <MoleculeDefinition Title=\"%s\"> \n") % mol_name;
	os << boost::format("      <StructureType>SINGLERES</StructureType> \n");
	os << boost::format("      <StructureDefinition> \n");

	os << boost::format("        <Instance Title=\"StateA\"> \n");
	os << boost::format("          <TopologySource>STRUCT</TopologySource> \n");
	os << boost::format("            <TopologyFiles> \n");
	os << boost::format("              <File Type=\"HIN\">%s</File>\n") % mol_a_hin_fname;
	os << boost::format("            </TopologyFiles> \n");
	os << boost::format("          <StructureFiles> \n");
	os << boost::format("            <File Type=\"HIN\">%s</File>\n") % mol_a_hin_fname;
	os << boost::format("          </StructureFiles> \n");
	os << boost::format("          <StructureSettings> \n");
	os << boost::format("            <GenerateChargeGroups>AUTO</GenerateChargeGroups> \n");
	os << boost::format("            <ConstrainBonds>false</ConstrainBonds> \n");
	os << boost::format("            <ConstrainAngles>false</ConstrainAngles> \n");
	os << boost::format("            <ConstrainTorsions>false</ConstrainTorsions> \n");
	os << boost::format("          </StructureSettings> \n");

	if (pos_restr_desc.size() > 0)
	{
		os << "          <RestrainingRules> \n";
		os << pos_restr_desc;
		os << "          </RestrainingRules> \n";
	}

	os << boost::format("        </Instance> \n");

	if (pmol->has_mut_atoms())
	{
		std::string mol_b_hin_fname = mol_name + "_b.hin";
		AtomSaveOptions opt;
		opt.SetSavedAtoms(*pmol);
		opt.alchemical_state = AlchemicalState::STATE_B;
		pmset->SaveHINFile(mol_b_hin_fname, opt);

		os << boost::format("        <Instance Title=\"StateB\"> \n");
		os << boost::format("          <TopologySource>STRUCT</TopologySource> \n");
		os << boost::format("            <TopologyFiles> \n");
		os << boost::format("              <File Type=\"HIN\">%s</File>\n") % mol_b_hin_fname;
		os << boost::format("            </TopologyFiles> \n");
		os << boost::format("          <StructureFiles> \n");
		os << boost::format("            <File Type=\"HIN\">%s</File>\n") % mol_b_hin_fname;
		os << boost::format("          </StructureFiles> \n");
		os << boost::format("          <StructureSettings> \n");
		os << boost::format("            <GenerateChargeGroups>AUTO</GenerateChargeGroups> \n");
		os << boost::format("            <ConstrainBonds>false</ConstrainBonds> \n");
		os << boost::format("            <ConstrainAngles>false</ConstrainAngles> \n");
		os << boost::format("            <ConstrainTorsions>false</ConstrainTorsions> \n");
		os << boost::format("          </StructureSettings> \n");
		os << boost::format("        </Instance> \n");
	}

	os << boost::format("      </StructureDefinition> \n");

	if (pmol->has_mut_atoms())
	{
		os << boost::format("      <TransformationRules> \n");
		os << boost::format("        <MutationRule> \n");
		os << boost::format("          <InstA>StateA</InstA> \n");
		os << boost::format("          <InstB>StateB</InstB> \n");
		os << boost::format("          <MutationMaps> \n");

		ResidueIteratorMolecule ritr(pmol);
		int ir_seq = 0;
		for (HaResidue* pres = ritr.GetFirstRes(); pres; pres = ritr.GetNextRes())
		{
			ir_seq++;
			if (pres->has_mut_atoms())
			{
				std::string res_name_a = pres->GetName();
				std::string res_name_b = pres->p_res_transform->res_name_b;
				std::string fname_map_a_b = std::string("map_") + res_name_a + "_" + res_name_b + ".xml";
				int res_no = ir_seq;

				pres->p_res_transform->SaveMutationMapArbalestFmt(fname_map_a_b);

				os << boost::format("            <Map ResID=\"%d\"> \n") % res_no;
				os << boost::format("              <ResA Title=\"%s\">%s</ResA> \n") % res_name_a % res_name_a;
				os << boost::format("              <ResB Title=\"%s\">%s</ResB> \n") % res_name_b % res_name_b;
				os << boost::format("              <MapFile>%s</MapFile> \n") % fname_map_a_b;
				os << boost::format("            </Map> \n");
			}
		}

		os << boost::format("          </MutationMaps> \n");
		os << boost::format("        </MutationRule> \n");

		os << boost::format("      </TransformationRules> \n");
	}
	os << boost::format("    </MoleculeDefinition> \n");

	return true;
}

bool MMDriverArbalest::SaveStdMolDefToStream(std::ostream& os , std::string mol_name )
{
	std::string mol_a_hin_fname = "Input/HIN/" + mol_name + ".hin";

	os << boost::format("    <MoleculeDefinition Title=\"%s\"> \n") % mol_name;
	os << boost::format("      <StructureType>SINGLERES</StructureType> \n");
	os << boost::format("      <StructureDefinition> \n");

	os << boost::format("        <Instance Title=\"StateA\"> \n");
	os << boost::format("          <TopologySource>STRUCT</TopologySource> \n");
	os << boost::format("            <TopologyFiles> \n");
	os << boost::format("              <File Type=\"HIN\">%s</File>\n") % mol_a_hin_fname;
	os << boost::format("            </TopologyFiles> \n");
	os << boost::format("          <StructureFiles> \n");
	os << boost::format("            <File Type=\"HIN\">%s</File>\n") % mol_a_hin_fname;
	os << boost::format("          </StructureFiles> \n");
	os << boost::format("          <StructureSettings> \n");
	os << boost::format("            <GenerateChargeGroups>AUTO</GenerateChargeGroups> \n");
	os << boost::format("            <ConstrainBonds>false</ConstrainBonds> \n");
	os << boost::format("            <ConstrainAngles>false</ConstrainAngles> \n");
	os << boost::format("            <ConstrainTorsions>false</ConstrainTorsions> \n");
	os << boost::format("          </StructureSettings> \n");
	os << boost::format("        </Instance> \n");
	os << boost::format("      </StructureDefinition> \n");
	os << boost::format("    </MoleculeDefinition> \n");

	return true;
}

std::vector<std::string> MMDriverArbalest::GetPosRestraintsDescAndList(std::set<HaAtom*>& atoms_saved)
{
	std::vector<std::string> restr_desc_and_list(2);
	std::ostringstream oss_desc;
	std::ostringstream oss_list;

	AtomGroup* restr_atoms = p_mm_model->GetRestrAtoms();
	if (restr_atoms && restr_atoms->GetNAtoms() > 0)
	{
		bool bres = SavePosRestraintsStream(oss_desc, oss_list, atoms_saved);
		if (bres)
		{
			restr_desc_and_list[0] = oss_desc.str();
			restr_desc_and_list[1] = oss_list.str();
		}
	}
	return restr_desc_and_list;
}


bool MMDriverArbalest::SavePosRestraintsStream(std::ostream& os_desc, std::ostream& os_list, std::set<HaAtom*>& atoms_saved)
{
	if (os_desc.fail()) return false;
	if (os_list.fail()) return false;

	AtomGroup* restr_atoms = p_mm_model->GetRestrAtoms();
	if (!restr_atoms || restr_atoms->GetNAtoms() == 0)
	{
		PrintLog("No Restrained Atoms defined \n");
		return false;
	}

	int num_restr_atoms = restr_atoms->GetNAtoms();
	if (p_mm_model->restr_ref_coords.size() != num_restr_atoms)
	{
		PrintLog("Error in MMDriverArbalest::SavePosRestraintsStream() \n");
		PrintLog("Dimension of Reference Coordinates array is not equal \n");
		PrintLog("To the size of restrained atoms group \n");
		return false;
	}

	//os_desc << "  <RestrainingRules> \n";
	//os_list << "  <Restraints> \n";

	std::string fc_str = (boost::format("%12.6f") % p_mm_model->GetAtomRestrForceConst()).str();

	AtomIteratorAtomGroup aitr(restr_atoms);
	HaAtom* aptr;
	bool use_chain = false;
	aptr = aitr.GetFirstAtom();
	if (aptr)
	{
		HaChain* pchain_fst = aptr->GetHostChain();
		for (; aptr; aptr = aitr.GetNextAtom())
		{
			HaChain* pchain = aptr->GetHostChain();
			if (pchain != pchain_fst)
			{
				use_chain = true;
				break;
			}
		}
	}

	MolSet* pmset = p_mm_model->GetMolSet();
	ResidueIteratorMolSet ritr(pmset);
	std::map<HaAtom*, int> atom_res_seq_map; // map atoms to residue sequential number 
	int n_res_seq = 0;
	for (HaResidue* rptr : ritr)
	{
		if (!rptr->HasAtomsInSet(atoms_saved)) continue;
		n_res_seq++;
		for (HaAtom* aptr : *rptr)
			atom_res_seq_map[aptr] = n_res_seq;
	}

	int idx_r = -1;
	for (aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		idx_r++;
		Vec3D crd_r = p_mm_model->restr_ref_coords[idx_r];

		if (atoms_saved.count(aptr) == 0) continue;

		HaResidue* pres = aptr->GetHostRes();
		HaChain* pchain = aptr->GetHostChain();
		HaMolecule* pmol = aptr->GetHostMol();

		std::string at_name = aptr->GetName();
		std::string res_name = pres->GetName();
		std::string name_mod = pres->GetNameModifier();
		std::string res_name_save = res_name;

		if (res_name == "HIS")
		{
			res_name_save = "HID";
			if (name_mod == "EPSILON") res_name_save = "HIE";
			if (name_mod == "PROT")    res_name_save = "HIP";
		}

		if (res_name == "CYS")
		{
			if (name_mod == "UNPROT")    res_name_save = "CYX";
		}

		//std::string res_n_str = harlem::ToString( pres->GetSerNo() );
		std::string res_n_str = std::to_string(atom_res_seq_map[aptr]); // use sequence number of the residue - it looks like Arbalest uses sequential residue numbers

		std::string atom_id_arb;
		std::string atom_restr_id;

		if (res_name.empty())
		{
			atom_id_arb = at_name;
			atom_restr_id = "PositionRestraint_" + at_name;
		}
		else
		{
			atom_id_arb = res_name_save + "(" + res_n_str + ")";
			if (use_chain)
			{
				atom_id_arb.push_back('\'');
				atom_id_arb.push_back(pchain->ident);
				atom_id_arb.push_back('\'');
			}
			atom_id_arb += ":" + at_name;
			atom_restr_id = "PositionRestraint_" + res_name_save + "_" + res_n_str;
			if (use_chain) atom_restr_id += std::string(1, pchain->ident);
			atom_restr_id += "_" + at_name;
			// PrintLog("use_chain = %d atom_restr_id = %s \n", use_chain, atom_restr_id.c_str());
		}
		if (!os_desc.fail())
		{
			os_desc << "            <RestraintDefinition IsRelativeCoords=\"false\" Title=\"" << atom_restr_id << "\"> \n";
			os_desc << "              <LinkToAtom_a>" << atom_id_arb << "</LinkToAtom_a> \n";
			os_desc << "                <Anchor_C X=\"" << boost::format("%8.3f") % crd_r.GetX_Ang() << "\"";
			os_desc << " Y=\"" << boost::format("%8.3f") % crd_r.GetY_Ang() << "\"";
			os_desc << " Z=\"" << boost::format("%8.3f") % crd_r.GetZ_Ang() << "\"" << " /> \n";
			os_desc << "                <ForceConstant_Ca>" << boost::trim_copy(fc_str) << "</ForceConstant_Ca> \n";
			os_desc << "                <EquilibriumLength_Ca>0</EquilibriumLength_Ca> \n";
			os_desc << "            </RestraintDefinition> \n";
		}
		if (!os_list.fail())
		{
			os_list << "          <Restraint MolGroup=\"\" Title=\"" << atom_restr_id << "\"  /> \n";
		}
	}
	//os_desc << "  </RestrainingRules> \n";
	//os_list << "  </Restraints> \n";
	return TRUE;
}

