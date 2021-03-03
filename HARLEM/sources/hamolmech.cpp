/*! \file hamolmech.cpp

    Classes to perform Molecular Mechanics calculations 
 
    \author Igor Kurnikov 
    \date 1999-2002
*/

#define HAMOLMECH_CPP

#include <assert.h>
#include <float.h>  
#include <math.h> 

#define HARLEM_MPI 1
#include <mpi.h>
#include <stdexcept>

#pragma warning (disable:4786)

#include <chrono>
#include <thread>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>

#include "hatypes.h"
#include "hamolmech.h"
#include "mm_elements.h"
#include "mm_model.h"

#include <wx/string.h>
#include <wx/filename.h>

#include "haio.h"
#include "tokens.h"
#include "hampi.h"
#include "harlemapp.h"

#include "haintermol.h"
#include "haatom.h"
#include "hamolset.h"
#include "moleditor.h"
#include "hamolview.h"
#include "hamolecule.h"
#include "haresdb.h"

#include "mm_traj_anal.h"
#include "mm_driver_amber.h"
#include "mm_driver_gromacs.h"

#ifndef _MSC_VER
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#endif

harlem::RunOptions HaMolMechMod::run_opt_default;

MMSysInfo::MMSysInfo(HaMolMechMod* p_mm_mod_new)
{
	p_mm_mod = p_mm_mod_new;
	clear();
}

MMSysInfo::~MMSysInfo()
{

}

void MMSysInfo::clear()
{
	nstep     = 0;
	time      = 0.0;   
	temp      = 0.0;   
	temp_solute = 0.0; 
	temp_solv   = 0.0;   
	press     = 0.0; 
	pres_x    = 0.0;
	pres_y    = 0.0;
	pres_z    = 0.0;
	press_scale_solute  = 0.0;
	press_scale_solvent = 0.0;

	tot_energy = 0.0;

	kin_ene  = 0.0;
	pot_ene  = 0.0;  
	bond_ene = 0.0; 
	vang_ene = 0.0;  
	dihed_ene = 0.0;  
	
	vdw_ene    = 0.0;  
	vdw_ene_14 = 0.0;    
	vdw_ene_nb = 0.0;    

	electr_ene = 0.0;  
	electr_ene_14 = 0.0;  
	electr_ene_nb = 0.0;   
	polar_ene = 0.0; 
	polar_dip_iter = 0.0;          
	polar_dip_rms  = 0.0;

	gb_ene    = 0.0;
	hbond_ene = 0.0; 

	constraints_ene = 0.0; 

	epol = 0.0;       
	e3body = 0.0;

	kin_ene_plus_half_dt  = 0.0;   
	kin_ene_minus_half_dt = 0.0;  
	kin_ene_pbs           = 0.0;

	kin_ene_com   = 0.0;
	kin_ene_com_x = 0.0;
	kin_ene_com_y = 0.0;
	kin_ene_com_z = 0.0;

	kin_ene_solute  = 0.0;  
	kin_ene_solvent = 0.0;

	virial_tot = 0.0;
	virial_x   = 0.0;
	virial_y   = 0.0;
	virial_z   = 0.0; 

	volume  = 0.0;
	density = 0.0;  

	dv_dlambda = 0.0;

	av_perm_moment = 0.0;
	av_ind_moment  = 0.0;
	av_tot_moment  = 0.0;
	
	pme_err_est = 0.0;
	rms_ene = 0.0;    
	grad_ene_max = 0.0;   
}


HaMolMechMod::HaMolMechMod(MolSet* new_phost_mset):
HaCompMod(COMP_MOD_MOLMECH,new_phost_mset)
{
	p_mm_model      = new MolMechModel(this);
	p_amber_driver  = new MMDriverAmber(this);
	p_mm_model->p_amber_model->p_amber_driver = p_amber_driver;

	p_gromacs_driver = new MMDriverGromacs(this);
	p_mm_info       = new MMSysInfo(this);
	p_ti_mod        = new TISimMod(this);
	p_md_mod        = new MDSimMod(this); 
	p_min_mod       = new MinEneMod(this); 
	p_traj_anal_mod = new MDTrajAnalMod(this);
	p_traj_io_agent = new MDTrajectoryIOAgent(p_amber_driver);

	period_bcond.SetMolMechMod(this);

	inter_model_comm = MPI_COMM_NULL;
	inter_model_rank = 0;
	single_job_comm  = MPI_COMM_NULL;
	single_job_rank  = 0;

	ctrl_thread = NULL;
	run_thread  = NULL;
	internal_mm_running = false;
	ctrl_thread_running = false;
#if !defined(HA_NOGUI)
	p_mm_dlg    = NULL;
#endif
	p_axx_mm_model = NULL;

	SetStdParams();
}

HaMolMechMod::~HaMolMechMod()
{
	delete p_mm_model;
	delete p_amber_driver;
	delete p_mm_info;
	delete p_ti_mod;
	delete p_md_mod;
	delete p_min_mod;
	delete p_traj_anal_mod;
	delete p_traj_io_agent;

	if( p_axx_mm_model != NULL) delete p_axx_mm_model;
}

int HaMolMechMod::SetStdParams()
{
	MolSet* pmset = GetMolSet();
	if( pmset != NULL ) SetPrefix( pmset->GetName() );

	to_init_simulations = TRUE;

	run_internal_flag = TRUE;
	run_ti    = FALSE;
	update_view_flag = TRUE;
	update_view_interval = 0.5;

	ext_proc_id = 0;

	wrap_coord = TRUE;
	restart_flag = FALSE;

// output params:

	wrt_log_freq = 10;
	wrt_rstrt_freq = 100;
	wrt_coord_freq = 0;
	wrt_vel_freq = 0;
	wrt_ener_freq = 10;
	wrt_constr_freq = 0;

    limit_wrt_atoms = 0;

// Potential function:

	nb_list_update_freq = 25;

//  Minimization parameters:
    	
	SetMaxNumMinimSteps( 500 );
	SetNumSteepDescentSteps ( 10 );
	SetInitMinStep( 0.01 );
	SetGradCnvrgVal( 0.0001 );
	SetZMatMin( false );

// MD parameters

	SetNumMDSteps( 1000 );
	SetRemoveInitRBMotion( FALSE );
	SetRemoveRBMotionFreq( 1000 );
	SetStartTime( 0.0 );
	SetMDTimeStep( 0.001 );

// Temperature regulation:

	ref_temp = 300.0;
	init_temp = 300.0;
	langevin_dump_const = 5.0;
	random_seed = 71277;
	rand_vel_freq = 1000;
	scale_init_vel = 0.0;
	last_solute_atom = 0;
	temp_deviation = 0.0;
	temp_relax_time_solute = 0.5;
	temp_relax_time_solvent = 0.5;
	vel_limit = 20.0; 

// Pressure regulation

	ref_pressure = 1.0;
	compressibility = 44.6;
	press_relax_time = 1.0;

// SHAKE bond length constraints:

	if ( p_mm_model->omit_interactions == p_mm_model->omit_interactions.OMIT_BONDS_H) 
	{
		shake_constr = shake_constr.H_ATOM_SHAKE;
	}
	shake_tol = 0.0005;

// Special Water treatment:

	solute_solvent_image_flag = 1;
	remove_solute_nb_cut_flag = 0;
    fast_water_method = NORMAL_FAST_WATER;

	return True;
}

typedef list<int> LIST_INT;


int HaMolMechMod::Initialize()
{
	InitMolMechModel();
	return TRUE;
}

int HaMolMechMod::InitMolMechModel(const ForceFieldType& ff_type)
{
	int ires = p_mm_model->InitModel(ff_type);
	return ires;
}

MolMechModel* HaMolMechMod::GetMolMechModel()
{
	return p_mm_model;
}

const MolMechModel* HaMolMechMod::GetMolMechModel() const
{
	return p_mm_model;
}

void HaMolMechMod::SetMMRunType( const MMRunType& run_type_new)
{
	run_type = run_type_new;
}

void HaMolMechMod::SetRunInternal( int run_internal_flag_new )
{
	run_internal_flag = run_internal_flag_new;
}

void HaMolMechMod::SetMMExternalProg( const MMExternalProg& ext_mm_prog_new)	
{
	run_internal_flag = FALSE;
	ext_mm_prog = ext_mm_prog_new;
}

MDSimMod* HaMolMechMod::GetMDSimMod()
{
	return p_md_mod;	
}

MinEneMod* HaMolMechMod::GetMinEneMod()
{
	return p_min_mod;
}

TISimMod* HaMolMechMod::GetTISimMod()
{
	return p_ti_mod;
}

MDTrajAnalMod* HaMolMechMod::GetTrajAnalMod()
{
	return p_traj_anal_mod;
}


void run_ctrl_thread(HaMolMechMod* p_mm_mod)
{
	if ( p_mm_mod == NULL) return;
	p_mm_mod->ctrl_thread_running = true;
	p_mm_mod->ControlCalc();
	p_mm_mod->ctrl_thread_running = false;
}

int HaMolMechMod::RunCtrlThread()
{
	if( ctrl_thread_running )
	{
		PrintLog(" Error in HaMolMechMod::RunCtrlThread() \n");
		PrintLog(" MM Control thread is already running \n");
		PrintLog(" Stop calculations before starting a new job \n");
		return FALSE;
	}
	std::thread ctrl_t(run_ctrl_thread, this);
	ctrl_t.detach();

	return TRUE;
}

void run_mm(HaMolMechMod* p_mm_mod)
{
	if (p_mm_mod == NULL) return;
	p_mm_mod->internal_mm_running = true;
	p_mm_mod->RunInternal();
	p_mm_mod->internal_mm_running = false;
}

int HaMolMechMod::Run( const harlem::HashMap* popt_par )
{
	const harlem::RunOptions* popt_c = dynamic_cast<const harlem::RunOptions*>(popt_par);
	std::auto_ptr<harlem::RunOptions> popt_auto( popt_c == NULL ? (harlem::RunOptions*) run_opt_default.clone() : (harlem::RunOptions*) popt_c->clone() );
	harlem::RunOptions* popt = popt_auto.get();

	PrintLog("\n HaMolMechMod::Run() pt 1   Current Dir: %s \n", boost::filesystem::current_path().string().c_str() );

	int ires = TRUE;
	if( to_init_simulations || p_mm_model->to_init_mm_model ) ires = InitMMSimulations();
	if(!ires) return TRUE;

	to_stop_simulations = FALSE;
	
	if( run_internal_flag )
	{
		if( !popt->ToRunSync() )
		{
			std::thread mm_t(run_mm, this);
			mm_t.detach();
		}
		else
		{
			RunInternal();
			return TRUE;
		}
	}
	else
	{
		RunExternal( popt );
		if( popt->ToRunSync() ) return TRUE;
	}
	RunCtrlThread();
	// PrintLog("\n HaMolMechMod::Run() pt end   Current Dir: %s \n", boost::filesystem::current_path().string().c_str());
	return TRUE;
}

int HaMolMechMod::RunMinEne( const harlem::HashMap* popt_par )
{
	const harlem::RunOptions* popt_c = dynamic_cast<const harlem::RunOptions*>(popt_par);
	std::auto_ptr<harlem::RunOptions> popt_auto( popt_c == NULL ? (harlem::RunOptions*) run_opt_default.clone() : (harlem::RunOptions*) popt_c->clone() );
	harlem::RunOptions* popt = popt_auto.get();

	SetMMRunType(MMRunType::MIN_RUN);
	return Run( popt );
}

int HaMolMechMod::RunMD( const harlem::HashMap* popt_par )
{
	const harlem::RunOptions* popt_c = dynamic_cast<const harlem::RunOptions*>(popt_par);
	std::auto_ptr<harlem::RunOptions> popt_auto( popt_c == NULL ? (harlem::RunOptions*) run_opt_default.clone() : (harlem::RunOptions*) popt_c->clone() );
	harlem::RunOptions* popt = popt_auto.get();

	SetMMRunType(MMRunType::MD_RUN);
	return Run( popt );
}

int HaMolMechMod::ControlCalc()
{
	typedef std::chrono::system_clock Time;
	auto update_time = Time::now();

	std::chrono::duration<long, std::milli> update_interval = std::chrono::milliseconds( (long) ( 1000 * update_view_interval) );

	for(;;)
	{
		if( update_time > Time::now() )
		{
		    std::this_thread::sleep_for(std::chrono::milliseconds(200));
			continue;
		}
		if(!run_internal_flag && to_stop_simulations)
		{
			PrintLog("Killing External MM Process: %ld \n", ext_proc_id);  
			HarlemApp::KillProc(ext_proc_id);
			to_stop_simulations = FALSE;
		}
		if(update_view_flag && update_time < Time::now() )
		{
			UpdateMolInfo();
			UpdateMolView();
			update_time = Time::now();
			update_time += update_interval;
		}

		if(run_internal_flag) 
		{
			if( !internal_mm_running )
			{
				UpdateMolInfo();
				UpdateMolView();
				PrintMessage("Internal Molecular Mechanics Execution Process has completed");
				break;
			}
		}
		else
		{
			int active_flag = HarlemApp::CheckProcIsActive(ext_proc_id);
			if(!active_flag)
			{
				ext_proc_id = 0;
				PrintMessage("External Molecular Mechanics Execution Process has completed");
				break;
			}
		}
	}
	return TRUE;
}

int HaMolMechMod::RunExternal( const harlem::HashMap* popt_par )
{
	const harlem::RunOptions* popt_c = dynamic_cast<const harlem::RunOptions*>(popt_par);
	std::auto_ptr<harlem::RunOptions> popt_auto( popt_c == NULL ? (harlem::RunOptions*) run_opt_default.clone() : (harlem::RunOptions*) popt_c->clone() );
	harlem::RunOptions* popt = popt_auto.get();

	bool sync = popt->ToRunSync();

	if( /* ext_mm_prog == ext_mm_prog.PMEMD_9 || ext_mm_prog == ext_mm_prog.SANDER_9 || 
		ext_mm_prog == ext_mm_prog.PMEMD_10 || */ ext_mm_prog == ext_mm_prog.PMEMD_12 || ext_mm_prog == ext_mm_prog.PMEMD_18)
	{
		p_amber_driver->RunAmberProg(sync);
	}
	return TRUE;
}

int HaMolMechMod::RunTI(MolMechModel* p_mm_model_2)
{
	this->run_type = MMRunType::MD_RUN;
	int ires = InitMixedHamSimulations(p_mm_model_2);
	if(!ires) return FALSE;

	RunInternal();
	
	return TRUE;
}

int HaMolMechMod::UpdateMolInfo()
{
	if(!run_internal_flag)
	{
		if( /* ext_mm_prog == ext_mm_prog.PMEMD_9 || ext_mm_prog == ext_mm_prog.SANDER_9 || 
			ext_mm_prog == ext_mm_prog.PMEMD_10 || */ ext_mm_prog == ext_mm_prog.PMEMD_12 || ext_mm_prog == ext_mm_prog.PMEMD_18)
		{
			p_amber_driver->LoadAmberRestartFile(p_amber_driver->amber_rst_file.c_str());
			p_amber_driver->LoadAmberMDInfoFile();
		}
	}	
	return TRUE;
}

int HaMolMechMod::UpdateMolView()
{
	char buf[256];
	MolSet* pmset= GetMolSet();
	pmset->info_str.clear();
	if( run_type == MMRunType::MD_RUN )
	{
		sprintf(buf," nstep = %d  time = %9.3f(ps) T = %9.3f(K) Press= %9.3f (atm)",p_mm_info->nstep, p_mm_info->time,p_mm_info->temp,p_mm_info->press);
		pmset->info_str.push_back(buf);
	}
	else if( run_type == MMRunType::MIN_RUN ) 
	{
		sprintf(buf," Minimization Run \n");
		sprintf(buf," nstep = %d  grad_max = %12.6e  RMS ene = %12.6e",p_mm_info->nstep,p_mm_info->grad_ene_max,p_mm_info->rms_ene);
		pmset->info_str.push_back(buf);
	}

	sprintf(buf," Tot ene = %14.4f  Pot Ene = %14.4f  Kin Ene = %14.4f \n", p_mm_info->tot_energy,p_mm_info->pot_ene,p_mm_info->kin_ene);
	pmset->info_str.push_back(buf);
	sprintf(buf," Bond ene = %14.4f  Vang Ene = %14.4f  Dihedral Ene = %14.4f \n", p_mm_info->bond_ene,p_mm_info->vang_ene,p_mm_info->dihed_ene);
	pmset->info_str.push_back(buf);
	sprintf(buf," VdW Ene = %14.4f  El Ene = %16.4f  VdW 1-4 = %12.4f  El 1-4 = %12.4f \n", 
		p_mm_info->vdw_ene,p_mm_info->electr_ene,p_mm_info->vdw_ene_14,p_mm_info->electr_ene_14);
	pmset->info_str.push_back(buf);

	pmset->AnnounceGeomChange();
	return TRUE;
}

int HaMolMechMod::InitMMSimulations()
{
	int ires;
	if(run_ti)
	{
		ires = InitMixedHamSimulations(p_axx_mm_model);
	}
	else
	{
		ires = InitSingleHamMMSimulations();
	}
	return ires;
}


int HaMolMechMod::InitSingleHamMMSimulations()
{
	int ires;
	try 
	{
		if(pApp->mpi_driver->nprocs > 1) 
		{
		    CallMMFunctionOnSlaves(MM_SET_MPI_COMM_ALL_PROCS);
			p_amber_driver->SetMPICommAllProcs();
		}
		p_amber_driver->InitParallelDatMod(); 
		
//		PrintLog(" HaMolMechMod::InitSingleHamMMSimulations() pt 1 \n");
		ires = p_mm_model->UpdateModel();
		if(!ires) throw std::runtime_error("Failed to update the MM Model ");
		
		p_amber_driver->InitCtrlParams();

		if(!run_internal_flag) return TRUE;

		p_amber_driver->SaveModelToFortran();
		p_amber_driver->SetAtomCrdToInternalArrays();
		if( p_amber_driver->numtasks > 1 ) CallMMFunctionOnSlaves(MM_INIT_SIMULATIONS_STEP_2);			
		p_amber_driver->InitSimulationsStep2();

		run_ti = FALSE;
		to_init_simulations = FALSE;
	}
	catch( std::exception& e)
	{
		PrintLog("Exception in HaMolMechMod::InitSingleHamMMSimulations() %s \n", e.what());
		return FALSE;
	}
	catch(...) 
	{
		PrintLog("Unidentified exception occured in HaMolMechMod::InitSingleHamMMSimulations() \n");
		return FALSE;
	}
	return TRUE;
}

int HaMolMechMod::InitMixedHamSimulations(MolMechModel* p_mm_model_2)
{
	int ires = TRUE;
	if( p_mm_model_2 == NULL ) return FALSE;
	ires = p_mm_model->UpdateModel();
	if(!ires) return FALSE;
	ires = p_mm_model_2->UpdateModel();
	if(!ires) return FALSE;

	if(!CheckModelsForTI(p_mm_model,p_mm_model_2)) return FALSE;
	this->p_axx_mm_model = p_mm_model_2;

	CallMMFunctionOnSlaves(MM_MOD_INIT_MIXED_HAMILTONIAN);
	ires = InitMixedHamSimulations_node(p_mm_model_2);
	return ires;
}

int HaMolMechMod::InitMixedHamSimulations_node(MolMechModel* p_mm_model_2)
{
	PrintLog(" HaMolMechMod::InitMixedHamSimulations_node() pt 1 \n");
	try 
	{
		SetMPICommSplit2();
		PrintLog(" HaMolMechMod::InitMixedHamSimulations_node() pt 2 \n");
		p_amber_driver->InitParallelDatMod();
		PrintLog(" HaMolMechMod::InitMixedHamSimulations_node() pt 3 \n");

		if( pApp->mpi_driver->myrank == 0)
		{
			p_amber_driver->InitCtrlParams();	
			p_amber_driver->SetAtomCrdToInternalArrays(); 
		}
		PrintLog(" HaMolMechMod::InitMixedHamSimulations_node() pt 4 \n");
		if( p_amber_driver->master )
		{
			PrintLog(" HaMolMechMod::InitMixedHamSimulations_node() pt 4_1 \n");
			BcastCtrlParams(inter_model_comm);
			PrintLog(" HaMolMechMod::InitMixedHamSimulations_node() pt 4_2 \n");
			if( inter_model_rank == 0 ) // MODEL 1
			{	
				PrintLog(" HaMolMechMod::InitMixedHamSimulations_node() pt 4_3 \n");
				p_mm_model_2->Bcast(inter_model_comm); // Broadcast Model 2 to processors with inter_model_rank != 0
				PrintLog(" HaMolMechMod::InitMixedHamSimulations_node() pt 4_4 \n");
				p_mm_model->p_amber_model->BcastAtmMass(inter_model_comm);  // Broadcast atom masses from  Model 1 (same masses should be used in TI)
			}
			else 
			{
				PrintLog(" HaMolMechMod::InitMixedHamSimulations_node() pt 4_5 \n");
				p_mm_model->Bcast(inter_model_comm);
				PrintLog(" HaMolMechMod::InitMixedHamSimulations_node() pt 4_6 \n");
				p_mm_model->p_amber_model->BcastAtmMass(inter_model_comm);
			}
			PrintLog(" HaMolMechMod::InitMixedHamSimulations_node() pt 4_7 \n");
			p_amber_driver->SaveModelToFortran();
			PrintLog(" HaMolMechMod::InitMixedHamSimulations_node() pt 4_8 \n");
			p_amber_driver->ResizeCrdVelFrcArrays(p_mm_model->p_amber_model->natom);
			p_amber_driver->BcastCrd(inter_model_comm);
			p_amber_driver->BcastVel(inter_model_comm);
			p_amber_driver->BcastPBox(inter_model_comm);
		}
		PrintLog(" HaMolMechMod::InitMixedHamSimulations_node() pt 5 \n");
		run_ti = TRUE;
		p_axx_mm_model = p_mm_model_2;
		p_amber_driver->InitSimulationsStep2();
		to_init_simulations = FALSE;
		to_stop_simulations = FALSE;
	}
	catch( std::exception& e)
	{
		PrintLog("Exception in HaMolMechMod::InitMixedHamSimulations_node %s \n", e.what());
		return FALSE;
	}
	catch(...) 
	{
		PrintLog("Unidentified exception occured in HaMolMechMod::InitMixedHamSimulations_node \n");
		return FALSE;
	}
	PrintLog(" HaMolMechMod::InitMixedHamSimulations_node() pt end \n");
	return TRUE;
}

void HaMolMechMod::BcastCtrlParams(MPI_Comm& comm)
{
	int ires;
#if defined(HARLEM_MPI)
	
	ires = run_type.Bcast(comm);
	ires = init_read_coord.Bcast(comm);
	ires = MPI_Bcast(&restart_flag,1,MPI_INT,0,comm);
	write_coord_format.Bcast(comm);

	ires = MPI_Bcast(&p_amber_driver->ntave,1,MPI_INT,0,comm);

	ires = MPI_Bcast(&wrt_rstrt_freq,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&wrap_coord,1,MPI_INT,0,comm);

	ires = MPI_Bcast(&wrt_coord_freq,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&wrt_constr_freq,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&wrt_vel_freq,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&wrt_ener_freq,1,MPI_INT,0,comm);

	ires = traj_wrt_format.Bcast(comm);

	ires = MPI_Bcast(&limit_wrt_atoms,1,MPI_INT,0,comm);
	ires = period_bcond.Bcast(comm);
	ires = MPI_Bcast(&p_amber_driver->nsnb,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&nb_list_update_freq,1,MPI_INT,0,comm);

	ires = MPI_Bcast(&p_amber_driver->nrespai,1,MPI_INT,0,comm);

	ires = MPI_Bcast(&max_num_minim_steps,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&num_steep_descent_steps,1,MPI_INT,0,comm);
	ires = min_type.Bcast(comm);
	ires = MPI_Bcast(&init_min_step,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&grad_cnvrg_val,1,MPI_DOUBLE,0,comm);

	ires = MPI_Bcast(&remove_rb_motion_freq,1,MPI_INT,0,comm);

	ires = MPI_Bcast(&start_time,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&md_time_step,1,MPI_DOUBLE,0,comm);

	ires = MPI_Bcast(&p_amber_driver->nrespa,1,MPI_INT,0,comm);

	ires = temp_control_method.Bcast(comm);

	ires = MPI_Bcast(&ref_temp,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&init_temp,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&random_seed,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&temp_relax_time_solute,1,MPI_DOUBLE,0,comm);

	ires = MPI_Bcast(&langevin_dump_const,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&rand_vel_freq,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&vel_limit,1,MPI_DOUBLE,0,comm);

	ires = MPI_Bcast(&p_amber_driver->ntp,1,MPI_INT,0,comm);
	ires = pressure_reg_method.Bcast(comm);

	ires = MPI_Bcast(&ref_pressure,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&press_relax_time,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&compressibility,1,MPI_DOUBLE,0,comm);

	ires = MPI_Bcast(&p_amber_driver->ntc,1,MPI_INT,0,comm);
	ires = shake_constr.Bcast(comm);

	ires = MPI_Bcast(&shake_tol,1,MPI_DOUBLE,0,comm);

	ires = MPI_Bcast(&p_amber_driver->jfastw,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&fast_water_method,1,MPI_INT,0,comm);

	ires = MPI_Bcast(&p_amber_driver->dbg_atom_redistribution,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&p_amber_driver->loadbal_verbose,1,MPI_INT,0,comm);
#endif
}

int HaMolMechMod::SaveXMLToStream(std::ostream& os, const harlem::SaveOptions* popt ) const
{
	os << "<mm_mod mset=\"" << GetMolSet()->GetName() << "\">" << std::endl;
	const MolMechModel* p_mm_model = GetMolMechModel();
	if( p_mm_model != NULL ) 
	{
		p_mm_model->SaveXMLToStream( os, popt );
	}
	os << "</mm_mod>" << std::endl;
	return TRUE;
}

int HaMolMechMod::RunInternal()
{
	if(pApp->mpi_driver->nprocs > 1 ) HaMolMechMod::CallMMFunctionOnSlaves(MM_DRIVER_AMBER_RUN_INTERNAL);
	RunInternal_node();

	return TRUE;
}

void HaMolMechMod::RunInternal_node()
{
	PrintLog(" HaMolMechMod::RunInternal_node() pt 1 \n");
	MolSet* pmset = this->GetMolSet();
	p_amber_driver->p_tm->InitTimers();

	if(pApp->mpi_driver->nprocs > 1 ) run_type.Bcast(single_job_comm);

	try
	{
		if ( run_type == MMRunType::MD_RUN )      // Do molecular dynamics:
		{
			p_amber_driver->RunMD(); 
		}
		else if( run_type == MMRunType::MIN_RUN ) // Do energy minimization:
		{
			if(p_amber_driver->master)
			{
				p_amber_driver->RunMinMaster();
			}
			else
			{
				p_amber_driver->RunMinSlave();
			}
			if(p_amber_driver->master) pmset->AnnounceGeomChange();
		}
		else if( run_type == MMRunType::ENER_RUN )
		{
			p_amber_driver->CalcCurrEne();
		}
	}
	catch( std::exception& e)
	{
		PrintLog("Exception in HaMolMechMod::RunInternal_node() %s \n", e.what());
	}   
	catch(...) 
	{
		PrintLog("Unidentified exception occured in HaMolMechMod::RunInternal_node() \n");
	}
	p_amber_driver->p_tm->EndRunTimers();
	if(run_type == MMRunType::MD_RUN || run_type == MMRunType::MIN_RUN) p_amber_driver->p_tm->PrintTimings();
	PrintLog(" HaMolMechMod::RunInternal_node() pt end \n");
}

int HaMolMechMod::CalcEnergy()
{
	int ires = TRUE;
	run_type = MMRunType::ENER_RUN;
	if( to_init_simulations || p_mm_model->to_init_mm_model) ires = InitMMSimulations();
	if(!ires) return FALSE;
	
//	p_amber_driver->SetAtomCrdToInternalArrays();

	RunInternal();

	return TRUE;
}

int HaMolMechMod::CheckModelsForTI(MolMechModel* p_mm_model_1, MolMechModel* p_mm_model_2)
{
	if(pApp->mpi_driver->nprocs < 2) 
	{
		PrintLog("Error in HaMolMechMod::CheckModelsForTI() \n");
		PrintLog("Number of processes is less than 2 \n");
		PrintLog("TI is available only in MPI runs \n");
		return FALSE;
	}
	if(pApp->mpi_driver->nprocs % 2 != 0) 
	{
		PrintLog("Error in HaMolMechMod::CheckModelsForTI() \n"); 
		PrintLog("Number of processes should be a multiple of 2 \n");
		return FALSE;
	}
	if(p_mm_model_1->Atoms.size() == 0 )
	{
		PrintLog("Error in HaMolMechMod::CheckModelsForTI() \n");
		PrintLog("No atoms in model 1 \n");
		return FALSE;
	}
	if(p_mm_model_2->Atoms.size() == 0 )
	{
		PrintLog("Error in HaMolMechMod::CheckModelsForTI() \n"); 
		PrintLog("No atoms in model 2 \n");
		return FALSE;
	}
	if(p_mm_model_1->Atoms.size() !=  p_mm_model_2->Atoms.size())
	{
		PrintLog("Error in HaMolMechMod::CheckModelsForTI() \n");
		PrintLog("Number of atoms in models 1 and 2 is different: %d and %d \n",p_mm_model_1->Atoms.size(),p_mm_model_2->Atoms.size());
		return FALSE;
	}
	return TRUE;
}

void HaMolMechMod::CallMMFunctionOnSlaves(int id)
{
//	PrintLog(" HaMolMechMod::CallMMFunctionOnSlaves() pt 1 id= %d \n", id);
	std::string msg = HaMPI::BuildXMLwxCmdEventBasic(HA_MOL_MECH_EVENT,id);
	pApp->mpi_driver->SendXmlMsgAllProc(msg.c_str());
}

int HaMolMechMod::SetMPICommSplit2()
{
	int ires, itmp;
#if defined(HARLEM_MPI)
	if(p_amber_driver->driver_mpi_comm  != MPI_COMM_NULL) MPI_Comm_free(&(p_amber_driver->driver_mpi_comm));
	if(p_amber_driver->driver_mpi_group != MPI_GROUP_NULL) MPI_Group_free(&(p_amber_driver->driver_mpi_group));
	
	int color = pApp->mpi_driver->myrank % 2;
	int key   = pApp->mpi_driver->myrank;

	ires = MPI_Comm_split(MPI_COMM_WORLD,color,key,&(p_amber_driver->driver_mpi_comm) );
	ires = MPI_Comm_group(p_amber_driver->driver_mpi_comm,&(p_amber_driver->driver_mpi_group));
	ires = MPI_Comm_size(p_amber_driver->driver_mpi_comm,&itmp); p_amber_driver->numtasks = itmp;
	ires = MPI_Comm_rank(p_amber_driver->driver_mpi_comm,&itmp); p_amber_driver->mytaskid = itmp;

	if(inter_model_comm != MPI_COMM_NULL) MPI_Comm_free(&inter_model_comm); 

	color = pApp->mpi_driver->myrank/2;
	key   = pApp->mpi_driver->myrank;

	ires = MPI_Comm_split(MPI_COMM_WORLD,color,key,&inter_model_comm);
	ires = MPI_Comm_rank(inter_model_comm,&inter_model_rank);
	
	if(single_job_comm != MPI_COMM_NULL) MPI_Comm_free(&single_job_comm);
	ires = MPI_Comm_dup(MPI_COMM_WORLD,&single_job_comm);
	ires = MPI_Comm_rank(single_job_comm,&single_job_rank);

//	PrintLog("HaMolMechMod::SetMPICommSplit2(): \n");
//	PrintLog(" pApp->mpi_driver->myrank = %d \n",pApp->mpi_driver->myrank);
//	PrintLog(" inter_model_rank = %d \n",inter_model_rank);
//	PrintLog(" p_amber_driver->mytaskid = %d \n",p_amber_driver->mytaskid);
#endif

	return TRUE;
}


int HaMolMechMod::StopCalc()
{
	to_stop_simulations = TRUE;
	return True;
}

double HaMolMechMod::GetEne() const
{
	return p_mm_info->pot_ene;
}

double HaMolMechMod::GetTotEne() const 
{ 
	return p_mm_info->tot_energy; 
}

double HaMolMechMod::GetPotEne() const 
{ 
	return p_mm_info->pot_ene; 
}
	
double HaMolMechMod::GetConstrEne() const 
{ 
	return p_mm_info->constraints_ene; 
}
	
double HaMolMechMod::GetUnConstrEne() const 
{ 
	return (p_mm_info->pot_ene - p_mm_info->constraints_ene); 
}

int HaMolMechMod::LoadAmberRestartFile(const std::string& rst_file_name)
{
	return p_amber_driver->LoadAmberRestartFile( rst_file_name.c_str() );
}

int HaMolMechMod::CalcEnergySimple()
{
	if(p_mm_model->to_init_mm_model || to_init_simulations ) InitMMSimulations();
	if(p_mm_model->build_nb_contact_list_flag) p_mm_model->BuildNonBondContactList();

	int np = p_mm_model->Atoms.size();
	p_mm_info->tot_energy = 0.0;
	if(np == 0) 
	{	
		return True;
	}
	
	if(np != p_mm_model->nonbond_contact_list.size())
	{
		ErrorInMod("HaMolMechMod::CalcEnergySimple()",
		" The number of Mat Atoms is not equal to the size of nonbonded contact list");
		return FALSE;
	}
	p_mm_info->clear();

	HaAtom* pt1;
	HaAtom* pt2;

	double vdw_at_ene;
	double el_at_ene;

	int nbonds = p_mm_model->MBonds.size();
	int i;
	double val;
	double delt;

	set<MMBond, less<MMBond> >::iterator mbitr = p_mm_model->MBonds.begin();

	for(; mbitr != p_mm_model->MBonds.end(); mbitr++)
	{
		MMBond& bond = (MMBond&) *mbitr;
		val = Vec3D::CalcDistance(bond.pt1,bond.pt2,ANGSTROM_U);
		delt = val - bond.r0;
		p_mm_info->bond_ene += bond.fc * delt * delt; 
	}

	set<MMValAngle, less<MMValAngle> >::iterator vaitr = p_mm_model->ValAngles.begin();

	for(; vaitr != p_mm_model->ValAngles.end(); vaitr++)
	{
		MMValAngle& ang = (MMValAngle&)(*vaitr);
		val = Vec3D::CalcAngle(ang.pt1,ang.pt2,ang.pt3);
		delt = val - ang.a0*DEG_TO_RAD;
		p_mm_info->vang_ene += ang.fc * delt * delt; 
	}

	for( i = 0; i < np; i++)
	{
		int nl_size = p_mm_model->nonbond_contact_list[i].size();
		set<HaAtom*, less<HaAtom*> >::iterator mitr;
		for( mitr = p_mm_model->nonbond_contact_list[i].begin();
		     mitr != p_mm_model->nonbond_contact_list[i].end(); mitr++)
		{
			pt1 = p_mm_model->Atoms[i];
			pt2 = *mitr;

			p_mm_model->CalcNonBondPt(pt1, pt2, vdw_at_ene, el_at_ene);

			p_mm_info->vdw_ene_nb += vdw_at_ene;
			p_mm_info->electr_ene_nb += el_at_ene;
		}
	}

	vector<MMDihedral>::iterator ditr;

	for( ditr = p_mm_model->Dihedrals.begin(); ditr != p_mm_model->Dihedrals.end(); ditr++)
	{
		pt1 = (*ditr).pt1;
		pt2 = (*ditr).pt4;
			
		p_mm_model->CalcNonBondPt(pt1, pt2, vdw_at_ene, el_at_ene);

		p_mm_info->vdw_ene_14 += vdw_at_ene;
		p_mm_info->electr_ene_14 += el_at_ene;
	}

	if(p_mm_model->GetNumHarmConstr() > 0)
	{
		double lc_ene;
		vector<AtomContact>::iterator citr;
		for(citr = p_mm_model->DistConstraints.begin(); citr != p_mm_model->DistConstraints.end(); citr++)
		{
			AtomContact& cnt = (*citr);
			if( cnt.cnt_type != AtomContactType::HARMONIC_CNT ) continue;

			double dist = Vec3D::CalcDistance(cnt.pt1,cnt.pt2,ANGSTROM_U);
			double fc   = cnt.GetHarmForceConst();
			double r0   = cnt.GetRMin();
			lc_ene = 0.5 *(dist - r0)*(dist - r0);
			p_mm_info->constraints_ene += lc_ene;
		}
	}

	p_mm_info->electr_ene_14 = p_mm_info->electr_ene_14/p_mm_model->scale_14_electr;
	p_mm_info->vdw_ene_14    = p_mm_info->vdw_ene_14/p_mm_model->scale_14_vdw;

	p_mm_info->vdw_ene    = p_mm_info->vdw_ene_nb    + p_mm_info->vdw_ene_14;
	p_mm_info->electr_ene = p_mm_info->electr_ene_nb + p_mm_info->electr_ene_14; 

	p_mm_info->tot_energy =  p_mm_info->bond_ene + p_mm_info->vang_ene + 
		                     p_mm_info->electr_ene + p_mm_info->vdw_ene + p_mm_info->constraints_ene;

	PrintLog(" Bond ene = %12.6f,  Val Angle Ene =  %12.6f (kcal/mol) \n", p_mm_info->bond_ene,p_mm_info->vang_ene);

	PrintLog(" Current VdW energy = %f kcal/mol = %f (1-4) + %f (non 1-4)  \n", 
		       p_mm_info->vdw_ene, p_mm_info->vdw_ene_14, p_mm_info->vdw_ene_nb );   

	PrintLog(" Current Electrostatic energy = %f kcal/mol %f (1-4) + %f (non 1-4)  \n",
		       p_mm_info->electr_ene, p_mm_info->electr_ene_14, p_mm_info->electr_ene_nb);
	
	if(!p_mm_model->DistConstraints.empty())
	{
		PrintLog(" Harmonic constraints energy = %12.6f kcal/mol \n", p_mm_info->constraints_ene);
	}

	PrintLog(" Current Total energy %f kcal/mol \n", p_mm_info->tot_energy );

	return True;
}

void HaMolMechMod::PrintEneStr(MMSysInfo& info,std::string& str_out)
{
	wxString str;

	str += wxString::Format(" NSTEP = %9d TIME(PS) =  %10.3f TEMP(K) = %7.2f PRESS = %7.1f\n",    
		info.nstep,info.time, info.temp, info.press);
	str += wxString::Format(" Etot   = %14.4f  EKtot   = %14.4f EPtot     = %14.4f\n",
		info.tot_energy,info.kin_ene,info.pot_ene);
	str += wxString::Format(" BOND   = %14.4f  ANGLE   = %14.4f DIHED     = %14.4f\n",
		info.bond_ene,info.vang_ene,info.dihed_ene);
	str += wxString::Format(" 1-4 NB = %14.4f  1-4 EEL = %14.4f VDWAALS   = %14.4f\n",
		info.vdw_ene_14,info.electr_ene_14,info.vdw_ene);

	if (p_mm_model->p_amber_model->using_pme_potential) 
	{
		str += wxString::Format(" EELEC = %14.4f  EHBOND = %14.4f RESTRAINT = %14.4f\n",
			info.electr_ene,info.hbond_ene,info.constraints_ene);
	}  
	else
	{
		str += wxString::Format(" EELEC = %14.4f  EGB   =  %14.4f RESTRAINT = %14.4f\n",
			info.electr_ene,info.gb_ene,info.constraints_ene);
	}
	if( fabs(info.constraints_ene) > 0)
	{
		str += wxString::Format(" EAMBER (non-restraint) = %14.4f\n",
			(info.pot_ene - info.constraints_ene) );
	}
	if( fabs(info.dv_dlambda) > 0)
	{
		str += wxString::Format(" DV/DL  = %14.4f\n", info.dv_dlambda );
	}
	if( info.volume != 0.0)
	{
		str += wxString::Format(" EKCMT  = %14.4f  VIRIAL  = %14.4f  VOLUME =   %14.4f\n",
			info.kin_ene_com,info.virial_tot,info.volume);
	}
	if( info.kin_ene_solvent != 0.0 )
	{
//		str += wxString::Format(" EK_solv = %14.4f \n", info.kin_ene_solvent );
	}
	if( info.epol != 0.0 || info.e3body != 0.0 )
	{
		str += wxString::Format(" EPOLZ  = %14.4f  E3BODY  = %14.4f\n",
			info.epol,info.e3body);
	}

	if(info.polar_ene != 0.0)
	{
		str += wxString::Format(" E_POLAR = %14.4f \n ", info.polar_ene );
	}
	if(info.density != 0.0)
	{
		str += wxString::Format("                                              Density    = %14.4f\n", info.density );
	}
	if(p_mm_model->p_amber_model->using_pme_potential)
	{
		str += wxString::Format(" Ewald error estimate: %12.4e\n", info.pme_err_est ); 
	}
	str += "----------------------------------------------------------------------------\n";
	//str += "Detailed ENE: \n";
	//str += wxString::Format(" NSTEP = %9d TIME(PS) =  %10.3f TEMP(K) = %7.2f PRESS = %16.9f\n",    
	//	info.nstep,info.time, info.temp, info.press);
	//str += wxString::Format(" Etot   = %16.9f  EKtot   = %16.9f EPtot     = %16.9f\n",
	//	info.tot_energy,info.kin_ene,info.pot_ene);
	//str += wxString::Format(" BOND   = %16.9f  ANGLE   = %16.9f DIHED     = %16.9f\n",
	//	info.bond_ene,info.vang_ene,info.dihed_ene);
	//str += wxString::Format(" 1-4 NB = %16.9f  1-4 EEL = %16.9f VDWAALS   = %16.9f\n",
	//	info.vdw_ene_14,info.electr_ene_14,info.vdw_ene);
	//if (p_mm_model->p_amber_model->using_pme_potential) 
	//{
	//	str += wxString::Format(" EELEC = %16.9f  EHBOND = %16.9f RESTRAINT = %16.9f\n",
	//		info.electr_ene,info.hbond_ene,info.constraints_ene);
	//}  
	//else
	//{
	//	str += wxString::Format(" EELEC = %16.9f  EGB   =  %16.9f RESTRAINT = %16.9f\n",
	//		info.electr_ene,info.gb_ene,info.constraints_ene);
	//}
	//str += "----------------------------------------------------------------------------\n";
	str_out = str;
}

void HaMolMechMod::PrintEneStrAccurate(MMSysInfo& info,std::string& str_out)
{
	wxString str;

	str += wxString::Format(" NSTEP = %9d TIME(PS) =  %10.3f TEMP(K) = %7.2f PRESS = %7.1f\n",    
		info.nstep,info.time, info.temp, info.press);
	str += wxString::Format(" Etot = %18.11f EKtot = %18.11f\n",
		info.tot_energy,info.kin_ene);
	str += wxString::Format(" EPtot = %18.11f\n",info.pot_ene);
	str += wxString::Format(" BOND = %18.11f  ANGLE = %18.11f\n",
		info.bond_ene,info.vang_ene);
	str += wxString::Format(" DIHED = %18.11f\n",info.dihed_ene);
	str += wxString::Format(" 1-4 NB = %18.11f  1-4 EEL = %18.11f \n",
		info.vdw_ene_14,info.electr_ene_14);
	str += wxString::Format(" VDWAALS = %18.11f\n",info.vdw_ene);

	if (p_mm_model->p_amber_model->using_pme_potential) 
	{
		str += wxString::Format(" EELEC = %18.11f EHBOND = %18.11f \n",
			info.electr_ene,info.hbond_ene);
	}  
	else
	{
		str += wxString::Format(" EELEC = %18.11f EGB    =  %18.11f\n",
			info.electr_ene,info.gb_ene);
	}
	str += wxString::Format(" RESTRAINT = %18.11f\n",info.constraints_ene);
	if( fabs(info.constraints_ene) > 0)
	{
		str += wxString::Format(" EAMBER (non-restraint) = %18.11f\n",
			(info.pot_ene - info.constraints_ene) );
	}
	if( fabs(info.dv_dlambda) > 0)
	{
		str += wxString::Format(" DV/DL  = %18.11f\n", info.dv_dlambda );
	}
	if( info.volume != 0.0)
	{
		str += wxString::Format(" EKCMT = %18.11f  VIRIAL  = %18.11f\n",
			info.kin_ene_com,info.virial_tot);
		str += wxString::Format(" VOLUME = %18.11f\n",info.volume);
	}
	if( info.kin_ene_solvent != 0.0 ) 
	{
//		str += wxString::Format(" EK_solv = %14.4f \n", info.kin_ene_solvent );
	}
	if( info.epol != 0.0 || info.e3body != 0.0 )
	{
		str += wxString::Format(" EPOLZ  = %18.11f  E3BODY  = %18.11f\n",
			info.epol,info.e3body);
	}
	if(info.polar_ene != 0.0)
	{
		str += wxString::Format(" E_POLAR = %18.11f \n ", info.polar_ene );
	}
	if(info.density != 0.0)
	{
		str += wxString::Format("                                  Density    = %18.11f\n", info.density );
	}
	if(p_mm_model->p_amber_model->using_pme_potential)
	{
		str += wxString::Format(" Ewald error estimate: %18.11e\n", info.pme_err_est ); 
	}
	str += "----------------------------------------------------------------------------\n";
	str_out = str;
}

void HaMolMechMod::PrintLogEne()
{
	std::string str;
	PrintEneStr((*this->p_mm_info),str);
	PrintLog("\n%s\n",str.c_str());
}



bool HaMolMechMod::Print_info(ostream& sout, const int level)
{
	if( level > 2)
	{
		sout << " Excluded Atom List " << endl;
		set<HaAtom*, less<HaAtom*> >::iterator si_itr;
		int nn = p_mm_model->Atoms.size();
		if( nn != p_mm_model->excluded_atom_list.size())
		{
			ErrorInMod("HaMolMechMod::Print_info()",
		" The size of excluded atom list not equal to the number of Atoms");
			return false;
		}
		int i; char buf[256];
		for( i = 0; i < nn; i++)
		{
			HaAtom* aptr = (HaAtom*) p_mm_model->Atoms[i];
			if(aptr == NULL) 
				continue;
			aptr->FillRef(buf);
			PrintLog("%5d %s ", i, buf);
			for(si_itr= p_mm_model->excluded_atom_list[i].begin();
			    si_itr != p_mm_model->excluded_atom_list[i].end (); si_itr++)
			{
				HaAtom* pt1 = *si_itr;
				pt1->FillRef(buf);
				PrintLog(" %s ",buf);
			}
			PrintLog("\n");
		}
	}
	return true;
}

int HaMolMechMod::ProcessEvent(int type, int id)
{
	MolSet* pmset = this->GetMolSet();
	//	PrintLog(" In:HaMolMechMod::ProcessEvent() \n");
	//	PrintLog(" MolSet Name=%s \n",pmset->GetName() );
	//	PrintLog(" Event ID =%d \n",id);

	try
	{
		switch (id)
		{
		case MM_INIT_SIMULATIONS_STEP_2:
			PrintLog("About to Execute InitSimulationsStep2() \n");
			this->p_amber_driver->InitSimulationsStep2();
			break;
		case MM_DRIVER_AMBER_RUN_INTERNAL:
			PrintLog("About to Execute RunInternal_node() \n");
			this->RunInternal_node();
			break;
		case MM_SET_MPI_COMM_ALL_PROCS:
			PrintLog("About to Execute SetMPICommAllProcs() \n");
			this->p_amber_driver->SetMPICommAllProcs();
			this->p_amber_driver->InitParallelDatMod();
			break;
		case MM_MOD_SET_MPI_COMM_SPLIT_2:
			PrintLog("About to Execute SetMPICommSplit2() \n");
			this->SetMPICommSplit2();
			this->p_amber_driver->InitParallelDatMod();
			break;
		case MM_MOD_INIT_MIXED_HAMILTONIAN:
			PrintLog("About to Execute HaMolMechMod::InitMixedHamSimulations_node()  \n");
			this->InitMixedHamSimulations_node(this->p_mm_model);
			break;
		case MM_UPDATE_CONSTR_2:
			PrintLog("About to Execute MolMechModel::UpdateConstraints_2() \n");
			this->GetMolMechModel()->UpdateConstraints_2();
			break;
		default:
			PrintLog("Error in HaMolMechMod::%s() Unknown ID= %d \n", __func__,id);
		}
	}
	catch (std::exception& e)
	{
		PrintLog(" Exception executing HaMolMechMod::%s() \n %s \n", __func__,e.what());
	}
	catch (...)
	{
		PrintLog(" Exception executing HaMolMechMod::::%s() \n",__func__);
	}
	return TRUE;
}

int HaMolMechMod::OnDelAtoms(AtomContainer& del_atoms)
{
	PtrSet pt_set;
	HaAtom* aptr;
	
	AtomIteratorGen aitr(&del_atoms);
	
	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if(aptr != NULL && !pt_set.IsMember(aptr)) pt_set.insert(aptr);
	}

	vector<HaAtom*>::iterator mpitr;
	for(mpitr = p_mm_model->Atoms.begin(); mpitr != p_mm_model->Atoms.end(); )
	{
		if( pt_set.IsMember( *mpitr ) )
		{
			mpitr = p_mm_model->Atoms.erase(mpitr);
		}
		else
			mpitr++;
	}

	set<MMBond,less<MMBond> >::iterator mbitr;
	set<MMBond,less<MMBond> >::iterator mbitr2;
	for( mbitr= p_mm_model->MBonds.begin(); mbitr != p_mm_model->MBonds.end();  mbitr++ )
	{
		if( pt_set.IsMember( (*mbitr).pt1 ) || pt_set.IsMember( (*mbitr).pt2 ) ) 
		{
			mbitr2 = mbitr;
			mbitr2++;
			p_mm_model->MBonds.erase( mbitr );
			mbitr = mbitr2;
		}
		else
			mbitr++;
	}

	set<MMValAngle, less<MMValAngle> >::iterator ang_itr;
	set<MMValAngle, less<MMValAngle> >::iterator ang_itr2;
	for( ang_itr= p_mm_model->ValAngles.begin(); ang_itr != p_mm_model->ValAngles.end(); )
	{
		if( pt_set.IsMember( (*ang_itr).pt1 ) || pt_set.IsMember( (*ang_itr).pt2 ) || pt_set.IsMember( (*ang_itr).pt3 ) )
		{
			ang_itr2 = ang_itr;
			ang_itr2++;
			p_mm_model->ValAngles.erase( ang_itr );
			ang_itr = ang_itr2;
		}
		else
			ang_itr++;
	}

	vector<AtomContact>::iterator vitr;
	for( vitr = p_mm_model->DistConstraints.begin(); vitr != p_mm_model->DistConstraints.end(); )
	{
		if( pt_set.IsMember( (*vitr).pt1 ) || pt_set.IsMember( (*vitr).pt2 ) )
		{
			vitr = p_mm_model->DistConstraints.erase( vitr );
		}
		else
			vitr++;
	}

	vector<MMDihedral>::iterator ditr;
	for( ditr = p_mm_model->Dihedrals.begin(); ditr != p_mm_model->Dihedrals.end(); )
	{
		if( pt_set.IsMember( (*ditr).pt1 ) || pt_set.IsMember( (*ditr).pt2 ) ||
			pt_set.IsMember( (*ditr).pt3 ) || pt_set.IsMember( (*ditr).pt4 ) )
		{
			ditr = p_mm_model->Dihedrals.erase( ditr );
		}
		else
			ditr++;
	}

	for( ditr = p_mm_model->ImprDihedrals.begin(); ditr != p_mm_model->ImprDihedrals.end(); )
	{
		if( pt_set.IsMember( (*ditr).pt1 ) || pt_set.IsMember( (*ditr).pt2 ) ||
			pt_set.IsMember( (*ditr).pt3 ) || pt_set.IsMember( (*ditr).pt4 ) )
		{
			ditr = p_mm_model->ImprDihedrals.erase( ditr );
		}
		else
			ditr++;
	}

	return TRUE;
}

void HaMolMechMod::SetPrefix(const std::string& prefix_par )
{
	std::string prefix_loc = prefix_par;
	boost::trim( prefix_loc );
	if( prefix_loc.empty() ) 
	{
		prefix = GetMolSet()->GetName();
	}
	else
	{
		prefix = prefix_loc;
	}

	p_amber_driver->amber_run_file        = prefix + ".run";
	p_amber_driver->amber_inp_file        = prefix + ".inp";
	p_amber_driver->amber_top_file        = prefix + ".top";
	p_amber_driver->amber_crd_file        = prefix + ".crd";
	p_amber_driver->amber_out_file        = prefix + ".out";
	p_amber_driver->amber_rst_file        = prefix + ".rst";
	p_amber_driver->amber_constr_crd_file = prefix + ".ref";
	p_amber_driver->amber_trj_coord_file  = prefix + ".mdcrd";
	p_amber_driver->amber_trj_vel_file    = prefix + ".mdvel";
	p_amber_driver->amber_trj_ene_file    = prefix + ".mden";

	constr_trj_fname                      = prefix + "_constr_trj.dat";
}
	
std::string HaMolMechMod::GetPrefix() const
{
	return prefix;
}

void HaMolMechMod::SetWrtLogFreq(int wrt_freq)
{
	wrt_log_freq = wrt_freq;
}

void HaMolMechMod::SetWrtRstrtFreq(int wrt_freq)
{
	wrt_rstrt_freq = wrt_freq;
}

void HaMolMechMod::SetWrtMDTrajFreq(int wrt_freq, int save_vel )
{
	wrt_coord_freq = wrt_freq;
	wrt_ener_freq  = wrt_freq;
	if( save_vel )
	{
		wrt_vel_freq = wrt_freq;
	}
	else
	{
		wrt_vel_freq = 0;
	}
}

void HaMolMechMod::SetWrtCoordFreq(int wrt_freq)
{
	wrt_coord_freq = wrt_freq;
}

void HaMolMechMod::SetWrtVelFreq(int wrt_freq)
{
	wrt_vel_freq = wrt_freq;
}

void HaMolMechMod::SetWrtEnerFreq(int wrt_freq)
{
	wrt_ener_freq  = wrt_freq;
}

void HaMolMechMod::SetWrtConstrFreq(int wrt_freq)
{
	wrt_constr_freq = wrt_freq;
}

void HaMolMechMod::SetRestrtFileFormat( const CrdFormatParam& format )
{
	write_coord_format = format;
}

void HaMolMechMod::SetMDtrajFileFormat( const CrdFormatParam& format)
{
	traj_wrt_format = format;
}

void HaMolMechMod::SetEneMinMethod( const EneMinMethod& method )
{
	min_type = method;
}

void HaMolMechMod::SetMaxNumMinimSteps( int max_num_minim_steps_new )
{
	max_num_minim_steps = max_num_minim_steps_new;
}

void HaMolMechMod::SetNumSteepDescentSteps( int nsteps )
{
	num_steep_descent_steps = nsteps;
}

void HaMolMechMod::SetInitMinStep( double init_min_step_new ) 
{
	init_min_step = init_min_step_new;
}

void HaMolMechMod::SetGradCnvrgVal(double grad_cnvrg_val_new )
{
	grad_cnvrg_val = grad_cnvrg_val_new;
}

void HaMolMechMod::SetZMatMin( bool set_par )
{
	zmat_min = set_par;
}

bool HaMolMechMod::IsZMatMin() const
{
	return zmat_min;
}

void HaMolMechMod::SetNumMDSteps( int num_md_steps_new )
{
	num_md_steps = num_md_steps_new;
}

void HaMolMechMod::SetRemoveInitRBMotion( int remove_init_motion )
{
	remove_init_rb_motion_flag = remove_init_motion;
}

void HaMolMechMod::SetRemoveRBMotionFreq( int freq )
{
	remove_rb_motion_freq = freq;
}

void HaMolMechMod::SetStartVelMethod( const StartVelMethod& start_vel_method_new)
{
	start_vel_method = start_vel_method_new;
}

void HaMolMechMod::SetStartTime( double start_time_new )
{
	start_time = start_time_new;
}

void HaMolMechMod::SetMDTimeStep( double md_time_step_new )
{
	md_time_step = md_time_step_new;
}

void HaMolMechMod::SetNBListUpdateFreq( int freq )
{
	nb_list_update_freq = freq;
}

void HaMolMechMod::SetPerBoundaryCondType( const PerBoundaryCondType& type)
{
	period_bcond = type;
}

void HaMolMechMod::SetRefTemp( double ref_temp_new )
{
	ref_temp = ref_temp_new;
}

void HaMolMechMod::SetInitTemp( double init_temp_new )
{
	init_temp = init_temp_new; 
}

void HaMolMechMod::SetLangevinDumpConst( double langevin_dump_const_new )
{
	langevin_dump_const = langevin_dump_const_new;
}

void HaMolMechMod::SetRandomSeed( int random_seed_new )
{
	random_seed = random_seed_new;
}

void HaMolMechMod::SetScaleInitVel( double scale_init_vel_new )
{
	scale_init_vel_new = scale_init_vel;
}

void HaMolMechMod::SetTempCtrlMethod( const TempCtrlMethod& method)
{
	temp_control_method = method;
}

void HaMolMechMod::SetTempDeviation( double temp_deviation_new )
{
	temp_deviation = temp_deviation_new;
}

void HaMolMechMod::SetTempRelaxTimeSolute( double temp_relax_time_solute_new )
{
	temp_relax_time_solute = temp_relax_time_solute_new;
}

void HaMolMechMod::SetTempRelaxTimeSolvent( double temp_relax_time_solvent_new )
{
	temp_relax_time_solvent = temp_relax_time_solvent_new;
}

void HaMolMechMod::SetVelLimit( double vel_limit_new )
{
	vel_limit = vel_limit_new;
}

void HaMolMechMod::SetPressureRegMethod( const PressureRegMethod& method)
{
	pressure_reg_method = method;
}

void HaMolMechMod::SetRefPressure( double ref_pressure_new )
{
	ref_pressure = ref_pressure_new;
}

void HaMolMechMod::SetCompressibility( double compressibility_new )
{
	compressibility = compressibility_new;
}
	
void HaMolMechMod::SetPressRelaxTime( double press_relax_time_new )
{
	press_relax_time = press_relax_time_new;
}

void HaMolMechMod::SetShakeConstr( const MMShakeParam& shake_constr_new)
{
	shake_constr = shake_constr_new;
}

void HaMolMechMod::SetShakeTol( double shake_tol_new )
{
	shake_tol = shake_tol_new;
}

#if defined(HA_NOGUI)
void HaMolMechMod::OnChangePeriodicity()
{
	period_bcond.SetCompatValue();
	p_mm_model->electr_method.SetCompatValue();
}
#endif

void HaMolMechMod::TestSaveAmoebaTopFile1()
{
	MolSet* pmset = GetCurMolSet();
	HaMolMechMod* p_mm_mod = pmset->GetMolMechMod(1);
	p_mm_mod->InitMolMechModel(ForceFieldType::AMOEBA);
	p_mm_mod->p_amber_driver->SaveAmberTopFile();
}

void HaMolMechMod::TestSaveAmoebaTopFile2()
{
	PrintLog(" HaMolMechMod::TestSaveAmoebaTopFile2()  Empty function \n");
//	MolSet* pmset = GetCurMolSet();
//	HaMolMechMod* p_mm_mod = pmset->GetMolMechMod(1);
//	p_mm_mod->InitMolMechModel(ForceFieldType::AMOEBA);
	
}


MinEneMod::MinEneMod(HaMolMechMod* p_mm_mod_new)
{
	p_mm_mod = p_mm_mod_new;
}
	
MinEneMod::~MinEneMod()
{

}

MDSimMod::MDSimMod(HaMolMechMod* p_mm_mod_new)
{
	p_mm_mod = p_mm_mod_new;
}
	
MDSimMod::~MDSimMod()
{

}


HaVec_int TISimMod::allowed_num_lmb;

TISimMod::TISimMod(HaMolMechMod* p_mm_mod_new)
{
	ti_sync_freq = 100;
	klambda_ti = 1.0;

	num_lmb = 3;
	cur_idx_lmb = 0;
	max_idx = -1;

	num_equilib_points = 2000;

	file_prefix = "";
	file_dvdl = NULL;
	p_mm_mod = p_mm_mod_new;

	allowed_num_lmb.clear();
	allowed_num_lmb.append(1);
	allowed_num_lmb.append(2);
	allowed_num_lmb.append(3);
	allowed_num_lmb.append(5);
	allowed_num_lmb.append(7);
	allowed_num_lmb.append(9);
	allowed_num_lmb.append(12);
}

TISimMod::~TISimMod()
{

}

double TISimMod::GetLambdaByIdx(int idx)
{
	double lmb = 0.0;
	
	if(idx < 0 || idx >= num_lmb)
	{
		PrintLog("Error in TISimMod::GetLambdaByIdx() \n");
		PrintLog("Invalid lambda index = %d for  num_lmb = %d \n",idx,num_lmb);
		return lmb;
	}
	int i;
	
	if(num_lmb == 1)
	{
		lmb = 0.5;
	}
	else if(num_lmb == 2)
	{
		if(idx == 0) lmb = 0.21133;
		if(idx == 1) lmb = 0.78867;
	}
	else if(num_lmb == 3)
	{
		if(idx == 0) lmb = 0.11271;
		if(idx == 1) lmb = 0.5;
		if(idx == 2) lmb = 0.88729;
	}
	else if(num_lmb == 5)
	{
		if(idx == 0) lmb = 0.04692;
		if(idx == 1) lmb = 0.23077;
		if(idx == 2) lmb = 0.5;
		if(idx == 3) lmb = 0.76923;
		if(idx == 4) lmb = 0.95308;
	}
	else if(num_lmb == 7)
	{
		if(idx == 0) lmb = 0.02545;
		if(idx == 1) lmb = 0.12924;
		if(idx == 2) lmb = 0.29708;
		if(idx == 3) lmb = 0.5;
		if(idx == 4) lmb = 0.70292;
		if(idx == 5) lmb = 0.87076;
		if(idx == 6) lmb = 0.97455;
	}
	else if(num_lmb == 9)
	{
		if(idx == 0) lmb = 0.01592;
		if(idx == 1) lmb = 0.08198;
		if(idx == 2) lmb = 0.19331;
		if(idx == 3) lmb = 0.33787;
		if(idx == 4) lmb = 0.5;
		if(idx == 5) lmb = 0.66213;
		if(idx == 6) lmb = 0.80669;
		if(idx == 7) lmb = 0.91802;
		if(idx == 8) lmb = 0.98408;
	}
	else if(num_lmb == 12)
	{
		if(idx == 0) lmb = 0.00922;
		if(idx == 1) lmb = 0.04794;
		if(idx == 2) lmb = 0.11505;
		if(idx == 3) lmb = 0.20634;
		if(idx == 4) lmb = 0.31608;
		if(idx == 5) lmb = 0.43738;
		if(idx == 6) lmb = 0.56262;
		if(idx == 7) lmb = 0.68392;
		if(idx == 8) lmb = 0.79366;
		if(idx == 9) lmb = 0.88495;
		if(idx == 10) lmb = 0.95206;
		if(idx == 11) lmb = 0.99078;
	}
	return lmb;
}

double TISimMod::GetIntegWtByIdx(int idx)
{
	double wt = 0.0;
	
	if(idx < 0 || idx >= num_lmb)
	{
		PrintLog("Error in TISimMod::GetIntegWeightByIdx() \n");
		PrintLog("Invalid lambda index = %d for  num_lmb = %d \n",idx,num_lmb);
		return wt;
	}
	int i;
	
	if(num_lmb == 1)
	{
		wt = 1.0;
	}
	else if(num_lmb == 2)
	{
		if(idx == 0) wt = 0.5;
		if(idx == 1) wt = 0.5;
	}
	else if(num_lmb == 3)
	{
		if(idx == 0) wt = 0.27777;
		if(idx == 1) wt = 0.44444;
		if(idx == 2) wt = 0.27777;
	}
	else if(num_lmb == 5)
	{
		if(idx == 0) wt = 0.11846;
		if(idx == 1) wt = 0.23931;
		if(idx == 2) wt = 0.28444;
		if(idx == 3) wt = 0.23931;
		if(idx == 4) wt = 0.11846;
	}
	else if(num_lmb == 7)
	{
		if(idx == 0) wt = 0.06474;
		if(idx == 1) wt = 0.13985;
		if(idx == 2) wt = 0.19091;
		if(idx == 3) wt = 0.20897;
		if(idx == 4) wt = 0.19091;
		if(idx == 5) wt = 0.13985;
		if(idx == 6) wt = 0.06474;
	}
	else if(num_lmb == 9)
	{
		if(idx == 0) wt = 0.04064;
		if(idx == 1) wt = 0.09032;
		if(idx == 2) wt = 0.13031;
		if(idx == 3) wt = 0.15617;
		if(idx == 4) wt = 0.16512;
		if(idx == 5) wt = 0.15617;
		if(idx == 6) wt = 0.13031;
		if(idx == 7) wt = 0.09032;
		if(idx == 8) wt = 0.04064;
	}
	else if(num_lmb == 12)
	{
		if(idx == 0) wt = 0.02359;
		if(idx == 1) wt = 0.05347;
		if(idx == 2) wt = 0.08004;
		if(idx == 3) wt = 0.10158;
		if(idx == 4) wt = 0.11675;
		if(idx == 5) wt = 0.12457;
		if(idx == 6) wt = 0.12457;
		if(idx == 7) wt = 0.11675;
		if(idx == 8) wt = 0.10158;
		if(idx == 9)  wt = 0.08004;
		if(idx == 10) wt = 0.05347;
		if(idx == 11) wt = 0.02359;
	}
	return wt;
}

double TISimMod::GetCurLambda()
{
	return GetLambdaByIdx(cur_idx_lmb);
}

void TISimMod::SetNumEqPoints(int num_eq_pt_new)
{
	num_equilib_points = num_eq_pt_new;
}

void TISimMod::SetNumLambda(int num_lmb_new)
{
	int not_allowed = TRUE;
	int i;
	for(i = 0; i < allowed_num_lmb.size(); i++)
	{
		if(num_lmb_new == allowed_num_lmb[i]) not_allowed = FALSE;
	}

	if(not_allowed)
	{
		PrintLog("Warning in TISimMod::SetNumLambda() \n");
		PrintLog("Invalid num_lmb = %d \n",num_lmb_new);
		PrintLog("Allowed values = 1,2,3,5,7,9 and 12 \n");
		PrintLog(" Set to standard value num_lmb = 5 \n");
		num_lmb = 5;
	}
	else
	{
		num_lmb = num_lmb_new;
	}
}

int  TISimMod::GetNumLambda()
{
	return num_lmb;
}

void TISimMod::SetCurIdxLambda(int cur_idx_lmb_new)
{
	cur_idx_lmb = cur_idx_lmb_new;
}

void TISimMod::SetMaxLambdaIdx(int max_idx_new)
{
	max_idx = max_idx_new;
}

void TISimMod::SetTI_OutputFileNames()
{
	MMDriverAmber* p_amber_driver = p_mm_mod->p_amber_driver;
	if(!p_amber_driver->master) return;

	wxString cur_prefix = GetCurFilePrefix();
	
	p_amber_driver->amber_out_file        = cur_prefix + ".out";
	p_amber_driver->amber_rst_file        = cur_prefix + ".rst";
	p_amber_driver->amber_trj_coord_file  = cur_prefix + ".mdcrd";
	p_amber_driver->amber_trj_vel_file    = cur_prefix + ".mdvel";
	p_amber_driver->amber_trj_ene_file    = cur_prefix + ".mden";
}

void TISimMod::CollectForceAndEneTI(HaVec_double& sys_info)
{
	int TAG_F = 301;
	int TAG_E = 302;
	MPI_Status status;
	int ires;
	int i,j;

	double ep_1;
	double ep_2;

	double scale;
	HaVec_double buf_arr;

	double lambda = GetCurLambda();

	MMDriverAmber* p_amber_driver = p_mm_mod->p_amber_driver;
	AmberMMModel*  p_amber_model  = p_amber_driver->p_amber_model;
	
	if( p_mm_mod->inter_model_rank == 0 )
	{
		if( (klambda_ti  - 1.0) < 0.01 )
		{
			scale = (1.0 - lambda );
		}
		else
		{
			scale = pow((1.0 - lambda), klambda_ti);
		}
	}
	else
	{
		if( (klambda_ti  - 1.0) < 0.01 )
		{
			scale = lambda;
		}
		else
		{
			scale = 1.0 - pow((1.0 - lambda),klambda_ti);
		}
	}

	int atm_lst_idx;
	for(atm_lst_idx = 0; atm_lst_idx < p_amber_driver->my_atm_cnt; atm_lst_idx++)
	{
		j = p_amber_driver->gbl_my_atm_lst[atm_lst_idx] - 1;
		p_amber_driver->atm_frc[3*j]   *= scale;
		p_amber_driver->atm_frc[3*j+1] *= scale;
		p_amber_driver->atm_frc[3*j+2] *= scale;
	}

	if( p_amber_driver->master )
	{
		if(p_mm_mod->inter_model_rank == 0) 
		{
			ep_1 = sys_info[MMDriverAmber::SI_POT_ENE]; 
			MPI_Bcast(&ep_2,1,MPI_DOUBLE,1,p_mm_mod->inter_model_comm);
		}
		if(p_mm_mod->inter_model_rank == 1) MPI_Bcast(&sys_info[MMDriverAmber::SI_POT_ENE],1,MPI_DOUBLE,1,p_mm_mod->inter_model_comm);
	}

//	PrintLogMDOUT(" pot_ene = %16.9f \n",sys_info[SI_POT_ENE]);
//	PrintLogMDOUT(" vdw_ene = %16.9f \n",sys_info[SI_VDW_ENE]);
//	PrintLogMDOUT(" elect_ene = %16.9f \n",sys_info[SI_ELECT_ENE]);
//	PrintLogMDOUT(" bond_ene = %16.9f \n",sys_info[SI_BOND_ENE]);
//	PrintLogMDOUT(" angle_ene = %16.9f \n",sys_info[SI_ANGLE_ENE]);
//	PrintLogMDOUT(" dihed_ene = %16.9f \n",sys_info[SI_DIHEDRAL_ENE]);
//	PrintLogMDOUT(" vdw_14_ene = %16.9f \n",sys_info[SI_VDW_14_ENE]);
//	PrintLogMDOUT(" elect_14_ene = %16.9f \n",sys_info[SI_ELECT_14_ENE]);

	sys_info[MMDriverAmber::SI_POT_ENE]       *= scale;
	sys_info[MMDriverAmber::SI_VDW_ENE]       *= scale;  
	sys_info[MMDriverAmber::SI_ELECT_ENE]     *= scale; 
	sys_info[MMDriverAmber::SI_HBOND_ENE]     *= scale; 
	sys_info[MMDriverAmber::SI_BOND_ENE]      *= scale; 
	sys_info[MMDriverAmber::SI_ANGLE_ENE]     *= scale; 
	sys_info[MMDriverAmber::SI_DIHEDRAL_ENE]  *= scale; 
	sys_info[MMDriverAmber::SI_VDW_14_ENE]    *= scale; 
	sys_info[MMDriverAmber::SI_ELECT_14_ENE]  *= scale; 
	sys_info[MMDriverAmber::SI_RESTRAINT_ENE] *= scale;
	sys_info[MMDriverAmber::SI_VIR_0]         *= scale;
	sys_info[MMDriverAmber::SI_VIR_1]         *= scale;
	sys_info[MMDriverAmber::SI_VIR_2]         *= scale;

	if( p_amber_driver->numtasks > 1 ) p_amber_driver->AllGatherVec(p_amber_driver->atm_frc);
	
	MPI_Barrier(MPI_COMM_WORLD);

	if( p_amber_driver->mytaskid == 0)
	{
		int n3 = 3*p_amber_model->natom;
		if(buf_arr.size() != n3) buf_arr.resize(n3);

		if(p_mm_mod->inter_model_rank == 1) 
		{
			ires = MPI_Bcast(p_amber_driver->atm_frc.v(),n3,MPI_DOUBLE,1,p_mm_mod->inter_model_comm);
			ires = MPI_Bcast(buf_arr.v(),n3,MPI_DOUBLE,0,p_mm_mod->inter_model_comm);
		}
		else
		{
			ires = MPI_Bcast(buf_arr.v(),n3,MPI_DOUBLE,1,p_mm_mod->inter_model_comm);
			ires = MPI_Bcast(p_amber_driver->atm_frc.v(),n3,MPI_DOUBLE,0,p_mm_mod->inter_model_comm);
		}
		p_amber_driver->atm_frc += buf_arr;
	}

	if( buf_arr.size() < MMDriverAmber::SI_CNT ) buf_arr.resize(MMDriverAmber::SI_CNT);

	MPI_Barrier(MPI_COMM_WORLD);

	if(p_mm_mod->inter_model_rank == 1) 
	{
		ires = MPI_Bcast(sys_info.v(),MMDriverAmber::SI_CNT,MPI_DOUBLE,1,p_mm_mod->inter_model_comm);
		ires = MPI_Bcast(buf_arr.v(), MMDriverAmber::SI_CNT,MPI_DOUBLE,0,p_mm_mod->inter_model_comm);
	}
	else
	{
		ires = MPI_Bcast(buf_arr.v(), MMDriverAmber::SI_CNT,MPI_DOUBLE,1,p_mm_mod->inter_model_comm);
		ires = MPI_Bcast(sys_info.v(),MMDriverAmber::SI_CNT,MPI_DOUBLE,0,p_mm_mod->inter_model_comm);
	}

	sys_info[MMDriverAmber::SI_POT_ENE]       += buf_arr[MMDriverAmber::SI_POT_ENE];
	sys_info[MMDriverAmber::SI_VDW_ENE]       += buf_arr[MMDriverAmber::SI_VDW_ENE]; 
	sys_info[MMDriverAmber::SI_ELECT_ENE]     += buf_arr[MMDriverAmber::SI_ELECT_ENE];
	sys_info[MMDriverAmber::SI_HBOND_ENE]     += buf_arr[MMDriverAmber::SI_HBOND_ENE];
	sys_info[MMDriverAmber::SI_BOND_ENE]      += buf_arr[MMDriverAmber::SI_BOND_ENE]; 
	sys_info[MMDriverAmber::SI_ANGLE_ENE]     += buf_arr[MMDriverAmber::SI_ANGLE_ENE]; 
	sys_info[MMDriverAmber::SI_DIHEDRAL_ENE]  += buf_arr[MMDriverAmber::SI_DIHEDRAL_ENE]; 
	sys_info[MMDriverAmber::SI_VDW_14_ENE]    += buf_arr[MMDriverAmber::SI_VDW_14_ENE]; 
	sys_info[MMDriverAmber::SI_ELECT_14_ENE]  += buf_arr[MMDriverAmber::SI_ELECT_14_ENE]; 
	sys_info[MMDriverAmber::SI_RESTRAINT_ENE] += buf_arr[MMDriverAmber::SI_RESTRAINT_ENE];
	sys_info[MMDriverAmber::SI_VIR_0]         += buf_arr[MMDriverAmber::SI_VIR_0];
	sys_info[MMDriverAmber::SI_VIR_1]         += buf_arr[MMDriverAmber::SI_VIR_1];
	sys_info[MMDriverAmber::SI_VIR_2]         += buf_arr[MMDriverAmber::SI_VIR_2];

	if( p_amber_driver->master && (p_mm_mod->inter_model_rank == 0) )
	{
		if( (klambda_ti  - 1.0) < 0.01 )
		{
			sys_info[MMDriverAmber::SI_DV_DLAMBDA] = ep_2 - ep_1;
		}
		else
		{
			sys_info[MMDriverAmber::SI_DV_DLAMBDA] = klambda_ti * (ep_2 - ep_1) * pow((1.0 - p_mm_mod->lambda_ti),(klambda_ti - 1));
		}
	}

	if(p_amber_driver->numtasks > 1)
	{
		p_amber_driver->BcastFrc(p_amber_driver->driver_mpi_comm);
//		parallel_mod_mp_distribute_crds_proxy_(&natom, &my_atm_cnt, atm_frc.v());
	}
}

void TISimMod::SincCrdAndVelTI()
{
	MMDriverAmber* p_amber_driver = p_mm_mod->p_amber_driver;
	AmberMMModel*  p_amber_model  = p_amber_driver->p_amber_model;

#if defined(HARLEM_MPI)
	MPI_Barrier(MPI_COMM_WORLD);
	if( p_amber_driver->numtasks > 1 && p_mm_mod->inter_model_rank == 0 && !p_amber_driver->all_crds_valid) 
	{
		p_amber_driver->AllGatherVec(p_amber_driver->atm_crd);
	}
	
	if( p_amber_driver->numtasks > 1 && p_mm_mod->inter_model_rank == 0 && !p_amber_driver->all_vels_valid) 
	{
		p_amber_driver->AllGatherVec(p_amber_driver->atm_vel);
	}
	
	if (p_amber_driver->ntp > 0 && p_amber_model->natc > 0)
	{
		p_amber_driver->AllGatherVec(p_amber_model->atm_xc);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	if( p_amber_driver->master )
	{
		p_amber_driver->BcastCrd(p_mm_mod->inter_model_comm);
		p_amber_driver->BcastVel(p_mm_mod->inter_model_comm);
		if (p_amber_driver->ntp > 0) p_amber_driver->BcastPBox(p_mm_mod->inter_model_comm);
		if (p_amber_driver->ntp > 0 && p_amber_model->natc > 0) 
			MPI_Bcast(p_amber_model->atm_xc.v(),3*p_amber_model->natom,MPI_DOUBLE,0,p_mm_mod->inter_model_comm);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(  p_amber_driver->numtasks > 1 && p_mm_mod->inter_model_rank == 1 )
	{
		p_amber_driver->BcastCrd(p_amber_driver->driver_mpi_comm);
		p_amber_driver->BcastVel(p_amber_driver->driver_mpi_comm);
		if( p_amber_driver->ntp > 0) p_amber_driver->BcastPBox(p_amber_driver->driver_mpi_comm);
		if ( p_amber_driver->ntp > 0 && p_amber_model->natc > 0) 
			MPI_Bcast(p_amber_model->atm_xc.v(),3*p_amber_model->natom,MPI_DOUBLE,0,p_amber_driver->driver_mpi_comm);
	}
	p_amber_driver->all_crds_valid = TRUE;
	p_amber_driver->all_vels_valid = TRUE;
#endif
}


int TISimMod::ComputeDvDlAvg()
{	
	if(num_lmb < 1) 
	{
		PrintLog("Error in ComputeDvDlAvg()  num_lmb = %d  < 1 \n",num_lmb);
		return FALSE;
	}

	dvdl_avg.resize(num_lmb);

	char buf[256];
	int i,j;
	int ires;
	for(i = 0; i < num_lmb; i++)
	{
		std::string dvdl_file_name = GetFilePrefixIdx(i);
		dvdl_file_name += "_dvdl.dat";
		FILE* finp = fopen(dvdl_file_name.c_str(),"r");
		if(finp == NULL)
		{
			PrintLog("Error in TISimMod::ComputeDvDlAvg() opening DvDl file \n");
			PrintLog(" lambda index = %d   File Name %s \n",num_lmb,dvdl_file_name.c_str());
			return FALSE;
		}
		int ip;
		double val;
		for(j = 0; j < num_equilib_points; j++)
		{
			ires = fscanf(finp,"%d %lf",&ip,&val);
			if( ires == EOF) 
			{
				fclose(finp);
				PrintLog("Error in TISimMod::ComputeDvDlAvg() skipping equilibration points \n");
				PrintLog(" lambda index = %d dvdl file = %s  assumed equilibration points number = %d \n",
					      i,dvdl_file_name.c_str(),num_equilib_points);
				return FALSE;
			}
		}
		int npt = 0;
		double dsum = 0.0;
		for(;;)
		{
			ires = fscanf(finp,"%d %lf",&ip,&val);
			if( ires == EOF) break;
			npt++;
			dsum += val;
		}
		fclose(finp);
		if( npt == 0 )
		{	
			PrintLog("Error in TISimMod::ComputeDvDlAvg() no points to compute average dvdl \n");
			PrintLog("lambda idx = %d   dvdl file name = %s \n",i,dvdl_file_name.c_str());
			return FALSE;
		}
		dvdl_avg[i] = dsum/npt;
	}
	return TRUE;
}

double TISimMod::CalcDeltaG(int recalc_dvdl_avg)
{
	int ires;
	if(recalc_dvdl_avg || dvdl_avg.size() != num_lmb )
	{
		ires = ComputeDvDlAvg();
		if( !ires )
		{
			PrintLog("Error in TISimMod::CalcDeltaG() computing dv/dl averages \n");
			return 0.0;
		}
	}
	int i;
	delta_g = 0.0;
	PrintLog(" idx     lambda       <dv/dl>           wt            DDG   \n");
	for( i = 0; i < num_lmb; i++)
	{
		double wt  = GetIntegWtByIdx(i);
		double lmb = GetLambdaByIdx(i); 
		double val = dvdl_avg[i];
		delta_g += wt*val;
		PrintLog(" %d  %12.6f  %12.6f    %12.6f  %12.6f \n",
			       i,lmb,val,wt,wt*val);
	}
	PrintLog(" TISimMod::CalcDeltaG() delta G = %12.6f (kcal/mol) \n",delta_g);
	return delta_g;
}

void TISimMod::ReduceDvDlData(int n_avg,const char* file_name)
{
	if(num_lmb < 1) 
	{
		PrintLog("Error in ReduceDvDlData()  num_lmb = %d  < 1 \n",num_lmb);
		return;
	}
	if(n_avg < 1) 
	{
		PrintLog("Error in ReduceDvDlData()  num_avg = %d  < 1  \n  Set to 1 \n",n_avg);
		n_avg = 1;
	}
	
	char buf[256];
	int i,j,ires;

	vector<FILE*>files_inp;
	files_inp.resize(num_lmb);

	int no_files = TRUE;
	for(i = 0; i < num_lmb; i++)
	{
		std::string dvdl_file_name = GetFilePrefixIdx(i);
		dvdl_file_name += "_dvdl.dat";
		files_inp[i] = fopen(dvdl_file_name.c_str(),"r");
		if( files_inp[i] != NULL) no_files = FALSE;
		
	}
	if( no_files )
	{
		PrintLog("Error in TISimMod::ReduceDvDlData() opening DvDl files \n");
		for(i = 0; i < num_lmb; i++)
		{
			std::string dvdl_file_name = GetFilePrefixIdx(i);
			dvdl_file_name += "_dvdl.dat";
			PrintLog(" lambda index = %d   File Name %s \n",i,dvdl_file_name.c_str());
		}
		return;
	}

	FILE* fout = fopen(file_name,"w");
	if(fout == NULL)
	{
		PrintLog("Error in TISimMod::ReduceDvDlData() opening output file %s \n",file_name);
		for(i= 0; i < num_lmb; i++)
		{
			if(files_inp[i] != NULL) fclose(files_inp[i]);
		}
		return;
	}
	
	fprintf(fout," lambda idx ");
	for( i = 0; i < num_lmb; i++)
	{
		fprintf(fout,"  %12d    ",i);
	}
	fprintf(fout,"\n");
	fprintf(fout,"  lambda  ");
	for( i = 0; i < num_lmb; i++)
	{
		fprintf(fout," %12.6f  ",GetLambdaByIdx(i));
	}
	fprintf(fout,"\n");
	fprintf(fout,"    n_pt \n");

	int ip;
	double val;
	HaVec_double val_avg;
	val_avg.resize(num_lmb);
	int n_pt = 0;
	for(;;)
	{
		val_avg = 0.0;
		int to_stop = FALSE; 
		for(j = 0; j < n_avg; j++)
		{
			n_pt++;
			int all_files_closed = TRUE;
			for(i = 0; i < num_lmb; i++)
			{
				if( files_inp[i] != NULL )
				{
					ires = fscanf(files_inp[i],"%d %lf",&ip,&val);
					if( ires == EOF) 
					{
						fclose(files_inp[i]);
						files_inp[i] = NULL;
					}
					else
					{
						val_avg[i] += val;
						all_files_closed = FALSE;
					}
				}
			}
			if( all_files_closed )
			{
				to_stop = TRUE;
				break;
			}
		}
		if(to_stop)
		{
			fclose(fout);
			return;
		}
		else
		{
			fprintf(fout," %9d ",n_pt);
			for(i = 0; i < num_lmb; i++)
			{
				if(files_inp[i] != NULL)
				{
					fprintf(fout," %12.6f ",val_avg[i]/n_avg);
				}
				else
				{
					fprintf(fout," %12.6f ",0.0);
				}
			}
			fprintf(fout,"\n");
		}
	}
}

std::string TISimMod::GetFilePrefixIdx(int idx)
{
	if(file_prefix.empty())
	{
		if( p_mm_mod->inter_model_rank == 0 )
		{
			file_prefix = (p_mm_mod->GetMolSet())->GetName();
		}
		else
		{
			file_prefix = "ti_simm_rank_1";
		}
	}

	std::string cur_prefix = file_prefix.c_str();
	wxString str_idx = wxString::Format("%d",idx);
	str_idx = str_idx.Strip(wxString::both);
	wxString str_num = wxString::Format("%d",num_lmb);
	str_num = str_num.Strip(wxString::both);
	cur_prefix += "_lmb_";
	cur_prefix += str_idx.c_str();
	cur_prefix += "_";
	cur_prefix += str_num.c_str();
	return cur_prefix;
}

std::string TISimMod::GetCurFilePrefix()
{
	std::string cur_prefix = GetFilePrefixIdx(cur_idx_lmb);
	return cur_prefix;
}

///////////////////////////////////////////////////////////////////////////////
MDTrajectory::MDTrajectory(MolSet* new_pmset)
{
	CrdFileName="";
	VelFileName="";
	EneFileName="";
	
	CrdFile=NULL;
	VelFile=NULL;
	EneFile=NULL;
#ifndef _MSC_VER
	xtctrj=NULL;
#endif
	format=FormatUnk;
	
	frame=0;
	pmset=new_pmset;
	
	if( pmset == NULL) 
	{
		PrintLog("Error in MDTrajectory::MDTrajectory() ");
		PrintLog("pmset == NULL");
		return;
	}

	if( pmset->per_bc->IsSet() )
	{
		AMBER_CRD_withBoxSize = true;
	}
	else
	{
		AMBER_CRD_withBoxSize = false;
	}

	HaAtom* aptr;

	AtomIteratorMolSet aitr(pmset);
	for(aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom())
	{
		Atoms.push_back(aptr);
	}

}
MDTrajectory::~MDTrajectory()
{
	Close();
}

int MDTrajectory::Open()
{
	PrintLog("MDTrajectory::Open\n");
	wxString wxCrd(CrdFileName.c_str());
	if(!wxCrd.IsEmpty())
	{
		wxFileName wxFN(wxCrd);
		wxString wxFileExt=wxFN.GetExt();
		if(wxFileExt.IsSameAs("mdcrd",false)||wxFileExt.IsSameAs("crd",false))
		{
			format=AMBER_CRD;
			if(OpenAMBER_CRD(wxCrd.mb_str())==FALSE) return FALSE;
		}
		else if(wxFileExt.IsSameAs("xtc",false))
		{
			format=GMX_XTC;
			if(OpenXTC(wxCrd.c_str())==FALSE) return FALSE;
		}
	}
	return TRUE;
}
int MDTrajectory::Close()
{
	if(CrdFile!=NULL)
		fclose(CrdFile);
	CrdFile=NULL;
	
	if(VelFile!=NULL)
		fclose(VelFile);
	VelFile=NULL;
	
	if(EneFile!=NULL)
		fclose(EneFile);
	EneFile=NULL;
	
	return TRUE;
}
int MDTrajectory::ReadNextFrame()
{
	if(format==AMBER_CRD)
		return ReadNextFrameAMBER_CRD();
	else if(format==GMX_XTC)
		return ReadNextFrameXTC();
	return FALSE;
}

int MDTrajectory::OpenAMBER_CRD(const char* filename)
{
	PrintLog("Open AMBER MDCRD trajectory (%s)\n",filename);
	CrdFile = fopen(filename,"r");
	if(CrdFile == NULL)
	{
		PrintLog("Can't find MD trajectory file %s\n",filename );
		return FALSE;
	}	

	int i;
	char buf[256];
	
	//read title or whats there
	fgets(buf,128,CrdFile);
	
	frame=0;
	return TRUE;
}
int MDTrajectory::ReadNextFrameAMBER_CRD()
{
	char buf[256];
	
	PrintLog("MDTrajectory::ReadNextFrameAMBER_CRD\n");
	
	int iread = FALSE;
	
	if(CrdFile != NULL)
	{
		iread = TRUE;
		int npt =Atoms.size();
		HaVec_double darr(3*npt);
		// Reading atom coordinate
		int i, ires;
		double box_x,box_y,box_z;
		
		for(i = 0 ; i < 3*npt; i++)
		{
			ires = fscanf(CrdFile,"%lf",&darr(i+1));
			if(ires == EOF)
			{
				PrintLog("End of AMBER Trajectory file \n");
				return FALSE;
			}
		}
		
		if( AMBER_CRD_withBoxSize )
		{
			ires = fscanf(CrdFile,"%lf %lf %lf",&box_x,&box_y,&box_z);
			if(ires == EOF)
			{
				PrintLog("End of AMBER Trajectory file \n");
				return FALSE;
			}

			pmset->per_bc->SetBox( box_x, box_y, box_z );
		}
		
		for(i = 0; i < npt; i++ )
		{
			Atoms[i]->SetX_Ang( darr[3*i]   );
			Atoms[i]->SetY_Ang( darr[3*i+1] );
			Atoms[i]->SetZ_Ang( darr[3*i+2] );
		}
		frame++;
		PrintLog("frame %d\n",frame);
	}
	else
	{
		PrintLog("MDTrajectory::ReadNextFrameAMBER_CRD there is no CrdFile\n");
	}
	return iread;
}
int MDTrajectory::RefreshAllViews()
{
	pmset->RefreshAllViews( RFRefresh | RFColour | RFApply );
	return TRUE;
}

int MDTrajectory::OpenXTC(const char* filename)
{
#ifndef _MSC_VER
	PrintLog("Open GROMACS XTC trajectory (%s)\n",filename);
	if(xtctrj!=NULL)
	{
		PrintLog("One XTC file is already open, close it first\n");
		return FALSE;
	}
	xtctrj=new XTCTraj();
	xtctrj->FileName=filename;
	if(xtctrj->Open()==FALSE)
	{
		PrintLog("Can't open MD trajectory file %s\n",filename );
		return FALSE;
	}
#else
	PrintLog("Not implemented\n");
	return FALSE;
#endif
	return TRUE;
}
int MDTrajectory::ReadNextFrameXTC()
{
#ifndef _MSC_VER
	if(xtctrj==NULL)
	{
		PrintLog("One XTC file is not open\n");
		return FALSE;
	}
	xtctrj->ReadNextFrame();
	int i;
	int npt =Atoms.size();
	for(i = 0; i < npt; i++ )
	{
		Atoms[i]->SetX( 10.0*xtctrj->r->GetVal_idx0(3*i)   );
		Atoms[i]->SetY( 10.0*xtctrj->r->GetVal_idx0(3*i+1) );
		Atoms[i]->SetZ( 10.0*xtctrj->r->GetVal_idx0(3*i+2) );
	}

	double px = 10.0*xtctrj->box[0][0];
	double py = 10.0*xtctrj->box[1][1];
	double pz = 10.0*xtctrj->box[2][2];
	
	pmset->per_bc->SetBox( px, py, pz );

	frame++;
	PrintLog("frame %d box [%g %g %g]\n",frame,
		pmset->per_bc->GetA(), pmset->per_bc->GetB(), pmset->per_bc->GetC());
#else
	PrintLog("Not implemented\n");
	return FALSE;
#endif
	return TRUE;
}

///////////////////////////////////////////////////////////////////////////////
#ifndef _MSC_VER
XTCTraj::XTCTraj()
{
	XTC_MAGIC=1995;
	fp=NULL;
	xdr=NULL;
	r=NULL;
}
XTCTraj::~XTCTraj()
{
	Close();
	if(r!=NULL)delete r;
}
int XTCTraj::XTC_CHECK(char *str,bool bResult)
{
	if (!bResult)
	{
		PrintLog("\nXTC error: read/write of %s failed\n",str);
		return 0;
	}
	return 1;
}
int XTCTraj::xdr_r2f_double(XDR *xdrs,double *r,bool bRead)
{
	float f;
	int   ret;
	
	if (!bRead)
		f = *r;
	ret = xdr_float(xdrs,&f);
	if (bRead)
		*r = f;
	
	return ret;
}
int XTCTraj::xdr_r2f_float(XDR *xdrs,float *r,bool bRead)
{
	float f;
	int   ret;
	ret = xdr_float(xdrs,&f);
	*r = f;
	return ret;
}
int XTCTraj::Open()
{
	xdr = new XDR;
	
	fp=fopen(FileName.c_str(), "rb");
	if(fp==NULL)
	{
		PrintLog("Can't find MD trajectory file %s\n",FileName.c_str() );
		return FALSE;
	}
	xdrstdio_create(xdr, fp, XDR_DECODE);
	return TRUE;
}
//ReadFirstFrame()
int XTCTraj::ReadNextFrame()
{
	xtc_header();
	check_xtc_magic();
	xtc_coord();
	return TRUE;
}
int XTCTraj::Close()
{
	if(xdr!=NULL)
	{
		xdr_destroy(xdr);
		delete xdr;
		xdr=NULL;
	}
	if(fp!=NULL)
	{
		fclose(fp);
		fp=NULL;
	}
	return TRUE;
}
int XTCTraj::check_xtc_magic()
{
	if (magic != XTC_MAGIC)
	{
		PrintLog("Magic Number Error in XTC file (read %d, should be %d)", magic, XTC_MAGIC);
		return FALSE;
	}
	return TRUE;
}
int XTCTraj::xtc_header()
{
	if (xdr_int(xdr,&magic) == 0)
	{
		PrintLog("xtc_header: can't read magic number\n");
		return FALSE;
	}
	if (xdr_int(xdr,&natoms) == 0)
	{
		PrintLog("xtc_header: can't read header: natoms\n");
		return FALSE;
	}
	if (xdr_int(xdr,&step) == 0)
	{
		PrintLog("xtc_header: can't read header: step\n");
		return FALSE;
	}
	if (xdr_r2f_float(xdr,&time,TRUE) == 0)
	{
		PrintLog("xtc_header: can't read header: time\n");
		return FALSE;
	}
	
	PrintLog("xtc_header: magic=%d natoms=%d step=%d time=%g\n", magic, natoms, step, time);
	
	
	return TRUE;
}
int XTCTraj::xtc_coord()
{
	if(r==NULL)
	{
		r = new HaVec_float(3*natoms);
	}
	//*r->begin()=0.0;
	//PrintLog("natoms=%d %d %g\n",natoms,r->size(),*r->begin());
	int result;
	int i,j;
	// box 
	result=1;
	for(i=0; ((i<3) && result); i++)
		for(j=0; ((j<3) && result); j++)
			result=XTC_CHECK("box",xdr_r2f_float(xdr,&(box[i][j]),TRUE));
	if (!result)
		return FALSE;
	
	float fprec;
	result=XTC_CHECK("x",xdr3dfcoord(r->begin(),&natoms,&fprec));
	
	//PrintLog("natoms=%d %d %g\n",natoms,r->size(),*r->begin());
	if (!result)
		return FALSE;
	return TRUE;
}
int XTCTraj::sizeofint(const int size) 
{
	unsigned int num = 1;
	int num_of_bits = 0;
    
	while (size >= num && num_of_bits < 32) {
		num_of_bits++;
		num <<= 1;
	}
	return num_of_bits;
}
int XTCTraj::sizeofints( const int num_of_ints, unsigned int sizes[]) {
	int i, num;
	unsigned int num_of_bytes, num_of_bits, bytes[32], bytecnt, tmp;
	num_of_bytes = 1;
	bytes[0] = 1;
	num_of_bits = 0;
	for (i=0; i < num_of_ints; i++) {	
		tmp = 0;
		for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++) {
			tmp = bytes[bytecnt] * sizes[i] + tmp;
			bytes[bytecnt] = tmp & 0xff;
			tmp >>= 8;
		}
		while (tmp != 0) {
			bytes[bytecnt++] = tmp & 0xff;
			tmp >>= 8;
		}
		num_of_bytes = bytecnt;
	}
	num = 1;
	num_of_bytes--;
	while (bytes[num_of_bytes] >= num) {
		num_of_bits++;
		num *= 2;
	}
	return num_of_bits + num_of_bytes * 8;

}
int XTCTraj::receivebits(int buf[], int num_of_bits) 
{

	int cnt, num; 
	unsigned int lastbits, lastbyte;
	unsigned char * cbuf;
	int mask = (1 << num_of_bits) -1;

	cbuf = ((unsigned char *)buf) + 3 * sizeof(*buf);
	cnt = buf[0];
	lastbits = (unsigned int) buf[1];
	lastbyte = (unsigned int) buf[2];
    
	num = 0;
	while (num_of_bits >= 8) {
		lastbyte = ( lastbyte << 8 ) | cbuf[cnt++];
		num |=  (lastbyte >> lastbits) << (num_of_bits - 8);
		num_of_bits -=8;
	}
	if (num_of_bits > 0) {
		if (lastbits < num_of_bits) {
			lastbits += 8;
			lastbyte = (lastbyte << 8) | cbuf[cnt++];
		}
		lastbits -= num_of_bits;
		num |= (lastbyte >> lastbits) & ((1 << num_of_bits) -1);
	}
	num &= mask;
	buf[0] = cnt;
	buf[1] = lastbits;
	buf[2] = lastbyte;
	return num; 
}
void XTCTraj::receiveints(int buf[], const int num_of_ints, int num_of_bits, unsigned int sizes[], int nums[])
{
	int bytes[32];
	int i, j, num_of_bytes, p, num;
	bytes[1] = bytes[2] = bytes[3] = 0;
	num_of_bytes = 0;
	while (num_of_bits > 8) {
		bytes[num_of_bytes++] = receivebits(buf, 8);
		num_of_bits -= 8;
	}
	if (num_of_bits > 0) {
		bytes[num_of_bytes++] = receivebits(buf, num_of_bits);
	}
	for (i = num_of_ints-1; i > 0; i--) {
		num = 0;
		for (j = num_of_bytes-1; j >=0; j--) {
			num = (num << 8) | bytes[j];
			p = num / sizes[i];
			bytes[j] = p;
			num = num - p * sizes[i];
		}
		nums[i] = num;
	}
	nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
}
/*____________________________________________________________________________
 |
 | xdr3dfcoord - read or write compressed 3d coordinates to xdr file.
 |
 | this routine reads or writes (depending on how you opened the file with
 | xdropen() ) a large number of 3d coordinates (stored in *fp).
 | The number of coordinates triplets to write is given by *size. On
 | read this number may be zero, in which case it reads as many as were written
 | or it may specify the number if triplets to read (which should match the
 | number written).
 | Compression is achieved by first converting all floating numbers to integer
 | using multiplication by *precision and rounding to the nearest integer.
 | Then the minimum and maximum value are calculated to determine the range.
 | The limited range of integers so found, is used to compress the coordinates.
 | In addition the differences between succesive coordinates is calculated.
 | If the difference happens to be 'small' then only the difference is saved,
 | compressing the data even more. The notion of 'small' is changed dynamically
 | and is enlarged or reduced whenever needed or possible.
 | Extra compression is achieved in the case of GROMOS and coordinates of
 | water molecules. GROMOS first writes out the Oxygen position, followed by
 | the two hydrogens. In order to make the differences smaller (and thereby
 | compression the data better) the order is changed into first one hydrogen
 | then the oxygen, followed by the other hydrogen. This is rather special, but
 | it shouldn't harm in the general case.
 |
 */
#define LASTIDX (sizeof(magicints) / sizeof(*magicints))
#define FIRSTIDX 9
static int magicints[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0,
 8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
 80, 101, 128, 161, 203, 256, 322, 406, 512, 645,
 812, 1024, 1290, 1625, 2048, 2580, 3250, 4096, 5060, 6501,
 8192, 10321, 13003, 16384, 20642, 26007, 32768, 41285, 52015, 65536,
 82570, 104031, 131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561,
 832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042,
 8388607, 10568983, 13316085, 16777216 };
int XTCTraj::xdr3dfcoord(float *fp, int *size, float *precision)
{
	static int *ip = NULL;
	static int oldsize;
	static int *buf;

	int minint[3], maxint[3], mindiff, *lip, diff;
	int lint1, lint2, lint3, oldlint1, oldlint2, oldlint3, smallidx;
	int minidx, maxidx;
	unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3, *luip;
	int flag, k;
	int smallnum, smaller, larger, i, is_small, is_smaller, run, prevrun;
	float *lfp, lf;
	int tmp, *thiscoord,  prevcoord[3];
	unsigned int tmpcoord[30];

	int bufsize, lsize;
	unsigned int bitsize;
	float inv_precision;
	int errval = 1;

	bitsizeint[0] = bitsizeint[1] = bitsizeint[2] = 0;
	prevcoord[0]  = prevcoord[1]  = prevcoord[2]  = 0;


	
		// xdr is open for reading
	
		if (xdr_int(xdr, &lsize) == 0) 
			return 0;
		if (*size != 0 && lsize != *size) {
			fprintf(stderr, "wrong number of coordinates in xdr3dfcoord; "
					"%d arg vs %d in file", *size, lsize);
		}
		*size = lsize;
		size3 = *size * 3;
		if (*size <= 9) {
			*precision = -1;
			return (xdr_vector(xdr, (char *) fp, (unsigned int)size3, 
							(unsigned int)sizeof(*fp), (xdrproc_t)xdr_float));
		}
		xdr_float(xdr, precision);
		if (ip == NULL) {
			ip = (int *)malloc((size_t)(size3 * sizeof(*ip)));
			if (ip == NULL) {
				fprintf(stderr,"malloc failed\n");
				exit(1);
			}
			bufsize = size3 * 1.2;
			buf = (int *)malloc((size_t)(bufsize * sizeof(*buf)));
			if (buf == NULL) {
				fprintf(stderr,"malloc failed\n");
				exit(1);
			}
			oldsize = *size;
		} else if (*size > oldsize) {
			ip = (int *)realloc(ip, (size_t)(size3 * sizeof(*ip)));
			if (ip == NULL) {
				fprintf(stderr,"malloc failed\n");
				exit(1);
			}
			bufsize = size3 * 1.2;
			buf = (int *)realloc(buf, (size_t)(bufsize * sizeof(*buf)));
			if (buf == NULL) {
				fprintf(stderr,"malloc failed\n");
				exit(1);
			}
			oldsize = *size;
		}
		buf[0] = buf[1] = buf[2] = 0;
	
		xdr_int(xdr, &(minint[0]));
		xdr_int(xdr, &(minint[1]));
		xdr_int(xdr, &(minint[2]));

		xdr_int(xdr, &(maxint[0]));
		xdr_int(xdr, &(maxint[1]));
		xdr_int(xdr, &(maxint[2]));
		
		sizeint[0] = maxint[0] - minint[0]+1;
		sizeint[1] = maxint[1] - minint[1]+1;
		sizeint[2] = maxint[2] - minint[2]+1;
	
		/* check if one of the sizes is to big to be multiplied */
		if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff) {
			bitsizeint[0] = sizeofint(sizeint[0]);
			bitsizeint[1] = sizeofint(sizeint[1]);
			bitsizeint[2] = sizeofint(sizeint[2]);
			bitsize = 0; /* flag the use of large sizes */
		} else {
			bitsize = sizeofints(3, sizeint);
		}
	
		if (xdr_int(xdr, &smallidx) == 0)	
			return 0;
		maxidx = MIN(LASTIDX, smallidx + 8) ;
		minidx = maxidx - 8; /* often this equal smallidx */
		smaller = magicints[MAX(FIRSTIDX, smallidx-1)] / 2;
		smallnum = magicints[smallidx] / 2;
		sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
		larger = magicints[maxidx];

		/* buf[0] holds the length in bytes */

		if (xdr_int(xdr, &(buf[0])) == 0)
			return 0;
		if (xdr_opaque(xdr, (char *)&(buf[3]), (unsigned int)buf[0]) == 0)
			return 0;
		buf[0] = buf[1] = buf[2] = 0;
	
		lfp = fp;
		inv_precision = 1.0 / * precision;
		run = 0;
		i = 0;
		lip = ip;
		while ( i < lsize ) {
			thiscoord = (int *)(lip) + i * 3;

			if (bitsize == 0) {
				thiscoord[0] = receivebits(buf, bitsizeint[0]);
				thiscoord[1] = receivebits(buf, bitsizeint[1]);
				thiscoord[2] = receivebits(buf, bitsizeint[2]);
			} else {
				receiveints(buf, 3, bitsize, sizeint, thiscoord);
			}
	    
			i++;
			thiscoord[0] += minint[0];
			thiscoord[1] += minint[1];
			thiscoord[2] += minint[2];
	    
			prevcoord[0] = thiscoord[0];
			prevcoord[1] = thiscoord[1];
			prevcoord[2] = thiscoord[2];
	    
	   
			flag = receivebits(buf, 1);
			is_smaller = 0;
			if (flag == 1) {
				run = receivebits(buf, 5);
				is_smaller = run % 3;
				run -= is_smaller;
				is_smaller--;
			}
			if (run > 0) {
				thiscoord += 3;
				for (k = 0; k < run; k+=3) {
					receiveints(buf, 3, smallidx, sizesmall, thiscoord);
					i++;
					thiscoord[0] += prevcoord[0] - smallnum;
					thiscoord[1] += prevcoord[1] - smallnum;
					thiscoord[2] += prevcoord[2] - smallnum;
					if (k == 0) {
			/* interchange first with second atom for better
						* compression of water molecules
			*/
						tmp = thiscoord[0]; thiscoord[0] = prevcoord[0];
						prevcoord[0] = tmp;
						tmp = thiscoord[1]; thiscoord[1] = prevcoord[1];
						prevcoord[1] = tmp;
						tmp = thiscoord[2]; thiscoord[2] = prevcoord[2];
						prevcoord[2] = tmp;
						*lfp++ = prevcoord[0] * inv_precision;
						*lfp++ = prevcoord[1] * inv_precision;
						*lfp++ = prevcoord[2] * inv_precision;
					} else {
						prevcoord[0] = thiscoord[0];
						prevcoord[1] = thiscoord[1];
						prevcoord[2] = thiscoord[2];
					}
					*lfp++ = thiscoord[0] * inv_precision;
					*lfp++ = thiscoord[1] * inv_precision;
					*lfp++ = thiscoord[2] * inv_precision;
				}
			} else {
				*lfp++ = thiscoord[0] * inv_precision;
				*lfp++ = thiscoord[1] * inv_precision;
				*lfp++ = thiscoord[2] * inv_precision;		
			}
			smallidx += is_smaller;
			if (is_smaller < 0) {
				smallnum = smaller;
				if (smallidx > FIRSTIDX) {
					smaller = magicints[smallidx - 1] /2;
				} else {
					smaller = 0;
				}
			} else if (is_smaller > 0) {
				smaller = smallnum;
				smallnum = magicints[smallidx] / 2;
			}
			sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
		}
	return 1;
}
#endif
