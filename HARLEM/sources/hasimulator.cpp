/*! \file hasimulator.cpp
 
    Basic Monte Carlo Simulator Module  

    \author Igor Kurnikov
    \date 1999-2009

*/

#include <memory>
#include "wx/thread.h"

#include "hamolset.h"
#include "hacoord.h"
#include "haenefunc.h"
#include "trajanal.h"
#include "hasimulator.h"

#include "randomip.h"

using namespace harlem;

MCSimulator::MCSimulator()
{
	num_mc_steps = 100;        
	stop_calc_flag = FALSE;
	
	p_crd = NULL; 

	delay_time = 0.0;           
	npt_begin = 0;         
	npt_step = 1; 
	npt_end = 9999999999;
	
	dont_calc_ene_flag = FALSE;  

	sim_thread = NULL;
	SetTemperature(300.0);
	p_ene_func = NULL; 
	p_rand_num_gen = NULL;
	pmset = NULL;
}
	
MCSimulator::~MCSimulator()
{

}

class SimulatorThread: public wxThread
{
public:
	SimulatorThread(MCSimulator* ptr_sim_new) 
	{ 
		ptr_sim = ptr_sim_new; 
	}

	virtual ExitCode Entry()
	{
		ptr_sim->RunMC();
		return 0;
	}
	MCSimulator* ptr_sim;
};

int MCSimulator::RunMCThread()
{
	sim_thread = new SimulatorThread(this);
	sim_thread->Create();
	sim_thread->Run();

	return TRUE;
}

int MCSimulator::SetEnergyFunc( HaEnergyFunc* p_ene_func_new)
{
	p_ene_func = p_ene_func_new;
	return TRUE;
}

int MCSimulator::InitEnergyFunc()
{
	return TRUE;
}

double MCSimulator::ComputeEnergy(harlem::Coord* pcrd)
{
	if( p_ene_func == NULL) return 0.0;
	return p_ene_func->ComputeEnergy(pcrd);
}

int MCSimulator::SetInitPoint(harlem::Coord* pcrd_new)
{
	p_crd = pcrd_new;
	return TRUE;
}

void MCSimulator::SetTemperature(double temp_new)
{
	double kT_MC = temp_new*BOLTZ * AVOG_NUM/(1000.0*CAL_TO_JOULES); 
}

double MCSimulator::GetTemperature() const
{
	return (kT_MC*1000.0*CAL_TO_JOULES)/(BOLTZ * AVOG_NUM);
}

int MCSimulator::RunMC() 
{	
	InitEnergyFunc();
	
	if( p_crd == NULL) SetInitPoint();
	
	int np_reject = 0;
	stop_calc_flag = FALSE;

	std::auto_ptr<harlem::Coord> p_crd_old( p_crd->clone() );  // to ensure deleting of pcrd upon exit from the function; 

	int istep = 0;
	double ene_cur = ComputeEnergy(p_crd);
	double ene_prev = ene_cur; 

	TrajPointInfo traj_pt_info;
	traj_pt_info.ipt = 0;
	traj_pt_info.pcrd = p_crd;
	traj_pt_info.is_accepted = NULL;

	InitTrajAnalysis(&traj_pt_info);

	try
	{
		for( istep = 0; istep < num_mc_steps ; istep++ )
		{
			if( stop_calc_flag )
			{
				stop_calc_flag = FALSE;
				break;
			}
			p_crd_old->SetFrom(p_crd);
			IncrementCrd(p_crd);
			ene_cur = ComputeEnergy(p_crd);			
			
			bool is_accept = false;	
			
			if( ene_cur < ene_prev ) 
			{
				is_accept = true;
			}
			else
			{
				double rand_dbl = p_rand_num_gen->GetRndNum(); 
				double dtest = exp( (ene_prev - ene_cur)/kT_MC );
				if( dtest > rand_dbl )
				{
					is_accept = true;
				}
			} 
			if( is_accept )
			{
				ene_prev = ene_cur;
			}

			traj_pt_info.Clear();
			traj_pt_info.ipt = istep+1;
			if( is_accept)
			{
				traj_pt_info.is_accepted = TRUE;
				traj_pt_info.pcrd = p_crd;
			}
			else
			{
				traj_pt_info.is_accepted = FALSE;
				traj_pt_info.pcrd          = p_crd_old.get();
				traj_pt_info.pcrd_rejected = p_crd; 
			}

			ComputePropTrajPoint(&traj_pt_info);

			if( !is_accept )
			{
				np_reject++;
				p_crd->SetFrom( p_crd_old.get() );
				SetCoord(p_crd);
			}
		}
	}
	catch(const char* msg)
	{
		PrintLog(msg);
	}

	FinalizeTrajAnalysis();

	if(istep > 0) PrintLog(" point acceptance ratio: %7.3f \n",    (istep - np_reject)/((double)(istep)) );  
	
	return True;
}

int MCSimulator::PauseMC()
{
	return TRUE;
}

int MCSimulator::ResumeMC()
{
	return TRUE;
}

int MCSimulator::StopMC()
{
	stop_calc_flag = TRUE;
	return TRUE;
}

int MCSimulator::AnalyzeTrajectory()
{
	TrajIOAgent* p_traj_io_ag = GetTrajectoryIOAgent();
	p_traj_io_ag->SetReadCoord(TRUE);
	p_traj_io_ag->SetWriteCoord(FALSE);
	
	TrajPointInfo pt_info;
	try
	{
		InitTrajAnalysis(&pt_info);
		int i;
		if(npt_begin > 0)
		{
			pt_info.do_skip = TRUE;
			ComputePropTrajPoint(&pt_info);
		}

		for(;;) // Start cycle on trajectory points 
		{
			pt_info.do_skip = FALSE;
			ComputePropTrajPoint(&pt_info);
			if(delay_time > 0)
			{
				wxThread::Sleep(delay_time);
			}	
			if(npt_step > 1)
			{
				for( i = 0; i < (npt_step - 1); i++)
				{
					pt_info.do_skip = TRUE;
					ComputePropTrajPoint(&pt_info);
				}
			}
		} // End cycle on trajectory points 		
	}
	catch( char* err_msg )
	{
		PrintLog(err_msg);
	}
	FinalizeTrajAnalysis();

	return TRUE;
}


int MCSimulator::InitTrajAnalysis(TrajPointInfo* ppt_info)
{
	int na = agents.size();
	int i;
	for( i = 0; i < na; i++)
	{
		TrajAnalAgent* p_ag = agents[i];
		p_ag->Init();
	}
	return TRUE;
}

int MCSimulator::ComputePropTrajPoint(TrajPointInfo* ppt_info )
{
	int na = agents.size();
	int i;
	for( i = 0; i < na; i++)
	{
		TrajAnalAgent* p_ag = agents[i];
		p_ag->AnalyzePt(ppt_info);
	}
	return TRUE;
}

int MCSimulator::FinalizeTrajAnalysis()
{
	int na = agents.size();
	int i;
	for( i = 0; i < na; i++)
	{
		TrajAnalAgent* p_ag = agents[i];
		p_ag->Finalize();
	}
	return TRUE;
}

int MCSimulator::DeleteTrajAnalAgent(TrajAnalAgent* p_ag)
{
	int i;
	int nag = agents.size();
	for( i= 0; i < nag; i++)
	{
		if( agents[i] == p_ag) 
		{
			agents.erase( agents.begin() + i);
			return TRUE;
		}
	}
	return FALSE;
}

TrajIOAgent* MCSimulator::GetTrajectoryIOAgent()
{
	int na = agents.size();
	int i;
	for( i = 0; i < na; i++)
	{
		TrajAnalAgent* p_ag = agents[i];
		if( p_ag->GetClassName() == "TrajIOAgent") return (TrajIOAgent*)p_ag;
	}
	TrajIOAgent* p_ene_ag = new TrajIOAgent(this);
	agents.push_back(p_ene_ag);
	return p_ene_ag;
}

TraceMolAgent* MCSimulator::GetTrajectoryTraceAgent(int create_agent)
{
	int na = agents.size();
	int i;
	for( i = 0; i < na; i++)
	{
		TrajAnalAgent* p_ag = agents[i];
		if( p_ag->GetClassName() == "TraceMolAgent") return (TraceMolAgent*)p_ag;
	}
	if( create_agent )
	{
		TraceMolAgent* p_tr_mol_ag = new TraceMolAgent(this->GetMolSet());
		agents.push_back(p_tr_mol_ag);
		return p_tr_mol_ag;
	}
	return NULL;
}

UpdateMolViewNotifyAgent*  MCSimulator::GetMolViewNotifyAgent(int create_agent)
{
	int na = agents.size();
	int i;
	for( i = 0; i < na; i++)
	{
		TrajAnalAgent* p_ag = agents[i];
		if( p_ag->GetClassName() == "UpdateMolViewNotifyAgent") return (UpdateMolViewNotifyAgent*)p_ag;
	}
	if( create_agent )
	{
		UpdateMolViewNotifyAgent* p_tr_mol_ag = new UpdateMolViewNotifyAgent(this->GetMolSet());
		agents.push_back(p_tr_mol_ag);
		return p_tr_mol_ag;
	}
	return NULL;

}


MolSet* MCSimulator::GetMolSet()
{
	return pmset;
}


