/*! \file trajectoryanalysis.cpp
 
    Classes for analysis of simulation trajectories

    \author Igor Kurnikov
    \date 2009

*/

#include <time.h>
#include "hamolset.h"
#include "hacoord.h"
#include "trajanal.h"
#include "hasimulator.h"
#include "hamolview.h"
#include "hamolecule.h"

using namespace harlem;

TrajPointInfo::TrajPointInfo()
{
	Clear();
}

TrajPointInfo::~TrajPointInfo()
{

}

void TrajPointInfo::Clear()
{
	do_skip = FALSE;
	ipt = 0;
	pcrd = NULL;   
	is_accepted = TRUE; 
	pcrd_rejected = NULL;
	double tot_energy = 0.0;      
	double energy_rejected = 0.0;
}


TrajIOAgent::TrajIOAgent(MCSimulator* p_sim_new)
{
	traj_ene_file_name = "ene_file.dat";
	traj_file_name     = "traj_file.dat"; 
	traj_all_pts_file_name = "all_mc_pts.dat";
	save_image_seq_gif  = FALSE;
	save_image_seq_pict = FALSE;
	output_rejected_points = FALSE; 
	p_sim = p_sim_new;

	npt = 0;
	average_ene = 0.0;

	active_flag = TRUE;
}

TrajIOAgent::~TrajIOAgent()
{

}

int TrajIOAgent::IsActive() const
{
	return active_flag;
}

void TrajIOAgent::SetActive(int active_flag_new)
{
	active_flag = active_flag_new;
}

int TrajIOAgent::Init(TrajPointInfo* ppt_info)
{
	if( IsReadCoord() )
	{
		traj_file.open(traj_file_name.c_str(), ios::in);
	}
	else if( IsWriteCoord() )
	{
		traj_file.open(traj_file_name.c_str(), ios::out);
		if(output_rejected_points) all_points_file.open(traj_all_pts_file_name.c_str(), ios::out);
	}
	if( traj_file.fail())
	{
		PrintLog("Error in TrajIOAgent::Init() \n Can't open file with trajectory \n");
		return FALSE;
	}

	if(IsReadEnergy())
	{
		ene_file.open(traj_ene_file_name.c_str(), ios::in);	
	}
	else if( IsWriteEnergy())
	{
		ene_file.open(traj_ene_file_name.c_str(), ios::out);		
	}

	npt = 0;
	average_ene = 0.0;

//	if( p_im_mod->freeze_first_mol )
//	{
//		for(i = 0; i < 6; i++)
//		{
//			p_crd->FreezeCrd(i);
//		}
//      HaVec_double crd_v = p_crd->AsVecDouble();
//		sprintf(buf,"%6d  %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f ",
//					-1, crd_v[0], crd_v[1], crd_v[2], crd_v[3], crd_v[4], crd_v[5]);
//		traj_file << buf << endl;
//	}
	return TRUE;
}


int TrajIOAgent::AnalyzePt(TrajPointInfo* ppt_info)
{
	char buf[256];
	
	MolSet* pmset = p_sim->GetMolSet();
	if( pmset) pmset->info_str.clear();

	if( IsReadCoord() )
	{
		traj_file >> ppt_info->ipt;
		LoadCrdOptions ropt;
		ropt.SetLoadAllCrd(true);
		p_sim->p_crd->LoadFromStream(traj_file, &ropt );
		ppt_info->pcrd = p_sim->p_crd;
		if( !ppt_info->do_skip )
		{
			p_sim->SetCoord(p_sim->p_crd);
		}
	}

	if( IsWriteCoord() )
	{
		traj_file << ppt_info->ipt;
		if(ppt_info->pcrd != NULL ) ppt_info->pcrd->SaveToStream(traj_file);

		if(output_rejected_points )
		{
			all_points_file << ppt_info->ipt;
			if( !(ppt_info->is_accepted) && ppt_info->pcrd_rejected != NULL )
				ppt_info->pcrd_rejected->SaveToStream(all_points_file);
			if( ppt_info->is_accepted && ppt_info->pcrd != NULL )
			{
				ppt_info->pcrd->SaveToStream(all_points_file);
				all_points_file << " accepted " << endl;
			}
		}
	}

	int ipt_ene;

	if(IsReadEnergy())
	{
		char ene_str[256];
		ene_file.getline(ene_str,255);
 		istrstream is_ene(ene_str);
		is_ene >> ipt_ene;
		is_ene >> ppt_info->tot_energy;
//		ppt_info->tot_energy = p_im_mod->cur_intermol_ene;
//		if( p_im_mod->calc_et_rate )
//		{
//			is_ene >> p_im_mod->add_eff_ene;
//			is_ene >> p_im_mod->add_eff_ene;
//		}
//		else
//		{
//			p_im_mod->add_eff_ene = 0.0;
//		}
	}

	if( IsWriteEnergy() && !ppt_info->do_skip )
	{
		npt++;
		average_ene += ppt_info->tot_energy; 
		ene_file << ppt_info->ipt << "  " << ppt_info->tot_energy;
		if( npt > 0 && (npt%100 == 0) ) PrintLog(" average energy: %12.6e \n", average_ene/npt );
	}

	if( !ppt_info->do_skip )
	{
		sprintf(buf,"point number %d",ppt_info->ipt );

		if( pmset) pmset->info_str.push_back(buf);
				
//		if(p_im_mod->calc_et_rate)
//		{
//			double real_intermol_ene = p_im_mod->cur_intermol_ene - p_im_mod->add_eff_ene;
//			if(!p_sim->dont_calc_ene_flag)
//			{
//				ene_file << "  " << real_intermol_ene << " " << p_im_mod->add_eff_ene;
//			}
//			sprintf(buf," intermolecular energy = %12.6f kcal/mol",real_intermol_ene);
//			pmset->info_str.push_back(buf);
//			sprintf(buf," ET effective energy = %12.6f kcal/mol", p_im_mod->add_eff_ene);
//			pmset->info_str.push_back(buf);
//		}
	}
	if((save_image_seq_gif || save_image_seq_pict) && !ppt_info->do_skip)
	{
		pmset->RefreshAllViews(RFApply | RFRefresh);
		HaMolView* pview = pmset->GetActiveMolView();
		if(save_image_seq_gif)
		{
			sprintf(buf,"traj%6d.gif",(100000 + ppt_info->ipt));
			pview->WriteGIFFile(buf);
		}
		if(save_image_seq_pict)
		{
			sprintf(buf,"traj%6d.pic",(100000 + ppt_info->ipt));
			pview->WritePICTFile(buf);
		}
	}
	
	if(!p_sim->dont_calc_ene_flag && !ppt_info->do_skip)
	{
		ene_file << endl;
	}
	return TRUE;
}

int TrajIOAgent::Finalize()
{
	if(ene_file.is_open()) ene_file.close();
	if((traj_io_mode & COORD_WRITE) && traj_file.is_open()) traj_file.close();
	if(all_points_file.is_open()) all_points_file.close();

	if( npt > 0 ) PrintLog(" average energy: %12.6e \n", average_ene/npt );
	return TRUE;
}

void  TrajIOAgent::SetReadCoord(int set_on)
{
	if( set_on )
	{
		traj_io_mode |= COORD_READ;
	}
	else
	{
		traj_io_mode &= (~COORD_READ);
	}
}

void  TrajIOAgent::SetWriteCoord(int set_on)
{
	if( set_on )
	{
		traj_io_mode |= COORD_WRITE;
	}
	else
	{
		traj_io_mode &= (~COORD_WRITE);
	}
}

void  TrajIOAgent::SetReadEnergy (int set_on)
{
	if( set_on )
	{
		traj_io_mode |= ENERGY_READ;
	}
	else
	{
		traj_io_mode &= (~ENERGY_READ);
	}
}

void  TrajIOAgent::SetWriteEnergy(int set_on)
{
	if( set_on )
	{
		traj_io_mode |= ENERGY_WRITE;
	}
	else
	{
		traj_io_mode &= (~ENERGY_WRITE);
	}
}

int TrajIOAgent::IsReadCoord()   const
{
	return (traj_io_mode & COORD_READ);
}

int TrajIOAgent::IsWriteCoord()  const
{
	return (traj_io_mode & COORD_WRITE);
}

int TrajIOAgent::IsReadEnergy()  const
{
	return (traj_io_mode & ENERGY_READ);
}

int TrajIOAgent::IsWriteEnergy() const
{
	return (traj_io_mode & ENERGY_WRITE);
}

TraceMolAgent::TraceMolAgent(MolSet* pmset_new)
{
	pmset = pmset_new;

	trace_mol = NULL;
	trace_chain = NULL;
	int itr_res = 0;

}

TraceMolAgent::~TraceMolAgent()
{

}

int TraceMolAgent::IsActive() const
{
	return active_flag;
}

void TraceMolAgent::SetActive(int active_flag_new)
{
	active_flag = active_flag_new;
}


int TraceMolAgent::Init(TrajPointInfo* ppt_info)
{
	trace_mol = pmset->GetMolByName("TRACE_MOL");
	if(trace_mol != NULL)
	{
		pmset->DeleteMol(trace_mol);
	}
	trace_mol = pmset->CreateMolecule();
	trace_mol->SetObjName("TRACE_MOL");
	trace_chain = trace_mol->AddChain(' ');
	trace_res   = NULL;

	return TRUE;
}

int TraceMolAgent::AnalyzePt(TrajPointInfo* ppt_info)
{
	itr_res++;
	HaResidue* prev_res = trace_res;
	trace_res = trace_chain->AddResidue(itr_res);
	AtomIteratorAtomGroup aitr(&traced_atoms);
	HaAtom* aptr_tr;
	HaAtom* aptr;
	int iat = 0;
	for(aptr_tr = aitr.GetFirstAtom(); aptr_tr; aptr_tr = aitr.GetNextAtom())
	{
		iat++;
		aptr = trace_res->AddNewAtom();
		aptr->SetElemNo(5+iat);
		char atn[10];
		sprintf(atn,"X%1i",iat);
		aptr->SetName(atn);
		aptr->SetX(aptr_tr->GetX());
		aptr->SetY(aptr_tr->GetY());
		aptr->SetZ(aptr_tr->GetZ());

		if(prev_res != NULL)
		{
			HaAtom* prev_aptr = prev_res->GetAtomByName(atn);
			if(prev_aptr != NULL && aptr != NULL)
			{
				//								pMol->AddBond(prev_aptr, aptr);
			}
		}
	}
	return TRUE;
}

int TraceMolAgent::Finalize()
{
	return TRUE;
}

UpdateMolViewNotifyAgent::UpdateMolViewNotifyAgent(MolSet* pmset_new)
{
	pmset = pmset_new;
	is_moved = FALSE;
	update_interval = 1.0;
	next_update_time = 0;

	active_flag = TRUE;
}

UpdateMolViewNotifyAgent::~UpdateMolViewNotifyAgent()
{

}

int UpdateMolViewNotifyAgent::IsActive() const
{
	return active_flag;
}

void UpdateMolViewNotifyAgent::SetActive(int active_flag_new)
{
	active_flag = active_flag_new;
}

int UpdateMolViewNotifyAgent::Init(TrajPointInfo* ppt_info)
{
	next_update_time = clock() + (clock_t)update_interval*CLOCKS_PER_SEC;
	is_moved = FALSE;
	return TRUE;
}

int UpdateMolViewNotifyAgent::AnalyzePt(TrajPointInfo* ppt_info)
{
	if( ppt_info->is_accepted ) 
	{
		is_moved = TRUE;
	}
	if( is_moved )
	{
		if( clock() > next_update_time )
		{
			pmset->AnnounceGeomChange();
			is_moved = FALSE;
			next_update_time = clock() + (clock_t)update_interval*CLOCKS_PER_SEC;
		}
	}
	return TRUE;
}

int UpdateMolViewNotifyAgent::Finalize()
{
	pmset->AnnounceGeomChange();
	return TRUE;
}

