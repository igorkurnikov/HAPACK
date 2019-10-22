/*! \file mm_traj_anal.cpp

    Classes to Analyze MD trajectories 
 
    \author Igor Kurnikov 
    \date 2010-
*/

#include <mpi.h>

#include <math.h>
#include <stdexcept>
#include <sstream>

#if !defined(HARLEM_PYTHON_NO)
#include "Python.h"
#endif

#include "haio.h"

#include <wx/process.h>
#include <wx/thread.h>

#include <boost/algorithm/string.hpp>

#include "FCMangle.h"

#include "harlemapp.h"
#include "vec3d.h"
#include "haatom.h"
#include "haatgroup.h"
#include "hamolset.h"
#include "moleditor.h"
#include "mm_traj_anal.h"
#include "hamolmech.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "mm_driver_amber.h"

extern "C" 
{
//extern void runfiles_mod_mp_mdwrit_(int* imin, int* ntwr, int* ntxo, int* ntb_par, int* nstep, int* atm_cnt, double* crd, double* box, double* vel, double* tt);
//extern void pbc_mod_mp_wrap_molecules_(int* nspm, int* nsp, double* crd);
//extern void pbc_mod_mp_wrap_to_(int* nspm, int* nsp, double* crd, double* box);
//extern void runfiles_mod_mp_corpac_(int* iend, double* crd, int* istart, int* nf);
extern void FC_FUNC_MODULE(runfiles_mod,mdwrit)(int* imin, int* ntwr, int* ntxo, int* ntb_par, int* nstep, int* atm_cnt, double* crd, double* box, double* vel, double* tt);
extern void FC_FUNC_MODULE(pbc_mod,wrap_molecules)(int* nspm, int* nsp, double* crd);
extern void FC_FUNC_MODULE(pbc_mod,wrap_to)(int* nspm, int* nsp, double* crd, double* box);
extern void FC_FUNC_MODULE(runfiles_mod,corpac)(int* iend, double* crd, int* istart, int* nf);
//extern void bintraj_mod_mp_write_binary_cell_dat_(double* box);
//extern void runfiles_mod_mp_mdeng_(int* nstep, double* time, double* sys_info, double* box);
extern void FC_FUNC_MODULE(runfiles_mod,mdeng)(int* nstep, double* time, double* sys_info, double* box);
// extern void bintraj_mod_mp_end_binary_frame_(int* nf);
//extern void nmr_calls_mod_mp_ndvptx_(double* crd, double* frc, double* mass, int* mdout);
//extern void gb_ene_mod_mp_print_born_radii_stat_(double* tspan);
extern void FC_FUNC_MODULE(nmr_calls_mod,ndvptxs)(double* crd, double* frc, double* mass, int* mdout);
extern void FC_FUNC_MODULE(gb_ene_mod,print_born_radii_stat)(double* tspan);
};

MDTrajAnalMod::MDTrajAnalMod(HaMolMechMod* p_mm_mod_new)
{
	p_mm_mod = p_mm_mod_new;
	p_mm_model = p_mm_mod->p_mm_model;
	p_amber_driver = p_mm_mod->p_amber_driver;
	pmset = p_mm_mod->GetMolSet();

	delay_time = 0.0;
    SetPtBegin(1);
	SetPtStep(1);
	SetPtEnd(999999999);
	ipt_curr  = 0;
	SetAlignToFstPt(false);
	SetAlignToCurrentCrd(false);
	SetReadPBox(PBOX_READ_IF_SET);
	SetWritePBox(PBOX_WRITE_IF_SET);
	SetWrapCrd(FALSE);
}

MDTrajAnalMod::~MDTrajAnalMod()
{

}

class AnalyzeMDTrajThread: public wxThread
{
public:
	AnalyzeMDTrajThread(HaMolMechMod* ptr_mm_mod_new): wxThread() 
	{
		PrintLog("AnalyzeMDTrajThread::MolMechThread() \n");
		ptr_mm_mod = ptr_mm_mod_new; 
	}
	virtual void* Entry()
	{
		PrintLog("AnalyzeMDTrajThread::Entry() \n");
		ptr_mm_mod->p_traj_anal_mod->AnalyzeTrajectoryInternal();
		return NULL;
	}
	virtual void OnExit()
	{
		ptr_mm_mod->run_thread = NULL;
		PrintLog("AnalyzeMDTrajThread::OnExit() \n");
	}

public:
	HaMolMechMod* ptr_mm_mod;
};

//class AnalMDCtrlThread: public wxThread
////!< MD analysis controlling thread
//{
//public:
//	AnalMDCtrlThread(HaMolMechMod* ptr_mm_mod_new): wxThread() 
//	{
//		PrintLog("AnalMDCtrlThread::AnalMDCtrlThread() \n");
//		ptr_mm_mod = ptr_mm_mod_new; 
//	}
//	virtual void* Entry()
//	{
//		ptr_mm_mod->p_traj_anal_mod->ControlCalc();
//		return NULL;
//	}
//	virtual void OnExit()
//	{
//		ptr_mm_mod->ctrl_thread = NULL;
//		PrintLog("AnalMDCtrlThread::OnExit() \n");
//	}
//
//public:
//	HaMolMechMod* ptr_mm_mod;
//};


int MDTrajAnalMod::AnalyzeTrajectory(int sync)
{
	int ires = TRUE;
	if( p_mm_mod->to_init_simulations || p_mm_mod->p_mm_model->to_init_mm_model ) ires = p_mm_mod->InitMMSimulations();
	if(!ires) return TRUE;

	if( !sync )
	{
		p_mm_mod->run_thread = new AnalyzeMDTrajThread(p_mm_mod);
		p_mm_mod->run_thread->Create();
		p_mm_mod->run_thread->Run();
	}
	else
	{
		AnalyzeTrajectoryInternal();
		return TRUE;
	}
	return TRUE;
}

int MDTrajAnalMod::AnalyzeTrajectoryInternal()
{
	p_mm_mod->to_stop_simulations = FALSE;
	bool files_opened = false;

	bool use_index = false;
	if( pt_pos.size() > 0 ) use_index = true;

	HaVec_double crd_init;
	pmset->SaveCrdToArray(crd_init);

	try
	{
		int ires = OpenAmberTrajFilesToRead();
		if(!ires) throw std::runtime_error("Error to open MD trajectory"); 
		files_opened = true;

		int i,j;
		for(i = 0; i < agents.size(); i++)
		{
			if(agents[i]->IsActive()) agents[i]->Init(NULL);
		}

		if(!traj_script.empty())
		{
			pApp->ExecuteScriptInString("script_status = SCRIPT_START");
			pApp->ExecuteScriptFromFile(traj_script.c_str());
		}

		clock_t update_time = clock();

		


		ipt_curr = 1;
		if( npt_begin > 1 ) 
		{
			if( use_index )
			{
				if( ipt_curr > pt_pos.size() ) throw std::runtime_error("index of initial point is larger than the number of points in the trajectory");
				ipt_curr = npt_begin;
			}
			else
			{
				for( int j = 0; j < (npt_begin - 1); j++)
				{
					ires = ReadTrajPoint();
					if( !ires ) throw std::runtime_error("Reading Error Skipping initial MD trajectory points");
					ipt_curr++;
				}
			}
		}
		
		for(;;)
		{
			if( use_index )
			{
				if( ipt_curr <= pt_pos.size() && ipt_curr <= npt_end ) ires = LoadCurrPt();
			}
			else
			{
				ires = ReadTrajPoint();
			}
			if( !ires ) p_mm_mod->to_stop_simulations = TRUE;
			
			if( ipt_curr > npt_end || ( use_index && (ipt_curr > pt_pos.size()) ) ) p_mm_mod->to_stop_simulations = TRUE;
			if( p_mm_mod->to_stop_simulations ) break;

			TrajPointInfo pt_info;
			pt_info.ipt = ipt_curr; 

			if(delay_time > 0.1) wxThread::Sleep( delay_time*1000 );

			for(i = 0; i < agents.size(); i++)
			{
				if(agents[i]->IsActive()) agents[i]->AnalyzePt(&pt_info);
			}
			if(!traj_script.empty())
			{
				std::string cmd_line = "script_status = SCRIPT_CONTINUE \n";
				cmd_line += "idx_curr_pt = ";
				cmd_line += harlem::ToString(ipt_curr);
				cmd_line += "\nipt_curr = ";
				cmd_line += harlem::ToString(ipt_curr);
				pApp->ExecuteScriptInString( cmd_line.c_str() );
				pApp->ExecuteScriptFromFile(traj_script.c_str());
			}
			if(!ires) break;

			if(p_mm_mod->update_view_flag && update_time < clock() )
			{
				p_mm_mod->UpdateMolView();
				update_time = clock() + CLOCKS_PER_SEC * p_mm_mod->update_view_interval;
			}

			if( npt_step > 1 )
			{
				if( use_index )
				{
					ipt_curr += npt_step;
				}
				else
				{
					for( int j = 0; j < (npt_step - 1); j++)
					{
						ires = ReadTrajPoint();
						if(!ires) 
						{
							p_mm_mod->to_stop_simulations = TRUE;
							break;
						}
						ipt_curr++;
					}
					ipt_curr++;
				}
			}
			else
			{
				ipt_curr++;
			}
		}
	}
	catch( const std::exception& ex )
	{
		p_mm_mod->to_stop_simulations = TRUE;
		if( files_opened ) CloseAmberTrajFiles();
		PrintLog("Error in MDTrajAnalMod::AnalyzeTrajectoryInternal() \n");
		PrintLog("%s\n",ex.what());
		pmset->SetCrdFromArray(crd_init);
		return FALSE;
	}
	for(int i = 0; i < agents.size(); i++)
	{
		if(agents[i]->IsActive()) agents[i]->Finalize();
	}

	if(!traj_script.empty())
	{
		pApp->ExecuteScriptInString("script_status = SCRIPT_STOP");
		pApp->ExecuteScriptFromFile(traj_script.c_str());
	}

	pmset->SetCrdFromArray(crd_init);
	p_mm_mod->to_stop_simulations = TRUE;
	if( files_opened ) CloseAmberTrajFiles();
	return TRUE;
}

int MDTrajAnalMod::OpenAmberTrajFilesToRead()
{
	p_mm_mod->p_amber_driver->trj_coord_fp = fopen(p_mm_mod->p_amber_driver->amber_trj_coord_file.c_str(),"rb");
	if(p_mm_mod->p_amber_driver->trj_coord_fp == NULL)
	{
		PrintLog("Can't find MD trajectory file %s\n",p_mm_mod->p_amber_driver->amber_trj_coord_file.c_str() );
		return FALSE;
	}	

	int i;
	char buf[256];
	
	fgets(buf,128,p_mm_mod->p_amber_driver->trj_coord_fp);
	
	p_mm_mod->p_amber_driver->trj_ene_fp = fopen(p_mm_mod->p_amber_driver->amber_trj_ene_file.c_str(),"rb");
	if(p_mm_mod->p_amber_driver->trj_ene_fp == NULL)
	{
		PrintLog("Can't find MD energy info file %s",p_mm_mod->p_amber_driver->amber_trj_ene_file.c_str() );
	}
	else
	{
		for(i = 0; i < 10; i++)
		{
			fgets(buf,128,p_mm_mod->p_amber_driver->trj_ene_fp);
		}
	}

	p_mm_mod->p_amber_driver->trj_vel_fp = fopen(p_mm_mod->p_amber_driver->amber_trj_vel_file.c_str(),"rb");
	return TRUE;
}

int MDTrajAnalMod::CloseAmberTrajFiles()
{
	if(p_mm_mod->p_amber_driver->trj_coord_fp != NULL)
	{
		fclose(p_mm_mod->p_amber_driver->trj_coord_fp);
		p_mm_mod->p_amber_driver->trj_coord_fp = NULL;
	}
	if(p_mm_mod->p_amber_driver->trj_ene_fp != NULL)
	{
		fclose(p_mm_mod->p_amber_driver->trj_ene_fp);
		p_mm_mod->p_amber_driver->trj_ene_fp = NULL;
	}
	if(p_mm_mod->p_amber_driver->trj_vel_fp != NULL)
	{
		fclose(p_mm_mod->p_amber_driver->trj_vel_fp);
		p_mm_mod->p_amber_driver->trj_vel_fp = NULL;
	}
	return TRUE;
}

int MDTrajAnalMod::BuildTrajIndex()
{
	char buf[256];
	char* cres;

	if( p_mm_mod->p_amber_driver->trj_coord_fp == NULL ) 
	{
		int ires = OpenAmberTrajFilesToRead();
		if( !ires ) return FALSE;
	}
	rewind( p_mm_mod->p_amber_driver->trj_coord_fp);

	int npt;
	npt = pmset->GetNAtoms();

	int ncrd = 3*npt;
	int nline = ncrd/10;
	if( ncrd % 10 != 0 ) nline++;

	std::list<fpos_t> pos_list;

	cres = fgets(buf,256,p_mm_mod->p_amber_driver->trj_coord_fp);
	if( cres == NULL ) return FALSE;

	int i;
	for(;;)
	{
		fpos_t pos;
		fgetpos( p_mm_mod->p_amber_driver->trj_coord_fp, &pos );

		int ch;
		for( i = 0 ; i < nline; i++ )
		{
			cres = fgets(buf,256,p_mm_mod->p_amber_driver->trj_coord_fp);
			if( cres == NULL ) break;
		}

		if( read_pbox == PBOX_READ || (read_pbox == PBOX_READ_IF_SET && pmset->per_bc->IsSet())  )
		{
			cres = fgets(buf,256,p_mm_mod->p_amber_driver->trj_coord_fp);		
		}
		if( cres == NULL ) break;
		pos_list.push_back(pos);
	}
	
	int np = pos_list.size();
	pt_pos.clear();
	pt_pos.reserve(np);
	std::list<fpos_t>::iterator itr;
	for( itr = pos_list.begin(); itr != pos_list.end(); itr++ )
	{
		pt_pos.push_back(*itr);
	}

	return TRUE;
}

int MDTrajAnalMod::ReadTrajPoint()
{
	char buf[256];

	PrintLog("MDTrajAnalMod::ReadTrajPoint() ipt = %d \n", ipt_curr );

    int iread = FALSE;

	if(p_mm_mod->p_amber_driver->trj_coord_fp != NULL)
	{
		iread = TRUE;
		int npt;

		npt = pmset->GetNAtoms();
		AtomIteratorMolSet aitr(pmset);
		HaVec_double darr(3*npt);

// Reading atom coordinate
		int i, ires;
		double box_x,box_y,box_z;

		for(i = 0 ; i < 3*npt; i++)
		{
			ires = fscanf(p_mm_mod->p_amber_driver->trj_coord_fp,"%lf",&darr[i]);
			if(ires == EOF )
			{
				PrintLog("End of AMBER Trajectory file \n");
				return FALSE;
			}
		}

		if( read_pbox == PBOX_READ || (read_pbox == PBOX_READ_IF_SET && pmset->per_bc->IsSet())  )
		{
		    ires = fscanf(p_mm_mod->p_amber_driver->trj_coord_fp,"%lf %lf %lf",&box_x,&box_y,&box_z);
		    if(ires == EOF )
			{
			   PrintLog("End of AMBER Trajectory file \n");
			   return FALSE;
			}
			pmset->per_bc->SetBox(box_x,box_y,box_z );
		}

		if( wrap_crd && pmset->per_bc->IsSet() )
		{
			MolEditor* p_mol_editor = pmset->GetMolEditor();
			int ires = p_mol_editor->WrapToUnitCell(pmset,pmset->per_bc);
		}

		HaAtom* aptr;
		i = 0;
		for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom() )
		{
			aptr->SetX( darr[3*i]   );
			aptr->SetY( darr[3*i+1] );
			aptr->SetZ( darr[3*i+2] );
			i++;
		}
	}
	if( p_amber_driver->trj_ene_fp != NULL)
	{
		fgets(buf,128,p_amber_driver->trj_ene_fp);
		sscanf(buf + 2, "%d %lf %lf %lf ", &p_mm_mod->p_mm_info->nstep, &p_mm_mod->p_mm_info->time, 
			            &p_mm_mod->p_mm_info->tot_energy, &p_mm_mod->p_mm_info->kin_ene);

		fgets(buf,128,p_amber_driver->trj_ene_fp);
		sscanf(buf + 2, "%lf %lf %lf %lf ", &p_mm_mod->p_mm_info->temp, &p_mm_mod->p_mm_info->temp_solute, 
			  &p_mm_mod->p_mm_info->temp_solv, &p_mm_mod->p_mm_info->press_scale_solute);

		fgets(buf,128,p_amber_driver->trj_ene_fp);
		double box_x,box_y,box_z;
		sscanf(buf + 2, "%lf %lf %lf %lf ", &p_mm_mod->p_mm_info->press_scale_solvent,&box_x,&box_y,&box_z);

//		pmset->per_bc->SetBox(box_x,box_y,box_z);

        fgets(buf,128,p_amber_driver->trj_ene_fp);
		sscanf(buf + 2, "%lf %lf %lf %lf ", &p_mm_mod->p_mm_info->volume,&p_mm_mod->p_mm_info->pres_x, 
			       &p_mm_mod->p_mm_info->pres_y, &p_mm_mod->p_mm_info->pres_z);

		fgets(buf,128,p_amber_driver->trj_ene_fp);
		sscanf(buf + 2,"%lf %lf %lf %lf ",&p_mm_mod->p_mm_info->press,&p_mm_mod->p_mm_info->kin_ene_com_x, 
			   &p_mm_mod->p_mm_info->kin_ene_com_y, &p_mm_mod->p_mm_info->kin_ene_com_z);

		fgets(buf,128,p_amber_driver->trj_ene_fp);
		sscanf(buf + 2,"%lf %lf %lf %lf ",&p_mm_mod->p_mm_info->kin_ene_com, &p_mm_mod->p_mm_info->virial_x, 
			       &p_mm_mod->p_mm_info->virial_y, &p_mm_mod->p_mm_info->virial_z);

		fgets(buf,128,p_amber_driver->trj_ene_fp);
		sscanf(buf + 2,"%lf %lf %lf %lf ",&p_mm_mod->p_mm_info->virial_tot,&p_mm_mod->p_mm_info->pot_ene,
			                 &p_mm_mod->p_mm_info->vdw_ene,&p_mm_mod->p_mm_info->electr_ene);

		fgets(buf,128,p_amber_driver->trj_ene_fp);
		sscanf(buf + 2,"%lf %lf %lf %lf ",&p_mm_mod->p_mm_info->hbond_ene, &p_mm_mod->p_mm_info->bond_ene,
			&p_mm_mod->p_mm_info->vang_ene,&p_mm_mod->p_mm_info->dihed_ene);

		fgets(buf,128,p_amber_driver->trj_ene_fp);
		sscanf(buf + 2,"%lf %lf %lf %lf ",&p_mm_mod->p_mm_info->vdw_ene_14,&p_mm_mod->p_mm_info->electr_ene_14,
			&p_mm_mod->p_mm_info->constraints_ene,&p_mm_mod->p_mm_info->polar_ene);
			
		fgets(buf,128,p_amber_driver->trj_ene_fp);
		sscanf(buf + 2,"%lf %lf %lf %lf ",&p_mm_mod->p_mm_info->av_perm_moment,&p_mm_mod->p_mm_info->av_ind_moment,
			&p_mm_mod->p_mm_info->av_tot_moment,&p_mm_mod->p_mm_info->density);     	
	}
	return iread;
}

int MDTrajAnalMod::LoadCurrPt()
{
	try
	{
		if( p_mm_mod->p_amber_driver->trj_coord_fp == NULL ) throw std::runtime_error(" MD trajectory is not open ");
		if( pt_pos.size() == 0 ) throw std::runtime_error(" Index of MD trajectory is not built ");
		if( ipt_curr < 1 && ipt_curr > pt_pos.size() ) throw std::runtime_error("  Current point index is out of range ");
		fpos_t& pos = pt_pos[ipt_curr - 1];
		fsetpos( p_mm_mod->p_amber_driver->trj_coord_fp, &pos );
		if( feof(p_mm_mod->p_amber_driver->trj_coord_fp) ) throw std::runtime_error(" Error to position file at current point index ");
		ReadTrajPoint();
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in MDTrajAnalMod::LoadCurrPt()  ipt_curr = \n", ipt_curr);
		PrintLog("%s\n",ex.what());
		return FALSE;

	}
	return TRUE;
}

int MDTrajAnalMod::ConvArbalestTrajToAmber( const std::string& md_traj_arbalest,  const std::string& md_traj_amber)
{
	std::ifstream ifs;

	int iunit = 120;
	bool amber_file_opened = false;
	char form_1[20]="(10F8.3)";
	char form_2[20]="(3F8.3)";

	try
	{
		ifs.open(md_traj_arbalest.c_str());
		if( ifs.fail() ) throw std::runtime_error((std::string)"Error to Open file " + md_traj_arbalest );

		std::string str;

		int ires;

		char traj_name_buf[120];
		strncpy(traj_name_buf,md_traj_amber.c_str(),119);

		char fstat[8];
		strcpy(fstat,"UNKNOWN");
		fstat[7] = 0;
	
		char file_format[10];
		strcpy(file_format,"FORMATTED");
		file_format[9] = 0;

		ires = openforf_(&iunit,traj_name_buf, strlen(traj_name_buf),
		                    file_format, strlen(file_format),
							fstat, strlen(fstat));

		if(!ires) throw std::runtime_error(" Unable to open AMBER trajectory for writing "); 

		amber_file_opened = true;

		bool fst_pt = true;
		HaVec_double crd;
		
		for(;;)
		{
			std::getline(ifs,str);
			if( ifs.fail() ) break;
			int ipt; 
			double time;
			int nat;
			int ncrd;
			HaVec_double box(3);

//			PrintLog(" Title string:\n,%s\n",str.c_str());
			std::istringstream ss(str);
			ss >> ipt;
			ss >> time;
			ss >> nat;
			ss >> box[0];
			ss >> box[1];
			ss >> box[2];

			PrintLog(" Ncrd = %d Point = %d time= %8.3f  box = %8.3f, %8.3f, %8.3f \n", ncrd, ipt, time, box[0], box[1], box[2]);    
			if( ss.fail() ) throw std::runtime_error((std::string) " Unable to process title line for MD point:  " + str );

			if( nat != pmset->GetNAtoms() ) throw std::runtime_error( (std::string)"Invalid number of atoms "  + harlem::ToString(nat) );
			
			ncrd = nat*3;
			crd.resize(ncrd);

			if (fst_pt )
			{
				char title_str[120]= "THIS IS ARBALEST MD TRAJECTORY FILE CONVERTED TO AMBER                                                   ";
				wrtfstr_(&iunit,title_str,80);
				fst_pt = false;
			}

			int i;
			for( i = 0; i < nat; i++)
			{
				std::getline(ifs,str);
				if( ifs.fail() ) break;

				std::istringstream ss(str);
				ss >> crd[3*i];
				ss >> crd[3*i+1];
				ss >> crd[3*i+2];
			}

			ires = wrtfarrd_(&iunit,crd.v(),&ncrd,form_1,strlen(form_1));
			if(!ires) std::runtime_error(" Unable write MD point into AMBER trajectory file ");

			int nbx = 3;
			ires = wrtfarrd_(&iunit,box.v(),&nbx,form_2,strlen(form_2));
			if(!ires) std::runtime_error(" Unable write Box info into AMBER trajectory file ");
		}
	}
	catch( const std::exception& ex )
	{
		PrintLog("Error in  MDTrajAnalMod::ConvArbalestTrajToAmber() \n");
		PrintLog(" %s\n",ex.what());
		if( amber_file_opened ) clsforf_(&iunit);
		return FALSE;
	}
	if( amber_file_opened ) clsforf_(&iunit);
	return TRUE;
}

int MDTrajAnalMod::ReduceAmberMDTraj(const char* traj_name_new, PointContainer* p_sub_group_arg )
//! Will use currently set npt_begin and npt_step factors
//!         
{
	PrintLog(" MDTrajAnalMod::ReduceAmberMDTraj() pt 1 \n");
	PointContainer* p_sub_group = this->pmset;
	if( p_sub_group_arg != NULL ) p_sub_group = p_sub_group_arg;

	int iunit = 120;
	int ires;

	char traj_name_buf[120];
	strncpy(traj_name_buf,traj_name_new,119);

	char fstat[8];
	strcpy(fstat,"UNKNOWN");
	fstat[7] = 0;
	
	char file_format[10];
	strcpy(file_format,"FORMATTED");
    file_format[9] = 0;

	ires = openforf_(&iunit,traj_name_buf, strlen(traj_name_buf),
		                    file_format, strlen(file_format),
							fstat, strlen(fstat));

	if(!ires)
	{	
		PrintLog("Error in HaMolMechMod::ReduceAmberMDTraj() \n");
		PrintLog("Unable to open new shortened MD trajectory file \n");
		return FALSE;
	}


	if(p_amber_driver->trj_coord_fp != NULL) 
	{
		CloseAmberTrajFiles(); 
	}
		
	ires = OpenAmberTrajFilesToRead();
	if(!ires)
	{
		clsforf_(&iunit);
		PrintLog("Error in HaMolMechMod::ReduceAmberMDTraj() \n");
		PrintLog("Unable to open MD trajectory file \n");
		return FALSE;
	}
	
	char title_str[120]= "THIS IS A SHORTENED MD TRAJECTORY FILE                                          ";

	wrtfstr_(&iunit,title_str,80);

	int npt = p_sub_group->GetNumPt();
	int ncrd = 3*npt;
	
	HaVec_double crd_grp(3*npt);
	HaVec_double box(3);

	char form_1[20]="(10F8.3)";
	char form_2[20]="(3F8.3)";

	double avx, avy, avz;

	HaMat_double rot_mat(3,3,0.0);
    rot_mat(1,1) = 1.0;  rot_mat(2,2) = 1.0; rot_mat(3,3) = 1.0; 
	HaVec_double trans_vec (3,0.0);
    double eps;

	Vec3DValArray ref_crd;
	
	if( align_to_first_pt || align_to_current_crd ) ref_crd.resize(npt);

	if( align_to_current_crd )
	{
		p_sub_group->SaveCrdToArray(crd_grp);
		ref_crd.SetCrdFromArray(crd_grp);
	}

	ipt_curr = 0;
	for(;;) 
	{
		ires = this->ReadTrajPoint();
		ipt_curr++;
		if(!ires) break;
		
		if( ipt_curr == 1 && align_to_first_pt && !align_to_current_crd )
		{
			p_sub_group->SaveCrdToArray(crd_grp);
			ref_crd.SetCrdFromArray(crd_grp);
		}

		if( ipt_curr < npt_begin ) continue;
		if( ipt_curr > npt_end ) break;
		
		if( npt_step > 1 && ( (ipt_curr - 1) % npt_step != 0)) continue;

		if( align_to_first_pt || align_to_current_crd )
		{
			ires = PointContainer::GetSuperimposeMat( ref_crd, *p_sub_group, rot_mat,  trans_vec, eps);
			p_sub_group->Transform(rot_mat,  trans_vec);
		}
		p_sub_group->SaveCrdToArray(crd_grp);

		ires = wrtfarrd_(&iunit,crd_grp.begin(),&ncrd,form_1,strlen(form_1));
		if(!ires) break;

		box[0] = 100.0; box[1] = 100.0; box[2] = 100.0;
		if( pmset->per_bc->IsSet() )
		{
			box[0] = pmset->per_bc->GetA();
			box[1] = pmset->per_bc->GetB();
			box[2] = pmset->per_bc->GetC();
		}

		if( (write_pbox == PBOX_WRITE_IF_SET && pmset->per_bc->IsSet()) ||  write_pbox == PBOX_WRITE )
		{
			int nbx = 3;
			ires = wrtfarrd_(&iunit,box.begin(),&nbx,form_2,strlen(form_2));
			if(!ires) break;
		}
	}
	
	clsforf_(&iunit);

	CloseAmberTrajFiles();
	return ires;
}

void MDTrajAnalMod::SetAmberMDCrdTraj(const std::string& md_crd_fname )
{
	p_mm_mod->p_amber_driver->amber_trj_coord_file = md_crd_fname;
}

void MDTrajAnalMod::SetAmberMDVelTraj(const std::string& md_vel_fname )
{
	p_mm_mod->p_amber_driver->amber_trj_vel_file = md_vel_fname;
}

void MDTrajAnalMod::SetAmberMDEneTraj(const std::string& md_ene_fname )
{
	p_mm_mod->p_amber_driver->amber_trj_ene_file = md_ene_fname;
}

void MDTrajAnalMod::SetAlignToFstPt( bool set_par )
{
	align_to_first_pt = set_par;
}

void MDTrajAnalMod::SetAlignToCurrentCrd( bool set_par )
{
	align_to_current_crd = set_par;
}

void MDTrajAnalMod::SetReadPBox( int set_par )
{
	read_pbox = set_par;
}

void MDTrajAnalMod::SetWritePBox( int set_par )
{
	write_pbox = set_par;
}

void MDTrajAnalMod::SetWrapCrd( int set_par )
{
	wrap_crd = set_par;
}

void MDTrajAnalMod::SetPtBegin(int npt_begin_par)
{
	npt_begin = npt_begin_par;
}

int MDTrajAnalMod::GetPtBegin() const
{
	return npt_begin;
}

void MDTrajAnalMod::SetPtStep(int npt_step_par)
{
	npt_step = npt_step_par;
}

int MDTrajAnalMod::GetPtStep() const
{
	return npt_step;
}

void MDTrajAnalMod::SetPtEnd(int npt_end_par)
{
	npt_end = npt_end_par;
}

int MDTrajAnalMod::GetPtEnd() const
{
	return npt_end;
}

int MDTrajAnalMod::GetCurrPtIdx() const 
{ 
	return ipt_curr; 
}

int MDTrajAnalMod::AddAgent( TrajAnalAgent* p_agent )
{
	if( p_agent == NULL ) return FALSE;
	std::vector<TrajAnalAgent*>::iterator aitr;
	for( aitr = agents.begin(); aitr != agents.end(); aitr++ )
	{
		if( (*aitr) == p_agent ) return FALSE;
	}
	agents.push_back( p_agent );
	return TRUE;
}

int MDTrajAnalMod::DeleteAgent( TrajAnalAgent* p_agent )
{
	if( p_agent == NULL ) return FALSE;
	std::vector<TrajAnalAgent*>::iterator aitr;
	for( aitr = agents.begin(); aitr != agents.end(); )
	{
		if( (*aitr) == p_agent )
		{
			aitr = agents.erase(aitr);
		}
		else
		{
			aitr++;
		}
	}
	agents.push_back( p_agent );
	return TRUE;
}

int MDTrajAnalMod::AddPythonAgent( PyObject* p_obj)
{
	if( p_obj == NULL ) return FALSE;
	std::vector<PyObject*>::iterator aitr;
	for( aitr = py_agents.begin(); aitr != py_agents.end(); aitr++ )
	{
		if( (*aitr) == p_obj ) return FALSE;
	}
	py_agents.push_back( p_obj );
	return TRUE;
}

void MDTrajAnalMod::PrintAgents()
{
	PrintLog(" Python Agents: \n");
	std::vector<PyObject*>::iterator aitr;
	for( aitr = py_agents.begin(); aitr != py_agents.end(); aitr++ )
	{
		PyObject* p_obj = (*aitr);

	}
}


TrajAnalAgent* MDTrajAnalMod::GetTrajAnalAgent(const char* agent_class_name, int create_flag)
{
	TrajAnalAgent* p_agent = NULL;
	int i;
	int n = agents.size();
	for(i = 0; i < n; i++)
	{
		if (agents[i]->GetClassName() ==  agent_class_name)
		{
			p_agent = agents[i];
			break;
		}
	}
	return p_agent;
}


RMSDAgent* MDTrajAnalMod::GetRMSDAgent( int create_flag )
{
	RMSDAgent* p_agent = (RMSDAgent*) GetTrajAnalAgent("RMSDAgent", FALSE );
	if( p_agent != NULL) return p_agent;
	p_agent = NULL;
	if( create_flag ) 
	{
		p_agent = new RMSDAgent();
		agents.push_back(p_agent);
		p_agent->SetMolSet(pmset);
	}
	return p_agent;
}

AtomCorrAgent*  MDTrajAnalMod::GetAtomCorrAgent( int create_flag )
{
	AtomCorrAgent* p_agent = (AtomCorrAgent*) GetTrajAnalAgent("AtomCorrAgent", FALSE );
	if( p_agent != NULL) return p_agent;
	p_agent = NULL;
	if( create_flag ) 
	{
		p_agent = new AtomCorrAgent(pmset);
		agents.push_back(p_agent);
	}
	return p_agent;
}

RMSDAgent::RMSDAgent()
{
	fname_rmsd_out = "rmsd_md.out";
	fname_rmsd_atom_out = "rmsd_atom.out";
	fname_rmsf_atom_out = "rmsf_atom.out";

	f_rmsd_out = NULL;
	npt = 0;
	pmset = NULL;
	pmset_ref = NULL;

	active_flag = TRUE;
	
	calc_rmsd_per_atom_flag = TRUE;
	calc_rmsf_per_atom_flag = TRUE;
	calc_avg_crd_flag  = FALSE;

	avg_crd_file_name  = "avg_coords.xyz";
}
	
RMSDAgent::~RMSDAgent()
{
	if(f_rmsd_out != NULL) fclose(f_rmsd_out);
}

int RMSDAgent::IsActive() const
{
	return active_flag;
}

void RMSDAgent::SetActive(int active_flag_new)
{
	active_flag = active_flag_new;
}


int RMSDAgent::Init(TrajPointInfo* ppt_info )
{
	try
	{
		if( f_rmsd_out != NULL)
		{
			fclose(f_rmsd_out);
			f_rmsd_out = NULL;
		}

		is_initiated_flag = FALSE;
		if( pmset == NULL ) throw std::runtime_error("Molecular Set is not set ");

		if( fit_atoms.empty() ) SetAtomsFit(*pmset);
		if( fit_atoms.size() != ref_coords_fit.size()) throw std::runtime_error("Size of fit atoms array doesn't match reference coordinate array size");

		if( rmsd_atoms.empty() ) 
		{
			rmsd_atoms = fit_atoms;
			ref_coords_rmsd = ref_coords_fit;
		}

		if( ref_coords_rmsd.size() != rmsd_atoms.size() ) throw std::runtime_error("Size of RMSD atoms array doesn't match reference coordinate array size");

		f_rmsd_out = fopen(fname_rmsd_out.c_str(),"w");
		if( f_rmsd_out == NULL) throw std::runtime_error("Error to open atom rmsd output file: " + fname_rmsd_out );

		rmsd_per_atom.resize( rmsd_atoms.size() );
		rmsd_per_atom = 0.0;

		rmsf_per_atom.resize( rmsd_atoms.size() );
		rmsf_per_atom = 0.0;

		sq_coords.resize( 3*rmsd_atoms.size() );
		sq_coords = 0.0;

		avg_coords.resize( 3*rmsd_atoms.size() );
		avg_coords = 0.0;

	}
	catch( const std::exception& ex )
	{
		PrintLog("Error in RMSDAgent::Init() \n");
		PrintLog("%s\n",ex.what());
		is_initiated_flag = FALSE;
		return FALSE;
	}

	npt = 0;
	is_initiated_flag = TRUE;
	return TRUE;
}

int RMSDAgent::AnalyzePt(TrajPointInfo* ppt_info) 
{
	if(!is_initiated_flag)
	{
		PrintLog("Error in RMSDAgent::AnalyzePt() \n");
		PrintLog("Agent is not initiated \n");
		return FALSE;
	}
	npt++;
	int idx_pt = npt;

	if( ppt_info != NULL && ppt_info->ipt > 0)
	{
		idx_pt = ppt_info->ipt; 
	}

	if( idx_pt == 1 )
	{
		if( ref_crd_fit_type == REFC_FIRST_PT )
		{
			HaVec_double crd_arr;
			fit_atoms.SaveCrdToArray( crd_arr );
			ref_coords_fit.SetCrdFromArray( crd_arr );
			rmsd_atoms.SaveCrdToArray( crd_arr );
			ref_coords_rmsd.SetCrdFromArray( crd_arr );
		}
		if( ref_crd_rmsd_type == REFC_FIRST_PT )
		{

		}
	}
	
	HaMat_double rot_mat(3,3);
	HaVec_double trans_vec(3);
	double rmsd;
	int i;
	int ires = PointContainer::GetSuperimposeMat( ref_coords_fit, fit_atoms, rot_mat, trans_vec, rmsd);
	int na = rmsd_atoms.size();
	if( ires )
	{
		Vec3DValArray rmsd_atoms_trcrd; // transformed atoms coordinates
	
		rmsd_atoms_trcrd.resize(na); 
		for(i = 0; i < na; i++)
		{
			rmsd_atoms_trcrd[i].SetCoordFrom(*rmsd_atoms[i]);
		}
		rmsd_atoms_trcrd.Transform(rot_mat,trans_vec);

		rmsd = 0.0;
		for(i = 0; i < na; i++)
		{
			Vec3D diff = rmsd_atoms_trcrd[i] - ref_coords_rmsd[i];
			double sqd_at = diff.length2();
			rmsd += sqd_at;
			if( calc_rmsd_per_atom_flag ) rmsd_per_atom[i] += sqd_at;

		}
		
		if( na > 0) rmsd = sqrt( rmsd/na );

		for(i = 0; i < na; i++)
		{
			avg_coords[3*i  ] += rmsd_atoms_trcrd[i].GetX();
			avg_coords[3*i+1] += rmsd_atoms_trcrd[i].GetY();
			avg_coords[3*i+2] += rmsd_atoms_trcrd[i].GetZ();

			sq_coords[3*i  ] += rmsd_atoms_trcrd[i].GetX()*rmsd_atoms_trcrd[i].GetX();
			sq_coords[3*i+1] += rmsd_atoms_trcrd[i].GetY()*rmsd_atoms_trcrd[i].GetY();
			sq_coords[3*i+2] += rmsd_atoms_trcrd[i].GetZ()*rmsd_atoms_trcrd[i].GetZ();
		}
	}
	else
	{
		rmsd = -1.0;
	}

	if(f_rmsd_out != NULL)
	{
		fprintf(f_rmsd_out,"%8d  %12.6f\n",idx_pt,rmsd);
	}
	
	return TRUE;
}

int RMSDAgent::Finalize()
{
	if(f_rmsd_out) 
	{
		fclose(f_rmsd_out);
	}
	int na = rmsd_atoms.size();

	if( calc_rmsd_per_atom_flag && npt > 0 && na > 0)
	{
		int i;
		FILE* fout_rmsd_atom = fopen(fname_rmsd_atom_out.c_str(),"w");
		if( fout_rmsd_atom != NULL )
		{
			for( i = 0; i < na; i++)
			{
				rmsd_per_atom[i] = sqrt(rmsd_per_atom[i]/npt);
				rmsd_atoms[i]->tempf = rmsd_per_atom[i];
				fprintf(fout_rmsd_atom," %s  %12.4f \n", (rmsd_atoms[i]->GetRef(HaAtom::ATOMREF_NO_MOL)).c_str(),rmsd_per_atom[i]); 
			}
			fclose(fout_rmsd_atom);
		}
	}

	if( npt > 0 && na > 0 ) avg_coords.scale(1.0/npt);
	if( npt > 0 && na > 0 ) sq_coords.scale(1.0/npt);

	if( calc_rmsf_per_atom_flag && npt > 0 && na > 0)
	{
		int i;
		FILE* fout_rmsf_atom = fopen(fname_rmsf_atom_out.c_str(),"w");
		if( fout_rmsf_atom != NULL )
		{
			for( i = 0; i < na; i++)
			{
				rmsf_per_atom[i]  = sq_coords[3*i] - avg_coords[3*i]*avg_coords[3*i];
				rmsf_per_atom[i] += sq_coords[3*i+1] - avg_coords[3*i+1]*avg_coords[3*i+1];
				rmsf_per_atom[i] += sq_coords[3*i+2] - avg_coords[3*i+2]*avg_coords[3*i+2];
				rmsf_per_atom[i] = sqrt( rmsf_per_atom[i] );
				fprintf(fout_rmsf_atom," %s  %12.4f \n", (rmsd_atoms[i]->GetRef(HaAtom::ATOMREF_NO_MOL)).c_str(),rmsf_per_atom[i]); 
			}
			fclose(fout_rmsf_atom);
		}
	}

	if( calc_avg_crd_flag && npt > 0 && na > 0)
	{
		HaVec_double crd_save;
		rmsd_atoms.SaveCrdToArray(crd_save);

		if( rmsd_atoms.SetCrdFromArray(avg_coords) )
		{
			AtomSaveOptions save_opt;
			save_opt.save_atom_ref = TRUE;
			save_opt.at_ref_type = HaAtom::ATOMREF_NO_MOL;
			rmsd_atoms.SaveXYZFile(avg_crd_file_name.c_str(), &save_opt);
			rmsd_atoms.SetCrdFromArray(crd_save);
		}
	}
	is_initiated_flag = FALSE;

	return TRUE;
}

int RMSDAgent::SetAtomsFit(const std::string& atom_group_name )
{
	int ires;
	try
	{
		fit_atoms.clear();
		if( pmset == NULL ) throw std::runtime_error(" Molecular set is not set");
		if( atom_group_name.empty() ) throw std::runtime_error(" atom group name is empty");

		if( boost::iequals(atom_group_name,"ALL") )
		{
			ires = SetAtomsFit(*pmset);
		}
		else
		{
			AtomGroup* pat_arr = pmset->GetAtomGroupByID(atom_group_name.c_str());
			if(pat_arr == NULL) throw std::runtime_error(" No atoms group " + atom_group_name );
			ires = SetAtomsFit(*pat_arr);
		}
	}
	catch( const std::exception& ex ) 
	{
		PrintLog(" RMSDAgent::SetFitAtoms() \n");
		PrintLog("%s\n",ex.what());
		return FALSE;
	}
	return ires;
}


int RMSDAgent::SetAtomsFit(AtomContainer& active_atoms_new)
{
	is_initiated_flag = FALSE;
	fit_atoms.clear();
	int na = active_atoms_new.GetNAtoms();

	if( na == 0) 
	{
		PrintLog("Error in RMSDAgent::SetFitAtoms() \n");
		PrintLog("Empty active atoms array \n");
		return FALSE;
	}
	
	fit_atoms.reserve(na);
	HaAtom* aptr;

	AtomIteratorGen aitr(&active_atoms_new);
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		fit_atoms.push_back(aptr);
	}

	SetRefCrdFit(fit_atoms);

	return TRUE;
}

int RMSDAgent::SetAtomsRMSD(const std::string& atom_group_name)
{
	int ires;
	try
	{
		rmsd_atoms.clear();
		if( pmset == NULL ) throw std::runtime_error(" Molecular set is not set");
		if( atom_group_name.empty() ) throw std::runtime_error(" atom group name is empty");

		if( boost::iequals(atom_group_name,"ALL") )
		{
			ires = SetAtomsFit(*pmset);
		}
		else
		{
			AtomGroup* pat_arr = pmset->GetAtomGroupByID(atom_group_name.c_str());
			if(pat_arr == NULL) throw std::runtime_error(" No atoms group " + atom_group_name );
			ires = SetAtomsRMSD(*pat_arr);
		}
	}
	catch( const std::exception& ex ) 
	{
		PrintLog(" RMSDAgent::SetAtomsRMSD() \n");
		PrintLog("%s\n",ex.what());
		return FALSE;
	}
	return ires;
}

int RMSDAgent::SetAtomsRMSD(AtomContainer& active_atoms_new)
{
	rmsd_atoms.clear();
	int na = active_atoms_new.GetNAtoms();

	if( na == 0) 
	{
		PrintLog("Error in RMSDAgent::SetAtomsRMSD() \n");
		PrintLog("Empty atoms array \n");
		return FALSE;
	}
	
	rmsd_atoms.reserve(na);
	HaAtom* aptr;

	AtomIteratorGen aitr(&active_atoms_new);
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		rmsd_atoms.push_back(aptr);
	}

	SetRefCrdRMSD(rmsd_atoms);

	return TRUE;
}

int RMSDAgent::SetRefCrdFit(PointContainer& ref_coords_new, int ref_crd_type_new )
{
	ref_coords_fit.clear();

	int na = ref_coords_new.GetNumPt();

	if( na == 0) 
	{
		PrintLog("Error in RMSDAgent::SetRefCrdFit() \n");
		PrintLog("Empty reference atom coordinates \n");
		return FALSE;
	}

	if( fit_atoms.size() != na)
	{
		PrintLog("Error in RMSDAgent::SetRefCrdFit() \n");
		PrintLog("Size active atoms array  = %d  not equal reference atom array = %d\n",
			     fit_atoms.size(), na );
		return FALSE;
	}
	
	ref_coords_fit.resize(na);
	Vec3D* pptr;
	int i = 0;

	PointIteratorGen pitr(ref_coords_new);
	for(pptr = pitr.GetFirstPt(); pptr; pptr = pitr.GetNextPt())
	{
		ref_coords_fit[i].SetCoordFrom(*pptr);
		i++;
	}
	if( ref_crd_type_new < 0 || ref_crd_type_new > REFC_SPECIAL ) 
	{
		ref_crd_fit_type = REFC_SPECIAL;
	}
	else
	{
		ref_crd_fit_type = (RefCrdType) ref_crd_type_new;
	}
	return TRUE;
}

int RMSDAgent::SetRefCrdRMSD(PointContainer& ref_coords_new, int ref_crd_type_new )
{
	ref_coords_rmsd.clear();

	int na = ref_coords_new.GetNumPt();

	if( na == 0) 
	{
		PrintLog("Error in RMSDAgent::SetRefCrdRMSD() \n");
		PrintLog("Empty reference atom coordinates \n");
		return FALSE;
	}

	if( rmsd_atoms.size() != na)
	{
		PrintLog("Error in RMSDAgent::SetRefCrdRMSD() \n");
		PrintLog("Size active atoms array  = %d  not equal reference atom array = %d\n",
			     rmsd_atoms.size(), na );
		return FALSE;
	}
	
	ref_coords_rmsd.resize(na);
	Vec3D* pptr;
	int i = 0;

	PointIteratorGen pitr(ref_coords_new);
	for(pptr = pitr.GetFirstPt(); pptr; pptr = pitr.GetNextPt())
	{
		ref_coords_rmsd[i].SetCoordFrom(*pptr);
		i++;
	}
	if( ref_crd_type_new < 0 || ref_crd_type_new > REFC_SPECIAL ) 
	{
		ref_crd_rmsd_type = REFC_SPECIAL;
	}
	else
	{
		ref_crd_rmsd_type = (RefCrdType) ref_crd_type_new;
	}
	return TRUE;
}


//int RMSDAgent::SetAtoms(AtomContainer& active_atoms_new, PointContainer& ref_coords_new)
//{
//	fit_atoms.clear();
//	ref_fit_coords.clear();
//
//	int na = active_atoms_new.GetNAtoms();
//
//	if( na == 0) 
//	{
//		PrintLog("Error in RMSDAgent::SetAtoms() \n");
//		PrintLog("Empty active atoms array \n");
//		return FALSE;
//	}
//	if( ref_coords_new.GetNumPt() != na)
//	{
//		PrintLog("Error in RMSDAgent::SetAtoms() \n");
//		PrintLog("Size active atoms array  = %d  not equal reference atom array = % d\n",
//			     na, ref_coords_new.GetNumPt() );
//		return FALSE;
//	}
//
//	fit_atoms.reserve(na);
//	HaAtom* aptr;
//
//	AtomIteratorGen aitr(&active_atoms_new);
//	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
//	{
//		fit_atoms.push_back(aptr);
//	}
//
//	int ires = SetRefCrd(ref_coords_new);
//
//	return ires;
//}

int RMSDAgent::SetMolSet(MolSet* pmset_new)
{
	pmset = pmset_new;
	return TRUE;
}

int RMSDAgent::SetRefCrdFitFromXYZFile( const std::string& ref_crd_file_name_new )
{
	std::ifstream refc_fs( ref_crd_file_name_new.c_str() );
	if( !refc_fs.good() )
	{	
		PrintLog(" Error in RMSDAgent::SetRefCrdFromXYZFile() \n");
		PrintLog(" Error to open file %s \n", ref_crd_file_name_new.c_str() );
		return FALSE;
	}

	char buf[256];
	int na;
	refc_fs.getline(buf,256);
	istrstream ss(buf);
	ss >> na;
	if(!ss)
	{
		PrintLog(" Error in RMSDAgent::SetRefCrdFromXYZFile() \n");
		PrintLog(" Error reating atom number from file %s \n", ref_crd_file_name_new.c_str() );
		return FALSE;
	}

	if( na != fit_atoms.size())
	{
		PrintLog(" Error in RMSDAgent::SetRefCrdFromXYZFile() \n");
		PrintLog(" Number of atoms = %d in file %s \n", na, ref_crd_file_name_new.c_str() );
		PrintLog(" Is different from the number of active atoms %d \n",fit_atoms.GetNAtoms()); 
		return FALSE;
	}
	ref_coords_fit.resize(na);
	int i;
	for(i = 0; i < na; i++)
	{
		refc_fs.getline(buf,256);
		istrstream line_s(buf);
		std::string id;
		int idx;
		double x,y,z;
		line_s >> idx;
		line_s >> id;
		line_s >> x;
		line_s >> y;
		line_s >> z;
		if( !line_s ) 
		{
			PrintLog(" Error in RMSDAgent::SetRefCrdFromXYZFile() \n");
			PrintLog(" Processing line %s \n", buf);
			ref_coords_fit.clear();
			return FALSE;
		}
		ref_coords_fit[i].SetX_Ang(x);
		ref_coords_fit[i].SetY_Ang(y);
		ref_coords_fit[i].SetZ_Ang(z);
	}

	return TRUE;
}

int RMSDAgent::SetRefCrdRMSDFromXYZFile( const std::string& ref_crd_file_name_new )
{
	std::ifstream refc_fs( ref_crd_file_name_new.c_str() );
	if( !refc_fs.good() )
	{	
		PrintLog(" Error in RMSDAgent::SetRefCrdRMSDFromXYZFile() \n");
		PrintLog(" Error to open file %s \n", ref_crd_file_name_new.c_str() );
		return FALSE;
	}

	char buf[256];
	int na;
	refc_fs.getline(buf,256);
	istrstream ss(buf);
	ss >> na;
	if(!ss)
	{
		PrintLog(" Error in RMSDAgent::SetRefCrdRMSDFromXYZFile() \n");
		PrintLog(" Error reating atom number from file %s \n", ref_crd_file_name_new.c_str() );
		return FALSE;
	}

	if( na != rmsd_atoms.size())
	{
		PrintLog(" Error in RMSDAgent::SetRefCrdRMSDFromXYZFile() \n");
		PrintLog(" Number of atoms = %d in file %s \n", na, ref_crd_file_name_new.c_str() );
		PrintLog(" Is different from the number of atoms to compute RMSD %d \n",rmsd_atoms.GetNAtoms()); 
		return FALSE;
	}

	ref_coords_rmsd.resize(na);
	int i;
	for(i = 0; i < na; i++)
	{
		refc_fs.getline(buf,256);
		istrstream line_s(buf);
		std::string id;
		int idx;
		double x,y,z;
		line_s >> idx;
		line_s >> id;
		line_s >> x;
		line_s >> y;
		line_s >> z;
		if( !line_s ) 
		{
			PrintLog(" Error in RMSDAgent::SetRefCrdRMSDFromXYZFile() \n");
			PrintLog(" Processing line %s \n", buf);
			ref_coords_rmsd.clear();
			return FALSE;
		}
		ref_coords_rmsd[i].SetX_Ang(x);
		ref_coords_rmsd[i].SetY_Ang(y);
		ref_coords_rmsd[i].SetZ_Ang(z);
	}

	return TRUE;
}

int RMSDAgent::SetRefCrdFitFromAtomGroup( const std::string& at_grp_id, MolSet* pmset_ref_new)
{
	if(pmset_ref_new == NULL) 
	{
		PrintLog(" Error in RMSDAgent::SetRefCrdFromAtomGroupID() \n");
		PrintLog(" pmset_ref_new == NULL \n");
		return FALSE;
	}
	
	AtomGroup* p_refc_arr = pmset_ref_new->GetAtomGroupByID(at_grp_id.c_str());
	if( p_refc_arr == NULL) 
	{
		PrintLog(" Error in RMSDAgent::SetRefCrdFromAtomGroupID() \n");
		PrintLog(" Atom Group with id = %s not found \n", at_grp_id.c_str() );
		return FALSE;
	}

	if( p_refc_arr->GetNAtoms() != fit_atoms.GetNAtoms() )
	{
		PrintLog(" Error in RMSDAgent::SetRefCrdFromAtomGroupID() \n");
		PrintLog(" number of Atoms in Ref Crd Array = %d not equal to atom numberes in active array %d = \n",
			       p_refc_arr->GetNAtoms(),fit_atoms.GetNAtoms());
		return FALSE;
	}

	SetRefCrdFit(*p_refc_arr,REFC_ATOM_ARRAY_ID);
	return TRUE;
}

int RMSDAgent::SetRefCrdRMSDFromAtomGroup( const std::string& at_grp_id, MolSet* pmset_ref_new)
{
	if(pmset_ref_new == NULL) 
	{
		PrintLog(" Error in RMSDAgent::SetRefCrdRMSDFromAtomGroup() \n");
		PrintLog(" pmset_ref_new == NULL \n");
		return FALSE;
	}
	
	AtomGroup* p_refc_arr = pmset_ref_new->GetAtomGroupByID(at_grp_id.c_str());
	if( p_refc_arr == NULL) 
	{
		PrintLog(" Error in RMSDAgent::SetRefCrdRMSDFromAtomGroup() \n");
		PrintLog(" Atom Group with id = %s not found \n",at_grp_id.c_str() );
		return FALSE;
	}

	if( p_refc_arr->GetNAtoms() != rmsd_atoms.GetNAtoms() )
	{
		PrintLog(" Error in RMSDAgent::SetRefCrdRMSDFromAtomGroup() \n");
		PrintLog(" number of Atoms in Ref Crd Array = %d not equal to atom numberes in active array %d = \n",
			       p_refc_arr->GetNAtoms(),rmsd_atoms.GetNAtoms());
		return FALSE;
	}

	SetRefCrdRMSD(*p_refc_arr,REFC_ATOM_ARRAY_ID);
	return TRUE;
}

AtomCorrAgent::AtomCorrAgent( MolSet* pmset_par )
{
	pmset = pmset_par;
	npt_proc = 0;
	SetDistRange(1.5, 6.0, 46);
	fname_out = "at_corr.out";
}

AtomCorrAgent::~AtomCorrAgent()
{
		
}


int AtomCorrAgent::IsActive() const
{
//	PrintLog(" AtomCorrAgent::IsActive() pt 1 at_grp_1.size() = %d \n", at_grp_1.size() );
	if( at_grp_1.empty() ) return FALSE;
//	PrintLog(" AtomCorrAgent::IsActive() pt 2 at_grp_2.size() = %d \n", at_grp_2.size()  );
	if( at_grp_2.empty() ) return FALSE;
//	PrintLog(" AtomCorrAgent::IsActive() pt 3 rmin = %9.3f  rmax = %9.3f \n", rmin, rmax );
	if( (rmax - rmin) < 0.1 ) return FALSE;
//	PrintLog(" AtomCorrAgent::IsActive() pt 4 nr = %d \n", nr );
	if( nr < 1 ) return FALSE;
	return TRUE;
}

	
void AtomCorrAgent::SetActive(int active_flag)
{
	
}

int AtomCorrAgent::Init(TrajPointInfo* ppt_info )
{
	PrintLog(" AtomCorrAgent::Init() pt 1 \n");
	npt_proc = 0;
	SetDistRange(rmin,rmax,nr);

	int na1 = at_grp_1.size();
	int na2 = at_grp_2.size();
	int npair = na1*na2;
	if( npair == 0 ) return FALSE;

	return TRUE;
}

int  AtomCorrAgent::AnalyzePt(TrajPointInfo* ppt_info )
{
//	PrintLog(" AtomCorrAgent::AnalyzePt() pt 1 \n");
	int na1 = at_grp_1.size();
	int na2 = at_grp_2.size();
	int npair = na1*na2;
	if( npair == 0 ) return FALSE;

	double dr = 0;
	if( (nr > 1) && (rmax > rmin) )
	{
		dr = (rmax - rmin)/(nr-1);
	}

	int i,j,k;
	for( i = 0; i < na1; i++ )
	{
		HaAtom* aptr1 = at_grp_1.at(i);
		HaVec_int nv(nr);
		nv = 0;
		for( j = 0; j < na2; j++ )
		{
			HaAtom* aptr2 = at_grp_2.at(j);
			double dist = Vec3D::CalcDistance(aptr1,aptr2);
			for( k = 0; k < nr; k++ )
			{
				double r = rmin + dr*k;
				if( dist < r ) nv[k]++;
			}
		}
		for( k = 0; k < nr; k++ )
		{
			nr_acc[k]   += nv[k]; 
			n2r_acc[k]  += nv[k]*nv[k];
		}
		npt_proc++;
	}

	return TRUE;
}
	
int  AtomCorrAgent::Finalize()
{
	PrintLog(" AtomCorrAgent::Finalize() pt 1 \n");
	char buf[256];

	if( npt_proc == 0 ) return FALSE;
	double dr = 0;
	if( (nr > 1) && (rmax > rmin) )
	{
		dr = (rmax - rmin)/(nr-1);
	}

	HaVec_double rcut =  GetRCut();
	HaVec_double r_gr(nr); r_gr = 0.0;


	double dens = 1.0;
	if( pmset == NULL ) return FALSE;
	if( pmset->per_bc->IsSet() )
	{
		double vol = pmset->per_bc->GetA() * pmset->per_bc->GetB() * pmset->per_bc->GetC();
		dens = at_grp_2.size()/vol;
	}

	int k;
	for( k = 0; k < nr; k++ )
	{
		navg[k]   = ((double)nr_acc[k] )/npt_proc;
		deltn2[k] = ((double)n2r_acc[k])/npt_proc - navg[k]*navg[k];  
		if( k > 0 && dr > 0.0 )
		{
			double r = rmin + dr*(k-0.5);
			r_gr[k-1] = r;
			gr[k-1] = (navg[k] - navg[k-1])/(4.0*PI*dens*r*r*dr);
		}
	}

	ofstream os(fname_out.c_str());
	if( os.is_open() )
	{
		os << "  i_pt        r           <n>    <n - <n>>^2  g(r) " << std::endl;
		for( k = 0; k < nr; k++ )
		{
			sprintf(buf," %5d  %9.4f  %9.4f  %9.4f  %9.4f ", k, rcut[k], navg[k], deltn2[k], gr[k]);
			os << buf << std::endl;
		}
	}

	return TRUE;
}


int AtomCorrAgent::SetAtGroup1ByExpr(const std::string& expr)
{
	at_grp_1.SetFromExprStr(expr.c_str(),pmset);
	return at_grp_1.size();
}

int AtomCorrAgent::SetAtGroup2ByExpr(const std::string& expr)
{
	at_grp_2.SetFromExprStr(expr.c_str(),pmset);
	return at_grp_2.size();
}

void AtomCorrAgent::SetDistRange( double rmin_par, double rmax_par, int nr_par )
{
	rmin = rmin_par;
	rmax = rmax_par;
	if( nr_par < 1 ) 
	{
		nr = (rmax - rmin )/0.1 + 1.000000001;
	}
	else
	{
		nr = nr_par;
	}

	if( nr < 1 ) nr = 50;
	ra.resize(nr);      ra = 0.0;
	nr_acc.resize(nr);  nr_acc = 0;
	n2r_acc.resize(nr); n2r_acc = 0;
	gr.resize(nr);      gr = 0.0;
	navg.resize(nr);    navg = 0;
	deltn2.resize(nr);  deltn2 = 0.0;
}

HaVec_double AtomCorrAgent::GetRCut() const
{
	double dr = 0;
	if( (nr > 1) && (rmax > rmin) )
	{
		dr = (rmax - rmin)/(nr-1);
	}

	HaVec_double rcut(nr);
	int k;
	for( k = 0; k < nr; k++ )
	{
		if( k == 0 ) 
		{
			rcut[k] = rmin;
		}
		else
		{
			rcut[k] = rcut[k-1] + dr;
		}
	}
	return rcut;
}

MDTrajectoryIOAgent::MDTrajectoryIOAgent(MMDriverAmber* p_mm_driver_new)
{
	p_mm_driver   = p_mm_driver_new;
	p_amber_model = p_mm_driver->p_amber_model;
	p_mm_mod      = p_mm_driver->p_mm_mod;
	p_mm_model    = p_mm_mod->p_mm_model;
}

MDTrajectoryIOAgent::~MDTrajectoryIOAgent()
{

}

std::string MDTrajectoryIOAgent::GetClassName() const 
{ 
	return "MDTrajectoryIOAgent";
}  

int MDTrajectoryIOAgent::Init()
{	
	sit.resize(MMDriverAmber::SI_CNT);
	sit = 0.0;
	sit_tmp.resize(MMDriverAmber::SI_CNT);
	sit_tmp = 0.0;
    sit2.resize(MMDriverAmber::SI_CNT);
	sit2 = 0.0;
	sit2_tmp.resize(MMDriverAmber::SI_CNT);
	sit2_tmp = 0.0;

	nvalid = 0;

	n_saved_crd = p_amber_model->natom * 3;
    if (p_mm_mod->limit_wrt_atoms > 0) n_saved_crd = p_mm_mod->limit_wrt_atoms * 3;

	if( p_mm_mod->wrt_constr_freq > 0 && p_mm_mod->p_mm_model->DistConstraints.size() > 0 )
	{
		constr_trj_fs.open(p_mm_mod->constr_trj_fname.c_str());
		if( !constr_trj_fs.is_open() )
		{
			PrintLog(" Error to open Constraints info file %s\n", p_mm_mod->constr_trj_fname.c_str());
		}	
	}

	return TRUE;
}	

int MDTrajectoryIOAgent::AnalyzePt(int nstep, HaVec_double& sys_info )
{
	char buf[256];

	int i_0 = 0;
	int i_1 = 1;
	int i_3 = 3;
	int i_7 = 7;

	int m;

	int prime_ham = TRUE; //!< indicate primary hamiltonian (in TI) - always true in simple MD 

	if( p_mm_mod->run_ti && p_mm_mod->inter_model_rank != 0 ) 
	{
		prime_ham = FALSE;
	}

	if( p_mm_driver->master )
	{
	// nvalid is the number of steps where all energies are calculated.

		if(p_mm_mod->run_ti)
		{
			if(p_mm_mod->inter_model_rank == 0)
			{
				fprintf(p_mm_mod->p_ti_mod->file_dvdl,"%8d %16.9f\n",nstep,sys_info[MMDriverAmber::SI_DV_DLAMBDA]);
			}
		}

		if ( (nstep % p_mm_driver->nrespa) == 0 )
		{
			nvalid++;

			for(m = 0; m < MMDriverAmber::SI_CNT; m++)
			{
				sit[m]  += sys_info[m];
				sit2[m] += sys_info[m] * sys_info[m];
			}
		}
	// Restrt:
		int write_restrt = FALSE;

		if (nstep == p_mm_mod->num_md_steps) 
		{
			write_restrt = TRUE;
		}
		else if (p_mm_mod->wrt_rstrt_freq != 0) 
		{
			if ( (nstep % p_mm_mod->wrt_rstrt_freq) == 0 && (nstep % p_mm_driver->nrespa) == 0 ) write_restrt = TRUE;
		}

		if(!prime_ham) write_restrt = FALSE;

		if (write_restrt) 
		{
			if (p_mm_mod->wrap_coord == 0) 
			{
				FC_FUNC_MODULE(runfiles_mod,mdwrit)(&p_mm_mod->run_type.value(),&p_mm_mod->wrt_rstrt_freq, &p_mm_mod->write_coord_format.value(), 
					                    &p_mm_mod->period_bcond.value(), &nstep, &p_amber_model->natom, 
					p_mm_driver->atm_crd.v(), p_mm_driver->pbc_box.v(), p_mm_driver->atm_vel.v(), &p_mm_driver->t);
			}
			else
			{
				HaVec_double crd_copy(p_mm_driver->atm_crd);
				if (p_mm_mod->period_bcond != p_mm_mod->period_bcond.NO_PERIODICITY )
				{
					FC_FUNC_MODULE(pbc_mod,wrap_molecules)(&p_amber_model->nspm, p_amber_model->atm_nsp.v(),crd_copy.v());
					if ( p_mm_model->GetMolSet()->per_bc->IsOctahedron() ) FC_FUNC_MODULE(pbc_mod,wrap_to)(&p_amber_model->nspm, p_amber_model->atm_nsp.v(), crd_copy.v(), p_mm_driver->pbc_box.v());
				}
				FC_FUNC_MODULE(runfiles_mod,mdwrit)(&p_mm_mod->run_type.value(),&p_mm_mod->wrt_rstrt_freq, &p_mm_mod->write_coord_format.value(), 
					                    &p_mm_mod->period_bcond.value(), &nstep, &p_amber_model->natom, 
					crd_copy.v(), p_mm_driver->pbc_box.v(), p_mm_driver->atm_vel.v(), &p_mm_driver->t);
			}
		}

		// Coordinate archive:

		if (p_mm_mod->wrt_coord_freq > 0 && prime_ham) 
		{
			if ( (nstep % p_mm_mod->wrt_coord_freq) == 0) 
			{
				if (p_mm_mod->wrap_coord == 0) 
				{
					 FC_FUNC_MODULE(runfiles_mod,corpac)(&n_saved_crd, p_mm_driver->atm_crd.v(),&i_1, &p_mm_driver->mdcrd);
				}
				else
				{
					HaVec_double crd_copy(p_mm_driver->atm_crd);
					if (p_mm_mod->period_bcond != p_mm_mod->period_bcond.NO_PERIODICITY)
					{
						FC_FUNC_MODULE(pbc_mod,wrap_molecules)(&p_amber_model->nspm, p_amber_model->atm_nsp.v(),crd_copy.v());
						if (p_mm_model->GetMolSet()->per_bc->IsOctahedron() && p_mm_mod->wrap_coord == 1)  
						{
							FC_FUNC_MODULE(pbc_mod,wrap_to)(&p_amber_model->nspm, p_amber_model->atm_nsp.v(), crd_copy.v(), 
								                            p_mm_driver->pbc_box.v() );
						}
					}
					FC_FUNC_MODULE(runfiles_mod,corpac)(&n_saved_crd, crd_copy.v(), &i_1, &p_mm_driver->mdcrd );
				}

				if ( p_mm_mod->period_bcond != p_mm_mod->period_bcond.NO_PERIODICITY ) 
				{
					if ( p_mm_mod->traj_wrt_format == p_mm_mod->traj_wrt_format.BINARY ) 
					{
						PrintLog(" No binary MD trajectory save in HARLEM \n");
						// bintraj_mod_mp_write_binary_cell_dat_(p_mm_driver->pbc_box.v());
					}
					else
					{
						FC_FUNC_MODULE(runfiles_mod,corpac)(&i_3, p_mm_driver->pbc_box.v(),&i_1, &p_mm_driver->mdcrd);
					} 
				}

				// For BINARY format, coordinate archive may also contain vels.

				if (p_mm_mod->traj_wrt_format == p_mm_mod->traj_wrt_format.BINARY && p_mm_driver->p_mm_mod->wrt_vel_freq < 0) 
				{
					FC_FUNC_MODULE(runfiles_mod,corpac)(&n_saved_crd, p_mm_driver->atm_vel.v(),&i_1, &p_mm_driver->mdcrd);
				}
			}
		} 
		// Velocity archive:
		if (p_mm_mod->wrt_vel_freq > 0 && prime_ham ) 
		{
			if ((nstep % p_mm_mod->wrt_vel_freq) == 0) 
			{
				FC_FUNC_MODULE(runfiles_mod,corpac)(&n_saved_crd, p_mm_driver->atm_vel.v(),&i_1, &p_mm_driver->mdvel);
			}
		}

		// Energy archive:

		if (p_mm_mod->wrt_ener_freq > 0 && prime_ham) 
		{
			if ( (nstep % p_mm_mod->wrt_ener_freq) == 0 && (nstep % p_mm_driver->nrespa) == 0 ) 
			{
				FC_FUNC_MODULE(runfiles_mod,mdeng)(&nstep, &p_mm_driver->t, sys_info.v(), p_mm_driver->pbc_box.v());
			}
		}

		if ( (nstep % p_mm_mod->wrt_log_freq) == 0 && (nstep % p_mm_driver->nrespa) == 0 ) 
		{
			p_mm_driver->PrintMDEneMDOUT(sys_info,nstep, p_mm_driver->t);
		}

		// Flush binary netCDF file(s) if necessary...

		if ( p_mm_mod->traj_wrt_format == p_mm_mod->traj_wrt_format.BINARY && prime_ham) 
		{
			if (p_mm_mod->wrt_coord_freq > 0) 
			{
				if ((nstep % p_mm_mod->wrt_coord_freq ) == 0) 
				{
					PrintLog(" No binary MD trajectory save in HARLEM \n");
					// bintraj_mod_mp_end_binary_frame_(&p_mm_driver->mdcrd);
				}
			}

			if (p_mm_mod->wrt_vel_freq > 0) 
			{
				if ((nstep % p_mm_mod->wrt_vel_freq) == 0) 
				{
					PrintLog(" No binary MD trajectory save in HARLEM \n");
					// bintraj_mod_mp_end_binary_frame_(&p_mm_driver->mdvel);
				}
			}
		}

		if( p_mm_mod->wrt_constr_freq > 0 && p_mm_mod->p_mm_model->DistConstraints.size() > 0 )
		{
			if( nstep % p_mm_mod->wrt_constr_freq == 0 && constr_trj_fs.is_open() )
			{
				int nc = p_mm_model->DistConstraints.size();
				int ic;
				
				sprintf(buf,"%10d ",nstep);
				constr_trj_fs << buf;

				for(ic = 0 ; ic < nc; ic++ )
				{
//					AtomContact* p_cnt = &p_mm_model->DistConstraints[ic];
//					sprintf(buf,"%16.9f ",Vec3D::CalcDistance(p_cnt->pt1,p_cnt->pt2) );

					int i1 = p_mm_model->p_amber_model->dist_constr_idx[3*ic] - 1;
					int i2 = p_mm_model->p_amber_model->dist_constr_idx[3*ic+1] - 1;

					double dx = p_mm_driver->atm_crd[3*i1]   - p_mm_driver->atm_crd[3*i2]; 
					double dy = p_mm_driver->atm_crd[3*i1+1] - p_mm_driver->atm_crd[3*i2+1];
					double dz = p_mm_driver->atm_crd[3*i1+2] - p_mm_driver->atm_crd[3*i2+2];

					double r = sqrt( dx*dx + dy*dy + dz*dz );
					sprintf(buf,"%16.9f ", r );

					constr_trj_fs << buf;
				}
				constr_trj_fs << std::endl;
			}
		}

		// Output running average:

		if (p_mm_driver->ntave > 0)
		{
			if ( (nstep % p_mm_driver->ntave) == 0 && (nstep % p_mm_driver->nrespa) == 0 )
			{
				p_mm_model->UpdateDataFromFort();
				p_mm_driver->PrintLogMDOUT("=============================================================================\n");

				double tspan = p_mm_driver->ntave/p_mm_driver->nrespa;

				// Coming into this loop, the _tmp variables hold the values of
				// sit, sit2 when this routine was last called (or 0.d0).  The _tmp
				// vars are used here as scatch space and then updated with the
				// current sit, sit2.

				for(m = 0; m < MMDriverAmber::SI_CNT; m++)
				{ 
					sit_tmp[m] = (sit[m] - sit_tmp[m]) / tspan;
					sit2_tmp[m] = (sit2[m] - sit2_tmp[m]) / tspan - sit_tmp[m] * sit_tmp[m];
					if (sit2_tmp[m] < 0.0) sit2_tmp[m] = 0.0;
					sit2_tmp[m] = sqrt(sit2_tmp[m]);
				}

				p_mm_driver->PrintLogMDOUT("     A V E R A G E S   O V E R   %7d    S T E P S \n", p_mm_driver->ntave/p_mm_driver->nrespa );

				p_mm_driver->PrintMDEneMDOUT(sit_tmp,nstep, p_mm_driver->t);
				
				p_mm_driver->PrintLogMDOUT("     R M S  F L U C T U A T I O N S   \n");
				p_mm_driver->PrintMDEneMDOUT(sit2_tmp,nstep, p_mm_driver->t);
				p_mm_driver->PrintLogMDOUT("=============================================================================\n");

				sit_tmp  = sit;
				sit2_tmp = sit2;
			}
		}	
	}
	return TRUE;
}

int MDTrajectoryIOAgent::Finalize(int nstep)
{
	int i_0 = 0;
	int i_1 = 1;

	int prime_ham = TRUE;

	if( p_mm_mod->run_ti && p_mm_mod->inter_model_rank != 0 ) 
	{
		prime_ham = FALSE;
	}

	if( constr_trj_fs.is_open() ) constr_trj_fs.close();

	if (p_mm_driver->master) 
	{
		double tspan = nvalid;

		if (nvalid > 0) 
		{
			for(int m = 0; m < MMDriverAmber::SI_CNT; m++)
			{
				sit[m] = sit[m]/tspan;
				sit2[m] = sit2[m]/tspan - sit[m] * sit[m];
				if (sit2[m] < 0.0) sit2[m] = 0.0;
				sit2[m] =  sqrt(sit2[m]);
			}

			p_mm_driver->PrintLogMDOUT("     A V E R A G E S   O V E R   %7d    S T E P S \n", nvalid );

			p_mm_driver->PrintMDEneMDOUT(sit,nstep, p_mm_driver->t);
			p_mm_driver->PrintLogMDOUT("     R M S  F L U C T U A T I O N S   \n");
			p_mm_driver->PrintMDEneMDOUT(sit2,nstep, p_mm_driver->t);

      // Print Born radii statistics:

			if (p_amber_model->using_gb_potential)
			{
				if (p_amber_model->rbornstat == 1) 
				{
					p_mm_driver->PrintLogMDOUT("   STATISTICS OF EFFECTIVE BORN RADII OVER  %7d  STEPS ",nstep);
					FC_FUNC_MODULE(gb_ene_mod,print_born_radii_stat)(&tspan);
				}
			}
		} // (nvalid > 0)
	}
	return TRUE;
}

