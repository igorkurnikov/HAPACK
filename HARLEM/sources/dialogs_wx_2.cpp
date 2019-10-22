/*! \file dialogs_wx_2.cpp

    wxWidgets Dialogs 2 in HARLEM implementation
 
    \author Igor Kurnikov  
    \date 2005-
*/

#include <mpi.h>

#include "wx/wxprec.h"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include "wx/valgen.h"
#include "wx/grid.h"
#include "wx/dir.h"
#include "wx/filename.h"
#include "wx/listctrl.h"
#include "wx/tokenzr.h"
#include "wx/notebook.h"
#include <wx/checkbox.h>
#include <wx/clipbrd.h>

#include "ctrl_wx.h"
#include "ha_wx_aux_1.h"
#include <fstream>
#include "hastring.h"
#include "canvas3d.h"
#include "dialogs_wx_1.h"
#include "dialogs_wx_2.h"
#include "ha_wx_res_2_wdr.h"
#include "haatgroup.h"
#include "trajanal.h"
#include "haintermol.h"
#include "hamolset.h"
#include "hamolecule.h"
#include "hamolview.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "haempirical.h"
#include "hamainframe_wx.h"

#include "ha_wx_res_mikola_wdr.h"//<mikola 26March08
#include "object3d.h"//<mikola 26March08
#include "hasurface.h"//<mikola 26March08
/////////////////////////////////////////////////////////////////////
//  InterMolDlgWX InterMolecular Interactions Module Dialog:

int InterMolDlgWX::dlg_open = FALSE;

InterMolDlgWX::InterMolDlgWX(HaInterMolMod* ptr_im_mod_new, wxWindow *parent ):
wxFrame( parent, -1, "Intermolecular Interactions Module")
{
	this->SetExtraStyle(wxWS_EX_VALIDATE_RECURSIVELY);
	ptr_im_mod = ptr_im_mod_new;
	if(ptr_im_mod) 
	{
		MolSet* pmset = ptr_im_mod->GetMolSet();
		ptr_emp_mod = pmset->GetEmpiricalMod(true);
		p_mc_sim =  ptr_im_mod->p_mc_sim;
		p_io_ag = p_mc_sim->GetTrajectoryIOAgent();
	}
    dlg_open = true;

	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);

	inter_mol_dlg( this, TRUE );
	OnInitDialog();
}

InterMolDlgWX::~InterMolDlgWX()
{
	dlg_open = FALSE;
}

void
InterMolDlgWX::OnClose(wxCloseEvent &event )
{
	event.Skip();
	dlg_open = FALSE;
}


BEGIN_EVENT_TABLE(InterMolDlgWX,wxFrame)
//	EVT_BUTTON(IDC_INTERMOL_INIT, InterMolDlgWX::OnIntermolInit)
	EVT_BUTTON(IDC_INTERMOL_EL_STAT, InterMolDlgWX::OnIntermolElStat)
	EVT_BUTTON(IDC_INTERMOL_TOT_ENE, InterMolDlgWX::OnIntermolTotEne)
	EVT_BUTTON(IDC_INTERMOL_SET_INT_COORD, InterMolDlgWX::OnIntermolSetIntCoord)
	EVT_BUTTON(IDC_INTERMOL_COMP_INT_COORD, InterMolDlgWX::OnIntermolCompIntCoord)
	EVT_BUTTON(IDC_MC_DOCK_RUN, InterMolDlgWX::OnMcDockRun)
	EVT_BUTTON(IDC_MC_DOCK_PLAYBACK_TRAJ, InterMolDlgWX::OnMcDockPlaybackTraj)
	EVT_BUTTON(IDC_MC_DOCK_PAUSE, InterMolDlgWX::OnMcDockPause)
	EVT_BUTTON(IDC_MC_DOCK_RESUME, InterMolDlgWX::OnMcDockResume)
	EVT_BUTTON(IDC_MC_DOCK_STOP, InterMolDlgWX::OnMcDockStop)
	EVT_CLOSE ( InterMolDlgWX::OnClose)
END_EVENT_TABLE()


void
InterMolDlgWX::OnInitDialog()
{
	wxTextCtrl* edit_ctrl;

    edit_ctrl = (wxTextCtrl*) FindWindow(IDC_MC_DOCK_TRAJ_FILE);
	edit_ctrl->SetValidator( StdStringValidator(&p_io_ag->traj_file_name) );
	
	DDX_Text_double(this,IDC_MC_DOCK_ANG_SCALE,p_mc_sim->ang_ratio,"%8.4f");
	DDX_Text_double(this,IDC_MC_DOCK_TR_SCALE,p_mc_sim->tr_ratio,"%8.4f");
	DDX_Text_long(this,IDC_MC_DOCK_NUM_STEPS,p_mc_sim->num_mc_steps);
	DDX_Text_double(this,IDC_MC_DOCK_DELAY,p_mc_sim->delay_time,"%9.3f");
	DDX_Text_long(this,IDC_MC_DOCK_SKIP_INIT,p_mc_sim->npt_begin);
	DDX_Text_long(this,IDC_MC_DOCK_SKIP_BETWEEN,p_mc_sim->npt_step);	
	DDX_Text_double(this,IDC_DEV_CONST_DIST,ptr_emp_mod->sigma_constr, "%9.3f");
	DDX_Text_double(this,IDC_WEIGHT_CONST_DIST, ptr_emp_mod->weight_constraints, "%3.1f");
	DDX_Text_double(this,IDC_PACK_DIST,ptr_emp_mod->pack_dist_com, "%9.3f");	
	DDX_Text_double(this,IDC_STDEV_PACK_DIST,ptr_emp_mod->sigma_pack_dist_com, "%9.3f");
    DDX_Text_double(this,IDC_WEIGHT_PACK_DIST,ptr_emp_mod->weight_pack_distance, "%3.1f");
	DDX_Text_double(this,IDC_PACK_ANGLE,ptr_emp_mod->pack_angle, "%3.1f");
	DDX_Text_double(this,IDC_STDEV_PACK_ANGLE, ptr_emp_mod->sigma_pack_angle, "%3.1f");
	DDX_Text_double(this,IDC_WEIGHT_PACK_ANGLE,ptr_emp_mod->weight_pack_angle, "%3.1f");
	DDX_Text_double(this,IDC_HEL_CONT,ptr_emp_mod->dist_contct_com, "%9.3f");
	DDX_Text_double(this,IDC_STDEV_HEL_CONT,ptr_emp_mod->sigma_dist_contct_com, "%9.3f"); 
	DDX_Text_double(this,IDC_WEIGHT_HEL_CONT,ptr_emp_mod->weight_num_contact, "%9.3f"); 
	DDX_Text_double(this,IDC_PACK_DEN,ptr_emp_mod->pack_dens, "%9.3f"); 
	DDX_Text_double(this,IDC_STDEV_PACK_DEN,ptr_emp_mod->sigma_pack_dens, "%9.3f");  
	DDX_Text_double(this,IDC_WEIGHT_PACK_DEN,ptr_emp_mod->weight_pack_dens, "%9.3f");
	DDX_Text_double(this,IDC_SYM,ptr_emp_mod->sigma_sym, "%9.3f");
	DDX_Text_double(this,IDC_WEIGHT_SYM,ptr_emp_mod->weight_sym, "%9.3f");
	DDX_Text_double(this,IDC_ANG_SA_FACE, ptr_emp_mod->face_up_bound, "%9.3f");
	DDX_Text_double(this,IDC_WEIGHT_SA_FACE, ptr_emp_mod->weight_sa, "%9.3f");
	DDX_Text_double(this,IDC_WEIGHT_VDW, ptr_emp_mod->weight_vdw, "%9.3f");
	DDX_Text_int(this,IDC_NUM_REPLICAS,ptr_im_mod->p_rex_sim->nreplicas);
	DDX_Text_int(this,IDC_REM_STEPS,ptr_im_mod->p_rex_sim->rem_steps);
	DDX_Text_double(this,IDC_HI_TEMP, ptr_im_mod->p_rex_sim->temperature_max, "%5.1f");
	DDX_Text_double(this,IDC_CUTOFF_VDW, ptr_emp_mod->dist_neighborhood, "%5.1f");

	TransferDataToWindow();
}

bool
InterMolDlgWX::TransferDataToWindow()
{
	wxCheckBox* check_box;
	
	check_box = (wxCheckBox*) FindWindow( IDC_INTERMOL_CHANGES_PK );
	if( ptr_im_mod->compute_pk) check_box->SetValue(true);
	else check_box->SetValue(false);

	check_box = (wxCheckBox*) FindWindow( IDC_MC_DOCK_ET_RATE2 );
	if( ptr_im_mod->calc_et_rate) check_box->SetValue(true);
	else check_box->SetValue(false);

	check_box = (wxCheckBox*) FindWindow( IDC_OPTIM_MD );
	if( ptr_im_mod->p_mc_sim->amber_flag) check_box->SetValue(true);
	else check_box->SetValue(false);	
	
	check_box = (wxCheckBox*) FindWindow( IDC_MC_DOCK_DONT_CALC_ENE );
	if( p_mc_sim->dont_calc_ene_flag) check_box->SetValue(true);
	else check_box->SetValue(false);
	
	check_box = (wxCheckBox*) FindWindow( IDC_MC_DOCK_SAVE_GIF );
	if( p_io_ag->save_image_seq_gif) check_box->SetValue(true);
	else check_box->SetValue(false);
	
	check_box = (wxCheckBox*) FindWindow( IDC_MC_DOCK_SAVE_PICT );
	if( p_io_ag->save_image_seq_pict) check_box->SetValue(true);
	else check_box->SetValue(false);
	
	TraceMolAgent* p_trace_mol_ag = p_mc_sim->GetTrajectoryTraceAgent(false);

	check_box = (wxCheckBox*) FindWindow( IDC_MC_DOCK_TRACE_ATOMS );
	if( p_trace_mol_ag != NULL && p_trace_mol_ag->traced_atoms.size() > 0) check_box->SetValue(true);
	else check_box->SetValue(false);

	check_box = (wxCheckBox*) FindWindow( IDC_MC_DOCK_DISCR_MOVES );
	if( p_mc_sim->IsDiscretizedMoves()) check_box->SetValue(true);
	else check_box->SetValue(false);

	check_box = (wxCheckBox*) FindWindow( IDC_SCORE);
	if( ptr_im_mod->empirical_flag) check_box->SetValue(true);
	else check_box->SetValue(false);	

	check_box = (wxCheckBox*) FindWindow( IDC_RUN_REM );
	if( ptr_im_mod->p_mc_sim->rex_flag ) check_box->SetValue(true);
	else check_box->SetValue(false);	

	check_box = (wxCheckBox*) FindWindow( IDC_VAR_TEMP );
	if( ptr_im_mod->p_rex_sim->vary_temperature_flag ) check_box->SetValue(true);
	else check_box->SetValue(false);
	
	int itemp;
	wxChoice* choice_ctrl;

	choice_ctrl = (wxChoice*) FindWindow( IDC_MC_DOCK_ELECTR_MODEL );
	itemp = ptr_im_mod->electr_model;
	choice_ctrl->SetSelection(itemp);

	wxTextCtrl* edit_mc_temp   = (wxTextCtrl*) FindWindow(IDC_MC_DOCK_TEMP);
	edit_mc_temp->SetValue( wxString::Format("%9.3f",p_mc_sim->GetTemperature()) );

	wxTextCtrl* edit_trace_at_ref   = (wxTextCtrl*) FindWindow(IDC_MC_DOCK_TRACE_AT_REF);
	
	MolSet* pmset = ptr_im_mod->GetMolSet();
	TraceMolAgent* tr_ag = p_mc_sim->GetTrajectoryTraceAgent(FALSE);
	wxString str;
	if(tr_ag != NULL)
	{
		AtomIteratorAtomGroup aitr(&tr_ag->traced_atoms);
		HaAtom* aptr;
		wxString str2;

		for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
			str2 = (aptr->GetRef()).c_str();
			str += str2;
			str += ";";
		}
	}
	edit_trace_at_ref->SetValue(str.c_str());

	return wxFrame::TransferDataToWindow();
}

bool
InterMolDlgWX::TransferDataFromWindow()
{
	wxCheckBox* check_box;

	check_box = (wxCheckBox*) FindWindow( IDC_INTERMOL_CHANGES_PK );
	ptr_im_mod->compute_pk = check_box->GetValue();

	check_box = (wxCheckBox*) FindWindow( IDC_MC_DOCK_ET_RATE2);
	ptr_im_mod->calc_et_rate = check_box->GetValue();

	check_box = (wxCheckBox*) FindWindow( IDC_OPTIM_MD );
	ptr_im_mod->p_mc_sim->amber_flag = check_box->GetValue();

	check_box = (wxCheckBox*) FindWindow( IDC_MC_DOCK_DONT_CALC_ENE );
	p_mc_sim->dont_calc_ene_flag = check_box->GetValue();

	check_box = (wxCheckBox*) FindWindow( IDC_MC_DOCK_SAVE_GIF );
	p_io_ag->save_image_seq_gif = check_box->GetValue();

	check_box = (wxCheckBox*) FindWindow( IDC_MC_DOCK_SAVE_PICT );
	p_io_ag->save_image_seq_pict = check_box->GetValue();

	check_box = (wxCheckBox*) FindWindow( IDC_MC_DOCK_DISCR_MOVES );
	if( check_box->GetValue()) p_mc_sim->SetDiscretizedMoves();

	check_box = (wxCheckBox*) FindWindow( IDC_SCORE );
	ptr_im_mod->empirical_flag = check_box->GetValue();

	check_box = (wxCheckBox*) FindWindow( IDC_RUN_REM );
	ptr_im_mod->p_mc_sim->rex_flag = check_box->GetValue();

	check_box = (wxCheckBox*) FindWindow( IDC_VAR_TEMP );
	ptr_im_mod->p_rex_sim->vary_temperature_flag = check_box->GetValue();

	int itemp;
	wxChoice* choice_ctrl;

	choice_ctrl = (wxChoice*) FindWindow( IDC_MC_DOCK_ELECTR_MODEL );
	itemp = choice_ctrl->GetSelection();
	ptr_im_mod->SetElectrModel(itemp);

	wxString str;
	double temp;

	wxTextCtrl* edit_mc_temp   = (wxTextCtrl*) FindWindow(IDC_MC_DOCK_TEMP);
	str = edit_mc_temp->GetValue();
	if(str.ToDouble(&temp)) p_mc_sim->SetTemperature(temp);

	TraceMolAgent* p_trace_mol_ag = NULL;

	check_box = (wxCheckBox*) FindWindow( IDC_MC_DOCK_TRACE_ATOMS );
	if( check_box->GetValue() )
	{
		p_trace_mol_ag = p_mc_sim->GetTrajectoryTraceAgent(TRUE);
		p_trace_mol_ag->traced_atoms.clear();
	}
	else
	{
		p_trace_mol_ag = p_mc_sim->GetTrajectoryTraceAgent(FALSE);
		p_mc_sim->DeleteTrajAnalAgent(p_trace_mol_ag);
	}

	if( check_box->GetValue() )
	{
		wxTextCtrl* edit_trace_at_ref   = (wxTextCtrl*) FindWindow(IDC_MC_DOCK_TRACE_AT_REF);

		MolSet* pmset = ptr_im_mod->GetMolSet();
		HaAtom* aptr;
	 
		wxString str;
		str = edit_trace_at_ref->GetValue();

		if( str.size() > 0 )
		{
			wxStringTokenizer tkz(str,";");
			int nat = 0;
			while ( tkz.HasMoreTokens())
			{
				wxString token = tkz.GetNextToken();
		
				aptr = pmset->GetAtomByRef( token.ToStdString() );
				if(aptr != NULL)
				{	
					p_trace_mol_ag->traced_atoms.InsertAtom(aptr);
					nat++;
				}
				PrintLog(" %d atoms are set to trace \n",nat);
			}
		}
	}

	return wxFrame::TransferDataFromWindow();
}


void InterMolDlgWX::OnIntermolInit(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	ptr_im_mod->Initialize();
}

void InterMolDlgWX::OnIntermolElStat(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	if( ptr_im_mod->electr_model == CHARGES_IN_FIELD_ELECTR ) 
	{
		ptr_im_mod->InitMolecularFields();
	}
	double el_inter_ene;
	el_inter_ene = ptr_im_mod->CalcElStaticInter();	

	PrintLog(" Electrostatic Intraction Energy %12.6f kcal/mol \n",el_inter_ene);
}


void InterMolDlgWX::OnIntermolTotEne(wxCommandEvent& event) 
{
	if( ptr_im_mod == NULL )
		return;
	TransferDataFromWindow();
	if( ptr_im_mod->electr_model == CHARGES_IN_FIELD_ELECTR ) 
	{
		ptr_im_mod->InitMolecularFields();
	}
	ptr_im_mod->CalcEffInterEne();
	PrintLog(" Total Intermolecular Energy      = %12.6f kcal/mol\n", ptr_im_mod->cur_intermol_ene);
	PrintLog(" Electrostatic Interaction Energy = %12.6f kcal/mol\n", ptr_im_mod->electr_inter_ene);
	PrintLog(" VdW Interaction Energy           = %12.6f kcal/mol\n", ptr_im_mod->vdw_inter_ene);
	if(ptr_im_mod->calc_et_rate)
	{
		PrintLog(" ET effective energy              = %12.6f kcal/mol\n", ptr_im_mod->add_eff_ene);
	}
	
}

void InterMolDlgWX::OnMcDockRun(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	ptr_im_mod->p_mc_sim->RunMCThread();	
}

void InterMolDlgWX::OnMcDockPlaybackTraj(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	p_mc_sim->AnalyzeTrajectory();	
}


void InterMolDlgWX::OnMcDockPause(wxCommandEvent& event) 
{
	p_mc_sim->PauseMC();	
}

void InterMolDlgWX::OnMcDockResume(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	p_mc_sim->ResumeMC();
}


void InterMolDlgWX::OnMcDockStop(wxCommandEvent& event) 
{
	p_mc_sim->StopMC();	
}

void InterMolDlgWX::OnIntermolSetIntCoord(wxCommandEvent& event) 
{
	wxTextCtrl* int_coord_1 = (wxTextCtrl*) FindWindow(IDC_INTERMOL_INT_COORD_1);
	wxTextCtrl* int_coord_2 = (wxTextCtrl*) FindWindow(IDC_INTERMOL_INT_COORD_2);

	wxString str;

	MolSet* pmset = ptr_im_mod->GetMolSet();
	HaMolecule* pMol1 = pmset->GetMolByIdx(0);
    HaMolecule* pMol2 = pmset->GetMolByIdx(1);

	if(pMol1 == NULL)
	{
		ErrorInMod("InterMolDlgWX::OnIntermolSetIntCoord()",
			       " No molecule 1");
		return;
	}
	str = int_coord_1->GetValue();
	pMol1->SetIntCoordFromStr(str.c_str());

	if(pMol2 == NULL)
	{
		ErrorInMod("InterMolDlgWX::OnIntermolSetIntCoord()",
			       " No molecule 2");
		return;
	}

	str = int_coord_2->GetValue();
	pMol2->SetIntCoordFromStr(str.c_str());

	pmset->RefreshAllViews(RFApply | RFRefresh);	
}

void InterMolDlgWX::OnIntermolCompIntCoord(wxCommandEvent& event) 
{
	wxTextCtrl* int_coord_1 = (wxTextCtrl*) FindWindow(IDC_INTERMOL_INT_COORD_1);
	wxTextCtrl* int_coord_2 = (wxTextCtrl*) FindWindow(IDC_INTERMOL_INT_COORD_2);

	wxString str;

	MolSet* pmset = ptr_im_mod->GetMolSet();
	HaMolecule* pMol1 = pmset->GetMolByIdx(0);
    HaMolecule* pMol2 = pmset->GetMolByIdx(1);

	if(pMol1 == NULL)
	{
		ErrorInMod("InterMolDlgWX::OnIntermolCompIntCoord()",
			       " No molecule 1");
		return;
	}
	
	double phi,cos_theta,psi;
	Vec3D trans;
	pMol1->GetPosEulerTrans( phi, cos_theta, psi, trans);
	str.Printf("%9.4f %9.4f %9.4f %9.4f %9.4f %9.4f",
		       trans[0], trans[1],trans[2],phi, cos_theta, psi);
	int_coord_1->SetValue(str);

	if(pMol2 == NULL)
	{
		ErrorInMod("InterMolDlgWX::OnIntermolCompIntCoord()",
			       " No molecule 2");
		return;
	}

	pMol2->GetPosEulerTrans( phi, cos_theta, psi, trans);
	str.Printf("%9.4f %9.4f %9.4f %9.4f %9.4f %9.4f",
		       trans[0], trans[1],trans[2],phi, cos_theta, psi);
	int_coord_2->SetValue(str);

}
/////////////////////////////////////////////////////////////////////
wxFieldPlaneView::wxFieldPlaneView(PlaneViewOfHaField3D* PlView, MolSet* _MolSet, wxWindow* parent,wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style, const wxString& name)
	: wxDialog(parent, id, title, pos, size, style, name)
		, m_MolSet(_MolSet)
{
	m_PlaneView=PlView;
	//this->SetExtraStyle(wxWS_EX_VALIDATE_RECURSIVELY);
	//this->SetExtraStyle(wxWS_EX_BLOCK_EVENTS);
	//dlg_open = TRUE;
	wxPlaneViewOfHaField3DDlgF(this,TRUE);
	
	wxCheckBox* CBHideZeroValues=(wxCheckBox*)FindWindow(IDCB_PV_HIDEZERO);
	CBHideZeroValues->SetValue(m_PlaneView->GetHideZeroValues());
	wxCheckBox* CBHide=(wxCheckBox*)FindWindow(IDCB_PV_HIDE);
	CBHide->SetValue(!m_PlaneView->IsDisplayed());

	OwnerOfView=0;
	
}
wxFieldPlaneView::~wxFieldPlaneView()
{
	//PrintLog("wxFieldPlaneView::~wxFieldPlaneView()");
	if(OwnerOfView==1)
	{
		m_MolSet->DeleteObject3D(m_PlaneView);
	}
	m_MolSet->RefreshAllViews(RFRefresh);
}
BEGIN_EVENT_TABLE(wxFieldPlaneView,wxDialog)
		EVT_COMMAND_SCROLL(IDS_PV_LEVEL,wxFieldPlaneView::OnScroll)
		EVT_CHOICE(IDC_PV_PLANE,wxFieldPlaneView::OnChoosePlane)
		EVT_TEXT_ENTER(IDT_PV_MIN,wxFieldPlaneView::OnEnter)
		EVT_TEXT_ENTER(IDT_PV_MAX,wxFieldPlaneView::OnEnter)
		EVT_SPINCTRL(IDSC_PV_LEVEL,wxFieldPlaneView::OnLevelChange)
		EVT_BUTTON(IDB_PV_CP1D1,wxFieldPlaneView::OnCopy1D1stVar)
		EVT_BUTTON(IDB_PV_CP1D2,wxFieldPlaneView::OnCopy1D2ndVar)
		EVT_BUTTON(IDB_PV_CP2D,wxFieldPlaneView::OnCopy2D)
		EVT_CHECKBOX(IDCB_PV_HIDEZERO, wxFieldPlaneView::OnHideZeroValue)
		EVT_CHECKBOX(IDCB_PV_HIDE, wxFieldPlaneView::OnHide)
		EVT_CLOSE(wxFieldPlaneView::OnClose)
END_EVENT_TABLE()
void wxFieldPlaneView::OnClose(wxCloseEvent& event )
{
	Destroy();
}
bool wxFieldPlaneView::TransferDataToWindow()
{
	wxDialog::TransferDataToWindow();
	if(m_PlaneView->GetHaField3D()==NULL)return true;
  
	float NMax,NMax2;
	int Plane=m_PlaneView->GetPlane();
  
	if(Plane==PlaneViewOfHaField3D::PlaneXY)
	{
		NMax=m_PlaneView->GetHaField3D()->GetNz();
		NMax2=m_PlaneView->GetHaField3D()->GetNy();
	}
	else if(Plane==PlaneViewOfHaField3D::PlaneZX)
	{
		NMax=m_PlaneView->GetHaField3D()->GetNy();
		NMax2=m_PlaneView->GetHaField3D()->GetNx();
	}
	else 
	{
		NMax=m_PlaneView->GetHaField3D()->GetNx();
		NMax2=m_PlaneView->GetHaField3D()->GetNy();
	}
	;
	wxString StrValueMin;
	StrValueMin<<m_PlaneView->GetValueMin();
	wxString StrValueMax;
	StrValueMax<<m_PlaneView->GetValueMax();
	((wxTextCtrl*)FindWindow(IDT_PV_MIN))->SetValue(StrValueMin);
	((wxTextCtrl*)FindWindow(IDT_PV_MAX))->SetValue(StrValueMax);
  //wxSpinCtrl
	wxSpinCtrl* SpinCtrl=(wxSpinCtrl*)FindWindow(IDSC_PV_LEVEL);
	SpinCtrl->SetRange(0,NMax-1);
	SpinCtrl->SetValue(m_PlaneView->GetLevel());
  
	SpinCtrl=(wxSpinCtrl*)FindWindow(IDSC_PV_LEVEL2);
	SpinCtrl->SetRange(0,NMax2-1);
  //wxChoice
	wxChoice* Choice=(wxChoice*)FindWindow(IDC_PV_PLANE);
	Choice->SetSelection(Plane);
	//Level Slider
	wxSlider* Slider=(wxSlider*)FindWindow(IDS_PV_LEVEL);
	Slider->SetValue(m_PlaneView->GetLevel());
	Slider->SetRange(0,NMax-1);
	
	return true;
}
bool wxFieldPlaneView::TransferDataFromWindow()
{
	wxDialog::TransferDataFromWindow();
  
//   wxString StrValueMin,StrValueMax;
//   StrValueMin=((wxTextCtrl*)FindWindow(IDT_PNP_MIN))->GetValue();
//   StrValueMax=((wxTextCtrl*)FindWindow(IDT_PNP_MAX))->GetValue();
//   double dtmp;
//   StrValueMin.ToDouble(&dtmp);
//   ValueMin=dtmp;
//   StrValueMax.ToDouble(&dtmp);
//   ValueMax=dtmp;
//   Plane=((wxChoice*)FindWindow(IDC_PNP_PLANE))->GetSelection();
//   ColoringSheme=((wxChoice*)FindWindow(IDC_PNP_COLSCH))->GetSelection();
	return true;
}
void wxFieldPlaneView::OnCopy1D1stVar( wxCommandEvent& event )
{
	int Plane=m_PlaneView->GetPlane();
  
	int Level=m_PlaneView->GetLevel();
	wxSpinCtrl* SpinCtrl1=(wxSpinCtrl*)FindWindow(IDSC_PV_LEVEL);
	int lev1=SpinCtrl1->GetValue();
	wxSpinCtrl* SpinCtrl2=(wxSpinCtrl*)FindWindow(IDSC_PV_LEVEL2);
	int lev2=SpinCtrl2->GetValue();
	
	int i,j,k,gpnt;
	int GS_X=m_PlaneView->GetHaField3D()->GetNx();
	int GS_Y=m_PlaneView->GetHaField3D()->GetNy();
	int GS_Z=m_PlaneView->GetHaField3D()->GetNz();
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_X*GS_Y*GS_Z;
	float *v=m_PlaneView->GetHaField3D()->GetFieldPtr();
	wxString TextOut("");
  
	if(Plane==m_PlaneView->PlaneXY)
	{
		printf("OnCopy1D Plane::XY [0 .. %d, %d, %d]\n",GS_X-1,lev2,Level);
		for(k=0;k<GS_X;k++)
		{
			gpnt=k+lev2*GS_X+Level*GS_XY;
			TextOut+=wxString::Format("%.6e\n",v[gpnt]);
      //TextOut<<v[gpnt]<<"\n";
		}
		//gpnt=GS_X-1+lev2*GS_X+Level*GS_XY;
		//TextOut+=wxString::Format("%.6e",v[gpnt]);
	}
	else if(Plane==m_PlaneView->PlaneZX)
	{
		printf("OnCopy1D Plane::ZX [%d, %d, 0 .. %d]\n",lev2,Level,GS_Z-1);
		for(k=0;k<GS_Z;k++)
		{
			gpnt=lev2+Level*GS_X+k*GS_XY;
			TextOut+=wxString::Format("%.6e\n",v[gpnt]);
      //TextOut<<v[gpnt]<<"\n";
		}
		//gpnt=lev2+Level*GS_X+(GS_Z-1)*GS_XY;
		//TextOut+=wxString::Format("%.6e",v[gpnt]);
	}
	else 
	{
		printf("OnCopy1D Plane::YZ [%d, 0 .. %d, %d]\n",Level,GS_Y-1,lev2);
		for(k=0;k<GS_Y;k++)
		{
			gpnt=Level+k*GS_X+lev2*GS_XY;
			TextOut+=wxString::Format("%.6e\n",v[gpnt]);
      //TextOut<<v[gpnt]<<"\n";
		}
		//gpnt=Level+(GS_Y-1)*GS_X+lev2*GS_XY;
		//TextOut+=wxString::Format("%.6e",v[gpnt]);
	}
	if (wxTheClipboard->Open())
	{
		wxTheClipboard->SetData( new wxTextDataObject(TextOut) );
		wxTheClipboard->Close();
	}
}
void wxFieldPlaneView::OnCopy1D2ndVar( wxCommandEvent& event )
{
	int Plane=m_PlaneView->GetPlane();
  
	int Level=m_PlaneView->GetLevel();
	wxSpinCtrl* SpinCtrl1=(wxSpinCtrl*)FindWindow(IDSC_PV_LEVEL);
	int lev1=SpinCtrl1->GetValue();
	wxSpinCtrl* SpinCtrl2=(wxSpinCtrl*)FindWindow(IDSC_PV_LEVEL2);
	int lev2=SpinCtrl2->GetValue();
	
	int i,j,k,gpnt;
	int GS_X=m_PlaneView->GetHaField3D()->GetNx();
	int GS_Y=m_PlaneView->GetHaField3D()->GetNy();
	int GS_Z=m_PlaneView->GetHaField3D()->GetNz();
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_X*GS_Y*GS_Z;
	float *v=m_PlaneView->GetHaField3D()->GetFieldPtr();
	
	
	wxString TextOut("");
  
	if(Plane==m_PlaneView->PlaneXY)
	{
		PrintLog("OnCopy1D Plane::XY [%d, 0 .. %d, %d]\n",lev1,GS_Y-1,Level);
		for(k=0;k<GS_Y;k++)
		{
			gpnt=lev1+k*GS_X+Level*GS_XY;
			TextOut+=wxString::Format("%.6e\n",v[gpnt]);
      //TextOut<<v[gpnt]<<"\n";
		}
		//gpnt=lev1+(GS_Y-1)*GS_X+Level*GS_XY;
		//TextOut+=wxString::Format("%.6e",v[gpnt]);
	}
	else if(Plane==m_PlaneView->PlaneZX)
	{
		PrintLog("OnCopy1D Plane::ZX [0 .. %d, %d, %d]\n",GS_X-1,lev1,Level);
		for(k=0;k<GS_X;k++)
		{
			gpnt=k+Level*GS_X+lev1*GS_XY;
			TextOut+=wxString::Format("%.6e\n",v[gpnt]);
      //TextOut<<v[gpnt]<<"\n";
		}
		//gpnt=GS_X-1+Level*GS_X+lev1*GS_XY;
		//TextOut+=wxString::Format("%.6e",v[gpnt]);
	}
	else 
	{
		PrintLog("OnCopy1D Plane::YZ [%d, %d, 0 .. %d]\n",Level,lev1,GS_Z-1);
		for(k=0;k<GS_Z;k++)
		{
			gpnt=Level+lev1*GS_X+k*GS_XY;
			TextOut+=wxString::Format("%.6e\n",v[gpnt]);
      //TextOut<<v[gpnt]<<"\n";
		}
		//gpnt=Level+lev1*GS_X+(GS_Z-1)*GS_XY;
		//TextOut+=wxString::Format("%.6e",v[gpnt]);
	}
	if (wxTheClipboard->Open())
	{
		wxTheClipboard->SetData( new wxTextDataObject(TextOut) );
		wxTheClipboard->Close();
	}
}
void wxFieldPlaneView::OnCopy2D( wxCommandEvent& event )
{
	int Plane=m_PlaneView->GetPlane();
  
	int lev1=m_PlaneView->GetLevel();
	int i,j,k,gpnt;
	int GS_X=m_PlaneView->GetHaField3D()->GetNx();
	int GS_Y=m_PlaneView->GetHaField3D()->GetNy();
	int GS_Z=m_PlaneView->GetHaField3D()->GetNz();
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_X*GS_Y*GS_Z;
	float *v=m_PlaneView->GetHaField3D()->GetFieldPtr();
	wxString TextOut("");
  
	if(Plane==m_PlaneView->PlaneXY)
	{
		PrintLog("OnCopy2D Plane::XY Z=%d\n",lev1);
		for(j=0;j<GS_Y-1;j++)
		{
			for(i=0;i<GS_X-1;i++)
			{
				gpnt=i+j*GS_X+lev1*GS_XY;
				TextOut+=wxString::Format("%.6e\t",v[gpnt]);
			}
			gpnt=GS_X-1+j*GS_X+lev1*GS_XY;
			TextOut+=wxString::Format("%.6e\n",v[gpnt]);
		}
		j=GS_Y-1;
		for(i=0;i<GS_X-1;i++)
		{
			gpnt=i+j*GS_X+lev1*GS_XY;
			TextOut+=wxString::Format("%.6e\t",v[gpnt]);
		}
		gpnt=GS_X-1+j*GS_X+lev1*GS_XY;
		TextOut+=wxString::Format("%.6e",v[gpnt]);
	}
	else if(Plane==m_PlaneView->PlaneZX)
	{
		PrintLog("OnCopy2D Plane::ZX %d\n",lev1);
		for(i=0;i<GS_X-1;i++)
		{
			for(k=0;k<GS_Z-1;k++)
			{
				gpnt=i+lev1*GS_X+k*GS_XY;
				TextOut+=wxString::Format("%.6e\t",v[gpnt]);
			}
			gpnt=i+lev1*GS_X+(GS_Z-1)*GS_XY;
			TextOut+=wxString::Format("%.6e\n",v[gpnt]);
		}
		i=GS_X-1;
		for(k=0;k<GS_Z-1;k++)
		{
			gpnt=i+lev1*GS_X+k*GS_XY;
			TextOut+=wxString::Format("%.6e\t",v[gpnt]);
		}
		gpnt=i+lev1*GS_X+(GS_Z-1)*GS_XY;
		TextOut+=wxString::Format("%.6e\n",v[gpnt]);
	}
	else 
	{
		PrintLog("OnCopy2D Plane::YZ %d\n",lev1);
		for(k=0;k<GS_Z-1;k++)
		{
			for(j=0;j<GS_Y-1;j++)
			{
				gpnt=lev1+j*GS_X+k*GS_XY;
				TextOut+=wxString::Format("%.6e\t",v[gpnt]);
			}
			gpnt=lev1+(GS_Y-1)*GS_X+k*GS_XY;
			TextOut+=wxString::Format("%.6e\n",v[gpnt]);
		}
		k=GS_Z-1;
		for(j=0;j<GS_Y-1;j++)
		{
			gpnt=lev1+j*GS_X+k*GS_XY;
			TextOut+=wxString::Format("%.6e\t",v[gpnt]);
		}
		gpnt=lev1+(GS_Y-1)*GS_X+k*GS_XY;
		TextOut+=wxString::Format("%.6e\n",v[gpnt]);
	}
	if (wxTheClipboard->Open())
	{
		wxTheClipboard->SetData( new wxTextDataObject(TextOut) );
		wxTheClipboard->Close();
	}
}
void wxFieldPlaneView::OnEnter(wxCommandEvent& event)
{
	wxString StrValueMin,StrValueMax;
	StrValueMin=((wxTextCtrl*)FindWindow(IDT_PV_MIN))->GetValue();
	StrValueMax=((wxTextCtrl*)FindWindow(IDT_PV_MAX))->GetValue();
	double dtmp1,dtmp2;
	StrValueMin.ToDouble(&dtmp1);
	StrValueMax.ToDouble(&dtmp2);
	m_PlaneView->SetMinMax(dtmp1,dtmp2);
	m_MolSet->RefreshAllViews(RFRefresh);
}
void wxFieldPlaneView::OnChoosePlane( wxCommandEvent& event )
{
	m_PlaneView->SetPlane(event.GetInt());
	TransferDataToWindow();
	m_MolSet->RefreshAllViews(RFRefresh);
}
void wxFieldPlaneView::OnLevelChange( wxSpinEvent& event )
{
	wxSpinCtrl* SpinCtrl=(wxSpinCtrl*)FindWindow(IDSC_PV_LEVEL);
	int l=SpinCtrl->GetValue();
	//if(l<0)l=0;
	//if(l>NMax-1)l=NMax-1;
	m_PlaneView->SetLevel(l);
	TransferDataToWindow();
	m_MolSet->RefreshAllViews(RFRefresh);
}
void wxFieldPlaneView::OnScroll( wxScrollEvent& event )
{
	wxSlider* Slider=(wxSlider*)FindWindow(IDS_PV_LEVEL);
	m_PlaneView->SetLevel(Slider->GetValue());
	TransferDataToWindow();
	m_MolSet->RefreshAllViews(RFRefresh);
}
void wxFieldPlaneView::OnHideZeroValue( wxCommandEvent& event )
{
	wxCheckBox* CB=(wxCheckBox*)FindWindow(IDCB_PV_HIDEZERO);
	m_PlaneView->SetHideZeroValues(CB->GetValue());
	m_MolSet->RefreshAllViews(RFRefresh);
}
void wxFieldPlaneView::OnHide( wxCommandEvent& event )
{
	wxCheckBox* CB=(wxCheckBox*)FindWindow(IDCB_PV_HIDE);
	m_PlaneView->SetDisplayed(!CB->GetValue());
	m_MolSet->RefreshAllViews(RFRefresh);
}
void CreatewxFieldPlaneView(PlaneViewOfHaField3D* PlView, MolSet* _MolSet,const char *title,int OwnerOfView)
{
	wxString Name(title);
	wxFieldPlaneView *PV=new wxFieldPlaneView(PlView,_MolSet,(wxWindow*)GetHaMainFrameWX(),-1,Name);
	PV->OwnerOfView=OwnerOfView;
	//wxFieldPlaneView *PV=new wxFieldPlaneView(PlView,_MolSet,parent,-1,Name);
	PV->Show(true);
}


