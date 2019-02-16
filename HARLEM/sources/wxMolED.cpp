/*! \file wxMolED.cpp

    Dialogs for Essential Dynamics and clustering Analysis
 
    \author Igor Kurnikov  
    \date 2011-
*/

#define WX_MOL_ED_CPP

#include <mpi.h>

#include "hastl.h"
#include <stdexcept>

#include "wx/wx.h"
#include "wx/notebook.h"
#include "wx/valgen.h"
#include "wx/filename.h"
#include "wx/app.h"

#include "ctrl_wx.h"
#include "ha_wx_aux_1.h"

#include "canvas3d.h"
#include "dialogs_wx_1.h"
#include "hamainframe_wx.h"
#include "harlemapp_wx.h"

#include "hachart_wx.h"
#include "hachart.h"

#include "haatgroup.h"
#include "hacompmod.h"
#include "haproteined.h"
#include "hamolmech.h"
#include "mm_traj_anal.h"
#include "haproteined.h"


#if defined(__WXMSW__)
// undefined MFC MACROS conflicting WX functions
#ifdef DrawText
#undef DrawText
#endif

#ifdef StartDoc
#undef StartDoc
#endif

#ifdef GetCharWidth
#undef GetCharWidth
#endif

#ifdef FindWindow
#undef FindWindow
#endif

#endif

#include "wxMolED.h"
#include <wx/filedlg.h>

#include "ha_wx_ed_wdr.h"

#include "hamolview.h"

int wxMolED::dlg_open = FALSE;

wxMolED::wxMolED( CollectCrdAnalMod* p_cluster_mod_new , wxWindow* parent)
:wxFrame( parent, -1, "Clustering Analysis" )
{
	this->SetExtraStyle(wxWS_EX_VALIDATE_RECURSIVELY);
	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);

	wxMenuBar* p_menu_bar = collect_crd_menu();
    this->SetMenuBar(p_menu_bar);    

	dlg_open = TRUE;
	cluster_anal_dlg(this, TRUE);
//	EssentialDynamicsDlg(this,TRUE);
	p_ccrd_mod = p_cluster_mod_new;
	pmset = p_ccrd_mod->GetMolSet();
	p_anim_view = NULL;
	this->parent = parent;
	shift_vec_val = 5.0;
	OnInitDialog();
}

wxMolED::~wxMolED()
{

}

BEGIN_EVENT_TABLE(wxMolED, wxFrame)
    EVT_MENU(IDM_ED_TEST_GRAPH, wxMolED::OnTestGraph)
	EVT_BUTTON(IDC_ED_CHOOSE_MD_TRAJ,    wxMolED::OnChooseMDTraj)
	EVT_BUTTON(IDC_ED_CHOOSE_PLATO_TRAJ, wxMolED::OnChoosePlatoTraj)
	EVT_BUTTON(IDC_ED_CHOOSE_ATGRP,   wxMolED::OnChooseAtGrpToAnal)
	EVT_BUTTON(IDC_ED_GEN_PLATO_TRAJ, wxMolED::OnGenPlatoTraj)
	EVT_BUTTON(IDC_ED_SAVE_PLATO_INPUT,    wxMolED::OnSavePlatoInputFiles )
	EVT_BUTTON(IDC_ED_RUN_PLATO,           wxMolED::OnRunPlato )
	EVT_BUTTON(IDC_ED_CHOOSE_PLATO_OUTPUT, wxMolED::OnChoosePlatoOutputFile )
	EVT_BUTTON(IDC_ED_LOAD_PLATO_OUTPUT,   wxMolED::OnLoadPlatoOutputFile )
	EVT_BUTTON(IDC_ED_RUN_VEC_ANIMATION,   wxMolED::OnRunVecAnimation )
	EVT_BUTTON(IDC_ED_STOP_VEC_ANIMATION,  wxMolED::OnStopVecAnimation )
	EVT_BUTTON(IDC_ED_SHIFT_ALONG_VEC,     wxMolED::OnShiftAlongVec )
	EVT_BUTTON(IDC_ED_PROJ_TRAJ,           wxMolED::OnProjTraj )
	EVT_BUTTON(IDC_ED_PLOT_TIME_PROJ,      wxMolED::OnPlotTimeProj )
	EVT_BUTTON(IDC_ED_PLOT_PROJ_2_VS_1,    wxMolED::OnPlotProj2vs1 )
	EVT_BUTTON(ID_BROWSE_PDB, wxMolED::OnTopFileBrowse)
	EVT_BUTTON(ID_BROWSE_EV, wxMolED::OnEigFileBrowse)
	EVT_BUTTON(ID_BROWSE_PROJ, wxMolED::OnPpjFileBrowse)
	EVT_BUTTON(ID_BUTTON_NEXT, wxMolED::OnNextButtonClick)
	EVT_CLOSE   (wxMolED::OnClose)
END_EVENT_TABLE()

void wxMolED::OnInitDialog()
{
	wxTextCtrl* edit_ctrl;

	edit_ctrl = (wxTextCtrl*) FindWindow( IDC_ED_MD_TRAJ_FILE ); 
	edit_ctrl->SetValidator( StdStringValidator(&p_ccrd_mod->p_md_traj->CrdFileName) );

	edit_ctrl = (wxTextCtrl*) FindWindow( IDC_ED_TRAJ_FILE_PLATO ); 
	edit_ctrl->SetValidator( StdStringValidator(&p_ccrd_mod->md_traj_fname_plato) );

	edit_ctrl = (wxTextCtrl*) FindWindow( IDC_ED_PLATO_INPUT_FILE ); 
	edit_ctrl->SetValidator( StdStringValidator(&p_ccrd_mod->plato_input_fname) );

	edit_ctrl = (wxTextCtrl*) FindWindow( IDC_ED_PLATO_OUTPUT_FILE ); 
	edit_ctrl->SetValidator( StdStringValidator(&p_ccrd_mod->plato_output_fname) );

	edit_ctrl = (wxTextCtrl*) FindWindow( IDC_ED_PLATO_RUN_FILE ); 
	edit_ctrl->SetValidator( StdStringValidator(&p_ccrd_mod->plato_run_fname) );

	edit_ctrl = (wxTextCtrl*) FindWindow( IDC_ED_TIME_PROJ_FILE ); 
	edit_ctrl->SetValidator( StdStringValidator(&p_ccrd_mod->time_proj_fname) );

	DDX_Text_int (this,IDC_ED_N_CLUST,     p_ccrd_mod->num_clusters);
	DDX_Text_int (this,IDC_ED_N_EIG_VEC,   p_ccrd_mod->num_eig_vec);
	DDX_Text_int (this,IDC_ED_N_TIME_PROJ, p_ccrd_mod->num_time_proj);

	DDX_Text_double(this,IDC_ED_SHIFT_VAL, shift_vec_val, "%9.3f" );

	wxCheckBox* check_box = (wxCheckBox*) FindWindow( IDC_ED_COMP_SIM_MATRIX );
	check_box->SetValidator( IntCheckBoxValidator(&p_ccrd_mod->sim_matrix_flag) );

	DDX_Choice_HaEnum(this,IDC_ED_PLATO_ANAL_TYPE,&p_ccrd_mod->anal_type );

	TransferDataToWindow();
}

bool wxMolED::TransferDataToWindow()
{
	wxString str;
	std::string atgrp_name = p_ccrd_mod->GetActiveAtomGroupName();
	wxTextCtrl* edit_ctrl = (wxTextCtrl*) FindWindow( IDC_ED_ATGRP_TO_ANAL );
	edit_ctrl->SetValue( atgrp_name.c_str() );

	wxListBox* eigval_list= (wxListBox*) FindWindow(IDC_ED_EIGVEC_LIST);
	eigval_list->Clear();
	int nv = p_ccrd_mod->eigen_vals.size();
	int i;
	for(i=1; i <= nv; i++)
	{
		str.Printf("%4d :%16.9f", i, p_ccrd_mod->eigen_vals[i-1]);
		eigval_list->Append(str);
	}
	wxListBox* proj_list= (wxListBox*) FindWindow(IDC_ED_TIME_PROJ_LIST);
	proj_list->Clear();

	int np = p_ccrd_mod->time_projections.size();
	for(i = 1; i <= np; i++)
	{
		str.Printf("%4d ", i);
		proj_list->Append(str);
	}

	return wxFrame::TransferDataToWindow();
}

bool wxMolED::TransferDataFromWindow()
{
	wxListBox* eigval_list= (wxListBox*) FindWindow(IDC_ED_EIGVEC_LIST);
	int nsel = eigval_list->GetSelections(sel_arr);
	
	return wxFrame::TransferDataFromWindow();
}

void wxMolED::OnClose(wxCloseEvent& event)
{
	dlg_open = FALSE;
	event.Skip();
}


void wxMolED::OnTestGraph(wxCommandEvent& event)
{
#if !defined HARLEM_WX_PLPLOT_NO
	HaChartFrame *frame = new HaChartFrame( _T( "wxPLplot demo" ) );
	frame->PlotTest();
    frame->Show( true );
#endif
}

void wxMolED::OnChooseMDTraj(wxCommandEvent& event)
{
	TransferDataFromWindow();

	wxString mdcrd_file_name = ::wxFileSelector("Choose MD Trajectory Coordinates File",
		::wxGetCwd(),p_ccrd_mod->p_md_traj->CrdFileName.c_str(),
		"mdcrd","*.mdcrd");

	if(!mdcrd_file_name.empty() )
	{
		wxFileName scr_fname(mdcrd_file_name);
		wxString cur_dir = ::wxGetCwd();
		PrintLog("Current Directory: %s \n",cur_dir.ToStdString().c_str());
        scr_fname.MakeRelativeTo(cur_dir);
		wxString mod_mdcrd_file_name = scr_fname.GetFullPath();
		p_ccrd_mod->p_md_traj->CrdFileName = mod_mdcrd_file_name.c_str();
		TransferDataToWindow();
	}
	else
	{
		PrintLog("Invalid MD Trajectory Coordinates File Name");
	}
	this->Raise();
}

void wxMolED::OnChoosePlatoTraj(wxCommandEvent& event)
{
	TransferDataFromWindow();

	wxString plato_file_name = ::wxFileSelector("Choose Trajectory File in PLATO format",
		::wxGetCwd(),p_ccrd_mod->md_traj_fname_plato.c_str(),
		"plato","*.plato");

	if(!plato_file_name.empty() )
	{
		wxFileName scr_fname(plato_file_name);
		wxString cur_dir = ::wxGetCwd();
		PrintLog("Current Directory: %s \n",cur_dir.ToStdString().c_str());
        scr_fname.MakeRelativeTo(cur_dir);
		wxString mod_plato_file_name = scr_fname.GetFullPath();
		p_ccrd_mod->md_traj_fname_plato = mod_plato_file_name.c_str();
		TransferDataToWindow();
	}
	else
	{
		PrintLog("Invalid PLATO Trajectory File Name");
	}
	this->Raise();
}

void wxMolED::OnGenPlatoTraj(wxCommandEvent& event)
{
	TransferDataFromWindow();
	p_ccrd_mod->ConvertMDTrajToPlato();
}


void wxMolED::OnChooseAtGrpToAnal(wxCommandEvent& event)
{
	EditGroupsDlg* edit_grp_dlg = new EditGroupsDlg( pmset, 2, NULL );
	edit_grp_dlg->ShowModal();
	AtomGroup* p_at_arr = edit_grp_dlg->GetSelGroup();
	if(p_at_arr == NULL) return;

	p_ccrd_mod->SetActiveAtomGroup( p_at_arr->GetID() );

	TransferDataToWindow();
	this->Raise();
}


void wxMolED::OnTopFileBrowse(wxCommandEvent& event)
{
	wxString cwd = ::wxGetCwd();
	wxFileDialog* browser = new wxFileDialog(parent, "Choose a PDB file", cwd, "", "*.pdb", wxFD_OPEN, wxDefaultPosition);
	browser->ShowModal();
	topFileName = browser->GetFilename();
	topFileDir = browser->GetDirectory();
	if (topFileName.IsEmpty())
		exit(1);
	// After getting the file name we want to store it in our little text box
	wxTextCtrl* pdbFileName = (wxTextCtrl*) FindWindow(ID_PDB_TEXTCTRL);
	pdbFileName->SetValue(topFileName);
	
}

void wxMolED::OnEigFileBrowse(wxCommandEvent& event)
{
	wxString cwd = ::wxGetCwd();
	wxFileDialog* browser = new wxFileDialog(parent, "Choose EigenValue File", cwd, "", "*.pev", wxFD_OPEN, wxDefaultPosition);
	browser->ShowModal();
	eigFileName = browser->GetFilename();
	eigFileDir = browser->GetDirectory();
	if (eigFileName.IsEmpty())
		exit(1);
	// After getting the file name we want to store it in our little text box
	wxTextCtrl* eigenvalFileName = (wxTextCtrl*) FindWindow(ID_EV_TEXTCTRL);
	eigenvalFileName->SetValue(eigFileName);
}

void wxMolED::OnPpjFileBrowse(wxCommandEvent& event)
{
	wxString cwd = ::wxGetCwd();
	wxFileDialog* browser = new wxFileDialog(parent, "Choose EigenValue File", cwd, "", "*.ppj", wxFD_OPEN, wxDefaultPosition);
	browser->ShowModal();
	ppjFileName = browser->GetFilename();
	ppjFileDir = browser->GetDirectory();
	if (ppjFileName.IsEmpty())
		exit(1);
	// After getting the file name we want to store it in our little text box
	wxTextCtrl* projFileName = (wxTextCtrl*) FindWindow(ID_PPJ_TEXTCTRL);
	projFileName->SetValue(ppjFileName);
}

void wxMolED::OnNextButtonClick(wxCommandEvent& event)
{
	// Now load the PDB file that was selected
	HaMainFrameWX* frame_main = GetHaMainFrameWX();
    
	wxFileName fname_full;
    fname_full.Assign(topFileDir,topFileName);
	wxString fullPathName = fname_full.GetFullPath();
	PrintLog("File Name is %s", fullPathName.ToStdString().c_str());

	if( pmset == NULL || pmset->canvas_wx == NULL)
    {
        pmset = new HaMolSet();
        pmset->canvas_wx = frame_main->CreateMolView(pmset);        
    }
	pmset->FetchFile(FormatPDB, fullPathName.c_str()); 

	pmset->canvas_wx->mol_view->InitialTransform();
    pmset->canvas_wx->mol_view->DefaultRepresentation();
    pmset->RefreshAllViews(RFRefresh | RFColour | RFApply);

	// Now close this dialog box and show the next dialog
	
	fname_full.Assign(eigFileDir, eigFileName);
	fullPathName = fname_full.GetFullPath();

}

void wxMolED::OnSavePlatoInputFiles(wxCommandEvent& event)
{
	TransferDataFromWindow();
	p_ccrd_mod->SavePlatoInputFiles();

}

void wxMolED::OnRunPlato(wxCommandEvent& event)
{
	wxBusyCursor wait;
	TransferDataFromWindow();
	p_ccrd_mod->RunPlato();
}

void wxMolED::OnChoosePlatoOutputFile(wxCommandEvent& event)
{
	TransferDataFromWindow();

	wxString plato_output_file_name = ::wxFileSelector("Choose PLATO output file",
		::wxGetCwd(),p_ccrd_mod->plato_output_fname.c_str(),
		"xml","*.xml");

	if(!plato_output_file_name.empty() )
	{
		wxFileName scr_fname(plato_output_file_name);
		wxString cur_dir = ::wxGetCwd();
		PrintLog("Current Directory: %s \n",cur_dir.ToStdString().c_str());
        scr_fname.MakeRelativeTo(cur_dir);
		wxString mod_plato_output_file_name = scr_fname.GetFullPath();
		p_ccrd_mod->plato_output_fname = mod_plato_output_file_name.c_str();
		TransferDataToWindow();
	}
	else
	{
		PrintLog("Invalid PLATO Output File Name");
	}
	this->Raise();
}

void wxMolED::OnLoadPlatoOutputFile(wxCommandEvent& event)
{
	TransferDataFromWindow();
	p_ccrd_mod->LoadPlatoOutputFile();
	TransferDataToWindow();
}

void wxMolED::OnRunVecAnimation(wxCommandEvent& event)
{
	TransferDataFromWindow();
	try
	{
		if( CurMolView == NULL ) throw std::runtime_error("Current Molecular View is not defined");
		p_anim_view = CurMolView;
		
		if(sel_arr.size() == 0 ) throw std::runtime_error("Empty eigen vector list or No selection");
		
		int isel = sel_arr[0];
		if( isel >= p_ccrd_mod->eigen_vecs.size() ) throw std::runtime_error("Index of selected eigen value is beyound bounds of eigenvalue array");
		p_anim_view->AnimateEigenVector( p_ccrd_mod->eigen_vecs[isel], p_ccrd_mod->GetActiveAtomGroup() );
	}
	catch( std::exception& ex )
	{
		PrintLog(" Error in wxMolED::OnRunVecAnimation() \n");
		PrintLog("%s\n",ex.what());
		p_anim_view = NULL;
	}
}

void wxMolED::OnStopVecAnimation(wxCommandEvent& event)
{
	if( p_anim_view != NULL )
	{
		p_anim_view->StopAnimation();
		p_anim_view = NULL;
	}
	else
	{
		PrintLog("Error in wxMolED::OnStopVecAnimation() \n");
		PrintLog("Molecular View to display animation is not defined \n");
	}
}

void wxMolED::OnProjTraj(wxCommandEvent& event)
{
	TransferDataFromWindow();
	p_ccrd_mod->CalcTimeProj();
	p_ccrd_mod->SaveTimeProjFile();
	TransferDataToWindow();
}


void wxMolED::OnShiftAlongVec(wxCommandEvent& event)
{
	TransferDataFromWindow();
	try
	{
		if(sel_arr.size() == 0 ) throw std::runtime_error("Empty eigen vector list or No selection");
		int isel = sel_arr[0];
		if( isel >= p_ccrd_mod->eigen_vecs.size() ) throw std::runtime_error("Index of selected eigen value is beyond bounds of eigenvalue array");

		p_ccrd_mod->ShiftAlongEigenVec(isel,shift_vec_val);
		pmset->RefreshAllViews(RFApply);
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in wxMolED::OnShiftAlongVec() \n");
		PrintLog(" %s\n",ex.what());
	}
}

void wxMolED::OnPlotTimeProj(wxCommandEvent& event)
{
#if !defined HARLEM_WX_PLPLOT_NO
	wxListBox* proj_list= (wxListBox*) FindWindow(IDC_ED_TIME_PROJ_LIST);

	wxArrayInt idx_sel;
    proj_list->GetSelections(idx_sel);
	int np = idx_sel.size();

	if( np == 0 )
	{
		wxLogMessage(" At least one time projection should be selected ");
		return;
	}
	
	HaChartFrame* plt_frame = new HaChartFrame("Collect Coord Time Projections");
	HaChartPanel* pplt = plt_frame->GetPlot();	
	HaChart2D* p2d = pplt->AddChart2D();

	int ip;
	for( ip = 0; ip < np; ip++)
	{
		int idx_p = idx_sel[ip];
		if( idx_p >= 0 && idx_p < p_ccrd_mod->time_projections.size() )
		{
			p2d->AddYData( p_ccrd_mod->time_projections[idx_p] );
		}
	}
	plt_frame->Plot();
	plt_frame->Show();
#endif
}

void wxMolED::OnPlotProj2vs1(wxCommandEvent& event)
{
#if !defined HARLEM_WX_PLPLOT_NO
    wxListBox* proj_list= (wxListBox*) FindWindow(IDC_ED_TIME_PROJ_LIST);

	wxArrayInt idx_sel;
    proj_list->GetSelections(idx_sel);
	int np = idx_sel.size();

	if( np != 2 && np != 3) 
	{
		wxLogMessage(" Two or three time projections should be selected ");
		return;
	}

	int idx_p1 = idx_sel[0];
	int idx_p2 = idx_sel[1];

	if( idx_p1 < 0 || idx_p1 >= p_ccrd_mod->time_projections.size() ||
		idx_p2 < 0 || idx_p2 >= p_ccrd_mod->time_projections.size() )
	{
		PrintLog(" Invalid projection indexes ");
		return;
	}

	int idx_p3 = -1;
	if( np > 2 ) 
	{
		idx_p3 = idx_sel[2];
		if( idx_p3 < 0 || idx_p3 >= p_ccrd_mod->time_projections.size())
		{
			PrintLog(" Invalid projection indexes ");
			return;
		}
	}

	HaChartFrame* plt_frame = new HaChartFrame("Collect Coord Time Projections");
	HaChartPanel* pplt = plt_frame->GetPlot();
	if( np == 2 )
	{
		HaChart2D* p2d = pplt->AddChart2D();
		p2d->AddXYData( p_ccrd_mod->time_projections[idx_p1], p_ccrd_mod->time_projections[idx_p2] );
	}
	else if( np == 3 )
	{
		HaChart3D* p3d = pplt->AddChart3D();
		p3d->AddXYZData( p_ccrd_mod->time_projections[idx_p1], p_ccrd_mod->time_projections[idx_p2],p_ccrd_mod->time_projections[idx_p3] );
	}

	plt_frame->Plot();
	plt_frame->Show();
#endif
}
