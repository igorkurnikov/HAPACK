/*! \file mm_dialogs_wx.cpp

    Dialogs for Molecular Mechanics and Related Modules
 
    \author Igor Kurnikov  
    \date 2010-
*/

#define MM_DIALOGS_WX_CPP

#include <mpi.h>

#include "hastl.h"
#include <cstdlib>

#include "wx/wx.h"
#include "wx/notebook.h"
#include "wx/valgen.h"
#include "wx/filename.h"

#include "hatypes.h"

#include "ctrl_wx.h"
#include "ha_wx_aux_1.h"
#include "canvas3d.h"
#include "dialogs_wx_1.h"

#include "harlemapp.h"
#include "hampi.h"
#include "hamolset.h"
#include "hamolview.h"
#include "haatgroup.h"
#include "hacompmod.h"
#include "hamolmech.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "mm_driver_amber.h"
#include "mm_driver_tinker.h"
#include "mm_driver_gromacs.h"
#include "mm_driver_arbalest.h"
#include "mm_traj_anal.h"
#include "mm_force_field.h"

#include "mm_dialogs_wx.h"
#include "ha_wx_res_wdr.h"

int MolMechDlgWX::dlg_open = FALSE;

MolMechDlgWX::MolMechDlgWX(HaMolMechMod* ptr_mm_mod_new, wxWindow *parent ):
wxFrame( parent, -1, "Molecular Mechanics Module")
{
	this->SetExtraStyle(wxWS_EX_VALIDATE_RECURSIVELY);
	ptr_mm_mod = ptr_mm_mod_new;
    dlg_open = true;

	cur_bond = NULL;
	cur_vang = NULL;
	cur_dih  = NULL;
	cur_impr_dih = NULL;

	p_mm_info_dlg = new MMInfoDlg(ptr_mm_mod, this);

	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);

	wxMenuBar* mm_menu_bar = mm_menu();
    this->SetMenuBar(mm_menu_bar);

	mol_mech_dlg( this, TRUE );

//    wxMenuBar* edit_groups_menu_bar = edit_groups_menu();
//    SetMenuBar(edit_groups_menu_bar); 

	ptr_mm_mod->p_mm_dlg = this;

	

	OnInitDialog();
}

MolMechDlgWX::~MolMechDlgWX()
{
     dlg_open = FALSE;
	 delete p_mm_info_dlg;
}


BEGIN_EVENT_TABLE(MolMechDlgWX,wxFrame)
    EVT_BUTTON(IDC_FFPAR_INIT_MOLMECH, MolMechDlgWX::OnInitMolMech)
	EVT_BUTTON(IDC_RESTRAIN_HBONDS, MolMechDlgWX::OnRestrainHBonds)
	EVT_BUTTON(IDC_MM_CHOOSE_RESTRAINED_ATOMS, MolMechDlgWX::ChooseRestrainedAtoms)
	EVT_BUTTON(IDC_MM_CHOOSE_MOVING_ATOMS,     MolMechDlgWX::ChooseMovingAtoms)
	EVT_BUTTON(IDC_MM_CHOOSE_RESTR_REF_CRD,  MolMechDlgWX::ChooseRestrRefCrd)
	EVT_BUTTON(IDC_MM_LOAD_HARM_CONSTR_FILE,  MolMechDlgWX::LoadAtomAtomRestrFile)
	EVT_BUTTON(IDC_SAVE_EXT_PROG_INP, MolMechDlgWX::OnSaveExtProgInp)
	EVT_BUTTON(IDC_MM_RUN_CALC, MolMechDlgWX::OnRunMMCalc)
	EVT_BUTTON(IDC_AMBER_LOAD_RESTART, MolMechDlgWX::OnAmberLoadRestart)
	EVT_BUTTON(IDC_MM_CALC_SP_ENE, MolMechDlgWX::OnCalcSinglePtEne)
	EVT_BUTTON(IDC_MM_LOAD_LOG_FILE, MolMechDlgWX::OnLoadLogFile)
	EVT_BUTTON(IDC_MM_CHOOSE_MDCRD_FILE, MolMechDlgWX::OnChooseMDCrdFile)
	EVT_BUTTON(IDC_MM_CHOOSE_MDVEL_FILE, MolMechDlgWX::OnChooseMDVelFile)
	EVT_BUTTON(IDC_MM_CHOOSE_MDENE_FILE, MolMechDlgWX::OnChooseMDEneFile)
	EVT_BUTTON(IDC_MM_CHOOSE_CONSTR_TRAJ_FILE, MolMechDlgWX::OnChooseConstrTrajFile)
	EVT_BUTTON(IDC_MM_STOP, MolMechDlgWX::OnMMStop)
	EVT_BUTTON(IDC_MM_PLAYBACK_TRJ, MolMechDlgWX::OnPlaybackTrj)
	EVT_BUTTON(IDC_MM_INDEX_TRAJ, MolMechDlgWX::OnIndexTrj)
	EVT_BUTTON(IDC_MM_SET_CURR_PT, MolMechDlgWX::OnSetCurrPt)
	EVT_BUTTON(IDC_MM_CHOOSE_MDANAL_SCRIPT, MolMechDlgWX::OnChooseMDAnalScript)
	EVT_BUTTON(IDC_MM_EDIT_MDANAL_SCRIPT, MolMechDlgWX::OnEditMDAnalScript)
	EVT_CHECKBOX(IDC_MM_CHK_RMSD_ANAL, MolMechDlgWX::OnChkRMSDAnal)
	EVT_CHECKBOX(IDC_MM_RUN_INT, MolMechDlgWX::OnChkMMRunInt)
	EVT_CHOICE(IDC_MM_RUN_TYPE, MolMechDlgWX::OnChoiceRunType )
	EVT_BUTTON(IDC_MM_CHOOSE_FIT_ATOMS, MolMechDlgWX::ChooseFitAtoms)
	EVT_BUTTON(IDC_MM_CHOOSE_RMSD_ATOMS, MolMechDlgWX::ChooseRMSDAtoms)
	EVT_BUTTON(IDC_MM_CHOOSE_REF_CRD_FIT,  MolMechDlgWX::ChooseRefCrdFit)
	EVT_BUTTON(IDC_MM_CHOOSE_REF_CRD_RMSD, MolMechDlgWX::ChooseRefCrdRMSD)
	EVT_BUTTON(IDC_MM_EDIT_AMBER_INP, MolMechDlgWX::OnEditAmberInp)
	EVT_BUTTON(IDC_MM_EDIT_AMBER_TOP, MolMechDlgWX::OnEditAmberTop)
	EVT_BUTTON(IDC_MM_EDIT_AMBER_RUN, MolMechDlgWX::OnEditAmberRun)
	EVT_BUTTON(IDC_MM_EDIT_AMBER_RST, MolMechDlgWX::OnEditAmberRst)
	EVT_BUTTON(IDC_AMBER_SAVE_RESTART, MolMechDlgWX::OnAmberSaveRestart)
	EVT_BUTTON(IDC_MM_SHOW_INFO, MolMechDlgWX::OnShowMMInfo)
	EVT_BUTTON(IDC_MM_EDIT_PERIODIC_BOX, MolMechDlgWX::OnEditPeriodicBox)
	EVT_COMBOBOX(IDC_MM_FF_TYPE_DEFAULT, MolMechDlgWX::OnSelChangeFFTypeDefault)
	EVT_BUTTON(IDC_MM_UPDATE_ELEM_LIST, MolMechDlgWX::OnUpdateElemList)
//    EVT_BUTTON(IDC_SELECT_FF_PARAM_FILE, MolMechDlgWX::OnSelectFFParamFile)
    EVT_RADIOBOX( IDC_RADIO_ELEMENTS, MolMechDlgWX::OnUpdateElemList)
	EVT_LISTBOX(IDC_MM_ELEM_LIST, MolMechDlgWX::OnChangeSelElem)
	EVT_BUTTON(IDC_RESPAR_SET_NEW_PAR, MolMechDlgWX::OnSetNewPar)
	EVT_BUTTON(IDC_DEL_IMPROPER_ANG, MolMechDlgWX::OnDelImprAng)
	EVT_BUTTON(IDC_SET_DNA_CS_PARS, MolMechDlgWX::OnSetDnaCSPars)
	EVT_MENU(IDC_SAVE_RESFF_FROM_MORT, MolMechDlgWX::OnSaveResffFromMort)
	EVT_MENU(IDM_SAVE_ATOM_IND_RESTR_ARB, MolMechDlgWX::OnSaveAtomIndRestrArb)
	EVT_NOTEBOOK_PAGE_CHANGING(ID_MM_DLG, MolMechDlgWX::OnChangingPage)
	EVT_NOTEBOOK_PAGE_CHANGED (ID_MM_DLG, MolMechDlgWX::OnChangePage)
	EVT_CLOSE(MolMechDlgWX::OnClose)
END_EVENT_TABLE()

void MolMechDlgWX::OnClose(wxCloseEvent& event)
{
	dlg_open = FALSE;
	event.Skip();
}

void  MolMechDlgWX::OnInitDialog()
{
	wxTextCtrl* edit_ctrl;
	wxStaticText* stext_ctrl;
	wxChoice* choice_ctrl;

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_AMBER_RUN_FILE);
	ext_prog_controls.push_back(edit_ctrl);

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_LOG_FILE);

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_INP_FILE);
	ext_prog_controls.push_back(edit_ctrl);

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_TOP_FILE);
	ext_prog_controls.push_back(edit_ctrl);

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_INIT_CRD_FILE);
	ext_prog_controls.push_back(edit_ctrl);

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_CONSTR_CRD_FILE);
	ext_prog_controls.push_back(edit_ctrl);

	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_MM_TRAJ_SCRIPT);
	edit_ctrl->SetValidator( StdStringValidator(&ptr_mm_mod->p_traj_anal_mod->traj_script) );
	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_MM_MAX_COMP_CYCLES);	      ene_min_controls.push_back(edit_ctrl);
	edit_ctrl->SetValidator( wxGenericValidator(&ptr_mm_mod->max_num_minim_steps) );
	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_MM_INIT_MIN_STEP);           ene_min_controls.push_back(edit_ctrl);
	edit_ctrl->SetValidator( wxDoubleValidator(&ptr_mm_mod->init_min_step,"%8.3f") );
	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_MM_NUM_STEP_STEEP_DESCENT);  ene_min_controls.push_back(edit_ctrl);
	edit_ctrl->SetValidator( wxGenericValidator(&ptr_mm_mod->num_steep_descent_steps) );
	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_MM_MAX_COMP_CYCLES);         ene_min_controls.push_back(edit_ctrl);
	edit_ctrl->SetValidator( wxGenericValidator(&ptr_mm_mod->max_num_minim_steps) );
	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_MM_MIN_CNVRG_CRITERIUM);     ene_min_controls.push_back(edit_ctrl);
	edit_ctrl->SetValidator( wxDoubleValidator(&ptr_mm_mod->grad_cnvrg_val,"%10.5f") );

	fit_atoms_grp_name = "ALL";
	rmsd_atoms_grp_name = "SAME AS FIT";
	
	stext_ctrl = (wxStaticText*) FindWindow( ID_TEXT_PARAM_MIN );         ene_min_controls.push_back(stext_ctrl);
	stext_ctrl = (wxStaticText*) FindWindow( ID_TEXT_MIN_TYPE );          ene_min_controls.push_back(stext_ctrl);
	stext_ctrl = (wxStaticText*) FindWindow( ID_TEXT_NUM_MIN_STEPS );     ene_min_controls.push_back(stext_ctrl);
	stext_ctrl = (wxStaticText*) FindWindow( ID_TEXT_NUM_STEEP_DEC_STEPS ); ene_min_controls.push_back(stext_ctrl);
	stext_ctrl = (wxStaticText*) FindWindow( ID_TEXT_CONV_CRT );          ene_min_controls.push_back(stext_ctrl);
	stext_ctrl = (wxStaticText*) FindWindow( ID_TEXT_KCAL_MOL_A );        ene_min_controls.push_back(stext_ctrl);
	stext_ctrl = (wxStaticText*) FindWindow( ID_TEXT_INIT_STEP );         ene_min_controls.push_back(stext_ctrl);

	stext_ctrl = (wxStaticText*) FindWindow( ID_TEXT_MD_PARAMS );         md_controls.push_back(stext_ctrl);
	stext_ctrl = (wxStaticText*) FindWindow( ID_TEXT_TEMP_CTRL_METH );    md_controls.push_back(stext_ctrl);
	stext_ctrl = (wxStaticText*) FindWindow( ID_TEXT_INIT_TEMP );         md_controls.push_back(stext_ctrl);
	stext_ctrl = (wxStaticText*) FindWindow( ID_TEXT_REF_TEMP );          md_controls.push_back(stext_ctrl);
	stext_ctrl = (wxStaticText*) FindWindow( ID_TEXT_LANG_DUMP_CONST );   md_controls.push_back(stext_ctrl);
	stext_ctrl = (wxStaticText*) FindWindow( ID_TEXT_MD_STEPS_NUM );      md_controls.push_back(stext_ctrl);
	stext_ctrl = (wxStaticText*) FindWindow( ID_TEXT_REM_COM_MOTION_FREQ );   md_controls.push_back(stext_ctrl);
	stext_ctrl = (wxStaticText*) FindWindow( ID_TEXT_START_VEL_METH );    md_controls.push_back(stext_ctrl);
	stext_ctrl = (wxStaticText*) FindWindow( ID_TEXT_MD_TIME_STEP );      md_controls.push_back(stext_ctrl);
	stext_ctrl = (wxStaticText*) FindWindow( ID_TEXT_PRESS_CTRL );        md_controls.push_back(stext_ctrl);
	stext_ctrl = (wxStaticText*) FindWindow( ID_TEXT_INIT_INFO_READ );    md_controls.push_back(stext_ctrl);


	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_MM_LENGTH_MD_RUN);           md_controls.push_back(edit_ctrl);
	edit_ctrl->SetValidator( wxGenericValidator(&ptr_mm_mod->num_md_steps) );
	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_MM_REMOVE_RB_MOTION_FREQ);   md_controls.push_back(edit_ctrl);
	edit_ctrl->SetValidator( wxGenericValidator(&ptr_mm_mod->remove_rb_motion_freq) );
	edit_ctrl= (wxTextCtrl*) FindWindow(IDC_CONSTRAIN_FRCCONST);
	edit_ctrl->SetValue("0.0");
	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_MM_INIT_TEMP);               md_controls.push_back(edit_ctrl);
	edit_ctrl->SetValidator( wxDoubleValidator(&ptr_mm_mod->init_temp,"%6.2f") );
	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_MM_REF_TEMP);                md_controls.push_back(edit_ctrl);
	edit_ctrl->SetValidator( wxDoubleValidator(&ptr_mm_mod->ref_temp,"%6.2f") );
	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_LANGEVIN_DUMP_CONST);        md_controls.push_back(edit_ctrl);
	edit_ctrl->SetValidator( wxDoubleValidator(&ptr_mm_mod->langevin_dump_const,"%6.2f") );
	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_MM_MD_TIME_STEP);            md_controls.push_back(edit_ctrl);
	edit_ctrl->SetValidator( wxDoubleValidator(&ptr_mm_mod->md_time_step,"%9.5f") );

	DDX_Text_int (this,IDC_MM_PRINT_FREQ,ptr_mm_mod->wrt_log_freq);
	DDX_Text_int (this,IDC_MM_WRT_RSTRT_FREQ,ptr_mm_mod->wrt_rstrt_freq);
    DDX_Text_int (this,IDC_MM_WRT_COORD_FREQ,ptr_mm_mod->wrt_coord_freq);
    DDX_Text_int (this,IDC_MM_WRT_VEL_FREQ,ptr_mm_mod->wrt_vel_freq);
	
	DDX_Text_int (this,IDC_MM_WRT_ENER_FREQ,ptr_mm_mod->wrt_ener_freq);
	DDX_Text_int (this,IDC_MM_WRT_CONSTR_FREQ,ptr_mm_mod->wrt_constr_freq);

	DDX_Text_double(this,IDC_MM_UPDATE_VIEW_INTERVAL,ptr_mm_mod->update_view_interval,"%9.3f");
	DDX_Text_double(this,IDC_MM_DELAY,ptr_mm_mod->p_traj_anal_mod->delay_time,"%9.3f");
	DDX_Text_int (this,IDC_MM_SKIP_INIT,ptr_mm_mod->p_traj_anal_mod->npt_begin);
	DDX_Text_int (this,IDC_MM_SKIP_BETWEEN,ptr_mm_mod->p_traj_anal_mod->npt_step);
	DDX_Text_int (this,IDC_MM_LAST_PT_IDX,ptr_mm_mod->p_traj_anal_mod->npt_end);
	DDX_Text_int (this,IDC_MM_CURR_PT,ptr_mm_mod->p_traj_anal_mod->ipt_curr);

	DDX_Choice_HaEnum(this,IDC_MM_FF_TYPE_DEFAULT,&MMForceField::ff_type_default);

	choice_ctrl_run_type = (wxChoice*) FindWindow( IDC_MM_RUN_TYPE ); 
	choice_ctrl_run_type->SetValidator( HaEnumValidator(&ptr_mm_mod->run_type) );

	choice_ctrl_per_bcond = (wxChoice*) FindWindow(IDC_MM_PERIOD_BCOND);
	choice_ctrl_per_bcond->SetValidator( HaEnumValidator(&ptr_mm_mod->period_bcond) );

	choice_ctrl_electr_method = (wxChoice*) FindWindow(IDC_MM_ELECTR_MODEL);
	choice_ctrl_electr_method->SetValidator( HaEnumValidator(&ptr_mm_mod->p_mm_model->electr_method) );

	choice_ctrl = (wxChoice*) FindWindow( IDC_MM_MIN_TYPE ); ene_min_controls.push_back(choice_ctrl);
	choice_ctrl->SetValidator( HaEnumValidator(&ptr_mm_mod->min_type) );

	choice_ctrl = (wxChoice*) FindWindow( IDC_TEMP_CONTROL_METHOD ); md_controls.push_back(choice_ctrl);
	choice_ctrl->SetValidator( HaEnumValidator(&ptr_mm_mod->temp_control_method) );

	choice_ctrl = (wxChoice*) FindWindow( IDC_MM_START_VEL_METHOD ); md_controls.push_back(choice_ctrl);
	choice_ctrl->SetValidator( HaEnumValidator(&ptr_mm_mod->start_vel_method) );

	choice_ctrl = (wxChoice*) FindWindow( IDC_MM_PRESSURE_REG_METHOD ); md_controls.push_back(choice_ctrl);
	choice_ctrl->SetValidator( HaEnumValidator(&ptr_mm_mod->pressure_reg_method) );

	choice_ctrl = (wxChoice*) FindWindow( IDC_MM_INIT_READ_COORD );   md_controls.push_back(choice_ctrl);
	choice_ctrl->SetValidator( HaEnumValidator(&ptr_mm_mod->init_read_coord) );

	choice_ctrl = (wxChoice*) FindWindow( IDC_MM_SHAKE_METHOD );     md_controls.push_back(choice_ctrl);
	choice_ctrl->SetValidator( HaEnumValidator(&ptr_mm_mod->shake_constr) );

	choice_ctrl_ext_prog = (wxChoice*) FindWindow( IDC_MM_EXT_PROG );  ext_prog_controls.push_back(choice_ctrl_ext_prog);
	choice_ctrl_ext_prog->SetValidator( HaEnumValidator(&ptr_mm_mod->ext_mm_prog) );
	choice_ctrl_ext_prog->Bind(wxEVT_CHOICE, &MolMechDlgWX::OnChangeExternalProg, this);

	btn_save_inp_files  = (wxButton*) FindWindow( IDC_SAVE_EXT_PROG_INP );
	ext_prog_controls.push_back(btn_save_inp_files);
	btn_mm_run_calc     = (wxButton*) FindWindow( IDC_MM_RUN_CALC );
	wxButton* btn = (wxButton*) FindWindow( IDC_MM_EDIT_AMBER_INP );
	ext_prog_controls.push_back(btn);
	btn = (wxButton*) FindWindow( IDC_MM_EDIT_AMBER_TOP );
	ext_prog_controls.push_back(btn);
	btn = (wxButton*) FindWindow( IDC_MM_EDIT_AMBER_RUN );
	ext_prog_controls.push_back(btn);
	
	wxCheckBox* check_box = (wxCheckBox*) FindWindow( IDC_MM_RESTART_FLAG );
	check_box->SetValidator( IntCheckBoxValidator(&ptr_mm_mod->restart_flag) );
	check_box = (wxCheckBox*) FindWindow( IDC_MM_USE_MORT );
	check_box->SetValidator( IntCheckBoxValidator(&ptr_mm_mod->p_mm_model->setup_params_from_mort_flag) );

	chk_run_internal = (wxCheckBox*) FindWindow( IDC_MM_RUN_INT );
	chk_run_internal->SetValidator( IntCheckBoxValidator(&ptr_mm_mod->run_internal_flag) );
	check_box = (wxCheckBox*) FindWindow( IDC_MM_UPDATE_MOL_VIEW );
	check_box->SetValidator( IntCheckBoxValidator(&ptr_mm_mod->update_view_flag) );
	check_box = (wxCheckBox*) FindWindow( IDC_MM_ANAL_RUN_IN_THREAD );
	check_box->SetValue(true);
	check_box = (wxCheckBox*) FindWindow( IDC_MM_REMOVE_INIT_MOTION );          md_controls.push_back(check_box);
	check_box->SetValidator( IntCheckBoxValidator(&ptr_mm_mod->remove_init_rb_motion_flag) );
	check_box = (wxCheckBox*) FindWindow( IDC_MM_WRAP_COORD );                  md_controls.push_back(check_box);
	check_box->SetValidator( IntCheckBoxValidator(&ptr_mm_mod->wrap_coord) );

	chk_rmsd_anal  = (wxCheckBox*) FindWindow( IDC_MM_CHK_RMSD_ANAL  ); 
	txt_fit_atoms  = (wxTextCtrl*) FindWindow( IDC_MM_FIT_ATOMS );      atom_superimpose_controls.push_back(txt_fit_atoms);
	txt_rmsd_atoms = (wxTextCtrl*) FindWindow( IDC_MM_RMSD_ATOMS );     atom_superimpose_controls.push_back(txt_rmsd_atoms);
	txt_rmsd_file_name = (wxTextCtrl*) FindWindow( IDC_MM_RMSD_FILE_NAME );     atom_superimpose_controls.push_back(txt_rmsd_file_name);
	txt_ref_crd_fit_fname = (wxTextCtrl*) FindWindow( IDC_MM_REF_CRD_FIT_FILE );   atom_superimpose_controls.push_back(txt_ref_crd_fit_fname);
	txt_ref_crd_rmsd_fname = (wxTextCtrl*) FindWindow( IDC_MM_REF_CRD_RMSD_FILE ); atom_superimpose_controls.push_back(txt_ref_crd_fit_fname);
	chk_rmsd_per_atom  =  (wxCheckBox*) FindWindow( IDC_MM_CHK_RMSD_ATOM );     atom_superimpose_controls.push_back(chk_rmsd_per_atom);
	chk_rmsf_per_atom  =  (wxCheckBox*) FindWindow( IDC_MM_CHK_RMSF_ATOM );     atom_superimpose_controls.push_back(chk_rmsf_per_atom);
	chk_avg_coord      =  (wxCheckBox*) FindWindow( IDC_MM_CHK_AVG_COORD );     atom_superimpose_controls.push_back(chk_avg_coord);
	txt_rmsd_per_atom_file = (wxTextCtrl*) FindWindow( IDC_MM_RMSD_ATOM_FILE ); atom_superimpose_controls.push_back(txt_rmsd_per_atom_file);
	txt_rmsf_per_atom_file = (wxTextCtrl*) FindWindow( IDC_MM_RMSF_ATOM_FILE ); atom_superimpose_controls.push_back(txt_rmsf_per_atom_file);
	txt_avg_coord_file     = (wxTextCtrl*) FindWindow( IDC_MM_AVG_COORD_FILE ); atom_superimpose_controls.push_back(txt_avg_coord_file);
	btn_choose_fit_at   = (wxButton*) FindWindow( IDC_MM_CHOOSE_FIT_ATOMS );     atom_superimpose_controls.push_back(btn_choose_fit_at);
	btn_choose_rmsd_at  = (wxButton*) FindWindow( IDC_MM_CHOOSE_RMSD_ATOMS );     atom_superimpose_controls.push_back(btn_choose_rmsd_at);
	btn_choose_ref_crd_fit  = (wxButton*) FindWindow( IDC_MM_CHOOSE_REF_CRD_FIT );   atom_superimpose_controls.push_back(btn_choose_ref_crd_fit);
	btn_choose_ref_crd_rmsd = (wxButton*) FindWindow( IDC_MM_CHOOSE_REF_CRD_RMSD );  atom_superimpose_controls.push_back(btn_choose_ref_crd_rmsd);
	btn_choose_rmsd_out_file = (wxButton*) FindWindow( IDC_MM_CHOOSE_RMSD_FILE ); atom_superimpose_controls.push_back(btn_choose_rmsd_out_file);
	btn_choose_rmsd_per_atom_file = (wxButton*) FindWindow( IDC_MM_CHOOSE_RMSD_ATOM_FILE ); atom_superimpose_controls.push_back(btn_choose_rmsd_per_atom_file);
	btn_choose_rmsf_per_atom_file = (wxButton*) FindWindow( IDC_MM_CHOOSE_RMSF_ATOM_FILE ); atom_superimpose_controls.push_back(btn_choose_rmsf_per_atom_file);
	btn_choose_avg_coord_file     = (wxButton*) FindWindow( IDC_MM_CHOOSE_AVG_COORD_FILE ); atom_superimpose_controls.push_back(btn_choose_avg_coord_file);

	combo_ref_crd_fit_type = (wxComboBox*) FindWindow( IDC_MM_REF_CRD_FIT_TYPE ); atom_superimpose_controls.push_back(combo_ref_crd_fit_type);
//	combo_ref_crd_fit_type->Append("Current Coordinates");
//	combo_ref_crd_fit_type->Append("First Trajectory Point");
//	combo_ref_crd_fit_type->Append("XYZ Coordinates File");
//	combo_ref_crd_fit_type->Append("Named Atom Group");
	combo_ref_crd_fit_type->SetSelection(0);
//	combo_ref_crd_fit_type->SetStringSelection("Current Coordinates");

	combo_ref_crd_rmsd_type = (wxComboBox*) FindWindow( IDC_MM_REF_CRD_RMSD_TYPE ); atom_superimpose_controls.push_back(combo_ref_crd_rmsd_type);
	combo_ref_crd_rmsd_type->Append("Current Coordinates");
	combo_ref_crd_rmsd_type->Append("First Trajectory Point");
	combo_ref_crd_rmsd_type->Append("XYZ Coordinates File");
//	combo_ref_crd_rmsd_type->Append("Named Atom Group");
//	combo_ref_crd_rmsd_type->SetSelection(0);
	combo_ref_crd_rmsd_type->SetStringSelection("Current Coordinates");

	combo_restr_ref_crd_type = (wxComboBox*) FindWindow( IDC_MM_RESTR_REF_CRD_TYPE ); 
	combo_restr_ref_crd_type->Append("Current Coordinates");
	combo_restr_ref_crd_type->Append("XYZ Coordinates File");
	combo_restr_ref_crd_type->SetSelection(0);

	p_si_ag = ptr_mm_mod->p_traj_anal_mod->GetRMSDAgent(FALSE);

	TransferDataToWindow();
	wxNotebook* noteb = (wxNotebook*)FindWindow(ID_MM_DLG);
	wxWindow* cur_page = noteb->GetCurrentPage();
	cur_page->Fit();
}

void MolMechDlgWX::OnChangingPage(wxNotebookEvent& event)
{
	bool bres;
	wxNotebook* noteb = (wxNotebook*) FindWindow( ID_MM_DLG );
//	PrintLog("\nMolMechDlgWX::OnChangingPage() \n");
	int np = event.GetSelection();
	int np_old = event.GetOldSelection();
//	PrintLog(" Selected Page Number = %d np_old = %d \n", np, np_old );
	wxString page_title = noteb->GetPageText(np);
	wxString page_title_old;
	if( np_old >= 0 ) page_title_old = noteb->GetPageText(np_old);
//	PrintLog(" Selected Page Title = %s \n", page_title.c_str());

	if( page_title_old == "Force Field Parameters" )
	{
		bres = TransferDataFromWindow();
		if(!bres) 
		{
			PrintLog("Invalid parameters on page \n");
			event.Veto();
		}
	}
}

void MolMechDlgWX::OnChangePage(wxNotebookEvent& event)
{
	bool bres;
	bres = TransferDataToWindow();
	if( !bres ) return; 
	wxNotebook* noteb = (wxNotebook*) FindWindow( ID_MM_DLG );
//	PrintLog("\nMolMechDlgWX::OnChangePage() \n");
	int np = event.GetSelection();
//	PrintLog(" Selected Page Number = %d \n", np);
	wxString page_title = noteb->GetPageText(np);
//	PrintLog(" Selected Page Title = %s \n", page_title.c_str());

	wxWindow* cur_page = noteb->GetCurrentPage();
	cur_page->Fit();

	if( page_title == "Edit MM Model")
	{
		OnUpdateElemList(event);
	}
}

bool MolMechDlgWX::TransferDataToWindow()
{
	int itemp;

	TransferExternalProgFileNames(true);

	wxTextCtrl* edit_ctrl = (wxTextCtrl*) FindWindow( IDC_MM_MOVING_ATOMS);
	if( edit_ctrl == NULL ) return false;
	edit_ctrl->SetValue(ptr_mm_mod->p_mm_model->moving_atoms.c_str());

	edit_ctrl = (wxTextCtrl*) FindWindow( IDC_RESTRAINED_ATOMS);
	edit_ctrl->SetValue(ptr_mm_mod->p_mm_model->restrained_atoms.c_str());

	edit_ctrl = (wxTextCtrl*) FindWindow( IDC_MM_FF_TYPE);
	edit_ctrl->SetValue( ptr_mm_mod->p_mm_model->ff_type.label()); 

	edit_ctrl= (wxTextCtrl*) FindWindow( IDC_MM_NONB_CUTOFF_DIST );
	edit_ctrl->SetValue( wxString::Format("%14.9f",ptr_mm_mod->p_mm_model->GetNBCutDist()) );

	edit_ctrl= (wxTextCtrl*) FindWindow( IDC_MM_SCALE_14_ELECTR );
	edit_ctrl->SetValue( wxString::Format("%6.2f",ptr_mm_mod->p_mm_model->GetScale14Electr()) );

	edit_ctrl= (wxTextCtrl*) FindWindow( IDC_MM_SCALE_14_VDW );
	edit_ctrl->SetValue( wxString::Format("%6.2f",ptr_mm_mod->p_mm_model->GetScale14VdW()) );

	edit_ctrl= (wxTextCtrl*) FindWindow( IDC_MM_DIEL_CONST );
	edit_ctrl->SetValue( wxString::Format("%6.2f",ptr_mm_mod->p_mm_model->GetDielConst()) );

	edit_ctrl= (wxTextCtrl*) FindWindow( IDC_MM_ION_STRENGTH );
	edit_ctrl->SetValue( wxString::Format("%6.2f",ptr_mm_mod->p_mm_model->GetIonStrength()) );

	edit_ctrl= (wxTextCtrl*) FindWindow( IDC_MM_THOLE_DUMP_CONST );
	edit_ctrl->SetValue( wxString::Format("%6.4f",ptr_mm_mod->p_mm_model->GetTholeExponCoef() ) );

	edit_ctrl= (wxTextCtrl*) FindWindow( IDC_MM_RESTR_FRC_CONST );
	edit_ctrl->SetValue( wxString::Format("%8.3f",ptr_mm_mod->p_mm_model->GetAtomRestrForceConst()) );

	edit_ctrl= (wxTextCtrl*) FindWindow( IDC_MM_NUM_HARM_CONSTR );
	edit_ctrl->SetValue( wxString::Format("%d",ptr_mm_mod->p_mm_model->GetNumHarmConstr()) );

	edit_ctrl= (wxTextCtrl*) FindWindow( IDC_MM_NPT_TRAJ );
	edit_ctrl->SetValue( wxString::Format("%d",ptr_mm_mod->p_traj_anal_mod->pt_pos.size() ) );

	TransferRunTypeDataToWindow();
	TransferExtProgDataToWindow();
	TransferAtomSuperimposeDataToWindow();
	return wxFrame::TransferDataToWindow();
}

bool MolMechDlgWX::TransferDataFromWindow()
{
	wxString str;
	int itemp;
	double dval;

	TransferExternalProgFileNames(false);

	wxTextCtrl* edit_ctrl = (wxTextCtrl*) FindWindow( IDC_MM_MOVING_ATOMS);
	str = edit_ctrl->GetValue();
	ptr_mm_mod->p_mm_model->SetMovingAtoms( str.ToStdString() );
	
	edit_ctrl = (wxTextCtrl*) FindWindow( IDC_RESTRAINED_ATOMS);
	str = edit_ctrl->GetValue();
	ptr_mm_mod->p_mm_model->SetRestrainedAtoms(str.ToStdString());

	edit_ctrl = (wxTextCtrl*) FindWindow( IDC_MM_NONB_CUTOFF_DIST );
	str = edit_ctrl->GetValue();
	if(!str.ToDouble(&dval))
	{
		edit_ctrl->SetValue( wxString::Format("%14.9f",ptr_mm_mod->p_mm_model->GetNBCutDist()) );
	}
	else
	{
		ptr_mm_mod->p_mm_model->SetNBCutDist(dval);
	}

	edit_ctrl = (wxTextCtrl*) FindWindow( IDC_MM_SCALE_14_ELECTR );
	str = edit_ctrl->GetValue();
	if(!str.ToDouble(&dval))
	{
		edit_ctrl->SetValue( wxString::Format("%6.2f",ptr_mm_mod->p_mm_model->GetScale14Electr()) );
	}
	else
	{	
		ptr_mm_mod->p_mm_model->SetScale14Electr(dval);
	}

	edit_ctrl = (wxTextCtrl*) FindWindow( IDC_MM_SCALE_14_VDW );
	str = edit_ctrl->GetValue();
	if(!str.ToDouble(&dval))
	{
		edit_ctrl->SetValue( wxString::Format("%6.2f",ptr_mm_mod->p_mm_model->GetScale14VdW()) );
	}
	else
	{
		ptr_mm_mod->p_mm_model->SetScale14VdW(dval);
	}

	edit_ctrl = (wxTextCtrl*) FindWindow( IDC_MM_DIEL_CONST );
	str = edit_ctrl->GetValue();
	if(!str.ToDouble(&dval))
	{
		edit_ctrl->SetValue( wxString::Format("%6.2f",ptr_mm_mod->p_mm_model->GetDielConst()) );
	}
	else
	{
		ptr_mm_mod->p_mm_model->SetDielConst(dval);
	}

	edit_ctrl = (wxTextCtrl*) FindWindow( IDC_MM_ION_STRENGTH );
	str = edit_ctrl->GetValue();
	if(!str.ToDouble(&dval))
	{
		edit_ctrl->SetValue( wxString::Format("%6.2f",ptr_mm_mod->p_mm_model->GetIonStrength()) );
	}
	else
	{
		ptr_mm_mod->p_mm_model->SetIonStrength(dval);
	}

	edit_ctrl = (wxTextCtrl*) FindWindow( IDC_MM_THOLE_DUMP_CONST );
	str = edit_ctrl->GetValue();
	if(!str.ToDouble(&dval))
	{
		edit_ctrl->SetValue( wxString::Format("%6.4f",ptr_mm_mod->p_mm_model->GetTholeExponCoef()) );
	}
	else
	{	
		ptr_mm_mod->p_mm_model->SetTholeExponCoef(dval);
	}

	edit_ctrl = (wxTextCtrl*) FindWindow( IDC_MM_RESTR_FRC_CONST );
	str = edit_ctrl->GetValue();
	if(!str.ToDouble(&dval))
	{
		edit_ctrl->SetValue( wxString::Format("%8.3f",ptr_mm_mod->p_mm_model->GetAtomRestrForceConst()) );
	}
	else
	{
		ptr_mm_mod->p_mm_model->SetAtomRestrForceConst(dval);
	}
	TransferAtomSuperimposeDataFromWindow();
	
	return wxFrame::TransferDataFromWindow();
}

void MolMechDlgWX::OnEditMDAnalScript(wxCommandEvent& event) 
{
	TransferDataFromWindow();	
	wxString fname_str = ptr_mm_mod->p_traj_anal_mod->traj_script.c_str();
	wxFileName fname(fname_str);

	if(!ptr_mm_mod->p_traj_anal_mod->traj_script.empty() && fname.FileExists() )
	{
		std::string cmd_line = pApp->word_editor;
		cmd_line += " ";
		cmd_line += ptr_mm_mod->p_traj_anal_mod->traj_script;

		std::system(cmd_line.c_str());
	}
	else
	{
		PrintMessage("MD Analysis Script not found");
	}
}

void MolMechDlgWX::OnChkRMSDAnal(wxCommandEvent& event)
{	
	p_si_ag = ptr_mm_mod->p_traj_anal_mod->GetRMSDAgent(FALSE);
	if( chk_rmsd_anal->GetValue() )
	{
		if( p_si_ag == NULL )
		{
			p_si_ag = ptr_mm_mod->p_traj_anal_mod->GetRMSDAgent(TRUE);
			if(p_si_ag == NULL) return;
			p_si_ag->SetActive(TRUE);
		}
		else
		{
			p_si_ag->SetActive(TRUE);
		}
	}
	else
	{
		if( p_si_ag != NULL )
		{
			p_si_ag->SetActive(FALSE);
		}
	}
	TransferAtomSuperimposeDataToWindow();
}

void MolMechDlgWX::OnChkMMRunInt(wxCommandEvent& event)
{	
	ptr_mm_mod->run_internal_flag = chk_run_internal->GetValue() ? TRUE : FALSE;
	TransferExtProgDataToWindow();
}

void MolMechDlgWX::OnChoiceRunType(wxCommandEvent& event)
{
	wxString str = choice_ctrl_run_type->GetStringSelection();
	ptr_mm_mod->SetRunType(str.ToStdString());
	TransferRunTypeDataToWindow();
	TransferDataToWindow();
	wxNotebook* noteb = (wxNotebook*)FindWindow(ID_MM_DLG);
	wxWindow* cur_page = noteb->GetCurrentPage();
	cur_page->Fit();
}

void MolMechDlgWX::ChooseMovingAtoms(wxCommandEvent& event)
{
	MolMechModel* p_mm_model = ptr_mm_mod->p_mm_model;
	if(!p_mm_model) return;

	MolSet* pmset = ptr_mm_mod->GetMolSet();
	EditGroupsDlg* edit_grp_dlg = new EditGroupsDlg( pmset, 2, NULL );
	edit_grp_dlg->ShowModal();
	AtomGroup* p_at_arr = edit_grp_dlg->GetSelGroup();
	if(p_at_arr == NULL) return;

	p_mm_model->SetMovingAtoms(p_at_arr->GetID());

	TransferDataToWindow();
	this->Raise();
}

void MolMechDlgWX::ChooseRestrainedAtoms(wxCommandEvent& event)
{
	MolMechModel* p_mm_model = ptr_mm_mod->p_mm_model;
	if(!p_mm_model) return;

	MolSet* pmset = ptr_mm_mod->GetMolSet();
	EditGroupsDlg* edit_grp_dlg = new EditGroupsDlg( pmset, 2, NULL );
	edit_grp_dlg->ShowModal();
	AtomGroup* p_at_arr = edit_grp_dlg->GetSelGroup();
	if(p_at_arr == NULL) return;

	p_mm_model->SetRestrainedAtoms(p_at_arr->GetID());

	TransferDataToWindow();
	this->Raise();
}

void MolMechDlgWX::ChooseFitAtoms(wxCommandEvent& event)
{
	if(!p_si_ag) return;

	MolSet* pmset = ptr_mm_mod->GetMolSet();
	EditGroupsDlg* edit_grp_dlg = new EditGroupsDlg( pmset, 2, NULL );
	edit_grp_dlg->ShowModal();
	AtomGroup* p_at_arr = edit_grp_dlg->GetSelGroup();
	if(p_at_arr == NULL) return;

	p_si_ag->SetAtomsFit(p_at_arr->GetID());
	fit_atoms_grp_name = p_at_arr->GetID();

	txt_fit_atoms->SetValue(fit_atoms_grp_name.c_str());

	this->Raise();
}

void MolMechDlgWX::ChooseRMSDAtoms(wxCommandEvent& event)
{
	if(!p_si_ag) return;

	MolSet* pmset = ptr_mm_mod->GetMolSet();
	EditGroupsDlg* edit_grp_dlg = new EditGroupsDlg( pmset, 2, NULL );
	edit_grp_dlg->ShowModal();
	AtomGroup* p_at_arr = edit_grp_dlg->GetSelGroup();
	if(p_at_arr == NULL) return;

	p_si_ag->SetAtomsRMSD(p_at_arr->GetID());
	rmsd_atoms_grp_name = p_at_arr->GetID();

	TransferDataToWindow();
	this->Raise();
}

void MolMechDlgWX::ChooseRefCrdFit(wxCommandEvent& event)
{
	if(!p_si_ag) return;

	if( combo_ref_crd_fit_type->GetSelection() == RMSDAgent::REFC_XYZ_CRD_FILE )
	{
		p_si_ag->ref_crd_fit_type = RMSDAgent::REFC_XYZ_CRD_FILE;
		wxString fname_inp = ::wxFileSelector("Select XYZ file to read Reference coordinates",
		                     ::wxGetCwd(),txt_ref_crd_fit_fname->GetValue(),"xyz","*.xyz");

		p_si_ag->SetRefCrdFitFromXYZFile( fname_inp.ToStdString() );
		ref_crd_fit_file_name = fname_inp.c_str();
	}
	
	TransferDataToWindow();
	this->Raise();
}

void MolMechDlgWX::ChooseRefCrdRMSD(wxCommandEvent& event)
{
	if(!p_si_ag) return;

	if( combo_ref_crd_rmsd_type->GetSelection() == RMSDAgent::REFC_XYZ_CRD_FILE )
	{
		p_si_ag->ref_crd_rmsd_type  = RMSDAgent::REFC_XYZ_CRD_FILE;
		wxString fname_inp = ::wxFileSelector("Select XYZ file to read Reference coordinates for RMSD calculations",
		                     ::wxGetCwd(),txt_ref_crd_rmsd_fname->GetValue(),"xyz","*.xyz");

		PrintLog(" MolMechDlgWX::ChooseRefCrdRMSD() fname_inp = %s \n",fname_inp.ToStdString().c_str());
		p_si_ag->SetRefCrdRMSDFromXYZFile( fname_inp.ToStdString() );
		ref_crd_rmsd_file_name = fname_inp.c_str();
	}
	
	TransferDataToWindow();
	this->Raise();
}

void MolMechDlgWX::ChooseRestrRefCrd(wxCommandEvent& event)
{
	PrintLog(" MolMechDlgWX::ChooseRestrRefCrd() \n");
	MolMechModel* p_mm_model = ptr_mm_mod->p_mm_model;
	if(!p_mm_model) return;

	if( combo_restr_ref_crd_type->GetSelection() == MolMechModel::RESTR_REFC_XYZ_CRD_FILE )
	{
		wxString fname_inp = ::wxFileSelector("Select XYZ file to read Restraints Reference coordinates",
		                     ::wxGetCwd(),"","xyz","*.xyz");

		p_mm_model->SetRestrRefCrdFromXYZFile(fname_inp.c_str());
	}
	TransferDataToWindow();
	this->Raise();
}

void MolMechDlgWX::LoadAtomAtomRestrFile(wxCommandEvent& event)
{
	PrintLog(" MolMechDlgWX::LoadAtomAtomRestraintFile() \n");
	MolMechModel* p_mm_model = ptr_mm_mod->p_mm_model;
	if(!p_mm_model) return;

	wxString fname_inp = ::wxFileSelector("Select Atom-Atom Distance Restraint File",
		                 ::wxGetCwd(),txt_ref_crd_fit_fname->GetValue(),"Atom-Atom Distance Restraint File","*.dat");
	
	p_mm_model->SetDistConstrFromFile(fname_inp.c_str());
	TransferDataToWindow();
	this->Raise();
}

void MolMechDlgWX::OnChangePeriodicity()
{
	if( choice_ctrl_per_bcond )
	{
		HaEnumValidator* pval = dynamic_cast<HaEnumValidator*>(choice_ctrl_per_bcond->GetValidator());
		if( pval ) pval->UpdateItems();
	}
	if( choice_ctrl_electr_method )
	{
		HaEnumValidator* pval = dynamic_cast<HaEnumValidator*>(choice_ctrl_electr_method->GetValidator());
		if( pval ) pval->UpdateItems();
	}
	TransferDataToWindow();
}

void HaMolMechMod::OnChangePeriodicity()
{
	if( p_mm_dlg ) p_mm_dlg->OnChangePeriodicity();	
}

void MolMechDlgWX::TransferRunTypeDataToWindow()
{
	int ic;
	int nc_md  =  md_controls.size();
	int nc_min =  ene_min_controls.size();

	if( ptr_mm_mod->run_type == ptr_mm_mod->run_type.MD_RUN )
	{
		for(ic = 0; ic < nc_min; ic++)
		{
	//		ene_min_controls[ic]->Disable();
			ene_min_controls[ic]->Show(false);
		}
		for(ic = 0; ic < nc_md; ic++)
		{
	//		md_controls[ic]->Enable();
			md_controls[ic]->Show(true);
		}
	}
	else if( ptr_mm_mod->run_type == ptr_mm_mod->run_type.MIN_RUN )
	{
		for(ic = 0; ic < nc_md; ic++)
		{
	//		md_controls[ic]->Disable();
			md_controls[ic]->Show(false);
		}
		for(ic = 0; ic < nc_min; ic++)
		{
	//		ene_min_controls[ic]->Enable();
			ene_min_controls[ic]->Show(true);
			HaEnumValidator* pval = dynamic_cast<HaEnumValidator*>(ene_min_controls[ic]->GetValidator());
			if (pval) pval->UpdateItems();
		}
	}
	else if( ptr_mm_mod->run_type == ptr_mm_mod->run_type.ENER_RUN )
	{
		for(ic = 0; ic < nc_md; ic++)
		{
//			md_controls[ic]->Disable();
			md_controls[ic]->Show(false);
		}
		for(ic = 0; ic < nc_min; ic++)
		{
//			ene_min_controls[ic]->Disable();
			ene_min_controls[ic]->Show(false);
		}
	}
}


void MolMechDlgWX::TransferExtProgDataToWindow()
{
	int ic;
	int nc = ext_prog_controls.size();
	if( ptr_mm_mod->run_internal_flag ) // Set to run using internal MD code
	{
		for(ic = 0; ic < nc; ic++)
		{
			ext_prog_controls[ic]->Disable();
		}
		btn_mm_run_calc->Enable();
	}
	else  // Set to run using external MD program
	{
		for(ic = 0; ic < nc; ic++)
		{
			ext_prog_controls[ic]->Enable();
		}
		if( /* ptr_mm_mod->ext_mm_prog == ptr_mm_mod->ext_mm_prog.PMEMD_9 ||
			ptr_mm_mod->ext_mm_prog == ptr_mm_mod->ext_mm_prog.SANDER_9 || 
			ptr_mm_mod->ext_mm_prog == ptr_mm_mod->ext_mm_prog.PMEMD_10 || 
			ptr_mm_mod->ext_mm_prog == ptr_mm_mod->ext_mm_prog.PMEMD_12 || */ 
			ptr_mm_mod->ext_mm_prog == ptr_mm_mod->ext_mm_prog.PMEMD_18)
		{
			if( ptr_mm_mod->p_amber_driver->to_save_input_files ) 
			{
				btn_mm_run_calc->Disable();
			}
			else
			{
				btn_mm_run_calc->Enable();
			}
		}
		else if(  ptr_mm_mod->ext_mm_prog == ptr_mm_mod->ext_mm_prog.TINKER_51 )
		{
			if( ptr_mm_mod->p_tinker_driver->to_save_input_files ) 
			{
				btn_mm_run_calc->Disable();
			}
			else
			{
				btn_mm_run_calc->Enable();
			}
		}
		else if(  ptr_mm_mod->ext_mm_prog == ptr_mm_mod->ext_mm_prog.GROMACS_51 )
		{
			if( ptr_mm_mod->p_gromacs_driver->to_save_input_files ) 
			{
				btn_mm_run_calc->Disable();
			}
			else
			{
				btn_mm_run_calc->Enable();
			}
		}
		else if (ptr_mm_mod->ext_mm_prog == ptr_mm_mod->ext_mm_prog.ARBALEST_25)
		{
			if (ptr_mm_mod->p_arbalest_driver->to_save_input_files)
			{
				btn_mm_run_calc->Disable();
			}
			else
			{
				btn_mm_run_calc->Enable();
			}
		}
	}
}


void MolMechDlgWX::TransferAtomSuperimposeDataToWindow()
{
	if(p_si_ag && p_si_ag->IsActive()) chk_rmsd_anal->SetValue(true);
	else chk_rmsd_anal->SetValue(false);

	if( p_si_ag )
	{
		txt_fit_atoms->SetValue( fit_atoms_grp_name );
		txt_rmsd_atoms->SetValue( rmsd_atoms_grp_name );
		txt_rmsd_file_name->SetValue( p_si_ag->fname_rmsd_out.c_str() );
		combo_ref_crd_fit_type->SetSelection( p_si_ag->ref_crd_fit_type );
		combo_ref_crd_rmsd_type->SetSelection( p_si_ag->ref_crd_rmsd_type );
		txt_rmsd_per_atom_file->SetValue( p_si_ag->fname_rmsd_atom_out.c_str() );
		txt_rmsf_per_atom_file->SetValue( p_si_ag->fname_rmsf_atom_out.c_str() );
		txt_avg_coord_file->SetValue( p_si_ag->avg_crd_file_name );
		wxFileName ref_crd_fit_wx_fname( ref_crd_fit_file_name.c_str() );
		ref_crd_fit_wx_fname.MakeRelativeTo("");
		txt_ref_crd_fit_fname->SetValue( ref_crd_fit_wx_fname.GetFullName() );
		wxFileName ref_crd_rmsd_wx_fname( ref_crd_rmsd_file_name.c_str() );
		ref_crd_rmsd_wx_fname.MakeRelativeTo("");
		txt_ref_crd_rmsd_fname->SetValue( ref_crd_rmsd_wx_fname.GetFullName() );
		if( p_si_ag->calc_rmsd_per_atom_flag ) chk_rmsd_per_atom->SetValue(true);
		else chk_rmsd_per_atom->SetValue(false);
		if( p_si_ag->calc_rmsf_per_atom_flag ) chk_rmsf_per_atom->SetValue(true);
		else chk_rmsf_per_atom->SetValue(false);
		if( p_si_ag->calc_avg_crd_flag ) chk_avg_coord->SetValue(true);
		else chk_avg_coord->SetValue(false);
	}
	else
	{
		txt_fit_atoms->SetValue("");
		txt_rmsd_atoms->SetValue("");
		txt_rmsd_file_name->SetValue("");
		combo_ref_crd_fit_type->SetSelection(0);
		combo_ref_crd_rmsd_type->SetSelection(0);
		txt_rmsd_per_atom_file->SetValue("");
		txt_rmsf_per_atom_file->SetValue("");
		txt_ref_crd_fit_fname->SetValue("");
		chk_rmsd_per_atom->SetValue(false);
		chk_avg_coord->SetValue(false);
		txt_avg_coord_file->SetValue("");
	}

	int ic;
	int nc = atom_superimpose_controls.size();
	if( chk_rmsd_anal->GetValue() )
	{
		for(ic = 0; ic < nc; ic++)
		{
			atom_superimpose_controls[ic]->Enable();
		}
	}
	else
	{
		for(ic = 0; ic < nc; ic++)
		{
			atom_superimpose_controls[ic]->Disable();
		}
	}
}

void MolMechDlgWX::TransferAtomSuperimposeDataFromWindow()
{
	if( chk_rmsd_anal->GetValue() )
	{
		if( p_si_ag == NULL )
		{
			p_si_ag = ptr_mm_mod->p_traj_anal_mod->GetRMSDAgent(TRUE);
			if(p_si_ag == NULL) return;
			p_si_ag->SetActive(TRUE);
			TransferAtomSuperimposeDataToWindow();
		}
		else
		{
			p_si_ag->SetActive(TRUE);
		}
	}
	else
	{
		if( p_si_ag != NULL ) p_si_ag->SetActive(FALSE);
		return;
	}

	fit_atoms_grp_name = txt_fit_atoms->GetValue(); 
	rmsd_atoms_grp_name = txt_rmsd_atoms->GetValue(); 
	p_si_ag->fname_rmsd_out      = txt_rmsd_file_name->GetValue();
	p_si_ag->ref_crd_fit_type    = (RMSDAgent::RefCrdType) combo_ref_crd_fit_type->GetSelection();
	p_si_ag->ref_crd_rmsd_type   = (RMSDAgent::RefCrdType) combo_ref_crd_rmsd_type->GetSelection();
	ref_crd_fit_file_name     = txt_ref_crd_fit_fname->GetValue();
	ref_crd_rmsd_file_name    = txt_ref_crd_rmsd_fname->GetValue();
	
	p_si_ag->fname_rmsd_atom_out   = txt_rmsd_per_atom_file->GetValue();
	p_si_ag->fname_rmsf_atom_out   = txt_rmsf_per_atom_file->GetValue();
	p_si_ag->calc_rmsd_per_atom_flag = ( chk_rmsd_per_atom->GetValue() ) ? TRUE : FALSE;
	p_si_ag->calc_rmsf_per_atom_flag = ( chk_rmsf_per_atom->GetValue() ) ? TRUE : FALSE;
	p_si_ag->avg_crd_file_name       = txt_avg_coord_file->GetValue();
	p_si_ag->calc_avg_crd_flag    = ( chk_avg_coord->GetValue() ) ? TRUE : FALSE;
}


void MolMechDlgWX::OnChooseMDAnalScript(wxCommandEvent& event) 
{
	wxString script_name = ::wxFileSelector("Choose PYTHON MD Analysis script","","",
		"py","*.py");

	TransferDataFromWindow();

	if(!script_name.empty() )
	{
		wxFileName scr_fname(script_name);
		wxString cur_dir = ::wxGetCwd();
		PrintLog("Current Directory: %s \n",cur_dir.ToStdString().c_str());
        scr_fname.MakeRelativeTo(cur_dir);
		wxString mod_script_name = scr_fname.GetFullPath();
		ptr_mm_mod->p_traj_anal_mod->traj_script = mod_script_name.c_str();
		TransferDataToWindow();
	}
	else
	{
		PrintLog("Invalid Script Name");
	}
}

void MolMechDlgWX::TransferExternalProgFileNames(bool to_window)
{
	// SANDER_12 = 5, PMEMD_12 = 6, PMEMD_18 = 7, TINKER_51 = 8, GROMACS_51 = 9

	wxTextCtrl* edit_ctrl;
	wxChoice* choice_ctrl;

	if (ptr_mm_mod->ext_mm_prog == MMExternalProg::SANDER_12 || ptr_mm_mod->ext_mm_prog ==  MMExternalProg::PMEMD_12 || ptr_mm_mod->ext_mm_prog ==  MMExternalProg::PMEMD_18)
	{
		MMDriverAmber* p_amber_driver = ptr_mm_mod->p_amber_driver;

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_AMBER_RUN_FILE);
		if (to_window) edit_ctrl->SetValue(p_amber_driver->amber_run_file);
		else p_amber_driver->amber_run_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_LOG_FILE);
		if (to_window) edit_ctrl->SetValue(p_amber_driver->amber_out_file);
		else p_amber_driver->amber_out_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_INP_FILE);
		if (to_window) edit_ctrl->SetValue(p_amber_driver->amber_inp_file);
		else p_amber_driver->amber_inp_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_TOP_FILE);
		if (to_window) edit_ctrl->SetValue(p_amber_driver->amber_top_file);
		else p_amber_driver->amber_top_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_INIT_CRD_FILE);
		if (to_window) edit_ctrl->SetValue(p_amber_driver->amber_crd_file);
		else p_amber_driver->amber_crd_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_CONSTR_CRD_FILE);
		if (to_window) edit_ctrl->SetValue(p_amber_driver->amber_constr_crd_file);
		else p_amber_driver->amber_constr_crd_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_RESTART_FILE);
		if (to_window) edit_ctrl->SetValue(p_amber_driver->amber_rst_file);
		else p_amber_driver->amber_rst_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_CRD_TRAJ_FILE);
		if (to_window) edit_ctrl->SetValue(p_amber_driver->amber_trj_coord_file);
		else p_amber_driver->amber_trj_coord_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_CRD_TRAJ_FILE_2);
		if (to_window) edit_ctrl->SetValue(p_amber_driver->amber_trj_coord_file);
		else p_amber_driver->amber_trj_coord_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_VEL_TRAJ_FILE);
		if (to_window) edit_ctrl->SetValue(p_amber_driver->amber_trj_vel_file);
		else p_amber_driver->amber_trj_vel_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_VEL_TRAJ_FILE_2);
		if (to_window) edit_ctrl->SetValue(p_amber_driver->amber_trj_vel_file);
		else p_amber_driver->amber_trj_vel_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_ENE_TRAJ_FILE);
		if (to_window) edit_ctrl->SetValue(p_amber_driver->amber_trj_ene_file);
		else p_amber_driver->amber_trj_ene_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_ENE_TRAJ_FILE_2);
		if (to_window) edit_ctrl->SetValue(p_amber_driver->amber_trj_ene_file);
		else p_amber_driver->amber_trj_ene_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_CONSTR_TRAJ_FILE);
		if (to_window) edit_ctrl->SetValue(ptr_mm_mod->constr_trj_fname);
		else ptr_mm_mod->constr_trj_fname = edit_ctrl->GetValue().ToStdString();
	}

	if (ptr_mm_mod->ext_mm_prog == MMExternalProg::GROMACS_51 )
	{
		MMDriverGromacs* p_gromacs_driver = ptr_mm_mod->p_gromacs_driver;

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_AMBER_RUN_FILE);
		if (to_window) edit_ctrl->SetValue(p_gromacs_driver->run_fname); 
		else p_gromacs_driver->run_fname = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_LOG_FILE);
		if (to_window) edit_ctrl->SetValue("");
		//else p_gromacs_driver->amber_out_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_INP_FILE);
		if (to_window) edit_ctrl->SetValue(p_gromacs_driver->inp_fname);
		else p_gromacs_driver->inp_fname = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_TOP_FILE);
		if (to_window) edit_ctrl->SetValue(p_gromacs_driver->top_fname);
		else p_gromacs_driver->top_fname = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_INIT_CRD_FILE);
		if (to_window) edit_ctrl->SetValue(p_gromacs_driver->init_crd_fname);
		else p_gromacs_driver->init_crd_fname = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_CONSTR_CRD_FILE);
		if (to_window) edit_ctrl->SetValue(p_gromacs_driver->restr_crd_fname);
		else p_gromacs_driver->restr_crd_fname = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_RESTART_FILE);
		if (to_window) edit_ctrl->SetValue("");
		//else p_amber_driver->amber_rst_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_CRD_TRAJ_FILE);
		if (to_window) edit_ctrl->SetValue(p_gromacs_driver->trj_fname);
		else p_gromacs_driver->trj_fname = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_CRD_TRAJ_FILE_2);
		if (to_window) edit_ctrl->SetValue(p_gromacs_driver->trj_fname);
		else p_gromacs_driver->trj_fname = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_VEL_TRAJ_FILE);
		if (to_window) edit_ctrl->SetValue("");
		//else p_amber_driver->amber_trj_vel_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_VEL_TRAJ_FILE_2);
		if (to_window) edit_ctrl->SetValue("");
		//else p_amber_driver->amber_trj_vel_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_ENE_TRAJ_FILE);
		if (to_window) edit_ctrl->SetValue(p_gromacs_driver->ene_fname);
		else p_gromacs_driver->ene_fname = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_ENE_TRAJ_FILE_2);
		if (to_window) edit_ctrl->SetValue(p_gromacs_driver->ene_fname);
		else p_gromacs_driver->ene_fname = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_CONSTR_TRAJ_FILE);
		if (to_window) edit_ctrl->SetValue("");
		//else ptr_mm_mod->constr_trj_fname = edit_ctrl->GetValue().ToStdString();
	}

	if (ptr_mm_mod->ext_mm_prog == MMExternalProg::ARBALEST_25)
	{
		MMDriverArbalest* p_arbalest_driver = ptr_mm_mod->p_arbalest_driver;

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_AMBER_RUN_FILE);
		if (to_window) edit_ctrl->SetValue(p_arbalest_driver->run_fname);
		else p_arbalest_driver->run_fname = edit_ctrl->GetValue().ToStdString();

		//edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_LOG_FILE);
		//if (to_window) edit_ctrl->SetValue("");
		//else p_gromacs_driver->amber_out_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_INP_FILE);
		if (to_window) edit_ctrl->SetValue(p_arbalest_driver->config_fname);
		else p_arbalest_driver->config_fname = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_TOP_FILE);
		if (to_window) edit_ctrl->SetValue("");
		//else p_gromacs_driver->top_fname = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_INIT_CRD_FILE);
		if (to_window) edit_ctrl->SetValue("");
		//else p_gromacs_driver->init_crd_fname = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_CONSTR_CRD_FILE);
		if (to_window) edit_ctrl->SetValue("");
		//else p_gromacs_driver->restr_crd_fname = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_RESTART_FILE);
		if (to_window) edit_ctrl->SetValue("");
		//else p_amber_driver->amber_rst_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_CRD_TRAJ_FILE);
		if (to_window) edit_ctrl->SetValue("");
		//else p_gromacs_driver->trj_fname = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_CRD_TRAJ_FILE_2);
		if (to_window) edit_ctrl->SetValue("");
		//else p_gromacs_driver->trj_fname = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_VEL_TRAJ_FILE);
		if (to_window) edit_ctrl->SetValue("");
		//else p_amber_driver->amber_trj_vel_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_VEL_TRAJ_FILE_2);
		if (to_window) edit_ctrl->SetValue("");
		//else p_amber_driver->amber_trj_vel_file = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_ENE_TRAJ_FILE);
		if (to_window) edit_ctrl->SetValue("");
		//else p_gromacs_driver->ene_fname = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_ENE_TRAJ_FILE_2);
		if (to_window) edit_ctrl->SetValue("");
		//else p_gromacs_driver->ene_fname = edit_ctrl->GetValue().ToStdString();

		edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_CONSTR_TRAJ_FILE);
		if (to_window) edit_ctrl->SetValue("");
		//else ptr_mm_mod->constr_trj_fname = edit_ctrl->GetValue().ToStdString();
	}

};

void MolMechDlgWX::OnChangeExternalProg(wxCommandEvent& event)
{
	int selection = event.GetSelection();
	std::string choiceText = event.GetString().ToStdString();
	PrintLog("%s\n", choiceText.c_str());

	ptr_mm_mod->ext_mm_prog.SetWithLabel(choiceText.c_str());

	TransferExtProgDataToWindow();
	TransferExternalProgFileNames(true);
}

void MolMechDlgWX::OnSaveExtProgInp(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	if( /* ptr_mm_mod->ext_mm_prog == ptr_mm_mod->ext_mm_prog.PMEMD_9 ||
		ptr_mm_mod->ext_mm_prog == ptr_mm_mod->ext_mm_prog.SANDER_9 || 
		ptr_mm_mod->ext_mm_prog == ptr_mm_mod->ext_mm_prog.PMEMD_10 || 
		ptr_mm_mod->ext_mm_prog == ptr_mm_mod->ext_mm_prog.SANDER_10 || 
		ptr_mm_mod->ext_mm_prog == ptr_mm_mod->ext_mm_prog.PMEMD_12 || 
		ptr_mm_mod->ext_mm_prog == ptr_mm_mod->ext_mm_prog.SANDER_12 || */ 
		ptr_mm_mod->ext_mm_prog == ptr_mm_mod->ext_mm_prog.PMEMD_18)
	{
		ptr_mm_mod->p_amber_driver->SaveAllInpFiles();
	}
	else if ( ptr_mm_mod->ext_mm_prog == ptr_mm_mod->ext_mm_prog.TINKER_51 )
	{
		ptr_mm_mod->p_tinker_driver->SaveAllInpFiles();
	}
	else if ( ptr_mm_mod->ext_mm_prog == ptr_mm_mod->ext_mm_prog.GROMACS_51 )
	{
		ptr_mm_mod->p_gromacs_driver->SaveAllInpFiles();
	}
	else if (ptr_mm_mod->ext_mm_prog == ptr_mm_mod->ext_mm_prog.ARBALEST_25)
	{
		ptr_mm_mod->p_arbalest_driver->SaveAllInpFiles();
	}

	TransferExtProgDataToWindow();
}

void MolMechDlgWX::OnRunMMCalc(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	harlem::RunOptions opt;
	opt.SetRunSync(false);
	ptr_mm_mod->Run(&opt);
}

void MolMechDlgWX::OnIndexTrj(wxCommandEvent& event)
{
	wxBusyCursor wait;
	TransferDataFromWindow();
	ptr_mm_mod->p_traj_anal_mod->BuildTrajIndex();
	TransferDataToWindow();
}

void MolMechDlgWX::OnSetCurrPt(wxCommandEvent& event)
{
	TransferDataFromWindow();
	ptr_mm_mod->p_traj_anal_mod->LoadCurrPt();
	ptr_mm_mod->GetMolSet()->RefreshAllViews(RFRefresh | RFApply );
	TransferDataToWindow();
}

void MolMechDlgWX::OnPlaybackTrj(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	wxCheckBox* check_box = (wxCheckBox*) FindWindow( IDC_MM_ANAL_RUN_IN_THREAD );
	int sync = check_box->IsChecked() ? FALSE : TRUE; 
	ptr_mm_mod->p_traj_anal_mod->AnalyzeTrajectory(sync);		
}

void MolMechDlgWX::OnAmberLoadRestart(wxCommandEvent& event) 
{
	wxString rst_file_name = ::wxFileSelector("Choose Restart File",  
		::wxGetCwd(), ptr_mm_mod->p_amber_driver->amber_rst_file.c_str() );
	if(!rst_file_name.empty())
	{
		ptr_mm_mod->p_amber_driver->LoadAmberRestartFile( rst_file_name.ToStdString() );	
		MolSet* pmset= ptr_mm_mod->GetMolSet();
		pmset->AnnounceGeomChange();
	}
}

void MolMechDlgWX::OnAmberSaveRestart(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	ptr_mm_mod->p_amber_driver->SaveAmberRstFile(ptr_mm_mod->p_amber_driver->amber_rst_file.c_str());	
}

void MolMechDlgWX::OnShowMMInfo(wxCommandEvent& event)
{
	p_mm_info_dlg->Show(TRUE);
	p_mm_info_dlg->TransferDataToWindow();
}

void MolMechDlgWX::OnEditPeriodicBox(wxCommandEvent& event)
{
	MolSet* pmset = ptr_mm_mod->GetMolSet();
	
	//AtomParamsDlgWX::ResetEditFlags();
	
	if(AtomParamsDlgWX::dlg_open) return;

	AtomParamsDlgWX* ptr_atom_params_dlg = new AtomParamsDlgWX( pmset, NULL );

	ptr_atom_params_dlg->Show(TRUE);	
}


void MolMechDlgWX::OnCalcSinglePtEne(wxCommandEvent& event)
{
	if(ptr_mm_mod == NULL )
		return;

	ptr_mm_mod->CalcEnergy();
	ptr_mm_mod->PrintLogEne();
}

void MolMechDlgWX::OnLoadLogFile(wxCommandEvent& event) 
{
	if(ptr_mm_mod == NULL ) return;

	std::string cmd_line = pApp->word_editor;
	cmd_line += " ";
	cmd_line += ptr_mm_mod->p_amber_driver->amber_out_file;

	std::system(cmd_line.c_str());
}

void MolMechDlgWX::OnChooseMDCrdFile(wxCommandEvent& event) 
{
	if(ptr_mm_mod == NULL )
		return;

	TransferDataFromWindow();

	wxString mdcrd_file_name = ::wxFileSelector("Choose MD Trajectory Coordinates File",
		::wxGetCwd(),ptr_mm_mod->p_amber_driver->amber_trj_coord_file.c_str(),
		"mdcrd","*.mdcrd");

	if(!mdcrd_file_name.empty() )
	{
		wxFileName scr_fname(mdcrd_file_name);
		wxString cur_dir = ::wxGetCwd();
		PrintLog("Current Directory: %s \n",cur_dir.ToStdString().c_str());
        scr_fname.MakeRelativeTo(cur_dir);
		wxString mod_mdcrd_file_name = scr_fname.GetFullPath();
		ptr_mm_mod->p_amber_driver->amber_trj_coord_file = mod_mdcrd_file_name.c_str();
		TransferDataToWindow();
	}
	else
	{
		PrintLog("Invalid MD Trajectory Coordinates File Name");
	}
	this->Raise();
}

void MolMechDlgWX::OnChooseMDVelFile(wxCommandEvent& event) 
{
	if(ptr_mm_mod == NULL )
		return;

	TransferDataFromWindow();

	wxString mdvel_file_name = ::wxFileSelector("Choose MD Trajectory Velocities File",
		::wxGetCwd(),ptr_mm_mod->p_amber_driver->amber_trj_vel_file.c_str(),
		"mdvel","*.mdvel");

	if(!mdvel_file_name.empty() )
	{
		wxFileName scr_fname(mdvel_file_name);
		wxString cur_dir = ::wxGetCwd();
		PrintLog("Current Directory: %s \n",cur_dir.ToStdString().c_str());
        scr_fname.MakeRelativeTo(cur_dir);
		wxString mod_mdvel_file_name = scr_fname.GetFullPath();
		ptr_mm_mod->p_amber_driver->amber_trj_vel_file = mod_mdvel_file_name.c_str();
		TransferDataToWindow();
	}
	else
	{
		PrintLog("Invalid MD Trajectory Velocities File Name");
	}
	this->Raise();
}

void MolMechDlgWX::OnChooseMDEneFile(wxCommandEvent& event) 
{
	TransferDataFromWindow();

	wxString mdene_file_name = ::wxFileSelector("Choose MD Trajectory Energies File",
		::wxGetCwd(),ptr_mm_mod->p_amber_driver->amber_trj_ene_file.c_str(),
		"mden","*.mden");

	if(!mdene_file_name.empty() )
	{
		wxFileName scr_fname(mdene_file_name);
		wxString cur_dir = ::wxGetCwd();
		PrintLog("Current Directory: %s \n",cur_dir.ToStdString().c_str());
        scr_fname.MakeRelativeTo(cur_dir);
		wxString mod_mdene_file_name = scr_fname.GetFullPath();
		ptr_mm_mod->p_amber_driver->amber_trj_ene_file = mod_mdene_file_name.c_str();
		TransferDataToWindow();
	}
	else
	{
		PrintLog("Invalid MD Trajectory Energy File Name");
	}
	this->Raise();
}

void MolMechDlgWX::OnChooseConstrTrajFile(wxCommandEvent& event)
{
	TransferDataFromWindow();

	wxString constr_traj_file_name = ::wxFileSelector("Choose Constraint Trajectory File",
		::wxGetCwd(),ptr_mm_mod->constr_trj_fname.c_str(),
		"dat","*.dat");

	if(!constr_traj_file_name.empty() )
	{
		wxFileName scr_fname(constr_traj_file_name);
		wxString cur_dir = ::wxGetCwd();
		PrintLog("Current Directory: %s \n",cur_dir.ToStdString().c_str());
        scr_fname.MakeRelativeTo(cur_dir);
		wxString mod_constr_traj_file_name = scr_fname.GetFullPath();
		ptr_mm_mod->constr_trj_fname = mod_constr_traj_file_name.c_str();
		TransferDataToWindow();
	}
	else
	{
		PrintLog("Invalid ConstraintTrajectory File Name");
	}
	this->Raise();	

}

void MolMechDlgWX::OnMMStop(wxCommandEvent& event) 
{
	ptr_mm_mod->StopCalc();
}

void MolMechDlgWX::OnEditAmberInp(wxCommandEvent& event) 
{
	TransferDataFromWindow();	
	wxString fname_str = ptr_mm_mod->p_amber_driver->amber_inp_file.c_str();
	wxFileName fname(fname_str);

	if(!ptr_mm_mod->p_amber_driver->amber_inp_file.empty() && fname.FileExists() )
	{
		std::string cmd_line = pApp->word_editor;
		cmd_line += " ";
		cmd_line += ptr_mm_mod->p_amber_driver->amber_inp_file;
		std::system(cmd_line.c_str());
	}
	else
	{
		PrintMessage("AMBER/SANDER input file not found");
	}
}

void MolMechDlgWX::OnEditAmberTop(wxCommandEvent& event) 
{
	TransferDataFromWindow();	
	wxString fname_str = ptr_mm_mod->p_amber_driver->amber_top_file.c_str();
	wxFileName fname(fname_str);

	if(!ptr_mm_mod->p_amber_driver->amber_top_file.empty() && fname.FileExists() )
	{
		std::string cmd_line = pApp->word_editor;
		cmd_line += " ";
		cmd_line += ptr_mm_mod->p_amber_driver->amber_top_file;
		std::system(cmd_line.c_str());
	}
	else
	{
		PrintMessage("AMBER top file (molecular topology/MM model) not found");
	}	
}

void MolMechDlgWX::OnEditAmberRun(wxCommandEvent& event) 
{
	TransferDataFromWindow();	
	wxString fname_str = ptr_mm_mod->p_amber_driver->amber_run_file.c_str();
	wxFileName fname(fname_str);

	if(!ptr_mm_mod->p_amber_driver->amber_run_file.empty() && fname.FileExists() )
	{
		std::string cmd_line = pApp->word_editor;
		cmd_line += " ";
		cmd_line += ptr_mm_mod->p_amber_driver->amber_run_file;

		std::system(cmd_line.c_str());
	}
	else
	{
		PrintMessage("AMBER/SANDER UNIX RUN script file not found");
	}		
}

void MolMechDlgWX::OnEditAmberRst(wxCommandEvent& event) 
{
	TransferDataFromWindow();	
	wxString fname_str = ptr_mm_mod->p_amber_driver->amber_rst_file.c_str();
	wxFileName fname(fname_str);

	if(!ptr_mm_mod->p_amber_driver->amber_rst_file.empty() && fname.FileExists() )
	{
		std::string cmd_line = pApp->word_editor;
		cmd_line += " ";
		cmd_line += ptr_mm_mod->p_amber_driver->amber_rst_file;

		std::system(cmd_line.c_str());
	}
	else
	{
		PrintMessage("AMBER/SANDER restart file not found");
	}		
}

//! Find a mode to print ref of atom 2 with no repeating info the same as in atom 1  
static int adjust_ref_mode(HaAtom* aptr1, HaAtom* aptr2)
{
	if( aptr1 == NULL || aptr2 == NULL) return HaAtom::ATOMREF_FULL;

	if( aptr1->GetHostRes() == aptr2->GetHostRes()) return HaAtom::ATOMREF_NO_RES;
	if( aptr1->GetHostMol() == aptr2->GetHostMol()) return HaAtom::ATOMREF_NO_MOL;
	return HaAtom::ATOMREF_FULL;
}

void MolMechDlgWX::OnSelectFFParamFile(wxCommandEvent& event)
{
//	wxFileName old_param_fname(MMForceField::ff_type_default.c_str());
//	wxString dir_name = old_param_fname.GetPath();
//	wxString path = ::wxFileSelector("Select Force Field Parameters File", dir_name, 
//		"", ".dat", "*.dat");
//	if ( !path ) return;
//
//	MMForceField* p_ff = MMForceField::GetMMForceField( MMForceField::ff_name_default,TRUE);
//	p_ff->LoadAmberParamFile( path.c_str() );
//	TransferDataToWindow();
}

void MolMechDlgWX::OnSelChangeFFTypeDefault( wxCommandEvent& event )
{
	wxComboBox* p_combo_box = (wxComboBox*) FindWindow(IDC_MM_FF_TYPE_DEFAULT);
	wxString label = p_combo_box->GetValue();
	MMForceField::ff_type_default.SetWithLabel(label.c_str());
}


void MolMechDlgWX::OnUpdateElemList(wxCommandEvent& event)
{
	wxRadioBox* elem_type = (wxRadioBox*) FindWindow(IDC_RADIO_ELEMENTS);
    wxString sel_type_str = elem_type->GetStringSelection();

//	int isel = elem_type->GetSelection();

      int i,nel;
	  char buf[256];

	  wxListBox* elem_list= (wxListBox*) FindWindow(IDC_MM_ELEM_LIST);
	  elem_list->Clear();	
      
	  int mref = HaAtom::ATOMREF_FULL;
		
	  if(sel_type_str == "Dihedrals")
	  {
		  int idx = 0;
          nel = ptr_mm_mod->p_mm_model->Dihedrals.size();
		  for( i=0; i < nel; i++)
		  {
			  MMDihedral* pdih = ptr_mm_mod->p_mm_model->Dihedrals[i].get();
			  if(pdih->pt2->Selected() && pdih->pt3->Selected())
			  {
				  mref = HaAtom::ATOMREF_FULL;
				  std::string rec;
				  pdih->pt1->FillRef(buf); rec += buf + (std::string)" - "; 
				  mref = adjust_ref_mode(pdih->pt1,pdih->pt2);
				  pdih->pt2->FillRef(buf,mref); rec += buf + (std::string)" - "; 
                  mref = adjust_ref_mode(pdih->pt2,pdih->pt3);
				  pdih->pt3->FillRef(buf,mref); rec += buf + (std::string)" - "; 
                  mref = adjust_ref_mode(pdih->pt3,pdih->pt4);
				  pdih->pt4->FillRef(buf,mref); rec += buf; 
				  elem_list->Append(rec.c_str());
				  elem_list->SetClientData( idx , (void*) pdih );
				  idx++;
			  }
		  }
	  }

	  if(sel_type_str == "Improper Dihedrals")
	  {
          nel = ptr_mm_mod->p_mm_model->ImprDihedrals.size();
		  int idx = 0;
		  for( i=0; i < nel; i++)
		  {
			  MMDihedral* pdih = ptr_mm_mod->p_mm_model->ImprDihedrals[i].get();
			  if(pdih->pt2->Selected() || pdih->pt3->Selected())
			  {
				  mref = HaAtom::ATOMREF_FULL;
				  std::string rec;
				  pdih->pt1->FillRef(buf,mref); rec += buf + (std::string)" - "; 
				  mref = adjust_ref_mode(pdih->pt1,pdih->pt2);
				  pdih->pt2->FillRef(buf,mref); rec += buf + (std::string)" - "; 
				  mref = adjust_ref_mode(pdih->pt2,pdih->pt3);
				  pdih->pt3->FillRef(buf,mref); rec += buf + (std::string)" - "; 
				  mref = adjust_ref_mode(pdih->pt3,pdih->pt4);
				  pdih->pt4->FillRef(buf,mref); rec += buf; 
				  elem_list->Append(rec.c_str());
				  elem_list->SetClientData( idx , (void*)pdih );
				  idx++;
			  }
		  }
	  }

	  if(sel_type_str == "Valence Bonds")
	  {
		  int idx = 0;
          set<MMBond, less<MMBond> >::iterator mbitr = ptr_mm_mod->p_mm_model->MBonds.begin();
          
		  for( ; mbitr != ptr_mm_mod->p_mm_model->MBonds.end(); mbitr++ )
		  {
			  MMBond* pbond = (MMBond*) &(*mbitr);
			  if(pbond->pt1->Selected() && pbond->pt2->Selected())
			  {
				  mref = HaAtom::ATOMREF_FULL;
				  std::string rec;
				  pbond->pt1->FillRef(buf,mref); rec += buf + (std::string)" - "; 
				  mref = adjust_ref_mode(pbond->pt1,pbond->pt2);
				  pbond->pt2->FillRef(buf,mref); rec += buf ; 
				  elem_list->Append(rec.c_str());
				  elem_list->SetClientData( idx , (void*)pbond );
				  idx++;
			  }
		  }
	  }	 

	  if(sel_type_str == "Valence Angles")
	  {
		  int idx = 0;
		  set<MMValAngle, less<MMValAngle> >::iterator vaitr = ptr_mm_mod->p_mm_model->ValAngles.begin();
		  for(; vaitr != ptr_mm_mod->p_mm_model->ValAngles.end(); vaitr++ )
		  {
			  MMValAngle* pval = (MMValAngle*) &(*vaitr);
			  if(pval->pt1->Selected() && pval->pt2->Selected() && pval->pt3->Selected())
			  {
				  mref = HaAtom::ATOMREF_FULL;
				  std::string rec;
				  pval->pt1->FillRef(buf,mref); rec += buf + (std::string)" - ";
				  mref = adjust_ref_mode(pval->pt1,pval->pt2);
				  pval->pt2->FillRef(buf,mref); rec += buf + (std::string)" - "; 
				  mref = adjust_ref_mode(pval->pt2,pval->pt3);
				  pval->pt3->FillRef(buf,mref); rec += buf ; 
				  elem_list->Append(rec.c_str());
				  elem_list->SetClientData( idx , (void*)pval );
				  idx++;
			  }
		  }
	  } 
}


void MolMechDlgWX::OnChangeSelElem(wxCommandEvent& event)
{
	  MolSet* pmset = ptr_mm_mod->GetMolSet();
	  pmset->ClearPickedAtoms();

	  cur_bond = NULL;
	  cur_vang = NULL;
	  cur_dih  = NULL;
	  cur_impr_dih = NULL;	  

      wxListBox* elem_list= (wxListBox*) FindWindow(IDC_MM_ELEM_LIST);

      wxAtomEdit* edit_at1 = (wxAtomEdit*) FindWindow(IDC_AT1);
      wxAtomEdit* edit_at2 = (wxAtomEdit*) FindWindow(IDC_AT2);
      wxAtomEdit* edit_at3 = (wxAtomEdit*) FindWindow(IDC_AT3);
      wxAtomEdit* edit_at4 = (wxAtomEdit*) FindWindow(IDC_AT4);

	  edit_at1->SetValue("");
	  edit_at2->SetValue("");
	  edit_at3->SetValue("");
	  edit_at4->SetValue("");

	  wxTextCtrl* edit_fconst =  (wxTextCtrl*) FindWindow(IDC_MM_FCONST);
	  wxTextCtrl* edit_eq_dist = (wxTextCtrl*) FindWindow(IDC_MM_EQ_DIST);
	  wxTextCtrl* edit_curr_dist = (wxTextCtrl*) FindWindow(IDC_MM_CUR_DIST);

	  edit_fconst->Clear();
	  edit_eq_dist->Clear();
	  edit_curr_dist->Clear();

	  wxTextCtrl* edit_at_symbol_1 = (wxTextCtrl*) FindWindow(IDC_MM_AT_SYMBOL_1);
	  wxTextCtrl* edit_at_symbol_2 = (wxTextCtrl*) FindWindow(IDC_MM_AT_SYMBOL_2);
	  wxTextCtrl* edit_at_symbol_3 = (wxTextCtrl*) FindWindow(IDC_MM_AT_SYMBOL_3);
	  wxTextCtrl* edit_at_symbol_4 = (wxTextCtrl*) FindWindow(IDC_MM_AT_SYMBOL_4);
 
	  edit_at_symbol_1->SetValue("");
	  edit_at_symbol_2->SetValue("");
	  edit_at_symbol_3->SetValue("");
	  edit_at_symbol_4->SetValue("");

	  wxRadioBox* elem_type = (wxRadioBox*) FindWindow(IDC_RADIO_ELEMENTS);
  	  wxString sel_type_str = elem_type->GetStringSelection();
	  
	 
	  wxComboBox* set_type_combo = (wxComboBox*) FindWindow(IDC_MM_SET_TYPE); 

	  int idx_sel = elem_list->GetSelection();
	  if(idx_sel == -1) 
		  return;

	  if(sel_type_str.IsEmpty()) return;

	  void* ptr_item = elem_list->GetClientData(idx_sel);
	  if( ptr_item == NULL) return;

	  char buf[256];

	  if( sel_type_str == "Dihedrals" || sel_type_str == "Improper Dihedrals" )
	  {
          MMDihedral* pdih = (MMDihedral*) ptr_item;

		  if(sel_type_str == "Dihedrals") cur_dih = pdih;
		  if(sel_type_str == "Improper Dihedrals") cur_impr_dih = pdih;

          pdih->pt1->FillRef(buf);
		  edit_at1->SetValue(buf);
          pdih->pt2->FillRef(buf);
		  edit_at2->SetValue(buf);
          pdih->pt3->FillRef(buf);
		  edit_at3->SetValue(buf);
          pdih->pt4->FillRef(buf);
		  edit_at4->SetValue(buf);

		  edit_at_symbol_1->SetValue(pdih->pt1->GetFFSymbol());
		  edit_at_symbol_2->SetValue(pdih->pt2->GetFFSymbol());
		  edit_at_symbol_3->SetValue(pdih->pt3->GetFFSymbol());
		  edit_at_symbol_4->SetValue(pdih->pt4->GetFFSymbol());

		  int nterm = pdih->GetNTerms();

		  buf[0] = 0;
		  int i;
		  int j = 0;
		  for( i=0; i < nterm; i++)
		  {
		     j += sprintf(buf+j,"%9.4f",pdih->pk[i]);
			 if( i != (nterm -1)) j += sprintf(buf+j,";");

		  }	

		  edit_fconst->SetValue(buf);
		  
		  buf[0] = 0;
		  j = 0;
		  for( i=0; i < nterm; i++)
		  {
		     j += sprintf(buf+j,"%9.4f",pdih->phase[i]);
			 if( i != (nterm -1)) j += sprintf(buf+j,";");
		  }
		  edit_eq_dist->SetValue(buf);

		  double cur_val = RAD_TO_DEG*Vec3D::CalcTorsion(pdih->pt1,pdih->pt2,pdih->pt3,pdih->pt4);
		  sprintf(buf,"%9.4f",cur_val);
		  edit_curr_dist->SetValue(buf);

		  set_type_combo->SetSelection(pdih->set_type);

		  pmset->picked_atoms.InsertAtom((HaAtom*)pdih->pt1);
		  pmset->picked_atoms.InsertAtom((HaAtom*)pdih->pt2);
		  pmset->picked_atoms.InsertAtom((HaAtom*)pdih->pt3);
		  pmset->picked_atoms.InsertAtom((HaAtom*)pdih->pt4);
	  }

	if(sel_type_str == "Valence Bonds")
	{
		MMBond* pbond = (MMBond*) ptr_item;
		if(pbond == NULL) return;
		cur_bond = pbond;

		pbond->pt1->FillRef(buf);
		edit_at1->SetValue(buf);
		pbond->pt2->FillRef(buf);
		edit_at2->SetValue(buf);

		edit_at_symbol_1->SetValue(pbond->pt1->GetFFSymbol());
		edit_at_symbol_2->SetValue(pbond->pt2->GetFFSymbol());
 
		sprintf(buf,"%9.4f",pbond->fc);
		edit_fconst->SetValue(buf);
		sprintf(buf,"%9.4f",pbond->r0);
		edit_eq_dist->SetValue(buf);

		double cur_val = Vec3D::CalcDistance(pbond->pt1,pbond->pt2);
		sprintf(buf,"%9.4f",cur_val);
		edit_curr_dist->SetValue(buf);

		set_type_combo->SetSelection(pbond->set_type);

		pmset->picked_atoms.InsertAtom((HaAtom*)pbond->pt1);
		pmset->picked_atoms.InsertAtom((HaAtom*)pbond->pt2);
	}

	  if(sel_type_str == "Valence Angles")
	  {
		  MMValAngle* pval = (MMValAngle*) ptr_item;
		  if(pval == NULL) return;
		  cur_vang = pval;

          pval->pt1->FillRef(buf);
		  edit_at1->SetValue(buf);
          pval->pt2->FillRef(buf);
		  edit_at2->SetValue(buf);
          pval->pt3->FillRef(buf);
		  edit_at3->SetValue(buf);

		  edit_at_symbol_1->SetValue(pval->pt1->GetFFSymbol());
		  edit_at_symbol_2->SetValue(pval->pt2->GetFFSymbol());
		  edit_at_symbol_3->SetValue(pval->pt3->GetFFSymbol());

		  sprintf(buf,"%9.4f",pval->fc);
		  edit_fconst->SetValue(buf);
	      sprintf(buf,"%9.4f",pval->a0);
		  edit_eq_dist->SetValue(buf);

		  double cur_val = RAD_TO_DEG*Vec3D::CalcAngle(pval->pt1,pval->pt2,pval->pt3);
		  sprintf(buf,"%9.4f",cur_val);
		  edit_curr_dist->SetValue(buf);

		  set_type_combo->SetSelection(pval->set_type);

		  pmset->picked_atoms.InsertAtom((HaAtom*)pval->pt1);
		  pmset->picked_atoms.InsertAtom((HaAtom*)pval->pt2);
		  pmset->picked_atoms.InsertAtom((HaAtom*)pval->pt3);
 	  }
	
	  pmset->RefreshAllViews(RFRefresh);

}

void MolMechDlgWX::OnSetNewPar(wxCommandEvent& event) 
{
	wxRadioBox* elem_type = (wxRadioBox*) FindWindow(IDC_RADIO_ELEMENTS);
	int isel_type = elem_type->GetSelection();
	wxString elem_type_str =  elem_type->GetStringSelection();

	wxTextCtrl* edit_fconst =  (wxTextCtrl*) FindWindow(IDC_MM_FCONST);
	wxTextCtrl* edit_eq_dist = (wxTextCtrl*) FindWindow(IDC_MM_EQ_DIST);
	wxComboBox* set_type_combo = (wxComboBox*) FindWindow(IDC_MM_SET_TYPE); 

	double fc;
	double eq_dist;
	bool bres;

	wxString str;
	str = edit_fconst->GetValue();
	bres = str.ToDouble(&fc);
	if(!bres)
	{
		PrintLog("Error in MolMechDlgWX::OnSetNewPar() \n");
		PrintLog(" Error reading force constant \n");
		return;
	}

	str = edit_eq_dist->GetValue();
	bres = str.ToDouble(&eq_dist);

	if(!bres)
	{
		PrintLog("Error in MolMechDlgWX::OnSetNewPar() \n");
		PrintLog(" Error reading equilibrium distance \n");
		return;
	}

	int changed = FALSE;

	if(elem_type_str == "Valence Bonds")
	{
		if(cur_bond == NULL) return;

		if( fabs(cur_bond->fc - fc) > 0.001) 
		{
			changed = TRUE;
			cur_bond->fc = fc;
		}

		if( fabs(cur_bond->r0 - eq_dist) > 0.001) 
		{
			changed = TRUE;
			cur_bond->r0 = eq_dist;
		}

		if(changed = TRUE)
		{
			cur_bond->set_type = MolMechModel::SET_SPEC;
			set_type_combo->SetSelection(MolMechModel::SET_SPEC);
		}
	}

	if(elem_type_str == "Valence Angles")
	{
		if(cur_vang == NULL) return;

		if( fabs(cur_vang->fc - fc) > 0.001) 
		{
			changed = TRUE;
			cur_vang->fc = fc;
		}

		if( fabs(cur_vang->a0 - eq_dist) > 0.001) 
		{
			changed = TRUE;
			cur_vang->a0 = eq_dist;
		}

		if(changed = TRUE)
		{
			cur_vang->set_type = MolMechModel::SET_SPEC;
			set_type_combo->SetSelection(MolMechModel::SET_SPEC);
		}
	}


	if(elem_type_str == "Dihedrals")
	{
		if(cur_dih == NULL) return;
	}
	
	if(elem_type_str == "Improper Dihedrals")
	{
		if(cur_impr_dih == NULL) return;
	}

	ptr_mm_mod->p_mm_model->p_amber_model->SetUpdateDataFlag();
}

void MolMechDlgWX::OnDelImprAng(wxCommandEvent& event) 
{
	ptr_mm_mod->p_mm_model->ImprDihedrals.clear();
	OnUpdateElemList(event);
}

void MolMechDlgWX::OnSendMPIMsg1(wxCommandEvent& event)
{
	std::string msg = "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>"; msg += "\n";
	msg+= "<wxCommandEvent>";             msg += "\n";
	msg+= wxString::Format("<type>%d</type>\n", HA_MOL_MECH_EVENT);
	msg+= wxString::Format("<id>%d</id>\n",MOL_MECH_ID_TEST1);            
	msg+= "</wxCommandEvent>";              msg += "\n";
	PrintLog("XML Message: \n%s",msg.c_str());

	pApp->mpi_driver->SendXmlMsgAllProc(msg.c_str());
}

void MolMechDlgWX::OnSetDnaCSPars(wxCommandEvent& event) 
{
	ptr_mm_mod->p_mm_model->SetCoarseGrainedDNAParams();
	TransferDataToWindow();
	OnUpdateElemList(event);
}

void MolMechDlgWX::OnSaveResffFromMort(wxCommandEvent& event)
{
	TransferDataFromWindow();

	std::string fname = "resff_test.xml";
	MolSet* pmset = ptr_mm_mod->GetMolSet();
	ForceFieldType ff_type = MMForceField::ff_type_default;

	MMForceField* p_ff = MMForceField::GetMMForceField( ff_type, TRUE );
	p_ff->SaveResFFTemplatesFromMort(fname.c_str(), pmset);
}

void MolMechDlgWX::OnSaveAtomIndRestrArb(wxCommandEvent& event)
{
	TransferDataFromWindow();

	MolMechModel* p_model = ptr_mm_mod->GetMolMechModel();
	p_model->SaveAtomRestrArbalestIndForm("atom_restr_rules.xml","atom_restr_list.xml");
}

void MolMechDlgWX::OnInitMolMech(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	ptr_mm_mod->Initialize();
	OnUpdateElemList(event);
	TransferDataToWindow();
}

void MolMechDlgWX::OnRestrainHBonds(wxCommandEvent& event) 
{
	wxTextCtrl* ptr_ffpar_frccnst= (wxTextCtrl*) FindWindow(IDC_RESTRAIN_HBONDS);
	wxString str = ptr_ffpar_frccnst->GetValue();
	double dval;
	bool bres;
	bres = str.ToDouble(&dval);

	if(bres)
	{
	   ptr_mm_mod->p_mm_model->SetHBondRestraints(dval);
	}
	else
	{
		ptr_ffpar_frccnst->SetValue("0.0");
		PrintLog("Error Reading HBOnd Constraint force constant");
	}
}

MMInfoDlg::MMInfoDlg(HaMolMechMod* ptr_mm_mod_new, wxWindow *parent ):
wxFrame( parent, -1, "MM Info")
{
	this->SetExtraStyle(wxWS_EX_VALIDATE_RECURSIVELY);
	ptr_mm_mod = ptr_mm_mod_new;
	p_mm_info  = ptr_mm_mod->p_mm_info;
	ene_format = "%12.4f";

	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);

	mm_info_dlg( this, TRUE );

//    wxMenuBar* edit_groups_menu_bar = edit_groups_menu();
//    SetMenuBar(edit_groups_menu_bar);    

	OnInitDialog();
}

MMInfoDlg::~MMInfoDlg()
{
     
}


BEGIN_EVENT_TABLE(MMInfoDlg,wxFrame)
	EVT_BUTTON(IDC_MM_INFO_UPDATE, MMInfoDlg::OnUpdate )
	EVT_BUTTON(IDC_MM_INFO_UPDATE_FORMAT, MMInfoDlg::OnUpdateFormat )
	EVT_BUTTON(IDC_MM_INFO_CLOSE, MMInfoDlg::OnCloseBtn )
	EVT_CLOSE(MMInfoDlg::OnClose)
END_EVENT_TABLE()

void  MMInfoDlg::OnInitDialog()
{
	wxTextCtrl* edit_ctrl;

	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_MM_INFO_DISPLAY_FORMAT); 
	edit_ctrl->SetValidator( wxGenericValidator(&ene_format) );

	DDX_Text_int (this, IDC_MM_INFO_STEP, p_mm_info->nstep );

	DDX_Text_double(this, IDC_MM_INFO_TOT_ENE, p_mm_info->tot_energy, ene_format );
	AddWindowToArray(this,IDC_MM_INFO_TOT_ENE, dbl_text_boxes);

	DDX_Text_double(this, IDC_MM_INFO_POT_ENE, p_mm_info->pot_ene, ene_format );
	AddWindowToArray(this,IDC_MM_INFO_POT_ENE, dbl_text_boxes);

	DDX_Text_double(this, IDC_MM_INFO_ELECTR_ENE, p_mm_info->electr_ene, ene_format );
	AddWindowToArray(this,IDC_MM_INFO_ELECTR_ENE, dbl_text_boxes);

	DDX_Text_double(this, IDC_MM_INFO_ELECTR_14_ENE, p_mm_info->electr_ene_14, ene_format );
	AddWindowToArray(this,IDC_MM_INFO_ELECTR_14_ENE, dbl_text_boxes);

	DDX_Text_double(this, IDC_MM_INFO_KIN_ENE, p_mm_info->kin_ene, ene_format );
	AddWindowToArray(this,IDC_MM_INFO_KIN_ENE, dbl_text_boxes);

	DDX_Text_double(this, IDC_MM_INFO_VDW_ENE, p_mm_info->vdw_ene, ene_format );
	AddWindowToArray(this,IDC_MM_INFO_VDW_ENE, dbl_text_boxes);

	DDX_Text_double(this, IDC_MM_INFO_VDW_14_ENE, p_mm_info->vdw_ene_14, ene_format );
	AddWindowToArray(this,IDC_MM_INFO_VDW_14_ENE, dbl_text_boxes);

	DDX_Text_double(this, IDC_MM_INFO_TEMP, p_mm_info->temp, ene_format );
	AddWindowToArray(this,IDC_MM_INFO_TEMP, dbl_text_boxes);

	DDX_Text_double(this, IDC_MM_INFO_BOND_ENE, p_mm_info->bond_ene, ene_format );
	AddWindowToArray(this,IDC_MM_INFO_BOND_ENE, dbl_text_boxes);

	DDX_Text_double(this, IDC_MM_INFO_VANG_ENE, p_mm_info->vang_ene, ene_format );
	AddWindowToArray(this,IDC_MM_INFO_VANG_ENE, dbl_text_boxes);

	DDX_Text_double(this, IDC_MM_INFO_DIH_ENE, p_mm_info->dihed_ene, ene_format );
	AddWindowToArray(this,IDC_MM_INFO_DIH_ENE, dbl_text_boxes);

	DDX_Text_double(this, IDC_MM_INFO_GB_ENE, p_mm_info->gb_ene, ene_format );
	AddWindowToArray(this,IDC_MM_INFO_GB_ENE, dbl_text_boxes);

	DDX_Text_double(this, IDC_MM_INFO_HBOND_ENE, p_mm_info->hbond_ene, ene_format );
	AddWindowToArray(this,IDC_MM_INFO_HBOND_ENE, dbl_text_boxes);

	DDX_Text_double(this, IDC_MM_INFO_PRESS, p_mm_info->press, ene_format );
	AddWindowToArray(this,IDC_MM_INFO_PRESS, dbl_text_boxes);

	DDX_Text_double(this, IDC_MM_INFO_HARM_CONSTR_ENE, p_mm_info->constraints_ene, ene_format );
	AddWindowToArray(this,IDC_MM_INFO_HARM_CONSTR_ENE, dbl_text_boxes);

	DDX_Text_double(this, IDC_MM_INFO_ELECTR_POLAR_ENE, p_mm_info->epol, ene_format );
	AddWindowToArray(this,IDC_MM_INFO_ELECTR_POLAR_ENE, dbl_text_boxes);

	DDX_Text_double(this, IDC_MM_INFO_DENSITY, p_mm_info->density, ene_format );
	AddWindowToArray(this,IDC_MM_INFO_DENSITY, dbl_text_boxes);


}

bool MMInfoDlg::TransferDataToWindow()
{
	return wxFrame::TransferDataToWindow();
}

bool MMInfoDlg::TransferDataFromWindow()
{
	return wxFrame::TransferDataFromWindow();
}

void MMInfoDlg::OnUpdateFormat(wxCommandEvent& event)
{
	wxTextCtrl* edit_ctrl;
	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_MM_INFO_DISPLAY_FORMAT);
	ene_format = edit_ctrl->GetValue();

	int n = dbl_text_boxes.size();
	int i;
	for(i = 0; i < n; i++)
	{
		 wxDoubleValidator* p_val = (wxDoubleValidator*) dbl_text_boxes[i]->GetValidator();
		 p_val->format = ene_format;
	}
	TransferDataToWindow();
}

void MMInfoDlg::OnUpdate(wxCommandEvent& event)
{
	TransferDataToWindow();
}

void MMInfoDlg::OnClose(wxCloseEvent& event)
{
	this->Show(FALSE);
}

void MMInfoDlg::OnCloseBtn(wxCommandEvent& event)
{
	this->Show(FALSE);
}
