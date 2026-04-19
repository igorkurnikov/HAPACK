/*! \file mm_dialogs_wx.cpp

    Dialogs for Molecular Mechanics and Related Modules
 
    \author Igor Kurnikov  
    \date 2010-
*/

#define MM_DIALOGS_WX_CPP

#include "hampi.h"
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

// ---------------------------------------------------------------------------
// Page builder functions relocated from ha_wx_res_wdr.cpp (wxDesigner-generated).
// Each function constructs the controls of one notebook tab.
// ---------------------------------------------------------------------------

wxSizer *mol_mech_dlg( wxWindow *parent, bool call_fit, bool set_sizer )
{
    wxBoxSizer *item0 = new wxBoxSizer( wxVERTICAL );

    wxBoxSizer *item1 = new wxBoxSizer( wxHORIZONTAL );

    wxNotebook *item3 = new wxNotebook( parent, ID_MM_DLG, wxDefaultPosition, wxDefaultSize, 0 );
#if !wxCHECK_VERSION(2,5,2)
    wxNotebookSizer *item2 = new wxNotebookSizer( item3 );
#else
    wxWindow *item2 = item3;
#endif

    wxPanel *item4 = new wxPanel( item3, -1 );
    mm_run_param_page( item4, FALSE );
    item3->AddPage( item4, wxT("MM Run parameters") );

    wxPanel *item5 = new wxPanel( item3, -1 );
    mm_par_setup_page( item5, FALSE );
    item3->AddPage( item5, wxT("Force Field Parameters") );

    wxPanel *item6 = new wxPanel( item3, -1 );
    mm_edit_model_page( item6, FALSE );
    item3->AddPage( item6, wxT("Edit MM Model") );

    wxPanel *item7 = new wxPanel( item3, -1 );
    mm_md_anal_page( item7, FALSE );
    item3->AddPage( item7, wxT("MD analysis") );

    item1->Add( item2, 1, wxGROW|wxALIGN_CENTER_HORIZONTAL|wxALL, 5 );

    item0->Add( item1, 1, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    if (set_sizer)
    {
        parent->SetSizer( item0 );
        if (call_fit)
            item0->SetSizeHints( parent );
    }

    return item0;
}

wxSizer *mm_edit_model_page( wxWindow *parent, bool call_fit, bool set_sizer )
{
    wxBoxSizer *item0 = new wxBoxSizer( wxVERTICAL );

    wxBoxSizer *item1 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item2 = new wxStaticText( parent, ID_TEXT, wxT("Molecular Mechanics Model Element List(Atoms, Bonds, Angles ):"), wxDefaultPosition, wxDefaultSize, 0 );
    item1->Add( item2, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5 );

    item0->Add( item1, 0, wxALIGN_CENTER|wxALL|wxSHAPED, 5 );

    wxBoxSizer *item3 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item4 = new wxStaticText( parent, ID_TEXT, wxT("Default Force Field Type:"), wxDefaultPosition, wxDefaultSize, 0 );
    item3->Add( item4, 0, wxALIGN_CENTER|wxALL, 5 );

    wxString *strs5 = (wxString*) NULL;
    wxComboBox *item5 = new wxComboBox( parent, IDC_MM_FF_TYPE_DEFAULT, wxT(""), wxDefaultPosition, wxSize(200,-1), 0, strs5, wxCB_DROPDOWN );
    item3->Add( item5, 0, wxALIGN_CENTER|wxALL, 5 );

    wxStaticText *item6 = new wxStaticText( parent, ID_TEXT, wxT("Current Force Field of the model:"), wxDefaultPosition, wxDefaultSize, 0 );
    item3->Add( item6, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item7 = new wxTextCtrl( parent, IDC_MM_FF_TYPE, wxT(""), wxDefaultPosition, wxSize(80,-1), wxTE_READONLY );
    item3->Add( item7, 1, wxALIGN_CENTER|wxALL, 5 );

    item0->Add( item3, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item8 = new wxBoxSizer( wxHORIZONTAL );

    wxBoxSizer *item9 = new wxBoxSizer( wxVERTICAL );

    wxArrayString strs10;
    strs10.Add(wxT("Valence Bonds")); 
    strs10.Add(wxT("Valence Angles"));
    strs10.Add(wxT("Dihedrals"));
    strs10.Add(wxT("Improper Dihedrals"));
    strs10.Add(wxT("Atoms"));
    wxRadioBox *item10 = new wxRadioBox( parent, IDC_RADIO_ELEMENTS, wxT("Elements type"), wxDefaultPosition, wxDefaultSize, strs10, 1, wxRA_SPECIFY_COLS );
    item9->Add( item10, 0, wxALIGN_CENTER|wxALL, 5 );

    item8->Add( item9, 0, wxGROW|wxALIGN_CENTER_HORIZONTAL|wxALL, 5 );

    wxString *strs11 = (wxString*) NULL;
    wxListBox *item11 = new wxListBox( parent, IDC_MM_ELEM_LIST, wxDefaultPosition, wxSize(400,150), 0, strs11, wxLB_SINGLE );
    item8->Add( item11, 1, wxGROW|wxALIGN_CENTER_HORIZONTAL|wxALL, 5 );

    item0->Add( item8, 1, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item12 = new wxBoxSizer( wxHORIZONTAL );

    wxBoxSizer *item13 = new wxBoxSizer( wxVERTICAL );

    wxStaticText *item14 = new wxStaticText( parent, ID_TEXT, wxT("Atoms of the selected element"), wxDefaultPosition, wxDefaultSize, 0 );
    item13->Add( item14, 0, wxALIGN_CENTER|wxALL, 5 );

    wxBoxSizer *item15 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item16 = new wxStaticText( parent, ID_TEXT, wxT("ATOM 1:"), wxDefaultPosition, wxDefaultSize, 0 );
    item15->Add( item16, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item17 = new wxTextCtrl( parent, IDC_AT1, wxT(""), wxDefaultPosition, wxSize(160,-1), 0 );
    item15->Add( item17, 1, wxALIGN_CENTER|wxALL, 5 );

    item13->Add( item15, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item18 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item19 = new wxStaticText( parent, ID_TEXT, wxT("ATOM 2:"), wxDefaultPosition, wxDefaultSize, 0 );
    item18->Add( item19, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item20 = new wxTextCtrl( parent, IDC_AT2, wxT(""), wxDefaultPosition, wxSize(160,-1), 0 );
    item18->Add( item20, 1, wxALIGN_CENTER|wxALL, 5 );

    item13->Add( item18, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item21 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item22 = new wxStaticText( parent, ID_TEXT, wxT("ATOM 3:"), wxDefaultPosition, wxDefaultSize, 0 );
    item21->Add( item22, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item23 = new wxTextCtrl( parent, IDC_AT3, wxT(""), wxDefaultPosition, wxSize(160,-1), 0 );
    item21->Add( item23, 1, wxALIGN_CENTER|wxALL, 5 );

    item13->Add( item21, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item24 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item25 = new wxStaticText( parent, ID_TEXT, wxT("ATOM 4:"), wxDefaultPosition, wxDefaultSize, 0 );
    item24->Add( item25, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item26 = new wxTextCtrl( parent, IDC_AT4, wxT(""), wxDefaultPosition, wxSize(160,-1), 0 );
    item24->Add( item26, 1, wxALIGN_CENTER|wxALL, 5 );

    item13->Add( item24, 1, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    item12->Add( item13, 1, wxGROW|wxALIGN_CENTER_HORIZONTAL|wxALL, 5 );

    wxBoxSizer *item27 = new wxBoxSizer( wxVERTICAL );

    wxStaticText *item28 = new wxStaticText( parent, ID_TEXT, wxT("Element FF params:"), wxDefaultPosition, wxDefaultSize, 0 );
    item27->Add( item28, 0, wxALIGN_CENTER|wxALL, 5 );

    wxBoxSizer *item29 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item30 = new wxStaticText( parent, ID_TEXT, wxT("force constant:"), wxDefaultPosition, wxDefaultSize, 0 );
    item29->Add( item30, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item31 = new wxTextCtrl( parent, IDC_MM_FCONST, wxT(""), wxDefaultPosition, wxDefaultSize, 0 );
    item29->Add( item31, 0, wxALIGN_CENTER|wxALL, 5 );

    wxStaticText *item32 = new wxStaticText( parent, IDC_FF_UNIT_LABEL, wxT("(kcal/A^2)"), wxDefaultPosition, wxDefaultSize, 0 );
    item29->Add( item32, 0, wxALIGN_CENTER|wxALL, 5 );

    item27->Add( item29, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item33 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item34 = new wxStaticText( parent, ID_TEXT, wxT("Eq dist (angle)"), wxDefaultPosition, wxDefaultSize, 0 );
    item33->Add( item34, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item35 = new wxTextCtrl( parent, IDC_MM_EQ_DIST, wxT(""), wxDefaultPosition, wxDefaultSize, 0 );
    item33->Add( item35, 0, wxALIGN_CENTER|wxALL, 5 );

    wxStaticText *item36 = new wxStaticText( parent, ID_TEXT, wxT("Ang"), wxDefaultPosition, wxDefaultSize, 0 );
    item33->Add( item36, 0, wxALIGN_CENTER|wxALL, 5 );

    item27->Add( item33, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item37 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item38 = new wxStaticText( parent, ID_TEXT, wxT("Curr Dist(Angle):"), wxDefaultPosition, wxDefaultSize, 0 );
    item37->Add( item38, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item39 = new wxTextCtrl( parent, IDC_MM_CUR_DIST, wxT(""), wxDefaultPosition, wxDefaultSize, 0 );
    item37->Add( item39, 0, wxALIGN_CENTER|wxALL, 5 );

    wxStaticText *item40 = new wxStaticText( parent, ID_TEXT, wxT("Ang"), wxDefaultPosition, wxDefaultSize, 0 );
    item37->Add( item40, 0, wxALIGN_CENTER|wxALL, 5 );

    item27->Add( item37, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item41 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item42 = new wxStaticText( parent, ID_TEXT, wxT("Param Set Type:"), wxDefaultPosition, wxDefaultSize, 0 );
    item41->Add( item42, 0, wxALIGN_CENTER|wxALL, 5 );

    wxArrayString strs43;
    strs43.Add(wxT("NOT_SET"));
    strs43.Add(wxT("DEFAULT VAL"));
    strs43.Add(wxT("FROM FF FIELD"));
    strs43.Add(wxT("FROM RES TEMPL"));
    strs43.Add(wxT("SPECIAL VALUE"));
    wxComboBox *item43 = new wxComboBox( parent, IDC_MM_SET_TYPE, wxT(""), wxDefaultPosition, wxDefaultSize, strs43, wxCB_DROPDOWN );
    item41->Add( item43, 0, wxALIGN_CENTER|wxALL, 5 );

    item27->Add( item41, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    item12->Add( item27, 0, wxGROW|wxALIGN_CENTER_HORIZONTAL|wxALL|wxSHAPED, 5 );

    wxBoxSizer *item44 = new wxBoxSizer( wxVERTICAL );

    wxStaticText *item45 = new wxStaticText( parent, ID_TEXT,
        wxT("FF Symbols of\n")
        wxT("atoms of the element:"),
        wxDefaultPosition, wxDefaultSize, wxALIGN_CENTRE );
    item44->Add( item45, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item46 = new wxTextCtrl( parent, IDC_MM_AT_SYMBOL_1, wxT(""), wxDefaultPosition, wxDefaultSize, 0 );
    item44->Add( item46, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item47 = new wxTextCtrl( parent, IDC_MM_AT_SYMBOL_2, wxT(""), wxDefaultPosition, wxDefaultSize, 0 );
    item44->Add( item47, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item48 = new wxTextCtrl( parent, IDC_MM_AT_SYMBOL_3, wxT(""), wxDefaultPosition, wxDefaultSize, 0 );
    item44->Add( item48, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item49 = new wxTextCtrl( parent, IDC_MM_AT_SYMBOL_4, wxT(""), wxDefaultPosition, wxDefaultSize, 0 );
    item44->Add( item49, 0, wxALIGN_CENTER|wxALL, 5 );

    item12->Add( item44, 0, wxGROW|wxALIGN_CENTER_HORIZONTAL|wxALL|wxSHAPED, 5 );

    item0->Add( item12, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item50 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item51 = new wxButton( parent, IDC_FFPAR_INIT_MOLMECH, wxT("Init MM Model"), wxDefaultPosition, wxDefaultSize, 0 );
    item50->Add( item51, 0, wxALIGN_CENTER|wxALL, 5 );

    wxStaticLine *item52 = new wxStaticLine( parent, ID_LINE_MM_1, wxDefaultPosition, wxSize(-1,20), wxLI_VERTICAL );
    item50->Add( item52, 0, wxALIGN_CENTER|wxALL, 5 );

    wxButton *item53 = new wxButton( parent, IDC_MM_UPDATE_ELEM_LIST, wxT("Update Model Elements List"), wxDefaultPosition, wxDefaultSize, 0 );
    item50->Add( item53, 0, wxALIGN_CENTER|wxALL, 5 );

    wxButton *item54 = new wxButton( parent, IDC_RESPAR_SET_NEW_PAR, wxT("Set New Params"), wxDefaultPosition, wxDefaultSize, 0 );
    item50->Add( item54, 0, wxALIGN_CENTER|wxALL, 5 );

    wxStaticLine *item55 = new wxStaticLine( parent, ID_LINE_SP_DNA, wxDefaultPosition, wxSize(-1,20), wxLI_VERTICAL );
    item50->Add( item55, 0, wxALIGN_CENTER|wxALL, 5 );

    item50->Add( 20, 20, 0, wxALIGN_CENTER|wxALL, 5 );

    item50->Add( 20, 20, 0, wxALIGN_CENTER|wxALL, 5 );

    wxButton *item56 = new wxButton( parent, IDC_SET_DNA_CS_PARS, wxT("Set DNA CG Pars"), wxDefaultPosition, wxDefaultSize, 0 );
    item50->Add( item56, 0, wxALIGN_CENTER|wxALL, 5 );

    item0->Add( item50, 0, wxALIGN_CENTER|wxALL, 5 );

    wxBoxSizer *item57 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item58 = new wxButton( parent, IDC_DEL_IMPROPER_ANG, wxT("Delete Improper Angles"), wxDefaultPosition, wxDefaultSize, 0 );
    item57->Add( item58, 0, wxALIGN_CENTER|wxALL, 5 );

    wxCheckBox *item59 = new wxCheckBox( parent, IDC_MM_USE_MORT, wxT("use MORT to set MM Model"), wxDefaultPosition, wxDefaultSize, 0 );
    item57->Add( item59, 0, wxALIGN_CENTER|wxALL, 5 );

    item0->Add( item57, 0, wxALIGN_CENTER|wxALL, 5 );

    if (set_sizer)
    {
        parent->SetSizer( item0 );
        if (call_fit)
            item0->SetSizeHints( parent );
    }

    return item0;
}

wxSizer *mm_par_setup_page( wxWindow *parent, bool call_fit, bool set_sizer )
{
    wxBoxSizer *item0 = new wxBoxSizer( wxVERTICAL );

    wxBoxSizer *item1 = new wxBoxSizer( wxHORIZONTAL );

    wxBoxSizer *item2 = new wxBoxSizer( wxVERTICAL );

    wxStaticText *item3 = new wxStaticText( parent, ID_TEXT, wxT("Electrostatics:"), wxDefaultPosition, wxDefaultSize, 0 );
    item2->Add( item3, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxString *strs4 = (wxString*) NULL;
    wxChoice *item4 = new wxChoice( parent, IDC_MM_ELECTR_MODEL, wxDefaultPosition, wxSize(200,-1), 0, strs4, 0 );
    item2->Add( item4, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item5 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item6 = new wxStaticText( parent, ID_TEXT, wxT("Dielectric Constant:"), wxDefaultPosition, wxDefaultSize, 0 );
    item5->Add( item6, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item7 = new wxTextCtrl( parent, IDC_MM_DIEL_CONST, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item5->Add( item7, 1, wxALIGN_CENTER|wxALL, 5 );

    item2->Add( item5, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item8 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item9 = new wxStaticText( parent, ID_TEXT, wxT("Non-Bond Cutoff (Ang):"), wxDefaultPosition, wxDefaultSize, 0 );
    item8->Add( item9, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item10 = new wxTextCtrl( parent, IDC_MM_NONB_CUTOFF_DIST, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item8->Add( item10, 1, wxALIGN_CENTER|wxALL, 0 );

    item2->Add( item8, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item11 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item12 = new wxStaticText( parent, ID_TEXT, wxT("Ionic Strength (M):"), wxDefaultPosition, wxDefaultSize, 0 );
    item11->Add( item12, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item13 = new wxTextCtrl( parent, IDC_MM_ION_STRENGTH, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item11->Add( item13, 1, wxALIGN_CENTER|wxALL, 5 );

    item2->Add( item11, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item14 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item15 = new wxStaticText( parent, ID_TEXT, wxT("Scale 1-4 Electr:"), wxDefaultPosition, wxDefaultSize, 0 );
    item14->Add( item15, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item16 = new wxTextCtrl( parent, IDC_MM_SCALE_14_ELECTR, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item14->Add( item16, 0, wxALIGN_CENTER|wxALL, 5 );

    item2->Add( item14, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item17 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item18 = new wxStaticText( parent, ID_TEXT, wxT("Scale 1-4 VdW:"), wxDefaultPosition, wxDefaultSize, 0 );
    item17->Add( item18, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item19 = new wxTextCtrl( parent, IDC_MM_SCALE_14_VDW, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item17->Add( item19, 0, wxALIGN_CENTER|wxALL, 5 );

    item2->Add( item17, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item20 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item21 = new wxStaticText( parent, ID_TEXT, wxT("AMOEBA params:"), wxDefaultPosition, wxDefaultSize, 0 );
    item21->SetFont( wxFont( 12, wxROMAN, wxNORMAL, wxBOLD ) );
    item20->Add( item21, 0, wxALIGN_CENTER|wxALL, 5 );

    item2->Add( item20, 0, wxALIGN_CENTER|wxALL, 0 );

    wxBoxSizer *item22 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item23 = new wxStaticText( parent, ID_TEXT, wxT("Thole's dumping const:"), wxDefaultPosition, wxDefaultSize, 0 );
    item22->Add( item23, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item24 = new wxTextCtrl( parent, IDC_MM_THOLE_DUMP_CONST, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item22->Add( item24, 1, wxALIGN_CENTER|wxALL, 5 );

    item2->Add( item22, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    item1->Add( item2, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5 );

    wxBoxSizer *item25 = new wxBoxSizer( wxVERTICAL );

    wxBoxSizer *item26 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item27 = new wxStaticText( parent, ID_TEXT, wxT("Atom Position Restraints:"), wxDefaultPosition, wxDefaultSize, 0 );
    item27->SetFont( wxFont( 12, wxROMAN, wxNORMAL, wxBOLD ) );
    item26->Add( item27, 1, wxALIGN_CENTER|wxALL, 5 );

    item25->Add( item26, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item28 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item29 = new wxButton( parent, IDC_MM_CHOOSE_MOVING_ATOMS, wxT("Moving Atoms:"), wxDefaultPosition, wxDefaultSize, 0 );
    item28->Add( item29, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item30 = new wxTextCtrl( parent, IDC_MM_MOVING_ATOMS, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item28->Add( item30, 1, wxALIGN_CENTER|wxALL, 5 );

    item25->Add( item28, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item31 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item32 = new wxButton( parent, IDC_MM_CHOOSE_RESTRAINED_ATOMS, wxT("Restrained Atoms:"), wxDefaultPosition, wxDefaultSize, 0 );
    item31->Add( item32, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item33 = new wxTextCtrl( parent, IDC_RESTRAINED_ATOMS, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item31->Add( item33, 1, wxALIGN_CENTER|wxALL, 5 );

    item25->Add( item31, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item34 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item35 = new wxButton( parent, IDC_MM_CHOOSE_RESTR_REF_CRD, wxT("Restrained Atoms Coordinates:"), wxDefaultPosition, wxDefaultSize, 0 );
    item34->Add( item35, 0, wxALIGN_CENTER|wxALL, 5 );

    wxString *strs36 = (wxString*) NULL;
    wxComboBox *item36 = new wxComboBox( parent, IDC_MM_RESTR_REF_CRD_TYPE, wxT(""), wxDefaultPosition, wxSize(100,-1), 0, strs36, wxCB_DROPDOWN );
    item34->Add( item36, 1, wxALIGN_CENTER|wxALL, 5 );

    item25->Add( item34, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item37 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item38 = new wxStaticText( parent, ID_TEXT, wxT("Restraints Force Constant:"), wxDefaultPosition, wxDefaultSize, 0 );
    item37->Add( item38, 0, wxALIGN_CENTER|wxALL, 10 );

    wxTextCtrl *item39 = new wxTextCtrl( parent, IDC_MM_RESTR_FRC_CONST, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item37->Add( item39, 0, wxALIGN_CENTER|wxALL, 0 );

    item25->Add( item37, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxStaticText *item40 = new wxStaticText( parent, ID_TEXT, wxT("Atom-Atom Distance Restraints:"), wxDefaultPosition, wxDefaultSize, 0 );
    item40->SetFont( wxFont( 12, wxROMAN, wxNORMAL, wxBOLD ) );
    item25->Add( item40, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item41 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item42 = new wxButton( parent, IDC_MM_LOAD_HARM_CONSTR_FILE, wxT("Load Distance Restraints File"), wxDefaultPosition, wxDefaultSize, 0 );
    item41->Add( item42, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxStaticText *item43 = new wxStaticText( parent, ID_TEXT, wxT("Num Restriants:"), wxDefaultPosition, wxDefaultSize, 0 );
    item41->Add( item43, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item44 = new wxTextCtrl( parent, IDC_MM_NUM_HARM_CONSTR, wxT(""), wxDefaultPosition, wxSize(80,-1), wxTE_READONLY );
    item41->Add( item44, 0, wxALIGN_CENTER|wxALL, 5 );

    item25->Add( item41, 0, wxALIGN_CENTER|wxALL, 5 );

    wxButton *item45 = new wxButton( parent, IDC_RESTRAIN_HBONDS, wxT("Set Distance Restraints for H-bonds"), wxDefaultPosition, wxDefaultSize, 0 );
    item25->Add( item45, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item46 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item47 = new wxStaticText( parent, ID_TEXT, wxT("Restraints Force Constant:"), wxDefaultPosition, wxDefaultSize, 0 );
    item46->Add( item47, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item48 = new wxTextCtrl( parent, IDC_CONSTRAIN_FRCCONST, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item46->Add( item48, 0, wxALIGN_CENTER|wxALL, 5 );

    item25->Add( item46, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxStaticText *item49 = new wxStaticText( parent, ID_TEXT, wxT("Valence Bonds Constraints:"), wxDefaultPosition, wxDefaultSize, 0 );
    item49->SetFont( wxFont( 12, wxROMAN, wxNORMAL, wxBOLD ) );
    item25->Add( item49, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item50 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item51 = new wxStaticText( parent, ID_TEXT, wxT("SHAKE Bond Constraints:"), wxDefaultPosition, wxDefaultSize, 0 );
    item50->Add( item51, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxString *strs52 = (wxString*) NULL;
    wxComboBox *item52 = new wxComboBox( parent, IDC_MM_SHAKE_METHOD, wxT(""), wxDefaultPosition, wxDefaultSize, 0, strs52, wxCB_DROPDOWN );
    item50->Add( item52, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    item25->Add( item50, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    item1->Add( item25, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5 );

    item0->Add( item1, 1, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    if (set_sizer)
    {
        parent->SetSizer( item0 );
        if (call_fit)
            item0->SetSizeHints( parent );
    }

    return item0;
}

wxSizer *mm_run_param_page( wxWindow *parent, bool call_fit, bool set_sizer )
{
    wxBoxSizer *item0 = new wxBoxSizer( wxVERTICAL );

    wxBoxSizer *item1 = new wxBoxSizer( wxHORIZONTAL );

    wxBoxSizer *item2 = new wxBoxSizer( wxVERTICAL );

    wxBoxSizer *item3 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item4 = new wxStaticText( parent, ID_TEXT, wxT("Type of Calculations"), wxDefaultPosition, wxDefaultSize, 0 );
    item3->Add( item4, 0, wxALIGN_CENTER|wxALL, 5 );

    wxString *strs5 = (wxString*) NULL;
    wxChoice *item5 = new wxChoice( parent, IDC_MM_RUN_TYPE, wxDefaultPosition, wxDefaultSize, 0, strs5, 0 );
    item3->Add( item5, 1, wxALIGN_CENTER|wxALL, 0 );

    item2->Add( item3, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxButton *item6 = new wxButton( parent, IDC_MM_CALC_SP_ENE, wxT("Calculate Single Point Energy"), wxDefaultPosition, wxDefaultSize, 0 );
    item2->Add( item6, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item7 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item8 = new wxButton( parent, IDC_MM_RUN_CALC, wxT("Run MM Calculations"), wxDefaultPosition, wxDefaultSize, 0 );
    item7->Add( item8, 1, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxCheckBox *item9 = new wxCheckBox( parent, IDC_MM_RUN_INT, wxT("Run Internally"), wxDefaultPosition, wxDefaultSize, 0 );
    item7->Add( item9, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    item2->Add( item7, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item10 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item11 = new wxStaticText( parent, ID_TEXT, wxT("External Program:"), wxDefaultPosition, wxDefaultSize, 0 );
    item10->Add( item11, 0, wxALIGN_CENTER|wxALL, 5 );

    wxString *strs12 = (wxString*) NULL;
    wxChoice *item12 = new wxChoice( parent, IDC_MM_EXT_PROG, wxDefaultPosition, wxSize(100,-1), 0, strs12, 0 );
    item10->Add( item12, 1, wxALIGN_CENTER|wxALL, 5 );

    item2->Add( item10, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxButton *item13 = new wxButton( parent, IDC_SAVE_EXT_PROG_INP, wxT("Save Input Files For External Program"), wxDefaultPosition, wxDefaultSize, 0 );
    item2->Add( item13, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxCheckBox *item14 = new wxCheckBox( parent, IDC_MM_RESTART_FLAG, wxT("Restart Job"), wxDefaultPosition, wxDefaultSize, 0 );
    item2->Add( item14, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item15 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item16 = new wxButton( parent, IDC_MM_STOP, wxT("Stop Calculations"), wxDefaultPosition, wxDefaultSize, 0 );
    item15->Add( item16, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxButton *item17 = new wxButton( parent, IDC_MM_SHOW_INFO, wxT("Show MM Info"), wxDefaultPosition, wxDefaultSize, 0 );
    item15->Add( item17, 0, wxALIGN_CENTER|wxALL, 5 );

    item2->Add( item15, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    item1->Add( item2, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5 );

    wxBoxSizer *item18 = new wxBoxSizer( wxVERTICAL );

    wxStaticText *item19 = new wxStaticText( parent, ID_TEXT_PARAM_MIN, wxT("Parameters for Minimization:"), wxDefaultPosition, wxDefaultSize, 0 );
    item19->SetFont( wxFont( 12, wxROMAN, wxNORMAL, wxBOLD ) );
    item18->Add( item19, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item20 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item21 = new wxStaticText( parent, ID_TEXT_MIN_TYPE, wxT("Minimization Type:"), wxDefaultPosition, wxDefaultSize, 0 );
    item20->Add( item21, 0, wxALIGN_CENTER|wxALL, 5 );

    wxString *strs22 = (wxString*) NULL;
    wxChoice *item22 = new wxChoice( parent, IDC_MM_MIN_TYPE, wxDefaultPosition, wxDefaultSize, 0, strs22, 0 );
    item20->Add( item22, 0, wxALIGN_CENTER|wxALL, 0 );

    item18->Add( item20, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item23 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item24 = new wxStaticText( parent, ID_TEXT_NUM_MIN_STEPS, wxT("Number of Minimization Steps:"), wxDefaultPosition, wxDefaultSize, 0 );
    item23->Add( item24, 0, wxALIGN_CENTER|wxALL, 0 );

    wxTextCtrl *item25 = new wxTextCtrl( parent, IDC_MM_MAX_COMP_CYCLES, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item23->Add( item25, 0, wxALIGN_CENTER|wxALL, 0 );

    wxStaticText *item26 = new wxStaticText( parent, ID_TEXT_NUM_STEEP_DEC_STEPS, wxT("Number of Steepest Descent Steps:"), wxDefaultPosition, wxDefaultSize, 0 );
    item23->Add( item26, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item27 = new wxTextCtrl( parent, IDC_MM_NUM_STEP_STEEP_DESCENT, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item23->Add( item27, 0, wxALIGN_CENTER|wxALL, 0 );

    item18->Add( item23, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item28 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item29 = new wxStaticText( parent, ID_TEXT_CONV_CRT, wxT("Convergence Criterion:"), wxDefaultPosition, wxDefaultSize, 0 );
    item28->Add( item29, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item30 = new wxTextCtrl( parent, IDC_MM_MIN_CNVRG_CRITERIUM, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item28->Add( item30, 0, wxALIGN_CENTER|wxALL, 0 );

    wxStaticText *item31 = new wxStaticText( parent, ID_TEXT_KCAL_MOL_A, wxT("kcal/mole A"), wxDefaultPosition, wxDefaultSize, 0 );
    item28->Add( item31, 0, wxALIGN_CENTER|wxALL, 5 );

    wxStaticText *item32 = new wxStaticText( parent, ID_TEXT_INIT_STEP, wxT("  Initial Step:"), wxDefaultPosition, wxDefaultSize, 0 );
    item28->Add( item32, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item33 = new wxTextCtrl( parent, IDC_MM_INIT_MIN_STEP, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item28->Add( item33, 0, wxALIGN_CENTER|wxALL, 0 );

    item18->Add( item28, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxStaticText *item34 = new wxStaticText( parent, ID_TEXT_MD_PARAMS, wxT("Molecular Dynamics parameters:"), wxDefaultPosition, wxDefaultSize, 0 );
    item34->SetFont( wxFont( 12, wxROMAN, wxNORMAL, wxBOLD ) );
    item18->Add( item34, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item35 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item36 = new wxStaticText( parent, ID_TEXT_MD_STEPS_NUM, wxT("MD steps number "), wxDefaultPosition, wxDefaultSize, 0 );
    item35->Add( item36, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item37 = new wxTextCtrl( parent, IDC_MM_LENGTH_MD_RUN, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item35->Add( item37, 0, wxALIGN_CENTER|wxALL, 0 );

    wxStaticText *item38 = new wxStaticText( parent, ID_TEXT_MD_TIME_STEP, wxT("MD time step (ps):"), wxDefaultPosition, wxDefaultSize, 0 );
    item35->Add( item38, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item39 = new wxTextCtrl( parent, IDC_MM_MD_TIME_STEP, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item35->Add( item39, 0, wxALIGN_CENTER|wxALL, 0 );

    item18->Add( item35, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item40 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item41 = new wxStaticText( parent, ID_TEXT_TEMP_CTRL_METH, wxT("Temperature Control Method:"), wxDefaultPosition, wxDefaultSize, 0 );
    item40->Add( item41, 0, wxALIGN_CENTER|wxALL, 5 );

    wxString *strs42 = (wxString*) NULL;
    wxChoice *item42 = new wxChoice( parent, IDC_TEMP_CONTROL_METHOD, wxDefaultPosition, wxSize(250,-1), 0, strs42, 0 );
    item40->Add( item42, 0, wxALIGN_CENTER|wxALL, 5 );

    item18->Add( item40, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item43 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item44 = new wxStaticText( parent, ID_TEXT_INIT_TEMP, wxT("Initial Temp (K):"), wxDefaultPosition, wxDefaultSize, 0 );
    item43->Add( item44, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item45 = new wxTextCtrl( parent, IDC_MM_INIT_TEMP, wxT(""), wxDefaultPosition, wxDefaultSize, 0 );
    item43->Add( item45, 0, wxALIGN_CENTER|wxALL, 5 );

    wxStaticText *item46 = new wxStaticText( parent, ID_TEXT_REF_TEMP, wxT("Reference Temp (K):"), wxDefaultPosition, wxDefaultSize, 0 );
    item43->Add( item46, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item47 = new wxTextCtrl( parent, IDC_MM_REF_TEMP, wxT(""), wxDefaultPosition, wxDefaultSize, 0 );
    item43->Add( item47, 0, wxALIGN_CENTER|wxALL, 5 );

    wxStaticText *item48 = new wxStaticText( parent, ID_TEXT_LANG_DUMP_CONST, wxT("Langevin Damping  (ps-1)"), wxDefaultPosition, wxDefaultSize, 0 );
    item43->Add( item48, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item49 = new wxTextCtrl( parent, IDC_LANGEVIN_DUMP_CONST, wxT(""), wxDefaultPosition, wxSize(60,-1), 0 );
    item43->Add( item49, 0, wxALIGN_CENTER|wxALL, 5 );

    item18->Add( item43, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item50 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item51 = new wxStaticText( parent, ID_TEXT_INIT_INFO_READ, wxT("Init Info to read:"), wxDefaultPosition, wxDefaultSize, 0 );
    item50->Add( item51, 0, wxALIGN_CENTER|wxALL, 10 );

    wxString *strs52 = (wxString*) NULL;
    wxChoice *item52 = new wxChoice( parent, IDC_MM_INIT_READ_COORD, wxDefaultPosition, wxDefaultSize, 0, strs52, 0 );
    item50->Add( item52, 0, wxALIGN_CENTER|wxALL, 5 );

    wxStaticText *item53 = new wxStaticText( parent, ID_TEXT_START_VEL_METH, wxT("Start velocity method:"), wxDefaultPosition, wxDefaultSize, 0 );
    item50->Add( item53, 0, wxALIGN_CENTER|wxALL, 10 );

    wxString *strs54 = (wxString*) NULL;
    wxChoice *item54 = new wxChoice( parent, IDC_MM_START_VEL_METHOD, wxDefaultPosition, wxDefaultSize, 0, strs54, 0 );
    item50->Add( item54, 0, wxALIGN_CENTER|wxALL, 0 );

    item18->Add( item50, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item55 = new wxBoxSizer( wxHORIZONTAL );

    wxCheckBox *item56 = new wxCheckBox( parent, IDC_MM_REMOVE_INIT_MOTION, wxT("Remove initial motion"), wxDefaultPosition, wxDefaultSize, 0 );
    item55->Add( item56, 0, wxALIGN_CENTER|wxALL, 5 );

    wxCheckBox *item57 = new wxCheckBox( parent, IDC_MM_WRAP_COORD, wxT("Wrap crd to box"), wxDefaultPosition, wxDefaultSize, 0 );
    item55->Add( item57, 0, wxALIGN_CENTER|wxALL, 5 );

    wxStaticText *item58 = new wxStaticText( parent, ID_TEXT_REM_COM_MOTION_FREQ, wxT("Remove COM motion freq:"), wxDefaultPosition, wxDefaultSize, 0 );
    item55->Add( item58, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item59 = new wxTextCtrl( parent, IDC_MM_REMOVE_RB_MOTION_FREQ, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item55->Add( item59, 0, wxALIGN_CENTER|wxALL, 0 );

    item18->Add( item55, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item60 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item61 = new wxStaticText( parent, ID_TEXT, wxT("Boundary Conditions:"), wxDefaultPosition, wxDefaultSize, 0 );
    item61->SetFont( wxFont( 12, wxROMAN, wxNORMAL, wxBOLD ) );
    item60->Add( item61, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxString *strs62 = (wxString*) NULL;
    wxChoice *item62 = new wxChoice( parent, IDC_MM_PERIOD_BCOND, wxDefaultPosition, wxSize(150,-1), 0, strs62, 0 );
    item60->Add( item62, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxStaticText *item63 = new wxStaticText( parent, ID_TEXT_PRESS_CTRL, wxT("Pressure control:"), wxDefaultPosition, wxDefaultSize, 0 );
    item60->Add( item63, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxString *strs64 = (wxString*) NULL;
    wxChoice *item64 = new wxChoice( parent, IDC_MM_PRESSURE_REG_METHOD, wxDefaultPosition, wxDefaultSize, 0, strs64, 0 );
    item60->Add( item64, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    item18->Add( item60, 1, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item65 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item66 = new wxButton( parent, IDC_MM_EDIT_PERIODIC_BOX, wxT("Edit Periodic Box"), wxDefaultPosition, wxDefaultSize, 0 );
    item65->Add( item66, 0, wxALIGN_CENTER|wxALL, 5 );

    item18->Add( item65, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    item1->Add( item18, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5 );

    item0->Add( item1, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item67 = new wxBoxSizer( wxHORIZONTAL );

    wxBoxSizer *item68 = new wxBoxSizer( wxVERTICAL );

    wxBoxSizer *item69 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item70 = new wxStaticText( parent, ID_TEXT, wxT("Input/Output files:"), wxDefaultPosition, wxDefaultSize, 0 );
    item70->SetFont( wxFont( 12, wxROMAN, wxNORMAL, wxBOLD ) );
    item69->Add( item70, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    item68->Add( item69, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item71 = new wxBoxSizer( wxHORIZONTAL );

    item71->Add( 80, 20, 2, wxALIGN_CENTER|wxALL, 0 );

    wxStaticText *item72 = new wxStaticText( parent, ID_TEXT, wxT("Input Files:    "), wxDefaultPosition, wxDefaultSize, wxALIGN_CENTRE );
    item71->Add( item72, 5, wxALIGN_CENTER|wxALL, 5 );

    wxStaticText *item73 = new wxStaticText( parent, ID_TEXT, wxT("Output files:        "), wxDefaultPosition, wxDefaultSize, wxALIGN_CENTRE );
    item71->Add( item73, 5, wxALIGN_CENTER|wxALL, 5 );

    wxStaticText *item74 = new wxStaticText( parent, ID_TEXT, wxT("Output Freq:"), wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT );
    item71->Add( item74, 3, wxALIGN_CENTER|wxALL, 5 );

    item68->Add( item71, 0, wxGROW|wxALIGN_CENTER_VERTICAL, 0 );

    wxBoxSizer *item75 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item76 = new wxButton( parent, IDC_MM_EDIT_AMBER_INP, wxT("Input Pars:"), wxDefaultPosition, wxDefaultSize, 0 );
    item75->Add( item76, 1, wxALIGN_CENTER|wxALL, 0 );

    wxTextCtrl *item77 = new wxTextCtrl( parent, IDC_MM_INP_FILE, wxT(""), wxDefaultPosition, wxSize(140,-1), 0 );
    item75->Add( item77, 3, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item78 = new wxTextCtrl( parent, IDC_MM_LOG_FILE, wxT(""), wxDefaultPosition, wxSize(140,-1), 0 );
    item75->Add( item78, 3, wxALIGN_CENTER|wxALL, 5 );

    wxButton *item79 = new wxButton( parent, IDC_MM_LOAD_LOG_FILE, wxT("Log File"), wxDefaultPosition, wxDefaultSize, 0 );
    item75->Add( item79, 1, wxALIGN_CENTER|wxALL, 0 );

    wxTextCtrl *item80 = new wxTextCtrl( parent, IDC_MM_PRINT_FREQ, wxT(""), wxDefaultPosition, wxDefaultSize, 0 );
    item75->Add( item80, 0, wxALIGN_CENTER|wxALL, 5 );

    item68->Add( item75, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item81 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item82 = new wxButton( parent, IDC_MM_EDIT_AMBER_TOP, wxT("Top file:"), wxDefaultPosition, wxDefaultSize, 0 );
    item81->Add( item82, 1, wxALIGN_CENTER|wxALL, 0 );

    wxTextCtrl *item83 = new wxTextCtrl( parent, IDC_MM_TOP_FILE, wxT(""), wxDefaultPosition, wxSize(140,-1), 0 );
    item81->Add( item83, 3, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item84 = new wxTextCtrl( parent, IDC_MM_CRD_TRAJ_FILE, wxT(""), wxDefaultPosition, wxSize(140,-1), 0 );
    item81->Add( item84, 3, wxALIGN_CENTER|wxALL, 5 );

    wxButton *item85 = new wxButton( parent, IDC_MM_CHOOSE_MDCRD_FILE, wxT("Coordinates"), wxDefaultPosition, wxDefaultSize, 0 );
    item81->Add( item85, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item86 = new wxTextCtrl( parent, IDC_MM_WRT_COORD_FREQ, wxT(""), wxDefaultPosition, wxDefaultSize, 0 );
    item81->Add( item86, 0, wxALIGN_CENTER|wxALL, 5 );

    item68->Add( item81, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item87 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item88 = new wxStaticText( parent, ID_TEXT, wxT("Init crd file:     "), wxDefaultPosition, wxDefaultSize, 0 );
    item87->Add( item88, 1, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item89 = new wxTextCtrl( parent, IDC_MM_INIT_CRD_FILE, wxT(""), wxDefaultPosition, wxSize(140,-1), 0 );
    item87->Add( item89, 3, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item90 = new wxTextCtrl( parent, IDC_MM_VEL_TRAJ_FILE, wxT(""), wxDefaultPosition, wxSize(140,-1), 0 );
    item87->Add( item90, 3, wxALIGN_CENTER|wxALL, 5 );

    wxButton *item91 = new wxButton( parent, IDC_MM_CHOOSE_MDVEL_FILE, wxT("Velocities"), wxDefaultPosition, wxDefaultSize, 0 );
    item87->Add( item91, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item92 = new wxTextCtrl( parent, IDC_MM_WRT_VEL_FREQ, wxT(""), wxDefaultPosition, wxDefaultSize, 0 );
    item87->Add( item92, 0, wxALIGN_CENTER|wxALL, 5 );

    item68->Add( item87, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item93 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item94 = new wxStaticText( parent, ID_TEXT, wxT("Restraints file:      "), wxDefaultPosition, wxDefaultSize, 0 );
    item93->Add( item94, 1, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item95 = new wxTextCtrl( parent, IDC_MM_CONSTR_CRD_FILE, wxT(""), wxDefaultPosition, wxSize(140,-1), 0 );
    item93->Add( item95, 3, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item96 = new wxTextCtrl( parent, IDC_MM_ENE_TRAJ_FILE, wxT(""), wxDefaultPosition, wxSize(140,-1), 0 );
    item93->Add( item96, 3, wxALIGN_CENTER|wxALL, 5 );

    wxButton *item97 = new wxButton( parent, IDC_MM_CHOOSE_MDENE_FILE, wxT("Energy"), wxDefaultPosition, wxDefaultSize, 0 );
    item93->Add( item97, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item98 = new wxTextCtrl( parent, IDC_MM_WRT_ENER_FREQ, wxT(""), wxDefaultPosition, wxDefaultSize, 0 );
    item93->Add( item98, 0, wxALIGN_CENTER|wxALL, 5 );

    item68->Add( item93, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item99 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item100 = new wxButton( parent, IDC_MM_EDIT_AMBER_RUN, wxT("UNIX Run file"), wxDefaultPosition, wxDefaultSize, 0 );
    item99->Add( item100, 1, wxALIGN_CENTER|wxALL, 0 );

    wxTextCtrl *item101 = new wxTextCtrl( parent, IDC_MM_AMBER_RUN_FILE, wxT(""), wxDefaultPosition, wxSize(140,-1), 0 );
    item99->Add( item101, 3, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item102 = new wxTextCtrl( parent, IDC_MM_RESTART_FILE, wxT(""), wxDefaultPosition, wxSize(140,-1), 0 );
    item99->Add( item102, 3, wxALIGN_CENTER|wxALL, 5 );

    wxButton *item103 = new wxButton( parent, IDC_MM_EDIT_AMBER_RST, wxT("Restart File"), wxDefaultPosition, wxDefaultSize, 0 );
    item99->Add( item103, 1, wxALIGN_CENTER|wxALL, 0 );

    wxTextCtrl *item104 = new wxTextCtrl( parent, IDC_MM_WRT_RSTRT_FREQ, wxT(""), wxDefaultPosition, wxDefaultSize, 0 );
    item99->Add( item104, 0, wxALIGN_CENTER|wxALL, 5 );

    item68->Add( item99, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item105 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item106 = new wxButton( parent, IDC_AMBER_LOAD_RESTART, wxT("Load Restart File"), wxDefaultPosition, wxDefaultSize, 0 );
    item105->Add( item106, 0, wxALIGN_CENTER|wxALL, 5 );

    wxButton *item107 = new wxButton( parent, IDC_AMBER_SAVE_RESTART, wxT("Save Restart File"), wxDefaultPosition, wxDefaultSize, 0 );
    item105->Add( item107, 0, wxALIGN_CENTER|wxALL, 5 );

    wxCheckBox *item108 = new wxCheckBox( parent, IDC_AMBER_WRITE_BINARY, wxT("Binary output"), wxDefaultPosition, wxDefaultSize, 0 );
    item105->Add( item108, 0, wxALIGN_CENTER|wxALL, 5 );

    item68->Add( item105, 0, wxALIGN_RIGHT|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item109 = new wxBoxSizer( wxHORIZONTAL );

    wxTextCtrl *item110 = new wxTextCtrl( parent, IDC_MM_CONSTR_TRAJ_FILE, wxT(""), wxDefaultPosition, wxSize(210,-1), 0 );
    item109->Add( item110, 1, wxALIGN_CENTER|wxALL, 5 );

    wxButton *item111 = new wxButton( parent, IDC_MM_CHOOSE_CONSTR_TRAJ_FILE, wxT("Constraints Traj"), wxDefaultPosition, wxDefaultSize, 0 );
    item109->Add( item111, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item112 = new wxTextCtrl( parent, IDC_MM_WRT_CONSTR_FREQ, wxT(""), wxDefaultPosition, wxSize(40,-1), 0 );
    item109->Add( item112, 0, wxALIGN_CENTER|wxALL, 5 );

    item68->Add( item109, 0, wxALIGN_RIGHT|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    item67->Add( item68, 0, wxALIGN_CENTER|wxALL, 5 );

    item0->Add( item67, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    if (set_sizer)
    {
        parent->SetSizer( item0 );
        if (call_fit)
            item0->SetSizeHints( parent );
    }

    return item0;
}

wxSizer *mm_md_anal_page( wxWindow *parent, bool call_fit, bool set_sizer )
{
    wxBoxSizer *item0 = new wxBoxSizer( wxHORIZONTAL );

    wxBoxSizer *item1 = new wxBoxSizer( wxVERTICAL );

    wxBoxSizer *item2 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item3 = new wxButton( parent, IDC_MM_INDEX_TRAJ, wxT("Index Trajectory"), wxDefaultPosition, wxDefaultSize, 0 );
    item2->Add( item3, 0, wxALIGN_CENTER|wxALL, 5 );

    wxStaticText *item4 = new wxStaticText( parent, ID_TEXT, wxT("N Pt:"), wxDefaultPosition, wxDefaultSize, 0 );
    item2->Add( item4, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item5 = new wxTextCtrl( parent, IDC_MM_NPT_TRAJ, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item2->Add( item5, 0, wxALIGN_CENTER|wxALL, 5 );

    item1->Add( item2, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item6 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item7 = new wxStaticText( parent, ID_TEXT, wxT("Curr Pt:"), wxDefaultPosition, wxDefaultSize, 0 );
    item6->Add( item7, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item8 = new wxTextCtrl( parent, IDC_MM_CURR_PT, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item6->Add( item8, 0, wxALIGN_CENTER|wxALL, 5 );

    wxButton *item9 = new wxButton( parent, IDC_MM_SET_CURR_PT, wxT("Set Curr Pt"), wxDefaultPosition, wxDefaultSize, 0 );
    item6->Add( item9, 0, wxALIGN_CENTER|wxALL, 5 );

    item1->Add( item6, 0, wxALIGN_CENTER|wxALL, 5 );

    wxSlider *item10 = new wxSlider( parent, ID_MM_SLIDER_MD, 0, 0, 100, wxDefaultPosition, wxSize(100,-1), wxSL_HORIZONTAL|wxSL_LABELS );
    item10->Enable( false );
    item1->Add( item10, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxButton *item11 = new wxButton( parent, IDC_MM_PLAYBACK_TRJ, wxT("Playback Trajectory"), wxDefaultPosition, wxDefaultSize, 0 );
    item1->Add( item11, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxButton *item12 = new wxButton( parent, IDC_MM_STOP, wxT("Stop Calculations"), wxDefaultPosition, wxDefaultSize, 0 );
    item1->Add( item12, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxStaticText *item13 = new wxStaticText( parent, ID_TEXT, wxT("Script to run on MD points:"), wxDefaultPosition, wxDefaultSize, 0 );
    item1->Add( item13, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxTextCtrl *item14 = new wxTextCtrl( parent, IDC_MM_TRAJ_SCRIPT, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item1->Add( item14, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item15 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item16 = new wxButton( parent, IDC_MM_CHOOSE_MDANAL_SCRIPT, wxT("Choose Script"), wxDefaultPosition, wxDefaultSize, 0 );
    item15->Add( item16, 1, wxALIGN_CENTER|wxALL, 5 );

    wxButton *item17 = new wxButton( parent, IDC_MM_EDIT_MDANAL_SCRIPT, wxT("Edit Script"), wxDefaultPosition, wxDefaultSize, 0 );
    item15->Add( item17, 1, wxALIGN_CENTER|wxALL, 5 );

    item1->Add( item15, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxStaticText *item18 = new wxStaticText( parent, ID_TEXT, wxT("MD playback control:"), wxDefaultPosition, wxDefaultSize, 0 );
    item1->Add( item18, 0, wxALIGN_CENTER|wxALL, 5 );

    wxBoxSizer *item19 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item20 = new wxStaticText( parent, ID_TEXT, wxT("Delay time (s)"), wxDefaultPosition, wxDefaultSize, 0 );
    item19->Add( item20, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item21 = new wxTextCtrl( parent, IDC_MM_DELAY, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item19->Add( item21, 1, wxALIGN_CENTER|wxALL, 5 );

    item1->Add( item19, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item22 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item23 = new wxStaticText( parent, ID_TEXT, wxT("Minimal update view time(s)"), wxDefaultPosition, wxDefaultSize, 0 );
    item22->Add( item23, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item24 = new wxTextCtrl( parent, IDC_MM_UPDATE_VIEW_INTERVAL, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item22->Add( item24, 1, wxALIGN_CENTER|wxALL, 5 );

    item1->Add( item22, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item25 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item26 = new wxStaticText( parent, ID_TEXT, wxT("Index of the First Point to analyze:"), wxDefaultPosition, wxDefaultSize, 0 );
    item25->Add( item26, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item27 = new wxTextCtrl( parent, IDC_MM_SKIP_INIT, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item25->Add( item27, 1, wxALIGN_CENTER|wxALL, 5 );

    item1->Add( item25, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item28 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item29 = new wxStaticText( parent, ID_TEXT, wxT(" Num Points to Step during analysis:"), wxDefaultPosition, wxDefaultSize, 0 );
    item28->Add( item29, 0, wxALIGN_CENTER|wxALL, 0 );

    wxTextCtrl *item30 = new wxTextCtrl( parent, IDC_MM_SKIP_BETWEEN, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item28->Add( item30, 1, wxALIGN_CENTER|wxALL, 5 );

    item1->Add( item28, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item31 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item32 = new wxStaticText( parent, ID_TEXT, wxT("Idx of the last point to analyze:"), wxDefaultPosition, wxDefaultSize, 0 );
    item31->Add( item32, 0, wxALIGN_CENTER|wxALL, 0 );

    wxTextCtrl *item33 = new wxTextCtrl( parent, IDC_MM_LAST_PT_IDX, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item31->Add( item33, 1, wxALIGN_CENTER|wxALL, 5 );

    item1->Add( item31, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxCheckBox *item34 = new wxCheckBox( parent, IDC_MM_UPDATE_MOL_VIEW, wxT("Update Molecular View "), wxDefaultPosition, wxDefaultSize, 0 );
    item1->Add( item34, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxCheckBox *item35 = new wxCheckBox( parent, IDC_MM_ANAL_RUN_IN_THREAD, wxT("Run in Separate Thread"), wxDefaultPosition, wxDefaultSize, 0 );
    item1->Add( item35, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    item0->Add( item1, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5 );

    wxBoxSizer *item36 = new wxBoxSizer( wxVERTICAL );

    wxBoxSizer *item37 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item38 = new wxStaticText( parent, ID_TEXT, wxT("MD Trajectory files:"), wxDefaultPosition, wxDefaultSize, 0 );
    item38->SetFont( wxFont( 12, wxROMAN, wxNORMAL, wxBOLD ) );
    item37->Add( item38, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    item36->Add( item37, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item39 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item40 = new wxButton( parent, IDC_MM_CHOOSE_MDCRD_FILE, wxT("Coordinates:"), wxDefaultPosition, wxDefaultSize, 0 );
    item39->Add( item40, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item41 = new wxTextCtrl( parent, IDC_MM_CRD_TRAJ_FILE_2, wxT(""), wxDefaultPosition, wxSize(250,-1), wxTE_READONLY );
    item39->Add( item41, 3, wxALIGN_CENTER|wxALL, 5 );

    item36->Add( item39, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item42 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item43 = new wxButton( parent, IDC_MM_CHOOSE_MDVEL_FILE, wxT("Velocities:"), wxDefaultPosition, wxDefaultSize, 0 );
    item42->Add( item43, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item44 = new wxTextCtrl( parent, IDC_MM_VEL_TRAJ_FILE_2, wxT(""), wxDefaultPosition, wxSize(140,-1), wxTE_READONLY );
    item42->Add( item44, 3, wxALIGN_CENTER|wxALL, 5 );

    item36->Add( item42, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item45 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item46 = new wxButton( parent, IDC_MM_CHOOSE_MDENE_FILE, wxT("Energy:"), wxDefaultPosition, wxDefaultSize, 0 );
    item45->Add( item46, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item47 = new wxTextCtrl( parent, IDC_MM_ENE_TRAJ_FILE_2, wxT(""), wxDefaultPosition, wxSize(140,-1), wxTE_READONLY );
    item45->Add( item47, 3, wxALIGN_CENTER|wxALL, 5 );

    item36->Add( item45, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item48 = new wxBoxSizer( wxHORIZONTAL );

    wxCheckBox *item49 = new wxCheckBox( parent, IDC_MM_CHK_RMSD_ANAL, wxT("Atom RMSD Analysis:"), wxDefaultPosition, wxSize(250,20), 0 );
    item49->SetFont( wxFont( 11, wxROMAN, wxNORMAL, wxBOLD ) );
    item48->Add( item49, 0, wxALIGN_CENTER|wxALL, 5 );

    item36->Add( item48, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item50 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item51 = new wxButton( parent, IDC_MM_CHOOSE_FIT_ATOMS, wxT("Atoms to Fit:"), wxDefaultPosition, wxDefaultSize, 0 );
    item50->Add( item51, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item52 = new wxTextCtrl( parent, IDC_MM_FIT_ATOMS, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item50->Add( item52, 1, wxALIGN_CENTER|wxALL, 5 );

    item36->Add( item50, 1, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item53 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item54 = new wxStaticText( parent, IDC_MM_REF_CRD_FIT_FILE_LBL, wxT("Ref Crd (Fit):"), wxDefaultPosition, wxDefaultSize, 0 );
    item53->Add( item54, 0, wxALIGN_CENTER|wxALL, 10 );

    wxArrayString strs55;
    strs55.Add(wxT("Current Coordinates"));
    strs55.Add(wxT("First Trajectory Point"));
    strs55.Add(wxT("XYZ Coordinates File"));
    wxComboBox *item55 = new wxComboBox( parent, IDC_MM_REF_CRD_FIT_TYPE, wxT(""), wxDefaultPosition, wxDefaultSize, strs55, wxCB_DROPDOWN );
    item53->Add( item55, 1, wxALIGN_CENTER|wxALL, 5 );

    item36->Add( item53, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item56 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item57 = new wxButton( parent, IDC_MM_CHOOSE_REF_CRD_FIT, wxT("Ref Crd (Fit) File:"), wxDefaultPosition, wxDefaultSize, 0 );
    item56->Add( item57, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item58 = new wxTextCtrl( parent, IDC_MM_REF_CRD_FIT_FILE, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item56->Add( item58, 1, wxALIGN_CENTER|wxALL, 5 );

    item36->Add( item56, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item59 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item60 = new wxButton( parent, IDC_MM_CHOOSE_RMSD_ATOMS, wxT("Atoms To Calc RMSD:"), wxDefaultPosition, wxDefaultSize, 0 );
    item59->Add( item60, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item61 = new wxTextCtrl( parent, IDC_MM_RMSD_ATOMS, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item59->Add( item61, 1, wxALIGN_CENTER|wxALL, 5 );

    item36->Add( item59, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item62 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText *item63 = new wxStaticText( parent, IDC_MM_REF_CRD_RMSD_FILE_LBL, wxT("Ref Crd (RMSD):"), wxDefaultPosition, wxDefaultSize, 0 );
    item62->Add( item63, 0, wxALIGN_CENTER|wxALL, 5 );

    wxString *strs64 = (wxString*) NULL;
    wxComboBox *item64 = new wxComboBox( parent, IDC_MM_REF_CRD_RMSD_TYPE, wxT(""), wxDefaultPosition, wxSize(100,-1), 0, strs64, wxCB_DROPDOWN );
    item62->Add( item64, 1, wxALIGN_CENTER|wxALL, 5 );

    item36->Add( item62, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item65 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item66 = new wxButton( parent, IDC_MM_CHOOSE_REF_CRD_RMSD, wxT("Ref Crd (RMSD) File:"), wxDefaultPosition, wxDefaultSize, 0 );
    item65->Add( item66, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item67 = new wxTextCtrl( parent, IDC_MM_REF_CRD_RMSD_FILE, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item65->Add( item67, 1, wxALIGN_CENTER|wxALL, 5 );

    item36->Add( item65, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item68 = new wxBoxSizer( wxHORIZONTAL );

    wxButton *item69 = new wxButton( parent, IDC_MM_CHOOSE_RMSD_FILE, wxT("File to Save RMSD:"), wxDefaultPosition, wxDefaultSize, 0 );
    item68->Add( item69, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item70 = new wxTextCtrl( parent, IDC_MM_RMSD_FILE_NAME, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item68->Add( item70, 1, wxALIGN_CENTER|wxALL, 5 );

    item36->Add( item68, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    wxBoxSizer *item71 = new wxBoxSizer( wxHORIZONTAL );

    wxCheckBox *item72 = new wxCheckBox( parent, IDC_MM_CHK_RMSD_ATOM, wxT("RMSD Per Atom:"), wxDefaultPosition, wxDefaultSize, 0 );
    item71->Add( item72, 0, wxALIGN_CENTER|wxALL, 5 );

    wxButton *item73 = new wxButton( parent, IDC_MM_CHOOSE_RMSD_ATOM_FILE, wxT("File:"), wxDefaultPosition, wxSize(40,-1), 0 );
    item71->Add( item73, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item74 = new wxTextCtrl( parent, IDC_MM_RMSD_ATOM_FILE, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item71->Add( item74, 1, wxALIGN_CENTER|wxALL, 5 );

    item36->Add( item71, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item75 = new wxBoxSizer( wxHORIZONTAL );

    wxCheckBox *item76 = new wxCheckBox( parent, IDC_MM_CHK_AVG_COORD, wxT("Average Coords:"), wxDefaultPosition, wxDefaultSize, 0 );
    item75->Add( item76, 0, wxALIGN_CENTER|wxALL, 5 );

    wxButton *item77 = new wxButton( parent, IDC_MM_CHOOSE_AVG_COORD_FILE, wxT("File:"), wxDefaultPosition, wxSize(40,-1), 0 );
    item75->Add( item77, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item78 = new wxTextCtrl( parent, IDC_MM_AVG_COORD_FILE, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item75->Add( item78, 1, wxALIGN_CENTER|wxALL, 5 );

    item36->Add( item75, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    wxBoxSizer *item79 = new wxBoxSizer( wxHORIZONTAL );

    wxCheckBox *item80 = new wxCheckBox( parent, IDC_MM_CHK_RMSF_ATOM, wxT("RMSF Per Atom:"), wxDefaultPosition, wxDefaultSize, 0 );
    item79->Add( item80, 0, wxALIGN_CENTER|wxALL, 5 );

    wxButton *item81 = new wxButton( parent, IDC_MM_CHOOSE_RMSF_ATOM_FILE, wxT("File:"), wxDefaultPosition, wxSize(40,-1), 0 );
    item79->Add( item81, 0, wxALIGN_CENTER|wxALL, 5 );

    wxTextCtrl *item82 = new wxTextCtrl( parent, IDC_MM_RMSF_ATOM_FILE, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
    item79->Add( item82, 1, wxALIGN_CENTER|wxALL, 5 );

    item36->Add( item79, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 0 );

    item0->Add( item36, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5 );

    if (set_sizer)
    {
        parent->SetSizer( item0 );
        if (call_fit)
            item0->SetSizeHints( parent );
    }

    return item0;
}


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

	// ptr_mm_mod->p_mm_dlg = this;

	// Event wiring (migrated from BEGIN_EVENT_TABLE / END_EVENT_TABLE).
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnInitMolMech, this, IDC_FFPAR_INIT_MOLMECH);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnRestrainHBonds, this, IDC_RESTRAIN_HBONDS);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::ChooseRestrainedAtoms, this, IDC_MM_CHOOSE_RESTRAINED_ATOMS);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::ChooseMovingAtoms, this, IDC_MM_CHOOSE_MOVING_ATOMS);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::ChooseRestrRefCrd, this, IDC_MM_CHOOSE_RESTR_REF_CRD);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::LoadAtomAtomRestrFile, this, IDC_MM_LOAD_HARM_CONSTR_FILE);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnSaveExtProgInp, this, IDC_SAVE_EXT_PROG_INP);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnRunMMCalc, this, IDC_MM_RUN_CALC);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnAmberLoadRestart, this, IDC_AMBER_LOAD_RESTART);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnCalcSinglePtEne, this, IDC_MM_CALC_SP_ENE);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnLoadLogFile, this, IDC_MM_LOAD_LOG_FILE);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnChooseMDCrdFile, this, IDC_MM_CHOOSE_MDCRD_FILE);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnChooseMDVelFile, this, IDC_MM_CHOOSE_MDVEL_FILE);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnChooseMDEneFile, this, IDC_MM_CHOOSE_MDENE_FILE);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnChooseConstrTrajFile, this, IDC_MM_CHOOSE_CONSTR_TRAJ_FILE);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnMMStop, this, IDC_MM_STOP);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnPlaybackTrj, this, IDC_MM_PLAYBACK_TRJ);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnIndexTrj, this, IDC_MM_INDEX_TRAJ);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnSetCurrPt, this, IDC_MM_SET_CURR_PT);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnChooseMDAnalScript, this, IDC_MM_CHOOSE_MDANAL_SCRIPT);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnEditMDAnalScript, this, IDC_MM_EDIT_MDANAL_SCRIPT);
	Bind(wxEVT_CHECKBOX, &MolMechDlgWX::OnChkRMSDAnal, this, IDC_MM_CHK_RMSD_ANAL);
	Bind(wxEVT_CHECKBOX, &MolMechDlgWX::OnChkMMRunInt, this, IDC_MM_RUN_INT);
	Bind(wxEVT_CHOICE, &MolMechDlgWX::OnChoiceRunType, this, IDC_MM_RUN_TYPE);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::ChooseFitAtoms, this, IDC_MM_CHOOSE_FIT_ATOMS);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::ChooseRMSDAtoms, this, IDC_MM_CHOOSE_RMSD_ATOMS);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::ChooseRefCrdFit, this, IDC_MM_CHOOSE_REF_CRD_FIT);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::ChooseRefCrdRMSD, this, IDC_MM_CHOOSE_REF_CRD_RMSD);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnEditAmberInp, this, IDC_MM_EDIT_AMBER_INP);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnEditAmberTop, this, IDC_MM_EDIT_AMBER_TOP);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnEditAmberRun, this, IDC_MM_EDIT_AMBER_RUN);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnEditAmberRst, this, IDC_MM_EDIT_AMBER_RST);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnAmberSaveRestart, this, IDC_AMBER_SAVE_RESTART);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnShowMMInfo, this, IDC_MM_SHOW_INFO);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnEditPeriodicBox, this, IDC_MM_EDIT_PERIODIC_BOX);
	Bind(wxEVT_COMBOBOX, &MolMechDlgWX::OnSelChangeFFTypeDefault, this, IDC_MM_FF_TYPE_DEFAULT);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnUpdateElemList, this, IDC_MM_UPDATE_ELEM_LIST);
	Bind(wxEVT_RADIOBOX, &MolMechDlgWX::OnUpdateElemList, this, IDC_RADIO_ELEMENTS);
	Bind(wxEVT_LISTBOX, &MolMechDlgWX::OnChangeSelElem, this, IDC_MM_ELEM_LIST);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnSetNewPar, this, IDC_RESPAR_SET_NEW_PAR);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnDelImprAng, this, IDC_DEL_IMPROPER_ANG);
	Bind(wxEVT_BUTTON, &MolMechDlgWX::OnSetDnaCSPars, this, IDC_SET_DNA_CS_PARS);
	Bind(wxEVT_MENU, &MolMechDlgWX::OnSaveResffFromMort, this, IDC_SAVE_RESFF_FROM_MORT);
	Bind(wxEVT_MENU, &MolMechDlgWX::OnSaveAtomIndRestrArb, this, IDM_SAVE_ATOM_IND_RESTR_ARB);
	Bind(wxEVT_NOTEBOOK_PAGE_CHANGING, &MolMechDlgWX::OnChangingPage, this, ID_MM_DLG);
	Bind(wxEVT_NOTEBOOK_PAGE_CHANGED, &MolMechDlgWX::OnChangePage, this, ID_MM_DLG);
	Bind(wxEVT_CLOSE_WINDOW, &MolMechDlgWX::OnClose, this);

	OnInitDialog();
}

MolMechDlgWX::~MolMechDlgWX()
{
     dlg_open = FALSE;
	 delete p_mm_info_dlg;
}


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
	TransferExtProgDataToWindow();
	TransferMMParamToWindow();
	TransferRunTypeDataToWindow();
	TransferAtomSuperimposeDataToWindow();
	return wxFrame::TransferDataToWindow();
}

bool MolMechDlgWX::TransferMMParamToWindow()
{
	wxTextCtrl* edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_MOVING_ATOMS);
	if (edit_ctrl == NULL) return false;
	edit_ctrl->SetValue(ptr_mm_mod->p_mm_model->moving_atoms.c_str());

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_RESTRAINED_ATOMS);
	edit_ctrl->SetValue(ptr_mm_mod->p_mm_model->restrained_atoms.c_str());

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_FF_TYPE);
	edit_ctrl->SetValue(ptr_mm_mod->p_mm_model->ff_type.label());

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_NONB_CUTOFF_DIST);
	edit_ctrl->SetValue(wxString::Format("%14.9f", ptr_mm_mod->p_mm_model->GetNBCutDist()));

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_SCALE_14_ELECTR);
	edit_ctrl->SetValue(wxString::Format("%6.2f", ptr_mm_mod->p_mm_model->GetScale14Electr()));

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_SCALE_14_VDW);
	edit_ctrl->SetValue(wxString::Format("%6.2f", ptr_mm_mod->p_mm_model->GetScale14VdW()));

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_DIEL_CONST);
	edit_ctrl->SetValue(wxString::Format("%6.2f", ptr_mm_mod->p_mm_model->GetDielConst()));

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_ION_STRENGTH);
	edit_ctrl->SetValue(wxString::Format("%6.2f", ptr_mm_mod->p_mm_model->GetIonStrength()));

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_THOLE_DUMP_CONST);
	edit_ctrl->SetValue(wxString::Format("%6.4f", ptr_mm_mod->p_mm_model->GetTholeExponCoef()));

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_RESTR_FRC_CONST);
	edit_ctrl->SetValue(wxString::Format("%8.3f", ptr_mm_mod->p_mm_model->GetAtomRestrForceConst()));

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_NUM_HARM_CONSTR);
	edit_ctrl->SetValue(wxString::Format("%d", ptr_mm_mod->p_mm_model->GetNumHarmConstr()));

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MM_NPT_TRAJ);
	edit_ctrl->SetValue(wxString::Format("%d", ptr_mm_mod->p_traj_anal_mod->pt_pos.size()));

	return true;
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
		std::string cmd_line = pApp->word_editor.string();
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
		if (pval)
		{
			pval->UpdateItems();
			pval->TransferToWindow();
		}
	}
	if( choice_ctrl_electr_method )
	{
		HaEnumValidator* pval = dynamic_cast<HaEnumValidator*>(choice_ctrl_electr_method->GetValidator());
		if (pval)
		{
			pval->UpdateItems();
			pval->TransferToWindow();
		}
	}
	//TransferDataToWindow();
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
			ptr_mm_mod->p_gromacs_driver->SetCompatibleParams();
			this->OnChangePeriodicity();
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

	std::string cmd_line = pApp->word_editor.string();
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
		"mdcrd","AMBER Trajectory (*.mdcrd)|*.mdcrd|NetCDF Trajectory (*.nc)|*.nc|All files (*.*)|*.*");

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
		std::string cmd_line = pApp->word_editor.string();
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
		std::string cmd_line = pApp->word_editor.string();
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
		std::string cmd_line = pApp->word_editor.string();
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
		std::string cmd_line = pApp->word_editor.string();
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
		        
		  for( auto spb : ptr_mm_mod->p_mm_model->MBonds )
		  {
			  MMBond* pbond = spb.get();
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
		  for( auto spa : ptr_mm_mod->p_mm_model->ValAngles )
		  {
			  MMValAngle* pval = spa.get();
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

#ifdef HARLEM_MPI
	pApp->mpi_driver->SendXmlMsgAllProc(msg.c_str());
#endif
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
