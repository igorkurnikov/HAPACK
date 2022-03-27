/*! \file edit_mut_map_dlg_wx.cpp

    Edit Mutation Map Dialog declarations

    \author Igor Kurnikov
    \date 2022-
*/

#define EDIT_MUT_MAP_DLG_WX_CPP

#include "wx/wxprec.h"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include "ha_wx_aux_1.h"

#include "hamolset.h"
#include "hamolview.h"
#include "ha_wx_res_wdr.h"
#include "edit_mut_map_dlg_wx.h"


EditMutMapDlg::EditMutMapDlg(MolSet* pmset_par, wxWindow* parent): wxFrame(parent, -1, "Edit Mutation Map")
{
    mol_shift = 5.0;
    pmset = pmset_par;

    wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
    SetBackgroundColour(back_colour);

    edit_mut_map_dlg(this, TRUE);
    OnInitDialog();
    dlg_open = TRUE; 
}

EditMutMapDlg::~EditMutMapDlg()
{
    dlg_open = FALSE;
    HaMolView* pView = pmset->GetActiveMolView();
    if (pView != NULL)
    {
        pView->UseMolShift = FALSE;
    }
    pmset->RefreshAllViews(RFApply);
}


void EditMutMapDlg::OnInitDialog()
{
    view_atom_map_chk = (wxCheckBox*)FindWindow(IDC_VIEW_ATOM_MAP);
    color_mappped_atoms_chk = (wxCheckBox*)FindWindow(IDC_COLOR_MAPPED_ATOMS);

    mol_shift_txt = (wxTextCtrl*)FindWindow(IDC_MOL_SHIFT);
    HaMolView* pView = pmset->GetActiveMolView();
    if (pView != NULL)
    {
        pView->UseMolShift = TRUE;
        mol_shift_txt->SetValidator(wxDoubleValidator(&(pView->mol_shift), "%6.2f"));
    }
    TransferDataToWindow();
    pmset->RefreshAllViews(RFApply);
}

//----------------------------------------------------------------------------
// EditMutMapDlg
//----------------------------------------------------------------------------

int EditMutMapDlg::dlg_open = FALSE;

void EditMutMapDlg::OnLoadMutMap(wxCommandEvent& event) noexcept
{

}

void EditMutMapDlg::OnSaveMutMap(wxCommandEvent& event) noexcept
{

}

void EditMutMapDlg::OnSelectAtomPair(wxListEvent& event) noexcept
{

}

void EditMutMapDlg::OnClose(wxCloseEvent& event)
{
    EditMutMapDlg::dlg_open = FALSE;
    event.Skip();
}

void EditMutMapDlg::OnSetMolShift(wxCommandEvent& event) noexcept
{
    TransferDataFromWindow();
    pmset->RefreshAllViews(RFApply);
}

bool EditMutMapDlg::TransferDataFromWindow()
{
    return wxFrame::TransferDataFromWindow();
}

bool EditMutMapDlg::TransferDataToWindow()
{

    return wxFrame::TransferDataToWindow();
}




// WDR: event table for EditMutMapDlg 
BEGIN_EVENT_TABLE(EditMutMapDlg, wxFrame)
EVT_BUTTON(IDC_SET_MOL_SHIFT,     EditMutMapDlg::OnSetMolShift)
EVT_BUTTON(IDC_LOAD_MUTATION_MAP, EditMutMapDlg::OnLoadMutMap)
EVT_BUTTON(IDC_SAVE_MUTATION_MAP, EditMutMapDlg::OnSaveMutMap)
EVT_LIST_ITEM_SELECTED(IDC_LIST_ATOM_PAIRS, EditMutMapDlg::OnSelectAtomPair)
EVT_CLOSE(EditMutMapDlg::OnClose)
END_EVENT_TABLE()




