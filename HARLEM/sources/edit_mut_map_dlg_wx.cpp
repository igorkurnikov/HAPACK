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

#include "hamolset.h"
#include "ha_wx_res_wdr.h"
#include "edit_mut_map_dlg_wx.h"


EditMutMapDlg::EditMutMapDlg(MolSet* pmset_par, wxWindow* parent): wxFrame(parent, -1, "Edit Mutation Map")
{
    mol_shift = 5.0;
    pmset = pmset_par;

    edit_mut_map_dlg(this, TRUE);
}


void EditMutMapDlg::OnInitDialog()
{
    TransferDataToWindow();
}

//----------------------------------------------------------------------------
// EditMutMapDlg
//----------------------------------------------------------------------------

int EditMutMapDlg::dlg_open = False;

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


// WDR: event table for EditMutMapDlg 
BEGIN_EVENT_TABLE(EditMutMapDlg, wxFrame)
EVT_BUTTON(IDC_LOAD_MUTATION_MAP, EditMutMapDlg::OnLoadMutMap)
EVT_BUTTON(IDC_SAVE_MUTATION_MAP, EditMutMapDlg::OnSaveMutMap)
EVT_LIST_ITEM_SELECTED(IDC_LIST_ATOM_PAIRS, EditMutMapDlg::OnSelectAtomPair)
EVT_CLOSE(EditMutMapDlg::OnClose)
END_EVENT_TABLE()




