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


#include "abstree.h"
#include "haatom.h"
#include "hamolecule.h"
#include "hamolset.h"
#include "hamolview.h"
#include "ha_wx_res_wdr.h"
#include "edit_mut_map_dlg_wx.h"
#include "dialogs_wx_1.h"


EditMutMapDlg::EditMutMapDlg(MolSet* pmset_par, wxWindow* parent): wxFrame(parent, -1, "Edit Mutation Map")
{
    mol_shift = 5.0;
    pmset = pmset_par;
    pat1 = NULL;
    pat2 = NULL;

    wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
    SetBackgroundColour(back_colour);

    edit_mut_map_dlg(this, TRUE);
    OnInitDialog();

    dlg_open = TRUE;
    active_dlg_ptr = this;
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
    color_mapped_atoms_chk = (wxCheckBox*)FindWindow(IDC_COLOR_MAPPED_ATOMS);

    mol_shift_txt = (wxTextCtrl*)FindWindow(IDC_MOL_SHIFT);
    HaMolView* pView = pmset->GetActiveMolView();
    if (pView != NULL)
    {
        pView->UseMolShift = TRUE;
        mol_shift_txt->SetValidator(wxDoubleValidator(&(pView->mol_shift), "%6.2f"));
    }
    atom_pair_list = (wxListBox*)FindWindow(IDC_LIST_ATOM_PAIRS);

    atom1_edt = (wxTextCtrl*)FindWindow(IDC_EDIT_ATOM_1);
    atom2_edt = (wxTextCtrl*)FindWindow(IDC_EDIT_ATOM_2);
    
    active_atom_edt = atom1_edt;

    TransferDataToWindow();
    pmset->RefreshAllViews(RFApply);
}

//----------------------------------------------------------------------------
// EditMutMapDlg
//----------------------------------------------------------------------------

int EditMutMapDlg::dlg_open = FALSE;
EditMutMapDlg* EditMutMapDlg::active_dlg_ptr = NULL;

void EditMutMapDlg::OnLoadMutMap(wxCommandEvent& event) noexcept
{
    if (pmset->GetNMol() != 2)
    {
        PrintLog("The Number of Molecules in MolSet = %5d  is not equal 2 \n", pmset->GetNMol());
        return;
    }

    wxString mut_map_fname_wx = ::wxFileSelector("Choose Mutation Map File",
        ::wxGetCwd(), "map_mutation.xml",
        "xml", "*.xml");

    std::string mut_map_fname = mut_map_fname_wx.ToStdString();

    if (!mut_map_fname.empty())
    {
        p_mut_map.reset( new MutationMap(pmset) );
        p_mut_map->LoadArbalestMutMap(mut_map_fname);
    }
    TransferDataToWindow();
}

void EditMutMapDlg::OnSaveMutMap(wxCommandEvent& event) noexcept
{

}

void EditMutMapDlg::OnSelectAtomPair(wxListEvent& event) noexcept
{

}

void EditMutMapDlg::OnAddAtomPair(wxCommandEvent& event) noexcept
{
    if (pat1 != NULL && pat2 != NULL)
    {
        p_mut_map->atom_atom_map[pat1] = pat2;
    }
    TransferDataToWindow();
}

void EditMutMapDlg::OnDeleteAtomPair(wxCommandEvent& event) noexcept
{
    int n = atom_pair_list->GetCount();
    for (int i = 0; i < n; i++)
    {
        if (!atom_pair_list->IsSelected(i)) continue;
        std::pair<HaAtom*, HaAtom*>* p_atom_pair = (std::pair<HaAtom*, HaAtom*>*) atom_pair_list->GetClientData(i);
        if (p_atom_pair == NULL) continue;
        HaAtom* aptr1 = (*p_atom_pair).first;
        HaAtom* aptr2 = (*p_atom_pair).second;

        auto mitr = p_mut_map->atom_atom_map.begin();
        for (; mitr != p_mut_map->atom_atom_map.end(); )
        {
            if ((*mitr).first == aptr1 && (*mitr).second == aptr2)
            {
                mitr = p_mut_map->atom_atom_map.erase(mitr);
            }
            else
            {
                mitr++;
            }
        }
    }
    TransferDataToWindow();
}

void EditMutMapDlg::OnClose(wxCloseEvent& event)
{
    EditMutMapDlg::dlg_open = FALSE;
    EditMutMapDlg::active_dlg_ptr = NULL;
    event.Skip();
}

void EditMutMapDlg::OnSetMolShift(wxCommandEvent& event) noexcept
{
    TransferDataFromWindow();
    pmset->RefreshAllViews(RFApply);
}

void EditMutMapDlg::OnAtomSelect(wxCommandEvent& event) noexcept
{
    char buf[256], buf2[256];
    if (PkAtom)
    {
        std::string sref = PkAtom->GetRef(HaAtom::ATOMREF_NO_MOL);
        
        if (active_atom_edt == atom1_edt)
        {
            pat1 = PkAtom;
            atom1_edt->SetValue(sref);
            active_atom_edt = atom2_edt;
            atom2_edt->SetFocus();
        }
        else if( active_atom_edt == atom2_edt )
        {
            pat2 = PkAtom;
            atom2_edt->SetValue(sref);
            active_atom_edt = atom1_edt;
            atom1_edt->SetFocus(); 
        }
    }
}

void EditMutMapDlg::OnChkColorMappedAtoms(wxCommandEvent& event) noexcept
{
    TransferDataToWindow();
}

bool EditMutMapDlg::TransferDataFromWindow()
{
    return wxFrame::TransferDataFromWindow();
}


bool EditMutMapDlg::TransferDataToWindow()
{
    HaMolView* pView = pmset->GetActiveMolView();
    atom_pair_list->Clear();

    if (p_mut_map)
    {
        HaAtom* aptr;
        if (p_mut_map->pmol1 != NULL)
        {
            AtomIteratorMolecule aitr1(p_mut_map->pmol1);
            for (aptr = aitr1.GetFirstAtom(); aptr; aptr = aitr1.GetNextAtom())
            {
                pView->ColorAtomCPK(aptr);
            }
        }
        if (p_mut_map->pmol1 != NULL)
        {
            AtomIteratorMolecule aitr2(p_mut_map->pmol2);
            for (aptr = aitr2.GetFirstAtom(); aptr; aptr = aitr2.GetNextAtom())
            {
                pView->ColorAtomCPK(aptr);
            }
        }
        auto mitr = p_mut_map->atom_atom_map.begin();
        for (; mitr != p_mut_map->atom_atom_map.end(); mitr++)
        {
            HaAtom* aptr1 = (*mitr).first;
            HaAtom* aptr2 = (*mitr).second;
            std::string s1 = aptr1->GetRef(HaAtom::ATOMREF_NO_MOL);
            std::string s2 = aptr2->GetRef(HaAtom::ATOMREF_NO_MOL);
            wxString lbl = s1 + " -> " + s2;
            auto p_atom_pair = new std::pair<HaAtom*, HaAtom*>(aptr1, aptr2);
            atom_pair_list->Append( lbl, p_atom_pair );

            if (color_mapped_atoms_chk->IsChecked())
            {
                pView->ColorAtom(aptr1, "YELLOW");
                pView->ColorAtom(aptr2, "YELLOW");
            }
        }
    }
    pmset->RefreshAllViews(RFRefresh | RFColour);

    return wxFrame::TransferDataToWindow();
}





// WDR: event table for EditMutMapDlg 
BEGIN_EVENT_TABLE(EditMutMapDlg, wxFrame)
EVT_BUTTON(IDC_SET_MOL_SHIFT,     EditMutMapDlg::OnSetMolShift)
EVT_BUTTON(IDC_ADD_ATOM_PAIR,     EditMutMapDlg::OnAddAtomPair)
EVT_BUTTON(IDC_DELETE_ATOM_PAIR,  EditMutMapDlg::OnDeleteAtomPair)
EVT_BUTTON(IDC_LOAD_MUTATION_MAP, EditMutMapDlg::OnLoadMutMap)
EVT_BUTTON(IDC_SAVE_MUTATION_MAP, EditMutMapDlg::OnSaveMutMap)
EVT_CHECKBOX(IDC_COLOR_MAPPED_ATOMS, EditMutMapDlg::OnChkColorMappedAtoms)
EVT_LIST_ITEM_SELECTED(IDC_LIST_ATOM_PAIRS, EditMutMapDlg::OnSelectAtomPair)
EVT_BUTTON(IDU_ATOM_PICK, EditMutMapDlg::OnAtomSelect)
EVT_CLOSE(EditMutMapDlg::OnClose)
END_EVENT_TABLE()




