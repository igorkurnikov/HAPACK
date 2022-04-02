/*! \file edit_mut_map_dlg_wx.h

	wxWidgets Dialogs 1 in HARLEM

	\author Igor Kurnikov
	\date 2003-
*/
#pragma once

#if !defined(EDIT_MUT_MAP_DLG_H)
#define EDIT_MUT_MAP_DLG_H

#if !defined MSC_VER_
#include <memory>
#endif

class MolSet;

class EditMutMapDlg : public wxFrame
{
public:
	// constructors and destructors
	EditMutMapDlg( MolSet* pmset, wxWindow* parent);
    virtual ~EditMutMapDlg();

	void OnInitDialog();

	void OnSetMolShift(wxCommandEvent& event) noexcept;
	void OnLoadMutMap(wxCommandEvent& event) noexcept;
	void OnSaveMutMap(wxCommandEvent& event) noexcept;
	void OnAddAtomPair(wxCommandEvent& event) noexcept;
	void OnDeleteAtomPair(wxCommandEvent& event) noexcept;
	void OnAtomSelect(wxCommandEvent& event) noexcept;
	void OnChkColorMappedAtoms(wxCommandEvent& event) noexcept;
	void OnSelectAtomPair(wxListEvent& event) noexcept;
	

	virtual bool TransferDataFromWindow();
	virtual bool TransferDataToWindow();

	static int dlg_open;
	static EditMutMapDlg* active_dlg_ptr;

	wxCheckBox* view_atom_map_chk;       //!< View Atom Map CheckBox
	wxCheckBox* color_mapped_atoms_chk; //!< Color Mapped Atoms CheckBox

	wxTextCtrl* mol_shift_txt; //!< Text Control for molecular shift
	wxListBox* atom_pair_list; //!< List of Atom Paits in the Mutation Map

	wxTextCtrl* atom1_edt;   //!< Text Cntrl for Atom1 
	wxTextCtrl* atom2_edt;   //!< Text Cntrl for Atom2

	wxTextCtrl* active_atom_edt; //!< Active Atom Text Cntrl

	HaAtom* pat1;  //!< Atom In Edit Text 1
	HaAtom* pat2;  //!< Atom in Edit Text 2

	double mol_shift; //!< Molecular shift in Ang 

protected:

	std::shared_ptr<MutationMap> p_mut_map;

	MolSet* pmset; 
	void OnClose(wxCloseEvent& event);

private:
	DECLARE_EVENT_TABLE()
};

#endif // if !defined(EDIT_MUT_MAP_DLG_H)
