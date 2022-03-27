/*! \file edit_mut_map_dlg_wx.h

	wxWidgets Dialogs 1 in HARLEM

	\author Igor Kurnikov
	\date 2003-
*/
#pragma once

#if !defined(EDIT_MUT_MAP_DLG_H)
#define EDIT_MUT_MAP_DLG_H

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
	void OnSelectAtomPair(wxListEvent& event) noexcept;

	virtual bool TransferDataFromWindow();
	virtual bool TransferDataToWindow();

	static int dlg_open;

	wxCheckBox* view_atom_map_chk;       //!< View Atom Map CheckBox
	wxCheckBox* color_mappped_atoms_chk; //!< Color Mapped Atoms CheckBox

	wxTextCtrl* mol_shift_txt; //!< Text Control for molecular shift

	std::string file_name;
	std::string dir_name;

	wxButton* choose_file_btn;

	double mol_shift; //!< Molecular shift in Ang 

protected:

	std::shared_ptr<MutationMap> p_mut_map;

	MolSet* pmset; 
	void OnClose(wxCloseEvent& event);

private:
	DECLARE_EVENT_TABLE()
};

#endif // if !defined(EDIT_MUT_MAP_DLG_H)