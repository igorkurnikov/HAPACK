/*! \file save_mol_file_dlg_base.h

    Base class for Save Molecular File dialog - generated layout.
    This class can be regenerated from save_mol_file_dlg.xrc via wxFormBuilder.

    \author Igor Kurnikov
    \date 2024-
*/

#if !defined(SAVE_MOL_FILE_DLG_BASE_H)
#define SAVE_MOL_FILE_DLG_BASE_H

#include "ha_file_dlg1.h"

class wxListCtrl;

///////////////////////////////////////////////////////////////////////////////
/// Class SaveMolFileDlgBase
///////////////////////////////////////////////////////////////////////////////

class SaveMolFileDlgBase : public HaFileDlg1
{
public:
    SaveMolFileDlgBase( wxWindow* parent, wxWindowID id = wxID_ANY,
        const wxString& title = wxT("Save Molecular File"),
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxDEFAULT_DIALOG_STYLE );

    ~SaveMolFileDlgBase();

protected:
    // Controls
    wxTextCtrl*   m_dirName;
    wxButton*     m_chooseDir;
    wxListCtrl*   m_fileList;
    wxTextCtrl*   m_fileName;
    wxButton*     m_saveFile;
    wxChoice*     m_fileType;
    wxCheckBox*   m_saveTransformed;
    wxCheckBox*   m_saveConnect;
    wxCheckBox*   m_saveFragments;
    wxCheckBox*   m_saveAmberPdb;
    wxCheckBox*   m_saveSepWatMol;

    // Virtual event handlers (override in derived class)
    virtual void OnSaveFileClicked( wxCommandEvent& event ) { event.Skip(); }
};

#endif // SAVE_MOL_FILE_DLG_BASE_H
