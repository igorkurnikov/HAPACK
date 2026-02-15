/*! \file save_mol_file_dlg_base.cpp

    Base class for Save Molecular File dialog - generated layout.
    This class can be regenerated from save_mol_file_dlg.xrc via wxFormBuilder.

    \author Igor Kurnikov
    \date 2024-
*/

#include "hampi.h"

#include "wx/wx.h"
#include <wx/listctrl.h>

#include "ha_wx_res_wdr.h"
#include "save_mol_file_dlg_base.h"

SaveMolFileDlgBase::SaveMolFileDlgBase( wxWindow* parent, wxWindowID id,
    const wxString& title, const wxPoint& pos, const wxSize& size, long style )
    : HaFileDlg1( parent, id, title, pos, size, style )
{
    wxBoxSizer* bSizer1 = new wxBoxSizer( wxVERTICAL );

    // Row 1: Directory
    wxBoxSizer* bSizer2 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText* staticText1 = new wxStaticText( this, ID_TEXT, wxT("Directory:"), wxDefaultPosition, wxDefaultSize, 0 );
    bSizer2->Add( staticText1, 0, wxALIGN_CENTER|wxALL, 5 );

    m_dirName = new wxTextCtrl( this, IDC_DIR_NAME, wxT(""), wxDefaultPosition, wxSize(250,-1), 0 );
    bSizer2->Add( m_dirName, 0, wxALIGN_CENTER|wxALL, 5 );

    m_chooseDir = new wxButton( this, IDC_CHOOSE_DIR, wxT("Choose Dir"), wxDefaultPosition, wxDefaultSize, 0 );
    bSizer2->Add( m_chooseDir, 0, wxALIGN_CENTER|wxALL, 5 );

    bSizer1->Add( bSizer2, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    // Row 2: File list
    m_fileList = new wxListCtrl( this, IDC_FILE_LIST, wxDefaultPosition, wxSize(400,120), wxLC_LIST|wxLC_NO_HEADER|wxSUNKEN_BORDER|wxLC_SORT_ASCENDING );
    bSizer1->Add( m_fileList, 0, wxALIGN_CENTER|wxALL, 5 );

    // Row 3: File name + Save button
    wxBoxSizer* bSizer3 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText* staticText2 = new wxStaticText( this, ID_TEXT, wxT("File Name"), wxDefaultPosition, wxDefaultSize, 0 );
    bSizer3->Add( staticText2, 0, wxALIGN_CENTER|wxALL, 5 );

    m_fileName = new wxTextCtrl( this, IDC_FILE_NAME, wxT(""), wxDefaultPosition, wxSize(250,-1), 0 );
    bSizer3->Add( m_fileName, 0, wxALIGN_CENTER|wxALL, 5 );

    m_saveFile = new wxButton( this, IDC_SAVE_FILE, wxT("Save File"), wxDefaultPosition, wxDefaultSize, 0 );
    bSizer3->Add( m_saveFile, 0, wxALIGN_CENTER|wxALL, 5 );

    bSizer1->Add( bSizer3, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    // Row 4: File type choice
    wxBoxSizer* bSizer4 = new wxBoxSizer( wxHORIZONTAL );

    wxStaticText* staticText3 = new wxStaticText( this, ID_TEXT, wxT("File Type"), wxDefaultPosition, wxDefaultSize, 0 );
    bSizer4->Add( staticText3, 0, wxALIGN_CENTER|wxALL, 5 );

    m_fileType = new wxChoice( this, IDC_FILE_TYPE, wxDefaultPosition, wxDefaultSize, 0, (wxString*)NULL, 0 );
    bSizer4->Add( m_fileType, 0, wxALIGN_CENTER|wxALL, 5 );

    bSizer1->Add( bSizer4, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    // Row 5: Checkboxes in two columns
    wxBoxSizer* bSizer5 = new wxBoxSizer( wxHORIZONTAL );

    // Left column
    wxBoxSizer* bSizer6 = new wxBoxSizer( wxVERTICAL );

    m_saveTransformed = new wxCheckBox( this, IDC_SAVE_TRANSFORMED, wxT("Save Transformed"), wxDefaultPosition, wxDefaultSize, 0 );
    bSizer6->Add( m_saveTransformed, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    m_saveConnect = new wxCheckBox( this, IDC_SAVE_CONNECT, wxT("Save Connectivity"), wxDefaultPosition, wxDefaultSize, 0 );
    bSizer6->Add( m_saveConnect, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    m_saveFragments = new wxCheckBox( this, IDC_SAVE_FRAGMENTS, wxT("Save Fragments"), wxDefaultPosition, wxDefaultSize, 0 );
    m_saveFragments->SetValue( TRUE );
    bSizer6->Add( m_saveFragments, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    bSizer5->Add( bSizer6, 0, wxALIGN_CENTER|wxALL, 5 );

    // Right column
    wxBoxSizer* bSizer7 = new wxBoxSizer( wxVERTICAL );

    m_saveAmberPdb = new wxCheckBox( this, IDC_SAVE_AMBER_PDB, wxT("Save AMBER PDB"), wxDefaultPosition, wxDefaultSize, 0 );
    bSizer7->Add( m_saveAmberPdb, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    m_saveSepWatMol = new wxCheckBox( this, IDC_SAVE_SEP_WAT_MOL, wxT("Save Separate Water Molecules"), wxDefaultPosition, wxDefaultSize, 0 );
    bSizer7->Add( m_saveSepWatMol, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    bSizer5->Add( bSizer7, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5 );

    bSizer1->Add( bSizer5, 0, wxGROW|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

    this->SetSizer( bSizer1 );
    bSizer1->SetSizeHints( this );

    // Set file_types_ch for HaFileDlg1 compatibility
    file_types_ch = m_fileType;

    // Bind events
    m_saveFile->Bind( wxEVT_BUTTON, &SaveMolFileDlgBase::OnSaveFileClicked, this );
}

SaveMolFileDlgBase::~SaveMolFileDlgBase()
{
    m_saveFile->Unbind( wxEVT_BUTTON, &SaveMolFileDlgBase::OnSaveFileClicked, this );
}
