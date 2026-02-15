///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version 4.2.1-115-g11c2dec8)
// http://www.wxformbuilder.org/
//
// PLEASE DO *NOT* EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

#include "MyProjectBase.h"

///////////////////////////////////////////////////////////////////////////

SaveMolFileDlgBase::SaveMolFileDlgBase( wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style ) : wxDialog( parent, id, title, pos, size, style )
{
	this->SetSizeHints( wxDefaultSize, wxDefaultSize );

	wxBoxSizer* ;
	= new wxBoxSizer( wxVERTICAL );

	wxBoxSizer* ;
	= new wxBoxSizer( wxHORIZONTAL );

	IDC_DIR_LABEL = new wxStaticText( this, wxID_ANY, _("Directory:"), wxDefaultPosition, wxDefaultSize, 0 );
	IDC_DIR_LABEL->Wrap( 0 );
	->Add( IDC_DIR_LABEL, 0, wxALIGN_CENTER|wxALL, 5 );

	IDC_DIR_NAME = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxSize( 250,-1 ), 0 );
	->Add( IDC_DIR_NAME, 0, wxALIGN_CENTER|wxALL, 5 );

	IDC_CHOOSE_DIR = new wxButton( this, wxID_ANY, _("Choose Dir"), wxDefaultPosition, wxDefaultSize, 0 );

	IDC_CHOOSE_DIR->SetDefault();
	IDC_CHOOSE_DIR->SetAuthNeeded();
	->Add( IDC_CHOOSE_DIR, 0, wxALIGN_CENTER|wxALL, 5 );


	->Add( , 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

	IDC_FILE_LIST = new wxListCtrl( this, wxID_ANY, wxDefaultPosition, wxSize( 400,120 ), wxLC_LIST|wxLC_NO_HEADER|wxSUNKEN_BORDER|wxLC_SORT_ASCENDING );
	->Add( IDC_FILE_LIST, 0, wxALIGN_CENTER|wxALL, 5 );

	wxBoxSizer* ;
	= new wxBoxSizer( wxHORIZONTAL );

	IDC_FILENAME_LABEL = new wxStaticText( this, wxID_ANY, _("File Name"), wxDefaultPosition, wxDefaultSize, 0 );
	IDC_FILENAME_LABEL->Wrap( 0 );
	->Add( IDC_FILENAME_LABEL, 0, wxALIGN_CENTER|wxALL, 5 );

	IDC_FILE_NAME = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxSize( 250,-1 ), 0 );
	->Add( IDC_FILE_NAME, 0, wxALIGN_CENTER|wxALL, 5 );

	IDC_SAVE_FILE = new wxButton( this, wxID_ANY, _("Save File"), wxDefaultPosition, wxDefaultSize, 0 );

	IDC_SAVE_FILE->SetDefault();
	IDC_SAVE_FILE->SetAuthNeeded();
	->Add( IDC_SAVE_FILE, 0, wxALIGN_CENTER|wxALL, 5 );


	->Add( , 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

	wxBoxSizer* ;
	= new wxBoxSizer( wxHORIZONTAL );

	IDC_FILETYPE_LABEL = new wxStaticText( this, wxID_ANY, _("File Type"), wxDefaultPosition, wxDefaultSize, 0 );
	IDC_FILETYPE_LABEL->Wrap( 0 );
	->Add( IDC_FILETYPE_LABEL, 0, wxALIGN_CENTER|wxALL, 5 );

	wxArrayString IDC_FILE_TYPEChoices;
	IDC_FILE_TYPE = new wxChoice( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, IDC_FILE_TYPEChoices, 0 );
	IDC_FILE_TYPE->SetSelection( 0 );
	->Add( IDC_FILE_TYPE, 0, wxALIGN_CENTER|wxALL, 5 );

	test = new wxButton( this, wxID_ANY, _("MyButton"), wxDefaultPosition, wxDefaultSize, 0 );
	->Add( test, 0, wxALL, 5 );


	->Add( , 0, wxEXPAND|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

	wxBoxSizer* ;
	= new wxBoxSizer( wxHORIZONTAL );

	wxBoxSizer* ;
	= new wxBoxSizer( wxVERTICAL );

	IDC_SAVE_TRANSFORMED = new wxCheckBox( this, wxID_ANY, _("Save Transformed"), wxDefaultPosition, wxDefaultSize, 0 );
	IDC_SAVE_TRANSFORMED->SetValue(true);
	->Add( IDC_SAVE_TRANSFORMED, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

	IDC_SAVE_CONNECT = new wxCheckBox( this, wxID_ANY, _("Save Connectivity"), wxDefaultPosition, wxDefaultSize, 0 );
	IDC_SAVE_CONNECT->SetValue(true);
	->Add( IDC_SAVE_CONNECT, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

	IDC_SAVE_FRAGMENTS = new wxCheckBox( this, wxID_ANY, _("Save Fragments"), wxDefaultPosition, wxDefaultSize, 0 );
	IDC_SAVE_FRAGMENTS->SetValue(true);
	->Add( IDC_SAVE_FRAGMENTS, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );


	->Add( , 0, wxALIGN_CENTER|wxALL, 5 );

	wxBoxSizer* ;
	= new wxBoxSizer( wxVERTICAL );

	IDC_SAVE_AMBER_PDB = new wxCheckBox( this, wxID_ANY, _("Save AMBER PDB"), wxDefaultPosition, wxDefaultSize, 0 );
	IDC_SAVE_AMBER_PDB->SetValue(true);
	->Add( IDC_SAVE_AMBER_PDB, 0, wxEXPAND|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

	IDC_SAVE_SEP_WAT_MOL = new wxCheckBox( this, wxID_ANY, _("Save Separate Water Molecules"), wxDefaultPosition, wxDefaultSize, 0 );
	IDC_SAVE_SEP_WAT_MOL->SetValue(true);
	->Add( IDC_SAVE_SEP_WAT_MOL, 0, wxEXPAND|wxALIGN_CENTER_VERTICAL|wxALL, 5 );


	->Add( , 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5 );


	->Add( , 0, wxEXPAND|wxALIGN_CENTER_VERTICAL|wxALL, 5 );


	this->SetSizer(  );
	this->Layout();
	->Fit( this );

	this->Centre( wxBOTH );
}

SaveMolFileDlgBase::~SaveMolFileDlgBase()
{
}
