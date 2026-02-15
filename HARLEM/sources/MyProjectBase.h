///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version 4.2.1-115-g11c2dec8)
// http://www.wxformbuilder.org/
//
// PLEASE DO *NOT* EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

#pragma once

#include <wx/artprov.h>
#include <wx/xrc/xmlres.h>
#include <wx/intl.h>
#include <wx/string.h>
#include <wx/stattext.h>
#include <wx/gdicmn.h>
#include <wx/font.h>
#include <wx/colour.h>
#include <wx/settings.h>
#include <wx/textctrl.h>
#include <wx/button.h>
#include <wx/bitmap.h>
#include <wx/image.h>
#include <wx/icon.h>
#include <wx/sizer.h>
#include <wx/listctrl.h>
#include <wx/choice.h>
#include <wx/checkbox.h>
#include <wx/dialog.h>

///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// Class SaveMolFileDlgBase
///////////////////////////////////////////////////////////////////////////////
class SaveMolFileDlgBase : public wxDialog
{
	private:

	protected:
		wxStaticText* IDC_DIR_LABEL;
		wxTextCtrl* IDC_DIR_NAME;
		wxButton* IDC_CHOOSE_DIR;
		wxListCtrl* IDC_FILE_LIST;
		wxStaticText* IDC_FILENAME_LABEL;
		wxTextCtrl* IDC_FILE_NAME;
		wxButton* IDC_SAVE_FILE;
		wxStaticText* IDC_FILETYPE_LABEL;
		wxChoice* IDC_FILE_TYPE;
		wxButton* test;
		wxCheckBox* IDC_SAVE_TRANSFORMED;
		wxCheckBox* IDC_SAVE_CONNECT;
		wxCheckBox* IDC_SAVE_FRAGMENTS;
		wxCheckBox* IDC_SAVE_AMBER_PDB;
		wxCheckBox* IDC_SAVE_SEP_WAT_MOL;

	public:

		SaveMolFileDlgBase( wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = _("Save Molecular File"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxDEFAULT_DIALOG_STYLE|wxRESIZE_BORDER );

		~SaveMolFileDlgBase();

};

