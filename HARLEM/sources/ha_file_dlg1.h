/*! \file ha_file_dlg1.h

    Base class for file dialogs in HARLEM (extracted from dialogs_wx_1.h)

    \author Igor Kurnikov
    \date 2003-
*/

#if !defined(HA_FILE_DLG1_H)
#define HA_FILE_DLG1_H

#include <wx/dialog.h>

class wxChoice;
class wxCommandEvent;
class wxCloseEvent;
class wxListEvent;

//----------------------------------------------------------------------------
// HaFileDlg1
//----------------------------------------------------------------------------

class HaFileDlg1: public wxDialog
{
public:
    // constructors and destructors
    HaFileDlg1( wxWindow *parent, wxWindowID id, const wxString &title,
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxDEFAULT_DIALOG_STYLE );

	//virtual ~HaFileDlg1();

	void OnClose(wxCloseEvent& event);

	wxChoice* file_types_ch; //!< choice box for file types choices

protected:
    // WDR: member variable declarations for HaFileDlg1
    int nsubdir; //!< number of subdirectories including .. and . listed in the ListCtrl

protected:
    // WDR: handler declarations for HaFileDlg1
    void OnChangeFileType( wxCommandEvent &event );
    void OnSelectFile( wxListEvent &event );
	void OnActivateFile( wxListEvent &event );
    void ChooseDir( wxCommandEvent &event );
    void OnLoadFile( wxCommandEvent &event );

    virtual void FillFileTypes(); //!< Fill File Types filters choice box
    void OnChangeDir();   //!< Fill file list in the current directory

private:
    DECLARE_EVENT_TABLE()
};

#endif // HA_FILE_DLG1_H
