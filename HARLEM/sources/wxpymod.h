/***************************************************************************
 *   Nikolay Simakov                                                       *
 *   nsimakov@andrew.cmu.edu                                               *
 *   Maria Kurnikova Research Group                                        *
 *   http://crete.chem.cmu.edu/                                            *
 *   Carnegie Mellon University, 2006                                      *
 ***************************************************************************/

#ifndef WXPYMOD_H
#define WXPYMOD_H
#include <wx/wx.h>
#include <wx/arrstr.h>
#include "hacompmod.h"
class wxTextCtrl;

/*! Class for Python Interactive
*/
class wxPyMod : public wxFrame
{
  public:
    wxPyMod(wxWindow* parent, wxWindowID id = -1, const wxString& name="Python Iteractive",const wxPoint& pos=wxDefaultPosition, const wxSize& size=wxDefaultSize);
    virtual ~wxPyMod();
    
    static int dlg_open;
    void OnClose(wxCloseEvent& event);
    
    void OnNewCmd(wxCommandEvent& event);
    void OnKeyUp(wxKeyEvent& event);
    void InsertCmdFromHistory();
  private:
    int LastPointOfEnter;
    wxTextCtrl* PyModTC;
    wxArrayString CmdHistory;
    int CmdHistoryCount;
    DECLARE_EVENT_TABLE()
};
#endif // end if !defined(PNPMOD_H)
