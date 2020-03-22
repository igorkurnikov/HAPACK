/***************************************************************************
 *   Nikolay Simakov                                                       *
 *   nsimakov@andrew.cmu.edu                                               *
 *   Maria Kurnikova Research Group                                        *
 *   http://crete.chem.cmu.edu/                                            *
 *   Carnegie Mellon University, 2006                                      *
 ***************************************************************************/

#define WXPYMOD_CPP

#include "wxpymod.h"
#include "ha_wx_res_mikola_wdr.h"
#include "harlemapp.h"

int wxPyMod::dlg_open = FALSE;

wxPyMod::wxPyMod(wxWindow* parent, wxWindowID id, const wxString& name,const wxPoint& pos, const wxSize& size)
    : wxFrame(parent, id, name, pos, size)
{
  this->SetExtraStyle(wxWS_EX_VALIDATE_RECURSIVELY);
  
  dlg_open = TRUE;
  
  wxMenuBar*  MenuBar = PyModMenuBarFunc();
  SetMenuBar(MenuBar);
  PyModDialogFunc(this,TRUE);
  SetSize(800,600);
  
  PyModTC = (wxTextCtrl*) FindWindow(IDT_PYMOD);
  PyModTC->AppendText(wxT(">>> "));
  LastPointOfEnter=4;
  CmdHistoryCount=0;
  
}

wxPyMod::~wxPyMod()
{
  fprintf(stdout,"wxPyMod::~wxPyMod\n");
  dlg_open = FALSE;
}

void
wxPyMod::OnClose(wxCloseEvent& event)
{
  wxPyMod::dlg_open = false;
  event.Skip();
}
void
wxPyMod::OnNewCmd(wxCommandEvent& event)
{
  if(PyModTC->GetRange(PyModTC->GetLastPosition()-1,PyModTC->GetLastPosition())!="\\")
  {
    wxString cmd = PyModTC->GetRange(LastPointOfEnter,PyModTC->GetLastPosition());
    int i=0;
    while(i<cmd.Length()-1)
    {
      if(cmd.GetChar(i)=='\\'&&cmd.GetChar(i+1)=='\n')
      {
        cmd.Remove(i,1);
      }
      else
      {
        i++;
      }
    }
    fprintf(stdout,"Wanted to execute: %s\n",cmd.ToStdString().c_str());
    CmdHistory.Add(cmd);
    CmdHistoryCount=0;
    pApp->cmd_pr.SetCmdLine( cmd.ToStdString() );
    pApp->ExecuteCommand();
    PyModTC->AppendText(wxT("\n>>> "));
    LastPointOfEnter=PyModTC->GetLastPosition();
  }
  else
  {
    PyModTC->AppendText(wxT("\n"));
  }
}
void
wxPyMod::InsertCmdFromHistory()
{
  if(CmdHistory.GetCount()>0)
  {
    int HisLine=CmdHistory.GetCount()+CmdHistoryCount;
    fprintf(stdout,"Up %d\n",HisLine);
    PyModTC->Remove(LastPointOfEnter,PyModTC->GetLastPosition());
    PyModTC->AppendText(CmdHistory[HisLine]);
  }
  else
  {
    CmdHistoryCount=0;
    PyModTC->Remove(LastPointOfEnter,PyModTC->GetLastPosition());
  }
}
void
wxPyMod::OnKeyUp(wxKeyEvent& event)
{
  if(PyModTC->IsEnabled())
  {
    long i,j,k=PyModTC->GetInsertionPoint();
    PyModTC->PositionToXY(k,&i,&j);
    switch (event.GetKeyCode())
    {
      case 317:
      {
        fprintf(stdout,"Up\n");
        if(PyModTC->GetNumberOfLines()>1)
        {
          PyModTC->SetInsertionPoint(PyModTC->XYToPosition(i,j+1));
          CmdHistoryCount--;
          if((-CmdHistoryCount)<=CmdHistory.GetCount())InsertCmdFromHistory();
          else CmdHistoryCount -= CmdHistory.GetCount();
        }
        else
        {
          PyModTC->SetInsertionPoint(LastPointOfEnter);
        }
        break;
      }
      case 319:
      {
        fprintf(stdout,"Down\n");
        if(j<PyModTC->GetNumberOfLines()-1)
        {
          PyModTC->SetInsertionPoint(PyModTC->XYToPosition(i,j-1));
          CmdHistoryCount++;
          if(CmdHistoryCount<0)InsertCmdFromHistory();
          else
          {
            CmdHistoryCount=0;
            PyModTC->Remove(LastPointOfEnter,PyModTC->GetLastPosition());
          }
        }
        break;
      }
      case 316:
      {
        fprintf(stdout,"left\n");
        if(k<LastPointOfEnter)
        {
          PyModTC->SetInsertionPoint(LastPointOfEnter);
        }
        break;
      }
    }
  }
}
BEGIN_EVENT_TABLE( wxPyMod, wxFrame )
    EVT_CLOSE   (wxPyMod::OnClose)
    EVT_KEY_UP(wxPyMod::OnKeyUp)
    EVT_TEXT_ENTER( IDT_PYMOD, wxPyMod::OnNewCmd )
END_EVENT_TABLE()
