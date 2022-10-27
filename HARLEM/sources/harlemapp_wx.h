/*! \file harlemapp_wx.h

    Main wxWidget classes for HARLEM application
    \author Igor Kurnikov  
    \date 2003-
*/

#if !defined(HARLEMAPP_WX_H)
#define HARLEMAPP_WX_H

#include <wx/app.h>
#include "hawx_add.h"
#include "harlemapp.h"

class HaMainFrameWX;

class HarlemAppWX: public wxApp, public HarlemApp
{ 
public:
  HarlemAppWX();
  virtual ~HarlemAppWX();

// Virtuals of wxApp: 

  virtual bool OnInit(void);
  virtual bool ProcessEvent(wxEvent& event); //!< overrides wxApp::ProcessEvent() function to add processing of Molecular Sets
  virtual int  OnExit();
  
  void OnExitApp(wxCommandEvent& event); //!< Exit Application

#if !defined(HA_NOGUI)
  void OnIdle(wxIdleEvent& event);

  HaMainFrameWX* GetMainFrame() { return m_mainFrame;}
  HaMainFrameWX* m_mainFrame;
  wxFrame*       m_cmd_win; 
#endif

  virtual int OnRun(); //!< Need to define for non gui applications

private:
  int LoadHaPyGUIModules();//!< Load harlem python GUI
  DECLARE_EVENT_TABLE();
};

extern const wxEventType wxEVT_HARLEM_APP;
const int HARLEM_APP_EXIT_ID = 20001;

DECLARE_APP(HarlemAppWX);


#endif // !defined(HARLEMAPP_WX_H)
