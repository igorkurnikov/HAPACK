#include "oglview.h"
//#include <osgProducer/Viewer>

#include <wx/wx.h>

#include <stdio.h>
#include <unistd.h>
#include <X11/X.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xmd.h>
#include <X11/cursorfont.h>

#include <gdk/gdkx.h>

#include <wx/defs.h>
#include <wx/dcclient.h>

#include <gtk/gtkwidget.h>
#include <gdk/gdktypes.h>

#include <osgDB/ReadFile>
#include <osgDB/WriteFile>
#include <osgUtil/Optimizer>
#include <osgProducer/Viewer>

#include <iostream>

using namespace std;
BEGIN_EVENT_TABLE(wxOsgProducerViewer,wxWindow)
    EVT_ERASE_BACKGROUND(wxOsgProducerViewer::OnEraseBackground)
    EVT_WINDOW_CREATE(wxOsgProducerViewer::OnCreate)
    EVT_SIZE(wxOsgProducerViewer::OnResize)
    EVT_PAINT(wxOsgProducerViewer::OnPaint)
END_EVENT_TABLE()
    wxOsgProducerViewer::wxOsgProducerViewer(wxWindow* parent, wxWindowID id, 
                                             const wxPoint& pos, const wxSize& size, 
                                             long style, const wxString& name)
  : wxWindow(parent,id,pos,size,
             style
                 |wxWS_EX_TRANSIENT
                 |wxWS_EX_PROCESS_UI_UPDATES
                 ,name)
{
  // we dont need it, as everything is painted by ourselves.
  SetAutoLayout(false);
  
  m_viewer.setUpViewer(osgProducer::Viewer::STANDARD_SETTINGS);
  GtkWidget *gtkwidget = this->GetHandle();
  GdkWindow *window = gtkwidget->window;
  
  Producer::RenderSurface * rs = m_viewer.getCamera(0)->getRenderSurface();
  
  rs->setDisplay( GDK_WINDOW_XDISPLAY(window) ); 
  rs->setParentWindow(  GDK_WINDOW_XID(window) ); 
  // The call to fullScreen seems to be bugous, and leads to occasional crashes
  // Relying on XFunctions instead for resizing, The viewport is somehow
  // updated correctly.  
  // rs->fullScreen(false);
  
  m_viewer.setSceneData(osgDB::readNodeFile("/home/mikola/MYPROG/dep/OSG_OP_OT/OpenSceneGraph-Data/SmokeBox.osg"));
   
  cout << "realizing" << endl;
  m_viewer.realize();
  cout << "finished OnInit" << endl;
}

void wxOsgProducerViewer::OnEraseBackground(wxEraseEvent& WXUNUSED(event)) {}


wxOsgProducerViewer::~wxOsgProducerViewer()
{
  cout << "wxOsgProducerViewer" << endl;
  m_viewer.setDone(true);
  m_viewer.sync();
}

void wxOsgProducerViewer::OnResize(wxSizeEvent& event) 
{
  unsigned int w, h;
  GetClientSize((int*)&w, (int*)&h);

  Producer::RenderSurface * rs = m_viewer.getCamera(0)->getRenderSurface();
  // The method: RenderSurface->setWindowRectangle misplaces the window 
  // (y-coord ist constantly changing) --> Directly call the X11 Method
  XResizeWindow(rs->getDisplay(), 
                rs->getWindow(), w, h);

  Refresh();
  Update();
}

void wxOsgProducerViewer::OnPaint(wxPaintEvent& WXUNUSED(event))
{
  m_viewer.sync();
  m_viewer.update();
  m_viewer.frame();
}

void wxOsgProducerViewer::OnCreate(wxWindowCreateEvent& WXUNUSED(event))
{
}
///////////////////////////////////////////////////////////////////
//
OGLViewFrame::OGLViewFrame(wxWindow* parent, wxWindowID id, const wxString& name,const wxPoint& pos, const wxSize& size)
  : wxFrame(parent, id, name, pos, size)
{
  //Create Menu
  wxMenuBar*  MenuBar = new wxMenuBar;
  wxMenu* item1 = new wxMenu;
  item1->Append( -1, wxT("sthm"), wxT("") );
  MenuBar->Append( item1, wxT("&File") );
  SetMenuBar(MenuBar);
  
  Show(true);
  ProducerViewer = new wxOsgProducerViewer(this, -1);
  
}
OGLViewFrame::~OGLViewFrame()
{
}
BEGIN_EVENT_TABLE(OGLViewFrame , wxFrame )
    
END_EVENT_TABLE()
