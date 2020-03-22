
#ifndef OGLVIEW_H
#define OGLVIEW_H
#include <wx/wx.h>
#include <wx/window.h>

#include <wx/glcanvas.h>
#include <osgProducer/Viewer>
/**
@author Nikolay Simakov
*/

class wxOsgProducerViewer:public wxWindow{
  public:
    wxOsgProducerViewer(wxWindow* parent, wxWindowID id, 
                        const wxPoint& pos = wxDefaultPosition, 
                        const wxSize& size = wxDefaultSize, 
                        long style = 0, const wxString& name = wxPanelNameStr);
    ~wxOsgProducerViewer();
    void OnCreate(wxWindowCreateEvent& WXUNUSED(event));
    void OnResize(wxSizeEvent& event);
    void OnPaint(wxPaintEvent& WXUNUSED(event));
    void OnEraseBackground(wxEraseEvent& WXUNUSED(event));
    osgProducer::Viewer* getViewer(){return &m_viewer;}

  protected:
    DECLARE_EVENT_TABLE()

    osgProducer::Viewer m_viewer;
};

class OGLViewFrame : public wxFrame
{
  public:
    OGLViewFrame(wxWindow* parent, wxWindowID id = -1, const wxString& name="OGLViewFrame",const wxPoint& pos=wxDefaultPosition, const wxSize& size=wxDefaultSize);
    ~OGLViewFrame();
    
    wxOsgProducerViewer* ProducerViewer;
    
  protected:
    
    DECLARE_EVENT_TABLE();
};



#endif
