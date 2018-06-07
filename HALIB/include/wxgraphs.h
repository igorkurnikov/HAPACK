//
// C++ Interface: wxgraph
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef WXGRAPHS_H
#define WXGRAPHS_H
#include <wx/wx.h>
#include <wx/scrolwin.h>
#include <wx/panel.h>
#include <wx/window.h>
#include <wxPlot.h>
#include <string>
#include <vector>

#ifndef HALIB_BASE
#include "halinalg.h"
#endif
const int IDB_PNP_SAVEXY=7001;
const int IDB_PNP_NBIN=7002;
const int IDB_PNP_DMIN=7003;
const int IDB_PNP_DMAX=7004;

/**
@author Nikolay Simakov
*/
class wxGraphXY : public wxWindow
{
  public:
    wxGraphXY(wxWindow* parent, wxWindowID id = -1, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxTAB_TRAVERSAL, const wxString& name = "panel");

    ~wxGraphXY();
    void Plot(wxDC& dc, int w_, int h_);
    void OnPaint(wxPaintEvent& event);
    void OnSize(wxSizeEvent& event);
    void OnSaveXY(wxCommandEvent& event);
    void OnSetNbin(wxCommandEvent& event);
    void OnSetDMin(wxCommandEvent& event);
    void OnSetDMax(wxCommandEvent& event);
    void OnRightUp(wxMouseEvent& event);
    void SetLabels(char* strT,char* strX,char* strY);
    void SetXY(std::vector<float> X,std::vector<float> Y);
    void SetXY(HaVec_float* X,HaVec_float* Y);
    void SetDataHist(std::vector<float> D);
    
    
  private:
    DECLARE_EVENT_TABLE();
    std::string Title;
	std::string XAxis, YAxis;
    int N;
    int Nd,nbin;

    PLFLT *x, *y, *d;
    PLFLT xmin,xmax,ymin,ymax,dmin,dmax;

    void FindMinMax();
    
    
};
/**
@author Nikolay Simakov
 */
class wxGraphXYFrame : public wxFrame
{
  public:
    wxGraphXYFrame(wxWindow* parent, wxWindowID id = -1, const wxString& name="wxPlanViewFrame",const wxPoint& pos=wxDefaultPosition, const wxSize& size=wxDefaultSize);
    ~wxGraphXYFrame();
    
    wxGraphXY* GraphXY;
  protected:
    
    DECLARE_EVENT_TABLE();
};
//!This Class is desighned speasially for calling from python script wxGraphs.GraphXY()
class wxGraphs
{
  public:
    wxGraphs(){}
    ~wxGraphs(){}
    static wxWindow* Parent;
    static void GraphXY();
    static void GraphXY(HaVec_float* Y);
};
#endif
