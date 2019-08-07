//
// C++ Implementation: wxgraph
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "wxgraphs.h"
#include <wx/dcbuffer.h>
#include <wx/filedlg.h>
#include <wx/textdlg.h>
#include <wx/string.h>
//#include <wxplot.h>

/* polar plot data */
#define PERIMETERPTS 100
/* Transformation function */

////////////////////////////////////////////////////////////////////////////////////////
//wxGraphXY
wxGraphXY::wxGraphXY(wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style, const wxString& name)
  : wxWindow(parent, id, pos, size, style, name)
{
  Title="Title";
  XAxis="x";
  YAxis="y";
  
  
  int i;
  N=101;
  Nd=0;
  nbin=100;
//#ifndef _WIN32
  x=new PLFLT[N];
  y=new PLFLT[N];
  for (i = 0; i < N; i++) {
    x[i] = 3.6 * i;
    y[i] = sin(x[i] * M_PI / 180.0);
  }
//#endif
  FindMinMax();
  SetBackgroundColour(*wxBLACK);
}


wxGraphXY::~wxGraphXY()
{
}
void wxGraphXY::OnPaint(wxPaintEvent& event)
{
  //wxBufferedPaintDC dc(this);
  wxPaintDC dc(this);
  //wxBrush brush(*wxBLACK_BRUSH);
  
  //dc.SetBackgroundMode(wxSOLID);
  //dc.SetBackground(*wxBLACK_BRUSH);

  Plot(dc,0,0);
}
void wxGraphXY::OnSize(wxSizeEvent& event)
{
  Refresh();
}
void wxGraphXY::Plot(wxDC& dc, int w_, int h_)
{
//#ifndef _WIN32
  int w,h;
  this->GetSize(&w,&h);

  if(w_!=0 && h_!=0)
  {
    w=w_;
    h=h_;
  }
  
  wxPlot *wxplot=new wxPlot(&dc, w, h);

  PLINT space0 = 0, mark0 = 0, space1 = 1500, mark1 = 1500;
  int i;

  wxplot->adv(0);
  
  wxplot->vsta();
  wxplot->scolbg(0,0,0);
  if(Nd>0)
  {
    wxplot->hist(Nd,d,dmin,dmax,nbin,0);
  }
  else
  {
    wxplot->wind(xmin,xmax,ymin,ymax);
    
    wxplot->col0(1);
    wxplot->box("bcnst", (xmax-xmin)/5.0, 5, "bcnst", (ymax-ymin)/5.0, 5);
    wxplot->styl(1, &mark1, &space1);
    wxplot->col0(2);
    //wxplot->box("g", 30.0, 0, "g", 0.2, 0);
    wxplot->box("g", (xmax-xmin)/25.0, 0, "g", (ymax-ymin)/25.0, 0);
    wxplot->styl(0, &mark0, &space0);
  
    wxplot->col0(3);
    wxplot->lab(XAxis.c_str(), YAxis.c_str(), Title.c_str());
  
  
    wxplot->col0(4);
    wxplot->line(101, x, y);
  }
  delete(wxplot);
//#endif
}
void wxGraphXY::FindMinMax()
{
//#ifndef _WIN32
  if(N>0)
  {
    xmin=x[0];
    xmax=x[0];
    ymin=y[0];
    ymax=y[0];
    for (int i = 0; i < N; i++)
    {
      if(x[i]<xmin)xmin=x[i];
      if(x[i]>xmax)xmax=x[i];
      if(y[i]<ymin)ymin=y[i];
      if(y[i]>ymax)ymax=y[i];
    }
  }
  else
  {
    xmin=0.0;
    xmax=0.0;
    ymin=0.0;
    ymax=0.0;
  }
  if(Nd>0)
  {
    dmin=d[0];
    dmax=d[0];
    for (int i = 0; i < Nd; i++)
    {
      if(d[i]<dmin)dmin=d[i];
      if(d[i]>dmax)dmax=d[i];
    }
  }
//#endif
}
void wxGraphXY::OnRightUp(wxMouseEvent& event)
{
  wxMenu menu(wxT("wxGraph:"));
  menu.Append(IDB_PNP_SAVEXY, wxT("Save (X,Y)"));
  if(Nd>0)
  {
    menu.Append(IDB_PNP_NBIN, wxT("Change nbin"));
    menu.Append(IDB_PNP_DMIN, wxT("Change min"));
    menu.Append(IDB_PNP_DMAX, wxT("Change max"));
  }
  PopupMenu(&menu,  event.GetPosition());
}
void wxGraphXY::OnSaveXY(wxCommandEvent& event)
{
  wxFileDialog dialog(this,
                      wxT("Save (X,Y)"),
                      wxEmptyString,
                      wxT("xy.dat"),
                      wxT("dat file  (*.dat)|*.dat"),
                      wxSAVE|wxOVERWRITE_PROMPT);
  if (dialog.ShowModal() == wxID_OK)
  {
    //wxString FileNameTree=dialog.GetFilename();
//#ifndef _WIN32
    FILE *out;
    out=fopen(dialog.GetFilename().c_str(),"w");
    int i;
    for (i = 0; i < N; i++)
    {
      fprintf(out,"%8.7e %8.7e\n",x[i],y[i]);
    }
    fclose(out);
//#endif
  }
}
void wxGraphXY::SetLabels(char* strT,char* strX,char* strY)
{
  Title=strT;
  XAxis=strX;
  YAxis=strY;
  Refresh();
}
void wxGraphXY::SetXY(HaVec_float* X,HaVec_float* Y)
{
//#ifndef _WIN32
  N=X->size();
  delete x;
  delete y;
  x=new PLFLT[N];
  y=new PLFLT[N];
  int i;
  for (i = 0; i < N; i++) {
    x[i] = (*X)[i];
    y[i] = (*Y)[i];
  }
  FindMinMax();
  Refresh();
//#endif
}
void wxGraphXY::SetXY(std::vector<float> X,std::vector<float> Y)
{
//#ifndef _WIN32
  N=X.size();
  delete x;
  delete y;
  x=new PLFLT[N];
  y=new PLFLT[N];
  int i;
  for (i = 0; i < N; i++) {
    x[i] = X[i];
    y[i] = Y[i];
  }
  FindMinMax();
  Refresh();
//#endif
}
void wxGraphXY::SetDataHist(std::vector<float> D)
{
//#ifndef _WIN32
  Nd=D.size();
  delete d;
  d=new PLFLT[Nd];
  for (int i = 0; i < Nd; i++)d[i]=D[i];
  FindMinMax();
  Refresh();
//#endif
}
void wxGraphXY::OnSetNbin(wxCommandEvent& event)
{
  wxTextEntryDialog dialog(this,"nbin","Histogram Stuff","100",wxOK | wxCANCEL);
  if (dialog.ShowModal() == wxID_OK)
  {
    wxString str=dialog.GetValue();
    long l;
    if(str.ToLong(&l))nbin=(int)l;
  }
  Refresh();
}
void wxGraphXY::OnSetDMin(wxCommandEvent& event)
{
  wxTextEntryDialog dialog(this,"min","Histogram Stuff","100",wxOK | wxCANCEL);
  if (dialog.ShowModal() == wxID_OK)
  {
    wxString str=dialog.GetValue();
    double dtmp;
    if(str.ToDouble(&dtmp))dmin=(float)dtmp;
  }
  Refresh();
}
void wxGraphXY::OnSetDMax(wxCommandEvent& event)
{
  wxTextEntryDialog dialog(this,"max","Histogram Stuff","100",wxOK | wxCANCEL);
  if (dialog.ShowModal() == wxID_OK)
  {
    wxString str=dialog.GetValue();
    double dtmp;
    if(str.ToDouble(&dtmp))dmax=(float)dtmp;
  }
  Refresh();
}
BEGIN_EVENT_TABLE(wxGraphXY, wxWindow)
    EVT_PAINT(wxGraphXY::OnPaint)
    EVT_SIZE (wxGraphXY::OnSize)
    EVT_MENU (IDB_PNP_SAVEXY,wxGraphXY::OnSaveXY)
    EVT_MENU (IDB_PNP_NBIN,wxGraphXY::OnSetNbin)
    EVT_MENU (IDB_PNP_DMIN,wxGraphXY::OnSetDMin)
    EVT_MENU (IDB_PNP_DMAX,wxGraphXY::OnSetDMax)
    EVT_RIGHT_UP(wxGraphXY::OnRightUp)
END_EVENT_TABLE()
    ///////////////////////////////////////////////////////////////////
//wxPlanViewFrame
    wxGraphXYFrame::wxGraphXYFrame(wxWindow* parent, wxWindowID id, const wxString& name,const wxPoint& pos, const wxSize& size)
  : wxFrame(parent, id, name, pos, size)
{
  //Create Menu
  wxMenuBar*  MenuBar = new wxMenuBar;
  wxMenu* item1 = new wxMenu;
  item1->Append( -1, wxT("sthm"), wxT("") );
  MenuBar->Append( item1, wxT("&File") );
  SetMenuBar(MenuBar);
  GraphXY=new wxGraphXY(this);
}
wxGraphXYFrame::~wxGraphXYFrame()
{
}
BEGIN_EVENT_TABLE(wxGraphXYFrame , wxFrame )
END_EVENT_TABLE()

///////////////////////////////////////////////////////////////////
wxWindow* wxGraphs::Parent=NULL;
void wxGraphs::GraphXY()
{
  wxGraphXYFrame* graph=new wxGraphXYFrame(wxGraphs::Parent);
  
  graph->Show(true);
}
void wxGraphs::GraphXY(HaVec_float* Y)
{
  HaVec_float X;
  X.newsize(Y->size());
  float i;
  for(i=0;i<Y->size();i++)X[i]=i;
  
  wxGraphXYFrame* graph=new wxGraphXYFrame(wxGraphs::Parent);
  graph->GraphXY->SetXY( &X,Y);
  graph->Show(true);
}
