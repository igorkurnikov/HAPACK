/*! \file hachart_wx.h

    wxWidgets Chart Classes in HARLEM Based on PLPlot package
 
    \author Igor Kurnikov  
    \date 2011-
*/

#if !defined(HACHART_WX_H)
#define HACHART_WX_H

#include <wx/window.h>
#include <wx/dcmemory.h>

#include "hachart.h"

#if WITH_PLPLOT
#include "plplot/wxPLplotstream.h"
#include "plplot/wxPLplotwindow.h"

#include "plplot/plstream.h"
#endif

#if WITH_PLPLOT
#if defined(PLPLOT_VERSION_MAJOR) and (PLPLOT_VERSION_MINOR > 12)
class HaChartWindow : public wxPLplotwindow<wxWindow>
#else
class HaChartWindow : public wxPLplotwindow
#endif
#else
class HaChartWindow : public wxFrame
#endif
{
public:
    HaChartWindow( wxFrame* frame, wxWindow* parent, wxWindowID id = -1, const wxPoint& pos = wxDefaultPosition,
                  const wxSize& size = wxDefaultSize, long style = 0,
		          bool useGraphicsContext = true);

    void OnChar( wxKeyEvent& event );

private:
    wxFrame* mframe;
};


class HaChartFrame  : public wxFrame
{
public:
    HaChartFrame( const wxString& title );
    void Plot();
	void PlotTest();

	HaChartPanel* GetPlot();

private:
    void OnQuit( wxCommandEvent& event );
    void OnAbout( wxCommandEvent& event );
    void OnBackgroundColor( wxCommandEvent& event );

private:
    HaChartWindow* plotwindow;
	HaChartPanel        plt;
    bool          bgcolor;
    int           m_backend;

    DECLARE_EVENT_TABLE()
};

#endif  // !defined HACHART_WX_H
