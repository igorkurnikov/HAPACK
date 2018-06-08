/*! \file hachart_wx.cpp

    wxWidgets Chart Classes in HARLEM Bases on PLPlot package implementation
 
    \author Igor Kurnikov  
    \date 2011-
*/

#define HACHART_WX_CPP

#include "wx/wxprec.h"
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include "hachart.h"
#include "hachart_wx.h" 

#include <cmath>

#define MAX( a, b )    ( ( a ) < ( b ) ? ( b ) : ( a ) )
#define MIN( a, b )    ( ( a ) < ( b ) ? ( a ) : ( b ) )

// Application icon as XPM
// This free icon was taken from http://2pt3.com/news/twotone-icons-for-free/
static const char *graph[] = {
// columns rows colors chars-per-pixel
    "16 16 4 2",
    "   c black",
    ".  c #BA1825",
    "X  c gray100",
    "UX c None",
// pixels
    "UX. . . . . . . . . . . . . . UX",
    ". . . . . . . . . . . . . . . . ",
    ". . . . . . . . . . . . . . . . ",
    ". . . . . . . . . . . X X . . . ",
    ". . . . . . . . . . . X X . . . ",
    ". . . . . . . . . . . X X . . . ",
    ". . . . . X X . . . . X X . . . ",
    ". . . . . X X . . . . X X . . . ",
    ". . . . . X X . X X . X X . . . ",
    ". . . . . X X . X X . X X . . . ",
    ". . . . . X X . X X . X X . . . ",
    ". . . . . X X . X X . X X . . . ",
    ". . . X X X X X X X X X X . . . ",
    ". . . . . . . . . . . . . . . . ",
    ". . . . . . . . . . . . . . . . ",
    "UX. . . . . . . . . . . . . . UX"
};


//--------------------------------------------------------------------------
// constants
//--------------------------------------------------------------------------
enum { wxPLplotDemo_Quit    = wxID_EXIT, wxPLplotDemo_About = wxID_ABOUT,
       wxPLplotDemo_BGColor = 10000 };

//--------------------------------------------------------------------------
// event tables and other macros for wxWidgets
//--------------------------------------------------------------------------
BEGIN_EVENT_TABLE( HaChartFrame, wxFrame )
EVT_MENU( wxPLplotDemo_Quit, HaChartFrame::OnQuit )
EVT_MENU( wxPLplotDemo_About, HaChartFrame::OnAbout )
EVT_MENU( wxPLplotDemo_BGColor, HaChartFrame::OnBackgroundColor )
END_EVENT_TABLE()


HaChartWindow::HaChartWindow( wxFrame* frame, wxWindow* parent, wxWindowID id, const wxPoint& pos,
                            const wxSize& size, long style, int pl_style ) :
    wxPLplotwindow( parent, id, pos, size, style, pl_style )
{
    mframe = frame;
}


void HaChartWindow::OnChar( wxKeyEvent& event )
{
    int keycode = event.GetKeyCode();

    if ( keycode == WXK_RETURN ||
         keycode == WXK_SPACE ||
         keycode == WXK_RIGHT ||
         keycode == WXK_ESCAPE )
        mframe->Close( true );
    else
        event.Skip();
}

HaChartFrame::HaChartFrame( const wxString& title ) : wxFrame( NULL, wxID_ANY, title )
{
    bgcolor = false;

    // add menu
    wxMenu *fileMenu = new wxMenu;
    fileMenu->Append( wxPLplotDemo_BGColor, _T( "&Change background color...\tAlt-C" ), _T( "Change background color" ) );
    fileMenu->Append( wxPLplotDemo_About, _T( "&About...\tF1" ), _T( "Show about dialog" ) );
    fileMenu->Append( wxPLplotDemo_Quit, _T( "E&xit\tAlt-X" ), _T( "Quit this program" ) );

    wxMenuBar *menuBar = new wxMenuBar();
    menuBar->Append( fileMenu, _T( "&File" ) );
    SetMenuBar( menuBar );
    SetIcon( wxIcon( graph ) );

    // add the wxPLplot
    wxPanel   * panel = new wxPanel( this );
    wxBoxSizer* box   = new wxBoxSizer( wxVERTICAL );

    plotwindow = new HaChartWindow( this, panel, -1, wxDefaultPosition, wxDefaultSize, wxWANTS_CHARS,
#if wxUSE_GRAPHICS_CONTEXT
        wxPLPLOT_BACKEND_GC | wxPLPLOT_DRAW_TEXT );
#else
        wxPLPLOT_BACKEND_AGG | wxPLPLOT_DRAW_TEXT );
#endif
    plotwindow->Connect( wxEVT_CHAR, wxKeyEventHandler( HaChartWindow::OnChar ) );
    box->Add( plotwindow, 1, wxALL | wxEXPAND, 0 );
    panel->SetSizer( box );
    SetSize( 640, 500 );      // set frame size
    SetSizeHints( 220, 150 ); // set minimum frame size

    wxString m_title = title;
    switch ( plotwindow->getBackend() )
    {
    case wxPLPLOT_BACKEND_DC:
        m_title += wxT( " (basic)" );
        break;
    case wxPLPLOT_BACKEND_GC:
        m_title += wxT( " (wxGC)" );
        break;
    case wxPLPLOT_BACKEND_AGG:
        m_title += wxT( " (AGG)" );
        break;
    default:
        break;
    }
    SetTitle( m_title );

    PlotTest();
}

void HaChartFrame::Plot()
{
	wxPLplotstream* pls = plotwindow->GetStream();
	plt.ToPLpStream(pls);
	plotwindow->RenewPlot();
}

HaChartPanel* HaChartFrame::GetPlot()
{
	return &plt;
}

void HaChartFrame::PlotTest()
{
	wxPLplotstream* pls = plotwindow->GetStream();

    const size_t  np = 500;
    PLFLT         x[np], y[np];
    PLFLT         xmin, xmax;
    PLFLT         ymin = 1e30, ymax = 1e-30;

    xmin = -2.0;
    xmax = 10.0;
    for ( size_t i = 0; i < np; i++ )
    {
        x[i] = ( xmax - xmin ) * i / np + xmin;
        y[i] = 1.0;
        if ( x[i] != 0.0 )
            y[i] = sin( x[i] ) / x[i];
        ymin = MIN( ymin, y[i] );
        ymax = MAX( ymax, y[i] );
    }

    pls->adv( 0 );
    if ( bgcolor )
    {
        pls->scol0( 0, 255, 255, 255 );
        pls->scol0( 15, 0, 0, 0 );
    }
    else
    {
        pls->scol0( 15, 255, 255, 255 );
        pls->scol0( 0, 0, 0, 0 );
    }
    pls->col0( 1 );
    pls->env( xmin, xmax, ymin, ymax, 0, 0 );
    pls->col0( 2 );
    pls->lab( "x", "y", "sin(x)/x" );

    pls->col0( 3 );
    pls->wid( 2 );
    pls->line( np, x, y );

    plotwindow->RenewPlot();
}


void HaChartFrame::OnQuit( wxCommandEvent& WXUNUSED( event ) )
{
    Close( true );
}


void HaChartFrame::OnBackgroundColor( wxCommandEvent& WXUNUSED( event ) )
{
    bgcolor = !bgcolor;
    Plot();
}


//! Show information if Menu entry About was choosen.
//
void HaChartFrame::OnAbout( wxCommandEvent& WXUNUSED( event ) )
{
    wxMessageBox( _T( "This is the About dialog of the wxPLplot demo.\n" ), _T( "About wxPLplot" ),
        wxOK | wxICON_INFORMATION, this );
}
