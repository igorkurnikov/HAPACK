/*! file harlemapp_wx.cpp

    wxWindows class for HARLEM application                 

    \author Igor Kurnikov 
    \date 2003-

*/

#define HARLEMAPP_WX_CPP

#include "haconst.h"

#include "mpi.h"
#if !defined(HARLEM_PYTHON_NO)
#include "Python.h"
#endif

#include "haconst.h"
#include "vec3d.h"

#if !defined(HA_NOGUI)
#	include <wx/wxprec.h> 
#	include <wx/wx.h>
#	include <wx/file.h>
#	include <wx/textfile.h>
#	include <wx/textctrl.h>
#	include "ha_wx_res_wdr.h"
#endif

#include "hamolview.h"
#include "hampi.h"
#include "harlemapp.h"

#include "hawx_add.h"
#include "harlemapp_wx.h"
#include "molset_evt_handler.h"

#if !defined(HA_NOGUI)
#	include "hamainframe_wx.h"
#	include "dialogs_wx_1.h"
#endif

#include "hasvnrev.h"
#if defined(_MSC_VER)
#	include <windows.h>
#endif

#if defined(HA_NOGUI)
//IMPLEMENT_APP_CONSOLE(HarlemAppWX)
//IMPLEMENT_APP_NO_MAIN(HarlemAppWX)
//
//int main(int argc, char **argv) 
//{
//	int ires = 0;
//	printf("main() pt 1 \n");
//	ires = wxEntry(argc, argv);
//	printf("main() pt 2 ires = %d \n\n",ires);
//	return ires;
//}
IMPLEMENT_APP(HarlemAppWX)
#else
IMPLEMENT_APP(HarlemAppWX)
#endif

const wxEventType wxEVT_HARLEM_APP = wxNewEventType();

BEGIN_EVENT_TABLE( HarlemAppWX,  wxApp)
#if !defined(HA_NOGUI)
    EVT_IDLE( HarlemAppWX::OnIdle)
#endif
    EVT_COMMAND( HARLEM_APP_EXIT_ID, wxEVT_HARLEM_APP, HarlemAppWX::OnExitApp)
END_EVENT_TABLE()

HarlemAppWX::HarlemAppWX()
{
#if !defined(HA_NOGUI)
	m_mainFrame = NULL;
	m_cmd_win   = NULL;
#endif
}

HarlemAppWX::~HarlemAppWX()
{
//	if( m_mainFrame != NULL ) delete m_mainFrame;
}

bool HarlemAppWX::ProcessEvent(wxEvent& event)
{
	if( wxApp::ProcessEvent(event)) return true;

	HaMolSet* pmset = GetCurMolSet();

	if(pmset != NULL) if( pmset->p_evt_h->ProcessEvent(event) ) return true; 

	int nm = molset_vec.size();
	int i;
	for( i = 0; i < nm; i++)
	{
		pmset = (HaMolSet*) molset_vec[i];
		if( pmset == GetCurMolSet() ) continue;
		if( pmset->p_evt_h->ProcessEvent(event) ) return true; 
	}
	return false;
}

#if !defined(HA_NOGUI)
void HarlemAppWX::OnIdle(wxIdleEvent& event)
{
//	int ret = PyRun_InteractiveOne(stdin, "<stdin>");  //Igor Python incorporation...
  
	//fprintf(stdout,"HaMolSet::OnIdle >\n");
  //PyRun_InteractiveOne(stdin,"none");
  //fprintf(stdout,"HaMolSet::OnIdle <\n");
      
  /*
  unsigned char uc[1024];
  //uc=fgetc(stdin);
  int err;
  size_t l=fread(uc,1024,1,stdin);
  if(l>0)
  {
    int i;
    for(i=0;i<l;i++)
      PyCmd[inxPyCmd+i]=uc[i];
    inxPyCmd=inxPyCmd+l;
    PyCmd[inxPyCmd]='\0';
    fprintf(stdout,"HaMolSet::OnIdle [%s]\n",PyCmd);
    
  }
  else
  {
    fprintf(stdout,"HaMolSet::OnIdle No input\n");
    clearerr(stdin);
  }*/
}
#endif

void HarlemAppWX::OnExitApp(wxCommandEvent& event)
{
	delete this;
}

#if !defined(HARLEM_PYTHON_NO)
class PythonThread: public wxThread
{
public:
	PythonThread() {}
	virtual ExitCode Entry()
	{
		int sts = PyRun_AnyFile(stdin, "<stdin>");		
		return 0;
	}
	void OnExit()
	{
		Py_Finalize();
	}
};
#endif

bool HarlemAppWX::OnInit(void)
{
	this->argc_loc = this->argc;
	this->argv_loc = this->argv;

// IGOR TMP:

	wxLog* p_log = new wxLogStderr();
	wxLog::SetActiveTarget(p_log);

#if !defined(HA_NOGUI)
	SetExitOnFrameDelete(true);   // Exit app when the top level frame is deleted
	if (gui_mode)
	{
		m_mainFrame = new HaMainFrameWX();

		m_mainFrame->CreateToolBar();
		MainToolBarFunc(m_mainFrame->GetToolBar());

		//   m_mainFrame->Centre(wxBOTH);
		m_mainFrame->Show(TRUE);

		SetTopWindow(m_mainFrame);
		//	    RedirectIOLogWindow();

		RedirectIOToConsole();
	}
#endif

	InitFirst();
	//run parallel python
	if (!mpi_py_script.empty()) ExecuteScriptFromFile(mpi_py_script.c_str());

	if (mpi_driver->myrank != 0) return true;

#if !defined(HA_NOGUI)
	if (gui_mode)
	{
		if (mpi_driver->nprocs == 1)
		{
			//			CreateCommandWindow();
			RedirectIOToConsole();
		}
//	    wxLog* p_log = new wxLogStderr();
//	    wxLog::SetActiveTarget(p_log);

//		cout << " HarlemApp::OnInit() test cout " << std::endl;
//	    printf("  HarlemApp::OnInit() test printf() \n");
//		PrintLog(" HarlemApp::OnInit() test PrintLog() \n");

		HaMolSet* pmset = GetCurMolSet();
		if( pmset != NULL)
		{
			pmset->canvas_wx = m_mainFrame->CreateMolView(pmset);
			pmset->canvas_wx->mol_view->InitialTransform();
			pmset->canvas_wx->mol_view->DefaultRepresentation();
			pmset->RefreshAllViews();
		}
		LoadHaPyGUIModules();
	}
#endif
	InitLast();
	
	//post init
	#ifdef __GNUG__
	  FILE *in;
	  std::string HarlemProfilePy=harlem_home_dir + (std::string)"scripts/harlem_postinit.py";
	  in=fopen(HarlemProfilePy.c_str(),"r");
	  if(in!=NULL)
	  {
		fclose(in);
		ExecuteScriptFromFile(HarlemProfilePy.c_str());
	  }
	#endif
	#ifdef _MSC_VER
	  FILE *in;
	  std::string HarlemProfilePy=harlem_home_dir + (std::string)"scripts\\harlem_postinit.py";
	  in=fopen(HarlemProfilePy.c_str(),"r");
	  if(in!=NULL)
	  {
		fclose(in);
		PrintLog("Loading %s\n",HarlemProfilePy.c_str());
		ExecuteScriptFromFile(HarlemProfilePy.c_str());
	  }
	#endif
	return TRUE;
}
int HarlemAppWX::LoadHaPyGUIModules()
{
#ifdef __GNUG__
	FILE *in;
	std::string HarlemProfilePy=harlem_home_dir + (std::string)"scripts/hapygui_init.py";
	in=fopen(HarlemProfilePy.c_str(),"r");
	if(in!=NULL)
	{
		fclose(in);
		return ExecuteScriptFromFile(HarlemProfilePy.c_str());
	}
#endif
#ifdef _MSC_VER
#if PY_VERSION_HEX >= 0x03000000
	// Initialize Harlem's wxPython modules
	PyRun_SimpleString(
		"try:\n"
		"    print('Initiating wxPython')\n"
		"    import wx\n"
		"    import hapygui_init\n"
		"except Exception as e :\n"
		"    print('Can not import hapygui_init module!')\n"
		"    print(str(e))\n"
		"    import traceback\n"
		"    traceback.print_exc()\n"
	);
#endif
#endif
	return 0;
}
int HarlemAppWX::OnExit()
{
	PrintLog("HarlemAppWX::OnExit() end \n");
	return TRUE;
}

#if !defined(HA_NOGUI)

#ifdef WXUSINGDLL
class HaLogFrame;
#else
class WXDLLIMPEXP_FWD_CORE HaLogFrame;
#endif

#ifdef WXUSINGDLL
class HaLogWindow : public wxLogPassThrough
#else
class WXDLLEXPORT HaLogWindow : public wxLogPassThrough
#endif
{
public:
    HaLogWindow(wxWindow *pParent);         // the parent frame (can be NULL)

    virtual ~HaLogWindow();

    // window operations
        // show/hide the log window
    void Show(bool bShow = true);
        // retrieve the pointer to the frame
    wxFrame *GetFrame() const;

    // overridables
        // called immediately after the log frame creation allowing for
        // any extra initializations
    virtual void OnFrameCreate(wxFrame *frame);
        // called if the user closes the window interactively, will not be
        // called if it is destroyed for another reason (such as when program
        // exits) - return true from here to allow the frame to close, false
        // to prevent this from happening
    virtual bool OnFrameClose(wxFrame *frame);
        // called right before the log frame is going to be deleted: will
        // always be called unlike OnFrameClose()
    virtual void OnFrameDelete(wxFrame *frame);

protected:
    virtual void DoLog(wxLogLevel level, const wxChar *szString, time_t t);
    virtual void DoLogString(const wxChar *szString, time_t t);

private:
    HaLogFrame *m_pLogFrame;   //!< the log frame
	std::streambuf *m_sbufOld; //!< cout stream buffer before redirection 

    DECLARE_NO_COPY_CLASS(HaLogWindow)
};

class HaLogFrame : public wxFrame
{
public:
    HaLogFrame(wxWindow *pParent, HaLogWindow *log, const wxString& szTitle);
    virtual ~HaLogFrame();

    // menu callbacks
    void OnClose(wxCommandEvent& event);
    void OnCloseWindow(wxCloseEvent& event);
    void OnSave (wxCommandEvent& event);
    void OnClear(wxCommandEvent& event);

    // accessors
    wxTextCtrl *TextCtrl() const { return m_pTextCtrl; }

private:
    // use standard ids for our commands!
    enum
    {
        Menu_Close = wxID_CLOSE,
        Menu_Save  = wxID_SAVE,
        Menu_Clear = wxID_CLEAR
    };

    // common part of OnClose() and OnCloseWindow()
    void DoClose();

    wxTextCtrl  *m_pTextCtrl;
    HaLogWindow *m_log;

    DECLARE_EVENT_TABLE()
    DECLARE_NO_COPY_CLASS(HaLogFrame)
};

BEGIN_EVENT_TABLE(HaLogFrame, wxFrame)
    // wxLogWindow menu events
    EVT_MENU(Menu_Close, HaLogFrame::OnClose)
    EVT_MENU(Menu_Save,  HaLogFrame::OnSave)
    EVT_MENU(Menu_Clear, HaLogFrame::OnClear)

    EVT_CLOSE(HaLogFrame::OnCloseWindow)
END_EVENT_TABLE()

HaLogFrame::HaLogFrame(wxWindow *pParent, HaLogWindow *log, const wxString& szTitle)
          : wxFrame(pParent, wxID_ANY, szTitle)
{
    m_log = log;

    m_pTextCtrl = new wxTextCtrl(this, wxID_ANY, wxEmptyString, wxDefaultPosition,
            wxDefaultSize,
            wxTE_MULTILINE  |
            wxHSCROLL       |
            // needed for Win32 to avoid 65Kb limit but it doesn't work well
            // when using RichEdit 2.0 which we always do in the Unicode build
            wxTE_RICH       |
            wxTE_READONLY);

    // create menu
    wxMenuBar *pMenuBar = new wxMenuBar;
    wxMenu *pMenu = new wxMenu;
    pMenu->Append(Menu_Save,  _("&Save..."), _("Save log contents to file"));
    pMenu->Append(Menu_Clear, _("C&lear"), _("Clear the log contents"));
    pMenu->AppendSeparator();
    pMenu->Append(Menu_Close, _("&Close"), _("Close this window"));
    pMenuBar->Append(pMenu, _("&Log"));
    SetMenuBar(pMenuBar);

    // status bar for menu prompts
    CreateStatusBar();

    m_log->OnFrameCreate(this);
}

void HaLogFrame::DoClose()
{
    if ( m_log->OnFrameClose(this) )
    {
        // instead of closing just hide the window to be able to Show() it
        // later
        Show(false);
    }
}

void HaLogFrame::OnClose(wxCommandEvent& WXUNUSED(event))
{
    DoClose();
}

void HaLogFrame::OnCloseWindow(wxCloseEvent& WXUNUSED(event))
{
    DoClose();
}

static int OpenLogFile(wxFile& file, wxString *pFilename, wxWindow *parent)
{
    // get the file name
    // -----------------
    wxString filename = wxSaveFileSelector(wxT("log"), wxT("txt"), wxT("log.txt"), parent);
    if ( !filename ) {
        // cancelled
        return -1;
    }

    // open file
    // ---------
    bool bOk wxDUMMY_INITIALIZE(false);
    if ( wxFile::Exists(filename) ) {
        bool bAppend = false;
        wxString strMsg;
        strMsg.Printf(_("Append log to file '%s' (choosing [No] will overwrite it)?"),
                      filename.c_str());
        switch ( wxMessageBox(strMsg, _("Question"),
                              wxICON_QUESTION | wxYES_NO | wxCANCEL) ) {
            case wxYES:
                bAppend = true;
                break;

            case wxNO:
                bAppend = false;
                break;

            case wxCANCEL:
                return -1;

            default:
                wxFAIL_MSG(_("invalid message box return value"));
        }

        if ( bAppend ) {
            bOk = file.Open(filename, wxFile::write_append);
        }
        else {
            bOk = file.Create(filename, true /* overwrite */);
        }
    }
    else {
        bOk = file.Create(filename);
    }

    if ( pFilename )
        *pFilename = filename;

    return bOk;
}


void HaLogFrame::OnSave(wxCommandEvent& WXUNUSED(event))
{
    wxString filename;
    wxFile file;
    int rc = OpenLogFile(file, &filename, this);
    if ( rc == -1 )
    {
        // cancelled
        return;
    }

    bool bOk = rc != 0;

    // retrieve text and save it
    // -------------------------
    int nLines = m_pTextCtrl->GetNumberOfLines();
    for ( int nLine = 0; bOk && nLine < nLines; nLine++ ) {
        bOk = file.Write(m_pTextCtrl->GetLineText(nLine) +
                         wxTextFile::GetEOL());
    }

    if ( bOk )
        bOk = file.Close();

    if ( !bOk ) {
        wxLogError(_("Can't save log contents to file."));
    }
    else {
        wxLogStatus(this, _("Log saved to the file '%s'."), filename.c_str());
    }
}

void HaLogFrame::OnClear(wxCommandEvent& WXUNUSED(event))
{
    m_pTextCtrl->Clear();
}

HaLogFrame::~HaLogFrame()
{
    m_log->OnFrameDelete(this);
}

// HaLogWindow
// -----------

HaLogWindow::HaLogWindow(wxWindow *pParent)
{
    PassMessages(true);
	wxString title("HARLEM CONSOLE");
    m_pLogFrame = new HaLogFrame(pParent, this, title);

    m_pLogFrame->Show();
	m_sbufOld = cout.rdbuf();
    cout.rdbuf(( std::streambuf *)m_pLogFrame->TextCtrl());
	cout << "HaLogWindow::HaLogWindow()" << endl;
	cout.flush();
}

HaLogWindow::~HaLogWindow()
{
	cout.rdbuf(m_sbufOld);
    delete m_pLogFrame;
}

void HaLogWindow::Show(bool bShow)
{
    m_pLogFrame->Show(bShow);
}

void HaLogWindow::DoLog(wxLogLevel level, const wxChar *szString, time_t t)
{
    // first let the previous logger show it
	wxLogRecordInfo info(__FILE__, __LINE__, "DoLog","HaLogWindow");
	info.timestamp = t;
    wxLogPassThrough::DoLogRecord(level, szString, info);

    if ( m_pLogFrame ) {
        switch ( level ) {
            case wxLOG_Status:
                // by default, these messages are ignored by wxLog, so process
                // them ourselves
                if ( !wxIsEmpty(szString) )
                {
                    wxString str;
                    str << _("Status: ") << szString;
                    DoLogString(str.wchar_str(), t);
                }
                break;

                // don't put trace messages in the text window for 2 reasons:
                // 1) there are too many of them
                // 2) they may provoke other trace messages thus sending a program
                //    into an infinite loop
            case wxLOG_Trace:
                break;

            default:
                // and this will format it nicely and call our DoLogString()
                wxLog::DoLogRecord(level, szString, info);
        }
    }
}

void HaLogWindow::DoLogString(const wxChar *szString, time_t WXUNUSED(t))
{
    // put the text into our window
    wxTextCtrl *pText = m_pLogFrame->TextCtrl();

    // remove selection (WriteText is in fact ReplaceSelection)
#ifdef _MSC_VER
    wxTextPos nLen = pText->GetLastPosition();
    pText->SetSelection(nLen, nLen);
#endif // Windows

    wxString msg;
    TimeStamp(&msg);
    msg << szString << wxT('\n');

    pText->AppendText(msg);
}

wxFrame *HaLogWindow::GetFrame() const
{
    return m_pLogFrame;
}

void HaLogWindow::OnFrameCreate(wxFrame * WXUNUSED(frame))
{
}

bool HaLogWindow::OnFrameClose(wxFrame * WXUNUSED(frame))
{
    return false;
}

void HaLogWindow::OnFrameDelete(wxFrame * WXUNUSED(frame))
{
    m_pLogFrame = (HaLogFrame *)NULL;
}

#endif

int HarlemApp::RedirectIOLogWindow()
{
	HarlemAppWX* pApp_wx = (HarlemAppWX*) this;
#if !defined(HA_NOGUI)
	wxLog* p_log = new HaLogWindow(pApp_wx->GetMainFrame());
#endif
	return TRUE;
}

int HarlemAppWX::OnRun()
{
//	PrintLog("HarlemAppWX::OnRun() pt 1 \n");
//	PrintLog( " nprocs = %d \n", mpi_driver->nprocs );
//	PrintLog( " myrank = %d \n", mpi_driver->myrank );
	if( mpi_driver->myrank == 0 )
	{
#if defined(HA_NOGUI)
		return TRUE;
#endif
		return wxApp::OnRun();
	}
	else
	{
//		RedirectIOLogFile();   // FOR MPI_DEBUG
		return mpi_driver->Listen();
	}
	return FALSE;
}

HarlemAppWX* GetAppWX()
{
	return (HarlemAppWX*)pApp;
}
