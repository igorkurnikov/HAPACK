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
IMPLEMENT_APP_NO_MAIN(HarlemAppWX)
#endif

#if defined(_MSC_VER)
extern "C" __declspec(dllexport) int __stdcall start_harlemappwx(int argc, char **argv)
{
	return wxEntry(argc, argv);
}
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
	//if( wxApp::ProcessEvent(event)) return true;

	//MolSet* pmset = GetCurMolSet();

	//if(pmset != NULL) if( pmset->p_evt_h->ProcessEvent(event) ) return true; 

	//int nm = molset_vec.size();
	//int i;
	//for( i = 0; i < nm; i++)
	//{
	//	pmset = (MolSet*) molset_vec[i];
	//	if( pmset == GetCurMolSet() ) continue;
	//	if( pmset->p_evt_h->ProcessEvent(event) ) return true; 
	//}
	return false;
}

#if !defined(HA_NOGUI)
void HarlemAppWX::OnIdle(wxIdleEvent& event)
{
//	int ret = PyRun_InteractiveOne(stdin, "<stdin>");  //Igor Python incorporation...
  
	//fprintf(stdout,"MolSet::OnIdle >\n");
  //PyRun_InteractiveOne(stdin,"none");
  //fprintf(stdout,"MolSet::OnIdle <\n");
      
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
    fprintf(stdout,"MolSet::OnIdle [%s]\n",PyCmd);
    
  }
  else
  {
    fprintf(stdout,"MolSet::OnIdle No input\n");
    clearerr(stdin);
  }*/
}
#endif

void HarlemAppWX::OnExitApp(wxCommandEvent& event)
{
	delete this;
}

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
		}
//	    wxLog* p_log = new wxLogStderr();
//	    wxLog::SetActiveTarget(p_log);

//		cout << " HarlemApp::OnInit() test cout " << std::endl;
//	    printf("  HarlemApp::OnInit() test printf() \n");
//		PrintLog(" HarlemApp::OnInit() test PrintLog() \n");

		MolSet* pmset = GetCurMolSet();
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
	return 0;
}
int HarlemAppWX::OnExit()
{
	PrintLog("HarlemAppWX::OnExit() end \n");
#if !defined(HA_NOGUI)
#if(_MSC_VER)
	//FreeConsole();
#endif
#endif
	return TRUE;
}

int HarlemAppWX::OnRun()
{
    PrintLog("HarlemAppWX::OnRun() pt 1 \n");
    PrintLog( " nprocs = %d \n", mpi_driver->nprocs );
    PrintLog( " myrank = %d \n", mpi_driver->myrank );
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


