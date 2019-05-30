/*! \file harlemapp.cpp
 
    Application class of HARLEM.

    \author Igor Kurnikov 
    \date 1999-

*/

#define HARLEMAPP_CPP

#include <mpi.h> 
#include "haconst.h"
#include <wchar.h>

#if !defined(HARLEM_PYTHON_NO)
#include <Python.h>
#endif

#include <boost/algorithm/string.hpp>

#if defined(_MSC_VER)
//#include <afx.h>
#endif

#include "wx/wxprec.h"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include <wx/stdpaths.h>
#include <wx/process.h>
#include <wx/filename.h>


#if !defined(_MSC_VER)
#include <unistd.h>
#endif


#include "hampi.h"
#include "harlemapp.h"
#include "harlemapp_wx.h"
#include "hamolview.h"
#include "hamolecule.h"
#include "hastring.h"
#include "haio.h"
#include "tokens.h"
#include "abstree.h"
#include "hamolset.h"
//#include "hamatdb.h"
#include "etcoupl.h"
#include "haqchem.h"
#include "hagaussian.h"
#include "hadalton.h"
#include "haresdb.h"
#include "mm_force_field.h"

//#include "memchk.h"

#if defined(linux)
int __argc_save;
char** __argv_save;
#endif

#if defined(HA_NOGUI)
int PrintMessage(const char* str)
{
	printf("%s\n",str);
	return TRUE;
}

#endif

void StartHarlemApp()
{
	HarlemApp *m_HarlemApp = new HarlemApp();
	m_HarlemApp->InitFirst();
	m_HarlemApp->InitLast();
}


HarlemApp::HarlemApp()
{
#if defined(HA_NOGUI)
	gui_mode = FALSE;
#else
    gui_mode = TRUE;
#endif
	cmd_prompt_mode = TRUE;
	mpi_py_script = "";

	argc_loc = 0;
	argv_loc = NULL;

	FileFormat = FormatPDB;
	only_rasmol_command = FALSE;

	FormatOpt["alchemy"]   =  FormatAlchemy;
	FormatOpt["biosym"]    =  FormatBiosym;
	FormatOpt["cif"]       =  FormatCIF;
	FormatOpt["charmm"]    =  FormatCharmm;
	FormatOpt["fdat"]      =  FormatFDAT;
	FormatOpt["gaussian"]  =  FormatGaussian;
	FormatOpt["macromodel"]=  FormatMacroMod;
	FormatOpt["mdl"]       =  FormatMDL;
	FormatOpt["mmdb"]      =  FormatMMDB;
	FormatOpt["mol2"]      =  FormatMol2;
	FormatOpt["mopac"]     =  FormatMOPAC;
	FormatOpt["nmrpdb"]    =  FormatNMRPDB;
	FormatOpt["pdb"]       =  FormatPDB; 
	FormatOpt["shelx"]     =  FormatSHELX; 
	FormatOpt["xyz"]       =  FormatXYZ; 
	FormatOpt["rwf"]       =  FormatRWF; 
	FormatOpt["hlm"]       =  FormatHarlem; 
	FormatOpt["in"]        =  FormatAmberPrep; 
	FormatOpt["top"]       =  FormatAmberTop;

	HaAtom::FillStdAtomTypes();
	HaResidue::InitStdResNames();
	HaResidue::InitResSynonym();

	mpi_driver    = NULL;
	python_thread = NULL;
	file_log            = NULL;
}

HarlemApp::~HarlemApp()
{
//   PrintLog("HarlemApp::~HarlemApp() pt 1 \n");
   int nm = molset_vec.size();
   if( nm > 0 )
   {
	  VecPtr molset_vec_axx = molset_vec;
	  int i;
	  for( i = 0; i < nm; i++)
	  {
		 HaMolSet* pmset = (HaMolSet*) molset_vec_axx[i];
		 this->DeleteMolSet(pmset);
	  }
   }
   if( file_log  != NULL)
   {
	   fclose(file_log);
	   file_log = NULL;
   }
   FinalizeXML();
//   if(mpi_driver->myrank == 0)
//   {
//		if( mpi_driver->nprocs > 1)
//		{
//			PrintLog("HarlemApp::~HarlemApp() pt 2 \n");
//			pApp->mpi_driver->SendKillAppMsgAllProc();
//			PrintLog("HarlemApp::~HarlemApp() pt 3 \n");
//		}
// }
   if( mpi_driver->nprocs > 1)
   {
	   MPI_Abort(MPI_COMM_WORLD, -1);
   }
   delete mpi_driver;
}


#if defined(_MSC_VER)

#include <io.h>
#include <fcntl.h>

#endif

int HarlemApp::InitFirst()
{
	pApp = this;

//	memory_manager = new Memory_Manager(100000);
//	memory_manager->on();
//	memory_manager->print(cout);

// Set HARLEM HOME directory
	wxString harlem_home_str;
	bool exist = wxGetEnv("HARLEM_HOME",&harlem_home_str);
    harlem_home_str.Trim();
#if(_MSC_VER)
//	auto path_obj = wxStandardPaths::Get();
	wxString exe_dir = wxStandardPaths::Get().GetDataDir();
	if( exist && !harlem_home_str.IsEmpty() )
	{
		harlem_home_dir = harlem_home_str;  
	}
	else
	{
//		harlem_home_dir = "C:\\HARLEM\\";
		harlem_home_dir = exe_dir.ToStdString() + "\\";
	}
	res_db_dir = harlem_home_dir + "residues_db\\";
	word_editor = harlem_home_dir + "scite.exe";
	html_browser = "c:\\Program Files\\Internet Explorer\\iexplore.exe";
#else
	if( exist && !harlem_home_str.IsEmpty() )
	{	
		if(harlem_home_str.Last() != '/') harlem_home_str += '/';
		harlem_home_dir = harlem_home_str;
	}
	else
	{
		harlem_home_dir = "/usr/local/lib/harlem/";
	}
	res_db_dir = harlem_home_dir + "residues_db/";
	word_editor = "scite";
	html_browser = "firefox";
#endif
	
	//char* henv= getenv("HARLEM_HOME");
	//if(henv != NULL) harlem_home_dir = henv;
	
	manual_main_page =  harlem_home_dir + "doc/prog_man_html/index.html";
//  manual_main_page =  harlem_home_dir + "manual/harlem_manual.chm";
	
	InitParallel();
	InitXML();
	ProcessOptions();

	if( mpi_driver->myrank != 0 && mpi_py_script.empty()) return TRUE;
	
//    PrintLog(" HarlemApp::InitFirst() pt 3 \n");

	InitCommand();

//	PrintLog("HARLEM: HAmiltonians to Research LargE Molecules \n");   
//	PrintLog("Igor Kurnikov, Nikolay Simakov, Kirill Speransky, Arvind Ramanathan, \n"); 
//	PrintLog("Maria Kurnikova 1997 - 2014 \n" );
//	PrintLog("Graphical Interface Based on RASMOL 2.6 of Roger Sayle \n");
//	PrintLog("and wxWidgets library \n");
//	PrintLog("Command line/scripting interface is based on PYTHON 2.7 \n");
//	PrintLog("Also linked to VFLIB, TINYXML, RAPIDXML, LAPACK, BLAS and other libraries \n");
//	PrintLog("Build date -  %s \n",__DATE__);               
//	PrintLog("Harlem Home Dir %s\n",harlem_home_dir.c_str());  

	InitRemoteComp();

	if( !finp_name.empty() ) // if filename has been specified on the command line: load the file
	{
		HaMolSet* pmset = new HaMolSet();
		wxString fname_str = finp_name.c_str();
		PrintLog(" input file name %s \n", fname_str.ToStdString().c_str());
		wxFileName fname_obj(fname_str);
		wxString path = fname_obj.GetPath();
		PrintLog(" path of the input file %s \n", path.ToStdString().c_str());
		path = path.Strip(wxString::both);
			 
		wxString cur_dir = ::wxGetCwd();
		PrintLog(" Current Working Directory %s \n",cur_dir.ToStdString().c_str());
		if(!path.IsEmpty())
		{
			wxFileName::SetCwd(path);
		}
		PrintLog(" Current Working Directory %s \n",cur_dir.ToStdString().c_str());
			
		std::string exten = harlem::GetExtFromFileName(finp_name);
		
		boost::to_upper(exten);
		if(exten == "HLM") FileFormat = FormatHarlem;
		else if( exten == "PDB" || exten == "ENT")
			FileFormat = FormatPDB;
		else if( exten == "XYZ" )
			FileFormat = FormatXYZ;
		else if( exten == "HIN" )
			FileFormat = FormatHIN;

		int result = pmset->FetchFile(FileFormat,finp_name.c_str());
		if(!result)
		{
			PrintLog(" HarlemAppWX::OnInit() ");
			PrintLog(" Error loading file %s \n",finp_name.c_str());
		}
	}
	return TRUE;
}

int HarlemApp::InitParallel()
{
	wxString cur_dir = ::wxGetCwd();	
	mpi_driver = new HaMPI();
	//Current Directory is changed by MPI_Init() - changing back...
	wxSetWorkingDirectory(cur_dir);
	return TRUE;
}

int HarlemApp::InitXML()
{
//  Commented initiation Xerces-c XML processing for now
//	try
//  {
//        XMLPlatformUtils::Initialize();
//    }
//    catch(const XMLException& ex)
//    {
//        PrintLog("Error in HarlemApp::InitXML() \n");
//		PrintLog("During Xerces-c Initialization \n");
//        PrintLog("  Exception message: %s \n",ex.getMessage());
//        return FALSE;
//    }
	return TRUE;
}

int HarlemApp::FinalizeXML()
{
// Commented initiation Xerces-c XML processing for now
// XMLPlatformUtils::Terminate();
	return TRUE;
}

int HarlemApp::InitRemoteComp()
{
	std::string comp_acc_file_name = harlem_home_dir;
	comp_acc_file_name +=  "/local/available_computers";
	
	FILE* comp_acc_file = fopen(comp_acc_file_name.c_str(),"r");

	if(comp_acc_file != NULL)
	{
		char buf[256];
		char name1[256];
		char name2[256];
		char name3[256];
		for(;;)
		{
			char* cres = fgets(buf,255,comp_acc_file);
			if(cres == NULL) break;
			sscanf(buf,"%s %s %s",name1,name2,name3);
			ComputerAccount new_acc;
			new_acc.acc_id = name1;
			new_acc.login_str = name2;
			new_acc.interm_acc_ID = name3;
			comp_accounts.push_back(new_acc);
		}
	}
	else
	{
//		PrintLog(" Computer Accounts File does not exist \n");
	}
//	int i;
//	int na = comp_accounts.size();
//	PrintLog("\n Available Computer Accounts: \n");
//	for(i = 0; i < na ; i++)
//	{
//		PrintLog(" %s  %s  %s \n",comp_accounts[i].acc_id.c_str(),comp_accounts[i].login_str.c_str(),
//			                      comp_accounts[i].interm_acc_ID.c_str());
//	}
	return TRUE;
}


int HarlemApp::RasMolCmd(const char* cmd)
{
	cmd_pr.SetCmdLine(cmd);
	only_rasmol_command = TRUE;
	this->ExecuteCommand();
	only_rasmol_command = FALSE;
	return TRUE;
}

int HarlemApp::ExecRasMolScript(const char* file_name)
{
	char buf[20000];
	ifstream sfile(file_name);
	for(;;)
	{
		sfile.getline(buf,19999);
		if(sfile.fail()) break;
		RasMolCmd(buf);
	}
	return TRUE;
}


int HarlemApp::ExecuteCommand()
{
    int option;
    int i;

	std::string fname;
	int ip;
	PrintLog("Executing Command: %s \n", cmd_pr.GetCmdLine() ); 

#if !defined(HARLEM_PYTHON_NO)
	if(!only_rasmol_command)
	{
		PyGILState_STATE gstate;
		gstate = PyGILState_Ensure();
		int ires = PyRun_SimpleString(cmd_pr.GetCmdLine());
		PyGILState_Release(gstate);
		
		if( ires == 0 )
			return TRUE;
		else
		{
//			  PyErr_Print();
			PrintLog("RASMOL type command \n");
		}
	}
#endif

	only_rasmol_command = FALSE;
	cmd_pr.ResetCursorPosition();
	
    if( !cmd_pr.FetchToken() )
    {   
		cmd_pr.cursor_pos = 0;
        return( FALSE );
    }

    switch( cmd_pr.CurToken )
    {   
	case(SelectTok): 
	case(RestrictTok):
	case(EditTok):
	case(LoadTok):  
	case(RenumTok): 
    case(StructureTok):
    case(OverlapMolTok):
	case(AlignOverlapMolTok):
	case(DefineTok):
	case(RefreshTok): 
	case(ConnectTok):
	case(GroupsTok):
	case(ETTok):
	case(QChemTok):
	case(InterMolTok):
	case(MolMechTok):
	case(MoleculeTok):
	case(ShowTok):
		if( GetCurMolSet() != NULL) 
			GetCurMolSet()->ExecuteCommand(cmd_pr);
		break;

	case(ColourTok):  
	case(WireframeTok):
	case(BackboneTok):
    case(CPKTok):
    case(SpacefillTok):
    case(DashTok): 
    case(SSBondTok):  
    case(HBondTok):   
    case(RibbonTok): 
    case(StrandsTok): 
    case(TraceTok):
    case(CartoonTok): 
    case(DotsTok):  
    case(MonitorTok): 
    case(SlabTok): 
    case(ZoomTok):  
    case(RotateTok):  
    case(TranslateTok):
    case(StereoTok):
    case(CentreTok):  
    case(ResizeTok):  
    case(ResetTok):
    case(SetTok):     
    case(LabelTok):
	case(BackgroundTok):
	case(PrintTok):
    case(ClipboardTok):
		if(CurMolView != NULL)
			CurMolView->ExecuteCommand(cmd_pr);
			break;

    case(EchoTok):    
		cmd_pr.FetchToken();
		PrintLog("\n");
			
		if( cmd_pr.CurToken == StringTok )
		{   
			PrintLog(cmd_pr.TokenIdent.c_str());
		} 
		else if( cmd_pr.CurToken )
		{
			PrintLog( cmd_pr.GetStartPosSubstr() );
		}
		PrintLog("\n");
		cmd_pr.CurToken = 0;
		break;

        case(SourceTok):
        case(ScriptTok):  
			cmd_pr.FetchToken();	  
			if( !cmd_pr.CurToken )
			{   
				PrintLog("Filename string expected\n");
				break;
			} 
			else if( cmd_pr.CurToken==StringTok )
			{      
				fname = cmd_pr.TokenIdent;
			} 
			else 
			{
				fname = cmd_pr.GetStartPosSubstr();
			}
			boost::trim(fname);
			ip = fname.find(' ');
			if(ip != -1) fname = fname.substr(ip);
			cmd_pr.CurToken = 0;
			
			ExecuteScriptFromFile(fname.c_str());
			break;

		break;

        case(ExitTok):    return( ExitTok );
        case(QuitTok):    return( QuitTok );
        default:          cmd_pr.CommandError("Unrecognised command");
                          break;
    }

    if( cmd_pr.CurToken )
	{
        if( cmd_pr.FetchToken() ) cmd_pr.CommandError("Warning: Ignoring rest of command");
	}
    return( False );
}


extern "C" {
#if PY_VERSION_HEX >= 0x03000000
	extern PyObject* PyInit__molset();
	extern PyObject* PyInit__halib();
	//extern PyObject* PyInit__llpnps();
#else
	extern void init_molset();
	extern void init_halib();
	//extern void init_llpnps();
#endif

}

int HarlemApp::Python_AppInit()
//! PYTHON initialization: 
{	
//Do not need it any more
#if 0
#if !defined(HARLEM_PYTHON_NO)
	int ires;

	char *command = NULL;
	char *filename = NULL;
	int stdin_is_interactive;
	FILE *fp = stdin;

	stdin_is_interactive = Py_FdIsInteractive(stdin, (char *)0);

//	PrintLog(" Python_AppInit pt 1 \n");

#ifdef _MSC_VER
		_setmode(fileno(stdin),  _O_BINARY);
		_setmode(fileno(stdout), _O_BINARY);
#endif

#ifdef HAVE_SETVBUF
		setvbuf(stdin,  (char *)NULL, _IONBF, BUFSIZ);
		setvbuf(stdout, (char *)NULL, _IONBF, BUFSIZ);
		setvbuf(stderr, (char *)NULL, _IONBF, BUFSIZ);
#else /* !HAVE_SETVBUF */
		setbuf(stdin,  (char *)NULL);
		setbuf(stdout, (char *)NULL);
		setbuf(stderr, (char *)NULL);
#endif /* !HAVE_SETVBUF */

//	PrintLog(" Python_AppInit pt 3 \n");
#if PY_VERSION_HEX >= 0x03000000
	wchar_t *prog_name = Py_DecodeLocale("HARLEM", NULL);

// Load Harlem Python extension modules  PYTHON 3
#if PY_VERSION_HEX >= 0x03000000
	PyImport_AppendInittab("_molset", &PyInit__molset );
	//PyImport_AppendInittab("_halib",  &PyInit__halib );
	//PyImport_AppendInittab("_llpnps", &PyInit__llpnps );
#endif

	Py_SetProgramName(prog_name);
//	PyMem_RawFree(prog_name);
#else
	Py_SetProgramName("HARLEM");
#endif

#ifdef _DEBUG
	char* env_p;
	if (env_p = std::getenv("PATH"))
		std::cout << "Your PATH is: " << env_p << '\n';
	if (env_p = std::getenv("PYTHONPATH"))
		std::cout << "Your PYTHONPATH is: " << env_p << '\n';
	std::cout << '\n';
#endif
	/* Initialize the Python interpreter.  Required. */
	Py_Initialize();	
//	PrintLog(" Python_AppInit pt 4 \n");

#if PY_VERSION_HEX >= 0x03000000

#else
	if (argc_loc > 1)
		PySys_SetArgv(argc_loc - 1, argv_loc + 1);
#endif

	PyObject *v;

	if( isatty(fileno(stdin))) 
	{
		v = PyImport_ImportModule("readline");
		if (v == NULL)
			PyErr_Clear();
		else
			Py_DECREF(v);
	}
	
// Load Harlem Python extension modules  PYTHON 2
#if PY_VERSION_HEX >= 0x03000000
//	PyImport_AppendInittab( "_molset", &PyInit__molset );
//	PyObject* halib_py =  PyInit__halib();
//	PyObject* llpnps = PyInit__llpnps();
#else
	init_molset();
	init_halib();
	//init_llpnps();
#endif
	//PrintLog(" Python_AppInit pt 5 \n")
	ires = PyRun_SimpleString(
		"import os\n"
		"import sys\n"
	);
	if( ires != 0 ) PrintLog("Error in HarlemApp::Python_AppInit() loading sys module \n");
#if defined(_MSC_VER)
	std::string set_harlem_home_dir="HARLEM_HOME_DIR = r'" + harlem_home_dir + "\\'[:-2]";
	ires = PyRun_SimpleString(set_harlem_home_dir.c_str());
	ires = PyRun_SimpleString(
		"sys.path.insert(1, HARLEM_HOME_DIR)\n"
		"sys.path.insert(1, os.path.join(HARLEM_HOME_DIR, 'residues_db'))\n"
	);

	// add path of harlem python modules
	ires = PyRun_SimpleString(
		"HAPACK_DIR = os.path.dirname(os.path.dirname(HARLEM_HOME_DIR))\n"
		"if HAPACK_DIR == 'C:\\MYPROG\\HAPACK':\n"
		"    print('Running from VS build directory, adding path to harlem python modules from source code.\\n')\n"
		"    sys.path.insert(1, os.path.join(HAPACK_DIR, 'HARLEM', 'scripts'))\n"
		"    sys.path.insert(1, os.path.join(HAPACK_DIR, 'PNPS', 'src'))\n"
		"    sys.path.insert(1, os.path.join(HAPACK_DIR, 'PNPS'))\n"
	);
#else
	std::string HaScriptDir=(std::string)"sys.path.insert(1,\"" + harlem_home_dir + (std::string)"scripts\")";
        //ires = PyRun_SimpleString(HaScriptDir.c_str());
	std::string HaDBDir=(std::string)"sys.path.insert(2,\"" + harlem_home_dir + (std::string)"residues_db\")";
	ires = PyRun_SimpleString(HaDBDir.c_str());
#endif

#ifdef _DEBUG
	ires = PyRun_SimpleString("print('sys.path = ' + str(sys.path))");
#endif

	ires = PyRun_SimpleString("import harlempy.halib");
	if( ires != 0 ) PrintLog("Error in HarlemApp::Python_AppInit() loading halib module \n");
	ires = PyRun_SimpleString("from harlempy.halib import *");

	ires = PyRun_SimpleString("import harlempy.molset");
	if( ires != 0 ) PrintLog("Error in HarlemApp::Python_AppInit() loading molset module \n");
	ires = PyRun_SimpleString("from harlempy.molset import *");

//	if( !ires ) PrintLog("Error in HarlemApp::Python_AppInit() loading harlem profile script \n");
	ires = PyRun_SimpleString("pApp = cvar.pApp");
  
	v = PySys_GetObject("ps1");
	if (v == NULL) {
#if PY_VERSION_HEX >= 0x03000000
		PySys_SetObject("ps1", v = PyUnicode_FromString(">>> "));
#else
		PySys_SetObject("ps1", v = PyString_FromString(">>> "));
#endif
		Py_XDECREF(v);
	}
	
	v = PySys_GetObject("ps2");
	if (v == NULL) {
#if PY_VERSION_HEX >= 0x03000000
		PySys_SetObject("ps2", v = PyUnicode_FromString("... "));
#else
		PySys_SetObject("ps2", v = PyString_FromString("... "));
#endif
		Py_XDECREF(v);
	}

//	Py_Finalize();
#endif
#endif
	return TRUE;
}

int HarlemApp::CreateCommandWindow()
{
	 int ires = PyRun_SimpleString("import wx");  if( ires ) return FALSE;
	 ires = PyRun_SimpleString("import wx.py");  if( ires ) return FALSE;
	 ires = PyRun_SimpleString("frm = wx.Frame(None, -1, \"HARLEM CONSOLE\")");   if( ires ) return FALSE;
	 ires = PyRun_SimpleString("sh = wx.py.shell.Shell(frm)");   if( ires ) return FALSE;
	 ires = PyRun_SimpleString("frm.Show()");   if( ires ) return FALSE;
//	 ires = PyRun_SimpleString("sh.redirectStderr()");   if( ires ) return FALSE;
//	 ires = PyRun_SimpleString("sh.redirectStdout()");   if( ires ) return FALSE;
//	 printf("hello from HarlemApp::CreateCommandWindow() \n" );
	 return TRUE;
}

int HarlemApp::RedirectIOLogFile(const std::string& fname_new )
{
	PrintLog("HarlemApp::RedirectIOLogFile() \n");
	std::string fname = boost::trim_copy(fname_new);

	if(fname.empty()) fname = "harlem_proc_" + harlem::ToString(mpi_driver->myrank) + ".log";

	PrintLog(" HarlemApp::RedirectIOLogFile()  fname = %s\n", fname.c_str() ); 
	
	if( file_log != NULL) 
	{
		fclose(file_log); file_log = NULL;
	}

	file_log = fopen(fname.c_str(),"w");
	if(file_log) 
	{
		*stdout = *file_log;
		setvbuf( stdout, NULL, _IONBF, 0 ); // set stdout unbuffered
		*stderr = *file_log;
		setvbuf( stderr, NULL, _IONBF, 0 ); // set stderr unbuffered
		
		ios::sync_with_stdio();

	    wxLog* p_log = new wxLogStderr();
        wxLog::SetActiveTarget(p_log);

		return TRUE;
	}
	return TRUE;
}


void HarlemApp::InitCommand()
{
	cmd_pr.InitKeywords();

// PYTHON initialization
 	if(!Python_AppInit() )
	{
		ErrorInMod("HarlemApp::InitCommand()",
			       "Error to Initialize PYTHON commands processing");
	}
	
}

int HarlemApp::InitLast()
{
	if( mpi_driver->myrank == 0 )
	{	
		LoadInitFile();
		FILE *fp;
		if( !script_name.empty() )
		{
			if( !(fp=fopen(script_name.c_str(),"r")) )
			{
				PrintLog("Error: File %s not found!\n",script_name.c_str());
			}
			else
			{
				fclose(fp);
				ExecuteScriptFromFile(script_name.c_str());
			}
		}
		if( !script_str.empty())
		{
			ExecuteScriptInString(script_str.c_str());
		}
		if( !cmd_prompt_mode )
		{
			PrintLog(" HARLEM EXIT 1 \n");
			delete this;
			exit(0);
		} 
		if( this->cmd_prompt_mode )
		{
			if( this->gui_mode )
			{
//				pApp->python_thread = new PythonThread();
//				((PythonThread*)pApp->python_thread)->Create();
//				((PythonThread*)pApp->python_thread)->Run();
			}
			else
			{
#if !defined(HARLEM_PYTHON_NO)
			    int sts = PyRun_AnyFile(stdin, "<stdin>");
#endif
			}
		}
		// Print wellcome message
		PrintLog("\n");
		PrintLog("HARLEM (HAmiltonians to Research LargE Molecules)\n");
		PrintLog("==============================================================================\n\n");
	}
	return TRUE;
}

static StrIntMap FormatOpt;


int HarlemApp::ProcessOptions()
//! Command line options
//!  -format  (where format = pdb, hlm, etc...) determine format of the file
//!                      to load, set FileFormat value
//! file_name - name of the file to load, set finp_name string
//! -script or script_fname  - load script file  (set script_name string)
//
{
	PyObject* sys_mod = PyImport_ImportModule("sys");
	PyObject* argv_obj = PyObject_GetAttrString(sys_mod, "argv" );
	int check_list = PyList_Check(argv_obj);
	Py_ssize_t argv_size = PyList_Size(argv_obj);
	
	Py_ssize_t i;
		
	for(i=1; i < argv_size; i++)
	{
		PyObject* arg_obj = PyList_GetItem(argv_obj, i);
		const char *arg_str = PyUnicode_AsUTF8(arg_obj);
		std::string option = arg_str;
		boost::trim(option);
		std::string option_next = "";
		if ( option.empty() ) continue;
		if ((i + 1) < argv_size)
		{
			arg_obj = PyList_GetItem(argv_obj, i+1);
			arg_str = PyUnicode_AsUTF8(arg_obj);
			option_next = arg_str;
			boost::trim(option_next);
		}
		
		if( option[0] == '-')
		{
			std::string trunc_opt = option.substr(1);
			boost::to_lower(trunc_opt);

			if( trunc_opt == "nogui" )
			{
				gui_mode = FALSE;
				continue;
			}

			if( trunc_opt == "noprompt" )
			{
				cmd_prompt_mode = FALSE;
				continue;
			}

			if( FormatOpt.count(trunc_opt.c_str()) > 0)
			{
				FileFormat = FormatOpt[trunc_opt];
				i++;
				if( i >= argv_size )
				{
					PrintLog("No Molecular Geometry Input File Name supplied in the command line\n");
					continue;
				}
				finp_name = option_next;
				continue;
			}

			if( trunc_opt == "script" )
			{
				i++;
				if( i >= argv_size )
				{
					PrintLog("No Script Name supplied in the command line\n");
					continue;
				}
				script_name = option_next;
				continue;
			}

			if( trunc_opt == "mpipy" )
			{
				gui_mode = FALSE;
				cmd_prompt_mode = FALSE;
				i++;
				if( i >= argv_size )
				{
					PrintLog("No Script Name supplied in the command line\n");
					continue;
				}
				mpi_py_script = option_next;
				continue;
			}

			if( trunc_opt == "wd" )
			{
				i++;
				if( i >= argv_size )
				{
					PrintLog("No Working Directory Name supplied for -wd option\n");
					continue;
				}
				wxString work_dir = option_next;
				PrintLog("Set working directory to %s \n",work_dir.ToStdString().c_str());
				::wxSetWorkingDirectory(work_dir);
				continue;
			}
		}
		else
		{
			if( !finp_name.empty() )
			{
				PrintLog("HarlemApp::ProcessOptions()", "Warning: Input File Name already set");
			}
			finp_name = option;
		}
	}

//	PrintLog(" ProcessOptions() \n");
//	PrintLog(" Command Line is : %s\n",cmd_line.c_str() ); 

#if 0
		else if( finp_name.empty() )
		{
			int in_quote = FALSE;
			while( *ptr )
			{
				if(*ptr == '"' && in_quote) break;
				if(*ptr == '"') 
				{
					in_quote = TRUE;
					ptr++;
				}
				if(*ptr == ' ' && !in_quote) break;
				finp_name += *ptr;
				ptr++;
			}
		}
#endif
	
	return( True );
}


static int PrefixString(register char* str1, register char* str2 )
{
	while( *str1 == *str2++ )
	if( *str1++ == '\0' )
		return( True );
	return( *str1 == '\0' );
}


void HarlemApp::LoadInitFile()
{
	char *src,*dst;
	FILE *initrc;
	char *fname;

	char buf[256];

	fname = "HARLEM.INI";
	initrc = fopen(fname,"r");
	if( !initrc && (src=(char*)getenv("HOME")) )
	{   
		dst = buf;
	    while( *src )
	       *dst++ = *src++;
	    *dst++ = '\\';

    	src = fname; fname = buf;
	    while( *dst++ = *src++ );
	    initrc = fopen(fname,"r");
    }

    if( initrc )
	{
		fclose(initrc);
		ExecuteScriptFromFile(fname);
	}
}

void HarlemApp::StartWait()
{

}

void HarlemApp::EndWait()
{

}

int HarlemApp::CheckProcIsActive(long proc_id)
{
//	PrintLog(" Process %d ",proc_id);
	bool bres = wxProcess::Exists(proc_id);
	if(bres) 
	{
//		PrintLog(" Running \n");
	}
	else
	{
//		PrintLog(" Does not exist \n");
	}
	if( bres ) return TRUE;
    return FALSE;
}

int HarlemApp::KillProc(long proc_id)
{
	int ires;
	wxKillError error;
//	int ires = wxKill(proc_id,wxSIGTERM,&error);
    ires = wxKill(proc_id,wxSIGKILL,&error);
	return ires;
}

int HarlemApp::SwitchThread()
{
//	SwitchToThread();
	wxThread::Sleep(100);

    return FALSE;
}

int HarlemApp::SleepThread(int ms_delay)
{
   wxThread::Sleep(ms_delay);
   return FALSE;
}

ComputerAccount* HarlemApp::GetAccountByID(const char* acc_id)
{
	ComputerAccount* pacc = NULL;
	int i;
	int na = comp_accounts.size();
	for( i = 0; i < na; i++)
	{
		if(comp_accounts[i].acc_id == acc_id)
		{
			pacc = &comp_accounts[i];
			break;
		}
	}
	return pacc;
}

int HarlemApp::ExecuteRemoteCmd(ComputerAccount* pacc, const char* cmd, StrVec& prog_output, int get_prog_output)
{
	if(pacc == NULL) return FALSE;
	
	StrVec prog_args;

	std::string loc_cmd_str;
	
	prog_args.push_back("ssh");

	if(pacc->interm_acc_ID == "DIRECT")
	{
		prog_args.push_back(pacc->login_str);
		prog_args.push_back(cmd);
	}
	else
	{
		ComputerAccount* interm_pacc = GetAccountByID(pacc->interm_acc_ID.c_str());
		if(interm_pacc == NULL)
		{
			PrintLog(" Error in HarlemApp::ExecuteRemoteCmd() \n");
			PrintLog(" No Intermediate Account with ID %s \n",
				       pacc->interm_acc_ID.c_str());
			return FALSE;
		}
		prog_args.push_back(interm_pacc->login_str);
		std::string cmd_str2;
		cmd_str2 = (std::string)"ssh " + pacc->login_str + " " + (std::string) cmd;
		
		prog_args.push_back(cmd_str2);
	}
	int i; 
	int narg = prog_args.size();
	
	for(i =0; i < narg; i++)
	{
		loc_cmd_str += " " + prog_args[i];
	}
	
//	PrintLog("Executing Remote Command \n %s \n",loc_cmd_str.c_str());
	
	RunExternalProgram(RUN_BACKGROUND,"ssh",prog_args, prog_output,get_prog_output);
	return TRUE;
}

int HarlemApp::ShowAccountsLoad()
{
	int i;
	int na = comp_accounts.size();

	for(i= 0; i < na; i++)
	{
		std::string cmd_str = "ps -eo ";
		cmd_str += '\"';
		cmd_str += "%p%C";
		cmd_str += '\"';
		StrVec prog_output;
		ExecuteRemoteCmd(&comp_accounts[i],cmd_str.c_str(),prog_output,TRUE);
		
		double tot_cpu = 0.0;
		double perc_cpu;
		int pid;

		for(int j=0; j < prog_output.size(); j++)
		{			
			if( i > 0)
			{
				sscanf(prog_output[j].c_str()," %d %lf ",&pid,&perc_cpu);
				tot_cpu += perc_cpu;
			}
		}
		PrintLog(" CPU load for Account %s  is %7.2f  \n",comp_accounts[i].acc_id.c_str(),tot_cpu);
	}

	return TRUE;
}

long HarlemApp::RunExternalProgram(RunMode rmode, const std::string& prog_name, StrVec& prog_args,
							  StrVec& prog_output, int get_prog_output )
{
	int narg = prog_args.size();
	const char** argv_p = (const char**) malloc( (narg + 1)* sizeof( const char* )); 
	
	long i, result;
	std::string cmd_line = prog_name;
	for( i =0; i < narg; i++)
	{
		argv_p[i] = prog_args[i].c_str();
		cmd_line += " " + prog_args[i];
	}

	argv_p[narg] = NULL;
	result = 0;

	long pid = 0;
 
	char buf[256];
	std::string cmd_output;

	int sync_wx = wxEXEC_ASYNC;
   
#if defined(_MSC_VER) 
	if(rmode == RUN_BACKGROUND)
	{
        sync_wx = wxEXEC_ASYNC;
	}
	else if(rmode == RUN_FOREGROUND )
	{
		sync_wx = wxEXEC_SYNC;
	}

	HANDLE hSaveStdout,hChildStdoutRdDup,hChildStdoutWr,hChildStdoutRd;

	if(get_prog_output)
	{
		hSaveStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
				
		SECURITY_ATTRIBUTES saAttr; 
        // Set the bInheritHandle flag so pipe handles are inherited. 
        saAttr.nLength = sizeof(SECURITY_ATTRIBUTES); 
        saAttr.bInheritHandle = TRUE; 
        saAttr.lpSecurityDescriptor = NULL; 
		
	    // Create a pipe for the child process's STDOUT. 
		
		if( !CreatePipe(&hChildStdoutRd, &hChildStdoutWr, &saAttr, 0) ) 
		{
			PrintLog("Stdout pipe creation failed\n"); 
			return FALSE;
		}
		
	    // Set a write handle to the pipe to be STDOUT. 
		
		if (!SetStdHandle(STD_OUTPUT_HANDLE, hChildStdoutWr))
		{
			PrintLog("Redirecting STDOUT failed \n"); 
			return FALSE;
		}
		// Create noninheritable read handle and close the inheritable read 
		// handle. 
		
		int fSuccess = DuplicateHandle(GetCurrentProcess(), hChildStdoutRd,
			GetCurrentProcess(), &hChildStdoutRdDup , 0,
			FALSE,
			DUPLICATE_SAME_ACCESS);
		if( !fSuccess )
		{
			PrintLog("DuplicateHandle failed \n");
			return FALSE;
		}
		CloseHandle(hChildStdoutRd);
	}

//	pid = _spawnvp(mode, prog_name.c_str(), argv_p);
    
	wxString cmd_line_wx = cmd_line.c_str();
    pid = ::wxExecute( cmd_line_wx, sync_wx);

	if(pid == -1)
	{
		PrintLog(" Error in submitting job %s\n",prog_name.c_str());;
		if(errno == E2BIG)
		{
			PrintLog(" Argument list exceeds 1024 bytes \n");
		}
		else if( errno == EINVAL)
		{
			PrintLog( " Mode argument is invalid \n");
		}
		else if( errno == ENOENT)
		{
			PrintLog( " File or path is not found \n");
		}
		else if( errno == ENOEXEC)
		{
			PrintLog( " Specified file is not executable or has invalid executable-file format \n");
		}
		else if( errno == ENOMEM )
		{
			PrintLog( " Not enough memory is available to execute new process \n");
		}
	}		

	if(get_prog_output)
	{		
		// Duplicate copy of original stdout back into stdout
		
		if (!SetStdHandle(STD_OUTPUT_HANDLE, hSaveStdout)) 
		{
			PrintLog("Re-redirecting Stdout failed\n"); 
			return FALSE;
		}
				
		if (!CloseHandle(hChildStdoutWr)) 
		{
			PrintLog("Closing handle failed \n"); 
			return FALSE;
		}
		
		prog_output.clear();

		std::string empty_str;
				
		if(pid)
		{
			prog_output.push_back(empty_str);
			std::string* pstr = &(prog_output.back());
			int isactive = TRUE;
			DWORD nOutRead;
			
			for (;;) 
			{ 	
				if( !ReadFile( hChildStdoutRdDup, buf, 255, &nOutRead, 
					NULL) || nOutRead == 0) break; 
				for(i=0 ; i < nOutRead; i++)
				{
					if(buf[i] == '\n')
					{
						prog_output.push_back(empty_str);
						pstr = &(prog_output.back());
					}
					else
					{
						(*pstr) += buf[i];
					}
				}
			}
		}
	}
	
	
#else
	PrintLog("Command line = :\n %s \n ",cmd_line.c_str());
	system(cmd_line.c_str());

#endif

	free(argv_p);
	return pid;
}

HaMolSet* HarlemApp::GetMolSetByName(const char* name)
{
	std::string ref_name = name;
	boost::trim(ref_name);
	boost::to_upper(ref_name);

	int nm = molset_vec.size();
	int i;

	for(i = 0; i < nm; i++)
	{
		HaMolSet* pmset = (HaMolSet*) molset_vec[i];
		std::string name2 = pmset->GetName();
		boost::trim(name2);
		boost::to_upper(name2);
		if(name2 == ref_name) return pmset;
	}
	return NULL;
}	

void HarlemApp::AddMolSet(HaMolSet* pmset)
{
	if(pmset != NULL)
	{
		int nm_old = molset_vec.size();
		molset_vec.push_back(pmset);
	}
}

void HarlemApp::DeleteMolSet(HaMolSet* pmset)
{
	if(pmset != NULL)
	{
		VecPtr::iterator mitr;
		mitr = molset_vec.begin();
		for(; mitr != molset_vec.end(); mitr++)
		{
			if( (*mitr) == pmset) 
			{
				molset_vec.erase(mitr);
				return;
			}
		}
	}
}

HaAtom* HarlemApp::GetAtomByRef(const char* at_ref)
{
	  std::string mset_name;
	  int i;
	  int len = strlen(at_ref);
	  HaAtom* aptr = NULL;
	  if(at_ref[0] == '!') 
	  {
	     for(i = 1; i < len; i++)
		 {
		    if( at_ref[i] == '!') break;
		    mset_name += at_ref[i];
		 }
		 HaMolSet* pmset = GetMolSetByName(mset_name.c_str());
		 if( pmset == NULL) return NULL;
		 if( i == len) return NULL;
		 aptr = pmset->GetAtomByRef(&at_ref[i+1]);
	  }
	  return aptr;
}

int HarlemApp::ExecuteScriptFromFile(const char* script_fname)
{
#if !defined(HARLEM_PYTHON_NO)
	FILE* finp;
	finp = fopen(script_fname,"r");
	
	if(finp != NULL )
	{
		char fname_var[256];
		strcpy(fname_var,script_fname);
		PyGILState_STATE gstate;
		gstate = PyGILState_Ensure();
		int ires = PyRun_SimpleFile(finp,fname_var);
		PyGILState_Release(gstate);
		fclose(finp);
		if( ires == 0 )
			return TRUE;
	}
	else
	{
		PrintLog("Can't find script file %s \n", script_fname);
	}
#endif
	
	return TRUE;
}

int HarlemApp::ExecuteScriptInString(const char* script_str)
{
#if !defined(HARLEM_PYTHON_NO)
	char script_str_int[50000];
	strcpy(script_str_int,script_str);
	PyGILState_STATE gstate;
	gstate = PyGILState_Ensure();
	int ires = PyRun_SimpleString(script_str_int);
	PyGILState_Release(gstate);
	if( ires == 0 )
		return True;	
#endif
	return True;
}
