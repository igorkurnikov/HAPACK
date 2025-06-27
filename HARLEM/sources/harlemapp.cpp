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
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>

#include "wx/wxprec.h"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include <wx/process.h>

#include <chrono>
#include <thread>

#if !defined(_MSC_VER)
#include <unistd.h>
#endif

#include "hampi.h"
#include "harlemapp.h"
#include "hamolview.h"
#include "hamolecule.h"
#include "hastring.h"
#include "haio.h"
#include "tokens.h"
#include "abstree.h"
#include "hamolset.h"
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

HarlemApp* HarlemApp::m_HarlemApp = NULL;

void StartHarlemApp()
{
	HarlemApp::m_HarlemApp = new HarlemApp();
	HarlemApp::m_HarlemApp->InitFirst();
	HarlemApp::m_HarlemApp->InitLast();
	pApp = HarlemApp::m_HarlemApp;
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
	only_rasmol_command = TRUE;

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
	file_log      = NULL;

	max_num_log_msg = 50;
}

HarlemApp::~HarlemApp()
{
   int nm = molset_vec.size();
   if( nm > 0 )
   {
	  VecPtr molset_vec_axx = molset_vec;
	  int i;
	  for( i = 0; i < nm; i++)
	  {
		 MolSet* pmset = (MolSet*) molset_vec_axx[i];
		 this->DeleteMolSet(pmset);
	  }
   }
   if( file_log  != NULL)
   {
	   fclose(file_log);
	   file_log = NULL;
   }
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
	std::string doc_dir;

// Set HARLEM HOME directory
	boost::filesystem::path harlem_home_path = ".";

	if( std::getenv("MOLSET_HOME") != NULL ) harlem_home_path = std::getenv("MOLSET_HOME");
	if (std::getenv("HARLEM_HOME") != NULL)  harlem_home_path = std::getenv("HARLEM_HOME");

	harlem_home_dir = harlem_home_path.string() + boost::filesystem::path::preferred_separator;
	res_db_dir = harlem_home_dir + "residues_db" + boost::filesystem::path::preferred_separator;
	doc_dir    = harlem_home_dir + "doc" + boost::filesystem::path::preferred_separator;
	//PrintLog("%s(): harlem_home_dir = %s \n", __func__, harlem_home_dir.c_str());
#if(_MSC_VER)
	word_editor = harlem_home_dir + "scite.exe";
#else
	word_editor = "scite";
#endif
	
	manual_main_page = doc_dir + "advanced_manual_html" + boost::filesystem::path::preferred_separator + "index.html" ;
	cmd_line_help_main_page = doc_dir +  "HARLEM_BeginnerUserManual.htm";
	
	InitParallel();

	ProcessOptions();

	if( mpi_driver->myrank != 0 && mpi_py_script.empty()) return TRUE;
	
	InitCommand();
	InitRemoteComp();

	//RedirectIOToMultipleFilesMPI("test.log");
	//RedirectIOLogWindow();
	//ios::sync_with_stdio();
	
	//cout << "HarlemApp::InitFirst() Check output to cout" << endl;
	//printf("HarlemApp::InitFirst()  Check output to stdout \n");
	//PrintLog(" HarlemApp::InitFirst() Check output to PrintLog \n");
	//cout.flush();

	if( !finp_name.empty() ) // if filename has been specified on the command line: load the file
	{
		MolSet* pmset = new MolSet();
		boost::filesystem::path finp_path = finp_name;

		// PrintLog(" input file name %s \n", finp_path.string().c_str());

		boost::filesystem::path cur_path = boost::filesystem::current_path();
		boost::filesystem::path finp_dir_path = finp_path.parent_path();
		// PrintLog(" HarlemApp::InitFirst() pt 1:  Current Working Directory %s \n",cur_path.string().c_str() );
		if(!finp_dir_path.empty() )
		{
			boost::filesystem::current_path(finp_dir_path);
		}
		cur_path = boost::filesystem::current_path();

		// PrintLog(" HarlemApp::InitFirst() pt 2:  Current Working Directory %s \n", cur_path.string().c_str());
			
		std::string exten = harlem::GetExtFromFileName(finp_name);
		
		boost::to_upper(exten);
		if(exten == "HLM") FileFormat = FormatHarlem;
		else if( exten == "PDB" || exten == "ENT")
			FileFormat = FormatPDB;
		else if( exten == "XYZ" )
			FileFormat = FormatXYZ;
		else if( exten == "HIN" )
			FileFormat = FormatHIN;
		else if (exten == "MOL2")
			FileFormat = FormatMol2;
		else if (exten == "NRG")
			FileFormat = FormatNRG;

		int result = pmset->FetchFile(FileFormat,finp_name.c_str());
		if(!result)
		{
			PrintLog(" HarlemApp::InitFirst() ");
			PrintLog(" Error loading file %s \n",finp_name.c_str());
		}
	}
	if (!mpi_py_script.empty()) ExecuteScriptFromFile(mpi_py_script.c_str());
	
	// PrintLog(" HarlemApp::InitFirst() pt 3:  mpi_driver->myrank = %d \n", mpi_driver->myrank );
	if ( mpi_driver->myrank != 0) return true;

	if (gui_mode)
	{
//		LoadHaPyGUIModules();
	}
	return TRUE;
}

int HarlemApp::InitParallel()
{
	boost::filesystem::path cur_path = boost::filesystem::current_path();
	mpi_driver = new HaMPI();
	if (mpi_driver->nprocs > 1 && mpi_driver->myrank > 0)
	{
		mpi_driver->Listen();
	}
	//Current Directory is changed by MPI_Init() - changing back...
	boost::filesystem::current_path(cur_path);
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
	int only_rasmol_command_old = only_rasmol_command;
	only_rasmol_command = TRUE;
	this->ExecuteCommand();
	only_rasmol_command = only_rasmol_command_old;
	return TRUE;
}

int HarlemApp::ExecRasMolScript(const char* file_name)
{
	char buf[20000];
	std::ifstream sfile(file_name);
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

	if(!only_rasmol_command)
	{
		PyGILState_STATE gstate;
		gstate = PyGILState_Ensure();
		if (PyErr_Occurred()) { // PyErr_Print(); 
		     PyErr_Clear(); }
		int ires = PyRun_SimpleString(cmd_pr.GetCmdLine());
		if (PyErr_Occurred()) { // PyErr_Print(); 
		     PyErr_Clear(); }
		PyGILState_Release(gstate);
		
		if( ires == 0 )
			return TRUE;
		else
		{
//			  PyErr_Print();
			PrintLog("RASMOL type command \n");
		}
	}

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
		
		std::ios::sync_with_stdio();

	    wxLog* p_log = new wxLogStderr();
        wxLog::SetActiveTarget(p_log);

		return TRUE;
	}
	return TRUE;
}


void HarlemApp::InitCommand()
{
	cmd_pr.InitKeywords();	
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
//		PrintLog(" Arg %d = %s\n", i, option.c_str());
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

			if( trunc_opt == "nogui" || trunc_opt == "-nogui")
			{
				gui_mode = FALSE;
				continue;
			}

			if( trunc_opt == "noprompt" || trunc_opt == "-noprompt" )
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

			if( trunc_opt == "script" || trunc_opt == "-script" || trunc_opt == "s" )
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

			if( trunc_opt == "mpipy" || trunc_opt == "-mpipy")
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

			if( trunc_opt == "wd" || trunc_opt == "-wd" )
			{
				i++;
				if( i >= argv_size )
				{
					PrintLog("No Working Directory Name supplied for --wd option\n");
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


static int PrefixString(char* str1, char* str2 )
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

HarlemApp* GetHarlemApp()
{
	if (pApp == NULL)
	{
		StartHarlemApp();
	}
	return pApp;
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
	std::this_thread::sleep_for(std::chrono::milliseconds(100));
    return FALSE;
}

int HarlemApp::SleepThread(int ms_delay)
{
	std::this_thread::sleep_for(std::chrono::milliseconds(ms_delay));
   return FALSE;
}

int HarlemApp::ProcessEvent(int type, int id)
{
	MolSet* pmset = GetCurMolSet();

	// printf(" HarlemApp::ProcessEvent()  myrank = %d \n", this->mpi_driver->myrank);
	// printf("pmset = %s \n", pmset->GetName());

	if (pmset != NULL && pmset->ProcessEvent(type,id)) return true;

	int nm = molset_vec.size();
	int i;
	for (i = 0; i < nm; i++)
	{
		pmset = (MolSet*)molset_vec[i];
		if (pmset == GetCurMolSet()) continue;
		if (pmset->ProcessEvent(type,id)) return true;
	}
	return false;
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

MolSet* HarlemApp::GetMolSetByName(const char* name)
{
	std::string ref_name = name;
	boost::trim(ref_name);
	boost::to_upper(ref_name);

	int nm = molset_vec.size();
	int i;

	for(i = 0; i < nm; i++)
	{
		MolSet* pmset = (MolSet*) molset_vec[i];
		std::string name2 = pmset->GetName();
		boost::trim(name2);
		boost::to_upper(name2);
		if(name2 == ref_name) return pmset;
	}
	return NULL;
}	

void HarlemApp::AddMolSet(MolSet* pmset)
{
	if(pmset != NULL)
	{
		int nm_old = molset_vec.size();
		molset_vec.push_back(pmset);
	}
}

void HarlemApp::DeleteMolSet(MolSet* pmset)
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
		 MolSet* pmset = GetMolSetByName(mset_name.c_str());
		 if( pmset == NULL) return NULL;
		 if( i == len) return NULL;
		 aptr = pmset->GetAtomByRef(&at_ref[i+1]);
	  }
	  return aptr;
}

int HarlemApp::ExecuteScriptFromFile(const char* script_fname)
{
	FILE* finp;
	finp = fopen(script_fname,"r");
	
	if(finp != NULL )
	{
		char fname_var[256];
		strcpy(fname_var,script_fname);
		PyGILState_STATE gstate;
		gstate = PyGILState_Ensure();
		if (PyErr_Occurred()) { // PyErr_Print(); 
		     PyErr_Clear(); }
		
		int ires = PyRun_SimpleString("from molset import *");
		ires = PyRun_SimpleFile(finp,fname_var);

		if (PyErr_Occurred()) { // PyErr_Print(); 
		    PyErr_Clear(); }
		PyGILState_Release(gstate);
		fclose(finp);
		if( ires == 0 )
			return TRUE;
	}
	else
	{
		PrintLog("Can't find script file %s \n", script_fname);
	}
	return TRUE;
}

int HarlemApp::ExecuteScriptInString(const char* script_str)
{
	char script_str_int[50000];
	strcpy(script_str_int,script_str);
	PyGILState_STATE gstate;
	gstate = PyGILState_Ensure();
	if (PyErr_Occurred()) { // PyErr_Print(); 
	      PyErr_Clear(); }
	int ires = PyRun_SimpleString(script_str_int);
	if (PyErr_Occurred()) { // PyErr_Print(); 
	      PyErr_Clear(); }
	PyGILState_Release(gstate);
	if( ires == 0 )
		return True;	
	return True;
}
