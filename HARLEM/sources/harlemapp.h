/*! \file harlemapp.h

    Application class of HARLEM.

    \author Igor Kurnikov
    \date 1999-2002

*/
#if !defined(HARLEMAPP_H)
#define HARLEMAPP_H

#include "hastring.h"
#include "command.h"

class HaMolView;
class MolSet;
class HaAtom;
class HaMPI;

//!class to describe Remote Computer Account
class ComputerAccount
{
public:
	std::string acc_id;            //!< Account ID to refer to
    std::string login_str;         //!< login string
	std::string interm_acc_ID;     //!< ID of the intermediate account
};

//! Class for the HARLEM application
class HarlemApp 
{
public:
	HarlemApp();
	virtual ~HarlemApp();

	int ProcessOptions();
	int InitFirst();
	int InitLast();
	int InitParallel();   //!< Init Parallel (MPI) Environment
	int InitRemoteComp();
	void InitCommand();   //!< Init RASMOL command Processing 

	int RedirectIOLogWindow(); //!<  Create Log Window (if necessary) and redirect stdout and stderr to it
	int RedirectIOLogFile(const std::string& fname = ""); //!< Redirect stdout and stderr to log file by default: with the name harlem_nproc.log (nproc - MPI rank)

	int CreateCommandWindow(); //!< Create PYTHON command window 
	int ExecuteCommand(); //!< Execute command in CurLine buffer
	int RasMolCmd(const char* cmd); //!< Execute command only in RASMOL command processor
	int ExecRasMolScript(const char* file_name); //!< Execute script consisting of RASMOl-like commands
	int ExecuteScriptFromFile(const char* script_fname);  //!< Execute script from file
	int ExecuteScriptInString(const char* script_str);    //!< Execute script in the string

	static long RunExternalProgram(RunMode rmode, const std::string& prog_name, StrVec& prog_args,
		                                          StrVec& prog_output, int get_prog_output = 0); //!< Run external program 
	static int CheckProcIsActive(long proc_id); //!< Check if the process specified by proc id is running
	static int KillProc(long proc_id);     //!< Terminate the process specified by proc id
	static int SwitchThread();             //!< Switch execution to an another thread 
	static int SleepThread(int ms_delay);  //!< Delay execution of the thread for ms_delay ms 

	int ProcessEvent(int type, int id);    //!< Process Signal (substitution for Event)

	std::vector<ComputerAccount> comp_accounts; //!< accessible computer accounts

	ComputerAccount* GetAccountByID(const char* acc_id);
	int ShowAccountsLoad(); //!< Show load of remote accounts
	int ExecuteRemoteCmd(ComputerAccount* pacc, const char* cmd, StrVec& prog_output, int get_prog_output); //!< Execute Command on a remote computer 

	virtual void Exit() {}

	void LoadInitFile();

	void StartWait(); //!< Display Wait cursor to indicator the program is busy
	void EndWait();   //!< End Wait, restore  normal cursor 

	MolSet* GetMolSetByName(const char* name); //!< Retrieve Molecular Set by name
	void AddMolSet(MolSet* pmset); //!< Add Molecular Set to the list of Molecular Sets
	void DeleteMolSet(MolSet* pmset); //!< Delete Molecular Set from the list of Molecular Sets
	HaAtom* GetAtomByRef(const char* at_ref); //!< Get Atom By Full Reference (including MolSet Name) 

	int gui_mode;        //!< flag for GUI mode of the application to be activated
	int cmd_prompt_mode;  //!< flag for python command mode of the application to be activated
	std::string mpi_py_script; //!<execute python script in parrallel, it's the script job to think about parrallelizm

	int argc_loc;
	char **argv_loc;

    StrIntMap FormatOpt;

	std::string finp_name;
	std::string script_name;
	std::string script_str;

	CmdParser cmd_pr; //!< Rasmol Command Processor
	int only_rasmol_command; //!< Flag to execute only RASMOL type command in CurLine buffer

	int FileFormat;

	std::string harlem_home_dir;   //!< HARLEM home directory
	std::string res_db_dir;        //!< Directory of residue templates
	std::string script_dir;        //!< HARLEM script directory
	std::string basis_dir;         //!< HARLEM quantum chemical gaussian basis set directory
	std::string word_editor;       //!< external word editor name 
	std::string manual_main_page;  //!< URL of the advanced manual page 
	std::string cmd_line_help_main_page;  //!< URL of beginners manual 

	HaMPI* mpi_driver; 
//	void* python_thread;  //!< Thread for Python execution
	
	FILE* file_log; //!< Log file

	std::map<int, int> log_msg_count; // map of counts of printed messages of a given type 
	int max_num_log_msg;              // maximal number of messages of a given type to print 

	VecPtr molset_vec;  //!< Vector of pointers to Molecular Sets in the application
	std::vector<std::shared_ptr<MolSet>> molset_vec_shared; //!< Vector of shared pointers to Molecular Sets in the application
	static HarlemApp* m_HarlemApp;
};

//! HarlemApp starter for python runs
void StartHarlemApp();

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

extern "C" {
#if defined(HARLEMAPP_CPP)
	HarlemApp* pApp;
	HarlemApp* GetHarlemApp();
#else
	extern HarlemApp* pApp;          //!< pointer to the current HARLEM application instance
	extern HarlemApp* GetHarlemApp();   //!< get pointer to the current HARLEM application instance
#endif
}

#endif // !defined(HARLEM_APP_MFC_H)
