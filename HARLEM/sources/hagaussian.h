/*! \file hagaussian.h

    Classes to provide interface to GAUSSIAN from HARLEM.

    \author Igor Kurnikov  
    \date 1998-2002

*/
#ifndef HAGAUSSIAN_H
#define HAGAUSSIAN_H

#include "hacompmod.h" 

namespace harlem
{
	class RunOptions;
};

class HaGaussMod: public HaCompMod
//!  Class to control quantum chemical calculations with the GAUSSIAN program   
{
public:

	HaGaussMod(MolSet* new_phost_mset);
	virtual ~HaGaussMod();

	friend class QChemParDlgWX;

	void SetStdFileNames();  //!< Set default names of Gaussian files
	void SetStdJobFlags(); //!< Set default job parameters
	bool SaveInpFile();  //!< Save Gaussian input files
	int LoadOutput();    //!< Load Information from Output Files
	int LoadOutFile(const std::string& out_fname); //!< Load calculation results from GAUSSIAN out file
	int LoadOutFromStream( std::istream& is ); //!< Load calculation results from GAUSSIAN out file contained in a stream 
	int LoadOutSummary( std::string summary_str ); //!< Load info from summary string of Gaussian output 

	int Run(const harlem::RunOptions* popt = NULL); //!< Run GAUSSIAN program ( if sync == TRUE - no return before GAUSSIAN process exited )
	int RunFormChk(const char* fname_chk, const char* fname_fchk); //!< Convert checkpoint(or rwf) file to Formatted form using formchk program

	static void PrintCurBcommon(); //!< Print Gaussian common/B/ 
    
	std::string gaussian_version;  //!< Version of Gaussian
	std::string gaussian_exe;      //!< Name of Gaussian Executable

	int SetFilePrefix(const char* prefix); //!< Set prefix for Gaussian input/output files
	std::string GetFilePrefix() const; //!< Get prefix for Gaussian input/output files

	std::string GetInpFileName() const;  //!< Get Input File name
	std::string GetCHKFileName() const;  //!< Get Checkpoint (*.chk)  File name
	std::string GetFCHKFileName() const; //!< Get Formatted Checkpoint (*.fchk)  File name
	std::string GetRWFFileName() const;  //!< Get Read-Write (*.rwf) file name
	std::string GetOutFileName() const;  //!< Get Output (*.out) file name

	void SetAddKWStr( const string& add_kw_str ); //!< Set string of additional GAUSSIAN keywords

protected:

	std::string inp_file_prefix; //!< Prefix for Gaussian input files
//	std::string inp_file;        //!< Gaussian input file name
//	std::string chk_file;        //!< Checkpoint file name
//	std::string fchk_file;       //!< Formatted checkpoint file name
//	std::string rwf_file;        //!< RWF file name
//	std::string out_file;        //!< out file name
	std::string add_kw_str;      //!< String with additional GAUSSIAN KEYWORDS to add 

public:

	void SetNumSharedMemCores(int n_sh_cores_new); //!< Set the number of Shared Memory Cores per MPI/LINDA process 
	void SetNumProc(int n_proc_new);               //!< Set the number of MPI/LINDA processors (each having several shared memory cores)
	void SetMaxMem(int max_mem );                  //!< Set maximal memory for GAUSSIAN calculations
	int  GetNumSharedMemCores() const;             //!< Get the number of Shared Memory Cores per MPI/LINDA process 
	int  GetNumProc() const;                       //!< Get the number of MPI/LINDA processors (each having several shared memory cores)
	int  GetMaxMem()  const;                       //!< Get Maximal memory for GAUSSIAN calculations (in Mb)

	void SetLoadNonOptGeom( bool set_par = true ); //!< Set Loading last non-optimal geometry for not converged minimization runs  
	void SetLoadGeomZMatOrient( bool set_par = true ); 
	void SetLoadGeomStdOrient ( bool set_par = true );
	void SetReadInitGeomChkFile( bool set_par = true );
	void SetReadHFGuessChkFile( bool set_par = true );
	void SetNoStdOrient( bool set_par = true );
	void SetSaveBasisSetGen( bool set_par = true ); //!< Set Saving Basis Set as generic 

	

protected:
 
	HaQCMod* p_qc_mod;     //!< The pointer to Quantum Chemical module associated with the Gaussian module 
	
	bool pseudo_pot_flag;     //!< flag to use pseudopotential in Gaussian calculations

	void FillSectionProcCommands(std::ostream& os) const; //!< Fill proc commands section of specification section of Gaussian input file
	void FillSectionJob(std::ostream& os) const;    //!< Fill job specification section of Gaussian input file
	void FillSectionCoord(std::ostream& os) const;  //!< Fill atom coordinate section of Gaussian input file
	void FillSectionBasis(std::ostream& os) const;  //!< Fill basis section of Gaussian input file
	void FillSectionExtCharges(std::ostream& os) const; //!< Fill external charges section of Gaussian input file
	void FillSectionInitMO(std::ostream& os) const; //!< Fill initial MOs section of Gaussian input file

	int n_sh_cores;       //!< Number of shared Memory Cores per MPI/LINDA process
	int n_proc;           //!< Number of MPI/LINDA processors (each having several shared memory cores)
	int max_mem;          //!< Maximal memory for GAUSSIAN calculations in MB   
	
	static harlem::RunOptions run_opt_default;
	
	bool no_std_orient;     //!< flag to prevent transformation to standard orientation and use of symmetry
	bool load_non_opt_geom; //!< flag to last non-optimal geometry for not converged minimization runs
	bool load_geom_zmat_orient; 
	bool load_geom_std_orient;
	bool read_init_geom_chk_file; //!< read initial geometry from checkpoint file
	bool read_hf_guess_chk_file;  //!< read inital guess of Hatree Fock functions from Checkpoint File 
	bool save_basis_set_gen; //!< save basis set to gaussian input file as generic 

};


#endif /* !HAGAUSSIAN_H */
