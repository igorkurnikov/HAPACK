/*! \file mm_dialogs_wx_1.h

    Dialogs for Molecular Mechanics and Related Modules
 
    \author Igor Kurnikov  
    \date 2010-
*/

#if !defined(MM_DIALOGS_WX_H)
#define MM_DIALOGS_WX_H

class MMBond;
class MMValAngle;
class MMDihedral;
class MolEditor;
class HaMolMechMod;
class MMSysInfo;
class MMInfoDlg;
class RMSDAgent;

class MolMechDlgWX : public wxFrame
{
public:
	MolMechDlgWX(HaMolMechMod* ptr_mm_mod_new, wxWindow *parent);   
    virtual ~MolMechDlgWX();

    static int dlg_open;

	virtual bool TransferDataToWindow();
	virtual bool TransferDataFromWindow();

	void OnInitDialog();

protected:
	HaMolMechMod* ptr_mm_mod;

	MMInfoDlg* p_mm_info_dlg;
	RMSDAgent* p_si_ag;

	MMBond* cur_bond;
	MMValAngle* cur_vang;
	MMDihedral* cur_dih;
	MMDihedral* cur_impr_dih;

   void TransferAtomSuperimposeDataToWindow();   //!< Set Dialog controls from RMSDAgent parameters
   void TransferAtomSuperimposeDataFromWindow(); //!< Set RMSDAgent parameters from Dialog controls
   void TransferExtProgDataToWindow();
   void TransferRunTypeDataToWindow(); //!< Set Dialog elements corresponding to the chosen run_type

public:
   void OnChangePeriodicity(); //!< Update Dialog on creation/deletion of periodical box of the system
protected:
// WDR: handler declarations for MolMechDlgWX 
    void OnChangeSelElem(wxCommandEvent& event);
	void OnSaveExtProgInp(wxCommandEvent& event);
	void OnRunMMCalc(wxCommandEvent& event);
	void OnAmberLoadRestart(wxCommandEvent& event);
	void OnEditchangeRunType(wxCommandEvent& event);
	void OnCalcSinglePtEne(wxCommandEvent& event);
	void OnLoadLogFile(wxCommandEvent& event);
	void ChooseRestrainedAtoms(wxCommandEvent& event);
	void ChooseMovingAtoms(wxCommandEvent& event);
	void OnChooseMDCrdFile(wxCommandEvent& event);
	void OnChooseMDVelFile(wxCommandEvent& event);
	void OnChooseMDEneFile(wxCommandEvent& event);
	void OnChooseConstrTrajFile(wxCommandEvent& event);
	void OnMMStop(wxCommandEvent& event);
	void ChooseRestrRefCrd(wxCommandEvent& event);
	void LoadAtomAtomRestrFile(wxCommandEvent& event);
	void OnPlaybackTrj(wxCommandEvent& event);
	void OnIndexTrj(wxCommandEvent& event);
	void OnSetCurrPt(wxCommandEvent& event);
	void OnChooseMDAnalScript(wxCommandEvent& event);
	void OnEditMDAnalScript(wxCommandEvent& event);
	void OnChkRMSDAnal(wxCommandEvent& event);
	void OnChkMMRunInt(wxCommandEvent& event);
	void OnChoiceRunType(wxCommandEvent& event);
	void ChooseFitAtoms(wxCommandEvent& event);
	void ChooseRMSDAtoms(wxCommandEvent& event);
	void ChooseRefCrdFit(wxCommandEvent& event);
	void ChooseRefCrdRMSD(wxCommandEvent& event);
	void OnEditAmberInp(wxCommandEvent& event);
	void OnEditAmberTop(wxCommandEvent& event);
	void OnEditAmberRun(wxCommandEvent& event);
	void OnEditAmberRst(wxCommandEvent& event);
	void OnAmberSaveRestart(wxCommandEvent& event);
	void OnShowMMInfo(wxCommandEvent& event);
	void OnEditPeriodicBox(wxCommandEvent& event);

	void OnInitMolMech(wxCommandEvent& event);
	void OnConstrainHBond(wxCommandEvent& event);

	void OnSelChangeFFTypeDefault( wxCommandEvent& event );
	void OnUpdateElemList(wxCommandEvent& event);
	void OnSelectFFParamFile(wxCommandEvent& event);
	void OnRadioListChange(wxCommandEvent& event);
	void OnSetNewPar(wxCommandEvent& event);
	void OnDelImprAng(wxCommandEvent& event);
	void OnSendMPIMsg1(wxCommandEvent& event);
	void OnSetDnaCSPars(wxCommandEvent& event);
	void OnSaveResffFromMort(wxCommandEvent& event);
	void OnSaveAtomIndRestrArb(wxCommandEvent& event);
	void OnClose(wxCloseEvent& event);
	void OnChangingPage(wxNotebookEvent& event);
	void OnChangePage(wxNotebookEvent& event);

private:
	DECLARE_EVENT_TABLE()

	wxCheckBox* chk_run_internal;     //!< CheckBox for HaMolMechMod::run_internal_flag 
	wxChoice*  choice_ctrl_ext_prog;  //!< Choice control for MM External Program
	wxChoice*  choice_ctrl_run_type;  //!< Choice control for MM External Program
	wxChoice*  choice_ctrl_per_bcond; //!< Choice control for Periodical Boundary conditions simulations
	wxChoice*  choice_ctrl_electr_method; //!< Choice control for Electrostatic interactions Method 

	wxButton*  btn_save_inp_files;   //!< Button to save input files for external MM program
	wxButton*  btn_mm_run_calc;      //!< Button to run MM calculations

// RMSDAgent setup parameters:

	wxCheckBox* chk_rmsd_anal;        //!< Check button to switch Atom Superimpose analysis on/off
	wxComboBox* combo_ref_crd_fit_type;  //!< Combo box RMSD Reference Coordinates type
	wxComboBox* combo_ref_crd_rmsd_type; //!< Combo box RMSD Reference Coordinates type
	wxComboBox* combo_restr_ref_crd_type; //!< Combo box Restrained Reference Coordinates type 
	wxTextCtrl* txt_fit_atoms;         //!< Text Control with ATOM GROUP ID for atoms to fit 
	wxTextCtrl* txt_rmsd_atoms;        //!< Text Control with ATOM GROUP ID for atoms to compute RMSD
	wxTextCtrl* txt_rmsd_file_name;       //!< Text Control with RMSD data output file name
	wxTextCtrl* txt_ref_crd_fit_fname;   //!< Text Control reference coordinate (to fit atoms) file name
	wxTextCtrl* txt_ref_crd_rmsd_fname;  //!< Text Control reference coordinate (to compute RMSD) file name
	wxCheckBox* chk_rmsd_per_atom;        //!< Check Button to switch on/off RMSD per atom analysis
	wxCheckBox* chk_rmsf_per_atom;        //!< Check Button to switch on/off RMSF per atom analysis
	wxCheckBox* chk_avg_coord;            //!< Check Button to switch on/off calculations of average coordinates
	wxTextCtrl* txt_rmsd_per_atom_file;   //!< Text Control with RMSD per atom output file name
	wxTextCtrl* txt_rmsf_per_atom_file;   //!< Text Control with RMSF per atom output file name
	wxTextCtrl* txt_avg_coord_file;       //!< Text Control with Averaged Coordinates output file name
	wxButton* btn_choose_fit_at;          //!< Button to Choose Atoms to fit      
	wxButton* btn_choose_rmsd_at;         //!< Button to Choose Atoms to compute RMSD
	wxButton* btn_choose_ref_crd_fit;     //!< Button to Choose Reference Coordinates To fit
	wxButton* btn_choose_ref_crd_rmsd;    //!< Button to Choose Reference Coordinates To compute RMSD
	wxButton* btn_choose_rmsd_out_file;   //!< Button to Choose RMSD output file
	wxButton* btn_choose_rmsd_per_atom_file; //!< Button to Choose RMSD per Atom output file
	wxButton* btn_choose_rmsf_per_atom_file; //!< Button to Choose RMSF per Atom output file
	wxButton* btn_choose_avg_coord_file;     //!< Button to Choose Average Coordinates output file

	std::vector<wxWindow*> md_controls;  //!< Controls for Molecular Dynamics simulations
	std::vector<wxWindow*> ene_min_controls; //!< Controls for Energy Minimization
	std::vector<wxWindow*> atom_superimpose_controls; //!< Controls for setup of RMSDAgent
	std::vector<wxWindow*> ext_prog_controls; //!< Controls to use with external MM program

// RMSD dialog:
	std::string fit_atoms_grp_name;
	std::string rmsd_atoms_grp_name;
	
//	int ref_crd_fit_type;
//	int ref_crd_rmsd_type;

	std::string ref_crd_fit_file_name;
	std::string ref_crd_rmsd_file_name;

};

class MMInfoDlg : public wxFrame
{
public:
	MMInfoDlg(HaMolMechMod* ptr_mm_mod_new, wxWindow *parent);   
    virtual ~MMInfoDlg();

	virtual bool TransferDataToWindow();
	virtual bool TransferDataFromWindow();

	void OnInitDialog();

	wxString ene_format; //!< Format string to display energies

	void OnUpdate(wxCommandEvent& event);
	void OnUpdateFormat(wxCommandEvent& event);
	void OnClose(wxCloseEvent& event);
	void OnCloseBtn(wxCommandEvent& event);

protected:
	DECLARE_EVENT_TABLE()

	HaMolMechMod* ptr_mm_mod;
	MMSysInfo*    p_mm_info;

	std::vector<wxWindow*> dbl_text_boxes;
};



#endif