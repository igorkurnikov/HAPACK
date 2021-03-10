/*! \file dialogs_wx_1.h

    wxWidgets Dialogs 1 in HARLEM 
 
    \author Igor Kurnikov  
    \date 2003-
*/

#if !defined(DIALOGS_WX_1_H)
#define DIALOGS_WX_1_H

class wxGrid;
class wxGridEvent;
class wxListEvent;
class wxSpinEvent;
class wxSpinButton;

const int IDU_ATOM_PICK = 50001;

// WDR: class declarations

class HaMat_double;
class ETCouplMod;
class HaQCMod;
class HaMolView;
class Object3D;
class MolSet;
class MolEditor;
class HaMolMechMod;
class ElectrostMod;
class HaInterMolMod;
class AtomGroup;
class HaMolecule;
class HaAtom;
class NuclAcidMod;
class HaEmpiricalMod;
class ProtonRedoxMod;
class CrdSnapshot;

//----------------------------------------------------------------------------
// HaFileDlg1
//----------------------------------------------------------------------------

class HaFileDlg1: public wxDialog
{
public:
    // constructors and destructors
    HaFileDlg1( wxWindow *parent, wxWindowID id, const wxString &title,
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxDEFAULT_DIALOG_STYLE );

	//virtual ~HaFileDlg1();

	void OnClose(wxCloseEvent& event);
   
	wxChoice* file_types_ch; //!< choice box for file types choices
    
protected:
    // WDR: member variable declarations for HaFileDlg1
    int nsubdir; //!< number of subdirectories including .. and . listed in the ListCtrl

protected:
    // WDR: handler declarations for HaFileDlg1
    void OnChangeFileType( wxCommandEvent &event );
    void OnSelectFile( wxListEvent &event );
	void OnActivateFile( wxListEvent &event );
    void ChooseDir( wxCommandEvent &event );
    void OnLoadFile( wxCommandEvent &event );

    virtual void FillFileTypes(); //!< Fill File Types filters choice box
    void OnChangeDir();   //!< Fill file list in the current directory

private:
    DECLARE_EVENT_TABLE()
};

//----------------------------------------------------------------------------
// ChooseMolFileDlg
//----------------------------------------------------------------------------

class ChooseMolFileDlg: public HaFileDlg1
{
public:
    // constructors and destructors
    ChooseMolFileDlg( wxWindow *parent, wxWindowID id, const wxString &title,
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxDEFAULT_DIALOG_STYLE );

	//virtual ~ChooseMolFileDlg();

	void OnInitDialog();

	wxBoxSizer *sizer_main_v; //!< Main vertical sizer of the dialog
	
	wxCheckBox* calc_bonds_chk; //!< Calc Bonds Option
	wxTextCtrl* file_name_txt; //!< Text Control for file name
	wxTextCtrl* dir_name_txt;  //!< Text Control for Directory name

	std::string file_name;
	std::string dir_name;

	wxButton* choose_file_btn;

	int file_format; //!< File Format

protected:

    void OnLoadFile( wxCommandEvent &event );
	void OnClose( wxCloseEvent& event );

    virtual void FillFileTypes(); //!< Fill File Types filters choice box
    
private:
    DECLARE_EVENT_TABLE()
};

//----------------------------------------------------------------------------
// SaveMolFileDlg
//----------------------------------------------------------------------------

class SaveMolFileDlg: public HaFileDlg1
{
public:
    // constructors and destructors
    SaveMolFileDlg( MolSet* pmset, wxWindow *parent, wxWindowID id, const wxString &title,
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxDEFAULT_DIALOG_STYLE );

	MolSet* pmset;

    // WDR: method declarations for SaveMolFileDlg
	virtual bool TransferDataToWindow();
	virtual bool TransferDataFromWindow();
    
protected:
    // WDR: member variable declarations for SaveMolFileDlg
    
protected:
    // WDR: handler declarations for SaveMolFileDlg
    void OnSaveFile( wxCommandEvent &event );

    virtual void FillFileTypes(); //!< Fill File Types filters choice box
    
private:
    DECLARE_EVENT_TABLE()
};

//----------------------------------------------------------------------------
// SaveImageFileDlg
//----------------------------------------------------------------------------

class SaveImageFileDlg: public HaFileDlg1
{
public:
    // constructors and destructors
    SaveImageFileDlg( HaMolView* pview, wxWindow *parent, wxWindowID id, const wxString &title,
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxDEFAULT_DIALOG_STYLE );
   

	HaMolView* mol_view;
    // WDR: method declarations for SaveImageFileDlg
    
protected:
    // WDR: member variable declarations for SaveImageFileDlg
    
protected:
    // WDR: handler declarations for SaveImageFileDlg
    void OnSaveFile( wxCommandEvent &event );

    virtual void FillFileTypes(); //!< Fill File Types filters choice box
    
private:
    DECLARE_EVENT_TABLE()
};

//----------------------------------------------------------------------------
// LoadScriptDlgWX
//----------------------------------------------------------------------------

class LoadScriptDlgWX: public HaFileDlg1
{
public:
    // constructors and destructors
    LoadScriptDlgWX( wxWindow *parent, wxWindowID id, const wxString &title,
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxDEFAULT_DIALOG_STYLE );
   

    // WDR: method declarations for LoadScriptDlgWX

public:
    // WDR: member variable declarations for LoadScriptDlgWX
	static int dlg_open;
	wxString script_dir;
    
protected:
    // WDR: handler declarations for LoadScriptDlgWX
	void OnSelectFile( wxListEvent &event );
	void OnExecScriptFile   ( wxCommandEvent &event );
    void OnExecScriptWindow ( wxCommandEvent &event );
    void OnSaveScriptToFile ( wxCommandEvent &event );
	void OnClose( wxCloseEvent& event );

    virtual void FillFileTypes(); //!< Fill File Types filters choice box
    
private:
    DECLARE_EVENT_TABLE()
};

class Object3DDlgWX : public wxDialog
{
// Construction
public:
	Object3DDlgWX(HaMolView* new_pview, wxWindow *parent);   // standard constructor
        virtual ~Object3DDlgWX();

    static int dlg_open;

	virtual bool TransferDataToWindow();

// WDR: method declarations for Object3DDlgWX 

protected:
// WDR: member variable declarations for Object3DDlgWX 
	HaMolView* pview;
	std::vector<Object3D*> obj_vec;

	void DDX_obj_list();

// WDR: handler declarations for Object3DDlgWX 
	void OnSetTransp(wxCommandEvent &event );
	void OnDelete(wxCommandEvent &event );
	void OnDisplay(wxCommandEvent &event );
	void OnUnDisplay(wxCommandEvent &event );
	void OnUpdate(wxCommandEvent &event );
	void OnClose(wxCloseEvent &event);

private:
	DECLARE_EVENT_TABLE()
};

class SolvateDlgWX : public wxDialog
{
// Construction
public:
	SolvateDlgWX(MolSet* new_pmset, wxWindow* parent );   // standard constructor
    virtual ~SolvateDlgWX();

	static int dlg_open;

	virtual bool TransferDataToWindow();
	virtual bool TransferDataFromWindow();

// WDR: member variable declarations for SolvateDlgWX	
	double	m_solv_buf;
	MolSet* pmset;
	MolEditor* p_mol_editor;

	void OnInitDialog(wxInitDialogEvent& event);

// WDR: handler declarations for SolvateDlgWX
	void OnSolvMol(wxCommandEvent &event);
	void OnClose(wxCloseEvent& event);
private:
	DECLARE_EVENT_TABLE()
};

class RedirectIODlgWX : public wxDialog
{
// Construction
public:
	RedirectIODlgWX(wxWindow* parent);   // standard constructor
    virtual ~RedirectIODlgWX();

    static int dlg_open;

	void OnInitDialog(wxInitDialogEvent& event);
protected:

// WDR: handler declarations for RedirectIODlgWX
	void OnStdoutFile(wxCommandEvent &event);
	void OnStdoutConsole(wxCommandEvent &event);
	void OnClose(wxCloseEvent& event);
private:
	DECLARE_EVENT_TABLE()
};


///////////////////////////////////////////////////////////////////////////////////////////
// Dialog to perform ET PATHWAYS calculations

class PathwaysDlgWX : public wxDialog
{
public:
	PathwaysDlgWX(ETCouplMod* new_etmod, wxWindow* parent);   
    virtual ~PathwaysDlgWX();

    static int dlg_open;

	virtual bool TransferDataToWindow();
	virtual bool TransferDataFromWindow();

	void OnInitDialog(wxInitDialogEvent& event);
	
// WDR: member variable declarations for PathwaysDlgWX
	double sel_thresh;
	ETCouplMod* etmod;

// WDR: handler declarations for PathwaysDlgWX
	void OnClose(wxCloseEvent& event);
	void OnCouplMap(wxCommandEvent &event);
	void OnBestPath(wxCommandEvent &event);
	void OnSelCoupled(wxCommandEvent &event);
	void SetModuleDataFromControls (wxCommandEvent &event);
	void OnDutton(wxCommandEvent &event);
	void OnColorMolSurfETCoupl( wxCommandEvent &event);
	void OnSHBond(wxCommandEvent& event);
    void OnHBondTerm(wxCommandEvent& event);

private:
	DECLARE_EVENT_TABLE()
};

class ElectrostDlgWX : public wxDialog
{
public:
	ElectrostDlgWX(ElectrostMod* new_electrost_mod, wxWindow* parent);
	virtual ~ElectrostDlgWX();
	
	static int dlg_open;
	
	virtual bool TransferDataToWindow();
	virtual bool TransferDataFromWindow();
	
	void OnInitDialog(wxInitDialogEvent& event);
	
// WDR: member variable declarations for ElectrostDlgWX
	ElectrostMod* electrost_mod;
	
// WDR: handler declarations for ElectrostDlgWX
	void OnColorSurfPot(wxCommandEvent& event);
	void OnElFieldLoad(wxCommandEvent& event);
	void OnElFieldSave(wxCommandEvent& event);
	void OnAvgPotDon(wxCommandEvent& event);
	void OnClose(wxCloseEvent& event);
	void OnSetParamMod(wxCommandEvent& event);
	void OnEditAtomParams(wxCommandEvent& event);
	void OnSaveInpFiles(wxCommandEvent& event);
	void OnETReorgEne(wxCommandEvent& event);
	void OnCalcRedoxPotShft(wxCommandEvent& event);
	void OnRun(wxCommandEvent& event);
	void OnCalcPotIsoLvl(wxCommandEvent& event);
	void OnCalcIndCharge(wxCommandEvent& event);
	void OnPlotIndCharge(wxCommandEvent& event);
	void OnCalcPotVdwDots(wxCommandEvent& event);
	void OnPotPlaneView(wxCommandEvent& event);
	void OnConcPlaneView(wxCommandEvent& event);
private:
	DECLARE_EVENT_TABLE();
};


///////////////////////////////////////////////////////////////////////////////////////////
// Class to set Parameters of Quantum Chemical Calculations

class QChemParDlgWX : public wxFrame
{
public:
	QChemParDlgWX(HaQCMod* new_phost_qcmod, wxWindow* parent);  
	virtual ~QChemParDlgWX();

    static int dlg_open;          

	void OnInitDialog();

	virtual bool TransferDataToWindow();
	virtual bool TransferDataFromWindow();

	void DDX_bas_name( bool from_window );
    void DDX_active_orb( bool from_window );
	void DDX_wave_fun_type( bool from_window );

// WDR: member variable declarations for QChemParDlgWX
	HaQCMod* p_qc_mod;

// WDR: handler declarations for QChemParDlgWX
	 void OnIpackTest1(wxCommandEvent &event);
	 void OnIpackTest2(wxCommandEvent &event);
	 void OnTestRandom(wxCommandEvent &event);
	 void OnRunCndo(wxCommandEvent &event);
	 void OnRunCndo2(wxCommandEvent &event);
	 void OnConvertChkToFchk(wxCommandEvent& event);
	 void OnStopCalc(wxCommandEvent &event);
	 void OnChangeMethod(wxCommandEvent& event);
	 void OnSetIntEngineIPACK(wxCommandEvent &event);
	 void OnSetIntEngGauss(wxCommandEvent &event);
	 void OnSelChangeBasisSet(wxCommandEvent& event);
	 void OnSelChangeActiveOrb(wxCommandEvent& event);
	 void OnClose(wxCloseEvent &event);
	 void OnSetPar(wxCommandEvent &event);
	 void OnSaveInpFile(wxCommandEvent &event);
	 void OnRunCalc(wxCommandEvent &event);
private:
	DECLARE_EVENT_TABLE()
};

/////////////////////////////////////////////////////////////////////////////
// ETEffHamDlgWX  dialog

class ETEffHamDlgWX : public wxFrame
{
public:
	ETEffHamDlgWX(ETCouplMod* new_ptr_et_mod, wxWindow* parent);
        virtual ~ETEffHamDlgWX();

    static int dlg_open;

	void OnInitDialog();
	
	virtual bool TransferDataToWindow();
	virtual bool TransferDataFromWindow();

// WDR: member variable declarations for ETEffHamDlgWX
	ETCouplMod* ptr_et_mod;

	void DDX_eig_val(int set_var);
	void DDX_tun_ener(int set_var);
	void DDX_da_field(int set_var);
	void DDX_list_loc_orb_gf(int set_var);

// WDR: handler declarations for ETEffHamDlgWX

	void OnSetParam(wxCommandEvent& event);
	void OnClearHeff(wxCommandEvent& event);
	void OnCalcHeff(wxCommandEvent& event);
	void OnCalcScanEne(wxCommandEvent& event);
	void OnDiagHeff(wxCommandEvent& event);
	void OnMinSplit(wxCommandEvent& event);
	void OnSaveDB(wxCommandEvent& event);
	void OnBuildDB(wxCommandEvent& event);
	void OnSaveHeffXml(wxCommandEvent& event);
	void OnLoadHeffXml(wxCommandEvent& event);
	void OnSaveRedoxOrbsXml(wxCommandEvent& event);
	void OnLoadDonorOrbsXml(wxCommandEvent& event);
	void OnLoadAccOrbsXml(wxCommandEvent& event);
	void OnChooseRedoxOrbs(wxCommandEvent& event);
	void OnPickDonorOrbBasis(wxCommandEvent& event);
	void OnPickAccOrbBasis(wxCommandEvent& event);
	void OnResetSrcLorb(wxCommandEvent& event);
	void OnResetTgtLorb(wxCommandEvent& event);
	void OnCalcGFfromMO(wxCommandEvent& event);
	void OnCalcGFfromHeff(wxCommandEvent& event);
	void OnPrintProtectMat(wxCommandEvent& event);
	void OnPlotEigenVec(wxCommandEvent& event);
	void OnCalcHDAfromGF(wxCommandEvent& event);
	void OnCalcHDAPert(wxCommandEvent& event);
	void OnPrintEigVec(wxCommandEvent& event);
	void OnZeroLong(wxCommandEvent& event);
	void OnPrintOvlpElem(wxCommandEvent& event);
	void OnPrintHeffElem(wxCommandEvent& event);
	void OnHeffEditDiag(wxCommandEvent& event);
	void OnSetDonorTruncMo(wxCommandEvent& event);
	void OnCalcHeffHuck(wxCommandEvent& event);
	void OnClose(wxCloseEvent& event);
//	void OnRightDown(wxMouseEvent& event);
	
	void FillListBoxHeffEigene();
private:
	DECLARE_EVENT_TABLE()
};

class SelectLocOrbDlgWX : public wxDialog
{ 
public:
	SelectLocOrbDlgWX(ETCouplMod* new_ptr_et_mod, int new_orb_type, wxWindow* parent); 
    virtual ~SelectLocOrbDlgWX();

    static int dlg_open;

	virtual bool TransferDataToWindow();

	ETCouplMod* ptr_et_mod;
	HaQCMod* ptr_qc_mod;
	wxWindow* phost_dialog;
	int orb_type;

	void DDX_lorb_list(int id);

	void OnClose(wxCloseEvent& event);
	void OnExtractSel(wxCommandEvent& event);

private:
	DECLARE_EVENT_TABLE()
};


class EditFragmDlgWX : public wxFrame
{
public:
	EditFragmDlgWX(MolSet* new_pmset, wxWindow* parent); 
        virtual ~EditFragmDlgWX();

	static int dlg_open;
	
	void OnInitDialog();
	
	virtual bool TransferDataToWindow();
	virtual bool TransferDataFromWindow();

// WDR: member variable declarations for EditFragmDlgWX    
	MolSet* pmset;

// WDR: handler declarations for EditFragmDlgWX 	
	void OnSaveFragm(wxCommandEvent& event);
	void OnCreateFragm(wxCommandEvent& event);
	void OnAddFragment(wxCommandEvent& event);
	void OnSelectAtomsMatchFrag(wxCommandEvent& event);
	void OnClose(wxCloseEvent& event);

private:
	DECLARE_EVENT_TABLE()
};

class BuildFilmDlgWX : public wxFrame
{
// Dialog Class to build a 2-D layer of absorbed molecules
// moldecules 
public:
	BuildFilmDlgWX(MolSet* new_pmset, wxWindow* parent); 
    virtual ~BuildFilmDlgWX();

    static int dlg_open;

// WDR: member variable declarations for BuildFilmDlgWX
	
	double dx; //  increment between molecules along x direction
	double dy; //  increment metween molecules along y direction
	int nx;    //  the number of molecules along X axis
	int ny;    //  the number of molecules along Y axis
	
	double alpha; // angle between x and y directions
	double tilt;  // inclination angle of the molecules to the surface
	
	int nunit; // the number of monomers in the adsorbed molecule

	int num_surf_layers; // the Number of Surface layers

	bool add_surf_below_flag;
	bool add_surf_top_flag;
	bool add_atom_top_flag;
	bool add_atom_below_flag;

	MolSet* pmset;

	void OnInitDialog();

// WDR: handler declarations for BuildFilmDlgWX 
	void OnCreateFilm(wxCommandEvent& event);
	void OnCreateSurf(wxCommandEvent& event);
	void OnCloseBtn(wxCommandEvent& event);
	void OnClose(wxCloseEvent& event);

private:
	DECLARE_EVENT_TABLE()
};

class CrdSnapshotDlg : public wxFrame
{
//! Dialog Class to manipulate Coordinate Snaphots 
public:
	CrdSnapshotDlg(MolSet* new_pmset, wxWindow* parent); 
    virtual ~CrdSnapshotDlg();

	virtual bool TransferDataFromWindow();
	virtual bool TransferDataToWindow();

    static int dlg_open;

	std::map<std::string, bool> prop_show_flags; //!< Show flags of Snapshot  Properties
	std::map<std::string, bool> prop_col_num;    //!< Column numbers of Snapshot  Properties

	MolSet* pmset;
	
	std::map <CrdSnapshot*, std::map<std::string,double> > snap_data ;
	std::map<std::string, CrdSnapshot*> snap_id_ptr_map;

	wxGrid*  snap_list;
	wxTextCtrl* sel_snap_id;
	wxTextCtrl* sel_snap_desc;
	wxSpinButton* move_snap_btn;

	void OnInitDialog();

	void SetColumns();
	void OnChangeProp(wxCommandEvent& event);

	void OnAddSnapshot(wxCommandEvent& event);
	void OnDelSnapshot(wxCommandEvent& event);
	void OnDelAllSnapshots(wxCommandEvent& event);
	void OnSetCrdFromSnapshot(wxCommandEvent& event);
	void OnSaveCrdToSnapshot(wxCommandEvent& event);
	void OnLoadSnapshotFromMolFile(wxCommandEvent& event);
	void OnLoadSnapshotsFromXMLFile(wxCommandEvent& event);
	void OnSaveSnapshotsToXMLFile(wxCommandEvent& event);
	void OnSnapMoveBtnUp( wxSpinEvent& event );
	void OnSnapMoveBtnDown( wxSpinEvent& event );
	void OnChangeSelSnapshot(wxCommandEvent& event);
	void OnEditSnapID( wxCommandEvent& event );
	void OnClose(wxCloseEvent& event);

	void OnChangeNumSnapshots(); //!< Function to call whem number of snapshots changed
	CrdSnapshot* GetSelSnapshot(); //!< Get Selected Snapshot

private:
	DECLARE_EVENT_TABLE()
};



class AtomParamsDlgWX : public wxFrame
{
public:
	AtomParamsDlgWX(MolSet* new_pmset, wxWindow* parent); 
    virtual ~AtomParamsDlgWX();

	virtual bool TransferDataFromWindow();
	virtual bool TransferDataToWindow();

    static int dlg_open;

	static bool coord_edit_flag;
	static bool element_edit_flag;
	static bool charge_edit_flag;
	static bool ff_symbol_edit_flag;
	static bool mass_edit_flag;
	static bool atname_edit_flag;
	static bool radius_edit_flag;
	static bool vdwrad_edit_flag;
	static bool vdwene_edit_flag;
	static bool temperature_edit_flag;
	static bool hybridization_edit_flag;
	static bool hb_status_edit_flag;
	static bool solv_access_area_flag;
       
	static int ResetEditFlags();

// Column numbers of Atom Properties:
	int n_x_coord;
	int n_y_coord;
	int n_z_coord;

	int n_element;
	int n_charge;
	int n_ff_symbol;
	int n_mass;
	int n_atname;
	int n_radius;
	int n_vdwrad;
	int n_vdwene;
	int n_temp;
	int n_hybrid;
	int n_hb_status;
	int n_solv_access;

	void OnInitDialog();

	MolSet* pmset;
	MolEditor* p_mol_editor;
	ProtonRedoxMod* p_prot_rdx_mod;
	wxGrid* m_atom_lctrl;
	VecPtr at_ptrs;

	void SetColumns();
	void FillAtomGroup();

// WDR: handler declarations for AtomParamsDlgWX
	void OnChangeProp(wxCommandEvent& event);
	void OnCheckPerBox(wxCommandEvent& event);
	void OnEdtAtAmberCh(wxCommandEvent& event);
	void OnEdtAtHBondStatus(wxCommandEvent& event);
	void OnCalcSolvAccessArea(wxCommandEvent& event);
	void OnCalcSolvAccessAreaGeoball(wxCommandEvent& event);
	void OnSaveCrdExclVolArb(wxCommandEvent& event);
	void OnCreateExclVolMol(wxCommandEvent& event);
	void OnEdtAtAmberFFSymb(wxCommandEvent& event);
	void OnClearAtomFFParams(wxCommandEvent& event);
	void OnEdtAtElemMass(wxCommandEvent& event);
	void OnSetAtElemFromTempl(wxCommandEvent& event);
	void OnCopyToChmap(wxCommandEvent& event);
	void OnSetFromChmap(wxCommandEvent& event);
	void OnEditAtPhCh(wxCommandEvent& event);
	void OnEditAtFormChTempl(wxCommandEvent& event);
	void OnMinDist(wxCommandEvent& event);
	void OnSetSpecRad(wxCommandEvent& event);
	void OnLoadCoords(wxCommandEvent& event);
	void OnSetStdZMat(wxCommandEvent& event);
	void OnPrintZMat(wxCommandEvent& event);
	void OnSetZMatIntCrdFromAtCrd(wxCommandEvent& event);
	void OnSetAtCrdFromZMat(wxCommandEvent& event);
	void OnCalcDipole(wxCommandEvent& event);
    void OnSetVdwRad(wxCommandEvent& event);
	void OnSetParseRad(wxCommandEvent& event);
	void OnSetHPPRad(wxCommandEvent& event);
	void OnClose(wxCloseEvent& event);
	void OnCloseBtn(wxCommandEvent& event);
	void OnSetStdCharges(wxCommandEvent& event);
	void OnSetZeroCharges(wxCommandEvent& event);
	void OnUpdateAtomGroup(wxCommandEvent& event);
	void OnTransferToWin(wxCommandEvent& event);
	void OnTransferFromWin(wxCommandEvent& event);
	void OnEndLabelEdit(wxGridEvent& event);

private:
	DECLARE_EVENT_TABLE()
};

class ResidueParamsDlgWX : public wxFrame
{
public:
	ResidueParamsDlgWX(MolSet* new_pmset, wxWindow* parent); 
    virtual ~ResidueParamsDlgWX();

	static int dlg_open;

	static bool	res_name_flag;
	static bool	res_num_flag;
	static bool	res_name_modifier_flag;
	static bool	chain_name_flag;

	static int ResetEditFlags();

// Column numbers of Residue Properties:

	int n_res_name;
	int n_res_name_modifier;
	int n_res_num;
	int n_chain_name;
	
	virtual bool TransferDataToWindow();

	void OnInitDialog();

protected:
	MolSet* pmset;
	wxGrid* m_residue_lctrl;
	VecPtr  res_ptrs;
	MolEditor* p_mol_editor;
	
	void SetColumns();
	void FillResidueList();

// WDR: handler declarations for ResidueParamsDlgWX 
	void OnChangeProp(wxCommandEvent& event);
	void OnClose(wxCloseEvent& event);
	void OnCloseBtn(wxCommandEvent& event);
	void OnCheckStruct(wxCommandEvent& event);
	void OnDelExtraAt(wxCommandEvent& event);
	void OnAddMissAt(wxCommandEvent& event);
	void OnFixBondsUsingTempl(wxCommandEvent& event);
	void OnOrderAtomsInRes(wxCommandEvent& event);
	void OnRenameAtomsToAmber(wxCommandEvent& event);
	void OnRenameAtomsToGromacs(wxCommandEvent& event);
	void OnConvertWaterArrowVB(wxCommandEvent& event);
	void OnUpdateResidueList(wxCommandEvent& event);
	void OnResidueRenumber(wxCommandEvent& event);
	void OnEndLabelEdit(wxGridEvent& event);

private:
	DECLARE_EVENT_TABLE()
};

class EditGeomDlgWX : public wxFrame
{
public:
	EditGeomDlgWX(MolSet* new_pmset, wxWindow* parent);
    virtual ~EditGeomDlgWX();

    static int dlg_open;
    static EditGeomDlgWX* active_dlg_ptr;

	void OnInitDialog();
protected:
// WDR: member variable declarations for EditGeomDlgWX

	MolSet* pmset;
	MolEditor* p_mol_editor;

	HaAtom* aptr1;
	HaAtom* aptr2;
	HaAtom* aptr3;
	HaAtom* aptr4;

	void GetAtoms();

	wxTextCtrl* edit_at1;
	wxTextCtrl* edit_at2;
	wxTextCtrl* edit_at3;
	wxTextCtrl* edit_at4;

	wxTextCtrl* active_edit_at;

// WDR: handler declarations for EditGeomDlgWX
	void OnClose(wxCloseEvent& event);
	void OnCalcGeom(wxCommandEvent& event);
	void OnSetGeom(wxCommandEvent& event);
	void OnIncremTors(wxCommandEvent& event);
	void OnDelBond(wxCommandEvent& event);
	void OnCreateBond(wxCommandEvent& event);
	void OnDelAtom(wxCommandEvent& event);
	void OnMergeMol(wxCommandEvent& event);
	void OnAttachMol(wxCommandEvent& event);
	void OnAtomSelect(wxCommandEvent& event);
	
	void SetIntFocusAt(wxFocusEvent& event);

private:
	DECLARE_EVENT_TABLE()
};

class EditGroupsDlg : public wxFrame
{
public:
	EditGroupsDlg(MolSet* new_pmset, int itype, wxWindow* parent);  
    virtual ~EditGroupsDlg();

    static int dlg_open;
    static EditGroupsDlg* active_dlg_ptr;
 
	virtual bool TransferDataToWindow();

	void OnInitDialog();
	
	int ShowModal(); //!< Run in modal mode to select atom group
protected:
	MolSet* pmset;

//	AtomGroup* sel_at_list; //!< Selected Atom List

	int modal_run_flag; //!< Flag to indicate modal run
	wxEventLoop* p_loc_event_loop; //!< Local Event loop while in modal mode

	int igrp_type; //!< type of the atom groups edited
	enum {CHEM_GRP_TYPE=0, REDOX_CNT_TYPE=1, NAMED_ATOM_GROUP_TYPE = 2 };
	
	void OnTestRadio(wxCommandEvent& event);
	void OnGroupType(wxCommandEvent& event);
    void OnPickLevel(wxCommandEvent& event);
	void OnSetSelFromGrp(wxCommandEvent& event);
	void OnAddSelToGrp(wxCommandEvent& event);
	void OnDelSelFromGrp(wxCommandEvent& event);
	void OnCalcDonAccDist(wxCommandEvent& event);
	void OnRenameGrp(wxCommandEvent& event);
	void OnCopyGrp(wxCommandEvent& event);
	void OnSetFromExpr(wxCommandEvent& event);
	void OnAddFromExpr(wxCommandEvent& event);
	void OnDelFromExpr(wxCommandEvent& event);
	void OnKeepOnlyExpr(wxCommandEvent& event);
	void OnCreateExtChMol(wxCommandEvent& event);
	void OnCloseBtn(wxCommandEvent& event);
	void OnClose(wxCloseEvent& event);
	void OnNewGroup(wxCommandEvent& event);
	void OnDelGroup(wxCommandEvent& event);
	void OnAtomSelect(wxCommandEvent& event);
	void OnResetMemb(wxCommandEvent& event);
	void OnDelMemb(wxCommandEvent& event);
	void DisplaySelectedGroup();
	void OnSetProt(wxCommandEvent& event);
	void OnSaveXYZFile(wxCommandEvent& event);
	void OnSaveNDXFile(wxCommandEvent& event);
	void OnSortGrpIdx(wxCommandEvent& event);
	void OnStdGroups(wxCommandEvent& event);
	void OnStdProteinGroups(wxCommandEvent& event);
	void OnRenumberGrp(wxCommandEvent& event);
	void OnChangeSelGroup2(wxCommandEvent& event);
	void OnChangeSelGroup();
	void OnColorRigidClusters(wxCommandEvent& event);
	
	void FillGroupList();
	void UpdateNumAtomCounter();
public:
	AtomGroup* GetSelGroup(); 
private:
	DECLARE_EVENT_TABLE()
};

class CompAccountsDlg : public wxFrame
{
public:
	CompAccountsDlg(wxWindow* parent); 
    virtual ~CompAccountsDlg();

    static int dlg_open;

	virtual bool TransferDataToWindow();
protected:

	void OnAccShowLoad(wxCommandEvent& event);
	void OnClose(wxCloseEvent& event);
	void OnCloseBtn(wxCommandEvent& event);
private:
	DECLARE_EVENT_TABLE()
};

class ResDBDlg : public wxFrame
{
// Construction
public:
	ResDBDlg(wxWindow* parent);   
    virtual ~ResDBDlg();

    static int dlg_open;
	

	void OnInitDialog();
protected:

	HaMolecule* sel_templ;

	void OnClose(wxCloseEvent& event);
    void OnCloseBtn(wxCommandEvent& event);
	void OnUpdateTemplList(wxCommandEvent& event);
	void UpdateTemplList();
	void OnResDBLoad(wxCommandEvent& event);
	void OnChangeSelTempl(wxCommandEvent& event);
private:
	DECLARE_EVENT_TABLE()
};




class MolViewParDlg : public wxFrame
{
public:
	MolViewParDlg(HaMolView* new_pview, wxWindow* parent);
    virtual ~MolViewParDlg();

    static int dlg_open;

	virtual bool TransferDataToWindow();
	virtual bool TransferDataFromWindow();

	void OnInitDialog();
protected:
// WDR: member variable declarations for MolViewParDlg

	HaMolView* pview;

// WDR: handler declarations for MolViewParDlg

	void OnClose(wxCloseEvent& event);
    void OnCloseBtn(wxCommandEvent& event);
	void OnResetViewRot(wxCommandEvent& event);
	void OnSetControlPar(wxCommandEvent& event);
	void OnSetViewPar(wxCommandEvent& event);
private:
	DECLARE_EVENT_TABLE()
};

class MolSetParDlg : public wxFrame
{
public:
	MolSetParDlg(MolSet* new_pmset, wxWindow* parent);   // standard constructor
    virtual ~MolSetParDlg();

    static int dlg_open;

	virtual bool TransferDataToWindow();
	virtual bool TransferDataFromWindow();

	void OnInitDialog();
protected:
// WDR: member variable declarations for MolSetParDlg
	MolSet* pmset;

// WDR: handler declarations for MolSetParDlg
	void OnClose(wxCloseEvent& event);
    void OnCloseBtn(wxCommandEvent& event);
	void OnMultMol(wxCommandEvent& event);
	void OnSetControlPar(wxCommandEvent& event);
	void OnSetMolPar(wxCommandEvent& event);
private:
	DECLARE_EVENT_TABLE()
};

class DValColorMap;

class AtomPropColorDlg : public wxFrame
{
public:
	AtomPropColorDlg(HaMolView* new_pview, wxWindow* parent);   // standard constructor
    virtual ~AtomPropColorDlg();

    static int dlg_open;

	virtual bool TransferDataToWindow();
	virtual bool TransferDataFromWindow();

	void OnInitDialog();
protected:

	HaMolView* pview;
	MolSet* pmset;

	DValColorMap color_map;
	void OnClose(wxCloseEvent& event);
	void OnColorAtomsProp(wxCommandEvent& event);
	void OnClearColorMap(wxCommandEvent& event);
	void OnAddColorByName(wxCommandEvent& event);
	void OnAddColorByRGB(wxCommandEvent& event);
	void OnSaveColorMapFile(wxCommandEvent& event);
	void OnLoadColorMapFile(wxCommandEvent& event);
	void OnGridCellChange(wxGridEvent& event);

	wxTextCtrl* p_txt_minval;
	wxTextCtrl* p_txt_num_interpol_colors;
	wxChoice* p_atom_prop;
	wxChoice* p_color_names;
	wxGrid* p_color_map_grid;

	wxTextCtrl* p_txt_r;
	wxTextCtrl* p_txt_g;
	wxTextCtrl* p_txt_b;

	long r_val;
	long g_val;
	long b_val;

	long num_interpol_colors;

	double min_prop_val;

	int n_col_mval;
	int n_col_cname;
	int n_col_r;
	int n_col_g;
	int n_col_b;
	
private:
	DECLARE_EVENT_TABLE()
};


class NuclAcidDlgWX : public wxFrame
{
// Construction
public:
	NuclAcidDlgWX(NuclAcidMod* new_jumna_mod, wxWindow* parent);   
    virtual ~NuclAcidDlgWX();

    static int dlg_open;
	
	virtual bool TransferDataToWindow();
	virtual bool TransferDataFromWindow();

	void OnInitDialog();
protected:
// WDR: member variable declarations for NuclAcidDlgWX
	enum DISPLAY_MODE { JUMNA_HLX_CRD = 1, BASE_BASE_CRD, GLOB_INTER_BASE_CRD, GLOB_INTER_BP_CRD,
		LOC_INTER_BASE_CRD, LOC_INTER_BP_CRD} display_mode;

	MolSet* pmset;
	NuclAcidMod* nucl_acid_mod;
	wxGrid* m_Grid;

// WDR: handler declarations for NuclAcidDlgWX
	void OnPageChange(wxNotebookEvent& event);
	void OnClose(wxCloseEvent& event);
    void OnCloseBtn(wxCommandEvent& event);
	void OnUpdateControls(wxCommandEvent& event);
	void OnSaveChanges(wxCommandEvent& event);
	void OnBuildNuclAcid(wxCommandEvent& event);
	void OnPdistToCnt(wxCommandEvent& event);
	void OnMinEne(wxCommandEvent& event);
	void OnCalcEne(wxCommandEvent& event);
	void OnGenComplStr(wxCommandEvent& event);
	void OnLocBaseStepPars(wxCommandEvent& event);
	void OnLocHlxParsBp(wxCommandEvent& event);
	void OnCalcAxis(wxCommandEvent& event);
	void OnGlobHlxPar(wxCommandEvent& event);
	void OnCalcBbCrd(wxCommandEvent& event);
	void OnUpdateContent();
	void OnBaseBasePars(wxCommandEvent& event);
	void OnSetFromJumna(wxCommandEvent& event);
    void OnSetToJumna(wxCommandEvent& event);
    void OnUpdateXYZ(wxCommandEvent& event);
	void OnGridEndEdit(wxGridEvent& event);
	void OnLockSelCrd(wxCommandEvent& event);
private:
	DECLARE_EVENT_TABLE()
};

class HaScatterMod;

class ScatterDlg  : public wxFrame
{
public:
	ScatterDlg (HaScatterMod* new_ptr_sc_mod, wxWindow* parent);
    virtual ~ScatterDlg ();

    static int dlg_open;

	virtual bool TransferDataToWindow();
	virtual bool TransferDataFromWindow();

	void OnInitDialog();
protected:
// WDR: member variable declarations for ScatterDlg 

	HaScatterMod* ptr_sc_mod;
	double x1,y1,z1;
	double x2,y2,z2;

// WDR: handler declarations for ScatterDlg 

	void OnClose(wxCloseEvent& event);
	void OnCalcPspXYZ(wxCommandEvent& event);
	void OnCalcPsPAO(wxCommandEvent& event);
	void OnCalcPsPGrid(wxCommandEvent& event);
	void OnSavePSPFile(wxCommandEvent& event);
	void OnGridEigVal(wxCommandEvent& event);
	void OnSetCoreHam(wxCommandEvent& event);
	void OnPSPGaussXYZ(wxCommandEvent& event);
	void OnKinEneGaussXYZ(wxCommandEvent& event);

private:
	DECLARE_EVENT_TABLE()
};

#endif // !defined(DIALOGS_WX_H)
