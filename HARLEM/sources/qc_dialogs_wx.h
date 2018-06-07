/*! \file qc_dialogs_wx_1.h

    Dialogs for Quantum Chemical calculations
 
    \author Igor Kurnikov  
    \date 2011-
*/

#if !defined(QC_DIALOGS_WX_H)
#define QC_DIALOGS_WX_H

//----------------------------------------------------------------------------
// EditHuckHamDlg
//----------------------------------------------------------------------------

class EditHuckHamDlg: public wxDialog
{
public:
    // constructors and destructors
    EditHuckHamDlg( ETCouplMod* et_mod, wxWindow *parent, wxWindowID id, const wxString &title,
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxDEFAULT_DIALOG_STYLE );

    virtual bool TransferDataToWindow(); //!< Overwrite fill controls method

    void SetGridData(); //!< Fill Grid with diagonal matrix elements

    wxGrid*  ene_grid; //!< grid with diagonal energies of the hamiltonian
    HaMat_double* pmat; //!< Pointer to the hamiltonian matrix
    HaMat_double* ssmat; //!< Pointer to the overlap matrix
    HaQCMod*      pqc_mod;  //!< Pointer to Quantum Chem Module
    
    // WDR: method declarations for EditHuckHamDlg
    
private:
    // WDR: member variable declarations for EditHuckHamDlg
    
private:
    // WDR: handler declarations for EditHuckHamDlg
    void RefreshView( wxCommandEvent &event );
    void OnGridCellChange( wxGridEvent &event );

private:
    DECLARE_EVENT_TABLE()
};

class LoadQCDatDlgWX : public wxDialog
{
// Construction
public:
	LoadQCDatDlgWX(HaQCMod* new_qcmod, wxWindow* parent);   // standard constructor
        virtual ~LoadQCDatDlgWX();
     
    static int dlg_open;

// Overrides

// WDR: method declarations for LoadQCDatDlgWX

	void OnInitDialog(wxInitDialogEvent& event);
protected:
// WDR: member variable declarations for LoadQCDatDlgWX
	HaQCMod* qcmod;

// WDR: handler declarations for LoadQCDatDlgWX
	void OnLoadData(wxCommandEvent &event);
	void OnClose(wxCloseEvent& event);
private:
	DECLARE_EVENT_TABLE()
};

class WaveFunAnalDlgWX : public wxDialog
{
public:
	WaveFunAnalDlgWX(HaQCMod* new_phost_qcmod, wxWindow* parent);  
    virtual ~WaveFunAnalDlgWX();

    static int dlg_open;

	virtual bool TransferDataToWindow();

	void OnInitDialog(wxInitDialogEvent& event);
	int GetSelMOIdx(); //!< Get Index of selected MO

// WDR: member variable declarations for WaveFunAnalDlgWX 	
	HaQCMod* p_qc_mod;

// WDR: handler declarations for WaveFunAnalDlgWX
	void OnMoveUp(wxCommandEvent &event);
	void OnMoveDown(wxCommandEvent &event);
	void OnClose(wxCloseEvent& event);
	void OnPlotMO(wxCommandEvent &event);
	void OnPrintMOCoef(wxCommandEvent &event);
	void OnMOrefresh(wxCommandEvent &event);
private:
	DECLARE_EVENT_TABLE()
};




#endif