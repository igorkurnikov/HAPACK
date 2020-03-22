/*! \file wx_prot_redox_dlg.h

    Dialogs for Protonation aand Redox Equilibrium Calculations
 
    \author Igor Kurnikov  
    \date 2010-
*/

#if !defined(WX_PROT_REDOX_DLG_H)
#define WX_PROT_REDOX_DLG_H

class ProtonRedoxDlg : public wxFrame
{
public:
	ProtonRedoxDlg( ProtonRedoxMod* p_prot_rdx_mod_new, wxWindow* parent); 
    virtual ~ProtonRedoxDlg();

	static int dlg_open;
	virtual bool TransferDataToWindow();
	void OnInitDialog();

protected:
	MolSet* pmset;
	wxGrid* p_alt_st_grid;
	MolEditor* p_mol_editor;
	ProtonRedoxMod* p_prot_rdx_mod;

	int nc_alt_st_idx;   //!< Column number of Alt State Index
	int nc_alt_st_type;  //!< Column number of Alt State Type
	int nc_alt_st_descr; //!< Column number of Alt State Description
	int nc_res_id;       //!< Column number of Residue ID
	int nc_pka_val;      //!< Column number of current  pKa values
	int nc_std_pka_val;  //!< Column number of standard pKa values
	int nc_alt_st_active; //!< Column number for Alt State active flag

	void SetColumns();
	void FillAltStateList();

	void OnChangeProp(wxCommandEvent& event);
	void OnClose(wxCloseEvent& event);
	void OnCloseBtn(wxCommandEvent& event);
	void OnStdPK(wxCommandEvent& event);
	void OnCalcPKa(wxCommandEvent& event);
	void OnCalcPKPNP(wxCommandEvent& event);
	void OnSetAltInactive(wxCommandEvent& event);
	void OnSetAltActive(wxCommandEvent& event);
	void OnStdPk1(wxCommandEvent& event);
	void OnStdPKEP(wxCommandEvent& event);
	void OnStdPKEP_1(wxCommandEvent& event);
	void OnUpdateAltStateList(wxCommandEvent& event);
	void OnEndLabelEdit(wxGridEvent& event);

private:
	DECLARE_EVENT_TABLE()
};


#endif // !defined(WX_PROT_REDOX_DLG_H)
