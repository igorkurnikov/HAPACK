/*! \file wx_prot_redox_dlg.cpp

    Dialogs for Protonation and Redox Equilibrium Calculations
 
    \author Igor Kurnikov  
    \date 2010-
*/

#include <mpi.h>

#include <wx/wx.h>
#include <wx/grid.h>

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include "ha_wx_aux_1.h"
#include "ha_wx_res_wdr.h"

#include "haatgroup.h"
#include "hamolset.h"
#include "moleditor.h"
#include "protonredox.h"
#include "wx_prot_redox_dlg.h"


//////////////////////////////////////////////////////////////////////
// ProtonRedoxDlg dialog

int  ProtonRedoxDlg::dlg_open = FALSE;

ProtonRedoxDlg::ProtonRedoxDlg(ProtonRedoxMod* p_prot_rdx_mod_new, wxWindow* parent):
wxFrame( parent, -1, "Protonation and Redox Equilibrium")
{
	p_prot_rdx_mod = p_prot_rdx_mod_new;
	pmset= p_prot_rdx_mod->GetMolSet();
	p_mol_editor = pmset->GetMolEditor(true);

	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);

	wxMenuBar* prot_redox_menu_bar = prot_redox_menu();
    this->SetMenuBar(prot_redox_menu_bar);    

	prot_redox_dlg(this,TRUE);
	p_alt_st_grid = (wxGrid*) FindWindow( IDC_ALT_ST_GRID );
	OnInitDialog();     
}

ProtonRedoxDlg::~ProtonRedoxDlg()
{
    dlg_open = FALSE;
}


void ProtonRedoxDlg::OnInitDialog()
{
	DDX_Choice_HaEnum(this,IDC_PKA_MULTI_SITE_METHOD,&p_prot_rdx_mod->multi_site_pop_method);
	DDX_Text_double (this,IDC_CURRENT_PH,p_prot_rdx_mod->ph,"%8.3f");
	DDX_Text_double (this,IDC_CURRENT_E0,p_prot_rdx_mod->e0,"%8.3f");
	DDX_Text_double (this,IDC_PH_MIN,p_prot_rdx_mod->ph_min,"%8.3f");
	DDX_Text_double (this,IDC_PH_MAX,p_prot_rdx_mod->ph_max,"%8.3f");
	DDX_Text_double (this,IDC_PH_STEP,p_prot_rdx_mod->ph_step,"%8.3f");
	DDX_Text_double (this,IDC_E0_MIN,p_prot_rdx_mod->e0_min,"%8.3f");
	DDX_Text_double (this,IDC_E0_MAX,p_prot_rdx_mod->e0_max,"%8.3f");
	DDX_Text_double (this,IDC_E0_STEP,p_prot_rdx_mod->e0_step,"%8.3f");
	
	wxCheckBox* check_box = (wxCheckBox*) FindWindow( IDC_READ_INTER_MATRIX );
	check_box->SetValidator( IntCheckBoxValidator(&p_prot_rdx_mod->read_alt_st_inter) );
	
	check_box = (wxCheckBox*) FindWindow( IDC_SAVE_INTER_MATRIX );
	check_box->SetValidator( IntCheckBoxValidator(&p_prot_rdx_mod->save_alt_st_inter) );

	DDX_Text_int (this,IDC_NUM_MC_STEPS,p_prot_rdx_mod->n_mc_cyc);

	SetColumns();	
	TransferDataToWindow();
	dlg_open = TRUE;
}

void ProtonRedoxDlg::SetColumns()
{
//	int tot_width = 550;
//	int height = 400;

	int tot_width;
	int height;

	p_alt_st_grid->GetClientSize(&tot_width, &height);

	int num_cols = 6;
	int col_width = (tot_width*2 /3)/num_cols;

	p_alt_st_grid->SetColLabelSize(21);
    
	int ncol_act = p_alt_st_grid->GetNumberCols();
	if(ncol_act != num_cols)
	{
		p_alt_st_grid->DeleteCols(0,ncol_act);
		p_alt_st_grid->InsertCols(0,num_cols);
	}

	int icol = -1;
     
	nc_alt_st_idx = ++icol;
	p_alt_st_grid->SetColLabelValue(nc_alt_st_idx, "Idx");
	nc_alt_st_type = ++icol;
	p_alt_st_grid->SetColLabelValue(nc_alt_st_type, "Type");
	nc_alt_st_descr = ++icol;
	p_alt_st_grid->SetColLabelValue(nc_alt_st_descr, "Description");
//	nc_res_id = ++icol;
//	p_alt_st_grid->SetColLabelValue(nc_res_id, "Residue ID");
	nc_pka_val = ++icol;
	p_alt_st_grid->SetColLabelValue(nc_pka_val, "pKa Val");
	nc_std_pka_val = ++icol;
	p_alt_st_grid->SetColLabelValue(nc_std_pka_val, "Std pKa Val");
	nc_alt_st_active = ++icol;
	p_alt_st_grid->SetColLabelValue(nc_alt_st_active, "Active");
}

void ProtonRedoxDlg::FillAltStateList()
{
	char buf[256];

	p_alt_st_grid->ClearGrid();

	if(pmset == NULL) return;
	
	HaResidue* pres;

	int ncol = p_alt_st_grid->GetNumberCols();
	int nrow = p_alt_st_grid->GetNumberRows();

	int nst = p_prot_rdx_mod->alt_chem_states.size();

	if(nrow != nst )
	{
		if( nrow > 0) p_alt_st_grid->DeleteRows(0, nrow);
        p_alt_st_grid->AppendRows(nst);
	}

	int max_lbl = 0;
	int ist;

	for(ist = 0; ist < nst; ist++)
	{
		AltChemState* p_alt_st =  p_prot_rdx_mod->alt_chem_states[ist];
		AtomGroup* p_grp = p_alt_st->GetHostAtomGroup();
		
		HaResidue* pres = dynamic_cast<HaResidue*>(p_grp);
		std::string res_lbl; 

		if( pres ) res_lbl = pres->GetRef();
		else res_lbl = p_grp->GetID();
		
		if(res_lbl.size() > max_lbl) max_lbl = res_lbl.size();
	
		p_alt_st_grid->SetRowLabelValue(ist, res_lbl.c_str() ); 
		p_alt_st_grid->SetCellValue(ist, nc_alt_st_idx, wxString::Format("%d",ist));
		p_alt_st_grid->SetCellValue(ist, nc_alt_st_type, p_alt_st->alt_state_type.label() );
		p_alt_st_grid->SetCellValue(ist, nc_alt_st_descr, p_alt_st->id.c_str());
		p_alt_st_grid->SetCellValue(ist, nc_pka_val, wxString::Format("%6.2f",p_alt_st->pk));
		p_alt_st_grid->SetCellValue(ist, nc_std_pka_val, wxString::Format("%6.2f",p_alt_st->std_pk));
		p_alt_st_grid->SetCellValue(ist, nc_alt_st_active, wxString::Format("%d",p_alt_st->active_flag));
	}
	p_alt_st_grid->SetRowLabelSize((int)(max_lbl*8.5));
	p_alt_st_grid->AutoSizeColumns();
}

bool ProtonRedoxDlg::TransferDataToWindow()
{
	FillAltStateList();
	return wxFrame::TransferDataToWindow();

}

BEGIN_EVENT_TABLE(ProtonRedoxDlg, wxFrame)
	EVT_BUTTON(IDC_EDTRES_UPDATE_RESLIST, ProtonRedoxDlg::OnUpdateAltStateList)
	EVT_BUTTON(IDC_CLOSE, ProtonRedoxDlg::OnCloseBtn)
	EVT_CLOSE (ProtonRedoxDlg::OnClose)
	EVT_MENU(IDC_PROT_RDX_CALC_PKA,ProtonRedoxDlg::OnCalcPKa)
	EVT_MENU(IDC_RESPAR_STD_PK,ProtonRedoxDlg::OnStdPK)
	EVT_MENU(IDC_SET_ALTST_INACTIVE, ProtonRedoxDlg::OnSetAltInactive)
	EVT_MENU(IDC_SET_ALTST_ACTIVE, ProtonRedoxDlg::OnSetAltActive)
	EVT_MENU(IDC_RESPAR_STD_PK_1, ProtonRedoxDlg::OnStdPk1)
	EVT_MENU(IDC_RESPAR_STD_PK_EP, ProtonRedoxDlg::OnStdPKEP)
	EVT_MENU(IDC_STD_PK_EP_1, ProtonRedoxDlg::OnStdPKEP_1)
	EVT_GRID_CELL_CHANGE( ProtonRedoxDlg::OnEndLabelEdit)
END_EVENT_TABLE()

void ProtonRedoxDlg::OnUpdateAltStateList(wxCommandEvent& event)
{
	FillAltStateList();
}

void ProtonRedoxDlg::OnCloseBtn(wxCommandEvent& event)
{
	Close();
}

void ProtonRedoxDlg::OnClose(wxCloseEvent& event)
{
	ProtonRedoxDlg::dlg_open = FALSE;
	event.Skip();

}

void ProtonRedoxDlg::OnStdPK(wxCommandEvent& event)
{
	p_prot_rdx_mod->set_std_redox_pot = FALSE;
	p_prot_rdx_mod->SetStdPKa();
	FillAltStateList();
}

void ProtonRedoxDlg::OnStdPKEP(wxCommandEvent& event) 
{
	p_prot_rdx_mod->set_std_redox_pot = TRUE;
	p_prot_rdx_mod->SetStdPKa();
	FillAltStateList();	
}

void ProtonRedoxDlg::OnStdPk1(wxCommandEvent& event) 
{
	p_prot_rdx_mod->set_std_redox_pot = FALSE;
	p_prot_rdx_mod->SetStdPKa_G1();
	FillAltStateList();	
}

void ProtonRedoxDlg::OnStdPKEP_1(wxCommandEvent& event) 
{
	p_prot_rdx_mod->set_std_redox_pot = TRUE;
	p_prot_rdx_mod->SetStdPKa_G1();
	FillAltStateList();	
}

void ProtonRedoxDlg::OnCalcPKa(wxCommandEvent& event)
{
	TransferDataFromWindow();
	p_prot_rdx_mod->CalcPKaForSelection();
	TransferDataToWindow();
}
void ProtonRedoxDlg::OnCalcPKPNP(wxCommandEvent& event)
{
	p_prot_rdx_mod->CalcPKaForSelection(true);
	TransferDataToWindow();
}


void ProtonRedoxDlg::OnChangeProp(wxCommandEvent& event)
{
	TransferDataFromWindow();
	SetColumns();
	FillAltStateList();
}

void ProtonRedoxDlg::OnSetAltInactive(wxCommandEvent& event) 
{
	p_prot_rdx_mod->SetAltStatesActive(FALSE);
	TransferDataToWindow();
}

void ProtonRedoxDlg::OnSetAltActive(wxCommandEvent& event) 
{
	p_prot_rdx_mod->SetAltStatesActive(TRUE);
	TransferDataToWindow();
}

void ProtonRedoxDlg::OnEndLabelEdit(wxGridEvent& event)
{
	PrintLog( " Enter ProtonRedoxDlg::OnEndLabelEdit() \n");
	
	int icol = event.GetCol();
	int irow = event.GetRow();
	
	std::string str_val = (p_alt_st_grid->GetCellValue(irow,icol)).ToStdString();

	boost::trim(str_val);
	
	if( irow >= p_prot_rdx_mod->alt_chem_states.size() )
	{
		PrintLog(" Error in ProtonRedoxDlg::OnEndLabelEdit() \n");
		PrintLog(" Invalid index %d  in the array of alternative chemical states \n",irow);
		return;
	}

	AltChemState* p_st = p_prot_rdx_mod->alt_chem_states[irow];

	double dval;
	int ival;
	if(icol == nc_pka_val)
	{
		try
		{
			dval = boost::lexical_cast<double>(str_val);
			p_st->pk = dval;
		}
		catch(boost::bad_lexical_cast&) 
		{
			PrintLog(" Error in ProtonRedoxDlg::OnEndLabelEdit() \n");
			PrintLog(" Invalid String for PKa value %s \n",str_val.c_str());
		}
		p_alt_st_grid->SetCellValue(irow, icol, wxString::Format("%6.2f",p_st->pk));
	}

	if(icol == nc_std_pka_val)
	{
		try
		{
			dval = boost::lexical_cast<double>(str_val);
			p_st->std_pk = dval;
		}
		catch(boost::bad_lexical_cast&) 
		{
			PrintLog(" Error in ProtonRedoxDlg::OnEndLabelEdit() \n");
			PrintLog(" Invalid String for STD PKa value %s \n",str_val.c_str());
		}
		p_alt_st_grid->SetCellValue(irow, icol, wxString::Format("%6.2f",p_st->std_pk));
	}
	if(icol == nc_alt_st_active)
	{
		try
		{
			ival = boost::lexical_cast<int>(str_val);
			p_st->active_flag = ival;
		}
		catch(boost::bad_lexical_cast&) 
		{
			PrintLog(" Error in ProtonRedoxDlg::OnEndLabelEdit() \n");
			PrintLog(" Invalid String for inactive flag value %s \n",str_val.c_str());
		}
		p_alt_st_grid->SetCellValue(irow, icol, wxString::Format("%d",p_st->active_flag));
	}
}
