/*! \file qc_dialogs_wx.cpp

    Dialogs for Quantum Chemical Calculations
 
    \author Igor Kurnikov  
    \date 2011-
*/

#define QC_DIALOGS_WX_CPP

#include <mpi.h>

#include "vec3d.h"
#include "g94_globals.h"
#include "g94_protos.h"

#include "hastl.h"

#include "wx/wx.h"
#include "wx/notebook.h"
#include "wx/valgen.h"
#include "wx/grid.h"
#include "wx/filename.h"

#include "ctrl_wx.h"
#include "ha_wx_aux_1.h"
#include "ha_wx_res_wdr.h"

#include "canvas3d.h"
#include "dialogs_wx_1.h"
#include "hamolset.h"
#include "hamolview.h"
#include "haatombasis.h"
#include "gaufile.h"
#include "haqchem.h"
#include "etcoupl.h"
#include "qc_dialogs_wx.h"

//----------------------------------------------------------------------------
// EditHuckHamDlg
//----------------------------------------------------------------------------

// WDR: event table for EditHuckHamDlg

BEGIN_EVENT_TABLE(EditHuckHamDlg,wxDialog)
    EVT_GRID_CELL_CHANGED( EditHuckHamDlg::OnGridCellChange )
    EVT_BUTTON( ID_EDH_REFRESH_GRID, EditHuckHamDlg::RefreshView )
END_EVENT_TABLE()

EditHuckHamDlg::EditHuckHamDlg( ETCouplMod* et_mod, wxWindow *parent, wxWindowID id, const wxString &title,
    const wxPoint &position, const wxSize& size, long style ) :
    wxDialog( parent, id, title, position, size, style )
{
    edit_ham_dlg( this, TRUE );

    pmat = &(et_mod->heff_mat);
    ssmat = &(et_mod->ssl);
    pqc_mod = et_mod->GetQCMod();
}

bool EditHuckHamDlg::TransferDataToWindow()
{
    ene_grid = (wxGrid*) FindWindow( ID_HUCK_HAM_GRID );
    ene_grid->SetColFormatFloat(1,9,4);
    bool bres = wxDialog::TransferDataToWindow();
    SetGridData();
    return bres;
}

void EditHuckHamDlg::SetGridData()
{
    int nrow_grid = ene_grid->GetNumberRows();
    int ncol_grid = ene_grid->GetNumberCols();
    if(pmat == NULL)
    {
         if( nrow_grid > 0 ) ene_grid->DeleteRows(0,nrow_grid);
         return;
    }
    int nrow_mat = pmat->num_rows();
    if(nrow_grid != nrow_mat )
    {
        if( nrow_grid > 0 ) ene_grid->DeleteRows(0,nrow_grid);
        ene_grid->AppendRows(nrow_mat);
    }

//  wxGridTableBase* grid_table = ene_grid->GetTable();
        
    int i;
    ene_grid->SetRowLabelSize(300);
    for( i = 0; i < nrow_mat; i++)
    {
        wxString str; 
        double dval = pmat->GetVal_idx0(i,i);
        str.Printf("%12.6f",dval);
        ene_grid->SetCellValue(i,0,str);
		std::string lbl = pqc_mod->ActBas->GetLabel(i);
        ene_grid->SetRowLabelValue(i,lbl.c_str());
    }
    ene_grid->AutoSizeColumns();
}


// WDR: handler implementations for EditHuckHamDlg
void EditHuckHamDlg::RefreshView( wxCommandEvent& WXUNUSED(event) )
{
    TransferDataToWindow();
}

void EditHuckHamDlg::OnGridCellChange( wxGridEvent &event )
{
    if(pmat == NULL) return;
    int row = event.GetRow();
    wxString sval = ene_grid->GetCellValue(row,0);
    
    double vold = pmat->GetVal_idx0(row,row);
    
    double vnew;
    bool bres = sval.ToDouble(&vnew);
    if(!bres)
    {
        PrintLog(" Error to read value in the cell \n");
        sval.Printf("%12.6f \n",vold);
        ene_grid->SetCellValue(row,0,sval);
        return;
    }

    double shift = vnew - vold;

    if( fabs(shift) > 0.0001) 
    {
        int n = ssmat->num_rows();
        int i,j;
        for(i = 0; i < n; i++)
        {
            for(j = 0; j < n; j++)
            {
                double s1 = ssmat->GetVal_idx0(row,i);
                double s2 = ssmat->GetVal_idx0(row,j);
                double vnew  = pmat->GetVal_idx0(i,j) + s1*s2*shift;
                pmat->SetVal_idx0(i,j,vnew);
            }
        }
//        pmat->SetVal_idx0(row,row,vnew);
        PrintLog(" Value of matrix at %d %d changed to %9.4f \n",row,row,vnew);
        TransferDataToWindow();
    }
}

//////////////////////////////////////////////////////////////////////
// LoadQCDatDlgWX dialog

int LoadQCDatDlgWX::dlg_open = FALSE;

LoadQCDatDlgWX::LoadQCDatDlgWX(HaQCMod* new_qcmod, wxWindow* parent):
wxDialog( parent, -1, "Load Quantum Chemical Data", wxDefaultPosition, wxDefaultSize, 
		   wxDEFAULT_DIALOG_STYLE )
{
	qcmod= new_qcmod;
    dlg_open = TRUE;
	load_qchem_data_dlg(this, true);
}

LoadQCDatDlgWX::~LoadQCDatDlgWX()
{
    dlg_open = FALSE;
}


void LoadQCDatDlgWX::OnInitDialog(wxInitDialogEvent& event)
{
	wxCheckBox* mo_load= (wxCheckBox*) FindWindow( IDC_QCDAT_MO );
	mo_load->SetValue(true);

	wxChoice* file_type= (wxChoice*) FindWindow(IDC_QCDAT_FILE_TYPE);
	file_type->Append("Gaussian binary CHK file");
	file_type->Append("Gaussian formated FCHK file");
    file_type->SetSelection(1);

	if(qcmod != NULL)
	{
		MolSet* pmset = qcmod->GetMolSet();
		if(pmset != NULL)
		{
			wxTextCtrl* dat_fname= (wxTextCtrl*) FindWindow(IDC_QCDAT_DAT_FNAME);
			wxString fname = pmset->GetName();
		    
			int ifile_type = file_type->GetSelection();
			if(ifile_type == 0)
				fname += ".chk";
			else if(ifile_type == 1)
				fname += ".fchk";
			dat_fname->SetValue(fname);
		}
	}
	event.Skip();
}


BEGIN_EVENT_TABLE(LoadQCDatDlgWX,wxDialog)
    EVT_INIT_DIALOG( LoadQCDatDlgWX::OnInitDialog )
    EVT_BUTTON( IDC_QCDAT_LOAD_DATA,    LoadQCDatDlgWX::OnLoadData )
    EVT_CLOSE( LoadQCDatDlgWX::OnClose)
END_EVENT_TABLE()


void LoadQCDatDlgWX::OnLoadData(wxCommandEvent &event)
{
	if(qcmod == NULL) return;

	std::string filetyp_opt;
	wxChoice* file_type= (wxChoice*) FindWindow( IDC_QCDAT_FILE_TYPE );
	
	wxTextCtrl* dat_fname= (wxTextCtrl*) FindWindow(IDC_QCDAT_DAT_FNAME);
	wxString fname = dat_fname->GetValue();
	
	int isel= file_type->GetSelection();
	
	wxBusyCursor wait;
	if(isel == 0)
	{
		wxCheckBox* mo_load= (wxCheckBox*) FindWindow( IDC_QCDAT_MO );
		if(mo_load->IsChecked() )
		{
		    GauFile gfile(fname.c_str(),1,"old");
			int result = gfile.open();
			if(result)
			{
				qcmod->InitMOs(gfile);
				gfile.close();
			}
		}
	}
	if(isel == 1)
	{
		qcmod->LoadDataFromFChk(fname.c_str());
	}
}

void LoadQCDatDlgWX::OnClose(wxCloseEvent& event)
{
	dlg_open = FALSE;
	event.Skip();
}

/////////////////////////////////////////////////////////////////////
// WaveFunAnalDlgWX Dialog:

int WaveFunAnalDlgWX::dlg_open = FALSE;

WaveFunAnalDlgWX::WaveFunAnalDlgWX(HaQCMod* new_phost_qcmod, wxWindow* parent):
wxDialog( parent, -1, "Wave Function Analysis", wxDefaultPosition, wxDefaultSize, 
		   wxDEFAULT_DIALOG_STYLE )
{
	p_qc_mod= new_phost_qcmod;
	wave_fun_anal_dlg(this,TRUE);
}

WaveFunAnalDlgWX::~WaveFunAnalDlgWX()
{
       dlg_open = FALSE;
}

void WaveFunAnalDlgWX::OnInitDialog(wxInitDialogEvent& event)
{
	dlg_open = TRUE;
	if(p_qc_mod == NULL) return;

	wxTextCtrl* edit_mo_isolvl = (wxTextCtrl*) FindWindow( IDC_WFANAL_MO_ISOLVL);
	edit_mo_isolvl->SetValue("0.1");

	wxTextCtrl* edit_grid_size = (wxTextCtrl*) FindWindow( IDC_WFANAL_GRID_SIZE);
	
	edit_grid_size->SetValidator(wxGenericValidator(&p_qc_mod->m_grid_size));

	event.Skip();
}

bool WaveFunAnalDlgWX::TransferDataToWindow()
{
	int i,n;
	wxListBox* mo_list = (wxListBox*) FindWindow( IDC_WFANAL_MO_LIST );
	
	mo_list->Clear();
	if(p_qc_mod == NULL)
		return false;
	if(p_qc_mod->MOene.size() == 0)
	{
		std::cerr << " Error in FillListBoxMOene(): " << std::endl;
		std::cerr << " No Molecular Orbital Energies are set for current QChem Module" << std::endl;
		return false;
	}
	n = p_qc_mod->MOene.size();
	wxString str;
	for(i=1; i <= n; i++)
	{
		double ene= p_qc_mod->MOene(i);
		str.Printf("%4d %16.9f",i,ene);
		mo_list->Append(str);
	}

	wxTextCtrl* text_nae = (wxTextCtrl*) FindWindow(IDC_NALPHA_EL);
	wxTextCtrl* text_nbe = (wxTextCtrl*) FindWindow(IDC_NBETA_EL);

	int nae_tmp = p_qc_mod->GetNumAlphaEl();
	int nbe_tmp = p_qc_mod->GetNumBetaEl();
	str.Printf("%d",nae_tmp);
	text_nae->SetValue(str);
	str.Printf("%d",nbe_tmp);
	text_nbe->SetValue(str);

	return wxDialog::TransferDataToWindow();
}


void WaveFunAnalDlgWX::OnClose(wxCloseEvent& event)
{
	dlg_open = FALSE;
	event.Skip();
}

void WaveFunAnalDlgWX::OnMOrefresh(wxCommandEvent &event)
{
	TransferDataToWindow();
}

void WaveFunAnalDlgWX::OnPlotMO(wxCommandEvent &event)
{
	TransferDataFromWindow();
	wxListBox* mo_list = (wxListBox*) FindWindow( IDC_WFANAL_MO_LIST );
	
	if(p_qc_mod->MO_coef.num_cols() == 0)
	{
		std::cerr << " WaveFunAnalDlgWX::OnPlotMO() " << std::endl;
		std::cerr << " Molecular Orbitals are not set " << std::endl;
		return;
	}
	int idx;
	idx = GetSelMOIdx();
	if(idx < 1) return;

	wxTextCtrl* edit_mo_isolvl = (wxTextCtrl*) FindWindow( IDC_WFANAL_MO_ISOLVL);
	wxString str = edit_mo_isolvl->GetValue();

	double flvl;
	bool bres = str.ToDouble(&flvl);

	if(!bres)
	{
		std::cerr << " WaveFunAnalDlgWX::OnPlotMO() " << std::endl;
		std::cerr << " Error Reading MO Isolevel value " << std::endl;
		return;
	}

	wxBusyCursor wait;
	bool result= p_qc_mod->CreateMOcontour(idx, flvl, p_qc_mod->m_grid_size);
	if(result)
		(p_qc_mod->GetMolSet())->RefreshAllViews(RFRefresh);
}

void WaveFunAnalDlgWX::OnPrintMOCoef(wxCommandEvent &event)
{
	int idx = GetSelMOIdx();

	if(idx < 1) return;

	PrintLog(" MO number: %d \n", idx);
	for(int i=1; i <= p_qc_mod->MO_coef.num_rows(); i++)
	{
		PrintLog("%5d %16.9f \n",i, p_qc_mod->MO_coef(i,idx));
	}

}

int WaveFunAnalDlgWX::GetSelMOIdx()
{
	wxListBox* mo_list = (wxListBox*) FindWindow( IDC_WFANAL_MO_LIST );
	
	if(p_qc_mod->MO_coef.num_cols() == 0)
	{
		PrintLog(" WaveFunAnalDlgWX::OnPrintMOCoef() \n"); 
		PrintLog(" Molecular Orbitals are not set \n"); 
		return -1;
	}
	int idx;
	idx = mo_list->GetSelection();
	if(idx < 0)
	{
		PrintLog(" WaveFunAnalDlgWX::OnPrintMOCoef() \n"); 
		PrintLog(" No current selection in MO list \n"); 
		return -1;
	}
	idx++;

	if( idx > p_qc_mod->MO_coef.num_cols() )
	{
		PrintLog(" WaveFunAnalDlgWX::OnPrintMOCoef() \n"); 
		PrintLog(" Selected MO index is larger than the number of MO \n"); 
		return -1;
	}
	return idx;
}

void WaveFunAnalDlgWX::OnMoveUp(wxCommandEvent &event) 
{
	wxListBox* mo_list = (wxListBox*) FindWindow( IDC_WFANAL_MO_LIST );
	int idx = GetSelMOIdx();
	if(idx < 1) return;

	int nb = p_qc_mod->MO_coef.num_rows();
	int nmo = p_qc_mod->MO_coef.num_cols();
	if(idx >= nmo) return;
	
	int i;
	double tmp;
	for(i=1; i <= nb; i++)
	{
		tmp = p_qc_mod->MO_coef(i,idx);
        p_qc_mod->MO_coef(i,idx) = p_qc_mod->MO_coef(i,idx+1);
        p_qc_mod->MO_coef(i,idx+1) = tmp;
	}
	
	tmp = p_qc_mod->MOene(idx);
	p_qc_mod->MOene(idx) = p_qc_mod->MOene(idx+1);
    p_qc_mod->MOene(idx+1) = tmp;
	TransferDataToWindow();
	mo_list->SetSelection(idx);
}

void WaveFunAnalDlgWX::OnMoveDown(wxCommandEvent &event) 
{
	wxListBox* mo_list = (wxListBox*) FindWindow( IDC_WFANAL_MO_LIST );
	int idx = GetSelMOIdx();
	if(idx <= 1) return;

	int nb = p_qc_mod->MO_coef.num_rows();
	int nmo = p_qc_mod->MO_coef.num_cols();
	
	int i;
	double tmp;
	for(i=1; i <= nb; i++)
	{
		tmp = p_qc_mod->MO_coef(i,idx);
        p_qc_mod->MO_coef(i,idx) = p_qc_mod->MO_coef(i,idx-1);
        p_qc_mod->MO_coef(i,idx-1) = tmp;
	}
	
	tmp = p_qc_mod->MOene(idx);
	p_qc_mod->MOene(idx) = p_qc_mod->MOene(idx-1);
    p_qc_mod->MOene(idx-1) = tmp;
	TransferDataToWindow();
	mo_list->SetSelection(idx-2);
}


BEGIN_EVENT_TABLE(WaveFunAnalDlgWX,wxDialog)
    EVT_INIT_DIALOG( WaveFunAnalDlgWX::OnInitDialog )
    EVT_BUTTON( IDC_MOVE_UP,    WaveFunAnalDlgWX::OnMoveUp )
	EVT_BUTTON( IDC_MOVE_DOWN,    WaveFunAnalDlgWX::OnMoveDown )
	EVT_BUTTON( IDC_WFANAL_PLOT_MO,    WaveFunAnalDlgWX::OnPlotMO )
	EVT_BUTTON( IDC_PRT_MO_COEF,    WaveFunAnalDlgWX::OnPrintMOCoef )
	EVT_BUTTON( IDC_WFANAL_MO_REFRESH,    WaveFunAnalDlgWX::OnMOrefresh )
    EVT_CLOSE( WaveFunAnalDlgWX::OnClose)
END_EVENT_TABLE()

