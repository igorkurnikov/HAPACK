// wxMolFlex.cpp: implementation of the wxMolFlex class.
//
//////////////////////////////////////////////////////////////////////
#include <mpi.h>
#include "wx/wxprec.h"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include "wx/valgen.h"
#include "wx/grid.h"
#include "wx/dir.h"
#include "wx/filename.h"
#include "wx/listctrl.h"
#include "wx/tokenzr.h"
#include "wx/notebook.h"
#include "wx/process.h"

#ifdef HAOGL
#include <wx/glcanvas.h>
#endif

#include "harlemapp_wx.h"
#include "ctrl_wx.h"
#include "ha_wx_aux_1.h"
#include "mm_dialogs_wx.h"
#include "halocorb.h"
#include "haqchem.h"
#include "etcoupl.h"
#include "hagaussian.h"
#include "hazindo.h"
#include "hadalton.h"
#include "harlemapp.h"
#include "hamolset.h"
#include "hamolview.h"
#include "hamolecule.h"
#include "moleditor.h"
#include "hamolmech.h"
#include "haintermol.h"
#include "haresdb.h"
#include "nuclacidmod.h"
#include "hascattermod.h"

#include "dialogs_wx_1.h"
//#include "ha_wx_res_wdr.h"


#if defined(__WXMSW__)
// undefined MFC MACROS conflicting WX functions
#ifdef DrawText
#undef DrawText
#endif

#ifdef StartDoc
#undef StartDoc
#endif

#ifdef GetCharWidth
#undef GetCharWidth
#endif

#ifdef FindWindow
#undef FindWindow
#endif

#endif

#include "wxMolFlex.h"
#include "ha_wx_res_molflex_wdr.h"


#if defined(__WXMSW__)
// undefined MFC MACROS conflicting WX functions
#ifdef DrawText
#undef DrawText
#endif

#ifdef StartDoc
#undef StartDoc
#endif

#ifdef GetCharWidth
#undef GetCharWidth
#endif

#ifdef FindWindow
#undef FindWindow
#endif

#endif

#include "haflexmod.h"
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
int wxMolFlex::dlg_open = FALSE;
wxMolFlex::wxMolFlex(wxWindow* parent)
:wxDialog( parent, -1, "Molecular Flexibility with FIRST", wxDefaultPosition, wxDefaultSize, 
					 wxDEFAULT_DIALOG_STYLE|wxRESIZE_BORDER )
{
	molFlexDlg(this,TRUE);
	pmset = GetCurMolSet();

	myHaFlexMod = pmset->GetFlexMod(true);
}

void wxMolFlex::OnInitDialog(wxInitDialogEvent& event)
{
	dlg_open = TRUE;

	wxTextCtrl* edit_ctrl;

	edit_ctrl= (wxTextCtrl*) FindWindow(ID_HBOND_CUTOFF_DIST);
	edit_ctrl->SetValidator( wxDoubleValidator(&myHaFlexMod->p_mol_editor->max_da_dist_no_acc,"%8.3f"));

	edit_ctrl= (wxTextCtrl*) FindWindow(ID_HB_MAX_ANGLE);
	edit_ctrl->SetValidator( wxDoubleValidator(&myHaFlexMod->p_mol_editor->max_hda_angle,"%8.3f",RAD_TO_DEG));

	edit_ctrl= (wxTextCtrl*) FindWindow(ID_PH_CUTOFF_DIST);
	edit_ctrl->SetValidator( wxDoubleValidator(&myHaFlexMod->hydrophob_dist_cutoff,"%8.3f"));

    edit_ctrl= (wxTextCtrl*) FindWindow(ID_HB_DUTY_CYCLE_CUTOFF);
	edit_ctrl->SetValidator( wxDoubleValidator(&myHaFlexMod->hb_duty_cycle_cutoff,"%7.1f"));

	edit_ctrl= (wxTextCtrl*) FindWindow(ID_HB_ENERGY_CUTOFF);
	edit_ctrl->SetValidator( wxDoubleValidator(&myHaFlexMod->hb_energy_cutoff,"%8.3f"));

    edit_ctrl= (wxTextCtrl*) FindWindow(ID_HPT_DUTY_CYCLE_CUTOFF);
	edit_ctrl->SetValidator( wxDoubleValidator(&myHaFlexMod->hpt_duty_cycle_cutoff,"%7.1f"));

	edit_ctrl= (wxTextCtrl*) FindWindow(ID_MD_TRAJ_NAME);
	edit_ctrl->SetValidator( StdStringValidator(&myHaFlexMod->md_traj_file_name) );

	wxListCtrl* hbGrid = (wxListCtrl*) FindWindow(ID_HB_LIST);
	hbGrid->ClearAll();

	hbGrid->InsertColumn(0,"Donor");
	hbGrid->SetColumnWidth(0, 120);
	hbGrid->InsertColumn(1,"Acceptor");
	hbGrid->SetColumnWidth(1, 120);
	hbGrid->InsertColumn(2,"Distance");
	hbGrid->SetColumnWidth(2, 80);
	hbGrid->InsertColumn(3, wxString("Angle"));
	hbGrid->SetColumnWidth(3, 80);
	hbGrid->InsertColumn(4, wxString("Duty Cycle"));
	hbGrid->SetColumnWidth(4, 80);

	wxListCtrl* phGrid = (wxListCtrl*) FindWindow(ID_PH_LIST);
	phGrid->ClearAll();

	phGrid->InsertColumn(0, wxString("Atom 1"));
	phGrid->SetColumnWidth(0, 120);
	phGrid->InsertColumn(1, wxString("Atom 2"));
	phGrid->SetColumnWidth(1, 120);
	phGrid->InsertColumn(2, wxString("Distance"));
	phGrid->SetColumnWidth(2, 80);
	phGrid->InsertColumn(3, wxString("Duty Cycle"));
	phGrid->SetColumnWidth(3, 80);

	event.Skip();
}

bool
wxMolFlex::TransferDataFromWindow()
{	
	bool bres = wxDialog::TransferDataFromWindow();
	return bres;
}

bool
wxMolFlex::TransferDataToWindow()
{	
	bool bres = wxDialog::TransferDataToWindow();
	UpdateGrids();
	return bres;
}

void wxMolFlex::UpdateGrids()
{
	int i;

	int n_hpt = myHaFlexMod->HPTArray.size();
	int nhb = myHaFlexMod->HBArray.size();

	wxListCtrl* hbGrid = (wxListCtrl*) FindWindow(ID_HB_LIST);
	hbGrid->DeleteAllItems();

	for( i = 0; i < nhb; i++) 
	{
		HBondAvg* p_hb =  myHaFlexMod->HBArray[i];
		HaAtom* p_don = p_hb->GetDonorAtom();
		HaAtom* p_acc = p_hb->GetAcceptorAtom();
		Vec3D h_coord;
		int h_found = p_hb->GetHCoord(h_coord);

		if(!p_don->Selected() || !p_acc->Selected()) continue;

		double dist = Vec3D::CalcDistance(p_hb->GetDonorAtom(),p_hb->GetAcceptorAtom(),ANGSTROM_U);
		if( p_hb->GetAvgDistance() > 0.1 ) dist = p_hb->GetAvgDistance();

		double angle = 0.0;
		if(h_found) angle = Vec3D::CalcAngle(&h_coord,p_don,p_acc)*RAD_TO_DEG;
		
		double duty_cycle = p_hb->GetDutyCycle(); 
		
		hbGrid->InsertItem(i, (p_don->GetRef()).c_str() );
		hbGrid->SetItemPtrData(i,(wxUIntPtr)p_hb);
		hbGrid->SetItem(i,1, (p_acc->GetRef()).c_str() );
		hbGrid->SetItem(i,2,wxString::Format("%8.3f", dist ));
		hbGrid->SetItem(i,3,wxString::Format("%8.3f", angle ));
		hbGrid->SetItem(i,4,wxString::Format("%7.3f", duty_cycle ));
	}
	
	wxListCtrl* phGrid = (wxListCtrl*) FindWindow(ID_PH_LIST);
	phGrid->DeleteAllItems();
	
	for(i = 0; i < n_hpt; i++) 
	{
		HaHydrophobicTether* p_hpt = myHaFlexMod->HPTArray[i];

		HaAtom* aptr1 = p_hpt->GetFirstAtom();
		HaAtom* aptr2 = p_hpt->GetSecondAtom();
		if(!aptr1->Selected() || !aptr2->Selected()) continue;

		double duty_cycle = p_hpt->GetDutyCycle();
		double dist = Vec3D::CalcDistance(aptr1,aptr2,ANGSTROM_U);
		if( p_hpt->GetAvgDistance() > 0.1 ) dist = p_hpt->GetAvgDistance();

		phGrid->InsertItem(i, (aptr1->GetRef()).c_str() );
		phGrid->SetItemPtrData(i,(wxUIntPtr)p_hpt);
		phGrid->SetItem(i,1, (aptr2->GetRef()).c_str() );
		phGrid->SetItem(i,2,wxString::Format("%8.3f", dist ));		
		phGrid->SetItem(i,3,wxString::Format("%8.3f", duty_cycle ));
	}
}

wxMolFlex::~wxMolFlex()
{

}

void wxMolFlex::OnClose(wxCloseEvent& event)
{
	dlg_open = FALSE;
	event.Skip();
}

void wxMolFlex::GetAllAtomGroups()
{
	MolSet* pmset = GetCurMolSet();
	// Have to ask  igor about methods that access the groups
	wxControlWithItems* atomsi = (wxControlWithItems *) FindWindow(ID_ATMGRP1_SELECT);
	wxControlWithItems* atomsj = (wxControlWithItems *) FindWindow(ID_ATMGRP2_SELECT);
	
	//Second, extract all groups defined in the file itself
	AtomGroupIteratorMolSet at_arr_itr(pmset);
	for(AtomGroup* myGroup = at_arr_itr.GetFirst(); myGroup != NULL; myGroup = at_arr_itr.GetNext())
		if(myGroup != NULL)
		{
			std::string str = myGroup->GetID();
			PrintLog("Name of Group: %s\n", str.c_str());
			if(atomsi->IsEnabled() == TRUE)
				atomsi->Append(str.c_str());
			if(atomsj->IsEnabled() == TRUE)
				atomsj->Append(str.c_str());
		}
	atomsi->Refresh();
	atomsj->Refresh();
}

void wxMolFlex::OpenMolMechDlg(wxCommandEvent& event)
{
	HaMolMechMod* myNewPtr = new HaMolMechMod(pmset);
	MolMechDlgWX* newWindow = new MolMechDlgWX( myNewPtr, this->GetParent() );
	newWindow->Show();
}


void wxMolFlex::OnSelectAtomGroup(wxCommandEvent& event)
{
	MolSet* pmset = myHaFlexMod->GetMolSet();
	wxControlWithItems* atmGrpComboBox = (wxControlWithItems*) FindWindow(ID_ATMGRP1_SELECT);

 	wxString mySelection = atmGrpComboBox->GetStringSelection();
	int myIntChoice = atmGrpComboBox->GetSelection();
	PrintLog("Selection: %s\n", mySelection.ToStdString().c_str());
	
/*	if(mySelection == "") // Take care of dirty business
		PrintLog("Nothing to select\n");
	else
	{
		AtomGroup* atg = pmset->GetAtomGroupByID(mySelection.c_str());
		if( atg != NULL) 
		{
			hbList = myHaFlexMod->getHydrogenBondsList(atg);
			phList = myHaFlexMod->getHydrophobesList(atg);
			UpdateGrids();
		}
	}
*/
}

void wxMolFlex::FindHBondsPH(wxCommandEvent& event)
{
	TransferDataFromWindow();
	myHaFlexMod->FindHydrogenBonds();
	myHaFlexMod->FindHydrophobicContacts();
	TransferDataToWindow();
}

//void wxMolFlex::OnUpdateChart(wxCommandEvent& event)
//{
//	MolSet* pmset = myHaFlexMod->GetMolSet();
//	wxControlWithItems* atmGrpComboBox = (wxControlWithItems*) FindWindow(ID_ATMGRP1_SELECT);

// 	wxString mySelection = atmGrpComboBox->GetStringSelection();
//	int myIntChoice = atmGrpComboBox->GetSelection();
//	PrintLog("Selection: %s\n", mySelection.c_str());
//	
//	if(mySelection == "") // Take care of dirty business
//		PrintLog("Nothing to select\n");
//	else
//	{
//		AtomGroup* atg = pmset->GetAtomGroupByID(mySelection.c_str());
//		if( atg != NULL) 
//		{
//			hbList = myHaFlexMod->getHydrogenBondsList(atg);
//			phList = myHaFlexMod->getHydrophobesList(atg);
//			UpdateGrids();
//		}
//	}
//
//}

void wxMolFlex::OnEditGroups(wxCommandEvent& event)
{
	EditGroupsDlg* newWindow = new EditGroupsDlg(pmset, 2, this->GetParent() );
	newWindow->Show();
}

void wxMolFlex::OnSaveFirstInputFiles(wxCommandEvent& event)
{
	TransferDataFromWindow();
	if(myHaFlexMod->HBArray.empty() || myHaFlexMod->HPTArray.empty())
		FindHBondsPH(event);
	myHaFlexMod->SaveFirstInpFiles();
}

void wxMolFlex::OnRunFirst(wxCommandEvent& event)
{
	MolSet* pmset = myHaFlexMod->GetMolSet();

	std::string cmdLine = "C:\\HARLEM\\FIRST.exe";
    cmdLine += " -non ";
    cmdLine += " -hbin ";
    cmdLine += " -phin ";
	cmdLine += myHaFlexMod->struct_file_name.c_str(); 
	PrintLog("cmdLine: %s\n", cmdLine.c_str());
	long myStatus = wxExecute(cmdLine, wxEXEC_SYNC, NULL);
		
	myHaFlexMod->ReadFirstDataFile();

	HaMolView* pview = pmset->GetActiveMolView();
	pview->RigidClusterColourAttrib();
	pview->UpdateThisView(RFRefresh);
}

void wxMolFlex::OnChangeSelectedHB(wxListEvent& event)
{
	MolSet* pmset = myHaFlexMod->GetMolSet();
	long idx = -1;
	wxListCtrl* hbGrid = (wxListCtrl*) FindWindow(ID_HB_LIST);
	PrintLog("Selected HB indexes: \n)");
	myHaFlexMod->selected_hb.clear();
	for(;;)
	{
		idx = hbGrid->GetNextItem(idx, wxLIST_NEXT_ALL,
                                       wxLIST_STATE_SELECTED);
		if( idx == -1) break;
		PrintLog(" %d \n", idx);
		HBondAvg* p_hb = (HBondAvg*) hbGrid->GetItemData(idx);
		myHaFlexMod->selected_hb.push_back( p_hb );
	}
	pmset->RefreshAllViews();
}

void
wxMolFlex::OnChangeSelectedPH(wxListEvent& event)
{
	MolSet* pmset = myHaFlexMod->GetMolSet();
	long idx = -1;
	wxListCtrl* phGrid = (wxListCtrl*) FindWindow(ID_PH_LIST);
	PrintLog("Selected HPT indexes: \n)");
	myHaFlexMod->selected_hpt.clear();
	for(;;)
	{
		idx = phGrid->GetNextItem(idx, wxLIST_NEXT_ALL,
                                       wxLIST_STATE_SELECTED);
		if( idx == -1) break;
		PrintLog(" %d \n", idx);
		HaHydrophobicTether* p_hpt = (HaHydrophobicTether*) phGrid->GetItemData(idx);
		myHaFlexMod->selected_hpt.push_back( p_hpt );
	}
	pmset->RefreshAllViews();
}

void wxMolFlex::OnDeleteSelectedHB(wxCommandEvent& event)
{
	myHaFlexMod->DeleteSelectedHBonds();
	
	MolSet* pmset = myHaFlexMod->GetMolSet();
	pmset->RefreshAllViews();
	TransferDataToWindow();
}

void wxMolFlex::OnDeleteSelectedPH(wxCommandEvent& event)
{
	myHaFlexMod->DeleteSelectedHPTethers();
	
	MolSet* pmset = myHaFlexMod->GetMolSet();
	pmset->RefreshAllViews();
	TransferDataToWindow();
}

void wxMolFlex::OnComputeDutyCycle(wxCommandEvent& event)
{
	TransferDataFromWindow();
	myHaFlexMod->ComputeBondsDutyCycleAlongMD();
	pmset->RefreshAllViews();
	TransferDataToWindow();
}

void wxMolFlex::OnChooseMDTrajectory(wxCommandEvent& event)
{
	wxString traj_file = wxFileSelector("Choose MD Trajectory File");
	myHaFlexMod->md_traj_file_name = traj_file.c_str();
	TransferDataToWindow();
}
 

BEGIN_EVENT_TABLE(wxMolFlex, wxDialog)
    EVT_INIT_DIALOG( wxMolFlex::OnInitDialog )
	EVT_BUTTON(ID_MD_STATS, wxMolFlex::OpenMolMechDlg)
	EVT_BUTTON(ID_FIND_HBOND_AND_PH,wxMolFlex::FindHBondsPH)
//	EVT_BUTTON(ID_UPDATE_CHART, wxMolFlex::OnUpdateChart)
	EVT_BUTTON(ID_EDIT_GRPS, wxMolFlex::OnEditGroups)
	EVT_BUTTON(ID_SAVE_FIRST_INPUT_FILES, wxMolFlex::OnSaveFirstInputFiles)
	EVT_BUTTON(ID_RUN_FIRST, wxMolFlex::OnRunFirst)
	EVT_BUTTON(ID_COMPUTE_DUTY_CYCLE, wxMolFlex::OnComputeDutyCycle)
	EVT_BUTTON(ID_DEL_SELECT_HB, wxMolFlex::OnDeleteSelectedHB)
	EVT_BUTTON(ID_DEL_SELECT_PH, wxMolFlex::OnDeleteSelectedPH)
	EVT_BUTTON(ID_CHOOSE_MD_TRAJ, wxMolFlex::OnChooseMDTrajectory)
	EVT_LIST_ITEM_SELECTED(ID_HB_LIST,wxMolFlex::OnChangeSelectedHB )
    EVT_LIST_ITEM_DESELECTED(ID_HB_LIST,wxMolFlex::OnChangeSelectedHB )
	EVT_LIST_ITEM_SELECTED(ID_PH_LIST,wxMolFlex::OnChangeSelectedPH )
    EVT_LIST_ITEM_DESELECTED(ID_PH_LIST,wxMolFlex::OnChangeSelectedPH )
	// EVT_CHOICE(ID_ATMGRP1_SELECT, wxMolFlex::OnSelectAtomGroup)
	// EVT_CHOICE(ID_ATMGRP2_SELECT, wxMolFlex::OnSelectAtomGroup)
	EVT_CLOSE   (wxMolFlex::OnClose)
END_EVENT_TABLE()
