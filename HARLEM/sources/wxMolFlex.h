// wxMolFlex.h: interface for the wxMolFlex class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(WXMOLFLEX_H)
#define WXMOLFLEX_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "hahbhp.h"
#include "hamolset.h"
#include "haflexmod.h"

#include <wx/dialog.h>

class wxMolFlex : public wxDialog  
{
public:
	static int dlg_open;
	wxMolFlex(wxWindow* parent);
	virtual ~wxMolFlex();

	void OnInitDialog(wxInitDialogEvent& event);

	bool TransferDataFromWindow();
	bool TransferDataToWindow();

	void OnClose(wxCloseEvent& event);
	void OnUpdateElementList(wxCommandEvent& event);
	void GetAllAtomGroups();
	void OnSelectAtomGroup(wxCommandEvent& event);
	void OnChangeSelectedHB(wxListEvent& event);
	void OnChangeSelectedPH(wxListEvent& event);
	void UpdateGrids();
	void FindHBondsPH(wxCommandEvent& event);
	void OpenMolMechDlg(wxCommandEvent& event);
//	void OnUpdateChart(wxCommandEvent& event);
	void OnEditGroups(wxCommandEvent& event);
	void OnSaveFirstInputFiles(wxCommandEvent& event);
	void OnRunFirst(wxCommandEvent& event);
	void OnComputeDutyCycle(wxCommandEvent& event);
	void OnChooseMDTrajectory(wxCommandEvent& event);
	void OnDeleteSelectedHB(wxCommandEvent& event);
	void OnDeleteSelectedPH(wxCommandEvent& event);

private:
  DECLARE_EVENT_TABLE();

  MolSet* pmset;
  HaFlexMod* myHaFlexMod;

};

#endif // !defined(WXMOLFLEX_H)
