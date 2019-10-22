/*! \file wxMolED.h

    Dialogs for Essential Gynamics and clustering analysis
 
    \author Igor Kurnikov  
    \date 2011-
*/

#if !defined(WX_MOL_ED_H)
#define WX_MOL_ED_H

#include "hamolset.h"

#include <wx/dialog.h>

class wxMolED: public wxFrame
{
public:
	wxMolED( CollectCrdAnalMod* p_cluster_mod_new, wxWindow* parent);
	virtual ~wxMolED();

	static int dlg_open;

	virtual bool TransferDataToWindow();
	virtual bool TransferDataFromWindow();

	void OnInitDialog();

	void OnTestGraph(wxCommandEvent& event);
	void OnChooseMDTraj(wxCommandEvent& event);
	void OnChoosePlatoTraj(wxCommandEvent& event);
	void OnChooseAtGrpToAnal(wxCommandEvent& event);
	void OnGenPlatoTraj(wxCommandEvent& event);
	void OnSavePlatoInputFiles(wxCommandEvent& event);
	void OnRunPlato(wxCommandEvent& event);
	void OnChoosePlatoOutputFile(wxCommandEvent& event);
	void OnLoadPlatoOutputFile(wxCommandEvent& event);
	void OnRunVecAnimation(wxCommandEvent& event);
	void OnStopVecAnimation(wxCommandEvent& event);
	void OnProjTraj(wxCommandEvent& event);
	void OnShiftAlongVec(wxCommandEvent& event);
	void OnPlotTimeProj(wxCommandEvent& event);
	void OnPlotProj2vs1(wxCommandEvent& event);
	void OnTopFileBrowse(wxCommandEvent& event);
	void OnEigFileBrowse(wxCommandEvent& event);
	void OnPpjFileBrowse(wxCommandEvent& event);
	void OnNextButtonClick(wxCommandEvent& event);
	void OnClose(wxCloseEvent& event);

protected:
	wxString topFileName;	//!< Topology File name
	wxString topFileDir;
	wxString eigFileName;	//!< Eigenvalue / vectors file name
	wxString eigFileDir;
	wxString ppjFileName;	//!< Projection file name
	wxString ppjFileDir;
	wxWindow* parent;
	
	MolSet* pmset;
	CollectCrdAnalMod* p_ccrd_mod;
	HaMolView* p_anim_view; //!< Molecular view to run eigenvector animation 

	double shift_vec_val; //!< Value to shift structure along selected eigen vector 
	wxArrayInt sel_arr;   //!< Array of indexes of selected eigenvectors  

private:
	DECLARE_EVENT_TABLE();
};

#endif


