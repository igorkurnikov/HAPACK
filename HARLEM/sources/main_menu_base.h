/*! \file main_menu_base.h

    Base class for HARLEM main frame menu bar - wxFormBuilder style.
    Builds the complete menu and provides virtual event handlers.

    \author Igor Kurnikov
    \date 2026-
*/

#if !defined(MAIN_MENU_BASE_H)
#define MAIN_MENU_BASE_H

#include <wx/mdi.h>

// New menu IDs not in wxDesigner
enum {
    IDM_LABELS_ATOM_CHARGES = 12001
};

///////////////////////////////////////////////////////////////////////////////
/// Class HaMainMenuBase
///////////////////////////////////////////////////////////////////////////////

class HaMainMenuBase : public wxMDIParentFrame
{
public:
    HaMainMenuBase( wxWindow* parent, wxWindowID id, const wxString& title,
        const wxPoint& pos, const wxSize& size, long style, const wxString& name );

    ~HaMainMenuBase();

protected:
    // Virtual event handlers (override in derived class)

    // File menu
    virtual void OnFileNew( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnFileOpen( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnMolSave( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnMolSaveAs( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoLoadScriptDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoRedirectIODialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void SaveOutputFile( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnSaveImageClipboard( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoCompAccountsDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnInfo( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnFileClose( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnPrint( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnSetup( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnExit( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnOpenMDTrajectory( wxCommandEvent& event ) { event.Skip(); }

    // Edit menu
    virtual void OnSelectAll( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoAtomParamsDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoResidueParamsDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoMolSetParamDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoAtomPropColorDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoEditGroupsDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoCrdSnapshotDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoEditGeomDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoEditMutationMapDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnFindHbonds( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoSolvateDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnWrapIntoUnitCell( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnCenterSolute( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnCenterMolInPBox( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnAddMissingAtoms( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnAddPolarHydrogens( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnAddHHybrid( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnAddHydrogens( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnDelSelAtoms( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnDelOvlpMols( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoEditFragmDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoBuildFilmDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnClearPicked( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnSelAtomsInBoundBox( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnRevertAtomSelection( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnExpandAtomSelectionBonded( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoNuclAcidDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnDescribeSecStruct( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnPrintHBonds( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnSetAlphaHelix( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnShowResDb( wxCommandEvent& event ) { event.Skip(); }

    // Edit->Measure submenu
    virtual void OnShowAtomID( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnMeasureDist( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnMeasureAngle( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnMeasureDihed( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnAddDistMonitor( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnShowAtomLabel( wxCommandEvent& event ) { event.Skip(); }

    // Display->Display Mode menu
    virtual void OnWireFrame( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnSticks( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnSpheres( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnBallStick( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnRibbons( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnStrands( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnCartoons( wxCommandEvent& event ) { event.Skip(); }

    // Display->Colours menu
    virtual void OnMono( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnCPK( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnShapely( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnColGroups( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnColResidues( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnColChain( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnColTemper( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnStruct( wxCommandEvent& event ) { event.Skip(); }

    // Display->Options menu
    virtual void OnSlab( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnHydrogen( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnHetero( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnSpecular( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnStereo( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnLabels( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnDisplayPBox( wxCommandEvent& event ) { event.Skip(); }

    // Display (other)
    virtual void DoViewParamDlg( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoObject3DDialog( wxCommandEvent& event ) { event.Skip(); }

    // Display->Labels menu
    virtual void OnLabelsAtomId( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnLabelsAtomNames( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnLabelsAtomSeqNum( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnLabelsAtomFFSymb( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnLabelsAtomCharges( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnLabelsOff( wxCommandEvent& event ) { event.Skip(); }

    // Display->Windows menu
    virtual void OnWindowsCascade( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnWindowsTileVert( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnCenterViewSel( wxCommandEvent& event ) { event.Skip(); }

    // ET menu
    virtual void OnEditRedox( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoPathwaysDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoETEffHamDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnResetETModule( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoScatterDialog( wxCommandEvent& event ) { event.Skip(); }

    // QChem menu
    virtual void DoQChemParamDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoLoadQCDatDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoWaveFunAnalDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnCalcPolarGcontr( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnSaveGrpOperMat( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnCalcPolarContrF( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnReadPolarContr( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnCalcBetaContr2idx( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnReadBetaContr2idx( wxCommandEvent& event ) { event.Skip(); }

    // Basic Models / Applications menu
    virtual void DoMolMechDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoInterMolDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoContElectrDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoPNPDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoSmlPNPDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoRFDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoFlexModDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoEDDialog( wxCommandEvent& event ) { event.Skip(); }
    virtual void DoProtRedoxDialog( wxCommandEvent& event ) { event.Skip(); }

    // Help menu
    virtual void OnAbout( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnHelp( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnLoadManual( wxCommandEvent& event ) { event.Skip(); }

    // Test menu
    virtual void OnTestOper1( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnTestOper2( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnTestQCMod1( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnCalcTMatr1( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnTestMap( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnTestMin1( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnTestGraph1( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnTestPython1( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnDumpMolInfo( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnDumpGaussBCommon( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnDumpOverlap( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnDumpOverlap2( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnOpenVectorFieldInHaMolView( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnOpenDielConstFromNindexInHaMolView( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnOpenDiffConstFromNindexInHaMolView( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnOpenQstFromNindexInHaMolView( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnTestAvgPopFunc( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnSaveAmoebaTopFile1( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnSaveAmoebaTopFile2( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnTestFFT1( wxCommandEvent& event ) { event.Skip(); }

    // Toolbar
    virtual void OnMolConnect( wxCommandEvent& event ) { event.Skip(); }
    virtual void OnWorldConnect( wxCommandEvent& event ) { event.Skip(); }
};

#endif // MAIN_MENU_BASE_H
