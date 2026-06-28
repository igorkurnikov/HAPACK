/*! \file main_menu_base.cpp

    Base class for HARLEM main frame menu bar - wxFormBuilder style.
    Builds the complete menu and binds events to virtual handlers.

    \author Igor Kurnikov
    \date 2026-
*/

#include "hampi.h"

#include "wx/wx.h"

#include "ha_wx_res_wdr.h"
#include "main_menu_base.h"

HaMainMenuBase::HaMainMenuBase( wxWindow* parent, wxWindowID id, const wxString& title,
    const wxPoint& pos, const wxSize& size, long style, const wxString& name )
    : wxMDIParentFrame( parent, id, title, pos, size, style, name )
{
    // ---------------------------------------------------------------
    // Build menu bar (replaces MainMenu() from ha_wx_res_wdr.cpp)
    // ---------------------------------------------------------------

    wxMenuBar* menuBar = new wxMenuBar;

    // === File Menu ===
    wxMenu* fileMenu = new wxMenu;
    fileMenu->Append( ID_FILE_NEW_WX, wxT("New"), wxT("") );
    fileMenu->Append( ID_FILE_OPEN_WX, wxT("Open Molecule\tCtrl-o"), wxT("") );
    fileMenu->Append( IDM_MOL_SAVE_WX, wxT("Save Molecule\tCtrl-s"), wxT("") );
    fileMenu->Append( IDM_MOL_SAVE_AS_WX, wxT("Save Molecule As"), wxT("") );
    fileMenu->AppendSeparator();
    fileMenu->Append( IDM_LOAD_SCRIPT_WX, wxT("Load Script\tAlt-Ctrl-s"), wxT("") );
    fileMenu->Append( IDM_REDIRECT_IO_WX, wxT("Redirect IO"), wxT("") );
    fileMenu->AppendSeparator();
    fileMenu->Append( IDM_SAVE_IMAGE_WX, wxT("Save Image to File"), wxT("") );
    fileMenu->Append( IDM_SAVE_IMAGE_CLIPBOARD_WX, wxT("Save Image to Clipboard"), wxT("") );
    fileMenu->AppendSeparator();
    fileMenu->Append( IDM_COMP_ACCOUNTS_WX, wxT("Computer Accounts"), wxT("") );
    fileMenu->Append( IDM_INFO_WX, wxT("Information"), wxT("") );
    fileMenu->Append( ID_FILE_CLOSE_WX, wxT("Close View"), wxT("") );
    fileMenu->AppendSeparator();
    fileMenu->Append( IDM_PRINT_WX, wxT("Print"), wxT("") );
    fileMenu->Append( IDM_SETUP_WX, wxT("Print Setup..."), wxT("") );
    fileMenu->AppendSeparator();
    fileMenu->Append( IDM_EXIT_WX, wxT("Exit"), wxT("") );
    menuBar->Append( fileMenu, wxT("File") );

    // === Edit Menu ===
    wxMenu* editMenu = new wxMenu;
    editMenu->Append( IDM_SELECT_WX, wxT("Select All"), wxT("") );
    editMenu->AppendSeparator();
    editMenu->Append( IDM_EDIT_ATOM_PARAM_WX, wxT("Edit Atom Parameters"), wxT("") );
    editMenu->Append( IDM_EDIT_RES_PARAM_WX, wxT("Edit Residue Parameters"), wxT("") );
    editMenu->Append( IDM_EDIT_MOLSETS, wxT("Edit Molecular Sets"), wxT("") );
    editMenu->Append( IDM_EDIT_GROUPS_WX, wxT("Edit Groups"), wxT("") );
    editMenu->Append( IDM_EDIT_MUT_MAP, wxT("Edit Mutation Map"), wxT("") );
    editMenu->Append( IDM_CRD_SNAPSHOT, wxT("Edit Coordinate Snapshots"), wxT("") );
    editMenu->Append( IDM_EDIT_GEOM_WX, wxT("Edit Geometry"), wxT("") );
    editMenu->Append( IDM_FIND_HBONDS_WX, wxT("Find Hydrogen Bonds"), wxT("") );
    editMenu->Append( IDM_SOLVATE_WX, wxT("Solvate Molecules"), wxT("") );
    editMenu->Append( IDM_WRAP_UNIT_CELL, wxT("Wrap Into Unit Cell"), wxT("") );
    editMenu->Append( IDM_CENTER_SOLUTE, wxT("Center Solute in Solvent"), wxT("") );
    editMenu->Append( IDM_CENTER_MOL_PBOX, wxT("Center Molecules in Periodical Box"), wxT("") );

    // Add and Delete Atoms submenu
    wxMenu* addDelAtomsMenu = new wxMenu;
    addDelAtomsMenu->Append( IDM_ADD_MISSING_ATOMS_WX, wxT("Add All Missing Atoms"), wxT("") );
    addDelAtomsMenu->Append( IDM_ADD_POLAR_HYDROGENS_WX, wxT("Add Polar Hydrogens"), wxT("") );
    addDelAtomsMenu->Append( IDM_ADD_H_HYBRID_WX, wxT("Add H Using Hybridization"), wxT("") );
    addDelAtomsMenu->Append( IDM_ADD_HYDROGENS_WX, wxT("Add Hydrogens"), wxT("") );
    addDelAtomsMenu->Append( IDM_DEL_SEL_ATOMS_WX, wxT("Delete Selected Atoms"), wxT("") );
    addDelAtomsMenu->Append( IDM_DEL_OVLP_MOLS, wxT("Delete Molecules Overlapping Selection"), wxT("") );
    editMenu->Append( ID_MENU, wxT("Add and Delete Atoms"), addDelAtomsMenu );

    editMenu->Append( IDM_EDIT_FRAGM_WX, wxT("Fragments"), wxT("") );

    // Construct Molecule submenu
    wxMenu* constructMolMenu = new wxMenu;
    constructMolMenu->Append( IDM_BUILD_FILM_WX, wxT("Build Molecular Film"), wxT("") );
    editMenu->Append( ID_MENU, wxT("Construct Molecule"), constructMolMenu );

    // Pick and Select submenu
    wxMenu* pickSelMenu = new wxMenu;
    pickSelMenu->Append( IDM_CLEAR_PICKED_WX, wxT("Clear Picked Atoms"), wxT("") );
    pickSelMenu->Append( IDM_SEL_ATOMS_IN_BOUND_BOX, wxT("Select Atoms Within Boundary Box"), wxT("") );
    pickSelMenu->Append( IDM_REVERT_SELECTION, wxT("Revert Atom Selection"), wxT("") );
    pickSelMenu->Append( IDM_EXPAND_SELECTION_BONDED, wxT("Expand Atom Selection with Bonded Atoms"), wxT("") );
    editMenu->Append( ID_MENU, wxT("Pick and Select"), pickSelMenu );

    // Biopolymers submenu
    wxMenu* biopolymersMenu = new wxMenu;
    biopolymersMenu->Append( IDM_NUCL_ACID_WX, wxT("Nucleic Acids Manipulation"), wxT("") );
    biopolymersMenu->Append( IDM_DESCRIBE_SEC_STRUCT, wxT("Describe Secondary Structure"), wxT("") );
    biopolymersMenu->Append( IDM_PRINT_HBONDS, wxT("Print Hydrogen Bonds"), wxT("") );
    biopolymersMenu->Append( IDM_SET_ALPHA_HELIX, wxT("Set Selection To Alpha Helix"), wxT("") );
    editMenu->Append( ID_MENU, wxT("Biopolymers"), biopolymersMenu );

    editMenu->Append( IDM_SHOW_RES_DB_WX, wxT("Show Residue DB"), wxT("") );

    // Measure submenu
    wxMenu* measureMenu = new wxMenu;
    measureMenu->Append( IDM_SHOW_ATOMID_WX, wxT("Atom ID"), wxT("") );
    measureMenu->Append( IDM_MEASURE_DISTANCE_WX, wxT("Distance"), wxT("") );
    measureMenu->Append( IDM_MEASURE_ANGLE_WX, wxT("Angle"), wxT("") );
    measureMenu->Append( IDM_MEASURE_DIHEDRAL_WX, wxT("Dihedral"), wxT("") );
    measureMenu->Append( IDM_ADD_DIST_MONITOR, wxT("Add Distance Monitor"), wxT("") );
    measureMenu->Append( IDM_SHOW_ATOM_LABEL, wxT("Show Atom Label"), wxT("") );
    editMenu->Append( ID_MENU, wxT("Measure"), measureMenu );

    menuBar->Append( editMenu, wxT("Edit") );

    // === Display Menu ===
    wxMenu* displayMenu = new wxMenu;

    // Display Mode submenu
    wxMenu* displayModeMenu = new wxMenu;
    displayModeMenu->Append( IDM_WIREFRAME_WX, wxT("Wireframe"), wxT("") );
    displayModeMenu->Append( IDM_STICKS_WX, wxT("Sticks"), wxT("") );
    displayModeMenu->Append( IDM_SPHERES_WX, wxT("Spacefill"), wxT("") );
    displayModeMenu->Append( IDM_BALLSTICK_WX, wxT("Ball & Stick"), wxT("") );
    displayModeMenu->Append( IDM_RIBBONS_WX, wxT("Ribbons"), wxT("") );
    displayModeMenu->Append( IDM_STRANDS_WX, wxT("Strands"), wxT("") );
    displayModeMenu->Append( IDM_CARTOONS_WX, wxT("Cartoons"), wxT("") );
    displayMenu->Append( ID_MENU, wxT("Display Mode"), displayModeMenu );

    // Colours submenu
    wxMenu* coloursMenu = new wxMenu;
    coloursMenu->Append( IDM_MONO_WX, wxT("Monochrome"), wxT("") );
    coloursMenu->Append( IDM_CPK_WX, wxT("CPK"), wxT("") );
    coloursMenu->Append( IDM_SHAPELY_WX, wxT("Shapely"), wxT("") );
    coloursMenu->Append( IDM_COL_GROUPS_WX, wxT("Groups"), wxT("") );
    coloursMenu->Append( IDM_COL_RESIDUES_WX, wxT("Residues"), wxT("") );
    coloursMenu->Append( IDM_CHAIN_WX, wxT("Chains"), wxT("") );
    coloursMenu->Append( IDM_TEMPER_WX, wxT("Temperature"), wxT("") );
    coloursMenu->Append( IDM_STRUCT_WX, wxT("Structure"), wxT("") );
    coloursMenu->Append( IDM_COL_DON_COUPL_WX, wxT("Coupling to Donor"), wxT("") );
    coloursMenu->Append( IDM_ATOM_PROP_COLORS, wxT("Atom Properties"), wxT("") );
    displayMenu->Append( ID_MENU, wxT("Colours"), coloursMenu );

    // Options submenu (checkable items)
    wxMenu* optionsMenu = new wxMenu;
    optionsMenu->Append( IDM_SLAB_WX, wxT("Slab Mode"), wxT(""), wxITEM_CHECK );
    optionsMenu->Append( IDM_HYDROGEN_WX, wxT("Hydrogens"), wxT(""), wxITEM_CHECK );
    optionsMenu->Append( IDM_HETERO_WX, wxT("Hetero Atoms"), wxT(""), wxITEM_CHECK );
    optionsMenu->Append( IDM_SPECULAR_WX, wxT("Specular"), wxT(""), wxITEM_CHECK );
    optionsMenu->Append( IDM_STEREO_WX, wxT("Stereo"), wxT(""), wxITEM_CHECK );
    optionsMenu->Append( IDM_LABELS_WX, wxT("Labels"), wxT(""), wxITEM_CHECK );
    optionsMenu->Append( IDM_DISPLAY_PBOX, wxT("Periodical Boundary Box"), wxT(""), wxITEM_CHECK );
    displayMenu->Append( ID_MENU, wxT("Options"), optionsMenu );

    displayMenu->Append( IDM_VIEW_PARAM_WX, wxT("Set View Parameters"), wxT("") );
    displayMenu->Append( IDM_OBJECT3D_WX, wxT("3D Objects"), wxT("") );

    // Labels submenu
    wxMenu* labelsMenu = new wxMenu;
    labelsMenu->Append( IDM_LABELS_ATOM_ID_WX, wxT("Atom IDs"), wxT("") );
    labelsMenu->Append( IDM_LABELS_ATOM_NAMES_WX, wxT("Atom Names"), wxT("") );
    labelsMenu->Append( IDM_LABELS_ATOM_SEQNUM_WX, wxT("Atom Seq Numbers"), wxT("") );
    labelsMenu->Append( IDM_LABELS_ATOM_FF_SYMB, wxT("Force Field Symbols"), wxT("") );
    labelsMenu->Append( IDM_LABELS_ATOM_CHARGES, wxT("Atomic Charges"), wxT("") );
    labelsMenu->Append( IDM_LABELS_OFF_WX, wxT("OFF"), wxT("") );
    displayMenu->Append( ID_MENU, wxT("Labels"), labelsMenu );

    // Windows submenu
    wxMenu* windowsMenu = new wxMenu;
    windowsMenu->Append( ID_WINDOW_CASCADE_WX, wxT("Cascade"), wxT("") );
    windowsMenu->Append( ID_WINDOW_TILE_VERT_WX, wxT("Tile"), wxT("") );
    displayMenu->Append( ID_MENU, wxT("Windows"), windowsMenu );

    displayMenu->Append( IDM_CENTER_VIEW_SEL_WX, wxT("Center View Sel Atoms"), wxT("") );
    menuBar->Append( displayMenu, wxT("Display") );

    // === Basic Models Menu ===
    wxMenu* basicModelsMenu = new wxMenu;
    basicModelsMenu->Append( IDM_MM_PARAM_WX, wxT("Molecular Mechanics"), wxT("") );
    basicModelsMenu->Append( IDM_CONT_ELECTR_WX, wxT("Continuum Electrostatics"), wxT("") );

    // Quantum Chemistry submenu
    wxMenu* qchemMenu = new wxMenu;
    qchemMenu->Append( IDM_QCHEM_PARAM_WX, wxT("Parameters"), wxT("") );
    qchemMenu->Append( IDM_QCHEM_LOAD_DATA_WX, wxT("Load Data"), wxT("") );
    qchemMenu->Append( IDM_WAVEFUN_ANAL_WX, wxT("WaveFun Anal"), wxT("") );

    // Optical Properties submenu
    wxMenu* opticalPropMenu = new wxMenu;
    opticalPropMenu->Append( IDM_CALC_POLAR_GCONTR_WX, wxT("Polar Grp Contr"), wxT("") );
    opticalPropMenu->Append( IDM_SAVE_GROUP_OPER_WX, wxT("Save Group Oper Mat"), wxT("") );
    opticalPropMenu->Append( IDM_POL_GRP_CONTR_F_WX, wxT("Calc Polarizability Group Contrib"), wxT("") );
    opticalPropMenu->Append( IDM_READ_POLAR_CONTR_WX, wxT("Read Polarizabilty Group Contr"), wxT("") );
    opticalPropMenu->Append( IDM_CALC_ROTANG_GRP_CONTR2_WX, wxT("Calc Rot Angle Group Contrib"), wxT("") );
    opticalPropMenu->Append( IDM_READ_ROTANG_GRP_CONTR2_WX, wxT("Read Rot Angle Group Contr"), wxT("") );
    qchemMenu->Append( ID_MENU, wxT("Optical Properties"), opticalPropMenu );

    basicModelsMenu->Append( ID_MENU, wxT("Quantum Chemistry"), qchemMenu );
    menuBar->Append( basicModelsMenu, wxT("Basic Models") );

    // === Applications Menu ===
    wxMenu* applicationsMenu = new wxMenu;
    applicationsMenu->Append( IDM_INTER_MOL_WX, wxT("Inter Molecular Interactions"), wxT("") );
    applicationsMenu->Append( IDM_PNP_SMLPNP, wxT("Poisson Nernst Planck"), wxT("") );
    applicationsMenu->Append( IDM_PNP_RF, wxT("Potential of Mean Force"), wxT("") );

    // ET submenu
    wxMenu* etMenu = new wxMenu;
    etMenu->Append( IDM_EDIT_REDOX_WX, wxT("Edit Redox Centers"), wxT("") );
    etMenu->Append( IDM_RUN_PATHWAYS_WX, wxT("Pathways Model"), wxT("") );
    etMenu->Append( IDM_ET_EFF_HAM_WX, wxT("Donor-Acceptor QChem Calculations"), wxT("") );
    etMenu->Append( IDM_RESET_ET_WX, wxT("Reset ET Module"), wxT("") );
    etMenu->Append( IDM_EL_SCATTER_WX, wxT("Electron Scattering Calc"), wxT("") );
    applicationsMenu->Append( ID_MENU, wxT("ET"), etMenu );

    applicationsMenu->Append( ID_FLEX_MOD, wxT("Flexibility Analysis (FIRST)"), wxT("") );
    applicationsMenu->Append( ID_MENU_ED, wxT("Essential Dynamics and Clustering Analysis"), wxT("") );
    applicationsMenu->Append( IDM_PROT_REDOX, wxT("Protonation and Redox Equilibrium"), wxT("") );
    menuBar->Append( applicationsMenu, wxT("Applications") );

    // === Help Menu ===
    wxMenu* helpMenu = new wxMenu;
    helpMenu->Append( IDM_ABOUT_WX, wxT("About HARLEM"), wxT("") );
    helpMenu->Append( IDM_HELP_WX, wxT("Beginner's Manual"), wxT("") );
    helpMenu->Append( IDM_MANUAL_WX, wxT("Advanced Manual"), wxT("") );
    menuBar->Append( helpMenu, wxT("Help") );

    // === Test Menu ===
    wxMenu* testMenu = new wxMenu;
    testMenu->Append( IDM_TEST_OPER_1_WX, wxT("Test Operators 1"), wxT("") );
    testMenu->Append( IDM_TEST_OPER_2_WX, wxT("Test Operators 2"), wxT("") );
    testMenu->Append( IDM_TEST_QCMOD_1_WX, wxT("Test HaQCMod 1"), wxT("") );
    testMenu->Append( IDM_CALC_MATR_1_WX, wxT("Calc Matricies 1"), wxT("") );
    testMenu->Append( IDM_TEST_MAP_WX, wxT("Test STD Lib Map"), wxT("") );
    testMenu->Append( IDM_TEST_MIN_1_WX, wxT("Test Min 1"), wxT("") );
    testMenu->Append( IDM_TEST_GRAPH_1_WX, wxT("test_graph_1"), wxT("") );
    testMenu->Append( IDM_PNP, wxT("PNP Solver\tAlt-Ctrl-p"), wxT("") );
    testMenu->Append( IDM_CONT_ELECTR_WX, wxT("Continuum Electrostatics"), wxT("") );

    // Transform submenu
    wxMenu* transformMenu = new wxMenu;
    transformMenu->Append( IDC_TRANSFORM_CONNECT_WX, wxT("Connect"), wxT("") );
    transformMenu->AppendSeparator();
    transformMenu->Append( IDC_POSITION_WX, wxT("Position"), wxT("") );
    transformMenu->Append( IDC_POSITION_RESET_WX, wxT("Reset"), wxT("") );
    transformMenu->Append( IDC_POSITION_STORE_WX, wxT("Store"), wxT("") );
    transformMenu->Append( IDC_POSITION_AXES_WX, wxT("Axis"), wxT("") );
    testMenu->Append( ID_MENU, wxT("Transform"), transformMenu );

    testMenu->Append( IDM_TEST_PYTHON_1, wxT("Test Python 1"), wxT("") );

    // Dump submenu
    wxMenu* dumpMenu = new wxMenu;
    dumpMenu->Append( IDM_MOL_INFO_WX, wxT("Molecule Info"), wxT("") );
    dumpMenu->Append( IDM_GAUSS_BCOMMON_WX, wxT("Gaussian basis common"), wxT("") );
    dumpMenu->Append( IDM_DUMP_OVERLAP, wxT("Overlap Matrix"), wxT("") );
    dumpMenu->Append( IDM_DUMP_OVERLAP2_WX, wxT("Overlap between bas sets"), wxT("") );
    testMenu->Append( ID_MENU, wxT("Dump"), dumpMenu );

    // PNPS submenu
    wxMenu* pnpsMenu = new wxMenu;
    pnpsMenu->Append( IDM_OPENVECFIELD_HAVIEW, wxT("Open VectorField3D in HMolView"), wxT("") );
    pnpsMenu->Append( IDM_OPEN_QST_NIND, wxT("Open Qstatic Map From NodeIndex File"), wxT("") );
    pnpsMenu->Append( IDM_OPEN_EPS_NIND, wxT("Open Diel.Conts. Maps From NodeIndex File"), wxT("") );
    pnpsMenu->Append( IDM_OPEN_DIFF_NIND, wxT("Open Diff.Conts. Maps From NodeIndex File"), wxT("") );
    testMenu->Append( ID_MENU, wxT("PNPS"), pnpsMenu );

    testMenu->Append( IDM_TEST_MOLSET_HANDLER, wxT("Test Molset Handler"), wxT("") );

    // MORT submenu
    wxMenu* mortMenu = new wxMenu;
    mortMenu->Append( IDM_SAVE_AMOEBA_TOP_FILE_1, wxT("Save Amoeba Top File 1"), wxT("") );
    mortMenu->Append( IDM_SAVE_AMOEBA_TOP_FILE_2, wxT("Save Amoeba Top File 2"), wxT("") );
    testMenu->Append( ID_MENU, wxT("MORT"), mortMenu );

    // MATH submenu
    wxMenu* mathMenu = new wxMenu;
    mathMenu->Append( IDM_TEST_FFT_1, wxT("Test FFT 1"), wxT("") );
    testMenu->Append( ID_MENU, wxT("MATH"), mathMenu );

    // PKA REDOX submenu
    wxMenu* pkaRedoxMenu = new wxMenu;
    pkaRedoxMenu->Append( IDM_TEST_AVG_POP_FUNC, wxT("Test Avg Population functions"), wxT("") );
    testMenu->Append( ID_MENU, wxT("PKA REDOX"), pkaRedoxMenu );

    menuBar->Append( testMenu, wxT("Test") );

    SetMenuBar( menuBar );

    // ---------------------------------------------------------------
    // Bind menu events to virtual handlers
    // ---------------------------------------------------------------

    // File menu
    Bind( wxEVT_MENU, &HaMainMenuBase::OnFileNew, this, ID_FILE_NEW_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnFileOpen, this, ID_FILE_OPEN_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnMolSave, this, IDM_MOL_SAVE_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnMolSaveAs, this, IDM_MOL_SAVE_AS_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoLoadScriptDialog, this, IDM_LOAD_SCRIPT_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoRedirectIODialog, this, IDM_REDIRECT_IO_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::SaveOutputFile, this, IDM_SAVE_IMAGE_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnSaveImageClipboard, this, IDM_SAVE_IMAGE_CLIPBOARD_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoCompAccountsDialog, this, IDM_COMP_ACCOUNTS_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnInfo, this, IDM_INFO_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnFileClose, this, ID_FILE_CLOSE_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnPrint, this, IDM_PRINT_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnSetup, this, IDM_SETUP_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnExit, this, IDM_EXIT_WX );

    // Edit menu
    Bind( wxEVT_MENU, &HaMainMenuBase::OnSelectAll, this, IDM_SELECT_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoAtomParamsDialog, this, IDM_EDIT_ATOM_PARAM_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoResidueParamsDialog, this, IDM_EDIT_RES_PARAM_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoMolSetParamDialog, this, IDM_EDIT_MOLSETS );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoAtomPropColorDialog, this, IDM_ATOM_PROP_COLORS );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoEditGroupsDialog, this, IDM_EDIT_GROUPS_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoCrdSnapshotDialog, this, IDM_CRD_SNAPSHOT );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoEditGeomDialog, this, IDM_EDIT_GEOM_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoEditMutationMapDialog, this, IDM_EDIT_MUT_MAP );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnFindHbonds, this, IDM_FIND_HBONDS_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoSolvateDialog, this, IDM_SOLVATE_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnWrapIntoUnitCell, this, IDM_WRAP_UNIT_CELL );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnCenterSolute, this, IDM_CENTER_SOLUTE );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnCenterMolInPBox, this, IDM_CENTER_MOL_PBOX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnAddMissingAtoms, this, IDM_ADD_MISSING_ATOMS_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnAddPolarHydrogens, this, IDM_ADD_POLAR_HYDROGENS_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnAddHHybrid, this, IDM_ADD_H_HYBRID_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnAddHydrogens, this, IDM_ADD_HYDROGENS_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnDelSelAtoms, this, IDM_DEL_SEL_ATOMS_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnDelOvlpMols, this, IDM_DEL_OVLP_MOLS );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoEditFragmDialog, this, IDM_EDIT_FRAGM_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoBuildFilmDialog, this, IDM_BUILD_FILM_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnClearPicked, this, IDM_CLEAR_PICKED_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnSelAtomsInBoundBox, this, IDM_SEL_ATOMS_IN_BOUND_BOX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnRevertAtomSelection, this, IDM_REVERT_SELECTION );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnExpandAtomSelectionBonded, this, IDM_EXPAND_SELECTION_BONDED );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoNuclAcidDialog, this, IDM_NUCL_ACID_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnDescribeSecStruct, this, IDM_DESCRIBE_SEC_STRUCT );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnPrintHBonds, this, IDM_PRINT_HBONDS );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnSetAlphaHelix, this, IDM_SET_ALPHA_HELIX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnShowResDb, this, IDM_SHOW_RES_DB_WX );

    // Edit->Measure submenu
    Bind( wxEVT_MENU, &HaMainMenuBase::OnShowAtomID, this, IDM_SHOW_ATOMID_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnMeasureDist, this, IDM_MEASURE_DISTANCE_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnMeasureAngle, this, IDM_MEASURE_ANGLE_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnMeasureDihed, this, IDM_MEASURE_DIHEDRAL_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnAddDistMonitor, this, IDM_ADD_DIST_MONITOR );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnShowAtomLabel, this, IDM_SHOW_ATOM_LABEL );

    // Display->Display Mode menu
    Bind( wxEVT_MENU, &HaMainMenuBase::OnWireFrame, this, IDM_WIREFRAME_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnSticks, this, IDM_STICKS_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnSpheres, this, IDM_SPHERES_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnBallStick, this, IDM_BALLSTICK_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnRibbons, this, IDM_RIBBONS_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnStrands, this, IDM_STRANDS_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnCartoons, this, IDM_CARTOONS_WX );

    // Display->Colours menu
    Bind( wxEVT_MENU, &HaMainMenuBase::OnMono, this, IDM_MONO_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnCPK, this, IDM_CPK_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnShapely, this, IDM_SHAPELY_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnColGroups, this, IDM_COL_GROUPS_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnColResidues, this, IDM_COL_RESIDUES_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnColChain, this, IDM_CHAIN_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnColTemper, this, IDM_TEMPER_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnStruct, this, IDM_STRUCT_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnColTemper, this, IDM_COL_DON_COUPL_WX );

    // Display->Options menu
    Bind( wxEVT_MENU, &HaMainMenuBase::OnSlab, this, IDM_SLAB_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnHydrogen, this, IDM_HYDROGEN_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnHetero, this, IDM_HETERO_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnSpecular, this, IDM_SPECULAR_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnStereo, this, IDM_STEREO_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnLabels, this, IDM_LABELS_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnDisplayPBox, this, IDM_DISPLAY_PBOX );

    // Display (other)
    Bind( wxEVT_MENU, &HaMainMenuBase::DoViewParamDlg, this, IDM_VIEW_PARAM_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoObject3DDialog, this, IDM_OBJECT3D_WX );

    // Display->Labels menu
    Bind( wxEVT_MENU, &HaMainMenuBase::OnLabelsAtomId, this, IDM_LABELS_ATOM_ID_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnLabelsAtomNames, this, IDM_LABELS_ATOM_NAMES_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnLabelsAtomSeqNum, this, IDM_LABELS_ATOM_SEQNUM_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnLabelsAtomFFSymb, this, IDM_LABELS_ATOM_FF_SYMB );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnLabelsAtomCharges, this, IDM_LABELS_ATOM_CHARGES );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnLabelsOff, this, IDM_LABELS_OFF_WX );

    // Display->Windows menu
    Bind( wxEVT_MENU, &HaMainMenuBase::OnWindowsCascade, this, ID_WINDOW_CASCADE_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnWindowsTileVert, this, ID_WINDOW_TILE_VERT_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnCenterViewSel, this, IDM_CENTER_VIEW_SEL_WX );

    // ET menu
    Bind( wxEVT_MENU, &HaMainMenuBase::OnEditRedox, this, IDM_EDIT_REDOX_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoPathwaysDialog, this, IDM_RUN_PATHWAYS_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoETEffHamDialog, this, IDM_ET_EFF_HAM_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnResetETModule, this, IDM_RESET_ET_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoScatterDialog, this, IDM_EL_SCATTER_WX );

    // QChem menu
    Bind( wxEVT_MENU, &HaMainMenuBase::DoQChemParamDialog, this, IDM_QCHEM_PARAM_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoLoadQCDatDialog, this, IDM_QCHEM_LOAD_DATA_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoWaveFunAnalDialog, this, IDM_WAVEFUN_ANAL_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnCalcPolarGcontr, this, IDM_CALC_POLAR_GCONTR_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnSaveGrpOperMat, this, IDM_SAVE_GROUP_OPER_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnCalcPolarContrF, this, IDM_POL_GRP_CONTR_F_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnReadPolarContr, this, IDM_READ_POLAR_CONTR_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnCalcBetaContr2idx, this, IDM_CALC_ROTANG_GRP_CONTR2_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnReadBetaContr2idx, this, IDM_READ_ROTANG_GRP_CONTR2_WX );

    // Basic Models / Applications menu
    Bind( wxEVT_MENU, &HaMainMenuBase::DoMolMechDialog, this, IDM_MM_PARAM_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoInterMolDialog, this, IDM_INTER_MOL_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoContElectrDialog, this, IDM_CONT_ELECTR_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoPNPDialog, this, IDM_PNP );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoSmlPNPDialog, this, IDM_PNP_SMLPNP );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoRFDialog, this, IDM_PNP_RF );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoFlexModDialog, this, ID_FLEX_MOD );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoEDDialog, this, ID_MENU_ED );
    Bind( wxEVT_MENU, &HaMainMenuBase::DoProtRedoxDialog, this, IDM_PROT_REDOX );

    // Help menu
    Bind( wxEVT_MENU, &HaMainMenuBase::OnAbout, this, IDM_ABOUT_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnHelp, this, IDM_HELP_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnLoadManual, this, IDM_MANUAL_WX );

    // Test menu
    Bind( wxEVT_MENU, &HaMainMenuBase::OnTestOper1, this, IDM_TEST_OPER_1_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnTestOper2, this, IDM_TEST_OPER_2_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnTestQCMod1, this, IDM_TEST_QCMOD_1_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnCalcTMatr1, this, IDM_CALC_MATR_1_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnTestMap, this, IDM_TEST_MAP_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnTestMin1, this, IDM_TEST_MIN_1_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnTestGraph1, this, IDM_TEST_GRAPH_1_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnTestPython1, this, IDM_TEST_PYTHON_1 );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnDumpMolInfo, this, IDM_MOL_INFO_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnDumpGaussBCommon, this, IDM_GAUSS_BCOMMON_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnDumpOverlap, this, IDM_DUMP_OVERLAP );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnDumpOverlap2, this, IDM_DUMP_OVERLAP2_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnOpenVectorFieldInHaMolView, this, IDM_OPENVECFIELD_HAVIEW );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnOpenDielConstFromNindexInHaMolView, this, IDM_OPEN_EPS_NIND );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnOpenDiffConstFromNindexInHaMolView, this, IDM_OPEN_DIFF_NIND );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnOpenQstFromNindexInHaMolView, this, IDM_OPEN_QST_NIND );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnTestAvgPopFunc, this, IDM_TEST_AVG_POP_FUNC );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnSaveAmoebaTopFile1, this, IDM_SAVE_AMOEBA_TOP_FILE_1 );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnSaveAmoebaTopFile2, this, IDM_SAVE_AMOEBA_TOP_FILE_2 );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnTestFFT1, this, IDM_TEST_FFT_1 );

    // Toolbar
    Bind( wxEVT_MENU, &HaMainMenuBase::OnMolConnect, this, IDM_MOL_CONNECT_WX );
    Bind( wxEVT_MENU, &HaMainMenuBase::OnWorldConnect, this, IDM_WORLD_CONNECT_WX );
}

HaMainMenuBase::~HaMainMenuBase()
{
}
