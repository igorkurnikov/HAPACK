# QRDB=None

def load_qr_db():
    # Load QRDB
    import harlempy.molset
    print("Loading QRDB...")
    QRDB = harlempy.molset.GetQRDB()
    # import sys
    # print sys.path
    from . import QR_AMBER
    from . import QR_OPLS88
    from . import QR_PARSE94
    QRDB.PrintFFsInfo()


def on_pnpgui(_):
    """
    Launch PNPS Frame as response to menu clicking
    """
    import wx
    import harlempy.molset
    import harlempy.hapygui

    pmset = harlempy.molset.GetCurMolSet()
    if pmset is None:
        dlg = wx.MessageDialog(harlempy.hapygui.GetMainFrame(), "MolSet is not loaded", "Message from PNPGUI", wx.OK)
        dlg.ShowModal()  # Show it
        dlg.Destroy()  # finally destroy it when finished.
    else:
        from . import pnpsgui
        pnpsgui.PNPFrame(pmset)


def ask_for_file(message, wildcard):
    import wx
    import harlempy.molset
    import harlempy.hapygui

    # Check if molset exists
    # @todo create new molset and view
    if harlempy.molset.GetCurMolSet() is None:
        wx.MessageBox(
            "There is no MolSet. Load molecule or create an empty molset.", "MolSet is Needed",
            style=wx.OK | wx.CENTRE | wx.ICON_EXCLAMATION)
        return None

    with wx.FileDialog(
            harlempy.hapygui.GetMainFrame(),
            message, wildcard=wildcard,
            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as fileDialog:

        if fileDialog.ShowModal() == wx.ID_CANCEL:
            return None

        return fileDialog.GetPath()


def view_hafield3d(field, title):
    """
    View HaField3D
    """
    import harlempy.molset
    import harlempy.hapygui
    # Initiate PlaneViewOfHaField3D
    plane_view = harlempy.molset.PlaneViewOfHaField3D(field, title)
    plane_view.SetHideZeroValues(True)
    plane_view.thisown = 0

    # Redraw
    pmset = harlempy.molset.GetCurMolSet()
    pmset.AddObject3D(plane_view)
    pView = pmset.GetActiveMolView()
    if pView is not None:
        pView.ReDrawFlag = pView.ReDrawFlag | harlempy.molset.RFInitial
        pView.InitialTransform()
    pmset.RefreshAllViews(harlempy.molset.RFRefresh)

    # Create wxFieldPlaneView to control plane
    harlempy.molset.CreatewxFieldPlaneView(plane_view, pmset, title, 1)


def on_open_vector_field(_):
    """
    View vector field potential from file
    """
    import harlempy.molset
    import harlempy.hapygui
    import pnps

    filename = ask_for_file(
        message="Open Vector Field File",
        wildcard="Vector Field3D (*.bin)|*.bin|Vector Field3D (*.gz)|*.gz")
    if filename is None:
        return

    # Load vector field
    vector_field = pnps.VectorField3D(filename, 1.0)

    for i in range(vector_field.GetNelem()):
        # Get zeroth field and initiate its boundaries
        field = harlempy.molset.HaField3D(
            vector_field.GetField(0), vector_field.GetNx(), vector_field.GetNy(), vector_field.GetNz(), True)
        field.SetCenterAsZero(vector_field.GridScale)
        # Set taken filed to None to avoid its deallocation
        vector_field.SetField(0, None)
        # Set thisown to avoid deallocation
        field.thisown = 0
        # View field
        view_hafield3d(field, "%s:%d" % (filename, i))

    del vector_field

def on_open_dielectric(_):
    """
    View dielectric map from NIndex file
    """
    import harlempy.molset
    import harlempy.hapygui
    import pnps

    filename = ask_for_file(
        message="Open Diel. Map from System Topology",
        wildcard="System Topology (*.systop)|*.systop")
    if filename is None:
        return

    # Load vector field
    nindex = pnps.NodeIndexing()
    nindex.ReadFromFile(filename)

    params_list = [
        ["X", pnps.NodeIndexing.DielConst, pnps.NodeIndexing.Epsilon0, 0.5 / nindex.GridScale, 0, 0],
        ["Y", pnps.NodeIndexing.DielConst, pnps.NodeIndexing.Epsilon1, 0, 0.5 / nindex.GridScale, 0],
        ["Z", pnps.NodeIndexing.DielConst, pnps.NodeIndexing.Epsilon2, 0, 0, 0.5 / nindex.GridScale]
    ]

    for name, field_type, mask, dx, dy, dz in params_list:
        title = "Epsilon" + name + "_" + filename
        # Get zeroth field and initiate its boundaries
        field = harlempy.molset.HaField3D(
            nindex.GetField(field_type, mask), nindex.GetNx(), nindex.GetNy(), nindex.GetNz(), True)
        field.SetCenterAsZero(nindex.GridScale)
        field.SetName(title)
        field.ShiftGridCorners(dx, dy, dz)
        # Set taken filed to None to avoid its deallocation
        # Set thisown to avoid deallocation
        field.thisown = 0
        # View field
        view_hafield3d(field, title)


def on_open_diffusion(_):
    """
    View diffusion from NIndex file
    """
    import harlempy.molset
    import harlempy.hapygui
    import pnps

    filename = ask_for_file(
        message="Open Diffusion Map from System Topology",
        wildcard="System Topology (*.systop)|*.systop")
    if filename is None:
        return

    # Load vector field
    nindex = pnps.NodeIndexing()
    nindex.ReadFromFile(filename)

    params_list = [
        ["0", pnps.NodeIndexing.DiffConst, pnps.NodeIndexing.Ion0],
        ["1", pnps.NodeIndexing.DiffConst, pnps.NodeIndexing.Ion1]
    ]

    for name, field_type, mask in params_list:
        title = "Epsilon" + name + "_" + filename
        # Get zeroth field and initiate its boundaries
        field = harlempy.molset.HaField3D(
            nindex.GetField(field_type, mask), nindex.GetNx(), nindex.GetNy(), nindex.GetNz(), True)
        field.SetCenterAsZero(nindex.GridScale)
        field.SetName(title)
        # Set taken filed to None to avoid its deallocation
        # Set thisown to avoid deallocation
        field.thisown = 0
        # View field
        view_hafield3d(field, title)


def init():
    """
    Add PNP GUI Menu to MainFrame
    """
    import wx
    import harlempy.molset
    import harlempy.hapygui

    # menu item for main GUI dialog
    menu_simple_gui = harlempy.hapygui.GetApplicationsMenu().Append(
        harlempy.hapygui.NewID(), "PNP GUI", "GUI for PNP/PNP-SR calculations")
    harlempy.hapygui.GetMainFrame().Bind(wx.EVT_MENU, on_pnpgui, menu_simple_gui)

    # utilities
    menu_pnps_utils = wx.Menu()
    open_potential = menu_pnps_utils.Append(
        harlempy.hapygui.NewID(), "Open Potential Field")
    open_concentration = menu_pnps_utils.Append(
        harlempy.hapygui.NewID(), "Open Concentration Field")
    open_dielectric = menu_pnps_utils.Append(
        harlempy.hapygui.NewID(), "Open Dielectric Map")
    open_diffusion = menu_pnps_utils.Append(
        harlempy.hapygui.NewID(), "Open Diffusion Map")
    open_vector_field = menu_pnps_utils.Append(
        harlempy.hapygui.NewID(), "Open Vector Field")

    harlempy.hapygui.GetApplicationsMenu().Append(
        harlempy.hapygui.NewID(), "PNP Utils", menu_pnps_utils)

    harlempy.hapygui.GetMainFrame().Bind(wx.EVT_MENU, on_open_vector_field, open_potential)
    harlempy.hapygui.GetMainFrame().Bind(wx.EVT_MENU, on_open_vector_field, open_concentration)
    harlempy.hapygui.GetMainFrame().Bind(wx.EVT_MENU, on_open_dielectric, open_dielectric)
    harlempy.hapygui.GetMainFrame().Bind(wx.EVT_MENU, on_open_diffusion, open_diffusion)
    harlempy.hapygui.GetMainFrame().Bind(wx.EVT_MENU, on_open_vector_field, open_vector_field)
