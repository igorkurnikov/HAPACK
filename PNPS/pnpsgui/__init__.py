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


def on_open_potential(_):
    import wx
    import harlempy.molset
    import harlempy.hapygui
    import pnps

    filename = ask_for_file(
        message="Open Potential Field File", wildcard="Vector Field3D (*.bin)|*.bin|Vector Field3D (*.gz)|*.gz")
    if filename is None:
        return

    # Load vector field
    vector_field = pnps.VectorField3D(filename, 1.0)

    # Get zeroth field and initiate its boundaries
    field = harlempy.molset.HaField3D(
        vector_field.GetField(0), vector_field.GetNx(), vector_field.GetNy(), vector_field.GetNz(), True)
    field.SetCenterAsZero(vector_field.GridScale)

    # Set taken filed to None to avoid its deallocation
    vector_field.SetField(0, None)
    del vector_field
    # Set thisown to avoid deallocation
    field.thisown = 0

    # Initiate PlaneViewOfHaField3D
    PlaneV = harlempy.molset.PlaneViewOfHaField3D(field, filename)
    PlaneV.SetHideZeroValues(True)
    PlaneV.thisown = 0

    # Redraw
    pmset = harlempy.molset.GetCurMolSet()
    pmset.AddObject3D(PlaneV)
    pView = pmset.GetActiveMolView()
    if pView is not None:
        pView.ReDrawFlag = pView.ReDrawFlag | harlempy.molset.RFInitial
        pView.InitialTransform()
    pmset.RefreshAllViews(harlempy.molset.RFRefresh)

    # Create wxFieldPlaneView to control plane
    harlempy.molset.CreatewxFieldPlaneView(PlaneV, pmset, filename, 1)


def on_open_concentration(_):
    pass


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
    #open_concentration = menu_pnps_utils.Append(
    #    harlempy.hapygui.NewID(), "Open Concentration Field")

    harlempy.hapygui.GetApplicationsMenu().Append(
        harlempy.hapygui.NewID(), "PNP Utils", menu_pnps_utils)

    harlempy.hapygui.GetMainFrame().Bind(wx.EVT_MENU, on_open_potential, open_potential)
    #harlempy.hapygui.GetMainFrame().Bind(wx.EVT_MENU, on_open_concentration, open_concentration)

