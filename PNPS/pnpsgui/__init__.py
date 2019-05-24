# QRDB=None

def load_qr_db():
    # Load QRDB
    print("Loading QRDB...")
    QRDB = GetQRDB()
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
        # import molset
        # @todo add ppp
        # pnpsgui.pApp = molset.pApp
        pnpsgui.PNPFrame(pmset)


def init():
    """
    Add PNP GUI Menu to MainFrame
    """
    import wx
    import harlempy.hapygui

    menu_simple_gui = harlempy.hapygui.GetApplicationsMenu().Append(harlempy.hapygui.NewID(), "PNP GUI",
                                                                    "GUI for PNP/PNP-SR calculations")
    harlempy.hapygui.GetMainFrame().Bind(wx.EVT_MENU, on_pnpgui, menu_simple_gui)
