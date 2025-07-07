
# QRDB=None

def load_qr_db():
    #Load QRDB
    print "Loading QRDB..."
    QRDB=GetQRDB()
    #import sys
    #print sys.path
    import QR_AMBER
    import QR_OPLS88
    import QR_PARSE94
    QRDB.PrintFFsInfo()


def on_pnpgui(_):
    # A message dialog box with an OK button. wx.OK is a standard ID in wxWidgets.
    import molset
    pmset=molset.GetCurMolSet()
    if pmset is None:
        dlg = wx.MessageDialog(hapygui.GetMainFrame(), "MolSet is not loaded", "Message from PNPGUI", wx.OK)
        dlg.ShowModal() # Show it
        dlg.Destroy() # finally destroy it when finished.
    else:
        import pnpgui
        import molset
        # @todo add ppp
        # pnpgui.pApp = molset.pApp
        pnpgui.PNPFrame(pmset)


def init():
    """Add PNP GUI Menu to MainFrame"""
    import wx
    import hapygui

    menu_simple_gui = hapygui.GetApplicationsMenu().Append(hapygui.NewID(), "PNP GUI","GUI for PNP/PNP-SR calculations")
    hapygui.GetMainFrame().Bind(wx.EVT_MENU, on_pnpgui, menu_simple_gui)
