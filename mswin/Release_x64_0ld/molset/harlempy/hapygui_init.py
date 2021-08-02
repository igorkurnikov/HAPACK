import wx
import harlempy.hapygui as hapygui
import molset

def OnCenterAtOrigin(e):
    print("Moving molecules geometric center to the coordinates origin ...")
    pmset=molset.GetCurMolSet()
    if pmset is None:
        print ("There is no MolSet selected")
        return
    moleditor=molset.MolEditor()
    moleditor.CenterAtOrigin(pmset)
    pmset.RefreshAllViews( molset.RFRefresh | molset.RFColour | molset.RFApply )
    print("    Done")

menuOnCenterAtOrigin = hapygui.GetEditMenu().Append(
    hapygui.NewID(), "Center Molecules to Origin", "Move molecules geometric center to the coordinates origin")
hapygui.GetMainFrame().Bind(wx.EVT_MENU, OnCenterAtOrigin, menuOnCenterAtOrigin)

try:
    import pnpsgui
    pnpsgui.init()
except Exception as e:
     print(" ")
#    print("Warning: no pnpgui loaded:")
#    print("Can not load pnpgui")
#    print(str(e))
#    import traceback
#    traceback.print_exc()
