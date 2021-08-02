import os

import wx

try:
    from harlempy import molset
except ImportError:
    molset = None

def FileBB_CheckFileExistence(FileBB, NoneIsValid=False):
    filename = str(FileBB.GetValue())
    if NoneIsValid and filename == "None":
        return True
    if os.path.exists(filename):
        return True
    return False


def FileBB_GetValue(FileBB):
    if FileBB is None:
        return None
    val = str(FileBB.GetValue())
    if val == "None":
        return None
    return val


def MolSet_RefreshAllViews(pmset):
    """Reset all views in `pmset` MolSet"""
    pview = pmset.GetActiveMolView()
    pview.ReDrawFlag = pview.ReDrawFlag | molset.RFInitial
    pview.InitialTransform()
    pmset.RefreshAllViews(molset.RFRefresh | molset.RFMagnify)


def ShowVectorField3D(filename, pmset, title):
    from .pnpsgui import main_frame, pnps, molset

    if not os.path.exists(filename):
        dlg = wx.MessageDialog(main_frame, "File %s does not exist" % filename, "Error In Input Values", wx.OK)
        dlg.ShowModal()  # Show it
        dlg.Destroy()  # finally destroy it when finished.
        return False

    VField = pnps.VectorField3D(filename)
    for i in range(VField.Nelem):
        field = VField.GetHaField3D(i)
        Name = "%s %d (%s , %s)" % (title, i, filename, pmset.GetName())
        field.SetName(Name)
        field.thisown = 0
        PlaneV = molset.PlaneViewOfHaField3D(field, Name, 1)
        PlaneV.SetHideZeroValues(1)
        PlaneV.thisown = 0
        pmset.AddObject3D(PlaneV)
        molset.CreatewxFieldPlaneView(PlaneV, pmset, Name, 1)
    pmset.RefreshAllViews(
        molset.RFRefresh | molset.RFColour | molset.RFApply | molset.RFMagnify | molset.RFRotate | molset.RFTransX | molset.RFTransY)
    pmset.RefreshAllViews(
        molset.RFRefresh | molset.RFColour | molset.RFApply | molset.RFMagnify | molset.RFRotate | molset.RFTransX | molset.RFTransY)