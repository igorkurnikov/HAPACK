"""
Here we define start_harlem() function that starts Harlem GUI.
"""
from . import test_wx

def start_harlem():
    """
    Start Harlem
    """
    import wx
    import wx.py
    import molset

    app = wx.App()
    molset.StartHarlemApp()
    molset.StartHaMainFrameWX()
    frm_main = molset.harlempy.GetMainFrame()

    test_wx.expand_file_menu()
    test_wx.add_tests()
    app.MainLoop()

"""
Launch as module, i.e. python3 -m harlempy
"""
if __name__ == "__main__":
    start_harlem()
