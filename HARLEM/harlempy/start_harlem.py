"""
Here we define start_harlem() function that starts Harlem.
"""

def start_harlem():
    """
    Start Harlem
    """
    import wx
    import harlempy.halib
    import harlempy.molset
    app = wx.App()
    harlempy.molset.StartHarlemApp()
    harlempy.molset.StartHaMainFrameWX()
    
    import harlempy.hapygui_init
    app.MainLoop()

"""
Launch as module, i.e. python3 -m harlempy.start_harlem
"""
if __name__ == "__main__":
    start_harlem()
