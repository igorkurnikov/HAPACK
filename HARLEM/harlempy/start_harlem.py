"""
Here we define start_harlem() function that starts Harlem.
"""

def start_harlem():
    """
    Start Harlem
    """
    import wx
    import molset

    app = wx.App()
    molset.StartHarlemApp()
    molset.StartHaMainFrameWX()
    
    app.MainLoop()

"""
Launch as module, i.e. python3 -m harlempy.start_harlem
"""
if __name__ == "__main__":
    start_harlem()
