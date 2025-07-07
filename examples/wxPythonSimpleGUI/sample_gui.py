"""
Example of combining wxPython GUI with main Harlem pure wxWidgets GUI.
Menu item added to existing menu and simple dialog is created on that 
menu item click.
"""
import wx
from hapygui import NewID,GetMainFrame,GetTestMenu

def OnSimpleGUI(e):
    # A message dialog box with an OK button. wx.OK is a standard ID in wxWidgets.
    dlg = wx.MessageDialog(GetMainFrame(), "This is a message", "Sample GUI", wx.OK)
    dlg.ShowModal() # Show it
    dlg.Destroy() # finally destroy it when finished.

# Setting up the menu.
# Create menu item
menuSimpleGUI = GetTestMenu().Append(NewID(), "Sample GUI","Call a sample message")
# Bind menu item event and OnSimpleGUI function
GetMainFrame().Bind(wx.EVT_MENU, OnSimpleGUI, menuSimpleGUI)

