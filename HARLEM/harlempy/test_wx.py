import hapygui
import wx

def print_test(e):
    print("TEST")

def add_tests():
    menu_print_test = hapygui.GetTestMenu().Append(hapygui.NewID(), "PRINT TEST", "PRINT TEST wxPython example") 
    hapygui.GetMainFrame().Bind(wx.EVT_MENU, print_test, menu_print_test)