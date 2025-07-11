import wx
from . import hapygui
from . import wxexamples

def print_test(e):
    print("TEST HO HO ")


def list_box_test(e):
    frame = hapygui.GetMainFrame()
    ex = wxexamples.ExampleListBox(frame)
    ex.Show()
    print("LIST BOX TEST")

def grid_test(e):
    frame = hapygui.GetMainFrame()
    ex = wxexamples.ExampleGrid2(frame)
    ex.Show()
    print("LIST BOX TEST")

def open_python_window(e):
    frm_main = hapygui.GetMainFrame()
    frm_py = wx.py.frame.Frame(frm_main)
    frm_py.Show(True)
    win = wx.py.shell.Shell(frm_py,-1)

def expand_file_menu():
    menu_python_window = hapygui.GetFileMenu().Insert(5,hapygui.NewID(), "Interactive Python Window", "Interactive Python Window") 
    hapygui.GetMainFrame().Bind(wx.EVT_MENU, open_python_window, menu_python_window)

def add_tests():
    menu_print_test = hapygui.GetTestMenu().Append(hapygui.NewID(), "PRINT TEST", "PRINT TEST wxPython example") 
    hapygui.GetMainFrame().Bind(wx.EVT_MENU, print_test, menu_print_test)
    menu_init_list_panel = hapygui.GetTestMenu().Append(hapygui.NewID(), "wxListBox test", "wxListBox test") 
    hapygui.GetMainFrame().Bind(wx.EVT_MENU, list_box_test, menu_init_list_panel)
    menu_init_grid_panel = hapygui.GetTestMenu().Append(hapygui.NewID(), "wxGrid test", "wxGrid test") 
    hapygui.GetMainFrame().Bind(wx.EVT_MENU, grid_test, menu_init_grid_panel)




