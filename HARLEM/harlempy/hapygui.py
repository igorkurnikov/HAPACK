""""helping routins for making wxPython GUI for Harlem"""
import wx

last_issued_id=26000

def NewID():
    """Return a new ID, hopefully unique, for use with wxPython, would not conflict with main wxWidgets GUI"""
    global last_issued_id
    last_issued_id=last_issued_id+1
    return last_issued_id


def GetMainFrame():
    return wx.FindWindowByName("HARLEM")


def GetMainMenu():
    pyMainFrame=GetMainFrame()
    return pyMainFrame.GetMenuBar()


def GetEditMenu():
    pyMainMenu=GetMainMenu()
    return pyMainMenu.GetMenu(pyMainMenu.FindMenu("Edit"))


def GetApplicationsMenu():
    pyMainMenu=GetMainMenu()
    return pyMainMenu.GetMenu(pyMainMenu.FindMenu("Applications"))


def GetTestMenu():
    pyMainMenu=GetMainMenu()
    return pyMainMenu.GetMenu(pyMainMenu.FindMenu("Test"))
