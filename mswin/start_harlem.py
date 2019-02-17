import sys

print(sys.path)
import wx
print("ASDASDAS")
import harlempy.molset
print("ASDASDAS2")
from harlempy.halib import *
print("ASDASDAS3")
from harlempy.molset import *
print("ASDASDAS4")
app = wx.App()
harlempy.molset.StartHarlemApp()
harlempy.molset.StartHaMainFrameWX()
app.MainLoop()



#GetCurMolSet' is not defined
#'mset_c' is not defined
