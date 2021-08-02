# import logging as log
# import re
# from math import fabs, isnan, isinf

import wx

try:
    import pnps
except ImportError:
    pnps = None

try:
    from harlempy import molset
except ImportError:
    molset = None

from . import pnpgui_wdr
from .pnpsbasepanel import PNPSBasePanel, add_script_header, format_cmd
from .utils import FileBB_GetValue


class PBSRPanel(PNPSBasePanel):
    """
    Panel which help with PB-SR Calculation
    """
    def __init__(self, parent):
        """
        `parent` should be Notebook of  PNPFrame
        """
        super(PBSRPanel, self).__init__(parent, pnpgui_wdr.pnpPBSR, pnpgui_wdr)
        self.script_interruption_is_possible = True

        # Set references to widgets
        # Grid Parameters
        self.TextCtrl_PBSR_Tolerance = self.FindWindowById(pnpgui_wdr.IDTCpnpPBSR_Tolerance)
        self.NumCtrl_PBSR_Niter = self.FindWindowById(pnpgui_wdr.IDTCpnpPBSR_Niter)
        self.NumCtrl_PBSR_Relax = self.FindWindowById(pnpgui_wdr.IDTCpnpPBSR_Relax)
        # Output
        assert self.FBB_SysTopFileName_In is not None
        assert self.FBB_SRFileName_In is not None
        assert self.FBB_ScriptFileName is not None
        assert self.FBB_PotFileName_Out is not None
        assert self.FBB_ConcFileName_Out is not None
        # Process Control
        assert self.Btn_Save is not None
        assert self.Btn_Run is not None
        assert self.Btn_Stop is not None
        assert self.Btn_Interrupt is not None
        # Results viewing
        assert self.Btn_ViewPot is not None
        assert self.Btn_ViewConc is not None

        # To validate user entry
        self.Validation_NumCtrls = [
            self.NumCtrl_PBSR_Niter,
            self.NumCtrl_PBSR_Relax
        ]
        # FBB_InputFiles/FBB_OutputFiles are list of (FileBrowseButton, NoneIsAllowed) values
        self.FBB_InputFiles = [
            (self.FBB_SysTopFileName_In, False),
            (self.FBB_SRFileName_In, True),
        ]

        self.FBB_OutputFiles = [
            (self.FBB_ScriptFileName, False),
            (self.FBB_PotFileName_Out, True),
            (self.FBB_ConcFileName_Out, True)
        ]

        # Other internal values

    def ValidateParameters(self):
        if not self._ValidateParameters():
            return False

        TolValCorrect = True
        try:
            strTol = str(self.TextCtrl_PBSR_Tolerance.GetValue())
            Tol = float(strTol)
            if Tol < 0.0:
                TolValCorrect = False
        except ValueError:
            TolValCorrect = False

        if not TolValCorrect:
            dlg = wx.MessageDialog(self, "Convergence(Tolerance) value is incorrect", "Error In Input Values",
                                   wx.OK)
            dlg.ShowModal()  # Show it
            dlg.Destroy()  # finally destroy it when finished.
            return False

        return True

    def GenerateRunScript(self):
        SysTopFileNameIn = FileBB_GetValue(self.FBB_SysTopFileName_In)
        SRFileNameIn = FileBB_GetValue(self.FBB_SRFileName_In)

        MaxIterations = int(self.NumCtrl_PBSR_Niter.GetValue())
        Tolerance = float(self.TextCtrl_PBSR_Tolerance.GetValue())
        Relaxation = float(self.NumCtrl_PBSR_Relax.GetValue())

        PotFileNameOut = FileBB_GetValue(self.FBB_PotFileName_Out)
        ConcFileNameOut = FileBB_GetValue(self.FBB_ConcFileName_Out)

        cmds = add_script_header()

        cmds += format_cmd(
            "# Load ContWorld from Maps\n"
            "contworld = pnps.LoadContWorld(\n"
            "    SysTop=\"{SysTop}\""
            "{PMF}"
            ")\n"
            "assert contworld is not None\n\n",
            SysTop=SysTopFileNameIn,
            PMF="" if SRFileNameIn is None else ",\n    PMF=\""+SRFileNameIn+"\"\n"
        )

        cmds += format_cmd(
            "# Do PBSR\n"
            "pnps.SolvePBSR(\n"
            "    contworld, MaxIterations={MaxIterations}, Tolerance={Tolerance}, Relaxation={Relaxation}\n"
            ")\n\n",
            MaxIterations=MaxIterations, Tolerance=Tolerance, Relaxation=Relaxation
        )

        cmds += format_cmd("# save results\n")
        if PotFileNameOut is not None:
            cmds += format_cmd("contworld.WritePotential(\"{}\")\n", PotFileNameOut)
        if ConcFileNameOut is not None:
            cmds += format_cmd("pnps.PNPUtil.ConvertPBLJresultsToDynamicCharge(contworld)\n")
            cmds += format_cmd("contworld.WriteDynamicCharge(\"{}\")\n", ConcFileNameOut)

        return cmds

    def AnalyseOutput(self, return_code):
        # E = contworld.SystemEnergy
        # del contworld
        # print "done running PB-SR"
        # print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
        #
        # if E > 1.0E13 or isnan(E) or isinf(E):
        #     dlg = wx.MessageDialog(self, "PBSR has diverged, try smaller relaxation", "Error In PBSR Execution",
        #                            wx.OK)
        #     dlg.ShowModal()  # Show it
        #     dlg.Destroy()  # finally destroy it when finished.
        #     return False
        self._AnalyseOutput(return_code, None)
