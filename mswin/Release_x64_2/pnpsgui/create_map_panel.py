import logging as log
import re
from math import fabs

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
from .utils import FileBB_GetValue, MolSet_RefreshAllViews


class CreateMapsPanel(PNPSBasePanel):
    """
    Panel which help with map creation
    """
    def __init__(self, parent):
        """
        `parent` should be Notebook of  PNPFrame
        """
        super(CreateMapsPanel, self).__init__(parent, pnpgui_wdr.pnpMapCreation, pnpgui_wdr)
        self.script_interruption_is_possible = False

        # Set references to widgets
        # Grid Parameters
        self.NumCtrl_GridSizeX = self.FindWindowById(pnpgui_wdr.IDCpnpGridSizeX)
        self.NumCtrl_GridSizeY = self.FindWindowById(pnpgui_wdr.IDCpnpGridSizeY)
        self.NumCtrl_GridSizeZ = self.FindWindowById(pnpgui_wdr.IDCpnpGridSizeZ)
        self.NumCtrl_GridScale = self.FindWindowById(pnpgui_wdr.IDCpnpGridScale)
        # Ions Parameters
        self.NumCtrl_IonsQ1 = self.FindWindowById(pnpgui_wdr.IDTCpnpIonsQ1)
        self.NumCtrl_IonsQ2 = self.FindWindowById(pnpgui_wdr.IDTCpnpIonsQ2)
        self.NumCtrl_IonsR1 = self.FindWindowById(pnpgui_wdr.IDTCpnpIonsR1)
        self.NumCtrl_IonsR2 = self.FindWindowById(pnpgui_wdr.IDTCpnpIonsR2)
        self.NumCtrl_Dbulk1 = self.FindWindowById(pnpgui_wdr.IDTCpnpIonsDbulk1)
        self.NumCtrl_Dbulk2 = self.FindWindowById(pnpgui_wdr.IDTCpnpIonsDbulk2)
        self.NumCtrl_Cbulk = self.FindWindowById(pnpgui_wdr.IDTCpnpIonsCbulk)
        # Solvent Parameters
        self.NumCtrl_SolDiel = self.FindWindowById(pnpgui_wdr.IDTCpnpSolDiel)
        # Protein Parameters
        self.NumCtrl_ProtDiel = self.FindWindowById(pnpgui_wdr.IDTCpnpProtDiel)
        self.FileBB_CrMap_AtomParFile = self.FindWindowById(pnpgui_wdr.IDFBBpnpCrMapAtomParFile)
        # Membrane Model: Slab with Cylindrical Hole
        self.NumCtrl_MemDiel = self.FindWindowById(pnpgui_wdr.IDTCpnpMemDiel)
        self.NumCtrl_MemZ1 = self.FindWindowById(pnpgui_wdr.IDTCpnpMemZ1)
        self.NumCtrl_MemZ2 = self.FindWindowById(pnpgui_wdr.IDTCpnpMemZ2)
        self.NumCtrl_MemX = self.FindWindowById(pnpgui_wdr.IDTCpnpMemX)
        self.NumCtrl_MemY = self.FindWindowById(pnpgui_wdr.IDTCpnpMemY)
        self.NumCtrl_MemR1 = self.FindWindowById(pnpgui_wdr.IDTCpnpMemR1)
        # Initial Diffusion Model
        self.RB_InitDiffMod = self.FindWindowById(pnpgui_wdr.IDRB_InitDiffMod)
        self.RB_InitDiffMod.SetSelection(1)
        # Other parameters
        self.NumCtrl_RprobeDiel = self.FindWindowById(pnpgui_wdr.IDTCpnpRprobeDiel)
        self.NumCtrl_RprobeDiff = self.FindWindowById(pnpgui_wdr.IDTCpnpRprobeDiff)
        # Output
        assert self.FBB_ScriptFileName is not None
        assert self.FBB_SysTopFileName_Out is not None
        assert self.FBB_DiffFileName_Out is not None
        assert self.FBB_SRFileName_Out is not None
        # Process Control
        assert self.Btn_Save is not None
        assert self.Btn_Run is not None
        assert self.Btn_Stop is not None
        # Results viewing
        assert self.Btn_ViewStaticCharge is not None
        assert self.Btn_ViewDiel is not None
        assert self.Btn_ViewSR is not None
        assert self.Btn_ViewPlainDiff is not None
        assert self.Btn_ViewDiff is not None
        # assert self.Btn_ViewPot is not None
        # assert self.Btn_ViewConc is not None

        # To validate user entry
        self.Validation_NumCtrls = [
            self.NumCtrl_GridSizeX,
            self.NumCtrl_GridSizeY,
            self.NumCtrl_GridSizeZ,
            self.NumCtrl_GridScale,
            self.NumCtrl_IonsQ1,
            self.NumCtrl_IonsQ2,
            self.NumCtrl_IonsR1,
            self.NumCtrl_IonsR2,
            self.NumCtrl_Dbulk1,
            self.NumCtrl_Dbulk2,
            self.NumCtrl_Cbulk,
            self.NumCtrl_SolDiel,
            self.NumCtrl_ProtDiel,
            self.NumCtrl_MemDiel,
            self.NumCtrl_MemZ1,
            self.NumCtrl_MemZ2,
            self.NumCtrl_MemX,
            self.NumCtrl_MemY,
            self.NumCtrl_MemR1,
            self.NumCtrl_RprobeDiel,
            self.NumCtrl_RprobeDiff]

        # FBB_InputFiles/FBB_OutputFiles are list of (FileBrowseButton, NoneIsAllowed) values
        self.FBB_InputFiles = [(self.FileBB_CrMap_AtomParFile, True)]

        self.FBB_OutputFiles = [
            (self.FBB_ScriptFileName, False),
            (self.FBB_SysTopFileName_Out, False),
            (self.FBB_DiffFileName_Out, True),
            (self.FBB_SRFileName_Out, True)
        ]

        # bind events
        self.Bind(wx.EVT_BUTTON, self.OnDrawSimBox, id=pnpgui_wdr.IDBpnpDrawSimBox)
        self.Bind(wx.EVT_BUTTON, self.OnDrawMembrane, id=pnpgui_wdr.IDBpnpDrawMem)
        self.Bind(wx.EVT_RADIOBOX, self.OnDiffusionModeChange, id=pnpgui_wdr.IDRB_InitDiffMod)

        # Other internal values
        self.DrawnSimBox_GridSize = None
        self.DrawnSimBox_GridScale = None
        self.DrawnSimBox_SimBox = None

        self.DrawnMembrane_GridSize = None
        self.DrawnMembrane_GridScale = None
        self.DrawnMembrane_MemZ = None
        self.DrawnMembrane_MemXY = None
        self.DrawnMembrane_MemR = None
        self.DrawnMembrane_Mem = None

    def OnDrawSimBox(self, _):
        if self.pmset is None:
            log.error("Can not draw sim box, there is no associated molset")
            return

        GridSize = [self.NumCtrl_GridSizeX.GetValue(),
                    self.NumCtrl_GridSizeY.GetValue(),
                    self.NumCtrl_GridSizeZ.GetValue()]
        GridScale = self.NumCtrl_GridScale.GetValue()
        if self.DrawnSimBox_SimBox is None:
            # Draw box
            ircenter0 = [GridSize[0] / 2, GridSize[1] / 2, GridSize[2] / 2]

            fr0 = [-ircenter0[0] / GridScale,
                   -ircenter0[1] / GridScale,
                   -ircenter0[2] / GridScale]
            fr1 = [(GridSize[0] - 1 - ircenter0[0]) / GridScale,
                   (GridSize[1] - 1 - ircenter0[1]) / GridScale,
                   (GridSize[2] - 1 - ircenter0[2]) / GridScale]

            SimBox = molset.BoxObj3D("PNPSimulationBox%s" % (self.pmset.GetName()), 128, 128, 255, fr0[0], fr0[1],
                                     fr0[2], fr1[0], fr1[1], fr1[2])
            SimBox.thisown = 0
            self.pmset.AddObject3D(SimBox)
            MolSet_RefreshAllViews(self.pmset)

            self.DrawnSimBox_GridSize = GridSize
            self.DrawnSimBox_GridScale = GridScale
            self.DrawnSimBox_SimBox = SimBox

        else:
            ValuesChanged = False
            if GridSize[0] != self.DrawnSimBox_GridSize[0]:
                ValuesChanged = True
            if GridSize[1] != self.DrawnSimBox_GridSize[1]:
                ValuesChanged = True
            if GridSize[2] != self.DrawnSimBox_GridSize[2]:
                ValuesChanged = True
            if fabs(GridScale - self.DrawnSimBox_GridScale) > 0.001:
                ValuesChanged = True

            if ValuesChanged:
                # redraw box
                ircenter0 = [GridSize[0] / 2, GridSize[1] / 2, GridSize[2] / 2]

                fr0 = [-ircenter0[0] / GridScale,
                       -ircenter0[1] / GridScale,
                       -ircenter0[2] / GridScale]
                fr1 = [(GridSize[0] - 1 - ircenter0[0]) / GridScale,
                       (GridSize[1] - 1 - ircenter0[1]) / GridScale,
                       (GridSize[2] - 1 - ircenter0[2]) / GridScale]

                self.DrawnSimBox_SimBox.SetBox(fr0[0], fr0[1], fr0[2], fr1[0], fr1[1], fr1[2])
                MolSet_RefreshAllViews(self.pmset)

                self.DrawnSimBox_GridSize = GridSize
                self.DrawnSimBox_GridScale = GridScale
            else:
                # delete box
                self.pmset.DeleteObject3D(self.DrawnSimBox_SimBox)
                MolSet_RefreshAllViews(self.pmset)
                self.DrawnSimBox_GridSize = None
                self.DrawnSimBox_GridScale = None
                self.DrawnSimBox_SimBox = None

    def OnDrawMembrane(self, _):
        if self.pmset is None:
            log.error("Can not draw sim box, there is no associated molset")
            return
        GridSize = [self.NumCtrl_GridSizeX.GetValue(), self.NumCtrl_GridSizeY.GetValue(),
                    self.NumCtrl_GridSizeZ.GetValue()]
        GridScale = self.NumCtrl_GridScale.GetValue()

        MemZ = [self.NumCtrl_MemZ1.GetValue(), self.NumCtrl_MemZ2.GetValue()]
        MemXY = [self.NumCtrl_MemX.GetValue(), self.NumCtrl_MemY.GetValue()]
        MemR = [self.NumCtrl_MemR1.GetValue(), max([GridSize[0] / GridScale, GridSize[1] / GridScale])]

        if self.DrawnMembrane_Mem is None:
            # Draw Membrane
            ircenter0 = [GridSize[0] / 2, GridSize[1] / 2, GridSize[2] / 2]

            fr0 = [-ircenter0[0] / GridScale,
                   -ircenter0[1] / GridScale,
                   -ircenter0[2] / GridScale]
            fr1 = [(GridSize[0] - 1 - ircenter0[0]) / GridScale,
                   (GridSize[1] - 1 - ircenter0[1]) / GridScale,
                   (GridSize[2] - 1 - ircenter0[2]) / GridScale]

            Mem = molset.TubeObj3D("TubeObj3D%s" % (self.pmset.GetName()), 255, 128, 128)

            Mem.SetBox(fr0[0], fr0[1], fr0[2], fr1[0], fr1[1], fr1[2])
            Mem.SetTube3d(MemZ[0], MemZ[1], MemXY[0], MemXY[1], MemR[0], MemR[1])
            Mem.thisown = 0
            self.pmset.AddObject3D(Mem)
            MolSet_RefreshAllViews(self.pmset)

            self.DrawnMembrane_GridSize = GridSize
            self.DrawnMembrane_GridScale = GridScale
            self.DrawnMembrane_MemZ = MemZ
            self.DrawnMembrane_MemXY = MemXY
            self.DrawnMembrane_MemR = MemR
            self.DrawnMembrane_Mem = Mem
        else:
            ValuesChanged = False
            if GridSize[0] != self.DrawnMembrane_GridSize[0]:
                ValuesChanged = True
            if GridSize[1] != self.DrawnMembrane_GridSize[1]:
                ValuesChanged = True
            if GridSize[2] != self.DrawnMembrane_GridSize[2]:
                ValuesChanged = True
            if fabs(GridScale - self.DrawnMembrane_GridScale) > 0.001:
                ValuesChanged = True
            if fabs(MemZ[0] - self.DrawnMembrane_MemZ[0]) > 0.001:
                ValuesChanged = True
            if fabs(MemZ[1] - self.DrawnMembrane_MemZ[1]) > 0.001:
                ValuesChanged = True
            if fabs(MemR[0] - self.DrawnMembrane_MemR[0]) > 0.001:
                ValuesChanged = True
            if fabs(MemR[1] - self.DrawnMembrane_MemR[1]) > 0.001:
                ValuesChanged = True
            if fabs(MemXY[0] - self.DrawnMembrane_MemXY[0]) > 0.001:
                ValuesChanged = True
            if fabs(MemXY[1] - self.DrawnMembrane_MemXY[1]) > 0.001:
                ValuesChanged = True

            if ValuesChanged:
                # redraw box
                ircenter0 = [GridSize[0] / 2, GridSize[1] / 2, GridSize[2] / 2]

                fr0 = [-ircenter0[0] / GridScale,
                       -ircenter0[1] / GridScale,
                       -ircenter0[2] / GridScale]
                fr1 = [(GridSize[0] - 1 - ircenter0[0]) / GridScale,
                       (GridSize[1] - 1 - ircenter0[1]) / GridScale,
                       (GridSize[2] - 1 - ircenter0[2]) / GridScale]
                self.DrawnMembrane_Mem.SetBox(fr0[0], fr0[1], fr0[2], fr1[0], fr1[1], fr1[2])
                self.DrawnMembrane_Mem.SetTube3d(MemZ[0], MemZ[1], MemXY[0], MemXY[1], MemR[0], MemR[1])
                MolSet_RefreshAllViews(self.pmset)

                self.DrawnMembrane_GridSize = GridSize
                self.DrawnMembrane_GridScale = GridScale
                self.DrawnMembrane_MemZ = MemZ
                self.DrawnMembrane_MemXY = MemXY
                self.DrawnMembrane_MemR = MemR
            else:
                # delete box
                self.pmset.DeleteObject3D(self.DrawnMembrane_Mem)
                MolSet_RefreshAllViews(self.pmset)
                self.DrawnMembrane_GridSize = None
                self.DrawnMembrane_GridScale = None
                self.DrawnMembrane_MemZ = None
                self.DrawnMembrane_MemXY = None
                self.DrawnMembrane_MemR = None
                self.DrawnMembrane_Mem = None

    def OnDiffusionModeChange(self, _):
        if self.RB_InitDiffMod.GetSelection() == 0:
            DiffusionModel = "Plain"
            self.FBB_DiffFileName_Out.Enable(False)
        else:
            DiffusionModel = "MD"
            self.FBB_DiffFileName_Out.Enable(True)
        log.debug("CreateMaps_OnInitDiffMod %s %s", self.RB_InitDiffMod.GetSelection(), DiffusionModel)

    def ValidateParameters(self):
        return self._ValidateParameters()

    def GenerateRunScript(self):
        # Get User Values
        GridSize = (self.NumCtrl_GridSizeX.GetValue(), self.NumCtrl_GridSizeY.GetValue(),
                    self.NumCtrl_GridSizeZ.GetValue())
        GridScale = self.NumCtrl_GridScale.GetValue()

        Qions = (self.NumCtrl_IonsQ1.GetValue(), self.NumCtrl_IonsQ2.GetValue())
        Dbulk = (self.NumCtrl_Dbulk1.GetValue(), self.NumCtrl_Dbulk2.GetValue())
        Rions = (self.NumCtrl_IonsR1.GetValue(), self.NumCtrl_IonsR2.GetValue())
        Cbulk = self.NumCtrl_Cbulk.GetValue()

        SolDiel = self.NumCtrl_SolDiel.GetValue()
        ProtDiel = self.NumCtrl_ProtDiel.GetValue()
        MemDiel = self.NumCtrl_MemDiel.GetValue()

        MemZ = (self.NumCtrl_MemZ1.GetValue(), self.NumCtrl_MemZ2.GetValue())
        MemXY = (self.NumCtrl_MemX.GetValue(), self.NumCtrl_MemY.GetValue())
        MemR = (self.NumCtrl_MemR1.GetValue(), max([GridSize[0] / GridScale, GridSize[1] / GridScale]))

        AtomParFile = str(self.FileBB_CrMap_AtomParFile.GetValue())

        RprobeDiel = self.NumCtrl_RprobeDiel.GetValue()
        RprobeDiff = self.NumCtrl_RprobeDiff.GetValue()

        if self.RB_InitDiffMod.GetSelection() == 0:
            DiffusionModel = "Plain"
        else:
            DiffusionModel = "MD"

        SysTopFileName = FileBB_GetValue(self.FBB_SysTopFileName_Out)
        DiffFileName = FileBB_GetValue(self.FBB_DiffFileName_Out)
        SRFileName = FileBB_GetValue(self.FBB_SRFileName_Out)

        cmds = add_script_header()
        cmds += format_cmd(
            "# Create Continuum representation\n"
            "contworld = pnps.ContWorld(\n"
            "    GridSize={GridSize},\n"
            "    GridScale={GridScale},\n"
            "    Qions={Qions}\n"
            ")\n\n",
            GridSize="[%d, %d, %d]" % GridSize,
            GridScale=GridScale,
            Qions="[%s, %s]" % Qions)

        # @todo add BC setup
        cmds += format_cmd(
            "# setup builder\n"
            "builder = pnps.GetWorldBuilder(\n"
            "    DielConstBulk={SolDiel},\n"
            "    DiffCoefBulk={Dbulk},\n"
            "    Cbulk={Cbulk},\n"
            "    RprobeDiel={RprobeDiel},\n"
            "    RprobeDiff={RprobeDiff},\n"
            "    Rions={Rions},\n"
            "    DiffusionModel=\"{DiffusionModel}\",\n"
            "    BoundaryCondition=\"ZeroBC\",\n"
            "    MakeDielectricMap=True,\n"
            "    MakeChargeMap=True,\n"
            "    MakeDiffusionMap=True,\n"
            "    MakeConcentrationMap=True,\n"
            "    MakeSoftRepultionMap=True\n"
            ")\n\n",
            SolDiel=SolDiel,
            Dbulk="[%s, %s]" % Dbulk,
            Cbulk=Cbulk,
            RprobeDiel=RprobeDiel,
            RprobeDiff=RprobeDiff,
            Rions="[%s, %s]" % Rions,
            DiffusionModel=DiffusionModel
        )

        # if AtomParFile != "None":
        if re.search(".*atmpar$", AtomParFile) is not None:
            cmds += format_cmd(
                "# add atoms\n"
                "builder.addAtoms(\n"
                "    DielConst={ProtDiel},\n"
                "    AtomicParam=\"{AtomParFile}\"\n"
                ")\n\n",
                ProtDiel=ProtDiel,
                AtomParFile=AtomParFile)

        # @todo add membrane options
        cmds += format_cmd(
            "# add membrane with hole, note that second value in R is large\n"
            "builder.addTube(\n"
            "    DielConst={MemDiel},\n"
            "    x={MemX}, y={MemY},\n"
            "    z={MemZ},\n"
            "    R={MemR}\n"
            ")\n\n",
            MemDiel=MemDiel,
            MemX=MemXY[0], MemY=MemXY[1],
            MemZ="[%s, %s]" % MemZ,
            MemR="[%s, %s]" % MemR
        )

        cmds += format_cmd(
            "# Do building\n"
            "builder.BuildContWorld(contworld)\n\n")

        cmds += format_cmd("# Write created maps to file system\n")

        cmds += format_cmd(
            "contworld.WriteNodeIndexing(\"{SysTopFileName}\")\n",
            SysTopFileName=SysTopFileName,
        )
        if SRFileName is not None:
            cmds += format_cmd(
                "contworld.WritePMF(\"{SRFileName}\")\n",
                SRFileName=SRFileName
            )

        if DiffusionModel != "Plain" and DiffFileName is not None:
            cmds += format_cmd(
                "contworld.WriteDiffusion(\"{DiffFileName}\")\n",
                DiffFileName=DiffFileName)

        cmds += format_cmd(
            "# Clean-up\n"
            "del builder\n"
            "del contworld\n"
        )
        return cmds
