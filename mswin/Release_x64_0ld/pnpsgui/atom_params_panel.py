import os
import logging as log
import wx
import wx.grid

try:
    import pnps
except ImportError:
    pnps = None

try:
    from harlempy import molset
except ImportError:
    molset = None

from . import pnpgui_wdr
from .pnpsbasepanel import PNPSBasePanel
from .utils import FileBB_GetValue


class AtomicParametersTable(wx.grid.PyGridTableBase):
    def __init__(self, m_wxGrid,
                 natoms=0,
                 AtomsSerNo=None,
                 ResidueSerNo=None,
                 AtomsCharge=None,
                 AtomsRadius=None,
                 AtomsName=None,
                 ResidueName=None,
                 AtomsCoorX=None,
                 AtomsCoorY=None,
                 AtomsCoorZ=None,
                 AtomsIER1=None,
                 AtomsIER2=None,
                 SR_A_K=None,
                 SR_N_K=None,
                 SR_A_Cl=None,
                 SR_N_Cl=None):
        wx.grid.PyGridTableBase.__init__(self)
        self.m_wxGrid = m_wxGrid
        self.natoms = natoms
        self.AtomsSerNo = AtomsSerNo
        self.ResidueSerNo = ResidueSerNo
        self.AtomsCharge = AtomsCharge
        self.AtomsRadius = AtomsRadius
        self.AtomsName = AtomsName
        self.ResidueName = ResidueName
        self.AtomsCoorX = AtomsCoorX
        self.AtomsCoorY = AtomsCoorY
        self.AtomsCoorZ = AtomsCoorZ
        self.AtomsIER1 = AtomsIER1
        self.AtomsIER2 = AtomsIER2
        self.SR_A_K = SR_A_K
        self.SR_N_K = SR_N_K
        self.SR_A_Cl = SR_A_Cl
        self.SR_N_Cl = SR_N_Cl
        self.currentRows = natoms
        self.currentCols = 15
        self.UserChangedSomeValues = False

    def SetAtomicParameters(self, natoms=0,
                            AtomsSerNo=None,
                            ResidueSerNo=None,
                            AtomsCharge=None,
                            AtomsRadius=None,
                            AtomsName=None,
                            ResidueName=None,
                            AtomsCoorX=None,
                            AtomsCoorY=None,
                            AtomsCoorZ=None,
                            AtomsIER1=None,
                            AtomsIER2=None,
                            SR_A_K=None,
                            SR_N_K=None,
                            SR_A_Cl=None,
                            SR_N_Cl=None):

        self.natoms = natoms
        self.AtomsSerNo = AtomsSerNo
        self.ResidueSerNo = ResidueSerNo
        self.AtomsCharge = AtomsCharge
        self.AtomsRadius = AtomsRadius
        self.AtomsName = AtomsName
        self.ResidueName = ResidueName
        self.AtomsCoorX = AtomsCoorX
        self.AtomsCoorY = AtomsCoorY
        self.AtomsCoorZ = AtomsCoorZ
        self.AtomsIER1 = AtomsIER1
        self.AtomsIER2 = AtomsIER2
        self.SR_A_K = SR_A_K
        self.SR_N_K = SR_N_K
        self.SR_A_Cl = SR_A_Cl
        self.SR_N_Cl = SR_N_Cl
        self.ResetView()
        self.currentRows = self.GetNumberRows()
        self.currentCols = self.GetNumberCols()
        self.UserChangedSomeValues = False

    def DeleteRows(*args, **kwargs):
        """DeleteRows(self, int pos=0, int numRows=1, bool updateLabels=True) -> bool"""
        return
        # return _grid.Grid_DeleteRows(*args, **kwargs)

    def ResetView(self):
        """Trim/extend the control's rows and update all values"""
        self.m_wxGrid.BeginBatch()
        for current, new, delmsg, addmsg in [
            (self.currentRows, self.GetNumberRows(), wx.grid.GRIDTABLE_NOTIFY_ROWS_DELETED,
             wx.grid.GRIDTABLE_NOTIFY_ROWS_APPENDED),
            (self.currentCols, self.GetNumberCols(), wx.grid.GRIDTABLE_NOTIFY_COLS_DELETED,
             wx.grid.GRIDTABLE_NOTIFY_COLS_APPENDED),
        ]:
            if new < current:
                msg = wx.grid.GridTableMessage(
                    self,
                    delmsg,
                    new,  # position
                    current - new,
                )
                self.m_wxGrid.ProcessTableMessage(msg)
            elif new > current:
                msg = wx.grid.GridTableMessage(
                    self,
                    addmsg,
                    new - current
                )
                self.m_wxGrid.ProcessTableMessage(msg)
        self.UpdateValues()
        self.m_wxGrid.EndBatch()

        # The scroll bars aren't resized (at least on windows)
        # Jiggling the size of the window rescales the scrollbars
        h, w = self.m_wxGrid.GetSize()
        self.m_wxGrid.SetSize((h + 1, w))
        self.m_wxGrid.SetSize((h, w))
        self.m_wxGrid.ForceRefresh()

    def UpdateValues(self):
        """Update all displayed values"""
        msg = wx.grid.GridTableMessage(self, wx.grid.GRIDTABLE_REQUEST_VIEW_GET_VALUES)
        self.m_wxGrid.ProcessTableMessage(msg)

    def GetRowLabelValue(self, row):
        return str(row + 1)

    def GetColLabelValue(self, col):
        if col == 0:
            return "ResNum"
        elif col == 1:
            return "ResName"
        elif col == 2:
            return "AtmNum"
        elif col == 3:
            return "AtmName"
        elif col == 4:
            return "x"
        elif col == 5:
            return "y"
        elif col == 6:
            return "z"
        elif col == 7:
            return "Q"
        elif col == 8:
            return "Rdiel"
        elif col == 9:
            return "Rier1"
        elif col == 10:
            return "Rier2"
        elif col == 11:
            return "SR-A1"
        elif col == 12:
            return "SR-N1"
        elif col == 13:
            return "SR-A2"
        elif col == 14:
            return "SR-N2"
        else:
            return "Unknown"

    # these five are the required methods
    def GetNumberRows(self):
        return self.natoms

    def GetNumberCols(self):
        return 15

    def IsEmptyCell(self, row, col):
        if col > 14: return False
        if row > self.natoms: return False
        return True

    def GetValue(self, row, col):
        if self.natoms > 0:
            if col == 0:
                return "%d" % (self.ResidueSerNo[row])
            elif col == 1:
                return self.ResidueName[row]
            elif col == 2:
                return "%d" % (self.AtomsSerNo[row])
            elif col == 3:
                return self.AtomsName[row]
            elif col == 4:
                return "%.4f" % (self.AtomsCoorX[row])
            elif col == 5:
                return "%.4f" % (self.AtomsCoorY[row])
            elif col == 6:
                return "%.4f" % (self.AtomsCoorZ[row])
            elif col == 7:
                return "%.4f" % (self.AtomsCharge[row])
            elif col == 8:
                return "%.4f" % (self.AtomsRadius[row])
            elif col == 9:
                return "%.4f" % (self.AtomsIER1[row])
            elif col == 10:
                return "%.4f" % (self.AtomsIER2[row])
            elif col == 11:
                return "%.4f" % (self.SR_A_K[row])
            elif col == 12:
                return "%.4f" % (self.SR_N_K[row])
            elif col == 13:
                return "%.4f" % (self.SR_A_Cl[row])
            elif col == 14:
                return "%.4f" % (self.SR_N_Cl[row])
            else:
                return "None"
        else:
            return "None"

    def SetValue(self, row, col, value):
        if self.natoms > 0:
            # if col==0: return "%d"%(self.ResidueSerNo[row])
            # elif col==1: return self.ResidueName[row].c_str()
            # elif col== 2: return "%d"%(self.AtomsSerNo[row])
            # elif col== 3: return self.AtomsName[row].c_str()

            if col == 4:
                self.AtomsCoorX[row] = float(value)
            elif col == 5:
                self.AtomsCoorY[row] = float(value)
            elif col == 6:
                self.AtomsCoorZ[row] = float(value)
            elif col == 7:
                self.AtomsCharge[row] = float(value)
            elif col == 8:
                self.AtomsRadius[row] = float(value)
            elif col == 9:
                self.AtomsIER1[row] = float(value)
            elif col == 10:
                self.AtomsIER2[row] = float(value)
            elif col == 11:
                self.SR_A_K[row] = float(value)
            elif col == 12:
                self.SR_N_K[row] = float(value)
            elif col == 13:
                self.SR_A_Cl[row] = float(value)
            elif col == 14:
                self.SR_N_Cl[row] = float(value)

            if 3 < col < 15:
                self.UserChangedSomeValues = True


class AtomParamsPanel(PNPSBasePanel):
    """
    Panel which help with PB-SR Calculation
    """
    def __init__(self, parent):
        """
        `parent` should be Notebook of  PNPFrame
        """
        super(AtomParamsPanel, self).__init__(parent, pnpgui_wdr.pnpPrep, pnpgui_wdr)
        self.pnpmod = None if self.pmset is None else self.pmset.GetPNPMod(True)

        # Set references to widgets
        self.ChoiceFF4QRdiel = self.FindWindowById(pnpgui_wdr.IDCpnpFF4QRdiel)
        self.ChoiceFF4SR = self.FindWindowById(pnpgui_wdr.IDCpnpFF4SR)

        self.NumCtrl_AddRion1 = self.FindWindowById(pnpgui_wdr.IDTCpnpAddRion1)
        self.NumCtrl_AddRion2 = self.FindWindowById(pnpgui_wdr.IDTCpnpAddRion2)

        self.GridAtomParam = self.FindWindowById(pnpgui_wdr.IDGpnpAtomParam)

        self.FBB_AtmParFileName_Out = self.FindWindowById(pnpgui_wdr.IDTCpnpAtmParFileName)
        # Output
        assert self.FBB_AtmParFileName_Out is not None
        # Process Control
        assert self.Btn_Run is None

        # To validate user entry
        self.Validation_NumCtrls = [
            self.NumCtrl_AddRion1,
            self.NumCtrl_AddRion2
        ]
        # FBB_InputFiles/FBB_OutputFiles are list of (FileBrowseButton, NoneIsAllowed) values
        self.FBB_InputFiles = [
        ]

        self.FBB_OutputFiles = [
            (self.FBB_AtmParFileName_Out, False)
        ]

        # bind events
        self.Bind(wx.EVT_BUTTON, self.OnSetParam, id=pnpgui_wdr.IDBpnpSetParam)
        self.Bind(wx.EVT_BUTTON, self.OnSaveAgainAtomParam, id=pnpgui_wdr.IDBpnpSaveAtmPar)
        self.Bind(wx.EVT_BUTTON, self.OnCalcStats, id=pnpgui_wdr.IDBpnpCalcStats)
        # Other internal values
        # Initiate controls
        self.AtmParTable = AtomicParametersTable(self.GridAtomParam)
        self.GridAtomParam.SetTable(self.AtmParTable, True)

        # init other staff which need both widgets and data
        self.PrepTab_InitDBforQRdiel()
        self.PrepTab_InitDBforSR()

    def PrepTab_InitDBforQRdiel(self):
        # FindWindowById
        # IDCpnpFF4QRdiel
        self.ChoiceFF4QRdiel.Append("None")
        self.ChoiceFF4QRdiel.SetSelection(0)

        if molset is None:
            log.error("Can not load QR DBs because executed without Harlem")
            return

        # Load QRDB
        print("Loading QRdielDB...")
        self.strHarlemSetQR = "Use Harlem Set Q and R"
        self.ChoiceFF4QRdiel.Append(self.strHarlemSetQR)
        self.ChoiceFF4QRdiel.SetSelection(1)

        QRdielDB = molset.GetQRDB()
        if QRdielDB.GetFF("AMBER94") is None:
            from .QR_AMBER import LoadQRDB_AMBER94
            LoadQRDB_AMBER94()

        if QRdielDB.GetFF("PARSE94") is None:
            from .QR_PARSE94 import LoadQRDB_PARSE94
            LoadQRDB_PARSE94()

        QRdielDB.PrintFFsInfo()
        for i in range(QRdielDB.NumFF()):
            ffname = QRdielDB.GetFFbyNum(i).Name
            print(ffname)
            self.ChoiceFF4QRdiel.Append(ffname)
            if ffname == "AMBER94":
                self.ChoiceFF4QRdiel.SetSelection(i + 1)
        print("done with loading QRdielDB.")

    def PrepTab_InitDBforSR(self):
        self.ChoiceFF4SR.Append("None")
        self.ChoiceFF4SR.SetSelection(0)

        if self.pnpmod is None:
            log.error("Can not load SR-MD DBs because executed without Harlem")
            return

        self.ChoiceFF4SR.Append("SR-MD")
        self.ChoiceFF4SR.SetSelection(1)

        from .pnpsgui import pnpgui_mod_dir
        res_db = os.path.join(pnpgui_mod_dir, "db_aar_sr1.pan")
        print(("Loading DB for SR-MD from %s ..." % res_db))
        self.pnpmod.ReadPANDB(res_db, True)
        print("done with Loading DB for SR-MD.")

    def ValidateParameters(self):
        return self._ValidateParameters()

    def GenerateRunScript(self):
        return

    def OnSetParam(self, _):
        if not self.ValidateParameters():
            return

        if molset is None or pnps is None:
            log.error("Can not executed this part without Harlem")
            return

        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print("Setting Parameters for Molecules...")
        # Process

        AtmParFileNameOut = FileBB_GetValue(self.FBB_AtmParFileName_Out)
        # Q and Rdiel
        QRdielDB = molset.GetQRDB()
        selectedFF = str(self.ChoiceFF4QRdiel.GetStringSelection())
        if selectedFF not in ("None", self.strHarlemSetQR):
            msg = "Atoms' charge and radii in Harlem MolSet will be overwritten from %s FF!\n" % selectedFF
            msg += "Continue?"
            log.warning(msg)
            dlg = wx.MessageDialog(self, msg, "Warning", wx.OK | wx.CANCEL)
            if dlg.ShowModal() != wx.ID_OK:
                dlg.Destroy()
                return False
            dlg.Destroy()
            QRdielDB.SetAtomsParam(self.pmset, selectedFF)
        # Init IER
        AddRion1 = float(self.NumCtrl_AddRion1.GetValue())
        AddRion2 = float(self.NumCtrl_AddRion2.GetValue())

        # Soft Repulsion
        if self.ChoiceFF4SR.GetStringSelection() == "SR-MD":
            self.pnpmod.AssignPAN(True)

        natoms = self.pnpmod.mSR_A_K.size()
        print(("Total number of atoms is %s" % natoms))
        # Output
        if self.GridAtomParam.GetNumberRows() > 0:
            self.GridAtomParam.DeleteRows(0, self.GridAtomParam.GetNumberRows())

        # self.GridAtomParam.AppendRows(natoms)

        # import time
        # start = time.time()
        pyacc = molset.PyAccMolSetProp(self.pmset)
        AtomsSerNo = pyacc.GetAtomsSerNoAsVec()
        ResidueSerNo = pyacc.GetResidueSerNoAsVec()
        AtomsCharge = pyacc.GetAtomsChargeAsVec()
        AtomsRadius = pyacc.GetAtomsRadiusAsVec()
        AtomsName = pyacc.GetAtomsNameAsVec()
        ResidueName = pyacc.GetResidueNameAsVec()
        AtomsCoorX = pyacc.GetAtomsCoorXAsVec()
        AtomsCoorY = pyacc.GetAtomsCoorYAsVec()
        AtomsCoorZ = pyacc.GetAtomsCoorZAsVec()

        AtomsIER1 = pyacc.GetAtomsIonExcludedRadiusAsVec(AddRion1)
        AtomsIER2 = pyacc.GetAtomsIonExcludedRadiusAsVec(AddRion2)

        if selectedFF == "None":
            for i in range(AtomsCharge.size()):
                AtomsCharge[i] = 0
                AtomsRadius[i] = 0
                AtomsIER1[i] = AddRion1
                AtomsIER2[i] = AddRion2

        # elapsed = (time.time() - start)
        # print "filling some values %f"%(elapsed)

        # start = time.time()

        self.AtmParTable.SetAtomicParameters(natoms=natoms,
                                             AtomsSerNo=AtomsSerNo,
                                             ResidueSerNo=ResidueSerNo,
                                             AtomsCharge=AtomsCharge,
                                             AtomsRadius=AtomsRadius,
                                             AtomsName=AtomsName,
                                             ResidueName=ResidueName,
                                             AtomsCoorX=AtomsCoorX,
                                             AtomsCoorY=AtomsCoorY,
                                             AtomsCoorZ=AtomsCoorZ,
                                             AtomsIER1=AtomsIER1,
                                             AtomsIER2=AtomsIER2,
                                             SR_A_K=self.pnpmod.mSR_A_K,
                                             SR_N_K=self.pnpmod.mSR_N_K,
                                             SR_A_Cl=self.pnpmod.mSR_A_Cl,
                                             SR_N_Cl=self.pnpmod.mSR_N_Cl)

        self.OnCalcStats(None)

        self.saveAtmParTable(AtmParFileNameOut)

    def saveAtmParTable(self, filename):
        q_tot = sum(self.AtmParTable.AtomsCharge)
        reminder = abs(round(q_tot)-q_tot)
        if reminder >= 0.01:
            msg = "Total charge is significantly non-integer (%s)!\n" % q_tot
            msg += "Save it to file anyway?"
            log.error(msg)
            dlg = wx.MessageDialog(self, msg, "Error In Atoms' Parameters", wx.OK | wx.CANCEL)
            if dlg.ShowModal() != wx.ID_OK:
                dlg.Destroy()
                return False
            dlg.Destroy()

        if abs(q_tot) >= 100.0:
            msg = "Total charge is large (%s)!\n" % q_tot
            msg += "Save it to file anyway?"
            log.error(msg)
            dlg = wx.MessageDialog(self, msg, "Error In Atoms' Parameters", wx.OK | wx.CANCEL)
            if dlg.ShowModal() != wx.ID_OK:
                dlg.Destroy()
                return False
            dlg.Destroy()

        pyacc = molset.PyAccMolSetProp(self.pmset)
        pyacc.WriteAtomParamFileForPNP(filename,
                                       self.AtmParTable.ResidueSerNo, self.AtmParTable.ResidueName,
                                       self.AtmParTable.AtomsSerNo, self.AtmParTable.AtomsName,
                                       self.AtmParTable.AtomsCoorX, self.AtmParTable.AtomsCoorY,
                                       self.AtmParTable.AtomsCoorZ,
                                       self.AtmParTable.AtomsCharge, self.AtmParTable.AtomsRadius,
                                       self.AtmParTable.AtomsIER1, self.AtmParTable.AtomsIER2,
                                       self.AtmParTable.SR_A_K, self.AtmParTable.SR_N_K, self.AtmParTable.SR_A_Cl,
                                       self.AtmParTable.SR_N_Cl)
        print("Done with Setting Parameters for Molecules")
        print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")

    def OnSaveAgainAtomParam(self, _):
        if not self.ValidateParameters():
            return
        # @todo make it work
        # do the trick
        AtmParFileName = str(self.FBB_AtmParFileName_Out.GetValue())
        print(("Saving Parameters for Molecules to %s ..." % AtmParFileName))
        self.OnCalcStats(None)
        self.saveAtmParTable(AtmParFileName)
        print("done")

    def OnCalcStats(self, _):
        if self.AtmParTable.natoms == 0:
            log.error("There is no atoms in ")
        print(("Total charge: %s" % sum(self.AtmParTable.AtomsCharge)))
        print(("x min,max: %s. %s" % (min(self.AtmParTable.AtomsCoorX),max(self.AtmParTable.AtomsCoorX))))
        print(("y min,max: %s. %s" % (min(self.AtmParTable.AtomsCoorY), max(self.AtmParTable.AtomsCoorY))))
        print(("z min,max: %s. %s" % (min(self.AtmParTable.AtomsCoorZ), max(self.AtmParTable.AtomsCoorZ))))



