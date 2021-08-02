import os
import inspect

from math import isnan, isinf

from .atom_params_panel import AtomParamsPanel
from .create_map_panel import CreateMapsPanel
from .pbsr_panel import PBSRPanel

from .utils import FileBB_CheckFileExistence, ShowVectorField3D
from .pnpgui_wdr import *

try:
    import pnps
    pnps_loaded = True
except ImportError:
    pnps = None
    pnps_loaded = False
try:
    from harlempy import molset
    harlem_module = True
except ImportError:
    molset = None
    harlem_module = False

# pointer to wx application
main_frame = None
pnpgui_mod_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))


#######################################################################
# Table for Map Patching
class cDielInUseDataTable(wx.grid.PyGridTableBase):
    def __init__(self):
        wx.grid.PyGridTableBase.__init__(self)

        self.colLabels = []
        self.dataTypes = []
        self.data = [[]]
        for i in range(1, 15):
            # self.colLabels.append("%d"%(i))
            self.colLabels.append("")
            self.dataTypes.append(wx.grid.GRID_VALUE_FLOAT + ':6,2')
            self.data[0].append(0.0)
        self.SetRowLabelValue(0, " ")
        self.data[0][0] = 80.0
        self.data[0][1] = 4.0

        self.PatchMapDielDataTable = None

    # --------------------------------------------------
    # required methods for the wxPyGridTableBase interface

    def GetNumberRows(self):
        return len(self.data)

    def GetNumberCols(self):
        return len(self.data[0])

    def IsEmptyCell(self, row, col):
        try:
            return not self.data[row][col]
        except IndexError:
            return True

    # Get/Set values in the table.  The Python version of these
    # methods can handle any data-type, (as long as the Editor and
    # Renderer understands the type too,) not just strings as in the
    # C++ version.
    def GetValue(self, row, col):
        try:
            return self.data[row][col]
        except IndexError:
            return ''

    def SetValue(self, row, col, value):
        try:
            self.data[row][col] = value
            if self.PatchMapDielDataTable != None:
                self.PatchMapDielDataTable.NewDielConst(self.data[0])
        except IndexError:
            pass

    # --------------------------------------------------
    # Some optional methods

    # Called when the grid needs to display labels
    def GetColLabelValue(self, col):
        return self.colLabels[col]

    # Called to determine the kind of editor/renderer to use by
    # default, doesn't necessarily have to be the same type used
    # natively by the editor/renderer if they know how to convert.
    def GetTypeName(self, row, col):
        return self.dataTypes[col]

    # Called to determine how the data can be fetched and stored by the
    # editor and renderer.  This allows you to enforce some type-safety
    # in the grid.
    def CanGetValueAs(self, row, col, typeName):
        colType = self.dataTypes[col].split(':')[0]
        if typeName == colType:
            return True
        else:
            return False

    def CanSetValueAs(self, row, col, typeName):
        return self.CanGetValueAs(row, col, typeName)

    def SetPatchMapDielDataTable(self, m_PatchMapDielDataTable):
        self.PatchMapDielDataTable = m_PatchMapDielDataTable

    def SetNewDielConstInUse(self, vec):
        for i in range(14):
            self.data[0][i] = 0.0
        m = 14
        if len(vec) < m:
            m = len(vec)
        for i in range(m):
            self.data[0][i] = vec[i]

        self.SetValue(0, 0, vec[0])
        msg = wx.grid.GridTableMessage(self, wx.grid.GRIDTABLE_REQUEST_VIEW_GET_VALUES)
        self.GetView().ProcessTableMessage(msg)


# ---------------------------------------------------------------------------
class cDielInUseTableGrid(wx.grid.Grid):
    def __init__(self, parent):
        wx.grid.Grid.__init__(self, parent, -1, style=0)
        table = cDielInUseDataTable()

        # The second parameter means that the grid is to take ownership of the
        # table and will destroy it when done.  Otherwise you would need to keep
        # a reference to it and call it's Destroy method later.
        self.SetTable(table, True)

        self.SetRowLabelSize(1)
        self.SetColLabelSize(1)
        self.EnableScrolling(False, False)
        self.SetMargins(1, 1)
        self.AutoSizeColumns(False)

        wx.grid.EVT_GRID_CELL_LEFT_DCLICK(self, self.OnLeftDClick)

    # I do this because I don't like the default behaviour of not starting the
    # cell editor on double clicks, but only a second click.
    def OnLeftDClick(self, evt):
        if self.CanEnableCellControl():
            self.EnableCellEditControl()


#######################################################################
# Table for Dielectric Constant Map Patching
class cPatchMapDielDataTable(wx.grid.PyGridTableBase):
    def GetStrRepOfDielConstInUse(self, m_DielConstInUse):
        m_DielConstInUseStrRep = []
        for i in range(len(m_DielConstInUse)):
            m_DielConstInUseStrRep.append("DielConst[%d]=%.2f" % (i + 1, m_DielConstInUse[i]))
        return m_DielConstInUseStrRep

    def GetStrRepOfJustDielConstInUse(self):
        m_DielConstInUseStrRep = []
        for i in range(len(self.DielConstInUse)):
            m_DielConstInUseStrRep.append("%.2f" % (self.DielConstInUse[i]))
        return m_DielConstInUseStrRep

    def GetDataTypeOfDielConstInUse(self, m_DielConstInUse):
        m_DielConstInUseStrRep = self.GetStrRepOfDielConstInUse(m_DielConstInUse)
        if len(m_DielConstInUseStrRep) >= 1:
            m_DielConstInUseDataType = wx.grid.GRID_VALUE_CHOICE + ':' + m_DielConstInUseStrRep[0]

            for i in range(1, len(m_DielConstInUseStrRep)):
                m_DielConstInUseDataType = m_DielConstInUseDataType + ',' + m_DielConstInUseStrRep[i]
            return m_DielConstInUseDataType
        else:
            return wx.grid.GRID_VALUE_CHOICE + ':empty'

    def GetJustDielStrRepr(self, FullStrValue):
        if self.DielConstInUseStrRep.count(FullStrValue) > 0:
            return "%.2f" % (self.DielConstInUse[self.DielConstInUseStrRep.index(FullStrValue)])
        else:
            return "None"

    def DielFloatStrRepr(self, FloatValue):
        return "%.2f" % (FloatValue)

    def GetDefaultElement(self):
        return [self.DielConstInUseStrRep[0], self.DielConstInUseStrRep[0], 0.0, 0.0, 0.0, 0.0, 0.0]

    def __init__(self):
        wx.grid.PyGridTableBase.__init__(self)
        # self.log = log
        self.DielConstInUse = [80.0, 40.0, 20.0]
        self.DielConstInUseStrRep = self.GetStrRepOfDielConstInUse(self.DielConstInUse)
        self.DielConstInUseDataType = self.GetDataTypeOfDielConstInUse(self.DielConstInUse)

        self.colLabels = ['Diel.Const. To Overwrite', 'New Diel.Const.', 'x', 'y', 'z1', 'z2', 'R']
        self.dataTypes = [self.DielConstInUseDataType,
                          self.DielConstInUseDataType,
                          wx.grid.GRID_VALUE_FLOAT + ':10,3',
                          wx.grid.GRID_VALUE_FLOAT + ':10,3',
                          wx.grid.GRID_VALUE_FLOAT + ':10,3',
                          wx.grid.GRID_VALUE_FLOAT + ':10,3',
                          wx.grid.GRID_VALUE_FLOAT + ':10,3']
        self.data = [
            self.GetDefaultElement()
        ]

    # --------------------------------------------------
    # required methods for the wxPyGridTableBase interface

    def GetNumberRows(self):
        return len(self.data)

    def GetNumberCols(self):
        return len(self.data[0])

    def IsEmptyCell(self, row, col):
        try:
            return not self.data[row][col]
        except IndexError:
            return True

    # Get/Set values in the table.  The Python version of these
    # methods can handle any data-type, (as long as the Editor and
    # Renderer understands the type too,) not just strings as in the
    # C++ version.
    def GetValue(self, row, col):
        try:
            return self.data[row][col]
        except IndexError:
            return ''

    def SetValue(self, row, col, value):
        def innerSetValue(row, col, value):
            try:
                self.data[row][col] = value
            except IndexError:
                pass

        innerSetValue(row, col, value)

        # --------------------------------------------------

    # Some optional methods

    # Called when the grid needs to display labels
    def GetColLabelValue(self, col):
        return self.colLabels[col]

    # Called to determine the kind of editor/renderer to use by
    # default, doesn't necessarily have to be the same type used
    # natively by the editor/renderer if they know how to convert.
    def GetTypeName(self, row, col):
        return self.dataTypes[col]

    # Called to determine how the data can be fetched and stored by the
    # editor and renderer.  This allows you to enforce some type-safety
    # in the grid.
    def CanGetValueAs(self, row, col, typeName):
        colType = self.dataTypes[col].split(':')[0]
        if typeName == colType:
            return True
        else:
            return False

    def CanSetValueAs(self, row, col, typeName):
        return self.CanGetValueAs(row, col, typeName)

    def NewDielConst(self, m_newDielConstInUse):
        # find last zero
        LastNotZero = 0
        for i in range(len(m_newDielConstInUse)):
            if m_newDielConstInUse[i] != 0.0:
                LastNotZero = i
        LastNotZero = LastNotZero + 1
        newDielConstInUse = m_newDielConstInUse[0:LastNotZero]
        newDielConstInUseStrRep = self.GetStrRepOfDielConstInUse(newDielConstInUse)
        newDielConstInUseDataType = self.GetDataTypeOfDielConstInUse(newDielConstInUse)

        # change values
        oldDielConstInUseStrRep = self.DielConstInUseStrRep
        for i in range(len(self.data)):
            if oldDielConstInUseStrRep.index(self.data[i][0]) < len(newDielConstInUseStrRep):
                self.data[i][0] = newDielConstInUseStrRep[oldDielConstInUseStrRep.index(self.data[i][0])]
            else:
                self.data[i][0] = newDielConstInUseStrRep[0]
            if oldDielConstInUseStrRep.index(self.data[i][1]) < len(newDielConstInUseStrRep):
                self.data[i][1] = newDielConstInUseStrRep[oldDielConstInUseStrRep.index(self.data[i][1])]
            else:
                self.data[i][1] = newDielConstInUseStrRep[0]

        self.DielConstInUse = newDielConstInUse
        self.DielConstInUseStrRep = newDielConstInUseStrRep
        self.DielConstInUseDataType = newDielConstInUseDataType
        self.dataTypes[0] = self.DielConstInUseDataType
        self.dataTypes[1] = self.DielConstInUseDataType

        # self.GetView().Update()
        msg = wx.grid.GridTableMessage(self, wx.grid.GRIDTABLE_REQUEST_VIEW_GET_VALUES)
        self.GetView().ProcessTableMessage(msg)

    def NewElement(self):
        # add a new row
        self.data.append(self.GetDefaultElement())

        # tell the grid we've added a row
        msg = wx.grid.GridTableMessage(self,  # The table
                                       wx.grid.GRIDTABLE_NOTIFY_ROWS_APPENDED,  # what we did to it
                                       1  # how many
                                       )

        self.GetView().ProcessTableMessage(msg)

    def NewElementWithValue(self, Value):
        # add a new row
        self.data.append(Value)

        # tell the grid we've added a row
        msg = wx.grid.GridTableMessage(self,  # The table
                                       wx.grid.GRIDTABLE_NOTIFY_ROWS_APPENDED,  # what we did to it
                                       1  # how many
                                       )

        self.GetView().ProcessTableMessage(msg)

    # def DeleteRows(self, pos, numRows):
    #    self.data.pop(pos)
    #    return True
    def DelElementByRow(self, row):
        # print "DelElement",self.GetView().GetSelectedRows()
        # self.GetView().DeleteRows(self.GetView().GetSelectedRows()[0],1,True)
        self.data.pop(row)
        # tell the grid we've delete a row
        msg = wx.grid.GridTableMessage(self,  # The table
                                       wx.grid.GRIDTABLE_NOTIFY_ROWS_DELETED,  # what we did to it
                                       len(self.data),  # how many
                                       1)

        self.GetView().ProcessTableMessage(msg)

    def DelElement(self):
        if len(self.GetView().GetSelectedRows()) > 0:
            # print "DelElement",self.GetView().GetSelectedRows()
            # self.GetView().DeleteRows(self.GetView().GetSelectedRows()[0],1,True)
            self.data.pop(self.GetView().GetSelectedRows()[0])
            # tell the grid we've delete a row
            msg = wx.grid.GridTableMessage(self,  # The table
                                           wx.grid.GRIDTABLE_NOTIFY_ROWS_DELETED,  # what we did to it
                                           len(self.data),  # how many
                                           1)

            self.GetView().ProcessTableMessage(msg)


# ---------------------------------------------------------------------------
class cPatchMapDielTableGrid(wx.grid.Grid):
    def __init__(self, parent):
        wx.grid.Grid.__init__(self, parent, -1)
        table = cPatchMapDielDataTable()

        # The second parameter means that the grid is to take ownership of the
        # table and will destroy it when done.  Otherwise you would need to keep
        # a reference to it and call it's Destroy method later.
        self.SetTable(table, True)

        self.SetColLabelSize(18)
        self.SetRowLabelSize(36)
        self.SetMargins(0, 0)
        self.AutoSizeColumns(False)

        wx.grid.EVT_GRID_CELL_LEFT_DCLICK(self, self.OnLeftDClick)

    # I do this because I don't like the default behaviour of not starting the
    # cell editor on double clicks, but only a second click.
    def OnLeftDClick(self, evt):
        if self.CanEnableCellControl():
            self.EnableCellEditControl()


#######################################################################
# Table for Diffusion Coefficient Map Patching
class cPatchMapDiffDataTable(wx.grid.PyGridTableBase):
    def GetDefaultElement(self):
        return [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    def __init__(self):
        wx.grid.PyGridTableBase.__init__(self)

        self.colLabels = ['Dscale1(K+)', 'Dscale2(K+)', 'Dscale1(Cl-)', 'Dscale2(Cl-)', 'x', 'y', 'z1', 'z2', 'R']
        self.dataTypes = [wx.grid.GRID_VALUE_FLOAT + ':10,3',
                          wx.grid.GRID_VALUE_FLOAT + ':10,3',
                          wx.grid.GRID_VALUE_FLOAT + ':10,3',
                          wx.grid.GRID_VALUE_FLOAT + ':10,3',
                          wx.grid.GRID_VALUE_FLOAT + ':10,3',
                          wx.grid.GRID_VALUE_FLOAT + ':10,3',
                          wx.grid.GRID_VALUE_FLOAT + ':10,3',
                          wx.grid.GRID_VALUE_FLOAT + ':10,3',
                          wx.grid.GRID_VALUE_FLOAT + ':10,3']
        self.data = [
            self.GetDefaultElement()
        ]

    # --------------------------------------------------
    # required methods for the wxPyGridTableBase interface

    def GetNumberRows(self):
        return len(self.data)

    def GetNumberCols(self):
        return len(self.data[0])

    def IsEmptyCell(self, row, col):
        try:
            return not self.data[row][col]
        except IndexError:
            return True

    # Get/Set values in the table.  The Python version of these
    # methods can handle any data-type, (as long as the Editor and
    # Renderer understands the type too,) not just strings as in the
    # C++ version.
    def GetValue(self, row, col):
        try:
            return self.data[row][col]
        except IndexError:
            return ''

    def SetValue(self, row, col, value):
        def innerSetValue(row, col, value):
            try:
                self.data[row][col] = value
            except IndexError:
                pass

        innerSetValue(row, col, value)

        # --------------------------------------------------

    # Some optional methods

    # Called when the grid needs to display labels
    def GetColLabelValue(self, col):
        return self.colLabels[col]

    # Called to determine the kind of editor/renderer to use by
    # default, doesn't necessarily have to be the same type used
    # natively by the editor/renderer if they know how to convert.
    def GetTypeName(self, row, col):
        return self.dataTypes[col]

    # Called to determine how the data can be fetched and stored by the
    # editor and renderer.  This allows you to enforce some type-safety
    # in the grid.
    def CanGetValueAs(self, row, col, typeName):
        colType = self.dataTypes[col].split(':')[0]
        if typeName == colType:
            return True
        else:
            return False

    def CanSetValueAs(self, row, col, typeName):
        return self.CanGetValueAs(row, col, typeName)

    def NewElement(self):
        # add a new row
        self.data.append(self.GetDefaultElement())

        # tell the grid we've added a row
        msg = wx.grid.GridTableMessage(self,  # The table
                                       wx.grid.GRIDTABLE_NOTIFY_ROWS_APPENDED,  # what we did to it
                                       1  # how many
                                       )

        self.GetView().ProcessTableMessage(msg)

    def NewElementWithValue(self, Value):
        # add a new row
        self.data.append(Value)

        # tell the grid we've added a row
        msg = wx.grid.GridTableMessage(self,  # The table
                                       wx.grid.GRIDTABLE_NOTIFY_ROWS_APPENDED,  # what we did to it
                                       1  # how many
                                       )

        self.GetView().ProcessTableMessage(msg)

    # def DeleteRows(self, pos, numRows):
    #    self.data.pop(pos)
    #    return True
    def DelElementByRow(self, row):
        # print "DelElement",self.GetView().GetSelectedRows()
        # self.GetView().DeleteRows(self.GetView().GetSelectedRows()[0],1,True)
        self.data.pop(row)
        # tell the grid we've delete a row
        msg = wx.grid.GridTableMessage(self,  # The table
                                       wx.grid.GRIDTABLE_NOTIFY_ROWS_DELETED,  # what we did to it
                                       len(self.data),  # how many
                                       1)

        self.GetView().ProcessTableMessage(msg)

    def DelElement(self):
        if len(self.GetView().GetSelectedRows()) > 0:
            self.data.pop(self.GetView().GetSelectedRows()[0])
            # tell the grid we've delete a row
            msg = wx.grid.GridTableMessage(self,  # The table
                                           wx.grid.GRIDTABLE_NOTIFY_ROWS_DELETED,  # what we did to it
                                           len(self.data),  # how many
                                           1)

            self.GetView().ProcessTableMessage(msg)


# ---------------------------------------------------------------------------
class cPatchMapDiffTableGrid(wx.grid.Grid):
    def __init__(self, parent):
        wx.grid.Grid.__init__(self, parent, -1)
        table = cPatchMapDiffDataTable()

        # The second parameter means that the grid is to take ownership of the
        # table and will destroy it when done.  Otherwise you would need to keep
        # a reference to it and call it's Destroy method later.
        self.SetTable(table, True)

        self.SetColLabelSize(18)
        self.SetRowLabelSize(36)
        self.SetMargins(0, 0)
        self.AutoSizeColumns(False)

        wx.grid.EVT_GRID_CELL_LEFT_DCLICK(self, self.OnLeftDClick)

    # I do this because I don't like the default behaviour of not starting the
    # cell editor on double clicks, but only a second click.
    def OnLeftDClick(self, evt):
        if self.CanEnableCellControl():
            self.EnableCellEditControl()


#######################################################################
# PNPFrame
class PNPFrame(wx.Frame):
    """
    Dialog for controlling PNP-SR calculations

    Each tab corresponds for a destinctive step of PNP-SR calculation, abreviated as {Step Name} within class:
        - Preporation of structure and atomic parameters
        - Map creation
        - Patches
        PBSR - Poisson Boltzmann with Soft Repulsion
        IAVRef - IAV Refinement
        - PNP-SR itself

    On Naming convention and similar withing this class
    Names of wx controls identification: ID{control type}pnp{Step Name}_{Button name}
        {control type}:
            B - button

        example:IDBpnpPBSR_Run
    """

    def __init__(self, pmset, parent="not set"):
        global main_frame
        main_frame = parent if parent != "not set" else wx.FindWindowByName("HARLEM")
        wx.Frame.__init__(self,
                          main_frame,
                          title="PNPGUI - %s" % ("None" if pmset is None else pmset.GetName()),
                          size=[640, 600])
        if main_frame is None:
            main_frame = self
        # init stuff
        self.pmset = pmset
        self.pnpmod = None if pmset is None else pmset.GetPNPMod(True)

        # init widgets
        set_sizer = True
        call_fit = True
        item0 = wx.GridSizer(0, 1, 0, 0)

        self.NotebookMain = wx.Notebook(self, IDpnpNOTEBOOK, wx.DefaultPosition, [640, 600], 0)

        self.atom_params_panel = AtomParamsPanel(self.NotebookMain)
        self.NotebookMain.AddPage(self.atom_params_panel, "Atomic Parameters")

        self.create_maps_panel = CreateMapsPanel(self.NotebookMain)
        self.NotebookMain.AddPage(self.create_maps_panel, "Maps Building")

        item5 = wx.Panel(self.NotebookMain, -1)
        pnpMapPatch(item5, False)
        self.NotebookMain.AddPage(item5, "Maps Patching")

        self.pbsr_panel = PBSRPanel(self.NotebookMain)
        self.NotebookMain.AddPage(self.pbsr_panel, "PB-SR")

        item7 = wx.Panel(self.NotebookMain, -1)
        pnpIAVRef(item7, False)
        self.NotebookMain.AddPage(item7, "IAV Refinement")

        item8 = wx.Panel(self.NotebookMain, -1)
        pnpPNPSR(item8, False)
        self.NotebookMain.AddPage(item8, "PNP-SR")

        item0.Add(self.NotebookMain, 0, wx.GROW | wx.ALL, 5)

        if set_sizer:
            self.SetSizer(item0)
            if call_fit:
                item0.SetSizeHints(self)

        #######################################################################
        # Initialize WX controls
        #   Preporation Tab
        #self.PrepTab_InitWX()
        #   Map Creation Interface
        # self.CreateMaps_InitWX()
        #   Map Patching Interface
        self.PatchMaps_InitWX()
        #   PB-SR Interface
        #self.PBSR_InitWX()
        #   IAVRefinement Interface
        self.IAVRef_InitWX()
        #   PNP-SR interface
        self.PNPSR_InitWX()
        #######################################################################

        # finally show dialog
        # self.Layout()
        self.Show(True)

    ###############################################################################################
    # Map Patching
    def PatchMaps_InitWX(self):

        NotebookDielOrDiffPatch = self.FindWindowById(IDNpnpPM_DielOrDiffPatch)
        # Patching Dielectric Constant Maps
        DielPatchPanel = wx.Panel(NotebookDielOrDiffPatch, -1)

        item0 = wx.BoxSizer(wx.VERTICAL)

        item1 = wx.StaticBox(DielPatchPanel, -1, "Dielectric Constant Values")
        item2 = wx.StaticBoxSizer(item1, wx.VERTICAL)
        #
        item3 = wx.StaticText(DielPatchPanel, -1,
                              "PNP module can handle only 14 descrite dielectric constant values. Following table shows values to use:",
                              wx.DefaultPosition, wx.DefaultSize, 0)
        item2.Add(item3, 0, wx.ALIGN_LEFT, 5)
        # Add Table for Diel InUse
        self.Table_PatchMap_DielInUse = cDielInUseTableGrid(DielPatchPanel)
        item2.Add(self.Table_PatchMap_DielInUse, 0, wx.GROW | wx.ALIGN_CENTER_VERTICAL, 5)
        item0.Add(item2, 0, wx.GROW | wx.ALIGN_CENTER_VERTICAL, 5)
        #
        item5 = wx.StaticBox(DielPatchPanel, -1, "Patching Dielectric Constant Maps")
        item6 = wx.StaticBoxSizer(item5, wx.VERTICAL)

        item7 = wx.BoxSizer(wx.HORIZONTAL)

        item7a = wx.Button(DielPatchPanel, -1, "New Entry", wx.DefaultPosition, wx.DefaultSize, 0)
        item7a.SetToolTip(wx.ToolTip("Append new entry to  the table"))
        item7.Add(item7a, 0, wx.ALIGN_CENTER, 5)

        item7b = wx.Button(DielPatchPanel, -1, "Delete Entry", wx.DefaultPosition, wx.DefaultSize, 0)
        item7b.SetToolTip(wx.ToolTip("Delete selected row from the table"))
        item7.Add(item7b, 0, wx.ALIGN_CENTER, 5)

        item7c = wx.Button(DielPatchPanel, -1, "Save Table", wx.DefaultPosition, wx.DefaultSize, 0)
        item7c.SetToolTip(wx.ToolTip("Save table as file"))
        item7.Add(item7c, 0, wx.ALIGN_CENTER, 5)

        item7d = wx.Button(DielPatchPanel, -1, "Load Table", wx.DefaultPosition, wx.DefaultSize, 0)
        item7d.SetToolTip(wx.ToolTip("Load table from file"))
        item7.Add(item7d, 0, wx.ALIGN_CENTER, 5)

        self.Bind(wx.EVT_BUTTON, self.PatchMaps_DielNewElement, item7a)
        self.Bind(wx.EVT_BUTTON, self.PatchMaps_DielDelElement, item7b)
        self.Bind(wx.EVT_BUTTON, self.PatchMaps_DielSaveTable, item7c)
        self.Bind(wx.EVT_BUTTON, self.PatchMaps_DielLoadTable, item7d)

        item6.Add(item7, 0, wx.GROW | wx.ALIGN_CENTER_VERTICAL, 5)
        # Add Table for Diel Patching
        self.Table_PatchMap_Diel = cPatchMapDielTableGrid(DielPatchPanel)
        self.Table_PatchMap_DielInUse.GetTable().SetPatchMapDielDataTable(self.Table_PatchMap_Diel.GetTable())
        self.Table_PatchMap_DielInUse.GetTable().SetValue(0, 3, 0.0)

        item6.Add(self.Table_PatchMap_Diel, 1, wx.GROW | wx.ALIGN_CENTER_VERTICAL, 5)

        item0.Add(item6, 1, wx.GROW | wx.ALIGN_CENTER_VERTICAL, 5)
        DielPatchPanel.SetSizer(item0)
        NotebookDielOrDiffPatch.AddPage(DielPatchPanel, "Patching Dielectric Constant Maps")
        #
        # Patching DiffusionCoefficient Maps
        DiffPatchPanel = wx.Panel(NotebookDielOrDiffPatch, -1)

        itemD0 = wx.BoxSizer(wx.VERTICAL)
        itemD5 = wx.StaticBox(DiffPatchPanel, -1, "Patching Diffusion Coefficient Maps")
        itemD6 = wx.StaticBoxSizer(itemD5, wx.VERTICAL)

        itemD7 = wx.BoxSizer(wx.HORIZONTAL)

        itemD7a = wx.Button(DiffPatchPanel, -1, "New Entry", wx.DefaultPosition, wx.DefaultSize, 0)
        itemD7a.SetToolTip(wx.ToolTip("Append new entry to  the table"))
        itemD7.Add(itemD7a, 0, wx.ALIGN_CENTER, 5)

        itemD7b = wx.Button(DiffPatchPanel, -1, "Delete Entry", wx.DefaultPosition, wx.DefaultSize, 0)
        itemD7b.SetToolTip(wx.ToolTip("Delete selected row from the table"))
        itemD7.Add(itemD7b, 0, wx.ALIGN_CENTER, 5)

        itemD7c = wx.Button(DiffPatchPanel, -1, "Save Table", wx.DefaultPosition, wx.DefaultSize, 0)
        itemD7c.SetToolTip(wx.ToolTip("Save table as file"))
        itemD7.Add(itemD7c, 0, wx.ALIGN_CENTER, 5)

        itemD7d = wx.Button(DiffPatchPanel, -1, "Load Table", wx.DefaultPosition, wx.DefaultSize, 0)
        itemD7d.SetToolTip(wx.ToolTip("Load table from file"))
        itemD7.Add(itemD7d, 0, wx.ALIGN_CENTER, 5)

        self.Bind(wx.EVT_BUTTON, self.PatchMaps_DiffNewElement, itemD7a)
        self.Bind(wx.EVT_BUTTON, self.PatchMaps_DiffDelElement, itemD7b)
        self.Bind(wx.EVT_BUTTON, self.PatchMaps_DiffSaveTable, itemD7c)
        self.Bind(wx.EVT_BUTTON, self.PatchMaps_DiffLoadTable, itemD7d)

        itemD6.Add(itemD7, 0, wx.GROW | wx.ALIGN_CENTER_VERTICAL, 5)
        # Add Table for Diff Patching
        self.Table_PatchMap_Diff = cPatchMapDiffTableGrid(DiffPatchPanel)
        itemD6.Add(self.Table_PatchMap_Diff, 1, wx.GROW | wx.ALIGN_CENTER_VERTICAL, 5)

        itemD0.Add(itemD6, 1, wx.GROW | wx.ALIGN_CENTER_VERTICAL, 5)
        DiffPatchPanel.SetSizer(itemD0)
        NotebookDielOrDiffPatch.AddPage(DiffPatchPanel, "Patching Diffusion Coefficient Maps")
        # Parameters
        # self.NumCtrl_PatchMap_DielValInUse=self.FindWindowById(IDTCpnpDielValInUse)
        # Input Maps
        self.FileBB_PatchMap_SysTopFileNameIN = self.FindWindowById(IDTCpnpPM_SysTopFileNameIN)
        self.FileBB_PatchMap_DiffFileNameIN = self.FindWindowById(IDTCpnpPM_DiffFileNameIN)
        # Output Maps
        self.FileBB_PatchMap_SysTopFileNameOUT = self.FindWindowById(IDTCpnpPM_SysTopFileNameOUT)
        self.FileBB_PatchMap_DiffFileNameOUT = self.FindWindowById(IDTCpnpPM_DiffFileNameOUT)
        # To validate user entry
        # self.PatchMap_NumCtrls4Validation=[
        #    self.NumCtrl_PatchMap_DielValInUse]
        # bind events
        self.Bind(wx.EVT_BUTTON, self.PatchMaps_LoadMaps, id=IDBpnpPM_LoadMaps)
        self.Bind(wx.EVT_BUTTON, self.PatchMaps_PatchMaps, id=IDBpnpPatchMaps)

        self.Bind(wx.EVT_BUTTON, self.PatchMaps_OnViewStaticCharge, id=IDBpnpPM_ViewStaticCharge)
        self.Bind(wx.EVT_BUTTON, self.PatchMaps_OnViewDiel, id=IDBpnpPM_ViewDiel)
        self.Bind(wx.EVT_BUTTON, self.PatchMaps_OnViewDiff, id=IDBpnpPM_ViewDiff)
        self.Bind(wx.EVT_BUTTON, self.PatchMaps_OnViewPlainDiff, id=IDBpnpPM_ViewPlainDiff)

    def PatchMaps_DielNewElement(self, event):
        self.Table_PatchMap_Diel.GetTable().NewElement()

    def PatchMaps_DielDelElement(self, event):
        self.Table_PatchMap_Diel.GetTable().DelElement()

    def PatchMaps_DielSaveTable(self, event):
        dlg = wx.FileDialog(
            self, message="Save table as ...", defaultDir=os.getcwd(),
            defaultFile="S1_DielPatchingTable.dat", wildcard="DAT files (*.dat)|*.dat|all(*)|*",
            style=wx.SAVE | wx.OVERWRITE_PROMPT
        )
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()

            DielPatchTable = self.Table_PatchMap_Diel.GetTable()
            Table_PatchMap_DielInUse = self.Table_PatchMap_DielInUse.GetTable()

            fOut = open(path, 'wt')
            fOut.write("#Used Dielectric Constants\n")
            for eps in DielPatchTable.GetStrRepOfJustDielConstInUse():
                fOut.write("%s " % (eps))
            fOut.write("\n")
            fOut.write("#Dielectric Constants Patching Table\n")
            fOut.write("#Diel.Const. To Overwrite;\tNew Diel.Const.;\tx;\ty;\tz1;\tz2;\tR;\n")
            for i in range(DielPatchTable.GetNumberRows()):
                epsToOver = DielPatchTable.GetJustDielStrRepr(DielPatchTable.GetValue(i, 0))
                epsNew = DielPatchTable.GetJustDielStrRepr(DielPatchTable.GetValue(i, 1))
                fOut.write("%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (epsToOver, epsNew,
                                                                       DielPatchTable.GetValue(i, 2),
                                                                       DielPatchTable.GetValue(i, 3),
                                                                       DielPatchTable.GetValue(i, 4),
                                                                       DielPatchTable.GetValue(i, 5),
                                                                       DielPatchTable.GetValue(i, 6)))

            fOut.close()

        dlg.Destroy()

    def PatchMaps_DielLoadTable(self, event):
        dlg = wx.FileDialog(
            self, message="Choose a file with table",
            defaultDir=os.getcwd(),
            defaultFile="S1_DielPatchingTable.dat",
            wildcard="DAT files (*.dat)|*.dat|all(*)|*",
            style=wx.OPEN | wx.CHANGE_DIR
        )

        # Show the dialog and retrieve the user response. If it is the OK response,
        # process the data.
        if dlg.ShowModal() == wx.ID_OK:
            # This returns a Python list of files that were selected.
            path = dlg.GetPath()
            DielPatchTable = self.Table_PatchMap_Diel.GetTable()
            Table_PatchMap_DielInUse = self.Table_PatchMap_DielInUse.GetTable()
            # delete all curent elements
            while DielPatchTable.GetNumberRows() > 0:
                DielPatchTable.DelElementByRow(0)
                print(DielPatchTable.GetNumberRows())
            fIn = open(path, 'rt')
            rawtab = fIn.readlines()
            fIn.close()
            tabstr = []
            for i in range(len(rawtab)):
                if rawtab[i][0] != '#':
                    tabstr.append(rawtab[i])

            dielStr = tabstr[0].split()
            for i in range(14):
                Table_PatchMap_DielInUse.SetValue(0, i, 0.0)
            dielStrFloat = []
            for i in range(len(dielStr)):
                Table_PatchMap_DielInUse.SetValue(0, i, float(dielStr[i]))
                dielStrFloat.append(DielPatchTable.DielFloatStrRepr(float(dielStr[i])))
            msg = wx.grid.GridTableMessage(Table_PatchMap_DielInUse, wx.grid.GRIDTABLE_REQUEST_VIEW_GET_VALUES)
            Table_PatchMap_DielInUse.GetView().ProcessTableMessage(msg)

            for i in range(1, len(tabstr)):
                (epsToOver, epsNew, x, y, z1, z2, R) = tabstr[i].split()
                epsToOver = DielPatchTable.DielConstInUseStrRep[
                    dielStrFloat.index(DielPatchTable.DielFloatStrRepr(float(epsToOver)))]
                epsNew = DielPatchTable.DielConstInUseStrRep[
                    dielStrFloat.index(DielPatchTable.DielFloatStrRepr(float(epsNew)))]
                val = [epsToOver, epsNew, float(x), float(y), float(z1), float(z2), float(R)]
                DielPatchTable.NewElementWithValue(val)

        dlg.Destroy()

    def PatchMaps_DiffNewElement(self, event):
        self.Table_PatchMap_Diff.GetTable().NewElement()

    def PatchMaps_DiffDelElement(self, event):
        self.Table_PatchMap_Diff.GetTable().DelElement()

    def PatchMaps_DiffSaveTable(self, event):
        dlg = wx.FileDialog(
            self, message="Save table as ...", defaultDir=os.getcwd(),
            defaultFile="S1_DiffPatchingTable.dat", wildcard="DAT files (*.dat)|*.dat|all(*)|*",
            style=wx.SAVE | wx.OVERWRITE_PROMPT
        )
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()

            DiffPatchTable = self.Table_PatchMap_Diff.GetTable()

            fOut = open(path, 'wt')
            fOut.write("#Dscale1(K+)\tDscale2(K+)\tDscale1(Cl-)\tDscale2(Cl-)\tx\ty\tz1\tz2\tR\n")
            for i in range(DiffPatchTable.GetNumberRows()):
                fOut.write("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (
                    DiffPatchTable.GetValue(i, 0), DiffPatchTable.GetValue(i, 1),
                    DiffPatchTable.GetValue(i, 2), DiffPatchTable.GetValue(i, 3),
                    DiffPatchTable.GetValue(i, 4), DiffPatchTable.GetValue(i, 5),
                    DiffPatchTable.GetValue(i, 6), DiffPatchTable.GetValue(i, 7),
                    DiffPatchTable.GetValue(i, 8)))

            fOut.close()

        dlg.Destroy()

    def PatchMaps_DiffLoadTable(self, event):
        dlg = wx.FileDialog(
            self, message="Choose a file with table",
            defaultDir=os.getcwd(),
            defaultFile="S1_DiffPatchingTable.dat",
            wildcard="DAT files (*.dat)|*.dat|all(*)|*",
            style=wx.OPEN | wx.CHANGE_DIR
        )

        # Show the dialog and retrieve the user response. If it is the OK response,
        # process the data.
        if dlg.ShowModal() == wx.ID_OK:
            # This returns a Python list of files that were selected.
            path = dlg.GetPath()
            DiffPatchTable = self.Table_PatchMap_Diff.GetTable()
            # delete all curent elements
            while DiffPatchTable.GetNumberRows() > 0:
                DiffPatchTable.DelElementByRow(0)
                print(DiffPatchTable.GetNumberRows())
            fIn = open(path, 'rt')
            rawtab = fIn.readlines()
            fIn.close()
            tabstr = []
            for i in range(len(rawtab)):
                if rawtab[i][0] != '#':
                    tabstr.append(rawtab[i])

            for i in range(0, len(tabstr)):
                (Dscale1K, Dscale2K, Dscale1Cl, Dscale2Cl, x, y, z1, z2, R) = tabstr[i].split()
                val = [float(Dscale1K), float(Dscale2K), float(Dscale1Cl), float(Dscale2Cl), float(x), float(y),
                       float(z1), float(z2), float(R)]
                DiffPatchTable.NewElementWithValue(val)

        dlg.Destroy()

    def PatchMaps_LoadMaps(self, event):
        print("Loading Map Information...")
        # check the existance of input map file
        SomeOfInputFilesDoNotExists = False
        SomeOfInputFilesDoNotExistsMessage = ""
        filebrowser = self.FileBB_PatchMap_SysTopFileNameIN
        if (not FileBB_CheckFileExistence(filebrowser)):
            SomeOfInputFilesDoNotExists = True
            SomeOfInputFilesDoNotExistsMessage = SomeOfInputFilesDoNotExistsMessage + "file %s do not exist.\n" % (
                filebrowser.GetValue())
        filebrowser = self.FileBB_PatchMap_DiffFileNameIN
        if (not FileBB_CheckFileExistence(filebrowser)) and (filebrowser.GetValue() != "None"):
            SomeOfInputFilesDoNotExists = True
            SomeOfInputFilesDoNotExistsMessage = SomeOfInputFilesDoNotExistsMessage + "file %s do not exist.\n" % (
                filebrowser.GetValue())
        if SomeOfInputFilesDoNotExists:
            SomeOfInputFilesDoNotExistsMessage = "Some of input files do not exists:" + SomeOfInputFilesDoNotExistsMessage
            dlg = wx.MessageDialog(self, SomeOfInputFilesDoNotExistsMessage, "Error In Input Values", wx.OK)
            dlg.ShowModal()  # Show it
            dlg.Destroy()  # finally destroy it when finished.
            return False
        # Read maps values
        NI = pnps.NodeIndexing()
        NI.ReadFromFile(self.FileBB_PatchMap_SysTopFileNameIN.GetValue())
        DielConstInUse = [0.0] * 14
        for i in range(14):
            DielConstInUse[i] = NI.GetDielConstInUse(i)
        del NI
        print(DielConstInUse)

        self.Table_PatchMap_DielInUse.GetTable().SetNewDielConstInUse(DielConstInUse)
        print("done with Loading Map Information")

    def PatchMaps_PatchMaps(self, event):
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print("Patching Maps...")
        # Check input values
        # check the existance of input map file
        SomeOfInputFilesDoNotExists = False
        SomeOfInputFilesDoNotExistsMessage = ""
        # "None"
        # for filebrowser in [self.FileBB_PatchMap_SysTopFileNameIN,self.FileBB_PatchMap_DiffFileNameIN]:
        filebrowser = self.FileBB_PatchMap_SysTopFileNameIN
        if not FileBB_CheckFileExistence(filebrowser):
            SomeOfInputFilesDoNotExists = True
            SomeOfInputFilesDoNotExistsMessage = SomeOfInputFilesDoNotExistsMessage + "file %s do not exist.\n" % (
                filebrowser.GetValue())
        filebrowser = self.FileBB_PatchMap_DiffFileNameIN
        if filebrowser.GetValue() != "None":
            if not FileBB_CheckFileExistence(filebrowser):
                SomeOfInputFilesDoNotExists = True
                SomeOfInputFilesDoNotExistsMessage = SomeOfInputFilesDoNotExistsMessage + "file %s do not exist.\n" % (
                    filebrowser.GetValue())
        if SomeOfInputFilesDoNotExists:
            SomeOfInputFilesDoNotExistsMessage = "Some of input files do not exists:" + SomeOfInputFilesDoNotExistsMessage
            dlg = wx.MessageDialog(self, SomeOfInputFilesDoNotExistsMessage, "Error In Input Values", wx.OK)
            dlg.ShowModal()  # Show it
            dlg.Destroy()  # finally destroy it when finished.
            return False
        # check the existance of these file and ask user for permition to rewrite
        SomeOfOutputFilesExists = False
        SomeOfOutputFilesExistsMessage = ""
        for filebrowser in [self.FileBB_PatchMap_SysTopFileNameOUT, self.FileBB_PatchMap_DiffFileNameOUT]:
            if filebrowser.GetValue() != "None":
                if FileBB_CheckFileExistence(filebrowser):
                    SomeOfOutputFilesExists = True
                    SomeOfOutputFilesExistsMessage = SomeOfOutputFilesExistsMessage + "file %s exist.\n" % (
                        filebrowser.GetValue())
        if SomeOfOutputFilesExists:
            SomeOfOutputFilesExistsMessage = SomeOfOutputFilesExistsMessage + "\nIf execution will be continued this(these) file(s) will be deleted."
            SomeOfOutputFilesExistsMessage = SomeOfOutputFilesExistsMessage + "\nContinue?"
            dlg = wx.MessageDialog(self, SomeOfOutputFilesExistsMessage, "Error In Input Values", wx.OK | wx.CANCEL)
            if dlg.ShowModal() != wx.ID_OK:
                dlg.Destroy()
                return False
            dlg.Destroy()

        if self.FileBB_PatchMap_SysTopFileNameOUT.GetValue() == "None":
            if self.FileBB_PatchMap_DiffFileNameOUT.GetValue() == "None":
                print("Nothing to do")
                print("done with Patching Maps")
                print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
                return
        # Do the magic
        if self.FileBB_PatchMap_DiffFileNameIN.GetValue() != "None":
            contworld = pnps.LoadContWorld(
                SysTop=self.FileBB_PatchMap_SysTopFileNameIN.GetValue(),
                Diffusion=self.FileBB_PatchMap_DiffFileNameIN.GetValue()
            )
        else:
            contworld = pnps.LoadContWorld(
                SysTop=self.FileBB_PatchMap_SysTopFileNameIN.GetValue()
            )
        if contworld == None:
            return False
        patcher = pnps.DielDiffMapsPatcher(contworld)

        DielPatchTable = self.Table_PatchMap_Diel.GetTable()
        Table_PatchMap_DielInUse = self.Table_PatchMap_DielInUse.GetTable()

        if self.FileBB_PatchMap_SysTopFileNameOUT.GetValue() != "None":
            # Set New Epsilons Values in Use
            patcher.setDielConstInUse(self.Table_PatchMap_DielInUse.GetTable().data[0])

            # To through table apply patches for dielectric constant
            patcher.GetIntDielMaps()
            for i in range(DielPatchTable.GetNumberRows()):
                epsToOver = DielPatchTable.DielConstInUseStrRep.index(DielPatchTable.GetValue(i, 0))
                epsNew = DielPatchTable.DielConstInUseStrRep.index(DielPatchTable.GetValue(i, 1))
                x = DielPatchTable.GetValue(i, 2)
                y = DielPatchTable.GetValue(i, 3)
                z1 = DielPatchTable.GetValue(i, 4)
                z2 = DielPatchTable.GetValue(i, 5)
                R = DielPatchTable.GetValue(i, 6)
                patcher.PatchDielMaps(epsToOver, epsNew, x, y, z1, z2, R)
            print("===============================================================================")
            patcher.PushIntDielMaps()

        if self.FileBB_PatchMap_DiffFileNameOUT.GetValue() != "None":
            # Patch Diffusion
            DiffPatchTable = self.Table_PatchMap_Diff.GetTable()
            patcher.InitNewDiff()
            for i in range(DiffPatchTable.GetNumberRows()):
                Dscale1K = DiffPatchTable.GetValue(i, 0)
                Dscale2K = DiffPatchTable.GetValue(i, 1)
                Dscale1Cl = DiffPatchTable.GetValue(i, 2)
                Dscale2Cl = DiffPatchTable.GetValue(i, 3)
                x = DiffPatchTable.GetValue(i, 4)
                y = DiffPatchTable.GetValue(i, 5)
                z1 = DiffPatchTable.GetValue(i, 6)
                z2 = DiffPatchTable.GetValue(i, 7)
                R = DiffPatchTable.GetValue(i, 8)
                patcher.PatchDiffMaps(Dscale1K, Dscale2K, Dscale1Cl, Dscale2Cl, x, y, z1, z2, R)
            print("===============================================================================")
            patcher.DelRefDiff()

        if self.FileBB_PatchMap_SysTopFileNameOUT.GetValue() != "None":
            contworld.WriteNodeIndexing(self.FileBB_PatchMap_SysTopFileNameOUT.GetValue())
        if self.FileBB_PatchMap_DiffFileNameOUT.GetValue() != "None":
            contworld.WriteDiffusion(self.FileBB_PatchMap_DiffFileNameOUT.GetValue())
        del patcher
        del contworld
        # Read maps values
        print("done with Patching Maps")
        print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")

    def PatchMaps_OnViewStaticCharge(self, event):
        if not FileBB_CheckFileExistence(self.FileBB_PatchMap_SysTopFileNameOUT):
            dlg = wx.MessageDialog(self, "File %s does not exist" % (self.FileBB_PatchMap_SysTopFileNameOUT.GetValue()),
                                   "Error In Input Values", wx.OK)
            dlg.ShowModal()  # Show it
            dlg.Destroy()  # finally destroy it when finished.
            return False

        NIndex = pnps.NodeIndexing()
        NIndex.ReadFromFile(self.FileBB_PatchMap_SysTopFileNameOUT.GetValue())

        field = NIndex.GetHaField3D(pnps.NodeIndexing.Charge, pnps.NodeIndexing.ChargeMask);
        Name = "Static Charge (%s , %s)" % (self.Str_CM_SysTopFileName.GetValue(), self.pmset.GetName())
        field.SetName(Name)
        field.thisown = 0
        PlaneV = molset.PlaneViewOfHaField3D(field, Name, 1)
        PlaneV.SetHideZeroValues(1);
        PlaneV.thisown = 0

        self.pmset.AddObject3D(PlaneV)
        self.pmset.RefreshAllViews(
            molset.RFRefresh | molset.RFColour | molset.RFApply | molset.RFMagnify | molset.RFRotate | molset.RFTransX | molset.RFTransY)
        self.pmset.RefreshAllViews(
            molset.RFRefresh | molset.RFColour | molset.RFApply | molset.RFMagnify | molset.RFRotate | molset.RFTransX | molset.RFTransY)

        molset.CreatewxFieldPlaneView(PlaneV, self.pmset, Name, 1)
        del NIndex

    def PatchMaps_OnViewDiel(self, event):
        if not FileBB_CheckFileExistence(self.FileBB_PatchMap_SysTopFileNameOUT):
            dlg = wx.MessageDialog(self, "File %s does not exist" % (self.FileBB_PatchMap_SysTopFileNameOUT.GetValue()),
                                   "Error In Input Values", wx.OK)
            dlg.ShowModal()  # Show it
            dlg.Destroy()  # finally destroy it when finished.
            return False

        NIndex = pnps.NodeIndexing()
        NIndex.ReadFromFile(self.FileBB_PatchMap_SysTopFileNameOUT.GetValue())

        Axes = ["X", "Y", "Z"]
        Mask = [pnps.NodeIndexing.Epsilon0, pnps.NodeIndexing.Epsilon1, pnps.NodeIndexing.Epsilon2]
        for i in range(3):
            field = NIndex.GetHaField3D(pnps.NodeIndexing.DielConst, Mask[i]);
            Name = "Diel.Const; Displaced in %s Charge (%s , %s)" % (
            Axes[i], self.Str_CM_SysTopFileName.GetValue(), self.pmset.GetName())
            field.SetName(Name)
            field.thisown = 0
            PlaneV = molset.PlaneViewOfHaField3D(field, Name, 1)
            PlaneV.SetHideZeroValues(1);
            PlaneV.thisown = 0

            self.pmset.AddObject3D(PlaneV)
            self.pmset.RefreshAllViews(molset.RFRefresh)

            molset.CreatewxFieldPlaneView(PlaneV, self.pmset, Name, 1)
        self.pmset.RefreshAllViews(
            molset.RFRefresh | molset.RFColour | molset.RFApply | molset.RFMagnify | molset.RFRotate | molset.RFTransX | molset.RFTransY)
        self.pmset.RefreshAllViews(
            molset.RFRefresh | molset.RFColour | molset.RFApply | molset.RFMagnify | molset.RFRotate | molset.RFTransX | molset.RFTransY)
        del NIndex

    def PatchMaps_OnViewPlainDiff(self, event):
        if not FileBB_CheckFileExistence(self.FileBB_PatchMap_SysTopFileNameOUT):
            dlg = wx.MessageDialog(self, "File %s does not exist" % (self.FileBB_PatchMap_SysTopFileNameOUT.GetValue()),
                                   "Error In Input Values", wx.OK)
            dlg.ShowModal()  # Show it
            dlg.Destroy()  # finally destroy it when finished.
            return False

        NIndex = pnps.NodeIndexing()
        NIndex.ReadFromFile(self.FileBB_PatchMap_SysTopFileNameOUT.GetValue())

        Mask = [pnps.NodeIndexing.Ion0, pnps.NodeIndexing.Ion1]
        for i in range(3):
            field = NIndex.GetHaField3D(pnps.NodeIndexing.DiffConst, Mask[i]);
            Name = "Plain Diffusion; Ion%d (%s , %s)" % (i, self.Str_CM_SysTopFileName.GetValue(), self.pmset.GetName())
            field.SetName(Name)
            field.thisown = 0
            PlaneV = molset.PlaneViewOfHaField3D(field, Name, 1)
            PlaneV.SetHideZeroValues(1);
            PlaneV.thisown = 0

            self.pmset.AddObject3D(PlaneV)
            self.pmset.RefreshAllViews(molset.RFRefresh)

            molset.CreatewxFieldPlaneView(PlaneV, self.pmset, Name, 1)
        self.pmset.RefreshAllViews(
            molset.RFRefresh | molset.RFColour | molset.RFApply | molset.RFMagnify | molset.RFRotate | molset.RFTransX | molset.RFTransY)
        self.pmset.RefreshAllViews(
            molset.RFRefresh | molset.RFColour | molset.RFApply | molset.RFMagnify | molset.RFRotate | molset.RFTransX | molset.RFTransY)
        del NIndex

    def PatchMaps_OnViewDiff(self, event):
        ShowVectorField3D(self.FileBB_PatchMap_DiffFileNameOUT.GetValue(), self.pmset, "Diffusion")

    ###############################################################################################
    # IAV Refinement Interface
    def IAVRef_InitWX(self):
        # Parameters
        self.TextCtrl_IAV_MaxdC = self.FindWindowById(IDTCpnpIAV_MaxdC)
        self.NumCtrl_IAV_MaxCycles = self.FindWindowById(IDTCpnpIAV_MaxCycles)
        self.NumCtrl_IAV_NP_Relax = self.FindWindowById(IDTCpnpIAV_NP_Relax)
        self.TextCtrl_IAV_PBSR_Tolerance = self.FindWindowById(IDTCpnpIAV_PBSR_Tolerance)
        self.NumCtrl_IAV_PBSR_Niter = self.FindWindowById(IDTCpnpIAV_PBSR_Niter)
        self.NumCtrl_IAV_PBSR_Relax = self.FindWindowById(IDTCpnpIAV_PBSR_Relax)

        # To validate user entry
        self.PBSR_NumCtrls4Validation = [
            self.NumCtrl_IAV_MaxCycles, self.NumCtrl_IAV_NP_Relax,
            self.NumCtrl_IAV_PBSR_Niter, self.NumCtrl_IAV_PBSR_Relax]

        # Input Maps
        self.FileBB_IAV_SysTopFileNameIN = self.FindWindowById(IDTCpnpIAV_SysTopFileNameIN)
        self.FileBB_IAV_DiffFileNameIN = self.FindWindowById(IDTCpnpIAV_DiffFileNameIN)
        self.FileBB_IAV_SRFileName = self.FindWindowById(IDTCpnpIAV_SRFileName)
        self.FileBB_IAV_PotFileNameIN = self.FindWindowById(IDTCpnpIAV_PotFileNameIN)
        self.FileBB_IAV_ConcFileNameIN = self.FindWindowById(IDTCpnpIAV_ConcFileNameIN)
        # To validate user entry
        self.IAV_FileIn = [self.FileBB_IAV_SysTopFileNameIN,
                           self.FileBB_IAV_PotFileNameIN, self.FileBB_IAV_ConcFileNameIN]
        self.IAV_FileInNoneAllowed = [self.FileBB_IAV_DiffFileNameIN, self.FileBB_IAV_SRFileName]
        # Output Maps
        self.FileBB_IAV_SysTopFileNameOUT = self.FindWindowById(IDTCpnpIAV_SysTopFileNameOUT)
        self.FileBB_IAV_DiffFileNameOUT = self.FindWindowById(IDTCpnpIAV_DiffFileNameOUT)
        self.FileBB_IAV_PotFileNameOUT = self.FindWindowById(IDTCpnpIAV_PotFileNameOUT)
        self.FileBB_IAV_ConcFileNameOUT = self.FindWindowById(IDTCpnpIAV_ConcFileNameOUT)
        # To validate user entry
        self.IAV_FileOUT = [self.FileBB_IAV_SysTopFileNameOUT,
                            self.FileBB_IAV_PotFileNameOUT, self.FileBB_IAV_ConcFileNameOUT]
        self.IAV_FileOUTNoneAllowed = [self.FileBB_IAV_DiffFileNameOUT]

        # bind events
        self.Bind(wx.EVT_BUTTON, self.IAVRef_OnRun, id=IDBpnpIAV_Run)

        self.Bind(wx.EVT_BUTTON, self.IAVRef_OnViewPlainDiff, id=IDBpnpIAV_ViewPlainDiff)
        self.Bind(wx.EVT_BUTTON, self.IAVRef_OnViewDiff, id=IDBpnpIAV_ViewDiff)
        self.Bind(wx.EVT_BUTTON, self.IAVRef_OnViewPot, id=IDBpnpIAV_ViewPot)
        self.Bind(wx.EVT_BUTTON, self.IAVRef_OnViewConc, id=IDBpnpIAV_ViewConc)

    def IAVRef_OnRun(self, event):
        # Check input values
        # check the existance of input map file
        SomeOfInputFilesDoNotExists = False
        SomeOfInputFilesDoNotExistsMessage = ""
        #
        for filebrowser in self.IAV_FileIn:
            if not FileBB_CheckFileExistence(filebrowser):
                SomeOfInputFilesDoNotExists = True
                SomeOfInputFilesDoNotExistsMessage = SomeOfInputFilesDoNotExistsMessage + "file %s do not exist.\n" % (
                    filebrowser.GetValue())
        for filebrowser in self.IAV_FileInNoneAllowed:
            if not FileBB_CheckFileExistence(filebrowser, True):
                SomeOfInputFilesDoNotExists = True
                SomeOfInputFilesDoNotExistsMessage = SomeOfInputFilesDoNotExistsMessage + "file %s do not exist.\n" % (
                    filebrowser.GetValue())
        if SomeOfInputFilesDoNotExists:
            SomeOfInputFilesDoNotExistsMessage = "Some of input files do not exists:" + SomeOfInputFilesDoNotExistsMessage
            dlg = wx.MessageDialog(self, SomeOfInputFilesDoNotExistsMessage, "Error In Input Values", wx.OK)
            dlg.ShowModal()  # Show it
            dlg.Destroy()  # finally destroy it when finished.
            return False
        # check the existance of these file and ask user for permition to rewrite
        SomeOfOutputFilesExists = False
        SomeOfOutputFilesExistsMessage = ""
        for filebrowser in self.IAV_FileOUT:
            if FileBB_CheckFileExistence(filebrowser):
                SomeOfOutputFilesExists = True
                SomeOfOutputFilesExistsMessage = SomeOfOutputFilesExistsMessage + "file %s exist.\n" % (
                    filebrowser.GetValue())
        for filebrowser in self.IAV_FileOUTNoneAllowed:
            if FileBB_CheckFileExistence(filebrowser, True):
                SomeOfOutputFilesExists = True
                SomeOfOutputFilesExistsMessage = SomeOfOutputFilesExistsMessage + "file %s exist.\n" % (
                    filebrowser.GetValue())
        if SomeOfOutputFilesExists:
            SomeOfOutputFilesExistsMessage = SomeOfOutputFilesExistsMessage + "\nIf execution will be continued this(these) file(s) will be deleted."
            SomeOfOutputFilesExistsMessage = SomeOfOutputFilesExistsMessage + "\nContinue?"
            dlg = wx.MessageDialog(self, SomeOfOutputFilesExistsMessage, "Error In Input Values", wx.OK | wx.CANCEL)
            if dlg.ShowModal() != wx.ID_OK:
                dlg.Destroy()
                return False
            dlg.Destroy()

        # Parameters
        # Valudate Other entries
        for ctrl in self.PBSR_NumCtrls4Validation:
            if not ctrl.IsInBounds():
                dlg = wx.MessageDialog(self,
                                       "Some input values are incorrect (some can be marked by yellow background)",
                                       "Error In Input Values", wx.OK)
                dlg.ShowModal()  # Show it
                dlg.Destroy()  # finally destroy it when finished.
                return False
        # TextCtrl_IAV_MaxdC
        TolValCorrect = True
        strTol = self.TextCtrl_IAV_MaxdC.GetValue()
        try:
            Tol = float(strTol)
            if Tol < 0.0:
                TolValCorrect = False
        except ValueError:
            TolValCorrect = False

        if not TolValCorrect:
            dlg = wx.MessageDialog(self, "Maximal Allowed Concentration change value is incorrect",
                                   "Error In Input Values", wx.OK)
            dlg.ShowModal()  # Show it
            dlg.Destroy()  # finally destroy it when finished.
            return False
        #
        TolValCorrect = True
        strTol = self.TextCtrl_IAV_PBSR_Tolerance.GetValue()
        try:
            Tol = float(strTol)
            if Tol < 0.0:
                TolValCorrect = False
        except ValueError:
            TolValCorrect = False

        if not TolValCorrect:
            dlg = wx.MessageDialog(self, "Convergence(Tolerance) value is incorrect", "Error In Input Values", wx.OK)
            dlg.ShowModal()  # Show it
            dlg.Destroy()  # finally destroy it when finished.
            return False
        # Get User Input
        # Parameters
        MaxdC = float(self.TextCtrl_IAV_MaxdC.GetValue())
        MaxCycles = self.NumCtrl_IAV_MaxCycles.GetValue()
        NP_Relax = self.NumCtrl_IAV_NP_Relax.GetValue()
        PBSR_Tolerance = float(self.TextCtrl_IAV_PBSR_Tolerance.GetValue())
        PBSR_Niter = self.NumCtrl_IAV_PBSR_Niter.GetValue()
        PBSR_Relax = self.NumCtrl_IAV_PBSR_Relax.GetValue()

        # Input Maps
        SysTopFileNameIN = self.FileBB_IAV_SysTopFileNameIN.GetValue()
        DiffFileNameIN = self.FileBB_IAV_DiffFileNameIN.GetValue()
        if DiffFileNameIN == "None":
            DiffFileNameIN = None
        SRFileName = self.FileBB_IAV_SRFileName.GetValue()
        if SRFileName == "None":
            SRFileName = None
        PotFileNameIN = self.FileBB_IAV_PotFileNameIN.GetValue()
        ConcFileNameIN = self.FileBB_IAV_ConcFileNameIN.GetValue()
        # Output Maps
        SysTopFileNameOUT = self.FileBB_IAV_SysTopFileNameOUT.GetValue()
        DiffFileNameOUT = self.FileBB_IAV_DiffFileNameOUT.GetValue()
        if DiffFileNameOUT == "None":
            DiffFileNameOUT = None
        PotFileNameOUT = self.FileBB_IAV_PotFileNameOUT.GetValue()
        ConcFileNameOUT = self.FileBB_IAV_ConcFileNameOUT.GetValue()
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print("Running IAV Refinement ...")
        # Load maps
        contworld = pnps.LoadContWorld(SysTop=SysTopFileNameIN, PMF=SRFileName,
                                         Potential=PotFileNameIN, Concentration=ConcFileNameIN,
                                         Diffusion=DiffFileNameIN)
        if contworld == None:
            print("Error: contworld was not allocated")
            return False

        pnps.PNPUtil.ConvertPBLJresultsToDynamicCharge(contworld)
        pnps.RefineIAV(contworld,
                         Relaxation=1.0,
                         RunPBLJwNPCriteriaCycle=True,
                         MaxdC=MaxdC,
                         MaxCycles=MaxCycles,
                         RunSuccessfullCycles=40,
                         PBSR_Param={"MaxIterations": PBSR_Niter, "Relaxation": PBSR_Relax, "Tolerance": PBSR_Tolerance,
                                     "Verbose": True},
                         NP_Param={"MaxIterations": 1, "Relaxation": NP_Relax, "Tolerance": 0.0, "ConvergenceCheck": 1}
                         )
        pnps.SingleWorldCalcController.CmdSetDZeroWithIndex(contworld)

        contworld.WriteNodeIndexing(SysTopFileNameOUT)
        contworld.WriteDiffusion(DiffFileNameOUT)
        contworld.WritePotential(PotFileNameOUT)
        contworld.WriteDynamicCharge(ConcFileNameOUT)
        # free allocations
        del contworld
        print("done IAV Refinement")
        print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")

    def IAVRef_OnViewPlainDiff(self, event):
        if not FileBB_CheckFileExistence(self.FileBB_IAV_SysTopFileNameOUT):
            dlg = wx.MessageDialog(self, "File %s does not exist" % (self.FileBB_IAV_SysTopFileNameOUT.GetValue()),
                                   "Error In Input Values", wx.OK)
            dlg.ShowModal()  # Show it
            dlg.Destroy()  # finally destroy it when finished.
            return False

        NIndex = pnps.NodeIndexing()
        NIndex.ReadFromFile(self.FileBB_IAV_SysTopFileNameOUT.GetValue())

        Mask = [pnps.NodeIndexing.Ion0, pnps.NodeIndexing.Ion1]
        for i in range(3):
            field = NIndex.GetHaField3D(pnps.NodeIndexing.DiffConst, Mask[i]);
            Name = "Plain Diffusion; Ion%d (%s , %s)" % (i, self.Str_CM_SysTopFileName.GetValue(), self.pmset.GetName())
            field.SetName(Name)
            field.thisown = 0
            PlaneV = molset.PlaneViewOfHaField3D(field, Name, 1)
            PlaneV.SetHideZeroValues(1);
            PlaneV.thisown = 0

            self.pmset.AddObject3D(PlaneV)
            self.pmset.RefreshAllViews(molset.RFRefresh)

            molset.CreatewxFieldPlaneView(PlaneV, self.pmset, Name, 1)
        self.pmset.RefreshAllViews(
            molset.RFRefresh | molset.RFColour | molset.RFApply | molset.RFMagnify | molset.RFRotate | molset.RFTransX | molset.RFTransY)
        self.pmset.RefreshAllViews(
            molset.RFRefresh | molset.RFColour | molset.RFApply | molset.RFMagnify | molset.RFRotate | molset.RFTransX | molset.RFTransY)
        del NIndex

    def IAVRef_OnViewDiff(self, event):
        ShowVectorField3D(self.FileBB_IAV_DiffFileNameIN.GetValue(), self.pmset, "Diffusion")

    def IAVRef_OnViewPot(self, event):
        ShowVectorField3D(self.FileBB_IAV_PotFileNameOUT.GetValue(), self.pmset, "Potential")

    def IAVRef_OnViewConc(self, event):
        ShowVectorField3D(self.FileBB_IAV_ConcFileNameOUT.GetValue(), self.pmset, "Concentrations")

    ###############################################################################################
    # PNP-SR Interface
    def PNPSR_InitWX(self):
        # Parameters
        # Mode for application of external potential mode
        self.NumCtrl_PNPSR_PotDiff = self.FindWindowById(IDTpnpPNPSR_PotDiff)

        self.NumCtrl_PNPSR_z0 = self.FindWindowById(IDTpnpPNPSR_z0)
        self.NumCtrl_PNPSR_z1 = self.FindWindowById(IDTpnpPNPSR_z1)
        self.NumCtrl_PNPSR_Pot0 = self.FindWindowById(IDTpnpPNPSR_Pot0)
        self.NumCtrl_PNPSR_Pot1 = self.FindWindowById(IDTpnpPNPSR_Pot1)

        self.PNPSR_AppPotMode_AloneMem_NumCtrls4Validation = [self.NumCtrl_PNPSR_PotDiff]
        self.PNPSR_Manual_NumCtrls4Validation = [
            self.NumCtrl_PNPSR_z0, self.NumCtrl_PNPSR_z1,
            self.NumCtrl_PNPSR_Pot0, self.NumCtrl_PNPSR_Pot1]

        # lables for AppPotMode
        self.Labels_PNPSR_AppPotMode_AloneMem = []
        self.Labels_PNPSR_AppPotMode_AloneMem.append(self.FindWindowById(IDTpnpPNPSR_AppPotMode_AloneMem1))
        self.Labels_PNPSR_AppPotMode_AloneMem.append(self.FindWindowById(IDTpnpPNPSR_AppPotMode_AloneMem2))
        self.Labels_PNPSR_AppPotMode_Manual = []
        self.Labels_PNPSR_AppPotMode_Manual.append(self.FindWindowById(IDTpnpPNPSR_AppPotMode_Manual1))
        self.Labels_PNPSR_AppPotMode_Manual.append(self.FindWindowById(IDTpnpPNPSR_AppPotMode_Manual2))
        self.Labels_PNPSR_AppPotMode_Manual.append(self.FindWindowById(IDTpnpPNPSR_AppPotMode_Manual3))
        self.Labels_PNPSR_AppPotMode_Manual.append(self.FindWindowById(IDTpnpPNPSR_AppPotMode_Manual4))
        self.Labels_PNPSR_AppPotMode_Manual.append(self.FindWindowById(IDTpnpPNPSR_AppPotMode_Manual5))

        # PNP param
        self.NumCtrl_PNPSR_MaxPNPCylces = self.FindWindowById(IDTpnpPNPSR_MaxPNPCylces)
        self.NumCtrl_PNPSR_PS_Niter = self.FindWindowById(IDTpnpPNPSR_PS_Niter)
        self.NumCtrl_PNPSR_PS_Relax = self.FindWindowById(IDTpnpPNPSR_PS_Relax)
        self.NumCtrl_PNPSR_NPS_Niter = self.FindWindowById(IDTpnpPNPSR_NPS_Niter)
        self.NumCtrl_PNPSR_NPS_Relax = self.FindWindowById(IDTpnpPNPSR_NPS_Relax)
        self.TextCtrl_PNPSR_CurrentTolerance = self.FindWindowById(IDTpnpPNPSR_CurrentTolerance)
        # To validate user entry
        self.PNPSR_NumCtrls4Validation = [
            self.NumCtrl_PNPSR_MaxPNPCylces, self.NumCtrl_PNPSR_PS_Niter,
            self.NumCtrl_PNPSR_PS_Relax, self.NumCtrl_PNPSR_NPS_Niter,
            self.NumCtrl_PNPSR_NPS_Relax]

        # Selection of potential application mode
        self.Radio_PNPSR_AppPotMode_AloneMem = self.FindWindowById(IDRpnpPNPSR_AppPotMode_AloneMem)
        self.Radio_PNPSR_AppPotMode_Manual = self.FindWindowById(IDRpnpPNPSR_AppPotMode_Manual)
        self.Bind(wx.EVT_RADIOBUTTON, self.PNPSR_OnChangeAppPotMode, self.Radio_PNPSR_AppPotMode_AloneMem)
        self.Bind(wx.EVT_RADIOBUTTON, self.PNPSR_OnChangeAppPotMode, self.Radio_PNPSR_AppPotMode_Manual)
        self.PNPSR_OnChangeAppPotMode(None)

        # Input Maps
        self.FileBB_PNPSR_SysTopFileNameIN = self.FindWindowById(IDTCpnpPNPSR_SysTopFileNameIN)
        self.FileBB_PNPSR_DiffFileNameIN = self.FindWindowById(IDTCpnpPNPSR_DiffFileNameIN)
        self.FileBB_PNPSR_SRFileName = self.FindWindowById(IDTCpnpPNPSR_SRFileName)
        self.FileBB_PNPSR_PotFileNameIN = self.FindWindowById(IDTCpnpPNPSR_PotFileNameIN)
        self.FileBB_PNPSR_ConcFileNameIN = self.FindWindowById(IDTCpnpPNPSR_ConcFileNameIN)
        # To validate user entry
        self.PNPSR_FileIn = [self.FileBB_PNPSR_SysTopFileNameIN]
        self.PNPSR_FileInNoneAllowed = [
            self.FileBB_PNPSR_DiffFileNameIN, self.FileBB_PNPSR_SRFileName,
            self.FileBB_PNPSR_PotFileNameIN, self.FileBB_PNPSR_ConcFileNameIN]
        # Output Maps
        self.FileBB_PNPSR_PotFileNameOUT = self.FindWindowById(IDTCpnpPNPSR_PotFileNameOUT)
        self.FileBB_PNPSR_ConcFileNameOUT = self.FindWindowById(IDTCpnpPNPSR_ConcFileNameOUT)
        # To validate user entry
        self.PNPSR_FileOUT = []
        self.PNPSR_FileOUTNoneAllowed = [self.FileBB_PNPSR_PotFileNameOUT, self.FileBB_PNPSR_ConcFileNameOUT]

        # bind events
        self.Bind(wx.EVT_BUTTON, self.PNPSR_OnRun, id=IDBpnpPNPSR_Run)

        self.Bind(wx.EVT_BUTTON, self.PNPSR_OnViewPot, id=IDBpnpPNPSR_ViewPot)
        self.Bind(wx.EVT_BUTTON, self.PNPSR_OnViewConc, id=IDBpnpPNPSR_ViewConc)

    def PNPSR_OnChangeAppPotMode(self, event):
        if self.Radio_PNPSR_AppPotMode_AloneMem.GetValue():
            for l in self.Labels_PNPSR_AppPotMode_AloneMem:
                l.Enable(True)
            for l in self.Labels_PNPSR_AppPotMode_Manual:
                l.Enable(False)
            self.NumCtrl_PNPSR_PotDiff.Enable(True)

            self.NumCtrl_PNPSR_z0.Enable(False)
            self.NumCtrl_PNPSR_z1.Enable(False)
            self.NumCtrl_PNPSR_Pot0.Enable(False)
            self.NumCtrl_PNPSR_Pot1.Enable(False)
        else:
            for l in self.Labels_PNPSR_AppPotMode_AloneMem:
                l.Enable(False)
            for l in self.Labels_PNPSR_AppPotMode_Manual:
                l.Enable(True)
            self.NumCtrl_PNPSR_PotDiff.Enable(False)

            self.NumCtrl_PNPSR_z0.Enable(True)
            self.NumCtrl_PNPSR_z1.Enable(True)
            self.NumCtrl_PNPSR_Pot0.Enable(True)
            self.NumCtrl_PNPSR_Pot1.Enable(True)

    def PNPSR_OnRun(self, event):
        # Check input values
        # check the existance of input map file
        SomeOfInputFilesDoNotExists = False
        SomeOfInputFilesDoNotExistsMessage = ""
        #
        for filebrowser in self.PNPSR_FileIn:
            if not FileBB_CheckFileExistence(filebrowser):
                SomeOfInputFilesDoNotExists = True
                SomeOfInputFilesDoNotExistsMessage = SomeOfInputFilesDoNotExistsMessage + "file %s do not exist.\n" % (
                    filebrowser.GetValue())
        for filebrowser in self.PNPSR_FileInNoneAllowed:
            if not FileBB_CheckFileExistence(filebrowser, True):
                SomeOfInputFilesDoNotExists = True
                SomeOfInputFilesDoNotExistsMessage = SomeOfInputFilesDoNotExistsMessage + "file %s do not exist.\n" % (
                    filebrowser.GetValue())
        if SomeOfInputFilesDoNotExists:
            SomeOfInputFilesDoNotExistsMessage = "Some of input files do not exists:" + SomeOfInputFilesDoNotExistsMessage
            dlg = wx.MessageDialog(self, SomeOfInputFilesDoNotExistsMessage, "Error In Input Values", wx.OK)
            dlg.ShowModal()  # Show it
            dlg.Destroy()  # finally destroy it when finished.
            return False

        # check the existance of these file and ask user for permition to rewrite
        SomeOfOutputFilesExists = False
        SomeOfOutputFilesExistsMessage = ""
        for filebrowser in self.PNPSR_FileOUT:
            if FileBB_CheckFileExistence(filebrowser):
                SomeOfOutputFilesExists = True
                SomeOfOutputFilesExistsMessage = SomeOfOutputFilesExistsMessage + "file %s exist.\n" % (
                    filebrowser.GetValue())
        for filebrowser in self.PNPSR_FileOUTNoneAllowed:
            if FileBB_CheckFileExistence(filebrowser, True):
                SomeOfOutputFilesExists = True
                SomeOfOutputFilesExistsMessage = SomeOfOutputFilesExistsMessage + "file %s exist.\n" % (
                    filebrowser.GetValue())
        if SomeOfOutputFilesExists:
            SomeOfOutputFilesExistsMessage = SomeOfOutputFilesExistsMessage + "\nIf execution will be continued this(these) file(s) will be deleted."
            SomeOfOutputFilesExistsMessage = SomeOfOutputFilesExistsMessage + "\nContinue?"
            dlg = wx.MessageDialog(self, SomeOfOutputFilesExistsMessage, "Error In Input Values", wx.OK | wx.CANCEL)
            if dlg.ShowModal() != wx.ID_OK:
                dlg.Destroy()
                return False
            dlg.Destroy()

        # Parameters
        # Valudate Other entries
        for ctrl in self.PNPSR_NumCtrls4Validation:
            if not ctrl.IsInBounds():
                dlg = wx.MessageDialog(self,
                                       "Some input values are incorrect (some can be marked by yellow background)",
                                       "Error In Input Values", wx.OK)
                dlg.ShowModal()  # Show it
                dlg.Destroy()  # finally destroy it when finished.
                return False

        if self.Radio_PNPSR_AppPotMode_AloneMem.GetValue():
            for ctrl in self.PNPSR_AppPotMode_AloneMem_NumCtrls4Validation:
                if not ctrl.IsInBounds():
                    dlg = wx.MessageDialog(self,
                                           "Some input values are incorrect (some can be marked by yellow background)",
                                           "Error In Input Values", wx.OK)
                    dlg.ShowModal()  # Show it
                    dlg.Destroy()  # finally destroy it when finished.
                    return False
        else:
            for ctrl in self.PNPSR_Manual_NumCtrls4Validation:
                if not ctrl.IsInBounds():
                    dlg = wx.MessageDialog(self,
                                           "Some input values are incorrect (some can be marked by yellow background)",
                                           "Error In Input Values", wx.OK)
                    dlg.ShowModal()  # Show it
                    dlg.Destroy()  # finally destroy it when finished.
                    return False
        #
        TolValCorrect = True
        strTol = self.TextCtrl_PNPSR_CurrentTolerance.GetValue()
        try:
            Tol = float(strTol)
            if Tol < 0.0:
                TolValCorrect = False
        except ValueError:
            TolValCorrect = False

        if not TolValCorrect:
            dlg = wx.MessageDialog(self, "Convergence(Tolerance) value is incorrect", "Error In Input Values", wx.OK)
            dlg.ShowModal()  # Show it
            dlg.Destroy()  # finally destroy it when finished.
            return False

        # Get User Input
        # Parameters
        MaxPNPCylces = self.NumCtrl_PNPSR_MaxPNPCylces.GetValue()
        PS_Niter = self.NumCtrl_PNPSR_PS_Niter.GetValue()
        PS_Relax = self.NumCtrl_PNPSR_PS_Relax.GetValue()
        NPS_Niter = self.NumCtrl_PNPSR_NPS_Niter.GetValue()
        NPS_Relax = self.NumCtrl_PNPSR_NPS_Relax.GetValue()
        CurrentTolerance = float(self.TextCtrl_PNPSR_CurrentTolerance.GetValue())

        AppPotMode_AloneMem = False
        if self.Radio_PNPSR_AppPotMode_AloneMem.GetValue():
            AppPotMode_AloneMem = True

        PotDiff = self.NumCtrl_PNPSR_PotDiff.GetValue()

        Pot_z0 = self.NumCtrl_PNPSR_z0.GetValue()
        Pot_z1 = self.NumCtrl_PNPSR_z1.GetValue()
        Pot0 = self.NumCtrl_PNPSR_Pot0.GetValue()
        Pot1 = self.NumCtrl_PNPSR_Pot1.GetValue()

        # Input Maps
        SysTopFileNameIN = self.FileBB_PNPSR_SysTopFileNameIN.GetValue()
        DiffFileNameIN = self.FileBB_PNPSR_DiffFileNameIN.GetValue()
        if DiffFileNameIN == "None":
            DiffFileNameIN = None
        SRFileName = self.FileBB_PNPSR_SRFileName.GetValue()
        if SRFileName == "None":
            SRFileName = None
        PotFileNameIN = self.FileBB_PNPSR_PotFileNameIN.GetValue()
        if PotFileNameIN == "None":
            PotFileNameIN = None
        ConcFileNameIN = self.FileBB_PNPSR_ConcFileNameIN.GetValue()
        if ConcFileNameIN == "None":
            ConcFileNameIN = None
        # Output Maps
        PotFileNameOUT = self.FileBB_PNPSR_PotFileNameOUT.GetValue()
        if PotFileNameOUT == "None":
            PotFileNameOUT = None
        ConcFileNameOUT = self.FileBB_PNPSR_ConcFileNameOUT.GetValue()
        if ConcFileNameOUT == "None":
            ConcFileNameOUT = None
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print("Running PNP-SR ...")
        # Load maps
        contworld = pnps.LoadContWorld(SysTop=SysTopFileNameIN, PMF=SRFileName,
                                         Potential=PotFileNameIN, Concentration=ConcFileNameIN,
                                         Diffusion=DiffFileNameIN)
        if contworld == None:
            print("Error: contworld was not allocated")
            return False

        # do PNP-SR
        LimitCurrentCalcZ = [Pot_z0, Pot_z1]
        if (AppPotMode_AloneMem):
            LimitCurrentCalcZ = contworld.AddPotentialAuto(PotDiff)
            if LimitCurrentCalcZ == None:
                print("Error: can not locate implicit membrane")
                del contworld
                return False
        else:
            contworld.AddPotential(Pot_z0, Pot_z1, Pot0, Pot1)

        pnps.SolvePNPSR(contworld, MaxIterations=MaxPNPCylces, Verbose=False, LimitCurrentCalcZ=LimitCurrentCalcZ,
                          P_Param={"MaxIterations": PS_Niter, "Relaxation": PS_Relax, "Tolerance": 0.0,
                                   "Verbose": False},
                          NP_Param={"MaxIterations": NPS_Niter, "Relaxation": NPS_Relax, "Tolerance": 0.0,
                                    "Verbose": False}
                          )
        # save maps
        contworld.WritePotential(PotFileNameOUT)
        contworld.WriteDynamicCharge(ConcFileNameOUT)
        # free allocations
        del contworld
        print("done PNP-SR")
        print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")

    def PNPSR_OnViewPot(self, event):
        ShowVectorField3D(str(self.FileBB_PNPSR_PotFileNameOUT.GetValue()), self.pmset, "Potential")

    def PNPSR_OnViewConc(self, event):
        ShowVectorField3D(str(self.FileBB_PNPSR_ConcFileNameOUT.GetValue()), self.pmset, "Concentrations")

#
