# Load benz_1.hlm, initialize QChem module and execute this script 
#
pmset = GetCurMolSet()
qcmod = pmset.GetQCMod(1)
fout = open("basis.xml","w")
qcmod.AtBasis.SaveXml(fout)
fout.close()
