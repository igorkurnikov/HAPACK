# Load tb5.hlm and execute this script 
#
pmset = GetCurMolSet()
fout = open("tb5.xml","w")
pmset.SaveXml(fout)
fout.close()
