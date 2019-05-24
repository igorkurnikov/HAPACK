from harlempy import halib
from harlempy.halib import *
from harlempy import molset
from harlempy.molset import *
import pnpmod
from pnpmod import *

#Load QRDB
print("Loading QRDB-AMBER...")
QRDB=GetQRDB()
ff=QRDB.NewFF("AMBER??","Ref to AMBER","Here is Some Notes")

r_S =2.0000
r_SH=2.0000
r_P =2.1000
r_IM=2.4700
r_Li=1.1370
r_IP=1.8680
r_K =2.6580
r_Rb=2.9560
r_Cs=3.3950
r_I =2.3500
r_F =1.7500
r_IB=5.0000
r_CL=2.2230
r_FE=1.2000
r_LO=1.6000
r_LC=1.8500
r_RU=1.2000
r_H_ST=1.4870
r_C_ST=1.9080
r_N_ST=1.8240
r_O_ST=1.6612
r_XXST=1.9080

r_H =0.6000
r_HO=0.0000
r_HS=0.6000
r_HC=1.4870
r_H1=1.3870
r_H2=1.2870
r_H3=1.1870
r_HP=1.1000
r_HA=1.4590
r_H4=1.4090
r_H5=1.3590
r_HW=0.0000

r_O =1.6612
r_O2=1.6612
r_OW=1.7683
r_OH=1.7210
r_OS=1.6837

r_N =1.8240
r_N3=r_NA=r_N2=r_NStr=r_NC=r_NB=r_N3=r_NP=r_NO=r_N4=r_N5=r_N6=r_N7=r_N8=r_N

r_C =1.9080
r_CT=r_CM=r_C
r_CStr=r_CA=r_CB=r_CC=r_CP=r_CG=r_CN=r_CM=r_CK=r_CQ=r_CW=r_CV=r_CR=r_CX=r_CY=r_CD=r_C




r=ff.NewRes("BKBN","")
r.SetAtom("N"   , -0.41570, r_N)
r.SetAtom("H"   ,  0.27190, r_H)
r.SetAtom("CA"  ,  0.03370, r_CA)
r.SetAtom("HA"  ,  0.08230, r_H1)
r.SetAtom("C"   ,  0.59730, r_C)
r.SetAtom("O"   , -0.56790, r_O)

r=ff.NewRes("GLY","")
r.CopyAtomsFromRes(ff.GetRes("BKBN",""))
r.DelAtom("HA")
r.SetAtom("CA"  , -0.02520, 2.0)
r.SetAtom("HA2" ,  0.06980, 2.0)
r.SetAtom("HA3" ,  0.06980, 2.0)


r=ff.NewRes("ALA","")
r.CopyAtomsFromRes(ff.GetRes("BKBN",""))
r.SetAtom("CB"  , -0.18250, r_CB)
r.SetAtom("HB1" ,  0.06030, r_HC)
r.SetAtom("HB2" ,  0.06030, r_HC)
r.SetAtom("HB3" ,  0.06030, r_HC)

ff.PrintEntries()
print()


