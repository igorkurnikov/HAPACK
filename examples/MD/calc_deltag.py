#
# This example script to run thermodynamic Integration 
#
pmset = GetCurMolSet()
pmset.SetName("mg_wat_1_ti")
mm_mod = pmset.GetMolMechMod(1)
ti_mod = mm_mod.GetTISimMod()
ti_mod.SetNumLambda(3)
ti_mod.SetNumEqPoints(2000)
delta_g = ti_mod.CalcDeltaG()
ti_mod.ReduceDvDlData(1000)
#print "Delta G = ",delta_g


