#bibl_db = BiblDB()
#bibl_db.Init()
for i in range(32041,32270):
  ref = bibl_db.GetRefByID(i)
  bibl_db.PrintRef1(ref) 
  bibl_db.LoadRefFromWeb(ref)



