db = HaMatDB()
db.open("ndo_pars_1.db","w")
db.AddFromFile("ndo_pars_1.dat")
db.ExtractAllToFile("out.dat")
db.close()


