db = HaMatDblDB()
db.open("test1.db","w")
db.AddFromFile("test_mat.dat")
db.ExtractAllToFile("out.dat")

