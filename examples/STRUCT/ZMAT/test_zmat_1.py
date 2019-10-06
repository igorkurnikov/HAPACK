mset = GetCurMolSet()
zm = mset.GetZMat()
zm.InitStdZMat()
print(zm.SaveToString())
