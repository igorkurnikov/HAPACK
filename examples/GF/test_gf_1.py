a = HaMat_double()
b = HaMat_double()
a.newsize(2,2)
b.newsize(2,1)
a.SetVal(1,1,0.4)
a.SetVal(1,2,0.5)
a.SetVal(2,1,0.2)
a.SetVal(2,2,0.3)
b.SetVal(1,1,1.0)
b.SetVal(2,1,1.0)
ir = HaMat_double_solv_lin_syst_1(a,b)
print "result = ",ir
print a.GetVal(1,1)
print b.GetVal(1,1)
print b.GetVal(2,1)
# Output: 
#
#  result =  1
#  0.4
#  -10.0
#   10.0

