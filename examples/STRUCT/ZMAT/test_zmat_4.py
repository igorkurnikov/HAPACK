mset = GetCurMolSet()
zm = mset.GetZMat()
zmat_str = """
Mg
O     1     rMO1
H     2     0.95          1     aMOH1
H     2     rOH12        3     aHOH1        1     dMHOH1
O     1     rMO2         2     aOMO2        3     dOMOH2
H     5     rOH21        1     aMOH2        2     dHOMO2
H     5     rOH22        6     aHOH2        1     dMHOH2
O     1     rMO3         2     aOMO3        3     dOMOH3
H     8     rOH31        1     aMOH3        2     dHOMO3
H     8     rOH32        9     aHOH3        1     dMHOH3
O     1     rMO4         5     aOMO4        6     dOMOH4
H     11     rOH41        1     aMOH4        5     dHOMO4
H     11     rOH42        12     aHOH4        1     dMHOH4
O     1     rMO5         2     aOMO5        3     dOMOH5
H     14     rOH51        1     aMOH5        2     dHOMO5
H     14     rOH52        15     aHOH5        1     dMHOH5
O     1     rMO6         2     aOMO6        3     dOMOH6
H     17     rOH61        1     aMOH6        2     dHOMO6
H     17     rOH62        18     aHOH6        1     dMHOH6
O     1     rMO7     2     aOMO7        3     dOMOH7
H     20     rOH71        1     aMOH7        2     dHOMO7
H     20     rOH72        21     aHOH7        1     dMHOH7

rMO1   2.11247309
rOH12    0.95341086
aMOH1   126.4789322
aHOH1   107.4375181
dMHOH1    -176.1937715
rMO2    2.11993303
rOH21    0.95291499
rOH22    0.95331212
aMOH2    126.0634331
aHOH2    106.9683016
dMHOH2   -167.3328815
aOMO2    86.53210581
dOMOH2   -175.2259583
dHOMO2    90.3834173
rMO3    2.09078869
rOH31    0.95248506
rOH32    0.96159493
aMOH3    128.8425032
aHOH3    108.2715025
dMHOH3   176.3897237
aOMO3    93.03391603
dOMOH3   -85.98541282
dHOMO3   3.92765029
rMO4   2.09467801
rOH41    0.95212312
rOH42    0.96103362
aMOH4   128.412064
aHOH4   107.0558795
dMHOH4   157.0655298
aOMO4    94.16458312
dOMOH4   -89.65408357
dHOMO4   -134.3273927
rMO5   2.12255818
rOH51    0.95345963
rOH52    0.95341795
aMOH5   124.8395297
aHOH5   106.9748438
dMHOH5   168.8434896
aOMO5    85.86304319
dOMOH5   5.80606768
dHOMO5   -103.4893761
rOH61   0.95391008
rOH62    0.95369514
aMOH6   123.9207953
aHOH6    106.8334524
dMHOH6   -174.5228904
aOMO6    89.62591725
dOMOH6   95.78766795
dHOMO6   -2.2583911
rOH71    0.9538636
rOH72    0.95369054
aMOH7    129.5063823
aHOH7    105.2902395
dMHOH7    -178.8382717
aOMO7    135.9682558
dOMOH7    -95.28277733
dHOMO7    99.61542168
rMO7   3.8

rMO6    2.2
"""
zm.LoadFromString(zmat_str)
aptr = mset.GetAtomByRef("HOH2.O")
elc = zm.GetRCrdForAtom(aptr)
elc.SetFrozen()
elc.SetValue(2.5)
aptr = mset.GetAtomByRef("HOH3.O")
elc = zm.GetRCrdForAtom(aptr)
elc.SetFrozen()
elc.SetValue(3.0)
print(zm.SaveToString())
