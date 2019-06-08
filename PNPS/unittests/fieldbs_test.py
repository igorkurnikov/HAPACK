import math
import pnpsll
from pnps import new_boolArray, boolArray_setitem, boolArray_getitem
from pnps import new_intArray, intArray_setitem, intArray_getitem
from pnps import new_floatArray, floatArray_setitem, floatArray_getitem


def gen_FieldBW(GS_X, GS_Y, GS_Z, Pot="Gen"):
    """
    Generate FieldBW instance.
    If Pot is "Gen", init values as ix+iy*10+iz*100
    return tuple(FieldBW,Potential)
    """
    potbw = pnpsll.pnpsll.FieldBW()

    # calculate strides
    GS_XY = GS_X * GS_Y
    GS_XYZ = GS_X * GS_Y * GS_Z

    # set grid size
    GridSize = new_intArray(3)
    intArray_setitem(GridSize, 0, GS_X)
    intArray_setitem(GridSize, 1, GS_Y)
    intArray_setitem(GridSize, 2, GS_Z)

    potbw.Init(GridSize)

    # initiate potential field
    PotIn = new_floatArray(GS_XYZ)
    for iz in range(GS_Z):
        for iy in range(GS_Y):
            for ix in range(GS_X):
                grd_pnt = ix + iy * GS_X + iz * GS_XY
                floatArray_setitem(PotIn, grd_pnt, ix + iy * 10 + iz * 100)

    potbw.SetFromField(PotIn)

    return potbw, PotIn


def test_FieldBW__check_field_setter_and_getter():
    """
    Test FieldBW::SetFromField and FieldBW::SetField
    """
    GS_X = 10
    GS_Y = 8
    GS_Z = 6
    GS_XYZ = GS_X * GS_Y * GS_Z

    potbw, PotIn = gen_FieldBW(GS_X, GS_Y, GS_Z, Pot="Gen")

    PotOut = new_floatArray(GS_XYZ)
    for grd_pnt in range(GS_XYZ):
        floatArray_setitem(PotOut, grd_pnt, -1)
    potbw.SetField(PotOut)

    sd = 0.0
    for grd_pnt in range(GS_XYZ):
        sd += (floatArray_getitem(PotOut, grd_pnt) - floatArray_getitem(PotIn, grd_pnt)) ** 2
    sd = math.sqrt(sd / GS_XYZ)
    # print "SD:", sd
    assert sd < 0.1

######################################################################################################################
# FieldBW::BorderExchange


def BorderExchangeTest(GS_X, GS_Y, GS_Z, PBC_X, PBC_Y, PBC_Z):
    GS_XY = GS_X * GS_Y
    GS_XYZ = GS_X * GS_Y * GS_Z

    pbc = (PBC_X, PBC_Y, PBC_Z)
    PBC = new_boolArray(3)
    boolArray_setitem(PBC, 0, pbc[0])
    boolArray_setitem(PBC, 1, pbc[1])
    boolArray_setitem(PBC, 2, pbc[2])

    potbw, PotIn = gen_FieldBW(GS_X, GS_Y, GS_Z, Pot="Gen")

    PotOut = new_floatArray(GS_XYZ)
    for grd_pnt in range(GS_XYZ):
        floatArray_setitem(PotOut, grd_pnt, -1)

    potbw.BorderExchange(PBC)

    potbw.SetField(PotOut)

    #         B  W  B   W   B  W
    #         W  B  W   B   W  B
    #         B  W  B   W   B  W
    #         W  B  W   B   W  B
    # index   0  0     LL   L   L
    # ix      0  1         G-2 G-1
    if pbc[0]:
        for iz in range(GS_Z):
            for iy in range(GS_Y):
                grd_pnt_0 = iy * GS_X + iz * GS_XY
                grd_pnt_L = GS_X - 1 + iy * GS_X + iz * GS_XY

                print(floatArray_getitem(PotOut, grd_pnt_0) , GS_X - 2 + iy * 10 + iz * 100)

                assert floatArray_getitem(PotOut, grd_pnt_0) == GS_X - 2 + iy * 10 + iz * 100
                assert floatArray_getitem(PotOut, grd_pnt_0 + 1) == 1 + iy * 10 + iz * 100
                assert floatArray_getitem(PotOut, grd_pnt_L) == 1 + iy * 10 + iz * 100
                assert floatArray_getitem(PotOut, grd_pnt_L - 1) == GS_X - 2 + iy * 10 + iz * 100
    if pbc[1]:
        for iz in range(GS_Z):
            for ix in range(GS_X):
                grd_pnt_0 = ix + iz * GS_XY
                grd_pnt_L = ix + (GS_Y - 1) * GS_X + iz * GS_XY



                assert floatArray_getitem(PotOut, grd_pnt_0) == ix + (GS_Y - 2) * 10 + iz * 100
                assert floatArray_getitem(PotOut, grd_pnt_0 + GS_X) == ix + 1 * 10 + iz * 100
                assert floatArray_getitem(PotOut, grd_pnt_L) == ix + 1 * 10 + iz * 100
                assert floatArray_getitem(PotOut, grd_pnt_L - GS_X) == ix + (GS_Y - 2) * 10 + iz * 100
    if pbc[2]:
        for iy in range(GS_Y):
            for ix in range(GS_X):
                grd_pnt_0 = ix + iy * GS_X
                grd_pnt_L = ix + iy * GS_X + (GS_Z - 1) * GS_XY

                assert floatArray_getitem(PotOut, grd_pnt_0) == ix + iy * 10 + (GS_Z - 2) * 100
                assert floatArray_getitem(PotOut, grd_pnt_0 + GS_XY) == ix + iy * 10 + 1 * 100
                assert floatArray_getitem(PotOut, grd_pnt_L) == ix + iy * 10 + 1 * 100
                assert floatArray_getitem(PotOut, grd_pnt_L - GS_XY) == ix + iy * 10 + (GS_Z - 2) * 100


def test_FieldBW__BorderExchange():
    """
    Test FieldBW::SetFromField and FieldBW::SetField
    """
    test_set = [
        (6, 8, 10, True, False, False),
        (7, 8, 10, True, False, False),
        (6, 7, 10, True, False, False),
        (6, 7, 11, True, False, False),
        (5, 7, 11, True, False, False),
        (6, 8, 10, False, True, False),
        (6, 7, 10, False, True, False),
        (6, 8, 11, False, True, False),
        (6, 7, 11, False, True, False),
        (5, 8, 10, False, True, False),
        (5, 7, 10, False, True, False),
        (6, 8, 10, False, False, True),
        (6, 8, 11, False, False, True),
        (6, 7, 10, False, False, True),
        (6, 7, 11, False, False, True),
        (5, 7, 10, False, False, True),
        (5, 7, 11, False, False, True)
    ]
    for t in test_set:
        BorderExchangeTest(*t)


if __name__ == "__main__":
    test_FieldBW__check_field_setter_and_getter()
    test_FieldBW__BorderExchange()
