def VwEpsilon0(filename):
    print("Load NodeIndexing from " + filename)
    NIndx = NodeIndexing()
    NIndx.ReadFromFile(filename)
    f3d = NIndx.GetHaField3D(NodeIndexing.DielConst, NodeIndexing.Epsilon0)

    Canv = GetmVueGlCanvas()
    PlaneV = Canv.CreatePlaneViewer()

    PlaneV.SetHaField3D(f3d)


def VwEpsilon1(filename):
    print("Load NodeIndexing from " + filename)
    NIndx = NodeIndexing()
    NIndx.ReadFromFile(filename)
    f3d = NIndx.GetHaField3D(NodeIndexing.DielConst, NodeIndexing.Epsilon1)

    Canv = GetmVueGlCanvas()
    PlaneV = Canv.CreatePlaneViewer()

    PlaneV.SetHaField3D(f3d)


def VwEpsilon2(filename):
    print("Load NodeIndexing from " + filename)
    NIndx = NodeIndexing()
    NIndx.ReadFromFile(filename)
    f3d = NIndx.GetHaField3D(NodeIndexing.DielConst, NodeIndexing.Epsilon2)

    Canv = GetmVueGlCanvas()
    PlaneV = Canv.CreatePlaneViewer()

    PlaneV.SetHaField3D(f3d)


def VwDiff0(filename):
    print("Load NodeIndexing from " + filename)
    NIndx = NodeIndexing()
    NIndx.ReadFromFile(filename)
    f3d = NIndx.GetHaField3D(NodeIndexing.DiffConst, NodeIndexing.Ion0)

    Canv = GetmVueGlCanvas()
    PlaneV = Canv.CreatePlaneViewer()

    PlaneV.SetHaField3D(f3d)


def VwDiff1(filename):
    print("Load NodeIndexing from " + filename)
    NIndx = NodeIndexing()
    NIndx.ReadFromFile(filename)
    f3d = NIndx.GetHaField3D(NodeIndexing.DiffConst, NodeIndexing.Ion1)

    Canv = GetmVueGlCanvas()
    PlaneV = Canv.CreatePlaneViewer()

    PlaneV.SetHaField3D(f3d)


def VwQ(filename):
    print("Load NodeIndexing from " + filename)
    NIndx = NodeIndexing()
    NIndx.ReadFromFile(filename)
    f3d = NIndx.GetHaField3D(NodeIndexing.Charge, NodeIndexing.ChargeMask)

    Canv = GetmVueGlCanvas()
    PlaneV = Canv.CreatePlaneViewer()

    PlaneV.SetHaField3D(f3d)


def CompQ(filename1, filename2):
    print("Load NodeIndexing from " + filename1)
    NIndx = NodeIndexing()
    NIndx.ReadFromFile(filename1)
    f3d1 = NIndx.GetHaField3D(NodeIndexing.Charge, NodeIndexing.ChargeMask)

    print("Load VectorField from " + filename2)
    VF = VectorField3D(filename2)

    f3d2 = VF.GetHaField3D(0)
    Nx = f3d1.GetNx()
    print("Nx" + `Nx`)
    rmsd = f3d1.CompareHaField3D(f3d2, 0.1)


def VwInitConc0(filename):
    print("Load NodeIndexing from " + filename)
    NIndx = NodeIndexing()
    NIndx.ReadFromFile(filename)
    f3d = NIndx.GetHaField3D(NodeIndexing.Conc, NodeIndexing.Ion0)

    Canv = GetmVueGlCanvas()
    PlaneV = Canv.CreatePlaneViewer()

    PlaneV.SetHaField3D(f3d)


def VwInitConc1(filename):
    print("Load NodeIndexing from " + filename)
    NIndx = NodeIndexing()
    NIndx.ReadFromFile(filename)
    f3d = NIndx.GetHaField3D(NodeIndexing.Conc, NodeIndexing.Ion1)

    Canv = GetmVueGlCanvas()
    PlaneV = Canv.CreatePlaneViewer()

    PlaneV.SetHaField3D(f3d)


def VwPot(filename, scale):
    print("Load Potential from " + filename)
    # f3d=HaField3D()
    # f3d.LoadGZ(filename)
    Canv = GetmVueGlCanvas()
    PlaneV = Canv.CreatePlaneViewer()
    PlaneV.SetHaField3DFromFileWithScale(filename, scale)
    # PlaneV.SetHaField3D(f3d)


def VwVectorField3D(filename, nelem):
    print("Load Concentration0 from " + filename)
    Conc = VectorField3D(filename)
    C0 = Conc.GetHaField3D(nelem)
    # Conc.__del__
    # f3d=HaField3D()
    # f3d.LoadGZ(filename)
    Canv = GetmVueGlCanvas()
    PlaneV = Canv.CreatePlaneViewer()
    PlaneV.SetHaField3D(C0)
    # PlaneV.SetHaField3D(f3d)


def VwVectorField3D2(filename):
    print("Load VectorField from " + filename)
    Conc = VectorField3D(filename)

    C0 = Conc.GetHaField3D(0)
    C0.SetName("VF[" + `0` + "]")
    Canv = GetmVueGlCanvas()
    PlaneV0 = Canv.CreatePlaneViewer()
    PlaneV0.SetHaField3D(C0)

    C1 = Conc.GetHaField3D(1)
    C1.SetName("VF[" + `1` + "]")
    Canv = GetmVueGlCanvas()
    PlaneV1 = Canv.CreatePlaneViewer()
    PlaneV1.SetHaField3D(C1)


def VwVectorField3D(filename):
    print("Load VectorField from " + filename)
    Conc = VectorField3D(filename)

    for x in range(0, Conc.Nelem, 1):
        C0 = Conc.GetHaField3D(x)
        # C0.SetName("VF["+`x`+"]")
        Canv = GetmVueGlCanvas()
        PlaneV0 = Canv.CreatePlaneViewer()
        PlaneV0.SetHaField3D(C0)
