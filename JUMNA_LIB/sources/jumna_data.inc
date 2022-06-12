      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      
c...Maximum number of rings/subunits per ligand
      parameter (n0=5)

c...Maximum number of atoms in the DNA oligomer
      parameter (n1=4000)

c...Maximum number of nucleotides
      parameter (n2=120)

c...Maximum number of additional distance constraints
      parameter (n3=300)

c...Maximum number of atoms per nucleotide
      parameter (n4=230)

c...Maximum number of types of nucleotide
      parameter (n5=20)

c...Maximum number of physical backbone variables
      parameter (n6=12000)

c...Maximum number of physical variables (lock only)
      parameter (n6a=12000)

c...Maximum number of independent, unlocked variables
      parameter (n7=4400)

c...Maximum number of physical variables per nucleotide
      parameter (n8=200)

c...Maximum number of non-bonded ligands
      parameter (n9=20)

c...Delta change of variables for numerical gradient test 
      parameter (deltv=1.d-7)

c...Conversion factor for Kcal/mole 
      parameter (convk=332.053986d0)

c...Conversion factor from degrees to radians 
      parameter (cdr=0.017453292519943295769d0)	

c...Conversion factor from radians to degrees 
      parameter (crd=57.295779513082320877d0)

c...Constant: Pi
      parameter (pi=3.14159265358979323846d0)

c...Limit distance for H-bonds (upper limit, angstroms)
      parameter (rhblim=3.d0)

c...Limit angle for H-bonds (upper limit, degrees)
      parameter (anglim=40.d0)

c...Limit energy for steric contacts (lower limit, kcal/mol)
      parameter (conlim=5.d0)
