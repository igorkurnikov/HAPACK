;Vij=((A/Rij)^N-1)kT:
[ ACE ]
  ;CH3  3.81 11.64  3.70 18.94
  CH3  3.81 11.31  3.72 18.18
    C  3.45 24.17  4.12 17.75
    O  2.51 36.48  4.49  8.53
[ NME ]
    N  4.33  7.77  3.18 33.21
 ; CH3  3.78 11.44  3.72 17.65
  CH3  3.81 11.31  3.72 18.18
[ GLY ]
; backbone average ["GLY","ALA","VAL","LEU","ILE"]
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
;    N  4.56  7.25  3.23 32.30
;   CA  4.40 11.43  3.84 23.00
;    C  3.50 20.30  4.08 21.25
;    O  2.50 46.66  4.48  7.65

;From first fit
;    N  4.36  8.01  3.16 44.69
;   CA  4.41  6.93  3.65 23.55
;    C  3.46 21.12  4.03 16.75
;    O  2.50 57.81  4.20 10.45
[ ALA ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
;    N  4.56  7.25  3.23 32.30
;   CA  4.40 11.43  3.84 23.00
;    C  3.50 20.30  4.08 21.25
;    O  2.50 46.66  4.48  7.65
;side chain "optimized"
   CB  3.83 10.69  3.77 19.19 ;from CH3 Avr ACE NME
;From first fit
;    N  4.66  7.01  3.20 36.41
;   CA  4.24 16.58  3.80 28.22
;   CB  3.86 12.16  3.73 19.00
;    C  3.51 18.25  4.03 20.34
;    O  2.50 29.26  4.64  6.69
[ VAL ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
;side chain "optimized"
;   CB  3.83 10.69  3.77 19.19 ;from CH3 Avr ACE NME
;   CB  5.00  5.98  4.30 10.58 ;from CB LEU, modifyed
   CB  4.47  9.99  4.12 13.67 ; VAL
  CG1  3.83 10.69  3.77 19.19 ;from CH3 Avr ACE NME
  CG2  3.83 10.69  3.77 19.19 ;from CH3 Avr ACE NME
;From first fit
;    N  4.51  8.45  3.34 21.28
;   CA  4.45 10.27  3.86 23.48
;   CB  4.47  9.99  4.12 13.67
;  CG1  3.61 11.73  3.85 14.73
;  CG2  3.94  7.74  3.72 37.44
;    C  3.57 16.05  4.09 23.21
;    O  2.49 34.11  4.70  5.96
[ LEU ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
;side chain "optimized"
;   CB  3.83 10.69  3.77 19.19 ;from CH3 Avr ACE NME
   CB  5.00  6.44  4.10  10.03 ;from CB LEU, modifyed
   CG  5.00  6.44  4.80  6.13 ;from CG LEU,modifyed
  CD1  3.83 10.69  3.77 19.19 ;from CH3 Avr ACE NME
  CD2  3.83 10.69  3.77 19.19 ;from CH3 Avr ACE NME
;From first fit
;    N  4.62  6.47  3.16 42.38
;   CA  4.43 12.24  3.88 22.13
;   CB  4.93  5.13  3.83 25.12
;   CG  4.88  5.59  4.32 12.23
;  CD1  3.73 11.89  3.95 12.60
;  CD2  3.77 11.36  3.92 13.32
;    C  3.50 18.44  4.11 22.28
;    O  2.49 30.08  4.46  8.99
[ ILE ] ;finish me
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
;   CA  4.57  9.16  3.65 23.55 ;from K-ILE Cl-GLY
;side chain "optimized"
   CB  5.00  6.44  4.10  10.03 ;from CB LEU, modifyed
  CG1  5.00  3.58  4.80  3.81 ;from hand fitting
  CG2  3.83 10.69  3.77 19.19 ;from CH3 Avr ACE NME
  CD1  3.83 10.69  3.77 19.19 ;from CH3 Avr ACE NME
;From first fit
;    N  4.66  6.34  3.24 35.80
;   CA  4.57  9.16  3.90 20.89
;   CB  4.59  7.77  3.93 26.11
;  CG1  4.25  7.39  3.91 15.90
;  CG2  3.53 17.35  3.87 14.25
;   CD  3.61 11.85  3.86 14.46
;    C  3.47 28.08  4.11 22.05
;    O  2.50 58.29  4.52  6.13
[ CYS ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
; backbone unique
   CA  4.41  6.93  3.78 13.84 ; CA-Cl==CB-Cl(w-fit)
;side chain "optimized"
;From first fit
;    N  4.63  7.11  3.22 33.00
;   CA  4.53  9.58  3.84 18.73
   CB  4.23  9.07  3.78 13.84
   SG  3.10 16.80  3.75 17.90
;    C  3.60 15.37  4.16 16.39
;    O  2.50 32.39  4.40  8.27
[ MET ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
;side chain "optimized"
;   CE  3.83 10.69  3.77 19.19 ;from CH3 Avr ACE NME, w-fit is slightly better
;From first fit
;    N  4.45  9.15  3.46 15.69
;   CA  4.54  9.63  3.84 24.74
   CB  4.42  7.38  3.88 22.14
   CG  4.06 12.57  3.80 20.40
   SD  3.01 21.49  4.20 14.37
   CE  4.01  7.69  3.84 12.61
;    C  3.51 17.94  4.04 38.21
;    O  2.50 56.46  4.43  8.20
[ SER ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
; backbone unique
;side chain "optimized"
;From first fit
;    N  4.24 11.58  3.19 38.51
;   CA  4.12 21.76  3.78 22.21
   CB  3.40 18.57  3.69 21.49
   OG  2.53 43.46  3.05 33.36
;    C  3.52 18.05  4.20 12.74
;    O  2.51 53.53  4.37  8.44
[ THR ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
; backbone unique
;side chain "optimized"
;From first fit
;    N  4.46 10.44  3.23 53.72
;   CA  4.48 10.19  3.71 26.46
   CB  3.46 15.91  3.76 23.33
  OG1  2.51 49.72  3.04 49.14
  CG2  3.55 13.32  3.77 17.56
;    C  3.54 22.94  3.94 37.32
;    O  2.47 32.10  4.21  8.78
[ PHE ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
; backbone unique
;side chain "optimized"
;   CG  3.71  5.68  4.38  7.98 ; avr over "CG","CD1","CD2","CE1","CE2","CZ", from w-fit-no11
;  CD1  3.71  5.68  4.38  7.98 ; avr over "CG","CD1","CD2","CE1","CE2","CZ", from w-fit-no11
;  CD2  3.71  5.68  4.38  7.98 ; avr over "CG","CD1","CD2","CE1","CE2","CZ", from w-fit-no11
;  CE1  3.71  5.68  4.38  7.98 ; avr over "CG","CD1","CD2","CE1","CE2","CZ", from w-fit-no11
;  CE2  3.71  5.68  4.38  7.98 ; avr over "CG","CD1","CD2","CE1","CE2","CZ", from w-fit-no11
;   CZ  3.71  5.68  4.38  7.98 ; avr over "CG","CD1","CD2","CE1","CE2","CZ", from w-fit-no11
   CG  3.18 19.40  4.17 14.87 ; avr over "CG","CD1","CD2","CE1","CE2","CZ", from w-fit-w11
  CD1  3.18 19.40  4.17 14.87 ; avr over "CG","CD1","CD2","CE1","CE2","CZ", from w-fit-w11
  CD2  3.18 19.40  4.17 14.87 ; avr over "CG","CD1","CD2","CE1","CE2","CZ", from w-fit-w11
  CE1  3.18 19.40  4.17 14.87 ; avr over "CG","CD1","CD2","CE1","CE2","CZ", from w-fit-w11
  CE2  3.18 19.40  4.17 14.87 ; avr over "CG","CD1","CD2","CE1","CE2","CZ", from w-fit-w11
   CZ  3.18 19.40  4.17 14.87 ; avr over "CG","CD1","CD2","CE1","CE2","CZ", from w-fit-w11
;From first fit
;    N  4.33 10.23  3.32 16.35
;   CA  4.33 12.05  3.86 17.40
   CB  4.02 13.25  3.64 33.64
;   CG  3.26 14.18  4.42 17.73
;  CD1  3.11 33.17  3.96 18.45
;  CD2  3.21 16.08  3.97 11.84
;  CE1  3.15 14.35  3.83 19.59
;  CE2  3.17 18.48  3.83 19.26
;   CZ  3.17 18.58  3.83 19.39
;    C  3.51 24.79  4.05 25.32
;    O  2.50 57.31  4.12 13.34
[ TYR ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
; backbone unique
;side chain "optimized"
;From first fit
;    N  4.47  8.94  3.25 47.91
;   CA  4.23 16.52  3.87 22.25
   CB  4.40  7.28  3.79 20.53
;   CG  3.43 27.71  3.83 22.89 ; CD1 CD2 Avr
   CG  4.54  4.85  4.76  8.58
;  CD1  3.43 21.91  3.82 25.91
;  CD2  3.44 31.02  3.84 18.59
;  CD1  3.43 27.71  3.83 22.89 ; CD1 CD2 Avr
;  CD2  3.43 27.71  3.83 22.89 ; CD1 CD2 Avr
  CD1  3.82 25.91  3.82 25.91 ; CD1-Cl All
  CD2  3.82 25.91  3.82 25.91 ; CD1-Cl All
;  CE1  3.22 21.27  3.83 15.79
;  CE2  3.34 21.04  3.70 27.76
  CE1  3.30 20.60  3.76 21.13 
  CE2  3.30 20.60  3.76 21.13
   CZ  3.05 27.06  3.84 25.56
   OH  2.53 43.22  2.99 36.38
;    C  3.52 23.39  4.08 23.28
;    O  2.48 29.30  4.68  6.51
[ TRP ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
; backbone unique
;side chain "optimized"
;From first fit
;    N  4.43  9.22  3.19 36.82
;   CA  4.48 11.75  3.81 26.75
   CB  4.57  5.87  3.76 22.55
   CG  3.42 23.07  4.45 14.09
  CD1  3.06 17.60  3.83 19.97
  CD2  3.24 29.62  4.76  9.62
  NE1  2.82 28.13  3.04 49.04
  CE2  3.05 43.57  3.88 31.54
  CE3  3.83  8.57  3.92 20.15
  CZ2  3.71  8.85  3.95 12.96
  CZ3  3.91  8.11  3.83 19.46
  CH2  3.92  7.99  3.94 12.89
;    C  3.61 14.67  4.11 21.91
;    O  2.49 32.93  4.35 10.18
[ HIS ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY ; fit is suitable but not very good but same as w-fit
; backbone unique
;side chain "optimized"
;From first fit
;    N  4.76  6.10  3.38 18.65
;   CA  4.36  9.10  3.90 16.54
   CB  4.11  7.49  3.99 14.22
   CG  5.00  3.74  4.11 22.50
  ND1  4.66  4.27  3.20 17.22
  CD2  3.56 16.57  3.55 22.41
  CE1  3.24 15.63  3.74 18.67
  NE2  2.64 26.45  4.05 13.42
;    C  3.61 12.37  4.13 20.96
;    O  2.54 24.04  4.56  8.10
[ ASN ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
; backbone unique
;side chain "optimized"
   CG  3.20 14.63  3.90 18.11
   CB  4.47  6.91  3.70 11.85
;From first fit
;    N  4.59  7.51  3.24 20.59
;   CA  4.59  8.30  3.87 23.19
;   CB  4.32  8.38  3.86 12.35
;   CG  3.36 15.75  4.07 19.63
  OD1  2.50 36.25  3.94 15.79
  ND2  3.47  9.01  3.20 18.46
;    C  3.44 32.48  4.06 25.44
;    O  2.48 30.13  4.68  6.13
[ GLN ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
; backbone unique
;side chain "optimized"
   CB  4.47  6.91  3.70 11.85 ; from CB ASN
   CG  4.20  6.71  3.80 12.20 
   CD  3.20 14.63  3.90 18.11 ; from CG ASN
  OE1  2.50 36.25  3.94 15.79 ; from OD1 ASN
  NE2  3.47  9.01  3.20 18.46 ; from ND2 ASN
;From first fit
;    N  4.52  8.68  3.25 19.60
;   CA  4.51  8.61  3.87 23.10
;   CB  4.47  6.91  3.82 15.86
;   CG  4.25  6.96  3.93 13.54
   CD  3.34 16.62  4.17 13.76
;  OE1  2.49 33.32  4.07 11.20
;  NE2  3.33 12.99  3.13 32.18
;    C  3.45 21.04  4.07 24.60
;    O  2.49 31.55  4.53  6.99
[ PRO ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
; backbone unique
;side chain "optimized"
;From first fit
;    N  4.48 12.05  4.93  6.72
;   CA  4.65  7.96  3.75 33.68
   CB  4.37  7.53  3.63 35.22
   CG  3.94 10.38  3.92 20.33
   CD  4.60  6.46  4.01 13.41
;    C  3.46 21.13  4.05 25.67
;    O  2.50 34.65  4.52 12.45
[ UGL ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
; backbone unique
;side chain "optimized"
  OE1  2.71 27.85  3.57 14.81 ; Avr OE1,OE2
  OE2  2.71 27.85  3.57 14.81 ; Avr OE1,OE2
;From first fit
;    N  4.52  7.68  3.18 38.45
   CA  4.74  7.24  4.07 15.08
   CB  4.67  5.66  3.64 34.11
   CG  4.22 13.92  3.83 15.37
   CD  3.91 11.49  3.86 17.91
;  OE1  2.72 27.19  3.51 18.15
;  OE2  2.70 29.61  3.65 11.46
;    C  3.59 12.46  4.21 14.47
;    O  2.50 36.06  4.41  8.36
[ GLU ]
; backbone "optimized"
;    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
;   CA  4.41  6.93  3.65 23.55 ;from GLY
; backbone unique
    N  3.60  9.90  4.70 7.03 ; m-fit
   CA  3.44  9.60  4.45 15.63 ; CD
;side chain "optimized"
   CB  3.44  9.60  4.45 15.63 ; CD
   CG  3.44  9.60  4.45 15.63 ; CD
   CD  3.10 16.23  4.45 15.63
  OE1  2.58 18.78  4.22 25.13 ; Avr OE1,OE2
  OE2  2.58 18.78  4.22 25.13 ; Avr OE1,OE2
;From first fit
;    N  4.50  6.54  5.00  4.08
;   CA  4.47  7.44  4.14 25.20
;   CB  4.68  5.28  3.76 30.60
;   CG  4.64  5.29  3.81 36.44
;   CD  3.44  9.60  4.45 15.63
;  OE1  2.59 17.79  4.16 34.12
;  OE2  2.57 20.05  4.32 16.43
;    C  3.74  9.75  4.93  8.44
;    O  2.53 26.08  4.79  6.12
[ ASP ]
; backbone "optimized"
;    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
;   CA  4.41  6.93  3.65 23.55 ;from GLY
; backbone unique
    N  3.60  9.90  4.70 7.03 ;from GLU m-fit
   CA  3.44  9.60  4.45 15.63 ;from GLU CD
;side chain "optimized"
   CB  3.44  9.60  4.45 15.63 ; from GLU CD
   CG  3.10 16.23  4.45 15.63 ; from GLU
  OD1  2.58 18.78  4.22 25.13 ; from GLU Avr OE1,OE2
  OD2  2.58 18.78  4.22 25.13 ; from GLU Avr OE1,OE2
;From first fit
;    N  4.72  5.92  7.63  2.37
;   CA  5.00 12.00  6.78  3.11
   CB  4.46  6.27  5.52  3.59
   CG  3.55  9.09  4.50 14.75
  OD1  2.61 16.73  4.44  8.73
  OD2  2.62 16.37  4.23 13.28
;    C  3.72 10.08  6.33  4.56
;    O  2.67 16.46  5.58  4.76
[ ARG ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
;    O  2.51 36.48  4.49  8.53 ;from Avr ACE
;   CA  4.41  6.93  3.65 23.55 ;from GLY
; backbone unique
    O  2.51 36.48  4.10  8.89 
   CA  4.00 15.29  3.65 23.55 
;side chain "optimized"
   CB  4.81  6.43  3.60 18.61 ;CD
   CG  4.81  6.43  3.60 18.61 ;CD
   CD  4.81  6.43  3.60 18.61
   NE  4.00  9.75  3.10 23.99
   CZ  4.70  6.46  3.55 23.94
  NH1  4.00  9.75  3.10 23.99
  NH2  4.00  9.75  3.10 23.99
;From first fit
;    N  4.47 10.36  3.23 22.85
;   CA  4.90  6.50  3.83 20.67
;   CB  4.75  8.21  4.05 10.52
;   CG  4.81  6.43  4.03 11.06
;   CD  4.66  4.62  4.04 10.64
;   NE  4.76  3.20  3.21 19.70
;   CZ  4.69  4.00  3.88 14.69
;  NH1  4.65  4.85  3.13 20.29
;  NH2  4.95  3.32  3.08 27.76
;    C  3.62 14.92  4.19 13.86
;    O  2.52 28.78  4.37  8.23
[ LYS ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
; backbone unique
;side chain "optimized"
   CB  4.86  5.93  3.60 13.63 ; from CE
   CG  4.86  5.93  3.70 14.04 ; from CE,mod
   CD  4.86  5.93  3.60 13.63 ; from CE
   CE  4.86  5.93  3.60 13.63
;From first fit
;    N  4.42  9.84  3.21 34.52
;   CA  4.43 10.91  3.90 22.24
;   CB  4.81  7.79  4.20  9.10
;   CG  4.95  6.85  4.19 10.39
;   CD  5.00  5.07  4.11 10.05
;   CE  4.46  8.45  3.85 11.68
   NZ  3.96 17.31  3.10 24.74
;    C  3.53 18.08  4.23 12.72
;    O  2.49 40.10  4.54  7.09

[ ULY ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
; backbone unique
;side chain "optimized"
;From first fit
;    N  4.62  6.57  3.20 35.45
;   CA  4.45 10.21  3.83 25.58
   CB  4.60  5.77  3.71 39.36
   CG  4.84  5.16  4.34 10.13
   CD  4.46  5.59  4.04 19.73
   CE  3.61 11.96  3.94 19.49
   NZ  2.96 17.98  3.46 21.18
;    C  3.45 30.87  4.12 21.81
;    O  2.49 34.48  4.50  7.80
[ UAR ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
    C  3.45 24.17  4.12 17.75 ;from Avr ACE
    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
; backbone unique
;side chain "optimized"
  NH1  2.76 24.13  3.39 17.47
  NH2  2.76 24.13  3.39 17.47
;From first fit
;    N  4.71  7.45  3.16 44.36
;   CA  4.56 10.09  3.82 20.03
   CB  4.61  7.61  3.87 14.90
   CG  4.01 13.29  3.85 24.70
   CD  3.45 21.17  4.34  8.57
   NE  3.04 14.50  3.21 35.64
   CZ  2.75 17.12  3.94 26.53
;  NH1  2.77 23.45  3.38 19.76
;  NH2  2.75 25.42  3.41 14.57
;    C  3.55 21.86  4.04 20.55
;    O  2.51 26.50  4.41  7.60
[ NTERMALA ]
; backbone "optimized"
    N  4.33  7.77  3.18 33.21 ;from Avr NME
;    C  3.45 24.17  4.12 17.75 ;from Avr ACE
;    O  2.51 36.48  4.49  8.53 ;from Avr ACE
   CA  4.41  6.93  3.65 23.55 ;from GLY
; backbone unique
    N  3.96 17.31  3.10 24.74 ;NZ of LYS
   CA  4.86  5.93  3.60 13.63 ;CE of LYS
    C  3.45 24.17  3.70 15.44 ;from Avr ACE
    O  2.51 36.48  3.50 13.21 ;from Avr ACE
;side chain "optimized"
   CB  3.83 10.69  3.77 19.19 ;from CH3 Avr ACE NME
;From first fit
;    N  4.02 39.37  3.13 22.17
;   CA  4.59 19.30  3.87 13.48
;   CB  4.58  7.55  3.93 11.05
;    C  3.50 35.99  4.24  9.17
;    O  2.48 27.97  4.06  8.24
[ CTERMALA ]
; backbone "optimized"
;    N  4.33  7.77  3.18 33.21 ;from Avr NME
;    C  3.45 24.17  4.12 17.75 ;from Avr ACE
;    O  2.51 36.48  4.49  8.53 ;from Avr ACE
;  OXT  2.51 36.48  4.49  8.53 ;from Avr ACE
;   CA  4.41  6.93  3.65 23.55 ;from GLY
; backbone unique
    N  3.60  9.90  4.70 7.03 ;from GLU m-fit
   CA  3.44  9.60  4.45 15.63 ;from GLU CD
;    O  2.60 17.81  4.19 19.36 ; Avr O OXT
;  OXT  2.60 17.81  4.19 19.36 ; Avr O OXT
    C  3.10 16.23  4.45 15.63 ; from GLU CG
    O  2.58 18.78  4.22 25.13 ; from GLU Avr OE1,OE2
  OXT  2.58 18.78  4.22 25.13 ; from GLU Avr OE1,OE2
;side chain "optimized"
   CB  3.83 10.69  3.77 19.19 ;from CH3 Avr ACE NME
;From first fit
;    N  4.69  5.17  5.00  3.70
;   CA  4.40  7.12  4.51 11.89
;   CB  4.26  5.91  4.01 13.25
;    C  3.56  8.96  4.55 13.32
;   O1  2.59 19.16  4.35 10.93
;   O2  2.61 16.97  4.11 27.84
[ K ]
     K  4.08 12.54  2.99 24.44
[ CL ]
    CL  2.99 24.44  4.08 12.54
[ K+ ]
     K+  4.08 12.54  2.99 24.44
[ CL- ]
    CL  2.99 24.44  4.08 12.54
