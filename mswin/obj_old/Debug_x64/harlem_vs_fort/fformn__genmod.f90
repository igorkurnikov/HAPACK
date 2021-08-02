        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:54 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FFORMN__genmod
          INTERFACE 
            SUBROUTINE FFORMN(IOUT,IPRINT,N,NATOMS,IAN,IATTYP,IU,LLIM,H,&
     &P1,P2,GSS,GSD,GDD,F1,ENERGY,G1,F2)
              INTEGER(KIND=4) :: IOUT
              INTEGER(KIND=4) :: IPRINT
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NATOMS
              INTEGER(KIND=4) :: IAN(*)
              INTEGER(KIND=4) :: IATTYP(*)
              INTEGER(KIND=4) :: IU(*)
              INTEGER(KIND=4) :: LLIM(*)
              REAL(KIND=8) :: H(*)
              REAL(KIND=8) :: P1(*)
              REAL(KIND=8) :: P2(*)
              REAL(KIND=8) :: GSS(*)
              REAL(KIND=8) :: GSD(*)
              REAL(KIND=8) :: GDD(*)
              REAL(KIND=8) :: F1(*)
              REAL(KIND=8) :: ENERGY
              REAL(KIND=8) :: G1(18)
              REAL(KIND=8) :: F2(18)
            END SUBROUTINE FFORMN
          END INTERFACE 
        END MODULE FFORMN__genmod
