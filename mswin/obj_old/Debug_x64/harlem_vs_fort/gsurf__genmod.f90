        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:53 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GSURF__genmod
          INTERFACE 
            SUBROUTINE GSURF(KSURF,RMIN,OFAC,RD,NDIV,ASS1,NATOM,GHOST,NP&
     &,VOL)
              INTEGER(KIND=4) :: KSURF
              REAL(KIND=4) :: RMIN
              REAL(KIND=4) :: OFAC
              REAL(KIND=4) :: RD
              INTEGER(KIND=4) :: NDIV
              LOGICAL(KIND=4) :: ASS1
              INTEGER(KIND=4) :: NATOM
              LOGICAL(KIND=4) :: GHOST
              INTEGER(KIND=4) :: NP
              REAL(KIND=8) :: VOL
            END SUBROUTINE GSURF
          END INTERFACE 
        END MODULE GSURF__genmod
