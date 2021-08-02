        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:20 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SETOUT__genmod
          INTERFACE 
            SUBROUTINE SETOUT(XX,MGRID,PHIMAP,PHIMAP1,PHIMAP2,PHIMAP3,  &
     &ATMCRG,IEPSMP,IEPSMP2,IDEBMAP,IOFF,XN2,RAD3,NATOM,EXRAD,RADPRB)
              INTEGER(KIND=4) :: MGRID
              REAL(KIND=4) :: XX(*)
              REAL(KIND=4) :: PHIMAP(MGRID,MGRID,*)
              REAL(KIND=4) :: PHIMAP1(*)
              REAL(KIND=4) :: PHIMAP2(*)
              REAL(KIND=4) :: PHIMAP3(*)
              REAL(KIND=4) :: ATMCRG(4,*)
              INTEGER(KIND=4) :: IEPSMP(MGRID,MGRID,MGRID,3)
              INTEGER(KIND=4) :: IEPSMP2(MGRID,MGRID,MGRID,3)
              INTEGER(KIND=4) :: IDEBMAP(MGRID,MGRID,MGRID)
              INTEGER(KIND=4) :: IOFF(3,*)
              REAL(KIND=4) :: XN2(3,50000)
              REAL(KIND=4) :: RAD3(50000)
              INTEGER(KIND=4) :: NATOM
              REAL(KIND=4) :: EXRAD
              REAL(KIND=4) :: RADPRB
            END SUBROUTINE SETOUT
          END INTERFACE 
        END MODULE SETOUT__genmod
