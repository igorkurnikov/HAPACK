        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:20 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MKEPS__genmod
          INTERFACE 
            SUBROUTINE MKEPS(XX,MGRID,PHIMAP,PHIMAP1,PHIMAP2,PHIMAP3,   &
     &ATMCRG,IEPSMP,IEPSMP2,IDEBMAP,RADPRB,EPSOUT,IBNUM,IBGRD,RAD3,XN2)
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
              REAL(KIND=4) :: RADPRB
              REAL(KIND=4) :: EPSOUT
              INTEGER(KIND=4) :: IBNUM
              INTEGER(KIND=4) :: IBGRD(4,*)
              REAL(KIND=4) :: RAD3(50000)
              REAL(KIND=4) :: XN2(3,50000)
            END SUBROUTINE MKEPS
          END INTERFACE 
        END MODULE MKEPS__genmod
