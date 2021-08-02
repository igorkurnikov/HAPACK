        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:20 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHRGUP__genmod
          INTERFACE 
            SUBROUTINE CHRGUP(XX,MGRID,PHIMAP,PHIMAP1,PHIMAP2,PHIMAP3,  &
     &ATMCRG,IEPSMP,IEPSMP2,IDEBMAP,NQGRD,CHRGV2,QVAL,IQPOS,FPOH,SIXEPS,&
     &DIFEPS,DEBFCT,GCHRGP,GCHRG,ICOUNT1A,ICOUNT1B)
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
              INTEGER(KIND=4) :: NQGRD
              REAL(KIND=4) :: CHRGV2(4,50000)
              REAL(KIND=4) :: QVAL(200000)
              INTEGER(KIND=4) :: IQPOS(200000)
              REAL(KIND=4) :: FPOH
              REAL(KIND=4) :: SIXEPS
              REAL(KIND=4) :: DIFEPS
              REAL(KIND=4) :: DEBFCT
              INTEGER(KIND=4) :: GCHRGP(3,200000)
              REAL(KIND=4) :: GCHRG(200000)
              INTEGER(KIND=4) :: ICOUNT1A
              INTEGER(KIND=4) :: ICOUNT1B
            END SUBROUTINE CHRGUP
          END INTERFACE 
        END MODULE CHRGUP__genmod
