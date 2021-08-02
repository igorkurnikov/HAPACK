        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:20 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SETBC__genmod
          INTERFACE 
            SUBROUTINE SETBC(XX,MGRID,PHIMAP,PHIMAP1,PHIMAP2,PHIMAP3,   &
     &ATMCRG,IEPSMP,IEPSMP2,IDEBMAP,IBCTYP,IPER,QPLUS,QMIN,CQPLUS,CQMIN,&
     &EPSOUT,DEBLEN,NQASS)
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
              INTEGER(KIND=4) :: IBCTYP
              LOGICAL(KIND=4) :: IPER(3)
              REAL(KIND=4) :: QPLUS
              REAL(KIND=4) :: QMIN
              REAL(KIND=4) :: CQPLUS(3)
              REAL(KIND=4) :: CQMIN(3)
              REAL(KIND=4) :: EPSOUT
              REAL(KIND=4) :: DEBLEN
              INTEGER(KIND=4) :: NQASS
            END SUBROUTINE SETBC
          END INTERFACE 
        END MODULE SETBC__genmod
