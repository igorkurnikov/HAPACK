        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:20 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DELPHI__genmod
          INTERFACE 
            SUBROUTINE DELPHI(XX,MGRID,PHIMAP,PHIMAP1,PHIMAP2,PHIMAP3,  &
     &ATMCRG,IEPSMP,IEPSMP2,IDEBMAP,SF1,SF2,QMAP1,QMAP2,DEBMAP1,DEBMAP2,&
     &BNDX1,BNDX2,BNDX3,BNDX4,BNDX,BNDY,BNDZ,DB,IBGRD,IDPOS,IOFF,IEPSV)
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
              REAL(KIND=4) :: SF1((MGRID*MGRID*MGRID+1)/2)
              REAL(KIND=4) :: SF2((MGRID*MGRID*MGRID+1)/2)
              REAL(KIND=4) :: QMAP1(*)
              REAL(KIND=4) :: QMAP2(*)
              REAL(KIND=4) :: DEBMAP1(*)
              REAL(KIND=4) :: DEBMAP2(*)
              REAL(KIND=4) :: BNDX1(*)
              REAL(KIND=4) :: BNDX2(*)
              REAL(KIND=4) :: BNDX3(*)
              REAL(KIND=4) :: BNDX4(*)
              REAL(KIND=4) :: BNDX(*)
              REAL(KIND=4) :: BNDY(*)
              REAL(KIND=4) :: BNDZ(*)
              REAL(KIND=4) :: DB(6,*)
              INTEGER(KIND=4) :: IBGRD(4,*)
              INTEGER(KIND=4) :: IDPOS(*)
              INTEGER(KIND=4) :: IOFF(*)
              INTEGER(KIND=4) :: IEPSV(*)
            END SUBROUTINE DELPHI
          END INTERFACE 
        END MODULE DELPHI__genmod
