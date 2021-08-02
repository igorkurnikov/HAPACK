        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:20 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NITIT__genmod
          INTERFACE 
            SUBROUTINE NITIT(XX,MGRID,PHIMAP,PHIMAP1,PHIMAP2,PHIMAP3,   &
     &ATMCRG,IEPSMP,IEPSMP2,IDEBMAP,QMAP1,QMAP2,DEBMAP1,DEBMAP2,BNDX1,  &
     &BNDX2,BNDX3,BNDX4,BNDX,BNDY,BNDZ,NLIT,NNIT,IPER,IDPOS,DB,SF1,SF2, &
     &IQPOS,QVAL,ICOUNT2A,ICOUNT2B,ICOUNT1A,ICOUNT1B,RIONST,SPEC,ICON1, &
     &ICON2,DEBFCT,EPSOUT)
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
              INTEGER(KIND=4) :: NLIT
              INTEGER(KIND=4) :: NNIT
              LOGICAL(KIND=4) :: IPER(3)
              INTEGER(KIND=4) :: IDPOS(*)
              REAL(KIND=4) :: DB(6,*)
              REAL(KIND=4) :: SF1(*)
              REAL(KIND=4) :: SF2(*)
              INTEGER(KIND=4) :: IQPOS(200000)
              REAL(KIND=4) :: QVAL(200000)
              INTEGER(KIND=4) :: ICOUNT2A
              INTEGER(KIND=4) :: ICOUNT2B
              INTEGER(KIND=4) :: ICOUNT1A
              INTEGER(KIND=4) :: ICOUNT1B
              REAL(KIND=4) :: RIONST
              REAL(KIND=4) :: SPEC
              INTEGER(KIND=4) :: ICON1
              INTEGER(KIND=4) :: ICON2
              REAL(KIND=4) :: DEBFCT
              REAL(KIND=4) :: EPSOUT
            END SUBROUTINE NITIT
          END INTERFACE 
        END MODULE NITIT__genmod
