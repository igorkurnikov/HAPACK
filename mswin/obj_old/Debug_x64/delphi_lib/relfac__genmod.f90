        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:20 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RELFAC__genmod
          INTERFACE 
            SUBROUTINE RELFAC(XX,MGRID,PHIMAP,PHIMAP1,PHIMAP2,PHIMAP3,  &
     &ATMCRG,IEPSMP,IEPSMP2,IDEBMAP,BNDX,BNDY,BNDZ,NLIT,NNIT,IPER,IDPOS,&
     &DB,SF1,SF2,ICOUNT2A,ICOUNT2B,RIONST,SPEC)
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
              INTEGER(KIND=4) :: ICOUNT2A
              INTEGER(KIND=4) :: ICOUNT2B
              REAL(KIND=4) :: RIONST
              REAL(KIND=4) :: SPEC
            END SUBROUTINE RELFAC
          END INTERFACE 
        END MODULE RELFAC__genmod
