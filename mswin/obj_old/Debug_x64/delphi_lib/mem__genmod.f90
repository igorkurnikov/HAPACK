        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:20 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MEM__genmod
          INTERFACE 
            SUBROUTINE MEM(XX,MGRID,PHIMAP,PHIMAP1,PHIMAP2,PHIMAP3,     &
     &ATMCRG,IEPSMP,IEPSMP2,IDEBMAP,ISLICE)
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
              INTEGER(KIND=4) :: ISLICE(2,*)
            END SUBROUTINE MEM
          END INTERFACE 
        END MODULE MEM__genmod
