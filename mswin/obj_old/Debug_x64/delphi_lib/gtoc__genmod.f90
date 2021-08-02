        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:20 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GTOC__genmod
          INTERFACE 
            SUBROUTINE GTOC(MGRID,OLDMID,SCALE,G,C)
              INTEGER(KIND=4) :: MGRID
              REAL(KIND=4) :: OLDMID(3)
              REAL(KIND=4) :: SCALE
              REAL(KIND=4) :: G(3)
              REAL(KIND=4) :: C(3)
            END SUBROUTINE GTOC
          END INTERFACE 
        END MODULE GTOC__genmod
