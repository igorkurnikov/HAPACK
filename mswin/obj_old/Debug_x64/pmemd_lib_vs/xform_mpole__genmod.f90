        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:40 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE XFORM_MPOLE__genmod
          INTERFACE 
            SUBROUTINE XFORM_MPOLE(MPOLE_XY,DIMXY,MPOLE_IN,MPOLE_OUT,   &
     &ORDER)
              INTEGER(KIND=4) :: DIMXY
              REAL(KIND=8) :: MPOLE_XY(DIMXY,DIMXY)
              REAL(KIND=8) :: MPOLE_IN(*)
              REAL(KIND=8) :: MPOLE_OUT(*)
              INTEGER(KIND=4) :: ORDER
            END SUBROUTINE XFORM_MPOLE
          END INTERFACE 
        END MODULE XFORM_MPOLE__genmod
