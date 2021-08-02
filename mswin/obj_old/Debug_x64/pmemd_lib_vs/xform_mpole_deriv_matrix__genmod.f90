        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:39 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE XFORM_MPOLE_DERIV_MATRIX__genmod
          INTERFACE 
            SUBROUTINE XFORM_MPOLE_DERIV_MATRIX(A_XY,DA_XY,DMPOLE_XY,   &
     &ORDER)
              INTEGER(KIND=4) :: ORDER
              REAL(KIND=8) :: A_XY(3,3)
              REAL(KIND=8) :: DA_XY(3,3)
              REAL(KIND=8) :: DMPOLE_XY(ORDER,ORDER)
            END SUBROUTINE XFORM_MPOLE_DERIV_MATRIX
          END INTERFACE 
        END MODULE XFORM_MPOLE_DERIV_MATRIX__genmod
