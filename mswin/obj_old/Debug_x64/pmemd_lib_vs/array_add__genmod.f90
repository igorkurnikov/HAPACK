        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:49 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ARRAY_ADD__genmod
          INTERFACE 
            SUBROUTINE ARRAY_ADD(A,B,NUM)
              REAL(KIND=8), INTENT(INOUT) :: A(*)
              REAL(KIND=8), INTENT(IN) :: B(*)
              INTEGER(KIND=4), INTENT(IN) :: NUM
            END SUBROUTINE ARRAY_ADD
          END INTERFACE 
        END MODULE ARRAY_ADD__genmod
