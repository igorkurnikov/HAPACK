        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:49 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ARRAY_COPY__genmod
          INTERFACE 
            SUBROUTINE ARRAY_COPY(A,B,NUM)
              INTEGER(KIND=4), INTENT(IN) :: NUM
              REAL(KIND=8), INTENT(INOUT) :: A(1:NUM)
              REAL(KIND=8), INTENT(IN) :: B(1:NUM)
            END SUBROUTINE ARRAY_COPY
          END INTERFACE 
        END MODULE ARRAY_COPY__genmod
