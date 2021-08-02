        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:21 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE AM_VAL_BCUCOF__genmod
          INTERFACE 
            SUBROUTINE AM_VAL_BCUCOF(Y,Y1,Y2,Y12,D1,D2,C)
              REAL(KIND=8), INTENT(IN) :: Y(4)
              REAL(KIND=8), INTENT(IN) :: Y1(4)
              REAL(KIND=8), INTENT(IN) :: Y2(4)
              REAL(KIND=8), INTENT(IN) :: Y12(4)
              REAL(KIND=8), INTENT(IN) :: D1
              REAL(KIND=8), INTENT(IN) :: D2
              REAL(KIND=8), INTENT(OUT) :: C(4,4)
            END SUBROUTINE AM_VAL_BCUCOF
          END INTERFACE 
        END MODULE AM_VAL_BCUCOF__genmod
