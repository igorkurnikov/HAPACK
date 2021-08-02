        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:51 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE VA13AD__genmod
          INTERFACE 
            SUBROUTINE VA13AD(IPTR,FUNC,N,X,F,G,SCALE,ACC,W)
              INTEGER(KIND=4) :: IPTR
              EXTERNAL FUNC
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: X(*)
              REAL(KIND=8) :: F
              REAL(KIND=8) :: G(*)
              REAL(KIND=8) :: SCALE(*)
              REAL(KIND=8) :: ACC
              REAL(KIND=8) :: W(*)
            END SUBROUTINE VA13AD
          END INTERFACE 
        END MODULE VA13AD__genmod
