        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:50 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE VA13CD__genmod
          INTERFACE 
            SUBROUTINE VA13CD(IPTR,FUNC,N,X,F,G,SCALE,ACC,H,D,W,XA,GA,XB&
     &,GB)
              INTEGER(KIND=4) :: IPTR
              EXTERNAL FUNC
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: X(*)
              REAL(KIND=8) :: F
              REAL(KIND=8) :: G(*)
              REAL(KIND=8) :: SCALE(*)
              REAL(KIND=8) :: ACC
              REAL(KIND=8) :: H(*)
              REAL(KIND=8) :: D(*)
              REAL(KIND=8) :: W(*)
              REAL(KIND=8) :: XA(*)
              REAL(KIND=8) :: GA(*)
              REAL(KIND=8) :: XB(*)
              REAL(KIND=8) :: GB(*)
            END SUBROUTINE VA13CD
          END INTERFACE 
        END MODULE VA13CD__genmod
