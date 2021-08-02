        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:23 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MINFOR__genmod
          INTERFACE 
            SUBROUTINE MINFOR(N,X,F,G,SCALE,ACC,H,D,W,XA,GA,XB,GB,MAXFUN&
     &,NFUN)
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
              INTEGER(KIND=4) :: MAXFUN
              INTEGER(KIND=4) :: NFUN
            END SUBROUTINE MINFOR
          END INTERFACE 
        END MODULE MINFOR__genmod
