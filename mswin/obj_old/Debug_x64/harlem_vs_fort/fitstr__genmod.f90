        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:51 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FITSTR__genmod
          INTERFACE 
            SUBROUTINE FITSTR(NAT,X1,X2,EPS,IROT,R,CC,W,IA)
              INTEGER(KIND=4) :: NAT
              REAL(KIND=8) :: X1(*)
              REAL(KIND=8) :: X2(*)
              REAL(KIND=8) :: EPS
              INTEGER(KIND=4) :: IROT
              REAL(KIND=8) :: R(3,3)
              REAL(KIND=8) :: CC(*)
              REAL(KIND=8) :: W(*)
              INTEGER(KIND=4) :: IA(*)
            END SUBROUTINE FITSTR
          END INTERFACE 
        END MODULE FITSTR__genmod
