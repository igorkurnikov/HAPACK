        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:55 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TORS_GS__genmod
          INTERFACE 
            SUBROUTINE TORS_GS(NOINT,I,J,K,L,B,IB,C)
              INTEGER(KIND=4) :: NOINT
              INTEGER(KIND=4) :: I
              INTEGER(KIND=4) :: J
              INTEGER(KIND=4) :: K
              INTEGER(KIND=4) :: L
              REAL(KIND=8) :: B(3,4,*)
              INTEGER(KIND=4) :: IB(4,*)
              REAL(KIND=8) :: C(1)
            END SUBROUTINE TORS_GS
          END INTERFACE 
        END MODULE TORS_GS__genmod
