        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:57 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FORMG_GS__genmod
          INTERFACE 
            SUBROUTINE FORMG_GS(NT,IB,B,G)
              INTEGER(KIND=4) :: NT
              INTEGER(KIND=4) :: IB(4,NT)
              REAL(KIND=8) :: B(3,4,NT)
              REAL(KIND=8) :: G(NT,NT)
            END SUBROUTINE FORMG_GS
          END INTERFACE 
        END MODULE FORMG_GS__genmod
