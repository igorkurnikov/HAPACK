        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:57 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TRANF_GS__genmod
          INTERFACE 
            SUBROUTINE TRANF_GS(NPARM,FX,F,IB,B,G)
              INTEGER(KIND=4) :: NPARM
              REAL(KIND=8) :: FX(3,*)
              REAL(KIND=8) :: F(*)
              INTEGER(KIND=4) :: IB(4,NPARM,*)
              REAL(KIND=8) :: B(3,4,NPARM,*)
              REAL(KIND=8) :: G(NPARM,NPARM)
            END SUBROUTINE TRANF_GS
          END INTERFACE 
        END MODULE TRANF_GS__genmod
