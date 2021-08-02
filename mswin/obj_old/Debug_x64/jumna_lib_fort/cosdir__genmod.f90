        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:23 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COSDIR__genmod
          INTERFACE 
            SUBROUTINE COSDIR(COR,NN,I,J,K,CDIR)
              INTEGER(KIND=4) :: NN
              REAL(KIND=8) :: COR(NN,3)
              INTEGER(KIND=4) :: I
              INTEGER(KIND=4) :: J
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: CDIR(9)
            END SUBROUTINE COSDIR
          END INTERFACE 
        END MODULE COSDIR__genmod
