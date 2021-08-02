        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:55 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SSZ__genmod
          INTERFACE 
            FUNCTION SSZ(NN1,NN2,LL1,LL2,M,AMU_X,BMU_X,FACT,BINCOE,LG)
              INTEGER(KIND=4) :: NN1
              INTEGER(KIND=4) :: NN2
              INTEGER(KIND=4) :: LL1
              INTEGER(KIND=4) :: LL2
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: AMU_X
              REAL(KIND=8) :: BMU_X
              REAL(KIND=8) :: FACT(30)
              REAL(KIND=8) :: BINCOE(*)
              INTEGER(KIND=4) :: LG
              REAL(KIND=8) :: SSZ
            END FUNCTION SSZ
          END INTERFACE 
        END MODULE SSZ__genmod
