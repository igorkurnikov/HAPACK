        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:50 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRTFARRD__genmod
          INTERFACE 
            FUNCTION WRTFARRD(NF,DARR,N,FORMSTR)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NF
              REAL(KIND=8) :: DARR(N)
              CHARACTER(*) :: FORMSTR
              INTEGER(KIND=4) :: WRTFARRD
            END FUNCTION WRTFARRD
          END INTERFACE 
        END MODULE WRTFARRD__genmod
