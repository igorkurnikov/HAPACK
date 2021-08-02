        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:50 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RDFARRD__genmod
          INTERFACE 
            FUNCTION RDFARRD(NF,DARR,N,FORM)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NF
              REAL(KIND=8) :: DARR(N)
              CHARACTER(*) :: FORM
              INTEGER(KIND=4) :: RDFARRD
            END FUNCTION RDFARRD
          END INTERFACE 
        END MODULE RDFARRD__genmod
