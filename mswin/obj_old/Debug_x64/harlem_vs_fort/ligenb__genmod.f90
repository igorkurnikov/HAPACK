        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:50 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LIGENB__genmod
          INTERFACE 
            SUBROUTINE LIGENB(A,VEC,EIG,IA,N,NDIM)
              INTEGER(KIND=4) :: NDIM
              REAL(KIND=8) :: A(*)
              REAL(KIND=8) :: VEC(NDIM,*)
              REAL(KIND=8) :: EIG(*)
              INTEGER(KIND=4) :: IA(*)
              INTEGER(KIND=4) :: N
            END SUBROUTINE LIGENB
          END INTERFACE 
        END MODULE LIGENB__genmod
