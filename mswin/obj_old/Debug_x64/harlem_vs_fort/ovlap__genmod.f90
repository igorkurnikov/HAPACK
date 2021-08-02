        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:55 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE OVLAP__genmod
          INTERFACE 
            SUBROUTINE OVLAP(N1,L1,AMU,N2,L2,BMU,R,S,SP,SD,SF,FACT,     &
     &BINCOE,PP)
              INTEGER(KIND=4) :: N1
              INTEGER(KIND=4) :: L1
              REAL(KIND=8) :: AMU
              INTEGER(KIND=4) :: N2
              INTEGER(KIND=4) :: L2
              REAL(KIND=8) :: BMU
              REAL(KIND=8) :: R
              REAL(KIND=8) :: S
              REAL(KIND=8) :: SP
              REAL(KIND=8) :: SD
              REAL(KIND=8) :: SF
              REAL(KIND=8) :: FACT(30)
              REAL(KIND=8) :: BINCOE(*)
              REAL(KIND=8) :: PP
            END SUBROUTINE OVLAP
          END INTERFACE 
        END MODULE OVLAP__genmod
