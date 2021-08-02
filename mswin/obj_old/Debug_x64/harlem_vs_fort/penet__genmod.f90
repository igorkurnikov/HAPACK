        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:55 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PENET__genmod
          INTERFACE 
            SUBROUTINE PENET(FACT,BINCOE,AMU,N,LA,BMU,M,LB,RR,PEN,PENPI,&
     &PEND,I)
              REAL(KIND=8) :: FACT(30)
              REAL(KIND=8) :: BINCOE(*)
              REAL(KIND=8) :: AMU
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: LA
              REAL(KIND=8) :: BMU
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: LB
              REAL(KIND=8) :: RR
              REAL(KIND=8) :: PEN
              REAL(KIND=8) :: PENPI
              REAL(KIND=8) :: PEND
              INTEGER(KIND=4) :: I
            END SUBROUTINE PENET
          END INTERFACE 
        END MODULE PENET__genmod
