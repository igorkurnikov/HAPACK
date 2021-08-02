        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:56 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FORMP_HA__genmod
          INTERFACE 
            SUBROUTINE FORMP_HA(INIT,NDIM,NBASIS,NE,A,P)
              INTEGER(KIND=4) :: NDIM
              LOGICAL(KIND=4) :: INIT
              INTEGER(KIND=4) :: NBASIS
              INTEGER(KIND=4) :: NE
              REAL(KIND=8) :: A(NDIM,1)
              REAL(KIND=8) :: P(1)
            END SUBROUTINE FORMP_HA
          END INTERFACE 
        END MODULE FORMP_HA__genmod
