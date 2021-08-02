        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:52 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CFUNC__genmod
          INTERFACE 
            SUBROUTINE CFUNC(IA,IB,IC,ID,IE,AMU,BMU,SNAG,BINCOE,A,B)
              INTEGER(KIND=4) :: IA
              INTEGER(KIND=4) :: IB
              INTEGER(KIND=4) :: IC
              INTEGER(KIND=4) :: ID
              INTEGER(KIND=4) :: IE
              REAL(KIND=8) :: AMU
              REAL(KIND=8) :: BMU
              REAL(KIND=8) :: SNAG
              REAL(KIND=8) :: BINCOE(*)
              REAL(KIND=8) :: A(*)
              REAL(KIND=8) :: B(*)
            END SUBROUTINE CFUNC
          END INTERFACE 
        END MODULE CFUNC__genmod
