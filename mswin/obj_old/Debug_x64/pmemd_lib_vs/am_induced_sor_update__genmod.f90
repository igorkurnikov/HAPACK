        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:46 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE AM_INDUCED_SOR_UPDATE__genmod
          INTERFACE 
            SUBROUTINE AM_INDUCED_SOR_UPDATE(ATM_CNT,IS_POLARIZABLE,    &
     &SOR_COEFF,IND_DIP_D,IND_DIP_P,OLD_DIP_D,OLD_DIP_P)
              INTEGER(KIND=4), INTENT(IN) :: ATM_CNT
              LOGICAL(KIND=4), INTENT(IN) :: IS_POLARIZABLE(*)
              REAL(KIND=8), INTENT(IN) :: SOR_COEFF
              REAL(KIND=8), INTENT(INOUT) :: IND_DIP_D(3,*)
              REAL(KIND=8), INTENT(INOUT) :: IND_DIP_P(3,*)
              REAL(KIND=8), INTENT(IN) :: OLD_DIP_D(3,*)
              REAL(KIND=8), INTENT(IN) :: OLD_DIP_P(3,*)
            END SUBROUTINE AM_INDUCED_SOR_UPDATE
          END INTERFACE 
        END MODULE AM_INDUCED_SOR_UPDATE__genmod
