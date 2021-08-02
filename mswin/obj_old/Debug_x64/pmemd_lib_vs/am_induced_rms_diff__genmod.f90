        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:46 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE AM_INDUCED_RMS_DIFF__genmod
          INTERFACE 
            SUBROUTINE AM_INDUCED_RMS_DIFF(ATM_CNT,IS_POLARIZABLE,RMS1, &
     &RMS2,IND_DIP_D,IND_DIP_P,OLD_DIP_D,OLD_DIP_P)
              INTEGER(KIND=4), INTENT(IN) :: ATM_CNT
              LOGICAL(KIND=4), INTENT(IN) :: IS_POLARIZABLE(*)
              REAL(KIND=8), INTENT(OUT) :: RMS1
              REAL(KIND=8), INTENT(OUT) :: RMS2
              REAL(KIND=8), INTENT(IN) :: IND_DIP_D(3,*)
              REAL(KIND=8), INTENT(IN) :: IND_DIP_P(3,*)
              REAL(KIND=8), INTENT(IN) :: OLD_DIP_D(3,*)
              REAL(KIND=8), INTENT(IN) :: OLD_DIP_P(3,*)
            END SUBROUTINE AM_INDUCED_RMS_DIFF
          END INTERFACE 
        END MODULE AM_INDUCED_RMS_DIFF__genmod
