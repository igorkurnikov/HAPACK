        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:49 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE AM_NONBOND_DIP_DIP_FIELDS__genmod
          INTERFACE 
            SUBROUTINE AM_NONBOND_DIP_DIP_FIELDS(ATM_CNT,CRD,IND_DIP_D, &
     &IND_DIP_P,ADJ_DIP_DIP_TENSOR,DIP_FIELD_D,DIP_FIELD_P)
              INTEGER(KIND=4), INTENT(IN) :: ATM_CNT
              REAL(KIND=8), INTENT(IN) :: CRD(3,*)
              REAL(KIND=8), INTENT(IN) :: IND_DIP_D(3,*)
              REAL(KIND=8), INTENT(IN) :: IND_DIP_P(3,*)
              REAL(KIND=8), INTENT(IN) :: ADJ_DIP_DIP_TENSOR(6,*)
              REAL(KIND=8), INTENT(OUT) :: DIP_FIELD_D(3,ATM_CNT)
              REAL(KIND=8), INTENT(OUT) :: DIP_FIELD_P(3,ATM_CNT)
            END SUBROUTINE AM_NONBOND_DIP_DIP_FIELDS
          END INTERFACE 
        END MODULE AM_NONBOND_DIP_DIP_FIELDS__genmod
