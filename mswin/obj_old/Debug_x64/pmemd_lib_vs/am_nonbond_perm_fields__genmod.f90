        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:49 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE AM_NONBOND_PERM_FIELDS__genmod
          INTERFACE 
            SUBROUTINE AM_NONBOND_PERM_FIELDS(ATM_CNT,IS_POLARIZABLE,CRD&
     &,DIP_FIELD_D,DIP_FIELD_P,ADJ_DIP_DIP_TENSOR,TRANVEC)
              INTEGER(KIND=4), INTENT(IN) :: ATM_CNT
              LOGICAL(KIND=4), INTENT(IN) :: IS_POLARIZABLE(*)
              REAL(KIND=8), INTENT(IN) :: CRD(3,*)
              REAL(KIND=8), INTENT(OUT) :: DIP_FIELD_D(3,ATM_CNT)
              REAL(KIND=8), INTENT(OUT) :: DIP_FIELD_P(3,ATM_CNT)
              REAL(KIND=8), INTENT(OUT) :: ADJ_DIP_DIP_TENSOR(6,*)
              REAL(KIND=8), INTENT(IN) :: TRANVEC(*)
            END SUBROUTINE AM_NONBOND_PERM_FIELDS
          END INTERFACE 
        END MODULE AM_NONBOND_PERM_FIELDS__genmod
