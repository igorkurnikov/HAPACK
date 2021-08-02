        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:46 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE AM_INDUCED_ADD_CART_TO_DIP__genmod
          INTERFACE 
            SUBROUTINE AM_INDUCED_ADD_CART_TO_DIP(ATM_CNT,IS_POLARIZABLE&
     &,CART_DIPOLE_FIELD,DIP_FIELD_D,DIP_FIELD_P)
              INTEGER(KIND=4), INTENT(IN) :: ATM_CNT
              LOGICAL(KIND=4), INTENT(IN) :: IS_POLARIZABLE(*)
              REAL(KIND=8), INTENT(IN) :: CART_DIPOLE_FIELD(3,*)
              REAL(KIND=8), INTENT(INOUT) :: DIP_FIELD_D(3,*)
              REAL(KIND=8), INTENT(INOUT) :: DIP_FIELD_P(3,*)
            END SUBROUTINE AM_INDUCED_ADD_CART_TO_DIP
          END INTERFACE 
        END MODULE AM_INDUCED_ADD_CART_TO_DIP__genmod
