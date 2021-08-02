        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:46 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE AM_INDUCED_FIELDS_TO_IND_DIPS__genmod
          INTERFACE 
            SUBROUTINE AM_INDUCED_FIELDS_TO_IND_DIPS(ATM_CNT,           &
     &IS_POLARIZABLE,DIP_FIELD_D,DIP_FIELD_P,POLARIZABILITY,HPOLAR,     &
     &POLARIZABILITY_CORR,IND_DIP_D,IND_DIP_P)
              INTEGER(KIND=4), INTENT(IN) :: ATM_CNT
              LOGICAL(KIND=4), INTENT(IN) :: IS_POLARIZABLE(*)
              REAL(KIND=8), INTENT(IN) :: DIP_FIELD_D(3,*)
              REAL(KIND=8), INTENT(IN) :: DIP_FIELD_P(3,*)
              REAL(KIND=8), INTENT(IN) :: POLARIZABILITY(*)
              REAL(KIND=8), INTENT(IN) :: HPOLAR(*)
              REAL(KIND=8), INTENT(IN) :: POLARIZABILITY_CORR(*)
              REAL(KIND=8), INTENT(OUT) :: IND_DIP_D(3,*)
              REAL(KIND=8), INTENT(OUT) :: IND_DIP_P(3,*)
            END SUBROUTINE AM_INDUCED_FIELDS_TO_IND_DIPS
          END INTERFACE 
        END MODULE AM_INDUCED_FIELDS_TO_IND_DIPS__genmod
