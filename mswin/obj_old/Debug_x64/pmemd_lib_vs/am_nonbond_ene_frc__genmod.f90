        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:49 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE AM_NONBOND_ENE_FRC__genmod
          INTERFACE 
            SUBROUTINE AM_NONBOND_ENE_FRC(ATM_CNT,CRD,ENE_PERM,ENE_IND, &
     &ENE_VDW,ENE_VDW_14,FRC,IMG_FRC,VIR_TENSOR,NETFRCS,ATM_OWNER_MAP,  &
     &TRANVEC)
              INTEGER(KIND=4), INTENT(IN) :: ATM_CNT
              REAL(KIND=8), INTENT(IN) :: CRD(3,*)
              REAL(KIND=8), INTENT(OUT) :: ENE_PERM
              REAL(KIND=8), INTENT(OUT) :: ENE_IND
              REAL(KIND=8), INTENT(OUT) :: ENE_VDW
              REAL(KIND=8), INTENT(OUT) :: ENE_VDW_14
              REAL(KIND=8), INTENT(INOUT) :: FRC(3,*)
              REAL(KIND=8), INTENT(INOUT) :: IMG_FRC(3,*)
              REAL(KIND=8), INTENT(INOUT) :: VIR_TENSOR(3,3)
              REAL(KIND=8), INTENT(OUT) :: NETFRCS(3)
              INTEGER(KIND=4), INTENT(INOUT) :: ATM_OWNER_MAP(ATM_CNT)
              REAL(KIND=8), INTENT(IN) :: TRANVEC(*)
            END SUBROUTINE AM_NONBOND_ENE_FRC
          END INTERFACE 
        END MODULE AM_NONBOND_ENE_FRC__genmod
