        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:50 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PME_FORCE_PROXY__genmod
          INTERFACE 
            SUBROUTINE PME_FORCE_PROXY(ATM_CNT,CRD,SAVED_CRD,BOX,       &
     &SAVED_BOX,VEL,FRC,MASS,NEW_LIST,ATM_JRC,ATM_XC,ATM_WEIGHT,IGROUP, &
     &NATC,SI,VIRIAL,EKCMT,PME_ERR_EST,ATM_OWNER_MAP,MY_ATM_LST,        &
     &MY_ATM_CNT,TRANVEC,IMIN_PAR)
              INTEGER(KIND=4) :: ATM_CNT
              REAL(KIND=8) :: CRD(3,ATM_CNT)
              REAL(KIND=8) :: SAVED_CRD(3,ATM_CNT)
              REAL(KIND=8) :: BOX(3)
              REAL(KIND=8) :: SAVED_BOX(3)
              REAL(KIND=8) :: VEL(3,ATM_CNT)
              REAL(KIND=8) :: FRC(3,ATM_CNT)
              REAL(KIND=8) :: MASS(ATM_CNT)
              INTEGER(KIND=4) :: NEW_LIST
              INTEGER(KIND=4) :: ATM_JRC(*)
              REAL(KIND=8) :: ATM_XC(3,*)
              REAL(KIND=8) :: ATM_WEIGHT(*)
              INTEGER(KIND=4) :: IGROUP(*)
              INTEGER(KIND=4) :: NATC
              REAL(KIND=8) :: SI(39)
              REAL(KIND=8) :: VIRIAL(3)
              REAL(KIND=8) :: EKCMT(3)
              REAL(KIND=8) :: PME_ERR_EST
              INTEGER(KIND=4) :: ATM_OWNER_MAP(ATM_CNT)
              INTEGER(KIND=4) :: MY_ATM_LST(*)
              INTEGER(KIND=4) :: MY_ATM_CNT
              REAL(KIND=8) :: TRANVEC(*)
              INTEGER(KIND=4) :: IMIN_PAR
            END SUBROUTINE PME_FORCE_PROXY
          END INTERFACE 
        END MODULE PME_FORCE_PROXY__genmod
