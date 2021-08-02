        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:46 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GB_FORCE_PROXY__genmod
          INTERFACE 
            SUBROUTINE GB_FORCE_PROXY(ATM_CNT,CRD,FRC,MASS,SI,NCALLS,   &
     &ATM_JRC,ATM_XC,ATM_WEIGHT,IGROUP,BELLY_ATM_CNT,NATC,ATM_OWNER_MAP,&
     &MY_ATM_LST,MY_ATM_CNT)
              INTEGER(KIND=4) :: ATM_CNT
              REAL(KIND=8) :: CRD(3,ATM_CNT)
              REAL(KIND=8) :: FRC(3,ATM_CNT)
              REAL(KIND=8) :: MASS(ATM_CNT)
              REAL(KIND=8) :: SI(39)
              INTEGER(KIND=4) :: NCALLS
              INTEGER(KIND=4) :: ATM_JRC(*)
              REAL(KIND=8) :: ATM_XC(3,*)
              REAL(KIND=8) :: ATM_WEIGHT(*)
              INTEGER(KIND=4) :: IGROUP(*)
              INTEGER(KIND=4) :: BELLY_ATM_CNT
              INTEGER(KIND=4) :: NATC
              INTEGER(KIND=4) :: ATM_OWNER_MAP(ATM_CNT)
              INTEGER(KIND=4) :: MY_ATM_LST(*)
              INTEGER(KIND=4) :: MY_ATM_CNT
            END SUBROUTINE GB_FORCE_PROXY
          END INTERFACE 
        END MODULE GB_FORCE_PROXY__genmod
