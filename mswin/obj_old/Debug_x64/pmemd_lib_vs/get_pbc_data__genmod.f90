        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:45 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_PBC_DATA__genmod
          INTERFACE 
            SUBROUTINE GET_PBC_DATA(IS_ORTHOG_N,PBC_ALPHA_N,PBC_BETA_N, &
     &PBC_GAMMA_N,RECIP_N,UCELL_N,CUT_FACTOR_N,RECLNG_N,UC_VOLUME_N,    &
     &UC_SPHERE_N)
              INTEGER(KIND=4) :: IS_ORTHOG_N
              REAL(KIND=8) :: PBC_ALPHA_N
              REAL(KIND=8) :: PBC_BETA_N
              REAL(KIND=8) :: PBC_GAMMA_N
              REAL(KIND=8) :: RECIP_N(3,3)
              REAL(KIND=8) :: UCELL_N(3,3)
              REAL(KIND=8) :: CUT_FACTOR_N(3)
              REAL(KIND=8) :: RECLNG_N(3)
              REAL(KIND=8) :: UC_VOLUME_N
              REAL(KIND=8) :: UC_SPHERE_N
            END SUBROUTINE GET_PBC_DATA
          END INTERFACE 
        END MODULE GET_PBC_DATA__genmod
