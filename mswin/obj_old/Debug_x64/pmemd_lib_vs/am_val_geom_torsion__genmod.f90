        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:21 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE AM_VAL_GEOM_TORSION__genmod
          INTERFACE 
            SUBROUTINE AM_VAL_GEOM_TORSION(CRD_ABCD,GRADPHI_ABCD,COSPHI,&
     &SINPHI)
              REAL(KIND=8), INTENT(IN) :: CRD_ABCD(12)
              REAL(KIND=8), INTENT(OUT) :: GRADPHI_ABCD(12)
              REAL(KIND=8), INTENT(OUT) :: COSPHI
              REAL(KIND=8), INTENT(OUT) :: SINPHI
            END SUBROUTINE AM_VAL_GEOM_TORSION
          END INTERFACE 
        END MODULE AM_VAL_GEOM_TORSION__genmod
