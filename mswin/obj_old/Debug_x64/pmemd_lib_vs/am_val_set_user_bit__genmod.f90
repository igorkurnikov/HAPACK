        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:49 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE AM_VAL_SET_USER_BIT__genmod
          INTERFACE 
            SUBROUTINE AM_VAL_SET_USER_BIT(DO_BOND,DO_UREYB,DO_REG_ANGLE&
     &,DO_TRIG_ANGLE,DO_OPBEND_ANGLE,DO_TORSIONS,DO_STR_TORSIONS,       &
     &DO_PITORSIONS,DO_STRETCH_BEND,DO_TORSION_TORSION)
              INTEGER(KIND=4), INTENT(IN) :: DO_BOND
              INTEGER(KIND=4), INTENT(IN) :: DO_UREYB
              INTEGER(KIND=4), INTENT(IN) :: DO_REG_ANGLE
              INTEGER(KIND=4), INTENT(IN) :: DO_TRIG_ANGLE
              INTEGER(KIND=4), INTENT(IN) :: DO_OPBEND_ANGLE
              INTEGER(KIND=4), INTENT(IN) :: DO_TORSIONS
              INTEGER(KIND=4), INTENT(IN) :: DO_STR_TORSIONS
              INTEGER(KIND=4), INTENT(IN) :: DO_PITORSIONS
              INTEGER(KIND=4), INTENT(IN) :: DO_STRETCH_BEND
              INTEGER(KIND=4), INTENT(IN) :: DO_TORSION_TORSION
            END SUBROUTINE AM_VAL_SET_USER_BIT
          END INTERFACE 
        END MODULE AM_VAL_SET_USER_BIT__genmod
