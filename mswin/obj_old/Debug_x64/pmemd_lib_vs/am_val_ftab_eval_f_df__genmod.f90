        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 14 06:12:20 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE AM_VAL_FTAB_EVAL_F_DF__genmod
          INTERFACE 
            SUBROUTINE AM_VAL_FTAB_EVAL_F_DF(NLIST,DEGREE,COEFF,ARG,FUNC&
     &,DFUNC_DARG)
              INTEGER(KIND=4), INTENT(IN) :: DEGREE
              INTEGER(KIND=4), INTENT(IN) :: NLIST
              REAL(KIND=8), INTENT(IN) :: COEFF(0:DEGREE)
              REAL(KIND=8), INTENT(IN) :: ARG(NLIST)
              REAL(KIND=8), INTENT(OUT) :: FUNC(NLIST)
              REAL(KIND=8), INTENT(OUT) :: DFUNC_DARG(NLIST)
            END SUBROUTINE AM_VAL_FTAB_EVAL_F_DF
          END INTERFACE 
        END MODULE AM_VAL_FTAB_EVAL_F_DF__genmod
