        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:49 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE OPENDX__genmod
          INTERFACE 
            SUBROUTINE OPENDX(LUDX,NAME,NELEM,NREC,STATUS,LRDX,NBDX,    &
     &OLDDX)
              INTEGER(KIND=4) :: LUDX
              CHARACTER(*) :: NAME
              INTEGER(KIND=4) :: NELEM
              INTEGER(KIND=4) :: NREC
              CHARACTER(*) :: STATUS
              INTEGER(KIND=4) :: LRDX
              INTEGER(KIND=4) :: NBDX
              LOGICAL(KIND=4) :: OLDDX
            END SUBROUTINE OPENDX
          END INTERFACE 
        END MODULE OPENDX__genmod
