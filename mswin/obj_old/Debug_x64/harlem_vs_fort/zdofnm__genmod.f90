        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:52 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZDOFNM__genmod
          INTERFACE 
            SUBROUTINE ZDOFNM(NATOMS,CZ,C,FXYZ)
              INTEGER(KIND=4) :: NATOMS
              REAL(KIND=8) :: CZ(*)
              REAL(KIND=8) :: C(3,*)
              REAL(KIND=8) :: FXYZ(3,*)
            END SUBROUTINE ZDOFNM
          END INTERFACE 
        END MODULE ZDOFNM__genmod
