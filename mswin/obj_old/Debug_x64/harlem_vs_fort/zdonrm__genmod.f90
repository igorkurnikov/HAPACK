        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:54 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZDONRM__genmod
          INTERFACE 
            FUNCTION ZDONRM(METHOD,NATOMS,CZ,IAN,C,GSS,GSD,GDD)
              INTEGER(KIND=4) :: NATOMS
              INTEGER(KIND=4) :: METHOD
              REAL(KIND=8) :: CZ(*)
              INTEGER(KIND=4) :: IAN(*)
              REAL(KIND=8) :: C(3,*)
              REAL(KIND=8) :: GSS(*)
              REAL(KIND=8) :: GSD(NATOMS,NATOMS)
              REAL(KIND=8) :: GDD(*)
              REAL(KIND=8) :: ZDONRM
            END FUNCTION ZDONRM
          END INTERFACE 
        END MODULE ZDONRM__genmod
