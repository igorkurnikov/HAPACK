        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:51 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LSQSTR__genmod
          INTERFACE 
            SUBROUTINE LSQSTR(NR,W,XP,X,E,IA,IROT,R)
              INTEGER(KIND=4) :: NR
              REAL(KIND=8) :: W(*)
              REAL(KIND=8) :: XP(*)
              REAL(KIND=8) :: X(*)
              REAL(KIND=8) :: E
              INTEGER(KIND=4) :: IA(*)
              INTEGER(KIND=4) :: IROT
              REAL(KIND=8) :: R(3,3)
            END SUBROUTINE LSQSTR
          END INTERFACE 
        END MODULE LSQSTR__genmod
