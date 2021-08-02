        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:49 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE OPENPM__genmod
          INTERFACE 
            SUBROUTINE OPENPM(NX,NY,NZ,XSTART,XEND,YSTART,YEND,ZSTART,  &
     &ZEND,FNAME)
              INTEGER(KIND=4) :: NX
              INTEGER(KIND=4) :: NY
              INTEGER(KIND=4) :: NZ
              REAL(KIND=4) :: XSTART
              REAL(KIND=4) :: XEND
              REAL(KIND=4) :: YSTART
              REAL(KIND=4) :: YEND
              REAL(KIND=4) :: ZSTART
              REAL(KIND=4) :: ZEND
              CHARACTER(*) :: FNAME
            END SUBROUTINE OPENPM
          END INTERFACE 
        END MODULE OPENPM__genmod
