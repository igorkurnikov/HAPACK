        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug 08 10:08:56 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EXTRDM__genmod
          INTERFACE 
            SUBROUTINE EXTRDM(IOUT,IPRINT,ICOUNT,N,NMAT,NTT,IFLAG,IEXTP,&
     &RMSDP,PCUR,SCR,DELTAP)
              INTEGER(KIND=4) :: NTT
              INTEGER(KIND=4) :: NMAT
              INTEGER(KIND=4) :: IOUT
              INTEGER(KIND=4) :: IPRINT
              INTEGER(KIND=4) :: ICOUNT
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IFLAG
              INTEGER(KIND=4) :: IEXTP
              REAL(KIND=8) :: RMSDP
              REAL(KIND=8) :: PCUR(NTT,NMAT)
              REAL(KIND=8) :: SCR(NTT,NMAT)
              REAL(KIND=8) :: DELTAP(NTT,NMAT,3)
            END SUBROUTINE EXTRDM
          END INTERFACE 
        END MODULE EXTRDM__genmod
