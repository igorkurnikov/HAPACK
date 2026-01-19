C       =============================================================
C       GEPOL93 (GEometria POLihedro,1993)
C       Version 8
C       Burjassot September 17, 1993.
C       =============================================================
C*** Written by
C
C       J.L. Pascual-Ahuir, E. Silla and I. Tunon
C       Departamento de Quimica Fisica.
C       Facultad de Quimica.
C       Universidad de Valencia.
C       C/Dr.Moliner 50.
C       Burjassot (Valencia) 46100
C       SPAIN
C       
C       Phone number: (34)-6-3864332
C       FAX   number: (34)-6-3864564
C       EMAIL   PASCUAL@EVALUN11
C                 SILLA@EVALUN11
C
	subroutine gsurf(KSURF,RMIN,OFAC,RD,NDIV,ASS1,NATOM,GHOST,NP,VOL)
      IMPLICIT NONE
      
      LOGICAL ASS1
      LOGICAL GHOST

      INTEGER*2 IUSE

      INTEGER*4 ISA,ISO,ITO
      INTEGER*4 MC,MV
      INTEGER*4 NATOM,NCOR,NDIV,NP
       
      REAL*4 AP
      REAL*4 DVEC
      REAL*4 OFAC
      REAL*4 RD,RE,RMIN
      REAL*4 XC1,XE,XP
      REAL*4 YC1,YE,YP
      REAL*4 ZC1,ZE,ZP

      REAL*8 CV
      REAL*8 STOT,VOL

      INTEGER*4  KSURF
      
      PARAMETER (MC=100000,MV=100000)
      
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
      COMMON/PUN/ITO(MV),ISO(MV),ISA(MV),XP(MV),YP(MV),ZP(MV),AP(MV)

C     KSURF= 0
C     ASS1=.FALSE.
C     RD=1.4E0
C     OFAC=0.80
C     RMIN=0.50E0
C     NDIV=3

      IF(KSURF.EQ.0)THEN

       WRITE(6,'(/A)')' The Van der Waals Surface is calculated'
       WRITE(6,'(A)')' ---------------------------------------'
       WRITE(6,'(A,I2/)')
     &  ' NDIV                                  =  ',NDIV

      ELSE IF(KSURF.EQ.1)THEN

C       WRITE(6,'(/A)')' The Accessible Surface is calculated'
C       WRITE(6,'(A)')' ------------------------------------'
C       WRITE(6,'(A,I2)')
C     &  ' NDIV                                  =  ',NDIV
C       WRITE(6,'(A,F10.5/)')
C     &  ' The radius of the solvent(RSOL) is    =',RD

      ELSE IF(KSURF.EQ.2)THEN

       WRITE(6,'(/A)')' The Solvent-Excluding Surface is calculated'
       WRITE(6,'(A)')' -------------------------------------------'
       WRITE(6,'(A,I2)')
     &  ' NDIV                                 =  ',NDIV
       WRITE(6,'(A,F10.5)')
     &  ' Radius of the solvent         (RSOL) =',RD
       WRITE(6,'(A,F10.5)')
     &  ' Minimum Radius for new sphere (RMIN) =',RMIN
       WRITE(6,'(A,F10.5/)')
     &  ' Overlapping factor            (OFAC) =',OFAC
 
      END IF
      IF (ASS1.AND.(KSURF.EQ.2))WRITE(6,'(A)')
     &' The new spheres will be assigned to the initials using ASSG1'
      
      NCOR=NATOM

C*****Make tesselation for sphere of radius 1.0*********

      CALL TES
      CALL DIVIDE(NDIV)

C*****Beginning  the calculations ******

C Create the new set of spheres if we are interested in the 
C Solvent-excluding Surface

      IF(KSURF.EQ.2)THEN

        IF(GHOST)CALL SHELL(NCOR,RD)
        CALL BULK(NATOM,NCOR,NDIV,OFAC,RD)
        CALL CLEAN5(NATOM,NCOR,NDIV)
        CALL CREA(NATOM,NCOR,RMIN,OFAC,RD,NDIV) 
        CALL CLEAN5(NATOM,NCOR,NDIV)
      END IF 

C
C  Compute the surface

      IF(KSURF.EQ.1)CALL SUM(NCOR,RD,'SUMA')
      CALL GEOCAV(NCOR,NP,NDIV,GHOST)
      IF(KSURF.EQ.2)THEN

         IF(ASS1)CALL ASSIGN1(NP,NATOM,GHOST)

      END IF
C
C Compute the area and volume
      CALL VOLARE(NP,GHOST,STOT,VOL)

      RETURN
      END
C 
      SUBROUTINE TES
C     --------------------------------------------------------------------
C     This computes the triangle vertex coordinates for a sphere of radius
C     one, projecting the pentakisdodecahedro onto it.
C     --------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER*4 I,II
      INTEGER*4 J

      REAL*4 XC1,YC1,ZC1

      REAL*8 CTH,CV
      REAL*8 FI,FIR,FIV
      REAL*8 STH
      REAL*8 TH,THEV
      
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
      
      DIMENSION THEV(6),FIV(6)
      
      DATA THEV/0.6523581397843682D0,1.1071487177940905D0,
     $          1.3820857960113345D0,1.7595068575784587D0,
     $          2.0344439357957027D0,2.4892345138054251D0/
      DATA FIV/0.6283185307179586D0,0.0 D0              ,
     $         0.6283185307179586D0,0.0 D0              ,
     $         0.6283185307179586D0,0.0 D0              /
      DATA FIR/1.2566370614359173 D0/

      CV(1,1)=0.D0
      CV(1,2)=0.D0
      CV(1,3)=1.D0
      CV(32,1)=0.D0
      CV(32,2)=0.D0
      CV(32,3)=-1.D0
      II=1
      DO 520 I=1,6
      TH=THEV(I)
      FI=FIV(I)
      CTH=DCOS(TH)
      STH=DSIN(TH)
      DO 521 J=1,5
      FI=FI+FIR
      IF(J.EQ.1) FI=FIV(I)
      II=II+1
      CV(II,1)=STH*DCOS(FI)
      CV(II,2)=STH*DSIN(FI)
      CV(II,3)=CTH
  521 CONTINUE
  520 CONTINUE
      RETURN
      END
C
      SUBROUTINE DIVIDE(NDIV)
C     ---------------------------------------------------------------
C     This divides the initial 60 spherical triangles to the level
C     indicated by NDIV
C     ---------------------------------------------------------------
      IMPLICIT NONE

      INTEGER*4 IJ
      INTEGER*4 J,J2,J3,J4,J5,JVT1,JVT2
      INTEGER*4 NDIV,NV1,NV2,NV21,NV22,NV23,NV3,NV31,NV32,NV33
      INTEGER*4 NV41,NV42,NV43,NV51,NV52,NV53
     
      REAL*4 XC1,YC1,ZC1
   
      REAL*8 CC,CV,CVN2,CVN3,CVN4,CVN5
      REAL*8 FOUR     
      REAL*8 PI
      REAL*8 XV1,XV2,XV3      
      REAL*8 YV1,YV2,YV3     
      REAL*8 ZERO,ZV1,ZV2,ZV3     

      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
      COMMON/PENTA/JVT1(3,60),JVT2(3,4)

      DIMENSION CVN2(6,3),CVN3(6,3),CVN4(6,3),CVN5(6,3),CC(3)
 
      DATA ZERO/0.0D0/
      DATA PI/3.1415926535897932D0/
      
      IJ=0
C*****Level 1****************************
      DO 10 J=1,60
        NV1=JVT1(1,J)
        NV2=JVT1(2,J)
        NV3=JVT1(3,J)
        XV1=CV(NV1,1)
        YV1=CV(NV1,2)
        ZV1=CV(NV1,3)
        XV2=CV(NV2,1)
        YV2=CV(NV2,2)
        ZV2=CV(NV2,3)
        XV3=CV(NV3,1)
        YV3=CV(NV3,2)
        ZV3=CV(NV3,3)
        IF(NDIV.GT.1) GO TO 20
        CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
        IJ=IJ+1
        XC1(IJ)=SNGL(CC(1))  ! save center of the initial triangles in XC1,YC1,ZC1
        YC1(IJ)=SNGL(CC(2))
        ZC1(IJ)=SNGL(CC(3))
        GO TO 10
C*****Level 2**********************
   20   CONTINUE
        CVN2(1,1)=XV1
        CVN2(1,2)=YV1
        CVN2(1,3)=ZV1
        CVN2(2,1)=XV2
        CVN2(2,2)=YV2
        CVN2(2,3)=ZV2
        CVN2(3,1)=XV3
        CVN2(3,2)=YV3
        CVN2(3,3)=ZV3
        CALL CALVER(CVN2)
        DO 21 J2=1,4
          NV21=JVT2(1,J2)
          NV22=JVT2(2,J2)
          NV23=JVT2(3,J2)
          XV1=CVN2(NV21,1)
          YV1=CVN2(NV21,2)
          ZV1=CVN2(NV21,3)
          XV2=CVN2(NV22,1)
          YV2=CVN2(NV22,2)
          ZV2=CVN2(NV22,3)
          XV3=CVN2(NV23,1)
          YV3=CVN2(NV23,2)
          ZV3=CVN2(NV23,3)
          IF(NDIV.GT.2) GO TO 30
          CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
          IJ=IJ+1
          XC1(IJ)=SNGL(CC(1))
          YC1(IJ)=SNGL(CC(2))
          ZC1(IJ)=SNGL(CC(3))
          GO TO 21
C*****Level 3**********************************
   30     CONTINUE
          CVN3(1,1)=XV1
          CVN3(1,2)=YV1
          CVN3(1,3)=ZV1
          CVN3(2,1)=XV2
          CVN3(2,2)=YV2
          CVN3(2,3)=ZV2
          CVN3(3,1)=XV3
          CVN3(3,2)=YV3
          CVN3(3,3)=ZV3
          CALL CALVER(CVN3)
          DO 31 J3=1,4
            NV31=JVT2(1,J3)
            NV32=JVT2(2,J3)
            NV33=JVT2(3,J3)
            XV1=CVN3(NV31,1)
            YV1=CVN3(NV31,2)
            ZV1=CVN3(NV31,3)
            XV2=CVN3(NV32,1)
            YV2=CVN3(NV32,2)
            ZV2=CVN3(NV32,3)
            XV3=CVN3(NV33,1)
            YV3=CVN3(NV33,2)
            ZV3=CVN3(NV33,3)
            IF(NDIV.GT.3) GO TO 40
            CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
            IJ=IJ+1
            XC1(IJ)=SNGL(CC(1))
            YC1(IJ)=SNGL(CC(2))
            ZC1(IJ)=SNGL(CC(3))
            GO TO 31
C*****Level 4******************************
   40       CONTINUE
            CVN4(1,1)=XV1
            CVN4(1,2)=YV1
            CVN4(1,3)=ZV1
            CVN4(2,1)=XV2
            CVN4(2,2)=YV2
            CVN4(2,3)=ZV2
            CVN4(3,1)=XV3
            CVN4(3,2)=YV3
            CVN4(3,3)=ZV3
            CALL CALVER(CVN4)
            DO 41 J4=1,4
              NV41=JVT2(1,J4)
              NV42=JVT2(2,J4)
              NV43=JVT2(3,J4)
              XV1=CVN4(NV41,1)
              YV1=CVN4(NV41,2)
              ZV1=CVN4(NV41,3)
              XV2=CVN4(NV42,1)
              YV2=CVN4(NV42,2)
              ZV2=CVN4(NV42,3)
              XV3=CVN4(NV43,1)
              YV3=CVN4(NV43,2)
              ZV3=CVN4(NV43,3)
              IF(NDIV.GT.4) GO TO 50
              CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
              IJ=IJ+1
              XC1(IJ)=SNGL(CC(1))
              YC1(IJ)=SNGL(CC(2))
              ZC1(IJ)=SNGL(CC(3))
              GO TO 41
C*****Level 5*************************************
   50         CONTINUE
              CVN5(1,1)=XV1
              CVN5(1,2)=YV1
              CVN5(1,3)=ZV1
              CVN5(2,1)=XV2
              CVN5(2,2)=YV2
              CVN5(2,3)=ZV2
              CVN5(3,1)=XV3
              CVN5(3,2)=YV3
              CVN5(3,3)=ZV3
              CALL CALVER(CVN5)
              DO 51 J5=1,4
                NV51=JVT2(1,J5)
                NV52=JVT2(2,J5)
                NV53=JVT2(3,J5)
                XV1=CVN5(NV51,1)
                YV1=CVN5(NV51,2)
                ZV1=CVN5(NV51,3)
                XV2=CVN5(NV52,1)
                YV2=CVN5(NV52,2)
                ZV2=CVN5(NV52,3)
                XV3=CVN5(NV53,1)
                YV3=CVN5(NV53,2)
                ZV3=CVN5(NV53,3)
                CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
                IJ=IJ+1
                XC1(IJ)=SNGL(CC(1))
                YC1(IJ)=SNGL(CC(2))
                ZC1(IJ)=SNGL(CC(3))
   51         CONTINUE
   41       CONTINUE
   31     CONTINUE
   21   CONTINUE
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE CALVER(CVN)
C     ---------------------------------------------------------------------
C     This divides one triangle into four.
C     ---------------------------------------------------------------------
C  AT the start CVN(I,J), I=1,2,3  J=1,2,3 - contain verticies of a triangle 
C  on a unit sphere  that is to be divided
C  points CVN(I,J), I=4,5,6 are first taken in the middle of bonds (3-1),(1-2) and (2-3)
C  and then scaled to be also put on a unit sphere
C
      IMPLICIT NONE
      REAL*8 CVN,XXX,YYY,ZZZ,RRR,FC
      INTEGER*4 N,N1,N2
      DIMENSION CVN(6,3)
      DO 7 N=1,3
        N2=N+3
        N1=N-1
        IF(N.EQ.1)N1=3
        XXX=(CVN(N,1)+CVN(N1,1))/2.0D0
        YYY=(CVN(N,2)+CVN(N1,2))/2.0D0
        ZZZ=(CVN(N,3)+CVN(N1,3))/2.0D0
        RRR=SQRT(XXX*XXX+YYY*YYY+ZZZ*ZZZ) 
        FC=1.0D0/RRR
        CVN(N2,1)=XXX*FC
        CVN(N2,2)=YYY*FC
        CVN(N2,3)=ZZZ*FC
    7 CONTINUE
      RETURN
      END
C
      SUBROUTINE CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
C     ---------------------------------------------------------------------
C     This computes the center of a spherical triangle.
C     ---------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC,XXX,YYY,ZZZ,RRR,FC
      DIMENSION CC(3)
      XXX=(XV1+XV2+XV3)/3.0D0
      YYY=(YV1+YV2+YV3)/3.0D0
      ZZZ=(ZV1+ZV2+ZV3)/3.0D0
      RRR=SQRT(XXX*XXX+YYY*YYY+ZZZ*ZZZ)
      FC=1.0D0/RRR
      CC(1)=XXX*FC
      CC(2)=YYY*FC
      CC(3)=ZZZ*FC
      RETURN
      END

C
      SUBROUTINE GEOCAV(NCOR,NP,NDIV,GHOST)
C     ------------------------------------------------------------------
C     This computes the surface.
C     ------------------------------------------------------------------
      IMPLICIT NONE

      LOGICAL GHOST

      INTEGER*2 IUSE

      INTEGER*4 I,IIJ,IJE,ITO,ISA,ISO
      INTEGER*4 J
      INTEGER*4 K 
      INTEGER*4 L1
      INTEGER*4 MC,MV
      INTEGER*4 N3,N4,NCOR,NDIV,NEJCI,NINF,NP,NSUP,NTRIAN,NTS
      INTEGER*4 UN3

      REAL*4 AP,ATP,ATS
      REAL*4 DD,DIJ2
      REAL*4 FC,FNDIV
      REAL*4 PI
      REAL*4 RE,REI,RREJ,RRR
      REAL*4 SRE2
      REAL*4 XC1,XE,XEI,XP,XPL,XSL,XSM
      REAL*4 YC1,YE,YEI,YP,YPL,YSL,YSM
      REAL*4 ZC1,ZE,ZEI,ZP,ZPL,ZSL,ZSM
      REAL*8 CV

      PARAMETER (MC=100000,MV=100000)
      
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
      COMMON/PUN/ITO(MV),ISO(MV),ISA(MV),XP(MV),YP(MV),ZP(MV),AP(MV)

      DIMENSION IJE(MC) ! will contain indexes of spheres ( tot number NEJCI) that overlap 
	                  ! with a given sphere

      DATA PI/3.141593E0/

      NP=0

C begin
      NTRIAN=4**(NDIV-1) ! the number of trianlges obtained by division of a main triangle
      FNDIV=60*NTRIAN    ! the total number of triangles on the sphere

C It selects one sphere
      DO 1 I=1,NCOR
      
        IF(IUSE(I).EQ.6) THEN ! if real sphere
          REI=RE(I)
          XEI=XE(I)
          YEI=YE(I)
          ZEI=ZE(I)
          ATS=4.0E0*PI*REI*REI/FNDIV  ! surface of a single axx triangle 

          IIJ=0
C
C  It determines which spheres are linked to sphere I
C  (Distance between centers of spheres is less then sum of radii
          DO J=1,NCOR

            IF(IUSE(J).GE.3) THEN

              IF(I.NE.J)THEN

                DIJ2=(XEI-XE(J))*(XEI-XE(J))+
     &               (YEI-YE(J))*(YEI-YE(J))+
     &               (ZEI-ZE(J))*(ZEI-ZE(J))
                SRE2=(REI+RE(J))*(REI+RE(J))

                IF(DIJ2.LT.SRE2) THEN
                  IIJ=IIJ+1
                  IJE(IIJ)=J
                  NEJCI=IIJ
                END IF

              END IF
            END IF
          END DO
C 
C It selects one main triangle.
          NSUP=0
          UN3=1
          DO 2 J=1,60
            XPL=0.0E0
            YPL=0.0E0
            ZPL=0.0E0
            NTS=0
            NINF=NSUP+1
            NSUP=NINF+NTRIAN-1
C
C It selects one secondary triangle.
            DO 3 K=NINF,NSUP
              XSL=XC1(K)*REI  
              YSL=YC1(K)*REI
              ZSL=ZC1(K)*REI
              XSM=XSL+XEI ! coordinates of the center of the secondary triangle 
              YSM=YSL+YEI
              ZSM=ZSL+ZEI
C
C It fixes if the secondary triangle is inside or outside.
              L1=UN3          ! set index sphere the same as for last sphere containing
	                        ! secondary triangle (for speed up apparently)
              DO N3=L1,NEJCI  ! cycle on overlapping spheres
                UN3=N3
                N4=IJE(N3)
                DD=(XSM-XE(N4))*(XSM-XE(N4))+
     &             (YSM-YE(N4))*(YSM-YE(N4))+
     &             (ZSM-ZE(N4))*(ZSM-ZE(N4))
                RREJ=RE(N4)*RE(N4)
                IF(DD.LT.RREJ) GO TO 3  ! if a secondary triangle inside one of the spheres
              END DO                    ! got to the next triangles

              DO N3=1,L1-1
                UN3=N3
                N4=IJE(N3)
                DD=(XSM-XE(N4))*(XSM-XE(N4))+
     &             (YSM-YE(N4))*(YSM-YE(N4))+
     &             (ZSM-ZE(N4))*(ZSM-ZE(N4))
                RREJ=RE(N4)*RE(N4)
                IF(DD.LT.RREJ) GO TO 3  ! skip a secondary triangle if inside one of the spheres
              END DO
C
C It prepares the coordinates for the main triangle
              XPL=XPL+XSL
              YPL=YPL+YSL
              ZPL=ZPL+ZSL
              NTS=NTS+1
    3       CONTINUE
C
C It reduces the secondary triangles to the main triangle.
            IF(NTS.EQ.0)GO TO 2 ! if free trianlges found 

            ATP=ATS*NTS
            XPL=XPL/NTS  ! average coordinates of the free triangles
            YPL=YPL/NTS
            ZPL=ZPL/NTS
            RRR=SQRT(XPL*XPL+YPL*YPL+ZPL*ZPL)
            FC=REI/RRR

            NP=NP+1
            XP(NP)=XPL*FC+XEI ! average position of the centers of secondary triangles
            YP(NP)=YPL*FC+YEI
            ZP(NP)=ZPL*FC+ZEI
            AP(NP)=ATP
            ITO(NP)=J 
            ISO(NP)=I
            ISA(NP)=I

    2     CONTINUE
        END IF
    1 CONTINUE

      RETURN
      END
C       
      SUBROUTINE ASSIGN1(NP,NATOM,GHOST)
C     --------------------------------------------------------------------
C     This assigns the tesserae of the new spheres to the initial ones.
C     --------------------------------------------------------------------
      IMPLICIT NONE

      LOGICAL GHOST

      INTEGER*2 IU,IUSE

      INTEGER*4 I,IASS,ISA,ISO,ITO
      INTEGER*4 J
      INTEGER*4 MC,MV
      INTEGER*4 NATOM,NI,NP

      REAL*4 AP
      REAL*4 DIS,DMI
      REAL*4 RE
      REAL*4 XE,XP
      REAL*4 YE,YP
      REAL*4 ZE,ZP

      PARAMETER (MC=100000,MV=100000)

      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
      COMMON/PUN/ITO(MV),ISO(MV),ISA(MV),XP(MV),YP(MV),ZP(MV),AP(MV)
      
      DO I=1,NP
 
           IF(ISO(I).GT.NATOM)THEN
                  DMI=99999.9
                  IASS=0

                  DO J=1,NATOM

                   IF(IUSE(J).NE.1)THEN
                   DIS= (XE(J)-XP(I))*(XE(J)-XP(I))+
     &                  (YE(J)-YP(I))*(YE(J)-YP(I))+
     &                  (ZE(J)-ZP(I))*(ZE(J)-ZP(I))
                   DIS=SQRT(DIS)-RE(J)

                     IF(DIS.LT.DMI)THEN
                     IASS=J
                     DMI=DIS
                     END IF

                   END IF

                  END DO

           ISA(I)=IASS
           END IF
      END DO
     
      IF(GHOST)THEN
        J=0
        DO I=1,NP
        NI=I-J
        ITO(NI)=ITO(I)
        ISO(NI)=ISO(I)
        ISA(NI)=ISA(I)
        XP(NI)=XP(I)
        YP(NI)=YP(I)
        ZP(NI)=ZP(I)
        AP(NI)=AP(I)
        IU=IUSE(ISA(I))        

         IF(IU.EQ.2)THEN
            IUSE(ISA(I))=6        
         ELSE IF(IU.EQ.3)THEN 
            J=J+1
         ELSE IF(IU.EQ.4)THEN 
            IUSE(ISA(I))=6        
         ELSE IF(IU.EQ.5)THEN 
            J=J+1
         END IF

        END DO
        NP=NP-J
      END IF

      RETURN
      END

C
      SUBROUTINE SHELL(NCOR,RD)
C     --------------------------------------------------------------------
C     This determines which ghost spheres are around the real ones.
C     --------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER*2 IUSE

      INTEGER*4 I
      INTEGER*4 J
      INTEGER*4 MC
      INTEGER*4 NCOR

      REAL*4 DD
      REAL*4 RD,RE
      REAL*4 TEST
      REAL*4 XE
      REAL*4 YE
      REAL*4 ZE

      PARAMETER (MC=100000)

      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)

      DO I=1,NCOR-1
      DO J=I+1,NCOR

        IF((IUSE(I).EQ.3).AND.(IUSE(J).EQ.6))THEN

              DD = (XE(I)-XE(J)) * (XE(I)-XE(J)) +
     &             (YE(I)-YE(J)) * (YE(I)-YE(J)) +
     &             (ZE(I)-ZE(J)) * (ZE(I)-ZE(J))
              TEST=(RE(I)+RE(J)+2*RD)*(RE(I)+RE(J)+2*RD)
              IF(DD.LT.TEST)THEN
              IUSE(I)=5
              END IF

        ELSE IF((IUSE(I).EQ.6).AND.(IUSE(J).EQ.3))THEN

              DD = (XE(I)-XE(J)) * (XE(I)-XE(J)) +
     &             (YE(I)-YE(J)) * (YE(I)-YE(J)) +
     &             (ZE(I)-ZE(J)) * (ZE(I)-ZE(J))
              TEST=(RE(I)+RE(J)+2*RD)*(RE(I)+RE(J)+2*RD)

              IF(DD.LT.TEST)THEN
              IUSE(J)=5
              END IF

        END IF

      END DO
      END DO

      RETURN 
      END
C
      SUBROUTINE VOLARE(NP,GHOST,STOT,VOL)
C     ---------------------------------------------------------------------
C     This calculates total area and volume.
C     ---------------------------------------------------------------------

      IMPLICIT NONE    

      LOGICAL GHOST

      INTEGER*2 IUSE

      INTEGER*4 I,ISA,ISO,ITO
      INTEGER*4 J,JU
      INTEGER*4 NP
      INTEGER*4 MV,MC

      REAL*4 AP
      REAL*4 DP
      REAL*4 RE,REJ
      REAL*4 VD,VN
      REAL*4 XE,XP,XEJ
      REAL*4 YE,YP,YEJ
      REAL*4 ZE,ZP,ZEJ

      REAL*8 STOT,VOL

      PARAMETER (MC=100000,MV=100000)

      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
      COMMON/PUN/ITO(MV),ISO(MV),ISA(MV),XP(MV),YP(MV),ZP(MV),AP(MV)

      STOT=0.0E0
      VOL=0.0E0

      DO I=1,NP
      STOT=STOT+AP(I)
      END DO


      IF(.NOT.GHOST)THEN
      
        JU=0
        DO I=1,NP
        J=ISO(I)

           IF(JU.NE.J)THEN
            JU=J
            XEJ=XE(J)
            YEJ=YE(J)
            ZEJ=ZE(J)
            REJ=RE(J)
           END IF

        VOL=VOL+((XP(I)-XEJ)*XP(I)+
     &           (YP(I)-YEJ)*YP(I)+
     &           (ZP(I)-ZEJ)*ZP(I)) *AP(I)/(3.0E0*REJ)


        END DO

      END IF


C      WRITE(6,'(//A)')' ----------------- RESULTS -----------------'
C      WRITE(6,'(A/)') ' -------------------------------------------'
C      WRITE(6,'(A)')  ' *******************************************'
C      WRITE(6,'(A)')  ' *******************************************'
C      WRITE(6,'(A,F17.3,A)')' ** Area             =',STOT,'   **'
C      IF(.NOT.GHOST)THEN
C      WRITE(6,'(A,F17.3,A)')' ** Volume           =',VOL,'   **'
C      END IF
C      WRITE(6,'(A,I11,A)')' ** Number of Points =  ',NP,'       **'
C      WRITE(6,'(A)')' *******************************************'
C      WRITE(6,'(A)')' *******************************************'

      RETURN 
      END
C
      SUBROUTINE SUM(NCOR,RD,OP)
C     --------------------------------------------------------------------
C     Add the solvent radius to every sphere raddi
C     Or subtract the solvent radius.
C     --------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER*2 IUSE

      INTEGER*4 I
      INTEGER*4 MC
      INTEGER*4 NCOR

      REAL*4 RD,RE
      REAL*4 XE
      REAL*4 YE
      REAL*4 ZE
      
      CHARACTER*4 OP
 
      PARAMETER (MC=100000)

      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)     

      IF(OP.EQ.'SUMA')THEN

        DO I=1,NCOR
         IF(IUSE(I).NE.1)RE(I)=RE(I)+RD
        END DO

      ELSE IF(OP.EQ.'REST')THEN

        DO I=1,NCOR
         IF(IUSE(I).NE.1)RE(I)=RE(I)-RD
        END DO

      ELSE

        WRITE(6,'(A)')' Variable OP badly defined'
        STOP

      END IF
 
      RETURN
      END

C
      SUBROUTINE CREA(NATOM,NCOR,RMIN,OFAC,RD,NDIV)
C     ----------------------------------------------------------------
C     This subroutine creates the new spheres
C     ----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER*2 IUSE
  
      INTEGER*4 I
      INTEGER*4 J
      INTEGER*4 K,KK,KG,KP,KENG
      INTEGER*4 LK
      INTEGER*4 MC
      INTEGER*4 NATOM,NCOR,NEK,NGE,NL1,NL2,NDIV,NENG
      
      REAL*4 DI
      REAL*4 FC,FC1
      REAL*4 OFAC,OFACT
      REAL*4 RD,RE,REK,REN,REG,REG2,REGD2,REND2A,REND2C
      REAL*4 REP,REP2,REPD2,RMIN,RMID2,RGN
      REAL*4 RIJ,RIJ2,RIK,RIN,RNK2
      REAL*4 TEST
      REAL*4 XE,XEN
      REAL*4 YE,YEN
      REAL*4 ZE,ZEN

      PARAMETER (MC=100000)
      
      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
 
      DIMENSION DI(MC),KENG(MC)

      OFACT=1.0E0-OFAC*2.0E0
      RMID2=(RMIN+RD)*(RMIN+RD)
      NGE=0

      NL1=2   
      NL2=NCOR

  600 CONTINUE

C Loop to select the first sphere of the pair.
      DO 602 I=NL1,NL2
      IF(IUSE(I).LE.4) GO TO 602

      DO  J=1,NCOR
          DI(J)= SQRT( (XE(I)-XE(J)) * (XE(I)-XE(J)) +
     &                 (YE(I)-YE(J)) * (YE(I)-YE(J)) +
     &                 (ZE(I)-ZE(J)) * (ZE(I)-ZE(J)) )
      END DO

C
C Loop to select the second sphere of the pair.

      DO 603 J=1,I-1
      IF(IUSE(J).LE.4)GO TO 603

      RIJ= DI(J)

C If the solvent can pass through the pair this pair is discarded
      TEST=RD+RD+RE(I)+RE(J)
      IF(RIJ.GE.TEST) GO TO 603


C Determine which is the largest and smallest sphere
      IF(RE(I).GT.RE(J))THEN
        REG=RE(I)
        REP=RE(J)
        KG=I
        KP=J
      ELSE
        REG=RE(J)
        REP=RE(I)
        KG=J
        KP=I
      END IF

      REG2=REG*REG
      REP2=REP*REP
      REGD2=(REG+RD)*(REG+RD)

C Determine whether the spheres are overlapped
      IF(RIJ.LE.(REP+REG))THEN

C Test of overlapping
         TEST=REG+REP*OFACT
         IF(RIJ.LT.TEST)GO TO 603

C Test of the small sphere
         RIJ2=RIJ*RIJ
         REPD2=(REP+RD)*(REP+RD)
         RGN=(RIJ-REP+REG)*0.5E0
         REND2A=REGD2+RGN*(RGN-(REGD2+RIJ2-REPD2)/RIJ)
         IF(REND2A.LE.RMID2) GO TO 603

C Spheres A
         FC=(RIJ-REP+REG)/(RIJ+REP-REG)
         FC1=FC+1.0E0
         XEN=(XE(KG)+FC*XE(KP))/FC1
         YEN=(YE(KG)+FC*YE(KP))/FC1
         ZEN=(ZE(KG)+FC*ZE(KP))/FC1
         REN=SQRT(REND2A)-RD
         RIN=(RIJ-RE(J)+RE(I))*0.5E0

      ELSE

C Test of the small sphere
         RIJ2=RIJ*RIJ
         REPD2=(REP+RD)*(REP+RD)
         REND2C=REGD2+REG2-(REG/RIJ)*(REGD2+RIJ2-REPD2)
         IF(REND2C.LE.RMID2) GO TO 603

C Calculate radius for sphere of kind B and separate B and C
         RGN=(RIJ-REP+REG)*0.5E0
         REND2A=REGD2+RGN*(RGN-(REGD2+RIJ2-REPD2)/RIJ)

         IF(REND2A.GT.RMID2)THEN

C Spheres B
            FC=(RIJ-REP+REG)/(RIJ+REP-REG)
            FC1=FC+1.0E0
            XEN=(XE(KG)+FC*XE(KP))/FC1
            YEN=(YE(KG)+FC*YE(KP))/FC1
            ZEN=(ZE(KG)+FC*ZE(KP))/FC1
            REN=SQRT(REND2A)-RD
            RIN=(RIJ-RE(J)+RE(I))*0.5E0

         ELSE

C Spheres C
            FC=REG/(RIJ-REG)
            FC1=FC+1.0E0
            XEN=(XE(KG)+FC*XE(KP))/FC1
            YEN=(YE(KG)+FC*YE(KP))/FC1
            ZEN=(ZE(KG)+FC*ZE(KP))/FC1
            REN=SQRT(REND2C)-RD
            IF(KG.EQ.I)THEN
               RIN=REG
            ELSE
               RIN=(RIJ-REG)
            END IF

         END IF

      END IF

C
C Test of overlapping for the new sphere
      NENG=0
      DO 604 K=1,NCOR

      RIK=DI(K)

      IF(RIK.GE.(RIN+REN+RE(K))) GO TO 604

      IF(IUSE(K).LE.3) GO TO 604

      RNK2 = (XEN-XE(K)) * (XEN-XE(K)) +
     &       (YEN-YE(K)) * (YEN-YE(K)) +
     &       (ZEN-ZE(K)) * (ZEN-ZE(K))
      
      REK=RE(K)

      IF(RNK2.GE.((REK+REN)*(REK+REN))) GO TO 604
      
      TEST=(REN-REK)*(REN-REK)
      IF(RNK2.LE.TEST)THEN

           IF(REN.GT.REK)THEN
              IF(TEST.LT.4.0E-4) GO TO 603
              NENG=NENG+1
              KENG(NENG)=K
              GO TO 604
           ELSE
              GO TO 603
           END IF

      END IF

      TEST=REK+REN*OFACT
      IF(TEST.LT.0.0E0) GO TO 604
      TEST=TEST*TEST
      IF(RNK2.LE.TEST)GO TO 603

  604 CONTINUE

C
C Mark spheres engulfed by the new sphere
C
         DO LK=1,NENG
         K=KENG(LK)
         IUSE(K)=2
         END DO
       
C
C Creates the new spheres
C
         NCOR=NCOR+1
         XE(NCOR)=XEN
         YE(NCOR)=YEN
         ZE(NCOR)=ZEN
         RE(NCOR)=REN
         IUSE(NCOR)=6
         DI(NCOR)=RIN

  603 CONTINUE

  602 CONTINUE


      IF(NCOR.GT.MC)THEN
        WRITE(6,'(A)')' %-ERROR-% DIMENSION. TOO MANY SPHERES CREATED'
        STOP
      END IF

      NGE=NGE+1

      IF(NCOR.NE.NL2) THEN 
C Mark spheres with area zero
      CALL MZERO5(NATOM,NCOR,NL2,NDIV)
        IF(NCOR.NE.NL2) THEN 
           NL1=NL2+1
           NL2=NCOR
           GO TO 600
        END IF
      END IF
      
C      WRITE(6,'(A,I10)')' Number of generations     =',NGE
     
      RETURN
      END
C
      SUBROUTINE BULK(NATOM,NCOR,NDIV,OFAC,RD)
C     ----------------------------------------------------------------
C     This subroutine creates new spheres 
C     ----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER*2 IUSE
  
      INTEGER*4 I
      INTEGER*4 J
      INTEGER*4 K,KK,KENG
      INTEGER*4 L,L1,LK,LL
      INTEGER*4 MC
      INTEGER*4 NATOM,NCOR,NEK,NGE,NDIV,NL1,NL2,NTRI,NENG
      INTEGER*4 SC
      INTEGER*4 ULL
      
      REAL*4 DD,DI
      REAL*4 FC,FC1
      REAL*4 OFAC,OFACT
      REAL*4 RD,RE,REI,REJ,REK,REN,RNK2
      REAL*4 RIJ,RIK,RIN
      REAL*4 TEST
      REAL*4 XC1,XE,XEN,XPN
      REAL*4 YC1,YE,YEN,YPN
      REAL*4 ZC1,ZE,ZEN,ZPN

      REAL*8 CV

      PARAMETER (MC=100000)
      
      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
 
      DIMENSION DI(MC),SC(MC),KENG(MC)
     
      OFACT=1.0E0-OFAC*2.0E0

      NTRI=60*4**(NDIV-1)
      NGE=0

      NL1=2   
      NL2=NCOR

  600 CONTINUE
C
C Loop to select the first sphere of the pair.
      DO 602 I=NL1,NL2
      IF(IUSE(I).LE.4) GO TO 602
      REI=RE(I)

      DO  J=1,NCOR
          DI(J)=SQRT( (XE(I)-XE(J)) * (XE(I)-XE(J)) +
     &                (YE(I)-YE(J)) * (YE(I)-YE(J)) +
     &                (ZE(I)-ZE(J)) * (ZE(I)-ZE(J)) )
      END DO
C
C Loop to select the second sphere of the pair.

      DO 603 J=1,I-1
      IF(IUSE(J).LE.4)GO TO 603
      REJ=RE(J)
      RIJ= DI(J)
C
C If the solvent can pass through the pair this pair is discarded

      IF(RIJ.GE.(RD+RD+REI+REJ)) GO TO 603
C
C Test of overlapping and assignment of the radius of the new sphere

      IF(REI.GT.REJ)THEN
          TEST=REI+REJ*OFACT
          IF(RIJ.LT.TEST)GO TO 603
          REN=REI
      ELSE
          TEST=REJ+REI*OFACT
          IF(RIJ.LT.TEST)GO TO 603
          REN=REJ
      END IF
C
C Computes coordinates of the new sphere 

      FC=(RIJ-REJ+REI)/(RIJ+REJ-REI)
      FC1=FC+1.0E0
      XEN=(XE(I)+FC*XE(J))/FC1
      YEN=(YE(I)+FC*YE(J))/FC1
      ZEN=(ZE(I)+FC*ZE(J))/FC1
      RIN=(RIJ-REJ+REI)*0.5E0

C
C Test of overlapping for the new sphere

      NENG=0
      DO 604 K=1,NCOR

      RIK=DI(K)

      IF(RIK.GE.(RIN+REN+RE(K))) GO TO 604

      IF(IUSE(K).LE.3) GO TO 604

      RNK2 = (XEN-XE(K)) * (XEN-XE(K)) +
     &       (YEN-YE(K)) * (YEN-YE(K)) +
     &       (ZEN-ZE(K)) * (ZEN-ZE(K))
      
      REK=RE(K)

      IF(RNK2.GE.((REK+REN)*(REK+REN))) GO TO 604
      
      TEST=(REN-REK)*(REN-REK)
      IF(RNK2.LE.TEST)THEN

           IF(REN.GT.REK)THEN
              IF(TEST.LT.4.0E-4) GO TO 603
              NENG=NENG+1
              KENG(NENG)=K
              GO TO 604
           ELSE
              GO TO 603
           END IF

      END IF

      TEST=REK+REN*OFACT
      IF(TEST.LT.0.0E0) GO TO 604
      TEST=TEST*TEST
      IF(RNK2.LE.TEST)GO TO 603

  604 CONTINUE
C
C Find spheres that are overlapped with the new sphere
C Use radius of the sphere plus RD
C
      NEK=0
      DO K=1,NATOM

         IF(IUSE(K).GE.2)THEN
            DD= (XEN-XE(K))*(XEN-XE(K)) +
     &          (YEN-YE(K))*(YEN-YE(K)) + 
     &          (ZEN-ZE(K))*(ZEN-ZE(K))

            TEST=(REN+RD+RD+RE(K))*(REN+RD+RD+RE(K))
               
                IF(DD.LT.TEST)THEN
                  NEK=NEK+1
                  SC(NEK)=K
                END IF

         END IF

      END DO


C
C Determine if the new sphere has accessible surface area zero.

      ULL=1
      DO 3 L=1,NTRI
      XPN=XC1(L)*(REN+RD)+XEN
      YPN=YC1(L)*(REN+RD)+YEN
      ZPN=ZC1(L)*(REN+RD)+ZEN

            L1=ULL
            DO LL=L1,NEK
            ULL=LL
            K=SC(LL)
            DD=(XPN-XE(K))*(XPN-XE(K)) +
     &         (YPN-YE(K))*(YPN-YE(K)) + 
     &         (ZPN-ZE(K))*(ZPN-ZE(K)) 
            TEST=(RE(K)+RD)*(RE(K)+RD)
            IF(DD.LT.TEST)GO TO 3
            END DO

            DO LL=1,L1-1
            ULL=LL
            K=SC(LL)
            DD=(XPN-XE(K))*(XPN-XE(K)) +
     &         (YPN-YE(K))*(YPN-YE(K)) + 
     &         (ZPN-ZE(K))*(ZPN-ZE(K)) 
            TEST=(RE(K)+RD)*(RE(K)+RD)
            IF(DD.LT.TEST)GO TO 3
            END DO
            
      GO TO 603
      
    3 CONTINUE

C
C Marks spheres that are engulfed by the new sphere
C
      DO LK=1,NENG
      K=KENG(LK)
      IUSE(K)=2
      END DO

C
C Save information about the new sphere
C
      NCOR=NCOR+1
      XE(NCOR)=XEN
      YE(NCOR)=YEN
      ZE(NCOR)=ZEN
      RE(NCOR)=REN
      IUSE(NCOR)=6
      DI(NCOR)=RIN

  603 CONTINUE

  602 CONTINUE

      IF(NCOR.GT.MC)THEN
        WRITE(6,'(A)')' %-ERROR-% DIMENSION. TOO MANY SPHERES CREATED'
        STOP
      END IF


C Compute number of generations
      NGE=NGE+1

C Check if there are new spheres
      IF(NCOR.NE.NL2) THEN 
C Mark spheres with area zero
        CALL MZERO5(NATOM,NCOR,NL2,NDIV)
        IF(NCOR.NE.NL2) THEN 
           NL1=NL2+1
           NL2=NCOR
           GO TO 600
        END IF
      END IF

      RETURN
      END
C
      SUBROUTINE MZERO5(NATOM,NCOR,NL2,NDIV)
C     ----------------------------------------------------------------
C     This discards spheres(IUSE=2) engulfed by another
C     Marks (IUSE=4) the spheres with total area zero.
C     ----------------------------------------------------------------
      IMPLICIT NONE
    
      LOGICAL CHECK

      INTEGER*2 IUSE
  
      INTEGER*4 I,J,K
      INTEGER*4 MC
      INTEGER*4 L,LL,L1
      INTEGER*4 NATOM,NCOR,NDIV,NEJ,NEW,NI,NL2,NTRI,NZERO
      INTEGER*4 SC
      INTEGER*4 ULL

      REAL*4 DD
      REAL*4 RE,REI
      REAL*4 TEST
      REAL*4 XE,XEI,XC1,XP
      REAL*4 YE,YEI,YC1,YP
      REAL*4 ZE,ZEI,ZC1,ZP

      REAL*8 CV

      PARAMETER (MC=100000)

      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)


      DIMENSION SC(MC),CHECK(MC)

C
C Discard Spheres that have been engulfed by another (iuse=2)
      J=0
      K=0
      DO I=NATOM+1,NCOR
      NI=I-J
      XE(NI)=XE(I)
      YE(NI)=YE(I)
      ZE(NI)=ZE(I)
      RE(NI)=RE(I)
      IUSE(NI)=IUSE(I)

        IF(IUSE(I).EQ.2) THEN
           J=J+1
           IF(I.LE.NL2)THEN
              K=K+1
           END IF
        END IF

      END DO

      NEW=NCOR-NL2
      NCOR=NCOR-J
      NL2=NL2-K

C Print information about spheres discarded
c      WRITE(6,'(A,I7)')' New spheres in last generation       =',NEW
      NEW=NEW-(J-K)
c      WRITE(6,'(A,I7)')' Final new spheres in last generation =',NEW
c      WRITE(6,'(A,I7)')' Spheres discarded (engulfed)         =',J
     

C Find spheres with zero area.
      NZERO=0
      NTRI=60*4**(NDIV-1)

      DO I=1,NL2
      CHECK(I)=.FALSE.
      END DO

      DO I=NL2+1,NCOR
      CHECK(I)=.FALSE.
      IF(IUSE(I).GT.4)CHECK(I)=.TRUE.
      END DO

C
C Start
C
      DO 5 I=NCOR,1,-1
      IF(CHECK(I))THEN

      XEI=XE(I)
      YEI=YE(I)
      ZEI=ZE(I)
      REI=RE(I)

C Find spheres that overlap sphere I.
C  
          NEJ=0
          DO 7 J=1,NCOR
          IF(IUSE(J).GE.3)THEN

              DD=(XEI-XE(J))*(XEI-XE(J)) +
     &           (YEI-YE(J))*(YEI-YE(J)) +
     &           (ZEI-ZE(J))*(ZEI-ZE(J)) 
              TEST=(REI+RE(J))*(REI+RE(J))

              IF(DD.LT.TEST)THEN
              IF(I.NE.J)THEN

                  NEJ=NEJ+1
                  SC(NEJ)=J

              END IF
              END IF

          END IF
     
    7     CONTINUE
C          
C Mark spheres that should be checked
C
          IF(I.GT.NL2)THEN
          DO LL=1,NEJ
          J=SC(LL)
          IF(IUSE(J).GT.4)CHECK(J)=.TRUE.
          END DO
          END IF

C
C Determine if the sphere I has area zero
          ULL=1
          DO 4 L=1,NTRI
          XP=XC1(L)*REI+XEI
          YP=YC1(L)*REI+YEI
          ZP=ZC1(L)*REI+ZEI

             L1=ULL
             DO LL=L1,NEJ 
             ULL=LL
             J=SC(LL)
             DD=(XP-XE(J))*(XP-XE(J)) +
     &          (YP-YE(J))*(YP-YE(J)) +
     &          (ZP-ZE(J))*(ZP-ZE(J)) 
             TEST=RE(J)*RE(J)
             IF(DD.LT.TEST) GO TO 4
             END DO

             DO LL=1,L1-1 
             ULL=LL
             J=SC(LL)
             DD=(XP-XE(J))*(XP-XE(J)) +
     &          (YP-YE(J))*(YP-YE(J)) +
     &          (ZP-ZE(J))*(ZP-ZE(J)) 
             TEST=RE(J)*RE(J)
             IF(DD.LT.TEST) GO TO 4
             END DO
          
             GO TO 5

    4     CONTINUE

      IUSE(I)=4
      NZERO=NZERO+1

      END IF
    5 CONTINUE

c      WRITE(6,'(A,I7)')' Spheres marked with zero area        =',NZERO

      RETURN
      END
C

C
      SUBROUTINE CLEAN5(NATOM,NCOR,NDIV)
C     ----------------------------------------------------------------
C     This discard the spheres of area zero that are not needed for
C     the computation of the surface.
C     ----------------------------------------------------------------
      IMPLICIT NONE
    
      LOGICAL USEFUL

      INTEGER*2 IUSE
  
      INTEGER*4 I,J
      INTEGER*4 MC
      INTEGER*4 L,LJ,L1
      INTEGER*4 NATOM,NCOR,NDIV,NEJ,NI,NTRI
      INTEGER*4 SC
      INTEGER*4 ULJ

      REAL*4 DD
      REAL*4 RE,REI
      REAL*4 TEST
      REAL*4 XE,XEI,XC1,XP
      REAL*4 YE,YEI,YC1,YP
      REAL*4 ZE,ZEI,ZC1,ZP

      REAL*8 CV

      PARAMETER (MC=100000)

      COMMON/CSFE/XE(MC),YE(MC),ZE(MC),RE(MC),IUSE(MC)
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)

      DIMENSION SC(MC),USEFUL(MC)

      NTRI=60*4**(NDIV-1)


      DO I=1,NATOM
         USEFUL(I)=.TRUE.   
      END DO

      DO I=NATOM+1,NCOR
          IF(IUSE(I).NE.4)THEN
            USEFUL(I)=.TRUE.   
          ELSE
            USEFUL(I)=.FALSE.   
          END IF
      END DO

      DO I=1,NCOR
      IF(IUSE(I).GE.4)THEN      
      XEI=XE(I)
      YEI=YE(I)
      ZEI=ZE(I)
      REI=RE(I)

C Find spheres that overlap sphere I.
          NEJ=0
          DO J=1,NCOR
          IF(IUSE(J).GE.3)THEN

              DD=(XEI-XE(J))*(XEI-XE(J)) +
     &           (YEI-YE(J))*(YEI-YE(J)) +
     &           (ZEI-ZE(J))*(ZEI-ZE(J)) 
              TEST=(REI+RE(J))*(REI+RE(J))

              IF(DD.LT.TEST)THEN
              IF(I.NE.J)THEN
                  NEJ=NEJ+1
                  SC(NEJ)=J
              END IF
              END IF

          END IF
          END DO

C Find spheres that are needed to discard the triangle L
          ULJ=1
          DO 40 L=1,NTRI
          XP=XC1(L)*REI+XEI
          YP=YC1(L)*REI+YEI
          ZP=ZC1(L)*REI+ZEI

C Among the useful spheres
             L1=ULJ
             DO LJ=L1,NEJ
             ULJ=LJ
             J=SC(LJ)
             IF(USEFUL(J))THEN
             DD=(XP-XE(J))*(XP-XE(J)) +
     &          (YP-YE(J))*(YP-YE(J)) +
     &          (ZP-ZE(J))*(ZP-ZE(J)) 
             TEST=RE(J)*RE(J)

                IF(DD.LT.TEST) GO TO 40

             END IF
             END DO

             DO LJ=1,L1-1
             ULJ=LJ
             J=SC(LJ)
             IF(USEFUL(J))THEN
             DD=(XP-XE(J))*(XP-XE(J)) +
     &          (YP-YE(J))*(YP-YE(J)) +
     &          (ZP-ZE(J))*(ZP-ZE(J)) 
             TEST=RE(J)*RE(J)

                IF(DD.LT.TEST) GO TO 40

             END IF
             END DO

C Among the not useful
             DO LJ=1,NEJ
             ULJ=LJ
             J=SC(LJ)
             IF(.NOT.USEFUL(J))THEN
             DD=(XP-XE(J))*(XP-XE(J)) +
     &          (YP-YE(J))*(YP-YE(J)) +
     &          (ZP-ZE(J))*(ZP-ZE(J)) 
             TEST=RE(J)*RE(J)

                IF(DD.LT.TEST)THEN
                USEFUL(J)=.TRUE.
                GO TO 40
                END IF
             END IF
             END DO
   40     CONTINUE
      END IF
      END DO

C Discard spheres totally inside of others
      J=0
      DO I=NATOM+1,NCOR
      NI=I-J
      XE(NI)=XE(I)
      YE(NI)=YE(I)
      ZE(NI)=ZE(I)
      RE(NI)=RE(I)
      IUSE(NI)=IUSE(I)

        IF(.NOT.USEFUL(I)) THEN
           J=J+1
        END IF

      END DO

      NCOR=NCOR-J
    
c      WRITE(6,'(A,I7)')' Spheres discarded (not useful)       =',J

      RETURN
      END
