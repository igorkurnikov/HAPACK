* *******************************************************************
* COPYRIGHT (c) 1975 Hyprotech UK
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of an HSL ARCHIVE
* Licence, see http://hsl.rl.ac.uk/archive/cou.html
*
* Please note that for an HSL ARCHIVE Licence:
*
* 1. The Package must not be copied for use by any other person.
*    Supply of any part of the library by the Licensee to a third party
*    shall be subject to prior written agreement between AEA
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
*######DATE 9 February 1994
C       Toolpack tool decs employed.
C       SAVE statements added.
C       Arg dimensions made *.
C	No. of sig. figs. reduced to 15 for printed output.
C
      SUBROUTINE VA13AD(IPTR,FUNC,N,X,F,G,SCALE,ACC,W)
      DOUBLE PRECISION ACC,F
      INTEGER N
      DOUBLE PRECISION G(*),SCALE(*),W(*),X(*)
      EXTERNAL FUNC
      INTEGER ND,NGA,NGB,NW,NXA,NXB
      EXTERNAL VA13CD,VA13DD
      ND = 1 + (N* (N+1))/2
      NW = ND + N
      NXA = NW + N
      NGA = NXA + N
      NXB = NGA + N
      NGB = NXB + N
      CALL VA13CD(IPTR,FUNC,N,X,F,G,SCALE,ACC,W,W(ND),W(NW),W(NXA),
     +            W(NGA),W(NXB),W(NGB))
      RETURN
      END
      BLOCK DATA VA13DD
      COMMON /VA13BD/IPRINT,LP,MAXFUN,MODE,NFUN
      INTEGER IPRINT,LP,MAXFUN,MODE,NFUN
      SAVE /VA13BD/
      DATA IPRINT/0/,LP/6/,MAXFUN/0/,MODE/1/
      END
      SUBROUTINE VA13CD(IPTR,FUNC,N,X,F,G,SCALE,ACC,H,D,W,XA,GA,XB,GB)
      DOUBLE PRECISION ACC,F
      INTEGER N
      DOUBLE PRECISION D(*),G(*),GA(*),GB(*),H(*),SCALE(*),W(*),X(*),
     +                 XA(*),XB(*)
      EXTERNAL FUNC
      DOUBLE PRECISION C,CC,DFF,DGA,DGB,FA,FB,FMIN,GL1,GL2,GMIN,STEP,
     +                 STEPBD,STEPLB,STMIN
      INTEGER I,IP,IPRA,IR,ISFV,ITR,K,NP
      EXTERNAL MC11AD,MC11BD,MC11ED
      INTRINSIC DABS,DMAX1,DMIN1,DSQRT,IABS,MIN0
      COMMON /VA13BD/IPRINT,LP,MAXFUN,MODE,NFUN
      INTEGER IPRINT,LP,MAXFUN,MODE,NFUN
      SAVE /VA13BD/
      IF (IPRINT.EQ.0) GO TO 10
      WRITE (LP,FMT=1000)
 1000 FORMAT ('1ENTRY TO VA13AD')
   10 CALL FUNC(IPTR,N,X,F,G)
      NFUN = 1
      ITR = 0
      NP = N + 1
      IF (MODE.GE.2) GO TO 60
   20 C = 0.D0
      DO 30 I = 1,N
   30 C = DMAX1(C,DABS(G(I)*SCALE(I)))
      IF (C.LE.0.D0) C = 1.D0
      K = (N*NP)/2
      DO 40 I = 1,K
   40 H(I) = 0.D0
      K = 1
      DO 50 I = 1,N
        H(K) = 0.01D0*C/SCALE(I)**2
   50 K = K + NP - I
      GO TO 100
   60 IF (MODE.GE.3) GO TO 80
      CALL MC11BD(H,N,K)
      IF (K.GE.N) GO TO 100
   70 WRITE (LP,FMT=1010)
 1010 FORMAT (/,5X,'BECAUSE THE HESSIAN GIVEN TO VA13AD IS NOT POS DEF,'
     +       ,/,5X,'IT HAS BEEN REPLACED BY A POSITIVE DIAGONAL MATRIX')
      GO TO 20
   80 K = 1
      DO 90 I = 1,N
        IF (H(K).LE.0D0) GO TO 70
   90 K = K + NP - I
  100 DFF = 0.D0
      IPRA = IABS(IPRINT)
      IP = IABS(IPRA-1)
  110 FA = F
      ISFV = 1
      DO 120 I = 1,N
        XA(I) = X(I)
  120 GA(I) = G(I)
  130 IP = IP + 1
      IF (IP.NE.IPRA) GO TO 140
      IP = 0
      WRITE (LP,FMT=1020) ITR,NFUN
 1020 FORMAT (/,5X,'ITERATION =',I5,5X,'FUNCTIONS =',I5)
      WRITE (LP,FMT=1030) FA
 1030 FORMAT (5X,'F =',D24.15)
      IF (IPRINT.LE.0) GO TO 140
      WRITE (LP,FMT=1040) (XA(I),I=1,N)
 1040 FORMAT (5X,'X(.) =', (5D24.15))
      WRITE (LP,FMT=1050) (GA(I),I=1,N)
 1050 FORMAT (5X,'G(.) =', (5D24.15))
  140 ITR = ITR + 1
      DO 150 I = 1,N
  150 D(I) = -GA(I)
      CALL MC11ED(H,N,D,W,N)
      C = 0.D0
      DGA = 0.D0
      DO 160 I = 1,N
        C = DMAX1(C,DABS(D(I)/SCALE(I)))
  160 DGA = DGA + GA(I)*D(I)
      IF (DGA.GE.0.D0) GO TO 240
      STMIN = 0.D0
      STEPBD = 0.D0
      STEPLB = ACC/C
      FMIN = FA
      GMIN = DGA
      STEP = 1.D0
      IF (DFF.LE.0D0) STEP = DMIN1(STEP,1D0/C)
      IF (DFF.GT.0D0) STEP = DMIN1(STEP, (DFF+DFF)/ (-DGA))
  170 C = STMIN + STEP
      IF (NFUN.EQ.MAXFUN) GO TO 250
      NFUN = NFUN + 1
      DO 180 I = 1,N
  180 XB(I) = XA(I) + C*D(I)
      CALL FUNC(IPTR,N,XB,FB,GB)
      ISFV = MIN0(2,ISFV)
      IF (FB.GT.F) GO TO 220
      IF (FB.LT.F) GO TO 200
      GL1 = 0.D0
      GL2 = 0.D0
      DO 190 I = 1,N
        GL1 = GL1 + (SCALE(I)*G(I))**2
  190 GL2 = GL2 + (SCALE(I)*GB(I))**2
      IF (GL2.GE.GL1) GO TO 220
  200 ISFV = 3
      F = FB
      DO 210 I = 1,N
        X(I) = XB(I)
  210 G(I) = GB(I)
  220 DGB = 0.D0
      DO 230 I = 1,N
  230 DGB = DGB + GB(I)*D(I)
      IF (FB-FA.LE.0.1D0*C*DGA) GO TO 280
      IF (STEP.GT.STEPLB) GO TO 270
  240 IF (ISFV.GE.2) GO TO 110
  250 IF (IPRINT.EQ.0) GO TO 260
      WRITE (LP,FMT=1070)
 1070 FORMAT (/,5X,'THE RESULTS FROM VA13AD ARE AS FOLLOWS')
      WRITE (LP,FMT=1020) ITR,NFUN
      WRITE (LP,FMT=1030) F
      WRITE (LP,FMT=1040) (X(I),I=1,N)
      WRITE (LP,FMT=1050) (G(I),I=1,N)
  260 RETURN
  270 STEPBD = STEP
      C = GMIN + DGB - 3.D0* (FB-FMIN)/STEP
      CC = DSQRT(C*C-GMIN*DGB)
      C = (C-GMIN+CC)/ (DGB-GMIN+CC+CC)
      STEP = STEP*DMAX1(0.1D0,C)
      GO TO 170
  280 STEPBD = STEPBD - STEP
      STMIN = C
      FMIN = FB
      GMIN = DGB
      STEP = 9.D0*STMIN
      IF (STEPBD.GT.0.D0) STEP = 0.5D0*STEPBD
      C = DGA + 3.D0*DGB - 4.D0* (FB-FA)/STMIN
      IF (C.GT.0D0) STEP = DMIN1(STEP,STMIN*DMAX1(1D0,-DGB/C))
      IF (DGB.LT.0.7D0*DGA) GO TO 170
      ISFV = 4 - ISFV
      IF (STMIN+STEP.LE.STEPLB) GO TO 240
      IR = -N
      DO 290 I = 1,N
        XA(I) = XB(I)
        XB(I) = GA(I)
        D(I) = GB(I) - GA(I)
  290 GA(I) = GB(I)
      CALL MC11AD(H,N,XB,1.D0/DGA,W,IR,1,0.D0)
      IR = -IR
      CALL MC11AD(H,N,D,1D0/ (STMIN* (DGB-DGA)),D,IR,0,0D0)
      IF (IR.LT.N) GO TO 250
      DFF = FA - FB
      FA = FB
      GO TO 130
      END
