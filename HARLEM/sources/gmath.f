c
c  Some fortran functions for GAUSSIAN interface 
c  
c
      INTEGER FUNCTION OPENFORF(NF,FILENM,FFORM,FSTAT)
      implicit none
      INTEGER NF
      CHARACTER*(*) FILENM, FSTAT, FFORM
      
      OPEN(NF, ERR=110, FILE=FILENM, FORM = FFORM, STATUS= FSTAT)
      
      OPENFORF = 1
      return
110      CONTINUE
      OPENFORF = 0
      return
      end
      
      INTEGER FUNCTION CLSFORF(NF)
      implicit none
      INTEGER NF
      
      CLOSE(NF, ERR=110)

      CLSFORF = 1
      return
110      continue
      CLSFORF = 0
      return
      END

      INTEGER FUNCTION WRTFSTR(NF,STR)
      implicit none
      INTEGER NF
      CHARACTER*(*) STR
      
      WRITE(NF,'(A)',ERR=110)STR
      
      WRTFSTR = 1
      return
110      continue
      wrtfstr = 0
      return
      end

      INTEGER FUNCTION WRTFARRD(NF,DARR,N,FORMSTR)
      implicit none
      INTEGER NF
      INTEGER N
      REAL*8 DARR(N)
      CHARACTER*(*) FORMSTR
      
      WRITE(NF,FORMSTR,ERR=110)DARR
      WRTFARRD = 1
      return
110   wrtfarrd = 0
      return
      end

      INTEGER FUNCTION RDFARRD(NF,DARR,N,FORM)
      implicit none
      INTEGER N
      INTEGER NF
      REAL*8 DARR(N)
      CHARACTER*(*) FORM
      
      READ(NF,FORM,ERR=110)DARR
      RDFARRD = 1
      return
110   RDFARRD = 0
      return
      end


      SUBROUTINE ASCALE(N,S,A,B)
      Implicit Real*8(A-H,O-Z)
      
      call dcopy(N,A,1,B,1)
      call dscal(N,S,B,1)
      

      return
      end
      SUBROUTINE AMOVE_GS(N,A,B)
      Implicit Real*8(A-H,O-Z)
      
      call dcopy(N,A,1,B,1)

      return
      end
      
      Subroutine AAdd(N,A,B,C)
      Implicit Real*8(A-H,O-Z)
C
C     ROUTINE TO DO VECTOR OPERATION
C
C     C = A + B
C
      DIMENSION A(1),B(1),C(1)
      IF(N.LT.1) RETURN
      DO 10 I=1,N
   10     C(I)=A(I)+B(I)
      RETURN
      END

      SUBROUTINE ACLEAR(N,A)
      Implicit Real*8(A-H,O-Z)
C
C     ROUTINE TO CLEAR N ELEMENTS IN ARRAY A.
C
      DIMENSION A(1)
      Save ZERO
      DATA ZERO/0.0D0/
      IF(N.LT.1) RETURN
      DO 10 I=1,N
   10     A(I)=ZERO
      RETURN
      END
      Subroutine APISq(N,S,A)
      Implicit Real*8(A-H,O-Z)
C
C     Add Fact times the unit matrix to square matrix A.
C
      Dimension A(N,N)
C
      Do 10 I = 1, N
   10   A(I,I) = A(I,I) + S
      Return
      End
      SUBROUTINE ASUB(N,A,B,C)
      Implicit Real*8(A-H,O-Z)
C
C     ROUTINE TO PERFORM THE VECTOR OPERATION
C
C     C = A - B
C
C     WHERE C, A, AND B ARE VECTORS OF LENGTH N.
C
      DIMENSION C(1),A(1),B(1)
      IF(N.LT.1) RETURN
      DO 1 I=1,N
    1   C(I)=A(I)-B(I)
      RETURN
      END

      FUNCTION VLEN(N,V)
      Implicit Real*8(A-H,O-Z)
C
C     A FUNCTION WHICH RETURNS THE LENGTH OF VECTOR V.
C
      DIMENSION V(1)

      VLEN = DNRM2(N,V,1)
      RETURN 
      END

      SUBROUTINE VPROD(VP,X,Y)
      Implicit Real*8(A-H,O-Z)
C
C     VP=X CROSS Y
C
      DIMENSION VP(3),X(3),Y(3)
C
      VP(1)=X(2)*Y(3)-X(3)*Y(2)
      VP(2)=X(3)*Y(1)-X(1)*Y(3)
      VP(3)=X(1)*Y(2)-X(2)*Y(1)
      RETURN
      END
      function gfloat (iarg)
      real*8 gfloat
c
c     integer to working-precision conversion
c
      gfloat=dfloat(iarg)
      return
      end
      SUBROUTINE ACASB(N,A,B,C,S)
      Implicit Real*8(A-H,O-Z)
C
C     ROUTINE TO PERFORM VECTOR OPEATION
C
C     C = A + S * B
C
C     WHERE A, B AND C ARE VECTORS OF LENGTH N AND S IS A SCALAR.
C
      DIMENSION A(1),B(1),C(1)
      IF(N.LT.1) RETURN
      DO 10 I=1,N
   10     C(I)=A(I)+S*B(I)
      RETURN
      END
      Function SCFTRC(A,B,N,NMat)
      Implicit Real*8(A-H,O-Z)
      Dimension A(1), B(1)
      Save Zero
      Data Zero/0.0D0/
C
      Sum = Zero
      Sum1 = Zero
      Len = (NMat*N*(N+1))/2
      Do 10 I = 1, Len
   10   Sum = Sum + A(I)*B(I)
      IJ = 0
      Do 20 IM = 1, NMat
        Do 20 I = 1, N
          IJ = IJ + I
   20     Sum1 = Sum1 + A(IJ)*B(IJ)
      SCFTrc = Sum + Sum - Sum1
      Return
      End
C
C Fortran subroutines for NDO calculations
C
      SUBROUTINE TRMATD(T,MAXL,E)
      Implicit Real*8(A-H,O-Z)
      DIMENSION T(9,9), E(3)
      Save ZERO, ONE, TWO, THREE, Cut1, Cut2
      DATA ZERO/0.0D0/, ONE/1.0D0/, TWO/2.0D0/, THREE/3.0D0/,
     $     CUT1/1.0D-10/, CUT2/1.0D-06/
C
      COST = E(3)
      SINT = ZERO
      IF((ONE-COST**2).GT.CUT1) SINT = Sqrt(ONE-COST**2)
      COSP = ONE
      SINP = ZERO
      IF(SINT.GT.CUT2) COSP = E(1) / SINT
      IF(SINT.GT.CUT2) SINP = E(2) / SINT
      CALL ACLEAR(81,T)
C
C     TRANSFormatION MATRIX ELEMENTS FOR D FUNCTIONS.
C
      IF(MAXL.LT.2) GOTO 100
      COS2T = COST**2 - SINT**2
      SIN2T = TWO*SINT*COST
      COS2P = COSP**2 - SINP**2
      SIN2P = TWO*SINP*COSP
      SQRT3 = Sqrt(THREE)
      T(5,5) = (THREE*COST**2-ONE)/TWO
      T(5,6) = -SQRT3*SIN2T/TWO
      T(5,8) = SQRT3*SINT**2/TWO
      T(6,5) = SQRT3*SIN2T*COSP/TWO
      T(6,6) = COS2T*COSP
      T(6,7) = -COST*SINP
      T(6,8) =-T(6,5)/SQRT3
      T(6,9) = SINT*SINP
      T(7,5) = SQRT3*SIN2T*SINP/TWO
      T(7,6) = COS2T*SINP
      T(7,7) = COST*COSP
      T(7,8) = -T(7,5)/SQRT3
      T(7,9) = -SINT*COSP
      T(8,5) = SQRT3   *SINT**2*COS2P/TWO
      T(8,6) = SIN2T*COS2P/TWO
      T(8,7) = -SINT*SIN2P
      T(8,8) = (ONE+COST**2)*COS2P/TWO
      T(8,9) = -COST*SIN2P
      T(9,5) = SQRT3   *SINT**2*SIN2P/TWO
      T(9,6) = SIN2T*SIN2P/TWO
      T(9,7) = SINT*COS2P
      T(9,8) = (ONE+COST**2)*SIN2P/TWO
      T(9,9) = COST*COS2P
C
C     TRANSFormatION MATRIX ELEMENTS FOR P FUNCTIONS.
C
  100 IF(MAXL.LT.1) GOTO 110
      T(2,2) = COST*COSP
      T(2,3) = -SINP
      T(2,4) = SINT*COSP
      T(3,2) = COST*SINP
      T(3,3) = COSP
      T(3,4) = SINT*SINP
      T(4,2) = -SINT
      T(4,4) = COST
C
C     TRANSFormatION MATRIX ELEMENTS FOR S FUNCTIONS.
C
  110 T(1,1) = ONE
      RETURN
      END
      Subroutine ExtrDM(IOut,IPrint,ICount,N,NMat,NTT,IFlag,IExtp,RMSDP,
     $                  PCur,Scr,DeltaP)
      Implicit Real*8(A-H,O-Z)
C
C     3 and 4 point extrapolation routine.  Arguments:
C     IOut   ... Printed output unit.
C     IPrint ... >0 for print.
C     ICount ... Flag.  Should be initialized to 0 by calling routine.
C     N      ... Number of basis functions.
C     NMat   ... Number of matrices (1 for RHF, 2 for UHF).
C     NTT    ... Used for dimensioning.
C     IFlag  ... Set 0 if no extrapolation this time, 1 for 3 point,
C                2 for 4 point, 3 if 3 point skipped, 4 if four-point
C                skipped.
C     IExtp  ... Flag to suppress extrapolation.
C                0 ... Do both 3 and 4 point.
C                1 ... No 3 point.
C                2 ... No 3 or 4 point.
C     PCur   ... Current NMat density matrix on input, extrapolated
C                density on output.
C     Scr    ... NMat*NTT scratch vector.
C     DeltaP ... 3*NMat*NTT array for saved vectors.
C
      Dimension PCur(NTT,NMat), Scr(NTT,NMat), DeltaP(NTT,NMat,3)
      Save Zero, Pt99, Pt995, One, OnePt9, Two, Four, Small
      Save SP12, SP22
      Data Zero/0.0d0/, Pt99/0.99d0/, Pt995/0.995d0/, One/1.0d0/,
     $  OnePt9/1.9d0/, Two/2.0d0/, Four/4.0d0/, Small/1.d-30/,
     $  SP12/0.0d0/, SP22/0.0d0/
C
C     Initialize.  If this is the first cycle or we just extrapolated,
C     there's nothing to do.
C
      NTTM = NMat * NTT
      IFlag = 0
      ICount = ICount + 1
      Loc1 = 3
      If(Mod(ICount,2).eq.1) Loc1 = 2
      Loc2 = 5 - Loc1
      RNMat = GFloat(NMat)
C
C     Form P(N)-P(N-1) in Scr;  find its length DP1.
C
      Call ASub(NTTM,PCur,DeltaP(1,1,1),Scr)
      SP11 = SCFTrc(Scr,Scr,N,NMat)
      DP1 = Sqrt(SP11/RNMat)
      RMSDP = DP1 / GFloat(N)
      If(ICount.eq.1) goto 260
      If(ICount.eq.2) Goto 240
      If(RMSDP.lt.1.d-12) goto 240
      If(ICount.lt.4) goto 140
          SP23 = SP12
          SP33 = SP22
          SP13 = SCFTrc(DeltaP(1,1,Loc1),Scr,N,NMat)
          DP3 = Sqrt(SP33/RNMat)
  140 SP12 = SCFTrc(DeltaP(1,1,Loc2),Scr,N,NMat)
      SP22 = SCFTrc(DeltaP(1,1,Loc2),DeltaP(1,1,Loc2),N,NMat)
      DP2 = Sqrt(SP22/RNMat)
      Denom = RNMat*DP1*DP2
      If(Denom.le.Small) goto 240
      CosPhi = SP12 / Denom
      If(ICount.eq.3) Goto 240
C
C     Find cosine of angle between DP3 and plane of DP1 and DP2.
C
      Denom = (SP11*SP22-SP12*SP12)
      If(Denom.le.Small) goto 240
      X = (SP13*SP22-SP12*SP23) / Denom
      Y = (SP23*SP11-SP12*SP13) / Denom
      If(DP3.le.Small.or.Abs(X).le.Small) goto 240
      CosPsi = Sqrt((X*X*SP11+Y*Y*SP22+Two*X*Y*SP12)/RNMat)/DP3
C
C     Do not extrapolate unless 4 consecutive points are nearly
C     coplanar.
C
      If(IExtp.ge.2) Goto 240
      If(CosPsi.le.Pt99) goto 240
C
C     EXPRESS VECTOR DP1 AS X*DP3(PROJECTED)+Y*DP2
C
      Y = -Y/X
      X = One/X
C     Test if 2*2 matrix has real eigenvalues between -.95 and +.95.
      XY = Y*Y+Four*X
      If(XY.lt.Zero) goto 190
      XY = Abs(Y) + Sqrt(XY)
      If(XY.le.OnePt9) goto 220
C
C     3-pt extrapolation if allowed and no 4-pt.
  190 If(IExtp.ge.1) goto 240
      If(Abs(CosPhi).le.Pt995) goto 240
      X = DP1/(DP2*CosPhi-DP1)
      Call ACasB(NTTM,PCur,Scr,PCur,X)
      IFlag = 1
      ICount = 0
      Goto 240
C
C     4-pt extrapolation if allowed.
C
  220 If(IExtp.ge.2) goto 240
      XXX = X/(One-X-Y)
      YYY = (X+Y)/(One-X-Y)
      Call ACasB(NTTM,PCur,DeltaP(1,1,Loc2),PCur,XXX)
      Call ACasB(NTTM,PCur,Scr,PCur,YYY)
      IFlag = 2
      ICount = 0
C
C     Store final matrices.
C
C     mikola changed Scr to Scr(1,1) to make iFort12 happy
  240 Call AMOVE_GS(NTTM,Scr(1,1),DeltaP(1,1,Loc1))
C     mikola changed PCur to PCur(1,1) to make iFort12 happy
  260 Call AMOVE_GS(NTTM,PCur(1,1),DeltaP(1,1,1))
      Return
      End

      SUBROUTINE TORS_GS(NOINT,I,J,K,L,B,IB,C)
C
C        ADAPTED FROM THE NORMAL COORDINATE ANALYSIS PROGRAM OF
C        SCHACHTSCHNEIDER, SHELL DEVELOPMENT .
C
      Implicit Real*8(A-H,O-Z)
      DIMENSION B(3,4,*),IB(4,*),C(1)
      DIMENSION RIJ(3),RJK(3),RKL(3),EIJ(3),EJK(3),EKL(3),CR1(3),CR2(3)
      Save Zero,One
      DATA ZERO,ONE/0.D0,1.D0/
C
      IAIND=3*(I-1)
      JAIND=3*(J-1)
      KAIND=3*(K-1)
      LAIND=3*(L-1)
      IB(1,NOINT)=I
      IB(2,NOINT)=J
      IB(3,NOINT)=K
      IB(4,NOINT)=L
      DIJSQ=ZERO
      DJKSQ=ZERO
      DKLSQ=ZERO
      DO 124 M=1,3
      RIJ(M)=C(M+JAIND)-C(M+IAIND)
      DIJSQ=DIJSQ+RIJ(M)**2
      RJK(M)=C(M+KAIND)-C(M+JAIND)
      DJKSQ=DJKSQ+RJK(M)**2
      RKL(M)=C(M+LAIND)-C(M+KAIND)
  124 DKLSQ=DKLSQ+RKL(M)**2
      DIJ = SQRT(DIJSQ)
      DJK = SQRT(DJKSQ)
      DKL = SQRT(DKLSQ)
      DO 136 M=1,3
      EIJ(M)=RIJ(M)/DIJ
      EJK(M)=RJK(M)/DJK
  136 EKL(M)=RKL(M)/DKL
      CR1(1)=EIJ(2)*EJK(3)-EIJ(3)*EJK(2)
      CR1(2)=EIJ(3)*EJK(1)-EIJ(1)*EJK(3)
      CR1(3)=EIJ(1)*EJK(2)-EIJ(2)*EJK(1)
      CR2(1)=EJK(2)*EKL(3)-EJK(3)*EKL(2)
      CR2(2)=EJK(3)*EKL(1)-EJK(1)*EKL(3)
      CR2(3)=EJK(1)*EKL(2)-EJK(2)*EKL(1)
      DOTPJ = -(EIJ(1)*EJK(1)+EIJ(2)*EJK(2)+EIJ(3)*EJK(3))
      DOTPK = -(EJK(1)*EKL(1)+EJK(2)*EKL(2)+EJK(3)*EKL(3))
      SINPJ = SQRT(ONE-DOTPJ**2)
      SINPK = SQRT(ONE-DOTPK**2)
      DO 164 M=1,3
      SMI=-CR1(M)/(DIJ*SINPJ*SINPJ)
      B(M,1,NOINT)=SMI
      F1=(CR1(M)*(DJK-DIJ*DOTPJ))/(DJK*DIJ*SINPJ*SINPJ)
      F2=(DOTPK*CR2(M))/(DJK*SINPK*SINPK)
      SMJ=F1-F2
      B(M,2,NOINT)=SMJ
      SML= CR2(M)/(DKL*SINPK*SINPK)
      B(M,4,NOINT)=SML
      B(M,3,NOINT)=(-SMI-SMJ-SML)
  164 CONTINUE
C
      RETURN
      END
      SUBROUTINE BEND_GS(NOINT,I,J,K,B,IB,C)
C
C        ADAPTED FROM THE NORMAL COORDINATE ANALYSIS PROGRAM OF
C        SCHACHTSCHNEIDER, SHELL DEVELOPMENT .
C
      Implicit Real*8(A-H,O-Z)
      DIMENSION B(3,4,*),IB(4,*),C(1)
      DIMENSION RJI(3),RJK(3),EJI(3),EJK(3)
      Save ZERO, ONE
      DATA ZERO,ONE/0.D0,1.D0/
C
      IAIND=3*(I-1)
      JAIND=3*(J-1)
      KAIND=3*(K-1)
      IB(1,NOINT)=I
      IB(2,NOINT)=J
      IB(3,NOINT)=K
      IB(4,NOINT)=0
      DJISQ=ZERO
      DJKSQ=ZERO
      DO 120 M=1,3
      RJI(M)=C(M+IAIND)-C(M+JAIND)
      RJK(M)=C(M+KAIND)-C(M+JAIND)
      DJISQ = DJISQ + RJI(M)**2
  120 DJKSQ = DJKSQ + RJK(M)**2
      DJI = SQRT(DJISQ)
      DJK = SQRT(DJKSQ)
      DOTJ = ZERO
      DO 132 M=1,3
      EJI(M)=RJI(M)/DJI
      EJK(M)=RJK(M)/DJK
  132 DOTJ=DOTJ+EJI(M)*EJK(M)
      SINJ = SQRT(ONE-DOTJ**2)
      DO 144 M=1,3
      B(M,3,NOINT)=(   (DOTJ*EJK(M)-EJI(M)))/(DJK*SINJ)
      B(M,1,NOINT)=(   (DOTJ*EJI(M)-EJK(M)))/(DJI*SINJ)
      B(M,2,NOINT)=-B(M,1,NOINT)-B(M,3,NOINT)
  144 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE STR_GS(NOINT,I,J,B,IB,C)
C
C        ADAPTED FROM THE NORMAL COORDINATE ANALYSIS PROGRAM OF
C        SCHACHTSCHNEIDER, SHELL DEVELOPMENT .
C
      Implicit Real*8(A-H,O-Z)
      DIMENSION B(3,4,*),IB(4,*),C(1),RIJ(3)
      Save ZERO
      DATA ZERO/0.D0/
C
      IAIND=3*(I-1)
      JAIND=3*(J-1)
      IB(1,NOINT)=I
      IB(2,NOINT)=J
      IB(3,NOINT)=0
      IB(4,NOINT)=0
      DIJSQ = ZERO
      DO 112 M=1,3
      RIJ(M)=C(M+JAIND)-C(M+IAIND)
  112 DIJSQ=DIJSQ+RIJ(M)**2
      DO 120 M=1,3
      B(M,1,NOINT) = -RIJ(M)/SQRT(DIJSQ)
      B(M,2,NOINT)=-B(M,1,NOINT)
  120 CONTINUE
C
      RETURN
      END
C
      Subroutine FORMG_GS(NT,IB,B,G)
      Implicit Real*8(A-H,O-Z)
C
C     Given the B-matrix, form the Wilson G-matrix.
C
      Dimension IB(4,NT), B(3,4,NT), G(NT,NT)
C     
      Call AClear(NT*NT,G)  
      Do 130 I = 1, NT
        Do 130 J = 1, I
          Do 120 I1 = 1, 4
            IBI = IB(I1,I)
            If(IBI.ne.0) then
              Do 110 J1 = 1, 4
                If(IBI.eq.IB(J1,J))
     $            G(J,I) = G(J,I) + B(1,I1,I)*B(1,J1,J)
     $              + B(2,I1,I)*B(2,J1,J) + B(3,I1,I)*B(3,J1,J)
  110           Continue
              endIf
  120       Continue
  130     G(I,J) = G(J,I)
      Return
      End
C
      Subroutine TRANF_GS(NParm,FX,F,IB,B,G)
      Implicit Real*8(A-H,O-Z)
C
C     Routine to transform first derivatives to internal coordinates.
C     Arguments:
C
C     NParm  ... Number of Z-Matrix degrees of freedom (3*NZ-6).
C     
C     FX     ... Input vector (length 3*NAtoms) containing Cartesian
C                derivatives.
C     F      ... Output vector of length NParm containing derivatives
C                over internal coordinates.
C     IB     ... Integer B-matrix as produced by FormBG.
C     B      ... B-matrix as produced by FormBG.
C     G      ... G-matrix produced by FormBG.  Only used if UseG is set.
C
      Dimension FX(3,*), F(*), IB(4,NParm,*), B(3,4,NParm,*),
     $  G(NParm,NParm)
      Save Zero
      Data Zero/0.d0/
C
      Ndim = 1
      Call AClear(NParm,F)
      Do 40 I = 1, NParm
        R = Zero
        Do 25 iDim = 1, NDim
          Do 20 K1 = 1, 4
            K = IB(K1,I,iDim)
            If(K.ne.0) R = R + B(1,K1,I,iDim)*FX(1,K)
     $        + B(2,K1,I,iDim)*FX(2,K) + B(3,K1,I,iDim)*FX(3,K)
   20       Continue
   25     Continue
          Do 30 J = 1, NParm
   30       F(J) = F(J) + G(J,I)*R
   40   Continue
      Return
      End