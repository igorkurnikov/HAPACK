* *******************************************************************
* COPYRIGHT (c) 1973 Hyprotech UK
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
*######DATE 1 Feb 1993
C       Toolpack tool decs employed.
C       Arg dimensions set to *.
C
      SUBROUTINE MC11AD(A,N,Z,SIG,W,IR,MK,EPS)
      DOUBLE PRECISION EPS,SIG
      INTEGER IR,MK,N
      DOUBLE PRECISION A(*),W(*),Z(*)
      DOUBLE PRECISION AL,B,GM,R,TI,TIM,V,Y
      INTEGER I,IJ,IP,J,MM,NP
      IF (N.GT.1) GO TO 1
      A(1) = A(1) + SIG*Z(1)**2
      IR = 1
      IF (A(1).GT.0.0D0) RETURN
      A(1) = 0.D0
      IR = 0
      RETURN
    1 CONTINUE
      NP = N + 1
      IF (SIG.GT.0.0D0) GO TO 40
      IF (SIG.EQ.0.0D0 .OR. IR.EQ.0) RETURN
      TI = 1.0D0/SIG
      IJ = 1
      IF (MK.EQ.0) GO TO 10
      DO 7 I = 1,N
        IF (A(IJ).NE.0.0D0) TI = TI + W(I)**2/A(IJ)
    7 IJ = IJ + NP - I
      GO TO 20
   10 CONTINUE
      DO 11 I = 1,N
   11 W(I) = Z(I)
      DO 15 I = 1,N
        IP = I + 1
        V = W(I)
        IF (A(IJ).GT.0.0D0) GO TO 12
        W(I) = 0.D0
        IJ = IJ + NP - I
        GO TO 15
   12   CONTINUE
        TI = TI + V**2/A(IJ)
        IF (I.EQ.N) GO TO 14
        DO 13 J = IP,N
          IJ = IJ + 1
   13   W(J) = W(J) - V*A(IJ)
   14   IJ = IJ + 1
   15 CONTINUE
   20 CONTINUE
      IF (IR.LE.0) GO TO 21
      IF (TI.GT.0.0D0) GO TO 22
      IF (MK-1) 40,40,23
   21 TI = 0.D0
      IR = -IR - 1
      GO TO 23
   22 TI = EPS/SIG
      IF (EPS.EQ.0.0D0) IR = IR - 1
   23 CONTINUE
      MM = 1
      TIM = TI
      DO 30 I = 1,N
        J = NP - I
        IJ = IJ - I
        IF (A(IJ).NE.0.0D0) TIM = TI - W(J)**2/A(IJ)
        W(J) = TI
   30 TI = TIM
      GO TO 41
   40 CONTINUE
      MM = 0
      TIM = 1.0D0/SIG
   41 CONTINUE
      IJ = 1
      DO 66 I = 1,N
        IP = I + 1
        V = Z(I)
        IF (A(IJ).GT.0.0D0) GO TO 53
        IF (IR.GT.0 .OR. SIG.LT.0.0D0 .OR. V.EQ.0.0D0) GO TO 52
        IR = 1 - IR
        A(IJ) = V**2/TIM
        IF (I.EQ.N) RETURN
        DO 51 J = IP,N
          IJ = IJ + 1
   51   A(IJ) = Z(J)/V
        RETURN
   52   CONTINUE
        TI = TIM
        IJ = IJ + NP - I
        GO TO 66
   53   CONTINUE
        AL = V/A(IJ)
        IF (MM) 54,54,55
   54   TI = TIM + V*AL
        GO TO 56
   55   TI = W(I)
   56   CONTINUE
        R = TI/TIM
        A(IJ) = A(IJ)*R
        IF (R.EQ.0.0D0) GO TO 70
        IF (I.EQ.N) GO TO 70
        B = AL/TI
        IF (R.GT.4.0D0) GO TO 62
        DO 61 J = IP,N
          IJ = IJ + 1
          Z(J) = Z(J) - V*A(IJ)
   61   A(IJ) = A(IJ) + B*Z(J)
        GO TO 64
   62   GM = TIM/TI
        DO 63 J = IP,N
          IJ = IJ + 1
          Y = A(IJ)
          A(IJ) = B*Z(J) + Y*GM
   63   Z(J) = Z(J) - V*Y
   64   CONTINUE
        TIM = TI
        IJ = IJ + 1
   66 CONTINUE
   70 CONTINUE
      IF (IR.LT.0) IR = -IR
      RETURN
      END
      SUBROUTINE MC11BD(A,N,IR)
      INTEGER IR,N
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION AA,V
      INTEGER I,II,IJ,IK,IP,JK,NI,NP
      IR = N
      IF (N.GT.1) GO TO 100
      IF (A(1).GT.0.0D0) RETURN
      A(1) = 0.D0
      IR = 0
      RETURN
  100 CONTINUE
      NP = N + 1
      II = 1
      DO 104 I = 2,N
        AA = A(II)
        NI = II + NP - I
        IF (AA.GT.0.0D0) GO TO 101
        A(II) = 0.D0
        IR = IR - 1
        II = NI + 1
        GO TO 104
  101   CONTINUE
        IP = II + 1
        II = NI + 1
        JK = II
        DO 103 IJ = IP,NI
          V = A(IJ)/AA
          DO 102 IK = IJ,NI
            A(JK) = A(JK) - A(IK)*V
  102     JK = JK + 1
  103   A(IJ) = V
  104 CONTINUE
      IF (A(II).GT.0.0D0) RETURN
      A(II) = 0.D0
      IR = IR - 1
      RETURN
      END
      SUBROUTINE MC11CD(A,N)
      INTEGER N
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION AA,V
      INTEGER II,IJ,IK,IP,JK,NI,NIP,NP
      IF (N.EQ.1) RETURN
      NP = N + 1
      II = N*NP/2
      DO 202 NIP = 2,N
        JK = II
        NI = II - 1
        II = II - NIP
        AA = A(II)
        IP = II + 1
        IF (AA.GT.0.0D0) GO TO 203
        DO 204 IJ = IP,NI
  204   A(IJ) = 0.D0
        GO TO 202
  203   CONTINUE
        DO 201 IJ = IP,NI
          V = A(IJ)*AA
          DO 200 IK = IJ,NI
            A(JK) = A(JK) + A(IK)*V
  200     JK = JK + 1
  201   A(IJ) = V
  202 CONTINUE
      RETURN
      END
      SUBROUTINE MC11DD(A,N,Z,W)
      INTEGER N
      DOUBLE PRECISION A(*),W(*),Z(*)
      DOUBLE PRECISION Y
      INTEGER I,II,IJ,IP,J,K,N1,NP
      IF (N.GT.1) GO TO 300
      Z(1) = Z(1)*A(1)
      W(1) = Z(1)
      RETURN
  300 CONTINUE
      NP = N + 1
      II = 1
      N1 = N - 1
      DO 303 I = 1,N1
        Y = Z(I)
        IF (A(II).EQ.0.0D0) GO TO 302
        IJ = II
        IP = I + 1
        DO 301 J = IP,N
          IJ = IJ + 1
  301   Y = Y + Z(J)*A(IJ)
  302   Z(I) = Y*A(II)
        W(I) = Z(I)
  303 II = II + NP - I
      Z(N) = Z(N)*A(II)
      W(N) = Z(N)
      DO 311 K = 1,N1
        I = N - K
        II = II - NP + I
        IF (Z(I).EQ.0.0D0) GO TO 311
        IP = I + 1
        IJ = II
        Y = Z(I)
        DO 310 J = IP,N
          IJ = IJ + 1
  310   Z(J) = Z(J) + A(IJ)*Z(I)
  311 CONTINUE
      RETURN
      END
      SUBROUTINE MC11ED(A,N,Z,W,IR)
      INTEGER IR,N
      DOUBLE PRECISION A(*),W(*),Z(*)
      DOUBLE PRECISION V
      INTEGER I,I1,II,IJ,IP,J,NIP,NP
      IF (IR.LT.N) RETURN
      W(1) = Z(1)
      IF (N.GT.1) GO TO 400
      Z(1) = Z(1)/A(1)
      RETURN
  400 CONTINUE
      DO 402 I = 2,N
        IJ = I
        I1 = I - 1
        V = Z(I)
        DO 401 J = 1,I1
          V = V - A(IJ)*Z(J)
  401   IJ = IJ + N - J
        W(I) = V
  402 Z(I) = V
      Z(N) = Z(N)/A(IJ)
      NP = N + 1
      DO 411 NIP = 2,N
        I = NP - NIP
        II = IJ - NIP
        V = Z(I)/A(II)
        IP = I + 1
        IJ = II
        DO 410 J = IP,N
          II = II + 1
  410   V = V - A(II)*Z(J)
  411 Z(I) = V
      RETURN
      END
      SUBROUTINE MC11FD(A,N,IR)
      INTEGER IR,N
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION AA,V
      INTEGER I,I1,II,IJ,IK,IP,J,JK,K,N1,NI,NIP,NP
      IF (IR.LT.N) RETURN
      A(1) = 1.0D0/A(1)
      IF (N.EQ.1) RETURN
      NP = N + 1
      N1 = N - 1
      II = 2
      DO 511 I = 2,N
        A(II) = -A(II)
        IJ = II + 1
        IF (I.EQ.N) GO TO 502
        DO 501 J = I,N1
          IK = II
          JK = IJ
          V = A(IJ)
          DO 500 K = I,J
            JK = JK + NP - K
            V = V + A(IK)*A(JK)
  500     IK = IK + 1
          A(IJ) = -V
  501   IJ = IJ + 1
  502   CONTINUE
        A(IJ) = 1.0D0/A(IJ)
        II = IJ + 1
        AA = A(IJ)
        IJ = I
        IP = I + 1
        NI = N - I
        DO 511 J = 2,I
          V = A(IJ)*AA
          IK = IJ
          K = IJ - IP + J
          I1 = IJ - 1
          NIP = NI + IJ
          DO 510 JK = K,I1
            A(JK) = A(JK) + V*A(IK)
  510     IK = IK + NIP - JK
          A(IJ) = V
  511 IJ = IJ + NP - J
      RETURN
      END
