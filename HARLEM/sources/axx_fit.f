	subroutine fitstr(nat,x1,x2,eps,irot,r,cc,w,ia)
      implicit real*8 (a-h,o-z)
C
C     ----- FITSTR moves two structures to their center of masses 
C     and then ROTATES THE ATOMS WITH COORDINATES X2 ABOUT THE
C           ORIGIN SUCH THAT THE FUNCTION 
C                  F  =  0.5*W*(X-XP)**2 
C         HAS A MINIMUM. REF: A.D. MCLACHLAN, J.MOL.BIOL. 128(1979) 49
C
C INPUT:
C Nat -number of atoms in structures
C IA(2*NR)- working array at least 2*NR dimension
C w(NR+3*NR+3*NR) - working array ( initialized in program 
C  to the array of weights snd temporal coordinate storing)
C X1(3,NR) - first molecule  (preserved)
C X2 (3,NR) - second molecule (preserved)
C
C
C OUTPUT:
C
C eps = 0.5*W*(X-XP)**2  - in matched orientation
C IROT = 1 - in a successful run
C      = 0 - in a failure
C R(3,3) - rotation matrix 
C CC(3) -translational vector
C such that to restore X1:
C X1(i,j) ~= \sum_k R(i,k) X2(k) +cc(i)
C


      logical error
      dimension x1(*),x2(*),r(3,3),w(*),ia(*),cc(*)
      dimension c1(3),c2(3), v1(3),v2(3)
      dimension addat1(6),addat2(6)

      small=1.0d-3
      natadd=0

      if(nat.lt.1)then
         Write(*,*)' error in fitstr: number of atoms less than 1'
         stop
      endif

      if(nat.eq.1)then
          Write(*,*)' warning in fitstr: number of atoms eq 1 ' 
          Write(*,*)' the transformation matrix is set to be unit'
          do i = 1, 3
	       do j = 1, 3
	          r(i,j) = 0.0d0
	       enddo
	    enddo
	    r(1,1) = 1.0d0
		r(1,1) = 1.0d0
		r(1,1) = 1.0d0
		cc(1) = x1(1) - x2(1)
		cc(2) = x1(2) - x2(2)
	    cc(3) = x1(3) - x2(3)
          return
      endif

      call chkpln(nat,x1,ierror,small)

      if(ierror.ne.2.and.ierror.ne.1 ) 
     &      call chkpln(nat,x2,ierror,small)

        if(ierror.eq.1) then
          write(*,*)' Error in fitstr: some atoms are coincide'
          write(*,*)' atom sets :'
          do i=1,nat
            write(*,'("  ", 4I10,3x,3(F8.3,2x))')
     &           i,x1(3*i-2),x1(3*i-1),x1(3*i)
          enddo
             write(*,*)"  "
          do i=1,nat
            write(*,'("  ", 4I10,3x,3(F8.3,2x))')
     &           i,x2(3*i-2),x2(3*i-1),x2(3*i)
          enddo
          
          return
        endif

      
      if(ierror.eq.2)then
C
C add two atoms to x1 and x2 if structure are linear:
C

       write(*,*)'Warning in fitstr: number of atoms in a fragment'
       write(*,*)'Less than 3' 
       write(*,*)'Orientation of the fragment cannot be determined' 

	    do i = 1, 3
	      c1(i) = x1(3+i) - x1(i)
	    enddo
		do i = 1, 3
            c2(i) = 0.0d0
	    enddo
          if( dabs(x1(1)).lt.dabs(x1(2)) )then 
            c2(1)=1.0d0
          else
            c2(2)=1.0d0
          endif
          call vprod(v1,c1,c2)
          ss= vlen(3,v1)
          ss= 1.0d0/ss
	    call dscal(3,ss,v1,1)
c          call ascale(3,ss,v1,v1)
          call aadd(3,x1(1),v1,addat1) 
          call vprod(v2,v1,c1)
          ss= vlen(3,v2)
          ss= 1.0d0/ss
          call aadd(3,x1(1),v2,addat1(4)) 

          call asub(3,x2(4),x2(1),c1)
          call aclear(3,c2)
            if( dabs(x2(1)).lt.dabs(x2(2)) )then 
              c2(1)=1.0d0
            else
              c2(2)=1.0d0
            endif
          call vprod(v1,c1,c2)
          ss= vlen(3,v1)
          ss= 1.0d0/ss
c          call ascale(3,ss,v1,v1)
      	call dscal(3,ss,v1,1)
          call aadd(3,x2(1),v1,addat2) 
          call vprod(v2,v1,c1)
          ss= vlen(3,v2)
          ss= 1.0d0/ss
          call aadd(3,x2(1),v2,addat2(4))
          
          natadd=2
   
      endif

      if(ierror.eq.3)then
C
C Add one atom to the structures if atoms are in plane
C

        call asub(3,x2(4),x2(1),c1)
        do i=3,nat
          isave=i
          call asub(3,x2(3*i-2),x2(1),c2)
          call vprod(v1,c1,c2)
          ss= vlen(3,v1)
          if(ss.gt.small) goto 20 
        enddo
20     continue
        ss=1.0d0/ss
c        call ascale(3,ss,v1,v1)
        call dscal(3,ss,v1,1)
        call aadd(3,x2(1),v1,addat2(1))
 	


        call asub(3,x1(4),x1(1),c1)
          call asub(3,x1(3*isave-2),x1(1),c2)
          call vprod(v1,c1,c2)
          ss= vlen(3,v1)
          ss=1/ss
c        call ascale(3,ss,v1,v1)
         call dscal(3,ss,v1,1)
         call aadd(3,x1(1),v1,addat1(1))

        natadd=1
      endif
       
      nr=nat+natadd
  
      do i=1,nr
        w(i)=1.0d0
      enddo
C
C Copy sets of coordinates to the working array:
C
      do i=1,3*nat
         w(nr+i)=x1(i)
         w(4*nr+i)=x2(i)
      enddo
C
C optionally add atoms to form 3-D structure 
C
      do i=1,natadd
         do j=1,3
           w(nr+ 3*nat+3*(i-1)+j) =addat1(3*(i-1)+j)
           w(4*nr+3*nat+3*(i-1)+j)=addat2(3*(i-1)+j)
         enddo
      enddo
          

      call rmsmov(w(nr+1),w(1),nr)
      do i=1,3
       cc(i)=x1(i)-w(nr+i)
      enddo

      call rmsmov(w(4*nr+1),w(1),nr)
      do i=1,3
       c2(i)=w(4*nr+i)-x2(i)
      enddo
 
      
      call lsqstr(nr,w(1),w(nr+1),w(4*nr+1),eps,ia,irot,r)

      do i=1,3 
         do j=1,3
           cc(i)=cc(i)+r(i,j)*c2(j)
         enddo
      enddo
    
      return
      end
C
C******************************************************************
      SUBROUTINE LSQSTR(NR,W,XP,X,E,IA,IROT,R)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C INPUT:
C
C IA(*)- working array at least 2*NR dimension
C NR -number of atoms in structures
C w(*) - array of weights
C XP(3,NR) - first molecule
C X (3,NR) - second molecule
C
C OUTPUT:
C
C E = 0.5*W*(X-XP)**2  - in matched orientation
C IROT = 1 - in a successful run
C      = 0 - in a failure
C R(3,3) - rotation matrix :
C
C XP(i,j) ~= \sum_k R(i,k) X(k)
C
C
C     ----- LSQSTR ROTATES THE ATOMS WITH COORDINATES X ABOUT THE
C           ORIGIN SUCH THAT THE FUNCTION 
C                  F  =  0.5*W*(X-XP)**2 
C           HAS A MINIMUM. REF: A.D. MCLACHLAN, J.MOL.BIOL. 128(1979) 49
C
      DIMENSION W(*),XP(*),X(*),IA(*)
      DIMENSION XJ(3),XPJ(3),U(3,3),COM(21),OM(6,6),VH(3,3),VK(3,3)
      DIMENSION R(3,3),EIG(6),CC(3)
      DATA SQT2 /1.414213562d0/
C
C INITIATE IA ARRAY
C
      DO 100 I = 1,10
         IA(I) = I*(I-1)/2
  100 CONTINUE
C
C     ----- CALCULATE THE MATRIX U AND ITS DETERMINANT -----
C
      DO 40 M2 = 1,3
         DO 38 M1 = 1,3
            U(M1,M2) = 0.0d0
   38    continue
   40 CONTINUE
      I = 0
      DO 70 J = 1,NR
         DO 50 M = 1,3
            XJ(M) = X(I+M)
            XPJ(M) = XP(I+M)
   50    continue
         DO 60 M2 = 1,3
            DO 58 M1 = 1,3
               U(M1,M2) = U(M1,M2)+W(J)*XJ(M1)*XPJ(M2)
   58       continue
   60    CONTINUE
         I = I+3
   70 CONTINUE
C
      DU =  U(1,1)*U(2,2)*U(3,3) + U(1,3)*U(2,1)*U(3,2)
     +     +U(1,2)*U(2,3)*U(3,1) - U(3,1)*U(2,2)*U(1,3)
     +     -U(3,3)*U(2,1)*U(1,2) - U(3,2)*U(2,3)*U(1,1)
      IF( ABS(DU).LT.1.0d-5) GOTO 1002
      SIGD = DU/ ABS(DU)
C
C     ----- CONSTRUCT OMEGA, DIAGONALIZE IT AND DETERMINE H AND K -----
C
      M = 0
      DO 130 M1 = 1,6
      DO 128 M2 = 1,M1
      M = M+1
      IF(M1.GT.3.AND.M2.LE.3) GOTO 126
      COM(M) = 0.0d0
      GOTO 128
  126 COM(M) = U(M2,M1-3)
  128 CONTINUE
  130 CONTINUE
C
      CALL LIGENB(COM,OM,EIG,IA,6,6)
      IF(DU.LT.0.0d0.AND. ABS(EIG(2)-EIG(3)).LT.1.0d-5) GOTO 1004
C
      DO 140 M2 = 1,3
      DO 138 M1 = 1,3
      VH(M1,M2) = SQT2*OM(M1,6-M2+1)
  138 VK(M1,M2) = SQT2*OM(M1+3,6-M2+1)
  140 CONTINUE
      SIG =  (VH(2,1)*VH(3,2)-VH(3,1)*VH(2,2))*VH(1,3)
     +      +(VH(3,1)*VH(1,2)-VH(1,1)*VH(3,2))*VH(2,3)
     +      +(VH(1,1)*VH(2,2)-VH(2,1)*VH(1,2))*VH(3,3)
      IF(SIG.GT.0.0d0) GOTO 160
      DO 150 M = 1,3
      VH(M,3) = -VH(M,3)
  150 VK(M,3) = -VK(M,3)
  160 CONTINUE
C
C     ----- DETERMINE R  -----
C
      DO 230 M2 = 1, 3
         DO 228 M1 = 1, 3
            R(M1,M2) = VK(M1,1)*VH(M2,1) +
     +                 VK(M1,2)*VH(M2,2) +
     +                 VK(M1,3)*VH(M2,3) * SIGD
  228    continue
  230 CONTINUE
C
      I = 0
      DO 270 J = 1, NR
         XPJ(1) = X(I+1)
         XPJ(2) = X(I+2)
         XPJ(3) = X(I+3)
         DO 260 M1 = 1,3
            XJ(M1) = 0.0d0
            DO 258 M2 = 1,3
               XJ(M1) = XJ(M1)+R(M1,M2)*XPJ(M2)
  258       continue
            X(I+M1) = XJ(M1)
  260    continue
         I = I+3
  270 CONTINUE
C
C     ----- CALCULATE E, WHEN REQUIRED -----
C
      E = 0.0d0
      I = 0
      DO 290 J = 1,NR
      DO 288 M = 1,3
      I = I+1
  288 E = E+W(J)*(X(I)-XP(I))**2
  290 CONTINUE
      E = E/2.0d0
      IROT = 1
      RETURN
 1002 WRITE(6,1102)
      IROT = 0
      RETURN
 1004 WRITE(6,1104)
      IROT = 0
      RETURN
 1102 FORMAT(/2x,' DETERMINANT OF U EQUALS ZERO')
 1104 FORMAT(/2x,' DETERMINANT OF U IS NEGATIVE AND CAPITAL OMEGA HAS',
     +' DEGENERATE EIGENVALUES')
      END
C******************************************************************
      SUBROUTINE LIGENB(A,VEC,EIG,IA,N,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
C
C     ----- A GIVENS HOUSHOLDER MATRIX DIAGONALIZATION
C           ROUTINE SAME AS EIGEN BUT WORKS WITH A
C           LINEAR ARRAY -----
C
      DIMENSION A(*),VEC(NDIM,*),EIG(*),IA(*)
      DIMENSION W(50),GAMA(50),BETA(50),BETASQ(50)
      DIMENSION P(50),Q(50),IPOSV(50),IVPOS(50),IORD(50)
      EQUIVALENCE (P(1),Q(1)),(P(1),BETA(1)),(P(1),IVPOS(1))
      EQUIVALENCE (BETASQ(1),IORD(1)),(IPOSV(1),GAMA(1))
C
      DATA ZERO,PT5,ONE,TWO /0.0E+00,0.5E+00,1.0E+00,2.0E+00/
      DATA RHOSQ /1.0E-24/
C
      IF(N .EQ. 0) GO TO 900
      N1 = N-1
      N2 = N-2
      GAMA(1) = A(1)
      IF(N2) 360,340,100
  100 DO 320 NR = 1,N2
      IK = IA(NR+1)+NR
      B = A(IK)
      S = ZERO
      DO 120 I = NR,N2
      IJ = IA(I+2)+NR
  120 S = S+A(IJ)**2
C
C     ----- PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION -----
C
      A(IK) = ZERO
      IF(S .LE. ZERO) GO TO 300
      S = S+B*B
      SGN = +ONE
      IF(B .GE. ZERO) GO TO 140
      SGN = -ONE
  140 SQRTS = SQRT(S)
      D = SGN/(SQRTS+SQRTS)
      TEMP = SQRT(PT5+B*D)
      W(NR) = TEMP
      A(IK) = TEMP
      D = D/TEMP
      B = -SGN*SQRTS
C
C     ----- -D- IS FACTOR OF PROPORTIONALITY. NOW
C           COMPUTE AND SAVE -W- VECTOR. EXTRA SINGLY
C           SUBSCRIPTED -W- VECTOR FOR SPEED -----
C
      DO 160 I = NR,N2
      IJ = IA(I+2)+NR
      TEMP = D*A(IJ)
      W(I+1) = TEMP
  160 A(IJ) = TEMP
C
C     ----- PREMULTIPLY VECTOR -W- BY MATRIX -A- TO
C           OBTAIN -P- VECTOR. SIMULTANEOUSLY ACCUMULATE
C           DOT PRODUCT -WP- -- SCALR -K- -----
C
      WTAW = ZERO
      DO 240 I = NR,N1
      SUM = ZERO
      II = IA(I+1)
      DO 180 J = NR,I
      IJ = II+J+1
  180 SUM = SUM+A(IJ)*W(J)
      I1 = I+1
      IF(N1 .LT. I1) GO TO 220
      DO 200 J = I1,N1
      IJ = IA(J+1)+I+1
  200 SUM = SUM+A(IJ)*W(J)
  220 P(I) = SUM
  240 WTAW = WTAW+SUM*W(I)
      DO 260 I = NR,N1
  260 Q(I) = P(I)-WTAW*W(I)
C
C     ----- NOW FORM -PAP- MATRIX, REQUIRED PART -----
C
      DO 280 J = NR,N1
      QJ = Q(J)
      WJ = W(J)
      JJ = J+1
      DO 280 I = J,N1
      IJ = IA(I+1)+JJ
  280 A(IJ) = A(IJ)-TWO*(W(I)*QJ+WJ*Q(I))
  300 BETA(NR) = B
      BETASQ(NR) = B*B
      IL = IK+1
  320 GAMA(NR+1) = A(IL)
  340 IJ = IA(N)+N-1
      B = A(IJ)
      BETA(N-1) = B
      BETASQ(N-1) = B*B
      IJ = IJ+1
      GAMA(N) = A(IJ)
  360 BETASQ(N) = ZERO
C
C     ----- ADJOIN AN IDENTYTY MATRIX TO BE POST-
C           MULTIPLIED BY ROTATIONS -----
C
      DO 400 I = 1,N
      DO 380 J = 1,N
  380 VEC(I,J) = ZERO
  400 VEC(I,I) = ONE
      M = N
      SUM = ZERO
      NPAS = 1
      GO TO 600
  420 SUM = SUM+SHIFT
      COSA = ONE
      G = GAMA(1)-SHIFT
      PP = G
      PPBS = PP*PP+BETASQ(1)
      PPBR = SQRT(PPBS)
      DO 540 J = 1,M
      COSAP = COSA
      IF(PPBS .NE. ZERO) GO TO 440
      SINA = ZERO
      SINA2 = ZERO
      COSA = ONE
      GO TO 500
  440 SINA = BETA(J)/PPBR
      SINA2 = BETASQ(J)/PPBS
      COSA = PP/PPBR
C
C     ----- POSTMULTIPLY IDENTITY BY -P- TRANSPOSE -----
C
      NT = J+NPAS
      IF(NT .LT. N) GO TO 460
      NT = N
  460 CONTINUE
      DO 480 I = 1,NT
      TEMP = COSA*VEC(I,J)+SINA*VEC(I,J+1)
      VEC(I,J+1) = -SINA*VEC(I,J)+COSA*VEC(I,J+1)
  480 VEC(I,J) = TEMP
  500 DIA = GAMA(J+1)-SHIFT
      U = SINA2*(G+DIA)
      GAMA(J) = G+U
      G = DIA-U
      PP = DIA*COSA-SINA*COSAP*BETA(J)
      IF(J .NE. M) GO TO 520
      BETA(J) = SINA*PP
      BETASQ(J) = SINA2*PP*PP
      GO TO 560
  520 PPBS = PP*PP+BETASQ(J+1)
      PPBR = SQRT(PPBS)
      BETA(J) = SINA*PPBR
  540 BETASQ(J) = SINA2*PPBS
  560 GAMA(M+1) = G
C
C     ----- TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT -----
C
      NPAS = NPAS+1
      IF(BETASQ(M) .GT. RHOSQ) GO TO 620
  580 EIG(M+1) = GAMA(M+1)+SUM
  600 BETA(M) = ZERO
      BETASQ(M) = ZERO
      M = M-1
      IF(M .EQ. 0) GO TO 660
      IF(BETASQ(M) .LE. RHOSQ) GO TO 580
C
C     ----- TAKE ROOT OF CORMER 2 BY 2 NEAREST TO
C           LOWER DIAGONAL IN VALUE AS ESTIMATE OF
C           EIGENVALUE TO USE FOR SHIFT -----
C
  620 A2 = GAMA(M+1)
      R2 = PT5*A2
      R1 = PT5*GAMA(M)
      R12 = R1+R2
      DIF = R1-R2
      TEMP = SQRT(DIF*DIF+BETASQ(M))
      R1 = R12+TEMP
      R2 = R12-TEMP
      DIF = ABS(A2-R1)- ABS(A2-R2)
      IF(DIF .LT. ZERO) GO TO 640
      SHIFT = R2
      GO TO 420
  640 SHIFT = R1
      GO TO 420
  660 EIG(1) = GAMA(1)+SUM
      DO 680 J = 1,N
      IPOSV(J) = J
      IVPOS(J) = J
  680 IORD(J) = J
      M = N
      GO TO 740
  700 DO 720 J = 1,M
      IF(EIG(J) .LE. EIG(J+1)) GO TO 720
      TEMP = EIG(J)
      EIG(J) = EIG(J+1)
      EIG(J+1) = TEMP
      ITEMP = IORD(J)
      IORD(J) = IORD(J+1)
      IORD(J+1) = ITEMP
  720 CONTINUE
  740 M = M-1
      IF(M .NE. 0) GO TO 700
      IF(N1 .EQ. 0) GO TO 800
      DO 780 L = 1,N1
      NV = IORD(L)
      NP = IPOSV(NV)
      IF(NP .EQ. L) GO TO 780
      LV = IVPOS(L)
      IVPOS(NP) = LV
      IPOSV(LV) = NP
      DO 760 I = 1,N
      TEMP = VEC(I,L)
      VEC(I,L) = VEC(I,NP)
  760 VEC(I,NP) = TEMP
  780 CONTINUE
C
C     ----- BACK TRANSFORM THE VECTORS OF THE TRIPLE
C           DIAGONAL MATRIX -----
C
  800 DO 880 NRR = 1,N
      K = N1
  820 K = K-1
      IF(K .LE. 0) GO TO 880
      SUM = ZERO
      DO 840 I = K,N1
      IJ = IA(I+1)+K
  840 SUM = SUM+VEC(I+1,NRR)*A(IJ)
      SUM = SUM+SUM
      DO 860 I = K,N1
      IJ = IA(I+1)+K
  860 VEC(I+1,NRR) = VEC(I+1,NRR)-SUM*A(IJ)
      GO TO 820
  880 CONTINUE
  900 CONTINUE
      RETURN
      END
C******************************************************************
      SUBROUTINE rmsmov(C,WT,NAT)
      Implicit real*8(a-h,o-z)
C
C     Rev A:  This is a modification of MOVECM from Amber 3.0.  Here we
c     translate the molecule to the center of geometry, based on
c     unit weights.   Weights may be zero for atoms to be ignored.
c     This routine does NOT use inverse atomic masses as MOVECM did.
C
      DIMENSION T(6),VEC(9),EIG(3)
      DIMENSION C(*),WT(*) 
C
      DATA ZERO/0.0D+00/
C
      SUMX = ZERO
      SUMY = ZERO
      SUMZ = ZERO
      SUM  = ZERO
      I3 = 0
      DO 100 I =1,NAT
      WTI = WT(I)
      SUMX = SUMX+C(I3+1)*WTI
      SUMY = SUMY+C(I3+2)*WTI
      SUMZ = SUMZ+C(I3+3)*WTI
      SUM = SUM+WTI
      I3 = I3+3
  100 CONTINUE
      XC = SUMX/SUM
      YC = SUMY/SUM
      ZC = SUMZ/SUM
      I3 = 0
      DO 120 I = 1,NAT
      C(I3+1) = C(I3+1)-XC
      C(I3+2) = C(I3+2)-YC
      C(I3+3) = C(I3+3)-ZC
      I3 = I3+3
  120 CONTINUE
C      WRITE(6,9008) XC,YC,ZC
      RETURN
 9008 FORMAT(5X,'CENTER OF MASS:',3F10.4)
      END
C******************************************************************
      Subroutine chkpln(nat,c,ierror,small)
      Implicit Real*8(a-h,o-z)
C
C Programs check if atoms in c are in the same plane 
C in one line or some of them coincide within a threshold small
C
C ierror=0 - atoms are not coincide and are not in one plane
C       
C
      logical Error
      dimension c(*)
      dimension c1(3),c2(3),v1(3),v2(3),v3(3)
      
      ierror=0
      if(nat.lt.2) then
        ierror=1
        return
      endif

      do i=1,(nat-1)
        do j=i+1,nat
           call asub(3,c(3*i-2),c(3*j-2),c1)
           if(vlen(3,c1).lt.small) then 
              ierror = 1
              return
           endif
         enddo
       enddo

       if(nat.lt.3)then
         ierror=2
         return
       endif
       
       call asub(3,c(4),c(1),c1)
 
       do i=3,nat
         call asub(3,c(3*i-2),c(1),c2)
         call vprod(v1,c1,c2)
         if(vlen(3,v1).gt.small) goto 20
       enddo
         ierror=2
         return
20     continue

       do i=3,nat
         call asub(3,c(3*i-2),c(1),c2)
         call vprod(v2,c1,c2)
         call vprod(v3,v1,v2)
         if(vlen(3,v3).gt.small)goto 30
       enddo
         ierror=3
         return
30     continue
       return
       end
      
