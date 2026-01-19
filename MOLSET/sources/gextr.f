c
c  Extracts from Gaussian
c
      Real*8 Function ZDONRM(Method,NAtoms,CZ,IAn,C,Gss,Gsd,Gdd)
      Implicit Real*8(A-H,O-Z)
C
C     This function computes the nuclear repulsion energy for NDO/1 and /S.
C
      Dimension CZ(*), IAn(*), C(3,*), V(3), Gss(*),
     $  Gsd(NAtoms,NAtoms), Gdd(*)
      Save Zero
      Data Zero/0.0D0/
C
      Rep = Zero
      If(Method.eq.4) then
        Do 10 I = 2, NAtoms
          II = (I*(I-1))/2
          Call ZNDWtCM(IAn(I),Wt1,Wt2,SI,PI,DI)
          SPI = SI + PI
          Do 10 J = 1, (I-1)
            Call ZNDWtCM(IAn(J),Wt1,Wt2,SJ,PJ,DJ)
            SPJ = SJ + PJ
   10       Rep = Rep + SPI*SPJ*Gss(II+J) + SPI*DJ*Gsd(I,J) +
     $        DI*SPJ*Gsd(J,I) + DI*DJ*Gdd(II+J)
      else
        Do 30 I = 2, NAtoms
          Do 20 J = 1, (I-1)
              Call ASub(3,C(1,I),C(1,J),V)
              R = Sqrt(V(1)*V(1) + V(2)*V(2) + V(3)*V(3))
              Rep = Rep + CZ(I)*CZ(J)/R
   20       Continue
   30     Continue
        endIf
      ZDONRM = Rep
      Return
      End

*Deck ZDOFN
      SUBROUTINE ZDOFNM(NAtoms,CZ,C,FXYZ)
      Implicit Real*8(A-H,O-Z)
C
C     ROUTINE TO COMPUTE THE CONTRIBUTION OF THE V(NN) TERM TO
C     THE FIRST DERIVATIVES.
C
C     ARGUMENTS
C
C     NAtoms ... NUMBER OF ATOMS.
C     CZ    ... ARRAY CONTAINING ATOMIC CHARGES.
C     C      ... ARRAY CONTAINING THE CARTESIAN COORDINATES OF THE
C                NAtoms CENTERS, STORED (X,Y,Z) FOR EACH ATOM.
C     FXYZ   ... OUTPUT ARRAY OF LENGTH 3*NAtoms CONTAINING
C                THE DERIVATIVE CONTRIBUTIONS.
C
      DIMENSION CZ(*), C(3,*), FXYZ(3,*), AB(3)
      Save ZERO
      DATA ZERO/0.D0/
C
      DO 30 I = 2, NAtoms
        I1 = I - 1
        DO 20 J = 1, I1
          R = ZERO
          DO 10 K = 1,3
            AB(K) = C(K,I) - C(K,J)
   10       R = R + AB(K)**2
          FN = CZ(I)*CZ(J)/(R*SQRT(R))
          DO 20 K = 1,3
            R = AB(K)*FN
            FXYZ(K,I) = FXYZ(K,I) - R
   20       FXYZ(K,J) = FXYZ(K,J) + R
   30     CONTINUE
      RETURN
      END
      Subroutine ZNDWtCM(IA,Wt1,Wt2,SEle,PEle,DEle)
      Implicit Real*8(A-H,O-Z)
C
C     Return the weighting and average populations for ZINDO for atomic
C     number IA.
C
      Parameter (MaxAn=48)
      Dimension Wt(MaxAn), LMNMp1(3,MaxAn), LMNMp2(3,MaxAn)
      Save Wt, LMNMp1, LMNMp2, One
      Data LMNMp1/1,0,0,2,0,0,1,0,0,2,0,0,2,1,0,2,2,0,2,3,0,2,4,0,
     $  2,5,0,2,6,0,1,0,0,2,0,0,2,1,0,2,2,0,2,3,0,2,4,0,2,5,0,2,6,0,
     $  1,0,0,2,0,0,2,0,1,2,0,2,2,0,3,2,0,4,2,0,5,2,0,6,2,0,7,2,0,8,
     $  2,0,9,2,0,10,2,1,0,2,2,0,2,3,0,2,4,0,2,5,0,2,6,0,
     $  1,0,0,2,0,0,2,0,1,2,0,2,2,0,3,2,0,4,2,0,5,1,0,7,1,0,8,1,0,9,
     $  1,0,10,2,0,10/
      Data LMNMp2/1,0,0,2,0,0,1,0,0,2,0,0,2,1,0,2,2,0,2,3,0,2,4,0,
     $  2,5,0,2,6,0,1,0,0,2,0,0,2,1,0,2,2,0,2,3,0,2,4,0,2,5,0,2,6,0,
     $  1,0,0,2,0,0,1,0,2,1,0,3,1,0,4,1,0,5,1,0,6,1,0,7,1,0,8,1,0,9,
     $  1,0,10,2,0,10,2,1,0,2,2,0,2,3,0,2,4,0,2,5,0,2,6,0,
     $  1,0,0,2,0,0,2,0,1,2,0,2,2,0,3,2,0,4,2,0,5,1,0,7,1,0,8,1,0,9,
     $  1,0,10,2,0,10/
      Data Wt/20*1.0d0,0.9399d0,0.9069d0,0.8390d0,0.7052d0,0.6652d0,
     $   0.3143d0,0.2065d0,0.1421d0,0.0956d0,1.0d0,18*0.d0/
      Data One/1.0d0/
C
      Wt1 = Wt(IA)
      Wt2 = One - Wt1
      SEle = Wt1*GFloat(LMNMp1(1,IA)) + Wt2*GFloat(LMNMp2(1,IA))
      PEle = Wt1*GFloat(LMNMp1(2,IA)) + Wt2*GFloat(LMNMp2(2,IA))
      DEle = Wt1*GFloat(LMNMp1(3,IA)) + Wt2*GFloat(LMNMp2(3,IA))
      Return
      End
c
c
      Subroutine FFormn(IOut,IPrint,N,NAtoms,IAn,IAtTyp,IU,LLIM,H,P1,P2,
     $  GSS,GSD,GDD,F1,Energy,G1,F2)
      Implicit Real*8(A-H,O-Z)
C
C     Given the mapping array IU, core hamiltonian H, density matrices
C     P1 and P2, and coulomb matrix G form the Fock matrix F1 and
C     compute the energy.
C
      Dimension IU(*), H(*), P1(*), P2(*), GSS(*), GSD(*), GDD(*),
     $	  F1(*), G1(18), F2(18), LLIM(*), IAn(*), IAtTyp(*)
      Save Half, Two, Three, Five, Six, TwFive
      Data Half/0.5D0/, Two/2.0D0/, Three/3.0D0/, Five/5.0D0/,
     $  Six/6.0D0/, TwFive/25.0D0/
      LInd(I,J) = (Max(I,J)*(Max(I,J)-1)/2) + Min(I,J)
C
C     Diagonal contributions.
C
      Do I = 1, N
        IA = IU(I)
	  LIDX = I-LLIM(IA)+1
        IndF = LInd(I,I)
	  IF( LIDX.LT.5) THEN   
          F1(IndF) = H(IndF) - P1(IndF) * GSS(LInd(IA,IA))  ! S and P orbitals
	  ELSE
	    F1(IndF) = H(IndF) - P1(IndF) * GDD(LInd(IA,IA))  ! D orbitals
	  ENDIF

        Do J = 1, N
          IndP = LInd(J,J)
          JA = IU(J)
	    LJDX = J-LLIM(JA)+1
		IF( LJDX.LT. 5 .AND. LIDX.LT.5 ) THEN
			F1(IndF) = F1(IndF) + (P1(IndP)+P2(IndP))*GSS(LInd(JA,IA))
		ELSE IF( LJDX.GT.5 .AND. LIDX.GT.5) THEN
			F1(IndF) = F1(IndF) + (P1(IndP)+P2(IndP))*GDD(LInd(JA,IA))
		ELSE
			F1(IndF) = F1(IndF) + (P1(IndP)+P2(IndP))*GSD(LInd(JA,IA))
		ENDIF	
C
C     Off diagonal contributions.
C
	    IF( I.GT.1 .AND. J.LT.I) THEN
            IndF2 = LInd(I,J)
		  IF( LJDX.LT.5 .AND. LIDX.LT.5 ) THEN
		 	F1(IndF2) = H(IndF2) - P1(IndF2) * GSS(LInd(IA,JA))
		  ELSE IF( LJDX.GT.5 .AND. LIDX.GT.5) THEN
			F1(IndF2) = H(IndF2) - P1(IndF2) * GDD(LInd(IA,JA))
		  ELSE
			F1(IndF2) = H(IndF2) - P1(IndF2) * GSD(LInd(IA,JA))
		  ENDIF	
	    ENDIF
	  ENDDO
      ENDDO 
C
C     INDO contributions.
C
      Do 40 II = 1, NAtoms
        If(IAtTyp(II).lt.0) goto 40
        K = IAn(II)
        If(K.gt.2) then
          I = LLim(II)
          LIndS = LInd(I,I)
          LIndX = LInd(I+1,I+1)
          LIndY = LInd(I+2,I+2)
          LIndZ = LInd(I+3,I+3)
          PAA = P1(LIndS) + P1(LIndX) + P1(LIndY) + P1(LIndZ)
          PAB = P2(LIndS) + P2(LIndX) + P2(LIndY) + P2(LIndZ)
          F1(LIndS) = F1(LIndS) - (PAA-P1(LIndS)) * G1(K) / Three
          Do 30 J = 1, 3
            LIndJ = LInd(I+J,I+J)
            LIndIJ = LInd(I,I+J)
            F1(LIndJ) = F1(LIndJ) +
     $        (P1(LIndJ)-(PAA-P1(LIndS)))*F2(K)/Five -
     $        P1(LIndS)*G1(K)/THREE +
     $        (SIX*P2(LIndJ)-TWO*(PAB-P2(LIndS)))*F2(K)/TwFive
            F1(LIndIJ) = F1(LIndIJ) +
     $        (P1(LIndIJ)+TWO*P2(LIndIJ))*G1(K)/Three
            Do 30 L = 1, 3
              If(J.gt.L) then
                LIndJL = LInd(I+L,I+J)
                F1(LIndJL) = F1(LIndJL) +
     $            (FIVE*P1(LIndJL)+SIX*P2(LIndJL))*F2(K)/TwFive
                endIf
   30         Continue
          endIf
   40   Continue
C
C     Contribution to the energy.
C
      Energy = Energy + Half*(SCFTRCM(P1,H,N,1)+SCFTRCM(P1,F1,N,1))
      Return
      End
c
c
      REAL*8 Function SCFTRCM(A,B,N,NMat)
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
      SCFTrcM = Sum + Sum - Sum1
      Return
      End
	
      SUBROUTINE AUXC(A,B,P)    
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                          
C                                                                       
C     THIS FUNCTION CALCULATES A(P) AND B(-P) FOR THE COULOMB FUNCTION. 
C                                                                                                  
C                                                                       
      DIMENSION A(35),B(35)                                             
C                                                                       
      X=DEXP(P)                                                         
      Y=DEXP(-P)                                                        
      A(1)=Y/P                                                          
      B(1)=(X-Y)/P                                                      
C                                                                       
C  ROUTINE ONLY CALLED FROM GINT1                                      
C  WHICH ONLY USES UPTO N=6                                             
C                                                                       
      DO 1 N=2,6                                                        
      G=DFLOAT(N-1)                                                     
      A(N)=(Y+G*A(N-1))/P                                               
    1 B(N)=-((-1.D0)**(N-1)*Y-X+G*B(N-1))/P                             
      RETURN                                                            
      END  
C	
C	                                                             
      SUBROUTINE AUX(X,T,A,B,FACT)
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                             
C	                                                                       
C     CALCULATES A(X) AND B(X)  for slater orbitals integrals
c                                                                                                                                                  
      DIMENSION A(35),B(35),FACT(35)
c                        
	kmax = 17

      call AIntgm(X,T*X,kmax,35,A)
	call BIntgm(T*X,kmax,35,B,FACT)

c      WRITE(*,* )" P,PT = ", X, T*X
c      WRITE(*,*) " A and B integrals"
c 	do i = 1, kmax
c	WRITE(*,*) i, A(i), B(i)
c	enddo

	return
	end

      Subroutine AIntgm(X,XB,K,MaxAB,A)
      Implicit Real*8(A-H,O-Z)
C
C     A integrals for use in evaluating integrals involving Slater functions.
C
      Dimension A(MaxAB)
      Save Zero
      Data Zero/0.0d0/
C
      EXX = Exp(DAbs(XB)-X)
      A(1) = EXX / X
      Do 10 I = 1, K
   10   A(I+1) = (A(I)*DFloat(I)+EXX) / X
      Do 20 I = (K+1), (MaxAB-1)
   20   A(I+1) = Zero
      Return
      End
      Subroutine BIntgm(X,K,MaxAB,B,fact)
      Implicit Real*8(A-H,O-Z)
C
C     Fills array of B-integrals.  Note that B(I) is B(I-1) in the usual notation.
C     For  X.GT.3                    Exponential formula is used.
C     For  2.LT.X.LE.3 AND K.LE.10   Exponential formula is used.
C     For  2.LT.X.LE.3 AND K.GT.10   15 term series is used.
C     For  1.LT.X.LE.2 AND K.LE.7    Exponential formula is used.
C     For  1.LT.X.LE.2 AND K.GT.7    12 term series is used.
C     For .5.LT.X.LE.1 AND K.LE.5    Exponential formula is used.
C     For .5.LT.X.LE.1 AND K.GT.5    7 term series is used.
C     For  X.LE..5                   6 term series is used.
C
      Dimension B(MaxAB),fact(35)
      Save Pt5, One, Two, Three, Zero
      Real*8 MDCutO
      Data Pt5/0.5d0/, One/1.0d0/, Two/2.0d0/, Three/3.0d0/, Zero/0.0d0/
C
      AbsX = DAbs(X)
      F = Exp(-AbsX)
      If(AbsX.gt.Three.or.(AbsX.gt.Two.and.K.le.10).or.
     $  (AbsX.gt.One.and.K.le.7).or.(AbsX.gt.Pt5.and.K.le.5)) then
        If(X.ge.Zero) then
          ExpX = One
          ExpMX = Exp(-Two*X)
        else
          ExpX = Exp(Two*X)
          ExpMX = One
          endIf
        B(1) = (ExpX-ExpMX)/X
        Do 10 I = 1, K
   10     B(I+1) = (DFloat(I)*B(I)+(-One)**I*ExpX-ExpMX)/X
      else
        If(X.gt.Two) then
          NTerm = 15
        else if(X.gt.One) then
          NTerm = 12
        else if(X.gt.Pt5) then
          NTerm = 7
        else
          NTerm = 6
          endIf
        Do 30 I = 0, K
          B(I+1) = Zero
          XM = One
          Do 20 M = 0, NTerm
            B(I+1) = B(I+1) + XM*(One-(-One)**(M+I+1)) /
     $        (fact(M)*DFloat(M+I+1))
   20       XM = -XM * X
   30     B(I+1) = F*B(I+1)
        endIf
      Do 40 I = (K+2), MaxAB
   40   B(I) = Zero
      Return
      End
c      real*8 Function factm(N)
c      Implicit Real*8(A-H,O-Z)
c      factm = 1.D0
c      Do 10 I = 1, N
c   10   factm = factm * I
c      Return
c      End
C
      SUBROUTINE CBINCF(BINCOE,FACT)
C
C  FORM BINOMIAL COEFFECIENT                                            
C                                                                       
C      ( N )                                                            
C      ( M )  INDEX = N*(N+1)/2 + M+1                                   
C 
C BUNCOE(465) - binominal coefficients returned
C FACT(I+1)   - factorial of I
C
	IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
	DIMENSION BINCOE(*),FACT(*)
      BINCOE(1)=1.D0                                                    
      IND=2                                                             
      DO 60 I=1,35                                                      
      I1=I+1                                                            
      DO 60 J=1,I1                                                      
      JJ=I1-J+1                                                         
      IF(J.LT.26)GO TO 66                                               
      BIN=1.D0                                                          
      IF(J.GE.JJ)THEN                                                   
      DO 62 K=J,I                                                       
62    BIN=BIN*DFLOAT(K)                                                 
      BINCOE(IND)=BIN/FACT(JJ)                                          
      ELSE                                                              
      DO 64 K=JJ,I                                                      
64    BIN=BIN*DFLOAT(K)                                                 
      BINCOE(IND)=BIN/FACT(J)                                           
      END IF                                                            
      GO TO 60                                                          
66    BINCOE(IND)=FACT(I1)/FACT(J)/FACT(JJ)                             
60    IND=IND+1                                                         
	
      RETURN
	END
C      
C 
      REAL*8 FUNCTION SSZ(NN1,NN2,LL1,LL2,M,AMU_X,BMU_X,FACT,BINCOE,LG)     
C  (MOLPAB)             
C                                                          
C     THIS FUNCTION CALCULATES TWO CENTERED OVERLAP INTEGRALS         
C     (NN1,LL1,M) (NN2,LL2,M) - quantum number of Slater orbitals 
C     AMU,BMU - orbital exponents multiplied by R 
C     FACT(I+1) - factorials of I
C
C     A,B - A and B integrals
C     BINCOE - array of binominal coefficients
C     
C     Very ugly need to rewrite later:
C     LG = 1 - call from coulomb integrals function - do not reverse N1-N2, L1 - L2 order
C                                                                                                                               
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)                                                                  
      Parameter (MaxAB=17)      
	DIMENSION A(35),B(35),FACT(30)                  
	DIMENSION BINCOE(*)                      
C                                                                       
      SSZ=0.0D0
	M1 = IAbs(M)
	M2 = IAbs(M)                                                                                                 
      STRAD=0.0D0                                                       

      N1 = NN1
      L1 = LL1
      N2 = NN2
      L2 = LL2
	AMU = AMU_X
	BMU = BMU_X
      P  = (AMU + BMU)/ 2.0D0
      Pt = (AMU - BMU)/ 2.0D0
      X = 0.0D0
C
C     Reverse quantum numbers if necessary.
C
      if(LG.EQ.1)then ! rewrite!!
	else
	  if(L2.lt.L1.or.(L2.eq.L1.and.N2.lt.N1)) then
          K = N1
          N1 = N2
          N2 = K
          K = L1
          L1 = L2
          L2 = K
	    EX = AMU
	    AMU = BMU
	    BMU = EX
          PT = -PT
        endif
	endif
C
C     Find A and B integrals.
C
      K = Mod((N1+N2-L1-L2),2)
C     kmax = N1 + N2
	kmax = 17
      Call AIntgm(P,PT,kmax,MaxAB,A)
C      Call BIntgm(PT,kmax,MaxAB,B,FACT(2)) ! factorial shifted as assumed FACT(N) = (N-1)!
                                                  
      F1=FACT(2*N1+1)                                                   
      F2=FACT(2*N2+1)                                                   
      IFF=L1-M1+1                                                       
      F3=FACT(IFF)                                                      
      IFF=L2-M2+1                                                       
      F4=FACT(IFF)                                                      
      IFF=L1+M1+1                                                       
      F5=FACT(IFF)                                                      
      IFF=L2+M2+1                                                       
      F6=FACT(IFF)                                                      
C      TERM=2.0D0**(N1+N2-L1-L2)*DSQRT(DFLOAT((2*L1+1)*(2*L2+1))*       
C     1   (F3/F1/F5)*(F4/F2/F6)*(AMU/BMU)**(2*N1+1))                    
      Q1=(AMU**N1/DSQRT(F1))                                  
      Q2=(BMU**N2/DSQRT(F2))                                   
      TERM=DSQRT(DFLOAT((2*L1+1)*(2*L2+1))*                             
     1  (F3/F5)*(F4/F6)*AMU*BMU)/2.D0**(L1+L2+1)                      
      TERM=Q1*TERM*Q2                                                   
      JEND=1+((L1-M1)/2)                                                
      KEND=1+((L2-M2)/2)                                                
      DO 50 J=1,JEND                                                    
      JU=J-1                                                            
      IFF=2*L1-2*JU+1                                                   
      F11=FACT(IFF)                                                     
      IFF=L1-M1-2*JU+1                                                  
      F13=FACT(IFF)                                                     
      F15=FACT(JU+1)                                                    
      IFF=L1-JU+1                                                       
      F17=FACT(IFF)                                                     
      DO 50 K=1,KEND                                                    
      KU=K-1                                                            
      IFF=2*L2-2*KU+1                                                   
      F12=FACT(IFF)                                                     
      IFF=L2-M2-2*KU+1                                                  
      F14=FACT(IFF)                                                     
      F16=FACT(KU+1)                                                    
      IFF=L2-KU+1                                                       
      F18=FACT(IFF)                                                     
      CALL CFUNC(N1-L1+2*JU,N2-L2+2*KU,L1-M1-2*JU,L2-M2-2*KU,M1,AMU,  
     1 BMU,VALUE,BINCOE,A,B)                                                      
   50 STRAD=STRAD+VALUE*(F11/F13/F15/F17)*(F12/F14/F16/F18)             
     1  *(-1)**(JU+KU)                                                  
      SSZ=TERM*STRAD                                                   
      RETURN                                                            
      END                                                               
C
      SUBROUTINE CFUNC(IA,IB,IC,ID,IE,AMU,BMU,SNAG,BINCOE,A,B)                  
C                                                                       
C     THIS SUBROUTINE CALCULATES THE C-FUNCTIONS OF ROOTHAAN FOR OVERLAP
C     INTEGRALS. ONLY C FNS WITH POSITIVE INDICES ARE NEEDED AND THEY   
C     WILL BE CALC. BY THE BINOMIAL THEOREM.                            
C              
C     BINCOE - binomial coefficients
C     A,B    - A and B integrals                                       
C                                                                                                                      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	DIMENSION BINCOE(*),A(*),B(*)                                                                                                       
   25 COUNT=0.0D0                                                       
      IAB=IA+1                                                          
      IBB=IB+1                                                          
      ICB=IC+1                                                          
      IDB=ID+1                                                          
      IEB=IE+1                                                          
      INDE= IEB*(IEB-1)/2      ! NIN(IEB)                                            
      INDD= IDB*(IDB-1)/2                                                     
      INDC= ICB*(ICB-1)/2                                                       
      INDB= IBB*(IBB-1)/2                                                     
      INDA= IAB*(IAB-1)/2                                                    
      DO 90 I6=1,IEB                                                    
      B6 = BINCOE(INDE+I6)                                              
      DO 90 I5=1,IEB                                                    
      B5 = BINCOE(INDE+I5)                                              
      DO 90 I4=1,IDB                                                    
      B4=BINCOE(INDD+I4)                                                
      DO 90 I3=1,ICB                                                    
      B3=BINCOE(INDC+I3)                                                
      DO 90 I2=1,IBB                                                    
      B2=BINCOE(INDB+I2)                                                
      DO 90 I1=1,IAB                                                    
      B1=BINCOE(INDA+I1)                                                
      TERM=B1*B2*B3*B4*B5*B6*(-1)**(I2+I5+I6+I4+IE+ID)                  
      IR=I1+I2-I3-I4+IE+IE-I6-I6+IC+ID+3                                
      IP=IA-I1+IB-I2+IE+IE-I5-I5+IC-I3+ID-I4+7                          
   90 COUNT=COUNT+A(IP)*B(IR)*TERM                                      
C      SNAG=COUNT*(BMU/2.0D0)**(IA+IB+IC+ID+IE+IE+1)                   
      SNAG=COUNT                                                        
   92 RETURN                                                            
      END                                                               
	REAL*8 FUNCTION GINT1(EXP1,EXP2,N1,N2,RR,FACT,BINCOE) 
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)
C
C 
C Compute 1 and 2 center coulomb integrals on slater orbitals
C 
C AMU - exponent of the first orbital
C BNU - exponent of the secons orbital
C
C N = N1+(N2*(N2-1))/2 assumed N1 < N2
C
C if (N1 == N2) AMU > BMU 
C
      DIMENSION FACT(30),AFAC(10),A(35),B(35)
	DIMENSION BINCOE(*)
      DATA AFAC/
     X 0.35355339059327D+00,0.13975424859374D+00,0.82679728470768D-01,
     X 0.57210228577644D-01,0.43066295528486D-01,0.34175618247351D-01,
     X 0.28119734619391D-01,0.23755299416522D-01,0.20475351001628D-01,
     X 0.17929397691743D-01/
	if( (N1.LT.N2) .or. ((N1.EQ.N2).and.(AMU.GT.BMU))) then
	  N = N1+(N2*(N2-1))/2
	  AMU = EXP1
	  BMU = EXP2
	else 
	  N = N2+(N1*(N1-1))/2
	  AMU = EXP2
	  BMU = EXP1
	endif 
	R = RR
      IF(DABS(R).LT.1.D-7) R=ZERO 
      TT=(AMU-BMU)/(AMU+BMU)
      T2=TT**2
      IF(T2.LE.1.0D-6)THEN
      TT=ZERO
      T2=ZERO
       AMU=BMU
      END IF
   22 CA=R*AMU
       CB=R*BMU
       C=0.5D0*(AMU+BMU)*R
      IF(AMU-BMU) 11,14,11
   11 D=(AMU**2+BMU**2)/(AMU**2-BMU**2)
   14 E=0.5D0*(AMU+BMU)
      IF(N.GE.10) GO TO 10
      GO TO(1,2,3,4,5,6,7,8,9,10),N
    1 IF(R) 12,13,12
   12 IF(AMU-BMU) 1111,1112,1111
 1111 COUL=(E/C)*(ONE-((ONE-D)**2*((TWO+D)+CA)*DEXP(-CA-CA)+(ONE+D)**2*
     1 ((TWO-D)+CB)*DEXP(-CB-CB))/FOUR)
      GO TO 500
 1112 COUL=(E/C)*(ONE-(ONE+11.0D0*C/EIGHT+0.75D0*C**2+C**3/SIX)*
     1            DEXP(-C- C))
      GO TO 500
   13 COUL=(ONE-T2)*(FIVE-T2)*E/EIGHT
      GO TO 500
    2 IF(R) 2002,2000,2002
 2000 COUL=(ONE-T2)*E*(14.0D0-SEVEN*TT-T2+THREE*T2*TT-T2*T2)/32.0D0
      GO TO 500
 2002 IF(T2.LE.0.0001) GO TO 7
      COUL=(E/C)*(ONE-((ONE-D)**3)*((ONE-FIVE*D-FOUR*D**2)/16.D0
     1      -0.125D0*D*CA)*DEXP(-CA-CA)- (ONE+D)**2*((15.D0
     2      -22.0D0*D+15.0D0*D**2-FOUR*D**3)/16.0D0+THREE*(THREE-
     3      THREE*D+D**2)*CB/EIGHT+0.25D0*(TWO-D)*CB**2+CB**3/12.0D0)
     4      *DEXP (-CB-CB))
      GO TO 500
    3 IF(R) 84,15,84
   15 COUL=(ONE-T2)*(93.0D0-47.0D0*T2+23.0D0*T2*T2-FIVE *T2**3)*E/
     1      256.0D0
      GO TO 500
   84 IF(AMU-BMU)16,17,16
   17 COUL=(E/C)*(ONE-(ONE+419.0D0*C/256.0D0+163.0D0*C**2/128.0D0+
     1     119.0D0*C**3/192.0D0+FIVE*C**4/24.0D0+C**5/20.0D0+C**6/
     2     120.0D0+C**7/1260.0D0)*DEXP (-C- C))
      GO TO 500
   16 COUL=(E/C)*(ONE-(ONE-  D)**3*((EIGHT-D-27.0D0*D**2-30.0D0*D**3
     1      -TEN*D**4)/16.0D0+(11.0D0-19.0D0*D-44.0D0*D**2-20.0D0*
     2      D**3)*CA/32.0D0+(ONE-FIVE*D-FOUR*D**2)*CA**2/16.0D0-D
     3      *CA**3/24.0D0)*DEXP(-CA-CA)-((ONE+D)**3)*((EIGHT+D-27.0D0
     4      *D**2+30.0D0*D**3-TEN*D**4)/16.0D0+(11.0D0+19.0D0*D
     5      -44.0D0*D**2+20.0D0*D**3)*CB/32.0D0+(ONE+FIVE*D-FOUR
     6      *D**2)*CB**2/16.0D0+D*CB**3/24.0D0)*DEXP(-CB-CB))
      GO TO 500
    4 IF(R) 2004,2003,2004
 2003 AA=(BMU/(AMU+BMU))**6
      COUL=BMU*(ONE-AA)/THREE-AMU*AA*BMU/(AMU+BMU)
      GO TO 500
 2004 CALL AUXC(A,B,CA)
      Z=720.0D0*(A(2)*B(1)-A(1)*B(2))
      P=TWO*C
      TT=-TT
C      CALL AUX(P,TT,A,B,FACT(2))
      X=CB
      Y=(X**5)*(A(7)*B(1)+FOUR*A(6)*B(2)+FIVE*A(5)*B(3)-FIVE*A(3)*B(5)
     1  -FOUR*A(2)*B(6)-A(1)*B(7))+TEN*(X**4)*(A(6)*B(1)+THREE*A(5)*
     2  B(2)+TWO*A(4)*B(3)-TWO*A(3)*B(4)-THREE*A(2)*B(5)-A(1)*B(6))+
     3  60.0D0*(X**3)*(A(5)*B(1)+TWO*A(4)*B(2)-TWO*A(2)*B(4)-A(1)*
     4  B(5))+240.0D0*(X**2)*(A(4)*B(1)+A(3)*B(2)-A(2)*B(3)-A(1)*
     5  B(4))+600.0D0*X*(A(3)*B(1)-A(1)*B(3))+720.0D0*(A(2)*B(1)-
     6  A(1)*B(2))
      COUL=( Z-Y)*(R**2)*(AMU**3)/1440.0D0
      GO TO 500
    5 IF(R) 2015,2010,2015
 2010 AA=BMU/(AMU+BMU)
      AA6=AA**6
      ABU=AMU/BMU
      COUL=(ONE-AA6)*BMU/THREE-0.50D0*AA6*AA*AMU*(28.0D0*ABU*ABU*AA*
     1      AA/THREE+SEVEN*ABU*AA+THREE)
      GO TO 500
 2015 CALL AUXC(A,B,CA)
      Z=720.D0*(A(4)*B(1)-THREE*A(3)*B(2)+THREE*A(2)*B(3)-A(1)*B(4))
      P=TWO*C
      TT=-TT
C      CALL AUX(P,TT,A,B,FACT(2))
      X=CB
      Y=-(X**5)*(A(9)*B(1)+TWO*A(8)*B(2)-TWO*A(7)*B(3)-SIX*A(6)*B(4)
     1  +SIX*A(4)*B(6)+TWO*A(3)*B(7)-TWO*A(2)*B(8)-A(1)*B(9))-TEN*
     2  (A(8)*B(1)+A(7)*B(2)-THREE*A(6)*B(3)-THREE*A(5)*B(4)+THREE*
     3  A(4)*B(5)+THREE*A(3)*B(6)-A(2)*B(7)-A(1)*B(8))*X**4-60.0D0*
     4  (A(7)*B(1)-THREE*A(5)*B(3)+THREE*A(3)*B(5)-A(1)*B(7))*X**3
      YY=  -240.0D0*(A(6)*B(1)-A(5)*B(2)-TWO*A(4)*B(3)+TWO*A(3)*B(4)
     5     +A(2)*B(5)-A(1)*B(6))*X**2-600.0D0*(A(5)*B(1)-TWO*A(4)*
     6     B(2)+TWO*A(2)*B(4)-A(1)*B(5))*X-720.0D0*(A(4)*B(1)-THREE*
     7     A(3)*B(2)+THREE*A(2)*B(3)-A(1)*B(4))
      Y=Y+YY
      COUL=(Z+Y)*(R**4)*(AMU**5)/17280.0D0
      GO TO 500
    6 IF(R) 18,19,18
   19 T=BMU+AMU
      COUL=((BMU**7)/(THREE*T**11))*(T**11/(BMU**6)-42.0D0*AMU**5-
     1     42.0D0*T*AMU**4 -28.0D0*AMU**3*(T**2)-14.0D0*(AMU**2)*
     2     (T**3)-FIVE*AMU*(T**4)-T**5)
      GO TO 500
   18 CALL AUXC(A,B,CB)
      Z=A(6)*B(1)-FIVE*A(5)*B(2)+TEN*A(4)*B(3)-TEN*A(3)*B(4)+FIVE*
     1  A(2)*B(5)-A(1)*B(6)
      P=TWO*C
C      CALL AUX(P,TT,A,B,FACT(2))
      X=CA
      U=X**3
      Y=(X**5)*(A(11)*B(1)-FIVE*A(9)*B(3)+TEN*A(7)*B(5)-TEN*A(5)*
     1  B(7)+FIVE*A(3)*B(9)-A(1)*B(11))+TEN*(X**4)*(A(10)*B(1)-A(9)*
     2  B(2)-FOUR*A(8)*B(3)+FOUR*A(7)*B(4)+SIX*A(6)*B(5)-SIX*A(5)*
     3  B(6)-FOUR*A(4)*B(7)+FOUR*A(3)*B(8)+A(2)*B(9)-A(1)*B(10))
      YY=60.0D0*U*(A(9)*B(1)-TWO*A(8)*B(2)-TWO*A(7)*B(3)+SIX*A(6)*
     4   B(4)-SIX*A(4)*B(6)+TWO*A(3)*B(7)+TWO*A(2)*B(8)-A(1)*B(9))+
     5   240.0D0*(X**2)*(A(8)*B(1)-THREE*A(7)*B(2)+A(6)*B(3)+FIVE*
     6   A(5)*B(4)-FIVE*A(4)*B(5)-A(3)*B(6)+THREE*A(2)*B(7)-A(1)*
     7   B(8))+600.0D0*X*(A(7)*B(1)-FOUR*A(6)*B(2)+FIVE*A(5)*B(3)-
     8   FIVE*A(3)*B(5)+FOUR*A(2)*B(6)-A(1)*B(7))+720.0D0*(A(6)*
     9   B(1)-FIVE*A(5)*B(2)+TEN*A(4)*B(3)-TEN*A(3)*B(4)+FIVE*A(2)*
     x   B(5)-A(1)*B(6))
      Y=Y+YY
      COUL=((BMU**7)*(R**6)/(518400.0D0))*(720.0D0*Z-Y)
      GO TO 500
    7 CONTINUE
    8 CONTINUE
    9 CONTINUE
   10 CONTINUE
      NS=N-15
      IF(N.LE.15) GO TO (80,81,81,82,82,82,83,83,83,83,90,90,90,90,90),N
      IF(N.GT.15) GO TO (92,92,92,92,92,92,94,94,94,94,94,94,94),NS
   80 ML=0
      NG=1
      GO TO 85
   81 ML=1
      NG=2
      GO TO 85
   82 ML=3
      NG=3
      GO TO 85
   83 ML=6
      NG=4
      GO TO 85
   90 ML=10
      NG=5
      GO TO 85
   92 ML=15
      NG=6
      GO TO 85
   94 ML=21
      NG=7
   85 NN=2*(N-ML)
      FNN=DFLOAT(NN)
      A2=2.0D0*AMU
      B2=2.0D0*BMU
      MM=2*NG
      IF(R.LT.1.D-5) GO TO 2040
      NGG=2*NG-1
	IONE = 1
      CALL PENET(FACT,BINCOE,BMU,NG,IZERO,BMU,NG,IZERO,RR,PEN,
     . PENPI,PEND,IONE)
C
C  SET PP TO -10.D0 SO AS TO INITIALISE OVLAP
C
      PP=-10.D0
      X= 0.0D0
      DO 100 I=1,NN
      NL=NN-I
      CALL OVLAP(NGG,IZERO,B2,NL,IZERO,A2,RR,S,SP,SD,SF,
     $           FACT,BINCOE,PP)
      FI=DFLOAT(I)/FNN
      AAA=DSQRT(FACT(2*NL+1    ))*FI/(FACT(NN-I+1)*(2.0D0)**(NN-I))
      X=X+S*AAA
  100 CONTINUE
      AAA=DSQRT(BMU**3/AMU)*AFAC(NG)
      COUL=PEN-X*AAA
      GO TO 500
C
C     GENERAL ONE CENTER INTEGRAL OF FORM (NS,NS/MS,MS), WHERE N=NN/2,
C     M=MM/2, AND AMU IS ASSOCIATED WITH N, BMU WITH M.
 2040 A2=TWO*AMU
       B2=TWO*BMU
      P=A2**(NN+1)*B2**(MM+1)/FACT(MM+1)
      COUL=(FACT(MM)/A2**(NN+1))*(1.0D0/(B2**MM)-1.0D0/((A2+B2)**MM))
      FNN=DFLOAT(NN)
      AA=A2*(A2+B2)**(MM+NN)
      DO 2050 I=2,NN
      II=I-1
      T=DFLOAT(II)/(FNN*FACT(NN-II+1))
      AA=AA*A2/(A2+B2)
 2050 COUL=COUL-T*FACT(MM+NN-II)/AA
      COUL=P*COUL
  500 GINT1=COUL
	RETURN
	END
      SUBROUTINE PENET(FACT,BINCOE,AMU,N,LA,BMU,M,LB,RR,PEN,
     $                 PENPI,PEND,I)      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
C     EVALUATES PENETRATION INTEGRALS BETWEEN S SYM. S.T.O.-S. I IS TYPE
C     OF PENET. INT. I=1 FOR INT. X(A)*Y(A)/R(B) WHERE X(A) IS ORB ON CE
C     NTER A, PRINC. QUANTUM NO. N AND EXPONENTIAL CONSTANT AMU. I=2 IS 
C     TYPE  INT.X(A)*X(B)/R(A), WHERE BMU AND M ARE ASSOCIATED WITH X(B)
C     RR IS SEPARATION BETWEEN ATOMIC CENTERS A AND B IN ANGSTROMS, FACT
C     (N+1) IS FACTORIAL N, AND PEN IS THE SIGMA COMPONENT OF THE ANSWER
C     PENPI THE PI COMPONENT, PEND THE DELTA COMPONENT IN HARTREES.     
C     LA, LB = 0,1,2,3                                                  
C              S,P,D,F                                                  
C     M. C. ZERNER   UF.  MODIFIED 1988.                                                                               
C                                                                                                                                                                                           
	DIMENSION A(35),B(35),DD(35),FACT(30)    
C
C DD - undefined - possible error
C              
	DIMENSION BINCOE(*)                                                                    
C                                                                       
      N2=N+M                                                            
      FN=DFLOAT(N2) 
	R = RR                                                                                                          
      GO TO (1,10),I                                                    
    1 CONTINUE                                                          
      PENPI=0.0D0                                                       
      PEND=0.0D0                                                        
C     (X(A)Y(A)/1/R(B))                                                 
      Z=(AMU+BMU)/2.0D0                                                 
      AR=Z*R                                                            
      IF (AR.GT.40.0D0) THEN                                            
         IF(LA.EQ.LB) THEN                                              
          PEN = 1.0D0/R                                                 
          IF(LA.GT.0) PENPI=PEN                                         
          IF(LA.GT.1) PEND =PEN                                         
         ELSE                                                           
          PEN=0.0D0                                                     
         ENDIF                                                          
        GO TO 500                                                       
      ENDIF 
	ONE = 1.0D0                                                            
      IF(R.LT.1.E-6) GO TO 7                                            
C      CALL AUX(AR,ONE,A,B,FACT(2))                                              
      X=0.0                                                             
      IF(LA.GE.LB) THEN                                                 
        IJ=(LA*(LA+1))/2 + LB+1                                         
      ELSE                                                              
        IJ=(LB*(LB+1))/2 + LA+1                                         
      ENDIF                                                             
      GO TO(101,102,103,104,105,106,107,108,109,110),IJ                 
C     S-S                                                               
  101 CONTINUE                                                          
      IBINC=((N2-1)*N2)/2                                               
      DO 5 J=1,N2                                                       
      IBINC=IBINC+1                                                     
      X=X+B(N2-J+1)*A(J)*BINCOE(IBINC)                                  
    5 CONTINUE                                                          
      FACFAC=DSQRT(FACT(N+N+1))*DSQRT(FACT(M+M+1))                      
      PEN=(R**N2)*(AMU**N)*(BMU**M)*DSQRT(AMU*BMU)                      
      PEN=PEN*X/FACFAC                                                  
      GO TO 500                                                         
C     S-P                                                               
  102 CONTINUE                                                          
      IN2=N2-1                                                          
      IBINC=((IN2-1)*IN2)/2                                             
      DO 52 J=1,IN2                                                     
      IBINC=IBINC+1                                                     
      X=X+(B(N2-J)*A(J)+B(N2-J+1)*A(J+1))*BINCOE(IBINC)                 
   52 CONTINUE                                                          
      FACFAC=DSQRT(FACT(N+N+1))*DSQRT(FACT(M+M+1))                      
      PEN=(R**N2)*(AMU**N)*(BMU**M)*DSQRT(AMU*BMU)*DSQRT(3.0D0)         
      PEN=PEN*X/FACFAC                                                  
      GO TO 500                                                         
  103 CONTINUE                                                          
C     P-P  X IS THE SIGMA PART, Y IS THE PI PART.                       
      IN2=N2-2                                                          
      IBINC=((IN2-1)*IN2)/2                                             
      Y=0.0D0                                                           
      DO 53 J=1,IN2                                                     
      IBINC=IBINC+1                                                     
      X=X+(B(N2-J-1)*A(J)+2.0D0*B(N2-J)*A(J+1)+B(N2-J+1)*A(J+2))        
     . *BINCOE(IBINC)                                                   
      Y=Y+(B(J)-B(J+2))*(A(N2-J+1)-A(N2-1-J))*BINCOE(IBINC)             
   53 CONTINUE                                                          
      FACFAC=DSQRT(FACT(N+N+1))*DSQRT(FACT(M+M+1))                      
      PEN=(R**N2)*(AMU**N)*(BMU**M)*DSQRT(AMU*BMU)*3.0D0                
      PENPI=PEN*Y/(2.0D0*FACFAC)                                        
      PEN=PEN*X/FACFAC                                                  
      GO TO 500                                                         
  104 CONTINUE                                                          
C     S-D                                                               
      IN2=N2-2                                                          
      IBINC=((IN2-1)*IN2)/2                                             
      DO 54 J=1,IN2                                                     
      IBINC=IBINC+1                                                     
      X=X+(A(J)*(3.0D0*B(N2-J-1)-B(N2-J+1))+A(J+2)*(3.0D0*B(N2-J+1)-    
     . B(N2-J-1))+4.0D0*A(J+1)*B(N2-J))*BINCOE(IBINC)                   
   54 CONTINUE                                                          
      FACFAC=DSQRT(FACT(N+N+1))*DSQRT(FACT(M+M+1))                      
      PEN=(R**N2)*(AMU**N)*(BMU**M)*DSQRT(AMU*BMU)*DSQRT(3.0D0)         
      PEN=PEN*X/FACFAC                                                  
      GO TO 500                                                         
  105 CONTINUE                                                          
C     P-D  X=SIGMA PART,  Y= PI PART                                    
      IN2=N2-3                                                          
      IBINC=((IN2-1)*IN2)/2                                             
      Y=0.0D0                                                           
      DO 55 J=1,IN2                                                     
      IBINC=IBINC+1                                                     
      X=X+(3.0D0*A(J)*B(IN2-J+1)+7.0D0*(A(J+2)*B(IN2-J+3)+              
     . A(J+1)*B(IN2-J+2))-A(J+2)*B(IN2-J+1)-A(J)*B(IN2-J+3)+            
     . 3.0D0*A(J+3)*B(IN2-J+4)-A(J+3)*B(IN2-J+2)-A(J+1)*B(IN2-J+4))     
     . *BINCOE(IBINC)                                                   
      Y=Y+((A(J+2)-A(J))*(B(IN2-J+1)-B(IN2-J+3))+(A(J+3)-A(J+1))*       
     . (B(IN2-J+2)-B(IN2-J+4)))*BINCOE(IBINC)                           
   55 CONTINUE                                                          
      FACFAC=DSQRT(FACT(N+N+1))*DSQRT(FACT(M+M+1))                      
      PEN=(R**N2)*(AMU**N)*(BMU**M)*DSQRT(AMU*BMU)                      
C     NORMALIZERS  1.936.. = SQRT(15/4), 6.708...= SQRT(45)             
      PENPI=6.708203932*PEN*Y/(2.0D0*FACFAC)                            
      PEN=1.936491673*PEN*X/FACFAC                                      
      GO TO 500                                                         
  106 CONTINUE                                                          
C     D - D   X = SIGMA, Y= PI, ZZ= DELTA PARTS.                        
  107 CONTINUE                                                          
      IN2=N2-4                                                          
      IBINC=((IN2-1)*IN2)/2                                             
      Y=0.0D0                                                           
      ZZ=0.0D0                                                          
      DO 57 J=1,IN2                                                     
      IBINC=IBINC+1                                                     
      X=X+(A(J+4)*(B(IN2-J+1)-6.0D0*B(IN2-J+3)+9.0D0*B(IN2-J+5))+       
     . A(J)*(B(IN2-J+5)+9.0D0*B(IN2-J+1)-6.0D0*B(IN2-J+3))+6.0D0*       
     . A(J+2)*(6.0D0*B(IN2-J+3)-B(IN2-J+5)-B(IN2-J+1))+8.0D0*           
     . A(J+1)*(3.0D0*B(IN2-J+2)-B(IN2-J+4))+8.0D0*A(J+3)*(3.0D0*        
     . B(IN2-J+4)-B(IN2-J+2)))*BINCOE(IBINC)                            
      Y=Y+((A(J+4)-A(J+2))*(B(IN2-J+3)-B(IN2-J+5))+2.0D0*(A(J+3)        
     . -A(J+1))*(B(IN2-J+2)-B(IN2-J+4))+(A(J+2)-A(J))*(B(IN2-J+1)-      
     . B(IN2-J+3)))*BINCOE(IBINC)                                       
      ZZ=ZZ+((A(J+4)-2.0D0*A(J+2)+A(J))*(B(IN2-J+1)-2.0D0*B(IN2-J+3)+   
     . B(IN2-J+5)))*BINCOE(IBINC)                                       
   57 CONTINUE                                                          
      FACFAC=DSQRT(FACT(N+N+1))*DSQRT(FACT(M+M+1))                      
      PEN=(R**N2)*(AMU**N)*(BMU**M)*DSQRT(AMU*BMU)                      
C     NORMALIZERS  3.75=15/4,  1.25 = 5/4                               
      PEND=3.75D0*PEN*ZZ/(2.0D0*FACFAC)                                 
      PENPI=15.0D0*PEN*Y/(2.0D0*FACFAC)                                 
      PEN=1.25D0*PEN*X/FACFAC                                           
      GO TO 500                                                         
  108 CONTINUE                                                          
  109 CONTINUE                                                          
C     YOU HAVE THESE FORMULA'S Z. AS FASCINATING AS THEY ARE            
C     DON'T WORK THEM OUT AGAIN.                                        
  110 CONTINUE                                                          
      WRITE (6,2020)                                                    
 2020 FORMAT (' ****ERROR IN PENET, INTEGRAL NOT AVAILABLE****')        
      CALL EXIT(3)                                                      
    7 IF(LA.EQ.LB) THEN                                                 
       IF (AMU.NE. BMU) GO TO 10                                        
       PEN=AMU*2.0D0/FN                                                 
       IF(LA.GT.0) PENPI=PEN                                            
       IF(LA.GT.1) PEND =PEN                                            
      ELSE                                                              
       PEN=0.0D0                                                        
      ENDIF                                                             
      GO TO 500                                                         
   10 AR=AMU*R                                                          
       BR=BMU*R                                                         
       M2=2*M+1                                                         
      N2=2*N+1                                                          
       P=(AMU+BMU)*R*0.5D0                                              
       T=(AMU-BMU)/(AMU+BMU)                                            
C      CALL AUX(P,T,A,B,FACT(2))                                                 
      T=FACT(M2)*FACT(N2)                                               
      P=DSQRT (AR**N2*BR**M2/T)/R                                       
      M2=M+N                                                            
       T=0.0D0                                                          
      DO 25 J=1,M2                                                      
      N2=M2-J+1                                                         
   25 T=T+DD(J)*B(J)*A(N2)                                              
      PEN=P*T                                                           
  500 CONTINUE                                                          
      RETURN                                                            
      END         
      SUBROUTINE OVLAP(N1,L1,AMU,N2,L2,BMU,R,S,SP,SD,SF,
     $                 FACT,BINCOE,PP) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)       
C                                                                       
C     EVALUATES OVERLAP BETWEEN TWO SLATER TYPE ORBITALS.               
C     N1,N2, ARE PRINC. Q.NO.-S, L1,L2, ARE SECONDARY Q.N.-S,, AMU AND  
C     BMU ARE THE EXPONENTIAL CONSTANTS. R IS THE SEPARATION IN BOHRS
C     S IS THE SIGMA, SP THE PI, SD THE DELTA AND SF THE F COMPOENTS OF 
C     THE OVERLAP, TO BE PUT TOGETHER USING LH AND SUB. GEOM.           
C     FACT(I+1)= FACTORIAL I.                                           
C                                                                                                                                                                                                                                   
      DIMENSION A(35),B(35)
      DIMENSION BINCOE(*)                                                                        
      DIMENSION FACT(30),SS(4)                                                                                           
C                                                                       
      IF(PP.EQ.-10.D0)THEN                                              
      TT= 0.0D0                                                           
      END IF                                                            
      RR=R                                                        
      P=(AMU+BMU)*RR/2.0D0                                               
      T=(AMU-BMU)/(AMU+BMU)                                             
      IF(R.LE.0.001) GO TO 300                                          
      IF(DABS (T).LE.0.0001) T=0.0D0                                     
      IF(P.EQ.PP.AND.T.EQ.TT) GO TO 5                                   
      PP=P                                                              
      TT=T                                                             
C      NM1=N1+N2+1 
C      WRITE(*,*) "(alpha + beta)/2 =",P, "(alpha - beta)/2 =",P*T                                              
C      CALL AUX(P,T,A,B,FACT(2))   	                                                
    5 CONTINUE                                                          
      SS(2)=0.0D0                                                        
      SS(3)=0.0D0                                                      
      SS(4)=0.0D0                                                        
      LMIN=L1                                                           
       IF(L2.LT.L1) LMIN=L2                                             
      LMIN=LMIN+1  
	LG = 1                                                     
      DO 10 II=1,LMIN                                                   
        I=II-1                      
	  AMU_R = AMU*RR
	  BMU_R = BMU*RR                                     
	  SS(II) = SSZ(N1,N2,L1,L2,I,AMU_R,BMU_R,FACT,BINCOE,LG)                                                               
   10 CONTINUE                                                          
      S=SS(1)                                                           
      SP=SS(2)                                                         
      SD=SS(3)                                                         
      SF=SS(4)                                                          
      GO TO 500                                                         
C     ONE CENTER OVERLAP.                                               
  300 S=0.0D0                                                            
      SP=0.0D0                                                          
      SD=0.0D0                                                         
      SF=0.0D0                                                          
      IF(L1.NE.L2) GO TO 500                                            
      NN=2*N1+1                                                         
       MM=2*N2+1                                                        
      F=(1.0D0-T)**MM*(1.0D0+T)**NN                                         
      S=FACT(N1+N2+1)*SQRT(F/(FACT(NN)*FACT(MM)))                      
      IF(L1.GT.0) SP=(-1.0D0)**(L1+1)*S                                   
      IF(L1.GT.1) SD=(-1.0D0)**(L1+2)*S                                   
      IF(L1.GT.2) SF=(-1.0D0)**(L1+3)*S                                   
      S=(-1.0D0)**L1*S                                                    
  500 RETURN                                                            
      END                                                               
      SUBROUTINE FORMP_HA(INIT,NDIM,NBASIS,NE,A,P)
      Implicit Real*8(A-H,O-Z)
C
C     CONVERT THE ORBITALS IN A INTO A DENSITY MATRIX IN P.
C
      LOGICAL INIT
      DIMENSION A(NDIM,1), P(1)
C
      IF(INIT) CALL ACLEAR(((NBASIS*(NBASIS+1))/2),P)
      If(NE.lt.1.or.NBasis.lt.1) Return
      IND = 0
      DO 10 I = 1, NBASIS
        DO 10 J = 1, I
          IND = IND + 1
          DO 10 K = 1, NE
   10       P(IND) = P(IND) + A(I,K)*A(J,K)
      Return
      End