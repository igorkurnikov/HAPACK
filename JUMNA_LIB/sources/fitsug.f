      subroutine fitsug(key)
      include 'jumna_data.inc'
      logical*2 lar,lock,kink,lthy
      dimension co(3,15),ct(2,9),dva(5),dvp(5),gs(5)
      character*4 seq*120,code*8,kode*8
      common/mnn/var(n7),gra(n7),con(n0*n2+n3),scl(n7),ncon,nvrc,
     1 ntba,nbac,nthe,nhel,ntki,nkin,ntri,ntot,nvar,nrin,
     1 lar(n7),lock(n6a)
      common/sst/phase(8,3),ampli(8,3),ph5(n2),am5(n2),kr5(n2)
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(120,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      data (co(1,j),j=1,15)/
     1  1.902762514,-0.0064873401,-0.127658642,
     1  0.000067789,-0.0077239115,-0.003000473,
     1  0.001272610,-0.0109947614, 0.054747910,
     1  0.001049814, 0.0035756702, 0.029464986,
     1 -0.001011133, 0.0105723382, 0.067594735/
     1     (co(2,j),j=1,15)/
     1  1.862110934,-0.0157206104,-0.060399117,
     1 -0.000034374, 0.0061854431, 0.009089240,
     1  0.000369412,-0.0036353439,-0.072562664,
     1 -0.001120066,-0.0041634023,-0.036666487,
     1 -0.000688633, 0.0067494541, 0.004192667/
     1     (co(3,j),j=1,15)/
     1  1.841838224,-0.0185462386,-0.076969336,
     1 -0.000489400, 0.0046177816,-0.008918594,
     1  0.000707169,-0.0061303810,-0.047148205,
     1 -0.000244070, 0.0087701119,-0.006696838,
     1  0.001091271,-0.0087170778,-0.017362905/
      data (ct(1,j),j=1,9)/
     1  0.000119607, 0.0003413446, 0.002945780,
     1  0.000407567,-0.8310368695, 0.009535682,
     1 -0.000672256, 0.5800252892,-0.013672097/
     1     (ct(2,j),j=1,9)/
     1  0.000019935,-0.0007499007, 0.002631154,
     1  0.000467024, 0.9764914255, 0.003757503,
     1  0.000101174, 0.0059428552,-0.001042681/
      do ir=1,nto
      k=kr5(ir)
      if(k.eq.0) goto 10
      if(key.eq.0) then
      ph=var(k)
      am=var(k+1)
      ph5(ir)=ph
      am5(ir)=am
         else
         ph=ph5(ir)
         am=am5(ir)
         var(k)=ph
         var(k+1)=am
         endif
      phr=ph*cdr
      amr=am*cdr
      cos1=cos(phr)
      sin1=sin(phr)
      cos2=cos(2*phr)
      sin2=sin(2*phr)
c-------------------------------------------calculate physical variables
      if(key.eq.0) then
      do i=1,3
      var(k+i-5)=(co(i,1)+co(i,2)*amr+co(i,3)*amr**2+co(i,4)*cos1+
     &            co(i,5)*amr*cos1+co(i,6)*amr**2*cos1+co(i,7)*cos2+
     &            co(i,8)*amr*cos2+co(i,9)*amr**2*cos2+co(i,10)*sin1+
     &            co(i,11)*amr*sin1+co(i,12)*amr**2*sin1+co(i,13)*sin2+
     &            co(i,14)*amr*sin2+co(i,15)*amr**2*sin2)*crd
      enddo
      var(k)  =(ct(1,1)+ct(1,2)*amr+ct(1,3)*amr**2+ct(1,4)*cos1+
     &          ct(1,5)*amr*cos1+ct(1,6)*amr**2*cos1+ct(1,7)*sin1+
     &          ct(1,8)*amr*sin1+ct(1,9)*amr**2*sin1)*crd
      var(k+1)=(ct(2,1)+ct(2,2)*amr+ct(2,3)*amr**2+ct(2,4)*cos1+
     &          ct(2,5)*amr*cos1+ct(2,6)*amr**2*cos1+ct(2,7)*sin1+
     &          ct(2,8)*amr*sin1+ct(2,9)*amr**2*sin1)*crd
c----------------------------------------------------calculate gradients
      else
      do i=1,3
      dva(i)= co(i,2)+co(i,3)*2*amr+co(i,5)*cos1+co(i,6)*2*amr*cos1+
     &        co(i,8)*cos2+co(i,9)*2*amr*cos2+co(i,11)*sin1+
     &        co(i,12)*2*amr*sin1+co(i,14)*sin2+co(i,15)*2*amr*sin2
      dvp(i)=-co(i,4)*sin1-co(i,5)*amr*sin1-co(i,6)*amr**2*sin1-
     &        co(i,7)*2*sin2-co(i,8)*amr*2*sin2-co(i,9)*2*amr**2*sin2+
     &        co(i,10)*cos1+co(i,11)*amr*cos1+co(i,12)*amr**2*cos1+
     &        co(i,13)*2*cos2+co(i,14)*amr*2*cos2+
     &        co(i,15)*2*amr**2*cos2
      enddo
      dva(4)=ct(1,2)+ct(1,3)*2*amr+ct(1,5)*cos1+ct(1,6)*amr*2*cos1+
     &       ct(1,8)*sin1+ct(1,9)*amr*2*sin1
      dva(5)=ct(2,2)+ct(2,3)*2*amr+ct(2,5)*cos1+ct(2,6)*amr*2*cos1+
     &       ct(2,8)*sin1+ct(2,9)*amr*2*sin1
      dvp(4)=ct(1,7)*cos1+ct(1,8)*amr*cos1+ct(1,9)*amr**2*cos1-
     &       ct(1,4)*sin1-ct(1,5)*amr*sin1-ct(1,6)*amr**2*sin1
      dvp(5)=ct(2,7)*cos1+ct(2,8)*amr*cos1+ct(2,9)*amr**2*cos1-
     &       ct(2,4)*sin1-ct(2,5)*amr*sin1-ct(2,6)*amr**2*sin1
         do j=1,5
         if(j.le.3) then
         gs(j)=gra(k+j-5)
         gra(k+j-5)=0.
            else
            gs(j)=gra(k+j-4)
            gra(k+j-4)=0.
            endif
         enddo
         do j=1,5
         gra(k)  =gra(k)  +gs(j)*dvp(j)
         gra(k+1)=gra(k+1)+gs(j)*dva(j)
         enddo
      endif
10    enddo
      return
      end
