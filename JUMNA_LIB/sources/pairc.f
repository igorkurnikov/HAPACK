      subroutine pairc(k)
      include 'jumna_data.inc'
      logical*2 lthy,kink,hst,bst,vst,lgi,lgj,lgu,ribose,cation
      character*4 mnam,munit,seq*120,code*8,kode*8
      integer*2 i23,i34,elim
      common/enf/for(n1,3),tor(n1,3),fot(n6),ener,elec,
     1 repl,disp,eang,etog,epen,i23(8*n1),i34(8*n1),elim(n1)
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/flx/sap(n6),refg,refb,refh,refv,iap(n8,n4,n5),
     1 nap(n8,7,n5),kap(3,n5)
      common/ind/iend(0:4),ise(n2),ienb(n2,2),kapt(n2+1),nsph
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/slg/hst(n2,6),bst(n2,n8),vst(n2,4),lgi(n1),lgj(n1),lgu(n2)
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
c------------------------------------------------------------set units to ignore
      do is=1,nto
      il=nuc(is-1)+iofs(is)
      ie=nuc(is-1)+iofe(is)
      do i=nuc(is-1)+1,nuc(is)
      elim(i)=0
      if(i.ge.il.and.i.le.ie) elim(i)=is
      enddo
      enddo
c-----------------------------------------------------------------------make i23
      k=0
      kk=0
      do is=1,nto
      do i=nuc(is-1)+1,nuc(is)
      if(lgi(i).or.lgj(i)) then
      k=k+1
      kold=k
      kk=kk+1
      kkold=kk
      do j=1,matd(i,7)
      im=abs(matd(i,j))
      if(im.gt.i) then
      k=k+1
      i23(k)=im
      endif
      do l=1,matd(im,7)
      imm=abs(matd(im,l))
      if(imm.gt.i) then
      k=k+1
      i23(k)=imm
      endif
      do ll=1,matd(imm,7)
      immm=abs(matd(imm,ll))
      if(immm.gt.i) then
      kk=kk+1
      i34(kk)=immm
      endif
      enddo
      enddo
      enddo
      i23(kold)=k-kold
      i34(kkold)=kk-kkold
      endif
      enddo
      enddo
      if(k.gt.8*n1) then
      write(6,*) '  ---- OVERFLOW OF I23 ----'
      stop
      endif
      if(kk.gt.8*n1) then
      write(6,*) '  ---- OVERFLOW OF I34 ----'
      stop
      endif
      return
      end
