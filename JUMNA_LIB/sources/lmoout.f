      subroutine lmoout
      include 'jumna_data.inc'
      logical*2 lar,lock,lthy,kink,ribose,cation,sup,rcom,homo,homo2,
     1 homo3,diep,sum,link,ecen,cyl,lcat,cent,autos,locr,tlp,quiet,
     1 convg,dloop,amber
      character*4 mnam,munit,seq*120,code*8,kode*8,mac*32,lmo*32,
     1 lib*32,libm*32,out*32,axe*32,noe*32,nol*32,axl*32,test*32,pdb*32,
     1 ins*32,bar*32,parm*32,lnam
      integer*4 opt,nal(n6,7),ial(n6,0:n4),kin(6),kan(6),moi(10),
     1 ord(33,10),nva(3),ipi(0:2),ips(-2:2),ipm(-2:2)
      common/cha/mac,lmo,lib,libm,out,axe,axl,noe,nol,test,pdb,ins,
     1 bar,parm
      common/datjm/acc,phos,epsi,epsr,plat,slope,rhbl,vfac,tfac,rfac,xfac,
     1 scale,rpiv,tpiv,fad,fan,damp,fnoew,fnoes,fnoea,df1,scnb,scee,
     1 catd,catr,catc,rad,pit,enit,nshel,maxn,opt,mhomo,mhomo2,mhomo3,
     1 nick,limit,nop,nrib,ncat,lig,naxo,sup,rcom,homo,homo2,homo3,diep,
     1 sum,link,ecen,cyl,lcat,cent,autos,amber
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/flx/sap(n6),refg,refb,refh,refv,iap(n8,n4,n5),
     1 nap(n8,7,n5),kap(3,n5)
      common/ind/iend(0:4),ise(n2),ienb(n2,2),kapt(n2+1),nsph
      common/lgd/rlig(n9,6),slig(n9,2),ilig(n9,6),lopt(n9),lpiv(n9),
     1 ntlg,nlgi,nlig,ntl,lnam(n9),locr(n2,n0)
      common/lop/dellp(2),stlp(2),gmx,indlp(2,2),lplow(2),lphig(2),
     1 nloop,icy,tlp(2),quiet,convg,dloop
      common/mnn/var(n7),gra(n7),con(n0*n2+n3),scl(n7),ncon,nvrc,
     1 ntba,nbac,nthe,nhel,ntki,nkin,ntri,ntot,nvar,nrin,
     1 lar(n7),lock(n6a)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(120,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      data ord/
     1 -1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-16,17,11,12,13,14,15,-18,-19,
     1 20,21,22,-23,-24,0,0,0,0,0,0,0,2*0,
     1 -1,-2,-3,-4,-5,-6,-7,-8,-9,-10,17,18,16,-19,20,11,12,13,14,
     1 15,-21,-22,23,24,25,-26,-27,28,29,0,0,0,0,
     1  1,2,3,4,5,6,7,8,9,10,17,18,16,19,20,11,12,13,14,15,21,22,23,
     1 24,25,26,27,28,29,0,0,0,0,
     1 -1,-2,-3,-4,-5,-6,7,8,9,10,-17,-18,-16,19,20,-11,12,13,-14,15,
     1 -21,-22,-23,-24,-25,26,27,-28,-29,0,0,0,0,
     1 -1,-2,-3,-4,-5,-6,-12,-13,-14,-7,8,9,-10,11,-15,-16,-17,-18,-19,
     1 -20,-21,0,0,0,0,0,0,0,0,0,0,2*0,
c -------ribose
     1 -1,-2,-3,-4,-5,-6,-7,-8,-9,-10,11,-17,18,12,13,14,15,16,-19,
     1 -20,21,22,23,-24,-25,-26,27,28,0,0,0,2*0,
     1 -1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-17,-18,11,-19,-20,21,12,13,14,15,
     1 16,-22,-23,24,25,26,-27,-28,29,30,-31,32,33,
     1  1,2,3,4,5,6,7,8,9,10,17,18,11,-19,20,21,12,13,14,15,16,22,23,24,
     1     25,26,27,28,29,30,31,32,33,
     1 -1,-2,-3,-4,-5,-6,7,8,9,10,17,18,11,19,20,21,-12,13,14,-15,16,
     1 -22,-23,-24,-25,-26,27,28,-29,-30,-31,32,33,
     1 -1,-2,-3,-4,-5,-6,7,-13,-14,15,-8,9,10,-11,12,-16,-17,-18,-19,
     1 -20,-21,-22,-23,24,25,0,0,0,0,0,0,2*0/
      data kin/12,15,9,13,16,10/,kan/17,20,14,18,20,15/,
     1 moi/3,3,0,14,10,3,3,0,15,11/,
     1 ipi/1,0,2/,ips/1,0,2,0,0/,ipm/0,0,2,-2,0/
         if(.not.link) then
         write(6,*) '  ---- .lmo output requires link=.T. ----'
c         stop
         link=.true.
         endif
      kfi=index(lmo,' ')-1
      open(unit=11,file=lmo(:kfi)//'.lmo',status='new')
c--------------------------------------------------------------loop over strands
      do ks=1,nst
      ist=iend(ks-1)+1
      ien=iend(ks)
      idir=idr(ks)
         if(idir.gt.0) then
         id1=1
         id2=2
         else
         id1=2
         id2=1
         endif
      nring=ien-ist+1
      ksm=nuc(ien)-nuc(ist-1)
      call joie(ist,ien,ksh)
      mof=nuc(ist-1)
      nof=ist-1
      isp=ist+(ien-ist)/2
      in=itr(isp)
      ipv=nuc(isp-1)+iofs(isp)-mof
         if(isp.eq.1) then
         write(6,*)
         write(6,*) '  --- No Lmoout for dinucleotide strands ---'
         stop
         endif
c----------------------------------------------------------------------variables
      n=0
      do ik=1,3
      do is=ist,ien
      iof=nuc(is-1)-mof
      isz=nuc(is)-nuc(is-1)
      if(is.ge.2)isz1=nuc(is-1)-nuc(is-2)
      if(idir.lt.0.and.is.lt.ien)isz1=nuc(is+1)-nuc(is)
      in=itr(is)
      inn=in
      if(inn.gt.3) inn=inn-3
      inp=inn
      if(inn.eq.2.and.is.eq.isp) inp=3
      if(inn.eq.2.and.is.gt.isp) inp=4
         if(idir.lt.0) then
         if(inp.eq.2) then
         inp=4
         else if(inp.eq.4) then
         inp=2
         endif
         endif
      if(inn.eq.3) inp=5
      if(ribose(is)) inp=inp+5
      if(ik.eq.1) then
      lo=1
      lu=kin(in)
      else if(ik.eq.2) then
      lo=kin(in)+1
      lu=kan(in)
      else
      lo=kan(in)+1
      lu=kap(2,in)
      endif
         if(idir.gt.0) then
         dup=-nuc(is-1)+nuc(is)
         if(is.gt.1) dow=-nuc(is-1)+nuc(is-2)
         else
         if(is.gt.1) dup=-nuc(is-1)+nuc(is-2)
         dow=-nuc(is-1)+nuc(is)
         endif
      do l=lo,lu
      n=n+1
      ln=abs(ord(l,inp))
      do j=1,7
      if(j.le.4) then
      jj=j
      if(ord(l,inp).lt.0) jj=5-j
      ix=nap(ln,5,in)
      del=iof
      if(nap(ln,4,in).ne.0) then
      if((ix.ge.1.and.jj.eq.4).or.(ix.eq.2.and.jj.ge.3)) del=iof+dup
      if((ix.le.-1.and.jj.eq.4).or.(ix.eq.-2.and.jj.ge.3)) del=iof+dow
      else
      if((ix.ge.1.and.jj.eq.3).or.(ix.eq.2.and.jj.ge.2)) del=iof+dup
      if((ix.le.-1.and.jj.eq.3).or.(ix.eq.-2.and.jj.ge.2)) del=iof+dow
      endif
      if(nap(ln,jj,in).eq.0) del=0
      nal(n,j)=nap(ln,jj,in)+del
c-------------------------------------------------------------------------------
      else
      nal(n,j)=nap(ln,j,in)
      endif
      enddo
      if(nal(n,1).eq.0) then
      do j=1,3
      nal(n,j)=nal(n,j+1)
      enddo
      nal(n,4)=0
      endif
         if(ik.eq.3) then
         do m=1,nva(2)
         if(nal(m,2).eq.nal(n,2).and.nal(m,3).eq.nal(n,3)
     1   .and.nal(m,4).ne.0) then
         nal(n,6)=-m
         goto 50
         endif
         enddo
 50      endif
c---------------------------------------------------------------------atom moves
      if(ik.le.2) then
      do j=1,iofs(is)-1
      ip=iap(ln,j,in)
      if(ord(l,inp).lt.0) then
      if(nap(ln,2,in).eq.1.and.nap(ln,4,in).eq.0) then
      ip=ips(ip)
      else if(ln.eq.moi(inp)) then
      ip=ipm(ip)
      else
      ip=ipi(ip)
      endif
      if(j.eq.nap(ln,2,in)) ip=0
      if(j.eq.nap(ln,3,in).and.nap(ln,4,in).ne.0) ip=0
      endif
      ial(n,j)=ip
      enddo

      if (ribose(is)) then
      if((ln.eq.18.and.inp.eq.6).or.(ln.eq.21.and.inp.eq.7)) then
         do j=1,iofs(is)-1
         ial(n,j)=-1
         enddo
         endif
      else
         if((ln.eq.17.and.inp.eq.1).or.(ln.eq.20.and.inp.eq.2)) then
         do j=1,iofs(is)-1
         ial(n,j)=-1
         enddo
         endif
      endif
      do j=iofs(is),isz
      ial(n,j)=ial(n,8)
      enddo
      if(ln.eq.4.and.inp.ne.3.and.inp.ne.8) then
       do j=iofs(is),isz
       ial(n,j)=1
       if (j.eq.nap(ln,2,in).or.j.eq.nap(ln,3,in)) ial(n,j)=0
       enddo
      else if(ln.eq.4.and.(inp.eq.3.or.inp.eq.8)) then
       do j=iofs(is),isz
       ial(n,j)=0
       enddo
      endif
      ial(n,isz+id1)=0
      ial(n,isz+id2)=0
      if(inp.eq.2.or.inp.eq.7) then
       ial(n,isz+id1)=ial(n,13)
       ial(n,isz+id2)=0
      else if(inp.eq.3.or.inp.eq.8) then
       ial(n,isz+id1)=ial(n,13)
       ial(n,isz+id2)=ial(n,16)
       if(ribose(is)) then
         if(ln.eq.20.or.ln.eq.21) ial(n,isz+id2)=1
       else
         if(ln.eq.19.or.ln.eq.20) ial(n,isz+id2)=1
       endif
       else if(inp.eq.4.or.inp.eq.9) then
         ial(n,isz+id1)=0
         ial(n,isz+id2)=ial(n,16)
        if (ribose(is)) then
         if(ln.eq.20.or.ln.eq.21) ial(n,isz+id2)=1
        else
         if(ln.eq.19.or.ln.eq.20) ial(n,isz+id2)=1
        endif
       endif
c-----------------------------
      if(inp.eq.2.or.inp.eq.3.or.inp.eq.7.or.inp.eq.8)then
       if(ribose(is)) then
        if(ln.eq.18)then
         ial(n,isz+id1)=1
         ial(n,isz+id2)=0
         ial(n,13)=2
         ial(n,14)=2
        else if(ln.eq.19)then
         do i=1,isz1
         ial(n,i)=1
         enddo
         ial(n,16)=0
         ial(n,isz1+id1)=1
         ial(n,isz1+id2)=0
        endif
       else
        if(ln.eq.17)then
         ial(n,isz+id1)=1
         ial(n,isz+id2)=0
         ial(n,13)=2
         ial(n,14)=2
        else if(ln.eq.18)then
         do i=1,isz1
         ial(n,i)=1
         enddo
         ial(n,16)=0
         ial(n,isz1+id1)=1
         ial(n,isz1+id2)=0
        endif
       endif ! of ribose
      else if(inp.eq.4.or.inp.eq.9) then
       if(ribose(is)) then
        if(ln.eq.18)then
         do i=1,isz
         ial(n,i)=1
         enddo
         ial(n,isz+id1)=0
         ial(n,isz+id2)=1
         ial(n,13)=2
         ial(n,14)=2
        else if(ln.eq.19)then
         do i=1,isz1
         ial(n,i)=0
         enddo
         ial(n,isz1+id2)=1
         ial(n,isz1+id1)=0
        endif
       else
        if(ln.eq.17)then
         ial(n,isz+id2)=1
         ial(n,isz+id1)=0
         do i=1,isz
         ial(n,i)=1
         enddo
         ial(n,13)=2
         ial(n,14)=2
        else if(ln.eq.18)then
         do i=1,isz1
         ial(n,i)=0
         enddo
         ial(n,isz1+id2)=1
         ial(n,isz1+id1)=0
        endif
       endif ! of ribose
      else if(inp.eq.5.or.inp.eq.10) then
       if(ribose(is)) then
        if(ln.eq.14)then
         ial(n,isz+id2)=0
         ial(n,isz+id1)=0
         do i=1,isz
         ial(n,i)=1
         enddo
         ial(n,13)=2
         ial(n,14)=2
         ial(n,6)=0
        else if(ln.eq.15)then
         do i=1,isz1
         ial(n,i)=0
         enddo
         ial(n,isz1+id2)=1
         ial(n,isz1+id1)=0
        endif
       else
        if(ln.eq.13)then
         do i=1,isz
         ial(n,i)=1
         enddo
         ial(n,13)=2
         ial(n,14)=2
         ial(n,6)=0
         ial(n,isz+id2)=0
         ial(n,isz+id1)=0
        else if(ln.eq.14)then
         ial(n,isz1+id2)=1
         ial(n,isz1+id1)=0
         do i=1,isz1
         ial(n,i)=0
         enddo
        endif
       endif ! of ribose
      endif
c-----------------------------
      ial(n,0)=isz+2
      if(inp.ge.2.and.inp.le.4.or.inp.ge.7.and.inp.le.9)then
       if(ribose(is).and.ln.eq.19)ial(n,0)=isz1+2
       if(.not.ribose(is).and.ln.eq.18)ial(n,0)=isz1+2
      else if(inp.eq.5.or.inp.eq.10)then
       if(ribose(is).and.ln.eq.15)ial(n,0)=isz1+2
       if(.not.ribose(is).and.ln.eq.14)ial(n,0)=isz1+2
      endif
      endif
      enddo
         if(ik.eq.1.and.lthy(is)) then
         n=n+1
         do j=1,4
         nal(n,j)=ithy(j)+iofs(is)+iof
         enddo
         nal(n,5)=0
         nal(n,6)=0
         nal(n,7)=nith
         do j=1,isz
         ial(n,j)=0
         enddo
         do j=4,6
         ial(n,ithy(j)+iofs(is))=1
         enddo
         ial(n,isz+id1)=0
         ial(n,isz+id2)=0
         ial(n,0)=isz+2
         endif
      enddo
      nva(ik)=n
      enddo
c-------------------------------------------------------------------------output
      write(11,16) ks
 16   format('STR',i1)
      write(11,18) ksm,ksh,0,0,0,nva(1),nva(3),nva(2),nring,ipv,2
 18   format(11i4)
      write(11,20) (matm(i),i=1,ksh)
 20   format(15i5)
      do is=ist,ien
      iof=iofs(is)+nuc(is-1)
      do i=nuc(is-1)+1,nuc(is)
      write(11,22) mnam(i),(corm(i,j),j=1,3),imch(i),imty(i),
     1 dmon(i),munit(iof),nunit(iof)-nof,icm(i),i-mof
 22   format(a4,3f10.5,2i2,f8.4,1x,a4,1x,i3,1x,2i4)
      enddo
      enddo
      write(11,24) (isr(1,2,is)+nuc(is-1)-mof,
     1 isr(1,1,is)+nuc(is-1)-mof,rsr(1,is),is=ist,ien)
 24   format(3(2i3,f6.3))
      write(11,28) (ial(i,0),i=1,nva(2))
      do j=1,7
      write(11,26) (nal(l,j),l=1,nva(3))
 26   format(20i4)
      enddo
      do j=1,nva(2)
      write(11,28) (ial(j,i),i=1,ial(j,0))
 28   format(25i3)
      enddo
      enddo
      close(11)
      return
      end
