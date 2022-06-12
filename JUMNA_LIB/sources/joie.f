      subroutine joie(ist,ien,l)
      include 'jumna_data.inc'
      logical*2 kink,lthy,ribose,cation,locr
      character*4 mnam,munit,snam,suni,sub,seq*120,lnam,
     1 mac*32,lmo*32,lib*32,libm*32,out*32,axe*32,noe*32,nol*32,axl*32,
     1 test*32,pdb*32,ins*32,bar*32,parm*32,code*8,kode*8
      common/cha/mac,lmo,lib,libm,out,axe,axl,noe,nol,test,pdb,ins,
     1 bar,parm
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/flx/sap(n6),refg,refb,refh,refv,iap(n8,n4,n5),
     1 nap(n8,7,n5),kap(3,n5)
      common/lgd/rlig(n9,6),slig(n9,2),ilig(n9,6),lopt(n9),lpiv(n9),
     1 ntlg,nlgi,nlig,ntl,lnam(n9),locr(n2,n0)
      common/moljm/sor(n4,n5,3),smon(n4,n5),snam(n4,n5),suni(n4,n5),
     1 nuni(n4,n5),sub(n5),isch(n4,n5),isty(n4,n5),ics(n4,n5),
     1 mats(3*n4,n5),kas(n5),khs(n5),ksub
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(120,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      mof=nuc(ist-1)
      l=0
      if(ist.gt.nto) goto 100
c--------------------------------------------------------------make linear b.mat
      do is=ist,ien
      in=itr(is)
      ioff=nuc(is-1)-mof
      iu=irec(is)
      istr=ilq(is,2)
      idir=sign(1,idr(istr))
      do m=1,khs(iu)
      l=l+1
      me=mats(m,iu)
      if(me.ge.10000) then
      matm(l)=me
      else if(me.ne.0) then
      idelt=0
      if((ribose(is).and.parm.eq.'Flex').and.in.le.6) idelt=1
      if(abs(me).lt.iofs(is)-idelt) idelt=0
      matm(l)=me+sign(ioff+idelt,me)
      else
      if(idir.eq.1)  matm(l)=nuc(is-2)-mof+16
      if(idir.eq.-1) matm(l)=nuc(is)-mof+16
      ml1=matm(l-1)
      ml2=matm(l)
      endif
      enddo
c---------------------------------------------------------------ribose extension
         if((ribose(is).and.parm.eq.'Flex').and.in.le.6) then
         l=l+1
         matm(l)=10000
         l=l+1
         matm(l)=ioff+10
         l=l+1
         matm(l)=ioff+iofs(is)-1
         endif
      l=l+1
      matm(l)=10000
      enddo
      matm(l)=0
      l=l-1
      return
c------------------------------------------------------------------------Ligands
 100  is=ist
      in=itr(is)
      ioff=nuc(is-1)-mof
      iu=irec(is)
      do m=1,khs(iu)
      l=l+1
      me=mats(m,iu)
      if(me.ge.10000) then
      matm(l)=me
      else
      matm(l)=me+sign(ioff,me)
      endif
      enddo
      return
      end
