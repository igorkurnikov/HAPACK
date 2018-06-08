      subroutine setvar
      include 'jumna_data.inc'
      character*4 mnam,munit,seq*120,code*8,kode*8,lnam
      logical*2 lthy,lar,lock,kink,ifhb,ribose,cation,locr,ihl,sup,
     1 rcom,homo,homo2,homo3,diep,link,ecen,cyl,lcat,cent,autos,
     1 amber,sum
      integer*4 opt
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
      common/mnn/var(n7),gra(n7),con(n0*n2+n3),scl(n7),ncon,nvrc,
     1 ntba,nbac,nthe,nhel,ntki,nkin,ntri,ntot,nvar,nrin,
     1 lar(n7),lock(n6a)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/parjm/aij(25,25),bij(25,25),ahb(15),bhb(15),vt(35),va(25),
     1 vo(25),delmax,time0,nzsh,nwdg,ivf(35),ifhb(25,25)
      common/sst/phase(8,3),ampli(8,3),ph5(n2),am5(n2),kr5(n2)
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      common/symjm/xsum,ysum,itwl(n2),nvs(n7),ihm(n2),ksym(3),nbrk(2),
     1 isym,isur,isup,ihl(n2)
      k=0
      ki=0
      do is=1,ntl
        kr5(is)=0
        in=itr(is)
        ioff=nuc(is-1)
        idir=0
        if(is.le.nto) idir=idr(ilq(is,2))
        krin=kap(1,in)+nsr(is)*5
        do l=1,kap(2,in)
          k=k+1
          ia=nap(l,1,in)+ioff
          ib=nap(l,2,in)+ioff
          ic=nap(l,3,in)+ioff
          id=nap(l,4,in)+ioff
          ix=nap(l,5,in)
          if(id.eq.ioff) then
c if coordinate is a valence angle
            if(idir.gt.0) then
              if(ix.ge.1) ic=ic-ioff+nuc(is)
              if(ix.eq.2) ib=ib-ioff+nuc(is)
              if(ix.le.-1) ic=ic-ioff+nuc(is-2)
              if(ix.eq.-2) ib=ib-ioff+nuc(is-2)
            else
              if(ix.ge.1) ic=ic-ioff+nuc(is-2)
              if(ix.eq.2) ib=ib-ioff+nuc(is-2)
              if(ix.le.-1) ic=ic-ioff+nuc(is)
              if(ix.eq.-2) ib=ib-ioff+nuc(is)
            endif
            sap(k)=ang(ia,ib,ic)
          else
c if coordinate is a dihedral angle 
		  if(l.gt.krin) then
			if(idir.gt.0) then
			  if(ix.ge.1) id=id-ioff+nuc(is)
			  if(ix.eq.2) ic=ic-ioff+nuc(is)
			  if(ix.le.-1) id=id-ioff+nuc(is-2)
			  if(ix.eq.-2) ic=ic-ioff+nuc(is-2)
			else
			  if(ix.ge.1) id=id-ioff+nuc(is-2)
			  if(ix.eq.2) ic=ic-ioff+nuc(is-2)
			  if(ix.le.-1) id=id-ioff+nuc(is)
			  if(ix.eq.-2) ic=ic-ioff+nuc(is)
			endif
		  endif
		  sap(k)=torp(ia,ib,ic,id)
	    endif
c-------------------------------------------------------------position variables
        if(.not.lock(k).and.l.le.kap(1,in)) then
          ki=ki+1
          if(ki.gt.n7) then
            write(6,*) '  ---- n7 too small ----'
            stop
          endif
          var(ki)=sap(k)
          lar(ki)=.false.
          if(id.eq.ioff) lar(ki)=.true.
          if(is.le.nto.and.l.eq.5) kr5(is)=ki
        endif
        enddo
        if(is.eq.nto) nsph=ki
      enddo
      ntba=k
      nbac=ki
c--------------------------------------------------------------helical variables
      do is=1,nto
        do l=1,6
          k=k+1
          if(.not.lock(k)) then
            ki=ki+1
            var(ki)=hel(is,l)
            if(l.gt.3) then
              lar(ki)=.true.
            else
              lar(ki)=.false.
            endif
          endif
        enddo
      enddo
      nthe=k
      nhel=ki
c--------------------------------------------------------------------------kinks
      do is=1,kseq
      if(kink(is)) then
      do l=1,4
      k=k+1
      if(.not.lock(k)) then
      ki=ki+1
      var(ki)=vkink(is,l)
      if(l.le.2) then
      lar(ki)=.false.
      else
      lar(ki)=.true.
      endif
      do kp=2,nst
      js=ieq(is,kp)
      enddo
      endif
      enddo
      endif
      enddo
      ntki=k
      nkin=ki
c-------------------------------------------------------------------ligand inter
      do i=1,nlig
      kal=nuc(i+nto)-nuc(i+nto-1)
      lint=6
      if(kal.eq.1) lint=3
      do j=1,lint
      k=k+1
      if(.not.lock(k)) then
      ki=ki+1
      var(ki)=rlig(i,j)
      lar(ki)=.false.
      if(j.gt.3) lar(ki)=.true.
      endif
      enddo
      enddo
      ntlg=k
      nlgi=ki
c-----------------------------------------------------------------thymine methyl
      kth=ntba
      do is=1,nto
        ks=(ilq(is,2)-1)*kseq+ilq(is,1)
        if(lthy(is)) then
          in=itr(is)
          ioff=nuc(is-1)
          k=k+1
          kth=kth+1
          ia=ithy(1)+iofs(is)+ioff
          ib=ithy(2)+iofs(is)+ioff
          ic=ithy(3)+iofs(is)+ioff
          id=ithy(4)+iofs(is)+ioff
          sap(kth)=torp(ia,ib,ic,id)
c-------------------------------------------------------------------reset methyl
          del=thy(is)-sap(kth)
          sap(kth)=thy(is)
          x0=corm(ib,1)
          y0=corm(ib,2)
          z0=corm(ib,3)
          rx=corm(ic,1)-x0
          ry=corm(ic,2)-y0
          rz=corm(ic,3)-z0
          rr=sqrt(rx*rx+ry*ry+rz*rz)
          rx=rx/rr
          ry=ry/rr
          rz=rz/rr
          ca=cos(cdr*(del))
          sa=sin(cdr*(del))
          do j=4,6
            jg=ithy(j)+iofs(is)+ioff
            xx=corm(jg,1)-x0
            yy=corm(jg,2)-y0
            zz=corm(jg,3)-z0
            corm(jg,1)=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     1          (rx*rz*(1-ca)+ry*sa)*zz+x0
            corm(jg,2)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1          (ry*rz*(1-ca)-rx*sa)*zz+y0
            corm(jg,3)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1          (rz*rz+(1-rz*rz)*ca)*zz+z0
          enddo
          if(.not.lock(k)) then
            ki=ki+1
            var(ki)=sap(kth)
            lar(ki)=.true.
          endif
      endif
      enddo
c---------------------------------------------------------------------superhelix
      if(nshel.eq.1.or.nshel.eq.3) then
        ki=ki+1
        var(ki)=rad
        lar(ki)=.false.
        isur=ki
      endif
      if(nshel.eq.2.or.nshel.eq.3) then
         ki=ki+1
         var(ki)=pit
         lar(ki)=.false.
         isup=ki
      endif
c-------------------------------------------------------------------------------
      ntot=k
      nvar=ki
      return
      end
