      subroutine energ91
      include 'jumna_data.inc'
      logical*2 ifhb,kink,lthy,hst,bst,vst,lgi,lgj,lgu,ribose,
     1 lock,sup,rcom,homo,homo2,homo3,diep,sum,link,ecen,cyl,lcat,
     1 cent,autos,locr,lar,cation,amber,mfhb,i14
      character*4 mnam,munit,seq*120,code*8,kode*8,lnam,atclas*2
      integer*2 i23,i34,elim,iitor(88)
      integer*4 opt,fold
      dimension etork(88)
      common/cnd/rcyl,dcyl,hcyl,tcyl,eneq,ucyl(3),pcyl(3),ecyl,
     1 ac(25),bc(25),lcyl(3)
      common/datjm/acc,phos,epsi,epsr,plat,slope,rhbl,vfac,tfac,rfac,xfac,
     1 scale,rpiv,tpiv,fad,fan,damp,fnoew,fnoes,fnoea,df1,scnb,scee,
     1 catd,catr,catc,rad,pit,enit,nshel,maxn,opt,mhomo,mhomo2,mhomo3,
     1 nick,limit,nop,nrib,ncat,lig,naxo,sup,rcom,homo,homo2,homo3,diep,
     1 sum,link,ecen,cyl,lcat,cent,autos,amber
      common/enf/for(n1,3),tor(n1,3),fot(n6),ener,elec,
     1 repl,disp,eang,etog,epen,i23(8*n1),i34(8*n1),elim(n1)
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/flx/sap(n6),refg,refb,refh,refv,iap(n8,n4,n5),
     1 nap(n8,7,n5),kap(3,n5)
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
      common/para/ro(57),reps(57),angc(191),fangc(191),btor(88),
     $ gtor(88),ftor(88),bdtor(30),gdtor(30),fdtor(30),amb(57,57),
     $ bmb(57,57),amh(9,11),bmh(9,11),mtor(88),itor(88,4),iadt(30,4),
     $ iangc(191,3),mof,nof,nang,ntor,ntors,nad,atclas(57)
      common/slg/hst(n2,6),bst(n2,n8),vst(n2,4),lgi(n1),lgj(n1),lgu(n2)
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      do i=1,ntor
      iitor(i)=0
      etork(i)=0.
      enddo
      do j=1,3
      do i=1,kam
      for(i,j)=0.
      tor(i,j)=0.
      enddo
      enddo
      do i=1,n6
      fot(i)=0.
      enddo
      repl=0.
      disp=0.
      elec=0.
      eang=0.
      etog=0.
      ehyd=0.
c      if(etog.eq.0) goto 600
      k=0
      kk=0
      do i=1,kam-1
      if(.not.lgj(i)) goto 100
      k=k+1
      kk=kk+1
      m=imty(i)
      dmi=dmon(i)
      ini=elim(i)
      do j=i+1,kam
      if(.not.lgj(j).or.(.not.lgi(i).and..not.lgi(j))) goto 200
      if(ini.ne.0.and.elim(j).eq.ini) goto 200
      do l=k+1,k+i23(k)
      if(j.eq.i23(l)) goto 200
      enddo
      face=1.
      facn=1.
      i14=.false.
      do l=kk+1,kk+i34(kk)
      if(j.eq.i34(l)) then
      face=scee
      facn=scnb
      i14=.true.
      endif
      enddo
      n=imty(j)
      dmj=dmon(j)
      mfhb=.false.
      if((m.eq.19.or.m.eq.24.or.m.eq.17).and.(n.eq.33.or.
     $ n.eq.34.or.n.ge.39.and.n.le.43))mfhb=.true.
      if((n.eq.19.or.n.eq.24.or.n.eq.17).and.(m.eq.33.or.
     $ m.eq.34.or.m.ge.39.and.m.le.43))mfhb=.true.
c-----------------------------------------------------------------electrostatics
      dx=corm(j,1)-corm(i,1)
      dy=corm(j,2)-corm(i,2)
      dz=corm(j,3)-corm(i,3)
      r2=dx*dx+dy*dy+dz*dz
      r=sqrt(r2)
      dx=dx/r2
      dy=dy/r2
      dz=dz/r2
      if(diep) then
      alr=slope*r
      ex=exp(-alr)
      alr2=alr*alr
      e=plat-(plat-enit)*(alr2+2.*alr+2.)*ex/2.
      g=(plat-enit)*slope*alr2*ex/(2.*e*e)
      re=r*e
      sgf=convk*dmi*dmj*g/face
      else if(epsr.ne.0.0) then
      re=epsr*r2
      srij=convk*dmi*dmj*2/(re*face)
      else
      re=r*epsi
      endif
      sij=convk*dmi*dmj/(re*face)
      elec=elec+sij
      if(diep) then
      ss=sij+sgf
      for(i,1)=for(i,1)-ss*dx
      for(i,2)=for(i,2)-ss*dy
      for(i,3)=for(i,3)-ss*dz
      for(j,1)=for(j,1)+ss*dx
      for(j,2)=for(j,2)+ss*dy
      for(j,3)=for(j,3)+ss*dz
      else if(epsr.ne.0.0) then
      for(i,1)=for(i,1)-srij*dx
      for(i,2)=for(i,2)-srij*dy
      for(i,3)=for(i,3)-srij*dz
      for(j,1)=for(j,1)+srij*dx
      for(j,2)=for(j,2)+srij*dy
      for(j,3)=for(j,3)+srij*dz
      else
      for(i,1)=for(i,1)-sij*dx
      for(i,2)=for(i,2)-sij*dy
      for(i,3)=for(i,3)-sij*dz
      for(j,1)=for(j,1)+sij*dx
      for(j,2)=for(j,2)+sij*dy
      for(j,3)=for(j,3)+sij*dz
      endif
c------------------------------------------------------------------lennard-jones
      r6=r**6
      if(m.le.0.or.m.gt.53) write(6,*)' m= ',m,' i=',i,mnam(i)
      if(n.le.0.or.n.gt.53) write(6,*)' n= ',n,' j=',j,mnam(j)
      if(mfhb.and.r.le.rhbl.and..not.i14) then
      r10=r6*r2*r2
      mm=min(m,n)-mof
      nn=max(m,n)-nof
      e10=bmh(mm,nn)/r10
      e12=amh(mm,nn)/(r6*r6)
      ee=12.*e12-10.*e10
      disp=disp-e10
      repl=repl+e12
      ehyd=ehyd-e10+e12
      for(i,1)=for(i,1)-ee*dx
      for(i,2)=for(i,2)-ee*dy
      for(i,3)=for(i,3)-ee*dz
      for(j,1)=for(j,1)+ee*dx
      for(j,2)=for(j,2)+ee*dy
      for(j,3)=for(j,3)+ee*dz
      else
c---------------------------------------------------------------no hydrogen bond
      e6=bmb(m,n)/(r6*facn)
      e12=amb(m,n)/(r6*r6*facn)
      ee=12.*e12-6.*e6
      disp=disp-e6
      repl=repl+e12
      for(i,1)=for(i,1)-ee*dx
      for(i,2)=for(i,2)-ee*dy
      for(i,3)=for(i,3)-ee*dz
      for(j,1)=for(j,1)+ee*dx
      for(j,2)=for(j,2)+ee*dy
      for(j,3)=for(j,3)+ee*dz
      endif
200   enddo
      k=k+i23(k)
      kk=kk+i34(kk)
100   enddo
c--------------------------------------------------torsional/angle strain energy
300   k=0
      ks=0
      kz=0
      do is=1,ntl
      in=itr(is)
      ioff=nuc(is-1)
      idir=0
      if(is.le.nto) idir=idr(ilq(is,2))
      krin=kap(1,in)+nsr(is)*5
      if(.not.lgu(is)) then
      k=k+kap(2,in)
      goto 500
      endif
      do l=1,kap(2,in)
      k=k+1
      ia=nap(l,1,in)+ioff
      ib=nap(l,2,in)+ioff
      ic=nap(l,3,in)+ioff
      id=nap(l,4,in)+ioff
      ik=nap(l,6,in)
      it=nap(l,7,in)
      if(id.eq.ioff) then
      ix=0
      if(l.gt.krin) then
      ix=nap(l,5,in)
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
      endif
      fcon=fangc(it)
      angk=ang(ia,ib,ic)
c      if((in.eq.2.and.ik.eq.-17).or.(in.eq.3.and.ik.eq.-13).or.
c     1   (in.eq.5.and.ik.eq.-18).or.(in.eq.6.and.ik.eq.-14)) then
c      angk=ang(ia,ib,ic)
      dang=(angk-angc(it))*cdr
c      else
c      if(abs(angk-sap(k)).gt.1.d-7) then
c      write(6,*)' k= ',k,' angk= ',angk,' sap ',sap(k)
c      write(6,*)mnam(ia),mnam(ib),mnam(ic),
c     *nunit(ia),nunit(ib),nunit(ic),
c     *(nap(l,ii,in),ii=1,7),in
c      endif
c      dang=(sap(k)-angc(it))*cdr
c      endif
      eang=eang+fcon*dang**2
      fcon=2.*fcon*dang
      call delval(fcon,ia,ib,ic)
      else
      fold=ftor(it)
      barr=btor(it)
      ix=0
      if(l.gt.krin) then
      ix=nap(l,5,in)
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
      tork=torp(ia,ib,ic,id)
c      if((in.eq.2.and.(ik.eq.-20.or.ik.eq.-16)).or.
c     1   (in.eq.1.and.ik.eq.-17).or.
c     1   (in.eq.5.and.(ik.eq.-21.or.ik.eq.-17)).or.
c     1   (in.eq.4.and.ik.eq.-18)) then
      etog=etog+
     $ barr*(1.+cos(cdr*(fold*tork-gtor(it))))/mtor(it)
      fcon=barr*fold*sin(cdr*(fold*tork-gtor(it)))/mtor(it)
c      else
c      if(abs(tork-sap(k)).gt.1.d-7) then
c      write(6,*)' k= ',k,' tork= ',tork,' sap ',sap(k)
c      write(6,*)mnam(ia),mnam(ib),mnam(ic),mnam(id),
c     *nunit(ia),nunit(ib),nunit(ic),nunit(id),
c     *(nap(l,ii,in),ii=1,7),in
c      endif
c      etog=etog+
c     $ barr*(1.+cos(cdr*(fold*sap(k)-gtor(it))))/mtor(it)
c      fcon=barr*fold*sin(cdr*(fold*sap(k)-gtor(it)))/mtor(it)
c      endif
      call deltor(tork,fcon,ia,ib,ic,id)
      endif
      enddo
500   ks=k
      enddo
c-----------------------------------------------------------------thymine methyl
c
c-------------------------------------------------------------------------------
600   call penalty
      if(cyl) call enecyl(cyl)
      ener=elec+repl+disp+eang+etog+epen
      return
      end
