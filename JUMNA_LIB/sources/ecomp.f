      subroutine ecomp
      include 'jumna_data.inc'
      parameter (nhbx=200,ncox=50,emin=5.d0)
      logical*2 ifhb,kink,lthy,lar,lock,hst,bst,vst,lgi,lgj,lgu,
     1 ribose,sup,rcom,homo,homo2,homo3,diep,link,ecen,cyl,lcat,
     1 autos,start,ihl,locr,itype,jtype,cation,icat,jcat,tlp,quiet,
     1 sum,cent,convg,dloop,amber
      character*4 mnam,munit,seq*120,code*8,kode*8,base*1,lnam
      integer*2 i23,i34,elim
      integer*4 opt,fold,lex(29),nun(n1)
      dimension ec(n2,n2,16,0:5),en(n2,n2,0:5),es(5,5,10,5),esv(10),
     1 ed(n2,n2,5),esn(10),esl(10),esh(5,5),index(16),ese(5,5,5),
     1 esta(16),ihb(nhbx,2),rhb(nhbx,2),ihc(ncox,2),rhc(ncox,2)
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
      common/ind/iend(0:4),ise(n2),ienb(n2,2),kapt(n2+1),nsph
      common/lop/dellp(2),stlp(2),gmx,indlp(2,2),lplow(2),lphig(2),
     1 nloop,icy,tlp(2),quiet,convg,dloop
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
      common/slg/hst(n2,6),bst(n2,n8),vst(n2,4),lgi(n1),lgj(n1),lgu(n2)
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      common/symjm/xsum,ysum,itwl(n2),nvs(n7),ihm(n2),ksym(3),nbrk(2),
     1 isym,isur,isup,ihl(n2)
      data lex/15*0,1,1,1,-1,-1,0,0,0,0,0,0,-1,1,1/
      data index/1,2,3,4,2,5,6,7,3,6,8,9,4,7,9,10/
      repl=0.
      disp=0.
      elec=0.
      eang=0.
      etog=0.
      do is=1,ntl
      do i=nuc(is-1)+1,nuc(is)
      nun(i)=is
      enddo
      enddo
c-------------------------------------------------------------------------------
      cang=cos(cdr*(anglim))
      nhb=0
      nco=0
      if(ecen) call pairs(mp)
      do i=1,ntl
      do j=1,ntl
      do k=1,16
      do l=0,5
      ec(i,j,k,l)=0.
      enddo
      enddo
      enddo
      enddo
      k=0
      do i=1,kam-1
      k=k+1
      m=imty(i)
      iun=nun(i)
      dmi=dmon(i)
      ini=elim(i)
      ind=3
      if(munit(i).eq.'PHOS'.or.munit(i).eq.'phos'.or.
     1 munit(i).eq.'PHOL') ind=1
      if(munit(i).eq.'SUCR'.or.munit(i).eq.'sucr'.or.
     1 munit(i).eq.'SUCL') ind=2
      if(i.gt.nuc(iun-1)+iofe(iun)) ind=4
      if(cation(iun).and.i.eq.nuc(iun)) ind=1
      if(iun.gt.nto) ind=4
      do j=i+1,kam
      if(ini.ne.0.and.elim(j).eq.ini) goto 200
      do l=k+1,k+i23(k)
      if(j.eq.i23(l)) goto 200
      enddo
      n=imty(j)
      dmj=dmon(j)
      jun=nun(j)
      jnd=3
      if(munit(j).eq.'PHOS'.or.munit(j).eq.'phos'.or.munit(j).
     1 eq.'PHOL') jnd=1
      if(munit(j).eq.'SUCR'.or.munit(j).eq.'sucr'.or.munit(j).
     1 eq.'SUCL') jnd=2
      if(j.gt.nuc(jun-1)+iofe(jun)) jnd=4
      if(cation(jun).and.j.eq.nuc(jun)) jnd=1
      if(jun.gt.nto) jnd=4
      knd=4*(ind-1)+jnd
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
      re=r*e
      else if (epsr.ne.0.0) then
      re=epsr*r2
      else
      re=r*epsi
      endif
      sij=convk*dmi*dmj/re
      ec(iun,jun,knd,2)=ec(iun,jun,knd,2)+sij
c------------------------------------------------------------------lennard-jones
      r6=r**6
      e6=aij(m,n)/r6
      e12=bij(m,n)/(r6*r6)
      ee=12.*e12-6.*e6
      if(ifhb(m,n).and.r.le.rhbl) then
      if(m.le.2) then
      ihy=i
      sid=1.
      iatype=n
      else
      ihy=j
      sid=-1.
      iatype=m
      endif
      mp=matd(ihy,1)
      idtype=imty(mp)
      if(iatype.eq.18) iatype=11
      if(idtype.eq.18) then
      ih=iatype-6
      else if(idtype.eq.8) then
      ih=iatype-1
      else
      ih=iatype+4
      endif
      qx=corm(ihy,1)-corm(mp,1)
      qy=corm(ihy,2)-corm(mp,2)
      qz=corm(ihy,3)-corm(mp,3)
      q=sqrt(qx*qx+qy*qy+qz*qz)
      cth=sid*(dx*qx+dy*qy+dz*qz)*r/q
      if(cth.lt.0.) then
      ec(iun,jun,knd,0)=ec(iun,jun,knd,0)+e12
      ec(iun,jun,knd,1)=ec(iun,jun,knd,1)-e6+e12
      if(nco.lt.50.and.e12-e6.gt.conlim) then
      nco=nco+1
      ihc(nco,1)=i
      ihc(nco,2)=j
      rhc(nco,1)=r
      rhc(nco,2)=e12-e6
      endif
      else
      adel=(ahb(ih)-aij(m,n))/r6
      bdel=(bhb(ih)-bij(m,n))/(r6*r6)
      if(nhb.lt.nhbx.and.cth.ge.cang.and.r.le.rhblim) then
      nhb=nhb+1
      ihb(nhb,1)=i
      ihb(nhb,2)=j
      rhb(nhb,1)=r
      rhb(nhb,2)=acos(cth)*crd
      endif
      ec(iun,jun,knd,0)=ec(iun,jun,knd,0)+e12+cth*bdel
      ec(iun,jun,knd,1)=ec(iun,jun,knd,1)-e6+e12+cth*(bdel-adel)
      endif
      else
c---------------------------------------------------------------no hydrogen bond
      ec(iun,jun,knd,0)=ec(iun,jun,knd,0)+e12
      ec(iun,jun,knd,1)=ec(iun,jun,knd,1)-e6+e12
      if(nco.lt.ncox.and.e12-e6.gt.conlim) then
      nco=nco+1
      ihc(nco,1)=i
      ihc(nco,2)=j
      rhc(nco,1)=r
      rhc(nco,2)=e12-e6
      endif
      endif
200   enddo
      k=k+i23(k)
      enddo
c--------------------------------------------------torsional/angle strain energy
      k=0
      do is=1,ntl
      in=itr(is)
      do l=1,kap(2,in)
      k=k+1
      ix=nap(l,2,in)+nuc(is-1)
      iy=nap(l,3,in)+nuc(is-1)
      it=nap(l,7,in)
      if((is.eq.nick.and.lex(l).eq.1).or.
     1  (is+1.eq.nick.and.lex(l).eq.-1)) goto 400
      iun=nun(ix)
      ind=3
      if(munit(ix).eq.'PHOS'.or.munit(ix).eq.'phos'.or.munit(ix).eq.
     1 'PHOL') ind=1
      if(munit(ix).eq.'SUCR'.or.munit(ix).eq.'sucr'.or.munit(ix).eq.
     1 'SUCL') ind=2
      if(ix.gt.nuc(iun-1)+iofe(iun)) ind=4
      if(iun.gt.nto) ind=4
      if(nap(l,4,in).eq.0) then
      fcon=va(it)
      dang=(sap(k)-vo(it))*cdr
      ea=fcon*dang**2
      knd=5*ind-4
      ec(iun,iun,knd,3)=ec(iun,iun,knd,3)+ea
      else
      jun=nun(iy)
      jnd=3
      if(munit(iy).eq.'PHOS'.or.munit(iy).eq.'phos'.or.
     1 munit(iy).eq.'PHOL') jnd=1
      if(munit(iy).eq.'SUCR'.or.munit(iy).eq.'sucr'.or.
     1 munit(iy).eq.'SUCL') jnd=2
      if(iy.gt.nuc(jun-1)+iofe(jun)) jnd=4
      if(jun.gt.nto) jnd=4
      fold=ivf(it)
      barr=vt(it)
      et=barr*(1.+sign(1,fold)*cos(cdr*(abs(fold)*sap(k))))/2.
      knd=4*(ind-1)+jnd
      if(iun.le.jun) then
      ec(iun,jun,knd,4)=ec(iun,jun,knd,4)+et
      else
      ec(jun,iun,knd,4)=ec(jun,iun,knd,4)+et
      endif
      endif
400   enddo
      enddo
c-----------------------------------------------------------------thymine methyl
      fold=ivf(nith)
      barr=vt(nith)
      do is=1,nto
      ks=(ilq(is,2)-1)*kseq+ilq(is,1)
      base=seq(ks:ks)
      if(lthy(is)) then
      k=k+1
      et=barr*(1.+sign(1,fold)*cos(cdr*(abs(fold)*sap(k))))/2.
      ec(is,is,11,4)=ec(is,is,11,4)+et
      endif
      enddo
c---------------------------------------------------------ecen nucleotide output
      do is=1,ntl
      do js=is,ntl
      do l=0,5
      en(is,js,l)=0.
      enddo
      do k=1,16
      tot=0.
      do l=0,4
      if(l.gt.0) tot=tot+ec(is,js,k,l)
      en(is,js,l)=en(is,js,l)+ec(is,js,k,l)
      enddo
      ec(is,js,k,5)=tot
      enddo
      do l=1,4
      en(is,js,5)=en(is,js,5)+en(is,js,l)
      enddo
      enddo
      enddo
         if(ecen) then
         do is=1,nto
         if(ilq(is,1).eq.1) goto 51
         do js=is,nto
         if(ilq(js,1).eq.1) goto 52
         if(.not.lgu(is).and..not.lgu(js)) goto 52
         repl=repl+en(is,js,0)
         disp=disp+en(is,js,1)-en(is,js,0)
         elec=elec+en(is,js,2)
         eang=eang+en(is,js,3)
         etog=etog+en(is,js,4)
52       enddo
51       enddo
         endif
      do is=1,ntl
      do js=is,ntl
      do l=1,5
      ed(is,js,l)=en(is,js,l)
      enddo
      enddo
      enddo
c      if(nto.gt.0.and..not.quiet) write(6,50)
c50    format(/2x,'Ecen by nucleotide ....',
c     1//4x,'Nuc',10x,'Elj ',7x,'Eel ',7x,'Eang',7x,'Etog',6x,'Total'/)
      do is=1,ntl
c      if(ise(is).lt.0.and.is.ne.1.and..not.quiet) write(6,8)
c8     format(6x,
c     1'-----------------------------------------------------------')
      ks=(ilq(is,2)-1)*kseq+ilq(is,1)
      if(ks.gt.0.and.ks.le.nto) base=seq(ks:ks)
      do k=0,5
      do js=1,ntl
      if(js.gt.is) then
      en(is,is,k)=en(is,is,k)+en(is,js,k)/2.
      else if(js.lt.is) then
      en(is,is,k)=en(is,is,k)+en(js,is,k)/2.
      endif
      enddo
      enddo
c      if(is.le.nto.and..not.quiet) write(6,60) is,base,
c     1 (en(is,is,l),l=1,5)
      enddo
      if(nst.gt.1) then
      if(.not.quiet) write(6,70)
70    format(/2x,'Ecen by nucleotide levels ....',
     1//4x,'Nuc',10x,'Elj ',7x,'Eel ',7x,'Eang',7x,'Etog',6x,'Total'/)
      do is=1,kseq
      base=seq(is:is)
      elj=0.
      ele=0.
      ean=0.
      eto=0.
      ene=0.
      kpn=0
      do kp=1,nst
      js=ieq(is,kp)
      if(js.ne.0) then
      kpn=kpn+1
      elj=elj+en(js,js,1)
      ele=ele+en(js,js,2)
      ean=ean+en(js,js,3)
      eto=eto+en(js,js,4)
      ene=ene+en(js,js,5)
      endif
      enddo
      if(.not.quiet) write(6,60) is,base,elj,ele,ean,eto,ene,kpn
60    format(2x,i3,':',a2,2x,5f11.3,:,' (',i1,' bases)')
      enddo
      endif
c-------------------------------------------------------------unit decomp output
      nt=nst
      if(nlig.gt.0) nt=nt+1
      do i=1,nt
      do j=i,nt
      esh(i,j)=0
      do k=1,10
      do l=1,5
      es(i,j,k,l)=0.
      enddo
      enddo
      enddo
      enddo
      do i=1,ntl
      if(i.le.nto) then
      ib=ilq(i,2)
      else
      ib=nst+1
      endif
      do j=i,ntl
      if(j.le.nto) then
      jb=ilq(j,2)
      else
      jb=nst+1
      endif
      do k=1,16
      ks=index(k)
      do l=1,5
      es(ib,jb,ks,l)=es(ib,jb,ks,l)+ec(i,j,k,l)
      enddo
      enddo
      enddo
      enddo
      write(6,10)
10    format(/2x,'DEC',5x,'PP',5x,'PS',5x,'PB',5x,'PL',5x,
     1 'SS',5x,'SB',5x,'SL',5x,'BB',5x,'BL',5x,'LL',5x,'Tot'/)
      etot=0.
      enuc=0.
      elig=0.
      do k=1,10
      esv(k)=0.
      esn(k)=0.
      esl(k)=0.
      do i=1,nt
      do j=i,nt
      esh(i,j)=esh(i,j)+es(i,j,k,5)
      if(i.ne.nt.and.j.ne.nt) then
      esn(k)=esn(k)+es(i,j,k,5)
      enuc=enuc+es(i,j,k,5)
      else if(i.ne.nt.or.j.ne.nt) then
      esl(k)=esl(k)+es(i,j,k,5)
      elig=elig+es(i,j,k,5)
      endif
      esv(k)=esv(k)+es(i,j,k,5)
      etot=etot+es(i,j,k,5)
      enddo
      enddo
      enddo
      do i=1,nst
      do j=i,nst
      write(6,20) i,j,(es(i,j,k,5),k=1,10),esh(i,j)
20    format(2x,i1,'-',i1,1x,11f7.1)
      enddo
      enddo
      if(nst.gt.1.and.nlig.ne.0) write(6,30) (esn(k),k=1,10),enuc
30    format(2x,'----',1x,77('-'),/2x,'Nuc',1x,11f7.1)
      if(nlig.gt.0) then
      if(nst.gt.1) then
      write(6,31)
31    format(2x,'----',1x,77('-'))
      do i=1,nst
      write(6,21) i,(es(i,nt,k,5),k=1,10),esh(i,nt)
21    format(2x,i1,'-','L',1x,11f7.1)
      enddo
      endif
      write(6,22) (esl(k),k=1,10),elig
22    format(2x,'----',1x,77('-'),/2x,'L-N',1x,11f7.1)
      write(6,32) (es(nt,nt,k,5),k=1,10),esh(nt,nt)
32    format(2x,'----',1x,77('-'),/2x,'Lig',1x,11f7.1)
      endif
      if(nst.gt.1.or.nlig.ne.0) write(6,33) (esv(k),k=1,10),etot
33    format(2x,'----',1x,77('-'),/2x,'Tot',1x,11f7.1)
c-----------------------------------------------------------energy decomp output
      do i=1,nt
      do j=i,nt
      do l=1,5
      tot=0.
      do k=1,10
      tot=tot+es(i,j,k,l)
      enddo
      ese(i,j,l)=tot
      enddo
      enddo
      enddo
      write(6,101)
101   format(/2x,'DEC',5x,'Lj',4x,'Elec',4x,'Ang',3x,'Torg',4x,'Tot'/)
      do l=1,5
      esn(l)=0.
      esl(l)=0.
      esv(l)=0.
      do i=1,nt
      do j=i,nt
      if(i.ne.nt.and.j.ne.nt) then
      esn(l)=esn(l)+ese(i,j,l)
      else if(i.ne.nt.or.j.ne.nt) then
      esl(l)=esl(l)+ese(i,j,l)
      endif
      esv(l)=esv(l)+ese(i,j,l)
      enddo
      enddo
      enddo
      do i=1,nst
      do j=i,nst
      write(6,20) i,j,(ese(i,j,l),l=1,5)
      enddo
      enddo
      if(nst.gt.1.and.nlig.ne.0) write(6,301) (esn(l),l=1,5)
301   format(2x,'----',1x,35('-'),/2x,'Nuc',1x,5f7.1)
      if(nlig.gt.0) then
      if(nst.gt.1) then
      write(6,311)
311   format(2x,'----',1x,35('-'))
      do i=1,nst
      write(6,211) i,(ese(i,nt,l),l=1,5)
211   format(2x,i1,'-','L',1x,5f7.1)
      enddo
      endif
      write(6,222) (esl(l),l=1,5)
222   format(2x,'----',1x,35('-'),/2x,'L-N',1x,5f7.1)
      write(6,322) (ese(nt,nt,l),l=1,5)
322   format(2x,'----',1x,35('-'),/2x,'Lig',1x,5f7.1)
      endif
      if(nst.gt.1.or.nlig.ne.0) write(6,333) (esv(l),l=1,5)
333   format(2x,'----',1x,35('-'),/2x,'Tot',1x,5f7.1)
c------------------------------------modified nucleotide and ligand interactions
      start=.true.
      do is=1,ntl
      in=itr(is)
      ino=ito(is)
      icat=.false.
      if(nuc(is)-nuc(is-1).eq.1) icat=.true.
      itype=.false.
      if(in.ne.ino.or.is.gt.nto) itype=.true.
      ia=nuc(is-1)+iofs(is)
      if(is.gt.nto) ia=nuc(is-1)+1
      if(itype.and.is.le.nto) ia=nuc(is-1)+iofe(is)+1
      do js=is,ntl
      jn=itr(js)
      jno=ito(js)
      jcat=.false.
      if(nuc(js)-nuc(js-1).eq.1) jcat=.true.
      jtype=.false.
      if(jn.ne.jno.or.js.gt.nto) jtype=.true.
      ja=nuc(js-1)+iofs(js)
      if(js.gt.nto) ja=nuc(js-1)+1
      if(jtype.and.js.le.nto) ja=nuc(js-1)+iofe(js)+1
      if((in.ne.ino.or.is.gt.nto.or.jn.ne.jno.or.js.gt.nto)
     1 .and.start) then
      start=.false.
      if(.not.quiet) write(6,36)
36    format(/2x,'Cov/Lig interactions ...',16x,'Elj',5x,'Eel',5x,
     1 'Eang',4x,'Etog',3x,'Total')
      endif
c------------------------------------------------------------------i & j special
      if(itype.and.jtype) then
      enuc=0.
      do k=4,16,4
      enuc=enuc+ec(is,js,k,5)
      enddo
      if(abs(enuc).gt.emin) then
      if(is.gt.nto.and.(is.ne.js.or..not.icat.or..not.jcat)) then
      if(.not.quiet) write(6,37) is,munit(ia),nunit(ia),munit(ja),
     1 nunit(ja),(ec(is,js,16,l),l=1,5)
37    format(/2x,i2,') ',a4,1x,i3,7x,' & ',a4,1x,i3,' Clig ',5f8.2)
      else
      if(.not.quiet) write(6,35) is,munit(ia),nunit(ia),
     1 munit(ja),nunit(ja),((ec(is,js,k,l),l=1,5),k=4,16,4)
35    format(/2x,i2,') ',a4,1x,i3,7x,' & ',a4,1x,i3,' Phos ',5f8.2,
     1  /32x,' Sucr ',5f8.2,/32x,' Base ',5f8.2,:/32x,' Clig ',5f8.2)
      endif
      endif
c----------------------------------------------------------------------i special
      else if(itype) then
      enuc=0.
      do k=13,16
      enuc=enuc+ec(is,js,k,5)
      enddo
      if(abs(enuc).gt.emin.and..not.quiet) write(6,35) is,munit(ia),
     1 nunit(ia),munit(ja),nunit(ja),((ec(is,js,k,l),l=1,5),k=13,15)
c----------------------------------------------------------------------j special
      else if(jtype) then
      enuc=0.
      do k=4,16,4
      enuc=enuc+ec(is,js,k,5)
      enddo
      if(abs(enuc).gt.emin.and..not.quiet) write(6,35) js,munit(ja),
     1 nunit(ja),munit(ia),nunit(ia),((ec(is,js,k,l),l=1,5),k=4,12,4)
      endif
      enddo
      enddo
c-----------------------------------------------------------h-bonds and stacking
      if(kseq.gt.0) write(6,110)
110   format(/2x,'H-bonding and Stacking ....',
     1//4x,'Nuc',5x,'H-bond',3x,'Stack',8x,'1-1',3x,'1-2',3x,'2-1',
     1 3x,'2-2'/)
      do i=1,kseq
      ks=(ilq(i,2)-1)*kseq+ilq(i,1)
      base=seq(ks:ks)
      est=0
      ehb=0
      kpn=0
      kst=0
      do k1=1,nst
      il=ieq(i,k1)
      if(il.ne.0) then
      kpn=kpn+1
      do k2=1,nst
      if(i.lt.kseq) then
      iu=ieq(i+1,k2)
      if(iu.ne.0) then
      est=est+ec(min(il,iu),max(il,iu),11,5)
      kst=kst+1
      esta(kst)=ec(min(il,iu),max(il,iu),11,5)
      endif
      endif
      in=0
      if(k2.gt.k1) in=ieq(i,k2)
      if(in.ne.0) then
      ehb=ehb+ec(il,in,11,5)
      endif
      enddo
      endif
      enddo
      if(i.gt.isym) goto 190
      if(nst.ne.2.or.est.eq.0.) then
      write(6,120) i,base,ehb,est,kpn
      else
      write(6,120) i,base,ehb,est,kpn,(esta(j),j=1,4)
      endif
120   format(2x,i3,':',a2,2x,2f8.2,' (',i1,'b)',4f6.2)
190   enddo
c-----------------------------------------------------------h-bonds and contacts
      if(nhb.gt.0) then
      write(6,12)
12    format(/2x,'Hydrogen bond geometry ...'/)
      do l=1,nhb
      i=ihb(l,1)
      j=ihb(l,2)
      ii=nunit(i)
      if(ilq(ii,1).gt.isym) goto 195
      write(6,14) l,munit(i),nunit(i),mnam(i),munit(j),nunit(j),mnam(j),
     1 rhb(l,1),rhb(l,2)
14    format(2x,i2,') ',a4,i3,1x,a4,'  &  ',a4,i3,1x,a4,
     1 '  r= ',f6.2,'  ang= ',f6.1)
195   enddo
      endif
      if(nco.gt.0) then
      write(6,16)
16    format(/2x,'Close contacts ...'/)
      do l=1,nco
      i=ihc(l,1)
      j=ihc(l,2)
      ii=nunit(i)
      if(ilq(ii,1).gt.isym) goto 197
      write(6,18) l,munit(i),nunit(i),mnam(i),munit(j),nunit(j),mnam(j),
     1 rhc(l,1),rhc(l,2)
18    format(2x,i2,') ',a4,i3,1x,a4,'  &  ',a4,i3,1x,a4,
     1 '  r= ',f6.2,'  Elj= ',f10.1)
197   enddo
      endif
c-------------------------------------------------------------------------------
      if(.not.ecen) then
      do is=1,ntl
      repl=repl+en(is,is,0)
      disp=disp+en(is,is,1)-en(is,is,0)
      elec=elec+en(is,is,2)
      eang=eang+en(is,is,3)
      etog=etog+en(is,is,4)
      enddo
      endif
      ener=elec+disp+repl+eang+etog+epen
      if(ecen) call pairc(mp)
      return
      end
