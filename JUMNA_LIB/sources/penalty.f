      subroutine penalty
      include 'jumna_data.inc'
      logical*2 lthy,lar,lock,kink,ifhb,sup,rcom,homo,homo2,homo3,diep,
     1 sum,link,ecen,cyl,lcat,cent,autos,locr,ribose,cation,jsk,angle,
     1 hst,bst,vst,lgi,lgj,lgu,ihl,amber
      character*4 mnam,munit,seq*120,code*8,kode*8,lnam,knam,atclas*2
      integer*2 i23,i34,elim
      integer*4 opt,w1,w2,w3,wh,fold
      dimension u(3),p(3),g(3),c(3),t(3),s1(3),s2(3),fco(n2*2,3),
     1 dtc(7,3),dgc(7,3),dg(7),dgt(7,3),dk(7,3),dcc(7,3),
     1 dtt(7,3),uu(7,3,3),xuu(7,3,3),dgu(7,3),dpp(7,3,3),dc(7),
     1 v1(3),v2(3),q1(3),q2(3),dv1(4,3),dv2(4,3),dt(4,3,3),
     1 drv1(4,3),drv2(4,3),dq1(4,3,3),dq2(4,3,3)
      common/cur/fb1x,fb1y,fb1z,fb2x,fb2y,fb2z,fbx,fby,fbz,ib1,ib2
      common/datjm/acc,phos,epsi,epsr,plat,slope,rhbl,vfac,tfac,rfac,xfac,
     1 scale,rpiv,tpiv,fad,fan,damp,fnoew,fnoes,fnoea,df1,scnb,scee,
     1 catd,catr,catc,rad,pit,enit,nshel,maxn,opt,mhomo,mhomo2,mhomo3,
     1 nick,limit,nop,nrib,ncat,lig,naxo,sup,rcom,homo,homo2,homo3,diep,
     1 sum,link,ecen,cyl,lcat,cent,autos,amber
      common/dcu/dpx,dpy,dpz,dux,duy,duz,crx,cry,crz,ctx,cty,ctz,
     1 csx,csy,csz,crs,clx,cly,crl
      common/enf/for(n1,3),tor(n1,3),fot(n6),ener,elec,
     1 repl,disp,eang,etog,epen,i23(8*n1),i34(8*n1),elim(n1)
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/flx/sap(n6),refg,refb,refh,refv,iap(n8,n4,n5),
     1 nap(n8,7,n5),kap(3,n5)
      common/hel/ua(n2,3),da(n2,3),ra(n2,3)
      common/ind/iend(0:4),ise(n2),ienb(n2,2),kapt(n2+1),nsph
      common/lgd/rlig(n9,6),slig(n9,2),ilig(n9,6),lopt(n9),lpiv(n9),
     1 ntlg,nlgi,nlig,ntl,lnam(n9),locr(n2,n0)
      common/mnn/var(n7),gra(n7),con(n0*n2+n3),scl(n7),ncon,nvrc,
     1 ntba,nbac,nthe,nhel,ntki,nkin,ntri,ntot,nvar,nrin,
     1 lar(n7),lock(n6a)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/nmrjm/rnoe(n3),bnoe(n3),fnoe(n3),snoe,knam(200,2),ktyp(200),
     1 icon(n2*2),inoe(n3,9),ih68(n2),jnoe(n3),nnoe,ndcs,ndiv(0:99),nlin
      common/parjm/aij(25,25),bij(25,25),ahb(15),bhb(15),vt(35),va(25),
     1 vo(25),delmax,time0,nzsh,nwdg,ivf(35),ifhb(25,25)
      common/para/ro(57),reps(57),angc(191),fangc(191),btor(88),
     1 gtor(88),ftor(88),bdtor(30),gdtor(30),fdtor(30),amb(57,57),
     1 bmb(57,57),amh(9,11),bmh(9,11),mtor(88),itor(88,4),iadt(30,4),
     1 iangc(191,3),mof,nof,nang,ntor,ntors,nad,atclas(57)
      common/sst/phase(8,3),ampli(8,3),ph5(n2),am5(n2),kr5(n2)
      common/slg/hst(n2,6),bst(n2,n8),vst(n2,4),lgi(n1),lgj(n1),lgu(n2)
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(120,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      common/symjm/xsum,ysum,itwl(n2),nvs(n7),ihm(n2),ksym(3),nbrk(2),
     1 isym,isur,isup,ihl(n2)
      c72=cos(cdr*72.d0)
      s72=sin(cdr*72.d0)
      c144=cos(cdr*144.d0)
      s144=sin(cdr*144.d0)
      epen=0.
c      if(epen.eq.0) goto 600
c===============================================================ring constraints
      k=0
      do is=1,ntl
        do lr=1,nsr(is)
          if(locr(is,lr).or.(rcom.and.is.le.nto.and.lr.eq.1
     1                      .and.kr5(is).ne.0)) goto 10
          k=k+1
          i1=nuc(is-1)+isr(lr,1,is)
          i2=nuc(is-1)+isr(lr,2,is)
          dx=corm(i1,1)-corm(i2,1)
          dy=corm(i1,2)-corm(i2,2)
          dz=corm(i1,3)-corm(i2,3)
          r=sqrt(dx*dx+dy*dy+dz*dz)
          dist=r-rsr(lr,is)
          con(k)=dist
          epen=epen+fad*dist**2
          fcon=2*fad*dist/r
          for(i1,1)=for(i1,1)-dx*fcon
          for(i1,2)=for(i1,2)-dy*fcon
          for(i1,3)=for(i1,3)-dz*fcon
          for(i2,1)=for(i2,1)+dx*fcon
          for(i2,2)=for(i2,2)+dy*fcon
          for(i2,3)=for(i2,3)+dz*fcon
10      enddo
      enddo
c===========================================================junction constraints
      koff=nto-nst
      kk=nrin
      ll=kk+koff
      mm=ll+koff
      do ij=2,nto
      if(ise(ij).lt.0) goto 200
      is=ij
      ip=ij-1
      idir=idr(ilq(ij,2))
      if(idir.lt.0) then
      ip=ij
      is=ij-1
      endif
      iprib=0
      if(ribose(ip)) iprib=1
      isrib=0
      if(ribose(is)) isrib=1
      in=itr(is)
      ino=ito(is)
      isdel=kap(1,in)-kap(1,ino)+(nsr(is)-1)*5
      inp=itr(ip)
      inpo=ito(ip)
      ipdel=kap(1,inp)-kap(1,inpo)+(nsr(ip)-1)*5
      w1=kapt(ip)+16+ipdel+iprib
      w2=kapt(ip)+17+ipdel+iprib
      wh=16+isdel+isrib
      if(inpo.eq.2.or.inpo.eq.5) then
      w1=w1+3
      w2=w2+3
      endif
      if(ino.eq.3.or.ino.eq.6) wh=wh-4
      w3=kapt(is)+wh
      i1=nuc(ip-1)+7
      i2=nuc(ip-1)+15
      i3=nuc(ip-1)+16
      i4=nuc(is-1)+6
      i5=nuc(is-1)+4
      i6=nuc(is-1)+3
      i7=nuc(is-1)+13
c-------------------------------------------------------------------c5'-h motion
c     i1: O3' / i2: P / i3: O5' / i4: C5' / i5: C4' / i6: C3' / i7: H5'
c-------------------------------------------------------------------tors and ang
      x0=corm(i4,1)
      y0=corm(i4,2)
      z0=corm(i4,3)
      rx=x0-corm(i5,1)
      ry=y0-corm(i5,2)
      rz=z0-corm(i5,3)
      rr=sqrt(rx*rx+ry*ry+rz*rz)
      rx=rx/rr
      ry=ry/rr
      rz=rz/rr
         bx=corm(i3,1)-x0
         by=corm(i3,2)-y0
         bz=corm(i3,3)-z0
         px=by*rz-bz*ry
         py=bz*rx-bx*rz
         pz=bx*ry-by*rx
         pp=sqrt(px*px+py*py+pz*pz)
         px=px/pp
         py=py/pp
         pz=pz/pp
      ah=0.
      bh=0.
      do j=1,nuc(is)-nuc(is-1)
      m=iap(wh,j,in)
      if(m.ne.0) then
      jg=j+nuc(is-1)
      dx=corm(jg,1)-x0
      dy=corm(jg,2)-y0
      dz=corm(jg,3)-z0
      tx=dy*for(jg,3)-dz*for(jg,2)+tor(jg,1)
      ty=dz*for(jg,1)-dx*for(jg,3)+tor(jg,2)
      tz=dx*for(jg,2)-dy*for(jg,1)+tor(jg,3)
      ah=ah+(rx*tx+ry*ty+rz*tz)
      bh=bh+(px*tx+py*ty+pz*tz)/2
      endif
      enddo
      kk=kk+1
      ll=ll+1
      mm=mm+1
c---------------------------------------------------------------distance o5'-C5'
      dx=corm(i4,1)-corm(i3,1)
      dy=corm(i4,2)-corm(i3,2)
      dz=corm(i4,3)-corm(i3,3)
      r=sqrt(dx*dx+dy*dy+dz*dz)
      con(kk)=r-refb
      if(ij.eq.nick) goto 200
      epen=epen+fad*con(kk)**2
      fcon=2*fad*con(kk)/r
      for(i3,1)=for(i3,1)+dx*fcon
      for(i3,2)=for(i3,2)+dy*fcon
      for(i3,3)=for(i3,3)+dz*fcon
      for(i4,1)=for(i4,1)-dx*fcon
      for(i4,2)=for(i4,2)-dy*fcon
      for(i4,3)=for(i4,3)-dz*fcon
c-----------------------------------------------angles p-o5'-C5' and o5'-C5'-c4'
      if(.not.amber) then
      con(ll)=ang(i2,i3,i4)-vo(20)
      call delval(-fot(w3+2),i2,i3,i4)
      con(mm)=ang(i3,i4,i5)-vo(2)
      fota=fot(w3+1)+bh
      call delval(-fota,i3,i4,i5)
         else
         con(ll)=ang(i2,i3,i4)-angc(175)
         con(mm)=ang(i3,i4,i5)-angc(74)
         call delval(-bh,i3,i4,i5)
         endif
c--------------------------------------------torsions p-o5', O5'-c5' AND C5'-c4'
      if(.not.amber) then
      call deltor(sap(w1),fot(w1),i1,i2,i3,i4)
      call deltor(sap(w2),fot(w2),i2,i3,i4,i5)
      endif
      fota=fot(w3)+ah
      call deltor(sap(w3),fota,i3,i4,i5,i6)
      sah=torp(i7,i4,i5,i6)
      call deltor(sah,-ah,i7,i4,i5,i6)
200   enddo
c===============================================================superhelical fix
600   if(cent) then
      klim=ksym(1)
      if(klim.gt.kseq) klim=kseq
      xsum=0.
      ysum=0.
      wdg=0.
      do is=1,klim
      x=hel(is,1)
      y=hel(is,2)
      wdg=wdg+hel(is,6)*cdr
      if(isur.ne.0) wdg=hel(is,6)*cdr
      cwd=cos(wdg)
      swd=sin(wdg)
      xsum=xsum+x*cwd+y*swd
      ysum=ysum+x*swd-y*cwd
      enddo
      epen=epen+fnoes*(xsum**2+ysum**2)
      endif
c=======================================================================dist con
      do k=1,nnoe
      jsk=.false.
      jnk=jnoe(k)
      if(jnk.lt.0) then
      jnk=-jnk
      jsk=.true.
      endif
      i=ndcs+k
      if(jnk.eq.1.or.jnk.eq.3.or.jnk.eq.11) then
      i1=inoe(k,1)
      i2=inoe(k,2)
      xx=corm(i2,1)-corm(i1,1)
      yy=corm(i2,2)-corm(i1,2)
      zz=corm(i2,3)-corm(i1,3)
      r=sqrt(xx*xx+yy*yy+zz*zz)
      if(jsk) rnoe(k)=r
      if(r.lt.rnoe(k).or.jnk.eq.1.or.jnk.eq.11) then
      dist=r-rnoe(k)
      else if(r.gt.bnoe(k)) then
      dist=r-bnoe(k)
      else
      goto 250
      endif
      con(i)=dist
      epen=epen+fnoe(k)*dist**2
      fpen=2*fnoe(k)*dist/r
      for(i1,1)=for(i1,1)+xx*fpen
      for(i1,2)=for(i1,2)+yy*fpen
      for(i1,3)=for(i1,3)+zz*fpen
      for(i2,1)=for(i2,1)-xx*fpen
      for(i2,2)=for(i2,2)-yy*fpen
      for(i2,3)=for(i2,3)-zz*fpen
c================================================================double dist con
      else if(jnk.eq.26.or.jnk.eq.28) then
      i1=inoe(k,1)
      i2=inoe(k,2)
      xx1=corm(i2,1)-corm(i1,1)
      yy1=corm(i2,2)-corm(i1,2)
      zz1=corm(i2,3)-corm(i1,3)
      r1=sqrt(xx1*xx1+yy1*yy1+zz1*zz1)
      i3=inoe(k,3)
      i4=inoe(k,4)
      xx2=corm(i4,1)-corm(i3,1)
      yy2=corm(i4,2)-corm(i3,2)
      zz2=corm(i4,3)-corm(i3,3)
      r2=sqrt(xx2*xx2+yy2*yy2+zz2*zz2)
      r=r1+inoe(k,5)*r2
      if(jsk) rnoe(k)=r
      if(r.lt.rnoe(k).or.jnk.eq.26) then
      dist=r-rnoe(k)
      else if(r.gt.bnoe(k)) then
      dist=r-bnoe(k)
      else
      goto 250
      endif
      con(i)=dist
      epen=epen+fnoe(k)*dist**2
      fpen1=2*fnoe(k)*dist/r1
      fpen2=inoe(k,5)*2*fnoe(k)*dist/r2
      for(i1,1)=for(i1,1)+xx1*fpen1
      for(i1,2)=for(i1,2)+yy1*fpen1
      for(i1,3)=for(i1,3)+zz1*fpen1
      for(i2,1)=for(i2,1)-xx1*fpen1
      for(i2,2)=for(i2,2)-yy1*fpen1
      for(i2,3)=for(i2,3)-zz1*fpen1
      for(i3,1)=for(i3,1)+xx2*fpen2
      for(i3,2)=for(i3,2)+yy2*fpen2
      for(i3,3)=for(i3,3)+zz2*fpen2
      for(i4,1)=for(i4,1)-xx2*fpen2
      for(i4,2)=for(i4,2)-yy2*fpen2
      for(i4,3)=for(i4,3)-zz2*fpen2
c====================================================================torsion con
      else if(jnk.eq.2.or.jnk.eq.4) then
      i1=inoe(k,1)
      i2=inoe(k,2)
      i3=inoe(k,3)
      i4=inoe(k,4)
      dih=torp(i1,i2,i3,i4)
      if(jsk) rnoe(k)=dih
      delt=dih-rnoe(k)
      if(jnk.eq.2) goto 201
      cx=cos(cdr*(bnoe(k)))
      sx=sin(cdr*(bnoe(k)))
      cm=cos(cdr*(rnoe(k)))
      sm=sin(cdr*(rnoe(k)))
      ct=cos(cdr*(dih))
      st=sin(cdr*(dih))
      v=cx*sm-sx*cm
      a=cm*st-sm*ct
      b=ct*sx-st*cx
      if((a.ge.0.and.b.ge.0).or.(v.gt.0.and.a*b.le.0)) goto 250
      if(a.gt.0.and.b.lt.0) then
      delt=dih-bnoe(k)
      else if(a.lt.0.and.b.lt.0) then
      a=ct*cm+st*sm
      b=cx*ct+sx*st
      if(a.lt.b) then
      delt=dih-bnoe(k)
      endif
      endif
201   if(abs(delt).gt.180.) delt=delt-sign(360.d0,delt)
      con(i)=delt
      delt=delt*cdr
      epen=epen+fnoe(k)*delt**2
      fpen=-2*fnoe(k)*delt
      call deltor(dih,fpen,i1,i2,i3,i4)
c=============================================================double torsion con
      else if(jnk.eq.14.or.jnk.eq.16) then
      i1=inoe(k,1)
      i2=inoe(k,2)
      i3=inoe(k,3)
      i4=inoe(k,4)
      dih1=torp(i1,i2,i3,i4)
      i5=inoe(k,5)
      i6=inoe(k,6)
      i7=inoe(k,7)
      i8=inoe(k,8)
      i9=inoe(k,9)
      dih2=torp(i5,i6,i7,i8)
      dih=dih1+sign(1,i9)*dih2
      if(jsk) rnoe(k)=dih
      delt=dih-rnoe(k)
      if(jnk.eq.14) goto 211
      cx=cos(cdr*(bnoe(k)))
      sx=sin(cdr*(bnoe(k)))
      cm=cos(cdr*(rnoe(k)))
      sm=sin(cdr*(rnoe(k)))
      ct=cos(cdr*(dih))
      st=sin(cdr*(dih))
      v=cx*sm-sx*cm
      a=cm*st-sm*ct
      b=ct*sx-st*cx
      if((a.ge.0.and.b.ge.0).or.(v.gt.0.and.a*b.le.0)) goto 250
      if(a.gt.0.and.b.lt.0) then
      delt=dih-bnoe(k)
      else if(a.lt.0.and.b.lt.0) then
      a=ct*cm+st*sm
      b=cx*ct+sx*st
      if(a.lt.b) then
      delt=dih-bnoe(k)
      endif
      endif
211   if(abs(delt).gt.180.) delt=delt-sign(360.d0,delt)
      con(i)=delt
      delt=delt*cdr
      epen=epen+fnoe(k)*delt**2
      fpen=-2*fnoe(k)*delt
      call deltor(dih1,fpen,i1,i2,i3,i4)
      call deltor(sign(1,i9)*dih2,fpen,i5,i6,i7,i8)
c================================================================sugar phase con
      else if(jnk.eq.5.or.jnk.eq.6) then
      koff=inoe(k,1)-1
      i2=nuc(koff)+1
      i3=i2+1
      i4=i3+1
      i5=i4+1
      i1=i5+1
      t1=torp(i1,i2,i3,i4)
      t2=torp(i2,i3,i4,i5)
      t3=torp(i3,i4,i5,i1)
      t4=torp(i4,i5,i1,i2)
      t5=torp(i5,i1,i2,i3)
      a=0.4*(t1*c144+t2+t3*c144+t4*c72+t5*c72)
      b=0.4*(t1*s144   -t3*s144+t4*s72-t5*s72)
      amp=sqrt(a**2+b**2)
      pha=acos(a/amp)*crd
      if(b.lt.0) pha=360.-pha
      fac=-crd/sqrt(1.-(a/amp)**2)
      if(b.lt.0) fac=-fac
      faa=-0.16*a/amp**3
      if(jsk) rnoe(k)=pha
      delt=pha-rnoe(k)
      if(jnk.eq.5) goto 202
      cx=cos(cdr*(bnoe(k)))
      sx=sin(cdr*(bnoe(k)))
      cm=cos(cdr*(rnoe(k)))
      sm=sin(cdr*(rnoe(k)))
      ct=cos(cdr*(pha))
      st=sin(cdr*(pha))
      v=cx*sm-sx*cm
      a=cm*st-sm*ct
      b=ct*sx-st*cx
      if((a.ge.0.and.b.ge.0).or.(v.gt.0.and.a*b.le.0)) goto 250
      if(a.gt.0.and.b.lt.0) then
      delt=pha-bnoe(k)
      else if(a.lt.0.and.b.lt.0) then
      a=ct*cm+st*sm
      b=cx*ct+sx*st
      if(a.lt.b) then
      delt=pha-bnoe(k)
      endif
      endif
202   if(abs(delt).gt.180.) delt=delt-sign(360.d0,delt)
      con(i)=delt
      delt=delt*cdr
      epen=epen+fnoe(k)*delt**2
      fpen=-2*fnoe(k)*delt
      dt1=fac*(faa*(t1+c72*(t3+t4)+c144*(t2+t5))+0.4*c144/amp)
      dt2=fac*(faa*(t2+c72*(t4+t5)+c144*(t1+t3))+0.4     /amp)
      dt3=fac*(faa*(t3+c72*(t1+t5)+c144*(t2+t4))+0.4*c144/amp)
      dt4=fac*(faa*(t4+c72*(t1+t2)+c144*(t3+t5))+0.4*c72 /amp)
      dt5=fac*(faa*(t5+c72*(t2+t3)+c144*(t1+t4))+0.4*c72 /amp)
      call deltor(t1,fpen*dt1,i1,i2,i3,i4)
      call deltor(t2,fpen*dt2,i2,i3,i4,i5)
      call deltor(t3,fpen*dt3,i3,i4,i5,i1)
      call deltor(t4,fpen*dt4,i4,i5,i1,i2)
      call deltor(t5,fpen*dt5,i5,i1,i2,i3)
c=============================================================double sugar phase
      else if(jnk.eq.18.or.jnk.eq.19) then
      koff=inoe(k,1)-1
      i2=nuc(koff)+1
      i3=i2+1
      i4=i3+1
      i5=i4+1
      i1=i5+1
      t1=torp(i1,i2,i3,i4)
      t2=torp(i2,i3,i4,i5)
      t3=torp(i3,i4,i5,i1)
      t4=torp(i4,i5,i1,i2)
      t5=torp(i5,i1,i2,i3)
      a=0.4*(t1*c144+t2+t3*c144+t4*c72+t5*c72)
      b=0.4*(t1*s144   -t3*s144+t4*s72-t5*s72)
      amp1=sqrt(a**2+b**2)
      pha1=acos(a/amp1)*crd
      if(b.lt.0) pha1=360.-pha1
      fac1=-crd/sqrt(1.-(a/amp1)**2)
      if(b.lt.0) fac1=-fac1
      faa1=-0.16*a/amp1**3
         koff=inoe(k,2)-1
         j2=nuc(koff)+1
         j3=j2+1
         j4=j3+1
         j5=j4+1
         j1=j5+1
         p1=torp(j1,j2,j3,j4)
         p2=torp(j2,j3,j4,j5)
         p3=torp(j3,j4,j5,j1)
         p4=torp(j4,j5,j1,j2)
         p5=torp(j5,j1,j2,j3)
         a=0.4*(p1*c144+p2+p3*c144+p4*c72+p5*c72)
         b=0.4*(p1*s144   -p3*s144+p4*s72-p5*s72)
         amp2=sqrt(a**2+b**2)
         pha2=acos(a/amp2)*crd
         if(b.lt.0) pha2=360.-pha2
         fac2=-crd/sqrt(1.-(a/amp2)**2)
         if(b.lt.0) fac2=-fac2
         faa2=-0.16*a/amp2**3
      i9=inoe(k,9)
      dpha=pha1+i9*pha2
      if(jsk) rnoe(k)=dpha
      if(dpha.lt.rnoe(k).or.jnk.eq.18) then
      delt=dpha-rnoe(k)
      else if(dpha.gt.bnoe(k)) then
      delt=dpha-bnoe(k)
      else
      delt=0.
      con(i)=0.
      goto 250
      endif
      con(i)=delt
      delt=delt*cdr
      epen=epen+fnoe(k)*delt**2
      fpen=-2*fnoe(k)*delt
      dt1=fac1*(faa1*(t1+c72*(t3+t4)+c144*(t2+t5))+0.4*c144/amp1)
      dt2=fac1*(faa1*(t2+c72*(t4+t5)+c144*(t1+t3))+0.4     /amp1)
      dt3=fac1*(faa1*(t3+c72*(t1+t5)+c144*(t2+t4))+0.4*c144/amp1)
      dt4=fac1*(faa1*(t4+c72*(t1+t2)+c144*(t3+t5))+0.4*c72 /amp1)
      dt5=fac1*(faa1*(t5+c72*(t2+t3)+c144*(t1+t4))+0.4*c72 /amp1)
         dp1=fac2*(faa2*(p1+c72*(p3+p4)+c144*(p2+p5))+0.4*c144/amp2)
         dp2=fac2*(faa2*(p2+c72*(p4+p5)+c144*(p1+p3))+0.4     /amp2)
         dp3=fac2*(faa2*(p3+c72*(p1+p5)+c144*(p2+p4))+0.4*c144/amp2)
         dp4=fac2*(faa2*(p4+c72*(p1+p2)+c144*(p3+p5))+0.4*c72 /amp2)
         dp5=fac2*(faa2*(p5+c72*(p2+p3)+c144*(p1+p4))+0.4*c72 /amp2)
      call deltor(t1,fpen*dt1,i1,i2,i3,i4)
      call deltor(t2,fpen*dt2,i2,i3,i4,i5)
      call deltor(t3,fpen*dt3,i3,i4,i5,i1)
      call deltor(t4,fpen*dt4,i4,i5,i1,i2)
      call deltor(t5,fpen*dt5,i5,i1,i2,i3)
         call deltor(i9*p1,fpen*dp1,j1,j2,j3,j4)
         call deltor(i9*p2,fpen*dp2,j2,j3,j4,j5)
         call deltor(i9*p3,fpen*dp3,j3,j4,j5,j1)
         call deltor(i9*p4,fpen*dp4,j4,j5,j1,j2)
         call deltor(i9*p5,fpen*dp5,j5,j1,j2,j3)
c============================================================sugar amplitude con
      else if(jnk.eq.7.or.jnk.eq.8) then
      koff=inoe(k,1)-1
      i2=nuc(koff)+1
      i3=i2+1
      i4=i3+1
      i5=i4+1
      i1=i5+1
      t1=torp(i1,i2,i3,i4)
      t2=torp(i2,i3,i4,i5)
      t3=torp(i3,i4,i5,i1)
      t4=torp(i4,i5,i1,i2)
      t5=torp(i5,i1,i2,i3)
      a=0.4*(t1*c144+t2+t3*c144+t4*c72+t5*c72)
      b=0.4*(t1*s144   -t3*s144+t4*s72-t5*s72)
      amp=sqrt(a**2+b**2)
      if(jsk) rnoe(k)=amp
      if(amp.lt.rnoe(k).or.jnk.eq.7) then
      delt=amp-rnoe(k)
      else if(amp.gt.bnoe(k)) then
      delt=amp-bnoe(k)
      else
      delt=0.
      con(i)=0.
      goto 250
      endif
      con(i)=delt
      delt=delt*cdr
      epen=epen+fnoe(k)*delt**2
      fpen=-2*fnoe(k)*delt
      faa=0.16/amp
      dt1=faa*(t1+c72*(t3+t4)+c144*(t2+t5))
      dt2=faa*(t2+c72*(t4+t5)+c144*(t1+t3))
      dt3=faa*(t3+c72*(t1+t5)+c144*(t2+t4))
      dt4=faa*(t4+c72*(t1+t2)+c144*(t3+t5))
      dt5=faa*(t5+c72*(t2+t3)+c144*(t1+t4))
      call deltor(t1,fpen*dt1,i1,i2,i3,i4)
      call deltor(t2,fpen*dt2,i2,i3,i4,i5)
      call deltor(t3,fpen*dt3,i3,i4,i5,i1)
      call deltor(t4,fpen*dt4,i4,i5,i1,i2)
      call deltor(t5,fpen*dt5,i5,i1,i2,i3)
c============================================================sugar amplitude con
      else if(jnk.eq.24.or.jnk.eq.25) then
      koff=inoe(k,1)-1
      i2=nuc(koff)+1
      i3=i2+1
      i4=i3+1
      i5=i4+1
      i1=i5+1
      t1=torp(i1,i2,i3,i4)
      t2=torp(i2,i3,i4,i5)
      t3=torp(i3,i4,i5,i1)
      t4=torp(i4,i5,i1,i2)
      t5=torp(i5,i1,i2,i3)
      a=0.4*(t1*c144+t2+t3*c144+t4*c72+t5*c72)
      b=0.4*(t1*s144   -t3*s144+t4*s72-t5*s72)
      amp1=sqrt(a**2+b**2)
      faa1=0.16/amp1
         koff=inoe(k,2)-1
         j2=nuc(koff)+1
         j3=j2+1
         j4=j3+1
         j5=j4+1
         j1=j5+1
         p1=torp(j1,j2,j3,j4)
         p2=torp(j2,j3,j4,j5)
         p3=torp(j3,j4,j5,j1)
         p4=torp(j4,j5,j1,j2)
         p5=torp(j5,j1,j2,j3)
         a=0.4*(p1*c144+p2+p3*c144+p4*c72+p5*c72)
         b=0.4*(p1*s144   -p3*s144+p4*s72-p5*s72)
         amp2=sqrt(a**2+b**2)
         faa2=0.16/amp2
      i9=inoe(k,9)
      famp=amp1+i9*amp2
      if(jsk) rnoe(k)=famp
      if(famp.lt.rnoe(k).or.jnk.eq.24) then
      delt=famp-rnoe(k)
      else if(famp.gt.bnoe(k)) then
      delt=famp-bnoe(k)
      else
      delt=0.
      con(i)=0.
      goto 250
      endif
      con(i)=delt
      delt=delt*cdr
      epen=epen+fnoe(k)*delt**2
      fpen=-2*fnoe(k)*delt
      dt1=faa1*(t1+c72*(t3+t4)+c144*(t2+t5))
      dt2=faa1*(t2+c72*(t4+t5)+c144*(t1+t3))
      dt3=faa1*(t3+c72*(t1+t5)+c144*(t2+t4))
      dt4=faa1*(t4+c72*(t1+t2)+c144*(t3+t5))
      dt5=faa1*(t5+c72*(t2+t3)+c144*(t1+t4))
         dp1=faa2*(p1+c72*(p3+p4)+c144*(p2+p5))
         dp2=faa2*(p2+c72*(p4+p5)+c144*(p1+p3))
         dp3=faa2*(p3+c72*(p1+p5)+c144*(p2+p4))
         dp4=faa2*(p4+c72*(p1+p2)+c144*(p3+p5))
         dp5=faa2*(p5+c72*(p2+p3)+c144*(p1+p4))
      call deltor(t1,fpen*dt1,i1,i2,i3,i4)
      call deltor(t2,fpen*dt2,i2,i3,i4,i5)
      call deltor(t3,fpen*dt3,i3,i4,i5,i1)
      call deltor(t4,fpen*dt4,i4,i5,i1,i2)
      call deltor(t5,fpen*dt5,i5,i1,i2,i3)
         call deltor(i9*p1,fpen*dp1,j1,j2,j3,j4)
         call deltor(i9*p2,fpen*dp2,j2,j3,j4,j5)
         call deltor(i9*p3,fpen*dp3,j3,j4,j5,j1)
         call deltor(i9*p4,fpen*dp4,j4,j5,j1,j2)
         call deltor(i9*p5,fpen*dp5,j5,j1,j2,j3)
c===============================================================sugar pucker con
      else if(jnk.eq.21.or.jnk.eq.22) then
        koff=inoe(k,1)-1
        i2=nuc(koff)+1
        i3=i2+1
        i4=i3+1
        i5=i4+1
        i1=i5+1
        ip=inoe(k,2)
        if(jnk.eq.21) then
          rnoe(k)=phase(ip,3)
        else
          rnoe(k)=phase(ip,1)
          bnoe(k)=phase(ip,2)
        endif
        t1=torp(i1,i2,i3,i4)
        t2=torp(i2,i3,i4,i5)
        t3=torp(i3,i4,i5,i1)
        t4=torp(i4,i5,i1,i2)
        t5=torp(i5,i1,i2,i3)
        a=0.4*(t1*c144+t2+t3*c144+t4*c72+t5*c72)
        b=0.4*(t1*s144   -t3*s144+t4*s72-t5*s72)
        amp=sqrt(a**2+b**2)
        pha=acos(a/amp)*crd
        if(b.lt.0) pha=360.-pha
          fac=-crd/sqrt(1.-(a/amp)**2)
        if(b.lt.0) fac=-fac
          faa=-0.16*a/amp**3
          delt=pha-rnoe(k)
        if(jnk.eq.21) goto 203
        cx=cos(cdr*(bnoe(k)))
        sx=sin(cdr*(bnoe(k)))
        cm=cos(cdr*(rnoe(k)))
        sm=sin(cdr*(rnoe(k)))
        ct=cos(cdr*(pha))
        st=sin(cdr*(pha))
        v=cx*sm-sx*cm
        a=cm*st-sm*ct
        b=ct*sx-st*cx
        if((a.ge.0.and.b.ge.0).or.(v.gt.0.and.a*b.le.0)) goto 204
        if(a.gt.0.and.b.lt.0) then
          delt=pha-bnoe(k)
        else if(a.lt.0.and.b.lt.0) then
          a=ct*cm+st*sm
          b=cx*ct+sx*st
          if(a.lt.b) then
            delt=pha-bnoe(k)
          endif
        endif
203     if(abs(delt).gt.180.) delt=delt-sign(360.d0,delt)
        con(i)=delt
        delt=delt*cdr
        epen=epen+fnoe(k)*delt**2
        fpen=-2*fnoe(k)*delt
        dt1=fac*(faa*(t1+c72*(t3+t4)+c144*(t2+t5))+0.4*c144/amp)
        dt2=fac*(faa*(t2+c72*(t4+t5)+c144*(t1+t3))+0.4     /amp)
        dt3=fac*(faa*(t3+c72*(t1+t5)+c144*(t2+t4))+0.4*c144/amp)
        dt4=fac*(faa*(t4+c72*(t1+t2)+c144*(t3+t5))+0.4*c72 /amp)
        dt5=fac*(faa*(t5+c72*(t2+t3)+c144*(t1+t4))+0.4*c72 /amp)
        call deltor(t1,fpen*dt1,i1,i2,i3,i4)
        call deltor(t2,fpen*dt2,i2,i3,i4,i5)
        call deltor(t3,fpen*dt3,i3,i4,i5,i1)
        call deltor(t4,fpen*dt4,i4,i5,i1,i2)
        call deltor(t5,fpen*dt5,i5,i1,i2,i3)
204     if(jnk.eq.21) then
          rnoe(k)=ampli(ip,3)
        else
          rnoe(k)=ampli(ip,1)
          bnoe(k)=ampli(ip,2)
        endif
         if(amp.lt.rnoe(k).or.jnk.eq.21) then
         delt=amp-rnoe(k)
         else if(amp.gt.bnoe(k)) then
         delt=amp-bnoe(k)
         else
         goto 250
         endif
         delt=delt*cdr
         epen=epen+fnoe(k)*delt**2
         fpen=-2*fnoe(k)*delt
         faa=0.16/amp
         dt1=faa*(t1+c72*(t3+t4)+c144*(t2+t5))
         dt2=faa*(t2+c72*(t4+t5)+c144*(t1+t3))
         dt3=faa*(t3+c72*(t1+t5)+c144*(t2+t4))
         dt4=faa*(t4+c72*(t1+t2)+c144*(t3+t5))
         dt5=faa*(t5+c72*(t2+t3)+c144*(t1+t4))
         call deltor(t1,fpen*dt1,i1,i2,i3,i4)
         call deltor(t2,fpen*dt2,i2,i3,i4,i5)
         call deltor(t3,fpen*dt3,i3,i4,i5,i1)
         call deltor(t4,fpen*dt4,i4,i5,i1,i2)
         call deltor(t5,fpen*dt5,i5,i1,i2,i3)
c==============================================================valence angle con
      else if(jnk.eq.9) then
      i1=inoe(k,1)
      i2=inoe(k,2)
      i3=inoe(k,3)
      val=ang(i1,i2,i3)
      if(jsk) rnoe(k)=val
      delt=val-rnoe(k)
      if(abs(delt).gt.180.) delt=delt-sign(360.d0,delt)
      con(i)=delt
      delt=delt*cdr
      epen=epen+fnoe(k)*delt**2
      fpen=2*fnoe(k)*delt
      call delval(fpen,i1,i2,i3)
c========================================================================opening
      else if(jnk.eq.29) then
      i1=inoe(k,1)
      i2=inoe(k,2)
      i3=inoe(k,3)
      i4=inoe(k,4)
      i5=inoe(k,5)
      i6=inoe(k,6)
      i7=inoe(k,7)
      do i=1,7
      dc(i)=0.
      dg(i)=0.
         do j=1,3
         dcc(i,j)=0.
         dgc(i,j)=0.
         dtt(i,j)=0.
            do m=1,3
            uu(i,j,m)=0.
            xuu(i,j,m)=0.
            enddo
         enddo
      enddo
      dc(1)=-1.
      dc(3)= 1.
      dg(1)=-1.
      dg(2)= 1.
      do j=1,3
      g(j)=corm(i2,j)-corm(i1,j)
      c(j)=corm(i3,j)-corm(i1,j)
      s1(j)=corm(i5,j)-corm(i4,j)
      s2(j)=corm(i7,j)-corm(i6,j)
      enddo
      rc2=c(1)*c(1)+c(2)*c(2)+c(3)*c(3)
      rc=sqrt(rc2)
      rc3=rc*rc2
      rc4=rc2*rc2
      dotgc=g(1)*c(1)+g(2)*c(2)+g(3)*c(3)
      rs1=sqrt(s1(1)*s1(1)+s1(2)*s1(2)+s1(3)*s1(3))
      rs2=sqrt(s2(1)*s2(1)+s2(2)*s2(2)+s2(3)*s2(3))
      t(1)=s1(1)/rs1+s2(1)/rs2
      t(2)=s1(2)/rs1+s2(2)/rs2
      t(3)=s1(3)/rs1+s2(3)/rs2
      dotc=t(1)*c(1)+t(2)*c(2)+t(3)*c(3)
      u(1)=t(1)-dotc*c(1)/rc2
      u(2)=t(2)-dotc*c(2)/rc2
      u(3)=t(3)-dotc*c(3)/rc2
      ru2=u(1)*u(1)+u(2)*u(2)+u(3)*u(3)
      ru=sqrt(ru2)
      dug=g(1)*u(1)+g(2)*u(2)+g(3)*u(3)
      p(1)=g(1)-dug*u(1)/ru2
      p(2)=g(2)-dug*u(2)/ru2
      p(3)=g(3)-dug*u(3)/ru2
      rp=sqrt(p(1)**2+p(2)**2+p(3)**2)
      dpu= p(1)*u(1)+p(2)*u(2)+p(3)*u(3)
      dpc=(p(1)*c(1)+p(2)*c(2)+p(3)*c(3))/(rp*rc)
      val=acos(dpc)*crd
         dx=c(2)*p(3)-c(3)*p(2)
         dy=c(3)*p(1)-c(1)*p(3)
         dz=c(1)*p(2)-c(2)*p(1)
         dot=(dx*u(1)+dy*u(2)+dz*u(3))*inoe(k,8)
         if(dot.lt.0) val=-val
      if(jsk) rnoe(k)=val
      delt=val-rnoe(k)
      con(i)=delt
      epen=epen+fnoe(k)*delt**2
      fpen=2*fnoe(k)*delt*sign(1.d0,dot)*crd/sqrt(1-dpc*dpc)
      do j=1,3
      dgc(1,j)=-g(j)-c(j)
      dgc(2,j)= c(j)
      dgc(3,j)= g(j)
      dgc(4,j)= 0.
      enddo
      do j=1,3
      dcc(1,j)=-c(j)
      dcc(3,j)= c(j)
      enddo
      r3s1=rs1**3
      r3s2=rs2**3
      do j=1,3
      dtt(4,j)=-1/rs1+1/r3s1*s1(j)*s1(j)
      dtt(5,j)= 1/rs1-1/r3s1*s1(j)*s1(j)
      dtt(6,j)=-1/rs2+1/r3s2*s2(j)*s2(j)
      dtt(7,j)= 1/rs2-1/r3s2*s2(j)*s2(j)
      enddo
c-----------------------------------d dotgt
      do j=1,3
      dgt(1,j)=-t(j)
      dgt(2,j)= t(j)
      dgt(3,j)= 0.
      do i=4,7
      dgt(i,j)=g(j)*dtt(i,j)
      enddo
      enddo
      dgt(4,1)=dgt(4,1)+(g(2)*s1(2)+g(3)*s1(3))*s1(1)/r3s1
      dgt(4,2)=dgt(4,2)+(g(1)*s1(1)+g(3)*s1(3))*s1(2)/r3s1
      dgt(4,3)=dgt(4,3)+(g(1)*s1(1)+g(2)*s1(2))*s1(3)/r3s1
      dgt(5,1)=dgt(5,1)-(g(2)*s1(2)+g(3)*s1(3))*s1(1)/r3s1
      dgt(5,2)=dgt(5,2)-(g(1)*s1(1)+g(3)*s1(3))*s1(2)/r3s1
      dgt(5,3)=dgt(5,3)-(g(1)*s1(1)+g(2)*s1(2))*s1(3)/r3s1
      dgt(6,1)=dgt(6,1)+(g(2)*s2(2)+g(3)*s2(3))*s2(1)/r3s2
      dgt(6,2)=dgt(6,2)+(g(1)*s2(1)+g(3)*s2(3))*s2(2)/r3s2
      dgt(6,3)=dgt(6,3)+(g(1)*s2(1)+g(2)*s2(2))*s2(3)/r3s2
      dgt(7,1)=dgt(7,1)-(g(2)*s2(2)+g(3)*s2(3))*s2(1)/r3s2
      dgt(7,2)=dgt(7,2)-(g(1)*s2(1)+g(3)*s2(3))*s2(2)/r3s2
      dgt(7,3)=dgt(7,3)-(g(1)*s2(1)+g(2)*s2(2))*s2(3)/r3s2
c-----------------------------------d dotct
      do j=1,3
      dtc(1,j)=-t(j)/rc+dotc/rc3*c(j)
      dtc(3,j)= t(j)/rc-dotc/rc3*c(j)
      dtc(2,j)=0.
      do i=4,7
      dtc(i,j)=c(j)*dtt(i,j)/rc
      enddo
      enddo
      dtc(4,1)=dtc(4,1)+(c(2)*s1(2)+c(3)*s1(3))*s1(1)/r3s1/rc
      dtc(4,2)=dtc(4,2)+(c(1)*s1(1)+c(3)*s1(3))*s1(2)/r3s1/rc
      dtc(4,3)=dtc(4,3)+(c(1)*s1(1)+c(2)*s1(2))*s1(3)/r3s1/rc
      dtc(5,1)=dtc(5,1)-(c(2)*s1(2)+c(3)*s1(3))*s1(1)/r3s1/rc
      dtc(5,2)=dtc(5,2)-(c(1)*s1(1)+c(3)*s1(3))*s1(2)/r3s1/rc
      dtc(5,3)=dtc(5,3)-(c(1)*s1(1)+c(2)*s1(2))*s1(3)/r3s1/rc
      dtc(6,1)=dtc(6,1)+(c(2)*s2(2)+c(3)*s2(3))*s2(1)/r3s2/rc
      dtc(6,2)=dtc(6,2)+(c(1)*s2(1)+c(3)*s2(3))*s2(2)/r3s2/rc
      dtc(6,3)=dtc(6,3)+(c(1)*s2(1)+c(2)*s2(2))*s2(3)/r3s2/rc
      dtc(7,1)=dtc(7,1)-(c(2)*s2(2)+c(3)*s2(3))*s2(1)/r3s2/rc
      dtc(7,2)=dtc(7,2)-(c(1)*s2(1)+c(3)*s2(3))*s2(2)/r3s2/rc
      dtc(7,3)=dtc(7,3)-(c(1)*s2(1)+c(2)*s2(2))*s2(3)/r3s2/rc
c------------------------------------------
      do j=1,3
      do m=1,3
      uu(4,j,m)= s1(j)*s1(m)/r3s1
      uu(5,j,m)=-s1(j)*s1(m)/r3s1
      uu(6,j,m)= s2(j)*s2(m)/r3s2
      uu(7,j,m)=-s2(j)*s2(m)/r3s2
      enddo
      uu(4,j,j)=uu(4,j,j)-1/rs1
      uu(5,j,j)=uu(5,j,j)+1/rs1
      uu(6,j,j)=uu(6,j,j)-1/rs2
      uu(7,j,j)=uu(7,j,j)+1/rs2
      enddo
c-----------------------------------------
      do i=1,7
      do j=1,3
      do m=1,3
      uu(i,j,m)=uu(i,j,m)-dtc(i,j)*c(m)/rc
     1 +dotc*c(j)*dcc(i,m)/rc4
      enddo
      uu(i,j,j)=uu(i,j,j)-dotc*dc(i)/rc2
      enddo
      enddo
c----------------------------------- u/|u|
      do i=1,7
      do j=1,3
      dk(i,j)=0.
      do m=1,3
      dk(i,j)=dk(i,j)+uu(i,j,m)*u(m)
      enddo
      do m=1,3
      xuu(i,j,m)=uu(i,j,m)/ru-1./(ru2*ru)*dk(i,j)*u(m)
      enddo
      enddo
      enddo
c------------------------------------- dpp
      do i=1,7
      do j=1,3
      dgu(i,j)=0.
      do m=1,3
      dgu(i,j)=dgu(i,j)+g(m)*xuu(i,j,m)
      enddo
      dgu(i,j)=dg(i)*u(j)/ru+dgu(i,j)
      do m=1,3
      dpp(i,j,m)=-dgu(i,j)*u(m)/ru-dug*xuu(i,j,m)/ru
      enddo
      dpp(i,j,j)=dpp(i,j,j)+dg(i)
      enddo
      enddo
c-----------------------------------------
      do i=1,7
      ia=inoe(k,i)
      do j=1,3
      dot=dpp(i,j,1)*p(1)+dpp(i,j,2)*p(2)+dpp(i,j,3)*p(3)
      for(ia,j)=for(ia,j)+fpen*(dgc(i,j)/(rp*rc)
     1 -dotgc*(dot/(rp**3*rc)+dcc(i,j)/(rc**3*rp)))
      enddo
      enddo
c========================================================================X twist
      else if(jnk.eq.30) then
         do i=1,(inoe(k,1)+1)*2
         do j=1,3
         fco(i,j)=0.
         enddo
         enddo
      ic=0
      vtot=0.
      do l=1,inoe(k,1)
      i1=icon(ic+1)
      i2=icon(ic+2)
      i3=icon(ic+3)
      i4=icon(ic+4)
      do j=1,3
      v1(j)=corm(i3,j)-corm(i1,j)
      v2(j)=corm(i4,j)-corm(i2,j)
      s1(j)=corm(i2,j)-corm(i1,j)
      s2(j)=corm(i4,j)-corm(i3,j)
      enddo
      rv1=sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
      rv2=sqrt(v2(1)**2+v2(2)**2+v2(3)**2)
      do j=1,3
      t(j)=v1(j)/rv1+v2(j)/rv2
      enddo
      rt=sqrt(t(1)*t(1)+t(2)*t(2)+t(3)*t(3))
      t(1)=t(1)/rt
      t(2)=t(2)/rt
      t(3)=t(3)/rt
      s1t=s1(1)*t(1)+s1(2)*t(2)+s1(3)*t(3)
      s2t=s2(1)*t(1)+s2(2)*t(2)+s2(3)*t(3)
      do j=1,3
      q1(j)=s1(j)-t(j)*s1t
      q2(j)=s2(j)-t(j)*s2t
      enddo
      rq1=sqrt(q1(1)**2+q1(2)**2+q1(3)**2)
      rq2=sqrt(q2(1)**2+q2(2)**2+q2(3)**2)
      dotq=q1(1)*q2(1)+q1(2)*q2(2)+q1(3)*q2(3)
      val1=dotq/(rq1*rq2)
      if(abs(val1).gt.1) val1=sign(1.d0,val1)
      val=acos(val1)*crd
         px=q1(2)*q2(3)-q1(3)*q2(2)
         py=q1(3)*q2(1)-q1(1)*q2(3)
         pz=q1(1)*q2(2)-q1(2)*q2(1)
         dpt=px*t(1)+py*t(2)+pz*t(3)
         if(dpt.lt.0) val=-val
      vtot=vtot+val
      if(abs(val1).eq.1.d0) then
      fpen=0.
      else
      fpen=sign(1.d0,dpt)*crd/sqrt(1-val1*val1)
      endif
      do j=1,3
      dv1(1,j)=-1.
      dv1(2,j)= 0.
      dv1(3,j)= 1.
      dv1(4,j)= 0.
      dv2(1,j)= 0.
      dv2(2,j)=-1.
      dv2(3,j)= 0.
      dv2(4,j)= 1.
      enddo
      do i=1,4
      do j=1,3
      drv1(i,j)=-dv1(i,j)*v1(j)/rv1**3
      drv2(i,j)=-dv2(i,j)*v2(j)/rv2**3
      enddo
      enddo
      do j=1,3
      do m=1,3
      dt(1,j,m)=v1(j)*drv1(1,m)
      dt(3,j,m)=v1(j)*drv1(3,m)
      dt(2,j,m)=v2(j)*drv2(2,m)
      dt(4,j,m)=v2(j)*drv2(4,m)
      enddo
      dt(1,j,j)=dt(1,j,j)+dv1(1,j)/rv1
      dt(3,j,j)=dt(3,j,j)+dv1(3,j)/rv1
      dt(2,j,j)=dt(2,j,j)+dv2(2,j)/rv2
      dt(4,j,j)=dt(4,j,j)+dv2(4,j)/rv2
      enddo
         do i=1,4
         do j=1,3
         dtm=0.
         do m=1,3
         dtm=dtm+dt(i,j,m)*t(m)
         enddo
         do m=1,3
         dt(i,j,m)=dt(i,j,m)/rt-t(m)*dtm/rt
         enddo
         enddo
         enddo
      do j=1,3
      dv1(1,j)=-1.
      dv1(2,j)= 1.
      dv1(3,j)= 0.
      dv1(4,j)= 0.
      dv2(1,j)= 0.
      dv2(2,j)= 0.
      dv2(3,j)=-1.
      dv2(4,j)= 1.
      enddo
      do i=1,4
      do j=1,3
      drv1(i,j)=0.
      drv2(i,j)=0.
      do m=1,3
      drv1(i,j)=drv1(i,j)+dt(i,j,m)*s1(m)
      drv2(i,j)=drv2(i,j)+dt(i,j,m)*s2(m)
      enddo
      drv1(i,j)=drv1(i,j)+dv1(i,j)*t(j)
      drv2(i,j)=drv2(i,j)+dv2(i,j)*t(j)
      enddo
      enddo
      do i=1,4
      do j=1,3
      do m=1,3
      dq1(i,j,m)=-(dt(i,j,m)*s1t+t(m)*drv1(i,j))
      dq2(i,j,m)=-(dt(i,j,m)*s2t+t(m)*drv2(i,j))
      enddo
      dq1(i,j,j)=dq1(i,j,j)+dv1(i,j)
      dq2(i,j,j)=dq2(i,j,j)+dv2(i,j)
      enddo
      enddo
      rq13=rq1**3*rq2
      rq23=rq2**3*rq1
      dotq=q1(1)*q2(1)+q1(2)*q2(2)+q1(3)*q2(3)
      do i=1,4
      do j=1,3
      dqq1=0.
      dqq2=0.
      dqq =0.
      do m=1,3
      dqq1=dqq1+dq1(i,j,m)*q1(m)
      dqq2=dqq2+dq2(i,j,m)*q2(m)
      dqq=dqq+dq1(i,j,m)*q2(m)+dq2(i,j,m)*q1(m)
      enddo
      fco(ic+i,j)=fco(ic+i,j)+fpen*(dqq/(rq1*rq2)-dotq*(dqq1/rq13
     1 +dqq2/rq23))
      enddo
      enddo
      ic=ic+2
      enddo
         if(jsk) rnoe(k)=vtot
         vlr=vtot*cdr
         vlc=cos(vlr)
         rnr=rnoe(k)*cdr
         delt=acos(cos(rnr)*cos(vlr)+sin(rnr)*sin(vlr))*crd
         cros=cos(rnr)*sin(vlr)-sin(rnr)*cos(vlr)
         if(cros.lt.0) delt=-delt
         fnr=fnoe(k)*cdr
         epen=epen+fnr*delt**2
         fpen=2*fnr*delt
         do i=1,(inoe(k,1)+1)*2
         ii=icon(i)
         do j=1,3
         for(ii,j)=for(ii,j)+fpen*fco(i,j)
         enddo
         enddo
c=============================================================cosine torsion con
      else if(jnk.eq.10) then
      i1=inoe(k,1)
      i2=inoe(k,2)
      i3=inoe(k,3)
      i4=inoe(k,4)
      dih=torp(i1,i2,i3,i4)
      con(i)=0.
      fold=nint(rnoe(k))
      barr=fnoe(k)
      epen=epen+barr*(1.+sign(1,fold)*cos(cdr*(abs(fold)*dih)))/2.
      fpen=barr*fold*sin(cdr*(abs(fold)*dih))/2.
      call deltor(dih,fpen,i1,i2,i3,i4)
c==================================================================variable sums
      else if(jnk.eq.15.or.jnk.eq.17) then
      sumv=0.
      angle=.false.
      do iv=1,jnk-13
      ia=inoe(k,iv)
      is=sign(1,ia)
      ia=abs(ia)
         do jv=1,nvar
         if(neq(jv).eq.ia) goto 222
         enddo
222   if(jv.le.nbac.or.lar(jv)) angle=.true.
      sumv=sumv+is*var(jv)
      enddo
      delt=sumv-rnoe(k)
      con(i)=delt
      if(angle) delt=delt*cdr
      epen=epen+fnoe(k)*delt**2
c=====================================================================total rise
      else if(jnk.eq.12) then
      i1=inoe(k,1)
      i2=inoe(k,2)
      zave=0.
      do ih=i1,i2
      zave=zave+hel(ih,3)
      enddo
      if(jsk) rnoe(k)=zave/(i2-i1+1)
      delt=zave-rnoe(k)*(i2-i1+1)
      con(i)=delt
      epen=epen+fnoe(k)*delt**2
c====================================================================total twist
      else if(jnk.eq.13) then
      i1=inoe(k,1)
      i2=inoe(k,2)
      wave=0.
      do ih=i1,i2
      wave=wave+hel(ih,6)
      enddo
      if(jsk) rnoe(k)=wave/(i2-i1+1)
      delt=wave-rnoe(k)*(i2-i1+1)
      con(i)=delt
      delt=delt*cdr
      epen=epen+fnoe(k)*delt**2
c======================================================================curvature
      else if(jnk.eq.20) then
      ib1=inoe(k,1)
      ib2=inoe(k,2)
      the=rnoe(k)*cdr
      ct=cos(the)
      st=sin(the)
      ct2=cos(the/2)
      st2=sin(the/2)
      phi=rnoe(k+1)*cdr
      cp=cos(phi)
      sp=sin(phi)
      x0=ra(ib1,1)
      y0=ra(ib1,2)
      z0=ra(ib1,3)
      sx=ra(ib2,1)-x0
      sy=ra(ib2,2)-y0
      sz=ra(ib2,3)-z0
      rs=sqrt(sx*sx+sy*sy+sz*sz)
      ux=ua(ib1,1)
      uy=ua(ib1,2)
      uz=ua(ib1,3)
      ex=uy*da(ib1,3)-uz*da(ib1,2)
      ey=uz*da(ib1,1)-ux*da(ib1,3)
      ez=ux*da(ib1,2)-uy*da(ib1,1)
      wx=da(ib1,1)*cp+ex*sp
      wy=da(ib1,2)*cp+ey*sp
      wz=da(ib1,3)*cp+ez*sp
      s0x=ux*ct2+wx*st2
      s0y=uy*ct2+wy*st2
      s0z=uz*ct2+wz*st2
         q1x=x0+s0x*rs
         q1y=y0+s0y*rs
         q1z=z0+s0z*rs
      v0x=ux*ct+wx*st
      v0y=uy*ct+wy*st
      v0z=uz*ct+wz*st
         q2x=q1x+v0x
         q2y=q1y+v0y
         q2z=q1z+v0z
      d1x=ra(ib2,1)-q1x
      d1y=ra(ib2,2)-q1y
      d1z=ra(ib2,3)-q1z
      d1=sqrt(d1x*d1x+d1y*d1y+d1z*d1z)
      d2x=ra(ib2,1)+ua(ib2,1)-q2x
      d2y=ra(ib2,2)+ua(ib2,2)-q2y
      d2z=ra(ib2,3)+ua(ib2,3)-q2z
      d2=sqrt(d2x*d2x+d2y*d2y+d2z*d2z)
      fac=fnoe(k)
      con(i)=d1+d2
      epen=epen+fac*(d1**2+d2**2)
      dot1=s0x*d1x+s0y*d1y+s0z*d1z
      dot2=s0x*d2x+s0y*d2y+s0z*d2z
      dot=2*fac*(dot1+dot2)/rs
      fb1x=dot*sx-2*fac*d1x
      fb1y=dot*sy-2*fac*d1y
      fb1z=dot*sz-2*fac*d1z
      fb2x=-2*fac*d2x
      fb2y=-2*fac*d2y
      fb2z=-2*fac*d2z
      fbx=fb1x+fb2x
      fby=fb1y+fb2y
      fbz=fb1z+fb2z
c===================================================================groove width
      else if(jnk.eq.31.or.jnk.eq.32) then
      do ks=1,inoe(k,9)
      rmin=10000.
      m=nuc(inoe(k,ks)-1)+12
      x0=corm(m,1)
      y0=corm(m,2)
      z0=corm(m,3)
         do ls=iend(1)+1,iend(2)
         l=nuc(ls-1)+12
         dx=corm(l,1)-x0
         dy=corm(l,2)-y0
         dz=corm(l,3)-z0
         r=sqrt(dx*dx+dy*dy+dz*dz)
            if(r.lt.rmin) then
            rmin=r
            lm=l
            xx=dx
            yy=dy
            zz=dz
            endif
         enddo
      r=rmin
      if(jsk) rnoe(k)=r
      if(r.lt.rnoe(k).or.jnk.eq.31) then
      dist=r-rnoe(k)
      else if(r.gt.bnoe(k)) then
      dist=r-bnoe(k)
      else
      goto 260
      endif
      con(i)=dist
      epen=epen+fnoe(k)*dist**2
      fpen=2*fnoe(k)*dist/r
      for(m,1)=for(m,1)+xx*fpen
      for(m,2)=for(m,2)+yy*fpen
      for(m,3)=for(m,3)+zz*fpen
      for(lm,1)=for(lm,1)-xx*fpen
      for(lm,2)=for(lm,2)-yy*fpen
      for(lm,3)=for(lm,3)-zz*fpen
260   enddo
c-------------------------------------------------------------------------------
      endif
250   enddo
      return
      end
