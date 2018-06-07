      subroutine assemb
c  Assemble gradients of independent variables
c
      include 'jumna_data.inc'
      logical*2 ifhb,lar,lock,kink,lthy,ribose,sup,rcom,homo,homo2,
     1 homo3,diep,link,sum,ecen,cyl,lcat,cent,autos,locr,cation,ihl,
     1 amber
      character*4 mnam,munit,seq*120,code*8,kode*8,knam,lnam
      integer*2 i23,i34,elim
      integer*4 opt,li(n2)
      real*8 ktl,gl(6)
      common/cur/fb1x,fb1y,fb1z,fb2x,fb2y,fb2z,fbx,fby,fbz,ib1,ib2
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
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      common/symjm/xsum,ysum,itwl(n2),nvs(n7),ihm(n2),ksym(3),nbrk(2),
     1 isym,isur,isup,ihl(n2)
      if(isur.ne.0) then
        pi2=pi*2
        root=pi2*rad
        if(isup.ne.0) root=sqrt((pi2*rad)**2+pit**2)
        fac=pi2/root
        rt2=root**2
        ftr=-pi2*pit/rt2
        ftp= pi2*rad/rt2
        fhr=-4*pi**2*rad/rt2
        fhp=-pit/rt2
      endif
c-------------------------------------------------------------------------------
      if(nzsh.ne.0) then
         iuz1=inoe(nzsh,1)
         iuz2=inoe(nzsh,2)
            else
            iuz1=0
            iuz2=0
            endif
         if(nwdg.ne.0) then
         iuw1=inoe(nwdg,1)
         iuw2=inoe(nwdg,2)
            else
            iuw1=0
            iuw2=0
            endif
c-------------------------------------------------------------------------------
      k=0
      ki=0
      do is=1,ntl
      ioff=nuc(is-1)
      in=itr(is)
      do lr=1,nsr(is)
      if(.not.locr(is,lr)) then
      ll=kap(1,in)+1+(lr-1)*5
      lu=ll+4
c---------------------------------------------------------ring dep contributions
      do l=ll,lu
      lk=l+k
      ia=nap(l,1,in)+ioff
      ib=nap(l,2,in)+ioff
      ic=nap(l,3,in)+ioff
      id=nap(l,4,in)+ioff
      if(id.eq.ioff) then
      ax=corm(ib,1)-corm(ia,1)
      ay=corm(ib,2)-corm(ia,2)
      az=corm(ib,3)-corm(ia,3)
      bx=corm(ic,1)-corm(ib,1)
      by=corm(ic,2)-corm(ib,2)
      bz=corm(ic,3)-corm(ib,3)
      rx=by*az-bz*ay
      ry=bz*ax-bx*az
      rz=bx*ay-by*ax
      else
      rx=corm(ic,1)-corm(ib,1)
      ry=corm(ic,2)-corm(ib,2)
      rz=corm(ic,3)-corm(ib,3)
      endif
      rr=sqrt(rx*rx+ry*ry+rz*rz)
      rx=rx/rr
      ry=ry/rr
      rz=rz/rr
      x0=corm(ib,1)
      y0=corm(ib,2)
      z0=corm(ib,3)
      exo=0.
      jl=nuc(is)-nuc(is-1)
      do j=1,jl
      m=iap(l,j,in)
      if(cation(is).and.j.eq.jl) m=iap(l,17,in)
      if(m.ne.0) then
      jg=j+ioff
      dx=corm(jg,1)-x0
      dy=corm(jg,2)-y0
      dz=corm(jg,3)-z0
      tx=dy*for(jg,3)-dz*for(jg,2)+tor(jg,1)
      ty=dz*for(jg,1)-dx*for(jg,3)+tor(jg,2)
      tz=dx*for(jg,2)-dy*for(jg,1)+tor(jg,3)
      exo=exo+(rx*tx+ry*ty+rz*tz)/m
      endif
      enddo
      if(id.eq.ioff) then
        if(amber) then
          call delval(-exo,ia,ib,ic)
        else
          call delval(-fot(lk)-exo,ia,ib,ic)
        endif
      else
      if(amber) then
        tork=torp(ia,ib,ic,id)
        call deltor(tork,exo,ia,ib,ic,id)
      else
         call deltor(sap(lk),fot(lk)+exo,ia,ib,ic,id)
      endif
      endif
      enddo
      endif
      enddo
c-------------------------------------------------------------independent derivs
      do l=1,kap(1,in)
      k=k+1
      if(.not.lock(k)) then
      ki=ki+1
      ia=nap(l,1,in)+ioff
      ib=nap(l,2,in)+ioff
      ic=nap(l,3,in)+ioff
      id=nap(l,4,in)+ioff
      if(id.eq.ioff) then
      ax=corm(ib,1)-corm(ia,1)
      ay=corm(ib,2)-corm(ia,2)
      az=corm(ib,3)-corm(ia,3)
      bx=corm(ic,1)-corm(ib,1)
      by=corm(ic,2)-corm(ib,2)
      bz=corm(ic,3)-corm(ib,3)
      rx=by*az-bz*ay
      ry=bz*ax-bx*az
      rz=bx*ay-by*ax
      else
      rx=corm(ic,1)-corm(ib,1)
      ry=corm(ic,2)-corm(ib,2)
      rz=corm(ic,3)-corm(ib,3)
      endif
      rr=sqrt(rx*rx+ry*ry+rz*rz)
      rx=rx/rr
      ry=ry/rr
      rz=rz/rr
      x0=corm(ib,1)
      y0=corm(ib,2)
      z0=corm(ib,3)
      gra(ki)=-fot(k)
      jl=nuc(is)-nuc(is-1)
      do j=1,jl
      m=iap(l,j,in)
      if(cation(is).and.j.eq.jl) m=iap(l,17,in)
      if(m.ne.0) then
      jg=j+ioff
      dx=corm(jg,1)-x0
      dy=corm(jg,2)-y0
      dz=corm(jg,3)-z0
      tx=dy*for(jg,3)-dz*for(jg,2)+tor(jg,1)
      ty=dz*for(jg,1)-dx*for(jg,3)+tor(jg,2)
      tz=dx*for(jg,2)-dy*for(jg,1)+tor(jg,3)
      gra(ki)=gra(ki)-(rx*tx+ry*ty+rz*tz)/m
      endif
      enddo
      endif
      enddo
      k=k+kap(2,in)-kap(1,in)
      enddo
c==============================================================helical variables
      ll=0
      wdg=0.
      klim=ksym(1)
      if(klim.gt.kseq) klim=kseq
      do is=1,nto
        ioff=nuc(is-1)
        in=itr(is)
        xsh=hel(is,1)
        ysh=hel(is,2)
        tlt=hel(is,4)
        ctl=cos(cdr*(tlt))
        stl=-sin(cdr*(tlt))
        ux=ua(is,1)
        uy=ua(is,2)
        uz=ua(is,3)
        rx=da(is,1)
        ry=da(is,2)
        rz=da(is,3)
        wx=uy*rz-uz*ry
        wy=uz*rx-ux*rz
        wz=ux*ry-uy*rx
        if(is.gt.kseq) then
          wx=-wx
          wy=-wy
          wz=-wz
        endif
	  dox=ra(is,1)+rx*xsh-wx*ysh
        doy=ra(is,2)+ry*xsh-wy*ysh
        doz=ra(is,3)+rz*xsh-wz*ysh
        px=(rx*rx+(1-rx*rx)*ctl)*wx+(rx*ry*(1-ctl)-rz*stl)*wy+
     1     (rx*rz*(1-ctl)+ry*stl)*wz
        py=(rx*ry*(1-ctl)+rz*stl)*wx+(ry*ry+(1-ry*ry)*ctl)*wy+
     1     (ry*rz*(1-ctl)-rx*stl)*wz
        pz=(rx*rz*(1-ctl)-ry*stl)*wx+(ry*rz*(1-ctl)+rx*stl)*wy+
     1     (rz*rz+(1-rz*rz)*ctl)*wz
c-------------------------------------------------------------------------------
      wdg=wdg+hel(is,6)*cdr
      if(isur.ne.0) wdg=hel(is,6)*cdr
      cwd=cos(wdg)
      swd=sin(wdg)
      do l=1,6
      k=k+1
      if(.not.lock(k)) then
      ki=ki+1
      gra(ki)=0.
      if(is.le.kseq) then
      if(l.eq.3.and.is.ge.iuz1.and.is.le.iuz2)
     1 gra(ki)=gra(ki)+2*fnoea*con(nzsh+ndcs)
      if(l.eq.6.and.is.ge.iuw1.and.is.le.iuw2)
     1 gra(ki)=gra(ki)+2*fnoea*con(nwdg+ndcs)*cdr
         if(l.eq.3.and.ib1.lt.is.and.ib2.ge.is)
     1   gra(ki)=gra(ki)-(ux*fbx+uy*fby+uz*fbz)
         if(l.eq.6.and.ib1.lt.is.and.ib2.ge.is) then
         dx=ra(ib2,1)-ra(is,1)
         dy=ra(ib2,2)-ra(is,2)
         dz=ra(ib2,3)-ra(is,3)
         ax=dy*fbz-dz*fby+ua(ib2,2)*fb2z-ua(ib2,3)*fb2y
         ay=dz*fbx-dx*fbz+ua(ib2,3)*fb2x-ua(ib2,1)*fb2z
         az=dx*fby-dy*fbx+ua(ib2,1)*fb2y-ua(ib2,2)*fb2x
         gra(ki)=gra(ki)-(ux*ax+uy*ay+uz*az)
         endif
      if(cent.and.is.le.klim) then
      fax=2*fnoes*xsum
      fay=2*fnoes*ysum
      if(l.eq.1) gra(ki)=gra(ki)+fax*cwd+fay*swd
      if(l.eq.2) gra(ki)=gra(ki)+fax*swd-fay*cwd
      if(l.eq.6) then
      sut=-fax*(xsh*swd-ysh*cwd)+fay*(xsh*cwd+ysh*swd)
      gra(ki)=gra(ki)+sut
         if(isur.eq.0) then
         ll=ll+1
         li(ll)=ki
         do m=1,ll-1
         ko=li(m)
         gra(ko)=gra(ko)+sut
         enddo
         endif
      endif
      endif
      endif
      ipos=ilq(is,1)
      istr=ilq(is,2)
      itl=ipos
      itu=ipos
      kpu=istr
      if(l.eq.3.or.l.eq.6) then
      itu=ilq(ihm(is),1)
      if(is.le.kseq) kpu=nst
      endif
      do kp=istr,kpu
         if(l.eq.3.or.l.eq.6) then
         itl=ipos
         if(kp.gt.istr) then
         js=0
         if(ipos.gt.1) js=ieq(ipos-1,kp)
c         if(js.eq.0) then
c         itl=ipos+1
c         else
c         itl=ilq(ihm(js),1)+1
         if(js.gt.0) itl=ilq(ihm(js),1)+1
c         endif
         endif
         endif
      do it=itl,itu
      ip=ieq(it,kp)
      if(ip.ne.0) then
      do jg=nuc(ip-1)+1,nuc(ip)
c-------------------------------------------------------------------------------
      fx=for(jg,1)
      fy=for(jg,2)
      fz=for(jg,3)
      if(l.eq.1) then
      gra(ki)=gra(ki)-(rx*fx+ry*fy+rz*fz)
      else if(l.eq.2) then
      gra(ki)=gra(ki)+(wx*fx+wy*fy+wz*fz)
      else if(l.eq.3) then
         if(isur.eq.0) then
         gra(ki)=gra(ki)-(ux*fx+uy*fy+uz*fz)
         else
         gra(ki)=gra(ki)-(corm(jg,3)*fx-(corm(jg,1)
     1   -rad)*fz+tor(jg,2)+fy*pit/pi2)*fac
         endif
      else if(l.eq.4) then
      dx=corm(jg,1)-dox
      dy=corm(jg,2)-doy
      dz=corm(jg,3)-doz
      tx=dy*fz-dz*fy+tor(jg,1)
      ty=dz*fx-dx*fz+tor(jg,2)
      tz=dx*fy-dy*fx+tor(jg,3)
      gra(ki)=gra(ki)+(rx*tx+ry*ty+rz*tz)
      else if(l.eq.5) then
      dx=corm(jg,1)-dox
      dy=corm(jg,2)-doy
      dz=corm(jg,3)-doz
      tx=dy*fz-dz*fy+tor(jg,1)
      ty=dz*fx-dx*fz+tor(jg,2)
      tz=dx*fy-dy*fx+tor(jg,3)
      gra(ki)=gra(ki)-(px*tx+py*ty+pz*tz)
      else if(isur.eq.0.or.(kp.eq.istr.and.it.eq.itl)) then
      dx=corm(jg,1)-ra(is,1)
      dy=corm(jg,2)-ra(is,2)
      dz=corm(jg,3)-ra(is,3)
      tx=dy*fz-dz*fy+tor(jg,1)
      ty=dz*fx-dx*fz+tor(jg,2)
      tz=dx*fy-dy*fx+tor(jg,3)
      gra(ki)=gra(ki)-(ux*tx+uy*ty+uz*tz)
      endif
      enddo
c------------------------------------------------------------------------ligands
      if(ip.le.kseq.and.nlig.ne.0) then
      do il=1,nlig
      if(ip.eq.ilig(il,1).and.(l.eq.3.or.l.eq.6)) then
      ih=nto+il
      do jg=nuc(ih-1)+1,nuc(ih)
      fx=for(jg,1)
      fy=for(jg,2)
      fz=for(jg,3)
      if(l.eq.3) then
         if(isur.eq.0) then
         gra(ki)=gra(ki)-(ux*fx+uy*fy+uz*fz)
         else
         gra(ki)=gra(ki)-(corm(jg,3)*fx-(corm(jg,1)-rad)*fz
     1   +tor(jg,2)+fy*pit/pi2)*fac
         endif
      else if(isur.eq.0.or.(kp.eq.istr.and.it.eq.itl)) then
      dx=corm(jg,1)-ra(is,1)
      dy=corm(jg,2)-ra(is,2)
      dz=corm(jg,3)-ra(is,3)
      tx=dy*fz-dz*fy+tor(jg,1)
      ty=dz*fx-dx*fz+tor(jg,2)
      tz=dx*fy-dy*fx+tor(jg,3)
      gra(ki)=gra(ki)-(ux*tx+uy*ty+uz*tz)
      endif
      enddo
      endif
      enddo
      endif
c-------------------------------------------------------------------------------
      endif
      enddo
      enddo
      endif
      enddo
      enddo
c=================================================================kink variables
      do is=1,kseq
      if(kink(is)) then
      ioff=nuc(is-1)
      dsx=vkink(is,1)
      dsy=vkink(is,2)
      ktl=vkink(is,3)
      ctl=cos(cdr*(ktl))
      stl=-sin(cdr*(ktl))
      ux=ua(is-1,1)
      uy=ua(is-1,2)
      uz=ua(is-1,3)
      rx=da(is-1,1)
      ry=da(is-1,2)
      rz=da(is-1,3)
      wx=uy*rz-uz*ry
      wy=uz*rx-ux*rz
      wz=ux*ry-uy*rx
      dox=ra(is-1,1)+rx*dsx-wx*dsy
      doy=ra(is-1,2)+ry*dsx-wy*dsy
      doz=ra(is-1,3)+rz*dsx-wz*dsy
      px=(rx*rx+(1-rx*rx)*ctl)*wx+(rx*ry*(1-ctl)-rz*stl)*wy+
     1  (rx*rz*(1-ctl)+ry*stl)*wz
      py=(rx*ry*(1-ctl)+rz*stl)*wx+(ry*ry+(1-ry*ry)*ctl)*wy+
     1  (ry*rz*(1-ctl)-rx*stl)*wz
      pz=(rx*rz*(1-ctl)-ry*stl)*wx+(ry*rz*(1-ctl)+rx*stl)*wy+
     1  (rz*rz+(1-rz*rz)*ctl)*wz
c-------------------------------------------------------------------------------
      do l=1,4
      k=k+1
      if(.not.lock(k)) then
      ki=ki+1
      gra(ki)=0.
c#####
         if(ib1.lt.is.and.ib2.ge.is) then
         if(l.eq.1) then
         gra(ki)=gra(ki)-(rx*fbx+ry*fby+rz*fbz)
         else if(l.eq.2) then
         gra(ki)=gra(ki)+(wx*fbx+wy*fby+wz*fbz)
         else if(l.eq.3) then
         dx=ra(ib2,1)-dox
         dy=ra(ib2,2)-doy
         dz=ra(ib2,3)-doz
         tx=dy*fbz-dz*fby+ua(ib2,2)*fb2z-ua(ib2,3)*fb2y
         ty=dz*fbx-dx*fbz+ua(ib2,3)*fb2x-ua(ib2,1)*fb2z
         tz=dx*fby-dy*fbx+ua(ib2,1)*fb2y-ua(ib2,2)*fb2x
         gra(ki)=gra(ki)+(rx*tx+ry*ty+rz*tz)
         else
         dx=ra(ib2,1)-dox
         dy=ra(ib2,2)-doy
         dz=ra(ib2,3)-doz
         tx=dy*fbz-dz*fby+ua(ib2,2)*fb2z-ua(ib2,3)*fb2y
         ty=dz*fbx-dx*fbz+ua(ib2,3)*fb2x-ua(ib2,1)*fb2z
         tz=dx*fby-dy*fbx+ua(ib2,1)*fb2y-ua(ib2,2)*fb2x
         gra(ki)=gra(ki)-(px*tx+py*ty+pz*tz)
         endif
         endif
c####
      do kp=1,nst
      do it=is,kseq
      ip=ieq(it,kp)
      if(ip.ne.0) then
      do jg=nuc(ip-1)+1,nuc(ip)
c-------------------------------------------------------------------------------
      fx=for(jg,1)
      fy=for(jg,2)
      fz=for(jg,3)
      if(l.eq.1) then
      gra(ki)=gra(ki)-(rx*fx+ry*fy+rz*fz)
      else if(l.eq.2) then
      gra(ki)=gra(ki)+(wx*fx+wy*fy+wz*fz)
      else if(l.eq.3) then
      dx=corm(jg,1)-dox
      dy=corm(jg,2)-doy
      dz=corm(jg,3)-doz
      tx=dy*fz-dz*fy+tor(jg,1)
      ty=dz*fx-dx*fz+tor(jg,2)
      tz=dx*fy-dy*fx+tor(jg,3)
      gra(ki)=gra(ki)+(rx*tx+ry*ty+rz*tz)
      else if(l.eq.4) then
      dx=corm(jg,1)-dox
      dy=corm(jg,2)-doy
      dz=corm(jg,3)-doz
      tx=dy*fz-dz*fy+tor(jg,1)
      ty=dz*fx-dx*fz+tor(jg,2)
      tz=dx*fy-dy*fx+tor(jg,3)
      gra(ki)=gra(ki)-(px*tx+py*ty+pz*tz)
      endif
      enddo
c------------------------------------------------------------------------ligands
      if(ip.le.kseq.and.nlig.ne.0) then
      do il=1,nlig
      if(ip.eq.ilig(il,1)) then
      ih=nto+il
      do jg=nuc(ih-1)+1,nuc(ih)
      fx=for(jg,1)
      fy=for(jg,2)
      fz=for(jg,3)
      if(l.eq.1) then
      gra(ki)=gra(ki)-(rx*fx+ry*fy+rz*fz)
      else if(l.eq.2) then
      gra(ki)=gra(ki)+(wx*fx+wy*fy+wz*fz)
      else if(l.eq.3) then
      dx=corm(jg,1)-dox
      dy=corm(jg,2)-doy
      dz=corm(jg,3)-doz
      tx=dy*fz-dz*fy+tor(jg,1)
      ty=dz*fx-dx*fz+tor(jg,2)
      tz=dx*fy-dy*fx+tor(jg,3)
      gra(ki)=gra(ki)+(rx*tx+ry*ty+rz*tz)
      else if(l.eq.4) then
      dx=corm(jg,1)-dox
      dy=corm(jg,2)-doy
      dz=corm(jg,3)-doz
      tx=dy*fz-dz*fy+tor(jg,1)
      ty=dz*fx-dx*fz+tor(jg,2)
      tz=dx*fy-dy*fx+tor(jg,3)
      gra(ki)=gra(ki)-(px*tx+py*ty+pz*tz)
      endif
      enddo
      endif
      enddo
      endif
c-------------------------------------------------------------------------------
      endif
      enddo
      enddo
      endif
      enddo
      endif
      enddo
c===================================================================superhelical
      if(isur.ne.0) then
      if(isur.gt.0) gra(isur)=0.
      if(isup.gt.0) gra(isup)=0.
      do ip=1,nto
      dox=ra(ip,1)
      doy=ra(ip,2)
      doz=ra(ip,3)
      the=acos((rad-dox)/rad)
      if(doz.lt.0) the=-the
      rx=(dox-rad)/rad
      rz=doz/rad
      do jg=nuc(ip-1)+1,nuc(ip)
      dx=corm(jg,1)-dox
      dy=corm(jg,2)-doy
      dz=corm(jg,3)-doz
      fx=for(jg,1)
      fy=for(jg,2)
      fz=for(jg,3)
      px=dy*fz-dz*fy+tor(jg,1)
      pz=dx*fy-dy*fx+tor(jg,3)
      dot=px*rx+pz*rz
      dop=(corm(jg,3)*fx-(corm(jg,1)-rad)*fz+tor(jg,2)+fy*pit/(pi2))*the
      dos=rx*fx+rz*fz
      if(isur.gt.0) gra(isur)=gra(isur)-dop*fhr-dos-dot*ftr
      if(isup.gt.0) gra(isup)=gra(isup)-dop*fhp-fy*the/(pi2)-dot*ftp
      enddo
c------------------------------------------------------------ligand contribution
      if(ip.le.kseq.and.nlig.ne.0) then
      do il=1,nlig
      if(ip.eq.ilig(il,1)) then
      ih=nto+il
      do jg=nuc(ih-1)+1,nuc(ih)
      dx=corm(jg,1)-dox
      dy=corm(jg,2)-doy
      dz=corm(jg,3)-doz
      fx=for(jg,1)
      fy=for(jg,2)
      fz=for(jg,3)
      px=dy*fz-dz*fy+tor(jg,1)
      pz=dx*fy-dy*fx+tor(jg,3)
      dot=px*rx+pz*rz
      dop=(corm(jg,3)*fx-(corm(jg,1)-rad)*fz+tor(jg,2)+fy*pit/(pi2))*the
      dos=rx*fx+rz*fz
      if(isur.gt.0) gra(isur)=gra(isur)-dop*fhr-dos-dot*ftr
      if(isup.gt.0) gra(isup)=gra(isup)-dop*fhp-fy*the/(pi2)-dot*ftp
      enddo
      endif
      enddo
      endif
      enddo
      endif
c===================================================================ligand inter
      do il=1,nlig
      is=nto+il
      i1=ilig(il,1)
      kal=nuc(is)-nuc(is-1)
      if(i1.eq.0) then
      rx=0.
      ry=0.
      rz=1.
      dx=-cos(cdr*slig(il,2))
      dy=-sin(cdr*slig(il,2))
      dz=0.
      else
      rx=ua(i1,1)
      ry=ua(i1,2)
      rz=ua(i1,3)
      dx=-da(i1,1)
      dy=-da(i1,2)
      dz=-da(i1,3)
      endif
      wx=ry*dz-dy*rz
      wy=rz*dx-dz*rx
      wz=rx*dy-dx*ry
      if(kal.eq.1) then
      i=nuc(is)
      lint=3
      fx=-for(i,1)
      fy=-for(i,2)
      fz=-for(i,3)
      gl(1)=fx*dx+fy*dy+fz*dz
      gl(2)=fx*wx+fy*wy+fz*wz
      gl(3)=fx*rx+fy*ry+fz*rz
      goto 100
      else
      lint=6
      do j=1,3
      gl(j)=0.
      enddo
      frx=0.
      fry=0.
      frz=0.
      ipv=nuc(is-1)+lpiv(il)
      x0=corm(ipv,1)
      y0=corm(ipv,2)
      z0=corm(ipv,3)
      do i=nuc(is-1)+1,nuc(is)
      fx=-for(i,1)
      fy=-for(i,2)
      fz=-for(i,3)
      gl(1)=gl(1)+fx*dx+fy*dy+fz*dz
      gl(2)=gl(2)+fx*wx+fy*wy+fz*wz
      gl(3)=gl(3)+fx*rx+fy*ry+fz*rz
      xx=corm(i,1)-x0
      yy=corm(i,2)-y0
      zz=corm(i,3)-z0
      frx=frx+yy*for(i,3)-zz*for(i,2)+tor(i,1)
      fry=fry+zz*for(i,1)-xx*for(i,3)+tor(i,2)
      frz=frz+xx*for(i,2)-yy*for(i,1)+tor(i,3)
      enddo
      czt=cos(cdr*(rlig(il,6)))
      szt=sin(cdr*(rlig(il,6)))
      cxt=cos(cdr*(rlig(il,4)))
      sxt=sin(cdr*(rlig(il,4)))
      dxt=(rx*rx+(1-rx*rx)*czt)*dx+(rx*ry*(1-czt)-rz*szt)*dy+
     1    (rx*rz*(1-czt)+ry*szt)*dz
      dyt=(rx*ry*(1-czt)+rz*szt)*dx+(ry*ry+(1-ry*ry)*czt)*dy+
     1    (ry*rz*(1-czt)-rx*szt)*dz
      dzt=(rx*rz*(1-czt)-ry*szt)*dx+(ry*rz*(1-czt)+rx*szt)*dy+
     1    (rz*rz+(1-rz*rz)*czt)*dz
      wx=ry*dzt-dyt*rz
      wy=rz*dxt-dzt*rx
      wz=rx*dyt-dxt*ry
      gl(4)=-(frx*dxt+fry*dyt+frz*dzt)
      gl(6)=-(frx*rx+fry*ry+frz*rz)
      rx=dxt
      ry=dyt
      rz=dzt
      wxt=(rx*rx+(1-rx*rx)*cxt)*wx+(rx*ry*(1-cxt)-rz*sxt)*wy+
     1    (rx*rz*(1-cxt)+ry*sxt)*wz
      wyt=(rx*ry*(1-cxt)+rz*sxt)*wx+(ry*ry+(1-ry*ry)*cxt)*wy+
     1    (ry*rz*(1-cxt)-rx*sxt)*wz
      wzt=(rx*rz*(1-cxt)-ry*sxt)*wx+(ry*rz*(1-cxt)+rx*sxt)*wy+
     1    (rz*rz+(1-rz*rz)*cxt)*wz
      gl(5)=-(frx*wxt+fry*wyt+frz*wzt)
      endif
100   do j=1,lint
      k=k+1
      if(.not.lock(k)) then
      ki=ki+1
      gra(ki)=gl(j)
      endif
      enddo
      enddo
c=================================================================thymine methyl
      kth=ntba
      do is=1,nto
      ks=(ilq(is,2)-1)*kseq+ilq(is,1)
      if(lthy(is)) then
      ioff=nuc(is-1)
      in=itr(is)
      k=k+1
      kth=kth+1
      if(.not.lock(k)) then
      ki=ki+1
      gra(ki)=-fot(kth)
      ib=ithy(2)+iofs(is)+ioff
      ic=ithy(3)+iofs(is)+ioff
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
      do j=4,6
      jg=ithy(j)+iofs(is)+ioff
      dx=corm(jg,1)-x0
      dy=corm(jg,2)-y0
      dz=corm(jg,3)-z0
      tx=dy*for(jg,3)-dz*for(jg,2)+tor(jg,1)
      ty=dz*for(jg,1)-dx*for(jg,3)+tor(jg,2)
      tz=dx*for(jg,2)-dy*for(jg,1)+tor(jg,3)
      gra(ki)=gra(ki)-(rx*tx+ry*ty+rz*tz)
      enddo
      endif
      endif
      enddo
c---------------------------------------------------------variable sum gradients
      do k=1,nnoe
      i=ndcs+k
      if(jnoe(k).eq.15.or.jnoe(k).eq.17) then
      fpen=2*fnoe(k)*con(i)
      do iv=1,jnoe(k)-13
      ia=inoe(k,iv)
      is=sign(1,ia)
      ia=abs(ia)
         do jv=1,nvar
         if(neq(jv).eq.ia) goto 222
         enddo
222   if(jv.le.nbac.or.lar(jv)) then
      gra(jv)=gra(jv)+is*fpen*cdr
         else
         gra(jv)=gra(jv)+is*fpen
         endif
      enddo
      endif
      enddo
c------------------------------------------------------------correct for degrees
      do k=1,nbac
      gra(k)=gra(k)*cdr
      enddo
      do k=nbac+1,nvar
      if(lar(k)) gra(k)=gra(k)*cdr
      enddo
      return
      end
