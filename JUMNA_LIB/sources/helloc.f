      subroutine helloc(nop)
      include 'jumna_data.inc'
      logical*2 kink,lthy
      character*4 mnam,munit,seq*120,code*8,kode*8
      real*8 ktl,ddx(20),ddy(20),ddz(20),rnx(20),rny(20),rnz(20),
     1 kpr,rex(20),rey(20),rez(20),dew(-1:1),dez(-1:1)
      common/dob/theta(20),ind(20,3),nba(20)
      common/hel/ua(n2,3),da(n2,3),ra(n2,3)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(120,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      data th1,th2,dis/132.193d0,125.488d0,4.5033d0/
      ca=cos(cdr*(th1))
      sa=sin(cdr*(th1))
      cb=cos(cdr*(th2))
      sb=sin(cdr*(th2))
c=================================================================locate ref pts
      do ib=1,nop
      is=nba(ib)
      i2=ind(ib,2)
      i3=ind(ib,3)
      i1=0
      do i=1,matd(i2,7)
      j=abs(matd(i2,i))
      if(mnam(j).eq.'C1''') i1=j
      enddo
      if(i1.eq.0.or.i2.eq.0.or.i3.eq.0) then
      write(6,*) '  ---- HELLOC: NECESSARY BASE ATOMS NOT FOUND ----'
      stop
      endif
c--------------------------------------------------------------------find normal
      x0=corm(i2,1)
      y0=corm(i2,2)
      z0=corm(i2,3)
      ax=corm(i1,1)-x0
      ay=corm(i1,2)-y0
      az=corm(i1,3)-z0
      rr=sqrt(ax*ax+ay*ay+az*az)
      bx=corm(i3,1)-x0
      by=corm(i3,2)-y0
      bz=corm(i3,3)-z0
      rx=ay*bz-az*by
      ry=az*bx-ax*bz
      rz=ax*by-ay*bx
      r=sqrt(rx*rx+ry*ry+rz*rz)
      rx=rx/r
      ry=ry/r
      rz=rz/r
      rnx(ib)=rx
      rny(ib)=ry
      rnz(ib)=rz
c------------------------------------------------------construct reference point
      fac=dis/rr
      xx=ax*fac
      yy=ay*fac
      zz=az*fac
      rex(ib)=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     1       (rx*rz*(1-ca)+ry*sa)*zz+x0
      rey(ib)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1       (ry*rz*(1-ca)-rx*sa)*zz+y0
      rez(ib)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1       (rz*rz+(1-rz*rz)*ca)*zz+z0
c------------------------------------------------------------construct dyad axis
      xx=ax/rr
      yy=ay/rr
      zz=az/rr
      ddx(ib)=(rx*rx+(1-rx*rx)*cb)*xx+(rx*ry*(1-cb)-rz*sb)*yy+
     1       (rx*rz*(1-cb)+ry*sb)*zz
      ddy(ib)=(rx*ry*(1-cb)+rz*sb)*xx+(ry*ry+(1-ry*ry)*cb)*yy+
     1       (ry*rz*(1-cb)-rx*sb)*zz
      ddz(ib)=(rx*rz*(1-cb)-ry*sb)*xx+(ry*rz*(1-cb)+rx*sb)*yy+
     1       (rz*rz+(1-rz*rz)*cb)*zz
      ux=ua(is,1)
      uy=ua(is,2)
      uz=ua(is,3)
      if(is.gt.kseq) then
      ux=-ux
      uy=-uy
      uz=-uz
      endif
      dot=ux*ddx(ib)+uy*ddy(ib)+uz*ddz(ib)
      wx=ddx(ib)-ux*dot
      wy=ddy(ib)-uy*dot
      wz=ddz(ib)-uz*dot
      r=sqrt(wx*wx+wy*wy+wz*wz)
      wx=wx/r
      wy=wy/r
      wz=wz/r
      da(is,1)=wy*uz-wz*uy
      da(is,2)=wz*ux-wx*uz
      da(is,3)=wx*uy-wy*ux
      dx=rex(ib)-ra(is,1)
      dy=rey(ib)-ra(is,2)
      dz=rez(ib)-ra(is,3)
      dot=ux*dx+uy*dy+uz*dz
      ra(is,1)=ra(is,1)+ux*dot
      ra(is,2)=ra(is,2)+uy*dot
      ra(is,3)=ra(is,3)+uz*dot
      enddo
c=========================================================find hel param changes
      do ib=1,nop
      is=nba(ib)
      im=is-1
      ip=is+1
      fx=rnx(ib)
      fy=rny(ib)
      fz=rnz(ib)
      dx=ddx(ib)
      dy=ddy(ib)
      dz=ddz(ib)
      ux=ua(is,1)
      uy=ua(is,2)
      uz=ua(is,3)
      if(is.gt.kseq) then
      ux=-ux
      uy=-uy
      uz=-uz
      endif
      vx=da(is,1)
      vy=da(is,2)
      vz=da(is,3)
      wx=uy*vz-uz*vy
      wy=uz*vx-ux*vz
      wz=ux*vy-uy*vx
      ptx=ra(is,1)
      pty=ra(is,2)
      ptz=ra(is,3)
c-------------------------------------------------------local helical parameters
      ex=rex(ib)-ptx
      ey=rey(ib)-pty
      ez=rez(ib)-ptz
      xsh= ex*vx+ey*vy+ez*vz
      ysh=-ex*wx-ey*wy-ez*wz
      dot=wx*dx+wy*dy+wz*dz
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      tlt=acos(dot)*crd
      tx=wy*dz-wz*dy
      ty=wz*dx-wx*dz
      tz=wx*dy-wy*dx
      if(tx*vx+ty*vy+tz*vz.gt.0) tlt=-tlt
      px=vy*dz-vz*dy
      py=vz*dx-vx*dz
      pz=vx*dy-vy*dx
      rp=sqrt(px*px+py*py+pz*pz)
      dot=(px*fx+py*fy+pz*fz)/rp
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      pro=acos(dot)*crd
      ox=py*fz-pz*fy
      oy=pz*fx-px*fz
      oz=px*fy-py*fx
      if(ox*dx+oy*dy+oz*dz.lt.0) pro=-pro
      hel(is,1)=xsh
      hel(is,2)=ysh
      hel(is,4)=tlt
      hel(is,5)=pro
c===================================================================first strand
      if(is.le.kseq) then
      uxm=ua(im,1)
      uym=ua(im,2)
      uzm=ua(im,3)
      vxm=da(im,1)
      vym=da(im,2)
      vzm=da(im,3)
      wxm=uym*vzm-uzm*vym
      wym=uzm*vxm-uxm*vzm
      wzm=uxm*vym-uym*vxm
c-------------------------------------------------------------kink wdg,zsh below
      if(kink(is)) then
      ca=cos(cdr*(vkink(is,3)))
      sa=sin(cdr*(-vkink(is,3)))
      rx=vxm
      ry=vym
      rz=vzm
      gx=(rx*rx+(1-rx*rx)*ca)*wxm+(rx*ry*(1-ca)-rz*sa)*wym+
     1  (rx*rz*(1-ca)+ry*sa)*wzm
      gy=(rx*ry*(1-ca)+rz*sa)*wxm+(ry*ry+(1-ry*ry)*ca)*wym+
     1  (ry*rz*(1-ca)-rx*sa)*wzm
      gz=(rx*rz*(1-ca)-ry*sa)*wxm+(ry*rz*(1-ca)+rx*sa)*wym+
     1  (rz*rz+(1-rz*rz)*ca)*wzm
      dot=gx*wx+gy*wy+gz*wz
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      wdg=acos(dot)*crd
      rx=gy*wz-gz*wy
      ry=gz*wx-gx*wz
      rz=gx*wy-gy*wx
      if(rx*ux+ry*uy+rz*uz.lt.0.) wdg=-wdg
      rx=ptx-ra(im,1)
      ry=pty-ra(im,2)
      rz=ptz-ra(im,3)
      hel(is,3)=(uxm*rx+uym*ry+uzm*rz)/(uxm*ux+uym*uy+uzm*uz)
      hel(is,6)=wdg
c-----------------------------------------------------------normal wdg,zsh below
      else
      dot=wxm*wx+wym*wy+wzm*wz
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      wdg=acos(dot)*crd
      rx=wym*wz-wzm*wy
      ry=wzm*wx-wxm*wz
      rz=wxm*wy-wym*wx
      if(rx*ux+ry*uy+rz*uz.lt.0.) wdg=-wdg
      rx=ptx-ra(im,1)
      ry=pty-ra(im,2)
      rz=ptz-ra(im,3)
      zsh=sqrt(rx*rx+ry*ry+rz*rz)
      if(rx*ux+ry*uy+rz*uz.lt.0.) zsh=-zsh
      hel(is,3)=zsh
      hel(is,6)=wdg
      endif
      uxp=ua(ip,1)
      uyp=ua(ip,2)
      uzp=ua(ip,3)
      vxp=da(ip,1)
      vyp=da(ip,2)
      vzp=da(ip,3)
      wxp=uyp*vzp-uzp*vyp
      wyp=uzp*vxp-uxp*vzp
      wzp=uxp*vyp-uyp*vxp
c---------------------------------------------------------------------kink above
      if(kink(ip)) then
      sx=uyp*vz-uzp*vy
      sy=uzp*vx-uxp*vz
      sz=uxp*vy-uyp*vx
      rs=sqrt(sx*sx+sy*sy+sz*sz)
      dot=(sx*wx+sy*wy+sz*wz)/rs
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      ktl=acos(dot)*crd
      x=wy*sz-wz*sy
      y=wz*sx-wx*sz
      z=wx*sy-wy*sx
      if(x*vx+y*vy+z*vz.gt.0.) ktl=-ktl
      ox=vy*sz-vz*sy
      oy=vz*sx-vx*sz
      oz=vx*sy-vy*sx
      ro=sqrt(ox*ox+oy*oy+oz*oz)
      dot=(ox*uxp+oy*uyp+oz*uzp)/ro
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      kpr=acos(dot)*crd
      x=oy*uzp-oz*uyp
      y=oz*uxp-ox*uzp
      z=ox*uyp-oy*uxp
      if(x*sx+y*sy+z*sz.lt.0.) kpr=-kpr
      rx=ra(ip,1)-ptx
      ry=ra(ip,2)-pty
      rz=ra(ip,3)-ptz
      fac=(ux*rx+uy*ry+uz*rz)/(ux*uxp+uy*uyp+uz*uzp)
      tx=rx-uxp*fac
      ty=ry-uyp*fac
      tz=rz-uzp*fac
      dsx=  tx*vx+ty*vy+tz*vz
      dsy=-(tx*wx+ty*wy+tz*wz)
      vkink(ip,1)=dsx
      vkink(ip,2)=dsy
      vkink(ip,3)=ktl
      vkink(ip,4)=kpr
c------------------------------------------------------------------zsh,wdg above
      ca=cos(cdr*(ktl))
      sa=sin(cdr*(-ktl))
      rx=vx
      ry=vy
      rz=vz
      gx=(rx*rx+(1-rx*rx)*ca)*wx+(rx*ry*(1-ca)-rz*sa)*wy+
     1  (rx*rz*(1-ca)+ry*sa)*wz
      gy=(rx*ry*(1-ca)+rz*sa)*wx+(ry*ry+(1-ry*ry)*ca)*wy+
     1  (ry*rz*(1-ca)-rx*sa)*wz
      gz=(rx*rz*(1-ca)-ry*sa)*wx+(ry*rz*(1-ca)+rx*sa)*wy+
     1  (rz*rz+(1-rz*rz)*ca)*wz
      dot=gx*wxp+gy*wyp+gz*wzp
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      wdg=acos(dot)*crd
      rx=gy*wzp-gz*wyp
      ry=gz*wxp-gx*wzp
      rz=gx*wyp-gy*wxp
      if(rx*uxp+ry*uyp+rz*uzp.lt.0.) wdg=-wdg
      hel(ip,3)=fac
      hel(ip,6)=wdg
c------------------------------------------------------------------zsh,wdg above
      else
      rx=ra(ip,1)-ptx
      ry=ra(ip,2)-pty
      rz=ra(ip,3)-ptz
      zsh=sqrt(rx*rx+ry*ry+rz*rz)
      dot=rx*uxp+ry*uyp+rz*uzp
      if(dot.lt.0.) zsh=-zsh
      wdg=acos(wx*wxp+wy*wyp+wz*wzp)*crd
      rx=wy*wzp-wz*wyp
      ry=wz*wxp-wx*wzp
      rz=wx*wyp-wy*wxp
      dot=rx*uxp+ry*uyp+rz*uzp
      if(dot.lt.0.) wdg=-wdg
      hel(ip,3)=zsh
      hel(ip,6)=wdg
      endif
c==================================================================second strand
      else
      js=ieq(is,1)
      do l=-1,1
      j=js+l
      i=is-l
      vxl=da(j,1)
      vyl=da(j,2)
      vzl=da(j,3)
      vxr=da(i,1)
      vyr=da(i,2)
      vzr=da(i,3)
      dot=vxl*vxr+vyl*vyr+vzl*vzr
      if(abs(dot).gt.1.) dot=sign(1.d0,dot)
      dew(l)=acos(dot)*crd
      x=vyl*vzr-vzl*vyr
      y=vzl*vxr-vxl*vzr
      z=vxl*vyr-vyl*vxr
      if(x*ua(i,1)+y*ua(i,2)+z*ua(i,3).lt.0.) dew(l)=-dew(l)
      pxl=ra(j,1)
      pyl=ra(j,2)
      pzl=ra(j,3)
      pxr=ra(i,1)
      pyr=ra(i,2)
      pzr=ra(i,3)
      x=pxr-pxl
      y=pyr-pyl
      z=pzr-pzl
      dez(l)=sqrt(x*x+y*y+z*z)
      if(x*ua(i,1)+y*ua(i,2)+z*ua(i,3).lt.0.) dez(l)=-dez(l)
      enddo
      hel(is,3)=hel(js,3)+dez(0)-dez(-1)
      hel(im,3)=hel(js+1,3)+dez(1)-dez(0)
      hel(is,6)=hel(js,6)+dew(0)-dew(-1)
      hel(im,6)=hel(js+1,6)+dew(1)-dew(0)
      endif
      enddo
      return
      end
