      subroutine ligaxe(hlig)
      include 'jumna_data.inc'
      logical*2 locr,kink,lthy
      character*4 mnam,munit,lnam,seq*120,code*8,kode*8
      dimension cob(6,3),cdir(9),hlig(n9,6)
      common/hel/ua(n2,3),da(n2,3),ra(n2,3)
      common/lgd/rlig(n9,6),slig(n9,2),ilig(n9,6),lopt(n9),lpiv(n9),
     1 ntlg,nlgi,nlig,ntl,lnam(n9),locr(n2,n0)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
c---------------------------------------------------------------find link params
      do il=1,nlig
      is=nto+il
      kal=nuc(is)-nuc(is-1)
      i1=ilig(il,1)
      do i=4,6
      ii=ilig(il,i)
      if(ii.ne.0) then
      do j=1,3
      cob(i,j)=corm(ii,j)
      enddo
      else
      do j=1,3
      cob(i,j)=cob(i-1,j)
      enddo
      if(i.eq.5) then
      cob(i,1)=cob(i,1)+1
      else if(i.eq.6) then
      cob(i,1)=cob(i,1)+1
      cob(i,2)=cob(i,2)+1
      endif
      endif
      enddo
      if(i1.eq.0) then
      do l=1,3
      cob(l,3)=slig(il,1)
      enddo
      sz=sin(cdr*slig(il,2))
      cz=cos(cdr*slig(il,2))
      cob(2,1)=-cz
      cob(2,2)=-sz
      cob(3,1)= sz
      cob(3,2)=-cz
      else
      dx=-da(i1,1)
      dy=-da(i1,2)
      dz=-da(i1,2)
      wx=ua(i1,2)*dz-dy*ua(i1,3)
      wy=ua(i1,3)*dx-dz*ua(i1,1)
      wz=ua(i1,1)*dy-dx*ua(i1,2)
      do j=1,3
      cob(1,j)=ra(i1,j)
      cob(2,j)=ra(i1,j)-da(i1,j)
      enddo
      cob(3,1)=ra(i1,1)+wx
      cob(3,2)=ra(i1,2)+wy
      cob(3,3)=ra(i1,3)+wz
      endif
      call cosdir(cob,6,1,2,3,cdir)
      do i=4,6
      x=cob(i,1)-cob(1,1)
      y=cob(i,2)-cob(1,2)
      z=cob(i,3)-cob(1,3)
      cob(i,1)=cdir(1)*x+cdir(4)*y+cdir(7)*z
      cob(i,2)=cdir(2)*x+cdir(5)*y+cdir(8)*z
      cob(i,3)=cdir(3)*x+cdir(6)*y+cdir(9)*z
      enddo
      x0=cob(4,1)
      y0=cob(4,2)
      z0=cob(4,3)
      hlig(il,1)=x0
      hlig(il,2)=y0
      hlig(il,3)=z0
      if(kal.eq.1) then
      hlig(il,4)=0.
      hlig(il,5)=0.
      hlig(il,6)=0.
      goto 100
      endif
      do i=4,6
      cob(i,1)=cob(i,1)-x0
      cob(i,2)=cob(i,2)-y0
      cob(i,3)=cob(i,3)-z0
      enddo
      ax=cob(5,1)
      ay=cob(5,2)
      az=cob(5,3)
      rr=sqrt(ax*ax+ay*ay+az*az)
c------------------------------------------------------------------two atom case
      if(kal.eq.2) then
      roy=acos(az/rr)*crd-90.
      rp=sqrt(ax*ax+ay*ay)
      roz=acos(ax/rp)*crd
      if(ay.lt.0.) roz=-roz
      hlig(il,4)=0.
      hlig(il,5)=roy
      hlig(il,6)=roz
      return
      endif
c-------------------------------------------------------------------------------
      bx=cob(6,1)
      by=cob(6,2)
      bz=cob(6,3)
      cx=ay*bz-az*by
      cy=az*bx-ax*bz
      cz=ax*by-ay*bx
      bx=cy*az-cz*ay
      by=cz*ax-cx*az
      bz=cx*ay-cy*ax
      rb=sqrt(bx*bx+by*by+bz*bz)
      px=-by
      py= bx
      pz= 0.
      rp=sqrt(px*px+py*py)
      dot=px*ax+py*ay
      if(dot.lt.0.) then
      dot=-dot
      px=-px
      py=-py
      endif
      val=dot/(rp*rr)
      if(abs(val).gt.1.) val=sign(1.d0,val)
      roy=acos(val)*crd
      crdot=py*az*bx-px*az*by+(px*ay-py*ax)*bz
      if(crdot.lt.0.) roy=-roy
      val=px/rp
      if(abs(val).gt.1.) val=sign(1.d0,val)
      roz=acos(val)*crd
      if(py.lt.0.) roz=-roz
      qx=-py
      qy= px
      qz= 0.
      rq=sqrt(qx*qx+qy*qy)
      dot=qx*bx+qy*by+qz*bz
      val=dot/(rq*rb)
      if(abs(val).gt.1.) val=sign(1.d0,val)
      rox=acos(val)*crd
      crdot=qy*bz*px-qx*bz*py
      if(crdot.lt.0.) rox=-rox
      hlig(il,4)=rox
      hlig(il,5)=roy
      hlig(il,6)=roz
100   enddo
      return
      end
