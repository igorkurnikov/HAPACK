      subroutine ligput(ipl)
      include 'jumna_data.inc'
      logical*2 locr,kink,lthy
      character*4 mnam,munit,lnam,seq*120,code*8,kode*8,
     1 mac*32,lmo*32,lib*32,libm*32,out*32,axe*32,noe*32,nol*32,axl*32,
     1 test*32,pdb*32,ins*32,bar*32,parm*32
      dimension rm(3),th(3),tau(3),cob(6,3),cdir(9),hlig(n9,6),ipl(n9)
      common/cha/mac,lmo,lib,libm,out,axe,axl,noe,nol,test,pdb,ins,
     1 bar,parm
      common/hel/ua(n2,3),da(n2,3),ra(n2,3)
      common/lgd/rlig(n9,6),slig(n9,2),ilig(n9,6),lopt(n9),lpiv(n9),
     1 ntlg,nlgi,nlig,ntl,lnam(n9),locr(n2,n0)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
c---------------------------------------------------------bring ligand to origin
      do il=1,nlig
      is=nto+il
      kal=nuc(is)-nuc(is-1)
      link=lopt(il)
      if(link.ne.2) then
      ipv=nuc(is-1)+lpiv(il)
      ilig(il,4)=ipv
      ilig(il,5)=abs(matd(ipv,1))
      ilig(il,6)=abs(matd(ipv,2))
      endif
      if(link.eq.0) goto 100
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
      call cosdir(cob,6,4,5,6,cdir)
      x0=cob(4,1)
      y0=cob(4,2)
      z0=cob(4,3)
      do i=nuc(is-1)+1,nuc(is)
      x=corm(i,1)-x0
      y=corm(i,2)-y0
      z=corm(i,3)-z0
      corm(i,1)=cdir(1)*x+cdir(4)*y+cdir(7)*z
      corm(i,2)=cdir(2)*x+cdir(5)*y+cdir(8)*z
      corm(i,3)=cdir(3)*x+cdir(6)*y+cdir(9)*z
      enddo
      do i=4,6
      x=cob(i,1)-x0
      y=cob(i,2)-y0
      z=cob(i,3)-z0
      cob(i,1)=cdir(1)*x+cdir(4)*y+cdir(7)*z
      cob(i,2)=cdir(2)*x+cdir(5)*y+cdir(8)*z
      cob(i,3)=cdir(3)*x+cdir(6)*y+cdir(9)*z
      enddo
c-------------------------------------------------------------link using 6 atoms
      if(link.eq.2) then
      rm(1)= rlig(il,1)
      th(1)= rlig(il,2)
      th(2)= rlig(il,3)
      tau(1)=rlig(il,4)
      tau(2)=rlig(il,5)
      tau(3)=rlig(il,6)
      if(kal.gt.1) then
      dx1=cob(4,1)-cob(5,1)
      dy1=cob(4,2)-cob(5,2)
      dz1=cob(4,3)-cob(5,3)
      rm(2)=sqrt(dx1**2+dy1**2+dz1**2)
      else
      rm(2)=1.
      th(2)=90.
      endif
      if(kal.gt.2) then
      dx2=cob(6,1)-cob(5,1)
      dy2=cob(6,2)-cob(5,2)
      dz2=cob(6,3)-cob(5,3)
      rm(3)=sqrt(dx2**2+dy2**2+dz2**2)
      th(3)=acos((dx1*dx2+dy1*dy2+dz1*dz2)/(rm(2)*rm(3)))*crd
      else
      rm(3)=1.
      th(3)=90.
      endif
c----------------------------------------------calculate position of link points
      do i=1,3
      ii=ilig(il,i)
      do j=1,3
      cob(i,j)=corm(ii,j)
      enddo
      enddo
      do i=4,6
      bl  =rm(i-3)
      angl=th(i-3)
      tar =tau(i-3)
      call cosdir(cob,6,i-1,i-2,i-3,cdir)
      xn= bl*cos(cdr*(angl))
      yn= bl*sin(cdr*(angl))*cos(cdr*(tar))
      zn=-bl*sin(cdr*(angl))*sin(cdr*(tar))
      cob(i,1)=cdir(1)*xn+cdir(2)*yn+cdir(3)*zn+cob(i-1,1)
      cob(i,2)=cdir(4)*xn+cdir(5)*yn+cdir(6)*zn+cob(i-1,2)
      cob(i,3)=cdir(7)*xn+cdir(8)*yn+cdir(9)*zn+cob(i-1,3)
      enddo
c------------------------------------------------------------------------helical
      else if(link.eq.1) then
      i1= ilig(il,1)
      rox=rlig(il,4)
      roy=rlig(il,5)
      roz=rlig(il,6)
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
      dz=-da(i1,3)
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
c---------------------------------------------------------------setup parameters
      cx=cos(cdr*(rox))
      sx=sin(cdr*(rox))
      cy=cos(cdr*(roy))
      sy=sin(cdr*(roy))
      cz=cos(cdr*(roz))
      sz=sin(cdr*(roz))
      do i=4,6
      x=cob(i,1)
      y=cob(i,2)
      cob(i,1)=x*cz-y*sz
      cob(i,2)=x*sz+y*cz
      enddo
      rx=cz
      ry=sz
      rz=0.
      do i=4,6
      xx=cob(i,1)
      yy=cob(i,2)
      zz=cob(i,3)
      cob(i,1)=(rx*rx+(1-rx*rx)*cx)*xx+(rx*ry*(1-cx)-rz*sx)*yy+
     1         (rx*rz*(1-cx)+ry*sx)*zz
      cob(i,2)=(rx*ry*(1-cx)+rz*sx)*xx+(ry*ry+(1-ry*ry)*cx)*yy+
     1         (ry*rz*(1-cx)-rx*sx)*zz
      cob(i,3)=(rx*rz*(1-cx)-ry*sx)*xx+(ry*rz*(1-cx)+rx*sx)*yy+
     1         (rz*rz+(1-rz*rz)*cx)*zz
      enddo
      rx=-sz*cx
      ry= cz*cx
      rz= sx
      do i=4,6
      xx=cob(i,1)
      yy=cob(i,2)
      zz=cob(i,3)
      cob(i,1)=(rx*rx+(1-rx*rx)*cy)*xx+(rx*ry*(1-cy)-rz*sy)*yy+
     1         (rx*rz*(1-cy)+ry*sy)*zz+rlig(il,1)
      cob(i,2)=(rx*ry*(1-cy)+rz*sy)*xx+(ry*ry+(1-ry*ry)*cy)*yy+
     1         (ry*rz*(1-cy)-rx*sy)*zz+rlig(il,2)
      cob(i,3)=(rx*rz*(1-cy)-ry*sy)*xx+(ry*rz*(1-cy)+rx*sy)*yy+
     1         (rz*rz+(1-rz*rz)*cy)*zz+rlig(il,3)
      enddo
      do i=4,6
      x=cob(i,1)
      y=cob(i,2)
      z=cob(i,3)
      cob(i,1)=cdir(1)*x+cdir(2)*y+cdir(3)*z+cob(1,1)
      cob(i,2)=cdir(4)*x+cdir(5)*y+cdir(6)*z+cob(1,2)
      cob(i,3)=cdir(7)*x+cdir(8)*y+cdir(9)*z+cob(1,3)
      enddo
      endif
c-------------------------------------------------------------------place ligand
      call cosdir(cob,6,4,5,6,cdir)
      do i=nuc(is-1)+1,nuc(is)
      x=corm(i,1)
      y=corm(i,2)
      z=corm(i,3)
      corm(i,1)=cdir(1)*x+cdir(2)*y+cdir(3)*z+cob(4,1)
      corm(i,2)=cdir(4)*x+cdir(5)*y+cdir(6)*z+cob(4,2)
      corm(i,3)=cdir(7)*x+cdir(8)*y+cdir(9)*z+cob(4,3)
      enddo
      if(link.eq.2) then
      ilig(il,1)=ipl(il)
      ipv=nuc(is-1)+lpiv(il)
      ilig(il,4)=ipv
      ilig(il,5)=abs(matd(ipv,1))
      ilig(il,6)=abs(matd(ipv,2))
      endif
100   enddo
      call ligaxe(hlig)
      do il=1,nlig
      do j=1,6
      rlig(il,j)=hlig(il,j)
      enddo
      enddo
      return
      end
