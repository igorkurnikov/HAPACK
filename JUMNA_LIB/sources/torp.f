      function torp(i1,i2,i3,i4)
      include 'jumna_data.inc'
      character*4 mnam,munit
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      dx1=corm(i2,1)-corm(i1,1)
      dy1=corm(i2,2)-corm(i1,2)
      dz1=corm(i2,3)-corm(i1,3)
      dx2=corm(i3,1)-corm(i2,1)
      dy2=corm(i3,2)-corm(i2,2)
      dz2=corm(i3,3)-corm(i2,3)
      ux1=dy1*dz2-dz1*dy2
      uy1=dz1*dx2-dx1*dz2
      uz1=dx1*dy2-dy1*dx2
      dx1=corm(i4,1)-corm(i3,1)
      dy1=corm(i4,2)-corm(i3,2)
      dz1=corm(i4,3)-corm(i3,3)
      ux2=dz1*dy2-dy1*dz2
      uy2=dx1*dz2-dz1*dx2
      uz2=dy1*dx2-dx1*dy2
      ctor=(ux1*ux2+uy1*uy2+uz1*uz2)/sqrt((ux1*ux1+
     1    uy1*uy1+uz1*uz1)*(ux2*ux2+uy2*uy2+uz2*uz2))
      if(abs(ctor).gt.1.) ctor=sign(1.d0,ctor)
      torp=acos(ctor)*crd
      if(ux1*(uy2*dz2-uz2*dy2)+uy1*(uz2*dx2-ux2*dz2)+uz1*(ux2*dy2-
     1  uy2*dx2).lt.0.) torp=-torp
      return
      end
