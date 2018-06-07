      function ang(i1,i2,i3)
      include 'jumna_data.inc'
      character*4 mnam,munit
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      dx=corm(i1,1)-corm(i2,1)
      dy=corm(i1,2)-corm(i2,2)
      dz=corm(i1,3)-corm(i2,3)
      rd=sqrt(dx*dx+dy*dy+dz*dz)
      bx=corm(i3,1)-corm(i2,1)
      by=corm(i3,2)-corm(i2,2)
      bz=corm(i3,3)-corm(i2,3)
      rb=sqrt(bx*bx+by*by+bz*bz)
      ang=acos((dx*bx+dy*by+dz*bz)/(rd*rb))*crd
      return
      end
