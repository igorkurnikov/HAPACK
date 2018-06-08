      function dis(i1,i2)
      include 'jumna_data.inc'
      character*4 mnam,munit
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      dx=corm(i1,1)-corm(i2,1)
      dy=corm(i1,2)-corm(i2,2)
      dz=corm(i1,3)-corm(i2,3)
      dis=sqrt(dx*dx+dy*dy+dz*dz)
      return
      end
