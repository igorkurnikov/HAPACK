      subroutine setd
      include 'jumna_data.inc'
      character*4 mnam,munit
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      do i=1,kam
        do j=1,7
          matd(i,j)=0
        enddo
      enddo
      do i=1,khm-1
        k=matm(i)
        kn=abs(k)
        if(k.lt.10000) then
          ik=matd(kn,7)
          kp=matm(i+1)
          kpn=abs(kp)
          if(kp.lt.10000) then
             ikp=matd(kpn,7)
             ik=ik+1
             matd(kn,ik)=kp
             matd(kn,7)=ik
             ikp=ikp+1
             matd(kpn,ikp)=sign(k,kp)
             matd(kpn,7)=ikp
          endif
          if(kp.eq.20000) then
             kpp=matm(i+2)
             kppn=abs(kpp)
             ikpp=matd(kppn,7)
             ik=ik+1
             matd(kn,ik)=kpp
             matd(kn,7)=ik
             ikpp=ikpp+1
             matd(kppn,ikpp)=sign(k,kpp)
             matd(kppn,7)=ikpp
          endif
        endif
      enddo
      return
      end
