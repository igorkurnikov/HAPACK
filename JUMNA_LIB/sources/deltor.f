      subroutine deltor(tora,fota,j1,j2,j3,j4)
      include 'jumna_data.inc'
      character*4 mnam,munit
      integer*2 i23,i34,elim
      common/enf/for(n1,3),tor(n1,3),fot(n6),ener,elec,
     1 repl,disp,eang,etog,epen,i23(8*n1),i34(8*n1),elim(n1)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      if(tora.eq.0) return
      x12=corm(j1,1)-corm(j2,1)
      y12=corm(j1,2)-corm(j2,2)
      z12=corm(j1,3)-corm(j2,3)
      x23=corm(j3,1)-corm(j2,1)
      y23=corm(j3,2)-corm(j2,2)
      z23=corm(j3,3)-corm(j2,3)
      x34=corm(j4,1)-corm(j3,1)
      y34=corm(j4,2)-corm(j3,2)
      z34=corm(j4,3)-corm(j3,3)
      ax=-(y12*z23-z12*y23)
      ay=-(z12*x23-x12*z23)
      az=-(x12*y23-y12*x23)
      bx=-(y34*z23-z34*y23)
      by=-(z34*x23-x34*z23)
      bz=-(x34*y23-y34*x23)
      a2=ax*ax+ay*ay+az*az
      b2=bx*bx+by*by+bz*bz
      dot=ax*bx+ay*by+az*bz
      a=sqrt(a2)
      b=sqrt(b2)
      cc=fota/sin(cdr*(tora))
      f1=cc/(a2*a*b)
      f2=cc/(a*b*b2)
      a1x=f1*(a2*bx-dot*ax)
      a1y=f1*(a2*by-dot*ay)
      a1z=f1*(a2*bz-dot*az)
      a2x=f2*(b2*ax-dot*bx)
      a2y=f2*(b2*ay-dot*by)
      a2z=f2*(b2*az-dot*bz)
      for(j1,1)=for(j1,1)+y23*a1z-z23*a1y
      for(j1,2)=for(j1,2)+z23*a1x-x23*a1z
      for(j1,3)=for(j1,3)+x23*a1y-y23*a1x
      for(j2,1)=for(j2,1)-(y23-y12)*a1z+(z23-z12)*a1y+(y34*a2z-z34*a2y)
      for(j2,2)=for(j2,2)-(z23-z12)*a1x+(x23-x12)*a1z+(z34*a2x-x34*a2z)
      for(j2,3)=for(j2,3)-(x23-x12)*a1y+(y23-y12)*a1x+(x34*a2y-y34*a2x)
      for(j3,1)=for(j3,1)-(y23+y34)*a2z+(z23+z34)*a2y-(y12*a1z-z12*a1y)
      for(j3,2)=for(j3,2)-(z23+z34)*a2x+(x23+x34)*a2z-(z12*a1x-x12*a1z)
      for(j3,3)=for(j3,3)-(x23+x34)*a2y+(y23+y34)*a2x-(x12*a1y-y12*a1x)
      for(j4,1)=for(j4,1)+y23*a2z-z23*a2y
      for(j4,2)=for(j4,2)+z23*a2x-x23*a2z
      for(j4,3)=for(j4,3)+x23*a2y-y23*a2x
      return
      end
