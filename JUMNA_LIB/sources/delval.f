      subroutine delval(fpo,j1,j2,j3)
      include 'jumna_data.inc'
      character*4 mnam,munit
      integer*2 i23,i34,elim
      common/enf/for(n1,3),tor(n1,3),fot(n6),ener,elec,
     1 repl,disp,eang,etog,epen,i23(8*n1),i34(8*n1),elim(n1)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      x12=corm(j1,1)-corm(j2,1)
      y12=corm(j1,2)-corm(j2,2)
      z12=corm(j1,3)-corm(j2,3)
      x23=corm(j3,1)-corm(j2,1)
      y23=corm(j3,2)-corm(j2,2)
      z23=corm(j3,3)-corm(j2,3)
      r1=sqrt(x12*x12+y12*y12+z12*z12)
      r2=sqrt(x23*x23+y23*y23+z23*z23)
      rr=r1*r2
      dot=(x12*x23+y12*y23+z12*z23)/rr
      a=fpo/sqrt(1-dot*dot)
      dot=dot/rr
      for(j1,1)=for(j1,1)+a*(x23/rr-dot*r2/r1*x12)
      for(j1,2)=for(j1,2)+a*(y23/rr-dot*r2/r1*y12)
      for(j1,3)=for(j1,3)+a*(z23/rr-dot*r2/r1*z12)
      for(j2,1)=for(j2,1)-a*((x12+x23)/rr-dot*(r2/r1*x12+r1/r2*x23))
      for(j2,2)=for(j2,2)-a*((y12+y23)/rr-dot*(r2/r1*y12+r1/r2*y23))
      for(j2,3)=for(j2,3)-a*((z12+z23)/rr-dot*(r2/r1*z12+r1/r2*z23))
      for(j3,1)=for(j3,1)+a*(x12/rr-dot*r1/r2*x23)
      for(j3,2)=for(j3,2)+a*(y12/rr-dot*r1/r2*y23)
      for(j3,3)=for(j3,3)+a*(z12/rr-dot*r1/r2*z23)
      return
      end
