      subroutine enecyl(cyl)
      include 'jumna_data.inc'
      logical*2 ifhb,lar,lock,cyl
      character*4 mnam,munit
      integer*2 i23,i34,elim
      common/cnd/rcyl,dcyl,hcyl,tcyl,eneq,ucyl(3),pcyl(3),ecyl,
     1 ac(25),bc(25),lcyl(3)
      common/enf/for(n1,3),tor(n1,3),fot(n6),ener,elec,
     1 repl,disp,eang,etog,epen,i23(8*n1),i34(8*n1),elim(n1)
      common/mnn/var(n7),gra(n7),con(n0*n2+n3),scl(n7),ncon,nvrc,
     1 ntba,nbac,nthe,nhel,ntki,nkin,ntri,ntot,nvar,nrin,
     1 lar(n7),lock(n6a)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/parjm/aij(25,25),bij(25,25),ahb(15),bhb(15),vt(35),va(25),
     1 vo(25),delmax,time0,nzsh,nwdg,ivf(35),ifhb(25,25)
      dx=cos(cdr*(tcyl))
      dy=sin(cdr*(tcyl))
      ecyl=0.
      ux=ucyl(1)
      uy=ucyl(2)
      uz=ucyl(3)
      gr1=0.
      gr2=0.
      gr3=0.
      do i=1,kam
      it=imty(i)
      x=corm(i,1)-pcyl(1)
      y=corm(i,2)-pcyl(2)
      z=corm(i,3)-pcyl(3)
      dot=x*ux+y*uy+z*uz
      rx=x-dot*ux
      ry=y-dot*uy
      rz=z-dot*uz
      r2=rx*rx+ry*ry+rz*rz
      r6=r2**3
      ecyl=ecyl-ac(it)/r6+bc(it)/(r6*r6)
      e=-6*ac(it)/(r6*r2)+12*bc(it)/(r6*r6*r2)
      fx=( rx*(1-ux*ux)-ry*uy*ux-rz*uz*ux)*e
      fy=(-rx*ux*uy+ry*(1-uy*uy)-rz*uz*uy)*e
      fz=(-rx*ux*uz-ry*uy*uz+rz*(1-uz*uz))*e
      for(i,1)=for(i,1)+fx
      for(i,2)=for(i,2)+fy
      for(i,3)=for(i,3)+fz
      tx=pcyl(1)+dot*ux
      ty=pcyl(2)+dot*uy
      gr1=gr1+(tx*fy-ty*fx)
      gr2=gr2+(dx*fx+dy*fy)
      gr3=gr3+fz
      enddo
c----------------------------------------------------------------------gradients
      if(lcyl(1).ne.0) gra(lcyl(1))=gr1*cdr
      if(lcyl(2).ne.0) gra(lcyl(2))=gr2
      if(lcyl(3).ne.0) gra(lcyl(3))=gr3
      epen=epen+ecyl
      return
      end
