      subroutine build(ctot)
      include 'jumna_data.inc'
      logical*2 ribose,cation,kink,lthy,locr,sup,rcom,homo,homo2,homo3,
     1 diep,sum,link,ecen,cyl,lcat,cent,autos,ihl,amber
      character*4 mnam,munit,snam,suni,sub,seq*120,knam,lnam,
     1 mac*32,lmo*32,lib*32,libm*32,out*32,axe*32,noe*32,nol*32,axl*32,
     1 test*32,pdb*32,ins*32,bar*32,parm*32,code*8,kode*8
      integer*4 opt
      real*8 ktl,kpr,ha(n2,9),cdir(9),drib(12,3)
      common/cha/mac,lmo,lib,libm,out,axe,axl,noe,nol,test,pdb,ins,
     1 bar,parm
      common/datjm/acc,phos,epsi,epsr,plat,slope,rhbl,vfac,tfac,rfac,xfac,
     1 scale,rpiv,tpiv,fad,fan,damp,fnoew,fnoes,fnoea,df1,scnb,scee,
     1 catd,catr,catc,rad,pit,enit,nshel,maxn,opt,mhomo,mhomo2,mhomo3,
     1 nick,limit,nop,nrib,ncat,lig,naxo,sup,rcom,homo,homo2,homo3,diep,
     1 sum,link,ecen,cyl,lcat,cent,autos,amber
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/hel/ua(n2,3),da(n2,3),ra(n2,3)
      common/ind/iend(0:4),ise(n2),ienb(n2,2),kapt(n2+1),nsph
      common/lgd/rlig(n9,6),slig(n9,2),ilig(n9,6),lopt(n9),lpiv(n9),
     1 ntlg,nlgi,nlig,ntl,lnam(n9),locr(n2,n0)
      common/moljm/sor(n4,n5,3),smon(n4,n5),snam(n4,n5),suni(n4,n5),
     1 nuni(n4,n5),sub(n5),isch(n4,n5),isty(n4,n5),ics(n4,n5),
     1 mats(3*n4,n5),kas(n5),khs(n5),ksub
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/nmrjm/rnoe(n3),bnoe(n3),fnoe(n3),snoe,knam(200,2),ktyp(200),
     1 icon(n2*2),inoe(n3,9),ih68(n2),jnoe(n3),nnoe,ndcs,ndiv(0:99),nlin
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      common/symjm/xsum,ysum,itwl(n2),nvs(n7),ihm(n2),ksym(3),nbrk(2),
     1 isym,isur,isup,ihl(n2)
	common/hacon/iret
      data rhy/0.97/,the/107.0/,tau/-60.0/
      data drib/ 0.3053, 0.1430, 0.0544, 0.1128,-0.3352,-0.1417,-0.2798,
     1          0.0652, 0.0617,-0.5131, 0.0643, 0.0589,
     1          0.3056, 0.1432, 0.0563, 0.1270,-0.3318,-0.0135,-0.2797,
     1          0.0652, 0.0617,-0.5131, 0.0645, 0.0606,
     1          0.3052, 0.1397, 0.1356, 0.1235,-0.3327,-0.0139,-0.5148,
     1          0.0651, 0.0613,-0.5139, 0.0608, 0.0602/
      equivalence (ha(1,1),ua(1,1))
      dphos=(1.+phos)/2.
      n=0
      nuc(0)=0
c==============================================================macfile matricies
      k=0
      lr=0
      lx=nto
c
c  kr - index of the strand
c  is - absolute residue number 
c
      do kr=1,nst
        idir=sign(1,idr(kr))
        do ir=1,kseq
          k=k+1
          is=ieq(ir,kr)
          if(is.eq.0) goto 50
          lr=lr+1
10        m=irec(lr)
          do i=1,kas(m)   ! Cycle over atoms of the template
            n=n+1
            if(n.gt.n1) then
               write(6,*) '---- n1 too small ----'
               return 
            endif
            if((ribose(lr).and.parm.eq.'Flex').and.i.eq.iofs(lr)) n=n+1
            corm(n,1)=sor(i,m,1)  ! initial cartesian coordinates
            corm(n,2)=sor(i,m,2)
            corm(n,3)=sor(i,m,3)
            mnam(n)=snam(i,m)     ! atom name
            munit(n)=suni(i,m)    ! subresidue name 
            if(i.le.iofe(lr)) then   ! Number of atoms in the residue ??
               nunit(n)=lr           ! residue number
            else
               nunit(n)=lr+nuni(i,m)  ! nuni(i,m) ??? 
            endif
            imch(n)=isch(i,m)      ! element number of the atom
            imty(n)=isty(i,m)      ! force field type
            dmon(n)=smon(i,m)      ! partial charge of the atom
            if(imty(n).eq.14) dmon(n)=catc
            icm(n)=ics(i,m)        ! ?? 
          enddo
          if(cation(lr)) then
            n=n+1
            nc=n
          endif
          if(kas(m).gt.iofe(lr)) lx=lx+nuni(kas(m),m)
          nuc(lr)=n
c------------------------------------------------------------------ribose create
          if((ribose(lr).and.parm.eq.'Flex').and.itr(lr).le.6) then
            nox=nuc(lr-1)+10
            nc1=nuc(lr-1)+1
            nc2=nuc(lr-1)+2
            noh=nuc(lr-1)+iofs(lr)
            dx=corm(nox,1)-corm(nc2,1)
            dy=corm(nox,2)-corm(nc2,2)
            dz=corm(nox,3)-corm(nc2,3)
            rf=1.42/sqrt(dx*dx+dy*dy+dz*dz)
            mnam(nox)='O2'''
            imch(nox)=8
            imty(nox)=8
            icm(nox)=-icm(nox)
            corm(nox,1)=corm(nc2,1)+dx*rf
            corm(nox,2)=corm(nc2,2)+dy*rf
            corm(nox,3)=corm(nc2,3)+dz*rf
            mnam(noh)='HO2'''
            munit(noh)=suni(10,m)
            nunit(noh)=lr
            imch(noh)=1
            imty(noh)=1
            icm(noh)=0
            dmon(noh)=0.3186
            call cosdir(corm,n1,nox,nc2,nc1,cdir)
            xn= rhy*cos(cdr*(the))
            yn= rhy*sin(cdr*(the))*cos(cdr*(tau))
            zn=-rhy*sin(cdr*(the))*sin(cdr*(tau))
            corm(noh,1)=cdir(1)*xn+cdir(2)*yn+cdir(3)*zn+corm(nox,1)
            corm(noh,2)=cdir(4)*xn+cdir(5)*yn+cdir(6)*zn+corm(nox,2)
            corm(noh,3)=cdir(7)*xn+cdir(8)*yn+cdir(9)*zn+corm(nox,3)
            loc=0
            do i=nuc(lr-1)+1,nuc(lr-1)+12
              loc=loc+1
              dmon(i)=drib(loc,itr(lr))
            enddo
         endif
c------------------------------------------------------------------cation create
         if(cation(lr)) then
           np=nuc(lr-1)+15
           no1=nuc(lr-1)+17 
           no2=nuc(lr-1)+18
           x0=corm(np,1)
           y0=corm(np,2)
           z0=corm(np,3)
           x1=corm(no1,1)-x0
           y1=corm(no1,2)-y0
           z1=corm(no1,3)-z0
           r1=sqrt(x1*x1+y1*y1+z1*z1)
           x2=corm(no2,1)-x0
           y2=corm(no2,2)-y0
           z2=corm(no2,3)-z0
           r2=sqrt(x2*x2+y2*y2+z2*z2)
           dx=x1/r1+x2/r2
           dy=y1/r1+y2/r2
           dz=z1/r1+z2/r2
           r=catd/sqrt(dx*dx+dy*dy+dz*dz)
           corm(nc,1)=x0+dx*r
           corm(nc,2)=y0+dy*r
           corm(nc,3)=z0+dz*r
           munit(nc)=suni(15,m)
           nunit(nc)=lr
           imch(nc)=11
           dmon(nc)=catc
           imty(nc)=14
           icm(nc)=0
           mnam(nc)='Na'
         endif
50      enddo
      enddo
c------------------------------------------------------------------------Ligands
      do is=nto+1,ntl
        m=irec(is)
        lr=lr+1
        do i=1,kas(m)
          n=n+1
          corm(n,1)=sor(i,m,1)
          corm(n,2)=sor(i,m,2)
          corm(n,3)=sor(i,m,3)
          mnam(n)=snam(i,m)
          munit(n)=suni(i,m)
          nunit(n)=lx+nuni(i,m)
          imch(n)=isch(i,m)
          imty(n)=isty(i,m)
          dmon(n)=smon(i,m)
          if(imty(n).eq.14) dmon(n)=catc
          icm(n)=ics(i,m)
        enddo
        lx=lx+nuni(kas(m),m)
        nuc(lr)=n
      enddo
      kam=n
      kcen=lr
c--------------------------------------------------------------make linear b.mat
      l=0
      do i=1,nto
        in=itr(i)
        ioff=nuc(i-1)
        iu=irec(i)       ! residue template number
        istr=ilq(i,2)    ! strand number
        idir=sign(1,idr(istr))  ! strand direction
        do m=1,khs(iu)   ! cycle on mats - internal coordinates(graph?) - read from res description
          l=l+1
          me=mats(m,iu)  
          if(me.ge.10000) then
            matm(l)=me
          else if(me.ne.0) then
            idelt=0
            if((ribose(i).and.parm.eq.'Flex').and.abs(me).ge.iofs(i)) idelt=idelt+1  ! if atom belongs to nucleotide take into account offset 
	                                                                               ! due to the presence of extra atom on ribose
            matm(l)=me+sign(ioff+idelt,me)      ! index of the atom (internal coord belongs to?) 
	                                          ! with some sign info? and account of extra atom in ribose
          else
            if(idir.eq.1)  matm(l)=nuc(i-2)+16  ! O5' Atom of the previous residue in strand
            if(idir.eq.-1) matm(l)=nuc(i)+16
          endif
        enddo
c---------------------------------------------------------------ribose extension
        if((ribose(i).and.parm.eq.'Flex').and.in.le.6) then
          l=l+1
          matm(l)=10000
          l=l+1
          matm(l)=ioff+10
          l=l+1
          matm(l)=ioff+iofs(i)
          itr(i)=itr(i)+3
          ito(i)=itr(i)
          iofs(i)=iofs(i)+1
          iofe(i)=iofe(i)+1
        endif
c---------------------------------------------------------------cation extension
        if(cation(i)) then
          l=l+1
          matm(l)=10000
          l=l+1
          matm(l)=nuc(i)
        endif
c-------------------------------------------------------------------------------
        l=l+1
        matm(l)=10000
      enddo
c------------------------------------------------------------------------Ligands
      do i=nto+1,ntl
      in=itr(i)
      ioff=nuc(i-1)
      iu=irec(i)
      do m=1,khs(iu)
      l=l+1
      me=mats(m,iu)
      if(me.ge.10000) then
      matm(l)=me
      else
      matm(l)=me+sign(ioff,me)
      endif
      enddo
      l=l+1
      matm(l)=10000
      enddo
      matm(l)=0
      khm=l-1
      if(khm.gt.3*n1) then
        write(6,100) khm-3*n1
100     format(/2x,'---- Overflow of Matm by ',i4,' elements ----'/)
        iret = 1
	  return
      endif
      call setd
c-------------------------------------------------------------------------------
      do is=1,nto
      i1=0
      i2=0
      do i=nuc(is-1)+1,nuc(is)
      if(mnam(i).eq.'O1P'.and.matd(i,7).eq.1) i1=i
      if(mnam(i).eq.'O2P'.and.matd(i,7).eq.1) i2=i
      enddo
      if(i1.ne.0.and.i2.ne.0) then
      dmon(i1)=dmon(i1)+dphos
      dmon(i2)=dmon(i2)+dphos
      endif
      enddo
      ctot=0.
      do i=1,kam
      ctot=ctot+dmon(i)
      enddo
c===============================================================helical geometry
      do is=1,nst
        wtot=0.0
        ztot=0.0
        do k=1,kseq
          ks=ieq(k,is)
          if(ks.ne.0) then
            xsh=hel(ks,1)
            ysh=hel(ks,2)
            ztot=ztot+hel(ks,3)
            tlt=hel(ks,4)
            pro=hel(ks,5)
            wtot=wtot+hel(ks,6)
            ctl=cos(cdr*(tlt))
            stl=sin(cdr*(tlt))
            cpr=cos(cdr*(pro))
            spr=sin(cdr*(pro))
            cwd=cos(cdr*(wtot))
            swd=sin(cdr*(wtot))
            ua(ks,1)=0.
            ua(ks,2)=0.
            ua(ks,3)=1.
            da(ks,1)=cwd
            da(ks,2)=swd
            da(ks,3)=0.
            ra(ks,1)=0.
            ra(ks,2)=0.
            ra(ks,3)=ztot
            do i=nuc(ks-1)+1,nuc(ks)
c-------------------------------------------------------------------first strand
              if(is.eq.1) then
                x0=corm(i,1)
                y0=corm(i,2)
                z0=corm(i,3)
                xp= x0*cpr+z0*spr
                yp= y0
                zp=-x0*spr+z0*cpr
                xt= xp+xsh
                yt= yp*ctl+zp*stl-ysh
                zt=-yp*stl+zp*ctl
                corm(i,1)=xt*cwd-yt*swd
                corm(i,2)=xt*swd+yt*cwd
                corm(i,3)=zt+ztot
              else
c-----------------------------------------------------------------second strands
                x0= corm(i,1)
                y0=-corm(i,2)
                z0=-corm(i,3)
                xp= x0*cpr-z0*spr
                yp= y0
                zp= x0*spr+z0*cpr
                xt= xp+xsh
                yt= yp*ctl+zp*stl+ysh
                zt=-yp*stl+zp*ctl
                corm(i,1)=xt*cwd-yt*swd
                corm(i,2)=xt*swd+yt*cwd
                corm(i,3)=zt+ztot
              endif
c-------------------------------------------------------------------------------
            enddo
          endif
        enddo
      enddo
c======================================================================kink axis
      do is=kseq,2,-1
       if(kink(is)) then
         dsx=vkink(is,1)
         dsy=vkink(is,2)
         ktl=vkink(is,3)
         kpr=vkink(is,4)
         ctl=cos(cdr*(ktl))
         stl=-sin(cdr*(ktl))
         cpr=cos(cdr*(kpr))
         spr=sin(cdr*(kpr))
         rx=da(is-1,1)
         ry=da(is-1,2)
         rz=da(is-1,3)
         wx=-ry
         wy= rx
         wz= 0.
         dox=rx*dsx-wx*dsy
         doy=ry*dsx-wy*dsy
         doz=ra(is-1,3)
         do ks=1,nst
           do js=is,kseq
             lr=ieq(js,ks)
             if(lr.ne.0) then
               do j=nuc(lr-1)+1,nuc(lr)
                 xx=corm(j,1)
                 yy=corm(j,2)
                 zz=corm(j,3)-doz
                 xx1=(wx*wx+(1-wx*wx)*cpr)*xx+(wx*wy*(1-cpr)-wz*spr)*yy+
     1               (wx*wz*(1-cpr)+wy*spr)*zz
                 yy1=(wx*wy*(1-cpr)+wz*spr)*xx+(wy*wy+(1-wy*wy)*cpr)*yy+
     1               (wy*wz*(1-cpr)-wx*spr)*zz
                 zz1=(wx*wz*(1-cpr)-wy*spr)*xx+(wy*wz*(1-cpr)+wx*spr)*yy+
     1               (wz*wz+(1-wz*wz)*cpr)*zz
                 corm(j,1)=(rx*rx+(1-rx*rx)*ctl)*xx1+(rx*ry*(1-ctl)-rz*stl)*yy1+
     1                    (rx*rz*(1-ctl)+ry*stl)*zz1+dox
                 corm(j,2)=(rx*ry*(1-ctl)+rz*stl)*xx1+(ry*ry+(1-ry*ry)*ctl)*yy1+
     1                     (ry*rz*(1-ctl)-rx*stl)*zz1+doy
                 corm(j,3)=(rx*rz*(1-ctl)-ry*stl)*xx1+(ry*rz*(1-ctl)+rx*stl)*yy1+
     1                     (rz*rz+(1-rz*rz)*ctl)*zz1+doz
               enddo
c---------------------------------------------------------update helical vectors
               do l=1,3
                 ls=(l-1)*3+1
                 xx=ha(lr,ls)
                 yy=ha(lr,ls+1)
                 zz=ha(lr,ls+2)
                 if(l.eq.3) zz=zz-doz
                 xx1=(wx*wx+(1-wx*wx)*cpr)*xx+(wx*wy*(1-cpr)-wz*spr)*yy+
     1               (wx*wz*(1-cpr)+wy*spr)*zz
                 yy1=(wx*wy*(1-cpr)+wz*spr)*xx+(wy*wy+(1-wy*wy)*cpr)*yy+
     1               (wy*wz*(1-cpr)-wx*spr)*zz
                 zz1=(wx*wz*(1-cpr)-wy*spr)*xx+(wy*wz*(1-cpr)+wx*spr)*yy+
     1               (wz*wz+(1-wz*wz)*cpr)*zz
                 ha(lr,ls)  =(rx*rx+(1-rx*rx)*ctl)*xx1+(rx*ry*(1-ctl)-rz*stl)*yy1+
     1                       (rx*rz*(1-ctl)+ry*stl)*zz1
                 ha(lr,ls+1)=(rx*ry*(1-ctl)+rz*stl)*xx1+(ry*ry+(1-ry*ry)*ctl)*yy1+
     1                       (ry*rz*(1-ctl)-rx*stl)*zz1
                 ha(lr,ls+2)=(rx*rz*(1-ctl)-ry*stl)*xx1+(ry*rz*(1-ctl)+rx*stl)*yy1+
     1                       (rz*rz+(1-rz*rz)*ctl)*zz1
                 if(l.eq.3) then
                     ha(lr,ls)  =ha(lr,ls)+dox
                     ha(lr,ls+1)=ha(lr,ls+1)+doy
                     ha(lr,ls+2)=ha(lr,ls+2)+doz
                 endif
               enddo
             endif
           enddo
         enddo
       endif
      enddo  ! end kink coordinates
c--------------------------------------------------------------superhelical bend
      if(isup.ne.0.and.isur.eq.0) then
        write(6,*) '  ---- Set superhelix rad as well as pit ----'
        iret = 1
	  return
      endif
      if(isur.ne.0) then
         if(homo.or.homo2.or.homo3) then
           write(6,*) '  ---- No homo with superhelicity ----'
           iret = 1
	     return  
         endif
         if(nbrk(1)+nbrk(2).gt.0) then
           write(6,*) '  ---- No breaks with superhelicity ----'
           iret = 1
	     return
         endif
         root=sqrt((2*pi*rad)**2+pit**2)
         ct=2*pi*rad/root
         st=pit/root
         klim=ksym(1)
         if(klim.eq.0) klim=kseq
c      xsum=0.
c      ysum=0.
         do ks=1,nst
           dz=0.
           dw=0.
           dt=0.
           do is=1,kseq
             lr=ieq(is,ks)
             if(lr.ne.0) then
               dz=dz+hel(lr,3)
               dw=dw+hel(lr,6)
               if(dw.ge.360) dw=dw-360
               the=2*pi*hel(lr,3)/root
               zed=pit*dz/root
               dt=dt+the
               cw=cos(dt)
               sw=sin(dt)
               do j=nuc(lr-1)+1,nuc(lr)
                 yy=corm(j,2)
                 zz=corm(j,3)-dz
                 yy1= ct*yy+st*zz
                 zz1=-st*yy+ct*zz
                 xx=corm(j,1)-rad
                 corm(j,1)=cw*xx+sw*zz1+rad
                 corm(j,2)=yy1+zed
                 corm(j,3)=-sw*xx+cw*zz1
               enddo
               hel(lr,6)=dw
c         if(ks.eq.1.and.is.le.klim) then
c         x=hel(is,1)
c         y=hel(is,2)
c         w=hel(is,6)*cdr
c         cwd=cos(w)
c         swd=sin(w)
c         xsum=xsum+x*cwd+y*swd
c         ysum=ysum+x*swd-y*cwd
c         endif
c---------------------------------------------------------update helical vectors
               do l=1,3
                 ls=(l-1)*3+1
                 xx=ha(lr,ls)
                 yy=ha(lr,ls+1)
                 zz=ha(lr,ls+2)
                 if(l.eq.3) then
                   xx=xx-rad
                   zz=zz-dz
                 endif
                 yy1= ct*yy+st*zz
                 zz1=-st*yy+ct*zz
                 ha(lr,ls)  = cw*xx+sw*zz1
                 ha(lr,ls+1)=yy1
                 ha(lr,ls+2)=-sw*xx+cw*zz1
                 if(l.eq.3) then
                   ha(lr,ls)=ha(lr,ls)+rad
                   ha(lr,ls+1)=ha(lr,ls+1)+zed
                 endif
               enddo
             endif
           enddo
         enddo
      endif
      if(amber)then
        do i=1,nto
          if(ribose(i)) then
            itr(i)=itr(i)+3
            ito(i)=itr(i)
            iofs(i)=iofs(i)+1
          endif
        enddo
        call kapa
      endif
      return
      end
