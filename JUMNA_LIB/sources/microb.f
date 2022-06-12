      subroutine microb
c
c Conformation update during minimization
c
      include 'jumna_data.inc'
      logical*2 lar,lock,kink,lthy,ribose,cation,locr,sup,rcom,homo,
     1 homo2,homo3,diep,link,ecen,cyl,lcat,cent,autos,ihl,amber,sum
      character*4 mnam,munit,knam,seq*120,code*8,kode*8,lnam
      integer*4 opt
      real*8 ktl,dlig(6)
      dimension wa(3),dh(6),dk(6),ha(n2,9),dol(3),don(3)
      common/datjm/acc,phos,epsi,epsr,plat,slope,rhbl,vfac,tfac,rfac,xfac,
     1 scale,rpiv,tpiv,fad,fan,damp,fnoew,fnoes,fnoea,df1,scnb,scee,
     1 catd,catr,catc,rad,pit,enit,nshel,maxn,opt,mhomo,mhomo2,mhomo3,
     1 nick,limit,nop,nrib,ncat,lig,naxo,sup,rcom,homo,homo2,homo3,diep,
     1 sum,link,ecen,cyl,lcat,cent,autos,amber
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/flx/sap(n6),refg,refb,refh,refv,iap(n8,n4,n5),
     1 nap(n8,7,n5),kap(3,n5)
      common/hel/ua(n2,3),da(n2,3),ra(n2,3)
      common/ind/iend(0:4),ise(n2),ienb(n2,2),kapt(n2+1),nsph
      common/lgd/rlig(n9,6),slig(n9,2),ilig(n9,6),lopt(n9),lpiv(n9),
     1 ntlg,nlgi,nlig,ntl,lnam(n9),locr(n2,n0)
      common/mnn/var(n7),gra(n7),con(n0*n2+n3),scl(n7),ncon,nvrc,
     1 ntba,nbac,nthe,nhel,ntki,nkin,ntri,ntot,nvar,nrin,
     1 lar(n7),lock(n6a)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/nmrjm/rnoe(n3),bnoe(n3),fnoe(n3),snoe,knam(200,2),ktyp(200),
     1 icon(n2*2),inoe(n3,9),ih68(n2),jnoe(n3),nnoe,ndcs,ndiv(0:99),nlin
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(120,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      common/symjm/xsum,ysum,itwl(n2),nvs(n7),ihm(n2),ksym(3),nbrk(2),
     1 isym,isur,isup,ihl(n2)
      equivalence(ha(1,1),ua(1,1))
         if(isur.ne.0) then
         root=sqrt((2*pi*rad)**2+pit**2)
         pi2=pi*2
         endif
      ki=0  ! index of non-locked independent variables in var
      k=0   ! index of all internal variables
c-------------------------------------------------------------loop over residues
      do is=1,ntl
        in=itr(is)
        ioff=nuc(is-1)
        khold=k
        kiold=ki
        nsv=nsr(is)*5
c=======================================================indep backbone variables
        do l=1,kap(1,in)
          k=k+1
          if(.not.lock(k)) then
             ki=ki+1
             del=var(ki)-sap(k)
             sap(k)=var(ki)
             do ld=kap(1,in)+nsv+1,kap(2,in)
               ip=-nap(ld,6,in)
               if(ip.eq.l) sap(khold+ld)=sap(khold+ld)+del
             enddo
             ia=nap(l,1,in)+ioff
             ib=nap(l,2,in)+ioff
             ic=nap(l,3,in)+ioff
             id=nap(l,4,in)+ioff
             if(id.eq.ioff) then  ! valence angle
                ax=corm(ib,1)-corm(ia,1)
                ay=corm(ib,2)-corm(ia,2)
                az=corm(ib,3)-corm(ia,3)
                bx=corm(ic,1)-corm(ib,1)
                by=corm(ic,2)-corm(ib,2)
                bz=corm(ic,3)-corm(ib,3)
                rx=by*az-bz*ay
                ry=bz*ax-bx*az
                rz=bx*ay-by*ax 
             else                ! torsional angle
                rx=corm(ic,1)-corm(ib,1)
                ry=corm(ic,2)-corm(ib,2)
                rz=corm(ic,3)-corm(ib,3)
             endif
             rr=sqrt(rx*rx+ry*ry+rz*rz)
             rx=rx/rr
             ry=ry/rr
             rz=rz/rr
             x0=corm(ib,1)
             y0=corm(ib,2)
             z0=corm(ib,3)
             jl=nuc(is)-nuc(is-1)
             do j=1,jl  ! cycle over residue atoms and move dependent
                m=iap(l,j,in)
                if(cation(is).and.j.eq.jl) m=iap(l,17,in)
                if(m.ne.0) then
                   jg=j+ioff
                   ca=cos(cdr*(del/m))
                   sa=sin(cdr*(del/m))
                   xx=corm(jg,1)-x0
                   yy=corm(jg,2)-y0
                   zz=corm(jg,3)-z0
                   corm(jg,1)=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     1                        (rx*rz*(1-ca)+ry*sa)*zz+x0
                   corm(jg,2)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1                        (ry*rz*(1-ca)-rx*sa)*zz+y0
                   corm(jg,3)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1                        (rz*rz+(1-rz*rz)*ca)*zz+z0
                endif
             enddo
          endif
        enddo
c--------------------------------------------------------------------------rings
        do lr=1,nsr(is) ! cycle over sugar rings of the residue
          if(.not.locr(is,lr)) then ! if sugar ring is not locked
             ll=kap(1,in)+1+(lr-1)*5
             lu=ll+4
c-------------------------------------------------------------dep ring variables
             do l=ll,lu  ! cycle over dependent angles of the sugar ring
               k=k+1
               ia=nap(l,1,in)+ioff
               ib=nap(l,2,in)+ioff
               ic=nap(l,3,in)+ioff
               id=nap(l,4,in)+ioff
               ir=nap(l,5,in)
               if(id.eq.ioff) then
                 del=ang(ia,ib,ic)-sap(k)
                 sap(k)=sap(k)+del
                 ax=corm(ib,1)-corm(ia,1)
                 ay=corm(ib,2)-corm(ia,2)
                 az=corm(ib,3)-corm(ia,3)
                 bx=corm(ic,1)-corm(ib,1)
                 by=corm(ic,2)-corm(ib,2)
                 bz=corm(ic,3)-corm(ib,3)
                 rx=by*az-bz*ay
                 ry=bz*ax-bx*az
                 rz=bx*ay-by*ax
               else
                 del=torp(ia,ib,ic,id)-sap(k)
                 sap(k)=sap(k)+del
                 do ld=kap(1,in)+nsv+1,kap(2,in)
                   ip=-nap(ld,6,in)
                   if(ip.eq.l) sap(khold+ld)=sap(khold+ld)+del
                 enddo
                 rx=corm(ic,1)-corm(ib,1)
                 ry=corm(ic,2)-corm(ib,2)
                 rz=corm(ic,3)-corm(ib,3)
               endif
               rr=sqrt(rx*rx+ry*ry+rz*rz)
               rx=rx/rr
               ry=ry/rr
               rz=rz/rr
               x0=corm(ib,1)
               y0=corm(ib,2)
               z0=corm(ib,3)
               jl=nuc(is)-nuc(is-1)
               do j=1,jl
                 m=iap(l,j,in)
                 if(cation(is).and.j.eq.jl) m=iap(l,17,in)
                 if(m.ne.0) then
                    jg=j+ioff
                    ca=cos(cdr*(del/m))
                    sa=sin(cdr*(del/m))
                    xx=corm(jg,1)-x0
                    yy=corm(jg,2)-y0
                    zz=corm(jg,3)-z0
                    corm(jg,1)=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     1                (rx*rz*(1-ca)+ry*sa)*zz+x0
                    corm(jg,2)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1                (ry*rz*(1-ca)-rx*sa)*zz+y0
                    corm(jg,3)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1                (rz*rz+(1-rz*rz)*ca)*zz+z0
                 endif
               enddo
             enddo
          else ! if sugar ring locked: 
            k=k+5
          endif
        enddo
        k=kapt(is+1) ! index of the last internal variable of the residue is
      enddo ! end loop of residues
c==============================================================helical variables
      do is=1,nto
        ioff=nuc(is-1)
        xsh=hel(is,1)
        ysh=hel(is,2)
        tlt=hel(is,4)
        do l=1,6
          k=k+1
          dh(l)=0.  ! array of changes of helical variables
          if(.not.lock(k)) then
            ki=ki+1
            dh(l)=var(ki)-hel(is,l)
            hel(is,l)=var(ki)
          endif
        enddo
c----------------------------------------------------------------nucleotide move
        wa(1)=ua(is,2)*da(is,3)-ua(is,3)*da(is,2)
        wa(2)=ua(is,3)*da(is,1)-ua(is,1)*da(is,3)
        wa(3)=ua(is,1)*da(is,2)-ua(is,2)*da(is,1)
        if(is.gt.kseq) then ! if residue belogn to 2 (or higher) chain
          wa(1)=-wa(1)
          wa(2)=-wa(2)
          wa(3)=-wa(3)
        endif
        ctl=cos(cdr*(dh(4)))
        stl=-sin(cdr*(dh(4)))
        cpr=cos(cdr*(dh(5)))
        spr=sin(cdr*(dh(5)))
        cwd=cos(cdr*(dh(6)))
        swd=sin(cdr*(dh(6)))
        ctt=cos(cdr*(tlt+dh(4)))
        stt=-sin(cdr*(tlt+dh(4)))
        do l=1,3
          dol(l)=ra(is,l)+da(is,l)*xsh-wa(l)*ysh
          don(l)=dol(l)+da(is,l)*dh(1)-wa(l)*dh(2)
        enddo
        rx=da(is,1)
        ry=da(is,2)
        rz=da(is,3)
        xx=wa(1)
        yy=wa(2)
        zz=wa(3)
        wx=(rx*rx+(1-rx*rx)*ctt)*xx+(rx*ry*(1-ctt)-rz*stt)*yy+
     1    (rx*rz*(1-ctt)+ry*stt)*zz
        wy=(rx*ry*(1-ctt)+rz*stt)*xx+(ry*ry+(1-ry*ry)*ctt)*yy+
     1    (ry*rz*(1-ctt)-rx*stt)*zz
        wz=(rx*rz*(1-ctt)-ry*stt)*xx+(ry*rz*(1-ctt)+rx*stt)*yy+
     1    (rz*rz+(1-rz*rz)*ctt)*zz
        do jg=ioff+1,nuc(is)  ! cycle over atoms of the residue
          xx=corm(jg,1)-dol(1)
          yy=corm(jg,2)-dol(2)
          zz=corm(jg,3)-dol(3)
          xx1=(rx*rx+(1-rx*rx)*ctl)*xx+(rx*ry*(1-ctl)-rz*stl)*yy+
     1        (rx*rz*(1-ctl)+ry*stl)*zz
          yy1=(rx*ry*(1-ctl)+rz*stl)*xx+(ry*ry+(1-ry*ry)*ctl)*yy+
     1        (ry*rz*(1-ctl)-rx*stl)*zz
          zz1=(rx*rz*(1-ctl)-ry*stl)*xx+(ry*rz*(1-ctl)+rx*stl)*yy+
     1        (rz*rz+(1-rz*rz)*ctl)*zz
          corm(jg,1)=(wx*wx+(1-wx*wx)*cpr)*xx1+(wx*wy*
     1        (1-cpr)-wz*spr)*yy1+(wx*wz*(1-cpr)+wy*spr)*zz1+don(1)
          corm(jg,2)=(wx*wy*(1-cpr)+wz*spr)*xx1+(wy*wy+
     1        (1-wy*wy)*cpr)*yy1+(wy*wz*(1-cpr)-wx*spr)*zz1+don(2)
          corm(jg,3)=(wx*wz*(1-cpr)-wy*spr)*xx1+(wy*wz*
     1        (1-cpr)+wx*spr)*yy1+(wz*wz+(1-wz*wz)*cpr)*zz1+don(3)
        enddo
        rx=ua(is,1)
        ry=ua(is,2)
        rz=ua(is,3)
        rax=ra(is,1)
        ray=ra(is,2)
        raz=ra(is,3)
        drx=rax+rx*dh(3)
        dry=ray+ry*dh(3)
        drz=raz+rz*dh(3)
        if(isur.ne.0) then
          dth=pi2*dh(3)/root
          zed=pit*dh(3)/root
          cw=cos(dth)
          sw=sin(dth)
          drx=rax
          dry=ray
          drz=raz
        endif
        istr=ilq(is,2)
        ipos=ilq(is,1)
        kpu=istr
        if(is.le.kseq) kpu=nst
        do kp=istr,kpu
          itl=ipos
          if(kp.gt.istr) then
            js=0
            if(ipos.gt.1) js=ieq(ipos-1,kp)
c           if(js.eq.0) then
c           itl=ipos+1
c           else
c           itl=ilq(ihm(js),1)+1
            if(js.gt.0) itl=ilq(ihm(js),1)+1
c           endif
          endif
          do it=itl,ilq(ihm(is),1)
            ip=ieq(it,kp)
            if(ip.ne.0) then
              do jg=nuc(ip-1)+1,nuc(ip)
                xx=corm(jg,1)-rax
                yy=corm(jg,2)-ray
                zz=corm(jg,3)-raz
                if(isur.eq.0.or.(kp.eq.istr.and.it.eq.itl)) then
                  corm(jg,1)=(rx*rx+(1-rx*rx)*cwd)*xx+(rx*ry*(1-cwd)-rz*swd)*yy+
     1                       (rx*rz*(1-cwd)+ry*swd)*zz+drx
                  corm(jg,2)=(rx*ry*(1-cwd)+rz*swd)*xx+(ry*ry+(1-ry*ry)*cwd)*yy+
     1                       (ry*rz*(1-cwd)-rx*swd)*zz+dry
                  corm(jg,3)=(rx*rz*(1-cwd)-ry*swd)*xx+(ry*rz*(1-cwd)+rx*swd)*yy+
     1                       (rz*rz+(1-rz*rz)*cwd)*zz+drz
                endif
                if(isur.ne.0) then  ! superhelical symmetry
                  xx=corm(jg,1)-rad
                  zz=corm(jg,3)
                  corm(jg,1)= cw*xx+sw*zz+rad
                  corm(jg,2)=corm(jg,2)+zed
                  corm(jg,3)=-sw*xx+cw*zz
                endif
              enddo
c---------------------------------------------------------update helical vectors
              do l=1,3
                ls=(l-1)*3+1
                xx=ha(ip,ls)
                yy=ha(ip,ls+1)
                zz=ha(ip,ls+2)
                if(isur.eq.0.or.(kp.eq.istr.and.it.eq.itl)) then
                  if(l.eq.3) then
                    xx=xx-rax
                    yy=yy-ray
                    zz=zz-raz
                  endif
                  ha(ip,ls)  =(rx*rx+(1-rx*rx)*cwd)*xx+(rx*ry*(1-cwd)-rz*swd)*yy+
     1                        (rx*rz*(1-cwd)+ry*swd)*zz
                  ha(ip,ls+1)=(rx*ry*(1-cwd)+rz*swd)*xx+(ry*ry+(1-ry*ry)*cwd)*yy+
     1                        (ry*rz*(1-cwd)-rx*swd)*zz
                  ha(ip,ls+2)=(rx*rz*(1-cwd)-ry*swd)*xx+(ry*rz*(1-cwd)+rx*swd)*yy+
     1                        (rz*rz+(1-rz*rz)*cwd)*zz
                  if(l.eq.3) then
                    ha(ip,ls)  =ha(ip,ls)+drx
                    ha(ip,ls+1)=ha(ip,ls+1)+dry
                    ha(ip,ls+2)=ha(ip,ls+2)+drz
                  endif
                endif
                if(isur.ne.0) then  ! superhelical symmetry
                  xx=ha(ip,ls)
                  zz=ha(ip,ls+2)
                  if(l.eq.3) xx=xx-rad
                  ha(ip,ls)  = cw*xx+sw*zz
                  ha(ip,ls+2)=-sw*xx+cw*zz
                  if(l.eq.3) then
                    ha(ip,ls)=ha(ip,ls)+rad
                    ha(ip,ls+1)=ha(ip,ls+1)+zed
                  endif
                endif
              enddo
c-----------------------------------------------------------------update ligands
              if(ip.le.kseq.and.nlig.ne.0) then
                do il=1,nlig
                  if(ip.eq.ilig(il,1)) then
                    ih=nto+il
                    do jg=nuc(ih-1)+1,nuc(ih)
                      xx=corm(jg,1)-rax
                      yy=corm(jg,2)-ray
                      zz=corm(jg,3)-raz
                      if(isur.eq.0.or.(kp.eq.istr.and.it.eq.itl)) then
                        corm(jg,1)=(rx*rx+(1-rx*rx)*cwd)*xx+(rx*ry*(1-cwd)-rz*swd)*yy+
     1                             (rx*rz*(1-cwd)+ry*swd)*zz+drx
                        corm(jg,2)=(rx*ry*(1-cwd)+rz*swd)*xx+(ry*ry+(1-ry*ry)*cwd)*yy+
     1                             (ry*rz*(1-cwd)-rx*swd)*zz+dry
                        corm(jg,3)=(rx*rz*(1-cwd)-ry*swd)*xx+(ry*rz*(1-cwd)+rx*swd)*yy+
     1                             (rz*rz+(1-rz*rz)*cwd)*zz+drz
                      endif
                      if(isur.ne.0) then
                        xx=corm(jg,1)-rad
                        zz=corm(jg,3)
                        corm(jg,1)=cw*xx+sw*zz+rad
                        corm(jg,2)=corm(jg,2)+zed
                        corm(jg,3)=-sw*xx+cw*zz
                      endif
                    enddo
                  endif
                enddo
              endif
c-------------------------------------------------------------------------------
            endif
          enddo
        enddo
      enddo    ! the end of cycle over residues for the helical variables
c=================================================================kink variables
      do is=1,kseq
        if(kink(is)) then
          dsx=vkink(is,1)
          dsy=vkink(is,2)
          ktl=vkink(is,3)
          do l=1,4
            k=k+1
            dk(l)=0.
            if(.not.lock(k)) then
              ki=ki+1
              dk(l)=var(ki)-vkink(is,l)
              vkink(is,l)=var(ki)
            endif
          enddo
          ux=ua(is-1,1)
          uy=ua(is-1,2)
          uz=ua(is-1,3)
          rx=da(is-1,1)
          ry=da(is-1,2)
          rz=da(is-1,3)
          wa(1)=uy*rz-uz*ry
          wa(2)=uz*rx-ux*rz
          wa(3)=ux*ry-uy*rx
          ctl=cos(cdr*(dk(3)))
          stl=-sin(cdr*(dk(3)))
          cpr=cos(cdr*(dk(4)))
          spr=sin(cdr*(dk(4)))
          ctt=cos(cdr*(ktl+dk(3)))
          stt=-sin(cdr*(ktl+dk(3)))
          do l=1,3
            dol(l)=ra(is-1,l)+da(is-1,l)*dsx-wa(l)*dsy
            don(l)=dol(l)+da(is-1,l)*dk(1)-wa(l)*dk(2)
          enddo
          px=(rx*rx+(1-rx*rx)*ctt)*wa(1)+(rx*ry*(1-ctt)-rz*stt)*wa(2)+
     1       (rx*rz*(1-ctt)+ry*stt)*wa(3)
          py=(rx*ry*(1-ctt)+rz*stt)*wa(1)+(ry*ry+(1-ry*ry)*ctt)*wa(2)+
     1       (ry*rz*(1-ctt)-rx*stt)*wa(3)
          pz=(rx*rz*(1-ctt)-ry*stt)*wa(1)+(ry*rz*(1-ctt)+rx*stt)*wa(2)+
     1       (rz*rz+(1-rz*rz)*ctt)*wa(3)
          do kp=1,nst
            do it=is,kseq
              ip=ieq(it,kp)
              if(ip.ne.0) then
                do jg=nuc(ip-1)+1,nuc(ip)
                  xx=corm(jg,1)-dol(1)
                  yy=corm(jg,2)-dol(2)
                  zz=corm(jg,3)-dol(3)
                  xx1=(rx*rx+(1-rx*rx)*ctl)*xx+(rx*ry*(1-ctl)-rz*stl)*yy+
     1                (rx*rz*(1-ctl)+ry*stl)*zz
                  yy1=(rx*ry*(1-ctl)+rz*stl)*xx+(ry*ry+(1-ry*ry)*ctl)*yy+
     1                (ry*rz*(1-ctl)-rx*stl)*zz
                  zz1=(rx*rz*(1-ctl)-ry*stl)*xx+(ry*rz*(1-ctl)+rx*stl)*yy+
     1                (rz*rz+(1-rz*rz)*ctl)*zz
                  corm(jg,1)=(px*px+(1-px*px)*cpr)*xx1+(px*py*(1-cpr)-pz*spr)*yy1+
     1                (px*pz*(1-cpr)+py*spr)*zz1+don(1)
                  corm(jg,2)=(px*py*(1-cpr)+pz*spr)*xx1+(py*py+(1-py*py)*cpr)*yy1+
     1                (py*pz*(1-cpr)-px*spr)*zz1+don(2)
                  corm(jg,3)=(px*pz*(1-cpr)-py*spr)*xx1+(py*pz*(1-cpr)+px*spr)*yy1+
     1                (pz*pz+(1-pz*pz)*cpr)*zz1+don(3)
                enddo
c---------------------------------------------------------update helical vectors
                do l=1,3
                  ls=(l-1)*3+1
                  xx=ha(ip,ls)
                  yy=ha(ip,ls+1)
                  zz=ha(ip,ls+2)
                  if(l.eq.3) then
                    xx=xx-dol(1)
                    yy=yy-dol(2)
                    zz=zz-dol(3)
                  endif
                  xx1=(rx*rx+(1-rx*rx)*ctl)*xx+(rx*ry*(1-ctl)-rz*stl)*yy+
     1                (rx*rz*(1-ctl)+ry*stl)*zz
                  yy1=(rx*ry*(1-ctl)+rz*stl)*xx+(ry*ry+(1-ry*ry)*ctl)*yy+
     1                (ry*rz*(1-ctl)-rx*stl)*zz
                  zz1=(rx*rz*(1-ctl)-ry*stl)*xx+(ry*rz*(1-ctl)+rx*stl)*yy+
     1                (rz*rz+(1-rz*rz)*ctl)*zz
                  ha(ip,ls)  =(px*px+(1-px*px)*cpr)*xx1+(px*py*(1-cpr)-pz*spr)*yy1+
     1                        (px*pz*(1-cpr)+py*spr)*zz1
                  ha(ip,ls+1)=(px*py*(1-cpr)+pz*spr)*xx1+(py*py+(1-py*py)*cpr)*yy1+
     1                        (py*pz*(1-cpr)-px*spr)*zz1
                  ha(ip,ls+2)=(px*pz*(1-cpr)-py*spr)*xx1+(py*pz*(1-cpr)+px*spr)*yy1+
     1                        (pz*pz+(1-pz*pz)*cpr)*zz1
                  if(l.eq.3) then
                    ha(ip,ls)  =ha(ip,ls)+don(1)
                    ha(ip,ls+1)=ha(ip,ls+1)+don(2)
                    ha(ip,ls+2)=ha(ip,ls+2)+don(3)
                  endif
                enddo
c-----------------------------------------------------------------update ligands
                if(ip.le.kseq.and.nlig.ne.0) then
                  do il=1,nlig
                    if(ip.eq.ilig(il,1)) then
                      ih=nto+il
                      do jg=nuc(ih-1)+1,nuc(ih)
                        xx=corm(jg,1)-dol(1)
                        yy=corm(jg,2)-dol(2)
                        zz=corm(jg,3)-dol(3)
                        xx1=(rx*rx+(1-rx*rx)*ctl)*xx+(rx*ry*(1-ctl)-rz*stl)*yy+
     1                      (rx*rz*(1-ctl)+ry*stl)*zz
                        yy1=(rx*ry*(1-ctl)+rz*stl)*xx+(ry*ry+(1-ry*ry)*ctl)*yy+
     1                      (ry*rz*(1-ctl)-rx*stl)*zz
                        zz1=(rx*rz*(1-ctl)-ry*stl)*xx+(ry*rz*(1-ctl)+rx*stl)*yy+
     1                      (rz*rz+(1-rz*rz)*ctl)*zz
                        corm(jg,1)=(px*px+(1-px*px)*cpr)*xx1+(px*py*(1-cpr)-pz*spr)*yy1+
     1                             (px*pz*(1-cpr)+py*spr)*zz1+don(1)
                        corm(jg,2)=(px*py*(1-cpr)+pz*spr)*xx1+(py*py+(1-py*py)*cpr)*yy1+
     1                             (py*pz*(1-cpr)-px*spr)*zz1+don(2)
                        corm(jg,3)=(px*pz*(1-cpr)-py*spr)*xx1+(py*pz*(1-cpr)+px*spr)*yy1+
     1                             (pz*pz+(1-pz*pz)*cpr)*zz1+don(3)
                      enddo
                    endif
                  enddo
                endif
c-------------------------------------------------------------------------------
              endif
            enddo
          enddo
        endif
      enddo    ! end of update of kink variables 
c=====================================================================superhelix
      if(isur.ne.0) then
        rado=rad
        if(isur.gt.0) rad=var(isur)
        pito=pit
        if(isup.gt.0) pit=var(isup)
        po=atan(pito/(pi2*rado))
        pn=atan(pit/(pi2*rad))
        cp=cos(pn-po)
        sp=sin(pn-po)
        roto=sqrt((2*pi*rado)**2+pito**2)
        root=sqrt((2*pi*rad)**2+pit**2)
        do ip=1,nto
          dox=ra(ip,1)
          doy=ra(ip,2)
          doz=ra(ip,3)
          th1=acos((rado-dox)/rado)
          if(doz.lt.0) th1=-th1
          th2=th1*roto/root
          cr=cos(th2-th1)
          sr=sin(th2-th1)
          ra(ip,1)=rad*(1-cos(th2))
          ra(ip,2)=th2*pit/(pi2)
          ra(ip,3)=rad*sin(th2)
          rx=(dox-rado)/rado
          ry=0.
          rz=doz/rado
          do jg=nuc(ip-1)+1,nuc(ip)
            xx=corm(jg,1)-dox
            yy=corm(jg,2)-doy
            zz=corm(jg,3)-doz
            xx1=(rx*rx+(1-rx*rx)*cp)*xx+(rx*ry*(1-cp)-rz*sp)*yy+
     1          (rx*rz*(1-cp)+ry*sp)*zz
            yy1=(rx*ry*(1-cp)+rz*sp)*xx+(ry*ry+(1-ry*ry)*cp)*yy+
     1          (ry*rz*(1-cp)-rx*sp)*zz
            zz1=(rx*rz*(1-cp)-ry*sp)*xx+(ry*rz*(1-cp)+rx*sp)*yy+
     1          (rz*rz+(1-rz*rz)*cp)*zz
            corm(jg,1)=cr*xx1+sr*zz1+ra(ip,1)
            corm(jg,2)=yy1+ra(ip,2)
            corm(jg,3)=-sr*xx1+cr*zz1+ra(ip,3)
          enddo
c---------------------------------------------------------update helical vectors
          do l=1,2
            ls=(l-1)*3+1
            xx=ha(ip,ls)
            yy=ha(ip,ls+1)
            zz=ha(ip,ls+2)
            xx1=(rx*rx+(1-rx*rx)*cp)*xx+(rx*ry*(1-cp)-rz*sp)*yy+
     1          (rx*rz*(1-cp)+ry*sp)*zz
            yy1=(rx*ry*(1-cp)+rz*sp)*xx+(ry*ry+(1-ry*ry)*cp)*yy+
     1          (ry*rz*(1-cp)-rx*sp)*zz
            zz1=(rx*rz*(1-cp)-ry*sp)*xx+(ry*rz*(1-cp)+rx*sp)*yy+
     1          (rz*rz+(1-rz*rz)*cp)*zz
            ha(ip,ls)  =cr*xx1+sr*zz1
            ha(ip,ls+1)=yy1
            ha(ip,ls+2)=-sr*xx1+cr*zz1
          enddo
c-----------------------------------------------------------------update ligands
          if(ip.le.kseq.and.nlig.ne.0) then
            do il=1,nlig
              if(ip.eq.ilig(il,1)) then
                ih=nto+il
                do jg=nuc(ih-1)+1,nuc(ih)
                  xx=corm(jg,1)-dox
                  yy=corm(jg,2)-doy
                  zz=corm(jg,3)-doz
                  xx1=(rx*rx+(1-rx*rx)*cp)*xx+(rx*ry*(1-cp)-rz*sp)*yy+
     1                (rx*rz*(1-cp)+ry*sp)*zz
                  yy1=(rx*ry*(1-cp)+rz*sp)*xx+(ry*ry+(1-ry*ry)*cp)*yy+
     1                (ry*rz*(1-cp)-rx*sp)*zz
                  zz1=(rx*rz*(1-cp)-ry*sp)*xx+(ry*rz*(1-cp)+rx*sp)*yy+
     1                (rz*rz+(1-rz*rz)*cp)*zz
                  corm(jg,1)=cr*xx1+sr*zz1+ra(ip,1)
                  corm(jg,2)=yy1+ra(ip,2)
                  corm(jg,3)=-sr*xx1+cr*zz1+ra(ip,3)
                enddo
              endif
            enddo
          endif
        enddo
      endif
c=============================================================reposition ligands
      do il=1,nlig
        is=nto+il
        kal=nuc(is)-nuc(is-1)
        i1=ilig(il,1)
        lint=6
        if(kal.eq.1) lint=3
        do j=1,lint
          k=k+1
          if(.not.lock(k)) then
            ki=ki+1
            dlig(j)=var(ki)-rlig(il,j)
            rlig(il,j)=var(ki)
          else
            dlig(j)=0.
          endif
        enddo
        if(i1.eq.0) then
          rx=0.
          ry=0.
          rz=1.
          dx=-cos(cdr*slig(il,2))
          dy=-sin(cdr*slig(il,2))
          dz=0.
        else
          rx=ua(i1,1)
          ry=ua(i1,2)
          rz=ua(i1,3)
          dx=-da(i1,1)
          dy=-da(i1,2)
          dz=-da(i1,3)
        endif
        wx=ry*dz-dy*rz
        wy=rz*dx-dz*rx
        wz=rx*dy-dx*ry
        delx=dlig(1)*dx+dlig(2)*wx+dlig(3)*rx
        dely=dlig(1)*dy+dlig(2)*wy+dlig(3)*ry
        delz=dlig(1)*dz+dlig(2)*wz+dlig(3)*rz
        if(kal.eq.1) then
          do i=nuc(is-1)+1,nuc(is)
            corm(i,1)=corm(i,1)+delx
            corm(i,2)=corm(i,2)+dely
            corm(i,3)=corm(i,3)+delz
          enddo
        else
          ipv=nuc(is-1)+lpiv(il)
          x0=corm(ipv,1)
          y0=corm(ipv,2)
          z0=corm(ipv,3)
          cx=cos(cdr*(dlig(4)))
          sx=sin(cdr*(dlig(4)))
          cxt=cos(cdr*(rlig(il,4)))
          sxt=sin(cdr*(rlig(il,4)))
          cy=cos(cdr*(dlig(5)))
          sy=sin(cdr*(dlig(5)))
          cz=cos(cdr*(dlig(6)))
          sz=sin(cdr*(dlig(6)))
          czt=cos(cdr*(rlig(il,6)))
          szt=sin(cdr*(rlig(il,6)))
          dxt=(rx*rx+(1-rx*rx)*czt)*dx+(rx*ry*(1-czt)-rz*szt)*dy+
     1        (rx*rz*(1-czt)+ry*szt)*dz
          dyt=(rx*ry*(1-czt)+rz*szt)*dx+(ry*ry+(1-ry*ry)*czt)*dy+
     1        (ry*rz*(1-czt)-rx*szt)*dz
          dzt=(rx*rz*(1-czt)-ry*szt)*dx+(ry*rz*(1-czt)+rx*szt)*dy+
     1        (rz*rz+(1-rz*rz)*czt)*dz
          wx=ry*dzt-dyt*rz
          wy=rz*dxt-dzt*rx
          wz=rx*dyt-dxt*ry
          rx=dxt
          ry=dyt
          rz=dzt
          wxt=(rx*rx+(1-rx*rx)*cxt)*wx+(rx*ry*(1-cxt)-rz*sxt)*wy+
     1        (rx*rz*(1-cxt)+ry*sxt)*wz
          wyt=(rx*ry*(1-cxt)+rz*sxt)*wx+(ry*ry+(1-ry*ry)*cxt)*wy+
     1        (ry*rz*(1-cxt)-rx*sxt)*wz
          wzt=(rx*rz*(1-cxt)-ry*sxt)*wx+(ry*rz*(1-cxt)+rx*sxt)*wy+
     1        (rz*rz+(1-rz*rz)*cxt)*wz
          do i=nuc(is-1)+1,nuc(is)
            xx=corm(i,1)-x0
            yy=corm(i,2)-y0
            zz=corm(i,3)-z0
            if(i1.eq.0) then
              rx=0.
              ry=0.
              rz=1.
            else
              rx=ua(i1,1)
              ry=ua(i1,2)
              rz=ua(i1,3)
            endif
            xx1=(rx*rx+(1-rx*rx)*cz)*xx+(rx*ry*(1-cz)-rz*sz)*yy+
     1          (rx*rz*(1-cz)+ry*sz)*zz
            yy1=(rx*ry*(1-cz)+rz*sz)*xx+(ry*ry+(1-ry*ry)*cz)*yy+
     1          (ry*rz*(1-cz)-rx*sz)*zz
            zz1=(rx*rz*(1-cz)-ry*sz)*xx+(ry*rz*(1-cz)+rx*sz)*yy+
     1          (rz*rz+(1-rz*rz)*cz)*zz
            rx=dxt
            ry=dyt
            rz=dzt
            xx2=(rx*rx+(1-rx*rx)*cx)*xx1+(rx*ry*(1-cx)-rz*sx)*yy1+
     1          (rx*rz*(1-cx)+ry*sx)*zz1
            yy2=(rx*ry*(1-cx)+rz*sx)*xx1+(ry*ry+(1-ry*ry)*cx)*yy1+
     1          (ry*rz*(1-cx)-rx*sx)*zz1
            zz2=(rx*rz*(1-cx)-ry*sx)*xx1+(ry*rz*(1-cx)+rx*sx)*yy1+
     1          (rz*rz+(1-rz*rz)*cx)*zz1
            rx=wxt
            ry=wyt
            rz=wzt
            corm(i,1)=(rx*rx+(1-rx*rx)*cy)*xx2+(rx*ry*(1-cy)-rz*sy)*yy2+
     1                (rx*rz*(1-cy)+ry*sy)*zz2+x0+delx
            corm(i,2)=(rx*ry*(1-cy)+rz*sy)*xx2+(ry*ry+(1-ry*ry)*cy)*yy2+
     1                (ry*rz*(1-cy)-rx*sy)*zz2+y0+dely
            corm(i,3)=(rx*rz*(1-cy)-ry*sy)*xx2+(ry*rz*(1-cy)+rx*sy)*yy2+
     1                (rz*rz+(1-rz*rz)*cy)*zz2+z0+delz
          enddo
        endif
      enddo
c=================================================================thymine methyl
      kth=ntba
      do is=1,nto
        ks=(ilq(is,2)-1)*kseq+ilq(is,1)
        if(lthy(is)) then
          in=itr(is)
          ioff=nuc(is-1)
          k=k+1
          kth=kth+1
          if(.not.lock(k)) then
            ki=ki+1
            del=var(ki)-sap(kth)
            sap(kth)=var(ki)
            ib=ithy(2)+iofs(is)+ioff
            ic=ithy(3)+iofs(is)+ioff
            x0=corm(ib,1)
            y0=corm(ib,2)
            z0=corm(ib,3)
            rx=corm(ic,1)-x0
            ry=corm(ic,2)-y0
            rz=corm(ic,3)-z0
            rr=sqrt(rx*rx+ry*ry+rz*rz)
            rx=rx/rr
            ry=ry/rr
            rz=rz/rr
            ca=cos(cdr*(del))
            sa=sin(cdr*(del))
            do j=4,6
              jg=ithy(j)+iofs(is)+ioff
              xx=corm(jg,1)-x0
              yy=corm(jg,2)-y0
              zz=corm(jg,3)-z0
              corm(jg,1)=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     1                   (rx*rz*(1-ca)+ry*sa)*zz+x0
              corm(jg,2)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1                   (ry*rz*(1-ca)-rx*sa)*zz+y0
              corm(jg,3)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1                   (rz*rz+(1-rz*rz)*ca)*zz+z0
            enddo
          endif
        endif
      enddo
c=========================================================junction dep variables
      k=0
      do is=1,nto
        in=itr(is)
        ino=ito(is)
        ioff=nuc(is-1)
        idir=idr(ilq(is,2))
        nsv=nsr(is)*5
        idel=kap(1,in)-kap(1,ino)+(nsr(is)-1)*5
        k=k+kap(1,in)+nsv
        do l=kap(1,in)+nsv+1,kap(2,in)
          k=k+1
          if(nap(l,6,in).eq.0) then
            ia=nap(l,1,in)+ioff
            ib=nap(l,2,in)+ioff
            ic=nap(l,3,in)+ioff
            id=nap(l,4,in)+ioff
            ix=nap(l,5,in)
            if(id.eq.ioff) then
              if(idir.gt.0) then
                if(ix.ge.1) ic=ic-ioff+nuc(is)
                if(ix.eq.2) ib=ib-ioff+nuc(is)
                if(ix.le.-1) ic=ic-ioff+nuc(is-2)
                if(ix.eq.-2) ib=ib-ioff+nuc(is-2)
              else
                if(ix.ge.1) ic=ic-ioff+nuc(is-2)
                if(ix.eq.2) ib=ib-ioff+nuc(is-2)
                if(ix.le.-1) ic=ic-ioff+nuc(is)
                if(ix.eq.-2) ib=ib-ioff+nuc(is)
              endif
              del=ang(ia,ib,ic)-sap(k)
            else 
              if(idir.gt.0) then
                if(ix.ge.1) id=id-ioff+nuc(is)
                if(ix.eq.2) ic=ic-ioff+nuc(is)
                if(ix.le.-1) id=id-ioff+nuc(is-2)
                if(ix.eq.-2) ic=ic-ioff+nuc(is-2)
              else
                if(ix.ge.1) id=id-ioff+nuc(is-2)
                if(ix.eq.2) ic=ic-ioff+nuc(is-2)
                if(ix.le.-1) id=id-ioff+nuc(is)
                if(ix.eq.-2) ic=ic-ioff+nuc(is)
              endif
              del=torp(ia,ib,ic,id)-sap(k)
            endif
            sap(k)=sap(k)+del
            do ld=kap(1,in)+nsv+1,kap(2,in) ! update variable that change synchroneously with given (l)
              ip=-nap(ld,6,in)
              if(ip.eq.l) sap(k-l+ld)=sap(k-l+ld)+del
            enddo
c------------------------------------------------------------reset c5' hydrogens
            if(nick.eq.is) goto 400
            lm=l-idel
            if((ino.eq.2.and.lm.eq.16).or.(ino.eq.3.and.lm.eq.12).or.
     1         (ino.eq.5.and.lm.eq.17).or.(ino.eq.6.and.lm.eq.13)) then
c torsion: 3-4-6-16 (C3'-C4'- C5'-O5'(prev res in strand)) for ino ==2, ==3,==5, ==6 ( normal and 3' end residue)
c rotate H5' and H5" to set to zero angle between average O5'-C5'-C4' - H5' and  O5'-C5'-C4' - H5"
c
               x0=corm(ib,1)
               y0=corm(ib,2)
               z0=corm(ib,3)
               rx=corm(ic,1)-x0
               ry=corm(ic,2)-y0
               rz=corm(ic,3)-z0
               rr=sqrt(rx*rx+ry*ry+rz*rz)
               rx=rx/rr
               ry=ry/rr
               rz=rz/rr
               t13=torp(id,ic,ib,ioff+13) !  torsion with H5'(13th atom) and 
               t14=torp(id,ic,ib,ioff+14) !  H5"(14th atom in the residue)
               den=(t13+t14)/2
               if(abs(t13-den).lt.90) den=den+180
               do j=1,nuc(is)-nuc(is-1)
                 m=iap(l,j,in)  ! this will be non-zero only for H5' and H5"
                 if(m.ne.0) then
                   jg=j+ioff
                   ca=cos(cdr*(den/m))
                   sa=sin(cdr*(den/m))
                   xx=corm(jg,1)-x0
                   yy=corm(jg,2)-y0
                   zz=corm(jg,3)-z0
                   corm(jg,1)=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     1                        (rx*rz*(1-ca)+ry*sa)*zz+x0
                   corm(jg,2)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1                        (ry*rz*(1-ca)-rx*sa)*zz+y0
                   corm(jg,3)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1                        (rz*rz+(1-rz*rz)*ca)*zz+z0
                 endif
               enddo
            endif
            if((ino.eq.2.and.lm.eq.17).or.(ino.eq.3.and.lm.eq.13).or.
     1         (ino.eq.5.and.lm.eq.18).or.(ino.eq.6.and.lm.eq.14)) then
c   valence angle 4-6-16 ( C4' - C5' - O5'(prev res in chain)) 
c    move H5' and H5" dur to changes of this angle 
c
              x0=corm(ib,1)
              y0=corm(ib,2)
              z0=corm(ib,3)
              ax=x0-corm(ia,1)
              ay=y0-corm(ia,2)
              az=z0-corm(ia,3)
              bx=corm(ic,1)-x0
              by=corm(ic,2)-y0
              bz=corm(ic,3)-z0
              rx=by*az-bz*ay
              ry=bz*ax-bx*az
              rz=bx*ay-by*ax
              rr=sqrt(rx*rx+ry*ry+rz*rz)
              rx=rx/rr
              ry=ry/rr
              rz=rz/rr
              do j=1,nuc(is)-nuc(is-1)
                m=iap(l,j,in)
                if(m.ne.0) then
                  jg=j+ioff
                  ca=cos(cdr*(del/m))
                  sa=sin(cdr*(del/m))
                  xx=corm(jg,1)-x0
                  yy=corm(jg,2)-y0
                  zz=corm(jg,3)-z0
                  corm(jg,1)=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     1                       (rx*rz*(1-ca)+ry*sa)*zz+x0
                  corm(jg,2)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1                       (ry*rz*(1-ca)-rx*sa)*zz+y0
                  corm(jg,3)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1                       (rz*rz+(1-rz*rz)*ca)*zz+z0
                endif
              enddo
            endif
400       endif
        enddo
      enddo
      return
      end
