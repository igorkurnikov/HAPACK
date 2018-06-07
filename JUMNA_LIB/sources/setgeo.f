      subroutine setgeo(axe,iabas)
c Set initial physical variable array
c 
      include 'jumna_data.inc'
      logical*2 lthy,lar,lock,kink,ribose,cation,hst,bst,vst,lgi,lgj,
     1 lgu,locr,logax,sup,rcom,homo,homo2,homo3,diep,link,sum,
     1 ecen,cyl,lcat,cent,autos,amber
      character*4 mnam,munit,seq*120,code*8,kode*8,axe*32,lnam
      integer*4 opt
      common/datjm/acc,phos,epsi,epsr,plat,slope,rhbl,vfac,tfac,rfac,xfac,
     1 scale,rpiv,tpiv,fad,fan,damp,fnoew,fnoes,fnoea,df1,scnb,scee,
     1 catd,catr,catc,rad,pit,enit,nshel,maxn,opt,mhomo,mhomo2,mhomo3,
     1 nick,limit,nop,nrib,ncat,lig,naxo,sup,rcom,homo,homo2,homo3,diep,
     1 sum,link,ecen,cyl,lcat,cent,autos,amber
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/flx/sap(n6),refg,refb,refh,refv,iap(n8,n4,n5),
     1 nap(n8,7,n5),kap(3,n5)
      common/ind/iend(0:4),ise(n2),ienb(n2,2),kapt(n2+1),nsph
      common/lgd/rlig(n9,6),slig(n9,2),ilig(n9,6),lopt(n9),lpiv(n9),
     1 ntlg,nlgi,nlig,ntl,lnam(n9),locr(n2,n0)
      common/mnn/var(n7),gra(n7),con(n0*n2+n3),scl(n7),ncon,nvrc,
     1 ntba,nbac,nthe,nhel,ntki,nkin,ntri,ntot,nvar,nrin,
     1 lar(n7),lock(n6a)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/slg/hst(n2,6),bst(n2,n8),vst(n2,4),lgi(n1),lgj(n1),lgu(n2)
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      logax=.false.
      if(axe.ne.' ') logax=.true.
c-------------------------------------------------------------loop over residues
      kapt(1)=0
      do is=2,ntl+1
        ip=is-1
        inp=itr(ip)
        kapt(is)=kapt(ip)+kap(2,inp)
      enddo
      if(kapt(ntl+1).gt.n6) then
        write(6,*) '  ---- n6 too small ----'
        stop
      endif
      k=0
      do is=1,ntl
        in=itr(is)
        ino=ito(is)
        ioff=nuc(is-1)
c--------------------------------------calculate dep ring variables before moves
        kt=k+kap(1,in)
        do l=kap(1,in)+1,kap(1,in)+nsr(is)*5
          kt=kt+1
          ia=nap(l,1,in)+ioff
          ib=nap(l,2,in)+ioff
          ic=nap(l,3,in)+ioff
          id=nap(l,4,in)+ioff
          if(id.eq.ioff) then
            sap(kt)=ang(ia,ib,ic)
          else
            sap(kt)=torp(ia,ib,ic,id)
          endif
        enddo
c----------------------------------------------------------------indep variables
        do l=1,kap(1,in)
          k=k+1
          ia=nap(l,1,in)+ioff
          ib=nap(l,2,in)+ioff
          ic=nap(l,3,in)+ioff
          id=nap(l,4,in)+ioff
          if(id.eq.ioff) then
             sap(k)=ang(ia,ib,ic)
          else
             sap(k)=torp(ia,ib,ic,id)
          endif
c------------------------------------------------------------------set variables
          if(bst(is,l).or.logax) then
             if(set(is,l).lt.999.) then
*------------------------------------------Abasic Grenoble
                if(is.eq.iabas.and.l.eq.4) then
                   sapg=torp(21+ioff,ib,ic,id)
                   del=set(is,l)-sapg
                   sap(k)=sap(k)+del
*---------------------------------------------------------
                else
                   del=set(is,l)-sap(k)
                   sap(k)=set(is,l)
                endif
             else
                del=0.
             endif
             if(id.eq.ioff) then
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
     1                 (rx*rz*(1-ca)+ry*sa)*zz+x0
                corm(jg,2)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1                 (ry*rz*(1-ca)-rx*sa)*zz+y0
                corm(jg,3)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1                 (rz*rz+(1-rz*rz)*ca)*zz+z0
               endif
             enddo
          endif
        enddo
c-------------------------------------------------------------dep ring variables
        do l=kap(1,in)+1,kap(1,in)+nsr(is)*5
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
     1           (rx*rz*(1-ca)+ry*sa)*zz+x0
             corm(jg,2)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1           (ry*rz*(1-ca)-rx*sa)*zz+y0
             corm(jg,3)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1           (rz*rz+(1-rz*rz)*ca)*zz+z0
            endif
          enddo
        enddo
        k=kapt(is+1)
      enddo
c---------------------------------------------calculate remaining dep and gauche
      k=0
      do is=1,ntl
        in=itr(is)
        ino=ito(is)
        ioff=nuc(is-1)
        idir=0
        if(is.le.nto) idir=idr(ilq(is,2))
        nsv=nsr(is)*5
        k=k+kap(1,in)+nsv
        do l=kap(1,in)+nsv+1,kap(2,in)
          k=k+1
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
             sap(k)=ang(ia,ib,ic)
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
               sap(k)=torp(ia,ib,ic,id)
          endif
c----------------------------------------------------------correct c5' hydrogens
          if(is.le.nto) then
            idel=kap(1,in)-kap(1,ino)+(nsr(is)-1)*5
            if(nick.eq.is) goto 400
            lm=l-idel
            if((ino.eq.2.and.lm.eq.16).or.(ino.eq.3.and.lm.eq.12).or.
     1         (ino.eq.5.and.lm.eq.17).or.(ino.eq.6.and.lm.eq.13)) then
              del=sap(k)-refh
              rx=corm(ic,1)-corm(ib,1)
              ry=corm(ic,2)-corm(ib,2)
              rz=corm(ic,3)-corm(ib,3)
              rr=sqrt(rx*rx+ry*ry+rz*rz)
              rx=rx/rr
              ry=ry/rr
              rz=rz/rr
              x0=corm(ib,1)
              y0=corm(ib,2)
              z0=corm(ib,3)
              t13=torp(id,ic,ib,ioff+13)
              t14=torp(id,ic,ib,ioff+14)
              den=(t13+t14)/2
              if(abs(t13-den).lt.90) den=den+180
              do j=1,nuc(is)-nuc(is-1)
                m=iap(l,j,in)
                if(m.ne.0) then
                  jg=j+ioff
                  ca=cos(cdr*(den/m))
                  sa=sin(cdr*(den/m))
                  xx=corm(jg,1)-x0
                  yy=corm(jg,2)-y0
                  zz=corm(jg,3)-z0
                  corm(jg,1)=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     1                       (rx*rz*(1-ca)+ry*sa)*zz+x0
                  corm(jg,2)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1                      (ry*rz*(1-ca)-rx*sa)*zz+y0
                  corm(jg,3)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1                      (rz*rz+(1-rz*rz)*ca)*zz+z0
                endif
              enddo
            endif
            if((ino.eq.2.and.lm.eq.17).or.(ino.eq.3.and.lm.eq.13).or.
     1         (ino.eq.5.and.lm.eq.18).or.(ino.eq.6.and.lm.eq.14)) then
               del=sap(k)-refv
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
     1                        (rx*rz*(1-ca)+ry*sa)*zz + x0
                   corm(jg,2)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1                        (ry*rz*(1-ca)-rx*sa)*zz + y0
                   corm(jg,3)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1                        (rz*rz+(1-rz*rz)*ca)*zz + z0
                 endif
               enddo
            endif
400       endif
        enddo
      enddo
      return
      end
