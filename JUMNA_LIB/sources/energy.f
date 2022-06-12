      subroutine energy
      include 'jumna_data.inc'
      logical*2 ifhb,kink,lthy,hst,bst,vst,lgi,lgj,lgu,ribose,
     1 lock,sup,rcom,homo,homo2,homo3,diep,link,ecen,cyl,lcat,
     1 sum,cent,autos,locr,lar,cation,amber
      character*4 mnam,munit,seq*120,code*8,kode*8,lnam
      integer*2 i23,i34,elim
      integer*4 opt,fold,lex(29)
      common/cnd/rcyl,dcyl,hcyl,tcyl,eneq,ucyl(3),pcyl(3),ecyl,
     1 ac(25),bc(25),lcyl(3)
      common/datjm/acc,phos,epsi,epsr,plat,slope,rhbl,vfac,tfac,rfac,xfac,
     1 scale,rpiv,tpiv,fad,fan,damp,fnoew,fnoes,fnoea,df1,scnb,scee,
     1 catd,catr,catc,rad,pit,enit,nshel,maxn,opt,mhomo,mhomo2,mhomo3,
     1 nick,limit,nop,nrib,ncat,lig,naxo,sup,rcom,homo,homo2,homo3,diep,
     1 sum,link,ecen,cyl,lcat,cent,autos,amber
      common/enf/for(n1,3),tor(n1,3),fot(n6),ener,elec,
     1 repl,disp,eang,etog,epen,i23(8*n1),i34(8*n1),elim(n1)
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/flx/sap(n6),refg,refb,refh,refv,iap(n8,n4,n5),
     1 nap(n8,7,n5),kap(3,n5)
      common/lgd/rlig(n9,6),slig(n9,2),ilig(n9,6),lopt(n9),lpiv(n9),
     1 ntlg,nlgi,nlig,ntl,lnam(n9),locr(n2,n0)
      common/mnn/var(n7),gra(n7),con(n0*n2+n3),scl(n7),ncon,nvrc,
     1 ntba,nbac,nthe,nhel,ntki,nkin,ntri,ntot,nvar,nrin,
     1 lar(n7),lock(n6a)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/parjm/aij(25,25),bij(25,25),ahb(15),bhb(15),vt(35),va(25),
     1 vo(25),delmax,time0,nzsh,nwdg,ivf(35),ifhb(25,25)
      common/slg/hst(n2,6),bst(n2,n8),vst(n2,4),lgi(n1),lgj(n1),lgu(n2)
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(120,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      data lex/15*0,1,1,1,-1,-1,0,0,0,0,0,0,-1,1,1/
      do j=1,3
        do i=1,kam
          for(i,j)=0.
          tor(i,j)=0.
        enddo
      enddo
      repl=0.
      disp=0.
      elec=0.
      eang=0.
      etog=0.
c      if(etog.eq.0) goto 600
      k=0
      do i=1,kam-1
        if(.not.lgj(i)) goto 100
        k=k+1
        m=imty(i)
        dmi=dmon(i)
        ini=elim(i)
        do j=i+1,kam
          if(.not.lgj(j).or.(.not.lgi(i).and..not.lgi(j))) goto 200
          if(ini.ne.0.and.elim(j).eq.ini) goto 200
          do l=k+1,k+i23(k)
             if(j.eq.i23(l)) goto 200
          enddo
          n=imty(j)
          dmj=dmon(j)
c-----------------------------------------------------------------electrostatics
          dx=corm(j,1)-corm(i,1)
          dy=corm(j,2)-corm(i,2)
          dz=corm(j,3)-corm(i,3)
          r2=dx*dx+dy*dy+dz*dz
          r=sqrt(r2)
          dx=dx/r2
          dy=dy/r2
          dz=dz/r2
          if(diep) then
            alr=slope*r
            ex=exp(-alr)
            alr2=alr*alr
            e=plat-(plat-enit)*(alr2+2.*alr+2.)*ex/2.
            g=(plat-enit)*slope*alr2*ex/(2.*e*e)
            re=r*e
            sgf=convk*dmi*dmj*g
          else if(epsr.ne.0.0) then
            re=epsr*r2
            srij=convk*dmi*dmj*2/re
          else
            re=r*epsi
          endif
          sij=convk*dmi*dmj/re
          elec=elec+sij
c	    write(6,'("cur elec =", F12.6)')elec
          if(diep) then
             ss=sij+sgf
             for(i,1)=for(i,1)-ss*dx
             for(i,2)=for(i,2)-ss*dy
             for(i,3)=for(i,3)-ss*dz
             for(j,1)=for(j,1)+ss*dx
             for(j,2)=for(j,2)+ss*dy
             for(j,3)=for(j,3)+ss*dz
          else if(epsr.ne.0.0) then
             for(i,1)=for(i,1)-srij*dx
             for(i,2)=for(i,2)-srij*dy
             for(i,3)=for(i,3)-srij*dz
             for(j,1)=for(j,1)+srij*dx
             for(j,2)=for(j,2)+srij*dy
             for(j,3)=for(j,3)+srij*dz
          else
             for(i,1)=for(i,1)-sij*dx
             for(i,2)=for(i,2)-sij*dy
             for(i,3)=for(i,3)-sij*dz
             for(j,1)=for(j,1)+sij*dx
             for(j,2)=for(j,2)+sij*dy
             for(j,3)=for(j,3)+sij*dz
          endif
c------------------------------------------------------------------lennard-jones
          r6=r**6
          e6=aij(m,n)/r6
          e12=bij(m,n)/(r6*r6)
          ee=12.*e12-6.*e6
          if(ifhb(m,n).and.r.le.rhbl) then
            if(m.le.2) then
              ihy=i
              sid=1.
              iatype=n
            else
              ihy=j
              sid=-1.
              iatype=m
            endif
            mp=matd(ihy,1)
            idtype=imty(mp)
            if(iatype.eq.18) iatype=11
            if(idtype.eq.18) then
               ih=iatype-6
            else if(idtype.eq.8) then
               ih=iatype-1
            else
               ih=iatype+4
            endif
            qx=corm(ihy,1)-corm(mp,1)
            qy=corm(ihy,2)-corm(mp,2)
            qz=corm(ihy,3)-corm(mp,3)
            q=sqrt(qx*qx+qy*qy+qz*qz)
            cth=sid*(dx*qx+dy*qy+dz*qz)*r/q
            if(cth.lt.0.) then
               disp=disp-e6
               repl=repl+e12
               for(i,1)=for(i,1)-ee*dx
               for(i,2)=for(i,2)-ee*dy
               for(i,3)=for(i,3)-ee*dz
               for(j,1)=for(j,1)+ee*dx
               for(j,2)=for(j,2)+ee*dy
               for(j,3)=for(j,3)+ee*dz
            else
               adel=(ahb(ih)-aij(m,n))/r6
               bdel=(bhb(ih)-bij(m,n))/(r6*r6)
               disp=disp-e6-cth*adel
               repl=repl+e12+cth*bdel
               fhb=(adel-bdel)*r/q
               fac=sid/r2
               px=qx*fac
               py=qy*fac
               pz=qz*fac
               tx=fhb*px-(7.*cth*adel-13.*cth*bdel-ee)*dx
               ty=fhb*py-(7.*cth*adel-13.*cth*bdel-ee)*dy
               tz=fhb*pz-(7.*cth*adel-13.*cth*bdel-ee)*dz
               for(i,1)=for(i,1)-tx
               for(i,2)=for(i,2)-ty
               for(i,3)=for(i,3)-tz
               for(j,1)=for(j,1)+tx
               for(j,2)=for(j,2)+ty
               for(j,3)=for(j,3)+tz
               tox=fhb*(dy*qz-dz*qy)
               toy=fhb*(dz*qx-dx*qz)
               toz=fhb*(dx*qy-dy*qx)
               if(ihy.eq.i) then
                  tor(i,1)=tor(i,1)-tox
                  tor(i,2)=tor(i,2)-toy
                  tor(i,3)=tor(i,3)-toz
               endif
               if(ihy.eq.j) then
                  tor(j,1)=tor(j,1)+tox
                  tor(j,2)=tor(j,2)+toy
                  tor(j,3)=tor(j,3)+toz
               endif
            endif
          else
c---------------------------------------------------------------no hydrogen bond
            disp=disp-e6
            repl=repl+e12
            for(i,1)=for(i,1)-ee*dx
            for(i,2)=for(i,2)-ee*dy
            for(i,3)=for(i,3)-ee*dz
            for(j,1)=for(j,1)+ee*dx
            for(j,2)=for(j,2)+ee*dy
            for(j,3)=for(j,3)+ee*dz
          endif
200     enddo
        k=k+i23(k)
100   enddo
c--------------------------------------------------torsional/angle strain energy
      k=0
      ks=0
      do is=1,ntl
      in=itr(is)
      if(.not.lgu(is)) then
      k=k+kap(2,in)
      goto 500
      endif
      do l=1,kap(2,in)
      k=k+1
      ix=nap(l,6,in)
      it=nap(l,7,in)
      if((is.eq.nick.and.lex(l).eq.1).or.
     1  (is+1.eq.nick.and.lex(l).eq.-1)) goto 400
      if(nap(l,4,in).eq.0) then
      fcon=va(it)
      dang=(sap(k)-vo(it))*cdr
      eang=eang+fcon*dang**2
      fot(k)=-2.*fcon*dang
      else
      fold=ivf(it)
      barr=vt(it)
      etog=etog+barr*(1.+sign(1,fold)*cos(cdr*(abs(fold)*sap(k))))/2.
      if(ix.ge.0) then
      fot(k)=barr*fold*sin(cdr*(abs(fold)*sap(k)))/2.
      else
      fot(ks-ix)=fot(ks-ix)+barr*fold*sin(cdr*(abs(fold)*sap(k)))/2.
      endif
      endif
400   enddo
500   ks=k
      enddo
c-----------------------------------------------------------------thymine methyl
      fold=ivf(nith)
      barr=vt(nith)
      do is=1,nto
      ks=(ilq(is,2)-1)*kseq+ilq(is,1)
      if(lgu(is).and.lthy(is)) then
      k=k+1
      et=barr*(1.+sign(1,fold)*cos(cdr*(abs(fold)*sap(k))))/2.
      fot(k)=barr*fold*sin(cdr*(abs(fold)*sap(k)))/2.
      etog=etog+et
      endif
      enddo
c-------------------------------------------------------------------------------
 600  call penalty
      if(cyl) call enecyl(cyl)
      ener=elec+repl+disp+eang+etog+epen
      return
      end
