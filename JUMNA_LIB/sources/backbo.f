      subroutine backbo
c
c List backbone variables and update C5' hydrogen position 
c
      include 'jumna_data.inc'
      logical*2 lar,lock,kink,lthy,ribose,cation,start,locr,
     1 ihl,ifhb,sup,rcom,homo,homo2,homo3,diep,link,
     1 sum,ecen,cyl,lcat,cent,autos,amber
      character*4 mnam,munit,sugt(10)*8,puck*8,seq*120,code*8,
     1 kode*8,lnam,snam,suni,sub,vin(n8)*1,vlo(n8)*1,vix(n8,4),
     1 rb*1,cat*1,knam,add*1
      dimension ts(5),tb(7),isg(5,6),isd(6),ibk(7,6),isg0(5,6),
     1 ibk0(7,6),ixv(n8),v1(3),v2(3),s1(3),s2(3),t(3),q1(3),q2(3)
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
      common/hel/ua(n2,3),da(n2,3),ra(n2,3)
      common/ind/iend(0:4),ise(n2),ienb(n2,2),kapt(n2+1),nsph
      common/lgd/rlig(n9,6),slig(n9,2),ilig(n9,6),lopt(n9),lpiv(n9),
     1 ntlg,nlgi,nlig,ntl,lnam(n9),locr(n2,n0)
      common/mnn/var(n7),gra(n7),con(n0*n2+n3),scl(n7),ncon,nvrc,
     1 ntba,nbac,nthe,nhel,ntki,nkin,ntri,ntot,nvar,nrin,
     1 lar(n7),lock(n6a)
      common/moljm/sor(n4,n5,3),smon(n4,n5),snam(n4,n5),suni(n4,n5),
     1 nuni(n4,n5),sub(n5),isch(n4,n5),isty(n4,n5),ics(n4,n5),
     1 mats(3*n4,n5),kas(n5),khs(n5),ksub
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/nmrjm/rnoe(n3),bnoe(n3),fnoe(n3),snoe,knam(200,2),ktyp(200),
     1 icon(n2*2),inoe(n3,9),ih68(n2),jnoe(n3),nnoe,ndcs,ndiv(0:99),nlin
      common/parjm/aij(25,25),bij(25,25),ahb(15),bhb(15),vt(35),va(25),
     1 vo(25),delmax,time0,nzsh,nwdg,ivf(35),ifhb(25,25)
      common/sst/phase(8,3),ampli(8,3),ph5(n2),am5(n2),kr5(n2)
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      common/symjm/xsum,ysum,itwl(n2),nvs(n7),ihm(n2),ksym(3),nbrk(2),
     1 isym,isur,isup,ihl(n2)
      data ((isg0(i,j),i=1,5),j=1,6)/ 5, 6,11,12,13,
     1                                5, 6,11,12,13,
     1                                5, 6, 7, 8, 9,
     1                                5, 6,12,13,14,
     1                                5, 6,12,13,14,
     1                                5, 6, 8, 9,10/
      data ((ibk0(i,j),i=1,7),j=1,6)/ 4, 0,11, 7, 8,16,17,
     1                                4,16,11, 7, 8,19,20,
     1                                4,12, 7, 0, 0, 0, 0,
     1                                4, 0,12, 7, 8,17,18,
     1                                4,17,12, 7, 8,20,21,
     1                                4,13, 8, 0, 0, 0, 0/
      data isd/0,0,0,11,11,7/
      data sugt/'C3''-endo','C4''-exo ','O1''-endo','C1''-exo ',
     1'C2''-endo','C3''-exo ','C4''-endo','O1''-exo ',
     1'C1''-endo','C2''-exo '/
c-------------------------------------------------------------------------------
      nt=ntba
      k=0
      ki=0
      if(nto.eq.0) goto 500
      write(6,14)
14      format(/2x,'Sugar   ',6x,'  C1-C2 ','  C2-C3 ','  Phase ',
     1    '  Ampli ',2x,'Pucker  ','  C2-O2 ','   Thy  '/)
c-------------------------------------------------------------loop over residues
      do is=1,nto
        ip=ilq(is,1)
        in=itr(is)
        ino=ito(is)
        nsv=nsr(is)*5
        idel=kap(1,in)-kap(1,ino)
        idels=idel+(nsr(is)-1)*5
        do j=1,6
          do i=1,5
            isg(i,j)=isg0(i,j)
            if(isg(i,j).gt.kap(1,ino)) isg(i,j)=isg(i,j)+idel
          enddo
        enddo
        ioff=nuc(is-1)
        idir=idr(ilq(is,2))
        io=ioff+iofs(is)
        do l=1,kap(2,in)
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
            del=ang(ia,ib,ic)-sap(k)
            sap(k)=sap(k)+del
          else
            if(l.gt.kap(1,in)+nsv) then
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
            endif
            del=torp(ia,ib,ic,id)-sap(k)
            sap(k)=sap(k)+del
            if(abs(sap(k)).gt.180.) sap(k)=sap(k)-sign(360.d0,sap(k))
          endif
c------------------------------------------------------------reset c5' hydrogens
          if(is.eq.nick) goto 400
          lm=l-idels
          if((ino.eq.2.and.lm.eq.16).or.(ino.eq.3.and.lm.eq.12).or.
     1       (ino.eq.5.and.lm.eq.17).or.(ino.eq.6.and.lm.eq.13)) then
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
     1                  (rx*rz*(1-ca)+ry*sa)*zz+x0
                  corm(jg,2)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1                  (ry*rz*(1-ca)-rx*sa)*zz+y0
                  corm(jg,3)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1                  (rz*rz+(1-rz*rz)*ca)*zz+z0
                endif
              enddo
          endif
          if((ino.eq.2.and.lm.eq.17).or.(ino.eq.3.and.lm.eq.13).or.
     1       (ino.eq.5.and.lm.eq.18).or.(ino.eq.6.and.lm.eq.14)) then
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
     1                  (rx*rz*(1-ca)+ry*sa)*zz+x0
                  corm(jg,2)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1                  (ry*rz*(1-ca)-rx*sa)*zz+y0
                  corm(jg,3)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1                  (rz*rz+(1-rz*rz)*ca)*zz+z0
                endif
              enddo
          endif
c----------------------------------------------------------------------reset var
400       if(l.le.kap(1,in).and..not.lock(k)) then
              ki=ki+1
              var(ki)=sap(k)
          endif
c----------------------------------------------------------------sugar variables
          do i=1,5
             if(isg(i,ino).eq.l) ts(i)=sap(k)
          enddo
          if(isd(ino).eq.l) rib=sap(k)
        enddo
c----------------------------------------------------------------pseudo-rotation
        a=0.
        b=0.
        do i=1,5
          j=i+1
          if(i.eq.5) j=1
          a=a+ts(j)*cos(cdr*(144.*(i-1)))
          b=b+ts(j)*sin(cdr*(144.*(i-1)))
        enddo
        a=a*2./5.
        b=-b*2./5.
        amp=sqrt(a*a+b*b)
        if(amp.gt.0.) then
          cp=a/amp
          sp=b/amp
          if(abs(cp).gt.1.) cp=sign(1.d0,cp)
            pha=acos(cp)*crd
            if(sp.lt.0.) pha=360.-pha
          else
            pha=0.
          endif
          puck=sugt(1+int(pha/36.))
          ph5(is)=pha
          am5(is)=amp
c-------------------------------------------------------------------------output
          if(ip.gt.isym) goto 190
          if(ise(is).lt.0.and.is.ne.1) write(6,8)
8         format(6x,
     1           '---------------------------------',
     1           '----------------------------------')
          cat=' '
          if(cation(is)) cat='*'
          if(lthy(is)) then
            nt=nt+1
            if(ribose(is)) then
               write(6,10) is,munit(io),nunit(io),'R',cat,ts(1),ts(2),pha,amp,
     1           puck,rib,sap(nt)
            else
               write(6,12) is,munit(io),nunit(io),' ',cat,ts(1),ts(2),pha,amp,
     1           puck,sap(nt)
            endif
          else
            if(ribose(is)) then
               write(6,10) is,munit(io),nunit(io),'R',cat,ts(1),ts(2),pha,amp,
     1                     puck,rib
            else
               write(6,10) is,munit(io),nunit(io),' ',cat,ts(1),ts(2),pha,amp,
     1                     puck
            endif
        endif
10      format(2x,i3,':',a4,i3,a1,a1,4f8.2,2x,a8,2f8.2)
12      format(2x,i3,':',a4,i3,a1,a1,4f8.2,2x,a8,8x,f8.2)
190   enddo
c-------------------------------------------------------------backbone variables
      write(6,20)
20    format(/2x,'Backbone',6x,'   Chi  ','  Gamma ','  Delta ',
     1              '  Epsil ','  Zeta  ','  Alpha ','  Beta  ',
     1              ' Eps-Zet',
     1      /16x,'  C1''-N ',' C5''-C4''',' C4''-C3''',' C3''-O3''',
     1            '  O3''-P ','  P-O5'' ',' O5''-C5''','  B1/B2 '/)
      k=0
      do is=1,nto
      ip=ilq(is,1)
      in=itr(is)
      ino=ito(is)
      idel=kap(1,in)-kap(1,ino)
      idels=idel+(nsr(is)-1)*5
      do j=1,6
      do i=1,7
      ib=ibk0(i,j)
      if(ib.gt.kap(1,ino).and.ib.le.kap(1,ino)+5) then
      ibk(i,j)=ib+idel
      else if(ib.gt.kap(1,ino)+5) then
      ibk(i,j)=ib+idels
      else
      ibk(i,j)=ib
      endif
      enddo
      enddo
      ioff=nuc(is-1)
      io=ioff+iofs(is)
      do l=1,kap(2,in)
      k=k+1
      do i=1,7
      if(ibk(i,ino).eq.l) tb(i)=sap(k)
      enddo
      tb(3)=torp(ioff+6,ioff+4,ioff+3,ioff+7)
      enddo
      rb=' '
      if(ribose(is)) rb='R'
      cat=' '
      if(cation(is)) cat='*'
      if(ise(is).lt.0.and.is.ne.1) write(6,28)
      if(ip.gt.isym+1) goto 195
      if(ino.eq.1.or.ino.eq.4) then
      write(6,31) is,munit(io),nunit(io),rb,cat,tb(1),(tb(i),i=3,7),
     1 tb(4)-tb(5)
      else if(ino.eq.2.or.ino.eq.5) then
      write(6,32) is,munit(io),nunit(io),rb,cat,(tb(i),i=1,7),
     1 tb(4)-tb(5)
      else
      write(6,33) is,munit(io),nunit(io),rb,cat,(tb(i),i=1,3)
      endif
28    format(6x,
     1'---------------------------------',
     1 '----------------------------------------')
31    format(2x,i3,':',a4,i3,a1,a1, f8.2,'  ......',6f8.2)
32    format(2x,i3,':',a4,i3,a1,a1,8f8.2)
33    format(2x,i3,':',a4,i3,a1,a1,3f8.2,'  ......',
     1 '  ......','  ......','  ......','  ......')
195   enddo
      call present
c----------------------------------------------------------constrained variables
500   start=.true.
      do k=1,nnoe
        i=ndcs+k
        j=jnoe(k)
        if(start.and.(j.ne.1.and.j.ne.3)) then
          write(6,87)
87        format(/2x,'Torsion/Valence/Sugar/Diff constraints ...'/)
          start=.false.
        endif
c-----------------------------------------------------------------------torsions
        if(j.eq.2) then
          i1=inoe(k,1)
          i2=inoe(k,2)
          i3=inoe(k,3)
          i4=inoe(k,4)
          dih=torp(i1,i2,i3,i4)
          delt=dih-rnoe(k)
          if(abs(delt).gt.180.) delt=delt-sign(360.d0,delt)
          write(6,81) dih,rnoe(k),delt,nunit(i1),mnam(i1),
     1      nunit(i2),mnam(i2),nunit(i3),mnam(i3),nunit(i4),mnam(i4)
81        format(2x,'T) Val: ',f7.2,' Req: ',6x,f7.2,6x,' Del: ',f7.2,
     1    i2,'/',a4,'-',i2,'/',a4,'-',i2,'/',a4,'-',i2,'/',a4)
        else if(j.eq.4) then
          i1=inoe(k,1)
          i2=inoe(k,2)
          i3=inoe(k,3)
          i4=inoe(k,4)
          dih=torp(i1,i2,i3,i4)
          delt=dih-rnoe(k)
          cx=cos(cdr*(bnoe(k)))
          sx=sin(cdr*(bnoe(k)))
          cm=cos(cdr*(rnoe(k)))
          sm=sin(cdr*(rnoe(k)))
          ct=cos(cdr*(dih))
          st=sin(cdr*(dih))
          v=cx*sm-sx*cm
          a=cm*st-sm*ct
          b=ct*sx-st*cx
          if((a.ge.0.and.b.ge.0).or.(v.gt.0.and.a*b.le.0)) then
             delt=0.
             goto 201
          endif
          if(a.gt.0.and.b.lt.0) then
             delt=dih-bnoe(k)
          else if(a.lt.0.and.b.lt.0) then
             a=ct*cm+st*sm
             b=cx*ct+sx*st
             if(a.lt.b) then
               delt=dih-bnoe(k)
             endif
          endif
201       if(abs(delt).gt.180.) delt=delt-sign(360.d0,delt)
             write(6,82) dih,rnoe(k),bnoe(k),delt,nunit(i1),mnam(i1),
     1         nunit(i2),mnam(i2),nunit(i3),mnam(i3),nunit(i4),mnam(i4)
82             format(2x,'t) Val: ',f7.2,' Req: ',f7.2,' <-> ',f7.2,
     1                      ' Del: ',f7.2,
     1         i2,'/',a4,'-',i2,'/',a4,'-',i2,'/',a4,'-',i2,'/',a4)
          endif
c----------------------------------------------------------------double torsions
        if(j.eq.14) then
          i1=inoe(k,1)
          i2=inoe(k,2)
          i3=inoe(k,3)
          i4=inoe(k,4)
          dih1=torp(i1,i2,i3,i4)
          i5=inoe(k,5)
          i6=inoe(k,6)
          i7=inoe(k,7)
          i8=inoe(k,8)
          dih2=torp(i5,i6,i7,i8)
          if(inoe(k,9).gt.0) then
            dih=dih1+dih2
            add='+'
          else
            dih=dih1-dih2
            add='-'
          endif
          delt=dih-rnoe(k)
          if(abs(delt).gt.180.) delt=delt-sign(360.d0,delt)
          write(6,71) dih,rnoe(k),delt,nunit(i1),mnam(i1),
     1      nunit(i2),mnam(i2),nunit(i3),mnam(i3),nunit(i4),mnam(i4),add,
     1      nunit(i5),mnam(i5),nunit(i6),mnam(i6),nunit(i7),mnam(i7),
     1      nunit(i8),mnam(i8)
71        format(2x,'T2 Val:',f7.2,' Req: ',6x,f7.2,6x,' Del: ',f7.2,
     1           i2,'/',a4,'-',i2,'/',a4,'-',i2,'/',a4,'-',i2,'/',a4,
     1      /52x,a1,1x,i2,'/',a4,'-',i2,'/',a4,'-',i2,'/',a4,'-',i2,'/',a4)
        else if(j.eq.16) then
          i1=inoe(k,1)
          i2=inoe(k,2)
          i3=inoe(k,3)
          i4=inoe(k,4)
          dih1=torp(i1,i2,i3,i4)
          i5=inoe(k,5)
          i6=inoe(k,6)
          i7=inoe(k,7)
          i8=inoe(k,8)
          dih2=torp(i5,i6,i7,i8)
          if(inoe(k,9).gt.0) then
            dih=dih1+dih2
            add='+'
          else
            dih=dih1-dih2
            add='-'
          endif
          delt=dih-rnoe(k)
          cx=cos(cdr*(bnoe(k)))
          sx=sin(cdr*(bnoe(k)))
          cm=cos(cdr*(rnoe(k)))
          sm=sin(cdr*(rnoe(k)))
          ct=cos(cdr*(dih))
          st=sin(cdr*(dih))
          v=cx*sm-sx*cm
          a=cm*st-sm*ct
          b=ct*sx-st*cx
          if((a.ge.0.and.b.ge.0).or.(v.gt.0.and.a*b.le.0)) then
            delt=0.
            goto 211
          endif
          if(a.gt.0.and.b.lt.0) then
            delt=dih-bnoe(k)
          else if(a.lt.0.and.b.lt.0) then
            a=ct*cm+st*sm
            b=cx*ct+sx*st
            if(a.lt.b) then
              delt=dih-bnoe(k)
            endif
          endif
211       if(abs(delt).gt.180.) delt=delt-sign(360.d0,delt)
          write(6,72) dih,rnoe(k),bnoe(k),delt,nunit(i1),mnam(i1),
     1      nunit(i2),mnam(i2),nunit(i3),mnam(i3),nunit(i4),mnam(i4),add,
     1      nunit(i5),mnam(i5),nunit(i6),mnam(i6),nunit(i7),mnam(i7),
     1      nunit(i8),mnam(i8)
72          format(2x,'t) Val: ',f7.2,' Req: ',f7.2,' <-> ',f7.2,
     1         ' Del: ',f7.2,
     1           i2,'/',a4,'-',i2,'/',a4,'-',i2,'/',a4,'-',i2,'/',a4,
     1      /52x,a1,1x,i2,'/',a4,'-',i2,'/',a4,'-',i2,'/',a4,'-',i2,'/',a4)
        endif
c---------------------------------------------------------------double distances
        m=0
        if(jnoe(k).eq.26.or.jnoe(k).eq.28) then
          i1=inoe(k,1)
          i2=inoe(k,2)
          xx1=corm(i2,1)-corm(i1,1)
          yy1=corm(i2,2)-corm(i1,2)
          zz1=corm(i2,3)-corm(i1,3)
          r1=sqrt(xx1*xx1+yy1*yy1+zz1*zz1)
          i3=inoe(k,3)
          i4=inoe(k,4)
          xx2=corm(i4,1)-corm(i3,1)
          yy2=corm(i4,2)-corm(i3,2)
          zz2=corm(i4,3)-corm(i3,3)
          r2=sqrt(xx2*xx2+yy2*yy2+zz2*zz2)
          r=r1+inoe(k,5)*r2
        endif
        if(jnoe(k).eq.26) then
          delt=r-rnoe(k)
          m=m+1
          if (inoe(k,5).gt.0) then
            add='+'
          else
            add='-'
          endif
          write(6,98) r,rnoe(k),delt,nunit(i1),mnam(i1),
     1       nunit(i2),mnam(i2),add,nunit(i3),mnam(i3),nunit(i4),
     1       mnam(i4)
98        format(2x,'D2 Val:',f7.2,' Req: ',6x,f7.2,6x,' Del: ',f7.2,'-',
     1        i2,'/',a4,'-',i2,'/',a4,/53x,a1,1x,i2,'/',a4,'-',i2,'/',a4)
        else if(j.eq.28) then
          if (inoe(k,5).gt.0) then
             add='+'
          else
             add='-'
          endif
          m=m+1
          if (r.ge.rnoe(k).and.r.le.bnoe(k)) then
            delt=0
          elseif (r.lt.rnoe(k)) then
            delt=r-rnoe(k)
          elseif (r.gt.bnoe(k)) then
            delt=r-bnoe(k)
          endif
          write(6,99) r,rnoe(k),bnoe(k),delt,nunit(i1),mnam(i1),
     1          nunit(i2),mnam(i2),add,
     1          nunit(i3),mnam(i3),nunit(i4),mnam(i4)
99        format(2x,'d2 Val: ',f7.2,' Req: ',f7.2,' <-> ',f7.2,
     1         ' Del: ',f7.2,'-',
     1      i2,'/',a4,'-',i2,'/',a4,/54x,a1,1x,i2,'/',a4,'-',i2,'/',a4)
        endif
c--------------------------------------------------------------------sugar phase
        if(j.eq.5.or.j.eq.21) then
           if(j.eq.21) then
             ip=inoe(k,2)
             rnoe(k)=phase(ip,3)
           endif
           is=inoe(k,1)
           delt=ph5(is)-rnoe(k)
           if(abs(delt).gt.180.) delt=delt-sign(360.d0,delt)
           write(6,83) ph5(is),rnoe(k),delt,is
83         format(2x,'P) Val: ',f7.2,' Req: ',6x,f7.2,6x,' Del: ',f7.2,
     1       ' Sug:',i3)
        endif
c-------------------------------------------------------------------groove width
        if(j.eq.31.or.j.eq.32) then
          do ks=1,inoe(k,9)
            rmin=10000.
            m=nuc(inoe(k,ks)-1)+12
            x0=corm(m,1)
            y0=corm(m,2)
            z0=corm(m,3)
            do ls=iend(1)+1,iend(2)
              l=nuc(ls-1)+12
              dx=corm(l,1)-x0
              dy=corm(l,2)-y0
              dz=corm(l,3)-z0
              r=sqrt(dx*dx+dy*dy+dz*dz)
              if(r.lt.rmin) then
                rmin=r
                lm=l
                xx=dx
                yy=dy
                zz=dz
              endif
            enddo
            r=rmin
            if(j.eq.31) then
              delt=rmin-rnoe(k)
              write(6,77) rmin,rnoe(k),delt,mnam(m),munit(m),
     1           nunit(m),mnam(lm),munit(lm),nunit(lm)
77            format(2x,'H) Val: ',f7.2,' Req: ',f7.2,
     1         ' Del: ',f7.2,2x,a4,1x,a4,i3,' : ',a4,1x,a4,i3)
            else
              if(rmin.lt.rnoe(k)) then
                 delt=rmin-rnoe(k)
              else if(rmin.gt.bnoe(k)) then
                 delt=rmin-bnoe(k)
              else
                 delt=0.
              endif
              write(6,78) rmin,rnoe(k),bnoe(k),delt,mnam(m),munit(m),
     1                    nunit(m),mnam(lm),munit(lm),nunit(lm)
78            format(2x,'h) Val: ',f7.2,' Req: ',f7.2,' <-> ',f7.2,
     1               ' Del: ',f7.2,2x,a4,1x,a4,i2,' : ',a4,1x,a4,i2)
            endif
          enddo
        endif
c--------------------------------------------------------------sugar phase range
        if(j.eq.6.or.j.eq.22) then
          if(j.eq.22) then
            ip=inoe(k,2)
            rnoe(k)=phase(ip,1)
            bnoe(k)=phase(ip,2)
          endif
          is=inoe(k,1)
          pha=ph5(is)
          delt=pha-rnoe(k)
          cx=cos(cdr*(bnoe(k)))
          sx=sin(cdr*(bnoe(k)))
          cm=cos(cdr*(rnoe(k)))
          sm=sin(cdr*(rnoe(k)))
          ct=cos(cdr*(pha))
          st=sin(cdr*(pha))
          v=cx*sm-sx*cm
          a=cm*st-sm*ct
          b=ct*sx-st*cx
          if((a.ge.0.and.b.ge.0).or.(v.gt.0.and.a*b.le.0)) then
             delt=0.
             goto 202
          endif
          if(a.gt.0.and.b.lt.0) then
            delt=pha-bnoe(k)
          else if(a.lt.0.and.b.lt.0) then
            a=ct*cm+st*sm
            b=cx*ct+sx*st
            if(a.lt.b) then
              delt=pha-bnoe(k)
            endif
          endif
202       if(abs(delt).gt.180.) delt=delt-sign(360.d0,delt)
            write(6,84) pha,rnoe(k),bnoe(k),delt,is
84          format(2x,'p) Val: ',f7.2,' Req: ',f7.2,' <-> ',f7.2,
     1       ' Del: ',f7.2,' Sug:',i3)
        endif
c-------------------------------------------------------------------double phase
        if(j.eq.18) then
          is1=inoe(k,1)
          is2=inoe(k,2)
          if(inoe(k,9).gt.0) then
            add='+'
            delt=ph5(is1)+ph5(is2)-rnoe(k)
          else
            add='-'
            delt=ph5(is1)-ph5(is2)-rnoe(k)
          endif
          if(abs(delt).gt.180.) delt=delt-sign(360.d0,delt)
          write(6,73) delt+rnoe(k),rnoe(k),delt,is1,add,is2
73        format(2x,'P2 Val: ',f7.2,' Req: ',6x,f7.2,6x,' Del: ',f7.2,
     1      ' Sug:',i3,1x,a1,i3)
        else if(j.eq.19) then
          is1=inoe(k,1)
          is2=inoe(k,2)
          if(inoe(k,9).gt.0) then
            add='+'
            dpha=ph5(is1)+ph5(is2)
          else
            add='-'
            dpha=ph5(is1)-ph5(is2)
          endif
          if(dpha.lt.rnoe(k)) then
            delt=dpha-rnoe(k)
          else if(dpha.gt.bnoe(k)) then
            delt=dpha-bnoe(k)
          else
            delt=0.
          endif
          if(abs(delt).gt.180.) delt=delt-sign(360.d0,delt)
          write(6,74) dpha,rnoe(k),bnoe(k),delt,is1,add,is2
74        format(2x,'p2 Val: ',f7.2,' Req: ',f7.2,' <-> ',f7.2,
     1      ' Del: ',f7.2,' Sug:',i3,1x,a1,i3)
        endif
c----------------------------------------------------------------sugar amplitude
        if(j.eq.7.or.j.eq.21) then
           if(j.eq.21) then
             ip=inoe(k,2)
             rnoe(k)=ampli(ip,3)
           endif
           is=inoe(k,1)
           write(6,85) am5(is),rnoe(k),am5(is)-rnoe(k),is
85         format(2x,'A) Val: ',f7.2,' Req: ',6x,f7.2,6x,' Del: ',f7.2,
     1         ' Sug:',i3)
        endif
c----------------------------------------------------------sugar amplitude range
      if(j.eq.8.or.j.eq.22) then
         if(j.eq.22) then
         ip=inoe(k,2)
         rnoe(k)=ampli(ip,1)
         bnoe(k)=ampli(ip,2)
         endif
      is=inoe(k,1)
      amp=am5(is)
      if(amp.lt.rnoe(k)) then
      delt=amp-rnoe(k)
      else if(amp.gt.bnoe(k)) then
      delt=amp-bnoe(k)
      else
      delt=0.
      endif
      write(6,86) amp,rnoe(k),bnoe(k),delt,is
86    format(2x,'a) Val: ',f7.2,' Req: ',f7.2,' <-> ',f7.2,
     1 ' Del: ',f7.2,' Sug:',i3)
      endif
c---------------------------------------------------------------double amplitude
      if(j.eq.24) then
      is1=inoe(k,1)
      is2=inoe(k,2)
      if(inoe(k,9).gt.0) then
      add='+'
      delt=am5(is1)+am5(is2)-rnoe(k)
      else
      add='-'
      delt=am5(is1)-am5(is2)-rnoe(k)
      endif
      if(abs(delt).gt.180.) delt=delt-sign(360.d0,delt)
      write(6,75) delt+rnoe(k),rnoe(k),delt,is1,add,is2
75    format(2x,'A2 Val: ',f7.2,' Req: ',6x,f7.2,6x,' Del: ',f7.2,
     1 ' Sug:',i3,1x,a1,i3)
      else if(j.eq.25) then
      is1=inoe(k,1)
      is2=inoe(k,2)
      if(inoe(k,9).gt.0) then
      add='+'
      dpha=am5(is1)+am5(is2)
      else
      add='-'
      dpha=am5(is1)-am5(is2)
      endif
      if(dpha.lt.rnoe(k)) then
      delt=dpha-rnoe(k)
      else if(dpha.gt.bnoe(k)) then
      delt=dpha-bnoe(k)
      else
      delt=0.
      endif
      if(abs(delt).gt.180.) delt=delt-sign(360.d0,delt)
      write(6,76) dpha,rnoe(k),bnoe(k),delt,is1,add,is2
76    format(2x,'a2 Val: ',f7.2,' Req: ',f7.2,' <-> ',f7.2,
     1 ' Del: ',f7.2,' Sug:',i3,1x,a1,i3)
      endif
c------------------------------------------------------------------valence angle
      if(j.eq.9) then
      i1=inoe(k,1)
      i2=inoe(k,2)
      i3=inoe(k,3)
      val=ang(i1,i2,i3)
      delt=val-rnoe(k)
      if(abs(delt).gt.180.) delt=delt-sign(360.d0,delt)
      write(6,89) val,rnoe(k),delt,nunit(i1),
     1 mnam(i1),nunit(i2),mnam(i2),nunit(i3),mnam(i3)
89    format(2x,'V) Val: ',f7.2,' Req: ',6x,f7.2,6x,' Del: ',f7.2,
     1 i2,'/',a4,'-',i2,'/',a4,'-',i2,'/',a4)
      endif
c------------------------------------------------------------------------opening
      if(j.eq.29) then
      i1=inoe(k,1)
      i2=inoe(k,2)
      i3=inoe(k,3)
      i4=inoe(k,4)
      i5=inoe(k,5)
      i6=inoe(k,6)
      i7=inoe(k,7)
      gx=corm(i2,1)-corm(i1,1)
      gy=corm(i2,2)-corm(i1,2)
      gz=corm(i2,3)-corm(i1,3)
      cx=corm(i3,1)-corm(i1,1)
      cy=corm(i3,2)-corm(i1,2)
      cz=corm(i3,3)-corm(i1,3)
      rc=sqrt(cx*cx+cy*cy+cz*cz)
      dotgc=gx*cx+gy*cy+gz*cz
      s1x=corm(i5,1)-corm(i4,1)
      s1y=corm(i5,2)-corm(i4,2)
      s1z=corm(i5,3)-corm(i4,3)
      rs1=sqrt(s1x*s1x+s1y*s1y+s1z*s1z)
      s2x=corm(i7,1)-corm(i6,1)
      s2y=corm(i7,2)-corm(i6,2)
      s2z=corm(i7,3)-corm(i6,3)
      rs2=sqrt(s2x*s2x+s2y*s2y+s2z*s2z)
      tx=s1x/rs1+s2x/rs2
      ty=s1y/rs1+s2y/rs2
      tz=s1z/rs1+s2z/rs2
      dtc=(tx*cx+ty*cy+tz*cz)/rc
      ux=tx-dtc*cx/rc
      uy=ty-dtc*cy/rc
      uz=tz-dtc*cz/rc
      ru2=ux*ux+uy*uy+uz*uz
      ru=sqrt(ru2)
      dug=(gx*ux+gy*uy+gz*uz)/ru2
      px=gx-dug*ux
      py=gy-dug*uy
      pz=gz-dug*uz
      rp=sqrt(px*px+py*py+pz*pz)
      dpu= px*ux+py*uy+pz*uz
      dpc=(px*cx+py*cy+pz*cz)/(rp*rc)
      val=acos(dpc)*crd
         dx=cy*pz-cz*py
         dy=cz*px-cx*pz
         dz=cx*py-cy*px
         dot=(dx*ux+dy*uy+dz*uz)*inoe(k,8)
         if(dot.lt.0) val=-val
      delt=val-rnoe(k)
      if(abs(delt).gt.180.) delt=delt-sign(360.d0,delt)
      write(6,97) val,rnoe(k),delt,nunit(i1),
     1 mnam(i1),nunit(i2),mnam(i2),nunit(i3),mnam(i3)
97    format(2x,'O) Val: ',f7.2,' Req: ',6x,f7.2,6x,' Del: ',f7.2,
     1 i2,'/',a4,'-',i2,'/',a4,'-',i2,'/',a4)
      endif
c------------------------------------------------------------------------X twist
      if(j.eq.30) then
      ic=0
      vtot=0.
      do i=1,inoe(k,1)
      i1=icon(ic+1)
      i2=icon(ic+2)
      i3=icon(ic+3)
      i4=icon(ic+4)
      do j=1,3
      v1(j)=corm(i3,j)-corm(i1,j)
      v2(j)=corm(i4,j)-corm(i2,j)
      s1(j)=corm(i2,j)-corm(i1,j)
      s2(j)=corm(i4,j)-corm(i3,j)
      enddo
      rv1=sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
      rv2=sqrt(v2(1)**2+v2(2)**2+v2(3)**2)
      do j=1,3
      t(j)=v1(j)/rv1+v2(j)/rv2
      enddo
      rt=sqrt(t(1)*t(1)+t(2)*t(2)+t(3)*t(3))
      t(1)=t(1)/rt
      t(2)=t(2)/rt
      t(3)=t(3)/rt
      s1t=s1(1)*t(1)+s1(2)*t(2)+s1(3)*t(3)
      s2t=s2(1)*t(1)+s2(2)*t(2)+s2(3)*t(3)
      do j=1,3
      q1(j)=s1(j)-t(j)*s1t
      q2(j)=s2(j)-t(j)*s2t
      enddo
      rq1=sqrt(q1(1)**2+q1(2)**2+q1(3)**2)
      rq2=sqrt(q2(1)**2+q2(2)**2+q2(3)**2)
      dotq=q1(1)*q2(1)+q1(2)*q2(2)+q1(3)*q2(3)
      val=acos(dotq/(rq1*rq2))*crd
         px=q1(2)*q2(3)-q1(3)*q2(2)
         py=q1(3)*q2(1)-q1(1)*q2(3)
         pz=q1(1)*q2(2)-q1(2)*q2(1)
         dpt=px*t(1)+py*t(2)+pz*t(3)
         if(dpt.lt.0) val=-val
      vtot=vtot+val
      ic=ic+2
      enddo
      vlr=vtot*cdr
      rnr=rnoe(k)*cdr
      delt=acos(cos(rnr)*cos(vlr)+sin(rnr)*sin(vlr))*crd
      cros=cos(rnr)*sin(vlr)-sin(rnr)*cos(vlr)
      if(cros.lt.0) delt=-delt
      add='X'
      if(inoe(k,1).gt.1) add='x'
      write(6,50) add,rnoe(k)+delt,rnoe(k),delt,nunit(i1),nunit(i3)
50    format(2x,a1,') Val: ',f7.2,' Req: ',6x,f7.2,6x,' Del: ',f7.2,
     1 4x,i2,' - ',i2)
      endif
c-----------------------------------------------------------cosine torsion angle
        if(j.eq.10) then
           i1=inoe(k,1)
      i2=inoe(k,2)
      i3=inoe(k,3)
      i4=inoe(k,4)
      dih=torp(i1,i2,i3,i4)
      write(6,90) dih,rnoe(k),nunit(i1),mnam(i1),
     1 nunit(i2),mnam(i2),nunit(i3),mnam(i3),nunit(i4),mnam(i4)
90    format(2x,'C) Val: ',f7.2,' Mul: ',6x,f7.2,19x,
     1 i2,'/',a4,'-',i2,'/',a4,'-',i2,'/',a4,'-',i2,'/',a4)
      endif
c------------------------------------------------------------------chemical bond
      if(j.eq.11) then
      i1=inoe(k,1)
      i2=inoe(k,2)
      dx=corm(i1,1)-corm(i2,1)
      dy=corm(i1,2)-corm(i2,2)
      dz=corm(i1,3)-corm(i2,3)
      r=sqrt(dx*dx+dy*dy+dz*dz)
      delt=r-rnoe(k)
      write(6,91) r,rnoe(k),delt,nunit(i1),
     1 mnam(i1),nunit(i2),mnam(i2)
91    format(2x,'B) Val: ',f7.2,' Req: ',6x,f7.2,6x,' Del: ',f7.2,
     1 i2,'/',a4,'-',i2,'/',a4)
      endif
c---------------------------------------------------------------------total rise
      if(j.eq.12) then
      iu1=inoe(k,1)
      iu2=inoe(k,2)
      zsum=0.
      do ih=iu1,iu2
      zsum=zsum+hel(ih,3)
      enddo
      zave=zsum/(iu2-iu1+1)
      delt=zave-rnoe(k)
      write(6,92) zave,rnoe(k),delt,iu1,iu2
92    format(2x,'R) Val: ',f7.2,' Req: ',6x,f7.2,6x,' Del: ',f7.2,
     1 ' Units: ',i3,' - ',i3)
      endif
c--------------------------------------------------------------------total twist
      if(j.eq.13) then
      iu1=inoe(k,1)
      iu2=inoe(k,2)
      wsum=0.
      do ih=iu1,iu2
      wsum=wsum+hel(ih,6)
      enddo
      wave=wsum/(iu2-iu1+1)
      delt=wave-rnoe(k)
      write(6,93) wave,rnoe(k),delt,iu1,iu2
93    format(2x,'W) Val: ',f7.2,' Req: ',6x,f7.2,6x,' Del: ',f7.2,
     1 ' Units: ',i3,' - ',i3)
      endif
c------------------------------------------------------------------variable sums
      if(j.eq.15.or.j.eq.17) then
      sumv=0.
      do iv=1,j-13
      ia=inoe(k,iv)
      is=sign(1,ia)
      ia=abs(ia)
         do jv=1,nvar
         if(neq(jv).eq.ia) goto 222
         enddo
      write(6,95) rnoe(k)
95    format(2x,'S) Val:   ---- ',' Req: ',6x,f7.2,6x,' Del:   ---- ')
      goto 223
222   sumv=sumv+is*var(jv)
      enddo
      delt=sumv-rnoe(k)
      write(6,94) sumv,rnoe(k),delt
94    format(2x,'S) Val: ',f7.2,' Req: ',6x,f7.2,6x,' Del: ',f7.2)
      endif
c----------------------------------------------------------------------curvature
      if(j.eq.20) then
        i1=inoe(k,1)
        i2=inoe(k,2)
      the=rnoe(k)
      ct=cos(the*cdr)
      st=sin(the*cdr)
      phi=rnoe(k+1)
      cp=cos(phi*cdr)
      sp=sin(phi*cdr)
      x0=ra(i1,1)
      y0=ra(i1,2)
      z0=ra(i1,3)
      sx=ra(i2,1)-x0
      sy=ra(i2,2)-y0
      sz=ra(i2,3)-z0
      rs=sqrt(sx*sx+sy*sy+sz*sz)
      ux=ua(i1,1)
      uy=ua(i1,2)
      uz=ua(i1,3)
      ax=da(i1,1)
      ay=da(i1,2)
      az=da(i1,3)
      ex=uy*az-uz*ay
      ey=uz*ax-ux*az
      ez=ux*ay-uy*ax
      wx=ax*cp+ex*sp
      wy=ay*cp+ey*sp
      wz=az*cp+ez*sp
      v0x=ux*ct+wx*st
      v0y=uy*ct+wy*st
      v0z=uz*ct+wz*st
      dot=ua(i2,1)*v0x+ua(i2,2)*v0y+ua(i2,3)*v0z
      if(abs(dot).gt.1.d0) dot=sign(1.d0,dot)
      dcr=acos(dot)*crd
      dot=ua(i2,1)*ux+ua(i2,2)*uy+ua(i2,3)*uz
      if(abs(dot).gt.1.d0) dot=sign(1.d0,dot)
      ben=sign(acos(dot)*crd,the)
      delt=ben-the
      if(ben.ne.0) then
      dot=ux*sx+uy*sy+uz*sz
      qx=sx-ux*dot
      qy=sy-uy*dot
      qz=sz-uz*dot
      rq=sqrt(qx*qx+qy*qy+qz*qz)
      dir=acos((ax*qx+ay*qy+az*qz)/rq)*crd
      px=ay*qz-az*qy
      py=az*qx-ax*qz
      pz=ax*qy-ay*qx
      if(ux*px+uy*py+uz*pz.lt.0) dir=-dir
      if(the.lt.0) dir=dir-sign(180.d0,dir)
         else
         dir=phi
         endif
      rcd=rs/(2*sin(abs(ben)*cdr/2))
      if(rcd.gt.9999) rcd=9999.
      deld=dir-phi
      write(6,88) ben,the,delt,i1,i2,dir,phi,deld,rcd,rs,dcr
88    format(2x,'U) Ang: ',f7.2,' Req: ',f7.2,' Del: ',f7.2,
     1 ' between units ',i2,' and ',i2,
     1      /2x,'   Dir: ',f7.2,' Req: ',f7.2,' Del: ',f7.2,
     1      /2x,'   Rad: ',f7.1,' Len: ',f7.2,' Dcr: ',f7.2)
      endif
223   enddo
c--------------------------------------------------modified nucleotide variables
      if(.not.sum) then
      start=.true.
      k=0
      do is=1,nto
      in=itr(is)
      ino=ito(is)
      io=nuc(is-1)+iofs(is)
      kl=k+kap(1,ino)+1
      ku=k+kap(1,in)
      if(ku.ge.kl) then
      if(start) then
      start=.false.
      write(6,22)
22    format(/2x,'Covlig',12x,'x1',6x,'x2',6x,'x3',6x,'x4',6x,
     1 'x5',6x,'x6',6x,'x7'/)
      endif
      write(6,24) is,munit(io),nunit(io),(sap(k),k=kl,ku)
24    format(2x,i3,':',a4,i3,2x,7f8.2,10(:/15x,7f8.2))
      endif
      k=kapt(is+1)
      enddo
c------------------------------------------summary modified nucleotide variables
      else
      k=0
      do is=1,nto
      in=itr(is)
      ino=ito(is)
      io=nuc(is-1)+iofs(is)
      ioff=nuc(is-1)
      kh=k
      lv=0
      do l=1,kap(2,in)
      k=k+1
      idel=kap(1,in)-kap(1,ino)
      iano=kap(2,in)-kap(2,ino)-idel-5*(nsr(is)-1)
      if((l.gt.kap(1,ino).and.l.le.kap(1,in)).or.
     1   (l.gt.kap(1,in)+5.and.l.le.kap(1,in)+5*nsr(is)).or.
     1    l.gt.kap(2,in)-iano) then
      lv=lv+1
      ia=nap(l,1,in)+ioff
      ib=nap(l,2,in)+ioff
      ic=nap(l,3,in)+ioff
      id=nap(l,4,in)+ioff
      ixv(lv)=k
      vix(lv,1)=mnam(ia)
      vix(lv,2)=mnam(ib)
      vix(lv,3)=mnam(ic)
      id=nap(l,4,in)+ioff
      vix(lv,4)='     '
      if(id.ne.ioff) vix(lv,4)='-'//mnam(id)
      vin(lv)='I'
      if(l.gt.kap(1,in)) vin(lv)='D'
      if(nap(l,6,in).lt.0) vin(lv)='A'
      if(l.gt.kap(3,in)) vin(lv)='A'
      vlo(lv)=' '
      if(lock(k)) vlo(lv)='L'
      endif
      enddo
c-------------------------------------------------------------------------output
      if(lv.gt.0) then
      write(6,15) is,munit(io),nunit(io)
15    format(/2x,'Covlig intra summary of ',i2,': ',a4,i3/)
      m=(lv+1)/2
      do i=1,m
      ip=i+m
      i1=ixv(i)
      i2=ixv(ip)
      if(ip.le.lv) then
      write(6,25) i1-kh,(vix(i,j),j=1,4),sap(i1),vin(i),vlo(i),
     1 i2-kh,(vix(ip,j),j=1,4),sap(i2),vin(ip),vlo(ip)
      else
      write(6,25) i1-kh,(vix(i,j),j=1,4),sap(i1),vin(i),vlo(i)
      endif
      enddo
      endif
      enddo
      endif
c---------------------------------------------------------intra ligand variables
      if(nlig.ne.0.and.kapt(ntl+1).gt.kapt(nto+1)) then
      if(.not.sum) then
      write(6,42)
42    format(/2x,'Ligand intra',5x,'L1',6x,'L2',6x,'L3',6x,'L4',6x,
     1 'L5',6x,'L6',6x,'L7'/)
      ki=nsph
      do il=1,nlig
      is=nto+il
      io=nuc(is-1)+1
      ioff=nuc(is-1)
      in=itr(is)
      kl=kapt(is)+1
      ku=kl+kap(1,in)-1
      do l=1,kap(2,in)
      k=kapt(is)+l
      ia=nap(l,1,in)+ioff
      ib=nap(l,2,in)+ioff
      ic=nap(l,3,in)+ioff
      id=nap(l,4,in)+ioff
      if(id.eq.ioff) then
      sap(k)=ang(ia,ib,ic)
      else
      sap(k)=torp(ia,ib,ic,id)
      endif
c----------------------------------------------------------------------reset var
      if(l.le.kap(1,in).and..not.lock(k)) then
      ki=ki+1
      var(ki)=sap(k)
      endif
      enddo
      if(ku.ne.kl) write(6,44) is,munit(io),nunit(io),(sap(k),k=kl,ku)
44    format(2x,i3,':',a4,i3,2x,7f8.2,10(:/15x,7f8.2))
      enddo
      else
c-------------------------------------------------------summary ligand variables
      ki=nsph
      k=kapt(nto+1)
      do il=1,nlig
      is=nto+il
      in=itr(is)
      ir=irec(is)
      ioff=nuc(is-1)
      lv=0
      kh=k
      do l=1,kap(2,in)
      k=k+1
      lv=lv+1
      ia=nap(l,1,in)+ioff
      ib=nap(l,2,in)+ioff
      ic=nap(l,3,in)+ioff
      id=nap(l,4,in)+ioff
      vix(lv,1)=mnam(ia)
      vix(lv,2)=mnam(ib)
      vix(lv,3)=mnam(ic)
      id=nap(l,4,in)+ioff
      if(id.eq.ioff) then
      sap(k)=ang(ia,ib,ic)
      vix(lv,4)='     '
      else
      sap(k)=torp(ia,ib,ic,id)
      vix(lv,4)='-'//mnam(id)
      endif
      vin(lv)='I'
      if(l.gt.kap(1,in)) vin(lv)='D'
      if(nap(l,6,in).lt.0) vin(lv)='A'
      vlo(lv)=' '
      if(lock(k)) vlo(lv)='L'
c----------------------------------------------------------------------reset var
      if(l.le.kap(1,in).and..not.lock(k)) then
      ki=ki+1
      var(ki)=sap(k)
      endif
      enddo
c-------------------------------------------------------------------------output
      if(lv.gt.0) then
      write(6,35) il,lnam(il),kas(ir)
35    format(/2x,'Ligand intra summary of ',i2,': ',a4,' (',
     1 i3,' atoms)'/)
      m=(lv+1)/2
      do i=1,m
      ip=i+m
      if(ip.le.lv) then
      write(6,25) i,(vix(i,j),j=1,4),sap(kh+i),vin(i),vlo(i),
     1 ip,(vix(ip,j),j=1,4),sap(kh+ip),vin(ip),vlo(ip)
      else
      write(6,25) i,(vix(i,j),j=1,4),sap(kh+i),vin(i),vlo(i)
      endif
      enddo
25    format(2x,'(',i2,') ',a4,'-',a4,'-',a4,a5,2x,f7.2,1x,2a1,
     1     :,5x,'(',i2,') ',a4,'-',a4,'-',a4,a5,2x,f7.2,1x,2a1)
      endif
      enddo
      endif
      endif
      return
      end
