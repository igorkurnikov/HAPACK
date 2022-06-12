      subroutine equim(newq,nspc,nhlc,nbcc,nklc,nlgc,istop)
c
c  check and list symmetry contractions
c
      include 'jumna_data.inc'
      logical*2 lar,lock,kink,lthy,locr,cation,fhel,
     1 sup,rcom,homo,homo2,homo3,diep,link,ecen,cyl,ribose,lcat,
     1 cent,angle,start,brk(n2),autos,ihl,amber,sum
      character*4 mnam,munit,seq*120,code*8,kode*8,lnam,rb*1,ct*1,ch(18)
      integer*4 opt,ipl(n8),ipr(n7),ient(n2),newq(n2,18),nu(n7)
      dimension chec(n7),ches(n7),win(n2)
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
      common/mnn/var(n7),gra(n7),con(n0*n2+n3),scl(n7),ncon,nvrc,ntba,
     1 nbac,nthe,nhel,ntki,nkin,ntri,ntot,nvar,nrin,lar(n7),lock(n6a)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/sst/phase(8,3),ampli(8,3),ph5(n2),am5(n2),kr5(n2)
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(120,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      common/symjm/xsum,ysum,itwl(n2),nvs(n7),ihm(n2),ksym(3),nbrk(2),
     1 isym,isur,isup,ihl(n2)
	common/hacon/iret
      data ch/'C1''','C2''','C3''','Gly','C12','C23','Eps','Zet',
     1'C3''','P  ','Rib','Thy','Xdis','Ydis','Rise','Inc','Tip','Twis'/
c------------------------------------------------------------------------
      k=0
      do i=1,nto
        in=itr(i)
        ino=ito(i)
        ienb(i,1)=k+1
        ienb(i,2)=k+kap(1,ino)
        ient(i)=k+kap(1,in)
        k=k+kap(2,in)
      enddo
      do j=1,kseq
        brk(j)=.false.
        if(j.eq.nbrk(1).or.j.eq.nbrk(2)) brk(j)=.true.
      enddo
c---------------------------------------------------------------------neq matrix
      do i=1,nspc
        ipr(i)=0
      enddo
      k=0
      ki=0
      nvax=0
      nmax=0
      do is=1,nto
        in=itr(is)
        ino=ito(is)
        do iv=1,11
          n=newq(is,iv)
          na=abs(n)
          if(n.ne.0) then
            ki=ki+1
            if(na.gt.nmax) nmax=na
              neq(ki)=n+sign(ipr(na),n)
              if(abs(neq(ki)).gt.nvax) nvax=abs(neq(ki))
            endif
        enddo
        kto=nvax
        k=k+kap(1,ino)
        do iv=1,kap(1,in)-kap(1,ino)
          k=k+1
          if(.not.lock(k)) then
            ki=ki+1
            kto=kto+1
            neq(ki)=kto
          endif
        enddo
      kdel=kto-nvax
      nvax=nvax+kdel
      do i=nmax+1,nspc
      ipr(i)=ipr(i)+kdel
      enddo
      k=k+kap(2,in)-kap(1,in)
      enddo
 
      k=kapt(nto+1)
      ki=nsph
      kto=nvax
      do il=1,nlig
      is=nto+il
      in=itr(is)
      do i=1,kap(1,in)
      k=k+1
      if(.not.lock(k)) then
      ki=ki+1
      kto=kto+1
      neq(ki)=kto
      endif
      enddo
      k=k+kap(2,in)-kap(1,in)
      enddo
      ndel=kto-nspc
 
      k=nbac
      do in=1,nto
      do iv=13,18
      n=newq(in,iv)
      if(n.ne.0) then
      k=k+1
      nd=n+sign(ndel,n)
      neq(k)=nd
      if(abs(nd).gt.kto) kto=abs(nd)
      endif
      enddo
      enddo
 
      k=nthe
      ki=nhel
      nhlc=kto
      do is=2,kseq
      if(kink(is)) then
      do l=1,4
      k=k+1
      if(.not.lock(k)) then
      ki=ki+1
      kto=kto+1
      neq(ki)=kto
      endif
      enddo
      endif
      enddo
 
      k=ntki
      ki=nkin
      do il=1,nlig
      is=nto+il
      kal=nuc(is)-nuc(is-1)
      jlim=6
      if(kal.eq.1) jlim=3
      do j=1,jlim
      k=k+1
      if(.not.lock(k)) then
      ki=ki+1
      kto=kto+1
      neq(ki)=kto
      endif
      enddo
      enddo
      ndel=ndel+kto-nhlc
 
      k=nlgi
      do in=1,nto
      n=newq(in,12)
      if(n.ne.0) then
      k=k+1
      neq(k)=n+sign(ndel,n)
      endif
      enddo
 
c---------------------------------------------------------------------superhelix
      if(nshel.eq.1.or.nshel.eq.3) then
      k=k+1
      neq(isur)=k
      endif
         if(nshel.ge.2) then
         k=k+1
         neq(isup)=k
         endif
c=================================================================compact sugars
      if(rcom) then
      do j=1,nvar
      nu(j)=1
      enddo
      do is=1,nto
      k=kr5(is)
      if(k.ne.0) then
      do j=k-4,k-2
      nu(abs(neq(j)))=0
      neq(j)=0
      enddo
      endif
      enddo
         k=0
         do j=1,nvar
         if(nu(j).ne.0) then
         k=k+1
         nu(j)=k
         endif
         enddo
      do j=1,nvar
      nj=abs(neq(j))
      if(nj.ne.0) neq(j)=sign(nu(nj),neq(j))
      enddo
      endif
c================================================================ligand symmetry
      do i=1,nvar
      j=abs(neq(i))
      if(j.gt.0) then
      if(nvs(j).ne.0) then
      neq(i)=nvs(j)
      endif
      endif
      enddo
      do i=1,nvar
      nvs(i)=0
      enddo
      do i=1,nvar
      ia=abs(neq(i))
      if(ia.gt.0) nvs(ia)=1
      enddo
      ic=1
      do i=1,nvar
      if(nvs(i).ne.0) then
      nvs(i)=ic
      ic=ic+1
      endif
      enddo
      do i=1,nvar
      ia=abs(neq(i))
      if(ia.gt.0) neq(i)=sign(nvs(ia),neq(i))
      enddo
c======================================================setup contracted counters
      k=0
      do i=1,nsph
      in=abs(neq(i))
      if(in.gt.k) k=in
      enddo
      nspc=k
      do i=nsph+1,nbac
      in=abs(neq(i))
      if(in.gt.k) k=in
      enddo
      nbcc=k
      do i=nbac+1,nhel
      in=abs(neq(i))
      if(in.gt.k) k=in
      enddo
      nhlc=k
      do i=nhel+1,nkin
      in=abs(neq(i))
      if(in.gt.k) k=in
      enddo
      nklc=k
      do i=nkin+1,nlgi
      in=abs(neq(i))
      if(in.gt.k) k=in
      enddo
      nlgc=k
      do i=nlgi+1,nvar
      in=abs(neq(i))
      if(in.gt.k) k=in
      enddo
      nvrc=k
c========================================================backbone/thy/rib output
      if(nto.eq.0) goto 500
      write(6,200)
200   format(/2x,'Backbone Locking and Symmetry ....',
     1 //7x,'Base',8x,'  C1''','  C2''','  C3''','  Gly','  C12',
     1 '  C23','  Eps','  Zet','  O3''','  P  ','  Rib','  Thy'/)
      ki=0
      m=ntlg
      mi=nlgi
      do is=1,nto
      ks=(ilq(is,2)-1)*kseq+ilq(is,1)
      ib=nuc(is-1)+iofs(is)
      in=itr(is)
      if(ise(is).lt.0.and.is.ne.1) write(6,210)
210   format(7x,
     1 '-------------------------------------------------------------')
      if(ilq(is,1).eq.nbrk(1).or.ilq(is,1).eq.nbrk(2)) write(6,211)
211   format(7x,
     1 '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
      do l=1,12
      ipr(l)=0
      enddo
      if(lthy(is)) then
      m=m+1
      if(lock(m)) then
      ipr(12)=0
      else
      mi=mi+1
      ipr(12)=neq(mi)
      endif
      endif
      l=0
      do k=ienb(is,1),ienb(is,2)
      l=l+1
      if(.not.lock(k)) then
      ki=ki+1
      ipr(l)=neq(ki)
      endif
      if(ito(is).eq.6) then
      ipr(11)=ipr(7)
      ipr(7)=0
      endif
      enddo
      rb=' '
      if(ribose(is)) rb='R'
      ct=' '
      if(cation(is)) ct='*'
      write(6,220) is,munit(ib),nunit(ib),rb,ct,(ipr(l),l=1,12)
220   format(2x,i3,') ',a4,i3,a1,a1,2x,12i5)
      do k=ienb(is,2)+1,ient(is)
      if(.not.lock(k)) ki=ki+1
      enddo
      enddo
c=====================================================modified nucleotide output
      start=.true.
      k=0
      ki=0
      do is=1,nto
      in=itr(is)
      ino=ito(is)
      io=nuc(is-1)+iofs(is)
      kl=k+kap(1,ino)+1
      ku=k+kap(1,in)
      do l=1,kap(1,ino)
      k=k+1
      if(.not.lock(k)) ki=ki+1
      enddo
      if(ku.ge.kl) then
      kn=ku-kl+1
      do l=1,kn
      ipl(l)=0
      enddo
      do l=1,kn
      k=k+1
      if(.not.lock(k)) then
      ki=ki+1
      ipl(l)=neq(ki)
      endif
      enddo
      if(start) then
      start=.false.
      write(6,201)
201   format(/2x,'Covlig Locking ....',
     1 //7x,'Covlig',8x,'x1',3x,'x2',3x,'x3',3x,'x4',3x,'x5',3x,'x6',
     1 3x,'x7',3x,'x8',3x,'x9',3x,'x10'/)
      endif
      ct=' '
      if(cation(is)) ct='*'
      write(6,202) is,munit(io),nunit(io),ct,(ipl(l),l=1,kn)
202   format(2x,i3,') ',a4,i3,1x,a1,2x,10i5,10(:/18x,10i5))
      endif
      k=k+kap(2,in)-kap(1,in)
      enddo
c============================================================intra ligand output
500   if(nlig.ne.0.and.nbac.gt.nsph) write(6,101)
101   format(/2x,'Ligand Intra Locking ....',
     1 //7x,'Ligand',8x,'L1',3x,'L2',3x,'L3',3x,'L4',3x,'L5',3x,'L6',
     1 3x,'L7',3x,'L8',3x,'L9',3x,'L10'/)
      k=kapt(nto+1)
      ki=nsph
      do il=1,nlig
      is=nto+il
      in=itr(is)
      ib=nuc(is-1)+1
      do i=1,kap(1,in)
      k=k+1
      if(lock(k)) then
      ipl(i)=0
      else
      ki=ki+1
      ipl(i)=neq(ki)
      endif
      enddo
      if(kap(1,in).ne.0) write(6,121) il,munit(ib),nunit(ib),
     1 (ipl(i),i=1,kap(1,in))
121   format(2x,i3,') ',a4,i3,4x,10i5,10(:/18x,10i5))
      enddo
c=================================================================helical output
      if(nto.eq.0) goto 600
      write(6,100)
100   format(/2x,'Helical Locking and Symmetry ....',
     1 //7x,'Base',8x,'Xdis',' Ydis',' Rise','  Inc','  Tip',
     1 ' Twist'/)
      k=ntba
      ki=nbac
      do is=1,nto
        itwl(is)=0
        ib=nuc(is-1)+iofs(is)
        if(ise(is).lt.0.and.is.ne.1) write(6,110)
110   format(7x,'-----------------------------------------')
        if(ilq(is,1).eq.nbrk(1).or.ilq(is,1).eq.nbrk(2)) write(6,112)
112   format(7x,'- - - - - - - - - - - - - - - - - - - - -')
        do l=1,6
          k=k+1
          if(lock(k)) then
            ipr(l)=0
          else
            ki=ki+1
            ipr(l)=neq(ki)
            if(l.eq.6) itwl(is)=ki
          endif
        enddo
        rb=' '
        if(ribose(is)) rb='R'
        ct=' '
        if(cation(is)) ct='*'
        write(6,120) is,munit(ib),nunit(ib),rb,ct,(ipr(l),l=1,6)
120     format(2x,i3,') ',a4,i3,a1,a1,2x,6i5)
      enddo
c====================================================================kink output
      k=nthe
      ki=nhel
      write(6,400)
400   format(/2x,'Kink variable locking ....',
     1 //7x,'Base',8x,'  Ax ','  Ay ',' Ainc',' Atip'/)
      do is=2,kseq
      ib=nuc(is-1)+iofs(is)
      do l=1,4
      ipr(l)=0
      enddo
      if(kink(is)) then
      do l=1,4
      k=k+1
      if(.not.lock(k)) then
      ki=ki+1
      ipr(l)=neq(ki)
      endif
      enddo
      endif
      if(is.eq.nbrk(1).or.is.eq.nbrk(2)) write(6,411)
411   format(7x,'- - - - - - - - - - - - - - - -')
      write(6,420) is,munit(ib),nunit(ib),(ipr(l),l=1,4)
420   format(2x,i3,') ',a4,i3,4x,4i5)
      enddo
c============================================================superhelical output
      if(isur.ne.0) then
      isu1=isur
      if(isur.gt.0) isu1=neq(isur)
      isu2=isup
      if(isup.gt.0) isu2=neq(isup)
      write(6,401) isu1,isu2
401   format(/2x,'Superhelix ....   Rad  Pit',//18x,2i5)
      endif
c============================================================inter ligand output
600   if(nlig.ne.0) write(6,102)
102   format(/2x,'Ligand Inter Locking ....',
     1 //7x,'Ligand',8x,'Xdisp',2x,'Ydisp',3x,'Rise',4x,'Inc',
     1 4x,'Tip',2x,'Twist'/)
      do il=1,nlig
        is=nto+il
        ib=nuc(is-1)+1
        kal=nuc(is)-nuc(is-1)
        jlim=6
        if(kal.eq.1) jlim=3
        do j=1,jlim
          k=k+1
          if(lock(k)) then
            ipr(j)=0
          else
            ki=ki+1
            ipr(j)=neq(ki)
          endif
        enddo
        write(6,122) il,munit(ib),nunit(ib),(ipr(j),j=1,jlim)
122   format(2x,i3,') ',a4,i3,4x,6i7)
      enddo
c=================================================================symmetry check
      do i=1,nvrc
        ipr(i)=0
      enddo
      start=.true.
	do i=1,n7
	   chec(i) = 0.0
	enddo
      do i=1,nvar
        ia=neq(i)
        if(ia.gt.0) chec(ia)=var(i)
      enddo
      do i=1,nvar
        fhel=.false.
        if(isur.ne.0) then
          do j=1,nto
             if(itwl(j).eq.i) fhel=.true.
          enddo
        endif
        ia=abs(neq(i))
        if(ia.gt.0) then
          if(abs(var(i)-chec(ia)).gt.5.d-2) then
            if(fhel) then
               if(abs(abs(var(i))+abs(chec(ia))-360).le.5.d-2) goto 160
            endif
            if(start) write(6,*)
            start=.false.
            if(.not.autos) then
               if(ipr(ia).eq.0) write(6,150) ia,var(i),chec(ia)
150            format(2x,'Sym group ',i3,' error: ',f8.3,
     1              ' / ',f8.3)
               istop=1
            endif
            ipr(ia)=1
          endif
        endif
160   enddo
c===============================================================symmetry forcing
      if(autos) then
      start=.true.
      do i=1,nvrc
      nvs(i)=0
      chec(i)=0.
      ches(i)=0.
      enddo
      do i=1,nvar
      angle=.true.
      if(i.gt.nbac.and.i.le.nlgi.and..not.lar(i)) angle=.false.
      ia=abs(neq(i))
      if(.not.angle) then
      chec(ia)=chec(ia)+var(i)
      else
      chec(ia)=chec(ia)+cos(var(i)*cdr)
      ches(ia)=ches(ia)+sin(var(i)*cdr)
      endif
      nvs(ia)=nvs(ia)+1
      enddo
      do i=1,nvar
      angle=.true.
      if(i.gt.nbac.and.i.le.nlgi.and..not.lar(i)) angle=.false.
      ia=abs(neq(i))
      if(.not.angle) then
      del=chec(ia)/nvs(ia)-var(i)
      var(i)=var(i)+del
      else
      r=sqrt(ches(ia)**2+chec(ia)**2)
      del=acos(chec(ia)/r)*crd-var(i)
      if(ches(ia).lt.0) del=-acos(chec(ia)/r)*crd-var(i)
      var(i)=var(i)+del
      endif
      if(neq(i).gt.0.and.ipr(ia).ne.0) then
      if(start) write(6,*)
      start=.false.
      write(6,190) ia,var(i),del
190   format(2x,'Sym group ',i3,' set to ',f8.3,'  del= ',f8.3)
      endif
      enddo
         if(isur.ne.0) then
         do k=1,nst
         ks=ksym(1)+1
         dw=hel(ks,6)-hel(1,6)
         if(abs(dw).gt.180) dw=dw-sign(360.d0,dw)
         dw=dw/(ks-1)
         do i=1,kseq
         if(i.le.ks) then
         win(i)=hel(i,6)-(i-1)*dw
         else
         win(i)=win(i-ks+1)
         endif
         l=itwl(ieq(i,k))
         if(l.gt.0) var(l)=win(i)
         enddo
         write(6,300) k,-dw
300      format(2x,'Sym suphe ',i3,' reset  ',f8.3,' (Twist)')
         enddo
         endif
c--------------------------------------------------exocyclic dihedral protection
      k=0
      ki=0
      do is=1,ntl
        in=itr(is)
        do l=1,kap(1,in)
          k=k+1
          if(.not.lock(k)) then
            ki=ki+1
            set(is,l)=var(ki)
          endif
        enddo
        k=k+kap(2,in)-kap(1,in)
      enddo
      call setbac
      endif
c===============================================================================
      if(istop.eq.1) iret = 1
      return
      end
