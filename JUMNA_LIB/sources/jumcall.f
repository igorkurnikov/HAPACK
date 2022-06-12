      subroutine jumcall
c*******************************************************************************
c**                                                                           **
c**  Junction Minimization of Nucleic Acids 10.1               R.L. Apr 1996  **
c**                                                            HP Version     **
c**  
c**                                                                           **
c*******************************************************************************
      include 'jumna_data.inc'
      logical*2 sup,rcom,homo,homo2,homo3,diep,sum,link,ecen,cyl,
     1 lar,lock,ifhb,ribose,recep,donor,kink,kint,hst,bst,vst,lgi,fabs,
     1 lgj,lgu,lcat,cent,autos,locr,lthy,mods,rst(n9),cation,avec(53),
     1 start,tlp,quiet,convg,ihl,modif(n2),dloop,amber
      logical*4 there
      character*4 mnam,munit,snam,suni,sub,symb*1,lnam,pab(25,25)*1,
     1 mac*32,lmo*32,lib*32,libm*32,out*32,axe*32,noe*32,nol*32,axl*32,
     1 test*32,pdb*32,ins*32,bar*32,parm*32,anam(n9,6),name(n2),
     1 knam,code*8,kode*8,dm*9,vcyl(0:1)*6,seq*120,line*100,
     1 atclas*2,atclasn*2,strand(4)*3,dirn(-1:1)*5,tran(15)*1,tram(15)*1
	character*80 libn
	character*132 out_axx
      integer*2 i23,i34,elim
      integer*4 opt,nri(n2)
      real*4 dtime,tarray(2)
      dimension kuc(n2),ifr(n1),ipl(n9),helt(n2,6),sett(n2,10),vkint(4),vdw(25)
	common/cha/mac,lmo,lib,libm,out,axe,axl,noe,nol,test,pdb,ins,
     1 bar,parm
	common/jmdirs/libn
      common/cnd/rcyl,dcyl,hcyl,tcyl,eneq,ucyl(3),pcyl(3),ecyl,
     1 ac(25),bc(25),lcyl(3)
      common/comjm/cut,icx(n2),icn(n2),ict(n2,10),ncut,nmx,mco,ncomb,icomb
      common/cur/fb1x,fb1y,fb1z,fb2x,fb2y,fb2z,fbx,fby,fbz,ib1,ib2
      common/datjm/acc,phos,epsi,epsr,plat,slope,rhbl,vfac,tfac,rfac,xfac,
     1 scale,rpiv,tpiv,fad,fan,damp,fnoew,fnoes,fnoea,df1,scnb,scee,
     1 catd,catr,catc,rad,pit,enit,nshel,maxn,opt,mhomo,mhomo2,mhomo3,
     1 nick,limit,nop,nrib,ncat,lig,naxo,sup,rcom,homo,homo2,homo3,diep,
     1 sum,link,ecen,cyl,lcat,cent,autos,amber
      common/dcu/dpx,dpy,dpz,dux,duy,duz,crx,cry,crz,ctx,cty,ctz,
     1 csx,csy,csz,crs,clx,cly,crl
      common/dob/theta(20),ind(20,3),nba(20)
      common/enf/for(n1,3),tor(n1,3),fot(n6),ener,elec,
     1 repl,disp,eang,etog,epen,i23(8*n1),i34(8*n1),elim(n1)
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/flx/sap(n6),refg,refb,refh,refv,iap(n8,n4,n5),
     1 nap(n8,7,n5),kap(3,n5)
      common/hel/ua(n2,3),da(n2,3),ra(n2,3)
      common/ind/iend(0:4),ise(n2),ienb(n2,2),kapt(n2+1),nsph
      common/lgd/rlig(n9,6),slig(n9,2),ilig(n9,6),lopt(n9),lpiv(n9),
     1 ntlg,nlgi,nlig,ntl,lnam(n9),locr(n2,n0)
      common/lop/dellp(2),stlp(2),gmx,indlp(2,2),lplow(2),lphig(2),
     1 nloop,icy,tlp(2),quiet,convg,dloop
      common/moljm/sor(n4,n5,3),smon(n4,n5),snam(n4,n5),suni(n4,n5),
     1 nuni(n4,n5),sub(n5),isch(n4,n5),isty(n4,n5),ics(n4,n5),
     1 mats(3*n4,n5),kas(n5),khs(n5),ksub
      common/mnn/var(n7),gra(n7),con(n0*n2+n3),scl(n7),ncon,nvrc,
     1 ntba,nbac,nthe,nhel,ntki,nkin,ntri,ntot,nvar,nrin,
     1 lar(n7),lock(n6a)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/nmrjm/rnoe(n3),bnoe(n3),fnoe(n3),snoe,knam(200,2),ktyp(200),
     1 icon(n2*2),inoe(n3,9),ih68(n2),jnoe(n3),nnoe,ndcs,ndiv(0:99),nlin
      common/parjm/aij(25,25),bij(25,25),ahb(15),bhb(15),vt(35),va(25),
     1 vo(25),delmax,time0,nzsh,nwdg,ivf(35),ifhb(25,25)
      common/para/ro(57),reps(57),angc(191),fangc(191),btor(88),
     $ gtor(88),ftor(88),bdtor(30),gdtor(30),fdtor(30),amb(57,57),
     $ bmb(57,57),amh(9,11),bmh(9,11),mtor(88),itor(88,4),iadt(30,4),
     $ iangc(191,3),mof,nof,nang,ntor,ntors,nad,atclas(57)
      common/paran/ron(53),repsn(53),angcn(188),fangcn(188),btorn(88),
     $ gtorn(88),ftorn(88),bdtorn(30),gdtorn(30),fdtorn(30),
     $ mtorn(88),itorn(88,4),iadtn(30,4),iangcn(188,3),mofn,nofn,nangn,
     $ ntorn,ntorsn,nadn,atclasn(53)
      common/slg/hst(n2,6),bst(n2,n8),vst(n2,4),lgi(n1),lgj(n1),lgu(n2)
      common/srs/corms(n1,3),saps(n6),vars(n7),gras(n7),hels(n2,6),
     1 vkis(n2,4),has(n2,9),rlis(n9,6),eref,rref,pref
      common/sst/phase(8,3),ampli(8,3),ph5(n2),am5(n2),kr5(n2)
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(120,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      common/symjm/xsum,ysum,itwl(n2),nvs(n7),ihm(n2),ksym(3),nbrk(2),
     1 isym,isur,isup,ihl(n2)
	common/hacon/iret
      data dm/'    -----'/,vcyl/'Locked','Free'/,
     1 strand/'1st','2nd','3rd','4th'/,dirn/'3''-5''',' ','5''-3'''/
      data vdw/1.2,1.2,1.6,1.6,1.5,1.5,1.5,1.4,1.4,1.4,1.9,1.2,1.2,
     &         0.95,0.65,1.33,0.99,1.85,1.85,0.0,1.35,1.8,3*0.0/
      data tran/'G','A','I','C','T','U',
     1          'H','B','J','D','S','V','O','Q','R'/
      data tram/'g','a','i','c','t','u',
     1          'h','b','j','d','s','v','o','q','r'/
      data refo,refp/119.13,101.39/
c-------------------------------------------------------------------default data
      iret = 0
	kk = 0
      if( kk.eq.20) then
      call flex
	call ampar91
	call ampar94
	endif
c	write(6,*) ' aij(1,1) = ',aij(1,1)
      cut=-1.0
      ncut=0
      mac=' '
      lmo=' '
      libm=' '
      out='jumnatest'
      axe=' '
      axl=' '
      noe=' '
      nol=' '
      test=' '
      pdb='jumnatest'
      ins=' '
      bar=' '
      sup=.true.
c      ecen=.false.
c      ecen = .true.
      cyl=.false.
      lcyl(1)=0
      lcyl(2)=0
      lcyl(3)=0
      rcom=.false.
c      rad=0.
c      pit=0.
      nshel=0
c      homo =.false.
c      homo = .true.
c      homo2=.false.
c      homo3=.false.
c      mhomo =-1
c      mhomo2=-1
c      mhomo3=-1
      nick=0
      nrib=0
      ncat=0
      lig=0
      naxo=0
      vfac=0.3
      tfac=3.0
      rfac=1.5
      xfac=0.1
      scale=0.1
      damp=1.
      fad=200.0
      fan=400.0
      limit=0
c      opt= 0
c      maxn=500
      acc=1.d-4
      df1=10.d0
      probe=0.0
      diep=.true.
      epsi=1.0
      epsr=0.0
c      slope=0.16
      plat=78.
      enit=1.
      rhbl=7.
c      phos=-0.5
      link=.true.
      sum=.false.
      rpiv=2.5
      tpiv=0.0
      nop=0
      fnoew=6.
      fnoes=12.
      fnoea=1000.
      nnoe=0.
      lcat=.false.
      catd=2.53
      catr=0.0
      catc=1.0
      cent=.false.
      autos=.false.
      nloop=0
      dloop=.false.
      quiet=.false.
      nzsh=0
      nwdg=0
      ib1=0
      ib2=0
      isur=0
      isup=0
      amber=.false.
c      parm='Flex'
      scnb=2.0
      scee=1.2
      snoe=10000.
      kfi=index(libn,' ')-1
      if(kfi.lt.0) then
        libn='.'
        kfi=1
      endif
      if(parm.eq.'Flex') then
         libn=libn(:kfi)//'/nucave'
      else if(parm.eq.'Amber91') then
         amber=.true.
         libn=libn(:kfi)//'/nucamb91'
         if(nrib.eq.-1) libn=libn(:kfi)//'/ribamb91'
         if(nrib.gt.0) then
           write(6,*) '  ---- DNA/RNA hybrids not allowed for parm91 ----'
           return
         endif
      else if(parm.eq.'Amber94') then
         amber=.true.
         libn=libn(:kfi)//'/nucamb94'
         if(nrib.eq.-1) libn=libn(:kfi)//'/ribamb94'
         if(nrib.gt.0) then
           write(6,*) '  ---- DNA/RNA hybrids not allowed for parm94 ----'
           return
         endif
      else
         write(6,*) '  ---- Error in ''parm'' ---'
         return
      endif
      fabs=.false.
      if(axe.ne.' ') then
      kfi=index(axe,' ')-1
      inquire(file=axe(:kfi)//'.axe',exist=there)
      if(.not.there) fabs=.true.
      endif
      if(noe.ne.' ') then
      kfi=index(noe,' ')-1
      inquire(file=noe(:kfi)//'.noe',exist=there)
      if(.not.there) fabs=.true.
      endif
      if(axl.ne.' ') then
      kfi=index(axl,' ')-1
      inquire(file=axl(:kfi)//'.axl',exist=there)
      if(.not.there) fabs=.true.
      endif
      if(nol.ne.' ') then
      kfi=index(nol,' ')-1
      inquire(file=nol(:kfi)//'.nol',exist=there)
      if(.not.there) fabs=.true.
      endif
         if(fabs) then
         write(6,*) '  ---- Missing I/P file : Check namelist ----'
         return
         endif
c----------------------------------------------------------first axl file to axe
         if(axl.ne.' ') then
            if(nol.ne.' ') then
            write(6,*) '  ---- Choose either Axl or Nol ! ----'
            return
            endif
         kfi=index(axl,' ')-1
         open(unit=1,file=axl(:kfi)//'.axl',status='old')
         read(1,51,end=100) line
         kfe=index(line,'axe')-2
         axe=line(:kfe)
         close(1)
         endif
c---------------------------------------------------------------setup parameters
      if(.not.amber) then
      do i=1,25
      recep=.false.
      donor=.false.
      if((i.ge.7.and.i.le.10).or.i.eq.18) recep=.true.
      if(i.le.2) donor=.true.
      do j=1,i
      pab(i,j)='-'
      if(aij(i,j).ne.0.and.bij(i,j).ne.0) pab(i,j)='*'
      pab(j,i)=pab(i,j)
      aij(j,i)=aij(i,j)
      bij(i,j)=bij(i,j)*1000
      bij(j,i)=bij(i,j)
      ifhb(i,j)=.false.
      if(recep.and.j.le.2) ifhb(i,j)=.true.
      if(donor.and.((j.ge.7.and.j.le.10).or.j.eq.18)) ifhb(i,j)=.true.
      ifhb(j,i)=ifhb(i,j)
      enddo
      enddo
      do i=1,29
      do j=21,n4
      do k=1,6
      iap(i,j,k)=0
      enddo
      enddo
      enddo
      else
c------------------------------------------------------------------------- amber
      if(parm.eq.'Amber94') then
        do i=1,53
          ro(i)=ron(i)
          reps(i)=repsn(i)
          atclas(i)=atclasn(i)
        enddo
        do i=1,nangn
          angc(i)=angcn(i)
          fangc(i)=fangcn(i)
          do j=1,3
            iangc(i,j)=iangcn(i,j)
          enddo
        enddo
        do i=1,ntorn
          btor(i)=btorn(i)
          gtor(i)=gtorn(i)
          ftor(i)=ftorn(i)
          mtor(i)=mtorn(i)
          do j=1,4
            itor(i,j)=itorn(i,j)
          enddo
        enddo
        do i=1,nadn
          bdtor(i)=bdtorn(i)
          gdtor(i)=gdtorn(i)
          fdtor(i)=fdtorn(i)
          do j=1,4
            iadt(i,j)=iadtn(i,j)
          enddo
        enddo
        nang=nangn
        ntor=ntorn
        ntors=ntorsn
        nad=nadn
      endif
c----------------------------------------------------------------------amb,bmb
        do i=1,53
          do j=1,53
            ee=sqrt(reps(i)*reps(j))
            r6=(ro(i)+ro(j))**6
            amb(i,j)=ee*r6*r6
            bmb(i,j)=2.*ee*r6
          enddo
        enddo
      endif
c------------------------------------------------------- end of amber parameters
      time0=dtime(tarray)
c      read(51,*)  (ksym(i),i=1,3),nbrk(1),nbrk(2)
      write(6,1) (ksym(i),i=1,3),nbrk(1),nbrk(2)
1     format(/2x,'Symmetry : ',3i4,
     1       /2x,'Break pts: ',2i4)
      isym=n2
      if(nbrk(1)+nbrk(2).eq.0.and.ksym(1).gt.0) isym=ksym(1)
      if(rad.ne.0) isur=-1
      if(pit.ne.0) isup=-1
c------------------------------------------------------------cation modification
      do j=1,25
      aa=aij(14,j)
      bb=bij(14,j)
      if(aa.ne.0..and.bb.ne.0.) then
      req=(2*bb/aa)**(1./6.)+catr
      eeq=-aa*aa/(4.*bb)
      aij(14,j)=-2*eeq*req**6
      bij(14,j)=-eeq*req**12
      aij(j,14)=aij(14,j)
      bij(j,14)=bij(14,j)
      endif
      enddo
      if(lcat.and.ncat.eq.0) then
      lcat=.false.
      write(6,*) '  ----Warning: Lcat set false because ncat=0 ----'
      endif
      fad=fad*damp
      fan=fan*damp
      nlig=lig
c-----------------------------force field parameters: torsion and valence angles
      if(bar.ne.' ') then
      kfi=index(bar,' ')-1
      write(6,112) bar(:kfi)//'.bar'
112   format(/2x,'New force field parameters from ',a,/)
      open(unit=1,file=bar(:kfi)//'.bar',status='old')
      read(1,*) (vt(i),i=1,35)
      read(1,*) (vo(i),i=1,15)
      read(1,*) (va(i),i=1,15)
      close(1)
      endif
c-------------------------------------------------------------------------------
      if(opt.lt.-1.or.opt.gt.1) then
      write(6,*) '  ---- ILLEGAL OPT VALUE ----'
      return
      endif
c      read(51,*) nst,(idr(j),j=1,4)
      kseq=idr(1)
      if(kseq.lt.0) then
      write(6,*) '  ---- FIRST STRAND MUST BE 5''-3'' ----'
      return
      endif
      ncheck=0
      do j=1,4
      if(idr(j).eq.0) goto 100
      ncheck=ncheck+1
      enddo
100   if(ncheck.ne.nst) then
      write(6,*) '  ---- ERROR IN 1ST DATA LINE ----'
      return
      endif
      nto=0
      do j=1,nst
      nto=nto+abs(idr(j))
      if(abs(idr(j)).gt.kseq) then
      write(6,*) '  ---- FIRST STRAND MUST BE THE LONGEST ----'
      return
      endif
      enddo
      ntl=nto+nlig
      if(nlig.gt.n9) then
      write(6,*) '  ---- n9 too small ----'
      return
      endif
      if(nto.gt.n2) then
      write(6,*) '  ---- n2 too small ----'
      return
      endif
      if(ecen.and.nlig.ne.0) then
      write(6,*) '  ---- NO LIGANDS ALLOWED WITH ECEN ----'
      return
      endif
      if(epsr.ne.0.and.diep) then
      write(6,*) '  ---- Warning EPSR ignored because DIEP true ----'
      return
      endif
      k=0
      jl=1
      ju=kseq
      do i=1,nst
c        read(51,5) (seq(j:j),j=jl,ju)
5       format(50a)
        do j=jl,ju
         l=j-jl+1
         if(seq(j:j).eq.'-') then
           ieq(l,i)=0
         else
           k=k+1
           ieq(l,i)=k
           ilq(k,1)=l
           ilq(k,2)=i
         endif
        enddo
        jl=jl+kseq
        ju=ju+kseq
      enddo

      if(k.ne.nto) then
      write(6,*) '  ---- NTO INCONSISTANT WITH I/P ----'
      return
      endif
      do i=1,kseq
      kink(i)=.false.
      kode(i)='--------'
      do j=1,4
      vkink(i,j)=0.
      enddo
      enddo
c----------------------------------------------------------------------find ends
      do i=0,4
        iend(i)=0
      enddo
      do i=1,nto
        ise(i)=0
      enddo
      il=1
      iu=0
c il - index of the first residue in the strand
c iu - index of the last residue in the strand
      do k=1,nst
        iu=iu+abs(idr(k))
        iend(k)=iu
        if(idr(k).gt.0) then
           ise(il)=-5
           ise(iu)=3
        else
           ise(il)=-3
           ise(iu)=5
        endif
        il=iu+1
      enddo
c--------------------------------------------------------------------------nicks
      if(nick.ne.0) then
      if(nick.lt.0) then
      write(6,*)  '  ---- Warning: Nick must be >0 ----'
      return
      endif
      if(ise(nick).ne.0.or.ise(nick-1).ne.0) then
      write(6,*)  '  ---- Warning: Nick includes terminal unit ----'
      return
      endif
      endif
c-----------------------------------------------------------itr for nucleic acid
c indicator of presence of modified bases
      mods=.false.  
      do is=1,ntl
        irec(is)=0
        lthy(is)=.false.
      enddo
      k=0
      lr=0
c lr - absolute index of the residue
c ir - index of the residue in the strand
c name(lr) - array of processed names with prefix X for 5',Y for 3' and Z for in strand, 
c            postfix 'r' - for ribonucleotide(if not Flex force field)
      koff=6
      do kr=1,nst
        idir=sign(1,idr(kr))
        do ir=1,kseq
          k=k+1
          is=ieq(ir,kr)
c if the position in the strand is empty get to the next position
          if(is.eq.0) goto 20   
          lr=lr+1
          itr(lr)=0
          do l=1,15
c standard bases
             if(seq(k:k).eq.tran(l)) then
                if(seq(k:k).eq.'T'.or.seq(k:k).eq.'t'.or.
     1             seq(k:k).eq.'S'.or.seq(k:k).eq.'s') lthy(lr)=.true.
                if(abs(ise(is)).eq.5) then
                   name(lr)='X'//tran(l)
                   itr(lr)=1
                   ito(lr)=1
                else if(abs(ise(is)).eq.3) then
                   name(lr)='Z'//tran(l)
                   itr(lr)=3
                   ito(lr)=3
                else
                   name(lr)='Y'//tran(l)
                   itr(lr)=2
                   ito(lr)=2
                endif
             else if(seq(k:k).eq.tram(l)) then
c modified bases
                if(seq(k:k).eq.'T'.or.seq(k:k).eq.'t'.or.
     1             seq(k:k).eq.'S'.or.seq(k:k).eq.'s') lthy(lr)=.true.
                mods=.true.
                if(abs(ise(is)).eq.5) then
                   name(lr)='X'//tram(l)
                else if(abs(ise(is)).eq.3) then
                   name(lr)='Z'//tram(l)
                else
                   name(lr)='Y'//tram(l)
                endif
                do m=1,lr-1
                   if(name(lr).eq.name(m)) then
                       itr(lr)=itr(m)
                       goto 64
                   endif
                enddo
                koff=koff+1
                itr(lr)=koff
             endif
64        enddo
          if(itr(lr).eq.0) then
             write(6,65) seq(k:k),kr,ir
65           format(/2x,'Nucleotide code ',a1,' in strand= ',i2,' at posn= ',
     1              i2,' was not translated'/)
             return
          endif
20      enddo
      enddo
c=========================================================optional i/p: cylinder
      if(cyl) read(51,*) lcyl(1),lcyl(2),lcyl(3),tcyl,dcyl,hcyl,rcyl,eneq
c-------------------------------------------------------------------------ribose
      do is=1,nto
        ribose(is)=.false.
      enddo
c if set all residues to ribonucleotides:
      if(nrib.lt.0) then
        do is=1,nto
          if(itr(is).le.6) ribose(is)=.true.
        enddo
      else if(nrib.gt.0) then
        read(51,*) (nri(i),i=1,nrib)
        do i=1,nrib
           in=nri(i)
           if(itr(in).gt.6) then
              write(6,*) '  ---- RIBOSE NOT ALLOWED FOR MOD_NUC, USE NCHEM ----'
              return
           endif
           ribose(in)=.true.
           if(parm.ne.'Flex') name(in)=name(in)(:2)//'r'
        enddo
      endif
c-------------------------------------------------------------------------cation
      do is=1,nto
        cation(is)=.false.
      enddo
      if(ncat.lt.0) then
        do is=1,nto
          if(abs(ise(is)).ne.3) cation(is)=.true.
        enddo
      else if(ncat.gt.0) then
        read(51,*) (nri(i),i=1,ncat)
        do i=1,ncat
          is=nri(i)
          if(is.gt.nto) then
             write(6,*) '  ---- NO CATIONS ALLOWED ON LIGANDS ----'
             return
          endif
          if(abs(ise(is)).eq.3) then
             write(6,*) ' ---- NO CATIONS ALLOWED ON 3'' END (no phos) ----'
             return
          else
             cation(is)=.true.
          endif
        enddo
      endif
c---------------------------------------------------------------------open bases
      if(nop.gt.0) then
        npr=0
        do i=1,nop
          read(51,*) nba(i),theta(i)
          if(nba(i).lt.npr) then
             write(6,*) '  ---- GIVE OPEN BASES IN INCR. ORDER ----'
             return
          endif
          npr=nba(i)
          if(ise(npr).ne.0) then
             write(6,*) '  ---- NO OPENING ALLOWED FOR TERMINAL BASES ----'
             return
          endif
        enddo
      endif
c-------------------------------------------------------------conformation input
c      k=0
c      do i=1,nto
c         if(k.gt.0) goto 55
c         read(51,5) line
c         k=index(line,'*')
c         if(k.gt.0) read(line(k+1:),*) k
c55       call liner(i,line)
c      if(index(code(i),'K').ne.0) then
c         if(i.eq.1.or.i.gt.kseq) then
c            write(6,*) '  ---- KINK BADLY PLACED ----'
c            return
c         endif
c         kink(i)=.true.
c      endif
c         if(kink(i)) call kline(i)
c         k=k-1
c      enddo
c--------------------------------------------------------non-bonded ligand input
      do il=1,nlig
      is=nto+il
      read(51,51) line
      code(is)=line(:8)
      k=index(line(10:),' ')+9
      lnam(il)=line(10:k-1)
      read(line(k:),*) lopt(il),ilig(il,1)
      if(lopt(il).eq.1) then
      if(ilig(il,1).gt.0) then
      read(line(k:),*) lopt(il),ilig(il,1),(rlig(il,j),j=1,6)
      else
      read(line(k:),*) lopt(il),ilig(il,1),(rlig(il,j),j=1,6),
     1 slig(il,1),slig(il,2)
      endif
      else if(lopt(il).eq.2) then
      read(line(k:),*) lopt(il),ipl(il),(anam(il,j),j=1,6),
     1 (rlig(il,j),j=1,6)
      endif
      if(lopt(il).eq.0) then
      rst(il)=.false.
      else
      rst(il)=.true.
      endif
      enddo
      call molin(name,lr,koff,mods)
	if(iret.eq.1)then
	  return
	endif
      if( il.eq.9999898) then
        call flex
	  call ampar91
	  call ampar94
	endif
c======================================================================axe input
      if(axe.ne.' ') then
        kfi=index(axe,' ')-1
        open(unit=7,file=axe(:kfi)//'.axe',status='old')
        read(7,51) line
51      format(a)
c------------------------------------------Abasic Grenoble
         iabas=0
         if(index(line,'*').ne.0) then
         do is=1,nto
           if(name(is).eq.'YO') iabas=is
         enddo
       endif
c---------------------------------------------------------
      if(index(line,'#').eq.0) then
      read(line,53) ns,id1,id2,id3,id4,rado,pito
53    format(5i4,2f8.3)
         if(rad.eq.0.and.rado.ne.0) then
         rad=rado
         isur=-1
         endif
         if(pit.eq.0.and.pito.ne.0) then
         pit=pito
         isup=-1
         endif
      else
      read(line,*) ns,id1,id2,id3,id4
      endif
      ntt=abs(id1)+abs(id2)+abs(id3)+abs(id4)
         if(nto.eq.0.and.ntt.gt.0.and.id1.eq.nlig) then
         do i=1,ntt
         read(7,52) slig(i,1),slig(i,2)
52       format(20x,f9.4,18x,f9.4,/)
         if(i.gt.1) then
         slig(i,1)=slig(i,1)+slig(i-1,1)
         slig(i,2)=slig(i,2)+slig(i-1,2)
         endif
         enddo
         ltt=ntt/7
         if(ltt*7.lt.ntt) ltt=ltt+1
         do i=1,2*ltt+id1-1
         read(7,*)
         enddo
         goto 700
         endif
         if(ns.ne.nst.or.id1.ne.idr(1).or.id2.ne.idr(2).or.
     1       id3.ne.idr(3).or.id4.ne.idr(4)) then
            write(6,*) '  ---- INCOMPATIBLE .axe ----'
            return
         endif
      do is=1,nto
        read(7,50) (helt(is,j),j=1,6),sett(is,9),sett(is,10),
     1             (sett(is,j),j=1,8)
        if(sett(is,9).le.0)  sett(is,9) =refo
        if(sett(is,10).le.0) sett(is,10)=refp
        helt(is,1)=-helt(is,1)
        helt(is,5)=-helt(is,5)
      enddo
      do is=1,nto
        do j=1,6
           if(.not.hst(is,j)) hel(is,j)=helt(is,j)
        enddo
        do j=1,8
           if(.not.bst(is,j)) set(is,j)=sett(is,j)
        enddo
        if(abs(ise(is)).ne.3) then
           set(is,9)=sett(is,9)
           set(is,10)=sett(is,10)
        endif
      enddo
50    format(2x,8f9.4,/2x,8f9.4)
      if(nto.gt.0) read(7,60) (thy(is),is=1,nto)
      if(nto.gt.0) read(7,60) (sett(is,1),is=1,nto)
60    format(2x,7f9.4)
      do is=1,nto
      if(abs(ise(is)).ne.3) then
      if(.not.bst(is,9)) set(is,11)=sett(is,1)
      else
      if(.not.bst(is,7)) set(is,7)=sett(is,1)
      endif
      enddo
      do is=2,kseq
      read(7,70) (vkint(j),j=1,4),kint
70    format(2x,4f9.4,l2)
      vkint(1)=-vkint(1)
      vkint(4)=-vkint(4)
      do j=1,4
      if(.not.vst(is,j)) vkink(is,j)=vkint(j)
      enddo
      if(.not.kink(is)) kink(is)=kint
      enddo
c-----------------------------------------------------------modified nucleotides
      do is=1,nto
        in=itr(is)
        ino=ito(is)
        idel=kap(1,in)-kap(1,ino)
        if(idel.gt.0) read(7,60) (set(is,j),j=kap(1,ino)+1,kap(1,in))
        if(kap(2,in).gt.n8) then
          write(6,*) '  ---- n8 too small ----'
          return
        endif
      enddo
c-------------------------------------------------------------non-bonded ligands
700   do il=1,nlig
      is=nto+il
      in=itr(is)
      koff=kapt(is)
      if(nto.eq.0.and.ntt.eq.0) then
      read(7,23,end=24) (sett(1,j),j=1,6),ilis,slig(il,1),slig(il,2)
      else
      read(7,23,end=24) (sett(1,j),j=1,6),ilis
      endif
23    format(2x,6f9.4,i3,2f9.4)
      if(.not.rst(il)) then
      lopt(il)=1
      if(nto.gt.0) ilig(il,1)=ilis
      do j=1,6
      rlig(il,j)=sett(1,j)
      enddo
      endif
      if(kap(1,in).gt.0) read(7,60) (set(is,j),j=1,kap(1,in))
      goto 26
24    do j=1,kap(1,in)
      set(is,j)=999.
      enddo
26    enddo
      close(7)
      endif
c----------------------------------------------------superhelical variable check
      if(nshel.gt.0.and.rad.eq.0) then
        write(6,*) '  ---- nshel allowed only for non-zero rad ----'
        return
      endif
      if(isur.ne.0) cent=.true.
c------------------------------------------------check non-covalent ligand links
      do il=1,nlig
      if(ilig(il,1).lt.0.or.ilig(il,1).gt.kseq) then
      write(6,*) '  ----LINK LIGAND TO 1ST STRAND NUCLEOTIDE ----'
      return
      endif
      enddo
c-------------------------------------------set/lock covlig and ligand variables
      k=0
      do is=1,ntl
      kuc(is)=k
      in=itr(is)
      ino=ito(is)
      modif(is)=.false.
      if(in.ne.ino.or.is.gt.nto) modif(is)=.true.
      do l=1,kap(2,in)
      k=k+1
      lock(k)=.false.
      enddo
      enddo
c58    read(51,51,end=500) line
c      symb=line(:1)
c      if(symb.eq.'S') then
c      read(line(2:),*) is,k,val
c      set(is,k)=val
c      bst(is,k)=.true.
c      else if(symb.eq.'L') then
c      read(line(2:),*) is,k
c      ino=ito(is)
c      if(is.le.nto.and.k.le.kap(1,ino)) then
c      write(6,57) is,k
c57    format(/2x,'---- Variable ',i2,' in nuc ',i2,
c     1 ' is not in covlig ----'/)
c      return
c      endif
c      lock(kuc(is)+k)=.true.
c      else if(symb.eq.'R') then
c      read(line(2:),*) is,k
c      if(k+1.gt.nsr(is)) then
c      write(6,97) is,k
c97    format(/2x,'---- Residue ',i2,' has no ring number ',i1,' ----'/)
c      return
c      endif
c     locr(is,k+1)=.true.
c      else if(symb.eq.'A') then
c      read(line(2:),*) is
c      in=itr(is)
c      if(is.le.nto) then
c      ino=ito(is)
c      kl=kap(1,ino)+1
c      else
c      kl=1
c      endif
c      do l=kl,kap(1,in)
c      lock(kuc(is)+l)=.true.
c      enddo
c      else
c      write(6,56) line(:1)
c56    format(/2x,'---- ',a1,' not known for covlig/ligand input ----'/)
c      return
c      endif
c      if(.not.modif(is)) then
c      write(6,59) is
c59    format(/2x,'---- Error in nuc no. ',i2,
c     1 ' for covlig/ligand data ----'/)
c      return
c      endif
c      goto 58
c--------------------------------------------------------------------------build
500   call build(ctot)
      if(iret.eq.1) return
      if(.not.amber) then
        do i=1,25
          avec(i)=.false.
        enddo
        do i=1,kam
          avec(imty(i))=.true.
        enddo
        do i=1,25
          do j=i,25
             if(i.eq.20.or.j.eq.20) goto 502
             if(pab(i,j).eq.'-'.and.(avec(i).and.avec(j))) then
                 write(6,501) i,j
501              format(/2x,'---- LJ params do not exist for classes ',i2,
     1                  '/',i2,' ----')
                 return
              endif
502       enddo
        enddo
      endif
      do il=1,nlig
      if(lopt(il).eq.2) then
      ip=ipl(il)
      is=nto+il
      kal=nuc(is)-nuc(is-1)
      jlim=6
      if(kal.le.2) jlim=3
      do j=1,3
      do i=nuc(ip-1)+1,nuc(ip)
      if(anam(il,j).eq.mnam(i)) ilig(il,j)=i
      enddo
      enddo
      if(jlim.eq.6) then
      do j=4,6
      do i=nuc(is-1)+1,nuc(is)
      if(anam(il,j).eq.mnam(i)) ilig(il,j)=i
      enddo
      enddo
      else
      ilig(il,4)=nuc(is-1)+lpiv(il)
      if(j.eq.5.and.kal.eq.2) then
      ilig(il,5)=nuc(is-1)+1
      if(ilig(il,5).eq.ilig(il,4)) ilig(il,5)=ilig(il,5)+1
      endif
      endif
      do j=1,jlim
      if(ilig(il,j).eq.0) then
      write(6,98) anam(il,j),lnam(il)
98    format(/2x,'---- Linking atom ',a4,' not found for ligand ',
     1 a4,' ----'/)
      return
      endif
      enddo
      endif
      enddo
      if(nop.gt.0) call openb(nop,rpiv,tpiv)
c-------------------------------------------------------------------------output
      write(6,206) ctot,kam,kcen
206   format(/2x,'ctot= ',f8.2,' kam =',i8,' kcen= ',i8)
      jl=1
      ju=kseq
      do k=1,nst
      if(k.eq.1) write(6,*)
      nk=abs(idr(k))
      ik=sign(1,idr(k))
      write(6,201) strand(k),nk,dirn(ik),seq(jl:ju)
201   format(2x,a3,' Strand (N= ',i2,' Dir= ',a5,') : ',a)
      jl=jl+kseq
      ju=ju+kseq
      enddo
c---------------------------------------------------------------------open bases
      if(nop.gt.0) then
        write(6,*) ' '
        do i=1,nop
          is=nuc(nba(i)-1)+iofs(nba(i))
          write(6,208) munit(is),nunit(is),theta(i)
        enddo
208     format(2x,'Open: ',a4,i3,' by ',f8.3,' degrees')
      endif
c-------------------------------------------------------------noe distance input
      do i=1,nto
        in=itr(i)
        ih68(i)=0
        do j=nuc(i-1)+iofs(i),nuc(i)
          if(mnam(j).eq.'HC6'.or.mnam(j).eq.'HC8') ih68(i)=j
        enddo
      enddo
      if(noe.ne.' ') then
        call renoe(noe)
	  if(iret.eq.1) return
        write(6,209) nnoe,ndiv(nlin),fnoes,fnoew,fnoea
209     format(/2x,'NOE data: No. of constraints= ',i3,
     1 ' No. of list distances= ',i3,
     1 /2x,'(Fnoe Strong = ',f6.1,' Weak= ',f6.1,' Angle= ',f7.1,')')
        start=.true.
        do k=1,nnoe
          if(jnoe(k).eq.11) then
            if(start) then
              write(6,*)
              start=.false.
            endif
            i1=inoe(k,1)
            i2=inoe(k,2)
            write(6,301) rnoe(k),fnoe(k),nunit(i1),mnam(i1),
     1        nunit(i2),mnam(i2)
301         format(2x,'B) r0= ',f6.2,' f= ',f6.2,2x,i2,'/',a4,'- ',i2,'/',a4)
          else if(jnoe(k).eq.9) then
            if(start) then
              write(6,*)
              start=.false.
            endif
            i1=inoe(k,1)
            i2=inoe(k,2)
            i3=inoe(k,3)
            write(6,302) rnoe(k),fnoe(k),nunit(i1),mnam(i1),nunit(i2),
     1                  mnam(i2),nunit(i3),mnam(i3)
302         format(2x,'V) v0= ',f6.2,' f= ',f6.2,2x,i2,'/',a4,'- ',i2,'/',a4,
     1                   '- ',i2,'/',a4)
          else if(jnoe(k).eq.10) then
            if(start) then
               write(6,*)
               start=.false.
            endif
            i1=inoe(k,1)
            i2=inoe(k,2)
            i3=inoe(k,3)
            i4=inoe(k,4)
            write(6,303) rnoe(k),fnoe(k),nunit(i1),mnam(i1),nunit(i2),
     1         mnam(i2),nunit(i3),mnam(i3),nunit(i4),mnam(i4)
303          format(2x,'C) t0= ',f6.2,' f= ',f6.2,2x,i2,'/',a4,'- ',i2,'/',a4,
     1         '- ',i2,'/',a4,'- ',i2,'/',a4)
          endif
        enddo
      endif
c-----------------------------------------------------------------------cylinder
      if(cyl) then
        write(6,115) rcyl,eneq,tcyl,vcyl(lcyl(1)),dcyl,vcyl(lcyl(2)),
     1               hcyl,vcyl(lcyl(3))
115     format(/2x,'Cylinder active, Rad= ',f7.3,' En(eq)= ',f7.3,
     1        //2x,'                 Theta = ',f7.3,2x,a6,
     1         /2x,'                 Dist  = ',f7.3,2x,a6,
     1         /2x,'                 Height= ',f7.3,2x,a6)
        ct=cos(cdr*(tcyl))
        st=sin(cdr*(tcyl))
        pcyl(1)=dcyl*ct
        pcyl(2)=dcyl*st
        pcyl(3)=hcyl
        ucyl(1)= st
        ucyl(2)=-ct
        ucyl(3)=0.
        do i=1,25
          req=vdw(i)+rcyl
          ac(i)=-2*eneq*req**6
          bc(i)=-eneq*req**12
        enddo
      endif
c------------------------------------------------minimisation control and output
      limit=limit*60
      out_axx = out
      call closejm(ipl,iabas)
	if(iret.eq.1) return
      if(nloop.eq.0) then
        if(out_axx.ne.' '.and.axl.eq.' '.and.nol.eq.' ') call axeout(out_axx)
        if(lmo.ne.' ') call lmoout
      endif
      return
      end
