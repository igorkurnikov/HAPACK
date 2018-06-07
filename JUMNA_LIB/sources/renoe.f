      subroutine renoe(noe)
      include 'jumna_data.inc'
      logical*2 sup,rcom,homo,homo2,homo3,diep,link,ecen,cyl,lcat,
     1 lar,ihl,lock,kink,ifhb,ribose,cation,lthy,blank,tlp,quiet,convg,
     1 double,autos,dloop,cent,amber,sum,beg
      character*4 mnam,munit,seq*120,code*8,kode*8,lett*1,
     1 noe*32,inpu*80,line*80,nam(8),list*9,knam,key*1
      integer*4 opt,iun(8),iwk(80),ich(50)
      real*8 ha(n2,9),rch(50)
      common/comjm/cut,icx(n2),icn(n2),ict(n2,10),ncut,nmx,mco,ncomb,icomb
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
      common/lop/dellp(2),stlp(2),gmx,indlp(2,2),lplow(2),lphig(2),
     1 nloop,icy,tlp(2),quiet,convg,dloop
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
      common/sst/phase(8,3),ampli(8,3),ph5(n2),am5(n2),kr5(n2)
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      common/symjm/xsum,ysum,itwl(n2),nvs(n7),ihm(n2),ksym(3),nbrk(2),
     1 isym,isur,isup,ihl(n2)
	common/hacon/iret
      equivalence (ha(1,1),ua(1,1))
      ml=0
      mco=0
      icmin=10000
      do i=1,n7
        nvs(i)=0
      enddo
      nch=0
c------------------------------------------------------------------read noe data
      nlis=0
      nlin=0
      ndiv(0)=0
      kfi=index(noe,' ')-1
      open(unit=2,file=noe(:kfi)//'.noe',status='old')
      kc=0
10    read(2,5,end=100) inpu
5     format(a)
      key=inpu(:1)
      if(inpu.eq.' '.or.key.eq.'*'.or.key.eq.'#') goto 10
      if(key.ne.'E'.and.key.ne.'F'.and.
     1   key.ne.'L'.and.key.ne.'R'.and.key.ne.'M') kc=kc+1  ! not adiabatic mapping,list dist or change charges
      if(kc.ge.n3) then
        write(6,*) '  ---- n3 too small ----'
        iret = 1
	  return
      endif
      double=.false.
      line=inpu(2:)
      if(key.eq.'A'.or.key.eq.'a') goto 70  ! sugar amplitude
      if(key.eq.'B') goto 25
      if(key.eq.'C') goto 90
      if(key.eq.'D'.or.key.eq.'d') goto 20  ! interatomic distance
      if(key.eq.'E') goto 130
      if(key.eq.'F') goto 120
      if(key.eq.'G') goto 110
      if(key.eq.'H'.or.key.eq.'h') goto 125
      if(key.eq.'K') goto 115
      if(key.eq.'L') goto 50
      if(key.eq.'M') goto 170
      if(key.eq.'O') go to 85
      if(key.eq.'P'.or.key.eq.'p') goto 60  ! sugar phase
      if(key.eq.'R') goto 180
      if(key.eq.'S'.or.key.eq.'s') goto 160
      if(key.eq.'T'.or.key.eq.'t') goto 30
      if(key.eq.'U') goto 140
      if(key.eq.'V') goto 80
      if(key.eq.'W') goto 45
      if(key.eq.'X'.or.key.eq.'x') goto 55
      if(key.eq.'Z') goto 40
      write(6,*) '  ---- Unknown NOE file key letter ----'
      iret = 1
	return
c----------------------------------------------------------------------distances
20    if (line(:1).eq.'2') then
      double=.true.
      line=inpu(3:)
      endif
      i=0
      m=1
81    i=i+1
      if(line(i:i).eq.' ') goto 81
      ibeg=i
82    i=i+1
      if(line(i:i).ne.' ') goto 82
      nam(m)=line(ibeg:i-1)
      m=m+1
      if(m.le.2) goto 81
      lup=2
      if(key.eq.'D') then
      read(line(i:),*) iun(1),iun(2),rnoe(kc)
      jnoe(kc)=1
      else
      read(line(i:79),*) iun(1),iun(2),rnoe(kc),bnoe(kc)
      jnoe(kc)=3
      if(bnoe(kc).lt.rnoe(kc)) then
      write(6,*) '  ---- BRACKET DISTANCES INCORRECT ----'
      iret = 1
	return
      endif
      endif
         if(double) then
         read(2,5) inpu
         inoe(kc,5)=1
         if(inpu(:1).eq.'-') inoe(kc,5)=-1
         line=inpu(2:)
         i=0
103      i=i+1
         if(line(i:i).eq.' ') goto 103
         ibeg=i
104      i=i+1
         if(line(i:i).ne.' ') goto 104
         nam(m)=line(ibeg:i-1)
         m=m+1
         if(m.le.4) goto 103
         read(line(i:),*) iun(3),iun(4)
         jnoe(kc)=jnoe(kc)+25
         lup=4
         endif
      do l=1,lup
      is=iun(l)
      in=itr(is)
      ioff=nuc(is-1)+iofs(is)
      do i=nuc(is-1)+1,nuc(is)
      if(mnam(i).eq.nam(l)) then
      inoe(kc,l)=i
      goto 27
      endif
      enddo
      if(nam(l).eq.'H68') then
      inoe(kc,l)=ih68(is)
      goto 27
      endif
      write(6,26) kc,nam(l)
26    format(/2x,'---- Noe input ',i3,' atom ',a4,' unknown ----'/)
      iret = 1
	return
27    enddo
      if(rnoe(kc).le.3.0.or.key.eq.'d'.or.double) then
      fnoe(kc)=fnoes
      else
      fnoe(kc)=fnoew
      endif
      goto 10
c-----------------------------------------------------------------chemical bonds
25    i=0
      m=1
801   i=i+1
      if(line(i:i).eq.' ') goto 801
      ibeg=i
802   i=i+1
      if(line(i:i).ne.' ') goto 802
      nam(m)=line(ibeg:i-1)
      m=m+1
      if(m.le.2) goto 801
      read(line(i:),*) iun(1),iun(2),rnoe(kc),fnoe(kc)
      jnoe(kc)=11
      do l=1,2
      is=iun(l)
      in=itr(is)
      ioff=nuc(is-1)+iofs(is)
      do i=nuc(is-1)+1,nuc(is)
      if(mnam(i).eq.nam(l)) then
      inoe(kc,l)=i
      goto 207
      endif
      enddo
      write(6,26) kc,nam(l)
      iret = 1
	return
207   enddo
c--------------------------update bonding matrix
          i1=inoe(kc,1)
          i2=inoe(kc,2)
          m1=matd(i1,7)+1
          m2=matd(i2,7)+1
          matd(i1,m1)=i2
          matd(i1,7)=m1
          matd(i2,m2)=i1
          matd(i2,7)=m2
c-----------------------------------------------
      goto 10
c-----------------------------------------------------------------------torsions
30    if(line(:1).eq.'2') then
      double=.true.
      line=inpu(3:)
      endif
      i=0
      m=1
91    i=i+1
      if(line(i:i).eq.' ') goto 91
      ibeg=i
92    i=i+1
      if(line(i:i).ne.' ') goto 92
      nam(m)=line(ibeg:i-1)
      m=m+1
      if(m.le.4) goto 91
      if(key.eq.'T') then
      read(line(i:),*) iun(1),iun(2),iun(3),iun(4),rnoe(kc)
      jnoe(kc)=2
      else
      read(line(i:),*) iun(1),iun(2),iun(3),iun(4),rnoe(kc),
     1 bnoe(kc)
      jnoe(kc)=4
      if(abs(bnoe(kc)).gt.180.) bnoe(kc)=bnoe(kc)-sign(360.d0,bnoe(kc))
      if(bnoe(kc).eq.90..or.bnoe(kc).eq.-90.) bnoe(kc)=bnoe(kc)+1
      endif
      if(abs(rnoe(kc)).gt.180.) rnoe(kc)=rnoe(kc)-sign(360.d0,rnoe(kc))
      if(rnoe(kc).eq.90..or.rnoe(kc).eq.-90.) rnoe(kc)=rnoe(kc)+1
      lup=4
         if(double) then
         read(2,5) inpu
         inoe(kc,9)=1
         if(inpu(:1).eq.'-') inoe(kc,9)=-1
         line=inpu(2:)
         i=0
101      i=i+1
         if(line(i:i).eq.' ') goto 101
         ibeg=i
102      i=i+1
         if(line(i:i).ne.' ') goto 102
         nam(m)=line(ibeg:i-1)
         m=m+1
         if(m.le.8) goto 101
         read(line(i:),*) iun(5),iun(6),iun(7),iun(8)
         jnoe(kc)=jnoe(kc)+12
         lup=8
         endif
      do l=1,lup
      is=iun(l)
      do i=nuc(is-1)+1,nuc(is)
      if(mnam(i).eq.nam(l)) then
      inoe(kc,l)=i
      goto 35
      endif
      enddo
      if(nam(l).eq.'H68') then
      inoe(kc,l)=ih68(is)
      goto 35
      endif
      write(6,26) kc,nam(l)
      iret = 1
	return
35    enddo
      fnoe(kc)=fnoea
      goto 10
c---------------------------------------------------------------------total rise
 40   read(line,*) iu1,iu2,rnoe(kc)
         if(iu1.eq.0.and.iu2.eq.0) then
           inoe(kc,1)=2
           inoe(kc,2)=kseq
         else if(iu1.ge.2.and.iu2.le.kseq.and.iu1.le.iu2) then
           inoe(kc,1)=iu1
           inoe(kc,2)=iu2
         else
           write(6,*) '  ---- Twist/Rise: 1 < iu1 =< iu2 =< kseq ----'
           iret = 1
	     return
         endif
         jnoe(kc)=12
         fnoe(kc)=fnoea
         nzsh=kc
         goto 10
c--------------------------------------------------------------------total twist
 45   read(line,*) iu1,iu2,rnoe(kc)
         if(iu1.eq.0.and.iu2.eq.0) then
         inoe(kc,1)=2
         inoe(kc,2)=kseq
         else if(iu1.ge.2.and.iu2.le.kseq.and.iu1.le.iu2) then
         inoe(kc,1)=iu1
         inoe(kc,2)=iu2
         else
         write(6,*) '  ---- Twist/Rise: 1 < iu1 =< iu2 =< kseq ----'
         iret = 1
	   return
         endif
      jnoe(kc)=13
      fnoe(kc)=fnoea
      nwdg=kc
      goto 10
c------------------------------------------------------------------------X twist
 55   read(line,*) iu1,iu2,rnoe(kc)
      jnoe(kc)=30
      jdel=iu2-iu1
      if(key.eq.'x') jdel=1
      fnoe(kc)=fnoea
      ic=0
      do i=iu1,iu2,jdel
        icon(ic+1)=nuc(i-1)+1
        level=ilq(i,1)
        ipair=ieq(level,2)
        if(ilq(i,2).ne.1) ipair=ieq(level,1)
        if(ipair.gt.0) then
          icon(ic+2)=nuc(ipair-1)+1
        else
          write(6,28)
          iret = 1
	    return
        endif
        ic=ic+2
      enddo
      inoe(kc,1)=ic/2-1
      goto 10
c--------------------------------------------------------------------sugar phase
60    if(line(:1).eq.'2') then
        double=.true.
        line=inpu(3:)
      endif
      if(key.eq.'P') then
        if(.not.double) then
          read(line,*) iu1,rnoe(kc)
          jnoe(kc)=5
        else
          read(line,*) iu1,iu2,rnoe(kc)
          inoe(kc,9)=sign(1,iu2)
          iu2=abs(iu2)
          jnoe(kc)=18
        endif
      else
        if(.not.double) then
          read(line,*) iu1,rnoe(kc),bnoe(kc)
          jnoe(kc)=6
        else
          read(line,*) iu1,iu2,rnoe(kc),bnoe(kc)
          inoe(kc,9)=sign(1,iu2)
          iu2=abs(iu2)
          jnoe(kc)=19
        endif
        if(bnoe(kc).lt.rnoe(kc)) then
          write(6,*) '  ---- BRACKET DISTANCES INCORRECT ----'
          iret = 1
	    return
        endif
      endif
      if(iu1.gt.nto) then
        write(6,22) iu1
 22     format(/2x,'---- ',i2,' has no sugar to constrain ----'/)
        iret = 1
	  return
      endif
      if(double) then
        if(iu2.gt.nto) then
            write(6,22) iu2
            iret = 1
	    return
        endif
         inoe(kc,2)=iu2
      endif
      inoe(kc,1)=iu1
      fnoe(kc)=fnoea
      goto 10
c----------------------------------------------------------------sugar amplitude
70    if(line(:1).eq.'2') then
        double=.true.
        line=inpu(3:)
      endif
      if(key.eq.'A') then
        if(.not.double) then
          read(line,*) iu1,rnoe(kc)
          jnoe(kc)=7
        else
          read(line,*) iu1,iu2,rnoe(kc)
          inoe(kc,9)=sign(1,iu2)
          iu2=abs(iu2)
          jnoe(kc)=24
        endif
      else
        if(.not.double) then
          read(line,*) iu1,rnoe(kc),bnoe(kc)
          jnoe(kc)=8
        else
          read(line,*) iu1,iu2,rnoe(kc),bnoe(kc)
          inoe(kc,9)=sign(1,iu2)
          iu2=abs(iu2)
          jnoe(kc)=25
        endif
        if(bnoe(kc).lt.rnoe(kc)) then
          write(6,*) '  ---- BRACKET DISTANCES INCORRECT ----'
          iret = 1
	    return
        endif
      endif
      if(iu1.gt.nto) then
        write(6,22) iu1
        iret = 1
	  return
      endif
      if(double) then
        if(iu2.gt.nto) then
          write(6,22) iu2
          iret = 1
	    return
        endif
        inoe(kc,2)=iu2
      endif
      inoe(kc,1)=iu1
      fnoe(kc)=fnoea
      goto 10
c-------------------------------------------------------------------sugar pucker
160   read(line,*) iu1,lett
      jnoe(kc)=21
      if(key.eq.'s') jnoe(kc)=22
      if(iu1.gt.nto) then
        write(6,22) iu1
        iret = 1
	  return
      endif
      inoe(kc,1)=iu1
      if(lett.eq.'S') then
        inoe(kc,2)=1
      else if(lett.eq.'Y') then
        inoe(kc,2)=2
      else if(lett.eq.'X') then
        inoe(kc,2)=3
      else if(lett.eq.'E') then
        inoe(kc,2)=4
      else if(lett.eq.'N') then
        inoe(kc,2)=5
      else if(lett.eq.'M') then
        inoe(kc,2)=6
      else if(lett.eq.'A') then
        inoe(kc,2)=7
      else if(lett.eq.'H') then
        inoe(kc,2)=8
      else
        write(6,23) lett
23      format(/2x,'---- ',a1,' is an unknown pucker code ----'/)
        iret = 1
	  return
      endif
      fnoe(kc)=fnoea
      goto 10
c-----------------------------------------------------------------valence angles
80    i=0
      m=1
901   i=i+1
      if(line(i:i).eq.' ') goto 901
      ibeg=i
902   i=i+1
      if(line(i:i).ne.' ') goto 902
      nam(m)=line(ibeg:i-1)
      m=m+1
      if(m.le.3) goto 901
      read(line(i:),*) iun(1),iun(2),iun(3),rnoe(kc),fnoe(kc)
      jnoe(kc)=9
      if(rnoe(kc).gt.180.or.rnoe(kc).lt.0) then
      write(6,903) kc
903   format(/2x,'---- V contraint ',i3,' must be >0 and <180 ----')
      istep = 1
	return
      endif
      do l=1,3
      is=iun(l)
      do i=nuc(is-1)+1,nuc(is)
      if(mnam(i).eq.nam(l)) then
      inoe(kc,l)=i
      goto 305
      endif
      enddo
      write(6,26) kc,nam(l)
      iret = 1
	return
305   enddo
      goto 10
c------------------------------------------------------------------------opening
85    read(line,*) iw,rnoe(kc)
      inoe(kc,1)=nuc(iw-1)+1
      inoe(kc,2)=nuc(iw-1)+iofs(iw)
      lev=ilq(iw,1)
      if(ilq(iw,2).eq.1) then
      iwb=ieq(lev-1,1)
      iwa=ieq(lev+1,1)
      ic =ieq(lev,2)
      icb=ieq(lev-1,2)
      ica=ieq(lev+1,2)
         else
         iwb=ieq(lev-1,2)
         iwa=ieq(lev+1,2)
         ic =ieq(lev,1)
         icb=ieq(lev-1,1)
         ica=ieq(lev+1,1)
         endif
      inoe(kc,8)=idr(ilq(iw,2))
      if(iwb.eq.0.or.iwa.eq.0.or.ic.eq.0.or.icb.eq.0.or.ica.eq.0) then
      write(6,28)
28    format(/2x,'ERROR: Open/Xtwi - missing paired/neighbour nucleot')
      iret =1
	return
      endif
      inoe(kc,3)=nuc(ic-1)+1
      inoe(kc,4)=nuc(iwb-1)+1
      inoe(kc,5)=nuc(iwa-1)+1
      inoe(kc,6)=nuc(icb-1)+1
      inoe(kc,7)=nuc(ica-1)+1
      jnoe(kc)=29
      fnoe(kc)=fnoea
      goto 10
c----------------------------------------------------------------cosine torsions
90    i=0
      m=1
904   i=i+1
      if(line(i:i).eq.' ') goto 904
      ibeg=i
905   i=i+1
      if(line(i:i).ne.' ') goto 905
      nam(m)=line(ibeg:i-1)
      m=m+1
      if(m.le.4) goto 904
      read(line(i:),*) iun(1),iun(2),iun(3),iun(4),rnoe(kc),fnoe(kc)
      jnoe(kc)=10
      do l=1,4
      is=iun(l)
      do i=nuc(is-1)+1,nuc(is)
      if(mnam(i).eq.nam(l)) then
      inoe(kc,l)=i
      goto 306
      endif
      enddo
      write(6,26) kc,nam(l)
      iret = 1
	return
306   enddo
      goto 10
c-------------------------------------------------------------------groove width
125   l=0
      beg=.false.
      do i=80,1,-1
      if(.not.beg.and.line(i:i).ne.' ') then
      beg=.true.
      l=l+1
         else if(beg.and.line(i:i).eq.' ') then
         beg=.false.
         endif
      enddo
      if(key.eq.'H') then
      if(l.gt.9) then
      write(6,*) '  ---- Too many H indices ----'
      iret = 1
	return
      endif
      read(line,*) (inoe(kc,m),m=1,l-1),rnoe(kc)
      inoe(kc,9)=l-1
      jnoe(kc)=31
         else
         if(l.gt.10) then
         write(6,*) '  ---- Too many H indices ----'
         iret = 1
	   return
         endif
         read(line,*) (inoe(kc,m),m=1,l-2),rnoe(kc),bnoe(kc)
         inoe(kc,9)=l-2
         jnoe(kc)=32
         if(bnoe(kc).lt.rnoe(kc)) then
         write(6,*) '  ---- GROOVE WIDTH BRACKET INCORRECT ----'
         iret = 1
	   return
         endif
      do m=1,inoe(kc,9)
      if(inoe(kc,m).gt.iend(1)) then
      write(6,*) '  ---- H indicies not in 1st strand ----'
      iret = 1
	return
      endif
      enddo
      endif
      fnoe(kc)=fnoes
      goto 10
c----------------------------------------------------------------ligand symmetry
110   blank=.true.
      ic=0
      do i=1,79
      if(blank.and.line(i:i).ne.' ') then
      blank=.false.
      ic=ic+1
      endif
      if(.not.blank.and.line(i:i).eq.' ') blank=.true.
      enddo
      read(line,*) (iwk(i),i=1,ic)
      ii=iwk(1)
      do i=2,ic
      nvs(iwk(i))=-ii
      enddo
      goto 10
c------------------------------------------------------------------variable sums
115   if(line(:1).eq.'2') then
      read(line(2:),*) (inoe(kc,i),i=1,2),rnoe(kc)
      jnoe(kc)=15
      fnoe(kc)=fnoea
      else if(line(:1).eq.'4') then
      read(line(2:),*) (inoe(kc,i),i=1,4),rnoe(kc)
      jnoe(kc)=17
      fnoe(kc)=fnoea
      else
      write(6,*) '  ---- Error in .noe S option use 2 or 4 ----'
      iret = 1
	return
      endif
      goto 10
c------------------------------------------------------------read list distances
50    i=0
51    i=i+1
      if(i.eq.78) then
      nlin=nlin+1
      ndiv(nlin)=nlis
      if(nlin.gt.99) then
      write(6,*) '  ---- Too many lines of NOE list data ----'
      iret = 1
	return
      endif
      if(ndiv(nlin)-ndiv(nlin-1).gt.7) then
      write(6,53) nlin
53    format(2x,'---- Too many NOE list data in line ',i2,' ----')
      iret = 1
	return
      endif
      goto 10
      endif
      if(line(i:i).eq.' ') goto 51
      ibeg=i
52    i=i+1
      if(line(i:i).ne.' ') goto 52
      nlis=nlis+1
      if(nlis.gt.200) then
      write(6,*) '  ---- Too many NOE list distances ----'
      iret = 1
	return
      endif
      list=line(ibeg:i-1)
      if(index(list,'-').ne.0) then
      k=index(list,'-')
      ktyp(nlis)=1
      else if(index(list,'^').ne.0) then
      k=index(list,'^')
      ktyp(nlis)=2
      else if(index(list,'v').ne.0) then
      k=index(list,'v')
      ktyp(nlis)=3
      else if(index(list,':').ne.0) then
      k=index(list,':')
      ktyp(nlis)=4
      else if(index(list,'/').ne.0) then
      k=index(list,'/')
      ktyp(nlis)=5
      else if(index(list,'\\').ne.0) then
      k=index(list,'\\')
      ktyp(nlis)=6
      else
      write(6,*) '  ---- Unrecognized NOE list symbol ----'
      iret = 1
	return
      endif
      knam(nlis,1)=list(:k-1)
      knam(nlis,2)=list(k+1:)
      goto 51
c------------------------------------------------------------------------F loops
120   if(nloop.lt.0) then
      write(6,*) '  ---- Cannot combine combinatorial and loops ----'
      iret = 1
	return
      endif
      nloop=nloop+1
      if(nloop.gt.2) then
      write(6,*) '  ---- Number of .noe For loops > 2 ----'
      iret = 1
	return
      endif
      i=index(line,':')
      if(i.gt.0) then
      read(line(:i-1),*) indlp(nloop,1)
      read(line(i+1:),*) indlp(nloop,2),dellp(nloop),lplow(nloop),
     1 lphig(nloop)
         else
         read(line,*) indlp(nloop,1),dellp(nloop),lplow(nloop),
     1   lphig(nloop)
         indlp(nloop,2)=indlp(nloop,1)
         endif
      if(lplow(nloop).gt.0.or.lphig(nloop).lt.0.or.
     1  -lplow(nloop)+lphig(nloop).lt.0) then
      write(6,*) '  ---- Error in lplow/lphig ----'
      iret = 1
	return
      endif
      tlp(nloop)=.true.
      goto 10
c------------------------------------------------------------------------E loops
130   if(nloop.lt.0) then
      write(6,*) '  ---- Cannot combine combinatorial and loops ----'
      iret = 1
	return
      endif
      nloop=nloop+1
      if(nloop.gt.2) then
      write(6,*) '  ---- Number of .noe For loops > 2 ----'
      iret = 1
	return
      endif
      i=index(line,':')
      if(i.gt.0) then
      read(line(:i-1),*) indlp(nloop,1)
      read(line(i+1:),*) indlp(nloop,2),dellp(nloop),lplow(nloop),
     1 lphig(nloop),stlp(nloop)
         else
         read(line,*) indlp(nloop,1),dellp(nloop),lplow(nloop),
     1   lphig(nloop),stlp(nloop)
         indlp(nloop,2)=indlp(nloop,1)
         endif
      if(lplow(nloop).gt.0.or.lphig(nloop).lt.0.or.
     1  -lplow(nloop)+lphig(nloop).lt.0) then
      write(6,*) '  ---- Error in lplow/lphig ----'
      iret = 1
	return
      endif
      tlp(nloop)=.false.
      goto 10
c-----------------------------------------------------------Combinatorial search
170   if(nloop.gt.0) then
      write(6,*) '  ---- Cannot combine combinatorial and loops ----'
      iret = 1
	return
      endif
      nloop=nloop-1
      il=index(line,'/')
         if(il.eq.0) then
         write(6,*) '  ---- Combinatorial symbol / missing ----'
         iret = 1
	   return
         endif
      i=0
921   i=i+1
      if(line(i:i).eq.' ') goto 921
      ibeg=i
922   i=i+1
      if(line(i:i).ne.' ') goto 922
      if(i.gt.il) goto 923
      mco=mco+1
      read(line(ibeg:i-1),*) icx(mco)
      goto 921
923   m=0
931   i=i+1
      if(i.gt.80) goto 933
      if(line(i:i).eq.' ') goto 931
      ibeg=i
932   i=i+1
      if(line(i:i).ne.' ') goto 932
      m=m+1
      read(line(ibeg:i-1),*) ict(mco,m)
      if(ict(mco,m).lt.icmin) icmin=ict(mco,m)
      goto 931
933   icn(mco)=m
      do i=ml+1,mco-1
      icn(i)=icn(mco)
      do j=1,m
      ict(i,j)=ict(mco,j)
      enddo
      enddo
      ml=mco
      goto 10
c----------------------------------------------------------------------Curvature
140   read(line,*) inoe(kc,1),inoe(kc,2),rnoe(kc),bnoe(kc)
      if(inoe(kc,2).gt.kseq) then
      write(6,*) '  ---- Units for curvature must be in 1st strand ----'
      iret = 1
	return
      endif
      if(inoe(kc,1).ge.inoe(kc,2)) then
      write(6,*) '  ---- Units for curvature must be i1<i2 ----'
      iret = 1
	return
      endif
      jnoe(kc)=20
      fnoe(kc)=fnoea
         kc=kc+1
         jnoe(kc)=0
         rnoe(kc)=bnoe(kc-1)
      goto 10
c-----------------------------------------------------------------modify charges
180   i=0
181   i=i+1
      if(line(i:i).eq.' ') goto 181
      ibeg=i
182   i=i+1
      if(line(i:i).ne.' ') goto 182
      nam(1)=line(ibeg:i-1)
      read(line(i:),*) is,charge
      in=itr(is)
      do i=nuc(is-1)+1,nuc(is)
      if(mnam(i).eq.nam(1)) goto 183
      enddo
      write(6,26) kc,nam(1)
      iret = 1
	return
183   dmon(i)=charge
      nch=nch+1
      ich(nch)=i
      rch(nch)=charge
      goto 10
c-------------------------------------------------------------------end of input
100   close(2)
      nnoe=kc
      if(mco.gt.0) then
      icomb=nnoe
      ncomb=1
      do i=1,mco
      ncomb=ncomb*icn(i)
      enddo
      write(6,114) ncomb,mco
114   format(/2x,'Combinatorial choices ...',i9,' for ',i2,' subunits'/)
      do i=1,mco
      write(6,112) icx(i),icn(i),(ict(i,j),j=1,icn(i))
112   format(2x,i3,'  (',i1,')  ',10i3)
      enddo
      nnoe=icmin-1
      endif
         if(nch.gt.0) then
         write(6,185) nch
185      format(/2x,'Modified ',i3,' charges ...'/)
         write(6,186) (mnam(ich(i)),nunit(ich(i)),rch(i),i=1,nch)
186      format(4(:,2x,a4,i3,f7.3,' /'))
         endif
      return
      end
