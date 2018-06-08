      subroutine closejm(ipl,iabas)
      include 'jumna_data.inc'
      logical*2 lar,lock,lthy,kink,ifhb,hst,bst,vst,lgi,lgj,lgu,
     1 ribose,cation,sup,rcom,homo,homo2,homo3,diep,link,ecen,sum,
     1 cyl,ihl,lcat,cent,autos,locr,tlp,quiet,convg,hom(3),ok,brk(n2),
     1 dloop,amber
      character*4 mnam,munit,seq*120,code*8,kode*8,type*8,
     1 mac*32,lmo*32,lib*32,libm*32,out*32,axe*32,noe*32,nol*32,axl*32,
     1 test*32,pdb*32,ins*32,bar*32,parm*32,lnam,knam
      integer*2 i23,i34,elim
      integer*4 opt
      dimension cvar(3),ipl(n9)
      common/cha/mac,lmo,lib,libm,out,axe,axl,noe,nol,test,pdb,ins,
     1 bar,parm
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
      common/ind/iend(0:4),ise(n2),ienb(n2,2),kapt(n2+1),nsph
      common/lgd/rlig(n9,6),slig(n9,2),ilig(n9,6),lopt(n9),lpiv(n9),
     1 ntlg,nlgi,nlig,ntl,lnam(n9),locr(n2,n0)
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
      common/slg/hst(n2,6),bst(n2,n8),vst(n2,4),lgi(n1),lgj(n1),lgu(n2)
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      common/symjm/xsum,ysum,itwl(n2),nvs(n7),ihm(n2),ksym(3),nbrk(2),
     1 isym,isur,isup,ihl(n2)
	common/hacon/iret
      equivalence (hom(1),homo)
c------------------------------------------------------------------setup locking
      k=0
      ki=0
      do is=1,nto
        in=itr(is)
        ino=ito(is)
        type=code(is)
        locr(is,1)=.false.
        if(index(type,'S').ne.0.or.index(type,'A').ne.0) locr(is,1)=.true. ! lock all backbone or sugar pars
        do l=1,kap(2,in)
          k=k+1
          if(index(type,'A').ne.0.and.l.le.kap(1,ino)) lock(k)=.true.
          if(index(type,'G').ne.0.and.l.eq.4) lock(k)=.true.
          if(index(type,'S').ne.0.and.nap(l,6,in).gt.0) lock(k)=.true.
          if(index(type,'R').ne.0.and.ribose(is).and.
     1       l.eq.kap(1,ino)) lock(k)=.true.
          if(index(type,'B').ne.0.and.ino.ne.3.and.ino.ne.6) then
             if(l.eq.7.or.l.eq.8) lock(k)=.true.
             if(nap(l,6,in).eq.-8) lock(k)=.true.
          endif
          if(index(type,'N').ne.0.and.ino.ne.3.and.ino.ne.6) then
            if(l.eq.9.or.l.eq.10) lock(k)=.true.
          endif
        enddo
      enddo
      nrin=0
      nsug=0
      do is=1,ntl
        if(is.le.nto.and..not.locr(is,1)) nsug=nsug+1
        do lr=1,nsr(is)
          if(.not.locr(is,lr)) nrin=nrin+1
        enddo
        in=itr(is)
        if(is.gt.nto) k=k+kap(2,in)
      enddo
c----------------------------------------------------------------helical locking
      do is=1,nto
        type=code(is)
        do l=1,6
          k=k+1
c         lock(k)=.false.
          if(is.eq.1.and.(l.eq.3.or.(isur.eq.0.and.l.eq.6))) lock(k)=.true.
          if(index(type,'H').ne.0) lock(k)=.true.
          if(index(type,'X').ne.0.and.l.eq.1) lock(k)=.true.
          if(index(type,'Y').ne.0.and.l.eq.2) lock(k)=.true.
          if(index(type,'Z').ne.0.and.l.eq.3) lock(k)=.true.
          if(index(type,'I').ne.0.and.l.eq.4) lock(k)=.true.
          if(index(type,'T').ne.0.and.l.eq.5) lock(k)=.true.
          if(index(type,'W').ne.0.and.l.eq.6) lock(k)=.true.
        enddo
      enddo
c--------------------------------------------------------------------------kinks
      do is=1,kseq
        if(kink(is)) then
          type=kode(is)
          do l=1,4
            k=k+1
            lock(k)=.false.
            if(index(type,'X').ne.0.and.l.eq.1) lock(k)=.true.
            if(index(type,'Y').ne.0.and.l.eq.2) lock(k)=.true.
            if(index(type,'I').ne.0.and.l.eq.3) lock(k)=.true.
            if(index(type,'T').ne.0.and.l.eq.4) lock(k)=.true.
          enddo
        endif
      enddo
c-------------------------------------------------------------------ligand inter
      do is=nto+1,ntl
        kal=nuc(is)-nuc(is-1)
        type=code(is)
        do l=1,3
          k=k+1
          lock(k)=.false.
          if(index(type,'H').ne.0) lock(k)=.true.
          if(index(type,'X').ne.0.and.l.eq.1) lock(k)=.true.
          if(index(type,'Y').ne.0.and.l.eq.2) lock(k)=.true.
          if(index(type,'Z').ne.0.and.l.eq.3) lock(k)=.true.
        enddo
        if(kal.gt.1) then
          do l=4,6
            k=k+1
            if(kal.eq.2.and.l.eq.4)  lock(k)=.true.
            if(index(type,'H').ne.0) lock(k)=.true.
            if(index(type,'I').ne.0.and.l.eq.4) lock(k)=.true.
            if(index(type,'T').ne.0.and.l.eq.5) lock(k)=.true.
            if(index(type,'W').ne.0.and.l.eq.6) lock(k)=.true.
          enddo
        endif
      enddo
c------------------------------------------------thymine methyl and ribose count
      numt=0
      numk=0
      do is=1,nto
        if(lthy(is)) then
          numt=numt+1
          type=code(is)
          k=k+1
          lock(k)=.true.
          if(index(type,'M').eq.0) lock(k)=.false.
        endif
        if(ribose(is)) numk=numk+1
      enddo
c--------------------------------------------------------------setup ihl and ihm
      do i=1,kseq
        ihl(i)=.false.
        brk(i)=.false.
        if(i.eq.nbrk(1).or.i.eq.nbrk(2)) brk(i)=.true.
      enddo
      do k=1,nst
        if(k.eq.1) goto 100
        jl=ilq(iend(k-1)+1,1)
        ok=.false.
        islice=1
        do j=1,kseq
          i=ieq(j,k)
          if(i.ne.0) ihl(i)=.false.
          if(i.eq.0.and.brk(j)) ok=.true.
          if(brk(j)) islice=islice+1
          if(i.ne.0) then
            if(.not.hom(islice).or.j.eq.jl.or.brk(j).or.ok) ihl(i)=.true.
            ok=.false.
          endif
        enddo
100     do i=iend(k-1)+1,iend(k)
          ihm(i)=iend(k)
          if(k.gt.1) then
            do j=i+1,iend(k)
              if(ihl(j)) then
                 ihm(i)=j-1
                 goto 200
              endif
            enddo
          endif
200     enddo
      enddo
c----------------------------------------------------2nd strand zsh,wdg to local
      do kp=2,nst
        zloc=0.
        wloc=0.
        do is=1,kseq
          js=ieq(is,kp)
          if(js.gt.0) then
             zloc=zloc+hel(js,3)-hel(is,3)
             wloc=wloc+hel(js,6)-hel(is,6)
             if(ihl(js)) then
                hel(js,3)=zloc
                if(isur.eq.0) hel(js,6)=wloc
             endif
          else
             zloc=zloc-hel(is,3)
             wloc=wloc-hel(is,6)
          endif
        enddo
      enddo
c------------------------------------------------set variables and check closure
      ncon=nrin+3*(nto-nst)
      ndcs=ncon
      ncon=ncon+nnoe
      call setgeo(axe,iabas)
      call ligput(ipl)
      call setvar
      call helix
      call backbo
      call equivjm(nspc,nbcc,nhlc,nklc,nlgc)
	if(iret.eq.1) return
      if(ecen) then
        call pairc(mp)
      else
        call pairs(mp)
      endif
c---------------------------------------------------------------check dimensions
      if(kam.gt.n1) write(6,112) kam-n1
112   format(/2x,'---- OVERFLOW OF MAX ATOMS, N1 BY ',i5,' ----'/)
      if(ntba.gt.n6) write(6,113) ntba-n6
113   format(/2x,'---- OVERFLOW OF MAX BACK PHYVAR N6 BY ',i5,' ----'/)
      if(nvar.gt.n7)  write(6,114) nvar-n7
114   format(/2x,'---- OVERFLOW OF MAX INDEP VAR N7 BY ',i5,' ----'/)
      if(ncon.gt.(4*n2+n3))  write(6,115) ncon-4*n2-n3
115   format(/2x,'---- OVERFLOW OF MAX CON (4N2+N3) BY ',i5,' ----'/)
      if(ntot.gt.n6a) write(6,116) ntot-n6a
116   format(/2x,'---- OVERFLOW OF MAX PHYVAR N6A BY ',i5,' ----'/)
      if(kam.gt.n1.or.ntba.gt.n6.or.nvar.gt.n7.or.
     1 ntot.gt.n6a.or.ncon.gt.(4*n2+n3))then
		iret = 1
	    return
	endif
c------------------------------------------------------------------cylinder mods
      if(cyl) then
        cvar(1)=tcyl
        cvar(2)=dcyl
        cvar(3)=hcyl
        do i=1,3
          if(lcyl(i).eq.1) then
            nvar=nvar+1
            nvrc=nvrc+1
            neq(nvar)=nvrc
            lcyl(i)=nvar
            var(nvar)=cvar(i)
          endif
        enddo
      endif
c-------------------------------------------------------------------minimization
      if(opt.gt.0) write(6,5) vfac,tfac,rfac,xfac,scale,maxn,acc
5        format(/2x,'Minimization: ',
     1              'Vfac= ',f6.2,' Tfac= ',f6.2,' Rfac= ',f6.2,' Xfac= ',f6.2,
     1   /16x,'Scal= ',f6.2,' Maxn= ',i6  ,' Acc = ',e9.2)
         write(6,10) numt,numk,nrin,ntba,nbac,nbcc,nthe-ntba,
     1               nhel-nbac,nhlc-nbcc,ntki-nthe,nkin-nhel,ncon,ntlg-ntki,
     1               nlgi-nkin,nlgc-nklc,ntot,nvar,nvrc
10    format(/2x,'Variables   : ',
     1      'Nthy= ',i6,' Nrib= ',i6,' Nrin= ',i6,
     1 /16x,'Ntba= ',i6,' Nbac= ',i6,' Nbcc= ',i6,
     1 /16x,'Nthe= ',i6,' Nhel= ',i6,' Nhlc= ',i6,
     1 /16x,'Nkti= ',i6,' Nkin= ',i6,' Ncon= ',i6,
     1 /16x,'Ntlg= ',i6,' Nlgi= ',i6,' Nlgc= ',i6,
     1 /16x,'Ntot= ',i6,' Nvar= ',i6,' Nvrc= ',i6)
      if(opt.le.0) call microb
      if(opt.gt.0) then
         if(nloop.gt.0) then
             call loops
             return
         else if(nloop.lt.0) then
             call combi
             return
         else if(nol.ne.' ') then
             call tjunl
             return
         else if(axl.ne.' ') then
             call tbasl
             return
         else
             call minim
         endif
      endif
c-------------------------------------------------------------------------------
      if(opt.ne.0) then
         if(axl.ne.' ') then
            call tbasl
            return
         endif
      call helix
      call backbo
      call penalty
      if(parm.eq.'Flex')then
        call ecomp
      else if(parm.eq.'Amber91') then
        call ecomp91
      else if(parm.eq.'Amber94') then
        call ecomp94
      endif
      write(6,20)
20    format(/2x,'Corrected final energy ...')
      write(6,21) repl,disp,repl+disp,epen,elec,eang,etog,ener-epen
21    format(
     1 /2x,'REPL= ',f11.3,' DISP= ',f11.3,' LJ  = ',f11.3,
     1 ' PEN = ',f11.3,/2x,'ELEC= ',f11.3,' ANGL= ',f11.3,
     1 ' TORG= ',f11.3,' TOT = ',f11.3/)
      endif
c--------------------------------------------------------------noe and gradtests
      if(noe.ne.' ') then
      if(nnoe.ne.0.and.opt.eq.0) call penalty
      call disth
      endif
      if(index(test,'S').ne.0) then
      call grads(nspc,nbcc,nhlc,nklc,nlgc,rcom)
      else if(test.ne.' ') then
      call gradt
      endif
      return
      end
