      subroutine combi
      include 'jumna_data.inc'
      parameter (n10=25)
      logical*2 lar,lock,lthy,kink,ifhb,hst,bst,vst,lgi,lgj,lgu,
     1 ribose,cation,sup,rcom,homo,homo2,homo3,diep,link,ecen,cyl,
     1 sum,lcat,cent,autos,locr,tlp,quiet,convg,cstop,dloop,amber
      character*4 mnam,munit,seq*120,code*8,kode*8,inaxe*40,
     1 mac*32,lmo*32,lib*32,libm*32,out*32,axe*32,noe*32,nol*32,axl*32,
     1 test*32,pdb*32,ins*32,bar*32,parm*32,lnam,knam,
     1 axnm(5000)*10
      logical*4 there
      integer*2 i23,i34,elim
      integer*4 opt
      dimension rsav(n3,3),isav(n3,10),ivec(n2)
      common/cha/mac,lmo,lib,libm,out,axe,axl,noe,nol,test,pdb,ins,
     1 bar,parm
      common/comjm/cut,icx(n2),icn(n2),ict(n2,10),ncut,nmx,mco,ncomb,icomb
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
      nmx=maxn
c----------------------------------------------------------------save combi data
      do i=nnoe+1,icomb
        rsav(i,1)=rnoe(i)
        rsav(i,2)=bnoe(i)
        rsav(i,3)=fnoe(i)
        do j=1,9
          isav(i,j)=inoe(i,j)
        enddo
        isav(i,10)=jnoe(i)
      enddo
c-------------------------------------------------------------setup combi vector
      icomb=nnoe
      do i=1,mco
        ivec(i)=1
        inoe(icomb+i,1)=icx(i)
      enddo
c---------------------------------------------------------------------combi loop
      cstop=.false.
      kfo=index(out,' ')-1
      km=1
      kc=0
1     kc=kc+1
      nnoe=icomb
      axnm(kc)=' '
      do i=1,mco
        k=ivec(i)
        write(axnm(kc)(i:i),'(i1)') k
        n=ict(i,k)
        nnoe=nnoe+1
        rnoe(nnoe)=rsav(n,1)
        bnoe(nnoe)=rsav(n,2)
        fnoe(nnoe)=rsav(n,3)
        do j=2,9
          inoe(nnoe,j)=isav(n,j)
        enddo
        jnoe(nnoe)=isav(n,10)
      enddo
c--------------------------------------------------------------check for restart
      inquire(file=out(:kfo)//'-'//axnm(kc)(:mco)//'.axe',exist=there)
      if(there) goto 12
c-------------------------------------------------------------reset and minimize
      if(kc.eq.1) then
         inaxe='initial'
         kfi=7
      else
         inaxe=out(:kfo)//'-'//axnm(kc-km)
         kfi=kfo+mco+1
         call reset(inaxe)
      endif
      cut=0.1
      ncut=0
      call minim
      icyo=icy
      nnoe=icomb
      cut=-1.0
      ncut=0
      maxn=nmx
      call minim
c-------------------------------------------------------------------------output
      write(6,14) inaxe(:kfi),out(:kfo)//'-'//axnm(kc)(:mco)
14    format(/2x,'%%%%%%%%%% From: ',a,' to ',a,' %%%%%%%%%%')
      if(convg) then
         write(6,20) icyo,icy,gmx
20       format(/2x,'Converged in ',i3,'+',i3,' cycles, gmax= ',f6.3)
      else
         write(6,22) icyo,icy,gmx
22       format(2x,'Failed to converge in ',i3,'+',i3,
     1   ' cycles, gmax= ',f6.3)
         cstop=.true.
      endif
      call helix
      call backbo
      call penalty
      call disth
      if(parm.eq.'Flex')then
        call ecomp
      else if(parm.eq.'Amber91') then
        call ecomp91
      else if(parm.eq.'Amber94') then
        call ecomp94
      endif
      write(6,21) repl,disp,repl+disp,epen,elec,eang,etog,ener-epen
21    format(/2x,'Corrected final energy ...'
     1//2x,'REPL= ',f11.3,' DISP= ',f11.3,' LJ  = ',f11.3,
     1' PEN = ',f11.3,/2x,'ELEC= ',f11.3,' ANGL= ',f11.3,
     1' TORG= ',f11.3,' TOT = ',f11.3)
c     if(cstop) stop
      call axeout(out(:kfo)//'-'//axnm(kc))
c------------------------------------------------------------------------step on
12    quiet=.true.
      ist=1
10    ivec(ist)=ivec(ist)+1
      if(ivec(ist).gt.icn(ist)) then
        ivec(ist)=1
        ist=ist+1
        if(ist.gt.mco) goto 100
        goto 10
      endif
      km=1
      do i=2,ist
        km=km*icn(i-1)
      enddo
      goto 1
100   return
      end
