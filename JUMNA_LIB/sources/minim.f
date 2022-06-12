      subroutine minim
      include 'jumna_data.inc'
      logical*2 lthy,ifhb,lar,lock,kink,ribose,cation,sup,rcom,homo,
     1 homo2,homo3,diep,link,ecen,cyl,lcat,cent,autos,locr,tlp,sum,
     1 quiet,convg,dloop,amber
      character*4 mnam,munit,seq*120,code*8,kode*8,lnam,mac*32,
     1 lmo*32,lib*32,libm*32,out*32,axe*32,noe*32,nol*32,axl*32,
     1 test*32,pdb*32,ins*32,bar*32,parm*32
      integer*2 i23,i34,elim
      integer*4 opt
      real*4 dtime,tarray(2)
      dimension w(n7*(n7+13)/2),vrc(n7),grc(n7)
      common/cha/mac,lmo,lib,libm,out,axe,axl,noe,nol,test,pdb,ins,
     1 bar,parm
      common/datjm/acc,phos,epsi,epsr,plat,slope,rhbl,vfac,tfac,rfac,xfac,
     1 scale,rpiv,tpiv,fad,fan,damp,fnoew,fnoes,fnoea,df1,scnb,scee,
     1 catd,catr,catc,rad,pit,enit,nshel,maxn,opt,mhomo,mhomo2,mhomo3,
     1 nick,limit,nop,nrib,ncat,lig,naxo,sup,rcom,homo,homo2,homo3,diep,
     1 sum,link,ecen,cyl,lcat,cent,autos,amber
      common/cnd/rcyl,dcyl,hcyl,tcyl,eneq,ucyl(3),pcyl(3),ecyl,
     1 ac(25),bc(25),lcyl(3)
      common/enf/for(n1,3),tor(n1,3),fot(n6),ener,elec,
     1 repl,disp,eang,etog,epen,i23(8*n1),i34(8*n1),elim(n1)
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/flx/sap(n6),refg,refb,refh,refv,iap(n8,n4,n5),
     1 nap(n8,7,n5),kap(3,n5)
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
      common/parjm/aij(25,25),bij(25,25),ahb(15),bhb(15),vt(35),va(25),
     1 vo(25),delmax,time0,nzsh,nwdg,ivf(35),ifhb(25,25)
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(120,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      convg=.false.
c------------------------------------------------------------------move cylinder
      if(cyl) then
        if(lcyl(1).ne.0) tcyl=var(lcyl(1))
        if(lcyl(2).ne.0) dcyl=var(lcyl(2))
        if(lcyl(3).ne.0) hcyl=var(lcyl(3))
        ct=cos(cdr*(tcyl))
        st=sin(cdr*(tcyl))
        pcyl(1)=dcyl*ct
        pcyl(2)=dcyl*st
        pcyl(3)=hcyl
        ucyl(1)= st
        ucyl(2)=-ct
        ucyl(3)=0.
      endif
c------------------------------------------------------------------------startup
      call microb
      if(parm.eq.'Flex') then
        call energy
      else if(parm.eq.'Amber91') then
        call energ91
      else if(parm.eq.'Amber94') then
        call energ94
      endif
      call assemb
      esum=ener
c-----------------------------------------------------------symmetry contraction
      if(rcom) call fitsug(1)
      do k=1,nvrc
        grc(k)=0.
      enddo
      do i=1,nvar
        if(i.le.nbac) then
          fac=tfac
          if(lar(i)) fac=vfac
        else if(i.le.nlgi) then
          fac=xfac
          if(lar(i)) fac=rfac
        else
        fac=tfac
        endif
        if(i.eq.lcyl(1)) fac=rfac
        if(i.eq.lcyl(2).or.i.eq.lcyl(3)) fac=xfac
        k=abs(neq(i))
        if(neq(i).gt.0) then
          vrc(k)=var(i)
          scl(k)=fac*scale
        endif
        if(k.ne.0) grc(k)=grc(k)+gra(i)
      enddo
c--------------------------------------------------------------------------cycle
      if(.not.quiet) then
        write(6,21) repl,disp,repl+disp,epen,elec,eang,
     1   etog,ener-epen
21    format(
     1 /2x,'REPL= ',f11.3,' DISP= ',f11.3,' LJ  = ',f11.3,
     1 ' PEN = ',f11.3,/2x,'ELEC= ',f11.3,' ANGL= ',f11.3,
     1 ' TORG= ',f11.3,' TOT = ',f11.3/)
      write(6,403) (grc(i),i=1,nvrc)
      write(6,*) ' '
      endif
c----------------------------------------------------------------------minimizer
      deltim=dtime(tarray)
      time0=time0+deltim
      if(limit.ne.0) limit=limit-deltim
      delmax=0.
      nd=1+(nvrc*(nvrc+1))/2
      nw=nd+nvrc
      nxa=nw+nvrc
      nga=nxa+nvrc
      nxb=nga+nvrc
      ngb=nxb+nvrc
      call minfor(nvrc,vrc,esum,grc,scl,acc,w,w(nd),w(nw),
     1 w(nxa),w(nga),w(nxb),w(ngb),maxn,nfun)
      if(nfun.lt.maxn) then
      convg=.true.
      if(.not.quiet) write(6,401)
401   format(/2x,'---- CONVERGENCE ----'/)
      else
      if(.not.quiet) write(6,402)
402   format(/2x,'---- STEP LIMIT ----'/)
      endif
         gmx=0.
         do i=1,nvrc
           if(abs(grc(i)).gt.gmx) gmx=abs(grc(i))
         enddo
      if(.not.quiet) write(6,403) (grc(i),i=1,nvrc)
403   format(2x,'G: ',8f8.3)
      return
      end
