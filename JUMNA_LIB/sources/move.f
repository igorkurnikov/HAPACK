      subroutine movejm(icyc,vrc,esum,grc)
      include 'jumna_data.inc'
      logical*2 lthy,ifhb,lar,lock,kink,ribose,cation,sup,rcom,homo,
     1 homo2,homo3,diep,link,ecen,cyl,lcat,cent,autos,tlp,quiet,
     1 convg,dloop,amber,sum
      character*4 seq*120,code*8,kode*8,knam,ichr*4,
     1 mac*32,lmo*32,lib*32,libm*32,out*32,axe*32,noe*32,nol*32,axl*32,
     1 test*32,pdb*32,ins*32,bar*32,parm*32
      integer*2 i23,i34,elim
      integer*4 opt
      real*8 delta
      real*4 dtime,tarray(2)
      dimension vrc(n7),grc(n7)
      common/cha/mac,lmo,lib,libm,out,axe,axl,noe,nol,test,pdb,ins,
     1 bar,parm
      common/cnd/rcyl,dcyl,hcyl,tcyl,eneq,ucyl(3),pcyl(3),ecyl,
     1 ac(25),bc(25),lcyl(3)
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
      common/lop/dellp(2),stlp(2),gmx,indlp(2,2),lplow(2),lphig(2),
     1 nloop,icy,tlp(2),quiet,convg,dloop
      common/mnn/var(n7),gra(n7),con(n0*n2+n3),scl(n7),ncon,nvrc,
     1 ntba,nbac,nthe,nhel,ntki,nkin,ntri,ntot,nvar,nrin,
     1 lar(n7),lock(n6a)
      common/nmrjm/rnoe(n3),bnoe(n3),fnoe(n3),snoe,knam(200,2),ktyp(200),
     1 icon(n2*2),inoe(n3,9),ih68(n2),jnoe(n3),nnoe,ndcs,ndiv(0:99),nlin
      common/parjm/aij(25,25),bij(25,25),ahb(15),bhb(15),vt(35),va(25),
     1 vo(25),delmax,time0,nzsh,nwdg,ivf(35),ifhb(25,25)
      common/srs/corms(n1,3),saps(n6),vars(n7),gras(n7),hels(n2,6),
     1 vkis(n2,4),has(n2,9),rlis(n9,6),eref,rref,pref
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
c-------------------------------------------------------------symmetry expansion
c      if(minz.eq.'m1qn3') then ! Newton minimization
c        do i=1,nvar
c          k=abs(neq(i))
c          if(k.ne.0) var(i)=vrc(k)/scl(k)
c        enddo
c      else  // if Newton minimization
        do i=1,nvar
          k=abs(neq(i))
          if(k.ne.0) var(i)=vrc(k)
        enddo
c      endif
      if(rcom) call fitsug(0)
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
c---------------------------------------------------------------------------step
      if(icyc.eq.1) eref=ener
      call microb
      if(parm.eq.'Flex') then
        call energy
      else if(parm.eq.'Amber91') then
        call energ91
      else if(parm.eq.'Amber94') then
        call energ94
      endif
      call assemb
      delta=ener-eref
      eref=ener
      esum=ener
c-----------------------------------------------------------symmetry contraction
      if(rcom) call fitsug(1)
      do k=1,nvrc
        grc(k)=0.
      enddo
c      if(minz.eq.'m1qn3') then ! Newton minimization routine
c        do i=1,nvar
c          k=abs(neq(i))
c          if(neq(i).gt.0) vrc(k)=var(i)*scl(k)
c          if(k.ne.0) grc(k)=grc(k)+gra(i)/scl(k)
c        enddo
c      else  
        do i=1,nvar
          k=abs(neq(i))
          if(neq(i).gt.0) vrc(k)=var(i)
          if(k.ne.0) grc(k)=grc(k)+gra(i)
        enddo
c      endif
c------------------------------------------------------------------max gradients
      gx=0.
      do k=1,nvrc
        g=abs(grc(k))
        if(g.gt.gx) gx=g
      enddo
      if(gx.lt.cut) then
         ncut=ncut+1
      else
         ncut=0
      endif
      if(ncut.eq.5) maxn=icyc
c-------------------------------------------------------------------------timing
      deltim=dtime(tarray)
      time0=time0+deltim
      if(deltim.gt.delmax) delmax=deltim
      if(limit.ne.0) then
        limit=limit-deltim
        if(3*delmax.gt.limit) maxn=icyc
      endif
c---------------------------------------------------------------------axe output
      if(nloop.eq.0.and.naxo.gt.0.and.mod(icyc,naxo).eq.0) then
        kfo=index(out,' ')-1
        ic=icyc/naxo
        in=4-int(log10(float(ic)))
        write(ichr,'(i4)') ic
        call axeout(out(:kfo)//'.'//ichr(in:)//' ')
      endif
c--------------------------------------------------------------------------cycle
      if(.not.sup.and..not.quiet) write(6,*) ' '
      if(.not.quiet) write(6,20) icyc,ener-epen,epen,delta,gx,deltim
20    format(2x,i3,') TOT= ',f11.3,'/',f11.3,' DEL= ',e11.5,' G= ',f6.1,
     1 ' T= ',f5.1)
      if(.not.sup.and..not.quiet) then
         write(6,22) repl+disp,elec,eang,etog,epen
22       format(2x,'L=',f11.3,' E=',f11.3,' A=',f11.3,
     1' T=',f11.3,' P=',f11.3)
      endif
      icy=icyc
      return
      end
