      subroutine loops
      include 'jumna_data.inc'
      parameter (n10=8000)
      logical*2 lar,lock,lthy,kink,ifhb,hst,bst,vst,lgi,lgj,lgu,
     1 ribose,cation,sup,rcom,homo,homo2,homo3,diep,link,ecen,cyl,
     1 lcat,cent,autos,locr,tlp,quiet,convg,ihl,dloop,amber,sum
      logical*4 there
      character*4 mnam,munit,seq*120,code*8,kode*8,inaxe*40,olaxe*10,
     1 mac*32,lmo*32,lib*32,libm*32,out*32,axe*32,noe*32,nol*32,axl*32,
     1 test*32,pdb*32,ins*32,bar*32,parm*32,lnam,knam,
     1 axnm(n10)*10,ic*3,jc*3
      integer*2 i23,i34,elim
      integer*4 opt
      dimension nvc(200,2),nl(2)
      common/cha/mac,lmo,lib,libm,out,axe,axl,noe,nol,test,pdb,ins,
     1 bar,parm
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
      olaxe=' '
      itws=0
      if(itwl(1).gt.0) itws=abs(neq(itwl(1)))
      np1=lphig(1)-lplow(1)+1
      lp1=lplow(1)
      np2=0
      lp2=0
         if(nloop.eq.2) then
         np2=lphig(2)-lplow(2)+1
         lp2=lplow(2)
         endif
c--------------------------------------------------------------------setup names
      do i=lplow(1),lphig(1)
      if(abs(i).le.9) then
      in=1
      write(ic,5) abs(i)
 5    format(i1)
      else if(abs(i).le.99) then
      in=2
      write(ic,10) abs(i)
 10   format(i2)
      else
      in=3
      write(ic,15) abs(i)
 15   format(i3)
      endif
      i0=np1*( -lp2)+(i-lp1)+1
      if(i.lt.0) then
      axnm(i0)='._'//ic(:in)
      in=in+2
      else
      axnm(i0)='.'//ic(:in)
      in=in+1
      endif
        if(nloop.eq.2) then
         do j=lplow(2),lphig(2)
         if(abs(j).le.9) then
         write(jc,5) abs(j)
         else if(abs(j).le.99) then
         write(jc,10) abs(j)
         else
         write(jc,15) abs(j)
         endif
         i0=np1*( -lp2)+(i-lp1)+1
         ij=np1*(j-lp2)+(i-lp1)+1
         if(j.lt.0) then
         axnm(ij)=axnm(i0)(:in)//'._'//jc
         else
         axnm(ij)=axnm(i0)(:in)//'.'//jc
         endif
         enddo
         endif
      enddo
c-----------------------------------------------------------setup neq compaction
      do k=1,nloop
      if(.not.tlp(k)) then
      l=0
      if(k.eq.2.and.indlp(2,1).gt.n) indlp(2,1)=indlp(2,1)-1
      n=abs(indlp(k,1))
      do i=1,nvar
      in=neq(i)
      ina=abs(in)
      if(ina.eq.n) then
      l=l+1
      nvc(l,k)=i
      neq(i)=0
         else if(ina.gt.n) then
         neq(i)=in-sign(1,in)
         endif
      enddo
      nl(k)=l
      nvrc=nvrc-1
      endif
      enddo
      start1=stlp(1)
      start2=stlp(2)
c--------------------------------------------------------------------setup loops
      del1=dellp(1)
      if(tlp(1)) then
      do nc1=indlp(1,1),indlp(1,2)
      k=jnoe(nc1)
      if(k.ne.20.and.k.ne.22) then
      if(nc1.eq.indlp(1,1)) start1=rnoe(nc1)
      else
      write(6,11)
 11   format(/2x,'---- Loop 1 applied to wrong noe type ----')
      stop
      endif
      enddo
      endif
         if(nloop.eq.2) then
         del2=dellp(2)
         if(tlp(2)) then
         do nc2=indlp(2,1),indlp(2,2)
         k=jnoe(nc2)
         if(k.ne.20.and.k.ne.22) then
         if(nc2.eq.indlp(2,1)) start2=rnoe(nc2)
         else
         write(6,18)
 18      format(/2x,'---- Loop 2 applied to wrong noe type ----')
         stop
         endif
         enddo
         endif
            else
            np2=1
            endif
      kfo=index(out,' ')-1
      kfi=kfo+8
c--------------------------------------------------------------------begin loops
      j=0
      jd=-1
      do k2=1,np2
         i=0
         id=-1
         do k1=1,np1
         if(k1.gt.1) then
         im=i-sign(1,i)
         jm=j
         else
            im=i
            jm=j-sign(1,j)
            endif
         r1=start1+i*del1
         if(nloop.eq.2) r2=start2+j*del2
         if(tlp(1)) then
         do nc1=indlp(1,1),indlp(1,2)
         rnoe(nc1)=r1
         enddo
         endif
            if(tlp(2)) then
            do nc2=indlp(2,1),indlp(2,2)
            rnoe(nc2)=r2
            enddo
            endif
c--------------------------------------------------------------check for restart
      ij=np1*(j-lp2)+(i-lp1)+1
      iax=index(axnm(ij),' ')-1
      inquire(file=out(:kfo)//axnm(ij)(:iax)//'.axe',exist=there)
      if(there) goto 12
c-------------------------------------------------------------------------------
         if(i.eq.0.and.j.eq.0) then
         inaxe='initial'
         else
         imjm=np1*(jm-lp2)+(im-lp1)+1
         inaxe=out(:kfo)//axnm(imjm)
         if(olaxe.ne.axnm(imjm)) call reset(inaxe)
         quiet=.true.
         endif
           if(.not.tlp(1)) then
               if(isur.ne.0.and.abs(indlp(1,1)).eq.itws) then
               del=r1-var(nvc(1,1))
               do m=1,nto
               mm=itwl(m)
               var(mm)=var(mm)+del
               enddo
               endif
            do m=1,nl(1)
            var(nvc(m,1))=r1
            enddo
            endif
         if(nloop.eq.2) then
            if(.not.tlp(2)) then
               if(isur.ne.0.and.abs(indlp(2,1)).eq.itws) then
               del=r2-var(nvc(1,2))
               do m=1,nto
               mm=itwl(m)
               var(mm)=var(mm)+del
               enddo
               endif
            do m=1,nl(2)
            var(nvc(m,2))=r2
            enddo
            endif
         endif
         call minim
         if(nloop.eq.1) then
         write(6,14) inaxe(:kfi),out(:kfo)//axnm(ij),r1
 14      format(/2x,'%%%%%%%%%% From: ',a,' to ',a,' Con: ',f8.3,
     1   ' %%%%%%%%%%')
         else
         write(6,16) inaxe(:kfi),out(:kfo)//axnm(ij),r1,r2
 16      format(/2x,'%%%%%%%%%% From: ',a,' to ',a,' Con: ',f8.3,
     1   '/',f8.3,' %%%%%%%%%%')
         endif
         if(convg) then
         write(6,20) icy,gmx
 20      format(/2x,'Converged in ',i4,' cycles, gmax= ',f6.3)
            else
            write(6,22) icy,gmx
 22         format(2x,'---- Failed to converge in ',i4,
     1                ' cycles, gmax= ',f6.3,' ----'/)
            stop
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
 21      format(/2x,'Corrected final energy ...'
     1   //2x,'REPL= ',f11.3,' DISP= ',f11.3,' LJ  = ',f11.3,
     1   ' PEN = ',f11.3,/2x,'ELEC= ',f11.3,' ANGL= ',f11.3,
     1   ' TORG= ',f11.3,' TOT = ',f11.3)
         call axeout(out(:kfo)//axnm(ij))
      olaxe=axnm(ij)
 12   i=i+id
      if(i.lt.lplow(1).or.i.gt.lphig(1)) then
      id=-id
      i=id
      endif
      enddo
         j=j+jd
         if(j.lt.lplow(2).or.j.gt.lphig(2)) then
         jd=-jd
         j=jd
         endif
         enddo
      return
      end
