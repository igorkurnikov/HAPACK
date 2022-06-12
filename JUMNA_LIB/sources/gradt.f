      subroutine gradt
      include 'jumna_data.inc'
      logical*2 ifhb,lar,lock,kink,lthy,ribose,cation,locr,sup,
     1 rcom,homo,homo2,homo3,diep,link,ecen,cyl,lcat,cent,autos,
     1 ihl,amber,sum
      character*4 seq*120,code*8,kode*8,han(6)*5,kan(4)*4,mac*32,
     1 lmo*32,lib*32,libm*32,out*32,axe*32,noe*32,nol*32,axl*32,
     1 test*32,pdb*32,ins*32,bar*32,parm*32,lnam,mnam,munit
      integer*2 i23,i34,elim
      integer*4 opt
      common/cha/mac,lmo,lib,libm,out,axe,axl,noe,nol,test,pdb,ins,
     1 bar,parm
      common/datjm/acc,phos,epsi,epsr,plat,slope,rhbl,vfac,tfac,rfac,xfac,
     1 scale,rpiv,tpiv,fad,fan,damp,fnoew,fnoes,fnoea,df1,scnb,scee,
     1 catd,catr,catc,rad,pit,enit,nshel,maxn,opt,mhomo,mhomo2,mhomo3,
     1 nick,limit,nop,nrib,ncat,lig,naxo,sup,rcom,homo,homo2,homo3,diep,
     1 sum,link,ecen,cyl,lcat,cent,autos,amber
      common/cnd/rcyl,dcyl,hcyl,tcyl,eneq,ucyl(3),pcyl(3),ecyl,
     1 ac(25),bc(25),lcyl(3)
      common/dcu/dpx,dpy,dpz,dux,duy,duz,crx,cry,crz,ctx,cty,ctz,
     1 csx,csy,csz,crs,clx,cly,crl
      common/enf/for(n1,3),tor(n1,3),fot(n6),ener,elec,
     1 repl,disp,eang,etog,epen,i23(8*n1),i34(8*n1),elim(n1)
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/flx/sap(n6),refg,refb,refh,refv,iap(n8,n4,n5),
     1 nap(n8,7,n5),kap(3,n5)
      common/ind/iend(0:4),ise(n2),ienb(n2,2),kapt(n2+1),nsph
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
      common/srs/corms(n1,3),saps(n6),vars(n7),gras(n7),hels(n2,6),
     1 vkis(n2,4),has(n2,9),rlis(n9,6),eref,rref,pref
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(120,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      common/symjm/xsum,ysum,itwl(n2),nvs(n7),ihm(n2),ksym(3),nbrk(2),
     1 isym,isur,isup,ihl(n2)
      data han/'Xdisp','Ydisp','Rise','Inc','Tip','Twist'/,
     1     kan/'Ax','Ay','Ainc','Atip'/
      koff=nto-nst
      if(nst.eq.1) then
      write(6,6)
6     format(/2x,'---- GRADTEST FOR SINGLE STRAND ----')
      else
      write(6,7)
7     format(/2x,'---- GRADTEST FOR MULTI-STRAND ----')
      endif
         if(rcom) then
         call fitsug(1)
         call fitsug(0)
         endif
      call microb
      if(parm.eq.'Flex') then
      call energy
      else if(parm.eq.'Amber91') then
      call energ91
      else if(parm.eq.'Amber94') then
      call energ94
      endif
      call assemb
         if(rcom) then
         call fitsug(1)
         call fitsug(0)
         endif
      call putbac(0)
      write(6,21) repl,disp,repl+disp,epen,elec,eang,etog,ener-epen
21    format(
     1 /2x,'REPL= ',f11.3,' DISP= ',f11.3,' LJ  = ',f11.3,
     1 ' PEN = ',f11.3,/2x,'ELEC= ',f11.3,' ANGL= ',f11.3,
     1 ' TORG= ',f11.3,' TOT = ',f11.3/)
c-------------------------------------------------------------test atomic forces
      if(index(test,'F').ne.0) then
      j=1
      if(index(test,'Y').ne.0) j=2
      if(index(test,'Z').ne.0) j=3
      do i=1,kam
      corm(i,j)=corm(i,j)+deltv
      if(parm.eq.'Flex') then
      call energy
      else if(parm.eq.'Amber91') then
      call energ91
      else if(parm.eq.'Amber94') then
      call energ94
      endif
      grd=(eref-ener)/deltv
      corm(i,j)=corm(i,j)-deltv
      if(abs(grd).gt.1.d-20.or.abs(for(i,j)).gt.1.d-20) write(6,25)
     1 munit(i),nunit(i),mnam(i),for(i,j),grd,grd-for(i,j)
25    format(2x,a4,i3,1x,a4,4x,3f14.2)
      enddo
      endif
c-------------------------------------------------------------test var gradients
      if(index(test,'B').eq.0) then
      k=kapt(nto+1)
      i=nsph
      goto 100
      endif
      write(6,20)
20    format(/2x,'-----------------------------------------------'/)
      k=0
      i=0
      do is=1,nto
      in=itr(is)
      do l=1,kap(1,in)
      k=k+1
      if(.not.lock(k)) then
      i=i+1
      if(.not.rcom) then
      var(i)=var(i)+deltv
         else
         if(neq(i).eq.0) goto 10
         call fitsug(1)
         var(i)=var(i)+deltv
         call fitsug(0)
         endif
      call microb
      if(parm.eq.'Flex') then
      call energy
      else if(parm.eq.'Amber91') then
      call energ91
      else if(parm.eq.'Amber94') then
      call energ94
      endif
      grd=(ener-eref)/deltv
      write(6,30) is,l,gras(i),grd,gras(i)-grd
30    format(2x,'K= ',i2,' L= ',i2,' GRA= ',3f10.4)
         if(rcom) then
         call fitsug(1)
         var(i)=var(i)-deltv
         call fitsug(0)
         endif
      call putbac(1)
      endif
10    enddo
      k=k+kap(2,in)-kap(1,in)
      enddo
c----------------------------------------------------test intra ligand gradients
100   if(index(test,'L').eq.0) then
      k=ntba
      i=nbac
      goto 110
      endif
      if(nlig.ne.0.and.ntba.ne.kapt(nto+1)) write(6,20)
      do il=1,nlig
      is=nto+il
      in=itr(is)
      do l=1,kap(1,in)
      k=k+1
      if(.not.lock(k)) then
      i=i+1
      var(i)=var(i)+deltv
      call microb
      if(parm.eq.'Flex') then
      call energy
      else if(parm.eq.'Amber91') then
      call energ91
      else if(parm.eq.'Amber94') then
      call energ94
      endif
      grd=(ener-eref)/deltv
      write(6,33) il,l,gras(i),grd,gras(i)-grd
33    format(2x,'A= ',i2,' L= ',i2,' GRA= ',3f10.4)
      call putbac(1)
      endif
      enddo
      k=k+kap(2,in)-kap(1,in)
      enddo
c---------------------------------------------------------test helical gradients
110   if(index(test,'H').eq.0) then
      k=nthe
      i=nhel
      goto 200
      endif
      write(6,20)
      do is=1,nto
      in=itr(is)
      do l=1,6
      k=k+1
      if(.not.lock(k)) then
      i=i+1
      var(i)=var(i)+deltv
      call microb
      if(parm.eq.'Flex') then
      call energy
      else if(parm.eq.'Amber91') then
      call energ91
      else if(parm.eq.'Amber94') then
      call energ94
      endif
      grd=(ener-eref)/deltv
      write(6,34) is,han(l),gras(i),grd,gras(i)-grd
34    format(2x,'K= ',i2,' L=', a5,' GRA= ',3f10.4)
      call putbac(1)
      endif
      enddo
      enddo
c---------------------------------------------------------------------test kinks
200   if(index(test,'K').eq.0) then
      k=ntki
      i=nkin
      goto 300
      endif
      write(6,20)
      do is=1,kseq
      if(kink(is)) then
      do l=1,4
      k=k+1
      if(.not.lock(k)) then
      i=i+1
      var(i)=var(i)+deltv
      call microb
      if(parm.eq.'Flex') then
      call energy
      else if(parm.eq.'Amber91') then
      call energ91
      else if(parm.eq.'Amber94') then
      call energ94
      endif
      grd=(ener-eref)/deltv
      write(6,38) is,kan(l),gras(i),grd,gras(i)-grd
38    format(2x,'K= ',i2,' L=', a4,' GRA= ',3f10.4)
      call putbac(1)
      endif
      enddo
      endif
      enddo
c-------------------------------------------------------------test inter ligands
300   if(index(test,'I').eq.0) goto 310
      if(nlig.ne.0) write(6,20)
      do il=1,nlig
      is=nto+il
      kal=nuc(is)-nuc(is-1)
      jlim=6
      if(kal.eq.1) jlim=3
      do j=1,jlim
      k=k+1
      if(.not.lock(k)) then
      i=i+1
      var(i)=var(i)+deltv
      call microb
      if(parm.eq.'Flex') then
      call energy
      else if(parm.eq.'Amber91') then
      call energ91
      else if(parm.eq.'Amber94') then
      call energ94
      endif
      grd=(ener-eref)/deltv
      write(6,39) il,han(j),gras(i),grd,gras(i)-grd
39    format(2x,'R= ',i2,' L=', a5,' GRA= ',3f10.4)
      call putbac(1)
      endif
      enddo
      enddo
c----------------------------------------------------------------test superhelix
310   if(index(test,'U').eq.0) goto 320
      write(6,20)
      if(isur.gt.0) then
      var(isur)=var(isur)+deltv
      call microb
      if(parm.eq.'Flex') then
      call energy
      else if(parm.eq.'Amber91') then
      call energ91
      else if(parm.eq.'Amber94') then
      call energ94
      endif
      grd=(ener-eref)/deltv
      write(6,37) 'Rad',isur,gras(isur),grd,gras(isur)-grd
37    format(2x,'S= ',a3,2x,i5,' GRA= ',3f10.4)
      call putbac(1)
      endif
      if(isup.gt.0) then
      var(isup)=var(isup)+deltv
      call microb
      if(parm.eq.'Flex') then
      call energy
      else if(parm.eq.'Amber91') then
      call energ91
      else if(parm.eq.'Amber94') then
      call energ94
      endif
      grd=(ener-eref)/deltv
      write(6,37) 'Pit',isup,gras(isup),grd,gras(isup)-grd
      call putbac(1)
      endif
c------------------------------------------------------------------test cylinder
320   if(index(test,'C').eq.0.or..not.cyl) return
      write(6,20)
      call microb
      do i=1,3
      iv=lcyl(i)
      if(iv.ne.0) then
      var(iv)=var(iv)+deltv
      if(i.eq.1) tcyl=var(iv)
      if(i.eq.2) dcyl=var(iv)
      if(i.eq.3) hcyl=var(iv)
      ct=cos(cdr*(tcyl))
      st=sin(cdr*(tcyl))
      pcyl(1)=dcyl*ct
      pcyl(2)=dcyl*st
      pcyl(3)=hcyl
      ucyl(1)= st
      ucyl(2)=-ct
      ucyl(3)=0.
      if(parm.eq.'Flex') then
      call energy
      else if(parm.eq.'Amber91') then
      call energ91
      else if(parm.eq.'Amber94') then
      call energ94
      endif
      grd=(ener-eref)/deltv
      write(6,48) i,gras(iv),grd,gras(iv)-grd
48    format( /2x,'CYL= ',i1,4x,' GRA= ',3f10.4)
      var(iv)=var(iv)-deltv
      if(i.eq.1) tcyl=var(iv)
      if(i.eq.2) dcyl=var(iv)
      if(i.eq.3) hcyl=var(iv)
      ct=cos(cdr*(tcyl))
      st=sin(cdr*(tcyl))
      pcyl(1)=dcyl*ct
      pcyl(2)=dcyl*st
      pcyl(3)=hcyl
      ucyl(1)= st
      ucyl(2)=-ct
      ucyl(3)=0.
      endif
      enddo
      return
      end
