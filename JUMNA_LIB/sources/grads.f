      subroutine grads(nspc,nbcc,nhlc,nklc,nlgc,rcom)
      include 'jumna_data.inc'
      logical*2 ifhb,lar,lock,kink,lthy,ribose,cation,locr,rcom
      character*4 seq*120,code*8,kode*8,lnam,
     1 mac*32,lmo*32,lib*32,libm*32,out*32,axe*32,noe*32,nol*32,axl*32,
     1 test*32,pdb*32,ins*32,bar*32,parm*32
      integer*2 i23,i34,elim
      dimension grc(n7)
      common/cha/mac,lmo,lib,libm,out,axe,axl,noe,nol,test,pdb,ins,
     1 bar,parm
      common/enf/for(n1,3),tor(n1,3),fot(n6),ener,elec,
     1 repl,disp,eang,etog,epen,i23(8*n1),i34(8*n1),elim(n1)
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/flx/sap(n6),refg,refb,refh,refv,iap(n8,n4,n5),
     1 nap(n8,7,n5),kap(3,n5)
      common/lgd/rlig(n9,6),slig(n9,2),ilig(n9,6),lopt(n9),lpiv(n9),
     1 ntlg,nlgi,nlig,ntl,lnam(n9),locr(n2,n0)
      common/mnn/var(n7),gra(n7),con(n0*n2+n3),scl(n7),ncon,nvrc,
     1 ntba,nbac,nthe,nhel,ntki,nkin,ntri,ntot,nvar,nrin,
     1 lar(n7),lock(n6a)
      common/parjm/aij(25,25),bij(25,25),ahb(15),bhb(15),vt(35),va(25),
     1 vo(25),delmax,time0,nzsh,nwdg,ivf(35),ifhb(25,25)
      common/srs/corms(n1,3),saps(n6),vars(n7),gras(n7),hels(n2,6),
     1 vkis(n2,4),has(n2,9),rlis(n9,6),eref,rref,pref
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      koff=nto-nst
      if(nst.eq.1) then
      write(6,6)
6     format(/2x,'---- CONTRACTED GRADTEST FOR SINGLE STRAND ----')
      else
      write(6,7)
7     format(/2x,'---- CONTRACTED GRADTEST FOR MULTI-STRAND ----')
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
      call putbac(0)
      write(6,21) repl,disp,repl+disp,epen,elec,eang,etog,ener-epen
21    format(
     1 /2x,'REPL= ',f11.3,' DISP= ',f11.3,' LJ  = ',f11.3,
     1 ' PEN = ',f11.3,/2x,'ELEC= ',f11.3,' ANGL= ',f11.3,
     1 ' TORG= ',f11.3,' TOT = ',f11.3/)
c-----------------------------------------------------------symmetry contraction
      if(rcom) call fitsug(1)
      do k=1,nvrc
        grc(k)=0.
      enddo
c      if(minz.eq.'m1qn3') then  ! Newton minimization:
c        do i=1,nvar
c          k=abs(neq(i))
c          if(k.ne.0) grc(k)=grc(k)+gra(i)/scl(k)
c        enddo
c      else   
        do i=1,nvar
          k=abs(neq(i))
          if(k.ne.0) grc(k)=grc(k)+gra(i)
        enddo
c      endif
c-------------------------------------------------------------test var gradients
      if(index(test,'B').eq.0) goto 100
      write(6,20)
20    format(/2x,'-----------------------------------------------'/)
      do k=1,nspc
        do i=1,nvar
          l=abs(neq(i))
          if(l.eq.k) then
c           if(minz.eq.'m1qn3') then ! Newton minimization:
c             var(i)=var(i)+deltv/scl(k)
c           else 
              var(i)=var(i)+deltv
c           endif
          endif
        enddo
        if(rcom) call fitsug(0)
        call microb
        if(parm.eq.'Flex') then
          call energy
        else if(parm.eq.'Amber91') then
          call energ91
        else if(parm.eq.'Amber94') then
          call energ94
      endif
      grd=(ener-eref)/deltv
      write(6,30) k,grc(k),grd,grc(k)-grd
30    format(2x,'K= ',i3,' GRA= ',3f10.4)
      if(rcom) call fitsug(1)
      call putbac(1)
      enddo
c----------------------------------------------------test intra ligand gradients
100   if(index(test,'L').eq.0) goto 110
      if(nspc.lt.nbcc) write(6,20)
      do k=nspc+1,nbcc
      do i=1,nvar
      l=abs(neq(i))
      if(l.eq.k) then  
c         if(minz.eq.'m1qn3') then  ! Newton minimization
c            var(i)=var(i)+deltv/scl(k)
c         else 
             var(i)=var(i)+deltv
c         endif
      endif
      enddo
      if(rcom) call fitsug(0)
      call microb
      if(parm.eq.'Flex') then
      call energy
      else if(parm.eq.'Amber91') then
      call energ91
      else if(parm.eq.'Amber94') then
      call energ94
      endif
      grd=(ener-eref)/deltv
      write(6,31) k,grc(k),grd,grc(k)-grd
31    format(2x,'A= ',i3,' GRA= ',3f10.4)
      if(rcom) call fitsug(1)
      call putbac(1)
      enddo
c---------------------------------------------------------test helical gradients
110   if(index(test,'H').eq.0) goto 200
      write(6,20)
      do k=nbcc+1,nhlc
      do i=1,nvar
      l=abs(neq(i))
      if(l.eq.k) then
c         if(minz.eq.'m1qn3') then  ! Newton mimimization
c           var(i)=var(i)+deltv/scl(k)
c         else 
            var(i)=var(i)+deltv
c         endif
      endif
      enddo
      if(rcom) call fitsug(0)
      call microb
      if(parm.eq.'Flex') then
        call energy
      else if(parm.eq.'Amber91') then
        call energ91
      else if(parm.eq.'Amber94') then
        call energ94
      endif
      grd=(ener-eref)/deltv
      write(6,34) k,grc(k),grd,grc(k)-grd
34    format(2x,'K= ',i3,' GRA= ',3f10.4)
      if(rcom) call fitsug(1)
      call putbac(1)
      enddo
c----------------------------------------------------test inter ligand gradients
200   if(index(test,'I').eq.0) goto 210
      if(nklc.lt.nlgc) write(6,20)
      do k=nklc+1,nlgc
      do i=1,nvar
      l=abs(neq(i))
      if(l.eq.k) then
c         if(minz.eq.'m1qn3') then  ! Newton minimization
c            var(i)=var(i)+deltv/scl(k)
c         else 
            var(i)=var(i)+deltv
c         endif
      endif
      enddo
      if(rcom) call fitsug(0)
      call microb
      if(parm.eq.'Flex') then
      call energy
      else if(parm.eq.'Amber91') then
      call energ91
      else if(parm.eq.'Amber94') then
      call energ94
      endif
      grd=(ener-eref)/deltv
      write(6,35) k,grc(k),grd,grc(k)-grd
35    format(2x,'R= ',i3,' GRA= ',3f10.4)
      if(rcom) call fitsug(1)
      call putbac(1)
      enddo
210   return
      end
