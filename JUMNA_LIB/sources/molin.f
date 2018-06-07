      subroutine molin(name,lr,koff,mods)
      include 'jumna_data.inc'
      logical*2 lthy,kink,ribose,cation,fault,locr,mods
      character*4 mac*32,lmo*32,lib*32,libm*32,out*32,axe*32,noe*32,
     1 nol*32,axl*32,pdb*32,ins*32,bar*32,parm*32,snam,suni,
     1 sub,seq*120,code*8,kode*8,name(n2),subnam,lnam,test*32
	character*80 libn
      integer*4 iofref(6),isrl(5,2)
      dimension rsrl(5)
      common/cha/mac,lmo,lib,libm,out,axe,axl,noe,nol,test,pdb,ins,
     1 bar,parm
	common/jmdirs/libn
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/flx/sap(n6),refg,refb,refh,refv,iap(n8,n4,n5),
     1 nap(n8,7,n5),kap(3,n5)
      common/ind/iend(0:4),ise(n2),ienb(n2,2),kapt(n2+1),nsph
      common/lgd/rlig(n9,6),slig(n9,2),ilig(n9,6),lopt(n9),lpiv(n9),
     1 ntlg,nlgi,nlig,ntl,lnam(n9),locr(n2,n0)
      common/moljm/sor(n4,n5,3),smon(n4,n5),snam(n4,n5),suni(n4,n5),
     1 nuni(n4,n5),sub(n5),isch(n4,n5),isty(n4,n5),ics(n4,n5),
     1 mats(3*n4,n5),kas(n5),khs(n5),ksub
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
	common/hacon/iret
      data iofref/20,19,16,21,20,17/
c-------------------------------------------------------------------------------
      if(nlig.ne.0) then
        mods=.true.
        koff=koff+1
        lrh=lr+1
        itr(lrh)=koff
        name(lrh)=lnam(1)
        do il=2,nlig
          is=nto+il
          subnam=lnam(il)
          do js=lrh,is-1
            if(name(js).eq.subnam) then
               itr(is)=itr(js)
               name(is)=subnam
               goto 101
            endif
          enddo
          koff=koff+1
          itr(is)=koff
          name(is)=subnam
101     enddo
      endif
c-------------------------------------------------------------------------------
      kfi=index(libn,' ')-1
      open(unit=3,file=libn(:kfi)//'.lre',status='old',err=600)
      k=0
      m=0
5     m=m+1
      read(3,10,end=100) subnam
      do l=1,nto
        if(subnam.eq.name(l)) then
           k=k+1
           sub(k)=subnam
           read(3,12) kas(k),khs(k)
           if(kas(k).gt.n4) then
             write(6,*) '  ---- n4 too small ----'
             return
           endif
           read(3,14) (mats(i,k),i=1,khs(k)) 
           read(3,16) (snam(i,k),sor(i,k,1),sor(i,k,2),sor(i,k,3),isch(i,k),
     1     isty(i,k),smon(i,k),suni(i,k),nuni(i,k),ics(i,k),i=1,kas(k))
10         format(a)
12         format(11i4)
14         format(15i5)
16         format(a4,3f10.5,2i2,f8.4,1x,a4,1x,i3,1x,i4)
           do j=1,nto
              if(name(l).eq.name(j)) then
                 irec(j)=k
                 iofs(j)=iofref(itr(j))
                 iofe(j)=kas(k)
                 nsr(j)=1
                 isr(1,1,j)=4
                 isr(1,2,j)=5
                 rsr(1,j)=refg
              endif
           enddo
           goto 5
         endif
      enddo
      read(3,12) kast,khst
      kl=khst/15
      if(kl*15.lt.khst) kl=kl+1
      do l=1,kl+kast
      read(3,*)
      enddo
      goto 5
100   close(3)
      if(.not.mods) goto 500
c-------------------------------------------------------------------------------
c  read modified nucleotide data base file
c 
      kfi=index(libm,' ')-1
      open(unit=4,file=libm(:kfi)//'.lmo',status='old')
      m=0
50    m=m+1
      read(4,10,end=200) subnam
      if(subnam.eq.' ') goto 200
      do l=1,ntl
        if(subnam.eq.name(l)) then
          k=k+1
          if(k.gt.n5) then
             write(6,*) '  ---- n5 too small ----'
             return
          endif
          kt=itr(l)
          sub(k)=subnam
          read(4,12) kas(k),khs(k),iofsl,iofel,itol,kap(1,kt),
     1         kap(2,kt),kap(3,kt),nsrl,ipvl,ifi4
          if(kas(k).gt.n4) then
               write(6,*) '  ---- n4 too small ----'
               return
           endif
           if(kap(2,kt).gt.n8) then
               write(6,*) '  ---- n8 too small ----'
               return
           endif
           if(nsrl.gt.n0) then
               write(6,*) '  ---- n0 too small ----'
               return
           endif
           read(4,14) (mats(i,k),i=1,khs(k))
           read(4,16) (snam(i,k),sor(i,k,1),sor(i,k,2),sor(i,k,3),isch(i,k),
     1         isty(i,k),smon(i,k),suni(i,k),nuni(i,k),ics(i,k),i=1,kas(k))
           if(nsrl.gt.0) read(4,19) (isrl(i,1),isrl(i,2),rsrl(i),i=1,nsrl)
19         format(5(2i3,f6.3))
           if(kap(2,kt).gt.0) then
               if(ifi4.eq.0) then
                  do i=1,7
                     read(4,18) (nap(j,i,kt),j=1,kap(2,kt))
18                   format(25i3)
                  enddo
               else
                  do i=1,7
                    read(4,118) (nap(j,i,kt),j=1,kap(2,kt))
118                 format(20i4)
                  enddo
               endif
           endif
           if(kap(3,kt).gt.0) then
             do j=1,kap(3,kt)
               read(4,18) (iap(j,i,kt),i=1,kas(k))
             enddo
           endif
           do i=1,kas(k)
             do j=1,4
               ic=ichar(suni(i,k)(j:j))
               if(ic.ge.65.and.ic.le.90) suni(i,k)(j:j)=char(ic+32)
             enddo
          enddo
          do j=1,ntl
            if(name(l).eq.name(j)) then
              irec(j)=k
              iofs(j)=iofsl
              iofe(j)=iofel
              ito(j)=itol
              if(l.le.nto) then
                nsr(j)=nsrl+1
                isr(1,1,j)=4
                isr(1,2,j)=5
                rsr(1,j)=refg
                do i=2,nsrl+1
                  isr(i,1,j)=isrl(i-1,1)
                  isr(i,2,j)=isrl(i-1,2)
                  rsr(i,j)=rsrl(i-1)
                  locr(j,i)=.false.
                enddo
              else
                lpiv(j-nto)=ipvl
                nsr(j)=nsrl
                do i=1,nsrl
                  isr(i,1,j)=isrl(i,1)
                  isr(i,2,j)=isrl(i,2)
                  rsr(i,j)=rsrl(i)
                  locr(j,i)=.false.
                enddo
              endif
              if(itol.gt.3) ribose(j)=.true.
            endif
          enddo
          goto 50
        endif
      enddo
      read(4,12) kast,khst,iofst,iofet,itot,kapt1,kapt2,kapt3,nsrt,ipvt
     *,ifi4
      kl=khst/15
      if(kl*15.lt.khst) kl=kl+1
      if(nsrt.gt.0) kl=kl+1
      if(ifi4.eq.0) then
      kla=kapt2/25
      if(kla*25.lt.kapt2) kla=kla+1
      else
      kla=kapt2/20
      if(kla*20.lt.kapt2) kla=kla+1
      endif
      klb=kast/25
      if(klb*25.lt.kast) klb=klb+1
      do l=1,kl+kast+7*kla+kapt3*klb
      read(4,*)
      enddo
      goto 50
200   close(4)
500   ksub=k
      fault=.false.
      do lr=1,ntl
      if(irec(lr).eq.0) then
      fault=.true.
      if(lr.le.nto) then
      write(6,65) name(lr)
65    format(/2x,'Nucleotide type ',a2,
     1 ' was not found in the libraries')
      else
      write(6,66) name(lr)
66    format(/2x,'Ligand type ',a4,
     1 ' was not found in the .lmo library')
      endif
      endif
      enddo
      if(fault) then
	  iret = 1
	endif
      return
600   write(6,'("Can not find nucleotide libary file ",A)')libn(:kfi)//'.lre'
      iret = 1
	return
      end
