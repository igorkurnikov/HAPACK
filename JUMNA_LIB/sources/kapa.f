      subroutine kapa
      include 'jumna_data.inc'
      character*4 mnam,munit,seq*120,code*8,kode*8,lnam,atclas*2
      logical*2 lthy,lar,lock,kink,ifhb,ribose,cation,locr,ihl,found,
     1 modk,found1
      dimension modk(n7),iv1(3),iv2(3)
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
      common/parjm/aij(25,25),bij(25,25),ahb(15),bhb(15),vt(35),va(20),
     1 vo(20),delmax,time0,nzsh,nwdg,ivf(35),ifhb(25,25)
      common/para/ro(57),reps(57),angc(191),fangc(191),btor(88),
     $ gtor(88),ftor(88),bdtor(30),gdtor(30),fdtor(30),amb(57,57),
     $ bmb(57,57),amh(9,11),bmh(9,11),mtor(88),itor(88,4),iadt(30,4),
     $ iangc(191,3),mof,nof,nang,ntor,ntors,nad,atclas(57)
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      common/symjm/xsum,ysum,itwl(n2),nvs(n7),ihm(n2),ksym(3),nbrk(2),
     1 isym,isur,isup,ihl(n2)
      data iv1/3,5,12/,iv2/13,14,19/
      k=0
      ki=0
      do i=1,n5
        modk(i)=.false. ! flag to identify that a particular topology type was processed
      enddo
      do is=1,nto   ! cycle on residues
        in=itr(is)  ! topology type
c      write(6,*)' kap ',(kap(i,in),i=1,3),modk(in)
        if(.not.modk(in)) then  
          do i=1,kap(2,in)  ! cycle over all variables of the residue(topology) type
            if(nap(i,6,in).lt.0) then
              ii=i
              goto 1200
            endif
          enddo
1200      nkap=ii-1
          idep=nkap  ! the number of all internal vars??
c      write(6,*)' nkap= ',nkap
          ioff=nuc(is-1)
          idir=0
          if(is.le.nto) idir=idr(ilq(is,2)) ! direction of the strand for a given residue
          krin=kap(1,in)+nsr(is)*5   ! nsr(is) - the number of dep variables
          do l=1,idep
            ia=nap(l,1,in)+ioff
            ib=nap(l,2,in)+ioff
            ic=nap(l,3,in)+ioff
            id=nap(l,4,in)+ioff
            if(id.eq.ioff) then  ! nap(l,4,in) == 0  - Valence angle
c---------------------------------------------------------------- angle
c---- angle itself
               icase=0
               ix=0
c      if(l.gt.krin) then
               ix=nap(l,5,in)
               if(idir.gt.0) then  ! set atom numbers to the next or previous residue
                 if(ix.ge.1) ic=ic-ioff+nuc(is)
                 if(ix.eq.2) ib=ib-ioff+nuc(is)
                 if(ix.le.-1) ic=ic-ioff+nuc(is-2)
                 if(ix.eq.-2) ib=ib-ioff+nuc(is-2)
                 icase=1
               else
                 if(ix.ge.1) ic=ic-ioff+nuc(is-2)
                 if(ix.eq.2) ib=ib-ioff+nuc(is-2)
                 if(ix.le.-1) ic=ic-ioff+nuc(is)
                 if(ix.eq.-2) ib=ib-ioff+nuc(is)
                 icase=2
               endif
c      endif
               ima=imty(ia) ! force field atom type index
               imb=imty(ib)
               imc=imty(ic)
c  set type of valence angle force field type
               do j=1,nang   ! cycle over valence angle types 
                 if(imb.eq.iangc(j,2)) then
                   if((ima.eq.iangc(j,1).and.imc.eq.iangc(j,3)).or.
     $                (imc.eq.iangc(j,1).and.ima.eq.iangc(j,3))) then
                      nap(l,7,in)=j
                      goto 99
                   endif
                 endif
               enddo
99             found=.false.
               do i=1,matd(ib,7)  ! cycle over atoms bonded to the middle atom of the valence angle
                 ii=abs(matd(ib,i))
                 if(ii.ne.ia.and.ii.ne.ic) then  ! if the bonded atom not in the valence angle 
                   imi=imty(ii) 
                   if(.not.found) then
                     ll=ii
                     iml=imty(ii)
                   endif
                   nkap=nkap+1
                   if(nkap.gt.n8) then
                      write(6,*)' Increase n8 '
                      stop
                   endif
                   nap(nkap,1,in)=nap(l,1,in)  ! add valence angle to nap array
                   nap(nkap,2,in)=nap(l,2,in)
                   nap(nkap,3,in)=ii-ioff
                   nap(nkap,4,in)=0
                   nap(nkap,5,in)=0
                   nap(nkap,6,in)=-l
c nie ma tego
c         if(icase.eq.1) then
c      if(ix.ge.1) nap(nkap,3,in)=ii-nuc(is)
c      if(ix.le.-1) nap(nkap,3,in)=ii-nuc(is-2)
c      else if(icase.eq.2) then
c      if(ix.ge.1) nap(nkap,3,in)=ii-nuc(is-2)
c      if(ix.le.-1) nap(nkap,3,in)=ii-nuc(is)
c      endif
                   nap(nkap,5,in)=0
                   do j=1,nang
                     if(imb.eq.iangc(j,2)) then
                       if((ima.eq.iangc(j,1).and.imi.eq.iangc(j,3)).or.
     $                    (imi.eq.iangc(j,1).and.ima.eq.iangc(j,3))) then
                          nap(nkap,7,in)=j
                          goto 100
                       endif
                     endif
                   enddo
                   write(6,*)' not found angle ',ima,imb,imi,' ',mnam(ia),
     $                         mnam(ib),mnam(ii)
c      if(ima.eq.39.or.imi.eq.39)nap(nkap,7,in)=79
c      if(ima.eq.43.or.imi.eq.43)nap(nkap,7,in)=81
100                nkap=nkap+1
                   if(nkap.gt.n8) then
                     write(6,*)' Increase n8 '
                     stop
                   endif
                   nap(nkap,1,in)=ii-ioff
                   nap(nkap,2,in)=nap(l,2,in)
                   nap(nkap,3,in)=nap(l,3,in)
                   nap(nkap,4,in)=0
                   nap(nkap,5,in)=0
                   nap(nkap,6,in)=-l
                   nap(nkap,5,in)=ix
                   do j=1,nang
                     if(imb.eq.iangc(j,2)) then
                        if((imc.eq.iangc(j,1).and.imi.eq.iangc(j,3)).or.
     $                    (imi.eq.iangc(j,1).and.imc.eq.iangc(j,3))) then
                          nap(nkap,7,in)=j
                          goto 101
                        endif
                     endif
                   enddo
                   write(6,*)' not found angle ',imi,imb,imc,' ',mnam(ii),
     $                        mnam(ib),mnam(ic)
c      if(imc.eq.39.or.imi.eq.39)nap(nkap,7,in)=79
c      if(imc.eq.43.or.imi.eq.43)nap(nkap,7,in)=81
101                if(found) then
                     nkap=nkap+1
                     if(nkap.gt.n8) then
                        write(6,*)' Increase n8 '
                        stop
                     endif
                     nap(nkap,1,in)=ii-ioff
                     nap(nkap,2,in)=nap(l,2,in)
                     nap(nkap,3,in)=ll-ioff
                     nap(nkap,4,in)=0
                     nap(nkap,5,in)=0
                     nap(nkap,6,in)=-l
                     nap(nkap,5,in)=0
                     do j=1,nang
                       if(imb.eq.iangc(j,2)) then
                         if((imi.eq.iangc(j,1).and.iml.eq.iangc(j,3)).or.
     $                      (iml.eq.iangc(j,1).and.imi.eq.iangc(j,3))) then
                             nap(nkap,7,in)=j
                             goto 102
                         endif
                       endif
                     enddo
                     write(6,*)' not found angle ',imi,imb,iml,' ',mnam(ii),
     $                         mnam(ib),mnam(ll)
c      if(imi.eq.39.or.iml.eq.39)nap(nkap,7,in)=79
c      if(imi.eq.43.or.iml.eq.43)nap(nkap,7,in)=81
                   endif
102                found=.true.
                 endif
               enddo
            else     ! nap(l,4,in) != 0 - torsion angle
c--------------------------------------------------------------- torsion
              icase=0
              ix=0
              if(l.gt.krin) then
                ix=nap(l,5,in)
                if(idir.gt.0) then
                  if(ix.ge.1) id=id-ioff+nuc(is)
                  if(ix.eq.2) ic=ic-ioff+nuc(is)
                  if(ix.le.-1) id=id-ioff+nuc(is-2)
                  if(ix.eq.-2) ic=ic-ioff+nuc(is-2)
                  icase=1
                else
                  if(ix.ge.1) id=id-ioff+nuc(is-2)
                  if(ix.eq.2) ic=ic-ioff+nuc(is-2)
                  if(ix.le.-1) id=id-ioff+nuc(is)
                  if(ix.eq.-2) ic=ic-ioff+nuc(is)
                  icase=2
                endif
              endif
              ima=imty(ia)
              imb=imty(ib)
              imc=imty(ic)
              imd=imty(id)
              found=.false.
c--------------------------------- new stuff
              do i=1,matd(ib,7)
                ii=abs(matd(ib,i))
                if(ii.ne.ic) then
                  imi=imty(ii)
                  do j=1,matd(ic,7)
                    jj=abs(matd(ic,j))
                    if(jj.ne.ib) then
                      imj=imty(jj)
                      found=.false.
                      found1=.false.
                      do k=1,ntors
                        if((imb.eq.itor(k,2).and.imc.eq.itor(k,3)).or.
     $                     (imb.eq.itor(k,3).and.imc.eq.itor(k,2))) then
                           if((imi.eq.itor(k,1).and.imj.eq.itor(k,4)).or.
     $                        (imj.eq.itor(k,1).and.imi.eq.itor(k,4))) then
                              if(ii.eq.ia.and.jj.eq.id.and..not.found.or.
     $                           ii.eq.id.and.jj.eq.ia.and..not.found ) then
                                 nap(l,7,in)=k
                                 found=.true.
                              else
                                 nkap=nkap+1
                                 if(nkap.gt.n8) then
                                   write(6,*)' Increase n8 '
                                   stop
                                 endif
                                 nap(nkap,1,in)=ii-ioff
                                 nap(nkap,2,in)=nap(l,2,in)
                                 nap(nkap,3,in)=nap(l,3,in)
                                 nap(nkap,4,in)=jj-ioff
                                 if(icase.eq.1) then
                                   if(ix.eq.2) nap(nkap,4,in)=jj-nuc(is)
                                   if(ix.eq.1.and.jj.eq.id)nap(nkap,4,in)=jj-nuc(is)
                                   if(ix.eq.-2) nap(nkap,4,in)=jj-nuc(is-2)
                                   if(ix.eq.-1.and.jj.eq.id)nap(nkap,4,in)=jj-nuc(is-2)
                                 else if(icase.eq.2) then
                                   if(ix.eq.2) nap(nkap,4,in)=jj-nuc(is-2)
                                   if(ix.eq.1.and.jj.eq.id)nap(nkap,4,in)=jj-nuc(is-2)
                                   if(ix.eq.-2) nap(nkap,4,in)=jj-nuc(is)
                                   if(ix.eq.-1.and.jj.eq.id)nap(nkap,4,in)=jj-nuc(is)
                                 endif
                                 nap(nkap,5,in)=ix
                                 if(abs(ix).le.1.and.nunit(jj).eq.nunit(ic))nap(nkap,5,in)=0
                                 if(found) then
                                    nap(nkap,6,in)=-l
                                 else if(found1) then
                                    nap(nkap,6,in)=-l
                                 else
                                    ll=nkap
                                    found1=.true.
                                    nap(nkap,6,in)=-l
                                 endif
                                 nap(nkap,7,in)=k
                                 found=.true.
                              endif
                           endif
                        endif
                      enddo
                      if(.not.found) then
                        do k=ntors+1,ntor
                          if((imb.eq.itor(k,2).and.imc.eq.itor(k,3)).or.
     $                     (imb.eq.itor(k,3).and.imc.eq.itor(k,2))) then
                           if(ii.eq.ia.and.jj.eq.id.or.ii.eq.id.and.jj.eq.ia) then
                              nap(l,7,in)=k
                              found=.true.
                           else
                              nkap=nkap+1
                              if(nkap.gt.n8) then
                                write(6,*)'  Increase n8 '
                                stop
                              endif
                              nap(nkap,1,in)=ii-ioff
                              nap(nkap,2,in)=nap(l,2,in)
                              nap(nkap,3,in)=nap(l,3,in)
                              nap(nkap,4,in)=jj-ioff
                              if(icase.eq.1) then
                                if(ix.eq.2) nap(nkap,4,in)=jj-nuc(is)
                                if(ix.eq.1.and.jj.eq.id)nap(nkap,4,in)=jj-nuc(is)
                                if(ix.eq.-2) nap(nkap,4,in)=jj-nuc(is-2)
                                if(ix.eq.-1.and.jj.eq.id)nap(nkap,4,in)=jj-nuc(is-2)
                              else if(icase.eq.2) then
                                if(ix.eq.2) nap(nkap,4,in)=jj-nuc(is-2)
                                if(ix.eq.1.and.jj.eq.id)nap(nkap,4,in)=jj-nuc(is-2)
                                if(ix.eq.-2) nap(nkap,4,in)=jj-nuc(is)
                                if(ix.eq.-1.and.jj.eq.id)nap(nkap,4,in)=jj-nuc(is)
                              endif
                              nap(nkap,5,in)=ix
                              if(abs(ix).le.1.and.nunit(jj).eq.nunit(ic))nap(nkap,5,in)=0
                              nap(nkap,6,in)=-l
                              nap(nkap,7,in)=k
                              found=.true.
                           endif
                          endif
                        enddo ! ok general
                      endif ! if not found specific
                      if(.not.found)write(6,*)' not found angle ',imi,imb,imc,imj,
     $                     ' ',mnam(ii),mnam(ib),mnam(ic),mnam(jj),l,nap(l,1,in),
     $                     nap(l,2,in),nap(l,3,in),nap(l,4,in)
                    endif
                  enddo ! jj
                endif
              enddo !ii
            endif ! if torsion
          enddo ! of l
          if(in.eq.1) then
            do i1=1,3
              do i2=1,3
                nkap=nkap+1
                nap(nkap,1,in)=iv1(i1)
                nap(nkap,2,in)=4
                nap(nkap,3,in)=6
                nap(nkap,4,in)=iv2(i2)
                nap(nkap,5,in)=0
                nap(nkap,6,in)=0
                ia=iv1(i1)+ioff
                ib=4+ioff
                ic=6+ioff
                id=iv2(i2)+ioff
                ima=imty(ia)
                imb=imty(ib)
                imc=imty(ic)
                imd=imty(id)
                found=.false.
                found1=.false.
                do k=1,ntors
                  if((imb.eq.itor(k,2).and.imc.eq.itor(k,3)).or.
     $              (imb.eq.itor(k,3).and.imc.eq.itor(k,2))) then
                     if((ima.eq.itor(k,1).and.imd.eq.itor(k,4)).or.
     $                 (imd.eq.itor(k,1).and.ima.eq.itor(k,4))) then
                       if(.not.found) then
                         nap(nkap,7,in)=k
                         found=.true.
                       else
                         nkap=nkap+1
                         if(nkap.gt.n8) then
                            write(6,*)' Increase n8 '
                            stop
                         endif
                         nap(nkap,1,in)=ia-ioff
                         nap(nkap,2,in)=ib-ioff
                         nap(nkap,3,in)=ic-ioff
                         nap(nkap,4,in)=id-ioff
                         nap(nkap,7,in)=k
                         found=.true.
                       endif
                     endif
                  endif
                enddo
                if(.not.found) then
                  do k=ntors+1,ntor
                     if((imb.eq.itor(k,2).and.imc.eq.itor(k,3)).or.
     $                  (imb.eq.itor(k,3).and.imc.eq.itor(k,2))) then
                         nap(nkap,7,in)=k
                         found=.true.
                     endif
                  enddo
                endif
                if(.not.found)write(6,*)' not found angle ',ima,imb,imc,imd,
     $             ' ',mnam(ia),mnam(ib),mnam(ic),mnam(id)
              enddo
            enddo
          endif ! if x unit
          kap(2,in)=nkap
          modk(in)=.true.
c      write(6,*)' is= ',is,' in= ',in,' nkap= ',nkap
c      do kk=1,nkap
c      write(6,*)kk,') ',(nap(kk,jj,in),jj=1,7)
c      enddo
        endif
      enddo
      return
      end
 
