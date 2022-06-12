      subroutine present
      include 'jumna_data.inc'
      logical*2 ifhb,lthy,kink,ribose,cation
      character*4 mnam,munit,seq*120,code*8,kode*8,rb*1,ct*1
      dimension ik(6),doc(n2),aoc(n2,4),sug(n2),refa(4)
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/flx/sap(n6),refg,refb,refh,refv,iap(n8,n4,n5),
     1 nap(n8,7,n5),kap(3,n5)
      common/ind/iend(0:4),ise(n2),ienb(n2,2),kapt(n2+1),nsph
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(120,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      common/parjm/aij(25,25),bij(25,25),ahb(15),bhb(15),vt(35),va(25),
     1 vo(25),delmax,time0,nzsh,nwdg,ivf(35),ifhb(25,25)
      refa(1)=vo(20)
      refa(2)=vo(21)
      refa(3)=vo(20)
      refa(4)=vo(2)
c----------------------------------------------------backbone angles and closure
      do i=1,nto
      ic=nuc(i-1)+4
      io=nuc(i-1)+5
      sug(i)=dis(ic,io)-refg
      enddo
      do i=1,nto-1
      if(ise(i+1).lt.0) goto 200
      if(idr(ilq(i,2)).gt.0) then
      ik(1)=nuc(i-1)+3
      ik(2)=nuc(i-1)+7
      ik(3)=nuc(i-1)+15
      ik(4)=nuc(i-1)+16
      ik(5)=nuc(i)+6
      ik(6)=nuc(i)+4
         else
         ik(1)=nuc(i)+3
         ik(2)=nuc(i)+7
         ik(3)=nuc(i)+15
         ik(4)=nuc(i)+16
         ik(5)=nuc(i-1)+6
         ik(6)=nuc(i-1)+4
         endif
      doc(i)=dis(ik(4),ik(5))-refb
      do j=1,4
      aoc(i,j)=ang(ik(j),ik(j+1),ik(j+2))-refa(j)
      enddo
200   enddo
c-------------------------------------------------------------------------output
      if(nto.eq.0) return
      write(6,10)
10    format(/2x,'Dist/Ang',5x,'   C4-O1','   C5-O5',
     1   '  C-O3-P','   O-P-O','  P-O5-C','  O-C5-C'/)
      do i=1,nto
      in=nuc(i-1)+iofs(i)
      rb=' '
      if(ribose(i)) rb='R'
      ct=' '
      if(cation(i)) ct='*'
      if(ise(i).gt.0) then
      write(6,20) i,munit(in),nunit(in),rb,ct,sug(i)
      else
      if(ise(i).lt.0.and.i.gt.1) write(6,15)
      write(6,20) i,munit(in),nunit(in),rb,ct,sug(i),doc(i),
     1 (aoc(i,j),j=1,4)
      endif
      enddo
15    format(6x,
     1 '---------------------------------------------------------')
20    format(2x,i3,':',a4,i3,a1,a1,6f8.2)
      return
      end
