      subroutine kline(l)
      include 'jumna_data.inc'
      logical*2 kink,lthy,hst,bst,vst,lgi,lgj,lgu,ribose,cation
      character seq*120,code*8,kode*8,line*80
      dimension vn(20)
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/slg/hst(n2,6),bst(n2,n8),vst(n2,4),lgi(n1),lgj(n1),lgu(n2)
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(120,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      read(51,5) line
5     format(a)
      kode(l)=line(:8)
c------------------------------------------------------------find number of data
      ic=0
      do k=1,8
      ik=ichar(line(k:k))
      if(ik.ge.49.and.ik.le.53) ic=ic+1
      enddo
      if(ic.gt.0) read(line(9:),*) (vn(k),k=1,ic)
c------------------------------------------------------assign junction variables
      do k=1,4
      vst(l,k)=.false.
      enddo
      m=0
      if(index(kode(l),'0').ne.0) then
      do k=1,4
      m=m+1
      vkink(l,k)=vn(m)
      if(k.eq.1.or.k.eq.4) vkink(l,k)=-vn(m)
      vst(l,k)=.true.
      enddo
      goto 200
      endif
      if(index(kode(l),'1').ne.0) then
      m=m+1
      vkink(l,1)=-vn(m)
      vst(l,1)=.true.
      endif
      if(index(kode(l),'2').ne.0) then
      m=m+1
      vkink(l,2)=vn(m)
      vst(l,2)=.true.
      endif
      if(index(kode(l),'4').ne.0) then
      m=m+1
      vkink(l,3)=vn(m)
      vst(l,3)=.true.
      endif
      if(index(kode(l),'5').ne.0) then
      m=m+1
      vkink(l,4)=-vn(m)
      vst(l,4)=.true.
      endif
200   return
      end
