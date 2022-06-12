      subroutine liner(l,line)
c
c  Set nucleotide coordinates from the line
c  l - nucleotide residue number
c
      include 'jumna_data.inc'
      logical*2 kink,lthy,hst,bst,vst,lgi,lgj,lgu,ribose,cation,
     1 sup,rcom,homo,homo2,homo3,diep,link,ecen,cyl,lcat,cent,
     1 sum,autos,amber
      character seq*120,code*8,kode*8,line*100
      integer*4 opt
      dimension vn(20)
      common/datjm/acc,phos,epsi,epsr,plat,slope,rhbl,vfac,tfac,rfac,xfac,
     1 scale,rpiv,tpiv,fad,fan,damp,fnoew,fnoes,fnoea,df1,scnb,scee,
     1 catd,catr,catc,rad,pit,enit,nshel,maxn,opt,mhomo,mhomo2,mhomo3,
     1 nick,limit,nop,nrib,ncat,lig,naxo,sup,rcom,homo,homo2,homo3,diep,
     1 sum,link,ecen,cyl,lcat,cent,autos,amber
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/flx/sap(n6),refg,refb,refh,refv,iap(n8,n4,n5),
     1 nap(n8,7,n5),kap(3,n5)
      common/ind/iend(0:4),ise(n2),ienb(n2,2),kapt(n2+1),nsph
      common/slg/hst(n2,6),bst(n2,n8),vst(n2,4),lgi(n1),lgj(n1),lgu(n2)
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(120,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      in=itr(l)
      ino=ito(l)
      do k=1,6
      hst(l,k)=.false.
      enddo
      do k=1,kap(1,in)
      bst(l,k)=.false.
      enddo
      code(l)=line(:8)
c------------------------------------------------------------find number of data
      ic=0
c indicate to set all helicoidal parameters
      if(index(line(:8),'0').ne.0) ic=6
c indicate to set all backbone parameters
      if(index(line(:8),'#').ne.0) ic=ic+5
      do k=1,8
      ik=ichar(line(k:k))
c if line(k:k) 1 to 7 - helicoidal params or glycosydic bond
      if(ik.ge.49.and.ik.le.55) ic=ic+1
c if line(k:k) 8 to 9 - sugar torsion or backbone torsions
      if(ik.ge.56.and.ik.le.57) ic=ic+2
      enddo
      if(ribose(l)) then
      if(index(line(:8),'#').ne.0) ic=ic+1
      if(index(line(:8),'8').ne.0) ic=ic+1
      endif
c read parameters from the line
      read(line(9:),*) (vn(k),k=1,ic)
c-------------------------------------------------------assign helical variables
      m=0
c if 
      if(index(code(l),'0').ne.0) then
      do k=1,6
      m=m+1
      hel(l,k)=vn(m)
      if(k.eq.1.or.k.eq.5) hel(l,k)=-vn(m)
      hst(l,k)=.true.
      enddo
      goto 200
      endif
      if(index(code(l),'1').ne.0) then
      m=m+1
      hel(l,1)=-vn(m)
      hst(l,1)=.true.
      endif
      if(index(code(l),'2').ne.0) then
      m=m+1
      hel(l,2)=vn(m)
      hst(l,2)=.true.
      endif
      if(index(code(l),'3').ne.0) then
      m=m+1
      hel(l,3)=vn(m)
      hst(l,3)=.true.
      endif
      if(index(code(l),'4').ne.0) then
      m=m+1
      hel(l,4)=vn(m)
      hst(l,4)=.true.
      endif
      if(index(code(l),'5').ne.0) then
      m=m+1
      hel(l,5)=-vn(m)
      hst(l,5)=.true.
      endif
      if(index(code(l),'6').ne.0) then
      m=m+1
      hel(l,6)=vn(m)
      hst(l,6)=.true.
      endif
c------------------------------------------------------assign backbone variables
200   if(index(code(l),'#').ne.0) then
        if(abs(ise(l)).ne.3) then
          do k=4,8
            m=m+1
            set(l,k)=vn(m)
            bst(l,k)=.true.
          enddo
          if(ribose(l)) then
            m=m+1
            set(l,11)=vn(m)
            bst(l,11)=.true.
          endif
        else
          do k=4,6
            m=m+1
            set(l,k)=vn(m)
            bst(l,k)=.true.
          enddo
          m=m+2
          if(ribose(l)) then
            m=m+1
            set(l,7)=vn(m)
            bst(l,7)=.true.
          endif
        endif
      goto 300
      endif
      if(index(code(l),'7').ne.0) then
      m=m+1
      set(l,4)=vn(m)
      bst(l,4)=.true.
      endif
      if(index(code(l),'8').ne.0) then
      m=m+1
      set(l,5)=vn(m)
      bst(l,5)=.true.
      m=m+1
      set(l,6)=vn(m)
      bst(l,6)=.true.
      endif
      if(index(code(l),'9').ne.0.and.abs(ise(l)).ne.3) then
      m=m+1
      set(l,7)=vn(m)
      bst(l,7)=.true.
      m=m+1
      set(l,8)=vn(m)
      bst(l,8)=.true.
      endif
      if(index(code(l),'8').ne.0.and.ribose(l)) then
        m=m+1
        if(abs(ise(l)).ne.3) then
          set(l,11)=vn(m)
          bst(l,11)=.true.
        else
          set(l,7)=vn(m)
          bst(l,7)=.true.
        endif
      endif
 300  return
      end
