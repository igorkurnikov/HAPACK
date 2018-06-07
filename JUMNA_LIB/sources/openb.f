      subroutine openb(nop,rpiv,tpiv)
      include 'jumna_data.inc'
      logical*2 lthy,kink,ribose,cation
      character*4 mnam,munit,seq*120,code*8,kode*8,tran(15)*1,perp(6),
     1 tram(15)
      dimension cp(3,3)
      common/dob/theta(20),ind(20,3),nba(20)
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/hel/ua(n2,3),da(n2,3),ra(n2,3)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      data perp/'C8','N9','C4','C6','N1','C2'/
      data tran/'G','A','I','C','T','U',
     1          'H','B','J','D','S','V','O','Q',
     1                      'R'/,
     1     tram/'g','a','i','c','t','u',
     1          'h','b','j','d','s','v','o','q',
     1                      'r'/
      do ib=1,nop
      is=nba(ib)
      ks=(ilq(is,2)-1)*kseq+ilq(is,1)
      do l=1,15
      if(seq(ks:ks).eq.tran(l).or.seq(ks:ks).eq.tram(l)) goto 200
      enddo
200   if(l.eq.13) then
      write(6,50) is
50    format(/2x,'---- No opening allowed for unit ',i2,' ----'/)
      stop
      endif
      loff=0
      if((l.ge.4.and.l.le.6).or.l.ge.10) loff=3
      ioff=iofs(is)
c-------------------------------------------------------pivot points and vectors
      do i=nuc(is-1)+ioff,nuc(is)
      do l=1,3
      if(mnam(i).eq.perp(l+loff)) then
      cp(l,1)=corm(i,1)
      cp(l,2)=corm(i,2)
      cp(l,3)=corm(i,3)
      if(l.eq.2) ind(ib,2)=i
      if(l.eq.3) ind(ib,3)=i
      endif
      enddo
      enddo
      ax=cp(1,1)-cp(2,1)
      ay=cp(1,2)-cp(2,2)
      az=cp(1,3)-cp(2,3)
      bx=cp(3,1)-cp(2,1)
      by=cp(3,2)-cp(2,2)
      bz=cp(3,3)-cp(2,3)
      rx=by*az-bz*ay
      ry=bz*ax-bx*az
      rz=bx*ay-by*ax
      rp=sqrt(rx*rx+ry*ry+rz*rz)
      rx=rx/rp
      ry=ry/rp
      rz=rz/rp
      ca=cos(cdr*(tpiv))
      sa=sin(cdr*(tpiv))
      x0=corm(nuc(is-1)+1,1)
      y0=corm(nuc(is-1)+1,2)
      z0=corm(nuc(is-1)+1,3)
      cx=x0-cp(2,1)
      cy=y0-cp(2,2)
      cz=z0-cp(2,3)
      fac=rpiv/sqrt(cx*cx+cy*cy+cz*cz)
      xx=cx*fac
      yy=cy*fac
      zz=cz*fac
      ptx=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     1   (rx*rz*(1-ca)+ry*sa)*zz+x0
      pty=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1   (ry*rz*(1-ca)-rx*sa)*zz+y0
      ptz=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1   (rz*rz+(1-rz*rz)*ca)*zz+z0
c---------------------------------------------------------------distort geometry
      ca=cos(cdr*(theta(ib)))
      sa=sin(cdr*(theta(ib)))
      do i=nuc(is-1)+1,nuc(is)
      xx=corm(i,1)-ptx
      yy=corm(i,2)-pty
      zz=corm(i,3)-ptz
      corm(i,1)=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     1         (rx*rz*(1-ca)+ry*sa)*zz+ptx
      corm(i,2)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1         (ry*rz*(1-ca)-rx*sa)*zz+pty
      corm(i,3)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1         (rz*rz+(1-rz*rz)*ca)*zz+ptz
      enddo
      enddo
      call helloc(nop)
      return
      end
