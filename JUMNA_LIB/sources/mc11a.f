      subroutine mc11a(a,n,z,sig,w,ir,mk,eps)
      include 'jumna_data.inc'
      dimension a(*),z(*),w(*)
      if(n.gt.1)goto1
      a(1)=a(1)+sig *z(1)**2
      ir=1
      if(a(1).gt.0.)return
      a(1)=0.
      ir=0
      return
1     continue
      np=n+1
      if(sig.gt.0.)goto40
      if(sig.eq.0..or.ir.eq.0)return
      ti=1./sig
      ij=1
      if(mk.eq.0)goto10
      do 7 i=1,n
      if(a(ij).ne.0.)ti=ti+w(i)**2/a(ij)
7     ij=ij+np-i
      goto20
10    continue
      do 11 i=1,n
11    w(i)=z(i)
      do 15 i=1,n
      ip=i+1
      v=w(i)
      if(a(ij).gt.0.)goto12
      w(i)=0.
      ij=ij+np-i
      goto15
12    continue
      ti=ti+v**2/a(ij)
      if(i.eq.n)goto14
      do 13 j=ip,n
      ij=ij+1
13    w(j)=w(j)-v*a(ij)
14    ij=ij+1
15    continue
20    continue
      if(ir.le.0 )goto21
      if(ti.gt.0.)goto22
      if(mk-1)40,40,23
21    ti=0.
      ir=-ir-1
      goto23
22    ti=eps/sig
      if(eps.eq.0.)ir=ir-1
23    continue
      mm=1
      tim=ti
      do 30 i=1,n
      j=np-i
      ij=ij-i
      if(a(ij).ne.0.)tim=ti-w(j)**2/a(ij)
      w(j)=ti
30    ti=tim
      goto41
40    continue
      mm=0
      tim=1./sig
41    continue
      ij=1
      do 66 i=1,n
      ip=i+1
      v=z(i)
      if(a(ij).gt.0.)goto53
      if(ir.gt.0 .or.sig.lt.0..or.v.eq.0.)goto52
      ir=1-ir
      a(ij)=v**2/tim
      if(i.eq.n)return
      do 51 j=ip,n
      ij=ij+1
51    a(ij)=z(j)/v
      return
52    continue
      ti=tim
      ij=ij+np-i
      goto66
53    continue
      al=v/a(ij)
      if(mm)54,54,55
54    ti=tim+v*al
      goto56
55    ti=w(i)
56    continue
      r=ti/tim
      a(ij)=a(ij)*r
      if(r.eq.0.)goto70
      if(i.eq.n)goto70
      b=al/ti
      if(r.gt.4.)goto62
      do 61 j=ip,n
      ij=ij+1
      z(j)=z(j)-v*a(ij)
61    a(ij)=a(ij)+b*z(j)
      goto64
62    gm=tim/ti
      do 63 j=ip,n
      ij=ij+1
      y=a(ij)
      a(ij)=b*z(j)+y*gm
63    z(j)=z(j)-v*y
64    continue
      tim=ti
      ij=ij+1
66    continue
70    continue
      if(ir.lt.0)ir=-ir
      return
c-----------------------------------multiply vector z by inverse of factors in a
      entry mc11e(a,n,z,w,ir)
      if(ir.lt.n)return
      w(1)=z(1)
      if(n.gt.1)goto400
      z(1)=z(1)/a(1)
      return
400   continue
      do 402 i=2,n
      ij=i
      i1=i-1
      v=z(i)
      do 401 j=1,i1
      v=v-a(ij)*z(j)
401   ij=ij+n-j
      w(i)=v
402   z(i)=v
      z(n)=z(n)/a(ij)
      np=n+1
      do 411 nip=2,n
      i=np-nip
      ii=ij-nip
      v=z(i)/a(ii)
      ip=i+1
      ij=ii
      do 410 j=ip,n
      ii=ii+1
410   v=v-a(ii)*z(j)
411   z(i)=v
      return
      end
