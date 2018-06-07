      subroutine disth
      include 'jumna_data.inc'
      logical*2 lar,lock,kink,lthy,ifhb,ribose,cation,
     1 klog(n2,50),dlog(n2,50)
      character*4 mnam,munit,seq*120,code*8,kode*8,knam,symb(6)*1,
     1 dirv(-1:1)*1,rb*1,line*80,lind(80)*1
      dimension dnoe(n2,50),knoe(n2,50,2),rms_row(n2),rms_col(50),kr(n2)
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/flx/sap(n6),refg,refb,refh,refv,iap(n8,n4,n5),
     1 nap(n8,7,n5),kap(3,n5)
      common/ind/iend(0:4),ise(n2),ienb(n2,2),kapt(n2+1),nsph
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
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      data dirv/'^',' ','v'/,symb/'-','^','v',':','/','\\'/,lind/80*'-'/
c-------------------------------------------------------calculate list distances
      nlis=ndiv(nlin)
      do is=1,nto
      ilev=ilq(is,1)
      istr=ilq(is,2)
      idir=sign(1,idr(istr))
      do k=1,nlis
      kt=ktyp(k)
      klog(is,k)=.false.
      knoe1=0
      knoe2=0
      do i=nuc(is-1)+1,nuc(is)
      if(mnam(i).eq.knam(k,1)) knoe1=i
      enddo
      if(knam(k,1).eq.'H68') knoe1=ih68(is)
      if(knoe1.ne.0) then
      iss=0
      if(kt.eq.1) then
      iss=is
      else if(kt.eq.2) then
      ils=ilev+idir
      if(ils.gt.0.and.ils.le.kseq) iss=ieq(ils,istr)
      else if(kt.eq.3) then
      ils=ilev-idir
      if(ils.gt.0.and.ils.le.kseq) iss=ieq(ils,istr)
      else if(kt.eq.4.and.istr.lt.nst) then
      iss=ieq(ilev,istr+1)
      else if(kt.eq.5.and.istr.lt.nst) then
      ils=ilev+idir
      if(ils.gt.0.and.ils.le.kseq) iss=ieq(ils,istr+1)
      else if(kt.eq.6.and.istr.lt.nst) then
      ils=ilev-idir
      if(ils.gt.0.and.ils.le.kseq) iss=ieq(ils,istr+1)
      endif
      if(iss.ne.0) then
      do i=nuc(iss-1)+1,nuc(iss)
      if(mnam(i).eq.knam(k,2)) knoe2=i
      enddo
      if(knam(k,2).eq.'H68') knoe2=ih68(iss)
      endif
      if(knoe2.ne.0) then
      klog(is,k)=.true.
      dnoe(is,k)=dis(knoe1,knoe2)
      knoe(is,k,1)=knoe1
      knoe(is,k,2)=knoe2
      endif
      endif
      enddo
      enddo
c----------------------------------------------------------output list distances
      if(nlis.ne.0) write(6,27)
27    format(/2x,'List NOE distances (angstroms) ...')
      do nl=1,nlin
      l1=ndiv(nl-1)+1
      l2=ndiv(nl)
      write(6,28) (knam(i,1)(2:),symb(ktyp(i)),
     1 knam(i,2)(2:),i=l1,l2),'Dir'
28    format(/4x,'Nucleotide',3x,8(:,1x,a3,:,a1,a3))
      write(6,*)
      do is=1,nto
      isa=nuc(is-1)+iofs(is)
      id=sign(1,idr(ilq(is,2)))
      rb=' '
      if(ribose(is)) rb='R'
      lim=13+(l2-l1+1)*8
      if(ise(is).lt.0.and.is.ne.1) write(6,29) (lind(j),j=1,lim)
      write(line,30) is,munit(isa),nunit(isa),rb
      l=16
      do k=l1,l2
      if(klog(is,k)) then
      write(line(l:),31) dnoe(is,k)
      else
      write(line(l:),32) '        '
      endif
      l=l+8
      enddo
      write(6,32) line(:l)//'   '//dirv(id)
29    format(7x,73a1)
30    format(2x,i3,') ',a4,i3,a1)
31    format(f8.1)
32    format(a)
      enddo
      enddo
c----------------------------------------------find noe correspondance distances
      kg=0
      if(nnoe.ne.0) then
      do is=1,nto
      rms_row(is)=0.
      kr(is)=0
      enddo
      do k=1,nlis
      m=0
      rms_col(k)=0.
      do is=1,nto
      dlog(is,k)=.false.
      k1=knoe(is,k,1)
      k2=knoe(is,k,2)
      do i=1,nnoe
      if(jnoe(i).eq.1.or.jnoe(i).eq.3) then
      i1=inoe(i,1)
      i2=inoe(i,2)
      if((i1.eq.k1.and.i2.eq.k2).or.(i1.eq.k2.and.i2.eq.k1)) then
      del=dnoe(is,k)-rnoe(i)
      del2=del*del
      dnoe(is,k)=del
      dlog(is,k)=.true.
      kg=kg+1
      m=m+1
      rms_col(k)=rms_col(k)+del2
      kr(is)=kr(is)+1
      rms_row(is)=rms_row(is)+del2
      endif
      endif
      enddo
      enddo
      if(m.ne.0) rms_col(k)=sqrt(rms_col(k))/m
      enddo
      do is=1,nto
      if(kr(is).ne.0) rms_row(is)=sqrt(rms_row(is))/kr(is)
      enddo
      endif
c--------------------------------------------output noe correspondance distances
      if(kg.ne.0) then
      if(nlis.ne.0) write(6,24)
24    format(/2x,'List NOE errors (angstroms) ...')
      do nl=1,nlin
      l1=ndiv(nl-1)+1
      l2=ndiv(nl)
      if(nl.ne.nlin) then
      write(6,28) (knam(i,1)(2:),symb(ktyp(i)),
     1 knam(i,2)(2:),i=l1,l2),'Dir'
      else
      write(6,28) (knam(i,1)(2:),symb(ktyp(i)),
     1 knam(i,2)(2:),i=l1,l2),'Rms'
      endif
      write(6,*)
      do is=1,nto
      isa=nuc(is-1)+iofs(is)
      id=sign(1,idr(ilq(is,2)))
      rb=' '
      if(ribose(is)) rb='R'
      lim=13+(l2-l1+1)*8
      if(nl.eq.nlin) lim=lim+1
      if(ise(is).lt.0.and.is.ne.1) write(6,29) (lind(j),j=1,lim)
      write(line,30) is,munit(isa),nunit(isa),rb
      l=16
      do k=l1,l2
      if(klog(is,k)) then
      if(dlog(is,k)) then
      write(line(l:),31) dnoe(is,k)
      else
      write(line(l:),32) '     ---'
      endif
      else
      write(line(l:),32) '        '
      endif
      l=l+8
      enddo
      if(l2.ne.nlis) then
      write(6,32) line(:l)//'   '//dirv(id)
      else
      write(line(l:),33) rms_row(is)
      write(6,32) line(:l+6)
      endif
      enddo
      write(6,34) (rms_col(k),k=l1,l2)
      enddo
33    format(2x,f4.1)
34    format(7x,'Rms: ',3x,7f8.1)
      endif
c----------------------------------------------------------------other distances
500   m=0
      do k=1,nnoe
      if(jnoe(k).eq.1.or.jnoe(k).eq.3) then
      i1=inoe(k,1)
      i2=inoe(k,2)
      xx=corm(i1,1)-corm(i2,1)
      yy=corm(i1,2)-corm(i2,2)
      zz=corm(i1,3)-corm(i2,3)
      r=sqrt(xx*xx+yy*yy+zz*zz)
      if(jnoe(k).eq.1) then
      delt=r-rnoe(k)
      m=m+1
      if(m.eq.1) write(6,52)
52    format(/2x,'Distance constraints ...',/)
      write(6,54) r,rnoe(k),delt,munit(i1),nunit(i1),mnam(i1),
     1 munit(i2),nunit(i2),mnam(i2)
54    format(2x,'D) Val: ',f7.2,' Req: ',6x,f7.2,6x,' Del: ',f7.2,
     1 2x,a4,i3,1x,a4,'-',a4,i3,1x,a4)
      else if(jnoe(k).eq.3) then
      m=m+1
      if(m.eq.1) write(6,52)
      if(r.lt.rnoe(k)) then
      delt=r-rnoe(k)
      else if(r.gt.bnoe(k)) then
      delt=r-bnoe(k)
      else
      delt=0.
      endif
      write(6,56) r,rnoe(k),bnoe(k),delt,munit(i1),nunit(i1),
     1 mnam(i1),munit(i2),nunit(i2),mnam(i2)
56    format(2x,'d) Val: ',f7.2,' Req: ',f7.2,' <-> ',f7.2,
     1 ' Del: ',f7.2,2x,a4,i3,1x,a4,'-',a4,i3,1x,a4)
      endif
      endif
      enddo
      return
      end
