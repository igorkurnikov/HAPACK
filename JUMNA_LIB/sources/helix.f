      subroutine helix
      include 'jumna_data.inc'
      logical*2 lar,lock,lthy,kink,ribose,sup,rcom,homo,homo2,homo3,
     1 diep,sum,link,ecen,cyl,lcat,cent,autos,locr,cation,ihl,amber
      character*4 mnam,munit,knam,seq*120,code*8,kode*8,
     1 dirv(-1:1)*1,rb*1,ct*1,lnam,rtl*3
      integer*4 opt
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
      common/lgd/rlig(n9,6),slig(n9,2),ilig(n9,6),lopt(n9),lpiv(n9),
     1 ntlg,nlgi,nlig,ntl,lnam(n9),locr(n2,n0)
      common/mnn/var(n7),gra(n7),con(n0*n2+n3),scl(n7),ncon,nvrc,
     1 ntba,nbac,nthe,nhel,ntki,nkin,ntri,ntot,nvar,nrin,
     1 lar(n7),lock(n6a)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/nmrjm/rnoe(n3),bnoe(n3),fnoe(n3),snoe,knam(200,2),ktyp(200),
     1 icon(n2*2),inoe(n3,9),ih68(n2),jnoe(n3),nnoe,ndcs,ndiv(0:99),nlin
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(120,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      common/symjm/xsum,ysum,itwl(n2),nvs(n7),ihm(n2),ksym(3),nbrk(2),
     1 isym,isur,isup,ihl(n2)
      data dirv/'^',' ','v'/
c-------------------------------------------------------------helical parameters
      if(nto.eq.0) goto 500
      write(6,28)
28    format(/2x,'Helical param: ',
     1 11x,'Xdisp',3x,'Ydisp',3x,'Rise',4x,'Inc',5x,'Tip',5x,
     1 'Twist',2x,'Dir',2x,'R/T'/)
      do i=1,nto
      ip=ilq(i,1)
      kp=ilq(i,2)
      if(ip.gt.isym+1) goto 190
      is=nuc(i-1)+iofs(i)
      id=sign(1,idr(ilq(i,2)))
      rb=' '
      if(ribose(i)) rb='R'
      ct=' '
      if(cation(i)) ct='*'
      if(ise(i).lt.0.and.i.ne.1) write(6,29)
29    format(7x,'--------------',
     1 '-------------------------------------------------------------')
c---------------------------------------------------2nd strand zsh,wdg to global
      if(i.le.kseq) then
      zglo=hel(i,3)
      wglo=hel(i,6)
      if(isur.ne.0.and.i.gt.1) wglo=hel(i,6)-hel(i-1,6)
      if(abs(wglo).gt.180) wglo=wglo-sign(360.d0,wglo)
      write(6,30) i,code(i),munit(is),nunit(is),rb,ct,-hel(i,1),
     1 hel(i,2),zglo,hel(i,4),-hel(i,5),wglo,dirv(id),'(V)'
30    format(2x,i3,2x,a8,') ',a4,i3,a1,a1,f7.2,5f8.2,3x,a1,3x,a3)
      else
      ip=ilq(i,1)
      i1=ieq(ip,1)
      if(ise(i).lt.0) then
      zdel=hel(i,3)
      wdel=hel(i,6)
      zglo=hel(i1,3)+zdel
      wglo=hel(i1,6)+wdel
      if(isur.ne.0) wglo=wdel
      else
         if(.not.ihl(i)) then
         zglo=hel(i,3)
         wglo=hel(i,6)
         zdel=zdel+zglo-hel(i1,3)
         wdel=wdel+wglo-hel(i1,6)
         else
         zglo=(hel(i1,3)+hel(i,3))-zdel
         wglo=hel(i1,6)+hel(i,6)-wdel
            ir=ip-1
            it=ilq(i,2)
35          if(ieq(ir,it).eq.0) then
            zglo=zglo+hel(i1-1,3)
            wglo=wglo+hel(i1-1,6)
            ir=ir-1
            goto 35
            endif
         zdel=hel(i,3)
         wdel=hel(i,6)
         endif
         if(isur.ne.0) wglo=hel(i,6)-hel(i-1,6)
      endif
      if(abs(wglo).gt.180) wglo=wglo-sign(360.d0,wglo)
      rtl='(H)'
      if(ihl(i)) rtl='(V)'
      write(6,30) i,code(i),munit(is),nunit(is),rb,ct,-hel(i,1),
     1 hel(i,2),zglo,hel(i,4),-hel(i,5),wglo,dirv(id),rtl
      endif
190   enddo
c---------------------------------------------------------------------superhelix
      if(isur.ne.0) write(6,27) rad,pit
27    format(/2x,'Suprhel param:  Rad= ',f8.3,'   Pit= ',f8.3)
c--------------------------------------------------------------------------kinks
      do i=1,kseq
      if(kink(i)) goto 490
      enddo
      goto 500
490   write(6,32)
32    format(/2x,'Kink params:',15x,'Ax',7x,'Ay',7x,'Ainc',5x,'Atip',
     1 1x,'State'/)
      do i=1,kseq
      is=nuc(i-1)+iofs(i)
      write(6,34) i,kode(i),munit(is),nunit(is),-vkink(i,1),vkink(i,2),
     1 vkink(i,3),-vkink(i,4),kink(i)
      enddo
34    format(2x,i3,2x,a8,') ',a4,i3,4f9.2,2x,l2)
c-------------------------------------------------------------------ligand inter
500   if(nlig.ne.0) then
      write(6,39)
39    format(/2x,'Ligand inter : ',
     1 11x,'Xdisp',3x,'Ydisp',3x,'Rise',4x,'Inc',5x,'Tip',
     1 5x,'Twist',2x,'Base'/)
      do il=1,nlig
      is=nto+il
      isa=nuc(is-1)+1
      kal=nuc(is)-nuc(is-1)
      i1=ilig(il,1)
      if(i1.eq.0) then
      if(kal.ne.1) then
      write(6,40) il,code(is),munit(isa),nunit(isa),
     1 (rlig(il,j),j=1,6),'----',0
      else
      write(6,42) il,code(is),munit(isa),nunit(isa),
     1 (rlig(il,j),j=1,3),'----',0
      endif
      else
      isb=nuc(i1-1)+iofs(i1)
      if(kal.ne.1) then
      write(6,40) il,code(is),munit(isa),nunit(isa),
     1 (rlig(il,j),j=1,6),munit(isb),nunit(isb)
40    format(2x,i3,2x,a8,') ',a4,i3,1x,6f8.2,1x,a4,i3)
      else
      write(6,42) il,code(is),munit(isa),nunit(isa),
     1 (rlig(il,j),j=1,3),munit(isb),nunit(isb)
42    format(2x,i3,2x,a8,') ',a4,i3,1x,3f8.2,25x,a4,i3)
      endif
      endif
      enddo
      endif
      return
      end
