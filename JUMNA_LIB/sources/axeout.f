      subroutine axeout(out)
      include 'jumna_data.inc'
      logical*2 lar,lock,kink,lthy,ribose,cation,locr,sup,rcom,homo,
     1 sum,homo2,homo3,diep,link,ecen,cyl,lcat,cent,autos,ihl,nuz,amber
      character*4 mnam,munit,knam,seq*120,code*8,kode*8,lnam
      character*132 out
      integer*4 opt
      dimension a(n2),vrib(n2),hcat(n2,3),icat(n2)
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
      kfi=index(out,' ')-1
      open(unit=8,file=out(:kfi)//'.axe',status='unknown')
      write(8,10) nst,(idr(j),j=1,4),rad,pit
10    format(5i4,2f8.3)
c---------------------------------------------------------------output variables
      do is=1,nto
c in - topology type
c ino - 5',3', or normal indicator
c koff - index of the last int coordinate in (sap) of the residue is
c ioff - last atom number of the previous residue 
c ilim - number of independent position variables, torsions in strjm_.set
        in=itr(is)
        ino=ito(is)
        koff=kapt(is)
        ioff=nuc(is-1)
        ilim=kap(1,ino)
        nuz=.true.
c if 3' end 
        if(ino.ne.3.and.ino.ne.6) then
           ilim=ilim-2
           nuz=.false.
        endif
        vrib(is)=0
        if(ribose(is)) then
           ilim=ilim-1
           vrib(is)=sap(koff+11)
           if(ino.eq.6) vrib(is)=sap(koff+7)
        endif
c---------------------------------------------------2nd strand zsh,wdg to global
c if first strand:
        if(is.le.kseq) then
           zglo=hel(is,3)
           wglo=hel(is,6)
           if(isur.ne.0.and.is.gt.1) wglo=hel(is,6)-hel(is-1,6)
           if(abs(wglo).gt.180) wglo=wglo-sign(360.d0,wglo)
c if not 3':
           if(nuz) then
              write(8,20) -hel(is,1),hel(is,2),zglo,hel(is,4),-hel(is,5),wglo
           else
              write(8,20) -hel(is,1),hel(is,2),zglo,hel(is,4),-hel(is,5),wglo,
     1        sap(koff+ilim+1),sap(koff+ilim+2)
           endif
           write(8,20) (sap(koff+l),l=1,ilim)
20         format(2x,8f9.4)
        else
c if second strand:
c ipos - index of the residue in the chain
c i1 - the number of the residue in the strand
c
           ipos=ilq(is,1)
           i1=ieq(ipos,1)
c first residue:
           if(ise(is).lt.0) then
             zdel=hel(is,3)
             wdel=hel(is,6)
             zglo=hel(i1,3)+zdel
             wglo=hel(i1,6)+wdel
             if(isur.ne.0) wglo=wdel
           else
             if(.not.ihl(is)) then
               zglo=hel(is,3)
               wglo=hel(is,6)
               zdel=zdel+zglo-hel(i1,3)
               wdel=wdel+wglo-hel(i1,6)
             else
               zglo=(hel(i1,3)+hel(is,3))-zdel
               wglo=hel(i1,6)+hel(is,6)-wdel
               ir=ipos-1
               it=ilq(is,2)
35             if(ieq(ir,it).eq.0) then
                 zglo=zglo+hel(i1-1,3)
                 wglo=wglo+hel(i1-1,6)
                 ir=ir-1
                 goto 35
               endif
               zdel=hel(is,3)
               wdel=hel(is,6)
             endif
             if(isur.ne.0) wglo=hel(is,6)-hel(is-1,6)
           endif
           if(abs(wglo).gt.180) wglo=wglo-sign(360.d0,wglo)
c      if not 3':
           if(nuz) then
              write(8,20) -hel(is,1),hel(is,2),zglo,hel(is,4),-hel(is,5),wglo
           else
              write(8,20) -hel(is,1),hel(is,2),zglo,hel(is,4),-hel(is,5),wglo,
     1              sap(koff+ilim+1),sap(koff+ilim+2)
           endif
           write(8,20) (sap(koff+l),l=1,ilim)
        endif
      enddo
c--------------------------------------------------------------------methyl data
      nt=ntba
      do is=1,nto
        a(is)=0.
        if(lthy(is)) then
          nt=nt+1
          a(is)=sap(nt)
        endif
      enddo
      if(nto.gt.0) write(8,30) (a(is),is=1,nto)
30    format(2x,7f9.4)
      if(nto.gt.0) write(8,30) (vrib(is),is=1,nto)
c--------------------------------------------------------------------------kinks
      if(kseq.gt.0) write(8,40) (-vkink(i,1),vkink(i,2),vkink(i,3),
     1 -vkink(i,4),kink(i),i=2,kseq)
40    format(2x,4f9.4,l2)
c-----------------------------------------------------------modified nucleotides
      do is=1,nto
        in=itr(is)
        ino=ito(is)
        koff=kapt(is)
        if(kap(1,in).ge.kap(1,ino)+1) write(8,22)
     1    (sap(koff+j),j=kap(1,ino)+1,kap(1,in))
22      format(2x,7f9.4)
      enddo
c-------------------------------------------------------------non-bonded ligands
      do il=1,nlig
        is=nto+il
        in=itr(is)
        koff=kapt(is)
        if(ilig(il,1).eq.0) then
          write(8,23) (rlig(il,j),j=1,6),ilig(il,1),slig(il,1),slig(il,2)
        else
          write(8,23) (rlig(il,j),j=1,6),ilig(il,1)
23        format(2x,6f9.4,i3,2f9.4)
        endif
        if(kap(1,in).gt.0) write(8,22) (sap(koff+j),j=1,kap(1,in))
      enddo
c-----------------------------------------------------------cations if lcat=true
      if(lcat) then
        call cataxe(hcat,icat)
        do il=1,ncat
          write(8,23) (hcat(il,j),j=1,3),0.,0.,0.,icat(il)
        enddo
      endif
      close(8)
      return
      end
