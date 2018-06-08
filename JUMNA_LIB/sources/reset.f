      subroutine reset(inaxe)
      include 'jumna_data.inc'
      character*4 mnam,munit,knam,seq*120,code*8,kode*8,lnam
      character*4 inaxe*32,line*80
      logical*2 lthy,lar,lock,kink,ribose,cation,locr,sup,rcom,homo,
     1 homo2,homo3,diep,link,ecen,cyl,lcat,cent,autos,ihl,amber,sum,new
      integer*4 opt
      dimension sett(n2),het(n2,6),rlit(n9,6),vkint(n2,4),tht(n2)
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
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      common/symjm/xsum,ysum,itwl(n2),nvs(n7),ihm(n2),ksym(3),nbrk(2),
     1 isym,isur,isup,ihl(n2)
      data refo,refp/119.13,101.39/
c======================================================================axe input
      kfi=index(inaxe,' ')-1
      open(unit=7,file=inaxe(:kfi)//'.axe',status='old')
      read(7,51) line
51    format(a)
      new=.true.
      if(index(line,'#').eq.0) then
        read(line,1) rado,pito
1       format(20x,2f8.3)
      else
         rado=0.
         pito=0.
      endif
      do is=1,nto
        read(7,50) (het(is,j),j=1,6),set(is,9),set(is,10)
        read(7,50) (set(is,j),j=1,8)
        if(set(is,9).le.0)  set(is,9) =refo
        if(set(is,10).le.0) set(is,10)=refp
50      format(2x,8f9.4)
          het(is,1)=-het(is,1)
          het(is,5)=-het(is,5)
      enddo
      if(nto.gt.0) read(7,60) (tht(is),is=1,nto)
      if(nto.gt.0) read(7,60) (sett(is),is=1,nto)
60    format(2x,7f9.4)
      do is=1,nto
        if(abs(ise(is)).ne.3) then
          set(is,11)=sett(is)
        else
          set(is,7)=sett(is)
        endif
      enddo
      do is=2,kseq
        read(7,70) (vkint(is,j),j=1,4)
70      format(2x,4f9.4)
        vkint(is,1)=-vkint(is,1)
        vkint(is,4)=-vkint(is,4)
      enddo
c-----------------------------------------------------------modified nucleotides
      do is=1,nto
		in=itr(is)
		ino=ito(is)
		idel=kap(1,in)-kap(1,ino)
		if(idel.gt.0) read(7,60) (set(is,j),j=kap(1,ino)+1,kap(1,in))
      enddo
c-------------------------------------------------------------non-bonded ligands
      do il=1,nlig
        is=nto+il
        in=itr(is)
        if(nto.eq.0) then
          read(7,23) (rlit(il,j),j=1,6),ilis,slig(il,1),slig(il,2)
23           format(2x,6f9.4,i3,2f9.4)
        else
          read(7,23) (rlit(il,j),j=1,6),ilis
        endif
        if(kap(1,in).gt.0) read(7,60) (set(is,j),j=1,kap(1,in))
      enddo
      close(7)
c-------------------------------------------------update superhelical rise/twist
      if(isur.ne.0) then
        do ks=1,nst
          dw=0.
          do is=1,kseq
            lr=ieq(is,ks)
            if(lr.ne.0) then
              dw=dw+het(lr,6)
              if(dw.ge.360) dw=dw-360
              het(lr,6)=dw
            endif
          enddo
        enddo
      endif
c----------------------------------------------------2nd strand zsh,wdg to local
      do kp=2,nst
        zloc=0.
        wloc=0.
        do is=1,kseq
          js=ieq(is,kp)
          if(js.gt.0) then
            zloc=zloc+het(js,3)-het(is,3)
            wloc=wloc+het(js,6)-het(is,6)
            if(ihl(js)) then
               het(js,3)=zloc
               if(isur.eq.0) het(js,6)=wloc
            endif
          else
            zloc=zloc-het(is,3)
            wloc=wloc-het(is,6)
          endif
        enddo
      enddo
c=============================================================position variables
      k=0
      ki=0
      do is=1,ntl
        in=itr(is)
        do l=1,kap(1,in)
          k=k+1
          if(.not.lock(k)) then
            ki=ki+1
            var(ki)=set(is,l)
          endif
        enddo
        k=k+kap(2,in)-kap(1,in)
      enddo
c--------------------------------------------------------------helical variables
      do is=1,nto
        do l=1,6
          k=k+1
          if(.not.lock(k)) then
             ki=ki+1
             var(ki)=het(is,l)
          endif
      enddo
      enddo
c--------------------------------------------------------------------------kinks
      do is=1,kseq
        if(kink(is)) then
          do l=1,4
            k=k+1
            if(.not.lock(k)) then
              ki=ki+1
              var(ki)=vkint(is,l)
            endif
          enddo
        endif
      enddo
c-------------------------------------------------------------------ligand inter
      do i=1,nlig
        kal=nuc(i+nto)-nuc(i+nto-1)
        lint=6
        if(kal.eq.1) lint=3
        do j=1,lint
          k=k+1
          if(.not.lock(k)) then
            ki=ki+1
            var(ki)=rlit(i,j)
          endif
        enddo
      enddo
c-----------------------------------------------------------------thymine methyl
      do is=1,nto
        if(lthy(is)) then
          k=k+1
          if(.not.lock(k)) then
            ki=ki+1
            var(ki)=tht(is)
          endif
        endif
      enddo
c---------------------------------------------------------------------superhelix
      if(isur.gt.0) var(isur)=rado
      if(isup.gt.0) var(isup)=pito
c-------------------------------------------------------------------------------
      call setbac
      call microb
      return
      end
