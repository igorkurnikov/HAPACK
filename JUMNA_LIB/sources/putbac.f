      subroutine putbac(key)
      include 'jumna_data.inc'
      logical*2 lar,lock,kink,lthy,locr,sup,rcom,homo,homo2,homo3,
     1 diep,sum,link,ecen,cyl,lcat,cent,autos,amber
      character*4 seq*120,code*8,kode*8,lnam,mnam,munit
      integer*2 i23,i34,elim
      integer*4 opt
      dimension ha(n2,9)
      common/datjm/acc,phos,epsi,epsr,plat,slope,rhbl,vfac,tfac,rfac,xfac,
     1 scale,rpiv,tpiv,fad,fan,damp,fnoew,fnoes,fnoea,df1,scnb,scee,
     1 catd,catr,catc,rad,pit,enit,nshel,maxn,opt,mhomo,mhomo2,mhomo3,
     1 nick,limit,nop,nrib,ncat,lig,naxo,sup,rcom,homo,homo2,homo3,diep,
     1 sum,link,ecen,cyl,lcat,cent,autos,amber
      common/enf/for(n1,3),tor(n1,3),fot(n6),ener,elec,
     1 repl,disp,eang,etog,epen,i23(8*n1),i34(8*n1),elim(n1)
      common/flx/sap(n6),refg,refb,refh,refv,iap(n8,n4,n5),
     1 nap(n8,7,n5),kap(3,n5)
      common/hel/ua(n2,3),da(n2,3),ra(n2,3)
      common/lgd/rlig(n9,6),slig(n9,2),ilig(n9,6),lopt(n9),lpiv(n9),
     1 ntlg,nlgi,nlig,ntl,lnam(n9),locr(n2,n0)
      common/mnn/var(n7),gra(n7),con(n0*n2+n3),scl(n7),ncon,nvrc,
     1 ntba,nbac,nthe,nhel,ntki,nkin,ntri,ntot,nvar,nrin,
     1 lar(n7),lock(n6a)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/srs/corms(n1,3),saps(n6),vars(n7),gras(n7),hels(n2,6),
     1 vkis(n2,4),has(n2,9),rlis(n9,6),eref,rref,pref
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(120,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      equivalence (ha(1,1),ua(1,1))
c--------------------------------------------------------------save conformation
      if(key.eq.0) then
        do i=1,kam
          do j=1,3
            corms(i,j)=corm(i,j)
          enddo
        enddo
        do i=1,ntot
          saps(i)=sap(i)
        enddo
        do i=1,nvar
          vars(i)=var(i)
          gras(i)=gra(i)
        enddo
        do i=1,nto
          do j=1,6
            hels(i,j)=hel(i,j)
          enddo
        enddo
        do i=1,kseq
          do j=1,4
            vkis(i,j)=vkink(i,j)
          enddo
        enddo
        do i=1,nlig
          do j=1,6
             rlis(i,j)=rlig(i,j)
          enddo
        enddo
        do i=1,nto
          do j=1,9
             has(i,j)=ha(i,j)
          enddo
        enddo
        rref=rad
        pref=pit
        eref=ener
c-----------------------------------------------------------restore conformation
      else
        do i=1,kam
          do j=1,3
             corm(i,j)=corms(i,j)
          enddo
        enddo
        do i=1,ntot
          sap(i)=saps(i)
        enddo
        do i=1,nvar
          var(i)=vars(i)
        enddo
        do i=1,nto
          do j=1,6
             hel(i,j)=hels(i,j)
          enddo
        enddo
        do i=1,kseq
          do j=1,4
             vkink(i,j)=vkis(i,j)
          enddo
        enddo
        do i=1,nlig
          do j=1,6
             rlig(i,j)=rlis(i,j)
          enddo
        enddo
        do i=1,nto
          do j=1,9
            ha(i,j)=has(i,j)
          enddo
        enddo
        rad=rref
        pit=pref
      endif
      return
      end
