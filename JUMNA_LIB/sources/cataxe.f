      subroutine cataxe(hcat,icat)
      include 'jumna_data.inc'
      logical*2 locr,kink,lthy,ribose,cation,sup,rcom,homo,
     1 sum,homo2,homo3,diep,link,ecen,cyl,lcat,cent,autos,amber
      character*4 mnam,munit,lnam,seq*120,code*8,kode*8
      integer*4 opt
      dimension cob(3,3),cdir(9),hcat(n2,3),icat(n2)
      common/datjm/acc,phos,epsi,epsr,plat,slope,rhbl,vfac,tfac,rfac,xfac,
     1 scale,rpiv,tpiv,fad,fan,damp,fnoew,fnoes,fnoea,df1,scnb,scee,
     1 catd,catr,catc,rad,pit,enit,nshel,maxn,opt,mhomo,mhomo2,mhomo3,
     1 nick,limit,nop,nrib,ncat,lig,naxo,sup,rcom,homo,homo2,homo3,diep,
     1 sum,link,ecen,cyl,lcat,cent,autos,amber
      common/extjm/thy(n2),rsr(n0,n2),iofs(n2),iofe(n2),ithy(6),nith,
     1 neq(n7),isr(n0,2,n2),nsr(n2),ribose(n2),cation(n2)
      common/hel/ua(n2,3),da(n2,3),ra(n2,3)
      common/lgd/rlig(n9,6),slig(n9,2),ilig(n9,6),lopt(n9),lpiv(n9),
     1 ntlg,nlgi,nlig,ntl,lnam(n9),locr(n2,n0)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(120,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
c---------------------------------------------------------------find link params
      k=0
      do is=1,nto
      if(cation(is)) then
      i1=ilq(is,1)
      k=k+1
      icat(k)=i1
      nc=nuc(is)
      dx=-da(i1,1)
      dy=-da(i1,2)
      dz=-da(i1,2)
      wx=ua(i1,2)*dz-dy*ua(i1,3)
      wy=ua(i1,3)*dx-dz*ua(i1,1)
      wz=ua(i1,1)*dy-dx*ua(i1,2)
      do j=1,3
      cob(1,j)=ra(i1,j)
      cob(2,j)=ra(i1,j)-da(i1,j)
      enddo
      cob(3,1)=ra(i1,1)+wx
      cob(3,2)=ra(i1,2)+wy
      cob(3,3)=ra(i1,3)+wz
      call cosdir(cob,3,1,2,3,cdir)
      x=corm(nc,1)-cob(1,1)
      y=corm(nc,2)-cob(1,2)
      z=corm(nc,3)-cob(1,3)
      hcat(k,1)=cdir(1)*x+cdir(4)*y+cdir(7)*z
      hcat(k,2)=cdir(2)*x+cdir(5)*y+cdir(8)*z
      hcat(k,3)=cdir(3)*x+cdir(6)*y+cdir(9)*z
      endif
      enddo
      ncat=k
      return
      end
