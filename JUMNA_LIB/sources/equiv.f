      subroutine equivjm(nspc,nbcc,nhlc,nklc,nlgc)
      include 'jumna_data.inc'
      logical*2 lar,lock,kink,lthy,hst,bst,vst,lgi,lgj,lgu,locr,cation,
     1 sup,rcom,homo,homo2,homo3,diep,link,ecen,cyl,ribose,lcat,
     1 cent,brk(n2+1),ok,autos,ihl,hom(3),hseq(3),amber,sum
      character*4 mnam,munit,seq*120,code*8,kode*8,lnam,rb*1,ct*1,ch(18)
      integer*4 opt,ipr(n2),ipt(2,3,n2),nbk(0:3),
     1 sy(n2+1),mho(3),nuch(n2),newq(n2,18)
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
      common/mnn/var(n7),gra(n7),con(n0*n2+n3),scl(n7),ncon,nvrc,ntba,
     1 nbac,nthe,nhel,ntki,nkin,ntri,ntot,nvar,nrin,lar(n7),lock(n6a)
      common/mrc/corm(n1,3),dmon(n1),mnam(n1),munit(n1),
     1 imch(n1),imty(n1),icm(n1),matm(3*n1),matd(n1,7),nunit(n1),
     1 nuc(0:n2),ncen(0:n0*n2),kam,khm,lkm,kcen
      common/slg/hst(n2,6),bst(n2,n8),vst(n2,4),lgi(n1),lgj(n1),lgu(n2)
      common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
     1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(120,4),ilq(n2,2),idr(4),
     1 kink(n2),lthy(n2)
      common/symjm/xsum,ysum,itwl(n2),nvs(n7),ihm(n2),ksym(3),nbrk(2),
     1 isym,isur,isup,ihl(n2)
	common/hacon/iret
      equivalence (hom(1),homo),(mho(1),mhomo)
      data ch/'C1''','C2''','C3''','Gly','C12','C23','Eps','Zet',
     1'O3''','P  ','Rib','Thy','Xdis','Ydis','Rise','Inc','Tip','Twis'/
      k=0
      do i=1,nto
        in=itr(i)
        ino=ito(i)
        ienb(i,1)=k+1
        ienb(i,2)=k+kap(1,ino)
        k=k+kap(2,in)
      do iv=1,18
        newq(i,iv)=0
      enddo
      enddo
c---------------------------------------------------------------update ksym,nbrk
      if(nbrk(1).gt.kseq.or.nbrk(2).gt.kseq) then
      write(6,509) kseq+1
509   format('   NBRK MUST BE LESS THAN',i4)
      istop=1
      endif
      nbk(0)=1
      do i=1,2
      if(nbrk(i).eq.0) nbrk(i)=kseq+1
      nbk(i)=nbrk(i)
      enddo
        nbk(3)=kseq+1
        if(ksym(1).eq.0) ksym(1)=2*kseq
        if(ksym(2).eq.0.and.nbrk(1).le.kseq) ksym(2)=2*kseq
        if(ksym(3).eq.0.and.nbrk(2).le.kseq) ksym(3)=2*kseq
      do j=1,kseq+1
        brk(j)=.false.
        if(j.eq.nbrk(1).or.j.eq.nbrk(2)) brk(j)=.true.
        sy(j)=1
        if(j.ge.nbrk(1)) sy(j)=2
        if(j.ge.nbrk(2)) sy(j)=3
      enddo
      nsli=sy(kseq)
      istop=0
      nsth=nst
      if(hom(1).or.hom(2).or.hom(3)) nsth=nst/2
c------------------------------------------------------------mhomo from sequence
      if(nsth.ne.nst) then
      do ns=1,nsli
        hseq(ns)=.false.
      if(hom(ns).and.mho(ns).eq.-1) then
        hseq(ns)=.true.
        mho(ns)=0
      endif
      enddo
        if(hseq(1).or.hseq(2).or.hseq(3)) write(6,*)
        if(nst.eq.4) goto 514
        ij=idr(1)*idr(2)
      do j=1,kseq
        ns=sy(j)
        if(.not.hseq(ns).or.ksym(ns).eq.2*kseq) goto 510
        l1=j
7     if(ij.gt.0) then
        j2=j+mho(ns)
        if(j2.ge.nbk(ns)) j2=j2-ksym(ns)
      elseif(ij.lt.0) then
        j2=nbk(ns-1)+nbk(ns)-1-j-mho(ns)
        if(j2.lt.nbk(ns-1)) j2=j2+ksym(ns)
      endif
        l2=kseq+j2
      if(seq(l1:l1).ne.seq(l2:l2)) then
        mho(ns)=mho(ns)+1
        if(mho(ns).lt.ksym(ns)) goto 7
      endif
510   enddo
      do ns=1,3
        if(hseq(ns).and.mho(ns).ge.ksym(ns)) mho(ns)=0
      enddo
514     if(hseq(1)) write(6,511) mho(1)
        if(hseq(2)) write(6,512) mho(2)
        if(hseq(3)) write(6,513) mho(3)
511     format(2x,'... mhomo  reset to ',i2,' by Equiv')
512     format(2x,'... mhomo2 reset to ',i2,' by Equiv')
513     format(2x,'... mhomo3 reset to ',i2,' by Equiv')
      endif
        write(6,*)
c----------------------------------------------missing nucleotides and nbrk test
      do i=1,kseq
        ipr(i)=i
      enddo
      do is=2,nst
        jl=ilq(iend(is-1)+1,1)
        ju=ilq(iend(is),1)
        ine=1
        ok=.false.
      do j=jl,ju
        in=ieq(j,is)
        ns=sy(j)
      if(in.eq.0) then
        if(ine.ne.0) j1=j
        if(ine.ne.0.and.ksym(ns).eq.2*kseq) miss=1
        if(ine.eq.0.and.ksym(ns).eq.2*kseq) miss=miss+1
        if(brk(j).or.brk(j+1)) ok=.true.
      endif
      if(in.ne.0.and.ine.eq.0) then
        if(miss.eq.j-j1) ok=.true.
      if(.not.ok) then
        j2=j
        write(6,30)(ipr(i),i=j1,j2)
30      format('   MISSING NUCLEOTIDE(S): NO SYMMETRY AT THEIR LEVELS,'
     1        /'  OR ONE SYMMETRY BREAK MUST BE GIVEN '
     1        /'  AT ANY LEVEL',10i3)
        istop=1
      endif
        ok=.false.
      endif
        ine=in
      enddo
      enddo
c-----------------------------------------------------------------error messages
      do ns=1,3
        if(ns.eq.1) ct=' '
        if(ns.gt.1) write(ct,516)ns
      if(mho(ns).ge.ksym(ns)) then
        write(6,517)ct,ns
        istop=1
      endif
      if(hom(ns).and.nbk(ns-1).gt.kseq) then
        write(6,518)ct,ns-1,kseq+1
        istop=1
      endif
      enddo
516   format(i1)
517   format('   MHOMO',a1,' MUST BE LESS THAN KSYM(',i1,')')
518   format('   IF HOMO',a1,'=.T. NBRK(',i1,') MUST BE LESS THAN',i4)
      if(ecen.and.nbrk(1).le.kseq) then
        write(6,*)'   SYMMETRY BREAK(S) NOT ALLOWED WITH ECEN'
        istop=1
      endif
      if(nbrk(1).gt.nbrk(2).and.nbrk(2).le.kseq) then
        write(6,*)'   NBRK(1) MUST BE LESS THAN NBRK(2)'
        istop=1
      endif
      if(nbrk(1).lt.0.or.nbrk(2).lt.0.or.
     1   nbrk(1).eq.1.or.nbrk(2).eq.1) then
        write(6,*)'   NBRK MUST BE  0 OR GREATER THAN  1'
        istop=1
      endif
      do ns=1,nsli
      if(ksym(ns).lt.2*kseq.and.
     1   ksym(ns).gt.nbk(ns)-nbk(ns-1)) then
        write(6,519) ns,nbk(ns)-nbk(ns-1)
        istop=1
      endif
      enddo
519   format('   KSYM(',i1,') CANNOT BE GREATER THAN',i4)
      do i=1,3
      if(ksym(i).lt.0) then
        write(6,*)'   KSYM MUST BE ZERO OR POSITIVE'
        istop=1
      endif
      enddo
      if(nbrk(1).gt.kseq.and.(ksym(2).gt.0.or.ksym(3).gt.0)) then
        write(6,*)'   KSYM(2) AND KSYM(3) MUST BE ZERO'
        istop=1
      endif
      if(nbrk(2).gt.kseq.and.ksym(3).gt.0) then
        write(6,*)'   KSYM(3) MUST BE ZERO'
        istop=1
      endif
      if(nsth.ne.nst) then
      if(nst.eq.1.or.nst.eq.3) then
        write(6,*)'   HOMO=.T. ALLOWED ONLY WITH 2 OR 4 STRANDS'
        istop=1
      endif
      do is=nsth+1,nst
        iss=is-nsth
        ij=idr(is)*idr(iss)
      do j=1,kseq
        ns=sy(j)
        if(.not.hom(ns)) goto 9
        if(ij.lt.0) js=nbk(ns-1)+nbk(ns)-1-j-mho(ns)
        if(ij.gt.0) js=j-mho(ns)
        if(js.lt.nbk(ns-1)) js=js+ksym(ns)
        if(js.gt.nbk(ns)) goto 9
        ct=seq(j+(is-1)*kseq:j+(is-1)*kseq)
        rb=seq(js+(iss-1)*kseq:js+(iss-1)*kseq)
      if((ct.eq.'-'.or.rb.eq.'-').and.ct.ne.rb) then
        write(6,*)'   MISSING NUCLEOTIDES MUST CORRESPOND'
        write(6,*)'  IN HOMOLOGOUS STRANDS (OR SECTIONS)'
        istop=1
        goto 10
      endif
9     enddo
      enddo
      endif
10    do j=1,kseq
        ns=sy(j)
      if(kink(j)) then
      if(.not.brk(j).and.ksym(ns).ne.2*kseq) then
        write(6,*) '   KINKS ARE ALLOWED ONLY IF NO SYMMETRY '
        write(6,*) '  OR AT SYMMETRY BREAKS'
        istop=1
      endif
      if(nsth.ne.nst) then
      if(.not.brk(j).and.mho(ns).gt.0) then
        write(6,*) '   KINKS ARE ALLOWED WITH HOMO=.T. ONLY IF MHOMO= 0'
        istop=1
      endif
      if((idr(nsth+1).lt.0.or.idr(2)*idr(4).lt.0).and.
     1 .not.kink(kseq+2-j)) then
        write(6,*) '   KINK(S) AND HOMO=.T. WITH ANTIPARALLEL STRANDS'
        write(6,*) '  ARE ALLOWED ONLY IF KINK(S) POSITION(S)'
        write(6,*) '  IS(ARE) COMPATIBLE WITH HOMO=.T.'
        istop=1
      endif
      endif
      endif
      enddo
c-----------------------------------------------------------------central energy
      write(6,*) 'ecen = ',ecen
      if(ecen) then
        do i=1,nto
          lgu(i)=.false.
        enddo
        do i=1,kam
          lgi(i)=.false.
          lgj(i)=.false.
        enddo
        do kp=1,nst
          do in=2,ksym(1)+1
            is=ieq(in,kp)
            if(is.gt.0) then
               lgu(is)=.true.
               do i=nuc(is-1)+1,nuc(is)
                 lgi(i)=.true.
                 lgj(i)=.true.
               enddo
            endif
          enddo
          do in=ksym(1)+2,kseq
            is=ieq(in,kp)
            do i=nuc(is-1)+1,nuc(is)
              lgj(i)=.true.
            enddo
          enddo
        enddo
      else
        do i=1,kcen
          lgu(i)=.true.
        enddo
        do i=1,kam
          lgi(i)=.true.
          lgj(i)=.true.
        enddo
      endif
c-----------------------------------------------search for homologous nucleotide
      do is=nsth+1,nst
        ij=idr(is)*idr(is-nsth)
      do j=1,kseq
        in=ieq(j,is)
        if(in.eq.0) goto 2
        ns=sy(j)
        if(.not.hom(ns)) goto 2
        if(ij.gt.0) js=j-mho(ns)
        if(ij.lt.0) js=nbk(ns-1)+nbk(ns)-1-j-mho(ns)
        if(js.lt.nbk(ns-1)) js=js+ksym(ns)
        ins=0
        if(js.lt.nbk(ns)) ins=ieq(js,is-nsth)
        nuch(in)=ins
2     enddo
      enddo
ccc----------------------------------------------------------backbone variables
        k=0
      do l=1,nsli
      do i=1,ksym(l)
        ipt(1,l,i)=0
        ipt(2,l,i)=0
      enddo
      enddo
      do is=1,nst
        ii=sign(1,idr(is))
        j5=1
        if(ii.lt.0) j5=kseq
        j3=j5+ii*(kseq-1)
        id=(1-ii)/2
        n=-1
        ipair=1
        if(nst-nsth.eq.2) ipair=mod(is-1,nsth)+1
      do j=j5,j3,ii
        if(brk(j+id)) n=-1
        in=ieq(j,is)
        if(in.eq.0) goto 19
        ien1=ienb(in,1)
        ien2=ienb(in,2)
        if(ribose(in)) ien2=ien2-1
c----------------------------------------------------------symmetric nucleotide
        n=n+1
        ns=sy(j)
        jm=mod(n,ksym(ns))+1
        ins=0
        if(n.ge.ksym(ns)) ins=in-ii*ksym(ns)
      if(hom(ns).and.is.gt.nsth) then
        ins=nuch(in)
        jm=mod(n-mho(ns),ksym(ns))+1
        if(jm.le.0) jm=jm+ksym(ns)
      endif
      if(ins.ne.0) then
        ies1=ienb(ins,1)
        ies2=ienb(ins,2)
        if(ribose(ins)) ies2=ies2-1
        if(ribose(in).neqv.ribose(ins)) write(6,24) ch(11),ins,in
      endif
c--------------------------------------------------------newq for each variable
        iv=0
      do i=ien1,ien2
        l=ies1+iv
        iv=iv+1
      if(ins.eq.0.or.l.gt.ies2.or.(brk(j+1-id).and.iv.ge.7)) then
      if(.not.lock(i)) then
        k=k+1
        newq(in,iv)=k
      endif
      else
      if(lock(i).neqv.lock(l)) then
        write(6,22) ch(iv),ins,in
        istop=1
      endif
        newq(in,iv)=-abs(newq(ins,iv))
      endif
      enddo
c----------------------------------------------------------newq for Rib variable
      if(ribose(in)) then
      if(.not.lock(ien2+1)) then
      if(ipt(ipair,ns,jm).eq.0) then
        k=k+1
        newq(in,11)=k
        ipt(ipair,ns,jm)=k
      else
        newq(in,11)=-ipt(ipair,ns,jm)
      endif
      endif
        if(ins.ne.0.and.ribose(ins).and.(lock(ies2+1).neqv.
     1  lock(ien2+1))) write(6,23) ch(11),ins,in
      endif
19    enddo
      do l=1,nsli
      do i=1,ksym(l)
        if(.not.hom(l)) ipt(ipair,l,i)=0
      enddo
      enddo
      enddo
      nspc=k
22    format(/'  LOCKING of variable ',a4,
     1  ' is NOT COMPATIBLE in identical nucleotides',i3,' and',i3)
23    format( '  locking of variable ',a4,
     1  ' is different in identical nucleotides',i3,' and',i3)
24    format( '  existence of variable ',a4,
     1  ' is different in identical nucleotides',i3,' and',i3)
ccc---------------------------------------------------------helicoidal variables
      do is=1,nst
        if(is.le.nsth) ij=0
        if(is.gt.nsth) ij=idr(is)*idr(is-nsth)
        iss=is
        n=-1
      do j=kseq,1,-1
        if(brk(j+1)) n=-1
        in=ieq(j,is)
        if(in.eq.0) goto 29
        kt=ntba+6*(in-1)
c-----------------------------------------------------------symmetric nucleotide
        n=n+1
        ns=sy(j)
        ins=0
        if(n.ge.ksym(ns)) ins=in+ksym(ns)
        if(hom(ns).and.ij.ne.0) ins=nuch(in)
        ks=ntba+6*(ins-1)
c---------------------------------------------------------newq for each variable
      do iv=1,5
        if(iv.eq.3) goto 27
      if(ins.eq.0.or.(hom(ns).and.ij.gt.0.and.iv.ne.2)) then
      if(.not.lock(kt+iv)) then
        k=k+1
        newq(in,iv+12)=k
      endif
      else
      if(lock(kt+iv).neqv.lock(ks+iv)) then
        write(6,22) ch(iv+12),ins,in
        istop=1
      endif
        newq(in,iv+12)=-abs(newq(ins,iv+12))
      endif
27    enddo
c-------------------------------------------------------newq for Rise and Twist
      if(hom(ns).and.ins.ne.0) then
        if(ij.lt.0) ins=ins+1
        if(ihl(in).or.ihl(ins)) ins=0
      endif
        if(ins.eq.1) ins=1+ksym(1)
        if(hom(ns).and.ij.ne.0) iss=is-nsth
        if(ins.gt.ieq(nbk(ns)-1,iss)) ins=0
        ks=ntba+6*(ins-1)
      do iv=3,6,3
      if(in.eq.1.and.(iv.eq.3.or.(isur.eq.0.and.iv.eq.6))) goto 28
      if(ins.eq.0.or.(is.eq.1.and.brk(j))) then
      if(.not.lock(kt+iv)) then
        k=k+1
        newq(in,iv+12)=k
      endif
      else
      if(lock(kt+iv).neqv.lock(ks+iv)) then
        write(6,25) ch(iv+12),ins,in
        istop=1
      endif
        newq(in,iv+12)=-abs(newq(ins,iv+12))
      endif
28    enddo
29    enddo
      enddo
25    format(/'  LOCKING of variable ',a4,
     1  ' is NOT COMPATIBLE in correspd. nucleotides',i3,' and',i3)
ccc---------------------------------------------------------------thymine methyl
        kt=ntlg
      do l=1,nsli
      do i=1,ksym(l)
        ipt(1,l,i)=0
        ipt(2,l,i)=0
      enddo
      enddo
      do is=1,nst
        ii=sign(1,idr(is))
        j5=1
        if(ii.lt.0) j5=kseq
        j3=j5+ii*(kseq-1)
        id=(1-ii)/2
        n=-1
        ipair=1
        if(nst-nsth.eq.2) ipair=mod(is-1,nsth)+1
      do j=j5,j3,ii
        if(brk(j+id)) n=-1
        in=ieq(j,is)
        if(in.eq.0) goto 39
c----------------------------------------------------------symmetric nucleotide
        n=n+1
        ns=sy(j)
        jm=mod(n,ksym(ns))+1
        ins=0
        if(n.ge.ksym(ns)) ins=in-ii*ksym(ns)
      if(hom(ns).and.is.gt.nsth) then
      ins=nuch(in)
      jm=mod(n-mho(ns),ksym(ns))+1
      if(jm.le.0) jm=jm+ksym(ns)
      endif
      if(ins.eq.0) goto 38
      if(lthy(in).neqv.lthy(ins)) write(6,24) ch(12),ins,in
c--------------------------------------------------------------------------newq
38    if(lthy(in)) then
      kt=kt+1
      ipr(in)=kt
      if(.not.lock(kt)) then
        if(ipt(ipair,ns,jm).eq.0) then
        k=k+1
        newq(in,12)=k
        ipt(ipair,ns,jm)=k
          else
          newq(in,12)=-ipt(ipair,ns,jm)
          endif
      endif
      if(ins.eq.0) goto 39
      if(lthy(ins).and.(lock(ipr(ins)).neqv.lock(kt))) write(6,23)
     1 ch(12),ins,in
      endif
39    enddo
      do l=1,nsli
      do i=1,ksym(l)
      if(.not.hom(l)) ipt(ipair,l,i)=0
      enddo
      enddo
      enddo
c------------------------------------------------------------
      call equim(newq,nspc,nhlc,nbcc,nklc,nlgc,istop)
      if(istop.eq.1) iret = 1
      return
      end
