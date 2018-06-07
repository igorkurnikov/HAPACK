	subroutine delphi(xx,mgrid,phimap,phimap1,phimap2,phimap3,atmcrg,
     &                  iepsmp,iepsmp2,idebmap,
     &                  sf1,sf2, qmap1,qmap2,debmap1,debmap2,
     &                  bndx1,bndx2,bndx3,bndx4,bndx,bndy,bndz,
     &                  db,ibgrd,idpos,ioff,iepsv)
c------------------------------------------------------------------
	include 'qdiffpar4.h'
c------------------------------------------------------------------
      character*8 hour
      character*9 day
	character*9 bclab(4)
	logical iper,iconc
c------------------------------------------------------------------
c------------------------------------------------------------------
	dimension iper(3)
	dimension xo(3)
	dimension xn(3),xn2(3,natmax),ibgrd(4,*),dbval(0:1,0:6,0:1)
	integer*4 gchrgp(3,ngcrg)
	dimension bndx1(*),bndx2(*),bndx3(*),bndx4(*),bndx(*),bndy(*),bndz(*)
	logical   it1(0:6)
	dimension rad3(natmax),chrgv2(4,ncrgmx),idpos(*),db(6,*)
	dimension sf1((mgrid*mgrid*mgrid+1)/2),sf2((mgrid*mgrid*mgrid+1)/2)
	dimension qmap1(*),qmap2(*),debmap1(*),debmap2(*)
	dimension qval(ngcrg),iqpos(ngcrg),iepsv(*)
	dimension sfd(5,0:1)
	dimension chrgv4(natmax)
	dimension ioff(*),gchrg(ngcrg)
	dimension cqplus(3),cqmin(3)   
	dimension islice(2,mgrid) ! TEMPORAL FIX IGOR change to alloacatable...
	character enc(4)
	logical   imem,iautocon
c	offset of midlines from grid points
c-----------------------------------------------------------------

c conversion to kT/charge at 25 celsius
c
			data epkt,fpi /561.0,12.566370614359172/
c
c converts between ionic st (M) and debye length
c
	data dfact / 3.047/	
c
c boundary condition types
c
	data bclab / 'zero     ','dipolar  ',
     &             'focussing','coulombic' /
c
c------------------------------------------------------------------
c------------------------------------------------------------------
      write(6,*)'  '
	write(6,*)' FDPB CALCULATIONS START '
	write(6,*)'  '
c
c initialization
c
	call rdprm2(iautocon,epsin,epsout,rionst,exrad
     &        ,radprb,ibctyp,iper,nlit,nnit,iconc,
     &        imem,icon1,icon2)
c
c
c convert ionic strength to debye length
c
	if(rionst.ne.0.) then
	  deblen = dfact/sqrt(rionst)
	else
	  deblen = 1.e6
	end if
c
c print title, description
c

	write(6,*)'  '
	write(6,*)'grid size                  :',igrid
	
	write(6,*)'inner,outer dielectrics    :',epsin,epsout
	write(6,*)'ionic strength (M)         :',rionst
	write(6,*)'debye length (A)           :',deblen
	write(6,*)'ion exclusion radius (A)   :',exrad
	write(6,*)'probe radius (A)           :',radprb
	write(6,*)'boundary conditions        : ',bclab(ibctyp)
	write(6,*)'x,y,z periodic bc. flags   :',iper
	if(iautocon) then
	if(gten.eq.0.) write(6,*)'# of linear iterations     : automatic'
	if(gten.ne.0.) write(6,*)'convergence by grid energy :',gten,' kt'
	else 
	write(6,*)'# of linear iterations     :',nlit
	end if 
	write(6,*)'# of non-linear iterations :',nnit
	write(6,*)'concentration map output   :',iconc
	write(6,*)'map file label             :'
c
	if((icon1.ne.10).or.(icon2.ne.1)) then
	write(6,*) 'convergence test interval is every',icon1,'loops'
	write(6,*) 'testing',100/icon2,'%'
	end if 
	write(6,*)' '
	
c
c scale dielectric so that potentials turn out as kT when 
c charges are in e's and distances are in angstroms
c
	epsin = epsin/epkT
	epsout = epsout/epkT
c
c set eps and deb map values to their exterior values
c prior to resetting interior values using atom coords and radii
c
	write(6,*)'   '
	write(6,*)'initializing dielectric and debye maps...'
	write(6,*)'   '
	  do 9008 k = 1,igrid
	    do 9009 j = 1,igrid
		  do 9010 i = 1,igrid
	        phimap(i,j,k) = 0.0
		    iepsmp(i,j,k,1) = 0
		    iepsmp(i,j,k,2) = 0
		    iepsmp(i,j,k,3) = 0
		    iepsmp2(i,j,k,1) = 0
		    iepsmp2(i,j,k,2) = 0
		    iepsmp2(i,j,k,3) = 0
		    idebmap(i,j,k)  = 1
9010		 continue
9009	    continue
9008	  continue
        lim1 = (igrid*igrid*igrid + 1)
	  do i = 1,lim1
	    phimap3(i) = 0.0
	  enddo
        lim2 = (igrid*igrid*igrid + 1)/2
        do i = 1, lim2
           phimap1(i) = 0.0
	     phimap2(i) = 0.0
	  enddo
c
c jump point for shortened version
c
c initialize a few things
c
	qnet  = 0.0
	qplus = 0.0
	qmin  = 0.0
	do 9014 k  = 1,3
	  cqplus(k) = 0.0
	  cqmin(k) = 0.0
9014	continue
	nqass = 0
	natom = 0
201   format(a)

      call setatq(xn2,rad3,chrgv4,natom,oldmid,scale)

	midg=(igrid+1)/2
	rmidg=midg
c 
c scale factors and oldmid 
c 
c	fpoh = 4.0*3.14159265358979323846*scale
	fpoh = fpi*scale
c
	do 3104 ix=1,natom
	xo(1)=xn2(1,ix)
	xo(2)=xn2(2,ix)
	xo(3)=xn2(3,ix)
	call ctog(mgrid,oldmid,scale,xo,xn)
	xn2(1,ix)=xn(1)
	xn2(2,ix)=xn(2)
	xn2(3,ix)=xn(3)
3104	continue 
c assign charges for boundary conditions
c 
	ic1=0
	do 2101 ix=1,natom
	if(chrgv4(ix).ne.0.) then
	ic1=ic1+1
	atmcrg(1,ic1)=xn2(1,ix)
	atmcrg(2,ic1)=xn2(2,ix)
	atmcrg(3,ic1)=xn2(3,ix)
	atmcrg(4,ic1)=chrgv4(ix)
	end if 
2101	continue 
c 
c ic = number of charges
c 
c find charge moments for dipole approximation
c 
	qnet=0.0
	qplus=0.0
	qmin=0.0
	do 2106 ix=1,3
	cqplus(ix)=0.0
	cqmin(ix)=0.0
2106	continue 
	do 2105 ix=1,ic1
	chrg=atmcrg(4,ix) 
	qnet=qnet + chrg  
	if(chrg.gt.0.) then
	qplus=qplus + chrg
	cqplus(1)=cqplus(1) + chrg*atmcrg(1,ix)
	cqplus(2)=cqplus(2) + chrg*atmcrg(2,ix)
	cqplus(3)=cqplus(3) + chrg*atmcrg(3,ix)
	else
	qmin=qmin + chrg
	cqmin(1)=cqmin(1) + chrg*atmcrg(1,ix)
	cqmin(2)=cqmin(2) + chrg*atmcrg(2,ix)
	cqmin(3)=cqmin(3) + chrg*atmcrg(3,ix)
	end if
2105    continue 
c
c divide by charge totals
c 
	if(qplus.ne.0.0) then
	  do 2110 k = 1,3
	    cqplus(k) = cqplus(k)/qplus
2110	  continue
	end if
	if(qmin.ne.0.0) then
	  do 2111 k = 1,3
	    cqmin(k) = cqmin(k)/qmin
2111	  continue
	end if
c
c select those charges which will be charging the grid
c 
	    nqass=ic1
	    rgrid=igrid
	    ic2=0
	do 2102 ix=1,ic1
	    if((atmcrg(1,ix).gt.1.).and.(atmcrg(1,ix).lt.rgrid)) then
	    if((atmcrg(2,ix).gt.1.).and.(atmcrg(2,ix).lt.rgrid)) then
	    if((atmcrg(3,ix).gt.1.).and.(atmcrg(3,ix).lt.rgrid)) then
	      ic2=ic2+1
	      chrgv2(1,ic2)=atmcrg(1,ix)
	      chrgv2(2,ic2)=atmcrg(2,ix)
	      chrgv2(3,ic2)=atmcrg(3,ix)
	      chrgv2(4,ic2)=atmcrg(4,ix)
	    end if
	    end if
	    end if
2102	continue 
	nqgrd=ic2 
c
c write details
c
	write(6,*)'  '
	write(6,*)'scale   (grids/A): ',scale
	write(6,*)'object centre (A): ',oldmid
	write(6,*)'  '
c
c
	write(6,*)'number of atom coordinates read  : ',natom
	write(6,*)'total number of charged atoms    : ',nqass
	write(6,*)'net assigned charge              : ',qnet
	write(6,*)'assigned positive charge         : ',qplus
	write(6,*)'centred at (gu) :',cqplus
	write(6,*)'assigned negative charge         : ',qmin
	write(6,*)'centred at (gu) :',cqmin
	write(6,*)'   '
c
	call setout(xx,mgrid,
     &     phimap,phimap1,phimap2,phimap3,atmcrg,
     &     iepsmp,iepsmp2,idebmap,
     &     ioff,xn2,rad3,natom,exrad,radprb)
c
c finish off dielectric map, producing solvent accessible
c volume defined by prob radius
c 
c check for boundary elements
c
	do i=1,3
	do iz=1,igrid
	do iy=1,igrid
	do ix=1,igrid
	if(iepsmp2(ix,iy,iz,i).ne.0) iepsmp(ix,iy,iz,i)=1
	end do
	end do
	end do
	end do
c
	it1(0)= .false.
	it1(6)= .false.
	do 487 ix=1,5
	it1(ix)= .true.
487	continue
	ibnum=0
	do 671 k=2,igrid-1
	  do 672 j=2,igrid-1
	    do 673 i=2,igrid-1
              ieps=   iepsmp(i,j,k,1) +
     &   iepsmp(i,j,k,2) +
     &   iepsmp(i,j,k,3) +
     &   iepsmp(i-1,j,k,1)+ 
     &   iepsmp(i,j-1,k,2) +
     &   iepsmp(i,j,k-1,3) 
	   if(it1(ieps)) then
	      ibnum=ibnum+1
		  ibgrd(1,ibnum)=i
		  ibgrd(2,ibnum)=j
		  ibgrd(3,ibnum)=k
	   end if
673	continue
672	continue
671	continue
c
c reset iepsmap using these boundary elements
	write(6,*) "number of grid points on expanded surface= ",ibnum
c
c	if(ihs) then
c	write(6,*) "writing expanded surface data file: hsurf.dat"
c	open(40,file="hsurf.dat")
c	sahen=1.2079/(scale**2)
c	do i=1,ibnum
c	xo(1)=ibgrd(1,i)
c	xo(2)=ibgrd(2,i)
c	xo(3)=ibgrd(3,i)
c	j=ibgrd(4,i)
c	call gtoc(mgrid,oldmid,scale,xo,xn)
c	write(40,414) j,xn,sahen
c	end do
c414	format(i5,3f8.3,f8.3)
c	close(40)
c	end if
c
	call mkeps(xx,mgrid,
     &   phimap,phimap1,phimap2,phimap3,atmcrg,	
     &   iepsmp,iepsmp2,idebmap,
     &   radprb,epsout,ibnum,ibgrd,rad3,xn2)
c
c have used iepsmp2, now reassign iepsmp2
c
	call setin(xx,mgrid,phimap,phimap1,phimap2,phimap3,atmcrg,
     &           iepsmp,iepsmp2,idebmap,
     &           ioff,xn2,rad3,natom)
c
c COMMENT OUT MEMBRANES FOR THE MOMENT
	 if(imem) call mem(xx,mgrid,phimap,phimap1,phimap2,phimap3,atmcrg,
     &                   iepsmp,iepsmp2,idebmap,islice )
c 
c recalculate boundary elements
c
	debfct = epsout/(deblen*scale)**2
	difeps=epsin-epsout
	sixeps=epsout*6.0
	sixth=1.0/6.0
	sixsalt=sixth*((1/(1+debfct/sixeps))-1.0)
	icgrid=igrid*igrid*igrid 
c
	if(rionst.ne.0.) then
	  do 491 iz=0,1
	    do 493 ix=1,3
	      denom= sixeps + ix*difeps + iz*debfct
	      dbval(0,ix,iz)= 0.0
	      dbval(1,ix,iz)= difeps/denom
	      sfd(ix,iz)=epsout/denom
493	    continue
491	  continue
	  do 492 iz=0,1
	    do 494 ix=4,5
	      denom= sixeps + ix*difeps + iz*debfct
	      dbval(1,ix,iz)= 0.0
	      dbval(0,ix,iz)= -difeps/denom
	      sfd(ix,iz)=epsin/denom
494	    continue
492	  continue
	else
	  do 591 iz=0,1
	    do 593 ix=1,5
	      denom= sixeps + ix*difeps 
	      dbval(0,ix,iz)= (epsout/denom) -sixth 
	      dbval(1,ix,iz)= (epsin/denom) - sixth
593	    continue
591	  continue
	end if
c
c dbval for db
c 
	ibnum=0
	isgrid=igrid**2
	do 771 k=2,igrid-1
	  do 772 j=2,igrid-1
	    do 773 i=2,igrid-1
            ieps =  (iepsmp(i,j,k,1) + 
     &	                iepsmp(i,j,k,2) + 
     &	                iepsmp(i,j,k,3) + 
     &	                iepsmp(i-1,j,k,1) + 
     &	                iepsmp(i,j-1,k,2) + 
     &	                iepsmp(i,j,k-1,3)) 
	     if(it1(ieps)) then
	  	     ibnum=ibnum+1
		     ibgrd(1,ibnum)=i
		     ibgrd(2,ibnum)=j
		     ibgrd(3,ibnum)=k
		     ioff(ibnum)=ieps
	     end if
773	continue
772	continue
771	continue
c
	ibnum1=0
	nsp2 = mgrid*mgrid*5
	ibnum2=nsp2
	do 774 ix=1,ibnum
	  i=ibgrd(1,ix)
	  j=ibgrd(2,ix)
	  k=ibgrd(3,ix)
	  ieps=ioff(ix)
	  iw=isgrid*(k-1) + igrid*(j-1) + i
	  iv=(iw+1)/2
	  deb=idebmap(i,j,k)
	  if(iw.ne.(2*iv)) then
		ibnum1=ibnum1+1
		ibnum3=ibnum1
	  else
		ibnum2=ibnum2+1
		ibnum3=ibnum2
	  end if
	  idpos(ibnum3) = iv
	  iepsv(ibnum3) = ieps
	  db(1,ibnum3)=dbval(iepsmp(i-1,j,k,1),ieps,deb)
	  db(2,ibnum3)=dbval(iepsmp(i,j,k,1),ieps,deb)
	  db(3,ibnum3)=dbval(iepsmp(i,j-1,k,2),ieps,deb)
	  db(4,ibnum3)=dbval(iepsmp(i,j,k,2),ieps,deb)
	  db(5,ibnum3)=dbval(iepsmp(i,j,k-1,3),ieps,deb)
	  db(6,ibnum3)=dbval(iepsmp(i,j,k,3),ieps,deb)
774	continue
c realign idpos and db
c
	icount2a=ibnum1
	icount2b=icount2a+ibnum2- nsp2
	itemp= nsp2
	do 781 ix=icount2a+1,icount2b
	  itemp=itemp+1
	  idpos(ix)=idpos(itemp)
	  iepsv(ix)=iepsv(itemp)
	  db(1,ix)=db(1,itemp)
	  db(2,ix)=db(2,itemp)
	  db(3,ix)=db(3,itemp)
	  db(4,ix)=db(4,itemp)
	  db(5,ix)=db(5,itemp)
	  db(6,ix)=db(6,itemp)
781	continue
c
c
	write(6,*) 'number of dielectric boundary points', icount2b 
c
	call chrgup(xx,mgrid,
     &    phimap,phimap1,phimap2,phimap3,atmcrg,
     &    iepsmp,iepsmp2,idebmap,
     &    nqgrd,chrgv2,qval,iqpos,fpoh,sixeps,difeps,
     &                 debfct,gchrgp,gchrg,icount1a,icount1b)
 
c	
	write(6,*) 'number of grid points assigned charge', icount1b
c
c set saltmaps 1 and 2
	if(rionst.ne.0.) then 
	iw=1
	do 841 iz=1,igrid
	  do 842 iy=1,igrid
	    do 843 ix=1,igrid
		  phimap3(iw)=sixth + idebmap(ix,iy,iz)*sixsalt  
	      iw=iw+1
843	continue
842	continue
841	continue
	iy=0
	sf1((icgrid+1)/2)=phimap3(icgrid)
	do 850 ix=1,icgrid-2,2
	iy=iy+1
	sf1(iy)=phimap3(ix)
	sf2(iy)=phimap3(ix+1) 
850	continue 
	do 844 ix=1,icount2a
	  i=1
	  if(sf1(idpos(ix)).eq.sixth) i=0
	  sf1(idpos(ix))=sfd(iepsv(ix),i)
844	continue
	do 845 ix=icount2a+1,icount2b
	  i=1
	  if(sf2(idpos(ix)).eq.sixth) i=0
	  sf2(idpos(ix))=sfd(iepsv(ix),i)
845	continue
	end if 
c
c
c calculate boundary conditions
c
c
	write(6,*)'  '
	write(6,*)'setting boundary conditions...'
	write(6,*)'  '
c
c
	call setbc(xx,mgrid,
     &          phimap,phimap1,phimap2,phimap3,atmcrg,
     &          iepsmp,iepsmp2,idebmap,
     &          ibctyp,iper,qplus,qmin,cqplus,cqmin,
     &          epsout,deblen,nqass)
c
c
c iterate
c
	call relfac(xx,mgrid,
     &   phimap,phimap1,phimap2,phimap3,atmcrg,
     &   iepsmp,iepsmp2,idebmap,
     &   bndx,bndy,bndz,
     &   nlit,nnit,iper,idpos,db,sf1,sf2,
     &            icount2a,icount2b,rionst,spec)
c 
	noit=int(7.8/log(1.0 + sqrt(1-spec)))
	write(6,*) 'estimated iterations to convergence',noit
	if(iautocon) nlit = noit 
c
c
	if(nnit.eq.0.or.rionst.eq.0.0) then
	call itit(xx,mgrid,
     & phimap,phimap1,phimap2,phimap3,atmcrg,
     & iepsmp,iepsmp2,idebmap,
     & nlit,nnit,iper,idpos,db,
     & bndx1,bndx2,bndx3,bndx4,bndx,bndy,bndz,
     & sf1,sf2,iqpos,qval,
     & icount2a,icount2b,icount1a,icount1b,rionst,spec,
     & icon1,icon2,epsin,epsout)
c
	else 
c
	call nitit(xx,mgrid,
     &  phimap,phimap1,phimap2,phimap3,atmcrg,
     &  iepsmp,iepsmp2,idebmap,
     &  qmap1,qmap2,debmap1,debmap2,
     &  bndx1,bndx2,bndx3,bndx4,bndx,bndy,bndz,
     &  nlit,nnit,iper,idpos,db,sf1,sf2,iqpos,qval,
     &            icount2a,icount2b,icount1a,icount1b,rionst,spec,
     &            icon1,icon2,debfct,epsout)
	end if 
c
	ergg=0.0
	do 587 i=1,icount1b
	ix=gchrgp(1,i)
	iy=gchrgp(2,i)
	iz=gchrgp(3,i)
	ergg=ergg + phimap(ix,iy,iz)*gchrg(i)
587	continue
	ergg=ergg/2.0
	write(6,*) ' '
	write(6,*) 'total energy (including grid energy): ',ergg,' kt'

      call geteneq(ergg)

	if(iconc) then
c
c convert potentials to concentrations
c
	  sixth = 1./6.
	  if(rionst.ne.0.0) then
	  write(6,*)'  '
	  write(6,*)'converting potentials to '
	  write(6,*)'net charge concentrations...'
	  write(6,*)'  '
	    do 9037 iz = 1,igrid
	      do 9038 iy = 1,igrid
		  do 9039 ix = 1,igrid
c
c use same number of terms in expansion
c of sinh as for iteration in itit.f
c
c		    temp = exp(-phimap(ix,iy,iz)) - exp(phimap(ix,iy,iz))
		    phi = phimap(ix,iy,iz)
		    phisq = phi**2
		    temp = phisq/120. + sixth
		    temp = phisq*temp + 1.
		    temp = 2.0*temp*phi
		    phimap(ix,iy,iz) = rionst*temp*idebmap(ix,iy,iz)
9039	    continue
9038	  continue
9037	continue
	  end if
	else
	end if
	
	return
	end
