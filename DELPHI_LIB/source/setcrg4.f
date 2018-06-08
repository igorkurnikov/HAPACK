	subroutine chrgup(xx,mgrid,
     &  phimap,phimap1,phimap2,phimap3,atmcrg,
     &  iepsmp,iepsmp2,idebmap,
     &  nqgrd,chrgv2,qval,iqpos,fpoh,
     & sixeps,difeps,debfct,gchrgp,gchrg,icount1a,icount1b)
c
	include 'qdiffpar4.h'
c 
	logical isum(3*ngrmax)
	dimension gchrg(ngcrg),gchrgd(ngcrg),chrgv2(4,ncrgmx)
	dimension qval(ngcrg),iqpos(ngcrg)
    	integer gchrgp(3,ngcrg)
	integer*4 gchrg2(ngcrg)
c
c assign fractional charges to grid points, using phimap
c as a dummy array (set to zero again later)
c 
	do 821 ig=1,nqgrd 
	  kb1=chrgv2(1,ig)
	  kb2=chrgv2(2,ig)
	  kb3=chrgv2(3,ig)
	  do 822 ix=0,1
	    i=kb1+ix
	    cg1=kb1-chrgv2(1,ig)+1-ix
	    do 823 iy=0,1
	      j=kb2+iy
	      cg2=kb2-chrgv2(2,ig)+1-iy
		  do 824 iz=0,1
		    k=kb3+iz
		    cg3=kb3-chrgv2(3,ig)+1-iz
	        phimap(i,j,k)=phimap(i,j,k)+abs(cg1*cg2*cg3)*chrgv2(4,ig)
824	      continue
823	    continue
822	  continue
821	continue
c
c set up odd/even logical array 
	do 10 i=1,3*igrid,2
10	isum(i)=.true.
	do 20 i=2,(3*igrid-1),2
20	isum(i)=.false. 
c
c find which grid points have charge assigned to them
c (will use this array later to calculate grid energy)
c 
	n=0
	do 100 k=2,igrid-1
	  do 110 j=2,igrid-1
	    do 120 i=2,igrid-1
		if(phimap(i,j,k).ne.0) then
		n=n+1
		gchrgp(1,n)=i
		gchrgp(2,n)=j
		gchrgp(3,n)=k
		gchrg(n)=phimap(i,j,k)
		phimap(i,j,k)=0.0
		end if
120	 continue
110	 continue
100	 continue 
	icount1b=n
c
c determine how many charged grid points are odd
c 
	icount1a=0
	do 200 i=1,n
	itemp=gchrgp(1,i)+gchrgp(2,i)+gchrgp(3,i)
	if(isum(itemp)) icount1a=icount1a+1
200	continue 
c
c set up odd/even pointer array, to be used in making qval
c and iqpos
c 
	i1=0
	i2=icount1a
	do 300 i=1,n
	itemp=gchrgp(1,i)+gchrgp(2,i)+gchrgp(3,i)
	if(isum(itemp)) then
	i1=i1+1
	gchrg2(i)=i1
	else
	i2=i2+1
	gchrg2(i)=i2
	end if
300	continue 
c
c determine denominator at all charged grid points
c 
	ib=0
	epsins6=sixeps+ 6.0*difeps 
	do 400 i=1,n
	  iz=gchrgp(3,i)
	  iy=gchrgp(2,i)
	  ix=gchrgp(1,i)
	  itemp1=iepsmp(ix,iy,iz,1)+iepsmp(ix-1,iy,iz,1)
	  itemp2=iepsmp(ix,iy,iz,2)+iepsmp(ix,iy-1,iz,2)
	  itemp3=iepsmp(ix,iy,iz,3)+iepsmp(ix,iy,iz-1,3)
	  itemp=itemp1+itemp2+itemp3
	  gchrgd(i)=itemp*difeps + debfct*idebmap(ix,iy,iz) + sixeps 
	  if(itemp.ne.6) then
	    ib=ib+1
c	    cgbp(1,ib)=ix 
c	    cgbp(2,ib)=iy
c	    cgbp(3,ib)=iz 
	    cgbp(1,ib)=gchrg(i)*fpoh/(epsins6+debfct*idebmap(ix,iy,iz))
	    cgbp(2,ib)=gchrg2(i)
	  end if 
400	continue 
	ibc=ib 
	write(6,*) '# grid points charged and at boundary=',ib
c
c make qval, fpoh term so potentials will be in kt/e
c 
	do 500 i=1,n
	j=gchrg2(i)
	qval(j)=gchrg(i)*fpoh/gchrgd(i)
	gval(j)=gchrg(i)
500	continue 
c
c make iqpos
c 
	isgrid=igrid**2
	do 600 i=1,n
	j=gchrg2(i)
	ix=gchrgp(1,i)
	iy=gchrgp(2,i)
	iz=gchrgp(3,i)
	iw=1+ix+igrid*(iy-1)+isgrid*(iz-1)
	iv=iw/2
	iqpos(j)=iv
600	continue 
c
c end of chrgup, return with qval,iqpos and gchrgp and gchrg
c also icount1a, icount1b 
	return
	end 
