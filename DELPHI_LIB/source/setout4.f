	subroutine setout(xx,mgrid,
     &  phimap,phimap1,phimap2,phimap3,atmcrg,
     &  iepsmp,iepsmp2,idebmap,
     &  ioff,xn2,rad3,natom,exrad,radprb)
	include 'qdiffpar4.h'
	dimension xn2(3,natmax),rad3(natmax),xn(3)
	integer ioff(3,*),ismin(3),ismax(3)
	dimension sq1(-15:15),sq2(-15:15),sq3(-15:15)
	dimension rad2a1(-15:15),rad2a2(-15:15),rad2a3(-15:15)
	logical itobig,itest2
c
	radmax2=0.0
	do 691 ix=1,natom
	radmax2=amax1(radmax2,rad3(ix))
691	continue
	temp=amax1(exrad,radprb)
	radmax2=scale*(radmax2+temp)
	lim= 1+ radmax2
c 
	limmax = 12
	if(lim.gt.limmax) itobig=.true.
c 
	if(itobig) goto 7878 
	radtest= (radmax2 + 0.5*sqrt(3.0))**2
	ibox=0
	do 692 ix=-lim,lim
	  do 693 iy=-lim,lim
	    do 694 iz=-lim,lim
	    dist=ix**2 + iy**2 + iz**2 
	    dist1=dist + ix + 0.25
	    dist2=dist + iy + 0.25
	    dist3=dist + iz + 0.25
	    itest=0
	    if(dist.lt.radtest) itest=1
	    if(dist1.lt.radtest) itest=1
	    if(dist2.lt.radtest) itest=1
	    if(dist3.lt.radtest) itest=1 
	    if(itest.eq.1) then 
	    ibox=ibox+1
	    ioff(1,ibox)=ix
	    ioff(2,ibox)=iy
	    ioff(3,ibox)=iz
	    end if
694	continue
693	continue
692	continue
7878	continue 
c
c 
c set interiors
c
c read data file
c
	do 607 iv=1, natom
c
c restore values
c
	rad= rad3(iv)
	xn(1)=xn2(1,iv)
	xn(2)=xn2(2,iv)
	xn(3)=xn2(3,iv)
	if(rad.eq.0.) goto 608
c
c scale radius to grid
c
	  iv10=10*iv+1
	  rad = rad*scale
	  rad5= (rad + 0.5)**2
	  radp = rad + exrad*scale
	  rad = rad + radprb*scale
	  rad4= (rad + 0.5)**2
	  rad2 = rad*rad
	  radp2 = radp*radp
c
c set dielectric map
c
c find upper and lower grid index limits
c that fully enclose sphere of atom
c ensure they lie within grid
c
c
	itest2=.false.
        do 9017 k = 1,3
          ismin(k) = (xn(k) - radmax2 - 1.)
	  itest1=ismin(k)
	    ismin(k) = min(ismin(k),igrid)
	    ismin(k) = max(ismin(k),1)
	    if(itest1.ne.ismin(k)) itest2=.true.
          ismax(k) = (xn(k) + radmax2 + 1.)
	  itest1=ismax(k)
	    ismax(k) = min(ismax(k),igrid)
	    ismax(k) = max(ismax(k),1)
	    if(itest1.ne.ismax(k)) itest2=.true.
9017	continue
c
c
	if(itest2.or.itobig) then
	num=num+1
	rad2a = rad2 - 0.25
          do 9019 iz =  ismin(3),ismax(3)
            do 9020 iy =  ismin(2),ismax(2)
              do 9021 ix =  ismin(1),ismax(1)
                  dxyz1 = (ix - xn(1))
                  dxyz2 = (iy - xn(2))
                  dxyz3 = (iz - xn(3))
                  distsq = dxyz1**2 +dxyz2**2 +dxyz3**2 
		  distsq1 = distsq + dxyz1
		  distsq2 = distsq + dxyz2
		  distsq3 = distsq + dxyz3
		if(distsq1.lt.rad2a) iepsmp2(ix,iy,iz,1)=iv
		if(distsq2.lt.rad2a) iepsmp2(ix,iy,iz,2)=iv
		if(distsq3.lt.rad2a) iepsmp2(ix,iy,iz,3)=iv
		if(distsq.lt.radp2)  idebmap(ix,iy,iz) =0
9021		continue
9020	    continue
9019	  continue
		else
	rad2a = rad2 - 0.25
	ixn1=nint(xn(1))
	iyn1=nint(xn(2))
	izn1=nint(xn(3))
	fxn1=ixn1-xn(1)
	fxn2=iyn1-xn(2)
	fxn3=izn1-xn(3)
	rad2ax=rad2a-fxn1
	rad2ay=rad2a-fxn2
	rad2az=rad2a-fxn3
	do 6020 ix=-lim,lim
	temp1=ix+fxn1
	temp2=ix+fxn2
	temp3=ix+fxn3
	sq1(ix)=temp1*temp1
	sq2(ix)=temp2*temp2
	sq3(ix)=temp3*temp3
	rad2a1(ix)=rad2a-temp1
	rad2a2(ix)=rad2a-temp2
	rad2a3(ix)=rad2a-temp3
6020	continue 
C$DIR NO_RECURRENCE
		  do 9024 i=1,ibox
		  i1= ioff(1,i)
		  i2= ioff(2,i)
		  i3= ioff(3,i)
		  ix=ixn1+ i1
		  iy=iyn1+ i2
		  iz=izn1+ i3
        distsq = sq1(i1) +sq2(i2) + sq3(i3)
	if(distsq.lt.rad2a1(i1)) iepsmp2(ix,iy,iz,1)=iv
	if(distsq.lt.rad2a2(i2)) iepsmp2(ix,iy,iz,2)=iv
	if(distsq.lt.rad2a3(i3)) iepsmp2(ix,iy,iz,3)=iv
        if(distsq.lt.radp2)   idebmap(ix,iy,iz)=0
9024	continue
		end if
c
608	continue
c
607	continue
c
	return
	end
