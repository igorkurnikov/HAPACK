	subroutine mkeps(xx,mgrid,
     & phimap,phimap1,phimap2,phimap3,atmcrg,	
     & iepsmp,iepsmp2,idebmap,
     & radprb,epsout,ibnum,ibgrd,rad3,xn2)
c
c fills up intersticies between  vanderWaals atom spheres
c that would otherwise be given a dielectric constant of water:
c if probe radius, radprb, not zero, then finds
c every inside grid point that has any outside eps
c points as neighbours, (ie all boundary points)
c and turns all eps points within
c probe radius back to outside dielectric- this generates an eps map
c that approximates the solvent accessible surface
c
	include 'qdiffpar4.h'
	dimension ibgrd(4,*),xn2(3,natmax),rad3(natmax),mv(6)
c
c----------------------------------------------------------------
c----------------------------------------------------------------
c
	if(radprb.eq.0.0) return
	write(6,*)' finishing off dielectric map...'
c
	ix1=igrid
	iy1=igrid
	iz1=igrid
	ix65=1
	iy65=1
	iz65=1
	do 543 i=1,ibnum
	ib1=ibgrd(1,i)
	ib2=ibgrd(2,i)
	ib3=ibgrd(3,i)
	ix1=min(ib1,ix1)
	iy1=min(ib2,iy1)
	iz1=min(ib3,iz1)
	ix65=max(ib1,ix65)
	iy65=max(ib2,iy65)
	iz65=max(ib3,iz65)
543	continue
	minp=min(ix1,iy1,iz1)
	maxp=max(ix65,iy65,iz65)
	write(6,*) "expanded surface maximum, minimum (g.u.)=",minp,maxp
c
c
	radp = radprb*scale
	radp2 = radp**2
	radp4= (radp-0.5)**2
	rgrid = igrid
	lim=int(radp+1.5)
c
	rad2a = radp2 - 0.25
	ivz=0
c
	do 50 i=1,ibnum
c
	ib1=ibgrd(1,i)
	ib2=ibgrd(2,i)
	ib3=ibgrd(3,i)
	xn=float(ib1)
	yn=float(ib2)
	zn=float(ib3)
c
	mv(1)=iepsmp2(ib1,ib2,ib3,1)
	mv(2)=iepsmp2(ib1,ib2,ib3,2)
	mv(3)=iepsmp2(ib1,ib2,ib3,3)
	mv(4)=iepsmp2(ib1-1,ib2,ib3,1)
	mv(5)=iepsmp2(ib1,ib2-1,ib3,2)
	mv(6)=iepsmp2(ib1,ib2,ib3-1,3)
	dism=1000.
	iat=100000
	do j=1,6
	  m=mv(j)
	  if((m.ne.0).and.(m.ne.iat)) then
	    dist=(xn2(1,m)-xn)**2 + (xn2(2,m)-yn)**2 + (xn2(3,m)-zn)**2
	    if(dist.lt.dism) then
	      iat=m
	      dism=dist
	    end if
	  end if
	end do
	iv=iat
c
	if(iv.eq.0) then
	ivz=ivz+1
	else
c
	xa=xn2(1,iv)
	ya=xn2(2,iv)
	za=xn2(3,iv)
c
	arad=(xn-xa)**2 + (yn-ya)**2 + (zn-za)**2
	rad=(rad3(iv)*scale+radp)
	sfact=rad/(sqrt(arad))
c
	xn=xa + sfact*(xn-xa)
	yn=ya + sfact*(yn-ya)
	zn=za + sfact*(zn-za)
	ib1=int(xn)
	ib2=int(yn)
	ib3=int(zn)
	end if
c
c set loop limits 
c
	lim1=max0(ib1-lim,1)
	lim2=min0(ib1+lim,igrid)
	lim3=max0(ib2-lim,1)
	lim4=min0(ib2+lim,igrid)
	lim5=max0(ib3-lim,1)
	lim6=min0(ib3+lim,igrid)
c
	  z=float(lim5)-1.0
	  y1=float(lim3)-1.0
	  x1=float(lim1)-1.0
c
          do 9019 iz = lim5,lim6 
	  z=z+1.0
	  y=y1
            do 9020 iy = lim3,lim4 
	    y=y+1.0
	    x=x1
              do 9021 ix = lim1,lim2 
	      x=x+1.0
                  dxyz1 = (x - xn)
                  dxyz2 = (y - yn)
                  dxyz3 = (z - zn)
                  distsq = dxyz1**2 +dxyz2**2 +dxyz3**2 
		  distsq1 = distsq + dxyz1
		  distsq2 = distsq + dxyz2
		  distsq3 = distsq + dxyz3
		if(distsq1.lt.rad2a) iepsmp(ix,iy,iz,1)=0
		if(distsq2.lt.rad2a) iepsmp(ix,iy,iz,2)=0
		if(distsq3.lt.rad2a) iepsmp(ix,iy,iz,3)=0
9021		continue
9020	    continue
9019	  continue
50	continue
	if(ivz.ne.0) write(6,*) "no. of surface points unassigned= ",ivz
	return
	end
