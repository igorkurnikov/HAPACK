	subroutine mem(xx,mgrid,phimap,phimap1,phimap2,phimap3,atmcrg,
     &               iepsmp,iepsmp2,idebmap,islice)
	include 'qdiffpar4.h'
	dimension off(3),goff(3),islice(2,*)
	write(6,*) ' '
	write(6,*) 'membrane layer being included...'
	write(6,*) ' ' 
	open(25,file='mem.dat')
	read(25,*) off1,off2,off3
	read(25,*) zext,hole
	off(1)=off1
	off(2)=off2
	off(3)=off3
	call ctog(mgrid,oldmid,scale,off,goff)
	off1=goff(1)
	off2=goff(2)
	off3=goff(3)
	ext=zext*scale
	iext1=off3-ext
	iext2=off3+ext 
c
c rmmin,rmmax are the upper and lower boundaries of the membrane in grid units
c they are global parameters
c
	rmmin=off3-ext
	rmmax=off3+ext
	write(6,*) "membrane runs from ",rmmin," to ",rmmax
	if(iext1.lt.1) then
	write(6,*) 'truncating membranes lower level to the z=1 plane'
	iext1=1
	end if
	if(iext2.gt.igrid) then
	write(6,*) 'truncating membranes upper level to the z=',igrid,'plane'
	iext2=igrid
	end if 
	hole=hole*scale 
	hole2=hole**2 
	n=0
	do 10 x=1,igrid
	do 20 y=1,igrid
	temp=(x-off1)**2 + (y-off2)**2
	if(temp.gt.hole2) then
	n=n+1
	islice(1,n)=x
	islice(2,n)=y
	end if
20	continue
10	continue 
	write(6,*) 'number of x-y grid points affected is',n 
	write(6,*) 'from z=',iext1,'    to     ',iext2 
	do 30 i=1,n
	ix=islice(1,i)
	iy=islice(2,i)
	do 40 iz=iext1,iext2
	iepsmp(ix,iy,iz,1)=1
	iepsmp(ix,iy,iz,2)=1
	iepsmp(ix,iy,iz,3)=1
	idebmap(ix,iy,iz)=0
40	continue
30	continue 
	do 31 i=1,n
	ix=islice(1,i)
	iy=islice(2,i)
	do 41 iz=iext1,(iext2+iext1)/2
	iepsmp2(ix,iy,iz,1)=-2
	iepsmp2(ix,iy,iz,2)=-2
	iepsmp2(ix,iy,iz,3)=-2
41	continue
	do 42 iz=(iext2+iext1+2)/2,iext2
	iepsmp2(ix,iy,iz,1)=-1
	iepsmp2(ix,iy,iz,2)=-1
	iepsmp2(ix,iy,iz,3)=-1
42	continue
31	continue 
	return
	end 
