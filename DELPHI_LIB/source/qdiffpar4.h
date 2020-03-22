c---------------------------------------------------
c parameters for qdiff4.f and subroutines 
c---------------------------------------------------
	   dimension xx(*)

	parameter (ncrgmx = 50000)
c ncrgmx = maximum number of charged atoms 
c	parameter (ngrid = 65)
	parameter (ngrmax = 1000)
c ngrid= maximum grid size 
c	parameter (ngp = (ngrid*ngrid*ngrid+1))
c ngp = ngrid**3 = number of grid points 
c	parameter (nhgp = ngp/2 )
c nhgp = half the previous number 
	parameter (natmax = 50000)
c natmax = maximum number of atoms 
	parameter (ngcrg =  200000)
c ngcrg = maximum number of grid points that may be assigned charge (was 4226)
c	parameter (nbgp = ngrmax*ngrmax*2 +1 )
c nbgp = number of points at one box boundary = ngrid**2
c	parameter (nsp = 5*nbgp)
c nsp = maximum number of dielectric boundary points, appox.= 5*nbgp 
c---------------------------------------------------
c	dimension atmcrg(4,ncrgmx)	
c	atom charges for boundary condition
c 	dimension phimap(ngrid,ngrid,ngrid),phimap1(nhgp),phimap2(nhgp)
c	dimension phimap3(ngp)
	dimension phimap(mgrid,mgrid,*),phimap1(*),phimap2(*)
	dimension phimap3(*)
	dimension atmcrg(4,*)
 	integer iepsmp(mgrid,mgrid,mgrid,3),iepsmp2(mgrid,mgrid,mgrid,3)
	integer idebmap(mgrid,mgrid,mgrid)
	dimension cgbp(2,ngcrg),gval(ngcrg)
	dimension oldmid(3)
c---------------------------------------------------
c---------------------------------------------------
      character*1 chn
      character*3 res
      character*4 rnum
      character*6 atm
c---------------------------------------------------
c---------------------------------------------------
	common
     &	/scaleq/ scale,oldmid,igrid,ibc,gten,cgbp,gval,rmmin,rmmax
     &  /hadel1/ iphi,iphi1,iphi2,iphi3
c---------------------------------------------------
c rmmin,rmmax are the upper and lower boundaries of the membrane in grid units
c---------------------------------------------------





