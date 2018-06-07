	subroutine setbc(xx,mgrid,
     &   phimap,phimap1,phimap2,phimap3,atmcrg,
     &   iepsmp,iepsmp2,idebmap,
     &   ibctyp,iper,qplus,qmin,cqplus,cqmin,
     &                 epsout,deblen,nqass)
c
c assigns either zero boundary conditions (ibctyp=1)
c quasi-coulombic based on debye dipole qplus at cqplus
c and qmin at cqmin (ibctyp=2)
c focussing (ibctyp=3)(continuance if scales are the same)
c full quasi-coulombic (ibctyp=4)
c or constant external filed (ibctyp=5)
c option 2 will be appropriately modified by periodic
c boundary condition flags (iper)
c--------------------------------------------------------------
	include 'qdiffpar4.h'
c--------------------------------------------------------------
	logical iper(3)
	dimension oldmid1(3)
	dimension cqplus(3),cqmin(3),g(3),go(3),c(3)
	character*20 toblbl
	character*10 label
	character*60 title
	character*80 filnam
	character*16 botlbl
c--------------------------------------------------------------
c
c zero option, clear boundary values
c
          do 9000 iz=1,igrid
                do 9002 ix=1,igrid,igrid-1
             do 9001 iy=1,igrid
                   phimap(ix,iy,iz) = 0.0
9001		 continue
9002		    continue
9000	    continue
          do 9003 iz=1,igrid
             do 9004 iy=1,igrid,igrid-1
                do 9005 ix=1,igrid
                   phimap(ix,iy,iz) = 0.0
9005		    continue
9004		 continue
9003	    continue
          do 9006 iz=1,igrid,igrid-1
             do 9007 iy=1,igrid
                do 9008 ix=1,igrid
                   phimap(ix,iy,iz) = 0.0
9008		    continue
9007		 continue
9006	    continue
c
c end of zero option
c
      if(ibctyp.eq.5) then
c
c constant field, 1 kT/e/grid unit, in the x direction
c
          do 9050 iz=1,igrid
             do 9051 iy=1,igrid
                do 9052 ix=1,igrid
                   phimap(ix,iy,iz) = ix
9052		    continue
9051		 continue
9050	    continue
	end if
c
c end external field option
c
      if(ibctyp.eq.2) then
c
c quasi coulombic dipole option
c
	do 9009 iz=1,igrid
	  do 9010 iy=1,igrid
            dist1 = (cqplus(1)-1)**2 + (cqplus(2)-iy)**2
     &	     + (cqplus(3)-iz)**2
 	      dist1 = sqrt(dist1)/scale
 	      tempp = qplus*exp(-dist1/deblen )/(dist1*epsout) 
            dist2 = (cqmin(1)-1)**2 + (cqmin(2)-iy)**2
     &	     + (cqmin(3)-iz)**2
 	      dist2 = sqrt(dist2)/scale
 	      tempn = qmin*exp(-dist2/deblen )/(dist2*epsout) 
            phimap(1,iy,iz) = phimap(1,iy,iz) + tempp + tempn
            dist3 = (cqplus(1)-igrid)**2 + (cqplus(2)-iy)**2
     &	     + (cqplus(3)-iz)**2
 	      dist3 = sqrt(dist3)/scale
 	      tempp = qplus*exp(-dist3/deblen )/(dist3*epsout) 
            dist4 = (cqmin(1)-igrid)**2 + (cqmin(2)-iy)**2
     &	     + (cqmin(3)-iz)**2
 	      dist4 = sqrt(dist4)/scale
 	      tempn = qmin*exp(-dist4/deblen )/(dist4*epsout) 
            phimap(igrid,iy,iz) = phimap(igrid,iy,iz) + tempp + tempn
9010		 continue
9009	    continue
	do 9012 iz=1,igrid
	  do 9013 iy=1,igrid,igrid-1
	    do 9014 ix=1,igrid
            dist = (cqplus(1)-ix)**2 + (cqplus(2)-iy)**2
     &	     + (cqplus(3)-iz)**2
 	      dist = sqrt(dist)/scale
 	      tempp = qplus*exp(-dist/deblen )/(dist*epsout) 
            dist = (cqmin(1)-ix)**2 + (cqmin(2)-iy)**2
     &	     + (cqmin(3)-iz)**2
 	      dist = sqrt(dist)/scale
 	      tempn = qmin*exp(-dist/deblen )/(dist*epsout) 
            phimap(ix,iy,iz) = phimap(ix,iy,iz) + tempp + tempn
9014		    continue
9013		 continue
9012	    continue
	do 9015 iz=1,igrid,igrid-1
	  do 9016 iy=1,igrid
	    do 9017 ix=1,igrid
            dist = (cqplus(1)-ix)**2 + (cqplus(2)-iy)**2
     &	     + (cqplus(3)-iz)**2
 	      dist = sqrt(dist)/scale
 	      tempp = qplus*exp(-dist/deblen )/(dist*epsout) 
            dist = (cqmin(1)-ix)**2 + (cqmin(2)-iy)**2
     &	     + (cqmin(3)-iz)**2
 	      dist = sqrt(dist)/scale
 	      tempn = qmin*exp(-dist/deblen )/(dist*epsout) 
            phimap(ix,iy,iz) = phimap(ix,iy,iz) + tempp + tempn
9017		    continue
9016		 continue
9015	    continue
      end if 
c
c end of quasi coulombic dipole option
c
      if(ibctyp.eq.4) then
c
c a summation of the potential resulted from each point of charge 
c
	do 9034 iz=1,igrid
	  do 9035 iy=1,igrid
	    do 9036 ix=1,igrid,igrid-1
		do 9037 ic = 1,nqass
              dist = (atmcrg(1,ic)-ix)**2 + (atmcrg(2,ic)-iy)**2
     &	     + (atmcrg(3,ic)-iz)**2
 	        dist = sqrt(dist)/scale
	        tempd = atmcrg(4,ic)*exp(-dist/deblen)/(dist*epsout)
              phimap(ix,iy,iz) = phimap(ix,iy,iz) + tempd 
9037		continue
9036	    continue
9035    continue
9034  continue
	do 9038 iz=1,igrid
	  do 9039 iy=1,igrid,igrid-1
	    do 9040 ix=1,igrid
		do 9041 ic = 1,nqass
              dist = (atmcrg(1,ic)-ix)**2 + (atmcrg(2,ic)-iy)**2
     &	     + (atmcrg(3,ic)-iz)**2
 	        dist = sqrt(dist)/scale
	        tempd = atmcrg(4,ic)*exp(-dist/deblen)/(dist*epsout)
              phimap(ix,iy,iz) = phimap(ix,iy,iz) + tempd 
9041		continue
9040	    continue
9039    continue
9038  continue
	do 9042 iz=1,igrid,igrid-1
	  do 9043 iy=1,igrid
	    do 9044 ix=1,igrid
		do 9045 ic = 1,nqass
              dist = (atmcrg(1,ic)-ix)**2 + (atmcrg(2,ic)-iy)**2
     &	     + (atmcrg(3,ic)-iz)**2
 	        dist = sqrt(dist)/scale
	        tempd = atmcrg(4,ic)*exp(-dist/deblen)/(dist*epsout)
              phimap(ix,iy,iz) = phimap(ix,iy,iz) + tempd 
9045		continue
9044	    continue
9043    continue
9042  continue
      end if 
c
c end of the option for the complete summation of potential
c
c
      if(ibctyp.eq.3) then
c 
c focussing option-bc's come from a previous phimap
c
	open(unit=18,status='old',
     1err=900,form='unformatted')
	filnam = ' '
	inquire(18,name = filnam)
	write(6,*)' '
	write(6,*)'focussing boundary condition '
	write(6,*)'read from file'
	write(6,*)filnam
	write(6,*)' '
c
c read in old potential map
c
	read(18)toplbl
      read(18)label,title
      read(18)(((phimap(i,j,k),i=1,mgrid),j=1,mgrid),k=1,mgrid)
      read(18)botlbl
	read(18)scale1,oldmid1
	close(18)
c check to see if this is a continuence
	
	if(scale1.eq.scale) then
	write(6,*) 'scales are the same.' 
	write(6,*) 'therefore assuming this to be a continuence'
	goto 511
	endif
	write(6,*)' '
      write(6,*)' focussing potential map:'
	write(6,200)title
	write(6,*)'original scale (grids/A)      : ',scale1
	write(6,*)'object centre at (A) : ',oldmid1
	write(6,*)' '
200	format(A60)
c
c check to see that new grid lies within old one that is going to
c provide bc's
c
      iout = 0
	goff = (igrid+1.)/2.
	do 9018 iz=1,igrid,igrid-1
        g(3) = iz
	  do 9019 iy=1,igrid,igrid-1
          g(2) = iy
	    do 9020 ix=1,igrid,igrid-1
            g(1) = ix
c
c for each new grid corner, calculate old grid coords
c
            call gtoc(mgrid,oldmid,scale,g,c)
            do 9021 i = 1,3
		  gold = (c(i) - oldmid1(i))*scale1 + goff
              if((gold.le.2.).or.(gold.ge.igrid-1.))  iout = 1
9021		continue
9020	    continue
9019	  continue
9018	continue
      if(iout.ne.0) then
	  write(6,*)'part of new grid lies outside old grid'
	  write(6,*)'check scaling of both grids'
	  write(6,*)'old scale:'
	  write(6,*)'scale (grids/A)      : ',scale1
	  write(6,*)'object centre at (A) : ',oldmid1
	  write(6,*)'new scale:'
	  write(6,*)'scale (grids/A)      : ',scale
	  write(6,*)'object centre at (A) : ',oldmid
	  stop
      end if
c
c for each boundary point
c convert to real coordinates
c convert to old grid coordinates
c interpolate potential
c note that can use same potential array for boundaries
c since old potentials at boundary are not used for new ones
c
c
c save new grid size, and set temporarily to 65
c
	isgrid = igrid
c	igrid = ngrid
	gmid = (isgrid + 1.)/2.
	write(6,*)igrid
      write(6,*)'pulling boundary values out of old potential map...'
	do 9022 iz=2,isgrid-1
        g(3) = iz
	  do 9023 iy=2,isgrid-1
          g(2) = iy
	    do 9024 ix=1,isgrid,isgrid-1
            g(1) = ix
c
c for each new grid side, calculate old grid coords
c
		do 9025 i = 1,3
		  c(i) = (g(i) - gmid)/scale + oldmid(i)
		  go(i) = (c(i) - oldmid1(i))*scale1 + goff

9025		continue
c
c find potential
c
            call phintp(xx,mgrid,phimap,phimap1,phimap2,phimap3,atmcrg,
     &                  iepsmp,iepsmp2,idebmap,go,phiv)
            phimap(ix,iy,iz) = phiv

9024		    continue
9023		 continue
9022	    continue

	do 9026 iz=2,isgrid-1
        g(3) = iz
	  do 9027 iy=1,isgrid,isgrid-1
          g(2) = iy
	    do 9028 ix=2,isgrid-1
            g(1) = ix
c
c for each new grid side, calculate old grid coords
c
		do 9029 i = 1,3
		  c(i) = (g(i) - gmid)/scale + oldmid(i)
		  go(i) = (c(i) - oldmid1(i))*scale1 + goff
9029		continue
c
c find potential
c
            call phintp(xx,mgrid,phimap,phimap1,phimap2,phimap3,atmcrg,
     &                  iepsmp,iepsmp2,idebmap,go,phiv)
            phimap(ix,iy,iz) = phiv

9028		    continue
9027		 continue
9026	    continue

	do 9030 iz=1,isgrid,isgrid-1
        g(3) = iz
	  do 9031 iy=2,isgrid-1
          g(2) = iy
	    do 9032 ix=2,isgrid-1
            g(1) = ix
c
c for each new grid side, calculate old grid coords
c
		do 9033 i = 1,3
		  c(i) = (g(i) - gmid)/scale + oldmid(i)
		  go(i) = (c(i) - oldmid1(i))*scale1 + goff
9033		continue
c
c find potential
c
            call phintp(xx,mgrid,phimap,phimap1,phimap2,phimap3,atmcrg,
     &                  iepsmp,iepsmp2,idebmap,go,phiv)
            phimap(ix,iy,iz) = phiv
c
9032		    continue
9031		 continue
9030	    continue
c restore new grid size
c
	igrid = isgrid
511   end if 
c
c end of focussing option
c

	return
900	write(6,*)' no potl map for focussing boundary conditions'
   	stop
	end
