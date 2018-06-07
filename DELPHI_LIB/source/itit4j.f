	subroutine itit(xx,mgrid,
     &  phimap,phimap1,phimap2,phimap3,atmcrg,
     &  iepsmp,iepsmp2,idebmap,
     &  nlit,nnit,iper,idpos,db,
     &  bndx1,bndx2,bndx3,bndx4,bndx,bndy,bndz,
     &  sf1,sf2,iqpos,qval,icount2a,icount2b,icount1a,icount1b,
     &  rionst,spec,icon1,icon2,epsin,epsout)
c
c
c NOTE THIS HAS BEEN ALTERED TO GIVE THE CONVERGENCE TO THE EXACT
C SOLUTION, WHICH IS INPUTTED VIA OPTION 3,(AND NOW WITH GAUSS
C SIEDEL IMPLEMENTED, WITH MAPPED SUBARRAYS.)
c
c  finally what we`ve all been waiting for-
c do the actual pb equation iteration
c first the linear, with relaxation parameter 1.
c then non-linear with relaxation parameter 0.6
c
c this is a modified iteration routine which runs at about
c three times that of the original version. this has been 
c accomplished by several features which will be mentioned
c at the appropiate juncture.
c note that the arrays slab,old and denom no longer exist
c and that epsmap is only used in setting up arrays
c in their place are several arrays, the two biggest being
c phimap2 and sf.
c however, the array space accessed in the main inner loop
c is no larger than before, i.e. memory requirements should
c not be much different.
c
c some notes on array sizes. there can be no more than 10000
c charges, or 12000 boundary elements, or 40000 high potential salt
c sites (for the nonlinear). also, for the latter the potential
c in salt should not exceed 10kt
c
	include 'qdiffpar4.h'
c 
	parameter ( nxran = 60 )
	parameter ( nyran = 20 )
	dimension rmsl(nxran),rmsn(nxran),rmaxl(nxran),rmaxn(nxran)
	dimension qval(ngcrg)
	dimension db(6,*),sf1(*),sf2(*)
	integer*4 iqpos(ngcrg),idpos(*)
	dimension bndx1(*),bndx2(*),bndx3(*),bndx4(*)
	dimension bndx(*),bndy(*),bndz(*),grdn(1000)
	character*8 hour
	character*1 symb
	character*70 title,iplot(nyran)
	logical iper(3),qstopper,resdat,igt
	integer*4 star,fin,sta1(ngrmax),sta2(ngrmax),fi1(ngrmax),fi2(ngrmax)
c

c
c some initialization
c
	icount2abac=icount2a
	icount2bbac=icount2b 
	epsrat=epsout/epsin 
	sixth = 1./6.
	th120 = 1./120.
	icgrid=igrid**3
	ihgd=(igrid+1)/2
	if(icon2.eq.0) then
	icon1=10
	icon2=1
	end if 
c	if(icon1.gt.nlit .and. nlit.gt.1) icon1=nlit 
c
C
	do 9000 i = 1,nxran
	  rmsl(i) = 0.0
	  rmsn(i) = 0.0
	  rmaxl(i) = 0.0
	  rmaxn(i) = 0.0
9000	continue
      npoint = (igrid-2)**3
c
c
c comment out for cray version 14, nov 88, kas
c
c ---------------------------------------------
c MAIN SET UP ROUTINE	
c ---------------------------------------------
c
	if(iper(1)) then
	n=0
	do 440 iz=2,igrid-1
	  iadd1=(iz-1)*igrid*igrid 
	  do 441 iy=2,igrid-1
	    iadd2=(iadd1+(iy-1)*igrid +2)/2
	    n=n+1
	    bndx(n)=iadd2
441	  continue
440     continue
	idif1x=(igrid-2)/2
	idif2x=idif1x+1
	inc1xa=1
	inc1xb=0
	inc2xa=0
	inc2xb=1
	end if 
	if(iper(2)) then
	n=0
	do 442 iz=2,igrid-1
	  iadd1=(iz-1)*igrid*igrid 
	  do 443 ix=2,igrid-1
	    iadd2=(iadd1+ix+1)/2
	    n=n+1
    	    bndy(n)=iadd2
443       continue
442     continue
	idif1y=igrid*(igrid-2)/2
	idif2y=idif1y+1
	inc1ya=(igrid/2)+1
	inc1yb=inc1ya-1
	inc2ya=inc1yb
	inc2yb=inc1ya
	end if
	if(iper(3)) then
	n=0
	do 444 ix=2,igrid-1
	  iadd1=ix+1
	  do  445 iy=2,igrid-1
	    iadd2=(iadd1+(iy-1)*igrid)/2
	    n=n+1
	    bndz(n)=iadd2
445	  continue
444	continue
	idif1z=igrid*igrid*(igrid-2)/2
	idif2z=idif1z+1
	inc1za=((igrid**2)/2)+1
	inc1zb=inc1za
	inc2za=inc1zb
	inc2zb=inc1za
	end if
c
c
c END OF SET UP	
c
	resdat = .false.
      res1=0.0
      res2=0.0
	res3=0.0
	res4=0.0
c      end if
	if(resdat) then 
	write(6,*) ' '
	write(6,*) 'linear resolution criteria are:',res1,res2
	if(nnit.ne.0) then
	write(6,*) 'non-linear resolution criteria are:',res3,res4
	end if
	write(6,*) ' '
	end if 
c
c
	write(6,*) ' '
	write(6,*) ' '
	if(gten.eq.0.)
     &        write(6,*) '  rms-change     max change       #iterations'
	if(gten.ne.0.)write(6,*) '  rms-change     max change  grid energy 
     &                #iterations'
c
c
c set up start and stop vectors
	sta1(2)=(igrid**2 + igrid +4)/2
	sta2(2)=sta1(2)-1
	fi1(2)=igrid**2 - (igrid+1)/2
	fi2(2)=fi1(2)
	itemp1=igrid + 2
	itemp2=igrid**2 -igrid -2
	do 225 i=3,igrid-1
    	  sta1(i)=fi1(i-1) + itemp1
	  sta2(i)=fi2(i-1) + itemp1
	  fi1(i)=sta1(i-1) + itemp2
	  fi2(i)=sta2(i-1) + itemp2
225   continue
c
c also
c
	lat1= (igrid-1)/2
	lat2= (igrid+1)/2
	long1= (igrid**2 - 1)/2
	long2= (igrid**2 + 1)/2
c 
c
      ires=0
c
c
	om2=2.0/(1.0 + sqrt(1 - spec))
	om1=1.0-om2
	do 801 ix=1,(icgrid+1)/2
	  sf1(ix)=sf1(ix)*om2
	  sf2(ix)=sf2(ix)*om2
801	continue
	do 802 ix=1,icount1b 
	  qval(ix)=qval(ix)*om2
802	continue 
	do 803 iy=1,6
	  do 804 ix=1,icount2b
	  db(iy,ix)=db(iy,ix)*om2
804	continue
803	continue 
	sixth=sixth*om2 
	i=1
c
      iw=1
	do 311 iz=1,igrid
	  do 312 iy=1,igrid
	    do 313 ix=1,igrid
	      phimap3(iw)=phimap(ix,iy,iz)
	      iw=iw+1
313       continue
312     continue
311   continue
c
c 
1689  do 314 ix=1,(icgrid+1)/2
	  iy=ix*2
	  phimap1(ix)=phimap3(iy-1)
	  phimap2(ix)=phimap3(iy)
314	continue
c
c       set boundary values 
c 
	star=(igrid+1)/2
	iy=(igrid*(igrid+1)/2) - igrid + 1
	fin=(igrid*(igrid-1)-2)/2
	ihgd2=ihgd-1
	do 9211 ix=star,fin
	  iy=iy+igrid
	  bndx1(ix)=phimap1(iy)
	  bndx2(ix)=phimap1(iy+ihgd2)
9211  continue
c 
	star=(igrid+2)/2
	iy=(igrid*(igrid+2)/2) - igrid +1
	fin=(igrid*(igrid-1)-1)/2
	do 8211 ix=star,fin
     	  iy=iy+igrid
	  bndx3(ix)=phimap2(iy)
	  bndx4(ix)=phimap2(iy+ihgd2)
8211  continue
c

1000      continue
c
c clear rms, max change
c
	  rmsch = 0.0
	  rmxch = 0.00

c
c if there is no salt then the main loop is executed without sf
c saving about 15% in execution time
c
c 
c
	if(rionst.ne.0.0) then
        do 9004 n = 2, igrid-1
	      star=sta1(n)
	      fin=fi1(n)
            do 9006 ix = star,fin
              temp1 = phimap2(ix) +
     &        phimap2(ix-1)
c
              temp2 = phimap2(ix+lat1) +
     &        phimap2(ix-lat2)
c
              temp3 = phimap2(ix+long1) +
     &        phimap2(ix-long2)
c
       	    phimap1(ix) =phimap1(ix)*om1 + (temp1+temp2+temp3)*sf1(ix)
9006	     continue
9004	  continue
c
c otherwise the main loop is as below:
c
      else
        do 9104 n = 2, igrid-1
	    star=sta1(n)
	    fin=fi1(n)
            do 9106 ix = star,fin
              temp1 = phimap2(ix) +
     &        phimap2(ix-1)
c
              temp2 = phimap2(ix+lat1) +
     &        phimap2(ix-lat2)
c
              temp3 = phimap2(ix+long1) +
     &        phimap2(ix-long2)
c
       	    phimap1(ix) = phimap1(ix)*om1 + (temp1+temp2+temp3)*sixth
9106	    continue
9104	  continue
      end if
c
c the above loops are about fourtimes faster than the original
c loop over all grid points for several reasons, the biggest being that
c we are only solving laplace's equation (unless salt is present), which
c numerically much simpler, hence faster. we put all we leave out, back
c in below, ending up with an equivalent calculation, but much faster.
c
c first we add back the dielectric boundary points, by recalculating them
c individually. note this is still vectorised by means of a gathering
c load by the compiler.
c
C$DIR NO_RECURRENCE 
      do 9010 k=1,icount2a
	  ix=idpos(k)
	  temp1=phimap2(ix-1)*db(1,k)+phimap2(ix)*db(2,k)
	  temp2=phimap2(ix-lat2)*db(3,k)+phimap2(ix+lat1)*db(4,k)
        temp3=phimap2(ix-long2)*db(5,k)+phimap2(ix+long1)*db(6,k)
        phimap1(ix)= phimap1(ix) + temp1+temp2+temp3
9010  continue
c
c next we add back an adjustment to all the charged grid points due to
c the charge assigned. the compiler directive just reassures the vector
c compiler that all is well as far as recurrence is concerned, i.e. it
c would think there is a recurrence below, where as in fact there is none.
c
c
c Now reset boundary values altered in above loops.
c
c 
	star=(igrid+1)/2
	iy=(igrid*(igrid+1)/2) - igrid +1
	fin=(igrid*(igrid-1)-2)/2
C$DIR NO_RECURRENCE
	do 9201 ix=star,fin
	  iy=iy+igrid
	  phimap1(iy)=bndx1(ix)
	  phimap1(iy+ihgd2)=bndx2(ix)
9201  continue
c
C$DIR NO_RECURRENCE 
	do 9011 k=1,icount1a
	  temp=qval(k)
	  ix=iqpos(k)
	  phimap1(ix)=phimap1(ix) + temp
9011  continue
c
c
c
c if periodic boundary condition option
c force periodicity using wrap around update of boundary values:
c 2nd slice-->last
c last-1 slice-->first
c
c z periodicity
c
      if(iper(3)) then
          do 9013 iz = 1,(igrid-1)**2,2
	      temp1=bndz(iz)
	      temp2=temp1+idif1z
	      temp3=temp2+inc1za
	      temp4=temp1+inc1zb
	      phimap1(temp1)=phimap2(temp2)
	      phimap1(temp3)=phimap2(temp4)
9013      continue
      end if
c
c y periodicity
c
      if(iper(2)) then
        do 9015 iy = 1,(igrid-1)**2,2
	      temp1=bndy(iy)
	      temp2=temp1+idif1y
	      temp3=temp2+inc1ya
	      temp4=temp1+inc1yb
	      phimap1(temp1)=phimap2(temp2)
	      phimap1(temp3)=phimap2(temp4)
9015	  continue
      end if
c
c x periodicity
c
      if(iper(1)) then
          do 9017 ix = 1,(igrid-1)**2,2
	      temp1=bndx(ix)
	      temp2=temp1+idif1x
	      temp3=temp2+inc1xa
	      temp4=temp1+inc1xb
	      phimap1(temp1)=phimap2(temp2)
	      phimap1(temp3)=phimap2(temp4)
9017      continue
      end if
c
c
c Next update phimap2 using the new phimap1
c
      if(rionst.ne.0.0) then	
        do 8004 n = 2, igrid-1
	    star=sta2(n)
	    fin=fi2(n)
          do 8006 ix = star,fin
            temp1 = phimap1(ix) +
     &      phimap1(ix+1)
c
            temp2 = phimap1(ix+lat2) +
     &      phimap1(ix-lat1)
c
            temp3 = phimap1(ix+long2) +
     &      phimap1(ix-long1)
c
       	  phimap2(ix) =phimap2(ix)*om1 + (temp1+temp2+temp3)*sf2(ix)
8006	    continue
8004	  continue
c
	else
c
        do 8104 n = 2, igrid-1
	    star=sta2(n)
	    fin=fi2(n)
            do 8106 ix = star,fin
              temp1 = phimap1(ix) +
     &        phimap1(ix+1)
c
              temp2 = phimap1(ix+lat2) +
     &        phimap1(ix-lat1)
c
              temp3 = phimap1(ix+long2) +
     &        phimap1(ix-long1)
c
       	    phimap2(ix) =phimap2(ix)*om1 + (temp1+temp2+temp3)*sixth
8106	      continue
8104	  continue
	end if
c
C$DIR NO_RECURRENCE 
      do 8010 k=icount2a+1,icount2b
	    ix=idpos(k)
	    temp1=phimap1(ix)*db(1,k)+phimap1(ix+1)*db(2,k)
	    temp2=phimap1(ix-lat1)*db(3,k)+phimap1(ix+lat2)*db(4,k)
          temp3=phimap1(ix-long1)*db(5,k)+phimap1(ix+long2)*db(6,k)
          phimap2(ix)=phimap2(ix) + temp1+temp2+temp3
8010  continue
c reset boundary condition
c
	star=(igrid+2)/2
	iy=(igrid*(igrid+2)/2) - igrid +1
	fin=(igrid*(igrid-1)-1)/2
C$DIR NO_RECURRENCE
	do 8201 ix=star,fin
	  iy=iy+igrid
	  phimap2(iy)=bndx3(ix)
	  phimap2(iy+ihgd2)=bndx4(ix)
8201  continue
c
c
c
C$DIR NO_RECURRENCE 
	do 8011 k=icount1a+1,icount1b
	  temp=qval(k)
	  ix=iqpos(k)
	  phimap2(ix)=phimap2(ix) + temp
8011  continue
c
c
c z periodicity
c
        if(iper(3)) then
          do 8013 iz = 2,(igrid-1)**2,2
	  temp1=bndz(iz)
	  temp2=temp1+idif2z
	  temp3=temp2+inc2za
	  temp4=temp1+inc2zb
	    phimap2(temp1)=phimap1(temp2)
	    phimap2(temp3)=phimap1(temp4)
8013      continue
        end if
c
c y periodicity
c
        if(iper(2)) then
          do 8015 iy = 2,(igrid-1)**2,2
	      temp1=bndy(iy)
	      temp2=temp1+idif2y
	      temp3=temp2+inc2ya
	      temp4=temp1+inc2yb
	      phimap2(temp1)=phimap1(temp2)
	      phimap2(temp3)=phimap1(temp4)
8015	    continue
        end if
c
c x periodicity
c
        if(iper(1)) then
          do 8017 ix = 2,(igrid-1)**2,2
	       temp1=bndx(ix)
	       temp2=temp1+idif2x
	       temp3=temp2+inc2xa
	       temp4=temp1+inc2xb
	       phimap2(temp1)=phimap1(temp2)
	       phimap2(temp3)=phimap1(temp4)
8017       continue
         end if
c
c we also save time by only checking convergence every ten
c iterations, rather than every single iteration.
c
	if(mod(i,icon1).eq.(icon1-1)) then
	do 8051 ix=2,(icgrid+1)/2,icon2
	  phimap3(ix)=phimap2(ix)
8051    continue
	end if
c
	if(gten.ne.0.0) then
	  grden=0.0
	  do ix=1,icount1a
	    iy=iqpos(ix)
	    grden=grden + phimap1(iy)*gval(ix)
	  end do
	  do ix=icount1a+1,icount1b
	    iy=iqpos(ix)
	    grden=grden+phimap2(iy)*gval(ix)
	  end do
	  grdn(i)=grden/2.0
	  if(i.gt.10) then
	    igt=.true.
	    do ix=i-4,i
	      do iy=i-4,i
	        if(abs(grdn(iy)-grdn(ix)).gt.gten) igt=.false.
	      end do
	    end do
	    if(igt)then
	       ires=1
	    end if
	  end if
c
	end if
c
	if((mod(i,icon1).eq.0).or.(ires.eq.1)) then
	  rnorm2=0
        do 8050 ix=2,(icgrid+1)/2,icon2
	    temp2=phimap3(ix)-phimap2(ix)
	    rnorm2=rnorm2+temp2**2
	    rmxch=amax1(rmxch,abs(temp2))
8050    continue
        rmsch = sqrt(float(icon2)*rnorm2/npoint)
	  rmsch2=rmsch
	  rmxch2=rmxch 
	  if(gten.eq.0.) write(6,*) rmsch2,rmxch2,' at ',i,'iterations'
	  if(gten.ne.0.) write(6,*) rmsch2,rmxch2,grden,' at ',i,'iterations'
        if(rmsch.le.res1.and.rmxch.le.res2) ires=1

      end if
c
c check to see if accuracy is sufficient
c
c LOOP
c
      i=i+1
      if(i.le.nlit.and.ires.eq.0) goto 1000
c
c	end of iteration loop
c
c       remap into phimap
c
	do 8701 iy=1,(icgrid-1)/2 
	  ix=iy*2 
	  phimap3(ix-1)= phimap1(iy)
	  phimap3(ix)=   phimap2(iy)
8701	continue 
	iw=1 
	do 8801 iz=1,igrid
	  do 8802 iy=1,igrid
	     do 8803 ix=1,igrid
	        phimap(ix,iy,iz)=phimap3(iw)
	        iw=iw+1
8803	     continue
8802    continue
8801  continue
	phimap(igrid,igrid,igrid)=phimap1((icgrid+1)/2)
c
	icount2a=icount2abac
	icount2b=icount2bbac
c 
	write(6,*)'finished qdiffx linear iterations'
	write(6,*)'# loops                  : ',(i-1)
	write(6,*)'mean,max change (kT/e)   : ',rmsch2,rmxch2
c
	return
	end
