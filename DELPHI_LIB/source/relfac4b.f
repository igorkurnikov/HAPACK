	subroutine relfac(xx,mgrid,
     &  phimap,phimap1,phimap2,phimap3,atmcrg,
     &  iepsmp,iepsmp2,idebmap,
     &  bndx,bndy,bndz,
     &  nlit,nnit,iper,idpos,db,
     &			sf1,sf2,icount2a,icount2b,
     &                  rionst,spec)
c
c
	include 'qdiffpar4.h'
c 

	dimension bndx(*),bndy(*),bndz(*)
	dimension sn1(ngrmax),sn2(ngrmax),sn3(ngrmax)
	dimension db(6,*),sf1(*),sf2(*)
	integer*4 idpos(*)
	logical iper(3)
	integer*4 star,fin,sta1(ngrmax),sta2(ngrmax),fi1(ngrmax),fi2(ngrmax)
c

c
c
c
	sixth = 1./6.
	icgrid=igrid**3
	ihgd=(igrid+1)/2
c
c
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
225     continue
c
c also
c
	lat1= (igrid-1)/2
	lat2= (igrid+1)/2
	long1= (igrid**2 - 1)/2
	long2= (igrid**2 + 1)/2
c 
c set up sn array for lowest eigenstate
c
      i=0
c
	sn1(1)=0.0
	sn1(igrid)=0.0
        sn2(1)=0.0
        sn2(igrid)=0.0
        sn3(1)=0.0
        sn3(igrid)=0.0
c
	do 550 ix=2,igrid-1
        temp=3.1415926536*float(ix-1)/float(igrid-1)
        sn1(ix)=sqrt(2.0)*sin(temp)/sqrt(float(igrid-1))
	sn2(ix)=sn1(ix)
	sn3(ix)=sn1(ix)
550	continue
	if(iper(1)) then
	do 571 ix=1,igrid
571     sn1(ix)=1.0/sqrt(float(igrid))
	end if
	if(iper(2)) then
	do 572 iy=1,igrid
572	sn2(ix)=1.0/sqrt(float(igrid))
	end if
	if(iper(3)) then
	do 573 iz=1,igrid
573	sn3(iz)=1.0/sqrt(float(igrid))
	end if
c
c
	iw=1
	do 301 iz=1,igrid
	  temp3=sn3(iz)
	  do 302 iy=1,igrid
	    temp2=temp3*sn2(iy)
	    do 303 ix=1,igrid
	       phimap3(iw)=temp2*sn1(ix)
	    iw=iw+1
303         continue
302       continue
301     continue
	temp=0.0
	do 500 ix=2,icgrid-1,2
	iy=ix/2
	phimap2(iy)=phimap3(ix)
	temp=temp + phimap3(ix)**2
500	continue 
c
	if(rionst.ne.0.0) then
        do 9004 n = 2, igrid-1
	star=sta1(n)
	fin=fi1(n)
            do 9006 ix = star,fin
            temp1 = phimap2(ix) +
     &      phimap2(ix-1)
c
            temp2 = phimap2(ix+lat1) +
     &      phimap2(ix-lat2)
c
            temp3 = phimap2(ix+long1) +
     &      phimap2(ix-long2)
c
       	phimap1(ix) = (temp1+temp2+temp3)*sf1(ix)
9006	    continue
9004	continue
c
c otherwise the main loop is as below:
c
        else
c
        do 9104 n = 2, igrid-1
	    star=sta1(n)
	    fin=fi1(n)
            do 9106 ix = star,fin
            temp1 = phimap2(ix) +
     &      phimap2(ix-1)
c
            temp2 = phimap2(ix+lat1) +
     &      phimap2(ix-lat2)
c
            temp3 = phimap2(ix+long1) +
     &      phimap2(ix-long2)
c
       	phimap1(ix) =  (temp1+temp2+temp3)*sixth
9106	    continue
9104	continue
        end if
c
c
C$DIR NO_RECURRENCE 
        do 9010 k=1,icount2a
	ix=idpos(k)
	temp1=phimap2(ix-1)*db(1,k)+phimap2(ix)*db(2,k)
	temp2=phimap2(ix-lat2)*db(3,k)+phimap2(ix+lat1)*db(4,k)
        temp3=phimap2(ix-long2)*db(5,k)+phimap2(ix+long1)*db(6,k)
        phimap1(ix)= phimap1(ix) + temp1+temp2+temp3
9010    continue
c
c Now reset boundary values altered in above loops.
c
c 
	star=igrid*(igrid+1)/2
	fin=igrid*(igrid*(igrid-1)-2)/2
C$DIR NO_RECURRENCE
	do 9201 ix=star,fin,igrid
	  phimap1(ix+1)=0.0
	  phimap1(ix+ihgd)=0.0
9201    continue
c
	temp=0.0
	do 1500 ix=1,(icgrid-1)/2
	temp=temp + phimap1(ix)*phimap3(2*ix-1)
1500	continue 
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
9015	    continue
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
9017       continue
         end if
c
c
c Next update phimap3 using the new phimap1
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
       	phimap3(ix) =(temp1+temp2+temp3)*sf2(ix)
8006	    continue
8004	continue
c
	else
c
        do 8104 n = 2, igrid-1
	    star=sta2(n)
	    fin=fi2(n)
            do 8106 ix = star,fin
            temp1 = phimap1(ix) +
     &      phimap1(ix+1)
c
            temp2 = phimap1(ix+lat2) +
     &      phimap1(ix-lat1)
c
            temp3 = phimap1(ix+long2) +
     &      phimap1(ix-long1)
c
       	phimap3(ix) = (temp1+temp2+temp3)*sixth
8106	    continue
8104	continue
	end if
c
C$DIR NO_RECURRENCE 
        do 8010 k=icount2a+1,icount2b
	ix=idpos(k)
	temp1=phimap1(ix)*db(1,k)+phimap1(ix+1)*db(2,k)
	temp2=phimap1(ix-lat1)*db(3,k)+phimap1(ix+lat2)*db(4,k)
        temp3=phimap1(ix-long1)*db(5,k)+phimap1(ix+long2)*db(6,k)
        phimap3(ix)=phimap3(ix) + temp1+temp2+temp3
8010    continue
c reset boundary condition
c
	star=(igrid+2)/2
	iy=(igrid*(igrid+2)/2) - igrid +1
	fin=(igrid*(igrid-1)-1)/2
	ihgd2=ihgd-1
C$DIR NO_RECURRENCE
	do 8201 ix=star,fin
	iy=iy+igrid
	phimap3(iy)=0.0
	phimap3(iy+ihgd2)=0.0
8201    continue
c
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
	    phimap3(temp1)=phimap1(temp2)
	    phimap3(temp3)=phimap1(temp4)
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
	    phimap3(temp1)=phimap1(temp2)
	    phimap3(temp3)=phimap1(temp4)
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
	    phimap3(temp1)=phimap1(temp2)
	    phimap3(temp3)=phimap1(temp4)
8017       continue
         end if
c
c
	temp=0.0
	do 2000 ix=1,(icgrid-1)/2
	temp=temp + (phimap3(ix)*phimap2(ix))
2000	continue 
	spec=(2.0*temp)
	write(6,*) ' '
	write(6,*) 'gauss-seidel spectral radius is',spec
	write(6,*) ' '
c
c
	return
	end
