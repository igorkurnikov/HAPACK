c	Delcx.f	
c
c	This file contains a library of fortran routines used to
c	compute the regular triangulation of a set of points in
c	3D space
c
c	Copyright (C) 2002 Patrice Koehl
c
c	This library is free software; you can redistribute it and/or
c	modify it under the terms of the GNU Lesser General Public
c	License as published by the Free Software Foundation; either
c	version 2.1 of the License, or (at your option) any later version.
c
c	This library is distributed in the hope that it will be useful,
c	but WITHOUT ANY WARRANTY; without even the implied warranty of
c	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
c	Lesser General Public License for more details.
c
c	You should have received a copy of the GNU Lesser General Public
c	License along with this library; if not, write to the Free Software
c	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
c
c	#include "defines.h"
c
c	Setup.f
c
c	Copyright (C) 2002 Patrice Koehl
c 
c	This program gets the coordinates of the N points considered, and stores
c	these into the structures (common blocks) used in all suite of programs
c	Regular3D.f
c
	subroutine setup(coord_sph,rad,nspheres)
c
c	Input:
c	*******
c
c	- nspheres:	number of points to be triangulated
c	- coord_sph:	coordinates of all points (in real*8)
c	- rad:   	weight of each point; this is the radius
c			of the sphere, while the "weight" usually
c			considered in regular triangulations
c			is the square of the radius
c
	integer	npointmax
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)

	parameter	(npointmax=MAX_POINT)
c
	integer		ndigit
	integer		nspheres, npoint,nvertex
	integer		i,j,k,ip,jp
c
	integer		redinfo(npointmax)
	integer		ranlist(npointmax)
c
	real*8		scale,eps,epsd,c_max,s
	real*8		x,xval,radius2
	real*8		y,z,w,xi,yi,zi,wi
	real*8		ranval(3*npointmax)
	real*8		coord_sph(3*npointmax),coord(3*npointmax)
	real*8		coord4(npointmax)
	real*8		radius(npointmax),rad(npointmax)
c
	common	/xyz_vertex/	coord,radius,coord4
	common  /vertex_zone/	npoint,nvertex,redinfo
	common /gmp_info/	scale,eps
c
	save
c
	parameter(EPS_APA= 0.01)
	parameter(SCALE_APA= 100000000.)
	parameter(NDIGIT_APA = 8)

	scale = SCALE_APA
	ndigit = NDIGIT_APA

c	write(6,*) 'SCALE_APA, EPS_APA, NDIGIT_APA', SCALE_APA, EPS_APA,NDIGIT_APA
c
c	Pre-processing:
c	***************
c
c	Eliminate all duplicates. If two points have the same coordinates
c	but different radii, keep the point with the largest radius
c
	npoint = nspheres
c
	c_max = 0
	do 100 i = 1,nspheres
		redinfo(i) = 0
		x = rad(i)
		call truncate_real(x,xval,ndigit)
		radius(i) = xval
		do 50 j = 1,3
			k = 3*(i-1)+j
			x = coord_sph(k)
			call truncate_real(x,xval,ndigit)
			coord(k) = xval
			if(abs(coord(k)).gt.c_max) c_max = abs(coord(k))
50		continue
100	continue
c
	c_max = max(100.d0,c_max)
c
c	Set eps
c
	epsd = 1.0d0
125	epsd = epsd / 2.0d0
	s = 1.0d0 + epsd
	if(s.gt.1.0d0) goto 125
	epsd = 2.0d0 * epsd
c	write(6,*) 'epsd =',epsd
c
	eps = 100d0* c_max**5 * epsd
	eps = max(eps,0.00001)

c
c	write(6,*) 'eps =',eps
c
c	Pre-compute all weights (stored in coord4):
c
	do 200 i = 1,npoint
		radius2 = radius(i)*radius(i)
		coord4(i) = -radius2
		do 150 j = 1,3
			k = 3*(i-1)+j
			x = coord(k)
			coord4(i) = coord4(i) + x*x
150 		continue
200	continue
c
c
c	Check for trivial redundant points: twice the same point
c
	do 300 i = 1,3*npoint
		ranval(i) = coord(i)
300	continue
	call hpsort_three(ranval,ranlist,npoint)
c
	jp = ranlist(1)
	x = coord(3*jp-2)
	y = coord(3*jp-1)
	z = coord(3*jp)
	w = radius(jp)
	do 400 i = 2,npoint
		ip = ranlist(i)
		xi = coord(3*ip-2)
		yi = coord(3*ip-1)
		zi = coord(3*ip)
		wi = radius(ip)
		if((xi-x)**2+(yi-y)**2+(zi-z)**2.le.100.d0*epsd) then
			if(wi.le.w) then
				redinfo(ip) = 1
			else
				redinfo(jp) = 1
				jp = ip
				w = wi
			endif
		else
			x = xi
			y = yi
			z = zi
			w = wi
			jp = ip
		endif
400	continue
c
	return
	end
c
c	Regular3D.f
c
c	Copyright (C) 2002 Patrice Koehl
c 
c	This program computes the regular triangulation of a set
c	of N weighted points in 3D, using the incremental flipping
c	algorithm of Edelsbrunner.
c
c	This implementation is based on the algorithm published in
c	H. Edelsbrunner and N.R. Shah, Algorithmica (1996) 15: 223-241
c
c	1) Algorithm:
c	*************
c
c	Briefly, the algorithm works as follows:
c
c	- first, a large tetrahedron initialises the program. All
c	four vertices of this tetrahedron are set at "infinite"
c
c	- All N points are added one by one.
c
c	- For each point:
c
c		- localize the tetrahedron in the current regular
c		triangulation that contains this point
c
c		- test if the point is redundant; if yes, remove
c
c		- If the point is not redundant, insert in the
c		tetrahedron : this is a "1-4" flip
c
c		- collect all "link facets" (i.e. all triangles
c		in tetrahedron containing the new point, that face
c		this new point) that are not regular.
c
c		- for each non-regular link facet, check if it
c		is "flippable". If yes, perform a "2-3", "3-2"
c		or "1-4" flip. Add new link facets in the list,
c		if needed.
c
c		- when link facet list if empty, move to next 
c		point
c
c	- Remove "infinite" tetrahedron, i.e. tetrahedron with
c	one vertice at "infinite"
c
c	- collect all remaining tetrahedron, define convex hull,
c	and exit.
c
c	2) Data structure:
c	******************
c
c	I maintain a minimal data structure that includes only 
c	the tetrahedrons of the triangulation (triangles,
c	edges and vertices are implicit).
c
c	For each tetrahedron, I store:
c
c	- the index of its four vertices
c	- pointers to its neighbours (4 maximum).
c		neighbor(i) is the tetrahedron that shares
c		all vertices of the tetrahedron considered, except i
c		(0 if the corresponding face is on the convex hull)
c	- its status: 1 "active" (i.e. part of the triangulation), 0 inactive
c	- its orientation
c
c	3) number representation:
c	**************************
c
c	I use double precision floating points. However, if one of
c	the geometricaly test becomes "imprecise", I switch to
c	arbitrary precision arithmetics (using the gmp package).
c
c	All primitives have therefore a floating point filter
c
	subroutine regular3(nredundant,list_redundant)
c
c	Input:
c	*******

c
c	- npoint  :	number of points to be triangulated
c	- coord   :	coordinates of all points (in real*8)
c	- coord4  :	weight of each point; this is the radius
c			of the sphere, while the "weight" usually
c			considered in regular triangulations
c			is the square of the radius
c
c	Output:
c	********
c
c	- ntetra:	number of tetrahedron in the final
c			regular triangulation
c	- tetra:	for each tetrahedron, gives the position of its
c			four vertices in ascending order
c	- hull:		for each tetrahedron on the convex hull,
c			gives the index (local index, i.e. in 1-4)
c			of the vertex NOT on the convex hull
c
c	Both input and output are exchanged via common blocks
c
c	First include all parameters defining the maximum size of
c	the arrays (in Fortran, all sizes are defined prior to
c	compilation). If the dimensions were set to small,
c	you need to edit this file, and recompile the program.
c
c	Now declare all variables
c
	integer	npointmax,ntetra_max,nfreemax
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter	(npointmax=MAX_POINT)
	parameter	(ntetra_max=MAX_TETRA)
	parameter	(nfreemax = MAX_FREE)
c
	integer		i,j
	integer		ival,iredundant,tetra_loc,tetra_in
	integer		npoint,nvertex
	integer		ntetra
	integer		iseed
	integer		nfree,nkill,n_new
	integer		nredundant,iweight
c
	integer		redinfo(npointmax)
	integer		list_redundant(npointmax)
c
	integer*1	tetra_status(ntetra_max)
	integer*1	tetra_orient(ntetra_max)
	integer*1	tetra_nindex(4,ntetra_max)
	integer*1	infpoint(npointmax)
c
	integer		list_new(ntetra_max)
	integer		tetra(4,ntetra_max)
	integer		tetra_neighbour(4,ntetra_max)
	integer		free(nfreemax),kill(nfreemax)
c
	integer		ranlist(npointmax)
c
	real*4		ran2
c
	real*8		scale,eps
c
	real*8		coord(3*npointmax)
	real*8		radius(npointmax)
	real*8		coord4(npointmax)
	real*8		ranval(npointmax)
c
c	for simplicity, information are stored in common blocks
c
	common	/xyz_vertex/	coord,radius,coord4
	common  /vertex_zone/	npoint,nvertex,redinfo
	common  /tetra_zone/	ntetra,tetra,tetra_neighbour
     	common /tetra_stat/	tetra_status,tetra_orient,
     1				tetra_nindex
	common /infinite/	infpoint
	common /freespace/	nfree,nkill,free,kill
	common /gmp_info/	scale,eps
	common /update/		n_new,list_new
	common /flags/		iweight
c
	save
c
c	Initialisation:
c	****************
c
c	Add four infinite points and initialize first tetrahedron
c
c	

	do 320 i = npoint,1,-1
		coord4(i+4) = coord4(i)
		redinfo(i+4) = redinfo(i)
		radius(i+4) = radius(i)
		infpoint(i+4) = 0
320     continue
c
        do 340 i = 3*npoint,1,-1
		coord(i+12) = coord(i)
340     continue
c
	nvertex = npoint + 4
c
	do 360 i = 1,12
		coord(i) = 0
360	continue
c
	do 380 i = 1,4
		radius(i) = 0
		coord4(i) = 0
380	continue
c

	call transfer_all_to_gmp(coord,radius,scale,nvertex)
	
c
	infpoint(1) = 1
	infpoint(2) = 1
	infpoint(3) = 1
	infpoint(4) = 1
c
	ntetra = 1
	tetra(1,ntetra) = 1
	tetra(2,ntetra) = 2
	tetra(3,ntetra) = 3
	tetra(4,ntetra) = 4
c
	tetra_neighbour(1,ntetra) = 0
	tetra_neighbour(2,ntetra) = 0
	tetra_neighbour(3,ntetra) = 0
	tetra_neighbour(4,ntetra) = 0
c
	tetra_status(ntetra) = 1
	tetra_orient(ntetra) = -1
c
c	Randomize order in which points are added
c
	iseed = -1
	do 300 i = 1,npoint
		ranval(i) = ran2(iseed)
300	continue
	call hpsort_key(ranval,ranlist,npoint)
c
c	Initialise "free" space to 0
c
	nfree = 0
	n_new = 0
c
c	Build regular triangulation
c	****************************
c
c	Now loop over all points:
c
	tetra_in = -1
	iweight = 1
c

	do 400 i = 1,npoint
c
		ival = ranlist(i) + 4
c
		n_new = 0
c
		if(redinfo(ival).eq.1) goto 400
c
c		first locate the point in the list of known tetrahedra
c
		iseed = -1
		tetra_loc = tetra_in
		call locate_jw(iseed,ival,tetra_loc,iredundant)
c
c		write(6,*) (coord(3*(ival-1)+j),j=1,3)
c		write(6,*) 'point ',ival,' found in tetrahedron # ',tetra_loc
c
c		If the point is redundant, move to next point
c
		if(iredundant.eq.1) then
			redinfo(ival) = 1
			goto 400
		endif
c
c		Otherwise, add point to tetrahedron : 1-4 flip
c
		call flipjw_1_4(ival,tetra_loc)
c
c		Now scan link_facet list, and flip till list is empty
c
		call flipjw
		
c
c		At this stage, I should have a regular triangulation
c		of the i+4 points (i-th real points+4 "infinite" points)
c		Now add another point
c
400	continue

c
c	call reorder_tetra(n_new,list_new)
c
c	I have the regular triangulation: I need to remove the
c	simplices including infinite points, and define the
c	convex hull
c

	call remove_inf
c
c	write(6,*) 'ntetra =',ntetra
c	do 410 i = 1,ntetra
c		write(6,*) (tetra(j,i),j=1,4),tetra_status(i)
c410	continue
c
c	Now I peel off flat tetrahedra at the boundary of the DT
c
c	call peel
c
c	I should be done now! go back to calling program
c
	nredundant = 0
	do 500 i = 1,npoint
		if(redinfo(i+4).eq.1) then
			nredundant = nredundant + 1
			list_redundant(nredundant) = i
		endif
500	continue
c
c
	return
	end
c
c	Locate_jw.f		Version 1 12/17/2001	Patrice Koehl
c
c	This subroutine locates the tetrahedron containing a new
c	point to be added in the triangulation
c
c	This implementation of the point location scheme
c	uses a "jump-and-walk" technique: first, N active
c	tetrahedra are chosen at random. The "distances" between
c	these tetrahedra and the point to be added are computed,
c	and the tetrahedron closest to the point is chosen as
c	a starting point. The program then "walks" from that tetrahedron
c	to the point, till we find a tetrahedron that contains
c	the point.
c	It also checks if the point is redundant in the current
c	tetrahedron. It it is, the search terminates.
c
	subroutine locate_jw(iseed,ival,tetra_loc,iredundant)
c
c	Input:
c	*******
c
c	- ival:	index of the points to be located
c
c	Output:
c	********
c
c	- tetra_loc:	tetrahedron containing the point
c	- iredundant:	flag for redundancy: 0 is not redundant,
c			1 otherwise
c
c	Define array size
c
	integer	ntetra_max
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter (ntetra_max = MAX_TETRA)
c
c	Declare variables
c
	integer	ival,itetra,iorient,idx,iseed
	integer	ntetra,tetra_in
	integer	tetra_loc,iredundant
	integer	a,b,c,d
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
	integer	tetra(4,ntetra_max)
	integer	tetra_neighbour(4,ntetra_max)
c
	logical	test_in,test_red
c
	common /tetra_zone/ ntetra,tetra,tetra_neighbour
	common /tetra_stat/ tetra_status,tetra_orient,
     1			    tetra_nindex
c
	save
c
c	Start at the root of the history dag: tetra(1)
c
	iredundant = 0
c
	if(ntetra.eq.1) then
		tetra_loc = 1
		return
	endif
c
	tetra_in = tetra_loc
c
c	write(6,*) 'before jump, iseed =',iseed
c
	if(tetra_in.lt.0) then
		call jump(iseed,ival,itetra)
	else
		itetra = tetra_in + 1
	endif
c	write(6,*) 'from jump, itetra =',itetra
c
100	continue
c
	a = tetra(1,itetra)
	b = tetra(2,itetra)
	c = tetra(3,itetra)
	d = tetra(4,itetra)
c
	iorient = tetra_orient(itetra)
c
	call inside_tetra_jw(ival,a,b,c,d,iorient,test_in,
     1		test_red,idx)
c
	
c	write(6,*) 'a,b,c,d,in,red :',a,b,c,d,test_in,test_red
c
	if(test_in) goto 200
c
	itetra = tetra_neighbour(idx,itetra)
	goto 100
c
200	continue
c
	tetra_loc = itetra
c	write(6,*) 'a,b,c,d,in,red :',a,b,c,d,test_in,test_red
c
c	Now that we have the tetrahedron (at a given layer of the
c	history dag), check if point is redundant
c
	if(test_red) iredundant = 1
c
	return
	end
c
c	Jump.f		Version 1 3/6/2002	Patrice Koehl
c
c	This subroutine picks N active tetrahedra at random,
c	computes the distance of the points to be inserted
c	to each of this tetrahedron, and selects the closest
c	one
c
	subroutine jump(iseed,ival,itetra)
c
	integer	npointmax,ntetra_max
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter	(npointmax=MAX_POINT)
	parameter	(ntetra_max=MAX_TETRA)
c
	integer	i,j,k,ival
	integer	Ntry,N,Nkeep
	integer	iseed
	integer	ntetra,itetra
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
	integer*1 infpoint(npointmax)
c
	integer	tetra(4,ntetra_max)
	integer	tetra_neighbour(4,ntetra_max)
	integer	list(100)
c
	real*4	ran2
c
	real*8	dist,distmin
	real*8	coord(3*npointmax),coord4(npointmax)
	real*8	radius(npointmax)
	real*8	xval(3),xa(3)
c
	common  /xyz_vertex/coord,radius,coord4
	common /tetra_zone/ ntetra,tetra,tetra_neighbour
	common /tetra_stat/ tetra_status,tetra_orient,
     1			    tetra_nindex
	common /infinite/   infpoint
c
	save
c
	Ntry = 20
c
	N = 200
c
50	continue
c
	Nkeep = 0
	do 100 i = 1,N
c
		j = int(ntetra*ran2(iseed)) + 1
		j = min(j,ntetra)
		if(tetra_status(j).ne.1) goto 100
		do 75 k = 1,Nkeep
			if(j.eq.list(k)) goto 100
75		continue
		Nkeep = Nkeep + 1
		list(Nkeep) = j
		if(nkeep.eq.Ntry) goto 125
c
100	continue
c
125	continue
c
	if(Nkeep.eq.0) goto 50
c
	do 200 i = 1,3
		xval(i) = coord(3*(ival-1)+i)
200	continue
c
	do 275 j = 1,4
		if(infpoint(tetra(j,list(1))).eq.1) goto 275
		do 250 i = 1,3
			xa(i) = coord(3*(tetra(j,list(1))-1)+i)
250		continue
		goto 285
275	continue
285	continue
c
	distmin = 0
	do 300 i = 1,3
		distmin = distmin + (xval(i)-xa(i))*
     1		(xval(i)-xa(i))
300	continue
c
	itetra = list(1)
c
	do 600 i = 2,Nkeep
c
		do 425 k = 1,4
			if(infpoint(tetra(k,list(i))).eq.1) goto 425
			do 400 j = 1,3
				xa(j) = coord(3*(tetra(k,list(i))-1)+j)
400			continue
			goto 450
425		continue
450		continue
c
		dist = 0
		do 500 j = 1,3
			dist = dist + (xval(j)-xa(j))*
     1			(xval(j)-xa(j))
500		continue
c
		if(dist.lt.distmin) then
			distmin = dist
			itetra = list(i)
		endif
600	continue
c
	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c	Inside_tetra_jw.f	Version 1 2/11/2002	Patrice Koehl
c
c	This subroutine tests if a point p is inside a tetrahedron
c	defined by four points (a,b,c,d) (with orientation "iorient")
c	If p is found inside the tetrahedron, it also checks if it
c	is redundant
c
c	Computation is done in floating point, but it is switched to
c	multiple precision if the result is imprecise
c
	subroutine inside_tetra_jw(p,a,b,c,d,iorient,is_in,redundant,
     1				ifail)
c
c	Input:
c		- p	:	index of point to be checked
c		- a,b,c,d:	the four vertices of the tetrahedron
c		- iorient:	orientation of the tetrahedron
c
c	Output:
c		- is_in:	"true" if p in the tetrahedron, "false"
c				otherwise
c		- redundant	"true" is p is redundant
c		- ifail:	In case p is not inside the tetrahedron,
c				ifail gives the index of the face that
c				fails the orientation test
c
	integer	npointmax
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter	(npointmax=MAX_POINT)
c
	integer	i,j,k,l,m
	integer	p,a,b,c,d
	integer	ia,ib,ic,id,ie,idx
	integer	ic1,ic5,ic1_k,ic1_l,sign,sign5,sign_k,sign_l
	integer	nswap,iswap,ninf
	integer	iorient,ifail
	integer	val
	integer iweight
c
	integer	list(4)
	integer inf4_1(4),sign4_1(4)
	integer inf4_2(4,4),sign4_2(4,4)
	integer	sign4_3(4)
	integer	inf5_2(4,4),sign5_2(4,4)
	integer inf5_3(4),sign5_3(4)
	integer	order1(3,4),order2(2,6),order3(2,6)
c
	integer*1	infpoint(npointmax)
c
	real*8	Sij_1,Sij_2,Sij_3,Skl_1,Skl_2,Skl_3
	real*8	det_pijk,det_pjil,det_pkjl,det_pikl,det_pijkl
	real*8	eps,scale
	real*8	detij(3)
	real*8	coordp(3),i_p(4),j_p(4),k_p(4),l_p(4)
	real*8	coord(3*npointmax)
        real*8	coord4(npointmax)
	real*8	radius(npointmax)
c
	logical	test_pijk,test_pjil,test_pkjl,test_pikl
	logical is_in,redundant
c
	common  /xyz_vertex/	coord,radius,coord4
	common  /gmp_info/	scale,eps
	common  /infinite/	infpoint
	common  /flags/		iweight
c
	save
c
	data	inf4_1 /2,2,1,1/
	data	sign4_1 /-1,1,1,-1/
c
	data	inf4_2 /0,2,3,3,
     1                  2,0,3,3,
     2                  3,3,0,1,
     3                  3,3,1,0/
c
        data	sign4_2 / 0,1,-1,1,
     1                   -1,0,1,-1,
     2                    1,-1,0,1,
     3                   -1,1,-1,0/
c
	data	sign4_3 /-1,1,-1,1/
c
	data inf5_2	/0,2,1,1,
     1 			2,0,1,1,
     2			1,1,0,1,
     3			1,1,1,0/
c
	data sign5_2    /0,-1,-1,1,
     1			1,0,-1,1,
     2			1,1,0,1,
     3			-1,-1,-1,0/
c
	data inf5_3	/1,1,3,3/
	data sign5_3	/1,1,-1,1/
c
	data order1 /3,2,4,1,3,4,2,1,4,1,2,3/
	data order2/3,4,4,2,2,3,1,4,3,1,1,2/
	data order3/1,2,1,3,1,4,2,3,2,4,3,4/
c
c	If (i,j,k,l) is the tetrahedron in positive orientation, we need
c	to test:
c		(p,i,j,k)
c		(p,j,i,l)
c		(p,k,j,l)
c		(p,i,k,l)
c	If all four are positive, than p is inside the tetrahedron.
c	All four tests relies on the sign of the corresponding 4x4
c	determinant. Interestingly, these four determinants share
c	some common lines, which can be used to speed up the computation.
c
c	Let us consider or example:
c
c	det(p,i,j,k) = | p(1) p(2) p(3) 1|
c		       | i(1) i(2) i(3) 1|
c		       | j(1) j(2) j(3) 1|
c		       | k(1) k(2) k(3) 1|
c
c	p appears in each determinant. The corresponding line can therefore
c	be substraced from all 3 other lines . Using the example above,
c	we find:
c
c	det(i,j,k,l) = - |ip(1) ip(2) ip(3)|
c		         |jp(1) jp(2) jp(3)|
c			 |kp(1) kp(2) kp(3)|
c
c	where :	xp(m) = x(m) - p(m) for x = i,j,k and m = 1,2,3
c
c	Now we notice that the first two lines of det(p,i,j,k) and
c	det(p,i,j,l) are the same.
c
c	Let us define: Sij_3 = |ip(1) ip(2)| Sij_2 = |ip(1) ip(3)| 
c			       |jp(1) jp(2)|         |jp(1) jp(3)|
c	and Sij_1 = |ip(2) ip(3)|
c		    |jp(2) jp(3)|
c
c	We find:
c		det(p,i,j,k) = - kp(1)*Sij_1 + kp(2)*Sij_2 - kp(3)*Sij_3
c	and:
c		det(p,j,i,l) =   lp(1)*Sij_1 - lp(2)*Sij_2 + lp(3)*Sij_3
c
c	Similarly, if we define: 
c
c	Skl_3 = |kp(1) kp(2)|	Skl_2 = |kp(1) kp(3)|	Skl_1 = |kp(2) kp(3)|
c		|lp(1) lp(2)|		|lp(1) lp(3)|		|lp(2) lp(3)|
c
c	We find:
c		det(p,k,j,l) = jp(1)*Skl_1 - jp(2)*Skl_2 + jp(3)*Skl_3
c	and:
c		det(p,i,k,l) = -ip(1)*Skl_1 + ip(2)*Skl_2 - ip(3)*Skl_3
c
c	Furthermore:
c
c	det(p,i,j,k,l) = -ip(4)*det(p,k,j,l)-jp(4)*det(p,i,k,l)
c			 -kp(4)*det(p,j,i,l)-lp(4)*det(p,i,j,k)
c
c	The equations above hold for the general case; special care is
c	required to take in account infinite points (see below)
c
	is_in = .false.
	redundant = .false.
	
	list(1) = a
	list(2) = b
	list(3) = c
	list(4) = d
c
	ninf = infpoint(a) + infpoint(b) + infpoint(c) + infpoint(d)
c
c	"General case" : no infinite point
c
	do 100 m = 1,3
		coordp(m) = coord(3*p-3+m)
100	continue
c
	if(ninf.eq.0) then
c
c		Define coordinates (with i=a, j=b, k=c and l=d)
c		(no need to change notation, just bad habit to use i,j,k,l
c		instead of a,b,c,d !)
c
		do 200 m = 1,3
			i_p(m) = coord(3*a-3+m) - coordp(m)
			j_p(m) = coord(3*b-3+m) - coordp(m)
			k_p(m) = coord(3*c-3+m) - coordp(m)
			l_p(m) = coord(3*d-3+m) - coordp(m)
200		continue
c
c		Now compute 2x2 determinants Sij and Skl
c
		Sij_1 = i_p(2)*j_p(3) - i_p(3)*j_p(2)
		Sij_2 = i_p(1)*j_p(3) - i_p(3)*j_p(1)
		Sij_3 = i_p(1)*j_p(2) - i_p(2)*j_p(1)
c
		Skl_1 = k_p(2)*l_p(3) - k_p(3)*l_p(2)
		Skl_2 = k_p(1)*l_p(3) - k_p(3)*l_p(1)
		Skl_3 = k_p(1)*l_p(2) - k_p(2)*l_p(1)
c
c	Now perform tests

c
c	Start with is_in = .false. :
c
		is_in = .false.
c
		det_pijk = -k_p(1)*Sij_1 + k_p(2)*Sij_2 - k_p(3)*Sij_3
		det_pijk = det_pijk*iorient
		test_pijk = abs(det_pijk).gt.eps
		if(test_pijk.and.det_pijk.gt.0) then
			ifail = 4
			return
		endif
c
c	We check all other four determinants
c
		det_pjil = l_p(1)*Sij_1 - l_p(2)*Sij_2 + l_p(3)*Sij_3
		det_pjil = det_pjil*iorient
		test_pjil = abs(det_pjil).gt.eps
		if(test_pjil.and.det_pjil.gt.0) then
			ifail = 3
			return
		endif
c
		det_pkjl = j_p(1)*Skl_1 - j_p(2)*Skl_2 + j_p(3)*Skl_3
		det_pkjl = det_pkjl*iorient
		test_pkjl = abs(det_pkjl).gt.eps
		if(test_pkjl.and.det_pkjl.gt.0) then
			ifail = 1
			return
		endif
c
		det_pikl = -i_p(1)*Skl_1 + i_p(2)*Skl_2 - i_p(3)*Skl_3
		det_pikl = det_pikl*iorient
		test_pikl = abs(det_pikl).gt.eps
		if(test_pikl.and.det_pikl.gt.0) then
			ifail = 2
			return
		endif
c
c		At this stage, either all four determinants are positive,
c		or one of the determinant is not precise enough, and
c		we need to switch the MP
c		In this case, since we may need SoS, we have to rank
c		the indices
c
		if(.not.test_pijk) then
			call valsort4(p,a,b,c,ia,ib,ic,id,nswap)
			call sos_minor4_gmp(ia,ib,ic,id,val)
			val = val*nswap*iorient
			if(val.eq.1) then
				ifail = 4
				return
			endif
		endif
c
		if(.not.test_pjil) then
			call valsort4(p,b,a,d,ia,ib,ic,id,nswap)
			call sos_minor4_gmp(ia,ib,ic,id,val)
			val = val*nswap*iorient
			if(val.eq.1) then
				ifail = 3
				return
			endif
		endif
c
		if(.not.test_pkjl) then
			call valsort4(p,c,b,d,ia,ib,ic,id,nswap)
			call sos_minor4_gmp(ia,ib,ic,id,val)
			val = val*nswap*iorient
			if(val.eq.1) then
				ifail = 1
				return
			endif
		endif
c
		if(.not.test_pikl) then
			call valsort4(p,a,c,d,ia,ib,ic,id,nswap)
			call sos_minor4_gmp(ia,ib,ic,id,val)
			val = val*nswap*iorient
			if(val.eq.1) then
				ifail = 2
				return
			endif
		endif
c
c		If we have gone that far, p is inside the tetrahedron
c
		is_in = .true.
		if(iweight.eq.0) return
c
c		Now we check if p is redundant
c
		i_p(4) = coord4(a) - coord4(p)
		j_p(4) = coord4(b) - coord4(p)
		k_p(4) = coord4(c) - coord4(p)
		l_p(4) = coord4(d) - coord4(p)
c
		det_pijkl = -i_p(4)*det_pkjl - j_p(4)*det_pikl
     1			   -k_p(4)*det_pjil - l_p(4)*det_pijk
c
c	No need to multiply by iorient, since all minors contains iorient...
c
		if(abs(det_pijkl).lt.eps) then
			call valsort5(p,a,b,c,d,ia,ib,ic,id,
     1			ie,nswap)
			call sos_minor5_gmp(ia,ib,ic,
     1			id,ie,val)
			det_pijkl = val*nswap*iorient
		endif
		redundant = det_pijkl.lt.0
c
	elseif(ninf.eq.1) then
c
c		We know that one of the 4 vertices a,b,c or d is 
c		infinite
c		To find which one it is, we use a map between
c		(inf(a),inf(b),inf(c),inf(d)) and X, where inf(i)
c		is 1 if i is infinite, 0 otherwise, and X = 1,2,3,4
c		if a,b,c or d are infinite, respectively.
c		A good mapping function is:
c		X = 3 - inf(a) - inf(a) -inf(b) + inf(d)
c
		idx=3-infpoint(a)-infpoint(a)-infpoint(b)+infpoint(d)
		l = list(idx)
c
c		The three finite points:
c
		i = list(order1(1,idx))
		j = list(order1(2,idx))
		k = list(order1(3,idx))
c
		ic1 = inf4_1(l)
		sign = sign4_1(l)
c
c	let us look at the four determinant we need to compute:
c
c	det_pijk	: unchanged
c	det_pjil	: 1 infinite point (l), becomes det3_pji
c			  where det3_pij = |p(ic1) p(ic2) 1|
c					   |i(ic1) i(ic2) 1|
c					   |j(ic1) j(ic2) 1|
c			  and ic1 and ic2 depends on which infinite
c			  (ic2 is always 3)
c			  point is considered
c	det_pkjl	: 1 infinite point (l), becomes det3_pkj
c	det_pikl	: 1 infinite point (l), becomes det3_pik
c
c	Get coordinates
c
		do 300 m = 1,3
			i_p(m) = coord(3*i-3+m) - coordp(m)
			j_p(m) = coord(3*j-3+m) - coordp(m)
			k_p(m) = coord(3*k-3+m) - coordp(m)
300		continue
c
		detij(1) = i_p(1)*j_p(3) - i_p(3)*j_p(1)
		detij(2) = i_p(2)*j_p(3) - i_p(3)*j_p(2)
		detij(3) = i_p(1)*j_p(2) - i_p(2)*j_p(1)
c
c	Now perform tests
c
c	Start with is_in = .false. :
c
		is_in = .false.
c
		det_pijk = -k_p(1)*detij(2)+k_p(2)*detij(1)
     1				- k_p(3)*detij(3)
		det_pijk = det_pijk*iorient
		test_pijk = abs(det_pijk).gt.eps
		if(test_pijk.and.det_pijk.gt.0) then
			ifail = idx
			return
		endif
c
		det_pjil = -detij(ic1)*sign*iorient
		test_pjil = abs(det_pjil).gt.eps
		if(test_pjil.and.det_pjil.gt.0) then
			ifail = order1(3,idx)
			return
		endif
c
		det_pkjl = k_p(ic1)*j_p(3) - k_p(3)*j_p(ic1)
		det_pkjl = sign*det_pkjl*iorient
		test_pkjl = abs(det_pkjl).gt.eps
		if(test_pkjl.and.det_pkjl.gt.0) then
			ifail = order1(1,idx)
			return
		endif
c
		det_pikl = i_p(ic1)*k_p(3) - i_p(3)*k_p(ic1)
		det_pikl = sign*det_pikl*iorient
		test_pikl = abs(det_pikl).gt.eps
		if(test_pikl.and.det_pikl.gt.0) then
			ifail = order1(2,idx)
			return
		endif
c
c		At this stage, either all four determinants are positive,
c		or one of the determinant is not precise enough, and
c		we need to switch the MP
c
c
		if(.not.test_pijk) then
			call valsort4(p,i,j,k,ia,ib,ic,id,nswap)
			call sos_minor4_gmp(ia,ib,ic,id,val)
			val = val*nswap*iorient
			if(val.eq.1) then
				ifail = idx
				return
			endif
		endif
c
		if(.not.test_pjil) then
			call valsort3(p,j,i,ia,ib,ic,nswap)
			call sos_minor3_gmp(ia,ib,ic,
     1				ic1,3,val)
			val = val*sign*nswap*iorient
			if(val.eq.1) then
				ifail = order1(3,idx)
				return
			endif
		endif
c
		if(.not.test_pkjl) then
			call valsort3(p,k,j,ia,ib,ic,nswap)
			call sos_minor3_gmp(ia,ib,ic,
     1				ic1,3,val)
			val = val*sign*nswap*iorient
			if(val.eq.1) then
				ifail = order1(1,idx)
				return
			endif
		endif
c
		if(.not.test_pikl) then
			call valsort3(p,i,k,ia,ib,ic,nswap)
			call sos_minor3_gmp(ia,ib,ic,
     1				ic1,3,val)
			val = val*sign*nswap*iorient
			if(val.eq.1) then
				ifail = order1(2,idx)
				return
			endif
		endif
c
c	If we have gone so far, p is inside the tetrahedron
c
		is_in = .true.
c
c		Now we check if p is redundant
c
c		since det_pijkl = det_pijk >1
c		p cannot be redundant !
c
		redundant = .false.
c
	elseif(ninf.eq.2) then
c
c		We know that two of the 4 vertices a,b,c or d are
c		infinite
c		To find which one it is, we use a map between
c		(inf(a),inf(b),inf(c),inf(d)) and X, where inf(i)
c		is 1 if i is infinite, 0 otherwise, and X = 1,2,3,4,5,6
c		if (a,b), (a,c), (a,d), (b,c), (b,d), or (c,d) are
c		infinite, respectively
c		A good mapping function is:
c		X = 3 - inf(a) - inf(a) +inf(c) + inf(d) + inf(d)
c
		idx = 3 -infpoint(a) -infpoint(a) +infpoint(c)
     1			+ infpoint(d) + infpoint(d)
c
c		The two infinite points :
c
		k = list(order3(1,idx))
		l = list(order3(2,idx))
c
c		The two finite points
c
		i = list(order2(1,idx))
		j = list(order2(2,idx))
c
		ic1_k = inf4_1(k)
		ic1_l = inf4_1(l)
		sign_k = sign4_1(k)
		sign_l = sign4_1(l)
		ic1 = inf4_2(k,l)
		sign = sign4_2(k,l)
c
c	Get coordinates
c
		do 400 m = 1,3
			i_p(m) = coord(3*i-3+m) - coordp(m)
			j_p(m) = coord(3*j-3+m) - coordp(m)
400		continue
c
c	Perform test; first set is_in .false.
c
		is_in = .false.
c
c	det_pijk is now det3_pij with k as infinite point
c
		det_pijk = i_p(ic1_k)*j_p(3)-i_p(3)*j_p(ic1_k)
		det_pijk = det_pijk*sign_k*iorient
		test_pijk = abs(det_pijk).gt.eps
		if(test_pijk.and.det_pijk.gt.0) then
			ifail = order3(2,idx)
			return
		endif
c
c	det_pjil is now det3_pji with l as infinite point
c
		det_pjil = i_p(3)*j_p(ic1_l)-i_p(ic1_l)*j_p(3)
		det_pjil = det_pjil*sign_l*iorient
		test_pjil = abs(det_pjil).gt.eps
		if(test_pjil.and.det_pjil.gt.0) then
			ifail = order3(1,idx)
			return
		endif
c
c	det_pkjl is now -det2_pj (k,l infinite)
c
		det_pkjl = j_p(ic1)*sign*iorient
		test_pkjl = abs(det_pkjl).gt.eps
		if(test_pkjl.and.det_pkjl.gt.0) then
			ifail = order2(1,idx)
			return
		endif
c
c	det_pikl is now det2_pi (k,l infinite)
c
		det_pikl = -i_p(ic1)*sign*iorient
		test_pikl = abs(det_pikl).gt.eps
		if(test_pikl.and.det_pikl.gt.0) then
			ifail = order2(2,idx)
			return
		endif
c
c		At this stage, either all four determinants are positive,
c		or one of the determinant is not precise enough, and
c		we need to switch the MP
c
		if(.not.test_pijk) then
			call valsort3(p,i,j,ia,ib,ic,nswap)
			call sos_minor3_gmp(ia,ib,ic,
     1				ic1_k,3,val)
			val = val*sign_k*nswap*iorient
			if(val.eq.1) then
				ifail = order3(2,idx)
				return
			endif
		endif
c
		if(.not.test_pjil) then
			call valsort3(p,j,i,ia,ib,ic,nswap)
			call sos_minor3_gmp(ia,ib,ic,
     1				ic1_l,3,val)
			val = val*sign_l*nswap*iorient
			if(val.eq.1) then
				ifail = order3(1,idx)
				return
			endif
		endif
c
		if(.not.test_pkjl) then
			call valsort2(p,j,ia,ib,nswap)
			call sos_minor2_gmp(ia,ib,ic1,val)
			val = -val*sign*nswap*iorient
			if(val.eq.1) then
				ifail = order2(1,idx)
				return
			endif
		endif
c
		if(.not.test_pikl) then
			call valsort2(p,i,ia,ib,nswap)
			call sos_minor2_gmp(ia,ib,ic1,val)
			val = val*sign*nswap*iorient
			if(val.eq.1) then
				ifail = order2(2,idx)
				return
			endif
		endif
c
c	Again, if we have gone so far, p is inside the tetrahedron
c
		is_in = .true.
		redundant = .false.
		if(iweight.eq.0) return
c
c		Now we check if p is redundant
c
c		det_pijkl becomes det3_pij
c
		ic5 = inf5_2(k,l)
		sign5 = sign5_2(k,l)
		det_pijkl = i_p(ic5)*j_p(3)-i_p(3)*j_p(ic5)
		if(abs(det_pijkl).lt.eps) then
			call valsort3(p,i,j,ia,ib,ic,nswap)
			call sos_minor3_gmp(ia,ib,ic,
     1			ic5,3,val)
			det_pijkl = val*nswap
		endif
		det_pijkl = det_pijkl*sign5*iorient
c
		redundant = det_pijkl.lt.0
c
	elseif(ninf.eq.3) then
c
c		We know that three of the 4 vertices a,b,c or d are
c		infinite
c		To find which one is finite, we use a map between
c		(inf(a),inf(b),inf(c),inf(d)) and X, where inf(i)
c		is 1 if i is infinite, 0 otherwise, and X = 1,2,3,4
c		if a,b,c or d are finite, respectively.
c		A good mapping function is:
c		X = 1 + inf(a) + inf(a) +inf(b) - inf(d)
c
		idx=1+infpoint(a)+infpoint(a)+infpoint(b)-infpoint(d)
		i = list(idx) 
		j = list(order1(1,idx))
		k = list(order1(2,idx))
		l = list(order1(3,idx))
c
c	Index of the "missing" infinite point (i.e. the fourth infinite
c	point)
c
		call missinf_sign(j,k,l,ie,iswap)
c
c	Get coordinates
c
		do 500 m = 1,3
			i_p(m) = coord(3*i-3+m) - coordp(m)
500		continue
c
c	Perform test; first set is_in to .false.
c
		is_in = .false.
c
c	det_pijk is now - det2_pi (missing j,k)
c
		det_pijk = i_p(inf4_2(j,k))*iorient*sign4_2(j,k)
		test_pijk = abs(det_pijk).gt.eps
		if(test_pijk.and.det_pijk.gt.0) then
			ifail = order1(3,idx)
			return
		endif
c
c	det_pjil is now det2_pi (missing j,l)
c
		det_pjil = -i_p(inf4_2(j,l))*iorient*sign4_2(j,l)
		test_pjil = abs(det_pjil).gt.eps
		if(test_pjil.and.det_pjil.gt.0) then
			ifail = order1(2,idx)
			return
		endif
c
c	det_pkjl is now det1_p
c
		det_pkjl = iorient*iswap*sign4_3(ie)
		if(det_pkjl.gt.0) then
			ifail = idx
			return
		endif
c
c	det_ikl is now - det2_pi (missing k,l)
c
		det_pikl = i_p(inf4_2(k,l))*iorient*sign4_2(k,l)
		test_pikl = abs(det_pikl).gt.eps
		if(test_pikl.and.det_pikl.gt.0) then
			ifail = order1(1,idx)
			return
		endif
c
c		At this stage, either all four determinants are positive,
c		or one of the determinant is not precise enough, and
c		we need to switch the MP
c
c
		if(.not.test_pijk) then
			call valsort2(p,i,ia,ib,nswap)
			call sos_minor2_gmp(ia,ib,
     1				inf4_2(j,k),val)
			val = -val*sign4_2(j,k)*iorient*nswap
			if(val.eq.1) then
				ifail = order1(3,idx)
				return
			endif
		endif
c
		if(.not.test_pjil) then
			call valsort2(p,i,ia,ib,nswap)
			call sos_minor2_gmp(ia,ib,
     1				inf4_2(j,l),val)
			val = val*sign4_2(j,l)*iorient*nswap
			if(val.eq.1) then
				ifail = order1(2,idx)
				return
			endif
		endif
c
		if(.not.test_pikl) then
			call valsort2(p,i,ia,ib,nswap)
			call sos_minor2_gmp(ia,ib,
     1				inf4_2(k,l),val)
			val = -val*sign4_2(k,l)*iorient*nswap
			if(val.eq.1) then
				ifail = order1(1,idx)
				return
			endif
		endif
c

		is_in = .true.
		redundant = .false.
		if(iweight.eq.0) return
c
c	Now check for redundancy
c
c		det_pijkl becomes -det2_pi
c
		ic1 = inf5_3(ie)
		sign5 = sign5_3(ie)
		det_pijkl = -i_p(ic1)
		if(abs(det_pijkl).lt.eps) then
			call valsort2(p,i,ia,ib,nswap)
			call sos_minor2_gmp(ia,ib,ic1,val)
			det_pijkl = val*nswap
		endif
		det_pijkl = - iorient*det_pijkl*sign5*iswap
c
		redundant = det_pijkl.lt.0
c
	else
c
c	In the case all four points ia,ib,ic, and id are infinite,
c	then is_in = .true. and redundant = .false.
c
		is_in = .true.
		redundant = .false.
c
	endif
c
	return
	end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCccccccccccccccccccccccccccccccccccccccccccccccccccCCCCCCCCCCCCCC


c	Regular_convex.f	Version 1 2/25/2002	Patrice Koehl
c
c	This subroutine checks if a link facet (a,b,c) is locally
c	regular, as well as if the union of two tetrahedra
c	(a,b,c,p) and (a,b,c,o) that connects to the facet is convex. 
c
c	Computation is done in floating point, but it is switched to
c	multiple precision if the result is imprecise
c
c	In floating point, there is no need to order points in lexicographic
c	order prior to computing a determinant. This simplification is no
c	more true if we need to apply SoS: consequently, as soon as
c	the program swtiches to GMP (i.e. multi-precision), I first
c	order the points, using a series of routines valsort*, where
c	* can be 2,3,4 or 5, all contained in the file valsort
c
	subroutine regular_convex(a,b,c,p,o,itest_abcp,
     1			regular,convex,test_abpo,test_bcpo,test_capo)
c
c	Input:
c		- a,b,c:	the three points defining the
c				link facet
c		- p:		the current point inserted in the
c				DT
c		- o:		the fourth point of the tetrahedron
c				that attaches to (a,b,c), opposite
c				to the tetrahedron (a,b,c,p,o)
c		- itest_abcp	orientation of the tetrahedron
c				(abcp)
c
c	Output:
c		- convex	"true" of (abcp)U(abco) is convex, "false"
c				otherwise
c		- regular	"true" if (abc) is locally regular, in which
c				case it does not matter if convex!
c
	integer	npointmax
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)

	parameter	(npointmax=MAX_POINT)
c
	integer	p,a,b,c,o,i,j,k,l,m
	integer	ia,ib,ic,id,ie
	integer	ninf,infp,info,iswap,iswap2,idx,val
	integer	icol1,sign1,icol2,sign2,icol4,sign4,icol5,sign5
	integer	itest_abcp
c
	integer	list(3)
	integer	order(2,3)
	integer inf4_1(4),sign4_1(4)
	integer inf4_2(4,4),sign4_2(4,4)
	integer sign4_3(4)
	integer	inf5_2(4,4),sign5_2(4,4)
	integer inf5_3(4),sign5_3(4)
	integer	order1(3,3)
c
	integer*1	infpoint(npointmax)
c
	real*8	eps,scale
	real*8	det_abpo,det_bcpo,det_capo,det_abcpo,det_abpc
	real*8	a_p(4),b_p(4),c_p(4),o_p(4)
	real*8	i_p(3),j_p(3)
	real*8	Mbo(3), Mca(3),Mjo(3),Mio(3)
	real*8	coordp(3)
	real*8	coord(3*npointmax),radius(npointmax)
        real*8	coord4(npointmax)
c
	logical convex,regular,test_abpo,test_bcpo,test_capo
	logical testc(3)
c
	common  /xyz_vertex/	coord,radius,coord4
	common  /gmp_info/	scale,eps
	common  /infinite/	infpoint
c
	save
c
	data	inf4_1 /2,2,1,1/
	data	sign4_1 /-1,1,1,-1/
c
	data	inf4_2 /0,2,3,3,
     1                  2,0,3,3,
     2                  3,3,0,1,
     3                  3,3,1,0/
c
        data	sign4_2 / 0,1,-1,1,
     1                   -1,0,1,-1,
     2                    1,-1,0,1,
     3                   -1,1,-1,0/
c
	data    sign4_3 /-1,1,-1,1/
c
	data inf5_2	/0,2,1,1,
     1 			2,0,1,1,
     2			1,1,0,1,
     3			1,1,1,0/
c
	data sign5_2    /0,-1,-1,1,
     1			1,0,-1,1,
     2			1,1,0,1,
     3			-1,-1,-1,0/
c
	data inf5_3     /1,1,3,3/
	data sign5_3    /1,1,-1,1/
c
	data order1	/1,2,3,3,1,2,2,3,1/
	data order	/2,3,3,1,1,2/
c
c	To test if the union of the two tetrahedron is convex, we check the
c	position of o with respect to the three faces (a,b,p), (b,c,p)
c	and (c,a,p) of (a,b,c,p). 
c	To do that, we evaluate the three determinants:
c		det(a,b,p,o)
c		det(b,c,p,o)
c		det(c,a,p,o)
c	If the three determinants are positive, and det(a,b,c,p) is negative,
c	then the union is convex
c	Also, if the three determinants are negative, and det(a,b,c,p) is 
c	positive, then the union is convex
c	In all other cases, the union is non convex
c
c	The regularity is tested by computing det(a,b,c,p,o) 
c
c	The computations required are very similar to those used in
c	inside_tetra.f . Look at comments there to decipher this subroutine
c
c	Let us first count how many infinite points we have:
c	(except o)
c	only a and/or b and/or c can be infinite:
c
	ninf = infpoint(a) + infpoint(b) + infpoint(c)
c
	list(1) = a
	list(2) = b
	list(3) = c
c
	do 100 m = 1,3
		coordp(m) = coord(3*p-3+m)
100	continue
c
c	"General case" : no infinite point
c
	if(ninf.eq.0) then
c
c		First, a simple case: if o is infinite, then
c		det(a,b,c,p,o) = -det(a,b,c,p) and consequently 
c		(a,b,c,p,o) is regular:nothing to do!
c
		if(infpoint(o).eq.1) then
			regular = .true.
			return
		endif
c
c	The three determinants det(a,b,p,o), det(b,c,p,o), and det(c,a,p,o)
c	are "real" 4x4 determinants. 
c	First, we substract the row corresponding to p from the other row,
c	and develop with respect to p. The determinants become:
c
c	det(a,b,p,o)= -| ap(1) ap(2) ap(3) |
c		       | bp(1) bp(2) bp(3) |
c		       | op(1) op(2) op(3) |
c
c	det(b,c,p,o)= -| bp(1) bp(2) bp(3) |
c		       | cp(1) cp(2) cp(3) |
c		       | op(1) op(2) op(3) |
c
c	det(c,a,p,o)= -| cp(1) cp(2) cp(3) |
c		       | ap(1) ap(2) ap(3) |
c		       | op(1) op(2) op(3) |
c
c	where ip(j) = i(j) - p(j) for all i in {a,b,c,o} and j in {1,2,3}
c
c	We compute two types of minors:
c
c		Mbo_ij = bp(i)op(j) - bp(j)op(i)
c	and
c		Mca_ij = cp(i)ap(j) - cp(j)op(i)
c
c	We store Mbo_12 in Mbo(3), Mbo_13 in Mbo(2),...
c
c	Get coordinates
c
		do 200 m = 1,3
			a_p(m) = coord(3*a-3+m) - coordp(m)
			b_p(m) = coord(3*b-3+m) - coordp(m)
			c_p(m) = coord(3*c-3+m) - coordp(m)
			o_p(m) = coord(3*o-3+m) - coordp(m)
200		continue
c
		a_p(4) = coord4(a) - coord4(p)
		b_p(4) = coord4(b) - coord4(p)
		c_p(4) = coord4(c) - coord4(p)
		o_p(4) = coord4(o) - coord4(p)
c
c	Now compute 2x2 determinants Mbo and Mca
c
		Mbo(1) = b_p(2)*o_p(3) - b_p(3)*o_p(2)
		Mbo(2) = b_p(1)*o_p(3) - b_p(3)*o_p(1)
		Mbo(3) = b_p(1)*o_p(2) - b_p(2)*o_p(1)
c
		Mca(1) = c_p(2)*a_p(3) - c_p(3)*a_p(2)
		Mca(2) = c_p(1)*a_p(3) - c_p(3)*a_p(1)
		Mca(3) = c_p(1)*a_p(2) - c_p(2)*a_p(1)
c
c	Now,
c
		det_abpo = - a_p(1)*Mbo(1)+a_p(2)*Mbo(2)
     1				-a_p(3)*Mbo(3)
		det_bcpo = c_p(1)*Mbo(1) - c_p(2)*Mbo(2) 
     1				+ c_p(3)*Mbo(3)
		det_capo = - o_p(1)*Mca(1) + o_p(2)*Mca(2) 
     1				- o_p(3)*Mca(3)
c
c
c	We also compute:
c
		det_abpc = - b_p(1)*Mca(1) + b_p(2)*Mca(2) 
     1				- b_p(3)*Mca(3)
c
c	Now we compute:
c		det(a,b,c,p,o) = | a(1) a(2) a(3) a(4) 1 |
c				 | b(1) b(2) b(3) b(4) 1 |
c				 | c(1) c(2) c(3) c(4) 1 |
c				 | p(1) p(2) p(3) p(4) 1 |
c				 | o(1) o(2) o(3) o(4) 1 |
c	We first substract row p :
c
c		det(a,b,c,p,o) = - | ap(1) ap(2) ap(3) ap(4) |
c				   | bp(1) bp(2) bp(3) bp(4) |
c				   | cp(1) cp(2) cp(3) cp(4) |
c			 	   | op(1) op(2) op(3) op(4) |
c
c	By developping with respect to the last column, we get:
c
		det_abcpo=-a_p(4)*det_bcpo-b_p(4)*det_capo 
     1			- c_p(4)*det_abpo + o_p(4)*det_abpc
c		write(6,*) 'det_abcpo :',det_abcpo
c
c	Test if (a,b,c,p,o) regular, in which case there is no need
c	to flip
c
		if(abs(det_abcpo).lt.eps) then
			call valsort5(a,b,c,p,o,ia,ib,ic,id,ie,iswap)
			call sos_minor5_gmp(ia,ib,ic,id,ie,val)
c			if(det_abcpo*val*iswap.lt.0) 
c     1				write(6,*) 'det_pabcpo,sos :',det_abcpo,val*iswap
c			write(6,*) 'Regular det_abcpo,sos :',det_abcpo,
c     1			val*iswap,itest_abcp
			det_abcpo = val*iswap
		endif
c
		if(det_abcpo*itest_abcp.lt.0) then
			regular = .true.
			return
		endif
		regular = .false.
c
c	If not regular, we test for convexity
c
		if(abs(det_abpo).lt.eps) then
			call valsort4(a,b,p,o,ia,ib,ic,id,iswap)
			call sos_minor4_gmp(ia,ib,ic,id,val)
			det_abpo = val*iswap
		endif
		if(abs(det_bcpo).lt.eps) then
			call valsort4(b,c,p,o,ia,ib,ic,id,iswap)
			call sos_minor4_gmp(ia,ib,ic,id,val)
			det_bcpo = val*iswap
		endif
		if(abs(det_capo).lt.eps) then
			call valsort4(c,a,p,o,ia,ib,ic,id,iswap)
			call sos_minor4_gmp(ia,ib,ic,id,val)
			det_capo = val*iswap
		endif
c
		test_abpo = det_abpo.gt.0
		test_bcpo = det_bcpo.gt.0
		test_capo = det_capo.gt.0
c
		convex = .false.
		if((itest_abcp*det_abpo).gt.0) return
		if((itest_abcp*det_bcpo).gt.0) return
		if((itest_abcp*det_capo).gt.0) return
		convex = .true.
c
c	Second case: either a,b or c is infinite:
c
	elseif(ninf.eq.1) then
c
c	Let us define as X the infinite point, and (i,j) the pair of finite
c	points.
c	If X = a, (i,j) = (b,c)
c	If X = b, (i,j) = (c,a)
c	If X = c, (i,j) = (a,b)
c	If we define inf(a) = 1 if a infinite, 0 otherwise,
c	then idx_X  = 2 - inf(a) + inf(c)
c
		idx = 2 -infpoint(a) + infpoint(c)
		infp = list(idx)
		i = list(order(1,idx))
		j = list(order(2,idx))
c
c	Get the coordinates
c
		do 300 m = 1,3
			i_p(m) = coord(3*i-3+m) - coordp(m)
			j_p(m) = coord(3*j-3+m) - coordp(m)
			o_p(m) = coord(3*o-3+m) - coordp(m)
300		continue
c
c	First case:	o is finite
c
		if(infpoint(o).eq.0) then
c
			icol1 = inf4_1(infp)
			sign1 = sign4_1(infp)
c
c	The three 4x4 determinants become:
c
c		-det(i,p,o) [X missing]
c		det(j,p,o) [X missing]
c		det(i,j,p,o)
c
c	And the 5x5 determinant becomes:
c
c		- det(i,j,p,o)
c
			Mjo(1) = j_p(1)*o_p(3) - j_p(3)*o_p(1)
			Mjo(2) = j_p(2)*o_p(3) - j_p(3)*o_p(2)
			Mjo(3) = j_p(1)*o_p(2) - j_p(2)*o_p(1)
c 
c	The correspondence between a,b,c and i,j is not essential
c	We use here the corresponce for a infinite; in the
c	two other cases (b infinite or c infinite), we would
c	have computed the same determinants, but they would
c	not come in the same order
c
			det_abpo = i_p(icol1)*o_p(3)-i_p(3)*o_p(icol1)
			if(abs(det_abpo).lt.eps) then
				call valsort3(i,p,o,ia,ib,ic,iswap)
				call sos_minor3_gmp(ia,ib,ic,icol1,
     1				3,val)
				det_abpo = -val*iswap
			endif
			det_abpo = det_abpo*sign1
			det_capo = - Mjo(icol1)
			if(abs(det_capo).lt.eps) then
				call valsort3(j,p,o,ia,ib,ic,iswap)
				call sos_minor3_gmp(ia,ib,ic,icol1,
     1				3,val)
				det_capo = val*iswap
			endif
			det_capo = det_capo*sign1
			det_bcpo =-i_p(1)*Mjo(2)+i_p(2)*Mjo(1)
     1					-i_p(3)*Mjo(3)
			if(abs(det_bcpo).lt.eps) then
				call valsort4(i,j,p,o,ia,ib,
     1					ic,id,iswap)
c				write(6,*) 'i,j,p,o:',i,j,p,o
				call sos_minor4_gmp(ia,ib,ic,id,val)
				det_bcpo = val*iswap
			endif
			det_abcpo = -det_bcpo
c			write(6,*) 'ninf = 1,det_abcpo :',det_abcpo
c
		else
c
c	Second case: o is infinite
c
			info = o
c
c	The three 4x4 determinants become:
c
c		-det(i,p) [o,X missing]
c		det(j,p) [o,X missing]
c		det(i,j,p) [o missing]
c
c	And the 5x5 determinant becomes:
c
c		det(i,j,p) [o,X missing]
c
			icol1 = inf4_2(info,infp)
			sign1 = sign4_2(info,infp)
c
			icol2 = inf4_1(info)
			sign2 = sign4_1(info)
c
			icol5 = inf5_2(info,infp)
			sign5 = sign5_2(info,infp)
c
			det_abpo =-i_p(icol1)*sign1
			if(abs(det_abpo).lt.eps) then
				call valsort2(i,p,ia,ib,iswap)
				call sos_minor2_gmp(ia,ib,icol1,val)
				det_abpo = -val*iswap*sign1
			endif
			det_capo =j_p(icol1)*sign1
			if(abs(det_capo).lt.eps) then
				call valsort2(j,p,ia,ib,iswap)
				call sos_minor2_gmp(ia,ib,icol1,val)
				det_capo = val*iswap*sign1
			endif
			det_bcpo =i_p(icol2)*j_p(3)
     1				-i_p(3)*j_p(icol2)
			if(abs(det_bcpo).lt.eps) then
				call valsort3(i,j,p,ia,ib,ic,iswap)
				call sos_minor3_gmp(ia,ib,ic,icol2,
     1				3,val)
				det_bcpo = val*iswap
			endif
			det_bcpo =det_bcpo*sign2
			det_abcpo=i_p(icol5)*j_p(3)
     1				-i_p(3)*j_p(icol5)
			if(abs(det_abcpo).lt.eps) then
				call valsort3(i,j,p,ia,ib,ic,iswap)
				call sos_minor3_gmp(ia,ib,ic,icol5,
     1					3,val)
				det_abcpo = val*iswap
			endif
			det_abcpo= det_abcpo*sign5
c			write(6,*) 'ninf = 1b,det_abcpo :',det_abcpo
c
		endif
c
c	Test if (a,b,c,p,o) regular, in which case there is no need
c	to flip
c
		if(det_abcpo*itest_abcp.lt.0) then
			regular = .true.
			return
		endif
		regular = .false.
c
c	If not regular, we test for convexity
c
		testc(1) = det_abpo.gt.0
		testc(2) = det_bcpo.gt.0
		testc(3) = det_capo.gt.0
		test_abpo = testc(order1(1,idx))
		test_bcpo = testc(order1(2,idx))
		test_capo = testc(order1(3,idx))
c
		convex = .false.
		if((itest_abcp*det_abpo).gt.0) return
		if((itest_abcp*det_bcpo).gt.0) return
		if((itest_abcp*det_capo).gt.0) return
		convex = .true.
c
c	Now we consider the case where two points are infinite
c
	elseif(ninf.eq.2) then
c
c	Let us define as (k,l) the two infinite points, and i the
c	point that is finite
c	If i = a, (k,l) = (b,c)
c	If i = b, (k,l) = (c,a)
c	If i = c, (k,l) = (a,b)
c
c	Again: i = 2 + inf(a) - inf(c)
c
		idx = 2 + infpoint(a) - infpoint(c)
		i = list(idx)
		k = list(order(1,idx))
		l = list(order(2,idx))
c
c	Get the coordinates
c
		do 400 m = 1,3
			i_p(m) = coord(3*i-3+m) - coordp(m)
			o_p(m) = coord(3*o-3+m) - coordp(m)
400		continue
c
c	First case: o is finite
c
		if(infpoint(o).eq.0) then
c
c	The three 4x4 determinants become:
c
c		det(i,p,o) [k missing]
c		-det(i,p,o) [l missing]
c		S*det(p,o) [k,l missing, with S =1 if k<l, -1 otherwise]
c
c	The 5x5 determinants become:
c
c		S*det(i,p,o) [k,l missing, with S=1 if k<l, -1 otherwise]
c
			icol1 = inf4_1(k)
			sign1 = sign4_1(k)
			icol2 = inf4_1(l)
			sign2 = sign4_1(l)
			icol4 = inf4_2(k,l)
			sign4 = sign4_2(k,l)
			icol5 = inf5_2(k,l)
			sign5 = sign5_2(k,l)
c
			Mio(1) = i_p(1)*o_p(3) - i_p(3)*o_p(1)
			Mio(2) = i_p(2)*o_p(3) - i_p(3)*o_p(2)
			Mio(3) = i_p(1)*o_p(2) - i_p(2)*o_p(1)
c 
c	The correspondence between a,b,c and i,j,k is not essential
c	We use here the correspondence for a finite; in the
c	two other cases (b finite or c finite), we would
c	have computed the same determinants, but they would
c	not come in the same order
c
			det_abpo = -Mio(icol1)*sign1
			if(abs(det_abpo).lt.eps) then
				call valsort3(i,p,o,ia,ib,ic,iswap)
				call sos_minor3_gmp(ia,ib,ic,icol1,
     1					3,val)
				det_abpo = val*iswap*sign1
			endif
			det_capo =  Mio(icol2)*sign2
			if(abs(det_capo).lt.eps) then
				call valsort3(i,p,o,ia,ib,ic,iswap)
				call sos_minor3_gmp(ia,ib,ic,icol2,
     1					3,val)
				det_capo = -val*iswap*sign2
			endif
			det_bcpo = -o_p(icol4)*sign4
			if(abs(det_bcpo).lt.eps) then
				call valsort2(p,o,ia,ib,iswap)
				call sos_minor2_gmp(ia,ib,icol4,val)
				det_bcpo = val*sign4*iswap
			endif
			det_abcpo = - Mio(icol5)*sign5
			if(abs(det_abcpo).lt.eps) then
				call valsort3(i,p,o,ia,ib,ic,iswap)
				call sos_minor3_gmp(ia,ib,ic,icol5,
     1				3,val)
				det_abcpo = val*iswap*sign5
			endif
c			write(6,*) 'ninf = 2,det_abcpo :',det_abcpo
c
		else
c
c	Second case: o is infinite
c
			info = o
c
c	The three 4x4 determinants become:
c
c		det(i,p) [o,k missing]
c		-det(i,p) [o,l missing]
c		Const [o,k,l missing]
c
c	The 5x5 determinants become:
c
c		Const*det(i,p) [o,k,l missing]
c	
			icol1 = inf4_2(info,k)
			sign1 = sign4_2(info,k)
			icol2 = inf4_2(info,l)
			sign2 = sign4_2(info,l)
c
			call missinf_sign(info,k,l,icol4,iswap)
c
			det_abpo = i_p(icol1)*sign1
			if(abs(det_abpo).lt.eps) then
				call valsort2(i,p,ia,ib,iswap2)
				call sos_minor2_gmp(ia,ib,icol1,val)
				det_abpo = val*iswap2*sign1
			endif
			det_capo = -i_p(icol2)*sign2
			if(abs(det_capo).lt.eps) then
				call valsort2(i,p,ia,ib,iswap2)
				call sos_minor2_gmp(ia,ib,icol2,val)
				det_capo = -val*iswap2*sign2
			endif
			det_bcpo = sign4_3(icol4)*iswap
			det_abcpo = sign5_3(icol4)*iswap
     1				*i_p(inf5_3(icol4))
			if(abs(det_abcpo).lt.eps) then
				call valsort2(i,p,ia,ib,iswap2)
				call sos_minor2_gmp(ia,ib,
     1				inf5_3(icol4),val)
				det_abcpo = val*iswap2*iswap*
     1					sign5_3(icol4)
			endif
c			write(6,*) 'ninf = 2b,det_abcpo :',det_abcpo
c
		endif
c
c	Test if (a,b,c,p,o) regular, in which case there is no need
c	to flip
c
		if(det_abcpo*itest_abcp.lt.0) then
			regular = .true.
			return
		endif
		regular = .false.
c
c	If not regular, we test for convexity
c
		testc(1) = det_abpo.gt.0
		testc(2) = det_bcpo.gt.0
		testc(3) = det_capo.gt.0
		test_abpo = testc(order1(1,idx))
		test_bcpo = testc(order1(2,idx))
		test_capo = testc(order1(3,idx))
c
		convex = .false.
		if((itest_abcp*det_abpo).gt.0) return
		if((itest_abcp*det_bcpo).gt.0) return
		if((itest_abcp*det_capo).gt.0) return
		convex = .true.
c
c	We cannot have all three points a,b,c infinite: in this case,
c	the facet a,b,c would be on the convex hull!
c
	elseif(ninf.eq.3) then
		write(6,*) 'This must be an error...'
	endif
c
	return
	end
c
c	Missinf_sign.f	Version1 3/1/2002	Patrice Koehl
c
c	This subroutine takes as input the indices of three "infinite"
c	points (between 1 and 4), finds the index of the "missing"
c	infinite point, and gives the signature of the permutation
c	required to put the three infinite point in order
c
	subroutine missinf_sign(i,j,k,l,sign)
c
c	Input:	i,j,k:		the three known infinite points
c
c	Output:	l		the "missing" infinite point
c		sign		the signature of the permutation
c				that orders i,j,k
c
	integer	i,j,k,l,sign
	integer	a,b,c,d
c
	save
c
	l = 10 -i -j -k
c
	a = i
	b = j
	c = k
c
	sign = 1
c
	if(a.gt.b) then
		d = a
		a = b
		b = d
		sign = -sign
	endif
c
	if(a.gt.c) then
		d = a
		a = c
		c = d
		sign = -sign
	endif
c
	if(b.gt.c) then
		sign = -sign
	endif
c
	return
	end
C	valsort.f	Version 1 3/2/2002	Patrice Koehl
c
c	This file contains a series of routines that sorts integer
c	values in ascending orders, and keep track of the number of
c	flip required (such as to define the signature of the permutation
c	that transforms the un-ordered data into an ordered array)
c
c	In all these routines, the input values are kept unaffected,
c	and new sorted output values are generated
c
c	1. Sort two numbers a and b
c
	subroutine valsort2(a,b,ia,ib,iswap)
c
	integer	a,b,ia,ib,iswap
c
	save
c
	iswap = 1
	if(a.gt.b) then
		ia = b
		ib = a
		iswap = -iswap
	else
		ia = a
		ib = b
	endif
c
	return
	end
c
c	2. Sort three numbers a, b and c
c
	subroutine valsort3(a,b,c,ia,ib,ic,iswap)
c
	integer	a,b,c,ia,ib,ic,iswap,temp
c
	save
c
	call valsort2(a,b,ia,ib,iswap)
c
	ic = c
c
	if(ib.gt.ic) then
		temp = ib
		ib = ic
		ic = temp
		iswap = -iswap
		if(ia.gt.ib) then
			temp = ia
			ia = ib
			ib = temp
			iswap = -iswap
		endif
	endif
c
	return
	end
c
c	3. Sort four numbers a, b, c and d
c
	subroutine valsort4(a,b,c,d,ia,ib,ic,id,iswap)
c
	integer	a,b,c,d,ia,ib,ic,id,iswap,temp
c
	save
c
	call valsort3(a,b,c,ia,ib,ic,iswap)
c
	id = d
c
	if(ic.gt.id) then
		temp = ic
		ic = id
		id = temp
		iswap = -iswap
		if(ib.gt.ic) then
			temp = ib
			ib = ic
			ic = temp
			iswap = -iswap
			if(ia.gt.ib) then
				temp = ia
				ia = ib
				ib = temp
				iswap = -iswap
			endif
		endif
	endif
c
	return
	end
c
c	4. Sort five numbers a, b, c, d and e
c
	subroutine valsort5(a,b,c,d,e,ia,ib,ic,id,ie,iswap)
c
	integer	a,b,c,d,e,ia,ib,ic,id,ie,iswap,temp
c
	save
c
	call valsort4(a,b,c,d,ia,ib,ic,id,iswap)
c
	ie = e
c
	if(id.gt.ie) then
		temp = id
		id = ie
		ie = temp
		iswap = -iswap
		if(ic.gt.id) then
			temp = ic
			ic = id
			id = temp
			iswap = -iswap
			if(ib.gt.ic) then
				temp = ib
				ib = ic
				ic = temp
				iswap = -iswap
				if(ia.gt.ib) then
					temp = ia
					ia = ib
					ib = temp
					iswap = -iswap
				endif
			endif
		endif
	endif
c
	return
	end
c	Flipjw.f	Version 1 12/17/2001	Patrice Koehl
c
c	After a point has been inserted, this subroutine goes over the 
c	link_facet list to restore regularity. When a link_facet is found
c	non_regular and "flippable" (see below), the program attempts
c	to flip it. If the flip is successful, new link_facets are added
c	on the queue.
c	The subroutine ends when the link facet is empty
c
c	This version of flip calls the "jw" flip subroutines, i.e. flip
c	subroutines that do not store a dag
c
c
	subroutine flipjw
c
c	Input:
c	********
c		- nlink_facet:  when this program is called, the
c				link_facet queue contains four triangles,
c				derived from the tetrahedron in which the
c				new point is added. Each triangle is defined
c				by two tetrahedra, defined by link_facet
c		- link_facet:	the four link facets
c
c	Include array dimensions
c
	integer	ntetra_max,nfreemax,nfacet_max
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter(MAX_FACET = 10000)

	parameter	(ntetra_max=MAX_TETRA)
	parameter	(nfreemax = MAX_FREE)
	parameter	(nfacet_max=MAX_FACET)
c
c	Define variables
c
	integer	j
	integer	p,o,a,b,c
	integer	ierr,ifind,nlink_facet
	integer	itetra,jtetra
	integer	tetra_ab,tetra_ac,tetra_bc
	integer	iorder,ntetra
	integer	ireflex,iflip
	integer	idx_p,idx_o,itest_abcp
	integer	idx_a,idx_b,idx_c
	integer	nfree,nkill,nkill_top,ns
c
	integer	ishare,jshare,kshare,lshare
	integer	idxi,idxj,idxk,idxl
	integer	ia,ib,ic,ii,ij
	integer	itetra_touch(3),jtetra_touch(3)
	integer	itetra_idx(3),jtetra_idx(3)
	integer	ktetra_touch(2,3)
	integer	ktetra_idx(2,3)
	integer	itouch(2),jtouch(2),ktouch(2)
	integer	i_idx(2),j_idx(2),k_idx(2)
	integer	tetra_flip(3),list_flip(3)
c
	integer table32(3,3),table32_2(2,3),table41(3,3)
	integer	table41_2(2,3)
	integer	vert_flip(5)
	integer free(nfreemax),kill(nfreemax)
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
	integer	link_facet(2,nfacet_max)
	integer	link_index(2,nfacet_max)
	integer	tetra(4,ntetra_max)
	integer	tetra_neighbour(4,ntetra_max)
c
	logical test,test_or(2,3),regular,convex
	logical test_abpo,test_abpc,test_capo,test_acpb
	logical test_bcpo,test_bcpa,test_acpo
c
	common  /tetra_zone/	ntetra,tetra,tetra_neighbour
	common  /tetra_stat/	tetra_status,tetra_orient,
     1				tetra_nindex
	common  /link_zone/	nlink_facet,link_facet,link_index
	common  /freespace/	nfree,nkill,free,kill
c
	save
c
	data table32 /1,2,3,1,3,2,3,1,2/
	data table32_2/1,2,1,3,2,3/
	data table41 /2,1,3,1,2,3,1,3,2/
	data table41_2/1,1,2,1,2,2/
c
	nkill_top = nint(nfreemax*0.9)
c
c	Loop over all link facets
c
	j = 0
c
100	if(j.eq.nlink_facet) goto 200
c
	if(nkill.ge.nkill_top) then
c		
		nkill =  nkill_top
c		ns = nfree
c		nfree = min(nfree+nkill,nkill_top)
c		do 150 j = ns+1,nfree
c			free(j) = kill(j-ns)
c150		continue
c		nkill = 0
c
	endif
c
c
c
c	write(6,*) 'j,nlink_facet :',j,nlink_facet
c
	j = j + 1
c
c	First defined the two tetrahedra that contains the link facet as
c	itetra and jtetra
c
	itetra = link_facet(1,j)
	jtetra = link_facet(2,j)
	idx_p  = link_index(1,j)
	idx_o  = link_index(2,j)
c
c
c	If the link facet is on the convect hull, discard
c
	if(itetra.eq.0.or.jtetra.eq.0) goto 100
c
c	If these tetrahedra have already been discarded, discard this
c	link facet
c
	if(tetra_status(itetra).ne.1) then
		if(tetra_status(jtetra).ne.1) then
			goto 100
		else
			itetra = tetra_neighbour(idx_o,jtetra)
			idx_p = tetra_nindex(idx_o,jtetra)
		endif
	endif
	if(tetra_status(jtetra).ne.1) then
		jtetra = tetra_neighbour(idx_p,itetra)
		idx_o = tetra_nindex(idx_p,itetra)
	endif
c
c	Let us define the vertices of the two tetrahedra:
c	itetra:		a,b,c,p
c	jtetra:		a,b,c,o
c
	a = tetra(1,itetra)
	b = tetra(2,itetra)
	c = tetra(3,itetra)
	p = tetra(4,itetra)
c
	o = tetra(idx_o,jtetra)
c
	itest_abcp = tetra_orient(itetra)
c
c	Check for local regularity (and convexity, at very little
c	extra cost)
c
	call regular_convex(a,b,c,p,o,itest_abcp,regular,convex,
     1	test_abpo,test_bcpo,test_capo)
c
c	write(6,*) 'regular, convex :',regular,convex
c
c	if the link facet is locally regular, discard
c
	if(regular) goto 100
c
c	Define neighbors of the facet on itetra and jtetra
c
	call define_facet(itetra,jtetra,idx_o,
     1		itetra_touch,itetra_idx,jtetra_touch,jtetra_idx)
c
	test_abpc = itest_abcp.ne.1
c
c	After discarding the trivial case, we now test if the tetrahedra
c	can be flipped. 
c
c	At this stage, I know that the link facet is not locally
c	regular. I still don't know if it is "flippable"
c
c	I first check if {itetra} U {jtetra} is convex. If it is, I
c	perform a 2-3 flip (this is the convexity test performed
c	at the same time as the regularity test)
c
	if(convex) then
		call flipjw_2_3(itetra,jtetra,p,o,
     1		itetra_touch,itetra_idx,jtetra_touch,jtetra_idx,
     2		test_abpo,test_bcpo,test_capo,ierr)
		goto 100
	endif
c
c	The union of the two tetrahedra is not convex...
c	I now check the edges of the triangle in the link facet, and
c	check if they are "reflexes" (see definition in Edelsbrunner and
c	Shah, Algorithmica (1996), 15:223-241)
c
	ireflex = 0
	iflip = 0
c
c	First check edge (ab): 
c		- (ab) is reflex iff o and c lies on opposite sides of
c		the hyperplane defined by (abp). We therefore test the
c		orientation of (abpo) and (abpc): if they differ (ab)
c		is reflex
c		- if (ab) is reflex, we test if it is of degree 3.
c		(ab) is of degree 3 if it is shared by 3 tetrahedra,
c		namely (abcp), (abco) and (abpo). The first two are itetra
c		and jtetra, so we only need to check if (abpo) exists.
c		since (abpo) contains p, (abp) should then be a link facet
c		of p, so we test all tetrahedra that define link facets
c
c
	if(test_abpo.neqv.test_abpc) then
c
		ireflex = ireflex + 1
c
		call find_tetra(itetra,3,a,b,o,ifind,
     1		tetra_ab,idx_a,idx_b)
c
		if(ifind.eq.1) then
			iflip = iflip + 1
			tetra_flip(iflip) = tetra_ab
			list_flip(iflip) = 1
			ktetra_touch(1,iflip)=tetra_neighbour(idx_a,
     1						tetra_ab)
			ktetra_touch(2,iflip)=tetra_neighbour(idx_b,
     1						tetra_ab)
			ktetra_idx(1,iflip)=tetra_nindex(idx_a,
     1						tetra_ab)
			ktetra_idx(2,iflip)=tetra_nindex(idx_b,
     1						tetra_ab)
			test_or(1,iflip) = test_bcpo
			test_or(2,iflip) = .not.test_capo
		endif
c
	endif
c
c	Now check edge (ac): 
c		- (ac) is reflex iff o and b lies on opposite sides of
c		the hyperplane defined by (acp). We therefore test the
c		orientation of (acpo) and (acpb): if they differ (ac)
c		is reflex
c		- if (ac) is reflex, we test if it is of degree 3.
c		(ac) is of degree 3 if it is shared by 3 tetrahedra,
c		namely (abcp), (abco) and (acpo). The first two are itetra
c		and jtetra, so we only need to check if (acpo) exists.
c		since (acpo) contains p, (acp) should then be a link facet
c		of p, so we test all tetrahedra that define link facets
c
	test_acpo = .not.test_capo
	test_acpb = .not.test_abpc
c
	if(test_acpo.neqv.test_acpb) then
c
		ireflex = ireflex + 1
c
		call find_tetra(itetra,2,a,c,o,ifind,
     1		tetra_ac,idx_a,idx_c)
c
		if(ifind.eq.1) then
			iflip = iflip + 1
			tetra_flip(iflip) = tetra_ac
			list_flip(iflip) = 2
			ktetra_touch(1,iflip)=tetra_neighbour(idx_a,
     1						tetra_ac)
			ktetra_touch(2,iflip)=tetra_neighbour(idx_c,
     1						tetra_ac)
			ktetra_idx(1,iflip)=tetra_nindex(idx_a,
     1						tetra_ac)
			ktetra_idx(2,iflip)=tetra_nindex(idx_c,
     1						tetra_ac)
			test_or(1,iflip) = .not.test_bcpo
			test_or(2,iflip) = test_abpo
		endif
c
	endif
c
c	Now check edge (bc): 
c		- (bc) is reflex iff o and a lies on opposite sides of
c		the hyperplane defined by (bcp). We therefore test the
c		orientation of (bcpo) and (bcpa): if they differ (bc)
c		is reflex
c		- if (bc) is reflex, we test if it is of degree 3.
c		(bc) is of degree 3 if it is shared by 3 tetrahedra,
c		namely (abcp), (abco) and (bcpo). The first two are itetra
c		and jtetra, so we only need to check if (bcpo) exists.
c		since (bcpo) contains p, (bcp) should then be a link facet
c		of p, so we test all tetrahedra that define link facets
c
	test_bcpa = test_abpc
c
	if(test_bcpo.neqv.test_bcpa) then
c
		ireflex = ireflex + 1
c
		call find_tetra(itetra,1,b,c,o,ifind,
     1		tetra_bc,idx_b,idx_c)
c
		if(ifind.eq.1) then
			iflip = iflip + 1
			tetra_flip(iflip) = tetra_bc
			list_flip(iflip) = 3
			ktetra_touch(1,iflip)=tetra_neighbour(idx_b,
     1						tetra_bc)
			ktetra_touch(2,iflip)=tetra_neighbour(idx_c,
     1						tetra_bc)
			ktetra_idx(1,iflip)=tetra_nindex(idx_b,
     1						tetra_bc)
			ktetra_idx(2,iflip)=tetra_nindex(idx_c,
     1						tetra_bc)
			test_or(1,iflip) = test_capo
			test_or(2,iflip) = .not.test_abpo
		endif
c
	endif
c
	if(ireflex.ne.iflip) goto 100
c
	if(iflip.eq.1) then
c
c		Only one edge is "flippable": we do a 3-2 flip
c
		iorder = list_flip(iflip)
		ia = table32(1,iorder)
		ib = table32(2,iorder)
		ic = table32(3,iorder)
		vert_flip(ia) = a
		vert_flip(ib) = b
		vert_flip(ic) = c
		vert_flip(4) = p
		vert_flip(5) = o
		ia = table32_2(1,iorder)
		ib = table32_2(2,iorder)
		itouch(1) = itetra_touch(ia)
		itouch(2) = itetra_touch(ib)
		i_idx(1)  = itetra_idx(ia)
		i_idx(2)  = itetra_idx(ib)
		jtouch(1) = jtetra_touch(ia)
		jtouch(2) = jtetra_touch(ib)
		j_idx(1)  = jtetra_idx(ia)
		j_idx(2)  = jtetra_idx(ib)
		ktouch(1) = ktetra_touch(1,iflip)
		ktouch(2) = ktetra_touch(2,iflip)
		k_idx(1)  = ktetra_idx(1,iflip)
		k_idx(2)  = ktetra_idx(2,iflip)
		call flipjw_3_2(itetra,jtetra,tetra_flip(1),vert_flip,
     1		itouch,i_idx,jtouch,j_idx,ktouch,k_idx,test_or(1,iflip),
     2		test_or(2,iflip),ierr)
c
	elseif(iflip.eq.2) then
c
c		In this case, one point is redundant: the point common to
c		the two edges that can be flipped. We then perform a 4-1
c		flip
c
		iorder = list_flip(1) + list_flip(2) - 2
		vert_flip(table41(1,iorder)) = a
		vert_flip(table41(2,iorder)) = b
		vert_flip(table41(3,iorder)) = c
		vert_flip(4) = p
		vert_flip(5) = o
		ii = table41_2(1,iorder)
		ij = table41_2(2,iorder)
		ishare = itetra_touch(iorder)
		idxi   = itetra_idx(iorder)
		jshare = jtetra_touch(iorder)
		idxj   = jtetra_idx(iorder)
		kshare = ktetra_touch(ii,1)
		idxk   = ktetra_idx(ii,1)
		lshare = ktetra_touch(ij,2)
		idxl   = ktetra_idx(ij,2)
		if(iorder.eq.1) then
			test = test_bcpo
		elseif(iorder.eq.2) then
			test = .not.test_capo
		else
			test = test_abpo
		endif
		call flipjw_4_1(itetra,jtetra,tetra_flip(1),
     1		tetra_flip(2),vert_flip,ishare,idxi,jshare,
     2		idxj,kshare,idxk,lshare,idxl,test,ierr)
c
	else
c
c	This case should not occur...
c
		write(6,*) 'Problem...three edges flippable!!'
c
	endif
c
	goto 100
c
200	continue
c
c	Add all "killed" tetrahedra in the free zone
c
	ns = nfree
	nfree = min(ns+nkill,nfreemax)
	do 300 j = ns+1,nfree
		free(j) = kill(j-ns)
300	continue
	nkill = 0
c
	return
	end
c	Flipjw_1_4.f		Version 1 12/17/2001	Patrice Koehl
c
c	This subroutine implements a 4->1 flip in 3D for regular triangulation
c
c	a 4->1 flip is a transformation in which a tetrahedron and a single
c	vertex included in the tetrahedron are transformed to 4 tetrahedra,
c	defined from the 4 four faces of the initial tetrahedron, connected
c	to the new point. Each of the faces are then called "link facet",
c	and stored on a queue
c
c	This version of flip_1_4 does not save the old tetrahedron
c	(i.e. no history dag) in order to save space. As a consequence,
c	it cannot be used with a point location scheme that uses the
c	history dag
c
	subroutine flipjw_1_4(ipoint,itetra)
c
c	Input:
c	******
c		- ipoint:	index of the point p to be included
c		- itetra:	index of the tetrahedra (a,b,c,d) considered
c
c	Output:
c	********
c		- nlink_facet:	4
c		- link_facet:	Add the four faces of the initial tetrahedron
c		- link_index:	A link_facet is a triangle defined from its
c				two neighboring tetrahedra. I store the position
c				of the vertex opposite to the triangle in each
c				tetrehedron in the array link_index
c
c	Include array dimensions
c
	integer	ntetra_max,nfreemax,nfacet_max
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter(MAX_FACET = 10000)


	parameter	(ntetra_max=MAX_TETRA)
	parameter	(nfreemax = MAX_FREE)
	parameter	(nfacet_max=MAX_FACET)
c
c	Define variables
c
	integer	i,j,k,ntetra,newtetra,n_new
	integer	ipoint,itetra,jtetra,nlink_facet
	integer	fact,idx
	integer	nfree,nkill
c
	integer	idx_list(3,4)
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
	integer	list_new(ntetra_max)
	integer	link_facet(2,nfacet_max)
	integer	link_index(2,nfacet_max)
	integer	tetra(4,ntetra_max)
	integer	tetra_neighbour(4,ntetra_max)
	integer	free(nfreemax),kill(nfreemax)
	integer	vertex(4),neighbour(4),nindex(4)
	integer	position(4)
c
	common  /tetra_zone/	ntetra,tetra,tetra_neighbour
	common  /tetra_stat/	tetra_status,tetra_orient,
     1				tetra_nindex
	common  /link_zone/	nlink_facet,link_facet,link_index
	common  /freespace/	nfree,nkill,free,kill
	common  /update/	n_new,list_new
c
	save
c
	data idx_list /1,1,1,1,2,2,2,2,3,3,3,3/
c
c	Store information about "old" tetrahedron
c
	do 50 i = 1,4
		vertex(i)    = tetra(i,itetra)
		neighbour(i) = tetra_neighbour(i,itetra)
		nindex(i)    = tetra_nindex(i,itetra)
50	continue
c
	fact = tetra_orient(itetra)
c
c	The four new tetrahedra are going to be stored
c	in : any free space in the tetrahedron list,
c	and at the end of the list of known tetrahedra
c
	k = 0
c
	do 100 i = nfree,max(nfree-3,1),-1
		k = k + 1
		position(k) = free(i)
100	continue
	nfree = max(nfree-4,0)
c
	do 150 i = k+1,4
		ntetra = ntetra + 1
		position(i) = ntetra
150	continue
c
c	itetra is set to 0, and added to the "kill" list
c
	tetra_status(itetra) = 0
	nkill = 1
	kill(nkill) = itetra
c
c	The tetrahedron is defined as (ijkl); four new tetrahedra are
c	created:	jklp, iklp, ijlp, and ijkp, where p is the new
c	point to be included
c
c	For each new tetrahedron, define all four neighbours:
c	For each neighbour, I store the index of the vertex opposite to 
c	the common face in array tetra_nindex
c
c	tetrahedron jklp : neighbours are iklp, ijlp, ijkp and neighbour
c			   of (ijkl) on face jkl
c	tetrahedron iklp : neighbours are jklp, ijlp, ijkp and neighbour
c			   of (ijkl) on face ikl
c	tetrahedron ijlp : neighbours are jklp, iklp, ijkp and neighbour
c			   of (ijkl) on face ijl
c	tetrahedron ijkp : neighbours are jklp, iklp, ijlp and neighbour
c			   of (ijkl) on face ijk
c
	do 250 i = 1,4
c
		newtetra = position(i)
		n_new = n_new + 1
		list_new(n_new) = newtetra
c
		k = 0
		do 200 j = 1,4
			if(j.eq.i) goto 200
			k = k+1
			tetra(k,newtetra)=vertex(j)
			tetra_neighbour(k,newtetra) = position(j)
			tetra_nindex(k,newtetra) = idx_list(k,i) 
200		continue
c
		jtetra = neighbour(i)
		idx = nindex(i)
		tetra(4,newtetra) = ipoint
		tetra_neighbour(4,newtetra) = jtetra
		tetra_nindex(4,newtetra) = idx
c
c		I must update the neighbors of the neighbour of itetra!
c
		if(jtetra.ne.0.and.idx.ne.0) then
			tetra_neighbour(idx,jtetra) = newtetra
			tetra_nindex(idx,jtetra) = 4
		endif
c
		tetra_status(newtetra) = 1
c
c		I also store the orientation of the tetrahedron
c		First, (jklp) and (ijlp) are cw, while (iklp) and (ijkp)
c		are ccw. 
c
		fact = -fact
		tetra_orient(newtetra) = fact
c
250	continue
c
c	Now add all fours faces of itetra in the link_facet queue.
c	Each link_facet (a triangle) is implicitly defined as the
c	intersection of two tetrahedra
c
c	link_facet:	jkl	tetrahedra:	jklp and neighbour of (ijkl)
c						on jkl
c	link_facet:	ikl	tetrahedra:	iklp and neighbour of (ijkl)
c						on ikl
c	link_facet:	ijl	tetrahedra:	ijlp and neighbour of (ijkl)
c						on ijl
c	link_facet:	ijk	tetrahedra:	ijkp and neighbour of (ijkl)
c						on ijk
c
	nlink_facet = 0
c
	do 300 i = 1,4
		newtetra = position(i)
		nlink_facet = nlink_facet + 1
		link_facet(1,nlink_facet) = newtetra
		link_facet(2,nlink_facet) = tetra_neighbour(4,newtetra)
		link_index(1,nlink_facet) = 4
		link_index(2,nlink_facet) = tetra_nindex(4,newtetra)
300	continue
c
c	end of flip_1_4
c
	return
	end
c	Flipjw_2_3.f		Version 1 12/17/2001	Patrice Koehl
c
c	This subroutine implements a 2->3 flip in 3D for regular triangulation
c
c	a 2->3 flip is a transformation in which two tetrahedrons are
c	flipped into three tetrahedra. The two tetrahedra (abcp) and
c	(abco) shares a triangle (abc) which is in the link_facet of the
c	current point p added to the triangulation. 
c	This flip is only possible if the union of the two tetrahedron is
c	convex, and if their shared triangle is not locally regular. 
c	We assume here that these tests have been performed and are
c	true. 
c	Once the flip has been performed, three new tetrahedra are added
c	and three new "link facet" are added to the link
c	facet queue
c
c	This version of flip_2_3 does not save the old tetrahedra
c	(i.e. no history dag) in order to save space. As a consequence,
c	it cannot be used with a point location scheme that uses a
c	history dag
c	
c
	subroutine flipjw_2_3(itetra,jtetra,p,o,itetra_touch,
     1	itetra_idx,jtetra_touch,jtetra_idx,test_abpo,test_bcpo,
     2	test_capo,ierr)
c
c	Input:
c	******
c		- itetra:	index of the tetrahedra (a,b,c,p) considered
c		- jtetra:	index of the tetrahedra (a,b,c,o) considered
c		- itetra_touch: the three tetrahedra that touches itetra on the
c				faces opposite to the 3 vertices a,b,c
c		- itetra_idx:	for the three tetrahedra defined by itetra_touch,
c				index of the vertex opposite to the face
c				common with itetra
c		- jtetra_touch: the three tetrahedra that touches jtetra on the
c				faces opposite to the 3 vertices a,b,c
c		- jtetra_idx:	for the three tetrahedra defined by jtetra_touch,
c				index of the vertex opposite to the face
c				common with jtetra
c		- test_abpo:	orientation of the four points a,b,p,o
c		- test_bcpo:	orientation of the four points b,c,p,o
c		- test_capo:	orientation of the four points c,a,p,o
c
c	Output:
c	********
c		- nlink_facet:	3 new link facets are added
c		- link_facet:	the three faces of the initial tetrahedron
c				(a,b,c,o) containing the vertex o are added
c				as link facets
c		- link_index:   A link_facet is a triangle defined from its
c				two neighbouring tetrahedra. I store the position
c				of the vertex opposite to the triangle in each
c				tetrehedron in the array link_index
c		- ierr:		1 if flip was not possible
c
c	Include array dimensions
c
	integer	ntetra_max,nfreemax,nfacet_max
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter(MAX_FACET = 10000)

	parameter	(ntetra_max=MAX_TETRA)
	parameter	(nfreemax = MAX_FREE)
	parameter	(nfacet_max=MAX_FACET)
c
c	Define variables
c
	integer	i,j,k,p,o
	integer	ierr,ntetra,nlink_facet
	integer	itetra,jtetra
	integer	it,jt,idx,jdx,pos,opos,o_place
	integer	nfree,nkill,newtetra,n_new
c
	integer	jtetra_touch(3),itetra_touch(3)
	integer	jtetra_idx(3),itetra_idx(3)
	integer	idx_list(2,3,4)
	integer	other_pos(2,3,4),o_pos(3,4)
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
	integer	list_new(ntetra_max)
	integer	link_facet(2,nfacet_max)
	integer	link_index(2,nfacet_max)
	integer	tetra(4,ntetra_max)
	integer	tetra_neighbour(4,ntetra_max)
	integer tests(3),position(3)
	integer	free(nfreemax),kill(nfreemax)
c
	logical test_abpo,test_bcpo,test_capo
c
	common  /tetra_zone/	ntetra,tetra,tetra_neighbour
	common  /tetra_stat/	tetra_status,tetra_orient,
     1				tetra_nindex
	common  /link_zone/	nlink_facet,link_facet,link_index
	common  /freespace/	nfree,nkill,free,kill
	common  /update/	n_new,list_new
c
	save
c
	data idx_list /2,2,2,3,3,3,1,1,2,3,3,3,1,1,1,2,3,3,1,1,1,2,2,2/
	data other_pos/2,3,2,3,2,3,2,3,1,3,1,3,1,3,1,3,1,2,1,2,1,2,1,2/
	data o_pos/1,1,1,1,2,2,2,2,3,3,3,3/
c
	ierr = 0
c
c	If itetra or jtetra are inactive, cannot flip
c
	if(tetra_status(itetra).ne.1.or.tetra_status(jtetra).ne.1) then
		ierr = 1
		return
	endif
c
c	The three new tetrahedra are going to be stored
c       in : any free space in the tetrahedron list,
c       and at the end of the list of known tetrahedra if needed
c
	k = 0
	do 50 i = nfree,max(nfree-2,1),-1
		k = k + 1
		position(k) = free(i)
50	continue
	nfree = max(nfree-3,0)
c
	do 100 i = k+1,3
		ntetra = ntetra + 1
		position(i) = ntetra
100	continue
c
c	Set itetra and jtetra to 0, and add them to kill list
c
	tetra_status(itetra) = 0
	tetra_status(jtetra) = 0
	kill(nkill+1) = itetra
	kill(nkill+2) = jtetra
	nkill = nkill + 2
c
c	I need :
c		- the vertices a,b,c are the first vertices of itetra
c		- the vertices p and o
c		- for each vertex in the triangle, define the opposing
c		faces in the two tetrahedra itetra and jtetra, and
c		the tetrahedra that share that faces with itetra and
c		jtetra, respectively. This information is stored
c		in two arrays, itetra_touch and jtetra_touch
c
c	These information are provided by the calling program
c
c	For bookkeeping reasons, I always store p as the last vertex
c
c	Now I define the three new tetrahedra: (bcop), (acop) and (abop)
c	as well as their neighbours
c
c	tetrahedron bcop : neighbours are acop, abop, neighbour of (abcp)
c			   on face bcp, and neighbour of (abco) on face bco
c	tetrahedron acop : neighbours are bcop, abop, neighbour of (abcp)
c			   on face acp, and neighbour of (abco) on face aco
c	tetrahedron abop : neighbours are bcop, acop, neighbour of (abcp)
c			   on face abp, and neighbour of (abco) on face abo
c
c
	tests(1) = 1
	if(test_bcpo) tests(1) = -1
	tests(2) = -1
	if(test_capo) tests(2) = 1
	tests(3) = 1
	if(test_abpo) tests(3) = -1
c
c	First look at the position of a with respect to a,b,c:
c
	if(o.lt.tetra(1,itetra)) then
		opos = 1
	elseif(o.lt.tetra(2,itetra)) then
		opos = 2
		tests(2) = -tests(2)
		tests(3) = -tests(3)
	elseif(o.lt.tetra(3,itetra)) then
		opos = 3
		tests(1) = -tests(1)
		tests(2) = -tests(2)
	else
		opos = 4
	endif
c		
	do 200 i = 1,3
c
		newtetra = position(i)
		n_new = n_new + 1
		list_new(n_new) = newtetra
c
		k = 0
		do 150 j = 1,3
			if(j.eq.i) goto 150
			k = k+1
			pos = other_pos(k,i,opos)
			tetra(pos,newtetra)=tetra(j,itetra)
			tetra_neighbour(pos,newtetra) = position(j)
			tetra_nindex(pos,newtetra) = idx_list(k,i,opos)
150		continue
c
		o_place = o_pos(i,opos)
		tetra(o_place,newtetra) = o
		it = itetra_touch(i)
		idx = itetra_idx(i)
		tetra_neighbour(o_place,newtetra) = it
		tetra_nindex(o_place,newtetra) = idx
		if(idx.ne.0.and.it.ne.0) then
			tetra_neighbour(idx,it) = newtetra
			tetra_nindex(idx,it) = o_place
		endif
c
		tetra(4,newtetra) = p
		jt = jtetra_touch(i)
		jdx = jtetra_idx(i)
		tetra_neighbour(4,newtetra) = jt
		tetra_nindex(4,newtetra) = jdx
		if(jdx.ne.0.and.jt.ne.0) then
			tetra_neighbour(jdx,jt) = newtetra
			tetra_nindex(jdx,jt) = 4
		endif
c
		tetra_status(newtetra) = 1
c
		tetra_orient(newtetra) = tests(i)
c
200	continue
c
c	if(opos.le.4) then
c		write(6,*) ' '
c		write(6,*) 'opos =',opos
c		write(6,*) 'New tetra 1'
c		write(6,*) 'Vertices :',(tetra(j,position(1)),j=1,4)
c		write(6,*) 'neighbours:',(tetra_neighbour(j,position(1)),j=1,4)
c		write(6,*) 'nindices:',(tetra_nindex(j,position(1)),j=1,4)
c		write(6,*) 'The four neighbours:'
c		write(6,*) (tetra(j,tetra_neighbour(1,position(1))),j=1,4)
c		write(6,*) (tetra(j,tetra_neighbour(2,position(1))),j=1,4)
c		write(6,*) (tetra(j,tetra_neighbour(3,position(1))),j=1,4)
c		if(tetra_neighbour(4,position(1)).ne.0) then
c		write(6,*) (tetra(j,tetra_neighbour(4,position(1))),j=1,4)
c		endif
c		write(6,*) ' '
c		write(6,*) 'New tetra 2'
c		write(6,*) 'Vertices :',(tetra(j,position(2)),j=1,4)
c		write(6,*) 'neighbours:',(tetra_neighbour(j,position(2)),j=1,4)
c		write(6,*) 'nindices:',(tetra_nindex(j,position(2)),j=1,4)
c		write(6,*) 'The four neighbours:'
c		write(6,*) (tetra(j,tetra_neighbour(1,position(2))),j=1,4)
c		write(6,*) (tetra(j,tetra_neighbour(2,position(2))),j=1,4)
c		write(6,*) (tetra(j,tetra_neighbour(3,position(2))),j=1,4)
c		if(tetra_neighbour(4,position(2)).ne.0) then
c		write(6,*) (tetra(j,tetra_neighbour(4,position(2))),j=1,4)
c		endif
c		write(6,*) ' '
c		write(6,*) 'New tetra 3'
c		write(6,*) 'Vertices :',(tetra(j,position(3)),j=1,4)
c		write(6,*) 'neighbours:',(tetra_neighbour(j,position(3)),j=1,4)
c		write(6,*) 'nindices:',(tetra_nindex(j,position(3)),j=1,4)
c		write(6,*) 'The four neighbours:'
c		write(6,*) (tetra(j,tetra_neighbour(1,position(3))),j=1,4)
c		write(6,*) (tetra(j,tetra_neighbour(2,position(3))),j=1,4)
c		write(6,*) (tetra(j,tetra_neighbour(3,position(3))),j=1,4)
c		if(tetra_neighbour(4,position(3)).ne.0) then
c		write(6,*) (tetra(j,tetra_neighbour(4,position(3))),j=1,4)
c		endif
c		write(6,*) ' '
cc		stop
c	endif
c
c	Now add all three faces of jtetra containing o in the link_facet queue.
c	Each link_facet (a triangle) is implicitly defined as the
c	intersection of two tetrahedra
c
c	link_facet:	bco	tetrahedra:	bcop and neighbour of (abco)
c						on bco
c	link_facet:	aco	tetrahedra:	acop and neighbour of (abco)
c						on aco
c	link_facet:	abo	tetrahedra:	abop and neighbour of (abco)
c						on abo
c
	do 250 i = 1,3
		newtetra = position(i)
		nlink_facet = nlink_facet + 1
		link_facet(1,nlink_facet) = newtetra
		link_facet(2,nlink_facet) = tetra_neighbour(4,newtetra)
		link_index(1,nlink_facet) = 4
		link_index(2,nlink_facet) = tetra_nindex(4,newtetra)
250	continue
c
	return
	end
c	Flipjw_3_2.f		Version 1 12/17/2001	Patrice Koehl
c
c	This subroutine implements a 3->2 flip in 3D for regular triangulation
c
c	a 3->2 flip is a transformation in which three tetrahedrons are
c	flipped into two tetrahedra. The two tetrahedra (abpo), (abcp) and
c	(abco) shares an edge (ab) which is in the link_facet of the
c	current point p added to the triangulation. 
c	This flip is only possible if the edge ab is reflex, with degree 3
c	We assume here that these tests have been performed and are
c	true. 
c	Once the flip has been performed, two new tetrahedra are added
c	and two new "link facet" are added to the link facet queue
c
c	This version of flip_3_2 does not save the old tetrahedron
c	(i.e. no history dag) in order to save space. As a consequence,
c	it cannot be used with a point location scheme that uses the
c	history dag
c
	subroutine flipjw_3_2(itetra,jtetra,ktetra,vertices,
     1			itetra_touch,itetra_idx,jtetra_touch,
     2			jtetra_idx,ktetra_touch,ktetra_idx,
     3			test_bcpo,test_acpo,ierr)
c
c	Input:
c	******
c		- itetra:	index of the tetrahedra (a,b,c,p) considered
c		- jtetra:	index of the tetrahedra (a,b,c,o) considered
c		- ktetra:	index of the tetrahedra (a,b,o,p) considered
c		- vertices:	the five vertices a,b,c,p,o
c		- itetra_touch:	indices of the two tetrahedra that share the
c				faces opposite to a and b in itetra,
c				respectively
c		- itetra_idx:   for the two tetrahedra defined by itetra_touch,
c				index position of the vertex opposite to the face
c				common with itetra
c		- jtetra_touch:	indices of the two tetrahedra that share the
c				faces opposite to a and b in jtetra,
c				respectively
c		- jtetra_idx:   for the two tetrahedra defined by jtetra_touch,
c				index position of the vertex opposite to the face
c				common with jtetra
c		- ktetra_touch:	indices of the two tetrahedra that share the
c				faces opposite to a and b in ktetra,
c				respectively
c		- ktetra_idx:   for the two tetrahedra defined by ktetra_touch,
c				index position of the vertex opposite to the face
c				common with ktetra
c		- test_bcpo   : orientation of the four points b,c,p,o
c		- test_acpo   : orientation of the four points a,c,p,o
c
c	Output:
c	********
c		- nlink_facet:	2 new link facets are added
c		- link_facet:	the two faces of the initial tetrahedron
c				(a,b,o,p) containing the edge op are added
c				as link facets
c		- link_index:   A link_facet is a triangle defined from its
c				two neighbouring tetrahedra. I store the position
c				of the vertex opposite to the triangle in each
c				tetrehedron in the array link_index
c		- ierr:		1 if flip was not possible
c
c	Include array dimensions
c
	integer	ntetra_max,nfreemax,nfacet_max
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter(MAX_FACET = 10000)


	parameter	(ntetra_max=MAX_TETRA)
	parameter	(nfreemax = MAX_FREE)
	parameter	(nfacet_max=MAX_FACET)
c
c	Define variables
c
	integer	i,j,k,p,o,c
	integer	ierr,ntetra,nlink_facet,n_new
	integer	itetra,jtetra,ktetra
	integer	it,jt,kt,idx,jdx,kdx
	integer	nfree,nkill,newtetra,nswap
c
	integer	edge(2),tests(2),vertices(5)
	integer	itetra_touch(2),jtetra_touch(2),ktetra_touch(2)
	integer	itetra_idx(2),jtetra_idx(2),ktetra_idx(2)
	integer	free(nfreemax),kill(nfreemax),position(2)
	integer vertex(3),neighbour(3),nindex(3),index(3)
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
	integer	list_new(ntetra_max)
	integer	link_facet(2,nfacet_max)
	integer	link_index(2,nfacet_max)
	integer	tetra(4,ntetra_max)
	integer	tetra_neighbour(4,ntetra_max)
c
	logical test_bcpo,test_acpo
c
	common  /tetra_zone/	ntetra,tetra,tetra_neighbour
	common /tetra_stat/	tetra_status,tetra_orient,
     1				tetra_nindex
	common /link_zone/	nlink_facet,link_facet,link_index
	common /freespace/	nfree,nkill,free,kill
	common /update/		n_new,list_new
c
	save
c
	tests(1) = 1
	if(test_bcpo) tests(1) = -1
	tests(2) = 1
	if(test_acpo) tests(2) = -1
c
	ierr = 0
c
c	If itetra, jtetra or ktetra are inactive, cannot flip
c
	if(tetra_status(itetra).ne.1.or.tetra_status(jtetra).ne.1
     1	.or.tetra_status(ktetra).ne.1) then
		ierr = 1
		return
	endif
c
	edge(1) = vertices(1)
	edge(2) = vertices(2)
	c       = vertices(3)
	p       = vertices(4)
	o       = vertices(5)
c
c	The two new tetrahedra are going to be stored "free" space, or 
c	at the end of the list
c
	k = 0
	do 50 i = nfree,max(nfree-1,1),-1
		k = k + 1
		position(k) = free(i)
50	continue
	nfree = max(nfree-2,0)
c
	do 100 i = k+1,2
		ntetra = ntetra + 1
		position(i) = ntetra
100	continue
c
c	itetra, jtetra and ktetra becomes "available"; they are added to the
c	"kill" list
c
	tetra_status(itetra) = 0
	tetra_status(jtetra) = 0
	tetra_status(ktetra) = 0
c
	kill(nkill+1) = itetra
	kill(nkill+2) = jtetra
	kill(nkill+3) = ktetra
	nkill = nkill + 3
c
c	I need :
c		- the two vertices that define their common edge (ab)
c		 these vertices are stored in the array edge
c		- the vertices c, p and o that form the new triangle
c		- for each vertex in the edge (ab), define the opposing
c		faces in the three tetrahedra itetra, jtetra and ktetra, and
c		the tetrahedron that share these faces with itetra, jtetra and
c		ktetra, respectively. This information is stored
c		in three arrays, itetra_touch, jtetra_touch and ktetra_touch
c
c	These information are given by the calling program
c
c	For bookkeeping reasons, I always set p to be the last vertex
c	of the new tetrahedra
c
c	Now I define the two new tetrahedra: (bcop) and (acop)
c	as well as their neighbours
c
c	tetrahedron bcop : neighbours are acop, neighbour of (abop)
c			   on face bpo, neighbour of (abcp) on face bcp
c			   and neighbour of (abco) on face (bco)
c	tetrahedron acop : neighbours are bcop, neighbour of (abop)
c			   on face apo, neighbour of (abcp) on face acp
c			   and neighbour of (abco) on face (aco)
c
	do 200 i = 1,2
c
		newtetra = position(i)
		n_new = n_new + 1
		list_new(n_new) = newtetra
c
		k = 0
		do 150 j = 1,2
			if(j.eq.i) goto 150
			k = k+1
			tetra(k,newtetra)=edge(j)
			tetra_neighbour(k,newtetra) = position(j)
			tetra_nindex(k,newtetra) = 1
150		continue
c
		tetra(2,newtetra) = c
		kt = ktetra_touch(i)
		kdx = ktetra_idx(i)
		tetra_neighbour(2,newtetra) = kt
		tetra_nindex(2,newtetra) = kdx
		if(kdx.ne.0.and.kt.ne.0) then
			tetra_neighbour(kdx,kt) = newtetra
			tetra_nindex(kdx,kt) = 2
		endif
c
		tetra(3,newtetra) = o
		it = itetra_touch(i)
		idx = itetra_idx(i)
		tetra_neighbour(3,newtetra) = it
		tetra_nindex(3,newtetra) = idx
		if(idx.ne.0.and.it.ne.0) then
			tetra_neighbour(idx,it) = newtetra
			tetra_nindex(idx,it) = 3
		endif
c
		tetra(4,newtetra) = p
		jt = jtetra_touch(i)
		jdx = jtetra_idx(i)
		tetra_neighbour(4,newtetra) = jt
		tetra_nindex(4,newtetra) = jdx
		if(jdx.ne.0.and.jt.ne.0) then
			tetra_neighbour(jdx,jt) = newtetra
			tetra_nindex(jdx,jt) = 4
		endif
c
		tetra_status(newtetra) = 1
c
		tetra_orient(newtetra) = tests(i)
c
200	continue
c
c	Now add the two faces of ktetra containing (co) in the link_facet 
c	queue.
c	Each link_facet (a triangle) is implicitly defined as the
c	intersection of two tetrahedra
c
c	link_facet:	bco	tetrahedra:	bcop and neighbour of (abco)
c						on bco
c	link_facet:	aco	tetrahedra:	acop and neighbour of (abco)
c						on aco
c
	do 600 k = 1,2
c
		i = position(k)
c
		do 300 j = 1,3
			vertex(j) = tetra(j,i)
300		continue
c
		call sort4_sign(vertex,index,nswap,3)
c
		do 400 j = 1,3
			neighbour(j)=tetra_neighbour(index(j),i)
			nindex(j)=tetra_nindex(index(j),i)
			if(neighbour(j).ne.0) then
				tetra_nindex(nindex(j),neighbour(j))=j
			endif
400		continue
c
		do 500 j = 1,3
			tetra(j,i)=vertex(j)
			tetra_neighbour(j,i)=neighbour(j)
			tetra_nindex(j,i) = nindex(j)
500		continue
c
		tetra_orient(i) = tetra_orient(i)*nswap
c
600	continue
c
	do 700 i = 1,2
		newtetra = position(i)
		nlink_facet = nlink_facet + 1
		link_facet(1,nlink_facet) = newtetra
		link_facet(2,nlink_facet) = tetra_neighbour(4,newtetra)
		link_index(1,nlink_facet) = 4
		link_index(2,nlink_facet) = tetra_nindex(4,newtetra)
700	continue
c
	return
	end
c	Flipjw_4_1.f		Version 1 12/17/2001	Patrice Koehl
c
c	This subroutine implements a 4->1 flip in 3D for regular triangulation
c
c	a 4->1 flip is a transformation in which four tetrahedra are
c	flipped into one tetrahedron. The two tetrahedra (abop), (bcop),
c	(abcp) and (abco) shares a vertex (b) which is in the link_facet of the
c	current point p added to the triangulation. After the flip, b
c	is set to redundant. 
c	This flip is only possible if the two edges (ab) and (bc)
c	are reflex of order 3.
c	We assume here that these tests have been performed and are
c	true. 
c	Once the flip has been performed, one tetrahedron is added
c	and one new "link facet" is added to the link facet queue
c
c       This version of flip_4_1 does not save the old tetrahedron
c       (i.e. no history dag) in order to save space. As a consequence,
c       it cannot be used with a point location scheme that uses the
c       history dag
c
	subroutine flipjw_4_1(itetra,jtetra,ktetra,ltetra,vertices,
     1			ishare,idx,jshare,jdx,kshare,kdx,lshare,ldx,
     2			test_acpo,ierr)
c
c	Input:
c	******
c		- itetra:	index of the tetrahedra (a,b,c,p) considered
c		- jtetra:	index of the tetrahedra (a,b,c,o) considered
c		- ktetra:	index of the tetrahedra (a,b,o,p) considered
c		- ltetra:	index of the tetrahedra (b,c,o,p) considered
c		- vertices:	index of a,b,c,p,o
c		- ishare:	index of tetrahedron sharing the face 
c				opposite to b in itetra
c		- idx		index of the vertex of ishare opposite to the
c				face of ishare shared with itetra
c		- jshare:	index of tetrahedron sharing the face 
c				opposite to b in jtetra
c		- jdx		index of the vertex of jshare opposite to the
c				face of jshare shared with jtetra
c		- kshare:	index of tetrahedron sharing the face 
c				opposite to b in ktetra
c		- kdx		index of the vertex of kshare opposite to the
c				face of kshare shared with ktetra
c		- lshare:	index of tetrahedron sharing the face 
c				opposite to b in ltetra
c		- ldx		index of the vertex of lshare opposite to the
c				face of lshare shared with ltetra
c		- test_acpo:	orientation of the 4 points (a,c,p,o)
c
c	Output:
c	********
c		- nlink_facet:	1 new link facet is added
c		- link_facet:	the face of the initial tetrahedron
c				(a,b,c,o) opposite to the vertex b is added
c				as link facet
c		- link_index:   A link_facet is a triangle defined from its
c				two neighbouring tetrahedra. I store the position
c				of the vertex opposite to the triangle in each
c				tetrehedron in the array link_index
c		- ierr:		1 if flip was not possible
c
c	Include array dimensions
c
	integer	ntetra_max,nfreemax,nfacet_max,npointmax
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter(MAX_FACET = 10000)

	parameter	(ntetra_max=MAX_TETRA)
	parameter	(nfreemax = MAX_FREE)
	parameter	(nfacet_max=MAX_FACET)
	parameter	(npointmax=MAX_POINT)
c
c	Define variables
c
	integer	i,j
	integer	p,o,a,b,c
	integer	ierr,npoints,nvertex,nlink_facet,n_new
	integer	ntetra,itetra,jtetra,ktetra,ltetra
	integer	ishare,jshare,kshare,lshare
	integer	idx,jdx,kdx,ldx
	integer	test1,newtetra
	integer	nfree,nkill,nswap
c
	integer	vertices(5)
c
	integer	redinfo(npointmax)
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
	integer	list_new(ntetra_max)
	integer	link_facet(2,nfacet_max)
	integer	link_index(2,nfacet_max)
	integer	tetra(4,ntetra_max)
	integer	tetra_neighbour(4,ntetra_max)
	integer	free(nfreemax),kill(nfreemax)
	integer vertex(3),neighbour(3),nindex(3),index(3)
c
	logical test_acpo
c
	common  /tetra_zone/	ntetra,tetra,tetra_neighbour
	common	/tetra_stat/	tetra_status,tetra_orient,
     1				tetra_nindex
	common /vertex_zone/	npoints,nvertex,redinfo
	common /link_zone/	nlink_facet,link_facet,link_index
	common /freespace/	nfree,nkill,free,kill
	common /update/		n_new,list_new
c
	save
c
	ierr = 0
c
	test1 = 1
	if(test_acpo) test1 = -1
c
c	If itetra, jtetra, ktetra, ltetra are inactive, cannot flip
c
	if(tetra_status(itetra).ne.1.or.tetra_status(jtetra).ne.1.
     1	or.tetra_status(ktetra).ne.1.or.tetra_status(ltetra).ne.1) then
		ierr = 1
		return
	endif
c
c	The new tetrahedron is going to be store in place of itetra
c
	if(nfree.ne.0) then
		newtetra = free(nfree)
		nfree = nfree -1
	else
		ntetra = ntetra + 1
		newtetra = ntetra
	endif
c
	n_new = n_new + 1
	list_new(n_new) = newtetra
c
c	jtetra, ktetra and ltetra become "available"; they
c	are added to the "kill" zone
c
	kill(nkill+1) = itetra
	kill(nkill+2) = jtetra
	kill(nkill+3) = ktetra
	kill(nkill+4) = ltetra
	nkill = nkill + 4
c
	tetra_status(itetra) = 0
	tetra_status(jtetra) = 0
	tetra_status(ktetra) = 0
	tetra_status(ltetra) = 0
c
c	I need :
c		- the vertex b that is shared by all 4 tetrahedra
c		- the vertices a, c, p and o
c		- for each tetrahedron, find neighbour attached to the face
c		oposite to b; this information is stored if *share,
c		where * can be i, j, k or l
c
c	These information are provided by the calling program
c
	a = vertices(1)
	b = vertices(2)
	c = vertices(3)
	p = vertices(4)
	o = vertices(5)
c
c	For bookkeeping reason, p is set to be the last vertex of the
c	new tetrahedron
c
c	Now I define the new tetrahedron: (acop)
c
c	tetrahedron acop : neighbor of (bcop) on face cpo, neighbor of (abop)
c			   on face apo, neighbor of (abcp) on face acp
c			   and neighbor of (abco) on face aco
c
	redinfo(b) = 1
c
	tetra(1,newtetra) = a
	tetra_neighbour(1,newtetra) = lshare
	tetra_nindex(1,newtetra) = ldx
	if(lshare.ne.0.and.ldx.ne.0) then
		tetra_neighbour(ldx,lshare) = newtetra
		tetra_nindex(ldx,lshare) = 1
	endif
c
	tetra(2,newtetra) = c
	tetra_neighbour(2,newtetra) = kshare
	tetra_nindex(2,newtetra) = kdx
	if(kshare.ne.0.and.kdx.ne.0) then
		tetra_neighbour(kdx,kshare) = newtetra
		tetra_nindex(kdx,kshare) = 2
	endif
c
	tetra(3,newtetra) = o
	tetra_neighbour(3,newtetra) = ishare
	tetra_nindex(3,newtetra) = idx
	if(ishare.ne.0.and.idx.ne.0) then
		tetra_neighbour(idx,ishare) = newtetra
		tetra_nindex(idx,ishare) = 3
	endif
c
	tetra(4,newtetra) = p
	tetra_neighbour(4,newtetra) = jshare
	tetra_nindex(4,newtetra) = jdx
	if(jshare.ne.0.and.jdx.ne.0) then
		tetra_neighbour(jdx,jshare) = newtetra
		tetra_nindex(jdx,jshare) = 4
	endif
c
	tetra_status(newtetra) = 1
c
	tetra_orient(newtetra) = test1
c
c	Now add one link facet : 
c
c	link_facet:	aco	tetrahedra:	acop and neighbour of (abco)
c						on aco
c
	nlink_facet = nlink_facet + 1
	link_facet(1,nlink_facet) = newtetra
	link_facet(2,nlink_facet) = jshare
	link_index(1,nlink_facet) = 4
	link_index(2,nlink_facet) = jdx
c
c	end of flip
c
	i = newtetra
c
	do 300 j = 1,3
		vertex(j) = tetra(j,i)
300	continue
c
	call sort4_sign(vertex,index,nswap,3)
c
	do 400 j = 1,3
		neighbour(j)=tetra_neighbour(index(j),i)
		nindex(j)=tetra_nindex(index(j),i)
		if(neighbour(j).ne.0) then
			tetra_nindex(nindex(j),neighbour(j))=j
		endif
400	continue
c
	do 500 j = 1,3
		tetra(j,i)=vertex(j)
		tetra_neighbour(j,i)=neighbour(j)
		tetra_nindex(j,i) = nindex(j)
500	continue
c
	tetra_orient(i) = tetra_orient(i)*nswap
c
	return
	end
c	Define_facet.f		Version 1 12/21/2001	Patrice Koehl
c
c	A triangle (or facet) is defined by the intersection of two 
c	tetrahedra itetra and jtetra
c	If we know the position of its three vertices in the first
c	tetrahedron (in fact the first three vertices a,b and c
c	of itetra), we need to find the indices of these vertices
c	in the second tetrahedron.
c	This routine also stores information about the neighbours
c	of the two tetrahedra considered
c
c	The vertices are called a,b,c,p, and o, where (abc) is the
c	common facet
c
	subroutine define_facet(itetra,jtetra,idx_o,
     1			itouch,idx,jtouch,jdx)
c
c	Input:
c	******
c		- itetra:	index of the tetrahedra (a,b,c,p) considered
c		- jtetra:	index of the tetrahedra (a,b,c,o) considered
c		- idx_o:	position of o in the vertices of jtetra
c
c	Output:
c	********
c		- itouch	itouch(i) is the tetrahedron sharing
c				the face opposite to i in tetrahedron itetra
c		- idx		idx(i) is the vertex of itouch(i) opposite
c				to the face shared with itetra
c		- jtouch	jtouch(i) is the tetrahedron sharing
c				the face opposite to i in tetrahedron jtetra
c		- jdx		jdx(i) is the vertex of jtouch(i) opposite
c				to the face shared with jtetra
c
c	Include array dimensions
c
	integer	ntetra_max
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter(MAX_FACET = 10000)


	parameter	(ntetra_max=MAX_TETRA)
c
c	Define variables
c
	integer	i,k
	integer	ia,ib,ie,if
	integer	ntetra,itetra,jtetra
	integer	idx_o
c
	integer	other(3,4),other2(2,4,4)
c
	integer corresp(3)
	integer	itouch(3),jtouch(3)
	integer	idx(3),jdx(3)
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
	integer	tetra(4,ntetra_max)
	integer	tetra_neighbour(4,ntetra_max)
c
	common  /tetra_zone/	ntetra,tetra,tetra_neighbour
	common  /tetra_stat/	tetra_status,tetra_orient,
     1				tetra_nindex
c
	save
c
	data other/2,3,4,
     1		   1,3,4,
     2		   1,2,4,
     3		   1,2,3/
c
	data other2/0,0,3,4,2,4,2,3,3,4,0,0,1,4,1,3,2,4,1,4,0,0,1,2,
     1			2,3,1,3,1,2,0,0/
c
c	I need to :
c		- find the three vertices that define their common face
c		 these vertices are stored in the array triangle
c		- find the vertices p and o
c
c	To define the common face of the two tetrahedra itetra and jtetra,
c	I look at the neighbours of itetra : one of them is jtetra!
c	This also provides p. The same procedure is repeated for jtetra,
c	to get o
c
	do 100 i = 1,3
		itouch(i) = tetra_neighbour(i,itetra)
		idx(i) = tetra_nindex(i,itetra)
100	continue
c
	ia = tetra(1,itetra)
	do 200 i = 1,3
		k = other(i,idx_o)
		ie = tetra(k,jtetra)
		if(ia.eq.ie) then
			corresp(1) = k
			goto 300
		endif
200	continue
c
300	continue
c
	ib = tetra(2,itetra)
	ie = other2(1,corresp(1),idx_o)
	if = other2(2,corresp(1),idx_o)
	if(ib.eq.tetra(ie,jtetra)) then
		corresp(2) = ie
		corresp(3) = if
	else
		corresp(2) = if
		corresp(3) = ie
	endif
c
	do 400 i = 1,3
		k = corresp(i)
		jtouch(i) = tetra_neighbour(k,jtetra)
		jdx(i) = tetra_nindex(k,jtetra)
400	continue
c
	return
	end
c	Find_tetra.f		Version 1 1/8/2002	Patrice Koehl
c
c	This subroutine tests if four given points form an existing
c	tetrahedron
c
	subroutine find_tetra(itetra,idx_c,a,b,o,ifind,tetra_loc,
     1		idx_a,idx_b)
c
c	Input:
c		itetra:		index of tetrahedra (abcp)
c		idx_c:		index of c in (abcp)
c		o:		index of vertex o
c
c	Output:
c		ifind:		1 if tetrahedron exists, 0 otherwise
c		tetra_loc:	index of existing tetrahedron, if it exists
c
c	We are testing in tetrahedron (abpo) exists. If it exists, it is
c	a neighbour of abcp, on the face opposite to vertex c.
c	We test that tetrahedron an see if it contains o
c
	integer	ntetra_max
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter(MAX_FACET = 10000)


	parameter	(ntetra_max=MAX_TETRA)
c
	integer	i,ifind,itetra,tetra_loc
	integer	ot,otx,otest
	integer	idx_c,idx_a,idx_b,o,a,b
	integer	ntetra
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
	integer tetra(4,ntetra_max),tetra_neighbour(4,ntetra_max)
c
	common  /tetra_zone/	ntetra,tetra,tetra_neighbour
	common  /tetra_stat/    tetra_status,tetra_orient,
     1				tetra_nindex
c
	save
c
	ot = tetra_neighbour(idx_c,itetra)
	otx = tetra_nindex(idx_c,itetra)
c
	otest = tetra(otx,ot)
c
	if(otest.eq.o) then
		ifind = 1
		tetra_loc = ot
c
c		We found the tetrahedron, let us define the position
c		of a and b in this tetrahedron
c
		do 100 i = 1,4
			if(tetra(i,tetra_loc).eq.a) then
				idx_a = i
			elseif(tetra(i,tetra_loc).eq.b) then
				idx_b = i
			endif
100		continue
c		
	else
		ifind = 0
	endif
c
	return
	end
c	Remove_inf.f		Version 1 1/16/2002	Patrice Koehl
c
c	This subroutine sets to 0 the status of tetrahedron that
c	contains infinite points
c
	subroutine remove_inf
c
	integer	ntetra_max,npointmax
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter(MAX_FACET = 10000)


	parameter	(npointmax=MAX_POINT)
	parameter	(ntetra_max=MAX_TETRA)
c
	integer	i,a,b,c,d
	integer	ntetra,ninf
	integer	npoint,nvertex
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 infpoint(npointmax)
	integer*1 tetra_nindex(4,ntetra_max)
c
	integer	redinfo(npointmax)
	integer tetra(4,ntetra_max)
	integer tetra_neighbour(4,ntetra_max)
c
	common	/tetra_zone/	ntetra,tetra,tetra_neighbour
	common	/tetra_stat/	tetra_status,tetra_orient,
     1				tetra_nindex
	common  /infinite/	infpoint
	common  /vertex_zone/	npoint,nvertex,redinfo
c
	save
c
	do 100 i = 1,ntetra
c
		if(tetra_status(i).ne.1) goto 100
c
		a = tetra(1,i)
		b = tetra(2,i)
		c = tetra(3,i)
		d = tetra(4,i)
c
		ninf = infpoint(a)+infpoint(b)+infpoint(c)+infpoint(d)
c
		if(ninf.ne.0) then
			tetra_status(i) = -3
			if(infpoint(a).eq.1) call mark_zero(i,1)
			if(infpoint(b).eq.1) call mark_zero(i,2)
			if(infpoint(c).eq.1) call mark_zero(i,3)
			if(infpoint(d).eq.1) call mark_zero(i,4)
		endif
c
100	continue
c
	do 200 i = 1,4
		redinfo(i) = 1
200	continue
c
	return
	end
c
c	Mark_zero.f	Version 1 3/11/2002	Patrice Koehl
c
c	This subroutine marks the tetrahedron that touches
c	a tetrahedron with infinite point as part of the
c	convex hull (i.e. one of its neighbor is 0)
c
	subroutine mark_zero(itetra,ivertex)
c
	integer	ntetra_max
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter(MAX_FACET = 10000)


	parameter	(ntetra_max=MAX_TETRA)
c
	integer	ntetra
	integer	itetra,ivertex
	integer	jtetra,jvertex
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
	integer tetra(4,ntetra_max)
	integer tetra_neighbour(4,ntetra_max)
c
	common	/tetra_zone/	ntetra,tetra,tetra_neighbour
	common	/tetra_stat/	tetra_status,tetra_orient,
     1				tetra_nindex
c
	save
c
	jtetra = tetra_neighbour(ivertex,itetra)
	jvertex = tetra_nindex(ivertex,itetra)
c
	if(jtetra.ne.0) then
		tetra_neighbour(jvertex,jtetra) = 0
	endif
c
	return
	end
c	Hpsort_key.f	Version 1 6/3/1995	Patrice Koehl
c
c	This subroutine rearranges an array in ascending order, and
c	provide an index of the ranked element
c
	subroutine hpsort_key(ra,index,n)
c
	integer	n,i,ir,j,l,idx
	integer	index(n)
c
	real*8	rra
	real*8	ra(n)
c
	save
c
	do 50 i = 1,n
		index(i) = i
50	continue
c
	if(n.lt.2) return
c
	l = n/2 + 1
	ir = n
100	continue
	if(l.gt.1) then
		l = l - 1
		rra = ra(l)
		idx = l
	else
		rra = ra(ir)
		idx = index(ir)
		ra(ir) = ra(1)
		index(ir) = index(1)
		ir = ir -1
		if(ir.eq.1) then
			ra(1) = rra
			index(1) = idx
			return
		endif
	endif
	i = l
	j = l+l
200	if(j.le.ir) then
		if(j.lt.ir) then
			if(ra(j).lt.ra(j+1)) j = j + 1
		endif
		if(rra.lt.ra(j)) then
			ra(i) = ra(j)
			index(i) = index(j)
			i = j
			j = j + j
		else
			j = ir + 1
		endif
		goto 200
	endif
	ra(i) = rra
	index(i) = idx
	goto 100
c
	end
C	Portable Random number generator.
C
C	From Numerical Recipes, W.H.Press, B.P.Flannery
C				S.A.Teukolsky, W.T.Vetterling
C		Cambridge Univ. Press
C
C	Function that returns a uniform random deviate between 0.0 and 1.0.
C	Set idum to any negative value to initialize the sequence.
C
	function ran2(idum)
C	-------------------
C
	integer	m,ia,ic
	real	rm,ran2
c
	parameter (m=714025,ia=1366,ic=150889,rm=1./m)
C
	integer	ir(97), iy,iff,idum,j
	data iff /0/
c
	save iy,ir
C
C	BEGIN.
C
	if(idum.lt.0.or.iff.eq.0) then
		iff=1
		idum=mod(ic-idum,m)
		do 11 j=1,97		! init the shuffle table
			idum=mod(ia*idum+ic,m)
			ir(j)=idum
11		continue
		idum=mod(ia*idum+ic,m)
		iy=idum
	endif
	j=1+(97*iy)/m
	if(j.gt.97.or.j.lt.1) then
		write(6,*) 'j=',j
		write(6,*) 'iy,m :',iy,m
		pause
	endif
	iy=ir(j)
C
C	RETURNED VALUE
C
	ran2=iy*rm
C
	idum=mod(ia*idum+ic,m)
	ir(j)=idum
C
	return
	end
c	Isort_indx.f		Version 1 1/10/2002	Patrice Koehl
c
c	This subroutine sorts 4 integers, gives the number of swaps
c	required for sorting, and provides the rank of the initial
c	numbers in the sorted array
c
	subroutine  isort_indx(list,idx,nswap,n)
c
	integer	i,j,a,n,nswap
	integer	list(n),idx(n)
c
	save
c
	do 100 i = 1,n
		idx(i) = i
100	continue
c
	nswap = 0
c
	do 300 i = 1,n-1
		do 200 j = i+1,n
			if(list(i).gt.list(j)) then
				a = list(i)
				list(i) = list(j)
				list(j) = a
				a = idx(i)
				idx(i) = idx(j)
				idx(j) = a
				nswap = nswap + 1
			endif
200		continue
300	continue
c
	return
	end
c
c	Isort_swap.f		Version 1 1/10/2002	Patrice Koehl
c
c	This subroutine sorts n integers, and gives the number of swaps
c	required for sorting
c
	subroutine  isort_swap(list,nswap,n)
c
	integer	i,j,a,n,nswap
c
	integer	list(n)
c
	save
c
	nswap = 0
c
	do 300 i = 1,n-1
		do 200 j = i+1,n
			if(list(i).gt.list(j)) then
				a = list(i)
				list(i) = list(j)
				list(j) = a
				nswap = nswap + 1
			endif
200		continue
300	continue
c
	return
	end
c
c	isort4_swap.f	Version 1 2/7/2002	Patrice Koehl
c
c	This subroutine sorts 4 numbers, and gives the number of swaps
c	required for sorting
c
	subroutine isort4_swap(a,b,c,d,nswap)
c
	integer	i,a,b,c,d,nswap
c
	save
c
	nswap = 0
c
	if(a.gt.b) then
		i = a
		a = b
		b = i
		nswap = 1
	endif
c
	if(a.gt.c) then
		i = a
		a = c
		c = i
		nswap = nswap + 1
	endif
c
	if(a.gt.d) then
		i = a
		a = d
		d = i
		nswap = nswap + 1
	endif
c
	if(b.gt.c) then
		i = b
		b = c
		c = i
		nswap = nswap + 1
	endif
c
	if(b.gt.d) then
		i = b
		b = d
		d = i
		nswap = nswap + 1
	endif
c
	if(c.gt.d) then
		i = c
		c = d
		d = i
		nswap = nswap + 1
	endif
c
	return
	end
c	Define_triangles.f	Version 1 3/11/2002	Patrice Koehl
c
c	This subroutine generates the list of all triangles
c	in the triangulation, using as input all tetrahedra
c
	subroutine define_triangles
c
	integer	ntetra_max,ntrig_max,nedge_max,npointmax
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter(MAX_FACET = 10000)


	parameter	(ntetra_max=MAX_TETRA)
	parameter	(ntrig_max=MAX_TRIG)
	parameter	(nedge_max=MAX_EDGE)
	parameter	(npointmax=MAX_POINT)
c
	integer	i,j,k,l
	integer	ntetra,ntrig
	integer	jtetra,jindex
	integer	nswap
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
	integer	vertice(4),neighbor(4),index(4)
c
	integer tetra(4,ntetra_max)
	integer tetra_neighbour(4,ntetra_max)
	integer	tetra_link(4,ntetra_max)
c
	integer	trig(3,ntrig_max),trig_link(3,ntrig_max)
	integer indx(4)
c
	integer*1 tetra_hull(ntetra_max),trig_hull(ntrig_max)
	integer*1 edge_hull(nedge_max),vertex_hull(npointmax)
c
	common	/tetra_zone/	ntetra,tetra,tetra_neighbour
	common	/tetra_stat/	tetra_status,tetra_orient,
     1				tetra_nindex
	common  /trig_zone/	ntrig,trig
	common  /links/		tetra_link,trig_link
	common  /hull/		tetra_hull,trig_hull,edge_hull,
     1				vertex_hull
c
	save
c
c	Loop over all tetrahedra, and generate all triangles
c
	ntrig = 0
c
	do 600 i = 1,ntetra
c
		if(tetra_status(i).ne.1) goto 600
c
		do 100 j = 1,4
			vertice(j) = tetra(j,i)
100		continue
c
		call isort_indx(vertice,indx,nswap,4)
c
		do 200 j = 1,4
			neighbor(j) = tetra_neighbour(indx(j),i)
			index(j) = tetra_nindex(indx(j),i)
200		continue
c
		if(mod(nswap,2).eq.1)tetra_orient(i)=-tetra_orient(i)
c
		do 300 j = 1,4
			tetra(j,i) = vertice(j)
			tetra_neighbour(j,i) = neighbor(j)
			tetra_nindex(j,i) = index(j)
			if(neighbor(j).ne.0) 
     1			tetra_nindex(index(j),neighbor(j)) = j
300		continue
c
		tetra_hull(i) = 0
c
c		Generate all triangles
c
		do 500 j = 1,4
c
			jtetra = tetra_neighbour(j,i)
c
			if(jtetra.eq.0.or.jtetra.gt.i) then
c
				ntrig = ntrig + 1
				l = 0
				do 400 k = 1,4
					if(k.eq.j) goto 400
					l = l + 1
					trig(l,ntrig) = tetra(k,i)
400				continue
c
				if(jtetra.eq.0) then
					tetra_hull(i) = 1
					trig_hull(ntrig) = 1
				else
					trig_hull(ntrig) = 0
				endif
c
				tetra_link(j,i) = ntrig
c
			else
c
				jindex = tetra_nindex(j,i)
				tetra_link(j,i) = tetra_link(jindex,
     1						jtetra)
c
			endif
c
500		continue
c
600	continue
c
	return
	end
c	Define_edges.f	Version 1 3/11/2002	Patrice Koehl
c
c	This subroutine generates the list of edges in the 
c	triangulation and link these edges to their "parent" triangles
c	It also labels each edge as 0 or 1 depending on 
c	the position of the edge with respect to the convex
c	hull
c
	subroutine define_edges
c
	integer	ntetra_max,ntrig_max,nedge_max,npointmax
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter(MAX_FACET = 10000)

	parameter	(ntetra_max=MAX_TETRA)
	parameter	(ntrig_max=MAX_TRIG)
	parameter	(nedge_max=MAX_EDGE)
	parameter	(npointmax=MAX_POINT)
c
	integer	i,idx
	integer	ipos,jpos
	integer	a,b,imask
	integer	trig1,trig2,triga,trigb,ipair
	integer	itrig1,itrig2
	integer	ktetra,trig_in,trig_out,npass,itetra
	integer	ntetra,ntrig,nedge
c
	integer	mask(6)
	integer	pair(2,6)
	integer	trig_info(2,6),trig_pos(2,6)
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
	integer*1 tetra_mask(ntetra_max)
c
	integer tetra(4,ntetra_max)
	integer tetra_neighbour(4,ntetra_max)
	integer	tetra_link(4,ntetra_max)
c
	integer	trig(3,ntrig_max)
	integer	trig_link(3,ntrig_max)
c
	integer	edge(2,nedge_max)
c
	integer*1 tetra_hull(ntetra_max),trig_hull(ntrig_max)
	integer*1 edge_hull(nedge_max),vertex_hull(npointmax)
c
	common	/tetra_zone/	ntetra,tetra,tetra_neighbour
	common	/tetra_stat/	tetra_status,tetra_orient,
     1				tetra_nindex
	common  /trig_zone/	ntrig,trig
	common  /edge_zone/	nedge,edge
	common  /links/		tetra_link,trig_link
	common  /hull/		tetra_hull,trig_hull,edge_hull,
     1				vertex_hull
c
	save
c
c	data pair/1,2,1,3,1,4,2,3,2,4,3,4/
c	data trig_info/3,4,2,4,2,3,1,4,1,3,1,2/
c	data trig_pos/3,3,3,2,2,2,3,1,2,1,1,1/
	data pair/3,4,2,4,2,3,1,4,1,3,1,2/
	data trig_info/1,2,1,3,1,4,2,3,2,4,3,4/
	data trig_pos/1,1,2,1,3,1,2,2,3,2,3,3/
	data mask/2,4,8,16,32,64/
c
	nedge = 0
c
	do 100 idx=1,ntetra
		tetra_mask(idx) = 0
100	continue
c
	do 600 idx = 1,ntetra
c
	   if(tetra_status(idx).ne.1) goto 600
c
	   imask = tetra_mask(idx)
c
	   do 500 i = 1,6
c
		if(iand(imask,mask(i)).ne.0) goto 500
c
		a = tetra(pair(1,i),idx)
		b = tetra(pair(2,i),idx)
c
		nedge = nedge + 1
		edge(1,nedge) = a
		edge(2,nedge) = b
c
		trig1 = trig_info(1,i)
		trig2 = trig_info(2,i)
c
		itrig1 = tetra_link(trig1,idx)
		itrig2 = tetra_link(trig2,idx)
c
		ipos = trig_pos(1,i)
		jpos = trig_pos(2,i)
c
		trig_link(ipos,itrig1) = nedge
		trig_link(jpos,itrig2) = nedge
c
		ktetra = idx
		npass = 1
		trig_out = trig1
c
c		Loop over all the tetrahedra in the link of
c		the edge: start in one direction, and if we 
c		hit the convex hull before cycling back to
c		the initial tetrahedron, restart in the other
c		direction
c
200		continue
c
		itetra = tetra_neighbour(trig_out,ktetra)
c
c		Leave this side of the link if we hit
c		the convex hull
c		
		if(itetra.eq.0) goto 300
c
c		Leave the loop completely if we have described
c		the full cycle
c
		if(itetra.eq.idx) goto 400
c
c		Identify the position of edge (ab) in
c		tetrahedron itetra
c
		if(a.eq.tetra(1,itetra)) then
			if(b.eq.tetra(2,itetra)) then
				ipair = 6
			elseif(b.eq.tetra(3,itetra)) then
				ipair=5
			else
				ipair=4
			endif
		elseif(a.eq.tetra(2,itetra)) then
			if(b.eq.tetra(3,itetra)) then
				ipair=3
			else
				ipair=2
			endif
		else
			ipair = 1
		endif
c
		tetra_mask(itetra)=tetra_mask(itetra)+mask(ipair)
c
		trig_in = tetra_nindex(trig_out,ktetra)
c
		triga = trig_info(1,ipair)
		trigb = trig_info(2,ipair)
c
		trig_out = triga
		ipos = trig_pos(1,ipair)
		if(trig_in.eq.triga) then
			trig_out = trigb
			ipos = trig_pos(2,ipair)
		endif
c
		trig_link(ipos,tetra_link(trig_out,itetra)) = nedge
c
		ktetra = itetra
c
		goto 200
c
300		continue
c
		if(npass.eq.2) goto 400
		npass = npass + 1
		ktetra = idx
		trig_out = trig2
c
		goto 200
c
400		continue
c
500	    continue
c
600	continue
c
c
	do 700 i = 1,ntrig
c
		if(trig_hull(i).eq.1) then
			edge_hull(trig_link(1,i)) = 1
			edge_hull(trig_link(2,i)) = 1
			edge_hull(trig_link(3,i)) = 1
			vertex_hull(trig(1,i)) = 1
			vertex_hull(trig(2,i)) = 1
			vertex_hull(trig(3,i)) = 1
		endif
c
700	continue
c
	return
	end
c	Peel.f	Version 1 4/1/2002	Patrice Koehl
c
c	This subroutine removes the flat tetrahedra at the boundary 
c	of the DT
c
	subroutine peel
c
	integer	ntetra_max,npointmax
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter(MAX_FACET = 10000)

	parameter (ntetra_max = MAX_TETRA)
	parameter (npointmax = MAX_POINT)
c
	integer	i,j,k,l
	integer	ia,ib,ic,id,val
	integer	ntetra
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
	integer tetra(4,ntetra_max)
	integer tetra_neighbour(4,ntetra_max)
c
	real*8  eps,scale
	real*8	vol
c
	real*8	coord(3*npointmax),radius(npointmax)
	real*8	coord4(npointmax)
c
	common  /xyz_vertex/	coord,radius,coord4
	common	/tetra_zone/	ntetra,tetra,tetra_neighbour
	common	/tetra_stat/	tetra_status,tetra_orient,
     1				tetra_nindex
	common  /gmp_info/	scale,eps
c
	save
c
c	Loop over all tetrahedra, and test the tetrahedra at
c	the boundary
c
	do 400 i = 1,ntetra
c
		if(tetra_status(i).ne.1) goto 400
c
		do 200 j = 1,4
			if(tetra_neighbour(j,i).eq.0) goto 300
200		continue
c
c		If we get here, the tetrahedron idx is interior, and
c		cannot be flat (see Edelsbrunner,...)
c
		goto 400
c
300		continue
c
c		This is a tetrahedron at the boundary: we test
c		if it is flat, i.e. if its volume is 0
c
		ia = tetra(1,i)
		ib = tetra(2,i)
		ic = tetra(3,i)
		id = tetra(4,i)
c
		call tetra_vol(coord,ia,ib,ic,id,vol,npointmax)
c
		if(abs(vol).lt.eps) then
			call minor4_gmp(ia,ib,ic,id,val)
			if(val.eq.0) tetra_status(i)=-1
		endif
c
400	continue
c
c	Now we remove that flat tetrahedra, and update the links
c	to their neighbours
c
	do 600 i = 1,ntetra
c
		if(tetra_status(i).eq.-1) then
			tetra_status(i) = 0
			do 500 j = 1,4
				k = tetra_neighbour(j,i)
				if(k.ne.0) then
					l = tetra_nindex(j,i)
					tetra_neighbour(l,k)=0
				endif
500			continue
c
		endif
c
600	continue
c
	return
	end
c
c	tetra_vol.f	Version 1 4/1/2002	Patrice Koehl
c
c	This subroutine computes the volume of a tetrahedron
c
	subroutine tetra_vol(coord,ia,ib,ic,id,vol,npointmax)
c
c	Input:
c		- coord:	array containing coordinates of all vertices
c		- ia,ib,ic,id	four vertices defining the tetrahedron
c		- npointmax	maximum number of vertices
c
c	Output:
c		- vol		volume of the tetrahedron (floating
c				point calculation)
c
	integer	i
	integer	ia,ib,ic,id
	integer	npointmax
c
	save
c
	real*8	vol
	real*8	ad(3),bd(3),cd(3)
	real*8	Sbcd(3)
	real*8	coord(3*npointmax)
c
c	The volume of the tetrahedron is proportional to:
c
c	vol = det | a(1)  a(2)  a(3)  1|
c		  | b(1)  b(2)  b(3)  1|
c		  | c(1)  c(2)  c(3)  1|
c		  | d(1)  d(2)  d(3)  1|
c
c	After substracting the last row from the first 3 rows, and
c	developping with respect to the last column, we obtain:
c
c	vol = det | ad(1)  ad(2)  ad(3) |
c		  | bd(1)  bd(2)  bd(3) |
c		  | cd(1)  cd(2)  cd(3) |
c
c	where ad(i) = a(i) - d(i), ...
c
	do 100 i = 1,3
		ad(i) = coord(3*(ia-1)+i) - coord(3*(id-1)+i)
		bd(i) = coord(3*(ib-1)+i) - coord(3*(id-1)+i)
		cd(i) = coord(3*(ic-1)+i) - coord(3*(id-1)+i)
100	continue
c
	Sbcd(3) = bd(1)*cd(2) - cd(1)*bd(2)
	Sbcd(2) = bd(1)*cd(3) - cd(1)*bd(3)
	Sbcd(1) = bd(2)*cd(3) - cd(2)*bd(3)
c
	vol = ad(1)*Sbcd(1) - ad(2)*Sbcd(2) + ad(3)*Sbcd(3)
c
	return
	end
c
c	tetra_order.f		Version 1 11/24/2000	Patrice Koehl
c
c	This subroutine organises the tetrahedron in a "depth" manner
c
	subroutine tetra_order
c
	integer	ntetra_max
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter(MAX_FACET = 10000)


	parameter (ntetra_max=MAX_TETRA)
c
	integer i,j,k,first
	integer	ntetra,norder,norder2
	integer idx
	integer	nvisited,pointer,pointer2
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
c	Information on the tetrahedra of the regular
c	triangulation
c
	integer tetra(4,ntetra_max)
	integer tetra_neighbour(4,ntetra_max)
c
	integer*1 tetra_visited(ntetra_max)
c
	integer	order(ntetra_max)
c
	common  /tetra_zone/    ntetra,tetra,tetra_neighbour
	common  /tetra_stat/    tetra_status,tetra_orient,
     1                          tetra_nindex
	common /tetra_print/	order
c
	save
c
c	Loop over all tetrahedra:
c
	nvisited = 0
	norder   = 0
	pointer  = 0
c
	first = 0
	do 25 i = 1,ntetra
		if(tetra_status(i).ne.1) then
			tetra_visited(i) = 1
			nvisited = nvisited + 1
		else
			if(first.eq.0) then
				k = i
				first = 1
			endif
		endif
25	continue
c
	if(k.eq.0) then
		return
	else
		norder = norder + 1
		order(norder) = k
		pointer = 1
		tetra_visited(k) = 1
		nvisited = nvisited + 1
	endif
c
50	continue
c
	norder2 = norder
c
	do 200 i = pointer, norder
c
		idx = order(i)
c
		do 100 j = 1,4
c
			k = tetra_neighbour(j,idx)
c
			if(k.eq.0)  goto 100
c
			if(tetra_visited(k).eq.1) goto 100
			norder2 = norder2 + 1
			order(norder2) = k
			tetra_visited(k) = 1
			nvisited = nvisited + 1
c
100		continue
c
200	continue
c
	pointer = norder + 1
	norder = norder2
c
	if(nvisited.eq.ntetra) goto 300
c
	goto 50
c
300	continue
c
	return
	end
c	sort4_sign.f	Version 1 9/30/2006	Patrice Koehl
c
c	This subroutine sorts the list of 4 numbers, and computes the
c	signature of the permutation
c
	subroutine sort4_sign(list,index,nswap,n)
c
	integer	i,j,k,nswap,n
	integer	list(n),index(n)
c
	save
c
	do 100 i = 1,n
		index(i) = i
100	continue
c
	nswap = 1
	do 300 i = 1,n-1
		do 200 j = i+1,n
			if(list(i).gt.list(j)) then
				k = list(i)
				list(i) = list(j)
				list(j) = k
				k = index(i)
				index(i) = index(j)
				index(j) = k
				nswap = -nswap
			endif
200		continue
300	continue
c
	return
	end
c
c	reorder_tetra.f	Version 1 9/30/2006	Patrice Koehl
c
c	This subroutine reorders the vertices of a list of tetrahedron,
c	such that now the indices are in increasing order
c
	subroutine reorder_tetra(new,list_tetra)
c
	integer	ntetra_max
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter(MAX_FACET = 10000)


	parameter (ntetra_max = MAX_TETRA)
c
	integer	i,j,idx,nswap,new
	integer	ntetra
	integer	index(4)
	integer	vertex(4),neighbour(4),nindex(4)
c
	integer		tetra(4,ntetra_max)
	integer		tetra_neighbour(4,ntetra_max)
c
	integer*1	tetra_status(ntetra_max)
	integer*1	tetra_orient(ntetra_max)
	integer*1	tetra_nindex(4,ntetra_max)
	integer		list_tetra(new)
c
	common  /tetra_zone/	ntetra,tetra,tetra_neighbour
     	common /tetra_stat/	tetra_status,tetra_orient,
     1				tetra_nindex
c
	save
c
	do 400 idx = 1,new
c
		i = list_tetra(idx)
		if(tetra_status(i).eq.1) then
c
			do 100 j = 1,4
				vertex(j) = tetra(j,i)
100			continue
c
			call sort4_sign(vertex,index,nswap,4)
c
			do 200 j = 1,4
				neighbour(j)=tetra_neighbour(index(j),i)
				nindex(j)=tetra_nindex(index(j),i)
				if(neighbour(j).ne.0) then
				tetra_nindex(nindex(j),neighbour(j))=j
				endif
200			continue
c
			do 300 j = 1,4
				tetra(j,i)=vertex(j)
				tetra_neighbour(j,i)=neighbour(j)
				tetra_nindex(j,i) = nindex(j)
300			continue
c
			tetra_orient(i) = tetra_orient(i)*nswap
c
		endif
c
400	continue
c
	return
	end
c	vertex_connect2.f	Version 1 3/11/2002	Patrice Koehl
c
c	This subroutine defines all tetrahedrons attached to a
c	vertex
c
	subroutine vertex_connect2(idx)
c
	integer	ntetra_max,nlinkmax
c
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter(MAX_FREE = 10000)
	parameter(MAX_FACET = 10000)
	parameter(MAX_LINK = 10000)

	parameter	(ntetra_max = MAX_TETRA)
	parameter	(nlinkmax = MAX_LINK)
c
	integer	i,j,k,idx
	integer	ntetra
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
	integer tetra(4,ntetra_max)
	integer tetra_neighbour(4,ntetra_max)
	integer	tetra_link(4,ntetra_max)
	integer	vertex_nlink, vertex_link(nlinkmax)
c
	common	/tetra_zone/	ntetra,tetra,tetra_neighbour
	common	/tetra_stat/	tetra_status,tetra_orient,
     1				tetra_nindex
	common  /vertex_tet/    vertex_nlink,vertex_link
c
	save
c
	vertex_nlink = 0
c
	do 600 i = 1,ntetra
c
		if(tetra_status(i).ne.1) goto 600
c
		do 300 j = 1,4
			k = tetra(j,i)
			if(k.eq.idx) then
				vertex_nlink = vertex_nlink + 1
				vertex_link(vertex_nlink) = i			
			endif
300		continue
c
600	continue
c
	return
	end
c	Hpsort_three.f	Version 1 6/3/1995	Patrice Koehl
c
c	This subroutine rearranges an array in ascending order, and
c	provide an index of the ranked element
c
	subroutine hpsort_three(ra,index,n)
c
	integer	n,i,ir,j,l,idx,comp3,k
	integer	index(n)
c
	save
c
	real*8	rra(3)
	real*8	ra(3,n)
c
	do 50 i = 1,n
		index(i) = i
50	continue
c
	if(n.lt.2) return
c
	l = n/2 + 1
	ir = n
100	continue
	if(l.gt.1) then
		l = l - 1
		do 10 k = 1,3
			rra(k) = ra(k,l)
10		continue
		idx = l
	else
		do 20 k = 1,3
			rra(k) = ra(k,ir)
20		continue
		idx = index(ir)
		do 30 k = 1,3
			ra(k,ir) = ra(k,1)
30		continue
		index(ir) = index(1)
		ir = ir -1
		if(ir.eq.1) then
			do 40 k = 1,3
				ra(k,1) = rra(k)
40			continue
			index(1) = idx
			return
		endif
	endif
	i = l
	j = l+l
200	if(j.le.ir) then
		if(j.lt.ir) then
			if(comp3(ra(1,j),ra(1,j+1)).eq.1) j = j + 1
		endif
		if(comp3(rra,ra(1,j)).eq.1) then
			do 210 k = 1,3
				ra(k,i) = ra(k,j)
210			continue
			index(i) = index(j)
			i = j
			j = j + j
		else
			j = ir + 1
		endif
		goto 200
	endif
	do 220 k = 1,3
		ra(k,i) = rra(k)
220	continue
	index(i) = idx
	goto 100
c
	end
c
c	Comp3.f		version 1 10/17/2006	Patrice Koehl
c
c	This function compares 2 arrays of 3 real numbers
c
	function comp3(a,b)
c
	integer	i,comp3,k
	real*8	a(3),b(3)
c
	save
c
	comp3 = 0
c
	do 100 k = 1,3
		if(a(k).lt.b(k)) then
			comp3 = 1
			return
		elseif(a(k).gt.b(k)) then
			return
		endif
100	continue
c
	return
	end
