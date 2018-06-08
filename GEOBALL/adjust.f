c	Adjust.f		Version 1 6/17/2005	Patrice Koehl
c
c	This file contains a set of routines that adjust the number of
c	spheres included in the calculation to be always at least 4,
c	such that the regular triangulation can be computed, and the
c	alpha shape derived
c
c	Copyright (C) 2005 Patrice Koehl
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
c	adjust_nsphere.f		Version 1 6/17/2005	Patrice Koehl
c
c	This subroutine adds spheres, if needed
c
	subroutine adjust_nsphere(coord_sph,rad,nsphere)
c
	integer	npointmax
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)

	parameter (npointmax = MAX_POINT)
c
	integer i,j,k
	integer	nsphere
	integer	sign(3,3)
c
c	Information on the vertices
c
	real*8	Rmax
	real*8	Dmax(3)
	real*8	rad(npointmax)
	real*8	coord_sph(3*npointmax)
c
	data ((sign(i,j),j=1,3),i=1,3) /1,1,1,1,-1,-1,-1,1,1/
c
	save
c
c
c	Do nothing if we already have at least 4 balls
c
	if(nsphere.ge.4) return
c
c	Get bounding box of the current set of spheres, as well as max
c	radius
c
	do 100 i = 1,3
		Dmax(i) = coord_sph(i)
100	continue
	Rmax = rad(1)
c
	do 300 i = 2,nsphere
		do 200 j = 1,3
			if(Dmax(j).lt.coord_sph(3*(i-1)+j)) then
				Dmax(j) = coord_sph(3*(i-1)+j)
			endif
200		continue
		if(Rmax.lt.rad(i)) Rmax = rad(i)
300	continue
c
c	Now add point(s) with center at Dmax + 3*Rmax, and radius Rmax/20
c
	do 500 i = nsphere+1,4
		j = i - nsphere
		do 400 k = 1,3
			coord_sph(3*(i-1)+k)=sign(i,k)*(Dmax(k)+2*Rmax)
400		continue
		rad(i) = Rmax/20
500	continue
c
	nsphere = 4
c
	return
	end
c	readjust_nsphere.f		Version 1 6/17/2005	Patrice Koehl
c
c	This subroutine removes the artificial spheres, if needed
c
	subroutine readjust_nsphere(nsphere,nred,listred)
c
	integer	npointmax
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)

	parameter (npointmax = MAX_POINT)
c
	integer i,j
	integer	npoints,nvertex,nsphere,nred
c
	integer	listred(nred)
c
c	Information on the vertices
c
	integer redinfo(npointmax)
c
	common  /vertex_zone/   npoints,nvertex,redinfo
c
	save
c
c	Do nothing if we already have at least 4 balls
c
	if(nsphere.ge.4) return
c
	do 100 i = nsphere+5,8
		redinfo(i) = 1
100	continue
c
	npoints = nsphere
	nvertex = npoints + 4
c
	j = 0
	do 200 i = 1,nred
		if(listred(i).le.nsphere) then
			j = j + 1
			listred(j) = listred(i)
		endif
200	continue
c
	return
	end
