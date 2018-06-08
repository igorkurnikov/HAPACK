c	prepare_deriv.f		Version 1 : 10/28/2002
c
c	Copyright (C) 2002 Patrice Koehl
c
c	This file contains a series of subroutines used to compute the
c	volume area of a union of balls, and optionally the derivatives
c	of the volume with respect to the coordinates of the centers of the ball.
c	As a side result, it also provides the surface area (and its derivatives)
c	of the union of balls.
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
c
c	prepare_deriv.f	
c
c	Copyright (C) 2002 Patrice Koehl
c
c	This program tests the "direction" of a tetrahedron (ccw or not),
c	and build the link lst of all triangle
c
	subroutine prepare_deriv
c
	integer ntetra_max,ntrig_max
c
	parameter (ntetra_max= 200000)
	parameter (ntrig_max = 400000)
c
	integer	i,j,k,l,idx
	integer	trig1,trig2,trig3,trig4
	integer	ntrig,ntetra
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
	integer tetra(4,ntetra_max)
	integer	trig(3,ntrig_max)
c
	integer	tetra_neighbour(4,ntetra_max)
c
	integer tetra_link(4,ntetra_max)
	integer trig_link(3,ntrig_max)
c
	integer*1	nlink_trig(ntrig_max)
	integer		link_trig(2,ntrig_max)
c
	common  /tetra_zone/	ntetra,tetra,tetra_neighbour
	common  /trig_zone/	ntrig,trig
	common  /tetra_stat/	tetra_status,tetra_orient,
     1				tetra_nindex
	common  /links/		tetra_link,trig_link
	common /trig_link/	nlink_trig,link_trig
c
	save
c
	do 100 idx = 1,ntrig
		nlink_trig(idx) = 0
100	continue
c
	do 300 idx = 1,ntetra
c
		if(tetra_status(idx).ne.1) goto 300
c
		i = tetra(1,idx)
		j = tetra(2,idx)
		k = tetra(3,idx)
		l = tetra(4,idx)
c
		trig1 = tetra_link(1,idx)
		trig2 = tetra_link(2,idx)
		trig3 = tetra_link(3,idx)
		trig4 = tetra_link(4,idx)
c
		nlink_trig(trig1) = nlink_trig(trig1) + 1
		link_trig(nlink_trig(trig1),trig1) = i
		nlink_trig(trig2) = nlink_trig(trig2) + 1
		link_trig(nlink_trig(trig2),trig2) = j
		nlink_trig(trig3) = nlink_trig(trig3) + 1
		link_trig(nlink_trig(trig3),trig3) = k
		nlink_trig(trig4) = nlink_trig(trig4) + 1
		link_trig(nlink_trig(trig4),trig4) = l
c
300	continue
c
	return
	end
