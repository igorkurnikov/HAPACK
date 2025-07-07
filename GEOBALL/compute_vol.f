c	compute_vol.f
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
c	Based on the alpha shape, this subroutine computes the
c	volume and surface of each atom of the protein
c	using the Inclusion-Exclusion formula
c
c	include "defines.h"
c
	subroutine compute_vol(coefasp,option)
c
	integer npointmax,nedge_max,ntrig_max,ntetra_max,natot,ncortot
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)
	parameter       (npointmax = MAX_POINT)
	parameter       (natot     = MAX_POINT)
	parameter       (ntetra_max = MAX_TETRA)
	parameter       (ntrig_max = MAX_TRIG)
	parameter       (nedge_max = MAX_EDGE)
	parameter       (ncortot   = 3*npointmax)
c
	integer isave
	integer	i,j,k,l,m
	integer	idx,pair1,pair2,pair3,pair4,pair5,pair6
	integer	trig1,trig2,trig3,trig4
	integer	flaga,flagb,flagc
	integer ntetra,ntrig,nedge
	integer npoints,nvertex
	integer	nlink,l1,l2
	integer	option,ivol
c
c	Information on tetrahedra
c
	integer*1 tetra_status(ntetra_max),tetra_orient(ntetra_max)
	integer*1 tetra_nindex(4,ntetra_max)
c
	integer tetra(4,ntetra_max)
	integer tetra_neighbour(4,ntetra_max)
	integer tetra_link(4,ntetra_max)
c
c	Information on triangles
c
	integer trig(3,ntrig_max)
	integer trig_link(3,ntrig_max)
	integer	link_trig(2,ntrig_max)
	integer*1 trig_coef(ntrig_max)
	integer*1 nlink_trig(ntrig_max)
c
c	Information on edges
c
	integer edge(2,nedge_max)
	integer*1 edge_coef(nedge_max)
c
c	Information on vertices
c
	integer redinfo(npointmax)
c
c
c       Information for alpha complex
c
	integer*1       tetra_alp(ntetra_max)
	integer*1       trig_alp(ntrig_max)
	integer*1       edge_alp(nedge_max)
	integer*1       vertex_alp(npointmax)
c
	real*8	pi,twopi
	real*8	voltot,surftot,wsurf,wvol
	real*8	distij,distij2
	real*8	eps,eps1,eps2,eps3,eps4
	real*8	ra,rb,rc,rd,ra2,rb2,rc2,rd2
	real*8	wa,wb,wc,wd,wl1,wl2
	real*8	surfa,surfb,surfc,surfd,vola,volb,volc,vold
	real*8	dist1,dist2,dist3,dist4,dist5,dist6
	real*8	d2_1,d2_2,d2_3,d2_4,d2_5,d2_6
	real*8	coefa,coefb,coefc,coefd
	real*8  sh_abc,sh_acb,sh_bca,sh_abd,sh_adb,sh_bda
	real*8  sh_acd,sh_adc,sh_cda,sh_bcd,sh_bdc,sh_cdb
	real*8	vol(natot),surf(natot)
	real*8	a(3),b(3),c(3),d(3)
	real*8  pabc(3),pacb(3),pabd(3),padb(3),padc(3),pacd(3)
	real*8  pbcd(3),pbdc(3)
	real*8	x(3),y(3)
	real*8	angles(3)
	real*8	dsurfa2(3,2),dsurfb2(3,2)
	real*8	dsurfa3(3,3),dsurfb3(3,3),dsurfc3(3,3)
	real*8	dsurfa4(3,4),dsurfb4(3,4),dsurfc4(3,4),dsurfd4(3,4)
	real*8	dvola2(3,2),dvolb2(3,2)
	real*8	dvola3(3,3),dvolb3(3,3),dvolc3(3,3)
	real*8	dvolda3(3,3),dvoldb3(3,3),dvoldc3(3,3)
	real*8	dvola4(3,4),dvolb4(3,4),dvolc4(3,4),dvold4(3,4)
	real*8	dsurf(3,natot),dvol(3,natot)
	real*8	trig_ang(3,ntrig_max)
	real*8	ang_abc(3),ang_abd(3),ang_acd(3),ang_bcd(3)
	real*8  trig_eps(ntrig_max),trig_sh(3,ntrig_max)
	real*8  trig_dual1(3,ntrig_max),trig_dual2(3,ntrig_max)
	real*8	coefasp(npointmax)
c
	real*8	distpair(nedge_max),distpair2(nedge_max)
	real*8  radius(npointmax),radius2(npointmax)
	real*8	weight(npointmax)
	real*8  coord(3*npointmax)

c
	logical test
c
	common	/xyz_vertex/	coord,radius,weight
	common	/tetra_zone/	ntetra,tetra,tetra_neighbour
	common	/tetra_stat/	tetra_status,tetra_orient,
     1				tetra_nindex
	common	/trig_zone/	ntrig,trig
	common	/trig_fact/	trig_coef
	common	/edge_zone/	nedge,edge
	common	/vertex_zone/	npoints,nvertex,redinfo
	common	/links/		tetra_link,trig_link
	common /trig_link/	nlink_trig,link_trig
c	
	common /constants/      pi,twopi
	common /results/	wsurf,wvol,surftot,voltot,surf,vol,
     1				dsurf,dvol
	common  /alpha_complex/ tetra_alp,trig_alp,edge_alp,vertex_alp
c
	data isave /0/
c
	save
c
	if(isave.eq.0) then
		pi = acos(-1.d0)
		twopi = 2.d0 * pi
		isave = 1
	endif
c
	do 100 i = 1,nvertex
		do 50 j = 1,3
			dsurf(j,i) = 0
			dvol(j,i) = 0
50		continue
100	continue
c
	do 150 i = 1,nvertex
		radius2(i) = radius(i)*radius(i)
		surf(i) = 0
		vol(i)  = 0
150	continue
c
	call prepare_deriv
c
	do 200 i = 1,nedge
		edge_coef(i) = edge_alp(i)
200	continue
c
	do 250 i = 1,ntrig
		if(trig_alp(i).eq.1) then
			pair1 = trig_link(1,i)
			pair2 = trig_link(2,i)
			pair3 = trig_link(3,i)
			edge_coef(pair1) = 0
			edge_coef(pair2) = 0
			edge_coef(pair3) = 0
		endif
250	continue
c
	do 300 i = 1,ntrig
		if(trig_coef(i).ne.0) then
			pair1 = trig_link(1,i)
			pair2 = trig_link(2,i)
			pair3 = trig_link(3,i)
			edge_coef(pair1) = 1
			edge_coef(pair2) = 1
			edge_coef(pair3) = 1
		endif
300	continue
c
	do 350 idx = 1,nvertex
c
		if(vertex_alp(idx).eq.0) goto 350
		surf(idx) = 4*pi*radius2(idx)
		vol(idx)  = surf(idx)*radius(idx)/3
c		write(6,*) 'i,surf(i) :',idx,surf(idx)
c
350	continue
c
	do 550 idx = 1,nedge
c
c		write(6,*) 'edge :,i,j, status :',edge(1,idx),
c     1			edge(2,idx),edge_alp(idx)
c
		if(edge_alp(idx).eq.0) goto 550
c
		i = edge(1,idx)
		j = edge(2,idx)
c
		coefa = coefasp(i)
		coefb = coefasp(j)
c
		call distance2(coord,i,j,distij2,ncortot)
		distij = sqrt(distij2)
c
		distpair(idx)  = distij
		distpair2(idx) = distij2
c
		do 400 k = 1,3
			a(k) = coord(3*(i-1)+k)
			b(k) = coord(3*(j-1)+k)
400		continue
c
		ra = radius(i)
		rb = radius(j)
		ra2 = radius2(i)
		rb2 = radius2(j)
c
		call twosphere_vol(a,b,ra,ra2,rb,rb2,distij,distij2,
     1		surfa,surfb,vola,volb,dsurfa2,dsurfb2,dvola2,dvolb2,
     2		option)
c
		surf(i) = surf(i) - surfa
		surf(j) = surf(j) - surfb
c
		vol(i)  = vol(i) - vola
		vol(j)  = vol(j) - volb
c
		if(option.eq.0) goto 550
c
		do 450 k = 1,3
			dsurf(k,i) = dsurf(k,i) - coefa*dsurfa2(k,1)
     1			-coefb*dsurfb2(k,1)
			dsurf(k,j) = dsurf(k,j) - coefa*dsurfa2(k,2)
     1			-coefb*dsurfb2(k,2)
450		continue
c
		do 500 k = 1,3
			dvol(k,i) = dvol(k,i) - coefa*dvola2(k,1)
     1			-coefb*dvolb2(k,1)
			dvol(k,j) = dvol(k,j) - coefa*dvola2(k,2)
     1			-coefb*dvolb2(k,2)
500		continue
c
550	continue
c
	do 1000 idx = 1,ntrig
c
		if(trig_alp(idx).eq.0) goto 1000
c
		i = trig(1,idx)
		j = trig(2,idx)
		k = trig(3,idx)
c
		do 600 m = 1,3
			a(m) = coord(3*(i-1)+m)
			b(m) = coord(3*(j-1)+m)
			c(m) = coord(3*(k-1)+m)
600		continue
c
		ra = radius(i)
		rb = radius(j)
		rc = radius(k)
c
		ra2 = radius2(i)
		rb2 = radius2(j)
		rc2 = radius2(k)
c
		wa = 0.5d0*weight(i)
		wb = 0.5d0*weight(j)
		wc = 0.5d0*weight(k)
c
		coefa = coefasp(i)
		coefb = coefasp(j)
		coefc = coefasp(k)
c
		pair1 = trig_link(3,idx)
		pair2 = trig_link(2,idx)
		pair3 = trig_link(1,idx)
c
		dist1 = distpair(pair1)
		dist2 = distpair(pair2)
		dist3 = distpair(pair3)
c
		d2_1 = distpair2(pair1)
		d2_2 = distpair2(pair2)
		d2_3 = distpair2(pair3)
c
		ivol = 0
		if(trig_coef(idx).eq.0) goto 775
c
		ivol = 1
		call threevol_dist(a,b,c,ra,rb,rc,ra2,rb2,rc2,
     1		wa,wb,wc,dist1,dist2,dist3,d2_1,d2_2,d2_3,
     2		surfa,surfb,surfc,vola,volb,volc,
     3		dsurfa3,dsurfb3,dsurfc3,dvola3,dvolb3,dvolc3,
     4		dvolda3,dvoldb3,dvoldc3,eps,pabc,pacb,angles,
     5		sh_abc,sh_acb,sh_bca,option)
c
		surf(i) = surf(i) + 0.5d0*surfa*trig_coef(idx)
		surf(j) = surf(j) + 0.5d0*surfb*trig_coef(idx)
		surf(k) = surf(k) + 0.5d0*surfc*trig_coef(idx)
c
		vol(i)  = vol(i) + 0.5d0*vola*trig_coef(idx)
		vol(j)  = vol(j) + 0.5d0*volb*trig_coef(idx)
		vol(k)  = vol(k) + 0.5d0*volc*trig_coef(idx)
c
		trig_ang(1,idx) = angles(1)
		trig_ang(2,idx) = angles(2)
		trig_ang(3,idx) = angles(3)
c
		trig_eps(idx) = eps
		trig_sh(1,idx) = sh_abc
		trig_sh(2,idx) = sh_acb
		trig_sh(3,idx) = sh_bca
c
		do 850 l = 1,3
			trig_dual1(l,idx) = pabc(l)
			trig_dual2(l,idx) = pacb(l)
850		continue
c
		if(option.eq.0) goto 1000
c
		do 900 l = 1,3
			dsurf(l,i) = dsurf(l,i) + 0.5d0*trig_coef(idx)
     1			*(coefa*dsurfa3(l,1)+ coefb*dsurfb3(l,1)
     2			+coefc*dsurfc3(l,1))
			dsurf(l,j) = dsurf(l,j) + 0.5d0*trig_coef(idx)
     1			*(coefa*dsurfa3(l,2)+ coefb*dsurfb3(l,2)
     2			+coefc*dsurfc3(l,2))
			dsurf(l,k) = dsurf(l,k) + 0.5d0*trig_coef(idx)
     1			*(coefa*dsurfa3(l,3)+ coefb*dsurfb3(l,3)+
     2			coefc*dsurfc3(l,3))
900		continue
c
		do 950 l = 1,3
			dvol(l,i) = dvol(l,i) + 0.5d0*trig_coef(idx)*
     1			(coefa*dvola3(l,1)+ coefb*dvolb3(l,1)
     2			+coefc*dvolc3(l,1))
			dvol(l,j) = dvol(l,j) + 0.5d0*trig_coef(idx)*
     1			(coefa*dvola3(l,2)+ coefb*dvolb3(l,2)
     2			+coefc*dvolc3(l,2))
			dvol(l,k) = dvol(l,k) +  0.5d0*trig_coef(idx)*
     1			(coefa*dvola3(l,3)+ coefb*dvolb3(l,3)
     2			+coefc*dvolc3(l,3))
			dvol(l,i) = dvol(l,i) + 
     1			(coefa*dvolda3(l,1)+ coefb*dvoldb3(l,1)
     2			+coefc*dvoldc3(l,1))
			dvol(l,j) = dvol(l,j) + 
     1			(coefa*dvolda3(l,2)+ coefb*dvoldb3(l,2)
     2			+coefc*dvoldc3(l,2))
			dvol(l,k) = dvol(l,k) + 
     1			(coefa*dvolda3(l,3)+ coefb*dvoldb3(l,3)
     2			+coefc*dvoldc3(l,3))
950		continue
c
775		continue
c
		if(ivol.eq.0) then
c
			call threevol_dir(a,b,c,ra,rb,rc,ra2,rb2,rc2,
     1			wa,wb,wc,dist1,dist2,dist3,d2_1,d2_2,d2_3,
     2			dvola3,dvolb3,dvolc3,angles,pabc,pacb,eps,
     3			sh_abc,sh_acb,sh_bca,option)
c
			trig_ang(1,idx) = angles(1)
			trig_ang(2,idx) = angles(2)
			trig_ang(3,idx) = angles(3)
c
			trig_eps(idx) = eps
			trig_sh(1,idx) = sh_abc
			trig_sh(2,idx) = sh_acb
			trig_sh(3,idx) = sh_bca
c
			do 650 l = 1,3
				trig_dual1(l,idx) = pabc(l)
				trig_dual2(l,idx) = pacb(l)
650			continue
c
			if(option.eq.0) goto 1000
c
			do 700 l = 1,3
				dvol(l,i) = dvol(l,i) + 
     1				(coefa*dvola3(l,1)+ coefb*dvolb3(l,1)
     2				+coefc*dvolc3(l,1))
				dvol(l,j) = dvol(l,j) + 
     1				(coefa*dvola3(l,2)+ coefb*dvolb3(l,2)
     2				+coefc*dvolc3(l,2))
				dvol(l,k) = dvol(l,k) + 
     1				(coefa*dvola3(l,3)+ coefb*dvolb3(l,3)
     2				+coefc*dvolc3(l,3))
700			continue
c
		endif
c
		if(edge_coef(pair1).eq.0.and.edge_coef(pair2)
     1		.eq.0.and.edge_coef(pair3).eq.0) goto 1000
c
		nlink= nlink_trig(idx)
		if(nlink.eq.0) then
			l1 = i
			l2 = i
		elseif(nlink.eq.1) then
			l1 = link_trig(1,idx)
			l2 = i
               	else
			l1 = link_trig(1,idx)
			l2 = link_trig(2,idx)
		endif
		wl1 = 0.5d0*weight(l1)
		wl2 = 0.5d0*weight(l2)
c
		flaga = edge_coef(pair1)
		flagb = edge_coef(pair2)
		flagc = edge_coef(pair3)
c
		do 750 l = 1,3
			x(l) = coord(3*l1-3+l)
			y(l) = coord(3*l2-3+l)
750		continue
c
		call threesurf_dir(a,b,c,x,y,
     1		nlink,pabc,pacb,eps,ra,rb,rc,
     2		ra2,wa,wb,wc,wl1,wl2,
     3		dist1,dist2,dist3,flaga,flagb,flagc,
     4		dsurfa3,dsurfb3,dsurfc3)
c
               	do 800 l = 1,3
                       	dsurf(l,i)=dsurf(l,i)+(coefa*
     1			dsurfa3(l,1)+coefb*dsurfb3(l,1)+
     2			coefc*dsurfc3(l,1))
                       	dsurf(l,j) = dsurf(l,j) + (coefa*
     1			dsurfa3(l,2)+ coefb*dsurfb3(l,2)+
     2			coefc*dsurfc3(l,2))
                       	dsurf(l,k) = dsurf(l,k) + (coefa*
     1			dsurfa3(l,3)+ coefb*dsurfb3(l,3)+
     2			coefc*dsurfc3(l,3))
800		continue
c
1000	continue
c
	do 1250 idx = 1,ntetra
c
		if(tetra_alp(idx).ne.1) goto 1250
c
		i = tetra(1,idx)
		j = tetra(2,idx)
		k = tetra(3,idx)
		l = tetra(4,idx)
c
		ra = radius(i)
		rb = radius(j)
		rc = radius(k)
		rd = radius(l)
c
		ra2 = radius2(i)
		rb2 = radius2(j)
		rc2 = radius2(k)
		rd2 = radius2(l)
c
		wa = 0.5d0*weight(i)
		wb = 0.5d0*weight(j)
		wc = 0.5d0*weight(k)
		wd = 0.5d0*weight(l)
c
		coefa = coefasp(i)
		coefb = coefasp(j)
		coefc = coefasp(k)
		coefd = coefasp(l)
c
		do 1050 m = 1,3
			a(m) = coord(3*(i-1)+m)
			b(m) = coord(3*(j-1)+m)
			c(m) = coord(3*(k-1)+m)
			d(m) = coord(3*(l-1)+m)
1050		continue
c
		trig1 = tetra_link(4,idx)
		trig2 = tetra_link(3,idx)
		trig3 = tetra_link(2,idx)
		trig4 = tetra_link(1,idx)
c
		pair1 = trig_link(3,trig1)
		pair2 = trig_link(2,trig1)
		pair4 = trig_link(1,trig1)
		pair3 = trig_link(2,trig2)
		pair5 = trig_link(1,trig2)
		pair6 = trig_link(1,trig3)
c
		dist1 = distpair(pair1)
		dist2 = distpair(pair2)
		dist3 = distpair(pair3)
		dist4 = distpair(pair4)
		dist5 = distpair(pair5)
		dist6 = distpair(pair6)
c
		d2_1 = distpair2(pair1)
		d2_2 = distpair2(pair2)
		d2_3 = distpair2(pair3)
		d2_4 = distpair2(pair4)
		d2_5 = distpair2(pair5)
		d2_6 = distpair2(pair6)
c
		eps1 = trig_eps(trig1)
		eps2 = trig_eps(trig2)
		eps3 = trig_eps(trig3)
		eps4 = trig_eps(trig4)
c
		sh_abc = trig_sh(1,trig1)
		sh_acb = trig_sh(2,trig1)
		sh_bca = trig_sh(3,trig1)
		sh_abd = trig_sh(1,trig2)
		sh_adb = trig_sh(2,trig2)
		sh_bda = trig_sh(3,trig2)
		sh_acd = trig_sh(1,trig3)
		sh_adc = trig_sh(2,trig3)
		sh_cda = trig_sh(3,trig3)
		sh_bcd = trig_sh(1,trig4)
		sh_bdc = trig_sh(2,trig4)
		sh_cdb = trig_sh(3,trig4)
c 
                do 1100 m = 1,3
			pabc(m)=trig_dual1(m,trig1)
			pacb(m)=trig_dual2(m,trig1)
			pabd(m)=trig_dual1(m,trig2)
			padb(m)=trig_dual2(m,trig2)
			pacd(m)=trig_dual1(m,trig3)
			padc(m)=trig_dual2(m,trig3)
			pbcd(m)=trig_dual1(m,trig4)
			pbdc(m)=trig_dual2(m,trig4)
1100		continue
c
		ang_abc(1) = trig_ang(1,trig1)
		ang_abc(2) = trig_ang(2,trig1)
		ang_abc(3) = trig_ang(3,trig1)
c
		ang_abd(1) = trig_ang(1,trig2)
		ang_abd(2) = trig_ang(2,trig2)
		ang_abd(3) = trig_ang(3,trig2)
c
		ang_acd(1) = trig_ang(1,trig3)
		ang_acd(2) = trig_ang(2,trig3)
		ang_acd(3) = trig_ang(3,trig3)
c
		ang_bcd(1) = trig_ang(1,trig4)
		ang_bcd(2) = trig_ang(2,trig4)
		ang_bcd(3) = trig_ang(3,trig4)
c
		test = tetra_orient(idx).eq.1
		if(test) then
			call foursphere_vol(a,b,c,d,ra,rb,rc,rd,
     1			ra2,rb2,rc2,rd2,dist1,dist2,dist3,dist4,
     2			dist5,dist6,d2_1,d2_2,d2_3,d2_4,d2_5,d2_6,
     3			wa,wb,wc,wd,eps1,eps2,eps3,eps4,
     4			sh_abc,sh_acb,sh_bca,sh_abd,sh_adb,sh_bda,
     5			sh_acd,sh_adc,sh_cda,sh_bcd,sh_bdc,sh_cdb,
     4			pacb,pabd,padc,pbcd,ang_abc,ang_abd,ang_acd,
     4			ang_bcd,surfa,surfb,surfc,surfd,vola,volb,volc,
     5			vold,dsurfa4,dsurfb4,dsurfc4,dsurfd4,
     6			dvola4,dvolb4,dvolc4,dvold4,option)
		else
			ang_acd(1) = trig_ang(2,trig3)
			ang_acd(2) = trig_ang(1,trig3)
			ang_acd(3) = trig_ang(3,trig3)
			ang_bcd(1) = trig_ang(2,trig4)
			ang_bcd(2) = trig_ang(1,trig4)
			ang_bcd(3) = trig_ang(3,trig4)
			call foursphere_vol(a,b,d,c,ra,rb,rd,rc,
     1			ra2,rb2,rd2,rc2,dist1,dist3,dist2,dist5,
     2			dist4,dist6,d2_1,d2_3,d2_2,d2_5,d2_4,d2_6,
     3			wa,wb,wd,wc,eps2,eps1,eps3,eps4,
     4			sh_abd,sh_adb,sh_bda,sh_abc,sh_acb,sh_bca,
     5			sh_adc,sh_acd,sh_cda,sh_bdc,sh_bcd,sh_cdb,
     4			padb,pabc,pacd,pbdc,ang_abd,ang_abc,ang_acd,
     7			ang_bcd,surfa,surfb,surfd,surfc,vola,volb,vold,
     8			volc,dsurfa4,dsurfb4,dsurfd4,dsurfc4,
     6			dvola4,dvolb4,dvold4,dvolc4,option)
		endif
c
		surf(i) = surf(i) - surfa
		surf(j) = surf(j) - surfb
		surf(k) = surf(k) - surfc
		surf(l) = surf(l) - surfd
c
		vol(i)  = vol(i) - vola
		vol(j)  = vol(j) - volb
		vol(k)  = vol(k) - volc
		vol(l)  = vol(l) - vold
c
		if(option.eq.0) goto 1250
c
		do 1150 m = 1,3
			dsurf(m,i) = dsurf(m,i) - coefa*dsurfa4(m,1) 
     1			-coefb*dsurfb4(m,1)-coefc*dsurfc4(m,1)-
     2			coefd*dsurfd4(m,1)
			dsurf(m,j) = dsurf(m,j) - coefa*dsurfa4(m,2) 
     1			-coefb*dsurfb4(m,2)-coefc*dsurfc4(m,2) 
     2			-coefd*dsurfd4(m,2)
			if(test) then
				dsurf(m,k) = dsurf(m,k) - coefa
     1				*dsurfa4(m,3)- coefb*dsurfb4(m,3)
     2				-coefc*dsurfc4(m,3)-coefd*dsurfd4(m,3)
				dsurf(m,l) = dsurf(m,l) - coefa
     1				*dsurfa4(m,4)- coefb*dsurfb4(m,4)
     2				-coefc*dsurfc4(m,4)-coefd*dsurfd4(m,4)
			else
				dsurf(m,k) = dsurf(m,k) - coefa
     1				*dsurfa4(m,4)- coefb*dsurfb4(m,4)
     2				-coefc*dsurfc4(m,4)-coefd*dsurfd4(m,4)
				dsurf(m,l) = dsurf(m,l) - coefa
     1				*dsurfa4(m,3)-coefb*dsurfb4(m,3)
     2				-coefc*dsurfc4(m,3)-coefd*dsurfd4(m,3)
			endif
1150		continue
c
		do 1200 m = 1,3
			dvol(m,i) = dvol(m,i) - coefa*dvola4(m,1) 
     1			-coefb*dvolb4(m,1)-coefc*dvolc4(m,1)
     2			-coefd*dvold4(m,1)
			dvol(m,j) = dvol(m,j) - coefa*dvola4(m,2) 
     1			-coefb*dvolb4(m,2)-coefc*dvolc4(m,2) 
     2			-coefd*dvold4(m,2)
			if(test) then
				dvol(m,k)=dvol(m,k)-coefa*dvola4(m,3)
     1				-coefb*dvolb4(m,3)-coefc*dvolc4(m,3)- 
     2				coefd*dvold4(m,3)
				dvol(m,l)=dvol(m,l)-coefa*dvola4(m,4)
     1				-coefb*dvolb4(m,4)-coefc*dvolc4(m,4)- 
     2				coefd*dvold4(m,4)
			else
				dvol(m,k) = dvol(m,k)-coefa*dvola4(m,4)
     1				- coefb*dvolb4(m,4)-coefc*dvolc4(m,4)- 
     2				coefd*dvold4(m,4)
				dvol(m,l) = dvol(m,l)-coefa*dvola4(m,3)
     1				-coefb*dvolb4(m,3)-coefc*dvolc4(m,3)- 
     2				coefd*dvold4(m,3)
			endif
1200		continue
c
1250	continue
c
	surftot = 0
	voltot  = 0
	wsurf   = 0
	wvol    = 0
c
	do 1300 i = 1,nvertex
		surftot = surftot + surf(i)
		wsurf  = wsurf + coefasp(i)*surf(i)
		voltot = voltot + vol(i)
		wvol   = wvol + coefasp(i)*vol(i)
1300	continue
c
c	Shift results, as first 4 points are the artificial infinite points
c
	do 1500 i = 1,npoints
		surf(i) = surf(i+4)
		vol(i) = vol(i+4)
		do 1400 j = 1,3
			dsurf(j,i) = dsurf(j,i+4)
			dvol(j,i) = dvol(j,i+4)
1400		continue
1500	continue
c
	return
	end
c
c	Alpha_result.f		Version 1 6/14/2005
c
c	This subroutine passes the volume and surface as arguments
c
	subroutine alpha_result(vtot,stot,wv,ws,volume,surface)
c
	integer	npointmax
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)

	parameter (npointmax=MAX_POINT)
c
	integer	i,npoints,nvertex
	integer	redinfo(npointmax)
c
	real*8	voltot,wvol,surftot,wsurf
	real*8	vtot,wv,stot,ws
	real*8	vol(npointmax),surf(npointmax)
	real*8	volume(npointmax),surface(npointmax)
	real*8	dvol(3,npointmax),dsurf(3,npointmax)
c
	common /results/	wsurf,wvol,surftot,voltot,surf,vol,
     1				dsurf,dvol
	common	/vertex_zone/	npoints,nvertex,redinfo
c
	vtot = voltot
	wv   = wvol
	stot = surftot
	ws   = wsurf
c
	do 100 i = 1,npoints
		volume(i) = vol(i)
		surface(i) = surf(i)

100	continue
c
	return
	end
c
c	Alpha_deriv.f		Version 1 6/14/2005
c
c	This subroutine passes the derivatives of surface area and volume
c
	subroutine alpha_deriv(dvolume,dsurface)
c
	integer	npointmax
	parameter(MAX_POINT = 20000)
	parameter(MAX_TETRA = 200000)
	parameter(MAX_TRIG = 400000)
	parameter(MAX_EDGE = 400000)

	parameter	(npointmax=MAX_POINT)
c
	integer	i,npoints,nvertex
	integer	redinfo(npointmax)
c
	real*8	voltot,wvol,surftot,wsurf
	real*8	vol(npointmax),surf(npointmax)
	real*8	dvolume(3*npointmax),dsurface(3*npointmax)
	real*8	dvol(3*npointmax),dsurf(3*npointmax)
c
	common /results/	wsurf,wvol,surftot,voltot,surf,vol,
     1				dsurf,dvol
	common	/vertex_zone/	npoints,nvertex,redinfo
c
	do 100 i = 1,3*npoints
		dvolume(i) = dvol(i)
		dsurface(i) = dsurf(i)
100	continue
c
	return
	end
