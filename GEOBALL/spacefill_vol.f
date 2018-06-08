c	spacefill_vol.f		Version 1 : 10/28/2002
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
c	center2.f
c
c	Copyright (C) 2002 Patrice Koehl
c
c	This subroutine finds the orthogonal center C of two spheres A and B
c	(if radii of A and B are equal, C is just the midpoint of the centers
c	of A and B)
c
	subroutine center2(a,b,ra2,rb2,rab2,c,lamda)
c
	integer	i
c
	real*8	ra2,rb2,rab2,lamda,uml
c
	real*8	a(3),b(3),c(3)
c
	lamda = 0.5d0 - (ra2-rb2)/(2*rab2)
	uml   = 1-lamda
c
	do 100 i = 1,3
		c(i) = lamda*a(i) + uml*b(i)
100	continue
c
	return
	end
c
c	Center3.f
c
c	Copyright (C) 2002 Patrice Koehl
c
c	This subroutine finds the orthogonal center Y of three spheres A and B
c	and C
c
	subroutine center3(a,b,c,i0,j0,k0,y)
c
	real*8	i0,j0,k0,a1,a2,a3,a4
	real*8	dx,dy,dz,d0
	real*8	val1,val2,val3
	real*8	a(3),b(3),c(3),y(3)
c
	a1=(b(2)-a(2))*(c(3)-a(3)) - (c(2)-a(2))*(b(3)-a(3))
	a2=(b(3)-a(3))*(c(1)-a(1)) - (c(3)-a(3))*(b(1)-a(1))
c
	val1=b(1)*c(2)-c(1)*b(2)
	val2=a(1)*c(2)-c(1)*a(2)
	val3=a(1)*b(2)-b(1)*a(2)
c
	a3 = val1-val2+val3
	a4 = a(3)*val1 - b(3)*val2 + c(3)*val3
c
	d0 = -a1*a1 - a2*a2 - a3*a3
c
	val1 = i0*(b(3)-c(3)) - j0*(a(3)-c(3)) + k0*(a(3)-b(3))
	val2 = i0*(b(2)-c(2)) - j0*(a(2)-c(2)) + k0*(a(2)-b(2))
	val3 = i0*(b(1)-c(1)) - j0*(a(1)-c(1)) + k0*(a(1)-b(1))
c
	dx = -a4*a1 + a2*val1 -a3*val2
	dy = -a1*val1 - a4*a2 + a3*val3
	dz = a1*val2 -a2*val3 -a4*a3
c
	y(1) = dx/d0
	y(2) = dy/d0
	y(3) = dz/d0
c
	return
	end
c	Center4.f
c
c	Copyright (C) 2002 Patrice Koehl
c
c	This subroutine finds the orthogonal center Y of four spheres A, B,
c	C and D
c
	subroutine center4(a,b,c,d,SA,SB,SC,SD,y)
c
	integer	i,j
c
	real*8	SA,SB,SC,SD,A1,B1,C1
	real*8	S(3),a(3),b(3),c(3),d(3),y(3),detval(4)
	real*8	ad(3),bd(3),cd(3)
c
	save
c
	S(1) = SA - SD
	S(2) = SB - SD
	S(3) = SC - SD
c
	do 200 i = 1,3
		ad(i) = a(i) - d(i)
		bd(i) = b(i) - d(i)
		cd(i) = c(i) - d(i)
200	continue
c
	A1 = bd(2)*cd(3)-cd(2)*bd(3)
	B1 = ad(2)*cd(3) - ad(3)*cd(2)
	C1 = ad(2)*bd(3) - ad(3)*bd(2)
c
	detval(1) = ad(1)*A1 - bd(1)*B1 + cd(1)*C1
	detval(2) = S(1)*A1 - S(2)*B1 + S(3)*C1
c
	A1 = S(2)*cd(1) - bd(1)*S(3)
	B1 = S(1)*cd(1) - S(3)*ad(1)
	C1 = S(1)*bd(1) - S(2)*ad(1)
c
	detval(3) = -ad(3)*A1 + bd(3)*B1 - cd(3)*C1
	detval(4) = ad(2)*A1 - bd(2)*B1 + cd(2)*C1
c
	if(detval(1).ne.0) then
c
		y(1) = detval(2)/detval(1)
		y(2) = detval(3)/detval(1)
		y(3) = detval(4)/detval(1)
c
	else
c
c		here we are in the special case that the four spheres
c		are co-planar, and the center3 of each subset of three
c		spheres are equal. We take this center3 as the "center4"
c
		call center3(a,b,c,SA,SB,SC,y)
c
	endif
c
	return
	end
c
c	det4_1.f
c
c	Copyright (C) 2002 Patrice Koehl
c
c	This function computes a 4 by 4 determinant, in which the last column
c	only contains the number 1
c
	function det4_1(mat4)
c
	real*8	det4_1,a,b,c,d
	real*8	val1,val2,val3,val4,val5,val6
	real*8 mat4(4,4)
c
	val1 = mat4(3,3) - mat4(4,3)
	val2 = mat4(2,3) - mat4(4,3)
	val3 = mat4(2,3) - mat4(3,3)
	val4 = mat4(1,3) - mat4(4,3)
	val5 = mat4(1,3) - mat4(3,3)
	val6 = mat4(1,3) - mat4(2,3)
c
	a = mat4(1,1)*
     1		(mat4(2,2)*val1
     2		- mat4(3,2)*val2
     3		+ mat4(4,2)*val3)
c
	b = mat4(2,1)*
     1		(mat4(1,2)*val1
     2		- mat4(3,2)*val4
     3		+ mat4(4,2)*val5)
c
	c = mat4(3,1)*
     1		(mat4(1,2)*val2
     2		-mat4(2,2)*val4
     3		+mat4(4,2)*val6)
c
	d = mat4(4,1)*
     1		(mat4(1,2)*val3
     2		- mat4(2,2)*val5
     3		+ mat4(3,2)*val6)
c
	det4_1 = a-b+c-d
c
	return
	end
c
c	triangle_dual.f
c
c	Copyright (C) 2002 Patrice Koehl
c
c	This subroutine computes the dual of a triangle (a,b,c), i.e.
c	the point P_ABC at which the 3 spheres centered at a, b and c
c	with radii ra, rb and rc intersect, and (a,b,c,P_ABC) has
c	positive orientation
c
	subroutine triangle_dual(a,b,c,center,eps2,ra2,d1,d2)
c
	integer	i
c
	real*8	ra2
	real*8	s2,s3,eps,eps2
	real*8	a(3),b(3),c(3),d1(3),d2(3)
	real*8	u1(3),u2(3),n(3),center(3)
c
	save
c
	call diffvect(b,a,u1)
	call diffvect(c,a,u2)
c
	call crossvect(u1,u2,n)
	call diffvect(center,a,u1)
c
	call dotvect(n,n,s2)
	call dotvect(u1,u1,s3)
c
	if( (ra2-s3).lt.0) then
	
		eps2 = sqrt(s3-ra2)
		
	else
		eps2 = sqrt(ra2-s3)
	
	end if

	eps = eps2/sqrt(s2)
c
	do 100 i = 1,3
		d1(i) = center(i) + eps*n(i)
		d2(i) = center(i) - eps*n(i)
100	continue

c
	return
	end
c
c	tetra_volume.f
c
c	Copyright (C) 2002 Patrice Koehl
c
c	This subroutine computes the volume of a tetrahedron
c
	function tetra_volume(a,b,c,d)
c
	integer	i
c
	real*8	det,det4_1
	real*8	tetra_volume
	real*8	a(3),b(3),c(3),d(3)
	real*8	mat4(4,4)
c
	save
c
	do 100 i = 1,4
		mat4(i,4) = 1.d0
100	continue
c
	do 200 i = 1,3
		mat4(1,i) = a(i)
		mat4(2,i) = b(i)
		mat4(3,i) = c(i)
		mat4(4,i) = d(i)
200	continue
c
	det = det4_1(mat4)
c
	tetra_volume = abs(det)/6
c
	return
	end
c
c	Distance2.f
c
c	Copyright (C) 2002 Patrice Koehl
c
c	This subroutine computes the square of the distance between two
c	sphere centers
c
	subroutine distance2(coord,n1,n2,dist,ncortot)
c
	integer	i,n1,n2,ncortot
c
	real*8	dist,val
	real*8	coord(ncortot)
c
	dist = 0
	do 100 i = 1,3
		val = coord(3*(n1-1)+i)-coord(3*(n2-1)+i)
		dist = dist + val*val
100	continue
c
	return
	end
c
c	twosphere_vol.f
c
c	This subroutine calculates the volume and surface of the
c	intersection of two spheres; it is only called when the
c	intersection exists
c	
c	Copyright (C) 2002 Patrice Koehl
c
c	It also provides the derivative of the surface and volume with
c	respect to the coordinates of the center of the atoms
c
	subroutine twosphere_vol(a,b,ra,ra2,rb,rb2,rab,rab2,
     1			surfa,surfb,vola,volb,
     2			dsurfa,dsurfb,dvola,dvolb,option)
c
c	Input:
c			a,b	: position of the centers of the 2 spheres
c			rab	: distance between the centers of the 2 spheres
c			rab2	: distance between the centers of the 2 spheres
c				  (squared)
c			ra,rb	: radii of sphere A and B, respectively
c			ra2 and rb2 are the squared of the quantities
c			above)
c			option	: 0 for surf and vol only, 1 if the derivatives
c				  are also computed
c	Output
c			surfa	: partial contribution of A to the total
c				  surface of the intersection
c			surfb	: partial contribution of B to the total
c				  surface of the intersection
c			vola	: partial contribution of A to the total
c				  volume of the intersection
c			volb	: partial contribution of B to the total
c				  volume of the intersection
c			dsurfa,dsurfb,dvola,dvolb: derivatives of surfa,
c						   surfb,vola and volb
c						   with respect to the
c						   coordinates of A and B
c
c
	integer	i,option
c
	real*8	ra,rb,surfa,surfb,vola,volb
	real*8	vala,valb,lamda
	real*8	ra2,rb2,rab,rab2,ha,hb,sa,ca,sb,cb
	real*8	dera,derb,der1,der2,coef1,coef2
	real*8	coefa,coefb,Aab
	real*8	pi,twopi
c
	real*8	a(3),b(3),c(3)
	real*8	u_ab(3),dsurfa(3,2),dsurfb(3,2)
	real*8	dvola(3,2),dvolb(3,2)
c
	common /constants/ pi,twopi
c
	save
c

c	Get "center" of the two spheres
c
	call center2(a,b,ra2,rb2,rab2,c,lamda)
c
	valb = lamda*rab
	vala = rab-valb
c
c	Get height of the cap of sphere A occluded by sphere B
c
	ha = ra - vala
c
c	same for sphere B ...
c
	hb = rb - valb
c
c	Get surfaces of intersection
c
	surfa = twopi*ra*ha
	surfb = twopi*rb*hb
c
c	Now get volume
c
	Aab = pi*(ra2-vala*vala)
c
	sa = ra*surfa
	ca = vala*Aab
c
	vola = (sa-ca)/3
c
	sb = rb*surfb
	cb = valb*Aab
c
	volb = (sb-cb)/3
c
c	Compute derivatives, if needed
c
	if(option.eq.0) return
c
	do 50 i = 1,3
		u_ab(i) = (a(i)-b(i))/rab
50	continue
c
	dera = -lamda
	derb = lamda-1
c
	der1 = twopi*ra*dera
	der2 = twopi*rb*derb
c
	do 100 i = 1,3
		coef1 = der1*u_ab(i)
		coef2 = der2*u_ab(i)
		dsurfa(i,1) = coef1
		dsurfa(i,2) = -coef1
		dsurfb(i,1) = coef2
		dsurfb(i,2) = -coef2
100	continue
c
	coefa = Aab*lamda
	coefb = Aab - coefa
c
	do 150 i = 1,3
		coef1 = -coefa*u_ab(i)
		coef2 = -coefb*u_ab(i)
		dvola(i,1) = coef1
		dvola(i,2) = -coef1
		dvolb(i,1) = coef2
		dvolb(i,2) = -coef2
150	continue
c
	return
	end
c
c	threevol_dist.f
c
c	Copyright (C) 2002 Patrice Koehl
c
c	This subroutine calculates the volume and surface of the intersection 
c	of three spheres; it is only called when the intersection exists
c
c	If needed, it also computes the derivatives of the surface area
c	and volume
c
	subroutine threevol_dist(a,b,c,ra,rb,rc,ra2,rb2,rc2,
     1		wa,wb,wc,rab,rac,rbc,rab2,rac2,rbc2,surfa,surfb,
     2		surfc,vola,volb,volc,dsurfa,dsurfb,dsurfc,
     3		dvola,dvolb,dvolc,dvolda,dvoldb,dvoldc,
     4		eps,pabc,pacb,angles,sh_abc,sh_acb,sh_bca,option)
c
c	Input:
c			ra,rb,rc  : radii of sphere A, B and C, respectively
c			ra2,rb2,rc2: radii squared
c			wa,wb,wc   : "weights" of A, B and C
c			(for example, wa = 0.5*(xa**2+ya**2+za**2-ra**2)
c			rab,rab2: distance between the centers of sphere A and B
c			rac,rac2: distance between the centers of sphere A and C
c			rbc,rbc2: distance between the centers of sphere B and C
c			option	: 0 for surface and vol only, 1 if derivatives
c				  are computed
c	Output
c			surfa,surfb,surfc : contribution of A, B and C to
c			the total surface of the intersection of A,B,C
c			vola,volb,volc : contribution of A, B and C to
c			the total volume of the intersection of A,B,C
c			dsurfa,dsurfb,dsurfc : derivatives of surfa,surfb,
c			and surfc with respect to the coordinates of A, B
c			and C
c			dvola,dvolb,dvolc : derivatives of vola,volb,
c			and volc with respect to the coordinates of A, B
c			and C
c
c	Other parameters are passed as output, and are then used by 
c	foursphere_vol
c
	integer	i,option
c
	real*8	pi,twopi
	real*8	surfa,surfb,surfc
	real*8	vola,volb,volc
	real*8	ra,rb,rc,rab,rac,rbc,rab2,rac2,rbc2
	real*8	ra2,rb2,rc2,wa,wb,wc
	real*8	a1,a2,a3,s2,c1,c2,eps
	real*8	seg_ang_abc,seg_ang_acb
	real*8	seg_ang_bca
	real*8	ang_abc,ang_acb,ang_bca
	real*8	cos_abc,cos_acb,cos_bca
	real*8	sin_abc,sin_acb,sin_bca
	real*8	ang_dih_abc,ang_dih_cab
	real*8	ang_dih_bac
	real*8	s_abc,s_acb,s_bca
	real*8	sh_abc,sh_acb,sh_bca
	real*8	val1,val2,val3,l1,l2,l3
	real*8	val1b,val2b,val3b
	real*8	coef,coef1,coef2,coef3
	real*8	rho_ab2,rho_ac2,rho_bc2
	real*8	der_ab,der_ac,der_bc,diff_ab,diff_ac,diff_bc
	real*8	dsurfa_ab,dsurfa_ac
	real*8	dsurfb_ab,dsurfb_bc
	real*8	dsurfc_ac,dsurfc_bc
	real*8	coef_ab,coef_ac,coef_bc
	real*8	a(3),b(3),c(3),center(3)
	real*8	c_ab(3),c_ac(3),c_bc(3)
	real*8	pabc(3),pacb(3)
	real*8	dsurfa(3,3),dsurfb(3,3),dsurfc(3,3)
	real*8	dvola(3,3),dvolb(3,3),dvolc(3,3)
	real*8	dvolda(3,3),dvoldb(3,3),dvoldc(3,3)
	real*8	u_ab(3),u_ac(3),u_bc(3)
	real*8	Vab(3),Vac(3),Vbc(3)
	real*8	cosine(3),sine(3),angles(3)
c
	common /constants/ pi,twopi
c
	save
c
	call center2(a,b,ra2,rb2,rab2,c_ab,l1)
	call center2(a,c,ra2,rc2,rac2,c_ac,l2)
	call center2(b,c,rb2,rc2,rbc2,c_bc,l3)
c
	val1 = l1*rab
	val2 = l2*rac
	val3 = l3*rbc
c
	val1b = rab-val1
	val2b = rac-val2
	val3b = rbc-val3
c
	call center3(a,b,c,wa,wb,wc,center)
c
	call triangle_dual(a,b,c,center,eps,ra2,pabc,pacb)
c
	call tetra3_noder(a,b,c,pabc,rab,rac,rbc,
     1	ra,rb,rc,seg_ang_abc,seg_ang_acb,
     2	seg_ang_bca,ang_dih_abc,ang_dih_bac,ang_dih_cab,cosine,
     3	sine)
c
	angles(1) = seg_ang_abc
	angles(2) = seg_ang_acb
	angles(3) = seg_ang_bca
c
	a1 = ra*(1-2*ang_dih_abc)
	a2 = 2*seg_ang_abc*val1b
	a3 = 2*seg_ang_acb*val2b
c
	surfa = twopi*ra*(a1 - a2 - a3)
c
	a1 = rb*(1-2*ang_dih_bac)
	a2 = 2*seg_ang_abc*val1
	a3 = 2*seg_ang_bca*val3b
c
	surfb = twopi*rb*(a1 - a2 - a3)
c
	a1 = rc*(1-2*ang_dih_cab)
	a2 = 2*seg_ang_acb*val2
	a3 = 2*seg_ang_bca*val3
c
	surfc = twopi*rc*(a1 - a2 - a3)
c
	
	ang_abc = twopi*seg_ang_abc
	ang_acb = twopi*seg_ang_acb
	ang_bca = twopi*seg_ang_bca
c
	cos_abc = cosine(1)
	sin_abc = sine(1)
	cos_acb = cosine(2)
	sin_acb = sine(2)
	cos_bca = cosine(3)
	sin_bca = sine(3)
c
	rho_ab2 = ra2 - val1b*val1b
	rho_ac2 = ra2 - val2b*val2b
	rho_bc2 = rb2 - val3b*val3b
c
	s_abc = rho_ab2*(ang_abc - sin_abc*cos_abc)
	s_acb = rho_ac2*(ang_acb - sin_acb*cos_acb)
	s_bca = rho_bc2*(ang_bca - sin_bca*cos_bca)
c
	s2 = ra*surfa
	c1 = val1b*s_abc
	c2 = val2b*s_acb
c
	vola = (s2 - c1 - c2)/3
c
	s2 = rb*surfb
	c1 = val1*s_abc
	c2 = val3b*s_bca
c
	volb = (s2 - c1 - c2)/3
c
	s2 = rc*surfc
	c1 = val2*s_acb
	c2 = val3*s_bca
c
	volc = (s2 - c1 - c2)/3
c
	sh_abc = eps*cos_abc/sin_abc
	sh_acb = eps*cos_acb/sin_acb
	sh_bca = eps*cos_bca/sin_bca
c
	if(option.eq.0) return
c
	do 50 i = 1,3
		u_ab(i) = (a(i)-b(i))/rab
		u_ac(i) = (a(i)-c(i))/rac
		u_bc(i) = (b(i)-c(i))/rbc
50	continue
c
	dsurfa_ab = -2*ra*ang_abc*l1 
	dsurfa_ac = - 2*ra*ang_acb*l2
c
	do 100 i = 1,3
		diff_ab = dsurfa_ab*u_ab(i)
		diff_ac = dsurfa_ac*u_ac(i)
		dsurfa(i,1) = diff_ab + diff_ac
		dsurfa(i,2) = -diff_ab
		dsurfa(i,3) = -diff_ac
100	continue
c
	dsurfb_ab = -2*rb*ang_abc*(1-l1)
	dsurfb_bc = -2*rb*ang_bca*l3
c
	do 200 i = 1,3
		diff_ab = dsurfb_ab*u_ab(i)
		diff_bc = dsurfb_bc*u_bc(i)
		dsurfb(i,1) = diff_ab 
		dsurfb(i,2) = -diff_ab + diff_bc
		dsurfb(i,3) = -diff_bc
200	continue
c
	dsurfc_ac = - 2*rc*ang_acb*(1-l2)
	dsurfc_bc = -2*rc*ang_bca*(1-l3)
c
	do 300 i = 1,3
		diff_ac = dsurfc_ac*u_ac(i)
		diff_bc = dsurfc_bc*u_bc(i)
		dsurfc(i,1) = diff_ac
		dsurfc(i,2) = diff_bc
		dsurfc(i,3) = -diff_ac - diff_bc
300	continue
c
	der_ab = -s_abc*l1
	der_ac = -s_acb*l2
	der_bc = 0.d0
c
	do 400 i = 1,3
		coef_ab = der_ab*u_ab(i)
		coef_ac = der_ac*u_ac(i)
		dvola(i,1) = coef_ab + coef_ac
		dvola(i,2) = -coef_ab
		dvola(i,3) = -coef_ac
400	continue
c
	der_ab = s_abc*(1-l1)
	der_ac = 0
	der_bc = -s_bca*l3
c
	do 500 i = 1,3
		coef_ab = der_ab*u_ab(i)
		coef_bc = der_bc*u_bc(i)
		dvolb(i,1) = -coef_ab
		dvolb(i,2) = coef_ab + coef_bc
		dvolb(i,3) = - coef_bc
500	continue
c
	der_ab = 0
	der_ac = s_acb*(1-l2)
	der_bc = s_bca*(1-l3)
c
	do 600 i = 1,3
		coef_ac = der_ac*u_ac(i)
		coef_bc = der_bc*u_bc(i)
		dvolc(i,1) = - coef_ac
		dvolc(i,2) = - coef_bc
		dvolc(i,3) = coef_ac + coef_bc
600	continue
c
	do 700 i = 1,3
		Vab(i) = center(i) - c_ab(i)
		Vac(i) = center(i) - c_ac(i)
		Vbc(i) = center(i) - c_bc(i)
700     continue
c
	coef = -2*(eps**3)/3
	coef1 = coef/(sh_abc*rab)
	coef2 = coef/(sh_acb*rac)
	coef3 = coef/(sh_bca*rbc)
c
	do 750 i = 1,3
		coef_ab = coef1*Vab(i)
		coef_ac = coef2*Vac(i)
		coef_bc = coef3*Vbc(i)
		dvolda(i,1) = coef_ab + coef_ac
		dvolda(i,2) = -coef_ab
		dvolda(i,3) = -coef_ac
		dvoldb(i,1) = -coef_ab
		dvoldb(i,2) = coef_ab + coef_bc
		dvoldb(i,3) = -coef_bc
		dvoldc(i,1) = -coef_ac
		dvoldc(i,2) = -coef_bc
		dvoldc(i,3) = coef_ac + coef_bc
750	continue
c
	return
	end
c
c	Foursphere_vol.f
c
c	This subroutine calculates the volume and surface area of the
c	intersection of four spheres; this intersection is supposed to exist
c
c	This routine assumes that the 4 points (a,b,c,d) are in ccw order
c
c	Copyright (C) 2002 Patrice Koehl
c
	subroutine foursphere_vol(a,b,c,d,ra,rb,rc,rd,
     1			ra2,rb2,rc2,rd2,rab,rac,rad,rbc,rbd,rcd,
     2			rab2,rac2,rad2,rbc2,rbd2,rcd2,wa,wb,wc,wd,
     3			eps1,eps3,eps5,
     4			eps7,shabc,shacb,shbca,shabd,shadb,shbda,
     6			shacd,shadc,shcda,shbcd,shbdc,shcdb,
     7			pacb,pabd,padc,pbcd,ang_abc,ang_abd,ang_acd,
     8			ang_bcd,surfa,surfb,surfc,surfd,vola,volb,volc,
     9			vold,dsurfa,dsurfb,dsurfc,dsurfd,
     9			dvola,dvolb,dvolc,dvold,option)
c
	integer	i,option
c
	real*8	pi,twopi
	real*8	ra,rb,rc,rd,ra2,rb2,rc2,rd2
	real*8	rab,rac,rad,rbc,rbd,rcd
	real*8	rab2,rac2,rad2,rbc2,rbd2,rcd2
	real*8	surfa,surfb,surfc,surfd
	real*8	wa,wb,wc,wd
	real*8	val_ab,val_ac,val_ad,val_bc,val_bd,val_cd
	real*8	val2_ab,val2_ac,val2_ad,val2_bc,val2_bd,val2_cd
	real*8	ang1,ang2,ang3,ang4,ang5,ang6
	real*8	l_ab,l_ac,l_ad,l_bc,l_bd,l_cd
	real*8	shabc,shabd,shacd,shbcd,shacb,shadb,shadc,shbdc
	real*8	shbca,shbda,shcda,shcdb
	real*8	dist1,dist3,dist5,dist7
	real*8	h1,h3,h5,h7
	real*8	eps1,eps3,eps5,eps7
	real*8	dab2,dac2,dad2,dbc2,dbd2,dcd2
	real*8	s1,t1,t2
	real*8	vola,volb,volc,vold
	real*8	coef
	real*8	der_ab,der_ac,der_ad,der_bc,der_bd,der_cd
	real*8	diff_ab,diff_ac,diff_ad,diff_bc,diff_bd,diff_cd
	real*8	cap_ab,cap_ac,cap_ad,cap_bc,cap_bd,cap_cd
	real*8  sin_ab,sin_ac,sin_ad,sin_bc,sin_bd,sin_cd
	real*8  cos_ab,cos_ac,cos_ad,cos_bc,cos_bd,cos_cd
c
	real*8	a(3),b(3),c(3),d(3)
	real*8	c_ab(3),c_ac(3),c_ad(3),c_bc(3),c_bd(3),c_cd(3)
	real*8	c_abcd(3)
	real*8	pacb(3),pabd(3),padc(3),pbcd(3)
	real*8  ang_abc(3),ang_abd(3),ang_acd(3),ang_bcd(3)
	real*8	dsurfa(3,4),dsurfb(3,4),dsurfc(3,4),dsurfd(3,4)
c
	real*8	dvola(3,4),dvolb(3,4),dvolc(3,4),dvold(3,4)
	real*8	u_ab(3),u_ac(3),u_ad(3),u_bc(3),u_bd(3),u_cd(3)
	real*8  v_abc(3),v_abd(3),v_acb(3),v_acd(3),v_adb(3),v_adc(3)
	real*8  v_bca(3),v_bcd(3),v_bda(3),v_bdc(3),v_cda(3),v_cdb(3)
	real*8  vect_ab(3),vect_ac(3),vect_ad(3),vect_bc(3),vect_bd(3)
	real*8  vect_cd(3)
	real*8  cof_ab(3),cof_ac(3),cof_ad(3),cof_bc(3),cof_bd(3)
	real*8  cof_cd(3)
c
	common /constants/ pi,twopi
c
	save
c
	call tetra_6dihed(a,b,c,d,ang1,ang2,ang4,ang3,ang5,ang6)
c
	call center2(a,b,ra2,rb2,rab2,c_ab,l_ab)
	call center2(a,c,ra2,rc2,rac2,c_ac,l_ac)
	call center2(a,d,ra2,rd2,rad2,c_ad,l_ad)
	call center2(b,c,rb2,rc2,rbc2,c_bc,l_bc)
	call center2(b,d,rb2,rd2,rbd2,c_bd,l_bd)
	call center2(c,d,rc2,rd2,rcd2,c_cd,l_cd)
c
	val_ab = l_ab*rab
	val_ac = l_ac*rac
	val_ad = l_ad*rad
	val_bc = l_bc*rbc
	val_bd = l_bd*rbd
	val_cd = l_cd*rcd
c
	val2_ab = rab - val_ab
	val2_ac = rac - val_ac
	val2_ad = rad - val_ad
	val2_bc = rbc - val_bc
	val2_bd = rbd - val_bd
	val2_cd = rcd - val_cd
c
	surfa = -0.5d0*ra + ang1*val2_ab + ang2*val2_ac +
     1		ang3*val2_ad
	surfa = twopi*ra*surfa
c
	surfb = -0.5d0*rb + ang1*val_ab + ang5*val2_bd +
     1		ang4*val2_bc
	surfb = twopi*rb*surfb
c
	surfc = -0.5d0*rc + ang2*val_ac + ang4*val_bc +
     1		ang6*val2_cd
	surfc = twopi*rc*surfc
c
	surfd = -0.5d0*rd + ang3*val_ad + ang6*val_cd +
     1		ang5*val_bd
	surfd = twopi*rd*surfd
c
c	Now computes volume
c
	call center4(a,b,c,d,wa,wb,wc,wd,c_abcd)
c
	dab2 = ra2 - val2_ab*val2_ab
	dac2 = ra2 - val2_ac*val2_ac
	dad2 = ra2 - val2_ad*val2_ad
	dbc2 = rb2 - val2_bc*val2_bc
	dbd2 = rb2 - val2_bd*val2_bd
	dcd2 = rc2 - val2_cd*val2_cd
c
	dist1 = 0
	dist3 = 0
	dist5 = 0
	dist7 = 0
	do 50 i = 1,3
		dist1 = dist1 + (pacb(i)-c_abcd(i))**2
		dist3 = dist3 + (pabd(i)-c_abcd(i))**2
		dist5 = dist5 + (padc(i)-c_abcd(i))**2
		dist7 = dist7 + (pbcd(i)-c_abcd(i))**2
50	continue
	dist1 = sqrt(dist1)
	dist3 = sqrt(dist3)
	dist5 = sqrt(dist5)
	dist7 = sqrt(dist7)
c
	h1 = dist1-eps1
	h3 = dist3-eps3
	h5 = dist5-eps5
	h7 = dist7-eps7
c
	s1 = -twopi*dab2*ang1
	t1 = shabc*h1
	t2 = shabd*h3
c
	cap_ab = s1 -t1 -t2
c
	s1 = -twopi*dac2*ang2
	t1 = shacd*h5
	t2 = shacb*h1
c
	cap_ac = s1 -t1 -t2
c
	s1 = -twopi*dad2*ang3
	t1 = shadb*h3
	t2 = shadc*h5
c
	cap_ad = s1 - t1 -t2
c
	s1 = -twopi*dbc2*ang4
	t1 = shbca*h1
	t2 = shbcd*h7
c
	cap_bc = s1 - t1 -t2
c
	s1 = -twopi*dbd2*ang5
	t1 = shbdc*h7
	t2 = shbda*h3
c
	cap_bd = s1 - t1 -t2
c
	s1 = -twopi*dcd2*ang6
	t1 = shcda*h5
	t2 = shcdb*h7
c
	cap_cd = s1 - t1 -t2
c
	vola = 2*ra*surfa-val2_ab*cap_ab-val2_ac*cap_ac-val2_ad*cap_ad
	vola = vola/6
c
	volb = 2*rb*surfb - val_ab*cap_ab-val2_bd*cap_bd-val2_bc*cap_bc
	volb = volb /6
c
	volc = 2*rc*surfc -val_ac*cap_ac-val_bc*cap_bc-val2_cd*cap_cd
	volc = volc/6
c
	vold = 2*rd*surfd-val_ad*cap_ad-val_bd*cap_bd-val_cd*cap_cd
	vold = vold/6
c
	if(option.eq.0) return
c
	do 70 i = 1,3
		u_ab(i) = (a(i)-b(i))/rab
		u_ac(i) = (a(i)-c(i))/rac
		u_ad(i) = (a(i)-d(i))/rad
		u_bc(i) = (b(i)-c(i))/rbc
		u_bd(i) = (b(i)-d(i))/rbd
		u_cd(i) = (c(i)-d(i))/rcd
70	continue
c
	coef = twopi*ra
c
	der_ab = coef*ang1*l_ab
	der_ac = coef*ang2*l_ac
	der_ad = coef*ang3*l_ad
c
	do 100 i = 1,3
		diff_ab = der_ab*u_ab(i)
		diff_ac = der_ac*u_ac(i)
		diff_ad = der_ad*u_ad(i)
		dsurfa(i,1) = diff_ab+diff_ac+ diff_ad
		dsurfa(i,2) = -diff_ab
		dsurfa(i,3) = -diff_ac
		dsurfa(i,4) = -diff_ad
100	continue
c
	coef = twopi*rb
	der_ab = coef*ang1*(1-l_ab)
	der_bc = coef*ang4*l_bc
	der_bd = coef*ang5*l_bd
c
	do 200 i = 1,3
		diff_ab = der_ab*u_ab(i)
		diff_bc = der_bc*u_bc(i)
		diff_bd = der_bd*u_bd(i)
		dsurfb(i,1) = diff_ab
		dsurfb(i,2) = -diff_ab+diff_bc+diff_bd
		dsurfb(i,3) = -diff_bc
		dsurfb(i,4) = -diff_bd
200	continue
c
	coef = twopi*rc
	der_ac = coef*ang2*(1-l_ac)
	der_bc = coef*ang4*(1-l_bc)
	der_cd = coef*ang6*l_cd
c
	do 300 i = 1,3
		diff_ac = der_ac*u_ac(i)
		diff_bc = der_bc*u_bc(i)
		diff_cd = der_cd*u_cd(i)
		dsurfc(i,1) = diff_ac
		dsurfc(i,2) = diff_bc
		dsurfc(i,3) = -diff_ac-diff_bc+diff_cd
		dsurfc(i,4) = -diff_cd
300	continue
c
	coef = twopi*rd
	der_ad = coef*ang3*(1-l_ad)
	der_bd = coef*ang5*(1-l_bd)
	der_cd = coef*ang6*(1-l_cd)
c
	do 400 i = 1,3
		diff_ad = der_ad*u_ad(i)
		diff_bd = der_bd*u_bd(i)
		diff_cd = der_cd*u_cd(i)
		dsurfd(i,1) = diff_ad
		dsurfd(i,2) = diff_bd
		dsurfd(i,3) = diff_cd
		dsurfd(i,4) = -diff_ad-diff_bd-diff_cd
400	continue
c
c	Now get volume derivatives
c
	do 500 i = 1,3
		v_abc(i) = (c_abcd(i)+pacb(i) - 2*c_ab(i))/2
		v_abd(i) = (c_abcd(i)+pabd(i) - 2*c_ab(i))/2
		v_acb(i) = (c_abcd(i)+pacb(i) - 2*c_ac(i))/2
		v_acd(i) = (c_abcd(i)+padc(i) - 2*c_ac(i))/2
		v_adb(i) = (c_abcd(i)+pabd(i) - 2*c_ad(i))/2
		v_adc(i) = (c_abcd(i)+padc(i) - 2*c_ad(i))/2
		v_bca(i) = (c_abcd(i)+pacb(i) - 2*c_bc(i))/2
		v_bcd(i) = (c_abcd(i)+pbcd(i) - 2*c_bc(i))/2
		v_bda(i) = (c_abcd(i)+pabd(i) - 2*c_bd(i))/2
		v_bdc(i) = (c_abcd(i)+pbcd(i) - 2*c_bd(i))/2
		v_cda(i) = (c_abcd(i)+padc(i) - 2*c_cd(i))/2
		v_cdb(i) = (c_abcd(i)+pbcd(i) - 2*c_cd(i))/2
		vect_ab(i) = (pacb(i)+pabd(i)-2*c_ab(i))/2
		vect_ac(i) = (pacb(i)+padc(i)-2*c_ac(i))/2
		vect_ad(i) = (pabd(i)+padc(i)-2*c_ad(i))/2
		vect_bc(i) = (pacb(i)+pbcd(i)-2*c_bc(i))/2
		vect_bd(i) = (pabd(i)+pbcd(i)-2*c_bd(i))/2
		vect_cd(i) = (padc(i)+pbcd(i)-2*c_cd(i))/2
500	continue
c
	sin_ab = sin(pi*(ang_abc(1)+ang_abd(1)-ang1))
	sin_ac = sin(pi*(ang_abc(2)+ang_acd(1)-ang2))
	sin_ad = sin(pi*(ang_abd(2)+ang_acd(2)-ang3))
	sin_bc = sin(pi*(ang_abc(3)+ang_bcd(1)-ang4))
	sin_bd = sin(pi*(ang_abd(3)+ang_bcd(2)-ang5))
	sin_cd = sin(pi*(ang_acd(3)+ang_bcd(3)-ang6))
c
	cos_ab = cos(pi*(ang_abc(1)+ang_abd(1)-ang1))
	cos_ac = cos(pi*(ang_abc(2)+ang_acd(1)-ang2))
	cos_ad = cos(pi*(ang_abd(2)+ang_acd(2)-ang3))
	cos_bc = cos(pi*(ang_abc(3)+ang_bcd(1)-ang4))
	cos_bd = cos(pi*(ang_abd(3)+ang_bcd(2)-ang5))
	cos_cd = cos(pi*(ang_acd(3)+ang_bcd(3)-ang6))
c
        do 600 i = 1,3
		cof_ab(i) = (shabc*dist1*v_abc(i) + shabd*dist3*
     1		v_abd(i)-2*dab2*sin_ab*vect_ab(i)/cos_ab)/(3*rab)
		cof_ac(i) = (shacb*dist1*v_acb(i) + shacd*dist5*
     1		v_acd(i)-2*dac2*sin_ac*vect_ac(i)/cos_ac)/(3*rac)
		cof_ad(i) = (shadb*dist3*v_adb(i) + shadc*dist5*
     1		v_adc(i)-2*dad2*sin_ad*vect_ad(i)/cos_ad)/(3*rad)
		cof_bc(i) = (shbca*dist1*v_bca(i) + shbcd*dist7*
     1		v_bcd(i)-2*dbc2*sin_bc*vect_bc(i)/cos_bc)/(3*rbc)
		cof_bd(i) = (shbda*dist3*v_bda(i) + shbdc*dist7*
     1		v_bdc(i)-2*dbd2*sin_bd*vect_bd(i)/cos_bd)/(3*rbd)
		cof_cd(i) = (shcda*dist5*v_cda(i) + shcdb*dist7*
     1		v_cdb(i)-2*dcd2*sin_cd*vect_cd(i)/cos_cd)/(3*rcd)
600     continue
c
	der_ab = -0.5d0*l_ab*cap_ab
	der_ac = -0.5d0*l_ac*cap_ac
	der_ad = -0.5d0*l_ad*cap_ad
c
	do 650 i = 1,3
		diff_ab = der_ab*u_ab(i)
		diff_ac = der_ac*u_ac(i)
		diff_ad = der_ad*u_ad(i)
		dvola(i,1) = diff_ab+diff_ac+ diff_ad
		dvola(i,2) = -diff_ab
		dvola(i,3) = -diff_ac
		dvola(i,4) = -diff_ad
650	continue
c
        do 700 i = 1,3
		dvola(i,1)=dvola(i,1)+cof_ab(i)+cof_ac(i)+cof_ad(i)
		dvola(i,2)=dvola(i,2)-cof_ab(i)
		dvola(i,3)=dvola(i,3)-cof_ac(i)
		dvola(i,4)=dvola(i,4)-cof_ad(i)
700     continue
c
	der_ab = -0.5d0*cap_ab*(1-l_ab)
	der_bc = -0.5d0*cap_bc*l_bc
	der_bd = -0.5d0*cap_bd*l_bd
c
	do 750 i = 1,3
		diff_ab = der_ab*u_ab(i)
		diff_bc = der_bc*u_bc(i)
		diff_bd = der_bd*u_bd(i)
		dvolb(i,1) = diff_ab
		dvolb(i,2) = -diff_ab+diff_bc+diff_bd
		dvolb(i,3) = -diff_bc
		dvolb(i,4) = -diff_bd
750	continue
c
        do 800 i = 1,3
		dvolb(i,1)=dvolb(i,1)-cof_ab(i)
		dvolb(i,2)=dvolb(i,2)+cof_ab(i)+cof_bc(i)+cof_bd(i)
		dvolb(i,3)=dvolb(i,3)-cof_bc(i)
		dvolb(i,4)=dvolb(i,4)-cof_bd(i)
800     continue
c
	der_ac = -0.5d0*cap_ac*(1-l_ac)
	der_bc = -0.5d0*cap_bc*(1-l_bc)
	der_cd = -0.5d0*cap_cd*l_cd
c
	do 850 i = 1,3
		diff_ac = der_ac*u_ac(i)
		diff_bc = der_bc*u_bc(i)
		diff_cd = der_cd*u_cd(i)
		dvolc(i,1) = diff_ac
		dvolc(i,2) = diff_bc
		dvolc(i,3) = -diff_ac-diff_bc+diff_cd
		dvolc(i,4) = -diff_cd
850	continue
c
        do 900 i = 1,3
                dvolc(i,1)=dvolc(i,1)-cof_ac(i)
                dvolc(i,2)=dvolc(i,2)-cof_bc(i)
                dvolc(i,3)=dvolc(i,3)+cof_ac(i)+cof_bc(i)+cof_cd(i)
                dvolc(i,4)=dvolc(i,4)-cof_cd(i)
900     continue
c
	der_ad = -0.5d0*cap_ad*(1-l_ad)
	der_bd = -0.5d0*cap_bd*(1-l_bd)
	der_cd = -0.5d0*cap_cd*(1-l_cd)
c
	do 950 i = 1,3
		diff_ad = der_ad*u_ad(i)
		diff_bd = der_bd*u_bd(i)
		diff_cd = der_cd*u_cd(i)
		dvold(i,1) = diff_ad
		dvold(i,2) = diff_bd
		dvold(i,3) = diff_cd
		dvold(i,4) = -diff_ad-diff_bd-diff_cd
950	continue
c
        do 1000 i = 1,3
		dvold(i,1)=dvold(i,1)-cof_ad(i)
		dvold(i,2)=dvold(i,2)-cof_bd(i)
		dvold(i,3)=dvold(i,3)-cof_cd(i)
		dvold(i,4)=dvold(i,4)+cof_ad(i)+cof_bd(i)+cof_cd(i)
1000     continue
c
	return
	end
c
c	threesurf_dir.f
c
c	Copyright (C) 2002 Patrice Koehl
c
c	This subroutine calculates the correction term for the derivative 
c	of the surface area
c
	subroutine threesurf_dir(a,b,c,d1,d2,nlink,pabc,pacb,eps,
     1		ra,rb,rc,ra2,wa,wb,wc,wd1,wd2,rab,rac,rbc,
     2		flag_ab,flag_ac,flag_bc,dsurfa,dsurfb,dsurfc)
c
c	Input:
c			rab: distance between the centers of 
c				sphere A and B
c			rac: distance between the centers of 
c				sphere A and C
c			rbc: distance between the centers of 
c				sphere B and C
c
	integer	i,j,nlink
	integer	flag_ab,flag_ac,flag_bc
	integer	flag_a,flag_b,flag_c
c
	real*8	pi,twopi
	real*8	rab,rac,rbc
	real*8	ra,rb,rc,ra2,wa,wb,wc,wd1,wd2
	real*8	eps,beta,beta1,beta2
	real*8	coef2,coef3,coef_ab,coef_ac,coef_bc
	real*8	c_ab_ac,c_ab_bc,c_ac_bc
	real*8	a(3),b(3),c(3),d1(3),d2(3),center(3)
	real*8	pabc(3),pacb(3)
	real*8	dsurfa(3,3),dsurfb(3,3),dsurfc(3,3)
	real*8	u_ab(3),u_ac(3),u_bc(3)
	real*8	e_abc(3),e_cab(3),e_bca(3)
	real*8	u_abc(3),u_cab(3),u_bca(3)
c
	common /constants/ pi,twopi
c
	save
c
	do 100 i =1,3
		do 50 j = 1,3
			dsurfa(j,i) = 0
			dsurfb(j,i) = 0
			dsurfc(j,i) = 0
50		continue
100	continue
c
	flag_a = ior(flag_ab,flag_ac)
	flag_b = ior(flag_ab,flag_bc)
	flag_c = ior(flag_ac,flag_bc)
c
	if(nlink.eq.2) then
		call center3(a,b,c,wa,wb,wc,center)
		call triangle_dual(a,b,c,center,eps,ra2,pabc,pacb)
	endif
c
	if(nlink.eq.1) then
		call segment_cover(a,pabc,pacb,d1,wa,wd1,beta)
		beta = (1-beta)
	elseif(nlink.eq.2) then
		call segment_cover(a,pabc,pacb,d1,wa,wd1,beta1)
		call segment_cover(a,pabc,pacb,d2,wa,wd2,beta2)
		beta = 1-(beta2+beta1)
	else
		beta = 1
	endif
c
	beta = beta*2*eps
c
	if(flag_a.eq.0) then
c
		do 150 i = 1,3
			u_ab(i) = (a(i)-b(i))/rab
			u_bc(i) = (b(i)-c(i))/rbc
150		continue
c
		call dotvect(u_ab,u_bc,c_ab_bc)
		do 200 i = 1,3
			e_bca(i) = -u_ab(i) + c_ab_bc*u_bc(i)
200		continue
		call unitvector(e_bca,u_bca)
c
		coef3 = beta*rb/rbc
c
		do 300 i = 1,3
			coef_bc = coef3*u_bca(i)
			dsurfb(i,1) = 0
			dsurfb(i,2) = coef_bc
			dsurfb(i,3) = -coef_bc
300		continue
c
		coef3 = beta*rc/rbc
c
		do 350 i = 1,3
			coef_bc = coef3*u_bca(i)
			dsurfc(i,1) = 0
			dsurfc(i,2) = -coef_bc
			dsurfc(i,3) = coef_bc
350		continue
c
		return
c
	endif
c
	if(flag_b.eq.0) then
c
		do 400 i = 1,3
			u_ac(i) = (a(i)-c(i))/rac
			u_bc(i) = (b(i)-c(i))/rbc
400		continue
c
		call dotvect(u_ac,u_bc,c_ac_bc)
		do 450 i = 1,3
			e_cab(i) = -u_bc(i) + c_ac_bc*u_ac(i)
450		continue
		call unitvector(e_cab,u_cab)
c
		coef3 = beta*ra/rac
c
		do 500 i = 1,3
			coef_ac = coef3*u_cab(i)
			dsurfa(i,1) = coef_ac
			dsurfa(i,2) = 0
			dsurfa(i,3) = -coef_ac
500		continue
c
		coef2 = beta*rc/rac
c
		do 550 i = 1,3
			coef_ac = coef2*u_cab(i)
			dsurfc(i,1) = -coef_ac
			dsurfc(i,2) = 0
			dsurfc(i,3) = coef_ac
550		continue
c
		return
c
	endif
c
	if(flag_c.eq.0) then
c
		do 600 i = 1,3
			u_ab(i) = (a(i)-b(i))/rab
			u_ac(i) = (a(i)-c(i))/rac
600		continue
c
		call dotvect(u_ab,u_ac,c_ab_ac)
		do 650 i = 1,3
			e_abc(i) = u_ac(i) - c_ab_ac*u_ab(i)
650		continue
		call unitvector(e_abc,u_abc)
c
		coef2 = beta*ra/rab
c
		do 700 i = 1,3
			coef_ab = coef2*u_abc(i)
			dsurfa(i,1) = coef_ab
			dsurfa(i,2) = -coef_ab
			dsurfa(i,3) = 0
700		continue
c
		coef2 = beta*rb/rab
c
		do 750 i = 1,3
			coef_ab = coef2*u_abc(i)
			dsurfb(i,1) = -coef_ab
			dsurfb(i,2) = coef_ab
			dsurfb(i,3) = 0
750		continue
c
		return
c
	endif
c
	do 800 i = 1,3
		u_ab(i) = (a(i)-b(i))/rab
		u_ac(i) = (a(i)-c(i))/rac
		u_bc(i) = (b(i)-c(i))/rbc
800	continue
c
	call dotvect(u_ab,u_ac,c_ab_ac)
	call dotvect(u_ab,u_bc,c_ab_bc)
	call dotvect(u_ac,u_bc,c_ac_bc)
	do 850 i = 1,3
		e_abc(i) = u_ac(i) - c_ab_ac*u_ab(i)
		e_bca(i) = -u_ab(i) + c_ab_bc*u_bc(i)
		e_cab(i) = -u_bc(i) + c_ac_bc*u_ac(i)
850	continue
	call unitvector(e_abc,u_abc)
	call unitvector(e_bca,u_bca)
	call unitvector(e_cab,u_cab)
c
	coef2 = flag_ab*beta*ra/rab
	coef3 = flag_ac*beta*ra/rac
c
	do 900 i = 1,3
		coef_ab = coef2*u_abc(i)
		coef_ac = coef3*u_cab(i)
		dsurfa(i,1) = coef_ab + coef_ac
		dsurfa(i,2) = -coef_ab
		dsurfa(i,3) = -coef_ac
900	continue
c
	coef2 = flag_ab*beta*rb/rab
	coef3 = flag_bc*beta*rb/rbc
c
	do 950 i = 1,3
		coef_ab = coef2*u_abc(i)
		coef_bc = coef3*u_bca(i)
		dsurfb(i,1) = -coef_ab
		dsurfb(i,2) = coef_ab + coef_bc
		dsurfb(i,3) = -coef_bc
950	continue
c
	coef2 = flag_ac*beta*rc/rac
	coef3 = flag_bc*beta*rc/rbc
c
	do 1000 i = 1,3
		coef_ac = coef2*u_cab(i)
		coef_bc = coef3*u_bca(i)
		dsurfc(i,1) = -coef_ac
		dsurfc(i,2) = -coef_bc
		dsurfc(i,3) = coef_ac +coef_bc
1000	continue
c
	return
	end
c
c	Segment_cover.f		Version 1 7/10/2001	Patrice Koehl
c
c	This subroutine computes the fraction of the line segment 
c	P1P2 in	the bisector of a sphere D
c
	subroutine segment_cover(a,p1,p2,d,wa,wd,beta)
c
	real*8	beta,wa,wd
	real*8	val1,val2,test
	real*8	a(3),p1(3),p2(3),d(3),ad(3),p1p2(3)
c
	save
c
	call diffvect(p1,p2,p1p2)
	call diffvect(a,d,ad)
	call dotvect(ad,p1p2,val1)
	call dotvect(ad,p2,val2)
c
	if(val1.eq.0) then
c		write(6,*) 'The problem is in segment_cover...'
		test = val1 - 2*(val2-wd)
		if(test.gt.0) then
			beta = 0
		else
			beta = 1
		endif
		return
	endif
	beta = (wa-wd+val2)/val1
c
	call dotvect(p1,p1,val1)
	call dotvect(p1,d,val2)
c
	test = val1 - 2*(val2-wd)
c
	if(test.lt.0) then
		beta = 1-beta
	endif
c
	return
	end
c
c	threevol_dir.f
c
c	Copyright (C) 2002 Patrice Koehl
c
c	This subroutine computes the contribution of distance preserving
c	motions to the derivative of the volume of a union of balls
c
	subroutine threevol_dir(a,b,c,ra,rb,rc,ra2,rb2,rc2,
     1		wa,wb,wc,rab,rac,rbc,rab2,rac2,rbc2,
     2		dvola,dvolb,dvolc,angles,pabc,pacb,eps,sh_abc,sh_acb,
     3		sh_bca,option)
c
c	Input:
c			ra,rb,rc  : radii of sphere A, B and C, respectively
c			ra2,rb2,rc2: radii squared
c			rab,rab2: distance between the centers of sphere A and B
c			rac,rac2: distance between the centers of sphere A and C
c			rbc,rbc2: distance between the centers of sphere B and C
c	Output
c			dvola,dvolb,dvolc : contribution of A, B and C to
c			the derivatives of the volume of the intersection of A,B,C
c
	integer	i,option
c
	real*8	pi,twopi
	real*8	ra,rb,rc,rab,rac,rbc,rab2,rac2,rbc2
	real*8	ra2,rb2,rc2,wa,wb,wc
	real*8	eps
	real*8	seg_ang_abc,seg_ang_acb
	real*8	seg_ang_bca
	real*8	cos_abc,cos_acb,cos_bca
	real*8	sin_abc,sin_acb,sin_bca
	real*8	ang_dih_abc,ang_dih_cab
	real*8	ang_dih_bac
	real*8	sh_abc,sh_acb,sh_bca
	real*8	l1,l2,l3
	real*8	coef_ab,coef_ac,coef_bc
	real*8	coef,coef1,coef2,coef3
	real*8	a(3),b(3),c(3),center(3)
	real*8	c_ab(3),c_ac(3),c_bc(3)
	real*8	pabc(3),pacb(3)
	real*8	dvola(3,3),dvolb(3,3),dvolc(3,3)
	real*8	Vab(3),Vac(3),Vbc(3)
	real*8	angles(3),cosine(3),sine(3)
c
	common /constants/ pi,twopi
c
	save
c
	call center2(a,b,ra2,rb2,rab2,c_ab,l1)
	call center2(a,c,ra2,rc2,rac2,c_ac,l2)
	call center2(b,c,rb2,rc2,rbc2,c_bc,l3)
c
	call center3(a,b,c,wa,wb,wc,center)
c
	call triangle_dual(a,b,c,center,eps,ra2,pabc,pacb)
c
	call tetra3_noder(a,b,c,pabc,rab,rac,rbc,ra,rb,rc,
     1	seg_ang_abc,seg_ang_acb,seg_ang_bca,ang_dih_abc,
     2	ang_dih_bac,ang_dih_cab,cosine,sine)
c
	angles(1) = seg_ang_abc
	angles(2) = seg_ang_acb
	angles(3) = seg_ang_bca
c
	cos_abc = cosine(1)
	sin_abc = sine(1)
	cos_acb = cosine(2)
	sin_acb = sine(2)
	cos_bca = cosine(3)
	sin_bca = sine(3)
c
	sh_abc = eps*cos_abc/sin_abc
	sh_acb = eps*cos_acb/sin_acb
	sh_bca = eps*cos_bca/sin_bca
c
	if(option.eq.0) return
c
	do 700 i = 1,3
		Vab(i) = center(i) - c_ab(i)
		Vac(i) = center(i) - c_ac(i)
		Vbc(i) = center(i) - c_bc(i)
700	continue
c
	coef = -2*(eps**3)/3
	coef1 = coef/(sh_abc*rab)
	coef2 = coef/(sh_acb*rac)
	coef3 = coef/(sh_bca*rbc)
c
	do 750 i = 1,3
		coef_ab = coef1*Vab(i)
		coef_ac = coef2*Vac(i)
		coef_bc = coef3*Vbc(i)
		dvola(i,1) = coef_ab + coef_ac
		dvola(i,2) = -coef_ab
		dvola(i,3) = -coef_ac
		dvolb(i,1) = -coef_ab 
		dvolb(i,2) = coef_ab + coef_bc
		dvolb(i,3) = -coef_bc
		dvolc(i,1) = -coef_ac 
		dvolc(i,2) = -coef_bc 
		dvolc(i,3) = coef_ac + coef_bc
750	continue
c
	return
	end
c
c	tetra_6dihed.f
c
c	Copyright (C) 2002 Patrice Koehl
c
c	Get all six dihedral angles
c
c	ang1 = angle_dihed(a,b,c,d) = angle_dihed(T3,T4)
c	ang2 = angle_dihed(a,c,b,d) = angle_dihed(T2,T4)
c	ang3 = angle_dihed(b,c,a,d) = angle_dihed(T1,T4)
c	ang4 = angle_dihed(a,d,b,c) = angle_dihed(T2,T3)
c	ang5 = angle_dihed(b,d,a,c) = angle_dihed(T1,T3)
c	ang6 = angle_dihed(c,d,a,b) = angle_dihed(T1,T2)
c
	subroutine tetra_6dihed(a,b,c,d,ang1,ang2,ang3,ang4,ang5,ang6)
c
	real*8	ang1,ang2,ang3,ang4,ang5,ang6
	real*8	cos_ang
	real*8	pi,twopi
c
	real*8  a(3),b(3),c(3),d(3)
	real*8	u_ab(3),u_ac(3),u_ad(3),u_bc(3),u_bd(3)
c
	real*8	u_abc(3),u_abd(3),u_acd(3),u_bcd(3)
	real*8	n_abc(3),n_abd(3),n_acd(3),n_bcd(3)
c
	common /constants/ pi,twopi
c
	save
c
	call diffvect(a,b,u_ab)
	call diffvect(a,c,u_ac)
	call diffvect(a,d,u_ad)
	call diffvect(b,c,u_bc)
	call diffvect(b,d,u_bd)
c
	call crossvect(u_ab,u_ac,u_abc)
	call crossvect(u_ab,u_ad,u_abd)
	call crossvect(u_ac,u_ad,u_acd)
	call crossvect(u_bc,u_bd,u_bcd)
c
	call unitvector(u_abc,n_abc)
	call unitvector(u_abd,n_abd)
	call unitvector(u_acd,n_acd)
	call unitvector(u_bcd,n_bcd)
c
	call dotvect(n_abc,n_abd,cos_ang)
	ang1 = acos(cos_ang)/twopi
c
	call dotvect(n_abc,n_acd,cos_ang)
	ang2 = acos(-cos_ang)/twopi
c
	call dotvect(n_abc,n_bcd,cos_ang)
	ang3 = acos(cos_ang)/twopi
c
	call dotvect(n_abd,n_acd,cos_ang)
	ang4 = acos(cos_ang)/twopi
c
	call dotvect(n_abd,n_bcd,cos_ang)
	ang5 = acos(-cos_ang)/twopi
c
	call dotvect(n_acd,n_bcd,cos_ang)
	ang6 = acos(cos_ang)/twopi
c
	return
	end
c
C	tetra3_noder.f
c
c	Copyright (C) 2002 Patrice Koehl
c
c	This subroutine computes the volume of a tetrahedron, its
c	six dihedral angles, as well as the derivatives of the six
c	dihedral angles with respect to 3 distances (the 3 other
c	distances are considered fixed)
c
	subroutine tetra3_noder(a,b,c,p,rab,rac,rbc,ra,rb,rc,
     1			ang1,ang2,ang3,ang4,ang5,ang6,cosine,sine)
c
	real*8	tetra_volume
c
	real*8	pi,twopi
	real*8	ra,rb,rc
	real*8	rab,rac,rbc
	real*8	ang1,ang2,ang3,ang4,ang5,ang6
	real*8	cos_ang1,cos_ang2,cos_ang3
	real*8	cos_ang4,cos_ang5,cos_ang6
	real*8	sin_ang1,sin_ang2,sin_ang3
	real*8	vol
	real*8	val1,val2
	real*8	sum
	real*8	s1_2,s1
	real*8	s2_2,s2
	real*8	s3_2,s3
	real*8	s4_2,s4,sum_s_2
	real*8	c1,c2,c3,c4
	real*8	a(3),b(3),c(3),p(3)
	real*8	cosine(3),sine(3)
c
	common /constants/pi,twopi
c
	save

c
c	Volume of the tetrahedron, based on all edge lengths
c	(in fact, this is the square of the volume, multiplied
c	by 288)
c	The derivative are really (derivative/vol) and are
c	computed from the fact that
c
c	vol2 = 2*ra2*(4*rb2*rc2 - val_bc*val_bc)
c     1		+ val_ab*(-2*rc2*val_ab - val_bc*val_ac)
c     2		- val_ac*(val_ab*val_bc + 2*rb2*val_ac)
c
c	where vol2 = 288*vol^2
c
	vol = tetra_volume(a,b,c,p)
c
c	Surfaces s1,s2,s3,s4 of the four faces of the tetrahedron,
c
c	(We use the fact that for a triangle T with side lengths
c	a,b,c, then
c
c	P (=perimeter) = (a+b+c)/2
c	Surf**2 = p*(p-a)*(p-b)*(p-c)
c
c	The four triangles considered are:
c
c		Triangle	Surface
c
c		T1 : BCP	s1
c		T2 : ACP	s2
c		T3 : ABP	s3
c		T4 : ABC	s4
c
	sum = (rb + rc + rbc)/2
	s1_2 = sum*(sum-rb)*(sum-rc)*(sum-rbc)
	s1   = sqrt(s1_2)
c
	sum = (ra + rc + rac)/2
	s2_2 = sum*(sum-ra)*(sum-rc)*(sum-rac)
	s2   = sqrt(s2_2)
c
	sum = (ra + rb + rab)/2
	s3_2 = sum*(sum-ra)*(sum-rb)*(sum-rab)
	s3   = sqrt(s3_2)
c
	sum = (rab + rac + rbc)/2
	s4_2 = sum*(sum-rab)*(sum-rac)*(sum-rbc)
	s4 = sqrt(s4_2)
c
	sum_s_2 = s1_2 + s2_2 + s3_2 + s4_2
c
c	Get all six dihedral angles
c
c	ang1 = angle_dihed(a,b,c,p) = angle_dihed(T3,T4)
c	ang2 = angle_dihed(a,c,b,p) = angle_dihed(T2,T4)
c	ang3 = angle_dihed(b,c,a,p) = angle_dihed(T1,T4)
c	ang4 = angle_dihed(a,p,b,c) = angle_dihed(T2,T3)
c	ang5 = angle_dihed(b,p,a,c) = angle_dihed(T1,T3)
c	ang6 = angle_dihed(c,p,a,b) = angle_dihed(T1,T2)
c
	call angle_dihed(b,p,a,c,ang5,cos_ang5)
	call angle_dihed(c,p,a,b,ang6,cos_ang6)
c
c	Get 4 other dihedral angle using cosine rule
c
	val1 = 2*s1*s3*cos_ang5
	val2 = 2*s1*s2*cos_ang6
c
	c1 = sum_s_2 - 2*s1_2
	c2 = sum_s_2 - 2*s2_2 - val1
	c3 = sum_s_2 - 2*s3_2 - val2
	c4 = sum_s_2 - 2*s4_2 - val1 - val2
c
	cos_ang1 = (c1+c2-c3-c4)/(4*s3*s4)
	cos_ang2 = (c1-c2+c3-c4)/(4*s2*s4)
	cos_ang3 = (-c1+c2+c3+c4)/(4*s1*s4)
	cos_ang4 = c4/(2*s2*s3)
c
	ang1 = acos(cos_ang1)/twopi
	ang2 = acos(cos_ang2)/twopi
	ang3 = acos(cos_ang3)/twopi
	ang4 = acos(cos_ang4)/twopi
c
	sin_ang1 = 1.5d0*rab*vol/(s3*s4)
	sin_ang2 = 1.5d0*rac*vol/(s2*s4)
	sin_ang3 = 1.5d0*rbc*vol/(s1*s4)
c	sin_ang4 = 1.5d0*ra*vol/(s2*s3)
c	sin_ang5 = 1.5d0*rb*vol/(s1*s3)
c	sin_ang6 = 1.5d0*rc*vol/(s1*s2)
c
	cosine(1) = cos_ang1
	cosine(2) = cos_ang2
	cosine(3) = cos_ang3
	sine(1) = sin_ang1
	sine(2) = sin_ang2
	sine(3) = sin_ang3
c
	return
	end
c
c	Angle_dihed.f
c
c	Copyright (C) 2002 Patrice Koehl
c
	subroutine angle_dihed(a,b,c,d,ang,cos_ang)
c
	real*8	ang,cos_ang
	real*8	pi,twopi
c
	real*8  a(3),b(3),c(3),d(3),u1(3),u2(3),m(3),n1(3),n2(3)
c
	common /constants/ pi,twopi
c
	save
c

	call diffvect(c,a,u1)
	call diffvect(c,b,u2)
c
	call crossvect(u1,u2,m)
	call unitvector(m,n1)
c
	call diffvect(d,a,u1)
	call diffvect(d,b,u2)
c
	call crossvect(u1,u2,m)
	call unitvector(m,n2)
c
	call dotvect(n1,n2,cos_ang)
c
	ang = acos(cos_ang)/twopi
c
	return
	end
