#include "copyright.i"

!*******************************************************************************!
! Function:  am_val_real_array_index
!
! Description: <TBS>
!
!*******************************************************************************

function am_val_real_array_index(val, array, num)

  implicit none

! Formal arguments:

  integer, intent(in)           :: num
  double precision, intent(in)  :: val, array(num)

! Local variables:

  integer                       :: am_val_real_array_index
  integer                       :: indhi, indlo, ind

! Bisection search for ind just before val (i.e. array(ind)<val<array(ind + 1)).
! Assume array is ordered in increasing order.

  indlo = 1 
  indhi = num

  do while (indhi - indlo .gt. 1)
    ind = (indhi + indlo) / 2
    if (array(ind) .gt. val) then
      indhi = ind
    else
      indlo = ind
    end if
  end do

  am_val_real_array_index = indlo

  return

end function am_val_real_array_index

!*******************************************************************************!
! Subroutine:  am_val_ftab_eval_f_df
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_val_ftab_eval_f_df(nlist, degree, coeff, arg, func, dfunc_darg)

  implicit none

! Formal arguments:

  integer, intent(in)                   :: nlist
  integer, intent(in)                   :: degree
  double precision, intent(in)          :: coeff(0:degree)
  double precision, intent(in)          :: arg(nlist)
  double precision, intent(out)         :: func(nlist)
  double precision, intent(out)         :: dfunc_darg(nlist)

! Local variables:

  integer                               :: n, deg
  double precision                      :: dx

! Set up in case of smaller degrees:

  if (degree .eq. 2) then
    do n = 1, nlist
      dx = arg(n)
      func(n) = coeff(0) + dx * (coeff(1) + dx * coeff(2))
      dfunc_darg(n) = coeff(1) + 2.d0 * dx * coeff(2)
    end do 
  else if (degree .eq. 3) then
    do n = 1, nlist
      dx = arg(n)
      func(n) = coeff(0) + dx * (coeff(1) + dx * (coeff(2) + dx * coeff(3)))
      dfunc_darg(n) = coeff(1) + dx * (2.d0 * coeff(2) + dx * 3.d0 * coeff(3))
    end do 
  else if (degree .eq. 4) then
    do n = 1, nlist

      dx = arg(n)

      func(n) = coeff(0) + dx * (coeff(1) + &
                           dx * (coeff(2) + &
                           dx * (coeff(3) + &
                           dx* coeff(4))))

      dfunc_darg(n) = coeff(1) + dx * (2.d0 * coeff(2) + dx * &
                      (3.d0 * coeff(3) + dx * 4.d0 * coeff(4)))

    end do 
  else 
    do n = 1, nlist
      dx = arg(n)
      func(n) = coeff(degree)
      dfunc_darg(n) = 0.d0
      do deg = degree - 1, 0, -1
        dfunc_darg(n) = func(n) + dx * dfunc_darg(n)
        func(n) = coeff(deg) + dx * func(n)
      end do
    end do
  end if

  return

end subroutine am_val_ftab_eval_f_df

!*******************************************************************************!
! Subroutine:  am_val_geom_torsion
!
! Description:
!
! Given coords of points a, b, c, d this routine calculates cosine and sine of
! torsion phi as well as gradient of phi with respect to coords of a, b, c, d.
! Units are in radians.
!
!*******************************************************************************

subroutine am_val_geom_torsion(crd_abcd, gradphi_abcd, cosphi, sinphi)

  implicit none

! Formal arguments:

  double precision, intent(in)  :: crd_abcd(12)
  double precision, intent(out) :: gradphi_abcd(12), cosphi, sinphi

! Local variables:

  double precision              :: rab(3)
  double precision              :: rcb(3)
  double precision              :: rdc(3)
  double precision              :: ucb(3)
  double precision              :: rcross(3)
  double precision              :: upab(3)
  double precision              :: upabc(3)
  double precision              :: upbcd(3)
  double precision              :: updc(3)
  double precision              :: S(3)
  double precision              :: siz, sizcb, sizpab, sizpdc
  double precision              :: dotp_ab_cb, dotp_dc_cb
  double precision              :: dot
  integer                       :: m

  do m = 1, 3
    rab(m) = crd_abcd(m) - crd_abcd(m + 3)
    rcb(m) = crd_abcd(m + 6) - crd_abcd(m + 3)
    rdc(m) = crd_abcd(m + 9) - crd_abcd(m + 6)
  end do

  sizcb = sqrt(rcb(1) * rcb(1) + rcb(2) * rcb(2) + rcb(3) * rcb(3))

  ucb(1) = rcb(1) / sizcb
  ucb(2) = rcb(2) / sizcb
  ucb(3) = rcb(3) / sizcb

  dotp_ab_cb = rab(1) * ucb(1) + rab(2) * ucb(2) + rab(3) * ucb(3)

! upab is unit vector along component rab perp to ucb

  dot = rab(1) * ucb(1) + rab(2) * ucb(2) + rab(3) * ucb(3)
  
  upab(1) = rab(1) - dot * ucb(1)
  upab(2) = rab(2) - dot * ucb(2)
  upab(3) = rab(3) - dot * ucb(3)
  
  sizpab = sqrt(upab(1) * upab(1) + upab(2) * upab(2) + upab(3) * upab(3))
  
  upab(1) = upab(1) / sizpab
  upab(2) = upab(2) / sizpab
  upab(3) = upab(3) / sizpab
  
  dotp_dc_cb = rdc(1) * ucb(1) + rdc(2) * ucb(2) + rdc(3) * ucb(3)

! updc is unit vector along component rdc perp to ucb

  dot = rdc(1) * ucb(1) + rdc(2) * ucb(2) + rdc(3) * ucb(3)

  updc(1) = rdc(1) - dot * ucb(1)
  updc(2) = rdc(2) - dot * ucb(2)
  updc(3) = rdc(3) - dot * ucb(3)

  sizpdc = sqrt(updc(1) * updc(1) + updc(2) * updc(2) + updc(3) * updc(3))

  updc(1) = updc(1) / sizpdc
  updc(2) = updc(2) / sizpdc
  updc(3) = updc(3) / sizpdc

! cosine of phi is given by dot product of upab and updc

  cosphi = upab(1) * updc(1) + upab(2) * updc(2) + upab(3) * updc(3)

! sine of phi is given by dot product of ucb and upab x updc

  rcross(1) = upab(2) * updc(3) - upab(3) * updc(2)
  rcross(2) = upab(3) * updc(1) - upab(1) * updc(3)
  rcross(3) = upab(1) * updc(2) - upab(2) * updc(1)
  sinphi = rcross(1) * ucb(1) + rcross(2) * ucb(2) + rcross(3) * ucb(3)
  
! gradient of phi wrt ra is perp to abc plane---movement of ra by dr perp
! to abc plane results in dphi of dr / sizpab
! perp to abc given by upab x ucb  (these are orthogonal unit vectors)

  upabc(1) = upab(2) * ucb(3) - upab(3) * ucb(2)
  upabc(2) = upab(3) * ucb(1) - upab(1) * ucb(3)
  upabc(3) = upab(1) * ucb(2) - upab(2) * ucb(1)

! grad of phi wrt rd is perp to bcd plane--calc sim to grad phi wrt ra
! perp given by updc x ucb or ucb x updc

  upbcd(1) = ucb(2) * updc(3) - ucb(3) * updc(2)
  upbcd(2) = ucb(3) * updc(1) - ucb(1) * updc(3)
  upbcd(3) = ucb(1) * updc(2) - ucb(2) * updc(1)

! now have enough for gradphi for a and d

  do m = 1, 3
    gradphi_abcd(m) = upabc(m) / sizpab
    gradphi_abcd(9 + m) = upbcd(m) / sizpdc
  end do

! following chap 5 of thesis of Bekker we have grad phi wrt b = -grad phi wrt a
! plus some vec S and rad phi wrt c = -grad phi wrt d - S
! S is perp to rcb; using simple torque rule and identity for 
! triple cross product he derives S (eqn 5.20)

  do m = 1, 3
    S(m) = (dotp_ab_cb / sizcb) * gradphi_abcd(m) + &
           (dotp_dc_cb / sizcb) * gradphi_abcd(m + 9)
    gradphi_abcd(m + 3) = S(m) - gradphi_abcd(m)
    gradphi_abcd(m + 6) = -S(m) - gradphi_abcd(m + 9)
  end do

  return

end subroutine am_val_geom_torsion

!*******************************************************************************!
! Function:  vec3d_unitperpto_unitvec
!
! Description:
!
! Removes component of v along unit vector u and returns length of new w;
! normalizes resulting w.
!
!*******************************************************************************

function vec3d_unitperpto_unitvec(v, u, w)

implicit none

! Formal arguments:

  double precision      :: v(3), u(3), w(3)

! Local variables:

  double precision      :: vec3d_unitperpto_unitvec
  double precision      :: siz, dot

  dot = v(1) * u(1) + v(2) * u(2) + v(3) * u(3)

  w(1) = v(1) - dot * u(1)
  w(2) = v(2) - dot * u(2)
  w(3) = v(3) - dot * u(3)

  siz = sqrt(w(1) * w(1) + w(2) * w(2) + w(3) * w(3))

  w(1) = w(1) / siz
  w(2) = w(2) / siz
  w(3) = w(3) / siz

  vec3d_unitperpto_unitvec = siz

  return

end function vec3d_unitperpto_unitvec

!*******************************************************************************!
! Subroutine:  am_val_vec3d_get_perp_to_vecs
!
! Description: 
!
! Given two vectors v1, v2, with associated unit vectors u1, u2, get the unit 
! cross product p = v1 x v2 / |v1 x v2|.
! A change dv1 in v1 within v1v2 plane has no effect on p;
! only a change dv1_p along p can change p --- to produce this change without
! changing v2 need to rotate about v2; this causes change in p along u2p
! where u2p is u2 x p; change dv1_p is equiv to rotation of size
! dv1_p / v1_dot_u2p about axis u2 = v2 / |v2| , where v1_dot_u2p is dot 
! product of v1 with u2p; change in p is therefore
! dp = (-dv1_p / v1_dot_u2p) * u2p thus
! dp_dv1_p = -u2p / v1_dot_u2p
!
!*******************************************************************************

subroutine am_val_vec3d_get_perp_to_vecs(v1, v2, p, dp_dv1_p, dp_dv2_p)

  implicit none

! Formal arguments:

  double precision :: v1(3), v2(3), p(3), dp_dv1_p(3), dp_dv2_p(3)

! Local variables:

  double precision :: siz, u1(3), u2(3), u1p(3), u2p(3), dotp
  
  siz = v1(1) * v1(1) + v1(2) * v1(2) + v1(3) * v1(3)

  u1(1) = v1(1) / siz
  u1(2) = v1(2) / siz
  u1(3) = v1(3) / siz

  siz = v2(1) * v2(1) + v2(2) * v2(2) + v2(3) * v2(3)

  u2(1) = v2(1) / siz
  u2(2) = v2(2) / siz
  u2(3) = v2(3) / siz

! p = u1 x u2 normed since u1 not perp to u2 necessarily

  p(1) = u1(2) * u2(3) - u1(3) * u2(2)
  p(2) = u1(3) * u2(1) - u1(1) * u2(3)
  p(3) = u1(1) * u2(2) - u1(2) * u2(1)

  siz = sqrt(p(1) * p(1) + p(2) * p(2) + p(3) * p(3))

  p(1) = p(1) / siz
  p(2) = p(2) / siz
  p(3) = p(3) / siz

! u1p is unit vec perp to u1 and p : p x u1

  u1p(1) = p(2) * u1(3) - p(3) * u1(2)
  u1p(2) = p(3) * u1(1) - p(1) * u1(3)
  u1p(3) = p(1) * u1(2) - p(2) * u1(1)

! u2p is unit vec perp to u2 and p : u2 x p

  u2p(1) = u2(2) * p(3) - u2(3) * p(2)
  u2p(2) = u2(3) * p(1) - u2(1) * p(3)
  u2p(3) = u2(1) * p(2) - u2(2) * p(1)

! change in v1 along p gives dp along -u2p

  dotp = v1(1) * u2p(1) + v1(2) * u2p(2) + v1(3) * u2p(3)
  dp_dv1_p(1) = -u2p(1) / dotp
  dp_dv1_p(2) = -u2p(2) / dotp
  dp_dv1_p(3) = -u2p(3) / dotp

! change in v2 along p gives dp along -u1p

  dotp = v2(1) * u1p(1) + v2(2) * u1p(2) + v2(3) * u1p(3)
  dp_dv2_p(1) = -u1p(1) / dotp
  dp_dv2_p(2) = -u1p(2) / dotp
  dp_dv2_p(3) = -u1p(3) / dotp

  return

end subroutine am_val_vec3d_get_perp_to_vecs

!*******************************************************************************!
! Subroutine:  am_val_bcuint1
!
! Description: 
!
! "bcuint1" performs a bicubic interpolation of the function
! value and gradient along the directions of a 2D spline grid
!
! contributed by Jay William Ponder
!
! literature reference:
!
! W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P.
! Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge
! University Press, 1992, Section 3.6
!
!*******************************************************************************

subroutine am_val_bcuint1(y, y1, y2, y12, x1l, x1u, x2l, x2u, x1, x2, &
                          ansy, ansy1, ansy2)

  implicit none

! Formal arguments:

  double precision, intent(in)  ::  y(4), y1(4), y2(4), y12(4)
  double precision, intent(in)  :: x1, x1l, x1u
  double precision, intent(in)  :: x2, x2l, x2u
  double precision, intent(out) :: ansy, ansy1, ansy2

! Local variables:

  integer                       :: i
  double precision              :: t, u, c(4,4)

! get coefficients, then perform bicubic interpolation

  call am_val_bcucof(y, y1, y2, y12, x1u - x1l, x2u - x2l, c)

  t = (x1 - x1l) / (x1u - x1l)
  u = (x2 - x2l) / (x2u - x2l)

  ansy = 0.d0
  ansy1 = 0.d0
  ansy2 = 0.d0

  do i = 4, 1, -1
     ansy = t * ansy + ((c(i,4) * u + c(i,3)) * u + c(i,2)) * u + c(i,1)
     ansy1 = u * ansy1 + (3.0d0 * c(4,i) * t + 2.0d0 * c(3,i)) * t + c(2,i)
     ansy2 = t * ansy2 + (3.0d0 * c(i,4) * u + 2.0d0 * c(i,3)) * u + c(i,2)
  end do

  ansy1 = ansy1 / (x1u - x1l)
  ansy2 = ansy2 / (x2u - x2l)

  return

end subroutine am_val_bcuint1

!*******************************************************************************!
! Subroutine:  am_val_bcucof
!
! Description:
!
!   "bcucof" determines the coefficient matrix needed for bicubic
!   interpolation of a function, gradients and cross derivatives.
!
!*******************************************************************************

subroutine am_val_bcucof(y, y1, y2, y12, d1, d2, c)

  implicit none

! Formal arguments:

  double precision, intent(in)  :: y(4), y1(4), y2(4), y12(4)
  double precision, intent(in)  :: d1, d2
  double precision, intent(out) :: c(4,4)

! Local variables:

  double precision              :: xx, d1d2
  double precision              ::  x(16), cl(16)
  double precision, save        ::  wt(16,16)
  integer                       ::  i, j, k
  
  data wt / 1.0d0, 0.0d0, -3.0d0, 2.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,  &
            -3.0d0, 0.0d0, 9.0d0, -6.0d0, 2.0d0, 0.0d0, -6.0d0, 4.0d0, &
             0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
             3.0d0, 0.0d0, -9.0d0, 6.0d0, -2.0d0, 0.0d0, 6.0d0, -4.0d0, &
             0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
             0.0d0, 0.0d0, 9.0d0, -6.0d0, 0.0d0, 0.0d0, -6.0d0, 4.0d0, &
             0.0d0, 0.0d0, 3.0d0, -2.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
             0.0d0, 0.0d0, -9.0d0, 6.0d0, 0.0d0, 0.0d0, 6.0d0, -4.0d0, &
             0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, -3.0d0, 2.0d0, &
            -2.0d0, 0.0d0, 6.0d0, -4.0d0, 1.0d0, 0.0d0, -3.0d0, 2.0d0, &
             0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
            -1.0d0, 0.0d0, 3.0d0, -2.0d0, 1.0d0, 0.0d0, -3.0d0, 2.0d0, &
             0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
             0.0d0, 0.0d0, -3.0d0, 2.0d0, 0.0d0, 0.0d0, 3.0d0, -2.0d0, &
             0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 3.0d0, -2.0d0, &
             0.0d0, 0.0d0, -6.0d0, 4.0d0, 0.0d0, 0.0d0, 3.0d0, -2.0d0, &
             0.0d0, 1.0d0, -2.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
             0.0d0, -3.0d0, 6.0d0, -3.0d0, 0.0d0, 2.0d0, -4.0d0, 2.0d0, &
             0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
             0.0d0, 3.0d0, -6.0d0, 3.0d0, 0.0d0, -2.0d0, 4.0d0, -2.0d0, &
             0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
             0.0d0, 0.0d0, -3.0d0, 3.0d0, 0.0d0, 0.0d0, 2.0d0, -2.0d0, &
             0.0d0, 0.0d0, -1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
             0.0d0, 0.0d0, 3.0d0, -3.0d0, 0.0d0, 0.0d0, -2.0d0, 2.0d0, &
             0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, -2.0d0, 1.0d0, &
             0.0d0, -2.0d0, 4.0d0, -2.0d0, 0.0d0, 1.0d0, -2.0d0, 1.0d0, &
             0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
             0.0d0, -1.0d0, 2.0d0, -1.0d0, 0.0d0, 1.0d0, -2.0d0, 1.0d0, &
             0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
             0.0d0, 0.0d0, 1.0d0, -1.0d0, 0.0d0, 0.0d0, -1.0d0, 1.0d0, &
             0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, -1.0d0, 1.0d0, &
             0.0d0, 0.0d0, 2.0d0, -2.0d0, 0.0d0, 0.0d0, -1.0d0, 1.0d0 /

! Pack a temporary vector of corner values

  d1d2 = d1 * d2

  do i = 1, 4
     x(i) = y(i)
     x(i + 4) = y1(i) * d1
     x(i + 8) = y2(i) * d2
     x(i + 12) = y12(i) * d1d2
  end do

! Matrix multiply by the stored weight table

  do i = 1, 16
     xx = 0.0d0
     do k = 1, 16
        xx = xx + wt(i,k) * x(k)
     end do
     cl(i) = xx
  end do

! Unpack the result into the coefficient table

  j = 0
  do i = 1, 4
     do k = 1, 4
        j = j + 1
        c(i,k) = cl(j)
     end do
  end do

  return

end subroutine am_val_bcucof

subroutine dummy_amoeba_valence
  return
end subroutine dummy_amoeba_valence
