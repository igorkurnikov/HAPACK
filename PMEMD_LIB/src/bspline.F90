#include "copyright.i"

!*******************************************************************************
!
! Module: bspline_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module bspline_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:  cubspl
!
! Description:
!              
!     ...this code is from netlib...
!
!     from  * a practical guide to splines *  by c. de boor    
!
!     ************************  input  ***************************
!     n = number of data points. assumed to be .ge. 2.
!     (tau(i), c(1,i), i=1,...,n) = abscissae and ordinates of the
!        data points. tau is assumed to be strictly increasing.
!     ibcbeg, ibcend = boundary condition indicators, and
!     c(2,1), c(2,n) = boundary condition information. specifically,
!        ibcbeg = 0  means no boundary condition at tau(1) is given.
!           in this case, the not-a-knot condition is used, i.e. the
!           jump in the third derivative across tau(2) is forced to
!           zero, thus the first and the second cubic polynomial pieces
!           are made to coincide.)
!        ibcbeg = 1  means that the slope at tau(1) is made to equal
!           c(2,1), supplied by input.
!        ibcbeg = 2  means that the second derivative at tau(1) is
!           made to equal c(2,1), supplied by input.
!        ibcend = 0, 1, or 2 has analogous meaning concerning the
!           boundary condition at tau(n), with the additional infor-
!           mation taken from c(2,n).
!
!     ***********************  output  **************************
!     c(j,i), j=1,...,4; i=1,...,l (= n-1) = the polynomial coefficients
!        of the cubic interpolating spline with interior knots (or
!        joints) tau(2), ..., tau(n-1). precisely, in the interval
!        (tau(i), tau(i+1)), the spline f is given by
!           f(x) = c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)/3.)/2.)
!        where h = x - tau(i). the function program *ppvalu* may be
!        used to evaluate f or its derivatives from tau,c, l = n-1,
!        and k=4.
!
!*******************************************************************************

subroutine cubspl(tau, c, n, ibcbeg, ibcend)

  implicit none

  integer               n
  double precision      tau(n)
  double precision      c(4, n)
  integer               ibcbeg
  integer               ibcend

  double precision      divdf1, divdf3, dtau, g
  integer               i, j, l, m

! A tridiagonal linear system for the unknown slopes s(i) of f 
! at tau(i), i=1,...,n, is generated and then solved by gauss elimination,
! with s(i) ending up in c(2,i), all i.
! c(3,.) and c(4,.) are used initially for temporary storage.

    l = n - 1

! Compute first differences of tau sequence and store in c(3,.).
! Also, compute first divided difference of data and store in c(4,.).

    do m = 2, n
      c(3, m) = tau(m) - tau(m - 1)
      c(4, m) = (c(1, m) - c(1, m - 1))/c(3, m)
    end do

! Construct first equation from the boundary condition, of the form
! c(4,1) * s(1) + c(3,1) * s(2) = c(2,1)

    if (ibcbeg - 1) 11, 15, 16

 11 if (n .gt. 2) goto 12

! No condition at left end and n = 2.

    c(4, 1) = 1.d0
    c(3, 1) = 1.d0
    c(2, 1) = 2.d0 * c(4, 2)

    goto 25

! Not-a-knot condition at left end and n .gt. 2.

 12 c(4, 1) = c(3, 3)
    c(3, 1) = c(3, 2) + c(3, 3)
    c(2, 1) = ((c(3, 2) + 2.d0 * c(3, 1)) * c(4, 2) * c(3, 3) + &
    c(3, 2) * c(3, 2) * c(4, 3))/c(3, 1)

    goto 19

! Slope prescribed at left end.

 15 c(4, 1) = 1.d0
    c(3, 1) = 0.d0

    goto 18

! Second derivative prescribed at left end.

 16 c(4, 1) = 2.d0
    c(3, 1) = 1.d0
    c(2, 1) = 3.d0 * c(4, 2) - c(3, 2)/2.d0 * c(2, 1)

 18 if (n .eq. 2) goto 25

! If there are interior knots, generate the corresp. equations and carry
! out the forward pass of gauss elimination, after which the m-th
! equation reads    c(4,m) * s(m) + c(3,m) * s(m+1) = c(2,m).

 19 do m = 2, l
      g = - c(3, m + 1) / c(4, m - 1)

      c(2, m) = g * c(2, m - 1) + 3.d0 * (c(3, m) * c(4, m + 1) + &
                c(3, m + 1) * c(4, m))

      c(4, m) = g * c(3, m - 1) + 2.d0 * (c(3, m) + c(3, m + 1))
    end do

! Construct last equation from the second boundary condition, of the form
! (-g * c(4,n-1)) * s(n-1) + c(4,n) * s(n) = c(2,n).
! If slope is prescribed at right end, one can go directly to back-
! substitution, since c array happens to be set up just right for it
! at this point.

    if (ibcend - 1) 21, 30, 24

 21 if (n .eq. 3 .and. ibcbeg .eq. 0) goto 22

! Not-a-knot and n .ge. 3, and either n .gt. 3 or also not-a-knot at
! left end point.

    g = c(3, n - 1) + c(3, n)

    c(2, n) = ((c(3, n) + 2.d0 * g) * c(4, n) * c(3, n - 1) + &
                c(3, n) * c(3, n) * (c(1, n - 1) - c(1, n - 2))/c(3, n - 1))/g

    g = - g/c(4, n - 1)
    c(4, n) = c(3, n - 1)

    goto 29

! Either (n=3 and not-a-knot also at left) or (n=2 and not not-a-
! knot at left end point).

 22 c(2, n) = 2.d0 * c(4, n)
    c(4, n) = 1.d0

    goto 28

! Second derivative prescribed at right endpoint.

 24 c(2, n) = 3.d0 * c(4, n) + c(3, n)/2.d0 * c(2, n)
    c(4, n) = 2.d0

    goto 28

 25 if (ibcend-1) 26, 30, 24

 26 if (ibcbeg .gt. 0) goto 22

! Not-a-knot at right endpoint and at left endpoint and n = 2.

    c(2, n) = c(4, n)

    goto 30

 28 g = - 1.d0 / c(4, n - 1)

! Complete forward pass of gauss elimination.

 29 c(4, n) = g * c(3, n - 1) + c(4, n)
    c(2, n) = (g * c(2, n - 1) + c(2, n))/c(4, n)

! carry out back substitution

 30 j = l 

 40 c(2, j) = (c(2, j) - c(3, j) * c(2, j + 1))/c(4, j)
    j = j - 1
    if (j .gt. 0) goto 40

! Generate cubic coefficients in each interval, i.e., the deriv.s
! at its left endpoint, from value and slope at its endpoints.

    do i = 2, n
      dtau = c(3, i)
      divdf1 = (c(1, i) - c(1, i - 1))/dtau
      divdf3 = c(2, i - 1) + c(2, i) - 2.d0 * divdf1
      c(3, i - 1) = 2.d0 * (divdf1 - c(2, i - 1) - divdf3)/dtau
      c(4, i - 1) = (divdf3/dtau) * (6.d0/dtau)
    end do

    return

end subroutine cubspl

!*******************************************************************************
!
! Subroutine:  fill_bspline_0
!
! Description:
!
! Use standard B-spline recursions: see doc file.
! w is fraction between 0 and 1; 
! order is the order of interpolation;
! array is the array of shifted bsplines;
! Using notation from Essmann et al; w = u-[u] and
! array(j) = M_n(w + order - j)  where n is order & w = u - [u]; 
!              
!*******************************************************************************

subroutine fill_bspline_0(w, order, array)

  implicit none

  double precision      w
  integer               order
  double precision      array(order)

  double precision      div
  integer               k

! Init order 2:

  array(2) = w
  array(1) = 1.d0 - w

  if (order .eq. 2) return

! One pass to order 3:

  array(3) = 0.5d0 * w * array(2)
  array(2) = 0.5d0 * ((w + 1.d0) * array(1) + (2.d0 - w) * array(2))
  array(1) = 0.5d0 * (1.d0 - w) * array(1)

  if (order .eq. 3) return

! One pass to order 4:

  div = 1.d0/3.d0
  array(4) = div * w * array(3)
  array(3) = div * ((w + 1.d0) * array(2) + (3.d0 - w) * array(3))
  array(2) = div * ((w + 2.d0) * array(1) + (2.d0 - w) * array(2))
  array(1) = div * (1.d0 - w) * array(1)

  if (order .eq. 4) return

  do k = 5, order
    call one_pass_bspline(array, w, k)
  end do

  return

end subroutine fill_bspline_0

!*******************************************************************************
!
! Subroutine:  fill_bspline_1
!
! Description:
!              
! Use standard B-spline recursions: see doc file.
! w is fraction between 0 and 1; 
! order is the order of interpolation;
! array is the array of shifted bsplines;
! darray is array of derivs;
! Using notation from Essmann et al; w = u-[u] and
! array(j) = M_n(w + order - j)  where n is order & w = u - [u]; 
!
!*******************************************************************************

subroutine fill_bspline_1(w, order, array, darray)

  implicit none

! Formal arguments:

  double precision      :: w
  integer               :: order
  double precision      :: array(order)
  double precision      :: darray(order)

! Local variables:

  integer                       :: k
  double precision, parameter   :: one_third = 1.d0 / 3.d0

! This routine works for orders of 3 or more!  There is no error checking here
! for order 2; it should be in input checking instead.

! Init order 2:

  array(1) = 1.d0 - w
  array(2) = w

  if (order .eq. 4) then

! One pass to order 3:
    
    array(3) = 0.5d0 * w * array(2)
    array(2) = 0.5d0 * ((w + 1.d0) * array(1) + (2.d0 - w) * array(2))
    array(1) = 0.5d0 * (1.d0 - w) * array(1)

! Diff to get darray:

    darray(1) = -array(1)
    darray(2) = array(1) - array(2)
    darray(3) = array(2) - array(3)
    darray(4) = array(3)

! One final pass to order 4:

    array(4) = one_third * w * array(3)
    array(3) = one_third * ((w + 1.d0) * array(2) + (3.d0 - w) * array(3))
    array(2) = one_third * ((w + 2.d0) * array(1) + (2.d0 - w) * array(2))
    array(1) = one_third * (1.d0 - w) * array(1)

  else if (order .gt. 4) then

! General order case:
! One pass to order 3:

    array(3) = 0.5d0 * w * array(2)
    array(2) = 0.5d0 * ((w + 1.d0) * array(1) + (2.d0 - w) * array(2))
    array(1) = 0.5d0 * (1.d0 - w) * array(1)

! Another pass to order 4:

    array(4) = one_third * w * array(3)
    array(3) = one_third * ((w + 1.d0) * array(2) + (3.d0 - w) * array(3))
    array(2) = one_third * ((w + 2.d0) * array(1) + (2.d0 - w) * array(2))
    array(1) = one_third * (1.d0 - w) * array(1)

! Compute standard b-spline recursion:

    do k = 5, order - 1
      call one_pass_bspline(array, w, k)
    end do

    call diff_bspline(array, darray, order)

! One more recursion:

    call one_pass_bspline(array, w, order)

  else ! order had better be 3...

! Deriv:

    darray(1) = -array(1)
    darray(2) = array(1) - array(2)
    darray(3) = array(2)

! One pass to order 3:

    array(3) = 0.5d0 * w * array(2)
    array(2) = 0.5d0 * ((w + 1.d0) * array(1) + (2.d0 - w) * array(2))
    array(1) = 0.5d0 * (1.d0 - w) * array(1)

  end if

  return

end subroutine fill_bspline_1


!*******************************************************************************
!
! Subroutine:  fill_bspline_1_3d
!
! Description:
!              
! Use standard B-spline recursions: see doc file.
! w is fraction between 0 and 1; 
! order is the order of interpolation;
! array is the 3d array of shifted bsplines;
! darray is array of derivs;
! Using notation from Essmann et al; w = u-[u] and
! array(j) = M_n(w + order - j)  where n is order & w = u - [u]; 
!
!*******************************************************************************

subroutine fill_bspline_1_3d(w, order, array, darray)

  implicit none

! Formal arguments:

  double precision      :: w(3)
  integer               :: order
  double precision      :: array(order, 3)
  double precision      :: darray(order, 3)

! Local variables:

  integer                       :: i, j
  double precision, parameter   :: one_third = 1.d0 / 3.d0

! This routine works for orders of 3 or more!  There is no error checking here
! for order 2; it should be in input checking instead.

! Init order 2:

  array(1, :) = 1.d0 - w(:)
  array(2, :) = w(:)

  if (order .eq. 4) then

    do i = 1, 3
! One pass to order 3:
    
      array(3, i) = 0.5d0 * w(i) * array(2, i)
      array(2, i) = 0.5d0 * ((w(i) + 1.d0) * array(1, i) + &
                    (2.d0 - w(i)) * array(2, i))
      array(1, i) = 0.5d0 * (1.d0 - w(i)) * array(1, i)

! Diff to get darray:

      darray(1, i) = -array(1, i)
      darray(2, i) = array(1, i) - array(2, i)
      darray(3, i) = array(2, i) - array(3, i)
      darray(4, i) = array(3, i)

! One final pass to order 4:

      array(4, i) = one_third * w(i) * array(3, i)
      array(3, i) = one_third * ((w(i) + 1.d0) * array(2, i) + &
                    (3.d0 - w(i)) * array(3, i))
      array(2, i) = one_third * ((w(i) + 2.d0) * array(1, i) + &
                    (2.d0 - w(i)) * array(2, i))
      array(1, i) = one_third * (1.d0 - w(i)) * array(1, i)
    end do

  else if (order .gt. 4) then

    do i = 1, 3
! General order case:
! One pass to order 3:

      array(3, i) = 0.5d0 * w(i) * array(2, i)
      array(2, i) = 0.5d0 * ((w(i) + 1.d0) * array(1, i) + &
                    (2.d0 - w(i)) * array(2, i))
      array(1, i) = 0.5d0 * (1.d0 - w(i)) * array(1, i)

! Another pass to order 4:

      array(4, i) = one_third * w(i) * array(3, i)
      array(3, i) = one_third * ((w(i) + 1.d0) * array(2, i) + &
                    (3.d0 - w(i)) * array(3, i))
      array(2, i) = one_third * ((w(i) + 2.d0) * array(1, i) + &
                    (2.d0 - w(i)) * array(2, i))
      array(1, i) = one_third * (1.d0 - w(i)) * array(1, i)

! Compute standard b-spline recursion:

      do j = 5, order - 1
        call one_pass_bspline(array(1, i), w(i), j)
      end do

      call diff_bspline(array(1, i), darray(1, i), order)

! One more recursion:

      call one_pass_bspline(array(1, i), w(i), order)

    end do

  else ! order had better be 3...

! Deriv:

    do i = 1, 3

      darray(1, i) = -array(1, i)
      darray(2, i) = array(1, i) - array(2, i)
      darray(3, i) = array(2, i)

! One pass to order 3:

      array(3, i) = 0.5d0 * w(i) * array(2, i)
      array(2, i) = 0.5d0 * ((w(i) + 1.d0) * array(1, i) + &
                    (2.d0 - w(i)) * array(2, i))
      array(1, i) = 0.5d0 * (1.d0 - w(i)) * array(1, i)

    end do

  end if

  return

end subroutine fill_bspline_1_3d

!*******************************************************************************
!
! Subroutine:  fill_bspline_2
!
! Description:
!
! Use standard B-spline recursions: see doc file.
! w is fraction between 0 and 1; 
! order is the order of interpolation;
! array is the array of shifted bsplines;
! darray is array of derivs;
! d2array is array of 2nd derivs;
! Using notation from Essmann et al; w = u-[u] and
! array(j) = M_n(w + order - j)  where n is order & w = u - [u]; 
!              
!*******************************************************************************

subroutine fill_bspline_2(w, order, array, darray, d2array)

  implicit none

  double precision      w
  integer               order
  double precision      array(order)
  double precision      darray(order)
  double precision      d2array(order)

  double precision      div
  integer               k

  if (order .eq. 4) then

! Init order 2:

    array(2) = w
    array(1) = 1.d0 - w

! In full:
!   order three deriv
!      darray(1) = -array(1)
!      darray(2) = array(1) - array(2)
!      darray(3) = array(2)
!   deriv of deriv
!      d2array(1) = -darray(1)
!      d2array(2) = darray(1) - darray(2)
!      d2array(3) = darray(2) - darray(3)
!      d2array(4) = darray(3)
! In short:

    d2array(1) = array(1)
    d2array(2) = array(2) - 2.d0 * array(1)
    d2array(3) = array(1) - 2.d0 * array(2)
    d2array(4) = array(2)

! One pass to order 3:

    array(3) = 0.5d0 * w * array(2)
    array(2) = 0.5d0 * ((w + 1.d0) * array(1) + (2.d0 - w) * array(2))
    array(1) = 0.5d0 * (1.d0 - w) * array(1)

! Order four deriv:

    darray(1) = - array(1)
    darray(2) = array(1) - array(2)
    darray(3) = array(2) - array(3)
    darray(4) = array(3)

! One pass to order 4:

    div = 1.d0/3.d0
    array(4) = div * w * array(3)
    array(3) = div * ((w + 1.d0) * array(2) + (3.d0 - w) * array(3))
    array(2) = div * ((w + 2.d0) * array(1) + (2.d0 - w) * array(2))
    array(1) = div * (1.d0 - w) * array(1)

    return

  else

! Init order 2:

    array(2) = w
    array(1) = 1.d0 - w

! One pass to order 3:

    array(3) = 0.5d0 * w * array(2)
    array(2) = 0.5d0 * ((w + 1.d0) * array(1) + (2.d0 - w) * array(2))
    array(1) = 0.5d0 * (1.d0 - w) * array(1)

! Compute standard b-spline recursion:

    do k = 4, order - 2
      call one_pass_bspline(array, w, k)
    end do

! Perform standard b-spline differentiation:

    call diff_bspline(array, darray, order - 1) 
! Deriv of deriv:

    call diff_bspline(darray, d2array, order)

! Perform standard b-spline differentiation:

    call one_pass_bspline(array, w, order - 1)
    call diff_bspline(array, darray, order)

! One more recursion:

    call one_pass_bspline(array, w, order)

    return

  end if

end subroutine fill_bspline_2

!*******************************************************************************
!
! Subroutine:  diff_bspline
!
! Description: 
!              
! Using notation from Essmann et al; w = u-[u] and
! array(j) = M_n(w + order - j)  where n is order
! DERIVATIVE:    d/dw M_n(w) = M_n-1(w) - M_n-1(w-1)
! i.e.   d/dw M_n(w+n-j) = M_n-1(w+n-j) - M_n-1(w+n-j-1)
! i.e.   new(j) = old(j-1) - old(j)
! where old is array before one_pass (thus n->n-1) and new is array afterwards.
!
!*******************************************************************************

subroutine diff_bspline(c, d, n)

  implicit none

  double precision c(*), d(*)
  integer n

  integer j

  d(1) = -c(1)

  do j = 2, n - 1
    d(j) = c(j - 1) - c(j)
  end do

  d(n) = c(n - 1)

  return

end subroutine diff_bspline

!*******************************************************************************
!
! Subroutine:  one_pass_bspline
!
! Description:
!              
! Using notation from Essmann et al; w = u-[u] and
! array(j) = M_n(w + order - j)  where n is order
! RECURSION:  M_n(w) = (w/(n-1))*M_n-1(w)+((n-w)/(n-1))*M_n-1(w-1)
! i.e.   M_n(w+n-j) = ((w+n-j)/(n-1))*M_n-1(w+n-j)+((j-w)/(n-1))*M_n-1(w+n-j-1)
! i.e.   new(j) = ((w+n-j)/(n-1))*old(j-1) + ((j-w)/(n-1))*old(j)
! where old is array before one_pass (thus n->n-1) and new is array afterwards.
! Write backwards to do it with one array.

!*******************************************************************************

subroutine one_pass_bspline(c, w, n)

  implicit none

  double precision      c(*), w
  integer               n

  double precision      div
  integer               j

  div = 1.d0 / (n - 1)
  c(n) = div * w * c(n - 1)

  do j = 1, n - 2
    c(n - j) = div * ((w + j) * c(n - j - 1) + (n - j - w) * c(n - j))
  end do

  c(1) = div * (1 - w) * c(1)

  return

end subroutine one_pass_bspline

end module bspline_mod
