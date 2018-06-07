#include "copyright.i"

!*******************************************************************************!
! Module: amoeba_multipoles_mod
!
! Description: <TBS>
!
!*******************************************************************************

module amoeba_multipoles_mod

  implicit none

  private

! Data that should be broadcast to slave processes from the master:

  integer, parameter    :: amoeba_multipoles_int_cnt = 4

  integer                               do_amoeba_multipoles_flag, &
                                        num_multipoles, num_chiral_frame_list, &
                                        num_frame_def_list

  common / amoeba_multipoles_int /      do_amoeba_multipoles_flag, &
                                        num_multipoles, num_chiral_frame_list, &
                                        num_frame_def_list

  save  :: / amoeba_multipoles_int /

  double precision, save, allocatable :: local_multipole(:,:)

  type chiral_frame
    integer :: frame_index
    integer :: fourth_atom
    integer :: chirality
  end type chiral_frame

  type(chiral_frame), save, allocatable         :: chiral_frame_list(:)

  type frame_def_list_entry
    integer :: frame_index        ! which atomic frame does this refer to
    integer :: frame_point_number ! which frame point (1 or 2)
    integer :: vector_tail_index  ! unit vector to add into point def
    integer :: vector_head_index
    integer :: num_vectors        ! no. of unit vec contribs to frame def point
  end type frame_def_list_entry

  type(frame_def_list_entry), save, allocatable :: frame_def_list(:)

  type frame
    logical             :: valid = .false.
    double precision    :: def_point1(3) = 0.d0
    double precision    :: def_point2(3) = 0.d0
    double precision    :: axis(3, 3)
  end type frame

  ! BUGBUG - Need to check if any of the 3 data structures below can
  !          be stack-allocated...

  type(frame), save, allocatable      :: frame_list(:)
  double precision, save, allocatable :: global_multipole(:,:)
  double precision, save, allocatable :: torque_field(:,:)

  integer, save :: frame_axis_order(3) = (/3,1,2/) ! gives axes with respect 
                                                   ! to def pts. Default is z
                                                   ! then x then y
  integer, parameter          :: MAXMP = 35
  ! From Tinker:
  double precision, parameter :: coulomb_const_kcal_per_mole = 332.05382d0

  public        global_multipole
  public        torque_field
  public        coulomb_const_kcal_per_mole
  public        am_mpole_torque_to_force
  public        am_mpole_local_to_global

contains

!*******************************************************************************!
! Subroutine:  am_mpole_rescale_multipoles
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_mpole_rescale_multipoles()

  implicit none

! Local variables:

  integer                       :: j, n
  double precision, parameter   :: bohr = 0.5291772083d0
  double precision, parameter   :: traced =  2.d0 / 3.d0

! traced is spherical harm expansion diff from Taylor's.
! See Stone's book.

! Now change units from Bohr to Angstroms

  do n = 1, num_multipoles

    do j = 2, 4
      local_multipole(j, n) = bohr * local_multipole(j, n)
    end do

    do j = 5, 10
      local_multipole(j, n) = traced * bohr * bohr * local_multipole(j, n)
    end do

  end do

  return

end subroutine am_mpole_rescale_multipoles

!*******************************************************************************!
! Subroutine:  am_mpole_local_to_global
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_mpole_local_to_global(crd)

  use amoeba_flags_mod
  use timers_mod

  implicit none

  ! Formal arguments:

  double precision, intent(in)  :: crd(3, *)

  ! Local variables:
  
  integer                       :: n

  if (.not. btest(do_amoeba_multipoles_flag, valid_bit)) return

  ! Clear frames

  do n = 1, num_multipoles
    frame_list(n)%def_point1 = 0.d0
    frame_list(n)%def_point2 = 0.d0
    frame_list(n)%axis = 0.d0
  end do

  if( num_frame_def_list > 0 ) call am_mpole_build_frame_def_pts(crd, frame_def_list, frame_list)

  call am_mpole_def_pts_to_axes(frame_axis_order, frame_list)

  if( num_chiral_frame_list > 0 ) call am_mpole_check_chirality(crd, frame_axis_order, chiral_frame_list, &
                                  frame_list, local_multipole)

  call am_mpole_rotate_multipole(frame_list, local_multipole, global_multipole)

  call update_pme_time(pme_misc_timer)

  return

end subroutine am_mpole_local_to_global

!*******************************************************************************!
! Subroutine:  am_mpole_torque_to_force
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_mpole_torque_to_force(atm_cnt, crd, frc, virial)

  use timers_mod

  implicit none

! Formal arguments:

  integer, intent(in)                   :: atm_cnt
  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in out)      :: frc(3, *)
  double precision, intent(in out)      :: virial(3, 3)

! Local variables:

  double precision                      :: de_drotsite(3, atm_cnt)
  double precision                      :: de_drotframe(3, atm_cnt)
  double precision                      :: de_ddefpt(3, 2, atm_cnt)

  call am_mpole_get_de_drot_mpoles(atm_cnt, global_multipole, torque_field, &
                                   de_drotsite)

  call am_mpole_accum_de_dframe_rot(atm_cnt, frame_axis_order, de_drotsite, &
                                    frame_list, de_drotframe)

  call am_mpole_accum_de_ddefpts(atm_cnt, frame_axis_order, de_drotframe, &
                                 frame_list, de_ddefpt)

  if( num_frame_def_list > 0 ) then
      call am_mpole_de_ddefpts_to_atoms(num_frame_def_list, frame_def_list, &
                                        de_ddefpt, crd, frc, virial)
  endif 
  
  call update_pme_time(pme_misc_timer)

  return

end subroutine am_mpole_torque_to_force

!*******************************************************************************!
! Subroutine:  am_mpole_build_frame_def_pts
!
! Description: routine needed by am_mpole_local_to_global,
!              am_mpole_torque_to_force
!
!*******************************************************************************

subroutine am_mpole_build_frame_def_pts(crd, fr_deflist, fr_list)

  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision, intent(in)                  :: crd(3, *)
  type(frame_def_list_entry), intent(in)        :: fr_deflist(*)
  type(frame), intent(in out)                   :: fr_list(*)

! Local variables:

  integer                                       :: n, i1, i2, j, k, m, p
  double precision                              :: dx, dy, dz, wt

  do n = 1, num_frame_def_list
    m = fr_deflist(n)%frame_index
    p = fr_deflist(n)%frame_point_number
    i1 = fr_deflist(n)%vector_tail_index
    i2 = fr_deflist(n)%vector_head_index
    k = fr_deflist(n)%num_vectors
    if (i1 .gt. 0 .and. i2 .gt. 0) then
      fr_list(m)%valid = .true.
      dx = crd(1, i2) - crd(1, i1)
      dy = crd(2, i2) - crd(2, i1)
      dz = crd(3, i2) - crd(3, i1)
      wt = k * sqrt(dx * dx + dy * dy + dz * dz)
      ! divide by length of i1i2 times num such
      if (p .eq. 1) then
        fr_list(m)%def_point1(1) = fr_list(m)%def_point1(1) + dx / wt
        fr_list(m)%def_point1(2) = fr_list(m)%def_point1(2) + dy / wt
        fr_list(m)%def_point1(3) = fr_list(m)%def_point1(3) + dz / wt
      else if (p .eq. 2) then
        fr_list(m)%def_point2(1) = fr_list(m)%def_point2(1) + dx / wt
        fr_list(m)%def_point2(2) = fr_list(m)%def_point2(2) + dy / wt
        fr_list(m)%def_point2(3) = fr_list(m)%def_point2(3) + dz / wt
      else
        error_msg =  'am_mpole_build_frame_def_pts: serious problem!'
        call mol_mech_error
      end if
    else
      fr_list(m)%valid = .false.
    end if
  end do

  return

end subroutine am_mpole_build_frame_def_pts

!*******************************************************************************!
! Subroutine:  am_mpole_def_pts_to_axes
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_mpole_def_pts_to_axes(axis_order, fr_list)

  implicit none

! Formal arguments:

  integer, intent(in)           :: axis_order(3)
  type(frame), intent(in out)   :: fr_list(*)

! Local variables:

  integer                       :: i, n, k1, k2, k3
  double precision              :: u(3), v(3), w(3), siz, dot

  k1 = axis_order(1)
  k2 = axis_order(2)
  k3 = axis_order(3)

  do n = 1, num_multipoles

    if (fr_list(n)%valid) then

  ! u is unit vector in direction of primary def pt

      u(1) = fr_list(n)%def_point1(1)
      u(2) = fr_list(n)%def_point1(2)
      u(3) = fr_list(n)%def_point1(3)
      siz = sqrt(u(1) * u(1) + u(2) * u(2) + u(3) * u(3))
      u(1) = u(1) / siz
      u(2) = u(2) / siz
      u(3) = u(3) / siz

  ! v is unit vector given by component of secondary pt orthog to u

      v(1) = fr_list(n)%def_point2(1)
      v(2) = fr_list(n)%def_point2(2)
      v(3) = fr_list(n)%def_point2(3)
      dot = u(1) * v(1) + u(2) * v(2) + u(3) * v(3)
      v(1) = v(1) - dot * u(1)
      v(2) = v(2) - dot * u(2)
      v(3) = v(3) - dot * u(3)
      siz = sqrt(v(1) * v(1) + v(2) * v(2) + v(3) * v(3))
      v(1) = v(1) / siz
      v(2) = v(2) / siz
      v(3) = v(3) / siz

  ! w is u cross v

      w(1) = u(2) * v(3) - u(3) * v(2)
      w(2) = u(3) * v(1) - u(1) * v(3)
      w(3) = u(1) * v(2) - u(2) * v(1)

  ! build axes

      do i = 1, 3
        fr_list(n)%axis(i, k1) = u(i)
        fr_list(n)%axis(i, k2) = v(i)
        fr_list(n)%axis(i, k3) = w(i)
      end do

    end if ! fr_list(n)%valid

  end do

  return

end subroutine am_mpole_def_pts_to_axes

!*******************************************************************************!
! Subroutine:  am_mpole_check_chirality
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_mpole_check_chirality(crd, axis_order, chiral_frlist, fr_list, &
                                    loc_mpole)

  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision, intent(in)          :: crd(3, *)
  integer, intent(in)                   :: axis_order(3)
  type(chiral_frame), intent(in out)    :: chiral_frlist(*)
  type(frame), intent(in)               :: fr_list(*)
  double precision, intent(in out)      :: loc_mpole(10, *)

! Local variables:

  integer                               :: n, j, k, k3
  double precision                      :: dx, dy, dz, dot

  k3 = axis_order(3)

  do n = 1, num_chiral_frame_list

    j = chiral_frlist(n)%frame_index

    if (.not. fr_list(j)%valid) then
      write(error_msg,*) 'chiral frame: serious problem in frame ', j
      call mol_mech_error
    end if

    k = chiral_frlist(n)%fourth_atom
    dx = crd(1, k) - crd(1, j)
    dy = crd(2, k) - crd(2, j)
    dz = crd(3, k) - crd(3, j)
    dot = dx * fr_list(j)%axis(1, k3)  + &
          dy * fr_list(j)%axis(2, k3)  + &
          dz * fr_list(j)%axis(3, k3)

    ! this should be negative if chirality is positive and vice versa

    if (((chiral_frlist(n)%chirality .eq. 1) .and. (dot .gt. 0)) .or. &
         ((chiral_frlist(n)%chirality .eq. -1) .and. (dot .lt. 0))) then

       chiral_frlist(n)%chirality = -chiral_frlist(n)%chirality

       if (k3 .eq. 1) then
         loc_mpole(2, n) = -loc_mpole(2, n)
         loc_mpole(8, n) = -loc_mpole(8, n)
         loc_mpole(9, n) = -loc_mpole(9, n)
       else if (k3 .eq. 2) then
         loc_mpole(3, n) = -loc_mpole(3, n)
         loc_mpole(8, n) = -loc_mpole(8, n)
         loc_mpole(10, n) = -loc_mpole(10, n)
       else
         loc_mpole(4, n) = -loc_mpole(4, n)
         loc_mpole(9, n) = -loc_mpole(9, n)
         loc_mpole(10, n) = -loc_mpole(10, n)
       end if

    end if

  end do

  return

end subroutine am_mpole_check_chirality

!*******************************************************************************!
! Subroutine:  am_mpole_rotate_multipole
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_mpole_rotate_multipole(fr_list, loc_mpole, glob_mpole)

  implicit none

! Formal arguments:

  type(frame), intent(in out)   :: fr_list(*)
  double precision, intent(in)  :: loc_mpole(10, *)
  double precision, intent(out) :: glob_mpole(10, *)

! Local variables:

  integer                       :: order, dimxy, n, j
  double precision              :: Mpole_xy(MAXMP * MAXMP)

  order = 10
  dimxy = 10

  do n = 1, num_multipoles

    if (fr_list(n)%valid ) then
      call xform_mpole_matrix(fr_list(n)%axis, Mpole_xy, order)
      call xform_mpole(Mpole_xy, dimxy, loc_mpole(:, n), glob_mpole(:, n), &
                       order)
    else
      do j = 1, 10
        glob_mpole(j, n) = loc_mpole(j, n) ! for charge only case e.g. ions
      end do
    end if
  end do

  return

end subroutine am_mpole_rotate_multipole

!*******************************************************************************!
! Subroutine:  am_mpole_get_de_drot_mpoles
!
! Description:
!
! Get the derivative of electrostatic energy with respect to infinitesmal 
! rotations of atomic frames about the x, y, z global axes
! Basic idea--electrostatic energy given by dot product of cartesian
! multipoles and the electrostatic potential and its cartesian derivatives
! i.e. ene = 1/2 * (q * phi + mux * dphidx + ...)
! Thus derivative obtained by rotating multipoles infinitesmally
!
!*******************************************************************************

subroutine am_mpole_get_de_drot_mpoles(atm_cnt, global_multipole, &
                                   torque_field, de_drotsite)

  implicit none

! Formal arguments:

  integer , intent(in)          :: atm_cnt
  double precision, intent(in)  :: global_multipole(10, *)
  double precision, intent(in)  :: torque_field(10, *)
  double precision, intent(out) :: de_drotsite(3, *)

! Local variables:

  integer                       :: i, j, k, n, dimxy, order
  double precision              :: DMP_x(MAXMP * MAXMP)
  double precision              :: DMP_y(MAXMP * MAXMP)
  double precision              :: DMP_z(MAXMP * MAXMP)
  double precision              :: A_xy(3, 3)
  double precision              :: DA_xy(3, 3)
  double precision              :: Tmp_x(10), Tmp_y(10), Tmp_z(10)

  order = 10

! to get de_drot we calculate the deriv of mpole wrt infinitesmal
! rotations about x, y and z axis

  do i = 1, 3
    do j = 1, 3
      A_xy(i, j) = 0.d0
    end do
    A_xy(i, i) = 1.d0
  end do
     
! do the maximal order
! x-axis rotation of dtheta

  do i = 1, 3
    do j = 1, 3
      DA_xy(i, j) = 0.d0
    end do
  end do

  DA_xy(3, 2) = 1.d0
  DA_xy(2, 3) = -1.d0

  call xform_mpole_deriv_matrix(A_xy, DA_xy, DMP_x, order)

! y-axis

  do i = 1, 3
    do j = 1, 3
      DA_xy(i, j) = 0.d0
    end do
  end do

  DA_xy(3, 1) = -1.d0
  DA_xy(1, 3) = 1.d0

  call xform_mpole_deriv_matrix(A_xy, DA_xy, DMP_y, order)

! z-axis

  do i = 1, 3
    do j = 1, 3
      DA_xy(i, j) = 0.d0
    end do
  end do

  DA_xy(2, 1) = 1.d0
  DA_xy(1, 2) = -1.d0

  call xform_mpole_deriv_matrix(A_xy, DA_xy, DMP_z, order)
 
  dimxy = order

  do n = 1, atm_cnt
    de_drotsite(1, n) = 0.d0
    de_drotsite(2, n) = 0.d0
    de_drotsite(3, n) = 0.d0
    call xform_mpole(DMP_x, dimxy, global_multipole(:, n), Tmp_x, order)
    call xform_mpole(DMP_y, dimxy, global_multipole(:, n), Tmp_y, order)
    call xform_mpole(DMP_z, dimxy, global_multipole(:, n), Tmp_z, order)
    do k = 1, order

      de_drotsite(1, n) = de_drotsite(1, n) +  &
        coulomb_const_kcal_per_mole * Tmp_x(k) * torque_field(k, n)
        
      de_drotsite(2, n) = de_drotsite(2, n) +  &
        coulomb_const_kcal_per_mole * Tmp_y(k) * torque_field(k, n)
        
      de_drotsite(3, n) = de_drotsite(3, n) +  &
        coulomb_const_kcal_per_mole * Tmp_z(k) * torque_field(k, n)
        
    end do
  end do

  return

end subroutine am_mpole_get_de_drot_mpoles

!*******************************************************************************!
! Subroutine:  am_mpole_accum_de_dframe_rot
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_mpole_accum_de_dframe_rot(atm_cnt, axis_order, de_drotsite, &
                                        fr_list, de_drotframe)

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  integer, intent(in)           :: axis_order(3)
  double precision, intent(in)  :: de_drotsite(3, atm_cnt)
  type(frame), intent(in)       :: fr_list(*)
  double precision, intent(out) :: de_drotframe(3, atm_cnt)

! Local variables:

  integer                       :: n, k1, k2, k3
  double precision              :: p2unit(3), siz

  k1 = axis_order(1)
  k2 = axis_order(2)
  k3 = axis_order(3)

! deriv of energy with respect to rotation about unit vectors along
! p1 p2 and their cross product
! note that unit vector along p1 corresponds to k1st frame axis
! and unit vector in p1 x p2 direction corresponds to k3rd frame axis
! the energy derivative with respect to rotation about any unit vector
! is for each mpole given by the dot product of de_drotpole 
! (which is negative of torque due to that mpole) with the unit vector
        
  do n = 1, atm_cnt
    ! initialize
    de_drotframe(1, n) = 0.d0
    de_drotframe(2, n) = 0.d0
    de_drotframe(3, n) = 0.d0
    if (fr_list(n)%valid) then
      siz = sqrt(fr_list(n)%def_point2(1) * fr_list(n)%def_point2(1) + &
                 fr_list(n)%def_point2(2) * fr_list(n)%def_point2(2) + &
                 fr_list(n)%def_point2(3) * fr_list(n)%def_point2(3))
      p2unit(1) = fr_list(n)%def_point2(1) / siz
      p2unit(2) = fr_list(n)%def_point2(2) / siz
      p2unit(3) = fr_list(n)%def_point2(3) / siz
      de_drotframe(k1, n) = de_drotframe(k1, n) +   &
                           de_drotsite(1, n) * fr_list(n)%axis(1, k1) +  &
                           de_drotsite(2, n) * fr_list(n)%axis(2, k1) +  &
                           de_drotsite(3, n) * fr_list(n)%axis(3, k1)
      de_drotframe(k2, n) = de_drotframe(k2, n) +   &
                           de_drotsite(1, n) * p2unit(1) +   &
                           de_drotsite(2, n) * p2unit(2) +   &
                           de_drotsite(3, n) * p2unit(3)
      de_drotframe(k3, n) = de_drotframe(k3, n) +    &
                           de_drotsite(1, n) * fr_list(n)%axis(1, k3) +  &
                           de_drotsite(2, n) * fr_list(n)%axis(2, k3) +  &
                           de_drotsite(3, n) * fr_list(n)%axis(3, k3)
    end if ! fr_list(n)%valid
  end do
      
  return

end subroutine am_mpole_accum_de_dframe_rot

!*******************************************************************************!
! Subroutine:  am_mpole_accum_de_ddefpts
!
! Description:
!
! get the derivs of energy with respect to movement of defpoints
! expressed in the local frame coord system
!
!*******************************************************************************

subroutine am_mpole_accum_de_ddefpts(atm_cnt, axis_order, &
                                     de_drotframe, fr_list, de_ddefpt)

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  integer, intent(in)           :: axis_order(3)
  double precision, intent(in)  :: de_drotframe(3, atm_cnt)
  type(frame), intent(in)       :: fr_list(*)
  double precision, intent(out) :: de_ddefpt(3, 2, atm_cnt)

! Local variables:

  integer                       :: n, k, j, k1, k2, k3
  double precision              :: p1(3), p2(3), p2unit(3)
  double precision              :: p2perp1(3), p1perp2(3)
  double precision              :: u(3), v(3), w(3)
  double precision              :: dotu, dotv, dotw
  double precision              :: sizp1perp2, sizp2perp1
  double precision              :: dot12, dot21
  double precision              :: sizp1, sizp2
  double precision              :: dedrotp1, dedrotp2
  double precision              :: dedu, dedv, dedw
  double precision              :: de_drotu, de_drotv, de_drotw

  k1 = axis_order(1)
  k2 = axis_order(2)
  k3 = axis_order(3)
  do n = 1, atm_cnt
    if (fr_list(n)%valid) then
      do j = 1, 3
        p1(j) = fr_list(n)%def_point1(j)
        p2(j) = fr_list(n)%def_point2(j)
        u(j) = fr_list(n)%axis(j, k1)
        v(j) = fr_list(n)%axis(j, k2)
        w(j) = fr_list(n)%axis(j, k3)
      end do
      de_drotu = de_drotframe(k1, n)
      de_drotv = de_drotframe(k2, n)
      de_drotw = de_drotframe(k3, n)

      sizp1 = sqrt(p1(1)**2 + p1(2)**2 + p1(3)**2)
      sizp2 = sqrt(p2(1)**2 + p2(2)**2 + p2(3)**2)
      do j = 1, 3
        p2unit(j) = p2(j) / sizp2
!       p1unit(j) = u(j) so no need to recalculate
      end do
      dot21 = u(1) * p2(1) + u(2) * p2(2) + u(3) * p2(3)
      dot12 = p1(1) * p2unit(1) + p1(2) * p2unit(2) + p1(3) * p2unit(3)
      do j = 1, 3
        p2perp1(j) = p2(j) - dot21 * u(j)
        p1perp2(j) = p1(j) - dot12 * p2unit(j)
      end do
      sizp2perp1 = sqrt(p2perp1(1)**2 + p2perp1(2)**2 + p2perp1(3)**2)
      sizp1perp2 = sqrt(p1perp2(1)**2 + p1perp2(2)**2 + p1perp2(3)**2)

! def point one is along axis one. movement du parallel to that axis does
! not rotate the frame..so deriv is zero
!       dedu = 0.d0
! movement dv in v-axis direction corresponds to rotation about local w-axis
! of dtheta = dv/sizp1; thus a change in energy of 
!    dE = dedrotw * dtheta = dedrotw * dv/sizp1
!    dE/dv = dedrotw /sizp1

      dedv = de_drotw / sizp1

! movement dw in w-axis direction corresponds to rotation about p2unit
! of dtheta = -dw/sizp1perp2 (clockwise rotation) 

      dedw = -de_drotv / sizp1perp2

! Now convert to derivs wrt x, y, z. u = p.u = x * u(1) + y * u(2) + z * u(3)
! thus dudx = u(1). similaryl dvdx = v(1)

      de_ddefpt(1, 1, n) = dedv * v(1) + dedw * w(1)
      de_ddefpt(2, 1, n) = dedv * v(2) + dedw * w(2)
      de_ddefpt(3, 1, n) = dedv * v(3) + dedw * w(3)
 
! for point 2..any movement in the local uv plane does not affect the frame
!       dedu = 0.d0
!       dedv = 0.d0
! movement dw in w direction corresponds to rotation about local u-axis 
! of dtheta = dw/sizp2perpu

      dedw = de_drotu/sizp2perp1
      de_ddefpt(1, 2, n) = dedw * w(1)
      de_ddefpt(2, 2, n) = dedw * w(2)
      de_ddefpt(3, 2, n) = dedw * w(3)
    else ! fr_list(n)%valid .eq. .false.
      do j = 1, 3
        de_ddefpt(j, 1, n) = 0.d0
        de_ddefpt(j, 2, n) = 0.d0
      end do
    end if ! fr_list(n)%valid
  end do ! n = 1, atm_cnt

  return

end subroutine am_mpole_accum_de_ddefpts

!*******************************************************************************!
! Subroutine:  am_mpole_de_ddefpts_to_atoms
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_mpole_de_ddefpts_to_atoms(num_fr_deflist, fr_deflist, &
                                        de_ddefpt, crd, frc, virial)
  implicit none

! Formal arguments:

  integer, intent(in)                    :: num_fr_deflist
  type(frame_def_list_entry), intent(in) :: fr_deflist(*)
  double precision, intent(in)           :: de_ddefpt(3, 2, *)
  double precision, intent(in)           :: crd(3, *)
  double precision, intent(in out)       :: frc(3, *)
  double precision, intent(in out)       :: virial(3, 3)

! Local variables:

  integer                                :: n, i, j, k, l, m, p
  double precision                       :: siz, siz2
  double precision                       :: dx, dy, dz
  double precision                       :: dux_dx, dux_dy, dux_dz
  double precision                       :: duy_dy, duy_dz, duz_dz
  double precision                       :: dedux, deduy, deduz
  double precision                       :: dedx, dedy, dedz
  double precision                       :: siz3inv, sizinv
  double precision                       :: vxx, vxy, vxz
  double precision                       :: vyx, vyy, vyz
  double precision                       :: vzx, vzy, vzz

  vxx = 0.d0
  vxy = 0.d0
  vxz = 0.d0
  vyx = 0.d0
  vyy = 0.d0
  vyz = 0.d0
  vzx = 0.d0
  vzy = 0.d0
  vzz = 0.d0

  do n = 1, num_fr_deflist
    i = fr_deflist(n)%vector_tail_index
    j = fr_deflist(n)%vector_head_index
    if (i .gt. 0 .and. j .gt. 0) then
      k = fr_deflist(n)%frame_index
      l = fr_deflist(n)%frame_point_number
      m = fr_deflist(n)%num_vectors
      dx = crd(1, j) - crd(1, i)
      dy = crd(2, j) - crd(2, i)
      dz = crd(3, j) - crd(3, i)
      siz2 = dx * dx + dy * dy + dz * dz
      siz = sqrt(siz2)
      siz3inv = 1.d0 / (siz2 * siz)
      sizinv = 1.d0 / siz

! ux, uy, uz are given by dx / siz, dy / siz, and dz / siz 

      dux_dx = sizinv - dx * dx * siz3inv
      dux_dy = -dx * dy * siz3inv   ! note duy_dx = dux_dy use this below
      dux_dz = -dx * dz * siz3inv
      duy_dy = sizinv - dy * dy * siz3inv
      duy_dz = -dy * dz * siz3inv
      duz_dz = sizinv - dz * dz * siz3inv

! the derivs of E wrt coordinates of unit vector in ij direction are given
! by (1/m) times the derivs of E wrt coords of def point (l, k)
! since the def point is the simple average of m of these unit vectors

      dedux = de_ddefpt(1, l, k) / m
      deduy = de_ddefpt(2, l, k) / m
      deduz = de_ddefpt(3, l, k) / m

! now apply chain rule, using symmetry e.g. dux_dy = duy_dx

      dedx = dedux * dux_dx + deduy * dux_dy + deduz * dux_dz
      dedy = dedux * dux_dy + deduy * duy_dy + deduz * duy_dz
      dedz = dedux * dux_dz + deduy * duy_dz + deduz * duz_dz

! finally apply forces. note force is negative of deriv of energy wrt position
! also note e.g. deriv of dx wrt x position of atoms i and j is -1,  + 1

      frc(1, i) = frc(1, i) + dedx
      frc(2, i) = frc(2, i) + dedy
      frc(3, i) = frc(3, i) + dedz
      frc(1, j) = frc(1, j) - dedx
      frc(2, j) = frc(2, j) - dedy
      frc(3, j) = frc(3, j) - dedz
      vxx = vxx + dedx * dx
      vxy = vxy + dedx * dy
      vxz = vxz + dedx * dz
      vyx = vyx + dedy * dx
      vyy = vyy + dedy * dy
      vyz = vyz + dedy * dz
      vzx = vzx + dedz * dx
      vzy = vzy + dedz * dy
      vzz = vzz + dedz * dz
    end if ! i > 0 .and. j > 0
  end do ! n = 1, num_fr_deflist

  virial(1, 1) = virial(1, 1) + vxx
  virial(1, 2) = virial(1, 2) + 0.5d0 * (vxy + vyx)
  virial(1, 3) = virial(1, 3) + 0.5d0 * (vxz + vzx)
  virial(2, 1) = virial(2, 1) + 0.5d0 * (vxy + vyx)
  virial(2, 2) = virial(2, 2) + vyy
  virial(2, 3) = virial(2, 3) + 0.5d0 * (vyz + vzy)
  virial(3, 1) = virial(3, 1) + 0.5d0 * (vxz + vzx)
  virial(3, 2) = virial(3, 2) + 0.5d0 * (vyz + vzy)
  virial(3, 3) = virial(3, 3) + vzz

  return

end subroutine am_mpole_de_ddefpts_to_atoms


subroutine set_local_multipoles_list(list_new,n)
	
	implicit none
	double precision, intent(in)  :: list_new(10,n)
	integer, intent(in)           :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	nold = size(local_multipole,2)
	num_multipoles = n
	
    if( nold .ne. n ) then
        if( nold .gt. 0) then 
            deallocate(local_multipole)
            deallocate(frame_list)
            deallocate(global_multipole)
            deallocate(torque_field)
        endif
        allocate( local_multipole(10,num_multipoles), &
                  frame_list(num_multipoles), &
                  global_multipole(10, num_multipoles), &
                  torque_field(10, num_multipoles), &
                  stat = alloc_failed)
        
        if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	if( num_multipoles .eq. 0) return 
	    
	local_multipole(:,1:n) = list_new(:,1:n) 
	
	call am_mpole_rescale_multipoles()
	
end subroutine set_local_multipoles_list

subroutine set_chiral_frame_list(list_new,n)
	
	implicit none
	integer, intent(in)         :: list_new(3,n)
	integer, intent(in)         :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	integer          ::   i
	
	nold = size(chiral_frame_list)
	num_chiral_frame_list = n
	
    if( nold .ne. n ) then
        if( nold .gt. 0) then 
            deallocate(chiral_frame_list)
        endif
        allocate( chiral_frame_list(num_chiral_frame_list), &
                  stat = alloc_failed)
        
        if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	if( num_chiral_frame_list .eq. 0) return 
	    
	do i = 1, num_chiral_frame_list
	    chiral_frame_list(i)%frame_index = list_new(1,i)
	    chiral_frame_list(i)%fourth_atom = list_new(2,i) 
	    chiral_frame_list(i)%chirality   = list_new(3,i)  
	enddo
	
end subroutine set_chiral_frame_list

subroutine set_reg_frame_list(list_new,n)
	
	implicit none
	integer, intent(in)         :: list_new(5,n)
	integer, intent(in)         :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	integer          ::   i
	
	nold = size(frame_def_list)
	num_frame_def_list = n
	
    if( nold .ne. n ) then
        if( nold .gt. 0) then 
            deallocate(frame_def_list)
        endif
        allocate( frame_def_list(num_frame_def_list), stat = alloc_failed)
        
        if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	if( num_frame_def_list .eq. 0) return 
	    
	do i = 1, num_frame_def_list
	    frame_def_list(i)%frame_index        = list_new(1,i)
	    frame_def_list(i)%frame_point_number = list_new(2,i)
	    frame_def_list(i)%vector_tail_index  = list_new(3,i)
	    frame_def_list(i)%vector_head_index  = list_new(4,i)
	    frame_def_list(i)%num_vectors        = list_new(5,i)
	enddo
	
end subroutine set_reg_frame_list

subroutine set_valid_bit(ival)

    use amoeba_flags_mod
    
    implicit none
	integer, intent(in)         :: ival
    
    if(ival .eq. 1) then
        do_amoeba_multipoles_flag = ibset(do_amoeba_multipoles_flag, valid_bit)
    else
        do_amoeba_multipoles_flag = ibclr(do_amoeba_multipoles_flag, valid_bit)
    endif
        
end subroutine set_valid_bit

subroutine get_global_multipole( natom, global_multipole_out )

implicit none

    integer,          intent(in)   :: natom
    double precision, intent(out)  :: global_multipole_out(*)

    integer          ::   i

    if ( 10*natom > size(global_multipole) ) return

	do i = 1, natom
        global_multipole_out(10*i - 9 ) = global_multipole(1,i)
        global_multipole_out(10*i - 8 ) = global_multipole(2,i)
        global_multipole_out(10*i - 7 ) = global_multipole(3,i)
        global_multipole_out(10*i - 6 ) = global_multipole(4,i)
        global_multipole_out(10*i - 5 ) = global_multipole(5,i)
        global_multipole_out(10*i - 4 ) = global_multipole(6,i)
        global_multipole_out(10*i - 3 ) = global_multipole(7,i)
        global_multipole_out(10*i - 2 ) = global_multipole(8,i)
        global_multipole_out(10*i - 1 ) = global_multipole(9,i)
        global_multipole_out(10*i)      = global_multipole(10,i)
    enddo 

end subroutine get_global_multipole

end module amoeba_multipoles_mod

!*******************************************************************************!
! Subroutine:   xform_mpole_matrix  
!
! Description: 
!
! A_xy is matrix of coefficients expanding y in terms of x
! i.e. y_i = A_i1 * x_1 + A_i2 * x_2 + A_i3 * x_3
! Mpole_xy is resulting matrix getting expansion of multipoles
! with basis y in terms of those with basis x
! order is order of multipoles 1 for charge, 4 for charge-dipoles
! 10 for charge, dipole, quadrupole, 20 for up to octupole and
! finally 35 for up to hexadecapole
!
! terms are obtained by examining taylor expansions of energies in terms
! of multipoles times field derivatives in coord system x or y,
! employing the chain rule to transform field derivatives as in
! subroutine xform_mpole_field_matrix below, reversing summation order and
! collecting terms. Note that symmetry applies, so that for example
! mpole_x(8) is the sum of the x1x2 and x2x1 multipole components
! this symmetry is applied in the matrix
!
!*******************************************************************************

subroutine xform_mpole_matrix(A_xy, Mpole_xy, order)
 
  implicit none

! Mpole 1 is charge... Mpole 10 is x_3, x_3 quadrupolar coefficients etc.

! Formal arguments:

  integer               :: order
  double precision      :: A_xy(3, 3)
  double precision      :: Mpole_xy(order, order)

! Local variables:

  integer               :: i, j, k, l, m, n, p, ind1, ind2, jj, kk
  integer               :: qind1(6), qind2(6)
  integer               :: oind1(10), oind2(10), oind3(10)
  integer               :: hind1(15), hind2(15), hind3(15), hind4(15)

  common /MxfInd/ qind1, qind2, oind1, oind2, oind3, hind1, hind2, hind3, hind4

  data qind1   /1, 2, 3, 1, 1, 2/
  data qind2   /1, 2, 3, 2, 3, 3/
  data oind1   /1, 2, 3, 1, 1, 2, 2, 3, 3, 1/
  data oind2   /1, 2, 3, 1, 1, 2, 2, 3, 3, 2/
  data oind3   /1, 2, 3, 2, 3, 1, 3, 1, 2, 3/
  data hind1   /1, 2, 3, 1, 1, 2, 2, 3, 3, 1, 1, 2, 1, 2, 3/
  data hind2   /1, 2, 3, 1, 1, 2, 2, 3, 3, 1, 1, 2, 1, 2, 3/
  data hind3   /1, 2, 3, 1, 1, 2, 2, 3, 3, 2, 3, 3, 2, 1, 1/
  data hind4   /1, 2, 3, 2, 3, 1, 3, 1, 2, 2, 3, 3, 3, 3, 2/
      

! CHARGE case

  Mpole_xy(1, 1) = 1.d0

  if (order .eq. 1) return

  do j = 2, order
   Mpole_xy(1, j) = 0.d0
   Mpole_xy(j, 1) = 0.d0
  end do

! DIPOLES
! d_yk = A_xy(k, 1) * d_x1 + A_xy(k, 2) * d_x2 + A_xy(k, 3) * d_x3
! D'_j

  do j = 2, 4
    do k = 2, 4
      Mpole_xy(j, k) = A_xy(j-1, k-1)
    end do
  end do

  if (order .eq. 4) return

  do j = 2, 4
    do k = 5, order
      Mpole_xy(j, k) = 0.d0
      Mpole_xy(k, j) = 0.d0
    end do
  end do

! QUADRUPOLES
! q_ykyl = sum_i, j (A_xy(k, i) * A_xy(l, j) + A_xy(k, j) * A_xy(l, i)) * 
! (q_xixj + q_xjxi)
! Mp(5) = q_y1y1, .., Mp(7) = q_y3y3, Mp(8) = q_y1y2 + q_y2y1, .,
! Mp(10) = q_y2y3 + q_y3y2
! Q'_kk

  do ind1 = 1, 3
    k = qind1(ind1)
    do ind2 = 1, 6
      i = qind1(ind2)
      j = qind2(ind2)
      Mpole_xy(ind1 + 4, ind2 + 4) = A_xy(k, i) * A_xy(k, j)
    end do
  end do

! Q'_kl

  do ind1 = 4, 6
    k = qind1(ind1)
    l = qind2(ind1)
    do ind2 = 1, 6
      i = qind1(ind2)
      j = qind2(ind2)
      Mpole_xy(ind1 + 4, ind2 + 4) = A_xy(k, i) * A_xy(l, j) + &
                                     A_xy(k, j) * A_xy(l, i)
    end do
  end do

  if (order .eq. 10) return

  do j = 5, 10
    do k = 11, order
      Mpole_xy(k, j) = 0.d0
      Mpole_xy(j, k) = 0.d0
    end do
  end do

! OCTUPOLES
! O'_lll

  do ind1 = 1, 3
    l = oind1(ind1)
    do ind2 = 1, 10
      i = oind1(ind2)
      j = oind2(ind2)
      k = oind3(ind2)
      Mpole_xy(ind1 + 10, ind2 + 10) = A_xy(l, i) * A_xy(l, j) * A_xy(l, k)
    end do
  end do

! O'_llm

  do ind1 = 4, 9
    l = oind1(ind1)
    m = oind3(ind1)
    do ind2 = 1, 9
      i = oind1(ind2)
      k = oind3(ind2)
      Mpole_xy(ind1 + 10, ind2 + 10) = A_xy(l, i) * A_xy(l, i) * A_xy(m, k) + &
        2.d0 * A_xy(l, i) * A_xy(m, i) * A_xy(l, k)
    end do
    Mpole_xy(ind1 + 10, 20) = A_xy(l, 1) * A_xy(l, 2) * A_xy(m, 3) + &
                              A_xy(l, 1) * A_xy(m, 2) * A_xy(l, 3) + &
                              A_xy(m, 1) * A_xy(l, 2) * A_xy(l, 3)
  end do

! O'_123

  Mpole_xy(20, 11) = 6.d0 * A_xy(1, 1) * A_xy(2, 1) * A_xy(3, 1)
  Mpole_xy(20, 12) = 6.d0 * A_xy(1, 2) * A_xy(2, 2) * A_xy(3, 2)
  Mpole_xy(20, 13) = 6.d0 * A_xy(1, 3) * A_xy(2, 3) * A_xy(3, 3)

  do ind2 = 4, 9
    i = oind1(ind2)
    k = oind3(ind2)
    Mpole_xy(20, 10 + ind2) = 2.d0 * (A_xy(1, i) * A_xy(2, i) * A_xy(3, k) + &
                              A_xy(1, i) * A_xy(2, k) * A_xy(3, i) + &
                              A_xy(1, k) * A_xy(2, i) * A_xy(3, i))
  end do

  Mpole_xy(20, 20) = A_xy(1, 1) * A_xy(2, 2) * A_xy(3, 3) + &
                     A_xy(1, 1) * A_xy(3, 2) * A_xy(2, 3) + &
                     A_xy(2, 1) * A_xy(1, 2) * A_xy(3, 3) + &
                     A_xy(2, 1) * A_xy(3, 2) * A_xy(1, 3) + &
                     A_xy(3, 1) * A_xy(1, 2) * A_xy(2, 3) + &
                     A_xy(3, 1) * A_xy(2, 2) * A_xy(1, 3) 

  if (order .eq. 20) return

  do j = 11, 20
    do k = 21, order
      Mpole_xy(k, j) = 0.d0
      Mpole_xy(j, k) = 0.d0
    end do
  end do

! HEXADECAPOLES
! H'_mmmm

  do ind1 = 1, 3
    m = hind1(ind1)
    do ind2 = 1, 15
      i = hind1(ind2)
      j = hind2(ind2)
      k = hind3(ind2)
      l = hind4(ind2)
      Mpole_xy(ind1 + 20, ind2 + 20) = &
        A_xy(m, i) * A_xy(m, j) * A_xy(m, k) * A_xy(m, l)
    end do
  end do
! H'_mmmn
  do ind1 = 4, 9
    m = hind1(ind1)
    n = hind4(ind1)
! can put iiii & iiij together
    do ind2 = 1, 9
      i = hind1(ind2)
      j = hind4(ind2)
      Mpole_xy(ind1 + 20, ind2 + 20) = &
        A_xy(m, i) * A_xy(m, i) * A_xy(m, i) * A_xy(n, j) + &
        3.d0 * A_xy(m, i) * A_xy(m, i) * A_xy(n, i) * A_xy(m, j) 
    end do
! can put iijj & iijk together
    do ind2 = 10, 15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Mpole_xy(ind1 + 20, ind2 + 20) = &
        A_xy(m, i) * A_xy(m, i) * A_xy(m, j) * A_xy(n, k) + &
        A_xy(m, i) * A_xy(m, i) * A_xy(m, k) * A_xy(n, j) + &
        2.d0 * A_xy(m, i) * A_xy(m, j) * A_xy(m, k) * A_xy(n, i) 
    end do
  end do

! H'_mmnn

  do ind1 = 10, 12
    m = hind1(ind1)
    n = hind3(ind1)
! can put iiii & iiij together
    do ind2 = 1, 9
      i = hind1(ind2)
      j = hind4(ind2)
      Mpole_xy(ind1 + 20, ind2 + 20) = &
        3.d0 * (A_xy(m, i) * A_xy(m, i) * A_xy(n, i) * A_xy(n, j) + &
        A_xy(m, i) * A_xy(m, j) * A_xy(n, i) * A_xy(n, i))
    end do
! can put iijj & iijk together
    do ind2 = 10, 15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Mpole_xy(ind1 + 20, ind2 + 20) = &
        A_xy(m, i) * A_xy(m, i) * A_xy(n, j) * A_xy(n, k) + &
        A_xy(m, j) * A_xy(m, k) * A_xy(n, i) * A_xy(n, i) + &
        2.d0 * A_xy(m, i) * A_xy(m, j) * A_xy(n, i) * A_xy(n, k) + &
        2.d0 * A_xy(m, i) * A_xy(m, k) * A_xy(n, i) * A_xy(n, j) 
    end do
  end do

! H'_mmnp

  do ind1 = 13, 15
    m = hind1(ind1)
    n = hind3(ind1)
    p = hind4(ind1)
! can put iiii & iiij together
    do ind2 = 1, 9
      i = hind1(ind2)
      j = hind4(ind2)
      Mpole_xy(ind1 + 20, ind2 + 20) = &
              3.d0 * A_xy(m, i) * A_xy(m, i) * A_xy(n, i) * A_xy(p, j) + &
              3.d0 * A_xy(m, i) * A_xy(m, i) * A_xy(n, j) * A_xy(p, i) + &
              6.d0 * A_xy(m, i) * A_xy(m, j) * A_xy(n, i) * A_xy(p, i) 
    end do

! can put iijj & iijk together

    do ind2 = 10, 15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Mpole_xy(ind1 + 20, ind2 + 20) = &
        A_xy(m, i) * A_xy(m, i) * A_xy(n, j) * A_xy(p, k) + &
        A_xy(m, i) * A_xy(m, i) * A_xy(n, k) * A_xy(p, j) + &
        2.d0 * (A_xy(m, i) * A_xy(m, j) * A_xy(n, i) * A_xy(p, k) + &
        A_xy(m, i) * A_xy(m, j) * A_xy(n, k) * A_xy(p, i) + &
        A_xy(m, i) * A_xy(m, k) * A_xy(n, i) * A_xy(p, j) + &
        A_xy(m, i) * A_xy(m, k) * A_xy(n, j) * A_xy(p, i) + &
        A_xy(m, j) * A_xy(m, k) * A_xy(n, i) * A_xy(p, i))

    end do
  end do

  return

end subroutine xform_mpole_matrix

!*******************************************************************************!
! Subroutine:  xform_mpole_deriv_matrix
!
! Description:
!
! calculate derivative of Mpole_xy as DMpole_xy
! A_xy is matrix of coefficients expanding y in terms of x
! i.e. y_i = A_i1 * x_1 + A_i2 * x_2 + A_i3 * x_3
! DA_xy is deriv of A_xy wrt to some parameter
! DMpole_xy is deriv of Mpole_xy wrt to the same parameter
! order is order of multipoles 1 for charge, 4 for charge-dipoles
! 10 for charge, dipole, quadrupole, 20 for up to octupole and
! finally 35 for up to hexadecapole
!
!*******************************************************************************

subroutine xform_mpole_deriv_matrix(A_xy, DA_xy, DMpole_xy, order)
 
! Mpole 1 is charge... Mpole 10 is x_3, x_3 quadrupolar coefficients etc.

  implicit none

! Formal arguments:

  integer               :: order
  double precision      :: A_xy(3, 3)
  double precision      :: DA_xy(3, 3)
  double precision      :: DMpole_xy(order, order)

! Local variables:

  integer               :: i, j, k, l, m, n, p, ind1, ind2, jj, kk
  integer               :: qind1(6), qind2(6)
  integer               :: oind1(10), oind2(10), oind3(10)
  integer               :: hind1(15), hind2(15), hind3(15), hind4(15)

  common /MxfInd/ qind1, qind2, oind1, oind2, oind3, hind1, hind2, hind3, hind4

! CHARGE case
!     Mpole_xy(1, 1) = 1.d0

  DMpole_xy(1, 1) = 0.d0
  if (order .eq. 1) return
  do j = 2, order
    DMpole_xy(1, j) = 0.d0
    DMpole_xy(j, 1) = 0.d0
  end do

! DIPOLES
! d_yk = A_xy(k, 1) * d_x1 + A_xy(k, 2) * d_x2 + A_xy(k, 3) * d_x3
! D'_j

  do j = 2, 4
    do k = 2, 4
!     Mpole_xy(j, k) = A_xy(j-1, k-1)
      DMpole_xy(j, k) = DA_xy(j-1, k-1)
    end do
  end do

  if (order .eq. 4) return

  do j = 2, 4
    do k = 5, order
      DMpole_xy(j, k) = 0.d0
      DMpole_xy(k, j) = 0.d0
    end do
  end do

! QUADRUPOLES
! q_ykyl = sum_i, j (A_xy(k, i) * A_xy(l, j) + A_xy(k, j) * A_xy(l, i)) * (q_xixj + q_xjxi)
! Mp(5) = q_y1y1, .., Mp(7) = q_y3y3, Mp(8) = q_y1y2 + q_y2y1, ., Mp(10) = q_y2y3 + q_y3y2
! Q'_kk

  do ind1 = 1, 3
    k = qind1(ind1)
    do ind2 = 1, 6
      i = qind1(ind2)
      j = qind2(ind2)
!     Mpole_xy(ind1 + 4, ind2 + 4) = A_xy(k, i) * A_xy(k, j)
      DMpole_xy(ind1 + 4, ind2 + 4) = DA_xy(k, i) * A_xy(k, j) + &
                                      A_xy(k, i) * DA_xy(k, j)
    end do
  end do

! Q'_kl

  do ind1 = 4, 6
    k = qind1(ind1)
    l = qind2(ind1)
    do ind2 = 1, 6
      i = qind1(ind2)
      j = qind2(ind2)
!     Mpole_xy(ind1 + 4, ind2 + 4) =   &
!            A_xy(k, i) * A_xy(l, j) + A_xy(k, j) * A_xy(l, i)
      DMpole_xy(ind1 + 4, ind2 + 4) =     &
            DA_xy(k, i) * A_xy(l, j) + DA_xy(k, j) * A_xy(l, i) +    &
             A_xy(k, i) * DA_xy(l, j) + A_xy(k, j) * DA_xy(l, i)
    end do
  end do

  if (order .eq. 10) return

  do j = 5, 10
    do k = 11, order
      DMpole_xy(k, j) = 0.d0
      DMpole_xy(j, k) = 0.d0
    end do
  end do

! OCTUPOLES
! O'_lll

  do ind1 = 1, 3
    l = oind1(ind1)
    do ind2 = 1, 10
      i = oind1(ind2)
      j = oind2(ind2)
      k = oind3(ind2)
!     Mpole_xy(ind1 + 10, ind2 + 10) = A_xy(l, i) * A_xy(l, j) * A_xy(l, k)
      DMpole_xy(ind1 + 10, ind2 + 10) =    &
                DA_xy(l, i) * A_xy(l, j) * A_xy(l, k) +    &
                A_xy(l, i) * DA_xy(l, j) * A_xy(l, k) +    &
                A_xy(l, i) * A_xy(l, j) * DA_xy(l, k) 
    end do
  end do

! O'_llm

  do ind1 = 4, 9
    l = oind1(ind1)
    m = oind3(ind1)
    do ind2 = 1, 9
      i = oind1(ind2)
      k = oind3(ind2)
!     Mpole_xy(ind1 + 10, ind2 + 10) =   &
!          A_xy(l, i) * A_xy(l, i) * A_xy(m, k) +    &
!          2.d0 * A_xy(l, i) * A_xy(m, i) * A_xy(l, k)
      DMpole_xy(ind1 + 10, ind2 + 10) = &
        DA_xy(l, i) * A_xy(l, i) * A_xy(m, k) + &
        2.d0 * DA_xy(l, i) * A_xy(m, i) * A_xy(l, k)

      DMpole_xy(ind1 + 10, ind2 + 10) = DMpole_xy(ind1 + 10, ind2 + 10) + &
        A_xy(l, i) * DA_xy(l, i) * A_xy(m, k) + &
        2.d0 * A_xy(l, i) * DA_xy(m, i) * A_xy(l, k)

      DMpole_xy(ind1 + 10, ind2 + 10) = DMpole_xy(ind1 + 10, ind2 + 10) + &
        A_xy(l, i) * A_xy(l, i) * DA_xy(m, k) + &
        2.d0 * A_xy(l, i) * A_xy(m, i) * DA_xy(l, k)
    end do
!   Mpole_xy(ind1 + 10, 20) =
!          A_xy(l, 1) * A_xy(l, 2) * A_xy(m, 3) +
!          A_xy(l, 1) * A_xy(m, 2) * A_xy(l, 3) +
!          A_xy(m, 1) * A_xy(l, 2) * A_xy(l, 3)
    DMpole_xy(ind1 + 10, 20) = &
      DA_xy(l, 1) * A_xy(l, 2) * A_xy(m, 3) + &
      DA_xy(l, 1) * A_xy(m, 2) * A_xy(l, 3) + &
      DA_xy(m, 1) * A_xy(l, 2) * A_xy(l, 3)

    DMpole_xy(ind1 + 10, 20) = DMpole_xy(ind1 + 10, 20) + &
      A_xy(l, 1) * DA_xy(l, 2) * A_xy(m, 3) + &
      A_xy(l, 1) * DA_xy(m, 2) * A_xy(l, 3) + &
      A_xy(m, 1) * DA_xy(l, 2) * A_xy(l, 3)

    DMpole_xy(ind1 + 10, 20) = DMpole_xy(ind1 + 10, 20) + &
      A_xy(l, 1) * A_xy(l, 2) * DA_xy(m, 3) +   &
      A_xy(l, 1) * A_xy(m, 2) * DA_xy(l, 3) +   &
      A_xy(m, 1) * A_xy(l, 2) * DA_xy(l, 3)

  end do

! O'_123

! Mpole_xy(20, 11) = 6.d0 * A_xy(1, 1) * A_xy(2, 1) * A_xy(3, 1)
  DMpole_xy(20, 11) = 6.d0 * DA_xy(1, 1) * A_xy(2, 1) * A_xy(3, 1)

  DMpole_xy(20, 11) = DMpole_xy(20, 11) + &
                      6.d0 * A_xy(1, 1) * DA_xy(2, 1) * A_xy(3, 1)

  DMpole_xy(20, 11) = DMpole_xy(20, 11) + &
                      6.d0 * A_xy(1, 1) * A_xy(2, 1) * DA_xy(3, 1)

! Mpole_xy(20, 12) = 6.d0 * A_xy(1, 2) * A_xy(2, 2) * A_xy(3, 2)

  DMpole_xy(20, 12) = 6.d0 * DA_xy(1, 2) * A_xy(2, 2) * A_xy(3, 2)

  DMpole_xy(20, 12) = DMpole_xy(20, 12) + &
                      6.d0 * A_xy(1, 2) * DA_xy(2, 2) * A_xy(3, 2)

  DMpole_xy(20,12) = DMpole_xy(20, 12) + &
                     6.d0 * A_xy(1, 2) * A_xy(2, 2) * DA_xy(3, 2)

! Mpole_xy(20, 13) = 6.d0 * A_xy(1, 3) * A_xy(2, 3) * A_xy(3, 3)

  DMpole_xy(20, 13) = 6.d0 * DA_xy(1, 3) * A_xy(2, 3) * A_xy(3, 3)

  DMpole_xy(20, 13) = DMpole_xy(20, 13) + &
                      6.d0 * A_xy(1, 3) * DA_xy(2, 3) * A_xy(3, 3)

  DMpole_xy(20, 13) = DMpole_xy(20, 13) + &
                      6.d0 * A_xy(1, 3) * A_xy(2, 3) * DA_xy(3, 3)
  do ind2 = 4, 9
    i = oind1(ind2)
    k = oind3(ind2)
!   Mpole_xy(20, 10 + ind2) = 2.d0*
!          (A_xy(1, i) * A_xy(2, i) * A_xy(3, k) +
!             A_xy(1, i) * A_xy(2, k) * A_xy(3, i) +
!             A_xy(1, k) * A_xy(2, i) * A_xy(3, i))
    DMpole_xy(20, 10 + ind2) = &
      2.d0 * (DA_xy(1, i) * A_xy(2, i) * A_xy(3, k) + &
      DA_xy(1, i) * A_xy(2, k) * A_xy(3, i) + &
      DA_xy(1, k) * A_xy(2, i) * A_xy(3, i))

    DMpole_xy(20, 10 + ind2) = &
      DMpole_xy(20, 10 + ind2) + &
      2.d0 * (A_xy(1, i) * DA_xy(2, i) * A_xy(3, k) + &
      A_xy(1, i) * DA_xy(2, k) * A_xy(3, i) +    &
      A_xy(1, k) * DA_xy(2, i) * A_xy(3, i))

    DMpole_xy(20, 10 + ind2) = &
      DMpole_xy(20, 10 + ind2) + &
      2.d0 * (A_xy(1, i) * A_xy(2, i) * DA_xy(3, k) + &
      A_xy(1, i) * A_xy(2, k) * DA_xy(3, i) +   &
      A_xy(1, k) * A_xy(2, i) * DA_xy(3, i))

  end do

!     Mpole_xy(20, 20) =
!             A_xy(1, 1) * A_xy(2, 2) * A_xy(3, 3) +
!             A_xy(1, 1) * A_xy(3, 2) * A_xy(2, 3) +
!             A_xy(2, 1) * A_xy(1, 2) * A_xy(3, 3) +
!             A_xy(2, 1) * A_xy(3, 2) * A_xy(1, 3) +
!             A_xy(3, 1) * A_xy(1, 2) * A_xy(2, 3) +
!             A_xy(3, 1) * A_xy(2, 2) * A_xy(1, 3)

  DMpole_xy(20, 20) = &
              DA_xy(1, 1) * A_xy(2, 2) * A_xy(3, 3) + &
              DA_xy(1, 1) * A_xy(3, 2) * A_xy(2, 3) + &
              DA_xy(2, 1) * A_xy(1, 2) * A_xy(3, 3) + &
              DA_xy(2, 1) * A_xy(3, 2) * A_xy(1, 3) + &
              DA_xy(3, 1) * A_xy(1, 2) * A_xy(2, 3) + &
              DA_xy(3, 1) * A_xy(2, 2) * A_xy(1, 3)

  DMpole_xy(20, 20) = DMpole_xy(20, 20) + &
              A_xy(1, 1) * DA_xy(2, 2) * A_xy(3, 3) + &
              A_xy(1, 1) * DA_xy(3, 2) * A_xy(2, 3) + &
              A_xy(2, 1) * DA_xy(1, 2) * A_xy(3, 3) + &
              A_xy(2, 1) * DA_xy(3, 2) * A_xy(1, 3) + &
              A_xy(3, 1) * DA_xy(1, 2) * A_xy(2, 3) + &
              A_xy(3, 1) * DA_xy(2, 2) * A_xy(1, 3)

  DMpole_xy(20, 20) = DMpole_xy(20, 20) + &
              A_xy(1, 1) * A_xy(2, 2) * DA_xy(3, 3) + &
              A_xy(1, 1) * A_xy(3, 2) * DA_xy(2, 3) + &
              A_xy(2, 1) * A_xy(1, 2) * DA_xy(3, 3) + &
              A_xy(2, 1) * A_xy(3, 2) * DA_xy(1, 3) + &
              A_xy(3, 1) * A_xy(1, 2) * DA_xy(2, 3) + &
              A_xy(3, 1) * A_xy(2, 2) * DA_xy(1, 3)

  if (order .eq. 20) return

  do j = 11, 20
    do k = 21, order
      DMpole_xy(k, j) = 0.d0
      DMpole_xy(j, k) = 0.d0
    end do
  end do

! HEXADECAPOLES
! H'_mmmm

  do ind1 = 1, 3
    m = hind1(ind1)
    do ind2 = 1, 15
      i = hind1(ind2)
      j = hind2(ind2)
      k = hind3(ind2)
      l = hind4(ind2)

!         Mpole_xy(ind1 + 20, ind2 + 20) =
!              A_xy(m, i) * A_xy(m, j) * A_xy(m, k) * A_xy(m, l)

      DMpole_xy(ind1 + 20, ind2 + 20) = &
               DA_xy(m, i) * A_xy(m, j) * A_xy(m, k) * A_xy(m, l)
      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) +  &
               A_xy(m, i) * DA_xy(m, j) * A_xy(m, k) * A_xy(m, l)
      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) +  &
               A_xy(m, i) * A_xy(m, j) * DA_xy(m, k) * A_xy(m, l)
      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) +  &
               A_xy(m, i) * A_xy(m, j) * A_xy(m, k) * DA_xy(m, l)
    end do
  end do

! H'_mmmn

  do ind1 = 4, 9
    m = hind1(ind1)
    n = hind4(ind1)
! can put iiii & iiij together
    do ind2 = 1, 9
      i = hind1(ind2)
      j = hind4(ind2)

!         Mpole_xy(ind1 + 20, ind2 + 20) =
!              A_xy(m, i) * A_xy(m, i) * A_xy(m, i) * A_xy(n, j) + 3.d0*
!              A_xy(m, i) * A_xy(m, i) * A_xy(n, i) * A_xy(m, j)

      DMpole_xy(ind1 + 20, ind2 + 20) = &
        DA_xy(m, i) * A_xy(m, i) * A_xy(m, i) * A_xy(n, j) + &
        3.d0 * DA_xy(m, i) * A_xy(m, i) * A_xy(n, i) * A_xy(m, j)

      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) + &
        A_xy(m, i) * DA_xy(m, i) * A_xy(m, i) * A_xy(n, j) + &
        3.d0 * A_xy(m, i) * DA_xy(m, i) * A_xy(n, i) * A_xy(m, j)

      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) + &
        A_xy(m, i) * A_xy(m, i) * DA_xy(m, i) * A_xy(n, j) + &
        3.d0 * A_xy(m, i) * A_xy(m, i) * DA_xy(n, i) * A_xy(m, j)

      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) + &
        A_xy(m, i) * A_xy(m, i) * A_xy(m, i) * DA_xy(n, j) + &
        3.d0 * A_xy(m, i) * A_xy(m, i) * A_xy(n, i) * DA_xy(m, j)

    end do

! can put iijj & iijk together

    do ind2 = 10, 15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)

!         Mpole_xy(ind1 + 20, ind2 + 20) =
!              A_xy(m, i) * A_xy(m, i) * A_xy(m, j) * A_xy(n, k) +
!              A_xy(m, i) * A_xy(m, i) * A_xy(m, k) * A_xy(n, j) + 2.d0*
!              A_xy(m, i) * A_xy(m, j) * A_xy(m, k) * A_xy(n, i)

      DMpole_xy(ind1 + 20, ind2 + 20) = &
        DA_xy(m, i) * A_xy(m, i) * A_xy(m, j) * A_xy(n, k) + &
        DA_xy(m, i) * A_xy(m, i) * A_xy(m, k) * A_xy(n, j) + &
        2.d0 * DA_xy(m, i) * A_xy(m, j) * A_xy(m, k) * A_xy(n, i)

      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) + &
        A_xy(m, i) * DA_xy(m, i) * A_xy(m, j) * A_xy(n, k) + &
        A_xy(m, i) * DA_xy(m, i) * A_xy(m, k) * A_xy(n, j) + &
        2.d0 * A_xy(m, i) * DA_xy(m, j) * A_xy(m, k) * A_xy(n, i)

      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) + &
        A_xy(m, i) * A_xy(m, i) * DA_xy(m, j) * A_xy(n, k) + &
        A_xy(m, i) * A_xy(m, i) * DA_xy(m, k) * A_xy(n, j) + &
        2.d0 * A_xy(m, i) * A_xy(m, j) * DA_xy(m, k) * A_xy(n, i)

      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) + &
        A_xy(m, i) * A_xy(m, i) * A_xy(m, j) * DA_xy(n, k) + &
        A_xy(m, i) * A_xy(m, i) * A_xy(m, k) * DA_xy(n, j) + &
        2.d0 * A_xy(m, i) * A_xy(m, j) * A_xy(m, k) * DA_xy(n, i)

    end do
  end do

! H'_mmnn

  do ind1 = 10, 12
    m = hind1(ind1)
    n = hind3(ind1)
! can put iiii & iiij together
    do ind2 = 1, 9
      i = hind1(ind2)
      j = hind4(ind2)

!        Mpole_xy(ind1 + 20, ind2 + 20) = 3.d0 * (
!              A_xy(m, i) * A_xy(m, i) * A_xy(n, i) * A_xy(n, j) +
!              A_xy(m, i) * A_xy(m, j) * A_xy(n, i) * A_xy(n, i))

      DMpole_xy(ind1 + 20, ind2 + 20) = &
        3.d0 * (DA_xy(m, i) * A_xy(m, i) * A_xy(n, i) * A_xy(n, j) + &
        DA_xy(m, i) * A_xy(m, j) * A_xy(n, i) * A_xy(n, i))

      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) + &
        3.d0 * (A_xy(m, i) * DA_xy(m, i) * A_xy(n, i) * A_xy(n, j) + &
        A_xy(m, i) * DA_xy(m, j) * A_xy(n, i) * A_xy(n, i))

      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) + &
        3.d0 * (A_xy(m, i) * A_xy(m, i) * DA_xy(n, i) * A_xy(n, j) + &
        A_xy(m, i) * A_xy(m, j) * DA_xy(n, i) * A_xy(n, i))

      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) + &
        3.d0 * (A_xy(m, i) * A_xy(m, i) * A_xy(n, i) * DA_xy(n, j) + &
        A_xy(m, i) * A_xy(m, j) * A_xy(n, i) * DA_xy(n, i))

    end do

! can put iijj & iijk together

    do ind2 = 10, 15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)

!         Mpole_xy(ind1 + 20, ind2 + 20) =
!              A_xy(m, i) * A_xy(m, i) * A_xy(n, j) * A_xy(n, k) +
!              A_xy(m, j) * A_xy(m, k) * A_xy(n, i) * A_xy(n, i) +  2.d0*
!              A_xy(m, i) * A_xy(m, j) * A_xy(n, i) * A_xy(n, k) +  2.d0*
!              A_xy(m, i) * A_xy(m, k) * A_xy(n, i) * A_xy(n, j)

      DMpole_xy(ind1 + 20, ind2 + 20) = &
        DA_xy(m, i) * A_xy(m, i) * A_xy(n, j) * A_xy(n, k) + &
        DA_xy(m, j) * A_xy(m, k) * A_xy(n, i) * A_xy(n, i) + &
        2.d0 * DA_xy(m, i) * A_xy(m, j) * A_xy(n, i) * A_xy(n, k) + &
        2.d0 * DA_xy(m, i) * A_xy(m, k) * A_xy(n, i) * A_xy(n, j)

      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) + &
        A_xy(m, i) * DA_xy(m, i) * A_xy(n, j) * A_xy(n, k) + &
        A_xy(m, j) * DA_xy(m, k) * A_xy(n, i) * A_xy(n, i) + &
        2.d0 * A_xy(m, i) * DA_xy(m, j) * A_xy(n, i) * A_xy(n, k) + &
        2.d0 * A_xy(m, i) * DA_xy(m, k) * A_xy(n, i) * A_xy(n, j)

      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) + &
        A_xy(m, i) * A_xy(m, i) * DA_xy(n, j) * A_xy(n, k) + &
        A_xy(m, j) * A_xy(m, k) * DA_xy(n, i) * A_xy(n, i) + &
        2.d0 * A_xy(m, i) * A_xy(m, j) * DA_xy(n, i) * A_xy(n, k) + &
        2.d0 * A_xy(m, i) * A_xy(m, k) * DA_xy(n, i) * A_xy(n, j)

      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) + &
        A_xy(m, i) * A_xy(m, i) * A_xy(n, j) * DA_xy(n, k) + &
        A_xy(m, j) * A_xy(m, k) * A_xy(n, i) * DA_xy(n, i) + &
        2.d0 * A_xy(m, i) * A_xy(m, j) * A_xy(n, i) * DA_xy(n, k) + &
        2.d0 * A_xy(m, i) * A_xy(m, k) * A_xy(n, i) * DA_xy(n, j)

    end do
  end do

! H'_mmnp

  do ind1 = 13, 15
    m = hind1(ind1)
    n = hind3(ind1)
    p = hind4(ind1)
! can put iiii & iiij together
    do ind2 = 1, 9
      i = hind1(ind2)
      j = hind4(ind2)

!         Mpole_xy(ind1 + 20, ind2 + 20) =
!             3.d0 * A_xy(m, i) * A_xy(m, i) * A_xy(n, i) * A_xy(p, j) +
!             3.d0 * A_xy(m, i) * A_xy(m, i) * A_xy(n, j) * A_xy(p, i) +
!             6.d0 * A_xy(m, i) * A_xy(m, j) * A_xy(n, i) * A_xy(p, i)

      DMpole_xy(ind1 + 20, ind2 + 20) = &
              3.d0 * DA_xy(m, i) * A_xy(m, i) * A_xy(n, i) * A_xy(p, j) + &
              3.d0 * DA_xy(m, i) * A_xy(m, i) * A_xy(n, j) * A_xy(p, i) + &
              6.d0 * DA_xy(m, i) * A_xy(m, j) * A_xy(n, i) * A_xy(p, i)

      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) + &
              3.d0 * A_xy(m, i) * DA_xy(m, i) * A_xy(n, i) * A_xy(p, j) + &
              3.d0 * A_xy(m, i) * DA_xy(m, i) * A_xy(n, j) * A_xy(p, i) + &
              6.d0 * A_xy(m, i) * DA_xy(m, j) * A_xy(n, i) * A_xy(p, i)

      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) + &
              3.d0 * A_xy(m, i) * A_xy(m, i) * DA_xy(n, i) * A_xy(p, j) + &
              3.d0 * A_xy(m, i) * A_xy(m, i) * DA_xy(n, j) * A_xy(p, i) + &
              6.d0 * A_xy(m, i) * A_xy(m, j) * DA_xy(n, i) * A_xy(p, i)

      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) + &
              3.d0 * A_xy(m, i) * A_xy(m, i) * A_xy(n, i) * DA_xy(p, j) + &
              3.d0 * A_xy(m, i) * A_xy(m, i) * A_xy(n, j) * DA_xy(p, i) + &
              6.d0 * A_xy(m, i) * A_xy(m, j) * A_xy(n, i) * DA_xy(p, i)

    end do

! can put iijj & iijk together

    do ind2 = 10, 15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
!         Mpole_xy(ind1 + 20, ind2 + 20) =
!              A_xy(m, i) * A_xy(m, i) * A_xy(n, j) * A_xy(p, k) +
!              A_xy(m, i) * A_xy(m, i) * A_xy(n, k) * A_xy(p, j) + 2.d0 * (
!              A_xy(m, i) * A_xy(m, j) * A_xy(n, i) * A_xy(p, k) +
!              A_xy(m, i) * A_xy(m, j) * A_xy(n, k) * A_xy(p, i) +
!              A_xy(m, i) * A_xy(m, k) * A_xy(n, i) * A_xy(p, j) +
!              A_xy(m, i) * A_xy(m, k) * A_xy(n, j) * A_xy(p, i) +
!              A_xy(m, j) * A_xy(m, k) * A_xy(n, i) * A_xy(p, i))

      DMpole_xy(ind1 + 20, ind2 + 20) = &
        DA_xy(m, i) * A_xy(m, i) * A_xy(n, j) * A_xy(p, k) + &
        DA_xy(m, i) * A_xy(m, i) * A_xy(n, k) * A_xy(p, j) + &
        2.d0 * (DA_xy(m, i) * A_xy(m, j) * A_xy(n, i) * A_xy(p, k) + &
        DA_xy(m, i) * A_xy(m, j) * A_xy(n, k) * A_xy(p, i) + &
        DA_xy(m, i) * A_xy(m, k) * A_xy(n, i) * A_xy(p, j) + &
        DA_xy(m, i) * A_xy(m, k) * A_xy(n, j) * A_xy(p, i) + &
        DA_xy(m, j) * A_xy(m, k) * A_xy(n, i) * A_xy(p, i))

      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) + &
        A_xy(m, i) * DA_xy(m, i) * A_xy(n, j) * A_xy(p, k) + &
        A_xy(m, i) * DA_xy(m, i) * A_xy(n, k) * A_xy(p, j) + &
        2.d0 * (A_xy(m, i) * DA_xy(m, j) * A_xy(n, i) * A_xy(p, k) + &
        A_xy(m, i) * DA_xy(m, j) * A_xy(n, k) * A_xy(p, i) + &
        A_xy(m, i) * DA_xy(m, k) * A_xy(n, i) * A_xy(p, j) + &
        A_xy(m, i) * DA_xy(m, k) * A_xy(n, j) * A_xy(p, i) + &
        A_xy(m, j) * DA_xy(m, k) * A_xy(n, i) * A_xy(p, i))

      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) + &
        A_xy(m, i) * A_xy(m, i) * DA_xy(n, j) * A_xy(p, k) + &
        A_xy(m, i) * A_xy(m, i) * DA_xy(n, k) * A_xy(p, j) + &
        2.d0 * (A_xy(m, i) * A_xy(m, j) * DA_xy(n, i) * A_xy(p, k) + &
        A_xy(m, i) * A_xy(m, j) * DA_xy(n, k) * A_xy(p, i) + &
        A_xy(m, i) * A_xy(m, k) * DA_xy(n, i) * A_xy(p, j) + &
        A_xy(m, i) * A_xy(m, k) * DA_xy(n, j) * A_xy(p, i) + &
        A_xy(m, j) * A_xy(m, k) * DA_xy(n, i) * A_xy(p, i))

      DMpole_xy(ind1 + 20, ind2 + 20) = DMpole_xy(ind1 + 20, ind2 + 20) + &
        A_xy(m, i) * A_xy(m, i) * A_xy(n, j) * DA_xy(p, k) + &
        A_xy(m, i) * A_xy(m, i) * A_xy(n, k) * DA_xy(p, j) + &
        2.d0 * (A_xy(m, i) * A_xy(m, j) * A_xy(n, i) * DA_xy(p, k) + &
        A_xy(m, i) * A_xy(m, j) * A_xy(n, k) * DA_xy(p, i) + &
        A_xy(m, i) * A_xy(m, k) * A_xy(n, i) * DA_xy(p, j) + &
        A_xy(m, i) * A_xy(m, k) * A_xy(n, j) * DA_xy(p, i) + &
        A_xy(m, j) * A_xy(m, k) * A_xy(n, i) * DA_xy(p, i))

    end do
  end do

  return

end subroutine xform_mpole_deriv_matrix

!*******************************************************************************!
! Subroutine:  xform_mpole_field_matrix
!
! Description:
!
! A_xy is matrix of first order terms connecting field in y-coord frame
! to those in the x-frame
! Field_xy expands this to get the higher order conversions between coord sys
! order is order of multipoles 1 for charge, 4 for charge-dipoles
! 10 for charge, dipole, quadrupole, 20 for up to octupole and
! finally 35 for up to hexadecapole
!
! Field terms in y-frame are dphi_dy1, dphi_dy2, dphi_dy3
!       those in x-frame are dphi_dx1, dphi_dx2, dphi_dx3
! Thus:  A_xy(1, 1) = dy1_dx1, A_xy(1, 2) = dy2_dx1, A_xy(1, 3) = dy3_dx1
!        A_xy(2, 1) = dy1_dx2, A_xy(2, 2) = dy2_dx2, A_xy(2, 3) = dy3_dx2
!        A_xy(3, 1) = dy1_dx3, A_xy(3, 2) = dy2_dx3, A_xy(3, 3) = dy3_dx3
! For example to get 2nd order derivs d^2_dx1_dx2 we expand
!  (d_dy1 * dy1_dx1 + d_dy2 * dy2_dx1 + d_dy3 * dy3_dx1)*
!  (d_dy1 * dy1_dx2 + d_dy2 * dy2_dx2 + d_dy3 * dy3_dx2)
! and collect terms involving each d^2_dyk_dyl
!
!*******************************************************************************

subroutine xform_mpole_field_matrix(A_xy, Field_xy, order)

  implicit none

! Formal arguments:

  integer               :: order
  double precision      :: A_xy(3, 3)
  double precision      :: Field_xy(order, order)

! Local variables:

  integer               :: i, j, k, l, m, n, p, q, ind1, ind2, jj, kk
  integer               :: qind1(6), qind2(6)
  integer               :: oind1(10), oind2(10), oind3(10)
  integer               :: hind1(15), hind2(15), hind3(15), hind4(15)

  common /MxfInd/ qind1, qind2, oind1, oind2, oind3, hind1, hind2, hind3, hind4

! CHARGE case

  Field_xy(1, 1) = 1.d0

  if (order .eq. 1) return

  do j = 2, order
    Field_xy(1, j) = 0.d0
    Field_xy(j, 1) = 0.d0
  end do

! DIPOLES
! d_yk = A_xy(k, 1) * d_x1 + A_xy(k, 2) * d_x2 + A_xy(k, 3) * d_x3
! D'_j

  do j = 2, 4
    do k = 2, 4
      Field_xy(j, k) = A_xy(j-1, k-1)
    end do
  end do

  if (order .eq. 4) return

  do j = 2, 4
    do k = 5, order
      Field_xy(j, k) = 0.d0
      Field_xy(k, j) = 0.d0
    end do
  end do

! QUADRUPOLES
! q_ykyl = sum_i, j (A_xy(k, i) * A_xy(l, j) + A_xy(k, j) * A_xy(l, i)) * (q_xixj + q_xjxi)
! Mp(5) = q_y1y1, .., Mp(7) = q_y3y3, Mp(8) = q_y1y2 + q_y2y1, ., Mp(10) = q_y2y3 + q_y3y2
! Q'_kk

  do ind1 = 1, 3
    l = qind1(ind1)
    do ind2 = 1, 3
      i = qind1(ind2)
      Field_xy(ind1 + 4, ind2 + 4) = A_xy(l, i) * A_xy(l, i)
    end do
    do ind2 = 4, 6
      i = qind1(ind2)
      j = qind2(ind2)
      Field_xy(ind1 + 4, ind2 + 4) = 2.d0 * A_xy(l, i) * A_xy(l, j)
    end do
  end do
  do ind1 = 4, 6
    l = qind1(ind1)
    m = qind2(ind1)
    do ind2 = 1, 3
      i = qind1(ind2)
      Field_xy(ind1 + 4, ind2 + 4) = A_xy(l, i) * A_xy(m, i)
    end do
    do ind2 = 4, 6
      i = qind1(ind2)
      j = qind2(ind2)
      Field_xy(ind1 + 4, ind2 + 4) = (A_xy(l, i) * A_xy(m, j) + &
                                      A_xy(m, i) * A_xy(l, j))
    end do
  end do

  if (order .eq. 10) return

  do j = 5, 10
    do k = 11, order
      Field_xy(k, j) = 0.d0
      Field_xy(j, k) = 0.d0
    end do
  end do

! OCTUPOLES level field
! F'_lll

  do ind1 = 1, 3
    l = oind1(ind1)
    do ind2 = 1, 3
      i = oind1(ind2)
      Field_xy(ind1 + 10, ind2 + 10) = A_xy(l, i)**3
    end do
    do ind2 = 4,9
      i = oind1(ind2)
      j = oind3(ind2)
      Field_xy(ind1 + 10, ind2 + 10) = 3.d0 * A_xy(l, i)**2 * A_xy(l, j)
    end do
    Field_xy(ind1 + 10, 20) = 6.d0 * A_xy(l, 1) * A_xy(l, 2) * A_xy(l, 3) 
  end do

! F'_llm

  do ind1 = 4, 9
    l = oind1(ind1)
    m = oind3(ind1)
    do ind2 = 1, 3
      i = oind1(ind2)
      Field_xy(ind1 + 10, ind2 + 10) = A_xy(l, i)**2 * A_xy(m, i)
    end do
    do ind2 = 4, 9
      i = oind1(ind2)
      j = oind3(ind2)
      Field_xy(ind1 + 10, ind2 + 10) = A_xy(l, i)**2 * A_xy(m, j) +   &
                                2.d0 * A_xy(l, i) * A_xy(l, j) * A_xy(m, i)
    end do

    Field_xy(ind1 + 10, 20) = 2.d0 * (A_xy(l, 1) * A_xy(l, 2) * A_xy(m, 3) + &
                                      A_xy(l, 1) * A_xy(m, 2) * A_xy(l, 3) + &
                                      A_xy(m, 1) * A_xy(l, 2) * A_xy(l, 3))
  end do

! F'_123

  do ind2 = 1, 3
    i = oind1(ind2)
    Field_xy(20, ind2 + 10) = A_xy(1, i) * A_xy(2, i) * A_xy(3, i)
  end do
  do ind2 = 4, 9
    i = oind1(ind2)
    j = oind3(ind2)
    Field_xy(20, ind2 + 10) = A_xy(1, i) * A_xy(2, i) * A_xy(3, j) +  &
                              A_xy(1, i) * A_xy(2, j) * A_xy(3, i) +  &
                              A_xy(1, j) * A_xy(2, i) * A_xy(3, i) 
  end do

  Field_xy(20, 20) = A_xy(1, 1) * A_xy(2, 2) * A_xy(3, 3) + &
                     A_xy(1, 1) * A_xy(2, 3) * A_xy(3, 2) + &
                     A_xy(1, 2) * A_xy(2, 1) * A_xy(3, 3) + &
                     A_xy(1, 2) * A_xy(2, 3) * A_xy(3, 1) + &
                     A_xy(1, 3) * A_xy(2, 1) * A_xy(3, 2) + &
                     A_xy(1, 3) * A_xy(2, 2) * A_xy(3, 1)

  if (order .eq. 20) return

  do j = 11, 20
    do k = 21, order
      Field_xy(k, j) = 0.d0
      Field_xy(j, k) = 0.d0
    end do
  end do

! HEXADECAPOLES
! F'_mmmm

  do ind1 = 1, 3
    m = hind1(ind1)
    do ind2 = 1, 3
      i = hind1(ind2)
      Field_xy(ind1 + 20, ind2 + 20) = A_xy(m, i)**4
    end do
    do ind2 = 4, 9
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1 + 20, ind2 + 20) = 4.d0 * A_xy(m, i)**3 * A_xy(m, j)
    end do
    do ind2 = 10, 12
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1 + 20, ind2 + 20) = 6.d0 * A_xy(m, i)**2 * A_xy(m, j)**2
    end do
    do ind2 = 13, 15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Field_xy(ind1 + 20, ind2 + 20) = 12.d0 * A_xy(m, i)**2 * &
                                       A_xy(m, j) * A_xy(m, k)
    end do
  end do
  do ind1 = 4, 9
    m = hind1(ind1)
    n = hind4(ind1)
    do ind2 = 1, 3
      i = hind1(ind2)
      Field_xy(ind1 + 20, ind2 + 20) = A_xy(m, i)**3 * A_xy(n, i)
    end do
    do ind2 = 4, 9
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1 + 20, ind2 + 20) = A_xy(m, i)**3 * A_xy(n, j) +  &
                                       3.d0 * A_xy(m, i)**2 * &
                                       A_xy(m, j) * A_xy(n, i)
    end do
    do ind2 = 10, 12
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1 + 20, ind2 + 20) =    &
              3.d0 * A_xy(m, i)**2 * A_xy(m, j) * A_xy(n, j) +    &
              3.d0 * A_xy(m, j)**2 * A_xy(m, i) * A_xy(n, i) 
    end do
    do ind2 = 13, 15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Field_xy(ind1 + 20, ind2 + 20) = &
        6.d0 * A_xy(m, i) * A_xy(m, j) * A_xy(m, k) * A_xy(n, i) + &
        3.d0 * A_xy(m, i)**2 * A_xy(m, j) * A_xy(n, k) + &
        3.d0 * A_xy(m, i)**2 * A_xy(m, k) * A_xy(n, j)
    end do
  end do
  do ind1 = 10, 12
    m = hind1(ind1)
    n = hind4(ind1)
    do ind2 = 1, 3
      i = hind1(ind2)
      Field_xy(ind1 + 20, ind2 + 20) = A_xy(m, i)**2 * A_xy(n, i)**2
    end do
    do ind2 = 4, 9
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1 + 20, ind2 + 20) =    &
        2.d0 * A_xy(m, i)**2 * A_xy(n, i) * A_xy(n, j) +   &
        2.d0 * A_xy(n, i)**2 * A_xy(m, i) * A_xy(m, j)
    end do
    do ind2 = 10, 12
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1 + 20, ind2 + 20) = &
        A_xy(m, i)**2 * A_xy(n, j)**2 + &
        A_xy(m, j)**2 * A_xy(n, i)**2 + &
        4.d0 *  A_xy(m, i) * A_xy(m, j) * A_xy(n, i) * A_xy(n, j)
    end do
    do ind2 = 13, 15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Field_xy(ind1 + 20, ind2 + 20) = &
        2.d0 * A_xy(m, i)**2 * A_xy(n, j) * A_xy(n, k) + &
        2.d0 * A_xy(n, i)**2 * A_xy(m, j) * A_xy(m, k) + &
        4.d0 * A_xy(m, i) * A_xy(m, j) * A_xy(n, i) * A_xy(n, k) + &
        4.d0 * A_xy(m, i) * A_xy(m, k) * A_xy(n, i) * A_xy(n, j) 
    end do
  end do
  do ind1 = 13, 15
    m = hind1(ind1)
    n = hind3(ind1)
    p = hind4(ind1)
    do ind2 = 1, 3
      i = hind1(ind2)
      Field_xy(ind1 + 20, ind2 + 20) = A_xy(m, i)**2 * A_xy(n, i) * A_xy(p, i)
    end do
    do ind2 = 4, 9
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1 + 20, ind2 + 20) = &
        2.d0 * A_xy(m, i) * A_xy(m, j) * A_xy(n, i) * A_xy(p, i) + &
        A_xy(m, i)**2 * (A_xy(n, i) * A_xy(p, j) + A_xy(n, j) * A_xy(p, i))
    end do
    do ind2 = 10, 12
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1 + 20, ind2 + 20) = &
        A_xy(m, i)**2 * A_xy(n, j) * A_xy(p, j) + &
        A_xy(m, j)**2 * A_xy(n, i) * A_xy(p, i) + &
        2.d0 * A_xy(m, i) * A_xy(m, j) * &
        (A_xy(n, i) * A_xy(p, j) + A_xy(n, j) * A_xy(p, i))
    end do
    do ind2 = 13, 15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Field_xy(ind1 + 20, ind2 + 20) = &
        A_xy(m, i)**2 * (A_xy(n, j) * A_xy(p, k) + A_xy(n, k) * A_xy(p, j)) + &
        2.d0 * A_xy(m, i) * A_xy(m, j) * (A_xy(n, i) * A_xy(p, k) + &
        A_xy(n, k) * A_xy(p, i)) + &
        2.d0 * A_xy(m, i) * A_xy(m, k) * (A_xy(n, i) * A_xy(p, j) + &
        A_xy(n, j) * A_xy(p, i)) + &
        2.d0 * A_xy(m, j) * A_xy(m, k) * A_xy(n, i) * A_xy(p, i)
    end do
  end do

  return

end subroutine xform_mpole_field_matrix

!*******************************************************************************!
! Subroutine:  xform_mpole
!
! Description: <TBS>
!
!*******************************************************************************

subroutine xform_mpole(Mpole_xy, dimxy, Mpole_in, Mpole_out, order)

  implicit none

! Formal arguments:

  integer               :: dimxy
  double precision      :: Mpole_xy(dimxy, dimxy)
  double precision      :: Mpole_in(*)
  double precision      :: Mpole_out(*)
  integer               :: order

! Local variables:

  integer               :: i, j

  if (order .eq. 0) return

  Mpole_out(1) = Mpole_xy(1, 1) * Mpole_in(1)

  if (order .eq. 1) return

! DIPOLES

  do i = 2, 4
    Mpole_out(i) = 0.d0
    do j = 2, 4
      Mpole_out(i) = Mpole_out(i) + Mpole_xy(i, j) * Mpole_in(j)
    end do
  end do

  if (order .eq. 4) return

! QUADRUPOLES

  do i = 5, 10
    Mpole_out(i) = 0.d0
    do j = 5, 10
      Mpole_out(i) = Mpole_out(i) + Mpole_xy(i, j) * Mpole_in(j)
    end do
  end do

  if (order .eq. 10) return

! OCTUPOLES

  do i = 11, 20
    Mpole_out(i) = 0.d0
    do j = 11, 20
      Mpole_out(i) = Mpole_out(i) + Mpole_xy(i, j) * Mpole_in(j)
    end do
  end do

  if (order .eq. 20) return

! HEXADECAPOLES

  do i = 21, 35
    Mpole_out(i) = 0.d0
    do j = 21, 35
      Mpole_out(i) = Mpole_out(i) + Mpole_xy(i, j) * Mpole_in(j)
    end do
  end do

  if (order .eq. 35) return

  return
  
end subroutine xform_mpole

!*******************************************************************************!
! Subroutine:   xform_field
!
! Description: <TBS>
!
!*******************************************************************************

subroutine xform_field(Field_xy, dimxy, Field_in, Field_out, order)

  implicit none

! Formal arguments:

  integer               :: dimxy
  double precision      :: Field_xy(dimxy, dimxy)
  double precision      :: Field_in(*)
  double precision      :: Field_out(*)
  integer               :: order

! Local variables:

  integer               :: i, j

  ! Field_out is xform of Field_in

  if (order .eq. 0) return

  Field_out(1) = Field_xy(1, 1) * Field_in(1)

  if (order .eq. 1) return

! DIPOLES

  do i = 2, 4
    Field_out(i) = 0.d0
    do j = 2, 4
      Field_out(i) = Field_out(i) + Field_xy(i, j) * Field_in(j)
    end do
  end do

  if (order .eq. 4) return

! QUADRUPOLES

  do i = 5, 10
    Field_out(i) = 0.d0
    do j = 5, 10
      Field_out(i) = Field_out(i) + Field_xy(i, j) * Field_in(j)
    end do
  end do

  if (order .eq. 10) return

! OCTUPOLES

  do i = 11, 20
    Field_out(i) = 0.d0
    do j = 11, 20
      Field_out(i) = Field_out(i) + Field_xy(i, j) * Field_in(j)
    end do
  end do

  if (order .eq. 20) return

! HEXADECAPOLES

  do i = 21, 35
    Field_out(i) = 0.d0
    do j = 21, 35
      Field_out(i) = Field_out(i) + Field_xy(i, j) * Field_in(j)
    end do
  end do

  if (order .eq. 35) return

  return

end subroutine xform_field
