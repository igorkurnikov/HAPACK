#include "copyright.i"

!*******************************************************************************
!
! Module:  dihedrals_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module dihedrals_mod

  use gbl_datatypes_mod

  implicit none

! The following are derived from prmtop dihedral info:

  integer, save                         :: cit_nphih, cit_nphia

  type(dihed_rec), allocatable, save    :: cit_dihed(:)

contains

!*******************************************************************************
!
! Subroutine:  dihedrals_setup
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine dihedrals_setup(use_atm_map, atm_owner_map)

  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer                       :: use_atm_map(natom)
  integer                       :: atm_owner_map(natom)

! Local variables:

  integer               :: alloc_failed
  type(dihed_rec)       :: diheds_copy(nphih + nphia)
  integer               :: atm_i, atm_j, atm_k, atm_l, diheds_idx
  integer               :: my_dihed_cnt

! This routine can handle reallocation, and thus can be called multiple
! times.

! Find all diheds for which this process owns either atom:

  my_dihed_cnt = 0

  do diheds_idx = 1, nphih

    atm_i = gbl_dihed(diheds_idx)%atm_i
    atm_j = gbl_dihed(diheds_idx)%atm_j
    atm_k = iabs(gbl_dihed(diheds_idx)%atm_k)
    atm_l = iabs(gbl_dihed(diheds_idx)%atm_l)

    if (atm_owner_map(atm_i) .eq. mytaskid) then
      my_dihed_cnt = my_dihed_cnt + 1
      diheds_copy(my_dihed_cnt) = gbl_dihed(diheds_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
    else if (atm_owner_map(atm_j) .eq. mytaskid) then
      my_dihed_cnt = my_dihed_cnt + 1
      diheds_copy(my_dihed_cnt) = gbl_dihed(diheds_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
    else if (atm_owner_map(atm_k) .eq. mytaskid) then
      my_dihed_cnt = my_dihed_cnt + 1
      diheds_copy(my_dihed_cnt) = gbl_dihed(diheds_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
    else if (atm_owner_map(atm_l) .eq. mytaskid) then
      my_dihed_cnt = my_dihed_cnt + 1
      diheds_copy(my_dihed_cnt) = gbl_dihed(diheds_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
    end if

  end do

  cit_nphih = my_dihed_cnt

  do diheds_idx = diheda_idx, diheda_idx + nphia - 1

    atm_i = gbl_dihed(diheds_idx)%atm_i
    atm_j = gbl_dihed(diheds_idx)%atm_j
    atm_k = iabs(gbl_dihed(diheds_idx)%atm_k)
    atm_l = iabs(gbl_dihed(diheds_idx)%atm_l)

    if (atm_owner_map(atm_i) .eq. mytaskid) then
      my_dihed_cnt = my_dihed_cnt + 1
      diheds_copy(my_dihed_cnt) = gbl_dihed(diheds_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
    else if (atm_owner_map(atm_j) .eq. mytaskid) then
      my_dihed_cnt = my_dihed_cnt + 1
      diheds_copy(my_dihed_cnt) = gbl_dihed(diheds_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
    else if (atm_owner_map(atm_k) .eq. mytaskid) then
      my_dihed_cnt = my_dihed_cnt + 1
      diheds_copy(my_dihed_cnt) = gbl_dihed(diheds_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
    else if (atm_owner_map(atm_l) .eq. mytaskid) then
      my_dihed_cnt = my_dihed_cnt + 1
      diheds_copy(my_dihed_cnt) = gbl_dihed(diheds_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
    end if

  end do

  cit_nphia = my_dihed_cnt - cit_nphih

  if (my_dihed_cnt .gt. 0) then
    if (allocated(cit_dihed)) then
      if (size(cit_dihed) .lt. my_dihed_cnt) then
        deallocate(cit_dihed)
        allocate(cit_dihed(my_dihed_cnt), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
      end if
    else
      allocate(cit_dihed(my_dihed_cnt), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
    end if
    cit_dihed(1:my_dihed_cnt) = diheds_copy(1:my_dihed_cnt)
  end if

!BEGIN DBG
! write(0,*)'task,cit_nphih,cit_nphia', mytaskid, &
!            cit_nphih, cit_nphia,
! write(0,*)'task,    nphih,    nphia', mytaskid, &
!                nphih,     nphia
!END DBG

  return

end subroutine dihedrals_setup

!*******************************************************************************
!
! Subroutine:  get_dihed_energy
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine get_dihed_energy(atm_owner_map, dihed_cnt, dihed, x, frc, ep, enbp, eelp, &
                            rel_crd, molvir, e14_vir )

  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer                       :: atm_owner_map(*)
  integer                       :: dihed_cnt
  type(dihed_rec)               :: dihed(*)
  double precision              :: x(3, *)
  double precision              :: frc(3, *)
  double precision              :: ep
  double precision              :: enbp
  double precision              :: eelp

  ! These can be left off for GB calls:

  double precision, optional    :: rel_crd(3, *)
  double precision, optional    :: molvir(3, 3)
  double precision, optional    :: e14_vir(3, 3)! e14_vir - must zero externally
                                                ! because routine called 2x.
                                                
              
! Local variables:

  double precision      :: ap
  double precision      :: cosnp
  double precision      :: cphi
  double precision      :: ct, ct0
  double precision      :: dc1, dc2, dc3, dc4, dc5, dc6
  double precision      :: df
  double precision      :: dfn
  double precision      :: dr1, dr2, dr3, dr4, dr5, dr6
  double precision      :: drx, dry, drz
  double precision      :: dums
  double precision      :: dx, dy, dz
  double precision      :: eelpl, enbpl, epl
  double precision      :: epw
  double precision      :: f1, f2
  double precision      :: fmuln
  double precision      :: fxi, fyi, fzi
  double precision      :: fxj, fyj, fzj
  double precision      :: fxk, fyk, fzk
  double precision      :: fxl, fyl, fzl
  double precision      :: g
  double precision      :: gmul(10)
  double precision      :: gx, gy, gz
  integer               :: ia1, ia2
  integer               :: ibig
  integer               :: ic, ic0
  integer               :: inc
  integer               :: i, j, k, kt, l, lt
  integer               :: jn
  double precision      :: one
  double precision      :: r1, r2, r6, r12
  double precision      :: s
  double precision      :: scnb0
  double precision      :: scee0
  double precision      :: sinnp
  double precision      :: six
  double precision      :: sphi
  double precision      :: tenm3
  double precision      :: tm06, tm24
  double precision      :: twelve
  double precision      :: xa, ya, za
  double precision      :: xij, yij, zij
  double precision      :: xkj, ykj, zkj
  double precision      :: xkl, ykl, zkl
  double precision      :: z1, z2
  double precision      :: z11, z12, z22
  double precision      :: zero
  double precision      :: vxx, vxy, vxz, vyy, vyz, vzz

  data gmul /0.d0, 2.d0, 0.d0, 4.d0, 0.d0, 6.d0, 0.d0, 8.d0, 0.d0, 10.d0/

  data tm24, tm06, tenm3/1.d-18, 1.d-06, 1.d-03/

  data zero, one, six, twelve /0.d0, 1.d0, 6.d0, 12.d0/

! e14 virials accumulators:

  vxx = zero
  vxy = zero
  vxz = zero
  vyy = zero
  vyz = zero
  vzz = zero

! Arrays gbl_gamc = gbl_pk * cos(phase) and gbl_gams = gbl_pk * sin(phase)

  epl = zero
  enbpl = zero
  eelpl = zero
  scnb0 = one / scnb
  scee0 = one / scee

! Grand loop for the dihedral stuff:

  do jn = 1, dihed_cnt

    i = dihed(jn)%atm_i
    j = dihed(jn)%atm_j
    kt = dihed(jn)%atm_k
    lt = dihed(jn)%atm_l
    k = iabs(kt)
    l = iabs(lt)

! Calculation of ij, kj, kl vectors:

    xij = x(1, i) - x(1, j)
    yij = x(2, i) - x(2, j)
    zij = x(3, i) - x(3, j)
    xkj = x(1, k) - x(1, j)
    ykj = x(2, k) - x(2, j)
    zkj = x(3, k) - x(3, j)
    xkl = x(1, k) - x(1, l)
    ykl = x(2, k) - x(2, l)
    zkl = x(3, k) - x(3, l)                                   

! Get the normal vector:

    dx = yij * zkj - zij * ykj
    dy = zij * xkj - xij * zkj
    dz = xij * ykj - yij * xkj
    gx = zkj * ykl - ykj * zkl
    gy = xkj * zkl - zkj * xkl
    gz = ykj * xkl - xkj * ykl

    fxi = sqrt(dx * dx + dy * dy + dz * dz + tm24)
    fyi = sqrt(gx * gx + gy * gy + gz * gz + tm24)
    ct = dx * gx + dy * gy + dz * gz

! Branch if linear dihedral:

    if (tenm3 .le. fxi) then
      z1 = one / fxi
    else
      z1 = zero
    endif

    if (tenm3 .le. fyi) then
      z2 = one / fyi
    else
      z2 = zero
    endif

    z12 = z1 * z2

    if (z12 .ne. zero) then
      fzi = one
    else
      fzi = zero
    end if

    s = xkj * (dz * gy - dy * gz) + ykj * (dx * gz - dz * gx) + &
        zkj * (dy * gx - dx * gy)

    ap = PI - sign(acos(max(-one, min(one, ct * z12))), s)

    cphi = cos(ap)
    sphi = sin(ap)

! Calculate the energy and the derivatives with respect to cosphi:

    ic = dihed(jn)%parm_idx
    inc = gbl_ipn(ic)
    ct0 = gbl_pn(ic) * ap
    cosnp = cos(ct0)
    sinnp = sin(ct0)
    epw = (gbl_pk(ic) + cosnp * gbl_gamc(ic) + sinnp * gbl_gams(ic)) * fzi
    
    dums = sphi + sign(tm24, sphi)

    if (tm06 .gt. abs(dums)) then
      df = fzi * gbl_gamc(ic) * (gbl_pn(ic) - gmul(inc) + gmul(inc) * cphi)
    else
      df = fzi * gbl_pn(ic) * (gbl_gamc(ic) * sinnp - gbl_gams(ic) * cosnp)/dums
    end if

! Now do torsional first derivatives:

! Now, set up array dc = 1st der. of cosphi w/respect to cartesian differences:

    z11 = z1 * z1
    z12 = z1 * z2
    z22 = z2 * z2
    dc1 = -gx * z12 - cphi * dx * z11
    dc2 = -gy * z12 - cphi * dy * z11
    dc3 = -gz * z12 - cphi * dz * z11
    dc4 =  dx * z12 + cphi * gx * z22
    dc5 =  dy * z12 + cphi * gy * z22
    dc6 =  dz * z12 + cphi * gz * z22

! Update the first derivative array:

    dr1 = df * ( dc3 * ykj - dc2 * zkj)
    dr2 = df * ( dc1 * zkj - dc3 * xkj)
    dr3 = df * ( dc2 * xkj - dc1 * ykj)
    dr4 = df * ( dc6 * ykj - dc5 * zkj)
    dr5 = df * ( dc4 * zkj - dc6 * xkj)
    dr6 = df * ( dc5 * xkj - dc4 * ykj)
    drx = df * (-dc2 * zij + dc3 * yij + dc5 * zkl - dc6 * ykl)
    dry = df * ( dc1 * zij - dc3 * xij - dc4 * zkl + dc6 * xkl)
    drz = df * (-dc1 * yij + dc2 * xij + dc4 * ykl - dc5 * xkl)
    fxi = -dr1
    fyi = -dr2
    fzi = -dr3
    fxj = -drx + dr1
    fyj = -dry + dr2
    fzj = -drz + dr3
    fxk = +drx + dr4
    fyk = +dry + dr5
    fzk = +drz + dr6
    fxl = -dr4
    fyl = -dr5
    fzl = -dr6

! End of a strip of dihedrals and start of 1-4 nb:

    ! NOTE that xij, yij, zij, ct are reassigned!

    xij = x(1, i) - x(1, l)
    yij = x(2, i) - x(2, l)
    zij = x(3, i) - x(3, l)

! Now loop over all the dihedrals (constant diel):

    ct = xij * xij + yij * yij + zij * zij

    fmuln = dble((2 + isign(1, kt) + isign(1, lt)) / 4) * gbl_fmn(ic)

    ia1 = atm_iac(i)
    ia2 = atm_iac(l)
    ibig = max0(ia1, ia2)
    ic0 = ibig * (ibig - 1) / 2 + min0(ia1, ia2)

! Calculate the 14-eel energy:

    r2 = fmuln / ct
    r1 = sqrt(r2)
    g = atm_qterm(i) * atm_qterm(l) * r1
    sphi = g
    r6 = r2 * r2 * r2
    r12 = r6 * r6
    f1 = gbl_cn1(ic0) * r12
    f2 = gbl_cn2(ic0) * r6
    cphi = f1 - f2
    
    dfn = ((-twelve * f1 + six * f2) * scnb0 - g * scee0) * r2

    xa = xij * dfn
    ya = yij * dfn 
    za = zij * dfn

    fxi = fxi - xa
    fyi = fyi - ya
    fzi = fzi - za
    fxl = fxl + xa
    fyl = fyl + ya
    fzl = fzl + za

! Sum up all the gradients:

! We use atm_i to determine who sums up the energy...
    if (atm_owner_map(i) .eq. mytaskid) then

      enbpl = enbpl + cphi
      eelpl = eelpl + sphi
      epl   = epl  + epw

      frc(1, i) = frc(1, i) + fxi
      frc(2, i) = frc(2, i) + fyi
      frc(3, i) = frc(3, i) + fzi


        vxx = vxx + xa * xij
        vxy = vxy + xa * yij
        vxz = vxz + xa * zij
        vyy = vyy + ya * yij
        vyz = vyz + ya * zij
        vzz = vzz + za * zij

        ! Correct the molecular virial:

        molvir(1, 1) = molvir(1, 1) - xa * rel_crd(1, i)
        molvir(2, 1) = molvir(2, 1) - ya * rel_crd(1, i)
        molvir(3, 1) = molvir(3, 1) - za * rel_crd(1, i)

        molvir(1, 2) = molvir(1, 2) - xa * rel_crd(2, i)
        molvir(2, 2) = molvir(2, 2) - ya * rel_crd(2, i)
        molvir(3, 2) = molvir(3, 2) - za * rel_crd(2, i)

        molvir(1, 3) = molvir(1, 3) - xa * rel_crd(3, i)
        molvir(2, 3) = molvir(2, 3) - ya * rel_crd(3, i)
        molvir(3, 3) = molvir(3, 3) - za * rel_crd(3, i)

    end if

    if (atm_owner_map(j) .eq. mytaskid) then
      frc(1, j) = frc(1, j) + fxj 
      frc(2, j) = frc(2, j) + fyj
      frc(3, j) = frc(3, j) + fzj
    end if

    if (atm_owner_map(k) .eq. mytaskid) then
      frc(1, k) = frc(1, k) + fxk
      frc(2, k) = frc(2, k) + fyk
      frc(3, k) = frc(3, k) + fzk
    end if

    if (atm_owner_map(l) .eq. mytaskid) then
      frc(1, l) = frc(1, l) + fxl
      frc(2, l) = frc(2, l) + fyl
      frc(3, l) = frc(3, l) + fzl

        ! Correct the molecular virial:

        molvir(1, 1) = molvir(1, 1) + xa * rel_crd(1, l)
        molvir(2, 1) = molvir(2, 1) + ya * rel_crd(1, l)
        molvir(3, 1) = molvir(3, 1) + za * rel_crd(1, l)

        molvir(1, 2) = molvir(1, 2) + xa * rel_crd(2, l)
        molvir(2, 2) = molvir(2, 2) + ya * rel_crd(2, l)
        molvir(3, 2) = molvir(3, 2) + za * rel_crd(2, l)

        molvir(1, 3) = molvir(1, 3) + xa * rel_crd(3, l)
        molvir(2, 3) = molvir(2, 3) + ya * rel_crd(3, l)
        molvir(3, 3) = molvir(3, 3) + za * rel_crd(3, l)

    end if
    
  end do

  enbp = enbpl * scnb0
  eelp = eelpl * scee0
  ep   = epl

    ! Save the virials.

    e14_vir(1, 1) = e14_vir(1, 1) + vxx
    e14_vir(1, 2) = e14_vir(1, 2) + vxy
    e14_vir(2, 1) = e14_vir(2, 1) + vxy
    e14_vir(1, 3) = e14_vir(1, 3) + vxz
    e14_vir(3, 1) = e14_vir(3, 1) + vxz
    e14_vir(2, 2) = e14_vir(2, 2) + vyy
    e14_vir(2, 3) = e14_vir(2, 3) + vyz
    e14_vir(3, 2) = e14_vir(3, 2) + vyz
    e14_vir(3, 3) = e14_vir(3, 3) + vzz

  return

end subroutine get_dihed_energy

end module dihedrals_mod
