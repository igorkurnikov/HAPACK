#ifdef DIRFRC_NOVEC

#ifdef BUILD_SHORT_ENE_NOVEC_2CUT

!*******************************************************************************
!
! Internal Subroutine:  short_ene_novec_2cut
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file, but all
!               versions here should only be invoked if DIRFRC_NOVEC is
!               defined.
!*******************************************************************************

subroutine short_ene_novec_2cut(img_frc, img, ef_tbl, eed_cub, ico, &
                                ipairs_sublst, img_iac, cn1, cn2, x_tran)

#else

!*******************************************************************************
!
! Internal Subroutine:  short_ene_novec
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file, but all
!               versions here should only be invoked if DIRFRC_NOVEC is
!               defined.
!*******************************************************************************

subroutine short_ene_novec(img_frc, img, ef_tbl, eed_cub, ico, ipairs_sublst, &
                           img_iac, cn1, cn2, x_tran)

#endif /* BUILD_SHORT_ENE_NOVEC_2CUT */

  implicit none

! Formal arguments:

  double precision              :: img_frc(3, *)
  type(img_rec), intent(in)     :: img(*)
  double precision, intent(in)  :: ef_tbl(*)
  double precision, intent(in)  :: eed_cub(*)
  integer, intent(in)           :: ico(*)
  integer                       :: ipairs_sublst(*)
  integer, intent(in)           :: img_iac(*)
  double precision, intent(in)  :: cn1(*), cn2(*)
  double precision, intent(in)  :: x_tran(1:3, 0:17)

! Local variables:

  integer, parameter            :: mask27 = Z"07FFFFFF"
  double precision, parameter   :: half = 1.d0/2.d0
  double precision, parameter   :: third = 1.d0/3.d0

  double precision      :: dumx, dumy, dumz

  double precision      :: cgi
  double precision      :: cgi_cgj
  double precision      :: b0, b1
  double precision      :: df
  double precision      :: dfx, dfy, dfz
  double precision      :: f6, r6, f12
  double precision      :: f10, r10

  double precision      :: du, du2, du3 ! 'u' is another name for delr2  ! Used only for DIRFRC_EFS
  double precision      :: del_efs                                       ! Used only for DIRFRC_EFS
  double precision      :: dens_efs                                      ! Used only for DIRFRC_EFS
  double precision      :: lowest_efs_u                                  ! Used only for DIRFRC_EFS

  ! Variables used with erfc switch table; name are historical:

  double precision      :: switch
  double precision      :: d_switch_dx
  double precision      :: x, dx, e3dx, e4dx2
  double precision      :: delr, delrinv

  integer               :: iaci
  integer               :: ic
  integer               :: ind
  integer               :: nxt_img_j, img_j
  integer               :: itran
  integer               :: sublst_idx
  integer               :: saved_pairlist_val

  double precision      :: nxt_delx, nxt_dely, nxt_delz
  double precision      :: delx, dely, delz, delr2, delr2inv

  if( dirfrc_efs ) then
    dens_efs = efs_tbl_dens
    del_efs = 1.d0 / dens_efs
    lowest_efs_u = lowest_efs_delr2
  end if

! First loop over the ee evaluation-only pairs:

  dumx = 0.d0
  dumy = 0.d0
  dumz = 0.d0

  cgi = img(img_i)%charge

  ! The pairlist must have one dummy end entry to cover reading past the
  ! end of the list...

  saved_pairlist_val = ipairs_sublst(ee_eval_cnt + full_eval_cnt + 1)

  ipairs_sublst(ee_eval_cnt + full_eval_cnt + 1) = &
    ipairs_sublst(ee_eval_cnt + full_eval_cnt)

#ifdef DIRFRC_COMTRANS
  if (common_tran .eq. 1) then
    nxt_img_j = ipairs_sublst(1)
    itran = 13
  else
#endif /* DIRFRC_COMTRANS */
    nxt_img_j = iand(ipairs_sublst(1), mask27)
    itran = ishft(ipairs_sublst(1), -27)
#ifdef DIRFRC_COMTRANS
  end if
#endif /* DIRFRC_COMTRANS */

  nxt_delx = img(nxt_img_j)%x + x_tran(1, itran)
  nxt_dely = img(nxt_img_j)%y + x_tran(2, itran)
  nxt_delz = img(nxt_img_j)%z + x_tran(3, itran)

  do sublst_idx = 2, ee_eval_cnt + 1

    img_j = nxt_img_j
    delx = nxt_delx
    dely = nxt_dely
    delz = nxt_delz

#ifdef DIRFRC_COMTRANS
    if (common_tran .eq. 1) then
      nxt_img_j = ipairs_sublst(sublst_idx)
      itran = 13
    else
#endif /* DIRFRC_COMTRANS */
      nxt_img_j = iand(ipairs_sublst(sublst_idx), mask27)
      itran = ishft(ipairs_sublst(sublst_idx), -27)
#ifdef DIRFRC_COMTRANS
    end if
#endif /* DIRFRC_COMTRANS */


    nxt_delx = img(nxt_img_j)%x + x_tran(1, itran)
    nxt_dely = img(nxt_img_j)%y + x_tran(2, itran)
    nxt_delz = img(nxt_img_j)%z + x_tran(3, itran)

    delr2 = delx * delx + dely * dely + delz * delz


#ifdef BUILD_SHORT_ENE_NOVEC_2CUT
    if (delr2 .lt. es_cut2) then
#else
    if (delr2 .lt. max_nb_cut2) then
#endif

      cgi_cgj = cgi * img(img_j)%charge

      if ( dirfrc_efs .and. (delr2 .ge. lowest_efs_u) ) then

        ! Do the Coulomb part of the direct sum using efs: 

        ind = int(dens_efs * delr2)
        du = delr2 - dble(ind) * del_efs
        du2 = du * du
        du3 = du * du2
        ind = ishft(ind, 3)             ! 8 * ind

        b0 = cgi_cgj * (ef_tbl(1 + ind) + du * ef_tbl(2 + ind) + &
             du2 * ef_tbl(3 + ind) + du3 * ef_tbl(4 + ind))

        df = cgi_cgj * (ef_tbl(5 + ind) + du * ef_tbl(6 + ind) + &
             du2 * ef_tbl(7 + ind) + du3 * ef_tbl(8 + ind))

        eed_stk = eed_stk + b0

        eedvir_stk = eedvir_stk - df * delr2

      else
      
        ! Do the Coulomb part of the direct sum using erfc spline table: 

        delr = sqrt(delr2)
        delrinv = 1.d0 / delr

        x = dxdr * delr
        ind = int(eedtbdns_stk * x)
        dx = x - dble(ind) * del
        ind = ishft(ind, 2)             ! 4 * ind

        e3dx  = dx * eed_cub(3 + ind)
        e4dx2 = dx * dx *  eed_cub(4 + ind)

        switch = eed_cub(1 + ind) + &
                 dx * (eed_cub(2 + ind) + half * (e3dx + third * e4dx2))

        d_switch_dx = eed_cub(2 + ind) + e3dx + half * e4dx2

        b0 = cgi_cgj * delrinv * switch
        b1 = (b0 - cgi_cgj * d_switch_dx * dxdr)

        eedvir_stk = eedvir_stk - b1
        eed_stk = eed_stk + b0
        df = b1 * delrinv * delrinv

      end if

      dfx = delx * df
      dfy = dely * df
      dfz = delz * df

      vxx = vxx - delx * dfx
      vxy = vxy - delx * dfy
      vxz = vxz - delx * dfz
      vyy = vyy - dely * dfy
      vyz = vyz - dely * dfz
      vzz = vzz - delz * dfz

      dumx = dumx + dfx
      dumy = dumy + dfy
      dumz = dumz + dfz

      img_frc(1, img_j) = img_frc(1, img_j) + dfx
      img_frc(2, img_j) = img_frc(2, img_j) + dfy
      img_frc(3, img_j) = img_frc(3, img_j) + dfz

    end if

  end do

  iaci = ntypes_stk * (img_iac(img_i) - 1)

  do sublst_idx = ee_eval_cnt + 2, ee_eval_cnt + full_eval_cnt + 1

    img_j = nxt_img_j
    delx = nxt_delx
    dely = nxt_dely
    delz = nxt_delz

#ifdef DIRFRC_COMTRANS
    if (common_tran .eq. 1) then
      nxt_img_j = ipairs_sublst(sublst_idx)
      itran = 13
    else
#endif /* DIRFRC_COMTRANS */
      nxt_img_j = iand(ipairs_sublst(sublst_idx), mask27)
      itran = ishft(ipairs_sublst(sublst_idx), -27)
#ifdef DIRFRC_COMTRANS
    end if
#endif /* DIRFRC_COMTRANS */


    nxt_delx = img(nxt_img_j)%x + x_tran(1, itran)
    nxt_dely = img(nxt_img_j)%y + x_tran(2, itran)
    nxt_delz = img(nxt_img_j)%z + x_tran(3, itran)

    delr2 = delx * delx + dely * dely + delz * delz

    if (delr2 .lt. max_nb_cut2) then

      ic = ico(iaci + img_iac(img_j))


#ifdef BUILD_SHORT_ENE_NOVEC_2CUT
      if (delr2 .lt. es_cut2) then
#endif

      cgi_cgj = cgi * img(img_j)%charge

      if ( dirfrc_efs .and. (delr2 .ge. lowest_efs_u) ) then

        ! Do the Coulomb part of the direct sum using efs: 

        ind = int(dens_efs * delr2)
        du = delr2 - dble(ind) * del_efs
        du2 = du * du
        du3 = du * du2
        ind = ishft(ind, 3)             ! 8 * ind

        b0 = cgi_cgj * (ef_tbl(1 + ind) + du * ef_tbl(2 + ind) + &
             du2 * ef_tbl(3 + ind) + du3 * ef_tbl(4 + ind))

        df = cgi_cgj * (ef_tbl(5 + ind) + du * ef_tbl(6 + ind) + &
             du2 * ef_tbl(7 + ind) + du3 * ef_tbl(8 + ind))

        eed_stk = eed_stk + b0

        eedvir_stk = eedvir_stk - df * delr2

        delr2inv = 1.d0 / delr2

      else
      
        ! Do the Coulomb part of the direct sum using erfc spline table: 

        delr = sqrt(delr2)
        delrinv = 1.d0 / delr

        x = dxdr * delr
        ind = int(eedtbdns_stk * x)
        dx = x - dble(ind) * del
        ind = ishft(ind, 2)             ! 4 * ind

        e3dx  = dx * eed_cub(3 + ind)
        e4dx2 = dx * dx *  eed_cub(4 + ind)

        switch = eed_cub(1 + ind) + &
                 dx * (eed_cub(2 + ind) + half * (e3dx + third * e4dx2))

        d_switch_dx = eed_cub(2 + ind) + e3dx + half * e4dx2

        b0 = cgi_cgj * delrinv * switch
        b1 = (b0 - cgi_cgj * d_switch_dx * dxdr)

        eed_stk = eed_stk + b0
        eedvir_stk = eedvir_stk - b1
        delr2inv = delrinv * delrinv
        df = b1 * delrinv * delrinv
        
      end if

#ifdef BUILD_SHORT_ENE_NOVEC_2CUT
      else
        delr2inv = 1.d0 / delr2
        df = 0.d0
      end if
#endif

      if (ic .gt. 0) then
        r6 = delr2inv * delr2inv * delr2inv
        f6 = cn2(ic) * r6
        f12 = cn1(ic) * (r6 * r6)
        evdw_stk = evdw_stk + f12 - f6
        df = df + (12.d0 * f12 - 6.d0 * f6) * delr2inv
      else
        ! This code allows 10-12 terms; in many (most?) (all?) cases, the
        ! only "nominal" 10-12 terms are on waters, where the asol and bsol
        ! parameters are always zero; hence we can skip this part.
        ic = - ic
        delr2inv = delrinv * delrinv
        r10 = delr2inv * delr2inv * delr2inv * delr2inv * delr2inv
        f10 = gbl_bsol(ic) * r10
        f12 = gbl_asol(ic) * (r10 * delr2inv)
        ehb_stk = ehb_stk + f12 - f10
        df = df + (12.d0 * f12 - 10.d0 * f10) * delr2inv
      end if

      dfx = delx * df
      dfy = dely * df
      dfz = delz * df

      vxx = vxx - delx * dfx
      vxy = vxy - delx * dfy
      vxz = vxz - delx * dfz
      vyy = vyy - dely * dfy
      vyz = vyz - dely * dfz
      vzz = vzz - delz * dfz

      dumx = dumx + dfx
      dumy = dumy + dfy
      dumz = dumz + dfz

      img_frc(1, img_j) = img_frc(1, img_j) + dfx
      img_frc(2, img_j) = img_frc(2, img_j) + dfy
      img_frc(3, img_j) = img_frc(3, img_j) + dfz

    end if

  end do

  img_frc(1, img_i) = img_frc(1, img_i) - dumx
  img_frc(2, img_i) = img_frc(2, img_i) - dumy
  img_frc(3, img_i) = img_frc(3, img_i) - dumz

  ipairs_sublst(ee_eval_cnt + full_eval_cnt + 1) = saved_pairlist_val

  return

#ifdef BUILD_SHORT_ENE_NOVEC_2CUT
end subroutine short_ene_novec_2cut
#else
end subroutine short_ene_novec
#endif

#endif /* DIRFRC_NOVEC */
