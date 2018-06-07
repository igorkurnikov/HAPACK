#ifndef DIRFRC_NOVEC

#ifdef BUILD_SHORT_ENE_VEC_2CUT

!*******************************************************************************
!
! Internal Subroutine:  short_ene_vec_2cut
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file, but all
!               versions here should only be invoked if DIRFRC_NOVEC is
!               undefined
!*******************************************************************************

subroutine short_ene_vec_2cut(img_frc, img, ef_tbl, eed_cub, ico, &
                              ipairs_sublst, img_iac, cn1, cn2, x_tran)

#else

!*******************************************************************************
!
! Internal Subroutine:  short_ene_vec
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file, but all
!               versions here should only be invoked if DIRFRC_NOVEC is
!               undefined
!*******************************************************************************

subroutine short_ene_vec(img_frc, img, ef_tbl, eed_cub, ico, ipairs_sublst, &
                         img_iac, cn1, cn2, x_tran)

#endif /* BUILD_SHORT_ENE_VEC_2CUT */

  implicit none

! Formal arguments:

  double precision              :: img_frc(3, *)
  type(img_rec), intent(in)     :: img(*)
  double precision, intent(in)  :: ef_tbl(*)
  double precision, intent(in)  :: eed_cub(*)
  integer, intent(in)           :: ico(*)
  integer, intent(in)           :: ipairs_sublst(*)
  integer, intent(in)           :: img_iac(*)
  double precision, intent(in)  :: cn1(*), cn2(*)
  double precision, intent(in)  :: x_tran(1:3, 0:17)

! Local variables:

  integer, parameter            :: mask27 = Z"07FFFFFF"
  double precision, parameter   :: half = 1.d0/2.d0
  double precision, parameter   :: third = 1.d0/3.d0
  integer, parameter            :: vec_size = 128

  double precision      :: dumx, dumy, dumz

  double precision      :: cgi
  double precision      :: cgi_cgj
  double precision      :: b0, b1
  double precision      :: df
  double precision      :: dfx, dfy, dfz
  double precision      :: f6, r6, f12
  double precision      :: f10, r10

  double precision      :: du, du2, du3 ! 'u' is another name for delr2   ! used only with DIRFRC_EFS
  double precision      :: del_efs                                        ! used only with DIRFRC_EFS
  double precision      :: dens_efs                                       ! used only with DIRFRC_EFS
  double precision      :: lowest_efs_u                                   ! used only with DIRFRC_EFS

  ! Variables used with erfc switch table; name are historical:

  double precision      :: switch
  double precision      :: d_switch_dx
  double precision      :: x, dx, e3dx, e4dx2
  double precision      :: delr, delrinv

  integer               :: iaci
  integer               :: ic
  integer               :: ind
  integer               :: img_j
  integer               :: itran
  integer               :: sublst_head, sublst_max
  integer               :: vec_idx, vec_max
#ifdef SLOW_INDIRECTVEC
  integer               :: enc_img, vec_cnt
#else
  integer               :: nxt_idx, nxt_cnt
#endif

  double precision      :: delr2, delr2inv
#ifdef SLOW_INDIRECTVEC
  double precision      :: delx, dely, delz
#else
  integer               :: nxt(vec_size)
#endif

  ! Arrays that seem to cut cache problems on some machines...

  integer               :: img_j_vec(vec_size)
  double precision      :: del_vec(3, vec_size)
  double precision      :: delr2_vec(vec_size)

  if( dirfrc_efs ) then
    dens_efs = efs_tbl_dens
    del_efs = 1.d0 / dens_efs
    lowest_efs_u = lowest_efs_delr2
  end if

! First loop over the ee evaluation-only pairs:

  dumx = 0.d0
  dumy = 0.d0
  dumz = 0.d0

  cgi = img(img_i)%qterm

  sublst_head = 0
  sublst_max = ee_eval_cnt

  do while (sublst_head .lt. sublst_max)

    vec_max = min(sublst_max - sublst_head, vec_size)

#ifdef SLOW_INDIRECTVEC

    vec_cnt = 0

#ifdef DIRFRC_COMTRANS
    if (common_tran .eq. 1) then

      do vec_idx = 1, vec_max
        img_j = ipairs_sublst(sublst_head + vec_idx)
        delx = img(img_j)%x + x_tran(1, 13)
        dely = img(img_j)%y + x_tran(2, 13)
        delz = img(img_j)%z + x_tran(3, 13)
        delr2 = delx * delx + dely * dely + delz * delz
#ifdef BUILD_SHORT_ENE_VEC_2CUT
        if (delr2 .lt. es_cut2) then
#else
        if (delr2 .lt. max_nb_cut2) then
#endif
          vec_cnt = vec_cnt + 1
          img_j_vec(vec_cnt) = img_j
          del_vec(1, vec_cnt) = delx
          del_vec(2, vec_cnt) = dely
          del_vec(3, vec_cnt) = delz
          delr2_vec(vec_cnt) = delr2
        end if
      end do

    else
#endif /* DIRFRC_COMTRANS */

      do vec_idx = 1, vec_max
        enc_img = ipairs_sublst(sublst_head + vec_idx)
        img_j = iand(enc_img, mask27)
        itran = ishft(enc_img, -27)
        delx = img(img_j)%x + x_tran(1, itran)
        dely = img(img_j)%y + x_tran(2, itran)
        delz = img(img_j)%z + x_tran(3, itran)
        delr2 = delx * delx + dely * dely + delz * delz
#ifdef BUILD_SHORT_ENE_VEC_2CUT
        if (delr2 .lt. es_cut2) then
#else
        if (delr2 .lt. max_nb_cut2) then
#endif
          vec_cnt = vec_cnt + 1
          img_j_vec(vec_cnt) = img_j
          del_vec(1, vec_cnt) = delx
          del_vec(2, vec_cnt) = dely
          del_vec(3, vec_cnt) = delz
          delr2_vec(vec_cnt) = delr2
        end if
      end do

#ifdef DIRFRC_COMTRANS
    end if
#endif /* DIRFRC_COMTRANS */

#else

#ifdef DIRFRC_COMTRANS
    if (common_tran .eq. 1) then

      do vec_idx = 1, vec_max
        img_j = ipairs_sublst(sublst_head + vec_idx)
        img_j_vec(vec_idx) = img_j
        del_vec(1, vec_idx) = img(img_j)%x + x_tran(1, 13)
        del_vec(2, vec_idx) = img(img_j)%y + x_tran(2, 13)
        del_vec(3, vec_idx) = img(img_j)%z + x_tran(3, 13)
        delr2_vec(vec_idx) = del_vec(1, vec_idx) * del_vec(1, vec_idx) + &
                             del_vec(2, vec_idx) * del_vec(2, vec_idx) + &
                             del_vec(3, vec_idx) * del_vec(3, vec_idx)
      end do

    else
#endif /* DIRFRC_COMTRANS */

      do vec_idx = 1, vec_max
        img_j = iand(ipairs_sublst(sublst_head + vec_idx), mask27)
        itran = ishft(ipairs_sublst(sublst_head + vec_idx), -27)
        img_j_vec(vec_idx) = img_j
        del_vec(1, vec_idx) = img(img_j)%x + x_tran(1, itran)
        del_vec(2, vec_idx) = img(img_j)%y + x_tran(2, itran)
        del_vec(3, vec_idx) = img(img_j)%z + x_tran(3, itran)
        delr2_vec(vec_idx) = del_vec(1, vec_idx) * del_vec(1, vec_idx) + &
                             del_vec(2, vec_idx) * del_vec(2, vec_idx) + &
                             del_vec(3, vec_idx) * del_vec(3, vec_idx)
      end do

#ifdef DIRFRC_COMTRANS
    end if
#endif /* DIRFRC_COMTRANS */

    nxt_cnt = 0
    do vec_idx = 1, vec_max
#ifdef BUILD_SHORT_ENE_VEC_2CUT
      if (delr2_vec(vec_idx) .lt. es_cut2) then
#else
      if (delr2_vec(vec_idx) .lt. max_nb_cut2) then
#endif
        nxt_cnt = nxt_cnt + 1
        nxt(nxt_cnt) = vec_idx
      end if 
    end do
#endif /* SLOW_INDIRECTVEC */

#ifdef SLOW_INDIRECTVEC
    do vec_idx = 1, vec_cnt
#else
    do nxt_idx = 1, nxt_cnt

      vec_idx = nxt(nxt_idx)
#endif /* SLOW_INDIRECTVEC */

      delr2 = delr2_vec(vec_idx)

      img_j = img_j_vec(vec_idx)
      cgi_cgj = cgi * img(img_j)%qterm

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

      dfx = del_vec(1, vec_idx) * df
      dfy = del_vec(2, vec_idx) * df
      dfz = del_vec(3, vec_idx) * df

      vxx = vxx - del_vec(1, vec_idx) * dfx
      vxy = vxy - del_vec(1, vec_idx) * dfy
      vxz = vxz - del_vec(1, vec_idx) * dfz
      vyy = vyy - del_vec(2, vec_idx) * dfy
      vyz = vyz - del_vec(2, vec_idx) * dfz
      vzz = vzz - del_vec(3, vec_idx) * dfz

      dumx = dumx + dfx
      dumy = dumy + dfy
      dumz = dumz + dfz

      img_frc(1, img_j) = img_frc(1, img_j) + dfx
      img_frc(2, img_j) = img_frc(2, img_j) + dfy
      img_frc(3, img_j) = img_frc(3, img_j) + dfz

    end do

    sublst_head = sublst_head + vec_max

  end do

  sublst_max = sublst_head + full_eval_cnt

  iaci = ntypes_stk * (img_iac(img_i) - 1)
  
  
  do while (sublst_head .lt. sublst_max)

    vec_max = min(sublst_max - sublst_head, vec_size)

#ifdef SLOW_INDIRECTVEC

    vec_cnt = 0

#ifdef DIRFRC_COMTRANS
    if (common_tran .eq. 1) then

      do vec_idx = 1, vec_max
        img_j = ipairs_sublst(sublst_head + vec_idx)
        delx = img(img_j)%x + x_tran(1, 13)
        dely = img(img_j)%y + x_tran(2, 13)
        delz = img(img_j)%z + x_tran(3, 13)
        delr2 = delx * delx + dely * dely + delz * delz
        if (delr2 .lt. max_nb_cut2) then
          vec_cnt = vec_cnt + 1
          img_j_vec(vec_cnt) = img_j
          del_vec(1, vec_cnt) = delx
          del_vec(2, vec_cnt) = dely
          del_vec(3, vec_cnt) = delz
          delr2_vec(vec_cnt) = delr2
        end if
      end do

    else
#endif /* DIRFRC_COMTRANS */

      do vec_idx = 1, vec_max
        enc_img = ipairs_sublst(sublst_head + vec_idx)
        img_j = iand(enc_img, mask27)
        itran = ishft(enc_img, -27)
        delx = img(img_j)%x + x_tran(1, itran)
        dely = img(img_j)%y + x_tran(2, itran)
        delz = img(img_j)%z + x_tran(3, itran)
        delr2 = delx * delx + dely * dely + delz * delz
        if (delr2 .lt. max_nb_cut2) then
          vec_cnt = vec_cnt + 1
          img_j_vec(vec_cnt) = img_j
          del_vec(1, vec_cnt) = delx
          del_vec(2, vec_cnt) = dely
          del_vec(3, vec_cnt) = delz
          delr2_vec(vec_cnt) = delr2
        end if
      end do

#ifdef DIRFRC_COMTRANS
    end if
#endif /* DIRFRC_COMTRANS */

#else

#ifdef DIRFRC_COMTRANS
    if (common_tran .eq. 1) then

      do vec_idx = 1, vec_max
        img_j = ipairs_sublst(sublst_head + vec_idx)
        img_j_vec(vec_idx) = img_j
        del_vec(1, vec_idx) = img(img_j)%x + x_tran(1, 13)
        del_vec(2, vec_idx) = img(img_j)%y + x_tran(2, 13)
        del_vec(3, vec_idx) = img(img_j)%z + x_tran(3, 13)
        delr2_vec(vec_idx) = del_vec(1, vec_idx) * del_vec(1, vec_idx) + &
                             del_vec(2, vec_idx) * del_vec(2, vec_idx) + &
                             del_vec(3, vec_idx) * del_vec(3, vec_idx)
      end do

    else
#endif /* DIRFRC_COMTRANS */

      do vec_idx = 1, vec_max
        img_j = iand(ipairs_sublst(sublst_head + vec_idx), mask27)
        itran = ishft(ipairs_sublst(sublst_head + vec_idx), -27)
        img_j_vec(vec_idx) = img_j
        del_vec(1, vec_idx) = img(img_j)%x + x_tran(1, itran)
        del_vec(2, vec_idx) = img(img_j)%y + x_tran(2, itran)
        del_vec(3, vec_idx) = img(img_j)%z + x_tran(3, itran)
        delr2_vec(vec_idx) = del_vec(1, vec_idx) * del_vec(1, vec_idx) + &
                             del_vec(2, vec_idx) * del_vec(2, vec_idx) + &
                             del_vec(3, vec_idx) * del_vec(3, vec_idx)
      end do

#ifdef DIRFRC_COMTRANS
    end if
#endif /* DIRFRC_COMTRANS */
    nxt_cnt = 0
    do vec_idx = 1, vec_max
      if (delr2_vec(vec_idx) .lt. max_nb_cut2) then
        nxt_cnt = nxt_cnt + 1
        nxt(nxt_cnt) = vec_idx
      end if 
    end do
#endif /* SLOW_INDIRECTVEC */

#ifdef SLOW_INDIRECTVEC
    do vec_idx = 1, vec_cnt
#else
    do nxt_idx = 1, nxt_cnt

      vec_idx = nxt(nxt_idx)
#endif /* SLOW_INDIRECTVEC */

      delr2 = delr2_vec(vec_idx)

      img_j = img_j_vec(vec_idx)
      ic = ico(iaci + img_iac(img_j))

#ifdef BUILD_SHORT_ENE_VEC_2CUT
      if (delr2 .lt. es_cut2) then
#endif
      cgi_cgj = cgi * img(img_j)%qterm

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

        eedvir_stk = eedvir_stk - b1
        eed_stk = eed_stk + b0
        delr2inv = delrinv * delrinv
        df = b1 * delr2inv

      end if

#ifdef BUILD_SHORT_ENE_VEC_2CUT
      else

        delr2inv = 1.d0 / delr2
        df = 0.d0

      end if
#endif /* BUILD_SHORT_ENE_VEC_2CUT */

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

      dfx = del_vec(1, vec_idx) * df
      dfy = del_vec(2, vec_idx) * df
      dfz = del_vec(3, vec_idx) * df

      vxx = vxx - del_vec(1, vec_idx) * dfx
      vxy = vxy - del_vec(1, vec_idx) * dfy
      vxz = vxz - del_vec(1, vec_idx) * dfz
      vyy = vyy - del_vec(2, vec_idx) * dfy
      vyz = vyz - del_vec(2, vec_idx) * dfz
      vzz = vzz - del_vec(3, vec_idx) * dfz

      dumx = dumx + dfx
      dumy = dumy + dfy
      dumz = dumz + dfz

      img_frc(1, img_j) = img_frc(1, img_j) + dfx
      img_frc(2, img_j) = img_frc(2, img_j) + dfy
      img_frc(3, img_j) = img_frc(3, img_j) + dfz

    end do

    sublst_head = sublst_head + vec_max

  end do

  img_frc(1, img_i) = img_frc(1, img_i) - dumx
  img_frc(2, img_i) = img_frc(2, img_i) - dumy
  img_frc(3, img_i) = img_frc(3, img_i) - dumz

  return

#ifdef BUILD_SHORT_ENE_VEC_2CUT
end subroutine short_ene_vec_2cut
#else
end subroutine short_ene_vec
#endif
#endif /* !DIRFRC_NOVEC */
