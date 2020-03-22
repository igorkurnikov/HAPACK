! include find_img_pairs_nocomtran_2cut

  do i = 0, 17
    i_tranvec(1, i) = tranvec(1, i) - x_i 
    i_tranvec(2, i) = tranvec(2, i) - y_i 
    i_tranvec(3, i) = tranvec(3, i) - z_i 
  end do

  z_loop4: &
  do z_bkts_idx = z_bkts_lo, z_bkts_hi
    z_bkt = z_bkts(z_bkts_idx)
    z_tran = z_trans(z_bkts_idx)
    y_loop4: &
    do y_bkts_idx = y_bkts_lo, y_bkts_hi
      yz_bkt = z_bkt + y_bkts(y_bkts_idx)
      yz_tran = z_tran + y_trans(y_bkts_idx)
      x_loop4: &
      do x_bkts_idx = x_bkts_lo, x_bkts_hi
        cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
        cur_tran = x_trans(x_bkts_idx) + yz_tran
        do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
          atm_j = img_atm_map(img_j)
          if (atm_j .gt. 0) then
            dx = img(img_j)%x + i_tranvec(1, cur_tran)
            dy = img(img_j)%y + i_tranvec(2, cur_tran)
            dz = img(img_j)%z + i_tranvec(3, cur_tran)
            r_sq = dx * dx + dy * dy + dz * dz
            if (r_sq .lt. cutlist_sq) then
              if (excl_img_flags(img_j) .eq. 0) then
                if (ico(iaci + img_iac(img_j)) .eq. 0) then
                  if (r_sq .lt. es_cut_sq) then
                    ee_eval_cnt = ee_eval_cnt + 1
                    img_j_ee_eval(ee_eval_cnt) = ior(img_j,ishft(cur_tran, 27))
                  end if
                else
                  full_eval_cnt = full_eval_cnt + 1
                  img_j_full_eval(full_eval_cnt) = ior(img_j,ishft(cur_tran,27))
                end if
              else
                excl_img_j_cnt = excl_img_j_cnt + 1
                excl_img_pairlst(excl_pairlst_len + excl_img_j_cnt) = img_j
              end if
            end if
          else
            atm_j = - atm_j
            if (img_iac(img_j) .ne. 0) then
              x_j = img(img_j)%x
              y_j = img(img_j)%y
              z_j = img(img_j)%z
            else
              if (is_orthog_stk .ne. 0) then
                x_j = fraction(1, atm_j) * x_box
                y_j = fraction(2, atm_j) * y_box
                z_j = fraction(3, atm_j) * z_box
              else
                f1 = fraction(1, atm_j)
                f2 = fraction(2, atm_j)
                f3 = fraction(3, atm_j)
                ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
                ! we can simplify the expression in this critical inner loop
                x_j = f1 * ucell_stk(1, 1) + f2 * ucell_stk(1, 2) + &
                      f3 * ucell_stk(1, 3)
                y_j = f2 * ucell_stk(2, 2) + f3 * ucell_stk(2, 3)
                z_j = f3 * ucell_stk(3, 3)
              end if
              img(img_j)%x = x_j
              img(img_j)%y = y_j
              img(img_j)%z = z_j
              img(img_j)%qterm = charge(atm_j)
              img_iac(img_j) = iac(atm_j)
            end if
            dx = x_j + i_tranvec(1, cur_tran)
            dy = y_j + i_tranvec(2, cur_tran)
            dz = z_j + i_tranvec(3, cur_tran)
            r_sq = dx * dx + dy * dy + dz * dz
            if (r_sq .lt. cutlist_sq) then
              if (excl_img_flags(img_j) .eq. 0) then
                if (ico(iaci + img_iac(img_j)) .eq. 0) then
                  if (r_sq .lt. es_cut_sq) then
                    ee_eval_cnt = ee_eval_cnt + 1
                    img_j_ee_eval(ee_eval_cnt) = ior(img_j,ishft(cur_tran,27))
                  else
                    cycle
                  end if
                else
                  full_eval_cnt = full_eval_cnt + 1
                  img_j_full_eval(full_eval_cnt) = ior(img_j,ishft(cur_tran,27))
                end if
              else
                excl_img_j_cnt = excl_img_j_cnt + 1
                excl_img_pairlst(excl_pairlst_len + excl_img_j_cnt) = img_j
              end if
              img_atm_map(img_j) = atm_j
            end if
          end if
        end do
        if (cur_bkt .eq. i_bkt) exit z_loop4 ! Done with cit buckets
                                             ! image i claims pairs from.
      end do x_loop4
    end do y_loop4
  end do z_loop4
