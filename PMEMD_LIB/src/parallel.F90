#include "copyright.i"

!*******************************************************************************
!
! Module: parallel_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module parallel_mod

  use parallel_dat_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:  get_send_atm_lst
!
! Description: <TBS>
!
!*******************************************************************************

subroutine get_send_atm_lst(atm_cnt, img_atm_map, atm_img_map, atm_owner_map, &
                            off_tbl, send_atm_lst, send_atm_cnts, &
                            used_img_lo, used_img_hi, used_img_range_wraps)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  integer               :: img_atm_map(atm_cnt)
  integer               :: atm_img_map(atm_cnt)
  integer               :: atm_owner_map(atm_cnt)
  integer               :: off_tbl(0:numtasks)
  integer               :: send_atm_lst(atm_cnt)
  integer               :: send_atm_cnts(0 : numtasks - 1)
  integer               :: used_img_lo
  integer               :: used_img_hi
  logical               :: used_img_range_wraps

! Local variables

  integer               :: atm_id, img_id, node
  integer               :: i

  send_atm_cnts(0: numtasks - 1) = 0

  if (used_img_range_wraps) then
    do img_id = used_img_lo, atm_cnt
      atm_id = img_atm_map(img_id)
      if (atm_id .gt. 0) then
        node = atm_owner_map(atm_id)
        send_atm_cnts(node) = send_atm_cnts(node) + 1
        send_atm_lst(off_tbl(node) + send_atm_cnts(node)) = atm_id
      end if
    end do
    do img_id = 1, used_img_hi
      atm_id = img_atm_map(img_id)
      if (atm_id .gt. 0) then
        node = atm_owner_map(atm_id)
        send_atm_cnts(node) = send_atm_cnts(node) + 1
        send_atm_lst(off_tbl(node) + send_atm_cnts(node)) = atm_id
      end if
    end do
  else
    do img_id = used_img_lo, used_img_hi
      atm_id = img_atm_map(img_id)
      if (atm_id .gt. 0) then
        node = atm_owner_map(atm_id)
        send_atm_cnts(node) = send_atm_cnts(node) + 1
        send_atm_lst(off_tbl(node) + send_atm_cnts(node)) = atm_id
      end if
    end do
  end if

  ! Okay, now add the extra used atoms, which are atoms we don't own but
  ! that we have to keep coordinates updated on because they are used in
  ! bond, angle, and dihedral calculations.  A little gotcha...
  ! There are very few of these (like 5-10 at most typically), but they
  ! really screw things up if you don't update their coordinates!
  ! We only need the coordinates, but because the send atoms mechanism is
  ! shared, forces will also get propagated, so it is important to zero the
  ! forces if they are not in use.

  do i = 1, extra_used_atm_cnt
    atm_id = gbl_extra_used_atms(i)
    if (img_atm_map(atm_img_map(atm_id)) .le. 0) then
      ! It is not already marked as used, so we include it in send atms list.
      node = atm_owner_map(atm_id)
      send_atm_cnts(node) = send_atm_cnts(node) + 1
      send_atm_lst(off_tbl(node) + send_atm_cnts(node)) = atm_id
    end if
  end do

  ! Sum up the atoms for which this process sends data elsewhere.  This is a
  ! measure of locality sent to the master and used in load balancing.  This
  ! will be adjusted to an average between load balancings.  The loadbalancing
  ! code controls when this info is actually used.

  do i = 0, numtasks - 1
    my_send_atm_cnts_total = my_send_atm_cnts_total + gbl_send_atm_cnts(i)
  end do

  my_send_atm_cnts_total = my_send_atm_cnts_total - &
                           gbl_send_atm_cnts(mytaskid)

  my_send_atm_cnts_sums = my_send_atm_cnts_sums + 1

! write(0,*)'DBG: task, gbl_send_atm_cnts()=', mytaskid, gbl_send_atm_cnts(:)

  return

end subroutine get_send_atm_lst

#ifdef SLOW_NONBLOCKING_MPI
! This is an inferior implementation for systems that seem unable to handle
! fully async transposes with good i/o overlap:
!*******************************************************************************
!
! Subroutine:  get_img_frc_distribution
!
! Description: <TBS>
!
!*******************************************************************************

subroutine get_img_frc_distribution(atm_cnt, my_atm_cnt, off_tbl, &
                                    recv_taskmap, send_taskmap, &
                                    send_atm_lst, send_atm_cnts, &
                                    recv_atm_lsts, recv_atm_cnts)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  integer               :: my_atm_cnt
  integer               :: off_tbl(0:numtasks)
  integer               :: recv_taskmap(numtasks - 1)
  integer               :: send_taskmap(numtasks - 1)
  integer               :: send_atm_lst(atm_cnt)
  integer               :: send_atm_cnts(0 : numtasks - 1)
  integer               :: recv_atm_lsts(my_atm_cnt, numtasks - 1)
  integer               :: recv_atm_cnts(0 : numtasks - 1)

! Local variables

  integer               :: recv_cnt
  integer               :: recv_stat(mpi_status_size)
  integer               :: recv_task
  integer               :: send_cnt, send_offset
  integer               :: send_task
  integer               :: send_req
  integer               :: send_stat(mpi_status_size)
  integer               :: taskmap_idx

  ! No need to receive from yourself, so there are numtasks - 1 bufs:

  ! Exchange atom lists...

  do taskmap_idx = 1, numtasks - 1

    send_task = send_taskmap(taskmap_idx)
    send_offset = off_tbl(send_task) + 1
    send_cnt = send_atm_cnts(send_task)

    call mpi_isend(send_atm_lst(send_offset), send_cnt, mpi_integer, &
                   send_task, gifd_tag, lib_mpi_comm, send_req, err_code_mpi)

    recv_task = recv_taskmap(taskmap_idx)

    call mpi_recv(recv_atm_lsts(1, taskmap_idx), my_atm_cnt, mpi_integer, &
                  recv_task, gifd_tag, lib_mpi_comm, recv_stat, err_code_mpi) 

    call mpi_get_count(recv_stat, mpi_integer, recv_cnt, err_code_mpi)

    recv_atm_cnts(recv_task) = recv_cnt

    call mpi_wait(send_req, send_stat, err_code_mpi)

  end do

  return

end subroutine get_img_frc_distribution
#else
! This is an implementation for where the n x n comm of the default
! implementation may be a problem - perhaps very high processor count...
#ifdef ALLTOALL_GETIMGFRCDIST
!*******************************************************************************
!
! Subroutine:  get_img_frc_distribution
!
! Description: <TBS>
!
!*******************************************************************************

subroutine get_img_frc_distribution(atm_cnt, my_atm_cnt, off_tbl, taskmap, &
                                    send_atm_lst, send_atm_cnts, &
                                    recv_atm_lsts, recv_atm_cnts)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  integer               :: my_atm_cnt
  integer               :: off_tbl(0:numtasks)
  integer               :: taskmap(numtasks - 1)
  integer               :: send_atm_lst(atm_cnt)
  integer               :: send_atm_cnts(0 : numtasks - 1)
  integer               :: recv_atm_lsts(my_atm_cnt, numtasks - 1)
  integer               :: recv_atm_cnts(0 : numtasks - 1)

! Local variables

  integer               :: sending_task_cnt, receiving_task_cnt
  integer               :: taskmap_idx
  integer               :: recv_task, send_task
  integer               :: recv_cnt, send_cnt, send_offset
  integer               :: stat(mpi_status_size, 2 * (numtasks - 1))
  integer               :: req(2 * (numtasks - 1))

  sending_task_cnt = 0
  receiving_task_cnt = 0

  ! Send counts of atoms that you need, and receive counts of atoms you
  ! should expect...

  call mpi_alltoall(send_atm_cnts, 1, mpi_integer, &
                    recv_atm_cnts, 1, mpi_integer, &
                    lib_mpi_comm, err_code_mpi)

  ! No need to receive from yourself, so there are numtasks - 1 bufs:

  ! Exchange atom lists...

  ! Post the asynchronous receives first:

  do taskmap_idx = 1, numtasks - 1

    recv_task = taskmap(taskmap_idx)
    recv_cnt = recv_atm_cnts(recv_task)

    if (recv_cnt .gt. 0) then
      sending_task_cnt = sending_task_cnt + 1
      call mpi_irecv(recv_atm_lsts(1, taskmap_idx), my_atm_cnt, mpi_integer, &
                     recv_task, gifd_tag, lib_mpi_comm, &
                     req(sending_task_cnt), err_code_mpi) 
    end if

  end do

  ! Now set up and post the asynchronous sends:

  do taskmap_idx = 1, numtasks - 1

    send_task = taskmap(taskmap_idx)
    send_offset = off_tbl(send_task) + 1
    send_cnt = send_atm_cnts(send_task)

    if (send_cnt .gt. 0) then
      receiving_task_cnt = receiving_task_cnt + 1
      call mpi_isend(send_atm_lst(send_offset), send_cnt, &
                     mpi_integer, send_task, gifd_tag, lib_mpi_comm, &
                     req(sending_task_cnt + receiving_task_cnt), err_code_mpi)
    end if

  end do

  ! Wait for all sends and receives to complete:

  call mpi_waitall(sending_task_cnt + receiving_task_cnt, req, stat, &
                   err_code_mpi)

  return

end subroutine get_img_frc_distribution
#else
! This is the default async implementation...
!*******************************************************************************
!
! Subroutine:  get_img_frc_distribution
!
! Description: <TBS>
!
!*******************************************************************************

subroutine get_img_frc_distribution(atm_cnt, my_atm_cnt, off_tbl, taskmap, &
                                    send_atm_lst, send_atm_cnts, &
                                    recv_atm_lsts, recv_atm_cnts)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  integer               :: my_atm_cnt
  integer               :: off_tbl(0:numtasks)
  integer               :: taskmap(numtasks - 1)
  integer               :: send_atm_lst(atm_cnt)
  integer               :: send_atm_cnts(0 : numtasks - 1)
  integer               :: recv_atm_lsts(my_atm_cnt, numtasks - 1)
  integer               :: recv_atm_cnts(0 : numtasks - 1)

! Local variables

  integer               :: wait_call
  integer               :: taskmap_idx
  integer               :: recv_task, send_task
  integer               :: recv_cnt, send_cnt, send_offset
  integer               :: irecv_stat(mpi_status_size)
  integer               :: isend_stat(mpi_status_size, numtasks - 1)
  integer               :: recv_req(numtasks - 1)
  integer               :: send_req(numtasks - 1)

  ! No need to receive from yourself, so there are numtasks - 1 bufs:

  ! Exchange atom lists...

  ! Post the asynchronous receives first:

  do taskmap_idx = 1, numtasks - 1

    recv_task = taskmap(taskmap_idx)

    call mpi_irecv(recv_atm_lsts(1, taskmap_idx), my_atm_cnt, mpi_integer, &
                   recv_task, gifd_tag, lib_mpi_comm, recv_req(taskmap_idx), &
                   err_code_mpi) 
  end do

  ! Now set up and post the asynchronous sends:

  do taskmap_idx = 1, numtasks - 1

    send_task = taskmap(taskmap_idx)

    send_offset = off_tbl(send_task) + 1
    send_cnt = send_atm_cnts(send_task)

    call mpi_isend(send_atm_lst(send_offset), send_cnt, &
                   mpi_integer, send_task, gifd_tag, lib_mpi_comm, &
                   send_req(taskmap_idx), err_code_mpi)
  end do

  ! Wait on and process the pending receive requests:

  do wait_call = 1, numtasks - 1

    call mpi_waitany(numtasks - 1, recv_req, taskmap_idx, irecv_stat, &
                     err_code_mpi)

    recv_task = taskmap(taskmap_idx)

    call mpi_get_count(irecv_stat, mpi_integer, recv_cnt, err_code_mpi)

    recv_atm_cnts(recv_task) = recv_cnt

  end do

  ! Wait for all sends to complete:

  call mpi_waitall(numtasks - 1, send_req, isend_stat, err_code_mpi)

  return

end subroutine get_img_frc_distribution
#endif /* ALLTOALL_GETIMGFRCDIST */
#endif /* SLOW_NONBLOCKING_MPI */

#ifdef SLOW_NONBLOCKING_MPI
!*******************************************************************************
!
! Subroutine:  distribute_img_frcs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine distribute_img_frcs(atm_cnt, img_frc, frc, atm_img_map, &
                               off_tbl, recv_taskmap, send_taskmap, &
                               send_atm_lst, send_atm_cnts, &
                               recv_atm_lsts, recv_atm_cnts, &
                               send_frc_lst, recv_frc_lst, my_atm_lst, my_atm_cnt)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: img_frc(3, atm_cnt)
  double precision      :: frc(3, atm_cnt)
  integer               :: atm_img_map(atm_cnt)
  integer               :: off_tbl(0:numtasks)
  integer               :: recv_taskmap(numtasks - 1)
  integer               :: send_taskmap(numtasks - 1)
  integer               :: send_atm_lst(atm_cnt)
  integer               :: send_atm_cnts(0 : numtasks - 1)
  integer               :: recv_atm_lsts(my_atm_cnt, numtasks - 1)
  integer               :: recv_atm_cnts(0 : numtasks - 1)
  double precision      :: send_frc_lst(3, atm_cnt)
  double precision      :: recv_frc_lst(3, my_atm_cnt)
  integer               :: my_atm_lst(atm_cnt)
  integer               :: my_atm_cnt

! Local variables

  integer               :: atm_id
  integer               :: atm_lst_idx
  integer               :: i, j
  integer               :: taskmap_idx
  integer               :: recv_task, send_task
  integer               :: recv_cnt, send_cnt
  integer               :: recv_stat(mpi_status_size)
  integer               :: send_req
  integer               :: send_stat(mpi_status_size)

  ! Zero your own nonbonded force entries:

  do atm_lst_idx = 1, my_atm_cnt
    i = my_atm_lst(atm_lst_idx)
    frc(:, i) = 0.d0
  end do

  ! Copy your own img_frc info into the frc buffer:

  do i = off_tbl(mytaskid) + 1, &
         off_tbl(mytaskid) + send_atm_cnts(mytaskid)
    atm_id = send_atm_lst(i)
    frc(:, atm_id) = img_frc(:, atm_img_map(atm_id))
  end do

  do taskmap_idx = 1, numtasks - 1

    send_task = send_taskmap(taskmap_idx)
    send_cnt = send_atm_cnts(send_task) * 3

    if (send_cnt .gt. 0) then
      j = 1
      do i = off_tbl(send_task) + 1, &
             off_tbl(send_task) + send_atm_cnts(send_task)
        atm_id = send_atm_lst(i)
        send_frc_lst(:, j) = img_frc(:, atm_img_map(atm_id))
        j = j + 1
      end do
      call mpi_isend(send_frc_lst(1, 1), send_cnt, mpi_double_precision, &
                     send_task, dif_tag, lib_mpi_comm, send_req, err_code_mpi)
    end if

    recv_task = recv_taskmap(taskmap_idx)
    recv_cnt = recv_atm_cnts(recv_task) * 3

    if (recv_cnt .gt. 0) then
      call mpi_recv(recv_frc_lst(1, 1), recv_cnt, mpi_double_precision, &
                    recv_task, dif_tag, lib_mpi_comm, recv_stat, err_code_mpi)

      do i = 1, recv_atm_cnts(recv_task)
        j = recv_atm_lsts(i, taskmap_idx)
        frc(:, j) = frc(:, j) + recv_frc_lst(:, i)
      end do
    end if

    ! Wait for the current send to complete:

    if (send_cnt .gt. 0) then
      call mpi_wait(send_req, send_stat, err_code_mpi)
    end if

  end do

  return

end subroutine distribute_img_frcs
#else
!*******************************************************************************
!
! Subroutine:  distribute_img_frcs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine distribute_img_frcs(atm_cnt, img_frc, frc, atm_img_map, &
                               off_tbl, recv_taskmap, send_taskmap, &
                               send_atm_lst, send_atm_cnts, &
                               recv_atm_lsts, recv_atm_cnts, &
                               send_frc_lst, recv_frc_lsts, my_atm_lst, my_atm_cnt)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: img_frc(3, atm_cnt)
  double precision      :: frc(3, atm_cnt)
  integer               :: atm_img_map(atm_cnt)
  integer               :: off_tbl(0:numtasks)
  integer               :: recv_taskmap(numtasks - 1)
  integer               :: send_taskmap(numtasks - 1)
  integer               :: send_atm_lst(atm_cnt)
  integer               :: send_atm_cnts(0 : numtasks - 1)
  integer               :: recv_atm_lsts(my_atm_cnt, numtasks - 1)
  integer               :: recv_atm_cnts(0 : numtasks - 1)
  double precision      :: send_frc_lst(3, atm_cnt)
  double precision      :: recv_frc_lsts(3, my_atm_cnt, numtasks - 1)
  integer               :: my_atm_lst(atm_cnt)
  integer               :: my_atm_cnt

  ! No need to receive from yourself, so there are numtasks - 1 bufs:

! Local variables:

  integer               :: atm_id, img_id, node, wait_call
  integer               :: atm_lst_idx
  integer               :: i, j
  integer               :: taskmap_idx
  integer               :: recvs_posted, sends_posted
  integer               :: recv_task, send_task
  integer               :: recv_cnt, send_cnt, send_offset
  integer               :: irecv_stat(mpi_status_size)
  integer               :: isend_stat(mpi_status_size, numtasks - 1)
  integer               :: recv_req(numtasks - 1)
  integer               :: send_req(numtasks - 1)

  recvs_posted = 0

  do taskmap_idx = 1, numtasks - 1

    recv_task = recv_taskmap(taskmap_idx)
    recv_cnt = recv_atm_cnts(recv_task) * 3

    if (recv_cnt .gt. 0) then
      call mpi_irecv(recv_frc_lsts(1, 1, taskmap_idx), recv_cnt, &
                     mpi_double_precision, recv_task, dif_tag, &
                     lib_mpi_comm, recv_req(taskmap_idx), err_code_mpi)
      recvs_posted = recvs_posted + 1
    else
      recv_req(taskmap_idx) = MPI_REQUEST_NULL
    end if

  end do

  sends_posted = 0

  do taskmap_idx = 1, numtasks - 1

    send_task = send_taskmap(taskmap_idx)
    send_cnt = send_atm_cnts(send_task) * 3

    if (send_cnt .gt. 0) then

!BEGIN DBG
!       if (off_tbl(send_task + 1) - off_tbl(send_task) .lt. &
!           send_atm_cnts(send_task)) then
!         write(0,*)'WARNING: send_cnt too big!!!'
!       end if
!END DBG

      send_offset = off_tbl(send_task)
      do i = send_offset + 1, send_offset + send_atm_cnts(send_task)
        atm_id = send_atm_lst(i)
        img_id = atm_img_map(atm_id)
        send_frc_lst(:, i) = img_frc(:, img_id)
      end do

      call mpi_isend(send_frc_lst(1, send_offset + 1), send_cnt, &
                     mpi_double_precision, send_task, dif_tag, lib_mpi_comm, &
                     send_req(taskmap_idx), err_code_mpi)

      sends_posted = sends_posted + 1
    else
      send_req(taskmap_idx) = MPI_REQUEST_NULL
    end if

  end do

  ! Zero your own nonbonded force entries:

  do atm_lst_idx = 1, my_atm_cnt
    i = my_atm_lst(atm_lst_idx)
    frc(:, i) = 0.d0
  end do


  ! Copy your own img_frc info into the frc buffer:

  send_offset = off_tbl(mytaskid)
  do i = 1, send_atm_cnts(mytaskid)
    atm_id = send_atm_lst(send_offset + i)
    img_id = atm_img_map(atm_id)
    frc(:, atm_id) = img_frc(:, img_id)
  end do

  ! Wait on and process any pending receive requests:

  do wait_call = 1, recvs_posted
    call mpi_waitany(numtasks - 1, recv_req, taskmap_idx, irecv_stat, &
                     err_code_mpi)

    recv_task = recv_taskmap(taskmap_idx)

! BEGIN DBG
!   call mpi_get_count(irecv_stat, mpi_double_precision, i, err_code_mpi)
!
!   if (i .ne. 3 * recv_atm_cnts(recv_task)) then
!     write(0,*)'WARNING: recv_atm_cnts value wrong!'
!   end if
!
!   if (i .gt. my_atm_cnt * 3) then
!     write(0,*)'WARNING: buffer overflow on recv!'
!   end if
! END DBG
    
    do i = 1, recv_atm_cnts(recv_task)
      j = recv_atm_lsts(i, taskmap_idx)
      frc(:, j) = frc(:, j) + recv_frc_lsts(:, i, taskmap_idx)
    end do
  end do

  ! Wait for all sends to complete:

  if (sends_posted .gt. 0) then
    call mpi_waitall(numtasks - 1, send_req, isend_stat, err_code_mpi)
  end if

  return

end subroutine distribute_img_frcs
#endif /* SLOW_NONBLOCKING_MPI */

#ifdef SLOW_NONBLOCKING_MPI
! This is an inferior implementation for systems that seem unable to handle
! fully async transposes with good i/o overlap:
!*******************************************************************************
!
! Subroutine:  distribute_crds
!
! Description: <TBS>
!
!*******************************************************************************

subroutine distribute_crds(atm_cnt, my_atm_cnt, crd, off_tbl, recv_taskmap, send_taskmap, &
                           recv_atm_lst, recv_atm_cnts, &
                           send_atm_lsts, send_atm_cnts, &
                           send_crd_lst, recv_crd_lst)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  integer               :: my_atm_cnt
  double precision      :: crd(3, atm_cnt)
  integer               :: off_tbl(0:numtasks)
  integer               :: recv_taskmap(numtasks - 1)
  integer               :: send_taskmap(numtasks - 1)
  integer               :: recv_atm_lst(atm_cnt)
  integer               :: recv_atm_cnts(0 : numtasks - 1)
  integer               :: send_atm_lsts(my_atm_cnt, numtasks - 1)
  integer               :: send_atm_cnts(0 : numtasks - 1)
  double precision      :: send_crd_lst(3, my_atm_cnt)
  double precision      :: recv_crd_lst(3, atm_cnt)

! Local variables

  integer               :: i
  integer               :: recv_cnt
  integer               :: recv_offset
  integer               :: recv_stat(mpi_status_size)
  integer               :: recv_task
  integer               :: send_cnt
  integer               :: send_req
  integer               :: send_stat(mpi_status_size)
  integer               :: send_task
  integer               :: taskmap_idx

  ! No need to send to yourself, so there are numtasks - 1 bufs.

  ! The global names gbl_send_atm_lst, cit_sent_atm_cnts, cit_recv_atm_lst,
  ! gbl_recv_atm_cnts all refer to the process required to distribute image
  ! forces to their owners.  Here we are gathering coordinates to be used with
  ! images from the atom owners, so the process is inverted (ie., we actually
  ! use send_atm_lst and send_atm_cnts to receive coordinates, and vice versa). 
  ! We therefore invert the names in the call list to make the names in this
  ! routine consistent.

  do taskmap_idx = 1, numtasks - 1

    send_task = send_taskmap(taskmap_idx)

    do i = 1, send_atm_cnts(send_task)
      send_crd_lst(:, i) = crd(:, send_atm_lsts(i, taskmap_idx))
    end do
  
    send_cnt = send_atm_cnts(send_task) * 3

    if (send_cnt .gt. 0) then
      call mpi_isend(send_crd_lst(1, 1), send_cnt, mpi_double_precision, &
                     send_task, dc_tag, lib_mpi_comm, send_req, err_code_mpi)
    end if

    recv_task = recv_taskmap(taskmap_idx)
    recv_cnt = recv_atm_cnts(recv_task) * 3

    if (recv_cnt .gt. 0) then
      call mpi_recv(recv_crd_lst(1, off_tbl(recv_task) + 1), recv_cnt, &
                    mpi_double_precision, recv_task, dc_tag, &
                    lib_mpi_comm, recv_stat, err_code_mpi)

      recv_offset = off_tbl(recv_task)
      do i = recv_offset + 1, recv_offset + recv_atm_cnts(recv_task)
        crd(:, recv_atm_lst(i)) = recv_crd_lst(:, i)
      end do
    end if

    if (send_cnt .gt. 0) then
      call mpi_wait(send_req, send_stat, err_code_mpi)
    end if

  end do

  return

end subroutine distribute_crds
#else
!*******************************************************************************
!
! Subroutine:  distribute_crds
!
! Description: <TBS>
!
!*******************************************************************************

subroutine distribute_crds(atm_cnt, my_atm_cnt, crd, off_tbl, taskmap, &
                               recv_atm_lst, recv_atm_cnts, &
                               send_atm_lsts, send_atm_cnts, &
                               send_crd_lsts, recv_crd_lst)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  integer               :: my_atm_cnt
  double precision      :: crd(3, atm_cnt)
  integer               :: off_tbl(0:numtasks)
  integer               :: taskmap(numtasks - 1)
  integer               :: recv_atm_lst(atm_cnt)
  integer               :: recv_atm_cnts(0 : numtasks - 1)
  integer               :: send_atm_lsts(my_atm_cnt, numtasks - 1)
  integer               :: send_atm_cnts(0 : numtasks - 1)
  double precision      :: send_crd_lsts(3, my_atm_cnt, numtasks - 1)
  double precision      :: recv_crd_lst(3, atm_cnt)


! Local variables

  integer               :: wait_call
  integer               :: i
  integer               :: taskmap_idx
  integer               :: recvs_posted, sends_posted
  integer               :: recv_task, send_task
  integer               :: recv_cnt, send_cnt, recv_offset
  integer               :: irecv_stat(mpi_status_size)
  integer               :: isend_stat(mpi_status_size, numtasks - 1)
  integer               :: recv_req(numtasks - 1)
  integer               :: send_req(numtasks - 1)

  ! No need to send to yourself, so there are numtasks - 1 bufs.

  ! The global names gbl_send_atm_lst, cit_sent_atm_cnts, cit_recv_atm_lst,
  ! gbl_recv_atm_cnts all refer to the process required to distribute image
  ! forces to their owners.  Here we are gathering coordinates to be used with
  ! images from the atom owners, so the process is inverted (ie., we actually
  ! use send_atm_lst and send_atm_cnts to receive coordinates, and vice versa). 
  ! We therefore invert the names in the call list to make the names in this
  ! routine consistent.

  ! Post the asynchronous receives first:

  recvs_posted = 0

  do taskmap_idx = 1, numtasks - 1

    recv_task = taskmap(taskmap_idx)
    recv_cnt = recv_atm_cnts(recv_task) * 3

    if (recv_cnt .gt. 0) then

      call mpi_irecv(recv_crd_lst(1, off_tbl(recv_task) + 1), recv_cnt, &
                     mpi_double_precision, recv_task, dc_tag, &
                     lib_mpi_comm, recv_req(taskmap_idx), err_code_mpi)
      recvs_posted = recvs_posted + 1
    else
      recv_req(taskmap_idx) = MPI_REQUEST_NULL
    end if
  end do

  ! Now set up and post the asynchronous sends:

  sends_posted = 0

  do taskmap_idx = 1, numtasks - 1

    send_task = taskmap(taskmap_idx)

    do i = 1, send_atm_cnts(send_task)
      send_crd_lsts(:, i, taskmap_idx) = crd(:, send_atm_lsts(i, taskmap_idx))
    end do
  
    send_cnt = send_atm_cnts(send_task) * 3

! BEGIN DBG
!       if (my_atm_cnt .lt. send_atm_cnts(send_task)) then
!         write(0,*)'WARNING: send_cnt too big, distribute_crds!!!'
!       end if
! END DBG

    if (send_cnt .gt. 0) then
      call mpi_isend(send_crd_lsts(1, 1, taskmap_idx), send_cnt, &
                     mpi_double_precision, send_task, dc_tag, lib_mpi_comm, &
                     send_req(taskmap_idx), err_code_mpi)
      sends_posted = sends_posted + 1
    else
      send_req(taskmap_idx) = MPI_REQUEST_NULL
    end if
  end do

  ! Wait on and process the pending receive requests:

  do wait_call = 1, recvs_posted

    call mpi_waitany(numtasks - 1, recv_req, taskmap_idx, irecv_stat, &
                     err_code_mpi)

    recv_task = taskmap(taskmap_idx)
    recv_offset = off_tbl(recv_task)

! BEGIN DBG
!   call mpi_get_count(irecv_stat, mpi_double_precision, i, err_code_mpi)
!
!   if (i .ne. 3 * recv_atm_cnts(recv_task)) then
!     write(0,*)'WARNING: recv_atm_cnts value wrong, distribute_crds!'
!   end if
!
!   if (recv_atm_cnts(recv_task) .gt. &
!       off_tbl(recv_task + 1) - off_tbl(recv_task)) then
!     write(0,*)'WARNING: buffer overflow on recv, distribute_crds!'
!   end if
! END DBG
    
    do i = recv_offset + 1, recv_offset + recv_atm_cnts(recv_task)
      crd(:, recv_atm_lst(i)) = recv_crd_lst(:, i)
    end do

  end do

  ! Wait for all sends to complete:

  if (sends_posted .gt. 0) then
    call mpi_waitall(numtasks - 1, send_req, isend_stat, err_code_mpi)
  end if

  return

end subroutine distribute_crds
#endif /* SLOW_NONBLOCKING_MPI */

!*******************************************************************************
!
! Subroutine:  mpi_allgathervec
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine mpi_allgathervec(atm_cnt, vec, atm_owner_map, my_atm_lst, my_atm_cnt)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: vec(3, atm_cnt)
  integer               :: atm_owner_map(atm_cnt)
  integer               :: my_atm_lst(atm_cnt)
  integer               :: my_atm_cnt
  
  call wrapped_mpi_allgathervec(atm_cnt, vec, dbl_mpi_send_buf, dbl_mpi_recv_buf, &
                                atm_owner_map, my_atm_lst, my_atm_cnt)

  return

end subroutine mpi_allgathervec

!*******************************************************************************
!
! Subroutine:  wrapped_mpi_allgathervec
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine wrapped_mpi_allgathervec(atm_cnt, vec, send_buf, recv_buf, atm_owner_map, my_atm_lst, my_atm_cnt)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: vec(3, atm_cnt)
  double precision      :: send_buf(3, *)
  double precision      :: recv_buf(3, *)
  integer               :: atm_owner_map(atm_cnt)
  integer               :: my_atm_lst(atm_cnt)
  integer               :: my_atm_cnt

! Local variables:

  integer               :: atm_lst_idx, atm_idx
  integer               :: my_sendoff
  integer               :: node
  integer               :: cur_node_off(0: numtasks - 1)

  ! We marshall, send/recv, and unmarshall from all processes including the
  ! current process in order to keep the current process from having to do a
  ! conditional across all the data in the unmarshalling.

  my_sendoff = gbl_atm_offsets(mytaskid)
  
  do atm_lst_idx = 1, my_atm_cnt
    atm_idx = my_atm_lst(atm_lst_idx)
    send_buf(:, my_sendoff + atm_lst_idx) = vec(:, atm_idx)
  end do

  call mpi_allgatherv(send_buf(1, my_sendoff + 1), &
                      gbl_vec_rcvcnts(mytaskid), &
                      mpi_double_precision, &
                      recv_buf, gbl_vec_rcvcnts, gbl_vec_offsets, &
                      mpi_double_precision, lib_mpi_comm, err_code_mpi)

  cur_node_off(0:numtasks - 1) = gbl_atm_offsets(0:numtasks - 1) + 1
  
  do atm_idx = 1, atm_cnt
    node = atm_owner_map(atm_idx)
    vec(:, atm_idx) = recv_buf(:, cur_node_off(node))
    cur_node_off(node) = cur_node_off(node) + 1
  end do

  return

end subroutine wrapped_mpi_allgathervec

!*******************************************************************************
!
! Subroutine:  mpi_gathervec
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine mpi_gathervec(atm_cnt, vec, atm_owner_map, my_atm_lst, my_atm_cnt)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: vec(3, atm_cnt)
  integer               :: atm_owner_map(atm_cnt)
  integer               :: my_atm_lst(atm_cnt)
  integer               :: my_atm_cnt
  
  call wrapped_mpi_gathervec(atm_cnt, vec, dbl_mpi_send_buf, dbl_mpi_recv_buf, & 
                             atm_owner_map, my_atm_lst, my_atm_cnt)

  return

end subroutine mpi_gathervec

!*******************************************************************************
!
! Subroutine:  wrapped_mpi_gathervec
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine wrapped_mpi_gathervec(atm_cnt, vec, send_buf, recv_buf, atm_owner_map, my_atm_lst, my_atm_cnt)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: vec(3, atm_cnt)
  double precision      :: send_buf(3, *)
  double precision      :: recv_buf(3, *)
  integer               :: atm_owner_map(atm_cnt)
  integer               :: my_atm_lst(atm_cnt)
  integer               :: my_atm_cnt

! Local variables:

  integer               :: i
  integer               :: atm_lst_idx, atm_idx
  integer               :: my_sendoff
  integer               :: node
  integer               :: cur_node_off(0: numtasks - 1)
  
  ! We marshall, send/recv, and unmarshall from all processes including root, in
  ! order to keep the root from having to do a conditional across all the
  ! data in the unmarshalling.

  my_sendoff = gbl_atm_offsets(mytaskid)
  
  do atm_lst_idx = 1, my_atm_cnt
    atm_idx = my_atm_lst(atm_lst_idx)
    send_buf(:, my_sendoff + atm_lst_idx) = vec(:, atm_idx)
  end do
  
  call mpi_gatherv(send_buf(1, my_sendoff + 1), &
                   gbl_vec_rcvcnts(mytaskid), &
                   mpi_double_precision, &
                   recv_buf, gbl_vec_rcvcnts, gbl_vec_offsets, &
                   mpi_double_precision, &
                   0, lib_mpi_comm, err_code_mpi)

  if (mytaskid .eq. 0) then
    cur_node_off(0:numtasks - 1) = gbl_atm_offsets(0:numtasks - 1) + 1
    do atm_idx = 1, atm_cnt
      node = atm_owner_map(atm_idx)
      vec(:, atm_idx) = recv_buf(:, cur_node_off(node))
      cur_node_off(node) = cur_node_off(node) + 1
    end do
  end if

  return

end subroutine wrapped_mpi_gathervec

!*******************************************************************************
!
! Subroutine:  set_minimum_mpi_bufs_size
!
! Description:  This routine grows dbl_mpi_send_buf and dbl_mpi_recv_buf as
!               required.  These buffers will always be the maximum value of
!               min_buf_size seen so far.  There is no shrinking, and the 
!               buffers are always the same size.
!
!*******************************************************************************

subroutine set_minimum_mpi_bufs_size(min_buf_size)

  implicit none

! Formal arguments:

  integer       :: min_buf_size

! Local variables:

  integer       :: alloc_failed

  if (min_buf_size .le. siz_dbl_mpi_bufs) return

  if (allocated(dbl_mpi_send_buf)) deallocate(dbl_mpi_send_buf)
 
  allocate(dbl_mpi_send_buf(min_buf_size), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error
  
  if (allocated(dbl_mpi_recv_buf)) deallocate(dbl_mpi_recv_buf)
  
  allocate(dbl_mpi_recv_buf(min_buf_size), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  siz_dbl_mpi_bufs = min_buf_size

  return

end subroutine set_minimum_mpi_bufs_size

!*******************************************************************************
!
! Subroutine:  gb_frcs_distrib
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine gb_frcs_distrib(atm_cnt, frc, recv_buf)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: frc(3 * atm_cnt)
  double precision      :: recv_buf(*)

! Local variables:

  integer               :: my_recvoff
  integer               :: my_recvcnt

  my_recvoff = gbl_vec_offsets(mytaskid)
  my_recvcnt = gbl_vec_rcvcnts(mytaskid)

  call mpi_reduce_scatter(frc, recv_buf, gbl_vec_rcvcnts, &
                          mpi_double_precision, mpi_sum, lib_mpi_comm, &
                          err_code_mpi)

  frc(my_recvoff + 1 : my_recvoff + my_recvcnt) = &
    recv_buf(1 : my_recvcnt)

  return

end subroutine gb_frcs_distrib

!*******************************************************************************
!
! Subroutine:  gb_mpi_allgathervec
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine gb_mpi_allgathervec(atm_cnt, vec)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: vec(3, atm_cnt)

! Local variables:

  call mpi_allgatherv(vec(1, gbl_atm_offsets(mytaskid) + 1), &
                      gbl_vec_rcvcnts(mytaskid), mpi_double_precision, &
                      vec, gbl_vec_rcvcnts, gbl_vec_offsets, &
                      mpi_double_precision, lib_mpi_comm, err_code_mpi)

  return

end subroutine gb_mpi_allgathervec

!*******************************************************************************
!
! Subroutine:  gb_mpi_gathervec
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine gb_mpi_gathervec(atm_cnt, vec)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: vec(3, atm_cnt)

! Local variables:

  call mpi_gatherv(vec(1, gbl_atm_offsets(mytaskid) + 1), &
                   gbl_vec_rcvcnts(mytaskid), mpi_double_precision, &
                   vec, gbl_vec_rcvcnts, gbl_vec_offsets, &
                   mpi_double_precision, 0, lib_mpi_comm, err_code_mpi)

  return

end subroutine gb_mpi_gathervec

subroutine distribute_crds_proxy(atm_cnt, my_atm_cnt, crd)

    implicit none
    
    ! Formal Arguments:
    
    integer            :: atm_cnt
    integer            :: my_atm_cnt
    double precision   :: crd(3,*)

#ifdef SLOW_NONBLOCKING_MPI
    call distribute_crds(atm_cnt, my_atm_cnt, crd, gbl_atm_offsets, gbl_inv_taskmap, &
                                       gbl_taskmap, &
                                       gbl_send_atm_lst, gbl_send_atm_cnts, &
                                       gbl_recv_atm_lsts, gbl_recv_atm_cnts, &
                                       dbl_mpi_send_buf, dbl_mpi_recv_buf)
#else
    call distribute_crds(atm_cnt, my_atm_cnt, crd, gbl_atm_offsets, gbl_taskmap, &
                                       gbl_send_atm_lst, gbl_send_atm_cnts, &
                                       gbl_recv_atm_lsts, gbl_recv_atm_cnts, &
                                       dbl_mpi_send_buf, dbl_mpi_recv_buf)
#endif
                                 
end subroutine distribute_crds_proxy
                                 
end module parallel_mod