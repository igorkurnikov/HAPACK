#include "copyright.i"

!*******************************************************************************
!
! Module: pme_fft_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module pme_fft_mod

  use fft1d_mod

  implicit none

! None of this stuff is broadcast; fft setup occurs in all tasks.

  ! alpha/beta coefficients
  double precision, allocatable, save, private  :: rc_coefs(:,:)

  ! Handles to the 1d fft's data:

  type(fft1d_rec), pointer, save        :: fft_x_hdl => null()
  type(fft1d_rec), pointer, save        :: fft_y_hdl => null()
  type(fft1d_rec), pointer, save        :: fft_z_hdl => null()

  ! Physical dimensions of fft data.  The data is considered complex, but
  ! stored as double precision in the form:
  ! double precision data(2, fft_x_dim, fft_y_dim, fft_z_dim)
  ! The above index ordering is the default input/output form.
  ! Intermediate data  may be reordered due to transpositions.

  integer, save         :: fft_x_dim
  integer, save         :: fft_y_dim
  integer, save         :: fft_z_dim

  integer, save         :: max_xy_slab_cnt
  integer, save         :: my_xy_slab_cnt
  integer, save         :: my_xy_slab_start
  integer, save         :: max_zx_slab_cnt
  integer, save         :: my_zx_slab_cnt
  integer, save         :: my_zx_slab_start
  integer, save         :: xy_slab_dbl_cnt
  integer, save         :: zx_slab_dbl_cnt

  integer, save                 :: siz_fft_mpi_buf      ! used only in MPI
  integer, allocatable, save    :: xy_slab_cnt(:)       ! used only in MPI
  integer, allocatable, save    :: xy_slab_start(:)     ! used only in MPI
  integer, allocatable, save    :: zx_slab_cnt(:)       ! used only in MPI
  integer, allocatable, save    :: zx_slab_start(:)     ! used only in MPI
  integer, save                 :: fft_taskmap_entries  ! used only in MPI
  integer, allocatable, save    :: fft_taskmap(:)       ! used only in MPI

  ! The FFT workload estimate is an INITIAL value used in loadbalancing.  It is
  ! chosen to typically underestimate fft workload in non-respa simulations;
  ! the number of recip force tasks will then be increased once actual
  ! workloads are known (erring on the underestimate side prevents a high
  ! fft transpose workload from screwing things up).
  
  double precision, parameter   :: fft_workload_estimate = 0.1d0  ! used only in MPI

contains

!*******************************************************************************
!
! Subroutine:  fft_setup
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine fft_setup(recip_numtasks, fft_slab_redist_enabled)

  use gbl_constants_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  ! recip_numtasks returns the initial number of tasks doing recip force calcs.

  integer, intent(out)          :: recip_numtasks
  logical, intent(in)           :: fft_slab_redist_enabled

! Local variables:

  integer               :: alloc_failed
  integer               :: i
  double precision      :: pi2n
  double precision      :: theta

  ! Set fft dimensions - these are physical array dimensions:

  if (mod(nfft1, 2) .ne. 0) then
    write(error_msg, *)'RCFFT requires that nfft1 is even!'
    call mol_mech_error
  end if

  fft_x_dim = (nfft1 + 2) / 2
  fft_y_dim = nfft2
  fft_z_dim = nfft3

  ! Create 1 dimensional fft handles, used to generically access all
  ! 1d fft implementations.

  call fft1d_create(fft_x_hdl, nfft1 / 2)
  call fft1d_create(fft_y_hdl, nfft2)
  call fft1d_create(fft_z_hdl, nfft3)

  ! Allocate and initialize arrays used in rc-cr fft conversions:
 
  if(allocated(rc_coefs)) deallocate(rc_coefs)

  allocate(rc_coefs(2, nfft1 / 2), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  pi2n = 2 * PI / nfft1

  theta = 0.d0
  do i = 1, nfft1 / 2
    rc_coefs(1, i) = cos(theta)
    rc_coefs(2, i) = sin(theta)
    theta = pi2n * i
  end do

  ! Set up multidimensional fft data structures:

  if( numtasks .gt. 1) then
      call  par_fft_setup(recip_numtasks, fft_slab_redist_enabled )
  else
      max_xy_slab_cnt = fft_z_dim
      max_zx_slab_cnt = fft_y_dim
      my_xy_slab_cnt = fft_z_dim
      my_xy_slab_start = 0
      my_zx_slab_cnt = fft_y_dim
      my_zx_slab_start = 0
      xy_slab_dbl_cnt = 2 * fft_x_dim * fft_y_dim
      zx_slab_dbl_cnt = 2 * fft_x_dim * fft_z_dim
  end if

  return

end subroutine fft_setup

!*******************************************************************************
!
! Subroutine:  par_fft_setup
!
! Description:  Setup routine for parallel 3d fft code.  This does the division
!               in an effort to get better overlap between fft slabs and
!               images.  It also allows for more processors than slabs.
!
!               Used only in MPI
!*******************************************************************************

subroutine par_fft_setup(recip_numtasks, fft_slab_redist_enabled)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  ! recip_numtasks returns the initial number of tasks doing recip force calcs.

  integer, intent(out)          :: recip_numtasks
  logical, intent(in)           :: fft_slab_redist_enabled

! Local variables:

  integer               :: alloc_failed

! Allocate space for the task arrays:
 
  if( allocated(xy_slab_cnt)) deallocate(xy_slab_cnt)
  if( allocated(xy_slab_start)) deallocate(xy_slab_start)
  if( allocated(zx_slab_cnt)) deallocate(zx_slab_cnt)
  if( allocated(zx_slab_start)) deallocate(zx_slab_start)
  if( allocated(fft_taskmap)) deallocate(fft_taskmap)

  allocate(xy_slab_cnt(0:numtasks - 1), &
           xy_slab_start(0:numtasks - 1), &
           zx_slab_cnt(0:numtasks - 1), &
           zx_slab_start(0:numtasks - 1), &
!           fft_taskmap(numtasks - 1), &
           fft_taskmap(numtasks), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

! xy_slab_dbl_cnt = total number of allocated spots in an xy slab.  This is
!                   the number of double values that has to be passed per slab
!                   to each task for 2D transforming.
! zx_slab_dbl_cnt = total number of allocated spots in a zx slab.

  xy_slab_dbl_cnt = fft_x_dim * fft_y_dim * 2
  zx_slab_dbl_cnt = fft_x_dim * fft_z_dim * 2

  if (fft_slab_redist_enabled) then
    recip_numtasks = ceiling(dble(numtasks) * fft_workload_estimate)
  else
    recip_numtasks = numtasks ! will be reduced if necessary in following call
  end if

!  recip_numtasks = numtasks ! TEMPORAL_IGOR_FIX
!  write(iunit_debug,*)" par_fft_setup() recip_numtasks = ",recip_numtasks
  call distribute_fft_slabs(recip_numtasks, fft_slab_redist_enabled, &
                            0 )

  return

end subroutine par_fft_setup

!*******************************************************************************
!
! Subroutine:  distribute_fft_slabs
!
! Description:  Distribute the fft slabs workload among the indicated number of
!               reciprocal tasks.
!               
!               Used only in MPI
!
! BUGBUG - UPDATE THE COMMENTS WHEN DONE WITH THE NEW IMPLEMENTATION!
!
!   Each task will get approx. the same number of xy slabs unless tasks
!   outnumber slabs.  Each processor's slabs will be contiguous.
!
!   nfft3  =  total number of xy slabs.
!
!   xy_slab_cnt(taskid) = number of slabs task taskid will have to do.
!
!   my_xy_slab_cnt = number of slabs this task will have to do.
!
!   zx_slab_start(taskid) = last slab of the previous task.
!
!   Same type of variables for the zx slabs:
!   (my_zx_slab_cnt, zx_slab_cnt(0:max), zx_slab_start(0:max))
!
! The grid nfft1 * nfft2 * nfft3 is divided into: 
!
!   x-y slabs nfft1 * nfft2 * xy_slab_cnt(taskid)
!   z-x slabs nfft1 * nfft3 * zx_slab_cnt(taskid)
!
! A subset of the processors get x-y slabs based on a rather complex heuristic
! that should always give slabs to the highest numbered task, should always
! give fewer or the same number of slabs to tasks with lower id's, and should
! always give slabs in chunks of at least max(bspl_order, ceiling(nfft3/nfft2))
! (note - this will change once dynamic fft loadbalancing is implemented!). The
! chunks formula cuts down on the number of places forces/crds must be sent; if
! individual slabs are assigned to a processor, then the processor actually
! processes that slab and the following bspl_order - 1 slabs; by specifying
! slabs in chunks we cut down on the amount of info that must be exchanged
! between tasks.  The ceiling(nfft3/nfft2) term insures that a zx slab can be
! assigned to every task that has an xy slab.  The tasks that have fft slabs are
! distributed roughly equally through the task id space, increasing the overlap
! with their direct force images.
! Note that the number of slabs for task 0 is NOT the maximum number of slabs
! distributed, as it is in the original sander code!
!
! Note that nfft1,2,3 above refer to logical grid sizes; the actual grids are
! dimensioned fft_[x,y,z]_dim, with an assumption that the leftmost index is
! a 2, allowing complex values to be represented as a pair of dbl prec values.
!*******************************************************************************

subroutine distribute_fft_slabs(recip_numtasks, fft_slab_redist_enabled, &
                                loadbal_step_ctr)

  use file_io_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use parallel_mod

  implicit none

! Formal arguments:

  integer, intent(in out)       :: recip_numtasks
  logical, intent(in)           :: fft_slab_redist_enabled
  integer, intent(in)           :: loadbal_step_ctr

! Local variables:

  integer               :: i
  integer               :: next_slab_start
  integer               :: slabs_left
  integer               :: taskid
  integer, save         :: distrib_cnt = 0
  
  distrib_cnt = distrib_cnt + 1 ! count of calls to this routine

  if (recip_numtasks .eq. 0) recip_numtasks = 1

  recip_numtasks = min(recip_numtasks, numtasks, fft_y_dim, fft_z_dim)

  xy_slab_cnt(0:numtasks - 1) = 0

  do i = 0, recip_numtasks - 1
    taskid = (numtasks - 1) - int(dble(i)*dble(numtasks)/dble(recip_numtasks))
    xy_slab_cnt(taskid) = 1
  end do

  slabs_left = fft_z_dim - recip_numtasks

  do
    if (slabs_left .eq. 0) exit
    do i = numtasks - 1, 0, -1
      if (slabs_left .eq. 0) exit
      if (xy_slab_cnt(i) .gt. 0) then
        xy_slab_cnt(i) = xy_slab_cnt(i) + 1
        slabs_left = slabs_left - 1
      end if
    end do
  end do

  next_slab_start = 0

  do i = 0, numtasks - 1
    if (xy_slab_cnt(i) .ne. 0) then
      xy_slab_start(i) = next_slab_start
      next_slab_start = next_slab_start + xy_slab_cnt(i)
    else
      xy_slab_start(i) = -1   ! Catch the bugs...
    end if
  end do

  my_xy_slab_cnt = xy_slab_cnt(mytaskid)
  my_xy_slab_start = xy_slab_start(mytaskid)
  
  ! Now assign zx slabs from the bottom up to processes that got 1 or
  ! more xy slabs.

  zx_slab_cnt(0:numtasks - 1) = 0

  slabs_left = fft_y_dim

  do
    if (slabs_left .eq. 0) exit
    do i = 0, numtasks - 1
      if (slabs_left .eq. 0) exit
      if (xy_slab_cnt(i) .gt. 0) then
        zx_slab_cnt(i) = zx_slab_cnt(i) + 1
        slabs_left = slabs_left - 1
      end if
    end do
  end do

  next_slab_start = 0

  do i = 0, numtasks - 1
    if (zx_slab_cnt(i) .ne. 0) then
      zx_slab_start(i) = next_slab_start
      next_slab_start = next_slab_start + zx_slab_cnt(i)
    else
      zx_slab_start(i) = -1   ! Catch the bugs...
    end if
  end do

  my_zx_slab_cnt = zx_slab_cnt(mytaskid)
  my_zx_slab_start = zx_slab_start(mytaskid)

  ! The slab counts are either both zero or both nonzero.

  i_do_recip = (my_xy_slab_cnt .gt. 0)

  max_xy_slab_cnt = xy_slab_cnt(0)
  max_zx_slab_cnt = zx_slab_cnt(0)
  
  do i = 1, numtasks - 1
    max_xy_slab_cnt = max(max_xy_slab_cnt, xy_slab_cnt(i))
    max_zx_slab_cnt = max(max_zx_slab_cnt, zx_slab_cnt(i))
  end do

  ! Set up the fft taskmap.  This excludes any tasks that don't have slabs.

  fft_taskmap_entries = 0
  fft_taskmap(:) = - 1  ! Force bugs.

  do i = 1, numtasks - 1
    taskid = gbl_taskmap(i)
    if (xy_slab_cnt(taskid) .ne. 0) then
      fft_taskmap_entries = fft_taskmap_entries + 1
      fft_taskmap(fft_taskmap_entries) = taskid
    end if
  end do

  ! Allocate space for multidimensional fft transposes.

  siz_fft_mpi_buf = 2 * fft_x_dim * max_xy_slab_cnt * max_zx_slab_cnt

#ifdef SLOW_NONBLOCKING_MPI
  call set_minimum_mpi_bufs_size(siz_fft_mpi_buf)
#else
  call set_minimum_mpi_bufs_size(siz_fft_mpi_buf * (numtasks - 1))
#endif

  if (master .and. idebug_flag .gt. 0) then

10 format(/, a, /)
12 format(/, a, i5, a, i7, a /)
20 format('  ', a, i4, a)
30 format('  ', a)
40 format('    ', 16i4)

    if( loadbal_verbose .le. 0 ) then
        !
    else if( loadbal_verbose .eq. 1 ) then

      if (distrib_cnt .eq. 1 .and. .not. fft_slab_redist_enabled) then

        ! Initial distribution info for run where dynamic slab distribution
        ! is turned off.

        write(iunit_debug, 10) 'Static FFT Slab Distribution:'
        write(iunit_debug, 20) 'FFT slabs assigned to ', recip_numtasks, ' tasks'
        write(iunit_debug, 20) 'Maximum of ', max_xy_slab_cnt, ' xy slabs per task'
        write(iunit_debug, 20) 'Maximum of ', max_zx_slab_cnt, ' zx slabs per task'
        write(iunit_debug, 30) 'Count of FFT xy slabs assigned to each task:'
        write(iunit_debug, 40) xy_slab_cnt(:)
        write(iunit_debug, 30) 'Count of FFT xz slabs assigned to each task:'
        write(iunit_debug, 40) zx_slab_cnt(:)
         
      else if (distrib_cnt .le. 2 .and. fft_slab_redist_enabled) then

        ! Initial distribution info for run where dynamic slab distribution
        ! is turned on; data is for first two dynamic distributions based on 
        ! the initial workload guess followed by the first actual workload
        ! data.

        if (distrib_cnt .eq. 1) then
          write(iunit_debug, 10) &
            'Initial FFT Slab Distribution Based on Workload Estimate:'
        else
          write(iunit_debug, 10) &
            'First FFT Slab Distribution Based on Actual Workload:'
        end if
        write(iunit_debug, 20) 'FFT slabs assigned to ', recip_numtasks, ' tasks'
        write(iunit_debug, 20) 'Maximum of ', max_xy_slab_cnt, ' xy slabs per task'
        write(iunit_debug, 20) 'Maximum of ', max_zx_slab_cnt, ' zx slabs per task'
        write(iunit_debug, 30) 'Count of FFT xy slabs assigned to each task:'
        write(iunit_debug, 40) xy_slab_cnt(:)
        write(iunit_debug, 30) 'Count of FFT xz slabs assigned to each task:'
        write(iunit_debug, 40) zx_slab_cnt(:)

      end if

    else

      ! Verbose mode.  Report fft slab redistribution for every call...

        write(iunit_debug, 12) &
          'FFT Slab Distribution No. ', distrib_cnt, ' at run step ', &
          loadbal_step_ctr, ':'
        write(iunit_debug, 20) 'FFT slabs assigned to ', recip_numtasks, ' tasks'
        write(iunit_debug, 20) 'Maximum of ', max_xy_slab_cnt, ' xy slabs per task'
        write(iunit_debug, 20) 'Maximum of ', max_zx_slab_cnt, ' zx slabs per task'
        write(iunit_debug, 30) 'Count of FFT xy slabs assigned to each task:'
        write(iunit_debug, 40) xy_slab_cnt(:)
        write(iunit_debug, 30) 'Count of FFT xz slabs assigned to each task:'
        write(iunit_debug, 40) zx_slab_cnt(:)

    end if

  end if

  return

end subroutine distribute_fft_slabs

!*******************************************************************************
!
! Subroutine:  fft3drc_forward
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine fft3drc_forward(xyz_data, zxy_data, x_dim, y_dim, z_dim)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer                               :: x_dim, y_dim, z_dim
  double precision, intent(in out)      :: xyz_data(2, x_dim, y_dim, *)
  double precision, intent(out)         :: zxy_data(2, z_dim, x_dim, *)

! Local variables:

  integer               :: i, j, k

! Each task should do their 2D ffts now:

  do k = 1, my_xy_slab_cnt
    call fft2drc_forward(x_dim, y_dim, xyz_data(1, 1, 1, k))
  end do

! Begin transpose...

  if( numtasks .gt. 1) then
#ifdef SLOW_NONBLOCKING_MPI
      call dist_xyz_zxy_transpose(x_dim, y_dim, z_dim, xyz_data, zxy_data)
#else
      call dist_xyz_zxy_transpose(x_dim, y_dim, z_dim, xyz_data, zxy_data, &
                              fft_taskmap, dbl_mpi_send_buf, dbl_mpi_recv_buf)
#endif
  else
      call xyz_zxy_transpose(x_dim, y_dim, z_dim, xyz_data, zxy_data)
  endif

! End transpose...

! Now do z-fft's:

  do k = 1, my_zx_slab_cnt
    do j = 1, x_dim
      call fft1d_forward(fft_z_hdl, zxy_data(1, 1, j, k))
    end do
  end do

! Results returned in zxy_data.

  return

end subroutine fft3drc_forward

!*******************************************************************************
!
! Subroutine:  xyz_zxy_transpose
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine xyz_zxy_transpose(x_dim, y_dim, z_dim, xyz_data, zxy_data)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: x_dim
  integer               :: y_dim
  integer               :: z_dim
  double precision      :: xyz_data(2 * x_dim * y_dim * z_dim)
  double precision      :: zxy_data(2 * z_dim * x_dim * y_dim)

! Local variables:

  integer               :: i, j, k
  integer               :: xyz_j_off, xyz_k_off
  integer               :: zxy_j_off, zxy_k_off
  integer               :: x_dim_2, z_dim_2

! What we are actually doing, but more efficiently, is:
!
! do j = 1, y_dim
!   do k = 1, z_dim
!     do i = 1, x_dim
!       zxy_data(:, k, i, j) = xyz_data(:, i, j, k)
!     end do
!   end do
! end do

  x_dim_2 = x_dim * 2
  z_dim_2 = z_dim * 2

  if( numtasks .gt. 1 ) then
      do k = 0, my_zx_slab_cnt - 1
        zxy_k_off = k * zx_slab_dbl_cnt + 2 * my_xy_slab_start + 1
        xyz_k_off = (my_zx_slab_start + k) * x_dim_2 + 1
        do j = 0, my_xy_slab_cnt - 1
          zxy_j_off = zxy_k_off + j * 2
          xyz_j_off = xyz_k_off + j * xy_slab_dbl_cnt
          do i = 0, x_dim - 1
            zxy_data(zxy_j_off + i*z_dim_2) = xyz_data(xyz_j_off + i*2)
            zxy_data(zxy_j_off + i*z_dim_2 + 1) = xyz_data(xyz_j_off + i*2 + 1)
          end do
        end do
      end do
  else
      do k = 0, y_dim - 1
        zxy_k_off = k * zx_slab_dbl_cnt + 1
        xyz_k_off = k * x_dim_2 + 1
        do j = 0, z_dim - 1
          zxy_j_off = zxy_k_off + j * 2
          xyz_j_off = xyz_k_off + j * xy_slab_dbl_cnt
          do i = 0, x_dim - 1
            zxy_data(zxy_j_off + i*z_dim_2) = xyz_data(xyz_j_off + i*2)
            zxy_data(zxy_j_off + i*z_dim_2 + 1) = xyz_data(xyz_j_off + i*2 + 1)
          end do
        end do
      end do
  endif

  return

end subroutine xyz_zxy_transpose

#ifdef SLOW_NONBLOCKING_MPI

! This is an inferior implementation for systems that seem unable to handle
! fully async transposes with good i/o overlap:

!*******************************************************************************
!
! Subroutine:  dist_xyz_zxy_transpose
!              used only in MPI
! Description:  <TBS>
!
!*******************************************************************************

subroutine dist_xyz_zxy_transpose(x_dim, y_dim, z_dim, xyz_data, zxy_data)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: x_dim, y_dim, z_dim
  double precision      :: xyz_data(*), zxy_data(*)

! Local variables:

  integer               :: i, j, k
  integer               :: xyz_j_off, xyz_k_off
  integer               :: buf_off
  integer               :: task_ctr, rtask, stask
  integer               :: ireq, isnd_stat(mpi_status_size)

  rtask = mytaskid

  do task_ctr = mytaskid + 1, mytaskid + numtasks - 1

    stask = mod(task_ctr, numtasks)
    if (stask .eq. mytaskid) cycle

    buf_off = 1

    do k = 0, zx_slab_cnt(stask) - 1
      xyz_k_off = (zx_slab_start(stask) + k) *  x_dim * 2 + 1
      do j = 0, my_xy_slab_cnt - 1
        xyz_j_off = xyz_k_off + j * xy_slab_dbl_cnt
        do i = 0, x_dim - 1
          dbl_mpi_send_buf(buf_off + i*2) = xyz_data(xyz_j_off + i*2)
          dbl_mpi_send_buf(buf_off + i*2 + 1) = xyz_data(xyz_j_off + i*2 + 1)
        end do
        buf_off = buf_off + x_dim * 2
      end do
    end do

    if (buf_off .gt. 1) then
      call mpi_isend(dbl_mpi_send_buf, buf_off - 1, mpi_double_precision, &
                     stask, xyzxt_tag, lib_mpi_comm, ireq, err_code_mpi)
    end if

    rtask = rtask - 1
    if (rtask .lt. 0) rtask = rtask + numtasks

    call xyz_zxy_trans_receive(x_dim, z_dim, zxy_data, rtask)

    if (buf_off .gt. 1) then
      call mpi_wait(ireq, isnd_stat, err_code_mpi)
    end if

  end do

  ! Now do your own transposes.

  call xyz_zxy_transpose(x_dim, y_dim, z_dim, xyz_data, zxy_data)

  return

end subroutine dist_xyz_zxy_transpose

!*******************************************************************************
!
! Subroutine:  xyz_zxy_trans_receive
!
! Description: SLOW_NONBLOCKING_MPI implementation!
!              
!*******************************************************************************

subroutine xyz_zxy_trans_receive(x_dim, z_dim, zxy_data, rtask)

  use parallel_dat_mod

  implicit none

  integer               :: x_dim, z_dim
  double precision      :: zxy_data(*)
  integer               :: rtask

  integer               :: i, j, k
  integer               :: zxy_j_off, zxy_k_off
  integer               :: numval, buf_off
  integer               :: z_dim_2
  integer               :: status(mpi_status_size)

  if (xy_slab_cnt(rtask) .eq. 0) return

  numval = 2 * x_dim * my_zx_slab_cnt * xy_slab_cnt(rtask)

  call mpi_recv(dbl_mpi_recv_buf, numval, mpi_double_precision, &
                rtask, xyzxt_tag, lib_mpi_comm, status, err_code_mpi) 

  z_dim_2 = z_dim * 2
  buf_off = 1
  do k = 0, my_zx_slab_cnt - 1
    zxy_k_off = k * zx_slab_dbl_cnt + 2*xy_slab_start(rtask)+1
    do j = 0, xy_slab_cnt(rtask) - 1
      zxy_j_off = zxy_k_off + j * 2
      do i = 0, x_dim - 1
        zxy_data(zxy_j_off + i*z_dim_2) = dbl_mpi_recv_buf(buf_off + i*2)
        zxy_data(zxy_j_off + i*z_dim_2+1) = dbl_mpi_recv_buf(buf_off + i*2+1)
      end do
      buf_off = buf_off + x_dim * 2
    end do
  end do

  return

end subroutine xyz_zxy_trans_receive

#else
!*******************************************************************************
!
! Subroutine:  dist_xyz_zxy_transpose
!
! Description: Fully asynchronous implementation!
!              
!*******************************************************************************

subroutine dist_xyz_zxy_transpose(x_dim, y_dim, z_dim, xyz_data, zxy_data, &
                                  taskmap, send_buf, recv_buf)
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: x_dim, y_dim, z_dim
  double precision      :: xyz_data(*), zxy_data(*)
  integer               :: taskmap(numtasks - 1)
  double precision      :: send_buf(siz_fft_mpi_buf * fft_taskmap_entries)
  double precision      :: recv_buf(siz_fft_mpi_buf, fft_taskmap_entries)

! Local variables:

  integer               :: i, j, k
  integer               :: xyz_j_off, xyz_k_off
  integer               :: numval, buf_base, buf_off
  integer               :: taskmap_idx, recv_task, send_task
  integer               :: irecv_stat(mpi_status_size)
  integer               :: isend_stat(mpi_status_size, fft_taskmap_entries)
  integer               :: recv_req(fft_taskmap_entries)
  integer               :: send_req(fft_taskmap_entries)

  ! Post the asynchronous receives first:

  do taskmap_idx = 1, fft_taskmap_entries

    recv_task = taskmap(taskmap_idx)
    numval = 2 * x_dim * my_zx_slab_cnt * xy_slab_cnt(recv_task)
    call mpi_irecv(recv_buf(1, taskmap_idx), numval, mpi_double_precision, &
                   recv_task, xyzxt_tag, lib_mpi_comm, &
                   recv_req(taskmap_idx), err_code_mpi) 
  end do

  ! Now set up and post the asynchronous sends:

  buf_base = 1

  do taskmap_idx = 1, fft_taskmap_entries

    send_task = taskmap(taskmap_idx)

    buf_off = buf_base

    do k = 0, zx_slab_cnt(send_task) - 1
      xyz_k_off = (zx_slab_start(send_task) + k) *  x_dim * 2 + 1
      do j = 0, my_xy_slab_cnt - 1
        xyz_j_off = xyz_k_off + j * xy_slab_dbl_cnt
        do i = 0, x_dim - 1
          send_buf(buf_off + i*2) = xyz_data(xyz_j_off + i*2)
          send_buf(buf_off + i*2 + 1) = xyz_data(xyz_j_off + i*2 + 1)
        end do
        buf_off = buf_off + x_dim * 2
      end do
    end do

!   write(0,*)'DBG: mpi_isend of ', numval * 8, ' bytes by task ', mytaskid

    call mpi_isend(send_buf(buf_base), buf_off - buf_base, &
                   mpi_double_precision, send_task, xyzxt_tag, lib_mpi_comm, &
                   send_req(taskmap_idx), err_code_mpi)

    buf_base = buf_base + siz_fft_mpi_buf

  end do

  ! Do your own transposes, overlapped with i/o.

  call xyz_zxy_transpose(x_dim, y_dim, z_dim, xyz_data, zxy_data)

  ! Wait on and process the pending receive requests:

  do j = 1, fft_taskmap_entries

    call mpi_waitany(fft_taskmap_entries, recv_req, taskmap_idx, irecv_stat, &
                     err_code_mpi)
    recv_task = taskmap(taskmap_idx)
    call xyz_zxy_trans_receive(x_dim, z_dim, zxy_data, &
                               recv_buf(1, taskmap_idx), recv_task)

  end do

  ! Wait for all sends to complete:

  call mpi_waitall(fft_taskmap_entries, send_req, isend_stat, err_code_mpi)

  return

end subroutine dist_xyz_zxy_transpose

!*******************************************************************************
!
! Subroutine:  xyz_zxy_trans_receive
!
! Description: Fully asynchronous implementation!
!              
!*******************************************************************************

subroutine xyz_zxy_trans_receive(x_dim, z_dim, zxy_data, recv_buf, recv_task)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: x_dim, z_dim
  double precision      :: zxy_data(*)
  double precision      :: recv_buf(*)
  integer               :: recv_task

! Local variables:

  integer               :: i, j, k
  integer               :: zxy_j_off, zxy_k_off
  integer               :: buf_off
  integer               :: z_dim_2

  z_dim_2 = z_dim * 2
  buf_off = 1
  do k = 0, my_zx_slab_cnt - 1
    zxy_k_off = k * zx_slab_dbl_cnt + 2*xy_slab_start(recv_task)+1
      do j = 0, xy_slab_cnt(recv_task) - 1
      zxy_j_off = zxy_k_off + j * 2
      do i = 0, x_dim - 1
        zxy_data(zxy_j_off + i*z_dim_2) = recv_buf(buf_off + i*2)
        zxy_data(zxy_j_off + i*z_dim_2+1) = recv_buf(buf_off + i*2+1)
      end do
      buf_off = buf_off + x_dim * 2
    end do
  end do

  return

end subroutine xyz_zxy_trans_receive
#endif /* NOT SLOW_NONBLOCKING_MPI */

!*******************************************************************************
!
! Subroutine:  fft2drc_forward
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine fft2drc_forward(x_dim, y_dim, xy_slab)

  implicit none

! Formal arguments:

  integer               :: x_dim
  integer               :: y_dim
  double precision      :: xy_slab(2, x_dim, y_dim)

! Local variables:

  double precision      :: a, b, c, d
  double precision      :: x_run(2, x_dim)
  double precision      :: y_run(2, y_dim)
  integer               :: i, j
  integer               :: x_lim       ! indexing limit

! real to complex:

! First the x direction, the data is already contiguous:

  x_lim = x_dim + 1

  do j = 1, y_dim

    do i = 1,  x_dim - 1
      x_run(:, i) = xy_slab(:, i, j)
    end do

    call fft1d_forward(fft_x_hdl, x_run)

    ! We have to put results in a temporary run to prevent overwriting of
    ! values needed later.

    do i = 2, x_dim - 1
      a =  (x_run(1, i) + x_run(1, x_lim - i)) ! Real F even * 2
      b =  (x_run(2, i) - x_run(2, x_lim - i)) ! Imag F even * 2
      c =  (x_run(2, i) + x_run(2, x_lim - i)) ! Real F odd * 2
      d = -(x_run(1, i) - x_run(1, x_lim - i)) ! Imag F odd * 2
      xy_slab(1, i, j) = .5d0 * (a + rc_coefs(1, i) * c + rc_coefs(2, i) * d)
      xy_slab(2, i, j) = .5d0 * (b + rc_coefs(1, i) * d - rc_coefs(2, i) * c)
    end do

! DC and nyquist:

    xy_slab(1, 1, j) = x_run(1, 1) + x_run(2, 1)
    xy_slab(2, 1, j) = 0.d0
    xy_slab(1, x_dim, j) = x_run(1, 1) - x_run(2, 1)
    xy_slab(2, x_dim, j) = 0.d0

  end do

! Now in the y direction, the data is in y now and we will put it into a
! contiguous 1D array first, transform, then put it back.

  do i = 1, x_dim

    do j = 1, y_dim
      y_run(:, j)   = xy_slab(:, i, j)
    end do            

    call fft1d_forward(fft_y_hdl, y_run)

    do j = 1, y_dim
      xy_slab(:, i, j) = y_run(:, j)
    end do            

  end do

  return
      
end subroutine fft2drc_forward

!*******************************************************************************
!
! Subroutine:  fft3drc_back
!
! Description:  <TBS>
!              
!*******************************************************************************

subroutine fft3drc_back(zxy_data, xyz_data, x_dim, y_dim, z_dim)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer                               :: x_dim, y_dim, z_dim
  double precision, intent(in out)      :: zxy_data(2, z_dim, x_dim, *)
  double precision, intent(out)         :: xyz_data(2, x_dim, y_dim, *)

! Local variables:

  integer               :: i, j, k
  integer               :: ja, ja00, jz, jz00, jidx
  integer               :: k0, ks

! Data enters z-x slabs in zxy_data...

! Do z ffts now:

  do k = 1, my_zx_slab_cnt
    do j = 1, x_dim
      call fft1d_back(fft_z_hdl, zxy_data(1, 1, j, k))
    end do
  end do

! Redistribute into xy slabs:

  if( numtasks .gt. 1) then
#ifdef SLOW_NONBLOCKING_MPI
      call dist_zxy_xyz_transpose(x_dim, y_dim, z_dim, zxy_data, xyz_data)
#else
      call dist_zxy_xyz_transpose(x_dim, y_dim, z_dim, zxy_data, xyz_data, &
                              fft_taskmap, dbl_mpi_send_buf, dbl_mpi_recv_buf)
#endif
  else
      call zxy_xyz_transpose(x_dim, y_dim, z_dim, zxy_data, xyz_data)
  endif

! Each task should do their 2D fft's now:

  do k = 1, my_xy_slab_cnt
    call fft2drc_back(x_dim, y_dim, xyz_data(1, 1, 1, k))
  end do

  return

end subroutine fft3drc_back

!*******************************************************************************
!
! Subroutine:  zxy_xyz_transpose
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine zxy_xyz_transpose(x_dim, y_dim, z_dim, zxy_data, xyz_data)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: x_dim
  integer               :: y_dim
  integer               :: z_dim
  double precision      :: zxy_data(2 * z_dim * x_dim * y_dim)
  double precision      :: xyz_data(2 * x_dim * y_dim * z_dim)

! Local variables:

  integer               :: i, j, k
  integer               :: xyz_j_off,  xyz_k_off
  integer               :: zxy_j_off, zxy_k_off
  integer               :: x_dim_2, z_dim_2

! What we are actually doing, but more efficiently, is:
!
! do j = 1, y_dim
!   do k = 1, z_dim
!     do i = 1, x_dim
!       xyz_data(:, i, j, k) = zxy_data(:, k, i, j)
!     end do
!   end do
! end do
!
  
  x_dim_2 = x_dim * 2
  z_dim_2 = z_dim * 2

  if( numtasks .gt. 1) then
      do k = 0, my_zx_slab_cnt - 1
        xyz_k_off = (my_zx_slab_start + k) * x_dim_2 + 1
        zxy_k_off = k * zx_slab_dbl_cnt + 2 * my_xy_slab_start + 1
        do j = 0, my_xy_slab_cnt - 1
          xyz_j_off = xyz_k_off + j * xy_slab_dbl_cnt
          zxy_j_off = zxy_k_off + j * 2
          do i = 0, x_dim - 1
            xyz_data(xyz_j_off + i*2) = zxy_data(zxy_j_off + i*z_dim_2)
            xyz_data(xyz_j_off + i*2 + 1) = zxy_data(zxy_j_off + i*z_dim_2 + 1)
          end do
        end do
      end do
  else
      do k = 0, y_dim - 1
        xyz_k_off = k * x_dim_2 + 1
        zxy_k_off = k * zx_slab_dbl_cnt + 1
        do j = 0, z_dim - 1
          xyz_j_off = xyz_k_off + j * xy_slab_dbl_cnt
          zxy_j_off = zxy_k_off + j * 2
          do i = 0, x_dim - 1
            xyz_data(xyz_j_off + i*2) = zxy_data(zxy_j_off + i*z_dim_2)
            xyz_data(xyz_j_off + i*2 + 1) = zxy_data(zxy_j_off + i*z_dim_2 + 1)
          end do
        end do
      end do
  end if ! ( numtasks > 1) 

  return

end subroutine zxy_xyz_transpose

#ifdef SLOW_NONBLOCKING_MPI

! This is an inferior implementation for systems that seem unable to handle
! fully async transposes with good i/o overlap:

!*******************************************************************************
!
! Subroutine:  dist_zxy_xyz_transpose
!
! Description: SLOW_NONBLOCKING_MPI implementation!
!              
!*******************************************************************************

subroutine dist_zxy_xyz_transpose(x_dim, y_dim, z_dim, zxy_data, xyz_data)

  use parallel_dat_mod

  implicit none

  ! Formal arguments:

  integer               :: x_dim, y_dim, z_dim
  double precision      :: zxy_data(*), xyz_data(*)

  ! Local variables:

  integer               :: i, j, k
  integer               :: zxy_j_off, zxy_k_off
  integer               :: buf_off
  integer               :: z_dim_2
  integer               :: task_ctr, rtask, stask
  integer               :: ireq, isnd_stat(mpi_status_size)

  rtask = mytaskid

  do task_ctr = mytaskid + 1, mytaskid + numtasks - 1

    stask = mod(task_ctr, numtasks)
    if (stask .eq. mytaskid) cycle

    z_dim_2 = z_dim * 2
    buf_off = 1

    do k = 0, my_zx_slab_cnt - 1
      zxy_k_off = k * zx_slab_dbl_cnt + 2 * xy_slab_start(stask) + 1
      do j = 0, xy_slab_cnt(stask) - 1
        zxy_j_off = zxy_k_off + j * 2
        do i = 0, x_dim - 1
          dbl_mpi_send_buf(buf_off + i*2) = zxy_data(zxy_j_off + i*z_dim_2)
          dbl_mpi_send_buf(buf_off + i*2+1) = zxy_data(zxy_j_off + i*z_dim_2+1)
        end do
        buf_off = buf_off + x_dim * 2
      end do
    end do

    if (buf_off .gt. 1) then
      call mpi_isend(dbl_mpi_send_buf, buf_off - 1, mpi_double_precision, &
                     stask, zxxyt_tag, lib_mpi_comm, ireq, err_code_mpi)
    end if

    rtask = rtask - 1
    if (rtask .lt. 0) rtask = rtask + numtasks

    call zxy_xyz_trans_receive(x_dim, z_dim, xyz_data, rtask)

    if (buf_off .gt. 1) then
      call mpi_wait(ireq, isnd_stat, err_code_mpi)
    end if

  end do

  ! Now do your own transposes.

  call zxy_xyz_transpose(x_dim, y_dim, z_dim, zxy_data, xyz_data)

  return

end subroutine dist_zxy_xyz_transpose

!*******************************************************************************
!
! Subroutine:  zxy_xyz_trans_receive
!
! Description: SLOW_NONBLOCKING_MPI implementation!
!              
!*******************************************************************************

subroutine zxy_xyz_trans_receive(x_dim, z_dim, xyz_data, rtask)

  use parallel_dat_mod

  implicit none

  integer               :: x_dim, z_dim
  double precision      :: xyz_data(*)
  integer               :: rtask

  integer               :: i, j, k
  integer               :: xyz_j_off, xyz_k_off
  integer               :: numval, buf_off
  integer               :: status(mpi_status_size)

  if (zx_slab_cnt(rtask) .eq. 0) return

  numval = 2 * x_dim * my_xy_slab_cnt * zx_slab_cnt(rtask)

  call mpi_recv(dbl_mpi_recv_buf, numval, mpi_double_precision, rtask, &
                zxxyt_tag, lib_mpi_comm, status, err_code_mpi) 

  buf_off = 1
  do k = 0, zx_slab_cnt(rtask) - 1
    xyz_k_off = (zx_slab_start(rtask) + k) * x_dim*2 + 1
    do j = 0, my_xy_slab_cnt - 1
      xyz_j_off = xyz_k_off + j * xy_slab_dbl_cnt
      do i = 0, x_dim - 1
        xyz_data(xyz_j_off + i*2) = dbl_mpi_recv_buf(buf_off + i*2)
        xyz_data(xyz_j_off + i*2+1) = dbl_mpi_recv_buf(buf_off + i*2+1)
      end do
      buf_off = buf_off + x_dim * 2
    end do
  end do

  return

end subroutine zxy_xyz_trans_receive

#else
!*******************************************************************************
!
! Subroutine:  dist_zxy_xyz_transpose
!
! Description: Fully asynchronous implementation!
!              
!*******************************************************************************

subroutine dist_zxy_xyz_transpose(x_dim, y_dim, z_dim, zxy_data, xyz_data, &
                                  taskmap, send_buf, recv_buf)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: x_dim, y_dim, z_dim
  double precision      :: zxy_data(*), xyz_data(*)
  integer               :: taskmap(numtasks - 1)
  double precision      :: send_buf(siz_fft_mpi_buf * fft_taskmap_entries)
  double precision      :: recv_buf(siz_fft_mpi_buf, fft_taskmap_entries)

! Local variables:

  integer               :: i, j, k
  integer               :: zxy_j_off, zxy_k_off
  integer               :: numval, buf_base, buf_off
  integer               :: z_dim_2
  integer               :: taskmap_idx, recv_task, send_task
  integer               :: irecv_stat(mpi_status_size)
  integer               :: isend_stat(mpi_status_size, fft_taskmap_entries)
  integer               :: recv_req(fft_taskmap_entries)
  integer               :: send_req(fft_taskmap_entries)

  ! Post the asynchronous receives first:

  do taskmap_idx = 1, fft_taskmap_entries

    recv_task = taskmap(taskmap_idx)
    numval = 2 * x_dim * my_xy_slab_cnt * zx_slab_cnt(recv_task)
    call mpi_irecv(recv_buf(1, taskmap_idx), numval, mpi_double_precision, &
                   recv_task, zxxyt_tag, lib_mpi_comm, &
                   recv_req(taskmap_idx), err_code_mpi)
  end do

  ! Now set up and post the asynchronous sends:

  buf_base = 1

  do taskmap_idx = 1, fft_taskmap_entries

    send_task = taskmap(taskmap_idx)

    z_dim_2 = z_dim * 2
    buf_off = buf_base

    do k = 0, my_zx_slab_cnt - 1
      zxy_k_off = k * zx_slab_dbl_cnt + 2 * xy_slab_start(send_task) + 1
      do j = 0, xy_slab_cnt(send_task) - 1
        zxy_j_off = zxy_k_off + j * 2
        do i = 0, x_dim - 1
          send_buf(buf_off + i*2) = zxy_data(zxy_j_off + i*z_dim_2)
          send_buf(buf_off + i*2+1) = zxy_data(zxy_j_off + i*z_dim_2+1)
        end do
        buf_off = buf_off + x_dim * 2
      end do
    end do

!   write(0,*)'DBG: mpi_isend of ', numval * 8, ' bytes by task ', mytaskid

    call mpi_isend(send_buf(buf_base), buf_off - buf_base, &
                   mpi_double_precision, send_task, zxxyt_tag, lib_mpi_comm, &
                   send_req(taskmap_idx), err_code_mpi)

    buf_base = buf_base + siz_fft_mpi_buf

  end do

  ! Do your own transposes, overlapped with i/o.

  call zxy_xyz_transpose(x_dim, y_dim, z_dim, zxy_data, xyz_data)

! Wait on and process the pending receive requests:

  do j = 1, fft_taskmap_entries

    call mpi_waitany(fft_taskmap_entries, recv_req, taskmap_idx, irecv_stat, &
                     err_code_mpi)
    recv_task = taskmap(taskmap_idx)
    call zxy_xyz_trans_receive(x_dim, z_dim, xyz_data, &
                               recv_buf(1, taskmap_idx), recv_task)

  end do

  ! Wait for all sends to complete:

  call mpi_waitall(fft_taskmap_entries, send_req, isend_stat, err_code_mpi)

  return

end subroutine dist_zxy_xyz_transpose

!*******************************************************************************
!
! Subroutine:  zxy_xyz_trans_receive
!
! Description: Fully asynchronous implementation!
!              
!*******************************************************************************

subroutine zxy_xyz_trans_receive(x_dim, z_dim, xyz_data, recv_buf, recv_task)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: x_dim, z_dim
  double precision      :: xyz_data(*)
  double precision      :: recv_buf(*)
  integer               :: recv_task

! Local variables:

  integer               :: i, j, k
  integer               :: xyz_j_off, xyz_k_off
  integer               :: buf_off

  buf_off = 1
  do k = 0, zx_slab_cnt(recv_task) - 1
    xyz_k_off = (zx_slab_start(recv_task) + k) * x_dim*2 + 1
    do j = 0, my_xy_slab_cnt - 1
      xyz_j_off = xyz_k_off + j * xy_slab_dbl_cnt
      do i = 0, x_dim - 1
        xyz_data(xyz_j_off + i*2) = recv_buf(buf_off + i*2)
        xyz_data(xyz_j_off + i*2+1) = recv_buf(buf_off + i*2+1)
      end do
      buf_off = buf_off + x_dim * 2
    end do
  end do

  return

end subroutine zxy_xyz_trans_receive
#endif /* NOT SLOW_NONBLOCKING_MPI */

!*******************************************************************************
!
! Subroutine:  fft2drc_back
!
! Description:
!
!*******************************************************************************

subroutine fft2drc_back(x_dim, y_dim, xy_slab)

  implicit none

! Formal arguments:

  integer               :: x_dim
  integer               :: y_dim
  double precision      :: xy_slab(2, x_dim, y_dim) ! an xy slab

! Local variables:

  double precision      :: a, b, c, d
  double precision      :: x_run(2, x_dim)    ! run of all x's for one y
  double precision      :: y_run(2, y_dim)    ! run of all y's for one x
  integer               :: i, j
  integer               :: x_lim              ! indexing limit

! complex to real:

! Now in the y direction, the data is in y now and we will put it into a
! contiguous 1D array first, transform, then put it back.

  x_lim = x_dim + 1

  do i = 1, x_dim

    do j = 1, y_dim
      y_run(:, j) = xy_slab(:, i, j)
    end do            

    call fft1d_back(fft_y_hdl, y_run)

    do j = 1, y_dim
      xy_slab(:, i, j) = y_run(:, j)
    end do            

  end do

  do j = 1, y_dim

    ! We have to put results in a temporary run to prevent overwriting of
    ! values needed later.

    do i = 2, x_dim - 1
      a = (xy_slab(1, i, j) + xy_slab(1, x_lim - i, j)) ! Real F even
      b = (xy_slab(2, i, j) - xy_slab(2, x_lim - i, j)) ! Imag F even
      c = (xy_slab(2, i, j) + xy_slab(2, x_lim - i, j)) ! F odd contrib
      d = (xy_slab(1, i, j) - xy_slab(1, x_lim - i, j)) ! F odd contrib
      x_run(1, i) = a - rc_coefs(1, i) * c - rc_coefs(2, i) * d
      x_run(2, i) = b + rc_coefs(1, i) * d - rc_coefs(2, i) * c
    end do

    x_run(1, 1) = xy_slab(1, 1, j) + xy_slab(1, x_dim, j)
    x_run(2, 1) = xy_slab(1, 1, j) - xy_slab(1, x_dim, j)

    call fft1d_back(fft_x_hdl, x_run)

    do i = 1, x_dim - 1
      xy_slab(:, i, j) = x_run(:, i)
    end do

  end do

  return

end subroutine fft2drc_back

end module pme_fft_mod
