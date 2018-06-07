#include "copyright.i"

!*******************************************************************************
!
! Module: parallel_dat_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module parallel_dat_mod

#ifdef INCLUDE_MPIF
  implicit none
# include <mpif.h>  
#else
  use mpi
  implicit none
#endif

  logical, save         :: master      ! Is this task the MPI master task?


! MAJOR NOTE:  Much of the stuff here is basically unnecessary under
!              Generalized Born; due to large cutoffs and small atom counts
!              in that environment, spatial decomposition is somewhat
!              impractical.  We DO use some of the data structures; these are
!              noted by a "VGB" sidebar - meaning "valid for GB".

! Global variables necessary for MPI implementation. These are NOT broadcast.

  integer, save         :: mytaskid     !  MPI rank of the process ( = 0 for master)
  integer, save         :: numtasks     !  number of MPI processes
  integer, save         :: lib_mpi_group   ! MPI processor group on which library calls are executed
  integer, save         :: lib_mpi_comm    ! MPI communicator on which library calls are executed
  integer, save         :: err_code_mpi
  integer, save         :: notdone                                         ! VGB
  logical, save         :: i_do_recip
  
  integer, save         :: iunit_debug     !  fortran unit to save debug information
  integer, save         :: idebug_flag = 1 !  flag to output debug info   

  integer, save         :: my_send_atm_cnts_total = 0
  integer, save         :: my_send_atm_cnts_sums = 0

  integer, parameter    :: dif_tag = 10
  integer, parameter    :: dc_tag = 11
  integer, parameter    :: gifd_tag = 12
  integer, parameter    :: xyzxt_tag = 13
  integer, parameter    :: zxxyt_tag = 14
  

  ! Used for async mpi sends/recvs in fft and force code.

  integer, allocatable, save    :: gbl_taskmap(:)
  integer, allocatable, save    :: gbl_inv_taskmap(:)

  integer, allocatable, save    :: gbl_img_div_tbl(:)  ! distribution of atoms on processors
                                                       ! idx from img_lo = gbl_img_div_tbl(task_id) + 1
                                                       ! img_hi = gbl_img_div_tbl(task_id + 1) 

  ! The following three arrays are for use in mpi routines defined in this
  ! module; the vectors they are used on (crd,frc,vel) DO NOT
  ! contain data that is segregated by owning task (ie., these are not
  ! really division tables).  Instead, data must be marshalled into and out of
  ! temporary buffers according to gbl_atm_owner_map().

  integer, save         :: extra_used_atm_cnt

  integer, allocatable, save    :: gbl_atm_offsets(:)      ! VGB
  integer, allocatable, save    :: gbl_vec_offsets(:)      ! VGB
  integer, allocatable, save    :: gbl_vec_rcvcnts(:)      ! allocated(0:numtasks) - number of coordiates per MPI process (task) ( num atoms per MPI process X 3)  

!  integer, allocatable, save    :: gbl_atm_owner_map(:)    ! allocated(natom), = taskid (rank) atom belong to  ( = 0 for serial run)

  integer, allocatable, save    :: gbl_send_atm_lst(:)
  integer, allocatable, save    :: gbl_send_atm_cnts(:)
  integer, allocatable, save    :: gbl_recv_atm_lsts(:,:)
  integer, allocatable, save    :: gbl_recv_atm_cnts(:)
  integer, allocatable, save    :: gbl_extra_used_atms(:)   ! Atom used when computing energies and forces but not owned by the processor 

  integer, save :: siz_dbl_mpi_bufs = 0                                     !VGB

  double precision, allocatable, save   :: dbl_mpi_send_buf(:)              !VGB
  double precision, allocatable, save   :: dbl_mpi_recv_buf(:)              !VGB

end module parallel_dat_mod

