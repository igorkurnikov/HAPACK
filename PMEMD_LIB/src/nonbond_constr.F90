#include "copyright.i"

!*******************************************************************************
!
! Module:  dist_constr_mod
!
! Description: Atom - Atom Distance constraints 
!              
!*******************************************************************************

module dist_constr_mod

  use gbl_datatypes_mod

  implicit none

! The following are derived from prmtop non-bond constrains info:

  integer, save                          :: cit_num_dist_constr
  
  type(dist_constr_rec), allocatable, save     :: cit_dist_constr(:)

contains

!*******************************************************************************
!
! Subroutine:  dist_constr_alloc
!
! Description:  setup atom-atom distance constraints
!
!*******************************************************************************
subroutine dist_constr_alloc(num_constr)

  use prmtop_dat_mod

implicit none

! Formal arguments:

  integer                       :: num_constr

! Local variables:

  integer               :: alloc_failed 
  
   if( allocated(gbl_dist_constr)) deallocate(gbl_dist_constr)

   allocate( gbl_dist_constr(num_constr),   &
             stat = alloc_failed )

   if (alloc_failed .ne. 0) call setup_alloc_error

return

end subroutine dist_constr_alloc


!*******************************************************************************
!
! Subroutine:  dist_constr_setup
!
! Description:  setup atom-atom distance constraints
!
!*******************************************************************************

subroutine dist_constr_setup(use_atm_map, atm_owner_map)

   use parallel_dat_mod
   use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer                       :: use_atm_map(natom)
  integer                       :: atm_owner_map(natom)

! Local variables:

  integer               :: alloc_failed
  integer               :: ic
  
  type(dist_constr_rec)   :: dist_constr_copy(num_dist_constr)

  ! This routine can handle reallocation, and thus can be called multiple
  ! times.

!  write(iunit_debug,*) " dist_constr_setup()  num_dist_constr = ",num_dist_constr

!  write(*,*) "  "
!  write(*,*) " dist_constr_setup() pt 1 num_dist_constr = ",num_dist_constr
!  write(*,*) "  "

  call find_my_dist_constr(num_dist_constr, gbl_dist_constr, cit_num_dist_constr, dist_constr_copy, use_atm_map, atm_owner_map)
  
!  do ic = 1,num_dist_constr
!     write(*,*) " dist_constr_setup() pt 2 gbl_dist_constr(ic)%atm_i ", ic, gbl_dist_constr(ic)%atm_i 
!     write(*,*) " dist_constr_setup() pt 2 gbl_dist_constr(ic)%atm_j ", ic, gbl_dist_constr(ic)%atm_j 
!     write(*,*) "  "
!     write(*,*) " dist_constr_setup() pt 2 dist_constr_copy(ic)%atm_i ", ic, dist_constr_copy(ic)%atm_i 
!     write(*,*) " dist_constr_setup() pt 2 dist_constr_copy(ic)%atm_i ", ic, dist_constr_copy(ic)%atm_j
!     write(*,*) "  "
!  enddo
  
  if (cit_num_dist_constr .gt. 0) then
    if (allocated(cit_dist_constr)) then
      if (size(cit_dist_constr) .lt. cit_num_dist_constr) then
        deallocate(cit_dist_constr)
        allocate(cit_dist_constr(cit_num_dist_constr), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
      end if
    else
      allocate(cit_dist_constr(cit_num_dist_constr), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
    end if
    cit_dist_constr(1:cit_num_dist_constr) = dist_constr_copy(1:cit_num_dist_constr)
  end if

  return

end subroutine dist_constr_setup

!*******************************************************************************
!
! Subroutine:  find_my_dist_constraints
!
! Description:  Set arrays of atom-atom distance constraints that belong to the process
!
!*******************************************************************************

subroutine find_my_dist_constr(dist_constr_cnt, dist_constr, my_dist_constr_cnt, my_dist_constr, &
                             use_atm_map, atm_owner_map)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: dist_constr_cnt
  type(dist_constr_rec)   :: dist_constr(dist_constr_cnt)
  integer               :: my_dist_constr_cnt
  type(dist_constr_rec)   :: my_dist_constr(*)
  integer               :: use_atm_map(*)
  integer               :: atm_owner_map(*)

! Local variables:

  integer               :: atm_i, atm_j, dist_constr_idx

! Find all nonbond constrains for which this process owns either atom:

  my_dist_constr_cnt = 0
   
  do dist_constr_idx = 1, dist_constr_cnt

    atm_i = dist_constr(dist_constr_idx)%atm_i
    atm_j = dist_constr(dist_constr_idx)%atm_j

    if ( (atm_owner_map(atm_i) .eq. mytaskid) .or. (atm_owner_map(atm_j) .eq. mytaskid) ) then
      my_dist_constr_cnt = my_dist_constr_cnt + 1
      my_dist_constr(my_dist_constr_cnt) = dist_constr(dist_constr_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      
!      write(iunit_debug,*)" find_my_dist_constr() my_dist_constr_cnt = ", my_dist_constr_cnt
!      write(iunit_debug,*)" my_dist_constr(my_dist_constr_cnt)%atm_i =",my_dist_constr(my_dist_constr_cnt)%atm_i
!      write(iunit_debug,*)" my_dist_constr(my_dist_constr_cnt)%atm_j =",my_dist_constr(my_dist_constr_cnt)%atm_j
!      write(iunit_debug,*)" my_dist_constr(my_dist_constr_cnt)%acoef =",my_dist_constr(my_dist_constr_cnt)%acoef
!      write(iunit_debug,*)" my_dist_constr(my_dist_constr_cnt)%bcoef =",my_dist_constr(my_dist_constr_cnt)%bcoef
!      write(iunit_debug,*)" my_dist_constr(my_dist_constr_cnt)%cnt_type =",my_dist_constr(my_dist_constr_cnt)%cnt_type
!      write(iunit_debug,*)" "
           
    end if
  end do
  
!  write(iunit_debug,*) " find_my_dist_constr()  my_dist_constr_cnt = ",my_dist_constr_cnt

  return

end subroutine find_my_dist_constr

!*******************************************************************************
!
! Subroutine:  get_dist_constr_energy
!
! Description:
!              
! Routine to compute atom-atom distance constraint energy:
! Harmonic constraints (type = ??)
!
! non-bonded constraints (type = ??) energy and forces for 
! nonbonded constrains of the type E = A/R^12 - B/R^6  and E = A/R^12 - B/R^10
!
!*******************************************************************************

subroutine get_dist_constr_energy(dist_constr_cnt, dist_constr, x, frc, dist_constr_energy, atm_owner_map)

  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: dist_constr_cnt
  type(dist_constr_rec) :: dist_constr(*)
  double precision      :: x(3, *)
  double precision      :: frc(3, *)
  double precision      :: dist_constr_energy
  integer               :: atm_owner_map(*)

! Local variables:

  double precision      :: da, df, de, evdw
  integer               :: i, j, ic, jn
  double precision      :: lcl_nb_constr_energy
  double precision      :: xa, ya, za
  double precision      :: rij
  double precision      :: xij, yij, zij
  double precision      :: rinv, r2inv, r6inv, r10inv
  double precision      :: f2, f6, f12, f10      
  double precision      :: r06inv           

  double precision      :: ene_constr

  lcl_nb_constr_energy = 0.0d0

! Grand loop for the atom-atom distance constraints:

!  write(*,*)"  get_dist_constr_energy() pt 1   dist_constr_cnt =", dist_constr_cnt
  
  do jn = 1, dist_constr_cnt

    i = dist_constr(jn)%atm_i
    j = dist_constr(jn)%atm_j
    
!    write(*,*)"  "
!    write(*,*)"  get_dist_constr_energy() pt 2   i =", i
!    write(*,*)"  get_dist_constr_energy() pt 2   j =", j

! Calculation of the non-bonded contact vector:

    xij = x(1, i) - x(1, j)
    yij = x(2, i) - x(2, j)
    zij = x(3, i) - x(3, j)

    rij = sqrt(xij * xij + yij * yij + zij * zij)
    
!    write(*,*)"  get_dist_constr_energy() pt 3   rij =", rij

! Calculation of the energy and deriv:

    if (dist_constr(jn)%cnt_type .eq. 0) then
          ! HARMONIC potential:
          
          da = rij - dist_constr(jn)%acoef    
          df = dist_constr(jn)%bcoef * da
          
 !         write(*,*)"  get_dist_constr_energy() pt 4   acoef =", dist_constr(jn)%acoef
 !         write(*,*)"  get_dist_constr_energy() pt 4   bcoef =", dist_constr(jn)%bcoef
          
          ene_constr = df*da 
          de = -(df + df) / rij
    else if (dist_constr(jn)%cnt_type .eq. 1) then
          ! 6-12 potential:
          
          rinv = 1.d0/rij;
          r2inv = rinv * rinv;

          r6inv = r2inv * r2inv * r2inv
          f6  = dist_constr(jn)%bcoef * r6inv
          f12 = dist_constr(jn)%acoef * (r6inv * r6inv)
                     
          ene_constr = (f12 - f6)
          de = (12.d0 * f12 - 6.d0 * f6) * r2inv
               
    else if( dist_constr(jn)%cnt_type .eq. 2) then
          ! 10-12 potential:
          rinv = 1.d0/rij;
          r2inv = rinv * rinv;
          
          r10inv = r2inv * r2inv * r2inv * r2inv * r2inv
          f10 = dist_constr(jn)%bcoef * r10inv
          f12 = dist_constr(jn)%acoef * r10inv * r2inv

          ene_constr = (f12 - f10)
          de = (12.d0 * f12 - 10.d0 * f10) * r2inv

	else if (dist_constr(jn)%cnt_type .eq. 3) then
          ! 6-12 potential + e_0 , when r < r_0 and = 0 when r > r_0
           
          rinv = 1.d0/rij;
          r2inv = rinv * rinv;

          r6inv = r2inv * r2inv * r2inv
		  r06inv = (dist_constr(jn)%bcoef/2.0d0*dist_constr(jn)%acoef); 
		  if( (r06inv/r6inv) .lt. 1.0d0) then

		   	 f6  = dist_constr(jn)%bcoef * r6inv
			 f12 = dist_constr(jn)%acoef * (r6inv * r6inv)
						
			 ene_constr = (f12 - f6) 
			 ene_constr = ene_constr - (dist_constr(jn)%acoef*r06inv*r06inv - dist_constr(jn)%bcoef*r06inv )
			 de = (12.d0 * f12 - 6.d0 * f6) * r2inv 
			
		  end if
	else
		  de = 0.0d0
    end if  ! (dist_constr(jn)%cnt_type .eq. 0)

    ! Forces:
         
      xa = de * xij
      ya = de * yij
      za = de * zij

    ! We use atm_i to determine who sums up the energy...
    if (atm_owner_map(i) .eq. mytaskid) then
      lcl_nb_constr_energy = lcl_nb_constr_energy + ene_constr
      frc(1, i) = frc(1, i) + xa
      frc(2, i) = frc(2, i) + ya
      frc(3, i) = frc(3, i) + za
    end if

    if (atm_owner_map(j) .eq. mytaskid) then
      frc(1, j) = frc(1, j) - xa
      frc(2, j) = frc(2, j) - ya
      frc(3, j) = frc(3, j) - za
    end if

  end do

  dist_constr_energy = lcl_nb_constr_energy

!  write(*,*) 'dist_constr_energy = ', dist_constr_energy

  return

end subroutine get_dist_constr_energy

end module dist_constr_mod
