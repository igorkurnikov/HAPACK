subroutine open_mdcrd_form
    use file_io_dat_mod
    use file_io_mod
    use prmtop_dat_mod
    
    call amopen(mdcrd, mdcrd_name, 'U', 'F', 'W')
    write(mdcrd, 1000) prmtop_ititl
1000 format(a80)
end subroutine open_mdcrd_form

subroutine open_mdvel_form
    use file_io_dat_mod
    use file_io_mod
    use prmtop_dat_mod
    
    call amopen(mdvel, mdvel_name, 'U', 'F', 'W')
    write(mdvel, 1000) prmtop_ititl
    
1000 format(a80)
end subroutine open_mdvel_form

subroutine open_mdcrd_bin4( ntb_par, ntwprt_par )
    use file_io_dat_mod
    use file_io_mod
    use prmtop_dat_mod
    use mdin_ctrl_dat_mod
    
    implicit none
    
    integer                       :: ntb_par
    integer                       :: ntwprt_par
    
    character(10), parameter      :: file_version = '9.00'
    integer                       :: box_flag
    character(100)  ::  bin4_title
    
    bin4_title = trim(mdcrd_name) // '.bin4'
    call amopen(mdcrd, bin4_title, 'U', 'U', 'W') 
    write(mdcrd) file_version
    write(mdcrd) prmtop_ititl

    if (ntb_par .gt. 0) then
        box_flag = 1
    else
        box_flag = 0
    end if

    if ( ntwprt_par .ne. 0) then
        write(mdcrd) ntwprt_par, box_flag
    else
        write(mdcrd) natom, box_flag
    end if
    
end subroutine open_mdcrd_bin4

subroutine open_mdvel_bin4( ntwprt_par )
    use file_io_dat_mod
    use file_io_mod
    use prmtop_dat_mod
    use mdin_ctrl_dat_mod
    
    implicit none
    
    integer                       :: ntwprt_par
    
    character(10), parameter      :: file_version = '9.00'
    integer                       :: box_flag
    character(100)  ::  bin4_title
    
    bin4_title = trim(mdvel_name) // '.bin4'
    call amopen(mdvel, bin4_title, 'U', 'U', 'W')   
    write(mdvel) file_version
    write(mdvel) prmtop_ititl

    box_flag = 0

    if (ntwprt_par .ne. 0) then
        write(mdvel) ntwprt_par, box_flag
    else
        write(mdvel) natom, box_flag
    end if
    
end subroutine open_mdvel_bin4

subroutine open_mden
    use file_io_dat_mod
    use file_io_mod
    
    call amopen(mden, mden_name, 'U', 'F', 'W')
    
end subroutine open_mden

subroutine open_restart_form
    use file_io_dat_mod
    use file_io_mod
    
    call amopen(restrt, restrt_name, 'U', 'F', 'W') 
    
end subroutine open_restart_form

subroutine open_restart_bin
    use file_io_dat_mod
    use file_io_mod
    
    call amopen(restrt, restrt_name, 'U', 'U', 'W') 
    
end subroutine open_restart_bin

subroutine close_mdcrd
    use file_io_dat_mod
    close(mdcrd)
end subroutine close_mdcrd

subroutine close_mdvel
    use file_io_dat_mod
    
    close(mdvel)
end subroutine close_mdvel

subroutine close_mden
    use file_io_dat_mod
    
    close(mden)
end subroutine close_mden

subroutine close_restart
    use file_io_dat_mod
    
    close(restrt)
end subroutine close_restart


