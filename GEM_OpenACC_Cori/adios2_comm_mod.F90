
module adios2_comm_module
    ! use gem_com,  only:myid
    use gem_com,  only: myid
    use mpi
#ifdef ADIOS2
    use adios2
#endif
    implicit none
#ifdef ADIOS2
    type(adios2_adios) :: adios2obj
    type(adios2_engine), allocatable :: list_engines(:)
    integer :: n_engines

    !! very simple single timer
    character(len=128) :: timer_name
    integer :: timer_comm
    integer :: timer_index
    real(kind=8) :: t_start

contains
    subroutine adios2_comm_init(initfile)
        !use input_module
        implicit none
        character(len=*), intent(in) :: initfile
        integer :: ierr
		
		!@effis-init xml=initfile, comm=sml_comm
		call adios2_init(adios2obj, initfile, mpi_comm_world, .true., ierr)
        allocate(list_engines(16))
        n_engines = 0
    end subroutine adios2_comm_init

    subroutine adios2_comm_finalize()
        implicit none
        integer :: ierr
        integer :: i

        do i = 1, n_engines
            if (myid.eq.0) print *, 'ADIOS2: close output ', trim(list_engines(i)%name)
            call adios2_close(list_engines(i), ierr)
        enddo
		!@effis-finalize
        call adios2_finalize(adios2obj, ierr)
    end subroutine adios2_comm_finalize

#endif
end module adios2_comm_module
