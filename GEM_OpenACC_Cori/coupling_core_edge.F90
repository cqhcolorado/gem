!author Junyi Cheng Oct 22 2018
!coupling
module coupling_core_edge
  use mpi
#ifdef ADIOS
  use adios_read_mod
#endif
#ifdef ADIOS2
  use adios2_comm_module
#endif
  implicit none

  integer :: cce_side,nphi,nwedge,coef_opt

  integer :: cce_first,cce_last,cce_all,cce_first_all_node,cce_last_all_node,cce_all_node_number

  integer :: cce_first_surface,cce_last_surface,cce_all_surface_number

  integer :: cce_first_node,cce_last_node,cce_all_surface_node_number

  integer :: cce_first_field,cce_last_field,cce_all_field_number

  integer :: cce_first_field_node,cce_last_field_node,cce_all_field_node_number

  integer :: cce_first_surface_coupling,cce_last_surface_coupling, cce_core_all
 
  real :: cce_in_weight,cce_out_weight,rhomax

  integer, dimension(:), allocatable :: cce_surface_first_node,cce_surface_last_node

  integer :: cce_step,cce_field_step

  real, dimension(:), allocatable :: cce_density,cce_field

  character(256) :: cce_folder,eq_filename,node_filename,surf_filename
  
  character(5) :: cce_my_side,cce_other_side

#ifdef ADIOS
  character (len=10) :: staging_read_method_name

  integer :: staging_read_method = ADIOS_READ_METHOD_BP !! internal value
#endif
contains

subroutine cce_initialize
  use mpi
#ifdef ADIOS
  use adios_read_mod
  use adios_write_mod
#endif
  use gem_com, only:myid

  integer :: i,ierr

  namelist /coupling/ cce_side,cce_folder,eq_filename,node_filename,surf_filename,   &
                      cce_first,cce_last,cce_first_surface,cce_last_surface,  &
                      cce_first_field,cce_last_field,      &
                      cce_first_surface_coupling,cce_last_surface_coupling, cce_core_all,&
                      cce_step,cce_field_step,nphi,nwedge,coef_opt,rhomax

  namelist /surfaces/ cce_surface_first_node,cce_surface_last_node

  open(unit=20,file='coupling.in', status='old',action='read')
  READ(20, NML=coupling)
  close(unit=20)

  if(cce_side==0)then
    cce_my_side='core'
    cce_other_side='edge'
  else
    cce_my_side='edge'
    cce_other_side='core'
  endif

  cce_all=cce_last-cce_first+1
 
  cce_all_surface_number=cce_last_surface-cce_first_surface+1

  cce_all_field_number=cce_last_field-cce_first_field+1

  allocate(cce_surface_first_node(cce_all))
  allocate(cce_surface_last_node(cce_all))

call mpi_barrier(mpi_comm_world,ierr)
if(myid==0)write(*,*)'after allocate cce 1'

  open(unit=20,file='coupling.in',status='old',action='read')
  read(20,NML=surfaces)
  close(unit=20)

  cce_first_node=cce_surface_first_node(cce_first_surface)
  cce_last_node=cce_surface_last_node(cce_last_surface)
  cce_all_surface_node_number=cce_last_node-cce_first_node+1

  cce_first_field_node=cce_surface_first_node(cce_first_field)
  cce_last_field_node=cce_surface_last_node(cce_last_field)
  cce_all_field_node_number=cce_last_field_node-cce_first_field_node+1

  cce_first_all_node=cce_surface_first_node(cce_first)
  cce_last_all_node=cce_surface_last_node(cce_last)
  cce_all_node_number=cce_last_all_node-cce_first_all_node+1

  !write(*,*)'field',cce_first_field_node,cce_last_field_node,cce_all_field_node_number
  !write(*,*)'all',cce_first_all_node,cce_last_all_node,cce_all_node_number
  allocate(cce_density(cce_all_surface_node_number))
  allocate(cce_field(cce_all_node_number))

call mpi_barrier(mpi_comm_world,ierr)
if(myid==0)write(*,*)'after allocate cce 2',cce_core_all,'cce_core_all',cce_first_surface_coupling,'cce_first_surface_coupling',cce_last_surface_coupling,'cce_last_surface_coupling'

  cce_in_weight=real(cce_first_surface_coupling-1)/real(cce_core_all-1)
  cce_out_weight=real(cce_last_surface_coupling-1)/real(cce_core_all-1)
call mpi_barrier(mpi_comm_world,ierr)
#ifdef ADIOS
if(myid==0)write(*,*)'after cce 2.0'
  staging_read_method_name = 'BP'//char(0)
  staging_read_method = adios_comm_get_read_method(trim(staging_read_method_name))
  if(myid==0) print *, 'staging_read_method_name=',trim(staging_read_method_name)
  if(myid==0) print *, 'staging_read_method=', staging_read_method
  call adios_init('adioscfg.xml'//char(0), mpi_comm_world, ierr)
  call adios_read_init_method(staging_read_method, mpi_comm_world, 'verbose=3', ierr)
if(myid==0)write(*,*)'adios init', staging_read_method_name
call mpi_barrier(mpi_comm_world,ierr)
if(myid==0)write(*,*)'after cce 3'
#endif

#ifdef ADIOS2
  call adios2_comm_init('adios2cfg.xml')
#endif 

end subroutine cce_initialize

subroutine cce_destroy
  use mpi
#ifdef ADIOS
  use adios_read_mod
  use adios_write_mod
#endif
  integer :: ierr

  if(allocated(cce_density))deallocate(cce_density)
  if(allocated(cce_field))deallocate(cce_field)
  if(allocated(cce_surface_first_node))deallocate(cce_surface_first_node)
  if(allocated(cce_surface_last_node))deallocate(cce_surface_last_node)
#ifdef ADIOS
  call adios_finalize(MPI_COMM_WORLD,ierr)
  call adios_read_finalize_method(staging_read_method,ierr)
#endif
end subroutine cce_destroy

subroutine cce_send_density(pot_density_3d, electron_on)
  use mpi
#ifdef ADIOS
  use adios_write_mod
#endif
  use gem_com, only:jmx,MyId,rho,gemout
  use gem_equil, only:nu,cn0i

  real, intent(in) :: pot_density_3d(cce_all_surface_node_number,nphi)
  integer, intent(in) :: electron_on

  character(5) :: cce_stepstr,planestr
  character(512) :: cce_filename

  integer*8 :: buf_id,buf_size,total_size
  integer :: i,j,k,n,ii,err,ierr
  
  real*8, dimension(:), allocatable :: arrtmp
 
#ifdef ADIOS2
  integer(8), dimension(2) :: gdims, goffset, ldims
  type(adios2_io), save :: io, io_electron
  type(adios2_engine), save :: engine, engine_electron
  logical, save :: init=.false.
  type(adios2_variable) :: varid
#endif

  cce_density=0.0
 
  !write(planestr,'(I0.5)')MyId
  write(cce_stepstr,'(I0.5)')cce_step
  if(cce_side .eq. 0)then
 
if(myid==0)then
    write(gemout,*)'cce_all_surface_node_number,nphi',cce_all_surface_node_number,nphi
    call flush(gemout)
  endif
 
#ifdef ADIOS2
    if(MyId==0)then
      if(.not. init)then
        gdims(1)=cce_all_surface_node_number
        gdims(2)=nphi
        goffset(1)=0
        goffset(2)=0
        ldims(1)=cce_all_surface_node_number
        ldims(2)=nphi

        call adios2_declare_io(io,adios2obj,'density',ierr)
        call adios2_define_variable(varid, io, "data",  adios2_type_dp, 2, gdims, goffset, ldims, adios2_constant_dims, ierr)
        call adios2_open(engine, io, trim(cce_folder) // '/' //'density.bp', adios2_mode_write, mpi_comm_self, ierr)
        call adios2_declare_io(io_electron,adios2obj,'electron_density',ierr)
        call adios2_define_variable(varid, io_electron, "data",  adios2_type_dp, 2, gdims, goffset, ldims, adios2_constant_dims, ierr)
        call adios2_open(engine_electron, io_electron, trim(cce_folder) // '/' //'electron_density.bp',adios2_mode_write, mpi_comm_self, ierr)
        init=.true.
      endif
    
      if(electron_on==0)then 
        write(*,*)'max density',maxval(pot_density_3d*nu)
        call adios2_begin_step(engine, adios2_step_mode_append, ierr)
        call adios2_put(engine, "data",  pot_density_3d*nu, ierr)
        call adios2_end_step(engine,ierr)
      else
        call adios2_begin_step(engine_electron, adios2_step_mode_append, ierr)
        call adios2_put(engine_electron, "data",  pot_density_3d*nu, ierr)
        call adios2_end_step(engine_electron,ierr)
      endif 
    endif
#elif ADIOS  
    if(MyId==0)then
      allocate(arrtmp(cce_all_surface_node_number))

      do i=1,nphi
        arrtmp=pot_density_3d(:,i)*nu
        ii=i-1
        write(planestr,'(I0.5)')ii
        !cce_filename=trim(cce_folder)//'/density_'//trim(cce_my_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.bp'
        cce_filename=trim(cce_folder)//'/density_'//trim(cce_my_side)//'_'//trim(planestr)//'.bp'
        call ADIOS_OPEN(buf_id,'coupling',cce_filename,'w',MPI_COMM_SELF,err)
        buf_size= 4*6 + 8*cce_all_surface_node_number  + 100 !last 100 is buffer

        call ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
        call ADIOS_WRITE(buf_id,"nphi",nphi,err)
        call ADIOS_WRITE(buf_id,"iphi",ii,err)
        call ADIOS_WRITE(buf_id,"first_node",cce_first_node,err)
        call ADIOS_WRITE(buf_id,"last_node",cce_last_node,err)
        call ADIOS_WRITE(buf_id,"node_number",cce_all_surface_node_number,err)
        call ADIOS_WRITE(buf_id,"cce_side",cce_side,err)
        call ADIOS_WRITE(buf_id,"cce_model",1,err)
        call ADIOS_WRITE(buf_id,"time",1.0,err)
        call ADIOS_WRITE(buf_id,"step",cce_step,err)
        !actual data
        call ADIOS_WRITE(buf_id,"data", arrtmp,err)

        call ADIOS_CLOSE(buf_id,err)

        !Create an unlock file
        cce_filename=trim(cce_folder)//'/density_'//trim(cce_my_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.unlock'
        open(20, file=cce_filename, status="new", action="write")
        close(20)
      enddo
      
      deallocate(arrtmp)
    endif
#endif   
  endif
  
end subroutine cce_send_density

subroutine cce_receive_density
#ifdef ADIOS
  use adios_read_mod
#endif
  use mpi
  use gem_com,only:jmx,MyId

  logical :: ex
  character(5)::cce_stepstr,planestr
  character(512) :: cce_filename
  integer*8 :: buf_id
  integer :: err
  integer*8 :: sel1=0
  integer :: stat
  real*8, dimension(:),allocatable :: arrtmp

  cce_density=0.0

  if(cce_side .eq. 1)then
    !write(cce_stepstr,'(I0.5)') cc_step
    !write(planestr,'(I0.5)') MyId

    !cce_filename=trim(cce_folder)//'/density_'//trim(cce_other_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.unlock'
    !err=-1
    !do while(err/=0)
    !  do while(.NOT.ex)
    !    inquire(file=cce_filename,EXIST=ex)
    !  end do
    !  open(20,file=cce_filename,iostat=stat,status="old")
    !  if(stat==0)close(20,status='delete')
    !  !cce_filename=trim(cce_folder)//'/density_'//trim(cce_other_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.bp'
    !  cce_filename=trim(cce_folder)//'/density_'//trim(cce_other_side)//'_'//trim(planestr)//'.bp'
    !  call adios_read_open_file (buf_id, cce_filename, 0, MPI_COMM_SELF, err)
    !  if(err/=0) then
    !    write(*,*)'coupling receive error: could not open file', cce_filename
    !  endif
    !enddo
    !allocate(arrtmp(cce_all_surface_node_number))
    !arrtmp=0.0
    !call adios_schedule_read (buf_id, sel1, 'data', 0, 1, arrtmp, err)
    !call adios_perform_reads (buf_id, err)
    !call adios_read_close (buf_id, err)

    !cce_density(1:cce_all_surface_node_number)=arrtmp(1:cce_all_surface_node_number)
    !if(MyId==0)write(*,*)'receive density', cce_density((1+cce_all_surface_number*jmx/2):(cce_all_surface_number*(jmx/2+1)))
    !deallocate(arrtmp)
  endif
  
end subroutine cce_receive_density

subroutine cce_process_density
  use mpi
  use gem_com,only:jmx,rho,MyId

  integer :: i,j,k,n
  real :: temp

  select case (cce_side)
    case (0)
    case (1)
      !if(MyId==0)then
      ! write(*,*)'before rho',rho(cce_first_field:cce_last_field,jmx/2,0)
      !endif
      
      !do i=cce_first_surface,cce_last_surface
      !   do j=0,jmx
      !      !do k=0,1
      !         n=i+1+cce_all_surface_number*j
      !         temp=cce_density(n)
      !         rho(i,j,0)=rho(i,j,0)+temp
      !         temp=cce_density(cce_all_surface_number*(jmx+1)+n)
      !         rho(i,j,1)=rho(i,j,1)+temp
            !enddo
      !   enddo
      !enddo
      
      !if(MyId==0)write(*,*)'after rho',rho(cce_first_field:cce_last_field,jmx/2,0)
   case default
     write(*,*)'unknown coupling model'
  end select

  cce_step=cce_step+1

end subroutine cce_process_density

subroutine cce_send_field
  use mpi
#ifdef ADIOS
  use adios_write_mod
#endif
  use gem_com,only:jmx,MyId,phi

  character(5) :: cce_stepstr,planestr
  character(512) :: cce_filename

  integer*8 :: buf_id,buf_size,total_size
  integer :: i,j,k,n,err

  real*8, dimension(:), allocatable :: arrtmp

  cce_field=0.0

  write(planestr,'(I0.5)')MyId
  write(cce_stepstr,'(I0.5)')cce_field_step

  if(cce_side .eq. 1)then
    !allocate(arrtmp(cce_all_field_node_number))
    !do i=cce_first_field,cce_last_field
    !   do j=0,jmx
    !      !do k=0,1
    !         n=i+1+cce_all_field_number*j
    !         arrtmp(n)=phi(i,j,0)
    !         arrtmp(cce_all_field_number*(jmx+1)+n)=phi(i,j,1)
    !      !enddo
    !   enddo
    !enddo

    !!cce_filename=trim(cce_folder)//'/field_'//trim(cce_my_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.bp'
    !cce_filename=trim(cce_folder)//'/field_'//trim(cce_my_side)//'_'//trim(planestr)//'.bp'

    !call ADIOS_OPEN(buf_id,'coupling_fields',cce_filename,'w',MPI_COMM_SELF,err)  !sml_comm_null,err)
    !buf_size= 4*6 + 8*cce_all_field_node_number  + 100 !last 100 is buffer

    !call ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
    !call ADIOS_WRITE(buf_id,'MyId',MyId,err)
    !call ADIOS_WRITE(buf_id,'first_field',cce_first_field,err)
    !call ADIOS_WRITE(buf_id,'last_field',cce_last_field,err)
    !call ADIOS_WRITE(buf_id,'cce_side',cce_side,err)
    !call ADIOS_WRITE(buf_id,'cce_field_step',cce_step,err)
    !call ADIOS_WRITE(buf_id,'cce_all_field_node_number',cce_all_field_node_number,err)
    !actual data
    !call ADIOS_WRITE(buf_id,'data', arrtmp,err)

    !call ADIOS_CLOSE(buf_id,err)
    !deallocate(arrtmp)

    !Create an unlock file
    !cce_filename=trim(cce_folder)//'/field_'//trim(cce_my_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.unlock'
    !open(20, file=cce_filename, status="new", action="write")
    !close(20)

  endif
end subroutine cce_send_field

subroutine cce_receive_field(pot_field_3d,opt_adi)
#ifdef ADIOS
  use adios_read_mod
#endif
  use mpi
  use gem_com, only:MyId,phi,gemout
  use gem_equil, only:Tu

  real, intent(out) :: pot_field_3d(cce_all_node_number,nphi)

#ifdef ADIOS2
  integer(8), dimension(2) :: starts,counts
  type(adios2_io), save :: io_field, io_field_electron
  type(adios2_engine), save :: engine, engine_electron
  type(adios2_variable) :: varid 
  logical, save :: init=.false.  
#endif
  logical :: ex
  character(5)::cce_stepstr,planestr
  character(512) :: cce_filename
  integer*8 :: buf_id
  integer ::  i,ii,err,stat,ierr
  integer*8 :: sel1=0
  real :: e,phi_env
  real*8, dimension(:),allocatable :: arrtmp
  integer :: opt_adi

  cce_field=0.0

  e=1.6022e-19
  phi_env=e/Tu


  if(myid==0)then
    write(gemout,*)'cce_all_field_node_number,nphi',cce_all_field_node_number,nphi
    call flush(gemout)
  endif
  if(cce_side .eq. 0)then
#ifdef ADIOS2
  if(MyId==0)then
    starts(1)=0
    counts(1)=cce_all_field_node_number
    starts(2)=0
    counts(2)=nphi
    if(.not. init)then
      call adios2_declare_io(io_field,adios2obj,'field',ierr)    
      call adios2_open(engine, io_field, trim(cce_folder) // '/' // 'field.bp', adios2_mode_read, mpi_comm_self, ierr)
      call adios2_declare_io(io_field_electron, adios2obj, 'electron_field', ierr)
      call adios2_open(engine_electron, io_field_electron, trim(cce_folder) // '/' // 'electron_field.bp', adios2_mode_read, mpi_comm_self, ierr)
      init=.true.
    endif
    if(opt_adi==0)then
      call adios2_begin_step(engine, adios2_step_mode_read, ierr)
      call adios2_inquire_variable(varid, io_field, "dpot1", ierr)
      call adios2_set_selection(varid, 2, starts, counts, ierr)
      call adios2_get(engine, varid, pot_field_3d(cce_first_field_node:cce_last_field_node,:), adios2_mode_deferred, ierr)
      call adios2_end_step(engine, ierr)  
    else
      call adios2_begin_step(engine_electron, adios2_step_mode_read, ierr)
      call adios2_inquire_variable(varid, io_field_electron, "dpot0", ierr)
      call adios2_set_selection(varid, 2, starts, counts, ierr)
      call adios2_get(engine_electron, varid, pot_field_3d(cce_first_field_node:cce_last_field_node,:), adios2_mode_deferred, ierr)
      call adios2_end_step(engine_electron, ierr)
    endif
    pot_field_3d=pot_field_3d*phi_env
  endif
#elif ADIOS
    write(cce_stepstr,'(I0.5)') cce_field_step

    if(MyId==0)then
      allocate(arrtmp(cce_all_field_node_number))

      do i=1,nphi
        ii=i-1
        write(planestr,'(I0.5)') ii

        cce_filename=trim(cce_folder)//'/field_'//trim(cce_other_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.unlock'
        err=-1
        !write(*,*)'before inquire',cce_filename,i
        do while(err/=0)
          do while(.NOT.ex)
            inquire(file=cce_filename,EXIST=ex)
          end do
          !write(*,*)'unlock open',cce_filename,i
          open(20,file=cce_filename,iostat=stat,status="old")
          if(stat==0)close(20,status='delete')
          !write(*,*)'before bp', cce_filename,i
          !cce_filename=trim(cce_folder)//'/field_'//trim(cce_other_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'.bp'
          
          cce_filename=trim(cce_folder)//'/field_'//trim(cce_other_side)//'_'//trim(planestr)//'.bp'
          !write(*,*)'after bp', cce_filename,i
          call adios_read_open_file (buf_id, cce_filename, 0, MPI_COMM_SELF, err)
          if(err/=0) then
            write(*,*)'coupling receive error: could not open file', cce_filename
          endif
        enddo
        arrtmp=0.0
        call adios_schedule_read (buf_id, sel1, 'dpot1', 0, 1, arrtmp, err)
        call adios_perform_reads (buf_id, err)
        call adios_read_close (buf_id, err)
        !write(*,*)'after bp reader', cce_filename, i
        !write(*,*)'iphi,arrtmp',i,arrtmp(63388)
        cce_field=0.0
        cce_field(cce_first_field_node:cce_last_field_node)=arrtmp(1:cce_all_field_node_number)*phi_env
        pot_field_3d(:,i)=cce_field
      enddo

      !write(*,*)'diag point pot3d',pot_field_3d(63388,1)
      deallocate(arrtmp)
    endif
#endif
  endif

end subroutine cce_receive_field

subroutine cce_process_field
  use mpi
  use gem_com,only:jmx,phi

  integer :: i,j,k,n
  real :: temp

  select case (cce_side)
    case (1)
    case (0)
      !do i=cce_first_surface,cce_last_surface
      !   do j=0,jmx
      !      !do k=0,1
      !         n=i+1+cce_all_field_number*j
      !         phi(i,j,0)=cce_field(n)
      !         phi(i,j,1)=cce_field(cce_all_field_number*(jmx+1)+n)
            !enddo
      !   enddo
      !enddo
   case default
     write(*,*)'unknown coupling model'
  end select

  cce_field_step=cce_field_step+1
  
end subroutine cce_process_field
#ifdef ADIOS
function adios_comm_get_read_method(method)
  implicit none
  integer :: adios_comm_get_read_method
  character(len=*), intent(in) :: method

  adios_comm_get_read_method = ADIOS_READ_METHOD_BP
  if (trim(method) .eq. 'FLEXPATH') then
    adios_comm_get_read_method = ADIOS_READ_METHOD_FLEXPATH
  else if (trim(method) .eq. 'DATASPACES') then
    adios_comm_get_read_method = ADIOS_READ_METHOD_DATASPACES
  else if (trim(method) .eq. 'DIMES') then
    adios_comm_get_read_method = ADIOS_READ_METHOD_DIMES
  end if
end function adios_comm_get_read_method
#endif
end module coupling_core_edge
