!Author: Junyi Cheng Nov.27 2018
!The 3d output of density, electrostatic field with XGC mesh
!The GEM mesh will be done later
!There are two condition:
!1: coupling 2: the GEM alone
!The same diagnosis point as XGC for the comparison,
!Should be on the point of the 2rd step of R-K method
subroutine XGC_3d_output(nstep,ndiag)
  use mpi
#ifdef COUPLING
  use gem_com, only: MyId,imx,jmx,kmx,den1
#else
  use gem_com, only: MyId,ntube,imx,jmx,kmx,den1,q,phi
  use gem_equil, only: nu,Tu
#endif 
  use coupling_core_edge, only: cce_first_node,cce_last_node,cce_all_surface_node_number,cce_first_all_node,cce_last_all_node,cce_all_node_number,nphi
  use mapping
#ifdef ADIOS  
  use adios_read_mod
  use adios_write_mod
#endif
#ifdef ADIOS2
#include "adios_macro.h"
  use adios2_comm_module
#endif

#ifdef ADIOS2
  type(adios2_engine), save :: engine_gem
  type(adios2_io), save :: io_gem
  type(adios2_variable) ::var
  integer(8), dimension(3) :: gdims_3d, goffset_3d, ldims_3d
  integer(8), dimension(2) :: gdims_2d, goffset_2d, ldims_2d
  integer, save :: init_gem=0
#endif

  integer, intent(in) :: nstep,ndiag

  integer :: stepid,err
  integer*8 :: buf_id,buf_size,total_size
  real :: e
  character(5) :: stepid_str
  character(512) :: cce_filename

  real :: den_tmp_3d(cce_all_node_number,nphi)
  real :: field_tmp_3d(cce_all_node_number,nphi)
#ifndef COUPLING
  real :: gem_3d(0:imx,0:jmx,0:kmx)
  real :: pot_3d(cce_all_surface_node_number,nphi)
  real :: gem_tmp(0:imx,0:jmx,0:1)  
#endif
  
  !The diagnosis point
  if(modulo(nstep,ndiag) .eq. 0)then
    stepid=int(nstep/ndiag)+1
    write(stepid_str,'(I0.5)')stepid
  else
    return
  endif
 
  den_tmp_3d=0.0
  field_tmp_3d=0.0
  
#ifndef COUPLING
  gem_3d=0.0
  pot_3d=0.0
  gem_tmp=0.0
#endif

#ifdef COUPLING
  !change the 3d output with the all points
  if(MyId==0)then 
    den_tmp_3d(cce_first_node:cce_last_node,:)=pot_density_3d(1:cce_all_surface_node_number,:)
    field_tmp_3d=pot_field_3d
  endif    
#else
  !construt the 3d ouput with all points:
  !1: mapping 2: change to the all points format
  
  !density part
  gem_tmp=den1(1,:,:,:)/real(q(1))*nu
#ifndef MAP_PARALLEL
  call mapping_GEM_XGC(gem_3d,pot_3d,gem_tmp,grid,mapping_GEM_XGC_linear_coef)
#else
  call mapping_GEM_XGC_new(gem_3d,pot_3d,gem_tmp,grid,mapping_GEM_XGC_linear_coef)
#endif
  if(MyId==0)then 
    den_tmp_3d(cce_first_node:cce_last_node,:)=pot_3d(1:cce_all_surface_node_number,:)
    density_3d=gem_3d
  endif

  !field part
  e=1.6022e-19
  gem_tmp=phi*Tu/e/real(ntube)
#ifndef MAP_PARALLEL
  call mapping_GEM_XGC(gem_3d,pot_3d,gem_tmp,grid,mapping_GEM_XGC_linear_coef)
#else
  call mapping_GEM_XGC_new(gem_3d,pot_3d,gem_tmp,grid,mapping_GEM_XGC_linear_coef)
#endif
  if(MyId==0)then
    field_tmp_3d(cce_first_node:cce_last_node,:)=pot_3d(1:cce_all_surface_node_number,:)
    field_3d=gem_3d
  endif
#endif

  !output the 3d data
  if(MyId==0)then
#ifdef ADIOS
    cce_filename='gem_xgc.3d.'//trim(stepid_str)//'.bp'
    call ADIOS_OPEN(buf_id,'gem_xgc_3d',cce_filename,'w',MPI_COMM_SELF,err)
    buf_size=4*5+cce_all_node_number*nphi*16+16*(imx+1)*(jmx+1)*(kmx+1)+100 !100 is buff
    
    call ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
    call ADIOS_WRITE(buf_id,"nphi",nphi,err)
    call ADIOS_WRITE(buf_id,"cce_all_node_number",cce_all_node_number,err)
    call ADIOS_WRITE(buf_id,"imx_1",imx+1,err)
    call ADIOS_WRITE(buf_id,"jmx_1",jmx+1,err)
    call ADIOS_WRITE(buf_id,"kmx_1",kmx+1,err)
    !actual data
    call ADIOS_WRITE(buf_id,"xgc_den_3d",den_tmp_3d,err)
    call ADIOS_WRITE(buf_id,"xgc_field_3d",field_tmp_3d,err)
    call ADIOS_WRITE(buf_id,"gem_den_3d",density_3d,err)
    call ADIOS_WRITE(buf_id,"gem_field_3d",field_3d,err)
    
    call ADIOS_CLOSE(buf_id,err) 
#elif ADIOS2
    gdims_3d(1)=imx+1
    gdims_3d(2)=jmx+1
    gdims_3d(3)=kmx+1
    goffset_3d(1)=0
    goffset_3d(2)=0
    goffset_3d(3)=0
    ldims_3d(1)=imx+1
    ldims_3d(2)=jmx+1
    ldims_3d(3)=kmx+1

    gdims_2d(1)=cce_all_node_number
    gdims_2d(2)=nphi
    goffset_2d(1)=0
    goffset_2d(2)=0
    ldims_2d(1)=cce_all_node_number
    ldims_2d(2)=nphi

    cce_filename='gem_xgc.3d.'//trim(stepid_str)//'.bp'
    if(init_gem==0)then
       call adios2_declare_io(io_gem, adios2obj, 'output_3d', err)
       call adios2_define_variable(var, io_gem, "gem_den_3d", adios2_type_dp, 3, gdims_3d, goffset_3d, ldims_3d, adios2_constant_dims, err)
       call adios2_define_variable(var, io_gem, "gem_field_3d", adios2_type_dp, 3, gdims_3d, goffset_3d, ldims_3d, adios2_constant_dims, err)
       call adios2_define_variable(var, io_gem, "xgc_den_3d", adios2_type_dp, 2, gdims_2d, goffset_2d, ldims_2d, adios2_constant_dims, err)
       call adios2_define_variable(var, io_gem, "xgc_field_3d", adios2_type_dp, 2, gdims_2d, goffset_2d, ldims_2d, adios2_constant_dims, err)
       init_gem=1
    endif

    call adios2_open(engine_gem, io_gem, cce_filename, adios2_mode_write, mpi_comm_self, err)
    call adios2_begin_step(engine_gem, adios2_step_mode_append, err)
    call adios2_put(engine_gem, "gem_den_3d", density_3d, err)
    call adios2_put(engine_gem, "gem_field_3d", field_3d, err)
    call adios2_put(engine_gem, "xgc_den_3d", den_tmp_3d, err)
    call adios2_put(engine_gem, "xgc_field_3d", field_tmp_3d, err)
    call adios2_end_step(engine_gem, err)
    call adios2_close(engine_gem, err)
     
#endif
  endif 
end subroutine XGC_3d_output
