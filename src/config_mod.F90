module config_mod

use mpi

implicit none

type config_type
  ! number of MPI tasks
  integer :: nproc
  integer :: nproc_A
  integer :: nproc_B
  
  ! global dimensions
  integer :: nx
  integer :: ny
  integer :: nfld
  integer :: mx
  integer :: my
  
  ! local values
  integer :: my_nx_l
  integer :: my_jxi  ! G 
  integer :: my_jxe  ! G 
  integer :: my_ny_l
  integer :: my_jyi  ! G, L
  integer :: my_jye  ! G, L
  integer :: my_nfld_l
  integer :: my_jfldi  ! L, M
  integer :: my_jflde  ! L, M
  integer :: my_mx_l
  integer :: my_kxi  ! M
  integer :: my_kxe  ! M
  integer :: my_nm_l  ! number of waves in M-space
  integer :: my_ns_l  ! number of waves in S-space
  integer :: my_jsi  ! S
  integer :: my_jse  ! S
  integer :: my_proc_A
  integer :: my_proc_B
  integer :: my_proc

  ! dimensions for different MPI tasks
  integer, pointer :: nx_l(:)
  integer, pointer :: ny_l(:)
  integer, pointer :: jxi_l(:)
  integer, pointer :: jxe_l(:)
  integer, pointer :: jyi_l(:)
  integer, pointer :: jye_l(:)
  integer, pointer :: jfldi_l(:)
  integer, pointer :: jflde_l(:)
  integer, pointer :: kxi_l(:)
  integer, pointer :: kxe_l(:)
  integer, pointer :: jsi_l(:)
  integer, pointer :: jse_l(:)
  integer, pointer :: nfld_l(:)
  integer, pointer :: mx_l(:)  ! number of x-wavenumbers on this proc
  
  ! M-space description
  integer, pointer :: ex(:)         ! elliptic truncation: maximum x-wavenumber for each y-wavenumber
  integer, pointer :: ey(:)         ! elliptic truncation: maximum y-wavenumber for each x-wavenumber
  integer, pointer :: nm_l(:)       ! number of waves (after truncation) for each proc in this B-set
  integer, pointer :: kx_m(:)       ! x-wavenumbers in M-space
  integer, pointer :: ky_m(:)       ! y-wavenumbers in M-space
  integer, pointer :: jproc_m(:)    ! processor where to send data from M to S

  ! S-space description
  integer, pointer :: ns_l(:)    ! number of waves for each proc in this A-set
  !integer, pointer :: my_kx(:)   ! x-wavenumbers
  !integer, pointer :: my_ky(:)   ! y-wavenumbers
  
  ! fftw plans
  integer*8 :: plan_fwd_single_x
  integer*8 :: plan_bwd_single_x
  integer*8 :: plan_fwd_single_y
  integer*8 :: plan_bwd_single_y
  
  ! fft992 data
  real, pointer :: trigx(:)
  integer, pointer :: facx(:)
  real, pointer :: trigy(:)
  integer, pointer :: facy(:)
  
  ! MPI groups
  integer :: mpi_comm_A
  integer :: mpi_comm_B
  
end type config_type

end module config_mod
