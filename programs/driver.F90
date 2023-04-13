program driver

use yomhook   ,only : lhook,   dr_hook, jphook
use config_mod, only : config_type
use edir_trans_mod, only : edir_trans
use einv_trans_mod, only : einv_trans

implicit none

type(config_type) :: config
real, allocatable :: fG(:,:,:)
real, allocatable :: fM(:,:,:)
integer :: nx, ny, nfld
integer :: nproc_A, nproc_B
integer :: truncation_order
integer :: jiter, niter
integer :: ierr
real(kind=jphook) :: zhook_handle, zhook_handle_t

call mpi_init(ierr)

if (lhook) call DR_HOOK('driver',0,zhook_handle)

nx=128
ny=128
nfld=10
niter=10
nproc_A=1
nproc_B=1
truncation_order=1

call eztrans_setup(config,nx,ny,nfld,nproc_A,nproc_B,truncation_order)

allocate(fG(config%my_nx_l,config%my_ny_l,config%nfld))
allocate(fM(config%my_mx_l,config%ny+2,config%my_nfld_l))
fG=0.
fM=0.

do jiter=1,niter
  ! forward transform
  call edir_trans(config,fG,fM)

  ! backward transform
  call einv_trans(config,fM,fG)
enddo

deallocate(fG)
deallocate(fM)

call eztrans_end(config)

if (lhook) call DR_HOOK('driver',1,zhook_handle)

call mpi_finalize(ierr)

end program driver
