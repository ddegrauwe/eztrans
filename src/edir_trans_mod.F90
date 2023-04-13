module edir_trans_mod

use mpi
use yomhook   ,only : lhook,   dr_hook, jphook

use config_mod, only: config_type
use trgtol_mod, only : trgtol
use trltom_mod, only : trltom

implicit none


contains

subroutine edir_trans(config,fG,fM)

! transform from gridpoint space to spectral space

! arguments
type(config_type), intent(in) :: config
real, intent(in)  :: fG(:,:,:)     ! (config%my_nx_l,config%my_ny_l,config%nfld)
real, intent(out) :: fM(:,:,:)   ! (config%my_mx_l,config%my,config%my_nfld_l)

! local variables
real :: fL(config%nx,config%my_ny_l,config%my_nfld_l)
real :: fLs(config%nx+2,config%my_ny_l,config%my_nfld_l)  ! fL in spectral x-space
real :: fMs(config%my_mx_l,config%ny,config%my_nfld_l) ! fM in spectral x-space

integer :: jfld, jx, jy
real(kind=jphook) :: zhook_handle,zhook_handle_t

if (lhook) call DR_HOOK('edir_trans',0,zhook_handle)

! transpose to x-pencils
call trgtol(config,fG,fL)

! x-fourier transform
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(jfld,jy,zhook_handle_t)
do jfld=1,config%my_nfld_l
  do jy=1,config%my_ny_l
    if (lhook) call DR_HOOK('edir_trans:fft_x',0,zhook_handle_t)
    call dfftw_execute_dft_r2c(config%plan_fwd_single_x, fL(:,jy,jfld), fLs(:,jy,jfld))
    if (lhook) call DR_HOOK('edir_trans:fft_x',1,zhook_handle_t)
  enddo
enddo
!$OMP END PARALLEL DO

! transpose to y-pencils
call trltom(config,fLs,fMs)

! y-fourier transform
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(jfld,jx,zhook_handle_t)
do jfld=1,config%my_nfld_l
  do jx=1,config%my_mx_l
    if (lhook) call DR_HOOK('edir_trans:fft_y',0,zhook_handle_t)
    call dfftw_execute_dft_r2c(config%plan_fwd_single_y, fMs(jx,:,jfld), fM(jx,:,jfld))
    if (lhook) call DR_HOOK('edir_trans:fft_y',1,zhook_handle_t)
  enddo
enddo
!$OMP END PARALLEL DO


if (lhook) call DR_HOOK('edir_trans',1,zhook_handle)

end subroutine edir_trans

end module edir_trans_mod
