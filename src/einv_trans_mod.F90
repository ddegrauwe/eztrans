module einv_trans_mod

use mpi
use yomhook   ,only : lhook,   dr_hook, jphook

use config_mod, only: config_type
use trmtol_mod, only : trmtol
use trltog_mod, only : trltog

implicit none

#define USE_FFT992

contains

subroutine einv_trans(config,fM,fG)

! transform from spectral space to gridpoint space

! arguments
type(config_type), intent(in) :: config
real, intent(in)  :: fM(config%my_mx_l,config%ny+2,config%my_nfld_l)
real, intent(out) :: fG(config%my_nx_l,config%my_ny_l,config%nfld)

! local variables
real :: fL(config%nx+2,config%my_ny_l,config%my_nfld_l)

integer :: jfld, jx, jy
real(kind=jphook) :: zhook_handle,zhook_handle_t

if (lhook) call DR_HOOK('einv_trans',0,zhook_handle)

! y-fourier transform
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(jfld,jx,zhook_handle_t)
do jfld=1,config%my_nfld_l
  do jx=1,config%my_mx_l
    if (lhook) call DR_HOOK('einv_trans:fft_y',0,zhook_handle_t)
#ifdef USE_FFTW
    call dfftw_execute_dft_c2r(config%plan_bwd_single_y, fM(jx,:,jfld), fM(jx,:,jfld))
#endif
#ifdef USE_FFT992
	! FFT992(A,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
    call fft992(fM(jx,1,jfld),config%trigy,config%facy,1,1,config%ny,1,1)   ! jump is wrong for batched setup
#endif
    if (lhook) call DR_HOOK('einv_trans:fft_y',1,zhook_handle_t)
  enddo
enddo
!$OMP END PARALLEL DO

! transpose to x-pencils
call trmtol(config,fM,fL)

! x-fourier transform
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(jfld,jy,zhook_handle_t)
do jfld=1,config%my_nfld_l
  do jy=1,config%my_ny_l
    if (lhook) call DR_HOOK('einv_trans:fft_x',0,zhook_handle_t)
#ifdef USE_FFTW
    call dfftw_execute_dft_c2r(config%plan_bwd_single_x, fL(:,jy,jfld), fL(:,jy,jfld))
#endif
#ifdef USE_FFT992
	! FFT992(A,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
    call fft992(fL(1,jy,jfld),config%trigx,config%facx,1,1,config%nx,1,1)   ! jump is wrong for batched setup
#endif
    if (lhook) call DR_HOOK('einv_trans:fft_x',1,zhook_handle_t)
  enddo
enddo
!$OMP END PARALLEL DO

! transpose to z-pencils
call trltog(config,fL,fG)

if (lhook) call DR_HOOK('einv_trans',1,zhook_handle)

end subroutine einv_trans

end module einv_trans_mod
