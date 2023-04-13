module tests_mod

implicit none

contains


subroutine test_ellips(mx,my)

use aux_mod

implicit none

integer, intent(in) :: mx,my
integer :: ex(my), ey(mx)
integer :: ex2(my), ey2(mx)
logical :: M(mx,my)
integer :: jx, jy
real :: kx, ky

write (*,*) '-------------------'
write (*,*)
write (*,*) 'Running test_ellips for mx = ',mx,'; my = ',my

call ellips(mx,my,ex,ey)

!write (*,*) 'elliptic truncation for mx = ',mx,', my = ',my,':'
!write (*,'(X,A,999I4)') 'ey = ',ey
!write (*,'(X,A,999I4)') 'ex = ',ex

M=.false.
do jy=1,my
  ky=real((jy-1)/2)/(my/2-1)
  do jx=1,mx
    kx=real((jx-1)/2)/(mx/2-1)
	M(jx,jy)= sqrt(kx**2+ky**2) < 1.+1.e-6
  enddo
enddo

do jx=1,mx
  do jy=1,my
    if ( M(jx,jy) ) then
	  ey2(jx)=jy
	  ex2(jy)=jx
	endif
  enddo
enddo

!write (*,'(X,A,999I4)') 'ey2 = ',ey2
!write (*,'(X,A,999I4)') 'ex2 = ',ex2

if ( maxval(abs(ex-ex2)) > 0 ) then
  write (*,*) '    difference on ex: ',maxval(abs(ex-ex2)), '    (small difference can be due to round-off errors'
endif
if ( maxval(abs(ey-ey2)) > 0 ) then
  write (*,*) '    difference on ey: ',maxval(abs(ey-ey2)), '    (small difference can be due to round-off errors'
endif

if ( sum(ex) /= sum(ey) ) then
  write (*,*) '  ERROR: inconsistency in x and y direction of elliptic truncation'
  call abort()
endif

write (*,*) '  All good'
write (*,*)

end subroutine test_ellips


subroutine test_fftw(nx,nfld)

use fftw3_mod

implicit none

integer, intent(in) :: nx, nfld

real :: fG(nx,nfld)
real :: fS(nx+2,nfld)
real, parameter :: PI=acos(-1.)
integer*8 :: plan_fwd, plan_bwd
integer :: ifld, ix
real :: errval

write (*,*) '-------------------'
write (*,*)
write (*,*) 'Running test_fftw for nx = ',nx,'; nfld = ',nfld


! fill with sensible values
do ifld=1,nfld
  do ix=1,nx
    fG(ix,ifld)=cos(2*PI*(ifld-1)*real(ix-1)/nx)
  enddo
enddo

! create plans
call dfftw_plan_dft_r2c_1d(plan_fwd,nx,fG,fS,FFTW_ESTIMATE)
call dfftw_plan_dft_c2r_1d(plan_bwd,nx,fS,fG,FFTW_ESTIMATE)

! forward transform
do ifld=1,nfld
  call dfftw_execute_dft_r2c(plan_fwd, fG(:,ifld), fS(:,ifld))
enddo

! verify result
errval=0.
do ifld=1,nfld
  do ix=1,nx
    if (ix == 2*ifld-1) then
	  if (ifld==1) then
		errval=max(errval,(fS(ix,ifld)-nx)**2)
	  else
	    errval=max(errval,(fS(ix,ifld)-nx/2)**2)
      endif
	else
      errval=max(errval,fS(ix,ifld)**2)
	endif
  enddo
enddo
write (*,*) '  maximum difference on fS : ',errval
if ( errval > 1.e-6) then
  call abort()
endif

! backward transform
do ifld=1,nfld
  call dfftw_execute_dft_c2r(plan_bwd, fS(:,ifld), fG(:,ifld))
enddo
fG=fG/nx ! rescale

! verify result
do ifld=1,nfld
  do ix=1,nx
    fG(ix,ifld)=fG(ix,ifld)-cos(2*PI*(ifld-1)*real(ix-1)/nx)
  enddo
enddo
errval=maxval(fG**2)
write (*,*) '  maximum difference on fG : ',errval
if ( errval > 1.e-6) then
  call abort()
endif

! destroy plans
call dfftw_destroy_plan(plan_fwd)
call dfftw_destroy_plan(plan_bwd)

write (*,*) '  All good'
write (*,*)

end subroutine test_fftw


subroutine test_fftw_batch(nx,nfld)

use fftw3_mod

implicit none

integer, intent(in) :: nx, nfld

real :: fG(nx,nfld)
real :: fS(nx+2,nfld)
real, parameter :: PI=acos(-1.)
real :: errval
integer*8 :: plan_fwd, plan_bwd
integer :: inembed,onembed
integer :: istride,idist,ostride,odist
integer :: ifld, ix

write (*,*) '-------------------'
write (*,*)
write (*,*) 'Running test_fftw_batch for nx = ',nx,'; nfld = ',nfld


! fill with sensible values
do ifld=1,nfld
  do ix=1,nx
    fG(ix,ifld)=cos(2*PI*(ifld-1)*real(ix-1)/nx)
  enddo
enddo

! create plans
inembed=0
istride=1
idist=nx
onembed=0
ostride=1
odist=nx/2+1
call dfftw_plan_many_dft_r2c(plan_fwd, 1, nx, nfld, &
 & fG, inembed, istride, idist, &
 & fS, onembed, ostride, odist, &
 & FFTW_ESTIMATE)
inembed=0
istride=1
idist=nx/2+1
onembed=0
ostride=1
odist=nx
call dfftw_plan_many_dft_c2r(plan_bwd, 1, nx, nfld, &
 & fS, inembed, istride, idist, &
 & fG, onembed, ostride, odist, &
 & FFTW_ESTIMATE)
 
 ! forward transform
call dfftw_execute_dft_r2c(plan_fwd, fG, fS)

! verify result
errval=0.
do ifld=1,nfld
  do ix=1,nx
    if (ix == 2*ifld-1) then
	  if (ifld==1) then
		errval=max(errval,(fS(ix,ifld)-nx)**2)
	  else
	    errval=max(errval,(fS(ix,ifld)-nx/2)**2)
      endif
	else
      errval=max(errval,fS(ix,ifld)**2)
	endif
  enddo
enddo
write (*,*) '  maximum difference on fS : ',errval
if ( errval > 1.e-6) then
  write (0,*) 'Difference too large!'
  !call abort()
endif

! backward transform
call dfftw_execute_dft_c2r(plan_bwd, fS, fG)
fG=fG/nx ! rescale

! verify result
do ifld=1,nfld
  do ix=1,nx
    fG(ix,ifld)=fG(ix,ifld)-cos(2*PI*(ifld-1)*real(ix-1)/nx)
  enddo
enddo
errval=maxval(fG**2)
write (*,*) '  maximum difference on fG : ',errval
if ( errval > 1.e-6) then
  write (0,*) 'Difference too large!'
  !call abort()
endif

! destroy plans
call dfftw_destroy_plan(plan_fwd)
call dfftw_destroy_plan(plan_bwd)

write (*,*) '  All good'
write (*,*)

end subroutine test_fftw_batch

subroutine test_drhook

use omp_lib
use yomhook   ,only : lhook,   dr_hook, jphook

integer :: i
real(kind=jphook) :: zhook_handle, zhook_handle_t
real :: x(100,100)

write (*,*) '-------------------'
write (*,*)
write (*,*) 'Running test_drhook'

if (lhook) call DR_HOOK('test_drhook',0,zhook_handle)
!$omp parallel do private(x,zhook_handle_t)
do i=1,50
  if (lhook) call DR_HOOK('test_drhook_1',0,zhook_handle_t)
  call random_number(x)
  x=exp(x)
  if (lhook) call DR_HOOK('test_drhook_1',1,zhook_handle_t)
enddo
!$omp end parallel do
if (lhook) call DR_HOOK('test_drhook',1,zhook_handle)

write (*,*) '  please check the drhook.prof.* files, but so far it looks ...'
write (*,*) "  if you don't get drhook.prof.* files, make sure to set export DR_HOOK=1; export DR_HOOK_OPT=prof"
write (*,*) '  All good'
write (*,*)

end subroutine test_drhook

subroutine test_fftw2d

use fftw3_mod
use yomhook   ,only : lhook,   dr_hook, jphook

integer :: ix, iy
integer, parameter :: nx=1024
integer, parameter :: ny=1024
real(kind=jphook) :: zhook_handle, zhook_handle_t
real :: fG(nx,ny), fT(nY,nx)
real :: ftemp(ny+2)
!real :: fS(nx+2,ny+2)
real :: fSx(nx+2,ny)
real :: fSy(nx,ny+2)
real :: trigsx(nx),trigsy(nx)
integer :: ifaxx(10), ifaxy(10)

integer :: iiter
integer, parameter :: niter=100
! fftw variables
integer*8 :: plan_fwd_x_loop, plan_fwd_x_batch, plan_fwd_y_loop, plan_fwd_y_batch
integer*8 :: plan_bwd_x_loop, plan_bwd_x_batch, plan_bwd_y_loop, plan_bwd_y_batch
integer :: inembed,onembed
integer :: istride,idist,ostride,odist
! fft992 variables
integer :: jbatch, lot
integer, parameter :: batch_size=1024

if (lhook) call DR_HOOK('test_fftw2d',0,zhook_handle)

write (*,*) '-------------------'
write (*,*)
write (*,*) 'Running test_fftw2d'

! create plans
inembed=0; istride=1; idist=nx; onembed=0; ostride=1; odist=nx/2+1
call dfftw_plan_many_dft_r2c(plan_fwd_x_batch, 1, nx, ny, &
 & fG, inembed, istride, idist, &
 & fSx, onembed, ostride, odist, &
 & FFTW_ESTIMATE)
call dfftw_plan_many_dft_r2c(plan_fwd_x_loop, 1, nx, 1, &
 & fG, inembed, istride, idist, &
 & fSx, onembed, ostride, odist, &
 & FFTW_ESTIMATE)
inembed=0; istride=1; idist=nx/2+1; onembed=0; ostride=1; odist=nx
call dfftw_plan_many_dft_c2r(plan_bwd_x_batch, 1, nx, ny, &
 & fSx, inembed, istride, idist, &
 & fG, onembed, ostride, odist, &
 & FFTW_ESTIMATE)
call dfftw_plan_many_dft_c2r(plan_bwd_x_loop, 1, nx, 1, &
 & fSx, inembed, istride, idist, &
 & fG, onembed, ostride, odist, &
 & FFTW_ESTIMATE)
inembed=0; istride=nx; idist=1; onembed=0; ostride=nx/2+1; odist=1
call dfftw_plan_many_dft_r2c(plan_fwd_y_batch, 1, ny, nx, &
 & fG, inembed, istride, idist, &
 & fSy, onembed, ostride, odist, &
 & FFTW_ESTIMATE)
call dfftw_plan_many_dft_r2c(plan_fwd_y_loop, 1, ny, 1, &
 & fG, inembed, istride, idist, &
 & fSy, onembed, ostride, odist, &
 & FFTW_ESTIMATE)
inembed=0; istride=nx/2+1; idist=1; onembed=0; ostride=nx; odist=1
call dfftw_plan_many_dft_c2r(plan_bwd_y_batch, 1, ny, nx, &
 & fSy, inembed, istride, idist, &
 & fG, onembed, ostride, odist, &
 & FFTW_ESTIMATE)
call dfftw_plan_many_dft_c2r(plan_bwd_y_loop, 1, ny, 1, &
 & fSy, inembed, istride, idist, &
 & fG, onembed, ostride, odist, &
 & FFTW_ESTIMATE)

fG=0.

if (lhook) call DR_HOOK('test_fftw2d:x_loop',0,zhook_handle_t)
do iiter=1,niter
  do iy=1,ny
    call dfftw_execute_dft_r2c(plan_fwd_x_loop, fG(1,iy), fSx(1,iy))
	call dfftw_execute_dft_c2r(plan_bwd_x_loop, fSx(1,iy), fG(1,iy))
  enddo
enddo
if (lhook) call DR_HOOK('test_fftw2d:x_loop',1,zhook_handle_t)

fG=0.

if (lhook) call DR_HOOK('test_fftw2d:x_batch',0,zhook_handle_t)
do iiter=1,niter
  call dfftw_execute_dft_r2c(plan_fwd_x_batch, fG, fSx)
  call dfftw_execute_dft_c2r(plan_bwd_x_batch, fSx, fG)
enddo
if (lhook) call DR_HOOK('test_fftw2d:x_batch',1,zhook_handle_t)

fG=0.

if (lhook) call DR_HOOK('test_fftw2d:y_loop',0,zhook_handle_t)
do iiter=1,niter
  do ix=1,nx
    
    call dfftw_execute_dft_r2c(plan_fwd_y_loop, fG(ix,1), ftemp)
	!fSy(ix,:)=ftemp
	!ftemp=fSy(ix,:)
	call dfftw_execute_dft_c2r(plan_bwd_y_loop, ftemp, fG(ix,1))
  enddo
enddo
if (lhook) call DR_HOOK('test_fftw2d:y_loop',1,zhook_handle_t)

fG=0.

if (lhook) call DR_HOOK('test_fftw2d:y_batch',0,zhook_handle_t)
do iiter=1,niter
  call dfftw_execute_dft_r2c(plan_fwd_y_batch, fG, fSy)
  call dfftw_execute_dft_c2r(plan_bwd_y_batch, fSy, fG)
enddo
if (lhook) call DR_HOOK('test_fftw2d:y_batch',1,zhook_handle_t)


if (lhook) call DR_HOOK('test_fftw2d:transpose',0,zhook_handle_t)
do iiter=1,niter
  do ix=1,nx
    do iy=1,ny
	  fT(iy,ix)=fG(ix,iy)
	enddo
  enddo
  do iy=1,ny
    do ix=1,nx
	  fG(ix,iy)=fT(iy,ix)
	enddo
  enddo
enddo
if (lhook) call DR_HOOK('test_fftw2d:transpose',1,zhook_handle_t)

call set99b(trigsx,ifaxx,nx)
fSx(1:nx,1:ny)=fG
if (lhook) call DR_HOOK('test_fftw2d:fft992_x',0,zhook_handle_t)
do iiter=1,niter
  ! SUBROUTINE FFT992(A,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
  do jbatch=1,ny,batch_size
    lot=min(jbatch+batch_size,ny+1)-jbatch
	!write (*,*) 'lot = ',lot
    call fft992(fSx(1,jbatch),trigsx,ifaxx,1,nx+2,nx,lot,-1)
    call fft992(fSx(1,jbatch),trigsx,ifaxx,1,nx+2,nx,lot,1)
  enddo
enddo
if (lhook) call DR_HOOK('test_fftw2d:fft992_x',1,zhook_handle_t)
fG=fSx(1:nx,1:ny)


call set99b(trigsy,ifaxy,ny)
if (lhook) call DR_HOOK('test_fftw2d:fft992_y',0,zhook_handle_t)
fSy(1:nx,1:ny)=fG
do iiter=1,niter
  ! SUBROUTINE FFT992(A,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
  do jbatch=1,nx,batch_size
    lot=min(jbatch+batch_size,nx+1)-jbatch
    !write (*,*) 'lot = ',lot
    call fft992(fSy(jbatch,1),trigsy,ifaxy,nx,1,ny,lot,-1)
    call fft992(fSy(jbatch,1),trigsy,ifaxy,nx,1,ny,lot,1)
  enddo
enddo
if (lhook) call DR_HOOK('test_fftw2d:fft992_y',1,zhook_handle_t)
fG=fSy(1:nx,1:ny)


if (lhook) call DR_HOOK('test_fftw2d',1,zhook_handle)

write (*,*) '  please check the drhook.prof.* files, but so far it looks ...'
write (*,*) "  if you don't get drhook.prof.* files, make sure to set export DR_HOOK=1; export DR_HOOK_OPT=prof"
write (*,*) '  All good'
write (*,*)

end subroutine test_fftw2d

end module tests_mod
