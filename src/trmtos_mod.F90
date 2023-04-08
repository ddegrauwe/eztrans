module trmtos_mod

use mpi
use yomhook   ,only : lhook,   dr_hook, jphook

use config_mod, only: config_type

implicit none


contains

subroutine trmtos(config,fM,fS)


! transpose from M-space (y-pencils) to S-space (z-pencils)

! arguments
type(config_type), intent(in) :: config
real, intent(in)  :: fM(:,:,:)   ! (config%my_mx_l,config%my,config%my_nfld_l)
real, intent(out) :: fS(:,:)     ! (config%my_ns_l,config%nfld)


! local variables
real :: send_buffer(config%my_nm_l*config%my_nfld_l)
real :: recv_buffer(config%my_ns_l*config%nfld)
integer :: sendcounts(config%nproc_A)
integer :: senddispls(config%nproc_A)
integer :: recvcounts(config%nproc_A)
integer :: recvdispls(config%nproc_A)
integer :: jproc, jfld, jx, jy, js, offset, jj
integer :: ierr
integer :: jproc_l(config%nfld)
integer :: jfld_l(config%nfld)
real(kind=jphook) :: zhook_handle, zhook_handle_t

if (lhook) call DR_HOOK('trmtos',0,zhook_handle)

! calculate send counts and displacements
do jproc=1,config%nproc_A
  sendcounts(jproc)=config%ns_l(jproc)*config%my_nfld_l
enddo
senddispls(1)=0
do jproc=2,config%nproc_A
  senddispls(jproc)=senddispls(jproc-1)+sendcounts(jproc-1)
enddo

! calculate recv counts and displacements
do jproc=1,config%nproc_A
  recvcounts(jproc)=config%my_ns_l*config%nfld_l(jproc)
enddo
recvdispls(1)=0
do jproc=2,config%nproc_A
  recvdispls(jproc)=recvdispls(jproc-1)+recvcounts(jproc-1)
enddo

! pack send buffer
!$OMP PARALLEL PRIVATE(jproc,jfld,offset,js,jx,jy,zhook_handle_t)
if (lhook) call DR_HOOK('trmtos:send_bufr',0,zhook_handle_t)
!$OMP DO COLLAPSE(2)
do jproc=1,config%nproc_A
  do jfld=1,config%my_nfld_l
    offset=senddispls(jproc)+(jfld-1)*config%ns_l(jproc)-config%jsi_l(jproc)+1
    do js=config%jsi_l(jproc),config%jse_l(jproc),4
	   jx=config%kx_m(js)
	   jy=config%ky_m(js)
	   send_buffer(offset+js+0)=fM(jx-config%my_kxi+1,jy,jfld)
	   send_buffer(offset+js+1)=fM(jx-config%my_kxi+2,jy,jfld)
	   send_buffer(offset+js+2)=fM(jx-config%my_kxi+1,jy+1,jfld)
	   send_buffer(offset+js+3)=fM(jx-config%my_kxi+2,jy+1,jfld)
	enddo
  enddo
enddo
!$OMP END DO
if (lhook) call DR_HOOK('trmtos:send_bufr',1,zhook_handle_t)
!$OMP END PARALLEL

! communications
call mpi_alltoallv(send_buffer, sendcounts, senddispls, MPI_FLOAT, &
 & recv_buffer, recvcounts, recvdispls, MPI_FLOAT, &
 & config%mpi_comm_A, ierr)

! fill arrays with local proc and fld indices; this is necessary because OpenMP collapsing doesn't work on non-rectangular loops
do jproc=1,config%nproc_A
  jproc_l(config%jfldi_l(jproc):config%jflde_l(jproc))=jproc
  jfld_l(config%jfldi_l(jproc):config%jflde_l(jproc))=(/ (jfld, jfld=1, config%nfld_l(jproc)) /)
enddo

! unpack recv buffer
!$OMP PARALLEL PRIVATE(jfld,offset,js,zhook_handle_t)
if (lhook) call DR_HOOK('trmtos:recv_bufr',0,zhook_handle_t)
!$OMP DO
do jfld=1,config%nfld
  offset=recvdispls(jproc_l(jfld))+(jfld_l(jfld)-1)*config%my_ns_l
  do js=1,config%my_ns_l
    fS(js,jfld)=recv_buffer(offset+js)
  enddo
enddo
!$OMP END DO
if (lhook) call DR_HOOK('trmtos:recv_bufr',1,zhook_handle_t)
!$OMP END PARALLEL

if (lhook) call DR_HOOK('trmtos',1,zhook_handle)

end subroutine trmtos

end module trmtos_mod
