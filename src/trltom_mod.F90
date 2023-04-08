module trltom_mod

use mpi
use yomhook   ,only : lhook,   dr_hook, jphook

use config_mod, only: config_type

implicit none


contains

subroutine trltom(config,fL,fM)


! transpose from L-space (x-pencils) to M-space (y-pencils)

! arguments
type(config_type), intent(in) :: config
real, intent(in)  :: fL(:,:,:)   ! (config%mx,config%my_ny_l,config%my_nfld_l)
real, intent(out) :: fM(:,:,:)   ! (config%my_mx_l,config%ny,config%my_nfld_l)


! local variables
real :: send_buffer(config%mx*config%my_ny_l*config%my_nfld_l)
real :: recv_buffer(config%my_mx_l*config%ny*config%my_nfld_l)
integer :: sendcounts(config%nproc_B)
integer :: senddispls(config%nproc_B)
integer :: recvcounts(config%nproc_B)
integer :: recvdispls(config%nproc_B)
integer :: jproc, jfld, jx, jy, offset, ierr, jjy
real(kind=jphook) :: zhook_handle, zhook_handle_t

if (lhook) call DR_HOOK('trltom',0,zhook_handle)

! calculate send counts and displacements
do jproc=1,config%nproc_B
  sendcounts(jproc)=config%mx_l(jproc)*config%my_ny_l*config%my_nfld_l
enddo
senddispls(1)=0
do jproc=2,config%nproc_B
  senddispls(jproc)=senddispls(jproc-1)+sendcounts(jproc-1)
enddo

! calculate recv counts and displacements
do jproc=1,config%nproc_B
  recvcounts(jproc)=config%my_mx_l*config%ny_l(jproc)*config%my_nfld_l
enddo
recvdispls(1)=0
do jproc=2,config%nproc_B
  recvdispls(jproc)=recvdispls(jproc-1)+recvcounts(jproc-1)
enddo

! pack send buffer
!$OMP PARALLEL PRIVATE(jproc,jfld,offset,jy,zhook_handle_t)
if (lhook) call DR_HOOK('trltom:send_buffer',0,zhook_handle_t)
!$OMP DO COLLAPSE(2)
do jproc=1,config%nproc_B
  do jfld=1,config%my_nfld_l
    offset=senddispls(jproc)+(jfld-1)*config%my_ny_l*config%mx_l(jproc)
    do jy=1,config%my_ny_l
	  send_buffer(offset+1:offset+config%mx_l(jproc))=fL(config%kxi_l(jproc):config%kxe_l(jproc),jy,jfld)
	  offset=offset+config%mx_l(jproc)
	enddo
  enddo
enddo
!$OMP END DO
if (lhook) call DR_HOOK('trltom:send_buffer',1,zhook_handle_t)
!$OMP END PARALLEL

! communications
call mpi_alltoallv(send_buffer, sendcounts, senddispls, MPI_FLOAT, &
 & recv_buffer, recvcounts, recvdispls, MPI_FLOAT, &
 & config%mpi_comm_B, ierr)
 
! unpack recv buffer
fM(:,:,:)=0.

!$OMP PARALLEL PRIVATE(jproc,jfld,offset,jy,zhook_handle_t)
if (lhook) call DR_HOOK('trltom:recv_buffer',0,zhook_handle_t)
!$OMP DO COLLAPSE(2)
do jproc=1,config%nproc_B
  do jfld=1,config%my_nfld_l
    offset=recvdispls(jproc)+(jfld-1)*config%ny_l(jproc)*config%my_mx_l
    do jy=config%jyi_l(jproc),config%jye_l(jproc)
	  fM(1:config%my_mx_l,jy,jfld)=recv_buffer(offset+1:offset+config%my_mx_l)
	  offset=offset+config%my_mx_l
	enddo
  enddo
enddo
!$OMP END DO
if (lhook) call DR_HOOK('trltom:recv_buffer',1,zhook_handle_t)
!$OMP END PARALLEL

if (lhook) call DR_HOOK('trltom',1,zhook_handle)

end subroutine trltom

end module trltom_mod
