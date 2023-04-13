module trltog_mod

use mpi
use yomhook   ,only : lhook,   dr_hook, jphook

use config_mod, only: config_type

implicit none


contains

subroutine trltog(config,fL,fG)


! transpose from L-space (x-pencils) to G-space (z-pencils)

! arguments
type(config_type), intent(in) :: config
real, intent(in)  :: fL(config%nx+2,config%my_ny_l,config%my_nfld_l)
real, intent(out) :: fG(config%my_nx_l,config%my_ny_l,config%nfld)


! local variables
real :: send_buffer(config%nx*config%my_ny_l*config%my_nfld_l)
!real :: recv_buffer(config%nx_l*config%ny_l*config%nfld)   ! fG can be used as recv buffer directly
integer :: sendcounts(config%nproc_A)
integer :: senddispls(config%nproc_A)
integer :: recvcounts(config%nproc_A)
integer :: recvdispls(config%nproc_A)
integer :: jproc, jfld, jx, jy, offset, ierr
real(kind=jphook) :: zhook_handle, zhook_handle_t

if (lhook) call DR_HOOK('trltog',0,zhook_handle)

! calculate send counts and displacements
do jproc=1,config%nproc_A
  sendcounts(jproc)=config%nx_l(jproc)*config%my_ny_l*config%my_nfld_l
enddo
senddispls(1)=0
do jproc=2,config%nproc_A
  senddispls(jproc)=senddispls(jproc-1)+sendcounts(jproc-1)
enddo

! calculate recv counts and displacements
do jproc=1,config%nproc_A
  recvcounts(jproc)=config%my_nx_l*config%my_ny_l*config%nfld_l(jproc)
enddo
recvdispls(1)=0
do jproc=2,config%nproc_A
  recvdispls(jproc)=recvdispls(jproc-1)+recvcounts(jproc-1)
enddo

! pack send buffer
!$OMP PARALLEL PRIVATE(jproc,jfld,offset,jy,zhook_handle_t)
if (lhook) call DR_HOOK('trltog:pack',0,zhook_handle_t)
!$OMP DO COLLAPSE(2)
do jproc=1,config%nproc_A
  do jfld=1,config%my_nfld_l
    offset=senddispls(jproc)+(jfld-1)*config%my_ny_l*config%nx_l(jproc)
    do jy=1,config%my_ny_l
      send_buffer(offset+1:offset+config%nx_l(jproc))=fL(config%jxi_l(jproc):config%jxe_l(jproc),jy,jfld)
	  offset=offset+config%nx_l(jproc)
	enddo
  enddo
enddo
!$OMP END DO
if (lhook) call DR_HOOK('trltog:pack',1,zhook_handle_t)
!$OMP END PARALLEL

! communications
call mpi_alltoallv(send_buffer, sendcounts, senddispls, MPI_DOUBLE_PRECISION, &
 & fG, recvcounts, recvdispls, MPI_DOUBLE_PRECISION, &
 & config%mpi_comm_A, ierr)
 
! no need to unpack recv buffer: written directly to fG

if (lhook) call DR_HOOK('trltog',1,zhook_handle)

end subroutine trltog

end module trltog_mod
