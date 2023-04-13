module trgtol_mod

use mpi
use yomhook   ,only : lhook,   dr_hook, jphook

use config_mod, only: config_type

implicit none

contains

subroutine trgtol(config,fG,fL)

! transpose from G-space (z-pencils) to L-space (x-pencils)

! arguments
type(config_type), intent(in) :: config
real, intent(in)  :: fG(config%my_nx_l,config%my_ny_l,config%nfld)
real, intent(out) :: fL(config%nx+2,config%my_ny_l,config%my_nfld_l)

! local variables
!real :: send_buffer(config%nx_l*config%ny_l*config%nfld)   ! fG can be used as send buffer directly
real :: recv_buffer(config%nx*config%my_ny_l*config%my_nfld_l)
integer :: sendcounts(config%nproc_A)
integer :: senddispls(config%nproc_A)
integer :: recvcounts(config%nproc_A)
integer :: recvdispls(config%nproc_A)
integer :: jproc, jfld, jx, jy, offset, ierr
real(kind=jphook) :: zhook_handle, zhook_handle_t

if (lhook) call DR_HOOK('trgtol',0,zhook_handle)

! calculate send counts and displacements
do jproc=1,config%nproc_A
  sendcounts(jproc)=config%my_nx_l*config%my_ny_l*config%nfld_l(jproc)
enddo
senddispls(1)=0
do jproc=2,config%nproc_A
  senddispls(jproc)=senddispls(jproc-1)+sendcounts(jproc-1)
enddo

! calculate recv counts and displacements
do jproc=1,config%nproc_A
  recvcounts(jproc)=config%nx_l(jproc)*config%my_ny_l*config%my_nfld_l
enddo
recvdispls(1)=0
do jproc=2,config%nproc_A
  recvdispls(jproc)=recvdispls(jproc-1)+recvcounts(jproc-1)
enddo

! communications
call mpi_alltoallv(fG, sendcounts, senddispls, MPI_DOUBLE_PRECISION, &
 & recv_buffer, recvcounts, recvdispls, MPI_DOUBLE_PRECISION, &
 & config%mpi_comm_A, ierr)
 
! unpack recv buffer
!$OMP PARALLEL PRIVATE(jproc,jfld,offset,jy,zhook_handle_t)
if (lhook) call DR_HOOK('trgtol:unpack',0,zhook_handle_t)
!$OMP DO COLLAPSE(2)
do jproc=1,config%nproc_A
  do jfld=1,config%my_nfld_l
    offset=recvdispls(jproc)+(jfld-1)*config%my_ny_l*config%nx_l(jproc)
    do jy=1,config%my_ny_l
      fL(config%jxi_l(jproc):config%jxe_l(jproc),jy,jfld)=recv_buffer(offset+1:offset+config%nx_l(jproc))
	  offset=offset+config%nx_l(jproc)
	enddo
  enddo
enddo
!$OMP END DO
if (lhook) call DR_HOOK('trgtol:unpack',1,zhook_handle_t)
!$OMP END PARALLEL

if (lhook) call DR_HOOK('trgtol',1,zhook_handle)

end subroutine trgtol

end module trgtol_mod
