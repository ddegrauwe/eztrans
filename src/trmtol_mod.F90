module trmtol_mod

use mpi

use config_mod, only: config_type

implicit none


contains

subroutine trmtol(config,fM,fL)


! transpose from M-space (y-pencils) to L-space (x-pencils)

! arguments
type(config_type), intent(in) :: config
real, intent(in)  :: fM(:,:,:)   ! (config%my_mx_l,config%ny,config%my_nfld_l)
real, intent(out) :: fL(:,:,:)   ! (config%mx,config%my_ny_l,config%my_nfld_l)


! local variables
real :: send_buffer(config%my_mx_l*config%ny*config%my_nfld_l)
real :: recv_buffer(config%mx*config%my_ny_l*config%my_nfld_l)
integer :: sendcounts(config%nproc_B)
integer :: senddispls(config%nproc_B)
integer :: recvcounts(config%nproc_B)
integer :: recvdispls(config%nproc_B)
integer :: jproc, jfld, jx, jy, offset, ierr, jjy

! calculate send counts and displacements
do jproc=1,config%nproc_B
  sendcounts(jproc)=config%my_mx_l*config%ny_l(jproc)*config%my_nfld_l
enddo
senddispls(1)=0
do jproc=2,config%nproc_B
  senddispls(jproc)=senddispls(jproc-1)+sendcounts(jproc-1)
enddo

! calculate recv counts and displacements
do jproc=1,config%nproc_B
  recvcounts(jproc)=config%mx_l(jproc)*config%my_ny_l*config%my_nfld_l
enddo
recvdispls(1)=0
do jproc=2,config%nproc_B
  recvdispls(jproc)=recvdispls(jproc-1)+recvcounts(jproc-1)
enddo

! pack send buffer
!$OMP PARALLEL DO PRIVATE(jproc,jfld,offset,jy) COLLAPSE(2)
do jproc=1,config%nproc_B
  do jfld=1,config%my_nfld_l
    offset=senddispls(jproc)+(jfld-1)*config%ny_l(jproc)*config%my_mx_l
    do jy=config%jyi_l(jproc),config%jye_l(jproc)
	  send_buffer(offset+1:offset+config%my_mx_l)=fM(1:config%my_mx_l,jy,jfld)
	  offset=offset+config%my_mx_l
	enddo
  enddo
enddo
!$OMP END PARALLEL DO

! communications
call mpi_alltoallv(send_buffer, sendcounts, senddispls, MPI_FLOAT, &
 & recv_buffer, recvcounts, recvdispls, MPI_FLOAT, &
 & config%mpi_comm_B, ierr)

! unpack recv buffer
fL(:,:,:)=0.
!$OMP PARALLEL DO PRIVATE(jproc,jfld,offset,jy) COLLAPSE(2)
do jproc=1,config%nproc_B
  do jfld=1,config%my_nfld_l
    offset=recvdispls(jproc)+(jfld-1)*config%my_ny_l*config%mx_l(jproc)
    do jy=1,config%my_ny_l
	  fL(config%kxi_l(jproc):config%kxe_l(jproc),jy,jfld)=recv_buffer(offset+1:offset+config%mx_l(jproc))
	  offset=offset+config%mx_l(jproc)
	enddo
  enddo
enddo
!$OMP END PARALLEL DO

end subroutine trmtol

end module trmtol_mod
