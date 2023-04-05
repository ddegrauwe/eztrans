module trstom_mod

use mpi

use config_mod, only: config_type

implicit none


contains

subroutine trstom(config,fS,fM)


! transpose from S-space (z-pencils) to M-space (y-pencils)

! arguments
type(config_type), intent(in) :: config
real, intent(in)  :: fS(:,:)     ! config%my_ns_l,config%nfld)
real, intent(out) :: fM(:,:,:)   ! (config%my_mx_l,config%my,config%my_nfld_l)


! local variables
real :: send_buffer(config%my_ns_l*config%nfld)
real :: recv_buffer(config%my_nm_l*config%my_nfld_l)
integer :: sendcounts(config%nproc_A)
integer :: senddispls(config%nproc_A)
integer :: recvcounts(config%nproc_A)
integer :: recvdispls(config%nproc_A)
integer :: jproc, jfld, jx, jy, js, offset, jj
integer :: ierr
integer :: jproc_l(config%nfld)
integer :: jfld_l(config%nfld)

! calculate send counts and displacements
do jproc=1,config%nproc_A
  sendcounts(jproc)=config%my_ns_l*config%nfld_l(jproc)
enddo
senddispls(1)=0
do jproc=2,config%nproc_A
  senddispls(jproc)=senddispls(jproc-1)+sendcounts(jproc-1)
enddo

! calculate recv counts and displacements
do jproc=1,config%nproc_A
  recvcounts(jproc)=config%ns_l(jproc)*config%my_nfld_l
enddo
recvdispls(1)=0
do jproc=2,config%nproc_A
  recvdispls(jproc)=recvdispls(jproc-1)+recvcounts(jproc-1)
enddo

! pack send buffer
! fill arrays with local proc and fld indices; this is necessary because OpenMP collapsing doesn't work on non-rectangular loops
do jproc=1,config%nproc_A
  jproc_l(config%jfldi_l(jproc):config%jflde_l(jproc))=jproc
  jfld_l(config%jfldi_l(jproc):config%jflde_l(jproc))=(/ (jfld, jfld=1, config%nfld_l(jproc)) /)
enddo

! unpack recv buffer
!$OMP PARALLEL DO PRIVATE(jfld,offset,js)
do jfld=1,config%nfld
  offset=senddispls(jproc_l(jfld))+(jfld_l(jfld)-1)*config%my_ns_l
  do js=1,config%my_ns_l
    send_buffer(offset+js)=fS(js,jfld)
  enddo
enddo
!$OMP END PARALLEL DO

! communications
call mpi_alltoallv(send_buffer, sendcounts, senddispls, MPI_FLOAT, &
 & recv_buffer, recvcounts, recvdispls, MPI_FLOAT, &
 & config%mpi_comm_A, ierr)

! unpack recv buffer
fM(:,:,:)=0.
!$OMP PARALLEL DO PRIVATE(jproc,jfld,offset,jx,jy,js) COLLAPSE(2)
do jproc=1,config%nproc_A
  do jfld=1,config%my_nfld_l
    offset=recvdispls(jproc)+(jfld-1)*config%ns_l(jproc)-config%jsi_l(jproc)+1
    do js=config%jsi_l(jproc),config%jse_l(jproc),4
	   jx=config%kx_m(js)
	   jy=config%ky_m(js)
	   fM(jx-config%my_kxi+1,jy+0,jfld)=recv_buffer(offset+js+0)
	   fM(jx-config%my_kxi+2,jy+0,jfld)=recv_buffer(offset+js+1)
	   fM(jx-config%my_kxi+1,jy+1,jfld)=recv_buffer(offset+js+2)
	   fM(jx-config%my_kxi+2,jy+1,jfld)=recv_buffer(offset+js+3)
	enddo
  enddo
enddo
!$OMP END PARALLEL DO

end subroutine trstom

end module trstom_mod
