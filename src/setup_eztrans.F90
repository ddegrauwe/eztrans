subroutine setup(config,nx,ny,nfld,nproc_A,nproc_B,truncation_order)

use mpi

use config_mod, only : config_type
use aux_mod, only : distribute, ellips

implicit none


! arguments
type(config_type), intent(inout) :: config
integer, intent(in) :: nx
integer, intent(in) :: ny
integer, intent(in) :: nfld
integer, intent(in) :: nproc_A
integer, intent(in) :: nproc_B
integer, intent(in) :: truncation_order

! local variables
integer :: ierr
integer :: jw
real, allocatable :: w(:)
integer :: jx, jy, jj, jproc, jj_l
integer :: nWM
integer :: kx_i, kx_e
character(len=1024) :: filename

! global dimensions
config%nx=nx
config%ny=ny
config%nfld=nfld
if ( truncation_order == 0 ) then
  ! no elliptic truncation
  config%mx=2*(config%nx/2)+2
  config%my=2*(config%ny/2)+2
elseif ( truncation_order == 1 ) then
  ! linear truncation
  config%mx=2*(config%nx/2)+2
  config%my=2*(config%ny/2)+2
elseif ( truncation_order == 2 ) then
  ! quadratic truncation
  config%mx=2*(config%nx/3)+2
  config%my=2*(config%ny/3)+2
elseif ( truncation_order == 3 ) then
  ! cubic truncation
  config%mx=2*(config%nx/4)+2
  config%my=2*(config%ny/4)+2
endif

! distribution parameters
config%nproc_A=nproc_A
config%nproc_B=nproc_B

! derived distribution parameters depending on runtime mpi
call mpi_comm_size(MPI_COMM_WORLD, config%nproc, ierr)
call mpi_comm_rank(MPI_COMM_WORLD, config%my_proc, ierr);
! small check
if ( config%nproc .ne. config%nproc_A*config%nproc_B ) then
  write (0,*) 'Mismatch between NPROC and NPROC_A*NPROC_B'
  call mpi_abort(MPI_COMM_WORLD,1,ierr);
endif
config%my_proc=config%my_proc+1
config%my_proc_A=modulo(config%my_proc-1,config%nproc_A)+1
config%my_proc_B=(config%my_proc-1)/config%nproc_A+1

! open output file for each task
write (filename,'(A,I2.2,A)') 'output_',config%my_proc,'.log'
open(unit=20,file=trim(filename))
write (20,*) 'Output from MPI task ',config%my_proc
write (20,*)
!write (20,*) '  nx     = ',config%nx
!write (20,*) '  ny     = ',config%ny
!write (20,*) '  mx     = ',config%mx
!write (20,*) '  my     = ',config%my

! allocations
allocate(config%nx_l(config%nproc_A))
allocate(config%jxi_l(config%nproc_A))
allocate(config%jxe_l(config%nproc_A))
allocate(config%ny_l(config%nproc_B))
allocate(config%jyi_l(config%nproc_B))
allocate(config%jye_l(config%nproc_B))
allocate(config%nfld_l(config%nproc_A))
allocate(config%jfldi_l(config%nproc_A))
allocate(config%jflde_l(config%nproc_A))
allocate(config%mx_l(config%nproc_B))
allocate(config%kxi_l(config%nproc_B))
allocate(config%kxe_l(config%nproc_B))
allocate(config%nm_l(config%nproc_B))
allocate(config%ns_l(config%nproc_A))
allocate(config%jsi_l(config%nproc_A))
allocate(config%jse_l(config%nproc_A))
allocate(config%ex(config%my))
allocate(config%ey(config%mx))


! distribution of G-space
call distribute(config%nx,config%nproc_A,config%nx_l)
config%my_nx_l=config%nx_l(config%my_proc_A)
do jproc=1,config%nproc_A
  config%jxi_l(jproc)=sum(config%nx_l(1:jproc-1))+1
  config%jxe_l(jproc)=sum(config%nx_l(1:jproc))
enddo
config%my_jxi=config%jxi_l(config%my_proc_A)
config%my_jxe=config%jxe_l(config%my_proc_A)

call distribute(config%ny,config%nproc_B,config%ny_l)
config%my_ny_l=config%ny_l(config%my_proc_B)
do jproc=1,config%nproc_B
  config%jyi_l(jproc)=sum(config%ny_l(1:jproc-1))+1
  config%jye_l(jproc)=sum(config%ny_l(1:jproc))
enddo
config%my_jyi=config%jyi_l(config%my_proc_B)
config%my_jye=config%jye_l(config%my_proc_B)

! distribution of L-space 
call distribute(config%nfld,config%nproc_A,config%nfld_l)
config%my_nfld_l=config%nfld_l(config%my_proc_A)
do jproc=1,config%nproc_A
  config%jfldi_l(jproc)=sum(config%nfld_l(1:jproc-1))+1
  config%jflde_l(jproc)=sum(config%nfld_l(1:jproc))
enddo
config%my_jfldi=config%jfldi_l(config%my_proc_A)
config%my_jflde=config%jflde_l(config%my_proc_A)

! distribution of M-space (uneven in mx)
allocate(w(config%mx/2))
do jw=1,config%mx/2
  w(jw)=0.1+sqrt(1.-(real(jw-1)/(config%mx/2-1))**2)   ! tunable !!!
enddo
call distribute(config%mx/2,config%nproc_B,config%mx_l,w)   ! factor 2 to ensure even and odd components are on the same proc
config%mx_l=2*config%mx_l
deallocate(w)
config%my_mx_l=config%mx_l(config%my_proc_B)
do jproc=1,config%nproc_B
  config%kxi_l(jproc)=sum(config%mx_l(1:jproc-1))+1
  config%kxe_l(jproc)=sum(config%mx_l(1:jproc))
enddo
config%my_kxi=config%kxi_l(config%my_proc_B)
config%my_kxe=config%kxe_l(config%my_proc_B)

! elliptic truncation
if ( truncation_order > 0) then
  call ellips(config%mx,config%my,config%ex,config%ey)
else
  config%ex(:)=config%mx
  config%ey(:)=config%my
endif

! starting and ending x-wavenumbers on this proc and total number of waves in M-space (nm_l)
kx_e=0
do jproc=1,nproc_B
  config%nm_l(jproc)=sum(config%ey(config%kxi_l(jproc):config%kxe_l(jproc)))
enddo
config%my_nm_l=config%nm_l(config%my_proc_B)

! distribution of S-space
call distribute(config%my_nm_l/4,config%nproc_A,config%ns_l) ! factor 4 to ensure even and odd components end up in the same task
config%ns_l=4*config%ns_l
config%my_ns_l=config%ns_l(config%my_proc_A)
do jproc=1,config%nproc_A
  config%jsi_l(jproc)=sum(config%ns_l(1:jproc-1))+1
  config%jse_l(jproc)=sum(config%ns_l(1:jproc))
enddo
config%my_jsi=config%jsi_l(config%my_proc_A)
config%my_jse=config%jse_l(config%my_proc_A)

! wavenumbers in M-space
allocate(config%kx_m(config%my_nm_l))
allocate(config%ky_m(config%my_nm_l))
allocate(config%jproc_m(config%my_nm_l))
if ( config%my_nm_l > 0 ) then
	jj_l=1
	jproc=1
	do jy=1,config%ey(config%my_kxi),2
	  do jx=config%my_kxi,min(config%my_kxe,config%ex(jy)),2
		config%kx_m(jj_l:jj_l+3)=jx
		config%ky_m(jj_l:jj_l+3)=jy
		jj_l=jj_l+4
	  enddo
	enddo
	jj_l=0
	do jproc=1,config%nproc_A
	  config%jproc_m(jj_l+1:jj_l+config%ns_l(jproc))=jproc
	  jj_l=jj_l+config%ns_l(jproc)
	enddo
endif

! define MPI groups
call MPI_Comm_split(MPI_COMM_WORLD, config%my_proc_B, config%my_proc_A, config%mpi_comm_A, ierr);
call MPI_Comm_split(MPI_COMM_WORLD, config%my_proc_A, config%my_proc_B, config%mpi_comm_B, ierr);

!write (20,*) '  ex     = ',config%ex
!write (20,*) '  ey     = ',config%ey
!write (20,'(X,A,9999I4)') '  kx_m   = ',config%kx_m(1:config%my_nm_l:4)
!write (20,'(X,A,9999I4)') '  ky_m   = ',config%ky_m(1:config%my_nm_l:4)
!write (20,'(X,A,9999I4)') '  jproc_m   = ',config%jproc_m(1:config%my_nm_l:4)

call flush(20)

end subroutine setup

