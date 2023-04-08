program driver

use yomhook   ,only : lhook,   dr_hook, jphook
use config_mod, only : config_type
use trgtol_mod, only : trgtol
use trltom_mod, only : trltom
use trmtos_mod, only : trmtos
use trstom_mod, only : trstom
use trmtol_mod, only : trmtol
use trltog_mod, only : trltog

implicit none

type(config_type) :: config
real, allocatable :: fG(:,:,:)
real, allocatable :: fL(:,:,:)
real, allocatable :: fM(:,:,:)
real, allocatable :: fS(:,:)
integer :: nx, ny, nfld
integer :: nproc_A, nproc_B
integer :: truncation_order
integer :: jfld, jx, jy
integer :: ierr
real(kind=jphook) :: zhook_handle, zhook_handle_t

call mpi_init(ierr)

if (lhook) call DR_HOOK('driver',0,zhook_handle)

nx=80
ny=60
nfld=10
nproc_A=1
nproc_B=1
truncation_order=1

call eztrans_setup(config,nx,ny,nfld,nproc_A,nproc_B,truncation_order)

allocate(fG(config%my_nx_l,config%my_ny_l,config%nfld))
allocate(fL(config%nx+2,config%my_ny_l,config%my_nfld_l))
allocate(fM(config%my_mx_l,config%ny+2,config%my_nfld_l))
allocate(fS(config%my_ns_l,config%nfld))
fG=0.
fL=0.
fM=0.
fS=0.

do jfld=1,config%nfld
  do jy=1,config%my_ny_l
    do jx=1,config%my_nx_l
	  fG(jx,jy,jfld)=10000*jfld+100*(sum(config%ny_l(1:config%my_proc_B-1))+jy)+sum(config%nx_l(1:config%my_proc_A-1))+jx
	enddo
  enddo
enddo

#undef DEBUG

#ifndef DEBUG
write (20,*) 'fG = '
do jfld=1,config%nfld
  do jy=1,config%my_ny_l
    write (20,'(999F8.0)') fG(:,jy,jfld)
  enddo
  write (20,*)
enddo
#endif

call trgtol(config,fG,fL)

#ifdef DEBUG
write (20,*) 'fL = '
do jfld=1,config%my_nfld_l
  do jy=1,config%my_ny_l
    write (20,'(999F8.0)') fL(:,jy,jfld)
  enddo
  write (20,*)
enddo
#endif

call trltom(config,fL,fM)

#ifdef DEBUG
write (20,*) 'fM = '
do jfld=1,config%my_nfld_l
  do jy=1,config%ny+2
    write (20,'(999F8.0)') fM(1:config%my_mx_l,jy,jfld)
  enddo
  write (20,*)
enddo
#endif

call trmtos(config,fM,fS)

#ifdef DEBUG
write (20,*) 'fS = '
do jfld=1,config%nfld
  write (20,'(999F8.0)') fS(:,jfld)
  write (20,*)
enddo
#endif

call trstom(config,fS,fM)

#ifdef DEBUG
write (20,*) 'fM = '
do jfld=1,config%my_nfld_l
  do jy=1,config%ny+2
    write (20,'(999F8.0)') fM(1:config%my_mx_l,jy,jfld)
  enddo
  write (20,*)
enddo
#endif

call trmtol(config,fM,fL)

#ifdef DEBUG
write (20,*) 'fL = '
do jfld=1,config%my_nfld_l
  do jy=1,config%my_ny_l
    write (20,'(999F8.0)') fL(:,jy,jfld)
  enddo
  write (20,*)
enddo
#endif

call trltog(config,fL,fG)

#ifndef DEBUG
write (20,*) 'fG = '
do jfld=1,config%nfld
  do jy=1,config%my_ny_l
    write (20,'(999F8.0)') fG(:,jy,jfld)
  enddo
  write (20,*)
enddo
#endif

deallocate(fG)
deallocate(fL)
deallocate(fM)
deallocate(fS)

call eztrans_end(config)

if (lhook) call DR_HOOK('driver',1,zhook_handle)

call mpi_finalize(ierr)

end program driver
