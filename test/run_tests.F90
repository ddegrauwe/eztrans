program run_tests

!#define TEST_FFTW
#define TEST_FFTW2D

use tests_mod

implicit none

integer :: mx, my, ierr

call mpi_init(ierr)

write (*,*)
write (*,*) 'EZTRANS test suite'
write (*,*)

#ifdef TEST_ELLIPS
! elliptic truncation test
mx=1202; my=802; call test_ellips(mx,my)
mx=82;   my=122; call test_ellips(mx,my)
mx=1202; my=802; call test_ellips(mx,my)
#endif

#ifdef TEST_FFTW
! fftw test
call test_fftw(nx=128,nfld=20)
call test_fftw_batch(nx=1024,nfld=1024)
#endif

#ifdef TEST_DRHOOK
! drhook test
call test_drhook()
#endif

#ifdef TEST_FFTW2D
call test_fftw2d()
#endif

call mpi_finalize(ierr)

end program run_tests
