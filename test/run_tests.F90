program run_tests

use tests_mod

implicit none

integer :: mx, my

write (*,*)
write (*,*) 'EZTRANS test suite'
write (*,*)

! elliptic truncation test
mx=1202; my=802; call test_ellips(mx,my)
mx=82;   my=122; call test_ellips(mx,my)
mx=1202; my=802; call test_ellips(mx,my)

! fftw test
call test_fftw(nx=128,nfld=20)
call test_fftw_batch(nx=128,nfld=20)


end program run_tests
