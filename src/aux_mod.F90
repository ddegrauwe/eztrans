module aux_mod

implicit none

contains

subroutine distribute(n_g,nProc,n_l,w)

implicit none

integer, intent(in) :: n_g
integer, intent(in) :: nProc
integer, intent(out) :: n_l(nProc)
real, optional, target, intent(in) :: w(n_g)

integer :: jj
real, pointer :: wp(:)
real :: wt, wa
real, parameter :: wtol=1.e-3

integer :: iproc
if (present(w)) then
  wp=>w
else
  allocate(wp(n_g))
  wp(:)=1.
endif

! total weight
wt=sum(wp)

! loop over points, adding to different procs as long as wa is smaller than 
wa=0.
iproc=1
n_l(:)=0
do jj=1,n_g,1
  n_l(iproc)=n_l(iproc)+1
  wa=wa+wp(jj)
  if ( wa .gt. iproc*wt/nProc - wtol ) then
	iproc=iproc+1
  endif
enddo


! verify if final distribution is okay
!jj=0
!write (*,*) 'weight distribution:'
!do iproc=1,nproc
!  write (*,*) sum(wp(jj+1:jj+n_l(iproc)))
!  jj=jj+n_l(iproc)
!enddo

if (.not. present(w) ) deallocate(wp)

end subroutine distribute

subroutine ellips(mx,my,ex,ey)
! arguments
integer, intent(in) :: mx
integer, intent(in) :: my
integer, intent(out) :: ex(my)
integer, intent(out) :: ey(mx)
! local variables
integer :: jx, jy, eys

do jx=1,mx/2
  ey(2*jx-1)=2*(floor((my/2-1)*sqrt(1.-(real(jx-1)/(mx/2-1))**2)+1))
  ey(2*jx)=ey(2*jx-1)
enddo
ey(1:2)=my

! alternative calculation enforcing consistency (otherwise round-off errors may lead to crashes)
!ex(:)=0
!do jx=1,mx
!  do jy=1,ey(jx)
!    ex(jy)=max(ex(jy),jx)
!  enddo
!enddo

! same, but avoiding double loop
eys=my
do jx=1,mx
  !write (*,*) 'jx = ',jx,'; eys = ',eys,' ey(jx) = ',ey(jx)
  if ( ey(jx) < eys ) then
    !write (*,*) 'BOINK! ','setting values ',ey(jx)+1,':',eys,' to ',jx-1
    ex(ey(jx)+1:eys)=jx-1
    eys=ey(jx)
  endif
enddo
ex(1:eys)=mx

!do jy=1,my/2
!  ex(2*jy-1)=2*(floor((mx/2-1)*sqrt(1.-(real(jy-1)/(my/2-1))**2)+1.e-6)+1)
!  ex(2*jy)=ex(2*jy-1)
!enddo


end subroutine ellips

end module aux_mod