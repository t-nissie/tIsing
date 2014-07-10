! tIsing.F - two dimensional Ising model
! Author: NISHIMATSU Takeshi
! Time-stamp: <2014-07-10 17:07:54 takeshi>
! References:
!  http://www.cp.tn.tudelft.nl/research/ising/
!  http://oscar.cacr.caltech.edu/Hrothgar/Ising/references.html
!  http://www-cms.phys.s.u-tokyo.ac.jp/~naoki/CIPINTRO/CIP/ising.html
!!
#define N_MEASUREMENTS 2
#define N_CORRELATION 5
program isingmodel
  implicit none
  integer :: seed, relax, samples, n
  double precision :: T,T_from,T_to,T_by   ! temperature (accurately, kT/J)
  integer, allocatable :: s(:)             ! size -n:n*n+n-1 array for n*n sigma
  double precision :: measurements(N_MEASUREMENTS+N_CORRELATION)
  double precision :: sum1(N_MEASUREMENTS+N_CORRELATION)
  double precision ::  avg(N_MEASUREMENTS+N_CORRELATION)
  double precision :: ran1, dummy
  external ran1
  integer i, j
  character(20) :: fmt
  write(fmt,'(a,i1,a)') '(f5.3,', N_MEASUREMENTS+N_CORRELATION, 'e14.6)'

  ! read a few parameters
  read(5,*) seed
  read(5,*) relax
  read(5,*) samples
  read(5,*) n
  read(5,*) T_from
  read(5,*) T_to
  read(5,*) T_by
  write(6,'(a)')           '#   seed, relax, samples, n, T_from, T_to, T_by'
  write(6,'(a,4i7,3f6.2)') '#', seed, relax, samples, n, T_from, T_to, T_by

  dummy = ran1(seed) ! From Numerical Recipes, 1986, p. 196

  allocate(s(-n:n*n+n*N_CORRELATION-1))
  s = 1.0   ! initial magnetization 1.0 or -1.0

  T = T_from
  do while (T<T_to)
     T = T + T_by
     do i = 1, relax
        call metropolis(n,T,s)
     end do
     sum1 = 0.0
     do i = 1, samples
        call metropolis(n,T,s)
        call evaluate(n,s,measurements)
        sum1 = sum1 + measurements
     end do
     avg = sum1/samples
     write(6,fmt) T,(avg(j),j=1,N_MEASUREMENTS+N_CORRELATION)
  end do
end program isingmodel

subroutine evaluate(n,s,measurements)
  implicit none
  integer :: n, s(-n:*)
  double precision :: measurements(N_MEASUREMENTS+N_CORRELATION)
  integer :: i,j,sum
  measurements=0.0
  do i = 0, n*n-1    ! magnetization
     measurements(1) = measurements(1) + s(i)
  end do
  measurements(1) = measurements(1)/n/n
  do i = 0, n*n-1   ! energy
     measurements(2) = measurements(2) - s(i)*(s(i+1)+s(i+n))
  end do
  measurements(2) = 2*measurements(2)/n/n
  do j = 1, N_CORRELATION   ! correlation functions
     do i = 0, n*n-1
        measurements(N_MEASUREMENTS+j) = measurements(N_MEASUREMENTS+j) + s(i)*(s(i+j)+s(i+j*n))
     end do
     measurements(N_MEASUREMENTS+j) = measurements(N_MEASUREMENTS+j)/4/n/n - measurements(1)**2/2
  end do
end subroutine evaluate

subroutine metropolis(n,T,s)
  implicit none
  integer :: n
  double precision :: T
  integer :: s(-n:*)
  double precision :: ran1
  external ran1
  double precision :: w(-4:4)
  integer :: i, j

  w(-4) = 1.0
  w(-2) = 1.0
  w( 0) = 1.0
  w( 2) = exp(-4/T)
  w( 4) = exp(-8/T)

  do i = 0, n*n-1
     j = s(i) * (s(i-n)+s(i-1)+s(i+1)+s(i+n))
     if (w(j).eq.1.0 .or. ran1(1).lt.w(j)) s(i)=-s(i)
  end do

  ! copy top and bottom edges
  do i = 0, n*N_CORRELATION-1
     s(n*n+i) = s(i)
  end do
  do i = 0, n-1
     s(-n+i)  = s(n*n-n+i)
  end do
end subroutine metropolis
