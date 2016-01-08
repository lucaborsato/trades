!
! --------------------------------------------
! Module of Quasi-Monte Carlo method utilities
! used by PSO: opti_pso.f90
! --------------------------------------------
!
!----------------------------------------------------------------------
module util_qmc
  use constants,only:dp
  private

  !============= PUBLIC ROUTINES =============
  public :: hammersley       ! hammersley point set
  public :: van_der_corput   ! van der Corput sequence


  !============= dpIVATE VARIABLES =============
  integer, parameter :: max_dim = 100   ! maximum dimension
  integer, parameter :: primes(1:max_dim) = (/ &
&                       2, 3, 5, 7, 11, 13, 17, 19, 23, 29, &
&                       31, 37, 41, 43, 47, 53, 59, 61, 67, 71, &
&                       73, 79, 83, 89, 97, 101, 103, 107, 109, 113, &
&                       127, 131, 137, 139, 149, 151, 157, 163, 167, 173, &
&                       179, 181, 191, 193, 197, 199, 211, 223, 227, 229, &
&                       233, 239, 241, 251, 257, 263, 269, 271, 277, 281, &
&                       283, 293, 307, 311, 313, 317, 331, 337, 347, 349, &
&                       353, 359, 367, 373, 379, 383, 389, 397, 401, 409, &
&                       419, 421, 431, 433, 439, 443, 449, 457, 461, 463, &
&                       467, 479, 487, 491, 499, 503, 509, 521, 523, 541 &
&                    /)

  contains
!======================================================================
!========================== PUBLIC ROUTINES ===========================

!----------------------------------------------------------------------
! Hammersley Point Set
!----------------------------------------------------------------------
subroutine hammersley(n, dim, x)
  integer, intent(in)  :: n
  integer, intent(in)  :: dim
  real(dp),    intent(out) :: x(:)
  integer :: i, j, jj
  if(dim > max_dim) then
    print *, 'hammersley: dimension too large'
    stop
  end if
  do j=0, n-1
    jj = j * dim
    x(jj + 1) = real(j,dp) / real(n,dp)
    do i=2, dim
      x(jj + i) = van_der_corput(j, primes(i - 1))
    end do
  end do
end subroutine hammersley

!----------------------------------------------------------------------
! van der Corput Sequence
!----------------------------------------------------------------------
function van_der_corput(i, base) result(r)
  integer, intent(in) :: i
  integer, intent(in) :: base
  real(dp)                :: r
  real(dp)    :: f, factor
  integer :: ii
  f = 1.0_dp / real(base,dp)
  factor = f
  ii = i
  r = 0.0_dp
  do while(ii > 0)
    r = r + real(mod(ii, base),dp) * factor
    ii = ii / base
    factor = factor * f
  end do
end function van_der_corput

end module util_qmc
