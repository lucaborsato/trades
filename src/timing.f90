module timing
    use constants
    !$ use omp_lib
    implicit none

  contains

  ! ------------------------------------------------------------------ !
  ! computes the time in hours minutes seconds between t2 and t1
  subroutine timer(t1,t2,h,m,s)
    real(dp), intent(in) :: t1,t2
    integer, intent(out) :: h,m
    real(dp),intent(out)::s
    real(dp) :: tots
    integer :: totm

    tots=(t2-t1)/ncpu
    totm=int(tots/60._dp)
    h=int(totm/60._dp)
    m=totm-h*int(60)
    s=tots-totm*60._dp

    return
  end subroutine timer
  ! ------------------------------------------------------------------ !

end module timing
