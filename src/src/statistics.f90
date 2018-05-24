module statistics
  use constants,only:dp
  implicit none

  contains

  ! calculates the mean of a vector
  function mean(x) result(xm)
    real(dp)::xm
    real(dp),dimension(:),intent(in)::x
    xm=sum(x)/size(x)

    return
  end function mean

  ! calculates the weighted mean
  function wmean(x,wx) result(xm)
    real(dp)::xm
    real(dp),dimension(:),intent(in)::x,wx
    xm=sum(x*wx)/sum(wx)

    return
  end function wmean

end module statistics
