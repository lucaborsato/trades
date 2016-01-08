module gaussian
  use constants,only:dp
  implicit none

  interface gaussdev
    module procedure gaussdev_s,gaussdev_v
  end interface gaussdev

  contains

  ! Based on Numerical Recipes
  ! it calculates a gaussian number of kind N(0,1)
  subroutine gaussdev_s(gauss)
    real(dp),intent(out)::gauss
    real(dp),dimension(2)::v1
    real(dp)::rsq,lrsq,lrsq2,slrsq2

    do
      call random_number(v1)
      v1 = 2._dp * v1 - 1._dp
      rsq = v1(1)**2 + v1(2)**2
      if( (rsq.gt.0._dp).and.(rsq.lt.1._dp) ) exit
    end do
    lrsq = log(rsq)
    lrsq2 = lrsq/rsq
    slrsq2 = sqrt( -2._dp * lrsq2 )
    gauss = v1(1) * slrsq2

    return
  end subroutine gaussdev_s

  ! Based on Numerical Recipes
  ! it calculates a gaussian vector of kind N(0,1)
  subroutine gaussdev_v(gauss)
    real(dp),dimension(:),intent(out)::gauss
    real(dp),dimension(size(gauss))::v1,v2,rsq,lrsq,lrsq2,slrsq2
    integer::n,ng,nn,m
    logical,dimension(size(gauss))::mask

    n=size(gauss)
    ng=1
    do
      if(ng.gt.n) exit
      call random_number(v1(ng:n))
      call random_number(v2(ng:n))
      v1(ng:n) = 2._dp * v1(ng:n) - 1._dp
      v2(ng:n) = 2._dp * v2(ng:n) - 1._dp
      rsq(ng:n) = v1(ng:n)**2 + v2(ng:n)**2
      mask(ng:n) = ( (rsq(ng:n).gt.0._dp) .and. (rsq(ng:n).lt.1._dp) )
      call array_copy(pack(v1(ng:n),mask(ng:n)),v1(ng:),nn,m)
      v2(ng:ng+nn-1) = pack(v2(ng:n),mask(ng:n))
      rsq(ng:ng+nn-1) = pack(rsq(ng:n),mask(ng:n))
      ng = ng+nn
    end do
    lrsq = log(rsq)
    lrsq2 = lrsq/rsq
    slrsq2 = sqrt(-2._dp * lrsq2)
    gauss = v1 * slrsq2

    return
  end subroutine gaussdev_v

  ! Numerical Recipes
  ! it copies a masked array into another, needed by gaussdev_v
  subroutine array_copy(src,dest,n_copied,n_not_copied)
    real(dp),dimension(:),intent(in)::src
    real(dp),dimension(:),intent(out)::dest
    integer,intent(out)::n_copied,n_not_copied

    n_copied = min(size(src), size(dest))
    n_not_copied = size(src) - n_copied
    dest(1:n_copied) = src(1:n_copied)

    return
  end subroutine array_copy

end module gaussian
