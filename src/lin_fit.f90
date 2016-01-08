module lin_fit
  use constants,only:dp
  implicit none

  interface linfit
    module procedure linfit_0,linfit_1
  end interface linfit

  contains

  ! ------------------------------------------------------------------ !
  ! linear fit with vector x and y and error on y (ey) => y = m*x + q
  subroutine linfit_0(x,y,ey,m,em,q,eq)
    integer,dimension(:),intent(in)::x
    real(dp),dimension(:),intent(in)::y,ey
    real(dp),intent(out)::m,em,q,eq
    integer::n
    real(dp),dimension(:),allocatable::xf,w,wyx,wy,wx,wx2
    real(dp)::Sw,Swyx,Swy,Swx,Swx2,delta,num,eq2,em2

    n=size(x)
    allocate(xf(n),w(n),wyx(n),wy(n),wx(n),wx2(n))

    xf=real(x,dp)
    w=1._dp/(ey*ey)
    Sw=sum(w)
    wyx=w*y*xf
    Swyx=sum(wyx)
    wy=w*y
    Swy=sum(wy)
    wx=w*xf
    Swx=sum(wx)
    wx2=w*(xf*xf)
    Swx2=sum(wx2)
    delta=Sw*Swx2-(Swx*Swx)
    num=Swy*Swx2-Swx*Swyx
    q=num/delta
    !write(*,*)"q => delta",delta,"num",num,"num/delta",q
    num=Sw*Swyx-Swy*Swx
    m=num/delta
    !write(*,*)"m => delta",delta,"num",num,"num/delta",m
    eq2=Swx2/delta
    eq=sqrt(eq2/(n-2))
    !write(*,*)"eq => eq2=Swx2/delta",eq2,"eq2/(n-2)",eq2/(n-2),"eq",eq
    em2=Sw/delta
    em=sqrt(em2/(n-2))
    !write(*,*)"em => em2=Sw/delta",em2,"em2/(n-2)",em2/(n-2),"em",em
    deallocate(xf,w,wyx,wy,wx,wx2)

    return
  end subroutine linfit_0

  ! linear fit with vector x and y and error on y (ey) => y = m*x + q
  subroutine linfit_1(x,y,ey,m,em,q,eq)
    real(dp),dimension(:),intent(in)::x
    real(dp),dimension(:),intent(in)::y,ey
    real(dp),intent(out)::m,em,q,eq
    integer::n
    real(dp),dimension(:),allocatable::w,wyx,wy,wx,wx2
    real(dp)::Sw,Swyx,Swy,Swx,Swx2,delta,num,eq2,em2

    n=size(x)
    allocate(w(n),wyx(n),wy(n),wx(n),wx2(n))
    w=1._dp/(ey*ey)
    Sw=sum(w)
    wyx=w*y*x
    Swyx=sum(wyx)
    wy=w*y
    Swy=sum(wy)
    wx=w*x
    Swx=sum(wx)
    wx2=w*(x*x)
    Swx2=sum(wx2)
    delta=Sw*Swx2-(Swx*Swx)
    num=Swy*Swx2-Swx*Swyx
    q=num/delta
    num=Sw*Swyx-Swy*Swx
    m=num/delta
    eq2=Swx2/delta
    eq=sqrt(eq2/(n-2))
    em2=Sw/delta
    em=sqrt(em2/(n-2))
    deallocate(w,wyx,wy,wx,wx2)

    return
  end subroutine linfit_1
  ! ------------------------------------------------------------------ !

end module lin_fit
