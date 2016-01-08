module numerical_integrator
  use constants,only:dp,zero,TOLERANCE
  use parameters
  use eq_motion
  implicit none

  contains

  ! ------------------------------------------------------------------ !
  ! RUNGE-KUTTA-CASH-KARP by Numerical Recepis
  ! rkck integrator that calls only eq. of motion in astrocentric coordinates
  !
  subroutine rkck_a(m,rin,drdt,h,rout,error)
    real(dp),dimension(:),intent(in)::m,rin,drdt
    real(dp),intent(in)::h
    real(dp),dimension(:),intent(out)::rout,error
    real(dp),dimension(:),allocatable::rtemp

    !coefficients k(6*NB)[s], s=1,6
    real(dp),dimension(:),allocatable::k2,k3,k4,k5,k6

    !coefficients A[s,s-1]
    !used in the intermediate state vector:
    !ks:=[t+Cs*h,y[t]+sum(A[s,s-1]*k[s-1])]
    real(dp),parameter::A21=1._dp/5._dp
    real(dp),parameter::A31=3._dp/40._dp,A32=9._dp/40._dp
    real(dp),parameter::A41=3._dp/10._dp,A42=-9._dp/10._dp,&
        &A43=6._dp/5._dp
    real(dp),parameter::A51=-11._dp/54._dp,A52=5._dp/2._dp,&
        &A53=-70._dp/27._dp,A54=35._dp/27._dp
    real(dp),parameter::A61=1631._dp/55296._dp,A62=175._dp/512._dp,&
        &A63=575._dp/13824._dp,A64=44275._dp/110592._dp,A65=253._dp/4096._dp

    !coefficients B[s]
    !used in the final state vector (4th order):
    !y[n+1]=y[n]+sum(B[s]*k[s])
    real(dp),parameter::B1=37._dp/378._dp,B3=250._dp/621._dp,&
        &B4=125._dp/594._dp,B6=512._dp/1771._dp

    !coefficients C[s]
    !used in the intermediate state vector:
    !ks:=[t+Cs*h,y[t]+sum(A[s,s-1]*k[s-1])]
    real(dp),parameter::C2=1._dp/5._dp,C3=3._dp/10._dp,C4=3._dp/5._dp,&
        &C6=7._dp/8._dp

    !coefficients E[s], s=1,6 for the error
    real(dp),parameter::E1=B1-2825._dp/27648._dp,&
        &E3=B3-18575._dp/48384._dp,E4=B4-13525._dp/55296._dp,&
        &E5=-277._dp/14336._dp,E6=B6-0.25_dp

    allocate(rtemp(NBDIM),k2(NBDIM),k3(NBDIM),k4(NBDIM),k5(NBDIM),k6(NBDIM))
    !Numerical Recipes algorithm
    rtemp=rin+A21*h*drdt
    call eqmastro(m,rtemp,k2)
    rtemp=rin+h*(A31*drdt+A32*k2)
    call eqmastro(m,rtemp,k3)
    rtemp=rin+h*(A41*drdt+A42*k2+A43*k3)
    call eqmastro(m,rtemp,k4)
    rtemp=rin+h*(A51*drdt+A52*k2+A53*k3+A54*k4)
    call eqmastro(m,rtemp,k5)
    rtemp=rin+h*(A61*drdt+A62*k2+A63*k3+A64*k4+A65*k5)
    call eqmastro(m,rtemp,k6)
    rout=rin+h*(B1*drdt+B3*k3+B4*k4+B6*k6)
    error=h*(E1*drdt+E3*k3+E4*k4+E5*k5+E6*k6)
    deallocate(rtemp,k2,k3,k4,k5,k6)

    return
  end subroutine rkck_a

  function get_emax(emold,err,rscal) result(emax)
    real(dp)::emax
    real(dp),intent(in)::emold
    real(dp),dimension(:),intent(in)::err,rscal
    real(dp),dimension(:),allocatable::ertemp
    real(dp)::maxtemp
    integer::i

    allocate(ertemp(size(err)))
    ertemp=zero
    ertemp(7:NBDIM)=abs(err(7:NBDIM))/rscal(7:NBDIM)
    do i=7,NBDIM
      !if(rscal(i).eq.zero)ertemp(i)=TOLERANCE
      if(abs(rscal(i)-zero).le.TOLERANCE)ertemp(i)=TOLERANCE
    end do
    maxtemp=maxval(ertemp(7:NBDIM))
    emax=max(emold,maxtemp)
    deallocate(ertemp)

    return
  end function get_emax

  ! it calls the integrator and select the right step for the integration
  subroutine int_rk_a(m,rin,drdt,h,hok,hnext,rout,err)
    real(dp),dimension(:),intent(in)::m,rin,drdt
    real(dp),intent(in)::h
    real(dp),intent(out)::hok,hnext
    real(dp),dimension(:),intent(out)::rout,err
    real(dp),dimension(:),allocatable::rscal
    real(dp)::emax,hh,scale_factor
    ! safety factor = sfac ;
    real(dp),parameter::sfac=0.9_dp
    real(dp),parameter::exp1=1._dp/5._dp


    hh=h !uses a temporary variable
    scale_factor=1._dp
    allocate(rscal(NBDIM))
    rscal=abs(rin)+abs(hh*drdt)
    sel: do
!       emax=0._dp
      emax=TOLERANCE
      call rkck_a(m,rin,drdt,hh,rout,err)
      emax=get_emax(emax,err,rscal)
      scale_factor=sfac*((tol_int/emax)**(exp1))
      hok=hh
      hh=hh*scale_factor
      if(tol_int.ge.emax)exit sel
    end do sel
    hnext=hh
    deallocate(rscal)

    return
  end subroutine int_rk_a
  ! ------------------------------------------------------------------ !

end module numerical_integrator
