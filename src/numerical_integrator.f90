module numerical_integrator
  use constants,only:dp,zero,one,two,three,TOLERANCE,TOL_dp
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
    real(dp),parameter::A21=one/5.0_dp
    real(dp),parameter::A31=three/40.0_dp,A32=9.0_dp/40.0_dp
    real(dp),parameter::A41=three/10.0_dp,A42=-9.0_dp/10.0_dp,&
        &A43=6.0_dp/5.0_dp
    real(dp),parameter::A51=-11.0_dp/54.0_dp,A52=5.0_dp/2.0_dp,&
        &A53=-70.0_dp/27.0_dp,A54=35.0_dp/27.0_dp
    real(dp),parameter::A61=1631.0_dp/55296.0_dp,A62=175.0_dp/512.0_dp,&
        &A63=575.0_dp/13824.0_dp,A64=44275.0_dp/110592.0_dp,A65=253.0_dp/4096.0_dp

    !coefficients B[s]
    !used in the final state vector (4th order):
    !y[n+1]=y[n]+sum(B[s]*k[s])
    real(dp),parameter::B1=37.0_dp/378.0_dp,B3=250.0_dp/621.0_dp,&
        &B4=125.0_dp/594.0_dp,B6=512.0_dp/1771.0_dp

    !coefficients C[s]
    !used in the intermediate state vector:
    !ks:=[t+Cs*h,y[t]+sum(A[s,s-1]*k[s-1])]
    real(dp),parameter::C2=one/5.0_dp,C3=three/10.0_dp,C4=three/5.0_dp,&
        &C6=7.0_dp/8.0_dp

    !coefficients E[s], s=1,6 for the error
    real(dp),parameter::E1=B1-2825.0_dp/27648.0_dp,&
        &E3=B3-18575.0_dp/48384.0_dp,E4=B4-13525.0_dp/55296.0_dp,&
        &E5=-277.0_dp/14336.0_dp,E6=B6-0.25_dp

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

  function get_emax(emold,error,rscal) result(emax)
    real(dp)::emax
    real(dp),intent(in)::emold
    real(dp),dimension(:),intent(in)::error,rscal
    real(dp),dimension(:),allocatable::ertemp
    real(dp)::maxtemp
    integer::i

    allocate(ertemp(size(error)))
    ertemp=zero
!     ertemp(7:NBDIM)=abs(error(7:NBDIM))/rscal(7:NBDIM)
    do i=7,NBDIM
      !if(rscal(i).eq.zero)ertemp(i)=TOL_dp
!       if(abs(rscal(i)-zero).le.TOL_dp)ertemp(i)=TOL_dp
      if(abs(rscal(i)).le.TOL_dp)then
        ertemp(i)=TOL_dp
      else
        ertemp(i)=abs(error(i))/rscal(i)
      end if
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
    real(dp),parameter::exp1=one/5.0_dp


    hh=h !uses a temporary variable
    scale_factor=one
    allocate(rscal(NBDIM))
    rscal=abs(rin)+abs(hh*drdt)
    sel: do
!       emax=0._dp
      emax=TOL_dp
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
  
  
  subroutine integrates_rk(m,rin,dt,rout)
    real(dp),dimension(:),intent(in)::m,rin
    real(dp),intent(in)::dt
    real(dp),dimension(:),intent(out)::rout
    
    real(dp)::itime,hw,hok,hnext
    real(dp),dimension(:),allocatable::r1,r2,drdt,error
    
    
    rout=zero
    allocate(r1(NBDIM),r2(NBDIM),drdt(NBDIM),error(NBDIM))
    r1=rin
    r2=rin
    
    hw=half*dt
    itime = zero
    
    
    loopint: do
      call eqmastro(m,r1,drdt)
      if(abs(itime+hw).gt.abs(dt))then

        hw=dt-itime
        call rkck_a(m,r1,drdt,hw,r2,error)
        itime=dt

      else
      
        call int_rk_a(m,r1,drdt,hw,hok,hnext,r2,error)
        itime=itime+hok
        hw=hnext
      
      end if
      if(abs(itime-dt).le.TOL_dp) exit loopint
      r1=r2
    end do loopint
    rout=r2
    deallocate(r1,r2,drdt,error)
    

    return
  end subroutine integrates_rk

end module numerical_integrator
