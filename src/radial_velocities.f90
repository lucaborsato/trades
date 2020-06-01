module radial_velocities
  use constants
  use custom_type
  use parameters
  implicit none

  interface check_RV
    module procedure check_RV_1,check_RV_2
  end interface check_RV

  
  contains

  ! ------------------------------------------------------------------ !
  ! uses the integrator to move to the obs RV jd and computed the RV (sim in m/s)
  subroutine calcRV(m,ri,drdt,stepRV,RVsim)
    use numerical_integrator,only:rkck_a
    use celestial_mechanics,only:barycenter
    real(dp),dimension(:),intent(in)::m
    real(dp),dimension(:),intent(in)::ri,drdt
    real(dp),intent(in)::stepRV
    real(dp),intent(out)::RVsim

    real(dp),dimension(:),allocatable::ro,err,rbarRV
    real(dp),dimension(:),allocatable::barRV

    allocate(barRV(6),rbarRV(NBDIM),ro(NBDIM),err(NBDIM))
    call rkck_a(m,ri,drdt,stepRV,ro,err) ! call integrator
    call barycenter(m,ro,barRV,rbarRV) ! astrocentric 2 barycentric
    RVsim=-rbarRV(6)*AU/s24h ! rv as -Zstar,bar m/s
    deallocate(barRV,rbarRV,ro,err)

    return
  end subroutine calcRV

  ! calculation of the trend of RV, order == gamma, non taken into account
  subroutine addRVtrend(t, coeff, trend)
    real(dp),dimension(:),intent(in)::t
    real(dp),dimension(:),intent(in)::coeff
    real(dp),dimension(:),intent(out)::trend

    real(dp),dimension(size(t))::w
    integer::o,order

    w = t - half*(maxval(t)-maxval(t))
    trend = zero
    order = size(coeff)
    do o=1,order
      trend = coeff(o)*w**o
    end do ! o -> order

    return
  end subroutine addRVtrend

  ! checks the RV, it uses the calcRV subroutine
!   subroutine check_RV_1(m,ri,drdt,ttemp,hok,cntRV,RV_stat,RV_sim)
subroutine check_RV_1(m,ri,drdt,ttemp,hok,simRV)
    real(dp),dimension(:),intent(in)::m,ri,drdt
    real(dp),intent(in)::ttemp
!     integer,intent(inout)::cntRV
!     real(dp),dimension(:),intent(inout)::RV_sim
!     integer,dimension(:),intent(inout)::RV_stat
    type(dataRV),intent(inout)::simRV
    real(dp),intent(inout)::hok
    
    real(dp)::stepRV,xRV
    integer::nRV,j

    nRV=obsData%obsRV%nRV
    xRV=zero
    ! RV check
    do j=1,nRV
      if(simRV%RV_stat(j).eq.0)then
        stepRV=obsData%obsRV%jd(j)-tepoch-ttemp
        if( stepRV*(stepRV-hok).le.zero )then
          call calcRV(m,ri,drdt,stepRV,xRV)
          simRV%RV(j)=xRV
          simRV%jd(j)=obsData%obsRV%jd(j)
          simRV%nRV=simRV%nRV+1
          simRV%RV_stat(j)=1
        end if
      end if
    end do

    return
  end subroutine check_RV_1

!   subroutine check_RV_2(m,ri,drdt,ttemp,hok,cntRV,tRV,RV_stat,RV_sim)
subroutine check_RV_2(m,ri,drdt,ttemp,hok,obsjd,simRV)
    real(dp),dimension(:),intent(in)::m,ri,drdt
    real(dp),intent(in)::ttemp
!     integer,intent(inout)::cntRV
!     real(dp),dimension(:),intent(in)::tRV
!     integer,dimension(:),intent(inout)::RV_stat
!     real(dp),dimension(:),intent(inout)::RV_sim
    real(dp),dimension(:),intent(in)::obsjd
    type(dataRV),intent(inout)::simRV
    real(dp),intent(inout)::hok
    real(dp)::stepRV
    integer::j,n_RV

    n_RV=size(obsjd)
    ! RV check
    do j=1,n_RV
      if(simRV%RV_stat(j).eq.0)then
        stepRV=obsjd(j)-tepoch-ttemp
        if(stepRV*(stepRV-hok).le.zero)then
          call calcRV(m,ri,drdt,stepRV,simRV%RV(j))
          simRV%jd=obsjd(j)
          simRV%nRV=simRV%nRV+1
          simRV%RV_stat(j)=1
        end if
      end if
    end do

    return
  end subroutine check_RV_2

  
end module radial_velocities
