module radial_velocities
  use constants
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
!     write(*,*)m
!     write(*,*)ri
!     write(*,*)drdt
!     write(*,*)stepRV
!     write(*,*)ro
    call barycenter(m,ro,barRV,rbarRV) ! astrocentric 2 barycentric
    RVsim=-rbarRV(6)*AU/s24h ! rv as -Zstar,bar m/s
!     write(*,*)barRV
!     write(*,*)rbarRV
!     write(*,*)RVsim
    deallocate(barRV,rbarRV,ro,err)

    return
  end subroutine calcRV

  ! checks the RV, it uses the calcRV subroutine
  subroutine check_RV_1(m,ri,drdt,ttemp,hok,cntRV,RV_stat,RV_sim)
    real(dp),dimension(:),intent(in)::m,ri,drdt
    real(dp),intent(in)::ttemp
    integer,intent(inout)::cntRV
    real(dp),dimension(:),intent(inout)::RV_sim
    integer,dimension(:),intent(inout)::RV_stat
    real(dp),intent(inout)::hok
    real(dp)::stepRV
    integer::j

    ! RV check
    do j=1,nRV
      if(RV_stat(j).eq.0)then
        stepRV=jdRV(j)-tepoch-ttemp
        if( stepRV*(stepRV-hok).le.zero )then
          call calcRV(m,ri,drdt,stepRV,RV_sim(j))
          cntRV=cntRV+1
          RV_stat(j)=1
        end if
      end if
    end do

    return
  end subroutine check_RV_1

  subroutine check_RV_2(m,ri,drdt,ttemp,hok,cntRV,tRV,RV_stat,RV_sim)
    real(dp),dimension(:),intent(in)::m,ri,drdt
    real(dp),intent(in)::ttemp
    integer,intent(inout)::cntRV
    real(dp),dimension(:),intent(in)::tRV
    integer,dimension(:),intent(inout)::RV_stat
    real(dp),dimension(:),intent(inout)::RV_sim
    real(dp),intent(inout)::hok
    real(dp)::stepRV
    integer::j,n_RV

    n_RV=size(tRV)
    ! RV check
    do j=1,n_RV
      if(RV_stat(j).eq.0)then
        stepRV=tRV(j)-tepoch-ttemp
        if(stepRV*(stepRV-hok).le.zero)then
          call calcRV(m,ri,drdt,stepRV,RV_sim(j))
          cntRV=cntRV+1
          RV_stat(j)=1
        end if
      end if
    end do

    return
  end subroutine check_RV_2

  
end module radial_velocities
