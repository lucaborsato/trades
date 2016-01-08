module eq_motion
  use constants,only:dp,zero,Giau
  use parameters,only:NB
  use celestial_mechanics,only:dist
  
  contains

!   ------------------------------------------------------------------ !
!   the eq. of motion in astrocentric coordinates
  subroutine eqmastro_1(m,r,drdt)
    real(dp),dimension(:),intent(in)::m,r
    real(dp),dimension(:),intent(out)::drdt
    real(dp),dimension(3)::AA,AB,BB
    real(dp),dimension(3)::ri,rj,rij
    integer::i,j,nci,ncj
    real(dp)::rimod,rjmod,rijmod

    drdt=zero

    do i=1,NB
      nci=(i-1)*6
      drdt(1+nci:3+nci)=r(4+nci:6+nci) !planets, vel.
    end do

    AA=zero
    AB=zero
    BB=zero

    do i=2,NB
      nci=(i-1)*6
      ri=r(1+nci:3+nci)
!       rimod=dist(r(1+nci:3+nci))
      rimod=dist(ri)
      AA=ri/(rimod**3)
      drdt(4+nci:6+nci)=-Giau*(m(1)+m(i))*AA
      do j=2,NB
        if(j.ne.i)then
          ncj=(j-1)*6
          rj=r(1+ncj:3+ncj)
          rij=ri-rj
          rijmod=dist(ri,rj)
          AB=rij/(rijmod**3)
          rjmod=dist(rj)
          BB=rj/(rjmod**3)
          drdt(4+nci:6+nci)=drdt(4+nci:6+nci)-Giau*m(j)*(AB+BB)
        end if
      end do
    end do

    return
  end subroutine eqmastro_1
  
  subroutine eqmastro_2(m,r,drdt)
    real(dp),dimension(:),intent(in)::m,r
    real(dp),dimension(:),intent(out)::drdt
    real(dp),dimension(NB,3)::AA
    real(dp),dimension(NB,NB,3)::AB
!     real(dp),dimension(:,:),allocatable::AA
!     real(dp),dimension(:,:,:),allocatable::AB
    real(dp),dimension(3)::ri,rj,rij
    integer::i,j,nci,ncj
    real(dp)::rimod,rijmod

    drdt=zero
!     allocate(AA(NB,3), AB(NB,NB,3))
    AA=zero
    AB=zero

    drdt(1:3)=r(4:6) !star velocity... zero
    do i=2,NB
      nci=(i-1)*6
      drdt(1+nci:3+nci)=r(4+nci:6+nci) !planets, vel.
      ri=r(1+nci:3+nci)
      rimod=dist(ri)
      AA(i,:)=ri/(rimod**3)
      drdt(4+nci:6+nci)=-Giau*(m(1)+m(i))*AA(i,:)
    end do

    
    do i=2,NB-1
      nci=(i-1)*6
      ri=r(1+nci:3+nci)
      do j=i+1,NB
        ncj=(j-1)*6
        rj=r(1+ncj:3+ncj)
        rij=ri-rj
        rijmod=dist(ri,rj)
        AB(i,j,:)=rij/(rijmod**3)
        AB(j,i,:)=-AB(i,j,:)
  !       drdt(4+nci:6+nci)=drdt(4+nci:6+nci)-Giau*m(j)*(AB+BB)
        drdt(4+nci:6+nci)=drdt(4+nci:6+nci)-Giau*m(j)*(AB(i,j,:)+AA(j,:))
        drdt(4+ncj:6+ncj)=drdt(4+ncj:6+ncj)-Giau*m(i)*(AB(j,i,:)+AA(i,:))
      end do
    end do
    
!     deallocate(AA, AB)
    
    return
  end subroutine eqmastro_2

  subroutine eqmastro(m,r,drdt)
    real(dp),dimension(:),intent(in)::m,r
    real(dp),dimension(:),intent(out)::drdt
    real(dp),dimension(3)::AA,BB
    real(dp),dimension(NB,NB,3)::AB
    real(dp),dimension(3)::ri,rj,rij
    integer::i,j,nci,ncj
    real(dp)::rimod,rjmod,rijmod

    AA=zero
    BB=zero
    AB=zero

    drdt=zero
!     drdt(1:3)=r(4:6) !star velocity... zero
    
    do i=2,NB,1
      nci=(i-1)*6
      drdt(1+nci:3+nci)=r(4+nci:6+nci)
      ri=r(1+nci:3+nci)
      rimod=dist(ri)
      AA=ri/(rimod**3)
      drdt(4+nci:6+nci)=drdt(4+nci:6+nci)-Giau*(m(1)+m(i))*AA
      do j=i+1,NB,1
        ncj=(j-1)*6
!         drdt(1+ncj:3+ncj)=r(4+ncj:6+ncj)
        rj=r(1+ncj:3+ncj)
        rij=ri-rj
        rjmod=dist(rj)
        rijmod=dist(ri,rj)
        AB(i,j,:)=rij/(rijmod**3)
!         AB(j,i,:)=-AB(i,j,:)
        BB=rj/(rjmod**3)
  !       drdt(4+nci:6+nci)=drdt(4+nci:6+nci)-Giau*m(j)*(AB+BB)
        drdt(4+nci:6+nci)=drdt(4+nci:6+nci)-Giau*m(j)*(AB(i,j,:)+BB)
!         drdt(4+ncj:6+ncj)=drdt(4+ncj:6+ncj)-Giau*m(i)*(AB(j,i,:)+AA)
        drdt(4+ncj:6+ncj)=drdt(4+ncj:6+ncj)-Giau*m(i)*(-AB(i,j,:)+AA)
      end do
    end do
    
!     deallocate(AA, AB)
    
    return
  end subroutine eqmastro

!   ------------------------------------------------------------------ !

end module eq_motion














