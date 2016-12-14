! ************************* !
! ROTATIONS                !
! ************************* !

module rotations
  use constants,only:dp,deg2rad,zero,one
  implicit none

  contains
  
  ! rotation around axis 1 (x)
  subroutine rotmat1(angle,matrix)
    real(dp),intent(in)::angle
    real(dp),dimension(:,:),intent(out)::matrix
    real(dp)::ang,cosa,sina

    ang=angle*deg2rad
    cosa=cos(ang)
    sina=sin(ang)
    matrix(1,1)=one
    matrix(1,2)=zero
    matrix(1,3)=zero
    matrix(2,1)=zero
    matrix(2,2)=cosa
    matrix(2,3)=sina
    matrix(3,1)=zero
    matrix(3,2)=-sina
    matrix(3,3)=cosa

    return
  end subroutine rotmat1

  ! rotation around axis 2 (y)
  subroutine rotmat2(angle,matrix)
    real(dp),intent(in)::angle
    real(dp),dimension(:,:),intent(out)::matrix
    real(dp)::ang,cosa,sina

    ang=angle*deg2rad
    cosa=cos(ang)
    sina=sin(ang)
    matrix(1,1)=cosa
    matrix(1,2)=zero
    matrix(1,3)=-sina
    matrix(2,1)=zero
    matrix(2,2)=one
    matrix(2,3)=zero
    matrix(3,1)=sina
    matrix(3,2)=zero
    matrix(3,3)=cosa

    return
  end subroutine rotmat2

  ! rotation around axis 3 (z)
  subroutine rotmat3(angle,matrix)
    real(dp),intent(in)::angle
    real(dp),dimension(:,:),intent(out)::matrix
    real(dp)::ang,cosa,sina

    ang=angle*deg2rad
    cosa=cos(ang)
    sina=sin(ang)
    matrix(1,1)=cosa
    matrix(1,2)=sina
    matrix(1,3)=zero
    matrix(2,1)=-sina
    matrix(2,2)=cosa
    matrix(2,3)=zero
    matrix(3,1)=zero
    matrix(3,2)=zero
    matrix(3,3)=one

    return
  end subroutine rotmat3

  ! successive rotations around axis: 3(z) 1(x) 3(z)
  subroutine rotmat313(angle1,angle2,angle3,CBA)
    real(dp),intent(in)::angle1,angle2,angle3
    real(dp),dimension(:,:),intent(out)::CBA
    real(dp),dimension(3,3)::A,B,C,BA

    call rotmat3(angle3,A)
    call rotmat1(angle2,B)
    call rotmat3(angle1,C)
    BA=matmul(B,A)
    CBA=matmul(C,BA)

    return
  end subroutine rotmat313

  ! rotate the orbital state vectors to observational state vectors
  ! User must pass -ang1, -ang2, -ang3 to obtain the right transformation
  subroutine orb2obs(rin,ang1,ang2,ang3,rout)
    use parameters,only:NB
    real(dp),dimension(:),intent(in)::rin
    real(dp),dimension(:),intent(out)::rout
    real(dp),dimension(:),intent(in)::ang1,ang2,ang3
    real(dp),dimension(3,3)::CBA
    integer::j,ncj

    rout=zero
    do j=2,NB
      ncj=(j-1)*6
      call rotmat313(ang1(j),ang2(j),ang3(j),CBA)
      rout(1+ncj:3+ncj)=matmul(CBA,rin(1+ncj:3+ncj))
      rout(4+ncj:6+ncj)=matmul(CBA,rin(4+ncj:6+ncj))
    end do

    return
  end subroutine orb2obs


end module rotations

