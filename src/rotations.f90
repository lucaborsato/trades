! ************************* !
! ROTATIONS                !
! ************************* !

module rotations
    use constants, only: dp, deg2rad, zero, one, TOLERANCE, TOL_dp
    implicit none

contains

    subroutine set_cossin(cosa,sina, tol)
        real(dp),intent(inout)::cosa,sina
        real(dp),intent(in)::tol

        if (abs(cosa) .le. tol) then
            cosa = zero
            sina = sign(one, sina)
        else if (abs(sina) .le. tol) then
            sina = zero
            cosa = sign(one, cosa)
        end if

        return
    end subroutine set_cossin

    ! rotation around axis 1 (x)
    subroutine rotmat1(angle, matrix)
        real(dp), intent(in)::angle
        real(dp), dimension(:, :), intent(out)::matrix
        real(dp)::ang, cosa, sina

        ang = angle*deg2rad
        cosa = cos(ang)
        sina = sin(ang)
        ! if (abs(cosa) .le. TOLERANCE) then
        !     cosa = zero
        !     sina = sign(one, sina)
        ! else if (abs(sina) .le. TOLERANCE) then
        !     sina = zero
        !     cosa = sign(one, cosa)
        ! end if
        call set_cossin(cosa,sina,TOLERANCE)

        !        (0 0     1   ) ?? (1 0 0)
        !mat(3x)=(0 cosa -sina)
        !        (0 sina  cosa)

        matrix(1, 1) = one
        matrix(1, 2) = zero
        matrix(1, 3) = zero
        matrix(2, 1) = zero
        matrix(2, 2) = cosa
!     matrix(2,3)=sina ! OLD
        matrix(2, 3) = -sina
        matrix(3, 1) = zero
!     matrix(3,2)=-sina ! OLD
        matrix(3, 2) = sina
        matrix(3, 3) = cosa

        return
    end subroutine rotmat1

    subroutine rotation_axis1(xyz_in, alpha_rad, xyz_out)
        ! Input
        real(dp),dimension(:),intent(in)::xyz_in
        real(dp),intent(in)::alpha_rad
        ! Output
        real(dp),dimension(:),intent(out)::xyz_out
        ! Local
        real(dp)::cosa,sina

        cosa = cos(alpha_rad)
        sina = sin(alpha_rad)
        ! if (abs(cosa) .le. TOLERANCE) then
        !     cosa = zero
        !     sina = sign(one, sina)
        ! else if (abs(sina) .le. TOLERANCE) then
        !     sina = zero
        !     cosa = sign(one, cosa)
        ! end if
        call set_cossin(cosa,sina,TOLERANCE)

        xyz_out(1) = xyz_in(1)
        xyz_out(2) = xyz_in(2)*cosa - xyz_in(3)*sina
        xyz_out(3) = xyz_in(2)*sina + xyz_in(3)*cosa

        return
    end subroutine rotation_axis1

!   ! rotation around axis 2 (y)
!   subroutine rotmat2(angle,matrix)
!     real(dp),intent(in)::angle
!     real(dp),dimension(:,:),intent(out)::matrix
!     real(dp)::ang,cosa,sina
!
!     ang=angle*deg2rad
!     cosa=cos(ang)
!     sina=sin(ang)
!     if(abs(cosa).le.TOLERANCE)then
!       cosa=zero
!       sina=sign(one,sina)
!     else if(abs(sina).le.TOLERANCE)then
!       sina=zero
!       cosa=sign(one,cosa)
!     end if
!     matrix(1,1)=cosa
!     matrix(1,2)=zero
!     matrix(1,3)=-sina
!     matrix(2,1)=zero
!     matrix(2,2)=one
!     matrix(2,3)=zero
!     matrix(3,1)=sina
!     matrix(3,2)=zero
!     matrix(3,3)=cosa
!
!     return
!   end subroutine rotmat2

    ! rotation around axis 3 (z)
    subroutine rotmat3(angle, matrix)
        real(dp), intent(in)::angle
        real(dp), dimension(:, :), intent(out)::matrix
        real(dp)::ang, cosa, sina

        ang = angle*deg2rad
        cosa = cos(ang)
        sina = sin(ang)
        ! if (abs(cosa) .le. TOLERANCE) then
        !     cosa = zero
        !     sina = sign(one, sina)
        ! else if (abs(sina) .le. TOLERANCE) then
        !     sina = zero
        !     cosa = sign(one, cosa)
        ! end if
        call set_cossin(cosa,sina,TOLERANCE)

        !        (cosa -sina 0)
        !mat(3x)=(sina  cosa 0)
        !        (0     0    1)

        matrix(1, 1) = cosa
!     matrix(1,2)=sina ! OLD
        matrix(1, 2) = -sina
        matrix(1, 3) = zero
!     matrix(2,1)=-sina ! OLD
        matrix(2, 1) = sina
        matrix(2, 2) = cosa
        matrix(2, 3) = zero
        matrix(3, 1) = zero
        matrix(3, 2) = zero
        matrix(3, 3) = one

        return
    end subroutine rotmat3

    subroutine rotation_axis3(xyz_in, alpha_rad, xyz_out)
        ! Input
        real(dp),dimension(:),intent(in)::xyz_in
        real(dp),intent(in)::alpha_rad
        ! Output
        real(dp),dimension(:),intent(out)::xyz_out
        ! Local
        real(dp)::cosa,sina

        cosa = cos(alpha_rad)
        sina = sin(alpha_rad)
        ! if (abs(cosa) .le. TOLERANCE) then
        !     cosa = zero
        !     sina = sign(one, sina)
        ! else if (abs(sina) .le. TOLERANCE) then
        !     sina = zero
        !     cosa = sign(one, cosa)
        ! end if
        call set_cossin(cosa,sina,TOLERANCE)

        xyz_out(1) = xyz_in(1)*cosa - xyz_in(2)*sina
        xyz_out(2) = xyz_in(1)*sina + xyz_in(2)*cosa
        xyz_out(3) = xyz_in(3)

        return
    end subroutine rotation_axis3

    ! successive rotations around axis: 3(z) 1(x) 3(z)
    subroutine rotmat313(angle1, angle2, angle3, CBA)
        real(dp), intent(in)::angle1, angle2, angle3
        real(dp), dimension(:, :), intent(out)::CBA
        real(dp), dimension(3, 3)::A, B, C, BA

        call rotmat3(angle3, A)
        call rotmat1(angle2, B)
        call rotmat3(angle1, C)
        BA = matmul(B, A)
        CBA = matmul(C, BA)

        return
    end subroutine rotmat313

    ! rotate the orbital state vectors to observational state vectors
    ! User must pass -ang1, -ang2, -ang3 to obtain the right transformation..TESTING
    subroutine orb2obs(rin, ang1, ang2, ang3, rout)
        use parameters, only: NB
        real(dp), dimension(:), intent(in)::rin
        real(dp), dimension(:), intent(out)::rout
        real(dp), dimension(:), intent(in)::ang1, ang2, ang3
        real(dp), dimension(3, 3)::CBA
        integer::j, ncj

        rout = zero
        do j = 2, NB
            ncj = (j-1)*6
            call rotmat313(ang1(j), ang2(j), ang3(j), CBA)
            rout(1+ncj:3+ncj) = matmul(CBA, rin(1+ncj:3+ncj))
            rout(4+ncj:6+ncj) = matmul(CBA, rin(4+ncj:6+ncj))
        end do

        return
    end subroutine orb2obs

    subroutine rotate_vector(xyz_in, argp, inc, longn, xyz_out)
        ! Input
        real(dp),dimension(:),intent(in)::xyz_in
        real(dp),intent(in)::argp,inc,longn
        ! Output
        real(dp),dimension(:),intent(out)::xyz_out
        ! Local
        real(dp)::argp_rad, inc_rad, longn_rad
        real(dp),dimension(3)::xyz_a1, xyz_a2

        argp_rad = argp*deg2rad
        inc_rad = inc*deg2rad
        longn_rad = longn*deg2rad

        call rotation_axis3(xyz_in, argp_rad, xyz_a1)
        call rotation_axis1(xyz_a1, inc_rad, xyz_a2)
        call rotation_axis3(xyz_a2, longn_rad, xyz_out)

        return
    end subroutine rotate_vector

end module rotations

