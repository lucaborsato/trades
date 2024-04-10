module eq_motion
    use constants, only: dp, zero, Giau, cm_au, cm_au_2
    use parameters, only: NB
    use celestial_mechanics, only: dist

contains


    subroutine eqmastro(m, r, drdt)
        real(dp), dimension(:), intent(in)::m, r
        real(dp), dimension(:), intent(out)::drdt

        real(dp), dimension(3)::AA, BB, AB
        real(dp), dimension(3)::ri, rj, rij
        integer::i, j, nci, ncj
        real(dp)::rimod, rjmod, rijmod

        drdt = zero

        do i = 2, NB, 1
            nci = (i-1)*6
            drdt(1+nci:3+nci) = r(4+nci:6+nci)
            ri = r(1+nci:3+nci)
            rimod = dist(ri)
            AA = ri/(rimod**3)
            drdt(4+nci:6+nci) = drdt(4+nci:6+nci)-Giau*(m(1)+m(i))*AA
            do j = i+1, NB, 1
                ncj = (j-1)*6
                rj = r(1+ncj:3+ncj)
                rij = ri-rj
                rjmod = dist(rj)
                rijmod = dist(ri, rj)
                AB = rij/((rijmod**2+cm_au_2)**1.5_dp)
                BB = rj/(rjmod**3)
                drdt(4+nci:6+nci) = drdt(4+nci:6+nci)-Giau*m(j)*(AB+BB)
                drdt(4+ncj:6+ncj) = drdt(4+ncj:6+ncj)-Giau*m(i)*(-AB+AA)
            end do
        end do

        return
    end subroutine eqmastro

!   ------------------------------------------------------------------ !

end module eq_motion

