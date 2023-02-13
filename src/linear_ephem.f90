module linear_ephem
    use constants
    use custom_type
    use parameters
    use lin_fit, only: linfit
    implicit none

    interface set_ephem
        module procedure set_ephem_noinput, set_ephem_winput_dataT0,&
          &set_ephem_winput_dataObs
    end interface set_ephem

contains

    ! ------------------------------------------------------------------ !
    ! given the T0 data it does a linear fit to the data and it finds the
    ! ephemeris T and P: tn = Tref + Pref*n
    subroutine set_ephem_noinput()

        ! write (*, '(a)') " COMPUTING LINEAR EPHEMERIS OF: "

        call set_ephem_winput_dataObs(obsData)

        return
    end subroutine set_ephem_noinput

    subroutine set_ephem_winput_dataT0(oT0)
        type(dataT0), intent(inout)::oT0

        real(dp)::Teph, eTeph, Peph, ePeph

        if (oT0%nT0 .gt. 0) then
            Peph = zero
            ePeph = zero
            Teph = zero
            eTeph = zero
            call linfit(oT0%epo, oT0%T0, oT0%eT0, Peph, ePeph, Teph, eTeph)
            oT0%Tephem = Teph
            oT0%eTephem = eTeph
            oT0%Pephem = Peph
            oT0%ePephem = ePeph
        end if

        return
    end subroutine set_ephem_winput_dataT0

    subroutine set_ephem_winput_dataObs(oDataIn)
        type(dataObs)::oDataIn

        integer::i_body

        do i_body = 1, NB-1
            call set_ephem_winput_dataT0(oDataIn%obsT0(i_body))
        end do

        return
    end subroutine set_ephem_winput_dataObs

    subroutine set_ephem_simT0(oT0)
        type(dataT0), intent(inout)::oT0

        real(dp)::Teph, Peph

        if (oT0%nT0 .gt. 0) then
            Peph = zero
            Teph = zero
            call linfit(oT0%epo, oT0%T0, Peph, Teph)
            oT0%Tephem = Teph
            oT0%Pephem = Peph
        end if

        return
    end subroutine set_ephem_simT0

    subroutine compute_oc_one_planet(oT0)
        type(dataT0), intent(inout)::oT0

        if (oT0%nT0 .gt. 0) then
            if (.not. allocated(oT0%oc)) allocate (oT0%oc(oT0%nT0))
            oT0%oc = oT0%T0-(oT0%Tephem+oT0%Pephem*oT0%epo)
        end if

        return
    end subroutine compute_oc_one_planet

    subroutine compute_oc(oT0s)
        type(dataT0), dimension(:), intent(inout)::oT0s

        integer::npl, ipl

        npl = size(oT0s)
        do ipl = 1, npl
            call compute_oc_one_planet(oT0s(ipl))
        end do

        return
    end subroutine compute_oc

end module linear_ephem
