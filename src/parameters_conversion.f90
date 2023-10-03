module parameters_conversion
    use constants
    use parameters
    use convert_type
    use statistics

contains

    subroutine set_all_parameter_names(kel_id)
        character(6), dimension(3:10), intent(in)::kel_id

        integer::i_par, body, cnt_kel
        allocate (all_names_list(npar)) ! all_names_list in 'parameters' module

        all_names_list(1) = trim(adjustl(kel_id(3)))//'1' ! mass of the star, id 1
        all_names_list(2) = trim(adjustl(kel_id(4)))//'1' ! radius of the star, id 1
        all_names_str = trim(adjustl(all_names_list(1)))//' '//trim(adjustl(all_names_list(2)))
        cnt_kel = 3
        do i_par = 3, npar
            body = int(int(i_par-3)/8)+2
            all_names_list(i_par) = trim(adjustl(kel_id(cnt_kel)))//trim(adjustl(string(body)))
            all_names_str = trim(adjustl(all_names_str))//' '//trim(adjustl(all_names_list(i_par)))
            cnt_kel = cnt_kel+1
            if (cnt_kel .eq. 11) cnt_kel = 3
        end do
!     write(*,'(a)')trim(adjustl(all_names_str))

        return
    end subroutine set_all_parameter_names

    subroutine idpar()
        integer::pos, cntid, j, body
        character(6), dimension(3:10)::elid
        data elid/"m", "R", "P", "e", "w", "mA", "i", "lN"/
        integer::nkep
        integer::order, iset

        ! write(*,*)"DEBUG: idpar()"
        ! write(*,*)

        if (.not. allocated(id)) allocate (id(nfit), idall(nfit), parid(nfit))
    
        pos = 0
        ! if Mstar to fit:
        j = 1
        if (tofit(j) .eq. 1)then
            pos = pos+1
            cntid = 3
            ! write(*,*)
            ! write(*,*)"DEBUG: start  j = ",j,"    pos = ",pos," cntid = ",cntid, " tofit(j)   = ",tofit(j)
            id(pos) = cntid
            idall(pos) = j
            parid(pos) = trim(adjustl(elid(cntid)))//trim(adjustl(string(1)))
            ! write(*,*)"DEBUG: update j = ",j,"    pos = ",pos,    " cntid = ",cntid, " tofit(j)   = ",tofit(j)
            ! write(*,*)"DEBUG: update j = ",j,"id(pos) = ",id(pos)," body  = ",1,  " parid(pos) = ",trim(parid(pos))
        end if
        ! if Rstar to fit:
        j = 2
        if (tofit(j) .eq. 1)then
            pos = pos+1
            cntid = 4
            ! write(*,*)
            ! write(*,*)"DEBUG: start  j = ",j,"    pos = ",pos," cntid = ",cntid, " tofit(j)   = ",tofit(j)
            id(pos) = cntid
            idall(pos) = j
            parid(pos) = trim(adjustl(elid(cntid)))//trim(adjustl(string(1)))
            ! write(*,*)"DEBUG: update j = ",j,"    pos = ",pos,    " cntid = ",cntid, " tofit(j)   = ",tofit(j)
            ! write(*,*)"DEBUG: update j = ",j,"id(pos) = ",id(pos)," body  = ",body,  " parid(pos) = ",trim(parid(pos))
        end if

        cntid = 2
        do j = 3, npar
            ! write(*,*)
            ! write(*,*)"DEBUG: start  j = ",j,"    pos = ",pos," cntid = ",cntid, " tofit(j)   = ",tofit(j)
            cntid = cntid+1

            if (tofit(j) .eq. 1) then
                pos = pos+1
                id(pos) = cntid
                idall(pos) = j
                body = int(int(j-3)/8)+2
                parid(pos) = trim(adjustl(elid(cntid)))//trim(adjustl(string(body)))
                parid(pos) = trim(adjustl(parid(pos)))
                ! write(*,*)"DEBUG: update j = ",j,"    pos = ",pos,    " cntid = ",cntid, " tofit(j)   = ",tofit(j)
                ! write(*,*)"DEBUG: update j = ",j,"id(pos) = ",id(pos)," body  = ",body,  " parid(pos) = ",trim(parid(pos))
            end if
            if (cntid .eq. 10) cntid = 2
        end do

        call set_all_parameter_names(elid)

        ! write(*,*)"DEBUG: idpar()"
        ! write(*,*)
        ! write(*,*)"DEBUG: all_names_list ",all_names_list

        ! nRVset = obsData%obsRV%nRVset ! global in parameters.f90
        nkep = nfit-rv_trend_order-2*nRVset
        nkel = nkep ! nkel global variable
        ! write(*,*)"nkep",nkep
        do j = 1, nkep

            body = int(int(idall(j)-3)/8)+2
            ! write(*,*)
            ! write(*,*)"DEBUG: 0) j = ",j," parid(j) = ",parid(j)


            if (id(j) .eq. 3) then ! Mp
                if (trim(parid(j)) .eq. "m1")then
                    parid(j) = "m1"
                else
                    parid(j) = "m"//trim(adjustl(string(body)))//"Ms"
                end if

            else if (id(j) .eq. 4) then ! Rp

                if (trim(parid(j)) .eq. "R1")then
                    parid(j) = "R1"
                else
                    parid(j) = "R"//trim(adjustl(string(body)))//"Rs"
                end if

            else if (id(j) .eq. 6) then
                if (j .lt. nkep) then
                    if (id(j+1) .eq. 7) then !e,w
                        parid(j) = "secosw"//trim(adjustl(string(body)))
                        parid(j+1) = "sesinw"//trim(adjustl(string(body)))
                    end if
                end if

            else if (id(j) .eq. 8) then ! mA
                parid(j) = "lambda"//trim(adjustl(string(body)))

            end if
            ! write(*,*)"DEBUG: 1) j = ",j," parid(j) = ",parid(j)

        end do ! end do j=1,nkep

        j = nkep

        if (nRVset .gt. 0) then
            ! jitter
            do iset = 1, nRVset
                j = j+1
                parid(j) = trim(adjustl("l2j_"//trim(adjustl(string(iset)))))
            end do
            ! gamma offset
            do iset = 1, nRVset
                j = j+1
                parid(j) = trim(adjustl("gamma_"//trim(adjustl(string(iset)))))
            end do
        end if

        if (rv_trend_order .gt. 0) then
            ! j=nkep # it has the last value from previous loop
            do order = 1, rv_trend_order
                j = j+1
                parid(j) = trim(adjustl("c_"//trim(adjustl(string(order)))))
            end do ! end j->rv_trend_order
        end if

        ! set keplerian elements id, for priors
        if (.not. allocated(mass_id)) allocate (mass_id(NB), radius_id(NB), period_id(NB), sma_id(NB))
        if (.not. allocated(ecc_id)) allocate (ecc_id(NB), argp_id(NB), meana_id(NB), inc_id(NB), longn_id(NB))
        do j = 1, NB
            mass_id(j) = "m"//trim(adjustl(string(j)))
            radius_id(j) = "R"//trim(adjustl(string(j)))
            period_id(j) = "P"//trim(adjustl(string(j)))
            sma_id(j) = "a"//trim(adjustl(string(j)))
            ecc_id(j) = "e"//trim(adjustl(string(j)))
            argp_id(j) = "w"//trim(adjustl(string(j)))
            meana_id(j) = "mA"//trim(adjustl(string(j)))
            inc_id(j) = "i"//trim(adjustl(string(j)))
            longn_id(j) = "lN"//trim(adjustl(string(j)))
        end do

        return
    end subroutine idpar

    ! given the id of the parameters to fit it creates a proper string
    subroutine set_parid_list()
        integer::i

        paridlist = ""
        sig_paridlist = ""
        do i = 1, nfit
            paridlist = trim(paridlist)//" "//trim(adjustl(parid(i)))
            sig_paridlist = trim(sig_paridlist)//" sig_"//trim(adjustl(parid(i)))
        end do

        return
    end subroutine set_parid_list

    ! it sets the boundaries
    subroutine set_minmax()
        integer::ifit, ipar !,nRVset
        logical, dimension(:), allocatable::done

        real(dp)::se
        real(dp)::vmin, vmax !, dRV
        ! real(dp)::cimin, cimax

        allocate (done(nfit))
        done = .false.
        if (.not. allocated(minpar)) allocate (minpar(nfit), maxpar(nfit))
        ifit = 0
        do ipar = 1, npar
            if (tofit(ipar) .eq. 1) then
                ifit = ifit+1

                if (.not. done(ifit)) then
                    ! fit m(3)==> mp/Ms
                    if (id(ifit) .eq. 3) then
                        minpar(ifit) = par_min(ipar)/(MR_star(1, 1)+MR_star(1,2))
                        maxpar(ifit) = par_max(ipar)/(MR_star(1, 1)-MR_star(1,2))
                        done(ifit) = .true.
                    
                    ! fit r(4)==> Rp/Rs
                    else if(id(ifit) .eq. 4) then
                        minpar(ifit) = par_min(ipar)/(MR_star(2, 1)+MR_star(2,2))
                        maxpar(ifit) = par_max(ipar)/(MR_star(2, 1)-MR_star(2,2))
                        done(ifit) = .true.

                        ! fit e(6) & w(7)==>(sqrt(e)cosw,sqrt(e)sinw) [-sqrt(e),+sqrt(e)],[-sqrt(e),+sqrt(e)]
                    else if (id(ifit) .eq. 6) then
                        if (ifit .lt. nfit) then
                            if (id(ifit+1) .eq. 7) then
                                se = sqrt(par_max(ipar))
                                minpar(ifit) = -se
                                minpar(ifit+1) = -se
                                maxpar(ifit) = se
                                maxpar(ifit+1) = se
                                done(ifit) = .true.
                                done(ifit+1) = .true.
                            end if
                        end if

                        ! fit mA(8) ==> lambda[0,360]
                    else if (id(ifit) .eq. 8) then
                        minpar(ifit) = zero !-circ
                        maxpar(ifit) = circ !720.0_dp
                        done(ifit) = .true.

                        ! 2022-07-29 fit inc and longN as they are (in physical unit)
                        ! fit i(9) & lN(10)==>(cosicoslN,cosisinlN) [min(cosi(),max(cosi)],[min(cosi(),max(cosi)]
                        ! else if (id(ifit) .eq. 9) then

                        ! if (ifit .lt. nfit) then
                        !     if (id(ifit+1) .eq. 10) then
                        !         ! ifit & ipar --> inc --> cos(inc)
                        !         cimin = minval(cos([par_min(ipar), par_max(ipar)]))
                        !         cimax = maxval(cos([par_min(ipar), par_max(ipar)]))

                        !         minpar(ifit)   = cimin  ! cosicoslN min
                        !         minpar(ifit+1) = cimin  ! cosisinlN min
                        !         maxpar(ifit)   = cimax  ! cosicoslN max
                        !         maxpar(ifit+1) = cimax  ! cosisinlN max

                        !         done(ifit) = .true.
                        !         done(ifit+1) = .true.
                        !     else
                        !         minpar(ifit) = par_min(ipar)
                        !         maxpar(ifit) = par_max(ipar)
                        !         done(ifit) = .true.
                        !     end if
                        ! end if

                    end if
                end if ! .not.done(ifit)

                if (.not. done(ifit)) then
                    ! default setting
                    minpar(ifit) = par_min(ipar)
                    maxpar(ifit) = par_max(ipar)
                    done(ifit) = .true.
                end if

            end if ! tofit(j)
        end do ! ipar

        ! check if nRVset > 0 and set minmax jitter to log2(1.0e-15, 100.0 m/s)
        ! nRVset=obsData%obsRV%nRVset
        if (nRVset .gt. 0) then
            ! jitter
            do ipar = 1, nRVset
                ifit = ifit+1
                minpar(ifit) = log(1.0e-15_dp)/log(two)
                maxpar(ifit) = log(100.0_dp)/log(two)
                done(ifit) = .true.
            end do ! ipar -> nRVset
            ! gamma offset
            do ipar = 1, nRVset
                ifit = ifit+1
                vmax = maxval(pack(obsData%obsRV%RV, obsData%obsRV%RVsetID == ipar))
                vmin = minval(pack(obsData%obsRV%RV, obsData%obsRV%RVsetID == ipar))
                ! dRV = vmax-vmin
                ! -+ 5000 m/s
                minpar(ifit) = vmin-5000.0_dp ! *dRV
                maxpar(ifit) = vmax+5000.0_dp ! *dRV
                done(ifit) = .true.
            end do ! ipar -> nRVset
        end if

        ! check if rv_trend_order and set minmax to (-1,1)
        if (rv_trend_order .gt. 0) then
            do ipar = 1, rv_trend_order
                ifit = ifit+1
                minpar(ifit) = -one
                maxpar(ifit) = one
                done(ifit) = .true.
            end do ! ipar -> rv_trend_order
        end if

        return
    end subroutine set_minmax
    ! ------------------------------------------------------------------ !

    ! ------------------------------------------------------------------ !
    ! defines limits of the integration as distance from the star
    subroutine sma_boundaries(radius, sma, smamin, smamax)
        real(dp), dimension(:), intent(in)::radius, sma
        real(dp), intent(out)::smamin, smamax

        smamin = radius(1)*RsunAU !smamin =  Rstar in AU
        smamax = 10.0_dp*maxval(sma) !smamax = 5 times the larger semi-major axis

        return
    end subroutine sma_boundaries

    ! puts the parameters in a big vector with dimension npar = 2+(NB-1)*8 (in tofit)
    subroutine set_all_param(mass, radius, period, ecc, argp, meanA, inc, longN, allpar)
        real(dp), dimension(:), intent(in)::mass, radius
        real(dp), dimension(:), intent(in)::period, ecc, argp, meanA, inc, longN
        real(dp), dimension(:), allocatable, intent(out)::allpar
        integer::j, j1

        if (.not. allocated(allpar)) allocate (allpar(npar))
        allpar(1) = mass(1)
        allpar(2) = radius(1)
        do j = 2, NB
            j1 = (j-2)*8
            allpar(3+j1) = mass(j)
            allpar(4+j1) = radius(j)
            allpar(5+j1) = period(j)
            allpar(6+j1) = ecc(j)
            allpar(7+j1) = argp(j)
            allpar(8+j1) = meanA(j)
            allpar(9+j1) = inc(j)
            allpar(10+j1) = longN(j)
        end do

        return
    end subroutine set_all_param

    ! set the par vector with parameters to be fitted
    subroutine init_param(allpar, par)
        real(dp), dimension(:), intent(in)::allpar
        real(dp), dimension(:), allocatable, intent(out)::par
        integer::j, ifit !,nRVset

        integer::nrvtemp
        real(dp), dimension(:), allocatable::rvtemp
        logical, dimension(:), allocatable::done

        allocate (done(nfit))
        done = .false.
        if (.not. allocated(par)) allocate (par(nfit))
        par = zero
        ifit = 0
        do j = 1, npar
            if (tofit(j) .eq. 1) then
                ifit = ifit+1

                if (.not. done(ifit)) then

                    ! m(3) ==> mp/Ms
                    if (id(ifit) .eq. 3) then
                        ! par(ifit) = allpar(j)/MR_star(1, 1)
                        par(ifit) = allpar(j)/allpar(1)
                        done(ifit) = .true.
                    
                    ! R(4) ==> Rp/Ms
                    else if (id(ifit) .eq. 4) then
                            ! par(ifit) = allpar(j)/MR_star(2, 2)
                            par(ifit) = allpar(j)/allpar(2)
                            done(ifit) = .true.

                    ! e(6),w(7)==>(sqrt(e)cosw,sqrt(e)sinw)
                    else if (id(ifit) .eq. 6) then
                        par(ifit) = allpar(j)
                        done(ifit) = .true.
                        if (ifit .lt. nfit) then
                            if (id(ifit+1) .eq. 7) then
                                par(ifit) = sqrt(allpar(j))*cos(allpar(j+1)*deg2rad)
                                par(ifit+1) = sqrt(allpar(j))*sin(allpar(j+1)*deg2rad)
                                done(ifit) = .true.
                                done(ifit+1) = .true.
                            end if
                        end if

                    ! mA(8) ==> lambda=mA(j)+w(j-1)+lN(j+2)
                    else if (id(ifit) .eq. 8) then
                        par(ifit) = mod(mod(allpar(j)+allpar(j-1)+allpar(j+2), circ)+circ, circ)
                        done(ifit) = .true.

                    else
                        par(ifit) = allpar(j)
                        done(ifit) = .true.
                    end if

                end if

                ! if (.not. done(ifit)) then
                !     par(ifit) = allpar(j)
                !     done(ifit) = .true.
                ! end if

            end if ! tofit(j)
        end do ! j

        ! check if nRVset > 0
        ! nRVset = obsData%obsRV%nRVset
        if (nRVset .gt. 0) then
            ! RV jitter
            do j = 1, nRVset
                ifit = ifit+1
                ! par(ifit) = log(1.0e-6_dp)/log(two)
                par(ifit) = log(half)/log(two)
                done(ifit) = .true.
            end do ! j -> nRVset
            ! RV gamma offset
            do j = 1, nRVset
                ifit = ifit+1
                rvtemp = pack(obsData%obsRV%RV, obsData%obsRV%RVsetID == j)
                nrvtemp = obsData%obsRV%nRVsingle(j)
                ! par(ifit) = sum(rvtemp)/nrvtemp
                par(ifit) = median_dp(rvtemp)
                done(ifit) = .true.
                deallocate (rvtemp)
            end do ! j -> nRVset
        end if

        ! Probably not needed because par = zero initially
        ! check if rv_trend_order > 0
        if (rv_trend_order .gt. 0) then
            do j = 1, rv_trend_order
                ifit = ifit+1
                par(ifit) = zero
                done(ifit) = .true.
            end do ! j -> rv_trend_order
        end if

        deallocate (done)

!     call set_minmax()

        return
    end subroutine init_param

    ! from keplerian orbital elements to the parameters needed by L-M
    ! calls some previous subroutines
    subroutine set_par(mass, radius, period, sma, ecc, argp, meanA, inc, longN, allpar, par)
        real(dp), dimension(:), intent(in)::mass, radius
        real(dp), dimension(:), intent(in)::period, sma, ecc, argp, meanA, inc, longN
        real(dp), dimension(:), allocatable, intent(out)::allpar, par

        call sma_boundaries(radius, sma, amin, amax) ! IT SETS SEMI-MAJOR AXIS BOUNDS
        call set_all_param(mass, radius, period, ecc, argp, meanA, inc, longN, allpar) ! IT DEFINES THE VECTOR ALLPAR FROM ORBITAL PARAMETERS
        call init_param(allpar, par) ! IT DEFINES THE PARAMETERS PAR TO BE FITTED

        return
    end subroutine set_par
    ! ------------------------------------------------------------------ !

    ! fix the system_parameters in case a parameter has been read with a value not in [par_min, par_max]

    subroutine fix_all_parameters(all_parameters)
        real(dp), dimension(:)::all_parameters
        integer::i

        do i = 1, npar
            if (tofit(i) .eq. 1) then
                if ((all_parameters(i) .lt. par_min(i)) .or.&
                    &(all_parameters(i) .gt. par_max(i))) then
                    all_parameters(i) = half*(par_min(i)+par_max(i))
                end if
            end if
        end do

        return
    end subroutine fix_all_parameters

    subroutine par2kel_fit(all_parameters, fit_parameters,&
        &mass, radius, period, sma, ecc, argp, meanA, inc, longN,&
        &checkpar)
        use celestial_mechanics, only: period_to_sma !get_semax
        real(dp), dimension(:), intent(in)::all_parameters, fit_parameters
        real(dp), dimension(:), intent(out)::mass, radius
        real(dp), dimension(:), intent(out)::period, sma, ecc, argp, meanA, inc, longN

        real(dp), dimension(:), allocatable::atemp
        real(dp), dimension(8)::temp_kel ! ==> mass, radius, period, ecc, argp, meanA, inc, longN
        integer::i_par, i_body
        logical, intent(out)::checkpar
        real(dp)::temp2

!     checkpar=.true.
        mass = zero
        radius = zero
        period = zero
        sma = zero
        ecc = zero
        argp = zero
        meanA = zero
        inc(1) = zero
        inc(2:NB) = 90.0_dp
        longN = zero
        allocate (atemp(npar))
        atemp = all_parameters

        ! nkel = nfit - rv_trend_order - nRVset ! it is global now
        if (nkel .gt. 0) then
            ! update atemp with fit_parameters
            do i_par = 1, nkel
                atemp(idall(i_par)) = fit_parameters(i_par)
            end do
        end if

        mass(1) = atemp(1)
        radius(1) = atemp(2)

        ! i_body = 2
        do i_body = 2, NB
            i_par = 3+((i_body-2)*8) ! first parameter id == 3 <-> Mass body 2, ..., 10 <-> lN body 2, 11 <-> Mass body 3, ...

            temp_kel = atemp(i_par:i_par+7)

            ! mp/Ms to mp
            if (tofit(i_par) .eq. 1) then
                ! temp_kel(1) = atemp(i_par)*MR_star(1, 1)
                temp_kel(1) = atemp(i_par)*mass(1)
            end if

            ! Rp/Ms to Rp
            if (tofit(i_par+1) .eq. 1) then
                temp_kel(2) = atemp(i_par+1)*radius(1)
            end if

            ! sqrt(e)cosw,sqrt(e)sinw to e,w
            if (tofit(i_par+3) .eq. 1 .and. tofit(i_par+4) .eq. 1) then
                temp2 = (atemp(i_par+3)*atemp(i_par+3))+(atemp(i_par+4)*atemp(i_par+4)) ! NEW: sqrt(e)cosw,sqrt(e)sinw ==> e,w
                if (temp2 .lt. e_bounds(1, i_body) .or. temp2 .gt. e_bounds(2, i_body)) then
                    ! write(*,*)"---- par2kel_fit body ", i_body, ": e_bounds ", e_bounds(:,i_body), " temp2 = ",temp2, " BAD"
                    checkpar = .false.
                end if

                ! temp_kel(4) = temp2
                ! temp_kel(5) = mod(rad2deg*atan2(atemp(i_par+4), atemp(i_par+3))+circ, circ)
                ! if (temp_kel(4) .le. TOLERANCE)then
                !     temp_kel(5) = 90.0_dp
                ! end if

                ! case ecc <= TOLERANCE ==> ecc = 0, argp = 90°
                if (temp2 .le. TOLERANCE) then
                    temp_kel(4) = zero
                    temp_kel(5) = 90.0_dp
                else
                    temp_kel(4) = temp2
                    temp_kel(5) = mod(rad2deg*atan2(atemp(i_par+4), atemp(i_par+3))+circ, circ)
                end if

                ! comment this check for testing
                ! if (abs(temp_kel(4)-atemp(i_par+3)) .le. TOL_dp) then
                !     if (abs(atemp(i_par+4)) .gt. TOL_dp) then
                !         checkpar = .false.
                !     end if
                ! else if (abs(temp_kel(4)-atemp(i_par+4)) .le. TOL_dp) then
                !     if (abs(atemp(i_par+3)) .gt. TOL_dp) then
                !         checkpar = .false.
                !     end if
                ! end if

            end if

            !inc fit
            if (tofit(i_par+6) .eq. 1) then

                ! 2022-07-29 fit inc and longN as they are (in physical unit)
                ! ! cosicoslN,cosisinlN to inc,lN
                ! if (tofit(i_par+7) .eq. 1) then
                !     temp_kel(7) = acos(sqrt(atemp(i_par+6)*atemp(i_par+6)+atemp(i_par+7)*atemp(i_par+7)))*rad2deg
                !     temp_kel(8) = mod(rad2deg*atan2(atemp(i_par+7), atemp(i_par+6))+circ, circ)
                ! end if
                temp_kel(7) = atemp(i_par+6)

                ! inc out of [0, 180] not allowed
                if (temp_kel(7) .le. zero .or. temp_kel(7) .ge. 180.0_dp) then
                    ! write(*,*)"---- par2kel_fit body ", i_body, ": inc out of [0, 180] not allowed: ", temp_kel(7)
                    checkpar = .false.
                end if

            end if

            ! longN fit
            if (tofit(i_par+7) .eq. 1) then
                temp_kel(8) = atemp(i_par+7)
            end if

            ! lambda to mA
            if (tofit(i_par+5) .eq. 1) then
                temp_kel(6) = mod(mod((atemp(i_par+5)-temp_kel(5)-temp_kel(8)), circ)+circ, circ) ! lambda to mA
            end if

            mass(i_body) = temp_kel(1)
            radius(i_body) = temp_kel(2)
            period(i_body) = temp_kel(3)
            ! sma(i_body) = get_semax(mass(1), temp_kel(1), temp_kel(3))
            call period_to_sma(mass(1), temp_kel(1), temp_kel(3), sma(i_body))
            ecc(i_body) = temp_kel(4)
            argp(i_body) = temp_kel(5)
            meanA(i_body) = temp_kel(6)
            inc(i_body) = temp_kel(7)
            longN(i_body) = temp_kel(8)

        end do

        deallocate (atemp)

        return
    end subroutine par2kel_fit

    function checkbounds_angle(val, min_val, max_val) result(check)
        logical::check
        real(dp), intent(in)::val, min_val, max_val

        ! real(dp)::temp

        check = .true.

        ! check if val is within min_val,max_val
        if ((val .lt. min_val) .or.&
            &(val .gt. max_val)) then
            check = .false.

            !! 2022-07-07
            !! Commenting this part, I decided that the boundaries has to be defined/provided,
            !! not checking circular values such as angles (some fitting lambdas are outside 0-360 bounds ... )
            ! ! check if argp -360 is within min_val,max_val
            ! if(.not.check)then
            !     temp=val-circ
            !     if ((temp .ge. min_val) .and.&
            !         &(temp .le. max_val)) then
            !         check = .true.
            !     end if
            ! end if

            ! ! check if mod(argp, 360) > 0 is within min_val,max_val
            ! if(.not.check)then
            !     temp = mod(mod(val, circ) +  circ, circ)
            !     if ((temp .ge. min_val) .and.&
            !     &(temp .le. max_val)) then
            !         check = .true.
            !     end if
            ! end if

        end if

        return
    end function checkbounds_angle

!   check the physical bounds of the fitted parameters
    function checkbounds_fit(fit_parameters) result(check)
        real(dp), dimension(:), intent(in)::fit_parameters
        logical::check

        real(dp)::temp
        character::p1
        character(2)::p2
        integer::ifit

        check = .true.

        if (nfit .gt. 0) then
            do ifit = 1, nfit
                !      body=int((idall(j)-3)/8)+2

                ! check if fit_parameters are within boundaries
                if ((fit_parameters(ifit) .lt. minpar(ifit)) .or.& ! fit < min
                    &(fit_parameters(ifit) .gt. maxpar(ifit))) then ! fit > max

                    check = .false.
                    if (.not. check) then
                        p1 = parid(ifit) (1:1)
                        p2 = parid(ifit) (1:2)
                        if ((p1 .eq. 'w') .or.&
                            &(p2 .eq. 'mA') .or.&
                            &(p2 .eq. 'la'))&
                            &then
                            temp = mod(mod(fit_parameters(ifit), circ)+circ, circ)
                            ! write(*,*)parid(ifit), temp, minpar(ifit), maxpar(ifit)
                            if ((temp .ge. minpar(ifit)) .and.&
                                &(temp .le. maxpar(ifit))) then
                                check = .true.
                            end if

                        end if
                    end if
                end if

                if (.not. check) return
            end do
        end if

        return
    end function checkbounds_fit

    function checkbounds_kel(mass, radius, period, ecc, argp, meanA, inc, longN) result(check)
        logical::check
        real(dp), dimension(:), intent(in)::mass, radius, period, ecc, argp, meanA, inc, longN

        integer::j, j1

        check = .true.

        ! write(*,*)'Ms'
        if ((mass(1) .lt. par_min(1)) .or.&
            &(mass(1) .gt. par_max(1))) then
            check = .false.
            return
        end if

        ! write(*,*)'Rs'
        if ((radius(1) .lt. par_min(2)) .or.&
            &(radius(1) .gt. par_max(2))) then
            check = .false.
            return
        end if

        do j = 2, NB
            j1 = (j-2)*8

            ! write(*,*)'Mp',j,': ',mass(j),par_min(j1+3),par_max(j1+3)
            ! mass
            if ((mass(j) .lt. par_min(j1+3)) .or.&
                &(mass(j) .gt. par_max(j1+3))) then
                check = .false.
                return
            end if

            ! write(*,*)'R',j,': ',radius(j),par_min(j1+4),par_max(j1+4)
            ! radius
            if ((radius(j) .lt. par_min(j1+4)) .or.&
                &(radius(j) .gt. par_max(j1+4))) then
                check = .false.
                return
            end if

            ! write(*,*)'P',j,': ',period(j),par_min(j1+5),par_max(j1+5)
            ! period
            if ((period(j) .lt. par_min(j1+5)) .or.&
                &(period(j) .gt. par_max(j1+5))) then
                check = .false.
                return
            end if

            ! write(*,*)'e',j,': ',ecc(j),par_min(j1+6),par_max(j1+6)
            ! eccentricity
            if ((ecc(j) .lt. par_min(j1+6)) .or.&
                &(ecc(j) .gt. par_max(j1+6))) then
                check = .false.
                return
            end if

            ! argument of pericentre
            ! write(*,*)'argp',j,': ',argp(j),par_min(j1+7),par_max(j1+7)
            check = checkbounds_angle(argp(j), par_min(j1+7), par_max(j1+7))
            ! write(*,*)"argp", check
            if (.not. check) return

            ! mean anomaly
            ! write(*,*)'meanA',j,': ',meanA(j),par_min(j1+8),par_max(j1+8)
            check = checkbounds_angle(meanA(j), par_min(j1+8), par_max(j1+8))
            ! write(*,*)"meanA", check
            if (.not. check) return

            ! write(*,*)'inc',j,': ',inc(j),par_min(j1+9),par_max(j1+9)
            ! inclination
            if ((inc(j) .lt. par_min(j1+9)) .or.&
                &(inc(j) .gt. par_max(j1+9))) then
                check = .false.
                return
            end if

            ! longitude of node
            ! write(*,*)'longN',j,': ',longN(j),par_min(j1+10),par_max(j1+10)
            check = checkbounds_angle(longN(j), par_min(j1+10), par_max(j1+10))
            ! write(*,*)"longN", check
            if (.not. check) return

        end do

        return
    end function checkbounds_kel

    ! ------------------------------------------------------------------ !
    subroutine convert_parameters(allpar, par,&
        &mass, radius, period, sma, ecc, argp, meanA, inc, longN, checkpar)
        real(dp), dimension(:), intent(in)::allpar, par
        real(dp), dimension(:), intent(out)::mass, radius
        real(dp), dimension(:), intent(out)::period, sma, ecc, argp, meanA, inc, longN
        logical, intent(inout)::checkpar

        checkpar = checkbounds_fit(par)
        if (checkpar) then
            call par2kel_fit(allpar, par, mass, radius,&
                &period, sma, ecc, argp, meanA, inc, longN, checkpar)
            if (checkpar) then
                checkpar = checkbounds_kel(mass, radius,&
                    &period, ecc, argp, meanA, inc, longN)
            end if
        end if

        return
    end subroutine convert_parameters

    subroutine only_convert_parameters(allpar, par,&
        &mass, radius, period, sma, ecc, argp, meanA, inc, longN, checkpar)
        real(dp), dimension(:), intent(in)::allpar, par
        real(dp), dimension(:), intent(out)::mass, radius
        real(dp), dimension(:), intent(out)::period, sma, ecc, argp, meanA, inc, longN
        logical, intent(inout)::checkpar

        call par2kel_fit(allpar, par, mass, radius,&
            &period, sma, ecc, argp, meanA, inc, longN, checkpar)

        return
    end subroutine only_convert_parameters
    ! ------------------------------------------------------------------ !

    ! subroutine needed to update the all_parameters from the fitted fit_parameters
    subroutine update_parameters_fit2all(fit_parameters, all_parameters)
        real(dp), dimension(:), intent(in)::fit_parameters
        real(dp), dimension(:), intent(inout)::all_parameters

        real(dp), dimension(:), allocatable::mass, radius
        real(dp), dimension(:), allocatable::period, sma, ecc, argp, meanA, inc, longN
        real(dp), dimension(:), allocatable::all_updating
        logical::check

        check = .true.
        allocate (mass(NB), radius(NB))
        allocate (period(NB), sma(NB), ecc(NB), argp(NB), meanA(NB), inc(NB), longN(NB))
        call par2kel_fit(all_parameters, fit_parameters, mass, radius,&
            &period, sma, ecc, argp, meanA, inc, longN, check)
        call set_all_param(mass, radius, period, ecc, argp, meanA, inc, longN, all_updating)
        all_parameters = all_updating
        deallocate (all_updating, mass, radius, period, sma, ecc, argp, meanA, inc, longN)

        return
    end subroutine update_parameters_fit2all

    ! function for pso/pik to initialize properly the first-generation population
    function check_only_boundaries(all_parameters, fit_parameters) result(check)
        logical::check
        real(dp), dimension(:), intent(in)::all_parameters, fit_parameters
        real(dp), dimension(:), allocatable::mass, radius,&
            &period, sma, ecc, argp, meanA, inc, longN

        check = .true.
        check = checkbounds_fit(fit_parameters)
        ! write(*,'(1x,a,l2)')' DEBUG: checkbounds_fit = ',check
        if (check) then
            allocate (mass(NB), radius(NB),&
                &period(NB), sma(NB), ecc(NB), argp(NB), meanA(NB), inc(NB), longN(NB))
            call par2kel_fit(all_parameters, fit_parameters,&
                &mass, radius,&
                &period, sma, ecc, argp, meanA, inc, longN,&
                &check)
            ! write(*,'(1x,a,l2)')' DEBUG: par2kel_fit = ',check
            if (check) then
                check = checkbounds_kel(mass, radius,&
                    &period, ecc, argp, meanA, inc, longN)
                ! write(*,'(1x,a,l2)')' DEBUG: checkbounds_kel = ',check
            end if
            deallocate (mass, radius, period, sma, ecc, argp, meanA, inc, longN)
        end if

        return
    end function check_only_boundaries

    ! boundaries check that output a scale parameter and not a logical .true./.false.
    ! the scale parameter will be multiplied to the weighted residuals to worsening the fitness
    !

    function get_fit_scale(value, minvalue, maxvalue) result(xscale)
        real(dp)::xscale
        real(dp), intent(in)::value, minvalue, maxvalue
        real(dp)::delta

        xscale = zero
        delta = abs(maxvalue-minvalue)

        ! check if fit_parameters are within boundaries
        if (value .lt. minvalue) then
            xscale = abs((value-minvalue)/delta)*100.0_dp
        else if (value .gt. maxvalue) then
            xscale = abs((value-maxvalue)/delta)*100.0_dp
        else
            xscale = zero
        end if

        return
    end function get_fit_scale

    function checkbounds_fit_scale(fit_parameters) result(bound_scale)
        real(dp)::bound_scale
        real(dp), dimension(:), intent(in)::fit_parameters
        integer::j
        real(dp)::xscale

        bound_scale = zero
        xscale = zero

        do j = 1, nfit
            ! check if fit_parameters are within boundaries
            xscale = get_fit_scale(fit_parameters(j), minpar(j), maxpar(j))
            bound_scale = bound_scale+xscale
        end do

        return
    end function checkbounds_fit_scale

    function checkbounds_kel_scale(mass, radius, period, ecc, argp, meanA, inc, longN) result(kel_scale)
        real(dp)::kel_scale
        real(dp), dimension(:), intent(in)::mass, radius, period, ecc, argp, meanA, inc, longN

        real(dp)::temp, xscale, xtemp1, xtemp2
        integer::j, j1

        kel_scale = zero
        xscale = zero

        ! Mstar
        xscale = get_fit_scale(mass(1), par_min(1), par_max(1))
        kel_scale = kel_scale+xscale

        ! Rstar
        xscale = get_fit_scale(radius(1), par_min(2), par_max(2))
        kel_scale = kel_scale+xscale

        do j = 2, NB
            j1 = (j-2)*8

            ! mass
            xscale = get_fit_scale(mass(j), par_min(j1+3), par_max(j1+3))
            kel_scale = kel_scale+xscale

            ! radius
            xscale = get_fit_scale(radius(j), par_min(j1+4), par_max(j1+4))
            kel_scale = kel_scale+xscale

            ! period
            xscale = get_fit_scale(period(j), par_min(j1+5), par_max(j1+5))
            kel_scale = kel_scale+xscale

            ! eccentricity
            xscale = get_fit_scale(ecc(j), par_min(j1+6), par_max(j1+6))
            kel_scale = kel_scale+xscale

            ! argument of pericentre
            temp = argp(j)-circ
            xtemp1 = get_fit_scale(temp, par_min(j1+7), par_max(j1+7))
            xtemp2 = get_fit_scale(argp(j), par_min(j1+7), par_max(j1+7))
            if ((xtemp1 .gt. zero) .and. (xtemp2 .gt. zero)) then
                xscale = max(xtemp1, xtemp2)
            else
                xscale = zero
            end if
            kel_scale = kel_scale+xscale

            ! mean anomaly
            temp = meanA(j)-circ
            xtemp1 = get_fit_scale(temp, par_min(j1+8), par_max(j1+8))
            xtemp2 = get_fit_scale(meanA(j), par_min(j1+8), par_max(j1+8))
            if ((xtemp1 .gt. zero) .and. (xtemp2 .gt. zero)) then
                xscale = max(xtemp1, xtemp2)
            else
                xscale = zero
            end if
            kel_scale = kel_scale+xscale

            ! inclination
            xscale = get_fit_scale(inc(j), par_min(j1+9), par_max(j1+9))
            kel_scale = kel_scale+xscale

            ! longitude of node
            temp = longN(j)-circ
            xtemp1 = get_fit_scale(temp, par_min(j1+10), par_max(j1+10))
            xtemp2 = get_fit_scale(longN(j), par_min(j1+10), par_max(j1+10))
            if ((xtemp1 .gt. zero) .and. (xtemp2 .gt. zero)) then
                xscale = max(xtemp1, xtemp2)
            else
                xscale = zero
            end if
            kel_scale = kel_scale+xscale

        end do

        return
    end function checkbounds_kel_scale

    subroutine convert_parameters_scale(all_parameters, fit_parameters,&
        &mass, radius, period, sma, ecc, argp, meanA, inc, longN, fit_scale)
        real(dp), dimension(:), intent(in)::all_parameters, fit_parameters
        real(dp), dimension(:), intent(out)::mass, radius
        real(dp), dimension(:), intent(out)::period, sma, ecc, argp, meanA, inc, longN
        real(dp), intent(out)::fit_scale
        real(dp)::bound_scale, kel_scale
        logical::checkpar

        checkpar = .true.
        fit_scale = zero
        bound_scale = checkbounds_fit_scale(fit_parameters)
        call par2kel_fit(all_parameters, fit_parameters,&
            &mass, radius, period, sma, ecc, argp, meanA, inc, longN, checkpar)
        kel_scale = checkbounds_kel_scale(mass, radius, period, ecc, argp, meanA, inc, longN)
        fit_scale = one+bound_scale+kel_scale ! min value of fit_scale has to be 1

        return
    end subroutine convert_parameters_scale

    function check_only_boundaries_scale(all_parameters, fit_parameters) result(fit_scale)
        real(dp)::fit_scale
        real(dp), dimension(:), intent(in)::all_parameters, fit_parameters

        real(dp), dimension(:), allocatable::mass, radius
        real(dp), dimension(:), allocatable::period, sma, ecc, argp, meanA, inc, longN
        real(dp)::bound_scale, kel_scale
        logical::checkpar

        checkpar = .true.
        bound_scale = checkbounds_fit_scale(fit_parameters)
        allocate (mass(NB), radius(NB))
        allocate (period(NB), sma(NB), ecc(NB), argp(NB), meanA(NB), inc(NB), longN(NB))
        call par2kel_fit(all_parameters, fit_parameters,&
            &mass, radius, period, sma, ecc, argp, meanA, inc, longN, checkpar)
        kel_scale = checkbounds_kel_scale(mass, radius, period, ecc, argp, meanA, inc, longN)
        deallocate (mass, radius, period, sma, ecc, argp, meanA, inc, longN)
        fit_scale = one+bound_scale+kel_scale

        return
    end function check_only_boundaries_scale
    ! ------------------------------------------------------------------ !
    ! given the computed parameters with the L-M it adjusts some parameters
    ! i.e. setting all the angles between 0 and 360 (no < 0 angles) etc.
    subroutine param_adj(par, sigpar)
        real(dp), dimension(:), intent(inout)::par, sigpar
        integer::j1
        real(dp), parameter::hcirc = 180._dp

        do j1 = 1, nfit

            if ((parid(j1) (1:1) .eq. 'w') .or.&
               &(parid(j1) (1:2) .eq. 'mA') .or.&
               &(parid(j1) (1:2) .eq. 'lN') .or.&
               &(parid(j1) (1:2) .eq. 'la')) then
                par(j1) = mod(par(j1), circ)
                sigpar(j1) = mod(sigpar(j1), circ)
            end if

        end do

        return
    end subroutine param_adj
    ! ------------------------------------------------------------------ !

    ! ------------------------------------------------------------------ !
    ! conversion from parameter boundaries [0, 1] --> [ParMin, ParMax]
    subroutine norm2par(norm, par)
        real(dp), dimension(:), intent(in)::norm
        real(dp), dimension(:), intent(out)::par
        real(dp)::dpar
        integer::j

        do j = 1, nfit
            dpar = abs(maxpar(j)-minpar(j))
            par(j) = minpar(j)+dpar*norm(j)
        end do

        return
    end subroutine norm2par

end module parameters_conversion
