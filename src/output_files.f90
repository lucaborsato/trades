module output_files
  use constants
  use parameters
  use parameters_conversion
  use init_trades,only:get_unit,state2string
  use convert_type,only:string
  use celestial_mechanics,only:barycenter
  implicit none

  interface write_RV
    module procedure write_RV_s,write_RV_f
  end interface write_RV
  
  interface write_T0
    module procedure write_T0_s,write_T0_f
  end interface write_T0

  interface write_par
    module procedure write_par_s,write_par_f
  end interface write_par

  contains

  ! THIS IS VERY OLD....
  ! it writes to screen and into file the inital combination of the grid search and the output
  ! of the LM, in particular:
  ! simID = number that identify the simulation
  ! LM_info = flag for LM convergence
  ! Chi2, Chi2r = Chi square and reduced chi square of the simulation with those initial parameters
  subroutine write_simlst(cpuid,Ngrid,grid,infos)
    integer,intent(in)::cpuid,Ngrid
    real(dp),dimension(:,:),intent(in)::grid
    integer,dimension(:),intent(in)::infos
    integer::ulst,j
    character(512)::fllst
    character(80)::fmt

    fllst=trim(path)//"sim.lst"
    fllst=trim(adjustl(fllst))
    ulst=get_unit(cpuid)
    open(ulst,file=trim(fllst))
    write(ulst,'(a,2x,a,8x,a,22x,a,22x,a,21x,a,25x,a,22x,a,20x,a)')&
        &"#    simID","LM_info","m[Msun]","P[d]","a[AU]","e","w[°]",&
        &"fitness","fitness*dof"
    write(*,'(a,2x,a,8x,a,22x,a,22x,a,21x,a,25x,a,22x,a,20x,a)')&
        "#     simID","LM_info","m[Msun]","P[d]","a[AU]","e","w[°]",&
        &"fitness","fitness*dof"
    fmt=adjustl("(i6,1x,i3,1x,1000("//trim(sprec)//",1x))")
    do j=1,Ngrid
      write(ulst,trim(fmt))j,infos(j),grid(:,j),(grid(6,j)*real(dof,dp))
      write(*,trim(fmt))j,infos(j),grid(:,j),(grid(6,j)*real(dof,dp))
    end do
    close(ulst)

    return
  end subroutine write_simlst
  
  
  ! GOOD GRID WRITE SUMMARY
  subroutine write_grid_summary(cpuid,n_grid,lm_flag,grid_summary)
    integer,intent(in)::cpuid,n_grid,lm_flag
    real(dp),dimension(:,:),intent(in)::grid_summary
    
    character(512)::output_file,header
    character(80)::fmt_wrt
    integer::u_wrt,i_sim
    
    u_wrt=get_unit(cpuid)
    
    fmt_wrt = "(i6,1x,1000("//trim(adjustl(sprec))//",1x))"
    output_file=trim(path)//"summary_grid_sims_"//&
      &trim(adjustl(string(lm_flag)))//".dat"
    header="# id_sim "//trim(adjustl(all_names_str))//" fitness_x_dof fitness"
    
    open(u_wrt,file=trim(output_file))
    write(u_wrt,'(a)')trim(header)
    do i_sim=1,n_grid
      write(u_wrt,trim(fmt_wrt))i_sim,grid_summary(i_sim,:)
    end do
    close(u_wrt)
    write(*,'(a,a)')" WRITTEN FILE:",trim(output_file)
    flush(6)
  
    return
  end subroutine write_grid_summary

  ! format string to write orbit file
  function fmtorbit() result(fmt)
    integer::ncol
    character(128)::fmt
    character(72)::col
    character(5)::nstr

    col="(1x,"//trim(sprec)//")"
    ncol=2+(NB*6)+1
    nstr=trim(adjustl(string(ncol)))
    fmt="("//trim(nstr)//trim(col)//")"
    fmt=trim(adjustl(fmt))

    return
  end function fmtorbit

  ! format string to write const. of motion file
  function fmtconst() result(fmt)
    implicit none
    character(128)::fmt

    fmt="(5(1x,"//trim(sprec)//"))"

    return
  end function fmtconst

  ! format string to write Keplerian elements file
  function fmtele() result(fmt)
    implicit none
    character(128)::fmt

    fmt="(10(1x,"//trim(sprec)//"))"

    return
  end function fmtele

  ! ------------------------------------------------------------------ !
  ! write to screen input parameters
  subroutine write_lm_inputpar(cpuid,par)
    integer,intent(in)::cpuid
    real(dp),dimension(:),intent(in)::par
    integer::i
    character(80)::fmt

    write(*,'(a,i5,a)')" CPU ",cpuid," PARAM: "
    fmt=adjustl("(1x,a6,"//trim(sprec)//")")
    do i=1,nfit
      write(*,trim(fmt))parid(i),par(i)
    end do

    return
  end subroutine write_lm_inputpar
  ! ------------------------------------------------------------------ !

  ! subroutine to write initial parameters of planets in a file
  ! during the call ode_out(...)
  subroutine outElements(isim,wrtid,m,R,P,a,e,w,mA,inc,lN)
    !$ use omp_lib
    integer,intent(in)::isim,wrtid
    real(dp),dimension(:),intent(in)::m,R,P,a,e,w,mA,inc,lN
    integer::cpuid,uwrt,ii
    character(512)::fwrt,fmtw

    fmtw=adjustl('(1000('//trim(sprec)//'))')
    fwrt=trim(path)//trim(adjustl(string(isim)))//"_"//&
      &trim(adjustl(string(wrtid)))//"_initialElements.dat"
    fwrt=trim(adjustl(fwrt))
    cpuid=1
    !$ cpuid=omp_get_thread_num()+1
    uwrt=get_unit(cpuid)
    open(uwrt,file=trim(fwrt))
    write(uwrt,'(a)')"# M_Msun R_Rsun P_d a_AU e w_deg mA_deg inc_deg lN_deg"
    !write(*,'(a)')""
    !write(*,'(a)')" Initial orbital parameters for each planet"
    !write(*,'(a)')"# M_Msun R_Rsun P_d a_AU e w_deg mA_deg inc_deg lN_deg"
    do ii=2,NB
      write(uwrt,trim(fmtw))m(ii),R(ii),P(ii),a(ii),e(ii),&
          &w(ii),mA(ii),inc(ii),lN(ii)
      !write(*,trim(fmtw))m(ii),R(ii),P(ii),a(ii),e(ii),&
      !     &w(ii),mA(ii),i(ii),lN(ii)
    end do
    close(uwrt)
    !write(*,'(a)')""

    return
  end subroutine outElements
  ! ------------------------------------------------------------------ !

  ! ------------------------------------------------------------------ !
  ! subroutine to store the orbital state vectors, ready to be written into file
  subroutine store_orb(pos,itime,m,rw,storeorb)
    integer,intent(in)::pos
    real(dp),intent(in)::itime
    real(dp),dimension(:),intent(in)::m,rw
    real(dp),dimension(:,:),intent(inout)::storeorb
    real(dp),dimension(6)::bar
    real(dp),dimension(:),allocatable::rbar

    allocate(rbar(NBDIM))
    storeorb(1,pos)=tepoch+itime
    call barycenter(m,rw,bar,rbar)
    storeorb(2,pos)=-rbar(3)/speedaud
    storeorb(3:NBDIM+2,pos)=rw
    storeorb(NBDIM+3,pos)=(-rbar(6)*AU/s24h)
    deallocate(rbar)

    return
  end subroutine store_orb

  ! subroutine to store constants of motion, ready to be written into file
  subroutine store_con(pos,itime,Etot,Eold,htot,hold,storecon)
    integer,intent(in)::pos
    real(dp),intent(in)::itime,Etot,Eold,htot,hold
    real(dp),dimension(:,:),intent(inout)::storecon

    storecon(1,pos)=tepoch+itime
    storecon(2,pos)=htot
    storecon(3,pos)=htot-hold
    storecon(4,pos)=Etot
    storecon(5,pos)=Etot-Eold

    return
  end subroutine store_con
  ! ------------------------------------------------------------------ !

  ! ------------------------------------------------------------------ !
  ! write orbit state vector into file
  subroutine write_orb(maxpos,uorb,fmorb,storeorb)
    integer,intent(in)::uorb,maxpos
    character(*),intent(in)::fmorb
    real(dp),dimension(:,:),intent(inout)::storeorb
    integer::j

    do j=1,maxpos
      write(uorb,trim(adjustl(fmorb)))storeorb(:,j)
    end do

    return
  end subroutine write_orb

  ! generic subroutine to write something into file
  subroutine write_file(maxpos,uwrt,fmwrt,store)
    integer,intent(in)::uwrt,maxpos
    character(*),intent(in)::fmwrt
    real(dp),dimension(:,:),intent(inout)::store
    integer::j

    do j=1,maxpos
      write(uwrt,trim(adjustl(fmwrt)))store(:,j)
    end do

    return
  end subroutine write_file

  ! write keplerian orbital elements into file
  subroutine write_elem(pos,uele,fmele,m,storeorb)
    use celestial_mechanics,only:elements
    integer,intent(in)::pos
    integer,dimension(:),intent(in)::uele
    character(*),intent(in)::fmele
    real(dp),dimension(:),intent(in)::m
    real(dp),dimension(:,:),intent(in)::storeorb
    real(dp),dimension(:),allocatable::P,a,e,i,mA,w,f,lN,dttau
    integer::j1,j2

    allocate(P(NB),a(NB),e(NB),i(NB),mA(NB),w(NB),f(NB),lN(NB),dttau(NB))
    j1loop: do j1=1,pos
      call elements(m,storeorb(3:NBDIM+2,j1),P,a,e,i,mA,w,f,lN,dttau)
      j2loop: do j2=2,NB
        write(uele(j2),trim(adjustl(fmele)))storeorb(1,j1),&
            &P(j2),a(j2),e(j2),i(j2),mA(j2),w(j2),lN(j2),&
            &f(j2),dttau(j2)
      end do j2loop
    end do j1loop
    deallocate(P,a,e,i,mA,w,f,lN,dttau)

    return
  end subroutine write_elem

  ! write transit time, light travel-time effect and state vector (at transit time)
  ! into file
  subroutine write_tra(pos,utra,stat_tra,storetra)
    integer,intent(in)::pos
    integer,dimension(:),intent(in)::utra
    integer,dimension(:,:),intent(in)::stat_tra
    real(dp),dimension(:,:),intent(in)::storetra
    integer::j1,j2
    character(80)::fmt

    fmt=adjustl("(1000("//trim(sprec)//",1x))") ! '(1000(es23.16,1x))'
    j1loop: do j1=1,pos
      j2loop: do j2=2,NB
        if(stat_tra(j2,j1).eq.1)then
          write(utra(j2),trim(fmt))storetra(:,j1)
        end if
      end do j2loop
    end do j1loop

    return
  end subroutine write_tra
  ! ------------------------------------------------------------------ !

    ! ------------------------------------------------------------------ !
  ! defines the unit and the name of orbit file
  subroutine set_file_orb(cpuid,isim,wrtid,uorb,florb)
    integer,intent(in)::cpuid,isim,wrtid
    integer,intent(inout)::uorb
    character(512),intent(inout)::florb
    character(512)::strState

    uorb=0
    uorb=get_unit(cpuid)
    florb=trim(path)//trim(adjustl(string(isim)))//"_"//&
      &trim(adjustl(string(wrtid)))//"_rotorbit.dat"
    florb=trim(adjustl(florb))
    open(uorb,file=trim(florb))
    strState=trim(adjustl(state2string(NB)))
    write(uorb,'(a,a,a)')"# Time_JD LTE_d ",trim(strState)," RV_ms^-1"

    return
  end subroutine set_file_orb

  ! defines the unit and the name of constants file
  subroutine set_file_con(cpuid,isim,wrtid,ucon,flcon)
    integer,intent(in)::cpuid,isim,wrtid
    integer,intent(inout)::ucon
    character(512),intent(inout)::flcon

    ucon=0
    ucon=get_unit(cpuid)
    flcon=trim(path)//trim(adjustl(string(isim)))//"_"//&
      &trim(adjustl(string(wrtid)))//"_constants.dat"
    flcon=trim(adjustl(flcon))
    open(ucon,file=trim(flcon))
    write(ucon,'(a)')"# Time htot dh Etot dE"

    return
  end subroutine set_file_con

  ! defines the unit and the name of keplerian orbital elements files
  subroutine set_file_elem(cpuid,isim,wrtid,uele,flele)
    integer,intent(in)::cpuid,isim,wrtid
    integer,dimension(:),allocatable,intent(inout)::uele
    character(512),dimension(:),allocatable,intent(inout)::flele
    character(512)::strElem
    integer::je

    allocate(uele(NB),flele(NB))
    uele=0
    do je=2,NB
      uele(je)=get_unit(cpuid)
      flele(je)=trim(path)//trim(adjustl(string(isim)))//"_"//&
        &trim(adjustl(string(wrtid)))//"_NB"//&
        &trim(adjustl(string(je)))//"_elements.dat"
      flele(je)=trim(adjustl(flele(je)))
      open(uele(je),file=trim(flele(je)))
      strElem="# Time "//&
          &"P"//trim(adjustl(string(je)))//"_d "//&
          &"a"//trim(adjustl(string(je)))//"_AU "//&
          &"e"//trim(adjustl(string(je)))//" "//&
          &"i"//trim(adjustl(string(je)))//"_deg "//&
          &"mA"//trim(adjustl(string(je)))//"_deg "//&
          &"w"//trim(adjustl(string(je)))//"_deg "//&
          &"lN"//trim(adjustl(string(je)))//"_deg "//&
          &"f"//trim(adjustl(string(je)))//"_deg "//&
          &"dtau"//trim(adjustl(string(je)))//"_d "
      write(uele(je),'(a)')trim(adjustl(strElem))
    end do

    return
  end subroutine set_file_elem

  ! closes and deallocates units and files for orbital elements
  subroutine close_elem(uele,flele)
    integer,dimension(:),allocatable,intent(inout)::uele
    character(512),dimension(:),allocatable,intent(inout)::flele
    integer::je

    do je=2,NB
      close(uele(je))
    end do
    deallocate(uele,flele)

    return
  end subroutine close_elem

  ! defines the unit and the name of transit files
  subroutine set_file_tra(cpuid,isim,wrtid,utra,fltra)
    integer,intent(in)::cpuid,isim,wrtid
    integer,dimension(:),allocatable,intent(inout)::utra
    character(512),dimension(:),allocatable,intent(inout)::fltra
    character(512)::strState,strTra
    integer::j

    strState=state2string(NB)
    allocate(utra(NB),fltra(NB))
    utra=0
    fltra=""
    do j=2,NB
      utra(j)=get_unit(cpuid)
      fltra=trim(path)//trim(adjustl(string(isim)))//"_"//&
        &trim(adjustl(string(wrtid)))//"_NB"//&
        &trim(adjustl(string(j)))//"_tra.dat"
      fltra(j)=trim(adjustl(fltra(j)))
      open(utra(j),file=trim(fltra(j)))
      strTra="# ttra_"//trim(adjustl(string(j)))//&
          &" LTE_"//trim(adjustl(string(j)))//&
          &" t1_"//trim(adjustl(string(j)))//&
          &" t2_"//trim(adjustl(string(j)))//&
          &" t3_"//trim(adjustl(string(j)))//&
          &" t4_"//trim(adjustl(string(j)))//&
          &" "//trim(adjustl(strState))
      write(utra(j),'(a)')trim(adjustl(strTra))
    end do

    return
  end subroutine set_file_tra

  ! closes and deallocates units and files for transit times
  subroutine close_tra(utra,fltra)
    integer,dimension(:),allocatable,intent(inout)::utra
    character(512),dimension(:),allocatable,intent(inout)::fltra
    integer::j

    do j=2,NB
      close(utra(j))
    end do
    deallocate(utra,fltra)

    return
  end subroutine close_tra
  ! ------------------------------------------------------------------ !

  subroutine setGammaWrt(gammaIn, gammaOut)
    real(dp),dimension(:,:),intent(in)::gammaIn
    real(dp),dimension(:,:),intent(out)::gammaOut
    integer::j,aRV,bRV
    aRV=0
    bRV=0
    do j=1,nRVset
      aRV=bRV+1
      bRV=bRV+nRVsingle(j)
      gammaOut(aRv:bRV,1)=gammaIn(j,1)
      gammaOut(aRv:bRV,2)=gammaIn(j,2)
    end do
    return
  end subroutine setGammaWrt
  
  subroutine setRVok(RV_sim,gamma,RVok)
    real(dp),dimension(:),intent(in)::RV_sim
    real(dp),dimension(:,:),intent(in)::gamma
    real(dp),dimension(:),intent(out)::RVok
    integer::j,aRV,bRV
    aRV=0
    bRV=0
    do j=1,nRVset
      aRV=bRV+1
      bRV=bRV+nRVsingle(j)
      RVok(aRv:bRV)=RV_sim(aRv:bRV)+gamma(j,1)
!       write(*,*)" SELECTED RV INTERVALS: ",aRV," -- ",bRV," adding gamma = ",gamma(j,1)
!       write(*,*)" RV_sim: "
!       write(*,*)RV_sim(aRv:bRV)
!       write(*,*)" RVok: "
!       write(*,*)RVok(aRv:bRV)
!       write(*,*)
    end do
    return
  end subroutine setRVok
  
  ! set RV data and gamma to write to screen and files
  subroutine setWriteRV(RV_sim,gamma,RV_simwrt,gamma_wrt)
    real(dp),dimension(:),intent(in)::RV_sim
    real(dp),dimension(:,:),intent(in)::gamma
    real(dp),dimension(:),intent(out)::RV_simwrt
    real(dp),dimension(:,:),intent(out)::gamma_wrt
    RV_simwrt=zero
    gamma_wrt=zero
    call setGammaWrt(gamma,gamma_wrt)!use this subroutine to set gamma_wrt
    call setRVok(RV_sim,gamma,RV_simwrt)
    return
  end subroutine setWriteRV
  
  ! ------------------------------------------------------------------ !
  ! write to screen the RV simulated compared to the RV observated
  subroutine write_RV_s(gamma,RV_sim,RV_stat)
    real(dp),dimension(:,:),intent(in)::gamma
    real(dp),dimension(:),intent(in)::RV_sim
    integer,dimension(:),intent(in)::RV_stat
    real(dp),dimension(:),allocatable::RV_simwrt
    real(dp),dimension(:,:),allocatable::gamma_wrt
    integer::j
    character(80)::fmt

    allocate(RV_simwrt(nRV),gamma_wrt(nRV,2))
    call setWriteRV(RV_sim,gamma,RV_simwrt,gamma_wrt)
!     write(*,'(2(a,es23.16))')" RV -> gamma = ",gamma(1),&
!         &" +/- ",gamma(2)
    write(*,'(a)')"# JD RVobs eRVobs rv_sim RV_sim gamma e_gamma RVsetID &
    &RV_stat"
    fmt=adjustl("(7("//trim(sprec)//",1x),i3,1x,i3)")
    do j=1,nRV
      write(*,trim(fmt))jdRV(j),RVobs(j),eRVobs(j),RV_sim(j),&
      &RV_simwrt(j),gamma_wrt(j,1),gamma_wrt(j,2),RVsetID(j),&
      &RV_stat(j)
    end do
    deallocate(RV_simwrt,gamma_wrt)

    return
  end subroutine write_RV_s

  ! write into file the RV simulated compared to the RV observated
  subroutine write_RV_f(cpuid,isim,wrtid,gamma,RV_sim,RV_stat)
    integer,intent(in)::cpuid,isim,wrtid ! wrtid = if 0 original parameters before LM, 1 parameters after LM
    real(dp),dimension(:,:),intent(in)::gamma
    real(dp),dimension(:),intent(in)::RV_sim
    integer,dimension(:),intent(in)::RV_stat
    real(dp),dimension(:),allocatable::RV_simwrt
    real(dp),dimension(:,:),allocatable::gamma_wrt
    integer::uRV,j
    character(512)::flRV
    character(80)::fmt

    allocate(RV_simwrt(nRV),gamma_wrt(nRV,2))
    call setWriteRV(RV_sim,gamma,RV_simwrt,gamma_wrt)

    uRV=get_unit(cpuid)
    flRV=""
    flRV=trim(path)//trim(adjustl(string(isim)))//"_"//&
      &trim(adjustl(string(wrtid)))//"_simRV.dat"
    flRV=trim(adjustl(flRV))
    open(uRV,file=trim(flRV))
!     fmt=adjustl("(a,28x,a,19x,a,18x,a,10x,2(a,"//trim(sprec)//"))")
    write(uRV,'(a)')"# JD RVobs eRVobs rv_sim RV_sim gamma e_gamma RVsetID &
    &RV_stat"
    fmt=""
    fmt=adjustl("(7("//trim(sprec)//",1x),i3,1x,i3)")
    do j=1,nRV
      write(uRV,trim(fmt))jdRV(j),RVobs(j),eRVobs(j),RV_sim(j),&
      &RV_simwrt(j),gamma_wrt(j,1),gamma_wrt(j,2),RVsetID(j),&
      &RV_stat(j)
    end do
    close(uRV)
    deallocate(RV_simwrt,gamma_wrt)
    write(*,'(a,a)')" WRITTEN RV INTO FILE: ",trim(flRV)

    return
  end subroutine write_RV_f
  ! ------------------------------------------------------------------ !

  ! ------------------------------------------------------------------ !
  ! write to screen the T0 simulated compared to the T0 observated
  subroutine write_T0_s(T0_sim,T0_stat)
    real(dp),dimension(:,:),intent(in)::T0_sim
    integer,dimension(:,:),intent(in)::T0_stat
    integer::j,j1
    character(80)::fmt

    fmt=adjustl("(i6,4("//trim(sprec)//",1x),i3)")
    do j=2,NB
      if(nT0(j).ne.0)then
        write(*,'(a,1x,a,19x,a,18x,a,18x,a,5x,a)')&
            "# epoT0obs "," T0obs "," eT0obs "," T0_sim ",&
            &" T0obs-T0_sim "," T0_stat "
        do j1=1,nT0(j)
          write(*,trim(fmt))epoT0obs(j1,j),T0obs(j1,j),eT0obs(j1,j),&
              &T0_sim(j1,j),(T0obs(j1,j)-T0_sim(j1,j)),T0_stat(j1,j)
        end do
        write(*,*)""
      end if
    end do

    return
  end subroutine write_T0_s

  ! write into file the T0 simulated compared to the T0 observated
  subroutine write_T0_f(cpuid,isim,wrtid,T0_sim,T0_stat)
    integer,intent(in)::cpuid,isim,wrtid ! wrtid = if 0 original parameters before LM, 1 parameters after LM
    real(dp),dimension(:,:),intent(in)::T0_sim
    integer,dimension(:,:),intent(in)::T0_stat
    character(512)::flT0
    integer::uT0,j,j1
    character(80)::fmt

    fmt=adjustl("(i6,4("//trim(sprec)//",1x),i3)")
    do j=2,NB
      if(nT0(j).ne.0)then
        uT0=get_unit(cpuid)
        flT0=""
        flT0=trim(path)//trim(adjustl(string(isim)))//"_"//&
          &trim(adjustl(string(wrtid)))//"_NB"//&
          &trim(adjustl(string(j)))//"_simT0.dat"
        flT0=trim(adjustl(flT0))
        open(uT0,file=trim(flT0))
        write(uT0,'(a,1x,a,19x,a,18x,a,18x,a,5x,a)')&
            &"# epoT0obs "," T0obs "," eT0obs "," T0_sim ",&
            &" T0obs-T0_sim "," T0_stat "
        do j1=1,nT0(j)
          write(uT0,trim(fmt))epoT0obs(j1,j),T0obs(j1,j),&
              &eT0obs(j1,j),T0_sim(j1,j),(T0obs(j1,j)-T0_sim(j1,j)),&
              &T0_stat(j1,j)
        end do
        flush(uT0)
        close(uT0)
        write(*,'(a,a)')" WRITTEN T0 INTO FILE: ",trim(flT0)
!         flush(6)
      end if
    end do

    return
  end subroutine write_T0_f
  ! ------------------------------------------------------------------ !

  ! ------------------------------------------------------------------ !
  ! adjust the parameters and writes them to screen
  subroutine write_par_s(cpuid,par,sigpar,resw)
    integer,intent(in)::cpuid
    real(dp),dimension(:),intent(in)::par,sigpar,resw
    real(dp),dimension(:),allocatable::par1,spar1
    integer::j
    character(80)::fmt

    allocate(par1(nfit),spar1(nfit))
    par1=par
    spar1=sigpar
    call param_adj(par1,spar1) ! adjusting parameters ... ex. angle with mod(angle,360) ...
    write(*,'(a,i4,a)')" CPU ",cpuid," ADJUSTED PARAMETERS "
    write(*,'(a,6x,a,11x,a,9x,a)')" PARAMETER ",&
        &" value "," +/- "," sigma "
    fmt=adjustl("(a7,1x,a,"//trim(sprec)//",1x,a,1x,"//trim(sprec)//")")
    do j=1,nfit
      write(*,trim(fmt))trim(parid(j))," = ",par1(j)," +/- ",spar1(j)
    end do
    write(*,*)""
    fmt=adjustl("(2(a,"//trim(sprec)//"),a,i6)")
    write(*,trim(fmt))" Fitness(Chi2r*k_chi2r + Chi2wr*k_chi2wr) = ",&
      &sum(resw*resw)," => Fitness*dof = ",(sum(resw*resw)*real(dof,dp)),&
      &" dof = ",dof
    write(*,*)""

    return
  end subroutine write_par_s

  ! adjust the parameters and writes them into file
  subroutine write_par_f(cpuid,isim,wrtid,par,sigpar,resw)
    integer,intent(in)::cpuid,isim,wrtid
    real(dp),dimension(:),intent(in)::par,sigpar,resw
    real(dp),dimension(:),allocatable::par1,spar1
    integer::upar,j
    character(512)::flpar
    character(80)::fmt

    allocate(par1(nfit),spar1(nfit))
    par1=par
    spar1=sigpar
    call param_adj(par1,spar1) ! adjusting parameters ... ex. angle with mod(angle,360) ...
    upar=get_unit(cpuid)
    flpar=""
    flpar=trim(path)//trim(adjustl(string(isim)))//"_"//&
      &trim(adjustl(string(wrtid)))//"_final"//&
      &trim(adjustl(string(nfit)))//"par.dat"
    flpar=trim(adjustl(flpar))
    open(upar,file=trim(flpar))
    write(*,'(a,i4,a,a)')" IN write_par_f => upar = ",upar," flpar = ",&
        &trim(flpar)
    write(upar,'(a,i3a)')"# CPU ",cpuid," ADJUSTED PARAMETERS "
    fmt=adjustl("(2(a,"//trim(sprec)//"),a,i6)")
    write(upar,trim(fmt))"# Fitness(Chi2r*k_chi2r + Chi2wr*k_chi2wr) = ",&
      &sum(resw*resw),&
      &" => Fitness*dof = ",(sum(resw*resw)*real(dof,dp)),&
      &" dof = ",dof
    write(upar,'(a,5x,a,24x,a)')"# parameter "," value "," sigma "
    fmt=adjustl("(a7,4x,"//trim(sprec)//",6x,"//trim(sprec)//")")
    do j=1,nfit
      write(upar,trim(fmt))trim(parid(j)),par1(j),spar1(j)
    end do
    flush(upar)
    close(upar)
    deallocate(par1,spar1)

    return
  end subroutine write_par_f
  ! ------------------------------------------------------------------ !

  ! ---
  ! adapted subroutine to write only final parameters
  subroutine write_parameters(cpuid,isim,wrtid,par,resw,fit_scale,gls_scale)
    integer,intent(in)::cpuid,isim,wrtid
    real(dp),dimension(:),intent(in)::par
    real(dp),dimension(:),intent(in)::resw
    real(dp),optional,intent(in)::fit_scale,gls_scale
    integer::upar,j
    character(512)::flpar
    character(80)::fmt
    real(dp)::fitness,bic,chi2,chi2r

    fitness = sum(resw*resw)
    if(present(fit_scale))then
      chi2r=fitness/((fit_scale*fit_scale)*(gls_scale*gls_scale))
    else
      chi2r=fitness
    end if
    chi2=chi2r*real(dof,dp)
!     bic=fitness*real(dof,dp) + real(nfit,dp)*log(real(ndata,dp))
    bic=chi2+real(nfit,dp)*log(real(ndata,dp))
    
    upar=get_unit(cpuid)
    flpar=""
    flpar=trim(path)//trim(adjustl(string(isim)))//"_"//&
      &trim(adjustl(string(wrtid)))//"_final"//&
      &trim(adjustl(string(nfit)))//"par.dat"
    flpar=trim(adjustl(flpar))
    open(upar,file=trim(flpar))
!     write(*,'(a,i4,a,a)')" IN write_par_f => upar = ",upar," flpar = ",&
!         &trim(flpar)
    write(upar,'(a,i3a)')"# CPU ",cpuid," PARAMETERS "
    
    fmt=adjustl("(2(a,"//trim(sprec)//"),a,i6)")
    if(present(fit_scale))then
      write(*,trim(fmt))"# Fitness(Chi2r*k_chi2r + Chi2wr*k_chi2wr) x fit_scale^2 x gls_scale^2 = ",&
        &fitness,&
        &" => Fitness_x_dof = ",(fitness*real(dof,dp)),&
        &" dof = ",dof
      write(upar,trim(fmt))"# Fitness(Chi2r*k_chi2r + Chi2wr*k_chi2wr) x fit_scale^2 x gls_scale^2 = ",&
        &fitness,&
        &" => Fitness_x_dof = ",(fitness*real(dof,dp)),&
        &" dof = ",dof
        
      write(*,'(2(a,es23.16))')"# fit_scale = ",fit_scale," gls_scale = ",gls_scale
      write(upar,'(2(a,es23.16))')"# fit_scale = ",fit_scale," gls_scale = ",gls_scale
    
    else
      write(*,trim(fmt))"# Fitness(Chi2r*k_chi2r + Chi2wr*k_chi2wr) = ",&
        &fitness,&
        &" => Fitness_x_dof = ",(fitness*real(dof,dp)),&
        &" dof = ",dof
      write(upar,trim(fmt))"# Fitness(Chi2r*k_chi2r + Chi2wr*k_chi2wr) = ",&
        &fitness,&
        &" => Fitness_x_dof = ",(fitness*real(dof,dp)),&
        &" dof = ",dof
    end if
    
    write(*,'(a,es23.16,a,i4,a,i5,a,es23.16)')"# BIC = Chi2 + nfit x ln(ndata) = ",chi2," + ",nfit," x ln(",ndata,") = ",bic
    write(upar,'(a,es23.16,a,i4,a,i5,a,es23.16)')"# BIC = Chi2 + nfit x ln(ndata) = ",chi2," + ",nfit," x ln(",ndata,") = ",bic
    
    write(*,'(a,es23.16)')"# LogLikelihood = - (dof/2)*ln(2pi) - (sum(ln(sigma**2)) - Fitness_x_dof / 2 = ", ln_err_const-half*fitness*real(dof,dp)
    write(upar,'(a,es23.16)')"# LogLikelihood = - (dof/2)*ln(2pi) - (sum(ln(sigma**2)) - Fitness_x_dof / 2 = ", ln_err_const-half*fitness*real(dof,dp)
    
    write(*,'(a,5x,a)')"# parameter "," value "
    write(upar,'(a,5x,a)')"# parameter "," value "
    
    fmt=adjustl("(a10,4x,"//trim(sprec)//")")
    do j=1,nfit
      write(*,trim(fmt))trim(parid(j)),par(j)
      write(upar,trim(fmt))trim(parid(j)),par(j)
    end do
    flush(6)
    flush(upar)
    close(upar)

    return
  end subroutine write_parameters
  ! ---
  
  ! write Fitness/Chi2r/Chi2wr to screen and file
  subroutine write_fitness_summary(cpuid,isim,wrtid,chi2r_RV,chi2wr_RV,chi2r_T0,chi2wr_T0,chi2r_oc,fitness,fit_scale,gls_scale)
    integer,intent(in)::cpuid,isim,wrtid
    real(dp),intent(in)::chi2r_RV,chi2r_T0,chi2wr_RV,chi2wr_T0,chi2r_oc,fitness
    real(dp),optional,intent(in)::fit_scale,gls_scale
    real(dp)::chi2wr_oc,chi2r,chi2wr,bic
    character(512)::summary_file
    integer::uwrt
    
    if(sum(nT0).gt.0)then
      chi2wr_oc=chi2r_oc*real(ndata,dp)/real(sum(nT0),dp)
      
      if(oc_fit.eq.2)then
        chi2r=chi2r_RV+half*(chi2r_T0+chi2r_oc)
        chi2wr=chi2wr_RV+half*(chi2wr_T0+chi2wr_oc)
      else if(oc_fit.eq.1)then
        chi2r=chi2r_RV+chi2r_oc
        chi2wr=chi2wr_RV+chi2wr_oc
      else
        chi2r=chi2r_RV+chi2r_T0
        chi2wr=chi2wr_RV+chi2wr_T0
      end if
    else
      chi2wr_oc=zero
      chi2r=chi2r_RV
      chi2wr=chi2wr_RV
    end if
    
    bic=chi2r*real(dof,dp) + real(nfit,dp)*log(real(ndata,dp))
    
    summary_file=trim(path)//trim(adjustl(string(isim)))//"_"//&
      &trim(adjustl(string(wrtid)))//"_fitness_summary.log"
    summary_file=trim(adjustl(summary_file))
    uwrt=get_unit(cpuid)
    open(uwrt,file=trim(summary_file))
    
    write(*,'(a,i2)')"# FITNESS SUMMARY: oc_fit == ",oc_fit
    write(uwrt,'(a,i2)')"# FITNESS SUMMARY: oc_fit == ",oc_fit
    
    if(ndata.gt.0)then
    
      ! print to screen
      write(*,'(3(a,i4))')" dof = ndata - nfit = ",ndata," - ",nfit," = ",dof
      
      if(present(fit_scale))then
        write(*,'(100(a,es23.16))')&
          &" Fitness (Chi2r*k_chi2r + Chi2wr*k_chi2wr) x fit_scale^2 x gls_scale^2 = ",fitness
        write(*,'(2(a,es23.16))')' fit_scale = ',fit_scale,' gls_scale = ',gls_scale
      else
        write(*,'(100(a,es23.16))')&
          &" Fitness (Chi2r*k_chi2r + Chi2wr*k_chi2wr) = ",fitness
      end if
      
      write(*,'(a,es23.16)')' Chi2 = ',chi2r*real(dof,dp)
      write(*,'(a,es23.16,a,i4,a,i5,a,es23.16)')" BIC = Chi2 + nfit x ln(ndata) = ",chi2r*real(dof,dp)," + ",nfit," x ln(",ndata,") = ",bic
      write(*,'(a,es23.16)')" LogLikelihood = - (dof/2)*ln(2pi) - (sum(ln(sigma**2))/2 - Fitness_x_dof/2 = ", ln_err_const-half*fitness*real(dof,dp)
      write(*,'(a)')''
      write(*,'(a,es23.16)')" k_chi2r   = ",k_chi2r
      write(*,'(a,es23.16)')" Chi2r     = ",chi2r
      write(*,'(a,es23.16)')" Chi2r_RV  = ",chi2r_RV
      write(*,'(a,es23.16)')" Chi2r_T0  = ",chi2r_T0
      write(*,'(a,es23.16)')" Chi2r_oc  = ",chi2r_oc
      write(*,'(a)')''
      write(*,'(a,es23.16)')" k_chi2wr  = ",k_chi2wr
      write(*,'(a,es23.16)')" Chi2wr    = ",chi2wr
      write(*,'(a,es23.16)')" Chi2wr_RV = ",chi2wr_RV
      write(*,'(a,es23.16)')" Chi2wr_T0 = ",chi2wr_T0
      write(*,'(a,es23.16)')" Chi2wr_oc = ",chi2wr_oc
      
      ! write into file
      write(uwrt,'(3(a,i4))')" dof = ndata - nfit = ",ndata," - ",nfit," = ",dof
      if(present(fit_scale))then
        write(uwrt,'(100(a,es23.16))')&
          &" Fitness (Chi2r*k_chi2r + Chi2wr*k_chi2wr) x fit_scale^2 x gls_scale^2 = ",fitness
        write(uwrt,'(2(a,es23.16))')' fit_scale = ',fit_scale,' gls_scale = ',gls_scale
        write(uwrt,'(a,es23.16)')' Chi2 = ',chi2r*real(dof,dp)
      else
        write(uwrt,'(100(a,es23.16))')&
          &" Fitness (Chi2r*k_chi2r + Chi2wr*k_chi2wr) = ",fitness
      end if
      
      write(uwrt,'(a,es23.16,a,i4,a,i5,a,es23.16)')" BIC = Chi2 + nfit x ln(ndata) = ",chi2r*real(dof,dp)," + ",nfit," x ln(",ndata,") = ",bic
      write(uwrt,'(a,es23.16)')" LogLikelihood = - (dof/2)*ln(2pi) - (sum(ln(sigma**2))/2 - Fitness_x_dof/2 = ", ln_err_const-half*fitness*real(dof,dp)
      write(uwrt,'(a)')''
      write(uwrt,'(a,es23.16)')" k_chi2r   = ",k_chi2r
      write(uwrt,'(a,es23.16)')" Chi2r     = ",chi2r
      write(uwrt,'(a,es23.16)')" Chi2r_RV  = ",chi2r_RV
      write(uwrt,'(a,es23.16)')" Chi2r_T0  = ",chi2r_T0
      write(uwrt,'(a,es23.16)')" Chi2r_oc  = ",chi2r_oc
      write(uwrt,'(a)')''
      write(uwrt,'(a,es23.16)')" k_chi2wr  = ",k_chi2wr
      write(uwrt,'(a,es23.16)')" Chi2wr    = ",chi2wr
      write(uwrt,'(a,es23.16)')" Chi2wr_RV = ",chi2wr_RV
      write(uwrt,'(a,es23.16)')" Chi2wr_T0 = ",chi2wr_T0
      write(uwrt,'(a,es23.16)')" Chi2wr_oc = ",chi2wr_oc
    else
      write(*,'(a)')" NOT ENOUGH DATA"
      write(uwrt,'(a)')" NOT ENOUGH DATA"
    end if
    
    flush(uwrt)
    close(uwrt)
    flush(6)
  
    return
  end subroutine write_fitness_summary

end module output_files
