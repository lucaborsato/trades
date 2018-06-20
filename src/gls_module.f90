! Luca conversion of the GLS (Zechmeister & Kuerster, 2009, A&A, 496, 577).
!
module gls_module
  use constants
  use parameters,only:NB,path
  use init_trades,only:get_unit
  use statistics
  use celestial_mechanics,only:calculate_true_anomaly ! ready to be used with TRADES
  implicit none

!   real(dp),parameter::delta_per=0.5_dp
  real(dp),parameter::delta_per=half
  
contains

  !     fit sine wave y=A_coeff*cosx+B_coeff*sinx+offset
  !     A_coeff,B_coeff,offset - fit parameter
  subroutine sine_fit(phase, data, weights, yy, A_coeff, B_coeff, offset, power_GLS, power_LS)
    real(dp),dimension(:),intent(in)::phase,data,weights
    real(dp),intent(in)::yy
    real(dp)::A_coeff,B_coeff,offset,power_GLS,power_LS

    integer::n_data
    real(dp),dimension(:),allocatable::cosx,sinx
    real(dp):: cc, ss, cs, c, s, yc, ys, d
    ! in GLS.f SineFit: COMMON /data/ v,z,ww,YY,N  ! v = phase, z = data, ww = weights
    ! in GLS.f main:    COMMON /data/ v,wy,ww,YY,N
    real(dp)::inv_yyd,inv_yy,inv_d
    
    n_data=size(data)
    allocate(cosx(n_data), sinx(n_data))
    
    cosx=cos(phase)
    sinx=sin(phase)
    cc=sum(weights*cosx*cosx)
    cs=sum(weights*cosx*sinx)
    c=sum(weights*cosx)
    s=sum(weights*sinx)
    yc=sum(data*cosx)
    ys=sum(data*sinx)
    
    ss=one-cc ! one in 'constants' module
    d=(cc*ss)-(cs*cs)
    
    inv_d=one/d
    inv_yy=one/yy
    inv_yyd=inv_yy*inv_d
    
  !   powLS = (SS*YC**2/D+CC*YS**2/D-2*CS*YC*YS/D) / YY ! original
  
!     power_LS=( ((ss*yc*yc)/d) + ((cc*ys*ys)/d) - ((two*cs*yc*ys)/d) ) / yy         ! Lomb-Scargle power

    power_LS= inv_yyd * ((ss*yc*yc) + (cc*ys*ys) - (two*cs*yc*ys))
    ! two in 'constants' module
!     power_LS= inv_yy * (inv_d*(ss*yc*yc) + inv_d*(cc*ys*ys) - inv_d*(two*cs*yc*ys))
    
    
    cc = cc - (c*c)
    ss = ss - (s*s)
    cs = cs - (c*s)
    d  = (cc*ss) - (cs*cs)
    inv_d=one/d
!     A_coeff = (yc*ss-ys*cs) / d
!     B_coeff = (ys*cc-yc*cs) / d
    A_coeff = ((yc*ss)-(ys*cs)) * inv_d
    B_coeff = ((ys*cc)-(yc*cs)) * inv_d
    offset = -(A_coeff*c)-(B_coeff*s)
    inv_yy=one/yy
    inv_yyd=inv_yy*inv_d
    
  !   pow= (A*YC+B*YS)/ YY
  !   pow = (SS*YC**2/D+CC*YS**2/D-2*CS*YC*YS/D) / YY
!     power_GLS = ( ((ss*yc*yc)/d)+((cc*ys*ys)/d)-((two*cs*yc*ys)/d) ) / yy           ! GLS power

    power_GLS = inv_yyd * ((ss*yc*yc) + (cc*ys*ys) - (two*cs*yc*ys))
!     power_GLS = inv_yy * (inv_d*(ss*yc*yc) + inv_d*(cc*ys*ys) - inv_d*(two*cs*yc*ys))
    
    deallocate(cosx,sinx)

    return
  end subroutine

  subroutine spectral_window(phase, n_data, ws_out) ! ws_out -> ws in GLS.f, it should mean window spectral ... maybe...
    real(dp),dimension(:),intent(in)::phase
    integer,intent(in)::n_data
    
    real(dp),intent(out)::ws_out
    
    real(dp)::wc,ws
    
    wc=sum(cos(phase))
    ws=sum(sin(phase))
    ws_out= (wc*wc + ws*ws) / (n_data*n_data)
    
    return
  end subroutine spectral_window
  
!   subroutine calculate_mass(mass_1,K,period,eccentricity,err_fac,mass_2,a1sini,a_axis)
!     real(dp),intent(in)::mass_1,K,period,eccentricity,err_fac
!     ! twopi,G,AU,Msun,Mjup in GLS.f from 'constants' module
!     real(dp),intent(out)::mass_2,a1sini,a_axis
!     
!     real(dp)::mass_func
!     
!     real(dp),parameter::sec_dpi=s24h/dpi !in 'constants' module
!     real(dp)::ecc_2
!     
!     ecc_2=eccentricity*eccentricity
!     mass_2=zero
!     if(eccentricity.le.TOLERANCE)then
!       write(*,'(a)')'Companion minimum mass Sine fit'
!     else
!       write(*,'(a)')'Companion minimum mass Kepler fit'
!     end if
!   
!     mass_func=sec_dpi/Gsi/Msun*period*(abs(K)*sqrt(one-ecc_2))**3
!     write(*,'(a,es23.16)')" mass function [Msun]:", mass_func
!   
! !     iteration to determine the companion mass m2
!     do i=1,10
!         mass_2 = (((mass_1+mass_2)**2)*mass_func)**onethird
!     end do
!     a1sini = sec_dpi*K*period*sqrt(one-ecc_2)
!     a_axis = ((mass_1+mass_2)/mass_2) * a1sini
!     
!     write(*,'(2(a,es23.16))')" Mmin  [Mjup]:  ", mass_2*Msunj,&
!       &" +/-", err_fac*mass_2*Msunj/3._dp
!     write(*,'(2(a,es23.16))')" Mmin  [Msun]:  ", mass_2,&
!       &" +/-", err_fac*mass_2/3._dp
!     write(*,'(2(a,es23.16))')" a1*sini [AU]:  ", a1sini/AU
!       &" +/-", err_fac*a1sini/AU
!     write(*,'(2(a,es23.16))')" a       [AU]:  ", a_axis/AU
!       &" +/-", err_fac*a_axis/AU
!     
!     mass_2 = mass_2*Msunj
! 
!     return
!   end subroutine calculate_mass
    
  subroutine phaser(time_jd,dpi_freq,phase)
    real(dp),dimension(:),intent(in)::time_jd
    real(dp),intent(in)::dpi_freq
    real(dp),dimension(:),intent(out)::phase
    
    phase=time_jd*dpi_freq
  
    return
  end subroutine phaser

  
  subroutine period_min_for_gls(jd,p_min)
    real(dp),dimension(:),intent(in)::jd
    real(dp),intent(out)::p_min
    
    real(dp),dimension(:),allocatable::djd
    integer::n_jd,i_jd
    
    p_min=zero
    n_jd=size(jd)
    allocate(djd(n_jd-1))
    do i_jd=2,n_jd
      djd(i_jd-1)=abs(jd(i_jd)-jd(i_jd-1))
    end do
    p_min=djd(nint(half*(n_jd-1)))
    deallocate(djd)
  
    return
  end subroutine period_min_for_gls
  
  subroutine gls(jd,rv,erv,freq,power_GLS,sp_window,power_LS)
    real(dp),dimension(:),intent(in)::jd,rv,erv
    real(dp),dimension(:),allocatable,intent(out)::freq,power_GLS,sp_window,power_LS
    
    integer::n_data
    real(dp),dimension(:),allocatable::time_jd,ww,wy,phase,res_sin
    real(dp),dimension(:),allocatable::inv_erv2
    real(dp)::wsum,time_base,inv_time_base,jd_min,rv_mean,yy

    real(dp)::A_coeff,B_coeff,C_offset
    
    real(dp)::per_min
    real(dp)::freq_start,freq_end,freq_step
    integer::n_freq,i_freq,freq_oversampling    
    
    n_data=size(jd)
    
    allocate(inv_erv2(n_data),ww(n_data),phase(n_data),res_sin(n_data))
    inv_erv2=one/(erv*erv)
    ww=inv_erv2 ! weights = 1/err^2
    
    wsum=sum(ww) ! sum of 1/err^2
    ww=ww/wsum ! normalized weights
    
    ! shift time series
    allocate(time_jd(n_data))
    jd_min=minval(jd)
    time_jd=jd-jd_min
    time_base=abs(maxval(jd)-jd_min)
    inv_time_base=one/time_base
    
    ! shift RV by mean and prepare Chi2_0 == YY
    rv_mean=sum(rv*ww) ! weighted mean with normalized weights
    wy=rv-rv_mean ! shifted RV by RV_mean
    yy=sum(wy*wy*ww)
    wy=wy*ww
    
    ! frequencies
    freq_start=one/time_base
    call period_min_for_gls(time_jd,per_min)
    freq_end=one/per_min
    freq_oversampling=20
    freq_step=inv_time_base/real(freq_oversampling,dp)        ! frequency step
    n_freq=int( ((freq_end-freq_start)/freq_step) + 1) ! number of steps
    
    allocate(freq(n_freq),power_GLS(n_freq),sp_window(n_freq),power_LS(n_freq))
    
    do i_freq=1,n_freq
      freq(i_freq)=freq_start+freq_step*real((i_freq-1),dp)
      phase=zero !reset phase
      call phaser(time_jd,dpi*freq(i_freq),phase)
      call sine_fit(phase,wy,ww,yy,&
        &A_coeff,B_coeff,C_offset,&
        &power_GLS(i_freq),power_LS(i_freq))
      call spectral_window(wy,n_data,sp_window(i_freq))
    end do
    
    deallocate(time_jd,inv_erv2,ww,phase,res_sin)
    
    return
  end subroutine gls
  
  subroutine calculates_power_max(periods_bounds,freq,power_GLS,power_max_within,power_max_outside,period_max,pl_max)
    real(dp),dimension(:,:),intent(in)::periods_bounds
    real(dp),dimension(:),intent(in)::freq,power_GLS
    real(dp),intent(out)::power_max_within,power_max_outside
    real(dp),intent(out)::period_max
    integer,intent(out)::pl_max
    
    logical,dimension(:),allocatable::power_stat

    integer::n_freq,i_nb,i_freq
    
    n_freq=size(freq)
    allocate(power_stat(n_freq))
    power_stat=.false.
    loopf: do i_freq=1,n_freq
      period_max=one/freq(i_freq)
      loopnb: do i_nb=2,NB
        if(period_max.ge.periods_bounds(1,i_nb))then
          if(period_max.le.periods_bounds(2,i_nb))then
            power_stat(i_freq)=.true.
            pl_max=i_nb
            exit loopnb
          end if
        end if
      end do loopnb
    end do loopf
!     power_max_within=maxval(power_GLS(power_stat)) ! not working
!     power_max_outside=maxval(power_GLS(.not.power_stat)) ! not working
    power_max_within=maxval(pack(power_GLS,power_stat))
    power_max_outside=maxval(pack(power_GLS,.not.power_stat))
    deallocate(power_stat)
    
    return
  end subroutine calculates_power_max
  
  subroutine check_periodogram(jd,rv,erv,periods,gls_check)
    real(dp),dimension(:),intent(in)::jd,rv,erv
    real(dp),dimension(:),intent(in)::periods
    logical,intent(out)::gls_check
    

    real(dp),dimension(:,:),allocatable::periods_bounds
    real(dp),dimension(:),allocatable::freq,power_GLS,sp_window,power_LS
    
    real(dp)::power_max_within,power_max_outside,power_threshold,period_max
    integer::pl_max
    
    gls_check=.true. ! if gls_check the fit is ok, it means no induced signals close to planetary periods
    
    call gls(jd,rv,erv,freq,power_GLS,sp_window,power_LS) ! run gls
    
    ! define boundaries of the periods to check, remember that the first 'element' is the star and it will not be taken into account
    allocate(periods_bounds(2,NB))
    periods_bounds(1,2:NB)=periods(2:NB)-delta_per
    periods_bounds(2,2:NB)=periods(2:NB)+delta_per
    call calculates_power_max(periods_bounds,freq,power_GLS,power_max_within,power_max_outside,period_max,pl_max)
    deallocate(periods_bounds)
    deallocate(freq,power_GLS,sp_window,power_LS)
  
    power_threshold=two*power_max_outside
    if(power_max_within.lt.power_threshold)then ! max peak with periods boudaries < 2*max peak outside --> GOOD
      gls_check=.true.
    else ! max peak with periods boudaries >= 2*max peak outside --> BAD
      gls_check=.false.
    end if
  
    return
  end subroutine check_periodogram

  subroutine check_periodogram_scale(jd,rv,erv,periods,gls_check,gls_scale)
    real(dp),dimension(:),intent(in)::jd,rv,erv
    real(dp),dimension(:),intent(in)::periods
    logical,intent(out)::gls_check
    real(dp),intent(out)::gls_scale
    
    real(dp),dimension(:,:),allocatable::periods_bounds
    real(dp),dimension(:),allocatable::freq,power_GLS,sp_window,power_LS
    
    real(dp)::power_max_within,power_max_outside,power_threshold,period_max,delta_max
    integer::pl_max

    gls_check=.true. ! if gls_check the fit is ok, it means no induced signals close to planetary periods
    gls_scale=zero
    
    call gls(jd,rv,erv,freq,power_GLS,sp_window,power_LS) ! run gls
    
    ! define boundaries of the periods to check, remember that the first 'element' is the star and it will not be taken into account
    allocate(periods_bounds(2,NB))
    periods_bounds(1,2:NB)=periods(2:NB)-delta_per
    periods_bounds(2,2:NB)=periods(2:NB)+delta_per
    call calculates_power_max(periods_bounds,freq,power_GLS,power_max_within,power_max_outside,period_max,pl_max)
    deallocate(freq,power_GLS,sp_window,power_LS)
  
    power_threshold=two*power_max_outside
    if(power_max_within.lt.power_threshold)then ! max peak with periods boudaries < 2*max peak outside --> GOOD
      gls_check=.true.
      gls_scale=one
    else ! max peak with periods boudaries >= 2*max peak outside --> BAD
      gls_check=.false.
      delta_max=max(abs(period_max-periods_bounds(1,pl_max)),abs(periods_bounds(2,pl_max)-period_max))
      if(delta_max.lt.TOL_dp) delta_max=TOL_dp
      gls_scale=one-log10(delta_max)
    end if
    deallocate(periods_bounds)
  
    return
  end subroutine check_periodogram_scale

    
  subroutine check_and_write_periodogram(cpuid,id_sim,wrt_id,jd,rv,erv,periods,gls_check,gls_scale)
    use convert_type,only:string
    integer,intent(in)::cpuid,id_sim,wrt_id
    real(dp),dimension(:),intent(in)::jd,rv,erv
    real(dp),dimension(:),intent(in)::periods
    logical,intent(out)::gls_check
    real(dp),optional,intent(out)::gls_scale
    
    real(dp),dimension(:,:),allocatable::periods_bounds
    real(dp),dimension(:),allocatable::freq,power_GLS,sp_window,power_LS
    
    real(dp)::power_max_within,power_max_outside,power_threshold,period_max,delta_max
    integer::pl_max
    character(512)::gls_file
    integer::i_freq,n_freq,u_gls
    
    gls_check=.true. ! if gls_check the fit is ok, it means no induced signals close to planetary periods
    
    call gls(jd,rv,erv,freq,power_GLS,sp_window,power_LS) ! run gls
    
    ! define boundaries of the periods to check, remember that the first 'element' is the star and it will not be taken into account
    allocate(periods_bounds(2,NB))
    periods_bounds(1,2:NB)=periods(2:NB)-delta_per
    periods_bounds(2,2:NB)=periods(2:NB)+delta_per
    call calculates_power_max(periods_bounds,freq,power_GLS,power_max_within,power_max_outside,period_max,pl_max)
    
    ! write file gls
    n_freq=size(freq)
    gls_file=trim(path)//trim(adjustl(string(id_sim)))//"_"//&
      &trim(adjustl(string(wrt_id)))//'_gls_output.dat'
    u_gls=get_unit(cpuid)
    open(u_gls,file=trim(gls_file))
    write(u_gls,'(a)')'# Period Frequency power_GLS sp_window power_LS'
    do i_freq=1,n_freq
      write(u_gls,'(1000(es23.16,1x))')one/freq(i_freq),freq(i_freq),power_GLS(i_freq),sp_window(i_freq),power_LS(i_freq)
    end do
    close(u_gls)
    deallocate(freq,power_GLS,sp_window,power_LS)
  
    power_threshold=two*power_max_outside
    if(power_max_within.lt.power_threshold)then ! max peak with periods boudaries < 2*max peak outside --> GOOD
      gls_check=.true.
      if(present(gls_scale))gls_scale=one
    else ! max peak with periods boudaries >= 2*max peak outside --> BAD
      gls_check=.false.
      if(present(gls_scale))then
        delta_max=max(abs(period_max-periods_bounds(1,pl_max)),abs(periods_bounds(2,pl_max)-period_max))
        if(delta_max.lt.TOL_dp) delta_max=TOL_dp
        gls_scale=one-log10(delta_max)
      end if
    end if
    deallocate(periods_bounds)
  
    return
  end subroutine check_and_write_periodogram
  
  
end module gls_module
