! module with all needed subroutines and function to initialize trades

module init_trades
  use constants
  use parameters
  use parameters_conversion
  use convert_type,only:string
  implicit none

!   interface init_random_seed
!     module procedure init_random_seed_one,init_random_seed_two
!   end interface init_random_seed
  
  contains

  ! function to check if a unit file is used or not
  function get_unit(cpu) result(newunit)
    integer::newunit
    integer,intent(in)::cpu
    integer,parameter::umin=1,umax=90
    logical::status
    integer::set
    set=umin
    uloop:do
      if(set.gt.90) set=umin
      inquire(unit=unit2(set,cpu),opened=status)
      if(.not.status)then
        newunit=unit2(set,cpu)
        exit uloop
      end if
      set=set+1
    end do uloop

    return
  end function get_unit

  ! ------------------------- !
  ! function to read a unit file and count the rows
  ! and rewind the file
  function get_rows(unit) result(nrows)
    integer::nrows
    integer,intent(in)::unit
    integer::stat
    character(128)::comment

    nrows=0
    count:do
      read(unit,*,IOSTAT=stat)comment
      comment=trim(adjustl(comment))
      if(IS_IOSTAT_END(stat)) exit count
      if(comment(1:1).ne."#") nrows=nrows+1
    end do count
    rewind(unit)

    return
  end function get_rows
  
  ! ------------------------------------------------------------------ !
  ! initialize the units needed to write files...even with openMP
  subroutine initu(nf,nc)
    integer,intent(in)::nf,nc
    integer::kc,kc1
    integer::kf
    integer::set
    if(allocated(unit2)) deallocate(unit2)
    allocate(unit2(nf,nc))
    unit2=0
    do kc=1,nc
      kc1=nf*(kc-1)
      do kf=1,nf
        set=10+kf+kc1
        unit2(kf,kc)=set
      end do
    end do

    return
  end subroutine initu
  ! ------------------------------------------------------------------ !

  ! ------------------------------------------------------------------ !
  ! writes an execution help
  subroutine help(arg)
    character(*),intent(in)::arg

    if(  (arg(1:2).eq."-h")&
        &.or.&
        &(arg(1:3).eq."--h")&
        &.or.&
        &(arg(1:5).eq."-help")&
        &.or.&
        &(arg(1:6).eq."--help"))then
      write(*,'(a)')""
      write(*,'(a)')" HELP: main_pik"
      write(*,'(a)')" a) execute: main_pik path/"
      write(*,'(a)')" with a / at the end"
      write(*,'(a)')" or"
      write(*,'(a)')" b) execute: main_pik"
      write(*,'(a)')" it will be used the path ./"
      write(*,'(a)')""
      stop
    end if

    return
  end subroutine help

  ! IT READS THE COMMAND ARGUMENT THAT DEFINE THE PATH OF THE FOLDER WITH THE FILES
  subroutine read_com_arg(arg1)
    character(*)::arg1
    integer::nargs,ilast
    character::last

    last=""
    nargs=command_argument_count()
    if(nargs.eq.0)then
      arg1="./"
    else
      call get_command_argument(1,arg1)
      arg1=trim(adjustl(arg1))
      ilast=len_trim(arg1)
      last=trim(arg1(ilast:ilast))
      if(last.ne."/") arg1=trim(arg1)//"/"
      arg1=trim(adjustl(arg1))
      ! HELP
      call help(arg1)
    end if

    return
  end subroutine read_com_arg
  ! ------------------------------------------------------------------ !

  ! ------------------------------------------------------------------ !
  ! initialize the vector for the selection of the parameter to be fitted
  subroutine tofit_init()

    npar=2+(NB-1)*8
    if(.not.allocated(tofit)) allocate(tofit(npar))

    return
  end subroutine tofit_init
  ! ------------------------------------------------------------------ !

  
  ! ------------------------------------------------------------------ !
  ! IT READS ARGUMENT/OPTION FOR THE PROGRAM FROM arg.in FILE
  subroutine read_arg(cpuid)
    integer,intent(in)::cpuid
    integer::uread
    logical::fstat
    character(512)::line
    integer::idx,idh,istat,it
    character(1)::hill_temp
    character(1),dimension(4)::list_true=(/'y', 'Y', 't', 'T'/)
    
    
    inquire(file=trim(path)//"arg.in",exist=fstat)
    if(fstat)then
      uread=get_unit(cpuid)
      open(uread,file=trim(path)//"arg.in",status='OLD')
      reading:do
        read(uread,'(a512)',IOSTAT=istat) line
        if(IS_IOSTAT_END(istat)) exit reading
        line = trim(adjustl(line))
        if(len(trim(line)).ne.0)then
          if(line(1:1).ne."#")then
            idh=index(line,"#")
            if(idh.eq.0) idh=len(trim(line))+1
            idx=index(line,"=")
            
            if(idx.ne.0)then
              if (line(1:idx-1).eq.'progtype')then
                read(line(idx+1:idh),*) progtype
              else if (line(1:idx-1).eq.'nboot')then
                read(line(idx+1:idh),*) nboot
              else if (line(1:idx-1).eq.'NB')then
                read(line(idx+1:idh),*) NB
              else if (line(1:idx-1).eq.'idtra')then
                read(line(idx+1:idh),*) idtra
              else if (line(1:idx-1).eq.'durcheck')then
                read(line(idx+1:idh),*) durcheck
              else if (line(1:idx-1).eq.'wrtorb')then
                read(line(idx+1:idh),*) wrtorb
              else if (line(1:idx-1).eq.'wrtconst')then
                read(line(idx+1:idh),*) wrtconst
              else if (line(1:idx-1).eq.'wrtel')then
                read(line(idx+1:idh),*) wrtel
              else if (line(1:idx-1).eq.'rvcheck')then
                read(line(idx+1:idh),*) rvcheck
              else if (line(1:idx-1).eq.'idpert')then
                read(line(idx+1:idh),*) idpert
              else if (line(1:idx-1).eq.'lmon')then
                read(line(idx+1:idh),*) lmon
              else if (line(1:idx-1).eq.'tepoch')then
                read(line(idx+1:idh),*) tepoch
              else if (line(1:idx-1).eq.'tstart')then
                read(line(idx+1:idh),*) tstart
              else if (line(1:idx-1).eq.'tint')then
                read(line(idx+1:idh),*) tint
              else if (line(1:idx-1).eq.'step')then
                read(line(idx+1:idh),*) step_0
              else if (line(1:idx-1).eq.'wrttime')then
                read(line(idx+1:idh),*) wrttime
              else if (line(1:idx-1).eq.'tol_int')then
                read(line(idx+1:idh),*) tol_int
              else if(line(1:idx-1).eq.'bootstrap_scaling')then
                read(line(idx+1:idh),*) bootstrap_scaling
              else if(line(1:idx-1).eq.'weight_chi_square')then
                read(line(idx+1:idh),*) k_chi2r
                if(k_chi2r.gt.one.or.k_chi2r.lt.zero) k_chi2r = one
                k_chi2wr = one - k_chi2r
              else if(line(1:idx-1).eq.'secondary_parameters')then
                read(line(idx+1:idh),*) secondary_parameters
              else if(line(1:idx-1).eq.'do_hill_check')then
                read(line(idx+1:idh),*) hill_temp
                do it=1,4
                  if(hill_temp.eq.list_true(it))then
                    do_hill_check=.true.
                    exit
                  end if
                end do
              else if(line(1:idx-1).eq.'oc_fit')then
                read(line(idx+1:idh),*) oc_fit
              end if
              
            end if
            
          end if
        end if
      end do reading
      close(uread)
    else
        write(*,'(a,a,a)')" CANNOT FIND ARGUMENT FILE ",trim(path),"arg.in"
        write(*,'(a,a,a)')""
        stop
    end if

    return
  end subroutine read_arg

  !IT READS BODY NAMES AND DETERMINES WHICH PARAMETER TO FIT
  subroutine read_list(cpuid)
    integer,intent(in)::cpuid
    character(128)::temp
    character::fitpar
    integer::i,i1,i2,j,pos,ulst,stat,ih,it
    logical::fstat

    pos=0
    call tofit_init()
    inquire(file=trim(path)//"bodies.lst",exist=fstat)
    if(fstat)then
      if(.not.allocated(bnames)) allocate(bnames(NB),bfiles(NB),do_transit(NB))
      do_transit(1)=.false. ! star not transiting ... ...
      do_transit(2:NB)=.true. ! set all planet to transit by default
      
      ulst=get_unit(cpuid)
      open(ulst,file=trim(path)//'bodies.lst',status='OLD')
!       do i=1,NB
      
      i=0
      readbody:do
        read(ulst,'(a128)', IOSTAT=stat) temp
        if(IS_IOSTAT_END(stat)) exit readbody
        
        if(temp(1:1).ne."#")then
          i=i+1 ! row
          i1=index(temp,' ')-1
          i2=i1+2
          bfiles(i)=trim(adjustl(temp(1:i1)))
          bnames(i)=trim(adjustl(temp(1:i1-4)))
        
          if(i.eq.1)then ! row 1 == star
        
            do j=1,2
              fitpar=temp(i2:i2)
              read(fitpar,*)tofit(j)
              i2=i2+2
            end do
        
          else ! row from 2 to NB == planets
        
            do j=1,8 ! read fitting parameters: 1 fit, 0 no fit
              pos=i+j+(i-2)*7
              fitpar=temp(i2:i2)
              read(fitpar,*)tofit(pos)
              i2=i2+2
            end do
            
            ih = index(temp(1:len(trim(adjustl(temp)))),'#')
            if(ih.eq.0)ih=len(trim(adjustl(temp)))
            ! then last column: planet should transit? T = True/Yes (or omit), F = False/No. Before # character.
            it=scan(temp(1:ih),'Ff')
            if(it.ne.0)then
              read(temp(it:it),*)do_transit(i)
            end if
        
          end if
        
        end if
      
      end do readbody
      
      close(ulst)
      nfit=sum(tofit)
    
    else
    
      write(*,'(a,a,a)')" CANNOT FIND ARGUMENT FILE ",trim(path),"bodies.lst"
      write(*,'(a,a,a)')""
      stop
    
    end if

    return
  end subroutine read_list

  
  ! initialize variables for Levenberg-Marquardt
  subroutine read_lm_opt(cpuid)
    !$ use omp_lib
    integer,intent(in)::cpuid
    integer::uread,ncpu
    logical::fstat

    inquire(file=trim(path)//"lm.opt",exist=fstat)
    if(fstat)then
      uread=get_unit(cpuid)
      open(uread,file=trim(path)//"lm.opt",status='OLD')
      read(uread,*)maxfev
      read(uread,*)ftol
      read(uread,*)xtol
      read(uread,*)gtol
      read(uread,*)epsfcn
      read(uread,*)nprint
      close(uread)
      if(maxfev.le.0) maxfev = 500*(nfit + 1)
      if(ftol.le.zero) ftol = TOLERANCE
      if(xtol.le.zero) xtol = TOLERANCE
      if(gtol.le.zero) gtol = TOLERANCE
      if(nprint.lt.0) nprint = 0
      if(epsfcn.le.zero) epsfcn=TOLERANCE
      ncpu = 1
      !$omp parallel
      !$omp master
      !$ ncpu = omp_get_num_threads()
      if(.not.allocated(lmtols)) allocate(lmtols(ncpu,4))
      lmtols(:,1)=xtol
      lmtols(:,2)=ftol
      lmtols(:,3)=gtol
      lmtols(:,4)=epsfcn
      !$omp end master
      !$omp end parallel
    else
      write(*,'(a,a,a)')" CANNOT FIND ARGUMENT FILE ",trim(path),"lm.opt"
      write(*,'(a,a,a)')""
      stop
    end if

    return
  end subroutine read_lm_opt

  ! allocation and zero initialization of Keplerian elements
  subroutine init_zero_par(m,R,P,a,e,w,mA,inc,lN,tau)
    real(dp),dimension(:),allocatable,intent(out)::m,R,P,a,e,w,mA,inc,lN
    real(dp),dimension(:),allocatable,intent(out)::tau

    if(.not.allocated(m)) allocate(m(NB),R(NB),P(NB),a(NB),e(NB),&
        &w(NB),mA(NB),inc(NB),lN(NB),tau(NB))
    m=zero
    m(2:NB)=9999._dp ! high values
    R=zero
    R(2:NB)=9999._dp ! high values
    P=zero
    P(2:NB)=9.e7_dp ! high values
    a=zero
    a(2:NB)=999._dp ! high values
    e=zero
    e(2:NB)=999._dp ! high values
    w=zero
    w(2:NB)=999._dp ! high values
    mA=zero
    mA(2:NB)=999._dp ! high values
    inc=zero
    inc(2:NB)=999._dp ! high values
    tau=zero
    tau(2:NB)=9.e8_dp ! high values
    lN=zero
    lN(2:NB)=999._dp ! high values

    return
  end subroutine init_zero_par

  ! reads the orbital elements from the files
!   subroutine read_par(cpuid,m,R,P,a,e,w,mA,i,lN)
!     use celestial_mechanics,only:semax,period
!     integer,intent(in)::cpuid
!     real(dp),dimension(:),allocatable,intent(out)::m,R,P,a,e,w,mA,i,lN
!     real(dp),dimension(:),allocatable::tau
!     character::ttype
!     integer::unit,j
!     logical::fstat
!     real(dp)::tempm
! 
!     inquire(file=trim(path)//trim(bfiles(1)),exist=fstat)
!     if(fstat)then
!       call init_zero_par(m,R,P,a,e,w,mA,i,lN,tau)
!       !read and convert star mass from Msun to Mjup
!       unit=get_unit(cpuid)
!       open(unit,file=trim(path)//bfiles(1),status='OLD')
!       read(unit,*) m(1) ! Msun
!       read(unit,*) R(1) ! Rsun
!       close(unit)
!       do j=2,NB
!         unit=get_unit(cpuid)
!         open(unit,file=trim(path)//trim(bfiles(j)),status='OLD')
!         read(unit,*) m(j)
!         read(unit,*) R(j)
!         read(unit,*) P(j)
!         read(unit,*) a(j)
!         read(unit,*) e(j)
!         read(unit,*) w(j)
!         if(w(j).ge.360._dp) w(j)=mod(w(j),360._dp)
!         read(unit,*) tau(j),tempm,ttype
!         if(ttype.eq.'t')then
!           mA(j)=mod((360._dp/P(j))*(tepoch-tau(j)),360._dp)
!         else
!           mA(j)=tau(j)
!         end if
!         read(unit,*) i(j)
!         read(unit,*) lN(j)
!         close(unit)
!         m(j)=m(j)*Mjups
!         R(j)=R(j)*Rjups
!         if(a(j).ge.999._dp)then
!           a(j)=semax(m(1),m(j),P(j))
!         else if(P(j).ge.9.e6_dp)then
!           P(j)=period(m(1),m(j),a(j))
!         end if
!       end do
!     else
!       write(*,'(a,a,a)')" CANNOT FIND STAR FILE ",trim(path),trim(bfiles(1))
!       write(*,'(a)')""
!       stop
!     end if
!     deallocate(tau)
! 
!     return
!   end subroutine read_par

  ! updated subroutine to read parameters from x.dat files
  ! WARNING: added tau row
  subroutine read_par(cpuid,m,R,P,a,e,w,mA,i,lN)
    use celestial_mechanics,only:semax,period,tau2mA
    integer,intent(in)::cpuid
    real(dp),dimension(:),allocatable,intent(out)::m,R,P,a,e,w,mA,i,lN
    
    real(dp),dimension(:),allocatable::tau
    integer::unit,j
    logical::fstat
    character(512)::temp_line
!     real(dp)::temp1,temp2
    
    MR_star = zero ! init MR_star to zero value: in 'parameters' module
    
    inquire(file=trim(path)//trim(bfiles(1)),exist=fstat)
    if(fstat)then
      call init_zero_par(m,R,P,a,e,w,mA,i,lN,tau)
      !read and convert star mass from Msun to Mjup
      unit=get_unit(cpuid)
      open(unit,file=trim(path)//bfiles(1),status='OLD')
!       read(unit,*) m(1) ! Msun
!       read(unit,*) R(1) ! Rsun
      
      read(unit,*) MR_star(1,1),temp_line ! Msun eMsun
      temp_line=trim(adjustl(temp_line))
      if(temp_line(1:1).ne.'#') read(temp_line,'(es23.16)') MR_star(1,2)
      
      read(unit,*) MR_star(2,1),temp_line ! Rsun eRsun
      temp_line=trim(adjustl(temp_line))
      temp_line=trim(adjustl(temp_line))
      if(temp_line(1:1).ne.'#') read(temp_line,'(es23.16)') MR_star(2,2)
      
      close(unit)
      
      m(1)=MR_star(1,1)
      R(1)=MR_star(2,1)
      
      do j=2,NB
        
        unit=get_unit(cpuid)
        open(unit,file=trim(path)//trim(bfiles(j)),status='OLD')
        read(unit,*) m(j)
        read(unit,*) R(j)
        read(unit,*) P(j)
        read(unit,*) a(j)
        read(unit,*) e(j)
        read(unit,*) w(j)
        read(unit,*) mA(j)
        read(unit,*) tau(j)
        read(unit,*) i(j)
        read(unit,*) lN(j)
        close(unit)
        
        ! adjust and select properly the parameters
        m(j)=m(j)*Mjups
        R(j)=R(j)*Rjups
        
        if(a(j).ge.999._dp)then
          a(j)=semax(m(1),m(j),P(j))
        else if(P(j).ge.9.e6_dp)then
          P(j)=period(m(1),m(j),a(j))
        end if

        if(w(j).ge.360._dp) w(j)=mod(w(j),360._dp)
        
        ! if mA >= 999 check if it has to use the pericenter time tau
        if(mA(j).ge.999._dp)then
          
          ! check also if tau < 9.e8
          if(tau(j).lt.9.e8_dp)then
            ! from tau to mA
!             mA(j)=mod((360._dp/P(j))*(tepoch-tau(j)),360._dp)
            mA(j)=tau2mA(tau(j),tepoch,P(j))
            
          else
            deallocate(tau)
            write(*,'(a,a)')' WARNING: MISSING BOTH MEAN ANOMAMLY AND TIME OF PERICENTER FOR PLANET INPUT FILE ',trim(bfiles(j))
            stop
          end if
          
        end if
        
      end do
      
    else
      deallocate(tau)
      write(*,'(a,a,a)')" CANNOT FIND STAR FILE ",trim(path),trim(bfiles(1))
      write(*,'(a)')""
      stop
    end if
    deallocate(tau)

    return
  end subroutine read_par

  ! it reads boundaries of the Keplerian elements for PIKAIA simulation
  subroutine read_par_boundaries(cpuid,m,R)
    use celestial_mechanics,only:period,semax
    integer,intent(in)::cpuid
    real(dp),dimension(:),intent(in)::m,R
    real(dp)::temp1,temp2
    real(dp),dimension(:),allocatable::Pvec
    integer::upar,j,j1

    if(.not.allocated(par_min)) allocate(par_min(npar),par_max(npar))

    par_min=zero
    par_max=zero

    allocate(Pvec(NB-1))
    par_min(1)=MR_star(1,1) ! Mstar
    par_max(1)=MR_star(1,1) ! Mstar
    par_min(2)=MR_star(2,1) ! Rstar
    par_max(2)=MR_star(2,1) ! Rstar


    do j=2,NB
      j1=(j-2)*8
      upar=get_unit(cpuid)
      open(upar,file=trim(path)//trim(bfiles(j)),status='OLD')
      
      ! Mass min & max of the planet [0-13]Mjup
      read(upar,*) par_min(3+j1),par_max(3+j1) !Mass
      if((par_max(3+j1).eq.0._dp).or.(par_max(3+j1).le.par_min(3+j1)))then
        ! Mmax =0 or Mmax <= Mmin --> Mmax = 2 * Mmin
        par_max(3+j1)=2._dp*par_min(3+j1)
      else if(par_min(3+j1).lt.0._dp)then
        ! Mmin < 0 --> Mmin = 0 , Mmax = 13 Mjup
        par_min(3+j1)=zero
        par_max(3+j1)=13._dp
!       else if(par_max(3+j1).gt.13._dp)then
        ! Mmax > 13 Mjups --> only alert!
!         write(*,'(a)')" WARNING YOU HAVE ENTERED A MASS GREATER THAN 13 MJup! "
      end if
      ! Mjup --> Msun
      par_min(3+j1)=par_min(3+j1)*Mjups
      par_max(3+j1)=par_max(3+j1)*Mjups

      ! Radius mjn = max = Radius planet
!       read(upar,*)temp !Radius
!       par_min(4+j1)=R(j)
!       par_max(4+j1)=R(j)
      read(upar,*)temp1,temp2 !Radius
      par_min(4+j1)=max(min(temp1,temp2),zero)*Rjups
      par_max(4+j1)=min(max(temp1,temp2),5._dp)*Rjups
      if(par_max(4+j1).lt.par_min(4+j1)) par_max(4+j1)=par_min(4+j1)

      ! Read Period min & max of the planet from file
      read(upar,*) par_min(5+j1),par_max(5+j1) !Period
      read(upar,*) temp1,temp2 !semi-major axis
      if(par_min(5+j1).ge.9e6_dp)then
        par_min(5+j1)=period(m(1),m(j),temp1)
        par_max(5+j1)=period(m(1),m(j),temp2)
      end if
      if((par_max(5+j1).eq.0._dp).or.(par_max(5+j1).le.par_min(5+j1)))then
        par_max(5+j1)=2._dp*par_min(5+j1)
!         write(*,*)" **WARNING** body ",j," has Period max not setted correctly"
!         write(*,*)" probably it has been setted to 0. or <= Period min"
!         write(*,*)" it will be set to 2*Pmin "
      end if
      Pvec(j-1)=par_max(5+j1)

      ! 2015-03-11
      ! TRADES WILL NOT FIT ANYMORE e AND w DISTINGUISHED,
      ! NOW IT WIL FIT A COMBINATION OF e AND w:
      ! AND THE MIN AND MAX WILL BE SET TO = [-1, 1] FOR BOTH
      !par_min(6+j1)=-one
      !par_max(6+j1)=one
      !par_min(7+j1)=-one
      !par_max(7+j1)=one
      !! skipping rows
      !read(upar,*)
      !read(upar,*)

      ! 2015-03-30
      ! eccentricity and arg. pericenter lines: set upper limit
      read(upar,*) temp1,temp2
      
! !       write(*,*)" e1 = ", e1," e2 = ", e2
! !       write(*,*)" e_bounds(:,j)[0] = ",e_bounds(:,j)
!       e_bounds(1,j) = max(min(e1,e2),TOLERANCE)
!       e_bounds(2,j) = min(max(e1,e2),1._dp-TOLERANCE)
!       !       write(*,*)" e_bounds(:,j)[1] = ",e_bounds(:,j)
!       par_min(6+j1)=-e_bounds(2,j)
!       par_max(6+j1)=e_bounds(2,j)
!       par_min(7+j1)=par_min(6+j1)
!       par_max(7+j1)=par_max(6+j1)
!       if((tofit(6+j1).eq.0).and.(tofit(7+j1).eq.1))then
!         par_min(7+j1)=zero
!         par_max(7+j1)=360._dp
!       end if
!       read(upar,*) ! skip line
      
      ! eccentricity := [0,1] max
      e_bounds(1,j) = max(min(temp1,temp2),zero)
      e_bounds(2,j) = min(max(temp1,temp2),one)
      par_min(6+j1)=e_bounds(1,j)
      par_max(6+j1)=e_bounds(2,j)
      ! argument of pericentre := [0,360] max
      read(upar,*)temp1,temp2
      par_min(7+j1)=min(temp1,temp2)
      par_max(7+j1)=max(temp1,temp2)
      
      
      ! OLD VERSION
!       eccentricity
!       read(upar,*) par_min(6+j1),par_max(6+j1)
!       if(par_min(6+j1).ge.par_max(6+j1))then
!         par_min(6+j1)=zero
!         par_max(6+j1)=0.95_dp
!       end if
!       if(par_min(6+j1).lt.0._dp)then
!         par_min(6+j1)=zero
!         par_max(6+j1)=0.95_dp
!       end if
!       ! argument of the pericenter
!       read(upar,*) par_min(7+j1),par_max(7+j1)
!       if(par_min(7+j1).ge.par_max(7+j1)) par_max(7+j1)=par_min(7+j1)+one

      ! mean anomaly
!       read(upar,*) par_min(8+j1),par_max(8+j1)
!       if(par_min(8+j1).ge.par_max(8+j1)) par_max(8+j1)=par_min(8+j1)+one
      read(upar,*)temp1,temp2
      par_min(8+j1)=min(temp1,temp2)
      par_max(8+j1)=max(temp1,temp2)
      
      ! time pericenter: skip
      read(upar,*)
      
      ! inclination
!       read(upar,*) par_min(9+j1),par_max(9+j1)
!       if(par_min(9+j1).ge.par_max(9+j1)) par_max(9+j1)=par_min(9+j1)+one
      read(upar,*)temp1,temp2
      par_min(9+j1)=max(min(temp1,temp2),zero)
      par_max(9+j1)=min(max(temp1,temp2),180._dp)
      
      ! longitude of node
!       read(upar,*) par_min(10+j1),par_max(10+j1)
!       if(par_min(10+j1).ge.par_max(10+j1)) par_max(10+j1)=par_min(10+j1)+one
      read(upar,*)temp1,temp2
      par_min(10+j1)=min(temp1,temp2)
      par_max(10+j1)=max(temp1,temp2)

      
      close(upar)
    end do

    amin=R(1)*RsunAU
    amax=5._dp*semax(m(1),zero,maxval(Pvec))
    deallocate(Pvec)
    
    call set_minmax() ! it modifies minpar and maxpar
    
    
    return
  end subroutine read_par_boundaries

  ! USING THIS SUBROUTINE
  ! subroutine that reads orbital elements from the files in bodies.lst
  ! but always assuming:
  ! col 1 == orbital parameters / initial guess
  ! col 2 == minimum value of the orbital elements
  ! col 3 == maximum value of the orbital elements
  ! in case of grid the files will be re-read in a different way
  subroutine read_fullpar(cpuid,m,R,P,a,e,w,mA,inc,lN,all_parameters)
    use celestial_mechanics,only:semax,period,tau2mA!,mA2tau
    integer,intent(in)::cpuid
    real(dp),dimension(:),allocatable,intent(out)::m,R,P,a,e,w,mA,inc,lN
    real(dp),dimension(:),allocatable,intent(out)::all_parameters
    
    real(dp),dimension(:),allocatable::tau
    real(dp),dimension(:),allocatable::Pvec
    integer::unit,j,j1
    logical::fstat,bstat
    character(512)::temp_line
    real(dp)::temp1,temp2
    
    MR_star = zero ! init MR_star to zero value: in 'parameters' module
    
    inquire(file=trim(path)//trim(bfiles(1)),exist=fstat)
    if(fstat)then
      call init_zero_par(m,R,P,a,e,w,mA,inc,lN,tau)
      !read and convert star mass from Msun to Mjup
      unit=get_unit(cpuid)
      open(unit,file=trim(path)//bfiles(1),status='OLD')
!       read(unit,*) m(1) ! Msun
!       read(unit,*) R(1) ! Rsun
      
      read(unit,*) MR_star(1,1),temp_line ! Msun eMsun
      temp_line=trim(adjustl(temp_line))
      if(temp_line(1:1).ne.'#') read(temp_line,'(es23.16)') MR_star(1,2)
      
      read(unit,*) MR_star(2,1),temp_line ! Rsun eRsun
      temp_line=trim(adjustl(temp_line))
      temp_line=trim(adjustl(temp_line))
      if(temp_line(1:1).ne.'#') read(temp_line,'(es23.16)') MR_star(2,2)
      
      close(unit)

      allocate(Pvec(NB-1))
      if(.not.allocated(all_parameters)) allocate(all_parameters(npar),par_min(npar),par_max(npar))
      all_parameters=zero
      par_min=zero
      par_max=zero
      
      ! Mstar
      m(1)=MR_star(1,1)
      all_parameters(1)=MR_star(1,1)
      par_min(1)=max(MR_star(1,1)-10._dp*MR_star(1,2),zero)
      par_max(1)=MR_star(1,1)+10._dp*MR_star(1,2)
      
      ! Rstar
      R(1)=MR_star(2,1)
      all_parameters(2)=MR_star(2,1)
      par_min(2)=max(MR_star(2,1)-10._dp*MR_star(2,2),zero)
      par_max(2)=MR_star(2,1)+10._dp*MR_star(2,2)
      
      
      readpar: do j=2,NB
        j1=(j-2)*8
        
        unit=get_unit(cpuid)
        inquire(file=trim(path)//trim(bfiles(j)),exist=bstat)
        
        if(bstat)then
          open(unit,file=trim(path)//trim(bfiles(j)),status='OLD')
          
          ! Mass
          read(unit,*) m(j),temp1,temp2
          m(j)=m(j)*Mjups
          all_parameters(3+j1)=m(j) ! Mjup to Msun
          par_min(3+j1)=max(min(temp1,temp2)*Mjups,TOLERANCE) ! Mjup to Msun
          par_max(3+j1)=min(max(temp1,temp2)*Mjups,one) ! Mjup to Msun
          if(par_max(3+j1).le.TOLERANCE) par_max(3+j1)=one ! 1 Msun
          
          ! Radius
          read(unit,*) R(j),temp1,temp2
          R(j)=R(j)*Rjups
          all_parameters(4+j1)=R(j)*Rjups ! Rjup to Rsun
          par_min(4+j1)=max(min(temp1,temp2)*Rjups,TOLERANCE)
          par_max(4+j1)=min(max(temp1,temp2),5._dp)*Rjups
          if(par_max(4+j1).le.TOLERANCE) par_max(4+j1)=5._dp*Rjups
          
          ! Period & semi-major axis
          read(unit,*) P(j),temp1,temp2
          if(P(j).ge.9.e6_dp)then ! set high value for Period --> using semi-major axis
            read(unit,*) a(j),temp1,temp2
            P(j)=period(m(1),m(j),a(j))
            temp1=period(m(1),m(j),temp1)
            temp2=period(m(1),m(j),temp2)
          else
            read(unit,*) ! skip semi-major axis row
            a(j)=semax(m(1),m(j),P(j))
          end if
          all_parameters(5+j1)=P(j)
          par_min(5+j1)=max(min(temp1,temp2),TOLERANCE)
          par_max(5+j1)=min(max(temp1,temp2),1000._dp*365.25_dp)
          if(par_max(5+j1).le.TOLERANCE) par_max(5+j1)=1000._dp*365.25_dp
          
          Pvec(j-1)=par_max(5+j1)
          
          ! eccentricity
          read(unit,*) e(j),temp1,temp2
          all_parameters(6+j1)=e(j)
          e_bounds(1,j) = max(min(temp1,temp2),zero)
          e_bounds(2,j) = min(max(temp1,temp2),one-TOLERANCE)
          if(e_bounds(2,j).le.TOLERANCE) e_bounds(2,j)=one-TOLERANCE
          par_min(6+j1)=e_bounds(1,j)
          par_max(6+j1)=e_bounds(2,j)
          
          
          ! argument of pericentre
          read(unit,*) w(j),temp1,temp2
          all_parameters(7+j1)=w(j)
          if(abs(temp1-temp2).le.TOLERANCE)then ! if min and max are equals, set to default range := [0,360] deg
            par_min(7+j1)=zero
            par_max(7+j1)=360._dp
          else
            par_min(7+j1)=min(temp1,temp2)
            par_max(7+j1)=max(temp1,temp2)
          end if
          
          ! mean Anomaly & time of the passage at pericentre
          read(unit,*) mA(j),temp1,temp2
          if(mA(j).ge.999._dp)then
            read(unit,*)tau(j),temp1,temp2
            
            if(tau(j).ge.9.e8_dp)then
              deallocate(tau)
              write(*,'(a,a)')' WARNING: MISSING BOTH MEAN ANOMAMLY AND TIME OF PERICENTER FOR PLANET INPUT FILE ',trim(bfiles(j))
              stop
            end if
            
            mA(j)=tau2mA(tau(j),tepoch,P(j))
            temp1=zero
            temp2=360._dp
            write(*,'(a)')' WARNING: USER PROVIDED TAU SO THE MIN-MAX OF THE MEAN ANOMALY WILL BE SET TO [0-360] DEG'
            flush(6)
            
          else
            read(unit,*) ! skip tau row
  !           tau(j)=mA2tau(mA(j),tepoch,P(j))
          end if
          all_parameters(8+j1)=mA(j)
          if(abs(temp1-temp2).le.TOLERANCE)then ! if min and max are equals, set to default range := [0,360] deg
            par_min(8+j1)=zero
            par_max(8+j1)=360._dp
          else
            par_min(8+j1)=min(temp1,temp2)
            par_max(8+j1)=max(temp1,temp2)
          end if
          
          ! inclination
          read(unit,*) inc(j),temp1,temp2
          all_parameters(9+j1)=inc(j)
          if(abs(temp1-temp2).le.TOLERANCE)then ! if min and max are equals, set to default range := [0,180] deg
            par_min(9+j1)=zero
            par_max(9+j1)=180._dp
          else
            par_min(9+j1)=max(min(temp1,temp2),zero)
            par_max(9+j1)=min(max(temp1,temp2),180._dp)
            if(par_max(9+j1).le.TOLERANCE) par_max(9+j1)=180._dp
          end if
          
          ! longitude of ascending node
          read(unit,*)lN(j),temp1,temp2
          all_parameters(10+j1)=lN(j)
          if(abs(temp1-temp2).le.TOLERANCE)then ! if min and max are equals, set to default range := [0,360] deg
            par_min(10+j1)=zero
            par_max(10+j1)=360._dp
          else
            par_min(10+j1)=min(temp1,temp2)
            par_max(10+j1)=max(temp1,temp2)
          end if

          close(unit)
          
          
        else
        
          deallocate(tau)
          write(*,'(a,a,a)')" CANNOT FIND BODY FILE ",trim(path),trim(bfiles(j))
          write(*,'(a)')""
          stop
      
        end if
        
      end do readpar
      
    else
    
      deallocate(tau)
      write(*,'(a,a,a)')" CANNOT FIND STAR FILE ",trim(path),trim(bfiles(1))
      write(*,'(a)')""
      stop
    end if
    
    deallocate(tau)
    
    amin=R(1)*RsunAU
    amax=5._dp*semax(m(1),zero,maxval(Pvec))
    deallocate(Pvec)
    
    call set_minmax() ! it modifies minpar and maxpar

    return
  end subroutine read_fullpar
  
  
!   ! fix the system_parameters in case a parameter has been read with a value not in [par_min, par_max]
! 
!   subroutine fix_system_parameters()
!     integer::i
!     
!     do i=1,npar
!       if( tofit(i).eq. 1)then
!         if( (system_parameters(i).lt.par_min(i)) .or. (system_parameters(i).gt.par_max(i)) )then
!           system_parameters(i) = par_min(i)
!         end if
!       end if
!     end do
!   
!     return
!   end subroutine fix_system_parameters

  
  ! read and initialize RV data
  ! FILE TYPE: JD RVobs err_RVobs
  subroutine read_RVobs(cpuid)
    integer,intent(in)::cpuid
    integer::j,urv,stat,nTemp,jSet
    logical::fstat
    character(512)::row

    nRV=0
    nRVset=1
    if(rvcheck.eq.1)then
      inquire(file=trim(path)//"obsRV.dat",exist=fstat)
      if(fstat)then
        urv=get_unit(cpuid)
        open(urv,file=trim(path)//"obsRV.dat",status='OLD')
        nRV=get_rows(urv)
        if(.not.allocated(jdRV)) &
        &allocate(jdRV(nRV),RVobs(nRV),eRVobs(nRV),RVsetID(nRV))
!         do j=1,nRV
        j=0
        RVdo:do
          read(urv,'(a512)',IOSTAT=stat)row
          if(IS_IOSTAT_END(stat)) exit RVdo
          row=trim(adjustl(row))
          if(row(1:1).ne."#")then
            j=j+1
            read(row,*)jdRV(j),RVobs(j),eRVobs(j),RVsetID(j)
            if(j.gt.1)then
              if(RVsetID(j).ne.RVsetID(j-1)) nRVset=nRVset+1
            end if
          end if
        end do RVdo
        close(urv)
        allocate(nRVsingle(nRVset))
        nTemp=1
        jSet=1
        do j=2,nRV
          if(RVsetID(j).ne.RVsetID(j-1))then
            nRVsingle(jSet)=nTemp
            jSet=jSet+1
            nTemp=0
          end if
          nTemp=nTemp+1
        end do
        nRVsingle(jSet)=nTemp
      else
        write(*,'(a,a,a)')" CANNOT FIND RV FILE ",trim(path),"obsRV.dat"
        write(*,'(a,a,a)')" CHECK YOUR RV FILE OR SET TO 0 THE RELATED OPTION IN ",&
        &trim(path),"arg.in"
        write(*,'(a)')""
        stop
      end if
    else
      write(*,'(a)')" SELECTED NO RV CHECK "
    end if

    return
  end subroutine read_RVobs

  ! read and initialize T0 data 
  !! TO DO: READ DURATION FROM FILE => nT0 x 2, allocation of variables, read file change
  ! FILE TYPE: N T_0,obs err_T_0,obs
  subroutine read_T0obs(cpuid)
    integer,intent(in)::cpuid
    integer::nTmax
    integer::j,j1,uread,stat
    logical::fstat
    character(512)::flt0,row

    flt0=""
    if(.not.allocated(nT0)) allocate(nT0(NB))
    nT0=0
    
    if(idtra.ne.0)then
      do j=2,NB
        flt0=trim(path)//'NB'//trim(adjustl(string(j)))//'_observations.dat'
        flt0=trim(adjustl(flt0))
        inquire(file=trim(flt0),exist=fstat)
        if(fstat)then
          uread=get_unit(cpuid)
          open(uread,file=trim(flt0),status='OLD')
          nT0(j)=get_rows(uread)
          close(uread)
        else
          write(*,'(a,a)')" CANNOT FIND T_0 FILE ",trim(flt0)
        end if
      end do
    
      nTmax=maxval(nT0)
      if(.not.allocated(T0obs))then
        allocate(T0obs(nTmax,NB),eT0obs(nTmax,NB),epoT0obs(nTmax,NB))
        T0obs=zero
        eT0obs=zero
        do j=2,NB
          flt0=""
          flt0=trim(path)//'NB'//trim(adjustl(string(j)))//'_observations.dat'
          flt0=trim(adjustl(flt0))
          inquire(file=trim(flt0),exist=fstat)
          if(fstat)then
            uread=get_unit(cpuid)
            open(uread,file=trim(flt0),status='OLD')
!             do j1=1,nT0(j)
            j1=0
            T0do:do
              read(uread,'(a512)',IOSTAT=stat)row
              if(IS_IOSTAT_END(stat)) exit T0do
              row=trim(adjustl(row))
              if(row(1:1).ne."#")then
                j1=j1+1
                read(row,*)epoT0obs(j1,j),T0obs(j1,j),eT0obs(j1,j)
              end if
            end do T0do
            close(uread)
          end if
        end do
      end if

    else
      write(*,'(a)')" SELECTED NO T_0 DATA "
    end if

    return
  end subroutine read_T0obs

  function get_ln_err_const(eRV,eT0) result(ln_const)
    real(dp)::ln_const
    real(dp),dimension(:),intent(in)::eRV
    real(dp),dimension(:,:),intent(in)::eT0
    
    real(dp)::ln_eRV,ln_eT0
    
    if(nRV.gt.0)then
      ln_eRV = sum(log(pack(eRV,eRV/=zero)*pack(eRV,eRV/=zero)))
    else
      ln_eRV = zero
    end if
    if(sum(nT0).gt.0)then
      ln_eT0 = sum(log(pack(eT0,eT0/=zero)*pack(eT0,eT0/=zero)))
    else
      ln_eT0 = zero
    end if
    ln_const = -(0.5_dp*real(dof,dp)*log(dpi))-(0.5_dp*(ln_eRV+ln_eT0))
!     write(*,'(a,es23.16)')' LN_ERR_CONST (SUBROUTINE) = ',ln_err_const
!     flush(6)
  
  end function get_ln_err_const
  
  
  ! it reads the parameters for PIKAIA simulation
  subroutine read_pik_opt(cpuid)
    integer,intent(in)::cpuid
    character(512)::flopt
    integer::uread,j
    logical::fstat


    uread=get_unit(cpuid)
    flopt=""
    flopt=trim(path)//"pikaia.opt"
    flopt=trim(adjustl(flopt))
    inquire(file=trim(flopt),exist=fstat)
    if(fstat)then
      open(uread,file=trim(flopt),status='OLD')
      do j=1,12
        read(uread,*)ctrl(j)
      end do
      read(uread,*)seed_pik
      if(seed_pik.le.0)seed_pik=123456
      read(uread,*)wrtAll
      read(uread,*)nGlobal
      close(uread)
      npop=int(ctrl(1))
      ngen=int(ctrl(2))
!       write(*,'(a,i3)')" Read nGlobal = ",nGlobal
!       stop
      if(nGlobal.lt.1) nGlobal=1
    else
      if(progtype.eq.3)then
        write(*,'(a,a,a)')" CANNOT FIND FILE ",trim(flopt)," FOR PIKAIA OPTIONS"
        stop
      end if
    end if

    return
  end subroutine read_pik_opt

  ! it reads the parameters for PSO
  subroutine read_pso_opt(cpuid)
    integer,intent(in)::cpuid
    character(512)::flopt
    integer::uread
    logical::fstat


    uread=get_unit(cpuid)
    flopt=""
    flopt=trim(path)//"pso.opt"
    flopt=trim(adjustl(flopt))
    inquire(file=trim(flopt),exist=fstat)
    if(fstat)then
      open(uread,file=trim(flopt),status='OLD')
      read(uread,*)np_pso
      read(uread,*)nit_pso
      read(uread,*)wrt_pso
      read(uread,*)wrtAll
      read(uread,*)nGlobal
      read(uread,*)seed_pso
      if(seed_pso.le.0) seed_pso=123456
      close(uread)
      if(nGlobal.lt.1) nGlobal=1
    else
      if(progtype.eq.3)then
        write(*,'(a,a,a)')" CANNOT FIND FILE ",trim(flopt)," FOR PSO OPTIONS"
        stop
      end if
    end if

    return
  end subroutine read_pso_opt

!   ! it sets the boundaries for the PSO simulation [not only]
!   subroutine set_minmax()
!     integer::ifit,ipar
! 
!     if(.not.allocated(minpar)) allocate(minpar(nfit),maxpar(nfit))
!     ifit=0
!     do ipar=1,npar
!       if(tofit(ipar).eq.1)then
!         ifit=ifit+1
!         minpar(ifit)=par_min(ipar)
!         maxpar(ifit)=par_max(ipar)
!       end if
!     end do
! 
!     return
!   end subroutine set_minmax
!   ! ------------------------------------------------------------------ !

  
  subroutine get_character_fields(line, variable)
    character(*),intent(in)::line
    real(8),intent(out),dimension(:),allocatable::variable
    character(512)::line_temp
    integer::n_len,n_space,n_elements
    integer::i
    
    line_temp = trim(adjustl(line))
    n_len = len_trim(line_temp)
    n_space = 0
    do i=1,n_len
      if(line_temp(i:i) .eq. ' ')then
        if(line_temp(i+1:i+1) .ne. ' ') n_space = n_space + 1
      end if
    end do
    n_elements = n_space + 1
    allocate(variable(n_elements))
    read(line_temp, *) variable

    return
  end subroutine get_character_fields
 
! read and set settings for PolyChord
  subroutine read_PC_opt(cpuid)
    integer,intent(in)::cpuid
    logical::fstat
    integer::upc,istat
    character(512)::line
    integer::idh,idx
    
    ! type(program_settings)::settings HAS BEEN DEFINED IN parameters.f90
    settings%nDims = nfit ! number of dimensions of PolyChord is the number of parameters to fit
    settings%nDerived = 0 ! none derived parameters to pass now
    
    inquire(file=trim(path)//'PolyChord.opt',exist=fstat)
    
    if(fstat)then
      upc=get_unit(cpuid)
      open(upc,file=trim(path)//'PolyChord.opt',status='OLD')
        
        reading:do
          read(upc,'(a512)',IOSTAT=istat)line
          if(IS_IOSTAT_END(istat)) exit reading
      
          if(len(trim(line)).ne.0)then ! look for non empty line
            
            if(line(1:1).ne."#")then ! if not commented line with a #
              idh=index(line,"#") ! find a # in the line, not at the beginning, further comment
              if(idh.eq.0) idh=len(trim(line))+1 ! if there are none # use the whole trimmed line
              
              idx=index(line,"=") ! find the = to distinguish keyword (left) and value (right)
              if(idx.ne.0)then
              
                if(line(1:idx-1).eq.'nlive')then
                  read(line(idx+1:idh),*) settings%nlive
                else if(line(1:idx-1).eq.'num_repeats')then
                  read(line(idx+1:idh),*) settings%num_repeats
                else if(line(1:idx-1).eq.'do_clustering')then
                  read(line(idx+1:idh),*) settings%do_clustering

                else if(line(1:idx-1).eq.'precision_criterion')then
                  read(line(idx+1:idh),*) settings%precision_criterion
                else if(line(1:idx-1).eq.'max_ndead')then
                  read(line(idx+1:idh),*) settings%max_ndead
                else if(line(1:idx-1).eq.'feedback')then
                  read(line(idx+1:idh),*) settings%feedback

                else if(line(1:idx-1).eq.'posteriors')then
                  read(line(idx+1:idh),*) settings%posteriors
                else if(line(1:idx-1).eq.'equals')then
                  read(line(idx+1:idh),*) settings%equals
                else if(line(1:idx-1).eq.'cluster_posteriors')then
                  read(line(idx+1:idh),*) settings%cluster_posteriors
                else if(line(1:idx-1).eq.'update_posterior')then
                  read(line(idx+1:idh),*) settings%update_posterior
                else if(line(1:idx-1).eq.'boost_posterior')then
                  read(line(idx+1:idh),*) settings%boost_posterior

                else if(line(1:idx-1).eq.'read_resume')then
                  read(line(idx+1:idh),*) settings%read_resume
                else if(line(1:idx-1).eq.'write_resume')then
                  read(line(idx+1:idh),*) settings%write_resume
                else if(line(1:idx-1).eq.'write_live')then
                  read(line(idx+1:idh),*) settings%write_live
                else if(line(1:idx-1).eq.'write_stats')then
                  read(line(idx+1:idh),*) settings%write_stats

                else if(line(1:idx-1).eq.'base_dir')then
                  read(line(idx+1:idh),*) settings%base_dir
                else if(line(1:idx-1).eq.'file_root')then
                  read(line(idx+1:idh),*) settings%file_root
                
                else if(line(1:idx-1).eq.'grade_frac')then
                  call get_character_fields(line(idx+1:idh), settings%grade_frac)
                
                end if
                      
              end if
              
            end if
          
          end if
          
        end do reading
      close(upc)
      
    else
    
      write(*,'(a,a,a)')" CANNOT FIND SETTINGS FILE ",trim(path),"PolyChord.opt"
      write(*,'(a,a,a)')" DEFAULT SETTINGS WILL BE USED"
    end if
    
    settings%update_resume = settings%nlive
    ! add absolute path 'path' to base_dir
    settings%base_dir=trim(path)//settings%base_dir
    ! create recursively the base_dir and the base_dir/'clusters' folder
    call execute_command_line('mkdir -p '//trim(settings%base_dir)//"/clusters")
  
    call set_minmax()
  
    return
  end subroutine read_PC_opt

  ! calls all the read subroutines necessary to initialize all the variables/vectors/arrays
  subroutine read_first(cpuid,m,R,P,a,e,w,mA,inc,lN)
    integer,intent(in)::cpuid
    real(dp),dimension(:),allocatable,intent(out)::m,R,P,a,e,w,mA,inc,lN
    integer::j
    character(80)::fmt
    integer,dimension(:),allocatable::nset

    ! IT DEFINES THE STRING TO WRITE THE REAL WITH RIGHT DECIMAL: PRECISION
!     sprec=trim(adjustl("g"//&
!     &trim(adjustl(string(2*prec)))//&
!     &"."//&
!     &trim(adjustl(string(prec)))))
    sprec='es23.16'
    ! prec, sprec in module 'constants'
    write(*,'(a,a,a,a)')" SELECTED PATH: ",trim(path),&
    &" FORMAT EDIT DESCRIPTOR: ",trim(sprec)
  
    ! IT READS THE ARGUMENTS OF INTEGRATION AND STORE IN COMMON MODULE PARAMETERS.
    ! THE VARIBLES WILL NOT BE MODIFIED FURTHERMORE.
    call read_arg(cpuid)
    write(*,'(a,a,a)')" READ ",trim(path),"arg.in"
    allocate(e_bounds(2,NB))
    e_bounds(1,:)=zero
    e_bounds(2,:)=one-TOLERANCE

    ! IT READS THE FILES AND THE NAMES OF THE BODIES AND DETERMINES THE PARAMETERS TO BE FITTED
    call read_list(cpuid)
    write(*,'(a,a,a)')" READ ",trim(path),"bodies.lst"
    write(*,'(a,a)')" NUMBER OF PARAMETERS TO FIT: nfit = ",trim(adjustl(string(nfit)))
    do j=1,NB
      write(*,'(a,1000a)')" BODY NAME: ",trim(bnames(j)),&
      &" FILE NAME: ",trim(path),trim(bfiles(j))
    end do
    
    call idpar() ! IT DEFINES THE ID OF THE PARAMETERS TO BE FITTED
     
    ! IT READS THE VARIABLES FOR THE LEVENBERG-MARQUARDT ALGORITHM
    call read_lm_opt(cpuid)
    write(*,'(a,a,a)')" READ ",trim(path),"lm.opt"

    ! IT READS THE PARAMETERS FROM THE FILES
!     call read_par(cpuid,m,R,P,a,e,w,mA,inc,lN)
    call read_fullpar(cpuid,m,R,P,a,e,w,mA,inc,lN,system_parameters)
    
    if(progtype.gt.1)then
      write(*,'(a)')'Initial-input Keplerian orbital elements'
      write(*,'(a, 1000(1x,es23.16))')"m   [Msun] = ", m(1),m(2:)
      write(*,'(a, 1000(1x,es23.16))')"R   [Rsun] = ", R
      write(*,'(a, 1000(1x,es23.16))')"P   [days] = ", P
      write(*,'(a, 1000(1x,es23.16))')"a   [au]   = ", a
      write(*,'(a, 1000(1x,es23.16))')"e          = ", e
      write(*,'(a, 1000(1x,es23.16))')"w   [deg]  = ", w
      write(*,'(a, 1000(1x,es23.16))')"mA  [deg]  = ", mA
      write(*,'(a, 1000(1x,es23.16))')"inc [deg]  = ", inc
      write(*,'(a, 1000(1x,es23.16))')"lN  [deg]  = ", lN
      write(*,'(a)')''
      
      write(*,'(a23,2(1x,a23))')'Parameters','( par_min ,','par_max )'
      do j=1,npar
        write(*,'(es23.16,2(a,es23.16),a)')system_parameters(j),' ( ',par_min(j),' , ',par_max(j),' )'
      end do
    end if
    
    nGlobal=1
!     if(progtype.ge.3)then
    
!       call read_par_boundaries(cpuid,m,R)
    if(progtype.eq.3)then
      ! IT READS THE VARIABLES FOR THE PIKAIA (GENETIC ALGORITHM) CODE
      call read_pik_opt(cpuid)
      write(*,'(a,a,a)')" READ ",trim(path),"pikaia.opt"
    else if(progtype.eq.4)then
      ! IT READS THE VARIABLES FOR THE PARTICLE SWARM CODE
      call read_pso_opt(cpuid)
      write(*,'(a,a,a)')" READ ",trim(path),"pso.opt"
!         call set_minmax() ! it modifies minpar and maxpar to be used with pso
    else if(progtype.eq.5)then
      call read_PC_opt(cpuid)
      write(*,'(a,a,a)')" READ ",trim(path),"PolyChord.opt"
    end if
      
!     end if

    ! IT READS RV DATA
    nRV=0
    call read_RVobs(cpuid)
    fmt=adjustl("(a,i4,a,i4,a,3i4))")
    if(rvcheck.eq.1) write(*,trim(fmt))" RV DATA: nRV = ",nRV,&
      &" in ",nRVset," set of RV: ",nRVsingle

    ! IT READS T0 DATA
    call read_T0obs(cpuid)
    if(idtra.ne.0) write(*,'(a,1000(i5,1x))') " T0 DATA: nT0 = ",nT0(2:)
    if(idtra.ne.0) write(*,'(a,1000(l2,1x))') ' DO TRANSIT FLAG: ',do_transit
    
    ! IT DETERMINES THE NDATA
    ndata=nRV+sum(nT0)
    dof=(ndata-nfit)
    inv_dof = one / real(dof,dp)
    write(*,'(a,i5)')" NUMBER OF DATA AVAILABLE: ndata = ",ndata
    write(*,'(a,i5)')" NUMBER OF PARAMETERS TO FIT: nfit = ",nfit
    write(*,'(a,i5)')&
        &" NUMBER OF DEGREES OF FREEDOM : dof = ndata - nfit = ",dof

    ! IT DETERMINES THE LN_ERR_CONST TO COMPUTE LOGLIKELIHOOD
    ln_err_const = get_ln_err_const(eRVobs,eT0obs)
!     write(*,'(a,es23.16)')' LN_ERR_CONST (init_trades) = ',ln_err_const
!     flush(6)
        
    call set_parid_list()

    ! set fitness parameters
    if(nRV.ne.0.and.sum(nT0).ne.0)then
      allocate(nset(2),k_b(2))
      nset(1)=nRV
      nset(2)=sum(nT0)
    else if(nRV.ne.0.and.sum(nT0).eq.0)then
      allocate(nset(1),k_b(1))
      nset(1)=nRV
    else if(nRV.eq.0.and.sum(nT0).ne.0)then
      allocate(nset(1),k_b(1))
      nset(1)=sum(nT0)
    else
      stop('No data-set available. Please check the files.')
    end if
    ! parameter to properly scale the residuals for the Chi2
    k_a = sqrt(k_chi2r*inv_dof)
    ! parameter to properly scale the residuals for the Chi2_weighted
!     k_b = sqrt(k_chi2wr/real(dof,dp))*(real(ndata,dp)/real(nset,dp))
    k_b = sqrt((k_chi2wr*inv_dof)*(real(ndata,dp)/real(nset,dp)))
    write(*,'(2(a,f7.4))')" k_chi2r = ",k_chi2r," k_chi2wr = ",k_chi2wr
    if(size(nset).eq.2)then
      write(*,'(a,i5,i5,a,f16.12,a,f16.12,f16.12)')" nset = ",nset," k_a = ",k_a," k_b = ",k_b
    else if(size(nset).eq.1)then
      write(*,'(a,i5,a,f16.12,a,f16.12)')" nset = ",nset," k_a = ",k_a," k_b = ",k_b
    end if
    deallocate(nset)
    
    return
  end subroutine read_first
  ! ------------------------------------------------------------------ !
  
  ! it creates string with all the state vector names, useful as header of orbit file
  function state2string(Nbody) result(out)
    character(512)::out
    integer,intent(in)::Nbody
    integer::i

    out=""
    do i=1,Nbody
      out=trim(adjustl(out))//&
          &" X_"//trim(adjustl(string(i)))//"_AU"//&
          &" Y_"//trim(adjustl(string(i)))//"_AU"//&
          &" Z_"//trim(adjustl(string(i)))//"_AU"//&
          &" VX_"//trim(adjustl(string(i)))//"_AU"//&
          &" VY_"//trim(adjustl(string(i)))//"_AU"//&
          &" VZ_"//trim(adjustl(string(i)))//"_AU"
    end do

    return
  end function state2string
  
end module init_trades

