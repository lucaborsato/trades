module parameters_conversion
  use constants
  use parameters
  use convert_type
  


  contains

  subroutine set_all_parameter_names(kel_id)
    character(6),dimension(3:10),intent(in)::kel_id
  
    integer::i_par,body,cnt_kel
    allocate(all_names_list(npar)) ! all_names_list in 'parameters' module
    
    all_names_list(1)=trim(adjustl(kel_id(3)))//'1' ! mass of the star, id 1
    all_names_list(2)=trim(adjustl(kel_id(4)))//'1' ! radius of the star, id 1
    all_names_str=trim(adjustl(all_names_list(1)))//' '//trim(adjustl(all_names_list(2)))
    cnt_kel=3
    do i_par=3,npar
      body=int(int(i_par-3)/8)+2
      all_names_list(i_par)=trim(adjustl(kel_id(cnt_kel)))//trim(adjustl(string(body)))
      all_names_str=trim(adjustl(all_names_str))//' '//trim(adjustl(all_names_list(i_par)))
      cnt_kel=cnt_kel+1
      if(cnt_kel.eq.11) cnt_kel=3
    end do
!     write(*,'(a)')trim(adjustl(all_names_str))
  
    return
  end subroutine set_all_parameter_names
  
!   ! determines the id of the parameters to be fitted ... only for helpful write
!   subroutine idpar()
!     integer::pos,cntid,j,body
!     character(5),dimension(3:10)::elid
!     data elid /"m", "R", "P", "e", "w", "mA", "i", "lN"/
! 
!     if(.not.allocated(id)) allocate(id(nfit),idall(nfit),parid(nfit))
!     pos=0
!     cntid=2
!     do j=1,npar
!       if(j.gt.2) cntid=cntid+1
!       if(tofit(j).eq.1)then
!         if(cntid.eq.6)then ! ecc
!           tofit(j+1)=1
!           tofit(j+2)=1
!         end if
!         pos=pos+1
!         id(pos)=cntid
!         idall(pos)=j
!         body=int(int(j-3)/8)+2
!         parid(pos)=trim(adjustl(elid(cntid)))//trim(adjustl(string(body)))
!         parid(pos)=trim(adjustl(parid(pos)))
!       end if
!       if(cntid.eq.10) cntid=2
!     end do
! 
!     call set_all_parameter_names(elid)
!     
!     return
!   end subroutine idpar
! 
!   subroutine idpar_fit()
!     integer::pos,cntid,j,body
!     character(5),dimension(3:10)::elid
!     data elid /"m", "R", "P", "ecosw", "esinw", "mA", "i", "lN"/
! 
!     if(.not.allocated(id)) allocate(id(nfit),idall(nfit),parid(nfit))
!     pos=0
!     cntid=2
!     do j=1,npar
!       if(j.gt.2) cntid=cntid+1
!       if(tofit(j).eq.1)then
!         pos=pos+1
!         id(pos)=cntid
!         idall(pos)=j
!         body=int(int(j-3)/8)+2
!         parid(pos)=trim(adjustl(elid(cntid)))//trim(adjustl(string(body)))
!         parid(pos)=trim(adjustl(parid(pos)))
!       end if
!       if(cntid.eq.10) cntid=2
!     end do
! 
!     call set_all_parameter_names(elid)
!     
!     return
!   end subroutine idpar_fit

  subroutine idpar()
    integer::pos,cntid,j,body
!     character(6),dimension(3:10)::elid_1,elid_2,elid_3,elid
    character(6),dimension(3:10)::elid
    data elid/"m", "R", "P", "e", "w", "mA", "i", "lN"/

    if(.not.allocated(id)) allocate(id(nfit),idall(nfit),parid(nfit))
    pos=0
    cntid=2
    do j=1,npar
      if(j.gt.2) cntid=cntid+1
      
      if(tofit(j).eq.1)then
        pos=pos+1
        id(pos)=cntid
        idall(pos)=j
        body=int(int(j-3)/8)+2
        parid(pos)=trim(adjustl(elid(cntid)))//trim(adjustl(string(body)))
        parid(pos)=trim(adjustl(parid(pos)))
      end if
      if(cntid.eq.10) cntid=2
    end do

    call set_all_parameter_names(elid)
    
    do j=1,nfit
      
      body=int(int(idall(j)-3)/8)+2
      
      if(id(j).eq.3)then ! Mp
        parid(j)="m"//trim(adjustl(string(body)))//"Ms"
        
      else if(id(j).eq.6)then
        if(j.lt.nfit)then
          if(id(j+1).eq.7)then !e,w
            parid(j)="ecosw"//trim(adjustl(string(body)))
            parid(j+1)="esinw"//trim(adjustl(string(body)))
          end if
        end if
        
      else if(id(j).eq.8)then ! mA
        parid(j)="lambda"//trim(adjustl(string(body)))
        
      else if(id(j).eq.9)then
        if(j.lt.nfit)then
          if(id(j+1).eq.10)then !i,lN
            parid(j)="icoslN"//trim(adjustl(string(body)))
            parid(j+1)="isinlN"//trim(adjustl(string(body)))
          end if
        end if
         
      end if
      
    end do
    
    return
  end subroutine idpar

  ! given the id of the parameters to fit it creates a proper string
  subroutine set_parid_list()
    integer::i

    paridlist=""
    sig_paridlist=""
    do i=1,nfit
      paridlist=trim(paridlist)//" "//trim(adjustl(parid(i)))
      sig_paridlist=trim(sig_paridlist)//" sig_"//trim(adjustl(parid(i)))
    end do
!     write(*,'(a)')" ID Parameters to fit:"
!     write(*,'(a)')trim(paridlist)
!     write(*,'(a)')trim(sig_paridlist)
!     write(*,*)

    return
  end subroutine set_parid_list

    ! it sets the boundaries for the PSO simulation [not only]
  subroutine set_minmax()
    integer::ifit,ipar
    logical,dimension(:),allocatable::done
    
    allocate(done(nfit))
    done=.false.
    if(.not.allocated(minpar)) allocate(minpar(nfit),maxpar(nfit))
    ifit=0
    do ipar=1,npar
      if(tofit(ipar).eq.1)then
        ifit=ifit+1
        
        
        if(.not.done(ifit))then
  !         ! default setting
  !         minpar(ifit)=par_min(ipar)
  !         maxpar(ifit)=par_max(ipar)
          
          ! fit m(3)==> mp/Ms
          if(id(ifit).eq.3)then
            minpar(ifit)=par_min(ipar)/MR_star(1,1)
            maxpar(ifit)=par_max(ipar)/MR_star(1,1)
            done(ifit)=.true.
            
          ! fit e(6) & w(7)==>(ecosw,esinw) [-e,+e],[-e,+e]
          else if(id(ifit).eq.6)then
  !           minpar(ifit)=par_min(ipar)
  !           maxpar(ifit)=par_max(ipar)
            if(ifit.lt.nfit)then
              if(id(ifit+1).eq.7)then
                minpar(ifit)  = -par_max(ipar)
                minpar(ifit+1)= -par_max(ipar)
                maxpar(ifit)  = par_max(ipar)
                maxpar(ifit+1)= par_max(ipar)
  !               write(*,*)id(ifit), id(ifit+1)
  !               write(*,*)minpar(ifit),maxpar(ifit)
  !               write(*,*)minpar(ifit+1),maxpar(ifit+1)
                done(ifit)=.true.
                done(ifit+1)=.true.
              end if
            end if

          ! fit mA(8) ==> lambda[0,360]
          else if(id(ifit).eq.8)then
!             minpar(ifit)=zero
!             maxpar(ifit)=360._dp
            minpar(ifit)=-360._dp
            maxpar(ifit)=720._dp
            done(ifit)=.true.
          
          ! fit i(9) & lN(10)==>(icoslN,isinlN) [-180,180],[-180,180]
          else if(id(ifit).eq.9)then
            
            if(ifit.lt.nfit)then
              if(id(ifit+1).eq.10)then
                minpar(ifit)  =-180._dp
                minpar(ifit+1)=-180._dp
                maxpar(ifit)  =+180._dp
                maxpar(ifit+1)=+180._dp
                done(ifit)=.true.
                done(ifit+1)=.true.
              else
                minpar(ifit) =par_min(ipar)
                maxpar(ifit) =par_max(ipar)
                done(ifit)=.true.
              end if
            end if
            
!           else if(ifit.gt.1)then
!             if((id(ifit-1).eq.6.and.id(ifit).eq.7).or.&
!               &(id(ifit-1).eq.9.and.id(ifit).eq.10))&
!               &cycle
! 
!           else
          end if
        end if ! .not.done(ifit)

        if(.not.done(ifit))then
          ! default setting
          minpar(ifit)=par_min(ipar)
          maxpar(ifit)=par_max(ipar)
          done(ifit)=.true.
        end if
      
      end if ! tofit(j)
    end do

    return
  end subroutine set_minmax
  ! ------------------------------------------------------------------ !

  
  ! ------------------------------------------------------------------ !
  ! defines limits of the integration as distance from the star
  subroutine sma_boundaries(R,a,smamin,smamax)
    real(dp),dimension(:),intent(in)::R,a
    real(dp),intent(out)::smamin,smamax

    smamin=R(1)*RsunAU !smamin =  Rstar in AU
    smamax=5._dp*maxval(a) !smamax = 5 times the larger semi-major axis

    return
  end subroutine sma_boundaries
  
  
      ! puts the parameters in a big vector with dimension npar = 2+(NB-1)*8 (in tofit)
  subroutine set_all_param(m,R,P,e,w,mA,inc,lN,allpar)
    real(dp),dimension(:),intent(in)::m,R,P,e,w,mA,inc,lN
    real(dp),dimension(:),allocatable,intent(out)::allpar
    integer::j,j1

    if(.not.allocated(allpar)) allocate(allpar(npar))
    allpar(1)=m(1)
    allpar(2)=R(1)
    do j=2,NB
      j1=(j-2)*8
      allpar(3+j1)=m(j)
      allpar(4+j1)=R(j)
      allpar(5+j1)=P(j)
      allpar(6+j1)=e(j)
      allpar(7+j1)=w(j)
      allpar(8+j1)=mA(j)
      allpar(9+j1)=inc(j)
      allpar(10+j1)=lN(j)
    end do

    return
  end subroutine set_all_param
  
   ! set the par vector with parameters to be fitted
  subroutine init_param(allpar,par)
    real(dp),dimension(:),intent(in)::allpar
    real(dp),dimension(:),allocatable,intent(out)::par
    integer::j,cnt
    logical,dimension(:),allocatable::done
    
    allocate(done(nfit))
    done=.false.
    if(.not.allocated(par)) allocate(par(nfit))
    cnt=0
    do j=1,npar
      if(tofit(j).eq.1)then
        cnt=cnt+1

        if(.not.done(cnt))then
          
          ! m(3) ==> mp/Ms
          if(id(cnt).eq.3)then
            par(cnt)=allpar(j)/MR_star(1,1)
            done(cnt)=.true.
          
          ! e(6),w(7)==>(ecosw,esinw)
          else if(id(cnt).eq.6)then
            par(cnt)=allpar(j)
            done(cnt)=.true.
            if(cnt.lt.nfit)then
              if(id(cnt+1).eq.7)then
                par(cnt)=allpar(j)*cos(allpar(j+1)*deg2rad)
                par(cnt+1)=allpar(j)*sin(allpar(j+1)*deg2rad)
                done(cnt)=.true.
                done(cnt+1)=.true.
              end if
            end if

          ! mA(8) ==> lambda=mA(j)+w(j-1)+lN(j+2)
          else if(id(cnt).eq.8)then
            par(cnt)=mod(mod(allpar(j)+allpar(j-1)+allpar(j+2),360._dp)+360._dp,360._dp)
            done(cnt)=.true.
          
          ! i(9),lN(10)==>(icoslN,isinlN)
          else if(id(cnt).eq.9)then
            par(cnt)=allpar(j)
            done(cnt)=.true.
            if(cnt.lt.nfit)then
              if(id(cnt+1).eq.10)then
                par(cnt)=allpar(j)*cos(allpar(j+1)*deg2rad)
                par(cnt+1)=allpar(j)*sin(allpar(j+1)*deg2rad)
                done(cnt)=.true.
                done(cnt+1)=.true.
              end if
            end if
          
          end if
        
        end if
        
        if(.not.done(cnt))then
          par(cnt)=allpar(j)
          done(cnt)=.true.
        end if
!         write(*,'("id = ",i3," => par(cnt=",i2,") = ",f23.8," = allpar(j=",i3,") = ",f23.8)')&
!           &id(cnt),cnt,par(cnt),j,allpar(j)
!         flush(6)
        
      end if ! tofit(j)
    end do
    
    deallocate(done)
    
    
!     call set_minmax()
    
    return
  end subroutine init_param
   
  ! from keplerian orbital elements to the parameters needed by L-M
  ! calls some previous subroutines
  subroutine set_par(m,R,P,a,e,w,mA,inc,lN,allpar,par)
    real(dp),dimension(:),intent(in)::m,R,P,a,e,w,mA,inc,lN
    real(dp),dimension(:),allocatable,intent(out)::allpar,par

    call sma_boundaries(R,a,amin,amax) ! IT SETS SEMI-MAJOR AXIS BOUNDS
    call set_all_param(m,R,P,e,w,mA,inc,lN,allpar) ! IT DEFINES THE VECTOR ALLPAR FROM ORBITAL PARAMETERS
    call init_param(allpar,par) ! IT DEFINES THE PARAMETERS PAR TO BE FITTED

    return
  end subroutine set_par
  ! ------------------------------------------------------------------ !
  
  ! fix the system_parameters in case a parameter has been read with a value not in [par_min, par_max]

  subroutine fix_all_parameters(all_parameters)
    real(dp),dimension(:)::all_parameters
    integer::i
    
    do i=1,npar
      if( tofit(i).eq. 1)then
        if( (all_parameters(i).lt.par_min(i)) .or. (all_parameters(i).gt.par_max(i)) )then
          all_parameters(i) = par_min(i)
        end if
      end if
    end do
  
    return
  end subroutine fix_all_parameters


  subroutine par2kel_fit(all_parameters,fit_parameters,m,R,P,a,e,w,mA,inc,lN,checkpar,wrt_info)
    use celestial_mechanics, only: semax
    real(dp),dimension(:),intent(in)::all_parameters,fit_parameters
    real(dp),dimension(:),intent(out)::m,R,P,a,e,w,mA,inc,lN
    real(dp),dimension(:),allocatable::atemp
    real(dp),dimension(8)::temp_kel ! ==> m,R,P,e,w,mA,inc,lN
    integer::j1,cnt
    logical,intent(out)::checkpar
    logical,intent(in),optional::wrt_info
    real(dp)::temp2
    real(dp),parameter::circ=360._dp
    
!     checkpar=.true.
    m=zero
    R=zero
    P=zero
    a=zero
    e=zero
    w=zero
    mA=zero
    inc=zero
    inc(2:)=90._dp
    lN=zero
    allocate(atemp(npar))
    atemp=all_parameters
    
    if(nfit.gt.0)then
      ! update atemp with fit_parameters
      do j1=1,nfit
          atemp(idall(j1))=fit_parameters(j1)
      end do
    end if
    
    m(1)=atemp(1)
    R(1)=atemp(2)
    cnt=2
    do cnt=2,NB
      j1=3+((cnt-2)*8) ! first parameter id == 3 <-> Mass body 2, ..., 10 <-> lN body 2, 11 <-> Mass body 3, ...
      
      temp_kel=atemp(j1:j1+7)
      
      if(tofit(j1).eq.1) temp_kel(1)=atemp(j1)*MR_star(1,1) ! mp/Ms to mp
      
      if(tofit(j1+3).eq.1.and.tofit(j1+4).eq.1)then ! ecosw,esinw to e,w
        temp2= sqrt((atemp(j1+3)*atemp(j1+3))+ (atemp(j1+4)*atemp(j1+4)))
        if(temp2.lt.e_bounds(1,cnt).or.temp2.gt.e_bounds(2,cnt))then
          if(present(wrt_info).and.wrt_info)then
            write(*,'(1x,a15,es23.16)')'ecc = ',temp2
            write(*,'(1x,2(a15,es23.16))')'emin = ',e_bounds(1,cnt),' emax = ',e_bounds(2,cnt)
            flush(6)
          end if
          checkpar=.false.
!           return
        end if
        temp_kel(4)=temp2
        temp_kel(5)=mod(rad2deg*atan2(atemp(j1+4),atemp(j1+3))+circ,circ)
        if(temp_kel(4).le.TOLERANCE) temp_kel(5)=90._dp
        if(abs(temp_kel(4)-atemp(j1+3)).le.TOLERANCE)then
          if(abs(atemp(j1+4)).gt.TOLERANCE) checkpar=.false.
        else if(abs(temp_kel(4)-atemp(j1+4)).le.TOLERANCE)then
          if(abs(atemp(j1+3)).gt.TOLERANCE) checkpar=.false.
        end if
      end if
      
      if(tofit(j1+6).eq.1)then !inc fit
        if(tofit(j1+7).eq.1)then
         temp_kel(7)=sqrt(atemp(j1+6)*atemp(j1+6)+atemp(j1+7)*atemp(j1+7)) ! icoslN,isinlN to inc,lN
         temp_kel(8)=mod(rad2deg*atan2(atemp(j1+7),atemp(j1+6))+circ,circ)
        end if
        if(temp_kel(7).le.zero.or.temp_kel(7).ge.180._dp)then
          if(present(wrt_info).and.wrt_info)write(*,'(1x,a15,es23.16)')'inc = ',temp2
          flush(6)
          checkpar=.false.
!           return
        end if
        
      end if
      
      if(tofit(j1+5).eq.1)temp_kel(6)=mod(mod((atemp(j1+5)-temp_kel(5)-temp_kel(8)),circ)+circ,circ) ! lambda to mA
      
      m(cnt)   =temp_kel(1)
      R(cnt)   =temp_kel(2)
      P(cnt)   =temp_kel(3)
      a(cnt)   =semax(m(1),temp_kel(1),temp_kel(3))
      e(cnt)   =temp_kel(4)
      w(cnt)   =temp_kel(5)
      mA(cnt)  =temp_kel(6)
      inc(cnt) =temp_kel(7)
      lN(cnt)  =temp_kel(8)
      
      if(present(wrt_info).and.wrt_info)then
        if(cnt.eq.2) write(*,'(1x,a15,3(1x,a23))')'name','min','val','max'
        write(*,'(1x,a15,3(1x,es23.16))')'m'//trim(string(cnt)), par_min(j1),  m(cnt),  par_max(j1)
        write(*,'(1x,a15,3(1x,es23.16))')'R'//trim(string(cnt)), par_min(j1+1),R(cnt),  par_max(j1+1)
        write(*,'(1x,a15,3(1x,es23.16))')'P'//trim(string(cnt)), par_min(j1+2),P(cnt),  par_max(j1+2)
!         write(*,'(1x,a15,3(1x,es23.16))')'a'//string(cnt), par_min(j1+1),a(cnt),  par_max(j1)
        write(*,'(1x,a15,3(1x,es23.16))')'e'//trim(string(cnt)), par_min(j1+3),e(cnt),  par_max(j1+3)
        write(*,'(1x,a15,3(1x,es23.16))')'w'//trim(string(cnt)), par_min(j1+4),w(cnt),  par_max(j1+4)
        write(*,'(1x,a15,3(1x,es23.16))')'mA'//trim(string(cnt)),par_min(j1+5),mA(cnt), par_max(j1+5)
        write(*,'(1x,a15,3(1x,es23.16))')'i'//trim(string(cnt)), par_min(j1+6),inc(cnt),par_max(j1+6)
        write(*,'(1x,a15,3(1x,es23.16))')'lN'//trim(string(cnt)),par_min(j1+7),lN(cnt), par_max(j1+7)
        flush(6)
      end if
      
    end do
    
    deallocate(atemp)

    return
  end subroutine par2kel_fit
  
!   check the physical bounds of the fitted parameters
  function checkbounds_fit(fit_parameters,wrt_info) result(check)
    real(dp),dimension(:),intent(in)::fit_parameters
    logical,intent(in),optional::wrt_info
    logical::check
    integer::j
    
    check = .true.

    if(nfit.gt.0)then
      do j=1,nfit
  !      body=int((idall(j)-3)/8)+2

        ! check if fit_parameters are within boundaries
        if(fit_parameters(j).lt.minpar(j).or.fit_parameters(j).gt.maxpar(j))then
          if(present(wrt_info))then
            if(wrt_info)then
              if(fit_parameters(j).lt.minpar(j)) write(*,'(a12,"( ",i3," )", F20.10, a, F20.10)')parid(j),j,fit_parameters(j),' < ',minpar(j)
              if(fit_parameters(j).gt.maxpar(j)) write(*,'(a12,"( ",i3," )", F20.10, a, F20.10)')parid(j),j,fit_parameters(j),' > ',maxpar(j)
            end if
          end if
          check=.false.
          return
        end if
      end do
    end if
    
    return
  end function checkbounds_fit
  

  function checkbounds_kel(m,R,P,e,w,mA,inc,lN) result(check)
    logical::check
    real(dp),dimension(:),intent(in)::m,R,P,e,w,mA,inc,lN
    
!     real(dp),dimension(:),allocatable::all_parameters
    real(dp)::temp
    integer::j,j1
    
    check=.true.
!     call set_all_param(m,R,P,e,w,mA,i,lN,all_parameters)
    
!     write(*,*)'Ms'
    if(m(1).lt.par_min(1).or.m(1).gt.par_max(1))then
      check=.false.
      return
    end if
    
!     write(*,*)'Rs'
    if(R(1).lt.par_min(2).or.R(1).gt.par_max(2))then
      check=.false.
      return
    end if
    
    do j=2,NB
      j1=(j-2)*8
      
!       write(*,*)'Mp',j,': ',m(j),par_min(j1+3),par_max(j1+3)
      ! mass
      if(m(j).lt.par_min(j1+3).or.m(j).gt.par_max(j1+3))then
        check=.false.
        return
      end if
      
!       write(*,*)'R',j,': ',R(j),par_min(j1+4),par_max(j1+4)
      ! radius
      if(R(j).lt.par_min(j1+4).or.R(j).gt.par_max(j1+4))then
        check=.false.
        return
      end if
      
!       write(*,*)'P',j,': ',P(j),par_min(j1+5),par_max(j1+5)
      ! period
      if(P(j).lt.par_min(j1+5).or.P(j).gt.par_max(j1+5))then
        check=.false.
        return
      end if
      
!       write(*,*)'e',j,': ',e(j),par_min(j1+6),par_max(j1+6)
      ! eccentricity
      if(e(j).lt.par_min(j1+6).or.e(j).gt.par_max(j1+6))then
        check=.false.
        return
      end if

!       write(*,*)'w',j,': ',w(j),par_min(j1+7),par_max(j1+7)
      ! argument of pericentre
      temp=w(j)-360._dp
      if(w(j).lt.par_min(j1+7).or.w(j).gt.par_max(j1+7))then
        check=.false.
!         temp=w(j)-360._dp
        if(temp.ge.par_min(j1+7).and.temp.le.par_max(j1+7))then
          check=.true.
        end if
        if(.not.check) return
      end if
!       write(*,*)'w = ',w(j),' w-360 = ',temp,' [min=',par_min(j1+7),' , max=', par_max(j1+7),'] ==> check = ',check

!       write(*,*)'mA',j,': ',mA(j),par_min(j1+8),par_max(j1+8)
      ! mean anomaly
      temp=mA(j)-360._dp
      if(mA(j).lt.par_min(j1+8).or.ma(j).gt.par_max(j1+8))then
        check=.false.
!         temp=mA(j)-360._dp
        if(temp.ge.par_min(j1+8).and.temp.le.par_max(j1+8))then
          check=.true.
        end if
        if(.not.check) return
      end if
!       write(*,*)'mA = ',mA(j),' mA-360 = ',temp,' [min=',par_min(j1+8),' , max=', par_max(j1+8),'] ==> check = ',check
      
!       write(*,*)'inc',j,': ',inc(j),par_min(j1+9),par_max(j1+9)
      ! inclination
      if(inc(j).lt.par_min(j1+9).or.inc(j).gt.par_max(j1+9))then
        check=.false.
        return
      end if
      
!       write(*,*)'lN',j,': ',lN(j),par_min(j1+10),par_max(j1+10)
      ! longitude of node
      temp=lN(j)-360._dp
      if(lN(j).lt.par_min(j1+10).or.lN(j).gt.par_max(j1+10))then
        check=.false.
!         temp=lN(j)-360._dp
        if(temp.ge.par_min(j1+10).and.temp.le.par_max(j1+10))then
          check=.true.
        end if
        if(.not.check) return
      end if
!       write(*,*)'lN = ',lN(j),' lN-360 = ',temp,' [min=',par_min(j1+10),' , max=', par_max(j1+10),'] ==> check = ',check

      
    end do
  
    return
  end function checkbounds_kel

    ! ------------------------------------------------------------------ !
  subroutine convert_parameters(allpar,par,m,R,P,a,e,w,mA,inc,lN,checkpar,wrt_info)
    real(dp),dimension(:),intent(in)::allpar,par
    real(dp),dimension(:),intent(out)::m,R,P,a,e,w,mA,inc,lN
    logical,intent(inout)::checkpar
    logical,intent(in),optional::wrt_info
    
    
!     checkpar=.true.
    checkpar=checkbounds_fit(par,wrt_info)
    if(present(wrt_info).and.wrt_info)then
      if(wrt_info) write(*,*)' checkbounds_fit: ',checkpar
      flush(6)
    end if
    if(checkpar)then
      call par2kel_fit(allpar,par,m,R,P,a,e,w,mA,inc,lN,checkpar,wrt_info)
      if(present(wrt_info))then
        if(wrt_info) write(*,*)' par2kel_fit: ',checkpar
        flush(6)
      end if
    end if
    if(checkpar)then
      checkpar=checkbounds_kel(m,R,P,e,w,mA,inc,lN)
      if(present(wrt_info))then
        if(wrt_info) write(*,*)' checkbounds_kel: ',checkpar
        flush(6)
      end if
    end if

    
    return
  end subroutine convert_parameters
  ! ------------------------------------------------------------------ !
  
  ! subroutine needed to update the all_parameters from the fitted fit_parameters
  subroutine update_parameters_fit2all(fit_parameters,all_parameters)
    real(dp),dimension(:),intent(in)::fit_parameters
    real(dp),dimension(:),intent(inout)::all_parameters
    
    real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,inc,lN
    real(dp),dimension(:),allocatable::all_updating
    logical::check
    
    check=.true.
    allocate(m(NB),R(NB),P(NB),a(NB),e(NB),w(NB),mA(NB),inc(NB),lN(NB))
    call par2kel_fit(all_parameters,fit_parameters,m,R,P,a,e,w,mA,inc,lN,check)
    call set_all_param(m,R,P,e,w,mA,inc,lN,all_updating)
    all_parameters=all_updating
    deallocate(all_updating,m,R,P,a,e,w,mA,inc,lN)

    return
  end subroutine update_parameters_fit2all
  
  ! function for pso/pik to initialize properly the first-generation population
  function check_only_boundaries(all_parameters,fit_parameters,wrt_info) result(check)
    logical::check
    real(dp),dimension(:),intent(in)::all_parameters,fit_parameters
    logical,intent(in),optional::wrt_info
    real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,inc,lN
    
    check=.true.
    check=checkbounds_fit(fit_parameters)
!     write(*,'(a,l2)')'checkbounds_fit = ',check
    if(check)then
      allocate(m(NB),R(NB),P(NB),a(NB),e(NB),w(NB),mA(NB),inc(NB),lN(NB))
      if(present(wrt_info).and.wrt_info)then
        call par2kel_fit(all_parameters,fit_parameters,m,R,P,a,e,w,mA,inc,lN,check,wrt_info)
      else
        call par2kel_fit(all_parameters,fit_parameters,m,R,P,a,e,w,mA,inc,lN,check)
      end if
!       write(*,'(a,l2)')'par2kel_fit = ',check
!       if(check) check=checkbounds_kel(m,R,P,e,w,mA,inc,lN)
        if(check)then
          check=checkbounds_kel(m,R,P,e,w,mA,inc,lN)
!           write(*,'(a,l2)')'checkbounds_kel = ',check
        end if
      deallocate(m,R,P,a,e,w,mA,inc,lN)
    end if
      
    return
  end function check_only_boundaries
  
  
  ! boundaries check that output a scale parameter and not a logical .true./.false.
  ! the scale parameter will be multiplied to the weighted residuals to worsening the fitness
  !
  
  function get_fit_scale(value,minvalue,maxvalue) result(xscale)
    real(dp)::xscale
    real(dp),intent(in)::value,minvalue,maxvalue
    real(dp)::delta
    
    xscale=zero
    delta=abs(maxvalue-minvalue)
    
  ! check if fit_parameters are within boundaries
    if(value.lt.minvalue)then
      xscale=abs((value-minvalue)/delta)*100._dp
    else if(value.gt.maxvalue)then
      xscale=abs((value-maxvalue)/delta)*100._dp
    else
      xscale=zero
    end if
    
  
    return
  end function get_fit_scale
  
  function checkbounds_fit_scale(fit_parameters) result(bound_scale)
    real(dp)::bound_scale
    real(dp),dimension(:),intent(in)::fit_parameters
    integer::j
    real(dp)::xscale
    
    bound_scale=zero
    xscale=zero
    
    do j=1,nfit
      ! check if fit_parameters are within boundaries
      xscale=get_fit_scale(fit_parameters(j),minpar(j),maxpar(j))
      bound_scale=bound_scale+xscale
    end do
    
    return
  end function checkbounds_fit_scale
  
    function checkbounds_kel_scale(m,R,P,e,w,mA,inc,lN) result(kel_scale)
    real(dp)::kel_scale
    real(dp),dimension(:),intent(in)::m,R,P,e,w,mA,inc,lN
    
    real(dp)::temp,xscale,xtemp1,xtemp2
    integer::j,j1
    
    kel_scale=zero
    xscale=zero

!     write(*,*)'Ms'
    xscale=get_fit_scale(m(1),par_min(1),par_max(1))
    kel_scale=kel_scale+xscale
    
!     write(*,*)'Rs'
    xscale=get_fit_scale(R(1),par_min(2),par_max(2))
    kel_scale=kel_scale+xscale
    
    do j=2,NB
      j1=(j-2)*8
      
      ! mass
      xscale=get_fit_scale(m(j),par_min(j1+3),par_max(j1+3))
      kel_scale=kel_scale+xscale
      
      ! radius
      xscale=get_fit_scale(R(j),par_min(j1+4),par_max(j1+4))
      kel_scale=kel_scale+xscale
      
      ! period
      xscale=get_fit_scale(P(j),par_min(j1+5),par_max(j1+5))
      kel_scale=kel_scale+xscale
      
      ! eccentricity
      xscale=get_fit_scale(e(j),par_min(j1+6),par_max(j1+6))
      kel_scale=kel_scale+xscale

      ! argument of pericentre
      temp=w(j)-360._dp
      xtemp1=get_fit_scale(temp,par_min(j1+7),par_max(j1+7))
      xtemp2=get_fit_scale(w(j),par_min(j1+7),par_max(j1+7))
      if(xtemp1.gt.zero.and.xtemp2.gt.zero)then
        xscale=max(xtemp1,xtemp2)
      else
        xscale=zero
      end if
      kel_scale=kel_scale+xscale
      
      ! mean anomaly
      temp=mA(j)-360._dp
      xtemp1=get_fit_scale(temp,par_min(j1+8),par_max(j1+8))
      xtemp2=get_fit_scale(mA(j),par_min(j1+8),par_max(j1+8))
      if(xtemp1.gt.zero.and.xtemp2.gt.zero)then
        xscale=max(xtemp1,xtemp2)
      else
        xscale=zero
      end if
      kel_scale=kel_scale+xscale
      
      ! inclination
      xscale=get_fit_scale(inc(j),par_min(j1+9),par_max(j1+9))
      kel_scale=kel_scale+xscale
      
      ! longitude of node
      temp=lN(j)-360._dp
      xtemp1=get_fit_scale(temp,par_min(j1+10),par_max(j1+10))
      xtemp2=get_fit_scale(lN(j),par_min(j1+10),par_max(j1+10))
      if(xtemp1.gt.zero.and.xtemp2.gt.zero)then
        xscale=max(xtemp1,xtemp2)
      else
        xscale=zero
      end if
      kel_scale=kel_scale+xscale
      
    end do
  
    return
  end function checkbounds_kel_scale
  
  subroutine convert_parameters_scale(all_parameters,fit_parameters,m,R,P,a,e,w,mA,inc,lN,fit_scale)
    real(dp),dimension(:),intent(in)::all_parameters,fit_parameters
    real(dp),dimension(:),intent(out)::m,R,P,a,e,w,mA,inc,lN
    real(dp),intent(out)::fit_scale
    real(dp)::bound_scale,kel_scale
    logical::checkpar
    
    checkpar=.true.
    fit_scale=zero
    bound_scale=checkbounds_fit_scale(fit_parameters)
    call par2kel_fit(all_parameters,fit_parameters,m,R,P,a,e,w,mA,inc,lN,checkpar)
    kel_scale=checkbounds_kel_scale(m,R,P,e,w,mA,inc,lN)
    fit_scale=one+bound_scale+kel_scale ! min value of fit_scale has to be 1
    
    
    return
  end subroutine convert_parameters_scale
  
  function check_only_boundaries_scale(all_parameters,fit_parameters) result(fit_scale)
    real(dp)::fit_scale
    real(dp),dimension(:),intent(in)::all_parameters,fit_parameters
    
    real(dp),dimension(:),allocatable::m,R,P,a,e,w,mA,inc,lN
    real(dp)::bound_scale,kel_scale
    logical::checkpar
    
    checkpar=.true.
    bound_scale=checkbounds_fit_scale(fit_parameters)
!     write(*,'(a,es23.16)')'checkbounds_fit = ',bound_scale
    allocate(m(NB),R(NB),P(NB),a(NB),e(NB),w(NB),mA(NB),inc(NB),lN(NB))
    call par2kel_fit(all_parameters,fit_parameters,m,R,P,a,e,w,mA,inc,lN,checkpar)
!     write(*,'(a,l2)')'par2kel_fit = ',checkpar
    kel_scale=checkbounds_kel_scale(m,R,P,e,w,mA,inc,lN)
    deallocate(m,R,P,a,e,w,mA,inc,lN)
!     write(*,'(a,es23.16)')'checkbounds_kel = ',kel_scale
    fit_scale=one+bound_scale+kel_scale
!     write(*,'(a,es23.16)')'fit_scale = ',fit_scale
      
    return
  end function check_only_boundaries_scale
    ! ------------------------------------------------------------------ !
  ! given the computed parameters with the L-M it adjusts some parameters
  ! i.e. setting all the angles between 0 and 360 (no < 0 angles) etc.
  subroutine param_adj(par,sigpar)
    real(dp),dimension(:),intent(inout)::par,sigpar
    integer::j1
    real(dp),parameter::circ=360._dp,hcirc=180._dp

    do j1=1,nfit

      if(parid(j1)(1:1).eq.'w'.or.parid(j1)(1:2).eq.'mA'.or.parid(j1)(1:2).eq.'lN'.or.parid(j1)(1:2).eq.'la')then
        par(j1)=mod(par(j1),circ)
        sigpar(j1)=mod(sigpar(j1),circ)
      end if
      
    end do

    return
  end subroutine param_adj
  ! ------------------------------------------------------------------ !
  
  ! ------------------------------------------------------------------ !
  ! conversion from parameter boundaries [0, 1] --> [ParMin, ParMax]
  subroutine norm2par(norm,par)
    real(dp),dimension(:),intent(in)::norm
    real(dp),dimension(:),intent(out)::par
    real(dp)::dpar
    integer::j

    do j=1,nfit
      dpar=abs(maxpar(j)-minpar(j))
      par(j)=minpar(j)+dpar*norm(j)
    end do

    return
  end subroutine norm2par

!   ! conversion from parameter boundaries [0, 1] --> [ParMin, ParMax]
!   subroutine norm2par_2(norm,par,allpar)
!     real(dp),dimension(:),intent(in)::norm
!     real(dp),dimension(:),intent(out)::par
! !     real(dp),dimension(:),intent(out)::allpar
!     real(dp),dimension(:)::allpar
!     real(dp)::dpar
!     integer::j,j1
! 
!     j1=0
!     do j=1,npar
!       if(tofit(j).eq.1)then
!         j1=j1+1
!         dpar=abs(par_max(j)-par_min(j))
!         par(j1)=par_min(j)+dpar*norm(j1)
!         allpar(j)=par(j1)
!       end if
!     end do
! 
!     return
!   end subroutine norm2par_2



end module parameters_conversion
