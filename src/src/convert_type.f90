module convert_type
  use constants,only:dp
  implicit none

  interface string
    module procedure i2s,r2s
  end interface string

  contains

  !function to read an integer to a character/string
  function i2s(number) result(out)
    character(512)::out
    integer,intent(in)::number

    write(out,*) number
    out=trim(adjustl(out))

    return
  end function i2s

  !function to read an real to a character/string
  function r2s(xxx) result(out)
    character(512)::out
    real(dp),intent(IN)::xxx
    character(128)::fmt
    integer::prec,dprec
    character(128)::c_prec,c_dprec

    prec = precision(xxx)
    dprec = 2 * prec
    write(c_prec,'(i3)') prec
    c_prec = trim(adjustl(c_prec))
    write(c_dprec,'(i3)') dprec
    c_dprec = trim(adjustl(c_dprec))
    fmt = adjustl('(f'//trim(c_dprec)//'.'//trim(c_prec)//')')
    write(out,trim(fmt)) xxx
    out = trim(adjustl(out))

    return
  end function r2s

end module convert_type
