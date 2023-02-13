module timing
    use constants
!$  use omp_lib
    implicit none

contains

    subroutine get_cpu_time(t1)
        real(dp), intent(out)::t1

        t1 = zero
        call cpu_time(t1)
!$      t1 = omp_get_wtime()

        return
    end subroutine get_cpu_time

    ! ------------------------------------------------------------------ !
    ! computes the time in hours minutes seconds between t2 and t1
    subroutine timer(t1, t2, h, m, s)
        real(dp), intent(in) :: t1, t2
        integer, intent(out) :: h, m
        real(dp), intent(out)::s
        real(dp) :: tots
        integer :: totm

        tots = (t2-t1)/ncpu
        totm = int(tots/60._dp)
        h = int(totm/60._dp)
        m = totm-h*int(60)
        s = tots-totm*60._dp

        return
    end subroutine timer

    ! provided the start time t1, it calculate automaticaly the end time t2
    ! and print the summarized elapsed time
    ! it uses also the openMP time
    subroutine elapsed_time(t1)
        real(dp), intent(in)::t1

        real(dp)::t2
        integer::h, m
        real(dp)::s

        call get_cpu_time(t2)

        call timer(t1, t2, h, m, s)

        return
    end subroutine elapsed_time

    ! ------------------------------------------------------------------ !

end module timing
