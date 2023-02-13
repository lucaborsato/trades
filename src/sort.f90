module swapping

    interface swap
        module procedure swap_i, swap_r, swap_rv
    end interface

contains

    ! ************************* !
    ! swapping subroutines by Numerical Recipes
    !
    subroutine swap_i(a, b)
        integer, intent(INOUT)::a, b
        integer::dum
        dum = a
        a = b
        b = dum
    end subroutine swap_i
    !BL
    subroutine swap_r(a, b)
        use constants, only: dp
        real(dp), intent(INOUT)::a, b
        real(dp)::dum
        dum = a
        a = b
        b = dum
    end subroutine swap_r
    !BL
    subroutine swap_rv(a, b)
        use constants, only: dp
        real(dp), dimension(:), intent(INOUT)::a, b
        real(dp), dimension(size(a))::dum
        dum = a
        a = b
        b = dum
    end subroutine swap_rv

end module swapping

module sorting

    interface sort
        module procedure sort1, sort2, sort3
    end interface

contains

    ! ************************* !
    ! sorting algorithm by Numerical Recipes
    !

    ! subroutine to make indexes vector (idx) with the position of arr to sort in the right order....arr will be updated
    subroutine indexx(arr, idx)
        use constants, only: dp
        use swapping
        implicit none
        real(dp), dimension(:), intent(in)::arr
        integer, dimension(:), intent(out)::idx
        integer, parameter::NN = 15, NSTACK = 50, NPAR_ARTH = 16, NPAR2_ARTH = 8
        real(dp)::a
        integer::n, k, i, j, indext, jstack, l, r
    !! integer::ii,temp
        integer, dimension(NSTACK)::istack

        n = size(arr)
        do i = 1, n
            idx(i) = i
        end do
        jstack = 0
        l = 1
        r = n
        bigdo: do
            if (r-l < NN) then
                jdo: do j = l+1, r
                    indext = idx(j)
                    a = arr(indext)
                    ido: do i = j-1, l, -1
                        if (arr(idx(i)) <= a) exit ido
                        idx(i+1) = idx(i)
                    end do ido
                    idx(i+1) = indext
                end do jdo
                if (jstack == 0) return
                r = istack(jstack)
                l = istack(jstack-1)
                jstack = jstack-2
            else
                k = (l+r)/2
                call swap(idx(k), idx(l+1))
                call icomp_xchg(idx(l), idx(r))
                call icomp_xchg(idx(l+1), idx(r))
                call icomp_xchg(idx(l), idx(l+1))
                i = l+1
                j = r
                indext = idx(l+1)
                a = arr(indext)
                indo: do
                    ado: do
                        i = i+1
                        if (arr(idx(i)) >= a) exit ado
                    end do ado
                    bdo: do
                        j = j-1
                        if (arr(idx(j)) <= a) exit bdo
                    end do bdo
                    if (j < i) exit indo
                    call swap(idx(i), idx(j))
                end do indo
                idx(l+1) = idx(j)
                idx(j) = indext
                jstack = jstack+2
                if (jstack > NSTACK) then
                    write (*, *) "indexx: NSTACK too small"
                    return
                end if
                if (r-i+1 >= j-l) then
                    istack(jstack) = r
                    istack(jstack-1) = i
                    r = j-1
                else
                    istack(jstack) = j-1
                    istack(jstack-1) = l
                    l = i
                end if
            end if
        end do bigdo

    contains

        subroutine icomp_xchg(i, j)
            integer, intent(inout)::i, j
            integer::swp
            if (arr(j) < arr(i)) then
                swp = i
                i = j
                j = swp
            end if
        end subroutine icomp_xchg

    end subroutine indexx

    subroutine sort1(arr)
        use constants, only: dp
        implicit none
        real(dp), dimension(:), intent(inout)::arr
        integer, dimension(size(arr))::idx
        call indexx(arr, idx)
        arr = arr(idx)
        return
    end subroutine sort1

    subroutine sort2(arr1, arr2)
        use constants, only: dp
        implicit none
        real(dp), dimension(:), intent(inout)::arr1, arr2
        integer, dimension(size(arr1))::idx
        call indexx(arr1, idx)
        arr1 = arr1(idx)
        arr2 = arr2(idx)
        return
    end subroutine sort2

    subroutine sort3(arr1, arr2, arr3)
        use constants, only: dp
        implicit none
        real(dp), dimension(:), intent(inout)::arr1, arr2, arr3
        integer, dimension(size(arr1))::idx
        call indexx(arr1, idx)
        arr1 = arr1(idx)
        arr2 = arr2(idx)
        arr3 = arr3(idx)
        return
    end subroutine sort3

end module sorting
