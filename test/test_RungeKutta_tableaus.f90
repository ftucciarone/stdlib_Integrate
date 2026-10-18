program test_RungeKutta_tableaus
    !! Algebraic checks of every method of the library, independent of
    !! any test problem: A strictly lower triangular, row sums c(i) =
    !! sum(A(i,:)), FSAL rows, and the Butcher order conditions (all rooted
    !! trees up to order 5): satisfied up to the declared order, and not all
    !! satisfied at the next one (so the declared order is not understated).
    !! The same for the embedded weights b_hat.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: rk_tableau_t, rk_tableau, rk_method_names
    implicit none

    integer, parameter :: MAX_ORDER = 5
    real(wp), parameter :: tol = 1.0e-13_wp
    character(len=:), allocatable :: names(:)
    type(rk_tableau_t) :: tab
    logical :: found, all_passed, method_passed
    integer :: m, i

    names = rk_method_names()
    all_passed = .true.

    do m = 1, size(names)
        tab = rk_tableau(trim(names(m)), found)
        method_passed = .true.
        call expect(found, "found by name")
        if (.not. found) cycle
        call expect(tab%name == trim(names(m)), "name matches")
        call expect(size(tab%b) == tab%stages .and. size(tab%c) == tab%stages, "sizes")
        do i = 1, tab%stages
            call expect(all(tab%A(i, i:) == 0.0_wp), "A strictly lower triangular")
            call expect(abs(tab%c(i) - sum(tab%A(i, :))) < tol, "c(i) = sum(A(i,:))")
        end do
        if (tab%fsal) call expect(all(abs(tab%A(tab%stages, :) - tab%b) < tol), "FSAL: last row of A = b")
        call check_order(tab%b, tab%order, "b")
        call expect(tab%linear_order >= tab%order, "linear_order >= order")
        if (allocated(tab%b_hat)) call check_order(tab%b_hat, tab%order_hat, "b_hat")
        if (method_passed) print '(a10, a, i0, a)', trim(tab%name), ": OK (order ", tab%order, ")"
    end do

    tab = rk_tableau("no_such_method", found)
    call expect(.not. found .and. tab%stages == 0, "unknown name not found")

    if (.not. all_passed) error stop "tableau checks FAILED"
    print *, "All tableau checks PASSED"

contains

    subroutine expect(condition, what)
        logical, intent(in) :: condition
        character(*), intent(in) :: what

        if (.not. condition) then
            print '(3a)', trim(tab%name), ": FAILED: ", what
            all_passed = .false.
            method_passed = .false.
        end if
    end subroutine expect

    ! Orders 1..order satisfied by weights w; order+1 not (if it is checked here).
    subroutine check_order(w, order, label)
        real(wp), intent(in) :: w(:)
        integer, intent(in) :: order
        character(*), intent(in) :: label
        integer :: p

        call expect(order >= 1 .and. order <= MAX_ORDER, label // ": order within 1..5")
        do p = 1, min(order, MAX_ORDER)
            call expect(maxval(abs(residuals(tab%A, w, tab%c, p))) < tol, &
                        label // ": order condition " // achar(iachar('0') + p) // " satisfied")
        end do
        if (order < MAX_ORDER) call expect(maxval(abs(residuals(tab%A, w, tab%c, order + 1))) > tol, &
                                           label // ": next order not satisfied (order understated?)")
    end subroutine check_order

    ! Residuals of the order conditions for the rooted trees of order p:
    ! sum(w * Phi(tree)) - 1/gamma(tree).
    function residuals(a, w, c, p) result(r)
        real(wp), intent(in) :: a(:,:), w(:), c(:)
        integer, intent(in) :: p
        real(wp), allocatable :: r(:)
        real(wp) :: a_c(size(c)), a_c2(size(c)), a_a_c(size(c))   ! matmul(a, c), matmul(a, c**2), matmul(a, a_c)

        a_c = matmul(a, c)
        a_c2 = matmul(a, c**2)
        a_a_c = matmul(a, a_c)
        select case (p)
        case (1)
            r = [sum(w) - 1]
        case (2)
            r = [sum(w*c) - 1/2.0_wp]
        case (3)
            r = [sum(w*c**2) - 1/3.0_wp, sum(w*a_c) - 1/6.0_wp]
        case (4)
            r = [sum(w*c**3) - 1/4.0_wp, sum(w*c*a_c) - 1/8.0_wp, &
                 sum(w*a_c2) - 1/12.0_wp, sum(w*a_a_c) - 1/24.0_wp]
        case (5)
            r = [sum(w*c**4) - 1/5.0_wp, sum(w*c**2*a_c) - 1/10.0_wp, &
                 sum(w*c*a_c2) - 1/15.0_wp, sum(w*c*a_a_c) - 1/30.0_wp, &
                 sum(w*a_c**2) - 1/20.0_wp, sum(w*matmul(a, c**3)) - 1/20.0_wp, &
                 sum(w*matmul(a, c*a_c)) - 1/40.0_wp, sum(w*matmul(a, a_c2)) - 1/60.0_wp, &
                 sum(w*matmul(a, a_a_c)) - 1/120.0_wp]
        end select
    end function residuals

end program test_RungeKutta_tableaus
