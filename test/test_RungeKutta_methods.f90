program test_RungeKutta_methods
    !! The schemes by name (rk_scheme): for every method of the library, the
    !! same result as the tableau engine (time_scheme_ERK_t) over 20 steps, bit
    !! for bit (for rk4, the hardcoded step against the engine on its
    !! tableau), the convergence order on dy/dt = -lambda y, and the tendency
    !! registers. Failures while building a scheme: test_RungeKutta_construction.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: status_t, time_scheme_c, rk_scheme, rk_tableau, rk_method_names, rk_tableau_t
    use stdlib_integrate_rk_steps, only: time_scheme_rk4_t
    use stdlib_integrate_rk_generic, only: time_scheme_ERK_t
    use test_decay_model, only: decay_model_t, decay_state_t
    implicit none

    real(wp), parameter :: lambda = 2.0_wp, t_end = 1.0_wp
    type(decay_model_t), target :: model
    character(len=:), allocatable :: names(:)
    class(time_scheme_c), allocatable :: scheme
    type(rk_tableau_t) :: tab
    logical :: all_passed
    integer :: m, n1
    real(wp) :: y_rk, y_erk, err1, err2, order

    model%lambda = lambda
    all_passed = .true.
    names = rk_method_names()

    do m = 1, size(names)
        tab = rk_tableau(trim(names(m)))

        ! Same as the tableau engine, bit for bit.
        call rk_scheme(scheme, trim(names(m)), 0.05_wp, model)
        y_rk = run(scheme, 20)
        call rk_scheme(scheme, tab, 0.05_wp, model)
        y_erk = run(scheme, 20)
        call expect(y_rk == y_erk, trim(names(m)) // ": same as time_scheme_ERK_t, bit for bit")

        ! Convergence order at the linear order (the test problem is linear).
        n1 = 50
        if (tab%linear_order >= 5) n1 = 20
        call rk_scheme(scheme, trim(names(m)), t_end/n1, model)
        err1 = abs(run(scheme, n1) - exp(-lambda*t_end))
        call rk_scheme(scheme, trim(names(m)), t_end/(2*n1), model)
        err2 = abs(run(scheme, 2*n1) - exp(-lambda*t_end))
        order = log(err1/err2) / log(2.0_wp)
        call expect(abs(order - tab%linear_order) < 0.15_wp, trim(names(m)) // ": convergence order")
        print '(a10, a, f6.3, a, i0, a)', trim(names(m)), ": order ", order, ", ", registers(scheme), " register(s)"
    end do

    ! Tendency registers: four for the hardcoded RK4, one per stage for the engine.
    call rk_scheme(scheme, "rk4", 0.1_wp, model)
    call expect(registers(scheme) == 4, "rk4: four tendency registers")
    select type (scheme)
    type is (time_scheme_rk4_t)
    class default
        call expect(.false., "rk4 by name: the hardcoded step")
    end select
    call rk_scheme(scheme, "rk3ws", 0.1_wp, model)
    call expect(registers(scheme) == 3, "rk3ws: the engine, three tendency registers")
    call rk_scheme(scheme, "dp54", 0.1_wp, model)
    call expect(registers(scheme) == 7, "dp54: the engine, seven tendency registers")

    if (.not. all_passed) error stop "scheme by name checks FAILED"
    print *, "All scheme by name checks PASSED"

contains

    subroutine expect(condition, what)
        logical, intent(in) :: condition
        character(*), intent(in) :: what

        if (.not. condition) then
            print '(2a)', "FAILED: ", what
            all_passed = .false.
        end if
    end subroutine expect

    ! y after n steps from y(0) = 1.
    real(wp) function run(scheme, n) result(y)
        class(time_scheme_c), intent(inout) :: scheme
        integer, intent(in) :: n
        type(decay_state_t) :: state
        type(status_t) :: status
        real(wp) :: time
        integer :: i

        state%y = 1.0_wp
        time = 0.0_wp
        do i = 1, n
            call scheme%advance(time, state, status)
            if (status%failed) error stop "step failed"
        end do
        y = state%y
    end function run

    integer function registers(scheme)
        class(time_scheme_c), intent(in) :: scheme

        registers = -1
        select type (scheme)
        type is (time_scheme_rk4_t)
            registers = size(scheme%k)
        type is (time_scheme_ERK_t)
            registers = size(scheme%k)
        end select
    end function registers

end program test_RungeKutta_methods
