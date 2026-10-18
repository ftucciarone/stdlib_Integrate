program test_RungeKutta_arrays
    !! The schemes on plain arrays (rk_array_scheme_t):
    !!   - for every method, by name and by tableau, the same steps as the
    !!     scheme built the same way on model classes, bit for bit, on a
    !!     nonlinear time-dependent system, and the same number of registers;
    !!   - the hardcoded RK4 and the tableau loop on RK4's tableau: the same
    !!     steps, bit for bit;
    !!   - convergence at each method's linear order on dy/dt = -lambda y, with
    !!     the right-hand side an internal procedure reading `lambda` from its
    !!     host;
    !!   - `boundaries` runs once per step plus once per stage state built, on
    !!     every state before it is evaluated, at its stage time;
    !!   - an unknown name, an invalid tableau, a size below 1, a scheme not
    !!     built, a state of the wrong size and a strided state are failures
    !!     returned in `status`.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: stage_t, status_t, time_scheme_c, rk_scheme, rk_method_names, &
                                   & rk_tableau, rk_tableau_t, rk_array_scheme_t
    use stdlib_integrate_rk_steps, only: time_scheme_rk4_t
    use stdlib_integrate_rk_generic, only: time_scheme_ERK_t
    use test_RungeKutta_arrays_model, only: array_model_t, array_state_t, rhs
    implicit none

    integer,  parameter :: n = 6, n_steps = 20
    real(wp), parameter :: dt = 0.05_wp

    character(len=8), allocatable :: names(:)
    real(wp) :: lambda, boundaries_at
    integer  :: n_boundaries, m
    logical  :: all_passed, boundaries_ok

    all_passed = .true.
    names = rk_method_names()

    do m = 1, size(names)
        call check_equivalence(trim(names(m)), by_tableau=.false.)
        call check_equivalence(trim(names(m)), by_tableau=.true.)
        call check_convergence(trim(names(m)))
        call check_boundaries(trim(names(m)))
    end do
    call check_rk4_paths()
    call check_failures()

    if (.not. all_passed) error stop "array scheme checks FAILED"
    print *, "All array scheme checks PASSED"

contains

    subroutine expect(condition, what)
        logical,      intent(in) :: condition
        character(*), intent(in) :: what

        if (.not. condition) then
            print '(2a)', "FAILED: ", what
            all_passed = .false.
        end if
    end subroutine expect

    function initial_state() result(y)
        real(wp) :: y(n)
        integer  :: i

        y = [(1.0_wp + 0.3_wp * i, i = 1, n)]
    end function initial_state

    !> n_steps steps of the array scheme and of the class scheme, built the
    !> same way (by name or by tableau), from the same state: the same values,
    !> bit for bit, the same time, and the same number of registers.
    subroutine check_equivalence(name, by_tableau)
        character(*), intent(in) :: name
        logical,      intent(in) :: by_tableau

        type(rk_array_scheme_t) :: array_scheme
        type(array_model_t), target :: model
        class(time_scheme_c), allocatable, target :: scheme
        type(array_state_t), target :: state
        real(wp) :: y(n), t_array, t_class
        integer  :: i, registers
        character(len=:), allocatable :: what

        model%n = n
        if (by_tableau) then
            call rk_scheme(array_scheme, rk_tableau(name), n, dt)
            call rk_scheme(scheme, rk_tableau(name), dt, model)
            what = name // " (by tableau)"
        else
            call rk_scheme(array_scheme, name, n, dt)
            call rk_scheme(scheme, name, dt, model)
            what = name // " (by name)   "
        end if
        registers = -1
        select type (scheme)
        type is (time_scheme_rk4_t)
            registers = size(scheme%k)
        type is (time_scheme_ERK_t)
            registers = size(scheme%k)
        end select
        call expect(registers == size(array_scheme%k, 2), what // ": the same registers as the class scheme")
        y = initial_state()
        state%y = initial_state()
        t_array = 0.3_wp
        t_class = 0.3_wp
        do i = 1, n_steps
            call array_scheme%advance(rhs, t_array, y)
            call scheme%advance(t_class, state)
        end do
        call expect(all(y == state%y) .and. t_array == t_class, what // ": same steps as the class scheme, bit for bit")
        if (any(y /= state%y)) print '(4x,a,es10.2)', "largest difference ", maxval(abs(y - state%y))
    end subroutine check_equivalence

    !> Observed order on dy/dt = -lambda y, the right-hand side an internal
    !> procedure that reads lambda from this program.
    subroutine check_convergence(name)
        character(*), intent(in) :: name

        type(rk_tableau_t) :: tableau
        real(wp) :: err_coarse, err_fine, observed

        lambda = 2.0_wp
        tableau = rk_tableau(name)
        err_coarse = abs(decay(name, 25) - exp(-lambda))
        err_fine   = abs(decay(name, 50) - exp(-lambda))
        observed = log(err_coarse / err_fine) / log(2.0_wp)
        print '(a10,a,f6.3)', name, ": observed order on arrays ", observed
        call expect(abs(observed - tableau%linear_order) < 0.15_wp, name // ": convergence order on arrays")
    end subroutine check_convergence

    !> y(1) after n_steps steps of dy/dt = -lambda y from y = 1 to t = 1.
    real(wp) function decay(name, n_steps) result(y_end)
        character(*), intent(in) :: name
        integer,      intent(in) :: n_steps

        type(rk_array_scheme_t) :: scheme
        real(wp) :: y(1), time
        integer  :: i

        call rk_scheme(scheme, name, 1, 1.0_wp / n_steps)
        y = 1.0_wp
        time = 0.0_wp
        do i = 1, n_steps
            call scheme%advance(decay_rhs, time, y)
        end do
        y_end = y(1)
    end function decay

    subroutine decay_rhs(t, y, dydt, stage)
        real(wp), intent(in)  :: t
        real(wp), intent(in)  :: y(:)
        real(wp), intent(out) :: dydt(:)
        type(stage_t), intent(in) :: stage

        dydt = -lambda * y
    end subroutine decay_rhs

    !> Calls to boundaries per step: one on the initial state, plus one per stage
    !> state built (stages after the first, but dp54's unused seventh).
    integer function expected_calls(name) result(calls)
        character(*), intent(in) :: name

        select case (name)
        case ("rk3ws");            calls = 3
        case ("rk4");              calls = 4
        case ("dp54");             calls = 6
        case default;              calls = -1
        end select
    end function expected_calls

    subroutine check_boundaries(name)
        character(*), intent(in) :: name

        type(rk_array_scheme_t) :: scheme
        real(wp) :: y(n), time
        integer  :: i

        call rk_scheme(scheme, name, n, dt)
        y = initial_state()
        time = 0.0_wp
        n_boundaries = 0
        boundaries_at = -1.0_wp
        boundaries_ok = .true.
        do i = 1, n_steps
            call scheme%advance(checked_rhs, time, y, boundaries=count_boundaries)
        end do
        call expect(n_boundaries == n_steps * expected_calls(name), name // ": calls to apply_boundaries per step")
        call expect(boundaries_ok, name // ": boundaries applied to every evaluated state, at its time")
    end subroutine check_boundaries

    subroutine count_boundaries(t, y, stage)
        real(wp), intent(in)    :: t
        real(wp), intent(inout) :: y(:)
        type(stage_t), intent(in) :: stage

        n_boundaries = n_boundaries + 1
        boundaries_at = t
    end subroutine count_boundaries

    subroutine checked_rhs(t, y, dydt, stage)
        real(wp), intent(in)  :: t
        real(wp), intent(in)  :: y(:)
        real(wp), intent(out) :: dydt(:)
        type(stage_t), intent(in) :: stage

        if (t /= boundaries_at) boundaries_ok = .false.
        call rhs(t, y, dydt, stage)
    end subroutine checked_rhs

    !> The hardcoded RK4 and the tableau loop on RK4's tableau: the same steps, bit for bit.
    subroutine check_rk4_paths()
        type(rk_array_scheme_t) :: hardcoded, loop
        real(wp) :: y1(n), y2(n), t1, t2
        logical  :: same
        integer  :: i

        call rk_scheme(hardcoded, "rk4", n, dt)
        call rk_scheme(loop, rk_tableau("rk4"), n, dt)
        y1 = initial_state()
        y2 = initial_state()
        t1 = 0.0_wp
        t2 = 0.0_wp
        same = .true.
        do i = 1, n_steps
            call hardcoded%advance(rhs, t1, y1)
            call loop%advance(rhs, t2, y2)
            same = same .and. all(y1 == y2) .and. t1 == t2
        end do
        call expect(same, "rk4: the hardcoded step and the tableau loop give the same steps, bit for bit")
    end subroutine check_rk4_paths

    subroutine check_failures()
        type(rk_array_scheme_t) :: scheme
        type(status_t) :: status
        real(wp) :: y(n), time

        call rk_scheme(scheme, "rk5", n, dt, status)
        call expect(status%failed, "an unknown name fails")
        if (status%failed) call expect(index(status%message, 'unknown Runge-Kutta method "rk5"') > 0, &
                                       "an unknown name: its message")

        status = status_t()
        call rk_scheme(scheme, rk_tableau_t(), n, dt, status)
        call expect(status%failed, "a tableau without stages fails")

        status = status_t()
        call rk_scheme(scheme, "rk4", 0, dt, status)
        call expect(status%failed, "a state size below 1 fails")

        status = status_t()
        y = initial_state()
        time = 0.0_wp
        call scheme%advance(rhs, time, y, status=status)
        call expect(status%failed, "advancing a scheme that was not built fails")

        status = status_t()
        call rk_scheme(scheme, "rk4", n + 1, dt)
        call scheme%advance(rhs, time, y, status=status)
        call expect(status%failed, "a state of the wrong size fails")
        if (status%failed) call expect(index(status%message, "size 6") > 0, "a state of the wrong size: its message")
        call expect(time == 0.0_wp .and. all(y == initial_state()), "a failed step leaves time and state unchanged")

        ! a strided state: refused, rather than copied in and out at every step
        block
            real(wp) :: wide(2*n)
            status = status_t()
            call rk_scheme(scheme, "rk4", n, dt)
            wide = 0.0_wp
            wide(1:2*n:2) = initial_state()
            call scheme%advance(rhs, time, wide(1:2*n:2), status=status)
            call expect(status%failed, "a strided state fails")
            if (status%failed) call expect(index(status%message, "contiguous") > 0, "a strided state: its message")
            call expect(all(wide(1:2*n:2) == initial_state()), "a strided state is left unchanged")
        end block
    end subroutine check_failures

end program test_RungeKutta_arrays
