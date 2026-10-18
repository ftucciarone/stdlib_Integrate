program test_ODEstage
    !! The stage information (stage_t) every engine passes to the model, for
    !! every method: rk_scheme by name (the hardcoded RK4, the generic engine
    !! for the others), the tableau engine (rk_scheme by tableau) and the array
    !! path (rk_array_scheme_t, by name and by tableau). Over two steps, each
    !! step must be,
    !! stage after stage of the method (stages nothing uses skipped):
    !! `apply_boundaries` then the evaluation, both with
    !!   - index = the stage, stages = the method's number of stages,
    !!     c = the tableau's c(index), dt = the step;
    !!   - first only on the first evaluation, last only on the last;
    !!   - time = t_n + c*dt;
    !! and the first evaluation on exactly the step's initial state y_n.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: time_scheme_c, rk_scheme, rk_method_names, &
                                   & rk_tableau, rk_tableau_t, rk_array_scheme_t
    use test_ODEstage_model, only: stage_model_t, stage_state_t, calls, n_calls, forget_calls, &
                                   & rec_rhs, rec_boundaries
    implicit none

    real(wp), parameter :: dt = 0.1_wp
    integer,  parameter :: n_steps = 2

    character(len=8), allocatable :: names(:)
    type(stage_model_t), target :: model
    type(rk_tableau_t) :: tab
    logical :: all_passed
    integer :: m

    all_passed = .true.
    names = rk_method_names()
    do m = 1, size(names)
        tab = rk_tableau(trim(names(m)))
        call check_class(trim(names(m)), by_tableau=.false.)
        call check_class(trim(names(m)), by_tableau=.true.)
        call check_arrays(trim(names(m)), by_tableau=.false.)
        call check_arrays(trim(names(m)), by_tableau=.true.)
    end do

    if (.not. all_passed) error stop "stage checks FAILED"
    print *, "All stage checks PASSED"

contains

    subroutine expect(condition, what)
        logical,      intent(in) :: condition
        character(*), intent(in) :: what

        if (.not. condition) then
            print '(2a)', "FAILED: ", what
            all_passed = .false.
        end if
    end subroutine expect

    ! Stage i is skipped when nothing uses its tendency.
    logical function unused(i)
        integer, intent(in) :: i
        unused = tab%b(i) == 0.0_wp .and. all(tab%A(i+1:, i) == 0.0_wp)
    end function unused

    ! The calls recorded during one step from (t_n, y_n), against the tableau `tab`.
    subroutine check_step(what, t_n, y_n)
        character(*), intent(in) :: what
        real(wp),     intent(in) :: t_n, y_n

        integer :: i, n, first, last
        character(len=:), allocatable :: at

        first = 0
        last = 0
        do i = tab%stages, 1, -1
            if (.not. unused(i)) first = i
        end do
        do i = 1, tab%stages
            if (.not. unused(i)) last = i
        end do

        n = 0
        do i = 1, tab%stages
            if (unused(i)) cycle
            at = what // ", stage " // achar(iachar('0') + i)
            if (n + 2 > n_calls) then
                call expect(.false., at // ": missing calls")
                return
            end if
            call expect(calls(n + 1)%kind == 'c' .and. calls(n + 2)%kind == 'f', at // ": apply_boundaries, then the evaluation")
            call check_call(n + 1, i, first, last, t_n, at // " (apply_boundaries)")
            call check_call(n + 2, i, first, last, t_n, at // " (evaluation)")
            if (i == first) call expect(calls(n + 2)%y == y_n, at // ": evaluated on y_n")
            n = n + 2
        end do
        call expect(n == n_calls, what // ": no other calls")
    end subroutine check_step

    subroutine check_call(j, i, first, last, t_n, at)
        integer,      intent(in) :: j, i, first, last
        real(wp),     intent(in) :: t_n
        character(*), intent(in) :: at

        associate (s => calls(j)%stage)
            call expect(s%index == i, at // ": index")
            call expect(s%stages == tab%stages, at // ": stages")
            call expect(s%c == tab%c(i), at // ": c")
            call expect(s%dt == dt, at // ": dt")
            call expect(s%first .eqv. (i == first), at // ": first")
            call expect(s%last .eqv. (i == last), at // ": last")
            call expect(calls(j)%time == t_n + s%c*s%dt, at // ": time = t_n + c*dt")
        end associate
    end subroutine check_call

    subroutine check_class(name, by_tableau)
        character(*), intent(in) :: name
        logical,      intent(in) :: by_tableau

        class(time_scheme_c), allocatable :: scheme
        type(stage_state_t) :: state
        real(wp) :: time, t_n, y_n
        integer  :: step
        character(len=:), allocatable :: what

        if (by_tableau) then
            call rk_scheme(scheme, tab, dt, model)
            what = "tableau engine " // name
        else
            call rk_scheme(scheme, name, dt, model)
            what = "by name " // name
        end if
        state%y = 1.0_wp
        time = 0.0_wp
        do step = 1, n_steps
            t_n = time
            y_n = state%y
            call forget_calls()
            call scheme%advance(time, state)
            call check_step(what, t_n, y_n)
        end do
    end subroutine check_class

    subroutine check_arrays(name, by_tableau)
        character(*), intent(in) :: name
        logical,      intent(in) :: by_tableau

        type(rk_array_scheme_t) :: scheme
        real(wp) :: y(1), time, t_n, y_n
        integer  :: step
        character(len=:), allocatable :: what

        if (by_tableau) then
            call rk_scheme(scheme, tab, 1, dt)
            what = "array tableau loop " // name
        else
            call rk_scheme(scheme, name, 1, dt)
            what = "array path " // name
        end if
        y = 1.0_wp
        time = 0.0_wp
        do step = 1, n_steps
            t_n = time
            y_n = y(1)
            call forget_calls()
            call scheme%advance(rec_rhs, time, y, boundaries=rec_boundaries)
            call check_step(what, t_n, y_n)
        end do
    end subroutine check_arrays

end program test_ODEstage
