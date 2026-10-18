program test_ODEmodel_boundaries
    !! The `apply_boundaries` hook (halos and boundary conditions), on every engine: each
    !! method by name (the hardcoded RK4, the generic engine for the others), and
    !! the tableau engine for every method. For each:
    !!   - every state a scheme evaluates has had its boundaries applied, at its stage time
    !!     (the test model raises an error otherwise), starting from an initial
    !!     state whose halos are wrong;
    !!   - `apply_boundaries` is called once per step plus once per stage state built;
    !!   - the steps match those of a model that fills halos on a scratch copy
    !!     instead, to round-off;
    !!   - a failing `apply_boundaries` returns its error in `status`.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: status_t, time_scheme_c, rk_scheme, rk_tableau, &
                                   & rk_method_names
    use test_ODEmodel_boundaries_model, only: halo_model_t, halo_state_t, N, MODE_HOOK, MODE_COPY, MODE_FAIL
    implicit none

    integer,  parameter :: n_steps = 5
    real(wp), parameter :: dt = 0.05_wp

    type(halo_model_t), target :: model
    class(time_scheme_c), allocatable, target :: scheme
    character(len=8), allocatable :: names(:)
    logical :: all_passed
    integer :: m

    all_passed = .true.
    names = rk_method_names()

    do m = 1, size(names)
        call check_engine(trim(names(m)), "by name", expected_calls(trim(names(m))))
        call check_engine(trim(names(m)), "by tableau", expected_calls(trim(names(m))))
    end do

    if (.not. all_passed) error stop "boundary hook checks FAILED"
    print *, "All boundary hook checks PASSED"

contains

    subroutine expect(condition, what)
        logical,      intent(in) :: condition
        character(*), intent(in) :: what

        if (.not. condition) then
            print '(2a)', "FAILED: ", what
            all_passed = .false.
        end if
    end subroutine expect

    !> Calls to `apply_boundaries` per step: one on the initial state, plus one per
    !> stage state built (stages after the first, but dp54's unused seventh).
    integer function expected_calls(name) result(n)
        character(*), intent(in) :: name

        select case (name)
        case ("rk3ws");            n = 3
        case ("rk4");              n = 4
        case ("dp54");             n = 6
        case default;              n = -1
        end select
    end function expected_calls

    !> Build the scheme for `name` with `engine` on the current model.
    subroutine build(name, engine)
        character(*), intent(in) :: name, engine

        select case (engine)
        case ("by name")
            call rk_scheme(scheme, name, dt, model)
        case ("by tableau")
            call rk_scheme(scheme, rk_tableau(name), dt, model)
        end select
    end subroutine build

    !> Run n_steps steps of the scheme built by `build`, from a smooth initial
    !> state whose halos are wrong; return the interior and the status.
    subroutine run(name, engine, y, status)
        character(*),   intent(in)    :: name, engine
        real(wp),       intent(out)   :: y(N)
        type(status_t), intent(inout) :: status

        type(halo_state_t), target :: state
        real(wp) :: time
        integer  :: i

        call build(name, engine)
        allocate(state%y(0:N+1))
        state%y = [(sin(0.7_wp * i) + 2.0_wp, i = 0, N+1)]
        state%y(0)   = -777.0_wp
        state%y(N+1) = -777.0_wp
        time = 0.0_wp
        do i = 1, n_steps
            call scheme%advance(time, state, status)
            if (status%failed) exit
        end do
        y = state%y(1:N)
    end subroutine run

    subroutine check_engine(name, engine, calls_per_step)
        character(*), intent(in) :: name, engine
        integer,      intent(in) :: calls_per_step

        type(status_t) :: status
        real(wp) :: y_hook(N), y_copy(N)
        character(len=:), allocatable :: label

        label = name // " (" // engine // ")"

        model%mode = MODE_HOOK
        model%n_boundaries = 0
        call run(name, engine, y_hook, status)
        call expect(.not. status%failed, label // ": boundaries applied to every evaluated state, at its time")
        if (status%failed) print '(4x,a)', status%message
        call expect(model%n_boundaries == n_steps * calls_per_step, label // ": calls to apply_boundaries per step")

        model%mode = MODE_COPY
        status = status_t()
        call run(name, engine, y_copy, status)
        call expect(.not. status%failed, label // ": scratch-copy model runs")
        call expect(maxval(abs(y_hook - y_copy)) <= 1.0e-14_wp * maxval(abs(y_copy)), &
                    label // ": same steps as filling halos on a scratch copy")

        model%mode = MODE_FAIL
        status = status_t()
        call run(name, engine, y_copy, status)
        call expect(status%failed, label // ": a failing apply_boundaries fails the step")
        if (status%failed) call expect(index(status%message, "apply_boundaries failed") > 0, &
                                       label // ": its message is returned")
    end subroutine check_engine

end program test_ODEmodel_boundaries
