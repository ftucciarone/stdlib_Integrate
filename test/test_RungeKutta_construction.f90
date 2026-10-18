program test_RungeKutta_construction
    !! Building a scheme in place (rk_scheme by name and by tableau, and the
    !! type-bound init of the tableau engine):
    !!   - a valid construction leaves status clear and the scheme usable;
    !!   - every failure (unknown name, invalid tableau, a model allocator
    !!     that allocates nothing) is returned in status with a message, and
    !!     rk_scheme leaves the scheme unallocated;
    !!   - by name, rk4 is the hardcoded step and every other method the
    !!     tableau engine; both need a stage state, so a model that cannot
    !!     allocate one cannot build them;
    !!   - init on an existing scheme rebuilds it (no stale buffers);
    !!   - the scheme refers to the model: a change to the model after
    !!     construction is seen by the next step.
    !! (A failure without status reports and stops, as advance does: not
    !! testable inside this program.)

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: status_t, time_scheme_c, rk_scheme, rk_tableau, rk_tableau_t
    use stdlib_integrate_rk_steps, only: time_scheme_rk4_t
    use stdlib_integrate_rk_generic, only: time_scheme_ERK_t
    use test_decay_model, only: decay_model_t, decay_state_t
    use test_RungeKutta_construction_models, only: faulty_model_t
    implicit none

    type(decay_model_t), target  :: model
    type(faulty_model_t), target :: faulty
    class(time_scheme_c), allocatable :: scheme
    type(time_scheme_ERK_t)   :: erk
    type(rk_tableau_t) :: tab
    type(status_t) :: status
    logical :: all_passed

    all_passed = .true.
    model%lambda = 2.0_wp

    ! ---- valid constructions ----
    status = status_t()
    call rk_scheme(scheme, "rk4", 0.1_wp, model, status)
    call expect(.not. status%failed .and. allocated(scheme), "by name: built, status clear")
    call expect(same_type_as(scheme, time_scheme_rk4_t()), "by name: rk4 is the hardcoded step")
    call expect(abs(one_step(scheme) - exp(-0.2_wp)) < 1.0e-5_wp, "by name: steps correctly")

    status = status_t()
    call rk_scheme(scheme, "dp54", 0.1_wp, model, status)
    call expect(.not. status%failed .and. allocated(scheme), "by name, dp54: built, status clear")
    if (allocated(scheme)) call expect(same_type_as(scheme, erk), "by name, dp54: the generic tableau engine")

    status = status_t()
    call rk_scheme(scheme, rk_tableau("rk4"), 0.1_wp, model, status)
    call expect(.not. status%failed .and. allocated(scheme), "by tableau: built, status clear")
    call expect(same_type_as(scheme, erk), "by tableau: the generic tableau engine")
    call expect(abs(one_step(scheme) - exp(-0.2_wp)) < 1.0e-5_wp, "by tableau: steps correctly")

    ! ---- unknown name ----
    status = status_t()
    call rk_scheme(scheme, "nope", 0.1_wp, model, status)
    call expect_failure("unknown name", 'unknown Runge-Kutta method "nope"')
    call expect(.not. allocated(scheme), "unknown name: scheme left unallocated")
    if (allocated(status%message)) &
        call expect(index(status%message, "rk4") > 0, "unknown name: message lists the known methods")

    ! ---- invalid tableaus ----
    status = status_t()
    call rk_scheme(scheme, rk_tableau("nope"), 0.1_wp, model, status)
    call expect_failure("empty tableau", "no stages")
    call expect(.not. allocated(scheme), "empty tableau: scheme left unallocated")

    tab = rk_tableau("rk3ws")
    tab%A(1, 1) = 1.0_wp
    status = status_t()
    call rk_scheme(scheme, tab, 0.1_wp, model, status)
    call expect_failure("implicit tableau", "strictly lower triangular")

    tab = rk_tableau("rk3ws")
    deallocate(tab%c)
    allocate(tab%c(4), source=0.0_wp)   ! rk3ws has 3 stages
    status = status_t()
    call rk_scheme(scheme, tab, 0.1_wp, model, status)
    call expect_failure("mis-sized tableau", "do not match the number of stages")

    status = status_t()
    call erk%init(rk_tableau("nope"), 0.1_wp, model, status)
    call expect_failure("ERK init, empty tableau", "no stages")

    ! ---- model allocators that allocate nothing ----
    faulty%no_tendencies = .true.
    status = status_t()
    call rk_scheme(scheme, "rk4", 0.1_wp, faulty, status)
    call expect_failure("no tendencies, by name", "allocate_tendencies")
    call expect(.not. allocated(scheme), "no tendencies, by name: scheme left unallocated")
    status = status_t()
    call rk_scheme(scheme, rk_tableau("rk4"), 0.1_wp, faulty, status)
    call expect_failure("no tendencies, by tableau", "allocate_tendencies")
    status = status_t()
    call rk_scheme(scheme, "rk3ws", 0.1_wp, faulty, status)
    call expect_failure("no tendencies, the engine by name", "allocate_tendencies")

    faulty%no_tendencies = .false.
    faulty%no_scratch = .true.
    status = status_t()
    call rk_scheme(scheme, "rk4", 0.1_wp, faulty, status)
    call expect_failure("no scratch", "allocate_scratch")

    faulty%no_scratch = .false.
    faulty%no_state = .true.
    status = status_t()
    call rk_scheme(scheme, "rk4", 0.1_wp, faulty, status)
    call expect_failure("no state, rk4 (needs a stage state)", "allocate_state")

    ! ---- init rebuilds an existing scheme ----
    call erk%init(rk_tableau("rk4"), 0.1_wp, model)
    call expect(size(erk%k) == 4, "ERK init: four registers for rk4")
    call erk%init(rk_tableau("rk3ws"), 0.1_wp, model)
    call expect(size(erk%k) == 3 .and. erk%tableau%stages == 3, "ERK init again: rebuilt for rk3ws")
    ! rk3ws is of order 3 on this linear problem: R(z) = 1 + z + z**2/2 + z**3/6, z = -0.2
    call expect(abs(one_step(erk) - (1.0_wp - 0.2_wp + 0.02_wp - 0.008_wp/6)) < 1.0e-14_wp, &
                "ERK init again: steps as rk3ws")

    ! ---- the scheme refers to the model ----
    call rk_scheme(scheme, "rk4", 0.1_wp, model)
    model%lambda = 1.0_wp
    ! rk4 on y' = -y, dt = 0.1: R(z) = 1 + z + z**2/2 + z**3/6 + z**4/24, z = -0.1
    call expect(abs(one_step(scheme) - (1.0_wp - 0.1_wp + 0.005_wp - 0.001_wp/6 + 0.0001_wp/24)) < 1.0e-14_wp, &
                "a change to the model is seen by the next step")
    model%lambda = 2.0_wp

    if (.not. all_passed) error stop "construction checks FAILED"
    print *, "All construction checks PASSED"

contains

    subroutine expect(condition, what)
        logical, intent(in) :: condition
        character(*), intent(in) :: what

        if (.not. condition) then
            print '(2a)', "FAILED: ", what
            all_passed = .false.
        end if
    end subroutine expect

    !> The last construction failed, with a message containing `fragment`.
    subroutine expect_failure(what, fragment)
        character(*), intent(in) :: what, fragment

        call expect(status%failed, what // ": fails")
        if (allocated(status%message)) then
            call expect(index(status%message, fragment) > 0, what // ": message mentions '" // fragment // "'")
        else
            call expect(.false., what // ": failure has a message")
        end if
    end subroutine expect_failure

    !> y after one step from y = 1.
    real(wp) function one_step(scheme) result(y)
        class(time_scheme_c), intent(inout) :: scheme
        type(decay_state_t) :: state
        type(status_t) :: step_status
        real(wp) :: time

        state%y = 1.0_wp
        time = 0.0_wp
        call scheme%advance(time, state, step_status)
        if (step_status%failed) error stop "step failed"
        y = state%y
    end function one_step

end program test_RungeKutta_construction
