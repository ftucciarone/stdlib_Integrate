module stdlib_integrate_rk_steps
    !! Explicit Runge-Kutta schemes on model classes, built by `rk_scheme`:
    !!
    !!   - by name, `rk_scheme(scheme, "rk4", dt, model, status)`: classic RK4
    !!     runs its hardcoded step (`time_scheme_rk4_t`), unrolled, with the
    !!     coefficients as constants; every other method of the library runs on
    !!     the generic engine, from its tableau;
    !!   - by tableau, `rk_scheme(scheme, tableau, dt, model, status)`: the
    !!     generic engine (`time_scheme_ERK_t`), for any explicit tableau,
    !!     including one built by hand.
    !!
    !! The hardcoded RK4 and the generic engine on RK4's tableau give the same
    !! steps, bit for bit (test_RungeKutta_methods).

    use stdlib_kinds, only: wp => dp
    use stdlib_integrate_status, only: status_t, fail_step
    use stdlib_integrate_stage, only: stage_t
    use stdlib_integrate_model, only: ODEmodel_c, ODEmodel_state_c, scratch_space_c, tendency_box_t
    use stdlib_integrate_scheme, only: time_scheme_c, allocate_buffers
    use stdlib_integrate_rk_tableaus, only: rk_tableau_t, rk_tableau, rk_method_names
    use stdlib_integrate_rk_generic, only: time_scheme_ERK_t
    implicit none
    private

    public :: rk_scheme
    public :: time_scheme_rk4_t

    !> Build a scheme in place, in `scheme` (`class(time_scheme_c),
    !> allocatable`): by method name, or by Butcher tableau. A failure (an
    !> unknown name, an invalid tableau, a model allocator that allocates
    !> nothing) leaves `scheme` unallocated, and is returned in `status`, or
    !> stops the program without it.
    interface rk_scheme
        module procedure rk_scheme_by_name
        module procedure rk_scheme_by_tableau
    end interface

    type, extends(time_scheme_c) :: time_scheme_rk4_t
        !! Classic RK4, hardcoded: four registers and one stage state.
        class(ODEmodel_c), pointer           :: model => null()  !! *reference* to model (not owned by us)
        class(scratch_space_c), allocatable  :: scratch
        class(ODEmodel_state_c), allocatable :: stage_state
        type(tendency_box_t), allocatable    :: k(:)
    contains
        procedure :: advance => rk4_advance
    end type

contains

    !> `scheme` = the method called `name` (see rk_method_names), with time
    !> step `dt`, for `model`: the hardcoded step for "rk4", the generic
    !> engine on the method's tableau otherwise.
    subroutine rk_scheme_by_name(scheme, name, dt, model, status)
        class(time_scheme_c), allocatable, intent(out)             :: scheme
        character(*),                      intent(in)              :: name   !! method name, e.g. "rk4"
        real(wp),                          intent(in)              :: dt     !! time step
        class(ODEmodel_c),                 intent(in),    target   :: model  !! model, kept by reference
        type(status_t),                    intent(inout), optional :: status !! failure, if present

        type(status_t)     :: local_status
        type(rk_tableau_t) :: tab
        logical            :: found

        tab = rk_tableau(name, found)
        if (.not. found) then
            call local_status%raise(unknown_method(name))
        else if (name == "rk4") then
            allocate(time_scheme_rk4_t :: scheme)
            select type (scheme)
            type is (time_scheme_rk4_t)
                scheme%description = tab%method_descriptor_t
                scheme%dt = dt
                scheme%model => model
                call allocate_buffers(model, local_status, scheme%scratch, 4, scheme%k, &
                                      stage_state=scheme%stage_state)
            end select
        else
            allocate(time_scheme_ERK_t :: scheme)
            select type (scheme)
            type is (time_scheme_ERK_t)
                call scheme%init(tab, dt, model, local_status)
            end select
        end if
        if (local_status%failed) then
            if (allocated(scheme)) deallocate(scheme)
            call fail_step(local_status, status, "rk_scheme", "cannot build the scheme")
        end if
    end subroutine rk_scheme_by_name

    !> `scheme` = the generic tableau engine for `tableau` (any explicit
    !> tableau, including one built by hand), with time step `dt`, for `model`.
    subroutine rk_scheme_by_tableau(scheme, tableau, dt, model, status)
        class(time_scheme_c), allocatable, intent(out)             :: scheme
        type(rk_tableau_t),                intent(in)              :: tableau !! Butcher tableau
        real(wp),                          intent(in)              :: dt      !! time step
        class(ODEmodel_c),                 intent(in),    target   :: model   !! model, kept by reference
        type(status_t),                    intent(inout), optional :: status  !! failure, if present

        type(status_t) :: local_status

        allocate(time_scheme_ERK_t :: scheme)
        select type (scheme)
        type is (time_scheme_ERK_t)
            call scheme%init(tableau, dt, model, local_status)
        end select
        if (local_status%failed) then
            deallocate(scheme)
            call fail_step(local_status, status, "rk_scheme", "cannot build the scheme")
        end if
    end subroutine rk_scheme_by_tableau

    !> The message for a method name that is not in the library.
    function unknown_method(name) result(message)
        character(*), intent(in)      :: name
        character(len=:), allocatable :: message

        character(len=:), allocatable :: names(:)
        integer :: i

        message = 'unknown Runge-Kutta method "' // name // '"; known:'
        names = rk_method_names()
        do i = 1, size(names)
            message = message // " " // trim(names(i))
        end do
    end function unknown_method

    subroutine rk4_advance(scheme, time, state, status)
        class(time_scheme_rk4_t), intent(inout), target   :: scheme
        real(wp),                 intent(inout)           :: time   !! time at beginning (in) or end (out) of the step
        class(ODEmodel_state_c),  intent(inout), target   :: state  !! state at beginning (in) or end (out) of the step
        type(status_t),           intent(inout), optional :: status !! failure; if absent, a failure stops

        type(status_t) :: local_status

        call step_rk4(scheme, time, state, local_status)
        if (local_status%failed) then
            call fail_step(local_status, status, "time_scheme_rk4%advance")
            return
        end if
        time = time + scheme%dt
    end subroutine rk4_advance

    !> Classic RK4 (order 4): Kutta (1901)
    subroutine step_rk4(scheme, time, state, status)
        class(time_scheme_rk4_t), intent(inout), target :: scheme
        real(wp),                 intent(in)            :: time
        class(ODEmodel_state_c),  intent(inout), target :: state
        type(status_t),           intent(inout)         :: status

        type(stage_t) :: st

        associate (model => scheme%model, scratch => scheme%scratch, dt => scheme%dt, k => scheme%k)
            ! The initial state: halos and boundary conditions, read by stage 1.
            st = stage_t(index=1, stages=4, c=0.0_wp, dt=dt, first=.true., last=.false.)
            call model%apply_boundaries(status, state, time, st)
            if (status%failed) return
            ! stage 1: k1 = f(y)  ->  register 1
            call model%compute_tendencies(status, state, time, scratch, k(1)%item, st)
            if (status%failed) return
            ! stage 2: k2 = f(y + dt*(1/2*k1))  ->  register 2
            st = stage_t(index=2, stages=4, c=1.0_wp/2.0_wp, dt=dt, first=.false., last=.false.)
            call scheme%stage_state%combine(status, state, [dt*(1.0_wp/2.0_wp)], &
                                            k(1:1))
            if (status%failed) return
            call model%apply_boundaries(status, scheme%stage_state, time + (1.0_wp/2.0_wp)*dt, st)
            if (status%failed) return
            call model%compute_tendencies(status, scheme%stage_state, time + (1.0_wp/2.0_wp)*dt, scratch, &
                                          k(2)%item, st)
            if (status%failed) return
            ! stage 3: k3 = f(y + dt*(1/2*k2))  ->  register 3
            st = stage_t(index=3, stages=4, c=1.0_wp/2.0_wp, dt=dt, first=.false., last=.false.)
            call scheme%stage_state%combine(status, state, [dt*(1.0_wp/2.0_wp)], &
                                            k(2:2))
            if (status%failed) return
            call model%apply_boundaries(status, scheme%stage_state, time + (1.0_wp/2.0_wp)*dt, st)
            if (status%failed) return
            call model%compute_tendencies(status, scheme%stage_state, time + (1.0_wp/2.0_wp)*dt, scratch, &
                                          k(3)%item, st)
            if (status%failed) return
            ! stage 4: k4 = f(y + dt*(1*k3))  ->  register 4
            st = stage_t(index=4, stages=4, c=1.0_wp, dt=dt, first=.false., last=.true.)
            call scheme%stage_state%combine(status, state, [dt*(1.0_wp)], &
                                            k(3:3))
            if (status%failed) return
            call model%apply_boundaries(status, scheme%stage_state, time + dt, st)
            if (status%failed) return
            call model%compute_tendencies(status, scheme%stage_state, time + dt, scratch, &
                                          k(4)%item, st)
            if (status%failed) return
            ! y = y + dt*(1/6*k1 + 1/3*k2 + 1/3*k3 + 1/6*k4)
            call state%accumulate(status, [dt*(1.0_wp/6.0_wp), dt*(1.0_wp/3.0_wp), dt*(1.0_wp/3.0_wp), &
                                           dt*(1.0_wp/6.0_wp)], k(1:4))
        end associate
    end subroutine step_rk4

end module stdlib_integrate_rk_steps
