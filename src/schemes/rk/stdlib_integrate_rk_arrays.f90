module stdlib_integrate_rk_arrays
    !! Explicit Runge-Kutta schemes on plain rank-1 arrays: no model, state or
    !! tendency classes, just a right-hand side `f(t, y, dydt, stage)`. Built
    !! by `rk_scheme`, as on the model classes (stdlib_integrate_rk_steps):
    !!
    !!   - by name, `rk_scheme(scheme, "rk4", n, dt, status)`: classic RK4 runs
    !!     its hardcoded step; every other method of the library runs a loop
    !!     over its tableau;
    !!   - by tableau, `rk_scheme(scheme, tableau, n, dt, status)`: the loop
    !!     over any explicit tableau, including one built by hand.
    !!
    !! Both run the operations of their class-path counterparts in the same
    !! order (the same stages skipped, every combination adding the same terms
    !! in the same order), so the two routes give the same steps, bit for bit.
    !!
    !! The right-hand side is passed to every `advance`, never stored: the
    !! scheme holds no reference to anything. It may be an internal procedure
    !! of the caller, which then reads the caller's variables (parameters,
    !! grids) by host association.

    use stdlib_kinds, only: wp => dp
    use stdlib_integrate_status, only: status_t, fail_step
    use stdlib_integrate_stage, only: stage_t
    use stdlib_integrate_descriptor, only: method_descriptor_t
    use stdlib_integrate_rk_tableaus, only: rk_tableau_t, rk_tableau, rk_method_names, check_tableau, &
                                            stage_unused
    implicit none
    private

    public :: rk_array_scheme_t
    public :: rk_rhs_i, rk_boundaries_i
    public :: rk_scheme

    !> Build an array scheme in place, as on the model classes
    !> (stdlib_integrate_rk_steps, which extends the same generic): the state
    !> size `n` in place of the model.
    interface rk_scheme
        module procedure rka_scheme_by_name
        module procedure rka_scheme_by_tableau
    end interface

    abstract interface
        !> Right-hand side: `dydt = f(t, y)`. `stage` tells which stage of the
        !> step this evaluation belongs to (see stdlib_integrate_stage).
        subroutine rk_rhs_i(t, y, dydt, stage)
            import wp, stage_t
            real(wp),      intent(in)  :: t       !! time
            real(wp),      intent(in)  :: y(:)    !! state
            real(wp),      intent(out) :: dydt(:) !! tendencies, size(y)
            type(stage_t), intent(in)  :: stage   !! stage this evaluation belongs to
        end subroutine rk_rhs_i

        !> Optional: fill the halos and apply the boundary conditions of `y`,
        !> in place, before it is evaluated (see ODEmodel_c's `apply_boundaries`).
        subroutine rk_boundaries_i(t, y, stage)
            import wp, stage_t
            real(wp),      intent(in)    :: t     !! time of `y`
            real(wp),      intent(inout) :: y(:)  !! state, its boundaries applied in place
            type(stage_t), intent(in)    :: stage !! stage they are applied for
        end subroutine rk_boundaries_i
    end interface

    type :: rk_array_scheme_t
        !! An explicit Runge-Kutta method on arrays of a fixed size, with its
        !! buffers. Build it with `rk_scheme(scheme, name | tableau, n, dt, status)`.
        logical                   :: built = .false.
        logical                   :: hardcoded_rk4 = .false. !! the hardcoded RK4 step, or the tableau loop
        integer                   :: n = 0                   !! size of the state
        real(wp)                  :: dt = 0.0_wp
        !> Tendency registers, k(:, slot): 4 for RK4, one per stage for the tableau loop
        real(wp), allocatable     :: k(:,:)
        !> Stage state
        real(wp), allocatable     :: ys(:)
        !> The tableau loop: the tableau, and what a step needs from it, computed when built
        type(rk_tableau_t)        :: tableau
        logical,  allocatable     :: used(:)      !! stage i's tendency is read by b or a later row of A
        logical,  allocatable     :: combines(:)  !! stage i builds a stage state (a nonzero A(i,1:i-1))
        integer                   :: last = 0     !! the last stage evaluated
        real(wp), allocatable     :: coef(:)      !! dt*A(i,:) or dt*b
        !> The scheme's method; see `descriptor`
        type(method_descriptor_t) :: description
    contains
        procedure :: advance => rka_advance
        !> The descriptor of the scheme's method, as time_scheme_c's `descriptor`
        procedure :: descriptor => rka_descriptor
        !> The time step, as time_scheme_c's `get_dt`
        procedure :: get_dt  => rka_get_dt
        !> Change the time step, as time_scheme_c's `set_dt`
        procedure :: set_dt  => rka_set_dt
    end type

contains

    !> Build the scheme of the method called `name` (see rk_method_names) for
    !> states of size `n`, with time step `dt` (any previous content is
    !> discarded): the hardcoded step for "rk4", the tableau loop otherwise. An
    !> unknown name or `n < 1` is a failure: returned in `status`, or stopping
    !> the program without it.
    subroutine rka_scheme_by_name(scheme, name, n, dt, status)
        type(rk_array_scheme_t), intent(out)             :: scheme
        character(*),            intent(in)              :: name   !! method name, e.g. "rk4"
        integer,                 intent(in)              :: n      !! size of the state
        real(wp),                intent(in)              :: dt     !! time step
        type(status_t),          intent(inout), optional :: status !! failure, if present

        type(status_t)     :: local_status
        type(rk_tableau_t) :: tab
        logical            :: found
        character(len=:), allocatable :: names(:)
        character(len=:), allocatable :: message
        integer :: i

        tab = rk_tableau(name, found)
        if (.not. found) then
            message = 'unknown Runge-Kutta method "' // name // '"; known:'
            names = rk_method_names()
            do i = 1, size(names)
                message = message // " " // trim(names(i))
            end do
            call local_status%raise(message)
        else if (n < 1) then
            call local_status%raise("the state size must be at least 1")
        end if
        if (local_status%failed) then
            call fail_step(local_status, status, "rk_scheme", "cannot build the scheme")
            return
        end if
        if (name == "rk4") then
            scheme%hardcoded_rk4 = .true.
            scheme%description = tab%method_descriptor_t
            allocate(scheme%k(n, 4), scheme%ys(n))
            scheme%n = n
            scheme%dt = dt
            scheme%built = .true.
        else
            call build_loop(scheme, tab, n, dt)
        end if
    end subroutine rka_scheme_by_name

    !> Build the tableau loop for `tableau` (any explicit tableau, including
    !> one built by hand), for states of size `n`, with time step `dt`. An
    !> invalid tableau or `n < 1` is a failure.
    subroutine rka_scheme_by_tableau(scheme, tableau, n, dt, status)
        type(rk_array_scheme_t), intent(out)             :: scheme
        type(rk_tableau_t),      intent(in)              :: tableau !! Butcher tableau
        integer,                 intent(in)              :: n       !! size of the state
        real(wp),                intent(in)              :: dt      !! time step
        type(status_t),          intent(inout), optional :: status  !! failure, if present

        type(status_t) :: local_status

        call check_tableau(tableau, local_status)
        if (.not. local_status%failed .and. n < 1) call local_status%raise("the state size must be at least 1")
        if (local_status%failed) then
            call fail_step(local_status, status, "rk_scheme", "cannot build the scheme")
            return
        end if
        call build_loop(scheme, tableau, n, dt)
    end subroutine rka_scheme_by_tableau

    !> The tableau loop, for a checked tableau: what every step needs from the
    !> tableau, once, as the class-path engine (time_scheme_ERK_t) computes it.
    subroutine build_loop(scheme, tableau, n, dt)
        type(rk_array_scheme_t), intent(inout) :: scheme
        type(rk_tableau_t),      intent(in)    :: tableau
        integer,                 intent(in)    :: n
        real(wp),                intent(in)    :: dt

        integer :: i, s

        s = tableau%stages
        scheme%tableau = tableau
        scheme%description = tableau%method_descriptor_t
        allocate(scheme%used(s), scheme%combines(s), scheme%coef(s))
        do i = 1, s
            scheme%used(i) = .not. stage_unused(tableau, i)
            scheme%combines(i) = any(tableau%A(i, 1:i-1) /= 0.0_wp)
            if (scheme%used(i)) scheme%last = i
        end do
        if (scheme%last == 0) scheme%last = 1     ! a tableau of unused stages still evaluates stage 1
        scheme%description%f_evals = count(scheme%used)
        allocate(scheme%k(n, s), scheme%ys(n))
        scheme%n = n
        scheme%dt = dt
        scheme%built = .true.
    end subroutine build_loop

    !> Advance `y` and `time` by one step of the scheme's method, evaluating the
    !> right-hand side `f`. If `boundaries` is present, it is called on `y` at the
    !> start of the step and on each stage state right after it is built. A
    !> scheme that was not built, or a `y` of the wrong size or not contiguous,
    !> is a failure: returned in `status`, or stopping the program without it.
    subroutine rka_advance(scheme, f, time, y, boundaries, status)
        class(rk_array_scheme_t), intent(inout)           :: scheme
        procedure(rk_rhs_i)                               :: f        !! right-hand side
        real(wp),                 intent(inout)           :: time     !! time at beginning (in) or end (out) of the step
        real(wp),                 intent(inout), target   :: y(:)     !! state at beginning (in) or end (out) of the step;
                                                                      !! contiguous (checked: a strided one would be copied)
        procedure(rk_boundaries_i),              optional :: boundaries !! halos and boundary conditions
        type(status_t),           intent(inout), optional :: status   !! failure; if absent, a failure stops

        type(status_t) :: local_status
        ! y as a contiguous pointer, once it is checked contiguous: the steps
        ! declare their y contiguous, and handing them an assumed-shape y costs
        ! the compiler's argument handling at every step (measured: +20 ns on a
        ! 23 ns step at n = 3, +23% at n = 1000); a contiguous pointer costs nothing.
        real(wp), pointer, contiguous :: yc(:)

        if (.not. scheme%built) then
            call local_status%raise("the scheme was not built (call rk_scheme first)")
        else if (size(y) /= scheme%n) then
            call local_status%raise("the state has size " // itoa(size(y)) // ", the scheme was built for " // &
                                    itoa(scheme%n))
        else if (.not. is_contiguous(y)) then
            ! The steps work on y in place, as a contiguous array: a strided one
            ! would be copied in and out at every step, behind the caller's back.
            call local_status%raise("the state must be a contiguous array, not a strided section")
        end if
        if (local_status%failed) then
            call fail_step(local_status, status, "rk_array_scheme%advance")
            return
        end if
        yc => y

        if (scheme%hardcoded_rk4) then
            call step_rk4(f, boundaries, time, scheme%dt, yc, scheme%k, scheme%ys)
        else
            call step_loop(scheme, f, boundaries, time, yc)
        end if
        time = time + scheme%dt
    end subroutine rka_advance

    !> Classic RK4 (order 4): Kutta (1901)
    subroutine step_rk4(f, boundaries, time, dt, y, k, ys)
        procedure(rk_rhs_i)                         :: f
        procedure(rk_boundaries_i),        optional :: boundaries
        real(wp),                 intent(in)        :: time, dt
        real(wp), contiguous,     intent(inout)     :: y(:)
        real(wp), contiguous,     intent(inout)     :: k(:,:), ys(:)

        type(stage_t) :: st

        st = stage_t(index=1, stages=4, c=0.0_wp, dt=dt, first=.true., last=.false.)
        if (present(boundaries)) call boundaries(time, y, st)
        ! stage 1: k1 = f(y)  ->  register 1
        call f(time, y, k(:, 1), st)
        ! stage 2: k2 = f(y + dt*(1/2*k1))  ->  register 2
        st = stage_t(index=2, stages=4, c=1.0_wp/2.0_wp, dt=dt, first=.false., last=.false.)
        ys = y + (dt*(1.0_wp/2.0_wp))*k(:, 1)
        if (present(boundaries)) call boundaries(time + (1.0_wp/2.0_wp)*dt, ys, st)
        call f(time + (1.0_wp/2.0_wp)*dt, ys, k(:, 2), st)
        ! stage 3: k3 = f(y + dt*(1/2*k2))  ->  register 3
        st = stage_t(index=3, stages=4, c=1.0_wp/2.0_wp, dt=dt, first=.false., last=.false.)
        ys = y + (dt*(1.0_wp/2.0_wp))*k(:, 2)
        if (present(boundaries)) call boundaries(time + (1.0_wp/2.0_wp)*dt, ys, st)
        call f(time + (1.0_wp/2.0_wp)*dt, ys, k(:, 3), st)
        ! stage 4: k4 = f(y + dt*(1*k3))  ->  register 4
        st = stage_t(index=4, stages=4, c=1.0_wp, dt=dt, first=.false., last=.true.)
        ys = y + (dt*(1.0_wp))*k(:, 3)
        if (present(boundaries)) call boundaries(time + dt, ys, st)
        call f(time + dt, ys, k(:, 4), st)
        ! y = y + dt*(1/6*k1 + 1/3*k2 + 1/3*k3 + 1/6*k4)
        y = y + (dt*(1.0_wp/6.0_wp))*k(:, 1) &
          & + (dt*(1.0_wp/3.0_wp))*k(:, 2) &
          & + (dt*(1.0_wp/3.0_wp))*k(:, 3) &
          & + (dt*(1.0_wp/6.0_wp))*k(:, 4)
    end subroutine step_rk4

    !> One step of the tableau loop: the operations of time_scheme_ERK_t's
    !> advance, in the same order (a stage state is `y`, then one term added
    !> per nonzero coefficient; the closing update adds one term per nonzero
    !> weight), on arrays.
    subroutine step_loop(scheme, f, boundaries, time, y)
        type(rk_array_scheme_t),    intent(inout) :: scheme
        procedure(rk_rhs_i)                       :: f
        procedure(rk_boundaries_i), optional      :: boundaries
        real(wp),                   intent(in)    :: time
        real(wp), contiguous,       intent(inout) :: y(:)

        type(stage_t) :: st
        integer       :: i, m

        associate (tab => scheme%tableau, dt => scheme%dt, k => scheme%k, ys => scheme%ys, &
                   coef => scheme%coef, last => scheme%last)
            ! The initial state: halos and boundary conditions, read by stage 1.
            st = stage_t(index=1, stages=tab%stages, c=tab%c(1), dt=dt, first=.true., last=(last == 1))
            if (present(boundaries)) call boundaries(time, y, st)

            do i = 1, tab%stages
                ! A stage nothing uses (only an embedded error estimate would): skip it.
                if (.not. scheme%used(i)) cycle
                st = stage_t(index=i, stages=tab%stages, c=tab%c(i), dt=dt, first=(i == 1), last=(i == last))
                if (scheme%combines(i)) then
                    coef(1:i-1) = dt*tab%A(i, 1:i-1)
                    ys = y
                    do m = 1, i - 1
                        if (coef(m) /= 0.0_wp) ys = ys + coef(m)*k(:, m)
                    end do
                    if (present(boundaries)) call boundaries(time + tab%c(i)*dt, ys, st)
                    call f(time + tab%c(i)*dt, ys, k(:, i), st)
                else
                    ! Stage state = y (always for stage 1): no copy.
                    call f(time + tab%c(i)*dt, y, k(:, i), st)
                end if
            end do

            ! y <- y + dt * sum_i b(i) k(i)
            coef(1:tab%stages) = dt*tab%b
            do m = 1, tab%stages
                if (coef(m) /= 0.0_wp) y = y + coef(m)*k(:, m)
            end do
        end associate
    end subroutine step_loop

    pure type(method_descriptor_t) function rka_descriptor(scheme) result(description)
        class(rk_array_scheme_t), intent(in) :: scheme
        description = scheme%description
    end function rka_descriptor

    pure real(wp) function rka_get_dt(scheme) result(dt)
        class(rk_array_scheme_t), intent(in) :: scheme
        dt = scheme%dt
    end function rka_get_dt

    subroutine rka_set_dt(scheme, dt)
        class(rk_array_scheme_t), intent(inout) :: scheme
        real(wp),                 intent(in)    :: dt
        scheme%dt = dt
    end subroutine rka_set_dt

    pure function itoa(i) result(s)
        integer, intent(in) :: i
        character(:), allocatable :: s
        character(len=16) :: buffer

        write(buffer, '(i0)') i
        s = trim(buffer)
    end function itoa

end module stdlib_integrate_rk_arrays
