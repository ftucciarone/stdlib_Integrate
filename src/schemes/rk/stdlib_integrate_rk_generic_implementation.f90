submodule(stdlib_integrate_rk_generic) stdlib_integrate_rk_generic_implementation
    use stdlib_integrate_status, only: fail_step
    use stdlib_integrate_stage, only: stage_t
    use stdlib_integrate_scheme, only: allocate_buffers
    use stdlib_integrate_rk_tableaus, only: check_tableau, stage_unused
    !
contains
    !
    module subroutine ERK_init(scheme, tableau, dt, model, status)
        class(time_scheme_ERK_t), intent(out)             :: scheme
        type(rk_tableau_t),       intent(in)              :: tableau
        real(wp),                 intent(in)              :: dt
        class(ODEmodel_c),        intent(in),    target   :: model
        type(status_t),           intent(inout), optional :: status
        !
        type(status_t)   :: local_status
        integer          :: i

        call check_tableau(tableau, local_status)
        if (.not. local_status%failed) then
            scheme%dt = dt
            scheme%model => model
            scheme%tableau = tableau
            scheme%description = tableau%method_descriptor_t
            ! What every step needs from the tableau, once
            allocate(scheme%used(tableau%stages), scheme%combines(tableau%stages), scheme%coef(tableau%stages))
            do i = 1, tableau%stages
                scheme%used(i) = .not. stage_unused(tableau, i)
                scheme%combines(i) = any(tableau%A(i, 1:i-1) /= 0.0_wp)
                if (scheme%used(i)) scheme%last = i
            end do
            if (scheme%last == 0) scheme%last = 1     ! a tableau of unused stages still evaluates stage 1
            ! The engine's own cost (also for a tableau built by hand): the stages it evaluates
            scheme%description%f_evals = count(scheme%used)
            call allocate_buffers(model, local_status, scheme%scratch, tableau%stages, scheme%k, &
                                  stage_state=scheme%stage_state)
        end if
        if (local_status%failed) &
            call fail_step(local_status, status, "time_scheme_ERK%init", "cannot build the scheme")
    end subroutine ERK_init

    module subroutine ERK_advance(scheme, time, state, status)
        class(time_scheme_ERK_t), intent(inout), target   :: scheme
        real(wp),                 intent(inout)           :: time
        class(ODEmodel_state_c),  intent(inout), target   :: state
        type(status_t),           intent(inout), optional :: status
        !
        type(status_t) :: local_status
        type(stage_t)  :: stage
        integer        :: i

        associate (      model => scheme%model,      &
                  &     scratch => scheme%scratch,    &
                  &          dt => scheme%dt,         &
                  &         tab => scheme%tableau,    &
                  &           k => scheme%k,          &
                  & stage_state => scheme%stage_state, &
                  &        coef => scheme%coef,       &
                  &        last => scheme%last        )

            ! ---- the initial state: halos and boundary conditions, read by stage 1 ----
            stage = stage_t(index=1, stages=tab%stages, c=tab%c(1), dt=dt, first=.true., last=(last == 1))
            call model%apply_boundaries(local_status, state, time, stage)
            if (fail()) return

            ! ---- stages: k(i) = f( state + dt * sum_j<i A(i,j) k(j),  time + c(i)*dt ) ----
            do i = 1, tab%stages
                ! A stage nothing uses (only the embedded error estimate would): skip it.
                if (.not. scheme%used(i)) cycle
                stage = stage_t(index=i, stages=tab%stages, c=tab%c(i), dt=dt, first=(i == 1), last=(i == last))

                if (scheme%combines(i)) then
                    ! The coefficients go into the scheme's buffer: no temporary per stage.
                    coef(1:i-1) = dt*tab%A(i, 1:i-1)
                    call stage_state%combine(local_status, state, coef(1:i-1), k(1:i-1))
                    if (fail()) return
                    call model%apply_boundaries(local_status, stage_state, time + tab%c(i)*dt, stage)
                    if (fail()) return
                    call model%compute_tendencies(local_status, stage_state, time + tab%c(i)*dt, scratch, k(i)%item, &
                                                  stage)
                else
                    ! Stage state = state (always for stage 1): no copy.
                    call model%compute_tendencies(local_status, state, time + tab%c(i)*dt, scratch, k(i)%item, stage)
                end if
                if (fail()) return
            end do

            ! ---- state <- state + dt * sum_i b(i) k(i), in one pass ----
            coef(1:tab%stages) = dt*tab%b
            call state%accumulate(local_status, coef(1:tab%stages), k)
            if (fail()) return

            time = time + dt
        end associate

    contains

        ! On failure: return it in `status`, or report it and stop.
        logical function fail()
            fail = local_status%failed
            if (fail) call fail_step(local_status, status, "time_scheme_ERK%advance")
        end function fail

    end subroutine ERK_advance

end submodule stdlib_integrate_rk_generic_implementation
