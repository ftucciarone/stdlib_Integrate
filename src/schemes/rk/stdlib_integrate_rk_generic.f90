module stdlib_integrate_rk_generic
    !! Generic explicit Runge-Kutta stepper, driven by a Butcher tableau
    !! from [[stdlib_integrate_rk_tableaus]]. A single `advance` implementation
    !! (in the [[stdlib_integrate_rk_generic_implementation]] submodule) serves every
    !! explicit RK method - adding a method means adding a tableau, not a type.

    use stdlib_kinds, only: wp => dp
    use stdlib_integrate_status, only: status_t
    use stdlib_integrate_model, only: ODEmodel_c, ODEmodel_state_c, ODEmodel_tendencies_c, scratch_space_c, &
                             & tendency_box_t
    use stdlib_integrate_scheme, only: time_scheme_c
    use stdlib_integrate_rk_tableaus, only: rk_tableau_t
    implicit none
    private

    type, public, extends(time_scheme_c) :: time_scheme_ERK_t
        !! Explicit Runge-Kutta scheme for an arbitrary Butcher tableau.
        !! Built in place by `rk_scheme(scheme, tableau, dt, model, status)`, or
        !! by name for a method without a hardcoded step.
        type(rk_tableau_t)                   :: tableau
        class(ODEmodel_c), pointer           :: model => null()  !! *reference* to model (not owned by us)
        class(scratch_space_c), allocatable  :: scratch
        class(ODEmodel_state_c), allocatable :: stage_state       !! reusable per-stage trial state buffer
        type(tendency_box_t), allocatable    :: k(:)              !! one tendency buffer per stage
        ! What the step needs from the tableau, computed once by `init` (so a
        ! step allocates nothing and rescans nothing); `tableau` must not be
        ! changed after `init`.
        logical,  allocatable                :: used(:)           !! stage i's tendency is read by b or a later row of A
        logical,  allocatable                :: combines(:)       !! stage i builds a stage state (a nonzero A(i,1:i-1))
        integer                              :: last = 0          !! the last stage evaluated
        real(wp), allocatable                :: coef(:)           !! dt*A(i,:) or dt*b, filled before each combine/accumulate
    contains
        procedure :: init    => ERK_init
        procedure :: advance => ERK_advance
    end type

    interface

        !> Build the scheme for `tableau`, time step `dt` and `model`, in place.
        !> Any previous content of `scheme` is discarded. `model` is kept by
        !> reference: it must be a `target` and outlive the scheme. An invalid
        !> tableau (no stages, missing or mis-sized coefficients, not explicit)
        !> or a model allocator that allocates nothing is a failure: it is
        !> returned in `status`, or reported and stops the program if `status`
        !> is absent; the scheme must not be used after a failure.
        module subroutine ERK_init(scheme, tableau, dt, model, status)
            class(time_scheme_ERK_t), intent(out)             :: scheme
            type(rk_tableau_t),       intent(in)              :: tableau !! Butcher tableau defining the method
            real(wp),                 intent(in)              :: dt      !! time step
            class(ODEmodel_c),        intent(in),    target   :: model   !! model
            type(status_t),           intent(inout), optional :: status  !! failure, if present
        end subroutine ERK_init

        module subroutine ERK_advance(scheme, time, state, status)
            class(time_scheme_ERK_t), intent(inout), target   :: scheme
            real(wp),                 intent(inout)           :: time
            class(ODEmodel_state_c),  intent(inout), target   :: state
            type(status_t),           intent(inout), optional :: status
        end subroutine ERK_advance

    end interface

end module stdlib_integrate_rk_generic
