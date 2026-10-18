module stdlib_integrate_scheme
    !! Defines the abstract time-scheme interface used to solve initial-value
    !! problems described by classes from [[stdlib_integrate_model]], and the
    !! allocation of a scheme's buffers shared by the concrete schemes.

    use stdlib_kinds, only: wp => dp
    use stdlib_integrate_status, only: status_t
    use stdlib_integrate_descriptor, only: method_descriptor_t
    use stdlib_integrate_model, only: ODEmodel_c, ODEmodel_state_c, scratch_space_c, tendency_box_t
    implicit none
    private

    public :: allocate_buffers

    type, abstract, public :: time_scheme_c
        !! Abstract class to be extended by concrete time schemes. A concrete
        !! scheme is built in place, by `rk_scheme(scheme, ...)` or by its
        !! type-bound `init`, never by a function returning it: a function
        !! result would be copied, buffers included, on assignment.
        real(wp) :: dt = 0.0_wp
            !! The time step: of every step at a fixed step; of the next step
            !! for an adaptive scheme, which updates it after every step.
        type(method_descriptor_t) :: description
            !! The descriptor of the scheme's method, set by its `init`; see `descriptor`.
    contains
        procedure(advance_i), deferred :: advance
        !! Advance model state in time by one time step.
        procedure :: descriptor => scheme_descriptor
        !! The descriptor of the scheme's method (name, order, cost, ...),
        !! whatever its family; f_evals is this scheme's cost per step.
        procedure :: get_dt => scheme_get_dt
        !! The time step (the next one, for an adaptive scheme).
        procedure :: set_dt => scheme_set_dt
        !! Change the time step: of every following step at a fixed step; of
        !! the next one for an adaptive scheme (e.g. to end exactly at a time:
        !! `call scheme%set_dt(min(scheme%get_dt(), t_end - time))`).
    end type

    abstract interface

        subroutine advance_i(scheme, time, state, status)
            import time_scheme_c, ODEmodel_state_c, wp, status_t
            class(time_scheme_c),    intent(inout), target   :: scheme   !! time integration scheme
            real(wp),                intent(inout)           :: time     !! time at beginning (in) or end (out) of time step
            class(ODEmodel_state_c), intent(inout), target   :: state    !! model state at beginning (in) or end (out) of time step
            type(status_t),          intent(inout), optional :: status   !! to raise errors
        end subroutine advance_i

    end interface

contains

    pure type(method_descriptor_t) function scheme_descriptor(scheme) result(description)
        class(time_scheme_c), intent(in) :: scheme
        description = scheme%description
    end function scheme_descriptor

    pure real(wp) function scheme_get_dt(scheme) result(dt)
        class(time_scheme_c), intent(in) :: scheme
        dt = scheme%dt
    end function scheme_get_dt

    subroutine scheme_set_dt(scheme, dt)
        class(time_scheme_c), intent(inout) :: scheme
        real(wp),             intent(in)    :: dt
        scheme%dt = dt
    end subroutine scheme_set_dt

    !> Allocate a scheme's buffers through the model: its scratch space,
    !> `n_registers` tendency registers in `k`, and, if `stage_state` is
    !> present, a stage state. Raises an error in `status` if one of the
    !> model's allocators returns without allocating.
    subroutine allocate_buffers(model, status, scratch, n_registers, k, stage_state)
        class(ODEmodel_c),                    intent(in)              :: model
        type(status_t),                       intent(inout)           :: status
        class(scratch_space_c),  allocatable, intent(inout)           :: scratch
        integer,                              intent(in)              :: n_registers
        type(tendency_box_t),    allocatable, intent(inout)           :: k(:)
        class(ODEmodel_state_c), allocatable, intent(inout), optional :: stage_state

        integer :: i

        call model%allocate_scratch(scratch)
        if (.not. allocated(scratch)) then
            call status%raise("the model's allocate_scratch returned an unallocated scratch space")
            return
        end if
        if (present(stage_state)) then
            call model%allocate_state(stage_state)
            if (.not. allocated(stage_state)) then
                call status%raise("the model's allocate_state returned an unallocated state")
                return
            end if
        end if
        allocate(k(n_registers))
        do i = 1, n_registers
            call model%allocate_tendencies(k(i)%item)
            if (.not. allocated(k(i)%item)) then
                call status%raise("the model's allocate_tendencies returned unallocated tendencies")
                return
            end if
        end do
    end subroutine allocate_buffers

end module stdlib_integrate_scheme
