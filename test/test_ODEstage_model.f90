module test_ODEstage_model
    !! Test fixture: scalar decay, dy/dt = -y, whose `compute_tendencies` and
    !! `apply_boundaries` (and the array-path `rhs` and `boundaries`) record every call:
    !! which of the two it was, the stage it received, its time and its state.
    !! test_ODEstage reads the record to check what each engine passes.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: stage_t, ODEmodel_c, ODEmodel_state_c, ODEmodel_tendencies_c, &
                                   & scratch_space_c, status_t
    implicit none
    private

    public :: record_t, calls, n_calls, forget_calls, rec_rhs, rec_boundaries

    !> One call: `kind` is 'f' (an evaluation) or 'c' (an `apply_boundaries`).
    type, public :: record_t
        character     :: kind = ' '
        type(stage_t) :: stage
        real(wp)      :: time = 0.0_wp
        real(wp)      :: y = 0.0_wp
    end type

    integer, parameter :: max_calls = 100
    type(record_t) :: calls(max_calls)
    integer        :: n_calls = 0

    type, extends(ODEmodel_c), public :: stage_model_t
    contains
        procedure :: allocate_state      => model_allocate_state
        procedure :: allocate_scratch    => model_allocate_scratch
        procedure :: allocate_tendencies => model_allocate_tendencies
        procedure :: compute_tendencies  => model_compute_tendencies
        procedure :: apply_boundaries            => model_apply_boundaries
    end type

    type, extends(ODEmodel_state_c), public :: stage_state_t
        real(wp) :: y = 0.0_wp
    contains
        procedure :: update => state_update
        procedure :: copy   => state_copy
    end type

    type, extends(ODEmodel_tendencies_c) :: stage_tendencies_t
        real(wp) :: dy(1) = 0.0_wp
    end type

contains

    subroutine forget_calls()
        n_calls = 0
    end subroutine forget_calls

    subroutine record(kind, stage, time, y)
        character,     intent(in) :: kind
        type(stage_t), intent(in) :: stage
        real(wp),      intent(in) :: time, y

        if (n_calls == max_calls) error stop "test_ODEstage_model: too many calls recorded"
        n_calls = n_calls + 1
        calls(n_calls) = record_t(kind, stage, time, y)
    end subroutine record

    ! ---- array path ----

    subroutine rec_rhs(t, y, dydt, stage)
        real(wp),      intent(in)  :: t
        real(wp),      intent(in)  :: y(:)
        real(wp),      intent(out) :: dydt(:)
        type(stage_t), intent(in)  :: stage

        call record('f', stage, t, y(1))
        dydt = -y
    end subroutine rec_rhs

    subroutine rec_boundaries(t, y, stage)
        real(wp),      intent(in)    :: t
        real(wp),      intent(inout) :: y(:)
        type(stage_t), intent(in)    :: stage

        call record('c', stage, t, y(1))
    end subroutine rec_boundaries

    ! ---- class path ----

    subroutine model_allocate_state(model, state)
        class(stage_model_t),                 intent(in)  :: model
        class(ODEmodel_state_c), allocatable, intent(out) :: state
        allocate(stage_state_t :: state)
    end subroutine model_allocate_state

    subroutine model_allocate_scratch(model, scratch)
        class(stage_model_t),                intent(in)  :: model
        class(scratch_space_c), allocatable, intent(out) :: scratch
        allocate(scratch_space_c :: scratch)
    end subroutine model_allocate_scratch

    subroutine model_allocate_tendencies(model, tendencies)
        class(stage_model_t),                      intent(in)  :: model
        class(ODEmodel_tendencies_c), allocatable, intent(out) :: tendencies
        allocate(stage_tendencies_t :: tendencies)
    end subroutine model_allocate_tendencies

    subroutine model_compute_tendencies(model, status, state, time, scratch, tendencies, stage)
        class(stage_model_t),         intent(inout) :: model
        class(status_t),              intent(inout) :: status
        real(wp),                     intent(in)    :: time
        class(ODEmodel_state_c),      intent(in)    :: state
        class(scratch_space_c),       intent(inout) :: scratch
        class(ODEmodel_tendencies_c), intent(inout) :: tendencies
        type(stage_t),                intent(in)    :: stage

        select type (state)
        type is (stage_state_t)
            call record('f', stage, time, state%y)
            select type (tendencies)
            type is (stage_tendencies_t)
                tendencies%dy(1) = -state%y
            end select
        end select
    end subroutine model_compute_tendencies

    subroutine model_apply_boundaries(model, status, state, time, stage)
        class(stage_model_t),    intent(inout) :: model
        class(status_t),         intent(inout) :: status
        class(ODEmodel_state_c), intent(inout) :: state
        real(wp),                intent(in)    :: time
        type(stage_t),           intent(in)    :: stage

        select type (state)
        type is (stage_state_t)
            call record('c', stage, time, state%y)
        end select
    end subroutine model_apply_boundaries

    subroutine state_update(state, status, scale, tendencies)
        class(stage_state_t),         intent(inout) :: state
        class(status_t),              intent(inout) :: status
        real(wp),                     intent(in)    :: scale
        class(ODEmodel_tendencies_c), intent(in)    :: tendencies

        select type (tendencies)
        type is (stage_tendencies_t)
            state%y = state%y + scale * tendencies%dy(1)
        end select
    end subroutine state_update

    subroutine state_copy(self, status, source)
        class(stage_state_t),    intent(inout) :: self
        class(status_t),         intent(inout) :: status
        class(ODEmodel_state_c), intent(in)    :: source

        select type (source)
        type is (stage_state_t)
            self%y = source%y
        end select
    end subroutine state_copy

end module test_ODEstage_model
