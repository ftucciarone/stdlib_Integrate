module test_decay_model
    !! Test fixture: scalar exponential decay, dy/dt = -lambda*y, with a known
    !! exact solution. A copy of examples/decay/decay_model.f90, kept separate so
    !! the tests do not depend on example code (the example may evolve).

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: stage_t, ODEmodel_c, ODEmodel_state_c, ODEmodel_tendencies_c, &
                                   & scratch_space_c, status_t, type_mismatch
    implicit none
    private

    public :: decay_model_t, decay_state_t

    type, extends(ODEmodel_c), public :: decay_model_t
        real(wp) :: lambda = 1.0_wp
    contains
        procedure :: allocate_state      => decay_allocate_state
        procedure :: allocate_scratch    => decay_allocate_scratch
        procedure :: allocate_tendencies => decay_allocate_tendencies
        procedure :: compute_tendencies  => decay_compute_tendencies
    end type

    type, extends(ODEmodel_state_c), public :: decay_state_t
        real(wp) :: y = 0.0_wp
    contains
        procedure :: update   => decay_state_update
        procedure :: copy     => decay_state_copy
    end type

    type, extends(ODEmodel_tendencies_c) :: decay_tendencies_t
        real(wp) :: dy = 0.0_wp
    end type

contains

    subroutine decay_allocate_state(model, state)
        class(decay_model_t), intent(in) :: model
        class(ODEmodel_state_c), allocatable, intent(out) :: state
        allocate(decay_state_t :: state)
    end subroutine decay_allocate_state

    subroutine decay_allocate_scratch(model, scratch)
        class(decay_model_t), intent(in) :: model
        class(scratch_space_c), allocatable, intent(out) :: scratch
        allocate(scratch_space_c :: scratch)
    end subroutine decay_allocate_scratch

    subroutine decay_allocate_tendencies(model, tendencies)
        class(decay_model_t), intent(in) :: model
        class(ODEmodel_tendencies_c), allocatable, intent(out) :: tendencies
        allocate(decay_tendencies_t :: tendencies)
    end subroutine decay_allocate_tendencies

    subroutine decay_compute_tendencies(model, status, state, time, scratch, tendencies, stage)
        class(decay_model_t),         intent(inout) :: model
        class(status_t),              intent(inout) :: status
        real(wp),                     intent(in)    :: time
        class(ODEmodel_state_c),      intent(in)    :: state
        class(scratch_space_c),       intent(inout) :: scratch
        class(ODEmodel_tendencies_c), intent(inout) :: tendencies
        type(stage_t),                intent(in)    :: stage

        select type (state)
        type is (decay_state_t)
            select type (tendencies)
            type is (decay_tendencies_t)
                tendencies%dy = -model%lambda * state%y
            end select
        end select
    end subroutine decay_compute_tendencies

    subroutine decay_state_update(state, status, scale, tendencies)
        class(decay_state_t),         intent(inout) :: state
        class(status_t),              intent(inout) :: status
        real(wp),                     intent(in)    :: scale
        class(ODEmodel_tendencies_c), intent(in)    :: tendencies

        select type (tendencies)
        type is (decay_tendencies_t)
            state%y = state%y + scale * tendencies%dy
        end select
    end subroutine decay_state_update

    subroutine decay_state_copy(self, status, source)
        class(decay_state_t),    intent(inout) :: self
        class(status_t),         intent(inout) :: status
        class(ODEmodel_state_c), intent(in)    :: source

        select type (source)
        type is (decay_state_t)
            self%y = source%y
        class default
            call status%raise(type_mismatch("decay_model", "decay_state_copy", "source", "decay_state_t"))
        end select
    end subroutine decay_state_copy

end module test_decay_model
