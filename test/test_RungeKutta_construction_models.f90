module test_RungeKutta_construction_models
    !! A decay model whose allocators can be told to allocate nothing, to
    !! check that building a scheme reports a model that does not allocate
    !! its buffers.
    use stdlib_Integrate, only: ODEmodel_state_c, ODEmodel_tendencies_c, scratch_space_c
    use test_decay_model, only: decay_model_t
    implicit none
    private

    type, extends(decay_model_t), public :: faulty_model_t
        logical :: no_state = .false.       !! allocate_state allocates nothing
        logical :: no_tendencies = .false.  !! allocate_tendencies allocates nothing
        logical :: no_scratch = .false.     !! allocate_scratch allocates nothing
    contains
        procedure :: allocate_state      => faulty_allocate_state
        procedure :: allocate_tendencies => faulty_allocate_tendencies
        procedure :: allocate_scratch    => faulty_allocate_scratch
    end type

contains

    subroutine faulty_allocate_state(model, state)
        class(faulty_model_t), intent(in) :: model
        class(ODEmodel_state_c), allocatable, intent(out) :: state
        if (.not. model%no_state) call model%decay_model_t%allocate_state(state)
    end subroutine faulty_allocate_state

    subroutine faulty_allocate_tendencies(model, tendencies)
        class(faulty_model_t), intent(in) :: model
        class(ODEmodel_tendencies_c), allocatable, intent(out) :: tendencies
        if (.not. model%no_tendencies) call model%decay_model_t%allocate_tendencies(tendencies)
    end subroutine faulty_allocate_tendencies

    subroutine faulty_allocate_scratch(model, scratch)
        class(faulty_model_t), intent(in) :: model
        class(scratch_space_c), allocatable, intent(out) :: scratch
        if (.not. model%no_scratch) call model%decay_model_t%allocate_scratch(scratch)
    end subroutine faulty_allocate_scratch

end module test_RungeKutta_construction_models
