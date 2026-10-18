module test_RungeKutta_arrays_model
    !! A nonlinear, time-dependent right-hand side on a plain array, `rhs`, and
    !! the same right-hand side as a model on the library's classes, so the
    !! array schemes and the class schemes can be compared step for step:
    !!   dy(i)/dt = -(1 + i/10) y(i) + sin(y(i+1)) / 2 + cos(t)   (i+1 wraps around)
    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: stage_t, ODEmodel_c, ODEmodel_state_c, ODEmodel_tendencies_c, scratch_space_c, &
                                   & status_t, type_mismatch
    implicit none
    private

    public :: rhs

    type, extends(ODEmodel_c), public :: array_model_t
        integer :: n = 0
    contains
        procedure :: allocate_state      => model_allocate_state
        procedure :: allocate_scratch    => model_allocate_scratch
        procedure :: allocate_tendencies => model_allocate_tendencies
        procedure :: compute_tendencies  => model_compute_tendencies
    end type

    type, extends(ODEmodel_state_c), public :: array_state_t
        real(wp), allocatable :: y(:)
    contains
        procedure :: update   => state_update
        procedure :: copy     => state_copy
    end type

    type, extends(ODEmodel_tendencies_c), public :: array_tendencies_t
        real(wp), allocatable :: dy(:)
    end type

contains

    subroutine rhs(t, y, dydt, stage)
        real(wp), intent(in)  :: t
        real(wp), intent(in)  :: y(:)
        real(wp), intent(out) :: dydt(:)
        type(stage_t), intent(in) :: stage
        integer :: i, n

        n = size(y)
        do i = 1, n
            dydt(i) = -(1.0_wp + i/10.0_wp) * y(i) + sin(y(mod(i, n) + 1)) / 2 + cos(t)
        end do
    end subroutine rhs

    subroutine model_allocate_state(model, state)
        class(array_model_t),                 intent(in)  :: model
        class(ODEmodel_state_c), allocatable, intent(out) :: state

        allocate(array_state_t :: state)
        select type (state)
        type is (array_state_t)
            allocate(state%y(model%n), source=0.0_wp)
        end select
    end subroutine model_allocate_state

    subroutine model_allocate_scratch(model, scratch)
        class(array_model_t),                intent(in)  :: model
        class(scratch_space_c), allocatable, intent(out) :: scratch

        allocate(scratch_space_c :: scratch)
    end subroutine model_allocate_scratch

    subroutine model_allocate_tendencies(model, tendencies)
        class(array_model_t),                      intent(in)  :: model
        class(ODEmodel_tendencies_c), allocatable, intent(out) :: tendencies

        allocate(array_tendencies_t :: tendencies)
        select type (tendencies)
        type is (array_tendencies_t)
            allocate(tendencies%dy(model%n), source=0.0_wp)
        end select
    end subroutine model_allocate_tendencies

    subroutine model_compute_tendencies(model, status, state, time, scratch, tendencies, stage)
        class(array_model_t),         intent(inout) :: model
        class(status_t),              intent(inout) :: status
        real(wp),                     intent(in)    :: time
        class(ODEmodel_state_c),      intent(in)    :: state
        class(scratch_space_c),       intent(inout) :: scratch
        class(ODEmodel_tendencies_c), intent(inout) :: tendencies
        type(stage_t),                intent(in)    :: stage

        select type (state)
        type is (array_state_t)
            select type (tendencies)
            type is (array_tendencies_t)
                call rhs(time, state%y, tendencies%dy, stage)
            end select
        end select
    end subroutine model_compute_tendencies

    subroutine state_update(state, status, scale, tendencies)
        class(array_state_t),         intent(inout) :: state
        class(status_t),              intent(inout) :: status
        real(wp),                     intent(in)    :: scale
        class(ODEmodel_tendencies_c), intent(in)    :: tendencies

        select type (tendencies)
        type is (array_tendencies_t)
            state%y = state%y + scale * tendencies%dy
        class default
            call status%raise(type_mismatch("arrays_model", "state_update", "tendencies", "array_tendencies_t"))
        end select
    end subroutine state_update

    subroutine state_copy(self, status, source)
        class(array_state_t),    intent(inout) :: self
        class(status_t),         intent(inout) :: status
        class(ODEmodel_state_c), intent(in)    :: source

        select type (source)
        type is (array_state_t)
            self%y = source%y
        class default
            call status%raise(type_mismatch("arrays_model", "state_copy", "source", "array_state_t"))
        end select
    end subroutine state_copy

end module test_RungeKutta_arrays_model
