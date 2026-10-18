module test_ODEmodel_other_state
    !! A second state type, to check that copy rejects a source of the wrong
    !! type.
    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: status_t, ODEmodel_state_c, ODEmodel_tendencies_c
    implicit none
    private

    type, extends(ODEmodel_state_c), public :: other_state_t
    contains
        procedure :: update => other_update
        procedure :: copy   => other_copy
    end type

contains

    subroutine other_update(state, status, scale, tendencies)
        class(other_state_t),         intent(inout) :: state
        class(status_t),              intent(inout) :: status
        real(wp),                     intent(in)    :: scale
        class(ODEmodel_tendencies_c), intent(in)    :: tendencies
    end subroutine other_update

    subroutine other_copy(self, status, source)
        class(other_state_t),    intent(inout) :: self
        class(status_t),         intent(inout) :: status
        class(ODEmodel_state_c), intent(in)    :: source
    end subroutine other_copy

end module test_ODEmodel_other_state
