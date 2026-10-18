module heat_model
    !! The heat equation u_t = u_xx on [0, 1], u = 0 at both ends, by the method
    !! of lines on the model classes: second-order central differences on the
    !! nodes x_j = j dx, j = 0, ..., n + 1, dx = 1/(n + 1). The state holds all
    !! nodes; the two boundary nodes are filled by `apply_boundaries` (u = 0), so the
    !! stencil reads a state whose boundaries are in place.
    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: stage_t, ODEmodel_c, ODEmodel_state_c, ODEmodel_tendencies_c, scratch_space_c, &
                                   & status_t, type_mismatch
    implicit none
    private

    type, extends(ODEmodel_c), public :: heat_model_t
        integer  :: n = 0                     !! interior nodes
        real(wp) :: dx = 0.0_wp
    contains
        procedure :: allocate_state      => model_allocate_state
        procedure :: allocate_scratch    => model_allocate_scratch
        procedure :: allocate_tendencies => model_allocate_tendencies
        procedure :: compute_tendencies  => model_compute_tendencies
        procedure :: apply_boundaries            => model_apply_boundaries
    end type

    type, extends(ODEmodel_state_c), public :: heat_state_t
        real(wp), allocatable :: u(:)         !! u(0:n+1); u(0), u(n+1) are the boundary nodes
    contains
        procedure :: update   => state_update
        procedure :: copy     => state_copy
    end type

    type, extends(ODEmodel_tendencies_c) :: heat_tendencies_t
        real(wp), allocatable :: du(:)        !! du(0:n+1), zero at the boundary nodes
    contains
    end type

    public :: make_heat_model

contains

    function make_heat_model(n) result(model)
        integer, intent(in) :: n
        type(heat_model_t) :: model

        model%n = n
        model%dx = 1.0_wp / (n + 1)
    end function make_heat_model

    subroutine model_allocate_state(model, state)
        class(heat_model_t),                  intent(in)  :: model
        class(ODEmodel_state_c), allocatable, intent(out) :: state

        allocate(heat_state_t :: state)
        select type (state)
        type is (heat_state_t)
            allocate(state%u(0:model%n+1), source=0.0_wp)
        end select
    end subroutine model_allocate_state

    subroutine model_allocate_scratch(model, scratch)
        class(heat_model_t),                 intent(in)  :: model
        class(scratch_space_c), allocatable, intent(out) :: scratch

        allocate(scratch_space_c :: scratch)
    end subroutine model_allocate_scratch

    subroutine model_allocate_tendencies(model, tendencies)
        class(heat_model_t),                       intent(in)  :: model
        class(ODEmodel_tendencies_c), allocatable, intent(out) :: tendencies

        allocate(heat_tendencies_t :: tendencies)
        select type (tendencies)
        type is (heat_tendencies_t)
            allocate(tendencies%du(0:model%n+1), source=0.0_wp)
        end select
    end subroutine model_allocate_tendencies

    !> Dirichlet boundary condition: the boundary nodes are zero.
    subroutine model_apply_boundaries(model, status, state, time, stage)
        class(heat_model_t),     intent(inout) :: model
        class(status_t),         intent(inout) :: status
        class(ODEmodel_state_c), intent(inout) :: state
        real(wp),                intent(in)    :: time
        type(stage_t),           intent(in)    :: stage

        select type (state)
        type is (heat_state_t)
            state%u(0) = 0.0_wp
            state%u(model%n + 1) = 0.0_wp
        class default
            call status%raise(type_mismatch("heat_model", "model_apply_boundaries", "state", "heat_state_t"))
        end select
    end subroutine model_apply_boundaries

    subroutine model_compute_tendencies(model, status, state, time, scratch, tendencies, stage)
        class(heat_model_t),          intent(inout) :: model
        class(status_t),              intent(inout) :: status
        real(wp),                     intent(in)    :: time
        class(ODEmodel_state_c),      intent(in)    :: state
        class(scratch_space_c),       intent(inout) :: scratch
        class(ODEmodel_tendencies_c), intent(inout) :: tendencies
        type(stage_t),                intent(in)    :: stage
        integer :: n

        n = model%n
        select type (state)
        type is (heat_state_t)
            select type (tendencies)
            type is (heat_tendencies_t)
                tendencies%du(1:n) = (state%u(0:n-1) - 2 * state%u(1:n) + state%u(2:n+1)) / model%dx**2
                tendencies%du(0) = 0.0_wp
                tendencies%du(n + 1) = 0.0_wp
            end select
        end select
    end subroutine model_compute_tendencies

    subroutine state_update(state, status, scale, tendencies)
        class(heat_state_t),          intent(inout) :: state
        class(status_t),              intent(inout) :: status
        real(wp),                     intent(in)    :: scale
        class(ODEmodel_tendencies_c), intent(in)    :: tendencies

        select type (tendencies)
        type is (heat_tendencies_t)
            state%u = state%u + scale * tendencies%du
        class default
            call status%raise(type_mismatch("heat_model", "state_update", "tendencies", "heat_tendencies_t"))
        end select
    end subroutine state_update

    subroutine state_copy(self, status, source)
        class(heat_state_t),     intent(inout) :: self
        class(status_t),         intent(inout) :: status
        class(ODEmodel_state_c), intent(in)    :: source

        select type (source)
        type is (heat_state_t)
            self%u = source%u
        class default
            call status%raise(type_mismatch("heat_model", "state_copy", "source", "heat_state_t"))
        end select
    end subroutine state_copy

end module heat_model
