module test_ODEmodel_boundaries_model
    !! Periodic 1-D diffusion, dy_i/dt = nu * (y(i-1) - 2 y(i) + y(i+1)), on N
    !! interior points stored with one halo point on each side, y(0:N+1). The
    !! stencil reads the halos, so they must be filled before every evaluation.
    !! `mode` selects who fills them:
    !!   MODE_HOOK  `apply_boundaries` fills the state's halos; `compute_tendencies` reads
    !!              the state directly and raises an error if a halo is not
    !!              filled, or if it evaluates a state at a time other than
    !!              the one its boundaries were applied at
    !!   MODE_COPY  `apply_boundaries` does nothing; `compute_tendencies` copies the state
    !!              into the scratch space and fills the copy's halos (what a
    !!              model does without the hook)
    !!   MODE_FAIL  `apply_boundaries` raises an error
    !! `n_boundaries` counts the calls to `apply_boundaries`. The tendencies' halos are set
    !! to GARBAGE, so every state a scheme builds has wrong halos until its boundaries are applied.
    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: stage_t, ODEmodel_c, ODEmodel_state_c, ODEmodel_tendencies_c, scratch_space_c, &
                                   & status_t, type_mismatch
    implicit none
    private

    integer,  parameter, public :: N = 8
    integer,  parameter, public :: MODE_HOOK = 1, MODE_COPY = 2, MODE_FAIL = 3
    real(wp), parameter, public :: GARBAGE = 1.0e3_wp

    type, extends(ODEmodel_c), public :: halo_model_t
        real(wp) :: nu = 1.0_wp
        integer  :: mode = MODE_HOOK
        integer  :: n_boundaries = 0
        real(wp) :: boundaries_at = -huge(1.0_wp)   !! time of the last `apply_boundaries`
    contains
        procedure :: allocate_state      => model_allocate_state
        procedure :: allocate_scratch    => model_allocate_scratch
        procedure :: allocate_tendencies => model_allocate_tendencies
        procedure :: compute_tendencies  => model_compute_tendencies
        procedure :: apply_boundaries            => model_apply_boundaries
    end type

    type, extends(ODEmodel_state_c), public :: halo_state_t
        real(wp), allocatable :: y(:)              !! y(0:N+1): interior 1..N, halos 0 and N+1
    contains
        procedure :: update   => state_update
        procedure :: copy     => state_copy
    end type

    type, extends(ODEmodel_tendencies_c), public :: halo_tendencies_t
        real(wp), allocatable :: dy(:)             !! dy(0:N+1); halos are GARBAGE
    contains
    end type

    type, extends(scratch_space_c) :: halo_scratch_t
        real(wp), allocatable :: w(:)              !! MODE_COPY: haloed working copy of y
    end type

contains

    subroutine model_allocate_state(model, state)
        class(halo_model_t),                  intent(in)  :: model
        class(ODEmodel_state_c), allocatable, intent(out) :: state

        allocate(halo_state_t :: state)
        select type (state)
        type is (halo_state_t)
            allocate(state%y(0:N+1), source=0.0_wp)
        end select
    end subroutine model_allocate_state

    subroutine model_allocate_scratch(model, scratch)
        class(halo_model_t),                 intent(in)  :: model
        class(scratch_space_c), allocatable, intent(out) :: scratch

        allocate(halo_scratch_t :: scratch)
        select type (scratch)
        type is (halo_scratch_t)
            allocate(scratch%w(0:N+1), source=0.0_wp)
        end select
    end subroutine model_allocate_scratch

    subroutine model_allocate_tendencies(model, tendencies)
        class(halo_model_t),                       intent(in)  :: model
        class(ODEmodel_tendencies_c), allocatable, intent(out) :: tendencies

        allocate(halo_tendencies_t :: tendencies)
        select type (tendencies)
        type is (halo_tendencies_t)
            allocate(tendencies%dy(0:N+1), source=0.0_wp)
        end select
    end subroutine model_allocate_tendencies

    !> Periodic halos: y(0) = y(N), y(N+1) = y(1).
    pure subroutine fill_halos(y)
        real(wp), intent(inout) :: y(0:)

        y(0)   = y(N)
        y(N+1) = y(1)
    end subroutine fill_halos

    subroutine model_apply_boundaries(model, status, state, time, stage)
        class(halo_model_t),     intent(inout) :: model
        class(status_t),         intent(inout) :: status
        class(ODEmodel_state_c), intent(inout) :: state
        real(wp),                intent(in)    :: time
        type(stage_t),           intent(in)    :: stage

        model%n_boundaries = model%n_boundaries + 1
        model%boundaries_at = time
        select case (model%mode)
        case (MODE_FAIL)
            call status%raise("halo_model: apply_boundaries failed")
        case (MODE_HOOK)
            select type (state)
            type is (halo_state_t)
                call fill_halos(state%y)
            class default
                call status%raise(type_mismatch("halo_model", "model_apply_boundaries", "state", "halo_state_t"))
            end select
        end select
    end subroutine model_apply_boundaries

    subroutine model_compute_tendencies(model, status, state, time, scratch, tendencies, stage)
        class(halo_model_t),          intent(inout) :: model
        class(status_t),              intent(inout) :: status
        real(wp),                     intent(in)    :: time
        class(ODEmodel_state_c),      intent(in)    :: state
        class(scratch_space_c),       intent(inout) :: scratch
        class(ODEmodel_tendencies_c), intent(inout) :: tendencies
        type(stage_t),                intent(in)    :: stage

        select type (state)
        type is (halo_state_t)
            select type (scratch)
            type is (halo_scratch_t)
                select type (tendencies)
                type is (halo_tendencies_t)
                    if (model%mode == MODE_COPY) then
                        scratch%w = state%y
                        call fill_halos(scratch%w)
                        call laplacian(model%nu, scratch%w, tendencies%dy)
                    else
                        if (state%y(0) /= state%y(N) .or. state%y(N+1) /= state%y(1)) then
                            call status%raise("halo_model: evaluated a state whose halos are not filled")
                        else if (time /= model%boundaries_at) then
                            call status%raise("halo_model: evaluated at a time other than its completion")
                        else
                            call laplacian(model%nu, state%y, tendencies%dy)
                        end if
                    end if
                    tendencies%dy(0)   = GARBAGE
                    tendencies%dy(N+1) = GARBAGE
                end select
            end select
        end select
    end subroutine model_compute_tendencies

    pure subroutine laplacian(nu, y, dy)
        real(wp), intent(in)    :: nu
        real(wp), intent(in)    :: y(0:)
        real(wp), intent(inout) :: dy(0:)

        dy(1:N) = nu * (y(0:N-1) - 2*y(1:N) + y(2:N+1))
    end subroutine laplacian

    subroutine state_update(state, status, scale, tendencies)
        class(halo_state_t),          intent(inout) :: state
        class(status_t),              intent(inout) :: status
        real(wp),                     intent(in)    :: scale
        class(ODEmodel_tendencies_c), intent(in)    :: tendencies

        select type (tendencies)
        type is (halo_tendencies_t)
            state%y = state%y + scale * tendencies%dy
        class default
            call status%raise(type_mismatch("halo_model", "state_update", "tendencies", "halo_tendencies_t"))
        end select
    end subroutine state_update

    subroutine state_copy(self, status, source)
        class(halo_state_t),     intent(inout) :: self
        class(status_t),         intent(inout) :: status
        class(ODEmodel_state_c), intent(in)    :: source

        select type (source)
        type is (halo_state_t)
            self%y = source%y
        class default
            call status%raise(type_mismatch("halo_model", "state_copy", "source", "halo_state_t"))
        end select
    end subroutine state_copy

end module test_ODEmodel_boundaries_model
