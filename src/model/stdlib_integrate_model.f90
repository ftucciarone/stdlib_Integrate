module stdlib_integrate_model
    !! Defines abstract classes representing models, model states and model tendencies.
    !! They are composed with the time schemes of stdlib_Integrate to solve initial-value problems.

    use stdlib_kinds, only: wp => dp

    use stdlib_integrate_status, only: status_t
    use stdlib_integrate_stage, only: stage_t

    implicit none
    private

    !============================== Model ===============================

    type, abstract, public :: ODEmodel_c
        !! Abstract class describing a model (initial-value problem, without the initial value).
        !! This class must be extended for a given model. It typically contains prescribed data
        !! (time-independent or time-dependent)
        !! needed to compute the time derivative of the model state.

    contains

        procedure(alloc_state_i), deferred :: allocate_state
            !! Allocates a `CLASS(ODEmodel_state_c)` object suitable for the model.
            !! Called by time schemes, typically once, at the time of their creation.

        procedure(alloc_scratch_i), deferred :: allocate_scratch
            !! Allocates a `CLASS(scratch_c)` object suitable for the model.
            !! Called by time schemes, typically once, at the time of their creation.

        procedure(alloc_tendencies_i), deferred :: allocate_tendencies
            !! Allocates a `CLASS(ODEmodel_tendencies_c)` object suitable for the model.
            !! Called by time schemes, typically a few times a the time of their creation.

        procedure(comp_tendencies_i), deferred :: compute_tendencies
            !! Given time, model state and scratch space, computes the model tendencies.
            !! Should check actual type of arguments via `SELECT TYPE`.
            !! Raises an error in `status` if arguments have the wrong types or another error occurs.
            !! Called by time schemes repeatedly during time integration.
            !! `model` is `INTENT(INOUT)` to allow reading external data during time integration.
            !! `stage` tells which stage of the step this evaluation belongs to (see stdlib_integrate_stage):
            !! a model may use it to compute some terms only at some stages.

        procedure :: apply_boundaries => default_apply_boundaries
            !! Makes a state ready to be read by `compute_tendencies`: fills its halos and
            !! applies its boundary conditions. Called by time schemes on every state they
            !! evaluate, before the evaluation: on the step's initial state at the start of
            !! each step, and on each stage state right after it is built. `stage` is the
            !! stage they are applied for. The default does nothing.

    end type

    !======================== Model state ===============================

    type, abstract, public :: ODEmodel_state_c
        !! Abstract class for the current state (degrees of freedom) of a model.
        !! This class must be extended for a given model or class of models sharing the same state space.
        !! It should contain only the time-dependent model state. Prescribed and time-independent data
        !! should be part of the model itself.

    contains

        procedure(update_i), deferred :: update
            !! Given a `CLASS(tendency_c)` and a scaling factor, adds the (scaled) tendencies to the model state.
            !! Should check actual type of arguments via `SELECT TYPE`
            !! and raise an error if arguments have the wrong types or another error occurs.
            !! Called by time schemes repeatedly during time integration.

        procedure(copy_i), deferred :: copy
            !! Copies another state of the same type into this one (`self = source`).
            !! Should check the actual type of `source` via `SELECT TYPE` and raise an error if it is wrong.
            !! Called by time schemes (and by the default `combine`) to build stage states.

        procedure :: combine => default_combine
            !! `self = x + sum(c(m) * k(m))`, the stage state of a Runge-Kutta step, `x` distinct from `self`.
            !! The default copies `x` and calls `update` once per nonzero coefficient.

        procedure :: accumulate => default_accumulate
            !! `self = self + sum(c(m) * k(m))`, in place: the final update of a step.
            !! The default calls `update` once per nonzero coefficient.

    end type

    !======================== Model tendencies ============================

    type, abstract, public :: ODEmodel_tendencies_c
        !! Abstract class for the tendencies (time derivatives) of a model state.
        !! This class must be extended for a given model or class of models sharing the same space for tendencies.
        !! This class is distinct from `ODEmodel_state_c` because it may contain fewer degrees of freedom (DOFs) if
        !! certain DOFs are known to have zero time derivatives.
    end type

    type, public :: tendency_box_t
        !! Wraps a polymorphic tendency object so an array of them can be stored:
        !! Fortran does not allow arrays of dissimilar CLASS(...) entities directly.
        class(ODEmodel_tendencies_c), allocatable :: item
    end type

    !======================== Scratch space ============================

    type, public :: scratch_space_c
        !! Type for pre-allocated scratch space used to compute model tendencies.
        !! This type should be extended for a given model or class of models sharing the same needs for scratch space.
        !! It can also be used as-is if no scratch space is needed.
    end type

    abstract interface

        subroutine alloc_state_i(model, state)
            import ODEmodel_c, ODEmodel_state_c
            class(ODEmodel_c),                    intent(in)  :: model !! model for which the state is allocated
            class(ODEmodel_state_c), allocatable, intent(out) :: state !! state object to be allocated
        end subroutine alloc_state_i

        subroutine alloc_scratch_i(model, scratch)
            import ODEmodel_c, scratch_space_c
            class(ODEmodel_c),                   intent(in)  :: model   !! model for which the scratch space is allocated
            class(scratch_space_c), allocatable, intent(out) :: scratch !! scratch space to be allocated
        end subroutine alloc_scratch_i

        subroutine alloc_tendencies_i(model, tendencies)
            import ODEmodel_c, ODEmodel_tendencies_c
            class(ODEmodel_c),                         intent(in)  :: model      !! model for which the tendencies are allocated
            class(ODEmodel_tendencies_c), allocatable, intent(out) :: tendencies !! tendencies object to be allocated
        end subroutine alloc_tendencies_i

        subroutine comp_tendencies_i(model, status, state, time, scratch, tendencies, stage)
            import ODEmodel_c, ODEmodel_state_c, scratch_space_c, ODEmodel_tendencies_c, wp, status_t, stage_t
            class(ODEmodel_c),            intent(inout) :: model      !! model
            class(status_t),              intent(inout) :: status     !! to raise an error
            real(wp),                     intent(in   ) :: time       !! current time
            class(ODEmodel_state_c),      intent(in   ) :: state      !! current model state
            class(scratch_space_c),       intent(inout) :: scratch    !! pre-allocated scratch space
            class(ODEmodel_tendencies_c), intent(inout) :: tendencies !! tendencies to be computed
            type(stage_t),                intent(in   ) :: stage      !! stage this evaluation belongs to
        end subroutine comp_tendencies_i

        subroutine update_i(state, status, scale, tendencies)
            import ODEmodel_state_c, ODEmodel_tendencies_c, wp, status_t
            class(ODEmodel_state_c),      intent(inout) :: state      !! model state
            class(status_t),              intent(inout) :: status     !! to raise an error
            real(wp),                     intent(in   ) :: scale      !! multiplicative factor
            class(ODEmodel_tendencies_c), intent(in   ) :: tendencies !! tendencies to be added to `state`, multiplied by `scale`
        end subroutine update_i

        subroutine copy_i(self, status, source)
            import ODEmodel_state_c, status_t
            class(ODEmodel_state_c), intent(inout) :: self   !! state to overwrite
            class(status_t),         intent(inout) :: status !! to raise an error
            class(ODEmodel_state_c), intent(in   ) :: source !! state to copy, of the same type
        end subroutine copy_i

    end interface

contains

    subroutine default_apply_boundaries(model, status, state, time, stage)
        !! No boundaries to apply: a model without halos or boundary conditions.
        class(ODEmodel_c),       intent(inout) :: model  !! model
        class(status_t),         intent(inout) :: status !! to raise an error
        class(ODEmodel_state_c), intent(inout) :: state  !! state, its boundaries applied in place
        real(wp),                intent(in   ) :: time   !! time of `state` (the stage time)
        type(stage_t),           intent(in   ) :: stage  !! stage they are applied for
    end subroutine default_apply_boundaries

    subroutine default_combine(self, status, x, c, k)
        !! `self = x + sum(c(m) * k(m))`, through `copy` and `update`. The
        !! `target` attributes are part of the interface an override must
        !! match, so that a one-pass default can be added without changing it.
        class(ODEmodel_state_c), intent(inout), target :: self   !! stage state to build
        class(status_t),         intent(inout)         :: status !! to raise an error
        class(ODEmodel_state_c), intent(in   ), target :: x      !! base state (not `self`)
        real(wp),                intent(in   )         :: c(:)   !! coefficients (dt included); zeros are skipped
        type(tendency_box_t),    intent(in   ), target :: k(:)   !! tendencies, size(k) == size(c)

        integer :: m

        call self%copy(status, x)
        do m = 1, size(c)
            if (status%failed) return
            if (c(m) /= 0.0_wp) call self%update(status, c(m), k(m)%item)
        end do
    end subroutine default_combine

    subroutine default_accumulate(self, status, c, k)
        !! `self = self + sum(c(m) * k(m))`, through `update`.
        class(ODEmodel_state_c), intent(inout), target :: self   !! state to update in place
        class(status_t),         intent(inout)         :: status !! to raise an error
        real(wp),                intent(in   )         :: c(:)   !! coefficients (dt included); zeros are skipped
        type(tendency_box_t),    intent(in   ), target :: k(:)   !! tendencies, size(k) == size(c)

        integer :: m

        do m = 1, size(c)
            if (status%failed) return
            if (c(m) /= 0.0_wp) call self%update(status, c(m), k(m)%item)
        end do
    end subroutine default_accumulate

end module stdlib_integrate_model
