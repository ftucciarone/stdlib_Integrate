module stdlib_integrate_stage
    !! Which stage of a time step an evaluation belongs to (`stage_t`). Every
    !! engine passes one to the model's `compute_tendencies` and `apply_boundaries`
    !! (class path) and to the right-hand side and `boundaries` (array path), so
    !! a model can vary what it computes from stage to stage:
    !!
    !!   - compute an expensive term only when `first` (on the step's initial
    !!     state), keep it in the scratch space and reuse it at later stages;
    !!     that term is then integrated by forward Euler (first order in time);
    !!   - for a stochastic forcing, draw the increment once per step, when
    !!     `first`, scaled with `sqrt(dt)`, and reuse it at every stage: all
    !!     stages of a step must see the same increment;
    !!   - compute diagnostics of the state at the start of the step, when
    !!     `first`: the state is then exactly the result of the previous step.
    !!
    !! Varying the right-hand side from stage to stage changes the method: the
    !! order of the tableau no longer applies, and is the model's
    !! responsibility. `last` marks the last evaluation of the step, which is on
    !! a stage state, not on the step's result: the result is never evaluated
    !! within the step that produces it.

    use stdlib_kinds, only: wp => dp
    implicit none
    private

    type, public :: stage_t
        !! The stage an evaluation belongs to. Passed `intent(in)`: read it, do
        !! not keep it.
        integer  :: index  = 0       !! stage of the method, 1..stages (stages nothing uses are skipped:
                                     !! dp54 evaluates 1..6 of 7)
        integer  :: stages = 0       !! number of stages of the method
        real(wp) :: c      = 0.0_wp  !! stage time fraction: the evaluation is at t_n + c*dt
        real(wp) :: dt     = 0.0_wp  !! the time step
        logical  :: first  = .false. !! first evaluation of the step: on the step's initial state, at t_n
        logical  :: last   = .false. !! last evaluation of the step (on a stage state, not on the result)
    end type

end module stdlib_integrate_stage
