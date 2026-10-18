module stdlib_integrate_descriptor
    !! The descriptor of a time-integration method, whatever its family: what
    !! describes it (name, order, error estimate, SSP coefficient, cost).
    !! The vocabulary of the library:
    !!   - a method is the mathematical object (RK4, Dormand-Prince 5(4));
    !!   - its descriptor, `method_descriptor_t`, is what describes it;
    !!   - its definition is the descriptor plus the coefficients that fully
    !!     define it: for an explicit Runge-Kutta method, its Butcher tableau
    !!     `rk_tableau_t`, which extends the descriptor;
    !!   - a scheme (`time_scheme_c`) is the method running on a model, with its
    !!     buffers and time step; `scheme%descriptor()` returns its descriptor.
    !!
    !! Data only: nothing here runs a method (a scheme does), and a field
    !! belongs here only if it means the same for every family.

    use stdlib_kinds, only: wp => dp
    implicit none
    private

    type, public :: method_descriptor_t
        !! Descriptor of a time-integration method. The values of a method
        !! of the library are set with its coefficients
        !! (stdlib_integrate_rk_tableaus); a method built by hand may leave
        !! the optional ones at their defaults (0: unknown or none).
        !> Short name ("rk4"), title ("Classic RK4") and source of the coefficients
        character(:), allocatable :: name
        character(:), allocatable :: title
        character(:), allocatable :: reference
        !> Order for general (nonlinear) problems
        integer  :: order = 0
        !> Order of the embedded error estimate; 0 if the method has none
        integer  :: order_hat = 0
        !> No equation to solve for the new state
        logical  :: explicit = .true.
        !> SSP coefficient C: strong-stability preserving for dt <= C*dt of
        !> forward Euler; 0 if not SSP
        real(wp) :: ssp_coef = 0.0_wp
        !> Tendency evaluations per step: the method's cost, comparable
        !> across families
        integer  :: f_evals = 0
    contains
        procedure, non_overridable :: ssp         => method_ssp
    end type

contains

    !> True if the method is strong-stability preserving (ssp_coef > 0).
    pure logical function method_ssp(self)
        class(method_descriptor_t), intent(in) :: self

        method_ssp = self%ssp_coef > 0.0_wp
    end function method_ssp

end module stdlib_integrate_descriptor
