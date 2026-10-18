module stdlib_Integrate
    !! Top-level interface module: re-exports the public entities of the status,
    !! model and time-scheme modules that make up stdlib_Integrate.
    use stdlib_integrate_status
    use stdlib_integrate_stage
    use stdlib_integrate_descriptor
    use stdlib_integrate_model
    use stdlib_integrate_scheme
    use stdlib_integrate_rk_tableaus
    use stdlib_integrate_rk_steps
    use stdlib_integrate_rk_arrays
    implicit none
    private

    ! stdlib_integrate_status
    public :: status_t, type_mismatch

    ! stdlib_integrate_stage, stdlib_integrate_descriptor
    public :: stage_t
    public :: method_descriptor_t

    ! stdlib_integrate_model
    public :: ODEmodel_c, ODEmodel_state_c, ODEmodel_tendencies_c, scratch_space_c, tendency_box_t

    ! stdlib_integrate_scheme
    public :: time_scheme_c

    ! stdlib_integrate_rk_tableaus
    public :: rk_tableau_t, rk_tableau, rk_method_names

    ! stdlib_integrate_rk_steps (model classes) and stdlib_integrate_rk_arrays (plain arrays)
    public :: rk_scheme
    public :: rk_array_scheme_t, rk_rhs_i, rk_boundaries_i

end module stdlib_Integrate
