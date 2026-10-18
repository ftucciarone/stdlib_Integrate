program test_RungeKutta_generic
    !! The generic tableau engine (time_scheme_ERK_t, built by rk_scheme from a
    !! tableau) reproduces the expected convergence order of every method of
    !! the library on dy/dt = -lambda*y, at its linear order since the
    !! problem is linear. Exact order conditions: test_RungeKutta_tableaus.f90.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: status_t, time_scheme_c, rk_scheme, rk_tableau_t, rk_tableau, &
                                   & rk_method_names
    use test_decay_model, only: decay_model_t, decay_state_t
    implicit none

    real(wp), parameter :: lambda = 2.0_wp, t_end = 1.0_wp
    integer,  parameter :: n_coarse = 50

    type(decay_model_t), target :: model
    type(rk_tableau_t) :: tab
    character(len=:), allocatable :: names(:)
    logical :: all_passed
    integer :: m

    model%lambda = lambda
    all_passed = .true.

    names = rk_method_names()
    do m = 1, size(names)
        tab = rk_tableau(trim(names(m)))
        call check_order(tab%title, tab, expected_order=real(tab%linear_order, wp), tol=0.15_wp)
    end do

    if (.not. all_passed) error stop "generic engine convergence checks FAILED"
    print *, "All generic engine convergence checks PASSED"

contains

    subroutine check_order(label, tableau, expected_order, tol)
        character(*),       intent(in) :: label
        type(rk_tableau_t), intent(in) :: tableau
        real(wp),           intent(in) :: expected_order, tol

        class(time_scheme_c), allocatable, target :: scheme
        real(wp) :: err_coarse, err_fine, observed_order
        integer :: n1

        ! High orders reach round-off at n_coarse steps: measure them on coarser steps.
        n1 = n_coarse
        if (tableau%linear_order >= 5) n1 = 20
        call rk_scheme(scheme, tableau, t_end/n1, model)
        err_coarse = final_error(scheme, n1)
        call rk_scheme(scheme, tableau, t_end/(2*n1), model)
        err_fine = final_error(scheme, 2*n1)
        observed_order = log(err_coarse / err_fine) / log(2.0_wp)

        if (abs(observed_order - expected_order) > tol) then
            print '(a28, a, f6.3, a, f6.3, a)', label, ": FAIL (observed order ", observed_order, &
                & ", expected ", expected_order, ")"
            all_passed = .false.
        else
            print '(a28, a, f6.3)', label, ": OK, observed order ", observed_order
        end if
    end subroutine check_order

    !> Error at t_end after n_steps steps from y(0) = 1.
    real(wp) function final_error(scheme, n_steps) result(err)
        class(time_scheme_c), intent(inout), target :: scheme
        integer,              intent(in)            :: n_steps

        type(decay_state_t), target :: state
        type(status_t) :: status
        real(wp) :: time
        integer  :: i

        state%y = 1.0_wp
        time = 0.0_wp
        do i = 1, n_steps
            call scheme%advance(time, state, status)
            if (status%failed) error stop "time step failed"
        end do
        err = abs(state%y - exp(-lambda * t_end))
    end function final_error

end program test_RungeKutta_generic
