program test_ODEmethod
    !! The descriptor of a method (method_descriptor_t), as its definition (the
    !! tableau) holds it and as every scheme returns it with `descriptor()`:
    !!   - a tableau's parent component is what its method's schemes return
    !!     (by name, generic tableau engine, plain arrays by name and by tableau);
    !!   - the values: tendency evaluations per step and SSP coefficients of
    !!     the methods of the library (the literature's), and ssp();
    !!   - the generic engines count their own cost, also for a tableau built
    !!     by hand, which leaves the rest of the description at its defaults.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: status_t, time_scheme_c, rk_scheme, rk_tableau, rk_tableau_t, &
                                   & rk_method_names, rk_array_scheme_t, method_descriptor_t
    use test_decay_model, only: decay_model_t
    implicit none

    character(*), parameter :: expected_names(3) = [character(8) :: "rk3ws", "rk4", "dp54"]
    integer, parameter :: expected_f_evals(3) = [3, 4, 6]
    real(wp), parameter :: expected_ssp(3) = [0, 0, 0]       ! no SSP method in the library

    type(decay_model_t), target :: model
    class(time_scheme_c), allocatable :: scheme
    type(rk_array_scheme_t) :: array_scheme
    type(rk_tableau_t) :: tab, by_hand
    type(method_descriptor_t) :: descriptor
    type(status_t) :: status
    character(len=:), allocatable :: names(:)
    logical :: all_passed
    integer :: m

    all_passed = .true.
    names = rk_method_names()
    call expect(size(names) == size(expected_names), "the table has the methods this test knows")

    do m = 1, min(size(names), size(expected_names))
        tab = rk_tableau(trim(names(m)))
        call expect(trim(names(m)) == trim(expected_names(m)), trim(names(m)) // ": table order")
        call expect(tab%f_evals == expected_f_evals(m), trim(names(m)) // ": f_evals")
        call expect(tab%ssp_coef == expected_ssp(m), trim(names(m)) // ": ssp_coef")
        call expect(tab%ssp() .eqv. expected_ssp(m) > 0, trim(names(m)) // ": ssp()")
        call expect(tab%explicit, trim(names(m)) // ": explicit")

        status = status_t()
        call rk_scheme(scheme, trim(names(m)), 0.1_wp, model, status)
        call expect(.not. status%failed, trim(names(m)) // ": scheme by name built")
        if (.not. status%failed) call expect_same(scheme%descriptor(), tab%method_descriptor_t, trim(names(m)) // ": by name")

        status = status_t()
        call rk_scheme(scheme, tab, 0.1_wp, model, status)
        call expect(.not. status%failed, trim(names(m)) // ": generic engine built")
        if (.not. status%failed) call expect_same(scheme%descriptor(), tab%method_descriptor_t, trim(names(m)) // ": generic")

        status = status_t()
        call rk_scheme(array_scheme, trim(names(m)), 3, 0.1_wp, status)
        call expect(.not. status%failed, trim(names(m)) // ": array scheme built")
        if (.not. status%failed) call expect_same(array_scheme%descriptor(), tab%method_descriptor_t, trim(names(m)) // ": arrays")

        status = status_t()
        call rk_scheme(array_scheme, tab, 3, 0.1_wp, status)
        call expect(.not. status%failed, trim(names(m)) // ": array tableau loop built")
        if (.not. status%failed) call expect_same(array_scheme%descriptor(), tab%method_descriptor_t, &
                                                  trim(names(m)) // ": array tableau loop")
    end do

    ! ---- a tableau built by hand: the engine counts its cost, the rest is unknown ----
    by_hand%stages = 3
    by_hand%A = reshape([0.0_wp, 0.5_wp, 0.0_wp,  0.0_wp, 0.0_wp, 1.0_wp,  0.0_wp, 0.0_wp, 0.0_wp], [3, 3])
    by_hand%b = [0.0_wp, 1.0_wp, 0.0_wp]     ! the midpoint method, plus a third stage nothing uses
    by_hand%c = [0.0_wp, 0.5_wp, 1.0_wp]
    status = status_t()
    call rk_scheme(scheme, by_hand, 0.1_wp, model, status)
    descriptor = scheme%descriptor()
    call expect(.not. status%failed .and. descriptor%f_evals == 2, "by hand: the engine evaluates 2 of 3 stages")
    call expect(descriptor%order == 0 .and. descriptor%ssp_coef == 0.0_wp .and. .not. allocated(descriptor%name), &
                "by hand: the description is left at its defaults")
    status = status_t()
    call rk_scheme(array_scheme, by_hand, 3, 0.1_wp, status)
    descriptor = array_scheme%descriptor()
    call expect(.not. status%failed .and. descriptor%f_evals == 2, "by hand, arrays: the loop evaluates 2 of 3 stages")

    if (.not. all_passed) error stop "method descriptor checks FAILED"
    print *, "All method descriptor checks PASSED"

contains

    subroutine expect(condition, what)
        logical, intent(in) :: condition
        character(*), intent(in) :: what

        if (.not. condition) then
            print '(2a)', "FAILED: ", what
            all_passed = .false.
        end if
    end subroutine expect

    subroutine expect_same(got, want, what)
        type(method_descriptor_t), intent(in) :: got, want
        character(*), intent(in) :: what

        call expect(allocated(got%name) .and. allocated(got%title) .and. allocated(got%reference), what // ": named")
        if (.not. allocated(got%name)) return
        call expect(got%name == want%name .and. got%title == want%title .and. got%reference == want%reference &
                    .and. got%order == want%order .and. got%order_hat == want%order_hat &
                    .and. (got%explicit .eqv. want%explicit) .and. got%ssp_coef == want%ssp_coef &
                    .and. got%f_evals == want%f_evals, what // ": same description as the tableau")
    end subroutine expect_same

end program test_ODEmethod
