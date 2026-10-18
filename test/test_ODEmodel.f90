program test_ODEmodel
    !! The state operations on a model that implements only update and copy
    !! (test_decay_model), so combine and accumulate are the defaults:
    !! results, skipped zero coefficients, and copy's type check; and a step
    !! failing inside a state operation returns its error in `status`.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: stage_t, status_t, scratch_space_c, tendency_box_t, time_scheme_c, &
                                   & rk_scheme, rk_tableau
    use test_decay_model, only: decay_model_t, decay_state_t
    use test_ODEmodel_other_state, only: other_state_t
    implicit none

    type(decay_model_t), target :: model
    class(time_scheme_c), allocatable :: scheme
    real(wp) :: time
    type(decay_state_t) :: x, y
    type(other_state_t) :: other
    type(tendency_box_t) :: k(3)
    class(scratch_space_c), allocatable :: scratch
    type(status_t) :: status
    real(wp), parameter :: tol = 1.0e-15_wp
    logical :: all_passed
    integer :: m

    all_passed = .true.
    model%lambda = 2.0_wp
    call model%allocate_scratch(scratch)

    ! Tendencies dy = -lambda * y of the states y = 1, 2, 3: k = [-2, -4, -6].
    do m = 1, 3
        call model%allocate_tendencies(k(m)%item)
        x%y = real(m, wp)
        call model%compute_tendencies(status, x, 0.0_wp, scratch, k(m)%item, stage_t())
    end do

    x%y = 10.0_wp
    call y%combine(status, x, [0.5_wp, 0.0_wp, 0.25_wp], k)
    call expect(.not. status%failed .and. abs(y%y - (10.0_wp - 1.0_wp - 1.5_wp)) < tol, "combine")
    call expect(x%y == 10.0_wp, "combine leaves x unchanged")

    call y%accumulate(status, [1.0_wp, 0.0_wp, 0.0_wp], k)
    call expect(.not. status%failed .and. abs(y%y - 5.5_wp) < tol, "accumulate")

    call y%copy(status, x)
    call expect(.not. status%failed .and. y%y == 10.0_wp, "copy")
    call y%copy(status, other)
    call expect(status%failed, "copy from another state type raises an error")

    ! A state of another type: the stage combine's copy raises a type
    ! mismatch, which advance returns in status (without status it would
    ! stop the program with that message).
    status = status_t()
    call rk_scheme(scheme, rk_tableau("rk4"), 0.1_wp, model)
    time = 0.0_wp
    call scheme%advance(time, other, status)
    call expect(status%failed, "advance on a state of the wrong type fails")
    if (allocated(status%message)) then
        call expect(index(status%message, "decay_state_copy") > 0, "message names the failing copy")
    else
        call expect(.false., "failure has a message")
    end if

    if (.not. all_passed) error stop "state operation checks FAILED"
    print *, "All state operation checks PASSED"

contains

    subroutine expect(condition, what)
        logical, intent(in) :: condition
        character(*), intent(in) :: what

        if (.not. condition) then
            print '(2a)', "FAILED: ", what
            all_passed = .false.
        end if
    end subroutine expect

end program test_ODEmodel
