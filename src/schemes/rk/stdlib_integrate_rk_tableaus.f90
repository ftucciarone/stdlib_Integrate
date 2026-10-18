module stdlib_integrate_rk_tableaus
    !! Butcher tableaux of the explicit Runge-Kutta methods of the library,
    !! looked up by name with `rk_tableau(name)`; `rk_method_names()` lists
    !! them. Adding a method is adding a constructor here and its name to the
    !! lookup: the generic engine runs any tableau.
    !!
    !! Methods: rk3ws, rk4, dp54.

    use stdlib_kinds, only: wp => dp
    use stdlib_integrate_descriptor, only: method_descriptor_t
    use stdlib_integrate_status, only: status_t
    implicit none
    private

    type, public, extends(method_descriptor_t) :: rk_tableau_t
        !! The definition of an s-stage explicit RK method: its descriptor
        !! (name, order, order_hat of b_hat, ssp_coef, f_evals: the parent,
        !! `tab%method_descriptor_t`) plus its Butcher tableau (A, b, c).
        !! `A` is strictly lower triangular. A tableau can also be built by
        !! hand (stages, A, b, c; the rest is optional).
        integer                   :: stages = 0
        real(wp), allocatable     :: A(:,:)
        real(wp), allocatable     :: b(:)
        real(wp), allocatable     :: c(:)
        !> Embedded weights for an error estimate (of order order_hat); unallocated if none
        real(wp), allocatable     :: b_hat(:)
        !> Order for linear autonomous problems (>= order)
        integer                   :: linear_order = 0
        !> Last stage = first stage of the next step (A(s,:) = b)
        logical                   :: fsal = .false.
    end type

    public :: rk_tableau
    public :: rk_method_names
    public :: check_tableau, stage_unused

contains

    !> Tableau of the method called `name` (see rk_method_names). For an
    !> unknown name, `found` is false and the result has stages = 0.
    function rk_tableau(name, found) result(tab)
        character(*), intent(in)       :: name
        logical, intent(out), optional :: found
        type(rk_tableau_t)             :: tab

        select case (name)
        case ("rk3ws")
            tab = rk3ws_tableau()
        case ("rk4")
            tab = rk4_tableau()
        case ("dp54")
            tab = dp54_tableau()
        end select
        if (present(found)) found = tab%stages > 0
    end function rk_tableau

    !> Names of all methods, in table order.
    pure function rk_method_names() result(names)
        character(len=5) :: names(3)

        names = [character(len=5) :: "rk3ws", "rk4", "dp54"]
    end function rk_method_names

    pure function rk3ws_tableau() result(tab)
        !! RK3 of Wicker and Skamarock: 3 stage(s), order 2 (3 for linear problems).
        !! Wicker & Skamarock (2002), Mon. Wea. Rev. 130, 2088
        type(rk_tableau_t) :: tab

        tab%name = "rk3ws"
        tab%title = "RK3 of Wicker and Skamarock"
        tab%reference = "Wicker & Skamarock (2002), Mon. Wea. Rev. 130, 2088"
        tab%stages = 3
        tab%order = 2
        tab%linear_order = 3
        tab%fsal = .false.
        tab%ssp_coef = 0.0_wp
        tab%f_evals = 3
        allocate(tab%A(3, 3), source=0.0_wp)
        tab%A(2, 1) = 1.0_wp/3.0_wp
        tab%A(3, 2) = 1.0_wp/2.0_wp
        tab%b = [0.0_wp, 0.0_wp, 1.0_wp]
        tab%c = [0.0_wp, 1.0_wp/3.0_wp, 1.0_wp/2.0_wp]
    end function rk3ws_tableau

    pure function rk4_tableau() result(tab)
        !! Classic RK4: 4 stage(s), order 4.
        !! Kutta (1901)
        type(rk_tableau_t) :: tab

        tab%name = "rk4"
        tab%title = "Classic RK4"
        tab%reference = "Kutta (1901)"
        tab%stages = 4
        tab%order = 4
        tab%linear_order = 4
        tab%fsal = .false.
        tab%ssp_coef = 0.0_wp
        tab%f_evals = 4
        allocate(tab%A(4, 4), source=0.0_wp)
        tab%A(2, 1) = 1.0_wp/2.0_wp
        tab%A(3, 2) = 1.0_wp/2.0_wp
        tab%A(4, 3) = 1.0_wp
        tab%b = [1.0_wp/6.0_wp, 1.0_wp/3.0_wp, 1.0_wp/3.0_wp, 1.0_wp/6.0_wp]
        tab%c = [0.0_wp, 1.0_wp/2.0_wp, 1.0_wp/2.0_wp, 1.0_wp]
    end function rk4_tableau

    pure function dp54_tableau() result(tab)
        !! Dormand-Prince 5(4): 7 stage(s), order 5, embedded order 4.
        !! Dormand & Prince (1980), J. Comput. Appl. Math. 6, 19; rklib rkdp54
        type(rk_tableau_t) :: tab

        tab%name = "dp54"
        tab%title = "Dormand-Prince 5(4)"
        tab%reference = "Dormand & Prince (1980), J. Comput. Appl. Math. 6, 19; rklib rkdp54"
        tab%stages = 7
        tab%order = 5
        tab%linear_order = 5
        tab%fsal = .true.
        tab%ssp_coef = 0.0_wp
        tab%f_evals = 6
        allocate(tab%A(7, 7), source=0.0_wp)
        tab%A(2, 1) = 1.0_wp/5.0_wp
        tab%A(3, 1) = 3.0_wp/40.0_wp
        tab%A(3, 2) = 9.0_wp/40.0_wp
        tab%A(4, 1) = 44.0_wp/45.0_wp
        tab%A(4, 2) = -56.0_wp/15.0_wp
        tab%A(4, 3) = 32.0_wp/9.0_wp
        tab%A(5, 1) = 19372.0_wp/6561.0_wp
        tab%A(5, 2) = -25360.0_wp/2187.0_wp
        tab%A(5, 3) = 64448.0_wp/6561.0_wp
        tab%A(5, 4) = -212.0_wp/729.0_wp
        tab%A(6, 1) = 9017.0_wp/3168.0_wp
        tab%A(6, 2) = -355.0_wp/33.0_wp
        tab%A(6, 3) = 46732.0_wp/5247.0_wp
        tab%A(6, 4) = 49.0_wp/176.0_wp
        tab%A(6, 5) = -5103.0_wp/18656.0_wp
        tab%A(7, 1) = 35.0_wp/384.0_wp
        tab%A(7, 3) = 500.0_wp/1113.0_wp
        tab%A(7, 4) = 125.0_wp/192.0_wp
        tab%A(7, 5) = -2187.0_wp/6784.0_wp
        tab%A(7, 6) = 11.0_wp/84.0_wp
        tab%b = [35.0_wp/384.0_wp, 0.0_wp, 500.0_wp/1113.0_wp, 125.0_wp/192.0_wp, -2187.0_wp/6784.0_wp, 11.0_wp/84.0_wp, 0.0_wp]
        tab%c = [0.0_wp, 1.0_wp/5.0_wp, 3.0_wp/10.0_wp, 4.0_wp/5.0_wp, 8.0_wp/9.0_wp, 1.0_wp, 1.0_wp]
        tab%b_hat = [5179.0_wp/57600.0_wp, 0.0_wp, 7571.0_wp/16695.0_wp, 393.0_wp/640.0_wp, -92097.0_wp/339200.0_wp,&
            & 187.0_wp/2100.0_wp, 1.0_wp/40.0_wp]
        tab%order_hat = 4
    end function dp54_tableau

    !> Raise an error in `status` unless `tableau` is usable by the explicit
    !> engine: at least one stage, A, b and c allocated with the stage count's
    !> sizes, and A strictly lower triangular.
    subroutine check_tableau(tableau, status)
        type(rk_tableau_t), intent(in)    :: tableau
        type(status_t),     intent(inout) :: status
        !
        integer :: s, i

        s = tableau%stages
        if (s < 1) then
            call status%raise("invalid tableau: no stages (an unknown name given to rk_tableau?)")
        else if (.not. (allocated(tableau%A) .and. allocated(tableau%b) .and. allocated(tableau%c))) then
            call status%raise("invalid tableau: A, b or c is not allocated")
        else if (any(shape(tableau%A) /= [s, s]) .or. size(tableau%b) /= s .or. size(tableau%c) /= s) then
            call status%raise("invalid tableau: A, b and c do not match the number of stages")
        else
            do i = 1, s
                if (any(tableau%A(i, i:) /= 0.0_wp)) then
                    call status%raise("invalid tableau: A is not strictly lower triangular (not an explicit method)")
                    return
                end if
            end do
        end if
    end subroutine check_tableau

    !> A stage nothing uses: neither b nor a later row of A reads its tendency.
    pure logical function stage_unused(tableau, i)
        type(rk_tableau_t), intent(in) :: tableau
        integer,            intent(in) :: i
        stage_unused = tableau%b(i) == 0.0_wp .and. all(tableau%A(i+1:, i) == 0.0_wp)
    end function stage_unused

end module stdlib_integrate_rk_tableaus
