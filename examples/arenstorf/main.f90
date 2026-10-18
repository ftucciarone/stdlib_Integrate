program main
    !! Arenstorf orbit (SciML prob_ode_threebody; Hairer, Norsett & Wanner,
    !! Solving ODEs I, section II.0): a light body in the plane of the Earth and
    !! the Moon, in the frame rotating with them (Earth at (-mu, 0), Moon at
    !! (1 - mu, 0)),
    !!   y1'' = y1 + 2 y2' - mu' (y1 + mu)/D1 - mu (y1 - mu')/D2,
    !!   y2'' = y2 - 2 y1' - mu' y2/D1 - mu y2/D2,
    !!   D1 = ((y1 + mu)^2 + y2^2)^(3/2),  D2 = ((y1 - mu')^2 + y2^2)^(3/2),  mu' = 1 - mu,
    !! with initial values on a periodic orbit of period T. The solution is not
    !! known in closed form, but it returns to its start after T: the closing
    !! error |u(T) - u(0)| is the global error. Tests:
    !!   - convergence: the closing error as the number of steps doubles;
    !!   - cost: the fewest steps (powers of 2 times n_start) with a closing
    !!     error below `tolerance`.
    !! Uses the plain-array schemes, state (y1, y2, y1', y2'). Writes
    !! examples/arenstorf/data/: parameters.tex, results.tex, reference.dat
    !! (dp54 with n_reference steps), coarse_<method>.dat and fine_<method>.dat
    !! (orbits with n_coarse and n_fine steps, n_points points each),
    !! history_<method>.dat (t, |y - y_reference| at n_fine steps),
    !! error_<method>.dat (dt, closing error), read by summary.tex. Run from the
    !! project root.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: rk_scheme, stage_t, rk_array_scheme_t, rk_method_names, rk_tableau, rk_tableau_t
    implicit none

    character(*), parameter :: out = "examples/arenstorf/data/"
    real(wp),     parameter :: mu = 0.012277471_wp, mu1 = 1 - mu
    real(wp),     parameter :: period = 17.0652165601579625588917206249_wp
    real(wp),     parameter :: u0(4) = [0.994_wp, 0.0_wp, 0.0_wp, -2.00158510637908252240537862224_wp]
    integer,      parameter :: n_coarse = 4000                 !! steps of the coarse solution figure
    integer,      parameter :: n_fine = 64000                  !! steps of the fine figure and the error history
    integer,      parameter :: n_reference = 10 * n_fine       !! steps of the reference (dp54)
    integer,      parameter :: n_points = 2000                 !! points written per orbit (TeX memory)
    integer,      parameter :: n_start = 1000, n_levels = 13   !! convergence: 1000, 2000, ..., 4096000 steps
    real(wp),     parameter :: tolerance = 1.0e-3_wp           !! closing error asked for in the cost test
    real(wp),     parameter :: round_off = 1.0e-8_wp           !! round-off accumulates to ~1e-9 over millions of steps

    character(len=8), allocatable :: names(:)
    type(rk_tableau_t) :: tableau
    real(wp), allocatable :: reference(:,:)
    real(wp) :: dt(n_levels), err(n_levels), observed
    integer  :: m, l, n_needed, table, dat

    names = rk_method_names()
    call execute_command_line("mkdir -p " // out)
    call write_parameters()

    ! Reference: dp54 with n_reference steps, kept at the n_fine step times.
    allocate(reference(4, 0:n_fine))
    call run("dp54", n_reference, err(1), keep=reference, every=n_reference / n_fine, orbit="reference.dat")
    print '(a, es10.3)', "reference closing error ", err(1)
    do m = 1, size(names)
        call run(trim(names(m)), n_coarse, err(1), orbit="coarse_" // trim(names(m)) // ".dat")
        call run(trim(names(m)), n_fine, err(1), orbit="fine_" // trim(names(m)) // ".dat", &
                 history="history_" // trim(names(m)) // ".dat")
    end do

    open(newunit=table, file=out // "results.tex", status="replace", action="write")
    write(table, '(a)') "% Written by examples/arenstorf/main.f90: do not edit."
    write(table, '(a)') "\begin{tabular}{lcccccc}"
    write(table, '(a)') "\toprule"
    write(table, '(a)') "name & order & observed & closing, " // itoa(n_fine) // " steps & " // &
                        "steps to " // sci(tolerance) // " & $f$ evaluations \\"
    write(table, '(a)') "\midrule"
    print '(a20, a10, a6, a9, a14, a10)', "method", "name", "order", "obs.", "closing", "steps"
    do m = 1, size(names)
        tableau = rk_tableau(trim(names(m)))
        open(newunit=dat, file=out // "error_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat, '(a)') "dt error"
        n_needed = 0
        do l = 1, n_levels
            dt(l) = period / (n_start * 2**(l - 1))
            call run(trim(names(m)), n_start * 2**(l - 1), err(l))
            if (err(l) > 1.0e-15_wp .and. err(l) < huge(1.0_wp)) write(dat, '(2es14.6)') dt(l), err(l)
            if (n_needed == 0 .and. err(l) < tolerance) n_needed = n_start * 2**(l - 1)
        end do
        close(dat)
        observed = slope(dt, err)
        call run(trim(names(m)), n_fine, err(1))
        print '(a20, a10, i6, f9.3, es14.3, i10)', tableau%title, trim(names(m)), tableau%order, observed, err(1), n_needed
        write(table, '("\texttt{", a, "} & ", i0, " & ", a, " & ", a, " & ", a, " & ", a, " \\")') &
            trim(names(m)), tableau%order, order_text(observed), sci(err(1)), steps_text(n_needed), &
            evaluations_text(n_needed, trim(names(m)))
    end do
    write(table, '(a)') "\bottomrule"
    write(table, '(a)') "\end{tabular}"
    close(table)

contains

    subroutine arenstorf(t, u, dudt, stage)
        real(wp), intent(in)  :: t
        real(wp), intent(in)  :: u(:)
        real(wp), intent(out) :: dudt(:)
        type(stage_t), intent(in) :: stage
        real(wp) :: d1, d2

        d1 = ((u(1) + mu)**2 + u(2)**2)**1.5_wp
        d2 = ((u(1) - mu1)**2 + u(2)**2)**1.5_wp
        dudt(1:2) = u(3:4)
        dudt(3) = u(1) + 2 * u(4) - mu1 * (u(1) + mu) / d1 - mu * (u(1) - mu1) / d2
        dudt(4) = u(2) - 2 * u(3) - mu1 * u(2) / d1 - mu * u(2) / d2
    end subroutine arenstorf

    !> One period in n steps: `err` = closing error max_i |u_i(T) - u_i(0)|
    !> (huge if the run blew up). Optional: keep every `every`-th state in
    !> `keep`, write the orbit (y1, y2), write the distance to the reference.
    subroutine run(name, n, err, keep, every, orbit, history)
        character(*), intent(in)  :: name
        integer,      intent(in)  :: n
        real(wp),     intent(out) :: err
        real(wp),     intent(inout), optional :: keep(:, 0:)
        integer,      intent(in),    optional :: every
        character(*), intent(in),    optional :: orbit, history

        type(rk_array_scheme_t) :: scheme
        real(wp) :: u(4), time
        integer  :: i, stride, u_orbit, u_history, orbit_stride

        stride = 1
        if (present(every)) stride = every
        orbit_stride = max(1, n / n_points)
        call rk_scheme(scheme, name, 4, period / n)
        u = u0
        time = 0.0_wp
        err = huge(1.0_wp)
        if (present(keep)) keep(:, 0) = u
        if (present(orbit)) then
            open(newunit=u_orbit, file=out // orbit, status="replace", action="write")
            write(u_orbit, '(a)') "x y"
            write(u_orbit, '(2es16.8)') u(1:2)
        end if
        if (present(history)) then
            open(newunit=u_history, file=out // history, status="replace", action="write")
            write(u_history, '(a)') "t error"
        end if
        do i = 1, n
            call scheme%advance(arenstorf, time, u)
            if (any(u /= u) .or. maxval(abs(u)) > 1.0e3_wp) exit           ! blew up
            if (present(keep) .and. mod(i, stride) == 0) keep(:, i / stride) = u
            if (present(orbit) .and. mod(i, orbit_stride) == 0) write(u_orbit, '(2es16.8)') u(1:2)
            if (present(history) .and. mod(i, orbit_stride) == 0) write(u_history, '(2es14.6)') time, &
                max(norm2(u(1:2) - reference(1:2, i)), 1.0e-16_wp)
        end do
        if (i > n) err = maxval(abs(u - u0))
        if (present(orbit)) close(u_orbit)
        if (present(history)) close(u_history)
    end subroutine run

    real(wp) function slope(dt, err) result(p)
        real(wp), intent(in) :: dt(:), err(:)
        integer :: l

        p = 0.0_wp
        do l = size(dt), 2, -1
            if (err(l) > round_off .and. err(l - 1) > round_off .and. err(l - 1) < 1.0e-1_wp) then
                p = log(err(l - 1) / err(l)) / log(dt(l - 1) / dt(l))
                return
            end if
        end do
    end function slope

    function order_text(p) result(s)
        real(wp), intent(in) :: p
        character(:), allocatable :: s
        character(len=8) :: buffer

        if (p == 0.0_wp) then
            s = "---"
        else
            write(buffer, '(f5.2)') p
            s = trim(adjustl(buffer))
        end if
    end function order_text

    function steps_text(n) result(s)
        integer, intent(in) :: n
        character(:), allocatable :: s

        if (n == 0) then
            s = "$> " // itoa(n_start * 2**(n_levels - 1)) // "$"
        else
            s = itoa(n)
        end if
    end function steps_text

    !> Evaluations of f for n steps: steps times the stages the step
    !> evaluates (dp54's seventh stage, used only by its error estimate, is skipped).
    function evaluations_text(n, name) result(s)
        integer,      intent(in) :: n
        character(*), intent(in) :: name
        character(:), allocatable :: s
        integer :: stages

        select case (name)
        case ("rk3ws");            stages = 3
        case ("rk4");              stages = 4
        case default;              stages = 6
        end select
        if (n == 0) then
            s = "---"
        else
            s = itoa(n * stages)
        end if
    end function evaluations_text

    subroutine write_parameters()
        integer :: u, k
        character(:), allocatable :: list

        list = ""
        do k = 1, size(names)
            if (k > 1) list = list // ", "
            list = list // "\texttt{" // trim(names(k)) // "}"
        end do
        open(newunit=u, file=out // "parameters.tex", status="replace", action="write")
        write(u, '(a)') "% Written by examples/arenstorf/main.f90: do not edit."
        write(u, '(a)') "\begin{tabular}{ll}"
        write(u, '(a)') "\toprule"
        write(u, '(a)') "mass ratio & $\mu = 0.012277471$ (Moon / (Earth + Moon)) \\"
        write(u, '(a)') "initial position & $y = (0.994, 0)$ \\"
        write(u, '(a)') "initial velocity & $\dot y = (0, -2.00158510637908252240537862224)$ \\"
        write(u, '(a)') "period & $T = 17.0652165601579625588917206249$ \\"
        write(u, '(a)') "reference & \texttt{dp54}, " // itoa(n_reference) // " steps \\"
        write(u, '(a)') "solution figures & " // itoa(n_coarse) // " and " // itoa(n_fine) // &
                        " steps per period (" // itoa(n_points) // " points drawn) \\"
        write(u, '(a)') "convergence & closing error $\max_i |u_i(T) - u_i(0)|$, " // itoa(n_start) // ", " // &
                        itoa(2 * n_start) // ", \dots, " // itoa(n_start * 2**(n_levels - 1)) // " steps \\"
        write(u, '(a)') "observed order & slope of the two finest levels with errors in $(10^{" // &
                        itoa(nint(log10(round_off))) // "}, 10^{-1})$ \\"
        write(u, '(a)') "cost & fewest such step count closing to " // sci(tolerance) // &
                        "; evaluations of $f$ \\"
        write(u, '(a)') "methods & " // list // " \\"
        write(u, '(a)') "interface & plain arrays, \texttt{rk\_array\_scheme\_t}, state $(y_1, y_2, \dot y_1, \dot y_2)$ \\"
        write(u, '(a)') "\bottomrule"
        write(u, '(a)') "\end{tabular}"
        close(u)
    end subroutine write_parameters

    function sci(x) result(s)
        real(wp), intent(in) :: x
        character(:), allocatable :: s
        character(len=32) :: mantissa
        integer :: ex

        if (x >= huge(1.0_wp)) then
            s = "---"
            return
        end if
        if (abs(x) >= 0.1_wp .and. abs(x) < 1000.0_wp) then
            write(mantissa, '(f8.2)') abs(x)
            s = "$" // trim(adjustl(mantissa)) // "$"
            return
        end if
        ex = floor(log10(abs(x)))
        write(mantissa, '(f4.2)') abs(x) / 10.0_wp**ex
        s = "$" // trim(mantissa) // "\times10^{" // itoa(ex) // "}$"
    end function sci

    function itoa(i) result(s)
        integer, intent(in) :: i
        character(:), allocatable :: s
        character(len=16) :: buffer

        write(buffer, '(i0)') i
        s = trim(buffer)
    end function itoa

end program main
