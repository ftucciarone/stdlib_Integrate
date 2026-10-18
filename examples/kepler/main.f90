program main
    !! Kepler problem: a body around a fixed mass (G M = 1),
    !!   q'' = -q / |q|^3,   q(0) = (1 - e, 0),   q'(0) = (0, sqrt((1 + e)/(1 - e))),
    !! an ellipse of eccentricity e, semi-major axis 1 and period 2 pi, starting
    !! at the pericentre. Exact solution from Kepler's equation E - e sin E = t:
    !!   q = (cos E - e, sqrt(1 - e^2) sin E),
    !!   q' = (-sin E, sqrt(1 - e^2) cos E) / (1 - e cos E).
    !! Energy H = |q'|^2/2 - 1/|q| = -1/2 and angular momentum L = q1 q2' - q2 q1'
    !! = sqrt(1 - e^2) are constant. Tests:
    !!   - convergence: the error at t = n_conv periods as the steps per period
    !!     double, at the method's order;
    !!   - long time: H and L over n_long periods at n_long_steps steps per
    !!     period: explicit Runge-Kutta methods keep neither.
    !! Uses the plain-array schemes, state (q1, q2, q1', q2'). Writes
    !! examples/kepler/data/: parameters.tex, results.tex, ellipse.dat (exact
    !! orbit), orbit_<method>.dat (q at n_plot_steps steps per period over
    !! n_plot periods), history_<method>.dat (t, |q - q_exact|), error_<method>.dat
    !! (dt, error), invariants_<method>.dat (t, running max |H - H0|, |L - L0|),
    !! read by summary.tex. Run from the project root.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: rk_scheme, stage_t, rk_array_scheme_t, rk_method_names, rk_tableau, rk_tableau_t
    implicit none

    character(*), parameter :: out = "examples/kepler/data/"
    real(wp),     parameter :: pi = acos(-1.0_wp)
    real(wp),     parameter :: e = 0.6_wp                       !! eccentricity
    integer,      parameter :: n_plot = 10, n_plot_steps = 200  !! solution figures: periods, steps per period
    integer,      parameter :: n_conv = 3                       !! convergence: periods
    integer,      parameter :: n_conv_steps = 100, n_levels = 8 !! convergence: 100, 200, ... steps per period
    integer,      parameter :: n_long = 100, n_long_steps = 400 !! long time: periods, steps per period
    real(wp),     parameter :: round_off = 1.0e-10_wp       !! round-off accumulates over the many steps

    character(len=8), allocatable :: names(:)
    type(rk_tableau_t) :: tableau
    real(wp) :: dt(n_levels), err(n_levels), observed, h_drift, l_drift, h0, l0, last
    integer  :: m, l, i, table, dat

    names = rk_method_names()
    h0 = energy(exact(0.0_wp))
    l0 = momentum(exact(0.0_wp))
    call execute_command_line("mkdir -p " // out)
    call write_parameters()

    open(newunit=dat, file=out // "ellipse.dat", status="replace", action="write")
    write(dat, '(a)') "x y"
    do i = 0, 400
        write(dat, '(2es16.8)') cos(2 * pi * i / 400) - e, sqrt(1 - e**2) * sin(2 * pi * i / 400)
    end do
    close(dat)
    do m = 1, size(names)
        call run(trim(names(m)), n_plot_steps, n_plot, err(1), h_drift, l_drift, last, &
                 "orbit_" // trim(names(m)) // ".dat", "history_" // trim(names(m)) // ".dat")
    end do

    open(newunit=table, file=out // "results.tex", status="replace", action="write")
    write(table, '(a)') "% Written by examples/kepler/main.f90: do not edit."
    write(table, '(a)') "\begin{tabular}{llcccc}"
    write(table, '(a)') "\toprule"
    write(table, '(a)') "method & name & order & observed & $\max|H - H_0|$ & $\max|L - L_0|$ \\"
    write(table, '(a)') "\midrule"
    print '(a20, a10, a6, a9, 2a13)', "method", "name", "order", "obs.", "max|H-H0|", "max|L-L0|"
    do m = 1, size(names)
        tableau = rk_tableau(trim(names(m)))
        open(newunit=dat, file=out // "error_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat, '(a)') "dt error"
        do l = 1, n_levels
            dt(l) = 2 * pi / (n_conv_steps * 2**(l - 1))
            call run(trim(names(m)), n_conv_steps * 2**(l - 1), n_conv, err(l), h_drift, l_drift, last)
            if (err(l) > 1.0e-15_wp .and. err(l) < huge(1.0_wp)) write(dat, '(2es14.6)') dt(l), err(l)
        end do
        close(dat)
        observed = slope(dt, err)
        call run(trim(names(m)), n_long_steps, n_long, err(1), h_drift, l_drift, last, &
                 invariants="invariants_" // trim(names(m)) // ".dat")
        print '(a20, a10, i6, f9.3, 2es13.3)', tableau%title, trim(names(m)), tableau%order, observed, h_drift, l_drift
        write(table, '(a, " & \texttt{", a, "} & ", i0, " & ", a, " & ", a, " & ", a, " \\")') &
            tableau%title, trim(names(m)), tableau%order, order_text(observed), sci(h_drift), sci(l_drift)
    end do
    write(table, '(a)') "\bottomrule"
    write(table, '(a)') "\end{tabular}"
    close(table)

contains

    subroutine kepler(t, u, dudt, stage)
        real(wp), intent(in)  :: t
        real(wp), intent(in)  :: u(:)
        real(wp), intent(out) :: dudt(:)
        type(stage_t), intent(in) :: stage
        real(wp) :: r3

        r3 = norm2(u(1:2))**3
        dudt(1:2) = u(3:4)
        dudt(3:4) = -u(1:2) / r3
    end subroutine kepler

    !> Exact state (q1, q2, q1', q2') at time t, from Kepler's equation
    !> E - e sin E = t solved by Newton's method.
    function exact(t) result(u)
        real(wp), intent(in) :: t
        real(wp) :: u(4), mean, ecc_anomaly, delta
        integer  :: it

        mean = modulo(t, 2 * pi)
        ecc_anomaly = mean + e * sin(mean)
        do it = 1, 50
            delta = (ecc_anomaly - e * sin(ecc_anomaly) - mean) / (1 - e * cos(ecc_anomaly))
            ecc_anomaly = ecc_anomaly - delta
            if (abs(delta) < 1.0e-16_wp) exit
        end do
        u(1) = cos(ecc_anomaly) - e
        u(2) = sqrt(1 - e**2) * sin(ecc_anomaly)
        u(3) = -sin(ecc_anomaly) / (1 - e * cos(ecc_anomaly))
        u(4) = sqrt(1 - e**2) * cos(ecc_anomaly) / (1 - e * cos(ecc_anomaly))
    end function exact

    pure real(wp) function energy(u)
        real(wp), intent(in) :: u(4)

        energy = (u(3)**2 + u(4)**2) / 2 - 1 / norm2(u(1:2))
    end function energy

    pure real(wp) function momentum(u)
        real(wp), intent(in) :: u(4)

        momentum = u(1) * u(4) - u(2) * u(3)
    end function momentum

    !> n_periods periods at steps_per_period: `err` = max-norm error of the
    !> state at the end, `h_drift`, `l_drift` = max |H - H0|, |L - L0| over the
    !> run, `last` = the time reached (earlier if the orbit escaped or
    !> collapsed). Optional files: the orbit (q1, q2), the position error
    !> history, the invariants' running maxima.
    subroutine run(name, steps_per_period, n_periods, err, h_drift, l_drift, last, orbit, history, invariants)
        character(*), intent(in)  :: name
        integer,      intent(in)  :: steps_per_period, n_periods
        real(wp),     intent(out) :: err, h_drift, l_drift, last
        character(*), intent(in), optional :: orbit, history, invariants

        type(rk_array_scheme_t) :: scheme
        real(wp) :: u(4), time, step
        integer  :: i, n, u_orbit, u_history, u_inv, next

        step = 2 * pi / steps_per_period
        n = steps_per_period * n_periods
        call rk_scheme(scheme, name, 4, step)
        u = exact(0.0_wp)
        time = 0.0_wp
        h_drift = 0.0_wp
        l_drift = 0.0_wp
        err = huge(1.0_wp)
        last = 0.0_wp
        next = 1
        if (present(orbit)) call open_file(orbit, "x y", u_orbit)
        if (present(history)) call open_file(history, "t error", u_history)
        if (present(invariants)) call open_file(invariants, "t energy momentum", u_inv)
        if (present(orbit)) write(u_orbit, '(2es16.8)') u(1:2)
        do i = 1, n
            call scheme%advance(kepler, time, u)
            if (any(u /= u) .or. norm2(u(1:2)) > 100 .or. norm2(u(1:2)) < 1.0e-6_wp) exit   ! escaped or collapsed
            last = time
            h_drift = max(h_drift, abs(energy(u) - h0))
            l_drift = max(l_drift, abs(momentum(u) - l0))
            if (present(orbit)) write(u_orbit, '(2es16.8)') u(1:2)
            if (present(history)) write(u_history, '(2es14.6)') time, max(norm2(u(1:2) - exact(time)), 1.0e-16_wp)
            if (present(invariants) .and. (i == next .or. i == n)) then
                write(u_inv, '(3es14.6)') time, max(h_drift, 1.0e-16_wp), max(l_drift, 1.0e-16_wp)
                next = max(next + 1, nint(next * 1.05_wp))
            end if
        end do
        if (i > n) err = maxval(abs(u - exact(time)))
        if (i <= n) then                                     ! did not finish: no meaningful drift either
            h_drift = huge(1.0_wp)
            l_drift = huge(1.0_wp)
        end if
        if (present(orbit)) close(u_orbit)
        if (present(history)) close(u_history)
        if (present(invariants)) close(u_inv)
    end subroutine run

    subroutine open_file(file, header, unit)
        character(*), intent(in)  :: file, header
        integer,      intent(out) :: unit

        open(newunit=unit, file=out // file, status="replace", action="write")
        write(unit, '(a)') header
    end subroutine open_file

    !> Slope of log(err) against log(dt) between the two finest levels whose
    !> errors are in (round_off, 1).
    real(wp) function slope(dt, err) result(p)
        real(wp), intent(in) :: dt(:), err(:)
        integer :: l

        p = 0.0_wp
        do l = size(dt), 2, -1
            if (err(l) > round_off .and. err(l - 1) > round_off .and. err(l - 1) < 1.0_wp) then
                p = log(err(l - 1) / err(l)) / log(dt(l - 1) / dt(l))
                return
            end if
        end do
    end function slope

    !> The observed order, or a dash when no two levels had errors in (round_off, 1).
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

    subroutine write_parameters()
        integer :: u, k
        character(:), allocatable :: list

        list = ""
        do k = 1, size(names)
            if (k > 1) list = list // ", "
            list = list // "\texttt{" // trim(names(k)) // "}"
        end do
        open(newunit=u, file=out // "parameters.tex", status="replace", action="write")
        write(u, '(a)') "% Written by examples/kepler/main.f90: do not edit."
        write(u, '(a)') "\begin{tabular}{ll}"
        write(u, '(a)') "\toprule"
        write(u, '(a)') "eccentricity & $e = " // fmt(e) // "$ (semi-major axis 1, period $2\pi$, $GM = 1$) \\"
        write(u, '(a)') "initial condition & $q = (" // fmt(1 - e) // ", 0)$, $\dot q = (0, " // &
                        fmt(sqrt((1 + e) / (1 - e))) // ")$: the pericentre \\"
        write(u, '(a)') "invariants & $H_0 = " // fmt(h0) // "$, $L_0 = \sqrt{1 - e^2} = " // fmt(l0) // "$ \\"
        write(u, '(a)') "solution figures & " // itoa(n_plot_steps) // " steps per period, " // itoa(n_plot) // &
                        " periods \\"
        write(u, '(a)') "convergence & error $\max_i |u_i - u_i^{\mathrm{exact}}|$ after " // itoa(n_conv) // &
                        " periods, " // itoa(n_conv_steps) // ", " // itoa(2 * n_conv_steps) // ", \dots, " // &
                        itoa(n_conv_steps * 2**(n_levels - 1)) // " steps per period \\"
        write(u, '(a)') "observed order & slope of the two finest levels with errors in $(10^{" // &
                        itoa(nint(log10(round_off))) // "}, 1)$ \\"
        write(u, '(a)') "long time & $\max|H - H_0|$, $\max|L - L_0|$ over " // itoa(n_long) // " periods, " // &
                        itoa(n_long_steps) // " steps per period \\"
        write(u, '(a)') "exact solution & Kepler's equation $E - e\sin E = t$, Newton's method \\"
        write(u, '(a)') "methods & " // list // " \\"
        write(u, '(a)') "interface & plain arrays, \texttt{rk\_array\_scheme\_t}, state $(q_1, q_2, \dot q_1, \dot q_2)$ \\"
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
        if (x == 0.0_wp) then
            s = "$0$"
            return
        end if
        ex = floor(log10(abs(x)))
        write(mantissa, '(f4.2)') abs(x) / 10.0_wp**ex
        s = "$" // trim(mantissa) // "\times10^{" // itoa(ex) // "}$"
    end function sci

    function fmt(x) result(s)
        real(wp), intent(in) :: x
        character(:), allocatable :: s
        character(len=24) :: buffer

        write(buffer, '(f14.4)') x
        s = trim(adjustl(buffer))
        do while (s(len(s):len(s)) == "0")
            s = s(1:len(s) - 1)
        end do
        if (s(len(s):len(s)) == ".") s = s(1:len(s) - 1)
    end function fmt

    function itoa(i) result(s)
        integer, intent(in) :: i
        character(:), allocatable :: s
        character(len=16) :: buffer

        write(buffer, '(i0)') i
        s = trim(buffer)
    end function itoa

end program main
