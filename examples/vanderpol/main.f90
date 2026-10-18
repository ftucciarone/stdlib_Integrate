program main
    !! Van der Pol oscillator (SciML prob_ode_vanderpol, mu = 1):
    !!   x' = y,   y' = mu ((1 - x^2) y - x),   (x, y)(0) = (sqrt(3), 0).
    !! With tau = sqrt(mu) t it is the classic x'' - eps (1 - x^2) x' + x = 0,
    !! eps = sqrt(mu): every solution but the origin tends to a limit cycle of
    !! amplitude about 2. No closed form: the reference is dp54 at a tiny step.
    !! Tests:
    !!   - convergence: the error at t = t_conv as the step halves;
    !!   - period of the limit cycle, from the times x crosses 0 upwards: each
    !!     crossing is found by solving for the step size that lands on x = 0
    !!     with the method itself (secant), so it is as accurate as the method;
    !!   - stiffness: on the slow branch near x = 2 the Jacobian has the
    !!     eigenvalue mu (1 - x^2) ~ -3 mu, so the largest stable step should
    !!     approach h*/(3 mu) for large mu, h* the method's real-axis limit.
    !! Uses the plain-array schemes. Writes examples/vanderpol/data/:
    !! parameters.tex, results.tex, field.dat (x, y, unit direction, speed),
    !! reference.dat (t, x, y), start_<k>.dat (trajectories from other initial
    !! points), cycle_<k>.dat (limit cycles for several mu, (t, x, y/sqrt(mu))),
    !! solution_<method>.dat (t, x at dt_plot), error_<method>.dat (dt, error),
    !! period_<method>.dat (dt, period error), stiff_<method>.dat (mu, mu dt_max),
    !! read by summary.tex. Run from the project root.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: rk_scheme, stage_t, rk_array_scheme_t, rk_method_names, rk_tableau, rk_tableau_t
    use vanderpol_stability, only: stability_r
    implicit none

    character(*), parameter :: out = "examples/vanderpol/data/"
    real(wp),     parameter :: mu1 = 1.0_wp                     !! the problem's mu
    real(wp),     parameter :: t_conv = 10.0_wp                 !! convergence: error at this time
    real(wp),     parameter :: dt_coarse = 0.1_wp
    integer,      parameter :: n_levels = 9                     !! dt = 0.1, 0.05, ..., 0.1/256
    real(wp),     parameter :: dt_reference = 1.0e-4_wp
    real(wp),     parameter :: dt_plot = 0.25_wp, t_plot = 20.0_wp
    integer,      parameter :: first_cycle = 10, last_cycle = 15 !! period from these crossings (after the transient)
    real(wp),     parameter :: mus(7) = [1.0_wp, 3.0_wp, 10.0_wp, 30.0_wp, 100.0_wp, 300.0_wp, 1000.0_wp]
    real(wp),     parameter :: t_stiff = 10.0_wp                !! stiffness runs: interval
    real(wp),     parameter :: floor = 1.0e-11_wp

    character(len=8), allocatable :: names(:)
    type(rk_tableau_t) :: tableau
    real(wp) :: mu, u_ref(2), u(2), period_ref, dt(n_levels), err(n_levels), perr(n_levels)
    real(wp) :: observed, observed_p, stiff_last, along_cycle
    real(wp), allocatable :: h_star(:)   ! per method of the table
    complex(wp), allocatable :: eigen(:)
    integer  :: m, l, k, table, dat

    names = rk_method_names()
    allocate(h_star(size(names)))
    do m = 1, size(names)
        h_star(m) = real_axis_limit(rk_tableau(trim(names(m))))
    end do
    call execute_command_line("mkdir -p " // out)
    mu = mu1

    ! Reference: state at t_conv and the period, with dp54 at dt_reference.
    u_ref = integrate("dp54", dt_reference, t_conv)
    period_ref = period("dp54", dt_reference)
    print '(a, f18.14)', "reference period ", period_ref
    call write_parameters()
    call write_figures()

    open(newunit=table, file=out // "results.tex", status="replace", action="write")
    write(table, '(a)') "% Written by examples/vanderpol/main.f90: do not edit."
    write(table, '(a)') "\begin{tabular}{lccccccl}"
    write(table, '(a)') "\toprule"
    write(table, '(a)') "& & \multicolumn{2}{c}{observed order} & \multicolumn{4}{c}{stiffness, " // &
                        "$\mu\,\Delta t_{\max}$ at $\mu = " // itoa(nint(mus(size(mus)))) // "$} \\"
    write(table, '(a)') "\cmidrule(lr){3-4}\cmidrule(lr){5-8}"
    write(table, '(a)') "name & order & solution & period & measured & $h^*/3$ & along the cycle & " // &
                        "first blows up \\"
    write(table, '(a)') "\midrule"
    eigen = cycle_eigenvalues(mus(size(mus)))
    print '(a10, a6, 2a9, 3a10)', "name", "order", "sol.", "period", "mu dt", "cycle", "h*/3"
    do m = 1, size(names)
        tableau = rk_tableau(trim(names(m)))
        mu = mu1
        open(newunit=dat, file=out // "error_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat, '(a)') "dt error"
        do l = 1, n_levels
            dt(l) = dt_coarse / 2**(l - 1)
            u = integrate(trim(names(m)), dt(l), t_conv)
            err(l) = maxval(abs(u - u_ref))
            if (err(l) > 1.0e-15_wp .and. err(l) < 1.0e3_wp) write(dat, '(2es14.6)') dt(l), err(l)
        end do
        close(dat)
        open(newunit=dat, file=out // "period_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat, '(a)') "dt error"
        do l = 1, n_levels
            perr(l) = abs(period(trim(names(m)), dt(l)) - period_ref)
            if (perr(l) > 1.0e-15_wp .and. perr(l) < 1.0e3_wp) write(dat, '(2es14.6)') dt(l), perr(l)
        end do
        close(dat)
        observed = slope(dt, err)
        observed_p = slope(dt, perr)
        open(newunit=dat, file=out // "stiff_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat, '(a)') "mu scaled"
        do k = 1, size(mus)
            mu = mus(k)
            stiff_last = mu * stable_step(trim(names(m)))
            write(dat, '(2es14.6)') mu, stiff_last
        end do
        close(dat)
        along_cycle = mus(size(mus)) * limit_along_cycle(tableau, eigen)
        print '(a10, i6, 2f9.3, 3f10.4)', trim(names(m)), tableau%order, observed, observed_p, stiff_last, &
            along_cycle, h_star(m) / 3
        write(table, '("\texttt{", a, "} & ", i0, " & ", a, " & ", a, " & ", f6.4, " & ", f6.4, " & ", f6.4, " & ", a, " \\")') &
            trim(names(m)), tableau%order, order_text(observed), order_text(observed_p), stiff_last, h_star(m) / 3, &
            along_cycle, failure_place(trim(names(m)), 1.05_wp * stiff_last / mus(size(mus)))
    end do
    write(table, '(a)') "\bottomrule"
    write(table, '(a)') "\end{tabular}"
    close(table)

contains

    subroutine vanderpol(t, u, dudt, stage)
        real(wp), intent(in)  :: t
        real(wp), intent(in)  :: u(:)
        real(wp), intent(out) :: dudt(:)
        type(stage_t), intent(in) :: stage

        dudt(1) = u(2)
        dudt(2) = mu * ((1 - u(1)**2) * u(2) - u(1))
    end subroutine vanderpol

    !> The state at t_end from (sqrt(3), 0), steps of size dt (t_end/dt integer).
    function integrate(name, dt, t_end) result(u)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: dt, t_end
        real(wp) :: u(2), time
        type(rk_array_scheme_t) :: scheme
        integer :: i

        call rk_scheme(scheme, name, 2, dt)
        u = [sqrt(3.0_wp), 0.0_wp]
        time = 0.0_wp
        do i = 1, nint(t_end / dt)
            call scheme%advance(vanderpol, time, u)
        end do
    end function integrate

    !> Period of the limit cycle: mean spacing of the upward zero crossings of x
    !> number first_cycle to last_cycle. A crossing inside a step is located by
    !> solving x(step of size s from u_n) = 0 for s with the secant method, so
    !> it has the accuracy of the method.
    real(wp) function period(name, dt) result(p)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: dt

        type(rk_array_scheme_t) :: scheme
        real(wp) :: u(2), u_old(2), time, crossings(last_cycle)
        integer  :: n_cross, n_steps

        call rk_scheme(scheme, name, 2, dt)
        u = [sqrt(3.0_wp), 0.0_wp]
        time = 0.0_wp
        n_cross = 0
        n_steps = 0
        do while (n_cross < last_cycle)
            u_old = u
            call scheme%advance(vanderpol, time, u)
            n_steps = n_steps + 1
            if (u_old(1) < 0 .and. u(1) >= 0) then
                n_cross = n_cross + 1
                ! (n_steps - 1) dt, not the accumulated `time`: summing dt
                ! 10^5 times would add ~1e-10 of round-off to the crossing time.
                crossings(n_cross) = (n_steps - 1) * dt + crossing_step(name, u_old, dt)
            end if
            if (time > 1.0e4_wp .or. any(u /= u)) then
                p = huge(1.0_wp)
                return
            end if
        end do
        p = (crossings(last_cycle) - crossings(first_cycle)) / (last_cycle - first_cycle)
    end function period

    !> The step s in (0, dt] after which one step of the method from u0 has x = 0.
    real(wp) function crossing_step(name, u0, dt) result(s)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: u0(2), dt
        real(wp) :: s0, s1, x0, x1
        integer  :: it

        s0 = 0.0_wp
        x0 = u0(1)
        s1 = dt
        x1 = x_after(name, u0, s1)
        do it = 1, 50
            s = s1 - x1 * (s1 - s0) / (x1 - x0)
            s0 = s1
            x0 = x1
            s1 = s
            x1 = x_after(name, u0, s1)
            if (abs(x1) < 1.0e-15_wp .or. abs(s1 - s0) < 1.0e-16_wp) exit
        end do
        s = s1
    end function crossing_step

    !> x after one step of size `step` of the method from u0.
    real(wp) function x_after(name, u0, step)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: u0(2), step
        type(rk_array_scheme_t) :: one
        real(wp) :: v(2), t

        v = u0
        t = 0.0_wp
        call rk_scheme(one, name, 2, step)
        call one%advance(vanderpol, t, v)
        x_after = v(1)
    end function x_after

    !> Where a run at step dt (just above the stable step) at the largest mu
    !> first leaves the cycle (|x| > 2.2 or |y| > 10 mu): on the slow branch
    !> (|y| < mu at the step before) or during a fast jump (|y| >= mu).
    function failure_place(name, dt) result(place)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: dt
        character(:), allocatable :: place
        type(rk_array_scheme_t) :: scheme
        real(wp) :: u(2), u_before(2), time
        integer  :: i

        mu = mus(size(mus))
        call rk_scheme(scheme, name, 2, dt)
        u = [sqrt(3.0_wp), 0.0_wp]
        time = 0.0_wp
        place = "---"
        do i = 1, ceiling(t_stiff / dt)
            u_before = u
            call scheme%advance(vanderpol, time, u)
            if (.not. (abs(u(1)) <= 2.2_wp .and. abs(u(2)) <= 10 * mu)) then
                if (abs(u_before(2)) < mu) then
                    place = "slow branch"
                else
                    place = "fast jump"
                end if
                exit
            end if
        end do
        mu = mu1
    end function failure_place

    !> Eigenvalues of the Jacobian [[0, 1], [-mu (2 x y + 1), mu (1 - x^2)]] along
    !> the limit cycle for `mu_cycle`, at the points where the exact dynamics
    !> contracts (both eigenvalues with negative real part).
    function cycle_eigenvalues(mu_cycle) result(lambda)
        real(wp), intent(in) :: mu_cycle
        complex(wp), allocatable :: lambda(:)
        type(rk_array_scheme_t) :: scheme
        complex(wp) :: root, pair(2)
        real(wp) :: u(2), time, a, b, step
        integer  :: i

        allocate(lambda(0))
        mu = mu_cycle
        step = 1.0e-3_wp / sqrt(mu)
        call rk_scheme(scheme, "dp54", 2, step)
        u = [2.0_wp, 0.0_wp]
        time = 0.0_wp
        do i = 1, nint(8.0_wp / step)
            call scheme%advance(vanderpol, time, u)
            if (time < 5.0_wp .or. mod(i, 10) /= 0) cycle          ! after the transient, every 10th step
            a = -mu * (2 * u(1) * u(2) + 1)                         ! J = [[0, 1], [a, b]]: lambda^2 - b lambda - a = 0
            b = mu * (1 - u(1)**2)
            root = sqrt(cmplx(b * b + 4 * a, 0.0_wp, wp))
            pair = [(b + root) / 2, (b - root) / 2]
            if (all(real(pair) < 0)) lambda = [lambda, pair]
        end do
        mu = mu1
    end function cycle_eigenvalues

    !> Largest h with |R(h lambda)| <= 1 for every lambda in the list (bisection
    !> after a scan, as for h*).
    real(wp) function limit_along_cycle(tableau, lambda) result(h)
        type(rk_tableau_t), intent(in) :: tableau
        complex(wp),        intent(in) :: lambda(:)
        real(wp) :: lo, hi, scale
        integer  :: it

        scale = maxval(abs(lambda))
        lo = 0.0_wp
        do
            hi = lo + 1.0e-3_wp / scale
            if (unstable_at(tableau, lambda, hi)) exit
            lo = hi
        end do
        do it = 1, 60
            h = (lo + hi) / 2
            if (unstable_at(tableau, lambda, h)) then
                hi = h
            else
                lo = h
            end if
        end do
        h = lo
    end function limit_along_cycle

    !> Whether |R(step lambda)| > 1 for some lambda of the list.
    logical function unstable_at(tableau, lambda, step)
        type(rk_tableau_t), intent(in) :: tableau
        complex(wp),        intent(in) :: lambda(:)
        real(wp),           intent(in) :: step
        integer :: k

        unstable_at = .false.
        do k = 1, size(lambda)
            if (abs(stability_r(tableau, step * lambda(k))) > 1.0_wp) then
                unstable_at = .true.
                return
            end if
        end do
    end function unstable_at

    !> h*: the first h > 0 with |R(-h)| > 1, R the stability polynomial of the
    !> tableau (scan, then bisection).
    real(wp) function real_axis_limit(tableau) result(h)
        type(rk_tableau_t), intent(in) :: tableau
        real(wp) :: lo, hi
        integer  :: it

        lo = 0.0_wp
        do
            hi = lo + 1.0e-3_wp
            if (abs(stability_r(tableau, cmplx(-hi, 0.0_wp, wp))) > 1.0_wp) exit
            lo = hi
        end do
        do it = 1, 60
            h = (lo + hi) / 2
            if (abs(stability_r(tableau, cmplx(-h, 0.0_wp, wp))) > 1.0_wp) then
                hi = h
            else
                lo = h
            end if
        end do
        h = lo
    end function real_axis_limit

    !> Largest step (bisection on a log scale) for which the run from (sqrt(3), 0)
    !> over 0 <= t <= t_stiff stays bounded (|x| < 10).
    real(wp) function stable_step(name) result(h)
        character(*), intent(in) :: name
        real(wp) :: lo, hi
        integer  :: it

        lo = 1.0e-3_wp / mu
        hi = 10.0_wp / mu
        do it = 1, 40
            h = sqrt(lo * hi)
            if (bounded(name, h)) then
                lo = h
            else
                hi = h
            end if
        end do
        h = lo
    end function stable_step

    logical function bounded(name, dt)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: dt
        type(rk_array_scheme_t) :: scheme
        real(wp) :: u(2), time
        integer  :: i

        call rk_scheme(scheme, name, 2, dt)
        u = [sqrt(3.0_wp), 0.0_wp]
        time = 0.0_wp
        bounded = .true.
        do i = 1, ceiling(t_stiff / dt)
            call scheme%advance(vanderpol, time, u)
            if (.not. (abs(u(1)) < 10)) then
                bounded = .false.
                return
            end if
        end do
    end function bounded

    !> Data of the solution figures, all at mu = 1 except the limit cycles.
    subroutine write_figures()
        real(wp), parameter :: starts(2, 5) = reshape([0.1_wp, 0.0_wp, 3.0_wp, 3.0_wp, -3.0_wp, -2.5_wp, &
                                                       -0.5_wp, 1.0_wp, 2.6_wp, -3.2_wp], [2, 5])
        real(wp), parameter :: cycle_mus(4) = [1.0_wp, 4.0_wp, 16.0_wp, 64.0_wp]
        type(rk_array_scheme_t) :: scheme
        real(wp) :: u(2), d(2), time, x, y
        integer  :: i, j, k, unit

        mu = mu1
        ! Direction field: unit vectors, coloured by speed.
        open(newunit=unit, file=out // "field.dat", status="replace", action="write")
        write(unit, '(a)') "x y u v speed"
        do i = 0, 20
            do j = 0, 20
                x = -3.0_wp + 6.0_wp * i / 20
                y = -3.6_wp + 7.2_wp * j / 20
                call vanderpol(0.0_wp, [x, y], d, stage_t())
                write(unit, '(5es14.6)') x, y, d / max(norm2(d), 1.0e-12_wp), norm2(d)
            end do
        end do
        close(unit)
        ! Reference trajectory and trajectories from other starting points.
        call trajectory("dp54", 0.01_wp, t_plot, [sqrt(3.0_wp), 0.0_wp], "reference.dat", 5)
        do k = 1, size(starts, 2)
            call trajectory("dp54", 0.01_wp, 15.0_wp, starts(:, k), "start_" // itoa(k) // ".dat", 5)
        end do
        ! Numerical solutions at dt_plot.
        do k = 1, size(names)
            call trajectory(trim(names(k)), dt_plot, t_plot, [sqrt(3.0_wp), 0.0_wp], &
                            "solution_" // trim(names(k)) // ".dat", 1)
        end do
        ! Limit cycles for growing mu, after the transient: (t, x, y/sqrt(mu)).
        do k = 1, size(cycle_mus)
            mu = cycle_mus(k)
            open(newunit=unit, file=out // "cycle_" // itoa(k) // ".dat", status="replace", action="write")
            write(unit, '(a)') "t x y"
            call rk_scheme(scheme, "dp54", 2, 1.0e-3_wp / sqrt(mu))
            u = [2.0_wp, 0.0_wp]
            time = 0.0_wp
            do i = 1, nint(20.0_wp / (1.0e-3_wp / sqrt(mu)))
                call scheme%advance(vanderpol, time, u)
                if (time > 12.0_wp .and. mod(i, 4) == 0) write(unit, '(3es14.6)') time, u(1), u(2) / sqrt(mu)
            end do
            close(unit)
        end do
        mu = mu1
    end subroutine write_figures

    subroutine trajectory(name, dt, t_end, start, file, every)
        character(*), intent(in) :: name, file
        real(wp),     intent(in) :: dt, t_end, start(2)
        integer,      intent(in) :: every
        type(rk_array_scheme_t) :: scheme
        real(wp) :: u(2), time
        integer  :: i, unit

        open(newunit=unit, file=out // file, status="replace", action="write")
        write(unit, '(a)') "t x y"
        call rk_scheme(scheme, name, 2, dt)
        u = start
        time = 0.0_wp
        write(unit, '(3es14.6)') time, u
        do i = 1, nint(t_end / dt)
            call scheme%advance(vanderpol, time, u)
            if (.not. all(abs(u) < 1.0e3_wp)) exit
            if (mod(i, every) == 0) write(unit, '(3es14.6)') time, u
        end do
        close(unit)
    end subroutine trajectory

    !> Observed order: slope of the two finest levels with errors in (floor, 0.1),
    !> kept only if within 0.3 of the next coarser slope; 0 (a dash) otherwise.
    real(wp) function slope(dt, err) result(p)
        real(wp), intent(in) :: dt(:), err(:)
        real(wp) :: coarser
        integer  :: l

        p = 0.0_wp
        do l = size(dt), 3, -1
            if (err(l) > floor .and. err(l - 1) > floor .and. err(l - 2) < 1.0e-1_wp) then
                p = log(err(l - 1) / err(l)) / log(dt(l - 1) / dt(l))
                coarser = log(err(l - 2) / err(l - 1)) / log(dt(l - 2) / dt(l - 1))
                if (abs(p - coarser) > 0.3_wp) p = 0.0_wp
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

    subroutine write_parameters()
        integer :: u, j
        character(:), allocatable :: list, mlist

        list = ""
        do j = 1, size(names)
            if (j > 1) list = list // ", "
            list = list // "\texttt{" // trim(names(j)) // "}"
        end do
        mlist = ""
        do j = 1, size(mus)
            if (j > 1) mlist = mlist // ", "
            mlist = mlist // itoa(nint(mus(j)))
        end do
        open(newunit=u, file=out // "parameters.tex", status="replace", action="write")
        write(u, '(a)') "% Written by examples/vanderpol/main.f90: do not edit."
        write(u, '(a)') "\begin{tabular}{l>{\raggedright\arraybackslash}p{12.4cm}}"
        write(u, '(a)') "\toprule"
        write(u, '(a)') "parameter & $\mu = 1$ (equivalently $\varepsilon = \sqrt\mu = 1$) \\"
        write(u, '(a)') "initial condition & $(x, y) = (\sqrt 3, 0)$ \\"
        write(u, '(a)') "reference & \texttt{dp54}, $\Delta t = 10^{" // itoa(nint(log10(dt_reference))) // &
                        "}$; period $T = " // fmt(period_ref, 10) // "$ \\"
        write(u, '(a)') "convergence & error $\max_i|u_i - u_i^{\mathrm{ref}}|$ at $t = " // fmt(t_conv, 1) // &
                        "$, $\Delta t = " // fmt(dt_coarse, 2) // ", \dots, " // fmt(dt_coarse, 2) // "/" // &
                        itoa(2**(n_levels - 1)) // "$ \\"
        write(u, '(a)') "period & mean spacing of upward zero crossings " // itoa(first_cycle) // " to " // &
                        itoa(last_cycle) // " of $x$; each crossing solved to the method's accuracy \\"
        write(u, '(a)') "observed order & two finest levels above $10^{" // itoa(nint(log10(floor))) // &
                        "}$; kept if within 0.3 of the next slope \\"
        write(u, '(a)') "stiffness & largest $\Delta t$ keeping $|x| < 10$ on $0 \le t \le " // fmt(t_stiff, 1) // &
                        "$, $\mu = " // mlist // "$ \\"
        write(u, '(a)') "figures & solutions at $\Delta t = " // fmt(dt_plot, 2) // "$, $0 \le t \le " // &
                        fmt(t_plot, 1) // "$; limit cycles for $\mu = 1, 4, 16, 64$ \\"
        write(u, '(a)') "methods & " // list // " \\"
        write(u, '(a)') "interface & plain arrays, \texttt{rk\_array\_scheme\_t} \\"
        write(u, '(a)') "\bottomrule"
        write(u, '(a)') "\end{tabular}"
        close(u)
    end subroutine write_parameters

    function fmt(x, digits) result(s)
        real(wp), intent(in) :: x
        integer,  intent(in) :: digits
        character(:), allocatable :: s
        character(len=32) :: buffer, form

        write(form, '("(f0.", i0, ")")') digits
        write(buffer, form) x
        s = trim(adjustl(buffer))
        if (s(1:1) == ".") s = "0" // s
    end function fmt

    function itoa(i) result(s)
        integer, intent(in) :: i
        character(:), allocatable :: s
        character(len=16) :: buffer

        write(buffer, '(i0)') i
        s = trim(buffer)
    end function itoa

end program main
