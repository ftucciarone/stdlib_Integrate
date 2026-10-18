program main
    !! Lorenz system (Lorenz 1963; SciML prob_ode_lorenz):
    !!   x' = sigma (y - x),  y' = x (rho - z) - y,  z' = x y - beta z,
    !!   (sigma, rho, beta) = (10, 28, 8/3),  (x, y, z)(0) = (1, 0, 0).
    !! Chaotic: errors grow like exp(lambda t), lambda ~ 0.906 the largest
    !! Lyapunov exponent. Reference: dp54 at a tiny step, trusted up to the time
    !! it departs from a run at half that step. Tests:
    !!   - convergence: the error at t = t_conv as the step halves;
    !!   - predictability: T_div, the first time a run departs from the reference
    !!     by more than `departure`. With an error C dt^p exp(lambda t),
    !!     T_div = (ln(departure / C) + p ln(1/dt)) / lambda: a straight line in
    !!     ln(1/dt) of slope p/lambda, so each method measures lambda;
    !!   - lambda directly, from two reference runs started 1e-10 apart.
    !! Uses the plain-array schemes. Writes examples/lorenz/data/:
    !! parameters.tex, results.tex, attractor.dat (t, x, y, z), reference.dat
    !! (t, x), solution_<method>.dat (t, x at dt_show), growth_<method>.dat (t,
    !! |u - u_ref| at dt_show), twins.dat (t, distance of the two nearby runs),
    !! error_<method>.dat (dt, error at t_conv), divergence_<method>.dat (ln(1/dt),
    !! T_div), read by summary.tex. Run from the project root.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: rk_scheme, stage_t, rk_array_scheme_t, rk_method_names, rk_tableau, rk_tableau_t
    implicit none

    character(*), parameter :: out = "examples/lorenz/data/"
    real(wp),     parameter :: sigma = 10.0_wp, rho = 28.0_wp, beta = 8.0_wp / 3.0_wp
    real(wp),     parameter :: u0(3) = [1.0_wp, 0.0_wp, 0.0_wp]
    real(wp),     parameter :: dt_base = 0.02_wp              !! coarsest step of the studies
    integer,      parameter :: n_levels = 8                   !! dt = 0.02, 0.01, ..., 0.02/128
    integer,      parameter :: ref_refine = 16                !! reference step = finest / ref_refine
    real(wp),     parameter :: t_window = 40.0_wp             !! the reference covers 0 <= t <= t_window
    real(wp),     parameter :: t_conv = 1.0_wp                !! convergence: error at this time
    real(wp),     parameter :: departure = 1.0_wp             !! T_div threshold
    real(wp),     parameter :: dt_show = 0.01_wp              !! step of the solution figures
    real(wp),     parameter :: floor = 1.0e-11_wp

    character(len=8), allocatable :: names(:)
    type(rk_tableau_t) :: tableau
    real(wp), allocatable :: reference(:,:)                   !! the reference at every finest step
    real(wp) :: dt_fine, t_trust, lambda_twins, dt(n_levels), err(n_levels), tdiv(n_levels)
    real(wp) :: observed, slope_div, lambda_method
    integer  :: n_fine, m, l, table, dat, n_fit

    names = rk_method_names()
    call execute_command_line("mkdir -p " // out)
    dt_fine = dt_base / 2**(n_levels - 1)
    n_fine = nint(t_window / dt_fine)

    ! Reference at every finest step, and how long it can be trusted.
    allocate(reference(3, 0:n_fine))
    call reference_run(dt_fine / ref_refine, reference, t_trust)
    lambda_twins = twins_lambda()
    print '(a, f8.3, a, f8.4)', "reference trusted to t = ", t_trust, ";  lambda from twins = ", lambda_twins
    call write_parameters()
    call write_figures()

    open(newunit=table, file=out // "results.tex", status="replace", action="write")
    write(table, '(a)') "% Written by examples/lorenz/main.f90: do not edit."
    write(table, '(a)') "\begin{tabular}{lcccccc}"
    write(table, '(a)') "\toprule"
    write(table, '(a)') "name & order & observed & $T_{\mathrm{div}}$ at $\Delta t = " // fmt(dt_show, 2) // &
                        "$ & slope of $T_{\mathrm{div}}$ & $p_{\mathrm{eff}}$ & $\lambda = p_{\mathrm{eff}}/\mathrm{slope}$ \\"
    write(table, '(a)') "\midrule"
    print '(a10, a6, a9, 3a10)', "name", "order", "obs.", "Tdiv", "slope", "lambda"
    do m = 1, size(names)
        tableau = rk_tableau(trim(names(m)))
        open(newunit=dat, file=out // "error_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat, '(a)') "dt error"
        do l = 1, n_levels
            dt(l) = dt_base / 2**(l - 1)
            call compare(trim(names(m)), l, err(l), tdiv(l))
            if (err(l) > 1.0e-15_wp .and. err(l) < 1.0e3_wp) write(dat, '(2es14.6)') dt(l), err(l)
        end do
        close(dat)
        observed = slope(dt, err)
        ! T_div against ln(1/dt): least-squares slope over the levels that
        ! diverge well inside the trusted window and after the first second.
        open(newunit=dat, file=out // "divergence_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat, '(a)') "logdt tdiv"
        do l = 1, n_levels
            if (tdiv(l) < t_trust) write(dat, '(2es14.6)') log(1 / dt(l)), tdiv(l)
        end do
        close(dat)
        call fit(log(1 / dt), tdiv, tdiv > 1.0_wp .and. tdiv < t_trust - 3, slope_div, n_fit)
        lambda_method = merge(effective_order(trim(names(m))) / slope_div, 0.0_wp, n_fit >= 3)
        call theory_line(trim(names(m)), log(1 / dt), tdiv, tdiv > 1.0_wp .and. tdiv < t_trust - 3, &
                         effective_order(trim(names(m))) / lambda_twins)
        print '(a10, i6, f9.3, 3f10.3)', trim(names(m)), tableau%order, observed, &
            divergence_time(trim(names(m)), dt_show), slope_div, lambda_method
        write(table, '("\texttt{", a, "} & ", i0, " & ", a, " & ", f5.2, " & ", a, " & ", i0, " & ", a, " \\")') &
            trim(names(m)), tableau%order, order_text(observed), divergence_time(trim(names(m)), dt_show), &
            value_text(slope_div, n_fit >= 3), effective_order(trim(names(m))), value_text(lambda_method, n_fit >= 3)
    end do
    write(table, '(a)') "\bottomrule"
    write(table, '(a)') "\end{tabular}"
    close(table)

contains

    subroutine lorenz(t, u, dudt, stage)
        real(wp), intent(in)  :: t
        real(wp), intent(in)  :: u(:)
        real(wp), intent(out) :: dudt(:)
        type(stage_t), intent(in) :: stage

        dudt(1) = sigma * (u(2) - u(1))
        dudt(2) = u(1) * (rho - u(3)) - u(2)
        dudt(3) = u(1) * u(2) - beta * u(3)
    end subroutine lorenz

    !> The reference (dp54 at step h) at every finest step, and t_trust: the
    !> first time it departs by 1e-6 from the same run at step 2h, plus the time
    !> its error needs to grow from 1e-6 to 1e-2 (ln(1e4)/0.9), a hundredth of
    !> `departure`: until then the reference decides T_div correctly.
    subroutine reference_run(h, keep, t_trust)
        real(wp), intent(in)  :: h
        real(wp), intent(out) :: keep(:, 0:), t_trust
        type(rk_array_scheme_t) :: fine, coarse
        real(wp) :: u(3), v(3), tu, tv
        integer  :: i, j

        call rk_scheme(fine, "dp54", 3, h)
        call rk_scheme(coarse, "dp54", 3, 2 * h)
        u = u0
        v = u0
        tu = 0.0_wp
        tv = 0.0_wp
        keep(:, 0) = u
        t_trust = t_window
        do i = 1, size(keep, 2) - 1
            do j = 1, ref_refine / 2
                call fine%advance(lorenz, tu, u)
                call fine%advance(lorenz, tu, u)
                call coarse%advance(lorenz, tv, v)
            end do
            keep(:, i) = u
            if (t_trust == t_window .and. norm2(u - v) > 1.0e-6_wp) &
                t_trust = min(i * dt_fine + log(1.0e4_wp) / 0.9_wp, t_window)
        end do
    end subroutine reference_run

    !> Largest Lyapunov exponent, by Benettin's method: two dp54 runs a distance
    !> d0 = 1e-8 apart; every time unit the growth ln(d/d0) is added up and the
    !> second run is pulled back to distance d0 along the separation. lambda is
    !> the mean growth rate over t_lyapunov (after one discarded unit). Also
    !> writes twins.dat: the distance of two runs started 1e-10 apart, never
    !> renormalised (the figure).
    real(wp) function twins_lambda() result(lam)
        real(wp), parameter :: d0 = 1.0e-8_wp, h = 1.0e-3_wp, t_lyapunov = 500.0_wp
        type(rk_array_scheme_t) :: s1, s2
        real(wp) :: u(3), v(3), t1, t2, d, total
        integer  :: i, j, unit

        call rk_scheme(s1, "dp54", 3, h)
        call rk_scheme(s2, "dp54", 3, h)
        u = u0
        v = u0 + [1.0e-10_wp, 0.0_wp, 0.0_wp]
        t1 = 0.0_wp
        t2 = 0.0_wp
        open(newunit=unit, file=out // "twins.dat", status="replace", action="write")
        write(unit, '(a)') "t distance"
        do i = 1, 40000
            call s1%advance(lorenz, t1, u)
            call s2%advance(lorenz, t2, v)
            if (mod(i, 20) == 0) write(unit, '(2es14.6)') i * h, norm2(u - v)
        end do
        close(unit)
        ! Benettin: start on the attractor (the state reached above).
        v = u + [d0, 0.0_wp, 0.0_wp]
        total = 0.0_wp
        do i = 1, nint(t_lyapunov) + 1
            do j = 1, nint(1.0_wp / h)
                call s1%advance(lorenz, t1, u)
                call s2%advance(lorenz, t2, v)
            end do
            d = norm2(v - u)
            if (i > 1) total = total + log(d / d0)
            v = u + (v - u) * (d0 / d)
        end do
        lam = total / t_lyapunov
    end function twins_lambda

    !> Run the method at level l (dt = dt_base / 2^(l-1)) against the reference:
    !> `err` = error at t_conv, `tdiv` = first time the distance exceeds
    !> `departure` (t_window if never).
    subroutine compare(name, l, err, tdiv)
        character(*), intent(in)  :: name
        integer,      intent(in)  :: l
        real(wp),     intent(out) :: err, tdiv
        type(rk_array_scheme_t) :: scheme
        real(wp) :: u(3), time, h
        integer  :: i, stride

        h = dt_base / 2**(l - 1)
        stride = 2**(n_levels - l)                      ! finest steps per step of this level
        call rk_scheme(scheme, name, 3, h)
        u = u0
        time = 0.0_wp
        err = huge(1.0_wp)
        tdiv = t_window
        do i = 1, nint(t_window / h)
            call scheme%advance(lorenz, time, u)
            if (i * stride == nint(t_conv / dt_fine)) err = maxval(abs(u - reference(:, i * stride)))
            if (norm2(u - reference(:, i * stride)) > departure) then
                tdiv = i * h
                exit
            end if
        end do
    end subroutine compare

    real(wp) function divergence_time(name, h) result(td)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: h
        real(wp) :: e
        integer  :: l

        l = nint(log(dt_base / h) / log(2.0_wp)) + 1
        call compare(name, l, e, td)
    end function divergence_time

    !> The order that governs error accumulation: the method's order, except
    !> rk3ws, whose second-order error term does not accumulate (effective
    !> order 3, see the summary).
    integer function effective_order(name) result(p)
        character(*), intent(in) :: name

        p = 0
        select case (name)
        case ("rk3ws");            p = 3
        case ("rk4");              p = 4
        case ("dp54");             p = 5
        end select
    end function effective_order

    !> divtheory_<name>.dat: the line of slope p_eff/lambda (lambda from
    !> Benettin) through the centroid of the selected points (x, y).
    subroutine theory_line(name, x, y, use, s)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: x(:), y(:), s
        logical,      intent(in) :: use(:)
        real(wp) :: mx, my, x_lo, x_hi
        integer  :: unit

        open(newunit=unit, file=out // "divtheory_" // name // ".dat", status="replace", action="write")
        write(unit, '(a)') "logdt tdiv"
        if (count(use) >= 2) then
            mx = sum(x, use) / count(use)
            my = sum(y, use) / count(use)
            x_lo = minval(x, use) - 0.3_wp
            x_hi = maxval(x, use) + 0.3_wp
            write(unit, '(2es14.6)') x_lo, my + s * (x_lo - mx)
            write(unit, '(2es14.6)') x_hi, my + s * (x_hi - mx)
        end if
        close(unit)
    end subroutine theory_line

    !> Least-squares slope of y against x over the selected points.
    subroutine fit(x, y, use, s, n)
        real(wp), intent(in)  :: x(:), y(:)
        logical,  intent(in)  :: use(:)
        real(wp), intent(out) :: s
        integer,  intent(out) :: n
        real(wp) :: mx, my

        n = count(use)
        s = 0.0_wp
        if (n < 2) return
        mx = sum(x, use) / n
        my = sum(y, use) / n
        s = sum((x - mx) * (y - my), use) / sum((x - mx)**2, use)
    end subroutine fit

    subroutine write_figures()
        type(rk_array_scheme_t) :: scheme
        real(wp) :: u(3), time
        integer  :: i, unit, m, l, stride

        ! The attractor: a long dp54 run.
        open(newunit=unit, file=out // "attractor.dat", status="replace", action="write")
        write(unit, '(a)') "t x y z"
        call rk_scheme(scheme, "dp54", 3, 1.0e-3_wp)
        u = u0
        time = 0.0_wp
        do i = 1, 60000
            call scheme%advance(lorenz, time, u)
            if (mod(i, 6) == 0) write(unit, '(4es14.6)') time, u
        end do
        close(unit)
        ! Reference x(t), every 0.01.
        open(newunit=unit, file=out // "reference.dat", status="replace", action="write")
        write(unit, '(a)') "t x"
        stride = nint(0.01_wp / dt_fine)
        do i = 0, n_fine, stride
            write(unit, '(2es14.6)') i * dt_fine, reference(1, i)
        end do
        close(unit)
        ! x(t) and the error growth of every method at dt_show.
        l = nint(log(dt_base / dt_show) / log(2.0_wp)) + 1
        stride = 2**(n_levels - l)
        do m = 1, size(names)
            call rk_scheme(scheme, trim(names(m)), 3, dt_show)
            u = u0
            time = 0.0_wp
            open(newunit=unit, file=out // "solution_" // trim(names(m)) // ".dat", status="replace", action="write")
            write(unit, '(a)') "t x error"
            do i = 1, nint(t_window / dt_show)
                call scheme%advance(lorenz, time, u)
                write(unit, '(3es14.6)') i * dt_show, u(1), max(norm2(u - reference(:, i * stride)), 1.0e-16_wp)
            end do
            close(unit)
        end do
    end subroutine write_figures

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

        s = value_text(p, p /= 0.0_wp)
    end function order_text

    function value_text(x, valid) result(s)
        real(wp), intent(in) :: x
        logical,  intent(in) :: valid
        character(:), allocatable :: s
        character(len=12) :: buffer

        if (.not. valid) then
            s = "---"
        else
            write(buffer, '(f6.2)') x
            s = trim(adjustl(buffer))
        end if
    end function value_text

    subroutine write_parameters()
        integer :: u, j
        character(:), allocatable :: list

        list = ""
        do j = 1, size(names)
            if (j > 1) list = list // ", "
            list = list // "\texttt{" // trim(names(j)) // "}"
        end do
        open(newunit=u, file=out // "parameters.tex", status="replace", action="write")
        write(u, '(a)') "% Written by examples/lorenz/main.f90: do not edit."
        write(u, '(a)') "\begin{tabular}{l>{\raggedright\arraybackslash}p{12.4cm}}"
        write(u, '(a)') "\toprule"
        write(u, '(a)') "parameters & $(\sigma, \rho, \beta) = (10, 28, 8/3)$ \\"
        write(u, '(a)') "initial condition & $(x, y, z) = (1, 0, 0)$ \\"
        write(u, '(a)') "reference & \texttt{dp54}, $\Delta t = " // fmt(dt_fine / ref_refine * 1.0e6_wp, 3) // &
                        "\times10^{-6}$, $0 \le t \le " // fmt(t_window, 0) // "$; trusted to $t = " // fmt(t_trust, 1) // &
                        "$ (departs by $10^{-6}$ from the run at twice the step, plus the time to grow to $10^{-2}$) \\"
        write(u, '(a)') "convergence & error $\max_i|u_i - u_i^{\mathrm{ref}}|$ at $t = " // fmt(t_conv, 0) // &
                        "$, $\Delta t = " // fmt(dt_base, 2) // ", \dots, " // fmt(dt_base, 2) // "/" // &
                        itoa(2**(n_levels - 1)) // "$; order from the two finest levels above $10^{" // &
                        itoa(nint(log10(floor))) // "}$, kept if within 0.3 of the next slope \\"
        write(u, '(a)') "predictability & $T_{\mathrm{div}}$: first $t$ with $|u - u^{\mathrm{ref}}| > " // &
                        fmt(departure, 0) // "$; slope: least squares of $T_{\mathrm{div}}$ against $\ln(1/\Delta t)$" // &
                        " over the levels with $1 < T_{\mathrm{div}} < t_{\mathrm{trust}} - 3$ (at least 3) \\"
        write(u, '(a)') "Lyapunov exponent & Benettin's method: two \texttt{dp54} runs ($\Delta t = 10^{-3}$), " // &
                        "separation $10^{-8}$ restored every time unit, mean growth rate over 500 units: " // &
                        "$\lambda = " // fmt(lambda_twins, 3) // "$ \\"
        write(u, '(a)') "figures & attractor: \texttt{dp54}, $0 \le t \le 60$; solutions at $\Delta t = " // &
                        fmt(dt_show, 2) // "$ \\"
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
        if (digits == 0 .and. s(len(s):len(s)) == ".") s = s(1:len(s) - 1)
    end function fmt

    function itoa(i) result(s)
        integer, intent(in) :: i
        character(:), allocatable :: s
        character(len=16) :: buffer

        write(buffer, '(i0)') i
        s = trim(buffer)
    end function itoa

end program main
