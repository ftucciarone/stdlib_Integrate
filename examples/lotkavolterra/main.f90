program main
    !! Lotka-Volterra predator-prey model (SciML prob_ode_lotkavolterra):
    !!   x' = a x - b x y,   y' = -c y + d x y,   (a, b, c, d) = (1.5, 1, 3, 1),  x(0) = y(0) = 1.
    !! No closed form in time, but V = d x - c ln x + b y - a ln y is constant,
    !! so the exact orbit is the closed level curve V = V0, computed here. Tests:
    !!   - convergence: max |V - V0| over 0 <= t <= t_conv as dt halves, at the
    !!     method's order (V is a smooth function of the solution);
    !!   - long time: |V - V0| over 0 <= t <= t_long at dt_long; Runge-Kutta
    !!     methods do not keep V, so the orbit drifts off the level curve.
    !! Uses the plain-array schemes. Writes examples/lotkavolterra/data/:
    !! parameters.tex, results.tex, orbit.dat (the exact orbit), reference.dat
    !! (t, x, y with dp54 at dt_reference), solution_<method>.dat (t, x, y at
    !! dt_plot), drift_<method>.dat (t, running max |V - V0| at dt_long), error_<method>.dat
    !! (dt, max |V - V0|), read by summary.tex. Run from the project root.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: rk_scheme, stage_t, rk_array_scheme_t, rk_method_names, rk_tableau, rk_tableau_t
    implicit none

    character(*), parameter :: out = "examples/lotkavolterra/data/"
    real(wp),     parameter :: a = 1.5_wp, b = 1.0_wp, c = 3.0_wp, d = 1.0_wp
    real(wp),     parameter :: x0 = 1.0_wp, y0 = 1.0_wp
    real(wp),     parameter :: t_plot = 20.0_wp, dt_plot = 0.1_wp        !! solution figure
    real(wp),     parameter :: dt_reference = 1.0e-4_wp                  !! reference curve (dp54)
    real(wp),     parameter :: t_conv = 10.0_wp, dt_coarse = 0.1_wp      !! convergence study
    integer,      parameter :: n_levels = 8
    real(wp),     parameter :: t_long = 1000.0_wp, dt_long = 0.05_wp     !! long-time drift
    real(wp),     parameter :: round_off = 1.0e-12_wp

    character(len=8), allocatable :: names(:)
    type(rk_tableau_t) :: tableau
    real(wp) :: v0, dt(n_levels), err(n_levels), observed, final_drift, survived
    integer  :: m, l, table, dat

    names = rk_method_names()
    v0 = invariant([x0, y0])
    call execute_command_line("mkdir -p " // out)
    call write_parameters()
    call write_orbit()
    call trajectory("dp54", dt_reference, t_plot, "reference.dat", every=500)
    do m = 1, size(names)
        call trajectory(trim(names(m)), dt_plot, t_plot, "solution_" // trim(names(m)) // ".dat", every=1)
    end do

    open(newunit=table, file=out // "results.tex", status="replace", action="write")
    write(table, '(a)') "% Written by examples/lotkavolterra/main.f90: do not edit."
    write(table, '(a)') "\begin{tabular}{llcccc}"
    write(table, '(a)') "\toprule"
    write(table, '(a)') "method & name & order & observed & $|V - V_0|$ at $t = " // fmt(t_long) // "$ & " // &
                        "positive until \\"
    write(table, '(a)') "\midrule"
    print '(a20, a10, a6, a9, 2a14)', "method", "name", "order", "obs.", "drift t_long", "survives to"
    do m = 1, size(names)
        tableau = rk_tableau(trim(names(m)))
        open(newunit=dat, file=out // "error_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat, '(a)') "dt error"
        do l = 1, n_levels
            dt(l) = dt_coarse / 2**(l - 1)
            call drift_run(trim(names(m)), dt(l), t_conv, err(l), survived)
            if (err(l) > 1.0e-15_wp .and. err(l) < huge(1.0_wp)) write(dat, '(2es14.6)') dt(l), err(l)  ! not if it left x, y > 0
        end do
        close(dat)
        observed = slope(dt, err)
        call drift_run(trim(names(m)), dt_long, t_long, final_drift, survived, "drift_" // trim(names(m)) // ".dat")
        print '(a20, a10, i6, f9.3, 2es14.4)', tableau%title, trim(names(m)), tableau%order, observed, final_drift, survived
        write(table, '(a, " & \texttt{", a, "} & ", i0, " & ", f5.2, " & ", a, " & ", a, " \\")') &
            tableau%title, trim(names(m)), tableau%order, observed, sci(final_drift), survived_text(survived)
    end do
    write(table, '(a)') "\bottomrule"
    write(table, '(a)') "\end{tabular}"
    close(table)

contains

    subroutine lotka_volterra(t, u, dudt, stage)
        real(wp), intent(in)  :: t
        real(wp), intent(in)  :: u(:)
        real(wp), intent(out) :: dudt(:)
        type(stage_t), intent(in) :: stage

        dudt(1) = a * u(1) - b * u(1) * u(2)
        dudt(2) = -c * u(2) + d * u(1) * u(2)
    end subroutine lotka_volterra

    !> V = d x - c ln x + b y - a ln y, constant along exact solutions.
    pure real(wp) function invariant(u)
        real(wp), intent(in) :: u(2)

        invariant = d * u(1) - c * log(u(1)) + b * u(2) - a * log(u(2))
    end function invariant

    !> Every `every` steps of size dt up to t_end, (t, x, y) to data/<file>.
    subroutine trajectory(name, dt, t_end, file, every)
        character(*), intent(in) :: name, file
        real(wp),     intent(in) :: dt, t_end
        integer,      intent(in) :: every

        type(rk_array_scheme_t) :: scheme
        real(wp) :: u(2), time
        integer  :: i, unit

        open(newunit=unit, file=out // file, status="replace", action="write")
        write(unit, '(a)') "t x y"
        call rk_scheme(scheme, name, 2, dt)
        u = [x0, y0]
        time = 0.0_wp
        write(unit, '(3es16.8)') time, u
        do i = 1, nint(t_end / dt)
            call scheme%advance(lotka_volterra, time, u)
            if (any(u <= 0)) exit                     ! left the positive quadrant
            if (mod(i, every) == 0) write(unit, '(3es16.8)') time, u
        end do
        close(unit)
    end subroutine trajectory

    !> Integrate to t_end with step dt: `drift` = max |V - V0| over the run,
    !> `survived` = time the solution stayed in x, y > 0 (t_end if always).
    !> With `file`, write (t, max |V - V0| so far) at times spread geometrically.
    subroutine drift_run(name, dt, t_end, drift, survived, file)
        character(*), intent(in)  :: name
        real(wp),     intent(in)  :: dt, t_end
        real(wp),     intent(out) :: drift, survived
        character(*), intent(in), optional :: file

        type(rk_array_scheme_t) :: scheme
        real(wp) :: u(2), time
        integer  :: i, n, unit, next

        n = nint(t_end / dt)
        call rk_scheme(scheme, name, 2, dt)
        u = [x0, y0]
        time = 0.0_wp
        drift = 0.0_wp
        survived = t_end
        next = 1
        if (present(file)) then
            open(newunit=unit, file=out // file, status="replace", action="write")
            write(unit, '(a)') "t drift"
        end if
        do i = 1, n
            call scheme%advance(lotka_volterra, time, u)
            if (any(u <= 0) .or. any(u /= u)) then    ! left the positive quadrant (or overflowed)
                survived = time
                drift = huge(1.0_wp)
                exit
            end if
            drift = max(drift, abs(invariant(u) - v0))
            if (present(file) .and. (i == next .or. i == n)) then
                write(unit, '(2es14.6)') time, max(drift, 1.0e-16_wp)       ! running maximum
                next = max(next + 1, nint(next * 1.05_wp))
            end if
        end do
        if (present(file)) close(unit)
    end subroutine drift_run

    !> The exact orbit, the level curve V(x, y) = V0: for x between the two
    !> roots of d x - c ln x = V0 - min_y(b y - a ln y), the two roots in y of
    !> b y - a ln y = V0 - d x + c ln x, one below and one above y = a/b.
    subroutine write_orbit()
        integer, parameter :: n = 400
        real(wp) :: x_lo, x_hi, x(0:n), y_lo(0:n), y_hi(0:n)
        integer  :: i, unit

        x_lo = root(hx, 1.0e-12_wp, c / d)
        x_hi = root(hx, c / d, 100.0_wp)
        do i = 0, n
            x(i) = x_lo + (x_hi - x_lo) * (1 - cos(acos(-1.0_wp) * i / n)) / 2   ! dense near the ends
            y_lo(i) = root_y(x(i), 1.0e-12_wp, a / b)
            y_hi(i) = root_y(x(i), a / b, 100.0_wp)
        end do
        open(newunit=unit, file=out // "orbit.dat", status="replace", action="write")
        write(unit, '(a)') "x y"
        do i = 0, n
            write(unit, '(2es16.8)') x(i), y_hi(i)
        end do
        do i = n, 0, -1
            write(unit, '(2es16.8)') x(i), y_lo(i)
        end do
        close(unit)
    end subroutine write_orbit

    !> d x - c ln x - (V0 - min_y(b y - a ln y)): zero at the orbit's extreme x.
    real(wp) function hx(s)
        real(wp), intent(in) :: s

        hx = d * s - c * log(s) - (v0 - (a - a * log(a / b)))
    end function hx

    !> Root in y of b y - a ln y = V0 - d xi + c ln xi on [lo, hi], by bisection.
    real(wp) function root_y(xi, lo, hi)
        real(wp), intent(in) :: xi, lo, hi
        real(wp) :: l, h, mid, rhs
        integer  :: it

        rhs = v0 - d * xi + c * log(xi)
        l = lo
        h = hi
        do it = 1, 200
            mid = (l + h) / 2
            if (sign(1.0_wp, b * l - a * log(l) - rhs) == sign(1.0_wp, b * mid - a * log(mid) - rhs)) then
                l = mid
            else
                h = mid
            end if
        end do
        root_y = (l + h) / 2
    end function root_y

    !> Bisection root of f on [lo, hi] (f changes sign there).
    real(wp) function root(f, lo, hi)
        interface
            real(wp) function f(s)
                import wp
                real(wp), intent(in) :: s
            end function f
        end interface
        real(wp), intent(in) :: lo, hi
        real(wp) :: l, h, mid
        integer  :: it

        l = lo
        h = hi
        do it = 1, 200
            mid = (l + h) / 2
            if (sign(1.0_wp, f(l)) == sign(1.0_wp, f(mid))) then
                l = mid
            else
                h = mid
            end if
        end do
        root = (l + h) / 2
    end function root

    !> Slope of log(err) against log(dt) between the two finest levels whose
    !> errors are both above round-off.
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

    function survived_text(t) result(s)
        real(wp), intent(in) :: t
        character(:), allocatable :: s

        if (t >= t_long) then
            s = "end"
        else
            s = "$t = " // fmt(t) // "$"
        end if
    end function survived_text

    subroutine write_parameters()
        integer :: u, k
        character(:), allocatable :: list

        list = ""
        do k = 1, size(names)
            if (k > 1) list = list // ", "
            list = list // "\texttt{" // trim(names(k)) // "}"
        end do
        open(newunit=u, file=out // "parameters.tex", status="replace", action="write")
        write(u, '(a)') "% Written by examples/lotkavolterra/main.f90: do not edit."
        write(u, '(a)') "\begin{tabular}{ll}"
        write(u, '(a)') "\toprule"
        write(u, '(a)') "parameters & $(a, b, c, d) = (" // fmt(a) // ", " // fmt(b) // ", " // fmt(c) // ", " // &
                        fmt(d) // ")$ \\"
        write(u, '(a)') "initial condition & $x(0) = " // fmt(x0) // "$, $y(0) = " // fmt(y0) // "$, $V_0 = " // &
                        fmt(v0) // "$ \\"
        write(u, '(a)') "equilibrium & $(c/d, a/b) = (" // fmt(c / d) // ", " // fmt(a / b) // ")$ \\"
        write(u, '(a)') "solution figure & $\Delta t = " // fmt(dt_plot) // "$, $0 \le t \le " // fmt(t_plot) // "$ \\"
        write(u, '(a)') "reference curve & \texttt{dp54}, $\Delta t = 10^{" // itoa(nint(log10(dt_reference))) // "}$ \\"
        write(u, '(a)') "convergence & $\max_t |V - V_0|$, $0 \le t \le " // fmt(t_conv) // "$, $\Delta t = " // &
                        fmt(dt_coarse) // ", \dots, " // fmt(dt_coarse) // "/" // itoa(2**(n_levels - 1)) // "$ \\"
        write(u, '(a)') "observed order & slope of the two finest levels with errors in $(10^{" // &
                        itoa(nint(log10(round_off))) // "}, 1)$ \\"
        write(u, '(a)') "long time & $|V - V_0|$ for $0 \le t \le " // fmt(t_long) // "$ at $\Delta t = " // &
                        fmt(dt_long) // "$ \\"
        write(u, '(a)') "methods & " // list // " \\"
        write(u, '(a)') "interface & plain arrays, \texttt{rk\_array\_scheme\_t} \\"
        write(u, '(a)') "\bottomrule"
        write(u, '(a)') "\end{tabular}"
        close(u)
    end subroutine write_parameters

    function sci(x) result(s)
        real(wp), intent(in) :: x
        character(:), allocatable :: s
        character(len=32) :: mantissa
        integer :: e

        if (x >= huge(1.0_wp)) then
            s = "---"
            return
        end if
        if (x == 0.0_wp) then
            s = "$0$"
            return
        end if
        e = floor(log10(abs(x)))
        write(mantissa, '(f4.2)') abs(x) / 10.0_wp**e
        s = "$" // trim(mantissa) // "\times10^{" // itoa(e) // "}$"
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
