program main
    !! Heat equation u_t = u_xx on [0, 1], u(0) = u(1) = 0, by the method of
    !! lines (heat_model: central differences on n interior nodes, the boundary
    !! filled by `apply_boundaries`), u(x, 0) = sin(pi x) + sin(4 pi x)/2 + sin(9 pi x)/4.
    !! The semi-discrete system u' = A u is solved exactly by its sine modes,
    !!   u_j(t) = sum_k a_k exp(lambda_k t) sin(k pi x_j),
    !!   lambda_k = -(4/dx^2) sin^2(k pi dx / 2),
    !! and the PDE by the same sum with -(k pi)^2 in place of lambda_k. Tests:
    !!   - time convergence at n = n_conv: the error against the exact
    !!     semi-discrete solution at t_conv, as dt halves;
    !!   - stability: the largest stable dt as the grid is refined, against
    !!     h*/|lambda_n| ~ h* dx^2/4, h* the method's real-axis limit;
    !!   - total error against the PDE at t_conv, with dt at 0.9 of each method's
    !!     limit: the space error O(dx^2) decides, not the method.
    !! Uses the model classes and rk_scheme (with `apply_boundaries`).
    !! Writes examples/heat/data/: parameters.tex, results.tex, spacetime.dat,
    !! profiles.dat, error_<method>.dat, stable_<method>.dat, pde_<method>.dat,
    !! read by summary.tex. Run from the project root.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: status_t, time_scheme_c, rk_scheme, rk_method_names, rk_tableau, rk_tableau_t
    use heat_model, only: heat_model_t, heat_state_t, make_heat_model
    use heat_stability, only: stability_r
    implicit none

    character(*), parameter :: out = "examples/heat/data/"
    real(wp),     parameter :: pi = acos(-1.0_wp)
    integer,      parameter :: modes(3) = [1, 4, 9]
    real(wp),     parameter :: amplitudes(3) = [1.0_wp, 0.5_wp, 0.25_wp]
    integer,      parameter :: n_conv = 63                     !! grid of the time convergence study
    real(wp),     parameter :: t_conv = 0.02_wp
    real(wp),     parameter :: dt_coarse = 1.0e-4_wp           !! stable for every method at n_conv
    integer,      parameter :: n_levels = 7
    integer,      parameter :: grids(5) = [15, 31, 63, 127, 255]
    integer,      parameter :: n_probe = 3000                   !! steps of a stability probe
    real(wp),     parameter :: floor = 1.0e-13_wp

    character(len=8), allocatable :: names(:)
    type(heat_model_t), target :: model
    type(rk_tableau_t) :: tableau
    real(wp) :: dt(n_levels), err(n_levels), observed, h_star, measured, predicted, worst_ratio
    real(wp) :: dx_g(size(grids)), pde_err(size(grids)), pde_order
    integer  :: m, l, g, table, dat, dat2

    names = rk_method_names()
    call execute_command_line("mkdir -p " // out)
    call write_parameters()
    call write_figures()

    open(newunit=table, file=out // "results.tex", status="replace", action="write")
    write(table, '(a)') "% Written by examples/heat/main.f90: do not edit."
    write(table, '(a)') "\begin{tabular}{lccccc}"
    write(table, '(a)') "\toprule"
    write(table, '(a)') "& & \multicolumn{2}{c}{time convergence, $n = " // itoa(n_conv) // "$} & stability & PDE \\"
    write(table, '(a)') "\cmidrule(lr){3-4}"
    write(table, '(a)') "name & $h^*$ & order & observed & measured / predicted & order in $\Delta x$ \\"
    write(table, '(a)') "\midrule"
    print '(a10, a8, a6, a9, 2a12)', "name", "h*", "order", "obs.", "meas/pred", "PDE order"
    do m = 1, size(names)
        tableau = rk_tableau(trim(names(m)))
        h_star = real_axis_limit(tableau)
        ! Time convergence against the semi-discrete solution.
        model = make_heat_model(n_conv)
        open(newunit=dat, file=out // "error_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat, '(a)') "dt error"
        do l = 1, n_levels
            dt(l) = dt_coarse / 2**(l - 1)
            err(l) = maxval(abs(integrate(trim(names(m)), dt(l), t_conv) - semi_discrete(t_conv)))
            if (err(l) > 1.0e-16_wp) write(dat, '(2es14.6)') dt(l), err(l)
        end do
        close(dat)
        observed = slope(dt, err)
        ! Stability against the grid, and the PDE error at 0.9 of the limit.
        open(newunit=dat, file=out // "stable_" // trim(names(m)) // ".dat", status="replace", action="write")
        open(newunit=dat2, file=out // "pde_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat, '(a)') "dx measured predicted"
        write(dat2, '(a)') "dx error"
        worst_ratio = 1.0_wp
        do g = 1, size(grids)
            model = make_heat_model(grids(g))
            predicted = h_star / abs(eigenvalue(grids(g)))
            measured = stable_step(trim(names(m)), predicted)
            if (abs(measured / predicted - 1) > abs(worst_ratio - 1)) worst_ratio = measured / predicted
            write(dat, '(3es14.6)') model%dx, measured, predicted
            dx_g(g) = model%dx
            pde_err(g) = maxval(abs(integrate(trim(names(m)), t_conv / ceiling(t_conv / (0.9_wp * predicted)), t_conv) &
                                    - pde(t_conv)))
            write(dat2, '(2es14.6)') model%dx, pde_err(g)
        end do
        close(dat)
        close(dat2)
        pde_order = log(pde_err(size(grids) - 1) / pde_err(size(grids))) / &
                    log(dx_g(size(grids) - 1) / dx_g(size(grids)))
        print '(a10, f8.4, i6, f9.3, f12.4, f12.3)', trim(names(m)), h_star, tableau%order, observed, worst_ratio, pde_order
        write(table, '("\texttt{", a, "} & ", f6.4, " & ", i0, " & ", a, " & ", f6.4, " & ", f5.2, " \\")') &
            trim(names(m)), h_star, tableau%order, order_text(observed), worst_ratio, pde_order
    end do
    write(table, '(a)') "\bottomrule"
    write(table, '(a)') "\end{tabular}"
    close(table)

contains

    !> The most negative eigenvalue of the semi-discrete Laplacian on n nodes.
    pure real(wp) function eigenvalue(n) result(lam)
        integer, intent(in) :: n
        real(wp) :: h

        h = 1.0_wp / (n + 1)
        lam = -(4 / h**2) * sin(n * pi * h / 2)**2
    end function eigenvalue

    !> The exact solution of the semi-discrete system at the interior nodes.
    function semi_discrete(t) result(u)
        real(wp), intent(in) :: t
        real(wp) :: u(model%n), x, lam
        integer  :: j, k

        do j = 1, model%n
            x = j * model%dx
            u(j) = 0.0_wp
            do k = 1, size(modes)
                lam = -(4 / model%dx**2) * sin(modes(k) * pi * model%dx / 2)**2
                u(j) = u(j) + amplitudes(k) * exp(lam * t) * sin(modes(k) * pi * x)
            end do
        end do
    end function semi_discrete

    !> The exact solution of the PDE at the interior nodes.
    function pde(t) result(u)
        real(wp), intent(in) :: t
        real(wp) :: u(model%n), x
        integer  :: j, k

        do j = 1, model%n
            x = j * model%dx
            u(j) = sum(amplitudes * exp(-(modes * pi)**2 * t) * sin(modes * pi * x))
        end do
    end function pde

    !> Interior values at t_end from u(x, 0), steps of size dt (t_end/dt integer).
    function integrate(name, dt, t_end) result(u)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: dt, t_end
        real(wp) :: u(model%n)
        class(time_scheme_c), allocatable, target :: scheme
        type(heat_state_t), target :: state
        type(status_t) :: status
        real(wp) :: time
        integer  :: i

        call rk_scheme(scheme, name, dt, model)
        state = initial_state()
        time = 0.0_wp
        do i = 1, nint(t_end / dt)
            call scheme%advance(time, state, status)
            if (status%failed) error stop status%message
        end do
        u = state%u(1:model%n)
    end function integrate

    function initial_state() result(state)
        type(heat_state_t) :: state
        integer :: j

        allocate(state%u(0:model%n+1))
        state%u = [(sum(amplitudes * sin(modes * pi * j * model%dx)), j = 0, model%n + 1)]
    end function initial_state

    !> Largest stable step near `guess`: bisection between 0.8 and 1.2 times it,
    !> stable if after n_probe steps the solution has not grown (the heat
    !> equation only decays). A grid-scale perturbation of 1e-6 makes sure the
    !> fastest mode is present.
    real(wp) function stable_step(name, guess) result(h)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: guess
        real(wp) :: lo, hi
        integer  :: it

        lo = 0.8_wp * guess
        hi = 1.2_wp * guess
        do it = 1, 30
            h = (lo + hi) / 2
            if (stays_bounded(name, h)) then
                lo = h
            else
                hi = h
            end if
        end do
        h = lo
    end function stable_step

    logical function stays_bounded(name, dt)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: dt
        class(time_scheme_c), allocatable, target :: scheme
        type(heat_state_t), target :: state
        type(status_t) :: status
        real(wp) :: time, top
        integer  :: i, j

        call rk_scheme(scheme, name, dt, model)
        state = initial_state()
        state%u(1:model%n) = state%u(1:model%n) + 1.0e-6_wp * [((-1.0_wp)**j, j = 1, model%n)]
        top = maxval(abs(state%u))
        time = 0.0_wp
        do i = 1, n_probe
            call scheme%advance(time, state, status)
        end do
        ! all(...): a NaN fails the comparison. maxval would skip NaNs and see only
        ! the boundary nodes, which `apply_boundaries` keeps at 0, after a blow-up.
        stays_bounded = all(abs(state%u) <= top)
    end function stays_bounded

    !> h*: the first h > 0 with |R(-h)| > 1, from the tableau.
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

    subroutine write_figures()
        class(time_scheme_c), allocatable, target :: scheme
        type(heat_state_t), target :: state
        type(status_t) :: status
        real(wp) :: time, dt, times(5)
        integer  :: i, j, k, unit

        model = make_heat_model(n_conv)
        ! Space-time map of the rk4 solution, 0 <= t <= 0.01.
        dt = 0.01_wp / 1000
        call rk_scheme(scheme, "rk4", dt, model)
        state = initial_state()
        time = 0.0_wp
        open(newunit=unit, file=out // "spacetime.dat", status="replace", action="write")
        write(unit, '(a)') "x t u"
        do i = 0, 1000
            if (mod(i, 10) == 0) then
                do j = 0, model%n + 1
                    write(unit, '(3es14.6)') j * model%dx, time, state%u(j)
                end do
                write(unit, '(a)') ""
            end if
            if (i < 1000) call scheme%advance(time, state, status)
        end do
        close(unit)
        ! Profiles: rk4 at a few times against the PDE solution.
        times = [0.0_wp, 0.002_wp, 0.005_wp, 0.01_wp, 0.03_wp]
        open(newunit=unit, file=out // "profiles.dat", status="replace", action="write")
        write(unit, '(a)') "x u1 u2 u3 u4 u5 e1 e2 e3 e4 e5"
        block
            real(wp) :: num(0:model%n+1, 5), ex(0:model%n+1, 5)
            do k = 1, 5
                if (times(k) > 0) then
                    num(1:model%n, k) = integrate("rk4", 1.0e-5_wp, times(k))     ! well inside rk4's limit
                else
                    num(1:model%n, k) = semi_discrete(0.0_wp)
                end if
                num(0, k) = 0.0_wp
                num(model%n + 1, k) = 0.0_wp
                ex(1:model%n, k) = pde(times(k))
                ex(0, k) = 0.0_wp
                ex(model%n + 1, k) = 0.0_wp
            end do
            do j = 0, model%n + 1
                write(unit, '(11es14.6)') j * model%dx, num(j, :), ex(j, :)
            end do
        end block
        close(unit)
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
        character(:), allocatable :: list, glist

        list = ""
        do j = 1, size(names)
            if (j > 1) list = list // ", "
            list = list // "\texttt{" // trim(names(j)) // "}"
        end do
        glist = ""
        do j = 1, size(grids)
            if (j > 1) glist = glist // ", "
            glist = glist // itoa(grids(j))
        end do
        open(newunit=u, file=out // "parameters.tex", status="replace", action="write")
        write(u, '(a)') "% Written by examples/heat/main.f90: do not edit."
        write(u, '(a)') "\begin{tabular}{l>{\raggedright\arraybackslash}p{12.4cm}}"
        write(u, '(a)') "\toprule"
        write(u, '(a)') "initial condition & $u(x, 0) = \sin\pi x + \tfrac12\sin 4\pi x + \tfrac14\sin 9\pi x$ \\"
        write(u, '(a)') "space & central differences, $n$ interior nodes, $\Delta x = 1/(n + 1)$; boundary nodes" // &
                        " set to 0 by \texttt{apply\_boundaries} \\"
        write(u, '(a)') "time convergence & $n = " // itoa(n_conv) // "$, error against the semi-discrete solution at" // &
                        " $t = " // fmt(t_conv) // "$, $\Delta t = 10^{-4}, \dots, 10^{-4}/" // &
                        itoa(2**(n_levels - 1)) // "$; order from the two finest levels above $10^{" // &
                        itoa(nint(log10(floor))) // "}$, kept if within 0.3 of the next slope \\"
        write(u, '(a)') "stability & $n = " // glist // "$; largest $\Delta t$ for which " // itoa(n_probe) // &
                        " steps do not let $\max|u|$ grow (initial value plus $10^{-6}(-1)^j$); bisection in" // &
                        " $[0.8, 1.2]$ times the prediction; table: the ratio farthest from 1 \\"
        write(u, '(a)') "PDE error & at $t = " // fmt(t_conv) // "$ with $\Delta t \approx 0.9$ of each method's limit;" // &
                        " order from the two finest grids \\"
        write(u, '(a)') "figures & \texttt{rk4}, $n = " // itoa(n_conv) // "$ \\"
        write(u, '(a)') "methods & " // list // " \\"
        write(u, '(a)') "interface & model classes, \texttt{rk\_scheme}, \texttt{apply\_boundaries} for the boundary \\"
        write(u, '(a)') "\bottomrule"
        write(u, '(a)') "\end{tabular}"
        close(u)
    end subroutine write_parameters

    function fmt(x) result(s)
        real(wp), intent(in) :: x
        character(:), allocatable :: s
        character(len=24) :: buffer

        write(buffer, '(f10.4)') x
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
