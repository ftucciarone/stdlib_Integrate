program main
    !! Brusselator in one dimension (SciML prob_ode_brusselator_1d; parameters
    !! of Hairer & Wanner, Solving ODEs II, section IV.1), periodic on [0, 1):
    !!   u_t = A - (B + 1) u + u^2 v + alpha u_xx,
    !!   v_t = B u - u^2 v + alpha v_xx,
    !!   A = 1, B = 3, alpha = 1/50,  u(x, 0) = 1 + sin(2 pi x), v(x, 0) = 3.
    !! Without diffusion the reaction has the fixed point (A, B/A), unstable for
    !! B > 1 + A^2, and oscillates on a limit cycle. Central differences on n
    !! points make the diffusion eigenvalues -(4 alpha/dx^2) sin^2(k pi dx), the
    !! largest in size 4 alpha/dx^2. Tests:
    !!   - convergence: the error at t_conv (n = n_conv) as dt halves, against a
    !!     dp54 reference;
    !!   - stability against the grid: the largest stable dt against the
    !!     diffusion prediction h* dx^2/(4 alpha), which the reaction perturbs on
    !!     coarse grids only.
    !! Uses the plain-array schemes, state (u(1:n), v(1:n)). Writes
    !! examples/brusselator/data/: parameters.tex, results.tex, spacetime.dat
    !! (x, t, u, v), probe.dat (t, u, v at x = 1/2), error_<method>.dat,
    !! stable_<method>.dat, read by summary.tex. Run from the project root.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: rk_scheme, stage_t, rk_array_scheme_t, rk_method_names, rk_tableau, rk_tableau_t
    use brusselator_stability, only: stability_r
    implicit none

    character(*), parameter :: out = "examples/brusselator/data/"
    real(wp),     parameter :: pi = acos(-1.0_wp)
    real(wp),     parameter :: a = 1.0_wp, b = 3.0_wp, alpha = 1.0_wp / 50
    integer,      parameter :: n_conv = 64                     !! grid of the convergence study and figures
    real(wp),     parameter :: t_conv = 10.0_wp
    real(wp),     parameter :: dt_coarse = 0.004_wp            !! stable for every method at n_conv
    integer,      parameter :: n_levels = 7
    real(wp),     parameter :: dt_reference = 1.0e-5_wp
    integer,      parameter :: grids(6) = [16, 32, 64, 128, 256, 512]
    real(wp),     parameter :: t_probe = 2.0_wp                !! stability probes run to this time
    real(wp),     parameter :: t_figure = 20.0_wp
    real(wp),     parameter :: floor = 1.0e-12_wp

    character(len=8), allocatable :: names(:)
    type(rk_tableau_t) :: tableau
    integer  :: n                                              !! current grid
    real(wp) :: dx
    real(wp), allocatable :: reference(:)
    real(wp) :: dt(n_levels), err(n_levels), observed, h_star, predicted, measured, ratio_fine, ratio_coarse
    integer  :: m, l, g, table, dat

    names = rk_method_names()
    call execute_command_line("mkdir -p " // out)
    call set_grid(n_conv)
    reference = integrate("dp54", dt_reference, t_conv)
    call write_parameters()
    call write_figures()

    open(newunit=table, file=out // "results.tex", status="replace", action="write")
    write(table, '(a)') "% Written by examples/brusselator/main.f90: do not edit."
    write(table, '(a)') "\begin{tabular}{lccccc}"
    write(table, '(a)') "\toprule"
    write(table, '(a)') "& & & & \multicolumn{2}{c}{stable $\Delta t$ / $h^*\Delta x^2/(4\alpha)$} \\"
    write(table, '(a)') "\cmidrule(lr){5-6}"
    write(table, '(a)') "name & $h^*$ & order & observed & $n = " // itoa(grids(1)) // "$ & $n = " // &
                        itoa(grids(size(grids))) // "$ \\"
    write(table, '(a)') "\midrule"
    print '(a10, a8, a6, a9, 2a10)', "name", "h*", "order", "obs.", "ratio n16", "ratio 512"
    do m = 1, size(names)
        tableau = rk_tableau(trim(names(m)))
        h_star = real_axis_limit(tableau)
        call set_grid(n_conv)
        open(newunit=dat, file=out // "error_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat, '(a)') "dt error"
        do l = 1, n_levels
            dt(l) = dt_coarse / 2**(l - 1)
            err(l) = maxval(abs(integrate(trim(names(m)), dt(l), t_conv) - reference))
            if (err(l) > 1.0e-16_wp .and. err(l) < 1.0e3_wp) write(dat, '(2es14.6)') dt(l), err(l)
        end do
        close(dat)
        observed = slope(dt, err)
        open(newunit=dat, file=out // "stable_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat, '(a)') "dx measured predicted"
        do g = 1, size(grids)
            call set_grid(grids(g))
            predicted = h_star * dx**2 / (4 * alpha)
            measured = stable_step(trim(names(m)), predicted)
            write(dat, '(3es14.6)') dx, measured, predicted
            if (g == 1) ratio_coarse = measured / predicted
            if (g == size(grids)) ratio_fine = measured / predicted
        end do
        close(dat)
        print '(a10, f8.4, i6, f9.3, 2f10.4)', trim(names(m)), h_star, tableau%order, observed, ratio_coarse, ratio_fine
        write(table, '("\texttt{", a, "} & ", f6.4, " & ", i0, " & ", a, " & ", f6.4, " & ", f6.4, " \\")') &
            trim(names(m)), h_star, tableau%order, order_text(observed), ratio_coarse, ratio_fine
    end do
    write(table, '(a)') "\bottomrule"
    write(table, '(a)') "\end{tabular}"
    close(table)

contains

    subroutine set_grid(points)
        integer, intent(in) :: points

        n = points
        dx = 1.0_wp / n
    end subroutine set_grid

    !> u(1:n) = y(1:n), v(1:n) = y(n+1:2n), periodic central differences.
    subroutine brusselator(t, y, dydt, stage)
        real(wp), intent(in)  :: t
        real(wp), intent(in)  :: y(:)
        real(wp), intent(out) :: dydt(:)
        type(stage_t), intent(in) :: stage
        real(wp) :: c
        integer  :: i, im, ip

        c = alpha / dx**2
        do i = 1, n
            im = modulo(i - 2, n) + 1
            ip = modulo(i, n) + 1
            dydt(i)     = a - (b + 1) * y(i) + y(i)**2 * y(n+i) + c * (y(im) - 2 * y(i) + y(ip))
            dydt(n + i) = b * y(i) - y(i)**2 * y(n+i) + c * (y(n+im) - 2 * y(n+i) + y(n+ip))
        end do
    end subroutine brusselator

    function initial() result(y)
        real(wp) :: y(2 * n)
        integer  :: i

        y(1:n) = [(1 + sin(2 * pi * (i - 1) * dx), i = 1, n)]
        y(n+1:2*n) = 3.0_wp
    end function initial

    function integrate(name, step, t_end) result(y)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: step, t_end
        real(wp) :: y(2 * n), time
        type(rk_array_scheme_t) :: scheme
        integer :: i

        call rk_scheme(scheme, name, 2 * n, step)
        y = initial()
        time = 0.0_wp
        do i = 1, nint(t_end / step)
            call scheme%advance(brusselator, time, y)
        end do
    end function integrate

    !> Largest step for which the run (to t_probe, and at least 2000 steps)
    !> stays bounded (all values finite and below 100): bisection between 0.3
    !> and 2 times `guess`.
    real(wp) function stable_step(name, guess) result(h)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: guess
        real(wp) :: lo, hi
        integer  :: it

        lo = 0.3_wp * guess
        hi = 2.0_wp * guess
        do it = 1, 30
            h = (lo + hi) / 2
            if (bounded(name, h)) then
                lo = h
            else
                hi = h
            end if
        end do
        h = lo
    end function stable_step

    logical function bounded(name, step)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: step
        type(rk_array_scheme_t) :: scheme
        real(wp) :: y(2 * n), time
        integer  :: i, j

        call rk_scheme(scheme, name, 2 * n, step)
        y = initial()
        y(1:n) = y(1:n) + 1.0e-6_wp * [((-1.0_wp)**j, j = 1, n)]     ! seed the fastest diffusion mode
        time = 0.0_wp
        bounded = .true.
        do i = 1, max(ceiling(t_probe / step), 2000)               ! enough steps for a slow growth to show
            call scheme%advance(brusselator, time, y)
            if (.not. all(abs(y) < 100)) then                         ! fails on NaN too
                bounded = .false.
                return
            end if
        end do
    end function bounded

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
        type(rk_array_scheme_t) :: scheme
        real(wp) :: y(2 * n), time, step
        integer  :: i, j, unit, up

        step = 1.0e-3_wp
        call rk_scheme(scheme, "rk4", 2 * n, step)
        y = initial()
        time = 0.0_wp
        open(newunit=unit, file=out // "spacetime.dat", status="replace", action="write")
        open(newunit=up, file=out // "probe.dat", status="replace", action="write")
        write(unit, '(a)') "x t u v"
        write(up, '(a)') "t u v"
        do i = 0, nint(t_figure / step)
            if (mod(i, 100) == 0) then
                do j = 1, n
                    write(unit, '(4es14.6)') (j - 1) * dx, time, y(j), y(n + j)
                end do
                write(unit, '(4es14.6)') 1.0_wp, time, y(1), y(n + 1)            ! close the period
                write(unit, '(a)') ""
            end if
            if (mod(i, 10) == 0) write(up, '(3es14.6)') time, y(n / 2 + 1), y(n + n / 2 + 1)
            if (i < nint(t_figure / step)) call scheme%advance(brusselator, time, y)
        end do
        close(unit)
        close(up)
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
        write(u, '(a)') "% Written by examples/brusselator/main.f90: do not edit."
        write(u, '(a)') "\begin{tabular}{l>{\raggedright\arraybackslash}p{12.4cm}}"
        write(u, '(a)') "\toprule"
        write(u, '(a)') "parameters & $A = 1$, $B = 3$, $\alpha = 1/50$ \\"
        write(u, '(a)') "initial condition & $u = 1 + \sin 2\pi x$, $v = 3$; periodic on $[0, 1)$ \\"
        write(u, '(a)') "space & central differences, $n$ points, $\Delta x = 1/n$; $2n$ unknowns \\"
        write(u, '(a)') "convergence & $n = " // itoa(n_conv) // "$, error at $t = " // fmt(t_conv) // &
                        "$ against \texttt{dp54} with $\Delta t = 10^{-5}$; $\Delta t = " // fmt(dt_coarse) // &
                        ", \dots, " // fmt(dt_coarse) // "/" // itoa(2**(n_levels - 1)) // "$; order from the two" // &
                        " finest levels above $10^{" // itoa(nint(log10(floor))) // "}$, kept if within 0.3 of the next \\"
        write(u, '(a)') "stability & $n = " // glist // "$; largest $\Delta t$ keeping every value finite and below" // &
                        " 100 up to $t = " // fmt(t_probe) // "$ and for at least 2000 steps (initial value plus" // &
                        " $10^{-6}(-1)^j$ in $u$); bisection" // &
                        " in $[0.3, 2]$ times $h^*\Delta x^2/(4\alpha)$ \\"
        write(u, '(a)') "figures & \texttt{rk4}, $n = " // itoa(n_conv) // "$, $\Delta t = 10^{-3}$, $0 \le t \le " // &
                        fmt(t_figure) // "$ \\"
        write(u, '(a)') "methods & " // list // " \\"
        write(u, '(a)') "interface & plain arrays, \texttt{rk\_array\_scheme\_t} \\"
        write(u, '(a)') "\bottomrule"
        write(u, '(a)') "\end{tabular}"
        close(u)
    end subroutine write_parameters

    function fmt(x) result(s)
        real(wp), intent(in) :: x
        character(:), allocatable :: s
        character(len=24) :: buffer

        write(buffer, '(f12.4)') x
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
