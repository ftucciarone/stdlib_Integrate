program main
    !! Nonlinear chemistry (SciML prob_ode_nonlinchem):
    !!   y1' = -y1,  y2' = y1 - y2^2,  y3' = y2^2,  y(0) = (1, 0, 0),  0 <= t <= 20,
    !! with an exact solution (nonlinchem_exact). Two tests:
    !!   - convergence: the error at t = 20 as dt halves. The problem is
    !!     nonlinear (y2^2), so the order that counts is the general one: rk3ws,
    !!     third order on linear problems, drops to 2 here;
    !!   - invariant: y1 + y2 + y3 = 1 for all t. Every Runge-Kutta method keeps
    !!     linear invariants exactly, so the sum must stay 1 to round-off.
    !! Uses the plain-array schemes. Writes examples/nonlinchem/data/:
    !! parameters.tex, results.tex, exact.dat and solution_<method>.dat (t, y1,
    !! y2, y3 at dt_plot), history_<method>.dat (t, error at dt_history),
    !! error_<method>.dat (dt, error at t_end), read by summary.tex. Run from the
    !! project root.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: rk_scheme, stage_t, rk_array_scheme_t, rk_method_names, rk_tableau, rk_tableau_t
    use nonlinchem_exact, only: exact_solution, bessel_i0, bessel_i1, bessel_k0, bessel_k1
    implicit none

    character(*), parameter :: out = "examples/nonlinchem/data/"
    real(wp),     parameter :: t_end = 20.0_wp
    real(wp),     parameter :: dt_coarse = 0.5_wp     !! coarsest step of the convergence study
    integer,      parameter :: n_levels = 8           !! dt = 0.5, 0.25, ..., 0.5/128
    real(wp),     parameter :: round_off = 1.0e-12_wp     !! below this, round-off pollutes the slope
    real(wp),     parameter :: dt_plot = 0.5_wp       !! step of the solution figure
    real(wp),     parameter :: dt_history = 0.25_wp   !! step of the error-history figure

    character(len=8), allocatable :: names(:)
    type(rk_tableau_t) :: tableau
    real(wp) :: dt(n_levels), err(n_levels), drift(n_levels), observed, check, worst_drift
    integer  :: m, l, i, table, dat

    names = rk_method_names()
    call execute_command_line("mkdir -p " // out)
    check = exact_check()
    print '(a, es9.2)', "exact solution against mpmath: largest relative difference ", check
    call write_parameters(check)

    ! The exact solution, the solutions at dt_plot, the error histories at dt_history.
    open(newunit=dat, file=out // "exact.dat", status="replace", action="write")
    write(dat, '(a)') "t y1 y2 y3"
    do i = 0, 400
        write(dat, '(4es16.8)') t_end * i / 400, exact_solution(t_end * i / 400)
    end do
    close(dat)
    do m = 1, size(names)
        call run(trim(names(m)), dt_plot, observed, worst_drift, "solution_" // trim(names(m)) // ".dat", .false.)
        call run(trim(names(m)), dt_history, observed, worst_drift, "history_" // trim(names(m)) // ".dat", .true.)
    end do

    ! Convergence and invariant.
    open(newunit=table, file=out // "results.tex", status="replace", action="write")
    write(table, '(a)') "% Written by examples/nonlinchem/main.f90: do not edit."
    write(table, '(a)') "\begin{tabular}{llccccc}"
    write(table, '(a)') "\toprule"
    write(table, '(a)') "method & name & order & linear order & observed & $\max|y_1 + y_2 + y_3 - 1|$ \\"
    write(table, '(a)') "\midrule"
    print '(a20, a10, 3a8, a14)', "method", "name", "order", "linear", "obs.", "max drift"
    do m = 1, size(names)
        tableau = rk_tableau(trim(names(m)))
        open(newunit=dat, file=out // "error_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat, '(a)') "dt error"
        do l = 1, n_levels
            dt(l) = dt_coarse / 2**(l - 1)
            call run(trim(names(m)), dt(l), err(l), drift(l))
            if (err(l) > 1.0e-15_wp) write(dat, '(2es14.6)') dt(l), err(l)
        end do
        close(dat)
        observed = slope(dt, err)
        print '(a20, a10, 2i8, f8.3, es14.2)', tableau%title, trim(names(m)), tableau%order, tableau%linear_order, &
            observed, maxval(drift)
        write(table, '(a, " & \texttt{", a, "} & ", i0, " & ", i0, " & ", f5.2, " & ", a, " \\")') &
            tableau%title, trim(names(m)), tableau%order, tableau%linear_order, observed, sci(maxval(drift))
    end do
    ! rk3ws with pre- and post-processing (effective order 3, see the summary):
    ! start from u0 - (h^2/24) f'f(u0), correct the result by + (h^2/24) f'f(u).
    open(newunit=dat, file=out // "error_rk3ws_processed.dat", status="replace", action="write")
    write(dat, '(a)') "dt error"
    do l = 1, n_levels
        dt(l) = dt_coarse / 2**(l - 1)
        call run("rk3ws", dt(l), err(l), drift(l), processed=.true.)
        if (err(l) > 1.0e-15_wp) write(dat, '(2es14.6)') dt(l), err(l)
    end do
    close(dat)
    observed = slope(dt, err)
    print '(a20, a10, 2i8, f8.3, es14.2)', "rk3ws, processed", "rk3ws", 2, 3, observed, maxval(drift)
    write(table, '(a)') "\midrule"
    write(table, '("\texttt{rk3ws}, pre- and post-processed & & & & ", f5.2, " & ", a, " \\")') &
        observed, sci(maxval(drift))
    write(table, '(a)') "\bottomrule"
    write(table, '(a)') "\end{tabular}"
    close(table)

contains

    subroutine chemistry(t, y, dydt, stage)
        real(wp), intent(in)  :: t
        real(wp), intent(in)  :: y(:)
        real(wp), intent(out) :: dydt(:)
        type(stage_t), intent(in) :: stage

        dydt(1) = -y(1)
        dydt(2) = y(1) - y(2)**2
        dydt(3) = y(2)**2
    end subroutine chemistry

    !> Integrate to t_end with step dt: `err` = max-norm error at t_end,
    !> `drift` = largest |y1 + y2 + y3 - 1| over the run. With `file`, write
    !> every step: (t, y1, y2, y3), or (t, max-norm error) if `errors`.
    !> f'(y) f(y) for the chemistry right-hand side: u'' along a solution.
    pure function second_derivative(y) result(a)
        real(wp), intent(in) :: y(3)
        real(wp) :: a(3)

        a(1) = y(1)
        a(2) = -y(1) - 2 * y(2) * (y(1) - y(2)**2)
        a(3) = 2 * y(2) * (y(1) - y(2)**2)
    end function second_derivative

    subroutine run(name, dt, err, drift, file, errors, processed)
        character(*), intent(in)  :: name
        real(wp),     intent(in)  :: dt
        real(wp),     intent(out) :: err, drift
        character(*), intent(in), optional :: file
        logical,      intent(in), optional :: errors
        logical,      intent(in), optional :: processed   !! pre- and post-process (rk3ws)

        type(rk_array_scheme_t) :: scheme
        real(wp) :: y(3), time
        integer  :: i, n, u

        n = nint(t_end / dt)
        call rk_scheme(scheme, name, 3, dt)
        y = [1.0_wp, 0.0_wp, 0.0_wp]
        if (present(processed)) y = y - dt**2 / 24 * second_derivative(y)
        time = 0.0_wp
        drift = 0.0_wp
        if (present(file)) then
            open(newunit=u, file=out // file, status="replace", action="write")
            if (errors) then
                write(u, '(a)') "t error"
            else
                write(u, '(a)') "t y1 y2 y3"
                write(u, '(4es16.8)') time, y
            end if
        end if
        do i = 1, n
            call scheme%advance(chemistry, time, y)
            drift = max(drift, abs(sum(y) - 1))
            if (present(file)) then
                if (errors) then
                    write(u, '(2es14.6)') time, max(maxval(abs(y - exact_solution(time))), 1.0e-16_wp)
                else
                    write(u, '(4es16.8)') time, y
                end if
            end if
        end do
        if (present(file)) close(u)
        if (present(processed)) y = y + dt**2 / 24 * second_derivative(y)
        err = maxval(abs(y - exact_solution(time)))
    end subroutine run

    !> Largest relative difference between this program's exact solution and
    !> Bessel functions and reference values computed with mpmath (30 digits).
    real(wp) function exact_check() result(worst)
        real(wp), parameter :: x(3) = [2.0_wp, 1.0_wp, 0.1_wp]
        real(wp), parameter :: i0(3) = [2.2795853023360673_wp, 1.2660658777520083_wp, 1.0025015629340956_wp]
        real(wp), parameter :: i1(3) = [1.5906368546373291_wp, 0.56515910399248503_wp, 0.050062526047092695_wp]
        real(wp), parameter :: k0(3) = [0.11389387274953344_wp, 0.42102443824070833_wp, 2.4270690247020166_wp]
        real(wp), parameter :: k1(3) = [0.13986588181652243_wp, 0.60190723019723457_wp, 9.8538447808706056_wp]
        real(wp), parameter :: t(3) = [2.0_wp, 10.0_wp, 20.0_wp]
        real(wp), parameter :: y(3, 3) = reshape([ &
            0.13533528323661269_wp, 0.47842176645076096_wp, 0.38624295031262635_wp, &
            4.5399929762484852e-5_wp, 0.11079059098117565_wp, 0.88916400908906187_wp, &
            2.0611536224385578e-9_wp, 0.052572280220485125_wp, 0.94742771771836125_wp], [3, 3])
        integer :: k

        worst = 0.0_wp
        do k = 1, 3
            worst = max(worst, abs(bessel_i0(x(k)) / i0(k) - 1), abs(bessel_i1(x(k)) / i1(k) - 1), &
                        abs(bessel_k0(x(k)) / k0(k) - 1), abs(bessel_k1(x(k)) / k1(k) - 1))
            worst = max(worst, maxval(abs(exact_solution(t(k)) / y(:, k) - 1)))
        end do
    end function exact_check

    !> Slope of log(err) against log(dt) between the two finest levels whose
    !> errors are both above the round-off floor.
    real(wp) function slope(dt, err) result(p)
        real(wp), intent(in) :: dt(:), err(:)
        integer :: l

        p = 0.0_wp
        do l = size(dt), 2, -1
            if (err(l) > round_off .and. err(l - 1) > round_off) then
                p = log(err(l - 1) / err(l)) / log(dt(l - 1) / dt(l))
                return
            end if
        end do
    end function slope

    !> The parameters of this run, as a LaTeX table.
    subroutine write_parameters(check)
        real(wp), intent(in) :: check
        integer :: u, k
        character(:), allocatable :: list

        list = ""
        do k = 1, size(names)
            if (k > 1) list = list // ", "
            list = list // "\texttt{" // trim(names(k)) // "}"
        end do
        open(newunit=u, file=out // "parameters.tex", status="replace", action="write")
        write(u, '(a)') "% Written by examples/nonlinchem/main.f90: do not edit."
        write(u, '(a)') "\begin{tabular}{ll}"
        write(u, '(a)') "\toprule"
        write(u, '(a)') "initial condition & $y(0) = (1, 0, 0)$ \\"
        write(u, '(a)') "interval & $0 \le t \le " // fmt(t_end) // "$ \\"
        write(u, '(a)') "steps, convergence & $\Delta t = " // fmt(dt_coarse) // ", " // fmt(dt_coarse / 2) // &
                        ", \dots, " // fmt(dt_coarse) // "/" // itoa(2**(n_levels - 1)) // "$ (" // itoa(n_levels) // &
                        " levels) \\"
        write(u, '(a)') "error & $\max_i |y_i - y_i^{\mathrm{exact}}|$ at $t = " // fmt(t_end) // "$ \\"
        write(u, '(a)') "observed order & slope of the two finest levels with errors above $10^{" // &
                        itoa(nint(log10(round_off))) // "}$ \\"
        write(u, '(a)') "invariant & $\max_n |y_1 + y_2 + y_3 - 1|$ over every run of the study \\"
        write(u, '(a)') "figures & solutions at $\Delta t = " // fmt(dt_plot) // "$; error histories at $\Delta t = " // &
                        fmt(dt_history) // "$ \\"
        write(u, '(a)') "methods & " // list // " \\"
        write(u, '(a)') "interface & plain arrays, \texttt{rk\_array\_scheme\_t} \\"
        write(u, '(a)') "exact solution & Bessel series (30 terms), within " // sci(check) // &
                        " (relative) of mpmath at 30 digits \\"
        write(u, '(a)') "\bottomrule"
        write(u, '(a)') "\end{tabular}"
        close(u)
    end subroutine write_parameters

    !> A number in LaTeX, e.g. $0.5$ or $2.3\times10^{-16}$.
    function sci(x) result(s)
        real(wp), intent(in) :: x
        character(:), allocatable :: s
        character(len=32) :: mantissa
        integer :: e

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
        character(len=16) :: buffer

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
