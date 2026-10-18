program main
    !! Stability limit on the negative real axis. On y' = -y a step multiplies y
    !! by R(-h), h = dt, so the run stays bounded exactly when |R(-h)| <= 1. For
    !! every method:
    !!   - theory: the largest h* with |R(-h)| <= 1 on (0, h*], R from the
    !!     tableau (stability_function, independent of the step routines);
    !!   - measured: the largest dt for which 2000 steps of the array scheme
    !!     keep |y| <= 1, found by bisection;
    !!   - the growth per step |y_N / y_0|^(1/N) measured at h = 0.05 ... 4,
    !!     against |R(-h)|.
    !! Writes examples/stability/data/: parameters.tex, results.tex,
    !! exact.dat and solution_<method>.dat (t, y at h_plot), edge_below.dat and
    !! edge_above.dat (t, |y| for rk4 just inside and outside its limit),
    !! theory_<method>.dat (h, |R|) and measured_<method>.dat (h, growth), read
    !! by summary.tex. Run from the project root.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: rk_scheme, stage_t, rk_array_scheme_t, rk_method_names, rk_tableau, rk_tableau_t
    use stability_polynomial, only: stability_r
    implicit none

    character(*), parameter :: out = "examples/stability/data/"
    integer,      parameter :: n_limit = 2000      !! steps of a run in the bisection
    integer,      parameter :: n_growth = 50       !! steps of a run measuring growth
    real(wp),     parameter :: h_max = 4.0_wp
    real(wp),     parameter :: h_plot = 1.5_wp     !! step of the solution figure
    integer,      parameter :: n_plot = 10         !! steps of the solution figure
    integer,      parameter :: n_edge = 300        !! steps of the runs at the edge
    real(wp),     parameter :: edge = 0.02_wp      !! edge runs at (1 -+ edge) h*

    character(len=8), allocatable :: names(:)
    type(rk_tableau_t) :: tableau
    real(wp) :: h_theory, h_measured, h
    integer  :: m, i, table, dat

    names = rk_method_names()
    call execute_command_line("mkdir -p " // out)
    call write_parameters()

    ! The solutions at h_plot against the exact one, and rk4 at its edge.
    open(newunit=dat, file=out // "exact.dat", status="replace", action="write")
    write(dat, '(a)') "t y"
    do i = 0, 300
        write(dat, '(2es14.6)') n_plot * h_plot * i / 300, exp(-n_plot * h_plot * i / 300)
    end do
    close(dat)
    do m = 1, size(names)
        call trajectory(trim(names(m)), h_plot, n_plot, "solution_" // trim(names(m)) // ".dat", .false.)
    end do
    h = limit_from_r(rk_tableau("rk4"))
    call trajectory("rk4", (1 - edge) * h, n_edge, "edge_below.dat", .true.)
    call trajectory("rk4", (1 + edge) * h, n_edge, "edge_above.dat", .true.)

    open(newunit=table, file=out // "results.tex", status="replace", action="write")
    write(table, '(a)') "% Written by examples/stability/main.f90: do not edit."
    write(table, '(a)') "\begin{tabular}{llcc}"
    write(table, '(a)') "\toprule"
    write(table, '(a)') "method & name & $h^*$ from $R$ & $h^*$ measured \\"
    write(table, '(a)') "\midrule"
    print '(a20, a10, 2a14)', "method", "name", "h* theory", "h* measured"

    do m = 1, size(names)
        tableau = rk_tableau(trim(names(m)))
        h_theory = limit_from_r(tableau)
        h_measured = limit_from_runs(trim(names(m)))
        print '(a20, a10, 2f14.5)', tableau%title, trim(names(m)), h_theory, h_measured
        write(table, '(a, " & \texttt{", a, "} & ", f7.4, " & ", f7.4, " \\")') &
            tableau%title, trim(names(m)), h_theory, h_measured

        open(newunit=dat, file=out // "theory_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat, '(a)') "h absR"
        do i = 0, 400
            h = h_max * i / 400
            write(dat, '(2es14.6)') h, abs(stability_r(tableau, cmplx(-h, 0.0_wp, wp)))
        end do
        close(dat)

        open(newunit=dat, file=out // "measured_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat, '(a)') "h growth"
        do i = 1, 80
            h = h_max * i / 80
            write(dat, '(2es14.6)') h, growth(trim(names(m)), h)
        end do
        close(dat)
    end do

    write(table, '(a)') "\bottomrule"
    write(table, '(a)') "\end{tabular}"
    close(table)

contains

    !> y' = -y
    subroutine decay(t, y, dydt, stage)
        real(wp), intent(in)  :: t
        real(wp), intent(in)  :: y(:)
        real(wp), intent(out) :: dydt(:)
        type(stage_t), intent(in) :: stage

        dydt = -y
    end subroutine decay

    !> |y| after n steps of size dt from y = 1.
    real(wp) function run(name, dt, n) result(y_end)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: dt
        integer,      intent(in) :: n

        type(rk_array_scheme_t) :: scheme
        real(wp) :: y(1), time
        integer  :: i

        call rk_scheme(scheme, name, 1, dt)
        y = 1.0_wp
        time = 0.0_wp
        do i = 1, n
            call scheme%advance(decay, time, y)
        end do
        y_end = abs(y(1))
    end function run

    !> Write (t, y), or (t, |y|), of n steps of size dt from y = 1 to data/<file>.
    subroutine trajectory(name, dt, n, file, absolute)
        character(*), intent(in) :: name, file
        real(wp),     intent(in) :: dt
        integer,      intent(in) :: n
        logical,      intent(in) :: absolute

        type(rk_array_scheme_t) :: scheme
        real(wp) :: y(1), time
        integer  :: i, u

        open(newunit=u, file=out // file, status="replace", action="write")
        write(u, '(a)') "t y"
        call rk_scheme(scheme, name, 1, dt)
        y = 1.0_wp
        time = 0.0_wp
        write(u, '(2es14.6)') time, y(1)
        do i = 1, n
            call scheme%advance(decay, time, y)
            write(u, '(2es14.6)') time, merge(abs(y(1)), y(1), absolute)
        end do
        close(u)
    end subroutine trajectory

    !> The parameters of this run, as a LaTeX table.
    subroutine write_parameters()
        integer :: u, k
        character(:), allocatable :: list

        list = ""
        do k = 1, size(names)
            if (k > 1) list = list // ", "
            list = list // "\texttt{" // trim(names(k)) // "}"
        end do
        open(newunit=u, file=out // "parameters.tex", status="replace", action="write")
        write(u, '(a)') "% Written by examples/stability/main.f90: do not edit."
        write(u, '(a)') "\begin{tabular}{ll}"
        write(u, '(a)') "\toprule"
        write(u, '(a)') "problem & $y' = -y$, $y(0) = 1$ \\"
        write(u, '(a)') "methods & " // list // " \\"
        write(u, '(a)') "interface & plain arrays, \texttt{rk\_array\_scheme\_t} \\"
        write(u, '(a)') "solution figure & $h = " // fmt(h_plot) // "$, " // itoa(n_plot) // " steps \\"
        write(u, '(a)') "edge runs & \texttt{rk4} at $(1 \mp " // fmt(edge) // ")\,h^*$, " // itoa(n_edge) // &
                        " steps \\"
        write(u, '(a)') "growth per step & $|y_N/y_0|^{1/N}$, $N = " // itoa(n_growth) // "$, at $h = " // &
                        fmt(h_max / 80) // ", " // fmt(2 * h_max / 80) // ", \dots, " // fmt(h_max) // "$ \\"
        write(u, '(a)') "$h^*$ from $R$ & first $h$ with $|R(-h)| > 1$: scan by $10^{-3}$, then 60 bisections \\"
        write(u, '(a)') "$h^*$ measured & largest $h$ with $|y_N| \le 1$ after $N = " // itoa(n_limit) // &
                        "$ steps: 50 bisections on $[" // fmt(0.5_wp) // ", " // fmt(h_max) // "]$ \\"
        write(u, '(a)') "\bottomrule"
        write(u, '(a)') "\end{tabular}"
        close(u)
    end subroutine write_parameters

    function fmt(x) result(s)
        real(wp), intent(in) :: x
        character(:), allocatable :: s
        character(len=16) :: buffer

        write(buffer, '(f8.3)') x
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

    !> Growth per step, |y_N / y_0|^(1/N), at step h.
    real(wp) function growth(name, h) result(g)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: h

        g = run(name, h, n_growth)**(1.0_wp / n_growth)
    end function growth

    !> Largest dt keeping 2000 steps bounded, by bisection (stable at 0.5,
    !> unstable at h_max for every method of the table).
    real(wp) function limit_from_runs(name) result(h)
        character(*), intent(in) :: name

        real(wp) :: lo, hi
        integer  :: it

        lo = 0.5_wp
        hi = h_max
        do it = 1, 50
            h = (lo + hi) / 2
            if (run(name, h, n_limit) <= 1.0_wp) then
                lo = h
            else
                hi = h
            end if
        end do
        h = lo
    end function limit_from_runs

    !> First h where |R(-h)| exceeds 1, scanning from 0, then bisected.
    real(wp) function limit_from_r(tableau) result(h)
        type(rk_tableau_t), intent(in) :: tableau

        real(wp) :: lo, hi
        integer  :: it

        lo = 0.0_wp
        do
            hi = lo + 1.0e-3_wp
            if (unstable(tableau, hi)) exit
            lo = hi
        end do
        do it = 1, 60
            h = (lo + hi) / 2
            if (unstable(tableau, h)) then
                hi = h
            else
                lo = h
            end if
        end do
        h = lo
    end function limit_from_r

    logical function unstable(tableau, h)
        type(rk_tableau_t), intent(in) :: tableau
        real(wp),           intent(in) :: h

        unstable = abs(stability_r(tableau, cmplx(-h, 0.0_wp, wp))) > 1.0_wp
    end function unstable

end program main
