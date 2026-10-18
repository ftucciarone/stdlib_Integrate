program main
    !! Harmonic oscillator x'' = -x, as x' = v, v' = -x, from (x, v) = (1, 0):
    !! exact solution x = cos t, v = -sin t, energy E = (x^2 + v^2)/2 constant.
    !! The system is y' = J y with eigenvalues +-i, so a step of size h multiplies
    !! the energy by |R(ih)|^2 and turns the phase by arg R(ih) instead of h. For
    !! every method:
    !!   - energy after n steps against |R(ih)|^(2n), over 100 periods at h = 0.1;
    !!   - phase error after 100 periods against 100*2*pi/h * (arg R(ih) - h);
    !!   - stability limit on the imaginary axis: the largest h with
    !!     |R(iy)| <= 1 on (0, h] from the tableau, against the largest h for
    !!     which 2000 steps do not increase the energy (bisection).
    !! Writes examples/oscillator/data/: parameters.tex, results.tex,
    !! tail_<method>.dat (t, x over the last n_tail_periods periods) and
    !! exact_tail.dat, portrait_<method>.dat (x, v every portrait_every steps),
    !! energy_<method>.dat (t, measured |E/E0 - 1|) and theory_<method>.dat
    !! (t, predicted), read by summary.tex. Run from the project root.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: rk_scheme, stage_t, rk_array_scheme_t, rk_method_names, rk_tableau, rk_tableau_t
    use oscillator_stability, only: stability_r
    implicit none

    character(*), parameter :: out = "examples/oscillator/data/"
    real(wp),     parameter :: pi = acos(-1.0_wp)
    real(wp),     parameter :: h = 0.1_wp
    integer,      parameter :: n_periods = 100
    integer,      parameter :: n_limit = 2000
    integer,      parameter :: n_tail_periods = 3   !! periods shown in the x(t) figure
    integer,      parameter :: portrait_every = 4   !! steps between points of the phase portrait

    character(len=8), allocatable :: names(:)
    type(rk_tableau_t) :: tableau
    complex(wp) :: r
    real(wp) :: y(2), time, e_ratio, e_predicted, phase, phase_predicted, lim_theory, lim_measured
    integer  :: m, n, n_steps, next_sample, table, dat_e, dat_t, dat_x, dat_p
    type(rk_array_scheme_t) :: scheme

    names = rk_method_names()
    n_steps = nint(n_periods * 2 * pi / h)
    call execute_command_line("mkdir -p " // out)
    call write_parameters()
    open(newunit=dat_x, file=out // "exact_tail.dat", status="replace", action="write")
    write(dat_x, '(a)') "t x"
    do n = 0, 600
        time = n_steps * h - n_tail_periods * 2 * pi * (1 - n / 600.0_wp)
        write(dat_x, '(2es14.6)') time, cos(time)
    end do
    close(dat_x)

    open(newunit=table, file=out // "results.tex", status="replace", action="write")
    write(table, '(a)') "% Written by examples/oscillator/main.f90: do not edit."
    write(table, '(a)') "\begin{tabular}{lcccccc}"
    write(table, '(a)') "\toprule"
    write(table, '(a)') "& \multicolumn{2}{c}{$E/E_0$ after 100 periods} & \multicolumn{2}{c}{phase error (rad)}" // &
                        " & \multicolumn{2}{c}{imaginary-axis limit} \\"
    write(table, '(a)') "\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}"
    write(table, '(a)') "name & predicted & measured & predicted & measured & from $R$ & measured \\"
    write(table, '(a)') "\midrule"
    print '(a10, 6a13)', "name", "E pred", "E meas", "phase pred", "phase meas", "lim R", "lim runs"

    do m = 1, size(names)
        tableau = rk_tableau(trim(names(m)))
        r = stability_r(tableau, cmplx(0.0_wp, h, wp))

        open(newunit=dat_e, file=out // "energy_" // trim(names(m)) // ".dat", status="replace", action="write")
        open(newunit=dat_t, file=out // "theory_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat_e, '(a)') "t drift"
        write(dat_t, '(a)') "t drift"
        open(newunit=dat_x, file=out // "tail_" // trim(names(m)) // ".dat", status="replace", action="write")
        open(newunit=dat_p, file=out // "portrait_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat_x, '(a)') "t x"
        write(dat_p, '(a)') "x v"
        write(dat_p, '(2es14.6)') 1.0_wp, 0.0_wp
        call rk_scheme(scheme, trim(names(m)), 2, h)
        y = [1.0_wp, 0.0_wp]
        time = 0.0_wp
        next_sample = 1
        do n = 1, n_steps
            call scheme%advance(oscillator, time, y)
            if (n == next_sample .or. n == n_steps) then          ! geometric sampling in time
                write(dat_e, '(2es14.6)') time, max(abs(energy(y) / 0.5_wp - 1.0_wp), 1.0e-16_wp)
                write(dat_t, '(2es14.6)') time, max(abs(abs(r)**(2 * n) - 1.0_wp), 1.0e-16_wp)
                next_sample = max(next_sample + 1, nint(next_sample * 1.25_wp))
            end if
            if (time >= n_steps * h - n_tail_periods * 2 * pi - h / 2) write(dat_x, '(2es14.6)') time, y(1)
            if (mod(n, portrait_every) == 0) write(dat_p, '(2es14.6)') y
        end do
        close(dat_e)
        close(dat_t)
        close(dat_x)
        close(dat_p)

        e_ratio = energy(y) / 0.5_wp
        e_predicted = abs(r)**(2 * n_steps)
        phase = wrap(atan2(-y(2), y(1)) - time)
        phase_predicted = wrap(n_steps * (atan2(aimag(r), real(r)) - h))
        lim_theory = limit_from_r(tableau)
        lim_measured = limit_from_runs(trim(names(m)))

        print '(a10, 6es13.4)', trim(names(m)), e_predicted, e_ratio, phase_predicted, phase, lim_theory, lim_measured
        write(table, '("\texttt{", a, "} & ", a, " & ", a, " & ", a, " & ", a, " & ", f6.4, " & ", f6.4, " \\")') &
            trim(names(m)), sci(e_predicted), sci(e_ratio), sci(phase_predicted), sci(phase), lim_theory, lim_measured
    end do

    write(table, '(a)') "\bottomrule"
    write(table, '(a)') "\end{tabular}"
    close(table)

contains

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
        write(u, '(a)') "% Written by examples/oscillator/main.f90: do not edit."
        write(u, '(a)') "\begin{tabular}{ll}"
        write(u, '(a)') "\toprule"
        write(u, '(a)') "frequency & $\omega = 1$ (period $2\pi$) \\"
        write(u, '(a)') "initial condition & $x(0) = 1$, $v(0) = 0$, $E_0 = 1/2$ \\"
        write(u, '(a)') "step & $h = \Delta t = " // fmt(h) // "$ \\"
        write(u, '(a)') "interval & " // itoa(n_periods) // " periods, " // itoa(n_steps) // " steps, $t_{\mathrm{end}} = " // &
                        fmt(n_steps * h) // "$ \\"
        write(u, '(a)') "methods & " // list // " \\"
        write(u, '(a)') "interface & plain arrays, \texttt{rk\_array\_scheme\_t} \\"
        write(u, '(a)') "figures & $x(t)$ over the last " // itoa(n_tail_periods) // &
                        " periods; $(x, v)$ every " // itoa(portrait_every) // " steps \\"
        write(u, '(a)') "phase measured & $\operatorname{atan2}(-v, x) - t$, wrapped to $(-\pi, \pi]$ \\"
        write(u, '(a)') "limit from $R$ & first $y$ with $|R(iy)| > 1$: scan by $10^{-3}$, then 60 bisections \\"
        write(u, '(a)') "limit measured & largest $h$ with $E_N \le E_0$ after $N = " // itoa(n_limit) // &
                        "$ steps: 50 bisections on $[0, 6]$ \\"
        write(u, '(a)') "\bottomrule"
        write(u, '(a)') "\end{tabular}"
        close(u)
    end subroutine write_parameters

    function fmt(x) result(s)
        real(wp), intent(in) :: x
        character(:), allocatable :: s
        character(len=16) :: buffer

        write(buffer, '(f10.3)') x
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

    !> x' = v, v' = -x
    subroutine oscillator(t, y, dydt, stage)
        real(wp), intent(in)  :: t
        real(wp), intent(in)  :: y(:)
        real(wp), intent(out) :: dydt(:)
        type(stage_t), intent(in) :: stage

        dydt(1) = y(2)
        dydt(2) = -y(1)
    end subroutine oscillator

    pure real(wp) function energy(y)
        real(wp), intent(in) :: y(2)

        energy = (y(1)**2 + y(2)**2) / 2
    end function energy

    !> Angle in (-pi, pi].
    pure real(wp) function wrap(a)
        real(wp), intent(in) :: a

        wrap = a - 2 * pi * nint(a / (2 * pi))
    end function wrap

    !> A number in LaTeX, e.g. $1.170081$, $-5.22\times10^{-4}$, $2.30\times10^{27}$.
    function sci(x) result(s)
        real(wp), intent(in) :: x
        character(:), allocatable :: s
        character(len=32) :: mantissa, expo
        character(len=1)  :: sign
        integer :: e

        sign = merge("-", " ", x < 0)
        if (abs(x) >= 0.01_wp .and. abs(x) < 1000.0_wp) then
            write(mantissa, '(f9.6)') abs(x)
            s = "$" // trim(sign) // trim(adjustl(mantissa)) // "$"
        else
            e = floor(log10(abs(x)))
            write(mantissa, '(f4.2)') abs(x) / 10.0_wp**e
            write(expo, '(i0)') e
            s = "$" // trim(sign) // trim(mantissa) // "\times10^{" // trim(expo) // "}$"
        end if
    end function sci

    !> Largest step with |R(iy)| <= 1 for all y in (0, h], from the tableau:
    !> scan from 0, then bisect; 0 if |R(iy)| > 1 for every small y.
    real(wp) function limit_from_r(tableau) result(lim)
        type(rk_tableau_t), intent(in) :: tableau

        real(wp) :: lo, hi
        integer  :: it

        lo = 0.0_wp
        do
            hi = lo + 1.0e-3_wp
            if (grows(tableau, hi) .or. hi > 6.0_wp) exit
            lo = hi
        end do
        if (lo == 0.0_wp) then
            lim = 0.0_wp
            return
        end if
        do it = 1, 60
            lim = (lo + hi) / 2
            if (grows(tableau, lim)) then
                hi = lim
            else
                lo = lim
            end if
        end do
        lim = lo
    end function limit_from_r

    logical function grows(tableau, y)
        type(rk_tableau_t), intent(in) :: tableau
        real(wp),           intent(in) :: y

        grows = abs(stability_r(tableau, cmplx(0.0_wp, y, wp))) > 1.0_wp
    end function grows

    !> Largest step for which n_limit steps do not increase the energy, by
    !> bisection on (0, 6). Below ~1e-4 the growth of a method unstable on the
    !> whole axis (|R(iy)|^2 = 1 + O(y^4)) is lost in round-off over n_limit steps.
    real(wp) function limit_from_runs(name) result(lim)
        character(*), intent(in) :: name

        real(wp) :: lo, hi
        integer  :: it

        lo = 0.0_wp
        hi = 6.0_wp
        do it = 1, 50
            lim = (lo + hi) / 2
            if (energy_after(name, lim) <= 0.5_wp) then
                lo = lim
            else
                hi = lim
            end if
        end do
        lim = lo
    end function limit_from_runs

    real(wp) function energy_after(name, dt) result(e)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: dt

        type(rk_array_scheme_t) :: s
        real(wp) :: z(2), t
        integer  :: i

        call rk_scheme(s, name, 2, dt)
        z = [1.0_wp, 0.0_wp]
        t = 0.0_wp
        do i = 1, n_limit
            call s%advance(oscillator, t, z)
        end do
        e = energy(z)
    end function energy_after

end program main
