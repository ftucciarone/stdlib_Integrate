program main
    !! Robertson's chemical kinetics (SciML prob_ode_rober), the classic stiff
    !! problem:
    !!   y1' = -0.04 y1 + 1e4 y2 y3,
    !!   y2' =  0.04 y1 - 1e4 y2 y3 - 3e7 y2^2,
    !!   y3' =  3e7 y2^2,                    y(0) = (1, 0, 0).
    !! y1 + y2 + y3 = 1 for all t. The Jacobian has the eigenvalue 0 (that
    !! invariant) and two more, from its trace T and the sum D of its principal
    !! 2x2 minors: lambda = (T +- sqrt(T^2 - 4D))/2, with
    !!   T = -0.04 - 1e4 y3 - 6e7 y2,  D = 2.4e6 y2 + 6e11 y2^2.
    !! The fast one, lambda_f, is of order -1e3 to -1e4 at all times, while the
    !! solution evolves on times of order 1 to 1e5: an explicit step is stable
    !! only while |lambda_f(t)| dt <= h*. Tests:
    !!   - the largest step that reaches T without blowing up, for T = 1, 10,
    !!     100, 1000, against h*/max_{t <= T} |lambda_f(t)|;
    !!   - the steps an explicit method needs to reach t_end;
    !!   - (figure) for a sweep of steps, the time the run blows up, against the
    !!     predicted onset, the first t with |lambda_f(t)| dt > h*.
    !! Uses the plain-array schemes. Writes examples/robertson/data/:
    !! parameters.tex, results.tex, reference.dat (t, y1, y2, y3, |lambda_f|),
    !! fixed_<k>.dat (t, y2 of rk4 at three steps), blowup_<method>.dat (dt,
    !! measured blow-up time), onset_<method>.dat (dt, predicted onset),
    !! reach_<method>.dat (T, measured and predicted largest step), read by
    !! summary.tex. Run from the project root.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: rk_scheme, stage_t, rk_array_scheme_t, rk_method_names, rk_tableau
    use robertson_stability, only: stability_r
    implicit none

    character(*), parameter :: out = "examples/robertson/data/"
    real(wp),     parameter :: t_end = 1000.0_wp
    real(wp),     parameter :: dt_reference = 5.0e-6_wp
    integer,      parameter :: n_samples = 400                !! reference samples, log-spaced in t
    integer,      parameter :: n_sweep = 14                   !! steps dt, log-spaced
    real(wp),     parameter :: dt_min = 1.0e-4_wp, dt_max = 1.0e-2_wp

    character(len=8), allocatable :: names(:)
    real(wp) :: ts(0:n_samples), ys(3, 0:n_samples), lf(0:n_samples)
    real(wp), parameter :: horizons(4) = [1.0_wp, 10.0_wp, 100.0_wp, 1000.0_wp]
    real(wp) :: h_star, dt, t_blow, t_onset, steps_needed, measured, predicted, ratio(size(horizons))
    integer  :: m, k, table, ub, uo, ur

    names = rk_method_names()
    call execute_command_line("mkdir -p " // out)
    call reference_run()
    call write_parameters()
    call write_fixed()

    open(newunit=table, file=out // "results.tex", status="replace", action="write")
    write(table, '(a)') "% Written by examples/robertson/main.f90: do not edit."
    write(table, '(a)') "\begin{tabular}{lcccccc}"
    write(table, '(a)') "\toprule"
    write(table, '(a)') "& & \multicolumn{4}{c}{largest step reaching $T$, measured / predicted} & \\"
    write(table, '(a)') "\cmidrule(lr){3-6}"
    write(table, '(a)') "name & $h^*$ & $T = 1$ & $T = 10$ & $T = 100$ & $T = 1000$ & steps to $t = " // &
                        fmt(t_end) // "$ \\"
    write(table, '(a)') "\midrule"
    print '(a10, a8, 4a9, a14)', "name", "h*", "T=1", "T=10", "T=100", "T=1000", "steps"
    do m = 1, size(names)
        h_star = real_axis_limit(trim(names(m)))
        open(newunit=ub, file=out // "blowup_" // trim(names(m)) // ".dat", status="replace", action="write")
        open(newunit=uo, file=out // "onset_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(ub, '(a)') "dt t"
        write(uo, '(a)') "dt t"
        do k = 0, n_sweep - 1
            dt = dt_min * (dt_max / dt_min)**(real(k, wp) / (n_sweep - 1))
            t_onset = onset(h_star, dt)
            t_blow = blow_up_time(trim(names(m)), dt)
            if (t_onset < t_end) write(uo, '(2es14.6)') dt, t_onset
            if (t_blow < t_end) write(ub, '(2es14.6)') dt, t_blow
        end do
        close(ub)
        close(uo)
        open(newunit=ur, file=out // "reach_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(ur, '(a)') "T measured predicted"
        do k = 1, size(horizons)
            predicted = h_star / maxval(lf, mask=ts <= horizons(k))
            measured = reaching_step(trim(names(m)), horizons(k), predicted)
            ratio(k) = measured / predicted
            write(ur, '(3es14.6)') horizons(k), measured, predicted
        end do
        close(ur)
        steps_needed = t_end * maxval(lf) / h_star
        print '(a10, f8.4, 4f9.4, es14.3)', trim(names(m)), h_star, ratio, steps_needed
        write(table, '("\texttt{", a, "} & ", f6.4, 4(" & ", f6.4), " & ", a, " \\")') &
            trim(names(m)), h_star, ratio, sci(steps_needed)
    end do
    write(table, '(a)') "\bottomrule"
    write(table, '(a)') "\end{tabular}"
    close(table)

contains

    subroutine robertson(t, y, dydt, stage)
        real(wp), intent(in)  :: t
        real(wp), intent(in)  :: y(:)
        real(wp), intent(out) :: dydt(:)
        type(stage_t), intent(in) :: stage

        dydt(1) = -0.04_wp * y(1) + 1.0e4_wp * y(2) * y(3)
        dydt(2) =  0.04_wp * y(1) - 1.0e4_wp * y(2) * y(3) - 3.0e7_wp * y(2)**2
        dydt(3) =  3.0e7_wp * y(2)**2
    end subroutine robertson

    !> |lambda_f|, the fast eigenvalue of the Jacobian, in size.
    pure real(wp) function fast(y)
        real(wp), intent(in) :: y(3)
        real(wp) :: tr, d

        tr = -0.04_wp - 1.0e4_wp * y(3) - 6.0e7_wp * y(2)
        d = 2.4e6_wp * y(2) + 6.0e11_wp * y(2)**2
        fast = abs((tr - sqrt(max(tr * tr - 4 * d, 0.0_wp))) / 2)
    end function fast

    !> dp54 at dt_reference to t_end, sampled at log-spaced times from 1e-6.
    subroutine reference_run()
        type(rk_array_scheme_t) :: scheme
        real(wp) :: y(3), time, next
        integer  :: i, k, unit

        call rk_scheme(scheme, "dp54", 3, dt_reference)
        y = [1.0_wp, 0.0_wp, 0.0_wp]
        time = 0.0_wp
        ts(0) = 0.0_wp
        ys(:, 0) = y
        lf(0) = fast(y)
        k = 1
        next = 1.0e-6_wp
        i = 0
        do while (k <= n_samples)
            call scheme%advance(robertson, time, y)
            i = i + 1
            if (i * dt_reference >= next) then
                ts(k) = i * dt_reference
                ys(:, k) = y
                lf(k) = fast(y)
                k = k + 1
                next = 1.0e-6_wp * (t_end / 1.0e-6_wp)**(real(k - 1, wp) / (n_samples - 1))
            end if
        end do
        open(newunit=unit, file=out // "reference.dat", status="replace", action="write")
        write(unit, '(a)') "t y1 y2 y3 lambda"
        do k = 1, n_samples
            write(unit, '(5es14.6)') ts(k), ys(:, k), lf(k)
        end do
        close(unit)
    end subroutine reference_run

    !> Predicted onset: the first sampled time with |lambda_f| dt > h*.
    real(wp) function onset(h_star, dt) result(t)
        real(wp), intent(in) :: h_star, dt
        integer :: k

        t = huge(1.0_wp)
        do k = 0, n_samples
            if (lf(k) * dt > h_star) then
                t = ts(k)
                return
            end if
        end do
    end function onset

    !> First time a run with step dt fails: a value not finite or outside
    !> [-0.01, 1.01], y2 below -1e-6, or (after t = 1e-3, once y2 has settled)
    !> y2 off the reference by more than half. The last catches an instability
    !> that grows without changing sign (R(-h) > 1, as for dp54): the stiff
    !> term 3e7 y2^2 then holds y2 at a wrong but bounded value.
    real(wp) function blow_up_time(name, dt, horizon) result(t)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: dt
        real(wp),     intent(in), optional :: horizon     !! run to this time (default t_end)
        type(rk_array_scheme_t) :: scheme
        real(wp) :: y(3), time, t_stop, y2_ref
        integer  :: k
        integer  :: i

        call rk_scheme(scheme, name, 3, dt)
        y = [1.0_wp, 0.0_wp, 0.0_wp]
        time = 0.0_wp
        t = huge(1.0_wp)
        t_stop = t_end
        if (present(horizon)) t_stop = horizon                ! not merge(): it would read an absent argument
        k = 1
        do i = 1, ceiling(t_stop / dt)
            call scheme%advance(robertson, time, y)
            if (.not. (all(y >= -0.01_wp .and. y <= 1.01_wp) .and. y(2) >= -1.0e-6_wp)) then
                t = i * dt
                return
            end if
            if (i * dt > 1.0e-3_wp) then
                do while (k < n_samples - 1 .and. ts(k + 1) < i * dt)
                    k = k + 1
                end do
                y2_ref = ys(2, k) + (ys(2, k + 1) - ys(2, k)) * (i * dt - ts(k)) / (ts(k + 1) - ts(k))
                if (abs(y(2) - y2_ref) > 0.5_wp * y2_ref) then
                    t = i * dt
                    return
                end if
            end if
        end do
    end function blow_up_time

    !> Largest step (bisection between 0.7 and 1.5 times `guess`) with which
    !> the run reaches `horizon` without blowing up.
    real(wp) function reaching_step(name, horizon, guess) result(h)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: horizon, guess
        real(wp) :: lo, hi
        integer  :: it

        lo = 0.7_wp * guess
        hi = 1.5_wp * guess
        do it = 1, 14
            h = (lo + hi) / 2
            if (blow_up_time(name, h, horizon) > horizon) then
                lo = h
            else
                hi = h
            end if
        end do
        h = lo
    end function reaching_step

    !> rk4 at three fixed steps: (t, y2) up to the blow-up, log-spaced output.
    subroutine write_fixed()
        real(wp), parameter :: steps(3) = [2.0e-3_wp, 7.0e-4_wp, 3.0e-4_wp]
        type(rk_array_scheme_t) :: scheme
        real(wp) :: y(3), time, next
        integer  :: i, k, unit

        do k = 1, 3
            call rk_scheme(scheme, "rk4", 3, steps(k))
            y = [1.0_wp, 0.0_wp, 0.0_wp]
            time = 0.0_wp
            next = steps(k)
            open(newunit=unit, file=out // "fixed_" // itoa(k) // ".dat", status="replace", action="write")
            write(unit, '(a)') "t y2"
            do i = 1, ceiling(t_end / steps(k))
                call scheme%advance(robertson, time, y)
                if (.not. (abs(y(2)) < 1.0e-3_wp)) exit
                if (i * steps(k) >= next) then
                    write(unit, '(2es14.6)') i * steps(k), max(y(2), 1.0e-12_wp)
                    next = next * 1.05_wp
                end if
            end do
            close(unit)
        end do
    end subroutine write_fixed

    real(wp) function real_axis_limit(name) result(h)
        character(*), intent(in) :: name
        real(wp) :: lo, hi
        integer  :: it

        lo = 0.0_wp
        do
            hi = lo + 1.0e-3_wp
            if (abs(stability_r(rk_tableau(name), cmplx(-hi, 0.0_wp, wp))) > 1.0_wp) exit
            lo = hi
        end do
        do it = 1, 60
            h = (lo + hi) / 2
            if (abs(stability_r(rk_tableau(name), cmplx(-h, 0.0_wp, wp))) > 1.0_wp) then
                hi = h
            else
                lo = h
            end if
        end do
        h = lo
    end function real_axis_limit

    subroutine write_parameters()
        integer :: u, j
        character(:), allocatable :: list

        list = ""
        do j = 1, size(names)
            if (j > 1) list = list // ", "
            list = list // "\texttt{" // trim(names(j)) // "}"
        end do
        open(newunit=u, file=out // "parameters.tex", status="replace", action="write")
        write(u, '(a)') "% Written by examples/robertson/main.f90: do not edit."
        write(u, '(a)') "\begin{tabular}{l>{\raggedright\arraybackslash}p{12.4cm}}"
        write(u, '(a)') "\toprule"
        write(u, '(a)') "rate constants & $k_1 = 0.04$, $k_2 = 3\times10^7$, $k_3 = 10^4$ \\"
        write(u, '(a)') "initial condition & $y = (1, 0, 0)$; interval $0 \le t \le " // fmt(t_end) // "$ \\"
        write(u, '(a)') "reference & \texttt{dp54}, $\Delta t = 5\times10^{-6}$; " // itoa(n_samples) // &
                        " samples, log-spaced from $t = 10^{-6}$ \\"
        write(u, '(a)') "fast eigenvalue & $|\lambda_f| = |T - \sqrt{T^2 - 4D}|/2$, " // &
                        "$T = -0.04 - 10^4 y_3 - 6\times10^7 y_2$, $D = 2.4\times10^6 y_2 + 6\times10^{11} y_2^2$ \\"
        write(u, '(a)') "step sweep & " // itoa(n_sweep) // " steps, log-spaced from $10^{-4}$ to $10^{-2}$ \\"
        write(u, '(a)') "failure & first $t$ with a value not finite or outside $[-0.01, 1.01]$, $y_2 < -10^{-6}$, or" // &
                        " (for $t > 10^{-3}$) $y_2$ off the reference by more than half \\"
        write(u, '(a)') "predicted onset & first sampled $t$ with $|\lambda_f(t)|\,\Delta t > h^*$ \\"
        write(u, '(a)') "reaching $T$ & largest $\Delta t$ whose run reaches $T$ without blowing up; bisection" // &
                        " (14 halvings) in $[0.7, 1.5]$ times $h^*/\max_{t \le T}|\lambda_f|$ \\"
        write(u, '(a)') "figures & \texttt{rk4} at $\Delta t = 2\times10^{-3}$, $7\times10^{-4}$, $3\times10^{-4}$ \\"
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
        integer :: ex

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
