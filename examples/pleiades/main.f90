program main
    !! Pleiades (SciML prob_ode_pleiades; Hairer, Norsett & Wanner, Solving
    !! ODEs I, section II.10): seven bodies in a plane, masses m_j = j (G = 1),
    !!   x_i'' = sum_{j /= i} m_j (x_j - x_i) / r_ij^3,  y_i'' likewise,
    !!   r_ij = ((x_i - x_j)^2 + (y_i - y_j)^2)^(1/2),   0 <= t <= 3.
    !! No closed form: the reference is dp54 with n_reference steps, its accuracy
    !! estimated against the run with half as many. Tests:
    !!   - convergence: the error at t = 1 and at t = 3 as the steps double
    !!     (before and after most close encounters);
    !!   - invariants: the total momentum P (linear: kept to round-off by every
    !!     Runge-Kutta method), the angular momentum L and the energy E (not
    !!     linear: not kept), over a run of n_invariants steps.
    !! Uses the plain-array schemes, state (x(7), y(7), x'(7), y'(7)). Writes
    !! examples/pleiades/data/: parameters.tex, results.tex, orbits.dat (t and
    !! the reference positions of the 7 bodies), encounter.dat (t, smallest
    !! distance between two bodies), history_<method>.dat (t, position error at
    !! n_invariants steps), energy_<method>.dat (t, |E - E0|), error1_<method>.dat
    !! and error_<method>.dat (dt, error at t = 1 and t = 3), read by
    !! summary.tex. Run from the project root.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: rk_scheme, stage_t, rk_array_scheme_t, rk_method_names, rk_tableau, rk_tableau_t
    implicit none

    character(*), parameter :: out = "examples/pleiades/data/"
    integer,      parameter :: nb = 7                          !! bodies
    real(wp),     parameter :: t_end = 3.0_wp
    integer,      parameter :: n_reference = 3 * 2**19         !! reference steps (dp54)
    integer,      parameter :: n_points = 3 * 2**9             !! points written per curve (t = 1 is point 512)
    integer,      parameter :: n_invariants = 3 * 2**10        !! steps of the invariant and history runs
    integer,      parameter :: n_start = 3 * 2**7, n_levels = 10  !! convergence: 384 ... 196608 steps
    real(wp),     parameter :: floor_1 = 3.0e-11_wp            !! error floor at t = 1 (reference ~1e-12)
    real(wp),     parameter :: floor_3 = 1.0e-8_wp             !! error floor at t = 3 (reference ~1e-10)

    character(len=8), allocatable :: names(:)
    type(rk_tableau_t) :: tableau
    real(wp) :: m(nb), u0(4 * nb), reference_end(4 * nb), check_end(4 * nb), ref_accuracy, u_mid(4 * nb)
    real(wp), allocatable :: reference(:,:)
    real(wp) :: dt(n_levels), err(n_levels), err1(n_levels), observed, observed1, e0, p0(2), l0, de, dp, dl
    integer  :: k, mth, l, table, dat, dat1

    m = [(real(k, wp), k = 1, nb)]
    u0 = 0.0_wp
    u0(1:nb)        = [3.0_wp, 3.0_wp, -1.0_wp, -3.0_wp, 2.0_wp, -2.0_wp, 2.0_wp]      ! x
    u0(nb+1:2*nb)   = [3.0_wp, -3.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, -4.0_wp, 4.0_wp]      ! y
    u0(2*nb+6)      = 1.75_wp                                                         ! x6'
    u0(2*nb+7)      = -1.5_wp                                                         ! x7'
    u0(3*nb+4)      = -1.25_wp                                                        ! y4'
    u0(3*nb+5)      = 1.0_wp                                                          ! y5'
    e0 = energy(u0)
    p0 = momentum(u0)
    l0 = angular(u0)
    names = rk_method_names()
    call execute_command_line("mkdir -p " // out)

    ! Reference, its accuracy, the orbits and the closest encounters.
    allocate(reference(4 * nb, 0:n_points))
    call run("dp54", n_reference, reference_end, de, dp, dl, keep=reference)
    call run("dp54", n_reference / 2, check_end, de, dp, dl)
    ref_accuracy = maxval(abs(reference_end - check_end))
    print '(a, es10.2)', "reference: difference to the run with half the steps ", ref_accuracy
    call write_reference()
    call write_parameters()

    do mth = 1, size(names)
        call run(trim(names(mth)), n_invariants, check_end, de, dp, dl, &
                 history="history_" // trim(names(mth)) // ".dat", energy_file="energy_" // trim(names(mth)) // ".dat")
    end do

    open(newunit=table, file=out // "results.tex", status="replace", action="write")
    write(table, '(a)') "% Written by examples/pleiades/main.f90: do not edit."
    write(table, '(a)') "\begin{tabular}{lccccccc}"
    write(table, '(a)') "\toprule"
    write(table, '(a)') "& & \multicolumn{2}{c}{observed order} & \multicolumn{3}{c}{largest drift, " // &
                        itoa(n_invariants) // " steps} \\"
    write(table, '(a)') "\cmidrule(lr){3-4}\cmidrule(lr){5-7}"
    write(table, '(a)') "name & order & $t = 1$ & $t = 3$ & momentum $|P - P_0|$ & angular $|L - L_0|$ & " // &
                        "energy $|E - E_0|$ \\"
    write(table, '(a)') "\midrule"
    print '(a10, a6, 2a9, 3a12)', "name", "order", "obs. t=1", "t=3", "|dP|", "|dL|", "|dE|"
    do mth = 1, size(names)
        tableau = rk_tableau(trim(names(mth)))
        open(newunit=dat, file=out // "error_" // trim(names(mth)) // ".dat", status="replace", action="write")
        open(newunit=dat1, file=out // "error1_" // trim(names(mth)) // ".dat", status="replace", action="write")
        write(dat, '(a)') "dt error"
        write(dat1, '(a)') "dt error"
        do l = 1, n_levels
            dt(l) = t_end / (n_start * 2**(l - 1))
            call run(trim(names(mth)), n_start * 2**(l - 1), check_end, de, dp, dl, u_mid=u_mid)
            err(l) = maxval(abs(check_end - reference_end))
            err1(l) = maxval(abs(u_mid - reference(:, n_points / 3)))
            if (err(l) /= err(l)) err(l) = huge(1.0_wp)
            if (err1(l) /= err1(l)) err1(l) = huge(1.0_wp)
            if (err(l) > 1.0e-15_wp .and. err(l) < 1.0e3_wp) write(dat, '(2es14.6)') dt(l), err(l)
            if (err1(l) > 1.0e-15_wp .and. err1(l) < 1.0e3_wp) write(dat1, '(2es14.6)') dt(l), err1(l)
        end do
        close(dat)
        close(dat1)
        observed = slope(dt, err, floor_3)
        observed1 = slope(dt, err1, floor_1)
        call run(trim(names(mth)), n_invariants, check_end, de, dp, dl)
        print '(a10, i6, 2f9.3, 3es12.3)', trim(names(mth)), tableau%order, observed1, observed, dp, dl, de
        write(table, '("\texttt{", a, "} & ", i0, " & ", a, " & ", a, " & ", a, " & ", a, " & ", a, " \\")') &
            trim(names(mth)), tableau%order, order_text(observed1), order_text(observed), sci(dp), sci(dl), sci(de)
    end do
    write(table, '(a)') "\bottomrule"
    write(table, '(a)') "\end{tabular}"
    close(table)

contains

    subroutine pleiades(t, u, dudt, stage)
        real(wp), intent(in)  :: t
        real(wp), intent(in)  :: u(:)
        real(wp), intent(out) :: dudt(:)
        type(stage_t), intent(in) :: stage
        real(wp) :: dx, dy, r3
        integer  :: i, j

        dudt(1:2*nb) = u(2*nb+1:4*nb)
        dudt(2*nb+1:4*nb) = 0.0_wp
        do i = 1, nb
            do j = 1, nb
                if (j == i) cycle
                dx = u(j) - u(i)
                dy = u(nb+j) - u(nb+i)
                r3 = (dx * dx + dy * dy)**1.5_wp
                dudt(2*nb+i) = dudt(2*nb+i) + m(j) * dx / r3
                dudt(3*nb+i) = dudt(3*nb+i) + m(j) * dy / r3
            end do
        end do
    end subroutine pleiades

    pure real(wp) function energy(u) result(e)
        real(wp), intent(in) :: u(:)
        integer :: i, j

        e = sum(m * (u(2*nb+1:3*nb)**2 + u(3*nb+1:4*nb)**2)) / 2
        do i = 1, nb
            do j = i + 1, nb
                e = e - m(i) * m(j) / hypot(u(i) - u(j), u(nb+i) - u(nb+j))
            end do
        end do
    end function energy

    pure function momentum(u) result(p)
        real(wp), intent(in) :: u(:)
        real(wp) :: p(2)

        p = [sum(m * u(2*nb+1:3*nb)), sum(m * u(3*nb+1:4*nb))]
    end function momentum

    pure real(wp) function angular(u)
        real(wp), intent(in) :: u(:)

        angular = sum(m * (u(1:nb) * u(3*nb+1:4*nb) - u(nb+1:2*nb) * u(2*nb+1:3*nb)))
    end function angular

    !> Smallest distance between two bodies, and which two.
    real(wp) function closest(u, pair)
        real(wp), intent(in)            :: u(:)
        integer,  intent(out), optional :: pair(2)
        real(wp) :: r
        integer  :: i, j

        closest = huge(1.0_wp)
        do i = 1, nb
            do j = i + 1, nb
                r = hypot(u(i) - u(j), u(nb+i) - u(nb+j))
                if (r < closest) then
                    closest = r
                    if (present(pair)) pair = [i, j]
                end if
            end do
        end do
    end function closest

    !> n steps to t_end: the final state, the largest drifts of E, P, L.
    !> Optional: keep n_points + 1 states (requires n divisible by n_points),
    !> write the position error against the kept reference, write |E - E0|.
    subroutine run(name, n, u_end, de, dp, dl, keep, history, energy_file, u_mid)
        character(*), intent(in)  :: name
        integer,      intent(in)  :: n
        real(wp),     intent(out) :: u_end(:), de, dp, dl
        real(wp),     intent(out), optional :: u_mid(:)        !! the state at t = t_end/3 (n divisible by 3)
        real(wp),     intent(inout), optional :: keep(:, 0:)
        character(*), intent(in),    optional :: history, energy_file

        type(rk_array_scheme_t) :: scheme
        real(wp) :: u(4 * nb), time
        integer  :: i, stride, uh, ue

        call rk_scheme(scheme, name, 4 * nb, t_end / n)
        u = u0
        time = 0.0_wp
        de = 0.0_wp
        dp = 0.0_wp
        dl = 0.0_wp
        stride = max(1, n / n_points)
        if (present(keep)) keep(:, 0) = u
        if (present(history)) then
            open(newunit=uh, file=out // history, status="replace", action="write")
            write(uh, '(a)') "t error"
        end if
        if (present(energy_file)) then
            open(newunit=ue, file=out // energy_file, status="replace", action="write")
            write(ue, '(a)') "t drift"
        end if
        do i = 1, n
            call scheme%advance(pleiades, time, u)
            if (any(u /= u)) exit
            if (present(u_mid) .and. i == n / 3) u_mid = u
            de = max(de, abs(energy(u) - e0))
            dp = max(dp, maxval(abs(momentum(u) - p0)))
            dl = max(dl, abs(angular(u) - l0))
            if (mod(i, stride) == 0) then
                if (present(keep)) keep(:, i / stride) = u
                if (present(history)) write(uh, '(2es14.6)') time, &
                    max(maxval(abs(u(1:2*nb) - reference(1:2*nb, i * n_points / n))), 1.0e-16_wp)
                if (present(energy_file)) write(ue, '(2es14.6)') time, max(abs(energy(u) - e0), 1.0e-16_wp)
            end if
        end do
        u_end = u
        if (present(history)) close(uh)
        if (present(energy_file)) close(ue)
    end subroutine run

    !> orbits.dat: t, x1..x7, y1..y7 of the reference; encounter.dat: t, the
    !> smallest distance between two bodies; encounters.dat: the close
    !> encounters (local minima of that distance below 1): t, distance, pair.
    subroutine write_reference()
        integer :: i, u, ue, un, pair(2)
        real(wp) :: d(0:n_points)

        open(newunit=u, file=out // "orbits.dat", status="replace", action="write")
        open(newunit=ue, file=out // "encounter.dat", status="replace", action="write")
        write(u, '(a)') "t x1 x2 x3 x4 x5 x6 x7 y1 y2 y3 y4 y5 y6 y7"
        write(ue, '(a)') "t distance"
        do i = 0, n_points
            d(i) = closest(reference(:, i), pair)
            write(u, '(15es15.7)') t_end * i / n_points, reference(1:2*nb, i)
            write(ue, '(2es14.6)') t_end * i / n_points, d(i)
        end do
        close(u)
        close(ue)
        open(newunit=un, file=out // "encounters.dat", status="replace", action="write")
        write(un, '(a)') "t distance pair"
        do i = 1, n_points - 1
            if (d(i) < d(i - 1) .and. d(i) <= d(i + 1) .and. d(i) < 1.0_wp) then
                d(i) = closest(reference(:, i), pair)
                write(un, '(2es14.6, 1x, a)') t_end * i / n_points, d(i), itoa(pair(1)) // "--" // itoa(pair(2))
            end if
        end do
        close(un)
    end subroutine write_reference

    !> The observed order: the slope of log(err) against log(dt) between the
    !> two finest levels with errors in (floor, 0.1), reported only if it has
    !> settled, agreeing within 0.3 with the slope of the next coarser pair;
    !> 0 (shown as a dash) otherwise.
    real(wp) function slope(dt, err, floor) result(p)
        real(wp), intent(in) :: dt(:), err(:), floor
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
        character(:), allocatable :: list

        list = ""
        do j = 1, size(names)
            if (j > 1) list = list // ", "
            list = list // "\texttt{" // trim(names(j)) // "}"
        end do
        open(newunit=u, file=out // "parameters.tex", status="replace", action="write")
        write(u, '(a)') "% Written by examples/pleiades/main.f90: do not edit."
        write(u, '(a)') "\begin{tabular}{l>{\raggedright\arraybackslash}p{12.4cm}}"
        write(u, '(a)') "\toprule"
        write(u, '(a)') "masses & $m_j = j$, $j = 1, \dots, 7$; $G = 1$ \\"
        write(u, '(a)') "positions & $(3,3), (3,-3), (-1,2), (-3,0), (2,0), (-2,-4), (2,4)$ \\"
        write(u, '(a)') "velocities & zero except $\dot x_6 = 1.75$, $\dot x_7 = -1.5$, $\dot y_4 = -1.25$, $\dot y_5 = 1$ \\"
        write(u, '(a)') "interval & $0 \le t \le 3$ \\"
        write(u, '(a)') "reference & \texttt{dp54}, " // itoa(n_reference) // " steps; " // sci(ref_accuracy) // &
                        " from the run with half the steps \\"
        write(u, '(a)') "convergence & error $\max_i |u_i - u_i^{\mathrm{ref}}|$ at $t = 1$, $3$; " // itoa(n_start) // ", " // &
                        itoa(2 * n_start) // ", \dots, " // itoa(n_start * 2**(n_levels - 1)) // " steps \\"
        write(u, '(a)') "observed order & two finest levels above " // sci(floor_1) // " ($t = 1$), " // &
                        sci(floor_3) // " ($t = 3$); kept if within 0.3 of the next slope \\"
        write(u, '(a)') "invariants, history & " // itoa(n_invariants) // " steps ($\Delta t = " // &
                        fmt(t_end / n_invariants) // "$); $E_0 = " // fmt(e0) // "$, $L_0 = " // fmt(l0) // &
                        "$, $P_0 = (" // fmt(p0(1)) // ", " // fmt(p0(2)) // ")$ \\"
        write(u, '(a)') "methods & " // list // " \\"
        write(u, '(a)') "interface & plain arrays, \texttt{rk\_array\_scheme\_t}, 28 unknowns \\"
        write(u, '(a)') "\bottomrule"
        write(u, '(a)') "\end{tabular}"
        close(u)
    end subroutine write_parameters

    function sci(x) result(s)
        real(wp), intent(in) :: x
        character(:), allocatable :: s
        character(len=32) :: mantissa
        integer :: ex

        if (x /= x .or. x >= huge(1.0_wp)) then
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
        if (s == "-0") s = "0"
    end function fmt

    function itoa(i) result(s)
        integer, intent(in) :: i
        character(:), allocatable :: s
        character(len=16) :: buffer

        write(buffer, '(i0)') i
        s = trim(buffer)
    end function itoa

end program main
