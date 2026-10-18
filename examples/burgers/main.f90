program main
    !! Strong stability preservation on Burgers' equation,
    !!   u_t + (u^2/2)_x = 0 on [0, 1), periodic,
    !!   u(x, 0) = 1 on 0.25 < x < 0.6, -0.5 elsewhere,
    !! a square wave: a shock and a rarefaction from each jump at once. Space:
    !! MUSCL reconstruction with the minmod limiter and the Rusanov flux on N
    !! cells: with forward Euler this scheme does not increase the total
    !! variation TV(u) = sum |u_{i+1} - u_i| (TVD) for CFL = max|u| dt/dx <= 1/2.
    !! A Runge-Kutta method with SSP coefficient C keeps that property for
    !! CFL <= C/2. Every method of the library (rk3ws, rk4, dp54) has C = 0:
    !! no guarantee at any step (an SSP method, such as SSPRK3 with C = 1, runs
    !! on the generic engine from its tableau). Tests, for CFL from
    !! 0.05 to 1.2: the largest increase of TV over one step, and the new
    !! extrema (max u - 1, -0.5 - min u), over 0 <= t <= t_end.
    !! Uses the plain-array schemes. Writes examples/burgers/data/:
    !! parameters.tex, results.tex, spacetime.dat (x, t, u of the reference),
    !! reference.dat and final_<method>.dat (x, u at t_end, CFL = cfl_show),
    !! tv_<method>.dat (CFL, largest TV increase), extrema_<method>.dat (CFL,
    !! largest new extremum), and at CFL = cfl_show early_<method>.dat (x, u
    !! after two steps), worst_<method>.dat (x, u at the largest new extremum),
    !! overshoot_<method>.dat (t, new extremum), read by summary.tex. Run from
    !! the project root.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: rk_scheme, stage_t, rk_array_scheme_t, rk_method_names
    implicit none

    character(*), parameter :: out = "examples/burgers/data/"
    real(wp),     parameter :: pi = acos(-1.0_wp)
    integer,      parameter :: nx = 200                       !! cells
    real(wp),     parameter :: t_end = 0.4_wp                 !! the shock and rarefaction interact
    real(wp),     parameter :: cfl_show = 0.6_wp              !! CFL of the solution figure
    real(wp),     parameter :: cfl_reference = 0.05_wp        !! CFL of the reference (rk4)
    integer,      parameter :: n_cfl = 24                     !! CFL = 0.05, 0.10, ..., 1.20
    real(wp),     parameter :: tolerance = 1.0e-12_wp         !! a TV increase above this counts

    character(len=8), allocatable :: names(:)
    real(wp) :: dx, x(nx), u(nx), cfl, tv_up, new_ext
    real(wp), allocatable :: tvd_limit(:), tv_at(:,:), ext_at(:,:)   ! per method of the table
    integer  :: m, k, unit, table

    names = rk_method_names()
    allocate(tvd_limit(size(names)), tv_at(size(names), n_cfl), ext_at(size(names), n_cfl))
    dx = 1.0_wp / nx
    x = [((k - 0.5_wp) * dx, k = 1, nx)]
    call execute_command_line("mkdir -p " // out)
    call write_parameters()

    ! Reference (rk4 at a small CFL) in space-time, and the solutions at cfl_show.
    call run("rk4", cfl_reference, tv_up, new_ext, final="reference.dat", spacetime="spacetime.dat")
    do m = 1, size(names)
        call run(trim(names(m)), cfl_show, tv_up, new_ext, final="final_" // trim(names(m)) // ".dat", &
                 worst="worst_" // trim(names(m)) // ".dat", history="overshoot_" // trim(names(m)) // ".dat", &
                 early="early_" // trim(names(m)) // ".dat")
    end do

    ! The CFL sweep.
    do m = 1, size(names)
        tvd_limit(m) = 0.0_wp
        open(newunit=unit, file=out // "tv_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(unit, '(a)') "cfl increase"
        do k = 1, n_cfl
            cfl = 0.05_wp * k
            call run(trim(names(m)), cfl, tv_at(m, k), ext_at(m, k))
            if (tv_at(m, k) < 1.0e3_wp) write(unit, '(2es14.6)') cfl, max(tv_at(m, k), 1.0e-16_wp)  ! not if it blew up
            if (tv_at(m, k) <= tolerance .and. tvd_limit(m) == 0.05_wp * (k - 1)) tvd_limit(m) = cfl
        end do
        close(unit)
        open(newunit=unit, file=out // "extrema_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(unit, '(a)') "cfl extremum"
        do k = 1, n_cfl
            if (ext_at(m, k) < 1.0e3_wp) write(unit, '(2es14.6)') 0.05_wp * k, max(ext_at(m, k), 1.0e-16_wp)
        end do
        close(unit)
    end do

    open(newunit=table, file=out // "results.tex", status="replace", action="write")
    write(table, '(a)') "% Written by examples/burgers/main.f90: do not edit."
    write(table, '(a)') "\begin{tabular}{lcccccc}"
    write(table, '(a)') "\toprule"
    write(table, '(a)') "& & & \multicolumn{2}{c}{largest TV increase per step} & " // &
                        "\multicolumn{2}{c}{largest new extremum} \\"
    write(table, '(a)') "\cmidrule(lr){4-5}\cmidrule(lr){6-7}"
    write(table, '(a)') "name & SSP coeff. & TVD up to CFL & at CFL 0.45 & at CFL 0.60 & " // &
                        "at CFL 0.45 & at CFL 0.60 \\"
    write(table, '(a)') "\midrule"
    print '(a10, a6, a10, 4a12)', "name", "C", "TVD to", "TV 0.45", "TV 0.60", "ext 0.45", "ext 0.60"
    do m = 1, size(names)
        print '(a10, i6, f10.2, 4es12.3)', trim(names(m)), ssp(trim(names(m))), tvd_limit(m), &
            tv_at(m, 9), tv_at(m, 12), ext_at(m, 9), ext_at(m, 12)
        write(table, '("\texttt{", a, "} & ", i0, " & ", a, " & ", a, " & ", a, " & ", a, " & ", a, " \\")') &
            trim(names(m)), ssp(trim(names(m))), limit_text(tvd_limit(m)), sci(tv_at(m, 9)), sci(tv_at(m, 12)), &
            sci(ext_at(m, 9)), sci(ext_at(m, 12))
    end do
    write(table, '(a)') "\bottomrule"
    write(table, '(a)') "\end{tabular}"
    close(table)

contains

    !> SSP coefficient of the method (Gottlieb, Ketcheson & Shu, 2011).
    integer function ssp(name)
        character(*), intent(in) :: name

        select case (name)
        case default;                     ssp = 0      ! no SSP method left in the table
        end select
    end function ssp

    !> Semi-discrete Burgers: MUSCL-minmod reconstruction, Rusanov flux.
    subroutine burgers(t, u, dudt, stage)
        real(wp), intent(in)  :: t
        real(wp), intent(in)  :: u(:)
        real(wp), intent(out) :: dudt(:)
        type(stage_t), intent(in) :: stage
        real(wp) :: slope(nx), flux(nx), ul, ur
        integer  :: i, ip, im

        do i = 1, nx
            ip = modulo(i, nx) + 1
            im = modulo(i - 2, nx) + 1
            slope(i) = minmod(u(i) - u(im), u(ip) - u(i))
        end do
        do i = 1, nx                                        ! flux through the right face of cell i
            ip = modulo(i, nx) + 1
            ul = u(i) + slope(i) / 2
            ur = u(ip) - slope(ip) / 2
            flux(i) = (ul**2 / 2 + ur**2 / 2) / 2 - max(abs(ul), abs(ur)) * (ur - ul) / 2
        end do
        do i = 1, nx
            im = modulo(i - 2, nx) + 1
            dudt(i) = -(flux(i) - flux(im)) / dx
        end do
    end subroutine burgers

    pure real(wp) function minmod(a, b)
        real(wp), intent(in) :: a, b

        minmod = 0.0_wp
        if (a * b > 0) minmod = sign(min(abs(a), abs(b)), a)
    end function minmod

    pure real(wp) function total_variation(u)
        real(wp), intent(in) :: u(:)

        total_variation = sum(abs(u - cshift(u, 1)))
    end function total_variation

    !> One run to t_end at the given CFL (dt = CFL dx / max|u0|, the last step
    !> shortened to land on t_end): the largest TV increase over one step and
    !> the largest new extremum (above 1 or below -0.5). Optional output files.
    subroutine run(name, cfl, tv_up, new_ext, final, spacetime, worst, history, early)
        character(*), intent(in)  :: name
        real(wp),     intent(in)  :: cfl
        real(wp),     intent(out) :: tv_up, new_ext
        character(*), intent(in), optional :: final, spacetime
        character(*), intent(in), optional :: worst     !! (x, u) when the new extremum is largest, its t in the header
        character(*), intent(in), optional :: history   !! (t, new extremum at that step)
        character(*), intent(in), optional :: early     !! (x, u) after two steps

        type(rk_array_scheme_t) :: scheme
        real(wp) :: u(nx), time, dt, tv_old, tv_new, u_worst(nx), t_worst, ext_now
        integer  :: i, n, unit, next_frame, uh, n_early

        n = ceiling(t_end / (cfl * dx))
        dt = t_end / n
        call rk_scheme(scheme, name, nx, dt)
        u = merge(1.0_wp, -0.5_wp, x > 0.25_wp .and. x < 0.6_wp)
        time = 0.0_wp
        tv_up = 0.0_wp
        new_ext = 0.0_wp
        tv_old = total_variation(u)
        next_frame = 0
        u_worst = u
        t_worst = 0.0_wp
        if (present(history)) then
            open(newunit=uh, file=out // history, status="replace", action="write")
            write(uh, '(a)') "t overshoot"
        end if
        if (present(spacetime)) then
            open(newunit=unit, file=out // spacetime, status="replace", action="write")
            write(unit, '(a)') "x t u"
        end if
        do i = 1, n
            if (present(spacetime)) then
                if (time >= next_frame * t_end / 80 - 1.0e-12_wp) then
                    call write_frame(unit, time, u)
                    next_frame = next_frame + 1
                end if
            end if
            call scheme%advance(burgers, time, u)
            if (any(u /= u) .or. maxval(abs(u)) > 1.0e3_wp) then
                tv_up = huge(1.0_wp)
                new_ext = huge(1.0_wp)
                exit
            end if
            tv_new = total_variation(u)
            tv_up = max(tv_up, tv_new - tv_old)
            ext_now = max(maxval(u) - 1, -0.5_wp - minval(u), 0.0_wp)
            if (ext_now > new_ext) then
                u_worst = u
                t_worst = time
            end if
            new_ext = max(new_ext, ext_now)
            if (present(history)) write(uh, '(2es14.6)') time, max(ext_now, 1.0e-16_wp)
            if (present(early) .and. i == 2) then
                open(newunit=unit, file=out // early, status="replace", action="write")
                write(unit, '(a)') "x u"
                write(unit, '(2es14.6)') (x(n_early), u(n_early), n_early = 1, nx)
                close(unit)
            end if
            tv_old = tv_new
        end do
        if (present(spacetime)) then
            call write_frame(unit, time, u)
            close(unit)
        end if
        if (present(history)) close(uh)
        if (present(worst)) then
            open(newunit=unit, file=out // worst, status="replace", action="write")
            write(unit, '(a)') "x u"
            do i = 1, nx
                write(unit, '(2es14.6)') x(i), u_worst(i)
            end do
            close(unit)
        end if
        if (present(final)) then
            open(newunit=unit, file=out // final, status="replace", action="write")
            write(unit, '(a)') "x u"
            do i = 1, nx
                write(unit, '(2es14.6)') x(i), u(i)
            end do
            close(unit)
        end if
    end subroutine run

    subroutine write_frame(unit, time, u)
        integer,  intent(in) :: unit
        real(wp), intent(in) :: time, u(:)
        integer :: i

        do i = 1, nx, 2
            write(unit, '(3es14.6)') x(i), time, u(i)
        end do
        write(unit, '(a)') ""
    end subroutine write_frame

    function limit_text(c) result(s)
        real(wp), intent(in) :: c
        character(:), allocatable :: s
        character(len=8) :: buffer

        if (c == 0.0_wp) then
            s = "---"
        else
            write(buffer, '(f4.2)') c
            s = trim(buffer)
        end if
    end function limit_text

    subroutine write_parameters()
        integer :: u, j
        character(:), allocatable :: list

        list = ""
        do j = 1, size(names)
            if (j > 1) list = list // ", "
            list = list // "\texttt{" // trim(names(j)) // "}"
        end do
        open(newunit=u, file=out // "parameters.tex", status="replace", action="write")
        write(u, '(a)') "% Written by examples/burgers/main.f90: do not edit."
        write(u, '(a)') "\begin{tabular}{l>{\raggedright\arraybackslash}p{12.4cm}}"
        write(u, '(a)') "\toprule"
        write(u, '(a)') "domain & $0 \le x < 1$, periodic; $" // itoa(nx) // "$ cells, $\Delta x = 1/" // itoa(nx) // "$ \\"
        write(u, '(a)') "initial condition & $u(x, 0) = 1$ on $0.25 < x < 0.6$, $-0.5$ elsewhere (a square wave) \\"
        write(u, '(a)') "interval & $0 \le t \le " // fmt(t_end) // "$ \\"
        write(u, '(a)') "space & MUSCL, minmod limiter, Rusanov flux: TVD with forward Euler for CFL $\le 1/2$ \\"
        write(u, '(a)') "time step & $\Delta t = \mathrm{CFL}\,\Delta x / \max|u_0|$, rounded down to land on $t = " // &
                        fmt(t_end) // "$ \\"
        write(u, '(a)') "reference & \texttt{rk4}, CFL $= " // fmt(cfl_reference) // "$ \\"
        write(u, '(a)') "CFL sweep & $0.05, 0.10, \dots, " // fmt(0.05_wp * n_cfl) // &
                        "$; largest $\mathrm{TV}(u_{n+1}) - \mathrm{TV}(u_n)$ and largest new extremum " // &
                        "$\max(\max u - 1, -0.5 - \min u)$ over the run; TVD up to the first CFL with an increase " // &
                        "above $10^{" // &
                        itoa(nint(log10(tolerance))) // "}$ \\"
        write(u, '(a)') "figures & solutions at CFL $= " // fmt(cfl_show) // "$ \\"
        write(u, '(a)') "methods & " // list // " \\"
        write(u, '(a)') "interface & plain arrays, \texttt{rk\_array\_scheme\_t}, " // itoa(nx) // " unknowns \\"
        write(u, '(a)') "\bottomrule"
        write(u, '(a)') "\end{tabular}"
        close(u)
    end subroutine write_parameters

    function sci(x) result(s)
        real(wp), intent(in) :: x
        character(:), allocatable :: s
        character(len=32) :: mantissa
        integer :: ex

        if (x >= huge(1.0_wp)) then
            s = "blow-up"
            return
        end if
        if (x <= tolerance) then
            s = "$\le 10^{" // itoa(nint(log10(tolerance))) // "}$"
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

end program main
