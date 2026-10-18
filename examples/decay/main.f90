program main
    !! Convergence study on exponential decay, dy/dt = -lambda*y, y(0) = 1, to
    !! t = 1: the error of every method against exp(-lambda) as the time step is
    !! halved, and the observed order (slope of log error versus log dt) against
    !! the order the method has on linear problems. Uses the model classes
    !! (decay_model) and rk_scheme by name.
    !! Writes examples/decay/data/: parameters.tex and results.tex (tables),
    !! exact.dat and solution_<method>.dat (t, y at dt_plot), error_<method>.dat
    !! (dt, error), read by examples/decay/summary.tex. Run from the project root.

    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: status_t, time_scheme_c, rk_scheme, rk_method_names, rk_tableau, &
                                   & rk_tableau_t
    use decay_model, only: decay_model_t, decay_state_t
    implicit none

    character(*), parameter :: out = "examples/decay/data/"
    real(wp),     parameter :: lambda = 2.0_wp, t_end = 1.0_wp
    integer,      parameter :: n_levels = 8          !! dt = 1/10, 1/20, ..., 1/1280
    real(wp),     parameter :: floor = 1.0e-12_wp    !! below this, round-off pollutes the slope
    real(wp),     parameter :: dt_plot = 0.25_wp      !! step of the solution figure

    type(decay_model_t), target :: model
    character(len=8), allocatable :: names(:)
    type(rk_tableau_t) :: tableau
    real(wp) :: dt(n_levels), err(n_levels), observed
    integer  :: m, l, table, dat

    model%lambda = lambda
    names = rk_method_names()
    call execute_command_line("mkdir -p " // out)
    call write_parameters()

    ! The solutions at dt_plot against the exact one.
    open(newunit=dat, file=out // "exact.dat", status="replace", action="write")
    write(dat, '(a)') "t y"
    do l = 0, 200
        write(dat, '(2es14.6)') t_end * l / 200, exp(-lambda * t_end * l / 200)
    end do
    close(dat)
    do m = 1, size(names)
        open(newunit=dat, file=out // "solution_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat, '(a)') "t y"
        observed = run(trim(names(m)), dt_plot, nint(t_end / dt_plot), dat)
        close(dat)
    end do

    open(newunit=table, file=out // "results.tex", status="replace", action="write")
    write(table, '(a)') "% Written by examples/decay/main.f90: do not edit."
    write(table, '(a)') "\begin{tabular}{llccc}"
    write(table, '(a)') "\toprule"
    write(table, '(a)') "method & name & order & linear order & observed \\"
    write(table, '(a)') "\midrule"

    print '(a20, a10, 3a10)', "method", "name", "order", "expected", "observed"
    do m = 1, size(names)
        tableau = rk_tableau(trim(names(m)))
        open(newunit=dat, file=out // "error_" // trim(names(m)) // ".dat", status="replace", action="write")
        write(dat, '(a)') "dt error"
        do l = 1, n_levels
            dt(l) = t_end / (10 * 2**(l - 1))
            err(l) = abs(run(trim(names(m)), dt(l), 10 * 2**(l - 1)) - exp(-lambda * t_end))
            if (err(l) > 1.0e-15_wp) write(dat, '(2es14.6)') dt(l), err(l)
        end do
        close(dat)
        observed = slope(dt, err)
        print '(a20, a10, i10, i10, f10.3)', tableau%title, trim(names(m)), tableau%order, tableau%linear_order, observed
        write(table, '(a, " & \texttt{", a, "} & ", i0, " & ", i0, " & ", f5.2, " \\")') &
            tableau%title, trim(names(m)), tableau%order, tableau%linear_order, observed
    end do

    write(table, '(a)') "\bottomrule"
    write(table, '(a)') "\end{tabular}"
    close(table)

contains

    !> y(t_end) after n_steps steps of size dt with the method called `name`;
    !> with `unit`, every (t, y) is written to it.
    real(wp) function run(name, dt, n_steps, unit) result(y)
        character(*), intent(in) :: name
        real(wp),     intent(in) :: dt
        integer,      intent(in) :: n_steps
        integer,      intent(in), optional :: unit

        class(time_scheme_c), allocatable, target :: scheme
        type(decay_state_t), target :: state
        type(status_t) :: status
        real(wp) :: time
        integer  :: i

        call rk_scheme(scheme, name, dt, model)
        state%y = 1.0_wp
        time = 0.0_wp
        if (present(unit)) write(unit, '(2es14.6)') time, state%y
        do i = 1, n_steps
            call scheme%advance(time, state, status)
            if (status%failed) error stop status%message
            if (present(unit)) write(unit, '(2es14.6)') time, state%y
        end do
        y = state%y
    end function run

    !> Slope of log(err) against log(dt) between the two finest levels whose
    !> errors are both above the round-off floor.
    real(wp) function slope(dt, err) result(p)
        real(wp), intent(in) :: dt(:), err(:)
        integer :: l

        p = 0.0_wp
        do l = size(dt), 2, -1
            if (err(l) > floor .and. err(l - 1) > floor) then
                p = log(err(l - 1) / err(l)) / log(dt(l - 1) / dt(l))
                return
            end if
        end do
    end function slope

    !> The parameters of this run, as a LaTeX table.
    subroutine write_parameters()
        integer :: u, k
        character(len=64) :: buffer
        character(:), allocatable :: list

        list = ""
        do k = 1, size(names)
            if (k > 1) list = list // ", "
            list = list // "\texttt{" // trim(names(k)) // "}"
        end do
        open(newunit=u, file=out // "parameters.tex", status="replace", action="write")
        write(u, '(a)') "% Written by examples/decay/main.f90: do not edit."
        write(u, '(a)') "\begin{tabular}{ll}"
        write(u, '(a)') "\toprule"
        write(buffer, '(f0.2)') lambda
        write(u, '(a)') "decay rate & $\lambda = " // trim(buffer) // "$ \\"
        write(u, '(a)') "initial condition & $y(0) = 1$ \\"
        write(buffer, '(f0.2)') t_end
        write(u, '(a)') "interval & $0 \le t \le " // trim(buffer) // "$ \\"
        write(buffer, '(i0)') 10 * 2**(n_levels - 1)
        write(u, '(a)') "steps, convergence & $\Delta t = 1/10, 1/20, \dots, 1/" // trim(buffer) // "$ (" // &
                        itoa(n_levels) // " levels) \\"
        write(buffer, '(f4.2)') dt_plot
        write(u, '(a)') "step, solution figure & $\Delta t = " // trim(buffer) // "$ (" // &
                        itoa(nint(t_end / dt_plot)) // " steps) \\"
        write(u, '(a)') "methods & " // list // " \\"
        write(u, '(a)') "interface & model classes, \texttt{rk\_scheme(scheme, name, dt, model)} \\"
        write(u, '(a)') "observed order & slope of the two finest levels with errors above $10^{" // &
                        itoa(nint(log10(floor))) // "}$ \\"
        write(u, '(a)') "\bottomrule"
        write(u, '(a)') "\end{tabular}"
        close(u)
    end subroutine write_parameters

    function itoa(i) result(s)
        integer, intent(in) :: i
        character(:), allocatable :: s
        character(len=16) :: buffer

        write(buffer, '(i0)') i
        s = trim(buffer)
    end function itoa

end program main
