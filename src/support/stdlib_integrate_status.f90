module stdlib_integrate_status
    !! Error status passed through the integrator's interfaces (models, states,
    !! time schemes), and what a time scheme does with a failure.
    !!
    !! A procedure that fails calls `status%raise(message)`; callers test
    !! `status%failed`. A time scheme's `advance` returns the status when the
    !! caller passed one; otherwise it stops the program with the failure's
    !! message (`fail_step`). Building a scheme (`rk_scheme`) follows the same
    !! rule.

    implicit none
    private

    type, public :: status_t
        !! Outcome of an operation: `failed` and, if so, a message.
        logical                       :: failed = .false.
        character(len=:), allocatable :: message
    contains
        procedure :: raise
    end type

    public :: type_mismatch
    public :: fail_step

contains

    !> Mark the status as failed, with `message`.
    subroutine raise(status, message)
        class(status_t), intent(inout) :: status
        character(*),    intent(in)    :: message

        status%failed = .true.
        status%message = message
    end subroutine raise

    !> Message for an argument of the wrong dynamic type, for `raise`.
    function type_mismatch(mod, proc, arg, typ) result(message)
        character(*), intent(in)      :: mod, proc, arg, typ
        character(len=:), allocatable :: message

        message = "type mismatch in " // trim(mod) // ": argument " // trim(arg) // " of " // &
                  trim(proc) // " should have type " // trim(typ)
    end function type_mismatch

    !> An operation (a time step, or building a scheme) failed with
    !> `local_status`: hand it to the caller's `status` if present, otherwise
    !> stop the program with its message. `procedure` names the failing
    !> procedure; `summary` is the message's first line after it (default
    !> "time step failed").
    subroutine fail_step(local_status, status, procedure, summary)
        type(status_t),           intent(in)              :: local_status
        type(status_t),           intent(inout), optional :: status
        character(*),             intent(in)              :: procedure
        character(*),             intent(in),    optional :: summary

        character(len=:), allocatable :: report

        if (present(status)) then
            status = local_status
            return
        end if
        if (present(summary)) then
            report = procedure // ": " // summary
        else
            report = procedure // ": time step failed"
        end if
        if (allocated(local_status%message)) report = report // new_line('a') // "  " // local_status%message
        error stop report
    end subroutine fail_step

end module stdlib_integrate_status
