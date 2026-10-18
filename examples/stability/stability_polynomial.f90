module stability_polynomial
    !! The stability function of an explicit Runge-Kutta method, from its
    !! tableau: one step on y' = lambda*y is y <- R(z) y with z = lambda*dt and
    !!   R(z) = 1 + z b^T (I - zA)^{-1} 1,
    !! computed by forward substitution (A is strictly lower triangular):
    !!   g(i) = 1 + z sum_{j<i} A(i,j) g(j),   R = 1 + z sum_i b(i) g(i).
    !! Independent of the step routines: the theory the runs are compared with.
    use stdlib_kinds, only: wp => dp
    use stdlib_Integrate, only: rk_tableau_t
    implicit none
    private

    public :: stability_r

contains

    complex(wp) function stability_r(tableau, z) result(r)
        type(rk_tableau_t), intent(in) :: tableau
        complex(wp),        intent(in) :: z

        complex(wp) :: g(tableau%stages)
        integer :: i

        do i = 1, tableau%stages
            g(i) = 1.0_wp + z * sum(tableau%A(i, 1:i-1) * g(1:i-1))
        end do
        r = 1.0_wp + z * sum(tableau%b * g)
    end function stability_r

end module stability_polynomial
