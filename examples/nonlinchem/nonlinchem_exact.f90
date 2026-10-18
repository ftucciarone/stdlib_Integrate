module nonlinchem_exact
    !! Exact solution of the nonlinear chemistry problem (SciML
    !! prob_ode_nonlinchem), y1' = -y1, y2' = y1 - y2^2, y3' = y2^2, y = (1, 0, 0):
    !!   y1 = exp(-t);
    !!   y2 = u'/u with u'' = exp(-t) u, which in s = 2 exp(-t/2) is the modified
    !!        Bessel equation of order 0, so u = I0(s) + c K0(s) and
    !!        y2 = -(s/2) (I1(s) - c K1(s)) / (I0(s) + c K0(s)), c = I1(2)/K1(2)
    !!        (from y2(0) = 0);
    !!   y3 = 1 - y1 - y2 (y1 + y2 + y3 is constant).
    !! The modified Bessel functions are summed from their power series, exact to
    !! round-off for 0 < s <= 2 (Abramowitz & Stegun 9.6.10-11).
    use stdlib_kinds, only: wp => dp
    implicit none
    private

    public :: exact_solution, bessel_i0, bessel_i1, bessel_k0, bessel_k1

    real(wp), parameter :: euler_gamma = 0.57721566490153286060651209_wp
    integer,  parameter :: n_terms = 30

contains

    !> (y1, y2, y3) at time t >= 0.
    function exact_solution(t) result(y)
        real(wp), intent(in) :: t
        real(wp) :: y(3), s, c

        s = 2 * exp(-t / 2)
        c = bessel_i1(2.0_wp) / bessel_k1(2.0_wp)
        y(1) = exp(-t)
        y(2) = -(s / 2) * (bessel_i1(s) - c * bessel_k1(s)) / (bessel_i0(s) + c * bessel_k0(s))
        y(3) = 1 - y(1) - y(2)
    end function exact_solution

    !> I0(x) = sum_k (x^2/4)^k / (k!)^2
    pure real(wp) function bessel_i0(x) result(r)
        real(wp), intent(in) :: x
        real(wp) :: term
        integer  :: k

        term = 1.0_wp
        r = term
        do k = 1, n_terms
            term = term * (x * x / 4) / (k * k)
            r = r + term
        end do
    end function bessel_i0

    !> I1(x) = (x/2) sum_k (x^2/4)^k / (k! (k+1)!)
    pure real(wp) function bessel_i1(x) result(r)
        real(wp), intent(in) :: x
        real(wp) :: term
        integer  :: k

        term = x / 2
        r = term
        do k = 1, n_terms
            term = term * (x * x / 4) / (k * (k + 1))
            r = r + term
        end do
    end function bessel_i1

    !> K0(x) = -(ln(x/2) + gamma) I0(x) + sum_{k>=1} H_k (x^2/4)^k / (k!)^2
    pure real(wp) function bessel_k0(x) result(r)
        real(wp), intent(in) :: x
        real(wp) :: term, harmonic
        integer  :: k

        r = -(log(x / 2) + euler_gamma) * bessel_i0(x)
        term = 1.0_wp
        harmonic = 0.0_wp
        do k = 1, n_terms
            term = term * (x * x / 4) / (k * k)
            harmonic = harmonic + 1.0_wp / k
            r = r + harmonic * term
        end do
    end function bessel_k0

    !> K1(x) = 1/x + ln(x/2) I1(x)
    !>         - (x/4) sum_k (psi(k+1) + psi(k+2)) (x^2/4)^k / (k! (k+1)!),
    !> psi(n) = -gamma + H_{n-1}.
    pure real(wp) function bessel_k1(x) result(r)
        real(wp), intent(in) :: x
        real(wp) :: term, h_k, sum
        integer  :: k

        term = 1.0_wp
        h_k = 0.0_wp
        sum = (-euler_gamma) + (-euler_gamma + 1.0_wp)
        do k = 1, n_terms
            term = term * (x * x / 4) / (k * (k + 1))
            h_k = h_k + 1.0_wp / k
            sum = sum + ((-euler_gamma + h_k) + (-euler_gamma + h_k + 1.0_wp / (k + 1))) * term
        end do
        r = 1 / x + log(x / 2) * bessel_i1(x) - (x / 4) * sum
    end function bessel_k1

end module nonlinchem_exact
