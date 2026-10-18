"""dp54 on the nonlinear chemistry problem in 40-digit arithmetic (mpmath), to
see its convergence below double-precision round-off. Its tableau is the one
of src/schemes/rk/stdlib_integrate_rk_tableaus.f90, as exact fractions; it
integrates to
t = 20 with dt = 0.5, 0.25, ..., 0.5/256 and writes data/dp54_mpmath.dat
(dt, error), plotted by summary.tex. Run from examples/nonlinchem/ (a few
seconds)."""
import pathlib

import mpmath as mp

mp.mp.dps = 40
here = pathlib.Path(__file__).resolve().parent

# Dormand-Prince 5(4): row i of A lists A(i, 1:i-1).
dp54 = dict(
    A=[[],
       ['1/5'],
       ['3/40', '9/40'],
       ['44/45', '-56/15', '32/9'],
       ['19372/6561', '-25360/2187', '64448/6561', '-212/729'],
       ['9017/3168', '-355/33', '46732/5247', '49/176', '-5103/18656'],
       ['35/384', '0', '500/1113', '125/192', '-2187/6784', '11/84']],
    b=['35/384', '0', '500/1113', '125/192', '-2187/6784', '11/84', '0'])


def coef(s):
    return mp.mpf(int(s.split("/")[0])) / int(s.split("/")[1]) if "/" in s else mp.mpf(s)


A = [[coef(x) for x in row] for row in dp54["A"]]
b = [coef(x) for x in dp54["b"]]
c_bessel = mp.besseli(1, 2) / mp.besselk(1, 2)


def exact(t):
    s = 2 * mp.e ** (-t / 2)
    y1 = mp.e ** (-t)
    y2 = -(s / 2) * (mp.besseli(1, s) - c_bessel * mp.besselk(1, s)) / (mp.besseli(0, s) + c_bessel * mp.besselk(0, s))
    return [y1, y2, 1 - y1 - y2]


def f(y):
    return [-y[0], y[0] - y[1] ** 2, y[1] ** 2]


def error(dt, t_end=20):
    y = [mp.mpf(1), mp.mpf(0), mp.mpf(0)]
    for _ in range(int(t_end / dt)):
        k = []
        for i in range(len(b)):
            k.append(f([y[j] + dt * sum(A[i][m] * k[m][j] for m in range(len(A[i]))) for j in range(3)]))
        y = [y[j] + dt * sum(b[i] * k[i][j] for i in range(len(b))) for j in range(3)]
    return max(abs(a - e) for a, e in zip(y, exact(t_end)))


(here / "data").mkdir(exist_ok=True)
with open(here / "data/dp54_mpmath.dat", "w") as out:
    out.write("dt error\n")
    for level in range(9):
        dt = mp.mpf(0.5) / 2 ** level
        e = error(dt)
        out.write(f"{mp.nstr(dt, 8)} {mp.nstr(e, 8)}\n")
        print(mp.nstr(dt, 6), mp.nstr(e, 5))
