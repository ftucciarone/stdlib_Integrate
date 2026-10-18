"""Boundaries of the linear stability regions |R(z)| <= 1 of explicit RK methods.

For an s-stage method of order s (s <= 4), R(z) = 1 + z + z^2/2! + ... + z^s/s!.
The boundary is the root locus { z : R(z) = exp(i*theta), 0 <= theta <= 2*pi }:
the s roots are found once (Durand-Kerner) and followed by Newton steps as theta
grows, one branch per root. Writes stability_rk<s>.dat (columns: x y; a blank
line separates branches) for pgfplots. Pure Python, no dependencies.

    python3 stability_regions.py
"""
import cmath
from math import factorial, pi

STEPS = 2000


def poly(s, z):
    """R(z) and R'(z) for the degree-s truncated exponential."""
    r = sum(z**k / factorial(k) for k in range(s + 1))
    dr = sum(z**(k - 1) / factorial(k - 1) for k in range(1, s + 1))
    return r, dr


def initial_roots(s, w):
    """All roots of R(z) - w = 0 (Durand-Kerner on the monic polynomial)."""
    monic = lambda z: (poly(s, z)[0] - w) * factorial(s)
    roots = [(0.4 + 0.9j) ** k for k in range(s)]
    for _ in range(500):
        new = []
        for i, z in enumerate(roots):
            den = 1
            for j, other in enumerate(roots):
                if j != i:
                    den *= z - other
            new.append(z - monic(z) / den)
        roots = new
    return roots


def newton(s, z, w):
    for _ in range(50):
        r, dr = poly(s, z)
        step = (r - w) / dr
        z -= step
        if abs(step) < 1e-14:
            break
    return z


for s in range(1, 5):
    branches = [[z] for z in initial_roots(s, 1.0)]
    for n in range(1, STEPS + 1):
        w = cmath.exp(2j * pi * n / STEPS)
        for branch in branches:
            branch.append(newton(s, branch[-1], w))
    with open(f"stability_rk{s}.dat", "w") as out:
        for branch in branches:
            for z in branch[::4] + [branch[-1]]:
                out.write(f"{z.real:.5f} {z.imag:.5f}\n")
            out.write("\n")
