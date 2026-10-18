![stdlib_Integrate: the word "Integrate", with the stability regions of Runge-Kutta methods of orders 1 to 4 behind its "I"](docs/tikz/stability_integrate.png)

# stdlib_Integrate

Time integration of initial-value problems, $dy/dt = f(t, y)$, $y(t_0) = y_0$, for Fortran, with explicit Runge-Kutta methods. Written as a candidate for the [Fortran standard library](https://github.com/fortran-lang/stdlib): **[proposal.md](proposal.md) is the full proposal**, and the place to start.

The program describes *what* is integrated; the library provides *how*. There are two ways to describe the problem, and they run the same methods with the same steps, bit for bit:

- **plain arrays**: the state is one rank-1 array, and the program writes only a right-hand side `f(t, y, dydt, stage)`;
- **model classes**: the program extends four abstract types (a model, its state, its tendencies, a scratch space), so that the state can be any data structure (several fields, multi-dimensional grids, halos, tiles), which the library never copies or reorganises. This is what models built on the method of lines need, and what existing Fortran ODE libraries do not offer.

```fortran
use stdlib_Integrate, only: rk_scheme, rk_array_scheme_t, stage_t

type(rk_array_scheme_t) :: scheme

call rk_scheme(scheme, "rk4", size(y), dt)     ! the method, by name
do it = 1, n_steps
    call scheme%advance(f, time, y)             ! y and time, one step forward
end do
```

## What is there

- three methods at fixed steps: Wicker-Skamarock RK3, classic RK4, Dormand-Prince 5(4);
- classic RK4 as a hardcoded step, and a generic engine that runs any explicit Butcher tableau, including one built by hand: every other method runs on it, and adding a method is adding a tableau;
- a hook to fill halos and apply boundary conditions before every evaluation.

Not yet: adaptive steps, low-storage forms, and other families of methods (multistep, implicit); the proposal says how they would fit.

Double precision only for now. Working and tested with fpm and gfortran.

## Building

With [fpm](https://fpm.fortran-lang.org):

```sh
fpm build
fpm test
fpm run --example decay     # 14 examples, in examples/
```

It depends only on [stdlib](https://github.com/fortran-lang/stdlib), which fpm fetches. To use it from another fpm project:

```toml
[dependencies]
stdlib_Integrate = { git = "https://github.com/ftucciarone/stdlib_Integrate" }
```

## Further reading

| | |
|---|---|
| [proposal.md](proposal.md) | the proposal: both routes from a first program, what a step does, the model contract, what would change for stdlib, open questions |
| `docs/usage.tex` | the user guide (LaTeX) |
| [examples/README.md](examples/README.md) | the examples: convergence, stability, invariants, stiff problems, PDEs |

## Licence

MIT, see [LICENSE](LICENSE).
