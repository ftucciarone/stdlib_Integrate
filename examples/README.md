# Examples

Runnable problems that check the methods against mathematics and physics, not
the code against its specification. `test/` asks "does the library do what it
says?" (buffers, failures, bit-for-bit equivalences); an example asks "does the
method behave the way the theory says it should?": convergence at the
predicted order, conservation of an invariant, the stability limit of the
stability region, the right phase or amplitude error.

## Conventions

- One folder per example, `examples/<name>/`:
  - `main.f90`, the program, and its model or helper modules;
  - `summary.tex`, a short report in three parts:
    1. **Problem**: the equations, initial condition and exact solution as
       displayed formulas, a sentence on what makes the problem interesting,
       and the table of every parameter the run used;
    2. **Solutions**: the numerical solutions of the methods plotted against
       the exact (or reference) solution, with the parameters of the plot;
    3. **Test**: what the theory predicts, the measured result (table and
       figure), and a short reading of it;
  - `data/`, everything the program writes: `parameters.tex` (the parameter
    table, written from the program's own constants so it cannot disagree
    with the run), `results.tex` (the result table) and `.dat` files (the
    curves). The program creates the folder; never edit these files, rerun
    the example.
- Examples share no code with each other or with `test/`: a module two of them
  need is copied (as `test/test_decay_model.f90` is from `decay/`), so each
  one can evolve freely. fpm builds all examples together, so **module names
  must be unique across examples**: prefix them with the example's name when
  they could clash (`oscillator_stability`).
- Each example is registered in `fpm.toml`:

  ```toml
  [[example]]
  name = "<name>"
  source-dir = "examples/<name>"
  main = "main.f90"
  ```

- Run it from the project root, then build its summary in its folder:

  ```sh
  fpm run --example <name>
  cd examples/<name> && pdflatex summary.tex
  ```

- What the theory predicts is computed independently of the step routines
  where possible (for example the stability polynomial from the Butcher
  tableau), so the comparison checks the steps against the mathematics.

## Catalogue

| Example | API | Question | Theory |
|---|---|---|---|
| `decay` | classes | Error at `t = 1` on `y' = -2y` as `dt` halves, every method | Slope = order on linear problems (rk3ws: 3, not 2) |
| `stability` | arrays | Largest stable `dt` on `y' = -y` | First `h` with `\|R(-h)\| > 1`: 2.51 (rk3ws), 2.79 (rk4), 3.31 (dp54) |
| `oscillator` | arrays | Energy and phase after 100 periods of `x'' = -x`; imaginary-axis limit | `\|R(ih)\|^(2n)`, `n (arg R(ih) - h)`; limits sqrt(3) (rk3ws), 2 sqrt(2) (rk4), 0.997 (dp54) |
| `nonlinchem` | arrays | SciML `prob_ode_nonlinchem`: error at `t = 20` as `dt` halves; `y1 + y2 + y3` | General order (rk3ws: 2, not 3); linear invariant kept to round-off. Exact solution via Bessel functions; `dp54_mpmath.py` reruns dp54 in 40 digits |
| `lotkavolterra` | arrays | SciML `prob_ode_lotkavolterra`, `(a,b,c,d) = (1.5,1,3,1)`: drift of `V = dx - c ln x + by - a ln y` as `dt` halves and over `t <= 1000` | Drift at the general order (rk3ws: 2); grows linearly in time, never stops; exact orbit = level curve `V = V0` |
| `kepler` | arrays | Ellipse `e = 0.6`: error after 3 periods as the steps per period double; energy and angular momentum over 100 periods | General order (rk3ws: 2); H and L drift; exact solution from Kepler's equation |
| `arenstorf` | arrays | SciML `prob_ode_threebody`: closing error after one period, 1000 to 4096000 steps | General order only once the Moon passages are resolved; cost to close to 1e-3: dp54 192k evaluations, rk4 512k, rk3ws more than 4096000 steps |
| `pleiades` | arrays | SciML `prob_ode_pleiades`, 7 bodies: error at t = 1 and t = 3; momentum, angular momentum, energy | Linear momentum kept to round-off by every method; rk3ws order 2 before the close encounters, apparently 3 after (effective order 3) |
| `vanderpol` | arrays | SciML `prob_ode_vanderpol`: error, limit-cycle period, stable step as mu grows to 1000 | Period 6.6632868593 (literature value matched); stable step ~ h*/(3 mu) for rk3ws, which fails on the slow branch; earlier for rk4 and dp54, which fail in the fast jumps |
| `lorenz` | arrays | Classic Lorenz system: error at t = 1; time to depart from the reference vs dt; Lyapunov exponent (Benettin) | lambda = 0.909; T_div grows like (p/lambda) ln(1/dt) in steps; rk3ws needs p = 3 (effective order) |
| `burgers` | arrays | Burgers, square wave, MUSCL-minmod-Rusanov: TV increase and new extrema vs CFL | No method of the library is SSP (C = 0): TVD lost at CFL 0.30 (dp54), 0.45 (rk4), 0.50 (rk3ws), inside the range an SSP method is guaranteed (0.5) |
| `heat` | classes (`apply_boundaries` for the boundary, pieces) | Heat equation by the method of lines: time error, stable step vs dx, PDE error | Stable step = h*/\|lambda_n\| within 0.25 %; PDE error O(dx^2) for every method |
| `brusselator` | arrays | SciML `prob_ode_brusselator_1d`, reaction-diffusion: error, stable step vs dx | Stable step -> h* dx^2/(4 alpha) within 0.4 % from 64 points |
| `robertson` | arrays | SciML `prob_ode_rober`, stiff: failure time vs dt, largest step reaching T | Failure at \|lambda_f\| dt = h*; largest step within ~1 % at T = 100, 1000; 2.0-2.7 million steps to t = 1000 |

## Building the summaries

The summaries with large figures (pleiades, vanderpol, lorenz, burgers, heat,
brusselator, robertson) need `lualatex` (pdfLaTeX runs out of memory); the
others build with either. Some examples also check their reference against an
independent computation: `nonlinchem/dp54_mpmath.py` reruns dp54 in 40 digits.

## Possible next

| Example | Question | Theory |
|---|---|---|
| forced decay | `y' = -y + cos t` | full order only with correct stage times |
| rk3ws, processed | rk3ws with pre- and post-processing on every example where it shows order 2 | effective order 3 (shown on `nonlinchem`) |
| adaptive steps | dp54 with its embedded error estimate on `arenstorf` | far fewer evaluations than the fixed-step cost table |
