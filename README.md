# l1fit

Fast **constrained L1 curve fitting** for Python and NumPy. `l1fit` solves

```math
\min_{x \in \mathbb{R}^n} \; \lVert A x - b \rVert_1
\quad \text{s.t.} \quad C x = d, \;\; E x \le f
```

with `A` of size `k x n`, optional equality constraints `C x = d` (`l` rows) and optional inequality
constraints `E x <= f` (`m` rows), and optional sign restrictions on `x` and on the residuals. The L1 norm
makes the fit robust against outliers: where a least-squares line is pulled towards them, the L1 line
ignores them.

The solver is the modified simplex method of Barrodale and Roberts (ACM TOMS algorithm 552, `CL1`), written in
modern Fortran. It works on a single dense tableau, so it needs no LP modelling layer and is much faster than
general-purpose solvers on this problem class (see [Benchmark](#benchmark-against-cvxpy)). The Python package is
a thin, typed [nanobind](https://nanobind.readthedocs.io) binding that releases the GIL.

## Installation

`l1fit` is not on PyPI yet. Wheels for Linux (x86_64, aarch64), macOS (arm64, x86_64) and Windows (x86_64) are built by
`.github/workflows/wheels.yml`. One wheel per platform covers Python 3.12 to 3.15 (stable ABI). They need
nothing but NumPy: the Fortran code uses no runtime library, so there is no `libgfortran` to bundle.

To build from source you need a Fortran compiler (`gfortran`; on Windows the MinGW-w64 one, next to MSVC for the
C++), a C++17 compiler and CMake. [uv](https://docs.astral.sh/uv/) builds the extension with scikit-build-core:

```bash
uv sync     # development environment (builds the extension in editable mode)
uv build    # wheel and sdist in dist/
```

## Usage

```python
import numpy as np
from l1fit import Status, solve_l1

# Fit a line through 7 points, one of them an outlier of size 10: the L1 fit ignores it.
t = np.arange(7.0)
b = 1.0 + 2.0 * t
b[3] += 10.0
A = np.column_stack([np.ones(7), t])

result = solve_l1(A, b)
assert result.status is Status.OPTIMAL
print(result.x)  # [1. 2.]
print(result.objective)  # 10.0 (up to rounding) = the sum of the absolute residuals

# Optional equality constraints C x = d, inequality constraints E x <= f and sign restrictions:
result = solve_l1(
    A,
    b,
    C=np.array([[1.0, 1.0]]),
    d=np.array([4.0]),  # x0 + x1 = 4
    E=np.array([[0.0, 1.0]]),
    f=np.array([1.5]),  # x1 <= 1.5
    x_sign=np.array([0, 1]),  # x1 >= 0
)
```

`solve_l1(A, b, C, d, E, f, *, x_sign, residual_sign, tol, max_iter)` takes `float64` numpy arrays (other
dtypes raise `TypeError`) and never modifies them. (The inequality matrix is called `E`, as in the
paper.) It returns a frozen `L1Result` with `x`, `residual` (`b - A x`), `equality_residual`,
`inequality_slack`, `objective`, `status` and `iterations`. Infeasibility, rounding errors and reaching
`max_iter` are reported through `result.status` (`Status.INFEASIBLE`, `ROUNDING_ERRORS`,
`MAX_ITERATIONS`), not as exceptions; invalid arguments raise `TypeError` or `ValueError`. NaN and infinity
are rejected. The binding releases the GIL, and the solver is reentrant, so threads can solve independent
problems in parallel.

The solver works in double precision. The default tolerance is `1e-10` and the default iteration limit
`10 * (k + l + m)`.

## Benchmark against CVXPY

CVXPY is the usual way to solve this problem in Python without a dedicated routine. The benchmark
(`benchmarks/bench_cvxpy.py`, instances from `tests/instances`) uses the natural formulation:

```python
x = cp.Variable(n)
problem = cp.Problem(cp.Minimize(cp.norm(b - A @ x, 1)), [C @ x == d, E @ x <= f])
problem.solve(solver=...)
```

CVXPY rewrites it as a linear program with one epigraph variable per row (`|b_i - A_i x| <= t_i`: `n + k`
variables and `2k + m` inequality rows) and passes it to a solver. Three were used:

* **HiGHS dual simplex** (`solver="HIGHS"`, option `solver="simplex"`) - the same algorithm family as `l1fit`, on
  the larger LP;
* **HiGHS interior point** (`solver="HIGHS"`, option `solver="ipm"`);
* **Clarabel** (`solver="CLARABEL"`), an interior-point solver for conic problems and CVXPY's default for LPs.

Wall-clock time of the whole call, minimum of 20 runs for `l1fit` and of 3 for CVXPY (a new `cp.Problem` for every
run, so the times include CVXPY's canonicalization). The factor is the time relative to `l1fit`:

| Instance | Size `k x n` | Equality / inequality rows | l1fit | CVXPY + HiGHS simplex | CVXPY + HiGHS IPM | CVXPY + Clarabel |
|---|---|---|--:|--:|--:|--:|
| `500x10` | 500 x 10 | 0 / 0 | **0.5 ms** | 29.6 ms (58x) | 29.1 ms (57x) | 11.4 ms (22x) |
| `2000x20` | 2000 x 20 | 0 / 0 | **8.6 ms** | 591 ms (68x) | 403 ms (47x) | 90.0 ms (10x) |
| `2000x20_l3_m10` | 2000 x 20 | 3 / 10 | **9.8 ms** | 520 ms (53x) | 363 ms (37x) | 95.7 ms (10x) |
| `300x100` | 300 x 100 | 0 / 0 | **9.9 ms** | 166 ms (17x) | 170 ms (17x) | 218 ms (22x) |
| `5000x30` | 5000 x 30 | 0 / 0 | **60.5 ms** | 5283 ms (87x) | 2673 ms (44x) | 399 ms (7x) |
| `cvxpy_script_1000x250` | 1000 x 250 | 10 / 200 | **162 ms** | 4514 ms (28x) | 3025 ms (19x) | 2586 ms (16x) |

Apple M1 Pro, Python 3.12.13, NumPy 2.5.3, CVXPY 1.9.3 (HiGHS 1.15.1 through `highspy`, Clarabel 0.11.1), default
solver tolerances, `l1fit` with its defaults (`tol=1e-10`). Raw numbers, versions and per-solver details are in
`benchmarks/results/cvxpy.json`; reproduce with

```bash
uv run --group bench python benchmarks/bench_cvxpy.py --json benchmarks/results/cvxpy.json
```

`l1fit` is 7 to 87 times faster. Even when only the time that the solver itself reports is compared (leaving
out CVXPY's modelling overhead), it stays 5 to 86 times faster: the dedicated simplex works on a `k x n` tableau,
while the general LP is much larger. The optimal objectives agree to `1e-12` (HiGHS) and `1e-9` (Clarabel, whose
interior-point tolerance is looser). The comparison is between a specialised and a general-purpose solver on
dense random problems (five with a planted feasible point and Student-t noise, `cvxpy_script_1000x250` with pure
Gaussian data); the factors will differ on other data, and simplex methods
slow down on very large or degenerate problems. CVXPY, of course, solves far more than this one problem class.

## Background: a modernized CALGO 552

The algorithm is not new: this is a refactored and modernized version of `CALGO552.f`, the Fortran 77 routine
`CL1` published in 1980 by I. Barrodale and F. D. K. Roberts (*Algorithm 552: Solution of the constrained l1
linear approximation problem*, ACM Transactions on Mathematical Software 6(2), 231-235). The original is a
single subroutine with 18 arguments, about 60 numeric labels, `GOTO`s and single precision. The
copy this project started from also had a few garbled lines that did not compile; they were repaired first.
What was done:

* **Source form and types:** fixed form to free-form Fortran 2018 modules (`implicit none (type, external)`,
  `private` by default); single precision to double precision (one kind parameter `wp`, the accumulator `dp`).
* **Control flow:** every `GOTO`, arithmetic `IF` and label removed. The single routine was split step by step
  into small procedures (`pivot_tableau`, `select_entering_column`, `select_leaving_row`,
  `decide_at_optimum`, `set_up_phase1_costs`, ...), each verified before the next one.
* **Interface:** assumed-shape `contiguous` arrays, so four dimension arguments (`KLMD`, `KLM2D`, `NKLMD`,
  `N2D`) and the three workspace arrays (`CU`, `IU`, `S`) are gone: the workspace is allocated inside, and 11
  arguments are left. Arguments are
  validated, and failures are reported as status codes (`0` optimal, `1` infeasible, `2` rounding errors, `3`
  iteration limit, `4` invalid input, `5` allocation failed) instead of undefined behaviour. A `bind(C)` entry
  point (`l1_cl1`) serves the Python extension.
* **Behaviour kept:** the algorithm itself was not changed. Compiled in single precision, the refactored code
  reproduced the original bit for bit (status, iterations, solution, residuals) on about 2600 test problems, and it
  runs at the same speed. Double precision is the only numerical difference: it makes the solver reliable
  on problems where the single-precision original stops with rounding errors (for example a 1000 x 250 problem with
  200 inequality rows).
* **Packaging:** the Python extension, wheels for all major platforms, and a test suite that checks every answer
  against an independent HiGHS reference.

The history of the repository contains the original and the intermediate steps.

## Development

```bash
uv sync                                  # rebuilds when pyproject.toml, CMakeLists.txt or src/ change
uv run pytest tests                      # tests (about 6 s), see tests/README.md
uv run ruff check && uv run ruff format --check
uv run pyrefly check
uv run python -m nanobind.stubgen -m l1fit._core -o python/l1fit/_core.pyi   # after changing src/_core.cpp
```

| Path | Content |
|------|---------|
| `python/l1fit/` | The typed Python package: `solve_l1`, `L1Result`, `Status`, the generated stub `_core.pyi` |
| `src/` | The native code: the Fortran solver (`l1_calgo552.f90` algorithm, `l1_c_api.f90` C entry point, `l1_precision.f90`) and the nanobind binding `_core.cpp` (about 60 lines: validates shapes, releases the GIL, calls Fortran) |
| `benchmarks/` | The CVXPY comparison and its recorded results |
| `scripts/` | `check_installed_wheel.py`, run on every built wheel |
| `tests/` | The test suite and the benchmark instances (`tests/README.md`) |

Ruff and pyrefly check the Python package only (`python/`); Fortitude lints the Fortran (`src/`); see
`pyproject.toml` and `fortitude.toml`.

## Releasing

Pushing a tag `vX.Y.Z` that matches `version` in `pyproject.toml` publishes to PyPI. The `Wheels` workflow builds
the sdist and the five wheels, runs the tests against every wheel, and only then runs the `publish` job, which
checks the tag and the file list and uploads through
[trusted publishing](https://docs.pypi.org/trusted-publishers/) (no API token is stored).

```bash
# 1. bump `version` in pyproject.toml, commit and push to main, wait for green checks
# 2.
git tag v0.1.0 && git push origin v0.1.0
```

One-time set-up: on pypi.org add a (pending) publisher for the project `l1fit` with owner `jhelgert`, repository
`l1fit`, workflow `wheels.yml` and environment `pypi`, and create the environment `pypi` in the repository settings
(required reviewers there make each release a manual approval). A published version cannot be replaced; to try the
pipeline first, use the same set-up on test.pypi.org.
