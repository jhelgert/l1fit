# l1fit

A tiny Python package for blazingly fast curve fitting by an tailored primal
simplex algorithm and tailored dual simplex algorithms. In detail, this package
can solve the unconstrained/constrained optimization problems:

```math
\min_{x \in \mathbb{R}^n} \|Ax - b\|_1 \text{ s.t .} Cx = d, \;Dx \leq f, \\
```

with matrices $A \in \mathbb{R}^{m \times n}$, $C \in \mathbb{R}^{p \times n}$,
$D \in \mathbb{R}^{q \times n}$, $f \in \mathbb{R}^q$ and vectors
$b \in \mathbb{R}^m$, $d \in \mathbb{R}^p$ and $f \in \mathbb{R}^q$. Note
that the equality and inequality constraints are optional.

## Installation

The package contains a Fortran solver, so installing from source needs a Fortran compiler (`gfortran`),
a C++17 compiler and CMake. [uv](https://docs.astral.sh/uv/) builds it with scikit-build-core and nanobind:

```bash
uv sync                # development environment (builds the extension in editable mode)
uv build               # wheel (stable ABI: one wheel for Python >= 3.12) and sdist in dist/
```

Python 3.12, 3.13, 3.14 and 3.15 are supported. Wheels for Linux (manylinux) and macOS are built by `.github/workflows/wheels.yml`; Windows is planned.

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

## Development

```bash
uv sync                                  # rebuilds when pyproject.toml, CMakeLists.txt or src/ change
uv run pytest tests                      # tests (about 15 s), see tests/README.md
uv run ruff check && uv run ruff format --check
uv run pyrefly check
uv run python -m nanobind.stubgen -m l1fit._core -o python/l1fit/_core.pyi   # after changing src/_core.cpp
```

| Path | Content |
|------|---------|
| `python/l1fit/` | The typed Python package: `solve_l1`, `L1Result`, `Status`, the generated stub `_core.pyi` |
| `src/` | The native code: the Fortran solver (`l1_calgo552.f90` algorithm, `l1_c_api.f90` C entry point, `l1_precision.f90`) and the nanobind binding `_core.cpp` (about 60 lines: validates shapes, releases the GIL, calls Fortran) |
| `legacy/` | The original and the first modernized version, frozen as test oracles |
| `tests/` | Tests, benchmarks and the comparison with the original (`tests/README.md`) |

Ruff and pyrefly check the Python package only (`python/`); see `pyproject.toml`.

## TODO

- Create a mid-large set of benchmark problem instances, solve all of them via CVXPY
and sort them by size and solver status etc. such that we don't need to rely
on too simply and straightforward benchmark problems.
