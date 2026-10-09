# Tests

The tests check the solver (`src/*.f90`) and the Python package (`l1fit`) against *independent* answers: a
HiGHS reference optimum, feasibility conditions and closed-form solutions. They do not compare with another
implementation of the algorithm.

## Layout

| Path | Purpose |
|------|---------|
| `instances.py` | Instance definitions (seeded) and `.npz` load/save. |
| `instances/*.npz`, `instances/manifest.json` | The 27 persisted instances, each with its HiGHS reference optimum (`ref_*`). |
| `generate_instances.py` | Regenerates the instance files and the manifest (only needed when `instances.py` changes). |
| `solvers.py` | The `Result` type, the adapter registry and the hang guard (`call_with_timeout`). |
| `fortran_src.py` | `ctypes` adapter of the Fortran sources through the C entry point `l1_cl1` (`src`), and the bounds-checked build (`SrcCheckedCL1`). Compiled with `gfortran` into `tests/_build/`. |
| `extension_adapter.py` | Adapter of the installed package: runs the instances through `l1fit.solve_l1` (`extension`). |
| `checks.py` | Solver-independent verification of a solution. |
| `test_l1_instances.py` | Correctness, closed-form and metamorphic tests, for both adapters on every instance. |
| `test_extension_api.py` | The Python interface and the contract of the binding (argument validation, copies, layouts, threads, the type stub). |
| `fortran/check_cl1_interface.f90`, `test_fortran_interface.py` | A Fortran program that checks the argument validation of `cl1` and oversized arrays (not reachable through the C interface); built with `-std=f2018 -Wall -Wextra -fcheck=all`. |
| `test_checked_build.py`, `run_checked_build.py` | All instances against `src/` built with `-fcheck=all -O0`, in a child process (a violated run-time check aborts the process). |
| `test_hang_guard.py` | Tests of the protection against solver calls that never return. |

The tests need `gfortran` (they compile the Fortran sources) and the installed package (`uv sync`).

## Commands (from the repository root)

```bash
uv sync                                  # builds the extension in editable mode; rebuilds when src/ changes
uv run pytest tests -q                   # about 6 s

uv run tests/generate_instances.py       # regenerate the instances; review the manifest diff

# strict compile and lint of the sources
uv run fortitude check --preview src/
gfortran -std=f2018 -Wall -Wextra -Wconversion -fcheck=all -g -c src/l1_precision.f90 src/l1_calgo552.f90 src/l1_c_api.f90
```

The `-Wcompare-reals` warnings of `l1_calgo552.f90` are deliberate: the algorithm compares floats for
(in)equality (`Q(KLM1,IN) /= XMAX`, the sign flags of `x` and `res`).

## Instances

* **Closed-form**: `median_1d`, `line_with_outlier`, `tiny_1x1`, `zero_rhs`, `square_exact`.
* **Edge cases**: `underdetermined_2x5` (non-unique x), `degenerate_duplicate_rows` (ties),
  `eq_fully_determined`, `infeasible_bounds` / `infeasible_equalities` (status 1), `iteration_limit`
  (status 3).
* **Sign restrictions (`KODE=1`)**: `kode1_x_signs`, `kode1_residual_signs`, `kode1_signs_and_constraints`.
* **Random**: heavy-tailed (Student-t, df=2) and Gaussian noise, with equality and inequality rows.
* **Benchmarks** (tag `benchmark`): 500x10 up to 5000x30, 300x100 (many columns) and
  `bench_cvxpy_script_1000x250`. They are also the instances of `benchmarks/bench_cvxpy.py`.

## What is asserted

1. `test_solution_is_correct` - for every adapter and instance: the expected status, the reported objective
   equals `||Ax-b||_1` recomputed from `x`, equals the HiGHS optimum (relative `1e-8`), `x` is feasible
   (equalities, inequalities, sign restrictions), the residual vector is consistent, and `x` equals the
   reference for unique optima.
2. Closed-form results, the iteration limit, and metamorphic properties (the optimum does not change when
   rows are permuted; scaling the whole problem by `s > 0` keeps `x` and scales the objective by `s`).
3. The Python interface (`test_extension_api.py`): validation of every argument, optional constraints and
   signs, inputs never modified, memory layouts, threads (the GIL is released: the test measures the
   speed-up), the generated type stub.
4. The Fortran interface: each invalid argument is rejected with `L1_INVALID_INPUT` and leaves `q`
   untouched; arrays larger than required give identical results; no hidden array temporaries.
5. The bounds-checked build solves all instances.

## Hang protection

A bug can make the solver loop forever, for example in the optimality test. An iteration limit cannot stop
that: `iter` is only incremented by the pivot step, so a loop that never pivots never counts. Instead:

* every adapter runs its foreign call in a worker thread with a wall-clock limit
  (`solvers.call_with_timeout`, default 30 s, `L1FIT_SOLVE_TIMEOUT` to change it; the slowest legitimate call
  takes about 1.5 s). A stuck call cannot be killed, so its thread is abandoned (a daemon thread) and
  `SolverHang` is raised, which is an ordinary test failure;
* after `MAX_HANGS` (3) hangs an adapter fails every further call immediately;
* the two tests that run subprocesses have their own limits (`test_fortran_interface.py`: 60 s,
  `test_checked_build.py`: 300 s);
* `test_hang_guard.py` tests the guard, including a call that blocks inside C.

## Adding another implementation

Write an adapter with `name`, `precision = "double"` and `solve(inst) -> Result`, and register it in
`solvers.available_solvers()`. Every correctness test then runs against it.

## Notes

* `uv sync` does not rebuild an editable install when only files in `src/` change, unless told to
  (`tool.uv.cache-keys` in `pyproject.toml`). Check that an experiment really runs the changed code.
* The Fortran code needs no runtime library: `L1FIT_FORTRAN_RUNTIME` in `CMakeLists.txt` is `NONE`
  (default), `SHARED` or `STATIC`, and the default build links neither `libgfortran` nor `libquadmath`. This
  is guarded at build time (`cmake/check_no_fortran_runtime.cmake`, which inspects the object files) and on
  the installed wheel (`scripts/check_installed_wheel.py`). It is what allows the Windows wheel to link
  MinGW-compiled Fortran into a module built with MSVC.
