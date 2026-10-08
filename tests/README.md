# Tests and benchmarks for the L1 solver

These tests pin down the behaviour of `CL1` (ACM TOMS algorithm 552, Barrodale & Roberts) so the
Fortran can be modernized and wrapped as a CPython extension without silently changing results.

## Source layout

| Path | Purpose |
|------|---------|
| `../src/l1_precision.f90` | Working precision `wp` (currently `c_double`) and accumulator precision `dp`. |
| `../src/l1_calgo552.f90` | Fortran interface `cl1(k,l,m,n,q,kode,toler,iter,x,res,error)`: assumed-shape `contiguous` arrays, validated arguments, local workspace. Status constants `L1_OPTIMAL ... L1_ALLOC_FAILED`. |
| `../src/l1_c_api.f90` | `bind(C)` wrapper `l1_cl1` with the same 11 arguments (explicit-shape arrays derived from `k,l,m,n`) for `ctypes` / the future CPython extension. |
| `../legacy/CALGO552.f` | **Frozen oracle #1**: the repaired original (fixed-form, single precision, 18-argument interface). Never refactor this file. |
| `../legacy/f90_double/` | **Frozen oracle #2**: the first free-form double-precision version (module + `cl1_` with the 18-argument interface), frozen before the interface was changed. Reference for bit-for-bit checks. Excluded from linting in `fortitude.toml`. |
| `../legacy/CALGO552_as_received.f90.txt` | The file exactly as received (9 garbled lines, did not compile). |

## Test layout

| Path | Purpose |
|------|---------|
| `instances.py` | Instance definitions (seeded) + `.npz` load/save. |
| `instances/*.npz`, `instances/manifest.json` | Persisted instances with HiGHS reference optimum (`ref_*`) and legacy golden output (`legacy_*`). |
| `generate_instances.py` | Regenerates the instance files and the manifest. |
| `solvers.py`, `legacy_cl1.py` | Solver adapters sharing one interface (`name`, `precision`, `is_oracle`, `matches_legacy`, `solve(inst) -> Result`). See below. |
| `checks.py` | Solver-independent solution verification (objective, feasibility, residual vector, x). |
| `test_l1_instances.py` | Correctness, closed-form, metamorphic and golden tests per adapter. |
| `test_differential.py` | Bit-for-bit differential fuzzing against the frozen oracles. |
| `fortran/check_cl1_interface.f90`, `test_fortran_interface.py` | Fortran program checking argument validation and oversized arrays (not reachable through the C ABI); built strictly with `-fcheck=all` and run by pytest. |
| `test_checked_build.py`, `run_checked_build.py` | All instances and fuzz cases against `src/` built with `-fcheck=all -O0`, in a child process (a violated run-time check aborts the whole process). |
| `benchmark.py`, `benchmarks/legacy_baseline.json` | Timing runner and the recorded baseline of the legacy oracle (single precision). |

### Adapters

| Name | Sources | C ABI | Precision | Role |
|------|---------|-------|-----------|------|
| `legacy` | `legacy/CALGO552.f` | `cl1_` (18 args) | single | oracle; golden values come from it |
| `legacy_double` | `legacy/f90_double/` | `cl1_` (18 args) | double | oracle for double-precision refactors |
| `src` | `src/*.f90` | `l1_cl1` (11 args) | read from `l1_precision.f90` | the code under development |
| `SrcCheckedCL1` | `src/*.f90`, `-fcheck=all -O0` | `l1_cl1` | idem | only used in the child process of `test_checked_build.py` |

All adapters call the library through `ctypes`; each is compiled with `gfortran` into
`tests/_build/` (cache keyed by source content and flags). The `src` adapter reads `wp` from
`l1_precision.f90`, so it uses float32 or float64 arrays and the matching default `TOLER`
(2e-5 single, 1e-10 double) automatically.

## Commands (from the repository root)

```bash
# all tests (about 10 s)
uv run --no-project --with numpy --with scipy --with pytest pytest tests -q

# benchmark (+ compare with the recorded baseline, same machine only)
uv run tests/benchmark.py --solver src --compare tests/benchmarks/legacy_baseline.json

# regenerate instances (only when instances.py changes; review the manifest diff)
uv run tests/generate_instances.py

# strict compile and lint of the sources
mkdir -p /tmp/l1mod
gfortran -std=f2018 -Wall -Wextra -Wconversion -fcheck=all -g -J /tmp/l1mod -I /tmp/l1mod -c src/l1_precision.f90 -o /tmp/l1mod/p.o
gfortran -std=f2018 -Wall -Wextra -Wconversion -fcheck=all -g -J /tmp/l1mod -I /tmp/l1mod -c src/l1_calgo552.f90 -o /tmp/l1mod/c.o
gfortran -std=f2018 -Wall -Wextra -Wconversion -fcheck=all -g -J /tmp/l1mod -I /tmp/l1mod -c src/l1_c_api.f90 -o /tmp/l1mod/a.o
fortitude check --preview src/
```

The three `-Wcompare-reals` warnings of `l1_calgo552.f90` are deliberate: the algorithm compares floats
for (in)equality (`Q(KLM1,IN) /= XMAX`, sign flags of `x`/`res`).

## Instances (27)

* **Closed-form**: `median_1d`, `line_with_outlier`, `tiny_1x1`, `zero_rhs`, `square_exact`.
* **Edge cases**: `underdetermined_2x5` (non-unique x), `degenerate_duplicate_rows` (ties),
  `eq_fully_determined`, `infeasible_bounds` / `infeasible_equalities` (`KODE=1`),
  `iteration_limit` (`KODE=3`).
* **`KODE=1` sign restrictions**: `kode1_x_signs`, `kode1_residual_signs`, `kode1_signs_and_constraints`.
* **Random**: heavy-tailed (Student-t, df=2) and Gaussian noise, with equality/inequality rows.
* **Benchmarks** (tag `benchmark`): 500x10 up to 5000x30, 300x100 (many columns), and
  `bench_cvxpy_script_1000x250` (the sizes of `test_cvxpy.py`).

## What the tests assert

1. `test_solution_is_correct` - for every adapter: expected `KODE`, reported objective equals
   `||Ax-b||_1` recomputed from `x`, equals the HiGHS optimum, `x` feasible, residual vector consistent,
   and `x` equals the reference for unique optima. Tolerances follow `solver.precision`
   (`single` ~1e-4, `double` ~1e-8).
2. Closed-form, iteration-limit and metamorphic tests (row permutation, scaling of the whole problem).
3. `test_oracle_reproduces_golden_output` - adapters flagged `matches_legacy` (single precision only)
   must reproduce the recorded legacy output exactly (`KODE`, iterations, `x`, objective).
4. `test_differential.py` - about 2600 generated problems; every non-oracle adapter is compared with the
   frozen oracle **of its own precision** and must give *identical* `KODE`, iteration count, `X`, `RES`
   and `ERROR` (`np.array_equal`). Four generators: general (ties, degeneracy, duplicates, infeasible
   systems, tiny iteration limits, random `TOLER`), sign-restriction heavy (`KODE=1`), tiny integer
   problems at exact boundaries (`TOLER=0`), and problems with one redundant equality row (`TOLER=0`)
   plus three pinned seeds (see "Mutation testing"). `fuzz_corpus()` lists every case and is shared with
   the bounds-checked child process.
5. Interface tests (Fortran): each invalid argument is rejected with `L1_INVALID_INPUT` and leaves `q`
   untouched (every case violates exactly one condition); arrays larger than required give bit-identical
   results; no hidden array temporaries (`-fcheck=all`).
6. Checked build: everything above against the bounds-checked build.

## Mutation testing (how much the suite can be trusted)

The suite is only useful if it fails when the code is wrong, so deliberate bugs were planted in `src/`
(and in the oracle copy while developing the harness) and the suite had to fail.

* **Algorithm (13 mutations)** - all killed: tie-break `>=`->`>`, `KFORCE` test, residual / x sign
  branches, both `<= TOLER` optimality tests, bypass-exit condition, infeasible `KODE`, iteration-limit
  `>=`->`>`, `SUM` index range, artificial-pivot guard, Gauss-Jordan swap condition. Three of them
  initially *survived* (a swapped residual-sign test and `<=` vs `<` against `TOLER`, twice); each exposed
  a gap in the fuzz generator, which is why the sign-restricted and boundary (`TOLER=0`) generators exist.
* **Swap helpers (`swap_rows`, `swap_columns`, 4 mutations)** - all killed, but one needed a new
  generator: leaving the label column out of the *artificial-pivot* row exchange survived 2300 cases.
  That branch is reached by only ~0.3% of the general fuzz cases (6 of 2327), and even then the label
  column only changes the reported residuals/objective (same `x`, `KODE`, iterations): 3 of 5000 cases
  that reach it. A generator with a redundant equality row and `TOLER=0` reaches it in ~8% of cases, and
  the three seeds that expose the label column are pinned (`PINNED_REDUNDANT_SEEDS`) so that thousands of
  random cases are not needed.
* **`compute_marginal_costs` (5 mutations)** - 3 killed (first-sign cost for flipped rows, skipped
  column, shortened correction loop). Two survive, for different reasons:
  * *Accumulating in `wp` instead of `dp`* is an equivalent mutant while `wp = c_double`. It is killed
    (63 failures) when `wp` is temporarily set to `c_float`; in that configuration the refactored `src`
    also reproduces the **original single-precision CALGO 552 bit-for-bit** (2775 passed), i.e. the switch to
    double precision is the only numerical difference left from the original.
  * *Using the first-sign cost for a flipped label in the **second** (correction) loop* is reached often
    (582 times in the corpus) but never observable: in all 582 events `cu(1,j) == cu(2,j)`. Empirical, not
    a proof; the sign test there looks like defensive symmetry in the original code. Not worth a test.
* **Phase 1 set-up (`set_up_phase1_costs`, `apply_sign_restrictions`, 12 mutations)** - all killed on the
  first attempt (wrong sign rows, shifted or shortened ranges, inverted conditions, forced-phase-1
  logic, the `kode` call-site test); the narrowest kill is a negative-sign residual restriction
  whose row label is positive (8 failing tests). No new generator was needed: the sign-restricted
  and boundary generators added earlier cover these branches.
* **Phase 2 set-up (`set_up_phase2_costs`, 9 mutations run)** - 7 killed (residual cost ranges and signs,
  restriction test of either sign, rows counted but not moved, rows moved but not counted). Two survive:
  *not resetting the cost to 0* of a restricted variable, for the positive or the negative label. Not a
  generator gap: on a throw-away instrumented copy of the frozen oracle the reset never changed a
  value in the whole corpus (2627 cases) nor in 3500 extra problems that were mostly residual sign
  restrictions. Empirical, not a proof; the reset looks like defensive symmetry in the original code.
  A tenth mutation (leaving `phase = 1` after the set-up) makes the solver loop forever and cannot be run.
* **Gauss-Jordan step (`pivot_tableau`, 11 mutations)** - 10 killed on the first failing test (wrong sign of
  the multiplier or of the pivot, skipped last column or objective row, pivot row modified, pivot cell not
  inverted, labels not exchanged, wrong label read by the caller). One is an *equivalent mutant*,
  established by reading the code: dropping the guard `j /= pivot_column` on the division of the pivot row
  changes the pivot cell, which is not read before it is overwritten by `1/pivot`.
* **Entering-column selection (`select_entering_column`, 14 mutations)** - all killed (first failing test
  each; run with `-x`): both tie-break directions, either restriction flag ignored or read from the wrong
  sign, the forced-original-variables rule and its boundary, the optimality threshold, either loop end,
  the arithmetic of both gains, and the exchange of the two gains for a flipped column.
* **Sign flip of the entering column (`orient_entering_column`, 6 mutations)** - all killed: inverted
  condition, column never negated, label row or first row left out of the negation, gain not stored,
  negated gain stored.
* **Artificial-pivot search (`pivot_on_restricted_row`, 10 mutations)** - 9 killed: both tie-break and
  threshold directions, signed instead of absolute value, search one row short, pivot row not moved, row
  count not decreased, stale pivot, `found` never reported, swap one column short (the label column; it
  is also pinned by `PINNED_REDUNDANT_SEEDS`). One is an *equivalent mutant*, established by reading the
  code: running the search in phase 1 as well. `iq` is only increased by `set_up_phase2_costs` (which
  sets phase 2) and only decreased by the search itself, so it is 0 throughout phase 1 and the helper
  returns immediately; the `iphase /= 1` guard at the call site is redundant but kept (it is the
  original logic).
* **New interface code (12 of 13 verified)** - all six validation conditions (`k<1`, `n<1`, rows, columns,
  `x`, `res`), the status constant, the shim argument order and the early-return reset were killed.
* **Known gaps**
  * Undersized workspaces (`cu`, `iu`, `s` one element short) corrupt memory. They are detected only
    indirectly: the `-fcheck=all` build aborts on the degenerate fuzz cases. That abort used to kill the
    whole pytest session; since the checked build moved into a child process
    (`test_checked_build.py`) it is one test failure, but this was **not re-run** after the change.
  * The allocation-failure branch (`L1_ALLOC_FAILED`) has no test (it needs an allocation to fail).
  * Never run undersized-workspace mutations against the unchecked `src` adapter: it hangs or crashes.
    Run them only against `tests/test_checked_build.py`, one at a time, with a time limit (macOS has no
    `timeout`; use the tool's own timeout or `gtimeout`).

## Known property of the legacy code

Single precision limits it: on `bench_cvxpy_script_1000x250` the default `TOLER` (2e-5) ends with
`KODE=2` (rounding errors) after 182 iterations; `TOLER=1e-3` solves it (`iter=936`, objective within
3e-7 of HiGHS). The instance carries `toler=1e-3` for single-precision builds only. The double-precision
`src` build solves it with the default tolerance (`KODE=0`, 921 iterations, objective within 2e-15 of
HiGHS), and with any `TOLER` from 1e-12 to 1e-6.

## Adding another implementation (e.g. the CPython extension)

Write an adapter with `name`, `precision`, `is_oracle = False`, `matches_legacy` and
`solve(inst) -> Result`, and register it in `solvers.available_solvers()`. Every correctness and
differential test and the benchmark then run against it automatically
(`uv run tests/benchmark.py --solver <name> --compare tests/benchmarks/legacy_baseline.json`).
