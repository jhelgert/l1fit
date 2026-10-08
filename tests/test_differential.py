"""Differential fuzz test: refactored sources vs. the frozen oracles, bit for bit.

Every non-oracle build is compared with the frozen oracle of *its own precision*:
``legacy/CALGO552.f`` (single) or ``legacy/f90_double/`` (double).  A refactoring must produce
*identical* output on arbitrary inputs: same ``KODE``, iteration count, objective, ``X`` and ``RES``.  The generator below deliberately stresses the unusual
branches of the simplex code (ties, degenerate vertices, infeasible systems, sign restrictions
(KODE=1), tiny iteration limits, loose/tight tolerances).
"""

from __future__ import annotations

import numpy as np
import pytest
from instances import Instance
from legacy_cl1 import LegacyCL1

N_CASES = 400
N_SIGN_CASES = 400  # cases dedicated to KODE=1 sign restrictions (found necessary by mutation testing)
N_REDUNDANT_CASES = 300  # redundant equality rows reach the artificial-pivot swap (found by mutation testing)
# Seeds of ``_redundant_equality_case`` for which the *label column* of the artificial-pivot row swap
# changes the reported residuals / objective (same x, iterations and KODE): only 3 of 5000 cases, found
# by comparing the frozen oracle with a mutant that skips that column.  Always run, so that damaging the
# swap is detected without needing thousands of random cases.
PINNED_REDUNDANT_SEEDS = (1511, 3003, 4455)
N_BOUNDARY_CASES = (
    1500  # tiny integer problems at exact-zero boundaries (TOLER=0), idem
)


def _random_case(seed: int) -> tuple[Instance, float]:
    rng = np.random.Generator(np.random.PCG64(1000 + seed))
    k = int(rng.integers(1, 41))
    n = int(rng.integers(1, 9))
    l = int(rng.choice([0, 0, 1, 2, 3]))
    m = int(rng.choice([0, 0, 2, 4, 8]))
    style = rng.choice(["gauss", "integer", "duplicate", "heavy"])

    def mat(rows: int) -> np.ndarray:
        if style == "integer":  # small integers -> exact ties and degenerate vertices
            return rng.integers(-3, 4, size=(rows, n)).astype(float)
        return rng.standard_normal((rows, n))

    A = mat(k)
    if style == "duplicate" and k >= 2:
        A[k // 2 :] = A[: k - k // 2]
    x0 = (
        rng.integers(-2, 3, size=n).astype(float)
        if style == "integer"
        else rng.standard_normal(n)
    )
    noise = rng.standard_t(2.0, size=k) if style == "heavy" else rng.standard_normal(k)
    if style == "integer":
        noise = rng.integers(-3, 4, size=k).astype(float)
    b = A @ x0 + noise
    C, E = mat(l), mat(m)
    plant = rng.random() < 0.7  # otherwise constraints may well be infeasible
    d = C @ x0 if plant else rng.standard_normal(l)
    f = (
        E @ x0 + (rng.integers(0, 3, size=m) if style == "integer" else rng.random(m))
        if plant
        else rng.standard_normal(m)
    )

    xsign = ressign = None
    if rng.random() < 0.35:  # KODE = 1
        if rng.random() < 0.8:
            xsign = rng.integers(-1, 2, size=n).astype(float)
        if rng.random() < 0.6:
            ressign = rng.integers(-1, 2, size=k).astype(float)
    max_iter = int(rng.integers(1, 6)) if rng.random() < 0.1 else None
    # TOLER = 0 is included on purpose: it exposes <= vs < slips at the exact-zero boundary
    toler = float(rng.choice([0.0, 1e-7, 2e-5, 1e-3, 1e-1]))
    inst = Instance(
        f"fuzz{seed}",
        "",
        A,
        b,
        C,
        d,
        E,
        f,
        xsign=xsign,
        ressign=ressign,
        max_iter=max_iter,
    )
    return inst, toler


def _sign_restricted_case(seed: int) -> tuple[Instance, float]:
    """Small problems dominated by KODE=1: random residual signs (any mix of -1/0/+1), x signs,
    optionally with constraint rows.  Targets the sign-restriction / phase-1 setup branches, which
    the general generator only reaches rarely."""
    rng = np.random.Generator(np.random.PCG64(50_000 + seed))
    k = int(rng.integers(2, 14))
    n = int(rng.integers(1, 5))
    l = int(rng.choice([0, 0, 1]))
    m = int(rng.choice([0, 0, 2, 3]))
    integer = seed % 2 == 1
    A = (
        rng.integers(-3, 4, size=(k, n)).astype(float)
        if integer
        else rng.standard_normal((k, n))
    )
    x0 = rng.standard_normal(n)
    b = A @ x0 + 2.0 * rng.standard_normal(k)
    C, E = rng.standard_normal((l, n)), rng.standard_normal((m, n))
    d, f = C @ x0, E @ x0 + rng.random(m)
    ressign = rng.integers(-1, 2, size=k).astype(float)
    xsign = rng.integers(-1, 2, size=n).astype(float) if rng.random() < 0.5 else None
    inst = Instance(f"sign{seed}", "", A, b, C, d, E, f, xsign=xsign, ressign=ressign)
    return inst, float(rng.choice([0.0, 2e-5, 1e-3]))


def _boundary_case(seed: int) -> tuple[Instance, float]:
    """Tiny integer problems with equality and inequality rows and ``TOLER = 0``.

    Integer data makes objective values *exactly* zero/equal at vertices, so comparisons against
    TOLER (``<=`` vs ``<``) and feasibility decisions at the phase-1 -> phase-2 hand-over are
    exercised on their exact boundary."""
    rng = np.random.Generator(np.random.PCG64(90_000 + seed))
    k, n = int(rng.integers(1, 8)), int(rng.integers(1, 4))
    l, m = int(rng.integers(0, 4)), int(rng.integers(0, 4))

    def ints(rows: int) -> np.ndarray:
        return rng.integers(-3, 4, size=(rows, n)).astype(float)

    def rhs(rows: int) -> np.ndarray:
        return rng.integers(-4, 5, size=rows).astype(float)

    inst = Instance(
        f"boundary{seed}", "", ints(k), rhs(k), ints(l), rhs(l), ints(m), rhs(m)
    )
    return inst, 0.0


def _redundant_equality_case(seed: int) -> tuple[Instance, float]:
    """Problems with one redundant (linearly dependent) equality row and ``TOLER = 0``.

    A redundant equality leaves an artificial variable in the basis at the end of phase 1; with a zero
    tolerance phase 2 then pivots it out through the "artificial pivot" row exchange, a branch that the
    other generators reach in only ~0.3% of their cases (found by mutation testing).  Here about 8% of
    the cases reach it."""
    rng = np.random.Generator(np.random.PCG64(95_000 + seed))
    n = int(rng.integers(2, 7))
    l = int(rng.integers(2, n + 1))  # noqa: E741
    k = int(rng.integers(n, 4 * n))
    m = int(rng.choice([0, 1, 2, 4]))
    integer = seed % 2 == 1

    def mat(rows: int) -> np.ndarray:
        if integer:
            return rng.integers(-3, 4, size=(rows, n)).astype(float)
        return rng.standard_normal((rows, n))

    A, E, C = mat(k), mat(m), mat(l - 1)
    C = np.vstack([C, C[0] + (C[1] if l > 2 else 0.0)])  # last equality = combination of earlier ones
    x0 = rng.integers(-2, 3, size=n).astype(float) if integer else rng.standard_normal(n)
    noise = rng.integers(-3, 4, size=k).astype(float) if integer else rng.standard_normal(k)
    slack = rng.integers(0, 3, size=m).astype(float) if integer else rng.random(m)
    inst = Instance(f"redundant{seed}", "", A, A @ x0 + noise, C, C @ x0, E, E @ x0 + slack)
    return inst, 0.0


def fuzz_corpus() -> list[tuple[str, int]]:
    """Every (generator name, seed) of the differential corpus, in one place (also used by the
    bounds-checked child process in ``run_checked_build.py`` so that the two cannot drift apart)."""
    return (
        [("random", s) for s in range(N_CASES)]
        + [("sign_restricted", s) for s in range(N_SIGN_CASES)]
        + [("boundary", s) for s in range(N_BOUNDARY_CASES)]
        + [("redundant_equality", s) for s in [*range(N_REDUNDANT_CASES), *PINNED_REDUNDANT_SEEDS]]
    )


GENERATORS = {
    "random": _random_case,
    "sign_restricted": _sign_restricted_case,
    "boundary": _boundary_case,
    "redundant_equality": _redundant_equality_case,
}


@pytest.fixture
def oracle(solver, all_solvers):
    """The frozen oracle with the same precision as the solver under test (None for oracles)."""
    if solver.is_oracle:
        pytest.skip("oracles are the reference")
    return all_solvers["legacy_double" if solver.precision == "double" else "legacy"]


@pytest.mark.parametrize("seed", range(N_CASES))
def test_matches_legacy_bit_for_bit(solver, oracle, seed):
    inst, toler = _random_case(seed)
    ref = oracle.solve(inst, toler=toler)
    got = solver.solve(inst, toler=toler)
    assert (got.kode, got.iterations) == (ref.kode, ref.iterations)
    assert np.array_equal(got.x, ref.x), "X differs"
    assert np.array_equal(got.res, ref.res), "RES differs"
    assert got.error == ref.error, "ERROR differs"


@pytest.mark.parametrize("seed", range(N_SIGN_CASES))
def test_sign_restricted_matches_legacy_bit_for_bit(solver, oracle, seed):
    inst, toler = _sign_restricted_case(seed)
    ref = oracle.solve(inst, toler=toler)
    got = solver.solve(inst, toler=toler)
    assert (got.kode, got.iterations) == (ref.kode, ref.iterations)
    assert np.array_equal(got.x, ref.x), "X differs"
    assert np.array_equal(got.res, ref.res), "RES differs"
    assert got.error == ref.error, "ERROR differs"


@pytest.mark.parametrize("seed", range(N_BOUNDARY_CASES))
def test_boundary_matches_legacy_bit_for_bit(solver, oracle, seed):
    inst, toler = _boundary_case(seed)
    ref = oracle.solve(inst, toler=toler)
    got = solver.solve(inst, toler=toler)
    assert (got.kode, got.iterations) == (ref.kode, ref.iterations)
    assert np.array_equal(got.x, ref.x), "X differs"
    assert np.array_equal(got.res, ref.res), "RES differs"
    assert got.error == ref.error, "ERROR differs"


@pytest.mark.parametrize("seed", [*range(N_REDUNDANT_CASES), *PINNED_REDUNDANT_SEEDS])
def test_redundant_equalities_match_legacy_bit_for_bit(solver, oracle, seed):
    inst, toler = _redundant_equality_case(seed)
    ref = oracle.solve(inst, toler=toler)
    got = solver.solve(inst, toler=toler)
    assert (got.kode, got.iterations) == (ref.kode, ref.iterations)
    assert np.array_equal(got.x, ref.x), "X differs"
    assert np.array_equal(got.res, ref.res), "RES differs"
    assert got.error == ref.error, "ERROR differs"


def test_fuzz_generator_covers_all_exit_codes(all_solvers):
    """The fuzz corpus is only useful if it reaches every KODE and both phases; guard that."""
    oracle = all_solvers["legacy"]
    kodes: dict[int, int] = {}
    kode1_cases = 0
    for seed in range(N_CASES):
        inst, toler = _random_case(seed)
        kode = oracle.solve(inst, toler=toler).kode
        kodes[kode] = kodes.get(kode, 0) + 1
        kode1_cases += inst.xsign is not None or inst.ressign is not None
    assert {0, 1, 3} <= set(kodes), kodes
    assert kode1_cases >= 40
