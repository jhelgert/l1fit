"""Correctness tests of the L1 solver on the persisted instances (``tests/instances``).

They run against the Fortran sources (``src``) and the installed package (``extension``), see
``conftest.py``. Nothing here compares with another implementation: the answers are checked against an
independent HiGHS reference optimum, against feasibility and against closed-form solutions.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest
from checks import TOLERANCES, check_solution
from instances import INSTANCE_DIR, KODE_OPTIMAL, load_instances


def test_instances_are_present():
    names = {i.name for i in load_instances()}
    assert len(names) >= 25, (
        "run `uv run tests/generate_instances.py` to create tests/instances"
    )
    assert {
        "median_1d",
        "infeasible_bounds",
        "iteration_limit",
        "bench_5000x30",
    } <= names
    assert (INSTANCE_DIR / "manifest.json").exists()


def test_solution_is_correct(solver, inst):
    """Exit code, objective, feasibility and x agree with the independent HiGHS reference."""
    check_solution(inst, solver.solve(inst), solver.precision)


# ---------------------------------------------------------------------------------------
# Hand-checkable problems with closed-form answers (independent of any solver / reference).
# ---------------------------------------------------------------------------------------
def _by_name(name):
    return next(i for i in load_instances() if i.name == name)


def test_median(solver):
    r = solver.solve(_by_name("median_1d"))
    assert r.kode == KODE_OPTIMAL
    assert r.x[0] == pytest.approx(3.0, abs=1e-5)
    assert r.error == pytest.approx(107.0, rel=1e-5)


def test_line_fit_ignores_outlier(solver):
    r = solver.solve(_by_name("line_with_outlier"))
    assert r.kode == KODE_OPTIMAL
    assert r.x == pytest.approx([1.0, 2.0], abs=1e-4)
    assert r.error == pytest.approx(10.0, rel=1e-5)


def test_tiny_problem(solver):
    r = solver.solve(_by_name("tiny_1x1"))
    assert (r.kode, r.x[0], r.error) == (
        KODE_OPTIMAL,
        pytest.approx(2.0),
        pytest.approx(0.0, abs=1e-6),
    )


def test_zero_rhs_gives_zero(solver):
    r = solver.solve(_by_name("zero_rhs"))
    assert r.kode == KODE_OPTIMAL
    assert np.allclose(r.x, 0.0, atol=1e-6) and r.error == pytest.approx(0.0, abs=1e-6)


def test_iteration_limit_is_respected(solver):
    inst = _by_name("iteration_limit")
    r = solver.solve(inst)
    assert r.kode == 3
    assert r.iterations <= inst.max_iter


# ---------------------------------------------------------------------------------------
# Metamorphic properties: cheap invariants of the problem that every correct solver must satisfy.
# ---------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    "name", ["rand_20x3", "rand_eq_ineq_60x6_l2_m8", "rand_300x15_l3_m20"]
)
def test_row_permutation_invariance(solver, name):
    inst = _by_name(name)
    rng = np.random.default_rng(0)
    perm = rng.permutation(inst.k)
    shuffled = dataclasses.replace(
        inst, A=inst.A[perm], b=inst.b[perm], ref_x=None, unique_x=False
    )
    r0, r1 = solver.solve(inst), solver.solve(shuffled)
    assert r0.kode == r1.kode == KODE_OPTIMAL
    assert r1.error == pytest.approx(
        r0.error, rel=TOLERANCES[solver.precision].objective
    )


@pytest.mark.parametrize(
    "name", ["rand_20x3", "rand_eq_ineq_60x6_l2_m8", "kode1_x_signs"]
)
def test_positive_scaling_of_the_whole_problem(solver, name):
    """Scaling (A,b,C,d,E,f) by s>0 leaves x unchanged and scales the objective by s."""
    inst = _by_name(name)
    s = 8.0  # power of two: scaling is exact in floating point
    scaled = dataclasses.replace(
        inst,
        A=s * inst.A,
        b=s * inst.b,
        C=s * inst.C,
        d=s * inst.d,
        E=s * inst.E,
        f=s * inst.f,
        ref_x=None,
        ref_objective=s * inst.ref_objective,
    )
    r0, r1 = solver.solve(inst), solver.solve(scaled)
    assert r0.kode == r1.kode == KODE_OPTIMAL
    assert r1.error == pytest.approx(
        s * r0.error, rel=TOLERANCES[solver.precision].objective
    )
    tol = TOLERANCES[solver.precision].x
    assert np.allclose(r1.x, r0.x, rtol=tol, atol=tol)
