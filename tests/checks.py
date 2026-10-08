"""Solver-independent verification of an L1 solution.

Nothing here knows how the answer was computed.  A solution is accepted when
  * the exit code is the expected one,
  * the reported objective equals ``||Ax-b||_1`` recomputed in float64 from ``x``,
  * that objective equals the independent HiGHS reference optimum,
  * ``x`` satisfies equalities, inequalities and sign restrictions,
  * the reported residual vector matches ``b-Ax``, ``d-Cx``, ``f-Ex``.
Tolerances depend on the working precision of the implementation under test.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from instances import KODE_OPTIMAL, Instance
from solvers import Result


@dataclass(frozen=True)
class Tolerances:
    objective: float  # relative, vs. the reference optimum
    feasibility: float  # relative to the row scale  1 + |row|.|x| + |rhs|
    x: float  # relative, vs. the reference x (only for unique optima)


TOLERANCES = {
    "single": Tolerances(objective=1e-4, feasibility=1e-3, x=2e-3),
    "double": Tolerances(objective=1e-8, feasibility=1e-7, x=1e-6),
}


def _row_scale(M: np.ndarray, x: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    return 1.0 + np.abs(M) @ np.abs(x) + np.abs(rhs)


def check_solution(inst: Instance, res: Result, precision: str) -> None:
    tol = TOLERANCES[precision]
    assert res.kode == inst.expected_kode, (
        f"{inst.name}: KODE={res.kode}, expected {inst.expected_kode}"
    )
    if inst.expected_kode != KODE_OPTIMAL:
        return

    x = res.x
    assert x.shape == (inst.n,) and np.all(np.isfinite(x)), (
        "x must be finite with shape (n,)"
    )
    ref = inst.ref_objective
    assert inst.ref_status == "optimal" and np.isfinite(ref)

    # --- objective ------------------------------------------------------------------
    residual = inst.b - inst.A @ x
    objective = float(np.abs(residual).sum())
    scale = max(1.0, abs(ref))
    assert abs(res.error - objective) <= tol.objective * scale, (
        f"reported error {res.error} != recomputed ||Ax-b||_1 {objective}"
    )
    assert abs(objective - ref) <= tol.objective * scale, (
        f"objective {objective} differs from reference optimum {ref}"
    )

    # --- feasibility ----------------------------------------------------------------
    if inst.l:
        viol = np.abs(inst.C @ x - inst.d) / _row_scale(inst.C, x, inst.d)
        assert viol.max() <= tol.feasibility, f"equality violation {viol.max():.2e}"
    if inst.m:
        viol = (inst.E @ x - inst.f) / _row_scale(inst.E, x, inst.f)
        assert viol.max() <= tol.feasibility, f"inequality violation {viol.max():.2e}"
    if inst.xsign is not None:
        xs = tol.feasibility * (1.0 + np.abs(x).max())
        assert np.all(x[inst.xsign > 0] >= -xs), "x >= 0 restriction violated"
        assert np.all(x[inst.xsign < 0] <= xs), "x <= 0 restriction violated"
    if inst.ressign is not None:
        rs = tol.feasibility * _row_scale(inst.A, x, inst.b)
        assert np.all(residual[inst.ressign > 0] >= -rs[inst.ressign > 0]), (
            "b-Ax >= 0 violated"
        )
        assert np.all(residual[inst.ressign < 0] <= rs[inst.ressign < 0]), (
            "b-Ax <= 0 violated"
        )

    # --- residual vector returned by the solver --------------------------------------
    k, l, m = inst.k, inst.l, inst.m
    assert res.res.shape == (k + l + m,)
    expected_res = np.concatenate([residual, inst.d - inst.C @ x, inst.f - inst.E @ x])
    row_scale = np.concatenate(
        [
            _row_scale(inst.A, x, inst.b),
            _row_scale(inst.C, x, inst.d),
            _row_scale(inst.E, x, inst.f),
        ]
    )
    err = np.abs(res.res - expected_res) / row_scale
    assert err.max() <= tol.feasibility, (
        f"RES inconsistent with b-Ax, d-Cx, f-Ex ({err.max():.2e})"
    )

    # --- x itself, when the optimum is unique -----------------------------------------
    if inst.unique_x and inst.ref_x is not None:
        dx = np.abs(x - inst.ref_x) / (1.0 + np.abs(inst.ref_x))
        assert dx.max() <= tol.x, f"x differs from reference solution ({dx.max():.2e})"
