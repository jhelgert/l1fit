"""Solver adapters used by the test-suite and the benchmark script.

Every adapter exposes the same tiny interface so that tests are written once
and run against the frozen legacy oracle today and against the modernized
Fortran / CPython extension later:

    adapter.name         -> str
    adapter.precision    -> "single" | "double"   (selects comparison tolerances)
    adapter.is_oracle    -> bool  (True only for the frozen legacy build)
    adapter.matches_legacy -> bool (True while exact reproduction of the legacy output is expected:
                            enables the golden and differential tests)
    adapter.solve(inst)  -> Result

To plug in a new implementation, write an adapter and register it in
``available_solvers``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Result:
    """Outcome of one L1 solve (all arrays float64, copies owned by the caller)."""

    kode: int
    """0 optimal, 1 infeasible, 2 rounding trouble, 3 iteration limit (CALGO 552 convention)."""
    iterations: int
    error: float
    """Objective value as reported by the solver (sum of |residuals|)."""
    x: np.ndarray
    """Solution vector, shape (n,)."""
    res: np.ndarray
    """Residuals b-Ax, d-Cx, f-Ex, shape (k+l+m,)."""


def available_solvers() -> dict[str, object]:
    """Return all solver adapters that can be built/imported in this environment."""
    from legacy_cl1 import FrozenDoubleCL1, LegacyCL1, SrcCL1

    solvers: dict[str, object] = {
        "legacy": LegacyCL1(),
        "legacy_double": FrozenDoubleCL1(),
        "src": SrcCL1(),
    }

    # The bounds-checked build (``SrcCheckedCL1``) is deliberately not registered here: a run-time
    # check failure aborts the process, so it runs in a child process, see test_checked_build.py.

    # Future hook: once the modernized extension exists, register it, e.g.
    #
    #   try:
    #       from l1fit_adapter import ExtensionCL1
    #       solvers["extension"] = ExtensionCL1()
    #   except ImportError:
    #       pass
    return solvers
