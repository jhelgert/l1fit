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

import os
import threading
from collections.abc import Callable
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


class SolverHang(RuntimeError):
    """A solver call did not return within the time limit (it is probably looping forever)."""


DEFAULT_SOLVE_TIMEOUT = 30.0
"""Seconds. The slowest legitimate call (the ``-O0 -fcheck=all`` build on the 1000x250 instance) takes
about 1.5 s. Override with the environment variable ``L1FIT_SOLVE_TIMEOUT``."""

MAX_HANGS = 3
"""After this many hangs an adapter fails every further call immediately instead of waiting again."""


def solve_timeout() -> float:
    return float(os.environ.get("L1FIT_SOLVE_TIMEOUT", DEFAULT_SOLVE_TIMEOUT))


def call_with_timeout(function: Callable[[], None], timeout: float) -> None:
    """Run ``function`` (a foreign call such as a ``ctypes`` call) and wait at most ``timeout`` seconds.

    ``ctypes`` releases the GIL during a call, so the calling thread stays responsive. A foreign call
    that does not return cannot be interrupted or killed from Python: its thread is abandoned (it is a
    daemon thread, so it cannot block interpreter exit, but it keeps a core busy until then) and
    :class:`SolverHang` is raised, so that a hang becomes an ordinary test failure.

    Note that an iteration limit cannot do this job: a loop that never pivots never increments the
    iteration counter.
    """
    outcome: dict[str, BaseException] = {}

    def target() -> None:
        try:
            function()
        except BaseException as exc:  # noqa: BLE001 - re-raised in the caller
            outcome["error"] = exc

    worker = threading.Thread(target=target, daemon=True, name="solver-call")
    worker.start()
    worker.join(timeout)
    if worker.is_alive():
        raise SolverHang(
            f"the solver did not return within {timeout:g} s (probably an infinite loop); "
            "the call is abandoned in a background thread"
        )
    if "error" in outcome:
        raise outcome["error"]


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
