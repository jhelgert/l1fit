"""Adapter that runs the test instances through the public Python interface ``l1fit.solve_l1``.

Unlike the ``ctypes`` adapters of ``fortran_src.py`` this exercises the whole stack: argument
validation, packing, the nanobind binding (``l1fit._core``) and the Fortran solver. It is registered in
``solvers.available_solvers`` only if the package is installed (``uv sync`` builds it), so the
Fortran-only tests keep working without it.
"""

from __future__ import annotations

import numpy as np

from fortran_src import DEFAULT_TOLER_DOUBLE
from solvers import MAX_HANGS, Result, SolverHang, call_with_timeout, solve_timeout


class ExtensionCL1:
    """The CPython extension ``l1fit`` (double precision), the same interface as the other adapters."""

    name = "extension"
    precision = "double"

    def __init__(self) -> None:
        import l1fit  # noqa: PLC0415 - optional dependency

        self._solve_l1 = l1fit.solve_l1
        self.toler = DEFAULT_TOLER_DOUBLE
        self._hangs = 0

    def solve(self, inst, toler: float | None = None, max_iter: int | None = None) -> Result:
        if self._hangs >= MAX_HANGS:
            raise SolverHang(f"{self.name}: {self._hangs} earlier calls did not return")
        if max_iter is None:
            max_iter = inst.max_iter  # None selects the default 10 * (k + l + m)
        outcome = []

        def call() -> None:
            outcome.append(
                self._solve_l1(
                    inst.A, inst.b, inst.C, inst.d, inst.E, inst.f,
                    x_sign=inst.xsign, residual_sign=inst.ressign,
                    tol=self.toler if toler is None else toler, max_iter=max_iter,
                )
            )

        try:
            call_with_timeout(call, solve_timeout())
        except SolverHang:
            self._hangs += 1
            raise
        r = outcome[0]
        return Result(
            kode=int(r.status),
            iterations=r.iterations,
            error=r.objective,
            x=r.x.astype(np.float64),
            res=np.concatenate([r.residual, r.equality_residual, r.inequality_slack]),
        )
