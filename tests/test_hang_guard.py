"""Tests of the protection against solver calls that never return.

A bug that makes the solver loop forever must fail a test instead of hanging the whole run. An
iteration limit cannot do that (a loop that never pivots never increments the iteration counter), so
every adapter runs its foreign call in a worker thread with a wall-clock limit, see
``solvers.call_with_timeout``.
"""

from __future__ import annotations

import ctypes
import time

import pytest

import legacy_cl1
import solvers
from instances import load_instances
from solvers import SolverHang, call_with_timeout

libc = ctypes.CDLL(None)  # a foreign call that blocks inside C (and releases the GIL), like a stuck Fortran loop


def test_returns_normally_and_propagates_exceptions():
    done = []
    call_with_timeout(lambda: done.append(1), timeout=5.0)
    assert done == [1]

    def fails():
        raise ValueError("boom")

    with pytest.raises(ValueError, match="boom"):
        call_with_timeout(fails, timeout=5.0)


def test_a_call_stuck_in_c_becomes_solver_hang_quickly():
    start = time.perf_counter()
    with pytest.raises(SolverHang, match="did not return"):
        call_with_timeout(lambda: libc.sleep(30), timeout=0.3)
    assert time.perf_counter() - start < 5.0, "the guard must not wait for the stuck call"


class StuckAdapter(legacy_cl1.LegacyCL1):
    """The frozen oracle whose foreign call never returns (it sleeps 30 s inside C)."""

    name = "stuck"

    def _invoke(self, *args) -> None:
        libc.sleep(30)


def test_adapter_reports_a_hang_then_fails_fast(monkeypatch):
    monkeypatch.setenv("L1FIT_SOLVE_TIMEOUT", "0.2")
    inst = next(i for i in load_instances() if i.name == "median_1d")
    adapter = StuckAdapter()

    for _ in range(solvers.MAX_HANGS):
        with pytest.raises(SolverHang, match="did not return"):
            adapter.solve(inst)

    start = time.perf_counter()
    with pytest.raises(SolverHang, match="failing immediately"):
        adapter.solve(inst)
    assert time.perf_counter() - start < 0.1, "after MAX_HANGS hangs further calls must fail at once"


def test_default_timeout_is_generous_and_configurable(monkeypatch):
    monkeypatch.delenv("L1FIT_SOLVE_TIMEOUT", raising=False)
    assert solvers.solve_timeout() == solvers.DEFAULT_SOLVE_TIMEOUT >= 10.0  # slowest legitimate call: ~1.5 s
    monkeypatch.setenv("L1FIT_SOLVE_TIMEOUT", "2.5")
    assert solvers.solve_timeout() == 2.5
