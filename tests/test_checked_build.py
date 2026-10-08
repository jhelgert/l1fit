"""Runs the correctness and differential tests against the bounds-checked build in a child process.

``SrcCheckedCL1`` is ``src/`` compiled with ``-fcheck=all -O0`` (array bounds, array temporaries,
uninitialised pointers, ...).  A violated check makes the Fortran run-time abort the *whole* process,
which would kill the pytest session, so the checked build only ever runs inside a child process and
its exit status becomes one ordinary test result here.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHILD_TIMEOUT = 300  # seconds; the child normally needs about 10 s


def test_bounds_checked_build_passes_all_instances_and_fuzz_cases():
    env = dict(os.environ, L1FIT_CHECKED_CHILD="1")
    try:
        child = subprocess.run(
            [sys.executable, str(HERE / "run_checked_build.py")],
            capture_output=True, text=True, env=env, timeout=CHILD_TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        raise AssertionError(
            f"the bounds-checked run did not finish within {CHILD_TIMEOUT} s (probably an infinite loop)"
        ) from None
    assert child.returncode == 0, (
        f"bounds-checked build failed (exit status {child.returncode}; a negative value is a signal, "
        f"i.e. the Fortran run-time aborted on a violated check)\n{child.stdout[-2000:]}\n{child.stderr[-3000:]}"
    )
    assert "checked: all cases passed" in child.stdout
