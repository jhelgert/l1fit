"""Builds and runs ``tests/fortran/check_cl1_interface.f90`` against the sources in ``src/``.

The Fortran interface ``cl1`` validates the shapes of its assumed-shape arrays, which the C ABI
(``l1_cl1``, explicit-shape arrays derived from the dimension arguments) can never get wrong, so
these checks are written in Fortran.  The program is compiled strictly (``-std=f2018 -Wall -Wextra``)
with run-time checking (``-fcheck=all``: bounds, array temporaries, ...); any failure of a check, any
run-time warning or any compiler warning other than the deliberate float comparisons fails the test.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

from fortran_src import BUILD_DIR, ROOT

PROGRAM = Path(__file__).resolve().parent / "fortran" / "check_cl1_interface.f90"
RUN_TIMEOUT = 60  # seconds; the program normally needs a fraction of a second
SOURCES = [ROOT / "src" / "l1_precision.f90", ROOT / "src" / "l1_calgo552.f90", PROGRAM]


def test_fortran_interface_checks():
    fc = os.environ.get("FC", "gfortran")
    out_dir = BUILD_DIR / "interface_checks"
    out_dir.mkdir(parents=True, exist_ok=True)
    exe = out_dir / "check_cl1_interface"
    build = subprocess.run(
        [
            fc,
            "-std=f2018",
            "-Wall",
            "-Wextra",
            "-fcheck=all",
            "-g",
            "-J",
            str(out_dir),
            "-I",
            str(out_dir),
            *map(str, SOURCES),
            "-o",
            str(exe),
        ],
        capture_output=True,
        text=True,
    )
    assert build.returncode == 0, build.stderr
    warnings = [
        line
        for line in build.stderr.splitlines()
        if "Warning" in line and "Wcompare-reals" not in line
    ]
    assert not warnings, "\n".join(warnings)

    try:
        run = subprocess.run([str(exe)], capture_output=True, text=True, timeout=RUN_TIMEOUT)
    except subprocess.TimeoutExpired:
        raise AssertionError(
            f"the Fortran interface check did not finish within {RUN_TIMEOUT} s (probably an infinite loop in cl1)"
        ) from None
    assert run.returncode == 0, run.stdout + run.stderr
    assert "all interface checks passed" in run.stdout
    assert "Fortran runtime warning" not in run.stderr, run.stderr
