"""Child process of ``test_checked_build.py``: every case against the bounds-checked ``src`` build.

Runs (a) all persisted instances through ``checks.check_solution`` and (b) all differential fuzz
cases bit-for-bit against the frozen oracle of the same precision.  Exits with status 0 and prints
``checked: all cases passed`` on success; any assertion failure raises, and a violated Fortran
run-time check aborts the process with a signal.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import test_differential as td  # noqa: E402
from checks import check_solution  # noqa: E402
from instances import load_instances  # noqa: E402
from legacy_cl1 import FrozenDoubleCL1, LegacyCL1, SrcCheckedCL1  # noqa: E402


def main() -> None:
    checked = SrcCheckedCL1()
    oracle = FrozenDoubleCL1() if checked.precision == "double" else LegacyCL1()

    n_instances = 0
    for inst in load_instances():
        check_solution(inst, checked.solve(inst), checked.precision)
        n_instances += 1

    generators = [
        (td._random_case, td.N_CASES),
        (td._sign_restricted_case, td.N_SIGN_CASES),
        (td._boundary_case, td.N_BOUNDARY_CASES),
    ]
    n_fuzz = 0
    for make_case, count in generators:
        for seed in range(count):
            inst, toler = make_case(seed)
            ref, got = oracle.solve(inst, toler=toler), checked.solve(inst, toler=toler)
            assert (got.kode, got.iterations) == (ref.kode, ref.iterations), (make_case.__name__, seed)
            assert np.array_equal(got.x, ref.x) and np.array_equal(got.res, ref.res), (make_case.__name__, seed)
            assert got.error == ref.error, (make_case.__name__, seed)
            n_fuzz += 1

    print(f"checked: all cases passed ({n_instances} instances, {n_fuzz} fuzz cases)")


if __name__ == "__main__":
    main()
