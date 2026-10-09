"""Child process of ``test_checked_build.py``: every persisted instance against the bounds-checked build.

``SrcCheckedCL1`` is ``src/`` compiled with ``-fcheck=all -O0``. Exits with status 0 and prints
``checked: all cases passed`` on success; a failed assertion raises, and a violated Fortran run-time
check aborts the process with a signal.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from checks import check_solution  # noqa: E402
from fortran_src import SrcCheckedCL1  # noqa: E402
from instances import load_instances  # noqa: E402


def main() -> None:
    checked = SrcCheckedCL1()
    n_instances = 0
    for inst in load_instances():
        check_solution(inst, checked.solve(inst), checked.precision)
        n_instances += 1
    print(f"checked: all cases passed ({n_instances} instances)")


if __name__ == "__main__":
    main()
