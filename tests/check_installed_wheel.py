"""Smoke test of an *installed* wheel (run by cibuildwheel after installing it, see pyproject.toml).

It checks what cannot be seen from the source tree:

* the package comes from site-packages, not from the checkout;
* the binary needs no compiler runtime: only system libraries (the wheels do not bundle any, because the
  Fortran solver does not use libgfortran, see CMakeLists.txt);
* a solve gives the right answer, and the type information ships with the wheel.

Standard library and numpy only, so that it runs in a bare test environment.
"""

from __future__ import annotations

import platform
import subprocess
import sys
from pathlib import Path

import numpy as np

import l1fit
from l1fit import _core

# Libraries that every system provides. Anything else (libgfortran, libquadmath, libgcc_s, libgomp...) would
# make the wheel fail on a machine without the compiler that built it.
ALLOWED_PREFIXES = {
    "Darwin": ("/usr/lib/libSystem", "/usr/lib/libc++", "/usr/lib/libc++abi"),
    "Linux": (
        "libc.so", "libm.so", "libdl.so", "libpthread.so", "librt.so", "libutil.so",
        "ld-linux", "linux-vdso", "libstdc++.so", "libgcc_s.so",
    ),
}


def needed_libraries(module: Path) -> list[str]:
    system = platform.system()
    if system == "Darwin":
        lines = subprocess.run(["otool", "-L", str(module)], capture_output=True, text=True, check=True).stdout
        return [line.split()[0] for line in lines.splitlines()[1:] if line.strip()]
    if system == "Linux":
        lines = subprocess.run(["readelf", "-d", str(module)], capture_output=True, text=True, check=True).stdout
        return [
            line.split("[")[1].rstrip("]")
            for line in lines.splitlines()
            if "(NEEDED)" in line and "[" in line
        ]
    return []


def main() -> None:
    module = Path(_core.__file__)
    assert "site-packages" in module.parts, f"not an installed wheel: {module}"

    libraries = needed_libraries(module)
    allowed = ALLOWED_PREFIXES.get(platform.system())
    if allowed is not None:
        unexpected = [lib for lib in libraries if not lib.startswith(allowed)]
        assert not unexpected, (
            f"{module.name} needs libraries that the wheel does not bundle: {unexpected} "
            "(a Fortran or compiler runtime must not be needed)"
        )

    # The median of five numbers: x = 3, objective 107.
    result = l1fit.solve_l1(np.ones((5, 1)), np.array([1.0, 2.0, 3.0, 10.0, 100.0]))
    assert result.status is l1fit.Status.OPTIMAL
    assert abs(result.x[0] - 3.0) < 1e-12 and abs(result.objective - 107.0) < 1e-9

    package = Path(l1fit.__file__).parent
    assert (package / "py.typed").exists() and (package / "_core.pyi").exists()

    print(f"ok: l1fit {l1fit.__version__} on Python {sys.version.split()[0]} ({platform.system()}),")
    print(f"    {module.name} needs only: {', '.join(libraries)}")


if __name__ == "__main__":
    main()
