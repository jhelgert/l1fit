"""ctypes adapters around the Fortran solver in ``src/`` (through the C entry point ``l1_cl1``).

* ``SrcCL1``        - the sources compiled as they are (``-O2``).
* ``SrcCheckedCL1`` - the same sources built with ``-fcheck=all -O0`` (bounds, array temporaries, ...).

Each is compiled into a shared library under ``tests/_build/`` (keyed by a hash of the source content and
the flags), so no f2py or numpy build machinery is needed. The adapters test the Fortran code on its own,
without the Python binding; ``extension_adapter.py`` tests the installed package.

Environment variable ``FC`` selects the Fortran compiler (default: gfortran).
"""

from __future__ import annotations

import ctypes
import hashlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from solvers import MAX_HANGS, Result, SolverHang, call_with_timeout, solve_timeout

ROOT = Path(__file__).resolve().parent.parent
BUILD_DIR = Path(__file__).resolve().parent / "_build"
# Compile order matters: modules before their users.
SRC_FILES = [
    ROOT / "src" / "l1_precision.f90",
    ROOT / "src" / "l1_calgo552.f90",
    ROOT / "src" / "l1_c_api.f90",
]
CHECK_FFLAGS = ["-O0", "-g", "-fcheck=all", "-ffp-contract=off", "-fPIC", "-w"]

# -ffp-contract=off: no fused multiply-add, so iteration counts are reproducible across
# optimisation levels and CPU architectures (arm64 would otherwise contract a*b+c).
FFLAGS = ["-O2", "-ffp-contract=off", "-fPIC", "-w"]

# CALGO 552 suggests TOLER = 10**(-D*2/3) for D decimal digits of accuracy: D ~ 16 in double precision.
DEFAULT_TOLER_DOUBLE = 1.0e-10


def _lib_suffix() -> str:
    if sys.platform == "darwin":
        return ".dylib"
    if sys.platform.startswith("win"):
        return ".dll"
    return ".so"


def build_library(sources: list[Path], fflags: list[str] | None = None) -> Path:
    """Compile ``sources`` (in order) into a shared library, reusing a build with identical content."""
    sources = [Path(s).resolve() for s in sources]
    BUILD_DIR.mkdir(exist_ok=True)
    fc = os.environ.get("FC", "gfortran")
    flags = list(fflags or FFLAGS)
    # The library name is keyed by the source *content* and the build flags, so a different source (or
    # an edited one) can never pick up a stale build.
    h = hashlib.sha256(" ".join([fc, *flags]).encode())
    for s in sources:
        h.update(s.read_bytes())
    digest = h.hexdigest()[:12]
    lib = BUILD_DIR / f"libcl1_{digest}{_lib_suffix()}"
    if lib.exists():
        return lib
    mod_dir = BUILD_DIR / f"mod_{digest}"
    mod_dir.mkdir(exist_ok=True)
    cmd = [fc, *flags, "-J", str(mod_dir), "-shared", *map(str, sources), "-o", str(lib)]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return lib


def load_library(path: Path) -> ctypes.CDLL:
    """Load a library built by ``build_library``.

    On Windows, Python does not search the PATH for the DLLs that a library depends on: the directory of
    the compiler (libgfortran, libgcc_s, libquadmath of a bounds-checked build) has to be added explicitly.
    """
    if sys.platform == "win32":
        compiler = shutil.which(os.environ.get("FC", "gfortran"))
        if compiler is not None:
            os.add_dll_directory(str(Path(compiler).resolve().parent))
    return ctypes.CDLL(str(path))


class SrcCL1:
    """The sources in ``src/`` through the C entry point ``l1_cl1`` (``l1_c_api.f90``, 11 arguments).

    Packs an instance into ``Q``/``X``/``RES``, calls the library and unpacks the result.
    """

    name = "src"
    precision = "double"
    fflags: list[str] | None = None

    def __init__(self, toler: float = DEFAULT_TOLER_DOUBLE) -> None:
        self.toler = toler
        self._hangs = 0  # number of calls that did not return in time
        self._lib = load_library(build_library(SRC_FILES, self.fflags))
        real_arr = np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS")
        int_ptr = ctypes.POINTER(ctypes.c_int)
        self._fn = self._lib.l1_cl1  # K L M N Q KODE TOLER ITER X RES ERROR
        self._fn.restype = None
        self._fn.argtypes = [
            *[int_ptr] * 4,
            real_arr,
            int_ptr,
            ctypes.POINTER(ctypes.c_double),
            int_ptr,
            real_arr,
            real_arr,
            ctypes.POINTER(ctypes.c_double),
        ]

    def solve(self, inst, toler: float | None = None, max_iter: int | None = None) -> Result:
        k, l, m, n = inst.k, inst.l, inst.m, inst.n
        klm = k + l + m

        q = np.zeros((klm + 2, n + 2), dtype=np.float64, order="F")
        q[:k, :n] = inst.A
        q[:k, n] = inst.b
        q[k : k + l, :n] = inst.C
        q[k : k + l, n] = inst.d
        q[k + l : klm, :n] = inst.E
        q[k + l : klm, n] = inst.f

        x = np.zeros(n, dtype=np.float64)
        res = np.zeros(klm, dtype=np.float64)
        kode = 0
        if inst.xsign is not None or inst.ressign is not None:
            kode = 1
            if inst.xsign is not None:
                x[:n] = inst.xsign
            if inst.ressign is not None:
                res[:k] = inst.ressign

        if max_iter is None:
            max_iter = inst.max_iter if inst.max_iter is not None else 10 * klm
        kode_c, iter_c = ctypes.c_int(kode), ctypes.c_int(max_iter)
        toler_c, error_c = ctypes.c_double(self.toler if toler is None else toler), ctypes.c_double(0.0)

        if self._hangs >= MAX_HANGS:
            raise SolverHang(
                f"{self.name}: {self._hangs} earlier calls did not return; failing immediately "
                "instead of waiting again"
            )
        try:
            call_with_timeout(
                lambda: self._invoke(k, l, m, n, q, kode_c, toler_c, iter_c, x, res, error_c),
                solve_timeout(),
            )
        except SolverHang:
            self._hangs += 1
            raise
        return Result(
            kode=int(kode_c.value),
            iterations=int(iter_c.value),
            error=float(error_c.value),
            x=x.copy(),
            res=res.copy(),
        )

    def _invoke(self, k, l, m, n, q, kode_c, toler_c, iter_c, x, res, error_c) -> None:
        dims = [ctypes.c_int(v) for v in (k, l, m, n)]
        self._fn(
            *[ctypes.byref(d) for d in dims],
            q,
            ctypes.byref(kode_c),
            ctypes.byref(toler_c),
            ctypes.byref(iter_c),
            x,
            res,
            ctypes.byref(error_c),
        )


class SrcCheckedCL1(SrcCL1):
    """``src/`` built with ``-fcheck=all -O0``: bounds, array temporaries, uninitialised use, ..."""

    name = "src_checked"
    fflags = CHECK_FFLAGS
