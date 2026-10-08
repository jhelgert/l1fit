"""ctypes adapters around the Fortran L1 solver.

Adapters (all share the packing logic of ``FortranCL1`` and the ``Result`` type):

* ``LegacyCL1``       - frozen oracle ``legacy/CALGO552.f``: the repaired original (fixed-form, single
  precision); 18-argument C ABI ``cl1_``.
* ``FrozenDoubleCL1`` - frozen oracle ``legacy/f90_double/``: first free-form, double-precision
  version (same algorithm, 18-argument ABI ``cl1_``). Reference for bit-for-bit checks of refactors.
* ``SrcCL1``          - the sources in ``src/``; slim 11-argument C ABI ``l1_cl1`` (``l1_c_api.f90``).
* ``SrcCheckedCL1``   - same sources built with ``-fcheck=all -O0`` (bounds, array temporaries, ...).

Each is compiled into a shared library under ``tests/_build/`` (keyed by a hash of the source
content and flags). No numpy/f2py build machinery is needed.

Environment variables:
    L1FIT_LEGACY_SRC   override the oracle source (used for mutation testing the test-suite)
    FC                 Fortran compiler (default: gfortran)
"""

from __future__ import annotations

import ctypes
import hashlib
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
from solvers import Result

ROOT = Path(__file__).resolve().parent.parent
BUILD_DIR = Path(__file__).resolve().parent / "_build"
DEFAULT_SRC = ROOT / "legacy" / "CALGO552.f"
# Compile order matters: modules before their users.
SRC_FILES = [
    ROOT / "src" / "l1_precision.f90",
    ROOT / "src" / "l1_calgo552.f90",
    ROOT / "src" / "l1_c_api.f90",
]
DOUBLE_ORACLE_FILES = [
    ROOT / "legacy" / "f90_double" / "l1_precision.f90",
    ROOT / "legacy" / "f90_double" / "CALGO552.f90",
]
CHECK_FFLAGS = ["-O0", "-g", "-fcheck=all", "-ffp-contract=off", "-fPIC", "-w"]

# -ffp-contract=off: no fused multiply-add, so iteration counts are reproducible across
# optimisation levels and CPU architectures (arm64 would otherwise contract a*b+c).
FFLAGS = ["-O2", "-ffp-contract=off", "-fPIC", "-w"]


def _is_free_form(src: Path) -> bool:
    return src.suffix.lower() in (".f90", ".f95", ".f03", ".f08")


def _source_flags(sources: list[Path]) -> list[str]:
    """Fixed-form legacy sources need ``-ffixed-form -std=legacy``; free-form ones need nothing."""
    if all(_is_free_form(s) for s in sources):
        return []
    assert not any(_is_free_form(s) for s in sources), (
        "do not mix fixed and free form sources"
    )
    return ["-ffixed-form", "-std=legacy"]


# Documented suggestion in the CALGO 552 header: TOLER = 10**(-D*2/3), D = decimal digits available.
# single precision -> D ~ 7 -> 2e-5;  double precision -> D ~ 16 -> 1e-10.
DEFAULT_TOLER = 2.0e-5
DEFAULT_TOLER_DOUBLE = 1.0e-10


def _wp_dtype_of_src() -> type:
    """Working precision (``wp``) declared in ``src/l1_precision.f90``: float32 or float64."""
    text = (ROOT / "src" / "l1_precision.f90").read_text()
    m = re.search(r"integer\s*,\s*parameter\s*::\s*wp\s*=\s*(\w+)", text, re.IGNORECASE)
    kinds = {
        "c_float": np.float32,
        "real32": np.float32,
        "c_double": np.float64,
        "real64": np.float64,
    }
    if not m or m.group(1).lower() not in kinds:
        raise RuntimeError(
            f"cannot determine wp from l1_precision.f90 ({m and m.group(1)})"
        )
    return kinds[m.group(1).lower()]


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
    flags = [*(fflags or FFLAGS), *_source_flags(sources)]
    # Library name is keyed by the source *content* and build flags, so a different source (or
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
    cmd = [
        fc,
        *flags,
        "-J",
        str(mod_dir),
        "-shared",
        *map(str, sources),
        "-o",
        str(lib),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return lib


class FortranCL1:
    """Base adapter: packs an instance into ``Q``/``X``/``RES``, calls the library, unpacks the result.

    Subclasses choose the sources, the dtype and the C ABI (``_bind`` / ``_invoke``).
    """

    name = "fortran"
    is_oracle = False

    def __init__(
        self,
        sources: list[Path],
        dtype: type = np.float32,
        toler: float | None = None,
        fflags: list[str] | None = None,
    ) -> None:
        self.dtype = np.dtype(dtype)
        single = self.dtype == np.float32
        self.precision = "single" if single else "double"
        # Exact reproduction of the legacy (single precision) output is only expected from
        # single precision builds; a double precision build is compared with the frozen *double*
        # oracle instead (see test_differential.py).
        self.matches_legacy = single
        self.toler = (
            toler
            if toler is not None
            else (DEFAULT_TOLER if single else DEFAULT_TOLER_DOUBLE)
        )
        self._c_real = ctypes.c_float if single else ctypes.c_double
        self._lib = ctypes.CDLL(str(build_library(sources, fflags)))
        self._real_arr = np.ctypeslib.ndpointer(dtype=self.dtype, flags="F_CONTIGUOUS")
        self._int_arr = np.ctypeslib.ndpointer(dtype=np.int32, flags="F_CONTIGUOUS")
        self._int_ptr = ctypes.POINTER(ctypes.c_int)
        self._real_ptr = ctypes.POINTER(self._c_real)
        self._bind()

    # -- C ABI hooks (overridden) -----------------------------------------------------------
    def _bind(self) -> None:
        """Declare ``cl1_``: K L M N KLMD KLM2D NKLMD N2D Q KODE TOLER ITER X RES ERROR CU IU S."""
        ip, ra, ia = self._int_ptr, self._real_arr, self._int_arr
        self._fn = self._lib.cl1_
        self._fn.restype = None
        self._fn.argtypes = [
            *[ip] * 8,
            ra,
            ip,
            self._real_ptr,
            ip,
            ra,
            ra,
            self._real_ptr,
            ra,
            ia,
            ia,
        ]

    def _invoke(self, k, l, m, n, q, kode_c, toler_c, iter_c, x, res, error_c) -> None:
        klm = k + l + m
        c = ctypes.c_int
        dims = [c(v) for v in (k, l, m, n, klm, klm + 2, n + klm, n + 2)]
        cu = np.zeros((2, n + klm), dtype=self.dtype, order="F")
        iu = np.zeros((2, n + klm), dtype=np.int32, order="F")
        s = np.zeros(klm, dtype=np.int32)
        self._fn(
            *[ctypes.byref(d) for d in dims],
            q,
            ctypes.byref(kode_c),
            ctypes.byref(toler_c),
            ctypes.byref(iter_c),
            x,
            res,
            ctypes.byref(error_c),
            cu,
            iu,
            s,
        )

    # -- common packing logic ----------------------------------------------------------------
    def solve(
        self, inst, toler: float | None = None, max_iter: int | None = None
    ) -> Result:
        k, l, m, n = inst.k, inst.l, inst.m, inst.n
        klm = k + l + m

        q = np.zeros((klm + 2, n + 2), dtype=self.dtype, order="F")
        q[:k, :n] = inst.A
        q[:k, n] = inst.b
        q[k : k + l, :n] = inst.C
        q[k : k + l, n] = inst.d
        q[k + l : klm, :n] = inst.E
        q[k + l : klm, n] = inst.f

        x = np.zeros(n + 2, dtype=self.dtype)  # legacy documentation asks for N+2 elements
        res = np.zeros(klm, dtype=self.dtype)
        kode = 0
        if inst.xsign is not None or inst.ressign is not None:
            kode = 1
            if inst.xsign is not None:
                x[:n] = inst.xsign
            if inst.ressign is not None:
                res[:k] = inst.ressign

        if max_iter is None:
            max_iter = (
                inst.max_iter if inst.max_iter is not None else 10 * klm
            )  # header suggestion
        if toler is None:
            # inst.toler is a workaround for the limits of single precision: only used there.
            use_inst = (
                self.precision == "single" and getattr(inst, "toler", None) is not None
            )
            toler = inst.toler if use_inst else self.toler
        kode_c, iter_c = ctypes.c_int(kode), ctypes.c_int(max_iter)
        toler_c, error_c = self._c_real(toler), self._c_real(0.0)

        self._invoke(k, l, m, n, q, kode_c, toler_c, iter_c, x, res, error_c)
        return Result(
            kode=int(kode_c.value),
            iterations=int(iter_c.value),
            error=float(error_c.value),
            x=x[:n].astype(np.float64),
            res=res[:klm].astype(np.float64),
        )


class LegacyCL1(FortranCL1):
    """The frozen oracle: original single-precision CALGO 552 code, never refactored."""

    name = "legacy"
    is_oracle = True

    def __init__(self) -> None:
        super().__init__(
            [Path(os.environ.get("L1FIT_LEGACY_SRC", DEFAULT_SRC))], np.float32
        )


class FrozenDoubleCL1(FortranCL1):
    """Frozen oracle #2: the first free-form double-precision version (18-argument interface)."""

    name = "legacy_double"
    is_oracle = True

    def __init__(self) -> None:
        super().__init__(DOUBLE_ORACLE_FILES, np.float64)


class SrcCL1(FortranCL1):
    """The sources in ``src/`` through the slim C ABI ``l1_cl1`` (11 arguments)."""

    name = "src"
    fflags: list[str] | None = None

    def __init__(self) -> None:
        super().__init__(SRC_FILES, _wp_dtype_of_src(), fflags=self.fflags)

    def _bind(self) -> None:
        """Declare ``l1_cl1``: K L M N Q KODE TOLER ITER X RES ERROR."""
        ip, ra = self._int_ptr, self._real_arr
        self._fn = self._lib.l1_cl1
        self._fn.restype = None
        self._fn.argtypes = [
            *[ip] * 4,
            ra,
            ip,
            self._real_ptr,
            ip,
            ra,
            ra,
            self._real_ptr,
        ]

    def _invoke(self, k, l, m, n, q, kode_c, toler_c, iter_c, x, res, error_c) -> None:
        c = ctypes.c_int
        dims = [c(v) for v in (k, l, m, n)]
        self._fn(
            *[ctypes.byref(d) for d in dims],
            q,
            ctypes.byref(kode_c),
            ctypes.byref(toler_c),
            ctypes.byref(iter_c),
            x[:n],  # the slim ABI takes exactly n elements
            res,
            ctypes.byref(error_c),
        )


class SrcCheckedCL1(SrcCL1):
    """``src/`` built with ``-fcheck=all -O0``: bounds, array temporaries, uninitialised use, ..."""

    name = "src_checked"
    fflags = CHECK_FFLAGS
