"""
Low-level binding of the Fortran L1 solver. Use l1fit.solve_l1 instead.
"""

from typing import Annotated

import numpy
from numpy.typing import NDArray


def cl1(k: int, l: int, m: int, n: int, q: Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='F', device='cpu')], kode: int, toler: float, iter: int, x: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C', device='cpu')], res: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C', device='cpu')]) -> tuple[int, int, float]:
    """
    Solve min ||A x - b||_1 s.t. C x = d, E x <= f in place.

    `q` (Fortran-contiguous float64, shape (k+l+m+2, n+2)) holds [A b; C d; E f] and is destroyed,
    `x` (shape (n,)) and `res` (shape (k+l+m,)) are float64 vectors that receive the solution and
    the residuals (and, if kode == 1, hold sign restrictions on entry). The arrays are never
    converted or copied: a wrong dtype, layout or read-only array raises TypeError.
    Returns (status, iterations, objective); `iter` is the maximum number of iterations.
    """
