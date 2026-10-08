"""Fast constrained L1 curve fitting.

A modernized Fortran implementation of ACM TOMS algorithm 552 (Barrodale and Roberts) that solves
``min ||A x - b||_1`` subject to ``C x = d`` and ``E x <= f``, with a typed numpy interface.
"""

from importlib.metadata import version

from ._solver import DEFAULT_TOLERANCE, FloatArray, L1Result, SignArray, Status, solve_l1

__version__ = version("l1fit")

__all__ = [
    "DEFAULT_TOLERANCE",
    "FloatArray",
    "L1Result",
    "SignArray",
    "Status",
    "__version__",
    "solve_l1",
]
