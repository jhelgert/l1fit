"""Typed Python interface of the L1 solver.

The numerical work is done by the Fortran implementation of ACM TOMS algorithm 552 (Barrodale and
Roberts), reached through the thin binding :mod:`l1fit._core`. This module validates the arguments,
packs the problem into the layout of the solver, and unpacks the result.
"""

from __future__ import annotations

import enum
import math
import operator
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from . import _core

type FloatArray = npt.NDArray[np.float64]
"""A float64 array. The solver works in double precision only."""

type SignArray = npt.NDArray[np.integer] | npt.NDArray[np.floating]
"""An array of sign restrictions: every entry is -1, 0 or 1."""

DEFAULT_TOLERANCE: float = 1e-10
"""Default tolerance: ``10**(-d*2/3)`` for ``d`` of about 16 decimal digits (rule of the paper)."""

_INT_MAX = 2**31 - 1
_STATUS_INVALID_INPUT = 4  # the Fortran statuses that are errors, not results
_STATUS_ALLOCATION_FAILED = 5


class Status(enum.IntEnum):
    """Outcome of a solve (the exit codes of the Fortran routine)."""

    OPTIMAL = 0
    """An optimal solution was found."""
    INFEASIBLE = 1
    """The constraints have no feasible solution."""
    ROUNDING_ERRORS = 2
    """The calculation was stopped because of rounding errors."""
    MAX_ITERATIONS = 3
    """The maximum number of iterations was reached."""


@dataclass(frozen=True, slots=True, eq=False)
class L1Result:
    """Result of :func:`solve_l1`.

    Unless ``status`` is :attr:`Status.OPTIMAL`, the arrays hold the state at which the solver
    stopped.
    """

    x: FloatArray
    """The solution, shape ``(n,)``."""
    residual: FloatArray
    """The residuals ``b - A x`` of the fitted equations, shape ``(k,)``."""
    equality_residual: FloatArray
    """The residuals ``d - C x`` of the equality constraints (zero), shape ``(l,)``."""
    inequality_slack: FloatArray
    """The slacks ``f - E x`` of the inequality constraints (non-negative), shape ``(m,)``."""
    objective: float
    """The minimum sum of the absolute values of ``residual``."""
    status: Status
    """Why the solver stopped."""
    iterations: int
    """The number of simplex iterations."""

    @property
    def success(self) -> bool:
        """Whether an optimal solution was found."""
        return self.status is Status.OPTIMAL


def _check_array(value: object, name: str, ndim: int) -> FloatArray:
    """Check that ``value`` is a finite float64 array with ``ndim`` dimensions."""
    if not isinstance(value, np.ndarray):
        raise TypeError(f"{name} must be a numpy.ndarray, got {type(value).__name__}")
    if value.dtype != np.float64:
        raise TypeError(f"{name} must have dtype float64, got {value.dtype}")
    if value.ndim != ndim:
        raise ValueError(f"{name} must have {ndim} dimension(s), got {value.ndim}")
    if not np.isfinite(value).all():
        raise ValueError(f"{name} must contain only finite values")
    return value


def _optional_constraints(
    matrix: FloatArray | None,
    rhs: FloatArray | None,
    names: tuple[str, str],
    n_cols: int,
) -> tuple[FloatArray, FloatArray]:
    """Check an optional pair (constraint matrix, right-hand side); ``None`` means no rows."""
    matrix_name, rhs_name = names
    if matrix is None and rhs is None:
        return np.zeros((0, n_cols)), np.zeros(0)
    if matrix is None or rhs is None:
        raise ValueError(f"{matrix_name} and {rhs_name} must be given together")
    checked_matrix = _check_array(matrix, matrix_name, 2)
    checked_rhs = _check_array(rhs, rhs_name, 1)
    if checked_matrix.shape[1] != n_cols:
        raise ValueError(
            f"{matrix_name} must have {n_cols} columns like A, got {checked_matrix.shape[1]}"
        )
    if checked_rhs.shape[0] != checked_matrix.shape[0]:
        raise ValueError(
            f"{rhs_name} must have {checked_matrix.shape[0]} entries like the rows of "
            f"{matrix_name}, got {checked_rhs.shape[0]}"
        )
    return checked_matrix, checked_rhs


def _check_signs(value: SignArray | None, name: str, length: int) -> FloatArray | None:
    """Check an optional array of sign restrictions and convert it to float64."""
    if value is None:
        return None
    if not isinstance(value, np.ndarray) or value.dtype.kind not in "iuf":
        raise TypeError(f"{name} must be a numpy.ndarray of integers or floats")
    if value.shape != (length,):
        raise ValueError(f"{name} must have shape ({length},), got {value.shape}")
    if not np.isin(value, (-1, 0, 1)).all():
        raise ValueError(f"{name} must only contain -1, 0 or 1")
    return value.astype(np.float64)


def _check_tolerance(tol: float | None) -> float:
    if tol is None:
        return DEFAULT_TOLERANCE
    value = float(tol)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"tol must be finite and not negative, got {tol}")
    return value


def _check_max_iter(max_iter: int | None, n_rows: int) -> int:
    if max_iter is None:
        return min(10 * n_rows, _INT_MAX)  # the suggestion of the paper
    if isinstance(max_iter, bool):
        raise TypeError("max_iter must be an integer, got bool")
    try:
        limit = operator.index(max_iter)  # accepts int and numpy integers, rejects floats
    except TypeError:
        raise TypeError(f"max_iter must be an integer, got {type(max_iter).__name__}") from None
    if not 0 <= limit <= _INT_MAX:
        raise ValueError(f"max_iter must be between 0 and {_INT_MAX}, got {max_iter}")
    return limit


def solve_l1(  # noqa: PLR0913
    A: FloatArray,
    b: FloatArray,
    C: FloatArray | None = None,
    d: FloatArray | None = None,
    E: FloatArray | None = None,
    f: FloatArray | None = None,
    *,
    x_sign: SignArray | None = None,
    residual_sign: SignArray | None = None,
    tol: float | None = None,
    max_iter: int | None = None,
) -> L1Result:
    """Solve a constrained L1 fitting problem.

    Minimizes ``||A x - b||_1`` subject to the equality constraints ``C x = d`` and the inequality
    constraints ``E x <= f`` with a modification of the simplex method (ACM TOMS algorithm 552).
    The input arrays are never modified.

    Parameters
    ----------
    A : numpy.ndarray
        The matrix of the equations to fit, float64, shape ``(k, n)`` with ``k, n >= 1``.
    b : numpy.ndarray
        The right-hand side, float64, shape ``(k,)``.
    C : numpy.ndarray, optional
        The matrix of the equality constraints, float64, shape ``(l, n)``. Needs ``d``.
    d : numpy.ndarray, optional
        The right-hand side of the equality constraints, float64, shape ``(l,)``.
    E : numpy.ndarray, optional
        The matrix of the inequality constraints, float64, shape ``(m, n)``. Needs ``f``.
    f : numpy.ndarray, optional
        The right-hand side of the inequality constraints, float64, shape ``(m,)``.
    x_sign : numpy.ndarray, optional
        Sign restrictions of the unknowns, shape ``(n,)``: -1 requires ``x[j] <= 0``, 0 means none,
        1 requires ``x[j] >= 0``. Cheaper than the same restriction as rows of ``E``.
    residual_sign : numpy.ndarray, optional
        Sign restrictions of the residuals ``b - A x``, shape ``(k,)``, with the same meaning.
    tol : float, optional
        Tolerance (finite, not negative): the solver cannot distinguish zero from a quantity whose
        magnitude does not exceed it. Default :data:`DEFAULT_TOLERANCE`.
    max_iter : int, optional
        Maximum number of simplex iterations. Default ``10 * (k + l + m)``.

    Returns
    -------
    L1Result
        The solution and the status. Infeasibility, rounding errors and reaching ``max_iter`` are
        reported through ``status``, not as exceptions.

    Raises
    ------
    TypeError
        If an array is not a numpy array or does not have dtype float64 (sign arrays may also be
        integer arrays), or ``max_iter`` is not an integer.
    ValueError
        If shapes do not match, an array contains NaN or infinity, ``C`` or ``E`` is given without
        its right-hand side, or a number is out of range.
    MemoryError
        If the solver cannot allocate its workspace.
    RuntimeError
        If the solver rejects the input although it was validated here (a bug).
    """
    matrix = _check_array(A, "A", 2)
    n_rows, n_cols = matrix.shape
    if n_rows < 1 or n_cols < 1:
        raise ValueError(f"A must have at least one row and one column, got shape {matrix.shape}")
    rhs = _check_array(b, "b", 1)
    if rhs.shape[0] != n_rows:
        raise ValueError(f"b must have {n_rows} entries like the rows of A, got {rhs.shape[0]}")
    eq_matrix, eq_rhs = _optional_constraints(C, d, ("C", "d"), n_cols)
    ineq_matrix, ineq_rhs = _optional_constraints(E, f, ("E", "f"), n_cols)
    unknown_signs = _check_signs(x_sign, "x_sign", n_cols)
    residual_signs = _check_signs(residual_sign, "residual_sign", n_rows)
    tolerance = _check_tolerance(tol)

    n_eq, n_ineq = eq_matrix.shape[0], ineq_matrix.shape[0]
    n_all = n_rows + n_eq + n_ineq
    iteration_limit = _check_max_iter(max_iter, n_all)

    # The solver works on one matrix [A b; C d; E f] with two extra rows and columns (destroyed).
    q = np.zeros((n_all + 2, n_cols + 2), dtype=np.float64, order="F")
    q[:n_rows, :n_cols] = matrix
    q[:n_rows, n_cols] = rhs
    q[n_rows : n_rows + n_eq, :n_cols] = eq_matrix
    q[n_rows : n_rows + n_eq, n_cols] = eq_rhs
    q[n_rows + n_eq : n_all, :n_cols] = ineq_matrix
    q[n_rows + n_eq : n_all, n_cols] = ineq_rhs

    x = np.zeros(n_cols, dtype=np.float64)
    res = np.zeros(n_all, dtype=np.float64)
    entry_flag = 0
    if unknown_signs is not None or residual_signs is not None:
        entry_flag = 1  # the solver then reads the sign restrictions from x and res
        if unknown_signs is not None:
            x[:] = unknown_signs
        if residual_signs is not None:
            res[:n_rows] = residual_signs

    status_code, iterations, objective = _core.cl1(
        n_rows, n_eq, n_ineq, n_cols, q, entry_flag, tolerance, iteration_limit, x, res
    )
    if status_code == _STATUS_INVALID_INPUT:
        raise RuntimeError("the solver rejected arguments that passed validation (this is a bug)")
    if status_code == _STATUS_ALLOCATION_FAILED:
        raise MemoryError("the solver could not allocate its workspace")

    return L1Result(
        x=x,
        residual=res[:n_rows].copy(),
        equality_residual=res[n_rows : n_rows + n_eq].copy(),
        inequality_slack=res[n_rows + n_eq :].copy(),
        objective=objective,
        status=Status(status_code),
        iterations=iterations,
    )
