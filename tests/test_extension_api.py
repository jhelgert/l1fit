"""Tests of the Python package ``l1fit``: the typed interface and the contract of the binding.

The numerical correctness and the bit-for-bit agreement with the frozen Fortran oracles is checked by
the adapter-based tests (``extension_adapter.py``, ``test_l1_instances.py``, ``test_differential.py``).
These tests cover what only the Python interface and the binding can get wrong: argument validation,
optional arguments, copies, layouts, threads and the type stub.
"""

from __future__ import annotations

import dataclasses
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pytest

l1fit = pytest.importorskip("l1fit")
from l1fit import L1Result, Status, _core, solve_l1  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent


def median_problem() -> tuple[np.ndarray, np.ndarray]:
    """Fit one constant to five numbers: the L1 optimum is their median 3, objective 107."""
    return np.ones((5, 1)), np.array([1.0, 2.0, 3.0, 10.0, 100.0])


def line_problem() -> tuple[np.ndarray, np.ndarray]:
    """A line through (0..6, 1 + 2 t) with one outlier: the L1 fit ignores it, objective 10."""
    t = np.arange(7.0)
    b = 1.0 + 2.0 * t
    b[3] += 10.0
    return np.column_stack([np.ones(7), t]), b


# ------------------------------------------------------------------------------------ results
def test_solves_the_median_problem() -> None:
    result = solve_l1(*median_problem())
    assert result.status is Status.OPTIMAL
    assert result.success
    assert result.x == pytest.approx([3.0])
    assert result.objective == pytest.approx(107.0)
    assert result.residual == pytest.approx([-2.0, -1.0, 0.0, 7.0, 97.0])
    assert result.equality_residual.shape == (0,)
    assert result.inequality_slack.shape == (0,)
    assert result.iterations >= 1


def test_result_is_immutable_and_typed() -> None:
    result = solve_l1(*line_problem())
    assert isinstance(result, L1Result)
    assert result.x == pytest.approx([1.0, 2.0])
    assert isinstance(result.status, Status)
    assert isinstance(result.iterations, int)
    assert isinstance(result.objective, float)
    with pytest.raises(dataclasses.FrozenInstanceError):
        result.objective = 0.0  # type: ignore[misc]


def test_equality_and_inequality_constraints() -> None:
    a, b = line_problem()
    # x0 + x1 == 4 pins the line to the constraint, and x1 <= 1.5 limits the slope
    result = solve_l1(
        a, b, C=np.array([[1.0, 1.0]]), d=np.array([4.0]), E=np.array([[0.0, 1.0]]), f=np.array([1.5])
    )
    assert result.status is Status.OPTIMAL
    assert result.x[0] + result.x[1] == pytest.approx(4.0)
    assert result.x[1] <= 1.5 + 1e-9
    assert result.equality_residual == pytest.approx([0.0], abs=1e-9)
    assert result.inequality_slack[0] == pytest.approx(1.5 - result.x[1])
    assert result.inequality_slack[0] >= -1e-9


def test_no_constraints_equals_zero_row_constraints() -> None:
    a, b = line_problem()
    plain = solve_l1(a, b)
    empty = solve_l1(a, b, np.zeros((0, 2)), np.zeros(0), np.zeros((0, 2)), np.zeros(0))
    assert np.array_equal(plain.x, empty.x)
    assert (plain.objective, plain.iterations) == (empty.objective, empty.iterations)


def test_sign_restrictions() -> None:
    a, b = line_problem()
    free = solve_l1(a, b)
    restricted = solve_l1(a, b, x_sign=np.array([0, -1]))  # slope <= 0, int signs are accepted
    assert free.x[1] == pytest.approx(2.0)
    assert restricted.x[1] <= 1e-12
    assert restricted.objective > free.objective
    from_floats = solve_l1(a, b, x_sign=np.array([0.0, -1.0]))
    assert np.array_equal(restricted.x, from_floats.x)
    below = solve_l1(a, b, residual_sign=np.ones(7))  # fit from below: b - A x >= 0
    assert below.status is Status.OPTIMAL
    assert (below.residual >= -1e-9).all()


def test_infeasible_and_iteration_limit_are_statuses_not_exceptions() -> None:
    a, b = line_problem()
    infeasible = solve_l1(a, b, E=np.array([[1.0, 0.0], [-1.0, 0.0]]), f=np.array([1.0, -2.0]))
    assert infeasible.status is Status.INFEASIBLE
    assert not infeasible.success
    limited = solve_l1(*line_problem(), max_iter=0)
    assert limited.status is Status.MAX_ITERATIONS
    assert limited.iterations == 0


def test_tolerance_zero_is_allowed_and_default_is_documented() -> None:
    assert solve_l1(*median_problem(), tol=0.0).status is Status.OPTIMAL
    assert l1fit.DEFAULT_TOLERANCE == 1e-10


def test_default_options_are_the_documented_ones(monkeypatch: pytest.MonkeyPatch) -> None:
    """``tol`` defaults to DEFAULT_TOLERANCE and ``max_iter`` to 10 * (k + l + m), constraints included."""
    calls = []
    real = _core.cl1

    def recording(*args):  # noqa: ANN002, ANN202
        calls.append(args)
        return real(*args)

    monkeypatch.setattr(_core, "cl1", recording)
    a, b = line_problem()  # k = 7
    solve_l1(a, b, C=np.array([[1.0, 1.0]]), d=np.array([4.0]), E=np.array([[0.0, 1.0]]), f=np.array([9.0]))
    solve_l1(a, b, tol=0.5, max_iter=3)
    (_k, _l, _m, _n, _q, _kode, tol_1, iter_1, _x, _res), second = calls
    assert (tol_1, iter_1) == (l1fit.DEFAULT_TOLERANCE, 10 * (7 + 1 + 1))
    assert (second[6], second[7]) == (0.5, 3)


@pytest.mark.parametrize(("code", "error"), [(4, RuntimeError), (5, MemoryError)])
def test_solver_errors_are_raised_not_returned(monkeypatch: pytest.MonkeyPatch, code: int, error: type) -> None:
    """Status 4 (invalid input) and 5 (allocation) are errors, not results: they need a mocked solver,
    since validated input never triggers them."""
    monkeypatch.setattr(_core, "cl1", lambda *args: (code, 0, 0.0))
    with pytest.raises(error):
        solve_l1(*median_problem())


# ------------------------------------------------------------------------------ input handling
def test_inputs_are_never_modified() -> None:
    a, b = line_problem()
    c, d = np.array([[1.0, 1.0]]), np.array([4.0])
    before = [array.copy() for array in (a, b, c, d)]
    solve_l1(a, b, c, d)
    for original, array in zip(before, (a, b, c, d), strict=True):
        assert np.array_equal(original, array)


def test_memory_layout_and_views_do_not_matter() -> None:
    a, b = line_problem()
    expected = solve_l1(a, b)
    c_order = np.ascontiguousarray(a)
    f_order = np.asfortranarray(a)
    strided = np.repeat(a, 2, axis=0)[::2]  # a non-contiguous view of the same values
    for variant in (c_order, f_order, strided):
        result = solve_l1(variant, b)
        assert np.array_equal(result.x, expected.x)
    read_only = a.copy()
    read_only.setflags(write=False)
    assert np.array_equal(solve_l1(read_only, b).x, expected.x)


@pytest.mark.parametrize(
    ("make_arguments", "error", "message"),
    [
        (lambda a, b: ([[1.0]], b), TypeError, "A must be a numpy.ndarray"),
        (lambda a, b: (a.astype(np.float32), b), TypeError, "A must have dtype float64"),
        (lambda a, b: (a, b.astype(np.int64)), TypeError, "b must have dtype float64"),
        (lambda a, b: (a[:, 0], b), ValueError, "A must have 2 dimension"),
        (lambda a, b: (a, b[:-1]), ValueError, "b must have 7 entries"),
        (lambda a, b: (a[:0], b[:0]), ValueError, "at least one row"),
        (lambda a, b: (a[:, :0], b), ValueError, "at least one row and one column"),
        (lambda a, b: (np.where(a == 1.0, np.nan, a), b), ValueError, "finite"),
        (lambda a, b: (a, np.where(b > 5.0, np.inf, b)), ValueError, "finite"),
    ],
)
def test_invalid_problem_data_is_rejected(make_arguments, error, message) -> None:  # noqa: ANN001
    a, b = line_problem()
    arguments = make_arguments(a, b)
    with pytest.raises(error, match=message):
        solve_l1(*arguments)


@pytest.mark.parametrize(
    ("keywords", "error", "message"),
    [
        ({"C": np.ones((1, 2))}, ValueError, "C and d must be given together"),
        ({"d": np.ones(1)}, ValueError, "C and d must be given together"),
        ({"E": np.ones((1, 2))}, ValueError, "E and f must be given together"),
        ({"C": np.ones((1, 3)), "d": np.ones(1)}, ValueError, "C must have 2 columns"),
        ({"C": np.ones((2, 2)), "d": np.ones(1)}, ValueError, "d must have 2 entries"),
        ({"E": np.ones((1, 2)), "f": np.ones((1, 1))}, ValueError, "f must have 1 dimension"),
        ({"E": np.ones((1, 2)), "f": np.array([np.nan])}, ValueError, "finite"),
        ({"x_sign": np.array([2, 0])}, ValueError, "only contain -1, 0 or 1"),
        ({"x_sign": np.zeros(3)}, ValueError, r"shape \(2,\)"),
        ({"x_sign": [0, 0]}, TypeError, "numpy.ndarray"),
        ({"x_sign": np.array(["a", "b"])}, TypeError, "integers or floats"),
        ({"residual_sign": np.zeros(2)}, ValueError, r"shape \(7,\)"),
        ({"tol": -1.0}, ValueError, "tol"),
        ({"tol": float("nan")}, ValueError, "tol"),
        ({"max_iter": -1}, ValueError, "max_iter"),
        ({"max_iter": 2**31}, ValueError, "max_iter"),
        ({"max_iter": 1.5}, TypeError, "integer"),
        ({"max_iter": True}, TypeError, "bool"),
    ],
)
def test_invalid_options_are_rejected(keywords, error, message) -> None:  # noqa: ANN001
    a, b = line_problem()
    with pytest.raises(error, match=message):
        solve_l1(a, b, **keywords)


def test_numpy_integer_max_iter_is_accepted() -> None:
    assert solve_l1(*median_problem(), max_iter=np.int64(50)).status is Status.OPTIMAL


# ------------------------------------------------------------------------------------ threads
def test_threads_can_solve_independent_problems_at_the_same_time() -> None:
    """The binding releases the GIL and the Fortran solver is reentrant: results must not mix."""
    rng = np.random.default_rng(7)
    problems = [(rng.standard_normal((300, 12)), rng.standard_normal(300)) for _ in range(24)]
    sequential = [solve_l1(a, b) for a, b in problems]
    with ThreadPoolExecutor(max_workers=8) as pool:
        parallel = list(pool.map(lambda problem: solve_l1(*problem), problems))
    for one, other in zip(sequential, parallel, strict=True):
        assert np.array_equal(one.x, other.x)
        assert (one.objective, one.iterations) == (other.objective, other.iterations)


@pytest.mark.skipif((os.cpu_count() or 1) < 4, reason="needs at least 4 CPUs to show a speed-up")
def test_the_gil_is_released_so_threads_run_in_parallel() -> None:
    """Parallel results alone would also hold with the GIL held; the benefit is the speed-up.

    Measured on a 10-core machine: 3.5x with 4 threads, 1.03x when the GIL is held. The threshold
    leaves a wide margin on both sides, and the best of several repetitions is compared.
    """
    rng = np.random.default_rng(3)
    a, b = rng.standard_normal((3000, 30)), rng.standard_normal(3000)  # about 8 ms per solve
    solve_l1(a, b)  # warm-up
    jobs = 8

    def timed(run) -> float:  # noqa: ANN001
        best = float("inf")
        for _ in range(3):
            start = time.perf_counter()
            run()
            best = min(best, time.perf_counter() - start)
        return best

    sequential = timed(lambda: [solve_l1(a, b) for _ in range(jobs)])

    def threaded() -> None:
        with ThreadPoolExecutor(max_workers=4) as pool:
            list(pool.map(lambda _: solve_l1(a, b), range(jobs)))

    assert sequential / timed(threaded) > 1.5


# ------------------------------------------------------------- the contract of the binding
class TestBindingContract:
    """The binding converts and copies nothing: a bad array must raise, never be fixed silently.

    The solver works in place, so a silent copy would lose the result without any error.
    """

    @staticmethod
    def arrays(k: int = 2, l: int = 0, m: int = 0, n: int = 1) -> tuple:  # noqa: E741
        klm = k + l + m
        return (
            np.zeros((klm + 2, n + 2), order="F"),
            np.zeros(n),
            np.zeros(klm),
        )

    def call(self, q, x, res, k=2, l=0, m=0, n=1) -> tuple:  # noqa: ANN001, E741
        return _core.cl1(k, l, m, n, q, 0, 1e-10, 10, x, res)

    def test_accepts_well_formed_arrays(self) -> None:
        q, x, res = self.arrays()
        status, _iterations, _objective = self.call(q, x, res)
        assert status in {0, 1, 2, 3}

    def test_results_are_written_into_the_callers_arrays(self) -> None:
        q, x, res = np.zeros((7, 3), order="F"), np.zeros(1), np.zeros(5)
        q[:5, 0] = 1.0
        q[:5, 1] = [1.0, 2.0, 3.0, 10.0, 100.0]
        status, _, objective = _core.cl1(5, 0, 0, 1, q, 0, 1e-10, 50, x, res)
        assert (status, x[0], objective) == (0, 3.0, 107.0)
        assert res.tolist() == [-2.0, -1.0, 0.0, 7.0, 97.0]

    @pytest.mark.parametrize("dtype", [np.float32, np.int64, np.complex128])
    def test_rejects_other_dtypes_instead_of_converting(self, dtype) -> None:  # noqa: ANN001
        q, x, res = self.arrays()
        for bad in ("q", "x", "res"):
            arrays = {"q": q, "x": x, "res": res}
            arrays[bad] = arrays[bad].astype(dtype, order="F" if bad == "q" else "C")
            with pytest.raises(TypeError):
                self.call(arrays["q"], arrays["x"], arrays["res"])

    def test_rejects_a_c_ordered_matrix_instead_of_copying_it(self) -> None:
        q, x, res = self.arrays()
        with pytest.raises(TypeError):
            self.call(np.ascontiguousarray(q), x, res)

    def test_rejects_non_contiguous_vectors(self) -> None:
        q, x, res = self.arrays(k=2, n=2)
        with pytest.raises(TypeError):
            self.call(q, np.zeros(4)[::2], res, n=2)

    def test_rejects_read_only_arrays(self) -> None:
        q, x, res = self.arrays()
        x.setflags(write=False)
        with pytest.raises(TypeError):
            self.call(q, x, res)

    def test_rejects_lists(self) -> None:
        q, x, res = self.arrays()
        with pytest.raises(TypeError):
            self.call(q, x.tolist(), res)

    @pytest.mark.parametrize(
        "bad",
        [
            lambda q, x, res: (q[:-1].copy(order="F"), x, res),  # one row short
            lambda q, x, res: (q[:, :-1].copy(order="F"), x, res),  # one column short
            lambda q, x, res: (np.zeros((q.shape[0] + 1, q.shape[1]), order="F"), x, res),  # one row long
            lambda q, x, res: (q, np.zeros(2), res),  # x too long
            lambda q, x, res: (q, np.zeros(0), res),  # x too short
            lambda q, x, res: (q, x, np.zeros(1)),  # res too short
            lambda q, x, res: (q, x, np.zeros(3)),  # res too long
        ],
    )
    def test_rejects_inconsistent_shapes(self, bad) -> None:  # noqa: ANN001
        """The Fortran entry point derives the shapes from k, l, m, n: only the binding can check."""
        q, x, res = bad(*self.arrays())
        with pytest.raises(ValueError, match="shape"):
            self.call(q, x, res)

    @pytest.mark.parametrize("dimensions", [(-1, 0, 0, 1), (2, -1, 0, 1), (2, 0, -1, 1), (2, 0, 0, -1)])
    def test_rejects_negative_dimensions(self, dimensions) -> None:  # noqa: ANN001
        q, x, res = self.arrays()
        with pytest.raises(ValueError, match="negative"):
            _core.cl1(*dimensions, q, 0, 1e-10, 10, x, res)

    def test_dimensions_that_would_overflow_are_rejected(self) -> None:
        q, x, res = self.arrays()
        with pytest.raises(ValueError, match="too large"):
            _core.cl1(2**31 - 1, 2**31 - 1, 0, 1, q, 0, 1e-10, 10, x, res)


# -------------------------------------------------------------------------- the type information
def test_the_stub_matches_the_compiled_module(tmp_path: Path) -> None:
    """``_core.pyi`` is generated (``python -m nanobind.stubgen``); it must not go stale."""
    pytest.importorskip("nanobind")
    generated = tmp_path / "_core.pyi"
    subprocess.run(
        [sys.executable, "-m", "nanobind.stubgen", "-m", "l1fit._core", "-o", str(generated)],
        check=True, capture_output=True,
    )
    committed = ROOT / "python" / "l1fit" / "_core.pyi"
    assert committed.read_text() == generated.read_text(), (
        "regenerate with: uv run python -m nanobind.stubgen -m l1fit._core -o python/l1fit/_core.pyi"
    )


def test_the_readme_example_runs() -> None:
    """The usage example of the README is executed, so that it cannot silently go stale."""
    import re  # noqa: PLC0415

    text = (ROOT / "README.md").read_text()
    match = re.search(r"## Usage\n\n```python\n(.*?)```", text, re.DOTALL)
    assert match, "the README has no usage example"
    namespace: dict[str, object] = {}
    exec(compile(match.group(1), "README.md", "exec"), namespace)  # noqa: S102 - our own README
    result = namespace["result"]
    assert isinstance(result, L1Result)
    assert result.status is Status.OPTIMAL


def test_package_is_marked_as_typed() -> None:
    assert (Path(l1fit.__file__).parent / "py.typed").exists()
    assert l1fit.__version__
