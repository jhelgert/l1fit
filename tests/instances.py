"""Deterministic L1-fit problem instances for tests and benchmarks.

Problem solved (CALGO 552 / ``CL1``)::

    min_x ||A x - b||_1   s.t.  C x = d,   E x <= f

optionally with sign restrictions (``KODE = 1``): ``xsign[j]`` in {-1, 0, 1} restricts x_j to
<=0 / free / >=0, ``ressign[i]`` restricts the residual b_i - A_i x the same way.

Instances are built with ``numpy.random.Generator(PCG64(seed))`` (stable stream) and persisted to
``tests/instances/<name>.npz`` by ``generate_instances.py`` together with *reference values*
(``ref_*``, computed independently with SciPy/HiGHS) and *golden values* (``legacy_*``, produced
by the frozen legacy oracle).  Tests only read the persisted files, so they stay valid even if
NumPy's RNG or the Fortran code changes.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

INSTANCE_DIR = Path(__file__).resolve().parent / "instances"

# CALGO 552 exit codes
KODE_OPTIMAL = 0
KODE_INFEASIBLE = 1
KODE_ROUNDING = 2
KODE_MAX_ITER = 3


@dataclass
class Instance:
    name: str
    description: str
    A: np.ndarray
    b: np.ndarray
    C: np.ndarray
    d: np.ndarray
    E: np.ndarray
    f: np.ndarray
    tags: tuple[str, ...] = ()
    xsign: np.ndarray | None = None
    ressign: np.ndarray | None = None
    max_iter: int | None = None
    """Iteration limit handed to the solver (None -> 10*(k+l+m), the header suggestion)."""
    toler: float | None = None
    """TOLER handed to the solver (None -> adapter default).  The single-precision legacy build needs
    a looser value (1e-3) on large ill-conditioned problems, see README."""
    unique_x: bool = True
    """False when the optimal x is not unique (only the objective value can be compared)."""
    expected_kode: int = KODE_OPTIMAL
    # Filled in by generate_instances.py ----------------------------------------------------
    ref_objective: float = float("nan")
    """Optimal objective from SciPy/HiGHS in float64 (nan if infeasible / not solved)."""
    ref_status: str = ""
    ref_x: np.ndarray | None = None
    legacy: dict = field(default_factory=dict)
    """Golden output of the frozen legacy oracle: kode, iterations, error, x."""

    # ------------------------------------------------------------------------------------
    @property
    def k(self) -> int:
        return self.A.shape[0]

    @property
    def l(self) -> int:
        return self.C.shape[0]

    @property
    def m(self) -> int:
        return self.E.shape[0]

    @property
    def n(self) -> int:
        return self.A.shape[1]

    # ------------------------------------------------------------------------------------
    def save(self, directory: Path = INSTANCE_DIR) -> Path:
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / f"{self.name}.npz"
        payload: dict[str, object] = {
            "name": self.name,
            "description": self.description,
            "tags": ",".join(self.tags),
            "A": self.A,
            "b": self.b,
            "C": self.C,
            "d": self.d,
            "E": self.E,
            "f": self.f,
            "max_iter": -1 if self.max_iter is None else self.max_iter,
            "toler": float("nan") if self.toler is None else self.toler,
            "unique_x": self.unique_x,
            "expected_kode": self.expected_kode,
            "ref_objective": self.ref_objective,
            "ref_status": self.ref_status,
        }
        if self.xsign is not None:
            payload["xsign"] = self.xsign
        if self.ressign is not None:
            payload["ressign"] = self.ressign
        if self.ref_x is not None:
            payload["ref_x"] = self.ref_x
        for key, value in self.legacy.items():
            payload[f"legacy_{key}"] = value
        np.savez_compressed(path, **payload)
        return path

    @classmethod
    def load(cls, path: Path) -> Instance:
        with np.load(path, allow_pickle=False) as z:
            max_iter = int(z["max_iter"])
            toler = float(z["toler"])
            legacy = {
                k[len("legacy_") :]: z[k][()] if z[k].ndim == 0 else z[k]
                for k in z.files
                if k.startswith("legacy_")
            }
            return cls(
                name=str(z["name"]),
                description=str(z["description"]),
                tags=tuple(t for t in str(z["tags"]).split(",") if t),
                A=z["A"],
                b=z["b"],
                C=z["C"],
                d=z["d"],
                E=z["E"],
                f=z["f"],
                xsign=z["xsign"] if "xsign" in z.files else None,
                ressign=z["ressign"] if "ressign" in z.files else None,
                max_iter=None if max_iter < 0 else max_iter,
                toler=None if np.isnan(toler) else toler,
                unique_x=bool(z["unique_x"]),
                expected_kode=int(z["expected_kode"]),
                ref_objective=float(z["ref_objective"]),
                ref_status=str(z["ref_status"]),
                ref_x=z["ref_x"] if "ref_x" in z.files else None,
                legacy=legacy,
            )


def load_instances(tag: str | None = None) -> list[Instance]:
    """Load persisted instances (optionally only those carrying ``tag``), sorted by name."""
    out = [Instance.load(p) for p in sorted(INSTANCE_DIR.glob("*.npz"))]
    return [i for i in out if tag is None or tag in i.tags]


# =====================================================================================
# Instance definitions
# =====================================================================================
def _empty(n: int) -> tuple[np.ndarray, np.ndarray]:
    return np.zeros((0, n)), np.zeros(0)


def _heavy_tailed(rng: np.random.Generator, size: int) -> np.ndarray:
    """Student-t (df=2) noise: the regime where an L1 fit differs markedly from least squares."""
    return rng.standard_t(2.0, size=size)


def _random_problem(
    name: str,
    description: str,
    seed: int,
    k: int,
    n: int,
    l: int = 0,
    m: int = 0,
    tags: tuple[str, ...] = ("test",),
    noise: str = "t",
    **kwargs,
) -> Instance:
    """Random feasible problem: constraints are generated around a hidden point ``x0``."""
    rng = np.random.Generator(np.random.PCG64(seed))
    A = rng.standard_normal((k, n))
    x0 = rng.standard_normal(n)
    eps = _heavy_tailed(rng, k) if noise == "t" else rng.standard_normal(k)
    b = A @ x0 + eps
    C, d = (rng.standard_normal((l, n)), None) if l else _empty(n)
    if l:
        d = C @ x0
    E, f = (rng.standard_normal((m, n)), None) if m else _empty(n)
    if m:
        f = E @ x0 + rng.uniform(0.0, 1.0, size=m)  # slack >= 0 keeps x0 feasible
    return Instance(name, description, A, b, C, d, E, f, tags=tags, **kwargs)


def build_all() -> list[Instance]:
    """Create every instance (without reference/golden values)."""
    inst: list[Instance] = []
    t = ("test",)

    # ---------------------------------------------------------------- hand-checkable cases
    # L1 estimate of a location parameter is the median: x = 3, error = 2+1+0+7+97 = 107.
    A = np.ones((5, 1))
    b = np.array([1.0, 2.0, 3.0, 10.0, 100.0])
    inst.append(
        Instance(
            "median_1d",
            "Median of 5 numbers (unique optimum x=3, error=107)",
            A,
            b,
            *_empty(1),
            *_empty(1),
            tags=t,
        )
    )

    # Exact line y = 2x+1 with one gross outlier: L1 ignores it, recovering (1, 2); error = 10.
    xs = np.arange(7.0)
    A = np.column_stack([np.ones(7), xs])
    b = 1.0 + 2.0 * xs
    b[3] += 10.0
    inst.append(
        Instance(
            "line_with_outlier",
            "Line fit, one outlier of size 10 (x=[1,2], error=10)",
            A,
            b,
            *_empty(2),
            *_empty(2),
            tags=t,
        )
    )

    A = np.array([[2.0]])
    b = np.array([4.0])
    inst.append(
        Instance(
            "tiny_1x1",
            "Smallest problem k=n=1 (x=2, error=0)",
            A,
            b,
            *_empty(1),
            *_empty(1),
            tags=t,
        )
    )

    rng = np.random.Generator(np.random.PCG64(11))
    A = rng.standard_normal((4, 4))
    b = rng.standard_normal(4)
    inst.append(
        Instance(
            "square_exact",
            "Square nonsingular system: error = 0",
            A,
            b,
            *_empty(4),
            *_empty(4),
            tags=t,
        )
    )

    rng = np.random.Generator(np.random.PCG64(12))
    A = rng.standard_normal((2, 5))
    b = rng.standard_normal(2)
    inst.append(
        Instance(
            "underdetermined_2x5",
            "k<n: error = 0, x not unique",
            A,
            b,
            *_empty(5),
            *_empty(5),
            tags=t,
            unique_x=False,
        )
    )

    rng = np.random.Generator(np.random.PCG64(13))
    A = rng.standard_normal((12, 3))
    inst.append(
        Instance(
            "zero_rhs",
            "b = 0: x = 0, error = 0",
            A,
            np.zeros(12),
            *_empty(3),
            *_empty(3),
            tags=t,
        )
    )

    # Degenerate: duplicated rows and rows scaled copies -> ties in the ratio test.
    rng = np.random.Generator(np.random.PCG64(14))
    base = rng.standard_normal((4, 3))
    rhs = rng.standard_normal(4)
    A = np.vstack([base, base, 2.0 * base])
    b = np.concatenate([rhs, rhs, 2.0 * rhs + 0.5])
    inst.append(
        Instance(
            "degenerate_duplicate_rows",
            "Duplicated / scaled rows (degenerate vertices)",
            A,
            b,
            *_empty(3),
            *_empty(3),
            tags=t,
            unique_x=False,
        )
    )

    # ----------------------------------------------------------------- constraint handling
    inst.append(
        _random_problem("rand_20x3", "Heavy-tailed noise, unconstrained", 21, 20, 3)
    )
    inst.append(
        _random_problem("rand_eq_30x4_l2", "2 equality constraints", 22, 30, 4, l=2)
    )
    inst.append(
        _random_problem("rand_ineq_30x4_m5", "5 inequality constraints", 23, 30, 4, m=5)
    )
    inst.append(
        _random_problem(
            "rand_eq_ineq_60x6_l2_m8",
            "Equality + inequality constraints",
            24,
            60,
            6,
            l=2,
            m=8,
        )
    )
    inst.append(
        _random_problem(
            "rand_gauss_40x5",
            "Gaussian noise, unconstrained",
            25,
            40,
            5,
            noise="normal",
        )
    )

    # Equality constraints that pin the solution completely (l = n): x = C^-1 d, error from A.
    rng = np.random.Generator(np.random.PCG64(26))
    A = rng.standard_normal((10, 3))
    b = rng.standard_normal(10)
    C = rng.standard_normal((3, 3))
    d = rng.standard_normal(3)
    inst.append(
        Instance(
            "eq_fully_determined",
            "l = n equality rows fix x",
            A,
            b,
            C,
            d,
            *_empty(3),
            tags=t,
        )
    )

    # Infeasible: x >= 2 (-x <= -2) and x <= 1 -> KODE 1.
    A = np.array([[1.0], [1.0], [1.0]])
    b = np.array([0.0, 1.0, 2.0])
    E = np.array([[-1.0], [1.0]])
    f = np.array([-2.0, 1.0])
    inst.append(
        Instance(
            "infeasible_bounds",
            "x>=2 and x<=1 -> KODE=1",
            A,
            b,
            *_empty(1),
            E,
            f,
            tags=t,
            expected_kode=KODE_INFEASIBLE,
        )
    )

    # Inconsistent equality constraints: x1 = 1 and x1 = 2.
    A = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    b = np.array([1.0, 2.0, 3.0])
    C = np.array([[1.0, 0.0], [1.0, 0.0]])
    d = np.array([1.0, 2.0])
    inst.append(
        Instance(
            "infeasible_equalities",
            "Contradictory equalities -> KODE=1",
            A,
            b,
            C,
            d,
            *_empty(2),
            tags=t,
            expected_kode=KODE_INFEASIBLE,
        )
    )

    # Iteration limit: a normal problem with max_iter = 1 -> KODE 3.
    p = _random_problem("iteration_limit", "max_iter=1 -> KODE=3", 27, 30, 4, tags=t)
    p.max_iter = 1
    p.expected_kode = KODE_MAX_ITER
    inst.append(p)

    # ---------------------------------------------------------- KODE = 1 (sign restrictions)
    rng = np.random.Generator(np.random.PCG64(31))
    A = rng.standard_normal((25, 3))
    b = A @ np.array([-1.0, 2.0, -3.0]) + _heavy_tailed(rng, 25)
    inst.append(
        Instance(
            "kode1_x_signs",
            "x1<=0, x2 free, x3>=0 (true optimum violates x3>=0)",
            A,
            b,
            *_empty(3),
            *_empty(3),
            tags=t,
            xsign=np.array([-1.0, 0.0, 1.0]),
        )
    )

    # Residual signs are made consistent with a hidden point (otherwise the problem is infeasible).
    rng = np.random.Generator(np.random.PCG64(32))
    A = rng.standard_normal((25, 3))
    eps = _heavy_tailed(rng, 25)
    rs = np.zeros(25)
    rs[:8] = 1.0
    rs[8:12] = -1.0
    eps = np.where(rs > 0, np.abs(eps), np.where(rs < 0, -np.abs(eps), eps))
    b = A @ rng.standard_normal(3) + eps
    inst.append(
        Instance(
            "kode1_residual_signs",
            "b-Ax >=0 for 8 rows, <=0 for 4 rows",
            A,
            b,
            *_empty(3),
            *_empty(3),
            tags=t,
            ressign=rs,
            unique_x=False,
        )
    )

    p = _random_problem(
        "kode1_signs_and_constraints",
        "KODE=1 combined with eq+ineq rows",
        33,
        40,
        4,
        l=1,
        m=4,
        tags=t,
    )
    p.xsign = np.array([1.0, 0.0, 0.0, -1.0])
    inst.append(p)

    # ---------------------------------------------------------------------- mid-size checks
    inst.append(_random_problem("rand_200x10", "Mid-size unconstrained", 41, 200, 10))
    inst.append(
        _random_problem(
            "rand_300x15_l3_m20", "Mid-size, 3 eq + 20 ineq", 42, 300, 15, l=3, m=20
        )
    )

    # ---------------------------------------------------------------------- benchmarks
    b_ = ("benchmark",)
    inst.append(
        _random_problem("bench_500x10", "Benchmark: 500 x 10", 101, 500, 10, tags=b_)
    )
    inst.append(
        _random_problem("bench_2000x20", "Benchmark: 2000 x 20", 102, 2000, 20, tags=b_)
    )
    inst.append(
        _random_problem(
            "bench_2000x20_l3_m10",
            "Benchmark: 2000 x 20, 3 eq, 10 ineq",
            103,
            2000,
            20,
            l=3,
            m=10,
            tags=b_,
        )
    )
    inst.append(
        _random_problem(
            "bench_300x100",
            "Benchmark: 300 x 100 (many columns)",
            104,
            300,
            100,
            tags=b_,
        )
    )
    inst.append(
        _random_problem("bench_5000x30", "Benchmark: 5000 x 30", 105, 5000, 30, tags=b_)
    )

    # Same sizes as tests/test_cvxpy.py (1000 x 250, 10 eq, 200 ineq), but a feasible problem:
    # the Gaussian data are the same, ``b`` is used (the script passes ``f`` as ``b`` by mistake).
    # No hidden feasible point is planted: with 10 equalities and 200 random inequalities in 250
    # variables the feasible set is non-empty with overwhelming probability (verified at generation).
    rs = np.random.RandomState(0)
    A = rs.randn(1000, 250)
    b = rs.randn(1000)
    C = rs.randn(10, 250)
    d = rs.randn(10)
    E = rs.randn(200, 250)
    f = rs.randn(200)
    inst.append(
        Instance(
            "bench_cvxpy_script_1000x250",
            "Sizes/data of tests/test_cvxpy.py (1000x250, 10 eq, 200 ineq)",
            A,
            b,
            C,
            d,
            E,
            f,
            tags=b_,
            toler=1e-3,  # legacy single precision breaks down (KODE=2) with the default 2e-5
        )
    )
    return inst
