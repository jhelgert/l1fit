"""Benchmark of ``l1fit.solve_l1`` against CVXPY on the persisted ``bench`` instances of ``tests/instances``.

    uv run --group bench python benchmarks/bench_cvxpy.py              # prints a Markdown table
    uv run --group bench python benchmarks/bench_cvxpy.py --json benchmarks/results/cvxpy.json

The CVXPY formulation of  min ||A x - b||_1  s.t.  C x = d,  E x <= f  is the natural one::

    x = cp.Variable(n)
    cp.Problem(cp.Minimize(cp.norm(b - A @ x, 1)), [C @ x == d, E @ x <= f])

CVXPY turns it into the standard LP with one epigraph variable t_i per row (|b_i - A_i x| <= t_i, i.e.
n + k variables and 2k + m inequality rows) and hands it to the solver:

* HiGHS dual simplex (``HIGHS`` with ``solver="simplex"``),
* HiGHS interior point (``HIGHS`` with ``solver="ipm"``),
* Clarabel, an interior-point solver for conic problems (``CLARABEL``).

Two times are reported for CVXPY: the whole ``cp.Problem(...)`` + ``solve()`` call (what a user waits for:
canonicalization, the solver, and unpacking) and the time that the solver itself reports. ``l1fit`` is timed
as the whole ``solve_l1`` call. All times are the minimum of several repetitions, single-threaded where the
solver allows it (HiGHS may use more).
"""

from __future__ import annotations

import argparse
import json
import platform
import statistics
import sys
import time
from pathlib import Path

import cvxpy as cp
import numpy as np

import l1fit

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tests"))
from instances import Instance, load_instances  # noqa: E402

CVXPY_SOLVERS: dict[str, tuple[str, dict]] = {
    "HiGHS simplex": ("HIGHS", {"highs_options": {"solver": "simplex"}}),
    "HiGHS ipm": ("HIGHS", {"highs_options": {"solver": "ipm"}}),
    "Clarabel": ("CLARABEL", {}),
}


def time_l1fit(inst: Instance, repeat: int) -> dict:
    times = []
    for _ in range(repeat):
        start = time.perf_counter()
        result = l1fit.solve_l1(inst.A, inst.b, inst.C, inst.d, inst.E, inst.f)
        times.append(time.perf_counter() - start)
    assert result.status is l1fit.Status.OPTIMAL
    return {"seconds": min(times), "objective": result.objective, "iterations": result.iterations, "x": result.x}


def build_problem(inst: Instance) -> tuple[cp.Problem, cp.Variable]:
    x = cp.Variable(inst.n)
    constraints = []
    if inst.l:
        constraints.append(inst.C @ x == inst.d)
    if inst.m:
        constraints.append(inst.E @ x <= inst.f)
    return cp.Problem(cp.Minimize(cp.norm(inst.b - inst.A @ x, 1)), constraints), x


def time_cvxpy(inst: Instance, solver: str, options: dict, repeat: int) -> dict:
    totals, internal = [], []
    for _ in range(repeat):
        start = time.perf_counter()
        problem, x = build_problem(inst)
        problem.solve(solver=solver, **{k: dict(v) if isinstance(v, dict) else v for k, v in options.items()})
        totals.append(time.perf_counter() - start)
        stats = problem.solver_stats
        internal.append(stats.solve_time if stats.solve_time is not None else float("nan"))
    return {
        "seconds": min(totals),
        "solver_seconds": min(internal),
        "objective": float(problem.value),
        "status": problem.status,
        "x": np.asarray(x.value),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repeat-l1fit", type=int, default=20)
    parser.add_argument("--repeat-cvxpy", type=int, default=3)
    parser.add_argument("--json", type=Path, help="write the results to this file")
    args = parser.parse_args()

    rows = []
    for inst in load_instances("benchmark"):
        ours = time_l1fit(inst, args.repeat_l1fit)
        row = {
            "instance": inst.name,
            "k": inst.k, "n": inst.n, "l": inst.l, "m": inst.m,
            "l1fit": {k: v for k, v in ours.items() if k != "x"},
            "cvxpy": {},
        }
        for label, (solver, options) in CVXPY_SOLVERS.items():
            r = time_cvxpy(inst, solver, options, args.repeat_cvxpy)
            r["objective_rel_diff"] = abs(r["objective"] - ours["objective"]) / max(1.0, abs(ours["objective"]))
            r["x_max_abs_diff"] = float(np.max(np.abs(r["x"] - ours["x"])))
            row["cvxpy"][label] = {k: v for k, v in r.items() if k != "x"}
        rows.append(row)
        print(f"done {inst.name}", file=sys.stderr)

    labels = list(CVXPY_SOLVERS)
    print("| instance (k x n, l eq, m ineq) | l1fit | " + " | ".join(f"{s} (total / solver)" for s in labels) + " |")
    print("|---|---|" + "---|" * len(labels))
    for r in rows:
        dims = f"{r['instance']} ({r['k']}x{r['n']}, {r['l']}, {r['m']})"
        cells = [f"{r['l1fit']['seconds'] * 1e3:.1f} ms"]
        for s in labels:
            c = r["cvxpy"][s]
            cells.append(f"{c['seconds'] * 1e3:.0f} / {c['solver_seconds'] * 1e3:.0f} ms ({c['seconds'] / r['l1fit']['seconds']:.0f}x)")
        print(f"| {dims} | " + " | ".join(cells) + " |")
    worst = max(c["objective_rel_diff"] for r in rows for c in r["cvxpy"].values())
    print(f"\nlargest relative objective difference to l1fit: {worst:.1e}", file=sys.stderr)

    if args.json:
        meta = {
            "machine": platform.platform(), "processor": platform.processor(), "python": sys.version.split()[0],
            "cvxpy": cp.__version__, "numpy": np.__version__, "l1fit": l1fit.__version__,
            "repeat_l1fit": args.repeat_l1fit, "repeat_cvxpy": args.repeat_cvxpy,
            "median_speedup_vs_best_cvxpy_total": statistics.median(
                min(c["seconds"] for c in r["cvxpy"].values()) / r["l1fit"]["seconds"] for r in rows
            ),
        }
        args.json.write_text(json.dumps({"meta": meta, "results": rows}, indent=2) + "\n")


if __name__ == "__main__":
    main()
