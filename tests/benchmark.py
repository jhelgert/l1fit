# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy", "scipy"]
# ///
"""Benchmark the available L1 solvers on the persisted ``benchmark`` instances.

    uv run tests/benchmark.py                       # table on stdout
    uv run tests/benchmark.py --repeat 20 --json tests/benchmarks/legacy_baseline.json
    uv run tests/benchmark.py --compare tests/benchmarks/legacy_baseline.json

Timings are wall-clock for the *whole* adapter call (data packing + Fortran), best of ``--repeat``
runs after one warm-up; they are machine specific, so compare only runs from the same machine.
Every run is also checked for correctness (objective vs. HiGHS reference), so a "fast" but wrong
refactor is reported as such.
"""

from __future__ import annotations

import argparse
import json
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from checks import check_solution
from instances import load_instances
from solvers import available_solvers


def _compiler_version() -> str:
    try:
        out = subprocess.run(
            ["gfortran", "--version"], capture_output=True, text=True
        ).stdout
        return out.splitlines()[0]
    except OSError:
        return "unknown"


def run(solver_name: str, repeat: int, tag: str) -> dict:
    solver = available_solvers()[solver_name]
    rows = []
    for inst in load_instances(tag):
        solver.solve(inst)  # warm-up
        times = []
        for _ in range(repeat):
            t0 = time.perf_counter()
            res = solver.solve(inst)
            times.append(time.perf_counter() - t0)
        try:
            check_solution(inst, res, solver.precision)
            correct = True
        except AssertionError as exc:
            correct = False
            print(f"  !! {inst.name}: {exc}", file=sys.stderr)
        rows.append(
            {
                "name": inst.name,
                "k": inst.k,
                "l": inst.l,
                "m": inst.m,
                "n": inst.n,
                "kode": res.kode,
                "iterations": res.iterations,
                "objective": res.error,
                "ref_objective": inst.ref_objective,
                "correct": correct,
                "best_ms": 1e3 * min(times),
                "median_ms": 1e3 * statistics.median(times),
            }
        )
    return {
        "solver": solver_name,
        "precision": solver.precision,
        "repeat": repeat,
        "machine": platform.platform(),
        "python": platform.python_version(),
        "fortran": _compiler_version(),
        "results": rows,
    }


def print_table(report: dict, baseline: dict | None = None) -> None:
    base = {r["name"]: r for r in baseline["results"]} if baseline else {}
    print(
        f"solver={report['solver']} ({report['precision']})  repeat={report['repeat']}  {report['fortran']}"
    )
    head = f"{'instance':32s} {'k x n':>10s} {'l':>3s} {'m':>4s} {'kode':>4s} {'iter':>6s} {'best ms':>9s} {'median ms':>10s} ok"
    print(head + ("   vs baseline" if base else ""))
    for r in report["results"]:
        line = (
            f"{r['name']:32s} {str(r['k']) + 'x' + str(r['n']):>10s} {r['l']:3d} {r['m']:4d} "
            f"{r['kode']:4d} {r['iterations']:6d} {r['best_ms']:9.3f} {r['median_ms']:10.3f} "
            f"{'yes' if r['correct'] else 'NO'}"
        )
        if r["name"] in base:
            line += f"   x{r['best_ms'] / base[r['name']]['best_ms']:.2f} time, {r['iterations'] - base[r['name']]['iterations']:+d} iter"
        print(line)
    total = sum(r["best_ms"] for r in report["results"])
    print(
        f"{'total (best)':32s} {'':>10s} {'':>3s} {'':>4s} {'':>4s} {'':>6s} {total:9.3f}"
    )


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--solver",
        default="legacy",
        help="adapter name from solvers.available_solvers()",
    )
    p.add_argument("--repeat", type=int, default=10)
    p.add_argument(
        "--tag",
        default="benchmark",
        help="instance tag to run ('test' for the small ones)",
    )
    p.add_argument("--json", type=Path, help="write the report to this file")
    p.add_argument(
        "--compare", type=Path, help="baseline report (json) to compare against"
    )
    a = p.parse_args()

    report = run(a.solver, a.repeat, a.tag)
    baseline = json.loads(a.compare.read_text()) if a.compare else None
    print_table(report, baseline)
    if a.json:
        a.json.parent.mkdir(parents=True, exist_ok=True)
        a.json.write_text(json.dumps(report, indent=2) + "\n")
        print(f"wrote {a.json}")
    if not all(r["correct"] for r in report["results"]):
        sys.exit(1)


if __name__ == "__main__":
    main()
