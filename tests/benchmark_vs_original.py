# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy", "scipy", "pytest"]
# ///
"""Compare the current implementation with the original CALGO 552 code.

    uv run tests/benchmark_vs_original.py                  # table + tests/benchmarks/vs_original.json
    uv run tests/benchmark_vs_original.py --quick          # fewer repetitions, no extra instances

Four builds, compiled by the same compiler with the same flags (``-O2 -ffp-contract=off`` and, as a
sensitivity check, ``-O3 -ffp-contract=off``; no fused multiply-add, so results stay reproducible):

  original       legacy/CALGO552.f            single precision, fixed form, 18 arguments
  current_single src/ with wp = c_float       (temporary copy; the repository is not touched)
  current_double src/ as it is                double precision
  frozen_double  legacy/f90_double/           the first double precision version, before the split

``current_single`` versus ``original`` isolates the effect of the refactoring (same precision, same
algorithm, so the results must be bit-identical); ``current_double`` versus ``frozen_double`` isolates the
refactoring in double precision; ``current_double`` versus ``original`` is what a user sees.

Only the Fortran call is timed (data packing in Python and the hang guard are excluded). The builds are
timed round-robin, rotating the order every repetition, so that drift of the machine (frequency,
other load) affects all of them alike. Run it on an otherwise idle machine.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import legacy_cl1 as L
import test_differential as td
from generate_instances import reference_solution
from instances import Instance, _random_problem, load_instances

FLAG_SETS = {
    "O2": ["-O2", "-ffp-contract=off", "-fPIC", "-w"],
    "O3": ["-O3", "-ffp-contract=off", "-fPIC", "-w"],
}
BUILDS = ["original", "current_single", "frozen_double", "current_double"]

SELECTED = [  # persisted instances: sizes from 500x10 to 5000x30, with and without constraints
    "bench_500x10",
    "bench_2000x20",
    "bench_2000x20_l3_m10",
    "bench_300x100",
    "bench_5000x30",
    "bench_cvxpy_script_1000x250",
    "rand_300x15_l3_m20",
    "kode1_signs_and_constraints",
]


class _Timed:
    """Wraps the ctypes function so that only the foreign call is timed."""

    def __init__(self, fn):
        self.fn, self.elapsed = fn, 0.0

    def __call__(self, *args):
        start = time.perf_counter()
        self.fn(*args)
        self.elapsed = time.perf_counter() - start


def make_single_sources(tmp: Path) -> list[Path]:
    """A temporary copy of src/ whose working precision is single (the repository is not modified)."""
    out = []
    for name in ("l1_precision.f90", "l1_calgo552.f90", "l1_c_api.f90"):
        text = (L.ROOT / "src" / name).read_text()
        if name == "l1_precision.f90":
            text = text.replace(
                "integer, parameter :: wp = c_double",
                "integer, parameter :: wp = c_float",
            )
            text = text.replace("only: c_double", "only: c_float, c_double")
            assert "wp = c_float" in text
        (tmp / name).write_text(text)
        out.append(tmp / name)
    return out


def make_builds(flags: list[str], tmp: Path) -> dict[str, L.FortranCL1]:
    single_sources = make_single_sources(tmp)

    class CurrentSingle(L.SrcCL1):
        name = "current_single"

        def __init__(self) -> None:
            L.FortranCL1.__init__(self, single_sources, np.float32, fflags=flags)

    class CurrentDouble(L.SrcCL1):
        name = "current_double"
        fflags = flags

    builds = {
        "original": L.FortranCL1([L.DEFAULT_SRC], np.float32, fflags=flags),
        "current_single": CurrentSingle(),
        "frozen_double": L.FortranCL1(L.DOUBLE_ORACLE_FILES, np.float64, fflags=flags),
        "current_double": CurrentDouble(),
    }
    for build in builds.values():
        build._fn = _Timed(build._fn)
    return builds


def extra_instances() -> list[Instance]:
    """Larger problems that are generated on the fly (seeded), for more stable timings."""
    tag = ("benchmark",)
    return [
        _random_problem("big_20000x20", "20000 x 20", 201, 20000, 20, tags=tag),
        _random_problem(
            "big_4000x60_l5_m30",
            "4000 x 60, 5 eq, 30 ineq",
            202,
            4000,
            60,
            l=5,
            m=30,
            tags=tag,
        ),
        _random_problem("wide_1500x150", "1500 x 150", 203, 1500, 150, tags=tag),
    ]


def objective(inst: Instance, x: np.ndarray) -> float:
    return float(np.abs(inst.b - inst.A @ x).sum())


def time_instance(
    inst: Instance, builds: dict[str, L.FortranCL1], repeat: int | None
) -> dict:
    results, times = {}, {name: [] for name in builds}
    for name, build in builds.items():  # warm-up, and the result to compare
        results[name] = build.solve(inst)
    if repeat is None:
        slowest = max(build._fn.elapsed for build in builds.values())
        repeat = int(min(300, max(15, 1.0 / max(slowest, 1e-5))))
    names = list(builds)
    for rep in range(repeat):
        for name in names[rep % len(names) :] + names[: rep % len(names)]:
            builds[name].solve(inst)
            times[name].append(builds[name]._fn.elapsed)
    row = {
        "name": inst.name,
        "k": inst.k,
        "l": inst.l,
        "m": inst.m,
        "n": inst.n,
        "repeat": repeat,
        "ref_objective": inst.ref_objective,
        "builds": {},
    }
    for name in builds:
        r = results[name]
        row["builds"][name] = {
            "kode": r.kode,
            "iterations": r.iterations,
            "objective": objective(inst, r.x),
            "reported_error": r.error,
            "min_ms": 1e3 * min(times[name]),
            "median_ms": 1e3 * statistics.median(times[name]),
        }
    a, b = results["original"], results["current_single"]
    c, d = results["current_double"], results["frozen_double"]
    same = lambda p, q: (
        p.kode == q.kode
        and p.iterations == q.iterations
        and np.array_equal(p.x, q.x)
        and np.array_equal(p.res, q.res)
        and p.error == q.error
    )
    row["original_vs_current_single_identical"] = bool(same(a, b))
    row["frozen_vs_current_double_identical"] = bool(same(d, c))
    row["max_dx_original_vs_current_double"] = float(
        np.abs(a.x - c.x).max() / (1.0 + np.abs(c.x).max())
    )
    return row


def corpus_identity(builds: dict[str, L.FortranCL1]) -> dict:
    """Bit-for-bit comparison over every persisted instance and the whole differential fuzz corpus."""
    cases = [(i, None) for i in load_instances()] + [
        td.GENERATORS[g](s) for g, s in td.fuzz_corpus()
    ]
    n_single = n_double = 0
    for inst, toler in cases:
        a, b = (
            builds["original"].solve(inst, toler=toler),
            builds["current_single"].solve(inst, toler=toler),
        )
        c, d = (
            builds["current_double"].solve(inst, toler=toler),
            builds["frozen_double"].solve(inst, toler=toler),
        )
        for p, q, which in ((a, b, "s"), (c, d, "d")):
            ok = (
                (p.kode, p.iterations, p.error) == (q.kode, q.iterations, q.error)
                and np.array_equal(p.x, q.x)
                and np.array_equal(p.res, q.res)
            )
            n_single += ok and which == "s"
            n_double += ok and which == "d"
    return {
        "cases": len(cases),
        "original_vs_current_single_identical": n_single,
        "frozen_vs_current_double_identical": n_double,
    }


def geomean(values: list[float]) -> float:
    return math.exp(sum(math.log(v) for v in values) / len(values))


def print_report(flag_name: str, rows: list[dict]) -> None:
    print(f"\n=== flags {flag_name}: {' '.join(FLAG_SETS[flag_name])} ===")
    print(
        "run time of the Fortran call, best of N (ms); ratios < 1 mean the current code is faster"
    )
    head = (
        f"{'instance':28s} {'k x n':>10s} {'l':>3s} {'m':>3s} | {'orig':>8s} {'cur.sgl':>8s} {'frozen':>8s} {'cur.dbl':>8s} | "
        f"{'sgl/orig':>8s} {'dbl/frozn':>9s} {'dbl/orig':>8s} | iterations (orig/cur.sgl/frozen/cur.dbl)"
    )
    print(head)
    ratios = {"sgl/orig": [], "dbl/frozen": [], "dbl/orig": []}
    for r in rows:
        b = r["builds"]
        t = {k: b[k]["min_ms"] for k in BUILDS}
        s, f, d = (
            t["current_single"] / t["original"],
            t["current_double"] / t["frozen_double"],
            t["current_double"] / t["original"],
        )
        ratios["sgl/orig"].append(s)
        ratios["dbl/frozen"].append(f)
        ratios["dbl/orig"].append(d)
        its = "/".join(str(b[k]["iterations"]) for k in BUILDS)
        print(
            f"{r['name']:28s} {str(r['k']) + 'x' + str(r['n']):>10s} {r['l']:3d} {r['m']:3d} | {t['original']:8.3f} {t['current_single']:8.3f} "
            f"{t['frozen_double']:8.3f} {t['current_double']:8.3f} | {s:8.2f} {f:9.2f} {d:8.2f} | {its}"
        )
    print(
        f"{'geometric mean':28s} {'':>10s} {'':>3s} {'':>3s} | {'':>8s} {'':>8s} {'':>8s} {'':>8s} | "
        f"{geomean(ratios['sgl/orig']):8.2f} {geomean(ratios['dbl/frozen']):9.2f} {geomean(ratios['dbl/orig']):8.2f} |"
    )


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--quick", action="store_true", help="few repetitions, no extra instances"
    )
    p.add_argument(
        "--flags", nargs="+", default=list(FLAG_SETS), choices=list(FLAG_SETS)
    )
    p.add_argument(
        "--json",
        type=Path,
        default=L.ROOT / "tests" / "benchmarks" / "vs_original.json",
    )
    a = p.parse_args()

    by_name = {i.name: i for i in load_instances()}
    instances = [by_name[n] for n in SELECTED]
    if not a.quick:
        for inst in extra_instances():
            inst.ref_status, inst.ref_objective, inst.ref_x = reference_solution(inst)
            instances.append(inst)

    fortran = subprocess.run(
        ["gfortran", "--version"], capture_output=True, text=True
    ).stdout.splitlines()[0]
    cpu = (
        subprocess.run(
            ["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True
        ).stdout.strip()
        or platform.processor()
    )
    report = {
        "machine": platform.platform(),
        "cpu": cpu,
        "fortran": fortran,
        "python": platform.python_version(),
        "runs": {},
    }
    print(f"{cpu} | {fortran} | numpy {np.__version__}")

    tmp = Path(tempfile.mkdtemp(prefix="l1fit_bench_"))
    try:
        for flag_name in a.flags:
            builds = make_builds(FLAG_SETS[flag_name], tmp)
            rows = [
                time_instance(inst, builds, 15 if a.quick else None)
                for inst in instances
            ]
            print_report(flag_name, rows)
            run = {"flags": FLAG_SETS[flag_name], "instances": rows}
            if flag_name == a.flags[0]:
                run["corpus_identity"] = corpus_identity(builds)
                print(
                    "\nbit-for-bit comparison over the instances and the whole fuzz corpus:",
                    run["corpus_identity"],
                )
            report["runs"][flag_name] = run
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    a.json.parent.mkdir(parents=True, exist_ok=True)
    a.json.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {a.json}")


if __name__ == "__main__":
    main()
