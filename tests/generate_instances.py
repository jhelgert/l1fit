# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy", "scipy"]
# ///
"""(Re)generate ``tests/instances/*.npz`` and ``tests/instances/manifest.json``.

For every instance this stores
  * the problem data,
  * an independent float64 *reference* optimum from SciPy/HiGHS (``ref_*``).

Every instance is also solved with the Fortran sources (``fortran_src.py``) to check that the expected exit
code and the reference optimum are consistent with the solver.

Run from the repository root::

    uv run tests/generate_instances.py

The persisted files are what the tests read; regenerate only when the instance definitions in
``instances.py`` change, and review the manifest diff when you do.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, str(Path(__file__).resolve().parent))

from instances import INSTANCE_DIR, Instance, build_all
from fortran_src import SrcCL1


def reference_solution(inst: Instance) -> tuple[str, float, np.ndarray | None]:
    """Solve the instance as an LP with HiGHS: variables (x, t), minimise sum(t), |Ax-b| <= t."""
    k, n = inst.k, inst.n
    c = np.concatenate([np.zeros(n), np.ones(k)])
    eye = np.eye(k)
    rows = [np.hstack([inst.A, -eye]), np.hstack([-inst.A, -eye])]
    rhs = [inst.b, -inst.b]
    if inst.m:
        rows.append(np.hstack([inst.E, np.zeros((inst.m, k))]))
        rhs.append(inst.f)
    if inst.ressign is not None:  # residual r = b - A x restricted in sign
        pos = np.flatnonzero(inst.ressign > 0)  # r >= 0  <=>  A x <= b
        neg = np.flatnonzero(inst.ressign < 0)  # r <= 0  <=>  -A x <= -b
        if pos.size:
            rows.append(np.hstack([inst.A[pos], np.zeros((pos.size, k))]))
            rhs.append(inst.b[pos])
        if neg.size:
            rows.append(np.hstack([-inst.A[neg], np.zeros((neg.size, k))]))
            rhs.append(-inst.b[neg])
    a_eq = np.hstack([inst.C, np.zeros((inst.l, k))]) if inst.l else None
    b_eq = inst.d if inst.l else None
    bounds = [(None, None)] * n + [(0, None)] * k
    if inst.xsign is not None:
        for j, s in enumerate(inst.xsign):
            bounds[j] = (None, 0) if s < 0 else (0, None) if s > 0 else (None, None)
    r = linprog(
        c,
        A_ub=np.vstack(rows),
        b_ub=np.concatenate(rhs),
        A_eq=a_eq,
        b_eq=b_eq,
        bounds=bounds,
        method="highs",
    )
    status = {0: "optimal", 2: "infeasible", 3: "unbounded"}.get(
        r.status, f"status{r.status}"
    )
    if r.status != 0:
        return status, float("nan"), None
    return status, float(r.fun), r.x[:n]


def main() -> None:
    solver = SrcCL1()
    manifest = []
    for inst in build_all():
        status, obj, x = reference_solution(inst)
        inst.ref_status, inst.ref_objective, inst.ref_x = status, obj, x

        expect_infeasible = inst.expected_kode == 1
        if expect_infeasible != (status == "infeasible"):
            raise SystemExit(
                f"{inst.name}: reference status '{status}' contradicts expected_kode={inst.expected_kode}"
            )

        t0 = time.perf_counter()
        res = solver.solve(inst)
        elapsed = time.perf_counter() - t0
        if res.kode != inst.expected_kode:
            raise SystemExit(f"{inst.name}: solver KODE={res.kode}, expected {inst.expected_kode}")
        inst.save()

        manifest.append(
            {
                "name": inst.name,
                "description": inst.description,
                "tags": list(inst.tags),
                "k": inst.k,
                "l": inst.l,
                "m": inst.m,
                "n": inst.n,
                "ref_status": status,
                "ref_objective": None if np.isnan(obj) else obj,
            }
        )
        rel = (
            abs(res.error - obj) / max(1.0, abs(obj))
            if not np.isnan(obj)
            else float("nan")
        )
        print(
            f"{inst.name:32s} {inst.k:5d}x{inst.n:<4d} l={inst.l:<3d} m={inst.m:<4d} "
            f"ref={obj:12.6f} solver={res.error:12.6f} (rel diff {rel:.1e}) "
            f"kode={res.kode} iter={res.iterations:6d} {elapsed * 1e3:8.1f} ms"
        )

    (INSTANCE_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"\nwrote {len(manifest)} instances to {INSTANCE_DIR}")


if __name__ == "__main__":
    main()
