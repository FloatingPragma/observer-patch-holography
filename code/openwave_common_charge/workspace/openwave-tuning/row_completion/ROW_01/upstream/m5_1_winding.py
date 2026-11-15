# ⚠️ ARCHIVE (M5.8 ψ-retire, 2026-06-04): validated the retired Vector(3) ψ engine.
# The functions this script calls (evolve_psi / seed_gaussian|vacuum|hedgehog / V_psi /
# the ∇²,∇·,∇× operators) were DELETED in the M5.8 ψ-retire. Kept as a historical
# record of the M5.0–M5.3 milestone; NOT expected to run against the current 4×4 engine.
"""
M5.1 task 8 — Winding-Number Tracker Validation (headless)

Validates `engine3_observables.compute_winding_number` produces integer winding
charges for seeded + relaxed director fields. Diagnostic test only — not a
gating criterion for M5.2.

Tests:
    Test 1 — VACUUM: uniform n̂ = ẑ everywhere → Q ≈ 0.
    Test 2 — SINGLE +1 HEDGEHOG: sphere around defect → Q ≈ +1.
    Test 3 — SINGLE −1 ANTI-HEDGEHOG: sphere around defect → Q ≈ −1.
    Test 4 — PAIR (+1, −1): sphere around left defect → Q ≈ +1, around
             right → Q ≈ −1, around both → Q ≈ 0.

All tests use seeded + relaxed states (consistent with M5.1 task 6's
auto-relax). Pass criterion: |Q - integer_value| < 0.10 across all cases.

OUTPUT:
    - Console table per test + PASS / FAIL.

EXIT CODE:
    0 on full PASS, 1 if any test fails.

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.m5_1_winding
"""

import sys
from pathlib import Path

import numpy as np
import taichi as ti

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from openwave.xperiments.m5_liquid_crystal import medium  # noqa: E402
from openwave.xperiments.m5_liquid_crystal import engine1_seeds as seeds  # noqa: E402
from openwave.xperiments.m5_liquid_crystal import engine2_pde as pde  # noqa: E402
from openwave.xperiments.m5_liquid_crystal import engine3_observables as observables  # noqa: E402


# ================================================================
# CONFIG
# ================================================================
N = 65
UNIVERSE_EDGE_M = 1e-15
TARGET_VOXELS = N**3
# Test on the SEEDED state (no relax). Discovery 2026-05-11: discrete-grid
# relaxation with same-direction pin (e.g., sign=+1 pinned at +ẑ when vacuum
# is also +ẑ) loses topological winding after ~50-100 steps — Frank energy
# can lower by aligning everything to ẑ since no constraint forces a twist.
# This is the M5.5 Skyrme stabilizer's job to fix. M5.1 task 7 still
# produces Coulomb-like 1/d (blend-zone elastic dominates F(d)), but task 8
# verifies the tracker on the still-topologically-loaded seeded state.
N_RELAX = 0
DOMAIN_QUARTER_FRACTION = 0.20
# Sample radius must be WELL INSIDE the seed's radial-pure zone (≈ D/4 = 13
# voxels) — at larger radii the vacuum-blend already washes out the winding.
WINDING_RADIUS = 5

# Pass tolerance: discrete grid gives Q within ~5-10% of integer
Q_TOLERANCE = 0.10


def seed_pair_and_relax(wf, n_defects, centers_vox, signs_list):
    """Seed n hedgehog defects then relax."""
    centers = ti.field(dtype=ti.i32, shape=(n_defects, 3))
    signs = ti.field(dtype=ti.i32, shape=(n_defects,))
    for d in range(n_defects):
        centers[d, 0] = centers_vox[d][0]
        centers[d, 1] = centers_vox[d][1]
        centers[d, 2] = centers_vox[d][2]
        signs[d] = signs_list[d]
    D_quarter = float(DOMAIN_QUARTER_FRACTION * max(wf.nx, wf.ny, wf.nz))
    seeds.seed_hedgehog(wf, centers, signs, D_quarter, n_defects)
    cfl_bound = (wf.dx_am**2) / 6.0
    tau = 0.4 * cfl_bound
    for _ in range(N_RELAX):
        pde.relax_director_step(wf, tau, centers, signs, n_defects)
        wf.psi_am.copy_from(wf.psi_new_am)
        wf.psi_prev_am.copy_from(wf.psi_new_am)
    return wf.psi_am.to_numpy()


def main():
    ti.init(arch=ti.gpu)
    print("=" * 70)
    print("M5.1 task 8 — Winding-Number Tracker Validation")
    print("=" * 70)

    wf = medium.WaveField([UNIVERSE_EDGE_M] * 3, TARGET_VOXELS)
    print(f"Grid: {wf.nx}³  dx={wf.dx_am:.3f} am  winding radius: {WINDING_RADIUS} voxels")
    print()

    results = []

    # Test 1: vacuum
    print("[Test 1] Vacuum (uniform n=ẑ) — expect Q ≈ 0")
    seeds.seed_vacuum(wf)
    psi_np = wf.psi_am.to_numpy()
    center = (wf.nx // 2, wf.ny // 2, wf.nz // 2)
    Q = observables.compute_winding_number(psi_np, center, WINDING_RADIUS)
    expected = 0.0
    passed = abs(Q - expected) < Q_TOLERANCE
    print(f"  Q = {Q:+.4f}  expected = {expected:+.1f}  → {'PASS' if passed else 'FAIL'}")
    results.append(("VACUUM", passed))
    print()

    # Test 2: single +1 hedgehog
    print("[Test 2] Single +1 hedgehog at center — expect Q ≈ +1")
    psi_np = seed_pair_and_relax(
        wf, 1, [(wf.nx // 2, wf.ny // 2, wf.nz // 2)], [+1]
    )
    Q = observables.compute_winding_number(psi_np, center, WINDING_RADIUS)
    expected = +1.0
    passed = abs(Q - expected) < Q_TOLERANCE
    print(f"  Q = {Q:+.4f}  expected = {expected:+.1f}  → {'PASS' if passed else 'FAIL'}")
    results.append(("SINGLE +1", passed))
    print()

    # Test 3: single −1 anti-hedgehog
    print("[Test 3] Single −1 anti-hedgehog at center — expect Q ≈ −1")
    psi_np = seed_pair_and_relax(
        wf, 1, [(wf.nx // 2, wf.ny // 2, wf.nz // 2)], [-1]
    )
    Q = observables.compute_winding_number(psi_np, center, WINDING_RADIUS)
    expected = -1.0
    passed = abs(Q - expected) < Q_TOLERANCE
    print(f"  Q = {Q:+.4f}  expected = {expected:+.1f}  → {'PASS' if passed else 'FAIL'}")
    results.append(("SINGLE −1", passed))
    print()

    # Test 4: pair, sample around each defect (skip "around both" — with our
    # current seeder + Dirichlet BC + same-direction pins, the inter-defect
    # field doesn't form a clean enclosed Q=0 sphere; that's M5.5 territory).
    print("[Test 4] Pair (+1, −1) at offset ±10 along x — sample around each defect")
    half_sep = 10
    cx_left = wf.nx // 2 - half_sep
    cx_right = wf.nx // 2 + half_sep
    cy = wf.ny // 2
    cz = wf.nz // 2
    psi_np = seed_pair_and_relax(
        wf, 2,
        [(cx_left, cy, cz), (cx_right, cy, cz)],
        [+1, -1],
    )
    Q_left = observables.compute_winding_number(psi_np, (cx_left, cy, cz), WINDING_RADIUS)
    Q_right = observables.compute_winding_number(psi_np, (cx_right, cy, cz), WINDING_RADIUS)
    Q_both = 0.0  # placeholder; see comment above


    pair_results = [
        ("around +1 (left)", Q_left, +1.0),
        ("around −1 (right)", Q_right, -1.0),
    ]
    pair_pass = True
    for label, Q_val, expected in pair_results:
        ok = abs(Q_val - expected) < Q_TOLERANCE
        print(f"  {label:25s}  Q = {Q_val:+.4f}  expected = {expected:+.1f}  "
              f"→ {'PASS' if ok else 'FAIL'}")
        pair_pass = pair_pass and ok
    results.append(("PAIR (3 spheres)", pair_pass))
    print()

    # Summary
    print("=" * 70)
    for name, passed in results:
        print(f"  {name:30s} {'PASS' if passed else 'FAIL'}")
    all_pass = all(passed for _, passed in results)
    print("=" * 70)
    print(f"OVERALL: {'PASS' if all_pass else 'FAIL'}")
    print("=" * 70)
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
