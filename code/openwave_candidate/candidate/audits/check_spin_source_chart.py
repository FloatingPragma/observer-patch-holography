"""Optional exact comparison with a public pinned RER checkout.

The portable package does not import RER. This audit checks its chart map.
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import sympy as s

PIN = "3bcacf1276ddb3ead35be2b9adcd63b6643e2d35"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rer", required=True, type=Path)
    args = parser.parse_args()
    if subprocess.check_output(["git", "-C", str(args.rer), "rev-parse", "HEAD"], text=True).strip() != PIN:
        raise ValueError("RER checkout differs from the specified public commit")
    subprocess.run(["git", "-C", str(args.rer), "diff", "--exit-code", PIN, "--", "model", "code"], check=True)
    sys.path.insert(0, str(args.rer.resolve()))
    from oph_openwave import spin as port
    from model.particle_primary_composite_spin import model as source
    from model.particle_primary_sector_threshold import model as sectors

    actual = source.local_multiplet((0, 2, 0, 0, 1))
    blocks = dict(sectors.native_blocks())
    weak, singlet = blocks[-3], blocks[6]
    mode_map = [15*row+i for row in (0, 1) for i in weak]+[15*row+singlet[0] for row in (0, 1)]
    masks, basis, _, _, spin, *_ = port.local_kernel()
    transform = s.zeros(12)
    where = {m: i for i, m in enumerate(actual["masks"])}
    for col, mask in enumerate(masks):
        ops = tuple((mode_map[i], True) for i in range(6) if mask & (1 << i))
        out, phase = source.apply_operators(source.OMEGA, ops)
        transform[where[out], col] = phase
    mapped = transform*basis
    projector = mapped*(mapped.H*mapped).inv()*mapped.H
    gauges_match = all(g*mapped == s.zeros(12, 6) for g in actual["gauge_matrices"])
    projector_matches = s.simplify(projector-actual["kernel_basis"]*actual["kernel_basis"].H) == s.zeros(12)
    spin_matches = [s.simplify(a*mapped-mapped*j) == s.zeros(12, 6)
                    for a, j in zip(actual["ambient_spin_generators"], spin)]
    # Reverse one exterior basis sign within a nontrivial kernel column.
    column = next(j for j in range(6) if sum(bool(basis[i, j]) for i in range(12)) > 1)
    row = next(i for i in range(12) if basis[i, column])
    broken = transform.copy()
    broken[:, row] = -broken[:, row]
    broken_rejected = any(g*broken*basis != s.zeros(12, 6) for g in actual["gauge_matrices"])
    if not (gauges_match and projector_matches and all(spin_matches) and broken_rejected):
        raise ArithmeticError("the source chart or its mutation check failed")
    print(json.dumps({"RER_commit": PIN, "mode_map": mode_map,
                      "all_twelve_native_gauges_annihilate_mapped_port_kernel": gauges_match,
                      "complete_exact_kernel_projector_matches": projector_matches,
                      "all_three_spin_intertwiners_match": spin_matches,
                      "changed_exterior_sign_rejected_by_native_gauss": broken_rejected}, indent=2))


if __name__ == "__main__":
    main()
