#!/usr/bin/env python3
"""Guard the live non-hadron predictive builders against reference injection."""

from __future__ import annotations

import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[1]

PREDICTIVE_BUILDERS = [
    ROOT / "particles" / "flavor" / "derive_charged_shared_absolute_scale_writeback.py",
    ROOT / "particles" / "flavor" / "derive_quark_sector_mean_split.py",
    ROOT / "particles" / "flavor" / "derive_quark_sector_descent.py",
    ROOT / "particles" / "flavor" / "build_forward_yukawas.py",
    ROOT / "particles" / "leptons" / "derive_lepton_ordered_shape_readout.py",
    ROOT / "particles" / "leptons" / "derive_lepton_log_spectrum_readout.py",
    ROOT / "particles" / "leptons" / "build_forward_charged_leptons.py",
    ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_pair.py",
    ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_readout.py",
    ROOT / "particles" / "calibration" / "derive_d11_forward_seed.py",
    ROOT / "particles" / "calibration" / "derive_d11_live_exact_split_pair_theorem.py",
    ROOT / "particles" / "calibration" / "derive_d11_critical_surface_readout.py",
]

FORBIDDEN_SNIPPETS = [
    "particle_reference_values.json",
    "REFERENCE_JSON",
    "current_reference_exact_witness",
    "exact_wz_candidate",
    "--exact-readout",
]


def main() -> int:
    failures: list[str] = []
    for path in PREDICTIVE_BUILDERS:
        text = path.read_text(encoding="utf-8")
        for snippet in FORBIDDEN_SNIPPETS:
            if snippet in text:
                failures.append(f"{path}: contains forbidden snippet `{snippet}`")
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
