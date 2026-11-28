"""Measure the mixed work required for the observed head odds to be thermally stationary."""

import argparse
from hashlib import sha256
import json

from mpmath import iv
import numpy as np

from . import thermal_pair as exp

HERE = exp.HERE


def run():
    energy = json.loads((HERE / "thermal_energy.json").read_bytes())
    pair = json.loads((HERE / "pair_contact.json").read_bytes())
    bands = exp.contact.overlap_bands(pair["source"])
    odds = iv.mpf(bands["reference_log_odds"])
    logr = iv.mpf(energy["log_trace_mixed"])
    beta = exp.kernel.BETA
    needed = (logr - odds) / beta
    # In the evacuated native source the complete two-head generator is the
    # sum of the calibrated physical head generators on these disjoint seams.
    mats = exp.contact.blocks()["full"]
    empty = []
    for a, b in ((0, 0), (0, 2), (2, 0), (2, 2)):
        A, B = mats[a][0], mats[b][0]
        H = np.kron(A, np.eye(len(B))) + np.kron(np.eye(len(A)), B)
        diag = H.diagonal() / 4
        assert min(diag) == max(diag)
        empty.append(float(diag[0]))
    mixed = empty[0] - empty[1] - empty[2] + empty[3]
    assert mixed == 0
    accepted = iv.mpf(exp.probe.log_factor(pair["source"]).tolist())
    q = iv.mpf(bands["connected_charge_nats"])
    heat = iv.mpf(energy["reversible_heat_into_ordinary"])
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (HERE / "pair_contact.json", HERE / "thermal_energy.json")
        },
        empty_head_energies=empty,
        mixed_empty_energy=mixed,
        measured_reference_log_odds=bands["reference_log_odds"],
        ordinary_log_partition_mixed=energy["log_trace_mixed"],
        fixed_geometry_free_energy_slope=exp.interval_output((odds - logr) / beta),
        required_mixed_source_work=exp.interval_output(needed),
        required_work_for_change=exp.interval_output(accepted * needed),
        joint_modular_increment_nats=exp.interval_output(q + beta * heat),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    path = HERE / "thermal_stationarity.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        exp.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), result
        )
    print(result)
