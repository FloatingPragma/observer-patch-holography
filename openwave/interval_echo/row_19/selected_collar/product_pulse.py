"""Pass the independent-head calibration through the same timing and head instruments."""

import argparse
import json
from hashlib import sha256

import numpy as np

from . import returned_head_lab as lab


def run():
    state = lab.heads.product_reference(lab.sources()["reference"])
    ports, history, err = lab.pulses(np.array([state]))
    rng = np.random.default_rng(61432879)
    probs = np.trace(history[0], axis1=-2, axis2=-1).real.sum(axis=2)
    tape = {
        "time": [lab.counts(np.r_[p, 1 - p.sum()], rng) for p in probs],
        "incoming": lab.pairs(state, rng, False),
        "returned": [lab.pairs(s, rng, True) for s in history[0].sum(axis=1)],
    }
    return {
        "protocol_sha256": sha256(
            (lab.HERE / "TRANSMITTED_REFERENCE.md").read_bytes()
        ).hexdigest(),
        "source_bank_sha256": sha256(
            (lab.HERE / "returned_head_lab.json").read_bytes()
        ).hexdigest(),
        "seed": 61432879,
        "shots_per_setting": lab.SHOTS,
        "ports": ports,
        "trace_effect_max_difference": err,
        "tape": tape,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    result = run()
    dest = lab.HERE / "product_pulse.json"
    if cli.parse_args().write:
        dest.write_text(
            json.dumps(result, allow_nan=False, separators=(",", ":")) + "\n"
        )
    else:
        stored = json.loads(dest.read_bytes())
        assert result["tape"] == stored["tape"] and result["ports"] == stored["ports"]
        for key in (
            "protocol_sha256",
            "source_bank_sha256",
            "shots_per_setting",
            "seed",
        ):
            assert result[key] == stored[key]
    print(result["ports"], result["trace_effect_max_difference"])
