"""Reproduce a bounded material clock survey at arbitrary-width addresses."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from level12_clocks import EPOCHS, LIMIT, survey, summarize
from strict_history import archive, wide_replay


def measure(level=202, digits=110, seed=archive.SEED):
    model = wide_replay.WideReplay(level, digits, seed)
    locations = [0] + [(face + 1) * 4**level - 1 for face in range(20)]
    stations = []
    for owner in dict.fromkeys(locations):
        history = model.history(owner)
        direct = []
        for port, peer, q in history.neighbors:
            incoming = model.links(12 * owner + port)
            back = model.links(12 * peer + q)
            closure = dict(incoming)[12 * peer + q]
            if dict(back)[12 * owner + port] != closure:
                raise ValueError("Source and receiver disagree on contact history")
            direct.append(
                dict(
                    port=port,
                    peer=str(peer),
                    peer_port=q,
                    local_load=2 * ((history.mask >> port) & 1),
                    peer_load=2 * ((model.history(peer).mask >> q) & 1),
                    settled_direct_echo_ticks=None if closure else 2,
                )
            )
        native = survey(model, owner)
        control = survey(model, owner, True)
        for row in native + control:
            row["peer"] = str(row["peer"])
        stations.append(
            dict(
                observer=str(owner),
                direct=direct,
                native=native,
                transparent=control,
                pulse_and_radar=summarize(native, (0, 16, 96)),
            )
        )
    return dict(
        level=level,
        sphere_observers=str(20 * 4**level),
        seed=seed,
        decimal_digits=digits,
        key_scheme=wide_replay.SCHEME,
        emissions=list(EPOCHS),
        deadline_per_leg_ticks=LIMIT,
        clock="Equal local pulse cadence. One tick per material transfer. Each challenge triggers an immediate echo; incoming pulses and challenges use separate matched trials.",
        enrollment="Observer zero and the central descendant of each base face. Every station probes port zero of its three immediate peers from port zero. Direct controls use the facing ports. Address strings are generator and matching keys.",
        sources_sha256={
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (
                Path(wide_replay.__file__),
                Path(__file__).with_name("level12_clocks.py"),
                archive.HERE / "local_replay.py",
                archive.HERE / "seams.json",
            )
        },
        paired_shifts={
            t: dict(
                Counter(
                    z for s in stations for z in s["pulse_and_radar"][t]["pulse_shifts"]
                )
            )
            for t in (0, 16, 96)
        },
        stations=stations,
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--level", type=int, default=202)
    ap.add_argument("--digits", type=int, default=110)
    ap.add_argument("--seed", type=int, default=archive.SEED)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    args.output.write_text(
        json.dumps(measure(args.level, args.digits, args.seed), indent=1) + "\n"
    )
