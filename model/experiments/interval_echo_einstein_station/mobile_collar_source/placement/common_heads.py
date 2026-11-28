"""Finish the two held-out finite populations of independently armed original heads."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import extended
from . import measure as m


def acquire():
    path = m.HERE / "extended.json"
    earlier = json.loads(path.read_text())
    radius = m.uncertainty()["covariance"]
    result = {}
    for level in (0, 2):
        b = m.source.federation(level)
        rec = earlier["boards"][str(b.cells)]
        rosters = list(m.routes(b, level == 0))
        charges = {}
        for q in (-1, 1):
            bank = rec["charges"][str(q)]
            known = {
                tuple(rec["minimum_routes"][i]): raw
                for i, raw in zip(
                    rec["standard_indices"], bank["standards"], strict=True
                )
            }
            transmitting = bank["survey_p"] if level == 0 else bank["census"]
            # Store only the genuinely new sharp-input probabilities.
            new = {}
            hydrated = []
            for i, r in enumerate(rosters):
                if r in known:
                    hydrated.append(known[r])
                    continue
                initial = m.initial_index(r, q)
                h = m.scaled_h(b, r, q)
                other = [
                    m.propagate_scaled(h, r, s).tolist()
                    for s in range(4)
                    if s != initial
                ]
                full = other.copy()
                full.insert(initial, transmitting[i])
                new[str(i)] = other
                hydrated.append(full)
            all_p = np.asarray(hydrated)
            summary = {
                str(j): extended.summarize(all_p[:, j], len(rosters), radius)
                for j in range(4)
            }
            summary["independent_fair"] = extended.summarize(
                all_p.mean(axis=1), len(rosters), radius
            )
            charges[str(q)] = dict(new_nontransmitting_inputs=new, summary=summary)
        result[str(b.cells)] = dict(total=len(rosters), charges=charges)
    return dict(
        extended_sha256=sha256(path.read_bytes()).hexdigest(),
        protocol_sha256=sha256((m.HERE / "COMMON_HEADS.md").read_bytes()).hexdigest(),
        boards=result,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    d = acquire()
    path = m.HERE / "common_heads.json"
    if args.write:
        path.write_text(json.dumps(d, separators=(",", ":"), allow_nan=False) + "\n")
    else:
        m.compare(json.loads(path.read_text()), d)
    for n, b in d["boards"].items():
        print(
            n,
            {
                q: {key: row[-1]["mean"] for key, row in v["summary"].items()}
                for q, v in b["charges"].items()
            },
        )
