"""Nested clock readings from the two previously commissioned endpoint banks."""

import json

from ..observer_banks import audit
from .kernels import CAPS


def read():
    iv = audit.iv
    iv.dps = 42
    bank = json.loads((audit.HERE / "records.json").read_text())
    audit.marginal_check(bank)
    panels = {}
    for cap in CAPS:
        # All starts survive in the denominator. Map every later return and
        # the full timeout onto the score for the new deadline.
        readings = []
        for endpoint in bank["counts"]:
            rows = []
            for raw in endpoint["clock"]:
                first = 2 * cap - 2
                row = list(map(int, raw))
                row[first] = sum(row[first:])
                row[first + 1 :] = [0] * (11 - first)
                assert len(row) == 12 and sum(row) == audit.N
                rows.append(row)
            # Confidence for the archived ticket probabilities. All 72 score
            # means are covered at once even though they share raw returns.
            width = (iv.mpf(cap) / 2 - iv.mpf(1) / 2) * iv.sqrt(
                iv.ln(14400) / (2 * audit.N)
            )
            readings.append(audit.geometry_box(rows, width))
        slopes = [(b - a) * 50000 for a, b in zip(*readings, strict=True)]
        panels[str(cap)] = {
            s: audit.rounded(v) for s, v in zip(("area", "volume"), slopes, strict=True)
        }
    return dict(denominator=10**6, intervals=panels)


if __name__ == "__main__":
    print(json.dumps(read(), indent=2))
