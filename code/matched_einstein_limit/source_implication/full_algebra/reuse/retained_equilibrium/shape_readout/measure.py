"""Reconstruct the measured edge response and certify flat realizability."""

import hashlib
import json
from pathlib import Path

from ..observer_banks import audit

HERE = Path(__file__).resolve().parent
BANK = audit.HERE / "records.json"
iv = audit.iv


def edge_logs(ranges):
    log_times = [iv.ln(x) for x in ranges]
    return [
        sum(p * t for p, t in zip(row, log_times, strict=True)) / 2
        for row in audit.range_powers()
    ]


def principal_minors(log_lengths):
    a, b, c, d, e, f = [iv.exp(2 * x) for x in log_lengths]
    p, q, r = (a + b - d) / 2, (a + c - e) / 2, (b + c - f) / 2
    return (
        a,
        a * b - p * p,
        a * b * c + 2 * p * q * r - a * r * r - b * q * q - c * p * p,
    )


def measure(precision=45):
    bank_bytes = BANK.read_bytes()
    bank = json.loads(bank_bytes)
    previous = iv.dps
    try:
        iv.dps = precision
        audit.marginal_check(bank)
        _, width = audit.instrument_widths()
        logs = [
            edge_logs(
                [t + iv.mpf([-1, 1]) * width for t in audit.clock_means(arm["clock"])]
            )
            for arm in bank["counts"]
        ]
        secants = [audit.rounded((b - a) * 50000) for a, b in zip(*logs, strict=True)]
        # Work on the published outward grid. This is a conservative lower bound.
        spread = max(lo for lo, _ in secants) - min(hi for _, hi in secants)
        return {
            "bank_sha256": hashlib.sha256(bank_bytes).hexdigest(),
            "units": "integer millionths, outward rounded",
            "edge_log_secants": {
                f"{a}-{b}": interval
                for (a, b), interval in zip(audit.EDGES, secants, strict=True)
            },
            "isotropic_max_error_lower": max(0, spread // 2),
            "gram_principal_minors": [
                [audit.rounded(x) for x in principal_minors(endpoint)]
                for endpoint in logs
            ],
        }
    finally:
        iv.dps = previous


if __name__ == "__main__":
    print(json.dumps(measure(), indent=2))
