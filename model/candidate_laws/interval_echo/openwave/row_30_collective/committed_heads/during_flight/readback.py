"""Read local coin-change and distant wall frequencies from complete count banks."""

import numpy as np

from ..readback import interval, subtract


def decode(banks, clocks, shots, sources):
    rows = {(b["kappa"], b["appointment"], b["menu"], b["source"]): b for b in banks}
    if len(rows) != len(banks):
        raise ValueError("duplicate detector bank")
    clocks_out = []
    last = 0
    for hand in clocks:
        hits = hand["counts"]
        estimates = 2 * np.array(hits) / shots - 1
        ci = 2 * np.array([interval(h, shots) for h in hits]) - 1
        angle = np.arctan2(estimates[1], estimates[0])
        angle += 2 * np.pi * round((last - angle) / (2 * np.pi))
        radius = np.linalg.norm(np.max(abs(ci - estimates[:, None]), axis=1))
        error = np.arcsin(radius / np.linalg.norm(estimates))
        clocks_out.append([float(angle - error), float(angle + error)])
        last = angle
    reports = []
    for k in (8, 0):
        for t in (2, 4, 6):
            responses = {}
            flip_rates = []
            conditional = {}
            for name in ("magnitude", "once", "duplicate", "spaced"):
                stands = []
                for source in sources:
                    bank = rows[k, t, name, source]
                    counts = np.array(bank["counts"], dtype=int).reshape(-1, 2)
                    if counts.sum() != shots or counts.min() < 0:
                        raise ValueError("incomplete count bank")
                    stands.append(interval(int(counts[:, 0].sum()), shots))
                    if name == "spaced":
                        flip_rates.append(interval(int(counts[[1, 2], :].sum()), shots))
                        for label, (yes, no) in zip(
                            bank["records"], counts, strict=True
                        ):
                            conditional.setdefault(label, []).append(
                                interval(int(yes), int(yes + no)) if yes + no else None
                            )
                    if name == "duplicate" and counts[[1, 2], :].sum():
                        raise ValueError("inconsistent immediate repeated record")
                responses[name] = subtract(*stands)
            reports.append(
                {
                    "kappa": k,
                    "appointment": t,
                    "source_responses": responses,
                    "spaced_minus_once": subtract(
                        responses["spaced"], responses["once"]
                    ),
                    "spaced_minus_half_once": subtract(
                        responses["spaced"], np.multiply(0.5, responses["once"])
                    ),
                    "coin_change_by_source": flip_rates,
                    "coin_change_difference": subtract(*flip_rates),
                    "conditional_spaced_response": {
                        key: subtract(*pair) if None not in pair else None
                        for key, pair in conditional.items()
                    },
                }
            )
    return {"clock_intervals": clocks_out, "comparisons": reports}
