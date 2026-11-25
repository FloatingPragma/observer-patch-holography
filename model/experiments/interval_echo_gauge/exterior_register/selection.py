"""Classify the full exterior census using only one station's decoded records."""

from collections import defaultdict


def census(received):
    if sorted(received) != list(range(32)):
        raise ValueError("The station requires each occupation word exactly once")
    groups = defaultdict(list)
    for word in received:
        colour, weak = word & 7, word >> 3
        groups[colour.bit_count(), weak.bit_count()].append((colour, weak))
    keys = sorted(
        (k for k in groups if sum(k) not in (0, 5)), key=lambda k: (sum(k), -k[0])
    )
    rows = []
    for c, w in keys:
        words = groups[c, w]
        cd, wd = len({a for a, _ in words}), len({b for _, b in words})
        charge = -2 * c + 3 * w
        rows.append(
            dict(
                colour_count=c,
                weak_count=w,
                colour_dim=cd,
                weak_dim=wd,
                charge=charge,
                conjugate=keys.index((3 - c, 2 - w)),
                occupations=sorted(a + 8 * b for a, b in words),
                anomalies=[
                    len(words) * charge,
                    wd * charge if cd == 3 else 0,
                    cd * charge if wd == 2 else 0,
                    len(words) * charge**3,
                ],
            )
        )
    return rows


def scan(rows):
    ballots = []
    for mask in range(1024):
        chosen = [r for j, r in enumerate(rows) if mask & (1 << j)]
        sums = [sum(r["anomalies"][t] for r in chosen) for t in range(4)]
        chiral = all(not mask & (1 << r["conjugate"]) for r in chosen)
        ballots.append([mask, int(bool(chosen)), int(chiral), *sums])
    return ballots


def survivors(ballots):
    return [
        mask
        for mask, nonempty, chiral, *sums in ballots
        if nonempty and chiral and not any(sums)
    ]
