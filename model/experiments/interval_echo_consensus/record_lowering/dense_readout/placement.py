"""Timing lower bounds which survive arbitrary common relabelings of the mounted owners."""

from .apparatus import SITES, squared_distance


def obstruction(arrivals):
    sites = range(len(SITES))
    target = sorted(sum(squared_distance(r, s) <= 4 for s in sites) for r in sites)
    rows = []
    for deadline in sorted({0, *(v for row in arrivals for v in row)}):
        measured = sorted(sum(v <= deadline for v in row) for row in arrivals)
        degree_error = sum(abs(a - b) for a, b in zip(measured, target, strict=True))
        asymmetric = sum(
            (arrivals[i][j] <= deadline) != (arrivals[j][i] <= deadline)
            for i in sites
            for j in range(i)
        )
        rows.append(
            dict(
                deadline=deadline,
                degree_error=degree_error,
                one_way_pairs=asymmetric,
                lower_bound=max(degree_error, asymmetric),
            )
        )
    return dict(
        target_sizes=target,
        all_deadlines=rows,
        best=min(rows, key=lambda r: (r["lower_bound"], r["deadline"])),
    )
