"""Compare regional source readings and the area of closed marker surveys."""

import numpy as np

from . import closed_balls


def adaptive_scores(words, mask_by_head, labels, refs):
    # The observer records which markers are in the ball and both oriented
    # boundary counts, as well as every enclosed load and enclosed head.
    mask = mask_by_head[words[:, -1]]
    tag = np.column_stack(
        (labels[words[:, -1]], np.where((mask[:, None] >> np.arange(7)) & 1, words, -1))
    )
    _, idx = np.unique(tag, axis=0, return_inverse=True)
    return np.array([-np.log(np.bincount(idx, weights=v)[idx]) for v in refs])


def adaptive(words, charts, sources, refs, hot, cold):
    head = words[:, -1]
    cache = {}
    rows = {}
    for (load, desc, r), panel in charts.items():
        by_head = panel["areas"]
        record = panel["records"]
        labels = panel["geometry_labels"]
        gain = hot - cold
        change = {}
        fixed = (
            panel["fixed_members"]
            & panel["fixed_boundary_count"]
            & (record == 127).all(0)
            & (panel["volumes"] < 240)
        )
        for div, (d, err) in sources.items():
            geom = np.tensordot(
                np.bincount(head, weights=d, minlength=11), by_head, axes=(0, 0)
            )
            survey = []
            for center in range(len(fixed)):
                y = labels[:, center]
                mask = record[:, center]
                key = (tuple(y), tuple(mask))
                if key not in cache:
                    cache[key] = adaptive_scores(words, mask, y, refs)
                scores = cache[key]
                tangent = scores @ d
                ds = scores[0] @ gain
                # Outward and inward counts are distinct ordinary instruments.
                da = gain @ by_head[head, center]
                det = da * tangent[0] - geom[center] * ds
                survey.append((*(tangent), *(det)))
            ts = np.array(survey)
            best = int(np.unravel_index(np.argmax(abs(ts[:, 2:])), ts[:, 2:].shape)[0])
            change[str(div)] = dict(
                max_area_change=np.max(abs(geom), axis=0).tolist(),
                sharp_contrast_range=[
                    float(min(ts[:, 0] - ts[:, 1])),
                    float(max(ts[:, 0] - ts[:, 1])),
                ],
                largest_ordinary_determinant=dict(
                    center=best,
                    values=ts[best, 2:].tolist(),
                    held_retained_entropy=float(ts[best, 0]),
                    held_area=geom[best].tolist(),
                ),
            )
        ids = np.flatnonzero(fixed)
        witness = None
        if len(ids):
            j = int(ids[0])
            witness = dict(
                center=j,
                volume=int(panel["volumes"][j]),
                boundary_out_in=by_head[0, j].tolist(),
            )
        rows[f"{load}/{'descending' if desc else 'ascending'}/{r}"] = dict(
            count_fixed_proper_balls_with_complete_source=int(sum(fixed)),
            count_source_enclosing=int(
                ((record == 127).all(0) & (panel["volumes"] < 240)).sum()
            ),
            fixed_witness=witness,
            rates=change,
        )
    return rows


def lapse(hops, lookup, sharp, coherent, clock, observed, sources, stat, words):
    board, _, _, terms, private, _ = closed_balls.energy.returned(words)
    out = {}
    for div, (d, error) in sources.items():
        runs = {}
        ordinary_energy = closed_balls.energy_at_ports(
            board, terms, private, 1 + len(d) * d
        )
        for label, operators in (("sharp", sharp), ("coherent", coherent)):
            predicted = np.zeros(128)
            bands = np.zeros((128, 2))
            for mask, v in operators.items():
                predicted[mask] = d @ v
                eps = error * max(abs(v)) + 2e-9
                bands[mask] = [predicted[mask] - eps, predicted[mask] + eps]
            # These two spectral inversions use numerical Gibbs reductions.
            # The independent sharp experiment below carries the statistical
            # and thermometer uncertainties through the same inverse.
            runs[label] = closed_balls.geometry_test(
                hops, lookup, bands, predicted, ordinary_energy
            )["comparisons"]
        if div in observed:
            d_lab = observed[div]
            bands = np.zeros((128, 2))
            scores = np.zeros(128)
            for mask, v in clock.items():
                mid = v.mean(1)
                width = (
                    stat[div] * (v[:, 1].max() - v[:, 0].min())
                    + abs(d_lab) @ (np.ptp(v, axis=1) / 2)
                    + 1e-11
                )
                scores[mask] = d_lab @ mid
                bands[mask] = [scores[mask] - width, scores[mask] + width]
            runs["independent_records"] = closed_balls.geometry_test(
                hops,
                lookup,
                bands,
                [float(d @ sharp[m]) if m in sharp else 0.0 for m in range(128)],
                ordinary_energy,
            )["comparisons"]
        out[str(div)] = runs
    return out
