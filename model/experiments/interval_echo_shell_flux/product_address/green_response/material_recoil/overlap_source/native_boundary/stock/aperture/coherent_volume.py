"""Test regional entropy stationarity with independent simulated interference records."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import coherent_chart as chart
from . import coherent_control as control
from . import regional_tomography as meter
from . import unlabelled_material

BETA = np.log(8 / 7)
SEED = 11076439


def acquire_material(board, edges, clamp, load, rng):
    cells = {}
    report = {}
    for c in (0, 1):
        twice, marks = chart.flight.apparatus(board, edges, load, clamp, c)
        for s, (a, b) in enumerate(chart.flight.SECTORS):
            h = twice / 2 + np.diag(a * marks[0] + b * marks[1])
            e, v = np.linalg.eigh(h)
            p = np.exp(-BETA * (e - e[0]))
            rho = (v * (p / p.sum())) @ v.T
            rec, error, receipt = meter.acquire(rho, rng)
            cells[c, s] = h, e, v, rho, rec, error, receipt["entropy_error"]
            report[f"{c},{s}"] = receipt
    return cells, report


def panel_geometry(survey, volume, reverse, rows):
    masks, ids, counts = [], {}, []
    members = []
    for c, _, _, f in rows:
        mask, area = survey[volume, reverse, c, f]
        key = mask.tobytes()
        if key not in ids:
            ids[key] = len(masks)
            masks.append(mask)
        members.append(ids[key])
        counts.append(area)
    return masks, np.array(members), np.array(counts)


def control_panel(
    cells, rows, exact, obs, weight_errors, masks, members, live, active, ports, rng
):
    stations = np.flatnonzero(active)
    x, y = [int(np.flatnonzero(stations == address)[0]) for address in ports]
    works = np.zeros((len(rows), 2, 2))
    base = np.empty(len(rows))
    data = np.empty((len(rows), 2, 2, len(control.AGES)))
    bound = np.empty(len(rows))
    cache, evolved = {}, {}
    a, b = np.where(np.triu(live.toarray()[np.ix_(active, active)], 1))
    for k, ((c, _, s, _), m) in enumerate(zip(rows, members, strict=True)):
        h, e, v, rho, rec, errors, _ = cells[c, s]
        mask = masks[m][active]
        cross = mask[a] ^ mask[b]
        for side, z in enumerate((x, y)):
            works[k, side] = control.work(rec, errors, h, z)
            if (c, s, z) not in evolved:
                evolved[c, s, z] = control.propagation(e, v, rho, z)
            key = (c, s, m, z)
            if key not in cache:
                cache[key] = control.readouts(
                    mask, (a[cross], b[cross]), h, rho, z, evolved[c, s, z], errors, rec
                )
            entropy, fringe, baseline, lower = cache[key]
            data[k, side] = entropy, fringe
            base[k] = baseline
            if side == 0:
                bound[k] = lower
    work = obs @ works[:, :, 0]
    width = obs @ works[:, :, 1] + weight_errors @ (
        abs(works[:, :, 0]) + works[:, :, 1]
    )
    assert np.all(work > width)
    rate = work[::-1] / sum(work)
    delivered = rate * work
    engine = np.einsum("k,s,ksta->sta", exact, rate, data)
    heat_residual = sum(rate * width) + abs(delivered[0] - delivered[1])
    entropy = engine[0, 0] - engine[1, 0]
    counts, fringe_bounds = control.measured_difference(
        exact @ base, engine[0, 1], engine[1, 1], rng
    )
    lower = rate[0] * (obs - weight_errors) @ bound
    assert entropy[0] >= lower > 0
    held = engine[0, 1] - engine[1, 1]
    assert np.all((fringe_bounds[:, 0] <= held) & (held <= fringe_bounds[:, 1]))
    return dict(
        ports=list(ports),
        rates=rate.tolist(),
        work=work.tolist(),
        work_error=width.tolist(),
        heat_mismatch_bound=heat_residual,
        entropy_slope_native=entropy.tolist(),
        initial_entropy_lower=lower,
        full_entropy_mismatch_bound=BETA * heat_residual,
        cut_fringe_counts=counts,
        cut_fringe_difference=fringe_bounds.tolist(),
        area_difference=0,
    )


def examine(board, old, edges, center, clamp, load, rows, true, rng):
    _, live, survey = chart.geometry(board, edges, clamp, center)
    assert all(mask[np.r_[clamp, np.ravel(edges)]].all() for mask, _ in survey.values())
    active = np.arange(board.n) != clamp
    controls = {
        n: control.ports(live, survey, n, clamp, edges) for n in chart.VOLUMES[:-1]
    }
    cells, receipts = acquire_material(board, edges, clamp, load, rng)
    # Independent banks for reference, ordinary perturbation and held re-pairing.
    counts = np.array([rng.multinomial(meter.N, w) for w in true])
    w = counts / meter.N
    err = meter.frequency_error(w)
    p, _ = chart.flight.sector_prior(old)
    delta = chart.flight.lab.probe.log_factor(old["source"])
    sharp_rng = np.random.default_rng([SEED, clamp, load, 1])
    sharp_counts = np.zeros(len(rows), dtype=np.int64)
    for head in range(4):
        selected = rows[:, 1] == head
        law = true[0, selected]
        sharp_counts[selected] = sharp_rng.multinomial(meter.N, law / law.sum())
    sharp = sharp_counts / meter.N
    mixed_weights = chart.flight.MIX[rows[:, 1]] * sharp
    panels = {}
    cache = {}
    for volume in chart.VOLUMES:
        for reverse in (False, True):
            masks, members, area = panel_geometry(survey, volume, reverse, rows)
            material, uncertainty = [], []
            for (c, _, s, _), m in zip(rows, members, strict=True):
                mask = masks[m][active]
                key = (c, s, mask.tobytes())
                if key not in cache:
                    rec = cells[c, s][4]
                    cache[key] = meter.entropy(rec[np.ix_(mask, mask)], mask.all())
                material.append(cache[key])
                uncertainty.append(cells[c, s][-1])
            scores, radii = meter.reference_scores(
                w[0], err[0], area, members, np.array(material), np.array(uncertainty)
            )
            bands = np.array(
                [
                    meter.shift_interval(w[k] - w[0], err[k] + err[0], scores, radii)
                    for k in (1, 2)
                ]
            )
            sharp_intervals = meter.shift_interval(
                mixed_weights, meter.frequency_error(sharp), scores, radii
            )
            forecasts = np.array(
                [meter.product_interval(delta, z) for z in sharp_intervals]
            )
            assert np.all(
                np.maximum(forecasts[:, 0], bands[1, :, 0])
                <= np.minimum(forecasts[:, 1], bands[1, :, 1])
            )
            if volume == 192 and not reverse and clamp != 0:
                np.testing.assert_array_equal(area, 52 - 8 * rows[:, 0])
            # Stationarity and the ordinary arm would convert the source-only
            # entropy forecast to this change in expected crossing count.
            implied_area = None
            if np.prod(bands[0, 1]) > 0:
                ratio = meter.product_interval(forecasts[1], 1 / bands[0, 1])
                implied_area = meter.product_interval(bands[0, 0], ratio).tolist()
            # Independent calibration of the source-to-area transfer with sharp heads.
            values = np.unique(area)
            sharp_reads = []
            for h in range(4):
                sel = rows[:, 1] == h
                probs = np.array(
                    [true[0, sel & (area == v)].sum() / p[h] for v in values]
                )
                np.testing.assert_allclose(probs.sum(), 1, atol=1e-14, rtol=0)
                sharp_reads.append(rng.multinomial(meter.N, probs / probs.sum()))
            mixes = (
                chart.flight.MIX
                @ np.array(sharp_reads)
                @ (values - values.min())
                / meter.N
            )
            mixed = mixes + np.array([-4, 4]) * meter.HOEFFDING * np.ptp(values)
            predicted = meter.product_interval(delta, mixed)
            assert max(predicted[0], bands[1, 0, 0]) <= min(
                predicted[1], bands[1, 0, 1]
            )

            def determinant(k):
                return meter.difference_interval(
                    meter.product_interval(bands[0, 0], bands[1, k]),
                    meter.product_interval(bands[1, 0], bands[0, k]),
                )

            out = dict(
                ordinary=bands[0].tolist(),
                overlap=bands[1].tolist(),
                source_forecasts=forecasts.tolist(),
                area_from_ordinary_entropy=implied_area,
                area_values=values.tolist(),
                sharp_area_counts=np.array(sharp_reads).tolist(),
                predicted_overlap_area=predicted.tolist(),
                determinant={
                    k: determinant(i).tolist()
                    for k, i in (("interference", 1), ("positions", 2))
                },
            )
            if volume in (192, chart.VOLUMES[-1]) and not reverse:
                # A separate receiver and random stream for the label-erasure control.
                unlabelled_rng = np.random.default_rng([SEED, clamp, load, volume, 2])
                out["decomposition"] = unlabelled_material.decompose(
                    rows,
                    w,
                    err,
                    np.asarray(material),
                    np.asarray(uncertainty),
                    cells,
                    BETA,
                    members,
                )
                out["without_source_labels"] = (
                    unlabelled_material.entropy_without_labels(
                        rows,
                        np.array([2 * true[0] - true[2], true[0], true[2]]),
                        masks,
                        members,
                        active,
                        cells,
                        unlabelled_rng,
                        delta,
                        chart.flight.MIX,
                    )
                )
            if volume != chart.VOLUMES[-1]:
                ports, deep = controls[volume]
                out["control"] = control_panel(
                    cells,
                    rows,
                    true[0],
                    w[0],
                    err[0],
                    masks,
                    members,
                    live,
                    active,
                    ports,
                    rng,
                )
                out["control"]["closed_neighborhood_inside_or_outside"] = deep
            panels[f"{volume},{int(reverse)}"] = out
    return dict(
        stratum_counts=counts.tolist(),
        sharp_counts=sharp_counts.tolist(),
        material_banks=receipts,
        panels=panels,
    )


def run():
    board, old, edges, _, first = chart.flight.chart()
    center, live, _ = chart.geometry(board, edges, first)
    chart.radar.certify_routing(board, live, edges)
    meter.verify_receiver(chart.phase.contact)
    selected, info = chart.ordinary_site(board, edges, center, first)
    rows, true = chart.strata(old)
    rng = np.random.default_rng(SEED)
    panels = {}
    for clamp in (first, selected):
        for species, load in (("vacancy", 0), ("excess", 2)):
            panels[f"{clamp},{species}"] = examine(
                board, old, edges, center, clamp, load, rows, true, rng
            )
    return dict(
        seed=SEED,
        shots_per_setting=meter.N,
        ages=list(control.AGES),
        sha256={
            f: sha256((chart.HERE / f).read_bytes()).hexdigest()
            for f in ("COHERENT_VOLUME.md", "pair_contact.json")
        },
        center=center,
        original_clamp=first,
        calibrated_clamp=selected,
        calibration=info,
        strata=rows.tolist(),
        results=panels,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    options = parser.parse_args()
    record = run()
    path = chart.HERE / "coherent_volume.json"
    if options.write:
        path.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    else:
        chart.flight.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), record
        )
    for key, value in record["results"].items():
        row = value["panels"]["192,0"]
        print(
            key,
            row["ordinary"],
            row["overlap"],
            row["determinant"],
            row["control"]["initial_entropy_lower"],
        )
