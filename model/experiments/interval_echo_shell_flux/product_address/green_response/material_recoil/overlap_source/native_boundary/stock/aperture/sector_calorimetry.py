"""Measure modular charge at fixed heat on shared native exchange records."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
import math

import numpy as np
from mpmath import iv

from . import peer_register_lab as lab
from . import sector_clock as clock
from . import sector_spectrum as spectrum
from . import sector_regional

thermal = spectrum.thermal
fibers = spectrum.fibers
HERE = lab.HERE


def interval(x):
    return [
        float(np.nextafter(float(x.a), -np.inf)),
        float(np.nextafter(float(x.b), np.inf)),
    ]


def ivq(x):
    q = Q(x)
    return iv.mpf(q.numerator) / q.denominator


def score_bound(z, orig, ids, price):
    d = np.array([sum(z[ids == j] - orig[ids == j]) for j in range(len(price))], object)
    return sum(ivq(v) * p for v, p in zip(d, price, strict=True))


def empirical(name, counts, ref, sid, cat, bank, expected, bands):
    # counts: one measured eight-bin spectrum histogram per independent register.
    # The stratified ordinary reference retains an equal count of each sharp word.
    shots, size = counts.shape[0], int(counts[0].sum())
    freq = counts.sum(0) / (shots * size) - ref / ref.sum()
    observed = counts @ expected / size
    values = np.array([expected[sid[cat == i][0]] for i in range(16)])
    pairs = fibers.BLANK_PAIRS
    diffs = np.array([values[j] - values[i] for i, j in pairs])
    order = np.argsort(diffs)
    vv = np.array([v for _, v in bank])
    totals = np.column_stack([vv[:, i] + vv[:, j] for i, j in pairs])
    blanks = sum(vv[:, j] for _, j in pairs)
    fixed = vv @ values - sum(
        vv[:, j] * change for (_, j), change in zip(pairs, diffs, strict=True)
    )
    limits = []
    for inds in (order, order[::-1]):
        remain = blanks.copy()
        score = fixed.copy()
        for i in inds:
            use = np.minimum(remain, totals[:, i])
            score += use * diffs[i]
            remain -= use
        assert not remain.any()
        limits.append(score / size)
    span = limits[1].max() - limits[0].min() + 1e-8
    # Empirical Bernstein for iid full-register scores. Constant protected
    # members cancel out of the range. Clock and source errors use a union bound.
    ell = math.log(2 / 0.001)
    width = (
        math.sqrt(2 * observed.var(ddof=1) * ell / shots)
        + 7 * span * ell / (3 * (shots - 1))
        + 1e-8
    )
    limits = np.array(bands)
    middle = limits.mean(1)
    estimate = float(freq @ middle)
    calibration = float(abs(freq) @ (np.ptp(limits, axis=1) / 2))
    return dict(
        observable=name,
        estimate=estimate,
        register_sampling_width=width,
        calibration_width=calibration,
        observed_register_sd=float(observed.std(ddof=1)),
        bounds=[estimate - width - calibration, estimate + width + calibration],
    )


def run():
    w, _, ids, blocks, ops = thermal.apparatus()
    _, _, _, cat, _, multiplicity, *_ = fibers.apparatus()
    rho = thermal.bath()
    native = fibers.source.records()[1]["original"]
    bank = []
    _, ens = lab.population.histogram(16, draw_bank=bank)
    z = ens["agreed"]
    d = z - native
    iv.dps = 50
    beta = iv.log(iv.mpf(8) / 7)
    _, _, _, shapes, sid, var, *_ = spectrum.design()
    by_shape = {shape: spectrum.thermal_price(shape, beta) for shape in set(shapes)}
    prices = [by_shape[s] for s in shapes]
    assert np.all(ops.reshape(6, len(w), len(w)).diagonal(axis1=1, axis2=2) @ d == 0)
    assert min(float(x.a) for x in prices) > 0
    C, V, _ = spectrum.reconstruction()
    delta_var = d @ var
    classes = sorted(set(shapes))
    reps = [shapes.index(v) for v in classes]
    coef = sum(ivq(t) * prices[j] for t, j in zip(C[:, 2], reps, strict=True))
    Qfull = coef * ivq(delta_var)
    direct = score_bound(z, native, ids, prices)
    assert float((direct - Qfull).a) <= 0 <= float((direct - Qfull).b)
    # Independent Ramsey calibrations and local comparative records.
    banks = clock.read()
    variance = clock.variance_meter()

    def combine(bands):
        out = []
        for shape in classes:
            indices = [j for j, s in enumerate(shapes) if s == shape]
            lo = max(bands[i][0] for i in indices)
            hi = min(bands[i][1] for i in indices)
            assert lo <= hi
            out.append([lo, hi])
        return out

    qbands = combine([r["centered_log_band"] for r in banks])
    vbands = combine([r["bounds"] for r in variance["groups"]])
    for j, p in enumerate(prices):
        lo, hi = qbands[classes.index(shapes[j])]
        assert lo <= float(p.a) <= float(p.b) <= hi
    reference_counts = {}
    var_class = np.array([float(var[sid == i][0]) for i in range(len(classes))])
    q_class = np.array([float(prices[i].mid) for i in reps])

    def receive(members, reference):
        counts = np.array(
            [np.bincount(sid[row], minlength=len(classes)) for row in members]
        )
        held = np.bincount(reference, minlength=len(w))
        quota = int(held.min())
        assert quota > 0
        ref = np.bincount(sid, minlength=len(classes)) * quota
        reference_counts.update(
            batches=len(members),
            physical_members=members.size,
            home_receivers=sector_regional.measure(members, banks),
            balanced_quota_per_word=quota,
            reference_returns=len(reference),
            reference_surplus=int(len(reference) - sum(ref)),
            native_spectrum_counts=counts.sum(0).tolist(),
            reference_spectrum_counts=ref.tolist(),
            modular=empirical(
                "sector modular charge", counts, ref, sid, cat, bank, q_class, qbands
            ),
            variance=empirical(
                "sector energy variance", counts, ref, sid, cat, bank, var_class, vbands
            ),
        )

    lab.run(observe=receive, pulse_seed=638172, register_seed=961332)
    coefficient_band = [
        8 / 4225 * (qbands[2][0] - qbands[3][1]),
        8 / 4225 * (qbands[2][1] - qbands[3][0]),
    ]
    predicted = np.outer(reference_counts["variance"]["bounds"], coefficient_band)
    measured = [float(predicted.min()), float(predicted.max())]
    assert measured[0] <= float(Qfull.a) <= float(Qfull.b) <= measured[1]
    assert (
        reference_counts["modular"]["bounds"][0]
        <= float(direct.a)
        <= float(direct.b)
        <= reference_counts["modular"]["bounds"][1]
    )

    # Construct two positive ensembles with identical heat, work, signed
    # overlap and blank readings by moving along an exact integer null direction.
    sign_changes = []
    for col in range(V.cols):
        u = [Q(s) for s in V[:, col]]
        if sum(u[i] * var[cat == i][0] for i in range(len(u))) == 0:
            continue
        step = np.array(
            [Q(0) if i < 0 else u[i] / int(multiplicity[i]) for i in cat], object
        )
        lo = max(-z[step > 0] / step[step > 0])
        hi = min(-z[step < 0] / step[step < 0])
        for t in (lo / 2, hi / 2):
            trial = z + t * step
            assert min(trial) > 0 and sum(trial) == 1
            f = np.vstack(
                (
                    w[:, :4].T,
                    (w[:, -1] == 0),
                    (w[:, 1] - 1) ** 2,
                    ops.reshape(6, len(w), len(w)).diagonal(axis1=1, axis2=2),
                    fibers.source.tidal.lab.energy.returned(w)[-1],
                )
            )
            assert np.all(f @ (trial - z) == 0) and np.all(trial[cat < 0] == z[cat < 0])
            sign_changes.append(
                dict(
                    scale=str(t),
                    direction=list(map(str, u)),
                    charge=interval(score_bound(trial, native, ids, prices)),
                )
            )
    assert (
        min(r["charge"][1] for r in sign_changes)
        < 0
        < max(r["charge"][0] for r in sign_changes)
    )

    step = 1 / 256
    warmer = thermal.bath(thermal.BETA - step)
    colder = thermal.bath(thermal.BETA + step)
    stiff = thermal.bath(stiffness=1 + step)
    soft = thermal.bath(stiffness=1 - step)
    regions = {}
    for family, ordinary in rho.items():
        panels = {}
        for start in range(4):
            for end in range(start, 4):
                panel, _, _ = thermal.panel(
                    np.array(d, float),
                    ordinary,
                    warmer[family] - colder[family],
                    stiff[family] - soft[family],
                    start,
                    end,
                )
                panels[f"{start}-{end}"] = panel
                if start == end:
                    assert (
                        abs(panel["sharp_modular_nats"] - panel["quantum_modular_nats"])
                        < 1e-12
                    )
                    lo, hi = reference_counts["home_receivers"][str(start)][
                        "receivers"
                    ][family]["inferred_nats"]
                    assert lo <= panel["quantum_modular_nats"] <= hi
        target = float(Qfull.mid) if family == "retained" else 0
        assert abs(panels["0-3"]["quantum_modular_nats"] - target) < 1e-10
        for start in range(4):
            for end in range(start + 1, 4):
                for kind in ("quantum_modular_nats", "sharp_modular_nats"):
                    q = panels[f"{start}-{end}"][kind]
                    q -= (
                        panels[f"{start + 1}-{end}"][kind]
                        + panels[f"{start}-{end - 1}"][kind]
                    )
                    if start + 1 <= end - 1:
                        q += panels[f"{start + 1}-{end - 1}"][kind]
                    panels[f"{start}-{end}"]["connected_" + kind] = q
        regions[family] = panels
    price0 = sum(
        iv.mpf(len(g)) * p / len(w) for g, p in zip(blocks, prices, strict=True)
    )
    upper = max(interval((p - price0) / beta)[1] for p in prices)
    lower = min(interval((p - price0) / beta)[0] for p in prices)
    grid = {}
    for k in (1, 3, 8, 11, 16):
        arm = z if k == 16 else lab.population.histogram(k)[1]["agreed"]
        grid[str(56 * k)] = dict(
            variance=str((arm - native) @ var),
            modular=interval(score_bound(arm, native, ids, prices)),
        )
    return dict(
        protocol_sha256=sha256(
            (HERE / "SECTOR_CALORIMETRY.md").read_bytes()
        ).hexdigest(),
        beta=float(beta.mid),
        full_charge_nats=interval(Qfull),
        full_equivalent_energy=interval(Qfull / beta),
        exchanged_charge_nats="0",
        heat_change="0",
        work_changes=["0"] * 3,
        variance_change=str(delta_var),
        variance_coefficient=interval(coef),
        registers=grid,
        twin_preparations=sign_changes,
        charge_from_measured_variance=measured,
        joint_measurement_coverage_at_least=0.994,
        joint_including_home_receivers_at_least=0.992,
        sector_log_enclosures=[interval(p) for p in prices],
        ordinary_clock_calibration=banks,
        variance_meter=variance,
        independent_source=reference_counts,
        regional=regions,
        thermal_control_step=step,
        matched_compact_energy_envelope=[lower, upper],
        spectral_classes=[
            dict(
                components=[j for j, s in enumerate(shapes) if s == shape],
                centered_characteristic_factors=[
                    dict(polynomial=list(map(str, p)), multiplicity_per_word=str(r))
                    for p, r in shape
                ],
            )
            for shape in classes
        ],
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    path = HERE / "sector_calorimetry.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        from .closed_balls import radar

        radar.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_text()), out
        )
    print(
        out["full_charge_nats"],
        out["independent_source"]["modular"]["bounds"],
        out["charge_from_measured_variance"],
    )
