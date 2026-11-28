"""Change the pointer at exactly zero mean energetic cost and test the held geometry."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import common_impulse_pointer as impulse
from . import flight_energy as energy

HERE, lab, reg, iv, d, mat = (
    impulse.HERE,
    impulse.lab,
    impulse.reg,
    impulse.iv,
    impulse.d,
    impulse.mat,
)
SEED = 28112124


def settings(record, fringes):
    phases = {}
    for t, w in fringes["windows"].items():
        subgroup = [z for k, z in record["histories"].items() if k.split("/")[0] == t]

        def quadrature(name):
            return (
                d.DOSE
                * sum(
                    energy.signed(v[name][:2], record["starts"], v["unit"])
                    for g in subgroup
                    for v in g["bonds"]
                )
                / len(subgroup)
            )

        work = quadrature("work")
        current = quadrature("incoming_current")
        length = iv.sqrt(work**2 + current**2)
        assert length.a > 0
        direction = [current / length, work / length]
        charge = impulse.apply_phase(
            w["measured_price_quadratures"]["connected"], direction
        )
        assert 0 not in charge
        orientation = 1 if charge.a > 0 else -1
        direction = [orientation * v for v in direction]
        center = complex(float(current.mid), float(work.mid))
        center *= orientation * (1 - 1e-14) / abs(center)
        assert all(
            v.a <= x <= v.b
            for x, v in zip((center.real, center.imag), direction, strict=True)
        )
        phases[t] = (direction, center, work, current)
    return phases


def stationary_bounds(fringe, frame, phase, full, volume):
    # Apply the ordinary stationary Lagrangian's measured Hessian bound
    # to the complete chord from dense to the calibrated energy-neutral source.
    reference = reg.src.observed(fringe["reference_head_counts"], lab.SHOTS, False)
    std = {
        k: reg.src.observed(v, lab.SHOTS, False)
        for k, v in frame["head_quadrature_counts"].items()
    }
    F = impulse.selected.fringe
    dense = tuple((u + v) / 2 for u, v in zip(std["plus"], std["minus"], strict=True))
    dense_offset = F.trace_bound(F.difference(dense, reference))
    chords = [
        F.trace_bound(F.difference(std[a], std[b])) / 2
        for a, b in (("plus", "minus"), ("left", "right"))
    ]
    length = sum(abs(v).b * m for v, m in zip(phase, chords, strict=True))
    radius = dense_offset + length
    error = (
        fringe["reference_gradient_error_upper"]
        + fringe["mixed_geometry_Hessian_upper"] * radius
    ) * length + 1e-12
    ordinary = energy.flow.dist.read("link_source_contraction")
    from pathlib import Path

    certificate = json.loads(
        Path(mat.c.dose.enclosure.__file__).with_name("certificate.json").read_bytes()
    )
    return (
        iv.mpf(ordinary["ordinary_volume_multiplier"]) * volume
        - full
        + iv.mpf([-error.b, error.b])
    ) / iv.mpf(certificate["gain_interval"])


def run():
    iv.dps = 35
    read = energy.flow.dist.read
    prior, fringe, cost, controls = (
        read(n)
        for n in (
            "charge_selected_flight",
            "source_interference",
            "density_material",
            "ordinary_service_distribution",
        )
    )
    phases = settings(read("flight_energy"), fringe)
    knobs = {k: [z[1].real, z[1].imag] for k, z in phases.items()}
    rng = np.random.default_rng(SEED)
    _, stocks, *_ = mat.c.trials()
    ports, bins = lab.six.survey.return_bins()
    mixed, _, active = d.history(
        stocks["reference"], ports, rng, {"p" + t: z[1] for t, z in phases.items()}
    )
    ref = reg.src.observed(prior["source_reference_heads"], lab.SHOTS, False)
    clock = 3 * np.eye(4) + np.einsum("rtbij,t->rbij", bins, reg.SCORE[:-1] - 3)
    labels = [f"clock/{i}/{j}" for i, j in lab.six.survey.ROUTES]
    anchors = [i for i, r in enumerate(lab.six.survey.ROUTES) if max(r) < 4]
    rows = {}
    for t, (phase, center, work, current) in phases.items():
        frame = prior["windows"][t]
        mod = fringe["windows"][t]
        q = {
            k: impulse.apply_phase(v, phase)
            for k, v in mod["measured_price_quadratures"].items()
        }
        heads = {
            k: reg.src.observed(z, lab.SHOTS, False)
            for k, z in frame["phase_standard_heads"].items()
        }
        baseline = tuple(
            (1 - d.DOSE) * r + d.DOSE * (a + b) / 2
            for r, a, b in zip(ref, heads["plus"], heads["minus"], strict=True)
        )
        base = [impulse.selected.fringe.read(e, baseline, False) for e in clock]
        changes = {
            k: impulse.apply_phase(y, phase)
            for k, y in frame["meter_quadratures"].items()
        }
        geom = impulse.selected.geometry_change(
            base, [changes[k] for k in labels], anchors
        )
        theorem = stationary_bounds(
            fringe, mod, phase, q["full"], iv.mpf(geom["volume"])
        )
        arms = ("p" + t, "dense")
        count = {k: lab.pairs(active[int(t)][k].sum(0), rng, False) for k in arms}
        src = [reg.src.observed(count[k], lab.SHOTS, False) for k in arms]
        observed_q = d.DOSE * d.price_contrast(ref, *src, "connected")
        radar = {
            k: [
                lab.counts(np.r_[v, 1 - v.sum()], rng)
                for v in mat.c.station.pair(bins, mixed[int(t)][k])
            ]
            for k in arms
        }
        means, _ = reg.timing_slopes(radar, lab.SHOTS, ((arms[0], arms[1], 1),))
        direct = impulse.selected.geometry_change(
            means["dense"],
            [a - b for a, b in zip(means[arms[0]], means["dense"], strict=True)],
            anchors,
        )
        mass = {
            k: lab.six.bounds(
                changes[k]
                * iv.mpf(cost["windows"][t]["ordinary_energy"])
                / iv.mpf(controls["windows"][t]["held_ordinary"][k])
            )
            for k in ("shift/-1", "shift/1")
        }
        rows[t] = {
            "phase": list(map(lab.six.bounds, phase)),
            "programmed_phase": knobs[t],
            "quadratures": {
                "work": lab.six.bounds(work),
                "current": lab.six.bounds(current),
            },
            "work_at_center": lab.six.bounds(
                center.real * work - center.imag * current
            ),
            "modular": {k: lab.six.bounds(v) for k, v in q.items()},
            "area_from_stationarity": lab.six.bounds(theorem),
            "source_geometry": geom,
            "ordinary_impulse_equivalents": mass,
            "head_counts": count,
            "clock_counts": radar,
            "direct_charge": lab.six.bounds(observed_q),
            "direct_geometry": direct,
        }
    calorimeter = energy.acquire(knobs, seed=SEED + 1)
    late = energy.decode(calorimeter)
    for t in rows:
        rows[t]["late_work_at_center"] = late["groups"][t]["late_minus_dense"][
            "overlap"
        ]["total"]
    return {
        "seed": SEED,
        "windows": rows,
        "calorimeter": calorimeter,
        "inputs": {
            n: sha256((HERE / n).read_bytes()).hexdigest()
            for n in (
                "ZERO_FLIGHT_WORK.md",
                "flight_energy.json",
                "charge_selected_flight.json",
                "source_interference.json",
                "density_material.json",
                "ordinary_service_distribution.json",
                "link_source_contraction.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    r = run()
    path = HERE / "zero_flight_work.json"
    if args.write:
        path.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), r)
    for t, z in r["windows"].items():
        print(
            t,
            z["direct_charge"],
            z["area_from_stationarity"],
            z["direct_geometry"]["area"],
            z["late_work_at_center"],
        )
