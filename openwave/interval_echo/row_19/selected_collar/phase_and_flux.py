"""Compare ordinary energy, entropy and two possible boundary readings at fixed volume."""

from math import log, sqrt

from mpmath import iv
import numpy as np
from scipy.linalg import expm

from . import coherent_balls as coherent

local = coherent.local
N = coherent.N
SEED = 28122924
DT = 1 / 64


def edge_read(m, i, j, rng, half):
    # A quarter relative phase and a balanced native exchange measure the
    # real symmetric and antisymmetric one-load modes, including absence.
    sector = m[np.ix_([i, j], [i, j])]
    unit = np.array([[1, 1j], [1j, 1]]) / sqrt(2) @ np.diag([1, -1j])
    probs = (unit @ sector @ unit.conj().T).diagonal().real
    read = coherent.draws([*probs, 1 - sum(probs)], rng)
    p = [iv.mpf([max(0, k / N - half), min(1, k / N + half)]) for k in read]
    return read, p


def flight(h, m, sites, rng, half):
    inside = np.zeros(len(m), bool)
    inside[sites] = True
    evolution = expm(-0.5j * DT * np.array(h, float))
    sample = []
    prob = []
    for sign in (1, -1):
        phases = np.where(inside, sign * 1j, 1)
        z = (
            evolution
            @ (phases[:, None] * m * phases.conjugate()[None, :])
            @ evolution.conj().T
        )
        chance = float(z.diagonal().real[inside].sum())
        counts = coherent.draws([chance, 1 - chance], rng)
        sample.append(counts)
        prob.append(iv.mpf([counts[0] / N - half, counts[0] / N + half]))
    return sample, (prob[0] - prob[1]) / (2 * DT)


def run():
    iv.dps = 45
    rng = np.random.default_rng(SEED)
    board = local.survey.bulk.federation(0)
    setup, _ = coherent.steady.symmetric_instrument(board)
    _, cells, ports, gates, partner, _, edges, _ = setup
    layouts, h = coherent.stationary_layouts(setup)
    law, _, _ = coherent.steady.parity(board, cells)
    mix = sum(float(w) * layouts[f] for f, w in law.items())
    quiet = layouts[0, 0]
    fresh = quiet + 25 / 36 * (mix - quiet)
    phase = np.ones(len(ports))
    phase[gates[0]] = -1
    ordinary = np.zeros_like(fresh)
    ordinary[partner[0], partner[0]] = 1
    cohorts = {
        "fresh": fresh,
        "shared": quiet + 101 / 144 * (mix - quiet),
        "ordinary": ordinary,
        "phase": phase[:, None] * fresh * phase[None, :],
    }
    row = local.inputs("restored_balls")["surveys"][0]
    q = local.inputs("coherent_balls")["comparisons"][0]
    assert (
        row["center"] == ports[gates[0]]
        and row["direction"] == "forward"
        and row["hop_radii"] == [0, 1]
    )
    balls = row["ball_sites"]
    cuts = [[e for e in edges if (e[0] in B) != (e[1] in B)] for B in balls]
    edge_set = sorted(set(cuts[0] + cuts[1]))
    assert len(cuts[0]) == 5 and all(gates[0] not in e for e in cuts[1])
    wait = float(
        coherent.F(
            local.inputs("stationary_polarization")["mirror_matched"][
                "combined_trace_error_upper"
            ]
        )
    )
    half = sqrt(log(2 * 10000 / 0.01) / (2 * N)) + wait + 5e-13
    archive = {}
    flows = {}
    beats = {}
    finite = {}
    phase_energy = {}
    for name, m in cohorts.items():
        raw = {}
        real = {}
        antisymmetric = {}
        for e in edge_set:
            raw[str(e)], p = edge_read(m, *e, rng, half)
            real[e] = (p[0] - p[1]) / 2
            antisymmetric[e] = p[1]
            if e == tuple(sorted(gates)):
                beats[name] = (p[0], p[1])
        flows[name] = [sum(real[e] for e in cut) for cut in cuts]
        phase_energy[name] = sum(antisymmetric[e] for e in cuts[0])
        shots = []
        finite[name] = []
        for B in balls:
            k, band = flight(h, m, B, rng, half)
            shots.append(k)
            finite[name].append(band)
        archive[name] = {
            "seam_counts": raw,
            "free_flight_counts": shots,
            "release_flux": [coherent.balls.bounds(z) for z in flows[name]],
            "timed_flux": [coherent.balls.bounds(z) for z in finite[name]],
        }
    # The external choice mixes two isospectral joint preparations.
    # Their full relative entropy bounds any measured classical divergence.
    # A pi shift exactly exchanges the two bridge outcomes. The third
    # outcome is fixed, including when the ordinary load is elsewhere.
    a, b = beats["fresh"]
    pmin = (a - b) * iv.log(a / b) / 2
    assert float(pmin.a) > 0
    gain = iv.mpf(q["ordinary_gain"])
    work = phase_energy["phase"] - phase_energy["fresh"]
    dflow = sum(x - y for x, y in zip(flows["phase"], flows["fresh"], strict=True)) / 2
    dtime = (
        sum(x - y for x, y in zip(finite["phase"], finite["fresh"], strict=True)) / 2
    )
    changes = {}
    for name, tag in (("shared", "record_source"), ("ordinary", "ordinary_standard")):
        beta = iv.mpf(row[tag]["public_probability_tangent"])
        delta = sum(x - y for x, y in zip(flows[name], flows["fresh"], strict=True)) / 2
        changes[name] = delta + beta * (flows["fresh"][1] - flows["fresh"][0])
    ds, os = [iv.mpf(q[k + "_entropy_tangent"]) for k in ("record", "ordinary")]
    determinant = changes["shared"] * os - changes["ordinary"] * ds
    return {
        "seed": SEED,
        "starts": N,
        "ports": ports,
        "cuts": [[list(e) for e in c] for c in cuts],
        "time": DT,
        "arms": archive,
        "phase_energy_change": coherent.balls.bounds(work),
        "phase_release_flux_change": coherent.balls.bounds(dflow),
        "phase_finite_flight_flux_change": coherent.balls.bounds(dtime),
        "joined_entropy_derivative_lower": float(pmin.a) - 1e-14,
        "required_area_derivative_upper_at_material_calibration": float(
            (-pmin / gain).b
        )
        + 1e-14,
        "phase_marker_volume_and_area_change": [0, 0],
        "overlap_and_ordinary_flux_tangents": [
            coherent.balls.bounds(changes[k]) for k in ("shared", "ordinary")
        ],
        "flux_common_coefficient_determinant": coherent.balls.bounds(determinant),
        "inputs": coherent.steady.exp.save_inputs(
            "PHASE_AND_FLUX.md",
            "coherent_balls.json",
            "restored_balls.json",
            "stationary_polarization.json",
        ),
    }


if __name__ == "__main__":
    z = run()
    coherent.steady.exp.emit("phase_and_flux", z)
    print({k: v for k, v in z.items() if k not in ("arms", "cuts", "ports", "inputs")})
