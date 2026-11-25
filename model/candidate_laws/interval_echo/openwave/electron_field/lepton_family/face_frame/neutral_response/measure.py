"""Two gauge-mounted observers test the current response of neutral preparations."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from .. import observers, shared_report
from ..field import transfer

HERE = Path(__file__).resolve().parent
MASKS = [0, 7, 8, 16, 15, 23, 24, 31]
DIM = len(MASKS) ** 2
SHOTS = 2**20
SEED = 473605
TAIL = 0.01 / (2 * 1024)
# These are apparatus settings on the declared effective target, in native
# comparison units. An independent charged standard measures each efficiency.
EXPOSURES = {"alice": {-1: 0.54, 1: 0.72}, "bob": {-1: 1.13, 1: 1.01}}


def ket(left, right):
    return np.eye(DIM)[:, len(MASKS) * MASKS.index(left) + MASKS.index(right)]


def sources():
    # All held states have Q_total=0. Use a common two-slot space, including
    # an empty first slot for the one-body preparations. The instrument gets
    # only the density matrix and its mounting, never this recipe or its ID.
    nu, anti, sterile = ket(0, 15), ket(0, 16), ket(0, 31)
    mimic, conjugate = ket(7, 8), ket(24, 23)
    charged_pair, neutral_pair = ket(8, 23), ket(16, 15)
    specs = [
        [(1.0, nu)],
        [(1.0, anti)],
        [(1.0, sterile)],
        [(1.0, mimic)],
        [(1.0, conjugate)],
        [(0.31, nu), (0.46, anti), (0.23, sterile)],
        [(1.0, charged_pair)],
        [(1.0, (charged_pair + neutral_pair) / np.sqrt(2))],
        [(1.0, (charged_pair - neutral_pair) / np.sqrt(2))],
        [(1.0, ket(7, 24))],
        [(0.4, charged_pair), (0.6, sterile)],
    ]
    return [sum(w * np.outer(v, v.conj()) for w, v in prep) for prep in specs]


def generators():
    occ = np.array([transfer(i, i).diagonal() for i in range(5)])
    q = np.diag((-occ[:3].sum(axis=0) / 3 + occ[3])[MASKS])
    lower = transfer(4, 3)[np.ix_(MASKS, MASKS)]
    eye = np.eye(len(MASKS))

    def lift(a):
        return np.kron(a, eye) + np.kron(eye, a)

    charge = lift(q)
    curr = {-1: lift(lower), 1: lift(lower.T)}
    for sign, generator in curr.items():
        np.testing.assert_allclose(
            charge @ generator - generator @ charge, sign * generator, atol=2e-14
        )
    return charge, curr


def apparatus(frame, charge, current, time):
    # The target spends a reference charge -sign. Its orthogonal read acts
    # on the target factor and commutes with the outgoing charge read.
    single = frame[np.ix_(MASKS, MASKS)]
    np.testing.assert_allclose(single.conj().T @ single, np.eye(len(MASKS)), atol=3e-14)
    G = np.kron(single, single)
    Q = G @ charge @ G.conj().T
    qs, axes = np.linalg.eigh(Q)
    projectors = [
        axes[:, np.isclose(qs, q)] @ axes[:, np.isclose(qs, q)].conj().T
        for q in (-1, 0, 1)
    ]
    projectors.append(np.eye(DIM) - sum(projectors))
    devices = {}
    for sign, generator in current.items():
        action = G @ generator @ G.conj().T
        zero = np.zeros_like(action)
        h = np.block([[zero, action.conj().T], [action, zero]])
        total = np.kron(np.eye(2), Q) + np.kron(np.diag([0, -sign]), np.eye(DIM))
        np.testing.assert_allclose(h @ total, total @ h, atol=1e-13)
        u = expm(-1j * time[sign] * h)[:, :DIM]
        ready, spent = np.split(u, 2)
        effects = [
            k.conj().T @ proj @ k for proj in projectors[:3] for k in (spent, ready)
        ]
        effects.append(sum(k.conj().T @ projectors[3] @ k for k in (spent, ready)))
        np.testing.assert_allclose(sum(effects), np.eye(DIM), atol=2e-13)
        devices[sign] = np.array(effects)
    return G, projectors, devices


def probabilities(effects, rho):
    p = np.einsum("kij,ji->k", effects, rho).real
    assert abs(1 - p.sum()) < 2e-12 and min(p) > -1e-13
    p[p < 5e-15] = 0
    return p / p.sum()


def interval(k):
    return np.array(observers.binomial(int(k), SHOTS, TAIL))


def overlap(a, b):
    return bool(a[0] <= b[1] and b[0] <= a[1])


def rate_calibrations(G, qs, device, rng):
    standards = {}
    for q, state in [(-1, ket(0, 23)), (0, ket(0, 31)), (1, ket(0, 8))]:
        v = G @ state
        counts = rng.multinomial(SHOTS, probabilities(qs, np.outer(v, v.conj())))
        assert counts[q + 1] == SHOTS
        standards[str(q)] = counts.tolist()
    rates = {}
    for sign in (-1, 1):
        # Both independent reference runs start CHARGED. Opposite sign is
        # established with the charge standards; the calibrated hit is spent
        # target plus a newly neutral outgoing body, never an input identity.
        v = G @ ket(0, 8 if sign == -1 else 23)
        counts = rng.multinomial(
            SHOTS, probabilities(device[sign], np.outer(v, v.conj()))
        )
        rates[str(sign)] = dict(
            counts=counts.tolist(),
            charged_standard_input=-sign,
            efficiency99=interval(counts[2]).tolist(),
        )
    return dict(charge_standards=standards, reciprocal_charge_to_neutral=rates)


def forecast(fitted, exposure_a, exposure_b):
    # Count data only. This is the single-doublet / sterile-mixture prediction.
    low, high = interval(fitted)
    weight = [low / exposure_a[1], high / exposure_a[0]]
    if weight[0] > 1:
        return weight, None
    return weight, [weight[0] * exposure_b[0], min(1.0, weight[1]) * exposure_b[1]]


def run():
    rng = np.random.default_rng(SEED)
    pre_rng = np.random.default_rng(SEED + 1)
    charge, current = generators()
    (_, ga), (_, gb) = observers.frames()
    a, qa, da = apparatus(ga, charge, current, EXPOSURES["alice"])
    b, qb, db = apparatus(gb, charge, current, EXPOSURES["bob"])
    wrong = apparatus(ga, charge, current, EXPOSURES["bob"])[2]
    canonical = apparatus(np.eye(len(ga)), charge, current, EXPOSURES["bob"])[2]
    link = b @ a.conj().T
    assert np.linalg.norm(link - np.eye(DIM)) > 1
    cal = {
        "alice": rate_calibrations(a, qa, da, rng),
        "bob": rate_calibrations(b, qb, db, rng),
    }
    # Use the established classical word only once each quantum measurement
    # is over. It carries the final charge digit and an independent spent bit.
    geo, channel = shared_report.channel()
    rows = []
    for j, rho in enumerate(sources()):
        src = a @ rho @ a.conj().T
        dst = link @ src @ link.conj().T
        neutral = [probabilities(q, r)[1] for q, r in [(qa, src), (qb, dst)]]
        np.testing.assert_allclose(neutral, 1, atol=2e-14)
        before = [
            pre_rng.multinomial(SHOTS, probabilities(q, r)).tolist()
            for q, r in ((qa, src), (qb, dst))
        ]
        banks = []
        for sign in (-1, 1):
            idx = 2 * (sign + 1)
            aa = probabilities(da[sign], src)
            bb = probabilities(db[sign], dst)
            np.testing.assert_allclose(
                bb, probabilities(canonical[sign], rho), atol=2e-13
            )
            # Keep Bob's dwell and state, mount the axes in Alice's frame.
            untouched = probabilities(wrong[sign], dst)
            aa_c = rng.multinomial(SHOTS, aa)
            bb_c = rng.multinomial(SHOTS, bb)
            assert bb_c[-1] == 0
            epsa = cal["alice"]["reciprocal_charge_to_neutral"][str(sign)][
                "efficiency99"
            ]
            epsb = cal["bob"]["reciprocal_charge_to_neutral"][str(sign)]["efficiency99"]
            weight, predicted = forecast(aa_c[idx], epsa, epsb)
            observed = interval(bb_c[idx]).tolist()
            joint = np.zeros((7, 7), dtype=np.int64)
            for action, count in zip(channel, bb_c[:6], strict=True):
                joint[action["bob"], action["alice"]] += count
            joint[6, 6] += bb_c[6]
            banks.append(
                dict(
                    sign=sign,
                    alice=aa_c.tolist(),
                    bob=bb_c.tolist(),
                    observed99=observed,
                    alice_hit=float(aa[idx]),
                    bob_hit=float(bb[idx]),
                    alice_doublet_weight99=weight,
                    forecast99=predicted,
                    classification=(
                        "excluded_at_alice"
                        if predicted is None
                        else "compatible"
                        if overlap(predicted, observed)
                        else "held_mismatch"
                    ),
                    same_event_material_joint=joint.tolist(),
                    unchanged_alice_device_tv=float(abs(untouched - bb).sum() / 2),
                )
            )
        rows.append(
            dict(
                blind_source=j,
                neutral_before_both=neutral,
                charge_preflight=before,
                alice_mixture_feasible=bool(
                    sum(bank["alice_doublet_weight99"][0] for bank in banks) <= 1
                ),
                banks=banks,
            )
        )
    return dict(
        shots=SHOTS,
        seed=SEED,
        tail=TAIL,
        source_calibration=cal,
        exposures=EXPOSURES,
        readout_order="outgoing charge -1,0,+1; spent then ready for each; other charge last",
        held_sources=rows,
        report_owners=[geo["bob"], geo["alice"]],
        source_key_after_scoring=[
            "active neutral",
            "conjugate active neutral",
            "sterile singlet",
            "positive active plus negative singlet",
            "negative active plus positive singlet",
            "31% active, 46% conjugate, 23% sterile",
            "oppositely charged active pair",
            "even coherent pair superposition",
            "odd coherent pair superposition",
            "opposite charged singlets",
            "40% opposite charged active pair, 60% sterile singlet",
        ],
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    result = run()
    if ap.parse_args().write:
        (HERE / "records.json").write_text(json.dumps(result, indent=2) + "\n")
    else:
        for row in result["held_sources"]:
            print(
                row["blind_source"],
                [(b["sign"], b["bob_hit"], b["classification"]) for b in row["banks"]],
            )
