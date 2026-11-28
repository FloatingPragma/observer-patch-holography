"""Invert exterior ranges for four occupied frequencies and their one visible joint moment."""

import math

import networkx as nx
import numpy as np
from sympy import Matrix

from . import sector_passage_charge_lab as charge_lab
from . import sector_passage_lab as lab
from . import sector_passage_source as source
from . import sector_passage_thermal as thermal


def basis(ports, adj, configs, rs, errs):
    # The unit current traverses the unique chain of articulation blocks
    # between two terminals. Off-chain branches carry zero current and the
    # chain's effective resistances add. Each edge belongs to one block.
    graph = nx.from_numpy_array(adj)
    blocks = list(nx.biconnected_components(graph))
    homes = np.searchsorted(ports, [8, 20, 32, 44])
    relevant = []
    for block in blocks:
        present = [j for j, h in enumerate(homes) if h in block]
        if present:
            relevant.append((present, len(block)))
    assert sorted(z[0] for z in relevant) == [[0, 3], [1], [2]]
    cuts = set(nx.articulation_points(graph))
    branch = next(block for block in blocks if homes[1] in block)
    assert ports[list(cuts & branch)].tolist() == [22]
    assert not ({25, 36} & set(ports[list(branch)]))
    occ = np.asarray(configs, dtype=int)
    F = np.column_stack((np.ones(len(occ), int), occ, occ[:, 0] * occ[:, 3]))
    ids = list(Matrix(F.T).rref()[1])
    inv = np.asarray(Matrix(F[ids]).inv(), float)
    coef = np.einsum("ij,jkl->ikl", inv, rs[ids])
    residual = float(np.max(abs(np.einsum("ij,jkl->ikl", F, coef) - rs)))
    err = float(max(errs) * abs(inv).sum(1).max() + 1e-12)
    assert residual < 7 * err
    return (
        F,
        coef,
        err,
        {
            "block_homes": [
                {"homes": ports[homes[z]].tolist(), "size": n} for z, n in relevant
            ],
            "ordinary_calibration_classes": ids,
            "held_classes": len(configs) - len(ids),
            "held_maximum_residual": residual,
            "coefficient_error_upper": err,
        },
    )


def baselines(D):
    # Fix an instrument using ordinary response maps alone. Enroll the
    # strongest mixed baseline; continue with independent normalized rows.
    normalized = D / np.linalg.norm(D, axis=0)
    chosen = [int(np.argmax(abs(D[:, -1])))]
    while len(chosen) < D.shape[1]:
        q, _ = np.linalg.qr(normalized[chosen].T)
        distances = np.linalg.norm(normalized - (normalized @ q) @ q.T, axis=1)
        chosen.append(int(np.argmax(distances)))
    return np.array(chosen)


def exterior_measurement(
    ports, i, j, picks, arms, group, steps, M, inverse, slop, confirmation
):
    measured = []
    n = 2**49 if confirmation else 2**47
    raw = []
    # Independent native excursions on each repeatable held preparation.
    # The exterior reconstruction uses the times at five station pairs and
    # the equal mixing operation. It receives no interior occupancy labels.
    for k, ix in enumerate(picks):
        ends = np.array([i[ix], j[ix]])
        measured_R = []
        bank = []
        for t, idx in enumerate(arms.ravel()):
            g = group[idx]
            counts = lab.acquire(
                steps[g],
                ends[0],
                ends,
                (713819 if confirmation else 933317) + 4 * k + t,
                n,
                excursion=True,
            )
            f = counts["first"][1] / n
            width = math.sqrt(math.log(2 * 20 / 0.004) / (2 * n))
            lo, hi = f - width, f + counts["pending"] / n + width
            if confirmation:
                # Simultaneous empirical Bernstein bounds at the fixed horizon
                # on both the accepted and potentially still pending fraction.
                a, b = f, f + counts["pending"] / n

                def margin(z):
                    ell = math.log(4 * 40 / 0.004)
                    return (
                        math.sqrt(2 * ell * z * (1 - z) / (n - 1))
                        + 7 * ell / (3 * (n - 1))
                        + 2e-9
                    )

                lo, hi = a - margin(a), b + margin(b)
            assert lo > 0
            measured_R.append([1 / (6 * hi), 1 / (6 * lo)])
            bank.append(counts)
        # E[R | material held] is averaged over the equal physical mixtures.
        r = np.array(measured_R)
        measured.append((r[:2].sum(0) - r[2:].sum(0)[::-1]) / 2)
        raw.append({"ports": ports[ends].tolist(), "native_excursions": bank})
    bands = np.array(measured)
    inferred = inverse @ bands.mean(1)
    widths = abs(inverse) @ (np.ptp(bands, axis=1) / 2) + slop
    return {
        "baselines": raw,
        "simultaneous_coverage": 0.996,
        "measured_R_contrast": bands.tolist(),
        "response": M.tolist(),
        "inverted": inferred.tolist(),
        "intervals": np.column_stack((inferred - widths, inferred + widths)).tolist(),
        "inverse_residual_upper": slop.tolist(),
    }


def experiment(
    ports,
    exterior,
    adj,
    configs,
    rs,
    errs,
    steps,
    group,
    weights,
    prepared,
    src_data,
    ordinary,
):
    F, C, err, cert = basis(ports, adj, configs, rs, errs)
    i, j = np.triu_indices(len(exterior), 1)
    i, j = exterior[i], exterior[j]
    D = C[1:, i, j].T
    mean = D.mean(0)
    ordinary_price = float(np.einsum("i,ij->j", ordinary, rs[:, i, j]).mean())
    standard_err = 2 * max(errs) + 1e-12
    assert abs(mean[1] - ordinary_price) < err + standard_err
    picks = baselines(D)
    M = D[picks]
    inverse = np.linalg.inv(M)
    slop = abs(np.eye(5) - inverse @ M).sum(1) + abs(inverse).sum(1) * 5 * err + 1e-10
    assert max(slop) < 1
    words, delta, subset, full, mats = src_data
    arms, ell, reference = source.reference(words, subset, full, mats)
    claim = source.all_rates(
        words, delta, subset, mats, mean, err, ordinary_price, standard_err, ell
    )
    pilot = exterior_measurement(
        ports, i, j, picks, arms, group, steps, M, inverse, slop, False
    )
    measured = exterior_measurement(
        ports, i, j, picks, arms, group, steps, M, inverse, slop, True
    )
    held = charge_lab.read(words, subset, full, mats, reference["log_odds_999"])
    features = F[group[arms[0]]].mean(0) - F[group[arms[1]]].mean(0)
    np.testing.assert_array_equal(features, [0, 0, 0, 0, 0, 0.5])
    rates = {}
    x, y = F[:, 1], F[:, 4]
    for div, (d, eps) in prepared.items():
        v = np.bincount(group, weights=d, minlength=len(F))
        z = v @ F[:, 1:]
        covariance = z[-1] - (weights @ x) * z[3] - (weights @ y) * z[0] - z[0] * z[3]
        rates[str(div)] = {
            "local_feature_contrasts": z.tolist(),
            "covariance_gain": float(covariance),
            "connected_modular_gain": -ell * float(covariance),
            "error_on_covariance": 4 * eps + 1e-11,
            "averaged_R_change": float(z @ mean),
            "averaged_R_connected_piece": float(covariance * mean[-1]),
            "equivalent_fills": float(z @ mean / ordinary_price),
            "connected_equivalent_fills": float(covariance * mean[-1] / ordinary_price),
        }
    joint = D[:, -1]
    return {
        "block_and_calibration": cert,
        "reference": reference,
        "feature_names": [
            "occupied_8",
            "occupied_20",
            "occupied_32",
            "occupied_44",
            "joint_8_44",
        ],
        "mean_coefficients": mean.tolist(),
        "ordinary_mean_response": ordinary_price,
        "correlation_response": {
            "minimum": float(min(joint)),
            "maximum": float(max(joint)),
            "positive": int(sum(joint > err)),
            "negative": int(sum(joint < -err)),
            "unresolved": int(sum(abs(joint) <= err)),
        },
        "measurement": measured,
        "pilot_measurement": pilot,
        "held_source": held,
        "thermal_source": thermal.compare(
            words, prepared, held, claim["minimum_covariance_gain_per_rate"]
        ),
        "rates": rates,
        "all_rates": claim,
    }
