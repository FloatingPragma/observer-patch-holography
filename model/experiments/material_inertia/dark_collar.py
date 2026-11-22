"""Read the stored energy, material contact and recoverability of one stationary collar.

Run at repository root with python -m model.experiments.material_inertia.dark_collar.
All waits, launches, preparations and detectors below are supplied apparatus.
"""

import argparse
import json
from math import log, pi
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo.candidate import KAPPA
from model.executable_model import federation
from .federation import propagate
from . import neutral_defect

ARCHIVE = Path(__file__).with_name("dark_collar.json")
WAITS = (0.25, 1, 2)
SHOTS = 131072
SEED = 26112219


def label(magnitude, parity, coin_parity=1):
    """One normalized evacuated collar mode in the native memory basis."""
    result = np.zeros(1 + law.N_LAB)
    for coin, coin_sign in ((0, 1), (1, coin_parity)):
        for d in {magnitude, -magnitude}:
            sign = parity if d < 0 else 1
            result[1 + law.LABELS.index((coin, d))] = coin_sign * sign
    return result / np.linalg.norm(result)


def native_readouts():
    """Verify the lifting of every binary two-port column to the actual 99-state gate."""
    exchange, contrast = law.local_term_parts()
    h = exchange + KAPPA * contrast
    blank = np.eye(law.DIM, law.N_SRC)
    loads = np.eye(law.N_SRC)
    base = np.kron(label(0, 1), loads[:, law.SRC.index((0, 0))])
    unitary = expm(-1j * pi / 128 * h)
    energies = {}
    modes = [(0, 1, 1), (2, 1, 1), (2, -1, 1), (2, 1, -1)]
    for magnitude, parity, coin_parity in modes:
        memory = label(magnitude, parity, coin_parity)
        vacuum = np.kron(memory, loads[:, law.SRC.index((0, 0))])
        gap = float(vacuum @ h @ vacuum)
        for i, (a, b) in enumerate(law.SRC):
            if a > 1 or b > 1:
                continue
            column = np.kron(memory, loads[:, i])
            overlap = coin_parity == 1 and a == b == 1
            np.testing.assert_allclose(
                h @ column, (0 if overlap else gap) * column, rtol=0, atol=3e-13
            )
        if magnitude == 0:
            continue
        start = (base + vacuum) / np.sqrt(2)
        x = abs(np.vdot(start, unitary @ start)) ** 2
        y = abs(np.vdot((base - 1j * vacuum) / np.sqrt(2), unitary @ start)) ** 2
        read = np.arctan2(2 * y - 1, 2 * x - 1) / (pi / 128)
        np.testing.assert_allclose(read, gap, rtol=0, atol=2e-12)
        energies[(parity, coin_parity)] = {
            "gap": float(read),
            "ramsey_plus": float(x),
            "ramsey_quadrature": float(y),
        }
    # Sources have no loads in the measured patch. A charged comparison has
    # one load coherently parked inside or outside it, at fixed total number.
    number = np.tile([a + b for a, b in law.SRC], 1 + law.N_LAB)
    gauge = np.exp(1j * pi * number)
    zero = np.kron(label(2, 1), loads[:, law.SRC.index((0, 0))])
    neutral = abs((1 + np.vdot(zero, gauge * zero)) / 2) ** 2
    inside = blank[:, law.SRC.index((1, 0))]
    parked = blank[:, law.SRC.index((0, 0))]
    # The latter arm has its unit on a distinct outside rail.
    charged = (
        abs((np.vdot(inside, gauge * inside) + np.vdot(parked, gauge * parked)) / 2)
        ** 2
    )
    np.testing.assert_allclose([neutral, charged], [1, 0], rtol=0, atol=1e-13)
    source_gate = expm(-1j * pi / 3 * h)
    error = abs(source_gate[law.SRC.index((0, 1)), law.SRC.index((1, 0))]) ** 2
    np.testing.assert_allclose(error, 1 / 4, rtol=0, atol=2e-14)
    return (
        energies,
        {
            "empty_source_return": float(neutral),
            "one_load_return": float(charged),
            "pulse_angle": pi,
            "observable": "patch load number",
        },
        float(error),
    )


def indirect_reader(energies, exposure_fraction=1):
    """Split and recombine one load next to a parked probe on the marked seam.

    Source (a,b) holds the written memory and a real neighboring blank seam
    (b,c) splits the load from b into spare c. All other terms are parked.
    The source-to-read interval uses the independently interrogated even gap.
    """
    kinetic, contrast = law.local_term_parts()
    h = kinetic + KAPPA * contrast
    t = exposure_fraction * pi / energies[(1, 1)]["gap"]
    source_gate = expm(-1j * t * h)
    rail_indices = [law.SRC.index((1, 0)), law.SRC.index((0, 1))]
    splitter = expm(-1j * pi / 2 * h)[np.ix_(rail_indices, rail_indices)]
    recombiner = expm(-1j * 3 * pi / 2 * h)[np.ix_(rail_indices, rail_indices)]
    rails = np.eye(2)[:, 0]
    entries = {}
    for a in (0, 1):
        records = {}
        for d, parity, coin in ((0, 1, 1), (2, 1, 1), (2, -1, 1), (2, 1, -1)):
            memory = label(d, parity, coin)
            # First arm holds the probe on b, second on c.
            lifts = [
                np.kron(memory, np.eye(law.N_SRC)[:, law.SRC.index((a, b))])
                for b in (1, 0)
            ]
            phases = np.array([np.vdot(v, source_gate @ v) for v in lifts])
            for v, phase in zip(lifts, phases, strict=True):
                np.testing.assert_allclose(
                    source_gate @ v, phase * v, atol=2e-13, rtol=0
                )
            out = recombiner @ (phases * (splitter @ rails))
            gap = 0 if d == 0 else energies[(parity, coin)]["gap"]
            relative = gap * t if a and coin == 1 else 0
            predicted = (1 + np.cos(relative)) / 2
            observed = abs(out[0]) ** 2
            np.testing.assert_allclose(observed, predicted, rtol=0, atol=3e-13)
            records[f"d{d}_sign{parity}_coin{coin}"] = float(observed)
        entries[str(a + 1)] = records
    return {"native_exposure": float(t), "return_by_probe_number": entries}


def common_pair_generator(sector, edge):
    """Replace the blank marked term by the written zero mode, preserving every other seam."""
    a, b = edge
    others = np.array([i for i in range(sector.n) if i not in edge])
    left = np.array([sector.index(a, int(i)) for i in others])
    right = np.array([sector.index(b, int(i)) for i in others])
    rows = np.r_[left, right, left, right]
    cols = np.r_[left, right, right, left]
    weights = np.repeat([KAPPA + 0.5, KAPPA + 0.5, -0.5, -0.5], len(others))
    blank_term = sparse.coo_matrix((weights, (rows, cols)), shape=(sector.P, sector.P))
    return (sector.hamiltonian(KAPPA) - blank_term).tocsr()


def patch_probs(sector, home, edge, base, energy):
    """Read both loads on the marked seam, another in-cell seam, or quiet."""
    selected = np.array(
        [sector.index(a, b) for a, b in sector.car.seams if a // 12 == b // 12 == home]
    )
    mark = sector.index(*edge)
    hit = sparse.coo_matrix(([-energy], ([mark], [mark])), shape=base.shape)
    h = (base + hit).tocsr()
    state = np.zeros(sector.P, complex)
    state[mark] = 1
    readings = []
    previous = 0
    for t in WAITS:
        state = propagate(h, state, t - previous)
        previous = t
        weights = abs(state) ** 2
        np.testing.assert_allclose(weights.sum(), 1, rtol=0, atol=3e-11)
        readings.append(
            [
                float(weights[mark]),
                float(weights[selected].sum() - weights[mark]),
                float(weights.sum() - weights[selected].sum()),
            ]
        )
    return np.array(readings)


def joint_history(reads, error):
    """Fair retained Z, noisy source C, and the native probe outcome B at C."""
    conditional = np.array([[1 - error, error], [error, 1 - error]])
    return conditional[:, :, None] * reads[None, :, :] / 2


def recovered_reference(joint):
    """Conditional recovery from B and signed K(actual) - K(recovered)."""
    bz = joint.sum(axis=1)
    bc = joint.sum(axis=0)
    pb = bc.sum(axis=0)
    if np.any(pb <= 0):
        raise ValueError("The chosen outcome alphabet contains a null bin")
    recovered = bz[:, None, :] * bc[None, :, :] / pb
    active = joint > 0
    if not np.array_equal(active, recovered > 0):
        raise ValueError("The reconstruction changes the support of the experiment")
    remainder = np.zeros_like(joint)
    remainder[active] = np.log(recovered[active] / joint[active])
    return recovered, remainder


def recovery(reads, error):
    """Petz recovery from the patch outcome B on its commuting record support."""
    joint = joint_history(reads, error)
    recovered, remainder = recovered_reference(joint)
    eta = 0.5 * np.abs(joint - recovered).sum()
    delta = float(np.sum((recovered - joint) * remainder))
    entropy = float(-np.sum(joint * remainder))
    capacity = 2 * log((1 - error) / error)
    if delta > capacity * eta + 1e-13 or delta < -1e-13:
        raise ValueError("Variation exceeded the bounded modular remainder")
    np.testing.assert_allclose(
        recovered.sum(axis=(0, 2)), joint.sum(axis=(0, 2)), rtol=0, atol=2e-14
    )
    return {
        "joint": joint.tolist(),
        "recovered": recovered.tolist(),
        "recovery_tv": float(eta),
        "conditional_mutual_information": entropy,
        "remainder_variation": delta,
        "uniform_bound": float(capacity * eta),
        "bare_collar_energy_variation": 0,
    }


def draw_history(reads, error, seed):
    """Generate all enrolled starts using three independent draws per trial."""
    rng = np.random.default_rng(seed)
    z = rng.integers(2, size=SHOTS)
    c = z ^ (rng.random(SHOTS) < error)
    dice = rng.random(SHOTS)
    b = (dice[:, None] > np.cumsum(reads, axis=1)[c, :-1]).sum(axis=1)
    index = (2 * z + c) * reads.shape[1] + b
    counts = np.bincount(index, minlength=4 * reads.shape[1]).reshape(2, 2, -1)
    if counts.sum() != SHOTS:
        raise ValueError("Lost a quiet trial")
    return counts.tolist()


def signed_history_charge(board, energies, full_reader):
    """Read K(actual reference) - K(recovered reference) on fresh preparations."""
    # The source and the probe use (0,1,5); disjoint rails retain Z and supply
    # the error bit. Other seams are parked during these addressed gates.
    source, analyzer, retained, flip = (0, 1), (1, 5), (2, 3), (4, 9)
    if not all(edge in board.seams for edge in (source, analyzer, retained, flip)):
        raise ValueError("A declared register has no native seam")
    exchange, contrast = law.local_term_parts()
    h = exchange + KAPPA * contrast
    left, right = (law.SRC.index(pair) for pair in ((1, 0), (0, 1)))
    rates = (1 / 8, 1 / 4, 3 / 8)
    angles = [2 * np.arcsin(np.sqrt(p)) for p in rates]
    for angle, rate in zip([pi / 2, pi, *angles], [0.5, 1, *rates], strict=True):
        # One sharp port read realizes a Bernoulli draw. A full swap can
        # return the parked error rail to its initial load after it is read.
        column = expm(-1j * angle * h)[:, left]
        np.testing.assert_allclose(
            abs(column[[left, right]]) ** 2, [1 - rate, rate], atol=3e-13, rtol=0
        )

    def gradient(delta, joint):
        # Reference part of d[sum delta log(P_zb P_cb / (P_b P_zcb))].
        def ratio(numerator, denominator):
            return np.divide(
                numerator,
                denominator,
                out=np.zeros_like(numerator),
                where=denominator > 0,
            )

        return (
            ratio(delta.sum(axis=1), joint.sum(axis=1))[:, None, :]
            + ratio(delta.sum(axis=0), joint.sum(axis=0))[None, :, :]
            - ratio(delta.sum(axis=(0, 1)), joint.sum(axis=(0, 1)))
            - ratio(delta, joint)
        )

    def variance(prob, observable, size):
        mean = np.sum(prob * observable)
        return np.sum(prob * (observable - mean) ** 2) / size

    rows = []
    for instrument, reader in enumerate((indirect_reader(energies, 0.5), full_reader)):
        returns = reader["return_by_probe_number"]["2"]
        q = np.array([returns[key] for key in ("d0_sign1_coin1", "d2_sign1_coin1")])
        # Structural zeros follow from the verified 0 and pi arm phase.
        reads = np.array([[1.0, 0.0], [0.5, 0.5] if instrument == 0 else [0.0, 1.0]])
        np.testing.assert_allclose(q, reads[:, 0], atol=3e-13, rtol=0)
        laws = [joint_history(reads, p) for p in rates]
        sigma = laws[1]
        _, exact_charge = recovered_reference(sigma)
        counts = np.array(
            [
                draw_history(reads, p, [SEED, 41, instrument, j])
                for j, p in enumerate(rates)
            ]
        )
        reference_counts = np.array(
            draw_history(reads, rates[1], [SEED, 42, instrument])
        )
        reference = reference_counts / SHOTS
        if not np.array_equal(reference > 0, sigma > 0):
            raise ValueError("The independent reference missed an occupied word")
        _, estimate_charge = recovered_reference(reference)
        readings = []
        for index in (0, 2):
            # Holding P(C,B) exactly fixed includes the source energy, the
            # material probe and every further common processing of that pair.
            np.testing.assert_allclose(
                laws[index].sum(axis=0), sigma.sum(axis=0), rtol=0, atol=1e-14
            )
            comparisons = []
            for collar in (None, 0, 1):
                mask = np.ones_like(sigma)
                if collar is not None:
                    mask[:, 1 - collar, :] = 0
                exact = (laws[index] - sigma) * mask / np.sum(sigma * mask)
                selected = counts[[index, 1]] * mask
                sizes = selected.sum(axis=(1, 2, 3))
                prep, baseline = selected / sizes[:, None, None, None]
                delta = prep - baseline
                observed = np.sum(delta * estimate_charge)
                expected = np.sum(exact * exact_charge)
                error = np.sqrt(
                    variance(reference, gradient(delta, reference), SHOTS)
                    + variance(prep, estimate_charge, sizes[0])
                    + variance(baseline, estimate_charge, sizes[1])
                )
                if abs(observed - expected) > 6 * error + 1e-14:
                    raise ValueError(
                        "Independent modular measurement missed the native prediction"
                    )
                comparisons.append(
                    dict(
                        collar=collar,
                        exact=float(expected),
                        sampled=float(observed),
                        first_order_se=float(error),
                    )
                )
            readings.append(dict(flip=rates[index], variations=comparisons))
        # A diagonal function of C,B gives the same value for both Z's in each
        # occupied bin. The minimax approximation chooses its midrange.
        cb_error = max(
            np.ptp(exact_charge[:, c, b][sigma[:, c, b] > 0]) / 2
            for c in range(2)
            for b in range(2)
            if sigma[:, c, b].sum()
        )
        rows.append(
            dict(
                native_contact_exposure=reader["native_exposure"],
                returns_by_collar=q.tolist(),
                reference_counts=reference_counts.tolist(),
                held_counts=counts.tolist(),
                reference_expectation=float((sigma * exact_charge).sum()),
                best_cb_only_operator_norm_error=float(cb_error),
                readings=readings,
            )
        )
    if not (
        rows[0]["readings"][0]["variations"][0]["exact"] < -0.09
        and rows[0]["readings"][1]["variations"][0]["exact"] > 0.09
        and abs(rows[1]["best_cb_only_operator_norm_error"]) < 1e-14
    ):
        raise ValueError("The native history and full-recovery controls disagree")
    return dict(
        home=0,
        source=source,
        analyzer=analyzer,
        retained_rail=retained,
        flip_rail=flip,
        rates=rates,
        native_flip_angles=angles,
        reference_flip=rates[1],
        history_order=["Z", "C", "B"],
        outcome_order=["return", "other"],
        instruments=rows,
    )


def run():
    energies, number_readout, error = native_readouts()
    board = federation(0)
    sector = PairSector(board)
    # Pick both a gateway-incident source and a held internal source using only incidence.
    sources = []
    for home, gateway in ((0, True), (7, False)):
        edge = next(
            e
            for e in board.seams
            if all(p // 12 == home for p in e)
            and any(board.degree[p] == 6 for p in e) == gateway
        )
        base = common_pair_generator(sector, edge)
        zero = patch_probs(sector, home, edge, base, 0)
        even = patch_probs(sector, home, edge, base, energies[(1, 1)]["gap"])
        odd = patch_probs(sector, home, edge, base, energies[(-1, 1)]["gap"])
        rows = []
        for n, t in enumerate(WAITS):
            reads = np.stack([zero[n], even[n]])
            assessment = recovery(reads, error)
            if np.max(abs(zero[n] - even[n])) < 1e-3:
                raise ValueError(
                    "The live local contact failed to resolve the populated collar"
                )
            rows.append(
                {
                    "exposure": t,
                    "written_zero": zero[n].tolist(),
                    "even_two": even[n].tolist(),
                    "odd_two": odd[n].tolist(),
                    "recovery": assessment,
                    "sampled_z_c_b": draw_history(reads, error, [SEED, home, n]),
                }
            )
        a, b = edge
        spare = next(p for p in board.neighbours[b] if p // 12 == home and p != a)
        sources.append(
            {"home": home, "edge": list(edge), "analyzer_spare": spare, "rows": rows}
        )
    blind = recovery(np.ones((2, 1)), error)
    resolved = recovery(np.eye(2), error)
    reader = indirect_reader(energies)
    np.testing.assert_allclose(
        [
            reader["return_by_probe_number"]["2"][k]
            for k in ("d0_sign1_coin1", "d2_sign1_coin1")
        ],
        [1, 0],
        rtol=0,
        atol=3e-13,
    )
    # The two ideal indirect return outcomes are exactly disjoint by the
    # pi phase difference above. Use this support in the logarithms.
    reader["ideal_indirect_recovery"] = resolved
    reader["sampled_z_c_return_word"] = draw_history(np.eye(2), error, [SEED, 2])
    np.testing.assert_allclose(
        [blind["recovery_tv"], resolved["recovery_tv"]], [1 / 4, 0], rtol=0, atol=1e-13
    )
    return {
        "protocol": {
            "seed": SEED,
            "shots_per_arm_and_time": SHOTS,
            "terminal_bins": ["marked_pair", "other_home_seam", "quiet"],
            "source_flip_rate": error,
            "cells": board.cells,
            "pair_states": sector.P,
        },
        "empty_collar_energy_reads": {
            f"sign_{p}_coin_{c}": v for (p, c), v in energies.items()
        },
        "load_u1_interferometer": number_readout,
        "marked_seam_addresses": {
            str(federation(i).cells): len(federation(i).seams) for i in (0, 1, 2)
        },
        "sources": sources,
        "blind_recovery": blind,
        "full_collar_recovery": resolved,
        "indirect_collar_read": reader,
        "signed_history_charge": signed_history_charge(board, energies, reader),
        "mobile_neutral_pair": neutral_defect.run(),
    }


def compare(old, actual):
    if isinstance(old, dict):
        if old.keys() != actual.keys():
            raise ValueError("Receipt keys changed")
        for key in old:
            compare(old[key], actual[key])
    elif isinstance(old, list):
        for a, b in zip(old, actual, strict=True):
            compare(a, b)
    elif isinstance(old, float):
        np.testing.assert_allclose(actual, old, atol=1e-11, rtol=1e-11)
    elif old != actual:
        raise ValueError("Receipt changed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write", action="store_true", help="regenerate the compact evidence"
    )
    result = run()
    if parser.parse_args().write:
        ARCHIVE.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        compare(json.loads(ARCHIVE.read_text()), result)
    for case in result["sources"]:
        row = case["rows"][-1]
        print(
            case["edge"],
            row["written_zero"][0],
            row["even_two"][0],
            row["recovery"]["recovery_tv"],
            row["recovery"]["remainder_variation"],
        )
