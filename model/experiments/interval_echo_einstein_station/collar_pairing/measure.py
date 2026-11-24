"""Regional modular variations and material returns on a formed two-head collar."""

import argparse
import json
from pathlib import Path

import numpy as np

from ..mobile_collar_source import measure as formation
from ..shared_geometry_probe.flight import prepared_states


RECEIPT = Path(__file__).with_name("results.json")
TIMES = (0.5, 1.0, 2.0)
BANK_SIZES = (2, 4, 16, 64)


def log_density(rho):
    weights, axes = np.linalg.eigh(rho)
    if weights[0] <= 0:
        raise ValueError("the regional logarithm requires a faithful state")
    return (axes * np.log(weights)) @ axes.conj().T


def relative_curvature(product, delta):
    """Tr D (d log at product)[D], and the relative upper perturbation."""
    weights, axes = np.linalg.eigh(product)
    if min(weights) <= 0:
        raise ValueError("singular product reference")
    d = axes.conj().T @ delta @ axes
    divided = np.empty((4, 4))
    for i, a in enumerate(weights):
        for j, b in enumerate(weights):
            low, high = sorted((a, b))
            gap = high - low
            divided[i, j] = np.log1p(gap / low) / gap if gap else 1 / low
    relative = d / np.sqrt(weights[:, None] * weights[None, :])
    upper = max(0.0, float(np.linalg.eigvalsh(relative)[-1]))
    return float(np.sum(abs(d) ** 2 * divided)), upper


def charge_slope(product, delta, probability):
    ref = product + probability * delta
    return -probability * np.trace(delta @ (log_density(ref) - log_density(product))).real


def range_reads(product, delta, full_slope, panels):
    rows = []
    for panel in panels:
        enrollments = []
        for size, number in ((1, panel["ball_count"]), (2, panel["doubled_ball_count"])):
            radius = size * panel["relative_ball_radius"]
            hit = number / (12 * panel["cells"])
            slope = charge_slope(product, delta, hit)
            enrollments.append(dict(radius=radius, hit=hit, unflagged_slope_nats=slope,
                                    flagged_slope_nats=hit * full_slope,
                                    unflagged_per_radius_fourth=slope / radius**4))
        inner, outer = enrollments
        rows.append(dict(cells=panel["cells"], enrollments=enrollments,
                         unflagged_doubling_exponent=np.log2(
                             outer["unflagged_slope_nats"] / inner["unflagged_slope_nats"]),
                         flagged_doubling_exponent=panel["doubling_dimension"]))
    return rows


def regional_charge(joint, product, panels):
    delta = joint - product
    curvature, upper = relative_curvature(product, delta)
    rows = []
    for n in (1, *BANK_SIZES):
        p = 1 / n
        slope = charge_slope(product, delta, p)
        rows.append(dict(bank_size=n, strengthening_slope_nats=slope,
                         scaled_magnitude=-slope / p**2,
                         lower=curvature / (1 + upper * p),
                         upper=curvature / (1 - p) if p < 1 else None))
    return dict(curvature_nats=curvature, relative_upper=upper,
                product_min_eigenvalue=float(np.linalg.eigvalsh(product)[0]),
                formed_min_eigenvalue=float(np.linalg.eigvalsh(joint)[0]), samples=rows,
                range_reads=range_reads(product, delta, rows[0]["strengthening_slope_nats"], panels))


def expectation(effect, state):
    return np.trace(effect @ state)


def origin_columns(h, origin):
    launch = np.zeros((h.shape[0], 4), complex)
    launch[4 * origin:4 * origin + 4] = np.eye(4)
    return launch


def material_read(h, ports, preparations):
    p, q = ports
    launched = np.column_stack([origin_columns(h, k) for k in ports])
    h_first = h @ launched[:, :4]
    mean = launched[:, :4].conj().T @ h_first
    waves = formation.flow(h, launched, TIMES)
    output = []
    for sign, pair in preparations.items():
        endpoints = []
        for state in pair[:2]:
            first = expectation(mean, state).real
            centered = h_first - first * launched[:, :4]
            variance = expectation(centered.conj().T @ centered, state).real
            records = []
            for time, wave in zip(TIMES, waves):
                stay = wave[4 * p:4 * p + 4, :4]
                outward = wave[4 * q:4 * q + 4, :4]
                inward = wave[4 * p:4 * p + 4, 4:]
                amplitude = expectation(stay, state)
                complete = inward @ outward
                reached = expectation(outward.conj().T @ outward, state).real
                returned = expectation(complete.conj().T @ complete, state).real
                records.append(dict(time=time, return_amplitude=[amplitude.real, amplitude.imag],
                                    relay_miss=1 - reached, return_miss=reached - returned,
                                    complete=returned))
            endpoints.append(dict(mean_energy=first,
                                  energy_variance=variance,
                                  appointments=records))
        output.append(dict(forming_sign=sign, formed=endpoints[0], independent=endpoints[1]))
    return output


def run():
    survey = RECEIPT.parent.parent / "proper_ball_ruler" / "coarse.json"
    panels = json.loads(survey.read_text())["panels"]
    boards = []
    for level in (0, 1):
        board = formation.federation(level)
        route = formation.region(board)
        states = {}
        for s in (-1, 1):
            joint = prepared_states(board, route, s)[0]
            np.testing.assert_allclose(np.trace(joint), 1, atol=2e-12, rtol=0)
            joint = joint / np.trace(joint).real
            states[s] = (joint, np.kron(*formation.marginals(joint)))
        modular = [dict(forming_sign=s, **regional_charge(*pair[:2], panels))
                   for s, pair in states.items()]
        reads = []
        for q in (-1, 1):
            h = formation.hamiltonian(board, route, q, 8)
            reads.append(dict(probe_sign=q, preparations=material_read(h, route[1:3], states)))
        boards.append(dict(cells=board.cells, ports=board.n, route=list(route),
                           regional=modular, material=reads))
    return {"boards": boards}


def compare(actual, archived):
    if isinstance(actual, dict):
        assert actual.keys() == archived.keys()
        for key, value in actual.items():
            compare(value, archived[key])
    elif isinstance(actual, list):
        assert len(actual) == len(archived)
        for a, b in zip(actual, archived):
            compare(a, b)
    elif actual is None or isinstance(actual, (int, str)):
        assert actual == archived
    else:
        np.testing.assert_allclose(actual, archived, rtol=1e-8, atol=3e-12)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = run()
    if args.write:
        RECEIPT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    else:
        compare(output, json.loads(RECEIPT.read_text()))
        print("Both collar banks and all material appointments agree with the receipt.")
