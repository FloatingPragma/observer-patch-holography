"""Read the two original native head coins at a blank load seam and test their sign."""

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.experiments.interval_echo_einstein_station.collar_pairing import (
    measure as source,
)
from . import measure


def sharp_masked(rho):
    probs = rho.diagonal().real
    return np.diag((probs + probs[::-1]) / 2)


def run():
    ex, contrast = candidate.native_law.local_term_parts()
    ids = [candidate.native_law.SRC.index(w) for w in product((0, 1), repeat=2)]
    h = (ex + 8 * contrast)[np.ix_(ids, ids)]
    flip = np.kron(np.eye(2), np.eye(2)[::-1])
    board = source.formation.federation(1)
    route = source.formation.region(board)
    readings = []
    noisy_ports = []
    steps, conductance = measure.green.measured_steps(None)
    visits = measure.green.poisson(measure.central_load(10), conductance)[
        1:-1, 1:-1, 1:-1
    ]
    admission_field = visits / measure.diagonal_visits(10, conductance[0])
    probes = []
    for sign in (-1, 1):
        joint = source.prepared_states(board, route, sign)[0]
        joint /= joint.trace().real
        reference = np.kron(*source.formation.marginals(joint))
        raw = joint.diagonal().real.reshape(2, 2)
        fresh = reference.diagonal().real.reshape(2, 2)
        for left, right in (
            ((0.1, 0.8), (0.1, 0.8)),
            ((0.8, 0.1), (0.8, 0.1)),
            ((0.1, 0.8), (0.8, 0.1)),
            ((0.5, 0.5), (0.1, 0.8)),
        ):
            lp, rp = np.array(left), np.array(right)
            lch, rch = np.array([1 - lp, lp]), np.array([1 - rp, rp])
            change = lch @ (fresh - raw) @ rch.T
            mismatch_change = change[0, 1] + change[1, 0]
            gain = (left[1] - left[0]) * (right[1] - right[0])
            predicted = 17 * gain * (raw[1, 1] - fresh[1, 1])
            np.testing.assert_allclose(8.5 * mismatch_change, predicted, atol=4e-15)
            noisy_ports.append(
                dict(
                    source_sign=sign,
                    left=left,
                    right=right,
                    erased_minus_formed_energy=float(8.5 * mismatch_change),
                )
            )
        for label, orient in (("matching", np.eye(4)), ("opposed", flip)):
            r, p = [
                sharp_masked(orient @ state @ orient.T) for state in (joint, reference)
            ]
            diff = p - r
            energy = float(np.trace(h @ diff))
            unlike = diff.diagonal()[[1, 2]].sum()
            np.testing.assert_allclose(energy, 8.5 * unlike, atol=2e-15)
            np.testing.assert_allclose(h @ diff, diff @ h, atol=2e-15)
            for hit in (0, 0.05, 0.5, 1):
                # Each original head is independently replaced by its own marginal on a miss.
                # All three failed-admission branches have exactly the same product state.
                received = hit**2 * r + (2 * (1 - hit) * hit + (1 - hit) ** 2) * p
                np.testing.assert_allclose(
                    np.trace(h @ (p - received)), hit**2 * energy, atol=3e-15
                )
            readings.append(
                dict(
                    source_sign=sign,
                    orientation=label,
                    formed_sharp_pairs=joint.diagonal().real.tolist(),
                    formed_unlike=float(r[1, 1] + r[2, 2]),
                    retained_energy=float(np.trace(h @ r)),
                    fully_erased_energy=float(np.trace(h @ p)),
                    already_independent_energy=float(
                        np.trace(h @ sharp_masked(orient @ reference @ orient.T))
                    ),
                    independent_unlike=float(p[1, 1] + p[2, 2]),
                    covariance_of_original_load_bits=float(
                        joint[3, 3].real
                        - sum(joint.diagonal().real[2:])
                        * sum(joint.diagonal().real[1::2])
                    ),
                    modular_change_fixed_formed=float(
                        -np.dot(diff.diagonal(), np.log(r.diagonal()))
                    ),
                    modular_change_subtract_product=float(
                        np.dot(diff.diagonal(), np.log(p.diagonal() / r.diagonal()))
                    ),
                    erased_minus_formed_energy=energy,
                    erased_minus_formed_ticket=energy / 512,
                )
            )
            if label == "matching":
                assert energy > 0
                acquire = measure.paired_acquisition(
                    energy / 512 * admission_field**2, steps
                )
                probes.append(dict(source_sign=sign, **acquire))
    thermal = []
    for beta in (0, 1 / 64, 1 / 8, 1):
        state = expm(-beta * h)
        state /= state.trace()
        product_state = np.kron(*source.formation.marginals(state))
        np.testing.assert_allclose(product_state, np.eye(4) / 4, atol=3e-15)
        u, v = np.exp(-8 * beta), np.exp(-9 * beta)
        dephased = sharp_masked(state)
        np.testing.assert_allclose(h @ dephased, dephased @ h)
        gaps = [float(np.trace(h @ (product_state - s))) for s in (state, dephased)]
        z = 2 + u + v
        np.testing.assert_allclose(
            gaps,
            [(15 * (1 - u) + 19 * (1 - v)) / (4 * z), 17 * (2 - u - v) / (4 * z)],
            atol=8e-15,
        )
        thermal.append(
            dict(
                beta=beta, coherent_erasure_energy=gaps[0], sharp_erasure_energy=gaps[1]
            )
        )
    return dict(
        cells=board.cells,
        formation_time=4,
        source=readings,
        noisy_ports=noisy_ports,
        passive_control=thermal,
        native_receiver_acquisitions=probes,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = Path(__file__).with_name("formed_heads.json")
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        saved = json.loads(path.read_text())
        assert (
            saved["cells"] == result["cells"]
            and saved["formation_time"] == result["formation_time"]
        )
        for old, new in zip(
            saved["native_receiver_acquisitions"],
            result["native_receiver_acquisitions"],
            strict=True,
        ):
            for key in (
                "source_sign",
                "proposals_per_address",
                "unmatched_releases",
                "visits",
                "last_drain",
            ):
                assert old[key] == new[key]
            for a, b in zip(old["panels"], new["panels"], strict=True):
                assert a["probe"] == b["probe"] and a["arms"] == b["arms"]
                for key in (
                    "signed_per_expected_excess",
                    "forecast",
                    "simultaneous_99_width",
                ):
                    np.testing.assert_allclose(a[key], b[key], rtol=3e-12)
        for group in ("source", "passive_control", "noisy_ports"):
            for old, new in zip(saved[group], result[group], strict=True):
                for key in new:
                    if isinstance(new[key], str):
                        assert old[key] == new[key]
                    else:
                        np.testing.assert_allclose(
                            old[key], new[key], rtol=1e-10, atol=2e-13
                        )
