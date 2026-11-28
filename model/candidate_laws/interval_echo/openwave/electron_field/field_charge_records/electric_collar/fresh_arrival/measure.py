"""Read the unchanged independent collar on the fresh local face-0 arrival."""

import argparse
from hashlib import sha256
from itertools import product
import json
from math import pi, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse import diags

from ....fan_path.fixture import FanPath, archived_fixture
from ....fan_path.local_launch.measure import PAULI
from ....fan_path.measure import conjugacy
from ...asymmetric_clock import STEP, UNDERFLOW_PAD, step_error, taylor_step
from ..apparatus import Collar
from ..measure import STEP as COLLAR_STEP
from ..measure import TICKS as COLLAR_TICKS
from ..measure import kernels

HERE = Path(__file__).resolve().parent
FIELD = HERE.parents[2]
SOURCE = FIELD / "fan_path" / "local_launch"
LABEL = [1045, 0, 334, 365]
ARRIVAL = 30.8125
WAIT = 2.0
PAD = 1e-11
LABELS = ("reference_ignored", "reference_y_plus", "reference_y_minus")


def complex_values(matrix):
    return np.stack([matrix.real, matrix.imag], axis=-1).tolist()


def source_receipt():
    saved = json.loads((SOURCE / "word_reads.json").read_text())
    arms = saved["cases"][-1]["stock_arms"]
    flight = next(
        entry
        for entry in next(
            arm["flights"] for arm in arms if arm["preparation"] == "clock_bypassed"
        )
        if entry["material_bridge_rate"] == 4 and entry["appointment"] == ARRIVAL
    )
    pair = next(
        row for row in flight["pair_readings"] if row["port_face_and_words"] == LABEL
    )
    return pair["receiver_reference_pauli_moments"]


def arrivals(path, gauge, paired_words):
    source = path.source() @ gauge
    np.testing.assert_array_equal(source.T @ source, np.eye(2))
    ticks = int(ARRIVAL / STEP)
    assert ticks * STEP == ARRIVAL
    _, _, truncation, roundoff = step_error(path.h)
    per_step_error = truncation + roundoff + UNDERFLOW_PAD
    assert ticks * per_step_error < 1 / 2
    wave = source.astype(complex)
    for _ in range(ticks):
        wave = taylor_step(path.h, wave)
    error = float(np.nextafter(float(4 * ticks * per_step_error), np.inf)) + PAD
    block = path.blocks[path.size - 1, 0]
    local_words = list(block.words[-1])
    indices = [local_words.index(w) for w in paired_words]
    assert paired_words[0] ^ paired_words[1] == path.field.masks[0]
    raw = wave[block.rows].reshape(-1, len(local_words), 2)
    actual = raw[:, indices, :]
    projected_code_then_face = np.einsum("ab,ebi->eai", path.code @ path.code.T, raw)[
        :, indices, :
    ]
    restricted = path.code[indices]
    np.testing.assert_array_equal(restricted.T @ restricted, np.eye(2) / 2)
    prior_port, hub = path.bridges[-1]
    assert hub == LABEL[0]
    # An isolated full bridge SWAP preserves every selected internal word and
    # every upstream label; it changes only the implicit incoming bridge head.
    for row in range(block.rows.start, block.rows.stop):
        port, words = path.basis[row]
        if words[-1] not in paired_words:
            continue
        assert port == hub
        partner = path.ids[prior_port, words]
        term = path.bridge_term.getrow(row)
        assert term.nnz == 2 and term[0, partner] == -4 and term[0, row] == 4
        assert path.bridge_term[partner, partner] == 4
    swap = expm((-1j * pi / 8) * np.array([[4, -4], [-4, 4]]))
    np.testing.assert_allclose(swap, [[0, 1], [1, 0]], atol=2e-15, rtol=0)
    joint = actual.reshape(-1, 4).T @ actual.reshape(-1, 4).conj() / 2
    moments = np.array(
        [[np.trace(joint @ np.kron(a, b)).real for b in PAULI] for a in PAULI]
    )
    vectors = np.array([[1, 1j], [1, -1j]]) / sqrt(2)
    conditioned = np.einsum("ewi,yi->yew", actual, vectors.conj()) / sqrt(2)
    sigmas = np.einsum("yew,yev->ywv", conditioned, conditioned.conj())
    ignored = np.einsum("ewi,evi->wv", actual, actual.conj()) / 2
    np.testing.assert_allclose(sigmas.sum(axis=0), ignored, atol=1e-16, rtol=0)
    density = np.concatenate((ignored[None, :, :], sigmas))
    return {
        "paired_material_words_in_positive_chart_order": paired_words,
        "admitted_flight_frobenius_error_upper": error,
        "field_density_in_real_imaginary_pairs_by_reference_read": complex_values(
            density
        ),
        "joint_receiver_reference_pauli_moments": moments.tolist(),
        "charge_and_face_acceptance_per_original": np.trace(
            density, axis1=1, axis2=2
        ).real.tolist(),
        "direct_face_amplitude_minus_code_then_face_frobenius_norm": float(
            np.linalg.norm(actual - projected_code_then_face) / sqrt(2)
        ),
        "code_then_face_probability_diagnostic_only": float(
            np.vdot(projected_code_then_face, projected_code_then_face).real / 2
        ),
    }, density


def collar_laws(path, signed_probe, paired_words, density, wave_err):
    collar = Collar(path, signed_probe)
    assert collar.hub == LABEL[0] and set(paired_words) <= set(collar.start_words)
    assert collar.verify_probe(path.board) == [5, 5]
    probabilities = np.trace(density, axis1=1, axis2=2).real
    cap_max = max(probabilities)
    capture_error = (2 * sqrt(cap_max) + wave_err) * wave_err + PAD
    hams, scalar = collar.generators()
    held = collar.terms["source_hop"] + collar.terms["source_degree"]
    hams["parked"] = hams["neutral"] + diags(held.diagonal())
    initial_states = [
        collar.ids[collar.hub, collar.probe[0], w | (b << len(collar.original))]
        for b in range(16)
        for w in paired_words
    ]
    columns = np.arange(16 * len(collar.start_words)).reshape(16, -1)
    choices = [collar.start_words.index(word) for word in paired_words]
    np.testing.assert_array_equal(columns[:, choices].ravel(), initial_states)
    end = np.array([p == collar.probe[1] for _, p, _ in collar.basis])
    raw = {}
    errors = {}
    assert COLLAR_TICKS[-1] * COLLAR_STEP == WAIT
    for name, h in hams.items():
        operator, kernel_err = kernels(h, end)[-1]
        selected = operator[:, :, choices]
        effect = np.einsum("dbi,dbj->bij", selected.conj(), selected)
        cells = np.einsum("bvw,ywv->by", effect, density).real
        amplitude_pad = wave_err + kernel_err * sqrt(cap_max) + PAD
        bound = (2 * sqrt(cap_max) + amplitude_pad) * amplitude_pad + PAD
        assert np.min(cells) >= -1e-14 and np.all(cells <= probabilities + bound)
        np.testing.assert_allclose(
            cells[:, 0], cells[:, 1:].sum(axis=1), atol=2e-15, rtol=0
        )
        raw[name] = cells
        errors[name] = float(np.nextafter(bound, np.inf))
    np.testing.assert_allclose(
        raw["mobile_cut"], raw["neutral_cut"], atol=2e-13, rtol=0
    )
    np.testing.assert_allclose(raw["parked"], raw["neutral"], atol=2e-13, rtol=0)
    return (
        collar,
        hams,
        {
            "probe_sign": signed_probe,
            "probe_endpoints": list(collar.probe),
            "probe_faces": [list(face) for face in collar.faces],
            "added_head_bits_in_order": [list(e) for e in collar.added],
            "parked_cost_at_hub": scalar,
            "phase_read_gates_per_original": probabilities.tolist(),
            "clock_refusal_per_original": 0,
            "source_write_refusal_per_original": 0,
            "arrival_refusal_per_original": float(1 - probabilities[0]),
            "probability_error_upper_by_law": errors,
            "gate_probability_error_upper": capture_error,
            "port_crossings_per_original_by_head_law_and_reference": {
                k: v.tolist() for k, v in raw.items()
            },
            "arrived_but_probe_did_not_cross": {
                k: (probabilities - v).tolist() for k, v in raw.items()
            },
        },
    )


def signed_comparison(arms):
    components = np.zeros((16, len(LABELS)))
    conditional = np.zeros_like(components)
    errors = np.zeros(len(LABELS))
    conditional_errors = np.zeros(len(LABELS))
    for (source, probe), (_, _, section) in arms.items():
        k = section["port_crossings_per_original_by_head_law_and_reference"]
        p = np.array(section["phase_read_gates_per_original"])
        q = section["gate_probability_error_upper"]
        n = sum(
            section["probability_error_upper_by_law"][name]
            for name in ("mobile", "neutral")
        )
        d = np.array(k["mobile"]) - np.array(k["neutral"])
        assert np.all(p > q)
        components += source * probe * d / 4
        conditional += source * probe * (d / p) / 4
        errors += n / 4
        conditional_errors += (n + q * np.max(abs(d / p), axis=0)) / (p - q) / 4
    assert np.all(conditional.min(axis=0) + conditional_errors < 0)
    assert np.all(conditional.max(axis=0) - conditional_errors > 0)
    assert np.all(conditional.mean(axis=0) + conditional_errors < 0)
    return {
        "part_odd_in_both_material_signs_per_original_by_head_and_reference": components.tolist(),
        "part_odd_in_both_material_signs_given_gate_by_head_and_reference": conditional.tolist(),
        "absolute_error_upper_per_original_by_reference": errors.tolist(),
        "absolute_error_upper_given_gate_by_reference": conditional_errors.tolist(),
        "uniform_background_mean_per_original_by_reference": components.mean(
            axis=0
        ).tolist(),
        "uniform_background_mean_given_gate_by_reference": conditional.mean(
            axis=0
        ).tolist(),
        "range_given_gate_by_reference": np.stack(
            [conditional.min(axis=0), conditional.max(axis=0)], axis=-1
        ).tolist(),
    }


def run():
    records, arms = {}, {}
    with archived_fixture() as archive:
        plus = FanPath(archive, 4, 1, asymmetric=True)
        minus = FanPath(archive, 4, -1, asymmetric=True)
        gauges = conjugacy(plus, minus)
        mask = (1 << len(plus.field.edges)) - 1
        for path, gauge, words in (
            (plus, np.eye(2), LABEL[-2:]),
            (minus, np.array(gauges[0]), [w ^ mask for w in LABEL[-2:]]),
        ):
            result, density = arrivals(path, gauge, words)
            records[path.sign] = result
            np.testing.assert_allclose(
                result["joint_receiver_reference_pauli_moments"],
                source_receipt(),
                rtol=0,
                atol=3e-12,
            )
            for probe in (-1, 1):
                arms[path.sign, probe] = collar_laws(
                    path,
                    probe,
                    words,
                    density,
                    result["admitted_flight_frobenius_error_upper"],
                )
    for probe in (-1, 1):
        a, m, x = arms[1, probe]
        b, n, y = arms[-1, -probe]
        reverse = (1 << len(a.edges)) - 1
        order = [b.ids[s, p, w ^ reverse] for s, p, w in a.basis]
        assert len(set(order)) == len(a.basis) == len(b.basis)
        for name in m:
            assert (m[name] != n[name][order][:, order]).nnz == 0
            bank = "port_crossings_per_original_by_head_law_and_reference"
            np.testing.assert_allclose(
                x[bank][name], y[bank][name][::-1], atol=4e-13, rtol=0
            )
    return {
        "source_word_receipt_sha256": sha256(
            (SOURCE / "word_reads.json").read_bytes()
        ).hexdigest(),
        "original_local_witness_sha256": sha256(
            (SOURCE / "witness.json").read_bytes()
        ).hexdigest(),
        "collar_recipe_sha256": sha256(
            (HERE.parent / "apparatus.py").read_bytes()
        ).hexdigest(),
        "arrival_label": LABEL,
        "arrival_wait": ARRIVAL,
        "probe_wait": WAIT,
        "reference_read_order": LABELS,
        "word_density_encoding": "each complex entry is [real, imaginary], no gate normalization",
        "original_start_arrival_records": records,
        "arms": [
            {"source_sign": s, **arms[s, p][2]} for s, p in product((-1, 1), repeat=2)
        ],
        "sign_comparison": signed_comparison(arms),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    record = run()
    text = json.dumps(record, indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x", encoding="utf-8") as file:
            file.write(text)
    else:
        from ....replay import compare

        compare(
            json.loads(text),
            json.loads((HERE / "results.json").read_text()),
            "fresh_arrival",
        )
        print("Reproduced the gated signed independent-probe records.")


if __name__ == "__main__":
    main()
