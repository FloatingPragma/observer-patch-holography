"""Independent arm origins and phases for the emitting and receiving observers."""

import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply

from ..apparatus import native_graph
from .instruments import material
from .measure import SHOTS, sample


def frame(origin, orientation, gradient):
    dest = []
    phases = []
    for mask in range(8):
        marked = [j for j in range(3) if (mask >> j) & 1]
        dest.append(sum(1 << ((origin + orientation * j) % 3) for j in marked))
        phases.append(np.exp(2j * np.pi * gradient * sum(marked) / 3))
    u = np.zeros((8, 8), complex)
    u[dest, np.arange(8)] = phases
    np.testing.assert_allclose(u.conj().T @ u, np.eye(8), atol=2e-14)
    return u


def run():
    h = native_graph()
    # For each physical port the word ordering is (receiver-head B, source-head A).
    frames = (frame(1, 1, 1), frame(1, -1, -1))
    two = np.kron(frames[1], frames[0])
    labels = np.argmax(abs(two), axis=0)
    phases = two[labels, np.arange(64)]
    n = h.shape[0]
    choice = 64 * (np.arange(n) // 64) + labels[np.arange(n) % 64]
    factor = phases[np.arange(n) % 64]
    raw = h.tocoo()
    changed = coo_matrix(
        (
            factor[raw.row] * raw.data * factor[raw.col].conj(),
            (choice[raw.row], choice[raw.col]),
        ),
        shape=h.shape,
    ).tocsr()
    local = (material("hub"), material("tip"))
    inputs = np.zeros((n, 12), complex)
    for endpoint in range(2):
        rows = np.arange(8) if endpoint == 0 else 64 + 8 * np.arange(8)
        inputs[rows, 6 * endpoint : 6 * (endpoint + 1)] = np.column_stack(
            [local[endpoint]["arm"], local[endpoint]["character"]]
        )
    transformed = np.zeros_like(inputs)
    transformed[choice] = factor[:, None] * inputs
    # Re-form the two actual emitters in their own head frames to check
    # the transformed input is not merely a relabelled output archive.
    for j in (0, 1):
        rows = np.arange(8) if j == 0 else 64 + 8 * np.arange(8)
        expected = frames[j] @ np.column_stack([local[j]["arm"], local[j]["character"]])
        np.testing.assert_allclose(
            transformed[rows, j * 6 : (j + 1) * 6], expected, atol=4e-13
        )
    state = expm_multiply(-32j * h, inputs, traceA=-32j * h.diagonal().sum())
    other = expm_multiply(
        -32j * changed, transformed, traceA=-32j * changed.diagonal().sum()
    )
    full_error = float(np.linalg.norm(other[choice] - factor[:, None] * state))
    assert full_error < 2e-11
    streams = [np.random.default_rng(113953) for _ in range(3)]
    totals = [dict(accepted=0, forbidden=0) for _ in range(3)]
    values = []
    maxerr = 0.0
    for endpoint in (0, 1):
        receiver = 1 - endpoint
        for basis_index, basis in enumerate(("arm", "character")):
            geometry = (64 * receiver, 64 * (1 + receiver))
            inputslice = slice(
                endpoint * 6 + basis_index * 3, endpoint * 6 + (basis_index + 1) * 3
            )
            slabs = [
                src[geometry[0] : geometry[1], inputslice].reshape(8, 8, 3)
                for src in (state, other, other)
            ]
            if receiver == 0:
                slabs = [z.transpose(1, 0, 2) for z in slabs]
            source = local[endpoint][basis]
            bases = local[receiver][basis]
            donors = np.column_stack([np.eye(8)[0], local[endpoint]["character"]])
            detectors = (bases, frames[receiver] @ bases, bases)
            meter = (donors, frames[endpoint] @ donors, frames[endpoint] @ donors)
            for i in range(3):
                src = np.linalg.norm(source[:, i]) ** 2
                for k in range(3):
                    for a in range(4) if basis == "character" else [-1]:
                        probabilities = []
                        counts = []
                        for j in range(3):
                            arrival = slabs[j][:, :, i]
                            wave = detectors[j][:, k].conj() @ arrival
                            hub = np.linalg.norm(arrival) ** 2
                            matched = np.linalg.norm(wave) ** 2
                            joint = (
                                abs(meter[j][:, a].conj() @ wave) ** 2
                                if a >= 0
                                else matched
                            )
                            prob = np.array(
                                [
                                    joint,
                                    matched - joint,
                                    hub - matched,
                                    src - hub,
                                    1 - src,
                                ]
                            )
                            bank = sample(prob, streams[j])
                            probabilities.append(prob)
                            counts.append(bank)
                            if basis == "character":
                                totals[j]["accepted"] += bank[0]
                                if (k + max(0, a - 1) - i) % 3:
                                    totals[j]["forbidden"] += bank[0]
                        assert counts[0] == counts[1]
                        maxerr = max(
                            maxerr,
                            float(np.max(abs(probabilities[0] - probabilities[1]))),
                        )
                        values.append([endpoint, basis, i, k, a, counts[0], counts[2]])
    assert totals[0]["forbidden"] == totals[1]["forbidden"] == 0
    assert totals[2]["forbidden"] > 0
    return dict(
        A=dict(origin=1, orientation=1, head_phase_gradient=1),
        B=dict(origin=1, orientation=-1, head_phase_gradient=-1),
        starts=SHOTS,
        seed=113953,
        settings=len(values),
        initial_charge_frames=["A", "B"],
        full_amplitude_error=full_error,
        max_bank_probability_error=maxerr,
        character_counts_shared_transported_untranslated=totals,
        record_columns=[
            "emitter",
            "basis",
            "emission",
            "receiver_projection",
            "donor_projection",
            "shared_and_transported",
            "untranslated_receiver",
        ],
        count_bins=[
            "joint",
            "receiver_only",
            "other_at_hub",
            "no_remote_charge",
            "source_refusal",
        ],
        rows=values,
    )


if __name__ == "__main__":
    result = run()
    vals = result.pop("rows")
    body = (
        json.dumps(result, indent=2)[:-2]
        + ',\n  "rows": [\n    '
        + ",\n    ".join(json.dumps(v) for v in vals)
        + "\n  ]\n}\n"
    )
    Path(__file__).with_name("frames.json").write_text(body)
