"""Recover a three-coordinate charged winding channel by independent terminal records."""

import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from ..apparatus import native_graph
from ..measure import ONE
from ..reciprocal import station_swap
from .instruments import BASES, extract, loan, material
from .spectrum import exact_moments, source_spectra, symmetry

HERE = Path(__file__).resolve().parent
SHOTS = 2**38
SEED = 811363


def interval(k, n=SHOTS):
    # Simultaneous over 512 signed and phase-reading marginals, both tails.
    p = 0.01 / (2 * 512)
    return [
        float(beta.ppf(p, k, n - k + 1)) if k else 0.0,
        float(beta.isf(p, k + 1, n - k)) if k < n else 1.0,
    ]


def pointer():
    result = {}
    for name, basis in BASES.items():
        out = np.zeros((8, 3), complex)
        out[ONE] = basis / np.sqrt(3)
        result[name] = out
    return result


def detection_probabilities(slab, sources, detectors):
    m = sources[ONE]
    after = np.einsum("rsi,ij->rsj", slab, m)
    amps = np.einsum("rk,rsj->ksj", detectors.conj(), after)
    return (
        np.sum(abs(sources) ** 2, axis=0),
        np.sum(abs(after) ** 2, axis=(0, 1)),
        np.sum(abs(amps) ** 2, axis=1),
        amps,
    )


def sample(probs, rng):
    assert min(probs) > -1e-13 and abs(sum(probs) - 1) < 1e-12
    p = np.maximum(probs, 0.0)
    return rng.multinomial(SHOTS, p / p.sum()).tolist()


def witness(rows):
    # Each cell is a physically independent binary accepted-outcome setting.
    # Sum the matched lower bounds and all upper bounds before conditioning.
    numbers = {}
    for b in BASES:
        same = 0
        all_hits = 0
        lower = 0.0
        upper = 0.0
        for basis, j, k, counts in rows:
            if basis != b:
                continue
            n = counts[0]
            all_hits += n
            lo, hi = interval(n)
            upper += hi + 2e-11
            if j == k:
                same += n
                lower += max(0.0, lo - 2e-11)
        numbers[b] = dict(
            match_fraction=same / all_hits,
            simultaneous_lower_99=lower / upper,
            matched=same,
            all_decoded=all_hits,
        )
    lower = sum(z["simultaneous_lower_99"] for z in numbers.values()) - 1
    return dict(measured_bases=numbers, entanglement_fidelity_lower_99=lower)


def crossed(slab, source, detector, rng):
    rows = []
    charge = []
    for basis in BASES:
        src, remote, hit, _ = detection_probabilities(
            slab, source[basis], detector[basis]
        )
        for i in range(3):
            # Independent occupation-only run, with the charge-projector effect
            # and its complementary outcomes. Its settings have no phase analyzer.
            charge.append(
                [basis, i, sample([remote[i], src[i] - remote[i], 1 - src[i]], rng)]
            )
            for k in range(3):
                p = hit[k, i]
                rows.append(
                    [
                        basis,
                        i,
                        k,
                        sample([p, remote[i] - p, src[i] - remote[i], 1 - src[i]], rng),
                    ]
                )
    return dict(bare_charge_counts=charge, counts=rows, **witness(rows))


def joint(slab, sources, detectors, donor, rng):
    src, remote, hit, amps = detection_probabilities(
        slab, sources["character"], detectors["character"]
    )
    local = [np.eye(8)[0], *donor["character"].T]
    rows = []
    allowed = []
    for i in range(3):
        for j, d in enumerate(local):
            for k in range(3):
                a = abs(d.conj() @ amps[k, :, i]) ** 2
                p = hit[k, i]
                counts = sample(
                    [a, p - a, remote[i] - p, src[i] - remote[i], 1 - src[i]], rng
                )
                label = "all_unwound" if j == 0 else j - 1
                match = (k + (0 if j == 0 else j - 1) - i) % 3 == 0
                if not match:
                    assert a < 2e-28
                else:
                    assert a > 1e-11
                allowed.append((match, counts[0], j))
                rows.append([i, label, k, float(a), counts])
    return dict(
        rows=rows,
        original_accepted_allowed=sum(n for m, n, j in allowed if m),
        forbidden_count=sum(n for m, n, j in allowed if not m),
        allowed_net_onehot=sum(1 for m, n, j in allowed if m and n and j),
    )


def fidelity(slab):
    ops = slab[ONE].transpose(1, 0, 2)
    s = np.linalg.norm(ops) ** 2 / 3
    f = sum(abs(np.trace(k)) ** 2 for k in ops) / 9
    values = {}
    for b, B in BASES.items():
        k = B.conj().T @ ops @ B
        values[b] = float(np.sum(abs(np.diagonal(k, axis1=1, axis2=2)) ** 2) / (3 * s))
    # Raw-state MUB upper inequality, without a channel model, and the bound
    # 2/3 on the overlap of Phi_3 with a mixed Schmidt-number-two state.
    phi = np.eye(3).ravel() / np.sqrt(3)
    projectors = []
    for B in BASES.values():
        vectors = np.column_stack([np.kron(v.conj(), v) for v in B.T])
        projectors.append(vectors @ vectors.conj().T)
    P = np.outer(phi, phi)
    np.testing.assert_allclose(projectors[0] @ projectors[1], P, atol=2e-14)
    assert np.linalg.eigvalsh(np.eye(9) + P - sum(projectors))[0] > -2e-14
    return dict(
        conditional_fidelity=float(f / s),
        original_instrument_yield=float(s / 9),
        complementary_matches=values,
        complementary_lower=sum(values.values()) - 1,
    )


def run():
    h = native_graph()
    group, perms = symmetry(h)
    station_swap(h)
    endpoints = [material("hub"), material("tip")]
    reference = pointer()
    y = np.zeros((h.shape[0], 6), complex)
    y[ONE, :3] = np.eye(3)
    y[64 + 8 * np.array(ONE), 3:] = np.eye(3)
    y = expm_multiply(-32j * h, y, traceA=-32j * h.diagonal().sum())
    left = y[64:128, :3].reshape(8, 8, 3)
    right = y[:64, 3:].reshape(8, 8, 3).transpose(1, 0, 2)
    discrepancy = np.linalg.norm(left - right)
    assert discrepancy < 5e-14
    rng = np.random.default_rng(SEED)
    benches = []
    splits = []
    for direction, slab_s in enumerate([left, right]):
        choices = [endpoints[direction], reference]
        other = [endpoints[1 - direction], reference]
        for a, b in ((0, 0), (1, 0), (0, 1)):
            bench = crossed(slab_s, choices[a], other[b], rng)
            assert bench["entanglement_fidelity_lower_99"] > 2 / 3
            benches.append(dict(direction=direction, source=a, receiver=b, **bench))
        splits.append(joint(slab_s, choices[0], other[0], choices[0], rng))
    # Full dephasing of emitted one-hot labels predicts a separable Choi
    # state; hold the subsequent native dynamics and receivers fixed.
    ops = left[ONE].transpose(1, 0, 2)
    deco = np.zeros((8, 24, 3), complex)
    for i in range(3):
        deco[:, 8 * i : 8 * (i + 1), i] = left[:, :, i]
    null = crossed(deco, endpoints[0], endpoints[1], rng)
    assert null["entanglement_fidelity_lower_99"] < 1 / 3
    physical = []
    for side, slab in zip(("tip", "hub"), (left, right), strict=True):
        matrix, ports, refused = extract(side)
        # Bare-hub check first. Before extraction all source heads are traced.
        # Removing the three calibrated local output phases defines the frame.
        out = np.einsum("ar,rsi->sai", matrix, slab)
        phases = (np.diag(matrix[:, ONE]) * np.sqrt(3)).conj()
        canonical = phases[None, :, None] * out * np.sqrt(3)
        np.testing.assert_allclose(canonical, slab[ONE].transpose(1, 0, 2), atol=5e-14)
        physical.append(
            dict(
                side=side,
                coherent_storage_ports=ports,
                rejected_ports=refused,
                relative_phases=[[float(z.real), float(z.imag)] for z in phases],
                extraction_probability_given_onehot=1 / 3,
            )
        )
    dephased = []
    for i in range(3):
        separated = ops.copy()
        separated[:, :, np.arange(3) != i] = 0
        dephased.extend(separated)
    dephased = np.asarray(dephased)
    s = np.linalg.norm(dephased) ** 2 / 3
    dephased_matches = [
        float(
            np.sum(abs(np.diagonal(B.conj().T @ dephased @ B, axis1=1, axis2=2)) ** 2)
            / (3 * s)
        )
        for B in BASES.values()
    ]
    return dict(
        symmetry=group,
        source_moment=exact_moments(h),
        visible_spectra=source_spectra(h, perms),
        probe_loan=loan(),
        time=32,
        starts_per_setting=SHOTS,
        seed=SEED,
        apparatus=["charge-path pulses", "independent three-path pointer"],
        direction_code=["4 to 160", "160 to 4"],
        crossed_count_order=[
            "phase hit",
            "other at remote hub",
            "no remote hole",
            "source refusal",
        ],
        bare_count_order=["remote hole", "not remote", "source refusal"],
        joint_count_order=[
            "two phase hits",
            "receiver only",
            "other at remote hub",
            "no remote hole",
            "source refusal",
        ],
        mirror_amplitude_error=float(discrepancy),
        exact=fidelity(left),
        dephased_matches=dephased_matches,
        counted_dephasing_control=null,
        common_coherent_extractors=physical,
        observed=benches,
        joint_character_read=splits,
    )


if __name__ == "__main__":
    # Keep exact individual bank counts on one line each in the receipt.
    data = run()
    joint_data = data.pop("joint_character_read")
    observed = data.pop("observed")
    body = json.dumps(data, indent=2)[:-2]
    for key, rows in [("observed", observed), ("joint_character_read", joint_data)]:
        body += (
            ",\n  "
            + json.dumps(key)
            + ": [\n    "
            + ",\n    ".join(json.dumps(row) for row in rows)
            + "\n  ]"
        )
    (HERE / "results.json").write_text(body + "\n}\n")
