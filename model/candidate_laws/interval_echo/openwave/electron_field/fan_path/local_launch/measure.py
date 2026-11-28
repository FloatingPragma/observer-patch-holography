"""Compare fresh local writes against an untouched reference after free flight."""

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import prod
from pathlib import Path
from runpy import run_path

import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh

from ...field_charge_records.asymmetric_clock import (
    STEP,
    TICKS,
    UNDERFLOW_PAD,
    checked_spectrum,
    interval,
    step_error,
    taylor_step,
)
from ..channel_bound import dyadic_vectors, sqrt_bounds
from ..fixture import FanPath, archived_fixture
from ..measure import conjugacy
from ..spectral import verified_inertia

HERE = Path(__file__).resolve().parent
PAULI = np.array([np.eye(2), [[0, 1], [1, 0]], [[0, -1j], [1j, 0]], np.diag([1, -1])])


def concentration(path, trial, norms, frame_error, cell, detector):
    """Outward bounds on both eigenvalues of the band's position effect."""
    block = path.blocks[cell, 0]
    data = trial[path.owner == cell] if detector == "cell" else trial[block.rows]
    data = data.astype(object)
    divisor = 1
    if detector == "code":
        local = (2 * path.code).astype(int)
        np.testing.assert_array_equal(local, 2 * path.code)
        amplitudes = data.reshape(*block.shape, 2)
        data = np.tensordot(local.T, amplitudes, axes=(1, cell)).reshape(-1, 2)
        divisor = 4
    gram = data.T @ data
    a, b = [Fraction(int(gram[i, i]), divisor * norms[i]) for i in range(2)]
    discriminant = (a - b) ** 2 + 4 * Fraction(
        int(gram[0, 1]) ** 2, divisor**2 * prod(norms)
    )
    gap_lo, gap_hi = sqrt_bounds(discriminant)
    pad = (2 + frame_error) * frame_error
    values = [
        ((a + b - gap_hi) / 2 - pad, (a + b - gap_lo) / 2 + pad),
        ((a + b + gap_lo) / 2 - pad, (a + b + gap_hi) / 2 + pad),
    ]
    return [
        [
            max(0.0, float(np.nextafter(float(lo), -np.inf))),
            min(1.0, float(np.nextafter(float(hi), np.inf))),
        ]
        for lo, hi in values
    ]


def band(path):
    cuts = (Fraction(17, 4), Fraction(145, 32))
    assert all(verified_inertia(path.h, cut)["count_below"] == 2 for cut in cuts)
    rng = np.random.default_rng(923)
    e, v = eigsh(
        path.h, k=4, which="SA", v0=rng.standard_normal(len(path.basis)), tol=1e-13
    )
    trial, norms, _, _, error = dyadic_vectors(
        path.h, e[:2], v[:, :2], cuts[0], bits=50
    )
    # Integer columns, with exact residual/gap enclosures, are used for the
    # numerical launch as well. Each conversion and normalization rounds once.
    frame = trial / np.sqrt(np.asarray(norms, dtype=float))
    ranges = {
        str(path.cells[c][0]): {
            name: concentration(path, trial, norms, error, c, name)
            for name in ("cell", "hub", "code")
        }
        for c in (0, path.size - 1)
    }
    return frame, float(error), ranges


def fresh_write(path, stock):
    """Keep the spectator density matrix after a successful source-code SWAP."""
    decoded = path.decode(stock, 0)
    spectators = decoded.reshape(len(decoded), -1) / np.sqrt(2)
    block = path.blocks[0, 0]
    launch = np.zeros((len(path.basis), 2 * spectators.shape[1]))
    launch[block.rows] = np.einsum("ai,ek->aeki", path.code, spectators).reshape(
        block.rows.stop - block.rows.start, -1
    )
    p = np.vdot(spectators, spectators).real
    inputs = launch.reshape(len(launch), -1, 2)
    np.testing.assert_allclose(
        np.einsum("bki,bkj->ij", inputs, inputs), p * np.eye(2), atol=5e-15
    )
    return launch, float(p)


def pauli_record(kraus, pad):
    """All nine reference/receiver analyzer laws via their Pauli coefficients."""
    flat = kraus.reshape(-1, 4)
    choi = flat.T @ flat.conj() / 2
    coefficients = np.array(
        [[np.trace(choi @ np.kron(a, b)).real for b in PAULI] for a in PAULI]
    )
    trace = coefficients[0, 0]
    partial = choi.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
    lo = float(checked_spectrum(partial)[0])
    # Z dephasing at the fresh reference keeps joint outcomes of a classical
    # source label; the resulting two diagonal reference blocks are positive.
    dephased = np.where(np.tile(np.eye(2), (2, 2)), partial, 0)
    dephased_min = float(checked_spectrum(dephased)[0])
    assert dephased_min >= -pad
    return {
        "reference_joint_success": float(trace),
        "partial_transpose_minimum_interval": interval(lo, pad),
        "conditional_negativity": float(max(0, -lo) / trace) if trace > 0 else None,
        "reference_dephased_minimum": dephased_min,
        "receiver_reference_pauli_moments": coefficients.tolist(),
    }


def terminal(path, wave, written, admitted, wave_error):
    block = path.blocks[path.size - 1, 0]
    # Receiver decoder leaves the four purifying alternatives of the source
    # herald unobserved, independently of which other stations report.
    tensor = path.decode(wave).reshape(*block.shape[:-1], 2, -1, 2)
    tensor = np.moveaxis(tensor, -2, -3)  # spectator words, old labels, out, in
    cell = np.linalg.norm(wave[path.owner == path.size - 1]) ** 2 / 2
    hub = np.linalg.norm(wave[block.rows]) ** 2 / 2
    decoded = np.linalg.norm(tensor) ** 2 / 2
    common = [
        1 - admitted,
        admitted - written,
        written - cell,
        cell - hub,
        hub - decoded,
    ]
    assert min(common) > -1e-11
    np.testing.assert_allclose(sum(common) + decoded, 1, rtol=0, atol=1e-12)
    # The purifying Bell vector has norm sqrt(written); projections contract.
    pad = (2 * np.sqrt(written) + wave_error) * wave_error + 1e-11
    checks = []
    minima = []
    indices = range(path.size - 1)
    for count in range(path.size):
        for subset in combinations(indices, count):
            quiet = tensor
            if subset:
                quiet = tensor.sum(axis=subset) / np.sqrt(
                    prod(block.shape[j] for j in subset)
                )
            record = pauli_record(quiet.reshape(-1, 2, 2), pad)
            selected = record["reference_joint_success"]
            lo, hi = record["partial_transpose_minimum_interval"]
            sites = [path.cells[i][0] for i in subset]
            if hi < 0 and not any(set(m) <= set(sites) for m in minima):
                minima.append(sites)
            verdict = (
                "entangled" if hi < 0 else "separable" if lo >= 0 else "unresolved"
            )
            if cell == 0:
                verdict = "zero support"
            checks.append(
                {
                    "vacuum_observers": sites,
                    "joint_record": verdict,
                    "conditional_success_per_fresh_write": selected / written,
                    "receiver_accepted_requested_vacuum_refused": max(
                        0.0, decoded - selected
                    ),
                    **record,
                }
            )
            np.testing.assert_allclose(
                sum(common) + decoded - selected + selected, 1, atol=1e-12
            )
    return {
        "complete_common_refusals": np.maximum(common, 0).tolist(),
        "signed_charge_in_receiver_cell": float(cell),
        "observable_absolute_error_upper": float(pad),
        "inclusion_minimal_entangled_comparisons": minima,
        "comparisons": checks,
    }


def run_case(path, negative, receive=terminal):
    gauge = conjugacy(path, negative)
    mask = sum(path.field.bits.values())
    for (c, x), block in path.blocks.items():
        opposite = negative.blocks[c, x]
        for left, right in zip(block.words, opposite.words, strict=True):
            assert {mask ^ w for w in left} == set(right)
        if c in (0, path.size - 1) and x == 0:
            local_map = [opposite.words[c].index(w ^ mask) for w in block.words[c]]
            local_gauge = gauge[0 if c == 0 else 1]
            np.testing.assert_array_equal(
                negative.code[local_map], path.code @ local_gauge
            )
    frame, error, ranges = band(path)
    end = path.blocks[path.size - 1, 0]
    for words in end.words[:-1]:
        ring = path.field.ring_laplacian(words)
        np.testing.assert_array_equal(ring.sum(axis=1), 0)
        assert connected_components(ring < 0, directed=False)[0] == 1
    source = path.source()
    eta = float(Fraction(2**41, 38**21))
    # Frobenius error for the purification of the original equal mixture.
    # The dense dot products have at most 3520 terms and unit-norm operands.
    preparation_error = (2 + error) * error + eta + 1e-11
    stocks = []
    launches = []
    for name, stock in (
        ("binary_clock", frame @ (frame.T @ source)),
        ("clock_bypassed", source),
    ):
        admitted = float(np.vdot(stock, stock).real / 2)
        fresh, written = fresh_write(path, stock)
        if name == "clock_bypassed":
            # The original stock is localized and has a product spectator vacuum.
            np.testing.assert_allclose(fresh[:, :2], source / np.sqrt(2), atol=2e-15)
            np.testing.assert_allclose(fresh[:, -2:], source / np.sqrt(2), atol=2e-15)
            np.testing.assert_allclose(fresh[:, 2:-2], 0, atol=2e-15)
            fresh = source.copy()
        assert 0 < written <= admitted + 2e-14
        read = (frame.T @ fresh).reshape(2, -1, 2)
        effect = np.einsum("eki,ekj->ij", read, read) / written
        stock_error = preparation_error if name == "binary_clock" else 1e-15
        # Conditioning here is only for this additional fresh-copy band test.
        # Normalization costs at most 2*stock_error/sqrt(written).
        retest_pad = 4 * stock_error / np.sqrt(written) + 4 * error + 1e-10
        stocks.append(
            {
                "preparation": name,
                "original_clock_admission": admitted,
                "fresh_write_per_original": written,
                "fresh_write_given_clock_acceptance": written / admitted,
                "source_amplitude_error_upper": stock_error,
                "band_retest_after_write_input_extrema_intervals": [
                    interval(value, retest_pad) for value in checked_spectrum(effect)
                ],
                "flights": [],
            }
        )
        launches.append(fresh)

    for rate in (4, 2):
        h = (path.parked + (rate / 4) * path.bridge_term).tocsr()
        _, _, truncation, roundoff = step_error(h)
        step_pad = float(truncation + roundoff + UNDERFLOW_PAD)
        wave = np.concatenate(launches, axis=1).astype(complex)
        batches = np.cumsum([0, *[a.shape[1] for a in launches]])
        for tick in range(TICKS[-1] + 1):
            if tick in (0, *TICKS):
                for i, row in enumerate(stocks):
                    sample = wave[:, batches[i] : batches[i + 1]]
                    written = row["fresh_write_per_original"]
                    np.testing.assert_allclose(
                        np.vdot(sample, sample).real, 2 * written, atol=5e-12, rtol=0
                    )
                    propagation_error = 4 * tick * step_pad * np.sqrt(written)
                    row["flights"].append(
                        {
                            "material_bridge_rate": rate,
                            "appointment": float(tick * STEP),
                            **receive(
                                path,
                                sample,
                                written,
                                row["original_clock_admission"],
                                row["source_amplitude_error_upper"] + propagation_error,
                            ),
                        }
                    )
            if tick < TICKS[-1]:
                wave = taylor_step(h, wave)
    parked = path.parked.tocoo()
    assert np.all(path.owner[parked.row] == path.owner[parked.col])
    return {
        "fan_count": path.size,
        "hub_addresses": [cell[0] for cell in path.cells],
        "sign_conjugation_source_receiver_gauges": gauge,
        "band_frame_error_upper": error,
        "position_effect_eigenvalue_intervals_in_two_line_band": ranges,
        "parked_nonlocal_material_matrix_entries": 0,
        "stock_arms": stocks,
    }


def run():
    with archived_fixture() as archived:
        cases = [
            run_case(FanPath(archived, n, 1, True), FanPath(archived, n, -1, True))
            for n in (3, 4)
        ]
    return {
        "design_sha256": sha256((HERE / "DESIGN.md").read_bytes()).hexdigest(),
        "pauli_order": ["I", "X", "Y", "Z"],
        "common_refusal_order": [
            "stock_clock_refused",
            "local_source_herald_refused",
            "charge_outside_receiving_fan",
            "charge_off_receiving_hub",
            "receiving_ring_code_refused",
        ],
        "cases": cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    record = run()
    if args.output:
        with args.output.open("x") as handle:
            handle.write(json.dumps(record, indent=2, allow_nan=False) + "\n")
    else:
        compare = run_path(str(HERE.parents[1] / "replay.py"))["compare"]
        compare(
            json.loads(json.dumps(record, allow_nan=False)),
            json.loads((HERE / "results.json").read_text()),
            "local_launch",
        )
        print(
            "Reproduced both delayed local-write preparations, flight controls and all station subsets."
        )


if __name__ == "__main__":
    main()
