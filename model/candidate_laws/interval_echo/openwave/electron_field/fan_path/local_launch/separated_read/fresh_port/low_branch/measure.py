"""Verify the low spectrum and the same local charged reference instruments."""

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from math import inf, nextafter, sqrt
from pathlib import Path
from runpy import run_path

import numpy as np
from scipy.sparse.linalg import eigsh

from ......field_charge_records.paths import (
    enclosing_records,
    joint_support,
    local_terms,
)
from .....channel_bound import sqrt_bounds
from .....spectral import verified_inertia
from ..measure import min_checked, operator_basis, partial_transpose
from ...measure import (
    FanPath,
    archived_fixture,
    evolution,
    interval,
    probability,
)

HERE = Path(__file__).resolve().parent
CUT = Fraction(81, 16)
COUNT = 14
BITS = 48
TIMES = (0, 433 / 16, 493 / 16, 497 / 16)
READ_PAD = 2e-12
PORTS = (98, 1010, 1011, 1045)
BANDS = {
    "ground": (0,),
    "lowest_pair": (0, 1),
    "closest_low_pair": (2, 3),
    "adjacent_encoded_pair": (4, 5),
    "five_line_cluster": (2, 3, 4, 5, 6),
    "complete_below_cut": tuple(range(COUNT)),
}


def certified_frame(h):
    assert h.shape == (3520, 3520)
    trial_energies, fitted = eigsh(
        h,
        k=COUNT + 2,
        which="SA",
        v0=np.random.default_rng(98).normal(size=h.shape[0]),
        ncv=60,
        tol=1e-14,
    )
    assert np.all(np.diff(trial_energies) > 0)
    fitted = fitted[:, :COUNT]
    anchors = np.sign(fitted[np.argmax(np.abs(fitted), axis=0), np.arange(COUNT)])
    scale = 2**BITS
    trial = np.rint(fitted * anchors * scale).astype(np.int64)
    integer_h = (2 * h).astype(np.int64)
    assert (integer_h != 2 * h).nnz == 0
    assert np.max(abs(integer_h).sum(axis=1)) * np.max(abs(trial)) < 2**60
    image = integer_h @ trial
    centers = [Fraction(round(float(x) * scale), scale) for x in trial_energies[:COUNT]]
    norm = [sum(int(w) ** 2 for w in line) for line in trial.T]
    radii = []
    for i, mu in enumerate(centers):
        n2 = sum(
            (mu.denominator * int(y) - 2 * mu.numerator * int(z)) ** 2
            for y, z in zip(image[:, i], trial[:, i], strict=True)
        )
        radii.append(
            sqrt_bounds(Fraction(n2, 4 * mu.denominator**2 * norm[i]), 2**62)[1]
        )
    for i in range(COUNT - 1):
        assert centers[i] + radii[i] < centers[i + 1] - radii[i + 1]
    assert centers[-1] + radii[-1] < CUT
    inertia = verified_inertia(h, CUT)
    assert inertia["count_below"] == COUNT
    bounds = []
    records = []
    for i, (mu, r) in enumerate(zip(centers, radii, strict=True)):
        gap = min(
            CUT - mu,
            *(abs(mu - nu) - radii[j] for j, nu in enumerate(centers) if j != i),
        )
        assert gap > r and gap > 0
        err = nextafter(float(2 * r / gap), inf)
        bounds.append(err)
        records.append(
            {
                "ordinal": i,
                "energy": [
                    nextafter(float(mu - r), -inf),
                    nextafter(float(mu + r), inf),
                ],
                "integer_trial_residual_radius_upper": nextafter(float(r), inf),
                "normalized_frame_column_error_upper": err,
                "multiplicity": 1,
            }
        )
    frame = trial / np.sqrt([float(x) for x in norm])
    return frame, np.array(list(map(float, centers))), radii, bounds, records, inertia


def detector(path, port):
    cell, block = next(
        (c, b) for (c, j), b in path.blocks.items() if path.cells[c][j] == port
    )
    pol = sum(path.field.bits.values()) if path.sign == -1 else 0
    words = block.words[cell]
    ordered = [w ^ pol for w in sorted(w ^ pol for w in words)]
    return cell, block, ordered


def tensor_at(vector, detector_chart):
    cell, block, order = detector_chart
    raw = np.moveaxis(
        vector[block.rows].reshape(*block.shape, vector.shape[1]), cell, -2
    )
    return raw.reshape(-1, len(order), vector.shape[1])[
        :, [block.words[cell].index(w) for w in order], :
    ]


def reduced_density(tensor):
    flat = tensor.reshape(len(tensor), -1)
    return flat.T @ flat.conj() / tensor.shape[-1]


def small_read(tensor, error, *, analyzer=False):
    rho = reduced_density(tensor)
    pad = (2 * np.linalg.norm(tensor) + error) * error / tensor.shape[-1] + READ_PAD
    p = float(rho.trace().real)
    out = {
        "joint_arrival_per_input": interval(p, pad),
        "absolute_record_error": float(pad),
    }
    if tensor.shape[-1] == 2:
        out["choi_minimum"] = interval(min_checked(rho), pad)
        partial = partial_transpose(rho)
        out["partial_transpose_minimum"] = interval(min_checked(partial), pad)
        if analyzer:
            bits = len(rho).bit_length() - 1
            assert 2**bits == len(rho)
            axes = tuple(product(range(4), repeat=bits))
            basis = [operator_basis(bits, setting) for setting in axes]
            out["local_Pauli_moments_lex_IXYZ_reference_last"] = [
                float(np.trace(moment @ rho).real) for moment in basis
            ]
            if out["partial_transpose_minimum"][1] < 0:
                _, modes = np.linalg.eigh(partial.real)
                vector = np.rint(2**20 * modes[:, 0]).astype(np.int64)
                if vector[np.argmax(abs(vector))] < 0:
                    vector = -vector
                pure = np.outer(vector, vector) / np.dot(vector, vector)
                measurable = partial_transpose(pure)
                expectation = float(np.trace(measurable @ rho).real)
                if expectation + pad < 0:
                    out["rounded_real_witness"] = {
                        "integer_vector": vector.tolist(),
                        "unnormalized_witness_expectation": interval(expectation, pad),
                        "Pauli_coefficients_same_axis_order": [
                            float(np.trace(measurable @ moment).real / len(rho))
                            for moment in basis
                        ],
                    }
    return out


def low_mode_resolution(columns, error, chart):
    raw = [tensor_at(columns[:, i : i + 1], chart) for i in range(2)]
    rho = [reduced_density(x) for x in raw]
    p = [float(x.trace().real) for x in rho]
    pads = [(2 * np.linalg.norm(x) + error) * error + READ_PAD for x in raw]
    delta = rho[0] / p[0] - rho[1] / p[1]
    d = float(np.linalg.norm(delta, ord="nuc") / 2)
    err = sum(pad / (pi - pad) for pi, pad in zip(p, pads, strict=True))
    lo, hi = interval(d, err)
    return {
        "arrival_probability_from_each_resolved_line": [
            interval(pi, pad) for pi, pad in zip(p, pads, strict=True)
        ],
        "equal_prior_conditional_trace_distance": [max(0, lo), min(1, hi)],
        "equal_prior_conditional_correct_binary_guess_upper": (1 + hi) / 2,
    }


def local_pulses(path):
    ring = path.field.masks[0]
    pol = sum(path.field.bits.values()) if path.sign < 0 else 0
    swap = []
    phase = []
    for i, (port, words) in enumerate(path.basis):
        change = words[2] ^ ring
        j = path.ids[port, (*words[:2], change, words[3])] if port == 1010 else i
        if port == 1010:
            clockwise = path.field.clockwise[0]
            assert words[2] & ring in (clockwise, clockwise ^ ring)
            assert path.parked[i, j] == -0.5
        swap.append(j)
        phase.append(-1 if port == 1010 and words[2] ^ pol == 101 else 1)
    swap = np.array(swap)
    phase = np.array(phase)
    np.testing.assert_array_equal(swap[swap], np.arange(len(swap)))
    return swap, phase


def physical_scope(path, source):
    from model.executable_model import cell
    from model.executable_model.symmetry import rotations

    idx, loads, heads, edges, _ = joint_support(
        path, source, path.sign, local_terms(path)
    )
    np.testing.assert_array_equal(idx, np.arange(len(path.basis)))
    enclosing_records(path, path.sign, loads, heads, edges)
    written = {
        tuple(sorted((path.cells[2][u] % 12, path.cells[2][v] % 12)))
        for u, v in path.field.edges
    }
    assert written <= set(cell().seams)
    gates = (1010 % 12, 1011 % 12)
    fixed = [r for r in rotations() if all(r[x] == x for x in gates)]
    assert fixed == [tuple(range(12))]
    for k, (left, right) in enumerate(path.bridges, 1):
        mask = (path.owner >= k).astype(int)
        i, j = path.h.nonzero()
        crossing = sorted(
            {
                tuple(sorted((path.basis[u][0], path.basis[v][0])))
                for u, v in zip(i, j)
                if mask[u] != mask[v]
            }
        )
        assert crossing == [tuple(sorted((left, right)))]
        assert set(np.ravel(path.h[i[mask[i] != mask[j]], j[mask[i] != mask[j]]])) == {
            -4
        }
    return {
        "proper_cell_poses": len(rotations()),
        "proper_poses_preserving_both_physical_aperture_addresses": len(fixed),
        "gauge_record_and_cut_current_both_checked_on_physical_sector": True,
        "support_size": len(idx),
    }


def pairs(path, wave, energy, radii, error, b, charts):
    gram_error = (2 + error) * error + READ_PAD
    gap = float(energy[1] - energy[0])
    gap_err = float(sum(radii)) + 1e-14
    stock_upper = float(
        np.linalg.eigvalsh(wave[path.owner == 0].T @ wave[path.owner == 0])[-1]
    )
    currents = []
    for i, bridge in enumerate(path.bridges, 1):
        cross = wave[path.owner >= i]
        coherence = float(abs((cross.T @ cross)[0, 1]))
        currents.append(
            {
                "bridge": list(bridge),
                "maximum_cut_current": interval(
                    gap * coherence, gap_err + (gap + gap_err) * gram_error
                ),
                "fixed_input_maximum_population_swing": interval(
                    2 * coherence, 2 * gram_error
                ),
            }
        )
    source_error = error + READ_PAD
    smin = max(0, float(np.linalg.svd(b, compute_uv=False)[-1]) - source_error)
    out = []
    for port, chart in charts.items():
        read = small_read(tensor_at(wave, chart), error, analyzer=port in (1010, 1011))
        pt = read["partial_transpose_minimum"]
        no_recovery = read["choi_minimum"][0] > 0
        out.append(
            {
                "port": port,
                "ordered_full_local_words": chart[2],
                "band_Bell_ideal_detector": read,
                "fine_energy_label_estimation_after_local_arrival": low_mode_resolution(
                    wave, error, chart
                ),
                "full_rank_choi_excludes_exact_heralded_local_reversal": no_recovery,
                "original_source_all_wait_negative_PT_magnitude_lower": (
                    max(0, nextafter(-pt[1] * smin * smin, -inf)) if pt[1] < 0 else 0
                ),
                "two_qubit_PPT_on_original_source_for_all_waits": (
                    len(chart[2]) == 2 and pt[0] > 0 and smin > 0
                ),
            }
        )
    swap, sign = local_pulses(path)
    pulse_waves = (wave[swap], wave * sign[:, None])
    leakage = []
    for name, u in zip(("1010_face_flip", "1010_head_phase"), pulse_waves, strict=True):
        overlap = wave.T @ u
        leakage.append(
            {
                "addressed_detector_operation": name,
                "outside_both_lines_probability_by_resolved_input": [
                    interval(float(1 - x), 3 * gram_error)
                    for x in np.sum(abs(overlap) ** 2, axis=0)
                ],
            }
        )
    return {
        "signed_load_measurement_current_uses_material_number": True,
        "energy_splitting": interval(gap, gap_err),
        "sender_cell_maximum_probability_any_normalized_band_input": interval(
            stock_upper, gram_error
        ),
        "complete_resolved_lines_each_rank_one": True,
        "always_zero_directed_cut_current_for_each_resolved_line": True,
        "maxima_of_native_signed_current_magnitudes": currents,
        "charged_local_detectors": out,
        "local_pulse_energy_refusals": leakage,
    }


def band_reads(path, source, v, energies, radii, errors, chart, raw):
    frame_error = float(np.linalg.norm(errors)) + READ_PAD
    coefficients = v.T @ source
    p = probability(coefficients)
    err_p = (
        2 * np.linalg.norm(coefficients) + frame_error
    ) * frame_error / 2 + READ_PAD
    charts = {port: detector(path, port) for port in PORTS}
    result = {
        "line_ordinals": list(chart),
        "ideal_coherent_global_band_admission_per_original_source": interval(p, err_p),
        "ideal_band_refusal_per_original_source": interval(1 - p, err_p),
        "source_matrix_minimum_singular_value": (
            interval(
                float(np.linalg.svd(coefficients, compute_uv=False)[-1]), frame_error
            )
            if len(energies) > 1
            else [0, 0]
        ),
    }
    if len(energies) == 2:
        result["rank_two_tests"] = pairs(
            path, v, energies, radii, frame_error, coefficients, charts
        )
    records = []
    for eta in TIMES:
        error = sqrt(2) * (
            (2 + frame_error) * frame_error
            + (1 + frame_error) ** 2 * eta * float(max(radii))
        )
        normalized = v @ (np.exp(-1j * eta * energies[:, None]) * coefficients)
        reads = []
        for port, mapping in charts.items():
            rec = small_read(tensor_at(normalized, mapping), error)
            seen = probability(normalized[mapping[1].rows])
            rec["port_refusal_after_band_pass_per_original"] = interval(
                p - seen, err_p + rec["absolute_record_error"]
            )
            if eta in raw:
                # A different instrument: charge read first, followed by a
                # global energy read. The face is not read before this pass.
                w, raw_error = raw[eta]
                accepted = w[mapping[1].rows]
                after = v[mapping[1].rows].T.conj() @ accepted
                ea = raw_error + frame_error * np.linalg.norm(accepted)
                if eta == 0 and port != 98:
                    assert np.count_nonzero(accepted) == 0
                    rec["local_charge_first_then_global_energy_per_original"] = [0, 0]
                else:
                    rec["local_charge_first_then_global_energy_per_original"] = (
                        interval(
                            probability(after),
                            (2 * np.linalg.norm(after) + ea) * ea / 2 + READ_PAD,
                        )
                    )
            reads.append({"port": port, **rec})
        records.append({"wait": eta, "band_first_then_local_charged_read": reads})
    result["original_source_and_diagnostic_reads"] = records
    return result


def case(path, frame, energies, radii, errors):
    source = path.source().astype(complex)
    physical = physical_scope(path, source)
    raw_waves = evolution(path.h, source, (0, 497), sqrt(2), scipy_ticks=(497,))
    raw = {k / 16: x for k, x in raw_waves.items()}
    original = []
    for wait, (wave, error) in raw.items():
        if wait == 0:
            for port in PORTS[1:]:
                assert np.count_nonzero(tensor_at(wave, detector(path, port))) == 0
        original.append(
            {
                "wait": wait,
                "local_charged_read_without_any_energy_preparation": [
                    {
                        "port": port,
                        **small_read(
                            tensor_at(wave, detector(path, port)),
                            error,
                            analyzer=port in (1010, 1011),
                        ),
                        **(
                            {"exact_zero_local_support": True}
                            if wait == 0 and port != 98
                            else {}
                        ),
                    }
                    for port in PORTS
                ],
            }
        )
    return {
        "material_sign": path.sign,
        "physical_and_signed_cut_scope": physical,
        "unfiltered_original": original,
        "diagnostic_bands": {
            label: band_reads(
                path,
                source,
                frame[:, ids],
                energies[list(ids)],
                [radii[i] for i in ids],
                [errors[i] for i in ids],
                ids,
                raw,
            )
            for label, ids in BANDS.items()
        },
    }


def run():
    with archived_fixture() as archive:
        pos, neg = [FanPath(archive, 4, sign, asymmetric=True) for sign in (1, -1)]
        v, energies, radii, errors, spectrum, inertia = certified_frame(pos.h)
        mask = sum(pos.field.bits.values())
        perm = [neg.ids[p, tuple(w ^ mask for w in words)] for p, words in pos.basis]
        assert (pos.h != neg.h[perm][:, perm]).nnz == 0
        np.testing.assert_array_equal(pos.source(), -neg.source()[perm])
        negative = np.empty_like(v)
        negative[perm] = v
        cases = [
            case(pos, v, energies, radii, errors),
            case(neg, negative, energies, radii, errors),
        ]
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "spectral_inertia_count_panel": inertia,
        "complete_spectrum_below_cut": spectrum,
        "exact_charge_conjugation": True,
        "signed_comparisons": cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.loads(json.dumps(run(), allow_nan=False))
    if args.output:
        with args.output.open("x") as out:
            out.write(json.dumps(result, allow_nan=False, indent=2) + "\n")
    else:
        compare = run_path(str(HERE.parents[4] / "replay.py"))["compare"]
        compare(result, json.loads((HERE / "records.json").read_text()), "low_branch")
        print("Replayed the counted low spectrum and signed local charged instruments.")


if __name__ == "__main__":
    main()
