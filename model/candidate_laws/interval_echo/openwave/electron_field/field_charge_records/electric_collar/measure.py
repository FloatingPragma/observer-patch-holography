"""Read a separate outside material probe after the original terminal capture."""

import argparse
from fractions import Fraction
from hashlib import sha256
from math import factorial, inf, nextafter, sin, sqrt
import json
from pathlib import Path
from runpy import run_path

import numpy as np
from scipy.sparse.linalg import eigsh

from ...fan_path.channel_bound import dyadic_vectors, sqrt_bounds
from ...fan_path.fixture import FanPath, archived_fixture
from ...fan_path.spectral import separated_reference, verified_inertia
from ..asymmetric_clock import gamma
from .apparatus import Collar


HERE = Path(__file__).resolve().parent
RELEASE_TIMES = (27.0625, 30.8125)
TICKS = (0, 8, 16, 32, 64)
STEP = Fraction(1, 32)
DEGREE = 22
PAD = 1e-11


def capture_records(path, bracket):
    start = np.random.default_rng(923).standard_normal(path.h.shape[0])
    values, modes = eigsh(path.h, k=4, which="SA", v0=start, tol=1e-13)
    _, norms, centers, radii, error = dyadic_vectors(
        path.h, values[:2], modes[:, :2], bracket[0]
    )
    frame = modes[:, :2]
    source = path.source()
    np.testing.assert_array_equal(source.T @ source, np.eye(2))
    coefficient = frame.T @ source
    block = path.blocks[path.size - 1, 0]
    # The sharp head word is unchanged by the exact inverse incoming transfer.
    entry, target = path.bridges[-1]
    for row in range(block.rows.start, block.rows.stop):
        port, words = path.basis[row]
        assert port == target
        partner = path.ids[entry, words]
        bridge = path.bridge_term.getrow(row)
        assert bridge.nnz == 2 and bridge[0, row] == 4 and bridge[0, partner] == -4
    projector = path.code @ path.code.T
    np.testing.assert_array_equal(projector @ projector, projector)
    output = []
    round_frame = 8 * sqrt(len(source)) / 2**40
    normalized_trial_error = sum(
        max(abs(1 - a), abs(1 - b))
        for n2 in norms
        for a, b in [sqrt_bounds(Fraction(n2, 2**80))]
    )
    quantization_error = sqrt_bounds(Fraction(len(source), 2**81))[1]
    assert quantization_error + normalized_trial_error < Fraction(round_frame)
    frame_error = float(error) + round_frame
    detuning = max(
        abs(Fraction(float(e)) - c) + r
        for e, c, r in zip(values[:2], centers, radii, strict=True)
    )
    for exposure in RELEASE_TIMES:
        wave = frame @ (np.exp(-1j * exposure * values[:2, None]) * coefficient)
        terminal = wave[block.rows].reshape(-1, 4, 2)
        words = np.einsum("ij,ejk->iek", projector, terminal)
        error_at_capture = (
            (2 + frame_error) * frame_error + exposure * float(detuning) + PAD
        )
        admission = np.sum(abs(coefficient) ** 2, axis=0)
        capture = np.sum(abs(words) ** 2, axis=(0, 1))
        output.append((words, float(error_at_capture), admission, capture))
    return output


def step_bound(h):
    """Taylor remainder and floating-point accumulation with exact dyadic row sums."""
    assert (h != h.T).nnz == 0
    twice = 2 * h
    np.testing.assert_array_equal(twice.data, np.rint(twice.data))
    assert max(abs(twice.data)) < 128
    terms = int(np.diff(h.indptr).max())
    assert 0 < terms <= 64 and h.shape[0] <= 4096
    x = STEP * Fraction(int(np.max(abs(twice).sum(axis=1))), 2)
    weights = [x**k / factorial(k) for k in range(DEGREE + 1)]
    factor = 1 + gamma(2 * terms + 4)
    assert x <= Fraction(5, 2)
    return (
        x ** (DEGREE + 1) / factorial(DEGREE + 1)
        + (1 + gamma(DEGREE + 1)) * sum(w * factor**k for k, w in enumerate(weights))
        - sum(weights)
        + Fraction(1, 2**920)
    )


def kernels(hamiltonian, end):
    n = hamiltonian.shape[0]
    width = 16 * 4
    initial = np.zeros((n, width), dtype=complex)
    initial[:width] = np.eye(width)
    output = [(initial[end].reshape(-1, 16, 4), 0)]
    bound = step_bound(hamiltonian)
    for tick in range(1, TICKS[-1] + 1):
        term = initial
        total = initial.copy()
        for k in range(1, DEGREE + 1):
            term = (-1j * float(STEP) / k) * (hamiltonian @ term)
            total += term
        initial = total
        if tick in TICKS:
            assert tick * bound < Fraction(1, 2)
            error = nextafter(float(2 * width * tick * bound), inf)
            output.append((initial[end].reshape(-1, 16, 4), error))
    return output


def report(path, probe_sign, captures):
    collar = Collar(path, probe_sign)
    probe_degrees = collar.verify_probe(path.board)
    gap = 8 * abs(probe_degrees[0] - probe_degrees[1])
    isolated_frequency = sqrt(4**2 + (gap / 2) ** 2)
    hams, clamp = collar.generators()
    end = np.array([p == collar.probe[1] for _, p, _ in collar.basis])
    histories = {label: kernels(h, end) for label, h in hams.items()}
    cases = []
    maximum = 0.0
    lower_resolved = 0.0
    strongest = None
    for launch, (words, prep_error, admission, capture) in zip(
        RELEASE_TIMES, captures, strict=True
    ):
        stages = []
        for step, tick in enumerate(TICKS):
            tests = {}
            uncertainty = {}
            for label, results in histories.items():
                kernel, local_error = results[step]
                joint = np.einsum("dbw,wei->dbei", kernel, words, optimize=True)
                probability = (abs(joint) ** 2).sum(axis=(0, 2))
                amplitude_error = (
                    prep_error + local_error * sqrt(float(np.max(capture))) + PAD
                )
                allowance = (
                    2 * sqrt(float(np.max(capture))) + amplitude_error
                ) * amplitude_error
                assert probability.min() >= 0 and np.all(
                    probability <= capture + allowance
                )
                uncertainty[label] = nextafter(float(allowance), inf)
                tests[label] = probability
            np.testing.assert_allclose(
                tests["mobile_cut"], tests["neutral_cut"], rtol=0, atol=2e-12
            )
            difference = tests["mobile"] - tests["neutral"]
            at = np.unravel_index(np.argmax(abs(difference)), difference.shape)
            observed = float(abs(difference[at]))
            resolved = observed - uncertainty["mobile"] - uncertainty["neutral"]
            if resolved > lower_resolved:
                maximum = observed
                lower_resolved = resolved
                strongest = {
                    "original_capture_appointment": launch,
                    "probe_exposure": float(tick * STEP),
                    "added_head_word": int(at[0]),
                    "input": int(at[1]),
                    "signed_difference_per_launch": float(difference[at]),
                    "absolute_difference_lower_per_launch": resolved,
                    "signed_difference_conditioned_on_capture": float(
                        difference[at] / capture[at[1]]
                    ),
                }
            stages.append(
                {
                    "probe_exposure": float(tick * STEP),
                    "per_launch_probe_crossings_by_added_word_and_input": {
                        k: v.tolist() for k, v in tests.items()
                    },
                    "calculated_absolute_per_launch_error_upper": uncertainty,
                }
            )
        cases.append(
            {
                "original_capture_appointment": launch,
                "energy_refusal_by_input": (1 - admission).tolist(),
                "band_admitted_but_terminal_charge_code_refused_by_input": (
                    admission - capture
                ).tolist(),
                "terminal_charge_and_code_by_input": capture.tolist(),
                "stages": stages,
            }
        )
    return (
        collar,
        hams,
        {
            "source_sign": path.sign,
            "probe_sign": probe_sign,
            "probe_physical_endpoints": list(collar.probe),
            "probe_port_degrees": probe_degrees,
            "isolated_probe_open_lane_probability_by_exposure": [
                (4 / isolated_frequency) ** 2
                * sin(float(tick * STEP) * isolated_frequency) ** 2
                for tick in TICKS
            ],
            "hub_to_probe_ranges": collar.ranges,
            "exterior_faces_in_order": [list(face) for face in collar.faces],
            "added_seams_bit_order": [list(edge) for edge in collar.added],
            "eligible_face_pairs": collar.eligible,
            "configuration_count": len(collar.basis),
            "parked_source_energy_scalar": clamp,
            "maximum_absolute_probe_response_per_launch": maximum,
            "maximum_resolved_absolute_probe_response_lower_per_launch": lower_resolved,
            "strongest_diagnostic_bank": strongest,
            "launches": cases,
        },
    )


def verify_conjugacy(positive, negative):
    a, pa, _ = positive
    b, pb, _ = negative
    full_word = (1 << len(a.edges)) - 1
    order = [b.ids[x, p, full_word ^ w] for x, p, w in a.basis]
    assert len(set(order)) == len(a.basis) == len(b.basis)
    for label in pa:
        assert (pa[label] != pb[label][order][:, order]).nnz == 0


def signed_response(arms):
    """Evaluate the part of the mobility signal odd in both calibrated signs."""
    selection = {(a["source_sign"], a["probe_sign"]): a for a in arms}
    low = high = None
    for launch in range(len(RELEASE_TIMES)):
        for local in range(len(TICKS)):
            combined = np.zeros((16, 2))
            uncertainty = 0
            for (source, probe), section in selection.items():
                stage = section["launches"][launch]["stages"][local]
                outputs = stage["per_launch_probe_crossings_by_added_word_and_input"]
                errors = stage["calculated_absolute_per_launch_error_upper"]
                conjugate = selection[-source, -probe]["launches"][launch]["stages"][
                    local
                ]
                counterpart = conjugate[
                    "per_launch_probe_crossings_by_added_word_and_input"
                ]
                for law, bank in outputs.items():
                    np.testing.assert_allclose(
                        bank, counterpart[law][::-1], atol=2e-12, rtol=0
                    )
                combined += (
                    source
                    * probe
                    * (np.array(outputs["mobile"]) - outputs["neutral"])
                    / 4
                )
                uncertainty += (errors["mobile"] + errors["neutral"]) / 4
            for operation in (np.argmin, np.argmax):
                b, j = np.unravel_index(operation(combined), combined.shape)
                value = float(combined[b, j])
                result = {
                    "per_launch": value,
                    "error_upper": nextafter(uncertainty, inf),
                    "added_head_word": int(b),
                    "original_input": int(j) + 1,
                    "original_capture_appointment": RELEASE_TIMES[launch],
                    "probe_exposure": float(TICKS[local] * STEP),
                }
                if low is None or value < low["per_launch"]:
                    low = result
                if high is None or value > high["per_launch"]:
                    high = result
    assert low["per_launch"] + low["error_upper"] < 0
    assert high["per_launch"] - high["error_upper"] > 0
    return {"minimum": low, "maximum": high}


def run():
    output = []
    with archived_fixture() as a:
        for size in (3, 4):
            p, m = [FanPath(a, size, q, asymmetric=True) for q in (1, -1)]
            _, bracket, _ = separated_reference(p)
            assert [verified_inertia(p.h, cut)["count_below"] for cut in bracket] == [
                2,
                2,
            ]
            mask = (1 << len(p.field.edges)) - 1
            order = [m.ids[x, tuple(w ^ mask for w in words)] for x, words in p.basis]
            assert (p.h != m.h[order][:, order]).nnz == 0
            results = {}
            for path in (p, m):
                launches = capture_records(path, bracket)
                for probe_sign in (1, -1):
                    results[path.sign, probe_sign] = report(path, probe_sign, launches)
            for probe_sign in (1, -1):
                verify_conjugacy(results[1, probe_sign], results[-1, -probe_sign])
            cases = [case for _, _, case in results.values()]
            output.append(
                {
                    "fan_count": size,
                    "terminal_hub": p.cells[-1][0],
                    "signed_mobility_extremes_per_launch": signed_response(cases),
                    "cases": cases,
                }
            )
            print(f"{size} fans: capture with separate probe resolved", flush=True)
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "outcome_notes": [
            "energy refusal",
            "energy pass but hub charge or code refusal",
            "terminal accepted and probe crosses",
            "terminal accepted and probe stays, given by terminal accepted minus probe crosses",
        ],
        "cases": output,
    }


def main():
    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument("--record", type=Path)
    args = command.parse_args()
    results = run()
    if args.record:
        with args.record.open("x") as stream:
            stream.write(
                json.dumps(results, separators=(",", ":"), allow_nan=False) + "\n"
            )
    else:
        compare = run_path(str(HERE.parents[1] / "replay.py"))["compare"]
        compare(
            results, json.loads((HERE / "results.json").read_text()), "electric_collar"
        )


if __name__ == "__main__":
    main()
