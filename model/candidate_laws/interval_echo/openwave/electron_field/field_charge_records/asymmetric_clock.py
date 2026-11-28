"""Independently enclose the asymmetric-guide receiver laws at two appointments."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial, inf, nextafter
from pathlib import Path
from runpy import run_path

import numpy as np
from scipy.sparse.csgraph import connected_components

from ..fan_path.fixture import FanPath, archived_fixture
from .measure import HERE


TICKS = (433, 493)
STEP = Fraction(1, 16)
DEGREE = 26
UNIT_ROUNDOFF = Fraction(1, 2**53)
UNDERFLOW_PAD = Fraction(1, 2**920)
READ_PAD = 1e-11


def gamma(operations):
    total = operations * UNIT_ROUNDOFF
    return total / (1 - total)


def step_error(h):
    """Outward Taylor error relative to the Frobenius norm of the input columns."""
    assert h.has_canonical_format and (h != h.T).nnz == 0
    twice = 2 * h
    np.testing.assert_array_equal(twice.data, np.rint(twice.data))
    assert max(abs(twice.data)) < 128
    width = int(np.max(np.diff(h.indptr)))
    assert width <= 64 and h.shape[0] <= 4096
    norm = Fraction(int(np.max(abs(twice).sum(axis=1))), 2)
    x = norm * STEP
    assert x <= Fraction(5, 2)
    # A real sparse row applied to a complex column has two independent
    # real dot products. The pad covers coefficient division and multiplication.
    term_factor = 1 + gamma(2 * width + 4)
    sum_factor = 1 + gamma(DEGREE + 1)
    weights = [x**k / factorial(k) for k in range(DEGREE + 1)]
    roundoff = sum_factor * sum(
        w * term_factor**k for k, w in enumerate(weights)
    ) - sum(weights)
    truncation = x ** (DEGREE + 1) / factorial(DEGREE + 1)
    return norm, width, truncation, roundoff


def taylor_step(h, source):
    term = source
    total = source.copy()
    for k in range(1, DEGREE + 1):
        term = (-1j * float(STEP) / k) * (h @ term)
        total += term
    return total


def source_columns(path):
    raw = path.source()
    block = path.blocks[0, 0]
    assert all(size == 4 for size in block.shape)
    assert np.count_nonzero(raw) == 2 * np.prod(block.shape)
    np.testing.assert_array_equal(abs(raw[block.rows]), 2.0**-path.size)
    np.testing.assert_array_equal(raw.T @ raw, np.eye(2))
    return raw.astype(complex)


def receiver(path, wave):
    block = path.blocks[path.size - 1, 0]
    np.testing.assert_array_equal(abs(path.code), 0.5)
    np.testing.assert_array_equal(path.code.T @ path.code, np.eye(2))
    for words in block.words[:-1]:
        ring = path.field.ring_laplacian(words)
        np.testing.assert_array_equal(ring.sum(axis=1), 0)
        assert connected_components(ring < 0, directed=False)[0] == 1
    raw = wave[block.rows].reshape(-1, 4, 2)
    kraus = np.einsum("do,edi->eoi", path.code, raw, optimize=False)
    assert len(kraus) <= 64
    vacuum = kraus.sum(axis=0, keepdims=True) / np.sqrt(len(kraus))
    return kraus, vacuum


def checked_spectrum(matrix):
    """Check a small dense congruence before using its approximate eigenvalues."""
    hermitian = (matrix + matrix.conj().T) / 2
    assert hermitian.shape[0] <= 4 and np.linalg.norm(hermitian, ord="fro") <= 2
    values, modes = np.linalg.eigh(hermitian)
    gram_defect = modes.conj().T @ modes - np.eye(len(values))
    difference = (modes * values) @ modes.conj().T - hermitian
    assert np.linalg.norm(gram_defect, ord="fro") < 1e-13
    assert np.linalg.norm(difference, ord="fro") < 1e-13
    return values


def interval(value, error):
    return [
        nextafter(float(value) - error, -inf),
        nextafter(float(value) + error, inf),
    ]


def read_bounds(kraus, wave_error):
    flat = kraus.reshape(-1, 4)
    assert np.linalg.norm(flat, ord="fro") < 1.5
    choi = flat.T @ flat.conj() / 2
    partial = choi.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
    effect = np.einsum("eoi,eoj->ij", kraus.conj(), kraus, optimize=False)
    minimum = checked_spectrum(partial)[0]
    low, high = checked_spectrum(effect)
    probability_pad = 3 * wave_error + READ_PAD
    choi_pad = 2 * wave_error + READ_PAD
    return {
        "success_by_input": [
            interval(x, probability_pad) for x in effect.diagonal().real
        ],
        "worst_input_success": interval(low, probability_pad),
        "best_input_success": interval(high, probability_pad),
        "reference_joint_success": interval(choi.trace().real, probability_pad),
        "reference_choi_partial_transpose_minimum": interval(minimum, choi_pad),
    }


def run_case(path):
    norm, width, truncation, roundoff = step_error(path.h)
    source = source_columns(path)
    saved = []
    for tick in range(1, TICKS[-1] + 1):
        source = taylor_step(path.h, source)
        if tick not in TICKS:
            continue
        local_error = truncation + roundoff + UNDERFLOW_PAD
        assert tick * local_error < Fraction(1, 2)
        wave_error = nextafter(float(4 * tick * local_error), inf)
        unrestricted, all_vacuum = receiver(path, source)
        local = read_bounds(unrestricted, wave_error)
        heralded = read_bounds(all_vacuum, wave_error)
        assert local["reference_choi_partial_transpose_minimum"][0] > 0
        assert heralded["reference_choi_partial_transpose_minimum"][1] < 0
        assert heralded["worst_input_success"][0] > 0
        saved.append(
            {
                "appointment": float(tick * STEP),
                "wave_frobenius_error_upper": wave_error,
                "unrestricted": local,
                "spectator_vacua_herald": heralded,
            }
        )
    return {
        "fan_count": path.size,
        "signed_load": path.sign,
        "physical_bridges": [list(bridge) for bridge in path.bridges],
        "exact_h_row_sum_maximum": str(norm),
        "maximum_terms_per_row": width,
        "step_truncation_error_upper": nextafter(float(truncation), inf),
        "step_roundoff_error_upper": nextafter(float(roundoff), inf),
        "reads": saved,
    }


def run():
    with archived_fixture() as a:
        cases = [run_case(FanPath(a, n, 1, asymmetric=True)) for n in (3, 4)]
    return {
        "certificate_sha256": sha256(
            (HERE / "ASYMMETRIC_CLOCK.md").read_bytes()
        ).hexdigest(),
        "step": str(STEP),
        "degree": DEGREE,
        "cases": cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as out:
            out.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        compare = run_path(str(HERE.parent / "replay.py"))["compare"]
        compare(
            result,
            json.loads((HERE / "asymmetric_clock.json").read_text()),
            "asymmetric_clock",
        )
        for case in result["cases"]:
            for read in case["reads"]:
                print(
                    case["fan_count"],
                    read["appointment"],
                    read["unrestricted"]["reference_choi_partial_transpose_minimum"],
                    read["spectator_vacua_herald"][
                        "reference_choi_partial_transpose_minimum"
                    ],
                )


if __name__ == "__main__":
    main()
