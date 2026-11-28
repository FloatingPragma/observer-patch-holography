"""Read the full connected modular increment with an independently counted reference."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import measure

AXES = np.asarray(measure.station.phase.AXES)
PRODUCTS = np.array([[np.kron(a, b) for b in AXES] for a in AXES])
SIGNS = np.array([1, -1])
SHOTS = (2**24, 2**38)


def emit(state, n, rng):
    tape = []
    for a in AXES[1:]:
        row = []
        for b in AXES[1:]:
            effects = [
                np.kron((AXES[0] + s * a) / 2, (AXES[0] + t * b) / 2)
                for s in SIGNS
                for t in SIGNS
            ]
            probs = np.einsum("kij,bji->bk", effects, state).real
            assert probs.min() >= 0 and abs(probs.sum() - 1) < 1e-12
            row.append(rng.multinomial(n, probs.ravel()).reshape(2, 2, 2).tolist())
        tape.append(row)
    return tape


def tomography(tape, n):
    counts = np.array(tape)
    if (
        counts.shape != (3, 3, 2, 2, 2)
        or n <= 0
        or not np.issubdtype(counts.dtype, np.integer)
        or counts.min() < 0
        or not np.all(counts.sum(axis=(2, 3, 4)) == n)
    ):
        raise ValueError(
            "Every axis pair must contain the stated number of signed parity outcomes"
        )
    p = counts / n
    moments = np.empty((2, 4, 4))
    moments[:, 0, 0] = p[0, 0].sum(axis=(1, 2))
    moments[:, 1:, 0] = np.einsum("ibst,s->bi", p[:, 0], SIGNS)
    moments[:, 0, 1:] = np.einsum("jbst,t->bj", p[0], SIGNS)
    moments[:, 1:, 1:] = np.einsum("ijbst,s,t->bij", p, SIGNS, SIGNS)
    return np.einsum("bij,ijkl->bkl", moments, PRODUCTS) / 4


def product_reference(ref):
    left, right = measure.marginals(ref)
    return np.array(
        [np.kron(a, b) / np.trace(a) for a, b in zip(left, right, strict=True)]
    )


def decode(tapes, n):
    # 32 Pauli coordinates in each of four states, at both precisions.
    # Every coordinate is a per-shot score in [-1,1], including unselected parity.
    delta = np.sqrt(2 * np.log(2 * 32 * 4 * len(SHOTS) / 0.01) / n)
    matrix_error = 2 * delta
    states = {k: tomography(v, n) for k, v in tapes.items()}
    diff = states["plus"] - states["minus"]
    # Each 4x4 block has Frobenius error <= 2 delta by Pauli orthogonality.
    # For a source difference, each block has trace norm error <= 4*matrix_error.
    size = sum(np.linalg.norm(x, ord="nuc") for x in diff) + 8 * matrix_error
    output = {}
    for name in ("equilibrium", "decorrelated"):
        ref = states[name]
        left, right = measure.marginals(ref)
        parity = np.trace(ref, axis1=1, axis2=2).real

        def log_error(smallest, error):
            if smallest <= error:
                raise ValueError("Reference faithfulness unresolved at this exposure")
            return error / (smallest - error)

        # Use the log resolvent integral and the lower bound on the spectra
        # of both calibration matrices along the line segment between them.
        op_error = max(
            log_error(np.linalg.eigvalsh(ref[b]).min(), matrix_error)
            + log_error(np.linalg.eigvalsh(left[b]).min(), 2 * matrix_error)
            + log_error(np.linalg.eigvalsh(right[b]).min(), 2 * matrix_error)
            + log_error(parity[b], delta)
            for b in range(2)
        )
        operator = measure.station.log_blocks(
            product_reference(ref)
        ) - measure.station.log_blocks(ref)
        value = measure.station.pair(operator, diff)
        coefficients = np.einsum("ijkl,blk->bij", PRODUCTS, operator).real / 4
        src_error = 2 * delta * abs(coefficients).sum()
        ref_error = op_error * size
        width = ref_error + src_error
        output[name] = dict(
            increment=value,
            source_error=src_error,
            reference_error=ref_error,
            interval=[value - width, value + width],
            positive=bool(value > width),
        )
    return output


def run():
    ref, _, _, _ = measure.jacobi.beam.aperture.source_and_changes()
    _, ends, _, _, _ = measure.select_source(ref)
    states = dict(
        equilibrium=ref,
        decorrelated=product_reference(ref),
        minus=ends[0],
        plus=ends[1],
    )
    rng = np.random.default_rng(81931)
    cases = {}
    for n in SHOTS:
        tapes = {k: emit(v, n, rng) for k, v in states.items()}
        cases[str(n)] = dict(tapes=tapes, readings=decode(tapes, n))
    return dict(
        seed=81931,
        protocol_sha256=sha256(
            (measure.HERE / "REFERENCE_READ.md").read_bytes()
        ).hexdigest(),
        order="left axis, right axis, parity, left sign, right sign",
        cases=cases,
    )


if __name__ == "__main__":
    cmd = argparse.ArgumentParser(description=__doc__)
    modes = cmd.add_mutually_exclusive_group()
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--counts-only", action="store_true")
    args = cmd.parse_args()
    path = measure.HERE / "reference_read.json"
    if args.counts_only:
        saved = json.loads(path.read_bytes())
        for n in SHOTS:
            case = saved["cases"][str(n)]
            readings = decode(case["tapes"], n)
            from cosmology.rar_deep_regime.lensing_profile import check

            check(case["readings"], readings)
            print(n, readings)
        raise SystemExit(0)
    out = json.loads(json.dumps(run(), allow_nan=False))
    if args.write:
        path.write_text(json.dumps(out, indent=2) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), out)
    for k, result in out["cases"].items():
        print(k, result["readings"])
