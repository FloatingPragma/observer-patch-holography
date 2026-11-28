"""Identify a two-head modular score and vary a declared overlap recovery channel."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from ..selected_collar import measure as head, matched_readings, observer_records
from ..screen_rulers import measure as rulers

HERE = Path(__file__).parent
RECEIPT = HERE / "results.json"
SHOTS = 2**42
SEED = 190827
KEEP = (0.0, 0.5, 1.0)
MODES = ("original", "flagged", "unflagged")


def norm1(state):
    return float(sum(abs(np.linalg.eigvalsh(state)).flat))


def product_reference(ref):
    left, right = head.marginals(ref)
    return np.array(
        [np.kron(a, b) / np.trace(a) for a, b in zip(left, right, strict=True)]
    )


def projected(operator, basis):
    return np.einsum("i,ibjk->bjk", head.station.pair(basis, operator), basis)


def calibration():
    ref, *_ = head.jacobi.beam.aperture.source_and_changes()
    direction, ends, basis, _, source = head.select_source(ref)
    product = product_reference(ref)
    logs = head.station.log_blocks
    score = logs(product) - logs(ref)
    visible = projected(score, basis)
    mins = np.linalg.eigvalsh(ref).min(), np.linalg.eigvalsh(product).min()
    matched = json.loads(matched_readings.RECEIPT.read_text())["source"]
    dmatch = np.einsum("i,ibjk->bjk", matched["selected_coefficients"], basis)
    held, _, _ = rulers.preparations()
    contrasts = {
        "selected": ends[1] - ends[0],
        "matched": 2 * matched["half_step"] * dmatch,
        "even_exchange": held["held_plus"] - held["held_minus"],
    }
    pair = head.station.pair
    tangents = {}
    for name, delta in contrasts.items():
        residual = max(
            abs(head.marginals(delta)[0]).max(), abs(head.marginals(delta)[1]).max()
        )
        np.testing.assert_allclose(pair(score - visible, delta), 0, atol=1e-13)
        tangents[name] = {
            "max_marginal_change": float(residual),
            "delta_full_score": float(pair(score, delta)),
            "delta_measured_projection": float(pair(visible, delta)),
            "delta_product_generator": float(pair(-logs(product), delta)),
        }
    recipes = {}
    for name, state in (("reference", ref), ("conditional_product", product)):
        try:
            rec = head.station.preparation(state)
        except ArithmeticError:
            recipes[name] = {"realizable": False}
        else:
            recipes[name] = {
                "realizable": min(rec["min_kept"], rec["min_dephased"]) >= 0
                and rec["max_error"] < 1e-11,
                **rec,
            }
    audit = {
        "selected_step": source["step"],
        "source_vector": source["basis_coefficients"],
        "reference_min_eigenvalue": float(mins[0]),
        "product_min_eigenvalue": float(mins[1]),
        "score_eigenvalues_by_parity": np.linalg.eigvalsh(score).tolist(),
        "score_norm": float(abs(np.linalg.eigvalsh(score)).max()),
        "omitted_operator_norm": float(abs(np.linalg.eigvalsh(score - visible)).max()),
        "log_commutator_norm": float(
            max(
                np.linalg.norm(m, ord=2)
                for m in logs(product) @ logs(ref) - logs(ref) @ logs(product)
            )
        ),
        "mean_full_on_reference": float(pair(score, ref)),
        "mean_projection_on_reference": float(pair(visible, ref)),
        "mean_full_on_product": float(pair(score, product)),
        "contrast_across_sources": tangents,
        "preparations": recipes,
    }
    return ref, product, basis, ends, score, audit


def read(rng, basis, ends, p):
    records = []
    for x in ends:
        settings = []
        for i, op in enumerate(basis):
            w = np.trace(x[i // 9]).real
            mean = head.station.pair(op, x)
            outcome = np.array([w / 2 + mean, w / 2 - mean, 1 - w])
            law = np.outer([1 - p, p], outcome)
            settings.append(rng.multinomial(SHOTS, law.ravel()).reshape(2, 3).tolist())
        records.append(settings)
    return records


def decode(records, coefficients, p):
    counts = np.asarray(records, dtype=np.int64)
    if (
        counts.shape != (2, 18, 2, 3)
        or counts.min() < 0
        or not np.all(counts.sum(axis=(2, 3)) == SHOTS)
    ):
        raise ValueError("Incomplete flag and joint-head outcome bank")
    if (p == 0 and np.any(counts[:, :, 1])) or (p == 1 and np.any(counts[:, :, 0])):
        raise ValueError("Outcome outside flag support")
    delta = counts[1] - counts[0]
    # Each joint signed product has magnitude 1/2. Original and unflagged
    # scores use the marginal counts, flagged uses only replacement branches.
    by_flag = (delta[:, :, 0] - delta[:, :, 1]) / (2 * SHOTS)
    # Bound both signs, marginalized or restricted to no flag, for every
    # endpoint, setting and retention. One Hoeffding event covers all scores.
    events = 2 * 2 * 2 * 18 * len(KEEP)
    radius = np.sqrt(np.log(2 * events / 0.01) / (2 * SHOTS))
    out = {}
    for mode, co in coefficients.items():
        co = np.asarray(co)
        z = by_flag[:, 0] if mode == "flagged" else by_flag.sum(axis=1)
        center = float(co @ z)
        half = 2 * abs(co).sum() * radius
        if mode == "flagged" and p == 1:
            half = 0
        out[mode] = {"estimate": center, "halfwidth": float(half)}
    return out


def run():
    ref, product, basis, ends, score, audit = calibration()
    delta = ends[1] - ends[0]
    base = head.station.pair(basis, score)
    m = min(audit["reference_min_eigenvalue"], audit["product_min_eigenvalue"])
    rng = np.random.default_rng(SEED)
    gates = {}
    for p, t in zip(KEEP, (0, np.pi / 2, np.pi), strict=True):
        h = observer_records.blank_clock(-1, 8)
        native = float(abs(expm(-1j * t * h)[1, 0]) ** 2)
        np.testing.assert_allclose(native, p, atol=2e-15, rtol=0)
        restored = p * ref + (1 - p) * product
        unflagged = head.station.log_blocks(restored) - head.station.log_blocks(ref)
        # On common direct-sum support the scalar log-weights cancel. The
        # retained branch has zero score; the erased branch has K_R-K_P.
        eta = 0.5 * (1 - p) * norm1(ref - product)
        np.testing.assert_allclose(0.5 * norm1(ref - restored), eta, atol=2e-15)
        truth = {
            "original": float(head.station.pair(score, delta)),
            "flagged": float((1 - p) * head.station.pair(score, delta)),
            "unflagged": float(head.station.pair(unflagged, delta)),
        }
        bound = 2 * norm1(delta) * eta / m
        if any(abs(truth[k]) > bound + 1e-14 for k in ("flagged", "unflagged")):
            raise ValueError("Finite logarithm continuity bound failed")
        co = {
            "original": base,
            "flagged": base,
            "unflagged": head.station.pair(basis, unflagged),
        }
        records = read(rng, basis, ends, p)
        co = {k: c.tolist() for k, c in co.items()}
        gates[str(p)] = {
            "native_keep_probability": native,
            "reference_trace_distance": eta,
            "active_flag_support": [
                k for k, prob in (("retained", p), ("erased", 1 - p)) if prob
            ],
            "finite_absolute_contrast_bound": bound,
            "expected": truth,
            "pauli_coefficients": co,
            "counts": records,
            "measured": decode(records, co, p),
        }
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "seed": SEED,
        "shots_per_setting": SHOTS,
        "operator_audit": audit,
        "gates": gates,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--counts-only", action="store_true")
    args = parser.parse_args()
    if args.counts_only:
        data = json.loads(RECEIPT.read_text())
        if (
            data["protocol_sha256"]
            != sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest()
            or data["shots_per_setting"] != SHOTS
            or set(data["gates"]) != {str(p) for p in KEEP}
        ):
            raise ValueError("Wrong instrument or support")
        for key, g in data["gates"].items():
            if (
                decode(g["counts"], g["pauli_coefficients"], float(key))
                != g["measured"]
            ):
                raise ValueError("Recomputed scores differ")
    else:
        data = run()
        if args.write:
            with RECEIPT.open("x") as f:
                json.dump(data, f, indent=2, allow_nan=False)
                f.write("\n")
        elif data != json.loads(RECEIPT.read_text()):
            raise ValueError("Acquisition differs from receipt")
    print(
        json.dumps(
            {
                k: {
                    "eta": g["reference_trace_distance"],
                    "expected": g["expected"],
                    "observed": g["measured"],
                }
                for k, g in data["gates"].items()
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
