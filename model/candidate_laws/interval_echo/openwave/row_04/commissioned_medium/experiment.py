"""Release the two row 4 loads on the archived row 30 spectroscopy materials."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.executable_model import federation
from model.executable_model.carrier import Carrier
from ...row_30_formed import experiment as donor
from .. import open_pair, open_pair_check as filed, open_pair_release as release
from ..corridor_pair import LAW, certify_native_collar, neutral_column
from ..corridor_release import monitored_source

HERE = Path(__file__).resolve().parent
DONOR_PIN = "44104ae607c93382f6115181337e26846ac868d9"
DONOR_SHA = "e6ba479b4ac47224bba800f9303851d5e07fa5e3d73116c827ada0e5be6669bf"
BANKS = (178121, 178122)
X, Y, Z = open_pair.BASIS
METER = (156, 161)
MARGIN = 0.94
FILED_DENSITY = float(filed.RETRY) * np.outer(filed.VECTOR, filed.VECTOR)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def materials(board):
    archive = donor.HERE / "results.json"
    if digest(archive) != DONOR_SHA:
        raise ValueError("Spectroscopy archive differs from the pinned material")
    studies = [
        s
        for s in json.loads(archive.read_text())["studies"]
        if s["cells"] == board.cells
    ]
    if tuple(s["seed"] for s in studies) != BANKS:
        raise ValueError("Incomplete or unexpected 20-cell spectroscopy roster")
    fair = donor.native()
    for study in studies:
        bank = study["seed"]
        history = study["formation"]
        regenerated = donor.formation(board, bank, 0, fair)[-1]
        if history != regenerated:
            raise ValueError(f"Bank {bank} fails the donor's own formation replay")
        packed = bytes.fromhex(history["coin_pairs_packed_little_bit_first_hex"])
        bits = np.unpackbits(
            np.frombuffer(packed, dtype=np.uint8),
            bitorder="little",
            count=2 * len(board.seams),
        )
        closed = bits.reshape(-1, 2).all(axis=1)
        written = {e for e, mark in zip(board.seams, closed, strict=True) if mark}
        if len(written) != history["written"]:
            raise ValueError("Decoded a different number of neutral writes")
        yield bank, written, sha256(packed).hexdigest()


def check_meter():
    """Check each selective occupation operation, including restoration of the memory."""
    certify_native_collar()  # Includes the terminal pair flag on blank and neutral seams.
    echo = LAW.echo()
    blank = np.arange(LAW.DIM) < LAW.N_SRC
    nonzero = np.repeat([False] + [d != 0 for _, d in LAW.LABELS], LAW.N_SRC)
    for occupied in (0, 1):
        word = (occupied, 0)
        sharp = np.eye(LAW.DIM)[:, LAW.SRC.index(word)]
        for state, click in ((sharp, nonzero), (neutral_column(word), ~blank)):
            for seen, flag in ((1, click), (0, ~click)):
                np.testing.assert_allclose(
                    echo.T @ (flag * (echo @ state)),
                    state * (occupied == seen),
                    rtol=0,
                    atol=2e-14,
                )


def acquire(board, written, kappa):
    """Compute the actual two addressed terms on X,Y,Z and monitor the center."""
    native = LAW.seam_term(kappa)
    block = np.zeros((3, 3))
    for i, pair in enumerate(open_pair.BASIS):
        for edge in (X, Z):
            if edge in written:
                continue
            word = tuple(int(v in pair) for v in edge)
            j = LAW.SRC.index(word)
            block[i, i] += native[j, j]
            if sum(word) == 1:
                after = tuple(sorted(set(pair) ^ set(edge)))
                block[open_pair.BASIS.index(after), i] += native[
                    LAW.SRC.index(word[::-1]), j
                ]
        block[i, i] += kappa * (6 - board.degree[156]) * (156 in pair)
    step = expm((-2j * np.pi / 3) * block)
    # A vacancy is Y; the reset to X exists exactly when the Z contact can swap.
    attempts = 1 if Z in written else release.ATTEMPTS
    rho, failed, reach, first = monitored_source(step, release.ROUNDS, attempts)
    accepted = rho[np.ix_([0, 2], [0, 2])]
    if X not in written and Z not in written:
        np.testing.assert_allclose(accepted, FILED_DENSITY, atol=2e-13, rtol=0)
    return accepted, {
        "first_vacancy_per_attempt": first,
        "attempt_reach": reach,
        "exhaustion": failed if Z not in written else 0.0,
        "failed_with_reset_unavailable": failed if Z in written else 0.0,
        "accepted_density_XZ_real": accepted.real.tolist(),
        "accepted_density_XZ_imag": accepted.imag.tolist(),
    }


def generator(board, written, home):
    live = Carrier(
        "row 4 flight contacts", board.n, [e for e in board.seams if e not in written]
    )
    pairs = PairSector(live)
    exchange = pairs.exchange()
    cost = live.degree + (6 - board.degree)
    contrast = cost[pairs.aa] + cost[pairs.bb] - 2 * pairs.adjacent
    launch = np.zeros((pairs.P, 2))
    launch[[pairs.index(*X), pairs.index(*Z)], [0, 1]] = 1
    detector = [pairs.index(*edge) for edge in home]
    return exchange, contrast, launch, detector


def home_bins(amplitude, rho):
    return np.einsum("ij,jk,ik->i", amplitude, rho, amplitude.conj()).real


def terminal(home):
    outcome = np.r_[home, 1 - home.sum()]
    if (
        not np.isfinite(outcome).all()
        or outcome.min() < -1e-12
        or abs(outcome.sum() - 1) > 2e-12
    ):
        raise ArithmeticError("Invalid home or quiet probability")
    return np.maximum(0.0, outcome).tolist()


def check_short_flight(exchange, contrast, kappa, launch, detector, amplitude):
    # The Chebyshev replay drops a common energy. Compare effects to remove its phase.
    h = exchange + sparse.diags(kappa * contrast, dtype=float, format="csr")
    full = expm_multiply(-1j * filed.TIMES[0] * h, launch)
    np.testing.assert_allclose(full.conj().T @ full, np.eye(2), rtol=0, atol=2e-11)
    direct = full[detector]
    expected = direct[:, :, None].conj() * direct[:, None, :]
    actual = amplitude[:, :, None].conj() * amplitude[:, None, :]
    np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-11)
    return float(np.max(abs(actual - expected)))


def study(board, home, bank, written, coin_hash):
    exchange, contrast, launch, detector = generator(board, written, home)
    stochastic = filed.stochastic_population(exchange, filed.TIMES, launch, detector)
    records, sources = [], []
    errors = []
    for kappa in filed.STRENGTHS:
        rho, source = acquire(board, written, kappa)
        sources.append({"kappa": kappa, **source})
        amplitudes = filed.flight_amplitudes(
            exchange, contrast, kappa, filed.TIMES, launch, detector
        )
        errors.append(
            check_short_flight(
                exchange, contrast, kappa, launch, detector, amplitudes[0]
            )
        )
        for i, (t, a) in enumerate(zip(filed.TIMES, amplitudes, strict=True)):
            bins = {
                "acquired": terminal(home_bins(a, rho)),
                "filed_input": terminal(home_bins(a, FILED_DENSITY)),
                "stochastic_filed_input": terminal(stochastic[i]),
            }
            best = float(np.linalg.eigvalsh(a.conj().T @ a)[-1])
            total = {arm: sum(v[:-1]) for arm, v in bins.items()}
            mass = rho.trace().real
            failed = source["exhaustion"] + source["failed_with_reset_unavailable"]
            if (
                abs(mass + failed - 1) > 2e-12
                or bins["acquired"][-1] + 2e-12 < failed
                or total["acquired"] > mass * best + 2e-12
                or best > 1 + 2e-12
            ):
                raise ArithmeticError(
                    "Home effect or preparation denominator is inconsistent"
                )
            records.append(
                {
                    "kappa": kappa,
                    "exposure": t,
                    "bins": bins,
                    "home": total,
                    "best_normalized_XZ_home": best,
                }
            )
    return {
        "bank": bank,
        "packed_coin_sha256": coin_hash,
        "writes": len(written),
        "pulse_seams_written_XZ": [X in written, Z in written],
        "meter_written": METER in written,
        "acquisitions": sources,
        "appointments": records,
        "largest_t1_effect_check_error": max(errors),
        "all_acquired_home_at_least_0_94": all(
            r["home"]["acquired"] >= MARGIN for r in records
        ),
    }


def controls(board, home):
    original = json.loads(filed.ARCHIVE.read_text())
    archived = filed.archive_bins(original, home)
    written = set(open_pair.geometry()[0])
    exchange, contrast, launch, detector = generator(board, written, home)
    worst = 0.0
    for kappa in filed.STRENGTHS:
        amplitudes = filed.flight_amplitudes(
            exchange, contrast, kappa, filed.TIMES, launch, detector
        )
        for t, a in zip(filed.TIMES, amplitudes, strict=True):
            bins = np.array(terminal(home_bins(a, FILED_DENSITY)))
            expected = archived[kappa, t]["coherent"]
            np.testing.assert_allclose(bins, expected, atol=2e-10, rtol=0)
            worst = max(worst, float(np.max(abs(bins - expected))))
    exchange, contrast, launch, detector = generator(board, set(), home)
    (a,) = filed.flight_amplitudes(exchange, contrast, 0.5, [37], launch, detector)
    blank = float(sum(home_bins(a, FILED_DENSITY)))
    np.testing.assert_allclose(blank, 0.440929105, rtol=0, atol=6e-10)
    return {
        "largest_filed_31_bin_error": worst,
        "blank_filed_input_home_k0_5_t37": blank,
    }


def run():
    board = federation(0)
    _, home = open_pair.geometry()
    if (
        digest(filed.BOARD) != filed.BOARD_HASH
        or METER not in board.seams
        or METER[1] in open_pair.PATH
    ):
        raise ValueError("Unexpected board or middle-port acquisition meter")
    check_meter()
    return {
        "protocol_sha256": digest(HERE / "PROTOCOL.md"),
        "board_sha256": filed.BOARD_HASH,
        "row_30_commit": DONOR_PIN,
        "row_30_receipt_sha256": DONOR_SHA,
        "source_order": [list(X), list(Z)],
        "meter": list(METER),
        "outcome_order": [*map(list, home), "quiet"],
        "controls": controls(board, home),
        "studies": [study(board, home, *material) for material in materials(board)],
    }


def check_record(actual, expected):
    if isinstance(actual, dict):
        if actual.keys() != expected.keys():
            raise ValueError("Archived record fields differ from replay")
        for key in actual:
            check_record(actual[key], expected[key])
    elif isinstance(actual, list):
        for value, saved in zip(actual, expected, strict=True):
            check_record(value, saved)
    elif isinstance(actual, float):
        if type(expected) not in (float, int) or not np.isclose(
            actual, expected, atol=2e-11, rtol=0
        ):
            raise ValueError(f"Saved value {expected} differs from replay {actual}")
    elif type(actual) is not type(expected) or actual != expected:
        raise ValueError(f"Saved metadata {expected} differs from replay {actual}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--verify", action="store_true", help="replay and compare the saved record"
    )
    args = parser.parse_args()
    result = run()
    if args.verify:
        check_record(result, json.loads((HERE / "results.json").read_text()))
        print(
            "Reproduced both materials, nine appointments per material and the filed controls."
        )
    else:
        print(json.dumps(result, indent=2, allow_nan=False))
