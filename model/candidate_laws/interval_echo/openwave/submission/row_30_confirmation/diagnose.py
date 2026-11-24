"""Read the frozen v2 pilot and validation counts; change only the inferred k.

Run from the repository root with the directory containing pilot/ and validation/.
The default checks diagnosis.json. No shots, new board settings or source
generator spectra enter this diagnostic.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from ...row_30_collective import holdout
from ...row_30_collective.confirmation import acquire, spatial
from ...row_30_collective.experiment import verify


def _gate_scores(board, wave_numbers, mass, speed2):
    modes = board["observer"]["modes"]
    rows = [
        dict(k=k, omega=r["omega"]) for r, k in zip(modes, wave_numbers, strict=True)
    ]
    measured = holdout.score(
        rows,
        board["observer"],
        board["transport"],
        board["changed_transport"],
        mass,
        speed2,
    )
    k = np.asarray(wave_numbers)
    omega = np.asarray([r["omega"] for r in rows])
    low = (k > 0) & (k <= 0.50)
    continuum = omega - np.sqrt(mass**2 + speed2 * k**2)
    max_low = float(np.max(abs(continuum[low]))) if np.any(low) else None
    # All other gates depend only on fixed meter, timing, commissioning or
    # transport records. Their frozen values were authenticated with the raw bank.
    gates = dict(board["gates"])
    gates.update(
        lattice=bool(np.all(omega > 0) and measured["lattice_max"] < 0.0048),
        shape=bool(
            measured["shape_rms_lattice"] < 0.5 * measured["shape_rms_parabola"]
        ),
        continuum=bool(sum(low) >= 2 and np.ptp(k[low]) >= 0.02 and max_low < 0.0048),
    )
    return gates, dict(
        lattice_max=measured["lattice_max"],
        low_k_continuum_max=max_low,
        low_k_count=int(sum(low)),
        lattice_over_parabola_rms=measured["shape_rms_lattice"]
        / measured["shape_rms_parabola"],
    )


def _pair_contributions(m, k, tests):
    residuals = spatial._residuals(m, k, tests).reshape(len(tests), m)
    pairs = {}
    for (a, b, _, _), residual in zip(tests, residuals, strict=True):
        pair = (a + 1, b + 1)
        pairs[pair] = pairs.get(pair, 0.0) + float(residual @ residual)
    return pairs


def board_diagnosis(board, folder, frozen, mass, speed2):
    n = len(board["route"])
    if acquire.audit_counts(folder, acquire.expected_counts(n)) != frozen["counts"]:
        raise ValueError("saved counts differ from the frozen inventory")
    original = np.array([r["k"] for r in board["observer"]["modes"]])
    # v3 exposes the two estimators separately. Pass the old wave numbers
    # explicitly and require every field in the original v2 spatial decision.
    baseline = spatial.score(board, folder, original)
    verify({k: baseline[k] for k in frozen["spatial"]}, frozen["spatial"])
    fit = np.array(baseline["even_trained_k"])
    old_gates, old_values = _gate_scores(board, original, mass, speed2)
    if old_gates != board["gates"] or old_gates != frozen["native_gates"]:
        raise ValueError("original verdict did not reproduce")
    new_gates, new_values = _gate_scores(board, fit, mass, speed2)
    m, _, _, _, odd, _ = spatial._coordinates(board, folder)
    old_pairs = _pair_contributions(m, original, odd)
    new_pairs = _pair_contributions(m, fit, odd)
    old_q, new_q = sum(old_pairs.values()), sum(new_pairs.values())
    np.testing.assert_allclose(
        [old_q, new_q], [baseline["candidate_q"], baseline["even_trained_q"]]
    )
    new_gates.update(
        profile=bool(new_q <= baseline["quadratic_cutoff"]),
        launch_label_rejected=baseline["launch_label_rejected"],
    )
    modes = []
    for i, (old, fitted) in enumerate(zip(original, fit, strict=True)):
        hybrid = original.copy()
        hybrid[i] = fitted
        gain = old_q - spatial._quadratic(m, hybrid, odd)
        modes.append(
            [
                i + 1,
                round(float(old), 10),
                round(float(fitted), 10),
                round(float(fitted - old), 10),
                round(float(gain), 5),
            ]
        )
    ranked = sorted(
        old_pairs, key=lambda pair: old_pairs[pair] - new_pairs[pair], reverse=True
    )
    return dict(
        level=board["level"],
        seed=board["seed"],
        dimers=n,
        frozen_status=frozen["status"],
        hypothetical_all_gates_pass=all(new_gates.values()),
        hypothetical_failed_gates=[g for g, ok in new_gates.items() if not ok],
        odd_q=dict(
            original=old_q,
            even_fit=new_q,
            launch_labels=baseline["launch_label_q"],
            cutoff=baseline["quadratic_cutoff"],
        ),
        dispersion=dict(original=old_values, even_fit=new_values),
        modes=modes,
        largest_pair_q_reductions=[
            dict(
                pair=list(pair),
                original=round(old_pairs[pair], 5),
                even_fit=round(new_pairs[pair], 5),
                reduction=round(old_pairs[pair] - new_pairs[pair], 5),
            )
            for pair in ranked[:3]
        ],
    )


def diagnose(root):
    v2_package = json.loads((Path(__file__).parent / "v2" / "PACKAGE.json").read_text())
    v2_seeds = {"pilot": [927016] * 3, "validation": [374621, 614357, 801173]}
    out = dict(
        scope="Post-hoc reanalysis of two existing local banks. Frozen v2 decisions are preserved.",
        method="Even-hold fitted k replaces the original full-recurrence k in the candidate profile, lattice, continuum and shape gates. Every other physical input is fixed.",
        mode_columns=[
            "launch",
            "original_k",
            "even_fit_k",
            "fit_minus_original_k",
            "odd_Q_reduction_replacing_only_this_mode",
        ],
        attribution="One-mode gains are conditional on keeping every other original k, so they need not add. Pair gains sum across both held receiver scans per pair.",
        gates=sorted(acquire.GATE_KEYS | {"profile", "launch_label_rejected"}),
        banks={},
    )
    for kind in ("pilot", "validation"):
        release = root / kind
        manifest = json.loads((release / "manifest.json").read_text())
        draw = json.loads((release / "draw.json").read_text())
        if draw != manifest["draw"] or draw["binding"]["mode"] != kind:
            raise ValueError("wrong release type or draw")
        binding = draw["binding"]
        if (
            draw["seeds"] != v2_seeds[kind]
            or binding["anchor"] is not None
            or binding["source_commit"] != v2_package["source_commit"]
            or binding["candidate_pin"] != v2_package["source_commit"]
            or binding["candidate_sha256"] != v2_package["candidate_sha256"]
            or binding["protocol_sha256"] != v2_package["task_sha256"]
        ):
            raise ValueError("local selectors or declaration differ")
        records = release / "records"
        for name, digest in manifest["metadata_sha256"].items():
            if hashlib.sha256((records / name).read_bytes()).hexdigest() != digest:
                raise ValueError("physical record metadata differ")
        raw = json.loads((records / "result.json").read_text())
        boards = []
        for level, (b, saved) in enumerate(
            zip(raw["boards"], manifest["boards"], strict=True)
        ):
            if level != b["level"] or b["seed"] != draw["seeds"][level]:
                raise ValueError("board ID differs")
            boards.append(
                board_diagnosis(
                    b, records / str(level), saved, raw["m"], raw["c_squared"]
                )
            )
        out["banks"][kind] = dict(
            source_commit=draw["binding"]["source_commit"],
            frozen_status=manifest["status"],
            boards=boards,
        )
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("release_root", type=Path)
    parser.add_argument("--output", type=Path, help="write a new diagnostic receipt")
    args = parser.parse_args()
    result = diagnose(args.release_root)
    if args.output:
        with args.output.open("x") as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        verify(
            result, json.loads(Path(__file__).with_name("diagnosis.json").read_text())
        )
        print("Reproduced all six retrospective board diagnostics.")
