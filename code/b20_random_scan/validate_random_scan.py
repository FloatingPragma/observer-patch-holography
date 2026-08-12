#!/usr/bin/env python3
"""Independent validator for the B20 random-scan preflight certificate.

Two verification layers, following the vendored-payload house pattern:

* **Offline algebra layer (always runs).**  The certificate vendors the
  26 committed state labels, the exact visit counts, and the realized
  256-row record/companion class table.  From those alone the validator
  recomputes, through its own code path (independent connectivity
  routine, elimination with a different pivot strategy than the
  producer), every per-subset battery entry: row-stochasticity, exact
  stationarity, the non-idempotence witness entries, the join component
  count, the fixed-space dimensions under both declared schedulers
  (which must also equal the component count — two independent routes to
  the same invariant), the protected-observable verdict, the pass flag,
  the subset enumeration, and the verdict/designation logic.
* **Custody layer (runs when the pinned run directory is present).**
  Re-derives the vendored objects from the pinned run: hash-checks every
  pinned input, recounts visits and the extra-step-field inventory from
  ``observer_views.jsonl``, rebuilds the record/companion classes from
  ``freezeout_fields.npz`` under the committed binning convention, and
  fails closed on any disagreement with the vendored copies.  Without
  the run directory the custody layer is reported as skipped; the
  pinned SHA-256 digests still bind the certificate to the run.

Semantic tamper classes rejected offline (exercised by the test suite):
a subset promoted to passing, a component count altered, a fixed-space
dimension altered (either scheduler), a witness entry altered, a
visit-count edit, a state-label edit, a verdict flip, a reference
fraction edit, and an emptied extra-field inventory.

Usage:

    python3 code/b20_random_scan/validate_random_scan.py [certificate]
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
RUN_DIR = Path(os.environ.get(
    "B20_RUN_DIR",
    str(REPO_ROOT.parent / "oph-physics-sim" / "runs" / "b12_prereg_16k_20260806"),
))
DEFAULT_CERT = HERE / "runtime" / "b20_preflight_certificate.json"

MAX_RECORD_CLASSES = 32
MAX_COMPANION_CLASSES = 16

PINNED_INPUT_NAMES = {
    "finite_repair_transition_matrix_report.json",
    "observer_views.jsonl",
    "conditional_resampling_realization_receipt.json",
    "freezeout_fields.npz",
    "git_commit.txt",
}

# SHA-256 of the load-bearing certificate prose.  A tampered grammar,
# route, or status sentence must fail closed, so the expected digests are
# pinned here; regenerating the certificate with changed prose requires a
# deliberate validator update.
EXPECTED_PROSE_SHA256 = {
    "no_go_grammar":
        "6fdc9a03469cd0ca3f4b4a981ba1b21ab2010974889dae17e53761dca551530a",
    "recorded_next_route":
        "97d6f89026c4dd10cfe7b35f23d3d73f1f2412279bdfd5ace01942f09d2f3cb7",
    "epistemic_status":
        "581a3fb5419a0710ac353e835a682d32c601f81c7078a01af5dd72ad1c59aaa2",
}


def fail(reason: str) -> None:
    print(f"VALIDATION FAILED: {reason}", file=sys.stderr)
    raise SystemExit(1)


def check(cond: bool, reason: str) -> None:
    if not cond:
        fail(reason)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def kernel(pi, fibre, n):
    mass = {}
    for x in range(n):
        mass[fibre[x]] = mass.get(fibre[x], Fraction(0)) + pi[x]
    return [
        [pi[y] / mass[fibre[x]] if fibre[y] == fibre[x] else Fraction(0)
         for y in range(n)]
        for x in range(n)
    ]


def components_multi(fibres, n) -> int:
    """Join components by iterative label propagation (not union-find)."""
    labels = list(range(n))
    changed = True
    while changed:
        changed = False
        group = {}
        for x in range(n):
            for tag, fibre in enumerate(fibres):
                key = (tag, fibre[x])
                if key in group:
                    if labels[x] != labels[group[key]]:
                        m = min(labels[x], labels[group[key]])
                        labels[x] = m
                        labels[group[key]] = m
                        changed = True
                else:
                    group[key] = x
        # propagate within groups
        rep = {}
        for x in range(n):
            for tag, fibre in enumerate(fibres):
                key = (tag, fibre[x])
                rep.setdefault(key, labels[x])
                if labels[x] != rep[key]:
                    m = min(labels[x], rep[key])
                    if labels[x] != m or rep[key] != m:
                        changed = True
                    labels[x] = m
                    rep[key] = m
    return len(set(labels))


def fixed_dimension(mix, n):
    """Kernel dimension of ``mix - I`` by elimination, pivoting from the
    bottom row upward (a different strategy than the producer's)."""
    rows = [
        [mix[i][j] - (1 if i == j else 0) for j in range(n)] for i in range(n)
    ]
    rank = 0
    for col in range(n - 1, -1, -1):
        pivot = None
        for r in range(n - 1 - rank, -1, -1):
            if rows[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            continue
        target = n - 1 - rank
        rows[target], rows[pivot] = rows[pivot], rows[target]
        inv = rows[target][col]
        rows[target] = [v / inv for v in rows[target]]
        for r in range(n):
            if r != target and rows[r][col] != 0:
                factor = rows[r][col]
                rows[r] = [
                    a - factor * b for a, b in zip(rows[r], rows[target])
                ]
        rank += 1
        if rank == n:
            break
    return n - rank


def analyse(pi, fibres, n):
    kernels = [kernel(pi, fibre, n) for fibre in fibres]
    k = len(kernels)
    mix = [
        [sum((ker[x][y] for ker in kernels), Fraction(0)) / k for y in range(n)]
        for x in range(n)
    ]
    stochastic = all(sum(row, Fraction(0)) == 1 for row in mix)
    stationary = all(
        sum((pi[x] * mix[x][y] for x in range(n)), Fraction(0)) == pi[y]
        for y in range(n)
    )
    witness = None
    for x in range(n):
        row2 = [
            sum((mix[x][z] * mix[z][y] for z in range(n)), Fraction(0))
            for y in range(n)
        ]
        for y in range(n):
            if row2[y] != mix[x][y]:
                witness = (x, y, mix[x][y], row2[y])
                break
        if witness:
            break
    comp = components_multi(fibres, n)
    dim = fixed_dimension(mix, n)

    blocks = [len(set(fibre)) for fibre in fibres]
    wtotal = sum(blocks)
    weights = [Fraction(b, wtotal) for b in blocks]
    mix2 = [
        [
            sum((w * ker[x][y] for w, ker in zip(weights, kernels)),
                Fraction(0))
            for y in range(n)
        ]
        for x in range(n)
    ]
    stationary2 = all(
        sum((pi[x] * mix2[x][y] for x in range(n)), Fraction(0)) == pi[y]
        for y in range(n)
    )
    dim2 = fixed_dimension(mix2, n)
    second = {
        "weights": [str(w) for w in weights],
        "stationary": stationary2,
        "dim": dim2,
    }
    return stochastic, stationary, witness, comp, dim, second


def check_entry(entry, pi, fibres, n, tag):
    stoch, stat, witness, comp, dim, second = analyse(pi, fibres, n)
    check(entry["scheduler"] == (
        f"uniform over the {len(fibres)} committed fields of the subset, "
        f"declared"), f"{tag}: scheduler declaration")
    check(entry["second_scheduler"]["rule"] == (
        "weights proportional to partition block counts, declared"),
        f"{tag}: second-scheduler declaration")
    check(entry["row_stochastic"] == stoch, f"{tag}: stochastic")
    check(entry["stationary_under_shared_reference"] == stat,
          f"{tag}: stationary")
    check(entry["non_idempotent"] == (witness is not None),
          f"{tag}: idempotence flag")
    if witness is not None:
        w = entry["non_idempotence_witness"]
        check(
            (w["x"], w["y"]) == (witness[0], witness[1])
            and w["mix_entry"] == str(witness[2])
            and w["mix_squared_entry"] == str(witness[3]),
            f"{tag}: witness values",
        )
    check(entry["join_component_count"] == comp, f"{tag}: components")
    check(entry["fixed_space_dimension"] == dim, f"{tag}: fixed dimension")
    check(dim == comp, f"{tag}: dimension/component cross-check")
    sched2 = entry["second_scheduler"]
    check(sched2["weights"] == second["weights"],
          f"{tag}: second-scheduler weights")
    check(sched2["stationary_under_shared_reference"] == second["stationary"],
          f"{tag}: second-scheduler stationarity")
    check(sched2["fixed_space_dimension"] == second["dim"],
          f"{tag}: second-scheduler dimension")
    check(second["dim"] == comp,
          f"{tag}: second-scheduler dimension/component cross-check")
    check(entry["nonconstant_protected_observable"] == (dim >= 2),
          f"{tag}: protected verdict")
    expected_pass = stoch and stat and (witness is not None) and dim >= 2
    check(entry["passes_all_requirements"] == expected_pass,
          f"{tag}: pass flag")
    return expected_pass


def offline_layer(cert):
    """Recompute every battery entry from the vendored objects alone."""
    fields = cert["field_order"]
    n = cert["state_count"]
    labels = [
        tuple((str(f), int(v)) for f, v in json.loads(s))
        for s in cert["state_labels"]
    ]
    check(len(labels) == n, "vendored label count disagrees with state count")
    check(len(set(labels)) == n, "vendored labels are not distinct")
    check(
        all(tuple(f for f, _v in label) == tuple(fields) for label in labels),
        "vendored label fields disagree with the committed field order",
    )
    visits = cert["visit_counts"]
    check(len(visits) == n and all(
        isinstance(v, int) and v > 0 for v in visits),
        "visit counts must be positive integers (faithful reference)")
    total = sum(visits)
    check(cert["transition_count"] == total,
          "transition count must equal the summed visit counts")
    check(isinstance(cert["observer_count"], int)
          and cert["observer_count"] >= 1, "observer count must be positive")
    pi = [Fraction(v, total) for v in visits]
    check([str(p) for p in pi] == cert["reference"], "reference mismatch")
    extra = cert["extra_step_fields_realized_values"]
    check(len(extra) >= 1 and all(
        len(vals) == 1 for vals in extra.values()),
        "extra step fields must be inventoried and constant "
        "(pinned to this run's step schema)")
    prose = {
        "no_go_grammar": cert["no_go_grammar"],
        "recorded_next_route": cert["recorded_next_route"],
        "epistemic_status": cert["provenance"]["epistemic_status"],
    }
    for field, text_value in prose.items():
        check(
            hashlib.sha256(text_value.encode()).hexdigest()
            == EXPECTED_PROSE_SHA256[field],
            f"load-bearing prose field tampered: {field}",
        )

    field_values = {f: [dict(label)[f] for label in labels] for f in fields}
    from itertools import combinations

    expected_subsets = [
        list(names)
        for size in range(2, len(fields) + 1)
        for names in combinations(fields, size)
    ]
    check([e["fields"] for e in cert["subset_results"]] == expected_subsets,
          "subset enumeration mismatch")

    expected_designated = None
    for entry in cert["subset_results"]:
        fibres = [field_values[f] for f in entry["fields"]]
        passed = check_entry(entry, pi, fibres, n, str(entry["fields"]))
        if passed and expected_designated is None:
            expected_designated = list(entry["fields"])

    rows = cert["arena2_states"]
    n2 = cert["arena2_state_count"]
    check(len(rows) == n2, "arena-2 state count disagrees with table")
    check(len({(r[0], r[1]) for r in rows}) == n2,
          "arena-2 class pairs are not distinct")
    check(all(int(r[2]) > 0 for r in rows),
          "arena-2 reference is not faithful")
    total2 = sum(int(r[2]) for r in rows)
    pi2 = [Fraction(int(r[2]), total2) for r in rows]
    rec_f = [int(r[0]) for r in rows]
    com_f = [int(r[1]) for r in rows]
    a2 = cert["arena2_result"]
    check(a2["fields"] == ["record_class", "companion_class"],
          "arena2 field names")
    passed2 = check_entry(a2, pi2, [rec_f, com_f], n2, "arena2")
    if passed2 and expected_designated is None:
        expected_designated = ["record_class", "companion_class"]

    check(cert["designated_subset"] == expected_designated,
          "designated subset disagrees with the committed designation rule")
    expected_verdict = "positive" if expected_designated else "negative"
    check(cert["verdict"] == expected_verdict, "verdict mismatch")
    return len(cert["subset_results"]), n2


def custody_layer(cert):
    """Re-derive the vendored objects from the pinned run directory."""
    pinned = cert["provenance"]["pinned_input_sha256"]
    check(set(pinned) == PINNED_INPUT_NAMES,
          "pinned input set disagrees with the required five run inputs")
    for name, digest in pinned.items():
        check(sha256_file(RUN_DIR / name) == digest,
              f"pinned input hash mismatch: {name}")
    check(
        cert["provenance"]["run_git_commit"]
        == (RUN_DIR / "git_commit.txt").read_text().strip(),
        "recorded run git commit disagrees with the pinned run",
    )

    report = json.loads(
        (RUN_DIR / "finite_repair_transition_matrix_report.json").read_text()
    )
    check(cert["field_order"] == [str(f) for f in report["packet_fields"]],
          "field order disagrees with the pinned report")
    check(cert["state_labels"] == list(report["state_labels"]),
          "vendored state labels disagree with the pinned report")

    fields = cert["field_order"]
    labels = [
        tuple((str(f), int(v)) for f, v in json.loads(s))
        for s in cert["state_labels"]
    ]
    index = {label: i for i, label in enumerate(labels)}
    n = len(labels)
    visits = [0] * n
    extra: dict[str, set] = {}
    with (RUN_DIR / "observer_views.jsonl").open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            view = json.loads(line)
            steps = (view.get("transition_history_descriptor") or {}).get("steps") or []
            if len(steps) < 2:
                continue
            for step in steps:
                for key in step:
                    if key not in fields:
                        extra.setdefault(key, set()).add(step[key])
            try:
                encoded = [
                    index[tuple((f, int(step[f])) for f in fields)]
                    for step in steps
                ]
            except KeyError:
                fail("a run step lies outside the vendored alphabet")
            for a, _b in zip(encoded, encoded[1:]):
                visits[a] += 1
    check(cert["visit_counts"] == visits, "visit counts mismatch")
    check(
        cert["extra_step_fields_realized_values"]
        == {key: sorted(vals) for key, vals in sorted(extra.items())},
        "extra step field inventory mismatch",
    )

    import numpy as np

    z = np.load(RUN_DIR / "freezeout_fields.npz")

    def bins(values, max_classes):
        flat = np.asarray(values).ravel()
        distinct = np.unique(flat)
        if distinct.size <= max_classes:
            lookup = {v: i for i, v in enumerate(distinct.tolist())}
            return [lookup[v] for v in flat.tolist()]
        edges = np.quantile(flat.astype(float),
                            np.linspace(0.0, 1.0, max_classes + 1)[1:-1])
        return np.digitize(flat.astype(float), edges).astype(int).tolist()

    rc = bins(z["record_signature"], MAX_RECORD_CLASSES)
    cc = bins(z["cumulative_repair_load"], MAX_COMPANION_CLASSES)
    pairs = {}
    for r, c in zip(rc, cc):
        pairs[(r, c)] = pairs.get((r, c), 0) + 1
    table = [[r, c, pairs[(r, c)]] for r, c in sorted(pairs)]
    check(cert["arena2_states"] == table,
          "arena-2 class table disagrees with the pinned run")


def main() -> None:
    args = [a for a in sys.argv[1:] if a != "--require-custody"]
    require_custody = "--require-custody" in sys.argv[1:]
    cert_path = Path(args[0]) if args else DEFAULT_CERT
    cert = json.loads(cert_path.read_text())
    check(cert.get("schema") == "oph.b20_random_scan_preflight.v2",
          "unexpected schema")

    n_subsets, n2 = offline_layer(cert)

    if RUN_DIR.is_dir():
        custody_layer(cert)
        custody = "custody layer verified against the pinned run"
    else:
        check(not require_custody,
              "custody layer required but the pinned run directory is absent")
        custody = ("custody layer skipped: pinned run directory absent; "
                   "SHA-256 digests still bind the certificate to the run")

    print("VALIDATION PASSED:", cert_path.name)
    print(f"  verdict {cert['verdict']}; {n_subsets} field subsets + arena2 "
          f"({n2} states) recomputed exactly from the vendored objects; "
          f"fixed-space dimensions match join components everywhere")
    print(f"  {custody}")


if __name__ == "__main__":
    main()
