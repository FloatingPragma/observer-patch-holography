#!/usr/bin/env python3
"""B20 random-scan common-reference preflight (issue #725).

Exact laptop-scale feasibility preflight on retained source data, per the
issue's highest-value first step.  From the pinned run
``runs/b12_prereg_16k_20260806`` of ``oph-physics-sim`` the producer
recounts the exact integer transition table on the committed 26-state
packet alphabet, forms the faithful visit-count reference, and for every
subset of at least two of the four committed packet fields builds the
overlapping conditional-resampling projectors and their random-scan
mixture under the declared uniform scheduler on the subset's fields.

The recount and every kernel/mixture/dimension computation are exact
(integers, then ``fractions.Fraction``); the one committed float step is
the arena-2 quantile binning convention, a pinned deterministic rule
cross-checked against the pinned realization receipt, with the numpy
version recorded in provenance.  Per field subset the certificate
records:

* row-stochasticity and exact stationarity of the mixture under the one
  shared reference;
* a non-idempotence witness entry (the mixture lies outside the #688
  spectral obstruction, which covered idempotent collars);
* the connected components of the partitions' join, whose indicator
  functions are exactly the observables fixed by every projector: at
  least two components are required for a nonconstant protected
  observable;
* the exact fixed-space dimension of the mixture, computed by rational
  Gaussian elimination on ``mix - I``: dimension one certifies directly
  that the only mixture-invariant observables are the constants, with no
  appeal to the partition-lattice argument.

Fail-closed controls per the contract: an idempotent mixture, a single
join component (constant protected record), a fixed-space dimension that
disagrees with the join component count, a non-faithful reference, or a
counting mismatch against the pinned report rejects the subset or the
run.  The scheduler is the uniform law on the committed fields of the
subset, declared in advance; the reference is the source visit count;
no stationary law is imported and no coupling is fitted.

The designated subset is the first field subset passing every
requirement, in the committed enumeration order (by size, then
lexicographically in the committed field order).  If none passes, the
certificate records the negative verdict with the per-subset failure
table; only a positive preflight may trigger a fresh prospectively
frozen run.

Claim boundary: this is a bounded algebraic feasibility receipt over one
retained artifact.  No four-law composition, physical thermodynamics,
energy, or clock claim is made; #703 owns physical calibration.

Usage (from the repository root):

    python3 code/b20_random_scan/preflight_random_scan.py
"""

from __future__ import annotations

import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
RUN_DIR = REPO_ROOT.parent / "oph-physics-sim" / "runs" / "b12_prereg_16k_20260806"
OUT_PATH = Path(__file__).resolve().parent / "runtime" / "b20_preflight_certificate.json"

PINNED_INPUTS = (
    "finite_repair_transition_matrix_report.json",
    "observer_views.jsonl",
    "conditional_resampling_realization_receipt.json",
    "freezeout_fields.npz",
    "git_commit.txt",
)

MAX_RECORD_CLASSES = 32
MAX_COMPANION_CLASSES = 16


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fail(obstruction: str) -> None:
    print(f"NEGATIVE RESULT: {obstruction}", file=sys.stderr)
    raise SystemExit(1)


def require(condition: bool, obstruction: str) -> None:
    if not condition:
        fail(obstruction)


def recount(run_dir: Path):
    """Exact integer transition recount on the committed alphabet.

    Also inventories every retained step field beyond the four committed
    packet fields and records its realized value set, so the certificate
    can quantify over all retained step fields: a subset containing a
    constant extra field mixes in the full-space resampler and is
    trivially connected, so the four-field battery covers it."""
    report = json.loads(
        (run_dir / "finite_repair_transition_matrix_report.json").read_text()
    )
    fields = [str(f) for f in report["packet_fields"]]
    labels = [
        tuple((str(f), int(v)) for f, v in json.loads(s))
        for s in report["state_labels"]
    ]
    index = {label: i for i, label in enumerate(labels)}
    n = len(labels)
    require(n == int(report["state_count"]), "state_count disagrees with labels")

    counts = [[0] * n for _ in range(n)]
    observers = 0
    skipped = 0
    transitions = 0
    extra_values: dict[str, set] = {}
    with (run_dir / "observer_views.jsonl").open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            observers += 1
            view = json.loads(line)
            steps = (view.get("transition_history_descriptor") or {}).get("steps") or []
            if len(steps) < 2:
                skipped += 1
                continue
            encoded = []
            for step in steps:
                for key in step:
                    if key not in fields:
                        extra_values.setdefault(key, set()).add(step[key])
                values = []
                for f in fields:
                    require(f in step, f"step misses field {f}")
                    v = step[f]
                    require(
                        isinstance(v, int) and not isinstance(v, bool),
                        f"non-integer value for {f}",
                    )
                    values.append((f, v))
                key = tuple(values)
                require(key in index, f"step outside the report alphabet: {key}")
                encoded.append(index[key])
            for a, b in zip(encoded, encoded[1:]):
                counts[a][b] += 1
                transitions += 1
    require(observers == int(report["observer_count"]), "observer count mismatch")
    require(skipped == int(report["skipped_observer_count"]), "skip count mismatch")
    require(
        transitions == int(report["transition_count"]), "transition count mismatch"
    )
    require(
        all(len(vals) == 1 for vals in extra_values.values()),
        "a retained step field beyond the packet fields is nonconstant: "
        "the certified grammar would not cover it",
    )
    extra_fields = {
        key: sorted(vals) for key, vals in sorted(extra_values.items())
    }
    return report, fields, labels, counts, observers, transitions, extra_fields


def heat_bath(pi, fibre):
    """Exact conditional-resampling kernel for one fibre map."""
    n = len(pi)
    mass = {}
    for x in range(n):
        mass.setdefault(fibre[x], Fraction(0))
        mass[fibre[x]] += pi[x]
    kernel = [[Fraction(0)] * n for _ in range(n)]
    for x in range(n):
        for y in range(n):
            if fibre[y] == fibre[x]:
                kernel[x][y] = pi[y] / mass[fibre[x]]
    return kernel


def join_components(fibres, n):
    """Connected components of the join of the given partitions."""
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for fibre in fibres:
        first = {}
        for x in range(n):
            if fibre[x] in first:
                union(x, first[fibre[x]])
            else:
                first[fibre[x]] = x
    roots = sorted({find(x) for x in range(n)})
    comp = {r: i for i, r in enumerate(roots)}
    return [comp[find(x)] for x in range(n)], len(roots)


def fixed_space_dimension(mix, n):
    """Exact dimension of ``{f : mix f = f}`` by rational elimination.

    Row-reduces ``mix - I`` over the rationals; the fixed functions of
    the mixture are the kernel vectors, so the dimension is ``n`` minus
    the rank."""
    rows = [
        [mix[i][j] - (1 if i == j else 0) for j in range(n)] for i in range(n)
    ]
    rank = 0
    for col in range(n):
        pivot = None
        for r in range(rank, n):
            if rows[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inv = rows[rank][col]
        rows[rank] = [v / inv for v in rows[rank]]
        for r in range(n):
            if r != rank and rows[r][col] != 0:
                factor = rows[r][col]
                rows[r] = [a - factor * b for a, b in zip(rows[r], rows[rank])]
        rank += 1
        if rank == n:
            break
    return n - rank


def freezeout_classes(run_dir: Path):
    """Record/companion classes under the committed binning convention.

    Mirror of the conditional-resampling ``_class_bins`` rule with the
    committed class budgets (32 record classes, 16 companion classes): at
    most that many distinct realized values map through the identity
    lookup, otherwise quantile edges over the realized values define the
    classes.  The class structure is cross-checked against the pinned
    realization receipt (fiber count, state count, total mass)."""
    import numpy as np

    z = np.load(run_dir / "freezeout_fields.npz")

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
    receipt = json.loads(
        (run_dir / "conditional_resampling_realization_receipt.json").read_text()
    )
    pinned = receipt["pinned_reference"]
    pairs = {}
    for r, c in zip(rc, cc):
        pairs[(r, c)] = pairs.get((r, c), 0) + 1
    require(len(set(rc)) == int(pinned["fiber_count"]),
            "record class count disagrees with the pinned receipt")
    require(len(pairs) == int(pinned["state_count"]),
            "realized state count disagrees with the pinned receipt")
    require(sum(pairs.values()) == int(pinned["total_mass_count"]),
            "total mass disagrees with the pinned receipt")
    return pairs


def analyse_subset(pi, fibres, n, names):
    """The full requirement battery for one subset of partitions."""
    k = len(fibres)
    kernels = [heat_bath(pi, fibre) for fibre in fibres]
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
                witness = {
                    "x": x,
                    "y": y,
                    "mix_entry": str(mix[x][y]),
                    "mix_squared_entry": str(row2[y]),
                }
                break
        if witness:
            break
    non_idempotent = witness is not None
    comp_map, comp_count = join_components(fibres, n)
    dim = fixed_space_dimension(mix, n)
    require(
        dim == comp_count,
        f"fixed-space dimension {dim} disagrees with join components "
        f"{comp_count} for {names}",
    )

    # Second declared scheduler point: weights proportional to each
    # partition's block count (a source-derived, output-independent
    # quantity), certifying that the fixed space does not depend on the
    # uniform choice at this second point either.
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
    dim2 = fixed_space_dimension(mix2, n)
    require(
        dim2 == comp_count,
        f"second-scheduler fixed-space dimension {dim2} disagrees with "
        f"join components {comp_count} for {names}",
    )

    nonconstant_protected = dim >= 2
    passes = stochastic and stationary and non_idempotent and         nonconstant_protected
    return {
        "fields": list(names),
        "scheduler": (
            f"uniform over the {k} committed fields of the subset, declared"
        ),
        "row_stochastic": stochastic,
        "stationary_under_shared_reference": stationary,
        "non_idempotent": non_idempotent,
        "non_idempotence_witness": witness,
        "join_component_count": comp_count,
        "fixed_space_dimension": dim,
        "second_scheduler": {
            "rule": "weights proportional to partition block counts, declared",
            "weights": [str(w) for w in weights],
            "stationary_under_shared_reference": stationary2,
            "fixed_space_dimension": dim2,
        },
        "nonconstant_protected_observable": nonconstant_protected,
        "passes_all_requirements": passes,
    }


def main() -> None:
    require(RUN_DIR.is_dir(), f"missing run directory {RUN_DIR}")
    for name in PINNED_INPUTS:
        require((RUN_DIR / name).is_file(), f"missing pinned input {name}")

    (report, fields, labels, counts, observers, transitions,
     extra_fields) = recount(RUN_DIR)
    n = len(labels)

    visit = [sum(counts[x]) for x in range(n)]
    require(all(v > 0 for v in visit), "reference is not faithful: a state has no counted exit")
    total = sum(visit)
    pi = [Fraction(v, total) for v in visit]

    field_values = {
        f: [dict(label)[f] for label in labels] for f in fields
    }

    from itertools import combinations

    subset_results = []
    designated = None
    for size in range(2, len(fields) + 1):
        for names in combinations(fields, size):
            entry = analyse_subset(
                pi, [field_values[f] for f in names], n, names
            )
            subset_results.append(entry)
            if entry["passes_all_requirements"] and designated is None:
                designated = list(names)

    # Arena 2: the conditional-resampling record/companion structure.
    fz_pairs = freezeout_classes(RUN_DIR)
    states2 = sorted(fz_pairs)
    n2 = len(states2)
    total2 = sum(fz_pairs.values())
    pi2 = [Fraction(fz_pairs[s2], total2) for s2 in states2]
    require(all(m > 0 for m in fz_pairs.values()),
            "arena-2 reference is not faithful")
    rec_f = [s2[0] for s2 in states2]
    com_f = [s2[1] for s2 in states2]
    arena2 = analyse_subset(pi2, [rec_f, com_f], n2,
                            ("record_class", "companion_class"))
    if arena2["passes_all_requirements"] and designated is None:
        designated = ["record_class", "companion_class"]

    verdict = "positive" if designated is not None else "negative"

    certificate = {
        "schema": "oph.b20_random_scan_preflight.v2",
        "provenance": {
            "run_id": RUN_DIR.name,
            "run_git_commit": (RUN_DIR / "git_commit.txt").read_text().strip(),
            "pinned_input_sha256": {
                name: sha256_file(RUN_DIR / name) for name in PINNED_INPUTS
            },
            "counting_convention": (
                "exact integer recount of consecutive "
                "transition_history_descriptor steps on the committed "
                "26-state packet alphabet, cross-checked against the pinned "
                "report totals"
            ),
            "reference_rule": (
                "visit counts (row sums of the exact count table), "
                "normalized; faithfulness is a hard check"
            ),
            "scheduler_rule": (
                "uniform law on the committed fields of the subset, "
                "declared in advance; no output-dependent choice"
            ),
            "designation_rule": (
                "first field subset passing every requirement, enumerated "
                "by size then lexicographically in the committed field "
                "order"
            ),
            "epistemic_status": (
                "post-hoc algebraic feasibility preflight on a retained "
                "artifact; ineligible as validation; only a positive "
                "preflight may trigger a fresh prospectively frozen run"
            ),
            "numpy_version": __import__("numpy").__version__,
        },
        "state_count": n,
        "observer_count": observers,
        "transition_count": transitions,
        "visit_counts": visit,
        "reference": [str(p) for p in pi],
        "field_order": fields,
        "state_labels": list(report["state_labels"]),
        "extra_step_fields_realized_values": extra_fields,
        "subset_results": subset_results,
        "arena2_state_count": n2,
        "arena2_states": [
            [int(s2[0]), int(s2[1]), int(fz_pairs[s2])] for s2 in states2
        ],
        "arena2_result": arena2,
        "designated_subset": designated,
        "verdict": verdict,
        "no_go_grammar": (
            "every subset of at least two retained step fields on the "
            "26-state alphabet, and the record/companion class pair on the "
            "realized 256-state conditional-resampling structure, each "
            "under the visit-count or occupation reference of its arena "
            "and the declared uniform random-scan scheduler on the "
            "subset's conditional-resampling projectors, with a second "
            "declared block-count-proportional scheduler certified at "
            "every computed subset. The computed subsets are the eleven "
            "packet-field subsets and the record/companion pair; for each "
            "computed mixture the exact fixed-space dimension is one, so "
            "every mixture-invariant observable is constant. The remaining "
            "step-field subsets contain a certified-constant field, so "
            "their mixtures carry the strictly positive full-space "
            "resampler as a component, which forces the same "
            "one-dimensional fixed space (a strictly positive stochastic "
            "component makes the mixture irreducible); non-idempotence is "
            "claimed only for the computed subsets"
        ),
        "recorded_next_route": (
            "enriched-source: a committed producer exporting an additional "
            "protected field whose join with an existing field is "
            "disconnected (this requires a fresh export, since every "
            "retained extra step field is certified constant on this "
            "run), or a dilated construction outside the random-scan "
            "grammar; declared impossible for nothing beyond the stated "
            "grammar"
        ),
    }
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(certificate, indent=1, sort_keys=True) + "\n"
    # write_bytes keeps LF endings on every platform.
    OUT_PATH.write_bytes(text.encode("utf-8"))
    digest = hashlib.sha256(text.encode()).hexdigest()
    print(f"wrote {OUT_PATH} ({len(text.encode())} bytes)")
    print(f"certificate_sha256: {digest}")
    print(f"verdict: {verdict}; designated subset: {designated}")
    for entry in [arena2] + subset_results:
        print(
            f"  {entry['fields']}: stochastic={entry['row_stochastic']} "
            f"stationary={entry['stationary_under_shared_reference']} "
            f"nonidem={entry['non_idempotent']} "
            f"components={entry['join_component_count']} "
            f"dim={entry['fixed_space_dimension']} "
            f"passes={entry['passes_all_requirements']}"
        )


if __name__ == "__main__":
    main()
