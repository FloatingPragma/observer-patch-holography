"""Compare material counterfactuals, certified register sets and direct writers."""

import argparse
import gzip
import hashlib
import json
from itertools import product
from pathlib import Path

from . import apparatus
from .apparatus import Board, decode, dispatch_record, final_record, registers

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
WORDS = tuple(product((0, 1), repeat=6))


def mask(reads):
    return sum(1 << i for i in reads)


def minimal_table(rows, diagnostic=False):
    """Enumerate every subcube containing a given initial word, using only readbacks."""
    if len(rows) != 64 or {tuple(r["source"]) for r in rows} != set(WORDS):
        raise ValueError("Incomplete physical intervention bank")
    reported = {}
    for row in rows:
        result = dispatch_record(decode(row["dispatch_raw"][0]))[0]
        if diagnostic:
            result = result, row["dispatcher_raw"]
        reported[mask(i for i, b in enumerate(row["source"]) if b)] = result
    certificates = []
    for base in range(64):
        diffs = [
            base ^ other for other in range(64) if reported[other] != reported[base]
        ]
        sufficient = [s for s in range(64) if all(d & s for d in diffs)]
        minima = [
            s for s in sufficient if all((s & k) != k for k in sufficient if k != s)
        ]
        sensitive = sum(
            1 << i for i in range(6) if reported[base ^ (1 << i)] != reported[base]
        )
        certificates.append(
            {"least": minima, "single_flip": sensitive, "sufficient": sufficient}
        )
    return certificates


def authenticated(register, support, parent_id, parent_writes, post, child_before):
    """Test the register named in a particular producer/consumer comparison."""
    return (
        register in support
        and register in parent_writes
        and child_before[register][1] == parent_id
        and child_before[register][0] == post[register][0]
    )


def lineage(run):
    snapshots = [registers(s) for s in run["audit_snapshots"]]
    operations = run["operations"]
    # Validate all mutations against the independently sensed frames and writer stamps.
    edges = []
    for i, (event, written) in enumerate(operations):
        before, after = snapshots[i : i + 2]
        if after[written][1] != event or any(
            before[r] != after[r] for r in range(8) if r != written
        ):
            raise ValueError("A sealed commit failed frame or stamp")
        support = (
            run["declared_offer_support"] if event == 8 else ([6] if event == 9 else [])
        )
        if event == 8:
            if decode(run["dispatcher_raw"]) != decode(run["audit_snapshots"][i][:24]):
                raise ValueError("The child did not read its audited source versions")
        if (
            event == 9
            and dispatch_record(decode(run["continuation_read"])) != before[6]
        ):
            raise ValueError(
                "The continuation did not read its audited dispatch version"
            )
        versions = [[r, *before[r]] for r in support]
        if json.dumps(run["certified_versions"][str(event)]) != json.dumps(versions):
            raise ValueError("Producer declared a different value or writer version")
        for reg in support:
            matches = [
                e
                for j, (e, w) in enumerate(operations[:i])
                if authenticated(reg, support, e, {w}, snapshots[j + 1], before)
            ]
            if len(matches) != 1:
                raise ValueError(
                    "Certified live register lacks a unique material writer"
                )
            edges.append([matches[0], event, reg])
    return edges


def selected_control(row):
    """Keep the three other conditions true while breaking one authentication condition."""
    i = next(i for i, (ev, _) in enumerate(row["operations"]) if ev == 8)
    versions = [registers(raw) for raw in row["audit_snapshots"]]
    before = versions[i]
    support = row["declared_offer_support"]
    reg = 2 + before[1][0]
    original = reg + 1
    if before[reg][1] != 7:
        raise ValueError("The controls need a replaced selected writer")
    excluded = 3 if reg == 2 else 2
    # In this held source both banks have the same payload.
    if before[excluded][0] != before[reg][0]:
        raise ValueError("Unequal controls")
    altered_post = versions[7].copy()
    altered_post[reg] = (1 - before[reg][0], 7)
    good = authenticated(reg, support, 7, {reg}, versions[7], before)
    # For a wrong declared write, the read's tag and value are deliberately kept.
    checks = {
        "honest": good,
        "replacement_with_wrong_write_set": authenticated(
            reg, support, 7, {excluded}, versions[7], before
        ),
        "outdated_equal_writer": authenticated(
            reg, support, original, {reg}, versions[original], before
        ),
        "changed_parent_value": authenticated(
            reg, support, 7, {reg}, altered_post, before
        ),
        "uncertified_equal_register": authenticated(
            excluded, support, excluded + 1, {excluded}, versions[excluded + 1], before
        ),
    }
    if not good or any(v for k, v in checks.items() if k != "honest"):
        raise ValueError("A separated provenance control passed")
    return checks


def report():
    board = Board()
    banks = {}
    for mode, replace in (
        ("native", None),
        ("native", 2),
        ("native", 3),
        ("hidden", None),
    ):
        rows = []
        for word in WORDS:
            row = board.run(word, replace, mode)
            row["reported_edges"] = lineage(row)
            rows.append(row)
        table = minimal_table(rows)
        diagnostic = minimal_table(rows, True)
        for row in rows:
            idx = mask(i for i, x in enumerate(row["source"]) if x)
            expected = mask(row["declared_offer_support"])
            resolved = table[idx]
            if mode == "hidden" and row["source"][0]:
                expected |= 1 << 5
            if resolved["least"] != [expected] or resolved["single_flip"] != expected:
                raise ValueError(
                    "Measured dependence differs from the predicted intervention"
                )
            if diagnostic[idx]["least"] != [63]:
                raise ValueError("Raw pointer tape failed to resolve an input")
        banks[f"{mode}_{replace}"] = {
            "runs": rows,
            "certificates_by_input_mask": table,
            "diagnostic_least_by_input_mask": [r["least"] for r in diagnostic],
        }
    held = next(
        r for r in banks["native_2"]["runs"] if r["source"] == [1, 0, 1, 1, 0, 1]
    )
    checks = selected_control(held)
    parked = board.run((1, 0, 1, 1, 0, 1), 2, "parked")
    if set(parked["dispatcher_raw"]) != {"0"} or registers(
        parked["independent_final_snapshot"]
    )[-1] != (None, 0):
        raise ValueError("The parked sensor reported a usable result")
    counts = {
        "nominal_accept": 0,
        "nominal_refuse": 0,
        "authenticated_edges": 0,
        "same_value_new_parent": 0,
        "unnecessary_spare_rejected": 0,
        "missing_onward_rejected": 0,
        "hidden_missing_spare": 0,
    }
    base = banks["native_None"]["runs"]
    for name, bank in banks.items():
        for reference, row in zip(base, bank["runs"], strict=True):
            src = row["source"]
            idx = mask(i for i, bit in enumerate(src) if bit)
            least = bank["certificates_by_input_mask"][idx]
            support = mask(row["declared_offer_support"])
            if name.startswith("hidden"):
                counts["hidden_missing_spare"] += support not in least["sufficient"]
                continue
            counts["authenticated_edges"] += len(row["reported_edges"])
            if not src[0]:
                counts["nominal_refuse"] += 1
                continue
            counts["nominal_accept"] += 1
            counts["unnecessary_spare_rejected"] += (support | 32) not in least["least"]
            counts["missing_onward_rejected"] += (support & ~16) not in least[
                "sufficient"
            ]
            if (
                row["dispatch_raw"] != reference["dispatch_raw"]
                or row["continuation_raw"] != reference["continuation_raw"]
            ):
                raise ValueError(
                    "Equal-valued source replacements changed public deliveries"
                )
            counts["same_value_new_parent"] += (
                row["reported_edges"] != reference["reported_edges"]
            )
    if tuple(counts.values()) != (96, 96, 480, 32, 96, 96, 32):
        raise ValueError(f"Unexpected intervention or ancestry counts: {counts}")
    files = [
        HERE / "PROTOCOL.md",
        HERE / "apparatus.py",
        Path(__file__),
        ROOT
        / "model/experiments/interval_echo_consensus/boolean_nodes/repair_ballots/measure.py",
        ROOT / "model/experiments/interval_echo_consensus/boolean_nodes/material.py",
        ROOT / "model/candidate_laws/interval_echo/candidate.py",
        ROOT / "model/candidate_laws/contrast_memory_spherical_v1/law.py",
        ROOT / "model/executable_model/carrier.py",
        ROOT / "model/executable_model/carriers/FEDERATION_LEVEL0.npz",
    ]
    return {
        "schema": "oph.versioned_dispatch.v1",
        "source_keys": apparatus.NAMES,
        "raw_digit": "two balanced-pointer clicks, most significant first; 2=zero, 1=one, 0 and 3 refuse",
        "provenance": "Deterministic sharp outcomes of routed native operators in simulation.",
        "sources": {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in files
        },
        "banks": banks,
        "auth_controls": checks,
        "counts": counts,
        "parked": parked,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = report()
    encoded = json.dumps(data, separators=(",", ":")).encode()
    archive = HERE / "records.json.gz"
    if args.write:
        archive.write_bytes(gzip.compress(encoded, mtime=0))
    if json.loads(encoded) != json.loads(gzip.decompress(archive.read_bytes())):
        raise ValueError("Material intervention receipt differs")
    for name, bank in data["banks"].items():
        accepted = sum(
            final_record(decode(r["continuation_raw"][0]))[0] is not None
            for r in bank["runs"]
        )
        print(name, len(bank["runs"]), accepted)
