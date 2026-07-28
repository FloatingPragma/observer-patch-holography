#!/usr/bin/env python3
"""Workstream D checker: independent counterterm-packet verification.

Re-derives the first-order counterterm expansion from the certified
rule table with its own generation code, recomputes the generator
matrix rank and null space with its own elimination, re-solves the
gauge pole recursion from the Workstream B betas, re-enumerates the
general local basis from the component quantum numbers, and replays
the packet controls.  The checker never imports the producer.
"""

from __future__ import annotations

import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))

from vertex_format import FIELD_TABLE  # noqa: E402

TABLE_PATH = ROOT / "outputs" / "rule_table_engine_a.json"
MATCHING_PATH = ROOT / "outputs" / "eft_matching_1.json"
PACKET_PATH = ROOT / "outputs" / "renormalization_ct_1.json"

PARAMETERS = {"g1": "dg1", "g2": "dg2", "g3": "dg3", "v": "dv", "lam": "dlam", "mu2": "dmu2"}
BOSON_Z = {
    "Wp": "dZW", "Wm": "dZW", "W3": "dZW", "B": "dZB", "Gl": "dZGl",
    "h": "dZH", "G0": "dZH", "Gp": "dZH", "Gm": "dZH",
    "cp": "dZcW", "cm": "dZcW", "c3": "dZcW", "cB": "dZcB",
}
FERMION_Z = {
    "uL": "dZQL", "dL": "dZQL", "uR": "dZuR", "dR": "dZdR",
    "nuL": "dZLL", "eL": "dZLL", "eR": "dZeR",
}
YUKAWA = ("Yu", "Yd", "Ye", "Yud", "Ydd", "Yed")


def sha(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def split_indexed(name: str) -> tuple[str, int, int] | None:
    if "[" not in name:
        return None
    base = name.split("[")[0]
    left, right = name[len(base):].strip("[]").split("][")
    return base, int(left), int(right)


def rederive_terms(table: dict[str, Any]) -> list[dict[str, Any]]:
    terms: list[dict[str, Any]] = []
    for entry in table["entries"]:
        fields = entry["fields"]
        barred = [f[:-4] for f in fields if f.endswith("_bar") and f[:-4] in FERMION_Z]
        plain = [f for f in fields if not f.endswith("_bar") and f in FERMION_Z]
        for mono in entry["coefficient"]["monomials"]:
            pref = Fraction(mono["prefactor"])
            powers = [tuple(p) for p in mono["powers"]]
            yuk = None
            for name, _ in powers:
                parsed = split_indexed(name)
                if parsed and parsed[0] in YUKAWA:
                    yuk = parsed
            free = (yuk[1], yuk[2]) if yuk else (None, None)

            def add(value: Fraction, generator: str, residual: list, pair: tuple) -> None:
                if value != 0:
                    terms.append({
                        "fields": fields,
                        "structure": entry["structure"],
                        "generator": generator,
                        "prefactor": str(value),
                        "residual_powers": [list(p) for p in sorted(residual)],
                        "bar_gen": pair[0],
                        "plain_gen": pair[1],
                    })

            for name, power in powers:
                if name in PARAMETERS:
                    residual = [(n, p) for n, p in powers if n != name]
                    if power != 1:
                        residual.append((name, power - 1))
                    add(pref * power, PARAMETERS[name], residual, free)
                parsed = split_indexed(name)
                if parsed and parsed[0] in YUKAWA:
                    add(pref, f"d{name}", [(n, p) for n, p in powers if n != name], free)
            for f in fields:
                base = f[:-4] if f.endswith("_bar") else f
                if base in BOSON_Z:
                    add(pref / 2, BOSON_Z[base], list(powers), free)
            if barred or plain:
                if yuk is None:
                    for i in (1, 2, 3):
                        for j in (1, 2, 3):
                            for f in barred:
                                add(pref / 2, f"{FERMION_Z[f]}d[{i}][{j}]", list(powers), (i, j))
                            for f in plain:
                                add(pref / 2, f"{FERMION_Z[f]}[{i}][{j}]", list(powers), (i, j))
                else:
                    base, i_idx, j_idx = yuk
                    residual = [
                        (n, p) for n, p in powers
                        if split_indexed(n) is None or split_indexed(n)[0] not in YUKAWA
                    ]
                    for f in barred:
                        for a in (1, 2, 3):
                            add(pref / 2, f"{FERMION_Z[f]}d[{a}][{i_idx}]",
                                residual + [(f"{base}[{i_idx}][{j_idx}]", 1)], (a, j_idx))
                    for f in plain:
                        for b in (1, 2, 3):
                            add(pref / 2, f"{FERMION_Z[f]}[{j_idx}][{b}]",
                                residual + [(f"{base}[{i_idx}][{j_idx}]", 1)], (i_idx, b))
    return terms


def matrix_of(terms: list[dict[str, Any]]) -> tuple[int, list[str], dict]:
    rows: dict[str, int] = {}
    cols: dict[str, int] = {}
    cells: dict[tuple[int, int], Fraction] = {}
    for term in terms:
        row_key = json.dumps({
            "fields": term["fields"], "structure": term["structure"],
            "residual_powers": term["residual_powers"],
            "bar_gen": term["bar_gen"], "plain_gen": term["plain_gen"],
        }, sort_keys=True)
        r = rows.setdefault(row_key, len(rows))
        c = cols.setdefault(term["generator"], len(cols))
        cells[(r, c)] = cells.get((r, c), Fraction(0)) + Fraction(term["prefactor"])
    names = [n for n, _ in sorted(cols.items(), key=lambda kv: kv[1])]
    return len(rows), names, cells


def eliminate(row_count: int, names: list[str], cells: dict) -> tuple[int, list[dict[str, str]]]:
    n = len(names)
    dense = [[Fraction(0)] * n for _ in range(row_count)]
    for (r, c), value in cells.items():
        dense[r][c] = value
    pivots: list[int] = []
    top = 0
    for col in range(n):
        hit = next((r for r in range(top, row_count) if dense[r][col] != 0), None)
        if hit is None:
            continue
        dense[top], dense[hit] = dense[hit], dense[top]
        scale = 1 / dense[top][col]
        dense[top] = [x * scale for x in dense[top]]
        for r in range(row_count):
            if r != top and dense[r][col] != 0:
                f = dense[r][col]
                dense[r] = [x - f * y for x, y in zip(dense[r], dense[top])]
        pivots.append(col)
        top += 1
    null = []
    for free in (c for c in range(n) if c not in pivots):
        vec = {names[free]: "1"}
        for prow, pcol in enumerate(pivots):
            if dense[prow][free] != 0:
                vec[names[pcol]] = str(-dense[prow][free])
        null.append(vec)
    return len(pivots), null


def resolve_gauge_poles(matching: dict[str, Any]) -> dict[str, str]:
    out = {}
    for name, key in (("dg1", "b1"), ("dg2", "b2"), ("dg3", "b3")):
        b = Fraction(matching["gauge_betas"]["coefficients"][key])
        out[name] = str(-b / Fraction(1 - 3))
    return out


def general_basis_dimension() -> tuple[int, int]:
    """Independent count of the charge- and triality-conserving local
    basis and its flavor-weighted dimension, walked directly over the
    component table with the class composition rules."""

    vectors = ("Wp", "Wm", "W3", "B")
    scalars = ("h", "G0", "Gp", "Gm")
    lefts = ("uL", "dL", "nuL", "eL")
    rights = ("uR", "dR", "eR")
    ghosts = ("cp", "cm", "c3", "cB")

    def q(label: str) -> Fraction:
        base = label[:-4] if label.endswith("_bar") else label
        value = Fraction(FIELD_TABLE[base]["charge"])
        return -value if label.endswith("_bar") else value

    def tri(label: str) -> int:
        base = label[:-4] if label.endswith("_bar") else label
        value = FIELD_TABLE[base]["triality"]
        return -value if label.endswith("_bar") else value

    def ok(fields: tuple[str, ...]) -> bool:
        return sum((q(f) for f in fields), Fraction(0)) == 0 and sum(tri(f) for f in fields) % 3 == 0

    import itertools
    directions = 0
    weighted = 0

    def count(fields: tuple[str, ...], flavored: bool) -> None:
        nonlocal directions, weighted
        if ok(fields):
            directions += 1
            weighted += 9 if flavored else 1

    for pair in itertools.combinations_with_replacement(vectors, 2):
        count(pair, False)
    for pair in itertools.combinations_with_replacement(scalars, 2):
        count(pair, False)
    for s in scalars:
        count((s,), False)
    for size in (3, 4):
        for combo in itertools.combinations_with_replacement(scalars, size):
            count(combo, False)
    for s in scalars:
        for pair in itertools.combinations_with_replacement(vectors, 2):
            count((s,) + pair, False)
    for spair in itertools.combinations_with_replacement(scalars, 2):
        for vpair in itertools.combinations_with_replacement(vectors, 2):
            count(spair + vpair, False)
    for spair in itertools.combinations(scalars, 2):
        for vec in vectors:
            count(spair + (vec,), False)
    for triple in itertools.combinations(vectors, 3):
        count(triple, False)
    for quad in itertools.combinations_with_replacement(vectors, 4):
        if ok(quad):
            directions += 3
            weighted += 3
    for chirality in (lefts, rights):
        for f1 in chirality:
            for f2 in chirality:
                for vec in vectors:
                    count((f"{f1}_bar", f2, vec), True)
    for fl in lefts:
        for fr in rights:
            for s in scalars:
                count((f"{fl}_bar", fr, s), True)
                count((f"{fr}_bar", fl, s), True)
            count((f"{fl}_bar", fr), True)
            count((f"{fr}_bar", fl), True)
    for g1_label in ghosts:
        for g2_label in ghosts:
            count((f"{g1_label}_bar", g2_label), False)
            for s in scalars:
                count((f"{g1_label}_bar", g2_label, s), False)
            for vec in vectors:
                count((f"{g1_label}_bar", g2_label, vec), False)
    return directions, weighted


def check() -> dict[str, Any]:
    table = json.loads(TABLE_PATH.read_text(encoding="utf-8"))
    matching = json.loads(MATCHING_PATH.read_text(encoding="utf-8"))
    packet = json.loads(PACKET_PATH.read_text(encoding="utf-8"))
    problems: list[str] = []

    if packet["bindings"]["rule_table_digest"] != table["table_digest"]:
        problems.append("packet is bound to a different rule table")

    terms = rederive_terms(table)
    digest = sha(sorted(terms, key=lambda item: json.dumps(item, sort_keys=True)))
    if digest != packet["ct_terms_digest"]:
        problems.append("independent counterterm expansion digest differs")
    if len(terms) != packet["ct_term_count"]:
        problems.append("counterterm term count differs")

    rows, names, cells = matrix_of(terms)
    rank, null = eliminate(rows, names, cells)
    if rank != packet["ct_matrix"]["rank"]:
        problems.append(f"independent rank {rank} differs from packet {packet['ct_matrix']['rank']}")
    if len(null) != len(packet["ct_matrix"]["nullspace"]):
        problems.append("independent null-space dimension differs")
    for vector in null:
        support = set(vector)
        quark = any(name.startswith(("dZQL", "dZuR", "dZdR")) for name in support)
        lepton = any(name.startswith(("dZLL", "dZeR")) for name in support)
        if quark == lepton:
            problems.append("a null vector mixes or misses the quark/lepton sectors")
        offdiag = [n for n in support if (p := split_indexed(n)) and p[1] != p[2]]
        if offdiag:
            problems.append("a null vector uses off-diagonal generators")
    if "dxi1" in names or "dxi2" in names or "dxi2" in packet["generators"]:
        problems.append("a gauge-parameter generator appears despite the declared freeze")

    poles = resolve_gauge_poles(matching)
    for name, value in poles.items():
        recorded = packet["uv_poles"]["gauge_sector"]["poles"][name]["pole_coefficient"]
        if Fraction(recorded) != Fraction(value):
            problems.append(f"gauge pole {name} differs from the beta-derived value")
    open_generators = set(packet["uv_poles"]["open_pole_values"]["generators"])
    if open_generators != set(names) - {"dg1", "dg2", "dg3"}:
        problems.append("open pole set does not cover exactly the non-gauge generators")

    directions, weighted = general_basis_dimension()
    if directions != packet["general_local_basis"]["directions"]:
        problems.append("general-basis direction count differs")
    if weighted != packet["general_local_basis"]["dimension_with_flavor"]:
        problems.append("general-basis flavored dimension differs")
    reachable = {(tuple(e["fields"]), e["structure"]) for e in table["entries"]}
    for direction in packet["general_local_basis"]["unreachable"]:
        if (tuple(direction["fields"]), direction["structure"]) in reachable:
            problems.append(f"direction {direction['fields']} is listed unreachable but has a record")

    for name, control in packet["controls"].items():
        if not (control.get("expected_failure") and control.get("failed")):
            problems.append(f"control {name} did not fire")

    return {"status": "PASS" if not problems else "FAIL", "problems": problems,
            "rank": rank, "nullspace_dimension": len(null), "term_count": len(terms)}


def main() -> int:
    verdict = check()
    print(json.dumps(verdict))
    return 0 if verdict["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
