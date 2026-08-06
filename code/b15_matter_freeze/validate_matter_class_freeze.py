#!/usr/bin/env python3
"""Deterministic validator for the B15 matter class freeze (issue #706).

Replays the frozen enumeration exactly and checks every count, every
admissible-irrep list, and every provenance hash recorded in
matter_class_freeze_v1.json. Fails closed on any mismatch.

The enumeration consumes only the frozen structural data:
  - the dimension cap (equal-budget clause: the twelve-port carrier dimension);
  - the weight-lattice rules (abstract integral lattice of the certified
    bracket; carrier-generated sublattice as a flagged refinement);
  - the Z/2 chirality grading counted up to the global parity flip;
  - the exact cubic anomaly functional, recorded as an observable.

No target data enters: no Standard Model multiplet table, no laboratory
quantum number, no mass or mixing input. Run:

    python3 validate_matter_class_freeze.py            # validate the freeze
    python3 validate_matter_class_freeze.py --emit     # print computed counts
"""

import hashlib
import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
FREEZE_PATH = os.path.join(HERE, "matter_class_freeze_v1.json")


# ---------------------------------------------------------------------------
# Frozen irrep grammars (structural data only)
# ---------------------------------------------------------------------------

def a1_irreps(cap):
    """sl(2) irreps labeled by dimension n (highest weight n-1), dim <= cap.

    Cubic invariant vanishes identically on A1. Carrier-lattice flag: the
    source carrier acts through the adjoint, whose weights generate the root
    lattice, i.e. even highest weights, i.e. odd dimension.
    """
    items = []
    for n in range(1, cap + 1):
        items.append(("A1(%d)" % n, n, 0, n % 2 == 1))
    return items


def a2_dim(p, q):
    return (p + 1) * (q + 1) * (p + q + 2) // 2


def a2_anomaly(p, q):
    """Exact cubic anomaly A(p,q) = (p-q)(p+2q+3)(2p+q+3) dim(p,q) / 60.

    Normalized so A(1,0) = 1. Integer-valued on every irrep; the validator
    checks integrality.
    """
    val = Fraction((p - q) * (p + 2 * q + 3) * (2 * p + q + 3) * a2_dim(p, q), 60)
    if val.denominator != 1:
        raise AssertionError("cubic anomaly non-integer at (%d,%d)" % (p, q))
    return int(val)


def a2_irreps(cap):
    """sl(3) irreps by Dynkin label (p,q), dim <= cap.

    Carrier-lattice flag: root-lattice weights, i.e. p == q mod 3.
    """
    items = []
    labels = []
    for p in range(0, cap):
        for q in range(0, cap):
            d = a2_dim(p, q)
            if d <= cap:
                labels.append((p, q))
    labels.sort(key=lambda pq: (a2_dim(*pq), pq))
    for (p, q) in labels:
        d = a2_dim(p, q)
        items.append(("A2(%d,%d)" % (p, q), d, a2_anomaly(p, q), (p - q) % 3 == 0))
    return items


def product_items(left, right, cap):
    """Irreps of a direct sum: outer products with dim <= cap.

    Anomaly is additive along the A2 factor only; the A1 cubic invariant
    vanishes, so the product anomaly is (anomaly of left) * (trivial count)
    resolved explicitly below: for (R, S) the cubic anomaly is
    A(R) * dim(S) when the cubic Casimir lives on the left factor and the
    right factor is anomaly-free, and symmetrically. Here every product in
    the frozen class has at most one anomaly-bearing factor.
    """
    items = []
    for (ll, ld, la, lc) in left:
        for (rl, rd, ra, rc) in right:
            d = ld * rd
            if d > cap:
                continue
            if la != 0 and ra != 0:
                raise AssertionError("two anomaly-bearing factors")
            anom = la * rd + ra * ld
            items.append(("%s*%s" % (ll, rl), d, anom, lc and rc))
    items.sort(key=lambda t: (t[1], t[0]))
    return items


def abelian_items(cap):
    """Abelian stratum: the certified charge lattice is {0} (constraint C3),
    so the trivial one-dimensional module is the only admissible irrep."""
    del cap
    return [("triv", 1, 0, True)]


def items_for_type(type_key, cap):
    if type_key == "abelian":
        return abelian_items(cap)
    if type_key == "A1":
        return a1_irreps(cap)
    if type_key == "A1A1":
        return product_items(a1_irreps(cap), a1_irreps(cap), cap)
    if type_key == "A2":
        return a2_irreps(cap)
    if type_key == "A2A1":
        return product_items(a2_irreps(cap), a1_irreps(cap), cap)
    raise AssertionError("unknown type key %r" % type_key)


# ---------------------------------------------------------------------------
# Enumeration (multisets of irreps under the budget)
# ---------------------------------------------------------------------------

def multisets(items, cap):
    """All multisets of items with total dimension <= cap, as tuples of
    multiplicities aligned with the item list. Deterministic recursion."""
    n = len(items)
    out = []

    def rec(idx, budget, acc):
        if idx == n:
            out.append(tuple(acc))
            return
        dim = items[idx][1]
        top = budget // dim
        for m in range(0, top + 1):
            acc.append(m)
            rec(idx + 1, budget - m * dim, acc)
            acc.pop()

    rec(0, cap, [])
    return out


def anomaly_zero_gradings(mults, anomalies):
    """Number of gradings (per-irrep even-copy counts e_i in 0..m_i) with
    sum_i (2 e_i - m_i) * A_i == 0, by exact dictionary DP."""
    table = {0: 1}
    for m, a in zip(mults, anomalies):
        if m == 0:
            continue
        nxt = {}
        for s, c in table.items():
            for e in range(0, m + 1):
                s2 = s + (2 * e - m) * a
                nxt[s2] = nxt.get(s2, 0) + c
        table = nxt
    return table.get(0, 0)


def count_class(items, cap):
    """Frozen counts for one bracket type.

    ungraded            nonempty multisets of admissible irreps, dim <= cap
    graded_up_to_flip   nonempty Z/2-graded such multisets, identified under
                        the global parity flip (Burnside per multiset)
    graded_anomaly_zero graded candidates (up to flip) whose signed cubic
                        anomaly sum vanishes (recorded observable, no
                        membership constraint)
    """
    anomalies = [it[2] for it in items]
    ungraded = 0
    graded = 0
    graded_zero = 0
    for mults in multisets(items, cap):
        if all(m == 0 for m in mults):
            continue
        ungraded += 1
        total_gradings = 1
        for m in mults:
            total_gradings *= (m + 1)
        fixed = 1 if all(m % 2 == 0 for m in mults) else 0
        assert (total_gradings + fixed) % 2 == 0
        graded += (total_gradings + fixed) // 2
        n0 = anomaly_zero_gradings(mults, anomalies)
        assert (n0 + fixed) % 2 == 0
        graded_zero += (n0 + fixed) // 2
    return {
        "ungraded": ungraded,
        "graded_up_to_flip": graded,
        "graded_anomaly_zero_up_to_flip": graded_zero,
    }


def counts_for_type(type_key, cap):
    items = items_for_type(type_key, cap)
    full = count_class(items, cap)
    carrier_items = [it for it in items if it[3]]
    carrier = count_class(carrier_items, cap)
    return items, full, carrier


TYPE_KEYS = ["abelian", "A1", "A1A1", "A2", "A2A1"]


def compute_all(cap):
    result = {}
    for tk in TYPE_KEYS:
        items, full, carrier = counts_for_type(tk, cap)
        result[tk] = {
            "admissible_irreps": [
                {
                    "label": lab,
                    "dim": dim,
                    "cubic_anomaly": anom,
                    "carrier_lattice": flag,
                }
                for (lab, dim, anom, flag) in items
            ],
            "admissible_irrep_count": len(items),
            "carrier_lattice_irrep_count": len([1 for it in items if it[3]]),
            "counts": full,
            "carrier_lattice_subclass_counts": carrier,
        }
    return result


# ---------------------------------------------------------------------------
# Freeze document checks
# ---------------------------------------------------------------------------

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def fail(msg):
    print("FREEZE INVALID: %s" % msg)
    sys.exit(1)


def validate(freeze):
    errors = []

    cap = freeze["budget"]["dimension_cap"]
    if cap != 12:
        errors.append("dimension cap %r differs from the frozen value 12" % cap)

    computed = compute_all(cap)

    # 1. Enumeration tables match byte-exactly.
    frozen_tables = freeze["enumeration"]["per_bracket_type"]
    if sorted(frozen_tables.keys()) != sorted(TYPE_KEYS):
        errors.append("bracket type keys differ: %s" % sorted(frozen_tables.keys()))
    for tk in TYPE_KEYS:
        if tk not in frozen_tables:
            continue
        if frozen_tables[tk] != computed[tk]:
            errors.append("enumeration mismatch for bracket type %s" % tk)

    # 2. Strata map onto computed types and mirror symmetry holds.
    strata = freeze["enumeration"]["strata"]
    for s in strata:
        if s["bracket_type"] not in TYPE_KEYS:
            errors.append("stratum %s names unknown type %s" % (s["id"], s["bracket_type"]))
    f2 = [s for s in strata if s["id"] == "F2"]
    g2 = [s for s in strata if s["id"] == "G2"]
    if not (f2 and g2 and f2[0]["bracket_type"] == g2[0]["bracket_type"]):
        errors.append("mirror strata F2/G2 do not share a bracket type")

    # 3. Provenance hashes.
    for pin in freeze["provenance"]["pinned_files"]:
        path = os.path.join(REPO, pin["path"])
        if not os.path.exists(path):
            errors.append("pinned file missing: %s" % pin["path"])
            continue
        got = sha256_file(path)
        if got != pin["sha256"]:
            errors.append("hash mismatch for %s: %s" % (pin["path"], got))

    # 4. Discipline flags.
    rule = freeze["decision_rule"]
    if rule.get("no_selection_is_valid_negative_exit") is not True:
        errors.append("negative exit declaration missing")
    if freeze.get("selection_performed") is not False:
        errors.append("freeze document must record selection_performed = false")
    barred = set(freeze["observables"]["barred_selection_inputs"])
    for required in [
        "lepton masses",
        "quark masses",
        "mixing angles",
        "laboratory quantum numbers",
        "Standard Model multiplet tables",
    ]:
        if required not in barred:
            errors.append("barred list lacks %r" % required)

    # 5. Anomaly formula spot checks (frozen normalization).
    checks = [((1, 0), 1), ((2, 0), 7), ((1, 1), 0), ((3, 0), 27), ((2, 1), 14)]
    for (p, q), want in checks:
        if a2_anomaly(p, q) != want:
            errors.append("anomaly formula check failed at (%d,%d)" % (p, q))

    if errors:
        for e in errors:
            print("FREEZE INVALID: %s" % e)
        sys.exit(1)

    print("FREEZE VALID: enumeration, provenance, and discipline checks pass")
    for tk in TYPE_KEYS:
        c = computed[tk]["counts"]
        cc = computed[tk]["carrier_lattice_subclass_counts"]
        print(
            "  %-7s irreps=%2d  ungraded=%5d graded=%6d anomaly0=%6d  "
            "[carrier: irreps=%2d ungraded=%4d graded=%5d anomaly0=%5d]"
            % (
                tk,
                computed[tk]["admissible_irrep_count"],
                c["ungraded"],
                c["graded_up_to_flip"],
                c["graded_anomaly_zero_up_to_flip"],
                computed[tk]["carrier_lattice_irrep_count"],
                cc["ungraded"],
                cc["graded_up_to_flip"],
                cc["graded_anomaly_zero_up_to_flip"],
            )
        )
    sys.exit(0)


def main():
    if "--emit" in sys.argv:
        print(json.dumps(compute_all(12), indent=1, sort_keys=True))
        return
    if not os.path.exists(FREEZE_PATH):
        fail("freeze document missing at %s" % FREEZE_PATH)
    with open(FREEZE_PATH) as fh:
        freeze = json.load(fh)
    validate(freeze)


if __name__ == "__main__":
    main()
