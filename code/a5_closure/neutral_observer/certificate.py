"""Exact signed-current and resolved-charge records in the exterior dictionary."""

import argparse
import json
from pathlib import Path
import sys

from sympy import I, Matrix, Rational, SparseMatrix, diag, eye, kronecker_product, zeros

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.face_frame.field import (  # noqa: E402
    transfer,
)


def generator(i, j):
    raw = transfer(i, j)
    values = raw.astype(int)
    assert (values == raw).all()
    return SparseMatrix(values.tolist())


def exterior():
    n = [generator(i, i) for i in range(5)]
    d = sum(n[:3], zeros(32))
    y6 = -2 * d + 3 * (n[3] + n[4])
    up = generator(3, 4)
    q = (y6 + 3 * (n[3] - n[4])) / 6
    color = [generator(i, j) for i in range(3) for j in range(3) if i != j]
    color += [n[0] - n[1], n[1] - n[2]]
    return q, {
        "T+": up,
        "T-": up.T,
        "6Y": y6,
        "2T3": n[3] - n[4],
        **{f"color_{j}": g for j, g in enumerate(color)},
    }


def supported(matrix):
    return [i for i in range(matrix.rows) if matrix[i, i] != 0]


def selection(q, value):
    return diag(*(int(q[i, i] == value) for i in range(q.rows)))


def signed_effects(q, generators):
    p0 = selection(q, 0)
    effects = []
    for sign, shift in ((-1, "T-"), (1, "T+")):
        j = generators[shift]
        effect = p0 * j.T * selection(q, sign) * j * p0
        assert effect * effect == effect and effect.T == effect
        effects.append(effect)
    minus, plus = effects
    dark = p0 - minus - plus
    assert dark**2 == dark and minus * plus == zeros(32)
    assert all(dark * e == zeros(32) for e in effects)
    # On electromagnetic neutrality the normalized Z amplitude is 2*T3.
    z = p0 * generators["2T3"] * p0
    assert z == minus - plus and z**2 == minus + plus
    return tuple(effects) + (dark,)


def restrict(g, ids):
    # Check invariance before using a block as a closed representation.
    outside = [i for i in range(g.rows) if i not in ids]
    assert g.extract(outside, ids) == zeros(len(outside), len(ids))
    return g.extract(ids, ids)


def matches(generators, left, right, direct):
    for g in generators.values():
        a, b = restrict(g, left), restrict(g, right)
        tensor = kronecker_product(a, eye(len(right))) + kronecker_product(
            eye(len(left)), b
        )
        assert tensor == restrict(g, direct)


def bound_controls(q, generators):
    rows = {}
    cases = [
        ("negative_daughter", [7], [8, 16], [15, 23], 0),
        ("positive_daughter", [24], [15, 23], [8, 16], 1),
        ("weak_dark", [7], [24], [0], 0),
    ]
    for label, left, right, direct, initial in cases:
        matches(generators, left, right, direct)
        word = [left[0], right[initial]]
        charges = [int(q[i, i]) for i in word]
        assert sorted(charges) == [-1, 1] and q[direct[initial], direct[initial]] == 0
        rows[label] = dict(
            resolved_masks=word,
            constituent_charges=charges,
            direct_mask=direct[initial],
            left_charge_squared=charges[0] ** 2,
            global_generators=len(generators),
        )
    return rows


def frame_and_flight(generators):
    # This SU(2) frame is exact and mixes both weak slots.
    w = Matrix([[3, 4 * I], [4 * I, 3]]) / 5
    assert w.H * w == eye(2) and w.det() == 1
    j = restrict(generators["T-"], [15, 23])
    incoming = Matrix([1, 0])
    final = diag(0, 1)
    effect = j.H * final * j
    transported = (w * j * w.H).H * (w * final * w.H) * (w * j * w.H)
    assert transported == w * effect * w.H
    moved = w * incoming
    matched = (moved.H * transported * moved)[0]
    stale = (moved.H * effect * moved)[0]
    # A supplied two-station unitary moves either common center with the
    # same response. Independent internal operations commute with this flight.
    flight = kronecker_product(Matrix([[3, -4], [4, 3]]) / 5, eye(2))
    arrival = kronecker_product(diag(0, 1), eye(2))
    v = kronecker_product(Matrix([1, 0]), incoming)
    assert flight.H * flight == eye(4)
    for g in generators.values():
        local = kronecker_product(eye(2), restrict(g, [15, 23]))
        assert local * flight == flight * local
    at_receiver = ((flight * v).H * arrival * flight * v)[0]
    return dict(matched=str(matched), stale=str(stale), center_arrival=str(at_receiver))


def flavor_reads():
    # A calibrated, supplied unitary relates a neutral basis and two
    # identified charged mass pointers. An unresolved third pointer is dark
    # only to this restricted instrument.
    mixing = Matrix([[3, 4 * I, 0], [4 * I, 3, 0], [0, 0, 5]]) / 5
    assert mixing.H * mixing == eye(3)
    readers = []
    for i in range(3):
        p = zeros(3)
        p[i, i] = 1
        readers.append(mixing.H * p * mixing)
    assert sum(readers, zeros(3)) == eye(3)
    assert all(a * a == a for a in readers)
    assert all(readers[i] * readers[j] == zeros(3) for i in range(3) for j in range(i))
    tags = [e[0, 0] for e in readers]
    assert tags == [Rational(9, 25), Rational(16, 25), 0]
    return dict(
        first_basis_tags=[str(x) for x in tags],
        third_tag_visible_to_first_two=int((readers[0] + readers[1])[2, 2]),
    )


def mixture_records(pair_controls, effects):
    # Project the signed transition strengths in each copy of the neutral
    # representation. Keep the two sterile exterior lines distinct as inputs.
    direct = [0, 15, 16, 31]
    masks = [*direct, *(r["direct_mask"] for r in pair_controls.values())]
    spectator = diag(*([0] * len(direct)), *([1] * len(pair_controls)))
    effects = [diag(*(e[m, m] for m in masks)) for e in effects]
    joint = [site * e for site in (eye(7) - spectator, spectator) for e in effects]
    assert sum(joint, zeros(7)) == eye(7)
    assert all(e**2 == e for e in joint)
    assert all(a * b == zeros(7) for i, a in enumerate(joint) for b in joint[:i])
    # Half bare active / half dark pair vs half active pair / half bare dark.
    first = diag(0, 1, 0, 0, 0, 0, 1) / 2
    second = diag(1, 0, 0, 0, 1, 0, 0) / 2

    def rates(state, readers):
        return [str((state * e).trace()) for e in readers]

    separate = [*effects, spectator]
    assert rates(first, separate) == rates(second, separate)
    assert rates(first, joint) != rates(second, joint)
    return dict(
        marginal_minus_plus_dark_pair=rates(first, separate),
        first_joint=rates(first, joint),
        second_joint=rates(second, joint),
    )


def calibration():
    efficiency = Matrix([[6, 1], [1, 4]]) / 8
    active = Matrix([Rational(1, 2), Rational(1, 3)])
    arrival = Rational(2, 5)
    observed = arrival * efficiency * active
    failed = 1 - sum(observed)
    solved = efficiency.inv() * observed / arrival
    assert solved == active
    # Rates from unreleased traffic and entered sterile sources are identical
    # without an independent arrival instrument, even at known efficiency.
    active_no_dark = active / sum(active)
    other_arrival = arrival * sum(active)
    assert other_arrival * efficiency * active_no_dark == observed
    # Discarding final signs also confounds signed sectors at equal efficiencies.
    assert Matrix([[1, 1]]) * Matrix([1, 0]) == Matrix([[1, 1]]) * Matrix([0, 1])
    return dict(
        calibration=[[str(x) for x in efficiency.row(i)] for i in range(2)],
        determinant=str(efficiency.det()),
        arrival=str(arrival),
        per_launch=[str(x) for x in (*observed, failed)],
        inferred_entering_minus_plus_dark=[str(x) for x in (*solved, 1 - sum(solved))],
        matching_no_dark_arrival=str(other_arrival),
    )


def material_source():
    receipt = ROOT / (
        "model/experiments/interval_echo_angular/charged_orientation/"
        "three_arm_spin/cyclic_charge/neutral_ladder/results.json"
    )
    data = json.loads(receipt.read_text())
    words = data["ternary_word_order"]
    result = []
    for row in data["experiments"]:
        if row["setting"] != "full":
            continue
        word = row["neutral_word"]
        per_port = [n - 1 for n in word]
        assert per_port == [row["sign"], -row["sign"], 0]
        read = words.index(word)
        hits = 0
        starts = 0
        for _, _, bank in row["load_banks"]:
            assert len(bank) == len(words) + 2
            assert sum(bank) == data["starts_per_bank"]
            assert sum(bank[: len(words)]) == bank[read]
            hits += bank[read]
            starts += sum(bank)
        result.append(
            dict(sign=row["sign"], hub_target_donor=word, hits=hits, starts=starts)
        )
    assert sorted(x["sign"] for x in result) == [-1, 1]
    return result


def run():
    q, g = exterior()
    em, ep, ed = signed_effects(q, g)
    pairs = bound_controls(q, g)
    return dict(
        neutral_masks=supported(selection(q, 0)),
        signed_projectors={
            "negative_daughter": supported(em),
            "positive_daughter": supported(ep),
            "dark": supported(ed),
        },
        bound_pair_controls=pairs,
        mixture_controls=mixture_records(pairs, (em, ep, ed)),
        transported_reads=frame_and_flight(g),
        signed_calibration=calibration(),
        declared_family_control=flavor_reads(),
        measured_material_full_setting=material_source(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = HERE / "results.json"
    computed = run()
    if args.write:
        receipt.write_text(json.dumps(computed, indent=2) + "\n")
    elif computed != json.loads(receipt.read_text()):
        raise SystemExit("neutral observer receipt differs")
    print("Exact exterior, instrument and stock checks pass.")
