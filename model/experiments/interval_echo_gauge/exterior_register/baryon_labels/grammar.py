"""Field records, multiset eligibility, and signed alternating contractions."""

from collections import Counter, defaultdict
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations_with_replacement, product
from pathlib import Path

from .. import material

# Reuse the declared contraction alphabet, pairings and rational rank algorithm.
ROOT = Path(__file__).resolve().parents[5]
spec = spec_from_file_location(
    "baryon_dimension_six_census",
    ROOT / "code/a5_closure/baryon_dimension_six_census.py",
)
formal = module_from_spec(spec)
spec.loader.exec_module(formal)

NAMES = formal.FIELD_ORDER
# This charge-to-name convention chooses which triplet is called Q and which
# fields are the adjoints of the left-Weyl even-degree rows.
CHARGES = (1, -3, 4, -2, -6, -1, 3, -4, 2, 6)


def field_census(packets, phases, side):
    fields = {}
    generators = defaultdict(list)
    for packet in packets:
        word, spin = (material.decode(s) for s in packet[side])
        if word in (0, 31):
            continue
        assert word is not None and spin in (0, 1)
        c, w = (word & 7).bit_count(), (word >> 3).bit_count()
        charge = -2 * c + 3 * w
        index = CHARGES.index(charge)
        triality = (0, 1, -1, 0)[c]
        left = int((c + w) % 2 == 0)
        entry = (
            charge,
            phases[word, "3B", side],
            phases[word, "L", side],
            triality,
            int(w == 1),
            left,
        )
        if index in fields:
            assert entry == fields[index]
        fields[index] = entry
        # The dual chart names a conjugate component by the complementary
        # positive word, including the dual weak index.
        positive = word if index < 5 else word ^ 31
        coloured = positive & 7
        colour = coloured.bit_length() - 1 if triality else 0
        weak = (positive >> 3).bit_length() - 1 if w == 1 else 0
        if triality:
            assert coloured.bit_count() == 1
        generators[index].append((2 * word + spin, colour, weak, spin))
    assert (
        len(fields) == 10
        and len({g[0] for group in generators.values() for g in group}) == 60
    )
    return fields, generators


def ballots(fields):
    records, admitted = [], []
    for row in combinations_with_replacement(range(10), 4):
        q, b, ell, colour, weak, left = map(
            sum, zip(*(fields[i] for i in row), strict=True)
        )
        tests = [b != 0, q == 0, colour % 3 == 0, weak % 2 == 0, left % 2 == 0]
        records.append(
            [list(row), [q, b, ell, colour, weak, left], [int(t) for t in tests]]
        )
        if all(tests):
            admitted.append(row)
    return records, admitted


def exterior_signs(trials, signed):
    signs = {}
    for order, setting, reported in trials:
        if setting == signed:
            result = material.decode(reported)
            assert result in (0, 16)
            signs[tuple(order)] = 1 if result == 0 else -1
    assert len(signs) == 24
    return signs


def contractions(row, fields, generators, signed, unsigned):
    coloured = [i for i, f in enumerate(row) if fields[f][3]]
    weak = [i for i, f in enumerate(row) if fields[f][4]]
    left = [i for i, f in enumerate(row) if fields[f][5]]
    right = [i for i in range(4) if i not in left]
    assert len(coloured) == 3
    signs = [signed, unsigned]
    outputs = []
    for weak_pairs, l_pairs, r_pairs in product(
        *(formal.perfect_pairings(slots) for slots in (weak, left, right))
    ):
        polys = [Counter(), Counter()]
        for components in product(*(generators[f] for f in row)):
            ids = tuple(c[0] for c in components)
            if len(set(ids)) != 4:
                continue
            factor = formal.eps3(*(components[i][1] for i in coloured))
            for a, b in weak_pairs:
                factor *= formal.eps2(components[a][2], components[b][2])
            for a, b in (*l_pairs, *r_pairs):
                factor *= formal.eps2(components[a][3], components[b][3])
            if factor:
                monomial = tuple(sorted(ids))
                order = tuple(monomial.index(n) for n in ids)
                for poly, sign in zip(polys, signs, strict=True):
                    poly[monomial] += sign[order] * factor
        outputs.append(
            dict(
                weak=weak_pairs,
                left=l_pairs,
                right=r_pairs,
                signed={k: v for k, v in polys[0].items() if v},
                unsigned={k: v for k, v in polys[1].items() if v},
            )
        )
    return outputs


def serialize_contractions(row, polys):
    def encode(poly):
        return [[*ids, n] for ids, n in sorted(poly.items())]

    return dict(
        row=list(row),
        signed_rank=formal.polynomial_rank([p["signed"] for p in polys]),
        unsigned_rank=formal.polynomial_rank([p["unsigned"] for p in polys]),
        pairings=[
            dict(
                weak=p["weak"],
                left=p["left"],
                right=p["right"],
                signed=encode(p["signed"]),
                unsigned=encode(p["unsigned"]),
            )
            for p in polys
        ],
    )
