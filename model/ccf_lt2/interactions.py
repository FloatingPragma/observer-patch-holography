#!/usr/bin/env python3
"""Full local CAR interaction additions for the declared CCF-2 source.

No decay projector or measured decay rate occurs in this Hamiltonian.
All coefficients are exact; the parent semileptonic contact is untouched.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as sy

from .core import (add, adjoint, clean, index, letter, mono, multiply,
                  native_charge_coordinates, normal, number, scale, substitute)

HERE = Path(__file__).resolve().parent

GL = sy.Rational(1, 4096)
GQ = sy.Rational(1, 4096)
UP_SHIFTS = (sy.Integer(0), sy.Integer(2), sy.Integer(8))
DOWN_SHIFTS = (sy.Integer(0), sy.Integer(1), sy.Integer(4))
FAMILY_ROTATION = sy.ImmutableMatrix([[sy.Rational(4, 5), sy.Rational(3, 5), 0],
                                    [-sy.Rational(3, 5), sy.Rational(4, 5), 0],
                                    [0, 0, 1]])
UP = (2, 5, 7)
DOWN = (3, 6, 8)
ANTI_UP = (4, 1, 0)
ANTI_DOWN = (14, 13, 12)
COLOUR_JOIN = (1, -1, 1)
LEPTON = (10, 11)
ANTI_LEPTON = (15, 9)
COLOURFUL = frozenset(UP + DOWN + ANTI_UP + ANTI_DOWN)
CHARGE_THIRDS = (-2, -2, 2, -1, -2, 2, -1, 2, -1, 3, 0, -3, 1, 1, 1, 0)
PAULI = (sy.ImmutableMatrix([[0, 1], [1, 0]]),
         sy.ImmutableMatrix([[0, -sy.I], [sy.I, 0]]),
         sy.ImmutableMatrix([[1, 0], [0, -1]]))
BAR_SIGMA = (sy.eye(2),) + tuple(-p for p in PAULI)
WEAK = tuple(p / 2 for p in PAULI)
ETA = (1, -1, -1, -1)


def _sum(terms):
    out = defaultdict(lambda: sy.Integer(0))
    for factor, poly in terms:
        if factor:
            for word, value in poly.items():
                out[word] += factor * value
    return clean(out)


def _lepton_field(v, s, weak, family, canonical=False):
    """Native left c, or its exact expression in canonical flavour modes.

    Neutral flavour modes are not asserted to be neutral mass eigenmodes.
    """
    mode = index(v, s, LEPTON[weak], family)
    if not canonical:
        return letter(mode)
    partner = index(v, 1-s, ANTI_LEPTON[weak], family)
    return scale(add(letter(mode), scale(letter(partner, True), -(1-2*s))),
                 1/sy.sqrt(2))


def _quark_field(v, s, weak, colour, family):
    """Full normalized left Nambu field on native quarks AND conjugates."""
    quark = (UP, DOWN)[weak][colour]
    antiquark = (ANTI_UP, ANTI_DOWN)[weak][colour]
    sign = (1-2*s)*COLOUR_JOIN[colour]
    return scale(add(letter(index(v, s, quark, family)),
                     scale(letter(index(v, 1-s, antiquark, family), True), sign)),
                 1/sy.sqrt(2))


def _quark_complement(v, s, weak, colour, family):
    quark = (UP, DOWN)[weak][colour]
    antiquark = (ANTI_UP, ANTI_DOWN)[weak][colour]
    sign = (1-2*s)*COLOUR_JOIN[colour]
    return scale(add(letter(index(v, 1-s, antiquark, family)),
                     scale(letter(index(v, s, quark, family), True), -sign)),
                 1/sy.sqrt(2))


def _currents(fields):
    """All twelve Hermitian currents, no positive-frequency restriction.

    fields[(spin, weak, channel)] denotes a CAR-linear left field; channels
    are the retained family labels and, for quarks, colour labels.
    """
    channels = sorted({key[2] for key in fields})
    currents = []
    for mu in range(4):
        for axis in range(3):
            terms = []
            for channel in channels:
                for s in range(2):
                    for z in range(2):
                        for w in range(2):
                            for u in range(2):
                                coefficient = BAR_SIGMA[mu][s, z]*WEAK[axis][w, u]
                                if coefficient:
                                    terms.append((coefficient, multiply(
                                        adjoint(fields[s, w, channel]), fields[z, u, channel])))
            currents.append(_sum(terms))
    return tuple(currents)


@lru_cache(maxsize=64)
def lepton_currents(v=0, canonical=False):
    fields = {(s, w, f): _lepton_field(v, s, w, f, canonical)
              for s in range(2) for w in range(2) for f in range(3)}
    return _currents(fields)


@lru_cache(maxsize=64)
def quark_currents(v=0, rotation=FAMILY_ROTATION):
    rotation = sy.ImmutableMatrix(rotation)
    if rotation.shape != (3, 3) or rotation*rotation.H != sy.eye(3):
        raise ValueError("Quark family frame must be exactly unitary")
    fields = {}
    for s in range(2):
        for colour in range(3):
            for f in range(3):
                fields[s, 0, (colour, f)] = _quark_field(v, s, 0, colour, f)
                fields[s, 1, (colour, f)] = _sum(
                    (rotation[f, g], _quark_field(v, s, 1, colour, g)) for g in range(3))
    return _currents(fields)


def _contact(currents, coupling):
    if coupling == 0:
        return {}
    return _sum((2*coupling*ETA[i//3], multiply(current, current))
                for i, current in enumerate(currents))


@lru_cache(maxsize=64)
def lepton_contact(v=0, coupling=GL, canonical=False):
    return _contact(lepton_currents(v, canonical), sy.sympify(coupling))


@lru_cache(maxsize=64)
def quark_contact(v=0, coupling=GQ, rotation=FAMILY_ROTATION):
    coupling = sy.sympify(coupling)
    rotation = sy.ImmutableMatrix(rotation)
    # Translation preserves CAR order because every mode gets one offset.
    # Assemble a local current square once, then retain its exact full word map.
    if v:
        return {tuple((i+96*v, creation) for i, creation in word): coefficient
                for word, coefficient in quark_contact(0, coupling, rotation).items()}
    return _contact(quark_currents(0, rotation), coupling)


def family_splitting(v=0, up_shifts=UP_SHIFTS, down_shifts=DOWN_SHIFTS):
    if len(up_shifts) != 3 or len(down_shifts) != 3:
        raise ValueError("All three family shifts are required")
    return _sum((sy.sympify(shifts[f]), number(index(v, s, coordinate, f)))
                for coordinates, shifts in ((UP+ANTI_UP, up_shifts),
                                            (DOWN+ANTI_DOWN, down_shifts))
                for coordinate in coordinates for s in range(2) for f in range(3))


def interaction_components(nvertices=2, *, g_leptonic=GL, g_nonleptonic=GQ,
                           up_shifts=UP_SHIFTS, down_shifts=DOWN_SHIFTS,
                           family_rotation=FAMILY_ROTATION):
    """Return complete ungauged native additions, local at every vertex."""
    if not isinstance(nvertices, int) or nvertices < 1:
        raise ValueError("nvertices must be a positive integer")
    return {
        "LL": add(*(lepton_contact(v, sy.sympify(g_leptonic)) for v in range(nvertices))),
        "quark_family": add(*(family_splitting(v, tuple(up_shifts), tuple(down_shifts))
                              for v in range(nvertices))),
        "QQ_nonleptonic": add(*(quark_contact(v, sy.sympify(g_nonleptonic),
                                             sy.ImmutableMatrix(family_rotation))
                                for v in range(nvertices))),
    }


def _word_change(word, weights):
    return sum((1 if creation else -1)*weights(mode) for mode, creation in word)


def _coordinate(mode):
    return (mode % 48)//3


def _baryon_thirds(mode):
    c = _coordinate(mode)
    return 1 if c in UP+DOWN else -1 if c in ANTI_UP+ANTI_DOWN else 0


def _lepton_number(mode):
    c = _coordinate(mode)
    return 1 if c in LEPTON else -1 if c in ANTI_LEPTON else 0


def audit_polynomial(poly, nvertices, charge_weight=None):
    charge_weight = charge_weight or (lambda m: CHARGE_THIRDS[_coordinate(m)])
    if adjoint(poly) != poly:
        raise AssertionError("Non-Hermitian interaction")
    counts = defaultdict(int)
    number_changing = 0
    family_changing = 0
    for word, coefficient in poly.items():
        if normal(word) != ((word, 1),) or coefficient == 0:
            raise AssertionError("Polynomial is not normalized")
        sites = {mode//96 for mode, _ in word}
        if len(sites) > 1 or any(site >= nvertices for site in sites):
            raise AssertionError("A local interaction has cross-site support")
        if _word_change(word, charge_weight):
            raise AssertionError("Electric charge is not preserved")
        if _word_change(word, _baryon_thirds):
            raise AssertionError("Baryon number is not preserved")
        if _word_change(word, _lepton_number):
            raise AssertionError("Dirac lepton number is not preserved")
        colourful = [(_coordinate(m) in COLOURFUL) for m, _ in word]
        if any(colourful) and not all(colourful):
            raise AssertionError("Module silently adds a semileptonic term")
        if all(colourful) and word:
            number_changing += bool(_word_change(word, lambda m: 1))
            family_changing += any(_word_change(word, lambda m, f=f: int(m % 3 == f))
                                   for f in range(3))
        counts[len(word)] += 1
    return {"monomials": len(poly), "degree_counts": dict(sorted(counts.items())),
            "colourful_number_changing_monomials": number_changing,
            "colourful_family_changing_monomials": family_changing,
            "hermiticity_charge_baryon_lepton_locality_exact": True,
            "no_added_semileptonic_monomials": True}


def act(poly, occupation):
    """Exact complete CAR action on one occupation; no target filtering."""
    result = defaultdict(lambda: sy.Integer(0))
    for word, coefficient in poly.items():
        state, sign = occupation, 1
        for mode, creation in reversed(word):
            if bool(state & (1 << mode)) == creation:
                sign = 0
                break
            sign *= -1 if (state & ((1 << mode)-1)).bit_count() % 2 else 1
            state ^= 1 << mode
        if sign:
            result[state] += sign*coefficient
    return {state: value for state, raw in result.items() if (value := sy.expand(raw)) != 0}


def _flavour_map(v=0):
    mapping = {}
    for s in range(2):
        for w in range(2):
            for f in range(3):
                i = index(v, s, LEPTON[w], f)
                j = index(v, 1-s, ANTI_LEPTON[w], f)
                sigma = 1-2*s
                mapping[i] = scale(add(letter(i), scale(letter(j, True), -sigma)), 1/sy.sqrt(2))
                mapping[j] = scale(add(letter(j), scale(letter(i, True), sigma)), 1/sy.sqrt(2))
    return mapping


def _car_gate(fields):
    for i, left in enumerate(fields):
        for j, right in enumerate(fields):
            if add(multiply(left, right), multiply(right, left)):
                raise AssertionError("Annihilator CAR failed")
            anti = add(multiply(left, adjoint(right)), multiply(adjoint(right), left))
            if anti != (mono() if i == j else {}):
                raise AssertionError("Mixed CAR failed")
    return len(fields)**2


def _colour_gate():
    """Read exact frozen exterior generators, checking their native basis."""
    artifact = HERE / "representation_snapshot.json"
    packet = json.loads(artifact.read_text())
    subsets = [tuple(row) for row in packet["native_basis_subsets"]]
    expected = [(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),
                (3,4),(0,1,2,3),(0,1,2,4),(0,1,3,4),(0,2,3,4),(1,2,3,4)]
    if subsets != expected:
        raise AssertionError("Native exterior order changed")
    generators = [sy.Matrix(row) for row in packet["colour_generators"]]
    M = sy.diag(*COLOUR_JOIN)
    checks = 0
    for lifted in generators[:8]:
        for q, a in ((UP, ANTI_UP), (DOWN, ANTI_DOWN)):
            if lifted.extract(q,q).T*M + M*lifted.extract(a,a) != sy.zeros(3):
                raise AssertionError("Quark/antiquark colour pairing failed")
            checks += 1
    return {"basis_entries": len(subsets), "exact_colour_intertwiner_checks": checks,
            "artifact_sha256": hashlib.sha256(artifact.read_bytes()).hexdigest()}


def checks():
    assert native_charge_coordinates() == CHARGE_THIRDS
    components = interaction_components(2)
    audit = {name: audit_polynomial(poly, 2) for name, poly in components.items()}
    assert audit["QQ_nonleptonic"]["colourful_number_changing_monomials"] > 0
    assert audit["QQ_nonleptonic"]["colourful_family_changing_monomials"] > 0
    assert _colour_gate()["exact_colour_intertwiner_checks"] == 16
    car_quark = _car_gate([factory(0, s, w, colour, f)
                           for factory in (_quark_field, _quark_complement)
                           for s in range(2) for w in range(2)
                           for colour in range(3) for f in range(3)])
    car_lepton = _car_gate(list(_flavour_map().values()))
    # Compare two independently composed expressions of the complete LL source.
    canonical = lepton_contact(0, GL, True)
    if substitute(lepton_contact(0), _flavour_map()) != canonical:
        raise AssertionError("Full native/canonical lepton polynomial differs")
    muon = 1 << index(0, 0, 11, 1)
    outgoing = act(canonical, muon)
    mu_channels = {}
    for es in range(2):
        for ns in range(2):
            for anti_s in range(2):
                mask = ((1 << index(0, es, 11, 0)) |
                        (1 << index(0, ns, 10, 1)) |
                        (1 << index(0, anti_s, 15, 0)))
                if mask in outgoing:
                    mu_channels[str(mask)] = str(outgoing[mask])
    if not mu_channels:
        raise AssertionError("No muon/electron/neutral-pair canonical vertex")
    assert not act(lepton_contact(0, sy.Integer(0), True), muon)
    # Exact canonical hadron occupation diagnostic: uds -> uud + d anti-u.
    initial = sum(1 << index(0, s, c, f) for s,c,f in
                  ((0, UP[0], 0), (1, DOWN[1], 0), (0, DOWN[2], 1)))
    qout = act(quark_contact(), initial)
    strange_channels = {}
    for mask, amplitude in qout.items():
        occupied = [mode for mode in range(96) if mask & (1 << mode)]
        if len(occupied) != 5 or any(mode % 3 != 0 for mode in occupied):
            continue
        species = [_coordinate(mode) for mode in occupied]
        if (sum(c in UP for c in species) == 2 and sum(c in DOWN for c in species) == 2
                and sum(c in ANTI_UP for c in species) == 1):
            strange_channels[str(mask)] = str(amplitude)
    if not strange_channels:
        raise AssertionError("No strange-to-light plus pair vertex")
    zero_rotation = act(quark_contact(0, GQ, sy.ImmutableMatrix(sy.eye(3))), initial)
    for mask in strange_channels:
        assert zero_rotation.get(int(mask), 0) == 0
    assert not quark_contact(0, sy.Integer(0))
    # Symmetry-relevant and integrity controls are reported separately.
    bad = dict(components["LL"])
    nonself = next(w for w in bad if adjoint({w: bad[w]}) != {w: bad[w]})
    bad.pop(nonself)
    try:
        audit_polynomial(bad, 2)
    except AssertionError:
        missing_adjoint_rejected = True
    else:
        raise AssertionError("Hermiticity control did not fail")
    try:
        audit_polynomial(components["LL"], 2, lambda m: CHARGE_THIRDS[_coordinate(m)]
                         + int(m == index(0,0,11,0)))
    except AssertionError:
        bad_charge_rejected = True
    else:
        raise AssertionError("Changed-charge control did not fail")
    # Native colourful number-conserving truncation kills every N3->N5 channel.
    truncation = {w:c for w,c in quark_contact().items()
                  if _word_change(w, lambda m: int(_coordinate(m) in COLOURFUL)) == 0}
    tout = act(truncation, initial)
    assert all(tout.get(int(mask), 0) == 0 for mask in strange_channels)
    canonical_truncated = {w:c for w,c in canonical.items()
                           if _word_change(w, lambda m: 1) == 0}
    lout = act(canonical_truncated, muon)
    assert all(lout.get(int(mask), 0) == 0 for mask in mu_channels)
    zero = interaction_components(2, g_leptonic=0, g_nonleptonic=0,
                                  up_shifts=(0,0,0), down_shifts=(0,0,0))
    assert not any(zero.values())
    return {
        "construction": "CCF-2 local leptonic/flavour/nonleptonic delta",
        "parameters": {"g_leptonic": str(GL), "g_nonleptonic": str(GQ),
                       "up_family_shifts": list(map(str, UP_SHIFTS)),
                       "down_family_shifts": list(map(str, DOWN_SHIFTS)),
                       "rotation": [[str(z) for z in row] for row in FAMILY_ROTATION.tolist()]},
        "full_two_site_native_polynomial": audit,
        "exact_colour_intertwining": _colour_gate(),
        "CAR_pair_checks": {"full_quark_Nambu_map": car_quark, "full_lepton_flavour_map": car_lepton},
        "complete_native_canonical_LL_identity": True,
        "muon_input_occupation": str(muon), "muon_full_action_outcomes": len(outgoing),
        "muon_to_e_neutral_pair_amplitudes": mu_channels,
        "strange_input_occupation": str(initial), "strange_full_action_outcomes": len(qout),
        "strange_to_light_plus_pair_amplitudes": strange_channels,
        "controls": {"zero_leptonic_contact_blocks_vertex": True,
                     "zero_nonleptonic_contact_blocks_vertex": True,
                     "identity_family_rotation_blocks_strange_channel": True,
                     "number_conserving_truncation_blocks_required_vertices": True,
                     "missing_adjoint_rejected": missing_adjoint_rejected,
                     "changed_charge_label_rejected": bad_charge_rejected,
                     "all_zero_deltas_return_zero": True},
        "uniform_extensive_lower_bound": {
            "LL_per_vertex": str(-216*abs(GL)),
            "QQ_per_vertex": str(-1944*abs(GQ)),
            "family_splitting_positive": True,
            "reason": "Each of 12 Hermitian weak currents has norm at most 3 (leptons) "
                      "or 9 (quarks). Triangle bounds on all current squares are "
                      "independent of total vertex count."},
        "scope": "Exact full-polynomial selection-rule and vertex checks, not decay probabilities. "
                 "Occupation probes are not Lambda/proton/pion eigenstates. Full interacting vacuum, "
                 "energy-conserving outgoing shells, physical masses and lifetimes are not established. "
                 "Parent semileptonic terms are retained unchanged by the source assembler.",
    }

