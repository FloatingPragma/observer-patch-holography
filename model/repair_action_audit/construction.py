"""A repair action needs the instrument, not just the classical endpoint rule.

All identities are exact in Q(sqrt(2))(i). The continuous parameter used for
generators is an explicitly optional unit-rate Poisson attempt parameter per
seam. It is neither the discrete attempt counter nor laboratory time.
"""
from __future__ import annotations

from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path

import numpy as np

from ..observer import native
from ..record_laws import family as f
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

S = f.S
I = S(0, 0, 1)
PHASES = {"0": S(1), "pi/2": I, "pi": S(-1), "3pi/2": -I}
ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = (
    "model/repair_action_audit/construction.py",
    "model/repair_action_audit/__init__.py",
    "model/repair_action_audit/__main__.py",
    "model/record_laws/family.py",
    "model/observer/native.py",
    "model/predictions/observer_instrument.py",
    "tools/strict_json.py",
)


def imaginary(value):
    """The real coefficient of i, retaining exact surds."""
    value = S.of(value)
    return S(value.c, value.d)


def unit_phase(value):
    value = S.of(value)
    if value*value.conj() != 1:
        raise ValueError("an exact unit-modulus phase is required")
    return value


def phase_instrument(phase):
    """Actual 78-state C instrument with a phase on waiting single loads.

    The original coin/current labels and every diagonal selective map remain
    unchanged. This is stronger than merely preserving the unread channel's
    classical marginal.
    """
    phase = unit_phase(phase)
    return f.coin_instrument(
        lambda key, coin, current: ("coin_current", coin, current),
        lambda key, coin: phase if sum(key) == 1 and f.current(key, coin) == 0 else 1,
    )


def phase_table(phase):
    phase = unit_phase(phase)
    return f.table(1, phase.conj()*F(1, 2), phase.conj()*F(1, 2), 0, F(1, 2))


def binary_channel(table):
    """Restrict an admissible table to the invariant binary two-port algebra."""
    if not f.positive_semidefinite(table):
        raise ValueError("an admissible complete record-law table is required")
    channel = np.full((16, 16), S(), dtype=object)
    multipliers = f.multipliers(table)
    for left, right in product(range(4), repeat=2):
        first, second = (left & 1, left >> 1), (right & 1, right >> 1)
        for current in f.CURRENTS:
            coefficient = multipliers.get((current, first, second), S())
            if coefficient:
                a = (first[0]-current)+2*(first[1]+current)
                b = (second[0]-current)+2*(second[1]+current)
                channel[4*a+b, 4*left+right] += coefficient
    return channel


def adjoint(matrix):
    return np.array([[S.of(value).conj() for value in row] for row in matrix.T], dtype=object)


def identity(size):
    return np.array([[S(i == j) for j in range(size)] for i in range(size)], dtype=object)


def native_binary_groups(phase, seam=(0, 1), ports=2):
    """Independent source-native Kraus construction on a small binary carrier."""
    phase = unit_phase(phase)
    if (type(ports) is not int or not 2 <= ports <= 3 or type(seam) is not tuple
            or len(seam) != 2 or any(type(p) is not int for p in seam)
            or not 0 <= seam[0] < seam[1] < ports):
        raise ValueError("an ordered seam on two or three binary ports is required")
    groups = {}
    for mask in range(1 << ports):
        loads = tuple((mask >> p) & 1 for p in range(ports))
        for coin in (0, 1):
            output, (_, current) = native.native(loads, seam, coin)
            target = sum(value << p for p, value in enumerate(output))
            waiting = loads[seam[0]]+loads[seam[1]] == 1 and current == 0
            groups.setdefault((coin, current), {})[mask] = (
                target, f.ROOT_HALF*(phase if waiting else 1))
    return groups


def group_channel(groups, dimension):
    channel = np.full((dimension**2, dimension**2), S(), dtype=object)
    for branch in groups.values():
        for left, right in product(branch, repeat=2):
            (a, u), (b, v) = branch[left], branch[right]
            channel[dimension*a+b, dimension*left+right] += u*v.conj()
    return channel


def spectrum(table):
    """The complete multiset in a declared orthogonal eigenbasis."""
    return ([S(1)]*3+[S(0)] + [table["k"], table["k"].conj()]
            + [table["g"]]*2+[table["g"].conj()]*2
            + [table["h"]]*2+[table["h"].conj()]*2+[table["m"]]*2)


def phase_cycle_defect(table):
    """Zero exactly when the Hilbert--Schmidt skew part is a commutator."""
    return imaginary(table["k"])+imaginary(table["h"])-imaginary(table["g"])


def commutator_generator(energies):
    """Matrix of -i[diag(energies), X] in row-major matrix-unit coordinates."""
    dim = len(energies)
    diagonal = [-I*(energies[a]-energies[b]) for a, b in product(range(dim), repeat=2)]
    result = np.full((dim*dim, dim*dim), S(), dtype=object)
    for i, value in enumerate(diagonal):
        result[i, i] = value
    return result


def local_certificate(table):
    channel = binary_channel(table)
    dagger = adjoint(channel)
    symmetric_action = identity(16)-(channel+dagger)*F(1, 2)
    skew = (channel-dagger)*F(1, 2)
    # Every off-diagonal matrix unit is an eigenvector. The only nondiagonal
    # block is the classical averaging of |10><10| and |01><01|.
    diagonal_indices = (0, 5, 10, 15)
    observed = []
    for index in range(16):
        if index in diagonal_indices:
            continue
        if any(channel[row, index] for row in range(16) if row != index):
            raise AssertionError("unexpected mixing of binary coherences")
        observed.append(channel[index, index])
    actual = [S(1)]*3+[S(0)]+observed
    if sorted(map(str, actual)) != sorted(map(str, spectrum(table))):
        raise AssertionError("binary channel spectrum differs")
    expected_population = np.array([[S(1), S(), S(), S()],
        [S(), S(F(1, 2)), S(F(1, 2)), S()],
        [S(), S(F(1, 2)), S(F(1, 2)), S()],
        [S(), S(), S(), S(1)]], dtype=object)
    if not np.array_equal(channel[np.ix_(diagonal_indices, diagonal_indices)], expected_population):
        raise AssertionError("classical population block changed")
    nonnegative = all((1-value.real()).sign() >= 0 for value in spectrum(table))
    if not nonnegative:
        raise AssertionError("negative symmetric action eigenvalue")
    selfadjoint = np.array_equal(channel, dagger)
    if selfadjoint != all(table[key].is_real() for key in ("k", "g", "h")):
        raise AssertionError("self-adjointness criterion failed")
    energies = (S(), imaginary(table["g"]), imaginary(table["g"]), imaginary(table["k"]))
    hamiltonian = np.array_equal(skew, commutator_generator(energies))
    if hamiltonian != (phase_cycle_defect(table) == 0):
        raise AssertionError("commutator criterion failed")
    # Rank of this explicit block decomposition gives its fixed dimension.
    fixed_dimension = 3+sum(value == 1 for value in observed)
    return {
        "binary_channel_spectrum": sorted(map(str, actual)),
        "self_adjoint_Liouville_action": selfadjoint,
        "symmetric_action_nonnegative": nonnegative,
        "fixed_operator_dimension": fixed_dimension,
        "skew_is_Hamiltonian_commutator": hamiltonian,
        "phase_cycle_defect": str(phase_cycle_defect(table)),
        "commutator_diagonal_if_exists": list(map(str, energies)) if hamiltonian else None,
    }, symmetric_action, skew


def source_selective_equality(phase):
    """Check every original input and actual coin/current branch, including failure."""
    base, changed = f.INSTRUMENTS["C"](), phase_instrument(phase)
    if set(base) != set(changed) or not f.completeness(changed):
        raise AssertionError("instrument labels or completeness changed")
    comparisons = 0
    for label in base:
        if set(base[label]) != set(changed[label]):
            raise AssertionError("branch support changed")
        for state in base[label]:
            (a, u), (b, v) = base[label][state], changed[label][state]
            if a != b or u*u.conj() != v*v.conj():
                raise AssertionError("a selective classical map changed")
            comparisons += 1
    if f.extract_table(changed) != phase_table(phase):
        raise AssertionError("declared phase table differs from source instrument")
    return comparisons


def network_certificate(phase):
    """Two overlapping native seams; a Hamiltonian interaction if Im phase != 0."""
    phase = unit_phase(phase)
    # Native carrier ports 0,1,5 contain this path; no screen or federation.
    carrier_ports = (0, 1, 5)
    seams = ((0, 1), (0, 2))
    source_edges = tuple(tuple(carrier_ports[p] for p in seam) for seam in seams)
    if source_edges != (native.E, native.F):
        raise AssertionError("chosen path is absent from the native carrier")
    channel_sum = sum((group_channel(native_binary_groups(phase, seam, 3), 8)
                       for seam in seams), np.full((64, 64), S(), dtype=object))
    skew = (channel_sum-adjoint(channel_sum))*F(1, 2)
    coefficient = -imaginary(phase)*F(1, 2)
    energies = tuple(coefficient*sum(((mask >> a) & 1) != ((mask >> b) & 1)
                                    for a, b in seams) for mask in range(8))
    if not np.array_equal(skew, commutator_generator(energies)):
        raise AssertionError("overlapping-seam Hamiltonian decomposition failed")
    # A nonzero mixed difference rules out a sum of independent site energies.
    interaction = energies[3]-energies[1]-energies[2]+energies[0]
    if interaction != imaginary(phase):
        raise AssertionError("wrong local interaction coefficient")
    return {"native_carrier_ports": list(carrier_ports),
            "native_carrier_seams": [list(e) for e in source_edges],
            "matrix_units_checked": 64,
            "exact_commutator_equality": True,
            "Hamiltonian_diagonal": list(map(str, energies)),
            "two_site_mixed_difference": str(interaction)}


def run():
    members = {}
    for label, phase in PHASES.items():
        table = phase_table(phase)
        certificate, _, _ = local_certificate(table)
        if not np.array_equal(binary_channel(table), group_channel(native_binary_groups(phase), 4)):
            raise AssertionError("independent native and table operators differ")
        certificate.update({
            "table": f.text(table),
            "unchanged_selective_classical_branches": source_selective_equality(phase),
            "same_C_fixed_operator_space": certificate["fixed_operator_dimension"] == 5,
            "wait_coherence_symmetric_action_coefficient": str(1-table["g"].real()),
            "coherent_probe_probabilities": {
                setting: str(f.plus_probability(f.apply_table(table, f.prepared(setting, f.oi.Z)), f.oi.Z))
                for setting in f.SETTINGS},
            "overlapping_native_seams": network_certificate(phase),
        })
        members[label] = certificate
    # A valid phase-recorded E law: rotate the full-output amplitude by -i.
    # This makes k=i, while g=h=m=n=0; it has a skew part which is not a derivation.
    obstruction_table = f.table(I, 0, 0, 0, 0)
    obstruction, _, skew = local_certificate(obstruction_table)
    # J(E_03 E_31) - J(E_03)E_31 - E_03 J(E_31) = -i E_01.
    defect = skew[1, 1]-skew[3, 3]-skew[13, 13]
    if defect != -I:
        raise AssertionError("skew-part Leibniz counterexample failed")
    return {
        "schema": "oph.repair_action_audit.v1",
        "domain": {"binary_local_Hilbert_dimension": 4,
                   "source_total_two_states_for_selective_checks": len(f.STATES),
                   "federation_required": False,
                   "generator_parameter": "optional unit-rate Poisson attempts per seam",
                   "discrete_transfer": "Phi, not exp(Phi-I)",
                   "compatibility": "fixed native endpoint law and actual coin/current selective diagonal maps"},
        "waiting_phase_family": members,
        "general_binary_action": {
            "self_adjoint_iff": "Im(k)=Im(g)=Im(h)=0",
            "symmetric_action": "I-(Phi+Phi*)/2 is positive semidefinite",
            "Hamiltonian_skew_iff": "Im(k)+Im(h)-Im(g)=0",
            "general_skew_counterexample": dict(obstruction, table=f.text(obstruction_table),
                                                Leibniz_defect_on_E03_E31=str(defect)),
        },
        "unavoidable_classical_mixing": {"input": "|10><10|", "input_purity": "1",
            "output": "(|10><10|+|01><01|)/2", "output_purity": "1/2",
            "applies_to_every_admissible_record_law": True},
        "scope": {"source_instrument_phase_is_extra_law_data": True,
            "same_classical_records_select_a_unique_action": False,
            "finite_native_reduced_process_is_unitary": False,
            "phase_family_introduces_no_new_classical_record_labels": True,
            "new_phase_law_selected_as_physics": False,
            "laboratory_clock_or_Lorentzian_action_derived": False,
            "binary_sector_result_not_claimed_on_collision_sector": True},
        "source_bindings": bound_pins({path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                            for path in SOURCE_PATHS}),
    }


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
