"""Complete native repair on the full twelve-mode binary Fock carrier.

The endpoint operation and raw records are the existing C/E law.  Its
positive lift is the hard-core occupation instrument; its negative lift is
the local, even CAR instrument.  A sign is an actual coherent operation,
not an annotation on a classical record.  Neither lift is selected by the
unchanged occupation histories.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from math import comb

from ..defect_field_bridge import full_domain
from ..observer import native, record_recovery
from ..observer.source_discrimination import circuit


MODES = 12
DIMENSION = 1 << MODES
FULL = DIMENSION-1
SEAMS = tuple(record_recovery.seams())
_UNSET = object()


def _integer(value, low, high, name):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(f'exact integer {name} in {low}..{high} required')
    return value


def _character(character):
    if type(character) is not int or character not in (-1, 1):
        raise ValueError('an exact integer exchange character +1 or -1 required')
    return character


def _label(law, seam, coin, raw=_UNSET):
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('the native C or E law required')
    if (type(seam) is not tuple or len(seam) != 2
            or any(type(p) is not int for p in seam) or seam not in SEAMS):
        raise ValueError('a sorted actual native seam required')
    _integer(coin, 0, 1, 'coin')
    if raw is not _UNSET:
        _integer(raw, 0, (1 << circuit.RECORD_BITS[law])-1, 'raw record')
    return law, seam, coin


def states(number=None):
    """All occupation addresses; number restriction never changes the law."""
    if number is not None:
        _integer(number, 0, MODES, 'particle number')
    return tuple(mask for mask in range(DIMENSION)
                 if number is None or mask.bit_count() == number)


def column(mask, law, seam, coin, character=-1):
    """(raw, output, phase, current), before the common Kraus factor1/sqrt2.

    Routing is delegated to the existing full-CAR native extension at one
    flavor per port.  In this binary domain no collision-policy choice is
    involved.  The positive control removes only its moving CAR phase.
    """
    _integer(mask, 0, FULL, 'occupation mask')
    _label(law, seam, coin)
    _character(character)
    raw, output, phase, current = full_domain.branch_action(
        mask, law, coin, sites=MODES, flavors=1, seam=seam)
    if (any(type(x) is not int for x in (raw, output, phase, current))
            or not 0 <= raw < (1 << circuit.RECORD_BITS[law])
            or not 0 <= output < DIMENSION or phase not in (-1, 1)
            or current not in (-1, 0, 1)):
        raise ArithmeticError('binary repair left the exact native output domain')
    return raw, output, phase if character == -1 else 1, current


@dataclass(frozen=True, slots=True)
class Branch:
    law: str
    seam: tuple[int, int]
    coin: int
    raw: int
    inputs: tuple[int, ...]
    outputs: tuple[int, ...]
    signs: tuple[int, ...]
    abort: bool


def _coin_branches(law, seam, coin, character, number):
    rows = [[] for _ in range(1 << circuit.RECORD_BITS[law])]
    for mask in states(number):
        raw, output, phase, _ = column(mask, law, seam, coin, character)
        rows[raw].append((mask, output, phase))
    return tuple(Branch(law, seam, coin, raw,
                        *(tuple(x[i] for x in entries) for i in range(3)),
                        raw not in circuit.record_values(law, coin))
                 for raw, entries in enumerate(rows))


def branch(law, seam, coin, raw, character=-1, *, number=None):
    """One complete raw slot, including valid zero and abort outputs."""
    _label(law, seam, coin, raw)
    _character(character)
    states(number)
    return _coin_branches(law, seam, coin, character, number)[raw]


def branches(law, character=-1, *, number=None):
    _label(law, SEAMS[0], 0)
    _character(character)
    states(number)
    return tuple(row for seam in SEAMS for coin in (0, 1)
                 for row in _coin_branches(law, seam, coin, character, number))


@dataclass(frozen=True, slots=True)
class History:
    inputs: tuple[int, ...]
    outputs: tuple[int, ...]
    signs: tuple[int, ...]
    squared_amplitude: Fraction


def history_map(history, character=-1, *, number=None):
    """Compose committed records without normalizing or merging histories."""
    if type(history) is not tuple:
        raise ValueError('a tuple of complete native record labels required')
    _character(character)
    image = {mask: (mask, 1) for mask in states(number)}
    for label in history:
        if type(label) is not tuple or len(label) != 4:
            raise ValueError('each label needs law, seam, coin and raw record')
        row = branch(*label, character, number=number)
        lookup = {i: (j, s) for i, j, s in zip(row.inputs, row.outputs, row.signs, strict=True)}
        image = {i: (lookup[j][0], s*lookup[j][1])
                 for i, (j, s) in image.items() if j in lookup}
    inputs = tuple(image)
    return History(inputs, tuple(image[i][0] for i in inputs),
                   tuple(image[i][1] for i in inputs), Fraction(1, 2**len(history)))


def _ladder(mask, port, creation):
    if bool(mask & (1 << port)) == creation:
        return None
    return mask ^ (1 << port), (-1)**((mask & ((1 << port)-1)).bit_count())


def ladder(mask, port, *, creation):
    """An actual CAR ladder column on all number sectors; None means zero.

    This is an operator, not a normalized preparation or a free odd physical
    service.  Operational addition/removal needs its declared reservoir.
    """
    _integer(mask, 0, FULL, 'occupation mask')
    _integer(port, 0, MODES-1, 'port')
    if type(creation) is not bool:
        raise ValueError('an explicit creation Boolean required')
    return _ladder(mask, port, creation)


def particle_hole(mask):
    """C|n> = (-1)^sum(p*n_p)|complement(n)>; C c_p = c_p^dagger C."""
    _integer(mask, 0, FULL, 'occupation mask')
    return FULL ^ mask, (-1)**sum(p for p in range(MODES) if mask >> p & 1)


def particle_hole_label(law, seam, coin, raw):
    _label(law, seam, coin, raw)
    code = raw & 3
    return law, seam, 1-coin, (raw & ~3)+({1: 2, 2: 1}.get(code, code))


def _product(mask, first, second):
    value = _ladder(mask, *first)
    if value is None:
        return None
    target = _ladder(value[0], *second)
    return None if target is None else (target[0], value[1]*target[1])


def car_certificate():
    """Full-carrier anticommutators, with no number-sector truncation."""
    checks = 0
    for a in range(MODES):
        for b in range(MODES):
            for create_b in (False, True):
                for mask in range(DIMENSION):
                    entries = Counter()
                    for value in (_product(mask, (b, create_b), (a, False)),
                                  _product(mask, (a, False), (b, create_b))):
                        if value is not None:
                            entries[value[0]] += value[1]
                    expected = {mask: 1} if create_b and a == b else {}
                    if {k: v for k, v in entries.items() if v} != expected:
                        raise ArithmeticError('full-carrier CAR anticommutator failed')
                    checks += 1
    return {'modes': MODES, 'complete_occupation_states': DIMENSION,
            'all_number_sectors': list(range(MODES+1)),
            'anticommutator_columns': checks,
            'identities': ['{c_a,c_b}=0', '{c_a,c_b^dagger}=delta_ab I'],
            'odd_ladders_are_not_free_physical_preparations': True}


def complete_certificate():
    """Bind both lifts to all native endpoints, records, signs and currents."""
    report = {}
    for law in ('C', 'E'):
        sector_counts = [0]*(MODES+1)
        raw_slots = zero_slots = abort_slots = moving = positive_checked = 0
        for seam in SEAMS:
            incoming, outgoing = [0]*DIMENSION, [0]*DIMENSION
            for coin in (0, 1):
                rows = _coin_branches(law, seam, coin, -1, None)
                for row in rows:
                    raw_slots += 1
                    zero_slots += not row.inputs
                    abort_slots += row.abort
                    if len(set(row.outputs)) != len(row.outputs):
                        raise ArithmeticError('native branch is not a partial injection')
                    if row.abort and row.inputs:
                        raise ArithmeticError('native abort acquired binary support')
                    for mask, output, phase in zip(row.inputs, row.outputs, row.signs, strict=True):
                        loads = tuple((mask >> p) & 1 for p in range(MODES))
                        target, (_, current) = native.native(loads, seam, coin)
                        expected = sum(n << p for p, n in enumerate(target))
                        raw = {0: 0, 1: 1, -1: 2}[current]
                        if law == 'E':
                            raw += 4*((loads[seam[0]]+loads[seam[1]]) % 2)
                        if (output, row.raw) != (expected, raw):
                            raise ArithmeticError('lift disagrees with the actual native endpoint/record law')
                        positive = column(mask, law, seam, coin, 1)
                        if (type(positive) is not tuple or len(positive) != 4
                                or any(type(x) is not int for x in positive)
                                or positive != (raw, expected, 1, current)):
                            raise ArithmeticError('positive lift disagrees with the native endpoint/record/phase law')
                        positive_checked += 1
                        if output.bit_count() != mask.bit_count():
                            raise ArithmeticError('native lift changed total particle number')
                        if current:
                            a, b = seam if current == 1 else tuple(reversed(seam))
                            if _product(mask, (a, False), (b, True)) != (output, phase):
                                raise ArithmeticError('moving column is not the local CAR bilinear')
                            moving += 1
                        elif phase != 1:
                            raise ArithmeticError('native hold acquired a phase')
                        if any(((output >> p) & 1)-loads[p] !=
                               current*(int(p == seam[1])-int(p == seam[0]))
                               for p in range(MODES)):
                            raise ArithmeticError('local number continuity failed')
                        incoming[mask] += 1
                        outgoing[output] += 1
                        sector_counts[mask.bit_count()] += 1
            if set(incoming) != {2} or set(outgoing) != {2}:
                raise ArithmeticError('complete fair-coin instrument lost TP or unitality')
        if sector_counts != [60*comb(MODES, n) for n in range(MODES+1)]:
            raise ArithmeticError('native sector decomposition changed')
        if positive_checked != sum(sector_counts):
            raise ArithmeticError('positive native column inventory changed')
        report[law] = {'raw_slots': raw_slots, 'zero_slots': zero_slots,
                       'abort_slots': abort_slots, 'columns': sum(sector_counts),
                       'moving_CAR_columns': moving,
                       'positive_columns_checked': positive_checked,
                       'columns_by_number': sector_counts,
                       'complete_effect_and_output_sums': '2 before squared Kraus factor1/2',
                       'positive_control_changes_only_moving_phases': True}
    return report


def spectator_certificate():
    """Every complete fermionic branch commutes with outside CAR ladders.

    The C raw labels and endpoint parity together determine E.  A spectator
    ladder changes neither, so this same exact column check covers both laws.
    It exposes why the independently positive all-number instrument is not
    the same CAR-local continuation.
    """
    checks = 0
    for seam in SEAMS:
        for coin in (0, 1):
            columns = [column(mask, 'C', seam, coin) for mask in range(DIMENSION)]
            for mask, (raw, output, phase, _) in enumerate(columns):
                for port in range(MODES):
                    if port in seam:
                        continue
                    for creation in (False, True):
                        left = _ladder(mask, port, creation)
                        right = _ladder(output, port, creation)
                        if left is None:
                            if right is not None:
                                raise ArithmeticError('spectator support changed')
                        else:
                            after = columns[left[0]]
                            if (right is None or after[0] != raw or after[1] != right[0]
                                    or left[1]*after[2] != phase*right[1]):
                                raise ArithmeticError('native continuation fails CAR spectator compatibility')
                        checks += 1
    # A positive-control failure is retained, not silently called CAR local.
    mask, seam, port, coin = 1, (0, 5), 1, 0
    added = _ladder(mask, port, True)
    before = column(mask, 'C', seam, coin, 1)
    after = column(added[0], 'C', seam, coin, 1)
    right = _ladder(before[1], port, True)
    if added[1]*after[2] == before[2]*right[1]:
        raise ArithmeticError('positive-lift CAR incompatibility control disappeared')
    return {'spectator_ladder_columns': checks,
            'complete_C_and_E_raw_branch_intertwining': True,
            'identity': 'K_h c_p^(dagger) = c_p^(dagger) K_h for p outside the addressed seam',
            'positive_control_failure': {'mask': mask, 'seam': list(seam), 'coin': coin,
                                         'created_spectator': port, 'raw': before[0],
                                         'left_sign': added[1]*after[2],
                                         'right_sign': before[2]*right[1]},
            'no_conditioned_branch_normalization': True}


def particle_hole_certificate():
    phases = [particle_hole(mask) for mask in range(DIMENSION)]
    for mask, (complement, phase) in enumerate(phases):
        if phases[complement] != (mask, phase):
            raise ArithmeticError('particle-hole map lost its involution')
        for port in range(MODES):
            left, right = _ladder(mask, port, False), _ladder(complement, port, True)
            if (left is None) != (right is None):
                raise ArithmeticError('particle-hole ladder support changed')
            if left is not None and (phases[left[0]][0], left[1]*phases[left[0]][1]) != (
                    right[0], phase*right[1]):
                raise ArithmeticError('particle-hole CAR conjugation failed')
    checks = 0
    for law in ('C', 'E'):
        for seam in SEAMS:
            for coin in (0, 1):
                for raw in range(1 << circuit.RECORD_BITS[law]):
                    other_raw = particle_hole_label(law, seam, coin, raw)[3]
                    if ((raw in circuit.record_values(law, coin)) !=
                            (other_raw in circuit.record_values(law, 1-coin))):
                        raise ArithmeticError('particle-hole conjugation changed the complete abort menu')
                for mask, (complement, phase) in enumerate(phases):
                    raw, output, sign, current = column(mask, law, seam, coin)
                    other = column(complement, law, seam, 1-coin)
                    expected_raw = particle_hole_label(law, seam, coin, raw)[3]
                    relative = -1 if current else 1
                    if (other[0] != expected_raw or other[1] != (FULL ^ output)
                            or other[3] != -current
                            or sign*phases[output][1] != relative*phase*other[2]):
                        raise ArithmeticError('fermionic particle-hole repair intertwiner failed')
                    checks += 1
    # Native triangle 0-5-11. One particle and one hole have opposite loop
    # phases, although their complemented occupation histories agree.
    triangle = (0, 5, 11, 0)
    loops = []
    for hole in (False, True):
        mask = FULL ^ 1 if hole else 1
        initial, total, labels = mask, 1, []
        for a, b in zip(triangle[:-1], triangle[1:], strict=True):
            seam = tuple(sorted((a, b)))
            wanted = mask ^ (1 << a) ^ (1 << b)
            matches = [(coin, column(mask, 'C', seam, coin)) for coin in (0, 1)]
            matches = [(coin, item) for coin, item in matches if item[1] == wanted]
            if len(matches) != 1:
                raise ArithmeticError('native triangle lost its unique moving record')
            coin, (raw, mask, sign, _) = matches[0]
            total *= sign
            labels.append(['C', list(seam), coin, raw])
        if mask != initial or total != (-1 if hole else 1):
            raise ArithmeticError('particle-hole odd-cycle phase control failed')
        loops.append({'number': MODES-1 if hole else 1, 'initial_mask': initial,
                      'history': labels, 'phase': total, 'history_probability': '1/8'})
    return {'all_native_column_checks': checks,
            'C_squared': 'I', 'ladder_conjugation': 'C c_p C = c_p^dagger',
            'native_intertwiner': 'C K_h = (-1)^abs(j_h) K_complement(h) C',
            'complemented_coin': '1-coin', 'complemented_current': '-current',
            'complete_zero_and_abort_menus_retained': True,
            'odd_triangle_loop_controls': loops,
            'uniform_moving_minus_removable_by_site_phases': False,
            'reason': 'the product of a site-phase coboundary around any triangle is +1, whereas three moving minus signs multiply to -1',
            'electromagnetic_charge_conjugation_identified': False}


def certificate():
    return {'schema': 'oph.native_exchange.many_body.v1',
            'carrier': {'modes': MODES, 'dimension': DIMENSION,
                        'sector_dimensions': [comb(MODES, n) for n in range(MODES+1)]},
            'native': complete_certificate(), 'CAR': car_certificate(),
            'creation_continuation': spectator_certificate(),
            'particle_hole': particle_hole_certificate(),
            'classical_record_and_terminal_population_laws_agree_for_both_lifts': True,
            'physical_statistics_or_particle_species_selected': False}
