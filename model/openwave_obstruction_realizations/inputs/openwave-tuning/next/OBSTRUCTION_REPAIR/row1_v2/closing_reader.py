"""Charge and a complete pointer effect from the thirty closing records alone.

No filling, initial pair, prefix, cumulative source tape or endpoint snapshot
is an observer input. The initial relative state is measured by the contrasts.
The instrument is conditional on a single-pair state at the start of this sweep.
"""
from copy import deepcopy
from fractions import Fraction
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[4]
PROVIDER = ROOT/'reverse-engineering-reality/code/openwave_row1_integrity'
sys.path.insert(0, str(PROVIDER))
import integrity as original

EDGES = original.EDGES


def need(condition, reason):
    if not condition:
        raise ValueError(reason)


def in_pair_sector(n):
    return (len(n) == 12 and sum(n) == 0 and all(type(x) is int and x in (-1, 0, 1) for x in n)
            and n.count(1) <= 1 and n.count(-1) <= 1)


def observe_closing(records, region):
    """Solve the connected difference equations for the sweep-start state."""
    region = original.region_set(region)
    need(isinstance(records, list) and len(records) == 30, 'exactly thirty closing records required')
    delta = [0]*12
    equations = []
    seen = set()
    for rec in records:
        need(isinstance(rec, dict) and rec.get('kind') == 'E', 'closing E record required')
        e, c, j, d = (rec.get(key) for key in ('seam', 'coin', 'current', 'contrast'))
        need(type(e) is int and 0 <= e < 30 and e not in seen, 'each seam must occur once')
        need(type(c) is int and c in (0, 1) and type(j) is int and type(d) is int,
             'integer native coin, current and contrast required')
        need(abs(j) <= 1 and abs(d) <= 2, 'record outside the finite single-pair alphabet')
        need(j == (d-(2*c-1)*(abs(d) % 2))//2, 'current contradicts native contrast law')
        seen.add(e)
        a, b = EDGES[e]
        equations.append((a, b, d-delta[a]+delta[b]))
        delta[a] -= j
        delta[b] += j
    values = {0: 0}
    pending = [0]
    while pending:
        v = pending.pop()
        for a, b, difference in equations:
            if a == v and b not in values:
                values[b] = values[a]-difference
                pending.append(b)
            elif b == v and a not in values:
                values[a] = values[b]+difference
                pending.append(a)
    need(len(values) == 12, 'disconnected observation graph')
    need(all(values[a]-values[b] == difference for a, b, difference in equations),
         'inconsistent closing contrasts')
    mean = Fraction(sum(values.values()), 12)
    initial = [Fraction(values[v])-mean for v in range(12)]
    need(all(x.denominator == 1 for x in initial), 'noninteger deviation from the mean')
    initial = [int(x) for x in initial]
    need(in_pair_sector(initial), 'measured initial state is outside the single-pair sector')
    state = initial.copy()
    for rec in records:
        a, b = EDGES[rec['seam']]
        need(state[a]-state[b] == rec['contrast'], 'reconstructed contrast mismatch')
        state[a] -= rec['current']
        state[b] += rec['current']
        need(in_pair_sector(state), 'native evolution left the declared sector')
    charge = sum(state[v] for v in region)
    need(charge in (-1, 0, 1), 'regional charge outside single-pair sector')
    return {'charge': charge, 'degree': -charge, 'regional_counter': (-charge) % 6,
            'initial_relative_state': initial, 'final_relative_state': state,
            'native_attempts': 30, 'absolute_filling_determined': False}


def admit_closing(records, region, pointer, retries=2):
    """Fine effect on the 133 relative source states x 12 pointer states.

    It is zero off the measured initial relative-state column; on that column
    its pointer diagonal is fine_pointer_effect. All failed outcomes remain.
    This is a sweep-start conditional experiment, not a prefix probability.
    """
    try:
        source = observe_closing(records, region)
        need(type(retries) is int and 1 <= retries <= 128, 'bounded retry count required')
        need(isinstance(pointer, dict) and isinstance(pointer.get('records'), list), 'pointer tape required')
        clean = deepcopy(pointer)
        clean.pop('final_loads', None)
        columns = []
        for x in range(12):
            loads = tuple(int(v == x) for v in range(12))
            for rec in clean['records']:
                need(isinstance(rec, dict) and isinstance(rec.get('seam'), (list, tuple))
                     and tuple(rec['seam']) in EDGES and type(rec.get('coin')) is int
                     and rec['coin'] in (0, 1), 'malformed pointer record')
                loads, _ = original.native(loads, tuple(rec['seam']), rec['coin'])
            candidate = deepcopy(clean)
            candidate['final_loads'] = list(loads)
            errors, _ = original.strict.validate_pointer_tape(candidate, x, source['degree'], retries)
            if not errors:
                columns.append(x)
        need(bool(columns), 'no compatible pointer initial state')
        attempts = 30+len(clean['records'])
        weight = Fraction(1, 2**attempts)
        return {'admitted': True, 'source_charge': source['charge'], 'source_degree': source['degree'],
                'source_initial_relative_column': source['initial_relative_state'],
                'compatible_pointer_initials': columns,
                'fine_pointer_effect': [str(weight if x in columns else 0) for x in range(12)],
                'effect_zero_on_other_source_columns': True,
                'success': clean['success'], 'reason': clean['reason'],
                'readout': source['charge'] if clean['success'] else None,
                'native_attempts': attempts, 'conditional_experiment': 'start of closing sweep'}
    except (ValueError, TypeError, KeyError, IndexError, OverflowError) as error:
        return {'admitted': False, 'reason': str(error)}
