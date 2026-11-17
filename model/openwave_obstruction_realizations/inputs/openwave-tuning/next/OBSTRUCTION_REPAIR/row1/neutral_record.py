"""Neutral pair preparation and a record-only native charge instrument.

The new preparation is a uniform integer filling and one recorded external
pair transport. Subsequent source and pointer events use the old native law.
Only a scalar filling, tapes and a region enter observer admission. No
terminal source loads or per-port preparation array are observer inputs.
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


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pair_sector(loads, filling):
    deviations = [v - filling for v in loads]
    return (sum(deviations) == 0 and all(v in (-1, 0, 1) for v in deviations)
            and deviations.count(1) <= 1 and deviations.count(-1) <= 1)


def execute_source(filling, pair_edge, pair_sign, prefix, sweep_coins):
    require(type(filling) is int and filling >= 1, 'positive integer neutral filling')
    require(type(pair_edge) is int and 0 <= pair_edge < 30 and type(pair_sign) is int and pair_sign in (-1, 1),
            'one recorded unit pair preparation')
    loads = [filling] * 12
    a, b = EDGES[pair_edge]
    loads[a] -= pair_sign
    loads[b] += pair_sign
    records = [{'kind': 'prepare_pair', 'seam': pair_edge, 'current': pair_sign}]
    order = tuple(range(7, 30)) + tuple(range(7))
    require(len(sweep_coins) == 30, 'thirty closing coins required')
    for kind, schedule in (('C', prefix), ('E', tuple(zip(order, sweep_coins)))):
        for e, coin in schedule:
            require(type(e) is int and 0 <= e < 30 and type(coin) is int and coin in (0, 1), 'native event')
            a, b = EDGES[e]
            d = loads[a] - loads[b]
            loads, (_, j) = original.native(tuple(loads), (a, b), coin)
            rec = {'kind': kind, 'seam': e, 'coin': coin, 'current': j}
            if kind == 'E':
                rec['contrast'] = d
            records.append(rec)
    return {'records': records, 'sweep_order': list(order)}, {'final_loads': list(loads)}


def observe_source(filling, tape, region):
    require(type(filling) is int and filling >= 1, 'positive integer neutral filling')
    region = original.region_set(region)
    require(isinstance(tape, dict) and isinstance(tape.get('records'), list), 'source records required')
    records = tape['records']
    require(len(records) >= 31, 'preparation plus closing sweep required')
    order = tape.get('sweep_order')
    require(isinstance(order, list) and all(type(e) is int for e in order)
            and sorted(order) == list(range(30)), 'complete closing sweep required')
    loads = [filling] * 12
    current = [0] * 30
    counters = [0] * 12
    for k, rec in enumerate(records):
        require(isinstance(rec, dict), 'event dictionary required')
        e, j = rec.get('seam'), rec.get('current')
        require(type(e) is int and 0 <= e < 30 and type(j) is int, 'integer current and seam required')
        a, b = EDGES[e]
        if k == 0:
            require(rec.get('kind') == 'prepare_pair' and j in (-1, 1), 'recorded unit pair intervention required')
        else:
            coin = rec.get('coin')
            require(type(coin) is int and coin in (0, 1), 'native coin required')
            d = loads[a] - loads[b]
            expected, (_, predicted) = original.native(tuple(loads), (a, b), coin)
            require(j == predicted, 'source current differs from native law')
            if k >= len(records)-30:
                require(rec.get('kind') == 'E' and e == order[k-(len(records)-30)], 'closing order mismatch')
                require(type(rec.get('contrast')) is int and rec['contrast'] == d, 'closing contrast mismatch')
            else:
                require(rec.get('kind') == 'C', 'unexpected source operation')
        loads = list(loads)
        loads[a] -= j
        loads[b] += j
        current[e] += j
        counters[a] = (counters[a] + j) % 6
        counters[b] = (counters[b] - j) % 6
        require(pair_sector(loads, filling), 'outside the single-pair invariant sector')
    residue = sum(counters[v] for v in region) % 6
    require(residue in (0, 1, 5), 'regional counter outside single-pair sector')
    degree = {0: 0, 1: 1, 5: -1}[residue]
    require(degree == sum(filling-loads[v] for v in region), 'counter/reconstruction inconsistency')
    return {'degree': degree, 'charge': -degree, 'regional_counter': residue,
            'current': current, 'native_attempts': len(records)-1,
            'preparation_interventions': 1}


def admit(filling, tape, region, pointer, retries=2):
    """Full fine-record effect, including every unsuccessful pointer outcome."""
    try:
        source = observe_source(filling, tape, region)
        require(type(retries) is int and 1 <= retries <= 128, 'bounded retry count required')
        require(isinstance(pointer, dict) and isinstance(pointer.get('records'), list), 'pointer tape required')
        clean = deepcopy(pointer)
        clean.pop('final_loads', None)
        columns = []
        for x in range(12):
            loads = tuple(int(v == x) for v in range(12))
            for rec in clean['records']:
                require(isinstance(rec, dict) and isinstance(rec.get('seam'), (list, tuple))
                        and tuple(rec['seam']) in EDGES and type(rec.get('coin')) is int
                        and rec['coin'] in (0, 1), 'malformed pointer record')
                loads, _ = original.native(loads, tuple(rec['seam']), rec['coin'])
            candidate = deepcopy(clean)
            candidate['final_loads'] = list(loads)
            errors, _ = original.strict.validate_pointer_tape(candidate, x, source['degree'], retries)
            if not errors:
                columns.append(x)
        require(bool(columns), 'no compatible unknown pointer initial state')
        weight = Fraction(1, 2 ** (source['native_attempts'] + len(clean['records'])))
        readout = source['charge'] if clean['success'] else None
        return {'admitted': True, 'readout': readout, 'success': clean['success'],
                'reason': clean['reason'], 'compatible_pointer_initials': columns,
                'fine_record_effect': [str(weight if x in columns else 0) for x in range(12)],
                'source_degree': source['degree'],
                'preparation_interventions': 1,
                'native_attempts': source['native_attempts'] + len(clean['records'])}
    except (ValueError, TypeError, KeyError, IndexError, OverflowError) as error:
        return {'admitted': False, 'reason': str(error)}
