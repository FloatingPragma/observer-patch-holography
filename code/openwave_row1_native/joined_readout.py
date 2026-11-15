"""Complete actual-source -> native-record counter -> native-pointer readout.

The source is scheduled first and idles during terminal pointer service.
The detector reports the retained source prefix's relative charge. The
controller never receives the terminal loads or a precomputed degree.
"""
from fractions import Fraction as F
import hashlib
from itertools import cycle
import json
from pathlib import Path

import native_cycle_meter as meter

HERE = Path(__file__).resolve().parent


def source_prefix(initial_G, schedule, coins):
    """Advance the actual source and its finite record-only observer."""
    initial_G = tuple(initial_G)
    if len(initial_G) != 12 or any(type(x) is not int or x not in (0, 1) for x in initial_G) or sum(initial_G) != 1:
        return {'admitted': False, 'reason': 'source_outside_one_load_zero_tape_domain',
                'source_initial': list(initial_G), 'records': [], 'counter': None}
    schedule = tuple(tuple(edge) for edge in schedule)
    coins = tuple(coins)
    if len(coins) != len(schedule) or any(type(c) is not int or c not in (0, 1) for c in coins):
        raise ValueError('one literal coin for every scheduled source attempt')
    if any(edge not in meter.EDGES for edge in schedule):
        raise ValueError('actual oriented source seams required')
    loads = initial_G
    counter = [0]*12
    records = []
    for event, (edge, coin) in enumerate(zip(schedule, coins)):
        loads, (_, current) = meter.native(loads, edge, coin)
        # This is the entire observer update. It reads no occupation value.
        counter[edge[0]] = (counter[edge[0]]+current) % 6
        counter[edge[1]] = (counter[edge[1]]-current) % 6
        records.append({'stage': 'source', 'event': event, 'seam': list(edge),
                        'coin': coin, 'current': current, 'counter_after': list(counter)})
    return {'admitted': True, 'initial_G': list(initial_G), 'counter': counter,
            'records': records, 'source_final_loads_for_audit_only': list(loads),
            'source_history_probability': str(F(1, 2**len(records)))}


def execute_joined(initial_G, schedule, source_coins, region_mask,
                   pointer_initial_port, pointer_coins, retries=16):
    """Single joined policy, with source and pointer tapes in one output."""
    if type(region_mask) is not int or not 0 <= region_mask < 1 << 12:
        raise ValueError('declared twelve-port region mask required')
    source = source_prefix(initial_G, schedule, source_coins)
    if not source['admitted']:
        return {'success': False, 'reason': source['reason'], 'source': source,
                'pointer_executed': False, 'relative_charge': None, 'spatial_degree': None,
                'combined_native_tape': [], 'combined_history_probability': '1'}
    # The terminal source loads and stored expected degrees are unavailable
    # to this policy. Only its accumulated current-derived counter is read.
    regional_class = sum(source['counter'][v] for v in range(12) if region_mask >> v & 1) % 6
    pointer = meter.execute_class(pointer_initial_port, regional_class, pointer_coins, retries)
    combined = list(source['records'])
    for offset, row in enumerate(pointer['records'], start=len(combined)):
        combined.append({**row, 'stage': 'pointer_'+row['stage'], 'event': offset})
    probability = F(source['source_history_probability'])*F(pointer.get('fine_history_probability', '1'))
    return {'success': pointer['success'], 'reason': pointer['reason'],
            'region_mask': region_mask, 'source': source,
            'observer_regional_class': regional_class,
            'pointer_executed': True, 'pointer': pointer,
            'spatial_degree': pointer['winding'],
            'pointer_winding': pointer['winding'],
            'relative_charge': -pointer['winding'] if pointer['success'] else None,
            'combined_native_tape': combined,
            'combined_history_probability': str(probability),
            'clock': {'source_prefix_attempts': len(source['records']),
                      'pointer_attempts': len(pointer['records']),
                      'total_native_attempts': len(combined),
                      'source_during_pointer': 'idle under the declared addressed scheduler',
                      'reported_charge_time': 'end of the retained source prefix; unchanged through terminal meter'}}


def discovery_hold_coins(initial_port, retries):
    for a, b in meter.EDGES:
        for _ in range(retries):
            yield 1 if initial_port == a else 0

