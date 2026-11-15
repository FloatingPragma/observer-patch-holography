"""Native occupancy service with equal conclusive efficiency for every source output.

The source is the frozen one-load, known-reference native process. Its
terminal apparatus uses six actual hops after reset for every class.
"""
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
FROZEN = ROOT / 'openwave-tuning/next/ROW1_NATIVE/routes/auditor_bridge'
PINS = {'native_cycle_meter.py': '4db8d707879c416ca5002b782a70e07aa71d4d9a7fe7111b76a64828c1a05329',
        'joined_readout.py': '1bf9bda370c20a6a23c253654f5fe60160bca20fd5d24b0548bac5e98ba6b467'}


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


for name, pin in PINS.items():
    if hashlib.sha256((FROZEN / name).read_bytes()).hexdigest() != pin:
        raise ValueError('frozen native readout dependency changed: ' + name)
meter = load('balanced_occupancy_meter', FROZEN / 'native_cycle_meter.py')
_prior_meter = sys.modules.get('native_cycle_meter')
sys.modules['native_cycle_meter'] = meter
joined = load('balanced_occupancy_source', FROZEN / 'joined_readout.py')
if _prior_meter is None:
    del sys.modules['native_cycle_meter']
else:
    sys.modules['native_cycle_meter'] = _prior_meter


def native_meter(initial_port, regional_class, coins, retries=16):
    """Complete equal-efficiency meter; zero follows three out-and-back pairs."""
    if regional_class not in (0, 1, 5):
        return {'success': False, 'reason': 'class_outside_one_load_domain',
                'records': [], 'winding': None, 'fine_history_probability': '1'}
    iterator = iter(coins)
    reset = meter.execute(initial_port, (), iterator, retries)
    if not reset['success']:
        return reset
    records = list(reset['records'])
    # State is simulator storage; the policy's current location after reset
    # is the certified fixed port zero. No branch reads hidden occupation.
    loads = tuple(reset['final_loads'])
    known_port = 0
    signed_steps = 0
    direction = {0: 0, 1: 1, 5: -1}[regional_class]
    path = ((1, 0, 1, 0, 1, 0) if direction == 0 else
            tuple(meter.CYCLE[(direction*k) % 6] for k in range(1, 7)))
    for target in path:
        edge = tuple(sorted((known_port, target)))
        if edge not in meter.EDGES:
            raise ValueError('pointer path is outside native seams')
        requested = 1 if meter.CYCLE[(meter.CYCLE.index(known_port)+1) % 6] == target else -1
        moved = False
        for _ in range(retries):
            coin = next(iterator)
            if type(coin) is not int or coin not in (0, 1):
                raise ValueError('literal native coin required')
            loads, (_, current) = meter.native(loads, edge, coin)
            records.append({'stage': 'balanced_measure', 'seam': list(edge),
                            'coin': coin, 'current': current,
                            'requested_cycle_direction': requested})
            if current:
                arrival = edge[1] if current > 0 else edge[0]
                if arrival != target:
                    raise ValueError('native pointer did not follow its requested seam')
                known_port = arrival
                signed_steps += requested
                moved = True
                break
        if not moved:
            return {'success': False, 'reason': 'measurement_timeout',
                    'records': records, 'winding': None, 'final_loads': list(loads),
                    'retained_signed_cycle_steps': signed_steps,
                    'fine_history_probability': str(F(1, 2**len(records)))}
    assert known_port == 0 and signed_steps % 6 == 0
    return {'success': True, 'reason': 'complete_balanced_native_winding',
            'records': records, 'winding': signed_steps // 6,
            'final_loads': list(loads), 'retained_signed_cycle_steps': signed_steps,
            'fine_history_probability': str(F(1, 2**len(records)))}


def execute_occupancy(initial_G, schedule, source_coins, region_mask,
                      pointer_initial_port, pointer_coins, retries=16):
    """Execute source, observed-current counter and balanced native apparatus."""
    if type(region_mask) is not int or not 0 <= region_mask < 1 << 12:
        raise ValueError('declared twelve-port region required')
    source = joined.source_prefix(initial_G, schedule, source_coins)
    if not source['admitted']:
        return {'success': False, 'reason': source['reason'], 'source': source,
                'pointer_executed': False, 'occupancy': None,
                'relative_charge': None, 'spatial_degree': None,
                'combined_native_tape': [], 'combined_history_probability': '1'}
    cls = sum(source['counter'][v] for v in range(12) if region_mask >> v & 1) % 6
    pointer = native_meter(pointer_initial_port, cls, pointer_coins, retries)
    combined = list(source['records'])
    for i, event in enumerate(pointer['records'], len(combined)):
        combined.append({**event, 'stage': 'pointer_' + event['stage'], 'event': i})
    degree = pointer['winding']
    initial_occupancy = sum(initial_G[v] for v in range(12) if region_mask >> v & 1)
    occupancy = initial_occupancy - degree if pointer['success'] else None
    assert occupancy in (None, 0, 1)
    return {'success': pointer['success'], 'reason': pointer['reason'],
            'source': source, 'pointer_executed': True, 'pointer': pointer,
            'region_mask': region_mask, 'observer_regional_class': cls,
            'spatial_degree': degree, 'relative_charge': -degree if pointer['success'] else None,
            'occupancy': occupancy, 'combined_native_tape': combined,
            'combined_history_probability': str(F(source['source_history_probability']) * F(pointer['fine_history_probability'])),
            'clock': {'source_prefix_attempts': len(source['records']),
                      'pointer_attempts': len(pointer['records']),
                      'total_native_attempts': len(combined),
                      'source_during_pointer': 'idle under the declared addressed scheduler',
                      'reported_occupancy_time': 'end of retained source prefix; unchanged during terminal meter'}}


def complete_law(pointer_initial_port, retries=16):
    """Same law for every admitted regional class, conditional on pointer input."""
    reset = meter.reset_success(pointer_initial_port, retries)
    success = reset * (1 - F(1, 2**retries))**6
    return {'success': success, 'timeout': 1 - success}
