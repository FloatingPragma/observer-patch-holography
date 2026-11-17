"""Independent exact comparisons, an actual native witness, and resource stops."""
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import random
import sys

from adaptive_phase import read_phase

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.dont_write_bytecode = True


def module(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def main():
    geometry_path = ROOT / 'openwave-tuning/row_completion/ROW_01/GEOMETRY.json'
    geom = json.loads(geometry_path.read_text())
    top = module('frozen_topology', 'openwave-tuning/row_completion/ROW_01/topology.py')
    native = module('frozen_native', 'reverse-engineering-reality/code/openwave_row1_integrity/native_kernel.py')
    witness_path = ROOT / 'openwave-push/b/G1_row1_strengthen/proposed_public/witness_four_turns.json'
    w = json.loads(witness_path.read_text())
    n = tuple(w['G'])
    current = [0] * 30
    schedule = w['events'] + list(zip(list(range(7, 30)) + list(range(7)), map(int, w['sweep_coins'])))
    for e, coin in schedule:
        n, (_, j) = native.native(n, tuple(geom['seams'][e]), coin)
        current[e] += j
    assert max(map(abs, current)) == 4
    exact = sum(w['G'][v] - n[v] for v in w['region'])
    assert top.discrete_curve_degree(geom, current, w['region']) == exact
    try:
        top.phase_degree(geom, current, w['region'], True)
    except ValueError:
        old_rejects = True
    else:
        old_rejects = False
    fixed = read_phase(geom, current, w['region'])
    assert old_rejects and fixed['status'] == 'read'
    assert abs(fixed['phase_degree'] - exact) < 1e-12
    rng = random.Random(711034)  # new diagnostic, never an acceptance seed
    max_error = 0.0
    counts = 0
    for amplitude in (0, 1, 3, 4, 8, 32, 128):
        for case in range(12):
            cochain = [rng.randint(-amplitude, amplitude) for _ in range(30)]
            region = tuple(v for v in range(12) if rng.getrandbits(1))
            q = sum(sign * cochain[e] for v in region for e, sign in geom['dual_port_boundaries'][v])
            for deform in (False, True):
                result = read_phase(geom, cochain, region, deformed=deform)
                assert result['status'] == 'read'
                err = abs(result['phase_degree'] - q)
                assert err < 2e-10
                assert result['max_principal_increment'] < math.pi/4 + 1e-12
                max_error = max(max_error, err)
                counts += 1
    # Rejecting controls: the historical bug and a low-resource readout
    # cannot silently become successful integer observations.
    stop = read_phase(geom, current, w['region'], sample_budget=0)
    assert stop['status'] == 'resource_stop' and stop['phase_degree'] is None
    huge = [10**100] + [0] * 29
    huge_stop = read_phase(geom, huge, tuple(range(12)))
    assert huge_stop['status'] == 'resource_stop'
    try:
        read_phase(geom, [.5] + [0] * 29, (0,))
    except ValueError:
        fraction_rejected = True
    else:
        fraction_rejected = False
    assert fraction_rejected
    pins = {}
    for p in (geometry_path, witness_path,
              ROOT/'openwave-tuning/row_completion/ROW_01/topology.py',
              ROOT/'reverse-engineering-reality/code/openwave_row1_integrity/native_kernel.py'):
        pins[str(p.relative_to(ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    out = {'schema': 'oph.obstruction_repair.adaptive_phase.v1',
           'problem_solved': 'fixed-resolution refusal on the honest four-turn witness',
           'source_law_changed': False, 'decoder_amendment': True,
           'physical_charge_identification_solved': False, 'official_rows_changed': [],
           'native_witness': {'current': current, 'final_loads': n, 'exact_degree': exact,
                              'old_rejects': old_rejects, 'amended_readout': fixed},
           'diagnostic_reads': counts, 'max_absolute_error': max_error,
           'controls': {'old_cap_rejects_witness': old_rejects,
                        'budget_stop_retained': stop['status'] == 'resource_stop',
                        'huge_integer_stops_before_float_conversion': huge_stop['status'] == 'resource_stop',
                        'fractional_input_rejected': fraction_rejected},
           'pins': pins}
    (HERE/'ADAPTIVE_RESULT.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({k: out[k] for k in ('problem_solved', 'diagnostic_reads', 'max_absolute_error', 'controls')}))


if __name__ == '__main__':
    main()
