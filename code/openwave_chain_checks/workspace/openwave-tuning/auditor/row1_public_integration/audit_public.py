"""Independent custody, execution and complete-law audit of the portable row-1 package."""
import ast
from collections import deque
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import importlib.util
from itertools import cycle
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PUBLIC = ROOT / 'reverse-engineering-reality/code/openwave_row1_native'
FROZEN = ROOT / 'openwave-tuning/next/ROW1_NATIVE/routes/auditor_bridge'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def functions(path):
    source = path.read_text()
    return {node.name: ast.get_source_segment(source, node)
            for node in ast.parse(source).body if isinstance(node, ast.FunctionDef)}


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@lru_cache(maxsize=None)
def independent_partition(native, edges, initial, retries, hops):
    """Direct native branches followed by an independently rebuilt route distance."""
    distances = {0: 0}
    queue = deque([0])
    while queue:
        x = queue.popleft()
        for a, b in edges:
            if x in (a, b):
                y = b if x == a else a
                if y not in distances:
                    distances[y] = distances[x] + 1
                    queue.append(y)
    state = tuple(int(v == initial) for v in range(12))
    undiscovered = F(1)
    discovered = {x: F() for x in range(12)}
    for edge in edges:
        for _ in range(retries):
            waiting = F()
            for coin in (0, 1):
                after, (_, current) = native(state, edge, coin)
                if current:
                    actual_arrival = after.index(1)
                    discovered[actual_arrival] += undiscovered / 2
                else:
                    assert after == state
                    waiting += undiscovered / 2
            undiscovered = waiting
    hop = 1 - F(1, 2**retries)
    reset = sum(p * hop**distances[x] for x, p in discovered.items())
    result = {'discovery_timeout': undiscovered,
              'routing_timeout': sum(discovered.values()) - reset,
              'measurement_timeout': reset * (1 - hop**hops),
              'success': reset * hop**hops}
    assert sum(result.values()) == 1 and all(p >= 0 for p in result.values())
    return result


def run():
    manifest = json.loads((PUBLIC / 'manifest.json').read_text())
    monitored = [PUBLIC / name for name in manifest['files_sha256']]
    monitored += [PUBLIC / 'manifest.json', FROZEN / 'joined_readout.py',
                  FROZEN / 'native_cycle_meter.py', FROZEN / 'JOINED_READOUT.json']
    packet = json.loads((PUBLIC / 'source_packet.json').read_text())
    monitored += [ROOT / name for name in packet['source_files_sha256']]
    before = {str(p.relative_to(ROOT)): sha(p) for p in monitored}
    for name, pin in manifest['files_sha256'].items():
        assert sha(PUBLIC / name) == pin, name
    for name, pin in packet['source_files_sha256'].items():
        assert sha(ROOT / name) == pin, name
    original_spatial = json.loads((ROOT / 'openwave-tuning/next/ROW1_NATIVE/routes/spatial/NATIVE_DUAL_WINDING.json').read_text())
    assert packet['spatial_result'] == original_spatial
    origins = json.loads((PUBLIC / 'provenance/ORIGINS.json').read_text())
    retained = {}
    for module in ('native_cycle_meter', 'joined_readout'):
        pub_functions = functions(PUBLIC / (module + '.py'))
        old_functions = functions(FROZEN / (module + '.py'))
        names = origins[module]['retained_function_sha256']
        assert set(pub_functions) == set(names)
        for name, pin in names.items():
            assert pub_functions[name] == old_functions[name], (module, name)
            assert hashlib.sha256(pub_functions[name].encode()).hexdigest() == pin
        retained[module] = sorted(names)
        assert (PUBLIC / ('provenance/' + module + '.original.py')).read_bytes() == (FROZEN / (module + '.py')).read_bytes()
    canonical_native = functions(ROOT / 'reverse-engineering-reality/model/observer/native.py')['native']
    assert canonical_native == functions(PUBLIC / 'native_kernel.py')['native']
    assert hashlib.sha256(canonical_native.encode()).hexdigest() == packet['native_function_sha256']

    old_meter = load('audit_frozen_meter', FROZEN / 'native_cycle_meter.py')
    sys.modules['native_cycle_meter'] = old_meter
    old_joined = load('audit_frozen_joined', FROZEN / 'joined_readout.py')
    sys.path.insert(0, str(PUBLIC))
    new_meter = load('audit_public_meter', PUBLIC / 'native_cycle_meter.py')
    sys.modules['native_cycle_meter'] = new_meter
    new_joined = load('audit_public_joined', PUBLIC / 'joined_readout.py')
    assert old_meter.EDGES == new_meter.EDGES == tuple(map(tuple, packet['seams']))
    assert old_meter.CYCLE == new_meter.CYCLE

    replays = 0
    reasons = set()
    law_digest = hashlib.sha256()
    laws = 0
    for initial in range(12):
        for retries in (1, 2, 16):
            for degree in (-1, 0, 1):
                partition = independent_partition(new_meter.native, new_meter.EDGES, initial, retries, 6 * abs(degree))
                reset = new_meter.reset_success(initial, retries)
                assert partition['success'] == reset * (1 - F(1, 2**retries))**(6 * abs(degree))
                law_digest.update(json.dumps({k: str(v) for k, v in partition.items()}, sort_keys=True).encode())
                laws += 1
    source_histories = 0
    joined_budget_count = 0
    for index, seed in enumerate(packet['spatial_result']['seed_records']):
        G = seed['initial_G']
        schedule = [tuple(event['seam']) for event in seed['actual_history']]
        source_laws = {port: {degree: F() for degree in (-1, 0, 1)} for port in range(12)}
        for branch in seed['complete_prescribed_schedule_histories']:
            coins = [event['coin'] for event in branch['events']]
            source_histories += 1
            for port in range(12):
                initial = (index + port) % 12
                pointer_coins = (0,) if port % 3 == 0 else (1,) if port % 3 == 1 else (0, 1)
                left = old_joined.execute_joined(G, schedule, coins, 1 << port, initial, cycle(pointer_coins), 2)
                right = new_joined.execute_joined(G, schedule, coins, 1 << port, initial, cycle(pointer_coins), 2)
                assert left == right
                assert F(right['combined_history_probability']) == F(1, 2**len(right['combined_native_tape']))
                assert right['clock']['source_during_pointer'] == 'idle under the declared addressed scheduler'
                assert right['clock']['total_native_attempts'] == len(right['combined_native_tape'])
                if not right['success']:
                    assert right['relative_charge'] is None and right['pointer_winding'] is None
                reasons.add(right['reason'])
                degree = G[port] - branch['final_loads'][port]
                source_laws[port][degree] += F(branch['probability'])
                replays += 1
        for port in range(12):
            assert sum(source_laws[port].values()) == 1
            for initial in range(12):
                outcome = {name: F() for name in ('discovery_timeout', 'routing_timeout', 'measurement_timeout', 'success')}
                for degree, weight in source_laws[port].items():
                    law = independent_partition(new_meter.native, new_meter.EDGES, initial, 16, 6 * abs(degree))
                    for name, probability in law.items():
                        outcome[name] += weight * probability
                assert sum(outcome.values()) == 1
                law_digest.update(json.dumps({k: str(v) for k, v in outcome.items()}, sort_keys=True).encode())
                joined_budget_count += 1
    assert source_histories == 456 and replays == 5472 and joined_budget_count == 19008

    frozen_examples = json.loads((FROZEN / 'JOINED_READOUT.json').read_text())['joined_execution_examples']
    example_checks = 0
    for archived in frozen_examples:
        if archived['reason'] == 'source_outside_one_load_zero_tape_domain':
            replay = new_joined.execute_joined([7] + [0]*11, (), (), 1, 0, iter(()))
        else:
            source = archived['source']
            initial = archived['pointer']['initial_port'] if 'initial_port' in archived['pointer'] else None
            # The unknown apparatus input is recovered for auditing by replaying
            # all admitted input columns against the retained native tape.
            candidates = []
            for x in range(12):
                try:
                    candidate = new_joined.execute_joined(source['initial_G'],
                        [tuple(e['seam']) for e in source['records']],
                        [e['coin'] for e in source['records']], archived['region_mask'], x,
                        [e['coin'] for e in archived['pointer']['records']])
                except StopIteration:
                    continue
                if candidate == {k: v for k, v in archived.items() if k != 'case'}:
                    candidates.append(candidate)
            assert candidates, archived['case']
            replay = candidates[0]
        assert replay == {k: v for k, v in archived.items() if k != 'case'}, archived['case']
        example_checks += 1
    assert example_checks == 11

    with tempfile.TemporaryDirectory(prefix='oph-row1-public-audit-') as temporary:
        isolated = Path(temporary)
        for name in manifest['files_sha256']:
            target = isolated / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(PUBLIC / name, target)
        shutil.copyfile(PUBLIC / 'manifest.json', isolated / 'manifest.json')
        env = {k: v for k, v in os.environ.items() if k not in ('PYTHONPATH', 'PYTHONHOME')}
        env['PYTHONDONTWRITEBYTECODE'] = '1'
        completed = subprocess.run([sys.executable, '-S', str(isolated / 'run.py'), '--output', str(isolated / 'result.json')],
                                   cwd=isolated, env=env, capture_output=True, text=True, check=True)
        isolated_result = json.loads((isolated / 'result.json').read_text())
        assert all(isolated_result['checks'].values())
        for template in isolated_result['apparatus_templates']:
            law = independent_partition(new_meter.native, new_meter.EDGES,
                template['initial_pointer_port'], 16, 6*abs(template['decoded_spatial_winding']))
            assert F(template['unconditional_law']['success']) == law['success']
            assert F(template['unconditional_law']['timeout']) == 1 - law['success']
        isolated_summary = {k: isolated_result[k] for k in ('checks', 'source_histories', 'source_events', 'source_noops', 'joined_executor_comparisons', 'worst_primary_meter_timeout')}
    after = {str(p.relative_to(ROOT)): sha(p) for p in monitored}
    assert before == after, 'reviewed source drift: ' + ', '.join(name for name in before if before[name] != after[name])
    return {'schema': 'oph.row1.public_integration_audit.v1',
            'checks': {'all_manifest_and_input_ancestry_pins_match': True,
                       'source_packet_matches_frozen_spatial_receipt': True,
                       'retained_function_bodies_and_native_kernel_exact': True,
                       'all5472_public_vs_frozen_joined_executions_equal': True,
                       'all19008_joint_laws_complete_with_distinct_timeouts': True,
                       'all11_frozen_examples_reproduced': True,
                       'isolated_standard_library_primary_passes': True,
                       'reviewed_sources_unchanged': True},
            'retained_functions': retained, 'actual_joined_comparisons': replays,
            'independent_meter_stage_laws': laws, 'complete_joint_laws': joined_budget_count,
            'observed_execution_reasons': sorted(reasons), 'frozen_example_replays': example_checks,
            'law_digest': law_digest.hexdigest(), 'isolated_primary': isolated_summary,
            'source_sha256': after, 'reviewer_sha256': sha(Path(__file__)),
            'scope': 'Declared one-load native source with retained initial reference, separate auxiliary pointer and addressed source-idle scheduler. Complete laws conditional on configured source/region/pointer input, with no input prior invented.',
            'official_rows_changed': []}


if __name__ == '__main__':
    result = run()
    (HERE / 'PUBLIC_INTEGRATION_AUDIT.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('checks', 'actual_joined_comparisons', 'complete_joint_laws', 'observed_execution_reasons')}, indent=2))
