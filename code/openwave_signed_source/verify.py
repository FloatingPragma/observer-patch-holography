"""Verify the signed-source evidence and recompute its three results.

Runs the custody check of the published experiment package
(code/openwave_static_signed_source, every evidence file against its manifest),
checks the sealed receipts, and recomputes the charge, rest-energy and frame
results from the carrier graph and the native repair rule (recompute.py).

    python3 code/openwave_signed_source/verify.py
"""
from __future__ import annotations

import importlib.util
import json
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent / 'openwave_static_signed_source'
SHARED = PACKAGE / 'evidence' / 'openwave' / 'shared'


def read(path):
    return json.loads(path.read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def custody():
    out = subprocess.run([sys.executable, '-B', str(PACKAGE / 'replay.py'), '--verify-only'],
                         capture_output=True, text=True, cwd=HERE.parents[1])
    require(out.returncode == 0, 'published evidence custody: ' + out.stderr[-400:])
    result = json.loads(out.stdout.strip().splitlines()[-1])
    require(result['verified_files'] == 316, 'published evidence file count')
    return result['verified_files']


def receipts():
    static = read(SHARED / 'static_candidate_join' / 'RECEIPT.json')
    require(static['all_pass'] and static['native_calls_total'] == 435080, 'native calls')
    require(static['charge_trials'] == 264 and static['native_conclusive'] == 264, 'charge trials')
    require(static['control_trials'] == 48 and static['controls_never_synchronize'], 'controls')
    plant = static['fractional_plant']
    require(plant['regions'] == {'0': '1/2', '1': '1/2'} and plant['constant_integer_decoder_rejected'],
            'planted half-current')
    moment = static['moment']
    require(moment['branch_alternatives'] == 720 and moment['proper_half_turn_tensor_return']
            and moment['polar_frame_reversal'] and moment['squared_frame_return'], 'half-turn')
    require(moment['fixed_seam_failed_inputs'] == 4 and moment['off_axis_plant_detected'], 'frame controls')
    require(math.isclose(static['pair_energy_keV'], 1024.4508487023415, rel_tol=0, abs_tol=1e-9), 'pair')
    compact = read(SHARED / 'compact_static_adjudication' / 'RECEIPT.json')
    require(compact['pair_all_columns_checked'] == 82575360, 'compact columns')
    return {'native_calls': static['native_calls_total'], 'charge_trials': static['charge_trials'],
            'changed_law_controls': static['control_trials'],
            'frame_branches': moment['branch_alternatives'],
            'compact_columns_checked': compact['pair_all_columns_checked']}


def verify():
    spec = importlib.util.spec_from_file_location('signed_source_recompute', HERE / 'recompute.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    graph, seams = module.load_graph()
    cell = [e for e in seams if e[1] < 12]
    recomputed = {'charge': module.charge(cell), 'rest_energy': module.energy(graph, seams),
                  'frame': module.frame(cell)}
    files = custody()
    sealed = receipts()
    require(recomputed['charge']['native_trials_absorbed'] == sealed['charge_trials'], 'trial count')
    require(recomputed['frame']['branches'] == sealed['frame_branches'], 'branch count')
    return {'published_evidence_files_verified': files,
            'sealed_native_calls': sealed['native_calls'],
            'full_replay': 'python3 -B code/openwave_static_signed_source/replay.py',
            'recomputed': recomputed}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
