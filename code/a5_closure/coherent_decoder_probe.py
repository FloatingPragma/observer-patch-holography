"""Finite fixed-template obstruction for a clean-register GV decoder.

No source execution, approximate decoding, or arbitrary residual-state
decoder is certified by the field obstruction.
"""
import argparse
import json
from pathlib import Path

import c12_current_attachment_probe as previous
from galois_port_frame_certificate import F, Z, encode, mm, transpose, require
from model.native_interface import tomography

ROOT = Path(__file__).resolve().parents[2]
BASELINE = 'eee3b6a2e0051e88d53053b5af1defbc0fb86988'
# Reviewed exact-value-preserving source update; retain the work baseline.
SOURCE_COMMIT = '83113a9cc33d82a56c84a8634445c8edd1e599a2'
REFERENCE = ROOT/'code/a5_closure/receipts/coherent_decoder_probe.json'
PATHS = sorted(set(previous.INPUTS + [
    'model/compiled_native_interface/construction.py',
    'model/native_interface/pointer.py',
    'model/native_quantum_history/construction.py',
    'model/observer/source_discrimination/circuit.py',
    'model/common_source/quantum_native.py',
    'model/common_source/native_routing.py',
]))


def build():
    pins = previous.custody(PATHS, SOURCE_COMMIT, 'FloatingPragma/observer-patch-holography')
    states = previous.native.states(binary=True)
    g = [[6*x[p]-1 for p in range(12)] for x in states]
    action_check = previous.slow_pair(g)
    b = mm([[F(x) for x in row] for row in g], previous.prior.marked.raw_frame())
    k = F(360)*F(10, 2)
    p = [[x/k for x in row] for row in mm(b, transpose(b))]
    require(sum((p[i][i] for i in range(66)), Z) == F(3), 'code rank')
    # A nonrational entry suffices: Q(sqrt(5)) intersects Q(i,sqrt(2)) in Q.
    i, j = next((i, j) for i in range(66) for j in range(66) if p[i][j].b)
    bound = max(2*len(tomography.route(a, b))+1
                for a in range(66) for b in range(a+1, 66))
    require(bound == 31, 'existing tomography resource envelope')
    return encode({
        'schema': 'oph.coherent-decoder-clean-register.v1',
        'baseline': BASELINE, 'source_custody': pins,
        'verdict': 'FINITE_GRAMMAR_GV_DECODER_OBSTRUCTED',
        'general_attachment_verdict': 'COHERENT_ACCESS_WITHOUT_ATTACHMENT',
        'grammar': {'qubits': 9, 'maximum_blocks': bound,
                    'alphabet_size': 9+9*8+2*9*8*7,
                    'templates': ['X', 'CNOT', 'mixer', 'y_analyzer'],
                    'placements': 'all ordered distinct wire placements',
                    'field': 'Q(i,sqrt(2))', 'free_angles': False},
        'target': {'receiver_bits': [0, 1], 'receiver_codes': [0, 1, 2],
                   'other_bits': 'fixed computational zero',
                   'equation': 'U J B/sqrt(k) = R',
                   'instrument': 'one unconditional unitary branch on all 512 addresses',
                   'postselection': False, 'arbitrary_residual_state': False},
        'projector': {'rank': 3, 'gram_scale': k, 'witness_index': [i, j],
                      'witness': p[i][j], 'sigma_witness': p[i][j].sigma(),
                      'sha256': previous.digest(p)},
        'marking': {'index': action_check['current_marking_index'],
                    'action': action_check['current_native_action'],
                    'actions_checked': action_check['actions_checked']},
        'matrix_units_required': 9, 'candidate_decoders': 0,
        'reference_extension_required': True,
        'approximation_lower_bound_claimed': False,
        'apparatus_service_constructed': False,
        'physical_golden_embedding_selected': False,
    })


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--write', type=Path)
    args = parser.parse_args()
    packet = build()
    if args.check:
        require(packet == json.loads(REFERENCE.read_text()), 'receipt replay')
    if args.write:
        args.write.write_text(json.dumps(packet, indent=2, sort_keys=True)+'\n')
    print(packet['verdict']+': fixed clean-register target only')
