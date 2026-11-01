"""Independent spectral-projector obstruction; never imports its producer."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

import verify_c12_current_attachment_independent as old
import verify_coherent_native_access_independent as access
from verify_galois_port_frame_independent import ZERO, scalar, sub, mul, div, total, sigma

ROOT = Path(__file__).resolve().parents[2]
BASELINE = 'eee3b6a2e0051e88d53053b5af1defbc0fb86988'
# Reviewed exact-value-preserving source update; retain the work baseline.
SOURCE_COMMIT = '83113a9cc33d82a56c84a8634445c8edd1e599a2'
REFERENCE = ROOT/'code/a5_closure/receipts/coherent_decoder_probe.json'
m = old.m
PATHS = sorted(set(old.INPUTS + [
    'model/compiled_native_interface/construction.py',
    'model/native_interface/pointer.py',
    'model/native_quantum_history/construction.py',
    'model/observer/source_discrimination/circuit.py',
    'model/common_source/quantum_native.py',
    'model/common_source/native_routing.py',
]))


def verify(packet, root=ROOT):
    pins = old.custody(PATHS, SOURCE_COMMIT, 'FloatingPragma/observer-patch-holography', root)
    _, distances = m.incidence(root)
    lap = [[scalar(5*int(i == j)-int(distances[i][j] == 1))
            for j in range(12)] for i in range(12)]
    slow = m.interpolate(lap, sub(scalar(5), m.S))
    states = old.configurations()
    g = [[scalar(6*x[p]-1) for p in range(12)] for x in states]
    # Unlike the producer's B B*/k, reconstruct from the spectral polynomial.
    p = [[div(x, scalar(360)) for x in row]
         for row in m.mm(m.mm(g, slow), old.prior.trans(g))]
    m.need(total(p[i][i] for i in range(66)) == scalar(3), 'rank three')
    i, j = next((i, j) for i in range(66) for j in range(66) if p[i][j][1])
    m.need(p[i][j] != sigma(p[i][j]), 'nonrational projector witness')
    frame, vertices = m.source_frame(root), m.response_vertices()
    mark = m.source_marked_action(root, vertices, [vertices.index(v) for v in frame])
    frozen = m.load(root/'code/a5_closure/manifests/marked_triplet_attachment_reference.json')
    m.need(frozen['actions'].index(list(mark)) == 4, 'fixed marker')
    b = m.mm(g, frame)
    k = (Q(2880), Q(1440))
    m.need(m.mm(old.prior.trans(b), b) == [[k if a == c else ZERO for c in range(3)] for a in range(3)], 'all three columns')
    lookup = {x: n for n, x in enumerate(states)}
    source = m.load(root/'code/source_selection_model/response.json')['closed_paths']
    m.need(len(source) == len(frozen['actions']) == 60, 'sixty actions')
    for action, row in zip(frozen['actions'], source, strict=True):
        u = [[(Q(x[0])-Q(x[1]), 2*Q(x[1])) for x in line]
             for line in row['observed_response'][0]]
        inverse = [action.index(n) for n in range(12)]
        for n, x in enumerate(states):
            image = lookup[tuple(x[inverse[t]] for t in range(12))]
            m.need(b[image] == [total(mul(u[a][c], b[n][c]) for c in range(3)) for a in range(3)], 'marked column intertwining')
    codes = access.source_codes()
    bound = max(2*len(access.route(codes, a, b))+1 for a in range(66) for b in range(a+1, 66))
    m.need(bound == 31, 'pre-existing route envelope')
    # Rebuild small source words and interval-check against exact field targets.
    access.raw_source_gates()
    expected = old.prior.encoded({
        'schema': 'oph.coherent-decoder-clean-register.v1',
        'baseline': BASELINE, 'source_custody': pins,
        'verdict': 'FINITE_GRAMMAR_GV_DECODER_OBSTRUCTED',
        'general_attachment_verdict': 'COHERENT_ACCESS_WITHOUT_ATTACHMENT',
        'grammar': {'qubits': 9, 'maximum_blocks': 31,
                    'alphabet_size': 1089,
                    'templates': ['X', 'CNOT', 'mixer', 'y_analyzer'],
                    'placements': 'all ordered distinct wire placements',
                    'field': 'Q(i,sqrt(2))', 'free_angles': False},
        'target': {'receiver_bits': [0, 1], 'receiver_codes': [0, 1, 2],
                   'other_bits': 'fixed computational zero',
                   'equation': 'U J B/sqrt(k) = R',
                   'instrument': 'one unconditional unitary branch on all 512 addresses',
                   'postselection': False, 'arbitrary_residual_state': False},
        'projector': {'rank': 3, 'gram_scale': k, 'witness_index': [i, j],
                      'witness': p[i][j], 'sigma_witness': sigma(p[i][j]),
                      'sha256': old.digest(p)},
        'marking': {'index': 4, 'action': list(mark), 'actions_checked': 60},
        'matrix_units_required': 9, 'candidate_decoders': 0,
        'reference_extension_required': True,
        'approximation_lower_bound_claimed': False,
        'apparatus_service_constructed': False,
        'physical_golden_embedding_selected': False,
    })
    m.need(packet == expected, 'complete decoder obstruction receipt')
    return {'verified': True, 'verdict': expected['verdict'],
            'scope': 'exact ideal fixed-template words; fixed zero residual only',
            'witness': expected['projector']['witness'], 'actions_checked': 60}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, default=REFERENCE)
    args = parser.parse_args()
    print(json.dumps(verify(m.load(args.receipt)), indent=2))
