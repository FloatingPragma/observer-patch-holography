"""Exact size-three rank and a rejected eight-bit compiler transfer."""
from pathlib import Path
import json
import sys

sys.dont_write_bytecode = True
import numpy as np
import frozen

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'vector_field'))
import cell_complex as cc


def modular_rank(matrix, prime=101):
    value = np.asarray(matrix, dtype=np.int64) % prime
    row = 0
    for col in range(value.shape[1]):
        nonzero = np.flatnonzero(value[row:, col])
        if not len(nonzero): continue
        pivot = row+int(nonzero[0])
        value[[row, pivot]] = value[[pivot, row]]
        value[row, col:] = value[row, col:]*pow(int(value[row, col]), -1, prime) % prime
        value[row+1:, col:] = (value[row+1:, col:]-value[row+1:, [col]]*value[[row], col:]) % prime
        row += 1
        if row == value.shape[0]: break
    return row


def run():
    before = frozen.verify_inputs()
    graph, B, A, periods, winding, labels = cc.build(3)
    assert B.shape == (324, 972) and A.shape == (972, 702)
    # build() verifies BA=0, periods*A=0 and periods*winding=I3 exactly.
    adjacency = [[] for _ in range(324)]
    for a, b in graph.seams: adjacency[a].append(b); adjacency[b].append(a)
    seen = {0}; queue = [0]
    for a in queue:
        for b in adjacency[a]:
            if b not in seen: seen.add(b); queue.append(b)
    assert len(seen) == 324
    # Incidence rank323 and three independent periods bound rank(A)<=646.
    # Modular rank gives an exact rational lower bound in the other direction.
    rank = modular_rank(A)
    assert rank == 646 == 972-(324-1)-3
    compiler = frozen._native().compiler
    rejected = False
    try: compiler.compile_diagonal(np.zeros(1024))
    except ValueError: rejected = True
    assert rejected
    count = 0; maxwire = 0
    for mask in range(1024):
        label = ''.join('Z' if mask >> j & 1 else 'I' for j in range(10))
        template = compiler.pauli_rotation(label, 0.)
        count += len(template); maxwire = max(maxwire, *(max(p['targets']) for p in template))
    assert count == (3*10-5)*1024+6 == 25606 and maxwire == 9
    assert before == frozen.verify_inputs()
    result = {'schema':'oph.photon_audit.size3_domain.v1',
              'checks':{'actual_size3_connected_graph':True,
                        'exact_mod101_rank_and_period_upper_bound':True,
                        'eight_bit_diagonal_compiler_rejects_size3_cap':True,
                        'ten_bit_Pauli_template_resource_count':True,
                        'producer_bytes_unchanged':True},
              'vertices':324, 'edges':972, 'faces':702, 'incidence_rank':323,
              'coexact_rank':646, 'harmonic_dimension':3,
              'cap1_physical_states':647, 'spectral_bits':10, 'padding_states':377,
              'clock_bits':1, 'joint_outcomes_per_C_E_law':2048,
              'diagonal_template_pulses':count, 'proposed_reset_counts':{'blank':11, 'pointer_ready':11},
              'cap1_vacuum_plus_excitation_error_factor':324,
              'scope':'Exact next-member geometry and compiler-domain/resource control only. No size3 spectrum, phase certificate, native experiment, harmonic discretization or interacting field word is supplied.',
              'official_rows_changed':[]}
    (HERE/'SPATIAL_DOMAIN.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__': run()
