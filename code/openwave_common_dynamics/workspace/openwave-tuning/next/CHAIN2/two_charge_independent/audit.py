"""Independent source-matrix replay of the two-charge phase instrument."""
from pathlib import Path
from functools import lru_cache
import hashlib
import json
import math
import sys

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
PRODUCER = HERE.parent / 'two_charge_response'
INPUT = PRODUCER / 'inputs'
RER = INPUT / 'reverse-engineering-reality'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def scalar(pair):
    from fractions import Fraction
    return float(Fraction(pair[0])) + math.sqrt(5) * float(Fraction(pair[1]))


def source_matrices():
    response = read(RER / 'code/source_selection_model/response.json')
    controls = read(RER / 'code/m1_source_realization/receipt.json')
    generators = []
    for blocks in response['generators']:
        matrix = np.zeros((6, 6), complex)
        for block, values in enumerate(blocks):
            for i, row in enumerate(values):
                for j, cell in enumerate(row):
                    matrix[3*block+i, 3*block+j] = scalar(cell[:2]) + 1j*scalar(cell[2:])
        generators.append(matrix)
    coefficients = {r['name']: np.array([scalar(x) for x in r['coefficients']])
                    for r in controls['evidence']['controls']['rows']}
    selected = {'X': coefficients['imaginary_01'], 'Y': coefficients['real_01'],
                'Z': coefficients['phase_0']-coefficients['phase_1'],
                'I': coefficients['phase_0']+coefficients['phase_1'],
                'J': coefficients['relative_phase']}
    selected['H'] = (selected['X']+selected['Z'])/math.sqrt(2)
    return {kind: sum(c*g for c, g in zip(row, generators))
            for kind, row in selected.items()}


GENERATORS = source_matrices()


@lru_cache(None)
def pulse(kind, angle):
    return expm(float(angle) * GENERATORS[kind])


@lru_cache(None)
def indices(bits, targets):
    # Independent integer bit gather; no producer tensor-axis implementation.
    all_indices = np.arange(1 << bits)
    mask = sum(1 << q for q in targets)
    bases = all_indices[(all_indices & mask) == 0]
    return np.array([bases + sum(((local >> (len(targets)-1-j)) & 1) << q
                                for j, q in enumerate(targets))
                     for local in range(1 << len(targets))])


def execute(words, bits, state):
    state = np.asarray(state, complex).copy()
    vector = state.ndim == 1
    state = state.reshape(1 << bits, -1)
    for word in words:
        targets = tuple(word['targets'])
        full = pulse(word['kind'], word['angle'])
        dimension = 1 << len(targets)
        assert np.linalg.norm(full[dimension:, :dimension]) < 2e-12
        places = indices(bits, targets)
        gathered = state[places]
        state[places] = np.einsum('ab,bcr->acr', full[:dimension, :dimension], gathered)
    return state[:, 0] if vector else state


def main():
    files = [PRODUCER/x for x in ('EDGE_ECHO.json', 'FULL_V_ECHO.json',
                                  'INPUT_SNAPSHOT.json', 'IMPORT_ANNEX.json')]
    files += [RER/'code/source_selection_model/response.json',
              RER/'code/m1_source_realization/receipt.json',
              INPUT/'openwave-tuning/common_dynamics/vector_field/FIELD_PROGRAM.json']
    before = {str(p.relative_to(HERE.parent)): sha(p) for p in files}
    edge = read(PRODUCER/'EDGE_ECHO.json')
    full = read(PRODUCER/'FULL_V_ECHO.json')
    alpha = edge['alpha']
    assert alpha == .17 and full['beta_face0'] == .41
    inputs = [((1 << left) | (1 << (2+right)) | (1 << 4))
              for left in (0, 1) for right in (0, 1)]
    columns = np.zeros((2048, 4), complex)
    columns[inputs, np.arange(4)] = 1
    phase = np.array([np.exp(-.5j*alpha*(1+left-right)**2)
                      for left in (0, 1) for right in (0, 1)])
    echo_word = edge['transport_word']+edge['electric_word']+edge['inverse_transport_word']
    actual = execute(echo_word, 11, columns)
    column_error = float(np.linalg.norm(actual-columns*phase))
    assert column_error < 3e-11
    blank = np.zeros(2048, complex); blank[0] = 1
    prepared = execute(edge['source_preparation_word'], 11, blank)
    final = execute(echo_word+edge['rephase_word']+edge['analysis_word'], 11, prepared)
    expected = np.array([(5+3*math.cos(alpha))/8]+[(1-math.cos(alpha))/8]*3)
    measured = abs(final[inputs])**2
    assert np.max(abs(measured-expected)) < 3e-11
    assert abs(np.vdot(final, final).real-1) < 3e-11
    # Every reference column is retained; Hilbert-Schmidt error upper-bounds
    # operator error on arbitrary input/reference supported on these columns.
    mutations = {}
    no_inverse = execute(edge['transport_word']+edge['electric_word'], 11, columns)
    mutations['omit_inverse_transport'] = float(np.linalg.norm(no_inverse-columns*phase)) > .5
    additive = phase.copy(); additive[3] = phase[1]*phase[2]/phase[0]
    separate = execute(edge['rephase_word']+edge['analysis_word'], 11, columns@additive/2)
    mutations['independent_electric_fields'] = abs(abs(separate[inputs[0]])**2-1) < 3e-11
    mutations['connected_phase_disappears_in_additive_control'] = abs(additive[0]*additive[3]/(additive[1]*additive[2])-1) < 1e-14
    mutations['nonzero_connected_phase_in_actual_source'] = abs(phase[0]*phase[3]/(phase[1]*phase[2])-1) > .1
    dephased = execute(echo_word+edge['rephase_word']+edge['analysis_word'], 11, columns)
    dephased_law = (abs(dephased[inputs])**2).sum(axis=1)/4
    mutations['dephase_prepared_qubits'] = float(np.abs(dephased_law-expected).sum()/2) > .7
    field_program = read(files[-1])
    seed = np.zeros(16, complex); seed[0] = 1
    field_actual = execute(field_program['plaquette_core_word']+3*field_program['electric_word'], 4, seed)
    # Independent dense circulant magnetic source on the full cycle orbit.
    shift = np.roll(np.eye(16), 1, axis=0)
    magnetic = expm(.5j*.41*(shift+shift.T))[:, 0]
    centered = np.array([q if q < 8 else q-16 for q in range(16)])
    field_expected = np.exp(-.5j*.17*3*centered**2)*magnetic
    field_error = float(np.linalg.norm(field_actual-field_expected))
    assert field_error < 3e-11
    full_expected = np.outer(abs(field_expected)**2, expected).reshape(-1)
    retained = np.array([r['probability'] for r in full['all64_joint_event_field_outcomes']])
    assert len(retained) == 64 and np.max(abs(retained-full_expected)) < 3e-13
    mutations['drop_nonzero_plaquette_labels'] = 1-float(abs(field_expected[0])**2) > .05
    assert all(mutations.values())
    maximum_error = max(r['full_native_TV_upper'] for r in full['native_complete_contracts'])
    signal = float(3*(1-math.cos(alpha))/8)
    assert signal > 100*maximum_error
    assert before == {str(p.relative_to(HERE.parent)): sha(p) for p in files}
    result = {'schema':'oph.two_charge.independent_source_matrix_audit.v1',
              'status':'PASS_SCOPED_SOURCE_PULSE_AND_FACTORIZED_RECORD_LAW',
              'column_operator_upper':column_error, 'common_field_factor_error':field_error,
              'record_law':measured.tolist(), 'connected_record_TV':signal,
              'all64_records_reproduced':True, 'reference_columns_retained':4,
              'mutation_controls':{k:bool(v) for k,v in mutations.items()}, 'source_pins':before,
              'scope':'Independent M6 matrix/bit-index pulse replay and 64-outcome factorized law. Full-source graph embedding and complete apparatus bound are separately audited dependencies; no independent 3203-qubit array execution.',
              'official_rows_changed':[], 'force_or_Coulomb_claim':False}
    (HERE/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_pins','scope')},indent=2),flush=True)


if __name__ == '__main__':
    main()
