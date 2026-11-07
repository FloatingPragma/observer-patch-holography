"""Reproduce the bounded native rank-one effect synthesis certificate."""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np
from tools import strict_json
from . import intervals as iv, synthesis as s
from model.native_effect_coverage import packet as coverage_packet

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_FILES = tuple('model/native_effect_synthesis/'+name for name in
                    ('__init__.py','intervals.py','synthesis.py','packet.py','test_synthesis.py','PROOF.md','README.md'))+(
    'Lean/EventAlgebra/OperationalBornBridge.lean','tools/strict_json.py','model/requirements.txt')
# Complete loaded local import closure in addition to the owning detector's
# declared source/data bindings. These are code inputs, never documentation
# notes or a dynamically changing list of modules from the caller's session.
IMPORT_FILES = (
    'model/anharmonic_transfer/__init__.py','model/common_source/composition.py',
    'model/common_source/operational.py','model/common_source/thermal_clock.py',
    'model/federated_observer/__init__.py','model/laboratory_clock/__init__.py',
    'model/laboratory_clock/event_geometry.py','model/local_repair_family/__init__.py',
    'model/native_causal_dataflow/__init__.py','model/native_interface/pointer.py',
    'model/native_modular_reference/__init__.py','model/native_modular_reference/cap_local_clock.py',
    'model/native_modular_reference/compatibility.py','model/native_modular_reference/completed_clock.py',
    'model/native_modular_reference/construction.py','model/native_modular_reference/full_flow.py',
    'model/native_modular_reference/local_history.py','model/native_program/__init__.py',
    'model/observer/record_recovery.py','model/observer/source_discrimination/__init__.py',
    'model/observer_reconstruction/__init__.py','model/observer_reconstruction/native_geometry.py',
    'model/predictions/protocol.py','model/predictions/repair_response.py',
    'model/primitive_observer/__init__.py','model/relativistic_scalar/__init__.py',
    'model/relativistic_scalar/lattice.py','model/repair_identification/__init__.py',
    'model/repair_identification/core.py','model/repair_scheduler/__init__.py',
    'model/scalar_observations/__init__.py','model/scalar_observations/experimental_readout/__init__.py',
    'model/scalar_observations/lattice.py','model/scalar_observations/source.py')


def source_bindings():
    paths = set(SOURCE_FILES)|set(IMPORT_FILES)|set(coverage_packet.source_bindings())
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def examples():
    sparse = tuple((F(3) if i == 0 else F(), F(4) if i == 65 else F()) for i in range(66))
    dense = tuple((F(i%5-2), F((3*i)%7-3)) for i in range(66))
    return (('complex_pair', sparse, F(3,5)), ('all_source_coordinates', dense, F(7,13)))


def word_digest(word):
    value = s.interpreter.source_word(word, 9)
    return hashlib.sha256(s.interpreter.canonical(value).encode()).hexdigest()


def produce():
    before = source_bindings(); rows = []
    s.certify_primitive.cache_clear()
    short = s.pair_word(0,65); checked = s.certify_word(short)
    row = s.selected_row(short,0)[list(s.tomo.codes())]
    desired = np.zeros(66,complex); desired[0]=.6; desired[65]=.8
    if np.linalg.norm(row-desired) > 1e-8:
        raise ArithmeticError('the actual inverse analyzer lost the unequal-weight effect')
    short_data = {'source_indices':[0,65],'target_pair_effect':[['9/25','12/25'],['12/25','16/25']],
                  'source_word_sha256':word_digest(short),'distinct_primitive_requests_checked':checked,
                  'complete_bounds':s.bounds(len(short),general=False),'outcomes':s.outcome_catalogue(short),
                  'actual_selected_row_within_diagnostic_gate':True}
    for name, vector, strength in examples():
        word = s.rank_one_word(vector,strength); checked = s.certify_word(word)
        norm = sum((x*x+y*y for x,y in vector),F())
        wanted = np.array([complex(float(x),-float(y)) for x,y in vector])*np.sqrt(float(strength/norm))
        row = s.selected_row(word)[list(s.tomo.codes())]
        if np.linalg.norm(row-wanted) > 1e-8:
            raise ArithmeticError('the literal negative source word has the wrong selected row')
        admitted = s.interpreter.source_word(word,9)
        stopped = s.interpreter.execute_word(admitted,np.zeros((1,1)))
        if stopped['status']!='RESOURCE_STOP' or stopped['pulses_executed']!=0:
            raise ArithmeticError('the generic word evaluator boundary changed')
        rows.append({'name':name,'rational_source_vector':[[str(x),str(y)] for x,y in vector],
                     'strength':str(strength),'source_word_sha256':word_digest(word),
                     'distinct_primitive_requests_checked':checked,'resources':admitted['declared_resources'],
                     'complete_bounds':s.bounds(len(word),general=True),'outcomes':s.outcome_catalogue(word),
                     'actual_selected_row_within_diagnostic_gate':True,
                     'generic_evaluator_status':stopped['status']})
    result = {'schema':'oph.native_effect_synthesis.v1','unequal_pair':short_data,'general_examples':rows,
              'construction':{'work_qubits':9,'work_dimension':512,'source_dimension':66,
                              'extra_qubits_or_new_primitive_services':0,'angle_grid_bits':iv.ANGLE_BITS,
                              'trigonometric_enclosure_bits':iv.BITS,'input_rational_bits_max':s.MAX_RATIONAL_BITS,
                              'amplitude_tree_pulses':12291,'full_diagonal_phase_pulses':11270,
                              'final_zero_to_source_pulses_max':18,'total_word_pulses_max':s.MAX_PULSES,
                              'source_word_admission_pulse_limit':s.interpreter.MAX_WORD,
                              'literal_selected_row_diagnostic_gate':'1/100000000',
                              'floating_residual_digits_serialized':False},
              'scope':{'fixed_precision_rank_one_effects_constructed':True,
                       'all_positive_effects_or_general_Naimark_isometry_compiled':False,
                       'universal_coexistent_refinement_or_Born_premises_closed':False,
                       'zero_error_or_arbitrarily_small_error_at_fixed_budget':False,
                       'rounded_maps_asserted_exactly_trace_preserving':False,
                       'Luders_update_inferred_from_effect':False,
                       'generic_host_limits_raised':False,'physical_QFT_or_particle_identification':False},
              'source_bindings':before}
    if source_bindings()!=before:
        raise ArithmeticError('source synthesis inputs changed during replay')
    return result


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path)))!=canonical(produce()):
        raise ValueError('native effect synthesis receipt differs from source replay')
    return True


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()));print('Native effect synthesis receipt written.')
    else:
        verify();print('Native effect synthesis receipt verified.')
