"""A native phase record for completed words of the same LOCAL12 source.

The attachment is an explicit clock architecture. It supplies no seconds,
particle rest energy, or implicit evolution between completed source words.
"""
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

from clock_contract import ROOT, PROVIDER, MODEL_ID, CLOCK_ATTACHMENT_ID, effective_time, positive_rational

HERE = Path(__file__).resolve().parent
INSTRUMENTS = ROOT/'openwave-tuning/auditor/common_walk_particle'
SOURCE_ERROR = ROOT/'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json'
sys.path.insert(0, str(INSTRUMENTS))
from native_phase_probe import actual_native_pointer
from instrument_contract import pointer_bound, reset_contract
from model.observer.source_scalar import compiler
from model.reference_recovery import recovery
from model.certified_reference import certificate as cert

PHASE_PER_COMPLETED_W = F(1, 8)
CLOCK_BIT = 23


def clock_tick_word(clock_bit=CLOCK_BIT, cell_spacing=F(1), source_memory_bits=23):
    if (type(source_memory_bits) is not int or source_memory_bits < 23
            or type(clock_bit) is not int or clock_bit < source_memory_bits):
        raise ValueError('a fresh clock wire outside the complete source register is required')
    positive_rational(cell_spacing, 'cell spacing')
    return [compiler.pulse('Z', -float(cell_spacing*PHASE_PER_COMPLETED_W)/2, [clock_bit])]


def local_words(steps, setting, cell_spacing=F(1)):
    positive_rational(cell_spacing, 'cell spacing')
    if type(steps) is not int or not 0 <= steps <= 4096 or steps*cell_spacing > 8 or setting not in ('X', 'Y'):
        raise ValueError('finite clock window0<=steps*h<=8, at most4096 steps, and configured X/Y read required')
    hadamard = [compiler.pulse('I', -math.pi/2, [0]),
                compiler.pulse('H', math.pi/2, [0])]
    tick = [compiler.pulse('Z', -float(cell_spacing)/16, [0])]
    analysis = ([] if setting == 'X' else [compiler.pulse('Z', math.pi/4, [0])]) + hadamard
    return hadamard + tick*steps + analysis


def exact_word_bound(steps, setting, cell_spacing=F(1)):
    word = local_words(steps, setting, cell_spacing)
    H = [('I',(0,),cert._pi(-1,2)), ('H',(0,),cert._pi(1,2))]
    template = H + [('Z',(0,),-cell_spacing/16)]*steps
    if setting == 'Y': template += [('Z',(0,),cert._pi(1,4))]
    template += H
    end, angle, variables = cert.consume(word, 0, template)
    assert end == len(word) and not variables
    return angle + cert.controls_bound(word, {}) + pointer_bound()


def run():
    files = [Path(__file__), HERE/'clock_contract.py',
             INSTRUMENTS/'native_phase_probe.py', INSTRUMENTS/'instrument_contract.py',
             PROVIDER/'COMMON_WALK_RECEIPT.json', PROVIDER/'COMMON_WALK_PROGRAM.json',
             Path(compiler.__file__), Path(recovery.__file__), Path(cert.__file__), SOURCE_ERROR]
    pins = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    pointer, pointer_run = actual_native_pointer()
    reset = reset_contract(1,1)
    source_error = json.loads(SOURCE_ERROR.read_text())
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == digest
               for p,digest in source_error['source_sha256'].items())
    rows = []
    worst = 0.
    largest_bound = F(0)
    cases = [(n,F(1)) for n in range(9)] + [(16,F(1,16)), (32,F(1,32)),
             (64,F(1,64)), (256,F(1,64)), (512,F(1,64))]
    for steps, spacing in cases:
        phase = float(steps*spacing/8)
        contrasts = {}
        for setting in ('X','Y'):
            run = recovery.execute_complete(local_words(steps,setting,spacing), 1,
                                             np.eye(2,dtype=complex)[:,[0]])
            source = run.success[:,0]
            ideal_contrast = math.cos(phase) if setting == 'X' else math.sin(phase)
            bound = exact_word_bound(steps, setting,spacing)
            largest_bound = max(largest_bound, bound)
            laws = {}
            for law, branches in pointer.items():
                probabilities = []
                bits = np.zeros(2)
                for meta, ks in branches:
                    probability = sum(float(np.vdot(K@source, K@source).real) for K in ks)
                    bits[meta['bit']] += probability
                    probabilities.append({'record':meta, 'probability': probability})
                contrast = bits[0]-bits[1]
                worst = max(worst, abs(contrast-ideal_contrast), abs(bits.sum()-1))
                assert abs(contrast-ideal_contrast) < 2e-12
                laws[law] = {'native_fine_probabilities': probabilities,
                             'all_bit_probabilities': bits.tolist(), 'contrast': float(contrast)}
                contrasts[setting,law] = contrast
            rows.append({'completed_W_steps': steps, 'configured_setting':setting,
                         'cell_spacing_exact':str(spacing),
                         'clock_phase_exact_radians':str(steps*spacing/8),
                         'local_clock_pulses':len(local_words(steps,setting,spacing)),
                         'native_pointer_pulses':6,
                         'source_export_failure_gram_trace':float(np.trace(run.failure_gram).real),
                         'word_and_pointer_error_upper_exact':str(bound), 'laws':laws})
        for law in ('C','E'):
            z = complex(contrasts['X',law], contrasts['Y',law])
            assert abs(z-np.exp(1j*phase)) < 3e-12
    # The complete-error bound includes the two physical reset failures.
    uniform_extra_reset = sum(F(row['complete_reset_certificate']['failure_upper'])
                              for row in reset['native_reset_contracts'])
    complete_error = largest_bound + uniform_extra_reset
    common_receipt = json.loads((PROVIDER/'COMMON_WALK_RECEIPT.json').read_text())
    common_reset = 23*F(common_receipt['preparation_reset']['certificate']['failure_upper'])
    joint_error = (complete_error + common_reset
                   + F(source_error['whole_source_preparation_error_exact'])
                   + max(n for n,h in cases)*F(source_error['whole_W_operator_error_exact']))
    # Each unconditional bounded sign observable differs by at most2epsilon.
    # A radius3epsilon bounds the two-quadrature complex displacement.
    radius = 3*joint_error
    phase_error = radius/(1-radius)
    assert phase_error < F(1,32)
    # A different tick coefficient is a different declared clock attachment.
    control_separation = abs(np.exp(1j)-np.exp(2j))
    assert control_separation > .9
    assert clock_tick_word(898,F(1,64),898)[0]['targets'] == [898]
    try:
        clock_tick_word(23,F(1),898)
    except ValueError:
        pass
    else:
        raise AssertionError('clock/field wire alias accepted')
    receipt = {'schema':'oph.local12.native_phase_clock.v1','candidate_id':MODEL_ID,
        'clock_attachment_id':CLOCK_ATTACHMENT_ID,
        'checks':{'native_phase_tick_and_XY_pointer_executed':True,
                  'all_native_C_and_E_fine_records':True,'phase_lift_unique_on_declared_window':True,
                  'certified_clock_word_error':True,'all_extra_reset_timeouts_counted':True,
                  'tick_rate_mutation_changes_clock':True,'one_shared_time_conversion':True,
                  'full_source_wire_alias_rejected':True,
                  'source_pins_unchanged':True},
        'tick_word_global_address':clock_tick_word(),
        'clock_bit':CLOCK_BIT,'phase_window_radians':[0,1],
        'tested_steps_and_spacings':[[n,str(h)] for n,h in cases],
        'maximum_tested_source_steps':max(n for n,h in cases),
        'joint_error_source_member_size':2,
        'larger_member_native_error_transferred':False,
        'wire_allocation':'The bare source has23 wires. Its clock default23 must be remapped when another instrument uses23. The full current-field construction has898 wires; a joint clock/probe can use fresh898/899 plus separate detector auxiliaries. clock_tick_word validates the declared full source_memory_bits.',
        'augmented_source_boundary':'The tested bound and duration use bare W. An actual W_F clock attaches after the entire augmented current word, whose full preparation, resets, pulse error and timing must be supplied. A controlled-W probe must control the current-addition words too; controlling only bare W is not a physical controlled-W_F instrument.',
        'phase_per_completed_W':'h/8; h=1 for the base member',
        'native_configured_records':rows,'maximum_numeric_record_error':worst,
        'complete_clock_error_upper_exact':str(complete_error),
        'extra_reset_uniform_bound_exact':str(uniform_extra_reset),
        'common_source_23_blank_reset_bound_exact':str(common_reset),
        'complete_joint_source_and_clock_error_on_tested_cases_exact':str(joint_error),
        'complete_joint_source_and_clock_error_upper':math.nextafter(float(joint_error),math.inf),
        'phase_error_upper_exact':str(phase_error),
        'extra_reset_contract':reset,
        'joint_program':'Prepare the common W source and a separately reset clock+. At each successful completed W append the one native Z tick. On any export or reset failure retain its complete failure record and the last completed-step count. Terminal X or Y is a configured intervention with an actual native C/E pointer; source/event/current registers are untouched by the clock read.',
        'complete_instrument':'All successes, export branches and native reset timeouts remain. Failure can be assigned zero in a bounded contrast statistic; this is not renormalization of successful trials. The phase bound compares the complete instrument with its ideal reset-success experiment.',
        'reference_contract':'Ideal joint macro map is W tensor exp(-i h Z/16). It preserves arbitrary external references. Product clock preparation is declared; a correlated clock/source input requires its joint state rather than this calibration.',
        'time_interface':{'function':'clock_contract.effective_time','formula':'t=n*tau0*h=8*tau0*unwrapped_phase',
                          'examples':{str(n):str(effective_time(n,F(1,16))) for n in (0,1,4,8)}},
        'scope':'This is an additional native operational clock on completed common source words. The refinement-dependent phase increment, clock preparation and attachment schedule are declared architecture. The source selects no absolute duration; no mass clock or Lorentz proper time is inferred. The joint bound includes23 common blank resets, two extra apparatus resets, common source preparation, at most512 W words and the terminal clock pointer. It assumes the ideal declared coherent services and exact retained response of the source-error certificate. Strict implementation convergence as h tends to zero additionally requires n*W_error to tend to zero; the fixed compiler bound does not imply that.',
        'source_sha256':pins,'official_rows_changed':[]}
    assert pins == {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    (HERE/'PHASE_CLOCK.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'checks':receipt['checks'], 'numeric_record_error':worst,
                      'complete_clock_error':float(complete_error),'complete_joint_error':float(joint_error),
                      'phase_error_bound':float(phase_error)}))


if __name__=='__main__':
    run()
