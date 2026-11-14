"""Size-three free radiative source, both field quadratures and common clock."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import sys

sys.dont_write_bytecode = True
import numpy as np
import frame
import providers
import walsh

HERE = Path(__file__).resolve().parent
nf = walsh.nf; compiler = nf.compiler; cert = nf.cert; iv = nf.iv
sys.path.insert(0, str(HERE.parent/'clock'))
import phase_clock
from model.native_clock import configuration
from model.observer.source_discrimination import circuit
from model.native_interface import pointer


@lru_cache(None)
def independent_executor():
    # Audited response-matrix executor, independent of the compiler executor.
    location = HERE.parent/'photon_audit/audit.py'
    sys.path.insert(0, str(location.parent))
    spec = importlib.util.spec_from_file_location('size3_opposing_matrix_executor', location)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module.execute


def targets(time=F(1)):
    if type(time) is not F or not 0 < time <= 8:
        raise ValueError('exact Fraction time in (0,8] required')
    eigenvalues, certificate = frame.spectrum()
    frequencies = [iv.sqrt_fraction(F(3, 2)*x) for x in eigenvalues]
    E0 = sum(frequencies, iv.Interval.point(0))/2
    energies = [E0]+[E0+x for x in frequencies]+[iv.Interval.point(0)]*377
    assert len(energies) == 1024
    return [-time*x for x in energies], frequencies, E0


def analysis(setting, bit):
    if setting not in ('X', 'Y'):
        raise ValueError('field and clock quadratures must be X or Y')
    return ([] if setting == 'X' else [compiler.pulse('Z', math.pi/4, [bit])])+compiler.remap(nf.H(0), [bit])


def full_word(diagonal, field_setting='X', clock_setting='X', time=F(1)):
    # A vacuum/one-photon superposition; all other modes start in vacuum.
    return nf.H(0)+diagonal+analysis(field_setting, 0)+compiler.remap(phase_clock.local_words(1, clock_setting, time), [10])


def probability_law(vector, effects):
    assert max(abs(effects[b][0, 1]) for b in (0, 1)) < 1e-13
    transition = np.array([[effects[b][x, x].real for x in (0, 1)] for b in (0, 1)])
    probabilities = abs(vector)**2
    for bit in range(11):
        base = np.arange(2048); lower = base[(base & (1 << bit)) == 0]
        positions = np.array([lower, lower+(1 << bit)])
        probabilities[positions] = transition@probabilities[positions]
    return probabilities


def run():
    before = providers.verify()
    frame_data, frozen_frame = frame.data()
    walsh_checks = walsh.controls()
    phase_targets, frequencies, E0 = targets()
    diagonal, phase_error = walsh.compile_intervals(phase_targets)
    assert len(diagonal) == 25606
    phase_mid = np.array([float(F(x.lo+x.hi, 2*iv.SCALE)) for x in phase_targets])
    execute = independent_executor()
    # Complete phase certificate is supplemented by physical/padding columns.
    selected = [0, 1, 2, 323, 646, 647, 1023]
    columns = np.eye(1024, dtype=complex)[:, selected]
    actual_columns = execute(diagonal, 10, columns)
    column_residual = float(np.linalg.norm(actual_columns-columns*np.exp(1j*phase_mid[selected])))
    assert column_residual < float(phase_error)+2e-11
    blank = np.zeros(1024, complex); blank[0] = 1
    prepared = execute(nf.H(0), 10, blank)[:, 0]
    evolved = execute(diagonal, 10, prepared)[:, 0]
    expected_evolved = np.zeros(1024, complex)
    expected_evolved[:2] = np.exp(1j*phase_mid[:2])/math.sqrt(2)
    assert np.linalg.norm(evolved-expected_evolved) < float(phase_error)+2e-11
    p0 = min(map(F, configuration.success_probabilities(4, 32, full=True, mode='reset', target=circuit.decode(0))))
    pr = min(map(F, configuration.success_probabilities(4, 32, full=True, mode='reset', target=circuit.decode(pointer.READY_CODE))))
    reset = 11*(1-p0)+11*(1-pr)
    htemplate = [('I', (0,), cert._pi(-1, 2)), ('H', (0,), cert._pi(1, 2))]
    herror = nf.root_certificate.exact_template_bound(nf.H(0), htemplate)
    extra_y = nf.root_certificate.exact_template_bound([compiler.pulse('Z', math.pi/4, [0])], [('Z', (0,), cert._pi(1, 4))])
    instruments, _ = phase_clock.actual_native_pointer()
    geometry = F(frozen_frame['certificate']['cap1_H_error_including_vacuum'])
    rows = []; programs = {}
    for field_setting in ('X', 'Y'):
        photon = execute(analysis(field_setting, 0), 10, evolved)[:, 0]
        for clock_setting in ('X', 'Y'):
            clock_word = phase_clock.local_words(1, clock_setting, F(1))
            clock = execute(clock_word, 1, np.array([1., 0.]))[:, 0]
            v = np.kron(clock, photon)
            # Disjoint source and clock factors prove this equals full11bit word.
            expected_photon = expected_evolved.copy()
            if field_setting == 'Y':
                expected_photon[0] *= np.exp(1j*math.pi/4); expected_photon[1] *= np.exp(-1j*math.pi/4)
            expected_photon[:2] = np.array([expected_photon[0]+expected_photon[1], expected_photon[0]-expected_photon[1]])/math.sqrt(2)
            angle = -1/16+(math.pi/4 if clock_setting == 'Y' else 0)
            expected = np.kron([math.cos(angle), 1j*math.sin(angle)], expected_photon)
            residual = float(np.linalg.norm(v-expected))
            native = phase_error+2*herror+(extra_y if field_setting == 'Y' else 0)
            native += phase_clock.exact_word_bound(1, clock_setting, F(1))+10*phase_clock.pointer_bound()
            complete = reset+native+geometry
            assert residual < float(native)+2e-11
            laws = {}
            for law, branches in instruments.items():
                effects = [np.zeros((2, 2), complex), np.zeros((2, 2), complex)]; fine = []
                for metadata, kraus in branches:
                    effect = sum((K.conj().T@K for K in kraus), np.zeros((2, 2), complex))
                    effects[metadata['bit']] += effect
                    fine.append({'record':metadata, 'effect_real':effect.real.tolist(), 'effect_imag':effect.imag.tolist()})
                probabilities = probability_law(v, effects)
                assert min(probabilities) > -1e-13 and abs(sum(probabilities)-np.vdot(v, v).real) < 3e-11
                contrast = sum(p*(-1 if j & 1 else 1) for j, p in enumerate(probabilities))
                frequency = float(F(frequencies[0].lo+frequencies[0].hi, 2*iv.SCALE))
                target = math.cos(frequency) if field_setting == 'X' else -math.sin(frequency)
                assert abs(contrast-target) < float(native)+3e-11
                laws[law] = {'all2048_joint_probabilities':probabilities.tolist(),
                    'radiative_quadrature_contrast':float(contrast), 'fine_one_bit_effects':fine,
                    'success_normalization':False,
                    'fine_rule':'unnormalized v* tensor(11 retained fine effects) v; first reset/export failures are separate outputs'}
            name = field_setting+clock_setting
            program = full_word(diagonal, field_setting, clock_setting)
            programs[name] = program
            rows.append({'field_quadrature':field_setting, 'clock_quadrature':clock_setting,
                'native_pulses':len(program), 'native_word_and_pointer_error_exact':str(native),
                'geometry_H_error_exact':str(geometry), 'reset_failure_bound_exact':str(reset),
                'complete_record_TV_bound_exact':str(complete), 'complete_record_TV_bound':float(complete),
                'native_state_residual':residual, 'laws':laws})
    # Both field quadratures identify the mode phase in the declared [0,pi] band.
    frequency = float(F(frequencies[0].lo+frequencies[0].hi, 2*iv.SCALE))
    assert 0 < frequency < math.pi
    X = rows[0]['laws']['C']['radiative_quadrature_contrast']
    Y = rows[2]['laws']['C']['radiative_quadrature_contrast']
    recovered = math.atan2(-Y, X)
    assert abs(recovered-frequency) < 3e-10
    no_frequency = max(abs(X-1), abs(Y))
    assert no_frequency > 20*max(r['complete_record_TV_bound'] for r in rows)
    erased = execute([compiler.pulse('Z', -frequency/2, [0])], 10, evolved)[:, 0]
    erasure_contrasts = {}
    for setting in ('X', 'Y'):
        measured = execute(analysis(setting, 0), 10, erased)[:, 0]
        erasure_contrasts[setting] = float(sum((-1 if j & 1 else 1)*abs(z)**2 for j, z in enumerate(measured)))
    assert abs(erasure_contrasts['X']-1) < 3e-10 and abs(erasure_contrasts['Y']) < 3e-10
    assert before == providers.verify()
    program_data = {'schema':'oph.quadratic.size3.native_program.v1', 'diagonal_word':diagonal,
        'experiments':programs, 'source_bits':list(range(10)), 'clock_bit':10,
        'pointer_auxiliary_bits':[11, 12], 'blank_resets':11, 'pointer_ready_resets':11,
        'physical_states':647, 'retained_padding':377,
        'encoder':'abstract exact coexact Fock frame U=polar(K Q Lambda^-1); no local electric-register encoder asserted',
        'readout':'X and Y quadratures of the declared radiative mode plus common-clock X/Y, all basis/failure outcomes retained'}
    (HERE/'PROGRAM.json').write_text(json.dumps(program_data, separators=(',', ':'))+'\n')
    result = {'schema':'oph.quadratic.size3.native_member.v1',
        'candidate_id':'native-local12-dirac-noncompact-quadratic-v4-proposed',
        'checks':{'exact_Walsh_controls':all(walsh_checks.values()), 'all1024phase_targets_certified':True,
            'native_response_matrix_columns_and_clock':True, 'all646modes_and377padding_retained':True,
            'both_radiative_field_quadratures':True, 'both_complete_native_pointer_laws':True,
            'all2048joint_outcomes_and_failures':True, 'source_vacuum_scalar_retained':True,
            'common_clock_and_exact_geometry_error_composed':True, 'frequency_deletion_fails':True,
            'producer_bytes_unchanged':True},
        'Walsh_checks':walsh_checks, 'diagonal_pulses':len(diagonal),
        'all_phase_operator_error_exact':str(phase_error), 'selected_native_column_residual':column_residual,
        'first_mode_frequency':frequency, 'native_recovered_phase_in_declared_band':recovered,
        'frequency_deletion_contrast_change':float(no_frequency), 'experiments':rows,
        'executed_native_counterphase_erasure_contrasts':erasure_contrasts,
        'vacuum_energy_interval':E0.json(), 'program_sha256':providers.sha(HERE/'PROGRAM.json'),
        'frame_certificate_sha256':providers.sha(HERE/'CERTIFICATE.json'), 'source_sha256':before,
        'failure_contract':'Rank-one conclusive reset theorem; retain first timeout and all native export complements, private labels, fine pointers and padding. No success renormalization.',
        'scope':'Free radiative cap1 member on actual size3 field. Canonical radiative matrix elements and traced harmonic reference are specified separately. No cap1 interacting Weyl accuracy.',
        'official_rows_changed':[]}
    (HERE/'NATIVE.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'checks':result['checks'], 'phase':frequency, 'TV':max(r['complete_record_TV_bound'] for r in rows),
                      'pulses':[r['native_pulses'] for r in rows], 'column_residual':column_residual}), flush=True)


if __name__ == '__main__': run()
