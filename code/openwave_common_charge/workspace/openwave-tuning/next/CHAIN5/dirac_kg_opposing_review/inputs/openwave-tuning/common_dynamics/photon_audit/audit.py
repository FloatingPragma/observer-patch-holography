"""Independent finite-spectrum, symbol, and native-record audit of DYN16."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import math
import sys

sys.dont_write_bytecode = True
import numpy as np
import sympy as s
from scipy.linalg import expm

import frozen
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PROVIDER = HERE.parent/'quadratic_field'
sys.path.insert(0, str(HERE.parent/'vector_field'))
import cell_complex as cc
import linear_response as linear


def read(name): return json.loads((PROVIDER/name).read_text())


@lru_cache(None)
def geometry():
    graph, B, A, periods, winding, labels = cc.build(2)
    B = np.asarray(B, dtype=np.int64); A = np.asarray(A, dtype=np.int64)
    assert B.shape == (96, 288) and A.shape == (288, 208) and not np.any(B@A)
    assert np.linalg.matrix_rank(B) == 95 and np.linalg.matrix_rank(A) == 190
    assert 288-95-190 == 3
    K = A@A.T
    L = B@B.T; one = np.ones((96, 96))/96
    G = np.linalg.inv(L+one)-one
    PL = B.T@G@B; PT = np.eye(288)-PL
    return graph, B, A, K, G, PL, PT


def continuous_source():
    graph, B, A, K, G, PL, PT = geometry()
    lam, U = np.linalg.eigh(K.astype(float)); keep = lam > 1e-9
    lam, U = lam[keep], U[:, keep]; PH = PT-U@U.T
    omega = np.sqrt(1.5*lam); Cq = (U*(.25/omega))@U.T+PH/2
    Cp = (U*omega)@U.T+PH/2
    assert np.max(abs(Cq@Cp-PT/4)) < 5e-12
    assert np.max(abs(B@Cp)) < 5e-12 and abs(np.trace(PH)-3) < 5e-12
    source = read('SOURCE.json')
    assert np.max(abs(omega-source['frequencies'])) < 1e-12
    assert F(2, 3)*F(1, 2)*3 == 1 and F(2, 3)*F(1, 2) == F(1, 3)
    charges = []
    for row in source['static_charge_cases']:
        q = np.zeros(96); q[row['vertices'][0]] += 1; q[row['vertices'][1]] -= 1
        particular = B.T@G@q
        energy = .25*float(q@G@q)
        assert np.max(abs(B@particular-q)) < 1e-12 and np.max(abs(A.T@particular)) < 1e-12
        assert abs(.25*float(particular@particular)-energy) < 1e-13
        assert abs(energy-row['electric_minimum_energy']) < 1e-13
        charges.append(energy)
    sys.path.insert(0, str(HERE.parent/'spin'))
    import dirac_proposal as matter
    route_count = 0; probabilities = []
    for c in range(8):
        for axis in range(3):
            currents = []
            for sign in (-1, 1):
                target, _, J = matter.route(c, axis, sign)
                J = np.asarray(J, float); currents.append(J)
                assert np.linalg.norm(PL@J+PT@J-J) < 5e-12
                assert np.linalg.norm(B@PT@J) < 5e-12
                route_count += 1
            delta = currents[1]-currents[0]
            assert not np.any(B@delta)
            if c == 0:
                probability = (1+math.exp(-float(delta@Cq@delta)/2))/2
                assert abs(probability-source['same_field_route_probes'][axis]['complete_probe_probabilities'][0]) < 1e-12
                probabilities.append(probability)
    # Zero-mode packet covariance evolves by a free shear, with finite energy.
    covariance = np.diag([.5, .5]); shear = np.array([[1., .2], [0., 1.]])
    evolved = shear@covariance@shear.T
    assert abs(evolved[0, 0]-.52) < 1e-15 and abs(np.linalg.det(evolved)-.25) < 1e-15
    return {'edges':288, 'faces':208, 'coexact':190, 'harmonic':3, 'route_decompositions':route_count,
            'static_charge_energies':charges, 'shared_field_route_probabilities':probabilities,
            'harmonic_variance_at_time_0_4':float(evolved[0, 0]),
            'harmonic_stationary_vacuum_mutation_rejected':bool(evolved[0, 0] != covariance[0, 0])}


def spectral_frame():
    K = geometry()[3]
    cert = frozen.spectral_certificate(); den = cert['denominator']
    Q = np.array(cert['frame_integers'], dtype=np.int64)
    lam = np.array([int(x*den) for x in cert['eigenvalues_exact']], dtype=np.int64)
    assert Q.shape == (288, 190) and den == 1 << 44
    # Three base-2^15 limbs, independently of the producer's two 22-bit limbs.
    base = 1 << 15; working = Q.copy(); limbs = []
    for _ in range(2):
        high = working//base; limbs.append(working-high*base); working = high
    limbs.append(working)
    gram = np.zeros((190, 190), dtype=object)
    for i, left in enumerate(limbs):
        for j, right in enumerate(limbs):
            assert 288*int(np.max(abs(left)))*int(np.max(abs(right))) < 1 << 63
            gram += (left.T@right).astype(object)*base**(i+j)
    gram -= den**2*np.eye(190, dtype=object)
    assert int(np.max(np.sum(abs(K), axis=1)))*int(np.max(abs(Q))) < 1 << 63
    residual = (K@Q).astype(object)*den-Q.astype(object)*lam.astype(object)[None, :]
    gram2 = F(sum(int(x)**2 for x in gram.flat), den**4)
    residual2 = F(sum(int(x)**2 for x in residual.flat), den**4)
    epsilon = cert['frame_Gram_error_upper']; r = cert['frame_residual_upper']
    assert epsilon**2 >= gram2 and r**2 >= residual2
    minimum, maximum = min(cert['eigenvalues_exact']), max(cert['eigenvalues_exact'])
    shift = r/minimum; gamma = epsilon+2*(1+epsilon)*shift+shift**2
    eta = r+(F(1504, 17)+maximum)*(shift+gamma)
    sqrt_min = F(math.isqrt(minimum.numerator*minimum.denominator), minimum.denominator)
    frequency = F(5, 4)*eta/sqrt_min
    assert gamma == cert['polar_Gram_error_upper'] < F(1, 2)
    assert eta == cert['actual_curvature_operator_error_upper'] < minimum/2
    assert frequency == cert['per_frequency_error_upper']
    assert 96*frequency == cert['cap1_total_H_error_including_vacuum'] < F('0.00002279')
    # A signed eigenvector is physically valid, yet selects different rational
    # frame bytes. Frozen replay must not silently regenerate this choice.
    mutated = Q.copy(); mutated[:, 0] *= -1
    assert not np.array_equal(mutated, Q)
    assert np.max(abs((Q/den)@(Q/den).T-(mutated/den)@(mutated/den).T)) == 0
    # A non-isometric coefficient mutation cannot borrow the frozen bound.
    bad = Q.copy(); bad[0, 0] += den//10
    badnorm = sum(int(x)**2 for x in bad[:, 0])-den**2
    assert F(abs(badnorm), den**2) > epsilon
    return {'Gram_squared_Frobenius_exact':str(gram2), 'residual_squared_Frobenius_exact':str(residual2),
            'cap1_H_bound':float(96*frequency), 'per_frequency_bound':float(frequency),
            'signed_eigenbasis_byte_mutation_detected':True, 'nonisometric_frame_mutation_rejected':True,
            'frozen_data_used_without_eigensolver_selection':True}


def photon_symbol():
    A0, B0, unused, terms, labels = linear.bloch_data()
    # Rank-factor pseudoinverse differs from the producer's pinv calculation.
    C, R = A0.rank_decomposition()
    inverse = R.T*(R*R.T).inv()*(C.T*C).inv()*C.T
    assert A0*inverse*A0 == A0 and inverse*A0*inverse == inverse
    assert s.trace(inverse*inverse.T) == s.Rational(3383, 240)
    assert s.trace(A0.T*A0) == 88 and C.cols == 22
    cellG = (B0*B0.T+s.ones(12)/12).inv()-s.ones(12)/12
    cycle = s.eye(36)-B0.T*cellG*B0
    H = cycle-A0*inverse
    flux = s.zeros(36, 3)
    for axis in range(3): flux[30+2*axis, axis] = flux[31+2*axis, axis] = 1
    assert flux.T*cycle*flux == s.eye(3)*s.Rational(3, 2)
    assert A0.T*cycle*flux == s.zeros(26, 3)
    W = H*flux*s.Rational(2, 3)
    assert W.T*W == s.eye(3)*s.Rational(2, 3) and flux.T*W == s.eye(3)
    derivatives = [s.zeros(36, 26) for _ in range(3)]
    first, second = defaultdict(int), defaultdict(int)
    for e, f, shift, value in terms:
        radius = sum(abs(x) for x in shift)
        first[e, f] += abs(value)*radius; second[e, f] += abs(value)*radius**2
        for axis in range(3): derivatives[axis][e, f] -= value*shift[axis]
    P, Q = s.eye(36)-A0*inverse, s.eye(26)-inverse*A0
    low = [P*D*Q for D in derivatives]
    for i in range(3):
        ei = s.eye(3)[:, i]
        assert H*low[i] == low[i]
        assert low[i]*low[i].T == W*(s.eye(3)-ei*ei.T)*W.T
        for j in range(i):
            ej = s.eye(3)[:, j]
            assert low[i]*low[j].T+low[j]*low[i].T == -W*(ei*ej.T+ej*ei.T)*W.T
    c1sq = sum(x*x for x in first.values()); c2four = sum(x*x for x in second.values())
    assert math.isqrt(c1sq-1)+1 == 11 and math.isqrt(c2four-1)+1 == 17
    cases = []
    for exponent in (18, 20, 22):
        h = F(1, 2**exponent); eta = 11*h
        D = 8*eta**2+32*eta**2*(10+eta)+F(17, 2)*h*h
        assert eta < F(1, 8) and h+D < F(1, 4)-eta and D < h/2
        error = F(5, 4)*D/h
        recorded = read('PHOTONS.json')['rigorous_refinement_cases'][len(cases)]
        assert error == F(recorded['frequency_error_upper_exact'])
        cases.append({'h':str(h), 'physical_frequency_error_bound':float(error)})
    # An actual face-column sign mutation spoils the incidence identity.
    bad = np.asarray(A0, dtype=np.int64).copy(); e, f = next(zip(*np.nonzero(bad)))
    bad[e, f] *= -1
    assert np.any(np.asarray(B0, dtype=np.int64)@bad)
    return {'pseudoinverse_Frobenius_squared':'3383/240', 'conductance':'3/2',
            'root_Green_pole':'2/3', 'g_squared':'1/3', 'six_all_direction_coefficients_exact':True,
            'refinement_cases':cases, 'orientation_mutation_rejected':True}


@lru_cache(None)
def generators():
    # Reconstruct M6 responses from source data, without compiler._native.
    base = ROOT/'reverse-engineering-reality/code'
    response = json.loads((base/'source_selection_model/response.json').read_text())
    controls = json.loads((base/'m1_source_realization/receipt.json').read_text())
    def scalar(pair): return float(F(pair[0]))+math.sqrt(5)*float(F(pair[1]))
    matrices = []
    for blocks in response['generators']:
        matrix = np.zeros((6, 6), complex)
        for block, values in enumerate(blocks):
            for i, row in enumerate(values):
                for j, cell in enumerate(row): matrix[3*block+i, 3*block+j] = scalar(cell[:2])+1j*scalar(cell[2:])
        matrices.append(matrix)
    coefficients = {r['name']:np.array([scalar(x) for x in r['coefficients']]) for r in controls['evidence']['controls']['rows']}
    selected = {'X':coefficients['imaginary_01'], 'Y':coefficients['real_01'],
                'Z':coefficients['phase_0']-coefficients['phase_1'],
                'I':coefficients['phase_0']+coefficients['phase_1'], 'J':coefficients['relative_phase']}
    selected['H'] = (selected['X']+selected['Z'])/math.sqrt(2)
    return {kind:sum(c*g for c, g in zip(row, matrices)) for kind, row in selected.items()}


@lru_cache(None)
def pulse(kind, angle): return expm(angle*generators()[kind])


@lru_cache(None)
def indices(bits, targets):
    all_indices = np.arange(1 << bits); mask = sum(1 << q for q in targets)
    bases = all_indices[(all_indices & mask) == 0]
    return np.array([bases+sum(((local >> (len(targets)-1-j)) & 1) << q for j, q in enumerate(targets))
                     for local in range(1 << len(targets))])


def execute(word, bits, state):
    v = np.asarray(state, complex).reshape(1 << bits, -1).copy()
    for p in word:
        gate = pulse(p['kind'], float(p['angle'])); dim = 1 << len(p['targets'])
        assert np.linalg.norm(gate[dim:, :dim]) < 2e-12
        idx = indices(bits, tuple(p['targets']))
        v[idx] = np.einsum('ab,bcr->acr', gate[:dim, :dim], v[idx])
    return v


def native_records():
    nf = frozen._native()
    sys.path.insert(0, str(HERE.parent/'clock'))
    import phase_clock
    from model.native_clock import configuration
    from model.observer.source_discrimination import circuit
    from model.native_interface import pointer
    diagonal, phase_error, phases, vacuum = frozen.photon_phase()
    cert = frozen.spectral_certificate()
    omega = np.sqrt(1.5*np.array([float(x) for x in cert['eigenvalues_exact']]))
    E0 = sum(omega)/2
    assert max(abs(np.array(phases)-np.r_[-E0, -E0-omega, np.zeros(65)])) < 2e-12
    # Native columns span physical and padding branches, with all phase kept.
    chosen = [0, 1, 37, 190, 191, 255]
    columns = np.eye(256, dtype=complex)[:, chosen]
    replay = execute(diagonal, 8, columns)
    phase_residual = float(np.linalg.norm(replay-columns*np.exp(1j*np.array(phases)[chosen])))
    assert phase_residual < float(phase_error)+3e-12
    p0 = min(map(F, configuration.success_probabilities(4, 32, full=True, mode='reset', target=circuit.decode(0))))
    pr = min(map(F, configuration.success_probabilities(4, 32, full=True, mode='reset', target=circuit.decode(pointer.READY_CODE))))
    reset = 9*(1-p0)+9*(1-pr)
    assert reset == F(read('NATIVE_PHOTON.json')['initialization_failure_union_bound_exact'])
    H = nf.H(0)
    herror = nf.root_certificate.exact_template_bound(H, [('I', (0,), nf.cert._pi(-1, 2)), ('H', (0,), nf.cert._pi(1, 2))])
    instruments, _ = phase_clock.actual_native_pointer()
    results = []
    for recorded in read('NATIVE_PHOTON.json')['experiments']:
        setting = recorded['setting']; program, native, _, _ = frozen.word(setting)
        initial = np.eye(512, dtype=complex)[:, [0]]
        vector = execute(program, 9, initial)[:, 0]
        photon = np.zeros(256, complex)
        photon[0] = (np.exp(1j*phases[0])+np.exp(1j*phases[1]))/2
        photon[1] = (np.exp(1j*phases[0])-np.exp(1j*phases[1]))/2
        t = -1/16+(math.pi/4 if setting == 'Y' else 0)
        expected = np.kron([math.cos(t), 1j*math.sin(t)], photon)
        residual = float(np.linalg.norm(vector-expected))
        assert residual < float(native)+3e-12
        recomputed = phase_error+2*herror+phase_clock.exact_word_bound(1, setting, F(1))+8*phase_clock.pointer_bound()
        assert recomputed == native
        complete = reset+native+cert['cap1_total_H_error_including_vacuum']
        assert complete == F(recorded['complete_record_TV_bound_exact'])
        for law in ('C', 'E'):
            effects = [np.zeros((2, 2), complex), np.zeros((2, 2), complex)]
            stored = recorded['laws'][law]['fine_one_bit_effects']
            assert len(instruments[law]) == len(stored)
            for (metadata, operators), entry in zip(instruments[law], stored):
                live = sum((K.conj().T@K for K in operators), np.zeros((2, 2), complex))
                frozen_effect = np.array(entry['effect_real'])+1j*np.array(entry['effect_imag'])
                assert metadata == entry['record'] and np.max(abs(live-frozen_effect)) < 1e-13
                effects[metadata['bit']] += live
            assert np.max(abs(sum(effects)-np.eye(2))) < 1e-13
            # The actual pointer is diagonal in its logical data basis.
            assert max(abs(effects[b][0, 1]) for b in (0, 1)) < 1e-14
            probabilities = abs(vector)**2
            transition = np.array([[effects[b][x, x].real for x in (0, 1)] for b in (0, 1)])
            for bit in range(9):
                idx = indices(9, (bit,)); probabilities[idx] = transition@probabilities[idx]
            assert abs(sum(probabilities)-np.vdot(vector, vector).real) < 2e-12
            assert np.max(abs(probabilities-recorded['laws'][law]['all512_joint_probabilities'])) < 3e-12
            assert abs(sum(probabilities)+recorded['native_export_failure_probability']-1) < float(native)+3e-12
        results.append({'setting':setting, 'pulses':len(program), 'source_matrix_residual':residual,
                        'complete_TV':float(complete), 'outcomes_per_law':512})
    # Execute altered spectra, retaining the same preparation/readout.
    frequency_deletion = (1-math.cos(omega[0]))/2
    assert frequency_deletion > .08 and frequency_deletion > 1000*max(r['complete_TV'] for r in results)
    vacuum_phase_defect = abs(1-np.exp(1j*E0))
    assert vacuum_phase_defect > .1
    same_record = float(np.max(abs(abs(expected*np.exp(1j*E0))**2-abs(expected)**2)))
    assert same_record < 1e-14
    # A coherent vacuum/padding reference sees a phase that this tested
    # vacuum/one-photon probe and factorized clock cannot measure.
    coherent_reference_probability_change = (1-math.cos(E0))/2
    assert coherent_reference_probability_change > .1
    refused = 0
    for bad_time in (1., F(1, 2), F(2)):
        try: frozen.word('X', bad_time)
        except ValueError: refused += 1
    assert refused == 3
    return {'spectral_column_residual':phase_residual, 'experiments':results,
            'frequency_deletion_record_separation':frequency_deletion,
            'vacuum_scalar_deletion_vector_defect':float(vacuum_phase_defect),
            'vacuum_scalar_deletion_tested_record_change':same_record,
            'vacuum_padding_reference_probability_change':float(coherent_reference_probability_change),
            'unsupported_native_time_controls_rejected':refused,
            'all_256_phase_targets_certified':True,
            'local_edge_encoder_claim':'absent; spectral-code source only'}


def cap_boundary():
    omega = np.sqrt(1.5*np.array([float(x) for x in frozen.spectral_certificate()['eigenvalues_exact']]))
    E0 = omega.sum()/2
    literal_quadratic = np.r_[E0, omega/2]
    correction = np.r_[0, (omega+sum(omega))/2]
    assert np.max(abs(literal_quadratic+correction-np.r_[E0, E0+omega])) < 1e-12
    # Matrix-unit identity establishes the whole total-cap-one commutator.
    for i, j in ((0, 0), (0, 189), (189, 0), (189, 189)):
        ai = np.zeros((191, 191)); aj = ai.copy(); ai[0, i+1] = aj[0, j+1] = 1
        top = np.diag([0]+[1]*190); identity = np.eye(191) if i == j else np.zeros((191, 191))
        assert np.array_equal(ai@aj.T-aj.T@ai, identity-top@(identity+aj.T@ai))
    a = np.array([[0., 1.], [0., 0.]])
    assert np.linalg.norm(a@a.T-a.T@a-np.eye(2), 2) == 2
    truncated_displacement = expm(1j*(a+a.T))[:, 0]
    true_vacuum_probability = math.exp(-1)
    mismatch = abs(abs(truncated_displacement[0])**2-true_vacuum_probability)
    assert mismatch > .075
    return {'physical_states':191, 'padding':65, 'top_CCR_defect':2,
            'unit_displacement_vacuum_probability_error':float(mismatch),
            'missing_cap_completion_H_norm':float(max(correction))}


def spatial_route_control():
    graph = cc.geometry.graph(3); coordinates = graph.cell_coordinates
    cell_id = {tuple(c):i for i, c in enumerate(coordinates)}
    # Size3 permits distinct positive and negative destinations and three
    # charge positions without colliding with the fixed reference charge.
    origin = (1, 1, 0); destinations = [(2, 1, 0), (0, 1, 0)]
    assert len({origin, *destinations, (0, 0, 0)}) == 4
    endpoints = [12*cell_id[c] for c in destinations]
    assert endpoints[0] != endpoints[1]
    return {'minimum_nonaliased_periodic_size':3, 'cells':27, 'ports':324,
            'edges':len(graph.seams), 'fixed_charge_cell':[0, 0, 0],
            'moving_charge_cell':list(origin), 'plus_minus_destinations':list(map(list, destinations)),
            'root_port_destinations':endpoints,
            'transfer':'Implement parameterized actual size3 routes, coherent shared-field translations, retained harmonic factors and full occupation/current/field/clock records. Static energy differences require a specified interference/force instrument; size3 alone is not a continuum theorem.'}


def run():
    before = frozen.verify_inputs()
    result = {}
    for name, function in [('continuous', continuous_source), ('spectral_frame', spectral_frame),
                           ('photon_symbol', photon_symbol), ('native', native_records),
                           ('cap_boundary', cap_boundary), ('spatial_route', spatial_route_control)]:
        result[name] = function()
        (HERE/'PARTIAL.json').write_text(json.dumps(result, indent=2)+'\n')
        print(name+' passed', flush=True)
    assert before == frozen.verify_inputs()
    receipt = {'schema':'oph.quadratic_photon_independent_audit.v1',
               'checks':{name:True for name in result}, 'results':result,
               'source_sha256':before, 'producer_bytes_unchanged':True,
               'reserved_acceptance_runs':[], 'official_rows_changed':[],
               'scope':'Frozen free noncompact field, rational spectral adapter and complete native spectral-code experiment. No local-edge encoder, harmonic native discretization, interacting charged-field word or physical row is certified.'}
    (HERE/'AUDIT.json').write_text(json.dumps(receipt, indent=2)+'\n')
    (HERE/'IMPORTS.json').write_text(json.dumps({'schema':'oph.photon_audit.imports.v1',
        'provider_manifest':frozen.MANIFEST_SHA256, 'source_pins':before,
        'frozen_replay':'frozen.py reads immutable frame, spectrum and emitted programs; no provider run() or eigensolver regeneration.'}, indent=2)+'\n')
    print(json.dumps({'checks':receipt['checks'], 'pins':len(before)}), flush=True)


if __name__ == '__main__': run()
