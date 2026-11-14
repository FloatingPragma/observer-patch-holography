"""Complete-observable continuum bounds for the declared LOCAL12 walk.

The finite-dimensional Bloch generator is computed from the same reflection
and routing. No continuum generator is an input to source execution.
"""
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np
import sympy as s
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'reverse-engineering-reality'))
sys.path.insert(0, str(ROOT / 'openwave-tuning/auditor/native_walk_bridge'))
from model.native_spatial_refinement import source as spatial
import common_walk as common

CLOCK_PATH = HERE.parent / 'clock/clock_contract.py'
spec = importlib.util.spec_from_file_location('local12_common_clock', CLOCK_PATH)
clock = importlib.util.module_from_spec(spec)
spec.loader.exec_module(clock)
MODEL = 'native-local12-current-v0-proposed'
SOURCE_ERROR_PATH = ROOT / 'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@lru_cache(None)
def rational_bloch_data():
    """Exact integer numerators; every multiplication stays below int64 bounds."""
    size = 144
    I = np.eye(size, dtype=np.int64)
    Rn = np.kron(np.eye(12, dtype=np.int64),
                 np.ones((12, 12), dtype=np.int64) - 6 * np.eye(12, dtype=np.int64))
    S0 = I.copy()
    for direction, matching in enumerate(common.INTERNAL_MATCHINGS):
        event = 2 * direction + 1
        for a, b in matching:
            x, y = 12 * a + event, 12 * b + event
            S0[x, x] = S0[y, y] = 0
            S0[x, y] = S0[y, x] = 1
    D = [np.zeros((size, size), dtype=np.int64) for _ in range(3)]
    for axis, row in enumerate(spatial.links()):
        for a, b in row:
            x, y = 12 * a + 11, 12 * b + 11
            S0[x, x] = S0[y, y] = 0
            S0[x, y] = S0[y, x] = 1
            D[axis][x, y], D[axis][y, x] = 1, -1
    assert np.array_equal(S0 @ S0, I) and np.array_equal(Rn @ Rn, 36 * I)
    Wn = S0 @ Rn
    Cn = Wn + Wn.T  # C=(W+W*)/2 has denominator12.
    assert spatial.exact_symbol_zero().eigenvals() == {s.Integer(0): 1,
        s.Integer(4): 3, s.Integer(6): 2, s.Integer(8): 6}
    Pn = I.copy()
    for shift in (-12, 8, 6, 4):
        factor = Cn - shift * I
        assert size * int(abs(Pn).max()) * int(abs(factor).max()) < 2**63
        Pn = Pn @ factor
    Pd = 4608
    divisor = math.gcd(Pd, int(np.gcd.reduce(Pn.reshape(-1))))
    Pn //= divisor
    Pd //= divisor
    assert size * int(abs(Pn).max())**2 < 2**63
    assert np.array_equal(Pn, Pn.T) and np.array_equal(Pn @ Pn, Pd * Pn)
    assert np.array_equal(Wn @ Pn, 6 * Pn) and int(np.trace(Pn)) == 26 * Pd
    one = np.ones(size, dtype=np.int64)
    assert np.array_equal(Pn @ one, Pd * one)
    H, nums, dens, vectors, displacement = [], [], [], [], []
    for d in D:
        current = Pn.copy()
        for factor in (d, Rn, Pn):
            assert size * int(abs(current).max()) * int(abs(factor).max()) < 2**63
            current = current @ factor
        denominator = 6 * Pd * Pd
        divisor = math.gcd(denominator, int(np.gcd.reduce(current.reshape(-1))))
        current //= divisor
        denominator //= divisor
        assert np.array_equal(current, current.T)
        oriented = d @ S0
        assert np.array_equal(oriented, np.diag(np.diag(oriented)))
        assert sorted(np.diag(oriented)) == [-1] * 2 + [0] * 140 + [1] * 2
        assert np.array_equal((Pn @ oriented @ Pn) * denominator, current * (Pd * Pd))
        # After S, arriving at a negative endpoint records a positive-cell
        # displacement. The actual terminal observable is -oriented.
        displacement.append(-oriented)
        gnum = current @ one
        assert int(gnum.sum()) == 0
        assert np.array_equal(144 * current, gnum[:, None] + gnum[None, :])
        H.append(current.astype(float) / denominator)
        nums.append(gnum)
        dens.append(denominator)
        vectors.append(gnum.astype(float) / (12 * denominator))
    for a in range(3):
        for b in range(3):
            dot = sum(int(x) * int(y) for x, y in zip(nums[a], nums[b]))
            assert dot == (3 * dens[a] * dens[b] if a == b else 0)
    return {'R': Rn.astype(float) / 6, 'S0': S0, 'D': D,
            'W0': Wn.astype(float) / 6, 'P': Pn.astype(float) / Pd,
            'H': H, 'scalar': np.ones(size) / 12, 'vectors': vectors,
            'native_terminal_displacement': displacement,
            'exact_projector_denominator': Pd, 'exact_gradient_matrix_denominators': dens}


def bloch_walk(k):
    data = rational_bloch_data()
    S = np.array(data['S0'], complex)
    for axis, row in enumerate(spatial.links()):
        for a, b in row:
            x, y = 12 * a + 11, 12 * b + 11
            S[x, y], S[y, x] = np.exp(1j * k[axis]), np.exp(-1j * k[axis])
    return S @ np.array(data['R'], float)


def effective_generator(p, tau0=F(1)):
    tau0 = F(tau0)
    if tau0 <= 0 or len(p) != 3:
        raise ValueError('three momenta and one positive shared clock scale required')
    data = rational_bloch_data()
    return 1j * sum(float(p[a]) * np.array(data['H'][a], float) for a in range(3)) / float(tau0)


def observer_error_bound(steps, cell_spacing, momentum_bound, tau0=F(1)):
    """Trace-distance/TV upper bound, including arbitrary spectator reference.

    Applies to initial states in the full26-dimensional W(0)=1 sector and
    Fourier support |p|_infinity<=K, followed by any one terminal instrument.
    Full-state dynamics and all terminal outcomes are retained.
    """
    if type(steps) is not int or steps < 0:
        raise ValueError('nonnegative integer coherent steps required')
    h, K, tau = F(cell_spacing), F(momentum_bound), F(tau0)
    if h <= 0 or K < 0 or tau <= 0 or h * K > F(1, 2):
        raise ValueError('positive h,tau0 and 0<=hK<=1/2 required')
    time = clock.effective_time(steps, h, tau)
    # ||X||<=sqrt(3/2)K<=5K/4. The exact effective generator has
    # ||K(p)||=|p|/sqrt48<=K/4. The one-step defect is bounded by
    # h²K²[5/4+1/2+1/32+5/16+(5/8)hK] <=(77/32)h²K².
    preparation = F(5, 2) * h * K
    evolution = F(77, 32) * steps * h * h * K * K
    return {'time': time, 'embedding_error': preparation, 'evolution_error': evolution,
            'state_vector_error': preparation + evolution,
            'trace_distance_and_terminal_TV_bound': min(F(1), preparation + evolution),
            'continuum_cell_bin_TV_bound': min(F(1), preparation + evolution + F(3, 4) * h * K),
            'native_displacement_mean_error_bound': 2 * min(F(1), preparation + evolution),
            'native_velocity_mean_error_bound': 2 * min(F(1), preparation + evolution) / tau,
            'clock_scale': tau, 'cell_spacing': h, 'momentum_bound': K}


def scalar_dispersion_bound(steps, cell_spacing, momentum_bound, tau0=F(1)):
    h, K, tau = F(cell_spacing), F(momentum_bound), F(tau0)
    observer_error_bound(steps, h, K, tau)
    r = h * K
    gap = 4 - r - r * r / 2
    eta = r * r / (gap * gap)
    coefficient = F(1, 16) + 1 / (576 * (1 - r * r / 48))
    # sqrt48 < 7 gives a rational uniform phase bound, including p=0.
    phase = 7 * steps * h**3 * K**3 * coefficient
    return {'optical_weight_bound': eta, 'phase_error_bound': phase,
            'normalized_low_band_trace_distance_bound': min(F(1), r / gap + phase),
            'normalized_low_band_state_vector_bound': min(F(2), 2 * r / gap + phase),
            'X_probe_TV_bound': min(F(1), eta + phase / 2),
            'Y_probe_probability_exact': F(1, 2),
            'speed_squared_per_common_time': 1 / (48 * tau * tau),
            'domain': 'uniform internal/event preparation, bandlimited envelope, no intervening source measurement'}


def compiled_size2_budget(steps, ideal_tv=F(0), instrument_error=F(0), reset_timeout=F(0), member_size=2):
    """Compose the frozen size2 source certificate, without family transfer."""
    if type(steps) is not int or steps < 0 or member_size != 2:
        raise ValueError('the available source-word certificate is only the size2 member')
    values = tuple(map(F, (ideal_tv, instrument_error, reset_timeout)))
    if any(v < 0 for v in values):
        raise ValueError('nonnegative complete-instrument error bounds required')
    certificate = json.loads(SOURCE_ERROR_PATH.read_text())
    assert certificate['candidate_id'] == MODEL and all(certificate['checks'].values())
    source = steps * F(certificate['whole_W_operator_error_exact'])
    preparation = F(certificate['whole_source_preparation_error_exact'])
    return {'member_size': 2, 'source_word_error_exact': str(source),
            'source_preparation_error_exact': str(preparation),
            'complete_record_TV_bound_exact': str(min(F(1), sum(values) + source + preparation)),
            'ideal_bound_input_exact': str(values[0]), 'instrument_error_input_exact': str(values[1]),
            'reset_timeout_input_exact': str(values[2]),
            'scope': 'Same size2 native source; no transfer to compiled words of other sizes.'}


def matrix_controls():
    data = rational_bloch_data()
    P = np.array(data['P'], float)
    W0 = np.array(data['W0'], float)
    Q = np.eye(144) - P
    # The inverse is on Q; P regularizes it without changing its Q action.
    inverse = np.linalg.inv(W0 - np.eye(144) + P) @ Q
    assert np.linalg.norm(inverse, 2) < 1.25
    records = []
    for p in ((1., 0., 0.), (1., 1., 0.), (1., -1., 1.)):
        K = max(abs(v) for v in p)
        W1 = 1j * sum(p[a] * np.array(data['D'][a] @ data['R'], float) for a in range(3))
        generator = effective_generator(p)
        X = -inverse @ Q @ W1 @ P
        assert np.linalg.norm(W0 @ X + W1 @ P - X - generator, 2) < 1e-12
        assert np.linalg.norm(generator + generator.conj().T, 2) < 1e-12
        for denominator in (8, 16, 32, 64):
            h = F(1, denominator)
            steps = denominator
            W = bloch_walk([float(h) * v for v in p])
            J = P + float(h) * X
            step_defect = np.linalg.norm(W @ J - J @ expm(float(h) * generator), 2)
            certified_defect = F(77, 32) * h * h * F(K)**2
            full_error = np.linalg.norm(np.linalg.matrix_power(W, steps) @ P
                                        - expm(steps * float(h) * generator) @ P, 2)
            bound = observer_error_bound(steps, h, F(K))
            assert step_defect <= float(certified_defect) + 1e-12
            assert full_error <= float(bound['state_vector_error']) + 1e-12
            records.append({'momentum': p, 'h_exact': str(h), 'steps': steps,
                'time_exact': str(bound['time']), 'full26sector_operator_error': float(full_error),
                'state_error_bound_exact': str(bound['state_vector_error']),
                'one_step_corrector_error': float(step_defect),
                'one_step_corrector_bound_exact': str(certified_defect),
                'all144_internal_event_coordinates_retained': True})
    return records


def development_geometry_controls():
    """Whole graph Fourier action, only existing development sizes2,3,4."""
    R = np.array(rational_bloch_data()['R'], float)
    rows = []
    for size in (2, 3, 4):
        graph = spatial.graph(size)
        coords = np.array(graph.cell_coordinates)
        for mode in ((0, 0, 0), (1, 0, 0), (1, 1, 0)):
            k = 2 * np.pi * np.array(mode) / size
            # Full graph action on one arbitrary144-component Bloch column;
            # independent direction assignment from actual graph neighbours.
            v = np.arange(1, 145, dtype=float) + 1j * np.arange(144, 0, -1)
            v /= np.linalg.norm(v)
            phased = np.exp(1j * coords @ k)
            source = (phased[:, None] * v[None, :]).reshape(graph.ports, 12)
            value = source @ R[:12, :12].T
            out = value.copy()
            internal_colours = {tuple(edge): d for d, matching in enumerate(common.INTERNAL_MATCHINGS)
                                for edge in matching}
            for a, b in graph.seams:
                direction = internal_colours[(a % 12, b % 12)] if a // 12 == b // 12 else 5
                event = 2 * direction + 1
                out[a, event], out[b, event] = value[b, event], value[a, event]
            target = phased[:, None] * (bloch_walk(k) @ v)[None, :]
            error = float(np.max(abs(out.reshape(graph.cells, 144) - target)))
            assert error < 1e-12
            rows.append({'development_size': size, 'mode': mode, 'full_graph_symbol_error': error})
    return rows


def run():
    paths = [Path(__file__), Path(spatial.__file__), Path(common.__file__), CLOCK_PATH,
             ROOT/'reverse-engineering-reality/model/native_spatial_refinement/PROOF.md',
             ROOT/'openwave-tuning/auditor/native_walk_bridge/COMMON_WALK_RECEIPT.json',
             ROOT/'openwave-tuning/auditor/common_walk_particle/MANIFEST.json',
             SOURCE_ERROR_PATH]
    pins = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    data = rational_bloch_data()
    print('exact26sector and rank4 scalar/vector generator checked', flush=True)
    matrices = matrix_controls()
    print('complete Bloch operator bounds checked', flush=True)
    graphs = development_geometry_controls()
    # Bounds on a unit torus use |2pi m|∞<7 for |m|∞<=1. These are
    # symbolic inequalities, not enumerations or acceptance runs at these L.
    refinement = []
    for size in (128, 256, 512, 1024, 2048):
        b = observer_error_bound(size, F(1, size), F(7))
        d = scalar_dispersion_bound(size, F(1, size), F(7))
        refinement.append({'formal_refinement_size': size, 'enumerated': False,
                           'whole_observable_bound': {k: str(v) for k, v in b.items()},
                           'scalar_probe_bound': {k: str(v) for k, v in d.items()}})
    mutations = []
    for kwargs in ({'steps': -1, 'cell_spacing': F(1, 16), 'momentum_bound': 1},
                   {'steps': 1, 'cell_spacing': 1, 'momentum_bound': 1},
                   {'steps': 1, 'cell_spacing': F(1, 16), 'momentum_bound': 1, 'tau0': 0}):
        try:
            observer_error_bound(**kwargs)
        except ValueError:
            mutations.append(True)
        else:
            mutations.append(False)
    assert all(mutations)
    try:
        compiled_size2_budget(8, member_size=4)
    except ValueError:
        different_member_rejected = True
    else:
        different_member_rejected = False
    assert different_member_rejected
    # A common nonzero phase sign is observable and cannot be removed.
    p = (1., 0., 0.)
    P = np.array(data['P'], float)
    h = F(1, 64)
    W = bloch_walk([float(h), 0., 0.])
    wrong_sign_error = float(np.linalg.norm(-W @ P - expm(float(h) * effective_generator(p)) @ P, 2))
    assert wrong_sign_error > 1.9
    assert pins == {str(p.relative_to(ROOT)): sha(p) for p in paths}
    result = {'schema': 'oph.local12.controlled_continuum.v1', 'model_id': MODEL,
        'checks': {'exact_zero_bloch_projector': True, 'effective_generator_source_derived': True,
                   'three_gradient_vectors_Gram_is_I_over48': True,
                   'native_signed_displacement_has_bilinear_amplitude_map': True,
                   'full26sector_operator_error_bounds': True, 'complete_graph_Bloch_binding_L234': True,
                   'both_acoustic_branches_and24_stationary_modes_retained': True,
                   'optical_preparation_weight_bounded': True, 'one_common_clock_imported': True,
                   'invalid_domains_and_global_sign_rejected': True,
                   'compiled_source_error_composed_without_family_substitution': True,
                   'source_pins_unchanged': True},
        'zero_phase_sector_dimension': 26, 'dynamic_union_dimension': 4,
        'exact_projector_denominator': data['exact_projector_denominator'],
        'exact_gradient_matrix_denominators': data['exact_gradient_matrix_denominators'],
        'generic_momentum_acoustic_eigenvalues': ['+i|p|/(sqrt48*tau0)', '-i|p|/(sqrt48*tau0)'],
        'generic_momentum_zero_modes': 24,
        'effective_field_equations': 'd_t phi=(1/(sqrt48*tau0)) div chi; d_t chi=(1/(sqrt48*tau0)) grad phi; chi is a vector amplitude and22 further amplitude components are stationary',
        'native_current_attachment': 'The actual terminal signed cell-displacement observable has compressed operator -H_a. For amplitudes(phi,chi), its first moment is -2/sqrt48 Re(conj(phi)*chi_a); dividing by tau0 gives continuum probability velocity flux. chi itself is not the native edge-current record or an electromagnetic field. The equality is a finite compressed observable identity, with complete-observable evolution error retained.',
        'complete_observable_theorem': 'For initial reference-entangled states in P0 and Fourier support |p|∞<=K, hK<=1/2, at t=n*tau0*h, state trace distance and any complete terminal-instrument TV are bounded by min(1,(5/2)hK+(77/32)n*h²*K²). No optical or flat output is discarded.',
        'continuum_position_bins': 'Comparing native cell-position probabilities with continuum density integrated over centered cells adds at most3hK/4 to the complete terminal bound; no equality between a discrete grid measure and an unbinned continuous measure is claimed.',
        'Bloch_controls': matrices, 'development_graph_binding': graphs,
        'formal_unit_torus_bounds': refinement,
        'frozen_size2_native_error_compositions': [compiled_size2_budget(n) for n in (1, 8, 64)],
        'implementation_refinement_condition': 'For fixed effective time, n is O(1/h). Uniform physical implementation convergence requires per-step error delta_W(h)=o(h), vanishing preparation/instrument errors and retained reset timeouts with total probability tending to zero. The fixed size2 certificate does not establish these family conditions.',
        'global_sign_mutation_operator_error': wrong_sign_error,
        'clock_contract': {'source': str(CLOCK_PATH.relative_to(ROOT)), 'api': 'effective_time(n,h,tau0)',
                           'physical_seconds_claimed': False},
        'instrument_domain': 'One complete terminal instrument, or a sequence whose intervening operations preserve Fourier bandwidth and P0 with certified error. Arbitrary sharp position interventions are excluded from that theorem; exact same-W full-state continuation is retained separately.',
        'native_execution_scope': 'The family Bloch/graph update is exact. Native pulse compilation is frozen for size2; no uniform-in-size physical pulse implementation or accumulated hardware error is inferred. Bandlimited input preparation is an explicit domain premise.',
        'field_scope': 'Free scalar/longitudinal first-order continuum on this one-load source. Stationary transverse components do not establish Maxwell radiation. No interacting manybody limit, species masses, coupling or all31 validation is inferred.',
        'source_sha256': pins, 'official_rows_changed': []}
    (HERE/'RESULT.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'checks': result['checks'], 'maximum_controlled_error': max(r['full26sector_operator_error'] for r in matrices)}))


if __name__ == '__main__':
    run()
