"""Complete finite model detector and exact original-H species attachment.

The 1024 physical CAR modes are not replaced by ten positional qubits.
The Gaussian one-particle execution represents the ideal comparison to a
separately certified physical-mode native network, not 2**1024 native rows.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from ..interacting_field.lattice import SIGMA
from ..contact_scattering.onshell import pi_interval
from ..relativistic_sector_audit.sector import sine_interval

SIDE = 8
VOLUME = SIDE**3
MODES = 2*VOLUME
SPEED = F(1, 8)


def address(site, spin):
    if (type(site) not in (tuple, list) or len(site) != 3 or
            any(type(x) is not int or not 0 <= x < SIDE for x in site) or
            type(spin) is not int or spin not in (0, 1)):
        raise ValueError('three side-eight coordinates and spin0 or1 required')
    return 2*((site[0]*SIDE+site[1])*SIDE+site[2])+spin


def decode(index):
    if type(index) is not int or not 0 <= index < MODES:
        raise ValueError('detector address in0..1023 required')
    site, spin = divmod(index, 2)
    xy, z = divmod(site, SIDE)
    x, y = divmod(xy, SIDE)
    return (x, y, z), spin


@lru_cache(maxsize=4)
def integer_operators(side=SIDE):
    """16H and16J_z at c=1/8, assembled from the original forward-bond rule."""
    if type(side) is not int or side not in (4, 8):
        raise ValueError('side4 source control or side8 detector required')
    sites = tuple(product(range(side), repeat=3))
    lookup = {x: j for j, x in enumerate(sites)}
    rows, cols, real, imag = [], [], [], []
    jr, jc, jv = [], [], []
    for x, site in enumerate(sites):
        for axis in range(3):
            y = list(site); y[axis] = (y[axis]+1) % side
            target = lookup[tuple(y)]
            for a, b in product(range(2), repeat=2):
                value = 1j*SIGMA[axis][a, b]
                if not value:
                    continue
                p, q = 2*x+a, 2*target+b
                rows.extend((p, q)); cols.extend((q, p))
                real.extend((int(value.real), int(value.real)))
                imag.extend((int(value.imag), -int(value.imag)))
                if axis == 2:
                    current = -int(SIGMA[2][a, b].real)
                    jr.extend((p, q)); jc.extend((q, p)); jv.extend((current, current))
    shape = (2*side**3,)*2
    R = sparse.coo_matrix((np.array(real, dtype=np.int64), (rows, cols)), shape=shape).tocsr()
    I = sparse.coo_matrix((np.array(imag, dtype=np.int64), (rows, cols)), shape=shape).tocsr()
    J = sparse.coo_matrix((np.array(jv, dtype=np.int64), (jr, jc)), shape=shape).tocsr()
    return R, I, J


def original_hamiltonian():
    real, imag, _ = integer_operators()
    return (real.astype(complex)+1j*imag)/16


def mode_catalog():
    """Sixteen measured modes: eight operational Weyl flavors, two signed branches."""
    result = []
    for bits in product((0, 1), repeat=3):
        signs = tuple((-1)**b for b in bits)
        chirality = -signs[0]*signs[1]*signs[2]
        labels = (4*bits[0], 4*bits[1], 4*bits[2]+1)
        for spin in (0, 1):
            sigma_z = 1-2*spin
            energy_sign = -signs[2]*sigma_z
            result.append({'flavor': ''.join(map(str, bits)),
                'valley_pi_coordinates': list(bits), 'momentum_grid_labels': list(labels),
                'momentum_pi_coefficients': [str(F(x, 4)) for x in labels],
                'detector_address': address(labels, spin), 'spin_occupation_label': spin,
                'model_number': 1, 'CAR_parity': -1,
                'jacobian_chirality': chirality,
                'energy_branch': '+' if energy_sign == 1 else '-',
                'energy_sqrt2_coefficient': str(F(energy_sign, 16)),
                'energy_display': energy_sign*np.sqrt(2)/16,
                'helicity': str(F(chirality*energy_sign, 2)),
                'helicity_convention': 'proper valley-frame cone spin; not bare CAR component or a physical spin assignment',
                'z_current_sqrt2_coefficient': str(F(energy_sign, 16)),
                'each_forward_z_bond_current_sqrt2_coefficient': str(F(energy_sign, 8192)),
                'continuum_mass': '0', 'electric_charge_assigned': False,
                'physical_species_name_assigned': False})
    return result


def algebraic_modes():
    """Integer A,B with normalized physical-mode columns psi=(A+sqrt2 B)/32."""
    catalog = mode_catalog()
    arrays = [np.zeros((MODES, len(catalog)), dtype=np.int64) for _ in range(4)]
    ar, ai, br, bi = arrays
    odd = {1: (1, 1), 3: (-1, 1), 5: (-1, -1), 7: (1, -1)}
    even = {0: (1, 0), 2: (0, 1), 4: (-1, 0), 6: (0, -1)}
    for column, mode in enumerate(catalog):
        labels, spin = mode['momentum_grid_labels'], mode['spin_occupation_label']
        for site in product(range(SIDE), repeat=3):
            phase = sum(x*k for x, k in zip(site, labels, strict=True)) % 8
            row = address(site, spin)
            if phase % 2:
                ar[row, column], ai[row, column] = odd[phase]
            else:
                br[row, column], bi[row, column] = even[phase]
    return tuple(arrays)


def exact_species_certificate():
    """Exact orthonormality, original-H spectral labels and nonzero current."""
    ar, ai, br, bi = algebraic_modes()
    real, imag, current = integer_operators()
    signs = np.array([int(F(row['energy_sqrt2_coefficient'])*16) for row in mode_catalog()], dtype=np.int64)
    def multiply(r, i): return real@r-imag@i, real@i+imag@r
    har, hai = multiply(ar, ai)
    hbr, hbi = multiply(br, bi)
    equations = (np.array_equal(har, 2*br*signs), np.array_equal(hai, 2*bi*signs),
                 np.array_equal(hbr, ar*signs), np.array_equal(hbi, ai*signs),
                 np.array_equal(current@ar, 2*br*signs), np.array_equal(current@ai, 2*bi*signs),
                 np.array_equal(current@br, ar*signs), np.array_equal(current@bi, ai*signs))
    if not all(equations):
        raise ArithmeticError('original-H energy or transported-current label is false')
    g0r = ar.T@ar+ai.T@ai+2*(br.T@br+bi.T@bi)
    g0i = ar.T@ai-ai.T@ar+2*(br.T@bi-bi.T@br)
    g1r = ar.T@br+ai.T@bi+br.T@ar+bi.T@ai
    g1i = ar.T@bi-ai.T@br+br.T@ai-bi.T@ar
    if not np.array_equal(g0r, 1024*np.eye(16, dtype=np.int64)) or any(np.any(x) for x in (g0i, g1r, g1i)):
        raise ArithmeticError('operational particle modes are not exactly orthonormal')
    commutator_real = real@current-current@real
    commutator_imag = imag@current-current@imag
    if not commutator_real.nnz and not commutator_imag.nnz:
        raise ArithmeticError('negative control: z-current unexpectedly conserved on every momentum')
    return {'number_of_modes': MODES, 'identified_orthonormal_modes': 16,
            'distinct_operational_Weyl_flavors': 8, 'signed_branches_per_flavor': 2,
            'source_integer_energy_and_current_equations': 8*MODES*16,
            'spectral_relation': 'H Pi_label = E_label Pi_label; Pi_label rank1',
            'outside_relation': '[H,I-sum_label Pi_label]=0',
            'unresolved_single_particle_outcomes': MODES-16,
            'contact_action_on_every_N1_input': 'zero for all source U>=0',
            'all_forward_z_bonds_have_nonzero_signed_flux': True,
            'z_current_stationary_on_selected_modes': True,
            'z_current_commutes_with_H_on_all_N1_states': False,
            'mode_normalization': 'psi=(A+sqrt(2)B)/32, Gaussian integer A,B',
            'CAR_exchange_statistics': 'fermionic by the supplied CAR algebra; not inferred from counting outcomes',
            'helicity_at_other_momenta_automatically_resolved': False}


def exact_mode_vectors():
    ar, ai, br, bi = algebraic_modes()
    return (ar+1j*ai+np.sqrt(2)*(br+1j*bi))/32


def fft_reference(state, inverse=False):
    """Independent numerical discrete Fourier transform, no source primitives."""
    array = np.asarray(state, dtype=complex)
    if array.ndim not in (1, 2) or array.shape[0] != MODES or not np.isfinite(array).all():
        raise ValueError('finite state columns on all1024 one-particle modes required')
    shape = (SIDE, SIDE, SIDE, 2)+array.shape[1:]
    operation = np.fft.ifftn if inverse else np.fft.fftn
    return operation(array.reshape(shape), axes=(0, 1, 2), norm='ortho').reshape(array.shape)


def refinement_table():
    """Exact N1 finite periodic family; no infinite-volume packet claim is reused."""
    pi_lo, pi_hi = pi_interval()
    result = []
    for n in range(5):
        a = F(1, 2**n)
        lower = SPEED/a*sine_interval(a*pi_lo/4)[0]
        upper = SPEED/a*sine_interval(a*pi_hi/4)[1]
        result.append({'side': 8*2**n, 'spacing': str(a), 'speed_coefficient': str(SPEED/a),
                       'physical_period_in_model_length': '8',
                       'fixed_physical_momentum': '(0,0,pi/4) about each valley',
                       'positive_energy_interval': [str(lower), str(upper)],
                       'positive_energy_display': float((lower+upper)/2),
                       'limiting_positive_energy_interval': [str(SPEED*pi_lo/4), str(SPEED*pi_hi/4)],
                       'limiting_energy': 'pi/32', 'continuum_mass': '0'})
    return result


def simulate_detector(duration=F(8)):
    """Native-network ideal Gaussian comparison; original H, every N1 output.

    Full literal-native predictions are certified separately with complete
    output-vector errors. This numerical calculation does not discard a
    literal-native off-number branch and claim to have executed it exactly.
    """
    from . import native_fft
    if type(duration) not in (int, F):
        raise ValueError('exact model duration required')
    catalog = mode_catalog()
    localized = np.zeros((MODES, len(catalog)), complex)
    for j, row in enumerate(catalog): localized[row['detector_address'], j] = 1
    prepared = native_fft.apply(localized, inverse=True)
    target = exact_mode_vectors()
    preparation_error = float(np.linalg.norm(prepared-target, 2))
    evolved = expm_multiply(-1j*float(duration)*original_hamiltonian(), prepared)
    measured = native_fft.apply(evolved)
    expected = np.array(localized)
    for j, row in enumerate(catalog): expected[:, j] *= np.exp(-1j*float(duration)*row['energy_display'])
    error = float(np.linalg.norm(measured-expected, 2))
    probabilities = np.abs(measured)**2
    gram_error = float(np.linalg.norm(measured.conj().T@measured-np.eye(len(catalog)), 2))
    if max(preparation_error, error, gram_error) > 2e-11:
        raise ArithmeticError('native-network comparison fails exact source mode or complete normalization')
    selected = [row['detector_address'] for row in catalog]
    return {'duration': str(F(duration)), 'ideal_preparation_difference_diagnostic': preparation_error,
            'complete_ideal_output_vector_difference_diagnostic': error,
            'complete_ideal_output_Gram_error_diagnostic': gram_error,
            'output_probabilities_all_1024_rows': probabilities.tolist(),
            'identified_outcome_confusion_matrix': probabilities[selected, :].tolist(),
            'comparison_kind': 'ideal Gaussian physical-mode network, not a literal2^1024-row native execution',
            'probabilities_normalized_or_postselected': False}


def sample_counts(seed=20261101, trials=65536):
    """Seeded ideal-model count illustration; literal-native errors stay separate."""
    if type(seed) is not int or seed < 0 or type(trials) is not int or trials != 65536:
        raise ValueError('nonnegative integer seed and the declared65536-trial campaign required')
    counts = np.random.Generator(np.random.PCG64(seed)).multinomial(trials, [1/16]*16)
    all_counts = [0]*MODES
    rows = []
    for mode, count in zip(mode_catalog(), counts, strict=True):
        all_counts[mode['detector_address']] = int(count)
        rows.append({'flavor': mode['flavor'], 'energy_branch': mode['energy_branch'],
                     'detector_address': mode['detector_address'], 'count': int(count)})
    if sum(all_counts) != trials:
        raise ArithmeticError('count illustration dropped an outcome')
    # Hoeffding plus a union bound: 2*16*exp(-2*n/64^2)=32*exp(-32).
    # The exact rational majorant uses e>8/3.
    statistical = F(1, 64)
    if any(abs(F(int(count), trials)-F(1, 16)) > statistical for count in counts):
        raise ArithmeticError('seeded illustration lies outside the stated simultaneous event')
    return {'sampler': 'PCG64', 'seed': seed, 'total_trials': trials,
            'input_distribution': 'uniform classical mixture of all16 source-prepared labels',
            'ideal_probabilities_exact': ['1/16']*16,
            'counts_by_identified_mode': rows, 'counts_all_1024_single_particle_addresses': all_counts,
            'vacuum_count': 0, 'unresolved_single_particle_count': 0,
            'multiple_occupation_or_source_abort_count': 0,
            'sampling_error_bound': str(statistical),
            'simultaneous_failure_probability_upper': str(32*F(3, 8)**32),
            'counts_are_ideal_model_simulation': True,
            'counts_are_literal_native_full_Fock_execution': False,
            'native_error_added_separately': True,
            'conditioned_on_one_particle_or_success': False}
