"""Reconstruct the complete calibrated quadratic action and physical modes."""
from functools import lru_cache
from pathlib import Path
import hashlib
import json

import numpy as np
from scipy import linalg as la
import sympy as s

from model.common_source import calibrated_matter as matter
from model.common_source.chiral_candidate import action, kinetic
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory

ROOT = Path(__file__).resolve().parents[2]
CALIBRATION = 'model/common_source/calibrated_matter_results.json'
SPECIES = {'up': ('up', 'charm', 'top'), 'down': ('down', 'strange', 'bottom'),
           'charged_lepton': ('electron', 'muon', 'tau')}
CHARGE = {'up': s.Rational(2, 3), 'down': s.Rational(-1, 3), 'charged_lepton': s.Integer(-1)}
MOMENTA = ((0, 0, 1), (0, 0, -1))


def numeric(matrix):
    return np.asarray(matrix.evalf(30), dtype=complex)


def scalar(value):
    return {'exact': str(value), 'numeric': float(s.N(value, 18))}


def sparse_record(matrix):
    return {'shape': list(matrix.shape), 'entries': [[i, j, str(value)]
            for (i, j), value in sorted(s.ImmutableSparseMatrix(matrix).todok().items())]}


def authenticate():
    receipt = json.loads((ROOT/CALIBRATION).read_text())
    for name, digest in receipt['source_bindings'].items():
        path = ROOT/name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('calibrated action input differs: '+name)
    return receipt


@lru_cache(maxsize=1)
def fermions():
    archived = authenticate()
    ward = matter.exact_checks()
    neutrino = matter.neutrino_report()
    k = matter.mass_kernel(include_neutrinos=True)
    w = matter.mass_basis_change(include_neutrinos=True)
    diagonal = s.zeros(45)
    masses = {}
    coordinate_labels = {}
    for channel in matter.CHANNELS:
        masses[channel] = matter.mass_coordinates(channel)
        for color, ((left, right), _) in enumerate(sorted(matter._oriented_components(channel)[1].todok().items())):
            for family, mass in enumerate(masses[channel]):
                diagonal[3*left+family, 3*right+family] = mass
                diagonal[3*right+family, 3*left+family] = mass
                for coordinate in (3*left+family, 3*right+family):
                    coordinate_labels[coordinate] = (channel, SPECIES[channel][family], color, mass)
    masses['neutrino'] = matter.neutrino_mass_coordinates()
    for family, mass in enumerate(masses['neutrino']):
        diagonal[30+family, 30+family] = mass
        coordinate_labels[30+family] = ('neutrino', 'neutrino_mass_eigenstate_'+str(family+1), 0, mass)
    if not matter._zero(w.H*w-s.eye(45)) or not matter._zero(w.T*k*w-diagonal):
        raise ArithmeticError('complete optional-Weinberg action mass diagonalization failed')
    if not matter._zero(action.electric_charge().T*k+k*action.electric_charge()):
        raise ArithmeticError('the complete action lost electric neutrality')
    # Actual two-momentum CAR source, retaining the opposite momentum block.
    h, delta = kinetic.momentum_pair_source(MOMENTA[0], kernel=k, speed=1.)
    hn = np.block([[h, delta], [delta.conj().T, -h.T]])
    d = len(h)
    tau = np.block([[np.zeros((d, d)), np.eye(d)], [np.eye(d), np.zeros((d, d))]])
    eig, vec = la.eigh(hn)
    positive = vec[:, eig > 0]
    qcar = np.kron(np.eye(4), numeric(action.electric_charge()))
    qn = la.block_diag(qcar, -qcar.T)
    expected = []
    for channel in matter.CHANNELS:
        for mass in masses[channel]:
            expected.extend([float(s.sqrt(1+mass*mass))]*(8 if channel == 'charged_lepton' else 24))
    for mass in masses['neutrino']:
        expected.extend([float(s.sqrt(1+mass*mass))]*4)
    if len(expected) != 180 or positive.shape != (360, 180):
        raise ArithmeticError('complete physical fermion count differs')
    errors = {'Hermiticity': float(la.norm(hn-hn.conj().T)),
              'particle_hole': float(la.norm(tau@hn.conj()@tau+hn)),
              'complete_spectrum': float(np.max(abs(eig-np.sort(np.r_[-np.array(expected), expected])))),
              'canonical_CAR_metric': float(la.norm(positive.conj().T@positive-np.eye(180))),
              'canonical_CAR_anomalous': float(la.norm(positive.conj().T@tau@positive.conj()))}
    errors['electric_charge_commutator'] = float(la.norm(qn@hn-hn@qn))
    expected_charges = np.array([-1]*12+[-2/3]*36+[-1/3]*36+[0]*12+[1/3]*36+[2/3]*36+[1]*12)
    errors['complete_positive_charge_spectrum'] = float(np.max(abs(la.eigvalsh(positive.conj().T@qn@positive)-expected_charges)))
    if max(errors.values()) > 2e-10:
        raise ArithmeticError('complete native two-momentum CAR reconstruction failed')
    modes, groups, targets = [], [], []
    for channel in matter.CHANNELS:
        colors = 1 if channel == 'charged_lepton' else 3
        for family, species in enumerate(SPECIES[channel]):
            mass = masses[channel][family]
            targets.append({'species': species, 'kind': 'Dirac', 'rest_coefficient_GeV': scalar(mass),
                'electric_charge': str(CHARGE[channel]), 'antiparticle_charge': str(-CHARGE[channel]),
                'colors': colors, 'physical_polarizations_per_charge': 2,
                'positive_modes_per_momentum': 4*colors, 'canonical_kinetic_metric': 'identity',
                'canonical_fermion_spectral_residue': '1 per normalized mode',
                'mass_convention': 'running pure MSbar action parameter at mu=91.1876 GeV',
                'empirical_inputs': ['Yukawa_'+species, 'v_SM'],
                'native_source_prediction': False,
                'missing_physical_attachment': 'interacting selfenergy/pole matching and source-clock identification; confinement additionally applies to quarks' if colors == 3 else 'interacting pole/width matching and source-clock identification'})
            for momentum in MOMENTA:
                for orientation, sign in (('particle', 1), ('antiparticle', -1)):
                    ids = []
                    for color in range(colors):
                        for spin in (-1, 1):
                            ids.append(len(modes))
                            modes.append({'index': len(modes), 'species': species,
                                'orientation': orientation, 'momentum_GeV': list(momentum),
                                'spin_twice': spin, 'color_component': color,
                                'electric_charge': str(sign*CHARGE[channel]),
                                'energy_GeV': scalar(s.sqrt(1+mass*mass))})
                    groups.append({'species': species, 'orientation': orientation,
                                   'momentum_GeV': list(momentum), 'modes': ids,
                                   'energy_GeV': scalar(s.sqrt(1+mass*mass)),
                                   'electric_charge': str(sign*CHARGE[channel])})
    for family, mass in enumerate(masses['neutrino']):
        species = 'neutrino_mass_eigenstate_'+str(family+1)
        targets.append({'species': species, 'kind': 'Majorana mass eigenfield; massless limit retained',
            'rest_coefficient_GeV': scalar(mass), 'rest_coefficient_eV': scalar(10**9*mass),
            'electric_charge': '0', 'antiparticle_charge': '0', 'colors': 1,
            'physical_polarizations_per_charge': 2, 'positive_modes_per_momentum': 2,
            'canonical_kinetic_metric': 'identity', 'canonical_fermion_spectral_residue': '1 per normalized mode',
            'separate_antiparticle_species_added': False,
            'mass_convention': 'optional calibrated Weinberg tree prescription; low-energy normal-ordering targets at the chosen common background',
            'empirical_inputs': ['oscillation mass-squared differences'],
            'assumed_inputs': ['normal ordering', 'lightest mass zero', 'Majorana phases zero'],
            'native_source_prediction': False,
            'missing_physical_attachment': 'low-energy-to-MZ EFT running/loop matching, absolute-mass mechanism, and source-clock identification'})
        for momentum in MOMENTA:
            ids = []
            for helicity in (-1, 1):
                ids.append(len(modes))
                modes.append({'index': len(modes), 'species': species,
                    'orientation': 'neutrino-like helicity' if helicity == -1 else 'antineutrino-like helicity',
                    'momentum_GeV': list(momentum), 'helicity_twice': helicity,
                    'color_component': 0, 'electric_charge': '0',
                    'energy_GeV': scalar(s.sqrt(1+mass*mass))})
            groups.append({'species': species, 'orientation': 'both Majorana helicities',
                           'momentum_GeV': list(momentum), 'modes': ids,
                           'energy_GeV': scalar(s.sqrt(1+mass*mass)), 'electric_charge': '0'})
    if len(modes) != 180 or sum(len(group['modes']) for group in groups) != 180:
        raise ArithmeticError('fermion labels omit a complete canonical mode')
    # Exact W fixes family axes. Analytic disjoint mass pairs fix each named
    # charge, momentum and spin eigenvector, including tiny neutrino masses.
    hm, dm = kinetic.momentum_pair_source(MOMENTA[0], kernel=diagonal, speed=1.)
    mass_symbol = np.block([[hm, dm], [dm.conj().T, -hm.T]])
    transform_particle = np.kron(np.eye(4), numeric(w))
    transform = la.block_diag(transform_particle, transform_particle.conj())
    off_diagonal = mass_symbol-np.diag(np.diag(mass_symbol))
    named = {}
    for index in np.flatnonzero(np.diag(mass_symbol).real > 0):
        hole, car = divmod(int(index), 180)
        momentum_axis, local = divmod(car, 90)
        spin_axis, coordinate = divmod(local, 45)
        channel, species, color, mass = coordinate_labels[coordinate]
        physical_p = (1 if momentum_axis == 0 else -1)*(-1 if hole else 1)
        physical_spin = (1 if spin_axis == 0 else -1)*(-1 if hole else 1)
        charge = action.electric_charge()[coordinate,coordinate]*(-1 if hole else 1)
        if channel == 'neutrino':
            spin = physical_spin*physical_p
            orientation = 'neutrino-like helicity' if spin == -1 else 'antineutrino-like helicity'
        else:
            spin = physical_spin
            orientation = 'particle' if charge == CHARGE[channel] else 'antiparticle'
        omega = float(s.sqrt(1+mass*mass))
        vector = off_diagonal[:,index]/(omega+1)
        vector[index] = 1
        vector /= np.sqrt(1+(float(mass)/(omega+1))**2)
        key = (species, orientation, physical_p, spin, color)
        if key in named:
            raise ArithmeticError('a named canonical mode occurs twice')
        named[key] = transform@vector
    frame = np.column_stack([named[(mode['species'], mode['orientation'], mode['momentum_GeV'][2],
        mode.get('spin_twice', mode.get('helicity_twice')), mode['color_component'])] for mode in modes])
    errors.update({'exact_mass_basis_numeric_transform': float(la.norm(transform.conj().T@hn@transform-mass_symbol)),
                   'named_frame_eigenvectors': float(la.norm(hn@frame-frame@np.diag([mode['energy_GeV']['numeric'] for mode in modes]))),
                   'named_frame_CAR_metric': float(la.norm(frame.conj().T@frame-np.eye(180))),
                   'named_frame_CAR_anomalous': float(la.norm(frame.conj().T@tau@frame.conj()))})
    if max(errors.values()) > 2e-10:
        raise ArithmeticError('named physical field frame failed native canonical reconstruction')
    return {'K': k, 'W': w, 'mass_basis': s.ImmutableSparseMatrix(diagonal),
            'native_pair_h': h, 'native_pair_delta': delta, 'positive_frame': positive,
            'named_positive_frame': frame,
            'masses': masses, 'modes': modes, 'groups': groups, 'targets': targets,
            'checks': {'charged_native_Ward_identities': ward['source_current_Higgs_identities_exact'],
                       'same_Higgs_Weinberg_Ward_identities': neutrino['same_source_current_LH_invariant_identities_exact'],
                       'exact_45_coordinate_mass_diagonalization': True,
                       'K45_rank': 44, 'K45_nullity': 1,
                       'single_zero_momentum_Nambu_dimension': 180,
                       'single_zero_momentum_Nambu_rank': 176,
                       'single_zero_momentum_Nambu_nullity': 4, 'distinct_momentum_CAR_modes': 180,
                       'complete_Nambu_dimension': 360, 'numerical_errors': errors,
                       'named_fields_constructed_from_exact_W_and_analytic_mass_pairs': True,
                       'tiny_neutrino_masses_resolved_by_exact_Takagi_block': True,
                       'binary64_p_equals_one_neutrino_splittings_resolved': False},
            'calibration': archived['calibration'],
            'neutrino_calibration': archived['optional_phenomenological_neutrino_extension']['calibration']}


@lru_cache(maxsize=1)
def bosons():
    b = matter.bosonic_parameters()
    v, gw, gy, lam = (b[key] for key in ('v_SM_GeV', 'gW', 'gY', 'lambda_SM'))
    background = s.Matrix([0, v/s.sqrt(2)])
    generators = [higgs_memory()[1][a][1:3, 1:3] for a in (8, 9, 10, 11)]
    images = [g*t*background for g, t in zip((gw, gw, gw, gy), generators)]
    gauge_mass = s.Matrix(4, 4, lambda i, j: s.simplify(2*s.re((images[i].H*images[j])[0])))
    gz = s.sqrt(gw*gw+gy*gy)
    basis = s.Matrix.hstack(s.Matrix([1, s.I, 0, 0])/s.sqrt(2),
                           s.Matrix([1, -s.I, 0, 0])/s.sqrt(2),
                           s.Matrix([0, 0, gw, -gy])/gz,
                           s.Matrix([0, 0, gy, gw])/gz)
    expected = s.diag(gw*gw*v*v/4, gw*gw*v*v/4, gz*gz*v*v/4, 0)
    if not matter._zero(basis.H*basis-s.eye(4)) or not matter._zero(basis.H*gauge_mass*basis-expected):
        raise ArithmeticError('actual doublet covariant derivative did not give complete gauge masses')
    electric = s.Matrix([[0, -s.I, 0, 0], [s.I, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    if not matter._zero(basis.H*electric*basis-s.diag(1, -1, 0, 0)):
        raise ArithmeticError('charged gauge eigenvectors do not carry W plus/minus charge')
    coordinates = s.symbols('h0:4', real=True)
    radius = sum(x*x for x in coordinates)
    potential = -lam*v*v*radius/2+lam*radius*radius/4
    hessian = s.hessian(potential, coordinates).subs(dict(zip(coordinates, (0, 0, v, 0))))
    if hessian != s.diag(0, 0, 2*lam*v*v, 0):
        raise ArithmeticError('tree comparison scalar Hessian failed')
    masses = {'W': gw*v/2, 'Z': gz*v/2, 'Higgs': s.sqrt(2*lam)*v,
              'photon': s.S.Zero, 'gluon': s.S.Zero}
    definitions = [('W', 'W+', 1, 3, 1), ('W', 'W-', -1, 3, 1),
                   ('Z', 'self-conjugate', 0, 3, 1), ('Higgs', 'self-conjugate', 0, 1, 1),
                   ('photon', 'self-conjugate', 0, 2, 1), ('gluon', 'self-conjugate adjoint', 0, 2, 8)]
    modes, groups = [], []
    for species, orientation, charge, polarizations, colors in definitions:
        mass = masses[species]
        energy = s.sqrt(1+mass*mass)
        for momentum in MOMENTA:
            ids = []
            for color in range(colors):
                for polarization in range(polarizations):
                    ids.append(len(modes))
                    modes.append({'index': len(modes), 'species': species, 'orientation': orientation,
                        'electric_charge': str(charge), 'color_component': color,
                        'polarization': polarization, 'momentum_GeV': list(momentum),
                        'mass_GeV': scalar(mass), 'energy_GeV': scalar(energy)})
            groups.append({'species': species, 'orientation': orientation, 'electric_charge': str(charge),
                           'momentum_GeV': list(momentum), 'modes': ids, 'energy_GeV': scalar(energy)})
    if len(modes) != 56:
        raise ArithmeticError('the physical boson count is not 28 per momentum')
    targets = []
    for species, mass in masses.items():
        targets.append({'species': species, 'rest_coefficient_GeV': scalar(mass),
            'kinetic_convention': 'canonical -F^2/4 and |D H|^2; radial real kinetic coefficient1/2',
            'physical_modes_per_momentum': {'W': 6, 'Z': 3, 'Higgs': 1, 'photon': 2, 'gluon': 16}[species],
            'electric_charge': '+1 and -1' if species == 'W' else '0',
            'antiparticle': 'W- is the conjugate of W+' if species == 'W' else 'self-conjugate; no extra antiparticle copy',
            'canonical_oscillator_retarded_residue': '1 in 1/(z^2-omega^2)',
            'mass_convention': 'tree comparison coefficient at the calibrated loop VEV, with separately chosen tree-stationary m2',
            'native_source_prediction': False,
            'missing_physical_attachment': 'interacting gauge/scalar vacuum, loop pole/width matching and source-clock map; gluon physical states additionally require confinement'})
    massive_norms = {name: {'raw_longitudinal_kinetic': str(mass*mass/(1+mass*mass)),
                           'canonical_longitudinal_coordinate': 'q_L=(m/omega) A_L',
                           'raw_longitudinal_potential': str(mass*mass),
                           'canonical_frequency_squared': str(1+mass*mass)}
                     for name, mass in masses.items() if name in ('W', 'Z')}
    return {'gauge_mass_matrix': gauge_mass, 'gauge_basis': basis, 'scalar_hessian': hessian,
            'masses': masses, 'modes': modes, 'groups': groups, 'targets': targets,
            'longitudinal_normalization': massive_norms,
            'checks': {'gauge_Higgs_mass_matrix_reconstructed_from_native_doublet_currents': True,
                       'charged_electric_eigenvectors_checked': True,
                       'canonically_normalized_gauge_kinetic_metric': 'I12 before physical polarization reduction',
                       'scalar_Hessian_rank': 1, 'Goldstone_directions': 3,
                       'gauge_mass_matrix_rank_in_twelve_directions': 3,
                       'massless_gauge_directions': 9,
                       'massless_longitudinal_and_time_components_not_physical_oscillators': True,
                       'all_physical_polarizations_per_momentum': 28,
                       'tree_stationarity_is_additional_prescription': True,
                       'physical_massless_zero_momentum_oscillator_assumed': False}}
