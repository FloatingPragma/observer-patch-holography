"""Concrete model energies and a conditional, dimensionally explicit mass map.

No numerical mass calibration is fitted. Exact intervals come from the
independent positive-measure certificate. Sparse ground values are labeled
diagnostics separately from a rigorous variational ground enclosure.
"""
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
import sympy as s

from ..cube_eta_protection import cube
from . import vacuum


ROOT = Path(__file__).resolve().parents[2]
CALIBRATION_SOURCES = (
    'model/clock_defect_energy/si_units.py',
    'model/defect_field_clock/README.md',
    'model/defect_field_clock/clock.py',
    'model/particle_comparison/README.md',
    'model/particle_comparison/observations.json',
    'model/moving_bound_pair/README.md',
)


def interval_strings(interval):
    return [str(x) for x in interval]


@lru_cache(maxsize=1)
def ground_certificate():
    """Rigorous bracket for the N6 minimum; separate sparse eigen diagnostics."""
    h = cube.restrict(cube.operators().hubbard_H.real, 1)
    if (h@h-3*sparse.eye(16)).nnz or np.any(h.diagonal()):
        raise ArithmeticError('the exact cube one-particle edge changed')
    up = [j for j, mask in enumerate(cube.weights(1)) if int(mask) & 0xAAAA == 0]
    down = [j for j in range(16) if j not in up]
    a = h[up][:, up]
    if h[up][:, down].nnz or (a@a-3*sparse.eye(8)).nnz or int(a.diagonal().sum()) != 0:
        raise ArithmeticError('the polarized variational sector changed')
    # h/16 has eight negative orbitals; six particles cannot have free
    # energy below -6sqrt3/16, and positive contact cannot lower that.
    # Six spin-up particles occupy four negative plus two positive orbitals,
    # have no double occupations, and realize exactly -2sqrt3/16.
    root = vacuum.root_interval((1, 0, -3), (1, 2))
    rigorous = (-3*root[1]/8, -root[0]/8)
    h6 = cube.restrict(cube.operators().hubbard_H.real, 6).astype(float)/16
    values, vectors = eigsh(h6, k=2, which='SA', tol=1e-12,
                            v0=np.random.default_rng(2026).normal(size=8008))
    order = np.argsort(values); values, vectors = values[order], vectors[:, order]
    residuals = np.linalg.norm(h6@vectors-vectors*values, axis=0)
    if max(residuals) > 1e-10 or not float(rigorous[0]) < values[0] < float(rigorous[1]):
        raise ArithmeticError('numerical N6 ground diagnostic exceeds its declared admission')
    return {'number_sector': 6, 'number_sector_dimension': 8008,
            'rigorous_ground_bounds_formula': '-3sqrt3/8 <= E0_N6 <= -sqrt3/8',
            'rigorous_ground_interval': interval_strings(rigorous),
            'lower_bound_reason': 'six free one-particle energies are bounded below by -sqrt3/16 each; contact is positive',
            'upper_bound_reason': 'six fully spin-polarized orbitals have exact energy -sqrt3/8 and zero contact',
            'numerical_lowest_eigenvalue_estimate': float(values[0]),
            'numerical_next_eigenvalue_estimate': float(values[1]),
            'numerical_eigenpair_residuals': residuals.tolist(),
            'numerical_eigenpair_admission': 1e-10,
            'sparse_estimate_is_an_interval_certified_ground': False,
            'full_Fock_ground_across_all_number_sectors_determined': False}


def mass_certificate():
    energy_unit, length_unit, hbar, gamma, momentum = s.symbols('E_unit L_unit hbar gamma p', positive=True)
    mass = 2*hbar**2/(energy_unit*length_unit**2)
    dimensionless_p = momentum*length_unit/hbar
    converted = s.expand(energy_unit*(dimensionless_p**2/4-gamma**2))
    if s.simplify(converted-(momentum**2/(2*mass)-energy_unit*gamma**2)) != 0:
        raise ArithmeticError('conditional inertial-mass dimensions changed')
    return {'model_dispersion': 'epsilon(P)=|P|^2/4-gamma^2',
            'generator_and_time_orientation': 'global reflected J_a=6/a^2-H_a; the original H_a curvature has the opposite sign',
            'conditional_momentum_map': 'p_phys=hbar P/L_unit',
            'conditional_energy_map': 'E_phys=E_unit epsilon+E_offset, E_unit>0',
            'required_time_map': 't_phys=T_unit t, E_unit T_unit=hbar',
            'conditional_inertial_mass': 'M=2 hbar^2/(E_unit L_unit^2)',
            'model_inertial_mass_hbar_model_equal_one': '2',
            'mass_ratio_to_single_edge_particle': '2',
            'conditional_positive_binding_energy': 'E_binding=E_unit gamma^2',
            'zero_momentum_energy': 'E_offset-E_unit gamma^2',
            'negative_binding_offset_is_rest_mass': False,
            'relativistic_rest_energy_offset_fixed': False,
            'SI_or_GeV_mass_identified': False}


def calibration_audit():
    observations = json.loads((ROOT/CALIBRATION_SOURCES[4]).read_text())
    assigned = observations['electron_calibration']
    if assigned['role'] != 'one common energy calibration input; not an absolute mass prediction':
        raise ArithmeticError('the inherited electron calibration role changed')
    return {'sources_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in CALIBRATION_SOURCES},
            'old_SI_helper': 'requires supplied seconds per model unit, signal speed and E=hbar omega; it derives no physical clock calibration',
            'old_defect_clock_scale': 'Omega_ref=1/4 GeV is a declared calibration of a separate co-driven calibrated field/clock design',
            'old_defect_clock_lab_calibration_proved': False,
            'old_dimer_scale': {'assigned_electron_energy': assigned['value'], 'unit': assigned['unit'],
                               'role': assigned['role']},
            'old_scale_inherited_by_cube_or_moving_pair': False,
            'same_cube_or_moving_pair_physical_clock_and_length_join_proved': False,
            'conditional_reuse_is_an_additional_hypothesis_not_a_prediction': True,
            'measured_mass_used_to_fit_new_candidates': False}


def certificate():
    spectral = vacuum.certificate(); ground = ground_certificate()
    ground_interval = tuple(map(F, ground['rigorous_ground_interval']))
    rows = []
    for row in spectral['spectral_rows']:
        z = tuple(map(F, row['z_interval']))
        energy = tuple(2+x/16 for x in z)
        gap = tuple(x-3 for x in energy)
        weight = tuple(map(F, row['positive_normalized_probe_weight_interval']))
        if list(map(str, gap)) != row['omega_interval']:
            raise ArithmeticError('the common reference energy was not retained')
        rows.append({'factor': row['factor'], 'z_interval': interval_strings(z),
                     'H_energy_interval': interval_strings(energy),
                     'H_energy_display_interval_12dp': interval_strings(vacuum.outward_decimal(energy)),
                     'H_energy_approximation': float(sum(energy)/2),
                     'signed_gap_to_R_interval': interval_strings(gap),
                     'signed_gap_to_R_display_interval_12dp': interval_strings(vacuum.outward_decimal(gap)),
                     'signed_gap_to_R_approximation': float(sum(gap)/2),
                     'normalized_probe_weight_interval': interval_strings(weight),
                     'normalized_probe_weight_approximation': float(sum(weight)/2),
                     'rigorous_gap_to_N6_ground_interval': interval_strings((energy[0]-ground_interval[1], energy[1]-ground_interval[0])),
                     'gap_to_numerically_estimated_N6_ground': float(sum(energy)/2)-ground['numerical_lowest_eigenvalue_estimate']})
    rows.sort(key=lambda row: F(row['H_energy_interval'][0]))
    minimum = tuple(map(F, rows[0]['H_energy_interval']))
    for index, row in enumerate(rows):
        energy = tuple(map(F, row['H_energy_interval']))
        difference = (F(), F()) if index == 0 else (energy[0]-minimum[1], energy[1]-minimum[0])
        row['gap_to_minimum_probe_energy_interval'] = interval_strings(difference)
    return {'model_energy_unit': 'the original cube contact coupling U=1; no physical joule or GeV calibration',
            'cube_reference_energy_H': '3', 'spectral_variable': 'E_H=2+z/16',
            'normalization_of_weights': 'normalized centered density probe; the local phase has15/16 times these weights and an additional1/16 elastic line at E=3',
            'exact_positive_probe_lines': rows,
            'energy_gap_and_weight_intervals_are_rationally_certified': True,
            'minimum_probe_energy_is_full_N6_ground': False,
            'ground_analysis': ground,
            'inertial_mass_and_binding': mass_certificate(),
            'physical_unit_calibration_audit': calibration_audit(),
            'physical_rest_masses_or_particle_names_assigned': False}
