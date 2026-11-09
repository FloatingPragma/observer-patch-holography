"""Exact logarithmic capacity bounds, with decimal values only for display."""
from decimal import Decimal, localcontext
from pathlib import Path
import hashlib

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_physical_boundary_state import boundary
from model.native_quantum_gauss import algebra as native

ROOT = Path(__file__).resolve().parents[2]
SOURCE_INPUTS = (
    'model/particle_primary_capacity/capacity.py',
    'model/particle_primary_capacity/PROOF.md',
    'model/particle_primary_capacity/test_capacity.py',
    'model/particle_primary_capacity/__main__.py',
    'model/particle_common_contract/PROOF.md',
    'model/particle_spatial_ground_state/PROOF.md',
    'model/particle_physical_boundary_state/PROOF.md',
    'model/observer/CAPACITY.md',
)


def log_expression(terms):
    """The integer coefficients of natural logarithms are the exact result."""
    if any(type(base) is not int or base < 2 or type(coefficient) is not int
           or coefficient < 0 for base, coefficient in terms):
        raise ValueError('positive integer bases and nonnegative integer powers required')
    with localcontext() as context:
        context.prec = 40
        value = sum((Decimal(base).ln()*coefficient for base, coefficient in terms), Decimal(0))
    return {'exact_log_terms': [{'base': base, 'coefficient': coefficient}
                                for base, coefficient in terms],
            'approximate_nats_for_display': str(value),
            'decimal_is_not_an_exact_capacity_or_error_certificate': True}


def singlet_rows(generators=None):
    """Each native M15 exterior top form is a gauge singlet, as is its vacuum."""
    generators = native.native_currents() if generators is None else tuple(generators)
    if not generators or any(g.shape != (15, 15) or g+g.H != s.zeros(15)
                             or s.trace(g) != 0 for g in generators):
        raise ValueError('the actual traceless anti-Hermitian M15 currents are required')
    return {'generators_checked': len(generators), 'top_wedge_generator_traces': ['0']*len(generators),
            'Dirac_rows_per_cell': 4, 'orthogonal_local_singlet_witnesses': 16,
            'constant_links_preserve_all_Gauss_constraints': True,
            'local_witness_span_claimed_Hamiltonian_reducing': False}


def member_capacity(member):
    contract.verify_member(member)
    control = singlet_rows()
    source = member['source']; volume = source['side']**3
    degree = source['per_link_word_degree']; links = 3*volume
    one = boundary.native_degree_one_control()
    # This value is independently reconstructed from the five native blocks.
    dimension = 1+2*sum(len(indices)**2 for _, indices in native.native_blocks())
    if dimension != 119:
        raise ArithmeticError('the complete degree-one link dimension differs')
    del one  # Its exact construction also validates all endpoint charges.
    ambient = log_expression([(2, 60*volume), (119, links)]) if degree == 1 else None
    upper = log_expression([(2, 60*volume), (451, degree*links)])
    bits = volume*(60+27*degree)
    if bits != member['preparation_plan']['source']['dimension_bits_upper']:
        raise ArithmeticError('capacity bound differs from the actual source plan')
    return {'schema': 'oph.primary_capacity_bound.v1', 'member_id': member['member_id'], **member['ids'],
            'level': source['level'], 'cells': volume, 'side': source['side'],
            'per_link_word_degree': degree, 'owned_links': links,
            'capacity_definition': 'N=ln(dim carrier), natural logarithm in nats',
            'source_spacing': source['spacing'],
            'physical_coordinate_extent': str(source['side']*s.Rational(source['spacing'])),
            'Hamiltonian_factor_from_inverse_spacing': 2**source['level'],
            'dimensionless_eigenvalues_can_change_with_volume_and_cutoff': True,
            'same_SI_clock_conversion_must_be_tracked_across_levels': True,
            'ambient_capacity_exact': ambient,
            'ambient_capacity_upper': upper,
            'physical_capacity_lower': log_expression([(2, 4*volume)]),
            'physical_capacity_upper': ambient if ambient is not None else upper,
            'exact_physical_Gauss_dimension_computed': False,
            'singlet_witness': control,
            'allocated_source_qubit_count_upper': bits,
            'allocated_source_plus_conjugate_reference_qubit_count_upper': 2*bits,
            'apparatus_and_retained_record_memory_included': False,
            'prepared_state_entropy_computed': False,
            'correctable_public_record_capacity_computed': False,
            'cosmological_capacity_identification_proved': False,
            'full_spatial_state_vector_executed': False,
            'mass_capacity_dependence_computed': False,
            'dimension_not_classical_statevector_memory': True}


def certificate():
    paths = provenance.definition_inputs(SOURCE_INPUTS)
    pins = provenance.definition_hashes(paths)
    result = {'schema': 'oph.primary_capacity_family.v1',
              'members': [member_capacity(contract.member(n)) for n in range(6)],
              'capacity_one_trillion_nats_selected': False,
              'mass_observables_must_retain_refinement_and_unit_provenance': True,
              'input_sha256': pins}
    if provenance.definition_hashes(paths) != pins:
        raise ValueError('capacity source changed during replay')
    return result
