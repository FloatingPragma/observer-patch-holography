"""Finite resource and complete-channel errors for existing native services."""
from fractions import Fraction as F


def generic_instrument_plan(spins, kraus_branches, error):
    """All computable finite Kraus menus on q spins through source Naimark words.

    The constructive proof specifies each primitive, all register bits and
    its exact finite synthesis bound. It uses the inherited refined dyadic
    control interface, not a limit of the binary64 host executor.
    """
    if (type(spins) is not int or spins < 1 or type(kraus_branches) is not int or kraus_branches < 1
            or type(error) is not F or not 0 < error < 1):
        raise ValueError('positive spin and fine Kraus-branch counts and exact error in (0,1) required')
    record_bits = (kraus_branches-1).bit_length()
    bits = spins+record_bits
    dimension = 1 << bits
    # Complex state preparation <=8D; each Householder <=20D;
    # at most D reflectors plus a diagonal phase, and spin code wrappers.
    pulses = 20*dimension*dimension+4*dimension+10*spins
    # The ideal unitary approximation and all source ports receive separate
    # half-trace budgets. Complete source operator error <=29N 2^-p.
    port_bits = 1
    while F(29*pulses, 1 << port_bits) > error/2:
        port_bits += 1
    angle_bits = 1
    while F(2*dimension*(bits+2), 1 << angle_bits) > error/4:
        angle_bits += 1
    return {'spin_bits': spins, 'retained_kraus_branch_count': kraus_branches,
            'coarse_classical_outcome_count_is_not_the_environment_dimension': True,
            'coherent_record_bits': record_bits,
            'decoded_system_and_record_bits': bits,
            'physical_code_plus_record_bits': 2*spins+record_bits,
            'unitary_dimension': dimension,
            'native_pulse_upper_bound': pulses,
            'normalized_factor_vector_error_budget': str(error/F(32*dimension)),
            'factor_phase_error_budget': str(error/F(8*(dimension+1))),
            'rank_one_tree_angle_bits': angle_bits,
            'rank_one_tree_operator_error_upper': str(F(2*dimension*(bits+2), 1 << angle_bits)),
            'integrated_port_fraction_bits': port_bits,
            'port_operator_error_upper': str(F(29*pulses, 1 << port_bits)),
            'complete_halfdiamond_error_upper': str(error),
            'arbitrary_untouched_reference': True,
            'coherent_addressed_source_service_required': True,
            'full_host_execution_at_arbitrary_size_claimed': False,
            'automatic_access_from_bare_repair_claimed': False}
