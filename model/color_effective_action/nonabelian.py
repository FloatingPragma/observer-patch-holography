"""Source normalization and exact sufficient coordinate-regulator bounds."""
from fractions import Fraction as F

from model.interacting_extension.higgs_join.interacting_carrier import NativeCarrier
from .action import Elimination, interval_packet


def fierz_certificate(carrier):
    import sympy as s
    rows = []
    for dimension, indices, generators in (
        (3, carrier.data['blocks'][0]['indices'],
         tuple(s.I*g/2 for g in carrier.anti[:7])+(s.I*carrier.anti[7]/(2*s.sqrt(3)),)),
        (2, carrier.data['blocks'][1]['indices'], tuple(s.I*g/2 for g in carrier.anti[8:11])),
    ):
        t = [g.extract(indices, indices) for g in generators]
        swap = s.zeros(dimension**2)
        for i in range(dimension):
            for j in range(dimension):
                swap[i*dimension+j, j*dimension+i] = 1
        tensor = sum((s.kronecker_product(g, g) for g in t), s.zeros(dimension**2))
        expected = swap/2-s.eye(dimension**2)/(2*dimension)
        if tensor.applyfunc(s.simplify) != expected:
            raise ArithmeticError('source-normalized Fierz contraction failed')
        casimir = sum((g*g for g in t), s.zeros(dimension)).applyfunc(s.simplify)
        value = s.Rational(dimension**2-1, 2*dimension)
        if casimir != value*s.eye(dimension):
            raise ArithmeticError('source fundamental Casimir failed')
        rows.append({'dimension':dimension, 'source_fundamental_Casimir':str(value),
                     'full_tensor_matrix_entries_checked':dimension**4,
                     'tensor_contraction':'swap/2-identity/(2*dimension)'})
    return rows


def highest_Q_power(n):
    if type(n) is not int or n < 0:
        raise ValueError('nonnegative integer tensor degree required')
    return {'tensor_degree':n, 'charge6':n,
            'representation_dimension':(n+1)**2*(n+2)//2,
            'native_Casimir':str(F(n*(7*n+18), 8))}


def growth_certificate():
    carrier = NativeCarrier()
    blocks = carrier.data['blocks']
    charges = [int(b['q6Y']) for b in blocks]
    dimensions = [b['dimension'] for b in blocks]
    if charges != [-4, -3, 1, 2, 6] or dimensions != [3, 2, 6, 3, 1]:
        raise ArithmeticError('source representation menu changed')
    if set(charges)&{-q for q in charges}:
        raise ArithmeticError('source conjugate absence changed')
    native_dimension = 1+sum(d*d for d in dimensions)
    conjugate_closed = 1+2*sum(d*d for d in dimensions)
    q_block = blocks[charges.index(1)]
    if q_block['dimension'] != highest_Q_power(1)['representation_dimension']:
        raise ArithmeticError('actual Q highest-weight seed differs')
    if {str(carrier.data['casimir'][i, i]) for i in q_block['indices']} != {'25/8'}:
        raise ArithmeticError('source Q Casimir differs')
    # On the conditional full coordinate realization, each native raising
    # matrix has norm sqrt(d). The two-edge trace has norm at most d².
    magnetic = 2*sum(b['dimension']**2 for b in blocks if b['color_casimir'])
    # Per edge:3families*2spins*(four doublet entries plus one singlet).
    # sqrt(2)<=3/2, two orientations and the adjoint are all retained.
    hopping = 2*2*3*2*(4*F(3, 2)+1)
    # Two Higgs components, two spin signs, three diagonal families, adjoint.
    if len(carrier.gamma) != 2 or any(abs(gamma) != 1 for _, _, gamma in carrier.gamma):
        raise ArithmeticError('literal Higgs tensor component count differs')
    higgs = 2*2*len(carrier.gamma)*3
    perturbation = magnetic+hopping+higgs+1
    if (native_dimension, conjugate_closed, magnetic, hopping, higgs, perturbation) != (60, 119, 108, 168, 24, 301):
        raise ArithmeticError('complete source coupling bound differs')
    return {'actual_charge_menu':charges, 'actual_dimensions':dimensions,
            'source_normalized_Fierz_identities':fierz_certificate(carrier),
            'native_link_dimension':native_dimension,
            'minimum_dimension_after_all_conjugates_added':conjugate_closed,
            'all_native_conjugate_blocks_absent':True,
            'zero_electric_energy_vacuum_coordinate_unitarity_defect':'1',
            'Q_tensor_highest_weight_family':[highest_Q_power(n) for n in range(9)],
            'conditional_same_control_rotor_bounds':{
                'two_edge_magnetic':str(magnetic), 'full_lepton_hopping':str(hopping),
                'Higgs_creation_and_adjoint':str(higgs), 'scalar_energy':'1',
                'total_bounded_perturbation':'301',
                'max_square_root_Casimir_change_per_link_per_insertion':'2'},
            'actual_native60_source_selects_this_completion':False,
            'capacity_selected_representation_cutoff_supplied':False}


def cutoff_bounds(cutoff_sqrt, initial_sqrt=F(0), mu=F(602), interaction=F(301), step=F(2)):
    values = (cutoff_sqrt, initial_sqrt, mu, interaction, step)
    if any(type(x) not in (int, F) for x in values):
        raise TypeError('exact rational Casimir-radius inputs required')
    radius, inner, mu, bound, step = map(F, values)
    if not 0 <= inner <= radius or not 0 < bound < mu or step <= 0:
        raise ValueError('positive regulator and resolvent margins required')
    # n insertions stay within radius inner+n*step. Exit needs strict >radius.
    count = (radius-inner)//step+1
    vector = 2/(mu-bound)*(bound/mu)**count
    local = 2/(mu-bound)*(bound/mu)**(2*count)
    return {'Casimir_cutoff':str(radius*radius), 'initial_Casimir_cutoff':str(inner*inner),
            'minimum_exit_insertions':count,
            'conditional_full_to_cutoff_resolvent_bound':str(vector),
            'conditional_local_action_resolvent_bound':str(local),
            'complete_coordinate_family_premise_required':True,
            'source_capacity_selection_inferred':False}


def coordinate_core_error(cutoff_sqrt, representation_sqrt, mean_electric):
    if any(type(x) not in (int, F) for x in (cutoff_sqrt, representation_sqrt, mean_electric)):
        raise TypeError('exact rational coordinate-radius inputs required')
    radius, rep, energy = map(F, (cutoff_sqrt, representation_sqrt, mean_electric))
    if not radius > rep >= 0 or energy < 0:
        raise ValueError('positive coordinate cutoff margin required')
    return min(F(1), energy/(radius-rep)**2)


def completion_response_witness(time=F(1, 2**20)):
    if type(time) not in (int, F) or time <= 0:
        raise ValueError('positive exact source time required')
    time = F(time)
    native = Elimination()
    if native.source.norm_bound >= 24:
        raise ArithmeticError('finite source norm admission changed')
    seed = {native.source.seed:F(1)}
    response = native.source.apply(seed)
    flux = {i for i, row in enumerate(native.source.labels)
            if row['link0_block'] != -1 or row['link1_block'] != -1}
    scalar = native.q
    native_flux = native.source.inner(response, response, flux)
    native_higgs = native.source.inner(response, response, scalar)
    if (native_flux, native_higgs, native.source.inner(response, response)) != (54, 12, 66):
        raise ArithmeticError('complete native blank response differs')
    rows = []
    for name, flux_coefficient, amplitude_upper, remainder in (
        ('native_raising_source', F(54), F(8), F(108)),
        ('full_coordinate_completion', F(108), F(11), F(1700)),
    ):
        flux_error = 2*amplitude_upper*remainder*time+remainder**2*time**2
        higgs_error = 8*remainder*time+remainder**2*time**2
        flo, fhi = flux_coefficient-flux_error, flux_coefficient+flux_error
        hlo, hhi = F(12)-higgs_error, F(12)+higgs_error
        if min(flo, hlo) <= 0:
            raise ValueError('source time exceeds the response comparison margin')
        lo, hi = flo/hhi, fhi/hlo
        rows.append({'source':name, 'flux_probability_over_time_squared':interval_packet(flux_coefficient, flux_error),
                     'Higgs_probability_over_time_squared':interval_packet(F(12), higgs_error),
                     'complete_flux_over_Higgs_probability_ratio':interval_packet((lo+hi)/2, (hi-lo)/2),
                     'state_vector_remainder_over_time_squared':str(remainder)})
    if F(rows[0]['complete_flux_over_Higgs_probability_ratio'][1]) >= F(rows[1]['complete_flux_over_Higgs_probability_ratio'][0]):
        raise ArithmeticError('finite complete records do not distinguish the attachments')
    return {'same_empty_input':True, 'same_default_coefficients':True, 'same_source_time':str(time),
            'original_flux_leading_coefficient':'54', 'completed_flux_leading_coefficient':'108',
            'Higgs_leading_coefficient_in_both':'12', 'finite_time_experiments':rows,
            'identical_complete_source_records':False,
            'different_low_energy_universality_class_proved':False}
